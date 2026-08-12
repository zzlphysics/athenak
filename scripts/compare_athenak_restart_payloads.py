#!/usr/bin/env python3
"""Compare topology metadata and opaque payloads of two closed AthenaK restarts.

The ParameterInput dumps at the beginning of the files are deliberately not
compared byte-for-byte.  Each dump is parsed by
``read_athenak_restart_metadata.py`` and only restart topology/state metadata is
compared.  Everything after each file's independently determined metadata
boundary is then hashed and compared in bounded chunks.  That opaque payload
contains serialized object state, data-size records, and field arrays.

Inputs must be closed files.  Path and open-file statistics are checked before,
during, and after inspection; a file that changes makes the comparison fail
closed with exit status 2.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import sys
from typing import BinaryIO, Iterable


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
if str(SCRIPT_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

from read_athenak_restart_metadata import (  # noqa: E402
    AMR_COUNTER_MAGIC,
    AMR_COUNTER_VERSION,
    EVENT_COUNTER_MAGIC,
    EVENT_COUNTER_VERSION,
    EVENT_SUM_COUNTER_COUNT,
    RestartMetadata,
    read_restart_metadata,
)


DEFAULT_CHUNK_BYTES = 1024 * 1024
MAX_CHUNK_BYTES = 64 * 1024 * 1024


@dataclass(frozen=True)
class FileSignature:
    """Identity and mutation-sensitive attributes of one local file."""

    device: int
    inode: int
    size: int
    mtime_ns: int
    ctime_ns: int


def _signature_from_stat(stat_result: os.stat_result) -> FileSignature:
    return FileSignature(
        device=stat_result.st_dev,
        inode=stat_result.st_ino,
        size=stat_result.st_size,
        mtime_ns=stat_result.st_mtime_ns,
        ctime_ns=stat_result.st_ctime_ns,
    )


def _file_signature(path: Path) -> FileSignature:
    return _signature_from_stat(path.stat())


def _assert_path_unchanged(
    path: Path, expected: FileSignature, checkpoint: str
) -> None:
    if _file_signature(path) != expected:
        raise RuntimeError(
            f"restart changed {checkpoint}: {path}; compare only closed files"
        )


def _assert_stream_unchanged(
    stream: BinaryIO, path: Path, expected: FileSignature, checkpoint: str
) -> None:
    if _signature_from_stat(os.fstat(stream.fileno())) != expected:
        raise RuntimeError(
            f"open restart changed {checkpoint}: {path}; compare only closed files"
        )


def _parameter_int(metadata: RestartMetadata, block: str, key: str) -> int:
    try:
        return int(metadata.parameters[block][key])
    except (KeyError, ValueError) as exc:
        raise ValueError(
            f"restart parameter <{block}>/{key} is missing or invalid"
        ) from exc


def _parameter_float(metadata: RestartMetadata, block: str, key: str) -> float:
    try:
        return float(metadata.parameters[block][key])
    except (KeyError, ValueError) as exc:
        raise ValueError(
            f"restart parameter <{block}>/{key} is missing or invalid"
        ) from exc


def _root_mesh_signature(metadata: RestartMetadata) -> tuple[object, ...]:
    """Return fixed-header mesh values without exposing ParameterInput text."""

    return (
        _parameter_int(metadata, "mesh", "nghost"),
        *(
            _parameter_int(metadata, "mesh", f"nx{axis}")
            for axis in range(1, 4)
        ),
        *(
            _parameter_int(metadata, "meshblock", f"nx{axis}")
            for axis in range(1, 4)
        ),
        *(
            _parameter_float(metadata, "mesh", f"x{axis}{bound}")
            for axis in range(1, 4)
            for bound in ("min", "max")
        ),
    )


def _compare_metadata(
    left: RestartMetadata, right: RestartMetadata
) -> dict[str, object]:
    """Compare restart topology metadata, excluding ParameterInput text."""

    fields = {
        "time": left.time == right.time,
        "last_dt": left.last_dt == right.last_dt,
        "cycle": left.cycle == right.cycle,
        "nmb_total": left.nmb_total == right.nmb_total,
        "root_level": left.root_level == right.root_level,
        "root_mesh": _root_mesh_signature(left) == _root_mesh_signature(right),
        "locations": left.locations == right.locations,
        "costs": left.costs == right.costs,
        "amr_cycle_counters": (
            left.amr_cycle_counters == right.amr_cycle_counters
        ),
        "event_counters": left.event_counters == right.event_counters,
        "byte_order": left.byte_order == right.byte_order,
        "real_bytes": left.real_bytes == right.real_bytes,
    }
    fields["abi"] = bool(fields["byte_order"] and fields["real_bytes"])
    return {
        "match": all(bool(value) for value in fields.values()),
        "fields": fields,
    }


def _read_expected(stream: BinaryIO, size: int, path: Path) -> bytes:
    if size == 0:
        return b""
    data = stream.read(size)
    if len(data) != size:
        raise RuntimeError(
            f"payload stream for {path} ended early: expected {size} bytes, "
            f"read {len(data)}"
        )
    return data


def _first_chunk_difference(left: bytes, right: bytes) -> int | None:
    common = min(len(left), len(right))
    if left[:common] != right[:common]:
        for index, (left_byte, right_byte) in enumerate(zip(left, right)):
            if left_byte != right_byte:
                return index
        raise AssertionError("unequal byte slices had no unequal element")
    if len(left) != len(right):
        return common
    return None


def _stream_payloads(
    left_path: Path,
    right_path: Path,
    left_signature: FileSignature,
    right_signature: FileSignature,
    left_offset: int,
    right_offset: int,
    chunk_bytes: int,
) -> dict[str, object]:
    """Hash and compare both payloads without retaining more than two chunks."""

    left_expected = left_signature.size - left_offset
    right_expected = right_signature.size - right_offset
    if left_expected < 0 or right_expected < 0:
        raise ValueError("restart metadata boundary lies beyond end of file")

    left_hash = hashlib.sha256()
    right_hash = hashlib.sha256()
    left_read = 0
    right_read = 0
    comparison_offset = 0
    first_difference: int | None = None

    with left_path.open("rb") as left_stream, right_path.open("rb") as right_stream:
        _assert_stream_unchanged(
            left_stream, left_path, left_signature, "before payload streaming"
        )
        _assert_stream_unchanged(
            right_stream, right_path, right_signature,
            "before payload streaming",
        )
        left_stream.seek(left_offset)
        right_stream.seek(right_offset)

        while left_read < left_expected or right_read < right_expected:
            left_size = min(chunk_bytes, left_expected - left_read)
            right_size = min(chunk_bytes, right_expected - right_read)
            left_chunk = _read_expected(left_stream, left_size, left_path)
            right_chunk = _read_expected(right_stream, right_size, right_path)
            left_hash.update(left_chunk)
            right_hash.update(right_chunk)
            if first_difference is None:
                local_difference = _first_chunk_difference(
                    left_chunk, right_chunk
                )
                if local_difference is not None:
                    first_difference = comparison_offset + local_difference
            left_read += len(left_chunk)
            right_read += len(right_chunk)
            comparison_offset += max(len(left_chunk), len(right_chunk))

        _assert_stream_unchanged(
            left_stream, left_path, left_signature, "after payload streaming"
        )
        _assert_stream_unchanged(
            right_stream, right_path, right_signature,
            "after payload streaming",
        )

    left_digest = left_hash.hexdigest()
    right_digest = right_hash.hexdigest()
    payload_equal = bool(
        left_read == right_read
        and first_difference is None
        and left_digest == right_digest
    )
    return {
        "payload_equal": payload_equal,
        "payload_bytes": {"left": left_expected, "right": right_expected},
        "sha256": {"left": left_digest, "right": right_digest},
        "first_difference_offset": first_difference,
        "field_data_read": True,
        "field_data_bytes_read": {"left": left_read, "right": right_read},
    }


def compare_restart_payloads(
    left_path: Path,
    right_path: Path,
    *,
    chunk_bytes: int = DEFAULT_CHUNK_BYTES,
) -> dict[str, object]:
    """Compare two stable single-file restart payloads and their topology."""

    if not 1 <= chunk_bytes <= MAX_CHUNK_BYTES:
        raise ValueError(
            f"chunk_bytes must be between 1 and {MAX_CHUNK_BYTES}"
        )
    left_path = left_path.resolve(strict=True)
    right_path = right_path.resolve(strict=True)
    left_signature = _file_signature(left_path)
    right_signature = _file_signature(right_path)

    left_metadata = read_restart_metadata(left_path)
    right_metadata = read_restart_metadata(right_path)
    _assert_path_unchanged(
        left_path, left_signature, "while metadata was being read"
    )
    _assert_path_unchanged(
        right_path, right_signature, "while metadata was being read"
    )

    metadata = _compare_metadata(left_metadata, right_metadata)
    payload = _stream_payloads(
        left_path,
        right_path,
        left_signature,
        right_signature,
        left_metadata.metadata_end,
        right_metadata.metadata_end,
        chunk_bytes,
    )

    _assert_path_unchanged(left_path, left_signature, "during comparison")
    _assert_path_unchanged(right_path, right_signature, "during comparison")
    result: dict[str, object] = {
        "kind": "athenak_restart_payload_comparison",
        "left": str(left_path),
        "right": str(right_path),
        "metadata": metadata,
        **payload,
    }
    result["match"] = bool(
        metadata["match"] and result["payload_equal"]
    )
    return result


def _chunk_bytes(text: str) -> int:
    try:
        value = int(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if not 1 <= value <= MAX_CHUNK_BYTES:
        raise argparse.ArgumentTypeError(
            f"must be between 1 and {MAX_CHUNK_BYTES}"
        )
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument(
        "--chunk-bytes",
        type=_chunk_bytes,
        default=DEFAULT_CHUNK_BYTES,
        help=f"streaming chunk size (default: {DEFAULT_CHUNK_BYTES})",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="emit JSON (the comparator's default output format)",
    )
    parser.add_argument(
        "--compact", action="store_true", help="emit compact rather than pretty JSON"
    )
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = compare_restart_payloads(
            args.left, args.right, chunk_bytes=args.chunk_bytes
        )
    except (OSError, RuntimeError, TypeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.compact:
        print(json.dumps(result, sort_keys=True, separators=(",", ":")))
    else:
        print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["match"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
