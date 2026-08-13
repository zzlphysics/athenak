#!/usr/bin/env python3
"""Stream-audit every stored Real in one closed AthenaK restart file."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import sys
from typing import Any, Sequence

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
if str(SCRIPT_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

from compare_athenak_restart_fields import _derive_layout  # noqa: E402
from output_integrity import (  # noqa: E402
    HASH_CHUNK_BYTES,
    _assert_closed,
    _assert_path_signature,
    _assert_stream_signature,
    _open_regular_nofollow,
)
from read_athenak_restart_metadata import read_restart_metadata_stream  # noqa: E402


def _audit_leaf_topology(metadata: Any) -> dict[str, Any]:
    """Prove that restart logical locations are exactly a dyadic AMR leaf domain."""

    try:
        mesh = metadata.parameters["mesh"]
        meshblock = metadata.parameters["meshblock"]
        num_levels = int(metadata.parameters["mesh_refinement"]["num_levels"])
        mesh_cells = tuple(int(mesh[f"nx{axis}"]) for axis in (1, 2, 3))
        block_cells = tuple(int(meshblock[f"nx{axis}"]) for axis in (1, 2, 3))
    except (KeyError, ValueError) as exc:
        raise ValueError("restart has invalid mesh dimensions for topology audit") from exc
    if (num_levels <= 0 or
            not all(value > 0 for value in mesh_cells + block_cells) or
            any(total % block for total, block in zip(mesh_cells, block_cells))):
        raise ValueError("restart mesh dimensions do not form integral root MeshBlocks")
    root_blocks = tuple(total // block
                        for total, block in zip(mesh_cells, block_cells))
    active = tuple(total > 1 for total in mesh_cells)
    topology: set[tuple[int, int, int, int]] = set()
    maximum_physical_level = 0
    for location in metadata.locations:
        physical_level = location.level - metadata.root_level
        if physical_level < 0 or physical_level >= num_levels:
            raise ValueError(
                "restart logical location lies outside configured physical levels")
        coordinates = (location.lx1, location.lx2, location.lx3)
        limits = tuple((root_blocks[axis] << physical_level)
                       if active[axis] else root_blocks[axis]
                       for axis in range(3))
        if not all(0 <= coordinate < limit
                   for coordinate, limit in zip(coordinates, limits)):
            raise ValueError("restart logical location is outside dyadic root domain")
        identity = (physical_level, *coordinates)
        if identity in topology:
            raise ValueError("restart contains duplicate logical locations")
        topology.add(identity)
        maximum_physical_level = max(maximum_physical_level, physical_level)
    for level, lx1, lx2, lx3 in topology:
        for ancestor_level in range(level):
            shift = level - ancestor_level
            ancestor_coordinates = tuple(
                coordinate >> shift if active[axis] else coordinate
                for axis, coordinate in enumerate((lx1, lx2, lx3)))
            if (ancestor_level, *ancestor_coordinates) in topology:
                raise ValueError(
                    "restart contains overlapping ancestor/descendant MeshBlocks")
    dimensions = sum(active)
    coverage_units = sum(
        1 << (dimensions * (maximum_physical_level - level))
        for level, *_ in topology)
    domain_units = math.prod(root_blocks) * (
        1 << (dimensions * maximum_physical_level))
    if coverage_units != domain_units:
        raise ValueError(
            "restart AMR leaf coverage is incomplete: "
            f"{coverage_units} != {domain_units}")
    return {
        "root_meshblocks": list(root_blocks),
        "active_dimensions": list(active),
        "maximum_physical_level": maximum_physical_level,
        "configured_physical_levels": num_levels,
        "logical_locations_unique": True,
        "ancestor_descendant_overlap": False,
        "coverage_units": coverage_units,
        "domain_units": domain_units,
        "complete_leaf_coverage": True,
    }


def _read_exact(stream: Any, size: int, description: str) -> bytes:
    pieces: list[bytes] = []
    remaining = size
    while remaining:
        piece = stream.read(min(remaining, HASH_CHUNK_BYTES))
        if not piece:
            break
        pieces.append(piece)
        remaining -= len(piece)
    if remaining:
        raise ValueError(
            f"restart ended early while reading {description}: expected {size} "
            f"bytes, got {size - remaining}"
        )
    return pieces[0] if len(pieces) == 1 else b"".join(pieces)


def audit_restart(path: os.PathLike[str] | str) -> dict[str, Any]:
    """Validate layout, EOF, finiteness, closure, identity, and SHA-256."""

    checked_path, stream, signature = _open_regular_nofollow(path)
    digest = hashlib.sha256()
    try:
        exempt = {(os.getpid(), stream.fileno())}
        _assert_closed(checked_path, signature, exempt)
        metadata = read_restart_metadata_stream(
            stream, signature["size"], str(checked_path)
        )
        topology = _audit_leaf_topology(metadata)
        layout = _derive_layout(metadata)
        size_record_bytes = struct.calcsize("Q")
        field_start = metadata.metadata_end + size_record_bytes
        expected_file_size = field_start + metadata.nmb_total * layout.data_size
        if signature["size"] != expected_file_size:
            raise ValueError(
                "restart size disagrees with derived field layout: "
                f"expected {expected_file_size}, found {signature['size']}"
            )
        if layout.data_size % layout.real_bytes:
            raise ValueError("per-MeshBlock data_size is not an integral Real array")

        stream.seek(0)
        header = _read_exact(stream, field_start, "header and data_size record")
        digest.update(header)
        byte_prefix = "<" if metadata.byte_order == "little" else ">"
        stored_data_size = struct.unpack(
            f"{byte_prefix}Q", header[-size_record_bytes:]
        )[0]
        if stored_data_size != layout.data_size:
            raise ValueError(
                "restart data_size mismatch: "
                f"stored {stored_data_size}, derived {layout.data_size}"
            )

        values_per_block = layout.data_size // layout.real_bytes
        finite_count = 0
        nonfinite_count = 0
        nan_count = 0
        positive_infinity_count = 0
        negative_infinity_count = 0
        max_abs = 0.0
        first_nonfinite: dict[str, int | str] | None = None
        for gid in range(metadata.nmb_total):
            block = _read_exact(
                stream, layout.data_size, f"MeshBlock GID {gid}"
            )
            digest.update(block)
            values = np.frombuffer(block, dtype=layout.real_dtype)
            if values.size != values_per_block:
                raise ValueError(
                    f"MeshBlock GID {gid} has {values.size} Reals; "
                    f"expected {values_per_block}"
                )
            finite = np.isfinite(values)
            block_finite = int(np.count_nonzero(finite))
            finite_count += block_finite
            block_nonfinite = values.size - block_finite
            nonfinite_count += block_nonfinite
            if block_finite:
                block_max = float(np.max(np.abs(values[finite])))
                if not math.isfinite(block_max):
                    raise ValueError("finite-value maximum unexpectedly overflowed")
                max_abs = max(max_abs, block_max)
            if block_nonfinite:
                nan_count += int(np.count_nonzero(np.isnan(values)))
                positive_infinity_count += int(np.count_nonzero(np.isposinf(values)))
                negative_infinity_count += int(np.count_nonzero(np.isneginf(values)))
                if first_nonfinite is None:
                    index = int(np.flatnonzero(~finite)[0])
                    value = values[index]
                    classification = (
                        "nan" if np.isnan(value)
                        else "positive_infinity" if value > 0
                        else "negative_infinity"
                    )
                    first_nonfinite = {
                        "gid": gid,
                        "real_index": index,
                        "classification": classification,
                    }

        if stream.read(1):
            raise ValueError("restart contains trailing bytes after the final MeshBlock")
        _assert_stream_signature(stream, checked_path, signature, "during restart audit")
        _assert_path_signature(checked_path, signature, "during restart audit")
        _assert_closed(checked_path, signature, exempt)
    finally:
        stream.close()

    _assert_path_signature(checked_path, signature, "after restart audit")
    _assert_closed(checked_path, signature)
    if nonfinite_count:
        raise ValueError(
            f"restart contains {nonfinite_count} non-finite stored Reals; "
            f"nan={nan_count}, +inf={positive_infinity_count}, "
            f"-inf={negative_infinity_count}, first={first_nonfinite}"
        )

    expected_values = metadata.nmb_total * values_per_block
    if finite_count != expected_values:
        raise ValueError(
            f"stored Real count mismatch: expected {expected_values}, "
            f"found {finite_count}"
        )
    return {
        "kind": "athenak_restart_audit",
        "valid": True,
        "path": str(checked_path),
        "signature": dict(signature),
        "sha256": digest.hexdigest(),
        "metadata": {
            "time": metadata.time,
            "last_dt": metadata.last_dt,
            "cycle": metadata.cycle,
            "nmb_total": metadata.nmb_total,
            "root_level": metadata.root_level,
            "real_bytes": metadata.real_bytes,
            "byte_order": metadata.byte_order,
            "metadata_end": metadata.metadata_end,
            "physical_levels": metadata.physical_level_rows(),
        },
        "layout": {
            "data_size": layout.data_size,
            "field_start": field_start,
            "expected_file_size": expected_file_size,
            "values_per_meshblock": values_per_block,
            "streaming_unit": "one MeshBlock",
        },
        "topology": topology,
        "stored_reals": {
            "count": expected_values,
            "finite_count": finite_count,
            "nonfinite_count": 0,
            "nan_count": 0,
            "positive_infinity_count": 0,
            "negative_infinity_count": 0,
            "max_abs": max_abs,
        },
        "closure_check": "linux_proc_fd",
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("restart", type=Path)
    parser.add_argument(
        "--compact", action="store_true", help="emit compact JSON"
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        result = audit_restart(args.restart)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            result,
            indent=None if args.compact else 2,
            sort_keys=True,
            allow_nan=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
