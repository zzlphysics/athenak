#!/usr/bin/env python3
"""Fail-closed integrity primitives for completed AthenaK output files.

The functions in this module are intentionally Linux-oriented.  A successful
checksum means that the input was a closed, non-symlink regular file whose
device, inode, size, modification time, and change time remained stable for the
whole inspection.  A JSON hand-off record is installed with a same-directory
hard link so an existing record is never silently replaced.
"""

from __future__ import annotations

import errno
import hashlib
import json
import math
import os
from pathlib import Path
import secrets
import stat
import sys
import time
from typing import IO, Any, Mapping


HASH_CHUNK_BYTES = 8 * 1024 * 1024
PROC_ROOT = Path("/proc")


def strict_json_loads(raw: bytes | str, source: str = "JSON") -> Any:
    """Decode JSON while rejecting duplicate keys and non-finite constants."""

    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"{source}: duplicate JSON key {key!r}")
            result[key] = value
        return result

    def no_constant(value: str) -> None:
        raise ValueError(f"{source}: non-finite JSON number {value}")

    text = raw.decode("utf-8") if isinstance(raw, bytes) else raw
    return json.loads(
        text, object_pairs_hook=no_duplicates, parse_constant=no_constant)


def _path(value: os.PathLike[str] | str) -> Path:
    """Return an absolute path without resolving its final symlink."""

    return Path(os.path.abspath(os.path.expanduser(os.fspath(value))))


def _signature(stat_result: os.stat_result) -> dict[str, int]:
    return {
        "device": stat_result.st_dev,
        "inode": stat_result.st_ino,
        "size": stat_result.st_size,
        "mtime_ns": stat_result.st_mtime_ns,
        "ctime_ns": stat_result.st_ctime_ns,
    }


def _assert_regular(stat_result: os.stat_result, path: Path) -> None:
    if not stat.S_ISREG(stat_result.st_mode):
        raise ValueError(f"not a regular file: {path}")


def _lstat_regular(path: Path) -> tuple[os.stat_result, dict[str, int]]:
    stat_result = path.lstat()
    if stat.S_ISLNK(stat_result.st_mode):
        raise ValueError(f"refusing symbolic link: {path}")
    _assert_regular(stat_result, path)
    return stat_result, _signature(stat_result)


def _assert_signature(
    actual_stat: os.stat_result,
    expected: Mapping[str, int],
    path: Path,
    checkpoint: str,
) -> None:
    _assert_regular(actual_stat, path)
    if _signature(actual_stat) != dict(expected):
        raise RuntimeError(f"file changed {checkpoint}: {path}")


def _open_regular_nofollow(
    path: os.PathLike[str] | str,
) -> tuple[Path, IO[bytes], dict[str, int]]:
    """Open *path* read-only and bind it to its lstat identity."""

    checked_path = _path(path)
    _, expected = _lstat_regular(checked_path)
    nofollow = getattr(os, "O_NOFOLLOW", None)
    if sys.platform != "linux" or nofollow is None:
        raise RuntimeError("O_NOFOLLOW integrity checks require Linux")
    flags = os.O_RDONLY | nofollow
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    try:
        descriptor = os.open(checked_path, flags)
    except OSError as exc:
        if exc.errno == errno.ELOOP:
            raise ValueError(f"refusing symbolic link: {checked_path}") from exc
        raise
    try:
        descriptor_stat = os.fstat(descriptor)
        _assert_signature(
            descriptor_stat, expected, checked_path, "between lstat and open"
        )
        stream = os.fdopen(descriptor, "rb", closefd=True)
    except BaseException:
        os.close(descriptor)
        raise
    return checked_path, stream, expected


def _assert_path_signature(
    path: Path, expected: Mapping[str, int], checkpoint: str
) -> None:
    try:
        stat_result = path.lstat()
    except FileNotFoundError as exc:
        raise RuntimeError(f"file disappeared {checkpoint}: {path}") from exc
    if stat.S_ISLNK(stat_result.st_mode):
        raise RuntimeError(f"file became a symbolic link {checkpoint}: {path}")
    _assert_signature(stat_result, expected, path, checkpoint)


def _assert_stream_signature(
    stream: IO[bytes], path: Path, expected: Mapping[str, int], checkpoint: str
) -> None:
    _assert_signature(os.fstat(stream.fileno()), expected, path, checkpoint)


def _process_owner(process: Path) -> int | None:
    try:
        return process.stat().st_uid
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        return None


def _open_holders(
    expected: Mapping[str, int], exempt: set[tuple[int, int]] | None = None
) -> list[tuple[int, int, str]]:
    """Find accessible Linux descriptors referring to *expected*'s inode.

    A disappearing process or descriptor is a normal race.  Inability to inspect
    descriptors owned by our effective uid is not: it makes closure unprovable and
    therefore aborts the audit.  Other-user descriptor directories may be hidden by
    the host's procfs policy; those cannot normally acquire a private campaign file
    and are skipped, while every accessible descriptor (including our own) is checked.
    """

    if sys.platform != "linux" or not (PROC_ROOT / "self" / "fd").is_dir():
        raise RuntimeError("cannot prove file closure because Linux /proc is unavailable")
    try:
        processes = list(PROC_ROOT.iterdir())
    except OSError as exc:
        raise RuntimeError("cannot enumerate Linux /proc to prove file closure") from exc

    exempt = exempt or set()
    holders: list[tuple[int, int, str]] = []
    effective_uid = os.geteuid()
    for process in processes:
        if not process.name.isdigit():
            continue
        pid = int(process.name)
        owner = _process_owner(process)
        descriptors = process / "fd"
        try:
            entries = list(descriptors.iterdir())
        except (FileNotFoundError, ProcessLookupError):
            continue
        except PermissionError as exc:
            if owner == effective_uid:
                raise RuntimeError(
                    f"cannot inspect same-user process {pid} while proving closure"
                ) from exc
            continue
        except OSError as exc:
            raise RuntimeError(
                f"cannot inspect process {pid} while proving file closure"
            ) from exc
        for descriptor in entries:
            try:
                fd = int(descriptor.name)
            except ValueError:
                continue
            if (pid, fd) in exempt:
                continue
            try:
                descriptor_stat = descriptor.stat()
            except (FileNotFoundError, ProcessLookupError):
                continue
            except PermissionError as exc:
                if owner == effective_uid:
                    raise RuntimeError(
                        f"cannot inspect same-user descriptor {pid}/{fd} "
                        "while proving closure"
                    ) from exc
                continue
            except OSError as exc:
                raise RuntimeError(
                    f"cannot inspect descriptor {pid}/{fd} while proving closure"
                ) from exc
            if (
                descriptor_stat.st_dev == expected["device"]
                and descriptor_stat.st_ino == expected["inode"]
            ):
                try:
                    command = (process / "comm").read_text(
                        encoding="utf-8"
                    ).strip()
                except (FileNotFoundError, PermissionError, ProcessLookupError):
                    command = "unknown"
                holders.append((pid, fd, command))
    return sorted(holders)


def _assert_closed(
    path: Path,
    expected: Mapping[str, int],
    exempt: set[tuple[int, int]] | None = None,
) -> None:
    holders = _open_holders(expected, exempt)
    if holders:
        detail = ", ".join(
            f"pid={pid}/fd={fd}({command})" for pid, fd, command in holders
        )
        raise RuntimeError(f"refusing open output file {path}: {detail}")


def _validate_min_age(min_age: float) -> float:
    value = float(min_age)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError("min_age must be finite and non-negative")
    return value


def stable_sha256(
    path: os.PathLike[str] | str, min_age: float = 0
) -> dict[str, Any]:
    """Stream-hash one proven-closed regular file and return its stable identity."""

    minimum_age = _validate_min_age(min_age)
    checked_path = _path(path)
    _, initial = _lstat_regular(checked_path)
    age_seconds = (time.time_ns() - initial["mtime_ns"]) / 1.0e9
    if age_seconds < minimum_age:
        raise RuntimeError(
            f"refusing possibly active file {checked_path}: age {age_seconds:.3f}s "
            f"< {minimum_age:.3f}s"
        )
    _assert_closed(checked_path, initial)

    checked_path, stream, expected = _open_regular_nofollow(checked_path)
    digest = hashlib.sha256()
    try:
        exempt = {(os.getpid(), stream.fileno())}
        _assert_closed(checked_path, expected, exempt)
        for chunk in iter(lambda: stream.read(HASH_CHUNK_BYTES), b""):
            digest.update(chunk)
        _assert_stream_signature(stream, checked_path, expected, "while hashing")
        _assert_path_signature(checked_path, expected, "while hashing")
        _assert_closed(checked_path, expected, exempt)
    finally:
        stream.close()

    _assert_path_signature(checked_path, expected, "after hashing")
    _assert_closed(checked_path, expected)
    return {
        "path": str(checked_path),
        **expected,
        "sha256": digest.hexdigest(),
        "age_seconds": age_seconds,
        "closure_check": "linux_proc_fd",
    }


def install_immutable_json(
    path: os.PathLike[str] | str, payload: Any
) -> dict[str, Any]:
    """Atomically install a new, read-only JSON record without replacement."""

    destination = _path(path)
    parent = destination.parent
    parent.mkdir(parents=True, exist_ok=True)
    parent_stat = parent.lstat()
    if stat.S_ISLNK(parent_stat.st_mode) or not stat.S_ISDIR(parent_stat.st_mode):
        raise ValueError(f"manifest parent is not a real directory: {parent}")
    encoded = (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")

    temporary: Path | None = None
    descriptor: int | None = None
    for _ in range(128):
        candidate = parent / (
            f".{destination.name}.{os.getpid()}.{secrets.token_hex(8)}.tmp"
        )
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        try:
            descriptor = os.open(candidate, flags, 0o600)
        except FileExistsError:
            continue
        temporary = candidate
        break
    if descriptor is None or temporary is None:
        raise RuntimeError(f"could not allocate manifest temporary in {parent}")

    installed = False
    try:
        with os.fdopen(descriptor, "wb", closefd=True) as stream:
            stream.write(encoded)
            stream.flush()
            os.fsync(stream.fileno())
            os.fchmod(stream.fileno(), 0o444)
            os.fsync(stream.fileno())
        try:
            os.link(temporary, destination, follow_symlinks=False)
        except FileExistsError as exc:
            raise FileExistsError(
                f"refusing to replace existing immutable JSON: {destination}"
            ) from exc
        installed = True
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        directory_fd = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)

    if not installed:
        raise RuntimeError(f"failed to install immutable JSON: {destination}")
    stat_result, signature = _lstat_regular(destination)
    if stat.S_IMODE(stat_result.st_mode) & 0o222:
        raise RuntimeError(f"installed JSON is unexpectedly writable: {destination}")
    return {
        "path": str(destination),
        **signature,
        "sha256": hashlib.sha256(encoded).hexdigest(),
    }
