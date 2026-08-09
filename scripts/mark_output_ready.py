#!/usr/bin/env python3
"""Atomically publish a checksum manifest for a completed AthenaK output segment."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative_file(root: Path, value: str) -> tuple[Path, Path]:
    candidate = Path(value)
    path = candidate if candidate.is_absolute() else root / candidate
    path = path.resolve(strict=True)
    try:
        relative = path.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"output is outside root: {path}") from exc
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"output is not a regular file: {path}")
    return path, relative


def processes_holding(path: Path) -> list[tuple[int, str]]:
    """Return Linux processes with an open descriptor for *path*.

    Size/mtime stability alone is not enough for production output: a writer may keep a
    large file open while it is temporarily idle.  Comparing device/inode pairs through
    ``/proc`` lets the publisher fail closed before a transfer manifest can expose such a
    file.  On non-Linux systems, where ``/proc`` is unavailable, publication is refused.
    """

    proc = Path("/proc")
    if not proc.is_dir():
        raise RuntimeError("cannot prove file closure because /proc is unavailable")
    target = path.stat()
    holders: list[tuple[int, str]] = []
    for process in proc.iterdir():
        if not process.name.isdigit() or int(process.name) == os.getpid():
            continue
        descriptors = process / "fd"
        try:
            entries = list(descriptors.iterdir())
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        held = False
        for descriptor in entries:
            try:
                descriptor_stat = descriptor.stat()
            except (FileNotFoundError, PermissionError, ProcessLookupError):
                continue
            if (descriptor_stat.st_dev, descriptor_stat.st_ino) == (
                target.st_dev,
                target.st_ino,
            ):
                held = True
                break
        if held:
            try:
                command = (process / "comm").read_text(encoding="utf-8").strip()
            except (FileNotFoundError, PermissionError, ProcessLookupError):
                command = "unknown"
            holders.append((int(process.name), command))
    return sorted(holders)


def assert_closed(path: Path) -> None:
    holders = processes_holding(path)
    if holders:
        detail = ", ".join(f"pid={pid}({command})" for pid, command in holders)
        raise SystemExit(f"refusing open output file {path}: {detail}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manifest-dir", required=True, type=Path)
    parser.add_argument("--segment", required=True)
    parser.add_argument("--min-age", type=float, default=120.0,
                        help="refuse files modified more recently than this many seconds")
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()

    root = args.root.expanduser().resolve(strict=True)
    manifest_dir = args.manifest_dir.expanduser().resolve()
    valid_segment_characters = (
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_."
    )
    if not args.segment or any(
        char not in valid_segment_characters for char in args.segment
    ):
        parser.error("segment may contain only letters, digits, '-', '_' and '.'")

    now = time.time()
    records = []
    seen: set[str] = set()
    for value in args.files:
        path, relative = relative_file(root, value)
        relative_text = relative.as_posix()
        if relative_text in seen:
            continue
        age = now - path.stat().st_mtime
        if age < args.min_age:
            raise SystemExit(
                f"refusing possibly open file {relative_text}: age {age:.1f}s "
                f"< {args.min_age}s"
            )
        assert_closed(path)
        stat_before = path.stat()
        checksum = file_sha256(path)
        stat_after = path.stat()
        if (stat_before.st_size, stat_before.st_mtime_ns) != (
                stat_after.st_size, stat_after.st_mtime_ns):
            raise SystemExit(f"file changed while hashing: {relative_text}")
        assert_closed(path)
        records.append({
            "path": relative_text,
            "size": stat_after.st_size,
            "sha256": checksum,
        })
        seen.add(relative_text)

    manifest_dir.mkdir(parents=True, exist_ok=True)
    destination = manifest_dir / f"{args.segment}.manifest.ready"
    temporary = manifest_dir / f".{args.segment}.{os.getpid()}.tmp"
    payload = {
        "schema": 1,
        "segment": args.segment,
        "created_unix": time.time(),
        "root": str(root),
        "files": records,
    }
    with temporary.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    # A ready manifest is an immutable hand-off record.  Installing it with a hard
    # link is atomic and, unlike os.replace(), cannot silently change the contents
    # behind an ACK that a receiver has already written.  Re-runs must use a new
    # segment identifier (or explicitly inspect/remove a stale unpublished file).
    try:
        os.link(temporary, destination)
    except FileExistsError as exc:
        raise SystemExit(
            f"refusing to replace existing ready manifest: {destination}"
        ) from exc
    finally:
        temporary.unlink(missing_ok=True)
    directory_fd = os.open(manifest_dir, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
