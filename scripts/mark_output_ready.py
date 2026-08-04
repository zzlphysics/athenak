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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manifest-dir", required=True, type=Path)
    parser.add_argument("--segment", required=True)
    parser.add_argument("--min-age", type=float, default=30.0,
                        help="refuse files modified more recently than this many seconds")
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()

    root = args.root.expanduser().resolve(strict=True)
    manifest_dir = args.manifest_dir.expanduser().resolve()
    if not args.segment or any(char not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_."
                               for char in args.segment):
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
                f"refusing possibly open file {relative_text}: age {age:.1f}s < {args.min_age}s")
        stat_before = path.stat()
        checksum = file_sha256(path)
        stat_after = path.stat()
        if (stat_before.st_size, stat_before.st_mtime_ns) != (
                stat_after.st_size, stat_after.st_mtime_ns):
            raise SystemExit(f"file changed while hashing: {relative_text}")
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
    os.replace(temporary, destination)
    print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
