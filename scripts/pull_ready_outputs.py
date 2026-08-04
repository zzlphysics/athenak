#!/usr/bin/env python3
"""Pull checksum-manifested cloud outputs into a writable NAS directory."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shlex
import subprocess
import sys
import time


def run(command: list[str], *, input_text: str | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, input=input_text, text=True, check=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def remote(host: str, command: list[str], *, input_text: str | None = None) -> str:
    remote_command = " ".join(shlex.quote(value) for value in command)
    return run(["ssh", host, remote_command], input_text=input_text).stdout


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_relative(value: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if path.is_absolute() or not path.parts or any(part in ("", ".", "..") for part in path.parts):
        raise ValueError(f"unsafe manifest path: {value!r}")
    return path


def assert_writable_directory(destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    probe = destination / f".athenak-write-test-{os.getpid()}"
    try:
        with probe.open("x", encoding="utf-8") as stream:
            stream.write("ok\n")
        probe.unlink()
    except OSError as exc:
        raise SystemExit(f"destination is not writable: {destination}: {exc}") from exc


def remote_usage_percent(host: str, remote_root: str) -> int:
    output = remote(host, ["df", "-P", remote_root])
    lines = [line for line in output.splitlines() if line.strip()]
    if len(lines) < 2:
        raise RuntimeError(f"unexpected df output: {output!r}")
    return int(lines[-1].split()[4].rstrip("%"))


def list_manifests(host: str, manifest_dir: str) -> list[str]:
    output = remote(host, ["find", manifest_dir, "-maxdepth", "1", "-type", "f",
                            "-name", "*.manifest.ready", "-printf", "%f\\n"])
    return sorted(line for line in output.splitlines() if line)


def load_manifest(host: str, manifest_dir: str, name: str) -> dict:
    if PurePosixPath(name).name != name or not name.endswith(".manifest.ready"):
        raise ValueError(f"unsafe manifest name: {name!r}")
    payload = remote(host, ["cat", f"{manifest_dir.rstrip('/')}/{name}"])
    manifest = json.loads(payload)
    if manifest.get("schema") != 1 or not isinstance(manifest.get("files"), list):
        raise ValueError(f"unsupported manifest schema in {name}")
    return manifest


def already_verified(path: Path, size: int, digest: str) -> bool:
    return path.is_file() and path.stat().st_size == size and sha256(path) == digest


def transfer_manifest(args: argparse.Namespace, name: str) -> None:
    manifest = load_manifest(args.host, args.remote_manifest_dir, name)
    segment = manifest.get("segment")
    if not isinstance(segment, str) or not segment:
        raise ValueError(f"manifest {name} has no segment")

    incoming_root = args.destination / ".incoming" / segment
    ack_root = args.destination / ".acks"
    transferred = []
    for record in manifest["files"]:
        relative = safe_relative(record["path"])
        size = int(record["size"])
        digest = str(record["sha256"])
        if size < 0 or len(digest) != 64:
            raise ValueError(f"invalid record in {name}: {record!r}")

        final_path = args.destination.joinpath(*relative.parts)
        if already_verified(final_path, size, digest):
            transferred.append(record)
            continue
        if final_path.exists():
            raise RuntimeError(
                f"refusing to overwrite existing file with different content: {final_path}")

        partial_path = incoming_root.joinpath(*relative.parts).with_name(relative.name + ".part")
        partial_path.parent.mkdir(parents=True, exist_ok=True)
        source = f"{args.host}:{args.remote_root.rstrip('/')}/{relative.as_posix()}"
        command = ["rsync", "--partial", "--append-verify", "--protect-args",
                   f"--bwlimit={args.bwlimit_kib}", source, str(partial_path)]
        run(command)
        if partial_path.stat().st_size != size:
            raise RuntimeError(f"size mismatch for {relative}")
        if sha256(partial_path) != digest:
            raise RuntimeError(f"SHA256 mismatch for {relative}")
        final_path.parent.mkdir(parents=True, exist_ok=True)
        os.replace(partial_path, final_path)
        transferred.append(record)

    ack_root.mkdir(parents=True, exist_ok=True)
    ack = {
        "schema": 1,
        "manifest": name,
        "segment": segment,
        "verified_unix": time.time(),
        "destination": str(args.destination),
        "files": transferred,
    }
    ack_text = json.dumps(ack, indent=2, sort_keys=True) + "\n"
    local_ack = ack_root / f"{name}.ack"
    temporary = ack_root / f".{name}.{os.getpid()}.tmp"
    temporary.write_text(ack_text, encoding="utf-8")
    os.replace(temporary, local_ack)

    remote_ack_dir = f"{args.remote_manifest_dir.rstrip('/')}/acks"
    remote(args.host, ["mkdir", "-p", remote_ack_dir])
    remote(args.host, ["sh", "-c", f"cat > {shlex.quote(remote_ack_dir + '/' + name + '.ack')}"] ,
           input_text=ack_text)
    print(f"verified {name}: {len(transferred)} files")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", required=True, help="SSH destination, e.g. root@gpu-host")
    parser.add_argument("--remote-root", required=True)
    parser.add_argument("--remote-manifest-dir", required=True)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--bwlimit-kib", type=int, default=9000)
    parser.add_argument("--poll-seconds", type=float, default=0.0,
                        help="0 performs one pass; positive values poll continuously")
    args = parser.parse_args()
    args.destination = args.destination.expanduser().resolve()
    assert_writable_directory(args.destination)

    while True:
        usage = remote_usage_percent(args.host, args.remote_root)
        level = "ok" if usage < 65 else "warning" if usage < 75 else "hold" if usage < 80 else "stop"
        print(f"remote disk: {usage}% ({level})", file=sys.stderr)
        for name in list_manifests(args.host, args.remote_manifest_dir):
            ack = args.destination / ".acks" / f"{name}.ack"
            if not ack.exists():
                transfer_manifest(args, name)
        if args.poll_seconds <= 0:
            return 0 if usage < 75 else 75
        time.sleep(args.poll_seconds)


if __name__ == "__main__":
    raise SystemExit(main())
