"""Tests for fail-closed cloud-output publication."""

from __future__ import annotations

import importlib.util
import json
import mmap
import os
from pathlib import Path
import stat
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
PUBLISHER = ROOT / "scripts" / "mark_output_ready.py"


def load_publisher():
    spec = importlib.util.spec_from_file_location("mark_output_ready", PUBLISHER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def command(root: Path, output: Path, manifest_dir: Path) -> list[str]:
    return [
        sys.executable,
        str(PUBLISHER),
        "--root",
        str(root),
        "--manifest-dir",
        str(manifest_dir),
        "--segment",
        "segment-1",
        "--min-age",
        "0",
        str(output.relative_to(root)),
    ]


def test_closed_file_is_published(tmp_path: Path) -> None:
    output = tmp_path / "state" / "closed.bin"
    output.parent.mkdir()
    output.write_bytes(b"closed-output")
    result = subprocess.run(
        command(tmp_path, output, tmp_path / "manifests"),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    manifest = tmp_path / "manifests" / "segment-1.manifest.ready"
    assert manifest.is_file()
    assert stat.S_IMODE(manifest.stat().st_mode) == 0o444
    assert manifest.stat().st_nlink == 1
    record = json.loads(manifest.read_text(encoding="utf-8"))["files"][0]
    assert record["path"] == "state/closed.bin"
    assert record["size"] == len(b"closed-output")


def test_open_file_is_refused(tmp_path: Path) -> None:
    output = tmp_path / "state" / "open.bin"
    output.parent.mkdir()
    descriptor = os.open(output, os.O_WRONLY | os.O_CREAT, 0o600)
    try:
        os.write(descriptor, b"still-open")
        result = subprocess.run(
            command(tmp_path, output, tmp_path / "manifests"),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    finally:
        os.close(descriptor)
    assert result.returncode != 0
    assert "refusing open output file" in result.stderr
    assert not (tmp_path / "manifests" / "segment-1.manifest.ready").exists()


def test_ready_manifest_is_immutable(tmp_path: Path) -> None:
    output = tmp_path / "state" / "closed.bin"
    output.parent.mkdir()
    output.write_bytes(b"first-content")
    first = subprocess.run(
        command(tmp_path, output, tmp_path / "manifests"),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    assert first.returncode == 0, first.stderr
    manifest = tmp_path / "manifests" / "segment-1.manifest.ready"
    original = manifest.read_bytes()

    output.write_bytes(b"different-content")
    second = subprocess.run(
        command(tmp_path, output, tmp_path / "manifests"),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    assert second.returncode != 0
    assert "refusing to replace existing ready manifest" in second.stderr
    assert manifest.read_bytes() == original


def test_symbolic_link_output_is_refused(tmp_path: Path) -> None:
    state = tmp_path / "state"
    state.mkdir()
    target = state / "target.bin"
    target.write_bytes(b"target")
    output = state / "link.bin"
    output.symlink_to(target.name)

    result = subprocess.run(
        command(tmp_path, output, tmp_path / "manifests"),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    assert result.returncode != 0
    assert "not a regular file" in result.stderr
    assert not (tmp_path / "manifests" / "segment-1.manifest.ready").exists()


def test_symbolic_link_parent_is_refused(tmp_path: Path) -> None:
    real_state = tmp_path / "real-state"
    real_state.mkdir()
    output = real_state / "output.bin"
    output.write_bytes(b"target")
    (tmp_path / "state").symlink_to(real_state.name, target_is_directory=True)

    result = subprocess.run(
        command(tmp_path, tmp_path / "state" / output.name, tmp_path / "manifests"),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    assert result.returncode != 0
    assert not (tmp_path / "manifests" / "segment-1.manifest.ready").exists()


def test_hard_link_output_is_refused(tmp_path: Path) -> None:
    state = tmp_path / "state"
    state.mkdir()
    output = state / "output.bin"
    output.write_bytes(b"target")
    os.link(output, state / "second-name.bin")

    result = subprocess.run(
        command(tmp_path, output, tmp_path / "manifests"),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    assert result.returncode != 0
    assert "exactly one hard link" in result.stderr
    assert not (tmp_path / "manifests" / "segment-1.manifest.ready").exists()


def test_path_replacement_during_hash_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    publisher = load_publisher()
    output = tmp_path / "state" / "output.bin"
    output.parent.mkdir()
    output.write_bytes(b"original-bytes")
    moved = output.with_name("original-moved.bin")
    original_hash = publisher.file_sha256_fd

    def replace_after_hash(descriptor: int) -> str:
        digest = original_hash(descriptor)
        output.rename(moved)
        output.write_bytes(b"replacement-ok")
        return digest

    monkeypatch.setattr(publisher, "file_sha256_fd", replace_after_hash)
    monkeypatch.setattr(
        sys,
        "argv",
        command(tmp_path, output, tmp_path / "manifests")[1:],
    )
    with pytest.raises(SystemExit, match="changed while hashing"):
        publisher.main()
    assert not (tmp_path / "manifests" / "segment-1.manifest.ready").exists()


def test_writable_mmap_is_refused_after_descriptor_close(tmp_path: Path) -> None:
    output = tmp_path / "state" / "mapped.bin"
    output.parent.mkdir()
    output.write_bytes(b"mapped-output")
    descriptor = os.open(output, os.O_RDWR)
    mapping = mmap.mmap(descriptor, 0, access=mmap.ACCESS_WRITE)
    os.close(descriptor)
    try:
        result = subprocess.run(
            command(tmp_path, output, tmp_path / "manifests"),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    finally:
        mapping.close()
    assert result.returncode != 0
    assert "writable-mmap" in result.stderr
    assert not (tmp_path / "manifests" / "segment-1.manifest.ready").exists()


def test_proc_permission_gap_fails_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    publisher = load_publisher()
    fake_proc = tmp_path / "proc"
    denied_fd = fake_proc / "987654" / "fd"
    denied_fd.mkdir(parents=True)
    (fake_proc / "987654" / "maps").write_text("", encoding="utf-8")
    original_iterdir = Path.iterdir

    def deny_descriptor_listing(path: Path):
        if path == denied_fd:
            raise PermissionError("test denial")
        return original_iterdir(path)

    monkeypatch.setattr(Path, "iterdir", deny_descriptor_listing)
    output = tmp_path / "output.bin"
    output.write_bytes(b"closed")
    descriptor = os.open(output, os.O_RDONLY)
    try:
        with pytest.raises(RuntimeError, match="permission denied inspecting"):
            publisher.processes_holding(os.fstat(descriptor), fake_proc)
    finally:
        os.close(descriptor)
