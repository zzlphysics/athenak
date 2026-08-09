"""Tests for fail-closed cloud-output publication."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
PUBLISHER = ROOT / "scripts" / "mark_output_ready.py"


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
    assert (tmp_path / "manifests" / "segment-1.manifest.ready").is_file()


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
