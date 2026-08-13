"""Tests for reusable fail-closed output-integrity primitives."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "output_integrity.py"
SPEC = importlib.util.spec_from_file_location("output_integrity", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
INTEGRITY = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = INTEGRITY
SPEC.loader.exec_module(INTEGRITY)


def test_stable_sha256_reports_full_stable_signature(tmp_path: Path) -> None:
    output = tmp_path / "closed.bin"
    content = b"athenak-closed-output" * 1024
    output.write_bytes(content)

    result = INTEGRITY.stable_sha256(output)
    observed = output.stat()

    assert result["sha256"] == hashlib.sha256(content).hexdigest()
    assert result["size"] == len(content)
    assert result["device"] == observed.st_dev
    assert result["inode"] == observed.st_ino
    assert result["mtime_ns"] == observed.st_mtime_ns
    assert result["ctime_ns"] == observed.st_ctime_ns
    assert result["closure_check"] == "linux_proc_fd"


def test_open_file_and_symbolic_link_are_refused(tmp_path: Path) -> None:
    output = tmp_path / "open.bin"
    output.write_bytes(b"open")
    descriptor = os.open(output, os.O_RDONLY)
    try:
        with pytest.raises(RuntimeError, match="refusing open output file"):
            INTEGRITY.stable_sha256(output)
    finally:
        os.close(descriptor)

    link = tmp_path / "link.bin"
    link.symlink_to(output.name)
    with pytest.raises(ValueError, match="symbolic link"):
        INTEGRITY.stable_sha256(link)


def test_minimum_age_and_missing_proc_are_fail_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "young.bin"
    output.write_bytes(b"young")
    with pytest.raises(RuntimeError, match="possibly active"):
        INTEGRITY.stable_sha256(output, min_age=3600)

    monkeypatch.setattr(INTEGRITY, "PROC_ROOT", tmp_path / "not-proc")
    with pytest.raises(RuntimeError, match="/proc is unavailable"):
        INTEGRITY.stable_sha256(output)


def test_immutable_json_is_atomic_read_only_and_never_replaced(
    tmp_path: Path,
) -> None:
    destination = tmp_path / "records" / "segment.ready.json"
    payload = {"schema": 1, "files": [{"path": "a.rst", "size": 7}]}

    result = INTEGRITY.install_immutable_json(destination, payload)

    assert json.loads(destination.read_text(encoding="utf-8")) == payload
    assert result["sha256"] == hashlib.sha256(destination.read_bytes()).hexdigest()
    assert stat.S_IMODE(destination.stat().st_mode) & 0o222 == 0
    original = destination.read_bytes()
    with pytest.raises(FileExistsError, match="refusing to replace"):
        INTEGRITY.install_immutable_json(destination, {"schema": 2})
    assert destination.read_bytes() == original
    assert not list(destination.parent.glob(".*.tmp"))


def test_json_rejects_nonfinite_before_installation(tmp_path: Path) -> None:
    destination = tmp_path / "bad.json"
    with pytest.raises(ValueError):
        INTEGRITY.install_immutable_json(destination, {"bad": float("nan")})
    assert not destination.exists()


def test_strict_json_rejects_duplicate_keys_and_nonfinite_constants() -> None:
    with pytest.raises(ValueError, match="duplicate JSON key"):
        INTEGRITY.strict_json_loads(b'{"schema":1,"schema":2}', "duplicate.json")
    with pytest.raises(ValueError, match="non-finite JSON number"):
        INTEGRITY.strict_json_loads(b'{"value":NaN}', "nan.json")
