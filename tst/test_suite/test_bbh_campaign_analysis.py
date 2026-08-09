"""Unit tests for verified BBH campaign analysis and tier alignment."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "vis" / "python"))

from analyze_bbh_grmhd_campaign import (  # noqa: E402
    load_verified_files,
    merge_histories,
    verify_file,
)
from compare_bbh_grmhd_convergence import analyze  # noqa: E402


def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def write_history(path: Path, rows: list[tuple[float, float]]) -> None:
    text = "# Athena++ history data\n# [0]=time [1]=mass\n"
    text += "".join(f"{time:.17e} {mass:.17e}\n" for time, mass in rows)
    path.write_text(text, encoding="utf-8")


def test_verified_file_change_is_rejected(tmp_path: Path) -> None:
    output = tmp_path / "state" / "run.hst"
    output.parent.mkdir()
    payload = b"closed-history\n"
    output.write_bytes(payload)
    ack_dir = tmp_path / ".acks"
    ack_dir.mkdir()
    manifest_name = "segment.manifest.ready"
    file_record = {
        "path": "state/run.hst",
        "size": len(payload),
        "sha256": digest(payload),
    }
    manifest_bytes = (
        json.dumps(
            {
                "schema": 1,
                "segment": "segment",
                "files": [file_record],
            }
        )
        + "\n"
    ).encode()
    manifest_dir = tmp_path / ".manifests"
    manifest_dir.mkdir()
    (manifest_dir / manifest_name).write_bytes(manifest_bytes)
    (ack_dir / "segment.manifest.ready.ack").write_text(
        json.dumps(
            {
                "schema": 1,
                "manifest": manifest_name,
                "manifest_sha256": digest(manifest_bytes),
                "segment": "segment",
                "files": [file_record],
            }
        ),
        encoding="utf-8",
    )
    record = load_verified_files(tmp_path)[0]
    verify_file(record, True)
    output.write_bytes(b"changed-history\n")
    try:
        verify_file(record, True)
    except RuntimeError as exc:
        assert "changed" in str(exc)
    else:
        raise AssertionError("changed ACK output was accepted")


def test_restart_history_merge_prefers_later_segment(tmp_path: Path) -> None:
    first = tmp_path / "first.hst"
    second = tmp_path / "second.hst"
    write_history(first, [(0.0, 1.0), (1.0, 2.0)])
    write_history(second, [(1.0, 20.0), (2.0, 30.0)])
    columns, merged = merge_histories([first, second])
    assert columns == ["time", "mass"]
    np.testing.assert_allclose(merged["time"], [0.0, 1.0, 2.0])
    np.testing.assert_allclose(merged["mass"], [1.0, 20.0, 30.0])


def test_empirical_second_order_difference_ratio() -> None:
    time = np.linspace(0.0, 10.0, 21)
    shape = np.sin(time / 10.0)
    truth = 0.02 * time
    series = []
    for label, error in (("L", 4.0e-3), ("M", 1.0e-3), ("H", 2.5e-4)):
        mass = 1.0 + truth + error * shape
        series.append((label, {"time": time, "mass": mass}))
    _, summaries, _ = analyze(series, ["mass_rel"], 2.0)
    assert abs(float(summaries["mass_rel"]["empirical_order"]) - 2.0) < 1.0e-12
