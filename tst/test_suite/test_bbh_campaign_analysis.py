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
    classify_binary,
    classify_history,
    load_verified_files,
    merge_histories,
    storage_projection,
    subcycling_work_model,
    summarize_event_logs,
    verify_file,
)
from compare_bbh_grmhd_convergence import analyze, automatic_diagnostics  # noqa: E402


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


def test_bbh_user_history_selects_native_convergence_diagnostics() -> None:
    time = np.asarray([0.0, 1.0])
    series = []
    for label in ("L", "M"):
        series.append(
            (
                label,
                {
                    "time": time,
                    "baryon_m": np.asarray([10.0, 9.9]),
                    "inner_D": np.asarray([1.0, 1.1]),
                    "emag_prp": np.asarray([0.1, 0.2]),
                    "sigma_D": np.asarray([0.01, 0.02]),
                },
            )
        )
    assert automatic_diagnostics(series) == [
        "baryon_mass_rel",
        "inner_mass_rel",
        "magnetic_energy_rel",
        "D_weighted_sigma",
    ]


def test_storage_budget_includes_forced_segment_restarts() -> None:
    streams = {
        "mhd_w_bcc": [
            {
                "time_M": 100.0,
                "bytes": 1_000_000_000,
                "configured_dt_M": 10.0,
            }
        ]
    }
    projection = storage_projection(
        streams,
        target_time=300.0,
        restart_sizes=[24_000_000_000],
        restart_dt=250.0,
        segment_span=100.0,
        root_dt=4.8,
        root_step_seconds=500.0,
        drain_mib_s=8.0,
    )
    assert projection["remaining_binary_bytes"] == 20_000_000_000
    assert projection["restarts"]["effective_archive_cadence_M"] == 100.0
    assert projection["restarts"]["remaining_restarts"] == 2
    assert projection["remaining_archive_bytes"] == 68_000_000_000
    assert projection["transfer_budget"]["average_drain_has_headroom"] is True


def test_subcycling_work_model_uses_strict_two_to_one_levels() -> None:
    work = subcycling_work_model({0: 8, 1: 4, 2: 2})
    assert work["subcycled_meshblock_updates_per_root_step"] == 24
    assert work["global_finest_dt_meshblock_updates_per_root_step"] == 56
    assert work["global_to_subcycled_update_ratio"] == 56 / 24


def test_native_grmhd_diagnostic_stream_is_classified() -> None:
    variables = ("gr_bsq", "gr_lorentz", "gr_sigma", "gr_beta_inv")
    assert classify_binary("run.mhd_gr_diagnostics.00001.bin", variables) == (
        "mhd_gr_diagnostics"
    )
    assert classify_binary("renamed-output.bin", variables) == "mhd_gr_diagnostics"


def test_generic_and_bbh_histories_are_kept_separate() -> None:
    assert classify_history("run.mhd.hst") == "mhd"
    assert classify_history("run.user.hst") == "user"
    assert classify_history("legacy.hst") == "other"


def test_event_log_merge_prefers_later_restart_segment(tmp_path: Path) -> None:
    header = (
        "# Athena event counter data\n"
        "#  cycle eos_dfloor eos_efloor eos_tfloor eos_vceil eos_fail c2p_it fofc\n"
    )
    first = tmp_path / "first.log"
    second = tmp_path / "second.log"
    first.write_text(
        header + "       5        1        2        3        4        5      6        7\n"
        "      10       10       20       30       40       50     60       70\n",
        encoding="utf-8",
    )
    second.write_text(
        header + "      10        2        3        4        5        6      7        8\n"
        "      20        3        4        5        6        7      8        9\n",
        encoding="utf-8",
    )
    summary = summarize_event_logs([first, second])
    assert summary["records_with_events"] == 3
    assert summary["cycle_min"] == 5
    assert summary["cycle_max"] == 20
    assert summary["c2p_iteration_max"] == 8
    assert summary["totals"] == {
        "eos_dfloor": 6,
        "eos_efloor": 9,
        "eos_tfloor": 12,
        "eos_vceil": 15,
        "eos_fail": 18,
        "fofc": 24,
    }
