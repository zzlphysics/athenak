"""Unit tests for verified BBH campaign analysis and tier alignment."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
import struct
import sys
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "vis" / "python"))

from analyze_bbh_grmhd_campaign import (  # noqa: E402
    classify_binary,
    classify_history,
    configured_output_dcycle,
    configured_output_dt,
    configured_output_region,
    is_athenak_binary_output,
    load_verified_files,
    merge_histories,
    meshblock_level_counts,
    storage_projection,
    subcycling_work_model,
    summarize_event_logs,
    verify_file,
)
from compare_bbh_grmhd_convergence import analyze, automatic_diagnostics  # noqa: E402
from plot_bbh_grmhd import (  # noqa: E402
    PANELS,
    SliceBlock,
    SliceData,
    panel_values,
    read_cell_face_interpolated_slice,
)


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


def test_storage_budget_uses_root_cycle_cadence_for_local_slices() -> None:
    streams = {
        "bbh_local_w": [
            {
                "time_M": 0.0,
                "bytes": 10_000_000,
                "configured_dt_M": None,
                "configured_root_dcycle": 1,
            }
        ]
    }
    projection = storage_projection(
        streams,
        target_time=9.6,
        restart_sizes=[],
        restart_dt=None,
        segment_span=None,
        root_dt=4.8,
        root_step_seconds=None,
        drain_mib_s=8.0,
    )
    local = projection["streams"]["bbh_local_w"]
    assert local["cadence_M"] == 4.8
    assert local["remaining_frames"] == 2
    assert local["remaining_bytes"] == 20_000_000


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
    assert classify_binary(
        "renamed-masked-output.bin", variables + ("gr_excision_mask",)
    ) == "mhd_gr_diagnostics"
    assert classify_binary("run.bbh_local_gr.00001.bin", variables) == "bbh_local_gr"
    assert classify_binary(
        "run.bbh_local_w.00001.bin", ("dens", "velx", "bcc1")
    ) == "bbh_local_w"


def test_athenak_binary_signature_rejects_unrelated_bin(tmp_path: Path) -> None:
    dump = tmp_path / "output.bin"
    dump.write_bytes(b"Athena binary output version=1.1\nheader\n")
    executable = tmp_path / "test_mpi_CXX.bin"
    executable.write_bytes(b"\x7fELF\x02\x01\x01\x00not-an-athenak-dump")

    assert is_athenak_binary_output(dump)
    assert not is_athenak_binary_output(executable)


def test_native_gr_panel_automatically_masks_excised_cells() -> None:
    block = SliceBlock(
        extent=(-1.0, 1.0, -1.0, 1.0),
        level=0,
        logical_location=(0, 0, 0),
        slice_shape=(1, 2),
        fields={
            "gr_sigma": np.asarray([[2.0, 1.0e8]]),
            "gr_excision_mask": np.asarray([[0.0, 1.0]]),
        },
    )
    values = panel_values(
        SimpleNamespace(blocks=[block]), PANELS["gr_sigma"], density_threshold=None
    )
    np.testing.assert_array_equal(np.isfinite(values[0]), [[True, False]])
    assert values[0][0, 0] == 2.0

    legacy_block = SliceBlock(
        extent=block.extent,
        level=block.level,
        logical_location=block.logical_location,
        slice_shape=block.slice_shape,
        fields={"gr_sigma": np.asarray([[2.0, 1.0e8]])},
    )
    legacy_values = panel_values(
        SimpleNamespace(blocks=[legacy_block]),
        PANELS["gr_sigma"],
        density_threshold=None,
    )
    assert np.isfinite(legacy_values[0]).all()


def test_cell_face_slice_averages_both_cell_centers() -> None:
    header = SimpleNamespace(
        parameters={"mesh": {"x3min": "-1", "x3max": "1", "nx3": "2"}}
    )
    common = {
        "extent": (-1.0, 1.0, -1.0, 1.0),
        "level": 0,
        "slice_shape": (1, 2),
    }
    lower = SliceData(
        header=header,
        plane="z",
        location=-1.0e-12,
        blocks=[
            SliceBlock(
                logical_location=(0, 0, 0),
                fields={"dens": np.asarray([[1.0, 3.0]])},
                **common,
            )
        ],
        level_counts=Counter({0: 2}),
        selected_level_counts=Counter({0: 1}),
        presliced=False,
    )
    upper = SliceData(
        header=header,
        plane="z",
        location=1.0e-12,
        blocks=[
            SliceBlock(
                logical_location=(0, 0, 1),
                fields={"dens": np.asarray([[5.0, 7.0]])},
                **common,
            )
        ],
        level_counts=lower.level_counts,
        selected_level_counts=lower.selected_level_counts,
        presliced=False,
    )
    with patch("plot_bbh_grmhd.read_binary_header", return_value=header), patch(
        "plot_bbh_grmhd.read_slice", side_effect=[lower, upper]
    ):
        interpolated = read_cell_face_interpolated_slice(
            Path("full-3d.bin"), ["dens"], "z", 0.0, 1.0
        )
    np.testing.assert_allclose(interpolated.blocks[0].fields["dens"], [[3.0, 5.0]])
    assert interpolated.location == 0.0
    assert not interpolated.presliced


def test_local_and_global_duplicate_variables_use_ids_for_cadence() -> None:
    header = SimpleNamespace(
        parameters={
            "output2": {
                "file_type": "bin",
                "variable": "mhd_w_bcc",
                "id": "mhd_w_bcc",
                "dt": "50",
            },
            "output7": {
                "file_type": "bin",
                "variable": "mhd_w_bcc",
                "id": "bbh_local_w",
                "dcycle": "1",
                "region_center": "bbh_com",
                "region_half_width": "40",
                "region_slice_axis": "3",
            },
        }
    )
    assert configured_output_dt(header, "mhd_w_bcc") == 50.0
    assert configured_output_dt(header, "bbh_local_w") is None
    assert configured_output_dcycle(header, "bbh_local_w") == 1
    assert configured_output_region(header, "bbh_local_w") == {
        "center": "bbh_com",
        "half_width_M": [40.0, 40.0, 40.0],
        "slice_axis": 3,
        "slice_offset_M": 0.0,
    }


def test_meshblock_scan_accepts_full_and_sliced_record_shapes(tmp_path: Path) -> None:
    path = tmp_path / "mixed-shape-records.bin"

    def record(shape: tuple[int, int, int], level: int) -> bytes:
        nx1, nx2, nx3 = shape
        indices = struct.pack("=6i", 0, nx1 - 1, 0, nx2 - 1, 0, nx3 - 1)
        logical = struct.pack("=4i", 0, 0, 0, level)
        geometry = struct.pack("=6d", 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
        data = bytes(2 * nx1 * nx2 * nx3 * 4)
        return indices + logical + geometry + data

    path.write_bytes(record((2, 2, 2), 0) + record((2, 2, 1), 1))
    header = SimpleNamespace(
        data_offset=0,
        variables=("a", "b"),
        variable_size=4,
        location_size=8,
    )
    assert meshblock_level_counts(path, header) == {0: 1, 1: 1}


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
    assert summary["records"] == 3
    assert summary["records_with_events"] == 3
    assert summary["records_without_corrective_events"] == 0
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
        "cons_adjust": 0,
        "mag_adjust": 0,
        "c2p_calls": 0,
        "fofc_tests": 0,
    }


def test_extended_event_log_uses_header_names_and_merges_with_legacy(
    tmp_path: Path,
) -> None:
    legacy = tmp_path / "legacy.log"
    extended = tmp_path / "extended.log"
    legacy.write_text(
        "# Athena event counter data\n"
        "# cycle eos_dfloor eos_efloor eos_tfloor eos_vceil eos_fail c2p_it fofc\n"
        "10 1 2 3 4 5 6 7\n",
        encoding="utf-8",
    )
    # Deliberately reorder both legacy and appended fields to make sure values
    # are associated through the header rather than fixed column positions.
    extended.write_text(
        "# Athena event counter data\n"
        "# cycle c2p_calls eos_fail eos_dfloor mag_adjust c2p_it fofc_tests "
        "eos_tfloor cons_adjust fofc eos_vceil eos_efloor\n"
        "10 100 50 10 90 80 110 30 70 60 40 20\n"
        "20 101 51 11 91 81 111 31 71 61 41 21\n",
        encoding="utf-8",
    )

    summary = summarize_event_logs([legacy, extended])

    assert summary["records"] == 2
    assert summary["records_with_events"] == 2
    assert summary["records_without_corrective_events"] == 0
    assert summary["cycle_min"] == 10
    assert summary["cycle_max"] == 20
    assert summary["c2p_iteration_max"] == 81
    assert summary["totals"] == {
        "eos_dfloor": 21,
        "eos_efloor": 41,
        "eos_tfloor": 61,
        "eos_vceil": 81,
        "eos_fail": 101,
        "fofc": 121,
        "cons_adjust": 141,
        "mag_adjust": 181,
        "c2p_calls": 201,
        "fofc_tests": 221,
    }


def test_dense_event_log_distinguishes_observation_from_correction(
    tmp_path: Path,
) -> None:
    event_log = tmp_path / "dense.log"
    event_log.write_text(
        "# Athena event counter data\n"
        "# cycle eos_dfloor eos_efloor eos_tfloor eos_vceil eos_fail c2p_it "
        "fofc cons_adjust mag_adjust c2p_calls fofc_tests\n"
        "1 0 0 0 0 0 3 0 0 0 100 80\n"
        "2 1 0 0 0 0 3 0 1 0 100 80\n",
        encoding="utf-8",
    )

    summary = summarize_event_logs([event_log])

    assert summary["records"] == 2
    assert summary["records_with_events"] == 1
    assert summary["records_without_corrective_events"] == 1
