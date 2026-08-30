"""Tests for restartable frozen direct-EMRI production planning."""

from __future__ import annotations

import math
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import prepare_frozen_direct_pilot as pilot_module  # noqa: E402
import prepare_frozen_production_campaign as production  # noqa: E402


def _pilot() -> dict[str, object]:
    overrides = [
        "job/basename=emri_frozen_direct_calibration",
        "time/nlim=15",
        "time/tlim=1e30",
        "time/subcycling=level",
        "time/root_dt_max=100",
        "adm/dynamic=true",
        "mesh_refinement/refinement=adaptive",
        "mesh_refinement/num_levels=10",
        "mesh_refinement/max_nmb_per_rank=3500",
        "mesh_refinement/ncycle_check=1",
        "mesh_refinement/refinement_interval=1",
        "mesh/nx1=128", "mesh/nx2=112", "mesh/nx3=80",
        "mesh/x1min=-2890.2711304431523",
        "mesh/x1max=1845.119124744866",
        "mesh/x2min=-2645.0832268188583",
        "mesh/x2max=1845.1185514675913",
        "mesh/x3min=-1318.1407462145046",
        "mesh/x3max=1320.1180539362786",
        "meshblock/nx1=16", "meshblock/nx2=16", "meshblock/nx3=16",
        "mesh/ix1_bc=outflow", "mesh/ox1_bc=user",
        "mesh/ix2_bc=outflow", "mesh/ox2_bc=user",
        "mesh/ix3_bc=user", "mesh/ox3_bc=outflow",
        "problem/user_hist=false",
        "problem/background_mode=full",
        "problem/primary_mass=99999.999999999985",
        "problem/secondary_mass=1",
        "problem/primary_chi=0.9375",
        "problem/orbital_radius=5600000",
        "problem/orbit_direction=1",
        "mhd/gamma=1.4444444444444444",
        "output1/dt=0", "output2/dt=0", "output3/dt=0",
        "output4/dt=0", "output5/dt=0",
    ]
    return {
        "classification": pilot_module.CLASSIFICATION,
        "case_id": "r56_p000_t5000",
        "source_state_sha256": "0" * 64,
        "athena_input": "inputs/emri/emri_windtunnel_smoke.athinput",
        "athena_overrides": overrides,
        "mesh": {
            "lower": [-2890.2711304431523, -2645.0832268188583,
                      -1318.1407462145046],
            "upper": [1845.119124744866, 1845.1185514675913,
                      1320.1180539362786],
            "root_dimensions": [128, 112, 80],
            "meshblock_dimensions": [16, 16, 16],
            "physical_refinement_levels": 9,
            "maximum_meshblocks_per_rank": 3500,
        },
        "bhl_plan": {
            "scales": {
                "capture_radius_factor_two_for_cost": 329.0412309432926,
                "capture_crossing_time": 4220.466758165481,
                "direct_settling_duration": 8440.933516330962,
            },
            "costs": {"direct": {
                "nested_refinement_levels": 9,
                "finest_steps": 350251,
                "ideal_level_subcycled_zone_updates": 802607171520,
            }},
        },
    }


def _calibration() -> dict[str, object]:
    return {
        "classification": production.CALIBRATION_CLASSIFICATION,
        "date": "2026-08-30",
        "source": {
            "solver_commit": "5" * 40,
            "campaign_case": "r56_p000_t5000",
            "source_state_sha256": "0" * 64,
        },
        "physical_case": {
            "root_dimensions": [128, 112, 80],
            "meshblock_dimensions": [16, 16, 16],
            "physical_refinement_levels": 9,
            "configured_maximum_meshblocks_per_rank": 3500,
        },
        "build": {"athena_sha256": "a" * 64},
        "result": {
            "passed": True,
            "final_cycle": 15,
            "final_time_in_secondary_masses": 75.51375,
            "steady_nine_level_root_cycle_seconds": 109.1235,
            "zone_cycles_per_second_reported": 10976700.0,
            "final_meshblocks": 2779,
            "peak_gpu_memory_MiB": 47621.0,
        },
    }


def _qualification() -> dict[str, object]:
    return {
        "classification": production.QUALIFICATION_CLASSIFICATION,
        "date": "2026-08-30",
        "case_id": "r56_p000_t5000",
        "source": {"commit": "1" * 40},
        "build": {"athena_sha256": "b" * 64},
        "result": {
            "passed": True,
            "final_meshblocks": 2779,
            "peak_gpu_memory_MiB": 47675.0,
            "final_restart_size_bytes": 7721877209,
            "terminal_output_bundle_size_bytes": 8177706197,
            "post_process_finalize_seconds_proxy": 0.7362776145935186,
        },
        "qualification_analysis": {
            "paired_full_topology_force_diagnostic": {
                "without_force_history_seconds": 108.965,
                "with_force_history_seconds": 119.5162,
            },
        },
        "downloaded_evidence": {"summary_sha256": "c" * 64},
    }


def _restart_qualification() -> dict[str, object]:
    return {
        "classification": production.RESTART_QUALIFICATION_CLASSIFICATION,
        "date": "2026-08-30",
        "case_id": "r56_p000_t5000",
        "source": {"commit": "2" * 40},
        "build": {"athena_sha256": "d" * 64},
        "result": {
            "passed": True,
            "source_checkpoint": {"meshblocks": 2779},
            "cold_restart": {
                "restart_load_tree_allocate_and_cache_rebuild_seconds": 23.90582,
            },
            "durable_sync_seconds": {
                "fresh_cycle10_outputs": 0.994044303894043,
                "resumed_cycle11_outputs": 0.7929956912994385,
                "continuous_cycle11_outputs": 0.8776090145111084,
            },
            "endpoint_comparison": {"all_stored_fields_match": True},
        },
        "qualification_analysis": {
            "production_budget": {
                "conservative_root_cycle_seconds": 120.1823,
                "recommended_operational_reserve_hours": 60.0,
            },
        },
        "downloaded_evidence": {"summary_sha256": "e" * 64},
    }


def test_real_calibration_produces_conservative_segment_and_output_policy() -> None:
    campaign = production.build_production_campaign(_pilot(), _calibration())

    runtime = campaign["runtime_projection"]
    assert runtime["uncorrected_ideal_hours"] == pytest.approx(20.3108800034)
    assert runtime["measured_to_planner_step_count_correction"] == pytest.approx(
        2.45145338629
    )
    assert runtime["measured_to_planner_topology_correction"] == pytest.approx(
        1.02093118259
    )
    assert runtime["empirical_steady_topology_hours"] == pytest.approx(
        50.83336375
    )
    assert runtime["empirical_with_production_diagnostics_hours"] \
        == pytest.approx(50.83336375)
    assert runtime["force_diagnostic_overhead_fraction"] is None
    assert runtime["empirical_steady_topology_hours"] \
        > runtime["cfl_corrected_lower_bound_hours"] \
        > runtime["uncorrected_ideal_hours"]

    segmentation = campaign["segmentation"]
    assert segmentation["root_steps_per_segment"] == 98
    assert segmentation["checkpoint_root_steps"] == 49
    assert segmentation["segments_if_root_dt_is_stationary"] == 18
    assert segmentation["wall_hours_per_segment_proxy"] == pytest.approx(2.97058417)

    scales = campaign["physical_scales"]
    assert scales["force_outer_radii_in_secondary_masses"] == pytest.approx(
        [164.5206154716463, 329.0412309432926, 658.0824618865852]
    )
    assert scales["force_outer_radii_in_secondary_masses"][-1] \
        < scales["conservative_domain_half_width_in_secondary_masses"]

    resources = campaign["resource_envelope"]
    assert resources["meshblock_capacity_headroom"] == 721
    assert resources["estimated_restart_GiB_at_calibrated_topology"] \
        == pytest.approx(7.19148302078)
    assert resources["minimum_working_disk_GiB"] >= 40

    overrides = campaign["fresh_overrides"]
    keys = [item.split("=", 1)[0] for item in overrides]
    assert len(keys) == len(set(keys))
    assert "problem/user_hist=true" in overrides
    assert "output1/variable=mhd_w_bcc" in overrides
    assert "output2/variable=mhd_divb" in overrides
    assert "output4/user_hist_only=true" in overrides
    assert any(item.startswith("problem/force_outer_radius_3=658.082")
               for item in overrides)
    assert campaign["stationarity_gate"]["status"] == "runtime_evidence_required"


def test_passed_production_qualification_controls_runtime_and_resources() -> None:
    campaign = production.build_production_campaign(
        _pilot(), _calibration(), qualification=_qualification()
    )

    runtime = campaign["runtime_projection"]
    assert runtime["budget_root_cycle_seconds"] == pytest.approx(119.5162)
    assert runtime["force_diagnostic_overhead_fraction"] == pytest.approx(
        0.09683109255
    )
    assert runtime["empirical_with_production_diagnostics_hours"] \
        == pytest.approx(55.6746298333)

    segmentation = campaign["segmentation"]
    assert segmentation["root_steps_per_segment"] == 90
    assert segmentation["checkpoint_root_steps"] == 45
    assert segmentation["segments_if_root_dt_is_stationary"] == 19
    assert segmentation["wall_hours_per_segment_proxy"] == pytest.approx(2.987905)

    resources = campaign["resource_envelope"]
    assert resources["qualification_passed"] is True
    assert resources["qualified_peak_gpu_memory_MiB"] == 47675.0
    assert resources["measured_restart_GiB_at_calibrated_topology"] \
        == pytest.approx(7.19155856315)
    assert resources["restart_size_estimate_relative_error"] \
        == pytest.approx(-1.05043110258e-5)
    assert campaign["source_qualification"]["evidence_summary_sha256"] == "c" * 64


def test_restart_qualification_closes_read_budget_and_shortens_segment() -> None:
    campaign = production.build_production_campaign(
        _pilot(),
        _calibration(),
        qualification=_qualification(),
        restart_qualification=_restart_qualification(),
    )

    runtime = campaign["runtime_projection"]
    assert runtime["budget_root_cycle_seconds"] == pytest.approx(120.1823)
    assert runtime["empirical_with_production_diagnostics_hours"] \
        == pytest.approx(55.9849214167)
    assert runtime["projected_restart_read_hours"] == pytest.approx(0.1195291)
    assert runtime["projected_durable_sync_hours"] \
        == pytest.approx(0.00524634494)
    assert runtime["qualified_nominal_hours"] == pytest.approx(56.1096968616)
    assert runtime["recommended_operational_reserve_hours"] == 60.0

    segmentation = campaign["segmentation"]
    assert segmentation["root_steps_per_segment"] == 89
    assert segmentation["checkpoint_root_steps"] == 45
    assert segmentation["segments_if_root_dt_is_stationary"] == 19
    assert segmentation["projected_cold_restart_reads"] == 18
    assert segmentation["wall_hours_per_segment_proxy"] \
        == pytest.approx(2.97117352778)

    resources = campaign["resource_envelope"]
    assert resources["restart_read_qualification_passed"] is True
    assert resources["qualified_cold_restart_seconds"] == pytest.approx(23.90582)
    assert resources["qualified_exact_endpoint_resume"] is True
    assert campaign["source_restart_qualification"][
        "evidence_summary_sha256"
    ] == "e" * 64


def test_identity_mismatch_and_unsafe_force_shell_fail_closed() -> None:
    calibration = _calibration()
    calibration["source"]["source_state_sha256"] = "1" * 64
    with pytest.raises(ValueError, match="source state hash"):
        production.build_production_campaign(_pilot(), calibration)

    pilot = _pilot()
    pilot["mesh"]["lower"] = [-500.0, -500.0, -500.0]
    pilot["mesh"]["upper"] = [500.0, 500.0, 500.0]
    with pytest.raises(ValueError, match="outer force shell"):
        production.build_production_campaign(pilot, _calibration())

    qualification = _qualification()
    qualification["case_id"] = "another-case"
    with pytest.raises(ValueError, match="qualification disagree on case id"):
        production.build_production_campaign(
            _pilot(), _calibration(), qualification=qualification
        )

    with pytest.raises(ValueError, match="requires its production I/O qualification"):
        production.build_production_campaign(
            _pilot(),
            _calibration(),
            restart_qualification=_restart_qualification(),
        )


def test_restart_layout_estimate_includes_faces_and_prescribed_adm() -> None:
    stored = 16 + 2 * 4
    expected_values = (5 + 17) * stored**3 + 3 * (stored + 1) * stored**2
    assert production._restart_bytes_per_meshblock(16) == 8 * expected_values
    assert math.isclose(
        production._restart_bytes_per_meshblock(16) / 2**20,
        2.64990234375,
    )
