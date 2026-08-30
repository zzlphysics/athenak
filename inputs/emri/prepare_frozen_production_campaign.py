#!/usr/bin/env python3
"""Turn a passed frozen direct pilot and GPU calibration into a run policy.

The resulting JSON is deliberately an execution manifest rather than a launcher.  It
contains the fresh-run command, the restart-only overrides for later segments, resource
projections, output cadences, and the force-stationarity gate that must be satisfied
before a history interval is treated as a science average.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re

import prepare_frozen_direct_pilot as direct_pilot


CLASSIFICATION = "athenak-emri-frozen-direct-production-campaign-v1"
CALIBRATION_CLASSIFICATION = "athenak-emri-frozen-direct-cloud-calibration-v1"
REAL_BYTES = 8
MHD_CONSERVED_VARIABLES = 5
PRESCRIBED_ADM_VARIABLES = 17
MHD_PRIMITIVE_OUTPUT_VARIABLES = 8


def _positive(value: float, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{label} must be finite and positive")
    return number


def _positive_integer(value: int, label: str) -> int:
    number = int(value)
    if number < 1:
        raise ValueError(f"{label} must be positive")
    return number


def _replace_overrides(
    source: list[str], replacements: dict[str, str]
) -> list[str]:
    """Replace an exact AthenaK block/name override without changing its order."""

    result: list[str] = []
    seen: set[str] = set()
    for item in source:
        if "=" not in item or "/" not in item.split("=", 1)[0]:
            raise ValueError(f"invalid AthenaK override: {item!r}")
        key = item.split("=", 1)[0]
        if key in seen:
            raise ValueError(f"duplicate AthenaK override: {key}")
        seen.add(key)
        result.append(f"{key}={replacements[key]}" if key in replacements else item)
    missing = set(replacements).difference(seen)
    if missing:
        raise ValueError(
            "pilot lacks required override(s): " + ", ".join(sorted(missing))
        )
    return result


def _restart_bytes_per_meshblock(meshblock_cells: int, nghost: int = 4) -> int:
    stored = meshblock_cells + 2 * nghost
    cell_centered = (
        MHD_CONSERVED_VARIABLES + PRESCRIBED_ADM_VARIABLES
    ) * stored**3
    face_centered = 3 * (stored + 1) * stored**2
    return REAL_BYTES * (cell_centered + face_centered)


def _safe_basename(case_id: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_-]+", "_", case_id).strip("_")
    if not cleaned:
        raise ValueError("case id cannot form a safe AthenaK basename")
    return f"emri_frozen_direct_{cleaned}_full_gradient"


def _validate_identity(
    pilot: dict[str, object], calibration: dict[str, object]
) -> None:
    if pilot.get("classification") != direct_pilot.CLASSIFICATION:
        raise ValueError("input pilot is not a frozen direct pilot")
    if calibration.get("classification") != CALIBRATION_CLASSIFICATION:
        raise ValueError("input calibration is not a frozen direct cloud calibration")
    result = calibration.get("result")
    if not isinstance(result, dict) or result.get("passed") is not True:
        raise ValueError("GPU calibration did not pass")
    source = calibration.get("source")
    physical = calibration.get("physical_case")
    mesh = pilot.get("mesh")
    if not isinstance(source, dict) or not isinstance(physical, dict) \
            or not isinstance(mesh, dict):
        raise ValueError("pilot/calibration identity metadata is incomplete")
    exact_pairs = (
        (pilot.get("case_id"), source.get("campaign_case"), "case id"),
        (pilot.get("source_state_sha256"), source.get("source_state_sha256"),
         "source state hash"),
        (mesh.get("root_dimensions"), physical.get("root_dimensions"),
         "root dimensions"),
        (mesh.get("meshblock_dimensions"), physical.get("meshblock_dimensions"),
         "MeshBlock dimensions"),
        (mesh.get("physical_refinement_levels"),
         physical.get("physical_refinement_levels"), "refinement levels"),
        (mesh.get("maximum_meshblocks_per_rank"),
         physical.get("configured_maximum_meshblocks_per_rank"),
         "MeshBlock capacity"),
    )
    for left, right, label in exact_pairs:
        if left != right:
            raise ValueError(f"pilot and calibration disagree on {label}")


def build_production_campaign(
    pilot: dict[str, object],
    calibration: dict[str, object],
    *,
    target_segment_wall_hours: float = 3.0,
    checkpoint_wall_hours: float = 1.5,
    history_root_steps: int = 1,
    field_outputs_per_crossing: int = 4,
    force_outer_capture_fractions: tuple[float, float, float] = (0.5, 1.0, 2.0),
    retained_restart_generations: int = 2,
) -> dict[str, object]:
    _validate_identity(pilot, calibration)
    target_segment_wall_hours = _positive(
        target_segment_wall_hours, "target segment wall hours"
    )
    checkpoint_wall_hours = _positive(
        checkpoint_wall_hours, "checkpoint wall hours"
    )
    history_root_steps = _positive_integer(history_root_steps, "history root steps")
    field_outputs_per_crossing = _positive_integer(
        field_outputs_per_crossing, "field outputs per crossing"
    )
    retained_restart_generations = _positive_integer(
        retained_restart_generations, "retained restart generations"
    )
    if len(force_outer_capture_fractions) != 3 or not all(
        math.isfinite(value) and value > 0.0
        for value in force_outer_capture_fractions
    ) or not all(
        left < right for left, right in zip(
            force_outer_capture_fractions[:-1], force_outer_capture_fractions[1:]
        )
    ):
        raise ValueError("force outer capture fractions must be three increasing positives")

    result = calibration["result"]
    physical = calibration["physical_case"]
    direct = pilot["bhl_plan"]["costs"]["direct"]
    scales = pilot["bhl_plan"]["scales"]
    mesh = pilot["mesh"]
    final_cycle = _positive_integer(result["final_cycle"], "calibration final cycle")
    final_time = _positive(result["final_time_in_secondary_masses"],
                           "calibration final time")
    measured_root_dt = final_time / final_cycle
    root_cycle_seconds = _positive(
        result["steady_nine_level_root_cycle_seconds"],
        "steady root-cycle seconds",
    )
    measured_rate = _positive(
        result["zone_cycles_per_second_reported"], "measured zone-cycle rate"
    )
    duration = _positive(scales["direct_settling_duration"], "target duration")
    crossing_time = _positive(scales["capture_crossing_time"], "crossing time")
    capture_radius = _positive(
        scales["capture_radius_factor_two_for_cost"], "capture radius"
    )
    levels = _positive_integer(direct["nested_refinement_levels"],
                               "physical refinement levels")
    finest_steps = _positive_integer(direct["finest_steps"], "planner finest steps")
    ideal_updates = _positive_integer(
        direct["ideal_level_subcycled_zone_updates"], "ideal zone updates"
    )

    projected_root_cycles = math.ceil(duration / measured_root_dt)
    planner_root_cycle_proxy = finest_steps / 2**levels
    step_count_correction = projected_root_cycles / planner_root_cycle_proxy
    measured_zone_cycles_per_root = root_cycle_seconds * measured_rate
    ideal_zone_cycles_per_root = ideal_updates / planner_root_cycle_proxy
    topology_correction = measured_zone_cycles_per_root / ideal_zone_cycles_per_root
    corrected_lower_seconds = ideal_updates * step_count_correction / measured_rate
    empirical_seconds = projected_root_cycles * root_cycle_seconds

    segment_root_steps = max(
        1, math.floor(target_segment_wall_hours * 3600.0 / root_cycle_seconds)
    )
    checkpoint_root_steps = max(
        1, round(checkpoint_wall_hours * 3600.0 / root_cycle_seconds)
    )
    checkpoint_root_steps = min(checkpoint_root_steps, segment_root_steps)
    field_root_steps = max(
        1, round(crossing_time / field_outputs_per_crossing / measured_root_dt)
    )
    history_dt = history_root_steps * measured_root_dt
    checkpoint_dt = checkpoint_root_steps * measured_root_dt
    field_dt = field_root_steps * measured_root_dt
    segment_duration_proxy = segment_root_steps * measured_root_dt
    segment_wall_hours_proxy = segment_root_steps * root_cycle_seconds / 3600.0

    force_radii = [capture_radius * value for value in force_outer_capture_fractions]
    conservative_half_width = min(
        min(abs(float(low)), abs(float(high)))
        for low, high in zip(mesh["lower"], mesh["upper"])
    )
    if force_radii[-1] >= conservative_half_width:
        raise ValueError("outer force shell does not fit inside the direct domain")

    meshblock_cells = int(mesh["meshblock_dimensions"][0])
    if mesh["meshblock_dimensions"] != [meshblock_cells] * 3:
        raise ValueError("restart estimate requires cubic MeshBlocks")
    restart_bytes_per_block = _restart_bytes_per_meshblock(meshblock_cells)
    final_meshblocks = _positive_integer(result["final_meshblocks"],
                                         "final MeshBlocks")
    capacity = _positive_integer(mesh["maximum_meshblocks_per_rank"],
                                 "MeshBlock capacity")
    if final_meshblocks >= capacity:
        raise ValueError("calibrated topology has no positive MeshBlock capacity headroom")
    restart_final_bytes = restart_bytes_per_block * final_meshblocks
    restart_capacity_bytes = restart_bytes_per_block * capacity
    active_cells_per_block = meshblock_cells**3
    primitive_final_bytes = (
        final_meshblocks * active_cells_per_block
        * MHD_PRIMITIVE_OUTPUT_VARIABLES * 4
    )
    divb_final_bytes = final_meshblocks * active_cells_per_block * 4
    field_count = math.floor(duration / field_dt) + 2
    retained_field_bytes = field_count * (primitive_final_bytes + divb_final_bytes)
    rolling_restart_bytes = retained_restart_generations * restart_capacity_bytes
    working_disk_bytes = rolling_restart_bytes + retained_field_bytes + 20 * 2**30

    basename = _safe_basename(str(pilot["case_id"]))
    replacements = {
        "job/basename": basename,
        "time/nlim": str(segment_root_steps),
        "time/tlim": f"{duration:.17g}",
        "problem/user_hist": "true",
        "output1/dt": f"{field_dt:.17g}",
        "output2/dt": f"{field_dt:.17g}",
        "output3/dt": f"{checkpoint_dt:.17g}",
        "output4/dt": f"{history_dt:.17g}",
        "output5/dt": "0",
    }
    fresh_overrides = _replace_overrides(pilot["athena_overrides"], replacements)
    # These parameters already exist in the smoke template, so command-line replacement
    # remains compatible with AthenaK's strict "known key only" override parser.
    fresh_overrides.extend([
        "output1/variable=mhd_w_bcc",
        "output2/variable=mhd_divb",
        "output4/user_hist_only=true",
        "problem/force_frame=source_tetrad",
        "problem/force_surface_radius=3",
        f"problem/force_outer_radius_1={force_radii[0]:.17g}",
        f"problem/force_outer_radius_2={force_radii[1]:.17g}",
        f"problem/force_outer_radius_3={force_radii[2]:.17g}",
        "problem/force_surface_nlevel=5",
        "problem/force_subtract_background=true",
    ])
    keys = [item.split("=", 1)[0] for item in fresh_overrides]
    if len(keys) != len(set(keys)):
        duplicates = sorted({key for key in keys if keys.count(key) > 1})
        raise ValueError("production overrides contain duplicates: " + ", ".join(duplicates))

    earliest_stationarity_time = 1.5 * crossing_time
    return {
        "classification": CLASSIFICATION,
        "case_id": pilot["case_id"],
        "source_state_sha256": pilot["source_state_sha256"],
        "source_pilot_classification": pilot["classification"],
        "source_calibration": {
            "classification": calibration["classification"],
            "date": calibration.get("date"),
            "solver_commit": calibration["source"]["solver_commit"],
            "athena_sha256": calibration["build"]["athena_sha256"],
        },
        "science_scope": (
            "full-primary, frozen-gradient direct BHL baseline; force stationarity and "
            "far-wake shell closure are runtime gates, not assumptions"
        ),
        "physical_scales": {
            "capture_radius_in_secondary_masses": capture_radius,
            "capture_crossing_time_in_secondary_masses": crossing_time,
            "target_crossings": duration / crossing_time,
            "target_duration_in_secondary_masses": duration,
            "force_outer_capture_fractions": list(force_outer_capture_fractions),
            "force_outer_radii_in_secondary_masses": force_radii,
            "conservative_domain_half_width_in_secondary_masses":
                conservative_half_width,
        },
        "runtime_projection": {
            "measured_root_dt_in_secondary_masses": measured_root_dt,
            "projected_root_cycles": projected_root_cycles,
            "planner_root_cycle_proxy": planner_root_cycle_proxy,
            "measured_to_planner_step_count_correction": step_count_correction,
            "measured_zone_cycles_per_steady_root_cycle":
                measured_zone_cycles_per_root,
            "planner_ideal_zone_cycles_per_root_proxy": ideal_zone_cycles_per_root,
            "measured_to_planner_topology_correction": topology_correction,
            "uncorrected_ideal_hours": ideal_updates / measured_rate / 3600.0,
            "cfl_corrected_lower_bound_hours": corrected_lower_seconds / 3600.0,
            "empirical_steady_topology_hours": empirical_seconds / 3600.0,
            "budget_rule": (
                "use empirical_steady_topology_hours before diagnostic, checkpoint, "
                "field-output, startup, and stationarity-extension overhead"
            ),
        },
        "segmentation": {
            "target_wall_hours": target_segment_wall_hours,
            "root_steps_per_segment": segment_root_steps,
            "duration_per_segment_proxy": segment_duration_proxy,
            "wall_hours_per_segment_proxy": segment_wall_hours_proxy,
            "segments_if_root_dt_is_stationary":
                math.ceil(projected_root_cycles / segment_root_steps),
            "checkpoint_root_steps": checkpoint_root_steps,
            "checkpoint_dt_in_secondary_masses": checkpoint_dt,
            "checkpoint_wall_hours_proxy":
                checkpoint_root_steps * root_cycle_seconds / 3600.0,
            "resume_rule": (
                "read the newest audited restart, set nlim=current_cycle+"
                f"{segment_root_steps}, retain the global tlim, and never infer the "
                "next restart filename from a failed segment"
            ),
        },
        "outputs": {
            "history_root_steps": history_root_steps,
            "history_dt_in_secondary_masses": history_dt,
            "field_root_steps": field_root_steps,
            "field_dt_in_secondary_masses": field_dt,
            "field_outputs_per_crossing_proxy": crossing_time / field_dt,
            "primitive_variable": "mhd_w_bcc",
            "constraint_variable": "mhd_divb",
            "restart_cadence_root_steps": checkpoint_root_steps,
        },
        "resource_envelope": {
            "gpu": "single A100-SXM4-80GB or a requalified equivalent",
            "calibrated_peak_gpu_memory_MiB_without_restart":
                result["peak_gpu_memory_MiB"],
            "calibrated_final_meshblocks": final_meshblocks,
            "configured_meshblock_capacity": capacity,
            "meshblock_capacity_headroom": capacity - final_meshblocks,
            "restart_layout_assumptions": {
                "real_bytes": REAL_BYTES,
                "nghost": 4,
                "mhd_conserved_variables": MHD_CONSERVED_VARIABLES,
                "prescribed_adm_variables": PRESCRIBED_ADM_VARIABLES,
                "bytes_per_meshblock": restart_bytes_per_block,
            },
            "estimated_restart_GiB_at_calibrated_topology":
                restart_final_bytes / 2**30,
            "estimated_restart_GiB_at_capacity": restart_capacity_bytes / 2**30,
            "retained_restart_generations": retained_restart_generations,
            "estimated_retained_field_GiB": retained_field_bytes / 2**30,
            "minimum_working_disk_GiB": math.ceil(working_disk_bytes / 2**30),
            "warning": (
                "restart output duplicates U/B and prescribed ADM arrays on the device; "
                "measure peak memory, file size, and write latency before a long run"
            ),
        },
        "stationarity_gate": {
            "status": "runtime_evidence_required",
            "discard_before_time_in_secondary_masses": crossing_time,
            "comparison_window_in_crossings": 0.25,
            "earliest_assessment_time_in_secondary_masses":
                earliest_stationarity_time,
            "minimum_stationary_average_in_crossings": 0.5,
            "required_observables": [
                "mdot_hat", "Ftotal1H", "Ftotal2H", "Ftotal3H",
                "dFrel21H", "dFrel32H",
            ],
            "provisional_gates": {
                "maximum_adjacent_window_vector_drift_fraction": 0.10,
                "maximum_adjacent_window_mean_difference_in_combined_SE": 2.0,
                "maximum_outer_shell_increment_fraction": 0.15,
                "minimum_blocks_per_window": 8,
            },
            "extension_rule": (
                "if any gate fails or the block length is shorter than the measured "
                "autocorrelation time, extend by 0.25 capture crossings and reassess"
            ),
        },
        "athena_input": pilot["athena_input"],
        "fresh_argv_template": [
            "{athena}", "-i", pilot["athena_input"], "-d", "{run_directory}",
            *fresh_overrides,
        ],
        "restart_argv_template": [
            "{athena}", "-r", "{audited_restart}", "-d", "{run_directory}",
            f"time/nlim={{current_cycle_plus_{segment_root_steps}}}",
            f"time/tlim={duration:.17g}",
        ],
        "fresh_overrides": fresh_overrides,
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target-segment-wall-hours", type=float, default=3.0)
    parser.add_argument("--checkpoint-wall-hours", type=float, default=1.5)
    parser.add_argument("--history-root-steps", type=int, default=1)
    parser.add_argument("--field-outputs-per-crossing", type=int, default=4)
    parser.add_argument(
        "--force-outer-capture-fractions", type=float, nargs=3,
        default=(0.5, 1.0, 2.0), metavar=("R1", "R2", "R3"),
    )
    parser.add_argument("--retained-restart-generations", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    pilot = json.loads(arguments.pilot.expanduser().resolve(strict=True).read_text())
    calibration = json.loads(
        arguments.calibration.expanduser().resolve(strict=True).read_text()
    )
    campaign = build_production_campaign(
        pilot,
        calibration,
        target_segment_wall_hours=arguments.target_segment_wall_hours,
        checkpoint_wall_hours=arguments.checkpoint_wall_hours,
        history_root_steps=arguments.history_root_steps,
        field_outputs_per_crossing=arguments.field_outputs_per_crossing,
        force_outer_capture_fractions=tuple(arguments.force_outer_capture_fractions),
        retained_restart_generations=arguments.retained_restart_generations,
    )
    output = arguments.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(campaign, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(output)
    runtime = campaign["runtime_projection"]
    segmentation = campaign["segmentation"]
    print(
        f"empirical_hours={runtime['empirical_steady_topology_hours']:.3f} "
        f"segments={segmentation['segments_if_root_dt_is_stationary']} "
        f"steps_per_segment={segmentation['root_steps_per_segment']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
