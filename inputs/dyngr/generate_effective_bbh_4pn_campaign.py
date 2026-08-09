#!/usr/bin/env python3
"""Validate a 4PN BBH trajectory and emit one gated L/M/H AthenaK input."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import struct
import sys
import tempfile
from typing import Iterable


HERE = Path(__file__).resolve().parent
DEFAULT_MATRIX = HERE / "effective_bbh_4pn_campaign.json"
POSITION = ((0, 1, 2), (3, 4, 5))
VELOCITY = ((6, 7, 8), (9, 10, 11))
SPIN = ((12, 13, 14), (15, 16, 17))
MASS = (18, 19)


class CampaignError(ValueError):
    """A fail-closed campaign validation error."""


@dataclass(frozen=True)
class Target:
    position: tuple[float, float, float]
    velocity: tuple[float, float, float]
    spin: tuple[float, float, float]
    mass: float


@dataclass(frozen=True)
class TrajectorySummary:
    rows: int
    first_time: float
    last_time: float
    state_at_start: tuple[float, ...]
    state_at_end: tuple[float, ...]
    targets_at_start: tuple[Target, ...]
    targets_at_end: tuple[Target, ...]
    maximum_guard_radius: float
    maximum_row_guard_radius: float
    maximum_middle_control_speed: float
    maximum_speed: float
    sha256: str


@dataclass(frozen=True)
class ProvenanceSummary:
    path: Path
    sha256: str
    model_labels: tuple[str, ...]
    source_revision: str
    source_artifact: Path
    source_sha256: str
    remnant_mass: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a resource-gated AthenaK input for one tier of the "
            "effective-BBH 4PN AMR convergence campaign."
        )
    )
    parser.add_argument("tier", choices=("L", "M", "H"))
    parser.add_argument("--trajectory", required=True, type=Path)
    parser.add_argument(
        "--trajectory-provenance",
        required=True,
        type=Path,
        help="v2 JSON sidecar emitted by the canonical single-remnant stitcher",
    )
    parser.add_argument(
        "--runtime-trajectory",
        default=None,
        help="path written into the input (default: absolute validation path)",
    )
    parser.add_argument(
        "--source-artifact",
        type=Path,
        default=None,
        help="local copy of trajectory_model.source_provenance.path to re-hash",
    )
    parser.add_argument(
        "--expected-source-revision",
        default=None,
        help="exact 40-hex revision; it must also match the audited matrix",
    )
    parser.add_argument(
        "--trajectory-time-offset",
        type=float,
        default=None,
        help="table time at simulation t=0 (default: first table time + metric_fd_step)",
    )
    parser.add_argument(
        "--tlim",
        type=float,
        default=None,
        help="simulation duration in M (default: table end minus FD padding)",
    )
    parser.add_argument(
        "--metric-fd-step",
        type=float,
        default=None,
        help="metric time-difference half-width (default from matrix)",
    )
    parser.add_argument(
        "--cfl-number",
        type=float,
        default=None,
        help="CFL for a temporal-convergence variant (default from matrix; may only decrease)",
    )
    parser.add_argument(
        "--real-precision",
        choices=("double", "single"),
        default="double",
        help="publication campaign requires double (default)",
    )
    parser.add_argument(
        "--allow-single-precision-systematic",
        action="store_true",
        help="permit a representable short single-precision non-certificate run",
    )
    parser.add_argument("--gpus", type=int, required=True)
    parser.add_argument("--gpu-memory-gib", type=float, default=40.0)
    parser.add_argument("--scratch-gib", type=float, required=True)
    parser.add_argument(
        "--streaming-drain",
        action="store_true",
        help="use the smaller segmented-run scratch gate",
    )
    parser.add_argument(
        "--drain-mib-s",
        type=float,
        default=None,
        help="sustained verified NAS drain rate; required with --streaming-drain",
    )
    parser.add_argument("--basename", default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def load_matrix(path: Path) -> dict[str, object]:
    try:
        with path.open("r", encoding="utf-8") as source:
            matrix = json.load(source)
    except (OSError, json.JSONDecodeError) as error:
        raise CampaignError(f"cannot load campaign matrix {path}: {error}") from error
    validate_matrix(matrix)
    return matrix


def validate_matrix(matrix: dict[str, object]) -> None:
    if matrix.get("schema_version") != 2:
        raise CampaignError("unsupported campaign matrix schema")
    common = matrix["common"]
    mesh = common["mesh"]
    time = common["time"]
    amr = common["amr"]
    if time["subcycling"] != "level" or time["refinement_time_ratio"] != 2:
        raise CampaignError("campaign must use strict per-level 2:1 subcycling")
    if time.get("publication_precision") != "double":
        raise CampaignError("publication campaign must require double precision")
    fd_step = float(time["metric_fd_step_M"])
    if tuple(float(value) for value in time["metric_fd_step_variants_M"]) != (
        0.5 * fd_step,
        fd_step,
        2.0 * fd_step,
    ):
        raise CampaignError("metric FD variants must be exactly h/2, h, and 2h")
    if amr["refinement"] != "adaptive" or amr["criterion"] != "track":
        raise CampaignError("campaign must use moving-hole adaptive refinement")
    if amr["ncycle_check"] != 1 or amr["refinement_interval"] != 1:
        raise CampaignError("moving-hole AMR must re-evaluate at every synchronized cycle")
    if common["resource_model"].get("empirical_memory_precision") != "double":
        raise CampaignError("campaign resource model must identify its double-precision basis")
    resource = common["resource_model"]
    rounding = int(resource["max_nmb_rounding"])
    rank_margin = float(resource["per_rank_partition_margin_fraction"])
    if rounding <= 0 or not 0.0 < rank_margin < 1.0:
        raise CampaignError("invalid per-rank partition rounding or safety margin")
    if resource.get("load_balance_partition_ordering") != (
        "contiguous MeshBlock GIDs (Z-order)"
    ):
        raise CampaignError("campaign partition audit must use contiguous MeshBlock GIDs")
    outputs = common["outputs"]
    output_bytes = outputs["estimated_bytes_per_meshblock"]
    output_specs = (
        ("full_state_dt_M", "full_state_mhd_w_bcc"),
        ("gr_diagnostics_dt_M", "gr_diagnostics"),
        ("divb_dt_M", "divb"),
        ("restart_dt_M", "restart_double_precision_dynadm"),
    )
    if not math.isfinite(float(outputs["history_dt_M"])) or float(
        outputs["history_dt_M"]
    ) <= 0.0:
        raise CampaignError("history output cadence must be finite and positive")
    for cadence_key, bytes_key in output_specs:
        if (
            not math.isfinite(float(outputs[cadence_key]))
            or float(outputs[cadence_key]) <= 0.0
            or int(output_bytes[bytes_key]) <= 0
        ):
            raise CampaignError(
                f"invalid output cadence or per-MeshBlock size for {cadence_key}"
            )
    root_dx = (
        2.0 * float(mesh["domain_half_extent_M"])
        / int(mesh["root_cells_per_dimension"])
    )
    if not math.isclose(root_dx, float(mesh["root_dx_M"]), rel_tol=0.0, abs_tol=1e-14):
        raise CampaignError("matrix root_dx_M is inconsistent with its domain and root mesh")

    common_outer: tuple[float, ...] | None = None
    for name, tier in matrix["tiers"].items():
        max_level = int(tier["finest_physical_level"])
        if int(tier["num_levels"]) != max_level + 1:
            raise CampaignError(f"tier {name}: num_levels must include root L0")
        radii = tier["shell_radii_M"]
        if set(radii) != {str(level) for level in range(1, max_level + 1)}:
            raise CampaignError(f"tier {name}: every physical-level radius is required")
        ordered = tuple(float(radii[str(level)]) for level in range(1, max_level + 1))
        if any(not math.isfinite(value) or value <= 0.0 for value in ordered):
            raise CampaignError(f"tier {name}: shell radii must be finite and positive")
        if any(finer > coarser for coarser, finer in zip(ordered, ordered[1:])):
            raise CampaignError(f"tier {name}: shell radii must be non-increasing")
        if common_outer is None:
            common_outer = ordered[:8]
        elif ordered[:8] != common_outer:
            raise CampaignError("L/M/H outer shell radii through L8 must be identical")
        expected_dx = root_dx / (2**max_level)
        if not math.isclose(
            expected_dx, float(tier["finest_dx_M"]), rel_tol=0.0, abs_tol=1e-14
        ):
            raise CampaignError(f"tier {name}: finest_dx_M is inconsistent")
        topology = tier["topology_estimate"]
        margin = float(common["resource_model"]["topology_margin_fraction"])
        gate_candidates = [
            int(topology["initial_preseed_meshblocks_athena_m"]),
            int(topology["initial_preseed_alignment_probe_max_meshblocks"]),
            int(topology["postmerger_2p5M_guard_meshblocks"]),
        ]
        observation = topology.get("rootcap_runtime_observation")
        if name in {"L", "M"} and not isinstance(observation, dict):
            raise CampaignError(
                f"tier {name}: a rootcap runtime observation is required"
            )
        if observation is not None:
            if not isinstance(observation, dict):
                raise CampaignError(
                    f"tier {name}: rootcap runtime observation must be an object"
                )
            observation_sha = observation.get("evidence_sha256")
            observation_commit = observation.get("commit")
            if not isinstance(observation_sha, str) or re.fullmatch(
                r"[0-9a-f]{64}", observation_sha
            ) is None:
                raise CampaignError(
                    f"tier {name}: runtime-observation evidence SHA-256 is invalid"
                )
            if not isinstance(observation_commit, str) or re.fullmatch(
                r"[0-9a-f]{40}", observation_commit
            ) is None:
                raise CampaignError(
                    f"tier {name}: runtime-observation commit is invalid"
                )
            observed_root_dt = float(observation["root_dt_max_M"])
            expected_root_dt = float(time["cfl_number"]) * root_dx
            if not math.isclose(
                observed_root_dt, expected_root_dt, rel_tol=0.0, abs_tol=1e-14
            ):
                raise CampaignError(
                    f"tier {name}: runtime observation used a different root-step cap"
                )
            observed_peak = int(observation["peak_meshblocks"])
            executed_peak = int(observation["executed_peak_meshblocks"])
            observed_ranks = int(observation["mpi_ranks"])
            observed_rank_peak = int(observation["peak_meshblocks_per_rank"])
            if (
                observed_peak <= 0
                or executed_peak <= 0
                or observed_peak < executed_peak
                or observed_ranks <= 0
                or observed_rank_peak < math.ceil(executed_peak / observed_ranks)
            ):
                raise CampaignError(
                    f"tier {name}: invalid rootcap runtime topology observation"
                )
            gate_candidates.append(observed_peak)
        required_gate = math.ceil(max(gate_candidates) * (1.0 + margin))
        if int(topology["campaign_meshblock_gate"]) != required_gate:
            raise CampaignError(
                f"tier {name}: MeshBlock gate must equal the larger measured topology "
                f"plus {margin:.0%}"
            )
        qualification = topology.get("rank_partition_qualification")
        if not isinstance(qualification, dict):
            raise CampaignError(f"tier {name}: per-rank partition qualification is required")
        artifact_sha = qualification.get("mesh_structure_sha256")
        if not isinstance(artifact_sha, str) or re.fullmatch(
            r"[0-9a-f]{64}", artifact_sha
        ) is None:
            raise CampaignError(f"tier {name}: partition artifact SHA-256 is invalid")
        qualified_ranks = int(qualification["mpi_ranks"])
        qualified_blocks = int(qualification["topology_meshblocks"])
        block_range = tuple(int(value) for value in qualification["meshblocks_per_rank_range"])
        cost_range = tuple(int(value) for value in qualification["weighted_cost_per_rank_range"])
        if qualified_ranks <= 0 or qualified_blocks <= 0:
            raise CampaignError(f"tier {name}: invalid partition qualification size")
        if len(block_range) != 2 or not 0 < block_range[0] <= block_range[1]:
            raise CampaignError(f"tier {name}: invalid per-rank MeshBlock range")
        if len(cost_range) != 2 or not 0 < cost_range[0] <= cost_range[1]:
            raise CampaignError(f"tier {name}: invalid per-rank weighted-cost range")
        if qualified_blocks not in {
            int(topology["initial_preseed_meshblocks_athena_m"]),
            int(topology["initial_preseed_alignment_probe_max_meshblocks"]),
        }:
            raise CampaignError(
                f"tier {name}: partition qualification is not tied to a measured topology"
            )
        hard_floor = int(qualification["max_nmb_per_rank_hard_floor"])
        required_floor = (
            math.ceil(block_range[1] * (1.0 + rank_margin) / rounding) * rounding
        )
        if hard_floor != required_floor:
            raise CampaignError(
                f"tier {name}: max_nmb per-rank floor must be the measured peak plus "
                f"{rank_margin:.0%}, rounded to {rounding}"
            )
        resource_gate = tier["resource_gate"]
        if int(resource_gate["partition_qualified_launch_gpus"]) != qualified_ranks:
            raise CampaignError(
                f"tier {name}: partition-qualified launch count must match the audit"
            )
        per_block_mib = float(resource["empirical_peak_MiB_per_meshblock"])
        usable_mib = 40.0 * 1024.0 * float(resource["gpu_memory_usable_fraction"])
        aggregate_lower_bound = math.ceil(required_gate / math.floor(usable_mib / per_block_mib))
        if int(resource_gate["aggregate_memory_lower_bound_A100_40G_gpus"]) != (
            aggregate_lower_bound
        ):
            raise CampaignError(
                f"tier {name}: aggregate A100-40G memory lower bound is inconsistent"
            )
        minimum_gpu_memory_gib = float(
            resource_gate["minimum_gpu_memory_GiB_at_qualified_count"]
        )
        qualified_capacity = math.floor(
            minimum_gpu_memory_gib
            * 1024.0
            * float(resource["gpu_memory_usable_fraction"])
            / per_block_mib
        )
        qualified_gate_per_rank = (
            math.ceil((required_gate / qualified_ranks) / rounding) * rounding
        )
        qualified_allocation = max(qualified_gate_per_rank, hard_floor)
        if qualified_allocation > qualified_capacity:
            raise CampaignError(
                f"tier {name}: declared qualified GPU memory cannot hold the "
                "campaign per-rank allocation"
            )
        reference_duration = float(resource["reference_duration_M"])
        segment_span = float(resource["streaming_segment_span_M"])
        restart_generations = int(resource["streaming_restart_generations"])
        streaming_margin = float(resource["streaming_safety_fraction"])
        if segment_span <= 0.0 or restart_generations < 2 or streaming_margin < 0.0:
            raise CampaignError("streaming storage policy is invalid")
        effective_restart_cadence = min(
            float(outputs["restart_dt_M"]), segment_span
        )
        projected_reference_gib = sum(
            (math.floor(reference_duration / float(cadence)) + 2)
            * required_gate
            * int(output_bytes[key])
            for cadence, key in (
                (outputs["full_state_dt_M"], "full_state_mhd_w_bcc"),
                (outputs["gr_diagnostics_dt_M"], "gr_diagnostics"),
                (outputs["divb_dt_M"], "divb"),
                (effective_restart_cadence, "restart_double_precision_dynadm"),
            )
        ) / 2**30
        reported_projection = float(
            resource_gate["projected_undrained_output_GiB_at_10000M"]
        )
        if not (
            projected_reference_gib <= reported_projection
            < projected_reference_gib + 0.1
        ):
            raise CampaignError(
                f"tier {name}: reference output projection is inconsistent"
            )
        minimum_undrained = int(
            resource_gate["minimum_undrained_scratch_GiB_at_10000M"]
        )
        if minimum_undrained < math.ceil(1.25 * projected_reference_gib):
            raise CampaignError(
                f"tier {name}: undrained scratch does not include the 25% margin"
            )
        segment_bytes = sum(
            (math.floor(segment_span / float(cadence)) + 2)
            * required_gate
            * int(output_bytes[key])
            for cadence, key in (
                (outputs["full_state_dt_M"], "full_state_mhd_w_bcc"),
                (outputs["gr_diagnostics_dt_M"], "gr_diagnostics"),
                (outputs["divb_dt_M"], "divb"),
            )
        )
        segment_bytes += (
            restart_generations
            * required_gate
            * int(output_bytes["restart_double_precision_dynadm"])
        )
        minimum_streaming = int(resource_gate["minimum_streaming_scratch_GiB"])
        required_streaming = math.ceil((1.0 + streaming_margin) * segment_bytes / 2**30)
        if minimum_streaming < required_streaming:
            raise CampaignError(
                f"tier {name}: streaming scratch cannot retain a guarded segment "
                "and the configured restart generations"
            )


def vector_norm(values: Iterable[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def close(first: float, second: float) -> bool:
    scale = max(1.0, abs(first), abs(second))
    return abs(first - second) <= 128.0 * sys.float_info.epsilon * scale


def tracking_targets(state: tuple[float, ...]) -> tuple[Target, ...]:
    targets = []
    for hole in range(2):
        mass = state[MASS[hole]]
        if mass <= 0.0:
            continue
        targets.append(
            Target(
                tuple(state[index] for index in POSITION[hole]),
                tuple(state[index] for index in VELOCITY[hole]),
                tuple(state[index] for index in SPIN[hole]),
                mass,
            )
        )
    if len(targets) != 2:
        return tuple(targets)
    coincident = all(
        close(first, second)
        for first_group, second_group in (
            (targets[0].position, targets[1].position),
            (targets[0].velocity, targets[1].velocity),
            (targets[0].spin, targets[1].spin),
        )
        for first, second in zip(first_group, second_group)
    )
    if not coincident:
        return tuple(targets)
    return (
        Target(
            tuple((a + b) / 2.0 for a, b in zip(targets[0].position, targets[1].position)),
            tuple((a + b) / 2.0 for a, b in zip(targets[0].velocity, targets[1].velocity)),
            tuple((a + b) / 2.0 for a, b in zip(targets[0].spin, targets[1].spin)),
            targets[0].mass + targets[1].mass,
        ),
    )


def horizon_guard_radius(target: Target, factor: float) -> float:
    speed2 = sum(value * value for value in target.velocity)
    if not speed2 < 1.0:
        raise CampaignError("trajectory contains a superluminal target")
    spin2 = sum(value * value for value in target.spin)
    horizon = target.mass + math.sqrt(max(target.mass * target.mass - spin2, 0.0))
    rest_enclosing = math.sqrt(spin2 + horizon * horizon)
    return factor * rest_enclosing / math.sqrt(1.0 - speed2)


def minimum_segment_norm(
    first: tuple[float, float, float], second: tuple[float, float, float]
) -> float:
    direction = tuple(b - a for a, b in zip(first, second))
    norm2 = sum(value * value for value in direction)
    weight = 0.0 if norm2 == 0.0 else max(
        0.0, min(1.0, -sum(a * d for a, d in zip(first, direction)) / norm2)
    )
    return vector_norm(a + weight * d for a, d in zip(first, direction))


def interval_guard_bound(
    lower_time: float,
    lower: tuple[float, ...],
    upper_time: float,
    upper: tuple[float, ...],
    hole: int,
    factor: float,
    merged: bool = False,
) -> tuple[float, float]:
    interval = upper_time - lower_time
    if merged:
        position0 = tuple(
            (lower[first] + lower[second]) / 2.0
            for first, second in zip(POSITION[0], POSITION[1])
        )
        position1 = tuple(
            (upper[first] + upper[second]) / 2.0
            for first, second in zip(POSITION[0], POSITION[1])
        )
        velocity0 = tuple(
            (lower[first] + lower[second]) / 2.0
            for first, second in zip(VELOCITY[0], VELOCITY[1])
        )
        velocity1 = tuple(
            (upper[first] + upper[second]) / 2.0
            for first, second in zip(VELOCITY[0], VELOCITY[1])
        )
    else:
        position0 = tuple(lower[index] for index in POSITION[hole])
        position1 = tuple(upper[index] for index in POSITION[hole])
        velocity0 = tuple(lower[index] for index in VELOCITY[hole])
        velocity1 = tuple(upper[index] for index in VELOCITY[hole])
    middle_control = tuple(
        3.0 * (p1 - p0) / interval - v0 - v1
        for p0, p1, v0, v1 in zip(position0, position1, velocity0, velocity1)
    )
    middle_speed = vector_norm(middle_control)
    if not math.isfinite(middle_speed) or not middle_speed < 1.0:
        raise CampaignError(
            f"trajectory interval [{lower_time:g}, {upper_time:g}] hole {hole + 1} "
            f"has non-subluminal Hermite middle control |v|={middle_speed:.17g}"
        )
    speed_bound = max(
        vector_norm(velocity0),
        middle_speed,
        vector_norm(velocity1),
    )
    if merged:
        mass0 = lower[MASS[0]] + lower[MASS[1]]
        mass1 = upper[MASS[0]] + upper[MASS[1]]
        spin0 = tuple((lower[a] + lower[b]) / 2.0 for a, b in zip(SPIN[0], SPIN[1]))
        spin1 = tuple((upper[a] + upper[b]) / 2.0 for a, b in zip(SPIN[0], SPIN[1]))
    else:
        mass0, mass1 = lower[MASS[hole]], upper[MASS[hole]]
        spin0 = tuple(lower[index] for index in SPIN[hole])
        spin1 = tuple(upper[index] for index in SPIN[hole])
    maximum_mass = max(mass0, mass1)
    if maximum_mass <= 0.0:
        return 0.0, middle_speed
    minimum_spin = minimum_segment_norm(spin0, spin1)
    maximum_spin = max(vector_norm(spin0), vector_norm(spin1))
    subextreme_spin = min(minimum_spin, maximum_mass)
    subextreme_radius = math.sqrt(
        2.0
        * maximum_mass
        * (
            maximum_mass
            + math.sqrt(max(maximum_mass * maximum_mass - subextreme_spin**2, 0.0))
        )
    )
    superextreme_radius = math.sqrt(maximum_spin**2 + maximum_mass**2)
    rest_radius_bound = max(subextreme_radius, superextreme_radius)
    guard_bound = factor * rest_radius_bound / math.sqrt(1.0 - speed_bound**2)
    return guard_bound, middle_speed


def terms_coincident(state: tuple[float, ...]) -> bool:
    return all(
        close(state[first], state[second])
        for first_group, second_group in (
            (POSITION[0], POSITION[1]),
            (VELOCITY[0], VELOCITY[1]),
            (SPIN[0], SPIN[1]),
        )
        for first, second in zip(first_group, second_group)
    )


def interpolate_state(
    lower: tuple[float, tuple[float, ...]],
    upper: tuple[float, tuple[float, ...]],
    time: float,
) -> tuple[float, ...]:
    t0, state0 = lower
    t1, state1 = upper
    if time == t0 or t0 == t1:
        return state0
    if time == t1:
        return state1
    interval = t1 - t0
    weight = (time - t0) / interval
    state = [a + weight * (b - a) for a, b in zip(state0, state1)]
    w2 = weight * weight
    w3 = w2 * weight
    h00 = 2.0 * w3 - 3.0 * w2 + 1.0
    h10 = w3 - 2.0 * w2 + weight
    h01 = -2.0 * w3 + 3.0 * w2
    h11 = w3 - w2
    for hole in range(2):
        for position, velocity in zip(POSITION[hole], VELOCITY[hole]):
            p0, p1 = state0[position], state1[position]
            v0, v1 = state0[velocity], state1[velocity]
            state[position] = h00 * p0 + h10 * interval * v0 + h01 * p1 + h11 * interval * v1
            state[velocity] = (
                ((6.0 * w2 - 6.0 * weight) / interval) * p0
                + (3.0 * w2 - 4.0 * weight + 1.0) * v0
                + ((-6.0 * w2 + 6.0 * weight) / interval) * p1
                + (3.0 * w2 - 2.0 * weight) * v1
            )
    return tuple(state)


def inspect_trajectory(
    path: Path, start_time: float, horizon_factor: float, guard_limit: float
) -> TrajectorySummary:
    if not path.is_file():
        raise CampaignError(f"trajectory does not exist: {path}")
    digest = hashlib.sha256()
    first: tuple[float, tuple[float, ...]] | None = None
    previous: tuple[float, tuple[float, ...]] | None = None
    lower: tuple[float, tuple[float, ...]] | None = None
    upper: tuple[float, tuple[float, ...]] | None = None
    maximum_row_guard = 0.0
    maximum_interval_guard = 0.0
    maximum_middle_control_speed = 0.0
    maximum_speed = 0.0
    rows = 0
    try:
        with path.open("rb") as source:
            for line_number, raw_line in enumerate(source, 1):
                digest.update(raw_line)
                payload = raw_line.split(b"#", 1)[0].strip()
                if not payload:
                    continue
                try:
                    values = tuple(float(token) for token in payload.split())
                except ValueError as error:
                    raise CampaignError(
                        f"{path}:{line_number}: non-numeric trajectory field"
                    ) from error
                if len(values) != 21:
                    raise CampaignError(
                        f"{path}:{line_number}: expected 21 columns, found {len(values)}"
                    )
                if not all(math.isfinite(value) for value in values):
                    raise CampaignError(f"{path}:{line_number}: non-finite value")
                time, state = values[0], values[1:]
                if previous is not None and not time > previous[0]:
                    raise CampaignError(f"{path}:{line_number}: times are not strictly increasing")
                if state[MASS[0]] < 0.0 or state[MASS[1]] < 0.0:
                    raise CampaignError(f"{path}:{line_number}: negative component mass")
                if not state[MASS[0]] + state[MASS[1]] > 0.0:
                    raise CampaignError(f"{path}:{line_number}: non-positive total mass")
                row = (time, state)
                if first is None:
                    first = row
                if time <= start_time:
                    lower = row
                if upper is None and time >= start_time:
                    upper = row
                for target in tracking_targets(state):
                    maximum_speed = max(maximum_speed, vector_norm(target.velocity))
                    row_guard = horizon_guard_radius(target, horizon_factor)
                    if row_guard > guard_limit:
                        raise CampaignError(
                            f"trajectory row {line_number} horizon guard {row_guard:.6g}M "
                            f"exceeds the audited {guard_limit:g}M bound"
                        )
                    maximum_row_guard = max(maximum_row_guard, row_guard)
                if previous is not None:
                    for hole in range(2):
                        guard, control_speed = interval_guard_bound(
                            previous[0], previous[1], time, state, hole, horizon_factor
                        )
                        maximum_interval_guard = max(maximum_interval_guard, guard)
                        if guard > guard_limit:
                            raise CampaignError(
                                f"trajectory interval [{previous[0]:g}, {time:g}] hole "
                                f"{hole + 1} has certified horizon guard {guard:.6g}M, "
                                f"above the audited {guard_limit:g}M bound"
                            )
                        maximum_middle_control_speed = max(
                            maximum_middle_control_speed, control_speed
                        )
                    if terms_coincident(previous[1]) and terms_coincident(state):
                        guard, _ = interval_guard_bound(
                            previous[0],
                            previous[1],
                            time,
                            state,
                            0,
                            horizon_factor,
                            merged=True,
                        )
                        maximum_interval_guard = max(maximum_interval_guard, guard)
                        if guard > guard_limit:
                            raise CampaignError(
                                f"coincident-target interval [{previous[0]:g}, {time:g}] "
                                f"has certified horizon guard {guard:.6g}M, above the "
                                f"audited {guard_limit:g}M bound"
                            )
                previous = row
                rows += 1
    except OSError as error:
        raise CampaignError(f"cannot read trajectory {path}: {error}") from error
    if rows < 2 or first is None or previous is None:
        raise CampaignError("trajectory must contain at least two data rows")
    if lower is None or upper is None:
        raise CampaignError("trajectory_time_offset lies outside the trajectory table")
    state_at_start = interpolate_state(lower, upper, start_time)
    targets_at_start = tracking_targets(state_at_start)
    if not targets_at_start:
        raise CampaignError("trajectory has no active target at simulation start")
    return TrajectorySummary(
        rows,
        first[0],
        previous[0],
        state_at_start,
        previous[1],
        targets_at_start,
        tracking_targets(previous[1]),
        max(maximum_row_guard, maximum_interval_guard),
        maximum_row_guard,
        maximum_middle_control_speed,
        maximum_speed,
        digest.hexdigest(),
    )


def first_trajectory_time(path: Path) -> float:
    try:
        with path.open("rb") as source:
            for line_number, raw_line in enumerate(source, 1):
                payload = raw_line.split(b"#", 1)[0].strip()
                if not payload:
                    continue
                try:
                    value = float(payload.split(maxsplit=1)[0])
                except ValueError as error:
                    raise CampaignError(
                        f"{path}:{line_number}: non-numeric first trajectory time"
                    ) from error
                if not math.isfinite(value):
                    raise CampaignError("first trajectory time is non-finite")
                return value
    except OSError as error:
        raise CampaignError(f"cannot read trajectory {path}: {error}") from error
    raise CampaignError("trajectory contains no data rows")


def validate_baseline(summary: TrajectorySummary, matrix: dict[str, object]) -> None:
    baseline = matrix["trajectory"]["baseline"]
    tolerance = float(baseline["normalization_tolerance"])
    targets = summary.targets_at_start
    if len(targets) != 2:
        raise CampaignError("baseline campaign must start before merger with two targets")
    total_mass = sum(target.mass for target in targets)
    expected_mass = float(baseline["initial_total_mass_M"])
    if not math.isclose(total_mass, expected_mass, rel_tol=tolerance, abs_tol=tolerance):
        raise CampaignError(f"initial total mass {total_mass:.17g} is not normalized to 1M")
    ratio = max(target.mass for target in targets) / min(target.mass for target in targets)
    if not math.isclose(
        ratio, float(baseline["mass_ratio"]), rel_tol=tolerance, abs_tol=tolerance
    ):
        raise CampaignError(f"initial mass ratio {ratio:.17g} is not the q=1 baseline")
    separation = vector_norm(
        first - second for first, second in zip(targets[0].position, targets[1].position)
    )
    if abs(separation - float(baseline["initial_separation_M"])) > float(
        baseline["separation_tolerance_M"]
    ):
        raise CampaignError(f"initial separation {separation:.17g}M is outside the 20M gate")
    for target in targets:
        if vector_norm(target.spin) > tolerance * max(1.0, target.mass):
            raise CampaignError("baseline trajectory is not initially nonspinning")
    center = tuple(
        sum(target.mass * target.position[axis] for target in targets) / total_mass
        for axis in range(3)
    )
    center_velocity = tuple(
        sum(target.mass * target.velocity[axis] for target in targets) / total_mass
        for axis in range(3)
    )
    if vector_norm(center) > 1.0e-4 or vector_norm(center_velocity) > 1.0e-4:
        raise CampaignError("baseline trajectory does not start in the center-of-mass frame")


def finite_vector(value: object, name: str) -> tuple[float, float, float]:
    if not isinstance(value, list) or len(value) != 3:
        raise CampaignError(f"provenance {name} must be a three-component array")
    vector = tuple(float(component) for component in value)
    if not all(math.isfinite(component) for component in vector):
        raise CampaignError(f"provenance {name} contains a non-finite value")
    return vector


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as source:
            for chunk in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as error:
        raise CampaignError(f"cannot hash source artifact {path}: {error}") from error
    return digest.hexdigest()


def load_provenance(
    path: Path,
    summary: TrajectorySummary,
    matrix: dict[str, object],
    source_artifact_override: Path | None,
    expected_revision_override: str | None,
) -> ProvenanceSummary:
    if not path.is_file():
        raise CampaignError(f"trajectory provenance sidecar does not exist: {path}")
    try:
        raw = path.read_bytes()
        provenance = json.loads(raw)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise CampaignError(f"cannot load provenance sidecar {path}: {error}") from error
    required_version = matrix["trajectory"]["required_provenance_schema_version"]
    if not isinstance(provenance, dict) or provenance.get("provenance_schema_version") != required_version:
        raise CampaignError("trajectory provenance must use schema version 2")
    sidecar_sha256 = hashlib.sha256(raw).hexdigest()
    if sidecar_sha256 != matrix["trajectory"]["expected_provenance_sha256"]:
        raise CampaignError("provenance sidecar SHA-256 is not the frozen campaign artifact")
    try:
        output = provenance["output"]
        model = provenance["trajectory_model"]
        remnant = provenance["remnant"]
        schema = provenance["schema"]
    except KeyError as error:
        raise CampaignError(f"provenance is missing required object {error.args[0]!r}") from error

    if not all(isinstance(value, dict) for value in (output, model, remnant, schema)):
        raise CampaignError("provenance required objects must be JSON objects")
    if output.get("sha256") != summary.sha256:
        raise CampaignError("provenance output.sha256 does not match the trajectory")
    try:
        output_rows = int(output.get("rows", -1))
    except (TypeError, ValueError) as error:
        raise CampaignError("provenance output.rows is unavailable") from error
    if output_rows != summary.rows:
        raise CampaignError("provenance output.rows does not match the trajectory")
    bounds = output.get("time_bounds")
    try:
        numeric_bounds = tuple(float(value) for value in bounds)
    except (TypeError, ValueError) as error:
        raise CampaignError("provenance output.time_bounds is unavailable") from error
    if len(numeric_bounds) != 2 or not all(math.isfinite(value) for value in numeric_bounds):
        raise CampaignError("provenance output.time_bounds is unavailable")
    if not math.isclose(numeric_bounds[0], summary.first_time, rel_tol=0.0, abs_tol=1e-12) or not math.isclose(
        numeric_bounds[1], summary.last_time, rel_tol=0.0, abs_tol=1e-12
    ):
        raise CampaignError("provenance output.time_bounds does not match the trajectory")

    labels_value = model.get("labels")
    if not isinstance(labels_value, list) or not all(
        isinstance(label, str) for label in labels_value
    ):
        raise CampaignError("provenance trajectory_model.labels is unavailable")
    labels = tuple(labels_value)
    required_labels = ("frozen-CBwaves", "local-instantaneous-4PN")
    if labels != required_labels:
        raise CampaignError("provenance trajectory_model labels are not the audited 4PN pair")
    if model.get("declaration_source") != "CLI":
        raise CampaignError("provenance trajectory_model.declaration_source must be CLI")
    source_revision = model.get("source_revision")
    revision_pattern = matrix["trajectory"]["required_source_revision_pattern"]
    if not isinstance(source_revision, str) or re.fullmatch(
        revision_pattern, source_revision.casefold()
    ) is None:
        raise CampaignError(
            "provenance trajectory_model.source_revision must be a full 40-hex commit"
        )
    matrix_revision = matrix["trajectory"]["expected_source_revision"]
    expected_revision = expected_revision_override or matrix_revision
    if expected_revision != matrix_revision:
        raise CampaignError("--expected-source-revision does not match the audited matrix")
    if source_revision.casefold() != expected_revision.casefold():
        raise CampaignError("provenance source revision does not match the audited campaign")
    source_provenance = model.get("source_provenance")
    if not isinstance(source_provenance, dict):
        raise CampaignError("provenance trajectory_model.source_provenance is unavailable")
    declared_source_path = source_provenance.get("path")
    declared_source_sha256 = source_provenance.get("sha256")
    if not isinstance(declared_source_path, str) or not declared_source_path.strip():
        raise CampaignError("provenance source artifact path is unavailable")
    if not isinstance(declared_source_sha256, str) or re.fullmatch(
        r"[0-9a-f]{64}", declared_source_sha256.casefold()
    ) is None:
        raise CampaignError("provenance source artifact SHA-256 is unavailable")
    if declared_source_sha256.casefold() != matrix["trajectory"]["expected_source_sha256"]:
        raise CampaignError("source_provenance SHA-256 is not the frozen campaign source")
    source_artifact = (
        source_artifact_override.resolve()
        if source_artifact_override is not None
        else Path(declared_source_path).resolve()
    )
    if not source_artifact.is_file():
        raise CampaignError(
            "source provenance cannot be re-hashed; provide its local copy with "
            "--source-artifact"
        )
    measured_source_sha256 = file_sha256(source_artifact)
    if measured_source_sha256 != declared_source_sha256.casefold():
        raise CampaignError("source artifact SHA-256 does not match source_provenance")

    if remnant.get("representation") != "canonical-single-term-1":
        raise CampaignError("provenance does not declare a canonical single-term remnant")
    try:
        remnant_mass = float(remnant["mass"])
    except (KeyError, TypeError, ValueError) as error:
        raise CampaignError("provenance remnant.mass is unavailable") from error
    if not math.isfinite(remnant_mass) or remnant_mass <= 0.0:
        raise CampaignError("provenance remnant.mass must be finite and positive")
    remnant_a = finite_vector(remnant.get("a"), "remnant.a")
    remnant_chi = finite_vector(remnant.get("chi"), "remnant.chi")
    remnant_kick = finite_vector(remnant.get("kick"), "remnant.kick")
    if remnant.get("a_to_chi_relation") != "a_length = Mf * chi; chi = a_length / Mf":
        raise CampaignError("provenance remnant a/chi convention is unavailable")
    for a_value, chi_value in zip(remnant_a, remnant_chi):
        if not math.isclose(a_value, remnant_mass * chi_value, rel_tol=1e-12, abs_tol=1e-12):
            raise CampaignError("provenance remnant a and chi fields are inconsistent")

    if schema.get("version") != 1 or schema.get("columns") != matrix["trajectory"]["columns"]:
        raise CampaignError("provenance trajectory schema does not match the 21-column matrix")
    if len(summary.targets_at_end) != 1 or not summary.state_at_end[MASS[0]] > 0.0 or not close(
        summary.state_at_end[MASS[1]], 0.0
    ):
        raise CampaignError("trajectory does not end in canonical single-term-1 form")
    final_target = summary.targets_at_end[0]
    if not math.isclose(final_target.mass, remnant_mass, rel_tol=1e-12, abs_tol=1e-12):
        raise CampaignError("provenance remnant.mass does not match the final trajectory row")
    for actual, declared, name in (
        (final_target.spin, remnant_a, "spin a"),
        (final_target.velocity, remnant_kick, "kick"),
    ):
        if any(
            not math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
            for a, b in zip(actual, declared)
        ):
            raise CampaignError(f"provenance remnant {name} does not match the final row")
    return ProvenanceSummary(
        path.resolve(),
        sidecar_sha256,
        labels,
        source_revision,
        source_artifact,
        measured_source_sha256,
        remnant_mass,
    )


def projected_output_gib(matrix: dict[str, object], tier: dict[str, object], tlim: float) -> float:
    outputs = matrix["common"]["outputs"]
    blocks = int(tier["topology_estimate"]["campaign_meshblock_gate"])
    sizes = outputs["estimated_bytes_per_meshblock"]
    restart_cadence = min(
        float(outputs["restart_dt_M"]),
        float(matrix["common"]["resource_model"]["streaming_segment_span_M"]),
    )
    dump_specs = (
        (float(outputs["full_state_dt_M"]), int(sizes["full_state_mhd_w_bcc"])),
        (float(outputs["gr_diagnostics_dt_M"]), int(sizes["gr_diagnostics"])),
        (float(outputs["divb_dt_M"]), int(sizes["divb"])),
        (restart_cadence, int(sizes["restart_double_precision_dynadm"])),
    )
    total = 0
    for cadence, bytes_per_block in dump_specs:
        count = math.floor(tlim / cadence) + 2
        total += count * blocks * bytes_per_block
    return total / 2**30


def validate_resources(
    args: argparse.Namespace,
    matrix: dict[str, object],
    tier: dict[str, object],
    tlim: float,
) -> tuple[int, int, int, float, float]:
    if args.gpus <= 0 or not math.isfinite(args.gpu_memory_gib) or args.gpu_memory_gib <= 0.0:
        raise CampaignError("GPU count and memory must be positive")
    if not math.isfinite(args.scratch_gib) or args.scratch_gib <= 0.0:
        raise CampaignError("scratch capacity must be positive")
    resource = matrix["common"]["resource_model"]
    per_block_mib = float(resource["empirical_peak_MiB_per_meshblock"])
    usable_mib = args.gpu_memory_gib * 1024.0 * float(resource["gpu_memory_usable_fraction"])
    capacity_per_rank = math.floor(usable_mib / per_block_mib)
    if capacity_per_rank <= 0:
        raise CampaignError("declared GPU memory cannot hold one empirical MeshBlock allocation")
    gate_blocks = int(tier["topology_estimate"]["campaign_meshblock_gate"])
    minimum_gpus = math.ceil(gate_blocks / capacity_per_rank)
    if args.gpus < minimum_gpus:
        raise CampaignError(
            f"resource gate requires at least {minimum_gpus} x {args.gpu_memory_gib:g}GiB "
            f"GPUs for {gate_blocks} gated MeshBlocks"
        )
    qualified_ranks = int(
        tier["topology_estimate"]["rank_partition_qualification"]["mpi_ranks"]
    )
    if args.gpus != qualified_ranks:
        raise CampaignError(
            f"tier {args.tier} is partition-qualified only for {qualified_ranks} MPI "
            f"ranks/GPUs, not {args.gpus}; run a real GID-ordered topology and peak-memory "
            f"qualification, then update the campaign matrix before using another count"
        )
    minimum_qualified_memory = float(
        tier["resource_gate"]["minimum_gpu_memory_GiB_at_qualified_count"]
    )
    if args.gpu_memory_gib < minimum_qualified_memory:
        raise CampaignError(
            f"tier {args.tier} requires at least {minimum_qualified_memory:g} GiB/GPU "
            "at its partition-qualified rank count"
        )
    rounding = int(resource["max_nmb_rounding"])
    gate_based_max_nmb = math.ceil((gate_blocks / args.gpus) / rounding) * rounding
    partition_floor = int(
        tier["topology_estimate"]["rank_partition_qualification"][
            "max_nmb_per_rank_hard_floor"
        ]
    )
    max_nmb_per_rank = max(gate_based_max_nmb, partition_floor)
    if max_nmb_per_rank > capacity_per_rank:
        raise CampaignError(
            f"per-rank MeshBlock allocation {max_nmb_per_rank} exceeds the empirical "
            f"GPU memory capacity {capacity_per_rank}"
        )

    projected_gib = projected_output_gib(matrix, tier, tlim)
    if args.streaming_drain:
        minimum_rate = float(resource["streaming_minimum_drain_MiB_per_second"])
        if args.drain_mib_s is None or args.drain_mib_s < minimum_rate:
            raise CampaignError(
                f"streaming mode requires a declared verified drain >= {minimum_rate:g} MiB/s"
            )
        required_scratch = float(tier["resource_gate"]["minimum_streaming_scratch_GiB"])
    else:
        if args.drain_mib_s is not None:
            raise CampaignError("--drain-mib-s requires --streaming-drain")
        required_scratch = math.ceil(1.25 * projected_gib)
    if args.scratch_gib < required_scratch:
        mode = "streaming" if args.streaming_drain else "undrained"
        raise CampaignError(
            f"{mode} storage gate requires {required_scratch:g} GiB; "
            f"only {args.scratch_gib:g} GiB was declared"
        )
    return minimum_gpus, max_nmb_per_rank, gate_based_max_nmb, projected_gib, required_scratch


def format_real(value: float) -> str:
    return f"{value:.17g}"


def float32(value: float) -> float:
    try:
        return struct.unpack("!f", struct.pack("!f", value))[0]
    except OverflowError as error:
        raise CampaignError("single-precision campaign time overflows float32") from error


def validate_precision(
    args: argparse.Namespace,
    matrix: dict[str, object],
    first_table_time: float,
    last_table_time: float,
    metric_fd_step: float,
) -> None:
    publication_precision = matrix["common"]["time"]["publication_precision"]
    if args.real_precision == "double":
        if args.allow_single_precision_systematic:
            raise CampaignError("--allow-single-precision-systematic requires --real-precision single")
        if publication_precision != "double":
            raise CampaignError("matrix publication precision is not double")
        represent = float
    else:
        if not args.allow_single_precision_systematic:
            raise CampaignError(
                "single precision is a non-certificate systematic; add "
                "--allow-single-precision-systematic explicitly"
            )
        if args.metric_fd_step is None:
            raise CampaignError("single-precision systematic requires an explicit --metric-fd-step")
        represent = float32
    step = represent(metric_fd_step)
    for table_time in (first_table_time, last_table_time):
        center = represent(table_time)
        earlier = represent(center - step)
        later = represent(center + step)
        if not earlier < center < later:
            raise CampaignError(
                f"metric_fd_step={metric_fd_step:g} is not representable around "
                f"table time {table_time:g} in {args.real_precision} precision"
            )


def refined_regions(
    matrix: dict[str, object], tier: dict[str, object], targets: tuple[Target, ...]
) -> str:
    mesh = matrix["common"]["mesh"]
    amr = matrix["common"]["amr"]
    half_extent = float(mesh["domain_half_extent_M"])
    floor_radius = float(amr["disk_floor_radius_M"])
    floor_level = int(amr["disk_floor_physical_level"])
    finest_level = int(tier["finest_physical_level"])
    finest_nominal = float(tier["shell_radii_M"][str(finest_level)])
    horizon_factor = float(amr["horizon_factor"])
    regions = [
        (
            floor_level,
            (-floor_radius, -floor_radius, -floor_radius),
            (floor_radius, floor_radius, floor_radius),
            "preseed the persistent central disk floor",
        )
    ]
    for index, target in enumerate(targets, 1):
        radius = max(finest_nominal, horizon_guard_radius(target, horizon_factor))
        lower = tuple(value - radius for value in target.position)
        upper = tuple(value + radius for value in target.position)
        if any(value <= -half_extent for value in lower) or any(
            value >= half_extent for value in upper
        ):
            raise CampaignError(f"initial target {index} seed touches the root boundary")
        regions.append(
            (finest_level, lower, upper, f"preseed finest moving-hole target {index}")
        )
    blocks = []
    for number, (level, lower, upper, description) in enumerate(regions, 1):
        blocks.extend(
            [
                f"# {description}",
                f"<refined_region{number}>",
                f"level = {level}",
                f"x1min = {format_real(lower[0])}",
                f"x1max = {format_real(upper[0])}",
                f"x2min = {format_real(lower[1])}",
                f"x2max = {format_real(upper[1])}",
                f"x3min = {format_real(lower[2])}",
                f"x3max = {format_real(upper[2])}",
                "",
            ]
        )
    return "\n".join(blocks).rstrip()


def render_input(
    args: argparse.Namespace,
    matrix: dict[str, object],
    tier: dict[str, object],
    summary: TrajectorySummary,
    provenance: ProvenanceSummary,
    runtime_trajectory: str,
    tlim: float,
    metric_fd_step: float,
    cfl_number: float,
    root_dt_max: float,
    max_nmb_per_rank: int,
    projected_gib: float,
    required_scratch: float,
) -> str:
    common = matrix["common"]
    mesh = common["mesh"]
    time = common["time"]
    amr = common["amr"]
    outputs = common["outputs"]
    half_extent = float(mesh["domain_half_extent_M"])
    finest_level = int(tier["finest_physical_level"])
    basename = args.basename or f"effective_bbh_4pn_{args.tier}"
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", basename):
        raise CampaignError("basename may contain only letters, digits, dot, underscore, and dash")
    seed_blocks = refined_regions(matrix, tier, summary.targets_at_start)
    radii = "\n".join(
        f"refinement_radius_level_{level} = "
        f"{format_real(float(tier['shell_radii_M'][str(level)]))}"
        for level in range(1, finest_level + 1)
    )
    return f"""# GENERATED FILE: effective-BBH 4PN multilevel-AMR convergence tier {args.tier}.
# Generator: inputs/dyngr/generate_effective_bbh_4pn_campaign.py
# Matrix: {matrix['campaign_id']} (schema {matrix['schema_version']})
# Trajectory SHA-256: {summary.sha256}
# Provenance sidecar: {provenance.path}
# Provenance SHA-256: {provenance.sha256}
# Audited trajectory model: {', '.join(provenance.model_labels)}
# Trajectory source revision: {provenance.source_revision}
# Source artifact SHA-256: {provenance.source_sha256}
# Source artifact re-hashed at: {provenance.source_artifact}
# REQUIRED BUILD PRECISION: {args.real_precision}
# RESOURCE CONTRACT: {args.gpus} MPI ranks/GPUs, {format_real(args.gpu_memory_gib)} GiB/GPU,
# max_nmb_per_rank={max_nmb_per_rank}, scratch={format_real(args.scratch_gib)} GiB.
# Estimated undrained output={projected_gib:.1f} GiB; active scratch gate={required_scratch:g} GiB.
# Do not change one tier independently: regenerate all L/M/H inputs from one matrix/table.

<comment>
problem = q1 nonspinning 4PN-to-remnant effective-BBH torus convergence tier {args.tier}

<job>
basename = {basename}

<mesh>
nghost = {int(mesh['nghost'])}
nx1 = {int(mesh['root_cells_per_dimension'])}
x1min = {-half_extent:g}
x1max = {half_extent:g}
ix1_bc = outflow
ox1_bc = outflow
nx2 = {int(mesh['root_cells_per_dimension'])}
x2min = {-half_extent:g}
x2max = {half_extent:g}
ix2_bc = outflow
ox2_bc = outflow
nx3 = {int(mesh['root_cells_per_dimension'])}
x3min = {-half_extent:g}
x3max = {half_extent:g}
ix3_bc = outflow
ox3_bc = outflow

<meshblock>
nx1 = {int(mesh['meshblock_cells_per_dimension'])}
nx2 = {int(mesh['meshblock_cells_per_dimension'])}
nx3 = {int(mesh['meshblock_cells_per_dimension'])}

<mesh_refinement>
refinement = adaptive
num_levels = {int(tier['num_levels'])}
max_nmb_per_rank = {max_nmb_per_rank}
ncycle_check = 1
refinement_interval = 1
prolong_primitives = false

<amr_criterion0>
method = user

{seed_blocks}

<time>
evolution = dynamic
integrator = {time['integrator']}
subcycling = level
root_dt_max = {format_real(root_dt_max)}
cfl_number = {format_real(cfl_number)}
nlim = -1
tlim = {format_real(tlim)}
ndiag = 1

<coord>
minkowski = false
excise = true
excision_scheme = lapse
excise_lapse = 0.25
dexcise = 1.0e-10
pexcise = 3.333333333333e-13

<adm>
dynamic = true

<mhd>
eos = ideal
dyn_eos = ideal
dyn_error = reset_floor
gamma = 1.3333333333333333
reconstruct = wenoz
rsolver = hlle
dfloor = 1.0e-10
tfloor = 3.333333333333e-13
dthreshold = 1.01
dyn_scratch = 1
fofc = true
enforce_maximum = false

<problem>
pgen_name = dynbbh
trajectory_mode = table
trajectory_file = {runtime_trajectory}
trajectory_time_offset = {format_real(args.trajectory_time_offset)}
separation = 20.0
mass_ratio = 1.0
chi1 = 0.0
chi2 = 0.0
theta1 = 0.0
theta2 = 0.0
phi1 = 0.0
phi2 = 0.0
metric_fd_step = {format_real(metric_fd_step)}
mass_scale1 = 1.0
mass_scale2 = 1.0
spin_buffer1 = 0.05
spin_buffer2 = 0.05
singularity_floor = 1.0e-3
amr_condition = track
alpha_threshold = 0.6
refinement_radius = 1.0
refinement_radius_ratio = 2.0
refinement_hysteresis = {format_real(float(amr['hysteresis']))}
refinement_horizon_factor = {format_real(float(amr['horizon_factor']))}
{radii}
refinement_floor_radius = {format_real(float(amr['disk_floor_radius_M']))}
refinement_floor_level = {int(amr['disk_floor_physical_level'])}
refinement_floor_center1 = {format_real(float(amr['disk_floor_center_M'][0]))}
refinement_floor_center2 = {format_real(float(amr['disk_floor_center_M'][1]))}
refinement_floor_center3 = {format_real(float(amr['disk_floor_center_M'][2]))}

initial_data = reference_fm_torus
torus_reference_mass = 1.0
torus_reference_center1 = 0.0
torus_reference_center2 = 0.0
torus_reference_center3 = 0.0
torus_reference_velocity1 = 0.0
torus_reference_velocity2 = 0.0
torus_reference_velocity3 = 0.0
torus_r_edge = 18.0
torus_r_peak = 29.0
torus_rho_max = 1.0
torus_rho_min = 1.0e-5
torus_rho_pow = -1.5
torus_pgas_min = 3.333333333333e-8
torus_pgas_pow = -2.5
torus_rho_atm = 1.0e-10
torus_temp_atm = 3.333333333333e-13
torus_pert_amp = 0.02
torus_seed = 1
torus_magnetic_field = single_loop
torus_potential_cutoff = 0.2
torus_mag_norm = density_weighted_beta
torus_mag_target = 100.0
torus_min_grid_peak_fraction = 0.9
torus_min_magnetic_cells = 64
torus_require_full_domain = true

<output1>
file_type = hst
data_format = %.17e
dt = {format_real(float(outputs['history_dt_M']))}

<output2>
file_type = bin
variable = mhd_w_bcc
dt = {format_real(float(outputs['full_state_dt_M']))}
ghost_zones = false

<output3>
file_type = bin
variable = mhd_divb
dt = {format_real(float(outputs['divb_dt_M']))}
ghost_zones = false

<output4>
file_type = rst
dt = {format_real(float(outputs['restart_dt_M']))}

<output5>
file_type = bin
variable = mhd_gr_diagnostics
dt = {format_real(float(outputs['gr_diagnostics_dt_M']))}
ghost_zones = false
"""


def publish(path: Path, content: str, force: bool) -> None:
    if not path.parent.is_dir():
        raise CampaignError(f"output directory does not exist: {path.parent}")
    if path.exists() and not force:
        raise CampaignError(f"refusing to replace existing output {path}; use --force")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as target:
            target.write(content)
            target.flush()
            os.fsync(target.fileno())
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def main() -> int:
    args = parse_args()
    try:
        matrix = load_matrix(args.matrix)
        tier = matrix["tiers"][args.tier]
        metric_fd_step = (
            float(matrix["common"]["time"]["metric_fd_step_M"])
            if args.metric_fd_step is None
            else args.metric_fd_step
        )
        if not math.isfinite(metric_fd_step) or metric_fd_step <= 0.0:
            raise CampaignError("metric FD step must be finite and positive")
        if args.trajectory_time_offset is None:
            args.trajectory_time_offset = math.nextafter(
                first_trajectory_time(args.trajectory) + metric_fd_step, math.inf
            )
        if not math.isfinite(args.trajectory_time_offset):
            raise CampaignError("trajectory time offset must be finite")
        baseline_cfl = float(matrix["common"]["time"]["cfl_number"])
        cfl_number = baseline_cfl if args.cfl_number is None else args.cfl_number
        if not math.isfinite(cfl_number) or not 0.0 < cfl_number <= baseline_cfl:
            raise CampaignError(f"CFL must lie in (0, {baseline_cfl:g}]")
        root_dt_max = cfl_number * float(matrix["common"]["mesh"]["root_dx_M"])
        if not math.isfinite(root_dt_max) or root_dt_max <= 0.0:
            raise CampaignError("derived root_dt_max must be finite and positive")
        summary = inspect_trajectory(
            args.trajectory,
            args.trajectory_time_offset,
            float(matrix["common"]["amr"]["horizon_factor"]),
            float(matrix["trajectory"]["maximum_horizon_guard_radius_M"]),
        )
        if summary.sha256 != matrix["trajectory"]["expected_ascii_sha256"]:
            raise CampaignError("trajectory SHA-256 is not the frozen campaign artifact")
        validate_baseline(summary, matrix)
        provenance = load_provenance(
            args.trajectory_provenance,
            summary,
            matrix,
            args.source_artifact,
            args.expected_source_revision,
        )
        horizon_limit = float(matrix["trajectory"]["maximum_horizon_guard_radius_M"])
        if summary.maximum_guard_radius > horizon_limit:
            raise CampaignError(
                f"trajectory horizon guard reaches {summary.maximum_guard_radius:.6g}M, "
                f"above the audited {horizon_limit:g}M resource bound"
            )
        tlim = (
            summary.last_time - args.trajectory_time_offset - metric_fd_step
            if args.tlim is None
            else args.tlim
        )
        if not math.isfinite(tlim) or tlim <= 0.0:
            raise CampaignError("tlim must be finite and positive")
        validate_precision(
            args,
            matrix,
            args.trajectory_time_offset,
            args.trajectory_time_offset + tlim,
            metric_fd_step,
        )
        required_first = args.trajectory_time_offset - metric_fd_step
        required_last = args.trajectory_time_offset + tlim + metric_fd_step
        tolerance = 32.0 * sys.float_info.epsilon * max(1.0, abs(required_last))
        if summary.first_time > required_first + tolerance:
            raise CampaignError(
                f"trajectory begins at {summary.first_time:g}, after required {required_first:g}"
            )
        if summary.last_time < required_last - tolerance:
            raise CampaignError(
                f"trajectory ends at {summary.last_time:g}, before required {required_last:g}"
            )
        (
            minimum_gpus,
            max_nmb,
            gate_based_max_nmb,
            projected_gib,
            required_scratch,
        ) = validate_resources(args, matrix, tier, tlim)
        runtime_trajectory = args.runtime_trajectory or str(args.trajectory.resolve())
        if any(character.isspace() for character in runtime_trajectory) or "#" in runtime_trajectory:
            raise CampaignError("runtime trajectory path cannot contain whitespace or '#'")
        report = {
            "tier": args.tier,
            "trajectory_sha256": summary.sha256,
            "provenance_path": str(provenance.path),
            "provenance_sha256": provenance.sha256,
            "trajectory_model_labels": provenance.model_labels,
            "trajectory_source_revision": provenance.source_revision,
            "source_artifact": str(provenance.source_artifact),
            "source_artifact_sha256": provenance.source_sha256,
            "baseline_validation": "q1-nonspinning-canonical-single-remnant",
            "trajectory_rows": summary.rows,
            "trajectory_time_range": [summary.first_time, summary.last_time],
            "trajectory_time_offset_M": args.trajectory_time_offset,
            "tlim_M": tlim,
            "cfl_number": cfl_number,
            "root_dt_max_M": root_dt_max,
            "real_precision": args.real_precision,
            "maximum_horizon_guard_radius_M": summary.maximum_guard_radius,
            "maximum_row_horizon_guard_radius_M": summary.maximum_row_guard_radius,
            "maximum_hermite_middle_control_speed": summary.maximum_middle_control_speed,
            "finest_dx_M": tier["finest_dx_M"],
            "campaign_meshblock_gate": tier["topology_estimate"]["campaign_meshblock_gate"],
            "rootcap_runtime_observation": tier["topology_estimate"].get(
                "rootcap_runtime_observation"
            ),
            "aggregate_memory_lower_bound_gpus": minimum_gpus,
            "qualified_launch_gpus": tier["topology_estimate"][
                "rank_partition_qualification"
            ]["mpi_ranks"],
            "declared_gpus": args.gpus,
            "declared_gpu_memory_GiB": args.gpu_memory_gib,
            "runtime_qualification_status": tier["resource_gate"][
                "runtime_qualification_status"
            ],
            "max_nmb_per_rank": max_nmb,
            "gate_based_max_nmb_per_rank": gate_based_max_nmb,
            "rank_partition_qualification": {
                "artifact": tier["topology_estimate"]["rank_partition_qualification"][
                    "topology_artifact"
                ],
                "mesh_structure_sha256": tier["topology_estimate"][
                    "rank_partition_qualification"
                ]["mesh_structure_sha256"],
                "mpi_ranks": tier["topology_estimate"]["rank_partition_qualification"][
                    "mpi_ranks"
                ],
                "meshblocks_per_rank_range": tier["topology_estimate"][
                    "rank_partition_qualification"
                ]["meshblocks_per_rank_range"],
                "weighted_cost_per_rank_range": tier["topology_estimate"][
                    "rank_partition_qualification"
                ]["weighted_cost_per_rank_range"],
                "max_nmb_per_rank_hard_floor": tier["topology_estimate"][
                    "rank_partition_qualification"
                ]["max_nmb_per_rank_hard_floor"],
                "declared_rank_count_matches_qualification": args.gpus
                == int(
                    tier["topology_estimate"]["rank_partition_qualification"][
                        "mpi_ranks"
                    ]
                ),
            },
            "projected_undrained_output_GiB": projected_gib,
            "required_scratch_GiB": required_scratch,
        }
        if not args.validate_only:
            if args.output is None:
                raise CampaignError("--output is required unless --validate-only is used")
            content = render_input(
                args,
                matrix,
                tier,
                summary,
                provenance,
                runtime_trajectory,
                tlim,
                metric_fd_step,
                cfl_number,
                root_dt_max,
                max_nmb,
                projected_gib,
                required_scratch,
            )
            publish(args.output, content, args.force)
            report["output"] = str(args.output)
        print(json.dumps(report, indent=2, sort_keys=True))
    except CampaignError as error:
        print(f"campaign validation failed: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
