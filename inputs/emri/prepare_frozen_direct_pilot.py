#!/usr/bin/env python3
"""Prepare a reduced direct-AMR AthenaK pilot from a frozen campaign case."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import build_frozen_taylor_campaign as frozen
import extract_static_taylor_worldtube as static
import plan_bhl_hierarchy as hierarchy


CLASSIFICATION = "athenak-emri-frozen-direct-pilot-v2"


def _case(campaign: dict[str, object], identifier: str) -> dict[str, object]:
    if campaign.get("classification") != frozen.CLASSIFICATION:
        raise ValueError("input is not a frozen Taylor campaign")
    matches = [case for case in campaign["cases"] if case.get("id") == identifier]
    if len(matches) != 1:
        raise ValueError(f"case {identifier!r} was not found exactly once")
    selected = matches[0]
    if not selected["assessment"]["local_taylor_model_passed"]:
        raise ValueError("selected case failed its local Taylor-model gates")
    return selected


def inflow_boundaries(spatial_four_velocity: list[float]) -> dict[str, str]:
    """Select fixed-wind faces from the signs of source-frame velocity."""

    if len(spatial_four_velocity) != 3 or not all(
        math.isfinite(float(value)) for value in spatial_four_velocity
    ):
        raise ValueError("spatial four-velocity must contain three finite values")
    result = {}
    for axis, value in enumerate(spatial_four_velocity, start=1):
        result[f"ix{axis}_bc"] = "user" if value > 0.0 else "outflow"
        result[f"ox{axis}_bc"] = "user" if value < 0.0 else "outflow"
    return result


def build_pilot(
    campaign: dict[str, object],
    identifier: str,
    *,
    cells_per_secondary_mass: float,
    outer_cells_per_capture_radius: float,
    flow_crossings: float,
    meshblock_cells: int,
    calibration_cycles: int,
    finest_refinement_radius: float,
    maximum_meshblocks_per_rank: int | None = None,
    parallel_ranks: int = 1,
) -> dict[str, object]:
    case = _case(campaign, identifier)
    if meshblock_cells < 4 or meshblock_cells % 2 != 0:
        raise ValueError("meshblock cells must be an even integer of at least four")
    if calibration_cycles < 1:
        raise ValueError("calibration cycles must be positive")
    if parallel_ranks < 1:
        raise ValueError("parallel ranks must be positive")
    if not math.isfinite(finest_refinement_radius) or finest_refinement_radius <= 0.0:
        raise ValueError("finest refinement radius must be finite and positive")
    reference = case["profiles"][0]
    parameters = reference["parameters"]
    thermodynamics = case["state_thermodynamics"]
    gamma = thermodynamics.get("adiabatic_index")
    if gamma is None:
        raise ValueError("selected case has no explicit source adiabatic index")
    settings = hierarchy.PlannerSettings(
        cells_per_secondary_mass=cells_per_secondary_mass,
        outer_cells_per_capture_radius=outer_cells_per_capture_radius,
        flow_crossings=flow_crossings,
        root_meshblock_cells_per_axis=meshblock_cells,
        max_global_timestep_zone_updates=10**18,
    ).validated()
    primary = campaign["primary"]
    mass_ratio = float(campaign["mass_ratio"])
    local_primary_mass = float(primary["mass_global"]) / mass_ratio
    local_orbital_radius = (
        float(case["orbit"]["boyer_lindquist_radius"]) / mass_ratio
    )
    velocity = [float(parameters[f"u{i}"]) for i in range(1, 4)]
    magnetic = [float(parameters[f"b{i}"]) for i in range(1, 4)]
    plan = hierarchy.build_plan(
        secondary_mass=1.0,
        primary_mass=local_primary_mass,
        orbital_radius=local_orbital_radius,
        rho=float(parameters["rho0"]),
        pgas=float(parameters["pgas0"]),
        gamma=float(gamma),
        spatial_four_velocity=velocity,
        magnetic_field=magnetic,
        coherence_scales=case["coherence_scales_local"],
        settings=settings,
    )
    direct = plan["costs"]["direct"]
    capture_radius = float(plan["scales"]["capture_radius_factor_two_for_cost"])
    lower = [
        capture_radius * float(value)
        for value in direct["domain_lower_in_capture_radii"]
    ]
    upper = [
        capture_radius * float(value)
        for value in direct["domain_upper_in_capture_radii"]
    ]
    dimensions = [int(value) for value in direct["base_grid_dimensions"]]
    levels = int(direct["nested_refinement_levels"])
    if levels < 1:
        raise ValueError("reduced pilot unexpectedly requires no refinement")
    refinement_radii = {
        str(level): finest_refinement_radius * 2.0 ** (levels - level)
        for level in range(1, levels + 1)
    }
    block_cells = meshblock_cells**3
    estimated_blocks = math.ceil(
        1.25 * int(direct["optimistic_nested_resident_cells"]) / block_cells
    )
    root_blocks = math.prod(value // meshblock_cells for value in dimensions)
    configured_blocks = (
        math.ceil(estimated_blocks / parallel_ranks)
        if maximum_meshblocks_per_rank is None
        else int(maximum_meshblocks_per_rank)
    )
    minimum_initial_blocks_per_rank = math.ceil(root_blocks / parallel_ranks)
    if configured_blocks < minimum_initial_blocks_per_rank:
        raise ValueError(
            "maximum meshblocks per rank cannot hold a balanced initial root grid: "
            f"{configured_blocks} < {minimum_initial_blocks_per_rank}"
        )
    boundaries = inflow_boundaries(velocity)
    root_block_widths = [
        (upper[axis] - lower[axis]) * meshblock_cells / dimensions[axis]
        for axis in range(3)
    ]
    user_boundary_clearances = {}
    for axis in range(3):
        inner = f"ix{axis + 1}_bc"
        outer = f"ox{axis + 1}_bc"
        if boundaries[inner] == "user":
            user_boundary_clearances[inner] = (
                abs(lower[axis]) - root_block_widths[axis]
            )
        if boundaries[outer] == "user":
            user_boundary_clearances[outer] = (
                abs(upper[axis]) - root_block_widths[axis]
            )
    coarsest_refinement_radius = refinement_radii["1"]
    minimum_user_clearance = min(user_boundary_clearances.values())
    if not coarsest_refinement_radius < minimum_user_clearance:
        raise ValueError(
            "the level-one refinement region can touch a user physical boundary; "
            "reduce the finest refinement radius or enlarge the domain"
        )
    overrides = [
        "job/basename=emri_frozen_direct_calibration",
        f"time/nlim={calibration_cycles}",
        "time/tlim=1e30",
        "time/subcycling=level",
        "time/root_dt_max=100",
        "adm/dynamic=true",
        "mesh_refinement/refinement=adaptive",
        f"mesh_refinement/num_levels={levels + 1}",
        f"mesh_refinement/max_nmb_per_rank={configured_blocks}",
        "mesh_refinement/ncycle_check=1",
        "mesh_refinement/refinement_interval=1",
        f"mesh/nx1={dimensions[0]}",
        f"mesh/nx2={dimensions[1]}",
        f"mesh/nx3={dimensions[2]}",
        f"mesh/x1min={lower[0]:.17g}",
        f"mesh/x1max={upper[0]:.17g}",
        f"mesh/x2min={lower[1]:.17g}",
        f"mesh/x2max={upper[1]:.17g}",
        f"mesh/x3min={lower[2]:.17g}",
        f"mesh/x3max={upper[2]:.17g}",
        f"meshblock/nx1={meshblock_cells}",
        f"meshblock/nx2={meshblock_cells}",
        f"meshblock/nx3={meshblock_cells}",
    ]
    overrides.extend(
        f"mesh/{name}={value}" for name, value in boundaries.items()
    )
    overrides.extend(
        [
            "problem/user_hist=false",
            "problem/background_mode=full",
            f"problem/primary_mass={local_primary_mass:.17g}",
            "problem/secondary_mass=1",
            f"problem/primary_chi={float(primary['dimensionless_spin']):.17g}",
            f"problem/orbital_radius={local_orbital_radius:.17g}",
            f"problem/orbit_direction={int(primary['orbit_direction'])}",
            f"mhd/gamma={float(gamma):.17g}",
            "mhd/reconstruct=ppm4",
            "output1/dt=0",
            "output2/dt=0",
            "output3/dt=0",
            "output4/dt=0",
            "output5/dt=0",
        ]
    )
    overrides.extend(
        f"problem/refinement_radius_level_{level}={radius:.17g}"
        for level, radius in ((int(key), value) for key, value in refinement_radii.items())
    )
    overrides.extend(
        f"problem/{name}={float(parameters[name]):.17g}"
        for name in static.PROFILE_PARAMETER_ORDER
    )
    return {
        "classification": CLASSIFICATION,
        "case_id": identifier,
        "purpose": "GPU throughput and memory calibration; not a settled-flow run",
        "source_state_sha256": case["state_sha256"],
        "settings": settings.__dict__,
        "bhl_plan": plan,
        "mesh": {
            "lower": lower,
            "upper": upper,
            "root_dimensions": dimensions,
            "meshblock_dimensions": [meshblock_cells] * 3,
            "physical_refinement_levels": levels,
            "refinement_radii": refinement_radii,
            "estimated_meshblocks_for_budget": estimated_blocks,
            "estimated_meshblocks_per_rank_for_budget": math.ceil(
                estimated_blocks / parallel_ranks
            ),
            "maximum_meshblocks_per_rank": configured_blocks,
            "parallel_ranks": parallel_ranks,
            "capacity_limited_calibration": (
                configured_blocks * parallel_ranks < estimated_blocks
            ),
            "boundary_conditions": boundaries,
            "root_meshblock_physical_widths": root_block_widths,
            "user_boundary_refinement_clearances": user_boundary_clearances,
            "minimum_user_boundary_refinement_clearance": minimum_user_clearance,
        },
        "calibration_cycles": calibration_cycles,
        "numerical_method": {
            "integrator": "rk2",
            "reconstruction": "ppm4",
            "riemann_solver": "hlle",
            "first_order_flux_correction": True,
            "subcycling": "level",
        },
        "athena_input": "inputs/emri/emri_windtunnel_smoke.athinput",
        "athena_overrides": overrides,
        "production_warning": (
            "force history and all field outputs are disabled; re-enable them only after "
            "the measured throughput, memory, AMR topology, and CT state pass"
        ),
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cells-per-secondary-mass", type=float, default=8.0)
    parser.add_argument("--outer-cells-per-capture-radius", type=float, default=8.0)
    parser.add_argument("--flow-crossings", type=float, default=2.0)
    parser.add_argument("--meshblock-cells", type=int, default=8)
    parser.add_argument("--calibration-cycles", type=int, default=15)
    parser.add_argument("--finest-refinement-radius", type=float, default=3.0)
    parser.add_argument("--maximum-meshblocks-per-rank", type=int)
    parser.add_argument("--parallel-ranks", type=int, default=1)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    campaign_path = arguments.campaign.expanduser().resolve(strict=True)
    campaign = json.loads(campaign_path.read_text(encoding="utf-8"))
    pilot = build_pilot(
        campaign,
        arguments.case,
        cells_per_secondary_mass=arguments.cells_per_secondary_mass,
        outer_cells_per_capture_radius=arguments.outer_cells_per_capture_radius,
        flow_crossings=arguments.flow_crossings,
        meshblock_cells=arguments.meshblock_cells,
        calibration_cycles=arguments.calibration_cycles,
        finest_refinement_radius=arguments.finest_refinement_radius,
        maximum_meshblocks_per_rank=arguments.maximum_meshblocks_per_rank,
        parallel_ranks=arguments.parallel_ranks,
    )
    output = arguments.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(pilot, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(output)
    direct = pilot["bhl_plan"]["costs"]["direct"]
    print(
        f"root={direct['base_grid_dimensions']} levels={direct['nested_refinement_levels']} "
        f"ideal_updates={direct['ideal_level_subcycled_zone_updates']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
