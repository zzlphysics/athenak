#!/usr/bin/env python3
"""Plan a computationally affordable EMRI BHL hierarchy.

This is a scale and cost auditor, not a replacement for an RMHD characteristic
calculation.  It turns one locally measured upstream state into conservative capture,
environment, runtime, and outer--inner matching estimates.  The emitted JSON also
records the data-transfer and force-partition contract needed by a later spatial
worldtube implementation.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Iterable


CLASSIFICATION = "athenak-emri-bhl-hierarchy-plan-v1"


def _positive(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _finite_vector(values: Iterable[float], name: str) -> tuple[float, float, float]:
    result = tuple(float(value) for value in values)
    if len(result) != 3 or not all(math.isfinite(value) for value in result):
        raise ValueError(f"{name} must contain three finite values")
    return result  # type: ignore[return-value]


def _dot(left: Iterable[float], right: Iterable[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def _ceil_ratio(value: float) -> int:
    return max(1, int(math.ceil(value - 64.0 * math.ulp(max(1.0, value)))))


def _refinement_levels(coarse_dx: float, target_dx: float) -> int:
    if coarse_dx <= target_dx:
        return 0
    return int(math.ceil(math.log2(coarse_dx / target_dx)))


@dataclass(frozen=True)
class PlannerSettings:
    cells_per_secondary_mass: float = 16.0
    outer_cells_per_capture_radius: float = 16.0
    cells_per_outer_sink_radius: float = 8.0
    upstream_capture_radii: float = 4.0
    downstream_capture_radii: float = 8.0
    transverse_capture_radii: float = 4.0
    flow_crossings: float = 5.0
    cfl: float = 0.3
    inner_asymptotic_radius_in_m: float = 20.0
    overlap_fraction: float = 0.25
    outer_sink_fraction_of_match: float = 0.5
    environment_validity_fraction: float = 0.3
    inner_window_crossings: float = 5.0
    inner_window_count: int = 4
    inner_boundary_cells_per_match_radius: float = 16.0
    max_resident_cells: int = 200_000_000
    max_finest_steps: int = 5_000_000
    max_refinement_levels: int = 14
    max_global_timestep_zone_updates: int = 10_000_000_000_000
    bytes_per_cell: float = 256.0

    def validated(self) -> "PlannerSettings":
        for name, value in asdict(self).items():
            if name in {
                "inner_window_count",
                "max_resident_cells",
                "max_finest_steps",
                "max_refinement_levels",
                "max_global_timestep_zone_updates",
            }:
                if int(value) <= 0:
                    raise ValueError(f"{name} must be positive")
            else:
                _positive(float(value), name)
        if self.overlap_fraction >= 1.0:
            raise ValueError("overlap_fraction must be smaller than one")
        if self.outer_sink_fraction_of_match >= 1.0:
            raise ValueError("outer_sink_fraction_of_match must be smaller than one")
        if self.environment_validity_fraction >= 1.0:
            raise ValueError("environment_validity_fraction must be smaller than one")
        return self


def characteristic_state(
    rho: float,
    pgas: float,
    gamma: float,
    spatial_four_velocity: Iterable[float],
    magnetic_field: Iterable[float],
) -> dict[str, float | list[float] | str | None]:
    """Return source-orthonormal-frame characteristic scale proxies.

    ``u`` is the Eulerian spatial four-velocity W*v and ``B`` is the Eulerian
    magnetic field, matching ``wind_frame=source_tetrad`` in the local problem.
    """

    rho = _positive(rho, "rho")
    pgas = _positive(pgas, "pgas")
    gamma = _positive(gamma, "gamma")
    if gamma <= 1.0:
        raise ValueError("gamma must be greater than one")
    u = _finite_vector(spatial_four_velocity, "spatial four-velocity")
    magnetic = _finite_vector(magnetic_field, "magnetic field")
    usq = _dot(u, u)
    lorentz_factor = math.sqrt(1.0 + usq)
    speed_squared = usq / (1.0 + usq)
    speed = math.sqrt(speed_squared)
    enthalpy = 1.0 + gamma * pgas / ((gamma - 1.0) * rho)
    sound_speed_squared = gamma * pgas / (rho * enthalpy)

    bsq_lab = _dot(magnetic, magnetic)
    bdotu = _dot(magnetic, u)
    comoving_b_squared = (bsq_lab + bdotu * bdotu) / (lorentz_factor**2)
    alfven_speed_squared = comoving_b_squared / (
        rho * enthalpy + comoving_b_squared
    )
    # The perpendicular relativistic fast-mode value is a useful conservative upper
    # proxy.  The exact directional GRMHD characteristics require the full dispersion
    # relation and are intentionally not claimed here.
    fast_proxy_squared = (
        sound_speed_squared
        + alfven_speed_squared
        - sound_speed_squared * alfven_speed_squared
    )
    effective_speed_squared = speed_squared + fast_proxy_squared
    crossing_speed = min(1.0, math.sqrt(effective_speed_squared))
    magnetization = comoving_b_squared / (rho * enthalpy)
    plasma_beta = (
        None if comoving_b_squared == 0.0 else 2.0 * pgas / comoving_b_squared
    )
    fast_mach_proxy = (
        math.inf if fast_proxy_squared == 0.0 else speed / math.sqrt(fast_proxy_squared)
    )
    return {
        "rho": rho,
        "pgas": pgas,
        "adiabatic_index": gamma,
        "spatial_four_velocity": list(u),
        "magnetic_field": list(magnetic),
        "lorentz_factor": lorentz_factor,
        "three_speed": speed,
        "specific_enthalpy": enthalpy,
        "sound_speed_squared": sound_speed_squared,
        "comoving_magnetic_b_squared": comoving_b_squared,
        "alfven_speed_squared": alfven_speed_squared,
        "perpendicular_fast_speed_proxy_squared": fast_proxy_squared,
        "fast_mach_proxy": fast_mach_proxy,
        "magnetization_sigma": magnetization,
        "plasma_beta": plasma_beta,
        "capture_effective_speed_squared": effective_speed_squared,
        "causal_crossing_speed_proxy": crossing_speed,
        "characteristic_warning": (
            "cost proxy only: the exact directional GRMHD fast characteristics are "
            "not evaluated"
        ),
    }


def _domain_cells(
    settings: PlannerSettings, cells_per_capture_radius: float
) -> tuple[int, tuple[int, int, int]]:
    dimensions = (
        _ceil_ratio(
            (settings.upstream_capture_radii + settings.downstream_capture_radii)
            * cells_per_capture_radius
        ),
        _ceil_ratio(2.0 * settings.transverse_capture_radii * cells_per_capture_radius),
        _ceil_ratio(2.0 * settings.transverse_capture_radii * cells_per_capture_radius),
    )
    return math.prod(dimensions), dimensions


def _cost_model(
    base_cells: int,
    base_dimensions: tuple[int, int, int],
    levels: int,
    duration: float,
    finest_dx: float,
    settings: PlannerSettings,
) -> dict[str, float | int | list[int]]:
    finest_steps = _ceil_ratio(duration / (settings.cfl * finest_dx))
    resident_cells = base_cells * (levels + 1)
    subcycle_factor = sum(0.5**offset for offset in range(levels + 1))
    return {
        "base_grid_dimensions": list(base_dimensions),
        "base_grid_cells": base_cells,
        "nested_refinement_levels": levels,
        "optimistic_nested_resident_cells": resident_cells,
        "optimistic_nested_memory_gib": (
            resident_cells * settings.bytes_per_cell / 1024.0**3
        ),
        "finest_dx": finest_dx,
        "finest_steps": finest_steps,
        "global_timestep_zone_updates": resident_cells * finest_steps,
        "ideal_level_subcycled_zone_updates": int(
            math.ceil(base_cells * finest_steps * subcycle_factor)
        ),
        "cost_warning": (
            "nested boxes are an optimistic lower-bound geometry; load imbalance, "
            "guard zones, solver iterations, and output are excluded"
        ),
    }


def _finite_scales(scales: dict[str, float]) -> dict[str, float]:
    result: dict[str, float] = {}
    for name, value in scales.items():
        if not name or any(character.isspace() for character in name):
            raise ValueError("coherence-scale names must be nonempty and contain no spaces")
        result[name] = _positive(value, f"coherence scale {name}")
    return result


def build_plan(
    *,
    secondary_mass: float,
    primary_mass: float,
    orbital_radius: float,
    rho: float,
    pgas: float,
    gamma: float,
    spatial_four_velocity: Iterable[float],
    magnetic_field: Iterable[float],
    coherence_scales: dict[str, float] | None = None,
    settings: PlannerSettings | None = None,
) -> dict[str, object]:
    """Build the hierarchy plan in one common set of geometrized units."""

    secondary_mass = _positive(secondary_mass, "secondary_mass")
    primary_mass = _positive(primary_mass, "primary_mass")
    orbital_radius = _positive(orbital_radius, "orbital_radius")
    settings = (settings or PlannerSettings()).validated()
    state = characteristic_state(
        rho, pgas, gamma, spatial_four_velocity, magnetic_field
    )
    effective_speed_squared = float(state["capture_effective_speed_squared"])
    if effective_speed_squared <= 0.0:
        raise ValueError("the capture-radius proxy is singular for a static cold state")

    capture_radius_one = secondary_mass / effective_speed_squared
    capture_radius = 2.0 * capture_radius_one
    mass_ratio = secondary_mass / primary_mass
    hill_radius = orbital_radius * (mass_ratio / 3.0) ** (1.0 / 3.0)
    user_scales = _finite_scales(coherence_scales or {})
    limiting_scales = {"hill_radius": hill_radius, **user_scales}
    limiting_name, limiting_scale = min(limiting_scales.items(), key=lambda item: item[1])
    capture_to_limit = capture_radius / limiting_scale

    crossing_speed = float(state["causal_crossing_speed_proxy"])
    capture_crossing_time = capture_radius / crossing_speed
    direct_duration = settings.flow_crossings * capture_crossing_time
    requested_inner_dx = secondary_mass / settings.cells_per_secondary_mass
    base_dx = capture_radius / settings.outer_cells_per_capture_radius
    direct_levels = _refinement_levels(base_dx, requested_inner_dx)
    direct_finest_dx = base_dx / 2.0**direct_levels
    base_cells, base_dimensions = _domain_cells(
        settings, settings.outer_cells_per_capture_radius
    )
    direct_cost = _cost_model(
        base_cells,
        base_dimensions,
        direct_levels,
        direct_duration,
        direct_finest_dx,
        settings,
    )
    uniform_dimensions = (
        _ceil_ratio(
            (settings.upstream_capture_radii + settings.downstream_capture_radii)
            * capture_radius
            / requested_inner_dx
        ),
        _ceil_ratio(
            2.0
            * settings.transverse_capture_radii
            * capture_radius
            / requested_inner_dx
        ),
        _ceil_ratio(
            2.0
            * settings.transverse_capture_radii
            * capture_radius
            / requested_inner_dx
        ),
    )
    direct_cost["uniform_horizon_resolving_dimensions"] = list(uniform_dimensions)
    direct_cost["uniform_horizon_resolving_cells"] = math.prod(uniform_dimensions)

    asymptotic_inner_radius = settings.inner_asymptotic_radius_in_m * secondary_mass
    match_radius = math.sqrt(asymptotic_inner_radius * capture_radius)
    inner_to_match = asymptotic_inner_radius / match_radius
    match_to_capture = match_radius / capture_radius
    clean_overlap = (
        inner_to_match <= settings.overlap_fraction
        and match_to_capture <= settings.overlap_fraction
    )
    outer_sink_radius = settings.outer_sink_fraction_of_match * match_radius
    outer_target_dx = outer_sink_radius / settings.cells_per_outer_sink_radius
    outer_levels = _refinement_levels(base_dx, outer_target_dx)
    outer_finest_dx = base_dx / 2.0**outer_levels
    outer_cost = _cost_model(
        base_cells,
        base_dimensions,
        outer_levels,
        direct_duration,
        outer_finest_dx,
        settings,
    )
    inner_window_duration = (
        settings.inner_window_crossings * match_radius / crossing_speed
    )
    inner_base_dimension = _ceil_ratio(
        2.0 * settings.inner_boundary_cells_per_match_radius
    )
    inner_base_dimensions = (
        inner_base_dimension,
        inner_base_dimension,
        inner_base_dimension,
    )
    inner_base_cells = math.prod(inner_base_dimensions)
    inner_coarse_dx = (
        match_radius / settings.inner_boundary_cells_per_match_radius
    )
    inner_levels = _refinement_levels(inner_coarse_dx, requested_inner_dx)
    inner_finest_dx = inner_coarse_dx / 2.0**inner_levels
    inner_window_cost = _cost_model(
        base_cells=inner_base_cells,
        base_dimensions=inner_base_dimensions,
        levels=inner_levels,
        duration=inner_window_duration,
        finest_dx=inner_finest_dx,
        settings=settings,
    )
    inner_steps_per_window = int(inner_window_cost["finest_steps"])
    continuous_inner_steps = _ceil_ratio(
        direct_duration / (settings.cfl * inner_finest_dx)
    )
    continuous_inner_zone_updates = (
        int(inner_window_cost["optimistic_nested_resident_cells"])
        * continuous_inner_steps
    )
    ensemble_inner_zone_updates = (
        settings.inner_window_count
        * int(inner_window_cost["global_timestep_zone_updates"])
    )

    direct_affordable = (
        int(direct_cost["optimistic_nested_resident_cells"])
        <= settings.max_resident_cells
        and int(direct_cost["finest_steps"]) <= settings.max_finest_steps
        and direct_levels <= settings.max_refinement_levels
        and int(direct_cost["global_timestep_zone_updates"])
        <= settings.max_global_timestep_zone_updates
    )
    matched_affordable = (
        int(outer_cost["optimistic_nested_resident_cells"])
        <= settings.max_resident_cells
        and int(outer_cost["finest_steps"]) <= settings.max_finest_steps
        and outer_levels <= settings.max_refinement_levels
        and inner_steps_per_window <= settings.max_finest_steps
        and int(outer_cost["global_timestep_zone_updates"])
        <= settings.max_global_timestep_zone_updates
        and int(inner_window_cost["optimistic_nested_resident_cells"])
        <= settings.max_resident_cells
        and inner_levels <= settings.max_refinement_levels
        and ensemble_inner_zone_updates
        <= settings.max_global_timestep_zone_updates
    )
    uniform_environment_valid = (
        capture_to_limit <= settings.environment_validity_fraction
    )
    nearly_spherical_reduced_option = (
        float(state["fast_mach_proxy"]) < 0.3
        and float(state["magnetization_sigma"]) < 1.0e-2
        and uniform_environment_valid
    )

    if not uniform_environment_valid:
        recommendation = "global_or_shearing_outer_with_relativistic_inner"
        rationale = (
            "the nominal capture radius is not small compared with the limiting "
            "disk/tidal coherence scale, so an enormous uniform BHL box would model "
            "the wrong upstream environment"
        )
    elif direct_affordable:
        recommendation = "direct_grmhd"
        rationale = "the configured direct nested-grid cell and finest-step budgets pass"
    elif clean_overlap and matched_affordable:
        recommendation = "matched_bhl_outer_inner"
        rationale = (
            "direct horizon-to-capture evolution exceeds the configured budget and a "
            "two-sided asymptotic overlap exists"
        )
    elif nearly_spherical_reduced_option:
        recommendation = "bondi_michel_reduced_outer_or_direct_amr"
        rationale = (
            "the flow is sub-fast and weakly magnetized but no affordable clean BHL "
            "matching hierarchy was found"
        )
    else:
        recommendation = "direct_amr_or_reduced_outer_without_clean_overlap"
        rationale = (
            "neither the direct budget nor the clean matched-asymptotic criteria pass; "
            "a sink/worldtube split would be uncontrolled without a convergence study"
        )

    overlap_low = (1.0 - settings.overlap_fraction) * match_radius
    overlap_high = (1.0 + settings.overlap_fraction) * match_radius
    transfer_contract = {
        "status": (
            "fixed-grid outer writer and inner CT magnetic replay implemented; "
            "moving/source-frame pullback and CT cochain projection implemented; "
            "cut-surface/AMR sampling pending; operational HLLE boundary replay and "
            "host seven-wave reference projector plus packed-data characteristic audit "
            "implemented; fixed-grid outer-to-inner closure audit implemented; "
            "mode-resolved reflection test and device seven-wave projection pending"
        ),
        "geometry": "secondary-centered cubical worldtube in source-tetrad axes",
        "cube_half_width": match_radius,
        "radial_buffer_interval": [overlap_low, overlap_high],
        "outer_sink_radius": outer_sink_radius,
        "required_outer_fields": [
            "time-centered fluid primitive or conserved states",
            "face-integrated magnetic fluxes on all six worldtube faces",
            "edge-integrated electric fields/EMFs at matching time levels",
            "metric/tetrad metadata and coordinate Jacobians",
        ],
        "regridding_rules": [
            "interpolate incoming fluid characteristic amplitudes; retain outgoing modes",
            "conservatively restrict/prolong normal face flux so child sums equal parent flux",
            "use one shared edge EMF to update every face incident on that edge",
            "interpolate in source proper time and preserve the outer integrator time centering",
        ],
        "force_partition": {
            "formula": (
                "F_total = -Fmom_inner + integral(w_inner*f_inner_rel) "
                "+ integral(w_outer*f_outer_grav), with w_inner+w_outer=1"
            ),
            "rule": (
                "discard the outer sink momentum force when the relativistic inner solution "
                "replaces it; count accreted momentum exactly once"
            ),
        },
        "mandatory_convergence_variations": [
            "move the worldtube/matching radius within the overlap",
            (
                "halve the outer sink radius or show the sink is causally downstream "
                "of the fast surface"
            ),
            "double worldtube spatial resolution",
            "halve output cadence and compare flux/EMF injection",
            "vary the partition-of-unity overlap width",
        ],
    }

    warnings = [
        str(state["characteristic_warning"]),
        (
            "the factor-two capture radius is used for conservative cost estimates; "
            "the factor-one convention is reported separately"
        ),
        (
            "continuous replay for the full outer settling interval retains the inner "
            "finest-step cost; steady means or separated stationary windows reduce cost "
            "but cannot recover low-frequency cross-correlations"
        ),
    ]
    if not uniform_environment_valid:
        warnings.append(
            f"capture_radius/{limiting_name}={capture_to_limit:.6g} exceeds the "
            f"configured {settings.environment_validity_fraction:.6g} validity fraction"
        )
    if not clean_overlap:
        warnings.append(
            "the proposed geometric-mean matching radius has no two-sided scale separation"
        )
    if float(state["fast_mach_proxy"]) <= 1.0:
        warnings.append(
            "the flow is not fast-supersonic by the proxy; upstream/downstream boundary "
            "causality and a Bondi/Michel-like alternative need explicit checks"
        )

    return {
        "classification": CLASSIFICATION,
        "recommendation": recommendation,
        "rationale": rationale,
        "units": "one common G=c=1 unit system for every supplied mass and length",
        "input": {
            "secondary_mass": secondary_mass,
            "primary_mass": primary_mass,
            "orbital_radius": orbital_radius,
            "mass_ratio": mass_ratio,
            "coherence_scales": user_scales,
            "settings": asdict(settings),
        },
        "upstream_state": state,
        "scales": {
            "capture_radius_factor_one": capture_radius_one,
            "capture_radius_factor_two_for_cost": capture_radius,
            "capture_radius_in_secondary_masses": capture_radius / secondary_mass,
            "hill_radius": hill_radius,
            "limiting_environment_scale_name": limiting_name,
            "limiting_environment_scale": limiting_scale,
            "capture_to_limiting_scale": capture_to_limit,
            "capture_crossing_time": capture_crossing_time,
            "direct_settling_duration": direct_duration,
            "inner_asymptotic_radius": asymptotic_inner_radius,
            "matching_radius": match_radius,
            "matching_radius_in_secondary_masses": match_radius / secondary_mass,
            "inner_asymptotic_to_match": inner_to_match,
            "match_to_capture": match_to_capture,
            "clean_overlap": clean_overlap,
        },
        "validity": {
            "uniform_bhl_environment": uniform_environment_valid,
            "direct_within_budget": direct_affordable,
            "matched_components_within_budget": matched_affordable,
            "weakly_magnetized_subfast_reduced_option": nearly_spherical_reduced_option,
        },
        "costs": {
            "direct": direct_cost,
            "matched_outer": outer_cost,
            "matched_inner": {
                **inner_window_cost,
                "duration_per_stationary_window": inner_window_duration,
                "finest_steps_per_stationary_window": inner_steps_per_window,
                "stationary_window_count": settings.inner_window_count,
                "ensemble_finest_steps": (
                    settings.inner_window_count * inner_steps_per_window
                ),
                "ensemble_global_timestep_zone_updates": ensemble_inner_zone_updates,
                "continuous_full_outer_duration_finest_steps": continuous_inner_steps,
                "continuous_full_outer_duration_zone_updates": (
                    continuous_inner_zone_updates
                ),
            },
        },
        "matching_contract": transfer_contract,
        "warnings": warnings,
    }


def _format_number(value: object) -> str:
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "infinity"
        return f"{value:.8g}"
    return str(value)


def markdown_report(plan: dict[str, object]) -> str:
    scales = plan["scales"]
    validity = plan["validity"]
    costs = plan["costs"]
    assert isinstance(scales, dict)
    assert isinstance(validity, dict)
    assert isinstance(costs, dict)
    direct = costs["direct"]
    outer = costs["matched_outer"]
    inner = costs["matched_inner"]
    assert isinstance(direct, dict) and isinstance(outer, dict) and isinstance(inner, dict)

    scale_rows = [
        ("factor-one capture radius", scales["capture_radius_factor_one"]),
        ("factor-two cost radius", scales["capture_radius_factor_two_for_cost"]),
        ("capture radius / m", scales["capture_radius_in_secondary_masses"]),
        ("Hill radius", scales["hill_radius"]),
        ("limiting environment scale", scales["limiting_environment_scale"]),
        ("capture / limiting scale", scales["capture_to_limiting_scale"]),
        ("matching radius / m", scales["matching_radius_in_secondary_masses"]),
        ("clean two-sided overlap", scales["clean_overlap"]),
    ]
    cost_rows = [
        (
            "direct",
            direct["optimistic_nested_resident_cells"],
            direct["nested_refinement_levels"],
            direct["finest_steps"],
        ),
        (
            "matched outer",
            outer["optimistic_nested_resident_cells"],
            outer["nested_refinement_levels"],
            outer["finest_steps"],
        ),
        (
            "inner one stationary window",
            inner["optimistic_nested_resident_cells"],
            inner["nested_refinement_levels"],
            inner["finest_steps_per_stationary_window"],
        ),
    ]
    lines = [
        "# EMRI BHL hierarchy plan",
        "",
        f"Recommendation: `{plan['recommendation']}`.",
        "",
        str(plan["rationale"]) + ".",
        "",
        "## Characteristic scales",
        "",
        "| quantity | estimate |",
        "|---|---:|",
    ]
    lines.extend(f"| {name} | {_format_number(value)} |" for name, value in scale_rows)
    lines.extend(
        [
            "",
            "## Optimistic cost floor",
            "",
            "| model | resident cells | refinement levels | finest steps |",
            "|---|---:|---:|---:|",
        ]
    )
    lines.extend(
        "| " + " | ".join(_format_number(value) for value in row) + " |"
        for row in cost_rows
    )
    lines.extend(
        [
            "",
            "The nested-box resident count is an optimistic geometry estimate. The direct",
            "uniform horizon-resolving cell count would be",
            f"`{_format_number(direct['uniform_horizon_resolving_cells'])}`.",
            "",
            "## Gates",
            "",
        ]
    )
    for name, value in validity.items():
        lines.append(f"- `{name}`: {_format_number(value)}")
    lines.extend(["", "## Required outer--inner contract", ""])
    contract = plan["matching_contract"]
    assert isinstance(contract, dict)
    lines.append(f"Status: {contract['status']}.")
    lines.extend(["", "Required outer data:", ""])
    for field in contract["required_outer_fields"]:
        lines.append(f"- {field}")
    lines.extend(["", "Transfer rules:", ""])
    for rule in contract["regridding_rules"]:
        lines.append(f"- {rule}")
    force_partition = contract["force_partition"]
    assert isinstance(force_partition, dict)
    lines.extend(
        [
            "",
            f"Force bookkeeping: `{force_partition['formula']}`.",
            "",
            str(force_partition["rule"]) + ".",
            "",
            "## Warnings",
            "",
        ]
    )
    for warning in plan["warnings"]:
        lines.append(f"- {warning}")
    lines.append("")
    return "\n".join(lines)


def _named_scale(specification: str) -> tuple[str, float]:
    if "=" not in specification:
        raise argparse.ArgumentTypeError("coherence scale must be NAME=VALUE")
    name, value_text = specification.split("=", 1)
    try:
        value = float(value_text)
        _positive(value, f"coherence scale {name}")
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error
    if not name or any(character.isspace() for character in name):
        raise argparse.ArgumentTypeError(
            "coherence-scale names must be nonempty and contain no spaces"
        )
    return name, value


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--secondary-mass", type=float, default=1.0)
    parser.add_argument("--primary-mass", type=float, required=True)
    parser.add_argument("--orbital-radius", type=float, required=True)
    parser.add_argument("--rho", type=float, required=True)
    parser.add_argument("--pgas", type=float, required=True)
    parser.add_argument("--adiabatic-index", type=float, default=4.0 / 3.0)
    parser.add_argument("--four-velocity", type=float, nargs=3, required=True,
                        metavar=("U1", "U2", "U3"))
    parser.add_argument("--magnetic-field", type=float, nargs=3, required=True,
                        metavar=("B1", "B2", "B3"))
    parser.add_argument(
        "--coherence-scale",
        action="append",
        type=_named_scale,
        default=[],
        metavar="NAME=VALUE",
        help="repeat for disk H and density/pressure/velocity/magnetic gradient scales",
    )
    parser.add_argument("--cells-per-secondary-mass", type=float, default=16.0)
    parser.add_argument("--outer-cells-per-capture-radius", type=float, default=16.0)
    parser.add_argument("--cells-per-outer-sink-radius", type=float, default=8.0)
    parser.add_argument("--flow-crossings", type=float, default=5.0)
    parser.add_argument("--inner-boundary-cells-per-match-radius", type=float,
                        default=16.0)
    parser.add_argument("--max-resident-cells", type=int, default=200_000_000)
    parser.add_argument("--max-finest-steps", type=int, default=5_000_000)
    parser.add_argument("--max-refinement-levels", type=int, default=14)
    parser.add_argument("--max-global-timestep-zone-updates", type=int,
                        default=10_000_000_000_000)
    parser.add_argument("--output-prefix", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    settings = PlannerSettings(
        cells_per_secondary_mass=arguments.cells_per_secondary_mass,
        outer_cells_per_capture_radius=arguments.outer_cells_per_capture_radius,
        cells_per_outer_sink_radius=arguments.cells_per_outer_sink_radius,
        flow_crossings=arguments.flow_crossings,
        inner_boundary_cells_per_match_radius=(
            arguments.inner_boundary_cells_per_match_radius
        ),
        max_resident_cells=arguments.max_resident_cells,
        max_finest_steps=arguments.max_finest_steps,
        max_refinement_levels=arguments.max_refinement_levels,
        max_global_timestep_zone_updates=(
            arguments.max_global_timestep_zone_updates
        ),
    )
    try:
        plan = build_plan(
            secondary_mass=arguments.secondary_mass,
            primary_mass=arguments.primary_mass,
            orbital_radius=arguments.orbital_radius,
            rho=arguments.rho,
            pgas=arguments.pgas,
            gamma=arguments.adiabatic_index,
            spatial_four_velocity=arguments.four_velocity,
            magnetic_field=arguments.magnetic_field,
            coherence_scales=dict(arguments.coherence_scale),
            settings=settings,
        )
    except ValueError as error:
        raise SystemExit(f"error: {error}") from error

    prefix = arguments.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    json_path = prefix.with_suffix(".json")
    markdown_path = prefix.with_suffix(".md")
    json_path.write_text(
        json.dumps(plan, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    markdown_path.write_text(markdown_report(plan), encoding="utf-8")
    print(f"wrote {json_path}")
    print(f"wrote {markdown_path}")
    print(f"recommendation: {plan['recommendation']}")


if __name__ == "__main__":
    main()
