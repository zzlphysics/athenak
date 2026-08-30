#!/usr/bin/env python3
"""Plan a resolved and affordable numerical-ADM EMRI inner replay.

This planner is intentionally fail closed.  Increasing the number of worldtube
face cells does not create source resolution, and satisfying the horizon gate
does not by itself make a direct BHL domain large enough.  The report therefore
keeps four independent questions separate:

* can one secondary mass satisfy both horizon resolution and boundary distance;
* does the cube contain the requested number of BHL capture radii;
* does the global GRMHD/ADM source actually resolve the requested local cells;
* do the fluid evolution and ADM-volume construction fit the declared budgets.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np

import extract_global_worldtube as extract
import plan_bhl_hierarchy as hierarchy
import run_adm_inner_replay_pilot as pilot
import worldtube_flux_emf as worldtube


CLASSIFICATION = "athenak-emri-adm-inner-replay-plan-v1"
REQUIRED_SOURCE_PROVENANCE = {
    "metric_content": "primary_only",
    "fluid_content": "global_grmhd",
    "secondary_backreaction": "absent",
}


def _finite_positive(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def _condition(
    relation: str, observed: float, threshold: float
) -> dict[str, float | str | bool]:
    if relation not in ("minimum", "maximum"):
        raise ValueError("condition relation must be minimum or maximum")
    return {
        "relation": relation,
        "observed": observed,
        "threshold": threshold,
        "passed": observed >= threshold if relation == "minimum" else observed <= threshold,
    }


def horizon_radius_per_mass(chi: float) -> float:
    spin = float(chi)
    if not math.isfinite(spin) or abs(spin) > 1.0:
        raise ValueError("secondary chi must lie in [-1,1]")
    return 1.0 + math.sqrt(max(0.0, 1.0 - spin * spin))


def mass_window(
    *,
    half_width: float,
    cells: int,
    secondary_chi: float,
    minimum_horizon_cells: float,
    minimum_boundary_horizon_radii: float,
    capture_radius_per_mass: float | None = None,
    minimum_boundary_capture_radii: float | None = None,
) -> dict[str, object]:
    """Return the secondary-mass interval allowed by geometric constraints."""

    half_width = _finite_positive(half_width, "half width")
    if not isinstance(cells, int) or isinstance(cells, bool) or cells < 1:
        raise ValueError("cells must be a positive integer")
    minimum_horizon_cells = _finite_positive(
        minimum_horizon_cells, "minimum horizon cells"
    )
    minimum_boundary_horizon_radii = _finite_positive(
        minimum_boundary_horizon_radii,
        "minimum boundary horizon radii",
    )
    radius_per_mass = horizon_radius_per_mass(secondary_chi)
    spacing = 2.0 * half_width / cells
    lower = minimum_horizon_cells * spacing / radius_per_mass
    horizon_upper = (
        half_width / (minimum_boundary_horizon_radii * radius_per_mass)
    )
    upper_constraints = {"worldtube_horizon_separation": horizon_upper}
    if (capture_radius_per_mass is None) != (
        minimum_boundary_capture_radii is None
    ):
        raise ValueError(
            "capture radius and minimum boundary capture radii must be supplied together"
        )
    if capture_radius_per_mass is not None:
        capture_radius_per_mass = _finite_positive(
            capture_radius_per_mass, "capture radius per mass"
        )
        minimum_boundary_capture_radii = _finite_positive(
            minimum_boundary_capture_radii,
            "minimum boundary capture radii",
        )
        upper_constraints["worldtube_capture_separation"] = (
            half_width
            / (minimum_boundary_capture_radii * capture_radius_per_mass)
        )
    upper_name, upper = min(upper_constraints.items(), key=lambda item: item[1])
    feasible = bool(lower <= upper * (1.0 + 128.0 * np.finfo(float).eps))
    recommended = math.sqrt(lower * upper) if feasible else None
    return {
        "feasible": feasible,
        "cell_spacing": spacing,
        "horizon_radius_per_secondary_mass": radius_per_mass,
        "minimum_secondary_mass": lower,
        "maximum_secondary_mass": upper,
        "limiting_upper_constraint": upper_name,
        "upper_constraints": upper_constraints,
        "balanced_secondary_mass": recommended,
    }


def minimum_cells_for_geometry(
    minimum_horizon_cells: float, minimum_boundary_horizon_radii: float
) -> int:
    value = 2.0 * _finite_positive(
        minimum_horizon_cells, "minimum horizon cells"
    ) * _finite_positive(
        minimum_boundary_horizon_radii, "minimum boundary horizon radii"
    )
    return int(math.ceil(value - 64.0 * np.finfo(float).eps * max(1.0, value)))


def minimum_cells_for_direct_capture(
    *,
    minimum_horizon_cells: float,
    minimum_boundary_capture_radii: float,
    capture_radius_per_mass: float,
    secondary_chi: float,
) -> int:
    value = (
        2.0
        * _finite_positive(minimum_horizon_cells, "minimum horizon cells")
        * _finite_positive(
            minimum_boundary_capture_radii, "minimum boundary capture radii"
        )
        * _finite_positive(capture_radius_per_mass, "capture radius per mass")
        / horizon_radius_per_mass(secondary_chi)
    )
    return int(math.ceil(value - 64.0 * np.finfo(float).eps * max(1.0, value)))


def mesh_friendly_cells(minimum_cells: int, minimum_meshblock_cells: int = 8) -> int:
    """Round upward until AthenaK can use a reasonably sized <=16-cell block."""

    if not isinstance(minimum_cells, int) or minimum_cells < 1:
        raise ValueError("minimum cells must be a positive integer")
    if not isinstance(minimum_meshblock_cells, int) or not 1 <= minimum_meshblock_cells <= 16:
        raise ValueError("minimum meshblock cells must lie in [1,16]")
    candidate = minimum_cells
    while True:
        divisor = next(
            value
            for value in range(min(candidate, 16), 0, -1)
            if candidate % value == 0
        )
        if divisor >= minimum_meshblock_cells:
            return candidate
        candidate += 1


def source_resolution_audit(
    *,
    source_spacing: Iterable[float],
    spatial_legs: Iterable[object],
    local_cell_spacing: float,
    half_width: float,
    minimum_source_cells_per_local_cell: float,
) -> dict[str, object]:
    """Compare mapped local-cell lengths with the source grid spacing.

    The smallest singular value of ``diag(1/dx_source) e^i_a`` is the
    least-resolved local direction in source-cell coordinates.  Time tilt and
    temporal cadence are deliberately not folded into this spatial audit.
    """

    source = np.asarray(tuple(source_spacing), dtype=np.float64)
    if source.shape != (3,) or not np.isfinite(source).all() or np.min(source) <= 0.0:
        raise ValueError("source spacing must contain three finite positive values")
    local_cell_spacing = _finite_positive(local_cell_spacing, "local cell spacing")
    half_width = _finite_positive(half_width, "half width")
    minimum_source_cells_per_local_cell = _finite_positive(
        minimum_source_cells_per_local_cell,
        "minimum source cells per local cell",
    )
    singular_values = []
    raw_singular_values = []
    for value in spatial_legs:
        legs = np.asarray(value, dtype=np.float64)
        if legs.shape != (4, 3) or not np.isfinite(legs).all():
            raise ValueError("each spatial-leg sample must have finite shape (4,3)")
        spatial = legs[1:, :]
        raw = np.linalg.svd(spatial, compute_uv=False)
        normalized = np.linalg.svd(spatial / source[:, None], compute_uv=False)
        if raw[-1] <= 0.0 or normalized[-1] <= 0.0:
            raise ValueError("frame spatial legs are singular")
        raw_singular_values.append(float(raw[-1]))
        singular_values.append(float(normalized[-1]))
    if not singular_values:
        raise ValueError("source resolution audit requires at least one frame")
    least_resolved = min(singular_values)
    least_raw = min(raw_singular_values)
    cells_per_local_cell = least_resolved * local_cell_spacing
    cells_across_diameter = least_resolved * (2.0 * half_width)
    maximum_isotropic_source_spacing = (
        least_raw * local_cell_spacing / minimum_source_cells_per_local_cell
    )
    return {
        "passed": bool(
            cells_per_local_cell >= minimum_source_cells_per_local_cell
        ),
        "minimum_source_cells_per_local_cell": cells_per_local_cell,
        "required_source_cells_per_local_cell": minimum_source_cells_per_local_cell,
        "minimum_source_cells_across_worldtube_diameter": cells_across_diameter,
        "maximum_isotropic_global_source_spacing_for_target": (
            maximum_isotropic_source_spacing
        ),
        "minimum_normalized_spatial_jacobian_singular_value": least_resolved,
        "minimum_raw_spatial_jacobian_singular_value": least_raw,
        "source_spacing_global_coordinates": source.tolist(),
        "scope": (
            "conservative spatial oversampling audit; temporal tilt and source cadence "
            "are assessed separately"
        ),
    }


def resource_estimate(
    *,
    fluid_cells: int,
    metric_cells: int,
    mesh_nghost: int,
    metric_halo: int,
    metric_samples: int,
    duration: float,
    cfl: float,
    half_width: float,
    fluid_bytes_per_allocated_cell: float,
) -> dict[str, object]:
    for name, value in (
        ("fluid_cells", fluid_cells),
        ("metric_cells", metric_cells),
        ("mesh_nghost", mesh_nghost),
        ("metric_halo", metric_halo),
        ("metric_samples", metric_samples),
    ):
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    duration = _finite_positive(duration, "duration")
    cfl = _finite_positive(cfl, "CFL number")
    half_width = _finite_positive(half_width, "half width")
    fluid_bytes_per_allocated_cell = _finite_positive(
        fluid_bytes_per_allocated_cell, "fluid bytes per allocated cell"
    )
    meshblock = next(
        value
        for value in range(min(fluid_cells, 16), 0, -1)
        if fluid_cells % value == 0
    )
    blocks_per_axis = fluid_cells // meshblock
    meshblocks = blocks_per_axis**3
    allocated_fluid_cells = meshblocks * (meshblock + 2 * mesh_nghost) ** 3
    active_fluid_cells = fluid_cells**3
    local_spacing = 2.0 * half_width / fluid_cells
    causal_step_proxy = int(math.ceil(duration / (cfl * local_spacing)))
    metric_nodes = metric_cells + 2 * metric_halo
    metric_node_count = metric_nodes**3
    gib = 1024.0**3
    binary_bytes = metric_samples * 16 * metric_node_count * 8
    device_two_slab_bytes = 2 * 16 * metric_node_count * 8
    # build_volume simultaneously retains g_mu_nu, K_ij, and the packed 16-field
    # output.  This is a strict array floor and excludes sampler/cache temporaries.
    builder_floor_bytes = metric_samples * (16 + 9 + 16) * metric_node_count * 8
    return {
        "fluid": {
            "active_cells": active_fluid_cells,
            "meshblock_cells_per_axis": meshblock,
            "meshblock_count": meshblocks,
            "ghost_inclusive_allocated_cells": allocated_fluid_cells,
            "resident_memory_proxy_gib": (
                allocated_fluid_cells * fluid_bytes_per_allocated_cell / gib
            ),
            "causal_cfl_step_proxy": causal_step_proxy,
            "rk3_active_zone_updates_proxy": 3 * causal_step_proxy * active_fluid_cells,
            "note": (
                "CFL estimate assumes unit coordinate signal speed; shift/lapse and "
                "solver overhead are excluded"
            ),
        },
        "adm_volume": {
            "active_cells_per_axis": metric_cells,
            "halo_nodes_per_side": metric_halo,
            "nodes_per_axis": metric_nodes,
            "sample_count": metric_samples,
            "binary_gib": binary_bytes / gib,
            "device_two_time_slab_gib": device_two_slab_bytes / gib,
            "python_builder_array_floor_gib": builder_floor_bytes / gib,
            "note": (
                "builder value is a lower bound for current uniform-volume generation; "
                "snapshot caches and temporary derivatives are excluded"
            ),
        },
    }


def convergence_cost_matrix(
    *,
    fluid_cells: int,
    resolution_factors: Iterable[float],
    mesh_nghost: int,
    metric_samples: int,
    duration: float,
    cfl: float,
    half_width: float,
    fluid_bytes_per_allocated_cell: float,
) -> dict[str, object]:
    """Estimate uniform ADM-volume cost without allocating any volume arrays."""

    result: dict[str, object] = {}
    for raw_factor in resolution_factors:
        factor = _finite_positive(raw_factor, "metric resolution factor")
        metric_cells = max(2, int(math.ceil(fluid_cells * factor)))
        halo = pilot.minimum_metric_halo(metric_cells, fluid_cells, mesh_nghost)
        estimate = resource_estimate(
            fluid_cells=fluid_cells,
            metric_cells=metric_cells,
            mesh_nghost=mesh_nghost,
            metric_halo=halo,
            metric_samples=metric_samples,
            duration=duration,
            cfl=cfl,
            half_width=half_width,
            fluid_bytes_per_allocated_cell=fluid_bytes_per_allocated_cell,
        )
        result[f"{factor:.8g}"] = {
            "metric_resolution_factor": factor,
            "metric_cells_per_axis": metric_cells,
            "minimum_metric_halo": halo,
            **estimate["adm_volume"],
        }
    return result


def _read_case(
    campaign_path: Path, requested_case: str | None
) -> tuple[Path, dict[str, object], str, dict[str, object]]:
    path = campaign_path.expanduser().resolve(strict=True)
    campaign = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(campaign, dict) or not isinstance(campaign.get("cases"), dict):
        raise ValueError("campaign summary has no cases object")
    name = requested_case or campaign.get("reference_case")
    if not isinstance(name, str) or name not in campaign["cases"]:
        raise ValueError("requested/default campaign case is unavailable")
    case = campaign["cases"][name]
    if not isinstance(case, dict):
        raise ValueError("campaign case is not an object")
    return path, campaign, name, case


def _provenance_assessment(manifest: dict[str, object]) -> dict[str, object]:
    supplied = manifest.get("source_provenance")
    conditions = {
        name: isinstance(supplied, dict) and supplied.get(name) == expected
        for name, expected in REQUIRED_SOURCE_PROVENANCE.items()
    }
    return {
        "passed": all(conditions.values()),
        "required": REQUIRED_SOURCE_PROVENANCE,
        "supplied": supplied,
        "conditions": conditions,
        "note": (
            "the numerical ADM volume adds the analytic secondary, so its global "
            "metric source must contain the primary only"
        ),
    }


def build_plan(arguments: argparse.Namespace) -> dict[str, object]:
    campaign_path, campaign, case_name, case = _read_case(
        arguments.campaign, arguments.case
    )
    worldtube_path = Path(case["worldtube"]).expanduser().resolve(strict=True)
    times, faces, metadata = worldtube.read_worldtube(worldtube_path)
    face = faces[worldtube.FACE_NAMES[0]]
    cells = int(face.cell_state.shape[-1])
    half_width = _finite_positive(float(metadata["half_width"]), "half width")
    local_spacing = 2.0 * half_width / cells

    boundary = case.get("boundary_state")
    if not isinstance(boundary, dict):
        raise ValueError("campaign case predates the boundary-state audit")
    means = boundary.get("initial_face_mean_by_variable")
    if not isinstance(means, dict):
        raise ValueError("campaign case has no initial boundary means")
    override = arguments.upstream_rho is not None
    rho = float(arguments.upstream_rho if override else means["rho"])
    pgas = float(arguments.upstream_pgas if override else means["pgas"])
    four_velocity = (
        tuple(float(value) for value in arguments.upstream_four_velocity)
        if override
        else tuple(float(means[f"u{axis}"]) for axis in (1, 2, 3))
    )
    magnetic_field = (
        tuple(float(value) for value in arguments.upstream_magnetic_field)
        if override
        else tuple(float(means[f"bcc{axis}"]) for axis in (1, 2, 3))
    )
    state = hierarchy.characteristic_state(
        rho,
        pgas,
        arguments.adiabatic_index,
        four_velocity,
        magnetic_field,
    )
    effective_speed_squared = float(state["capture_effective_speed_squared"])
    capture_radius_per_mass = 2.0 / effective_speed_squared

    geometry_window = mass_window(
        half_width=half_width,
        cells=cells,
        secondary_chi=arguments.secondary_chi,
        minimum_horizon_cells=arguments.minimum_horizon_cells,
        minimum_boundary_horizon_radii=arguments.minimum_boundary_horizon_radii,
    )
    direct_window = mass_window(
        half_width=half_width,
        cells=cells,
        secondary_chi=arguments.secondary_chi,
        minimum_horizon_cells=arguments.minimum_horizon_cells,
        minimum_boundary_horizon_radii=arguments.minimum_boundary_horizon_radii,
        capture_radius_per_mass=capture_radius_per_mass,
        minimum_boundary_capture_radii=arguments.minimum_boundary_capture_radii,
    )

    manifest_path = Path(case["manifest"]).expanduser().resolve(strict=True)
    scan = extract.scan_snapshot_manifest(manifest_path)
    manifest = scan.document
    frame_document = extract._frame_document(manifest.get("frame"), scan.path.parent)
    frames = extract.AffineFrameSeries.from_document(frame_document)
    legs = [
        frames.evaluate(float(time))[1].spatial_legs for time in frames.times
    ]
    source_audit = source_resolution_audit(
        source_spacing=scan.descriptors[0].spacing,
        spatial_legs=legs,
        local_cell_spacing=local_spacing,
        half_width=half_width,
        minimum_source_cells_per_local_cell=(
            arguments.minimum_source_cells_per_local_cell
        ),
    )
    provenance = _provenance_assessment(manifest)

    metric_cells = arguments.metric_cells_per_axis or cells
    try:
        metric_indices = pilot.selected_time_indices(
            times.size, arguments.metric_cadence_stride
        )
        cadence_error = None
    except ValueError as error:
        metric_indices = []
        cadence_error = str(error)
    duration = (
        float(times[-1] - times[0]) if arguments.tlim is None else arguments.tlim
    )
    if duration > float(times[-1] - times[0]) * (
        1.0 + 128.0 * np.finfo(float).eps
    ):
        raise ValueError("requested duration exceeds the worldtube time table")
    resources = resource_estimate(
        fluid_cells=cells,
        metric_cells=metric_cells,
        mesh_nghost=arguments.mesh_nghost,
        metric_halo=arguments.metric_halo,
        metric_samples=max(1, len(metric_indices)),
        duration=duration,
        cfl=arguments.cfl,
        half_width=half_width,
        fluid_bytes_per_allocated_cell=arguments.fluid_bytes_per_allocated_cell,
    )
    required_halo = pilot.minimum_metric_halo(
        metric_cells, cells, arguments.mesh_nghost
    )

    magnetization = (
        sum(value * value for value in magnetic_field) / rho
        if override
        else float(boundary["maximum_b_squared_over_density_proxy"])
    )
    conditions = {
        "geometry_mass_window": {
            "passed": bool(geometry_window["feasible"]),
            "minimum_secondary_mass": geometry_window["minimum_secondary_mass"],
            "maximum_secondary_mass": geometry_window["maximum_secondary_mass"],
        },
        "direct_bhl_mass_window": {
            "passed": bool(direct_window["feasible"]),
            "minimum_secondary_mass": direct_window["minimum_secondary_mass"],
            "maximum_secondary_mass": direct_window["maximum_secondary_mass"],
        },
        "source_spatial_resolution": {
            "passed": bool(source_audit["passed"]),
            "observed_source_cells_per_local_cell": source_audit[
                "minimum_source_cells_per_local_cell"
            ],
            "required_source_cells_per_local_cell": (
                arguments.minimum_source_cells_per_local_cell
            ),
        },
        "boundary_magnetization_proxy": _condition(
            "maximum",
            magnetization,
            arguments.maximum_boundary_magnetization_proxy,
        ),
        "metric_time_samples": {
            "passed": cadence_error is None,
            "selected_indices": metric_indices,
            "error": cadence_error,
        },
        "metric_halo": _condition(
            "minimum", float(arguments.metric_halo), float(required_halo)
        ),
        "adm_builder_memory": _condition(
            "maximum",
            float(resources["adm_volume"]["python_builder_array_floor_gib"]),
            arguments.maximum_adm_builder_gib,
        ),
        "source_provenance": provenance,
    }
    implementation_ready = all(
        bool(conditions[name]["passed"])
        for name in (
            "geometry_mass_window",
            "source_spatial_resolution",
            "boundary_magnetization_proxy",
            "metric_time_samples",
            "metric_halo",
            "adm_builder_memory",
        )
    )
    science_ready = (
        implementation_ready
        and bool(conditions["direct_bhl_mass_window"]["passed"])
        and bool(provenance["passed"])
    )
    magnetic_amplitude_factor = (
        1.0
        if magnetization <= arguments.maximum_boundary_magnetization_proxy
        else math.sqrt(
            arguments.maximum_boundary_magnetization_proxy / magnetization
        )
    )
    minimum_geometry_cells = minimum_cells_for_geometry(
        arguments.minimum_horizon_cells,
        arguments.minimum_boundary_horizon_radii,
    )
    minimum_direct_cells = minimum_cells_for_direct_capture(
        minimum_horizon_cells=arguments.minimum_horizon_cells,
        minimum_boundary_capture_radii=arguments.minimum_boundary_capture_radii,
        capture_radius_per_mass=capture_radius_per_mass,
        secondary_chi=arguments.secondary_chi,
    )
    suggested_cells = mesh_friendly_cells(
        max(minimum_geometry_cells, minimum_direct_cells)
    )
    suggested_spacing = 2.0 * half_width / suggested_cells
    suggested_source_spacing = (
        float(source_audit["minimum_raw_spatial_jacobian_singular_value"])
        * suggested_spacing
        / arguments.minimum_source_cells_per_local_cell
    )

    def recommended_design(
        design_cells: int, include_capture: bool
    ) -> dict[str, object]:
        design_window = mass_window(
            half_width=half_width,
            cells=design_cells,
            secondary_chi=arguments.secondary_chi,
            minimum_horizon_cells=arguments.minimum_horizon_cells,
            minimum_boundary_horizon_radii=(
                arguments.minimum_boundary_horizon_radii
            ),
            capture_radius_per_mass=(
                capture_radius_per_mass if include_capture else None
            ),
            minimum_boundary_capture_radii=(
                arguments.minimum_boundary_capture_radii
                if include_capture
                else None
            ),
        )
        design_spacing = 2.0 * half_width / design_cells
        matrix = convergence_cost_matrix(
            fluid_cells=design_cells,
            resolution_factors=arguments.metric_resolution_factors,
            mesh_nghost=arguments.mesh_nghost,
            metric_samples=max(1, len(metric_indices)),
            duration=duration,
            cfl=arguments.cfl,
            half_width=half_width,
            fluid_bytes_per_allocated_cell=(
                arguments.fluid_bytes_per_allocated_cell
            ),
        )
        return {
            "face_cells_per_axis": design_cells,
            "local_cell_spacing": design_spacing,
            "secondary_mass_window": design_window,
            "maximum_isotropic_global_source_spacing": (
                float(source_audit[
                    "minimum_raw_spatial_jacobian_singular_value"
                ])
                * design_spacing
                / arguments.minimum_source_cells_per_local_cell
            ),
            "metric_convergence_costs": matrix,
            "fluid_cost": resource_estimate(
                fluid_cells=design_cells,
                metric_cells=design_cells,
                mesh_nghost=arguments.mesh_nghost,
                metric_halo=pilot.minimum_metric_halo(
                    design_cells, design_cells, arguments.mesh_nghost
                ),
                metric_samples=max(1, len(metric_indices)),
                duration=duration,
                cfl=arguments.cfl,
                half_width=half_width,
                fluid_bytes_per_allocated_cell=(
                    arguments.fluid_bytes_per_allocated_cell
                ),
            )["fluid"],
        }
    return {
        "classification": CLASSIFICATION,
        "campaign": str(campaign_path),
        "case": case_name,
        "worldtube": str(worldtube_path),
        "manifest": str(manifest_path),
        "current_design": {
            "face_cells_per_axis": cells,
            "half_width": half_width,
            "local_cell_spacing": local_spacing,
            "duration": duration,
            "secondary_chi": arguments.secondary_chi,
        },
        "upstream_coordinate_proxy": {
            **state,
            "origin": "command_line_override" if override else "campaign_boundary_mean",
            "capture_radius_per_secondary_mass": capture_radius_per_mass,
            "warning": (
                "boundary means are treated as orthonormal components only for this "
                "cost proxy; final BHL scales require the source tetrad and local metric"
            ),
        },
        "geometry_only_mass_window": geometry_window,
        "direct_bhl_mass_window": direct_window,
        "source_resolution": source_audit,
        "source_provenance": provenance,
        "resources": resources,
        "requirements": {
            "minimum_geometry_face_cells": minimum_geometry_cells,
            "minimum_direct_bhl_face_cells_proxy": minimum_direct_cells,
            "minimum_horizon_cells": arguments.minimum_horizon_cells,
            "minimum_boundary_horizon_radii": (
                arguments.minimum_boundary_horizon_radii
            ),
            "minimum_boundary_capture_radii": (
                arguments.minimum_boundary_capture_radii
            ),
            "maximum_boundary_magnetization_proxy": (
                arguments.maximum_boundary_magnetization_proxy
            ),
        },
        "recommended_new_source": {
            "face_cells_per_axis": suggested_cells,
            "maximum_isotropic_global_source_spacing": suggested_source_spacing,
            "maximum_magnetic_amplitude_factor_at_fixed_density": (
                magnetic_amplitude_factor
            ),
            "required_provenance": REQUIRED_SOURCE_PROVENANCE,
            "note": (
                "generate these cells in the global/AMR source itself; resampling an "
                "existing coarse worldtube does not satisfy the source-resolution gate"
            ),
        },
        "recommended_designs": {
            "geometry_only_near_zone": recommended_design(
                mesh_friendly_cells(minimum_geometry_cells), False
            ),
            "direct_bhl_cube_proxy": recommended_design(
                mesh_friendly_cells(minimum_direct_cells), True
            ),
        },
        "assessment": {
            "implementation_geometry_pilot_ready": implementation_ready,
            "direct_bhl_science_ready": science_ready,
            "conditions": conditions,
            "interpretation": (
                "geometry-pilot readiness does not claim a settled BHL wake; direct-BHL "
                "readiness additionally requires capture-radius coverage and the "
                "declared primary-only global provenance contract"
            ),
        },
        "limitations": [
            "the capture radius uses an orthonormalized coordinate-state proxy",
            "the source spatial audit does not replace source-resolution convergence",
            "the ADM builder memory is an array lower bound rather than peak RSS",
            "settling time and force-history convergence remain runtime observables",
        ],
    }


def _format(value: object) -> str:
    if isinstance(value, bool):
        return "PASS" if value else "FAIL"
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float):
        return f"{value:.8g}"
    return str(value)


def markdown_report(plan: dict[str, object]) -> str:
    current = plan["current_design"]
    requirements = plan["requirements"]
    resources = plan["resources"]
    assessment = plan["assessment"]
    geometry = plan["geometry_only_mass_window"]
    direct = plan["direct_bhl_mass_window"]
    source = plan["source_resolution"]
    designs = plan["recommended_designs"]
    assert all(isinstance(value, dict) for value in (
        current,
        requirements,
        resources,
        assessment,
        geometry,
        direct,
        source,
        designs,
    ))
    conditions = assessment["conditions"]
    lines = [
        "# Numerical-ADM EMRI inner replay plan",
        "",
        f"Case: `{plan['case']}`.",
        "",
        "| gate | result |",
        "|---|---:|",
    ]
    for name, condition in conditions.items():
        lines.append(f"| `{name}` | {_format(condition['passed'])} |")
    lines.extend(
        [
            "",
            "## Scale windows",
            "",
            "| quantity | value |",
            "|---|---:|",
            f"| current face cells | {_format(current['face_cells_per_axis'])} |",
            "| geometry minimum face cells | "
            f"{_format(requirements['minimum_geometry_face_cells'])} |",
            "| direct-BHL minimum face cells (proxy) | "
            f"{_format(requirements['minimum_direct_bhl_face_cells_proxy'])} |",
            f"| geometry mass lower bound | {_format(geometry['minimum_secondary_mass'])} |",
            f"| geometry mass upper bound | {_format(geometry['maximum_secondary_mass'])} |",
            f"| direct-BHL mass upper bound | {_format(direct['maximum_secondary_mass'])} |",
            "| source cells per requested local cell | "
            f"{_format(source['minimum_source_cells_per_local_cell'])} |",
            "",
            "## Cost floor",
            "",
            "| quantity | estimate |",
            "|---|---:|",
            f"| active fluid cells | {_format(resources['fluid']['active_cells'])} |",
            "| ghost-inclusive fluid cells | "
            f"{_format(resources['fluid']['ghost_inclusive_allocated_cells'])} |",
            f"| causal CFL steps | {_format(resources['fluid']['causal_cfl_step_proxy'])} |",
            f"| ADM binary GiB | {_format(resources['adm_volume']['binary_gib'])} |",
            "| ADM builder array floor GiB | "
            f"{_format(resources['adm_volume']['python_builder_array_floor_gib'])} |",
            "",
            "Implementation geometry pilot ready: "
            f"**{_format(assessment['implementation_geometry_pilot_ready'])}**.",
            "",
            f"Direct BHL science ready: **{_format(assessment['direct_bhl_science_ready'])}**.",
            "",
            str(assessment["interpretation"]),
        ]
    )
    lines.extend(
        [
            "",
            "## Recommended design costs",
            "",
            "| design | cells | balanced mass | RK3 zone updates | "
            "finest ADM factor | builder GiB |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for name, design in designs.items():
        matrix = design["metric_convergence_costs"]
        finest = max(
            matrix.values(), key=lambda entry: entry["metric_resolution_factor"]
        )
        mass = design["secondary_mass_window"]["balanced_secondary_mass"]
        lines.append(
            f"| `{name}` | {_format(design['face_cells_per_axis'])} | "
            f"{_format(mass)} | "
            f"{_format(design['fluid_cost']['rk3_active_zone_updates_proxy'])} | "
            f"{_format(finest['metric_resolution_factor'])} | "
            f"{_format(finest['python_builder_array_floor_gib'])} |"
        )
    return "\n".join(lines) + "\n"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--case")
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--secondary-chi", type=float, default=0.0)
    parser.add_argument("--adiabatic-index", type=float, default=4.0 / 3.0)
    parser.add_argument("--upstream-rho", type=float)
    parser.add_argument("--upstream-pgas", type=float)
    parser.add_argument(
        "--upstream-four-velocity", type=float, nargs=3, metavar=("U1", "U2", "U3")
    )
    parser.add_argument(
        "--upstream-magnetic-field", type=float, nargs=3, metavar=("B1", "B2", "B3")
    )
    parser.add_argument("--minimum-horizon-cells", type=float, default=4.0)
    parser.add_argument("--minimum-boundary-horizon-radii", type=float, default=5.0)
    parser.add_argument("--minimum-boundary-capture-radii", type=float, default=8.0)
    parser.add_argument(
        "--minimum-source-cells-per-local-cell", type=float, default=1.0
    )
    parser.add_argument(
        "--maximum-boundary-magnetization-proxy", type=float, default=1.0e3
    )
    parser.add_argument("--metric-cells-per-axis", type=int)
    parser.add_argument("--metric-cadence-stride", type=int, default=1)
    parser.add_argument("--metric-halo", type=int, default=4)
    parser.add_argument("--mesh-nghost", type=int, default=4)
    parser.add_argument("--tlim", type=float)
    parser.add_argument("--cfl", type=float, default=0.02)
    parser.add_argument("--fluid-bytes-per-allocated-cell", type=float, default=512.0)
    parser.add_argument("--maximum-adm-builder-gib", type=float, default=16.0)
    parser.add_argument(
        "--metric-resolution-factors", type=float, nargs="+", default=(0.5, 1.0, 2.0, 4.0)
    )
    arguments = parser.parse_args()
    positive = (
        "adiabatic_index",
        "minimum_horizon_cells",
        "minimum_boundary_horizon_radii",
        "minimum_boundary_capture_radii",
        "minimum_source_cells_per_local_cell",
        "maximum_boundary_magnetization_proxy",
        "cfl",
        "fluid_bytes_per_allocated_cell",
        "maximum_adm_builder_gib",
    )
    for name in positive:
        value = float(getattr(arguments, name))
        if not math.isfinite(value) or value <= 0.0:
            parser.error(f"--{name.replace('_', '-')} must be finite and positive")
    if arguments.adiabatic_index <= 1.0:
        parser.error("--adiabatic-index must be greater than one")
    upstream = (
        arguments.upstream_rho,
        arguments.upstream_pgas,
        arguments.upstream_four_velocity,
        arguments.upstream_magnetic_field,
    )
    if any(value is not None for value in upstream) and not all(
        value is not None for value in upstream
    ):
        parser.error(
            "upstream override requires rho, pgas, four-velocity, and magnetic-field"
        )
    if arguments.upstream_rho is not None:
        for name in ("upstream_rho", "upstream_pgas"):
            value = getattr(arguments, name)
            if not math.isfinite(value) or value <= 0.0:
                parser.error(f"--{name.replace('_', '-')} must be finite and positive")
        for name in ("upstream_four_velocity", "upstream_magnetic_field"):
            if not all(math.isfinite(value) for value in getattr(arguments, name)):
                parser.error(f"--{name.replace('_', '-')} must contain finite values")
    if not math.isfinite(arguments.secondary_chi) or abs(arguments.secondary_chi) > 1.0:
        parser.error("--secondary-chi must lie in [-1,1]")
    if arguments.metric_cells_per_axis is not None and arguments.metric_cells_per_axis < 2:
        parser.error("--metric-cells-per-axis must be at least two")
    for name in ("metric_cadence_stride", "metric_halo", "mesh_nghost"):
        if getattr(arguments, name) < 1:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if arguments.tlim is not None and (
        not math.isfinite(arguments.tlim) or arguments.tlim <= 0.0
    ):
        parser.error("--tlim must be finite and positive")
    if (
        not arguments.metric_resolution_factors
        or any(
            not math.isfinite(value) or value <= 0.0
            for value in arguments.metric_resolution_factors
        )
        or len(set(arguments.metric_resolution_factors))
        != len(arguments.metric_resolution_factors)
    ):
        parser.error("--metric-resolution-factors must be distinct and positive")
    return arguments


def main() -> int:
    arguments = parse_arguments()
    plan = build_plan(arguments)
    prefix = arguments.output_prefix.expanduser().resolve()
    json_path = prefix.with_suffix(".json")
    markdown_path = prefix.with_suffix(".md")
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(plan, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    markdown_path.write_text(markdown_report(plan), encoding="utf-8")
    assessment = plan["assessment"]
    print(
        "implementation_geometry_pilot_ready="
        f"{assessment['implementation_geometry_pilot_ready']}"
    )
    print(f"direct_bhl_science_ready={assessment['direct_bhl_science_ready']}")
    print(json_path)
    print(markdown_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
