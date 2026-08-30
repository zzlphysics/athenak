"""Tests for the resolved numerical-ADM inner-replay planner."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import plan_adm_inner_replay as planner  # noqa: E402


def test_geometry_window_proves_forty_cell_lower_bound() -> None:
    exact = planner.mass_window(
        half_width=4.0,
        cells=40,
        secondary_chi=0.0,
        minimum_horizon_cells=4.0,
        minimum_boundary_horizon_radii=5.0,
    )
    assert exact["feasible"]
    assert np.isclose(exact["minimum_secondary_mass"], 0.4)
    assert np.isclose(exact["maximum_secondary_mass"], 0.4)
    assert planner.minimum_cells_for_geometry(4.0, 5.0) == 40

    unresolved = planner.mass_window(
        half_width=4.0,
        cells=39,
        secondary_chi=0.0,
        minimum_horizon_cells=4.0,
        minimum_boundary_horizon_radii=5.0,
    )
    assert not unresolved["feasible"]


def test_capture_radius_expands_direct_bhl_requirement() -> None:
    # r_capture/m=4, r_H/m=2, four horizon cells and eight capture radii
    # require N >= 2*4*8*4/2 = 128 cells per cube edge.
    assert planner.minimum_cells_for_direct_capture(
        minimum_horizon_cells=4.0,
        minimum_boundary_capture_radii=8.0,
        capture_radius_per_mass=4.0,
        secondary_chi=0.0,
    ) == 128
    assert planner.mesh_friendly_cells(40) == 40
    assert planner.mesh_friendly_cells(301) == 304
    window = planner.mass_window(
        half_width=4.0,
        cells=127,
        secondary_chi=0.0,
        minimum_horizon_cells=4.0,
        minimum_boundary_horizon_radii=5.0,
        capture_radius_per_mass=4.0,
        minimum_boundary_capture_radii=8.0,
    )
    assert not window["feasible"]


def test_source_resolution_rejects_target_oversampling() -> None:
    identity_legs = np.asarray(
        (
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
        )
    )
    resolved = planner.source_resolution_audit(
        source_spacing=(0.2, 0.2, 0.2),
        spatial_legs=(identity_legs,),
        local_cell_spacing=0.2,
        half_width=4.0,
        minimum_source_cells_per_local_cell=1.0,
    )
    assert resolved["passed"]
    assert np.isclose(resolved["minimum_source_cells_per_local_cell"], 1.0)

    oversampled = planner.source_resolution_audit(
        source_spacing=(1.0, 1.0, 1.0),
        spatial_legs=(identity_legs,),
        local_cell_spacing=0.2,
        half_width=4.0,
        minimum_source_cells_per_local_cell=1.0,
    )
    assert not oversampled["passed"]
    assert np.isclose(oversampled["minimum_source_cells_per_local_cell"], 0.2)


def test_resource_estimate_separates_fluid_and_adm_costs() -> None:
    estimate = planner.resource_estimate(
        fluid_cells=40,
        metric_cells=20,
        mesh_nghost=4,
        metric_halo=4,
        metric_samples=4,
        duration=0.08,
        cfl=0.02,
        half_width=4.0,
        fluid_bytes_per_allocated_cell=512.0,
    )
    assert estimate["fluid"]["active_cells"] == 40**3
    assert estimate["fluid"]["meshblock_count"] == 64
    assert estimate["adm_volume"]["nodes_per_axis"] == 28
    assert estimate["adm_volume"]["binary_gib"] > 0.0
    assert (
        estimate["adm_volume"]["python_builder_array_floor_gib"]
        > estimate["adm_volume"]["binary_gib"]
    )
    matrix = planner.convergence_cost_matrix(
        fluid_cells=40,
        resolution_factors=(0.5, 1.0, 2.0),
        mesh_nghost=4,
        metric_samples=4,
        duration=0.08,
        cfl=0.02,
        half_width=4.0,
        fluid_bytes_per_allocated_cell=512.0,
    )
    assert matrix["0.5"]["metric_cells_per_axis"] == 20
    assert matrix["2"]["minimum_metric_halo"] == 8
    assert matrix["2"]["binary_gib"] > matrix["1"]["binary_gib"]


def test_primary_only_provenance_is_mandatory() -> None:
    missing = planner._provenance_assessment({})
    assert not missing["passed"]
    verified = planner._provenance_assessment(
        {"source_provenance": dict(planner.REQUIRED_SOURCE_PROVENANCE)}
    )
    assert verified["passed"]
