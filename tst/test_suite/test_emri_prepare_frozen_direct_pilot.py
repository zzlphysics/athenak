"""Tests for reduced direct frozen-snapshot pilot preparation."""

import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import prepare_frozen_direct_pilot as pilot  # noqa: E402


def test_inflow_faces_follow_oblique_velocity_signs() -> None:
    boundaries = pilot.inflow_boundaries([-0.8, -0.6, 0.001])
    assert boundaries == {
        "ix1_bc": "outflow",
        "ox1_bc": "user",
        "ix2_bc": "outflow",
        "ox2_bc": "user",
        "ix3_bc": "user",
        "ox3_bc": "outflow",
    }


def test_reduced_pilot_has_runnable_mesh_and_source_eos() -> None:
    parameters = {
        name: 0.0 for name in pilot.static.PROFILE_PARAMETER_ORDER
    }
    parameters.update(
        {"rho0": 1.0, "pgas0": 0.003, "u1": -0.04, "u2": -0.03, "b3": 0.001}
    )
    case = {
        "id": "selected",
        "state_sha256": "0" * 64,
        "state_thermodynamics": {"adiabatic_index": 13.0 / 9.0},
        "orbit": {"boyer_lindquist_radius": 56.0},
        "profiles": [{"parameters": parameters}],
        "coherence_scales_local": {"disk_H": 100000.0},
        "assessment": {"local_taylor_model_passed": True},
    }
    campaign = {
        "classification": pilot.frozen.CLASSIFICATION,
        "primary": {
            "mass_global": 1.0,
            "dimensionless_spin": 0.9375,
            "orbit_direction": 1,
        },
        "mass_ratio": 1.0e-5,
        "cases": [case],
    }
    result = pilot.build_pilot(
        campaign,
        "selected",
        cells_per_secondary_mass=8.0,
        outer_cells_per_capture_radius=8.0,
        flow_crossings=2.0,
        meshblock_cells=8,
        calibration_cycles=15,
        finest_refinement_radius=1.0,
    )
    mesh = result["mesh"]
    assert all(value % 8 == 0 for value in mesh["root_dimensions"])
    assert mesh["physical_refinement_levels"] >= 1
    assert len(mesh["refinement_radii"]) == mesh["physical_refinement_levels"]
    assert (
        mesh["refinement_radii"]["1"]
        < mesh["minimum_user_boundary_refinement_clearance"]
    )
    assert any(value == "mhd/gamma=1.4444444444444444"
               for value in result["athena_overrides"])
    assert math.isfinite(result["bhl_plan"]["scales"]["capture_radius_factor_two_for_cost"])
