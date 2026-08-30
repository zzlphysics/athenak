"""Unit tests for the EMRI BHL hierarchy planner."""

import json
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import plan_bhl_hierarchy as planner  # noqa: E402


def _slow_wind(**overrides: object) -> dict[str, object]:
    inputs: dict[str, object] = {
        "secondary_mass": 1.0,
        "primary_mass": 100_000.0,
        "orbital_radius": 1_000_000.0,
        "rho": 1.0,
        "pgas": 1.0e-4,
        "gamma": 4.0 / 3.0,
        "spatial_four_velocity": (0.05, 0.0, 0.0),
        "magnetic_field": (0.01, 0.0, 0.0),
        "coherence_scales": {"disk_H": 100_000.0},
    }
    inputs.update(overrides)
    return inputs


def test_characteristic_state_uses_spatial_four_velocity() -> None:
    state = planner.characteristic_state(
        1.0, 0.1, 4.0 / 3.0, (3.0, 0.0, 0.0), (2.0, 0.0, 0.0)
    )
    assert math.isclose(float(state["lorentz_factor"]), math.sqrt(10.0))
    assert math.isclose(float(state["three_speed"]), 3.0 / math.sqrt(10.0))
    assert math.isclose(float(state["specific_enthalpy"]), 1.4)
    assert math.isclose(float(state["comoving_magnetic_b_squared"]), 4.0)


def test_large_coherent_capture_scale_selects_matching() -> None:
    plan = planner.build_plan(**_slow_wind())
    assert plan["recommendation"] == "matched_bhl_outer_inner"
    assert plan["scales"]["clean_overlap"] is True
    assert plan["validity"]["uniform_bhl_environment"] is True
    assert plan["validity"]["direct_within_budget"] is False
    assert plan["validity"]["matched_components_within_budget"] is True
    assert (
        plan["costs"]["matched_outer"]["finest_steps"]
        < plan["costs"]["direct"]["finest_steps"]
    )


def test_disk_scale_caps_uniform_bhl_before_cost_decision() -> None:
    plan = planner.build_plan(
        **_slow_wind(coherence_scales={"disk_H": 100.0})
    )
    assert (
        plan["recommendation"]
        == "global_or_shearing_outer_with_relativistic_inner"
    )
    assert plan["scales"]["limiting_environment_scale_name"] == "disk_H"
    assert plan["validity"]["uniform_bhl_environment"] is False


def test_fast_compact_case_remains_direct() -> None:
    plan = planner.build_plan(
        **_slow_wind(
            pgas=0.01,
            spatial_four_velocity=(1.0, 0.0, 0.0),
            magnetic_field=(0.0, 0.0, 0.0),
            coherence_scales={},
        )
    )
    assert plan["recommendation"] == "direct_grmhd"
    assert plan["validity"]["direct_within_budget"] is True


def test_oblique_wind_expands_fixed_axis_domain() -> None:
    settings = planner.PlannerSettings()
    aligned = planner.axis_aligned_domain_envelope(settings, (1.0, 0.0, 0.0))
    oblique = planner.axis_aligned_domain_envelope(settings, (1.0, 1.0, 0.0))
    assert aligned["widths_in_capture_radii"] == [12.0, 8.0, 8.0]
    assert oblique["widths_in_capture_radii"][0] > 12.0
    assert oblique["widths_in_capture_radii"][1] > 8.0
    assert oblique["widths_in_capture_radii"][2] == 8.0
    assert math.isclose(
        sum(value * value for value in oblique["wind_unit_vector"]), 1.0
    )


def test_oblique_root_grid_rounds_to_meshblock_multiple() -> None:
    plan = planner.build_plan(
        **_slow_wind(spatial_four_velocity=(0.05, 0.05, 0.0))
    )
    dimensions = plan["costs"]["direct"]["base_grid_dimensions"]
    assert all(value % 8 == 0 for value in dimensions)


def test_zero_field_plan_is_strict_json() -> None:
    plan = planner.build_plan(
        **_slow_wind(magnetic_field=(0.0, 0.0, 0.0))
    )
    assert plan["upstream_state"]["plasma_beta"] is None
    encoded = json.dumps(plan, allow_nan=False)
    assert planner.CLASSIFICATION in encoded
    report = planner.markdown_report(plan)
    assert "face-integrated magnetic flux" in report
    assert "count accreted momentum exactly once" in report
