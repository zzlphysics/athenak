"""Unit tests for frozen real-GRMHD Taylor campaign planning."""

import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import build_frozen_taylor_campaign as campaign  # noqa: E402


def _parameters() -> dict[str, float]:
    result = {name: 0.0 for name in campaign.static.PROFILE_PARAMETER_ORDER}
    result.update(
        {
            "rho0": 1.0,
            "pgas0": 0.01,
            "u1": 0.1,
            "b2": 0.02,
            "dlnrho_dxh1": 0.001,
            "dlnpgas_dxh2": 0.002,
            "du1_dxh1": 0.0001,
            "db2_dxh3": 0.0002,
        }
    )
    return result


def test_case_orbit_is_circular_and_outside_isco() -> None:
    orbit = campaign.case_orbit(1.0, 0.9375, 14.0, 1, math.pi / 4.0)
    anchor = orbit["anchor"]
    velocity = orbit["source_velocity"]
    assert math.isclose(math.hypot(anchor[0], anchor[1]),
                        orbit["cartesian_equatorial_radius"])
    assert math.isclose(anchor[0] * velocity[0] + anchor[1] * velocity[1], 0.0,
                        abs_tol=1.0e-15)
    assert orbit["boyer_lindquist_radius"] > orbit["kerr_isco_radius"]


def test_gradient_coherence_scales_and_source_gamma() -> None:
    scales, state = campaign.gradient_coherence_scales(
        _parameters(), 13.0 / 9.0, 1.0e-4
    )
    assert math.isclose(scales["L_rho"], 1000.0)
    assert math.isclose(scales["L_pgas"], 500.0)
    assert math.isclose(scales["L_velocity"], 1000.0)
    assert math.isclose(scales["L_magnetic"], 100.0)
    assert state["adiabatic_index"] == 13.0 / 9.0
    assert scales["disk_H_hydrostatic_proxy"] > 0.0


def test_fitting_scale_comparison_uses_joint_vector_norms() -> None:
    reference = _parameters()
    changed = dict(reference)
    changed["b2"] *= 1.1
    maximum, detail = campaign.fitting_scale_comparison(
        [
            {"fit_radius_source": 1.0, "parameters": reference},
            {"fit_radius_source": 2.0, "parameters": changed},
        ]
    )
    magnetic = detail["r2"]["groups"]["magnetic_field"]
    assert maximum == magnetic["symmetric_relative_l2"]
    assert math.isclose(maximum, 0.1 / 1.05)


def test_request_rejects_duplicate_case_ids() -> None:
    request = {
        "classification": campaign.REQUEST_CLASSIFICATION,
        "primary": {"mass": 1.0, "dimensionless_spin": 0.5},
        "mass_ratio": 1.0e-5,
        "density_renormalization": 1.0e10,
        "cases": [
            {
                "id": "same",
                "state": __file__,
                "orbital_radius": 10.0,
                "phase": 0.0,
                "fit_radii": [0.5, 1.0],
            },
            {
                "id": "same",
                "state": __file__,
                "orbital_radius": 12.0,
                "phase": 1.0,
                "fit_radii": [0.5, 1.0],
            },
        ],
    }
    try:
        campaign._validated_request(request, Path(__file__).resolve())
    except ValueError as error:
        assert "duplicated" in str(error)
    else:
        raise AssertionError("duplicate case ids were accepted")
