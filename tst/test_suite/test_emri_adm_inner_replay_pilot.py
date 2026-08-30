"""Tests for the fail-closed ADM inner replay pilot."""

from argparse import Namespace
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import run_adm_inner_replay_pilot as pilot  # noqa: E402
import worldtube_flux_emf as worldtube  # noqa: E402


def _affine_faces(cells: int = 4) -> tuple[dict, dict]:
    names = ["rho", "u1", "u2", "u3", "pgas", "bcc1", "bcc2", "bcc3"]
    center = np.zeros(3)
    half_width = 2.0
    faces = {}
    for name in worldtube.FACE_NAMES:
        positions = pilot._face_positions(cells, center, half_width, name)
        state = np.empty((2, len(names), cells, cells))
        values = np.column_stack(
            (
                np.exp(0.2 + positions @ np.asarray((0.1, -0.2, 0.05))),
                0.3 + positions @ np.asarray((0.02, 0.01, -0.03)),
                -0.1 + positions @ np.asarray((0.04, 0.0, 0.01)),
                positions @ np.asarray((0.0, -0.02, 0.03)),
                np.exp(-2.0 + positions @ np.asarray((-0.03, 0.02, 0.04))),
                np.full(positions.shape[0], 0.01),
                np.full(positions.shape[0], -0.02),
                np.full(positions.shape[0], 0.03),
            )
        )
        state[0] = values.T.reshape(len(names), cells, cells)
        state[1] = state[0]
        faces[name] = worldtube.FaceData(
            state,
            np.zeros((2, cells, cells)),
            np.zeros((1, cells + 1, cells)),
            np.zeros((1, cells, cells + 1)),
        )
    return faces, {
        "state_variables": names,
        "center": center.tolist(),
        "half_width": half_width,
    }


def test_affine_profile_fit_recovers_boundary_model() -> None:
    faces, metadata = _affine_faces()
    fit = pilot.fit_initial_affine_profile(faces, metadata)
    assert np.isclose(fit["rho0"], np.exp(0.2))
    assert np.isclose(fit["pgas0"], np.exp(-2.0))
    np.testing.assert_allclose(fit["u"], (0.3, -0.1, 0.0), atol=1.0e-15)
    np.testing.assert_allclose(
        fit["log_density_gradient"], (0.1, -0.2, 0.05), atol=1.0e-15
    )
    assert max(
        entry["relative_l2"] for entry in fit["fit_residuals"].values()
    ) < 1.0e-14
    floors = pilot._floor_controls(
        fit,
        {
            "minimum_density": 1.0,
            "maximum_density": 2.0,
            "minimum_pressure": 0.1,
        },
        2.0,
    )
    assert 0.0 < floors["density_floor"] < 1.0e-4
    assert 0.0 < floors["temperature_floor"] < 1.0e-5


def test_structural_assessment_requires_resolution_and_separation() -> None:
    arguments = Namespace(
        secondary_mass=0.04,
        secondary_chi=0.0,
        minimum_horizon_cells=4.0,
        minimum_boundary_horizon_radii=5.0,
        maximum_boundary_magnetization_proxy=1000.0,
        maximum_initial_volume_flux_divergence=1.0e-10,
    )
    case = {
        "boundary_state": {"maximum_b_squared_over_density_proxy": 2.0},
        "inner_validation": {
            "initial_volume_flux": {
                "maximum_relative_cell_flux_divergence": 1.0e-14
            }
        },
    }
    assessment = pilot.structural_assessment(arguments, case, 64, 0.5)
    assert assessment["passed"]
    unresolved = pilot.structural_assessment(arguments, case, 16, 0.5)
    assert not unresolved["passed"]
    assert not unresolved["conditions"]["secondary_horizon_resolution"]["passed"]
