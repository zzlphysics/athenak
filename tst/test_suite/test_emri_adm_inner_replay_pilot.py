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
    coframes = np.broadcast_to(np.eye(4), (2, 4, 4)).copy()
    coframes[:, 1, 1] = 0.5
    volume = pilot.adm_volume.ADMVolume(
        np.asarray((0.0, 1.0)),
        np.zeros(3),
        np.ones(3),
        np.zeros((2, len(pilot.adm_volume.FIELD_NAMES), 2, 2, 2)),
        {},
        coframes,
    )
    pilot.apply_numerical_coframe_assessment(
        assessment, volume, arguments, 0.5
    )
    assert np.isclose(assessment["maximum_secondary_coordinate_stretch"], 2.0)
    assert not assessment["conditions"]["worldtube_separation"]["passed"]
    assert pilot.minimum_metric_halo(4, 4, 3) == 3
    assert pilot.minimum_metric_halo(8, 4, 3) == 6
    assert pilot.selected_time_indices(4, 2) == [0, 2, 3]
    try:
        pilot.selected_time_indices(4, 3)
    except ValueError as error:
        assert "fewer than three" in str(error)
    else:
        raise AssertionError("two-sample K_ij table was accepted")


def test_expected_adm_fields_reproduces_time_interpolation_and_decomposition() -> None:
    cells = 2
    halo = 1
    nodes = cells + 2 * halo
    times = np.asarray((5.0, 7.0))
    fields = np.zeros((2, len(pilot.adm_volume.FIELD_NAMES), nodes, nodes, nodes))
    metric_left = np.asarray(
        (
            (-1.8, 0.2, -0.1, 0.05),
            (0.2, 1.5, 0.1, 0.0),
            (-0.1, 0.1, 1.2, 0.02),
            (0.05, 0.0, 0.02, 0.9),
        )
    )
    metric_right = metric_left.copy()
    metric_right[0, 0] -= 0.4
    metric_right[0, 1] += 0.1
    metric_right[1, 0] += 0.1
    metric_right[2, 2] += 0.2
    for time_index, metric in enumerate((metric_left, metric_right)):
        for field, (left, right) in enumerate(pilot.adm_volume.METRIC_COMPONENTS):
            fields[time_index, field] = metric[left, right]
        for offset in range(len(pilot.adm_volume.CURVATURE_COMPONENTS)):
            fields[time_index, 10 + offset] = time_index + 0.1 * offset
    volume = pilot.adm_volume.ADMVolume(
        times,
        np.full(3, -1.5),
        np.ones(3),
        fields,
        {},
    )
    expected = pilot._expected_adm_fields(volume, 6.0, cells, 1.0)
    midpoint = 0.5 * (metric_left + metric_right)
    decomposed, _ = pilot.adm_volume.decompose_four_metric(midpoint)
    np.testing.assert_allclose(expected["adm_alpha"], decomposed["alpha"])
    np.testing.assert_allclose(expected["adm_betax"], decomposed["beta"][0])
    np.testing.assert_allclose(expected["adm_gzz"], midpoint[3, 3])
    np.testing.assert_allclose(expected["adm_Kyz"], 0.9)
    assert all(expected[name].shape == (cells, cells, cells) for name in expected)
