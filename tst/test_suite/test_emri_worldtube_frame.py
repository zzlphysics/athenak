"""Tests for moving/source-frame EMRI worldtube geometry and CT projection."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import worldtube_flux_emf as worldtube  # noqa: E402
import worldtube_frame as frame  # noqa: E402


def _static_frame(
    velocity: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> frame.AffineFrame:
    tangent = np.asarray((1.0, *velocity))
    legs = np.zeros((4, 3))
    legs[1:, :] = np.eye(3)
    return frame.AffineFrame(tangent, legs, np.zeros((4, 3)))


def _constant_cube(cells: int, times: np.ndarray) -> dict[str, worldtube.FaceData]:
    magnetic = np.asarray((0.3, -0.2, 0.1))
    electric = np.asarray((0.04, -0.03, 0.02))
    area = (2.0 / cells) ** 2
    length = 2.0 / cells
    faces = {}
    for name in worldtube.FACE_NAMES:
        orientation = worldtube.ORIENTATIONS[name]
        faces[name] = worldtube.FaceData(
            cell_state=np.ones((times.size, 5, cells, cells)),
            normal_flux=np.full(
                (times.size, cells, cells),
                orientation.normal_sign
                * magnetic[orientation.normal_axis]
                * area,
            ),
            emf_u=np.full(
                (times.size - 1, cells + 1, cells),
                orientation.u_sign * electric[orientation.u_axis] * length,
            ),
            emf_v=np.full(
                (times.size - 1, cells, cells + 1),
                orientation.v_sign * electric[orientation.v_axis] * length,
            ),
        )
    return faces


def test_faraday_pullback_contains_motional_emf() -> None:
    electric = np.asarray((0.2, -0.1, 0.05))
    magnetic = np.asarray((0.4, 0.3, -0.2))
    velocity = np.asarray((0.07, -0.11, 0.09))
    pulled_electric, pulled_magnetic = frame.pullback_electric_magnetic(
        electric, magnetic, _static_frame(tuple(velocity))
    )
    np.testing.assert_allclose(
        pulled_electric, electric + np.cross(velocity, magnetic), atol=2.0e-15
    )
    np.testing.assert_allclose(pulled_magnetic, magnetic, atol=2.0e-15)


def test_rotating_affine_frame_uses_position_dependent_surface_velocity() -> None:
    angular_velocity = np.asarray((0.0, 0.0, 0.17))
    leg_derivative = np.zeros((4, 3))
    for column in range(3):
        leg_derivative[1:, column] = np.cross(
            angular_velocity, np.eye(3)[:, column]
        )
    rotating = frame.AffineFrame(
        np.asarray((1.0, 0.0, 0.0, 0.0)),
        np.vstack((np.zeros(3), np.eye(3))),
        leg_derivative,
    )
    position = np.asarray((0.8, -0.3, 0.2))
    magnetic = np.asarray((0.2, -0.4, 0.7))
    pulled_electric, _ = frame.pullback_electric_magnetic(
        np.zeros(3), magnetic, rotating, position
    )
    surface_velocity = np.cross(angular_velocity, position)
    np.testing.assert_allclose(
        pulled_electric, np.cross(surface_velocity, magnetic), atol=2.0e-15
    )


def test_four_vector_round_trip_through_spacetime_jacobian() -> None:
    affine = _static_frame((0.2, -0.1, 0.05))
    jacobian = affine.jacobian()
    local = np.asarray((1.7, -0.2, 0.3, 0.4))
    global_vector = jacobian @ local
    np.testing.assert_allclose(
        frame.transform_contravariant_vector(global_vector, jacobian),
        local,
        atol=2.0e-15,
    )


def test_frame_audit_separates_static_moving_and_cut_surface_cases() -> None:
    static = frame.audit_frame(_static_frame())
    moving = frame.audit_frame(_static_frame((0.1, 0.0, 0.0)))
    angle = 0.2
    rotated_legs = np.zeros((4, 3))
    rotated_legs[1:, :] = np.asarray(
        (
            (np.cos(angle), -np.sin(angle), 0.0),
            (np.sin(angle), np.cos(angle), 0.0),
            (0.0, 0.0, 1.0),
        )
    )
    cut = frame.audit_frame(
        frame.AffineFrame(
            np.asarray((1.0, 0.0, 0.0, 0.0)),
            rotated_legs,
            np.zeros((4, 3)),
        )
    )
    assert static["capability"] == "exact_static_cube_relabel"
    assert static["current_fixed_writer_supported"] is True
    assert moving["capability"] == "moving_axis_aligned_ale_required"
    assert moving["motional_emf_required"] is True
    assert cut["capability"] == "moving_cut_surface_required"
    assert cut["current_fixed_writer_supported"] is False


def test_surface_complex_has_boundary_of_boundary_zero() -> None:
    complex_ = frame.CubeSurfaceComplex(5)
    generator = np.random.default_rng(8162)
    edge_values = generator.normal(size=complex_.edge_count)
    circulation = complex_.curl(edge_values)
    assert abs(float(np.sum(circulation))) < 2.0e-14


def test_projection_repairs_sampled_edges_without_changing_flux() -> None:
    times = np.asarray((0.0, 0.2, 0.55))
    faces = _constant_cube(5, times)
    original_flux = {
        name: face_data.normal_flux.copy() for name, face_data in faces.items()
    }
    generator = np.random.default_rng(9374)
    for face_data in faces.values():
        face_data.emf_u += 3.0e-3 * generator.normal(size=face_data.emf_u.shape)
        face_data.emf_v += 3.0e-3 * generator.normal(size=face_data.emf_v.shape)
    projected, diagnostics = frame.project_moving_samples(times, faces)
    validated = worldtube.validate_worldtube(times, projected)
    for name in worldtube.FACE_NAMES:
        np.testing.assert_array_equal(projected[name].normal_flux, original_flux[name])
    assert validated["maximum_shared_edge_emf_residual"] == 0.0
    assert max(validated["maximum_faraday_residual_by_face"].values()) < 2.0e-14
    assert diagnostics["intervals"][0]["maximum_edge_correction"] > 0.0
    assert diagnostics["intervals"][0]["sampled_maximum_edge"] > 0.0
    expected_relative = (
        diagnostics["intervals"][0]["maximum_edge_correction"]
        / diagnostics["intervals"][0]["sampled_maximum_edge"]
    )
    assert np.isclose(
        diagnostics["intervals"][0]["relative_maximum_edge_correction"],
        expected_relative,
    )
    assert diagnostics["intervals"][0]["final_maximum_faraday_residual"] < 2.0e-14


def test_projection_rejects_changing_closed_surface_flux() -> None:
    times = np.asarray((0.0, 0.3))
    faces = _constant_cube(3, times)
    faces["x1p"].normal_flux[1, 0, 0] += 0.01
    try:
        frame.project_moving_samples(times, faces)
    except ValueError as error:
        assert "closed-surface flux changes" in str(error)
    else:
        raise AssertionError("an incompatible endpoint flux change was accepted")


def test_endpoint_projection_removes_only_closed_flux_mean() -> None:
    times = np.asarray((0.0, 0.3, 0.8))
    faces = _constant_cube(4, times)
    faces["x1p"].normal_flux[0, 1, 2] += 0.03
    faces["x2m"].normal_flux[1, 0, 0] -= 0.02
    original_states = {
        name: values.cell_state.copy() for name, values in faces.items()
    }
    projected, diagnostics = frame.project_closed_surface_fluxes(times, faces)
    for endpoint in range(times.size):
        net = sum(
            np.sum(projected[name].normal_flux[endpoint])
            for name in worldtube.FACE_NAMES
        )
        assert abs(float(net)) < 3.0e-16
    for name in worldtube.FACE_NAMES:
        np.testing.assert_array_equal(projected[name].cell_state, original_states[name])
        delta = projected[name].normal_flux - faces[name].normal_flux
        for endpoint in range(times.size):
            np.testing.assert_allclose(
                delta[endpoint], delta[endpoint, 0, 0], rtol=0.0, atol=5.0e-18
            )
    assert diagnostics["endpoints"][0]["raw_closed_surface_flux"] != 0.0
