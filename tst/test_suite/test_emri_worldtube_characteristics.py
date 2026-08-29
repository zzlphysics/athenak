"""Tests for the local-orthonormal EMRI worldtube characteristic projector."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import worldtube_characteristics as characteristics  # noqa: E402
import worldtube_flux_emf as worldtube  # noqa: E402


def test_static_srmhd_conserved_state_and_flux() -> None:
    density = 1.2
    pressure = 0.15
    gamma = 4.0 / 3.0
    magnetic = np.asarray((0.3, -0.2, 0.1))
    primitive = np.asarray(
        (density, pressure, 0.0, 0.0, 0.0, magnetic[1], magnetic[2])
    )
    conserved, flux = characteristics.srmhd_conserved_flux(
        primitive, magnetic[0], gamma
    )
    expected_tau = pressure / (gamma - 1.0) + 0.5 * np.dot(magnetic, magnetic)
    expected_momentum_flux = -magnetic[0] * magnetic
    expected_momentum_flux[0] += pressure + 0.5 * np.dot(magnetic, magnetic)
    np.testing.assert_allclose(conserved[0], density, atol=2.0e-15)
    np.testing.assert_allclose(conserved[1:4], 0.0, atol=2.0e-15)
    np.testing.assert_allclose(conserved[4], expected_tau, atol=2.0e-15)
    np.testing.assert_allclose(flux[1:4], expected_momentum_flux, atol=2.0e-15)
    np.testing.assert_allclose(flux[[0, 4, 5, 6]], 0.0, atol=2.0e-15)


def test_unmagnetized_characteristics_recover_relativistic_sound_speeds() -> None:
    density = 1.0
    pressure = 0.1
    gamma = 4.0 / 3.0
    normal_u = 0.3
    primitive = np.asarray((density, pressure, normal_u, 0.0, 0.0, 0.0, 0.0))
    basis = characteristics.characteristic_basis(primitive, 0.0, gamma)
    lorentz = np.sqrt(1.0 + normal_u**2)
    velocity = normal_u / lorentz
    enthalpy_density = density + gamma * pressure / (gamma - 1.0)
    sound_speed = np.sqrt(gamma * pressure / enthalpy_density)
    expected_minus = (velocity - sound_speed) / (1.0 - velocity * sound_speed)
    expected_plus = (velocity + sound_speed) / (1.0 + velocity * sound_speed)
    np.testing.assert_allclose(basis.speeds[0], expected_minus, atol=2.0e-8)
    np.testing.assert_allclose(basis.speeds[-1], expected_plus, atol=2.0e-8)
    np.testing.assert_allclose(basis.speeds[1:-1], velocity, atol=2.0e-8)
    assert basis.jacobian_residual < 2.0e-13


def test_superfast_inflow_sets_exterior_and_outflow_keeps_interior() -> None:
    interior_in = np.asarray((1.0, 0.1, -3.0, 0.05, -0.02, 0.2, 0.1))
    exterior_in = np.asarray((1.1, 0.12, -3.2, 0.02, -0.01, 0.18, 0.13))
    boundary_in, diagnostics_in = characteristics.project_incoming_characteristics(
        interior_in, exterior_in, 0.15, 4.0 / 3.0
    )
    assert np.all(diagnostics_in.incoming)
    np.testing.assert_allclose(boundary_in, exterior_in, atol=2.0e-14)

    interior_out = interior_in.copy()
    exterior_out = exterior_in.copy()
    interior_out[2] *= -1.0
    exterior_out[2] *= -1.0
    boundary_out, diagnostics_out = characteristics.project_incoming_characteristics(
        interior_out, exterior_out, 0.15, 4.0 / 3.0
    )
    assert not np.any(diagnostics_out.incoming)
    np.testing.assert_array_equal(boundary_out, interior_out)


def test_mixed_boundary_replaces_only_incoming_linearized_amplitudes() -> None:
    reference = np.asarray((1.0, 0.1, 0.3, 0.05, -0.02, 0.2, 0.1))
    difference = np.asarray((0.0, 0.0, 0.02, -0.01, 0.015, 0.008, -0.006))
    interior = reference - 0.5 * difference
    exterior = reference + 0.5 * difference
    boundary, diagnostics = characteristics.project_incoming_characteristics(
        interior, exterior, 0.15, 4.0 / 3.0
    )
    basis = characteristics.characteristic_basis(reference, 0.15, 4.0 / 3.0)
    applied_amplitudes = basis.left_eigenvectors @ (boundary - interior)
    target_amplitudes = basis.left_eigenvectors @ difference
    np.testing.assert_allclose(
        applied_amplitudes[diagnostics.incoming],
        target_amplitudes[diagnostics.incoming],
        atol=2.0e-14,
    )
    np.testing.assert_allclose(
        applied_amplitudes[~diagnostics.incoming], 0.0, atol=2.0e-14
    )
    assert 0 < np.count_nonzero(diagnostics.incoming) < 7


def test_all_cubical_face_rotations_round_trip_and_preserve_normal_field() -> None:
    state = np.asarray((1.0, 0.3, -0.2, 0.1, 0.08, 0.4, 0.05, -0.07))
    normal_magnetic = 0.23
    for name in worldtube.FACE_NAMES:
        primitive = characteristics.state_to_face_primitive(
            state, name, normal_magnetic
        )
        recovered = characteristics.face_primitive_to_state(
            primitive, name, normal_magnetic
        )
        basis = characteristics.face_basis(name)
        np.testing.assert_allclose(recovered[:5], state[:5], atol=2.0e-15)
        np.testing.assert_allclose(
            basis[0] @ recovered[5:8], normal_magnetic, atol=2.0e-15
        )
        np.testing.assert_allclose(
            (basis @ recovered[5:8])[1:], (basis @ state[5:8])[1:], atol=2.0e-15
        )


def test_source_state_wrapper_uses_outward_incoming_sign() -> None:
    interior = np.asarray((1.0, 3.0, 0.05, -0.02, 0.1, 0.15, 0.2, 0.1))
    exterior = np.asarray((1.1, 3.2, 0.02, -0.01, 0.12, 0.15, 0.18, 0.13))
    # On x1m the outward normal is -x1, so positive source u1 is superfast inflow.
    boundary, diagnostics = characteristics.characteristic_boundary_state(
        interior, exterior, "x1m", 0.15, 4.0 / 3.0
    )
    assert np.all(diagnostics.incoming)
    expected = exterior.copy()
    expected[5] = -0.15
    np.testing.assert_allclose(boundary, expected, atol=3.0e-14)
