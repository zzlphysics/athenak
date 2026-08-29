"""Tests for numerical Kerr ephemerides and tetrad transport."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import build_kerr_circular_frame as circular  # noqa: E402
import build_kerr_ephemeris_frame as ephemeris_frame  # noqa: E402


def _circular_ephemeris(
    spacing: float, final_time: float = 0.2
) -> ephemeris_frame.HermiteEphemeris:
    mass = 1.0
    spin = 0.5
    radius = 8.0
    coordinate_radius = np.sqrt(radius**2 + spin**2)
    omega = circular.circular_kerr_omega(mass, spin, radius, 1)
    times = np.arange(0.0, final_time + 0.5 * spacing, spacing)
    phase = omega * times
    positions = np.column_stack(
        (
            coordinate_radius * np.cos(phase),
            coordinate_radius * np.sin(phase),
            np.zeros_like(phase),
        )
    )
    velocities = np.column_stack(
        (
            -omega * coordinate_radius * np.sin(phase),
            omega * coordinate_radius * np.cos(phase),
            np.zeros_like(phase),
        )
    )
    return ephemeris_frame.HermiteEphemeris(times, positions, velocities)


def _maximum_proper_acceleration(
    ephemeris: ephemeris_frame.HermiteEphemeris,
) -> float:
    result = 0.0
    for interval in range(ephemeris.global_times.size - 1):
        time = 0.5 * (
            ephemeris.global_times[interval]
            + ephemeris.global_times[interval + 1]
        )
        position, velocity, acceleration = ephemeris.evaluate(time, interval)
        kinematics = ephemeris_frame.worldline_kinematics(
            position,
            velocity,
            acceleration,
            np.zeros(3),
            primary_mass=1.0,
            primary_spin=0.5,
            metric_fd_step=1.0e-5,
        )
        result = max(result, kinematics.proper_acceleration)
    return result


def test_hermite_ephemeris_recovers_a_cubic_and_its_acceleration() -> None:
    times = np.asarray((0.0, 1.0))
    positions = np.column_stack((times**3, 2.0 * times**3, -times**3))
    velocities = np.column_stack((3.0 * times**2, 6.0 * times**2, -3.0 * times**2))
    ephemeris = ephemeris_frame.HermiteEphemeris(
        times, positions, velocities
    )
    position, velocity, acceleration = ephemeris.evaluate(0.3)
    np.testing.assert_allclose(position, (0.027, 0.054, -0.027), atol=2.0e-16)
    np.testing.assert_allclose(velocity, (0.27, 0.54, -0.27), atol=3.0e-16)
    np.testing.assert_allclose(acceleration, (1.8, 3.6, -1.8), atol=5.0e-16)

    scale = 1.0e-6
    discontinuous_acceleration = ephemeris_frame.HermiteEphemeris(
        np.asarray((0.0, 1.0, 2.0)),
        np.zeros((3, 3)),
        np.asarray(((0.0, 0.0, 0.0), (scale, 0.0, 0.0), (0.0, 0.0, 0.0))),
    )
    jump = discontinuous_acceleration.acceleration_jump_diagnostics()
    np.testing.assert_allclose(
        jump["maximum_coordinate_acceleration_jump"], 8.0 * scale
    )
    np.testing.assert_allclose(
        jump["maximum_relative_coordinate_acceleration_jump"], 2.0
    )


def test_numerical_kerr_connection_is_symmetric_and_metric_compatible() -> None:
    metric, connection, derivative = ephemeris_frame.kerr_connection(
        np.asarray((8.0, 1.0, 0.5)),
        np.zeros(3),
        primary_mass=1.0,
        primary_spin=0.4,
        metric_fd_step=1.0e-5,
    )
    np.testing.assert_allclose(
        connection, connection.transpose(0, 2, 1), rtol=0.0, atol=2.0e-18
    )
    residual = np.empty((4, 4, 4))
    for lam in range(4):
        for mu in range(4):
            for nu in range(4):
                residual[lam, mu, nu] = (
                    derivative[lam, mu, nu]
                    - connection[:, lam, mu] @ metric[:, nu]
                    - connection[:, lam, nu] @ metric[mu, :]
                )
    assert np.max(np.abs(residual)) < 3.0e-17


def test_fermi_walker_keeps_a_stationary_kerr_observer_frame_constant() -> None:
    ephemeris = ephemeris_frame.HermiteEphemeris(
        np.asarray((0.0, 1.0, 2.0)),
        np.asarray(((10.0, 0.0, 0.0),) * 3),
        np.zeros((3, 3)),
    )
    scale = 1.0e5
    document = ephemeris_frame.build_frame_document(
        ephemeris,
        primary_mass=1.0,
        dimensionless_spin=0.0,
        global_length_in_local_units=scale,
        integration_substeps_per_interval=8,
    )
    physical_legs = scale * np.asarray(document["spatial_legs"])
    np.testing.assert_allclose(
        physical_legs,
        np.broadcast_to(physical_legs[0], physical_legs.shape),
        rtol=0.0,
        atol=3.0e-15,
    )
    diagnostics = document["diagnostics"]
    expected_acceleration = 1.0 / (100.0 * np.sqrt(0.8))
    np.testing.assert_allclose(
        diagnostics["maximum_dimensionless_proper_acceleration"],
        expected_acceleration,
        rtol=2.0e-10,
    )
    assert diagnostics["maximum_interpolated_tetrad_gram_error"] < 3.0e-15
    assert diagnostics["maximum_coarse_fine_spatial_leg_difference"] < 3.0e-15
    assert document["times"] == [0.0, 1.0e5, 2.0e5]


def test_parallel_transport_rejects_an_accelerated_worldline() -> None:
    ephemeris = ephemeris_frame.HermiteEphemeris(
        np.asarray((0.0, 1.0)),
        np.asarray(((10.0, 0.0, 0.0),) * 2),
        np.zeros((2, 3)),
    )
    try:
        ephemeris_frame.build_frame_document(
            ephemeris,
            primary_mass=1.0,
            dimensionless_spin=0.0,
            transport_mode="parallel",
            integration_substeps_per_interval=2,
        )
    except ValueError as error:
        assert "parallel transport requires a geodesic ephemeris" in str(error)
    else:
        raise AssertionError("parallel transport accepted an accelerated observer")


def test_circular_ephemeris_geodesic_residual_converges_with_cadence() -> None:
    coarse = _maximum_proper_acceleration(_circular_ephemeris(0.1))
    fine = _maximum_proper_acceleration(_circular_ephemeris(0.05))
    assert fine < coarse / 3.8
    assert fine < 6.0e-9


def test_parallel_transport_accepts_a_resolved_circular_geodesic() -> None:
    ephemeris = _circular_ephemeris(0.02, final_time=0.1)
    document = ephemeris_frame.build_frame_document(
        ephemeris,
        primary_mass=1.0,
        dimensionless_spin=0.5,
        transport_mode="parallel",
        integration_substeps_per_interval=4,
        geodesic_acceleration_tolerance=1.0e-8,
    )
    diagnostics = document["diagnostics"]
    assert diagnostics["maximum_dimensionless_proper_acceleration"] < 2.0e-9
    assert diagnostics["maximum_interpolated_tetrad_gram_error"] < 4.0e-12
    assert diagnostics["maximum_coarse_fine_spatial_leg_difference"] < 2.0e-11
