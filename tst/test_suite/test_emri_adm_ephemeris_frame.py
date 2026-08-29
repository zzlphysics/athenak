"""Tests for numerical-ADM tetrad transport along an EMRI ephemeris."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import build_adm_ephemeris_frame as adm_frame  # noqa: E402
import build_kerr_ephemeris_frame as ephemeris_frame  # noqa: E402
import extract_global_worldtube as extract  # noqa: E402


def _adm_snapshot(
    time: float,
    spatial_scale: float = 1.0,
    lapse_gradient: float = 0.0,
) -> extract.UniformSnapshot:
    shape = (8, 8, 8)
    lower = np.asarray((-2.0, -2.0, -2.0))
    spacing = np.asarray((0.5, 0.5, 0.5))
    values = np.zeros((len(extract.ADM_VARIABLES), *shape))
    values[0] = spatial_scale
    values[3] = spatial_scale
    values[5] = spatial_scale
    x = lower[0] + (np.arange(shape[0]) + 0.5) * spacing[0]
    values[6] = 1.0 + lapse_gradient * x[None, None, :]
    return extract.UniformSnapshot(
        time=time,
        cycle=round(10 * time),
        lower=lower,
        spacing=spacing,
        shape_xyz=shape,
        values=values,
        state_path=Path(f"state-{time}.bin"),
        adm_path=Path(f"adm-{time}.bin"),
    )


def _stationary_ephemeris(
    times: tuple[float, ...] = (0.0, 1.0),
) -> ephemeris_frame.HermiteEphemeris:
    return ephemeris_frame.HermiteEphemeris(
        np.asarray(times),
        np.asarray(((0.75, 0.0, 0.0),) * len(times)),
        np.zeros((len(times), 3)),
    )


def test_adm_metric_time_derivative_matches_finite_difference() -> None:
    values = np.asarray((1.2, 0.1, -0.03, 0.9, 0.04, 1.1, 0.8, 0.02, -0.01, 0.03))
    derivative = np.asarray(
        (0.2, -0.04, 0.01, 0.1, -0.02, 0.3, -0.05, 0.03, 0.02, -0.01)
    )
    analytic = adm_frame.adm_metric_time_derivative(values, derivative)
    step = 1.0e-6
    numerical = (
        adm_frame.adm_metric(values + step * derivative)
        - adm_frame.adm_metric(values - step * derivative)
    ) / (2.0 * step)
    np.testing.assert_allclose(analytic, numerical, rtol=2.0e-9, atol=2.0e-10)


def test_source_snapshot_knot_selects_one_sided_adm_time_derivative() -> None:
    snapshots = extract.SnapshotSeries(
        (
            _adm_snapshot(0.0, 1.0),
            _adm_snapshot(0.5, 1.1),
            _adm_snapshot(1.0, 1.4),
        )
    )
    geometry = adm_frame.ADMSnapshotGeometry(
        snapshots, spatial_fd_step=0.1, audit=adm_frame.GeometryAudit()
    )
    position = np.asarray((0.75, 0.0, 0.0))
    _, _, left = geometry.metric_connection(0.5, position, 0)
    _, _, right = geometry.metric_connection(0.5, position, 1)
    np.testing.assert_allclose(np.diag(left[0])[1:], 0.2, atol=2.0e-15)
    np.testing.assert_allclose(np.diag(right[0])[1:], 0.6, atol=2.0e-15)


def test_parallel_transport_in_time_dependent_homogeneous_adm_metric() -> None:
    snapshots = extract.SnapshotSeries(
        (
            _adm_snapshot(0.0, 1.0),
            _adm_snapshot(0.5, 1.1),
            _adm_snapshot(1.0, 1.2),
        )
    )
    scale = 10.0
    document = adm_frame.build_frame_document(
        _stationary_ephemeris(),
        snapshots,
        transport_mode="parallel",
        global_length_in_local_units=scale,
        metric_fd_step=0.1,
        integration_substeps_per_interval=16,
        geodesic_acceleration_tolerance=1.0e-10,
    )
    assert document["times"] == [0.0, 5.0, 10.0]
    physical_legs = scale * np.asarray(document["spatial_legs"])
    expected_norms = 1.0 / np.sqrt(np.asarray((1.0, 1.1, 1.2)))
    np.testing.assert_allclose(
        physical_legs[:, 1, 0], expected_norms, rtol=2.0e-12, atol=2.0e-14
    )
    np.testing.assert_allclose(
        physical_legs[:, 2, 1], expected_norms, rtol=2.0e-12, atol=2.0e-14
    )
    np.testing.assert_allclose(
        physical_legs[:, 3, 2], expected_norms, rtol=2.0e-12, atol=2.0e-14
    )
    diagnostics = document["diagnostics"]
    assert diagnostics["inserted_source_snapshot_knot_count"] == 1
    assert diagnostics["maximum_dimensionless_proper_acceleration"] < 2.0e-15
    assert diagnostics["maximum_coarse_fine_spatial_leg_difference"] < 2.0e-10
    assert diagnostics["maximum_interpolated_tetrad_gram_error"] < 5.0e-6
    assert diagnostics["maximum_half_step_connection_absolute_difference"] < 2.0e-15


def test_fermi_walker_adm_frame_recovers_rindler_proper_acceleration() -> None:
    gradient = 0.1
    snapshots = extract.SnapshotSeries(
        (
            _adm_snapshot(0.0, lapse_gradient=gradient),
            _adm_snapshot(1.0, lapse_gradient=gradient),
        )
    )
    document = adm_frame.build_frame_document(
        _stationary_ephemeris(),
        snapshots,
        transport_mode="fermi_walker",
        metric_fd_step=0.1,
        integration_substeps_per_interval=8,
    )
    expected = gradient / (1.0 + 0.75 * gradient)
    np.testing.assert_allclose(
        document["diagnostics"]["maximum_dimensionless_proper_acceleration"],
        expected,
        rtol=2.0e-13,
    )
    physical_legs = np.asarray(document["spatial_legs"])
    np.testing.assert_allclose(
        physical_legs,
        np.broadcast_to(physical_legs[0], physical_legs.shape),
        rtol=0.0,
        atol=3.0e-15,
    )


def test_parallel_adm_transport_rejects_accelerated_worldline() -> None:
    snapshots = extract.SnapshotSeries(
        (
            _adm_snapshot(0.0, lapse_gradient=0.1),
            _adm_snapshot(1.0, lapse_gradient=0.1),
        )
    )
    try:
        adm_frame.build_frame_document(
            _stationary_ephemeris(),
            snapshots,
            transport_mode="parallel",
            metric_fd_step=0.1,
            integration_substeps_per_interval=2,
        )
    except ValueError as error:
        assert "geodesic ephemeris in the sampled ADM geometry" in str(error)
    else:
        raise AssertionError("parallel ADM transport accepted accelerated worldline")
