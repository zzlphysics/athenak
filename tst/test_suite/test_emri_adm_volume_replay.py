"""Tests for numerical-background plus secondary ADM volume extraction."""

from pathlib import Path
import sys
import tempfile

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import adm_volume_replay as replay  # noqa: E402


def _minkowski_series(times: np.ndarray, nodes: int = 3) -> np.ndarray:
    metric = np.zeros((times.size, nodes, nodes, nodes, 4, 4))
    metric[..., 0, 0] = -1.0
    metric[..., 1, 1] = 1.0
    metric[..., 2, 2] = 1.0
    metric[..., 3, 3] = 1.0
    return metric


def test_extrinsic_curvature_vanishes_for_minkowski() -> None:
    times = np.asarray((0.0, 0.5, 1.0))
    metric = _minkowski_series(times)
    curvature, diagnostics = replay.extrinsic_curvature(
        metric, times, np.ones(3)
    )
    assert np.max(np.abs(curvature)) == 0.0
    assert diagnostics["minimum_lapse"] == 1.0
    assert diagnostics["minimum_spatial_metric_eigenvalue"] == 1.0


def test_extrinsic_curvature_recovers_linear_scale_factor() -> None:
    times = np.asarray((0.0, 0.5, 1.0))
    metric = _minkowski_series(times)
    scale = 1.0 + 0.1 * times
    for index, value in enumerate(scale):
        for axis in range(1, 4):
            metric[index, ..., axis, axis] = value**2
    curvature, _ = replay.extrinsic_curvature(metric, times, np.ones(3))
    for index, value in enumerate(scale):
        expected = -0.1 * value
        for axis in range(3):
            np.testing.assert_allclose(
                curvature[index, ..., axis, axis], expected, atol=2.0e-15
            )


def test_metric_derivative_decomposition_recovers_linear_scale_factor() -> None:
    times = np.asarray((0.0, 0.5, 1.0))
    metric = _minkowski_series(times)
    derivatives = np.zeros(metric.shape[:-2] + (4, 4, 4))
    scale = 1.0 + 0.1 * times
    for index, value in enumerate(scale):
        for axis in range(1, 4):
            metric[index, ..., axis, axis] = value**2
            derivatives[index, ..., 0, axis, axis] = 0.2 * value
    _, curvature = replay.decompose_metric_derivatives(metric, derivatives)
    for index, value in enumerate(scale):
        expected = -0.1 * value
        for axis in range(3):
            np.testing.assert_allclose(
                curvature[index, ..., axis, axis], expected, atol=2.0e-15
            )


def test_spatial_metric_derivatives_recover_linear_metric() -> None:
    times = np.asarray((0.0, 1.0))
    metric = _minkowski_series(times, nodes=4)
    z, y, x = np.meshgrid(
        np.arange(4.0), np.arange(4.0), np.arange(4.0), indexing="ij"
    )
    metric[..., 0, 0] -= x + 2.0 * y + 3.0 * z
    derivatives = replay.spatial_metric_derivatives(metric, np.ones(3))
    np.testing.assert_allclose(derivatives[..., 0, 0, 0], -1.0)
    np.testing.assert_allclose(derivatives[..., 1, 0, 0], -2.0)
    np.testing.assert_allclose(derivatives[..., 2, 0, 0], -3.0)


def test_secondary_kerr_term_is_symmetric_and_tapered() -> None:
    positions = np.asarray(((0.0, 0.0, 0.0), (3.0, 0.2, -0.1)))
    perturbation = replay.secondary_kerr_perturbation(
        positions, 1.0, 0.3, np.eye(4)
    )
    assert np.max(np.abs(perturbation[0])) == 0.0
    np.testing.assert_allclose(
        perturbation, perturbation.transpose(0, 2, 1), atol=0.0
    )
    assert np.max(np.abs(perturbation[1])) > 0.0


def test_adm_volume_binary_round_trip_and_checksum() -> None:
    times = np.asarray((0.0, 0.5, 1.0))
    metric = _minkowski_series(times)
    curvature, _ = replay.extrinsic_curvature(metric, times, np.ones(3))
    fields = np.empty((3, len(replay.FIELD_NAMES), 3, 3, 3))
    for field, (left, right) in enumerate(replay.METRIC_COMPONENTS):
        fields[:, field] = metric[..., left, right]
    for offset, (left, right) in enumerate(replay.CURVATURE_COMPONENTS):
        fields[:, len(replay.METRIC_COMPONENTS) + offset] = curvature[
            ..., left, right
        ]
    volume = replay.ADMVolume(
        times,
        np.asarray((-1.0, -1.0, -1.0)),
        np.ones(3),
        fields,
        {"classification": replay.CLASSIFICATION},
        np.broadcast_to(np.eye(4), (times.size, 4, 4)).copy(),
    )
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "adm.bin"
        replay.write_binary(path, volume)
        loaded = replay.read_binary(path)
        np.testing.assert_array_equal(loaded.times, times)
        np.testing.assert_array_equal(loaded.fields, fields)
        np.testing.assert_array_equal(
            loaded.secondary_coframes, volume.secondary_coframes
        )
        assert loaded.metadata["binary_version"] == replay.BINARY_VERSION
        corrupted = bytearray(path.read_bytes())
        corrupted[-1] ^= 1
        bad = Path(directory) / "bad.bin"
        bad.write_bytes(corrupted)
        try:
            replay.read_binary(bad)
        except ValueError as error:
            assert "checksum" in str(error)
        else:
            raise AssertionError("corrupt ADM volume binary was accepted")

        legacy_path = Path(directory) / "legacy.bin"
        legacy = replay.ADMVolume(
            times,
            volume.lower,
            volume.spacing,
            fields,
            volume.metadata,
        )
        replay.write_binary(legacy_path, legacy)
        loaded_legacy = replay.read_binary(legacy_path)
        assert loaded_legacy.secondary_coframes is None
        assert loaded_legacy.metadata["binary_version"] == replay.LEGACY_BINARY_VERSION


def test_hybrid_adm_volume_binary_round_trip() -> None:
    times = np.asarray((0.0, 0.5, 1.0))
    fields = np.zeros((3, len(replay.HYBRID_FIELD_NAMES), 3, 3, 3))
    fields[:, 0] = -1.0
    for field in (4, 7, 9):
        fields[:, field] = 1.0
    parameters = np.asarray((0.01, 0.3, 0.05, 1.0e-3, 5.0e-7))
    volume = replay.ADMVolume(
        times,
        np.asarray((-1.0, -1.0, -1.0)),
        np.ones(3),
        fields,
        {"classification": replay.CLASSIFICATION},
        np.broadcast_to(np.eye(4), (times.size, 4, 4)).copy(),
        parameters,
    )
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "hybrid-adm.bin"
        sidecar = replay.write_binary(path, volume)
        loaded = replay.read_binary(path)
        np.testing.assert_array_equal(loaded.times, times)
        np.testing.assert_array_equal(loaded.fields, fields)
        np.testing.assert_array_equal(loaded.hybrid_parameters, parameters)
        assert sidecar["classification"] == replay.HYBRID_BINARY_CLASSIFICATION
        assert loaded.metadata["binary_version"] == replay.HYBRID_BINARY_VERSION
        assert loaded.metadata["field_names"] == list(replay.HYBRID_FIELD_NAMES)
