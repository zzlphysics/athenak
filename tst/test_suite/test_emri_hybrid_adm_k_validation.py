"""Tests for the manufactured hybrid ADM-v3 K_ij validation driver."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import adm_volume_replay as adm_volume  # noqa: E402
import run_hybrid_adm_k_validation as validation  # noqa: E402


def test_minkowski_manufactured_volume_has_exact_zero_primary_derivatives() -> None:
    volume = validation.manufactured_volume(
        np.asarray((0.0, 0.1)),
        half_width=0.5,
        metric_cells=4,
        metric_halo=2,
        secondary_mass=0.05,
        secondary_chi=0.4,
        fd_ratio=5.0e-5,
        affine_primary=False,
    )
    assert volume.fields.shape == (2, len(adm_volume.HYBRID_FIELD_NAMES), 8, 8, 8)
    np.testing.assert_array_equal(volume.fields[:, 0], -1.0)
    for field in (4, 7, 9):
        np.testing.assert_array_equal(volume.fields[:, field], 1.0)
    np.testing.assert_array_equal(volume.fields[:, 10:], 0.0)
    np.testing.assert_allclose(
        volume.secondary_coframes,
        np.broadcast_to(np.eye(4), volume.secondary_coframes.shape),
        atol=0.0,
    )
    assert volume.hybrid_parameters[4] == 2.5e-6


def test_affine_primary_stored_derivatives_match_direct_difference() -> None:
    coordinate = np.asarray((0.12, -0.08, 0.04))
    time = 0.03
    _, derivative = validation._primary_metric(
        time, coordinate[0], coordinate[1], coordinate[2], True
    )
    step = 1.0e-5
    for direction in range(3):
        lower = coordinate.copy()
        upper = coordinate.copy()
        lower[direction] -= step
        upper[direction] += step
        lower_metric, _ = validation._primary_metric(
            time, lower[0], lower[1], lower[2], True
        )
        upper_metric, _ = validation._primary_metric(
            time, upper[0], upper[1], upper[2], True
        )
        numerical = (upper_metric - lower_metric) / (2.0 * step)
        np.testing.assert_allclose(
            derivative[direction], numerical, rtol=6.0e-10, atol=6.0e-12
        )


def test_group_error_uses_joint_tensor_norm() -> None:
    reference = {"a": np.asarray((3.0, 4.0)), "b": np.asarray((0.0,))}
    candidate = {"a": np.asarray((3.0, 5.0)), "b": np.asarray((0.0,))}
    error = validation._group_error(reference, candidate, ("a", "b"))
    np.testing.assert_allclose(error["relative_l2"], 0.2)
    np.testing.assert_allclose(error["relative_linf"], 0.25)
    np.testing.assert_allclose(error["maximum_absolute"], 1.0)
