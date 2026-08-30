"""Tests for Taylor fitting-scale and cadence convergence analysis."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import analyze_taylor_profile_convergence as convergence  # noqa: E402


def test_group_error_uses_joint_vector_normalization() -> None:
    reference = np.asarray(((3.0, 4.0), (3.0, 4.0)))
    candidate = np.asarray(((3.0, 5.0), (3.0, 5.0)))
    error = convergence.group_error(
        reference, candidate, np.asarray((0, 1))
    )
    np.testing.assert_allclose(error["absolute_rms"], 1.0)
    np.testing.assert_allclose(error["reference_rms"], 5.0)
    np.testing.assert_allclose(error["relative_l2"], 0.2)


def test_cadence_selection_retains_last_sample_and_interpolates() -> None:
    times = np.asarray((0.0, 1.0, 2.0, 3.0, 4.0, 5.0))
    values = np.column_stack((times**2, times))
    selected, candidate = convergence.cadence_candidate(times, values, 2)
    np.testing.assert_array_equal(selected, (0, 2, 4, 5))
    np.testing.assert_array_equal(candidate[selected], values[selected])
    np.testing.assert_allclose(candidate[1, 0], 2.0)
    np.testing.assert_allclose(candidate[3, 0], 10.0)
