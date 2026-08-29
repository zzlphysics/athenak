"""Tests for mode-resolved worldtube reflection diagnostics."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import analyze_worldtube_reflection as reflection  # noqa: E402
import worldtube_characteristics as characteristics  # noqa: E402


def test_opposite_speed_eigenmode_is_measured_as_reflection() -> None:
    background = np.asarray((4.0, 1.0, -0.106, 0.318, -0.053, 1.8, -1.2))
    normal_field = 2.5
    gamma = 4.0 / 3.0
    basis = characteristics.characteristic_basis(background, normal_field, gamma)
    phase = np.sin(np.linspace(-0.8, 0.9, 40))
    source_mode = 6
    reflected_mode = 0
    reference = background[:, None] + basis.right_eigenvectors[:, source_mode, None] * phase
    candidate = (
        reference
        + 0.1 * basis.right_eigenvectors[:, reflected_mode, None] * phase
    )
    report = reflection.measure_mode_content(
        reference, candidate, background, normal_field, gamma, source_mode
    )
    np.testing.assert_allclose(report["reflected_amplitude_coefficient"], 0.1)
    np.testing.assert_allclose(report["mode_error_coefficients"][0], 0.1)
    assert report["reference_non_source_leakage"] < 2.0e-15


def test_same_direction_error_is_not_counted_as_reflection() -> None:
    background = np.asarray((4.0, 1.0, -0.106, 0.318, -0.053, 1.8, -1.2))
    normal_field = 2.5
    gamma = 4.0 / 3.0
    basis = characteristics.characteristic_basis(background, normal_field, gamma)
    phase = np.sin(np.linspace(-1.1, 0.7, 32))
    reference = background[:, None] + basis.right_eigenvectors[:, 6, None] * phase
    candidate = reference + 0.25 * basis.right_eigenvectors[:, 5, None] * phase
    report = reflection.measure_mode_content(
        reference, candidate, background, normal_field, gamma, 6
    )
    assert report["reflected_amplitude_coefficient"] < 2.0e-15
    np.testing.assert_allclose(report["mode_error_coefficients"][5], 0.25)
