"""Unit tests for the time-dependent EMRI Taylor-table builder."""

from pathlib import Path
import sys
import tempfile

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import build_taylor_worldtube_series as series  # noqa: E402
import extract_static_taylor_worldtube as static  # noqa: E402


def _parameters(scale: float) -> dict[str, float]:
    return {
        name: scale * (index + 1)
        for index, name in enumerate(static.PROFILE_PARAMETER_ORDER)
    }


def test_taylor_table_round_trip() -> None:
    samples = [
        {"time": 2.0, "parameters": _parameters(0.25)},
        {"time": 3.5, "parameters": _parameters(-0.5)},
    ]
    samples[0]["parameters"]["rho0"] = 1.0
    samples[0]["parameters"]["pgas0"] = 0.1
    samples[1]["parameters"]["rho0"] = 2.0
    samples[1]["parameters"]["pgas0"] = 0.2
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "profiles.dat"
        series.write_table(
            path,
            samples,
            source_coordinate_radius_local_units=20.0,
            source_coordinate_angular_frequency_local_units=0.02,
        )
        text = path.read_text(encoding="utf-8")
        times, values = series.read_table(path)
    assert "# athenak-emri-taylor-series-v2" in text
    assert "# source_coordinate_radius_local_units: 20" in text
    np.testing.assert_array_equal(times, (2.0, 3.5))
    expected = np.asarray(
        [
            [sample["parameters"][name] for name in static.PROFILE_PARAMETER_ORDER]
            for sample in samples
        ]
    )
    np.testing.assert_array_equal(values, expected)


def test_circular_orbit_contract() -> None:
    samples = [
        {
            "anchor_global": [10.0, 0.0, 0.0],
            "source_coordinate_velocity": [0.0, 2.0, 0.0],
        },
        {
            "anchor_global": [0.0, 10.0, 0.0],
            "source_coordinate_velocity": [-2.0, 0.0, 0.0],
        },
    ]
    diagnostics = series.circular_orbit_diagnostics(
        samples, np.zeros(3), np.asarray((0.0, 0.0, 1.0))
    )
    series.validate_circular_orbit(diagnostics, 1.0e-13)
    np.testing.assert_allclose(
        diagnostics["mean_coordinate_angular_frequency"], 0.2, atol=1.0e-15
    )


def test_noncircular_orbit_is_rejected() -> None:
    samples = [
        {
            "anchor_global": [10.0, 0.0, 0.0],
            "source_coordinate_velocity": [0.1, 2.0, 0.0],
        },
        {
            "anchor_global": [0.0, 10.1, 0.0],
            "source_coordinate_velocity": [-2.0, 0.0, 0.0],
        },
    ]
    diagnostics = series.circular_orbit_diagnostics(
        samples, np.zeros(3), np.asarray((0.0, 0.0, 1.0))
    )
    try:
        series.validate_circular_orbit(diagnostics, 1.0e-4)
    except RuntimeError as error:
        assert "circular-equatorial" in str(error)
    else:
        raise AssertionError("noncircular worldline was accepted")
