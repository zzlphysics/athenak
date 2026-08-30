"""Tests for state-only static-Kerr worldline manifest generation."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import build_static_kerr_worldline as worldline  # noqa: E402


def test_circular_samples_follow_kerr_orbit_and_sort_dump_times() -> None:
    records = [
        (Path("late.bin"), 2.0, 20),
        (Path("early.bin"), 0.0, 10),
        (Path("middle.bin"), 1.0, 15),
    ]
    samples, orbit = worldline.circular_samples(
        records,
        primary_mass=1.0,
        dimensionless_spin=0.5,
        orbital_radius=10.0,
        direction=1,
        initial_phase=0.2,
    )
    assert [sample["time"] for sample in samples] == [0.0, 1.0, 2.0]
    assert [sample["cycle"] for sample in samples] == [10, 15, 20]
    radius = orbit["cartesian_equatorial_radius"]
    for sample in samples:
        position = np.asarray(sample["anchor"])
        velocity = np.asarray(sample["source_velocity"])
        np.testing.assert_allclose(np.linalg.norm(position), radius)
        np.testing.assert_allclose(position @ velocity, 0.0, atol=2.0e-16)
        np.testing.assert_allclose(
            np.linalg.norm(velocity),
            abs(orbit["coordinate_angular_frequency"]) * radius,
        )
