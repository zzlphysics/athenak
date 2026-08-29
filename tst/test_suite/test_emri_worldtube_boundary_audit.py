"""Tests for characteristic/cadence audits of an EMRI worldtube."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import audit_worldtube_boundary as audit  # noqa: E402
import worldtube_flux_emf as worldtube  # noqa: E402


def _uniform_worldtube() -> tuple[np.ndarray, dict[str, worldtube.FaceData], dict]:
    times = np.asarray((0.0, 0.2))
    cells = 2
    state = np.zeros((2, 8, cells, cells))
    state[:, 0] = 1.0
    state[:, 1] = 0.3
    state[:, 2] = 0.05
    state[:, 3] = -0.02
    state[:, 4] = 0.1
    state[:, 5] = 0.15
    state[:, 6] = 0.2
    state[:, 7] = 0.1
    faces = {}
    for name in worldtube.FACE_NAMES:
        orientation = worldtube.ORIENTATIONS[name]
        normal_field = orientation.normal_sign * state[0, 5 + orientation.normal_axis, 0, 0]
        faces[name] = worldtube.FaceData(
            cell_state=state.copy(),
            normal_flux=np.full((2, cells, cells), normal_field),
            emf_u=np.zeros((1, cells + 1, cells)),
            emf_v=np.zeros((1, cells, cells + 1)),
        )
    return times, faces, {"half_width": 1.0}


def test_audit_classifies_all_faces_and_reports_hll_proxy() -> None:
    times, faces, metadata = _uniform_worldtube()
    report = audit.audit_worldtube(times, faces, metadata, 4.0 / 3.0)
    assert report["sample_count"] == 48
    assert report["eigensystem_failures"] == 0
    assert set(report["faces"]) == set(worldtube.FACE_NAMES)
    assert 0.0 < report["mixed_fan_fraction"] <= 1.0
    assert report["maximum_worldtube_cadence_courant"] < 0.2
    assert report["linear_hll_flux_gain_error"]["maximum"] > 0.0


def test_audit_stride_reduces_spatial_and_temporal_samples() -> None:
    times, faces, metadata = _uniform_worldtube()
    report = audit.audit_worldtube(
        times, faces, metadata, 4.0 / 3.0, sample_stride=2
    )
    assert report["sample_count"] == 6
