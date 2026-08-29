"""Tests for conservative EMRI cubical flux/EMF worldtube transfer."""

import json
from pathlib import Path
import sys
import tempfile

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import worldtube_flux_emf as worldtube  # noqa: E402


def _constant_cube(
    cells: int,
    times: np.ndarray,
    magnetic: tuple[float, float, float] = (0.3, -0.2, 0.1),
    electric: tuple[float, float, float] = (0.04, -0.03, 0.02),
) -> dict[str, worldtube.FaceData]:
    faces = {}
    area = (2.0 / cells) ** 2
    length = 2.0 / cells
    for name in worldtube.FACE_NAMES:
        orientation = worldtube.ORIENTATIONS[name]
        state = np.ones((times.size, 2, cells, cells))
        normal_value = (
            orientation.normal_sign * magnetic[orientation.normal_axis] * area
        )
        normal_flux = np.full((times.size, cells, cells), normal_value)
        emf_u = np.full(
            (times.size - 1, cells + 1, cells),
            orientation.u_sign * electric[orientation.u_axis] * length,
        )
        emf_v = np.full(
            (times.size - 1, cells, cells + 1),
            orientation.v_sign * electric[orientation.v_axis] * length,
        )
        faces[name] = worldtube.FaceData(state, normal_flux, emf_u, emf_v)
    return faces


def test_face_regrid_commutes_with_discrete_faraday() -> None:
    generator = np.random.default_rng(1847)
    times = np.asarray((0.0, 0.17))
    nv, nu = 4, 6
    initial_flux = generator.normal(size=(nv, nu))
    emf_u = generator.normal(size=(nv + 1, nu))
    emf_v = generator.normal(size=(nv, nu + 1))
    final_flux = worldtube.faraday_update(
        initial_flux, emf_u, emf_v, times[1] - times[0]
    )
    face = worldtube.FaceData(
        cell_state=generator.normal(size=(2, 3, nv, nu)),
        normal_flux=np.stack((initial_flux, final_flux)),
        emf_u=emf_u[None, ...],
        emf_v=emf_v[None, ...],
    )
    transferred = worldtube.resample_face(face, target_nv=7, target_nu=5)
    residual = worldtube.faraday_residuals(transferred, times)
    assert float(np.max(np.abs(residual))) < 2.0e-14
    np.testing.assert_allclose(
        np.sum(transferred.normal_flux, axis=(1, 2)),
        np.sum(face.normal_flux, axis=(1, 2)),
        rtol=0.0,
        atol=2.0e-14,
    )


def test_constant_cube_has_shared_edges_and_zero_closed_flux() -> None:
    times = np.asarray((1.0, 1.25, 1.75))
    faces = _constant_cube(4, times)
    diagnostics = worldtube.validate_worldtube(times, faces)
    assert diagnostics["maximum_shared_edge_emf_residual"] == 0.0
    assert diagnostics["maximum_closed_surface_flux"] < 1.0e-15
    assert max(diagnostics["maximum_faraday_residual_by_face"].values()) == 0.0


def test_worldtube_round_trip_and_noninteger_regrid() -> None:
    times = np.asarray((0.0, 0.2))
    faces = _constant_cube(3, times)
    with tempfile.TemporaryDirectory() as directory:
        source = Path(directory) / "source.npz"
        target = Path(directory) / "target.npz"
        worldtube.write_worldtube(source, times, faces, {"case": "constant"})
        loaded_times, loaded_faces, metadata = worldtube.read_worldtube(source)
        transferred = worldtube.resample_worldtube(loaded_times, loaded_faces, 7)
        worldtube.write_worldtube(target, loaded_times, transferred, metadata)
        final_times, final_faces, final_metadata = worldtube.read_worldtube(target)
    np.testing.assert_array_equal(final_times, times)
    assert final_metadata["case"] == "constant"
    diagnostics = worldtube.validate_worldtube(final_times, final_faces)
    assert diagnostics["maximum_closed_surface_flux"] < 2.0e-15


def test_mismatched_duplicate_cube_edge_is_rejected() -> None:
    times = np.asarray((0.0, 0.5))
    faces = _constant_cube(4, times)
    faces["x1p"].emf_u[0, 0, 0] += 0.1
    try:
        worldtube.validate_worldtube(times, faces)
    except ValueError as error:
        assert "edge" in str(error)
    else:
        raise AssertionError("inconsistent duplicated cube edge was accepted")


def test_outer_stream_binary_manifest_round_trip() -> None:
    times = np.asarray((2.0, 2.25))
    cells = 3
    faces = _constant_cube(cells, times)
    with tempfile.TemporaryDirectory() as directory_name:
        directory = Path(directory_name)
        times_file = directory / "times.bin"
        np.asarray(times, dtype="<f8").tofile(times_file)
        file_table = {}
        for name, face in faces.items():
            file_table[name] = {}
            for field_name in ("cell_state", "normal_flux", "emf_u", "emf_v"):
                filename = f"{name}.{field_name}.bin"
                np.asarray(getattr(face, field_name), dtype="<f8").tofile(
                    directory / filename
                )
                file_table[name][field_name] = filename
        manifest = {
            "classification": worldtube.OUTER_STREAM_CLASSIFICATION,
            "target_classification": worldtube.CLASSIFICATION,
            "complete": True,
            "binary_dtype": "<f8",
            "times_file": times_file.name,
            "nt": 2,
            "ninterval": 1,
            "nvar": 2,
            "cells_per_face_axis": cells,
            "state_variables": ["a", "b"],
            "faces": file_table,
        }
        manifest_path = directory / "stream.manifest.json"
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        loaded_times, loaded_faces, metadata = worldtube.read_outer_stream(
            manifest_path
        )
    np.testing.assert_array_equal(loaded_times, times)
    assert metadata["state_variables"] == ["a", "b"]
    diagnostics = worldtube.validate_worldtube(loaded_times, loaded_faces)
    assert diagnostics["maximum_shared_edge_emf_residual"] == 0.0
