"""Tests for outer-subvolume versus inner-replay closure comparison."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import compare_worldtube_closure as closure  # noqa: E402


def _output(array: np.ndarray, bounds: tuple[float, ...], blocks: int) -> dict:
    nz, ny, nx = array.shape
    variables = closure.DEFAULT_VARIABLES
    if blocks == 1:
        geometries = np.asarray((bounds,))
        values = [array]
    elif blocks == 8:
        geometries = []
        values = []
        for lz in range(2):
            for ly in range(2):
                for lx in range(2):
                    xmid = 0.5 * (bounds[0] + bounds[1])
                    ymid = 0.5 * (bounds[2] + bounds[3])
                    zmid = 0.5 * (bounds[4] + bounds[5])
                    geometries.append(
                        (
                            (bounds[0], xmid)[lx],
                            (xmid, bounds[1])[lx],
                            (bounds[2], ymid)[ly],
                            (ymid, bounds[3])[ly],
                            (bounds[4], zmid)[lz],
                            (zmid, bounds[5])[lz],
                        )
                    )
                    values.append(
                        array[
                            lz * nz // 2 : (lz + 1) * nz // 2,
                            ly * ny // 2 : (ly + 1) * ny // 2,
                            lx * nx // 2 : (lx + 1) * nx // 2,
                        ]
                    )
        geometries = np.asarray(geometries)
    return {
        "time": 0.2,
        "cycle": 3,
        "Nx1": nx,
        "Nx2": ny,
        "Nx3": nz,
        "x1min": bounds[0],
        "x1max": bounds[1],
        "x2min": bounds[2],
        "x2max": bounds[3],
        "x3min": bounds[4],
        "x3max": bounds[5],
        "mb_geometry": geometries,
        "mb_data": {name: [value.copy() for value in values] for name in variables},
    }


def test_closure_extracts_subvolume_and_assembles_eight_blocks() -> None:
    reference_array = np.arange(8**3, dtype=np.float64).reshape(8, 8, 8) + 1.0
    reference = _output(reference_array, (-4.0, 4.0, -4.0, 4.0, -4.0, 4.0), 1)
    candidate_array = reference_array[2:6, 2:6, 2:6].copy()
    candidate_array[1, 1, 1] += 0.5
    candidate = _output(candidate_array, (-2.0, 2.0, -2.0, 2.0, -2.0, 2.0), 8)
    report = closure.compare_loaded_outputs(reference, candidate)
    assert report["candidate_shape"] == [4, 4, 4]
    expected_l2 = 0.5 / np.linalg.norm(reference_array[2:6, 2:6, 2:6])
    np.testing.assert_allclose(
        report["variables"]["dens"]["relative_l2"], expected_l2
    )
    np.testing.assert_allclose(
        report["vector_groups"]["velocity"]["absolute_linf"], 0.5
    )
