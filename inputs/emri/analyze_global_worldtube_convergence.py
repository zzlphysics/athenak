#!/usr/bin/env python3
"""Compare an offline global-snapshot worldtube with an online RK reference."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np

import worldtube_flux_emf as worldtube


CLASSIFICATION = "athenak-emri-global-worldtube-convergence-v1"


def _error_norms(
    reference: np.ndarray, candidate: np.ndarray
) -> dict[str, float | None]:
    left = np.asarray(reference, dtype=np.float64)
    right = np.asarray(candidate, dtype=np.float64)
    if (
        left.shape != right.shape
        or not np.isfinite(left).all()
        or not np.isfinite(right).all()
    ):
        raise ValueError("worldtube comparison arrays are incompatible or non-finite")
    difference = right - left
    reference_l1 = float(np.sum(np.abs(left)))
    reference_l2 = float(np.linalg.norm(left.ravel()))
    reference_linf = float(np.max(np.abs(left)))
    absolute_l1 = float(np.sum(np.abs(difference)))
    absolute_l2 = float(np.linalg.norm(difference.ravel()))
    absolute_linf = float(np.max(np.abs(difference)))
    return {
        "relative_l1": absolute_l1 / reference_l1 if reference_l1 > 0.0 else None,
        "relative_l2": absolute_l2 / reference_l2 if reference_l2 > 0.0 else None,
        "relative_linf": (
            absolute_linf / reference_linf if reference_linf > 0.0 else None
        ),
        "absolute_l1": absolute_l1,
        "absolute_l2": absolute_l2,
        "absolute_linf": absolute_linf,
    }


def match_endpoint_indices(
    reference_times: Iterable[float], target_times: Iterable[float]
) -> np.ndarray:
    reference = worldtube.validate_times(reference_times)
    target = worldtube.validate_times(target_times)
    indices = []
    for time in target:
        insertion = int(np.searchsorted(reference, time))
        candidates = [
            index
            for index in (insertion - 1, insertion)
            if 0 <= index < reference.size
        ]
        if not candidates:
            raise ValueError("candidate time lies outside the reference time range")
        index = min(candidates, key=lambda candidate: abs(reference[candidate] - time))
        scale = max(1.0, abs(float(time)), abs(float(reference[index])))
        if not math.isclose(
            float(reference[index]),
            float(time),
            rel_tol=0.0,
            abs_tol=128.0 * np.finfo(float).eps * scale,
        ):
            raise ValueError(
                f"candidate time {time:.17g} is absent from the online reference"
            )
        indices.append(index)
    result = np.asarray(indices, dtype=np.int64)
    if result[0] != 0 or result[-1] != reference.size - 1:
        raise ValueError("candidate must retain the full reference time range")
    if np.any(np.diff(result) <= 0):
        raise ValueError("matched candidate times do not increase strictly")
    return result


def _concatenate(
    faces: dict[str, worldtube.FaceData], field: str
) -> np.ndarray:
    return np.concatenate(
        [np.asarray(getattr(faces[name], field)).ravel() for name in worldtube.FACE_NAMES]
    )


def compare_worldtubes(
    reference_times: Iterable[float],
    reference_faces: dict[str, worldtube.FaceData],
    candidate_times: Iterable[float],
    candidate_faces: dict[str, worldtube.FaceData],
    state_variables: Iterable[str] | None = None,
) -> dict[str, object]:
    """Compare at candidate endpoints after exact reference EMF aggregation."""

    reference_times_array = worldtube.validate_times(reference_times)
    candidate_times_array = worldtube.validate_times(candidate_times)
    indices = match_endpoint_indices(reference_times_array, candidate_times_array)
    coarsened_times, coarsened = worldtube.coarsen_worldtube_time(
        reference_times_array, reference_faces, indices
    )
    if not np.array_equal(coarsened_times, candidate_times_array):
        raise RuntimeError("matched reference times changed during coarsening")
    checked_candidate = {
        name: worldtube.validate_face(
            candidate_faces[name], candidate_times_array, name
        )
        for name in worldtube.FACE_NAMES
    }
    reference_shape = coarsened[worldtube.FACE_NAMES[0]].cell_state.shape
    for name in worldtube.FACE_NAMES:
        if coarsened[name].cell_state.shape != checked_candidate[name].cell_state.shape:
            raise ValueError("candidate and reference face-state shapes differ")
        if coarsened[name].normal_flux.shape != checked_candidate[name].normal_flux.shape:
            raise ValueError("candidate and reference face-flux shapes differ")

    variable_count = reference_shape[1]
    names = list(state_variables or ())
    if names and len(names) != variable_count:
        raise ValueError("state variable metadata does not match the stored state")
    if not names:
        names = [f"state_{index}" for index in range(variable_count)]
    state_by_variable = {}
    for variable, name in enumerate(names):
        reference = np.concatenate(
            [
                coarsened[face].cell_state[:, variable].ravel()
                for face in worldtube.FACE_NAMES
            ]
        )
        candidate = np.concatenate(
            [
                checked_candidate[face].cell_state[:, variable].ravel()
                for face in worldtube.FACE_NAMES
            ]
        )
        state_by_variable[name] = _error_norms(reference, candidate)

    interval_errors = []
    for interval, dt in enumerate(np.diff(candidate_times_array)):
        reference_emf = np.concatenate(
            [
                coarsened[name].emf_u[interval].ravel()
                for name in worldtube.FACE_NAMES
            ]
            + [
                coarsened[name].emf_v[interval].ravel()
                for name in worldtube.FACE_NAMES
            ]
        )
        candidate_emf = np.concatenate(
            [
                checked_candidate[name].emf_u[interval].ravel()
                for name in worldtube.FACE_NAMES
            ]
            + [
                checked_candidate[name].emf_v[interval].ravel()
                for name in worldtube.FACE_NAMES
            ]
        )
        interval_errors.append(
            {
                "interval": interval,
                "time_range": [
                    float(candidate_times_array[interval]),
                    float(candidate_times_array[interval + 1]),
                ],
                "dt": float(dt),
                "emf": _error_norms(reference_emf, candidate_emf),
            }
        )

    reference_emf = np.concatenate(
        (_concatenate(coarsened, "emf_u"), _concatenate(coarsened, "emf_v"))
    )
    candidate_emf = np.concatenate(
        (
            _concatenate(checked_candidate, "emf_u"),
            _concatenate(checked_candidate, "emf_v"),
        )
    )
    return {
        "classification": CLASSIFICATION,
        "reference_sample_count": int(reference_times_array.size),
        "candidate_sample_count": int(candidate_times_array.size),
        "reference_endpoint_indices": indices.tolist(),
        "time_range": [
            float(candidate_times_array[0]),
            float(candidate_times_array[-1]),
        ],
        "face_cells_per_axis": int(reference_shape[-1]),
        "state": _error_norms(
            _concatenate(coarsened, "cell_state"),
            _concatenate(checked_candidate, "cell_state"),
        ),
        "state_by_variable": state_by_variable,
        "normal_flux": _error_norms(
            _concatenate(coarsened, "normal_flux"),
            _concatenate(checked_candidate, "normal_flux"),
        ),
        "emf": _error_norms(reference_emf, candidate_emf),
        "intervals": interval_errors,
        "candidate_validation": worldtube.validate_worldtube(
            candidate_times_array, checked_candidate
        ),
    }


def read_worldtube_any(
    path: Path,
) -> tuple[np.ndarray, dict[str, worldtube.FaceData], dict[str, object]]:
    candidate = path.expanduser().resolve(strict=True)
    if candidate.suffix == ".json":
        document = json.loads(candidate.read_text(encoding="utf-8"))
        if document.get("classification") == worldtube.OUTER_STREAM_CLASSIFICATION:
            return worldtube.read_outer_stream(candidate)
    return worldtube.read_worldtube(candidate)


def analyze_files(reference: Path, candidate: Path) -> dict[str, object]:
    reference_times, reference_faces, reference_metadata = read_worldtube_any(reference)
    candidate_times, candidate_faces, candidate_metadata = read_worldtube_any(candidate)
    state_variables = candidate_metadata.get(
        "state_variables", reference_metadata.get("state_variables", [])
    )
    report = compare_worldtubes(
        reference_times,
        reference_faces,
        candidate_times,
        candidate_faces,
        state_variables=state_variables,
    )
    report["reference"] = str(reference.expanduser().resolve())
    report["candidate"] = str(candidate.expanduser().resolve())
    report["candidate_sampling_diagnostics"] = candidate_metadata.get(
        "sampling_diagnostics"
    )
    return report


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    report = analyze_files(arguments.reference, arguments.candidate)
    encoded = json.dumps(report, indent=2, sort_keys=True)
    if arguments.output is None:
        print(encoded)
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(encoded + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
