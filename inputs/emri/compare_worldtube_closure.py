#!/usr/bin/env python3
"""Compare a replayed EMRI inner cube with the same subvolume of an outer run."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
VIS_PYTHON = ROOT / "vis" / "python"
if str(VIS_PYTHON) not in sys.path:
    sys.path.insert(0, str(VIS_PYTHON))

import bin_convert  # noqa: E402


DEFAULT_VARIABLES = (
    "dens",
    "velx",
    "vely",
    "velz",
    "press",
    "bcc1",
    "bcc2",
    "bcc3",
)
VECTOR_GROUPS = {
    "velocity": ("velx", "vely", "velz"),
    "magnetic": ("bcc1", "bcc2", "bcc3"),
}


def _axis_spacing(data: dict[str, Any]) -> tuple[float, float, float]:
    result = (
        (float(data["x1max"]) - float(data["x1min"])) / int(data["Nx1"]),
        (float(data["x2max"]) - float(data["x2min"])) / int(data["Nx2"]),
        (float(data["x3max"]) - float(data["x3min"])) / int(data["Nx3"]),
    )
    if not all(math.isfinite(value) and value > 0.0 for value in result):
        raise ValueError("binary output has invalid uniform-grid spacing")
    return result


def assemble_uniform_grid(data: dict[str, Any], variable: str) -> np.ndarray:
    """Assemble fixed-level MeshBlock arrays into one ``(z,y,x)`` array."""

    if variable not in data["mb_data"]:
        raise ValueError(f"binary output does not contain {variable}")
    shape = (int(data["Nx3"]), int(data["Nx2"]), int(data["Nx1"]))
    result = np.full(shape, np.nan, dtype=np.float64)
    coverage = np.zeros(shape, dtype=np.uint8)
    spacing = _axis_spacing(data)
    global_minimum = (
        float(data["x1min"]),
        float(data["x2min"]),
        float(data["x3min"]),
    )
    geometries = np.asarray(data["mb_geometry"], dtype=np.float64)
    arrays = data["mb_data"][variable]
    if len(arrays) != geometries.shape[0]:
        raise ValueError("MeshBlock geometry and data counts differ")
    for geometry, values in zip(geometries, arrays, strict=True):
        block = np.asarray(values, dtype=np.float64)
        if block.ndim != 3:
            raise ValueError("closure comparison requires three-dimensional outputs")
        starts = []
        for axis in range(3):
            lower = float(geometry[2 * axis])
            upper = float(geometry[2 * axis + 1])
            cell_count = block.shape[2 - axis]
            block_spacing = (upper - lower) / cell_count
            if not math.isclose(
                block_spacing, spacing[axis], rel_tol=2.0e-12, abs_tol=2.0e-14
            ):
                raise ValueError("closure comparison does not support AMR output")
            index = round((lower - global_minimum[axis]) / spacing[axis])
            if not math.isclose(
                lower,
                global_minimum[axis] + index * spacing[axis],
                rel_tol=2.0e-12,
                abs_tol=2.0e-14,
            ):
                raise ValueError("MeshBlock is not aligned with the declared root grid")
            starts.append(index)
        slices = tuple(
            slice(starts[axis], starts[axis] + block.shape[2 - axis])
            for axis in (2, 1, 0)
        )
        if np.any(coverage[slices] != 0):
            raise ValueError("uniform output MeshBlocks overlap")
        result[slices] = block
        coverage[slices] = 1
    if np.any(coverage != 1) or not np.isfinite(result).all():
        raise ValueError("uniform output does not cover its declared domain exactly")
    return result


def assemble_physical_grid(
    data: dict[str, Any],
    variable: str,
    adiabatic_index: float = 4.0 / 3.0,
) -> np.ndarray:
    """Assemble a canonical primitive, converting ``eint`` to pressure if needed."""

    if variable != "press" or "press" in data["mb_data"]:
        return assemble_uniform_grid(data, variable)
    if "eint" not in data["mb_data"]:
        raise ValueError("binary output contains neither press nor eint")
    gamma = float(adiabatic_index)
    if not math.isfinite(gamma) or gamma <= 1.0:
        raise ValueError("adiabatic_index must be finite and greater than one")
    return (gamma - 1.0) * assemble_uniform_grid(data, "eint")


def reference_subvolume(
    reference: dict[str, Any],
    candidate: dict[str, Any],
    variable: str,
    adiabatic_index: float = 4.0 / 3.0,
) -> np.ndarray:
    """Extract the candidate-aligned subvolume from a uniform reference output."""
    reference_spacing = _axis_spacing(reference)
    candidate_spacing = _axis_spacing(candidate)
    starts = []
    for axis, label in enumerate(("x1", "x2", "x3")):
        if not math.isclose(
            reference_spacing[axis],
            candidate_spacing[axis],
            rel_tol=2.0e-12,
            abs_tol=2.0e-14,
        ):
            raise ValueError("reference and candidate cell spacings differ")
        lower = float(candidate[f"{label}min"])
        reference_lower = float(reference[f"{label}min"])
        index = round((lower - reference_lower) / reference_spacing[axis])
        if not math.isclose(
            lower,
            reference_lower + index * reference_spacing[axis],
            rel_tol=2.0e-12,
            abs_tol=2.0e-14,
        ):
            raise ValueError("candidate domain is not aligned inside the reference")
        starts.append(index)
    array = assemble_physical_grid(reference, variable, adiabatic_index)
    counts = (int(candidate["Nx1"]), int(candidate["Nx2"]), int(candidate["Nx3"]))
    slices = tuple(
        slice(starts[axis], starts[axis] + counts[axis]) for axis in (2, 1, 0)
    )
    result = array[slices]
    expected = (counts[2], counts[1], counts[0])
    if result.shape != expected:
        raise ValueError("candidate domain extends beyond the reference output")
    return result


def _error_norms(reference: np.ndarray, candidate: np.ndarray) -> dict[str, float | None]:
    difference = candidate - reference
    reference_l1 = float(np.sum(np.abs(reference)))
    reference_l2 = float(np.linalg.norm(reference.ravel()))
    reference_linf = float(np.max(np.abs(reference)))
    return {
        "relative_l1": (
            float(np.sum(np.abs(difference))) / reference_l1
            if reference_l1 > 0.0
            else None
        ),
        "relative_l2": (
            float(np.linalg.norm(difference.ravel())) / reference_l2
            if reference_l2 > 0.0
            else None
        ),
        "absolute_linf": float(np.max(np.abs(difference))),
        "relative_linf": (
            float(np.max(np.abs(difference))) / reference_linf
            if reference_linf > 0.0
            else None
        ),
    }


def compare_loaded_outputs(
    reference: dict[str, Any],
    candidate: dict[str, Any],
    variables: tuple[str, ...] = DEFAULT_VARIABLES,
    adiabatic_index: float = 4.0 / 3.0,
) -> dict[str, Any]:
    """Compare a candidate cube with its coincident reference subvolume."""

    time_scale = max(abs(float(reference["time"])), abs(float(candidate["time"])), 1.0)
    if not math.isclose(
        float(reference["time"]),
        float(candidate["time"]),
        rel_tol=2.0e-12,
        abs_tol=2.0e-14 * time_scale,
    ):
        raise ValueError("reference and candidate output times differ")
    reference_arrays = {
        name: reference_subvolume(
            reference, candidate, name, adiabatic_index
        )
        for name in variables
    }
    candidate_arrays = {
        name: assemble_physical_grid(candidate, name, adiabatic_index)
        for name in variables
    }
    variable_errors = {
        name: _error_norms(reference_arrays[name], candidate_arrays[name])
        for name in variables
    }
    group_errors = {}
    for group, names in VECTOR_GROUPS.items():
        if not set(names).issubset(variables):
            continue
        reference_vector = np.stack([reference_arrays[name] for name in names])
        candidate_vector = np.stack([candidate_arrays[name] for name in names])
        group_errors[group] = _error_norms(reference_vector, candidate_vector)
    return {
        "classification": "athenak-emri-worldtube-closure-v1",
        "time": float(reference["time"]),
        "reference_cycle": int(reference["cycle"]),
        "candidate_cycle": int(candidate["cycle"]),
        "candidate_shape": [
            int(candidate["Nx3"]),
            int(candidate["Nx2"]),
            int(candidate["Nx1"]),
        ],
        "variables": variable_errors,
        "vector_groups": group_errors,
    }


def compare_files(
    reference: Path,
    candidate: Path,
    adiabatic_index: float = 4.0 / 3.0,
) -> dict[str, Any]:
    return compare_loaded_outputs(
        bin_convert.read_binary(str(reference)),
        bin_convert.read_binary(str(candidate)),
        adiabatic_index=adiabatic_index,
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--gamma", type=float, default=4.0 / 3.0)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    report = compare_files(
        arguments.reference, arguments.candidate, arguments.gamma
    )
    encoded = json.dumps(report, indent=2, sort_keys=True)
    if arguments.output is None:
        print(encoded)
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(encoded + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
