#!/usr/bin/env python3
"""Measure reflected RMHD characteristic content in a worldtube closure run.

The outer solution restricted to the inner cube is the no-interface reference.  The
difference between the replayed inner solution and that reference is projected onto a
single primitive-variable RMHD eigenbasis.  Modes whose x1 characteristic speed has
the opposite sign to the injected wave are reported as reflected content.

This is an amplitude-like closure coefficient in a Euclidean-normalized primitive
eigenbasis.  It is not an energy reflection coefficient; an energy statement would
require a physical symmetrizer/energy norm.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import compare_worldtube_closure as closure
import worldtube_characteristics as characteristics


MODE_NAMES = (
    "left_fast",
    "left_alfven",
    "left_slow",
    "entropy",
    "right_slow",
    "right_alfven",
    "right_fast",
)
PRIMITIVE_OUTPUTS = ("dens", "press", "velx", "vely", "velz", "bcc2", "bcc3")


def _rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values, dtype=np.float64))))


def _combined_rms(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum(np.square(values, dtype=np.float64), axis=0))))


def primitive_arrays(
    output: dict[str, Any],
    candidate: dict[str, Any] | None = None,
    adiabatic_index: float = 4.0 / 3.0,
) -> np.ndarray:
    """Return ``rho,p,u1,u2,u3,B2,B3`` arrays, optionally on a subvolume."""

    if candidate is None:
        arrays = [
            closure.assemble_physical_grid(output, name, adiabatic_index)
            for name in PRIMITIVE_OUTPUTS
        ]
    else:
        arrays = [
            closure.reference_subvolume(
                output, candidate, name, adiabatic_index
            )
            for name in PRIMITIVE_OUTPUTS
        ]
    return np.stack(arrays)


def background_from_reference(
    reference: dict[str, Any], adiabatic_index: float = 4.0 / 3.0
) -> tuple[np.ndarray, float]:
    """Infer the constant state by averaging the full-wavelength outer domain."""

    primitive = primitive_arrays(reference, adiabatic_index=adiabatic_index)
    background = np.mean(primitive, axis=(1, 2, 3))
    normal_magnetic = float(np.mean(closure.assemble_uniform_grid(reference, "bcc1")))
    if not np.isfinite(background).all() or not math.isfinite(normal_magnetic):
        raise ValueError("reference mean state is not finite")
    return background, normal_magnetic


def measure_mode_content(
    reference_primitive: object,
    candidate_primitive: object,
    background_primitive: object,
    normal_magnetic: float,
    adiabatic_index: float,
    source_mode: int,
    speed_tolerance: float = 1.0e-8,
) -> dict[str, Any]:
    """Project closure error and return opposite-speed amplitude coefficients."""

    reference = np.asarray(reference_primitive, dtype=np.float64)
    candidate = np.asarray(candidate_primitive, dtype=np.float64)
    background = np.asarray(background_primitive, dtype=np.float64)
    if reference.shape != candidate.shape or reference.ndim < 2 or reference.shape[0] != 7:
        raise ValueError("reference and candidate primitives must share shape (7,...)")
    if background.shape != (7,) or not np.isfinite(background).all():
        raise ValueError("background_primitive must be a finite seven-vector")
    if not 0 <= source_mode < 7:
        raise ValueError("source_mode must be in [0,6]")
    if not math.isfinite(speed_tolerance) or speed_tolerance < 0.0:
        raise ValueError("speed_tolerance must be finite and nonnegative")
    if not np.isfinite(reference).all() or not np.isfinite(candidate).all():
        raise ValueError("primitive arrays must be finite")

    basis = characteristics.characteristic_basis(
        background, normal_magnetic, adiabatic_index
    )
    points = int(np.prod(reference.shape[1:]))
    reference_delta = reference.reshape(7, points) - background[:, None]
    closure_error = candidate.reshape(7, points) - reference.reshape(7, points)
    reference_amplitude = basis.left_eigenvectors @ reference_delta
    error_amplitude = basis.left_eigenvectors @ closure_error
    demeaned_error = closure_error - np.mean(closure_error, axis=1, keepdims=True)
    demeaned_amplitude = basis.left_eigenvectors @ demeaned_error

    source_rms = _rms(reference_amplitude[source_mode])
    if not math.isfinite(source_rms) or source_rms <= 0.0:
        raise ValueError("reference source-mode RMS amplitude is zero")
    source_speed = float(basis.speeds[source_mode])
    if abs(source_speed) <= speed_tolerance:
        raise ValueError("source mode is stationary at the requested speed tolerance")
    if source_speed > 0.0:
        reflected = np.flatnonzero(basis.speeds < -speed_tolerance)
    else:
        reflected = np.flatnonzero(basis.speeds > speed_tolerance)
    if reflected.size == 0:
        raise ValueError("background state has no opposite-speed characteristic modes")

    error_rms = np.asarray([_rms(error_amplitude[index]) for index in range(7)])
    demeaned_rms = np.asarray(
        [_rms(demeaned_amplitude[index]) for index in range(7)]
    )
    non_source = np.asarray([index for index in range(7) if index != source_mode])
    reflected_rms = _combined_rms(error_amplitude[reflected])
    reflected_demeaned_rms = _combined_rms(demeaned_amplitude[reflected])
    reference_leakage = _combined_rms(reference_amplitude[non_source]) / source_rms
    return {
        "classification": "athenak-emri-worldtube-reflected-modes-v1",
        "normalization": (
            "RMS primitive-eigenmode closure error divided by reference source-mode RMS"
        ),
        "energy_interpretation": (
            "amplitude-like only; do not square as an energy coefficient without a "
            "physical symmetrizer"
        ),
        "source_mode": source_mode,
        "source_mode_name": MODE_NAMES[source_mode],
        "source_speed": source_speed,
        "source_rms_amplitude": source_rms,
        "characteristic_speeds": [float(value) for value in basis.speeds],
        "reflected_mode_indices": [int(value) for value in reflected],
        "reflected_mode_names": [MODE_NAMES[int(value)] for value in reflected],
        "reflected_amplitude_coefficient": reflected_rms / source_rms,
        "demeaned_reflected_amplitude_coefficient": reflected_demeaned_rms / source_rms,
        "mode_error_coefficients": [float(value / source_rms) for value in error_rms],
        "demeaned_mode_error_coefficients": [
            float(value / source_rms) for value in demeaned_rms
        ],
        "reference_non_source_leakage": float(reference_leakage),
        "eigenvector_condition_number": basis.condition_number,
        "eigensystem_residual": basis.jacobian_residual,
    }


def analyze_loaded_outputs(
    reference: dict[str, Any],
    candidate: dict[str, Any],
    source_mode: int,
    adiabatic_index: float = 4.0 / 3.0,
    speed_tolerance: float = 1.0e-8,
) -> dict[str, Any]:
    """Analyze one outer-reference/inner-replay output pair."""

    time_scale = max(abs(float(reference["time"])), abs(float(candidate["time"])), 1.0)
    if not math.isclose(
        float(reference["time"]),
        float(candidate["time"]),
        rel_tol=2.0e-12,
        abs_tol=2.0e-14 * time_scale,
    ):
        raise ValueError("reference and candidate output times differ")
    background, normal_magnetic = background_from_reference(
        reference, adiabatic_index
    )
    reference_primitive = primitive_arrays(
        reference, candidate, adiabatic_index
    )
    candidate_primitive = primitive_arrays(
        candidate, adiabatic_index=adiabatic_index
    )
    result = measure_mode_content(
        reference_primitive,
        candidate_primitive,
        background,
        normal_magnetic,
        adiabatic_index,
        source_mode,
        speed_tolerance,
    )
    result.update(
        {
            "time": float(reference["time"]),
            "reference_cycle": int(reference["cycle"]),
            "candidate_cycle": int(candidate["cycle"]),
            "candidate_shape": [
                int(candidate["Nx3"]),
                int(candidate["Nx2"]),
                int(candidate["Nx1"]),
            ],
            "background_primitive": [float(value) for value in background],
            "background_normal_magnetic": normal_magnetic,
            "closure": closure.compare_loaded_outputs(
                reference, candidate, adiabatic_index=adiabatic_index
            ),
        }
    )
    return result


def analyze_files(
    reference: Path,
    candidate: Path,
    source_mode: int,
    adiabatic_index: float = 4.0 / 3.0,
    speed_tolerance: float = 1.0e-8,
) -> dict[str, Any]:
    return analyze_loaded_outputs(
        closure.bin_convert.read_binary(str(reference)),
        closure.bin_convert.read_binary(str(candidate)),
        source_mode,
        adiabatic_index,
        speed_tolerance,
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--source-mode", type=int, required=True, choices=range(7))
    parser.add_argument("--gamma", type=float, default=4.0 / 3.0)
    parser.add_argument("--speed-tolerance", type=float, default=1.0e-8)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    report = analyze_files(
        arguments.reference,
        arguments.candidate,
        arguments.source_mode,
        arguments.gamma,
        arguments.speed_tolerance,
    )
    encoded = json.dumps(report, indent=2, sort_keys=True)
    if arguments.output is None:
        print(encoded)
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(encoded + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
