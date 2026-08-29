#!/usr/bin/env python3
"""Audit GRMHD characteristic content and HLL risk in an EMRI worldtube file.

The input may be a packed ``.npz`` worldtube, an outer-stream manifest, or the strict
inner replay binary.  The report is deliberately local: it diagnoses whether the
boundary is super-fast or has a mixed characteristic fan, whether the numerical
eigenbasis is usable, how finely the fluid data are sampled in time, and how strongly
the two-wave HLL flux changes individual linear RMHD modes.  It does not claim to
replace a propagated wave-packet reflection measurement.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import worldtube_characteristics as characteristics
import worldtube_flux_emf as worldtube


def _read_input(
    path: Path,
) -> tuple[np.ndarray, dict[str, worldtube.FaceData], dict[str, object]]:
    if path.name.endswith(".manifest.json"):
        return worldtube.read_outer_stream(path)
    if path.suffix == ".bin":
        return worldtube.read_inner_binary(path)
    return worldtube.read_worldtube(path)


def _maximum_finite(values: list[float], default: float = 0.0) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return max(finite, default=default)


def _gain_error_summary(values: list[float]) -> dict[str, float]:
    finite = np.asarray([value for value in values if math.isfinite(value)])
    if finite.size == 0:
        return {"median": 0.0, "p95": 0.0, "maximum": 0.0}
    return {
        "median": float(np.median(finite)),
        "p95": float(np.percentile(finite, 95.0)),
        "maximum": float(np.max(finite)),
    }


def _failure_category(error: Exception) -> str:
    message = str(error)
    if "ill-conditioned" in message:
        return "ill_conditioned_eigenbasis"
    if "not numerically real" in message or "non-hyperbolic" in message:
        return "non_real_eigenbasis"
    if "singular" in message:
        return "singular_primitive_jacobian"
    if isinstance(error, ValueError):
        return "invalid_state"
    return "other_numerical_failure"


def audit_worldtube(
    times: object,
    faces: dict[str, worldtube.FaceData],
    metadata: dict[str, object],
    adiabatic_index: float,
    sample_stride: int = 1,
    speed_tolerance: float = 1.0e-8,
) -> dict[str, Any]:
    """Return a JSON-serializable local-characteristic audit."""

    checked_times = worldtube.validate_times(times)
    if sample_stride < 1:
        raise ValueError("sample_stride must be at least one")
    if not math.isfinite(speed_tolerance) or speed_tolerance < 0.0:
        raise ValueError("speed_tolerance must be finite and nonnegative")
    gamma = float(adiabatic_index)
    if not math.isfinite(gamma) or gamma <= 1.0:
        raise ValueError("adiabatic_index must be finite and greater than one")
    checked_faces = {
        name: worldtube.validate_face(faces[name], checked_times, name)
        for name in worldtube.FACE_NAMES
    }
    first_shape = checked_faces[worldtube.FACE_NAMES[0]].cell_state.shape
    cells = first_shape[-1]
    if first_shape[-2] != cells or first_shape[1] < 8:
        raise ValueError(
            "characteristic audit requires square faces and ideal-MHD state with "
            "rho,u1,u2,u3,pgas,...,bcc1,bcc2,bcc3"
        )
    state_variables = metadata.get("state_variables", [])
    if not state_variables and isinstance(metadata.get("source_metadata"), dict):
        state_variables = metadata["source_metadata"].get("state_variables", [])
    thermodynamic = "pgas"
    if isinstance(state_variables, list) and len(state_variables) == first_shape[1]:
        thermodynamic = str(state_variables[4])
    if thermodynamic not in ("pgas", "eint"):
        raise ValueError(
            "characteristic audit requires state variable 5 to be pgas or eint"
        )
    half_width = float(metadata.get("half_width", math.nan))
    if not math.isfinite(half_width) or half_width <= 0.0:
        outer = metadata.get("outer_stream_manifest")
        if isinstance(outer, dict):
            half_width = float(outer.get("half_width", math.nan))
    if not math.isfinite(half_width) or half_width <= 0.0:
        raise ValueError("worldtube metadata must provide a positive half_width")
    dx = 2.0 * half_width / cells
    area = dx * dx
    maximum_speed = 0.0
    face_reports: dict[str, object] = {}
    total_samples = 0
    total_failures = 0
    total_mixed = 0
    all_gain_errors: list[float] = []

    for face_name in worldtube.FACE_NAMES:
        face = checked_faces[face_name]
        incoming_histogram: Counter[int] = Counter()
        conditions: list[float] = []
        residuals: list[float] = []
        gain_errors: list[float] = []
        speeds: list[float] = []
        failure_reasons: Counter[str] = Counter()
        face_failures = 0
        face_samples = 0
        for time_index in range(0, checked_times.size, sample_stride):
            for v in range(0, cells, sample_stride):
                for u in range(0, cells, sample_stride):
                    face_samples += 1
                    state_raw = face.cell_state[time_index, :, v, u]
                    state = np.concatenate((state_raw[:5], state_raw[-3:]))
                    if thermodynamic == "eint":
                        state[4] *= gamma - 1.0
                    normal_field = float(
                        face.normal_flux[time_index, v, u] / area
                    )
                    try:
                        primitive = characteristics.state_to_face_primitive(
                            state, face_name, normal_field
                        )
                        basis = characteristics.characteristic_basis(
                            primitive, normal_field, gamma
                        )
                        hll = characteristics.linear_hll_gains_from_speeds(
                            basis.speeds,
                            speed_tolerance=speed_tolerance,
                        )
                    except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
                        face_failures += 1
                        failure_reasons[_failure_category(error)] += 1
                        continue
                    incoming_histogram[
                        int(np.count_nonzero(basis.speeds < -speed_tolerance))
                    ] += 1
                    conditions.append(basis.condition_number)
                    residuals.append(basis.jacobian_residual)
                    speeds.extend(float(value) for value in basis.speeds)
                    gain_errors.extend(
                        float(value)
                        for value in hll.gain_error
                        if math.isfinite(float(value))
                    )
                    maximum_speed = max(
                        maximum_speed, float(np.max(np.abs(basis.speeds)))
                    )
        successful = face_samples - face_failures
        mixed = sum(
            count
            for incoming, count in incoming_histogram.items()
            if 0 < incoming < 7
        )
        face_reports[face_name] = {
            "sample_count": face_samples,
            "eigensystem_failures": face_failures,
            "eigensystem_failure_reasons": dict(sorted(failure_reasons.items())),
            "incoming_mode_histogram": {
                str(key): incoming_histogram[key] for key in sorted(incoming_histogram)
            },
            "mixed_fan_fraction": (mixed / successful) if successful else None,
            "minimum_characteristic_speed": min(speeds, default=None),
            "maximum_characteristic_speed": max(speeds, default=None),
            "maximum_eigenvector_condition_number": _maximum_finite(conditions),
            "maximum_eigensystem_residual": _maximum_finite(residuals),
            "linear_hll_flux_gain_error": _gain_error_summary(gain_errors),
        }
        total_samples += face_samples
        total_failures += face_failures
        total_mixed += mixed
        all_gain_errors.extend(gain_errors)

    successful_samples = total_samples - total_failures
    maximum_dt = float(np.max(np.diff(checked_times)))
    return {
        "classification": "athenak-emri-worldtube-characteristic-audit-v1",
        "adiabatic_index": gamma,
        "speed_tolerance": speed_tolerance,
        "sample_stride": sample_stride,
        "input_thermodynamic_variable": thermodynamic,
        "sample_count": total_samples,
        "eigensystem_failures": total_failures,
        "mixed_fan_fraction": (
            total_mixed / successful_samples if successful_samples else None
        ),
        "maximum_worldtube_cadence_courant": maximum_speed * maximum_dt / dx,
        "linear_hll_flux_gain_error": _gain_error_summary(all_gain_errors),
        "interpretation": {
            "cadence_courant": "max(|lambda| dt_worldtube / dx_face)",
            "hll_gain_error": (
                "local linear flux-response mismatch; it is a reflection-risk proxy, "
                "not a propagated wave-packet reflection coefficient"
            ),
        },
        "faces": face_reports,
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--gamma", type=float, default=4.0 / 3.0)
    parser.add_argument("--sample-stride", type=int, default=1)
    parser.add_argument("--speed-tolerance", type=float, default=1.0e-8)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    times, faces, metadata = _read_input(arguments.input)
    report = audit_worldtube(
        times,
        faces,
        metadata,
        arguments.gamma,
        sample_stride=arguments.sample_stride,
        speed_tolerance=arguments.speed_tolerance,
    )
    encoded = json.dumps(report, indent=2, sort_keys=True)
    if arguments.output is None:
        print(encoded)
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(encoded + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
