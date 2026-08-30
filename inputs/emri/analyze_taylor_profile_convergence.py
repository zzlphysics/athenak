#!/usr/bin/env python3
"""Assess fitting-scale and snapshot-cadence sensitivity of Taylor winds."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import build_taylor_worldtube_series as series
import extract_static_taylor_worldtube as static


CLASSIFICATION = "athenak-emri-real-taylor-convergence-v1"
GROUPS = {
    "density": ("rho0",),
    "pressure": ("pgas0",),
    "velocity": ("u1", "u2", "u3"),
    "magnetic_field": ("b1", "b2", "b3"),
    "density_gradient": tuple(f"dlnrho_dxh{i}" for i in range(1, 4)),
    "pressure_gradient": tuple(f"dlnpgas_dxh{i}" for i in range(1, 4)),
    "velocity_gradient": tuple(
        f"du{i}_dxh{j}" for i in range(1, 4) for j in range(1, 4)
    ),
    "magnetic_gradient": tuple(
        f"db{i}_dxh{j}" for i in range(1, 4) for j in range(1, 4)
    ),
}


def _group_indices() -> dict[str, np.ndarray]:
    order = list(static.PROFILE_PARAMETER_ORDER)
    return {
        group: np.asarray([order.index(name) for name in names], dtype=np.int64)
        for group, names in GROUPS.items()
    }


def group_error(
    reference: np.ndarray, candidate: np.ndarray, indices: np.ndarray
) -> dict[str, float | None]:
    ref = np.asarray(reference, dtype=np.float64)[:, indices]
    cand = np.asarray(candidate, dtype=np.float64)[:, indices]
    if ref.shape != cand.shape or ref.ndim != 2 or not (
        np.isfinite(ref).all() and np.isfinite(cand).all()
    ):
        raise ValueError("profile comparison arrays are incompatible")
    difference = cand - ref
    absolute_rms = math.sqrt(float(np.mean(np.sum(difference**2, axis=1))))
    reference_rms = math.sqrt(float(np.mean(np.sum(ref**2, axis=1))))
    candidate_rms = math.sqrt(float(np.mean(np.sum(cand**2, axis=1))))
    scale = max(reference_rms, np.finfo(float).tiny)
    symmetric_scale = max(
        0.5 * (reference_rms + candidate_rms), np.finfo(float).tiny
    )
    return {
        "absolute_rms": absolute_rms,
        "reference_rms": reference_rms,
        "relative_l2": absolute_rms / scale,
        "symmetric_relative_l2": absolute_rms / symmetric_scale,
        "maximum_absolute": float(np.max(np.abs(difference))),
    }


def selected_indices(count: int, stride: int) -> np.ndarray:
    if count < 2 or stride < 1:
        raise ValueError("cadence selection requires at least two samples and stride >= 1")
    result = list(range(0, count, stride))
    if result[-1] != count - 1:
        result.append(count - 1)
    return np.asarray(result, dtype=np.int64)


def cadence_candidate(
    times: np.ndarray, values: np.ndarray, stride: int
) -> tuple[np.ndarray, np.ndarray]:
    selected = selected_indices(times.size, stride)
    candidate = np.column_stack(
        [np.interp(times, times[selected], values[selected, field])
         for field in range(values.shape[1])]
    )
    return selected, candidate


def _condition_maximum(observed: float, threshold: float) -> dict[str, object]:
    return {
        "relation": "maximum",
        "observed": observed,
        "threshold": threshold,
        "passed": math.isfinite(observed) and observed <= threshold,
    }


def _condition_minimum(observed: float, threshold: float) -> dict[str, object]:
    return {
        "relation": "minimum",
        "observed": observed,
        "threshold": threshold,
        "passed": math.isfinite(observed) and observed >= threshold,
    }


def _load(path: Path) -> tuple[dict[str, object], np.ndarray, np.ndarray]:
    resolved = path.expanduser().resolve(strict=True)
    document = json.loads(resolved.read_text(encoding="utf-8"))
    if document.get("classification") != series.TABLE_CLASSIFICATION:
        raise ValueError(f"{resolved} is not a Taylor-series manifest")
    times, values = series.read_table(
        Path(document["table_file"]).expanduser().resolve(strict=True)
    )
    if list(document.get("parameter_order", ())) != list(
        static.PROFILE_PARAMETER_ORDER
    ):
        raise ValueError(f"{resolved} has an incompatible parameter order")
    document["_path"] = str(resolved)
    return document, times, values


def analyze(arguments: argparse.Namespace) -> dict[str, object]:
    loaded = [_load(path) for path in arguments.series]
    loaded.sort(key=lambda entry: float(entry[0]["fit_radius_source"]))
    radii = [float(entry[0]["fit_radius_source"]) for entry in loaded]
    if len(set(radii)) != len(radii):
        raise ValueError("fit radii must be unique")
    reference_document, reference_times, reference_values = loaded[0]
    invariant_keys = (
        "source_manifest_sha256",
        "global_length_in_local_units",
        "density_renormalization",
        "metric_source",
        "parameter_order",
    )
    for document, times, _ in loaded[1:]:
        if any(document.get(key) != reference_document.get(key) for key in invariant_keys):
            raise ValueError("Taylor-series manifests do not share one source contract")
        if not np.array_equal(times, reference_times):
            raise ValueError("Taylor-series manifests do not share exact sample times")
    indices = _group_indices()
    fitting_scale = {}
    for document, _, values in loaded:
        label = f"r{float(document['fit_radius_source']):g}"
        fitting_scale[label] = {
            "fit_radius_source": float(document["fit_radius_source"]),
            "reference": document is reference_document,
            "groups": {
                group: group_error(reference_values, values, selected)
                for group, selected in indices.items()
            },
            "maximum_fit_residuals": {
                name: max(float(sample["fit_diagnostics"][name])
                          for sample in document["samples"])
                for name in (
                    "density_log_weighted_rms",
                    "pressure_log_weighted_rms",
                    "velocity_relative_rms",
                    "magnetic_relative_rms",
                )
            },
            "minimum_fit_samples": min(
                int(sample["fit_diagnostics"]["sample_count"])
                for sample in document["samples"]
            ),
            "maximum_magnetic_gradient_trace": max(
                abs(float(sample["fit_diagnostics"]["magnetic_gradient_trace"]))
                for sample in document["samples"]
            ),
            "maximum_cell_volume_ratio": max(
                float(sample["fit_diagnostics"]["maximum_cell_volume"])
                / float(sample["fit_diagnostics"]["minimum_cell_volume"])
                for sample in document["samples"]
            ),
        }
    cadence = {}
    for stride in arguments.cadence_strides:
        selected, candidate = cadence_candidate(
            reference_times, reference_values, stride
        )
        cadence[f"s{stride}"] = {
            "stride": stride,
            "reference": stride == 1,
            "selected_indices": selected.tolist(),
            "selected_times": reference_times[selected].tolist(),
            "groups": {
                group: group_error(reference_values, candidate, columns)
                for group, columns in indices.items()
            },
        }
    reference_fit = fitting_scale[f"r{radii[0]:g}"]
    nonreference_spatial = [
        entry for entry in fitting_scale.values() if not entry["reference"]
    ]
    nonreference_cadence = [
        entry for entry in cadence.values() if not entry["reference"]
    ]
    maximum_spatial = max(
        float(group["symmetric_relative_l2"])
        for entry in nonreference_spatial
        for group in entry["groups"].values()
    )
    maximum_cadence = max(
        float(group["symmetric_relative_l2"])
        for entry in nonreference_cadence
        for group in entry["groups"].values()
    )
    residuals = reference_fit["maximum_fit_residuals"]
    conditions = {
        "minimum_fit_samples": _condition_minimum(
            float(reference_fit["minimum_fit_samples"]),
            float(arguments.minimum_fit_samples),
        ),
        "density_log_fit_rms": _condition_maximum(
            float(residuals["density_log_weighted_rms"]),
            arguments.maximum_density_log_fit_rms,
        ),
        "pressure_log_fit_rms": _condition_maximum(
            float(residuals["pressure_log_weighted_rms"]),
            arguments.maximum_pressure_log_fit_rms,
        ),
        "velocity_fit_relative_rms": _condition_maximum(
            float(residuals["velocity_relative_rms"]),
            arguments.maximum_velocity_fit_relative_rms,
        ),
        "magnetic_fit_relative_rms": _condition_maximum(
            float(residuals["magnetic_relative_rms"]),
            arguments.maximum_magnetic_fit_relative_rms,
        ),
        "magnetic_gradient_trace": _condition_maximum(
            float(reference_fit["maximum_magnetic_gradient_trace"]),
            arguments.maximum_magnetic_gradient_trace,
        ),
        "cell_volume_ratio": _condition_maximum(
            float(reference_fit["maximum_cell_volume_ratio"]),
            arguments.maximum_cell_volume_ratio,
        ),
        "fitting_scale_sensitivity": _condition_maximum(
            maximum_spatial, arguments.maximum_fitting_scale_sensitivity
        ),
        "cadence_sensitivity": _condition_maximum(
            maximum_cadence, arguments.maximum_cadence_sensitivity
        ),
    }
    passed = all(bool(condition["passed"]) for condition in conditions.values())
    return {
        "classification": CLASSIFICATION,
        "run_status": "completed",
        "assessment": {
            "passed": passed,
            "science_ready": False,
            "conditions": conditions,
            "science_blockers": [
                (
                    "fit-radius variants share one source grid and do not establish "
                    "independent source-resolution convergence"
                ),
                (
                    "local times between global dumps are far longer than a short "
                    "EMRI wind-tunnel relaxation and should be treated as a "
                    "frozen-snapshot ensemble"
                ),
            ],
        },
        "reference_fit_radius_source": radii[0],
        "fit_radii_source": radii,
        "cadence_strides": arguments.cadence_strides,
        "sample_times_local": reference_times.tolist(),
        "fitting_scale_sensitivity": fitting_scale,
        "cadence_sensitivity": cadence,
        "source": {
            key: reference_document[key]
            for key in invariant_keys
        },
        "input_manifests": [
            {
                "path": document["_path"],
                "sha256": static.file_sha256(Path(document["_path"])),
            }
            for document, _, _ in loaded
        ],
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--series", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cadence-strides", type=int, nargs="+", default=(1, 2))
    parser.add_argument("--minimum-fit-samples", type=int, default=16)
    parser.add_argument("--maximum-density-log-fit-rms", type=float, default=0.03)
    parser.add_argument("--maximum-pressure-log-fit-rms", type=float, default=0.03)
    parser.add_argument("--maximum-velocity-fit-relative-rms", type=float, default=0.03)
    parser.add_argument("--maximum-magnetic-fit-relative-rms", type=float, default=0.05)
    parser.add_argument("--maximum-magnetic-gradient-trace", type=float, default=1e-12)
    parser.add_argument("--maximum-cell-volume-ratio", type=float, default=1.000001)
    parser.add_argument("--maximum-fitting-scale-sensitivity", type=float, default=0.15)
    parser.add_argument("--maximum-cadence-sensitivity", type=float, default=0.1)
    parser.add_argument("--fail-on-gate", action="store_true")
    arguments = parser.parse_args()
    if len(arguments.series) < 2:
        parser.error("at least two fitting radii are required")
    if (
        not arguments.cadence_strides
        or 1 not in arguments.cadence_strides
        or any(value < 1 for value in arguments.cadence_strides)
        or len(set(arguments.cadence_strides)) != len(arguments.cadence_strides)
    ):
        parser.error("cadence strides must be unique, positive, and include one")
    arguments.cadence_strides = sorted(arguments.cadence_strides)
    if arguments.minimum_fit_samples < 4:
        parser.error("minimum fit samples must be at least four")
    for name, value in vars(arguments).items():
        if name.startswith("maximum_") and (
            not math.isfinite(value) or value <= 0.0
        ):
            parser.error(f"--{name.replace('_', '-')} must be finite and positive")
    return arguments


def main() -> int:
    arguments = parse_arguments()
    report = analyze(arguments)
    output = arguments.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    state = "PASS" if report["assessment"]["passed"] else "FAIL"
    print(f"Taylor_convergence={state}")
    print(f"science_ready={report['assessment']['science_ready']}")
    print(output)
    return 1 if arguments.fail_on_gate and state == "FAIL" else 0


if __name__ == "__main__":
    raise SystemExit(main())
