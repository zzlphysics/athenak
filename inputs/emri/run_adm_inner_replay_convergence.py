#!/usr/bin/env python3
"""Run numerical-ADM spatial/cadence convergence for EMRI force histories."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np

import analyze_force_history as force_history
import extract_global_worldtube as extract
import run_adm_inner_replay_pilot as pilot
import worldtube_flux_emf as worldtube


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "inputs" / "emri" / "emri_windtunnel_smoke.athinput"
CLASSIFICATION = "athenak-emri-adm-inner-replay-convergence-v1"
HISTORY_GROUPS = {
    "mdot": ("mdot",),
    "momentum_flux": ("Fmom_x", "Fmom_y", "Fmom_z"),
    "newtonian_volume_force": ("Fnewt_x", "Fnewt_y", "Fnewt_z"),
    "relativistic_volume_force_r1": ("Frel1_x", "Frel1_y", "Frel1_z"),
    "relativistic_volume_force_r2": ("Frel2_x", "Frel2_y", "Frel2_z"),
    "relativistic_volume_force_r3": ("Frel3_x", "Frel3_y", "Frel3_z"),
}
ADM_GROUPS = {
    "metric": (
        "adm_alpha",
        "adm_betax",
        "adm_betay",
        "adm_betaz",
        "adm_gxx",
        "adm_gxy",
        "adm_gxz",
        "adm_gyy",
        "adm_gyz",
        "adm_gzz",
    ),
    "extrinsic_curvature": (
        "adm_Kxx",
        "adm_Kxy",
        "adm_Kxz",
        "adm_Kyy",
        "adm_Kyz",
        "adm_Kzz",
    ),
}


@dataclass(frozen=True)
class ConvergenceCase:
    resolution_factor: int
    cadence_stride: int

    @property
    def name(self) -> str:
        return f"mr{self.resolution_factor}_s{self.cadence_stride}"


def convergence_cases(
    resolution_factors: list[int], cadence_strides: list[int]
) -> tuple[ConvergenceCase, ...]:
    cases = tuple(
        ConvergenceCase(factor, stride)
        for factor in resolution_factors
        for stride in cadence_strides
    )
    if len({case.name for case in cases}) != len(cases):
        raise ValueError("convergence case names are not unique")
    return cases


def reference_case(cases: tuple[ConvergenceCase, ...]) -> ConvergenceCase:
    if not cases:
        raise ValueError("convergence matrix is empty")
    target = (
        max(case.resolution_factor for case in cases),
        min(case.cadence_stride for case in cases),
    )
    matches = [
        case
        for case in cases
        if (case.resolution_factor, case.cadence_stride) == target
    ]
    if len(matches) != 1:
        raise RuntimeError("convergence matrix has no unique finest reference")
    return matches[0]


def _write_json(path: Path, document: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(document, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _interpolated_value(
    times: np.ndarray, values: np.ndarray, target: float
) -> np.ndarray:
    if target <= times[0]:
        return values[0]
    if target >= times[-1]:
        return values[-1]
    right = int(np.searchsorted(times, target, side="right"))
    left = right - 1
    fraction = (target - times[left]) / (times[right] - times[left])
    return (1.0 - fraction) * values[left] + fraction * values[right]


def compare_series(
    reference_times: object,
    reference_values: object,
    candidate_times: object,
    candidate_values: object,
) -> dict[str, object]:
    ref_times = np.asarray(reference_times, dtype=np.float64)
    cand_times = np.asarray(candidate_times, dtype=np.float64)
    ref_values = np.asarray(reference_values, dtype=np.float64)
    cand_values = np.asarray(candidate_values, dtype=np.float64)
    if (
        ref_times.ndim != 1
        or cand_times.ndim != 1
        or ref_times.size < 2
        or cand_times.size < 2
        or ref_values.shape[0] != ref_times.size
        or cand_values.shape[0] != cand_times.size
        or ref_values.shape[1:] != cand_values.shape[1:]
    ):
        raise ValueError("sampled series dimensions are incompatible")
    start = max(float(ref_times[0]), float(cand_times[0]))
    stop = min(float(ref_times[-1]), float(cand_times[-1]))
    if not stop > start:
        raise ValueError("sampled series have no positive-duration overlap")
    knots = {start, stop}
    knots.update(float(value) for value in ref_times if start < value < stop)
    knots.update(float(value) for value in cand_times if start < value < stop)
    times = np.asarray(sorted(knots), dtype=np.float64)
    difference_power = np.empty(times.size)
    reference_power = np.empty(times.size)
    maximum_absolute = 0.0
    for index, sample_time in enumerate(times):
        reference = _interpolated_value(ref_times, ref_values, float(sample_time))
        candidate = _interpolated_value(
            cand_times, cand_values, float(sample_time)
        )
        difference = candidate - reference
        difference_power[index] = float(np.sum(difference * difference))
        reference_power[index] = float(np.sum(reference * reference))
        maximum_absolute = max(
            maximum_absolute, float(np.max(np.abs(difference)))
        )
    widths = np.diff(times)
    difference_integral = float(
        np.sum(0.5 * widths * (difference_power[:-1] + difference_power[1:]))
    )
    reference_integral = float(
        np.sum(0.5 * widths * (reference_power[:-1] + reference_power[1:]))
    )
    duration = stop - start
    absolute_rms = math.sqrt(max(difference_integral / duration, 0.0))
    reference_rms = math.sqrt(max(reference_integral / duration, 0.0))
    if reference_rms > np.finfo(float).tiny:
        relative_l2: float | None = absolute_rms / reference_rms
    else:
        relative_l2 = 0.0 if absolute_rms == 0.0 else None
    return {
        "time_range": [start, stop],
        "comparison_knots": int(times.size),
        "absolute_rms": absolute_rms,
        "reference_rms": reference_rms,
        "relative_l2": relative_l2,
        "maximum_absolute": maximum_absolute,
    }


def _history_groups(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    history = force_history.ForceHistory(path.stem, path)
    if history.force_is_source_tetrad:
        raise ValueError("numerical ADM convergence requires coordinate force history")
    times = np.asarray(history.times, dtype=np.float64)
    result = {}
    for group, quantities in HISTORY_GROUPS.items():
        values = np.column_stack(
            [np.asarray(history.values[quantity], dtype=np.float64)
             for quantity in quantities]
        )
        result[group] = (times, values)
    return result


def _adm_groups(
    report: dict[str, object]
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    comparison = report.get("adm_replay_comparison")
    if not isinstance(comparison, dict):
        raise ValueError("pilot report has no ADM audit")
    samples = comparison.get("samples")
    if not isinstance(samples, list) or len(samples) < 2:
        raise ValueError("pilot ADM audit has fewer than two samples")
    ordered = sorted(samples, key=lambda entry: float(entry["athena_time"]))
    times = np.asarray(
        [float(entry["athena_time"]) for entry in ordered], dtype=np.float64
    )
    result = {}
    for group, fields in ADM_GROUPS.items():
        snapshots = []
        for entry in ordered:
            data = extract.bin_convert.read_binary(str(entry["path"]))
            snapshots.append(
                np.stack(
                    [extract.closure.assemble_uniform_grid(data, field)
                     for field in fields]
                )
            )
        result[group] = (times, np.stack(snapshots))
    return result


def _observable_groups(
    report: dict[str, object]
) -> dict[str, dict[str, tuple[np.ndarray, np.ndarray]]]:
    return {
        "history": _history_groups(Path(report["history_summary"]["path"])),
        "adm": _adm_groups(report),
    }


def compare_observable_groups(
    reference: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]],
    candidate: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]],
) -> dict[str, object]:
    return {
        "history": {
            group: compare_series(
                *reference["history"][group], *candidate["history"][group]
            )
            for group in HISTORY_GROUPS
        },
        "adm": {
            group: compare_series(*reference["adm"][group], *candidate["adm"][group])
            for group in ADM_GROUPS
        },
    }


def compare_case_reports(
    reference: dict[str, object], candidate: dict[str, object]
) -> dict[str, object]:
    return compare_observable_groups(
        _observable_groups(reference), _observable_groups(candidate)
    )


def _condition_maximum(observed: float, threshold: float) -> dict[str, object]:
    finite = math.isfinite(observed)
    return {
        "relation": "maximum",
        "observed": observed if finite else None,
        "threshold": threshold,
        "passed": finite and observed <= threshold,
    }


def _condition_minimum(observed: float, threshold: float) -> dict[str, object]:
    finite = math.isfinite(observed)
    return {
        "relation": "minimum",
        "observed": observed if finite else None,
        "threshold": threshold,
        "passed": finite and observed >= threshold,
    }


def _relative_l2(entry: dict[str, object]) -> float:
    value = entry["relative_l2"]
    return math.inf if value is None else float(value)


def observed_spatial_orders(
    comparisons: dict[str, dict[str, object]],
    resolution_factors: list[int],
    cadence_stride: int,
) -> dict[str, dict[str, list[float]]]:
    factors = sorted(resolution_factors)
    result: dict[str, dict[str, list[float]]] = {"history": {}, "adm": {}}
    if len(factors) < 3:
        return result
    for family, groups in (("history", HISTORY_GROUPS), ("adm", ADM_GROUPS)):
        for group in groups:
            errors = [
                _relative_l2(
                    comparisons[f"mr{factor}_s{cadence_stride}"]["comparison"][
                        family
                    ][group]
                )
                for factor in factors[:-1]
            ]
            orders = []
            for left, right, left_factor, right_factor in zip(
                errors[:-1], errors[1:], factors[:-2], factors[1:-1]
            ):
                if left > 0.0 and right > 0.0 and all(
                    math.isfinite(value) for value in (left, right)
                ):
                    orders.append(
                        math.log(left / right)
                        / math.log(right_factor / left_factor)
                    )
            result[family][group] = orders
    return result


def assess_convergence(
    case_documents: dict[str, dict[str, object]],
    comparisons: dict[str, dict[str, object]],
    resolution_factors: list[int],
    cadence_strides: list[int],
    arguments: argparse.Namespace,
) -> dict[str, object]:
    nonreference = [entry for entry in comparisons.values() if not entry["reference"]]
    finest_cadence = min(cadence_strides)
    spatial_orders = observed_spatial_orders(
        comparisons, resolution_factors, finest_cadence
    )
    flattened_orders = [
        order
        for family in spatial_orders.values()
        for orders in family.values()
        for order in orders
    ]
    conditions = {
        "metric_relative_l2": _condition_maximum(
            max(_relative_l2(entry["comparison"]["adm"]["metric"])
                for entry in nonreference),
            arguments.maximum_metric_relative_l2,
        ),
        "extrinsic_curvature_relative_l2": _condition_maximum(
            max(
                _relative_l2(
                    entry["comparison"]["adm"]["extrinsic_curvature"]
                )
                for entry in nonreference),
            arguments.maximum_k_relative_l2,
        ),
        "mdot_relative_l2": _condition_maximum(
            max(_relative_l2(entry["comparison"]["history"]["mdot"])
                for entry in nonreference),
            arguments.maximum_mdot_relative_l2,
        ),
        "force_relative_l2": _condition_maximum(
            max(
                _relative_l2(group)
                for entry in nonreference
                for name, group in entry["comparison"]["history"].items()
                if name != "mdot"
            ),
            arguments.maximum_force_relative_l2,
        ),
    }
    if len(resolution_factors) >= 3:
        conditions["minimum_observed_spatial_order"] = _condition_minimum(
            min(flattened_orders) if flattened_orders else -math.inf,
            arguments.minimum_observed_spatial_order,
        )
    completed = all(
        document.get("run_status") == "completed"
        for document in case_documents.values()
    )
    structural = completed and all(
        bool(document["structural_assessment"]["passed"])
        for document in case_documents.values()
    )
    runtime = completed and all(
        bool(document["runtime_assessment"]["passed"])
        for document in case_documents.values()
    )
    production_matrix = len(resolution_factors) >= 3 and len(cadence_strides) >= 2
    convergence_passed = all(entry["passed"] for entry in conditions.values())
    blockers = []
    if not production_matrix:
        blockers.append(
            "production claim requires at least three resolutions and two cadences"
        )
    if not structural:
        blockers.append("one or more pilot structural gates failed")
    if not runtime:
        blockers.append("one or more pilot runtime gates failed")
    if not convergence_passed:
        blockers.append("one or more observable convergence thresholds failed")
    return {
        "passed": completed and structural and runtime and convergence_passed,
        "observable_convergence_passed": convergence_passed,
        "conditions": conditions,
        "observed_spatial_orders": spatial_orders,
        "all_cases_completed": completed,
        "all_structural_gates_passed": structural,
        "all_runtime_gates_passed": runtime,
        "production_matrix": production_matrix,
        "science_ready": not blockers,
        "science_blockers": blockers,
    }


def _case_summary_matches(
    document: dict[str, object], expected_control: dict[str, object]
) -> bool:
    return (
        document.get("classification") == pilot.CLASSIFICATION
        and document.get("run_status") == "completed"
        and document.get("convergence_control") == expected_control
    )


def _run_case(
    arguments: argparse.Namespace,
    case: ConvergenceCase,
    fluid_cells: int,
    directory: Path,
) -> dict[str, object]:
    metric_cells = fluid_cells * case.resolution_factor
    metric_halo = pilot.minimum_metric_halo(
        metric_cells, fluid_cells, arguments.mesh_nghost
    ) + arguments.metric_halo_padding
    control = {
        "resolution_factor": case.resolution_factor,
        "cadence_stride": case.cadence_stride,
        "metric_active_cells_per_axis": metric_cells,
        "metric_halo": metric_halo,
        "mesh_nghost": arguments.mesh_nghost,
        "secondary_mass": arguments.secondary_mass,
        "secondary_chi": arguments.secondary_chi,
        "hybrid_primary_adm": arguments.hybrid_primary_adm,
        "secondary_metric_fd_step": arguments.secondary_metric_fd_step,
        "fluid_cells_per_axis": arguments.fluid_cells_per_axis,
        "meshblock_cells_per_axis": arguments.meshblock_cells_per_axis,
        "amr_levels": arguments.amr_levels,
        "refinement_radius": arguments.refinement_radius,
        "history_samples": arguments.history_samples,
        "adm_audit_samples": arguments.adm_audit_samples,
        "tlim": arguments.tlim,
        "cfl": arguments.cfl,
    }
    summary_path = directory / "summary.json"
    if directory.exists():
        if arguments.resume and summary_path.exists():
            document = json.loads(summary_path.read_text(encoding="utf-8"))
            if _case_summary_matches(document, control):
                return document
        if arguments.resume and not summary_path.exists():
            suffix = 1
            archive = directory.with_name(directory.name + f".incomplete{suffix}")
            while archive.exists():
                suffix += 1
                archive = directory.with_name(
                    directory.name + f".incomplete{suffix}"
                )
            directory.rename(archive)
        else:
            raise FileExistsError(
                f"refusing to overwrite incomplete or incompatible case {directory}"
            )
    case_arguments = argparse.Namespace(
        campaign=arguments.campaign,
        case=arguments.case,
        athena=arguments.athena,
        input=arguments.input,
        workdir=directory,
        secondary_mass=arguments.secondary_mass,
        secondary_chi=arguments.secondary_chi,
        tlim=arguments.tlim,
        cfl=arguments.cfl,
        history_samples=arguments.history_samples,
        adm_audit_samples=arguments.adm_audit_samples,
        numerical_adm_volume=True,
        hybrid_primary_adm=arguments.hybrid_primary_adm,
        secondary_metric_fd_step=arguments.secondary_metric_fd_step,
        fluid_cells_per_axis=arguments.fluid_cells_per_axis,
        meshblock_cells_per_axis=arguments.meshblock_cells_per_axis,
        amr_levels=arguments.amr_levels,
        refinement_radius=arguments.refinement_radius,
        metric_cells_per_axis=metric_cells,
        metric_cadence_stride=case.cadence_stride,
        metric_halo=metric_halo,
        mesh_nghost=arguments.mesh_nghost,
        minimum_horizon_cells=arguments.minimum_horizon_cells,
        minimum_boundary_horizon_radii=arguments.minimum_boundary_horizon_radii,
        maximum_boundary_magnetization_proxy=(
            arguments.maximum_boundary_magnetization_proxy
        ),
        maximum_initial_volume_flux_divergence=(
            arguments.maximum_initial_volume_flux_divergence
        ),
        maximum_fallback_fraction=arguments.maximum_fallback_fraction,
        maximum_boundary_flux_residual=arguments.maximum_boundary_flux_residual,
        maximum_divb=arguments.maximum_divb,
        maximum_adm_replay_relative_error=(
            arguments.maximum_adm_replay_relative_error
        ),
        allow_unsafe_structural_smoke=arguments.allow_unsafe_structural_smoke,
        fail_on_gate=False,
    )
    document = pilot.run_pilot(case_arguments)
    document["convergence_control"] = control
    _write_json(summary_path, document)
    return document


def run_convergence(arguments: argparse.Namespace) -> dict[str, object]:
    source_path, source = pilot._read_campaign(arguments.campaign)
    case_name, source_case = pilot._case_document(source, arguments.case)
    worldtube_path = Path(source_case["worldtube"]).expanduser().resolve(strict=True)
    _, faces, _ = worldtube.read_worldtube(worldtube_path)
    source_fluid_cells = faces[worldtube.FACE_NAMES[0]].cell_state.shape[-1]
    fluid_cells = (
        source_fluid_cells
        if arguments.fluid_cells_per_axis is None
        else arguments.fluid_cells_per_axis
    )
    cases = convergence_cases(
        arguments.metric_resolution_factors, arguments.metric_cadence_strides
    )
    reference = reference_case(cases)
    workdir = arguments.workdir.expanduser().resolve()
    if workdir.exists() and not arguments.resume:
        raise FileExistsError(f"refusing to reuse convergence directory {workdir}")
    workdir.mkdir(parents=True, exist_ok=arguments.resume)
    cases_directory = workdir / "cases"
    cases_directory.mkdir(exist_ok=arguments.resume)
    documents: dict[str, dict[str, object]] = {}
    case_index: dict[str, dict[str, object]] = {}
    failures = {}
    for case in cases:
        start = time.perf_counter()
        try:
            document = _run_case(
                arguments, case, fluid_cells, cases_directory / case.name
            )
            documents[case.name] = document
            case_index[case.name] = {
                "resolution_factor": case.resolution_factor,
                "cadence_stride": case.cadence_stride,
                "metric_active_cells_per_axis": (
                    fluid_cells * case.resolution_factor
                ),
                "summary": str(cases_directory / case.name / "summary.json"),
                "run_status": document["run_status"],
                "wall_seconds": document.get("wall_seconds"),
            }
        except (KeyError, OSError, RuntimeError, ValueError) as error:
            failures[case.name] = str(error)
        _write_json(
            workdir / "progress.json",
            {
                "classification": CLASSIFICATION,
                "completed_cases": sorted(documents),
                "failed_cases": failures,
                "last_case_wall_seconds": time.perf_counter() - start,
            },
        )
    refused = {
        name: str(document.get("refusal_reason", "pilot did not complete"))
        for name, document in documents.items()
        if document.get("run_status") != "completed"
    }
    failures.update(refused)
    if failures or len(documents) != len(cases):
        return {
            "classification": CLASSIFICATION,
            "run_status": "incomplete",
            "campaign": str(source_path),
            "case": case_name,
            "cases": case_index,
            "failures": failures,
            "science_ready": False,
        }
    observables = {
        name: _observable_groups(document) for name, document in documents.items()
    }
    reference_observables = observables[reference.name]
    comparisons = {}
    for case in cases:
        comparison = compare_observable_groups(
            reference_observables, observables[case.name]
        )
        comparisons[case.name] = {
            "reference": case == reference,
            "comparison": comparison,
        }
    spatial_sensitivity = {}
    for stride in arguments.metric_cadence_strides:
        slice_reference = ConvergenceCase(
            max(arguments.metric_resolution_factors), stride
        )
        spatial_sensitivity[f"s{stride}"] = {
            case.name: {
                "reference_case": slice_reference.name,
                "comparison": compare_observable_groups(
                    observables[slice_reference.name], observables[case.name]
                ),
            }
            for case in cases
            if case.cadence_stride == stride
        }
    cadence_sensitivity = {}
    finest_cadence = min(arguments.metric_cadence_strides)
    for factor in arguments.metric_resolution_factors:
        slice_reference = ConvergenceCase(factor, finest_cadence)
        cadence_sensitivity[f"mr{factor}"] = {
            case.name: {
                "reference_case": slice_reference.name,
                "comparison": compare_observable_groups(
                    observables[slice_reference.name], observables[case.name]
                ),
            }
            for case in cases
            if case.resolution_factor == factor
        }
    assessment = assess_convergence(
        documents,
        comparisons,
        arguments.metric_resolution_factors,
        arguments.metric_cadence_strides,
        arguments,
    )
    return {
        "classification": CLASSIFICATION,
        "run_status": "completed",
        "campaign": str(source_path),
        "case": case_name,
        "source_worldtube": str(worldtube_path),
        "fluid_cells_per_axis": fluid_cells,
        "source_fluid_cells_per_axis": source_fluid_cells,
        "secondary_mass": arguments.secondary_mass,
        "secondary_chi": arguments.secondary_chi,
        "metric_resolution_factors": arguments.metric_resolution_factors,
        "metric_cadence_strides": arguments.metric_cadence_strides,
        "reference_case": reference.name,
        "cases": case_index,
        "comparisons_to_reference": comparisons,
        "spatial_sensitivity_at_fixed_cadence": spatial_sensitivity,
        "cadence_sensitivity_at_fixed_resolution": cadence_sensitivity,
        "assessment": assessment,
        "limitations": [
            (
                "the fluid AMR hierarchy is shared by all cases"
                if arguments.amr_levels > 1
                else "the fluid mesh is fixed so this isolates metric representation error"
            ),
            "one-way coupling remains in force",
            "threshold passage is not a substitute for a long-duration physical run",
        ],
    }


def _unique_positive(values: list[int], name: str) -> list[int]:
    if not values or any(value < 1 for value in values):
        raise ValueError(f"{name} must contain positive integers")
    if len(set(values)) != len(values):
        raise ValueError(f"{name} contains duplicates")
    return sorted(values)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--case")
    parser.add_argument("--athena", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--secondary-mass", type=float, required=True)
    parser.add_argument("--secondary-chi", type=float, default=0.0)
    parser.add_argument("--hybrid-primary-adm", action="store_true")
    parser.add_argument("--secondary-metric-fd-step", type=float)
    parser.add_argument("--fluid-cells-per-axis", type=int)
    parser.add_argument("--meshblock-cells-per-axis", type=int)
    parser.add_argument("--amr-levels", type=int, default=1)
    parser.add_argument("--refinement-radius", type=float)
    parser.add_argument(
        "--metric-resolution-factors", type=int, nargs="+", default=(1, 2, 4)
    )
    parser.add_argument(
        "--metric-cadence-strides", type=int, nargs="+", default=(1, 2)
    )
    parser.add_argument("--metric-halo-padding", type=int, default=0)
    parser.add_argument("--mesh-nghost", type=int, default=4)
    parser.add_argument("--history-samples", type=int, default=32)
    parser.add_argument("--adm-audit-samples", type=int, default=4)
    parser.add_argument("--tlim", type=float)
    parser.add_argument("--cfl", type=float, default=0.02)
    parser.add_argument("--maximum-metric-relative-l2", type=float, default=5.0e-2)
    parser.add_argument("--maximum-k-relative-l2", type=float, default=5.0e-2)
    parser.add_argument("--maximum-mdot-relative-l2", type=float, default=1.0e-1)
    parser.add_argument("--maximum-force-relative-l2", type=float, default=1.0e-1)
    parser.add_argument(
        "--minimum-observed-spatial-order", type=float, default=5.0e-1
    )
    parser.add_argument("--minimum-horizon-cells", type=float, default=4.0)
    parser.add_argument("--minimum-boundary-horizon-radii", type=float, default=5.0)
    parser.add_argument(
        "--maximum-boundary-magnetization-proxy", type=float, default=1.0e3
    )
    parser.add_argument(
        "--maximum-initial-volume-flux-divergence", type=float, default=1.0e-10
    )
    parser.add_argument("--maximum-fallback-fraction", type=float, default=5.0e-4)
    parser.add_argument("--maximum-boundary-flux-residual", type=float, default=1.0e-10)
    parser.add_argument("--maximum-divb", type=float, default=1.0e-10)
    parser.add_argument(
        "--maximum-adm-replay-relative-error", type=float, default=1.0e-5
    )
    parser.add_argument("--allow-unsafe-structural-smoke", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--fail-on-gate", action="store_true")
    arguments = parser.parse_args()
    try:
        arguments.metric_resolution_factors = _unique_positive(
            list(arguments.metric_resolution_factors), "metric resolution factors"
        )
        arguments.metric_cadence_strides = _unique_positive(
            list(arguments.metric_cadence_strides), "metric cadence strides"
        )
    except ValueError as error:
        parser.error(str(error))
    if 1 not in arguments.metric_cadence_strides:
        parser.error("--metric-cadence-strides must include one")
    if (
        len(arguments.metric_resolution_factors)
        * len(arguments.metric_cadence_strides)
        < 2
    ):
        parser.error("convergence matrix must contain at least two cases")
    if arguments.metric_halo_padding < 0 or arguments.mesh_nghost < 1:
        parser.error("metric halo padding must be nonnegative and nghost positive")
    if arguments.history_samples < 1 or arguments.adm_audit_samples < 1:
        parser.error("history and ADM audit sample counts must be positive")
    finite_positive = (
        "secondary_mass",
        "cfl",
        "maximum_metric_relative_l2",
        "maximum_k_relative_l2",
        "maximum_mdot_relative_l2",
        "maximum_force_relative_l2",
        "minimum_observed_spatial_order",
        "minimum_horizon_cells",
        "minimum_boundary_horizon_radii",
        "maximum_boundary_magnetization_proxy",
        "maximum_initial_volume_flux_divergence",
        "maximum_boundary_flux_residual",
        "maximum_divb",
        "maximum_adm_replay_relative_error",
    )
    for name in finite_positive:
        value = getattr(arguments, name)
        if not math.isfinite(value) or value <= 0.0:
            parser.error(f"--{name.replace('_', '-')} must be finite and positive")
    if not math.isfinite(arguments.secondary_chi) or abs(arguments.secondary_chi) > 1.0:
        parser.error("--secondary-chi must lie in [-1,1]")
    if arguments.tlim is not None and (
        not math.isfinite(arguments.tlim) or arguments.tlim <= 0.0
    ):
        parser.error("--tlim must be finite and positive")
    if arguments.secondary_metric_fd_step is not None and (
        not math.isfinite(arguments.secondary_metric_fd_step)
        or arguments.secondary_metric_fd_step <= 0.0
    ):
        parser.error("--secondary-metric-fd-step must be finite and positive")
    if arguments.fluid_cells_per_axis is not None and (
        arguments.fluid_cells_per_axis < 1
    ):
        parser.error("--fluid-cells-per-axis must be positive")
    if arguments.meshblock_cells_per_axis is not None and (
        arguments.meshblock_cells_per_axis < 1
    ):
        parser.error("--meshblock-cells-per-axis must be positive")
    if arguments.amr_levels < 1:
        parser.error("--amr-levels must be positive")
    if arguments.amr_levels > 1 and not arguments.hybrid_primary_adm:
        parser.error("AMR requires --hybrid-primary-adm")
    if arguments.amr_levels > 1 and arguments.meshblock_cells_per_axis is None:
        parser.error("AMR requires --meshblock-cells-per-axis")
    if arguments.refinement_radius is not None and (
        not math.isfinite(arguments.refinement_radius)
        or arguments.refinement_radius <= 0.0
    ):
        parser.error("--refinement-radius must be finite and positive")
    if (
        not math.isfinite(arguments.maximum_fallback_fraction)
        or not 0.0 <= arguments.maximum_fallback_fraction <= 1.0
    ):
        parser.error("--maximum-fallback-fraction must lie in [0,1]")
    arguments.athena = arguments.athena.expanduser().resolve(strict=True)
    arguments.input = arguments.input.expanduser().resolve(strict=True)
    arguments.campaign = arguments.campaign.expanduser().resolve(strict=True)
    return arguments


def main() -> int:
    arguments = parse_arguments()
    report = run_convergence(arguments)
    output = arguments.workdir.expanduser().resolve() / "summary.json"
    _write_json(output, report)
    print(f"run_status={report['run_status']}")
    if report["run_status"] == "completed":
        state = "PASS" if report["assessment"]["passed"] else "FAIL"
        print(f"convergence_gates={state}")
        print(f"science_ready={report['assessment']['science_ready']}")
    print(output)
    if report["run_status"] != "completed":
        return 2
    if arguments.fail_on_gate and not report["assessment"]["passed"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
