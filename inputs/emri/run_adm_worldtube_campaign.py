#!/usr/bin/env python3
"""Build and audit an ADM-driven EMRI worldtube production campaign.

The driver builds numerical-ADM ephemeris frames, trims unsafe endpoint times
with the exact metadata preflight, extracts a convergence matrix, validates the
cubical CT cochains, and prepares bounded-memory inner replay binaries.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import re
import time
from typing import Callable

import numpy as np

import analyze_global_worldtube_convergence as comparison
import build_adm_ephemeris_frame as adm_frame
import build_kerr_ephemeris_frame as ephemeris_frame
import extract_global_worldtube as extract
import preflight_global_worldtube as preflight
import worldtube_flux_emf as worldtube


CLASSIFICATION = "athenak-emri-adm-worldtube-campaign-v1"
TEMPORAL_COVERAGE_PATTERN = re.compile(
    r"mapped global time ([^ ]+) lies outside snapshot range \[([^,]+), ([^\]]+)\]"
)


def _float_slug(value: float) -> str:
    return (
        f"{value:.8g}".replace("-", "m").replace(".", "p").replace("+", "")
    )


@dataclass(frozen=True)
class CampaignCase:
    metric_fd_multiplier: float
    cadence_stride: int
    face_cells: int
    quadrature_order: int

    @property
    def name(self) -> str:
        return (
            f"fd{_float_slug(self.metric_fd_multiplier)}_s{self.cadence_stride}_"
            f"n{self.face_cells}_q{self.quadrature_order}"
        )


def campaign_cases(
    metric_fd_multipliers: list[float],
    cadence_strides: list[int],
    face_cells: list[int],
    quadrature_orders: list[int],
) -> tuple[CampaignCase, ...]:
    result = tuple(
        CampaignCase(multiplier, stride, cells, order)
        for multiplier in metric_fd_multipliers
        for stride in cadence_strides
        for cells in face_cells
        for order in quadrature_orders
    )
    names = [case.name for case in result]
    if len(names) != len(set(names)):
        raise ValueError("campaign case names are not unique")
    return result


def reference_case(cases: tuple[CampaignCase, ...]) -> CampaignCase:
    if not cases:
        raise ValueError("campaign requires at least one case")
    target = (
        min(case.metric_fd_multiplier for case in cases),
        min(case.cadence_stride for case in cases),
        max(case.face_cells for case in cases),
        max(case.quadrature_order for case in cases),
    )
    matches = [
        case
        for case in cases
        if (
            case.metric_fd_multiplier,
            case.cadence_stride,
            case.face_cells,
            case.quadrature_order,
        )
        == target
    ]
    if len(matches) != 1:
        raise RuntimeError("campaign matrix does not contain one finest reference")
    return matches[0]


def selected_indices(sample_count: int, stride: int) -> list[int]:
    if sample_count < 2 or stride < 1:
        raise ValueError("sample_count must be at least two and stride positive")
    result = list(range(0, sample_count, stride))
    if result[-1] != sample_count - 1:
        result.append(sample_count - 1)
    return result


def subsample_scan(
    scan: extract.SnapshotManifestScan, stride: int
) -> tuple[extract.SnapshotManifestScan, list[int]]:
    indices = selected_indices(len(scan.descriptors), stride)
    return (
        extract.SnapshotManifestScan(
            path=scan.path,
            document=scan.document,
            entries=tuple(scan.entries[index] for index in indices),
            source_level=scan.source_level,
            snapshot_cache_size=scan.snapshot_cache_size,
            hash_source_files=scan.hash_source_files,
            descriptors=tuple(scan.descriptors[index] for index in indices),
        ),
        indices,
    )


def _temporal_failure_side(message: str) -> str | None:
    match = TEMPORAL_COVERAGE_PATTERN.search(message)
    if match is None:
        return None
    mapped, lower, upper = (float(value) for value in match.groups())
    if mapped < lower:
        return "leading"
    if mapped > upper:
        return "trailing"
    raise RuntimeError("temporal coverage error reports an in-range mapped time")


def trim_safe_sample_times(
    sample_times: object,
    probe: Callable[[np.ndarray], list[str]],
) -> tuple[np.ndarray, dict[str, object]]:
    requested = worldtube.validate_times(sample_times)
    retained = requested.copy()
    leading = []
    trailing = []
    iterations = 0
    while True:
        failures = probe(retained)
        if not failures:
            break
        sides = []
        for failure in failures:
            side = _temporal_failure_side(failure)
            if side is None:
                raise RuntimeError(
                    "worldtube preflight failed for a non-temporal reason: " + failure
                )
            sides.append(side)
        drop_leading = "leading" in sides
        drop_trailing = "trailing" in sides
        minimum_count = 2 + int(drop_leading) + int(drop_trailing)
        if retained.size < minimum_count:
            raise RuntimeError(
                "temporal tetrad tilt leaves fewer than two safe sample times; "
                "extend the source snapshots beyond the requested ephemeris window"
            )
        if drop_leading:
            leading.append(float(retained[0]))
            retained = retained[1:]
        if drop_trailing:
            trailing.append(float(retained[-1]))
            retained = retained[:-1]
        iterations += 1
        if iterations > requested.size:
            raise RuntimeError("safe-time trimming failed to converge")
    return retained, {
        "requested_sample_count": int(requested.size),
        "retained_sample_count": int(retained.size),
        "requested_local_time_range": [
            float(requested[0]),
            float(requested[-1]),
        ],
        "retained_local_time_range": [
            float(retained[0]),
            float(retained[-1]),
        ],
        "dropped_leading_sample_times": leading,
        "dropped_trailing_sample_times": trailing,
        "preflight_trimming_iterations": iterations,
    }


def _manifest_document(
    scan: extract.SnapshotManifestScan,
    frame_path: Path,
    sample_times: np.ndarray,
    center: np.ndarray,
    half_width: float,
    face_cells: int,
    quadrature_order: int,
    length_scale: float,
    density_renormalization: float,
) -> dict[str, object]:
    document: dict[str, object] = {
        "classification": extract.INPUT_CLASSIFICATION,
        "snapshots": [
            {
                "time": descriptor.time,
                "state": str(descriptor.state_path),
                "adm": str(descriptor.adm_path),
            }
            for descriptor in scan.descriptors
        ],
        "snapshot_cache_size": scan.snapshot_cache_size,
        "frame": str(frame_path),
        "sample_times": sample_times.tolist(),
        "worldtube": {
            "center": center.tolist(),
            "half_width": half_width,
            "cells_per_edge": face_cells,
        },
        "quadrature_order": quadrature_order,
        "hash_source_files": False,
        "global_length_in_local_units": length_scale,
        "density_renormalization": density_renormalization,
    }
    if scan.source_level is not None:
        document["source_level"] = scan.source_level
    return document


def _write_json(path: Path, document: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _projection_summary(metadata: dict[str, object]) -> dict[str, float]:
    diagnostics = metadata["sampling_diagnostics"]
    flux_endpoints = diagnostics["closed_flux_projection"]["endpoints"]
    intervals = diagnostics["edge_emf_projection"]["intervals"]
    return {
        "maximum_relative_closed_flux_correction": max(
            float(entry["relative_maximum_flux_correction"])
            for entry in flux_endpoints
        ),
        "maximum_raw_faraday_residual": max(
            float(entry["raw_maximum_faraday_residual"])
            for entry in intervals
        ),
        "maximum_edge_correction": max(
            float(entry["maximum_edge_correction"]) for entry in intervals
        ),
        "maximum_relative_edge_correction": max(
            float(entry["relative_maximum_edge_correction"])
            for entry in intervals
        ),
        "maximum_final_faraday_residual": max(
            float(entry["final_maximum_faraday_residual"])
            for entry in intervals
        ),
    }


def _compact_preflight(report: dict[str, object]) -> dict[str, object]:
    return {
        "passed": report["passed"],
        "sampled_local_event_count": report["sampled_local_event_count"],
        "mapped_global_time_range": report["mapped_global_time_range"],
        "minimum_cell_center_envelope_margin_cells": report[
            "minimum_cell_center_envelope_margin_cells"
        ],
        "minimum_additional_stencil_halo_cells": report[
            "minimum_additional_stencil_halo_cells"
        ],
        "minimum_scaled_jacobian_absolute_determinant": report[
            "minimum_scaled_jacobian_absolute_determinant"
        ],
        "maximum_jacobian_condition_number": report[
            "maximum_jacobian_condition_number"
        ],
        "warnings": report["warnings"],
    }


def _condition_maximum(observed: object, threshold: float) -> dict[str, object]:
    valid = observed is not None and math.isfinite(float(observed))
    return {
        "relation": "maximum",
        "observed": observed,
        "threshold": threshold,
        "passed": bool(valid and float(observed) <= threshold),
    }


def _condition_minimum(observed: object, threshold: float) -> dict[str, object]:
    valid = observed is not None and math.isfinite(float(observed))
    return {
        "relation": "minimum",
        "observed": observed,
        "threshold": threshold,
        "passed": bool(valid and float(observed) >= threshold),
    }


def _comparison_maximum(
    cases: dict[str, dict[str, object]], field: str
) -> float | None:
    values = []
    for case in cases.values():
        report = case.get("comparison")
        if not isinstance(report, dict):
            continue
        value = report[field]["relative_l2"]
        if value is not None:
            values.append(float(value))
    return max(values) if values else None


def assess_campaign(
    arguments: argparse.Namespace,
    frames: dict[tuple[float, int], dict[str, object]],
    cases: dict[str, dict[str, object]],
) -> dict[str, object]:
    conditions: dict[str, dict[str, object]] = {}
    for (multiplier, stride), frame_document in frames.items():
        diagnostics = frame_document["diagnostics"]
        prefix = f"frame_fd_{multiplier:.8g}_stride_{stride}"
        conditions[f"{prefix}_gram_error"] = _condition_maximum(
            diagnostics["maximum_interpolated_tetrad_gram_error"],
            arguments.maximum_frame_gram_error,
        )
        conditions[f"{prefix}_rk_difference"] = _condition_maximum(
            diagnostics["maximum_coarse_fine_spatial_leg_difference"],
            arguments.maximum_rk_leg_difference,
        )
        conditions[f"{prefix}_metric_fd_difference"] = _condition_maximum(
            diagnostics["maximum_half_step_connection_relative_difference"],
            arguments.maximum_metric_fd_relative_difference,
        )
    for name, case in cases.items():
        case_preflight = case["preflight"]
        projection = case["projection"]
        conditions[f"{name}_stencil_halo"] = _condition_minimum(
            case_preflight["minimum_additional_stencil_halo_cells"],
            arguments.minimum_stencil_halo_cells,
        )
        conditions[f"{name}_closed_flux_correction"] = _condition_maximum(
            projection["maximum_relative_closed_flux_correction"],
            arguments.maximum_relative_closed_flux_correction,
        )
        conditions[f"{name}_edge_correction"] = _condition_maximum(
            projection["maximum_relative_edge_correction"],
            arguments.maximum_relative_edge_correction,
        )
        conditions[f"{name}_faraday_residual"] = _condition_maximum(
            projection["maximum_final_faraday_residual"],
            arguments.maximum_faraday_residual,
        )
    convergence = {
        field: _comparison_maximum(cases, field)
        for field in ("state", "normal_flux", "emf")
    }
    optional_thresholds = {
        "state": arguments.maximum_state_relative_change,
        "normal_flux": arguments.maximum_flux_relative_change,
        "emf": arguments.maximum_emf_relative_change,
    }
    for field, threshold in optional_thresholds.items():
        if threshold is not None:
            conditions[f"convergence_{field}_relative_l2"] = _condition_maximum(
                convergence[field], threshold
            )
    return {
        "passed": all(condition["passed"] for condition in conditions.values()),
        "conditions": conditions,
        "maximum_relative_l2_change_from_reference": convergence,
        "note": (
            "field comparisons across metric-FD variants include the small change "
            "in the transported local frame and physical sampling surface"
        ),
    }


def _base_metric_fd_step(
    arguments: argparse.Namespace,
    scan: extract.SnapshotManifestScan,
    ephemeris_document: dict[str, object],
) -> float:
    if arguments.base_metric_fd_step is not None:
        return arguments.base_metric_fd_step
    default = 0.25 * min(
        float(np.min(descriptor.spacing)) for descriptor in scan.descriptors
    )
    return ephemeris_frame._document_float(
        ephemeris_document, "adm_metric_fd_step_global_units", default
    )


def _build_frames(
    arguments: argparse.Namespace,
    scans: dict[int, extract.SnapshotManifestScan],
    base_metric_fd_step: float,
    directory: Path,
) -> tuple[
    dict[tuple[float, int], dict[str, object]],
    dict[tuple[float, int], Path],
]:
    documents: dict[tuple[float, int], dict[str, object]] = {}
    paths: dict[tuple[float, int], Path] = {}
    shared_snapshots: dict[str, dict[str, object]] | None = None
    hashes_recorded = False
    build_index = 0
    ordered_strides = [1] + [
        stride for stride in arguments.cadence_strides if stride != 1
    ]
    for stride in ordered_strides:
        scan = scans[stride]
        stride_times = None
        for multiplier in arguments.metric_fd_multipliers:
            document = adm_frame.build_frame_from_scan(
                scan,
                arguments.ephemeris,
                metric_fd_step_override=base_metric_fd_step * multiplier,
                hash_source_files_override=(
                    scan.hash_source_files if build_index == 0 else False
                ),
            )
            generator = document["generator"]
            if build_index == 0:
                shared_snapshots = {
                    str(entry["state"]): entry
                    for entry in generator["source_snapshots"]
                }
                hashes_recorded = bool(generator["source_file_hashes_recorded"])
            else:
                if shared_snapshots is None:
                    raise RuntimeError("source provenance was not initialized")
                generator["source_snapshots"] = [
                    shared_snapshots[str(descriptor.state_path)]
                    for descriptor in scan.descriptors
                ]
                generator["source_file_hashes_recorded"] = hashes_recorded
            key = (multiplier, stride)
            path = directory / (
                f"fd_{_float_slug(multiplier)}_stride_{stride}.json"
            )
            _write_json(path, document)
            documents[key] = document
            paths[key] = path
            times = np.asarray(document["times"], dtype=np.float64)
            if stride_times is None:
                stride_times = times
            elif not np.array_equal(times, stride_times):
                raise RuntimeError(
                    "metric-FD frame variants have different knot times"
                )
            build_index += 1
    return documents, paths


def _safe_time_probe(
    arguments: argparse.Namespace,
    scans: dict[int, extract.SnapshotManifestScan],
    frame_paths: dict[tuple[float, int], Path],
    length_scale: float,
    times: np.ndarray,
    directory: Path,
) -> list[str]:
    failures = []
    probe_cells = min(arguments.face_cells)
    for (multiplier, stride), frame_path in frame_paths.items():
        scan = scans[stride]
        for order in arguments.quadrature_orders:
            path = directory / (
                f"fd{multiplier:.8g}_s{stride}_n{probe_cells}_q{order}.json"
            )
            document = _manifest_document(
                scan,
                frame_path,
                times,
                arguments.center,
                arguments.half_width,
                probe_cells,
                order,
                length_scale,
                arguments.density_renormalization,
            )
            _write_json(path, document)
            try:
                preflight.audit_manifest(path)
            except (OSError, RuntimeError, ValueError) as error:
                failures.append(str(error))
    return failures


def _run_extraction_case(
    arguments: argparse.Namespace,
    scan: extract.SnapshotManifestScan,
    case: CampaignCase,
    frame_path: Path,
    safe_times: np.ndarray,
    length_scale: float,
    directory: Path,
) -> dict[str, object]:
    directory.mkdir()
    indices = selected_indices(safe_times.size, case.cadence_stride)
    times = safe_times[indices]
    manifest_path = directory / "manifest.json"
    manifest = _manifest_document(
        scan,
        frame_path,
        times,
        arguments.center,
        arguments.half_width,
        case.face_cells,
        case.quadrature_order,
        length_scale,
        arguments.density_renormalization,
    )
    _write_json(manifest_path, manifest)
    preflight_report = preflight.audit_manifest(manifest_path)
    start = time.perf_counter()
    extracted_times, faces, metadata, diagnostics = extract.extract_manifest(
        manifest_path
    )
    elapsed = time.perf_counter() - start
    output = directory / "worldtube.npz"
    validation = worldtube.write_worldtube(output, extracted_times, faces, metadata)
    inner = directory / "inner.bin"
    inner_validation = worldtube.write_inner_binary(
        inner, extracted_times, faces, metadata
    )
    _write_json(directory / "extraction.json", diagnostics)
    return {
        "metric_fd_multiplier": case.metric_fd_multiplier,
        "metric_fd_step_global_units": (
            arguments.base_metric_fd_step_resolved * case.metric_fd_multiplier
        ),
        "cadence_stride": case.cadence_stride,
        "face_cells_per_edge": case.face_cells,
        "quadrature_order": case.quadrature_order,
        "selected_safe_time_indices": indices,
        "sample_times": extracted_times.tolist(),
        "manifest": str(manifest_path),
        "frame": str(frame_path),
        "worldtube": str(output),
        "inner_binary": str(inner),
        "extraction_wall_seconds": elapsed,
        "preflight": _compact_preflight(preflight_report),
        "projection": _projection_summary(metadata),
        "raw_sampling": diagnostics["raw_sampling"],
        "snapshot_loading": diagnostics["snapshot_loading"],
        "worldtube_validation": validation,
        "inner_validation": inner_validation,
    }


def run_campaign(arguments: argparse.Namespace) -> dict[str, object]:
    arguments.manifest = arguments.manifest.expanduser().resolve(strict=True)
    arguments.ephemeris = arguments.ephemeris.expanduser().resolve(strict=True)
    workdir = arguments.workdir.expanduser().resolve()
    if workdir.exists():
        raise FileExistsError(f"refusing to reuse existing campaign directory {workdir}")
    scan = extract.scan_snapshot_manifest(arguments.manifest)
    ephemeris_document, _ = ephemeris_frame.read_ephemeris_document(
        arguments.ephemeris
    )
    length_scale = ephemeris_frame._document_float(
        ephemeris_document, "global_length_in_local_units", 1.0
    )
    base_fd_step = _base_metric_fd_step(arguments, scan, ephemeris_document)
    if not math.isfinite(base_fd_step) or base_fd_step <= 0.0:
        raise ValueError("resolved base ADM metric finite-difference step is invalid")
    arguments.base_metric_fd_step_resolved = base_fd_step
    cases_tuple = campaign_cases(
        arguments.metric_fd_multipliers,
        arguments.cadence_strides,
        arguments.face_cells,
        arguments.quadrature_orders,
    )
    reference = reference_case(cases_tuple)
    scans = {}
    source_snapshot_indices = {}
    for stride in arguments.cadence_strides:
        scans[stride], source_snapshot_indices[stride] = subsample_scan(
            scan, stride
        )
    workdir.mkdir(parents=True)
    frame_documents, frame_paths = _build_frames(
        arguments, scans, base_fd_step, workdir / "frames"
    )
    requested_times = np.asarray(
        frame_documents[(reference.metric_fd_multiplier, 1)]["times"],
        dtype=np.float64,
    )
    probe_directory = workdir / "safe_time_probes"
    probe_directory.mkdir()

    def probe(times: np.ndarray) -> list[str]:
        return _safe_time_probe(
            arguments,
            scans,
            frame_paths,
            length_scale,
            times,
            probe_directory,
        )

    safe_times, safe_time_report = trim_safe_sample_times(requested_times, probe)
    cases: dict[str, dict[str, object]] = {}
    cases_directory = workdir / "cases"
    cases_directory.mkdir()
    for case in cases_tuple:
        result = _run_extraction_case(
            arguments,
            scans[case.cadence_stride],
            case,
            frame_paths[(case.metric_fd_multiplier, case.cadence_stride)],
            safe_times,
            length_scale,
            cases_directory / case.name,
        )
        cases[case.name] = result
        _write_json(cases_directory / case.name / "summary.json", result)

    reference_path = Path(cases[reference.name]["worldtube"])
    reference_times, reference_faces, reference_metadata = worldtube.read_worldtube(
        reference_path
    )
    resampled_references: dict[int, dict[str, worldtube.FaceData]] = {}
    for case in cases_tuple:
        result = cases[case.name]
        candidate_path = Path(result["worldtube"])
        candidate_times, candidate_faces, candidate_metadata = (
            worldtube.read_worldtube(candidate_path)
        )
        if case.face_cells not in resampled_references:
            if case.face_cells == reference.face_cells:
                resampled_references[case.face_cells] = reference_faces
            else:
                resampled_references[case.face_cells] = worldtube.resample_worldtube(
                    reference_times, reference_faces, case.face_cells
                )
        report = comparison.compare_worldtubes(
            reference_times,
            resampled_references[case.face_cells],
            candidate_times,
            candidate_faces,
            state_variables=candidate_metadata.get(
                "state_variables", reference_metadata.get("state_variables", [])
            ),
        )
        report["reference_case"] = reference.name
        report["reference_worldtube"] = str(reference_path)
        report["candidate_worldtube"] = str(candidate_path)
        report["reference_resampled_face_cells"] = case.face_cells
        result["comparison"] = report
        _write_json(cases_directory / case.name / "comparison.json", report)
        _write_json(cases_directory / case.name / "summary.json", result)

    assessment = assess_campaign(arguments, frame_documents, cases)
    report = {
        "classification": CLASSIFICATION,
        "source_manifest": str(scan.path),
        "source_manifest_sha256": extract.static.file_sha256(scan.path),
        "source_snapshot_count": len(scan.descriptors),
        "source_level": scan.descriptors[0].source_level,
        "source_storage": scan.descriptors[0].source_storage,
        "source_time_range_global_units": [
            scan.descriptors[0].time,
            scan.descriptors[-1].time,
        ],
        "ephemeris": str(arguments.ephemeris),
        "ephemeris_sha256": extract.static.file_sha256(arguments.ephemeris),
        "global_length_in_local_units": length_scale,
        "density_renormalization": arguments.density_renormalization,
        "worldtube": {
            "center": arguments.center.tolist(),
            "half_width": arguments.half_width,
            "face_cells": arguments.face_cells,
        },
        "base_metric_fd_step_global_units": base_fd_step,
        "metric_fd_multipliers": arguments.metric_fd_multipliers,
        "cadence_strides": arguments.cadence_strides,
        "source_snapshot_indices_by_cadence_stride": {
            str(stride): indices
            for stride, indices in source_snapshot_indices.items()
        },
        "quadrature_orders": arguments.quadrature_orders,
        "safe_time_window": safe_time_report,
        "reference_case": reference.name,
        "frames": {
            f"fd_{multiplier:.17g}_stride_{stride}": {
                "path": str(frame_paths[(multiplier, stride)]),
                "metric_fd_step_global_units": base_fd_step * multiplier,
                "cadence_stride": stride,
                "diagnostics": document["diagnostics"],
            }
            for (multiplier, stride), document in frame_documents.items()
        },
        "cases": cases,
        "assessment": assessment,
        "limitations": [
            "one-way coupling; the local accretor does not modify the global disk",
            (
                "snapshot interpolation remains trilinear in space and linear in "
                "time; the campaign measures but does not remove this error"
            ),
            (
                "comparisons across metric-FD variants include the corresponding "
                "change in the transported frame and sampling surface"
            ),
            (
                "inner binaries are prepared and structurally validated here; a "
                "separate AthenaK replay is still required before science use"
            ),
        ],
    }
    return report


def _unique_positive_ints(values: list[int], name: str) -> list[int]:
    if not values or any(value < 1 for value in values):
        raise ValueError(f"{name} must contain positive integers")
    if len(values) != len(set(values)):
        raise ValueError(f"{name} contains duplicate values")
    return values


def _unique_positive_floats(values: list[float], name: str) -> list[float]:
    if not values or any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise ValueError(f"{name} must contain finite positive values")
    if len(values) != len(set(values)):
        raise ValueError(f"{name} contains duplicate values")
    return values


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--ephemeris", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--half-width", type=float, required=True)
    parser.add_argument("--center", type=float, nargs=3, default=(0.0, 0.0, 0.0))
    parser.add_argument("--face-cells", type=int, nargs="+", default=(32,))
    parser.add_argument(
        "--metric-fd-multipliers", type=float, nargs="+", default=(1.0, 0.5)
    )
    parser.add_argument("--cadence-strides", type=int, nargs="+", default=(1, 2))
    parser.add_argument(
        "--quadrature-orders", type=int, nargs="+", default=(2, 3)
    )
    parser.add_argument("--base-metric-fd-step", type=float)
    parser.add_argument("--density-renormalization", type=float, default=1.0)
    parser.add_argument("--minimum-stencil-halo-cells", type=int, default=1)
    parser.add_argument("--maximum-frame-gram-error", type=float, default=1.0e-8)
    parser.add_argument("--maximum-rk-leg-difference", type=float, default=1.0e-8)
    parser.add_argument(
        "--maximum-metric-fd-relative-difference", type=float, default=1.0e-1
    )
    parser.add_argument(
        "--maximum-relative-closed-flux-correction", type=float, default=1.0e-10
    )
    parser.add_argument(
        "--maximum-relative-edge-correction", type=float, default=5.0e-2
    )
    parser.add_argument("--maximum-faraday-residual", type=float, default=1.0e-10)
    parser.add_argument("--maximum-state-relative-change", type=float)
    parser.add_argument("--maximum-flux-relative-change", type=float)
    parser.add_argument("--maximum-emf-relative-change", type=float)
    parser.add_argument("--fail-on-gate", action="store_true")
    arguments = parser.parse_args()
    try:
        arguments.face_cells = _unique_positive_ints(
            list(arguments.face_cells), "face_cells"
        )
        arguments.metric_fd_multipliers = _unique_positive_floats(
            list(arguments.metric_fd_multipliers), "metric_fd_multipliers"
        )
        arguments.cadence_strides = _unique_positive_ints(
            list(arguments.cadence_strides), "cadence_strides"
        )
        arguments.quadrature_orders = _unique_positive_ints(
            list(arguments.quadrature_orders), "quadrature_orders"
        )
    except ValueError as error:
        parser.error(str(error))
    if min(arguments.face_cells) < 2:
        parser.error("--face-cells entries must be at least two")
    if max(arguments.quadrature_orders) > 8:
        parser.error("--quadrature-orders cannot exceed eight")
    if 1 not in arguments.cadence_strides:
        parser.error("--cadence-strides must include one for the reference")
    positive_scalars = (
        arguments.half_width,
        arguments.density_renormalization,
    )
    if any(not math.isfinite(value) or value <= 0.0 for value in positive_scalars):
        parser.error("half width and density renormalization must be positive")
    if arguments.base_metric_fd_step is not None and (
        not math.isfinite(arguments.base_metric_fd_step)
        or arguments.base_metric_fd_step <= 0.0
    ):
        parser.error("--base-metric-fd-step must be finite and positive")
    if arguments.minimum_stencil_halo_cells < 0:
        parser.error("--minimum-stencil-halo-cells cannot be negative")
    threshold_names = (
        "maximum_frame_gram_error",
        "maximum_rk_leg_difference",
        "maximum_metric_fd_relative_difference",
        "maximum_relative_closed_flux_correction",
        "maximum_relative_edge_correction",
        "maximum_faraday_residual",
        "maximum_state_relative_change",
        "maximum_flux_relative_change",
        "maximum_emf_relative_change",
    )
    for name in threshold_names:
        value = getattr(arguments, name)
        if value is not None and (not math.isfinite(value) or value < 0.0):
            parser.error(f"--{name.replace('_', '-')} must be finite and nonnegative")
    center = np.asarray(arguments.center, dtype=np.float64)
    if center.shape != (3,) or not np.isfinite(center).all():
        parser.error("--center must contain three finite values")
    arguments.center = center
    return arguments


def main() -> int:
    arguments = parse_arguments()
    report = run_campaign(arguments)
    output = arguments.workdir.expanduser().resolve() / "summary.json"
    _write_json(output, report)
    assessment = report["assessment"]
    print(f"reference_case={report['reference_case']}")
    print(
        "safe_local_time_range="
        f"{report['safe_time_window']['retained_local_time_range']}"
    )
    print(f"gates={'PASS' if assessment['passed'] else 'FAIL'}")
    print(output)
    return 1 if arguments.fail_on_gate and not assessment["passed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
