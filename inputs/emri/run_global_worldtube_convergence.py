#!/usr/bin/env python3
"""Run snapshot-cadence and quadrature convergence for offline EMRI worldtubes."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import time
from typing import Any

import numpy as np

import analyze_global_worldtube_convergence as analysis
import compare_worldtube_closure as closure
import extract_global_worldtube as extract
import worldtube_flux_emf as worldtube


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "inputs" / "emri" / "emri_windtunnel_smoke.athinput"
CLASSIFICATION = "athenak-emri-global-worldtube-convergence-campaign-v1"
SMOKE_GATE_THRESHOLDS = {
    "maximum_state_relative_l2": 1.0e-6,
    "maximum_normal_flux_relative_l2": 2.0e-4,
    "maximum_emf_relative_l2": 1.0e-2,
    "maximum_relative_closed_flux_correction": 1.0e-12,
    "maximum_relative_edge_correction": 2.0e-2,
    "maximum_final_faraday_residual": 1.0e-12,
    "minimum_flux_spatial_order": 1.5,
    "minimum_emf_spatial_order": 1.5,
    "maximum_flux_cadence_sensitivity": 3.0e-1,
    "maximum_emf_cadence_sensitivity": 2.0e-2,
    "maximum_replay_fallback_fraction": 5.0e-4,
    "maximum_replay_boundary_flux_residual": 1.0e-12,
    "maximum_replay_density_relative_l2": 1.0e-4,
    "maximum_replay_pressure_relative_l2": 5.0e-4,
    "maximum_replay_velocity_relative_l2": 2.0e-4,
    "maximum_replay_magnetic_relative_l2": 2.0e-4,
}


def _run(command: list[str], directory: Path, log_path: Path) -> float:
    start = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            command,
            cwd=directory,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    elapsed = time.perf_counter() - start
    if completed.returncode != 0:
        tail = log_path.read_text(encoding="utf-8", errors="replace").splitlines()[
            -50:
        ]
        raise RuntimeError(
            f"AthenaK failed with exit code {completed.returncode}:\n"
            + "\n".join(tail)
        )
    return elapsed


def _largest_divisor_at_most(value: int, limit: int) -> int:
    for candidate in range(min(value, limit), 0, -1):
        if value % candidate == 0:
            return candidate
    raise RuntimeError("positive integer unexpectedly has no divisor")


def selected_indices(sample_count: int, stride: int) -> list[int]:
    if sample_count < 2 or stride < 1:
        raise ValueError("sample_count must be at least two and stride positive")
    result = list(range(0, sample_count, stride))
    if result[-1] != sample_count - 1:
        result.append(sample_count - 1)
    return result


def _mesh_arguments(cells: int, half_width: float, meshblock: int) -> list[str]:
    return [
        f"mesh/nx1={cells}",
        f"mesh/nx2={cells}",
        f"mesh/nx3={cells}",
        f"mesh/x1min={-half_width:.17g}",
        f"mesh/x1max={half_width:.17g}",
        f"mesh/x2min={-half_width:.17g}",
        f"mesh/x2max={half_width:.17g}",
        f"mesh/x3min={-half_width:.17g}",
        f"mesh/x3max={half_width:.17g}",
        f"meshblock/nx1={meshblock}",
        f"meshblock/nx2={meshblock}",
        f"meshblock/nx3={meshblock}",
    ]


def _evolution_arguments(tlim: float, cfl: float, output_every_step: bool) -> list[str]:
    output_dt = 1.0e-30 if output_every_step else tlim
    return [
        "time/integrator=rk3",
        "time/nlim=-1",
        f"time/tlim={tlim:.17g}",
        f"time/cfl_number={cfl:.17g}",
        "problem/user_hist=false",
        "output1/variable=mhd_w_bcc",
        f"output1/dt={output_dt:.17g}",
        "output2/variable=adm",
        f"output2/dt={output_dt:.17g}",
        "output3/dt=0",
        "output4/dt=0",
    ]


def _read_dump_metadata(path: Path) -> tuple[float, int]:
    data = extract.bin_convert.read_binary(str(path), variables=())
    return float(data["time"]), int(data["cycle"])


def _paired_snapshots(directory: Path) -> list[dict[str, object]]:
    state_files = sorted((directory / "bin").glob("outer.mhd_w_bcc.*.bin"))
    adm_files = sorted((directory / "bin").glob("outer.adm.*.bin"))
    if len(state_files) < 2 or len(state_files) != len(adm_files):
        raise RuntimeError("outer run did not produce paired fluid and ADM time series")
    result = []
    for state_path, adm_path in zip(state_files, adm_files, strict=True):
        state_time, state_cycle = _read_dump_metadata(state_path)
        adm_time, adm_cycle = _read_dump_metadata(adm_path)
        scale = max(1.0, abs(state_time), abs(adm_time))
        if state_cycle != adm_cycle or not math.isclose(
            state_time,
            adm_time,
            rel_tol=0.0,
            abs_tol=128.0 * np.finfo(float).eps * scale,
        ):
            raise RuntimeError("outer fluid and ADM dump sequences are not co-temporal")
        result.append(
            {
                "time": state_time,
                "cycle": state_cycle,
                "state": state_path,
                "adm": adm_path,
            }
        )
    return result


def identity_frame_document(times: list[float]) -> dict[str, object]:
    if len(times) < 2 or np.any(np.diff(times) <= 0.0):
        raise ValueError("identity frame requires increasing times")
    legs = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
    zero_derivative = [[0.0, 0.0, 0.0] for _ in range(4)]
    return {
        "classification": extract.FRAME_SERIES_CLASSIFICATION,
        "times": times,
        "worldline": [[value, 0.0, 0.0, 0.0] for value in times],
        "worldline_tangent": [[1.0, 0.0, 0.0, 0.0] for _ in times],
        "spatial_legs": [legs for _ in times],
        "spatial_leg_derivative": [zero_derivative for _ in times],
    }


def _write_extraction_manifest(
    path: Path,
    snapshots: list[dict[str, object]],
    indices: list[int],
    half_width: float,
    face_cells: int,
    quadrature_order: int,
) -> None:
    selected = [snapshots[index] for index in indices]
    times = [float(entry["time"]) for entry in selected]
    document = {
        "classification": extract.INPUT_CLASSIFICATION,
        "snapshots": [
            {"state": str(entry["state"]), "adm": str(entry["adm"])}
            for entry in selected
        ],
        "frame": identity_frame_document(times),
        "sample_times": times,
        "worldtube": {
            "center": [0.0, 0.0, 0.0],
            "half_width": half_width,
            "cells_per_edge": face_cells,
        },
        "quadrature_order": quadrature_order,
        "hash_source_files": False,
        "global_length_in_local_units": 1.0,
        "density_renormalization": 1.0,
    }
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


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
            float(entry["raw_maximum_faraday_residual"]) for entry in intervals
        ),
        "maximum_edge_correction": max(
            float(entry["maximum_edge_correction"]) for entry in intervals
        ),
        "maximum_relative_edge_correction": max(
            float(entry["relative_maximum_edge_correction"])
            for entry in intervals
        ),
        "maximum_final_faraday_residual": max(
            float(entry["final_maximum_faraday_residual"]) for entry in intervals
        ),
    }


def _characteristic_counts(log_path: Path) -> tuple[int, int] | None:
    match = re.search(
        r"characteristic GR boundary: projections=(\d+), fallbacks=(\d+)",
        log_path.read_text(encoding="utf-8", errors="replace"),
    )
    if match is None:
        return None
    return int(match.group(1)), int(match.group(2))


def _flux_residual(log_path: Path) -> float | None:
    matches = re.findall(
        r"maximum normalized boundary-flux residual=([0-9.eE+-]+)",
        log_path.read_text(encoding="utf-8", errors="replace"),
    )
    return None if not matches else max(float(value) for value in matches)


def _final_output(directory: Path, stem: str) -> Path:
    candidates = sorted((directory / "bin").glob(f"{stem}.mhd_w_bcc.*.bin"))
    if not candidates:
        raise RuntimeError(f"no mhd_w_bcc output found below {directory}")
    return candidates[-1]


def _run_inner_replay(
    arguments: argparse.Namespace,
    directory: Path,
    replay: Path,
    face_cells: int,
    reference: Path,
) -> dict[str, object]:
    directory.mkdir()
    meshblock = _largest_divisor_at_most(face_cells, max(1, face_cells // 2))
    log_path = directory.with_suffix(".log")
    force_radii = (0.55, 0.65, 0.75, 0.875)
    command = [
        str(arguments.athena),
        "-i",
        str(arguments.input),
        "-d",
        str(directory),
        "job/basename=inner",
        *_mesh_arguments(face_cells, arguments.half_width, meshblock),
        *_evolution_arguments(arguments.tlim, arguments.cfl, False),
        f"problem/force_surface_radius={force_radii[0] * arguments.half_width:.17g}",
        f"problem/force_outer_radius_1={force_radii[1] * arguments.half_width:.17g}",
        f"problem/force_outer_radius_2={force_radii[2] * arguments.half_width:.17g}",
        f"problem/force_outer_radius_3={force_radii[3] * arguments.half_width:.17g}",
        "emri_worldtube/enabled=true",
        "emri_worldtube/mode=inner",
        f"emri_worldtube/file={replay}",
        "emri_worldtube/fluid_boundary=characteristic_gr",
        "emri_worldtube/fluid_state_frame=inner_coordinate",
    ]
    wall_seconds = _run(command, ROOT, log_path)
    candidate = _final_output(directory, "inner")
    counts = _characteristic_counts(log_path)
    return {
        "wall_seconds": wall_seconds,
        "candidate_output": str(candidate),
        "projection_count": None if counts is None else counts[0],
        "fallback_count": None if counts is None else counts[1],
        "maximum_boundary_flux_residual": _flux_residual(log_path),
        "closure": closure.compare_files(reference, candidate, arguments.gamma),
    }


def _effective_cadence(times: np.ndarray) -> dict[str, float]:
    intervals = np.diff(times)
    return {
        "minimum_dt": float(np.min(intervals)),
        "maximum_dt": float(np.max(intervals)),
        "mean_dt": float(np.mean(intervals)),
    }


def _comparison_error(case: dict[str, Any], field: str) -> float:
    value = case["comparison"][field]["relative_l2"]
    if value is None or not math.isfinite(float(value)):
        raise ValueError(f"case has no finite relative L2 error for {field}")
    return float(value)


def _observed_orders(spacings: list[float], errors: list[float]) -> list[float]:
    if len(spacings) != len(errors) or len(spacings) < 2:
        raise ValueError("observed order requires equally sized spatial series")
    result = []
    for coarse_h, fine_h, coarse_error, fine_error in zip(
        spacings[:-1], spacings[1:], errors[:-1], errors[1:], strict=True
    ):
        if not coarse_h > fine_h > 0.0 or not coarse_error > 0.0 or fine_error <= 0.0:
            raise ValueError("spacings must decrease and errors must be positive")
        result.append(
            math.log(coarse_error / fine_error) / math.log(coarse_h / fine_h)
        )
    return result


def _relative_spread(values: list[float]) -> float:
    if not values or values[0] <= 0.0:
        raise ValueError("relative spread requires a positive reference value")
    return max(abs(value - values[0]) for value in values) / values[0]


def _gate(
    conditions: dict[str, dict[str, object]],
    name: str,
    observed: float,
    threshold: float,
    relation: str,
) -> None:
    if relation not in ("maximum", "minimum"):
        raise ValueError("gate relation must be maximum or minimum")
    conditions[name] = {
        "relation": relation,
        "observed": observed,
        "threshold": threshold,
        "passed": (
            observed <= threshold if relation == "maximum" else observed >= threshold
        ),
    }


def _spatial_summary(
    cases: dict[str, Any], cells: list[int], strides: list[int], orders: list[int]
) -> dict[str, object]:
    spatial: dict[str, object] = {}
    for stride in strides:
        for order in orders:
            selected = [cases[f"n{cells_}_s{stride}_q{order}"] for cells_ in cells]
            spacings = [float(case["source_spacing"]) for case in selected]
            entry: dict[str, object] = {
                "outer_cells": cells,
                "source_spacing": spacings,
            }
            for field in ("normal_flux", "emf"):
                errors = [_comparison_error(case, field) for case in selected]
                entry[field] = {
                    "relative_l2": errors,
                    "observed_orders": (
                        _observed_orders(spacings, errors)
                        if len(spacings) >= 2
                        else []
                    ),
                }
            spatial[f"s{stride}_q{order}"] = entry
    return spatial


def _cadence_summary(
    cases: dict[str, Any], cells: list[int], strides: list[int], orders: list[int]
) -> dict[str, object]:
    cadence: dict[str, object] = {}
    for cells_ in cells:
        for order in orders:
            selected = [cases[f"n{cells_}_s{stride}_q{order}"] for stride in strides]
            entry = {}
            for field in ("normal_flux", "emf"):
                errors = [_comparison_error(case, field) for case in selected]
                entry[field] = {
                    "relative_l2": errors,
                    "relative_spread_from_finest_cadence": _relative_spread(errors),
                }
            cadence[f"n{cells_}_q{order}"] = entry
    return cadence


def _efficiency_summary(
    cases: dict[str, Any], cells: list[int], stride: int, orders: list[int]
) -> dict[str, object]:
    efficiency = {}
    for cells_ in cells:
        efficiency[str(cells_)] = {}
        for order in orders:
            case = cases[f"n{cells_}_s{stride}_q{order}"]
            efficiency[str(cells_)][str(order)] = {
                "extraction_wall_seconds": float(case["extraction_wall_seconds"]),
                "normal_flux_relative_l2": _comparison_error(case, "normal_flux"),
                "emf_relative_l2": _comparison_error(case, "emf"),
            }
    return efficiency


def _add_replay_gates(
    conditions: dict[str, dict[str, object]], replay_cases: list[dict[str, Any]]
) -> dict[str, float | int] | None:
    if not replay_cases:
        return None
    thresholds = SMOKE_GATE_THRESHOLDS
    projections = sum(int(case["projection_count"]) for case in replay_cases)
    fallbacks = sum(int(case["fallback_count"]) for case in replay_cases)
    fallback_fraction = fallbacks / projections
    _gate(
        conditions,
        "replay_fallback_fraction",
        fallback_fraction,
        thresholds["maximum_replay_fallback_fraction"],
        "maximum",
    )
    _gate(
        conditions,
        "replay_boundary_flux_residual",
        max(float(case["maximum_boundary_flux_residual"]) for case in replay_cases),
        thresholds["maximum_replay_boundary_flux_residual"],
        "maximum",
    )
    closure_paths = {
        "density": ("variables", "dens"),
        "pressure": ("variables", "press"),
        "velocity": ("vector_groups", "velocity"),
        "magnetic": ("vector_groups", "magnetic"),
    }
    for name, path in closure_paths.items():
        observed = max(
            float(case["closure"][path[0]][path[1]]["relative_l2"])
            for case in replay_cases
        )
        _gate(
            conditions,
            f"replay_{name}_relative_l2",
            observed,
            thresholds[f"maximum_replay_{name}_relative_l2"],
            "maximum",
        )
    return {
        "projection_count": projections,
        "fallback_count": fallbacks,
        "fallback_fraction": fallback_fraction,
    }


def assess_campaign(report: dict[str, Any]) -> dict[str, object]:
    """Summarize convergence and apply smooth-smoke regression gates."""

    cases = report["cases"]
    cells = sorted(int(value) for value in report["outer_cells"])
    strides = sorted(int(value) for value in report["cadence_strides"])
    orders = sorted(int(value) for value in report["quadrature_orders"])
    spatial = _spatial_summary(cases, cells, strides, orders)
    cadence = _cadence_summary(cases, cells, strides, orders)
    efficiency = _efficiency_summary(cases, cells, strides[0], orders)
    conditions: dict[str, dict[str, object]] = {}
    thresholds = SMOKE_GATE_THRESHOLDS
    for field in ("state", "normal_flux", "emf"):
        _gate(
            conditions,
            f"{field}_relative_l2",
            max(_comparison_error(case, field) for case in cases.values()),
            thresholds[f"maximum_{field}_relative_l2"],
            "maximum",
        )
    for name in (
        "maximum_relative_closed_flux_correction",
        "maximum_relative_edge_correction",
        "maximum_final_faraday_residual",
    ):
        _gate(
            conditions,
            name,
            max(float(case["projection"][name]) for case in cases.values()),
            thresholds[name],
            "maximum",
        )

    reference_series = spatial[f"s{strides[0]}_q{orders[-1]}"]
    flux_orders = reference_series["normal_flux"]["observed_orders"]
    emf_orders = reference_series["emf"]["observed_orders"]
    if flux_orders:
        _gate(
            conditions,
            "flux_spatial_order",
            min(flux_orders),
            thresholds["minimum_flux_spatial_order"],
            "minimum",
        )
        _gate(
            conditions,
            "emf_spatial_order",
            min(emf_orders),
            thresholds["minimum_emf_spatial_order"],
            "minimum",
        )
    for field in ("normal_flux", "emf"):
        maximum_spread = max(
            float(entry[field]["relative_spread_from_finest_cadence"])
            for entry in cadence.values()
        )
        threshold_name = "flux" if field == "normal_flux" else field
        _gate(
            conditions,
            f"{field}_cadence_sensitivity",
            maximum_spread,
            thresholds[f"maximum_{threshold_name}_cadence_sensitivity"],
            "maximum",
        )

    replay_cases = [
        case["inner_replay"] for case in cases.values() if "inner_replay" in case
    ]
    replay_summary = _add_replay_gates(conditions, replay_cases)
    return {
        "spatial_convergence": spatial,
        "cadence_sensitivity": cadence,
        "quadrature_efficiency_at_finest_cadence": efficiency,
        "replay": replay_summary,
        "smoke_gates": {
            "scope": (
                "regression gates for this smooth short-duration problem; physical "
                "production runs require observable-level convergence"
            ),
            "passed": all(bool(entry["passed"]) for entry in conditions.values()),
            "conditions": conditions,
        },
    }


def run_campaign(arguments: argparse.Namespace) -> dict[str, object]:
    arguments.athena = arguments.athena.expanduser().resolve(strict=True)
    arguments.input = arguments.input.expanduser().resolve(strict=True)
    workdir = arguments.workdir.expanduser().resolve()
    if workdir.exists():
        raise FileExistsError(f"refusing to reuse existing campaign directory {workdir}")
    workdir.mkdir(parents=True)
    cases: dict[str, Any] = {}
    outer_runs: dict[str, Any] = {}
    for outer_cells in arguments.outer_cells:
        face_cells = outer_cells // 2
        source_spacing = 4.0 * arguments.half_width / outer_cells
        outer_directory = workdir / f"outer_n{outer_cells}"
        outer_directory.mkdir()
        meshblock = _largest_divisor_at_most(outer_cells, 16)
        tube_basename = outer_directory / "tube"
        outer_command = [
            str(arguments.athena),
            "-i",
            str(arguments.input),
            "-d",
            str(outer_directory),
            "job/basename=outer",
            *_mesh_arguments(
                outer_cells, 2.0 * arguments.half_width, meshblock
            ),
            *_evolution_arguments(arguments.tlim, arguments.cfl, True),
            "emri_worldtube/enabled=true",
            "emri_worldtube/mode=outer",
            "emri_worldtube/overwrite=true",
            f"emri_worldtube/half_width={arguments.half_width:.17g}",
            "emri_worldtube/dcycle=1",
            f"emri_worldtube/file_basename={tube_basename}",
        ]
        outer_log = workdir / f"outer_n{outer_cells}.log"
        outer_wall = _run(outer_command, ROOT, outer_log)
        manifests = sorted(outer_directory.glob("tube.cycle*.manifest.json"))
        if len(manifests) != 1:
            raise RuntimeError("outer run did not produce exactly one worldtube manifest")
        reference_times, reference_faces, reference_metadata = (
            worldtube.read_outer_stream(manifests[0])
        )
        snapshots = _paired_snapshots(outer_directory)
        snapshot_times = np.asarray([entry["time"] for entry in snapshots])
        if not np.array_equal(snapshot_times, reference_times):
            raise RuntimeError("binary dump and online worldtube time axes differ")
        reference_output = _final_output(outer_directory, "outer")
        outer_runs[str(outer_cells)] = {
            "wall_seconds": outer_wall,
            "source_spacing": source_spacing,
            "face_cells_per_axis": face_cells,
            "sample_count": int(reference_times.size),
            "time_cadence": _effective_cadence(reference_times),
            "worldtube_manifest": str(manifests[0]),
            "reference_output": str(reference_output),
        }
        if reference_metadata.get("state_variables") != [
            "rho",
            "u1",
            "u2",
            "u3",
            "pgas",
            "bcc1",
            "bcc2",
            "bcc3",
        ]:
            raise RuntimeError("online DynGRMHD worldtube has unexpected state variables")

        for stride in arguments.cadence_strides:
            indices = selected_indices(reference_times.size, stride)
            selected_times = reference_times[indices]
            for quadrature_order in arguments.quadrature_orders:
                case_name = f"n{outer_cells}_s{stride}_q{quadrature_order}"
                case_directory = workdir / case_name
                case_directory.mkdir()
                manifest_path = case_directory / "extract.json"
                _write_extraction_manifest(
                    manifest_path,
                    snapshots,
                    indices,
                    arguments.half_width,
                    face_cells,
                    quadrature_order,
                )
                start = time.perf_counter()
                candidate_times, candidate_faces, metadata, diagnostics = (
                    extract.extract_manifest(manifest_path)
                )
                extraction_seconds = time.perf_counter() - start
                candidate_path = case_directory / "worldtube.npz"
                worldtube.write_worldtube(
                    candidate_path, candidate_times, candidate_faces, metadata
                )
                comparison = analysis.compare_worldtubes(
                    reference_times,
                    reference_faces,
                    candidate_times,
                    candidate_faces,
                    state_variables=metadata["state_variables"],
                )
                case: dict[str, Any] = {
                    "outer_cells_per_axis": outer_cells,
                    "source_spacing": source_spacing,
                    "face_cells_per_axis": face_cells,
                    "cadence_stride": stride,
                    "quadrature_order": quadrature_order,
                    "selected_endpoint_indices": indices,
                    "selected_sample_count": len(indices),
                    "selected_time_cadence": _effective_cadence(selected_times),
                    "extraction_wall_seconds": extraction_seconds,
                    "worldtube": str(candidate_path),
                    "comparison": comparison,
                    "projection": _projection_summary(metadata),
                    "raw_sampling": diagnostics["raw_sampling"],
                    "snapshot_loading": diagnostics["snapshot_loading"],
                }
                if (
                    not arguments.skip_replay
                    and quadrature_order == arguments.replay_quadrature
                ):
                    replay = case_directory / "inner.bin"
                    worldtube.write_inner_binary(
                        replay, candidate_times, candidate_faces, metadata
                    )
                    case["inner_replay"] = _run_inner_replay(
                        arguments,
                        case_directory / "inner",
                        replay,
                        face_cells,
                        reference_output,
                    )
                cases[case_name] = case
                (case_directory / "summary.json").write_text(
                    json.dumps(case, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )

    report = {
        "classification": CLASSIFICATION,
        "athena": str(arguments.athena),
        "input": str(arguments.input),
        "half_width": arguments.half_width,
        "tlim": arguments.tlim,
        "cfl_number": arguments.cfl,
        "outer_cells": arguments.outer_cells,
        "cadence_strides": arguments.cadence_strides,
        "quadrature_orders": arguments.quadrature_orders,
        "replay_quadrature": arguments.replay_quadrature,
        "outer_runs": outer_runs,
        "cases": cases,
        "caveat": (
            "identity/static extraction isolates interpolation and snapshot-EMF "
            "errors; moving-frame accuracy is covered by manufactured tensor tests"
        ),
    }
    report["assessment"] = assess_campaign(report)
    return report


def _positive_ints(values: list[int], name: str) -> list[int]:
    if not values or any(value < 1 for value in values):
        raise ValueError(f"{name} must contain positive integers")
    if len(values) != len(set(values)):
        raise ValueError(f"{name} contains duplicate values")
    return values


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--athena", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--outer-cells", type=int, nargs="+", default=(16, 24, 32))
    parser.add_argument("--cadence-strides", type=int, nargs="+", default=(1, 2, 4))
    parser.add_argument("--quadrature-orders", type=int, nargs="+", default=(1, 2, 3))
    parser.add_argument("--replay-quadrature", type=int, default=2)
    parser.add_argument("--skip-replay", action="store_true")
    parser.add_argument("--fail-on-gate", action="store_true")
    parser.add_argument("--half-width", type=float, default=4.0)
    parser.add_argument("--tlim", type=float, default=0.08)
    parser.add_argument("--cfl", type=float, default=0.02)
    parser.add_argument("--gamma", type=float, default=4.0 / 3.0)
    arguments = parser.parse_args()
    arguments.outer_cells = _positive_ints(arguments.outer_cells, "outer_cells")
    if any(value < 8 or value % 2 != 0 for value in arguments.outer_cells):
        parser.error("--outer-cells entries must be even and at least eight")
    arguments.cadence_strides = _positive_ints(
        arguments.cadence_strides, "cadence_strides"
    )
    arguments.quadrature_orders = _positive_ints(
        arguments.quadrature_orders, "quadrature_orders"
    )
    if any(value > 8 for value in arguments.quadrature_orders):
        parser.error("--quadrature-orders cannot exceed eight")
    if arguments.replay_quadrature not in arguments.quadrature_orders:
        parser.error("--replay-quadrature must be one of --quadrature-orders")
    if any(
        not math.isfinite(value) or value <= 0.0
        for value in (arguments.half_width, arguments.tlim, arguments.cfl)
    ):
        parser.error("--half-width, --tlim, and --cfl must be finite and positive")
    if not math.isfinite(arguments.gamma) or arguments.gamma <= 1.0:
        parser.error("--gamma must be finite and greater than one")
    return arguments


def main() -> int:
    arguments = parse_arguments()
    report = run_campaign(arguments)
    output = arguments.workdir.expanduser().resolve() / "summary.json"
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for name, case in report["cases"].items():
        comparison = case["comparison"]
        print(
            f"{name}: state={comparison['state']['relative_l2']:.6e}, "
            f"flux={comparison['normal_flux']['relative_l2']:.6e}, "
            f"emf={comparison['emf']['relative_l2']:.6e}, "
            f"edge_correction={case['projection']['maximum_edge_correction']:.6e}"
        )
    gates = report["assessment"]["smoke_gates"]
    print(f"smoke_gates={'PASS' if gates['passed'] else 'FAIL'}")
    print(output)
    return 1 if arguments.fail_on_gate and not gates["passed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
