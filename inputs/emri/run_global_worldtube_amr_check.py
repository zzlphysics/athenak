#!/usr/bin/env python3
"""Validate fixed-leaf-level extraction against an equivalent uniform source."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time
from typing import Any

import numpy as np

import analyze_global_worldtube_convergence as analysis
import extract_global_worldtube as extract
import run_global_worldtube_convergence as campaign
import worldtube_flux_emf as worldtube


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "inputs" / "emri" / "emri_windtunnel_smoke.athinput"
CLASSIFICATION = "athenak-emri-global-worldtube-amr-check-v1"


def _source_command(
    arguments: argparse.Namespace,
    directory: Path,
    cells: int,
    meshblock: int,
    refinement: str,
) -> list[str]:
    return [
        str(arguments.athena),
        "-i",
        str(arguments.input),
        "-d",
        str(directory),
        "job/basename=outer",
        *campaign._mesh_arguments(cells, arguments.source_half_width, meshblock),
        *campaign._evolution_arguments(arguments.tlim, arguments.cfl, True),
        f"mesh_refinement/refinement={refinement}",
        "emri_worldtube/enabled=false",
    ]


def _topology(path: Path) -> dict[str, object]:
    data = extract.bin_convert.read_binary(str(path), variables=())
    levels, counts = np.unique(data["mb_logical"][:, 3], return_counts=True)
    return {
        "leaf_meshblock_count": int(data["n_mbs"]),
        "leaf_level_counts": {
            str(int(level)): int(count)
            for level, count in zip(levels, counts, strict=True)
        },
    }


def _write_manifest(
    path: Path,
    snapshots: list[dict[str, object]],
    arguments: argparse.Namespace,
    source_level: int | None,
) -> None:
    times = [float(entry["time"]) for entry in snapshots]
    document: dict[str, Any] = {
        "classification": extract.INPUT_CLASSIFICATION,
        "snapshots": [
            {"state": str(entry["state"]), "adm": str(entry["adm"])}
            for entry in snapshots
        ],
        "frame": campaign.identity_frame_document(times),
        "sample_times": times,
        "worldtube": {
            "center": [0.0, 0.0, 0.0],
            "half_width": arguments.half_width,
            "cells_per_edge": arguments.face_cells,
        },
        "quadrature_order": arguments.quadrature_order,
        "hash_source_files": False,
        "global_length_in_local_units": 1.0,
        "density_renormalization": 1.0,
    }
    if source_level is not None:
        document["source_level"] = source_level
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def _extract(
    manifest: Path, output: Path
) -> tuple[
    np.ndarray,
    dict[str, worldtube.FaceData],
    dict[str, object],
    float,
]:
    start = time.perf_counter()
    times, faces, metadata, _ = extract.extract_manifest(manifest)
    elapsed = time.perf_counter() - start
    worldtube.write_worldtube(output, times, faces, metadata)
    return times, faces, metadata, elapsed


def _maximum(condition: dict[str, object], observed: float, threshold: float) -> None:
    condition.update(
        {
            "relation": "maximum",
            "observed": observed,
            "threshold": threshold,
            "passed": observed <= threshold,
        }
    )


def _assessment(
    comparison: dict[str, object],
    projection: dict[str, float],
    replay: dict[str, object],
) -> dict[str, object]:
    conditions: dict[str, dict[str, object]] = {}
    thresholds = {
        "state_relative_l2": 1.0e-6,
        "normal_flux_relative_l2": 1.0e-6,
        "emf_relative_l2": 1.0e-3,
        "relative_edge_correction": 1.0e-2,
        "boundary_flux_residual": 1.0e-12,
        "fallback_fraction": 1.0e-6,
        "density_closure": 1.0e-5,
        "pressure_closure": 1.0e-5,
        "velocity_closure": 1.0e-5,
        "magnetic_closure": 1.0e-5,
    }
    for field in ("state", "normal_flux", "emf"):
        _maximum(
            conditions.setdefault(f"{field}_relative_l2", {}),
            float(comparison[field]["relative_l2"]),
            thresholds[f"{field}_relative_l2"],
        )
    _maximum(
        conditions.setdefault("relative_edge_correction", {}),
        projection["maximum_relative_edge_correction"],
        thresholds["relative_edge_correction"],
    )
    projections = int(replay["projection_count"])
    fallbacks = int(replay["fallback_count"])
    _maximum(
        conditions.setdefault("fallback_fraction", {}),
        fallbacks / projections,
        thresholds["fallback_fraction"],
    )
    _maximum(
        conditions.setdefault("boundary_flux_residual", {}),
        float(replay["maximum_boundary_flux_residual"]),
        thresholds["boundary_flux_residual"],
    )
    closure = replay["closure"]
    closure_paths = {
        "density_closure": ("variables", "dens"),
        "pressure_closure": ("variables", "press"),
        "velocity_closure": ("vector_groups", "velocity"),
        "magnetic_closure": ("vector_groups", "magnetic"),
    }
    for name, path in closure_paths.items():
        _maximum(
            conditions.setdefault(name, {}),
            float(closure[path[0]][path[1]]["relative_l2"]),
            thresholds[name],
        )
    return {
        "passed": all(bool(value["passed"]) for value in conditions.values()),
        "conditions": conditions,
    }


def run_check(arguments: argparse.Namespace) -> dict[str, object]:
    arguments.athena = arguments.athena.expanduser().resolve(strict=True)
    arguments.input = arguments.input.expanduser().resolve(strict=True)
    workdir = arguments.workdir.expanduser().resolve()
    if workdir.exists():
        raise FileExistsError(f"refusing to reuse existing check directory {workdir}")
    workdir.mkdir(parents=True)
    fine_cells = arguments.root_cells * (1 << arguments.source_level)
    source_spacing = 2.0 * arguments.source_half_width / fine_cells
    target_spacing = 2.0 * arguments.half_width / arguments.face_cells
    if not math.isclose(source_spacing, target_spacing, rel_tol=0.0, abs_tol=1.0e-14):
        raise ValueError("target spacing must equal the selected AMR source spacing")

    source_results = {}
    source_snapshots = {}
    configurations = {
        "amr": (
            arguments.root_cells,
            arguments.amr_meshblock,
            "static",
        ),
        "uniform": (
            fine_cells,
            arguments.uniform_meshblock,
            "none",
        ),
    }
    for name, (cells, meshblock, refinement) in configurations.items():
        directory = workdir / f"source_{name}"
        directory.mkdir()
        wall = campaign._run(
            _source_command(
                arguments, directory, cells, meshblock, refinement
            ),
            ROOT,
            workdir / f"source_{name}.log",
        )
        snapshots = campaign._paired_snapshots(directory)
        final_output = campaign._final_output(directory, "outer")
        source_snapshots[name] = snapshots
        source_results[name] = {
            "wall_seconds": wall,
            "final_output": str(final_output),
            "topology": _topology(final_output),
        }

    products = {}
    loaded = {}
    for name in ("amr", "uniform"):
        directory = workdir / f"extract_{name}"
        directory.mkdir()
        manifest = directory / "manifest.json"
        output = directory / "worldtube.npz"
        _write_manifest(
            manifest,
            source_snapshots[name],
            arguments,
            arguments.source_level if name == "amr" else None,
        )
        times, faces, metadata, wall = _extract(manifest, output)
        loaded[name] = (times, faces, metadata)
        products[name] = {
            "manifest": str(manifest),
            "worldtube": str(output),
            "wall_seconds": wall,
            "source_level": metadata["source_level"],
            "snapshot_loading": metadata["sampling_diagnostics"][
                "snapshot_loading"
            ],
            "projection": campaign._projection_summary(metadata),
        }

    uniform_times, uniform_faces, _ = loaded["uniform"]
    amr_times, amr_faces, amr_metadata = loaded["amr"]
    comparison = analysis.compare_worldtubes(
        uniform_times,
        uniform_faces,
        amr_times,
        amr_faces,
        state_variables=amr_metadata["state_variables"],
    )
    replay_path = workdir / "extract_amr" / "inner.bin"
    worldtube.write_inner_binary(
        replay_path, amr_times, amr_faces, amr_metadata
    )
    replay = campaign._run_inner_replay(
        arguments,
        workdir / "inner_replay",
        replay_path,
        arguments.face_cells,
        Path(source_results["uniform"]["final_output"]),
    )
    assessment = _assessment(
        comparison, products["amr"]["projection"], replay
    )
    amr_levels = source_results["amr"]["topology"]["leaf_level_counts"]
    multilevel = len(amr_levels) >= 2 and str(arguments.source_level) in amr_levels
    assessment["conditions"]["multilevel_source_topology"] = {
        "relation": "required",
        "observed": amr_levels,
        "required_source_level": arguments.source_level,
        "passed": multilevel,
    }
    assessment["passed"] = bool(assessment["passed"] and multilevel)
    return {
        "classification": CLASSIFICATION,
        "athena": str(arguments.athena),
        "input": str(arguments.input),
        "source_level": arguments.source_level,
        "source_spacing": source_spacing,
        "source_runs": source_results,
        "products": products,
        "comparison": comparison,
        "inner_replay": replay,
        "assessment": assessment,
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--athena", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--source-half-width", type=float, default=8.0)
    parser.add_argument("--root-cells", type=int, default=16)
    parser.add_argument("--source-level", type=int, default=1)
    parser.add_argument("--amr-meshblock", type=int, default=4)
    parser.add_argument("--uniform-meshblock", type=int, default=8)
    parser.add_argument("--half-width", type=float, default=3.0)
    parser.add_argument("--face-cells", type=int, default=12)
    parser.add_argument("--quadrature-order", type=int, default=2)
    parser.add_argument("--tlim", type=float, default=1.0e-4)
    parser.add_argument("--cfl", type=float, default=0.02)
    parser.add_argument("--gamma", type=float, default=4.0 / 3.0)
    parser.add_argument("--fail-on-gate", action="store_true")
    arguments = parser.parse_args()
    integer_values = (
        arguments.root_cells,
        arguments.amr_meshblock,
        arguments.uniform_meshblock,
        arguments.face_cells,
        arguments.quadrature_order,
    )
    if any(value < 1 for value in integer_values) or arguments.source_level < 0:
        parser.error("cell counts and quadrature must be positive; level nonnegative")
    real_values = (
        arguments.source_half_width,
        arguments.half_width,
        arguments.tlim,
        arguments.cfl,
    )
    if any(not math.isfinite(value) or value <= 0.0 for value in real_values):
        parser.error("lengths, tlim, and cfl must be finite and positive")
    return arguments


def main() -> int:
    arguments = parse_arguments()
    report = run_check(arguments)
    output = arguments.workdir.expanduser().resolve() / "summary.json"
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    comparison = report["comparison"]
    print(
        f"state={comparison['state']['relative_l2']:.6e}, "
        f"flux={comparison['normal_flux']['relative_l2']:.6e}, "
        f"emf={comparison['emf']['relative_l2']:.6e}"
    )
    passed = bool(report["assessment"]["passed"])
    print(f"amr_gates={'PASS' if passed else 'FAIL'}")
    print(output)
    return 1 if arguments.fail_on_gate and not passed else 0


if __name__ == "__main__":
    raise SystemExit(main())
