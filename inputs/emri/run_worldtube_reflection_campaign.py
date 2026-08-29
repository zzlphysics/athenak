#!/usr/bin/env python3
"""Run reproducible seven-wave SRMHD outer-to-inner worldtube closure cases."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess
import time
from typing import Any

import numpy as np

import analyze_worldtube_reflection as reflection
import worldtube_flux_emf as worldtube


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "inputs" / "emri" / "worldtube_linear_wave.athinput"


def _tag(value: float) -> str:
    return format(value, ".8g").replace("-", "m").replace(".", "p")


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
        tail = log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-30:]
        raise RuntimeError(
            f"AthenaK failed with exit code {completed.returncode}:\n"
            + "\n".join(tail)
        )
    return elapsed


def _final_output(directory: Path, stem: str) -> Path:
    candidates = sorted((directory / "bin").glob(f"{stem}.mhd_w_bcc.*.bin"))
    if not candidates:
        raise RuntimeError(f"no final mhd_w_bcc output found below {directory}")
    return candidates[-1]


def _common_arguments(
    mode: int,
    amplitude: float,
    cfl: float,
    mesh_cells: tuple[int, int, int],
    meshblock_cells: tuple[int, int, int],
) -> list[str]:
    nx1, nx2, nx3 = mesh_cells
    bx1, bx2, bx3 = meshblock_cells
    return [
        f"mesh/nx1={nx1}",
        f"mesh/nx2={nx2}",
        f"mesh/nx3={nx3}",
        f"meshblock/nx1={bx1}",
        f"meshblock/nx2={bx2}",
        f"meshblock/nx3={bx3}",
        f"time/cfl_number={cfl:.17g}",
        f"problem/wave_flag={mode}",
        f"problem/amp={amplitude:.17g}",
        "output1/dt=100.0",
        "output2/dt=100.0",
    ]


def run_case(
    athena: Path,
    input_file: Path,
    campaign_directory: Path,
    mode: int,
    inner_cells: int,
    amplitude: float,
    cfl: float,
    dcycle: int,
    adiabatic_index: float,
    fluid_boundary: str,
) -> dict[str, Any]:
    """Run one mode and return its JSON-serializable reflection report."""

    case = campaign_directory / f"mode{mode}"
    if case.exists():
        raise FileExistsError(f"refusing to reuse existing case directory {case}")
    outer_directory = case / "outer"
    inner_directory = case / "inner"
    outer_directory.mkdir(parents=True)
    inner_directory.mkdir()
    outer_cells = (4 * inner_cells, 2 * inner_cells, 2 * inner_cells)
    block = min(inner_cells, 16)
    if inner_cells % block != 0:
        raise ValueError("inner_cells must be divisible by min(inner_cells,16)")

    tube_basename = case / "tube"
    outer_command = [
        str(athena),
        "-i",
        str(input_file),
        "-d",
        str(outer_directory),
        "job/basename=outer",
        *_common_arguments(
            mode, amplitude, cfl, outer_cells, (block, block, block)
        ),
        "emri_worldtube/enabled=true",
        "emri_worldtube/mode=outer",
        "emri_worldtube/overwrite=true",
        f"emri_worldtube/dcycle={dcycle}",
        f"emri_worldtube/file_basename={tube_basename}",
    ]
    outer_seconds = _run(outer_command, ROOT, case / "outer.log")
    manifests = sorted(case.glob("tube.cycle*.manifest.json"))
    if len(manifests) != 1:
        raise RuntimeError(f"expected one complete outer manifest, found {len(manifests)}")
    times, faces, metadata = worldtube.read_outer_stream(manifests[0])
    packed = case / "tube.npz"
    worldtube.write_worldtube(packed, times, faces, metadata)
    replay = case / "inner.bin"
    transfer_diagnostics = worldtube.write_inner_binary(
        replay, times, faces, metadata
    )

    inner_command = [
        str(athena),
        "-i",
        str(input_file),
        "-d",
        str(inner_directory),
        "job/basename=inner",
        *_common_arguments(
            mode,
            amplitude,
            cfl,
            (inner_cells, inner_cells, inner_cells),
            (block, block, block),
        ),
        "mesh/x1min=-0.25",
        "mesh/x1max=0.25",
        "mesh/x2min=-0.25",
        "mesh/x2max=0.25",
        "mesh/x3min=-0.25",
        "mesh/x3max=0.25",
        "mesh/ix1_bc=outflow",
        "mesh/ox1_bc=outflow",
        "mesh/ix2_bc=outflow",
        "mesh/ox2_bc=outflow",
        "mesh/ix3_bc=outflow",
        "mesh/ox3_bc=outflow",
        "emri_worldtube/enabled=true",
        "emri_worldtube/mode=inner",
        f"emri_worldtube/file={replay}",
        f"emri_worldtube/fluid_boundary={fluid_boundary}",
    ]
    inner_seconds = _run(inner_command, ROOT, case / "inner.log")
    reference_path = _final_output(outer_directory, "outer")
    candidate_path = _final_output(inner_directory, "inner")
    report = reflection.analyze_files(
        reference_path,
        candidate_path,
        mode,
        adiabatic_index,
    )
    report["campaign"] = {
        "inner_cells_per_axis": inner_cells,
        "outer_cells": list(outer_cells),
        "amplitude": amplitude,
        "cfl_number": cfl,
        "worldtube_dcycle": dcycle,
        "fluid_boundary": fluid_boundary,
        "worldtube_samples": int(times.size),
        "maximum_worldtube_dt": float(np.max(np.diff(times))),
        "outer_wall_seconds": outer_seconds,
        "inner_wall_seconds": inner_seconds,
        "reference_output": str(reference_path),
        "candidate_output": str(candidate_path),
        "transfer_maximum_closed_surface_flux": transfer_diagnostics[
            "maximum_closed_surface_flux"
        ],
        "transfer_maximum_faraday_residual": max(
            transfer_diagnostics["maximum_faraday_residual_by_face"].values()
        ),
    }
    (case / "reflection.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--athena", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--modes", type=int, nargs="+", default=list(range(7)))
    parser.add_argument("--inner-cells", type=int, default=8)
    parser.add_argument("--amplitude", type=float, default=1.0e-3)
    parser.add_argument("--cfl", type=float, default=0.3)
    parser.add_argument("--dcycle", type=int, default=1)
    parser.add_argument("--gamma", type=float, default=4.0 / 3.0)
    parser.add_argument(
        "--fluid-boundary",
        choices=("riemann", "characteristic_sr"),
        default="riemann",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    athena = arguments.athena.resolve()
    input_file = arguments.input.resolve()
    if not athena.is_file() or not input_file.is_file():
        raise FileNotFoundError("--athena and --input must name existing files")
    if arguments.inner_cells < 4 or arguments.dcycle < 1:
        raise ValueError("inner-cells must be >=4 and dcycle must be >=1")
    if (
        not math.isfinite(arguments.amplitude)
        or arguments.amplitude <= 0.0
        or not math.isfinite(arguments.cfl)
        or arguments.cfl <= 0.0
    ):
        raise ValueError("amplitude and cfl must be finite and positive")
    modes = list(dict.fromkeys(arguments.modes))
    if any(mode < 0 or mode > 6 for mode in modes):
        raise ValueError("all modes must lie in [0,6]")
    campaign = arguments.workdir.resolve() / (
        f"n{arguments.inner_cells}_cfl{_tag(arguments.cfl)}_"
        f"dc{arguments.dcycle}_amp{_tag(arguments.amplitude)}_"
        f"fb{arguments.fluid_boundary}"
    )
    if campaign.exists():
        raise FileExistsError(f"refusing to reuse existing campaign {campaign}")
    campaign.mkdir(parents=True)
    reports: list[dict[str, Any]] = []
    for mode in modes:
        report = run_case(
            athena,
            input_file,
            campaign,
            mode,
            arguments.inner_cells,
            arguments.amplitude,
            arguments.cfl,
            arguments.dcycle,
            arguments.gamma,
            arguments.fluid_boundary,
        )
        reports.append(report)
        print(
            f"mode {mode}: R_amp={report['reflected_amplitude_coefficient']:.6e}, "
            f"demeaned={report['demeaned_reflected_amplitude_coefficient']:.6e}",
            flush=True,
        )
    summary = {
        "classification": "athenak-emri-worldtube-reflection-campaign-v1",
        "athena": str(athena),
        "input": str(input_file),
        "modes": modes,
        "cases": reports,
    }
    (campaign / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(campaign / "summary.json")


if __name__ == "__main__":
    main()
