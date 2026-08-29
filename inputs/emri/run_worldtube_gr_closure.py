#!/usr/bin/env python3
"""Run a reproducible curved-spacetime Riemann/characteristic worldtube closure test."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import time
from typing import Any

import compare_worldtube_closure as closure
import worldtube_flux_emf as worldtube


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "inputs" / "emri" / "emri_windtunnel_smoke.athinput"


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
        tail = log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-40:]
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


def _evolution_arguments(tlim: float, cfl: float, root_dt_max: float) -> list[str]:
    return [
        "time/integrator=rk3",
        "time/nlim=-1",
        f"time/tlim={tlim:.17g}",
        f"time/cfl_number={cfl:.17g}",
        f"time/root_dt_max={root_dt_max:.17g}",
        "problem/user_hist=false",
        "output1/variable=mhd_w_bcc",
        f"output1/dt={tlim:.17g}",
        "output2/dt=100",
        "output3/dt=0",
        "output4/dt=0",
    ]


def _projection_counts(log_path: Path) -> tuple[int, int] | None:
    match = re.search(
        r"characteristic GR boundary: projections=(\d+), fallbacks=(\d+)",
        log_path.read_text(encoding="utf-8", errors="replace"),
    )
    if match is None:
        return None
    return int(match.group(1)), int(match.group(2))


def run_campaign(arguments: argparse.Namespace) -> dict[str, Any]:
    athena = arguments.athena.resolve()
    input_file = arguments.input.resolve()
    workdir = arguments.workdir.resolve()
    if workdir.exists():
        raise FileExistsError(f"refusing to reuse existing campaign directory {workdir}")
    if not athena.is_file() or not input_file.is_file():
        raise FileNotFoundError("--athena and --input must name existing files")
    if arguments.inner_cells < 4 or arguments.inner_cells % 2 != 0:
        raise ValueError("inner-cells must be even and at least four")
    if arguments.dcycle < 1:
        raise ValueError("dcycle must be positive")
    if any(
        not math.isfinite(value) or value <= 0.0
        for value in (
            arguments.half_width,
            arguments.tlim,
            arguments.cfl,
            arguments.root_dt_max,
        )
    ):
        raise ValueError("half-width, tlim, cfl, and root-dt-max must be positive")
    if not math.isfinite(arguments.gamma) or arguments.gamma <= 1.0:
        raise ValueError("gamma must be finite and greater than one")
    workdir.mkdir(parents=True)
    outer_directory = workdir / "outer"
    outer_directory.mkdir()
    outer_cells = 2 * arguments.inner_cells
    outer_half_width = 2.0 * arguments.half_width
    outer_block = min(outer_cells, 16)
    inner_block = arguments.inner_cells // 2
    tube_basename = workdir / "tube"
    common_evolution = _evolution_arguments(
        arguments.tlim, arguments.cfl, arguments.root_dt_max
    )
    outer_command = [
        str(athena),
        "-i",
        str(input_file),
        "-d",
        str(outer_directory),
        "job/basename=outer",
        *_mesh_arguments(outer_cells, outer_half_width, outer_block),
        *common_evolution,
        "emri_worldtube/enabled=true",
        "emri_worldtube/mode=outer",
        "emri_worldtube/overwrite=true",
        f"emri_worldtube/half_width={arguments.half_width:.17g}",
        f"emri_worldtube/dcycle={arguments.dcycle}",
        f"emri_worldtube/file_basename={tube_basename}",
    ]
    outer_seconds = _run(outer_command, ROOT, workdir / "outer.log")
    manifests = sorted(workdir.glob("tube.cycle*.manifest.json"))
    if len(manifests) != 1:
        raise RuntimeError(f"expected one complete outer manifest, found {len(manifests)}")
    times, faces, metadata = worldtube.read_outer_stream(manifests[0])
    packed = workdir / "tube.npz"
    worldtube.write_worldtube(packed, times, faces, metadata)
    replay = workdir / "inner.bin"
    transfer = worldtube.write_inner_binary(replay, times, faces, metadata)
    reference = _final_output(outer_directory, "outer")

    cases: dict[str, Any] = {}
    force_radii = (0.55, 0.65, 0.75, 0.875)
    for boundary in arguments.fluid_boundaries:
        directory = workdir / f"inner_{boundary}"
        directory.mkdir()
        log_path = workdir / f"inner_{boundary}.log"
        inner_command = [
            str(athena),
            "-i",
            str(input_file),
            "-d",
            str(directory),
            "job/basename=inner",
            *_mesh_arguments(
                arguments.inner_cells, arguments.half_width, inner_block
            ),
            *common_evolution,
            f"problem/force_surface_radius={force_radii[0]*arguments.half_width:.17g}",
            f"problem/force_outer_radius_1={force_radii[1]*arguments.half_width:.17g}",
            f"problem/force_outer_radius_2={force_radii[2]*arguments.half_width:.17g}",
            f"problem/force_outer_radius_3={force_radii[3]*arguments.half_width:.17g}",
            "emri_worldtube/enabled=true",
            "emri_worldtube/mode=inner",
            f"emri_worldtube/file={replay}",
            f"emri_worldtube/fluid_boundary={boundary}",
            "emri_worldtube/fluid_state_frame=inner_coordinate",
        ]
        inner_seconds = _run(inner_command, ROOT, log_path)
        candidate = _final_output(directory, "inner")
        comparison = closure.compare_files(reference, candidate, arguments.gamma)
        counts = _projection_counts(log_path)
        cases[boundary] = {
            "wall_seconds": inner_seconds,
            "candidate_output": str(candidate),
            "projection_count": None if counts is None else counts[0],
            "fallback_count": None if counts is None else counts[1],
            "closure": comparison,
        }

    return {
        "classification": "athenak-emri-worldtube-gr-closure-campaign-v1",
        "athena": str(athena),
        "input": str(input_file),
        "reference_output": str(reference),
        "outer_wall_seconds": outer_seconds,
        "inner_cells_per_axis": arguments.inner_cells,
        "outer_cells_per_axis": outer_cells,
        "half_width": arguments.half_width,
        "tlim": arguments.tlim,
        "cfl_number": arguments.cfl,
        "root_dt_max": arguments.root_dt_max,
        "worldtube_dcycle": arguments.dcycle,
        "worldtube_samples": int(times.size),
        "transfer_maximum_closed_surface_flux": transfer[
            "maximum_closed_surface_flux"
        ],
        "transfer_maximum_faraday_residual": max(
            transfer["maximum_faraday_residual_by_face"].values()
        ),
        "cases": cases,
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--athena", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--inner-cells", type=int, default=8)
    parser.add_argument("--half-width", type=float, default=4.0)
    parser.add_argument("--tlim", type=float, default=0.02)
    parser.add_argument("--cfl", type=float, default=0.15)
    parser.add_argument("--root-dt-max", type=float, default=0.004)
    parser.add_argument("--dcycle", type=int, default=1)
    parser.add_argument("--gamma", type=float, default=4.0 / 3.0)
    parser.add_argument(
        "--fluid-boundaries",
        nargs="+",
        choices=("riemann", "characteristic_gr"),
        default=("riemann", "characteristic_gr"),
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    report = run_campaign(arguments)
    output = arguments.workdir.resolve() / "summary.json"
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for boundary, case in report["cases"].items():
        comparison = case["closure"]
        print(
            f"{boundary}: density L2={comparison['variables']['dens']['relative_l2']:.6e}, "
            f"pressure L2={comparison['variables']['press']['relative_l2']:.6e}, "
            f"velocity L2={comparison['vector_groups']['velocity']['relative_l2']:.6e}"
        )
    print(output)


if __name__ == "__main__":
    main()
