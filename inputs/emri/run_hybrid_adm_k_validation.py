#!/usr/bin/env python3
"""Validate hybrid ADM-v3 ``K_ij`` reconstruction with AthenaK.

The validation is independent of a global GRMHD snapshot.  It first compares a
Minkowski-primary v3 replay against the existing isolated analytic-secondary
metric at identical finite-difference steps.  It then replays a time-dependent
affine primary metric under dynamic AMR and compares every ADM output cell with
the host-side reconstruction.  The latter exercises primary spatial
derivatives, slab time derivatives, a time-dependent secondary coframe, online
secondary derivatives, and direct evaluation on refined MeshBlocks.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import adm_volume_replay as adm_volume
import extract_global_worldtube as extract
import run_adm_inner_replay_pilot as pilot
import worldtube_flux_emf as worldtube


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "inputs" / "emri" / "emri_windtunnel_smoke.athinput"
CLASSIFICATION = "athenak-emri-hybrid-adm-k-validation-v1"
K_FIELDS = tuple(name for name in pilot.ATHENA_ADM_FIELDS if "_K" in name)
NON_K_FIELDS = tuple(name for name in pilot.ATHENA_ADM_FIELDS if name not in K_FIELDS)


def _positive(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def _constant_worldtube(
    path: Path, times: np.ndarray, cells: int, half_width: float
) -> dict[str, object]:
    names = ("rho", "u1", "u2", "u3", "pgas", "bcc1", "bcc2", "bcc3")
    state = np.zeros((times.size, len(names), cells, cells), dtype=np.float64)
    state[:, 0] = 1.0
    state[:, 1] = 0.15
    state[:, 4] = 1.0e-2
    faces = {}
    for name in worldtube.FACE_NAMES:
        faces[name] = worldtube.FaceData(
            cell_state=state.copy(),
            normal_flux=np.zeros((times.size, cells, cells)),
            emf_u=np.zeros((times.size - 1, cells + 1, cells)),
            emf_v=np.zeros((times.size - 1, cells, cells + 1)),
        )
    metadata = {
        "classification": CLASSIFICATION,
        "construction": "constant unmagnetized manufactured boundary",
        "state_variables": list(names),
        "center": [0.0, 0.0, 0.0],
        "half_width": half_width,
    }
    return worldtube.write_inner_binary(path, times, faces, metadata)


def _primary_metric(
    time: float, x: np.ndarray, y: np.ndarray, z: np.ndarray, affine: bool
) -> tuple[np.ndarray, np.ndarray]:
    shape = x.shape
    metric = np.zeros(shape + (4, 4), dtype=np.float64)
    derivative = np.zeros(shape + (3, 4, 4), dtype=np.float64)
    metric[..., 0, 0] = -1.0
    metric[..., 1, 1] = 1.0
    metric[..., 2, 2] = 1.0
    metric[..., 3, 3] = 1.0
    if not affine:
        return metric, derivative

    # A deliberately non-diagonal, time-dependent, affine Lorentzian metric.
    # Trilinear replay and the stored spatial derivatives are exact for it.
    metric[..., 0, 0] = -1.25 + 0.006 * time + 0.010 * x - 0.004 * z
    metric[..., 0, 1] = 0.020 + 0.002 * time + 0.008 * x - 0.003 * y
    metric[..., 0, 2] = -0.015 + 0.001 * time + 0.004 * y
    metric[..., 0, 3] = 0.010 - 0.0015 * time + 0.003 * z
    metric[..., 1, 1] = 1.20 + 0.004 * time + 0.020 * x
    metric[..., 1, 2] = 0.010 + 0.003 * x - 0.002 * y
    metric[..., 1, 3] = -0.005 + 0.002 * z
    metric[..., 2, 2] = 1.10 - 0.003 * time - 0.015 * y
    metric[..., 2, 3] = 0.004 + 0.002 * x
    metric[..., 3, 3] = 1.05 + 0.002 * time + 0.010 * z
    for left in range(4):
        for right in range(left):
            metric[..., left, right] = metric[..., right, left]

    # direction, row, column, value
    coefficients = (
        (0, 0, 0, 0.010),
        (2, 0, 0, -0.004),
        (0, 0, 1, 0.008),
        (1, 0, 1, -0.003),
        (1, 0, 2, 0.004),
        (2, 0, 3, 0.003),
        (0, 1, 1, 0.020),
        (0, 1, 2, 0.003),
        (1, 1, 2, -0.002),
        (2, 1, 3, 0.002),
        (1, 2, 2, -0.015),
        (0, 2, 3, 0.002),
        (2, 3, 3, 0.010),
    )
    for direction, left, right, value in coefficients:
        derivative[..., direction, left, right] += value
        if left != right:
            derivative[..., direction, right, left] += value
    return metric, derivative


def manufactured_volume(
    times: np.ndarray,
    half_width: float,
    metric_cells: int,
    metric_halo: int,
    secondary_mass: float,
    secondary_chi: float,
    fd_ratio: float,
    affine_primary: bool,
) -> adm_volume.ADMVolume:
    spacing = 2.0 * half_width / metric_cells
    nodes = metric_cells + 2 * metric_halo
    lower = np.full(3, -half_width - (metric_halo - 0.5) * spacing)
    coordinates = lower[0] + spacing * np.arange(nodes)
    z, y, x = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    fields = np.empty(
        (times.size, len(adm_volume.HYBRID_FIELD_NAMES), nodes, nodes, nodes),
        dtype=np.float64,
    )
    coframes = np.empty((times.size, 4, 4), dtype=np.float64)
    minimum_lapse = math.inf
    minimum_eigenvalue = math.inf
    for time_index, time in enumerate(times):
        metric, derivative = _primary_metric(float(time), x, y, z, affine_primary)
        for field, (left, right) in enumerate(adm_volume.METRIC_COMPONENTS):
            fields[time_index, field] = metric[..., left, right]
        for direction in range(3):
            offset = len(adm_volume.METRIC_COMPONENTS) * (direction + 1)
            for field, (left, right) in enumerate(adm_volume.METRIC_COMPONENTS):
                fields[time_index, offset + field] = derivative[
                    ..., direction, left, right
                ]
        origin_metric, _ = _primary_metric(
            float(time), np.asarray(0.0), np.asarray(0.0), np.asarray(0.0), affine_primary
        )
        _, coframe = adm_volume._source_tetrad(origin_metric)
        coframes[time_index] = coframe
        positions = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
        composed = metric.reshape(-1, 4, 4) + adm_volume.secondary_kerr_perturbation(
            positions,
            secondary_mass,
            secondary_chi,
            coframe,
        )
        _, audit = adm_volume.decompose_four_metric(composed)
        minimum_lapse = min(minimum_lapse, float(np.min(audit["lapse"])))
        minimum_eigenvalue = min(
            minimum_eigenvalue,
            float(np.min(audit["minimum_spatial_eigenvalue"])),
        )
    parameters = np.asarray(
        (secondary_mass, secondary_chi, 0.05, 1.0e-3, fd_ratio * secondary_mass),
        dtype=np.float64,
    )
    return adm_volume.ADMVolume(
        times=np.asarray(times, dtype=np.float64),
        lower=lower,
        spacing=np.full(3, spacing),
        fields=fields,
        metadata={
            "classification": CLASSIFICATION,
            "primary": "affine" if affine_primary else "Minkowski",
            "minimum_composed_lapse": minimum_lapse,
            "minimum_composed_spatial_metric_eigenvalue": minimum_eigenvalue,
        },
        secondary_coframes=coframes,
        hybrid_parameters=parameters,
    )


def _output_path(directory: Path, basename: str, variable: str) -> Path:
    matches = sorted((directory / "bin").glob(f"{basename}.{variable}.*.bin"))
    if not matches:
        raise RuntimeError(f"AthenaK produced no {variable} output in {directory}")
    return matches[-1]


def _read_adm(path: Path) -> tuple[list[tuple[float, ...]], dict[str, np.ndarray]]:
    data = extract.bin_convert.read_binary(str(path))
    geometries = [tuple(float(value) for value in row) for row in data["mb_geometry"]]
    order = sorted(range(len(geometries)), key=geometries.__getitem__)
    fields = {
        name: np.concatenate(
            [np.asarray(data["mb_data"][name][index]).ravel() for index in order]
        ).astype(np.float64)
        for name in pilot.ATHENA_ADM_FIELDS
    }
    return [geometries[index] for index in order], fields


def _group_error(
    reference: dict[str, np.ndarray], candidate: dict[str, np.ndarray], names: tuple[str, ...]
) -> dict[str, float]:
    left = np.concatenate([reference[name] for name in names])
    right = np.concatenate([candidate[name] for name in names])
    difference = right - left
    tiny = np.finfo(float).tiny
    return {
        "relative_l2": float(np.linalg.norm(difference) / max(np.linalg.norm(left), tiny)),
        "relative_linf": float(
            np.max(np.abs(difference)) / max(np.max(np.abs(left)), tiny)
        ),
        "maximum_absolute": float(np.max(np.abs(difference))),
    }


def compare_adm_outputs(reference_path: Path, candidate_path: Path) -> dict[str, object]:
    reference_geometry, reference = _read_adm(reference_path)
    candidate_geometry, candidate = _read_adm(candidate_path)
    if reference_geometry != candidate_geometry:
        raise ValueError("ADM comparison outputs have different MeshBlock geometry")
    return {
        "reference": str(reference_path),
        "candidate": str(candidate_path),
        "meshblock_count": len(reference_geometry),
        "K_ij": _group_error(reference, candidate, K_FIELDS),
        "non_K_ADM": _group_error(reference, candidate, NON_K_FIELDS),
        "fields": {
            name: _group_error(reference, candidate, (name,))
            for name in pilot.ATHENA_ADM_FIELDS
        },
    }


def _base_command(
    arguments: argparse.Namespace,
    workdir: Path,
    basename: str,
    worldtube_path: Path,
    cells: int,
    meshblock_cells: int,
    half_width: float,
    duration: float,
    fd_step: float,
) -> list[str]:
    output = workdir / "run"
    output.mkdir()
    command = [
        str(arguments.athena.expanduser().resolve(strict=True)),
        "-i",
        str(arguments.input.expanduser().resolve(strict=True)),
        "-d",
        str(output),
        f"job/basename={basename}",
        *pilot._mesh_arguments(cells, half_width, meshblock_cells),
        f"mesh/nghost={arguments.mesh_nghost}",
        "mesh/ix1_bc=user",
        "mesh/ox1_bc=user",
        "mesh/ix2_bc=user",
        "mesh/ox2_bc=user",
        "mesh/ix3_bc=user",
        "mesh/ox3_bc=user",
        "time/integrator=rk3",
        f"time/cfl_number={arguments.cfl:.17g}",
        "time/nlim=-1",
        f"time/tlim={duration:.17g}",
        "problem/background_mode=isolated",
        "problem/secondary_embedding=tangent_tetrad",
        "problem/primary_mass=1",
        f"problem/secondary_mass={arguments.secondary_mass:.17g}",
        "problem/primary_chi=0",
        f"problem/secondary_chi={arguments.secondary_chi:.17g}",
        "problem/orbital_radius=10",
        "problem/require_stable_orbit=false",
        "problem/user_hist=false",
        "problem/wind_frame=normal_frame",
        "problem/rho0=1",
        "problem/pgas0=0.01",
        "problem/u1=0.15",
        "problem/u2=0",
        "problem/u3=0",
        "problem/b1=0",
        "problem/b2=0",
        "problem/b3=0",
        f"problem/metric_fd_step={fd_step:.17g}",
        "mhd/dfloor=1e-10",
        "mhd/tfloor=1e-10",
        "emri_worldtube/enabled=true",
        "emri_worldtube/mode=inner",
        f"emri_worldtube/file={worldtube_path}",
        # The metric test keeps CT time synchronization and topology rebuilding,
        # while analytic user ghost states isolate it from characteristic fallbacks.
        "emri_worldtube/fluid_boundary=off",
        "emri_worldtube/fluid_state_frame=inner_coordinate",
        "output1/variable=mhd_w_bcc",
        f"output1/dt={duration:.17g}",
        "output2/variable=mhd_divb",
        f"output2/dt={duration:.17g}",
        "output3/dt=0",
        "output4/dt=0",
        "output5/variable=adm",
        f"output5/dt={duration:.17g}",
    ]
    return command


def _run_case(command: list[str], workdir: Path) -> tuple[float, Path]:
    log = workdir / "athena.log"
    wall_seconds = pilot._run(command, log)
    return wall_seconds, log


def _parity_case(
    arguments: argparse.Namespace,
    parent: Path,
    fd_ratio: float,
) -> tuple[dict[str, object], dict[str, np.ndarray]]:
    label = f"fd_{fd_ratio:.8g}".replace(".", "p").replace("-", "m")
    case_dir = parent / label
    case_dir.mkdir(parents=True)
    times = np.asarray((0.0, arguments.duration), dtype=np.float64)
    worldtube_path = case_dir / "inner.bin"
    _constant_worldtube(
        worldtube_path, times, arguments.cells, arguments.half_width
    )
    volume = manufactured_volume(
        times,
        arguments.half_width,
        arguments.metric_cells,
        arguments.metric_halo,
        arguments.secondary_mass,
        arguments.secondary_chi,
        fd_ratio,
        affine_primary=False,
    )
    metric_path = case_dir / "adm.bin"
    adm_volume.write_binary(metric_path, volume)
    fd_step = fd_ratio * arguments.secondary_mass

    analytic_dir = case_dir / "analytic"
    analytic_dir.mkdir()
    analytic_command = _base_command(
        arguments,
        analytic_dir,
        "analytic",
        worldtube_path,
        arguments.cells,
        arguments.meshblock_cells,
        arguments.half_width,
        arguments.duration,
        fd_step,
    )
    analytic_wall, analytic_log = _run_case(analytic_command, analytic_dir)

    hybrid_dir = case_dir / "hybrid"
    hybrid_dir.mkdir()
    hybrid_command = _base_command(
        arguments,
        hybrid_dir,
        "hybrid",
        worldtube_path,
        arguments.cells,
        arguments.meshblock_cells,
        arguments.half_width,
        arguments.duration,
        fd_step,
    )
    hybrid_command.extend(
        (
            "adm/dynamic=true",
            "emri_adm_replay/enabled=true",
            f"emri_adm_replay/file={metric_path}",
        )
    )
    hybrid_wall, hybrid_log = _run_case(hybrid_command, hybrid_dir)
    analytic_adm = _output_path(analytic_dir / "run", "analytic", "adm")
    hybrid_adm = _output_path(hybrid_dir / "run", "hybrid", "adm")
    comparison = compare_adm_outputs(analytic_adm, hybrid_adm)
    _, hybrid_fields = _read_adm(hybrid_adm)
    report = {
        "fd_ratio_to_secondary_mass": fd_ratio,
        "fd_step": fd_step,
        "analytic_wall_seconds": analytic_wall,
        "hybrid_wall_seconds": hybrid_wall,
        "analytic_log": str(analytic_log),
        "hybrid_log": str(hybrid_log),
        "comparison": comparison,
    }
    return report, hybrid_fields


def _amr_case(arguments: argparse.Namespace, parent: Path) -> dict[str, object]:
    case_dir = parent / "affine_primary_amr"
    case_dir.mkdir()
    times = np.linspace(0.0, arguments.duration, 3)
    worldtube_path = case_dir / "inner.bin"
    _constant_worldtube(
        worldtube_path, times, arguments.amr_cells, arguments.half_width
    )
    middle_ratio = arguments.fd_ratios[len(arguments.fd_ratios) // 2]
    volume = manufactured_volume(
        times,
        arguments.half_width,
        arguments.metric_cells,
        arguments.metric_halo,
        arguments.secondary_mass,
        arguments.secondary_chi,
        middle_ratio,
        affine_primary=True,
    )
    metric_path = case_dir / "adm.bin"
    adm_volume.write_binary(metric_path, volume)
    run_input = pilot._amr_pilot_input(
        arguments.input.expanduser().resolve(strict=True), case_dir / "amr.athinput"
    )
    amr_arguments = argparse.Namespace(**vars(arguments))
    amr_arguments.input = run_input
    command = _base_command(
        amr_arguments,
        case_dir,
        "amr",
        worldtube_path,
        arguments.amr_cells,
        arguments.amr_meshblock_cells,
        arguments.half_width,
        arguments.duration,
        middle_ratio * arguments.secondary_mass,
    )
    maximum_blocks = (
        (arguments.amr_cells // arguments.amr_meshblock_cells) ** 3
        * 8 ** (arguments.amr_levels - 1)
    )
    command.extend(
        (
            "adm/dynamic=true",
            "emri_adm_replay/enabled=true",
            f"emri_adm_replay/file={metric_path}",
            "mesh_refinement/refinement=adaptive",
            f"mesh_refinement/num_levels={arguments.amr_levels}",
            f"mesh_refinement/max_nmb_per_rank={maximum_blocks}",
            "mesh_refinement/ncycle_check=1",
            "mesh_refinement/refinement_interval=1",
            "mesh_refinement/prolong_primitives=false",
            "amr_criterion0/method=user",
        )
    )
    for level in range(1, arguments.amr_levels):
        command.append(
            f"problem/refinement_radius_level_{level}="
            f"{arguments.refinement_radius / 2.0 ** (level - 1):.17g}"
        )
    wall_seconds, log = _run_case(command, case_dir)
    adm_paths = sorted((case_dir / "run" / "bin").glob("amr.adm.*.bin"))
    replay = pilot._adm_replay_comparison(
        adm_paths, volume, arguments.amr_cells, arguments.half_width
    )
    maximum_k_relative_linf = max(
        sample["fields"][name]["relative_linf"]
        for sample in replay["samples"]
        for name in K_FIELDS
    )
    divergence_path = _output_path(case_dir / "run", "amr", "mhd_divb")
    divergence = pilot._field_summary(divergence_path)["divb"]
    maximum_divb = max(abs(float(divergence["minimum"])), abs(float(divergence["maximum"])))
    log_diagnostics = pilot._log_diagnostics(log)
    return {
        "wall_seconds": wall_seconds,
        "log": str(log),
        "fd_ratio_to_secondary_mass": middle_ratio,
        "adm_replay_comparison": replay,
        "maximum_K_relative_linf": maximum_k_relative_linf,
        "maximum_divb": maximum_divb,
        "log_diagnostics": log_diagnostics,
    }


def run_validation(arguments: argparse.Namespace) -> dict[str, object]:
    workdir = arguments.workdir.expanduser().resolve()
    if workdir.exists():
        raise FileExistsError(f"refusing to reuse validation directory {workdir}")
    workdir.mkdir(parents=True)
    parity_reports = []
    hybrid_fields = []
    for fd_ratio in arguments.fd_ratios:
        report, fields = _parity_case(arguments, workdir / "parity", fd_ratio)
        parity_reports.append(report)
        hybrid_fields.append(fields)

    reference_index = len(arguments.fd_ratios) // 2
    step_sensitivity = []
    for index, fd_ratio in enumerate(arguments.fd_ratios):
        step_sensitivity.append(
            {
                "fd_ratio_to_secondary_mass": fd_ratio,
                "relative_to_middle": {
                    "K_ij": _group_error(
                        hybrid_fields[reference_index], hybrid_fields[index], K_FIELDS
                    ),
                    "non_K_ADM": _group_error(
                        hybrid_fields[reference_index], hybrid_fields[index], NON_K_FIELDS
                    ),
                },
            }
        )
    amr = None if arguments.skip_amr else _amr_case(arguments, workdir)

    conditions = {
        "analytic_hybrid_K_parity": {
            "relation": "maximum",
            "observed": max(
                item["comparison"]["K_ij"]["relative_l2"] for item in parity_reports
            ),
            "threshold": arguments.maximum_parity_relative_l2,
        },
        "fd_step_K_sensitivity": {
            "relation": "maximum",
            "observed": max(
                item["relative_to_middle"]["K_ij"]["relative_l2"]
                for item in step_sensitivity
            ),
            "threshold": arguments.maximum_fd_sensitivity,
        },
    }
    if amr is not None:
        conditions.update(
            {
                "amr_ADM_replay": {
                    "relation": "maximum",
                    "observed": amr["adm_replay_comparison"]["maximum_relative_linf"],
                    "threshold": arguments.maximum_adm_replay_relative_error,
                },
                "amr_K_replay": {
                    "relation": "maximum",
                    "observed": amr["maximum_K_relative_linf"],
                    "threshold": arguments.maximum_adm_replay_relative_error,
                },
                "amr_divB": {
                    "relation": "maximum",
                    "observed": amr["maximum_divb"],
                    "threshold": arguments.maximum_divb,
                },
                "amr_topology_rebuild": {
                    "relation": "minimum",
                    "observed": amr["log_diagnostics"][
                        "amr_boundary_topology_rebuild_count"
                    ],
                    "threshold": 1,
                },
            }
        )
    for condition in conditions.values():
        if condition["relation"] == "maximum":
            condition["passed"] = condition["observed"] <= condition["threshold"]
        else:
            condition["passed"] = condition["observed"] >= condition["threshold"]
    return {
        "classification": CLASSIFICATION,
        "athena": str(arguments.athena.expanduser().resolve()),
        "input": str(arguments.input.expanduser().resolve()),
        "secondary_mass": arguments.secondary_mass,
        "secondary_chi": arguments.secondary_chi,
        "fd_ratios": arguments.fd_ratios,
        "parity_cases": parity_reports,
        "fd_step_sensitivity": step_sensitivity,
        "affine_primary_amr": amr,
        "assessment": {
            "passed": all(condition["passed"] for condition in conditions.values()),
            "conditions": conditions,
        },
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--athena", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--secondary-mass", type=float, default=0.05)
    parser.add_argument("--secondary-chi", type=float, default=0.4)
    parser.add_argument(
        "--fd-ratios", type=float, nargs="+", default=(2.5e-5, 5.0e-5, 1.0e-4)
    )
    parser.add_argument("--half-width", type=float, default=0.5)
    parser.add_argument("--duration", type=float, default=2.0e-3)
    parser.add_argument("--cells", type=int, default=32)
    parser.add_argument("--meshblock-cells", type=int, default=8)
    parser.add_argument("--metric-cells", type=int, default=4)
    parser.add_argument("--metric-halo", type=int, default=4)
    parser.add_argument("--mesh-nghost", type=int, default=4)
    parser.add_argument("--cfl", type=float, default=0.02)
    parser.add_argument("--skip-amr", action="store_true")
    parser.add_argument("--amr-cells", type=int, default=32)
    parser.add_argument("--amr-meshblock-cells", type=int, default=8)
    parser.add_argument("--amr-levels", type=int, default=2)
    parser.add_argument("--refinement-radius", type=float, default=0.20)
    parser.add_argument("--maximum-parity-relative-l2", type=float, default=2.0e-6)
    parser.add_argument("--maximum-fd-sensitivity", type=float, default=1.0e-2)
    parser.add_argument(
        "--maximum-adm-replay-relative-error", type=float, default=1.0e-5
    )
    parser.add_argument("--maximum-divb", type=float, default=1.0e-10)
    parser.add_argument("--fail-on-gate", action="store_true")
    arguments = parser.parse_args()
    for name in (
        "secondary_mass",
        "half_width",
        "duration",
        "cfl",
        "refinement_radius",
        "maximum_parity_relative_l2",
        "maximum_fd_sensitivity",
        "maximum_adm_replay_relative_error",
        "maximum_divb",
    ):
        try:
            _positive(getattr(arguments, name), name)
        except ValueError as error:
            parser.error(str(error))
    if not math.isfinite(arguments.secondary_chi) or abs(arguments.secondary_chi) > 1.0:
        parser.error("secondary_chi must lie in [-1,1]")
    arguments.fd_ratios = sorted(float(value) for value in arguments.fd_ratios)
    if len(arguments.fd_ratios) < 3 or any(
        not math.isfinite(value) or value <= 0.0 for value in arguments.fd_ratios
    ):
        parser.error("fd-ratios requires at least three finite positive values")
    for cells_name, block_name in (
        ("cells", "meshblock_cells"),
        ("amr_cells", "amr_meshblock_cells"),
    ):
        cells = getattr(arguments, cells_name)
        block = getattr(arguments, block_name)
        if cells < 1 or block < 1 or cells % block:
            parser.error(f"{block_name} must be a positive divisor of {cells_name}")
    if arguments.metric_cells < 2 or arguments.metric_halo < 1 or arguments.mesh_nghost < 1:
        parser.error("metric_cells, metric_halo, and mesh_nghost are invalid")
    if arguments.amr_levels < 2 and not arguments.skip_amr:
        parser.error("AMR validation requires at least two levels")
    return arguments


def main() -> int:
    arguments = parse_arguments()
    report = run_validation(arguments)
    summary = arguments.workdir.expanduser().resolve() / "summary.json"
    summary.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(f"K_validation={'PASS' if report['assessment']['passed'] else 'FAIL'}")
    print(summary)
    return int(arguments.fail_on_gate and not report["assessment"]["passed"])


if __name__ == "__main__":
    raise SystemExit(main())
