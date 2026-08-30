#!/usr/bin/env python3
"""Run a fail-closed inner replay pilot from an ADM worldtube campaign.

The default structural mode retains the isolated analytic secondary metric.  With
``--numerical-adm-volume``, the pilot instead pulls back the numerical primary ADM
background and adds the analytic secondary Kerr perturbation in its tangent frame.
Force history remains disabled in that mode until its off-grid metric sampler uses
the same replayed metric.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import time

import numpy as np

import adm_volume_replay as adm_volume
import analyze_force_history as force_history
import extract_global_worldtube as extract
import run_adm_worldtube_campaign as campaign
import worldtube_flux_emf as worldtube


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "inputs" / "emri" / "emri_windtunnel_smoke.athinput"
CLASSIFICATION = "athenak-emri-adm-inner-replay-pilot-v2"
ATHENA_ADM_FIELDS = (
    "adm_gxx",
    "adm_gxy",
    "adm_gxz",
    "adm_gyy",
    "adm_gyz",
    "adm_gzz",
    "adm_Kxx",
    "adm_Kxy",
    "adm_Kxz",
    "adm_Kyy",
    "adm_Kyz",
    "adm_Kzz",
    "adm_psi4",
    "adm_alpha",
    "adm_betax",
    "adm_betay",
    "adm_betaz",
)


def _read_campaign(path: Path) -> tuple[Path, dict[str, object]]:
    resolved = path.expanduser().resolve(strict=True)
    document = json.loads(resolved.read_text(encoding="utf-8"))
    if document.get("classification") != campaign.CLASSIFICATION:
        raise ValueError("input is not an ADM worldtube campaign summary")
    return resolved, document


def _case_document(
    document: dict[str, object], requested: str | None
) -> tuple[str, dict[str, object]]:
    name = document["reference_case"] if requested is None else requested
    cases = document.get("cases")
    if not isinstance(cases, dict) or name not in cases:
        raise ValueError(f"campaign has no case named {name}")
    selected = cases[name]
    if not isinstance(selected, dict):
        raise ValueError("selected campaign case is not an object")
    return str(name), selected


def _face_positions(
    cells: int, center: np.ndarray, half_width: float, name: str
) -> np.ndarray:
    spacing = 2.0 * half_width / cells
    coordinates = -half_width + (np.arange(cells) + 0.5) * spacing
    u, v = np.meshgrid(coordinates, coordinates, indexing="xy")
    geometry = extract.CubeGeometry(center, half_width, cells)
    return geometry.face_positions(
        name, u, v, exterior_offset=0.5 * spacing
    ).reshape(-1, 3)


def fit_initial_affine_profile(
    faces: dict[str, worldtube.FaceData], metadata: dict[str, object]
) -> dict[str, object]:
    names = list(metadata.get("state_variables", ()))
    required = ("rho", "u1", "u2", "u3", "pgas", "bcc1", "bcc2", "bcc3")
    if any(name not in names for name in required):
        raise ValueError("worldtube does not contain the eight DynGRMHD boundary fields")
    center = np.asarray(metadata.get("center", ()), dtype=np.float64)
    half_width = float(metadata.get("half_width"))
    cells = faces[worldtube.FACE_NAMES[0]].cell_state.shape[-1]
    positions = np.concatenate(
        [
            _face_positions(cells, center, half_width, name)
            for name in worldtube.FACE_NAMES
        ]
    )
    state = np.concatenate(
        [
            faces[name].cell_state[0].reshape(len(names), -1).T
            for name in worldtube.FACE_NAMES
        ]
    )
    design = np.column_stack((np.ones(positions.shape[0]), positions))
    coefficients: dict[str, np.ndarray] = {}
    residuals = {}
    for name in required:
        values = state[:, names.index(name)]
        target = np.log(values) if name in ("rho", "pgas") else values
        fit, _, _, _ = np.linalg.lstsq(design, target, rcond=None)
        prediction = design @ fit
        coefficients[name] = fit
        scale = max(float(np.sqrt(np.mean(target**2))), np.finfo(float).tiny)
        residuals[name] = {
            "relative_l2": float(np.sqrt(np.mean((prediction - target) ** 2)))
            / scale,
            "maximum_absolute": float(np.max(np.abs(prediction - target))),
        }
    return {
        "rho0": float(math.exp(coefficients["rho"][0])),
        "pgas0": float(math.exp(coefficients["pgas"][0])),
        "u": [float(coefficients[f"u{axis}"][0]) for axis in (1, 2, 3)],
        "b": [float(coefficients[f"bcc{axis}"][0]) for axis in (1, 2, 3)],
        "log_density_gradient": coefficients["rho"][1:].tolist(),
        "log_pressure_gradient": coefficients["pgas"][1:].tolist(),
        "velocity_gradient": [
            coefficients[f"u{component}"][1:].tolist()
            for component in (1, 2, 3)
        ],
        "fit_residuals": residuals,
        "sample_count": int(positions.shape[0]),
    }


def _condition(
    relation: str, observed: float, threshold: float
) -> dict[str, object]:
    passed = observed >= threshold if relation == "minimum" else observed <= threshold
    return {
        "relation": relation,
        "observed": observed,
        "threshold": threshold,
        "passed": passed,
    }


def structural_assessment(
    arguments: argparse.Namespace,
    case: dict[str, object],
    cells: int,
    half_width: float,
) -> dict[str, object]:
    spin = arguments.secondary_chi
    horizon_radius = arguments.secondary_mass * (
        1.0 + math.sqrt(max(0.0, 1.0 - spin * spin))
    )
    spacing = 2.0 * half_width / cells
    boundary_state = case.get("boundary_state")
    if not isinstance(boundary_state, dict):
        raise ValueError("campaign case predates the required boundary-state audit")
    inner_validation = case.get("inner_validation")
    if not isinstance(inner_validation, dict):
        raise ValueError("campaign case has no inner-binary validation")
    initial_flux = inner_validation.get("initial_volume_flux")
    if not isinstance(initial_flux, dict):
        raise ValueError("campaign inner binary predates volume-flux initialization")
    conditions = {
        "secondary_horizon_resolution": _condition(
            "minimum",
            horizon_radius / spacing,
            arguments.minimum_horizon_cells,
        ),
        "worldtube_separation": _condition(
            "minimum",
            half_width / horizon_radius,
            arguments.minimum_boundary_horizon_radii,
        ),
        "boundary_magnetization_proxy": _condition(
            "maximum",
            float(boundary_state["maximum_b_squared_over_density_proxy"]),
            arguments.maximum_boundary_magnetization_proxy,
        ),
        "initial_volume_flux_divergence": _condition(
            "maximum",
            float(initial_flux["maximum_relative_cell_flux_divergence"]),
            arguments.maximum_initial_volume_flux_divergence,
        ),
    }
    return {
        "passed": all(entry["passed"] for entry in conditions.values()),
        "conditions": conditions,
        "secondary_horizon_radius": horizon_radius,
        "cell_spacing": spacing,
        "metric_model": "isolated analytic secondary Kerr",
        "science_ready": False,
        "science_blocker": (
            "the volume does not yet replay the transformed numerical ADM metric"
        ),
    }


def _mesh_arguments(cells: int, half_width: float) -> list[str]:
    meshblock = next(
        value for value in range(min(cells, 16), 0, -1) if cells % value == 0
    )
    result = []
    for axis in (1, 2, 3):
        result.extend(
            (
                f"mesh/nx{axis}={cells}",
                f"mesh/x{axis}min={-half_width:.17g}",
                f"mesh/x{axis}max={half_width:.17g}",
                f"meshblock/nx{axis}={meshblock}",
            )
        )
    return result


def _profile_arguments(profile: dict[str, object]) -> list[str]:
    result = [
        f"problem/rho0={profile['rho0']:.17g}",
        f"problem/pgas0={profile['pgas0']:.17g}",
        "problem/wind_frame=normal_frame",
    ]
    for component, value in enumerate(profile["u"], start=1):
        result.append(f"problem/u{component}={value:.17g}")
    for component, value in enumerate(profile["b"], start=1):
        result.append(f"problem/b{component}={value:.17g}")
    for direction, value in enumerate(profile["log_density_gradient"], start=1):
        result.append(f"problem/dlnrho_dxh{direction}={value:.17g}")
    for direction, value in enumerate(profile["log_pressure_gradient"], start=1):
        result.append(f"problem/dlnpgas_dxh{direction}={value:.17g}")
    for component, row in enumerate(profile["velocity_gradient"], start=1):
        for direction, value in enumerate(row, start=1):
            result.append(f"problem/du{component}_dxh{direction}={value:.17g}")
    return result


def _floor_controls(
    profile: dict[str, object], boundary_state: dict[str, object], half_width: float
) -> dict[str, float]:
    density_gradient = np.asarray(
        profile["log_density_gradient"], dtype=np.float64
    )
    pressure_gradient = np.asarray(
        profile["log_pressure_gradient"], dtype=np.float64
    )
    fitted_minimum_density = float(profile["rho0"]) * math.exp(
        -half_width * float(np.sum(np.abs(density_gradient)))
    )
    fitted_maximum_density = float(profile["rho0"]) * math.exp(
        half_width * float(np.sum(np.abs(density_gradient)))
    )
    fitted_minimum_pressure = float(profile["pgas0"]) * math.exp(
        -half_width * float(np.sum(np.abs(pressure_gradient)))
    )
    minimum_density = min(
        float(boundary_state["minimum_density"]), fitted_minimum_density
    )
    minimum_temperature = min(
        float(boundary_state["minimum_pressure"])
        / float(boundary_state["maximum_density"]),
        fitted_minimum_pressure / fitted_maximum_density,
    )
    tiny = np.finfo(float).tiny
    return {
        "density_floor": max(tiny, 1.0e-4 * minimum_density),
        "temperature_floor": max(tiny, 1.0e-4 * minimum_temperature),
    }


def _run(command: list[str], log_path: Path) -> float:
    start = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as stream:
        completed = subprocess.run(
            command,
            cwd=ROOT,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    elapsed = time.perf_counter() - start
    if completed.returncode != 0:
        tail = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
        raise RuntimeError(
            f"AthenaK replay failed with exit code {completed.returncode}:\n"
            + "\n".join(tail[-50:])
        )
    return elapsed


def _log_diagnostics(path: Path) -> dict[str, object]:
    content = path.read_text(encoding="utf-8", errors="replace")
    count = re.search(
        r"characteristic GR boundary: projections=(\d+), fallbacks=(\d+)",
        content,
    )
    residuals = [
        float(value)
        for value in re.findall(
            r"maximum normalized boundary-flux residual=([0-9.eE+-]+)", content
        )
    ]
    projections = 0 if count is None else int(count.group(1))
    fallbacks = 0 if count is None else int(count.group(2))
    total = projections + fallbacks
    return {
        "projection_count": projections,
        "fallback_count": fallbacks,
        "fallback_fraction": fallbacks / total if total else None,
        "maximum_boundary_flux_residual": max(residuals) if residuals else None,
    }


def _field_summary(path: Path) -> dict[str, object]:
    data = extract.bin_convert.read_binary(str(path))["mb_data"]
    result = {}
    for name, blocks in data.items():
        values = np.concatenate([np.asarray(block).ravel() for block in blocks])
        result[name] = {
            "minimum": float(np.min(values)),
            "maximum": float(np.max(values)),
            "finite": bool(np.isfinite(values).all()),
        }
    return result


def _expected_adm_fields(
    volume: adm_volume.ADMVolume,
    table_time: float,
    cells: int,
    halo: int,
) -> dict[str, np.ndarray]:
    """Reproduce the C++ time interpolation and ADM decomposition on active cells."""

    tolerance = 128.0 * np.finfo(float).eps * max(
        abs(table_time), abs(float(volume.times[0])), abs(float(volume.times[-1])), 1.0
    )
    if (
        table_time < volume.times[0] - tolerance
        or table_time > volume.times[-1] + tolerance
    ):
        raise ValueError("Athena output time lies outside the ADM replay table")
    interval = int(np.searchsorted(volume.times, table_time, side="right") - 1)
    interval = min(max(interval, 0), volume.times.size - 2)
    fraction = (table_time - volume.times[interval]) / (
        volume.times[interval + 1] - volume.times[interval]
    )
    fraction = min(max(float(fraction), 0.0), 1.0)
    fields = (
        (1.0 - fraction) * volume.fields[interval]
        + fraction * volume.fields[interval + 1]
    )
    active = fields[
        :,
        halo : halo + cells,
        halo : halo + cells,
        halo : halo + cells,
    ]
    if active.shape != (len(adm_volume.FIELD_NAMES), cells, cells, cells):
        raise ValueError("ADM volume does not contain the requested active-cell cube")
    metric = np.zeros((cells, cells, cells, 4, 4), dtype=np.float64)
    for field, (left, right) in enumerate(adm_volume.METRIC_COMPONENTS):
        metric[..., left, right] = active[field]
        metric[..., right, left] = active[field]
    adm, _ = adm_volume.decompose_four_metric(metric)
    gamma = adm["gamma"]
    curvature = active[len(adm_volume.METRIC_COMPONENTS) :]
    expected = {
        "adm_gxx": gamma[..., 0, 0],
        "adm_gxy": gamma[..., 0, 1],
        "adm_gxz": gamma[..., 0, 2],
        "adm_gyy": gamma[..., 1, 1],
        "adm_gyz": gamma[..., 1, 2],
        "adm_gzz": gamma[..., 2, 2],
        "adm_Kxx": curvature[0],
        "adm_Kxy": curvature[1],
        "adm_Kxz": curvature[2],
        "adm_Kyy": curvature[3],
        "adm_Kyz": curvature[4],
        "adm_Kzz": curvature[5],
        "adm_psi4": np.cbrt(np.linalg.det(gamma)),
        "adm_alpha": adm["alpha"],
        "adm_betax": adm["beta"][..., 0],
        "adm_betay": adm["beta"][..., 1],
        "adm_betaz": adm["beta"][..., 2],
    }
    return expected


def _adm_replay_comparison(
    paths: list[Path],
    volume: adm_volume.ADMVolume,
    cells: int,
    halo: int,
) -> dict[str, object]:
    if not paths:
        raise ValueError("numerical ADM pilot produced no ADM output")
    samples = []
    maximum_relative_error = 0.0
    maximum_absolute_error = 0.0
    for path in paths:
        data = extract.bin_convert.read_binary(str(path))
        athena_time = float(data["time"])
        table_time = float(volume.times[0]) + athena_time
        expected = _expected_adm_fields(volume, table_time, cells, halo)
        field_errors = {}
        for name in ATHENA_ADM_FIELDS:
            observed = extract.closure.assemble_uniform_grid(data, name)
            difference = observed - expected[name]
            absolute = float(np.max(np.abs(difference)))
            scale = max(float(np.max(np.abs(expected[name]))), np.finfo(float).tiny)
            relative = absolute / scale
            field_errors[name] = {
                "maximum_absolute": absolute,
                "relative_linf": relative,
            }
            maximum_absolute_error = max(maximum_absolute_error, absolute)
            maximum_relative_error = max(maximum_relative_error, relative)
        samples.append(
            {
                "path": str(path),
                "athena_time": athena_time,
                "replay_table_time": table_time,
                "fields": field_errors,
            }
        )
    return {
        "samples": samples,
        "maximum_absolute_error": maximum_absolute_error,
        "maximum_relative_linf": maximum_relative_error,
    }


def _history_summary(paths: list[Path]) -> dict[str, object]:
    if not paths:
        return {
            "valid": False,
            "error": "no force/accretion history file was produced",
        }
    try:
        history = force_history.ForceHistory("pilot", paths[-1])
    except (OSError, ValueError) as error:
        return {"valid": False, "error": str(error)}
    return {
        "valid": True,
        "path": str(paths[-1]),
        "sample_count": len(history.times),
        "time_range": [history.times[0], history.times[-1]],
        "accretion_rate_column": history.mdot_name,
        "force_frame": (
            "source_tetrad" if history.force_is_source_tetrad else "coordinate"
        ),
    }


def run_pilot(arguments: argparse.Namespace) -> dict[str, object]:
    source_path, source = _read_campaign(arguments.campaign)
    case_name, case = _case_document(source, arguments.case)
    worldtube_path = Path(case["worldtube"]).expanduser().resolve(strict=True)
    times, faces, metadata = worldtube.read_worldtube(worldtube_path)
    center = np.asarray(metadata.get("center", ()), dtype=np.float64)
    if center.shape != (3,) or not np.allclose(center, 0.0, rtol=0.0, atol=1.0e-14):
        raise ValueError(
            "isolated-secondary pilot requires a cube centered at the origin"
        )
    cells = faces[worldtube.FACE_NAMES[0]].cell_state.shape[-1]
    half_width = float(metadata["half_width"])
    profile = fit_initial_affine_profile(faces, metadata)
    assessment = structural_assessment(arguments, case, cells, half_width)
    if arguments.numerical_adm_volume:
        assessment.update(
            {
                "metric_model": (
                    "transformed numerical primary ADM plus analytic secondary Kerr"
                ),
                "science_ready": False,
                "science_blocker": (
                    "the off-grid force/accretion diagnostic is not yet connected "
                    "to the replayed numerical metric"
                ),
            }
        )
    boundary_state = case["boundary_state"]
    if not isinstance(boundary_state, dict):
        raise ValueError("campaign case has no boundary-state audit")
    floor_controls = _floor_controls(profile, boundary_state, half_width)
    workdir = arguments.workdir.expanduser().resolve()
    if workdir.exists():
        raise FileExistsError(f"refusing to reuse existing pilot directory {workdir}")
    workdir.mkdir(parents=True)
    replay_path = workdir / "inner.bin"
    inner_validation = worldtube.write_inner_binary(
        replay_path, times, faces, metadata
    )
    report: dict[str, object] = {
        "classification": CLASSIFICATION,
        "campaign": str(source_path),
        "case": case_name,
        "source_worldtube": str(worldtube_path),
        "replay_binary": str(replay_path),
        "sample_times": times.tolist(),
        "cells_per_axis": cells,
        "half_width": half_width,
        "secondary_mass": arguments.secondary_mass,
        "secondary_chi": arguments.secondary_chi,
        "initial_affine_profile": profile,
        "floor_controls": floor_controls,
        "inner_validation": inner_validation,
        "structural_assessment": assessment,
        "adm_volume_status": (
            "pending" if arguments.numerical_adm_volume else "disabled"
        ),
        "run_status": "refused",
    }
    if not assessment["passed"] and not arguments.allow_unsafe_structural_smoke:
        report["refusal_reason"] = (
            "structural replay gates failed; use --allow-unsafe-structural-smoke "
            "only for implementation debugging"
        )
        return report

    duration = float(times[-1] - times[0])
    tlim = duration if arguments.tlim is None else arguments.tlim
    if tlim > duration * (1.0 + 128.0 * np.finfo(float).eps):
        raise ValueError("pilot tlim extends beyond the replay time table")
    metric_volume = None
    metric_path = workdir / "adm.bin"
    if arguments.numerical_adm_volume:
        try:
            manifest_path = Path(case["manifest"]).expanduser().resolve(strict=True)
            metric_volume = adm_volume.build_volume(
                manifest_path,
                times,
                half_width,
                cells,
                arguments.metric_halo,
                arguments.secondary_mass,
                arguments.secondary_chi,
            )
            metric_validation = adm_volume.write_binary(metric_path, metric_volume)
        except (KeyError, OSError, RuntimeError, ValueError) as error:
            report.update(
                {
                    "adm_volume_status": "failed",
                    "adm_volume_error": str(error),
                    "refusal_reason": (
                        "numerical ADM volume extraction or validation failed"
                    ),
                }
            )
            return report
        report.update(
            {
                "adm_volume_status": "validated",
                "adm_volume": str(metric_path),
                "adm_volume_validation": metric_validation,
                "metric_halo": arguments.metric_halo,
                "mesh_nghost": arguments.mesh_nghost,
            }
        )
    horizon = float(assessment["secondary_horizon_radius"])
    force_radii = (
        max(1.5 * horizon, 0.25 * half_width),
        0.55 * half_width,
        0.70 * half_width,
        0.875 * half_width,
    )
    if not all(left < right for left, right in zip(force_radii[:-1], force_radii[1:])):
        raise ValueError("secondary is too large for the configured force shells")
    output = workdir / "run"
    output.mkdir()
    log_path = workdir / "athena.log"
    history_dt = tlim / arguments.history_samples
    command = [
        str(arguments.athena.expanduser().resolve(strict=True)),
        "-i",
        str(arguments.input.expanduser().resolve(strict=True)),
        "-d",
        str(output),
        "job/basename=inner",
        *_mesh_arguments(cells, half_width),
        f"mesh/nghost={arguments.mesh_nghost}",
        *_profile_arguments(profile),
        "time/integrator=rk3",
        f"time/cfl_number={arguments.cfl:.17g}",
        "time/nlim=-1",
        f"time/tlim={tlim:.17g}",
        f"mhd/dfloor={floor_controls['density_floor']:.17g}",
        f"mhd/tfloor={floor_controls['temperature_floor']:.17g}",
        "problem/background_mode=isolated",
        "problem/secondary_embedding=tangent_tetrad",
        "problem/primary_mass=1",
        f"problem/secondary_mass={arguments.secondary_mass:.17g}",
        "problem/primary_chi=0",
        f"problem/secondary_chi={arguments.secondary_chi:.17g}",
        "problem/orbital_radius=10",
        "problem/require_stable_orbit=false",
        f"problem/user_hist={'false' if metric_volume is not None else 'true'}",
        "problem/force_frame=coordinate",
        f"problem/force_surface_radius={force_radii[0]:.17g}",
        f"problem/force_outer_radius_1={force_radii[1]:.17g}",
        f"problem/force_outer_radius_2={force_radii[2]:.17g}",
        f"problem/force_outer_radius_3={force_radii[3]:.17g}",
        "emri_worldtube/enabled=true",
        "emri_worldtube/mode=inner",
        f"emri_worldtube/file={replay_path}",
        "emri_worldtube/fluid_boundary=characteristic_gr",
        "emri_worldtube/fluid_state_frame=inner_coordinate",
        "output1/variable=mhd_w_bcc",
        f"output1/dt={tlim:.17g}",
        "output2/variable=mhd_divb",
        f"output2/dt={tlim:.17g}",
        "output3/dt=0",
        f"output4/dt={history_dt:.17g}" if metric_volume is None else "output4/dt=0",
        "output4/user_hist_only=true",
        "output5/variable=adm",
        f"output5/dt={tlim:.17g}" if metric_volume is not None else "output5/dt=0",
    ]
    if metric_volume is not None:
        command.extend(
            (
                "adm/dynamic=true",
                "emri_adm_replay/enabled=true",
                f"emri_adm_replay/file={metric_path}",
            )
        )
    wall_seconds = _run(command, log_path)
    states = sorted((output / "bin").glob("inner.mhd_w_bcc.*.bin"))
    divergences = sorted((output / "bin").glob("inner.mhd_divb.*.bin"))
    if not states or not divergences:
        raise RuntimeError("pilot did not produce final state and divB outputs")
    state_summary = _field_summary(states[-1])
    divergence_summary = _field_summary(divergences[-1])["divb"]
    log_diagnostics = _log_diagnostics(log_path)
    boundary_residual = log_diagnostics["maximum_boundary_flux_residual"]
    fallback_fraction = log_diagnostics["fallback_fraction"]
    runtime_conditions = {
        "characteristic_fallback_fraction": _condition(
            "maximum",
            math.inf if fallback_fraction is None else float(fallback_fraction),
            arguments.maximum_fallback_fraction,
        ),
        "boundary_flux_residual": _condition(
            "maximum",
            math.inf if boundary_residual is None else float(boundary_residual),
            arguments.maximum_boundary_flux_residual,
        ),
        "maximum_divb": _condition(
            "maximum",
            max(
                abs(float(divergence_summary["minimum"])),
                abs(float(divergence_summary["maximum"])),
            ),
            arguments.maximum_divb,
        ),
    }
    finite_state = all(entry["finite"] for entry in state_summary.values())
    positive_state = (
        float(state_summary["dens"]["minimum"]) > 0.0
        and float(state_summary["press"]["minimum"]) > 0.0
    )
    runtime_conditions["finite_positive_fluid_state"] = {
        "relation": "required",
        "observed": finite_state and positive_state,
        "threshold": True,
        "passed": finite_state and positive_state,
    }
    histories = sorted(output.rglob("*.hst"))
    if metric_volume is None:
        history_summary = _history_summary(histories)
        runtime_conditions["force_accretion_history"] = {
            "relation": "required",
            "observed": history_summary["valid"],
            "threshold": True,
            "passed": history_summary["valid"],
        }
        adm_comparison = None
    else:
        history_summary = {
            "valid": False,
            "disabled_reason": (
                "force/accretion diagnostic still samples the analytic off-grid metric"
            ),
        }
        adm_outputs = sorted((output / "bin").glob("inner.adm.*.bin"))
        adm_comparison = _adm_replay_comparison(
            adm_outputs,
            metric_volume,
            cells,
            arguments.metric_halo,
        )
        runtime_conditions["adm_endpoint_replay"] = _condition(
            "maximum",
            float(adm_comparison["maximum_relative_linf"]),
            arguments.maximum_adm_replay_relative_error,
        )
    report.update(
        {
            "run_status": "completed",
            "command": command,
            "wall_seconds": wall_seconds,
            "log": str(log_path),
            "final_state": str(states[-1]),
            "final_divb": str(divergences[-1]),
            "history_files": [str(path) for path in histories],
            "history_summary": history_summary,
            "adm_replay_comparison": adm_comparison,
            "log_diagnostics": log_diagnostics,
            "final_state_summary": state_summary,
            "final_divb_summary": divergence_summary,
            "runtime_assessment": {
                "passed": all(item["passed"] for item in runtime_conditions.values()),
                "conditions": runtime_conditions,
            },
        }
    )
    return report


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--case")
    parser.add_argument("--athena", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--secondary-mass", type=float, required=True)
    parser.add_argument("--secondary-chi", type=float, default=0.0)
    parser.add_argument("--tlim", type=float)
    parser.add_argument("--cfl", type=float, default=0.02)
    parser.add_argument("--history-samples", type=int, default=16)
    parser.add_argument("--numerical-adm-volume", action="store_true")
    parser.add_argument("--metric-halo", type=int, default=4)
    parser.add_argument("--mesh-nghost", type=int, default=4)
    parser.add_argument("--minimum-horizon-cells", type=float, default=4.0)
    parser.add_argument(
        "--minimum-boundary-horizon-radii", type=float, default=5.0
    )
    parser.add_argument(
        "--maximum-boundary-magnetization-proxy", type=float, default=1.0e3
    )
    parser.add_argument(
        "--maximum-initial-volume-flux-divergence", type=float, default=1.0e-10
    )
    parser.add_argument("--maximum-fallback-fraction", type=float, default=5.0e-4)
    parser.add_argument(
        "--maximum-boundary-flux-residual", type=float, default=1.0e-10
    )
    parser.add_argument("--maximum-divb", type=float, default=1.0e-10)
    parser.add_argument(
        "--maximum-adm-replay-relative-error", type=float, default=1.0e-5
    )
    parser.add_argument("--allow-unsafe-structural-smoke", action="store_true")
    parser.add_argument("--fail-on-gate", action="store_true")
    arguments = parser.parse_args()
    positive = (
        "secondary_mass",
        "cfl",
        "minimum_horizon_cells",
        "minimum_boundary_horizon_radii",
        "maximum_boundary_magnetization_proxy",
        "maximum_initial_volume_flux_divergence",
        "maximum_boundary_flux_residual",
        "maximum_divb",
        "maximum_adm_replay_relative_error",
    )
    for name in positive:
        value = getattr(arguments, name)
        if not math.isfinite(value) or value <= 0.0:
            parser.error(f"--{name.replace('_', '-')} must be finite and positive")
    if not math.isfinite(arguments.secondary_chi) or abs(arguments.secondary_chi) > 1.0:
        parser.error("--secondary-chi must lie in [-1,1]")
    if arguments.tlim is not None and (
        not math.isfinite(arguments.tlim) or arguments.tlim <= 0.0
    ):
        parser.error("--tlim must be finite and positive")
    if arguments.history_samples < 1:
        parser.error("--history-samples must be positive")
    if arguments.metric_halo < 1 or arguments.mesh_nghost < 1:
        parser.error("--metric-halo and --mesh-nghost must be positive")
    if arguments.metric_halo < arguments.mesh_nghost:
        parser.error("--metric-halo must be at least --mesh-nghost")
    if (
        not math.isfinite(arguments.maximum_fallback_fraction)
        or not 0.0 <= arguments.maximum_fallback_fraction <= 1.0
    ):
        parser.error("--maximum-fallback-fraction must lie in [0,1]")
    return arguments


def main() -> int:
    arguments = parse_arguments()
    report = run_pilot(arguments)
    output = arguments.workdir.expanduser().resolve() / "summary.json"
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"run_status={report['run_status']}")
    print(
        "structural_gates="
        f"{'PASS' if report['structural_assessment']['passed'] else 'FAIL'}"
    )
    if report["run_status"] == "completed":
        runtime = report["runtime_assessment"]
        print(f"runtime_gates={'PASS' if runtime['passed'] else 'FAIL'}")
    print(output)
    if report["run_status"] == "refused":
        return 2
    if arguments.fail_on_gate and not report["runtime_assessment"]["passed"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
