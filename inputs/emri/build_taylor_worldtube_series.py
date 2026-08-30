#!/usr/bin/env python3
"""Build a time-dependent EMRI Taylor-profile table from global GRMHD dumps.

The input manifest supplies AthenaK ``mhd_w_bcc`` dumps with the secondary
worldline position and global coordinate velocity.  The primary geometry can be
a co-temporal numerical ``adm`` dump per sample or one exact, stationary Kerr
metric declared at manifest level.  Every sample is reduced with
``extract_static_taylor_worldtube``; the resulting source-tetrad profiles are
written in the strict column order consumed by ``emri_windtunnel``.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np

import extract_static_taylor_worldtube as static


INPUT_CLASSIFICATION = "athenak-emri-worldline-v1"
TABLE_CLASSIFICATION = "athenak-emri-taylor-series-v2"


def _finite_vector(values: Iterable[float], name: str) -> np.ndarray:
    result = np.asarray(tuple(values), dtype=float)
    if result.shape != (3,) or not np.isfinite(result).all():
        raise ValueError(f"{name} must contain three finite values")
    return result


def _resolve_dump(path: str, manifest_directory: Path) -> Path:
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        candidate = manifest_directory / candidate
    return candidate.resolve(strict=True)


def _cotemporal(state: dict, adm: dict) -> None:
    state_time = float(state["time"])
    adm_time = float(adm["time"])
    tolerance = 64.0 * np.finfo(float).eps * max(
        1.0, abs(state_time), abs(adm_time)
    )
    if int(state["cycle"]) != int(adm["cycle"]) or not math.isclose(
        state_time, adm_time, rel_tol=0.0, abs_tol=tolerance
    ):
        raise RuntimeError("state and ADM dumps are not co-temporal")


def extract_entry(
    entry: dict,
    manifest_directory: Path,
    primary_center: np.ndarray,
    disk_normal: np.ndarray,
    fit_radius: float,
    metric_fit_radius: float,
    global_length_in_local_units: float,
    density_renormalization: float,
    analytic_kerr: tuple[float, float] | None = None,
) -> dict[str, object]:
    required = {"state", "anchor", "source_velocity"}
    if analytic_kerr is None:
        required.add("adm")
    missing = sorted(required.difference(entry))
    if missing:
        raise ValueError(f"worldline sample is missing: {', '.join(missing)}")
    state_path = _resolve_dump(str(entry["state"]), manifest_directory)
    adm_path = (
        None
        if analytic_kerr is not None
        else _resolve_dump(str(entry["adm"]), manifest_directory)
    )
    anchor = _finite_vector(entry["anchor"], "anchor")
    source_velocity = _finite_vector(entry["source_velocity"], "source velocity")
    state, state_thermodynamics = static.read_state_dump(state_path)
    if adm_path is not None:
        adm = static.bin_convert.read_binary(str(adm_path))
        static._check_dump(adm, static.ADM_VARIABLES, "ADM")
        _cotemporal(state, adm)
        adm_order = static._aligned_adm_blocks(state, adm)
        gamma, alpha, beta, metric_errors = static.fit_anchor_adm(
            adm, anchor, metric_fit_radius
        )
        metric_source = {
            "classification": "athenak-emri-numerical-adm-primary-v1"
        }
    else:
        assert analytic_kerr is not None
        adm = None
        adm_order = None
        gamma, alpha, beta = static.analytic_kerr_anchor_adm(
            anchor, *analytic_kerr
        )
        metric_errors = {"analytic_primary_metric": 0.0}
        metric_source = static.analytic_kerr_metric_source(*analytic_kerr)
    if "time" in entry:
        requested_time = float(entry["time"])
        tolerance = 64.0 * np.finfo(float).eps * max(
            1.0, abs(requested_time), abs(float(state["time"]))
        )
        if not math.isclose(
            requested_time, float(state["time"]), rel_tol=0.0, abs_tol=tolerance
        ):
            raise RuntimeError("manifest time does not match dump time")

    spatial_basis = static.canonical_spatial_basis(
        anchor, primary_center, disk_normal
    )
    metric = static.spacetime_metric_from_adm(gamma, alpha, beta)
    tetrad, coframe = static.build_source_tetrad(
        metric, source_velocity, spatial_basis
    )
    cloud = static.collect_local_samples(
        state,
        adm,
        adm_order,
        anchor,
        coframe,
        fit_radius,
        analytic_kerr=analytic_kerr,
    )
    parameters, diagnostics = static.fit_static_profile(
        cloud, coframe, fit_radius
    )
    parameters = static.rescale_profile_parameters(
        parameters, global_length_in_local_units, density_renormalization
    )
    diagnostics = static.rescale_profile_diagnostics(
        diagnostics, global_length_in_local_units, density_renormalization
    )
    global_time = float(state["time"])
    result: dict[str, object] = {
        "time": global_time * global_length_in_local_units,
        "time_global_units": global_time,
        "cycle": int(state["cycle"]),
        "state_file": str(state_path),
        "state_sha256": static.file_sha256(state_path),
        "metric_source": metric_source,
        "state_thermodynamics": state_thermodynamics,
        "anchor_global": anchor.tolist(),
        "source_coordinate_velocity": source_velocity.tolist(),
        "source_tetrad_contravariant": tetrad.tolist(),
        "source_tetrad_coframe": coframe.tolist(),
        "parameters": parameters,
        "fit_diagnostics": diagnostics,
        "anchor_metric_relative_fit_errors": metric_errors,
    }
    if adm_path is not None:
        result["adm_file"] = str(adm_path)
        result["adm_sha256"] = static.file_sha256(adm_path)
    return result


def circular_orbit_diagnostics(
    samples: list[dict[str, object]],
    primary_center: np.ndarray,
    disk_normal: np.ndarray,
) -> dict[str, float]:
    radii = []
    heights = []
    radial_speeds = []
    vertical_speeds = []
    angular_frequencies = []
    tangential_speeds = []
    for sample in samples:
        separation = np.asarray(sample["anchor_global"]) - primary_center
        velocity = np.asarray(sample["source_coordinate_velocity"])
        height = float(separation @ disk_normal)
        in_plane = separation - height * disk_normal
        radius = float(np.linalg.norm(in_plane))
        if not radius > 0.0:
            raise RuntimeError("worldline intersects the primary/disk axis")
        radial = in_plane / radius
        tangent = np.cross(disk_normal, radial)
        radii.append(radius)
        heights.append(height)
        radial_speeds.append(float(velocity @ radial))
        vertical_speeds.append(float(velocity @ disk_normal))
        tangential_speed = float(velocity @ tangent)
        tangential_speeds.append(tangential_speed)
        angular_frequencies.append(tangential_speed / radius)

    mean_radius = float(np.mean(radii))
    mean_omega = float(np.mean(angular_frequencies))
    speed_scale = max(float(np.max(np.abs(tangential_speeds))), 1.0)
    return {
        "mean_coordinate_radius": mean_radius,
        "maximum_fractional_radius_variation": float(
            np.max(np.abs(np.asarray(radii) - mean_radius)) / mean_radius
        ),
        "maximum_fractional_height": float(
            np.max(np.abs(heights)) / mean_radius
        ),
        "mean_coordinate_angular_frequency": mean_omega,
        "maximum_fractional_angular_frequency_variation": float(
            np.max(np.abs(np.asarray(angular_frequencies) - mean_omega))
            / max(abs(mean_omega), np.finfo(float).tiny)
        ),
        "maximum_relative_radial_speed": float(
            np.max(np.abs(radial_speeds)) / speed_scale
        ),
        "maximum_relative_vertical_speed": float(
            np.max(np.abs(vertical_speeds)) / speed_scale
        ),
    }


def validate_circular_orbit(diagnostics: dict[str, float], tolerance: float) -> None:
    checked = (
        "maximum_fractional_radius_variation",
        "maximum_fractional_height",
        "maximum_fractional_angular_frequency_variation",
        "maximum_relative_radial_speed",
        "maximum_relative_vertical_speed",
    )
    violations = [name for name in checked if diagnostics[name] > tolerance]
    if violations:
        values = ", ".join(f"{name}={diagnostics[name]:.6g}" for name in violations)
        raise RuntimeError(
            "worldline is incompatible with the current circular-equatorial local "
            f"metric ({values}; tolerance={tolerance:.6g})"
        )


def write_table(
    path: Path,
    samples: list[dict[str, object]],
    global_length_in_local_units: float = 1.0,
    density_renormalization: float = 1.0,
    source_coordinate_radius_local_units: float | None = None,
    source_coordinate_angular_frequency_local_units: float | None = None,
    orbit_tolerance: float = 1.0e-6,
) -> None:
    lines = [
        f"# {TABLE_CLASSIFICATION}",
        "# global_length_in_local_units: "
        + static._format_real(global_length_in_local_units),
        "# density_renormalization: "
        + static._format_real(density_renormalization),
    ]
    if source_coordinate_radius_local_units is not None:
        lines.append(
            "# source_coordinate_radius_local_units: "
            + static._format_real(source_coordinate_radius_local_units)
        )
    if source_coordinate_angular_frequency_local_units is not None:
        lines.append(
            "# source_coordinate_angular_frequency_local_units: "
            + static._format_real(
                source_coordinate_angular_frequency_local_units
            )
        )
    lines.append("# orbit_tolerance: " + static._format_real(orbit_tolerance))
    lines.append("# columns: time " + " ".join(static.PROFILE_PARAMETER_ORDER))
    for sample in samples:
        values = [float(sample["time"])]
        parameters = sample["parameters"]
        values.extend(float(parameters[name]) for name in static.PROFILE_PARAMETER_ORDER)
        lines.append(" ".join(static._format_real(value) for value in values))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_table(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = []
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        values = np.fromstring(line, sep=" ")
        expected = 1 + len(static.PROFILE_PARAMETER_ORDER)
        if values.shape != (expected,) or not np.isfinite(values).all():
            raise RuntimeError(f"invalid Taylor table row {line_number}")
        rows.append(values)
    if len(rows) < 2:
        raise RuntimeError("Taylor profile table requires at least two rows")
    table = np.asarray(rows)
    if np.any(np.diff(table[:, 0]) <= 0.0):
        raise RuntimeError("Taylor profile times must increase strictly")
    return table[:, 0], table[:, 1:]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--primary-center", type=float, nargs=3, default=(0, 0, 0))
    parser.add_argument("--disk-normal", type=float, nargs=3, default=(0, 0, 1))
    parser.add_argument("--fit-radius", type=float, required=True)
    parser.add_argument("--metric-fit-radius", type=float)
    parser.add_argument("--global-length-in-local-units", type=float, default=1.0)
    parser.add_argument("--density-renormalization", type=float)
    parser.add_argument("--orbit-tolerance", type=float, default=1.0e-6)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest_path = args.manifest.expanduser().resolve(strict=True)
    source_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if source_manifest.get("classification") != INPUT_CLASSIFICATION:
        raise RuntimeError(
            f"manifest classification must be {INPUT_CLASSIFICATION!r}"
        )
    entries = source_manifest.get("samples")
    if not isinstance(entries, list) or len(entries) < 2:
        raise RuntimeError("worldline manifest requires at least two samples")
    try:
        analytic_kerr = static.parse_analytic_kerr_metric_source(
            source_manifest.get("metric_source")
        )
    except (TypeError, ValueError) as error:
        raise RuntimeError(str(error)) from error
    if analytic_kerr is not None and any("adm" in entry for entry in entries):
        raise RuntimeError(
            "analytic-Kerr worldline samples must not also declare ADM dumps"
        )
    if not args.fit_radius > 0.0 or not math.isfinite(args.fit_radius):
        raise SystemExit("--fit-radius must be finite and positive")
    metric_fit_radius = args.metric_fit_radius or args.fit_radius
    if not metric_fit_radius > 0.0 or not math.isfinite(metric_fit_radius):
        raise SystemExit("--metric-fit-radius must be finite and positive")
    if not args.orbit_tolerance > 0.0 or not math.isfinite(args.orbit_tolerance):
        raise SystemExit("--orbit-tolerance must be finite and positive")
    if (
        not args.global_length_in_local_units > 0.0
        or not math.isfinite(args.global_length_in_local_units)
    ):
        raise SystemExit("--global-length-in-local-units must be finite and positive")
    density_renormalization = args.density_renormalization
    if density_renormalization is None:
        if args.global_length_in_local_units != 1.0:
            raise SystemExit(
                "--density-renormalization must be chosen explicitly when global and "
                "local length units differ"
            )
        density_renormalization = 1.0
    if not density_renormalization > 0.0 or not math.isfinite(
        density_renormalization
    ):
        raise SystemExit("--density-renormalization must be finite and positive")
    primary_center = _finite_vector(args.primary_center, "primary center")
    disk_normal = static._normalize(
        _finite_vector(args.disk_normal, "disk normal"), "disk normal"
    )

    samples = [
        extract_entry(
            entry,
            manifest_path.parent,
            primary_center,
            disk_normal,
            args.fit_radius,
            metric_fit_radius,
            args.global_length_in_local_units,
            density_renormalization,
            analytic_kerr,
        )
        for entry in entries
    ]
    times = np.asarray([sample["time"] for sample in samples])
    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0.0):
        raise RuntimeError("worldline dump times must increase strictly")
    global_times = np.asarray([sample["time_global_units"] for sample in samples])
    local_deltas = np.diff(times)
    global_deltas = np.diff(global_times)
    time_sampling = {
        "sample_count": len(samples),
        "minimum_interval_local_units": float(np.min(local_deltas)),
        "maximum_interval_local_units": float(np.max(local_deltas)),
        "median_interval_local_units": float(np.median(local_deltas)),
        "minimum_interval_global_units": float(np.min(global_deltas)),
        "maximum_interval_global_units": float(np.max(global_deltas)),
        "median_interval_global_units": float(np.median(global_deltas)),
    }
    orbit_diagnostics = circular_orbit_diagnostics(
        samples, primary_center, disk_normal
    )
    validate_circular_orbit(orbit_diagnostics, args.orbit_tolerance)
    orbit_diagnostics["mean_coordinate_radius_local_units"] = (
        orbit_diagnostics["mean_coordinate_radius"]
        * args.global_length_in_local_units
    )
    orbit_diagnostics["mean_coordinate_angular_frequency_local_units"] = (
        orbit_diagnostics["mean_coordinate_angular_frequency"]
        / args.global_length_in_local_units
    )

    prefix = args.output_prefix.expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    table_path = prefix.with_suffix(".dat")
    output_manifest_path = prefix.with_suffix(".json")
    write_table(
        table_path,
        samples,
        args.global_length_in_local_units,
        density_renormalization,
        orbit_diagnostics["mean_coordinate_radius_local_units"],
        orbit_diagnostics["mean_coordinate_angular_frequency_local_units"],
        args.orbit_tolerance,
    )
    read_table(table_path)
    output_manifest = {
        "schema": 1,
        "classification": TABLE_CLASSIFICATION,
        "source_manifest": str(manifest_path),
        "source_manifest_sha256": static.file_sha256(manifest_path),
        "table_file": str(table_path),
        "table_sha256": static.file_sha256(table_path),
        "parameter_order": list(static.PROFILE_PARAMETER_ORDER),
        "primary_center_global": primary_center.tolist(),
        "disk_normal_global": disk_normal.tolist(),
        "fit_radius_source": args.fit_radius,
        "metric_fit_radius_coordinate": metric_fit_radius,
        "global_length_in_local_units": args.global_length_in_local_units,
        "density_renormalization": density_renormalization,
        "metric_source": (
            static.analytic_kerr_metric_source(*analytic_kerr)
            if analytic_kerr is not None
            else {"classification": "athenak-emri-numerical-adm-primary-v1"}
        ),
        "time_sampling": time_sampling,
        "orbit_tolerance": args.orbit_tolerance,
        "orbit_diagnostics": orbit_diagnostics,
        "samples": samples,
        "limitations": [
            "one-way first-order Taylor replay; no response of the global disk",
            "current local metric requires a circular equatorial worldline",
            "time interpolation occurs between independently constrained snapshots",
        ],
    }
    output_manifest_path.write_text(
        json.dumps(output_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(table_path)
    print(output_manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
