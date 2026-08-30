#!/usr/bin/env python3
"""Build an audited affine frame for a circular equatorial Kerr orbit.

The orbit follows the test-particle angular frequency already used by the
``emri_windtunnel`` problem generator.  Its spatial legs are the orthonormal
radial, prograde-tangential, and vertical source tetrad in Cartesian
Kerr-Schild coordinates.  Analytic rotation derivatives complete the cubic
Hermite contract consumed by ``extract_global_worldtube.py``.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import extract_global_worldtube as extract
import extract_static_taylor_worldtube as static
import kerr_schild_background as kerr_background


GENERATOR_CLASSIFICATION = "athenak-emri-kerr-circular-frame-generator-v1"


def circular_kerr_omega(
    primary_mass: float,
    primary_spin: float,
    orbital_radius: float,
    direction: int,
) -> float:
    """Coordinate angular frequency in aligned-spin Kerr-Schild time."""

    if primary_mass <= 0.0 or orbital_radius <= 0.0 or direction not in (-1, 1):
        raise ValueError("circular Kerr orbit parameters are invalid")
    root_mass = math.sqrt(primary_mass)
    denominator = (
        orbital_radius * math.sqrt(orbital_radius)
        + direction * primary_spin * root_mass
    )
    if denominator <= 0.0:
        raise ValueError("circular Kerr angular-frequency denominator is non-positive")
    return direction * root_mass / denominator


def kerr_isco(primary_mass: float, dimensionless_spin: float, direction: int) -> float:
    """Equatorial Kerr ISCO radius in Boyer-Lindquist/Kerr-Schild r."""

    if (
        primary_mass <= 0.0
        or abs(dimensionless_spin) > 1.0
        or direction not in (-1, 1)
    ):
        raise ValueError("Kerr ISCO parameters are invalid")
    effective_spin = direction * dimensionless_spin
    z1 = 1.0 + np.cbrt(1.0 - effective_spin**2) * (
        np.cbrt(1.0 + effective_spin) + np.cbrt(1.0 - effective_spin)
    )
    z2 = math.sqrt(3.0 * effective_spin**2 + z1**2)
    signed_root = 0.0
    if effective_spin != 0.0:
        signed_root = math.copysign(
            math.sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2)),
            effective_spin,
        )
    return primary_mass * (3.0 + z2 - signed_root)


def kerr_schild_metric(
    position: object, primary_mass: float, primary_spin: float
) -> np.ndarray:
    """Unregularized aligned-spin Cartesian Kerr-Schild covariant metric."""

    return kerr_background.covariant_metric(
        position, primary_mass, primary_spin
    )


def _checked_times(values: object) -> np.ndarray:
    times = np.asarray(values, dtype=np.float64)
    if (
        times.ndim != 1
        or times.size < 2
        or not np.isfinite(times).all()
        or np.any(np.diff(times) <= 0.0)
    ):
        raise ValueError("frame times must be finite and strictly increasing")
    return times


def _exact_sample(
    local_time: float,
    primary_mass: float,
    primary_spin: float,
    coordinate_radius: float,
    omega: float,
    initial_phase: float,
    phase_reference_time: float,
    global_length_in_local_units: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    global_time = local_time / global_length_in_local_units
    phase = initial_phase + omega * (global_time - phase_reference_time)
    cosine = math.cos(phase)
    sine = math.sin(phase)
    radial = np.asarray((cosine, sine, 0.0))
    tangential = np.asarray((-sine, cosine, 0.0))
    position = coordinate_radius * radial
    coordinate_velocity = omega * coordinate_radius * tangential
    metric = kerr_schild_metric(position, primary_mass, primary_spin)
    basis = np.column_stack((radial, tangential, np.asarray((0.0, 0.0, 1.0))))
    tetrad, _ = static.build_source_tetrad(metric, coordinate_velocity, basis)
    tangent = np.concatenate(([1.0], coordinate_velocity))
    tangent /= global_length_in_local_units
    legs = tetrad[1:].T / global_length_in_local_units
    derivative = np.zeros_like(legs)
    local_omega = omega / global_length_in_local_units
    derivative[1, :] = -local_omega * legs[2, :]
    derivative[2, :] = local_omega * legs[1, :]
    return position, tangent, legs, derivative


def _audit_frame(
    document: dict[str, object],
    primary_mass: float,
    primary_spin: float,
    coordinate_radius: float,
    omega: float,
    initial_phase: float,
    phase_reference_time: float,
    global_length_in_local_units: float,
    subsamples: int,
) -> dict[str, float | int]:
    series = extract.AffineFrameSeries.from_document(document)
    audit_times = []
    for left, right in zip(series.times[:-1], series.times[1:], strict=True):
        audit_times.extend(np.linspace(left, right, subsamples + 1)[:-1])
    audit_times.append(series.times[-1])
    maximum_position_error = 0.0
    maximum_tangent_error = 0.0
    maximum_leg_error = 0.0
    maximum_gram_error = 0.0
    maximum_spatial_leg_time_row_norm = 0.0
    minimum_tangent_timelike_margin = math.inf
    minimum_jacobian_determinant = math.inf
    maximum_jacobian_condition = 0.0
    minkowski = np.diag((-1.0, 1.0, 1.0, 1.0))
    for time in audit_times:
        worldline, instantaneous = series.evaluate(float(time))
        exact_position, exact_tangent, exact_legs, _ = _exact_sample(
            float(time),
            primary_mass,
            primary_spin,
            coordinate_radius,
            omega,
            initial_phase,
            phase_reference_time,
            global_length_in_local_units,
        )
        exact_worldline = np.concatenate(
            ([float(time) / global_length_in_local_units], exact_position)
        )
        tangent = instantaneous.worldline_tangent
        metric = kerr_schild_metric(worldline[1:], primary_mass, primary_spin)
        norm2 = static.metric_inner(metric, tangent, tangent)
        if norm2 >= 0.0:
            raise RuntimeError("Hermite-interpolated circular worldline is not timelike")
        unit_tangent = tangent / math.sqrt(-norm2)
        tetrad = np.vstack(
            (
                unit_tangent,
                global_length_in_local_units * instantaneous.spatial_legs.T,
            )
        )
        gram_error = float(np.max(np.abs(tetrad @ metric @ tetrad.T - minkowski)))
        jacobian = instantaneous.jacobian(np.zeros(3))
        maximum_position_error = max(
            maximum_position_error,
            float(np.linalg.norm(worldline - exact_worldline)) / coordinate_radius,
        )
        maximum_tangent_error = max(
            maximum_tangent_error,
            float(np.linalg.norm(tangent - exact_tangent)),
        )
        maximum_leg_error = max(
            maximum_leg_error,
            float(np.max(np.abs(instantaneous.spatial_legs - exact_legs))),
        )
        maximum_gram_error = max(maximum_gram_error, gram_error)
        maximum_spatial_leg_time_row_norm = max(
            maximum_spatial_leg_time_row_norm,
            float(np.linalg.norm(instantaneous.spatial_legs[0])),
        )
        minimum_tangent_timelike_margin = min(
            minimum_tangent_timelike_margin,
            -norm2 * global_length_in_local_units**2,
        )
        minimum_jacobian_determinant = min(
            minimum_jacobian_determinant, abs(float(np.linalg.det(jacobian)))
        )
        maximum_jacobian_condition = max(
            maximum_jacobian_condition, float(np.linalg.cond(jacobian))
        )
    return {
        "audit_event_count": len(audit_times),
        "audit_subsamples_per_interval": subsamples,
        "maximum_worldline_position_error_over_coordinate_radius": (
            maximum_position_error
        ),
        "maximum_worldline_tangent_absolute_error": maximum_tangent_error,
        "maximum_spatial_leg_absolute_error": maximum_leg_error,
        "maximum_interpolated_tetrad_gram_error": maximum_gram_error,
        "maximum_spatial_leg_time_row_norm": (
            maximum_spatial_leg_time_row_norm
        ),
        "minimum_worldline_timelike_margin": minimum_tangent_timelike_margin,
        "minimum_origin_jacobian_absolute_determinant": (
            minimum_jacobian_determinant
        ),
        "minimum_scaled_origin_jacobian_absolute_determinant": (
            minimum_jacobian_determinant * global_length_in_local_units**4
        ),
        "maximum_origin_jacobian_condition_number": maximum_jacobian_condition,
    }


def build_frame_document(
    global_times: object,
    primary_mass: float,
    dimensionless_spin: float,
    orbital_radius: float,
    direction: int = 1,
    initial_phase: float = 0.0,
    phase_reference_time: float | None = None,
    require_outside_isco: bool = True,
    audit_subsamples: int = 4,
    global_length_in_local_units: float = 1.0,
) -> dict[str, object]:
    """Return a complete affine-frame-series JSON document."""

    checked_global_times = _checked_times(global_times)
    if (
        not math.isfinite(primary_mass)
        or primary_mass <= 0.0
        or not math.isfinite(dimensionless_spin)
        or abs(dimensionless_spin) > 1.0
        or not math.isfinite(orbital_radius)
        or orbital_radius <= 0.0
        or direction not in (-1, 1)
        or not math.isfinite(initial_phase)
        or not isinstance(audit_subsamples, int)
        or isinstance(audit_subsamples, bool)
        or audit_subsamples < 1
        or not math.isfinite(global_length_in_local_units)
        or global_length_in_local_units <= 0.0
    ):
        raise ValueError("Kerr circular-frame parameters are invalid")
    reference_time = (
        float(checked_global_times[0])
        if phase_reference_time is None
        else float(phase_reference_time)
    )
    if not math.isfinite(reference_time):
        raise ValueError("phase reference time must be finite")
    primary_spin = dimensionless_spin * primary_mass
    isco = kerr_isco(primary_mass, dimensionless_spin, direction)
    if require_outside_isco and orbital_radius < isco:
        raise ValueError(
            f"orbital radius {orbital_radius:.17g} lies inside Kerr ISCO "
            f"{isco:.17g}"
        )
    omega = circular_kerr_omega(
        primary_mass, primary_spin, orbital_radius, direction
    )
    coordinate_radius = math.sqrt(orbital_radius**2 + primary_spin**2)
    coordinate_speed = abs(omega * coordinate_radius)
    if coordinate_speed >= 1.0:
        raise ValueError(
            "circular trajectory is superluminal in Cartesian KS coordinates"
        )
    local_times = checked_global_times * global_length_in_local_units
    worldline = np.empty((local_times.size, 4))
    tangent = np.empty((local_times.size, 4))
    legs = np.empty((local_times.size, 4, 3))
    leg_derivative = np.empty_like(legs)
    for index, local_time in enumerate(local_times):
        position, sample_tangent, sample_legs, sample_derivative = _exact_sample(
            float(local_time),
            primary_mass,
            primary_spin,
            coordinate_radius,
            omega,
            initial_phase,
            reference_time,
            global_length_in_local_units,
        )
        worldline[index] = np.concatenate(
            ([local_time / global_length_in_local_units], position)
        )
        tangent[index] = sample_tangent
        legs[index] = sample_legs
        leg_derivative[index] = sample_derivative
    document: dict[str, object] = {
        "classification": extract.FRAME_SERIES_CLASSIFICATION,
        "times": local_times.tolist(),
        "worldline": worldline.tolist(),
        "worldline_tangent": tangent.tolist(),
        "spatial_legs": legs.tolist(),
        "spatial_leg_derivative": leg_derivative.tolist(),
        "generator": {
            "classification": GENERATOR_CLASSIFICATION,
            "coordinate_system": "Cartesian Kerr-Schild",
            "tetrad_transport": "corotating radial-prograde-tangential-vertical",
            "primary_mass": primary_mass,
            "primary_dimensionless_spin": dimensionless_spin,
            "primary_spin_length": primary_spin,
            "orbital_radius_boyer_lindquist": orbital_radius,
            "orbital_coordinate_radius": coordinate_radius,
            "orbit_direction": direction,
            "coordinate_angular_frequency": omega,
            "coordinate_orbital_speed": coordinate_speed,
            "coordinate_angular_frequency_local_units": (
                omega / global_length_in_local_units
            ),
            "global_length_in_local_units": global_length_in_local_units,
            "source_snapshot_times_global_units": checked_global_times.tolist(),
            "kerr_isco_radius": isco,
            "initial_phase": initial_phase,
            "phase_reference_time": reference_time,
        },
    }
    document["diagnostics"] = _audit_frame(
        document,
        primary_mass,
        primary_spin,
        coordinate_radius,
        omega,
        initial_phase,
        reference_time,
        global_length_in_local_units,
        audit_subsamples,
    )
    return document


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    time_group = parser.add_mutually_exclusive_group(required=True)
    time_group.add_argument("--snapshot-manifest", type=Path)
    time_group.add_argument("--times", type=float, nargs="+")
    parser.add_argument("--primary-mass", type=float, required=True)
    parser.add_argument("--primary-chi", type=float, required=True)
    parser.add_argument("--orbital-radius", type=float, required=True)
    parser.add_argument("--orbit-direction", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--initial-phase", type=float, default=0.0)
    parser.add_argument("--phase-reference-time", type=float)
    parser.add_argument("--audit-subsamples", type=int, default=4)
    parser.add_argument("--global-length-in-local-units", type=float)
    parser.add_argument("--allow-inside-isco", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    if arguments.snapshot_manifest is not None:
        scan = extract.scan_snapshot_manifest(arguments.snapshot_manifest)
        times = [descriptor.time for descriptor in scan.descriptors]
        source_manifest = str(scan.path)
        manifest_scale = float(
            scan.document.get("global_length_in_local_units", 1.0)
        )
        if arguments.global_length_in_local_units is not None and not math.isclose(
            arguments.global_length_in_local_units,
            manifest_scale,
            rel_tol=0.0,
            abs_tol=64.0
            * np.finfo(float).eps
            * max(abs(manifest_scale), 1.0),
        ):
            raise ValueError(
                "command-line length scale disagrees with snapshot manifest"
            )
        length_scale = manifest_scale
    else:
        times = arguments.times
        source_manifest = None
        length_scale = (
            1.0
            if arguments.global_length_in_local_units is None
            else arguments.global_length_in_local_units
        )
    document = build_frame_document(
        times,
        arguments.primary_mass,
        arguments.primary_chi,
        arguments.orbital_radius,
        arguments.orbit_direction,
        arguments.initial_phase,
        arguments.phase_reference_time,
        not arguments.allow_inside_isco,
        arguments.audit_subsamples,
        length_scale,
    )
    if source_manifest is not None:
        document["generator"]["source_snapshot_manifest"] = source_manifest
    output = arguments.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(document["diagnostics"], indent=2, sort_keys=True))
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
