#!/usr/bin/env python3
"""Transport an affine source tetrad along a numerical Kerr ephemeris.

The input supplies global Cartesian Kerr-Schild times, positions, and coordinate
velocities.  Cubic Hermite interpolation provides a continuous worldline and
piecewise-continuous coordinate acceleration.  Spatial tetrad legs are evolved
with Fermi-Walker transport, or with parallel transport after verifying that the
worldline is geodesic to a requested four-acceleration tolerance.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np

import build_kerr_circular_frame as circular
import extract_global_worldtube as extract
import extract_static_taylor_worldtube as static


INPUT_CLASSIFICATION = "athenak-emri-kerr-ephemeris-v1"
GENERATOR_CLASSIFICATION = "athenak-emri-kerr-ephemeris-frame-generator-v1"
TRANSPORT_MODES = ("fermi_walker", "parallel")


def _finite_array(values: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if result.shape != shape or not np.isfinite(result).all():
        raise ValueError(f"{name} must be a finite array with shape {shape}")
    return result


@dataclass(frozen=True)
class HermiteEphemeris:
    global_times: np.ndarray
    positions: np.ndarray
    coordinate_velocities: np.ndarray

    def __post_init__(self) -> None:
        times = np.asarray(self.global_times, dtype=np.float64)
        if (
            times.ndim != 1
            or times.size < 2
            or not np.isfinite(times).all()
            or np.any(np.diff(times) <= 0.0)
        ):
            raise ValueError("ephemeris times must be finite and strictly increasing")
        count = times.size
        positions = _finite_array(self.positions, (count, 3), "ephemeris positions")
        velocities = _finite_array(
            self.coordinate_velocities,
            (count, 3),
            "ephemeris coordinate velocities",
        )
        object.__setattr__(self, "global_times", times)
        object.__setattr__(self, "positions", positions)
        object.__setattr__(self, "coordinate_velocities", velocities)

    def interval(self, time: float) -> int:
        checked = float(time)
        tolerance = 64.0 * np.finfo(float).eps * max(
            1.0, abs(checked), abs(self.global_times[0]), abs(self.global_times[-1])
        )
        if (
            checked < self.global_times[0] - tolerance
            or checked > self.global_times[-1] + tolerance
        ):
            raise ValueError("requested time lies outside the ephemeris")
        checked = min(max(checked, self.global_times[0]), self.global_times[-1])
        interval = int(np.searchsorted(self.global_times, checked, side="right") - 1)
        return min(max(interval, 0), self.global_times.size - 2)

    def evaluate(
        self, time: float, interval: int | None = None
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if interval is None:
            interval = self.interval(time)
        if interval < 0 or interval >= self.global_times.size - 1:
            raise ValueError("ephemeris interval index is out of range")
        left = float(self.global_times[interval])
        right = float(self.global_times[interval + 1])
        tolerance = 64.0 * np.finfo(float).eps * max(
            1.0, abs(left), abs(right), abs(float(time))
        )
        if time < left - tolerance or time > right + tolerance:
            raise ValueError("time is outside the requested ephemeris interval")
        width = right - left
        coordinate = (min(max(float(time), left), right) - left) / width
        s = coordinate
        p0 = self.positions[interval]
        p1 = self.positions[interval + 1]
        v0 = self.coordinate_velocities[interval]
        v1 = self.coordinate_velocities[interval + 1]
        position = (
            (2.0 * s**3 - 3.0 * s**2 + 1.0) * p0
            + (s**3 - 2.0 * s**2 + s) * width * v0
            + (-2.0 * s**3 + 3.0 * s**2) * p1
            + (s**3 - s**2) * width * v1
        )
        velocity = (
            (6.0 * s**2 - 6.0 * s) * p0 / width
            + (3.0 * s**2 - 4.0 * s + 1.0) * v0
            + (-6.0 * s**2 + 6.0 * s) * p1 / width
            + (3.0 * s**2 - 2.0 * s) * v1
        )
        acceleration = (
            (12.0 * s - 6.0) * p0 / width**2
            + (6.0 * s - 4.0) * v0 / width
            + (-12.0 * s + 6.0) * p1 / width**2
            + (6.0 * s - 2.0) * v1 / width
        )
        return position, velocity, acceleration

    def acceleration_jump_diagnostics(self) -> dict[str, float]:
        maximum_absolute = 0.0
        maximum_relative = 0.0
        for knot in range(1, self.global_times.size - 1):
            time = float(self.global_times[knot])
            left = self.evaluate(time, knot - 1)[2]
            right = self.evaluate(time, knot)[2]
            difference = float(np.linalg.norm(right - left))
            scale = max(float(np.linalg.norm(left)), float(np.linalg.norm(right)))
            maximum_absolute = max(maximum_absolute, difference)
            if scale > 0.0:
                maximum_relative = max(maximum_relative, difference / scale)
        return {
            "maximum_coordinate_acceleration_jump": maximum_absolute,
            "maximum_relative_coordinate_acceleration_jump": maximum_relative,
        }


def kerr_connection(
    position: np.ndarray,
    primary_center: np.ndarray,
    primary_mass: float,
    primary_spin: float,
    metric_fd_step: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``g``, ``Gamma^mu_ab``, and ``partial_lambda g_ab``."""

    relative = np.asarray(position, dtype=np.float64) - primary_center
    metric = circular.kerr_schild_metric(relative, primary_mass, primary_spin)
    derivative = np.zeros((4, 4, 4), dtype=np.float64)
    for axis in range(3):
        offset = np.zeros(3)
        offset[axis] = metric_fd_step
        lower = circular.kerr_schild_metric(
            relative - offset, primary_mass, primary_spin
        )
        upper = circular.kerr_schild_metric(
            relative + offset, primary_mass, primary_spin
        )
        derivative[axis + 1] = (upper - lower) / (2.0 * metric_fd_step)
    inverse = np.linalg.inv(metric)
    connection = np.empty((4, 4, 4), dtype=np.float64)
    for mu in range(4):
        for alpha in range(4):
            for beta in range(4):
                value = 0.0
                for nu in range(4):
                    value += 0.5 * inverse[mu, nu] * (
                        derivative[alpha, nu, beta]
                        + derivative[beta, nu, alpha]
                        - derivative[nu, alpha, beta]
                    )
                connection[mu, alpha, beta] = value
    return metric, connection, derivative


@dataclass(frozen=True)
class WorldlineKinematics:
    metric: np.ndarray
    connection: np.ndarray
    coordinate_tangent: np.ndarray
    proper_time_rate: float
    four_velocity: np.ndarray
    four_acceleration: np.ndarray
    proper_acceleration: float
    four_velocity_normalization_error: float
    acceleration_orthogonality_error: float


def worldline_kinematics(
    position: np.ndarray,
    coordinate_velocity: np.ndarray,
    coordinate_acceleration: np.ndarray,
    primary_center: np.ndarray,
    primary_mass: float,
    primary_spin: float,
    metric_fd_step: float,
) -> WorldlineKinematics:
    metric, connection, metric_derivative = kerr_connection(
        position,
        primary_center,
        primary_mass,
        primary_spin,
        metric_fd_step,
    )
    tangent = np.concatenate(([1.0], coordinate_velocity))
    tangent_derivative = np.concatenate(([0.0], coordinate_acceleration))
    norm2 = float(tangent @ metric @ tangent)
    if norm2 >= 0.0:
        raise ValueError("ephemeris worldline is not timelike")
    proper_time_rate = math.sqrt(-norm2)
    metric_along_worldline = np.einsum(
        "l,lmn->mn", tangent, metric_derivative
    )
    norm2_derivative = float(
        tangent @ metric_along_worldline @ tangent
        + 2.0 * tangent_derivative @ metric @ tangent
    )
    proper_time_rate_derivative = -norm2_derivative / (2.0 * proper_time_rate)
    four_velocity = tangent / proper_time_rate
    four_velocity_coordinate_derivative = (
        tangent_derivative / proper_time_rate
        - tangent * proper_time_rate_derivative / proper_time_rate**2
    )
    four_acceleration = (
        four_velocity_coordinate_derivative / proper_time_rate
        + np.einsum("mab,a,b->m", connection, four_velocity, four_velocity)
    )
    acceleration2 = float(four_acceleration @ metric @ four_acceleration)
    roundoff = 512.0 * np.finfo(float).eps * max(
        float(np.max(np.abs(four_acceleration))) ** 2,
        1.0,
    )
    if acceleration2 < -roundoff:
        raise RuntimeError("computed four-acceleration is not spacelike")
    proper_acceleration = math.sqrt(max(acceleration2, 0.0))
    return WorldlineKinematics(
        metric=metric,
        connection=connection,
        coordinate_tangent=tangent,
        proper_time_rate=proper_time_rate,
        four_velocity=four_velocity,
        four_acceleration=four_acceleration,
        proper_acceleration=proper_acceleration,
        four_velocity_normalization_error=abs(
            float(four_velocity @ metric @ four_velocity) + 1.0
        ),
        acceleration_orthogonality_error=abs(
            float(four_velocity @ metric @ four_acceleration)
        ),
    )


def _orthonormalize_spatial_legs(
    metric: np.ndarray, four_velocity: np.ndarray, legs: np.ndarray
) -> tuple[np.ndarray, float]:
    result = np.empty((4, 3), dtype=np.float64)
    for leg in range(3):
        candidate = legs[:, leg].copy()
        candidate += float(candidate @ metric @ four_velocity) * four_velocity
        for previous in range(leg):
            candidate -= (
                float(candidate @ metric @ result[:, previous])
                * result[:, previous]
            )
        norm2 = float(candidate @ metric @ candidate)
        if norm2 <= 0.0:
            raise RuntimeError("transported spatial tetrad became degenerate")
        result[:, leg] = candidate / math.sqrt(norm2)
    correction = float(np.max(np.abs(result - legs)))
    return result, correction


@dataclass
class TransportAudit:
    rhs_evaluation_count: int = 0
    maximum_four_velocity_normalization_error: float = 0.0
    maximum_acceleration_orthogonality_error: float = 0.0
    maximum_dimensionless_proper_acceleration: float = 0.0
    maximum_reorthonormalization_correction: float = 0.0

    def update(self, kinematics: WorldlineKinematics, primary_mass: float) -> None:
        self.rhs_evaluation_count += 1
        self.maximum_four_velocity_normalization_error = max(
            self.maximum_four_velocity_normalization_error,
            kinematics.four_velocity_normalization_error,
        )
        self.maximum_acceleration_orthogonality_error = max(
            self.maximum_acceleration_orthogonality_error,
            kinematics.acceleration_orthogonality_error,
        )
        self.maximum_dimensionless_proper_acceleration = max(
            self.maximum_dimensionless_proper_acceleration,
            primary_mass * kinematics.proper_acceleration,
        )

    def document(self) -> dict[str, float | int]:
        return {
            "transport_rhs_evaluation_count": self.rhs_evaluation_count,
            "maximum_four_velocity_normalization_error": (
                self.maximum_four_velocity_normalization_error
            ),
            "maximum_acceleration_orthogonality_error": (
                self.maximum_acceleration_orthogonality_error
            ),
            "maximum_dimensionless_proper_acceleration": (
                self.maximum_dimensionless_proper_acceleration
            ),
            "maximum_reorthonormalization_correction": (
                self.maximum_reorthonormalization_correction
            ),
        }


@dataclass
class TetradTransporter:
    ephemeris: HermiteEphemeris
    primary_center: np.ndarray
    primary_mass: float
    primary_spin: float
    metric_fd_step: float
    transport_mode: str
    geodesic_acceleration_tolerance: float
    audit: TransportAudit

    def kinematics(
        self, time: float, interval: int | None = None
    ) -> WorldlineKinematics:
        position, velocity, acceleration = self.ephemeris.evaluate(time, interval)
        result = worldline_kinematics(
            position,
            velocity,
            acceleration,
            self.primary_center,
            self.primary_mass,
            self.primary_spin,
            self.metric_fd_step,
        )
        self.audit.update(result, self.primary_mass)
        if (
            self.transport_mode == "parallel"
            and self.primary_mass * result.proper_acceleration
            > self.geodesic_acceleration_tolerance
        ):
            raise ValueError(
                "parallel transport requires a geodesic ephemeris: "
                f"M|a|={self.primary_mass * result.proper_acceleration:.6e} "
                f"exceeds {self.geodesic_acceleration_tolerance:.6e}"
            )
        return result

    def rhs(self, time: float, legs: np.ndarray, interval: int) -> np.ndarray:
        kinematics = self.kinematics(time, interval)
        connection_term = -np.einsum(
            "mab,a,bj->mj",
            kinematics.connection,
            kinematics.coordinate_tangent,
            legs,
        )
        if self.transport_mode == "parallel":
            return connection_term
        acceleration_covector = (
            kinematics.metric @ kinematics.four_acceleration
        )
        velocity_covector = kinematics.metric @ kinematics.four_velocity
        fermi_walker = np.outer(
            kinematics.four_velocity, acceleration_covector
        ) - np.outer(kinematics.four_acceleration, velocity_covector)
        return (
            connection_term
            + kinematics.proper_time_rate * (fermi_walker @ legs)
        )

    def integrate_interval(
        self,
        interval: int,
        initial_legs: np.ndarray,
        substeps: int,
    ) -> np.ndarray:
        left = float(self.ephemeris.global_times[interval])
        right = float(self.ephemeris.global_times[interval + 1])
        step = (right - left) / substeps
        legs = initial_legs.copy()
        for substep in range(substeps):
            time = left + substep * step
            k1 = self.rhs(time, legs, interval)
            k2 = self.rhs(time + 0.5 * step, legs + 0.5 * step * k1, interval)
            k3 = self.rhs(time + 0.5 * step, legs + 0.5 * step * k2, interval)
            k4 = self.rhs(time + step, legs + step * k3, interval)
            candidate = legs + step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            final_kinematics = self.kinematics(time + step, interval)
            legs, correction = _orthonormalize_spatial_legs(
                final_kinematics.metric,
                final_kinematics.four_velocity,
                candidate,
            )
            self.audit.maximum_reorthonormalization_correction = max(
                self.audit.maximum_reorthonormalization_correction,
                correction,
            )
        return legs


def _knot_derivative(
    transporter: TetradTransporter,
    knot: int,
    legs: np.ndarray,
) -> np.ndarray:
    time = float(transporter.ephemeris.global_times[knot])
    if knot == 0:
        return transporter.rhs(time, legs, 0)
    if knot == transporter.ephemeris.global_times.size - 1:
        return transporter.rhs(time, legs, knot - 1)
    left = transporter.rhs(time, legs, knot - 1)
    right = transporter.rhs(time, legs, knot)
    return 0.5 * (left + right)


def _integrate_knots(
    transporter: TetradTransporter,
    initial_legs: np.ndarray,
    substeps: int,
) -> np.ndarray:
    count = transporter.ephemeris.global_times.size
    result = np.empty((count, 4, 3), dtype=np.float64)
    result[0] = initial_legs
    for interval in range(count - 1):
        result[interval + 1] = transporter.integrate_interval(
            interval, result[interval], substeps
        )
    return result


def _audit_frame_document(
    document: dict[str, object],
    ephemeris: HermiteEphemeris,
    primary_center: np.ndarray,
    primary_mass: float,
    primary_spin: float,
    length_scale: float,
    audit_subsamples: int,
) -> dict[str, float | int]:
    series = extract.AffineFrameSeries.from_document(document)
    audit_times = []
    for left, right in zip(series.times[:-1], series.times[1:], strict=True):
        audit_times.extend(np.linspace(left, right, audit_subsamples + 1)[:-1])
    audit_times.append(series.times[-1])
    maximum_position_error = 0.0
    maximum_tangent_error = 0.0
    maximum_gram_error = 0.0
    maximum_time_row_norm = 0.0
    minimum_scaled_determinant = math.inf
    maximum_condition = 0.0
    minimum_timelike_margin = math.inf
    minkowski = np.diag((-1.0, 1.0, 1.0, 1.0))
    for local_time in audit_times:
        global_time = float(local_time) / length_scale
        expected_position, expected_velocity, _ = ephemeris.evaluate(global_time)
        worldline, instantaneous = series.evaluate(float(local_time))
        expected_worldline = np.concatenate(([global_time], expected_position))
        expected_tangent = np.concatenate(([1.0], expected_velocity)) / length_scale
        metric = circular.kerr_schild_metric(
            worldline[1:] - primary_center, primary_mass, primary_spin
        )
        tangent = instantaneous.worldline_tangent
        norm2 = float(tangent @ metric @ tangent)
        if norm2 >= 0.0:
            raise RuntimeError("interpolated ephemeris worldline is not timelike")
        tetrad = np.vstack(
            (tangent / math.sqrt(-norm2), length_scale * instantaneous.spatial_legs.T)
        )
        jacobian = instantaneous.jacobian(np.zeros(3))
        maximum_position_error = max(
            maximum_position_error,
            float(np.max(np.abs(worldline - expected_worldline))),
        )
        maximum_tangent_error = max(
            maximum_tangent_error,
            float(np.max(np.abs(tangent - expected_tangent))),
        )
        maximum_gram_error = max(
            maximum_gram_error,
            float(np.max(np.abs(tetrad @ metric @ tetrad.T - minkowski))),
        )
        maximum_time_row_norm = max(
            maximum_time_row_norm,
            float(np.linalg.norm(instantaneous.spatial_legs[0])),
        )
        minimum_scaled_determinant = min(
            minimum_scaled_determinant,
            abs(float(np.linalg.det(jacobian))) * length_scale**4,
        )
        maximum_condition = max(maximum_condition, float(np.linalg.cond(jacobian)))
        minimum_timelike_margin = min(
            minimum_timelike_margin, -norm2 * length_scale**2
        )
    return {
        "audit_event_count": len(audit_times),
        "audit_subsamples_per_interval": audit_subsamples,
        "maximum_worldline_interpolation_absolute_error": maximum_position_error,
        "maximum_worldline_tangent_interpolation_absolute_error": (
            maximum_tangent_error
        ),
        "maximum_interpolated_tetrad_gram_error": maximum_gram_error,
        "maximum_spatial_leg_time_row_norm": maximum_time_row_norm,
        "minimum_scaled_origin_jacobian_absolute_determinant": (
            minimum_scaled_determinant
        ),
        "maximum_origin_jacobian_condition_number": maximum_condition,
        "minimum_worldline_timelike_margin": minimum_timelike_margin,
    }


def build_frame_document(
    ephemeris: HermiteEphemeris,
    primary_mass: float,
    dimensionless_spin: float,
    primary_center: object = (0.0, 0.0, 0.0),
    disk_normal: object = (0.0, 0.0, 1.0),
    initial_spatial_basis: object | None = None,
    transport_mode: str = "fermi_walker",
    global_length_in_local_units: float = 1.0,
    metric_fd_step: float | None = None,
    integration_substeps_per_interval: int = 16,
    geodesic_acceleration_tolerance: float = 1.0e-8,
    audit_subsamples_per_interval: int = 4,
    run_coarse_convergence_check: bool = True,
) -> dict[str, object]:
    if (
        isinstance(primary_mass, bool)
        or isinstance(dimensionless_spin, bool)
        or isinstance(global_length_in_local_units, bool)
        or isinstance(geodesic_acceleration_tolerance, bool)
        or not math.isfinite(primary_mass)
        or primary_mass <= 0.0
        or not math.isfinite(dimensionless_spin)
        or abs(dimensionless_spin) > 1.0
        or transport_mode not in TRANSPORT_MODES
        or not math.isfinite(global_length_in_local_units)
        or global_length_in_local_units <= 0.0
        or not isinstance(integration_substeps_per_interval, int)
        or isinstance(integration_substeps_per_interval, bool)
        or integration_substeps_per_interval < 1
        or not isinstance(audit_subsamples_per_interval, int)
        or isinstance(audit_subsamples_per_interval, bool)
        or audit_subsamples_per_interval < 1
        or not math.isfinite(geodesic_acceleration_tolerance)
        or geodesic_acceleration_tolerance < 0.0
        or not isinstance(run_coarse_convergence_check, bool)
    ):
        raise ValueError("Kerr ephemeris-frame parameters are invalid")
    center = _finite_array(primary_center, (3,), "primary center")
    normal = _finite_array(disk_normal, (3,), "disk normal")
    normal_norm = float(np.linalg.norm(normal))
    if normal_norm <= 0.0:
        raise ValueError("disk normal cannot vanish")
    normal /= normal_norm
    fd_step = (
        1.0e-5 * primary_mass if metric_fd_step is None else float(metric_fd_step)
    )
    if (
        isinstance(metric_fd_step, bool)
        or not math.isfinite(fd_step)
        or fd_step <= 0.0
    ):
        raise ValueError("metric_fd_step must be finite and positive")
    primary_spin = dimensionless_spin * primary_mass
    first_position = ephemeris.positions[0]
    first_velocity = ephemeris.coordinate_velocities[0]
    if initial_spatial_basis is None:
        spatial_basis = static.canonical_spatial_basis(
            first_position, center, normal
        )
        basis_source = "cylindrical_radial-prograde_tangential-disk_normal"
    else:
        spatial_basis = _finite_array(
            initial_spatial_basis, (3, 3), "initial spatial basis"
        )
        if abs(float(np.linalg.det(spatial_basis))) <= 64.0 * np.finfo(float).eps:
            raise ValueError("initial spatial basis is singular")
        basis_source = "user_supplied_columns"
    initial_metric = circular.kerr_schild_metric(
        first_position - center, primary_mass, primary_spin
    )
    initial_tetrad, _ = static.build_source_tetrad(
        initial_metric, first_velocity, spatial_basis
    )
    audit = TransportAudit()
    transporter = TetradTransporter(
        ephemeris=ephemeris,
        primary_center=center,
        primary_mass=primary_mass,
        primary_spin=primary_spin,
        metric_fd_step=fd_step,
        transport_mode=transport_mode,
        geodesic_acceleration_tolerance=geodesic_acceleration_tolerance,
        audit=audit,
    )
    count = ephemeris.global_times.size
    initial_legs = initial_tetrad[1:].T
    physical_legs = _integrate_knots(
        transporter, initial_legs, integration_substeps_per_interval
    )
    coarse_substeps = max(1, integration_substeps_per_interval // 2)
    coarse_fine_difference = None
    if run_coarse_convergence_check and integration_substeps_per_interval > 1:
        coarse_transporter = TetradTransporter(
            ephemeris=ephemeris,
            primary_center=center,
            primary_mass=primary_mass,
            primary_spin=primary_spin,
            metric_fd_step=fd_step,
            transport_mode=transport_mode,
            geodesic_acceleration_tolerance=geodesic_acceleration_tolerance,
            audit=TransportAudit(),
        )
        coarse_legs = _integrate_knots(
            coarse_transporter, initial_legs, coarse_substeps
        )
        coarse_fine_difference = float(
            np.max(np.abs(physical_legs - coarse_legs))
        )
    physical_derivative = np.empty_like(physical_legs)
    for knot in range(count):
        physical_derivative[knot] = _knot_derivative(
            transporter, knot, physical_legs[knot]
        )
    length_scale = global_length_in_local_units
    local_times = ephemeris.global_times * length_scale
    worldline = np.column_stack((ephemeris.global_times, ephemeris.positions))
    tangent = np.column_stack(
        (np.ones(count), ephemeris.coordinate_velocities)
    ) / length_scale
    legs = physical_legs / length_scale
    leg_derivative = physical_derivative / length_scale**2
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
            "primary_mass": primary_mass,
            "primary_dimensionless_spin": dimensionless_spin,
            "primary_spin_length": primary_spin,
            "primary_spin_axis": [0.0, 0.0, 1.0],
            "primary_center": center.tolist(),
            "disk_normal": normal.tolist(),
            "initial_spatial_basis": spatial_basis.tolist(),
            "initial_spatial_basis_source": basis_source,
            "transport_mode": transport_mode,
            "global_length_in_local_units": length_scale,
            "metric_fd_step_global_units": fd_step,
            "metric_derivative_method": "centered second-order spatial difference",
            "ephemeris_interpolation": "piecewise cubic Hermite",
            "integration_method": "fixed-substep RK4 with Gram projection",
            "integration_substeps_per_interval": (
                integration_substeps_per_interval
            ),
            "coarse_convergence_check_enabled": run_coarse_convergence_check,
            "geodesic_acceleration_tolerance_Ma": (
                geodesic_acceleration_tolerance
            ),
            "source_ephemeris_times_global_units": (
                ephemeris.global_times.tolist()
            ),
        },
    }
    diagnostics: dict[str, object] = {}
    diagnostics.update(audit.document())
    diagnostics["coarse_transport_substeps_per_interval"] = (
        coarse_substeps if coarse_fine_difference is not None else None
    )
    diagnostics["maximum_coarse_fine_spatial_leg_difference"] = (
        coarse_fine_difference
    )
    diagnostics.update(ephemeris.acceleration_jump_diagnostics())
    diagnostics.update(
        _audit_frame_document(
            document,
            ephemeris,
            center,
            primary_mass,
            primary_spin,
            length_scale,
            audit_subsamples_per_interval,
        )
    )
    document["diagnostics"] = diagnostics
    return document


def read_ephemeris_document(path: Path) -> tuple[dict[str, object], HermiteEphemeris]:
    resolved = path.expanduser().resolve(strict=True)
    document = json.loads(resolved.read_text(encoding="utf-8"))
    if (
        not isinstance(document, dict)
        or document.get("classification") != INPUT_CLASSIFICATION
    ):
        raise ValueError("Kerr ephemeris classification is unsupported")
    times = np.asarray(document.get("global_times"), dtype=np.float64)
    if times.ndim != 1:
        raise ValueError("ephemeris global_times must be one-dimensional")
    ephemeris = HermiteEphemeris(
        times,
        _finite_array(
            document.get("positions"), (times.size, 3), "ephemeris positions"
        ),
        _finite_array(
            document.get("coordinate_velocities"),
            (times.size, 3),
            "ephemeris coordinate velocities",
        ),
    )
    return document, ephemeris


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ephemeris", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def _document_float(
    document: dict[str, object], key: str, default: float | None = None
) -> float:
    if key not in document:
        if default is None:
            raise ValueError(f"ephemeris is missing numeric field {key}")
        return default
    value = document[key]
    if isinstance(value, bool):
        raise ValueError(f"ephemeris field {key} must be numeric")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"ephemeris field {key} must be numeric") from error
    if not math.isfinite(result):
        raise ValueError(f"ephemeris field {key} must be finite")
    return result


def _document_int(
    document: dict[str, object], key: str, default: int
) -> int:
    value = document.get(key, default)
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f"ephemeris field {key} must be an integer")
    return value


def _document_bool(
    document: dict[str, object], key: str, default: bool
) -> bool:
    value = document.get(key, default)
    if not isinstance(value, bool):
        raise ValueError(f"ephemeris field {key} must be true or false")
    return value


def main() -> int:
    arguments = parse_arguments()
    source_path = arguments.ephemeris.expanduser().resolve(strict=True)
    source, ephemeris = read_ephemeris_document(source_path)
    document = build_frame_document(
        ephemeris,
        primary_mass=_document_float(source, "primary_mass"),
        dimensionless_spin=_document_float(
            source, "primary_dimensionless_spin", 0.0
        ),
        primary_center=source.get("primary_center", (0.0, 0.0, 0.0)),
        disk_normal=source.get("disk_normal", (0.0, 0.0, 1.0)),
        initial_spatial_basis=source.get("initial_spatial_basis"),
        transport_mode=str(source.get("transport_mode", "fermi_walker")),
        global_length_in_local_units=_document_float(
            source, "global_length_in_local_units", 1.0
        ),
        metric_fd_step=(
            None
            if source.get("metric_fd_step_global_units") is None
            else _document_float(source, "metric_fd_step_global_units")
        ),
        integration_substeps_per_interval=_document_int(
            source, "integration_substeps_per_interval", 16
        ),
        geodesic_acceleration_tolerance=_document_float(
            source, "geodesic_acceleration_tolerance_Ma", 1.0e-8
        ),
        audit_subsamples_per_interval=_document_int(
            source, "audit_subsamples_per_interval", 4
        ),
        run_coarse_convergence_check=_document_bool(
            source, "run_coarse_convergence_check", True
        ),
    )
    document["generator"]["source_ephemeris"] = str(source_path)
    document["generator"]["source_ephemeris_sha256"] = static.file_sha256(
        source_path
    )
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
