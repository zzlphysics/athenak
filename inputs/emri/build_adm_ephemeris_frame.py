#!/usr/bin/env python3
"""Transport an affine ephemeris frame in a numerical ADM snapshot series.

The same fixed-level, trilinear-space and linear-time interpolation contract as
``extract_global_worldtube.py`` is used.  Spatial metric derivatives are centered
finite differences.  The time derivative is analytic for the piecewise-linear
ADM interpolation.  Transport intervals are split at both ephemeris and source
snapshot knots so no Runge--Kutta step silently crosses a derivative jump.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np

import build_kerr_ephemeris_frame as ephemeris_frame
import extract_global_worldtube as extract
import extract_static_taylor_worldtube as static


GENERATOR_CLASSIFICATION = "athenak-emri-adm-ephemeris-frame-generator-v1"


def _adm_geometry(values: object) -> tuple[np.ndarray, float, np.ndarray]:
    source = ephemeris_frame._finite_array(
        values, (len(extract.ADM_VARIABLES),), "interpolated ADM state"
    )
    gamma = np.asarray(
        (
            (source[0], source[1], source[2]),
            (source[1], source[3], source[4]),
            (source[2], source[4], source[5]),
        )
    )
    alpha = float(source[6])
    beta = source[7:10]
    if (
        alpha <= 0.0
        or float(np.linalg.det(gamma)) <= 0.0
        or float(np.min(np.linalg.eigvalsh(gamma))) <= 0.0
    ):
        raise RuntimeError("interpolated ADM geometry is invalid")
    return gamma, alpha, beta


def adm_metric(values: object) -> np.ndarray:
    gamma, alpha, beta = _adm_geometry(values)
    return static.spacetime_metric_from_adm(gamma, alpha, beta)


def adm_metric_time_derivative(values: object, derivatives: object) -> np.ndarray:
    source = ephemeris_frame._finite_array(
        values, (len(extract.ADM_VARIABLES),), "interpolated ADM state"
    )
    change = ephemeris_frame._finite_array(
        derivatives, (len(extract.ADM_VARIABLES),), "ADM time derivative"
    )
    gamma, alpha, beta = _adm_geometry(source)
    gamma_change = np.asarray(
        (
            (change[0], change[1], change[2]),
            (change[1], change[3], change[4]),
            (change[2], change[4], change[5]),
        )
    )
    alpha_change = float(change[6])
    beta_change = change[7:10]
    beta_lower = gamma @ beta
    beta_lower_change = gamma_change @ beta + gamma @ beta_change
    result = np.zeros((4, 4), dtype=np.float64)
    result[1:, 1:] = gamma_change
    result[0, 1:] = beta_lower_change
    result[1:, 0] = beta_lower_change
    result[0, 0] = (
        -2.0 * alpha * alpha_change
        + 2.0 * float(beta_change @ beta_lower)
        + float(beta @ gamma_change @ beta)
    )
    return result


def connection_from_metric_derivative(
    metric: np.ndarray, derivative: np.ndarray
) -> np.ndarray:
    inverse = np.linalg.inv(metric)
    connection = np.empty((4, 4, 4), dtype=np.float64)
    for mu in range(4):
        for alpha in range(4):
            for beta in range(4):
                connection[mu, alpha, beta] = 0.5 * sum(
                    inverse[mu, nu]
                    * (
                        derivative[alpha, nu, beta]
                        + derivative[beta, nu, alpha]
                        - derivative[nu, alpha, beta]
                    )
                    for nu in range(4)
                )
    return connection


@dataclass
class GeometryAudit:
    evaluation_count: int = 0
    sampled_event_count: int = 0
    minimum_lapse: float = math.inf
    minimum_spatial_metric_eigenvalue: float = math.inf
    maximum_connection_lower_index_asymmetry: float = 0.0
    maximum_metric_compatibility_residual: float = 0.0

    def update(
        self,
        adm_values: np.ndarray,
        metric: np.ndarray,
        derivative: np.ndarray,
        connection: np.ndarray,
        event_count: int,
    ) -> None:
        gamma, alpha, _ = _adm_geometry(adm_values)
        residual = np.empty((4, 4, 4), dtype=np.float64)
        for lam in range(4):
            for mu in range(4):
                for nu in range(4):
                    residual[lam, mu, nu] = (
                        derivative[lam, mu, nu]
                        - connection[:, lam, mu] @ metric[:, nu]
                        - connection[:, lam, nu] @ metric[mu, :]
                    )
        self.evaluation_count += 1
        self.sampled_event_count += event_count
        self.minimum_lapse = min(self.minimum_lapse, alpha)
        self.minimum_spatial_metric_eigenvalue = min(
            self.minimum_spatial_metric_eigenvalue,
            float(np.min(np.linalg.eigvalsh(gamma))),
        )
        self.maximum_connection_lower_index_asymmetry = max(
            self.maximum_connection_lower_index_asymmetry,
            float(np.max(np.abs(connection - connection.transpose(0, 2, 1)))),
        )
        self.maximum_metric_compatibility_residual = max(
            self.maximum_metric_compatibility_residual,
            float(np.max(np.abs(residual))),
        )

    def document(self) -> dict[str, float | int]:
        return {
            "adm_geometry_evaluation_count": self.evaluation_count,
            "adm_sampled_spacetime_event_count": self.sampled_event_count,
            "minimum_interpolated_adm_lapse": self.minimum_lapse,
            "minimum_interpolated_adm_spatial_metric_eigenvalue": (
                self.minimum_spatial_metric_eigenvalue
            ),
            "maximum_connection_lower_index_asymmetry": (
                self.maximum_connection_lower_index_asymmetry
            ),
            "maximum_metric_compatibility_residual": (
                self.maximum_metric_compatibility_residual
            ),
        }


@dataclass
class ADMSnapshotGeometry:
    snapshots: extract.SourceSeries
    spatial_fd_step: float
    audit: GeometryAudit

    def __post_init__(self) -> None:
        if not math.isfinite(self.spatial_fd_step) or self.spatial_fd_step <= 0.0:
            raise ValueError("ADM metric finite-difference step must be positive")

    def metric_connection(
        self,
        time: float,
        position: object,
        source_interval: int,
        spatial_fd_step: float | None = None,
        record_audit: bool = True,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        center = ephemeris_frame._finite_array(position, (3,), "ADM sample position")
        step = self.spatial_fd_step if spatial_fd_step is None else spatial_fd_step
        if not math.isfinite(step) or step <= 0.0:
            raise ValueError("ADM metric finite-difference step must be positive")
        positions = [center]
        for axis in range(3):
            offset = np.zeros(3)
            offset[axis] = step
            positions.extend((center - offset, center + offset))
        events = np.column_stack(
            (
                np.full(len(positions), float(time)),
                np.asarray(positions),
            )
        )
        intervals = np.full(len(positions), source_interval, dtype=np.int64)
        values, time_derivatives = self.snapshots.sample_with_time_derivative(
            events, intervals
        )
        if values.shape != (7, len(extract.ADM_VARIABLES)):
            raise RuntimeError("ADM-only snapshot series returned an invalid field set")
        metrics = np.asarray([adm_metric(row) for row in values])
        derivative = np.empty((4, 4, 4), dtype=np.float64)
        derivative[0] = adm_metric_time_derivative(
            values[0], time_derivatives[0]
        )
        for axis in range(3):
            derivative[axis + 1] = (
                metrics[2 * axis + 2] - metrics[2 * axis + 1]
            ) / (2.0 * step)
        connection = connection_from_metric_derivative(metrics[0], derivative)
        if record_audit:
            self.audit.update(
                values[0], metrics[0], derivative, connection, len(positions)
            )
        return metrics[0], connection, derivative

    def metric(
        self, time: float, position: object, source_interval: int
    ) -> np.ndarray:
        center = ephemeris_frame._finite_array(position, (3,), "ADM sample position")
        event = np.asarray(((float(time), *center),))
        values, _ = self.snapshots.sample_with_time_derivative(
            event, np.asarray((source_interval,))
        )
        if values.shape != (1, len(extract.ADM_VARIABLES)):
            raise RuntimeError("ADM-only snapshot series returned an invalid field set")
        return adm_metric(values[0])


def worldline_kinematics(
    geometry: ADMSnapshotGeometry,
    ephemeris: ephemeris_frame.HermiteEphemeris,
    time: float,
    ephemeris_interval: int,
    source_interval: int,
) -> ephemeris_frame.WorldlineKinematics:
    position, velocity, acceleration = ephemeris.evaluate(
        time, ephemeris_interval
    )
    metric, connection, metric_derivative = geometry.metric_connection(
        time, position, source_interval
    )
    tangent = np.concatenate(([1.0], velocity))
    tangent_derivative = np.concatenate(([0.0], acceleration))
    norm2 = float(tangent @ metric @ tangent)
    if norm2 >= 0.0:
        raise ValueError("ephemeris worldline is not timelike in the ADM snapshots")
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
    coordinate_derivative = (
        tangent_derivative / proper_time_rate
        - tangent * proper_time_rate_derivative / proper_time_rate**2
    )
    four_acceleration = (
        coordinate_derivative / proper_time_rate
        + np.einsum("mab,a,b->m", connection, four_velocity, four_velocity)
    )
    acceleration2 = float(four_acceleration @ metric @ four_acceleration)
    roundoff = 512.0 * np.finfo(float).eps * max(
        float(np.max(np.abs(four_acceleration))) ** 2, 1.0
    )
    if acceleration2 < -roundoff:
        raise RuntimeError("computed ADM four-acceleration is not spacelike")
    return ephemeris_frame.WorldlineKinematics(
        metric=metric,
        connection=connection,
        coordinate_tangent=tangent,
        proper_time_rate=proper_time_rate,
        four_velocity=four_velocity,
        four_acceleration=four_acceleration,
        proper_acceleration=math.sqrt(max(acceleration2, 0.0)),
        four_velocity_normalization_error=abs(
            float(four_velocity @ metric @ four_velocity) + 1.0
        ),
        acceleration_orthogonality_error=abs(
            float(four_velocity @ metric @ four_acceleration)
        ),
    )


@dataclass(frozen=True)
class TransportSegment:
    left: float
    right: float
    ephemeris_interval: int
    source_interval: int


def _merged_transport_knots(
    ephemeris: ephemeris_frame.HermiteEphemeris, source_times: object
) -> np.ndarray:
    snapshots = np.asarray(source_times, dtype=np.float64)
    if (
        snapshots.ndim != 1
        or snapshots.size < 2
        or not np.isfinite(snapshots).all()
        or np.any(np.diff(snapshots) <= 0.0)
    ):
        raise ValueError("ADM snapshot times must be finite and strictly increasing")
    start = float(ephemeris.global_times[0])
    stop = float(ephemeris.global_times[-1])
    tolerance = 128.0 * np.finfo(float).eps * max(1.0, abs(start), abs(stop))
    if start < snapshots[0] - tolerance or stop > snapshots[-1] + tolerance:
        raise ValueError("ephemeris time range lies outside the ADM snapshot series")
    candidates = np.concatenate(
        (
            ephemeris.global_times,
            snapshots[(snapshots > start + tolerance) & (snapshots < stop - tolerance)],
        )
    )
    candidates.sort()
    merged = [float(candidates[0])]
    for candidate in candidates[1:]:
        local_tolerance = 128.0 * np.finfo(float).eps * max(
            1.0, abs(float(candidate)), abs(merged[-1])
        )
        if float(candidate) - merged[-1] > local_tolerance:
            merged.append(float(candidate))
    return np.asarray(merged)


def _transport_segments(
    knots: np.ndarray,
    ephemeris: ephemeris_frame.HermiteEphemeris,
    source_times: np.ndarray,
) -> tuple[TransportSegment, ...]:
    result = []
    for left, right in zip(knots[:-1], knots[1:], strict=True):
        midpoint = 0.5 * (float(left) + float(right))
        ephemeris_interval = ephemeris.interval(midpoint)
        source_interval = int(np.searchsorted(source_times, midpoint) - 1)
        source_interval = min(max(source_interval, 0), source_times.size - 2)
        result.append(
            TransportSegment(
                float(left), float(right), ephemeris_interval, source_interval
            )
        )
    return tuple(result)


@dataclass
class ADMTetradTransporter:
    ephemeris: ephemeris_frame.HermiteEphemeris
    geometry: ADMSnapshotGeometry
    transport_mode: str
    acceleration_scale: float
    geodesic_acceleration_tolerance: float
    audit: ephemeris_frame.TransportAudit

    def kinematics(
        self, time: float, segment: TransportSegment
    ) -> ephemeris_frame.WorldlineKinematics:
        result = worldline_kinematics(
            self.geometry,
            self.ephemeris,
            time,
            segment.ephemeris_interval,
            segment.source_interval,
        )
        self.audit.update(result, self.acceleration_scale)
        scaled_acceleration = self.acceleration_scale * result.proper_acceleration
        if (
            self.transport_mode == "parallel"
            and scaled_acceleration > self.geodesic_acceleration_tolerance
        ):
            raise ValueError(
                "parallel transport requires a geodesic ephemeris in the sampled "
                f"ADM geometry: scale*|a|={scaled_acceleration:.6e} exceeds "
                f"{self.geodesic_acceleration_tolerance:.6e}"
            )
        return result

    def rhs(
        self, time: float, legs: np.ndarray, segment: TransportSegment
    ) -> np.ndarray:
        kinematics = self.kinematics(time, segment)
        connection_term = -np.einsum(
            "mab,a,bj->mj",
            kinematics.connection,
            kinematics.coordinate_tangent,
            legs,
        )
        if self.transport_mode == "parallel":
            return connection_term
        acceleration_covector = kinematics.metric @ kinematics.four_acceleration
        velocity_covector = kinematics.metric @ kinematics.four_velocity
        fermi_walker = np.outer(
            kinematics.four_velocity, acceleration_covector
        ) - np.outer(kinematics.four_acceleration, velocity_covector)
        return (
            connection_term
            + kinematics.proper_time_rate * (fermi_walker @ legs)
        )

    def integrate_segment(
        self,
        segment: TransportSegment,
        initial_legs: np.ndarray,
        substeps: int,
    ) -> np.ndarray:
        step = (segment.right - segment.left) / substeps
        legs = initial_legs.copy()
        for substep in range(substeps):
            time = segment.left + substep * step
            k1 = self.rhs(time, legs, segment)
            k2 = self.rhs(
                time + 0.5 * step, legs + 0.5 * step * k1, segment
            )
            k3 = self.rhs(
                time + 0.5 * step, legs + 0.5 * step * k2, segment
            )
            k4 = self.rhs(time + step, legs + step * k3, segment)
            candidate = legs + step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
            final_kinematics = self.kinematics(time + step, segment)
            legs, correction = ephemeris_frame._orthonormalize_spatial_legs(
                final_kinematics.metric,
                final_kinematics.four_velocity,
                candidate,
            )
            self.audit.maximum_reorthonormalization_correction = max(
                self.audit.maximum_reorthonormalization_correction, correction
            )
        return legs


def _transport_and_difference_diagnostics(
    transporter: ADMTetradTransporter,
    coarse_transporter: ADMTetradTransporter | None,
    segments: tuple[TransportSegment, ...],
    initial_legs: np.ndarray,
    fine_substeps: int,
    coarse_substeps: int,
    run_fd_check: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None, dict[str, object]]:
    """Transport chronologically so a two-snapshot cache is sufficient."""

    count = len(segments) + 1
    fine_legs = np.empty((count, 4, 3), dtype=np.float64)
    fine_legs[0] = initial_legs
    coarse_legs = None
    if coarse_transporter is not None:
        coarse_legs = np.empty_like(fine_legs)
        coarse_legs[0] = initial_legs
    left_derivatives = np.empty((len(segments), 4, 3), dtype=np.float64)
    right_derivatives = np.empty_like(left_derivatives)
    maximum_fd_absolute = 0.0
    maximum_fd_relative = 0.0
    maximum_jump_absolute = 0.0
    maximum_jump_relative = 0.0
    fd_count = 0
    jump_count = 0
    geometry = transporter.geometry
    ephemeris = transporter.ephemeris
    for index, segment in enumerate(segments):
        if (
            index > 0
            and segments[index - 1].source_interval != segment.source_interval
        ):
            position = ephemeris.evaluate(segment.left)[0]
            _, left_connection, _ = geometry.metric_connection(
                segment.left,
                position,
                segments[index - 1].source_interval,
                record_audit=False,
            )
            _, right_connection, _ = geometry.metric_connection(
                segment.left,
                position,
                segment.source_interval,
                record_audit=False,
            )
            difference = float(
                np.max(np.abs(right_connection - left_connection))
            )
            scale = max(
                float(np.max(np.abs(left_connection))),
                float(np.max(np.abs(right_connection))),
            )
            maximum_jump_absolute = max(maximum_jump_absolute, difference)
            if scale > 0.0:
                maximum_jump_relative = max(
                    maximum_jump_relative, difference / scale
                )
            jump_count += 1
        if run_fd_check:
            time = 0.5 * (segment.left + segment.right)
            position = ephemeris.evaluate(time, segment.ephemeris_interval)[0]
            _, connection, _ = geometry.metric_connection(
                time, position, segment.source_interval, record_audit=False
            )
            _, refined, _ = geometry.metric_connection(
                time,
                position,
                segment.source_interval,
                spatial_fd_step=0.5 * geometry.spatial_fd_step,
                record_audit=False,
            )
            difference = float(np.max(np.abs(refined - connection)))
            scale = max(
                float(np.max(np.abs(refined))),
                float(np.max(np.abs(connection))),
            )
            maximum_fd_absolute = max(maximum_fd_absolute, difference)
            if scale > 0.0:
                maximum_fd_relative = max(
                    maximum_fd_relative, difference / scale
                )
            fd_count += 1
        left_derivatives[index] = transporter.rhs(
            segment.left, fine_legs[index], segment
        )
        fine_legs[index + 1] = transporter.integrate_segment(
            segment, fine_legs[index], fine_substeps
        )
        right_derivatives[index] = transporter.rhs(
            segment.right, fine_legs[index + 1], segment
        )
        if coarse_transporter is not None and coarse_legs is not None:
            coarse_legs[index + 1] = coarse_transporter.integrate_segment(
                segment, coarse_legs[index], coarse_substeps
            )
    derivatives = np.empty_like(fine_legs)
    derivatives[0] = left_derivatives[0]
    derivatives[-1] = right_derivatives[-1]
    for knot in range(1, len(segments)):
        derivatives[knot] = 0.5 * (
            right_derivatives[knot - 1] + left_derivatives[knot]
        )
    diagnostics: dict[str, object] = {
        "metric_fd_convergence_check_enabled": run_fd_check,
        "metric_fd_convergence_sample_count": fd_count,
        "maximum_half_step_connection_absolute_difference": (
            maximum_fd_absolute if run_fd_check else None
        ),
        "maximum_half_step_connection_relative_difference": (
            maximum_fd_relative if run_fd_check else None
        ),
        "source_time_derivative_jump_count": jump_count,
        "maximum_source_knot_connection_absolute_jump": maximum_jump_absolute,
        "maximum_source_knot_connection_relative_jump": maximum_jump_relative,
    }
    return fine_legs, derivatives, coarse_legs, diagnostics


def _audit_frame_document(
    document: dict[str, object],
    ephemeris: ephemeris_frame.HermiteEphemeris,
    geometry: ADMSnapshotGeometry,
    segments: tuple[TransportSegment, ...],
    length_scale: float,
    audit_subsamples: int,
) -> dict[str, float | int]:
    series = extract.AffineFrameSeries.from_document(document)
    maximum_position_error = 0.0
    maximum_tangent_error = 0.0
    maximum_gram_error = 0.0
    maximum_time_row_norm = 0.0
    minimum_scaled_determinant = math.inf
    maximum_condition = 0.0
    minimum_timelike_margin = math.inf
    event_count = 0
    minkowski = np.diag((-1.0, 1.0, 1.0, 1.0))
    for segment in segments:
        times = np.linspace(
            segment.left, segment.right, audit_subsamples + 1
        )[:-1]
        for global_time in times:
            expected_position, expected_velocity, _ = ephemeris.evaluate(
                float(global_time), segment.ephemeris_interval
            )
            local_time = float(global_time) * length_scale
            worldline, instantaneous = series.evaluate(local_time)
            expected_worldline = np.concatenate(
                ([float(global_time)], expected_position)
            )
            expected_tangent = np.concatenate(
                ([1.0], expected_velocity)
            ) / length_scale
            metric = geometry.metric(
                float(global_time), worldline[1:], segment.source_interval
            )
            tangent = instantaneous.worldline_tangent
            norm2 = float(tangent @ metric @ tangent)
            if norm2 >= 0.0:
                raise RuntimeError("interpolated frame worldline is not timelike")
            tetrad = np.vstack(
                (
                    tangent / math.sqrt(-norm2),
                    length_scale * instantaneous.spatial_legs.T,
                )
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
            maximum_condition = max(
                maximum_condition, float(np.linalg.cond(jacobian))
            )
            minimum_timelike_margin = min(
                minimum_timelike_margin, -norm2 * length_scale**2
            )
            event_count += 1

    final_segment = segments[-1]
    final_global_time = final_segment.right
    final_local_time = final_global_time * length_scale
    expected_position, expected_velocity, _ = ephemeris.evaluate(
        final_global_time, final_segment.ephemeris_interval
    )
    worldline, instantaneous = series.evaluate(final_local_time)
    metric = geometry.metric(
        final_global_time, worldline[1:], final_segment.source_interval
    )
    tangent = instantaneous.worldline_tangent
    norm2 = float(tangent @ metric @ tangent)
    if norm2 >= 0.0:
        raise RuntimeError("final interpolated frame worldline is not timelike")
    tetrad = np.vstack(
        (
            tangent / math.sqrt(-norm2),
            length_scale * instantaneous.spatial_legs.T,
        )
    )
    jacobian = instantaneous.jacobian(np.zeros(3))
    maximum_position_error = max(
        maximum_position_error,
        float(
            np.max(
                np.abs(
                    worldline
                    - np.concatenate(([final_global_time], expected_position))
                )
            )
        ),
    )
    maximum_tangent_error = max(
        maximum_tangent_error,
        float(
            np.max(
                np.abs(
                    tangent
                    - np.concatenate(([1.0], expected_velocity)) / length_scale
                )
            )
        ),
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
    event_count += 1
    return {
        "audit_event_count": event_count,
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
    ephemeris: ephemeris_frame.HermiteEphemeris,
    snapshots: extract.SourceSeries,
    primary_center: object = (0.0, 0.0, 0.0),
    disk_normal: object = (0.0, 0.0, 1.0),
    initial_spatial_basis: object | None = None,
    transport_mode: str = "fermi_walker",
    global_length_in_local_units: float = 1.0,
    acceleration_scale_global_units: float = 1.0,
    metric_fd_step: float = 1.0e-3,
    integration_substeps_per_interval: int = 16,
    geodesic_acceleration_tolerance: float = 1.0e-8,
    audit_subsamples_per_interval: int = 4,
    run_coarse_convergence_check: bool = True,
    run_metric_fd_convergence_check: bool = True,
) -> dict[str, object]:
    scalar_values = (
        global_length_in_local_units,
        acceleration_scale_global_units,
        metric_fd_step,
        geodesic_acceleration_tolerance,
    )
    if (
        any(isinstance(value, bool) for value in scalar_values)
        or not all(math.isfinite(float(value)) for value in scalar_values)
        or global_length_in_local_units <= 0.0
        or acceleration_scale_global_units <= 0.0
        or metric_fd_step <= 0.0
        or geodesic_acceleration_tolerance < 0.0
        or transport_mode not in ephemeris_frame.TRANSPORT_MODES
        or not isinstance(integration_substeps_per_interval, int)
        or isinstance(integration_substeps_per_interval, bool)
        or integration_substeps_per_interval < 1
        or not isinstance(audit_subsamples_per_interval, int)
        or isinstance(audit_subsamples_per_interval, bool)
        or audit_subsamples_per_interval < 1
        or not isinstance(run_coarse_convergence_check, bool)
        or not isinstance(run_metric_fd_convergence_check, bool)
    ):
        raise ValueError("ADM ephemeris-frame parameters are invalid")
    center = ephemeris_frame._finite_array(primary_center, (3,), "primary center")
    normal = ephemeris_frame._finite_array(disk_normal, (3,), "disk normal")
    normal_norm = float(np.linalg.norm(normal))
    if normal_norm <= 0.0:
        raise ValueError("disk normal cannot vanish")
    normal /= normal_norm
    source_times = np.asarray(snapshots.times, dtype=np.float64)
    knots = _merged_transport_knots(ephemeris, source_times)
    segments = _transport_segments(knots, ephemeris, source_times)
    if not segments:
        raise ValueError("ADM transport requires at least one time interval")
    geometry_audit = GeometryAudit()
    geometry = ADMSnapshotGeometry(snapshots, metric_fd_step, geometry_audit)
    first_position, first_velocity, _ = ephemeris.evaluate(
        float(knots[0]), segments[0].ephemeris_interval
    )
    if initial_spatial_basis is None:
        spatial_basis = static.canonical_spatial_basis(
            first_position, center, normal
        )
        basis_source = "cylindrical_radial-prograde_tangential-disk_normal"
    else:
        spatial_basis = ephemeris_frame._finite_array(
            initial_spatial_basis, (3, 3), "initial spatial basis"
        )
        if abs(float(np.linalg.det(spatial_basis))) <= 64.0 * np.finfo(float).eps:
            raise ValueError("initial spatial basis is singular")
        basis_source = "user_supplied_columns"
    initial_metric, _, _ = geometry.metric_connection(
        float(knots[0]), first_position, segments[0].source_interval
    )
    initial_tetrad, _ = static.build_source_tetrad(
        initial_metric, first_velocity, spatial_basis
    )
    transport_audit = ephemeris_frame.TransportAudit()
    transporter = ADMTetradTransporter(
        ephemeris=ephemeris,
        geometry=geometry,
        transport_mode=transport_mode,
        acceleration_scale=acceleration_scale_global_units,
        geodesic_acceleration_tolerance=geodesic_acceleration_tolerance,
        audit=transport_audit,
    )
    initial_legs = initial_tetrad[1:].T
    coarse_substeps = max(1, integration_substeps_per_interval // 2)
    coarse_transporter = None
    if run_coarse_convergence_check and integration_substeps_per_interval > 1:
        coarse_transporter = ADMTetradTransporter(
            ephemeris=ephemeris,
            geometry=geometry,
            transport_mode=transport_mode,
            acceleration_scale=acceleration_scale_global_units,
            geodesic_acceleration_tolerance=geodesic_acceleration_tolerance,
            audit=ephemeris_frame.TransportAudit(),
        )
    (
        physical_legs,
        physical_derivative,
        coarse_legs,
        difference_diagnostics,
    ) = _transport_and_difference_diagnostics(
        transporter,
        coarse_transporter,
        segments,
        initial_legs,
        integration_substeps_per_interval,
        coarse_substeps,
        run_metric_fd_convergence_check,
    )
    coarse_fine_difference = None
    if coarse_legs is not None:
        coarse_fine_difference = float(
            np.max(np.abs(physical_legs - coarse_legs))
        )
    positions = []
    velocities = []
    for time in knots:
        position, velocity, _ = ephemeris.evaluate(float(time))
        positions.append(position)
        velocities.append(velocity)
    positions_array = np.asarray(positions)
    velocities_array = np.asarray(velocities)
    length_scale = global_length_in_local_units
    local_times = knots * length_scale
    worldline = np.column_stack((knots, positions_array))
    tangent = np.column_stack((np.ones(knots.size), velocities_array)) / length_scale
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
            "coordinate_system": "source ADM Cartesian coordinates",
            "primary_center": center.tolist(),
            "disk_normal": normal.tolist(),
            "initial_spatial_basis": spatial_basis.tolist(),
            "initial_spatial_basis_source": basis_source,
            "transport_mode": transport_mode,
            "global_length_in_local_units": length_scale,
            "acceleration_scale_global_units": acceleration_scale_global_units,
            "metric_fd_step_global_units": metric_fd_step,
            "metric_spatial_derivative_method": "centered second-order difference",
            "metric_time_derivative_method": (
                "analytic derivative of piecewise-linear ADM interpolation"
            ),
            "ephemeris_interpolation": "piecewise cubic Hermite",
            "integration_method": "fixed-substep RK4 with Gram projection",
            "integration_substeps_per_merged_interval": (
                integration_substeps_per_interval
            ),
            "coarse_convergence_check_enabled": run_coarse_convergence_check,
            "metric_fd_convergence_check_enabled": (
                run_metric_fd_convergence_check
            ),
            "geodesic_acceleration_tolerance_scaled": (
                geodesic_acceleration_tolerance
            ),
            "source_snapshot_times_global_units": source_times.tolist(),
            "source_ephemeris_times_global_units": (
                ephemeris.global_times.tolist()
            ),
            "merged_transport_times_global_units": knots.tolist(),
        },
    }
    diagnostics: dict[str, object] = {}
    diagnostics.update(transport_audit.document())
    diagnostics.update(geometry_audit.document())
    diagnostics["merged_transport_interval_count"] = len(segments)
    diagnostics["inserted_source_snapshot_knot_count"] = int(
        knots.size - ephemeris.global_times.size
    )
    diagnostics["coarse_transport_substeps_per_merged_interval"] = (
        coarse_substeps if coarse_fine_difference is not None else None
    )
    diagnostics["maximum_coarse_fine_spatial_leg_difference"] = (
        coarse_fine_difference
    )
    diagnostics.update(ephemeris.acceleration_jump_diagnostics())
    diagnostics.update(difference_diagnostics)
    diagnostics.update(
        _audit_frame_document(
            document,
            ephemeris,
            geometry,
            segments,
            length_scale,
            audit_subsamples_per_interval,
        )
    )
    document["diagnostics"] = diagnostics
    return document


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--ephemeris", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def _default_acceleration_scale(document: dict[str, object]) -> float:
    if "acceleration_scale_global_units" in document:
        return ephemeris_frame._document_float(
            document, "acceleration_scale_global_units"
        )
    return ephemeris_frame._document_float(document, "primary_mass", 1.0)


def build_frame_from_scan(
    scan: extract.SnapshotManifestScan,
    ephemeris_path: Path,
    metric_fd_step_override: float | None = None,
    hash_source_files_override: bool | None = None,
) -> dict[str, object]:
    ephemeris_path = ephemeris_path.expanduser().resolve(strict=True)
    source, ephemeris = ephemeris_frame.read_ephemeris_document(ephemeris_path)

    def load_snapshot(index: int) -> extract.Snapshot:
        return extract._load_adm_snapshot(
            scan.entries[index], scan.path.parent, scan.source_level
        )

    snapshots = extract.LazySnapshotSeries(
        scan.descriptors, load_snapshot, cache_size=scan.snapshot_cache_size
    )
    default_fd_step = 0.25 * min(
        float(np.min(descriptor.spacing)) for descriptor in scan.descriptors
    )
    metric_fd_step = (
        ephemeris_frame._document_float(
            source, "adm_metric_fd_step_global_units", default_fd_step
        )
        if metric_fd_step_override is None
        else float(metric_fd_step_override)
    )
    document = build_frame_document(
        ephemeris,
        snapshots,
        primary_center=source.get("primary_center", (0.0, 0.0, 0.0)),
        disk_normal=source.get("disk_normal", (0.0, 0.0, 1.0)),
        initial_spatial_basis=source.get("initial_spatial_basis"),
        transport_mode=str(source.get("transport_mode", "fermi_walker")),
        global_length_in_local_units=ephemeris_frame._document_float(
            source, "global_length_in_local_units", 1.0
        ),
        acceleration_scale_global_units=_default_acceleration_scale(source),
        metric_fd_step=metric_fd_step,
        integration_substeps_per_interval=ephemeris_frame._document_int(
            source, "integration_substeps_per_interval", 16
        ),
        geodesic_acceleration_tolerance=ephemeris_frame._document_float(
            source, "geodesic_acceleration_tolerance_scaled", 1.0e-8
        ),
        audit_subsamples_per_interval=ephemeris_frame._document_int(
            source, "audit_subsamples_per_interval", 4
        ),
        run_coarse_convergence_check=ephemeris_frame._document_bool(
            source, "run_coarse_convergence_check", True
        ),
        run_metric_fd_convergence_check=ephemeris_frame._document_bool(
            source, "run_metric_fd_convergence_check", True
        ),
    )
    generator = document["generator"]
    generator["source_manifest"] = str(scan.path)
    generator["source_manifest_sha256"] = static.file_sha256(scan.path)
    generator["source_ephemeris"] = str(ephemeris_path)
    generator["source_ephemeris_sha256"] = static.file_sha256(ephemeris_path)
    generator["source_level"] = scan.descriptors[0].source_level
    generator["source_snapshot_storage"] = scan.descriptors[0].source_storage
    hash_source_files = (
        scan.hash_source_files
        if hash_source_files_override is None
        else hash_source_files_override
    )
    if not isinstance(hash_source_files, bool):
        raise ValueError("hash_source_files_override must be true or false")
    source_snapshots = []
    for descriptor in scan.descriptors:
        provenance: dict[str, object] = {
            "time": descriptor.time,
            "cycle": descriptor.cycle,
            "state": str(descriptor.state_path),
            "adm": str(descriptor.adm_path),
            "source_level": descriptor.source_level,
            "source_storage": descriptor.source_storage,
            "selected_leaf_meshblocks": descriptor.source_meshblock_count,
        }
        if hash_source_files:
            provenance["state_sha256"] = static.file_sha256(descriptor.state_path)
            provenance["adm_sha256"] = static.file_sha256(descriptor.adm_path)
        source_snapshots.append(provenance)
    generator["source_snapshots"] = source_snapshots
    generator["source_file_hashes_recorded"] = hash_source_files
    document["diagnostics"]["snapshot_loading"] = snapshots.loading_document()
    return document


def main() -> int:
    arguments = parse_arguments()
    scan = extract.scan_snapshot_manifest(arguments.manifest)
    document = build_frame_from_scan(scan, arguments.ephemeris)
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
