#!/usr/bin/env python3
"""Stitch a normalized BBH inspiral table to one canonical Kerr remnant.

The input and output use AthenaK's 21-column trajectory schema.  During the
transition the two source positions are blended to one kicked remnant worldline.
The velocities are the analytic derivatives of those blended positions.  After
the transition, term 1 contains the full remnant and term 2 is exactly zero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Sequence

import numpy as np


COLUMN_NAMES = (
    "t",
    "x1",
    "y1",
    "z1",
    "x2",
    "y2",
    "z2",
    "vx1",
    "vy1",
    "vz1",
    "vx2",
    "vy2",
    "vz2",
    "a1x",
    "a1y",
    "a1z",
    "a2x",
    "a2y",
    "a2z",
    "m1_full",
    "m2_full",
)
SCHEMA_VERSION = 1
PROVENANCE_SCHEMA_VERSION = 2
TOOL_VERSION = 2
WEIGHT_FORMULA_VERSION = "cinf-bump-ratio-v1"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line options."""
    parser = argparse.ArgumentParser(
        description=(
            "Stitch a normalized 21-column inspiral trajectory to a canonical "
            "single-term Kerr remnant."
        )
    )
    parser.add_argument("input", type=Path, help="normalized 21-column input table")
    parser.add_argument("output", type=Path, help="new stitched output table")
    parser.add_argument(
        "--source-model-label",
        action="append",
        required=True,
        help=(
            "auditable source-model label supplied by the caller; repeat once "
            "for every applicable label"
        ),
    )
    parser.add_argument(
        "--source-revision",
        required=True,
        help="non-empty source code/data revision supplied by the caller",
    )
    parser.add_argument(
        "--source-provenance",
        type=Path,
        default=None,
        help="optional source manifest or provenance artifact to hash and record",
    )
    start = parser.add_mutually_exclusive_group(required=True)
    start.add_argument(
        "--transition-start-time",
        type=float,
        help="transition start in input-table time units",
    )
    start.add_argument(
        "--transition-start-separation",
        type=float,
        help=(
            "select the first input row whose coordinate separation is no larger "
            "than this value"
        ),
    )
    parser.add_argument(
        "--transition-end-time",
        type=float,
        required=True,
        help="transition end in input-table time units",
    )
    parser.add_argument(
        "--final-mass",
        type=float,
        required=True,
        help="final remnant mass Mf in normalized code units",
    )
    parser.add_argument(
        "--final-a",
        type=float,
        nargs=3,
        metavar=("AX", "AY", "AZ"),
        required=True,
        help="length-valued final Kerr a-vector (not dimensionless chi)",
    )
    parser.add_argument(
        "--kick",
        type=float,
        nargs=3,
        metavar=("VX", "VY", "VZ"),
        required=True,
        help="constant post-merger kick velocity in units with c=1",
    )
    parser.add_argument(
        "--postmerger-duration",
        type=float,
        required=True,
        help="duration to append after the transition endpoint",
    )
    parser.add_argument(
        "--postmerger-dt",
        type=float,
        default=None,
        help=(
            "post-merger table cadence; default is the median of the last up to "
            "32 input intervals ending before the transition endpoint"
        ),
    )
    parser.add_argument(
        "--max-relative-dxdt-v-error",
        type=float,
        default=None,
        help=(
            "reject the finished table if the maximum full-table relative "
            "finite-difference |dx/dt-v| error exceeds this value; by default "
            "record the diagnostic without enforcing a limit"
        ),
    )
    parser.add_argument(
        "--provenance",
        type=Path,
        default=None,
        help="JSON provenance path (default: OUTPUT.provenance.json)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace an existing output table and provenance file",
    )
    return parser.parse_args(argv)


def file_sha256(path: Path) -> str:
    """Return a lowercase SHA-256 digest for *path*."""
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_metadata(
    args: argparse.Namespace, output_path: Path, provenance_path: Path
) -> dict:
    """Validate caller declarations without inferring a physical source model."""
    labels = [label.strip() for label in args.source_model_label]
    if not labels or any(not label for label in labels):
        raise ValueError("every --source-model-label must be non-empty")
    if len(set(labels)) != len(labels):
        raise ValueError("--source-model-label values must be unique")
    source_revision = args.source_revision.strip()
    if not source_revision:
        raise ValueError("--source-revision must be non-empty")

    source_provenance = None
    if args.source_provenance is not None:
        source_path = args.source_provenance.resolve()
        if source_path in (output_path, provenance_path):
            raise ValueError(
                "--source-provenance must differ from output artifact paths"
            )
        if not source_path.is_file():
            raise ValueError(
                f"--source-provenance is not a readable file: {source_path}"
            )
        source_provenance = {
            "path": str(source_path),
            "sha256": file_sha256(source_path),
        }
    return {
        "declaration_source": "CLI",
        "labels": labels,
        "source_provenance": source_provenance,
        "source_revision": source_revision,
    }


def load_table(path: Path) -> np.ndarray:
    """Load and strictly validate the normalized 21-column schema."""
    try:
        table = np.loadtxt(path, dtype=np.float64, ndmin=2)
    except (OSError, ValueError) as exc:
        message = f"could not read a numeric trajectory from {path}: {exc}"
        raise ValueError(message) from exc
    if table.ndim != 2 or table.shape[1] != len(COLUMN_NAMES):
        raise ValueError(
            "input trajectory must have exactly 21 numeric columns in the order "
            + " ".join(COLUMN_NAMES)
        )
    if table.shape[0] < 2:
        raise ValueError("input trajectory must contain at least two rows")
    if not np.all(np.isfinite(table)):
        raise ValueError("input trajectory contains a non-finite value")
    if not np.all(np.diff(table[:, 0]) > 0.0):
        raise ValueError("input trajectory times must be strictly increasing")
    _validate_physical_state(table, "input trajectory")
    return table


def _state_views(table: np.ndarray) -> tuple[np.ndarray, ...]:
    """Return position, velocity, Kerr-a, and mass views with hole axes."""
    positions = np.stack((table[:, 1:4], table[:, 4:7]), axis=1)
    velocities = np.stack((table[:, 7:10], table[:, 10:13]), axis=1)
    spins = np.stack((table[:, 13:16], table[:, 16:19]), axis=1)
    masses = table[:, 19:21]
    return positions, velocities, spins, masses


def _validate_physical_state(table: np.ndarray, label: str) -> None:
    """Reject states that cannot represent subluminal, subextremal Kerr terms."""
    positions, velocities, spins, masses = _state_views(table)
    del positions
    if not np.all(np.isfinite(table)):
        raise ValueError(f"{label} contains a non-finite value")
    speed = np.linalg.norm(velocities, axis=2)
    if np.any(speed >= 1.0):
        index = np.argwhere(speed >= 1.0)[0]
        raise ValueError(
            f"{label} contains a superluminal velocity at row {index[0]}, "
            f"term {index[1] + 1}"
        )
    if np.any(masses < 0.0):
        raise ValueError(f"{label} contains a negative individual mass")
    if np.any(np.sum(masses, axis=1) <= 0.0):
        raise ValueError(f"{label} must have positive total mass in every row")
    spin_magnitude = np.linalg.norm(spins, axis=2)
    tolerance = 64.0 * np.finfo(np.float64).eps * np.maximum(1.0, masses)
    if np.any(spin_magnitude > masses + tolerance):
        index = np.argwhere(spin_magnitude > masses + tolerance)[0]
        raise ValueError(
            f"{label} contains a superextremal Kerr term at row {index[0]}, "
            f"term {index[1] + 1}: |a|={spin_magnitude[tuple(index)]:.17g} > "
            f"m={masses[tuple(index)]:.17g}"
        )


def cinf_step(
    query_time: np.ndarray, start_time: float, end_time: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return a C-infinity step and its analytic time derivative.

    For 0 < u < 1, W = f(u)/(f(u)+f(1-u)) with f(u)=exp(-1/u).
    Values indistinguishable from an endpoint in binary64 are deliberately
    rounded to exactly zero or one together with a zero derivative.
    """
    query_time = np.asarray(query_time, dtype=np.float64)
    weight = np.zeros_like(query_time)
    derivative = np.zeros_like(query_time)
    weight[query_time >= end_time] = 1.0
    interior = (query_time > start_time) & (query_time < end_time)
    u = (query_time[interior] - start_time) / (end_time - start_time)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        logistic_argument = 1.0 / u - 1.0 / (1.0 - u)
    local_weight = np.empty_like(u)
    positive = logistic_argument >= 0.0
    exp_negative = np.exp(-logistic_argument[positive])
    local_weight[positive] = exp_negative / (1.0 + exp_negative)
    exp_positive = np.exp(logistic_argument[~positive])
    local_weight[~positive] = 1.0 / (1.0 + exp_positive)

    local_derivative = np.zeros_like(u)
    resolved = np.abs(logistic_argument) < 350.0
    ur = u[resolved]
    slope_u = 1.0 / (ur * ur) + 1.0 / ((1.0 - ur) * (1.0 - ur))
    wr = local_weight[resolved]
    local_derivative[resolved] = (
        wr * (1.0 - wr) * slope_u / (end_time - start_time)
    )
    weight[interior] = local_weight
    derivative[interior] = local_derivative
    return weight, derivative


def _interval_coordinates(
    source_time: np.ndarray, query_time: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Locate query points and return left indices, normalized position, and dt."""
    index = np.searchsorted(source_time, query_time, side="right") - 1
    index = np.clip(index, 0, source_time.size - 2)
    interval = source_time[index + 1] - source_time[index]
    fraction = (query_time - source_time[index]) / interval
    return index, fraction, interval


def _sample_inspiral(
    table: np.ndarray, query_time: np.ndarray
) -> tuple[np.ndarray, ...]:
    """Sample source data, using position/velocity cubic Hermite pairs."""
    source_time = table[:, 0]
    positions, velocities, spins, masses = _state_views(table)
    index, fraction, interval = _interval_coordinates(source_time, query_time)
    s = fraction[:, None, None]
    h = interval[:, None, None]
    p0 = positions[index]
    p1 = positions[index + 1]
    v0 = velocities[index]
    v1 = velocities[index + 1]

    h00 = 2.0 * s**3 - 3.0 * s**2 + 1.0
    h10 = s**3 - 2.0 * s**2 + s
    h01 = -2.0 * s**3 + 3.0 * s**2
    h11 = s**3 - s**2
    sampled_positions = h00 * p0 + h10 * h * v0 + h01 * p1 + h11 * h * v1

    dh00 = (6.0 * s**2 - 6.0 * s) / h
    dh10 = 3.0 * s**2 - 4.0 * s + 1.0
    dh01 = (-6.0 * s**2 + 6.0 * s) / h
    dh11 = 3.0 * s**2 - 2.0 * s
    sampled_velocities = dh00 * p0 + dh10 * v0 + dh01 * p1 + dh11 * v1

    linear_fraction = fraction[:, None, None]
    sampled_spins = (
        (1.0 - linear_fraction) * spins[index]
        + linear_fraction * spins[index + 1]
    )
    mass_fraction = fraction[:, None]
    sampled_masses = (
        (1.0 - mass_fraction) * masses[index]
        + mass_fraction * masses[index + 1]
    )
    return sampled_positions, sampled_velocities, sampled_spins, sampled_masses


def _select_start_time(table: np.ndarray, args: argparse.Namespace) -> tuple[float, dict]:
    """Select a start time explicitly or at the first separation threshold row."""
    if args.transition_start_time is not None:
        start_time = float(args.transition_start_time)
        selection = {"mode": "time", "requested_time": start_time}
    else:
        threshold = float(args.transition_start_separation)
        if not np.isfinite(threshold) or threshold <= 0.0:
            raise ValueError("--transition-start-separation must be finite and positive")
        positions, _, _, _ = _state_views(table)
        separation = np.linalg.norm(positions[:, 0] - positions[:, 1], axis=1)
        matches = np.flatnonzero(separation <= threshold)
        if matches.size == 0:
            raise ValueError(
                "input separation never reaches --transition-start-separation="
                f"{threshold:.17g}"
            )
        row = int(matches[0])
        start_time = float(table[row, 0])
        selection = {
            "mode": "first-row-at-or-below-separation",
            "requested_separation": threshold,
            "selected_row": row,
            "selected_separation": float(separation[row]),
        }
    return start_time, selection


def _validate_parameters(
    table: np.ndarray,
    start_time: float,
    end_time: float,
    final_mass: float,
    final_a: np.ndarray,
    kick: np.ndarray,
    postmerger_duration: float,
    postmerger_dt: float | None,
    max_relative_dxdt_v_error: float | None,
) -> None:
    """Validate all construction parameters before creating an artifact."""
    scalars = (start_time, end_time, final_mass, postmerger_duration)
    if not all(np.isfinite(value) for value in scalars):
        raise ValueError("all scalar construction parameters must be finite")
    if not np.all(np.isfinite(final_a)) or not np.all(np.isfinite(kick)):
        raise ValueError("--final-a and --kick must contain only finite values")
    if start_time < table[0, 0] or start_time >= table[-1, 0]:
        raise ValueError("transition start must lie in the input time range")
    if end_time <= start_time or end_time > table[-1, 0]:
        raise ValueError(
            "transition end must be later than the start and no later than the "
            "last input-table time"
        )
    if final_mass <= 0.0:
        raise ValueError("--final-mass must be finite and positive")
    if np.linalg.norm(final_a) > final_mass:
        raise ValueError("--final-a must satisfy |a_f| <= Mf")
    if np.linalg.norm(kick) >= 1.0:
        raise ValueError("--kick must have magnitude strictly below the speed of light")
    if postmerger_duration < 0.0:
        raise ValueError("--postmerger-duration must be non-negative")
    if postmerger_dt is not None and (
        not np.isfinite(postmerger_dt) or postmerger_dt <= 0.0
    ):
        raise ValueError("--postmerger-dt must be finite and positive")
    if max_relative_dxdt_v_error is not None and (
        not np.isfinite(max_relative_dxdt_v_error)
        or max_relative_dxdt_v_error <= 0.0
    ):
        raise ValueError(
            "--max-relative-dxdt-v-error must be finite and positive"
        )


def _default_postmerger_dt(table: np.ndarray, end_time: float) -> float:
    """Choose a deterministic local input cadence."""
    retained_time = table[table[:, 0] <= end_time, 0]
    if retained_time.size < 2:
        retained_time = table[:2, 0]
    intervals = np.diff(retained_time)[-32:]
    cadence = float(np.median(intervals))
    if not np.isfinite(cadence) or cadence <= 0.0:
        raise ValueError("could not infer a finite positive post-merger cadence")
    return cadence


def _postmerger_times(
    end_time: float, duration: float, cadence: float
) -> np.ndarray:
    """Return deterministic cadence samples plus the exact requested endpoint."""
    if duration == 0.0:
        return np.empty(0, dtype=np.float64)
    final_time = end_time + duration
    count = int(np.floor(duration / cadence))
    samples = end_time + cadence * np.arange(1, count + 1, dtype=np.float64)
    scale = max(1.0, abs(final_time))
    samples = samples[samples < final_time - 16.0 * np.finfo(float).eps * scale]
    return np.concatenate((samples, np.asarray([final_time], dtype=np.float64)))


def stitch_table(
    table: np.ndarray,
    start_time: float,
    end_time: float,
    final_mass: float,
    final_a: np.ndarray,
    kick: np.ndarray,
    postmerger_duration: float,
    postmerger_dt: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Construct and validate a canonical single-remnant output table."""
    before = table[table[:, 0] < start_time].copy()
    source_transition_times = table[
        (table[:, 0] >= start_time) & (table[:, 0] <= end_time), 0
    ]
    transition_time = np.unique(
        np.concatenate(
            (
                np.asarray([start_time]),
                source_transition_times,
                np.asarray([end_time]),
            )
        )
    )
    positions, velocities, spins, masses = _sample_inspiral(
        table, transition_time
    )
    start_positions, _, _, start_masses = _sample_inspiral(
        table, np.asarray([start_time])
    )
    total_start_mass = float(np.sum(start_masses[0]))
    com_start = np.sum(
        start_masses[0, :, None] * start_positions[0], axis=0
    ) / total_start_mass

    weight, weight_derivative = cinf_step(
        transition_time, start_time, end_time
    )
    remnant_position = com_start + (
        transition_time[:, None] - start_time
    ) * kick
    w_vector = weight[:, None, None]
    wd_vector = weight_derivative[:, None, None]
    remnant_for_terms = remnant_position[:, None, :]
    kick_for_terms = kick[None, None, :]
    blended_positions = (
        (1.0 - w_vector) * positions + w_vector * remnant_for_terms
    )
    blended_velocities = (
        (1.0 - w_vector) * velocities
        + w_vector * kick_for_terms
        + wd_vector * (remnant_for_terms - positions)
    )

    target_spins = np.zeros_like(spins)
    target_spins[:, 0, :] = final_a
    blended_spins = (1.0 - w_vector) * spins + w_vector * target_spins
    target_masses = np.zeros_like(masses)
    target_masses[:, 0] = final_mass
    w_scalar = weight[:, None]
    blended_masses = (1.0 - w_scalar) * masses + w_scalar * target_masses

    transition = np.column_stack(
        (
            transition_time,
            blended_positions[:, 0, :],
            blended_positions[:, 1, :],
            blended_velocities[:, 0, :],
            blended_velocities[:, 1, :],
            blended_spins[:, 0, :],
            blended_spins[:, 1, :],
            blended_masses,
        )
    )

    appended_time = _postmerger_times(
        end_time, postmerger_duration, postmerger_dt
    )
    postmerger = np.empty((appended_time.size, len(COLUMN_NAMES)))
    if appended_time.size:
        postmerger[:, 0] = appended_time
        path = com_start + (appended_time[:, None] - start_time) * kick
        postmerger[:, 1:4] = path
        postmerger[:, 4:7] = path
        postmerger[:, 7:10] = kick
        postmerger[:, 10:13] = kick
        postmerger[:, 13:16] = final_a
        postmerger[:, 16:19] = 0.0
        postmerger[:, 19] = final_mass
        postmerger[:, 20] = 0.0

    output = np.vstack((before, transition, postmerger))
    if output.shape[1] != len(COLUMN_NAMES):
        raise RuntimeError("internal error: stitched table has the wrong schema")
    if not np.all(np.diff(output[:, 0]) > 0.0):
        raise RuntimeError("internal error: stitched output times are not increasing")
    _validate_physical_state(output, "stitched trajectory")
    return output, com_start


def _time_sampling_summary(time: np.ndarray) -> dict:
    """Return reproducible interval diagnostics for one table-time subset."""
    intervals = np.diff(np.asarray(time, dtype=np.float64))
    if intervals.size == 0:
        return {
            "interval_count": 0,
            "interval_sha256_float64_le": None,
            "max_dt": None,
            "median_dt": None,
            "min_dt": None,
            "row_count": int(np.size(time)),
            "uniform_within_1e-12": True,
        }
    little_endian = np.asarray(intervals, dtype="<f8")
    return {
        "interval_count": int(intervals.size),
        "interval_sha256_float64_le": hashlib.sha256(
            little_endian.tobytes()
        ).hexdigest(),
        "max_dt": float(np.max(intervals)),
        "median_dt": float(np.median(intervals)),
        "min_dt": float(np.min(intervals)),
        "row_count": int(np.size(time)),
        "uniform_within_1e-12": bool(
            np.allclose(intervals, intervals[0], rtol=1.0e-12, atol=0.0)
        ),
    }


def velocity_consistency_diagnostics(table: np.ndarray) -> dict:
    """Compare stored velocities with full-table finite differences of position."""
    time = table[:, 0]
    positions, velocities, _, _ = _state_views(table)
    edge_order = 2 if time.size >= 3 else 1
    numerical_velocity = np.gradient(
        positions, time, axis=0, edge_order=edge_order
    )
    absolute_error = np.linalg.norm(numerical_velocity - velocities, axis=2)
    normalization = np.maximum(
        np.maximum(
            np.linalg.norm(numerical_velocity, axis=2),
            np.linalg.norm(velocities, axis=2),
        ),
        1.0e-12,
    )
    relative_error = absolute_error / normalization
    if not np.all(np.isfinite(absolute_error)) or not np.all(
        np.isfinite(relative_error)
    ):
        raise ValueError("dx/dt-v consistency diagnostic produced a non-finite value")

    def statistics(values: np.ndarray, term_offset: int = 0) -> dict:
        row, term = np.unravel_index(int(np.argmax(values)), values.shape)
        return {
            "max": float(values[row, term]),
            "max_row": int(row),
            "max_term": int(term + term_offset + 1),
            "max_time": float(time[row]),
            "rms": float(np.sqrt(np.mean(values * values))),
        }

    return {
        "absolute_l2": statistics(absolute_error),
        "edge_order": edge_order,
        "normalization": "max(|v_finite_difference|, |v_table|, 1e-12)",
        "relative_l2": statistics(relative_error),
        "row_count": int(time.size),
        "scheme": "numpy.gradient over the full output table",
        "term_1": {
            "absolute_l2": statistics(absolute_error[:, 0:1]),
            "relative_l2": statistics(relative_error[:, 0:1]),
        },
        "term_2": {
            "absolute_l2": statistics(absolute_error[:, 1:2], term_offset=1),
            "relative_l2": statistics(relative_error[:, 1:2], term_offset=1),
        },
    }


def _stage_table(target: Path, table: np.ndarray) -> Path:
    """Write and fsync a deterministic table in the target directory."""
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        header = (
            "canonical_single_term_remnant = true; schema_version = 1; "
            f"weight_formula = {WEIGHT_FORMULA_VERSION}\n"
            + " ".join(COLUMN_NAMES)
        )
        np.savetxt(temporary, table, fmt="%.17e", header=header)
        with temporary.open("rb") as artifact:
            os.fsync(artifact.fileno())
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return temporary


def _stage_json(target: Path, document: dict) -> Path:
    """Write and fsync deterministic JSON in the target directory."""
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            json.dump(document, stream, indent=2, sort_keys=True, ensure_ascii=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return temporary


def _backup_link(target: Path) -> Path:
    """Create a same-filesystem hard-link backup without modifying *target*."""
    descriptor, backup_name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".rollback", dir=target.parent
    )
    os.close(descriptor)
    backup = Path(backup_name)
    backup.unlink()
    try:
        os.link(target, backup)
    except BaseException:
        backup.unlink(missing_ok=True)
        raise
    return backup


def _publish_transaction(
    staged_targets: Sequence[tuple[Path, Path]], force: bool
) -> None:
    """Publish a table/provenance set and roll back every partial installation.

    Each individual installation is atomic.  With ``force``, hard-link backups
    preserve the exact old inodes until every replacement succeeds.  Without
    ``force``, hard links provide atomic no-clobber creation.  Any failure after
    the first installation restores all pre-call targets before propagating the
    original exception.
    """
    targets = [target for _, target in staged_targets]
    if len(set(targets)) != len(targets):
        raise ValueError("transaction targets must be distinct")
    if not force:
        for target in targets:
            if target.exists():
                raise FileExistsError(
                    f"refusing to replace existing output {target}"
                )

    backups: dict[Path, Path | None] = {target: None for target in targets}
    if force:
        try:
            for target in targets:
                if target.exists():
                    backups[target] = _backup_link(target)
        except BaseException:
            for backup in backups.values():
                if backup is not None:
                    backup.unlink(missing_ok=True)
            raise

    installed: list[Path] = []
    try:
        for staged, target in staged_targets:
            if force:
                os.replace(staged, target)
                installed.append(target)
            else:
                try:
                    os.link(staged, target)
                except FileExistsError:
                    raise FileExistsError(
                        f"refusing to replace existing output {target}"
                    )
                installed.append(target)
                staged.unlink()
    except BaseException as publish_error:
        rollback_errors: list[str] = []
        for target in reversed(installed):
            backup = backups[target]
            try:
                if backup is None:
                    target.unlink(missing_ok=True)
                else:
                    os.replace(backup, target)
                    backups[target] = None
            except OSError as rollback_error:
                rollback_errors.append(f"{target}: {rollback_error}")

        if rollback_errors:
            surviving_backups = [
                str(path) for path in backups.values() if path is not None
            ]
            detail = "; ".join(rollback_errors)
            recovery = ", ".join(surviving_backups) or "none"
            raise RuntimeError(
                "output transaction failed and rollback was incomplete; "
                f"rollback errors: {detail}; retained backups: {recovery}"
            ) from publish_error

        for backup in backups.values():
            if backup is not None:
                backup.unlink(missing_ok=True)
        raise
    else:
        for backup in backups.values():
            if backup is not None:
                backup.unlink(missing_ok=True)


def run(args: argparse.Namespace) -> tuple[Path, Path]:
    """Validate inputs, construct artifacts, and atomically publish them."""
    input_path = args.input.resolve()
    output_path = args.output.resolve()
    provenance_path = (
        args.provenance.resolve()
        if args.provenance is not None
        else Path(str(output_path) + ".provenance.json")
    )
    if output_path == input_path:
        raise ValueError("input and output paths must be different")
    if provenance_path in (input_path, output_path):
        raise ValueError("provenance path must differ from input and output paths")
    source_metadata = _source_metadata(args, output_path, provenance_path)
    if not args.force:
        for target in (output_path, provenance_path):
            if target.exists():
                raise FileExistsError(f"refusing to replace existing output {target}")

    table = load_table(input_path)
    start_time, start_selection = _select_start_time(table, args)
    end_time = float(args.transition_end_time)
    final_mass = float(args.final_mass)
    final_a = np.asarray(args.final_a, dtype=np.float64)
    kick = np.asarray(args.kick, dtype=np.float64)
    duration = float(args.postmerger_duration)
    requested_dt = args.postmerger_dt
    max_relative_error = args.max_relative_dxdt_v_error
    _validate_parameters(
        table,
        start_time,
        end_time,
        final_mass,
        final_a,
        kick,
        duration,
        requested_dt,
        max_relative_error,
    )
    cadence = (
        _default_postmerger_dt(table, end_time)
        if requested_dt is None
        else float(requested_dt)
    )
    stitched, com_start = stitch_table(
        table,
        start_time,
        end_time,
        final_mass,
        final_a,
        kick,
        duration,
        cadence,
    )
    velocity_diagnostics = velocity_consistency_diagnostics(stitched)
    measured_relative_error = velocity_diagnostics["relative_l2"]["max"]
    velocity_diagnostics["threshold"] = {
        "enforced": max_relative_error is not None,
        "limit": max_relative_error,
        "metric": "relative_l2.max",
        "passed": (
            None
            if max_relative_error is None
            else measured_relative_error <= max_relative_error
        ),
    }
    if (
        max_relative_error is not None
        and measured_relative_error > max_relative_error
    ):
        raise ValueError(
            "full-table relative |dx/dt-v| error "
            f"{measured_relative_error:.6e} exceeds "
            f"--max-relative-dxdt-v-error={max_relative_error:.6e}"
        )

    transition_time = stitched[
        (stitched[:, 0] >= start_time) & (stitched[:, 0] <= end_time), 0
    ]
    transition_sampling = _time_sampling_summary(transition_time)
    final_chi = final_a / final_mass

    output_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_table: Path | None = None
    temporary_provenance: Path | None = None
    try:
        temporary_table = _stage_table(output_path, stitched)
        output_digest = file_sha256(temporary_table)
        provenance = {
            "provenance_schema_version": PROVENANCE_SCHEMA_VERSION,
            "input": {
                "path": str(input_path),
                "rows": int(table.shape[0]),
                "sha256": file_sha256(input_path),
                "time_bounds": [float(table[0, 0]), float(table[-1, 0])],
            },
            "output": {
                "path": str(output_path),
                "rows": int(stitched.shape[0]),
                "sha256": output_digest,
                "time_bounds": [float(stitched[0, 0]), float(stitched[-1, 0])],
            },
            "remnant": {
                "a": final_a.tolist(),
                "a_length": final_a.tolist(),
                "a_magnitude": float(np.linalg.norm(final_a)),
                "a_to_chi_relation": "a_length = Mf * chi; chi = a_length / Mf",
                "chi": final_chi.tolist(),
                "chi_magnitude": float(np.linalg.norm(final_chi)),
                "kick": kick.tolist(),
                "mass": final_mass,
                "representation": "canonical-single-term-1",
            },
            "trajectory_model": source_metadata,
            "sampling": {
                "postmerger_dt": cadence,
                "postmerger_dt_source": (
                    "local-input-median" if requested_dt is None else "explicit"
                ),
                "postmerger_duration": duration,
            },
            "schema": {
                "columns": list(COLUMN_NAMES),
                "version": SCHEMA_VERSION,
            },
            "tool": {
                "name": Path(__file__).name,
                "version": TOOL_VERSION,
            },
            "transition": {
                "com_position_at_start": com_start.tolist(),
                "end_time": end_time,
                "sampling": transition_sampling,
                "start_selection": start_selection,
                "start_time": start_time,
                "weight": {
                    "endpoint_definition": "W=0 for u<=0; W=1 for u>=1",
                    "formula": (
                        "u=(t-t0)/(t1-t0); f(s)=exp(-1/s); "
                        "W=f(u)/(f(u)+f(1-u)) for 0<u<1"
                    ),
                    "formula_version": WEIGHT_FORMULA_VERSION,
                    "regularity": "C-infinity",
                },
            },
            "velocity_consistency": velocity_diagnostics,
        }
        temporary_provenance = _stage_json(provenance_path, provenance)
        _publish_transaction(
            (
                (temporary_table, output_path),
                (temporary_provenance, provenance_path),
            ),
            args.force,
        )
        temporary_table = None
        temporary_provenance = None
    finally:
        if temporary_table is not None:
            temporary_table.unlink(missing_ok=True)
        if temporary_provenance is not None:
            temporary_provenance.unlink(missing_ok=True)
    return output_path, provenance_path


def main(argv: Sequence[str] | None = None) -> None:
    """CLI entry point."""
    output_path, provenance_path = run(parse_args(argv))
    print(f"wrote trajectory: {output_path}")
    print(f"wrote provenance: {provenance_path}")


if __name__ == "__main__":
    main()
