#!/usr/bin/env python3
"""Convert AnalyticalBBH trajectory HDF5 data to AthenaK's portable table format."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import warnings

import h5py
import numpy as np


STATE_NAMES = (
    "x1", "y1", "z1", "x2", "y2", "z2",
    "vx1", "vy1", "vz1", "vx2", "vy2", "vz2",
    "a1x", "a1y", "a1z", "a2x", "a2y", "a2z",
    "m1_full", "m2_full",
)

MANIFEST_SCHEMA_VERSION = 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert an AnalyticalBBH HDF5 trajectory to 21-column ASCII."
    )
    parser.add_argument("input", type=Path, help="input HDF5 trajectory")
    parser.add_argument("output", type=Path, help="output ASCII trajectory")
    parser.add_argument(
        "--mass-unit",
        type=float,
        default=None,
        help=(
            "input mass used as one AthenaK length/time unit; by default use the "
            "initial total mass"
        ),
    )
    parser.add_argument(
        "--assume-missing-spin-zero",
        action="store_true",
        help="treat missing legacy s1*/s2* spin components as zero",
    )
    parser.add_argument(
        "--velocity-tolerance",
        type=float,
        default=5.0e-2,
        help=(
            "warn when finite-difference dx/dt disagrees fractionally with v "
            "(default: 0.05)"
        ),
    )
    parser.add_argument(
        "--strict-velocity-mismatch",
        "--strict-velocity",
        dest="strict_velocity_mismatch",
        action="store_true",
        help=(
            "reject the trajectory instead of warning when |dx/dt-v| exceeds "
            "--velocity-tolerance"
        ),
    )
    parser.add_argument(
        "--validation-profile",
        choices=("none", "q1-nonspinning"),
        default="none",
        help=(
            "optional physical baseline checks; q1-nonspinning requires an "
            "equal-mass, centered, initially nonspinning binary and an "
            "unambiguous final remnant representation"
        ),
    )
    parser.add_argument(
        "--validation-tolerance",
        type=float,
        default=1.0e-8,
        help=(
            "absolute/relative tolerance in normalized code units for the "
            "selected validation profile (default: 1e-8)"
        ),
    )
    manifest_group = parser.add_mutually_exclusive_group()
    manifest_group.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help=(
            "JSON manifest path (default: OUTPUT.manifest.json); the manifest "
            "is deterministic and records endpoint diagnostics and hashes"
        ),
    )
    manifest_group.add_argument(
        "--no-manifest",
        action="store_true",
        help="do not create the default JSON manifest",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace an existing output file",
    )
    return parser.parse_args()


def _sha256(path: Path) -> str:
    """Return the lowercase SHA-256 digest of one file."""
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _temporary_path(target: Path) -> Path:
    """Allocate a private staging path on the target filesystem."""
    descriptor, name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
    )
    os.close(descriptor)
    return Path(name)


def _backup_link(target: Path) -> Path:
    """Create a same-filesystem hard-link backup without moving the target."""
    descriptor, name = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".rollback", dir=target.parent
    )
    os.close(descriptor)
    backup = Path(name)
    backup.unlink()
    try:
        os.link(target, backup)
    except BaseException:
        backup.unlink(missing_ok=True)
        raise
    return backup


def _publish_staged(staged_targets: list[tuple[Path, Path]], force: bool) -> None:
    """Publish staged artifacts with atomic installation and pair rollback.

    With ``force``, hard-link backups preserve old targets until every replacement
    succeeds.  Without ``force``, linking the staged inode into its final name is an
    atomic no-clobber operation, including when another process creates a target after
    the preliminary existence check.  A failed later installation rolls back every
    target installed by this call without deleting a concurrently created target.
    """
    try:
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
                    except FileExistsError as error:
                        raise FileExistsError(
                            f"refusing to replace existing output {target}"
                        ) from error
                    installed.append(target)
                    staged.unlink()

            for directory in {target.parent for target in targets}:
                descriptor = os.open(directory, os.O_RDONLY)
                try:
                    os.fsync(descriptor)
                finally:
                    os.close(descriptor)
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
    finally:
        for staged, _ in staged_targets:
            staged.unlink(missing_ok=True)


def _endpoint_summary(
    row: int, time: np.ndarray, state: dict[str, np.ndarray]
) -> dict[str, object]:
    """Return JSON-safe physical diagnostics for one endpoint."""
    positions = np.asarray(
        [[state[f"{component}{hole}"][row] for component in "xyz"]
         for hole in (1, 2)],
        dtype=np.float64,
    )
    velocities = np.asarray(
        [[state[f"v{component}{hole}"][row] for component in "xyz"]
         for hole in (1, 2)],
        dtype=np.float64,
    )
    spin_a = np.asarray(
        [[state[f"a{hole}{component}"][row] for component in "xyz"]
         for hole in (1, 2)],
        dtype=np.float64,
    )
    masses = np.asarray(
        [state["m1_full"][row], state["m2_full"][row]], dtype=np.float64
    )
    return {
        "time": float(time[row]),
        "masses": masses.tolist(),
        "total_mass": float(np.sum(masses)),
        "positions": positions.tolist(),
        "separation": float(np.linalg.norm(positions[0] - positions[1])),
        "velocities": velocities.tolist(),
        "speeds": np.linalg.norm(velocities, axis=1).tolist(),
        "spin_a": spin_a.tolist(),
    }


def _hermite_velocity_certificate(
    time: np.ndarray, state: dict[str, np.ndarray]
) -> dict[str, object]:
    """Mirror dynbbh's rigorous subluminal cubic-Hermite certificate.

    On an interval, the Hermite velocity is a quadratic Bezier curve whose
    controls are ``v0``, ``3*(p1-p0)/dt-v0-v1``, and ``v1``.  The open unit ball
    is convex, so checking every endpoint velocity and middle control certifies
    the complete interpolated worldline rather than only sampled table rows.
    """
    endpoint_records: list[dict[str, object]] = []
    middle_records: list[dict[str, object]] = []
    hole_records: dict[str, dict[str, object]] = {}
    interval_duration = np.diff(time)

    for hole in (1, 2):
        position = np.column_stack(
            tuple(state[f"{component}{hole}"] for component in "xyz")
        )
        velocity = np.column_stack(
            tuple(state[f"v{component}{hole}"] for component in "xyz")
        )
        endpoint_speed_squared = np.einsum("ij,ij->i", velocity, velocity)
        endpoint_row = int(np.argmax(endpoint_speed_squared))
        endpoint_record = {
            "hole": hole,
            "row_index": endpoint_row,
            "time": float(time[endpoint_row]),
            "velocity": velocity[endpoint_row].tolist(),
            "speed_squared": float(endpoint_speed_squared[endpoint_row]),
            "speed": float(np.sqrt(endpoint_speed_squared[endpoint_row])),
        }
        endpoint_records.append(endpoint_record)
        invalid_endpoints = np.flatnonzero(
            ~np.isfinite(endpoint_speed_squared) | (endpoint_speed_squared >= 1.0)
        )
        if invalid_endpoints.size:
            row = int(invalid_endpoints[0])
            raise ValueError(
                "trajectory contains a superluminal black-hole velocity at "
                f"row {row}, hole {hole}: |v|^2="
                f"{endpoint_speed_squared[row]:.17g}"
            )

        with np.errstate(over="ignore", invalid="ignore"):
            middle_control = (
                3.0 * np.diff(position, axis=0) / interval_duration[:, None]
                - velocity[:-1]
                - velocity[1:]
            )
            middle_speed_squared = np.einsum(
                "ij,ij->i", middle_control, middle_control
            )
        interval_index = int(np.argmax(middle_speed_squared))
        middle_record = {
            "hole": hole,
            "interval_index": interval_index,
            "time_start": float(time[interval_index]),
            "time_end": float(time[interval_index + 1]),
            "middle_control": middle_control[interval_index].tolist(),
            "speed_squared": float(middle_speed_squared[interval_index]),
            "speed": float(np.sqrt(middle_speed_squared[interval_index])),
        }
        middle_records.append(middle_record)
        invalid_intervals = np.flatnonzero(
            ~np.isfinite(middle_speed_squared) | (middle_speed_squared >= 1.0)
        )
        if invalid_intervals.size:
            interval_index = int(invalid_intervals[0])
            raise ValueError(
                f"trajectory interval [{time[interval_index]:.17g}, "
                f"{time[interval_index + 1]:.17g}] for hole {hole} cannot "
                "certify subluminal Hermite interpolation; middle velocity "
                f"control has |v|^2={middle_speed_squared[interval_index]:.17g}"
            )

        hole_records[str(hole)] = {
            "max_endpoint": {
                key: value for key, value in endpoint_record.items() if key != "hole"
            },
            "max_middle_control": {
                key: value for key, value in middle_record.items() if key != "hole"
            },
        }

    worst_endpoint = max(
        endpoint_records, key=lambda record: float(record["speed_squared"])
    )
    worst_middle = max(
        middle_records, key=lambda record: float(record["speed_squared"])
    )
    return {
        "formula_version": "dynbbh-quadratic-bezier-controls-v1",
        "certified_subluminal": True,
        "interval_indexing": "zero-based intervals [row i, row i+1]",
        "endpoint_max": worst_endpoint,
        "middle_control_max": worst_middle,
        "holes": hole_records,
    }


def _scaled_close_norm(value: np.ndarray, scale: float, tolerance: float) -> bool:
    """Apply one transparent absolute/relative norm bound."""
    return bool(np.linalg.norm(value) <= tolerance * max(1.0, float(scale)))


def _validate_q1_nonspinning(
    state: dict[str, np.ndarray], tolerance: float
) -> dict[str, object]:
    """Validate the documented equal-mass, nonspinning reference baseline."""
    initial_positions = np.asarray(
        [[state[f"{component}{hole}"][0] for component in "xyz"]
         for hole in (1, 2)]
    )
    initial_velocities = np.asarray(
        [[state[f"v{component}{hole}"][0] for component in "xyz"]
         for hole in (1, 2)]
    )
    initial_spin = np.asarray(
        [[state[f"a{hole}{component}"][0] for component in "xyz"]
         for hole in (1, 2)]
    )
    initial_masses = np.asarray(
        [state["m1_full"][0], state["m2_full"][0]]
    )
    initial_total_mass = float(np.sum(initial_masses))
    initial_separation = float(
        np.linalg.norm(initial_positions[0] - initial_positions[1])
    )
    mass_scale = max(1.0, initial_total_mass)
    if abs(float(initial_masses[0] - initial_masses[1])) > tolerance * mass_scale:
        raise ValueError(
            "q1-nonspinning validation failed: initial component masses are unequal"
        )
    center_of_mass = (
        initial_masses[:, None] * initial_positions
    ).sum(axis=0) / initial_total_mass
    center_velocity = (
        initial_masses[:, None] * initial_velocities
    ).sum(axis=0) / initial_total_mass
    if not _scaled_close_norm(center_of_mass, initial_separation, tolerance):
        raise ValueError(
            "q1-nonspinning validation failed: initial center of mass is not at "
            "the coordinate origin"
        )
    if not _scaled_close_norm(center_velocity, 1.0, tolerance):
        raise ValueError(
            "q1-nonspinning validation failed: initial center-of-mass velocity "
            "is nonzero"
        )
    if not _scaled_close_norm(initial_spin, 1.0, tolerance):
        raise ValueError(
            "q1-nonspinning validation failed: initial spin is nonzero"
        )

    final_positions = np.asarray(
        [[state[f"{component}{hole}"][-1] for component in "xyz"]
         for hole in (1, 2)]
    )
    final_velocities = np.asarray(
        [[state[f"v{component}{hole}"][-1] for component in "xyz"]
         for hole in (1, 2)]
    )
    final_spin = np.asarray(
        [[state[f"a{hole}{component}"][-1] for component in "xyz"]
         for hole in (1, 2)]
    )
    final_masses = np.asarray([state["m1_full"][-1], state["m2_full"][-1]])
    final_mass_scale = max(1.0, float(np.sum(final_masses)))
    inactive_component = bool(np.any(final_masses <= tolerance * final_mass_scale))
    duplicate_remnant = (
        _scaled_close_norm(
            final_positions[0] - final_positions[1],
            max(np.linalg.norm(final_positions[0]),
                np.linalg.norm(final_positions[1])),
            tolerance,
        )
        and _scaled_close_norm(
            final_velocities[0] - final_velocities[1], 1.0, tolerance
        )
        and _scaled_close_norm(
            final_spin[0] - final_spin[1],
            max(np.linalg.norm(final_spin[0]), np.linalg.norm(final_spin[1])),
            tolerance,
        )
    )
    if not inactive_component and not duplicate_remnant:
        raise ValueError(
            "q1-nonspinning validation failed: final remnant is neither two "
            "identical position/velocity/spin terms nor a single active term"
        )
    return {
        "profile": "q1-nonspinning",
        "tolerance": float(tolerance),
        "initial_center_of_mass": center_of_mass.tolist(),
        "initial_center_of_mass_velocity": center_velocity.tolist(),
        "remnant_representation": (
            "single-active-term" if inactive_component else "duplicate-identical-terms"
        ),
    }


def main() -> None:
    args = parse_args()
    manifest_path = None
    if not args.no_manifest:
        manifest_path = (
            args.manifest
            if args.manifest is not None
            else Path(f"{args.output}.manifest.json")
        )
    targets = [args.output] + ([] if manifest_path is None else [manifest_path])
    resolved_targets = [target.resolve() for target in targets]
    if len(set(resolved_targets)) != len(resolved_targets):
        raise ValueError("output and manifest paths must be distinct")
    if args.input.resolve() in resolved_targets:
        raise ValueError("input, output, and manifest paths must be distinct")
    for target in targets:
        if target.exists() and not args.force:
            raise FileExistsError(f"refusing to replace existing output {target}")
    if args.mass_unit is not None and (
        not np.isfinite(args.mass_unit) or args.mass_unit <= 0.0
    ):
        raise ValueError("--mass-unit must be finite and positive")
    if not np.isfinite(args.velocity_tolerance) or args.velocity_tolerance <= 0.0:
        raise ValueError("--velocity-tolerance must be finite and positive")
    if (
        not np.isfinite(args.validation_tolerance)
        or args.validation_tolerance <= 0.0
    ):
        raise ValueError("--validation-tolerance must be finite and positive")

    input_sha256 = _sha256(args.input)

    with h5py.File(args.input, "r") as source:
        if "t" not in source:
            raise ValueError("missing required dataset 't'")
        time = np.asarray(source["t"][()], dtype=np.float64)
        if time.ndim != 1 or time.size < 2:
            raise ValueError(
                "dataset 't' must be a one-dimensional array with >=2 rows"
            )
        if not np.all(np.isfinite(time)) or not np.all(np.diff(time) > 0.0):
            raise ValueError("trajectory times must be finite and strictly increasing")

        def series(name: str, *, required: bool = True) -> np.ndarray:
            if name not in source:
                if required:
                    raise ValueError(f"missing required dataset {name!r}")
                return np.zeros_like(time)
            values = np.asarray(source[name][()], dtype=np.float64)
            if values.ndim == 0:
                values = np.full_like(time, float(values))
            if values.shape != time.shape:
                raise ValueError(
                    f"dataset {name!r} has shape {values.shape}, expected {time.shape}"
                )
            if not np.all(np.isfinite(values)):
                raise ValueError(f"dataset {name!r} contains a non-finite value")
            return values

        state: dict[str, np.ndarray] = {}
        for name in STATE_NAMES[:12]:
            state[name] = series(name)

        full_mass_names = ("m1_full", "m2_full")
        if any(name in source for name in full_mass_names):
            if not all(name in source for name in full_mass_names):
                raise ValueError("m1_full and m2_full must be provided together")
            state["m1_full"] = series("m1_full")
            state["m2_full"] = series("m2_full")
        else:
            state["m1_full"] = series("m1")
            state["m2_full"] = series("m2")

        new_spin_names = STATE_NAMES[12:18]
        if any(name in source for name in new_spin_names):
            if not all(name in source for name in new_spin_names):
                raise ValueError(
                    "an a^i trajectory must provide all six a1*/a2* datasets"
                )
            for name in new_spin_names:
                state[name] = series(name)
        else:
            # Older pre-merger tables store dimensionless chi vectors as s1*/s2*.
            # Convert them to the length-valued Kerr parameter a^i = chi^i M.
            legacy_spin_names = tuple(
                f"s{hole}{component}" for hole in (1, 2) for component in "xyz"
            )
            if not all(name in source for name in legacy_spin_names):
                missing = [name for name in legacy_spin_names if name not in source]
                if not args.assume_missing_spin_zero:
                    raise ValueError(
                        "legacy spin data are incomplete; missing "
                        + ", ".join(missing)
                        + "; pass --assume-missing-spin-zero only after inspection"
                    )
            for component in "xyz":
                state[f"a1{component}"] = (
                    series(
                        f"s1{component}", required=not args.assume_missing_spin_zero
                    )
                    * state["m1_full"]
                )
                state[f"a2{component}"] = (
                    series(
                        f"s2{component}", required=not args.assume_missing_spin_zero
                    )
                    * state["m2_full"]
                )

    initial_total_mass = state["m1_full"][0] + state["m2_full"][0]
    mass_unit = initial_total_mass if args.mass_unit is None else args.mass_unit
    if not np.isfinite(mass_unit) or mass_unit <= 0.0:
        raise ValueError(
            "the selected trajectory mass unit must be finite and positive"
        )
    time = time / mass_unit
    dimensional_names = STATE_NAMES[:6] + STATE_NAMES[12:]
    for name in dimensional_names:
        state[name] = state[name] / mass_unit

    hermite_velocity_certificate = _hermite_velocity_certificate(time, state)

    hole_mismatches: dict[str, dict[str, object]] = {}
    velocity_mismatch: dict[str, object] = {
        "definition": "fractional Euclidean norm |dx/dt-v|/max(|dx/dt|,|v|,1e-12)",
        "tolerance": float(args.velocity_tolerance),
        "strict": bool(args.strict_velocity_mismatch),
        "row_indexing": "zero-based ASCII data rows; comments/header excluded",
        "evaluated_rows": (
            "interior rows; endpoint one-sided differences excluded"
            if time.size > 2
            else "all rows"
        ),
        "holes": hole_mismatches,
    }
    mismatch_records: list[dict[str, object]] = []
    for hole in (1, 2):
        position = np.column_stack(
            tuple(state[f"{component}{hole}"] for component in "xyz")
        )
        velocity = np.column_stack(
            tuple(state[f"v{component}{hole}"] for component in "xyz")
        )
        edge_order = 2 if time.size > 2 else 1
        numerical_velocity = np.column_stack(
            tuple(
                np.gradient(position[:, axis], time, edge_order=edge_order)
                for axis in range(3)
            )
        )
        scale = np.maximum(
            np.maximum(np.linalg.norm(velocity, axis=1),
                       np.linalg.norm(numerical_velocity, axis=1)),
            1.0e-12,
        )
        mismatch = np.linalg.norm(velocity - numerical_velocity, axis=1) / scale
        candidate_rows = (
            np.arange(1, mismatch.size - 1, dtype=np.int64)
            if mismatch.size > 2
            else np.arange(mismatch.size, dtype=np.int64)
        )
        candidate_mismatch = mismatch[candidate_rows]
        local_max = int(np.argmax(candidate_mismatch))
        row_index = int(candidate_rows[local_max])
        maximum = float(candidate_mismatch[local_max])
        record = {
            "hole": hole,
            "max_fractional": maximum,
            "row_index": row_index,
            "time": float(time[row_index]),
        }
        mismatch_records.append(record)
        hole_mismatches[str(hole)] = {
            key: value for key, value in record.items() if key != "hole"
        }
        if maximum > args.velocity_tolerance and args.strict_velocity_mismatch:
            raise ValueError(
                f"hole {hole}: max fractional |dx/dt-v| mismatch {maximum:.3e} "
                f"at row {row_index} exceeds tolerance "
                f"{args.velocity_tolerance:.3e}"
            )
        if maximum > args.velocity_tolerance:
            warnings.warn(
                f"hole {hole}: max fractional |dx/dt-v| mismatch is "
                f"{maximum:.3e} at row {row_index}; inspect trajectory units and "
                "merger smoothing",
                stacklevel=1,
            )

    overall_mismatch = max(
        mismatch_records, key=lambda record: float(record["max_fractional"])
    )
    velocity_mismatch["overall_max"] = {
        key: value for key, value in overall_mismatch.items()
    }

    values = np.column_stack((time, *(state[name] for name in STATE_NAMES)))
    if np.any(values[:, 19] < 0.0) or np.any(values[:, 20] < 0.0):
        raise ValueError("individual trajectory masses must be non-negative")
    if np.any(values[:, 19] + values[:, 20] <= 0.0):
        raise ValueError("the total trajectory mass must be positive in every row")

    validation: dict[str, object]
    if args.validation_profile == "q1-nonspinning":
        validation = _validate_q1_nonspinning(state, args.validation_tolerance)
    else:
        validation = {"profile": "none"}

    header = (
        f"input_mass_unit = {mass_unit:.17e}; dimensional columns normalized "
        "by this value\n"
        + "t "
        + " ".join(STATE_NAMES)
    )
    for target in targets:
        target.parent.mkdir(parents=True, exist_ok=True)

    staged_output = _temporary_path(args.output)
    staged_manifest = None
    try:
        with staged_output.open("w", encoding="utf-8", newline="\n") as output_file:
            np.savetxt(output_file, values, fmt="%.17e", header=header)
            output_file.flush()
            os.fsync(output_file.fileno())
        output_sha256 = _sha256(staged_output)

        staged_targets = [(staged_output, args.output)]
        if manifest_path is not None:
            manifest = {
                "schema_version": MANIFEST_SCHEMA_VERSION,
                "artifact": "athenak-bbh-trajectory-ascii",
                "input_sha256": input_sha256,
                "output_sha256": output_sha256,
                "columns": ["t", *STATE_NAMES],
                "row_count": int(time.size),
                "time_range": {
                    "start": float(time[0]),
                    "end": float(time[-1]),
                    "unit": "normalized mass time M",
                },
                "mass_unit": {
                    "input_value": float(mass_unit),
                    "selection": (
                        "initial_total_mass"
                        if args.mass_unit is None
                        else "explicit"
                    ),
                    "output_convention": "G=c=M_unit=1",
                },
                "initial_state": _endpoint_summary(0, time, state),
                "final_state": _endpoint_summary(-1, time, state),
                "velocity_mismatch": velocity_mismatch,
                "hermite_velocity_certificate": hermite_velocity_certificate,
                "validation": validation,
            }
            staged_manifest = _temporary_path(manifest_path)
            with staged_manifest.open(
                "w", encoding="utf-8", newline="\n"
            ) as manifest_file:
                json.dump(
                    manifest,
                    manifest_file,
                    sort_keys=True,
                    indent=2,
                    allow_nan=False,
                )
                manifest_file.write("\n")
                manifest_file.flush()
                os.fsync(manifest_file.fileno())
            staged_targets.append((staged_manifest, manifest_path))

        _publish_staged(staged_targets, args.force)
    except BaseException:
        for staged in (staged_output, staged_manifest):
            if staged is not None:
                try:
                    staged.unlink()
                except FileNotFoundError:
                    pass
        raise
    print(f"wrote {time.size} rows to {args.output}")
    if manifest_path is not None:
        print(f"wrote manifest to {manifest_path}")


if __name__ == "__main__":
    main()
