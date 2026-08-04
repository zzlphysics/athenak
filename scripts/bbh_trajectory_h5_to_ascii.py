#!/usr/bin/env python3
"""Convert AnalyticalBBH trajectory HDF5 data to AthenaK's portable table format."""

from __future__ import annotations

import argparse
from pathlib import Path
import warnings

import h5py
import numpy as np


STATE_NAMES = (
    "x1", "y1", "z1", "x2", "y2", "z2",
    "vx1", "vy1", "vz1", "vx2", "vy2", "vz2",
    "a1x", "a1y", "a1z", "a2x", "a2y", "a2z",
    "m1_full", "m2_full",
)


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
        "--force",
        action="store_true",
        help="replace an existing output file",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output.exists() and not args.force:
        raise FileExistsError(f"refusing to replace existing output {args.output}")
    if args.mass_unit is not None and (
        not np.isfinite(args.mass_unit) or args.mass_unit <= 0.0
    ):
        raise ValueError("--mass-unit must be finite and positive")
    if not np.isfinite(args.velocity_tolerance) or args.velocity_tolerance <= 0.0:
        raise ValueError("--velocity-tolerance must be finite and positive")

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
        interior = mismatch[1:-1] if mismatch.size > 2 else mismatch
        if np.max(interior) > args.velocity_tolerance:
            warnings.warn(
                f"hole {hole}: max fractional |dx/dt-v| mismatch is "
                f"{np.max(interior):.3e}; inspect trajectory units and "
                "merger smoothing",
                stacklevel=1,
            )

    values = np.column_stack((time, *(state[name] for name in STATE_NAMES)))
    speed1_squared = sum(values[:, index] ** 2 for index in (7, 8, 9))
    speed2_squared = sum(values[:, index] ** 2 for index in (10, 11, 12))
    if np.any(speed1_squared >= 1.0) or np.any(speed2_squared >= 1.0):
        raise ValueError("trajectory contains a superluminal black-hole velocity")
    if np.any(values[:, 19] < 0.0) or np.any(values[:, 20] < 0.0):
        raise ValueError("individual trajectory masses must be non-negative")
    if np.any(values[:, 19] + values[:, 20] <= 0.0):
        raise ValueError("the total trajectory mass must be positive in every row")

    header = (
        f"input_mass_unit = {mass_unit:.17e}; dimensional columns normalized "
        "by this value\n"
        + "t "
        + " ".join(STATE_NAMES)
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(args.output, values, fmt="%.17e", header=header)
    print(f"wrote {time.size} rows to {args.output}")


if __name__ == "__main__":
    main()
