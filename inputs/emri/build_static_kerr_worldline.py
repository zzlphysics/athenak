#!/usr/bin/env python3
"""Build a circular worldline manifest for state-only static-Kerr GRMHD dumps."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import build_kerr_circular_frame as circular
import extract_static_taylor_worldtube as static


CLASSIFICATION = "athenak-emri-worldline-v1"


def circular_samples(
    records: list[tuple[Path, float, int]],
    primary_mass: float,
    dimensionless_spin: float,
    orbital_radius: float,
    direction: int,
    initial_phase: float,
    phase_reference_time: float | None = None,
) -> tuple[list[dict[str, object]], dict[str, float | int]]:
    if len(records) < 2:
        raise ValueError("a worldline manifest requires at least two state dumps")
    if (
        not math.isfinite(primary_mass)
        or primary_mass <= 0.0
        or not math.isfinite(dimensionless_spin)
        or abs(dimensionless_spin) > 1.0
        or not math.isfinite(orbital_radius)
        or orbital_radius <= 0.0
        or direction not in (-1, 1)
        or not math.isfinite(initial_phase)
    ):
        raise ValueError("circular static-Kerr worldline parameters are invalid")
    ordered = sorted(records, key=lambda entry: entry[1])
    times = np.asarray([entry[1] for entry in ordered], dtype=np.float64)
    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0.0):
        raise ValueError("state dump times must increase strictly")
    reference_time = (
        float(times[0])
        if phase_reference_time is None
        else float(phase_reference_time)
    )
    if not math.isfinite(reference_time):
        raise ValueError("phase reference time must be finite")
    primary_spin = primary_mass * dimensionless_spin
    isco = circular.kerr_isco(primary_mass, dimensionless_spin, direction)
    if orbital_radius <= isco:
        raise ValueError(
            f"orbital radius {orbital_radius:.17g} is not outside Kerr ISCO "
            f"{isco:.17g}"
        )
    omega = circular.circular_kerr_omega(
        primary_mass, primary_spin, orbital_radius, direction
    )
    coordinate_radius = math.sqrt(orbital_radius**2 + primary_spin**2)
    samples = []
    for path, time, cycle in ordered:
        phase = initial_phase + omega * (time - reference_time)
        radial = np.asarray((math.cos(phase), math.sin(phase), 0.0))
        tangent = np.asarray((-math.sin(phase), math.cos(phase), 0.0))
        samples.append(
            {
                "state": str(path),
                "time": float(time),
                "cycle": int(cycle),
                "anchor": (coordinate_radius * radial).tolist(),
                "source_velocity": (omega * coordinate_radius * tangent).tolist(),
            }
        )
    return samples, {
        "primary_mass": primary_mass,
        "dimensionless_spin": dimensionless_spin,
        "orbit_direction": direction,
        "boyer_lindquist_radius": orbital_radius,
        "cartesian_equatorial_radius": coordinate_radius,
        "coordinate_angular_frequency": omega,
        "kerr_isco_radius": isco,
        "initial_phase": initial_phase,
        "phase_reference_time": reference_time,
        "phase_advance_per_median_dump": omega * float(np.median(np.diff(times))),
    }


def _metadata(path: Path) -> tuple[Path, float, int, dict[str, object]]:
    resolved = path.expanduser().resolve(strict=True)
    data = static.bin_convert.read_binary(str(resolved), variables=())
    required = {"dens", "velx", "vely", "velz", "bcc1", "bcc2", "bcc3"}
    missing = sorted(required.difference(data["var_names"]))
    if missing or not ({"press", "eint"} & set(data["var_names"])):
        detail = ", ".join(missing) if missing else "press or eint"
        raise RuntimeError(f"state dump {resolved} is missing {detail}")
    try:
        general_rel = static._header_value(
            data["header"], "<coord>", "general_rel"
        ).lower()
        stored_spin = float(static._header_value(data["header"], "<coord>", "a"))
    except (KeyError, ValueError) as error:
        raise RuntimeError(
            f"state dump {resolved} has no static Kerr coordinate contract"
        ) from error
    if general_rel not in ("true", "1"):
        raise RuntimeError(f"state dump {resolved} is not general relativistic")
    return resolved, float(data["time"]), int(data["cycle"]), {
        "stored_spin": stored_spin,
        "available_leaf_levels": sorted(
            set(int(value) for value in data["mb_logical"][:, 3])
        ),
        "root_shape": [int(data[name]) for name in ("Nx1", "Nx2", "Nx3")],
        "meshblock_shape": [
            int(data[name]) for name in ("nx1_mb", "nx2_mb", "nx3_mb")
        ],
        "thermodynamic_variable": (
            "press" if "press" in data["var_names"] else "eint"
        ),
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--states", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--primary-mass", type=float, default=1.0)
    parser.add_argument("--primary-chi", type=float, required=True)
    parser.add_argument("--orbital-radius", type=float, required=True)
    parser.add_argument("--orbit-direction", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--initial-phase", type=float, default=0.0)
    parser.add_argument("--phase-reference-time", type=float)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    metadata = [_metadata(path) for path in arguments.states]
    stored_spins = np.asarray(
        [float(entry[3]["stored_spin"]) for entry in metadata]
    )
    expected_spin = arguments.primary_mass * arguments.primary_chi
    tolerance = 64.0 * np.finfo(float).eps * max(1.0, abs(expected_spin))
    if np.max(np.abs(stored_spins - expected_spin)) > tolerance:
        raise RuntimeError(
            "requested primary spin does not match the state-dump Kerr header"
        )
    topologies = [
        {
            key: entry[3][key]
            for key in (
                "available_leaf_levels",
                "root_shape",
                "meshblock_shape",
                "thermodynamic_variable",
            )
        }
        for entry in metadata
    ]
    if any(topology != topologies[0] for topology in topologies[1:]):
        raise RuntimeError("state-dump topology or primitive convention changes")
    records = [(entry[0], entry[1], entry[2]) for entry in metadata]
    samples, orbit = circular_samples(
        records,
        arguments.primary_mass,
        arguments.primary_chi,
        arguments.orbital_radius,
        arguments.orbit_direction,
        arguments.initial_phase,
        arguments.phase_reference_time,
    )
    document = {
        "classification": CLASSIFICATION,
        "source_provenance": {
            "metric_content": "primary_only",
            "fluid_content": "global_grmhd",
            "secondary_backreaction": "absent",
        },
        "metric_source": static.analytic_kerr_metric_source(
            arguments.primary_mass, arguments.primary_chi
        ),
        "state_topology": topologies[0],
        "orbit": orbit,
        "samples": samples,
    }
    output = arguments.output.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
