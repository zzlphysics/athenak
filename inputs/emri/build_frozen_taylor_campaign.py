#!/usr/bin/env python3
"""Build and audit a frozen-snapshot EMRI Taylor-wind campaign.

Each case selects one global GRMHD dump, a circular Kerr orbit radius and an
azimuth.  The dump is loaded only once even when several fitting radii or cases
reuse it.  The output combines local Taylor-model gates, physical coherence
scales, a BHL hierarchy/cost decision, and exact AthenaK command-line overrides.
It deliberately does not launch an evolution.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

import analyze_taylor_profile_convergence as convergence
import build_kerr_circular_frame as circular
import extract_static_taylor_worldtube as static
import plan_bhl_hierarchy as hierarchy


REQUEST_CLASSIFICATION = "athenak-emri-frozen-taylor-request-v1"
CLASSIFICATION = "athenak-emri-frozen-taylor-campaign-v1"
CASE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
DEFAULT_GATES = {
    "minimum_fit_samples": 16,
    "maximum_density_log_fit_rms": 0.03,
    "maximum_pressure_log_fit_rms": 0.03,
    "maximum_velocity_fit_relative_rms": 0.03,
    "maximum_magnetic_fit_relative_rms": 0.05,
    "maximum_magnetic_gradient_trace": 1.0e-12,
    "maximum_cell_volume_ratio": 1.000001,
    "maximum_fitting_scale_sensitivity": 0.15,
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite_positive(value: object, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def _finite(value: object, name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _maximum_condition(observed: float, threshold: float) -> dict[str, object]:
    return {
        "relation": "maximum",
        "observed": observed,
        "threshold": threshold,
        "passed": math.isfinite(observed) and observed <= threshold,
    }


def _minimum_condition(observed: float, threshold: float) -> dict[str, object]:
    return {
        "relation": "minimum",
        "observed": observed,
        "threshold": threshold,
        "passed": math.isfinite(observed) and observed >= threshold,
    }


def case_orbit(
    primary_mass: float,
    dimensionless_spin: float,
    boyer_lindquist_radius: float,
    direction: int,
    phase: float,
) -> dict[str, object]:
    """Return a frozen circular-orbit anchor and coordinate velocity."""

    primary_mass = _finite_positive(primary_mass, "primary mass")
    dimensionless_spin = _finite(dimensionless_spin, "primary spin")
    radius = _finite_positive(boyer_lindquist_radius, "orbital radius")
    phase = _finite(phase, "orbital phase")
    if abs(dimensionless_spin) > 1.0 or direction not in (-1, 1):
        raise ValueError("Kerr spin or orbit direction is invalid")
    isco = circular.kerr_isco(primary_mass, dimensionless_spin, direction)
    if radius <= isco:
        raise ValueError(
            f"orbital radius {radius:.17g} is not outside Kerr ISCO {isco:.17g}"
        )
    spin = primary_mass * dimensionless_spin
    omega = circular.circular_kerr_omega(
        primary_mass, spin, radius, direction
    )
    coordinate_radius = math.sqrt(radius**2 + spin**2)
    radial = np.asarray((math.cos(phase), math.sin(phase), 0.0))
    tangent = np.asarray((-math.sin(phase), math.cos(phase), 0.0))
    return {
        "boyer_lindquist_radius": radius,
        "cartesian_equatorial_radius": coordinate_radius,
        "phase": phase,
        "coordinate_angular_frequency": omega,
        "kerr_isco_radius": isco,
        "anchor": (coordinate_radius * radial).tolist(),
        "source_velocity": (omega * coordinate_radius * tangent).tolist(),
    }


def gradient_coherence_scales(
    parameters: dict[str, float],
    gamma: float,
    angular_frequency_local: float,
) -> tuple[dict[str, float], dict[str, object]]:
    """Derive conservative local Taylor and disk-height scale proxies."""

    def vector(prefix: str) -> np.ndarray:
        return np.asarray(
            [parameters[f"{prefix}{index}"] for index in range(1, 4)],
            dtype=np.float64,
        )

    velocity = vector("u")
    magnetic = vector("b")
    density_gradient = vector("dlnrho_dxh")
    pressure_gradient = vector("dlnpgas_dxh")
    velocity_gradient = np.asarray(
        [
            [parameters[f"du{i}_dxh{j}"] for j in range(1, 4)]
            for i in range(1, 4)
        ],
        dtype=np.float64,
    )
    magnetic_gradient = np.asarray(
        [
            [parameters[f"db{i}_dxh{j}"] for j in range(1, 4)]
            for i in range(1, 4)
        ],
        dtype=np.float64,
    )
    scales: dict[str, float] = {}

    def inverse_norm(name: str, values: np.ndarray) -> None:
        norm = float(np.linalg.norm(values))
        if norm > 0.0:
            scales[name] = 1.0 / norm

    def amplitude_over_gradient(
        name: str, amplitude: np.ndarray, gradient: np.ndarray
    ) -> None:
        numerator = float(np.linalg.norm(amplitude))
        denominator = float(np.linalg.norm(gradient))
        if numerator > 0.0 and denominator > 0.0:
            scales[name] = numerator / denominator

    inverse_norm("L_rho", density_gradient)
    inverse_norm("L_pgas", pressure_gradient)
    amplitude_over_gradient("L_velocity", velocity, velocity_gradient)
    amplitude_over_gradient("L_magnetic", magnetic, magnetic_gradient)
    state = hierarchy.characteristic_state(
        parameters["rho0"],
        parameters["pgas0"],
        gamma,
        velocity,
        magnetic,
    )
    omega = abs(_finite(angular_frequency_local, "local angular frequency"))
    if omega > 0.0:
        scales["disk_H_hydrostatic_proxy"] = (
            math.sqrt(float(state["sound_speed_squared"])) / omega
        )
    return scales, state


def fitting_scale_comparison(
    profiles: list[dict[str, object]],
) -> tuple[float, dict[str, object]]:
    if len(profiles) < 2:
        raise ValueError("each frozen case requires at least two fitting radii")
    order = list(static.PROFILE_PARAMETER_ORDER)
    values = np.asarray(
        [[float(profile["parameters"][name]) for name in order]
         for profile in profiles],
        dtype=np.float64,
    )
    indices = convergence._group_indices()
    detail: dict[str, object] = {}
    maximum = 0.0
    for candidate, profile in enumerate(profiles):
        label = f"r{float(profile['fit_radius_source']):g}"
        groups = {
            group: convergence.group_error(
                values[:1], values[candidate:candidate + 1], columns
            )
            for group, columns in indices.items()
        }
        detail[label] = {
            "reference": candidate == 0,
            "fit_radius_source": float(profile["fit_radius_source"]),
            "groups": groups,
        }
        if candidate > 0:
            maximum = max(
                maximum,
                *(float(group["symmetric_relative_l2"])
                  for group in groups.values()),
            )
    return maximum, detail


def _validated_gates(document: dict[str, object]) -> dict[str, float]:
    requested = document.get("gates", {})
    if not isinstance(requested, dict):
        raise ValueError("gates must be an object")
    unknown = sorted(set(requested).difference(DEFAULT_GATES))
    if unknown:
        raise ValueError(f"unknown frozen-campaign gates: {', '.join(unknown)}")
    result = dict(DEFAULT_GATES)
    result.update(requested)
    for name, value in result.items():
        result[name] = _finite_positive(value, name)
    return result


def _validated_request(
    document: object, request_path: Path
) -> dict[str, object]:
    if not isinstance(document, dict) or document.get(
        "classification"
    ) != REQUEST_CLASSIFICATION:
        raise ValueError("request classification is unsupported")
    primary = document.get("primary")
    if not isinstance(primary, dict):
        raise ValueError("request primary must be an object")
    primary_mass = _finite_positive(primary.get("mass"), "primary mass")
    primary_chi = _finite(primary.get("dimensionless_spin"), "primary spin")
    if abs(primary_chi) > 1.0:
        raise ValueError("primary dimensionless spin must have magnitude <= 1")
    direction = int(primary.get("orbit_direction", 1))
    if direction not in (-1, 1):
        raise ValueError("orbit direction must be -1 or 1")
    mass_ratio = _finite_positive(document.get("mass_ratio"), "mass ratio")
    if mass_ratio >= 1.0:
        raise ValueError("EMRI mass ratio must be smaller than one")
    density_renormalization = _finite_positive(
        document.get("density_renormalization"), "density renormalization"
    )
    cases = document.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("request cases must be a nonempty list")
    normalized_cases = []
    identifiers: set[str] = set()
    for raw in cases:
        if not isinstance(raw, dict):
            raise ValueError("every campaign case must be an object")
        identifier = str(raw.get("id", ""))
        if CASE_ID.fullmatch(identifier) is None or identifier in identifiers:
            raise ValueError(f"case id is invalid or duplicated: {identifier!r}")
        identifiers.add(identifier)
        state_text = raw.get("state")
        if not isinstance(state_text, str) or not state_text:
            raise ValueError(f"case {identifier} has no state path")
        state_path = Path(state_text).expanduser()
        if not state_path.is_absolute():
            state_path = request_path.parent / state_path
        state_path = state_path.resolve(strict=True)
        radii = raw.get("fit_radii")
        if not isinstance(radii, list) or len(radii) < 2:
            raise ValueError(f"case {identifier} requires at least two fit radii")
        fit_radii = [
            _finite_positive(value, f"case {identifier} fit radius")
            for value in radii
        ]
        if fit_radii != sorted(set(fit_radii)):
            raise ValueError(
                f"case {identifier} fit radii must be unique and increasing"
            )
        normalized_cases.append(
            {
                "id": identifier,
                "state": state_path,
                "orbital_radius": _finite_positive(
                    raw.get("orbital_radius"),
                    f"case {identifier} orbital radius",
                ),
                "phase": _finite(raw.get("phase"), f"case {identifier} phase"),
                "fit_radii": fit_radii,
            }
        )
    return {
        "primary_mass": primary_mass,
        "primary_chi": primary_chi,
        "orbit_direction": direction,
        "mass_ratio": mass_ratio,
        "density_renormalization": density_renormalization,
        "gates": _validated_gates(document),
        "cases": normalized_cases,
    }


def _fit_one_case(
    case: dict[str, object],
    state: dict,
    thermodynamics: dict[str, object],
    state_sha256: str,
    request: dict[str, object],
) -> dict[str, object]:
    primary_mass = float(request["primary_mass"])
    primary_chi = float(request["primary_chi"])
    direction = int(request["orbit_direction"])
    mass_ratio = float(request["mass_ratio"])
    length_scale = primary_mass / mass_ratio
    density_renormalization = float(request["density_renormalization"])
    orbit = case_orbit(
        primary_mass,
        primary_chi,
        float(case["orbital_radius"]),
        direction,
        float(case["phase"]),
    )
    anchor = np.asarray(orbit["anchor"], dtype=np.float64)
    source_velocity = np.asarray(orbit["source_velocity"], dtype=np.float64)
    primary_center = np.zeros(3)
    disk_normal = np.asarray((0.0, 0.0, 1.0))
    basis = static.canonical_spatial_basis(anchor, primary_center, disk_normal)
    gamma_adm, alpha, beta = static.analytic_kerr_anchor_adm(
        anchor, primary_mass, primary_chi
    )
    metric = static.spacetime_metric_from_adm(gamma_adm, alpha, beta)
    tetrad, coframe = static.build_source_tetrad(
        metric, source_velocity, basis
    )
    profiles = []
    for fit_radius in case["fit_radii"]:
        cloud = static.collect_local_samples(
            state,
            None,
            None,
            anchor,
            coframe,
            float(fit_radius),
            analytic_kerr=(primary_mass, primary_chi),
        )
        parameters, diagnostics = static.fit_static_profile(
            cloud, coframe, float(fit_radius)
        )
        parameters = static.rescale_profile_parameters(
            parameters, length_scale, density_renormalization
        )
        diagnostics = static.rescale_profile_diagnostics(
            diagnostics, length_scale, density_renormalization
        )
        profiles.append(
            {
                "fit_radius_source": float(fit_radius),
                "parameters": parameters,
                "fit_diagnostics": diagnostics,
            }
        )
    maximum_sensitivity, sensitivity = fitting_scale_comparison(profiles)
    reference = profiles[0]
    diagnostics = reference["fit_diagnostics"]
    gates = request["gates"]
    assert isinstance(gates, dict)
    conditions = {
        "minimum_fit_samples": _minimum_condition(
            min(float(profile["fit_diagnostics"]["sample_count"])
                for profile in profiles),
            float(gates["minimum_fit_samples"]),
        ),
        "density_log_fit_rms": _maximum_condition(
            float(diagnostics["density_log_weighted_rms"]),
            float(gates["maximum_density_log_fit_rms"]),
        ),
        "pressure_log_fit_rms": _maximum_condition(
            float(diagnostics["pressure_log_weighted_rms"]),
            float(gates["maximum_pressure_log_fit_rms"]),
        ),
        "velocity_fit_relative_rms": _maximum_condition(
            float(diagnostics["velocity_relative_rms"]),
            float(gates["maximum_velocity_fit_relative_rms"]),
        ),
        "magnetic_fit_relative_rms": _maximum_condition(
            float(diagnostics["magnetic_relative_rms"]),
            float(gates["maximum_magnetic_fit_relative_rms"]),
        ),
        "magnetic_gradient_trace": _maximum_condition(
            max(abs(float(profile["fit_diagnostics"]["magnetic_gradient_trace"]))
                for profile in profiles),
            float(gates["maximum_magnetic_gradient_trace"]),
        ),
        "cell_volume_ratio": _maximum_condition(
            max(
                float(profile["fit_diagnostics"]["maximum_cell_volume"])
                / float(profile["fit_diagnostics"]["minimum_cell_volume"])
                for profile in profiles
            ),
            float(gates["maximum_cell_volume_ratio"]),
        ),
        "fitting_scale_sensitivity": _maximum_condition(
            maximum_sensitivity,
            float(gates["maximum_fitting_scale_sensitivity"]),
        ),
    }
    local_model_passed = all(
        bool(condition["passed"]) for condition in conditions.values()
    )
    source_gamma = thermodynamics.get("adiabatic_index")
    if source_gamma is None:
        raise RuntimeError(
            "frozen-campaign planning requires an explicit source adiabatic index"
        )
    gamma = float(source_gamma)
    omega_local = float(orbit["coordinate_angular_frequency"]) / length_scale
    parameters = reference["parameters"]
    coherence, characteristic = gradient_coherence_scales(
        parameters, gamma, omega_local
    )
    bhl_plan = hierarchy.build_plan(
        secondary_mass=1.0,
        primary_mass=primary_mass / mass_ratio,
        orbital_radius=float(case["orbital_radius"]) / mass_ratio,
        rho=float(parameters["rho0"]),
        pgas=float(parameters["pgas0"]),
        gamma=gamma,
        spatial_four_velocity=(parameters["u1"], parameters["u2"], parameters["u3"]),
        magnetic_field=(parameters["b1"], parameters["b2"], parameters["b3"]),
        coherence_scales=coherence,
    )
    overrides = [
        "adm/dynamic=true",
        "time/subcycling=level",
        "problem/background_mode=full",
        f"problem/primary_mass={primary_mass / mass_ratio:.17g}",
        "problem/secondary_mass=1",
        f"problem/primary_chi={primary_chi:.17g}",
        f"problem/orbital_radius={float(case['orbital_radius']) / mass_ratio:.17g}",
        f"problem/orbit_direction={direction}",
        f"mhd/gamma={gamma:.17g}",
    ]
    overrides.extend(
        f"problem/{name}={float(parameters[name]):.17g}"
        for name in static.PROFILE_PARAMETER_ORDER
    )
    hierarchy_recommendation = str(bhl_plan["recommendation"])
    execution_ready = local_model_passed and hierarchy_recommendation == "direct_grmhd"
    return {
        "id": case["id"],
        "state_file": str(case["state"]),
        "state_sha256": state_sha256,
        "time_global": float(state["time"]),
        "cycle": int(state["cycle"]),
        "state_thermodynamics": thermodynamics,
        "orbit": orbit,
        "global_length_in_local_units": length_scale,
        "density_renormalization": density_renormalization,
        "source_tetrad_contravariant": tetrad.tolist(),
        "source_tetrad_coframe": coframe.tolist(),
        "profiles": profiles,
        "fitting_scale_sensitivity": sensitivity,
        "coherence_scales_local": coherence,
        "characteristic_state": characteristic,
        "bhl_hierarchy_plan": bhl_plan,
        "athena_overrides": overrides,
        "assessment": {
            "local_taylor_model_passed": local_model_passed,
            "execution_ready": execution_ready,
            "conditions": conditions,
            "hierarchy_recommendation": hierarchy_recommendation,
            "blocking_reason": (
                None
                if execution_ready
                else (
                    "local Taylor gates failed"
                    if not local_model_passed
                    else "recommended outer-inner route is not yet production-ready"
                )
            ),
        },
    }


def build_campaign(request: dict[str, object]) -> dict[str, object]:
    loaded: dict[Path, tuple[dict, dict[str, object], str]] = {}
    cases = []
    expected_spin = float(request["primary_mass"]) * float(request["primary_chi"])
    for case in request["cases"]:
        path = case["state"]
        assert isinstance(path, Path)
        if path not in loaded:
            state, thermodynamics = static.read_state_dump(path)
            try:
                stored_spin = float(
                    static._header_value(state["header"], "<coord>", "a")
                )
            except (KeyError, ValueError) as error:
                raise RuntimeError(f"state dump {path} has no Kerr spin header") from error
            tolerance = 64.0 * np.finfo(float).eps * max(1.0, abs(expected_spin))
            if abs(stored_spin - expected_spin) > tolerance:
                raise RuntimeError(
                    f"state dump {path} Kerr spin does not match the request"
                )
            loaded[path] = (state, thermodynamics, _sha256(path))
        state, thermodynamics, digest = loaded[path]
        cases.append(
            _fit_one_case(case, state, thermodynamics, digest, request)
        )
    return {
        "classification": CLASSIFICATION,
        "run_status": "completed",
        "science_ready": False,
        "primary": {
            "mass_global": request["primary_mass"],
            "dimensionless_spin": request["primary_chi"],
            "orbit_direction": request["orbit_direction"],
        },
        "mass_ratio": request["mass_ratio"],
        "density_renormalization": request["density_renormalization"],
        "gates": request["gates"],
        "cases": cases,
        "summary": {
            "case_count": len(cases),
            "local_taylor_pass_count": sum(
                bool(case["assessment"]["local_taylor_model_passed"])
                for case in cases
            ),
            "execution_ready_count": sum(
                bool(case["assessment"]["execution_ready"])
                for case in cases
            ),
            "recommendation_counts": {
                name: sum(
                    case["assessment"]["hierarchy_recommendation"] == name
                    for case in cases
                )
                for name in sorted(
                    {case["assessment"]["hierarchy_recommendation"] for case in cases}
                )
            },
        },
        "science_blockers": [
            "all fit-radius variants share their parent global source resolution",
            "the frozen one-way Taylor model omits secondary feedback on the global disk",
            "an outer-inner recommendation is not executable science until characteristic fluid replay and matching-radius convergence are qualified",
        ],
    }


def _write_overrides(output_directory: Path, campaign: dict[str, object]) -> None:
    overrides = output_directory / "overrides"
    overrides.mkdir()
    for case in campaign["cases"]:
        path = overrides / f"{case['id']}.txt"
        path.write_text("\n".join(case["athena_overrides"]) + "\n", encoding="utf-8")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    request_path = arguments.request.expanduser().resolve(strict=True)
    document = json.loads(request_path.read_text(encoding="utf-8"))
    request = _validated_request(document, request_path)
    output = arguments.output_directory.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.mkdir(parents=True)
    campaign = build_campaign(request)
    campaign["request_file"] = str(request_path)
    campaign["request_sha256"] = _sha256(request_path)
    campaign_path = output / "campaign.json"
    campaign_path.write_text(
        json.dumps(campaign, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    _write_overrides(output, campaign)
    print(campaign_path)
    print(json.dumps(campaign["summary"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
