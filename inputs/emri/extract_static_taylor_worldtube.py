#!/usr/bin/env python3
"""Fit a static source-tetrad EMRI wind profile from AthenaK GRMHD dumps.

This is the first, deliberately low-order global-to-local coupling path.  It reads
co-temporal ``mhd_w_bcc`` and ``adm`` leaf-MeshBlock dumps, constructs the secondary
source tetrad at a requested worldline point, and volume-weight fits the local value and
first spatial derivative of the GRMHD state.  The magnetic fit is constrained to be
divergence-free.  Outputs are an AthenaK ``<problem>`` fragment and a provenance JSON
manifest; the profile is consumed by ``emri_windtunnel``'s analytic-gradient interface.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
VIS_PYTHON = ROOT / "vis" / "python"
if str(VIS_PYTHON) not in sys.path:
    sys.path.insert(0, str(VIS_PYTHON))

import bin_convert  # noqa: E402


STATE_VARIABLES = (
    "dens",
    "velx",
    "vely",
    "velz",
    "press",
    "bcc1",
    "bcc2",
    "bcc3",
)
ADM_VARIABLES = (
    "adm_gxx",
    "adm_gxy",
    "adm_gxz",
    "adm_gyy",
    "adm_gyz",
    "adm_gzz",
    "adm_alpha",
    "adm_betax",
    "adm_betay",
    "adm_betaz",
)
PROFILE_PARAMETER_ORDER = (
    "rho0",
    "pgas0",
    "u1",
    "u2",
    "u3",
    "b1",
    "b2",
    "b3",
    *(f"dlnrho_dxh{i}" for i in range(1, 4)),
    *(f"dlnpgas_dxh{i}" for i in range(1, 4)),
    *(f"du{i}_dxh{j}" for i in range(1, 4) for j in range(1, 4)),
    *(f"db{i}_dxh{j}" for i in range(1, 4) for j in range(1, 4)),
)


@dataclass(frozen=True)
class FitResult:
    intercept: np.ndarray
    gradient: np.ndarray
    weighted_rms: float
    relative_rms: float


@dataclass(frozen=True)
class SampleCloud:
    global_position: np.ndarray
    local_position: np.ndarray
    cell_volume: np.ndarray
    primitive: dict[str, np.ndarray]
    adm: dict[str, np.ndarray]


def _vector(values: Iterable[float], name: str) -> np.ndarray:
    result = np.asarray(tuple(values), dtype=float)
    if result.shape != (3,) or not np.isfinite(result).all():
        raise ValueError(f"{name} must contain three finite values")
    return result


def _normalize(vector: np.ndarray, name: str) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if not norm > 0.0 or not math.isfinite(norm):
        raise ValueError(f"cannot normalize {name}")
    return vector / norm


def canonical_spatial_basis(
    anchor: np.ndarray, primary_center: np.ndarray, disk_normal: np.ndarray
) -> np.ndarray:
    """Return global Cartesian columns (radial, prograde tangent, vertical)."""

    vertical = _normalize(disk_normal, "disk normal")
    separation = anchor - primary_center
    radial = separation - np.dot(separation, vertical) * vertical
    radial = _normalize(radial, "projected source-primary separation")
    tangent = _normalize(np.cross(vertical, radial), "prograde tangent")
    basis = np.column_stack((radial, tangent, vertical))
    if np.linalg.det(basis) <= 0.0:
        raise RuntimeError("canonical local spatial basis is not right-handed")
    return basis


def spacetime_metric_from_adm(
    gamma: np.ndarray, alpha: float, beta: np.ndarray
) -> np.ndarray:
    metric = np.zeros((4, 4), dtype=float)
    metric[1:, 1:] = gamma
    beta_lower = gamma @ beta
    metric[0, 1:] = beta_lower
    metric[1:, 0] = beta_lower
    metric[0, 0] = -alpha * alpha + float(beta @ beta_lower)
    return metric


def metric_inner(metric: np.ndarray, left: np.ndarray, right: np.ndarray) -> float:
    return float(left @ metric @ right)


def build_source_tetrad(
    metric: np.ndarray, source_velocity: np.ndarray, spatial_basis: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Construct contravariant tetrad legs and their dual coframe at the anchor."""

    tangent = np.concatenate(([1.0], source_velocity))
    norm2 = metric_inner(metric, tangent, tangent)
    if not norm2 < 0.0:
        raise ValueError("requested source coordinate velocity is not timelike")
    tetrad = np.zeros((4, 4), dtype=float)
    tetrad[0] = tangent / math.sqrt(-norm2)
    for leg in range(3):
        candidate = np.concatenate(([0.0], spatial_basis[:, leg])).astype(float)
        candidate += metric_inner(metric, candidate, tetrad[0]) * tetrad[0]
        for previous in range(leg):
            candidate -= (
                metric_inner(metric, candidate, tetrad[previous + 1])
                * tetrad[previous + 1]
            )
        candidate_norm2 = metric_inner(metric, candidate, candidate)
        if not candidate_norm2 > 0.0:
            raise RuntimeError("failed to construct source spatial tetrad")
        tetrad[leg + 1] = candidate / math.sqrt(candidate_norm2)

    coframe = np.empty((4, 4), dtype=float)
    for leg, signature in enumerate((-1.0, 1.0, 1.0, 1.0)):
        coframe[leg] = signature * metric @ tetrad[leg]
    gram = tetrad @ metric @ tetrad.T
    expected = np.diag((-1.0, 1.0, 1.0, 1.0))
    if not np.allclose(gram, expected, rtol=0.0, atol=2.0e-11):
        raise RuntimeError("source tetrad failed its orthonormality check")
    if not np.allclose(coframe @ tetrad.T, np.eye(4), rtol=0.0, atol=2.0e-11):
        raise RuntimeError("source coframe is not dual to the tetrad")
    return tetrad, coframe


def _block_key(row: np.ndarray) -> tuple[int, int, int, int]:
    return tuple(int(value) for value in row)


def _aligned_adm_blocks(state: dict, adm: dict) -> list[int]:
    state_logical = np.asarray(state["mb_logical"])
    adm_logical = np.asarray(adm["mb_logical"])
    adm_index = {_block_key(row): index for index, row in enumerate(adm_logical)}
    if len(adm_index) != len(adm_logical):
        raise RuntimeError("ADM dump contains duplicate logical MeshBlocks")
    aligned = []
    for state_index, row in enumerate(state_logical):
        key = _block_key(row)
        if key not in adm_index:
            raise RuntimeError(f"ADM dump is missing state MeshBlock {key}")
        candidate = adm_index[key]
        if not np.allclose(
            np.asarray(state["mb_geometry"])[state_index],
            np.asarray(adm["mb_geometry"])[candidate],
            rtol=0.0,
            atol=1.0e-12,
        ):
            raise RuntimeError(f"state/ADM geometry differs for MeshBlock {key}")
        aligned.append(candidate)
    if len(aligned) != len(adm_logical):
        raise RuntimeError("ADM dump has MeshBlocks absent from the state dump")
    return aligned


def _check_dump(dump: dict, variables: Iterable[str], label: str) -> None:
    missing = sorted(set(variables).difference(dump["mb_data"]))
    if missing:
        raise RuntimeError(f"{label} dump is missing fields: {', '.join(missing)}")
    if int(dump["Nx1"]) < 2 or int(dump["Nx2"]) < 2 or int(dump["Nx3"]) < 2:
        raise RuntimeError("static worldtube extraction requires a three-dimensional dump")


def _block_cell_centers(
    geometry: np.ndarray, shape: tuple[int, int, int]
) -> tuple[np.ndarray, float]:
    nz, ny, nx = shape
    x_faces = np.linspace(float(geometry[0]), float(geometry[1]), nx + 1)
    y_faces = np.linspace(float(geometry[2]), float(geometry[3]), ny + 1)
    z_faces = np.linspace(float(geometry[4]), float(geometry[5]), nz + 1)
    x = 0.5 * (x_faces[:-1] + x_faces[1:])
    y = 0.5 * (y_faces[:-1] + y_faces[1:])
    z = 0.5 * (z_faces[:-1] + z_faces[1:])
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    positions = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    volume = (
        (x_faces[1] - x_faces[0])
        * (y_faces[1] - y_faces[0])
        * (z_faces[1] - z_faces[0])
    )
    return positions, float(volume)


def _distance_squared_to_box(point: np.ndarray, geometry: np.ndarray) -> float:
    lower = np.asarray((geometry[0], geometry[2], geometry[4]), dtype=float)
    upper = np.asarray((geometry[1], geometry[3], geometry[5]), dtype=float)
    offset = np.maximum(np.maximum(lower - point, 0.0), point - upper)
    return float(offset @ offset)


def fit_anchor_adm(
    adm: dict, anchor: np.ndarray, coordinate_radius: float
) -> tuple[np.ndarray, float, np.ndarray, dict[str, float]]:
    positions = []
    volumes = []
    fields = {name: [] for name in ADM_VARIABLES}
    for block, geometry in enumerate(np.asarray(adm["mb_geometry"])):
        if _distance_squared_to_box(anchor, geometry) > coordinate_radius**2:
            continue
        reference = np.asarray(adm["mb_data"]["adm_alpha"][block])
        block_positions, volume = _block_cell_centers(geometry, reference.shape)
        radius2 = np.sum((block_positions - anchor) ** 2, axis=1)
        selected = radius2 <= coordinate_radius**2
        if not np.any(selected):
            continue
        positions.append(block_positions[selected])
        volumes.append(np.full(np.count_nonzero(selected), volume))
        for name in ADM_VARIABLES:
            fields[name].append(np.asarray(adm["mb_data"][name][block]).ravel()[selected])
    if not positions:
        raise RuntimeError("anchor metric fitting sphere contains no ADM cells")
    position = np.concatenate(positions)
    volume = np.concatenate(volumes)
    displacement = position - anchor
    radius2 = np.sum(displacement**2, axis=1)
    weights = volume * np.exp(-radius2 / coordinate_radius**2)
    intercepts = {}
    errors = {}
    for name in ADM_VARIABLES:
        result = weighted_affine_fit(displacement, np.concatenate(fields[name]), weights)
        intercepts[name] = float(result.intercept)
        errors[name] = result.relative_rms
    gamma = np.asarray(
        (
            (intercepts["adm_gxx"], intercepts["adm_gxy"], intercepts["adm_gxz"]),
            (intercepts["adm_gxy"], intercepts["adm_gyy"], intercepts["adm_gyz"]),
            (intercepts["adm_gxz"], intercepts["adm_gyz"], intercepts["adm_gzz"]),
        )
    )
    alpha = intercepts["adm_alpha"]
    beta = np.asarray(
        (intercepts["adm_betax"], intercepts["adm_betay"], intercepts["adm_betaz"])
    )
    if not np.all(np.linalg.eigvalsh(gamma) > 0.0) or not alpha > 0.0:
        raise RuntimeError("fitted anchor ADM state is not positive definite")
    return gamma, alpha, beta, errors


def collect_local_samples(
    state: dict,
    adm: dict,
    adm_order: list[int],
    anchor: np.ndarray,
    coframe: np.ndarray,
    fit_radius: float,
) -> SampleCloud:
    coordinate_inverse = np.linalg.inv(coframe[1:, 1:])
    coordinate_bound = fit_radius * float(np.linalg.norm(coordinate_inverse, ord=2))
    global_positions = []
    local_positions = []
    volumes = []
    primitive = {name: [] for name in STATE_VARIABLES}
    adm_values = {name: [] for name in ADM_VARIABLES}
    geometries = np.asarray(state["mb_geometry"])
    for state_block, geometry in enumerate(geometries):
        if _distance_squared_to_box(anchor, geometry) > coordinate_bound**2:
            continue
        reference = np.asarray(state["mb_data"]["dens"][state_block])
        block_position, volume = _block_cell_centers(geometry, reference.shape)
        displacement = block_position - anchor
        local = displacement @ coframe[1:, 1:].T
        selected = np.sum(local**2, axis=1) <= fit_radius**2
        if not np.any(selected):
            continue
        global_positions.append(block_position[selected])
        local_positions.append(local[selected])
        volumes.append(np.full(np.count_nonzero(selected), volume))
        adm_block = adm_order[state_block]
        for name in STATE_VARIABLES:
            if np.asarray(state["mb_data"][name][state_block]).shape != reference.shape:
                raise RuntimeError(
                    f"state field {name} has inconsistent shape in MeshBlock "
                    f"{_block_key(state['mb_logical'][state_block])}"
                )
            primitive[name].append(
                np.asarray(state["mb_data"][name][state_block]).ravel()[selected]
            )
        for name in ADM_VARIABLES:
            if np.asarray(adm["mb_data"][name][adm_block]).shape != reference.shape:
                raise RuntimeError(
                    f"state/ADM cell shape differs in MeshBlock "
                    f"{_block_key(state['mb_logical'][state_block])}"
                )
            adm_values[name].append(
                np.asarray(adm["mb_data"][name][adm_block]).ravel()[selected]
            )
    if not global_positions:
        raise RuntimeError("local fitting sphere contains no state cells")
    return SampleCloud(
        global_position=np.concatenate(global_positions),
        local_position=np.concatenate(local_positions),
        cell_volume=np.concatenate(volumes),
        primitive={name: np.concatenate(values) for name, values in primitive.items()},
        adm={name: np.concatenate(values) for name, values in adm_values.items()},
    )


def weighted_affine_fit(
    position: np.ndarray, values: np.ndarray, weights: np.ndarray
) -> FitResult:
    position = np.asarray(position, dtype=float)
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if position.ndim != 2 or position.shape[1] != 3:
        raise ValueError("affine-fit positions must have shape (N,3)")
    if values.shape[0] != position.shape[0] or weights.shape != (position.shape[0],):
        raise ValueError("affine-fit arrays have incompatible shapes")
    design = np.column_stack((np.ones(position.shape[0]), position))
    sqrt_weight = np.sqrt(weights)
    weighted_design = design * sqrt_weight[:, None]
    weighted_values = values * sqrt_weight[(slice(None),) + (None,) * (values.ndim - 1)]
    coefficients, _, rank, singular = np.linalg.lstsq(
        weighted_design, weighted_values, rcond=None
    )
    if rank < 4 or singular[-1] <= 0.0:
        raise RuntimeError("affine fit is rank deficient")
    prediction = design @ coefficients
    residual2 = np.sum(weights.reshape((-1,) + (1,) * (values.ndim - 1))
                       * (prediction - values) ** 2)
    normalization = float(np.sum(weights)) * max(1, int(np.prod(values.shape[1:])))
    weighted_rms = math.sqrt(float(residual2) / normalization)
    signal2 = np.sum(weights.reshape((-1,) + (1,) * (values.ndim - 1)) * values**2)
    signal_rms = math.sqrt(float(signal2) / normalization)
    relative = weighted_rms / max(signal_rms, np.finfo(float).tiny)
    return FitResult(
        intercept=np.asarray(coefficients[0]),
        gradient=np.asarray(coefficients[1:]).T,
        weighted_rms=weighted_rms,
        relative_rms=relative,
    )


def fit_divergence_free_field(
    position: np.ndarray, field: np.ndarray, weights: np.ndarray
) -> FitResult:
    """Jointly fit B=B0+G x using eleven parameters with tr(G)=0."""

    position = np.asarray(position, dtype=float)
    field = np.asarray(field, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if field.shape != position.shape or position.ndim != 2 or position.shape[1] != 3:
        raise ValueError("magnetic fit requires position and field arrays of shape (N,3)")
    rows = np.zeros((3 * position.shape[0], 11), dtype=float)
    target = field.reshape(-1)
    for sample, (x, y, z) in enumerate(position):
        row = 3 * sample
        rows[row, 0] = 1.0
        rows[row, 3:6] = (x, y, z)
        rows[row + 1, 1] = 1.0
        rows[row + 1, 6:9] = (x, y, z)
        rows[row + 2, 2] = 1.0
        rows[row + 2, 9:11] = (x, y)
        rows[row + 2, 3] = -z
        rows[row + 2, 7] = -z
    repeated_sqrt_weight = np.repeat(np.sqrt(weights), 3)
    coefficients, _, rank, singular = np.linalg.lstsq(
        rows * repeated_sqrt_weight[:, None],
        target * repeated_sqrt_weight,
        rcond=None,
    )
    if rank < 11 or singular[-1] <= 0.0:
        raise RuntimeError("divergence-free magnetic fit is rank deficient")
    intercept = coefficients[:3]
    gradient = np.asarray(
        (
            coefficients[3:6],
            coefficients[6:9],
            (coefficients[9], coefficients[10], -coefficients[3] - coefficients[7]),
        )
    )
    prediction = intercept + position @ gradient.T
    residual2 = np.sum(weights[:, None] * (prediction - field) ** 2)
    normalization = 3.0 * float(np.sum(weights))
    weighted_rms = math.sqrt(float(residual2) / normalization)
    signal_rms = math.sqrt(float(np.sum(weights[:, None] * field**2)) / normalization)
    return FitResult(
        intercept=intercept,
        gradient=gradient,
        weighted_rms=weighted_rms,
        relative_rms=weighted_rms / max(signal_rms, np.finfo(float).tiny),
    )


def project_grmhd_to_source(
    cloud: SampleCloud, coframe: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Project AthenaK Wv^i and densitized B^i into the fixed source tetrad."""

    count = cloud.local_position.shape[0]
    gamma = np.empty((count, 3, 3), dtype=float)
    gamma[:, 0, 0] = cloud.adm["adm_gxx"]
    gamma[:, 0, 1] = gamma[:, 1, 0] = cloud.adm["adm_gxy"]
    gamma[:, 0, 2] = gamma[:, 2, 0] = cloud.adm["adm_gxz"]
    gamma[:, 1, 1] = cloud.adm["adm_gyy"]
    gamma[:, 1, 2] = gamma[:, 2, 1] = cloud.adm["adm_gyz"]
    gamma[:, 2, 2] = cloud.adm["adm_gzz"]
    alpha = cloud.adm["adm_alpha"]
    beta = np.column_stack(
        (cloud.adm["adm_betax"], cloud.adm["adm_betay"], cloud.adm["adm_betaz"])
    )
    primitive_velocity = np.column_stack(
        (
            cloud.primitive["velx"],
            cloud.primitive["vely"],
            cloud.primitive["velz"],
        )
    )
    lorentz = np.sqrt(
        1.0
        + np.einsum("ni,nij,nj->n", primitive_velocity, gamma, primitive_velocity)
    )
    four_velocity = np.empty((count, 4), dtype=float)
    four_velocity[:, 0] = lorentz / alpha
    four_velocity[:, 1:] = primitive_velocity - four_velocity[:, :1] * beta

    determinant = np.linalg.det(gamma)
    if np.any(determinant <= 0.0) or np.any(alpha <= 0.0):
        raise RuntimeError("sample cloud contains invalid ADM geometry")
    densitized_field = np.column_stack(
        (cloud.primitive["bcc1"], cloud.primitive["bcc2"], cloud.primitive["bcc3"])
    )
    field = densitized_field / np.sqrt(determinant)[:, None]
    velocity_lower = np.einsum("nij,nj->ni", gamma, primitive_velocity)
    field_velocity = np.einsum("ni,ni->n", field, velocity_lower)
    magnetic_four = np.empty((count, 4), dtype=float)
    magnetic_four[:, 0] = field_velocity / alpha
    magnetic_four[:, 1:] = (
        field + (alpha * magnetic_four[:, 0])[:, None] * four_velocity[:, 1:]
    ) / lorentz[:, None]

    source_velocity = four_velocity @ coframe.T
    source_magnetic_four = magnetic_four @ coframe.T
    source_field = (
        source_velocity[:, :1] * source_magnetic_four[:, 1:]
        - source_magnetic_four[:, :1] * source_velocity[:, 1:]
    )
    if not np.isfinite(source_velocity).all() or not np.isfinite(source_field).all():
        raise RuntimeError("source-tetrad GRMHD projection produced non-finite data")
    return source_velocity[:, 1:], source_field


def fit_static_profile(
    cloud: SampleCloud, coframe: np.ndarray, fit_radius: float
) -> tuple[dict[str, float], dict[str, object]]:
    density = cloud.primitive["dens"]
    pressure = cloud.primitive["press"]
    valid = (
        np.isfinite(density)
        & np.isfinite(pressure)
        & (density > 0.0)
        & (pressure > 0.0)
    )
    if np.count_nonzero(valid) < 16:
        raise RuntimeError("fewer than sixteen positive finite GRMHD samples remain")
    position = cloud.local_position[valid]
    radius2 = np.sum(position**2, axis=1)
    weights = cloud.cell_volume[valid] * np.exp(-radius2 / fit_radius**2)
    filtered = SampleCloud(
        global_position=cloud.global_position[valid],
        local_position=position,
        cell_volume=cloud.cell_volume[valid],
        primitive={name: values[valid] for name, values in cloud.primitive.items()},
        adm={name: values[valid] for name, values in cloud.adm.items()},
    )
    source_velocity, source_field = project_grmhd_to_source(filtered, coframe)
    density_fit = weighted_affine_fit(position, np.log(density[valid]), weights)
    pressure_fit = weighted_affine_fit(position, np.log(pressure[valid]), weights)
    velocity_fit = weighted_affine_fit(position, source_velocity, weights)
    magnetic_fit = fit_divergence_free_field(position, source_field, weights)
    parameters: dict[str, float] = {
        "rho0": math.exp(float(density_fit.intercept)),
        "pgas0": math.exp(float(pressure_fit.intercept)),
    }
    for direction in range(3):
        parameters[f"dlnrho_dxh{direction + 1}"] = float(
            density_fit.gradient[direction]
        )
        parameters[f"dlnpgas_dxh{direction + 1}"] = float(
            pressure_fit.gradient[direction]
        )
    for component in range(3):
        parameters[f"u{component + 1}"] = float(velocity_fit.intercept[component])
        parameters[f"b{component + 1}"] = float(magnetic_fit.intercept[component])
        for direction in range(3):
            parameters[f"du{component + 1}_dxh{direction + 1}"] = float(
                velocity_fit.gradient[component, direction]
            )
            parameters[f"db{component + 1}_dxh{direction + 1}"] = float(
                magnetic_fit.gradient[component, direction]
            )
    diagnostics = {
        "sample_count": int(position.shape[0]),
        "density_log_weighted_rms": density_fit.weighted_rms,
        "density_log_relative_rms": density_fit.relative_rms,
        "pressure_log_weighted_rms": pressure_fit.weighted_rms,
        "pressure_log_relative_rms": pressure_fit.relative_rms,
        "velocity_weighted_rms": velocity_fit.weighted_rms,
        "velocity_relative_rms": velocity_fit.relative_rms,
        "magnetic_weighted_rms": magnetic_fit.weighted_rms,
        "magnetic_relative_rms": magnetic_fit.relative_rms,
        "magnetic_gradient_trace": float(np.trace(magnetic_fit.gradient)),
        "minimum_cell_volume": float(np.min(filtered.cell_volume)),
        "maximum_cell_volume": float(np.max(filtered.cell_volume)),
    }
    return parameters, diagnostics


def rescale_profile_parameters(
    parameters: dict[str, float],
    global_length_in_local_units: float,
    density_renormalization: float,
) -> dict[str, float]:
    """Convert a profile from global code units to local code units.

    If one global length unit contains ``L`` local length units, then
    ``x_local=L*x_global`` and ``t_local=L*t_global``.  Density and pressure have
    dimensions length^-2 and magnetic field length^-1.  The optional density
    renormalization preserves pressure/rho and magnetization by scaling B with
    its square root.
    """

    length_scale = float(global_length_in_local_units)
    density_scale = float(density_renormalization)
    if not length_scale > 0.0 or not math.isfinite(length_scale):
        raise ValueError("global length in local units must be finite and positive")
    if not density_scale > 0.0 or not math.isfinite(density_scale):
        raise ValueError("density renormalization must be finite and positive")
    result = dict(parameters)
    thermodynamic_factor = density_scale / length_scale**2
    magnetic_factor = math.sqrt(density_scale) / length_scale
    result["rho0"] *= thermodynamic_factor
    result["pgas0"] *= thermodynamic_factor
    for component in range(1, 4):
        result[f"b{component}"] *= magnetic_factor
    for direction in range(1, 4):
        result[f"dlnrho_dxh{direction}"] /= length_scale
        result[f"dlnpgas_dxh{direction}"] /= length_scale
        for component in range(1, 4):
            result[f"du{component}_dxh{direction}"] /= length_scale
            result[f"db{component}_dxh{direction}"] *= (
                magnetic_factor / length_scale
            )
    return result


def rescale_profile_diagnostics(
    diagnostics: dict[str, object],
    global_length_in_local_units: float,
    density_renormalization: float,
) -> dict[str, object]:
    result = dict(diagnostics)
    magnetic_factor = (
        math.sqrt(density_renormalization) / global_length_in_local_units
    )
    result["magnetic_weighted_rms"] = (
        float(result["magnetic_weighted_rms"]) * magnetic_factor
    )
    result["profile_amplitude_units"] = "local"
    result["cell_volume_units"] = "global_length_cubed"
    return result


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _format_real(value: float) -> str:
    return format(float(value), ".17g")


def write_outputs(
    prefix: Path,
    state_path: Path,
    adm_path: Path,
    state: dict,
    anchor: np.ndarray,
    primary_center: np.ndarray,
    disk_normal: np.ndarray,
    source_velocity: np.ndarray,
    tetrad: np.ndarray,
    coframe: np.ndarray,
    fit_radius: float,
    metric_fit_radius: float,
    global_length_in_local_units: float,
    density_renormalization: float,
    parameters: dict[str, float],
    diagnostics: dict[str, object],
    metric_errors: dict[str, float],
) -> tuple[Path, Path]:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fragment_path = prefix.with_suffix(".athinput")
    manifest_path = prefix.with_suffix(".json")
    lines = [
        "# Generated static first-order EMRI worldtube profile.",
        "# Copy these values into the <problem> block of an emri_windtunnel input.",
        "# global_length_in_local_units = "
        + _format_real(global_length_in_local_units),
        "# density_renormalization = " + _format_real(density_renormalization),
    ]
    lines.extend(
        f"{name:<18} = {_format_real(parameters[name])}"
        for name in PROFILE_PARAMETER_ORDER
    )
    fragment_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    overrides = [
        f"problem/{name}={_format_real(parameters[name])}"
        for name in PROFILE_PARAMETER_ORDER
    ]
    manifest = {
        "schema": 1,
        "classification": "athenak-emri-static-taylor-worldtube",
        "state_file": str(state_path),
        "state_sha256": file_sha256(state_path),
        "adm_file": str(adm_path),
        "adm_sha256": file_sha256(adm_path),
        "time": float(state["time"]),
        "time_global_units": float(state["time"]),
        "time_local_units": float(state["time"])
        * global_length_in_local_units,
        "cycle": int(state["cycle"]),
        "anchor_global": anchor.tolist(),
        "primary_center_global": primary_center.tolist(),
        "disk_normal_global": disk_normal.tolist(),
        "source_coordinate_velocity": source_velocity.tolist(),
        "fit_radius_source": fit_radius,
        "metric_fit_radius_coordinate": metric_fit_radius,
        "global_length_in_local_units": global_length_in_local_units,
        "density_renormalization": density_renormalization,
        "source_tetrad_contravariant": tetrad.tolist(),
        "source_tetrad_coframe": coframe.tolist(),
        "parameters": parameters,
        "athena_overrides": overrides,
        "fit_diagnostics": diagnostics,
        "anchor_metric_relative_fit_errors": metric_errors,
        "limitations": [
            "one-way static first-order fit; no response of the global disk",
            "fixed anchor tetrad across the fitting sphere",
            "cell-centered magnetic data are replaced by a trace-free linear fit",
        ],
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return fragment_path, manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", type=Path, required=True, help="mhd_w_bcc binary dump")
    parser.add_argument("--adm", type=Path, required=True, help="co-temporal adm binary dump")
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--anchor", type=float, nargs=3, required=True)
    parser.add_argument("--primary-center", type=float, nargs=3, default=(0.0, 0.0, 0.0))
    parser.add_argument("--disk-normal", type=float, nargs=3, default=(0.0, 0.0, 1.0))
    velocity_group = parser.add_mutually_exclusive_group(required=True)
    velocity_group.add_argument(
        "--source-velocity", type=float, nargs=3, help="global coordinate dx^i/dt"
    )
    velocity_group.add_argument(
        "--orbital-omega",
        type=float,
        help="derive dx/dt=Omega disk_normal cross (anchor-primary_center)",
    )
    parser.add_argument("--fit-radius", type=float, required=True)
    parser.add_argument("--metric-fit-radius", type=float)
    parser.add_argument("--global-length-in-local-units", type=float, default=1.0)
    parser.add_argument("--density-renormalization", type=float)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    state_path = args.state.expanduser().resolve(strict=True)
    adm_path = args.adm.expanduser().resolve(strict=True)
    anchor = _vector(args.anchor, "anchor")
    primary_center = _vector(args.primary_center, "primary center")
    disk_normal = _normalize(_vector(args.disk_normal, "disk normal"), "disk normal")
    if not args.fit_radius > 0.0 or not math.isfinite(args.fit_radius):
        raise SystemExit("--fit-radius must be finite and positive")
    metric_fit_radius = args.metric_fit_radius or args.fit_radius
    if not metric_fit_radius > 0.0 or not math.isfinite(metric_fit_radius):
        raise SystemExit("--metric-fit-radius must be finite and positive")
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
    if args.source_velocity is not None:
        source_velocity = _vector(args.source_velocity, "source velocity")
    else:
        source_velocity = args.orbital_omega * np.cross(
            disk_normal, anchor - primary_center
        )
    spatial_basis = canonical_spatial_basis(anchor, primary_center, disk_normal)

    state = bin_convert.read_binary(str(state_path))
    adm = bin_convert.read_binary(str(adm_path))
    _check_dump(state, STATE_VARIABLES, "state")
    _check_dump(adm, ADM_VARIABLES, "ADM")
    state_time = float(state["time"])
    adm_time = float(adm["time"])
    time_tolerance = 64.0 * np.finfo(float).eps * max(
        1.0, abs(state_time), abs(adm_time)
    )
    if int(state["cycle"]) != int(adm["cycle"]) or not math.isclose(
        state_time, adm_time, rel_tol=0.0, abs_tol=time_tolerance
    ):
        raise RuntimeError("state and ADM dumps are not co-temporal")
    adm_order = _aligned_adm_blocks(state, adm)
    gamma, alpha, beta, metric_errors = fit_anchor_adm(
        adm, anchor, metric_fit_radius
    )
    metric = spacetime_metric_from_adm(gamma, alpha, beta)
    tetrad, coframe = build_source_tetrad(metric, source_velocity, spatial_basis)
    cloud = collect_local_samples(
        state, adm, adm_order, anchor, coframe, args.fit_radius
    )
    parameters, diagnostics = fit_static_profile(cloud, coframe, args.fit_radius)
    parameters = rescale_profile_parameters(
        parameters,
        args.global_length_in_local_units,
        density_renormalization,
    )
    diagnostics = rescale_profile_diagnostics(
        diagnostics,
        args.global_length_in_local_units,
        density_renormalization,
    )
    fragment, manifest = write_outputs(
        args.output_prefix.expanduser().resolve(),
        state_path,
        adm_path,
        state,
        anchor,
        primary_center,
        disk_normal,
        source_velocity,
        tetrad,
        coframe,
        args.fit_radius,
        metric_fit_radius,
        args.global_length_in_local_units,
        density_renormalization,
        parameters,
        diagnostics,
        metric_errors,
    )
    print(fragment)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
