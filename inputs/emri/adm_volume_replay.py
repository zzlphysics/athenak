#!/usr/bin/env python3
"""Extract and validate a time-dependent local ADM volume for inner replay.

The source numerical-ADM metric is pulled back through the same affine frame used
by the fluid worldtube.  An analytic secondary Kerr-Schild perturbation is then
embedded in the orthonormal tangent frame of that numerical background.  The
resulting effective four-metric can be stored directly, or the smooth primary
four-metric and its spatial derivatives can be stored separately for hybrid
runtime composition with an analytic secondary on the fluid/AMR grid.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import struct
import zlib

import numpy as np

import build_adm_ephemeris_frame as adm_frame
import extract_global_worldtube as extract
import extract_static_taylor_worldtube as static
import worldtube_flux_emf as worldtube


CLASSIFICATION = "athenak-emri-adm-volume-replay-v1"
LEGACY_BINARY_CLASSIFICATION = "athenak-emri-adm-volume-binary-v1"
BINARY_CLASSIFICATION = "athenak-emri-adm-volume-binary-v2"
HYBRID_BINARY_CLASSIFICATION = "athenak-emri-adm-volume-binary-v3"
BINARY_MAGIC = b"AEMRIADMVOL001\x00\x00"
LEGACY_BINARY_VERSION = 1
BINARY_VERSION = 2
HYBRID_BINARY_VERSION = 3
BINARY_HEADER = struct.Struct("<16sIIIIII6dQ")
FIELD_NAMES = (
    "g00",
    "g01",
    "g02",
    "g03",
    "g11",
    "g12",
    "g13",
    "g22",
    "g23",
    "g33",
    "K11",
    "K12",
    "K13",
    "K22",
    "K23",
    "K33",
)
METRIC_COMPONENTS = (
    (0, 0),
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 1),
    (1, 2),
    (1, 3),
    (2, 2),
    (2, 3),
    (3, 3),
)
CURVATURE_COMPONENTS = (
    (0, 0),
    (0, 1),
    (0, 2),
    (1, 1),
    (1, 2),
    (2, 2),
)
HYBRID_FIELD_NAMES = tuple(
    list(FIELD_NAMES[: len(METRIC_COMPONENTS)])
    + [
        f"d{direction}_{FIELD_NAMES[field]}"
        for direction in ("x1", "x2", "x3")
        for field in range(len(METRIC_COMPONENTS))
    ]
)
HYBRID_PARAMETER_NAMES = (
    "secondary_mass",
    "secondary_chi",
    "secondary_spin_buffer",
    "secondary_singularity_floor",
    "secondary_metric_fd_step",
)


@dataclass(frozen=True)
class ADMVolume:
    times: np.ndarray
    lower: np.ndarray
    spacing: np.ndarray
    fields: np.ndarray
    metadata: dict[str, object]
    secondary_coframes: np.ndarray | None = None
    hybrid_parameters: np.ndarray | None = None


def _finite_positive(value: object, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def _global_metric(values: object) -> np.ndarray:
    source = np.asarray(values, dtype=np.float64)
    if source.shape != (len(extract.ADM_VARIABLES),):
        raise ValueError("ADM sample has the wrong field count")
    return adm_frame.adm_metric(source)


class LocalMetricSampler:
    def __init__(
        self,
        snapshots: extract.SourceSeries,
        frames: extract.AffineFrameSeries,
        length_scale: float,
    ) -> None:
        self.snapshots = snapshots
        self.frames = frames
        self.length_scale = _finite_positive(length_scale, "length scale")
        self.sampled_event_count = 0
        self.minimum_lapse = math.inf
        self.minimum_spatial_eigenvalue = math.inf

    def sample(self, local_time: float, local_positions: object) -> np.ndarray:
        positions = np.asarray(local_positions, dtype=np.float64)
        if (
            positions.ndim != 2
            or positions.shape[1] != 3
            or not np.isfinite(positions).all()
        ):
            raise ValueError("local metric positions must have finite shape (N,3)")
        worldline, instantaneous = self.frames.evaluate(float(local_time))
        events = worldline[None, :] + positions @ instantaneous.spatial_legs.T
        source = self.snapshots.sample(events)
        metrics = np.empty((positions.shape[0], 4, 4), dtype=np.float64)
        for index, position in enumerate(positions):
            jacobian = instantaneous.jacobian(position)
            global_metric = _global_metric(source[index])
            metrics[index] = self.length_scale**2 * (
                jacobian.T @ global_metric @ jacobian
            )
        _, audit = decompose_four_metric(metrics)
        self.sampled_event_count += positions.shape[0]
        self.minimum_lapse = min(
            self.minimum_lapse, float(np.min(audit["lapse"]))
        )
        self.minimum_spatial_eigenvalue = min(
            self.minimum_spatial_eigenvalue,
            float(np.min(audit["minimum_spatial_eigenvalue"])),
        )
        return metrics

    def diagnostics(self) -> dict[str, object]:
        return {
            "sampled_spacetime_event_count": self.sampled_event_count,
            "minimum_background_lapse": self.minimum_lapse,
            "minimum_background_spatial_metric_eigenvalue": (
                self.minimum_spatial_eigenvalue
            ),
            "snapshot_loading": self.snapshots.loading_document(),
        }


def decompose_four_metric(
    metric: object,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    values = np.asarray(metric, dtype=np.float64)
    if values.shape[-2:] != (4, 4) or not np.isfinite(values).all():
        raise ValueError("four-metric array must be finite with trailing shape (4,4)")
    symmetry = np.max(np.abs(values - np.swapaxes(values, -1, -2)))
    if symmetry > 1.0e-11:
        raise ValueError("four-metric array is not symmetric")
    gamma = values[..., 1:, 1:]
    beta_lower = values[..., 0, 1:]
    eigenvalues = np.linalg.eigvalsh(gamma)
    if float(np.min(eigenvalues)) <= 0.0:
        raise ValueError("spatial metric is not positive definite")
    beta = np.linalg.solve(gamma, beta_lower[..., None])[..., 0]
    lapse_squared = np.einsum("...i,...i->...", beta_lower, beta) - values[
        ..., 0, 0
    ]
    if float(np.min(lapse_squared)) <= 0.0:
        raise ValueError("four-metric has non-positive lapse squared")
    return {
        "alpha": np.sqrt(lapse_squared),
        "beta": beta,
        "gamma": gamma,
    }, {
        "lapse": np.sqrt(lapse_squared),
        "minimum_spatial_eigenvalue": eigenvalues[..., 0],
    }


def _source_tetrad(background_metric: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return static.build_source_tetrad(
        background_metric, np.zeros(3), np.eye(3)
    )


def secondary_kerr_perturbation(
    positions: object,
    mass: float,
    chi: float,
    coframe: object,
    spin_buffer: float = 0.05,
    singularity_floor: float = 1.0e-3,
) -> np.ndarray:
    """Return the regularized secondary Kerr-Schild term in local coordinates."""

    coordinates = np.asarray(positions, dtype=np.float64)
    basis = np.asarray(coframe, dtype=np.float64)
    mass = _finite_positive(mass, "secondary mass")
    spin_buffer = _finite_positive(spin_buffer, "secondary spin buffer")
    singularity_floor = _finite_positive(
        singularity_floor, "secondary singularity floor"
    )
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("secondary positions must have shape (N,3)")
    if basis.shape != (4, 4) or not np.isfinite(basis).all():
        raise ValueError("secondary coframe must have finite shape (4,4)")
    if not math.isfinite(chi) or abs(chi) > 1.0:
        raise ValueError("secondary dimensionless spin must lie in [-1,1]")
    displacement = np.column_stack((np.zeros(coordinates.shape[0]), coordinates))
    rest = displacement @ basis[1:].T
    spin = np.asarray((0.0, 0.0, chi * mass))
    spin_squared = float(spin @ spin)
    radius_squared = np.einsum("ni,ni->n", rest, rest)
    spin_dot_position = rest @ spin
    radial_term = radius_squared - spin_squared
    kerr_radius_squared = 0.5 * (
        radial_term
        + np.sqrt(radial_term**2 + 4.0 * spin_dot_position**2)
    )
    horizon_radius = mass + math.sqrt(max(mass * mass - spin_squared, 0.0))
    ring_radius = math.sqrt(spin_squared) * spin_buffer + mass * singularity_floor
    regularization_radius = max(horizon_radius, ring_radius)
    outer_squared = regularization_radius**2
    inner_squared = 0.25 * outer_squared
    active = kerr_radius_squared > inner_squared
    taper = np.zeros(coordinates.shape[0])
    taper[kerr_radius_squared >= outer_squared] = 1.0
    transition = active & (kerr_radius_squared < outer_squared)
    u = (kerr_radius_squared[transition] - inner_squared) / (
        outer_squared - inner_squared
    )
    taper[transition] = u**3 * (10.0 - 15.0 * u + 6.0 * u**2)

    result_tetrad = np.zeros((coordinates.shape[0], 4, 4))
    safe = active & (kerr_radius_squared > 0.0)
    if np.any(safe):
        r2 = kerr_radius_squared[safe]
        radius = np.sqrt(r2)
        adotx = spin_dot_position[safe]
        denominator = r2**2 + adotx**2
        hfun = taper[safe] * mass * r2 * radius / denominator
        cross = np.cross(rest[safe], spin)
        spatial_null = (
            radius[:, None] * rest[safe]
            + cross
            + (adotx / radius)[:, None] * spin[None, :]
        ) / (r2 + spin_squared)[:, None]
        null = np.column_stack((np.ones(radius.size), spatial_null))
        result_tetrad[safe] = (
            2.0 * hfun[:, None, None] * null[:, :, None] * null[:, None, :]
        )
    return np.einsum("am,nab,bj->nmj", basis, result_tetrad, basis)


def extrinsic_curvature(
    metric: object, times: object, spacing: object
) -> tuple[np.ndarray, dict[str, object]]:
    """Compute K_ij from the sampled four-metric and ADM evolution identity."""

    values = np.asarray(metric, dtype=np.float64)
    sample_times = worldtube.validate_times(times)
    delta = np.asarray(spacing, dtype=np.float64)
    if values.ndim != 6 or values.shape[0] != sample_times.size:
        raise ValueError("metric series must have shape (nt,nz,ny,nx,4,4)")
    if values.shape[-2:] != (4, 4) or min(values.shape[1:4]) < 3:
        raise ValueError("metric volume needs at least three nodes on every axis")
    if delta.shape != (3,) or not np.isfinite(delta).all() or np.min(delta) <= 0.0:
        raise ValueError("metric volume spacing must contain three positive values")
    adm, audit = decompose_four_metric(values)
    gamma = adm["gamma"]
    beta = adm["beta"]
    beta_lower = values[..., 0, 1:]
    time_order = 2 if sample_times.size >= 3 else 1
    gamma_time = np.gradient(
        gamma, sample_times, axis=0, edge_order=time_order
    )
    curvature = np.empty_like(gamma)
    spatial_axes = (2, 1, 0)
    maximum_identity_residual = 0.0
    for time_index in range(sample_times.size):
        spatial_gamma = []
        spatial_beta_lower = []
        for direction, array_axis in enumerate(spatial_axes):
            spatial_gamma.append(
                np.gradient(
                    gamma[time_index],
                    delta[direction],
                    axis=array_axis,
                    edge_order=2,
                )
            )
            spatial_beta_lower.append(
                np.gradient(
                    beta_lower[time_index],
                    delta[direction],
                    axis=array_axis,
                    edge_order=2,
                )
            )
        for i in range(3):
            for j in range(i, 3):
                covariant_sum = (
                    spatial_beta_lower[i][..., j]
                    + spatial_beta_lower[j][..., i]
                )
                connection_contraction = np.zeros_like(covariant_sum)
                for lower in range(3):
                    connection_contraction += beta[time_index, ..., lower] * (
                        spatial_gamma[i][..., lower, j]
                        + spatial_gamma[j][..., lower, i]
                        - spatial_gamma[lower][..., i, j]
                    )
                numerator = (
                    covariant_sum
                    - connection_contraction
                    - gamma_time[time_index, ..., i, j]
                )
                value = numerator / (2.0 * adm["alpha"][time_index])
                curvature[time_index, ..., i, j] = value
                curvature[time_index, ..., j, i] = value
                reconstructed = (
                    covariant_sum
                    - connection_contraction
                    - 2.0 * adm["alpha"][time_index] * value
                )
                maximum_identity_residual = max(
                    maximum_identity_residual,
                    float(
                        np.max(
                            np.abs(
                                reconstructed
                                - gamma_time[time_index, ..., i, j]
                            )
                        )
                    ),
                )
    return curvature, {
        "minimum_lapse": float(np.min(audit["lapse"])),
        "minimum_spatial_metric_eigenvalue": float(
            np.min(audit["minimum_spatial_eigenvalue"])
        ),
        "maximum_adm_evolution_identity_residual": maximum_identity_residual,
        "temporal_derivative": (
            f"nonuniform finite difference with edge_order={time_order}"
        ),
        "spatial_derivative": "second-order finite difference on the volume grid",
    }


def spatial_metric_derivatives(metric: object, spacing: object) -> np.ndarray:
    """Differentiate a metric volume along local x1, x2, and x3."""

    values = np.asarray(metric, dtype=np.float64)
    delta = np.asarray(spacing, dtype=np.float64)
    if values.ndim != 6 or values.shape[-2:] != (4, 4):
        raise ValueError("metric series must have shape (nt,nz,ny,nx,4,4)")
    if min(values.shape[1:4]) < 3:
        raise ValueError("metric volume needs at least three nodes on every axis")
    if delta.shape != (3,) or not np.isfinite(delta).all() or np.min(delta) <= 0.0:
        raise ValueError("metric volume spacing must contain three positive values")
    result = np.empty(values.shape[:-2] + (3, 4, 4), dtype=np.float64)
    # The stored volume is indexed (z,y,x), while derivative directions are
    # ordered (x1,x2,x3).
    for direction, array_axis in enumerate((3, 2, 1)):
        result[..., direction, :, :] = np.gradient(
            values, delta[direction], axis=array_axis, edge_order=2
        )
    return result


def decompose_metric_derivatives(
    metric: object, metric_derivatives: object
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Return ADM variables and K_ij from g_ab and partial_c g_ab.

    ``metric_derivatives[..., c, a, b]`` uses c=0 for time and c=1,2,3
    for x1,x2,x3.  Reconstructing K from the composed metric avoids assuming
    that primary and secondary extrinsic curvatures add independently.
    """

    values = np.asarray(metric, dtype=np.float64)
    derivatives = np.asarray(metric_derivatives, dtype=np.float64)
    if derivatives.shape != values.shape[:-2] + (4, 4, 4):
        raise ValueError("metric derivatives have incompatible dimensions")
    if not np.isfinite(derivatives).all():
        raise ValueError("metric derivatives must be finite")
    adm, _ = decompose_four_metric(values)
    beta = adm["beta"]
    curvature = np.empty(values.shape[:-2] + (3, 3), dtype=np.float64)
    for i in range(3):
        for j in range(i, 3):
            covariant_sum = (
                derivatives[..., i + 1, 0, j + 1]
                + derivatives[..., j + 1, 0, i + 1]
            )
            connection_contraction = np.zeros(values.shape[:-2], dtype=np.float64)
            for lower in range(3):
                connection_contraction += beta[..., lower] * (
                    derivatives[..., i + 1, lower + 1, j + 1]
                    + derivatives[..., j + 1, lower + 1, i + 1]
                    - derivatives[..., lower + 1, i + 1, j + 1]
                )
            value = (
                covariant_sum
                - connection_contraction
                - derivatives[..., 0, i + 1, j + 1]
            ) / (2.0 * adm["alpha"])
            curvature[..., i, j] = value
            curvature[..., j, i] = value
    return adm, curvature


def _source_series(
    manifest_path: Path,
) -> tuple[extract.SourceSeries, dict[str, object], Path]:
    scan = extract.scan_snapshot_manifest(manifest_path)

    def load_snapshot(index: int) -> extract.Snapshot:
        return extract._load_adm_snapshot(
            scan.entries[index], scan.path.parent, scan.source_level
        )

    series = extract.LazySnapshotSeries(
        scan.descriptors, load_snapshot, cache_size=scan.snapshot_cache_size
    )
    return series, scan.document, scan.path


def build_volume(
    manifest_path: Path,
    times: object,
    half_width: float,
    cells: int,
    halo: int,
    secondary_mass: float,
    secondary_chi: float = 0.0,
    chunk_size: int = 65536,
    spin_buffer: float = 0.05,
    singularity_floor: float = 1.0e-3,
    hybrid_primary: bool = False,
    secondary_metric_fd_step: float | None = None,
) -> ADMVolume:
    sample_times = worldtube.validate_times(times)
    half_width = _finite_positive(half_width, "volume half width")
    secondary_mass = _finite_positive(secondary_mass, "secondary mass")
    if secondary_metric_fd_step is None:
        secondary_metric_fd_step = 5.0e-5 * secondary_mass
    secondary_metric_fd_step = _finite_positive(
        secondary_metric_fd_step, "secondary metric finite-difference step"
    )
    if not isinstance(cells, int) or isinstance(cells, bool) or cells < 2:
        raise ValueError("volume cells must be an integer of at least two")
    if not isinstance(halo, int) or isinstance(halo, bool) or halo < 1:
        raise ValueError("volume halo must be a positive integer")
    if not isinstance(chunk_size, int) or chunk_size < 1:
        raise ValueError("chunk size must be positive")
    snapshots, document, resolved_manifest = _source_series(manifest_path)
    frame_document = extract._frame_document(
        document.get("frame"), resolved_manifest.parent
    )
    frames = extract.AffineFrameSeries.from_document(frame_document)
    length_scale = _finite_positive(
        document.get("global_length_in_local_units", 1.0),
        "global_length_in_local_units",
    )
    sampler = LocalMetricSampler(snapshots, frames, length_scale)
    spacing = 2.0 * half_width / cells
    nodes = cells + 2 * halo
    lower = np.full(3, -half_width - (halo - 0.5) * spacing)
    coordinates = lower[0] + spacing * np.arange(nodes)
    z, y, x = np.meshgrid(coordinates, coordinates, coordinates, indexing="ij")
    positions = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    metric = np.empty((sample_times.size, nodes, nodes, nodes, 4, 4))
    secondary_coframes = np.empty((sample_times.size, 4, 4))
    coframe_diagnostics = []
    minimum_composed_lapse = math.inf
    minimum_composed_spatial_eigenvalue = math.inf
    for time_index, local_time in enumerate(sample_times):
        background_origin = sampler.sample(
            float(local_time), np.zeros((1, 3))
        )[0]
        tetrad, coframe = _source_tetrad(background_origin)
        secondary_coframes[time_index] = coframe
        gram = tetrad @ background_origin @ tetrad.T
        dual = coframe @ tetrad.T
        coframe_diagnostics.append(
            {
                "time": float(local_time),
                "maximum_orthonormal_gram_error": float(
                    np.max(np.abs(gram - np.diag((-1.0, 1.0, 1.0, 1.0))))
                ),
                "maximum_dual_identity_error": float(
                    np.max(np.abs(dual - np.eye(4)))
                ),
            }
        )
        flattened = metric[time_index].reshape(-1, 4, 4)
        for start in range(0, positions.shape[0], chunk_size):
            stop = min(start + chunk_size, positions.shape[0])
            background = sampler.sample(float(local_time), positions[start:stop])
            if hybrid_primary:
                flattened[start:stop] = background
                secondary = secondary_kerr_perturbation(
                    positions[start:stop],
                    secondary_mass,
                    secondary_chi,
                    coframe,
                    spin_buffer,
                    singularity_floor,
                )
                _, composed_audit = decompose_four_metric(background + secondary)
                minimum_composed_lapse = min(
                    minimum_composed_lapse,
                    float(np.min(composed_audit["lapse"])),
                )
                minimum_composed_spatial_eigenvalue = min(
                    minimum_composed_spatial_eigenvalue,
                    float(np.min(composed_audit["minimum_spatial_eigenvalue"])),
                )
            else:
                secondary = secondary_kerr_perturbation(
                    positions[start:stop],
                    secondary_mass,
                    secondary_chi,
                    coframe,
                    spin_buffer,
                    singularity_floor,
                )
                flattened[start:stop] = background + secondary
    hybrid_parameters = None
    if hybrid_primary:
        derivatives = spatial_metric_derivatives(metric, np.full(3, spacing))
        fields = np.empty(
            (sample_times.size, len(HYBRID_FIELD_NAMES), nodes, nodes, nodes),
            dtype=np.float64,
        )
        for field, (left, right) in enumerate(METRIC_COMPONENTS):
            fields[:, field] = metric[..., left, right]
            for direction in range(3):
                offset = len(METRIC_COMPONENTS) * (direction + 1) + field
                fields[:, offset] = derivatives[..., direction, left, right]
        _, background_audit = decompose_four_metric(metric)
        adm_diagnostics = {
            "minimum_lapse": minimum_composed_lapse,
            "minimum_spatial_metric_eigenvalue": (
                minimum_composed_spatial_eigenvalue
            ),
            "minimum_primary_lapse": float(np.min(background_audit["lapse"])),
            "minimum_primary_spatial_metric_eigenvalue": float(
                np.min(background_audit["minimum_spatial_eigenvalue"])
            ),
            "temporal_derivative": (
                "exact derivative of piecewise-linear endpoint composition at runtime"
            ),
            "spatial_derivative": (
                "stored second-order primary derivative plus online analytic secondary"
            ),
        }
        hybrid_parameters = np.asarray(
            (
                secondary_mass,
                secondary_chi,
                spin_buffer,
                singularity_floor,
                secondary_metric_fd_step,
            ),
            dtype=np.float64,
        )
        field_names = HYBRID_FIELD_NAMES
        representation = "primary_metric_derivatives_plus_online_secondary"
    else:
        curvature, adm_diagnostics = extrinsic_curvature(
            metric, sample_times, np.full(3, spacing)
        )
        fields = np.empty(
            (sample_times.size, len(FIELD_NAMES), nodes, nodes, nodes),
            dtype=np.float64,
        )
        for field, (left, right) in enumerate(METRIC_COMPONENTS):
            fields[:, field] = metric[..., left, right]
        for offset, (left, right) in enumerate(CURVATURE_COMPONENTS):
            fields[:, len(METRIC_COMPONENTS) + offset] = curvature[..., left, right]
        field_names = FIELD_NAMES
        representation = "composed_metric_and_extrinsic_curvature"
    metadata = {
        "classification": CLASSIFICATION,
        "source_manifest": str(resolved_manifest),
        "source_frame_classification": frame_document.get("classification"),
        "global_length_in_local_units": length_scale,
        "half_width": half_width,
        "active_cells_per_axis": cells,
        "halo_nodes_per_side": halo,
        "secondary_mass": secondary_mass,
        "secondary_chi": secondary_chi,
        "secondary_embedding": "numerical-background tangent tetrad",
        "secondary_spin_buffer": spin_buffer,
        "secondary_singularity_floor": singularity_floor,
        "field_names": list(field_names),
        "representation": representation,
        "sampling": sampler.diagnostics(),
        "coframe_diagnostics": coframe_diagnostics,
        "adm_diagnostics": adm_diagnostics,
    }
    return ADMVolume(
        sample_times,
        lower,
        np.full(3, spacing),
        fields,
        metadata,
        secondary_coframes,
        hybrid_parameters,
    )


def write_binary(path: Path, volume: ADMVolume) -> dict[str, object]:
    times = worldtube.validate_times(volume.times)
    lower = np.asarray(volume.lower, dtype=np.float64)
    spacing = np.asarray(volume.spacing, dtype=np.float64)
    fields = np.asarray(volume.fields, dtype=np.float64)
    coframes = None
    if volume.secondary_coframes is not None:
        coframes = np.asarray(volume.secondary_coframes, dtype=np.float64)
    hybrid_parameters = None
    if volume.hybrid_parameters is not None:
        hybrid_parameters = np.asarray(volume.hybrid_parameters, dtype=np.float64)
    if lower.shape != (3,) or spacing.shape != (3,):
        raise ValueError("ADM volume lower and spacing must contain three values")
    if not np.isfinite(lower).all() or not np.isfinite(spacing).all():
        raise ValueError("ADM volume geometry is non-finite")
    if np.min(spacing) <= 0.0:
        raise ValueError("ADM volume spacing must be positive")
    field_names = HYBRID_FIELD_NAMES if hybrid_parameters is not None else FIELD_NAMES
    if fields.ndim != 5 or fields.shape[:2] != (times.size, len(field_names)):
        raise ValueError("ADM volume fields have incompatible dimensions")
    if min(fields.shape[2:]) < 2 or not np.isfinite(fields).all():
        raise ValueError("ADM volume fields are invalid")
    if coframes is not None and (
        coframes.shape != (times.size, 4, 4) or not np.isfinite(coframes).all()
    ):
        raise ValueError("secondary coframes must have finite shape (nt,4,4)")
    if hybrid_parameters is not None:
        if coframes is None:
            raise ValueError("hybrid ADM volume requires secondary coframes")
        if (
            hybrid_parameters.shape != (len(HYBRID_PARAMETER_NAMES),)
            or not np.isfinite(hybrid_parameters).all()
        ):
            raise ValueError("hybrid ADM parameters have incompatible dimensions")
        if (
            hybrid_parameters[0] <= 0.0
            or abs(hybrid_parameters[1]) > 1.0
            or hybrid_parameters[2] < 0.0
            or hybrid_parameters[3] <= 0.0
            or hybrid_parameters[4] <= 0.0
        ):
            raise ValueError("hybrid ADM parameters are outside their physical ranges")
        version = HYBRID_BINARY_VERSION
        classification = HYBRID_BINARY_CLASSIFICATION
    elif coframes is not None:
        version = BINARY_VERSION
        classification = BINARY_CLASSIFICATION
    else:
        version = LEGACY_BINARY_VERSION
        classification = LEGACY_BINARY_CLASSIFICATION
    nz, ny, nx = fields.shape[2:]
    payload = times.astype("<f8", copy=False).tobytes(order="C")
    if coframes is not None:
        payload += coframes.astype("<f8", copy=False).tobytes(order="C")
    if hybrid_parameters is not None:
        payload += hybrid_parameters.astype("<f8", copy=False).tobytes(order="C")
    payload += fields.astype("<f8", copy=False).tobytes(order="C")
    checksum = zlib.crc32(payload) & 0xFFFFFFFF
    header = BINARY_HEADER.pack(
        BINARY_MAGIC,
        version,
        times.size,
        len(field_names),
        nx,
        ny,
        nz,
        *lower,
        *spacing,
        checksum,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + payload)
    sidecar = dict(volume.metadata)
    sidecar.update(
        {
            "classification": classification,
            "binary_file": path.name,
            "binary_version": version,
            "payload_crc32": f"{checksum:08x}",
            "field_names": list(field_names),
            "times": times.tolist(),
            "grid_lower_node": lower.tolist(),
            "grid_spacing": spacing.tolist(),
            "grid_shape_xyz": [nx, ny, nz],
        }
    )
    if hybrid_parameters is not None:
        sidecar["hybrid_parameter_names"] = list(HYBRID_PARAMETER_NAMES)
        sidecar["hybrid_parameters"] = hybrid_parameters.tolist()
    sidecar_path = path.with_suffix(path.suffix + ".json")
    sidecar_path.write_text(
        json.dumps(sidecar, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return sidecar


def read_binary(path: Path) -> ADMVolume:
    payload = path.read_bytes()
    if len(payload) < BINARY_HEADER.size:
        raise ValueError("ADM volume binary is shorter than its header")
    header = BINARY_HEADER.unpack_from(payload)
    magic, version, nt, nvar, nx, ny, nz, *remaining = header
    lower = np.asarray(remaining[:3])
    spacing = np.asarray(remaining[3:6])
    expected_checksum = remaining[6]
    if magic != BINARY_MAGIC or version not in (
        LEGACY_BINARY_VERSION,
        BINARY_VERSION,
        HYBRID_BINARY_VERSION,
    ):
        raise ValueError("ADM volume binary magic or version is unsupported")
    field_names = HYBRID_FIELD_NAMES if version == HYBRID_BINARY_VERSION else FIELD_NAMES
    if nt < 2 or nvar != len(field_names) or min(nx, ny, nz) < 2:
        raise ValueError("ADM volume binary dimensions are invalid")
    binary_payload = payload[BINARY_HEADER.size :]
    checksum = zlib.crc32(binary_payload) & 0xFFFFFFFF
    if checksum != expected_checksum:
        raise ValueError("ADM volume binary checksum mismatch")
    coframe_values = nt * 16 if version >= BINARY_VERSION else 0
    parameter_values = (
        len(HYBRID_PARAMETER_NAMES) if version == HYBRID_BINARY_VERSION else 0
    )
    expected_values = (
        nt + coframe_values + parameter_values + nt * nvar * nx * ny * nz
    )
    values = np.frombuffer(binary_payload, dtype="<f8")
    if values.size != expected_values or not np.isfinite(values).all():
        raise ValueError("ADM volume binary payload dimensions are invalid")
    times = worldtube.validate_times(values[:nt])
    cursor = nt
    coframes = None
    if version >= BINARY_VERSION:
        coframes = np.asarray(values[cursor : cursor + coframe_values].reshape(nt, 4, 4))
        cursor += coframe_values
    hybrid_parameters = None
    if version == HYBRID_BINARY_VERSION:
        hybrid_parameters = np.asarray(values[cursor : cursor + parameter_values])
        cursor += parameter_values
        if (
            hybrid_parameters[0] <= 0.0
            or abs(hybrid_parameters[1]) > 1.0
            or hybrid_parameters[2] < 0.0
            or hybrid_parameters[3] <= 0.0
            or hybrid_parameters[4] <= 0.0
        ):
            raise ValueError("hybrid ADM binary parameters are invalid")
    fields = np.asarray(values[cursor:].reshape(nt, nvar, nz, ny, nx))
    metadata: dict[str, object] = {
        "classification": (
            HYBRID_BINARY_CLASSIFICATION
            if version == HYBRID_BINARY_VERSION
            else (
                BINARY_CLASSIFICATION
                if version >= BINARY_VERSION
                else LEGACY_BINARY_CLASSIFICATION
            )
        ),
        "binary_version": version,
        "payload_crc32": f"{checksum:08x}",
        "field_names": list(field_names),
    }
    if hybrid_parameters is not None:
        metadata["hybrid_parameter_names"] = list(HYBRID_PARAMETER_NAMES)
        metadata["hybrid_parameters"] = hybrid_parameters.tolist()
    sidecar_path = path.with_suffix(path.suffix + ".json")
    if sidecar_path.exists():
        sidecar = json.loads(sidecar_path.read_text(encoding="utf-8"))
        if (
            sidecar.get("classification") != metadata["classification"]
            or sidecar.get("payload_crc32") != f"{checksum:08x}"
            or sidecar.get("field_names") != list(field_names)
            or (
                hybrid_parameters is not None
                and (
                    sidecar.get("hybrid_parameter_names")
                    != list(HYBRID_PARAMETER_NAMES)
                    or sidecar.get("hybrid_parameters")
                    != hybrid_parameters.tolist()
                )
            )
        ):
            raise ValueError("ADM volume sidecar does not match its binary")
        metadata.update(sidecar)
    return ADMVolume(
        times, lower, spacing, fields, metadata, coframes, hybrid_parameters
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--times", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--half-width", type=float, required=True)
    parser.add_argument("--cells", type=int, required=True)
    parser.add_argument("--halo", type=int, default=3)
    parser.add_argument("--secondary-mass", type=float, required=True)
    parser.add_argument("--secondary-chi", type=float, default=0.0)
    parser.add_argument("--chunk-size", type=int, default=65536)
    parser.add_argument("--hybrid-primary", action="store_true")
    parser.add_argument("--secondary-metric-fd-step", type=float)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    volume = build_volume(
        arguments.manifest,
        arguments.times,
        arguments.half_width,
        arguments.cells,
        arguments.halo,
        arguments.secondary_mass,
        arguments.secondary_chi,
        arguments.chunk_size,
        hybrid_primary=arguments.hybrid_primary,
        secondary_metric_fd_step=arguments.secondary_metric_fd_step,
    )
    report = write_binary(arguments.output, volume)
    print(arguments.output)
    print(
        "minimum_lapse="
        f"{report['adm_diagnostics']['minimum_lapse']:.8e}"
    )
    print(
        "minimum_spatial_metric_eigenvalue="
        f"{report['adm_diagnostics']['minimum_spatial_metric_eigenvalue']:.8e}"
    )


if __name__ == "__main__":
    main()
