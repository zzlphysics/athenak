#!/usr/bin/env python3
"""Extract a moving inner-coordinate CT worldtube from global GRMHD dumps.

The source is a time series of co-temporal AthenaK ``mhd_w_bcc`` and ``adm``
binary outputs.  Uniform grids are assembled densely; AMR/SMR outputs are
sampled sparsely on one explicitly selected leaf level and reject interpolation
stencils that cross a refinement boundary.  A cubic-Hermite affine frame maps
local events to the global chart.  Fluid four-velocity and the ideal MHD
Faraday two-form are reconstructed in the global ADM geometry, sampled with
linear spacetime interpolation, and transformed into the inner coordinates.

Face fluxes and edge EMFs are integrated directly as pulled-back differential
forms.  The endpoint face cochain and interval edge cochain are then projected
onto the closed cubical CT complex.  Projection corrections are diagnostics:
they must converge toward zero with source, target, and quadrature resolution.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
from dataclasses import dataclass, field
import json
import math
from pathlib import Path
import sys
from typing import Any, Callable, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
VIS_PYTHON = ROOT / "vis" / "python"
if str(VIS_PYTHON) not in sys.path:
    sys.path.insert(0, str(VIS_PYTHON))

import bin_convert  # noqa: E402
import compare_worldtube_closure as closure  # noqa: E402
import extract_static_taylor_worldtube as static  # noqa: E402
import worldtube_flux_emf as worldtube  # noqa: E402
import worldtube_frame as frame  # noqa: E402


INPUT_CLASSIFICATION = "athenak-emri-global-worldtube-series-v1"
FRAME_SERIES_CLASSIFICATION = "athenak-emri-affine-worldtube-frame-series-v1"
OUTPUT_CLASSIFICATION = "athenak-emri-global-to-local-worldtube-v1"
STATE_VARIABLES = ("dens", "velx", "vely", "velz", "press", "bcc1", "bcc2", "bcc3")
ADM_VARIABLES = static.ADM_VARIABLES
ALL_VARIABLES = STATE_VARIABLES + ADM_VARIABLES


def _finite_array(values: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if result.shape != shape or not np.isfinite(result).all():
        raise ValueError(f"{name} must be a finite array with shape {shape}")
    return result


def _finite_positive(value: object, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def _resolve_path(value: object, directory: Path, name: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a nonempty path")
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = directory / candidate
    return candidate.resolve(strict=True)


@dataclass(frozen=True)
class HermiteSeries:
    """Cubic-Hermite values with analytic first derivatives."""

    times: np.ndarray
    values: np.ndarray
    derivatives: np.ndarray
    name: str

    def __post_init__(self) -> None:
        times = np.asarray(self.times, dtype=np.float64)
        values = np.asarray(self.values, dtype=np.float64)
        derivatives = np.asarray(self.derivatives, dtype=np.float64)
        if (
            times.ndim != 1
            or times.size < 2
            or not np.isfinite(times).all()
            or np.any(np.diff(times) <= 0.0)
        ):
            raise ValueError(f"{self.name} times must be finite and strictly increasing")
        if (
            values.shape[0] != times.size
            or derivatives.shape != values.shape
            or not np.isfinite(values).all()
            or not np.isfinite(derivatives).all()
        ):
            raise ValueError(f"{self.name} Hermite arrays have incompatible shapes")
        object.__setattr__(self, "times", times)
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "derivatives", derivatives)

    def evaluate(self, time: float) -> tuple[np.ndarray, np.ndarray]:
        checked_time = float(time)
        tolerance = 64.0 * np.finfo(float).eps * max(
            1.0, abs(checked_time), abs(float(self.times[0])), abs(float(self.times[-1]))
        )
        if (
            checked_time < self.times[0] - tolerance
            or checked_time > self.times[-1] + tolerance
        ):
            raise ValueError(f"{self.name} time {checked_time:.17g} is outside its table")
        checked_time = min(max(checked_time, float(self.times[0])), float(self.times[-1]))
        interval = int(np.searchsorted(self.times, checked_time, side="right") - 1)
        interval = min(max(interval, 0), self.times.size - 2)
        left = float(self.times[interval])
        width = float(self.times[interval + 1] - self.times[interval])
        coordinate = (checked_time - left) / width
        s = coordinate
        h00 = 2.0 * s**3 - 3.0 * s**2 + 1.0
        h10 = s**3 - 2.0 * s**2 + s
        h01 = -2.0 * s**3 + 3.0 * s**2
        h11 = s**3 - s**2
        value = (
            h00 * self.values[interval]
            + h10 * width * self.derivatives[interval]
            + h01 * self.values[interval + 1]
            + h11 * width * self.derivatives[interval + 1]
        )
        derivative = (
            (6.0 * s**2 - 6.0 * s) * self.values[interval] / width
            + (3.0 * s**2 - 4.0 * s + 1.0) * self.derivatives[interval]
            + (-6.0 * s**2 + 6.0 * s) * self.values[interval + 1] / width
            + (3.0 * s**2 - 2.0 * s) * self.derivatives[interval + 1]
        )
        return value, derivative


@dataclass(frozen=True)
class AffineFrameSeries:
    times: np.ndarray
    worldline: HermiteSeries
    spatial_legs: HermiteSeries

    @classmethod
    def from_document(cls, document: dict[str, object]) -> "AffineFrameSeries":
        if document.get("classification") != FRAME_SERIES_CLASSIFICATION:
            raise ValueError("frame-series classification is missing or unsupported")
        times = np.asarray(document.get("times"), dtype=np.float64)
        if times.ndim != 1 or times.size < 2:
            raise ValueError("frame series requires at least two times")
        count = times.size
        worldline = HermiteSeries(
            times,
            _finite_array(document.get("worldline"), (count, 4), "worldline"),
            _finite_array(
                document.get("worldline_tangent"),
                (count, 4),
                "worldline_tangent",
            ),
            "worldline",
        )
        legs = HermiteSeries(
            times,
            _finite_array(
                document.get("spatial_legs"), (count, 4, 3), "spatial_legs"
            ),
            _finite_array(
                document.get("spatial_leg_derivative"),
                (count, 4, 3),
                "spatial_leg_derivative",
            ),
            "spatial legs",
        )
        result = cls(times=worldline.times, worldline=worldline, spatial_legs=legs)
        for sample_time in times:
            _, instantaneous = result.evaluate(float(sample_time))
            instantaneous.jacobian()
        return result

    def evaluate(self, time: float) -> tuple[np.ndarray, frame.AffineFrame]:
        position, tangent = self.worldline.evaluate(time)
        legs, leg_derivative = self.spatial_legs.evaluate(time)
        return position, frame.AffineFrame(tangent, legs, leg_derivative)


@dataclass(frozen=True)
class UniformSnapshot:
    time: float
    cycle: int
    lower: np.ndarray
    spacing: np.ndarray
    shape_xyz: tuple[int, int, int]
    values: np.ndarray
    state_path: Path
    adm_path: Path
    source_level: int = 0
    source_meshblock_count: int = 1
    available_leaf_levels: tuple[int, ...] = (0,)

    @property
    def source_storage(self) -> str:
        return "dense_uniform"

    def sample(self, positions: object) -> np.ndarray:
        points = np.asarray(positions, dtype=np.float64)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("snapshot positions must have finite shape (N,3)")
        count = points.shape[0]
        lower_index = np.empty((count, 3), dtype=np.int64)
        fraction = np.empty((count, 3), dtype=np.float64)
        for axis, cells in enumerate(self.shape_xyz):
            coordinate = (points[:, axis] - self.lower[axis]) / self.spacing[axis] - 0.5
            tolerance = 128.0 * np.finfo(float).eps * np.maximum(1.0, np.abs(coordinate))
            outside = (coordinate < -tolerance) | (coordinate > cells - 1 + tolerance)
            if np.any(outside):
                first = int(np.flatnonzero(outside)[0])
                raise ValueError(
                    "worldtube sample lies outside the source cell-center envelope: "
                    f"point={points[first].tolist()}, axis={axis + 1}, "
                    f"coordinate_index={coordinate[first]:.17g}"
                )
            coordinate = np.clip(coordinate, 0.0, cells - 1.0)
            index = np.floor(coordinate).astype(np.int64)
            index = np.minimum(index, cells - 2)
            lower_index[:, axis] = index
            fraction[:, axis] = coordinate - index

        ix, iy, iz = lower_index[:, 0], lower_index[:, 1], lower_index[:, 2]
        wx, wy, wz = fraction[:, 0], fraction[:, 1], fraction[:, 2]
        result = np.zeros((count, self.values.shape[0]), dtype=np.float64)
        for dz in (0, 1):
            weight_z = wz if dz else 1.0 - wz
            for dy in (0, 1):
                weight_y = wy if dy else 1.0 - wy
                for dx in (0, 1):
                    weight_x = wx if dx else 1.0 - wx
                    weight = weight_x * weight_y * weight_z
                    corner = self.values[:, iz + dz, iy + dy, ix + dx].T
                    result += weight[:, None] * corner
        return result


@dataclass(frozen=True)
class FixedLevelSnapshot:
    """Sparse same-level leaf blocks supporting cross-block interpolation."""

    time: float
    cycle: int
    lower: np.ndarray
    spacing: np.ndarray
    shape_xyz: tuple[int, int, int]
    block_shape_xyz: tuple[int, int, int]
    block_logical: np.ndarray
    values: np.ndarray
    state_path: Path
    adm_path: Path
    source_level: int
    available_leaf_levels: tuple[int, ...]
    block_lookup: dict[tuple[int, int, int], int] = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        lower = _finite_array(self.lower, (3,), "fixed-level lower bound")
        spacing = _finite_array(self.spacing, (3,), "fixed-level spacing")
        if np.any(spacing <= 0.0):
            raise ValueError("fixed-level spacing must be positive")
        if (
            len(self.shape_xyz) != 3
            or len(self.block_shape_xyz) != 3
            or min(self.shape_xyz) < 2
            or min(self.block_shape_xyz) < 1
        ):
            raise ValueError(
                "fixed-level grid and block shapes must be three-dimensional"
            )
        if not isinstance(self.source_level, int) or self.source_level < 0:
            raise ValueError("source_level must be a nonnegative integer")
        if self.source_level not in self.available_leaf_levels:
            raise ValueError("source_level is absent from available_leaf_levels")
        logical = np.asarray(self.block_logical)
        values = np.asarray(self.values)
        expected = (
            logical.shape[0],
            len(ALL_VARIABLES),
            self.block_shape_xyz[2],
            self.block_shape_xyz[1],
            self.block_shape_xyz[0],
        )
        if (
            logical.ndim != 2
            or logical.shape[1:] != (3,)
            or not np.issubdtype(logical.dtype, np.integer)
            or values.shape != expected
            or not np.isfinite(values).all()
        ):
            raise ValueError("fixed-level logical blocks or values are incompatible")
        lookup = {
            tuple(int(value) for value in row): index
            for index, row in enumerate(logical)
        }
        if len(lookup) != logical.shape[0]:
            raise ValueError("fixed-level snapshot contains duplicate logical blocks")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "spacing", spacing)
        object.__setattr__(self, "block_logical", logical.astype(np.int64))
        object.__setattr__(self, "values", np.asarray(values, dtype=np.float64))
        object.__setattr__(self, "block_lookup", lookup)

    @property
    def source_storage(self) -> str:
        return "sparse_fixed_leaf_level"

    @property
    def source_meshblock_count(self) -> int:
        return int(self.block_logical.shape[0])

    def sample(self, positions: object) -> np.ndarray:
        points = np.asarray(positions, dtype=np.float64)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("snapshot positions must have finite shape (N,3)")
        count = points.shape[0]
        lower_index = np.empty((count, 3), dtype=np.int64)
        fraction = np.empty((count, 3), dtype=np.float64)
        for axis, cells in enumerate(self.shape_xyz):
            coordinate = (points[:, axis] - self.lower[axis]) / self.spacing[axis] - 0.5
            tolerance = 128.0 * np.finfo(float).eps * np.maximum(
                1.0, np.abs(coordinate)
            )
            outside = (coordinate < -tolerance) | (
                coordinate > cells - 1 + tolerance
            )
            if np.any(outside):
                first = int(np.flatnonzero(outside)[0])
                raise ValueError(
                    "worldtube sample lies outside the source cell-center envelope: "
                    f"point={points[first].tolist()}, axis={axis + 1}, "
                    f"coordinate_index={coordinate[first]:.17g}"
                )
            coordinate = np.clip(coordinate, 0.0, cells - 1.0)
            index = np.floor(coordinate).astype(np.int64)
            index = np.minimum(index, cells - 2)
            lower_index[:, axis] = index
            fraction[:, axis] = coordinate - index

        result = np.zeros((count, self.values.shape[1]), dtype=np.float64)
        for dz in (0, 1):
            weight_z = fraction[:, 2] if dz else 1.0 - fraction[:, 2]
            for dy in (0, 1):
                weight_y = fraction[:, 1] if dy else 1.0 - fraction[:, 1]
                for dx in (0, 1):
                    weight_x = fraction[:, 0] if dx else 1.0 - fraction[:, 0]
                    cell = lower_index + np.asarray((dx, dy, dz))
                    block = cell // np.asarray(self.block_shape_xyz)
                    local = cell % np.asarray(self.block_shape_xyz)
                    unique, inverse = np.unique(block, axis=0, return_inverse=True)
                    block_ids = np.empty(unique.shape[0], dtype=np.int64)
                    for item, logical in enumerate(unique):
                        key = tuple(int(value) for value in logical)
                        if key not in self.block_lookup:
                            first = int(np.flatnonzero(inverse == item)[0])
                            raise ValueError(
                                "fixed-level interpolation stencil is not covered by "
                                f"a level-{self.source_level} leaf MeshBlock: "
                                f"point={points[first].tolist()}, logical_block={key}; "
                                "move the worldtube away from the AMR interface or "
                                "choose another uniformly covering source_level"
                            )
                        block_ids[item] = self.block_lookup[key]
                    corner = np.empty((count, self.values.shape[1]))
                    for item, block_id in enumerate(block_ids):
                        selected = inverse == item
                        corner[selected] = self.values[
                            block_id,
                            :,
                            local[selected, 2],
                            local[selected, 1],
                            local[selected, 0],
                        ]
                    weight = weight_x * weight_y * weight_z
                    result += weight[:, None] * corner
        return result


Snapshot = UniformSnapshot | FixedLevelSnapshot


@dataclass(frozen=True)
class SnapshotDescriptor:
    """Lightweight source metadata retained for a long snapshot series."""

    time: float
    cycle: int
    lower: np.ndarray
    spacing: np.ndarray
    shape_xyz: tuple[int, int, int]
    state_path: Path
    adm_path: Path
    source_level: int
    source_meshblock_count: int
    available_leaf_levels: tuple[int, ...]
    source_storage: str

    def __post_init__(self) -> None:
        lower = _finite_array(self.lower, (3,), "snapshot-descriptor lower bound")
        spacing = _finite_array(self.spacing, (3,), "snapshot-descriptor spacing")
        if np.any(spacing <= 0.0):
            raise ValueError("snapshot-descriptor spacing must be positive")
        if (
            not math.isfinite(self.time)
            or len(self.shape_xyz) != 3
            or min(self.shape_xyz) < 2
            or not isinstance(self.source_level, int)
            or self.source_level < 0
            or self.source_meshblock_count < 1
        ):
            raise ValueError("snapshot descriptor has invalid scalar metadata")
        if self.source_level not in self.available_leaf_levels:
            raise ValueError("descriptor source_level is absent from leaf levels")
        if self.source_storage not in (
            "dense_uniform",
            "sparse_fixed_leaf_level",
        ):
            raise ValueError("snapshot descriptor has unsupported source storage")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "spacing", spacing)


SnapshotMetadata = Snapshot | SnapshotDescriptor


def _validate_snapshot_sequence(snapshots: tuple[SnapshotMetadata, ...]) -> None:
    if len(snapshots) < 2:
        raise ValueError("global worldtube extraction requires at least two snapshots")
    times = np.asarray([snapshot.time for snapshot in snapshots])
    if not np.isfinite(times).all() or np.any(np.diff(times) <= 0.0):
        raise ValueError("global snapshot times must increase strictly")
    reference = snapshots[0]
    for snapshot in snapshots[1:]:
        if (
            snapshot.shape_xyz != reference.shape_xyz
            or snapshot.source_level != reference.source_level
            or not np.allclose(snapshot.lower, reference.lower, rtol=0.0, atol=2.0e-13)
            or not np.allclose(
                snapshot.spacing, reference.spacing, rtol=0.0, atol=2.0e-13
            )
        ):
            raise ValueError(
                "global snapshot source level or Cartesian geometry changes in time"
            )


def _sample_snapshot_events(
    times: np.ndarray,
    snapshot_at: Callable[[int], Snapshot],
    events: object,
) -> np.ndarray:
    points = np.asarray(events, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 4 or not np.isfinite(points).all():
        raise ValueError("global events must have finite shape (N,4)")
    tolerance = 128.0 * np.finfo(float).eps * np.maximum(1.0, np.abs(points[:, 0]))
    outside = (points[:, 0] < times[0] - tolerance) | (
        points[:, 0] > times[-1] + tolerance
    )
    if np.any(outside):
        first = int(np.flatnonzero(outside)[0])
        raise ValueError(
            f"mapped global time {points[first, 0]:.17g} lies outside snapshot "
            f"range [{times[0]:.17g}, {times[-1]:.17g}]"
        )
    event_times = np.clip(points[:, 0], times[0], times[-1])
    intervals = np.searchsorted(times, event_times, side="right") - 1
    intervals = np.clip(intervals, 0, times.size - 2)
    result = np.empty((points.shape[0], len(ALL_VARIABLES)), dtype=np.float64)
    for interval in np.unique(intervals):
        selected = intervals == interval
        left = snapshot_at(int(interval)).sample(points[selected, 1:])
        right = snapshot_at(int(interval) + 1).sample(points[selected, 1:])
        fraction = (event_times[selected] - times[interval]) / (
            times[interval + 1] - times[interval]
        )
        result[selected] = left + fraction[:, None] * (right - left)
    return result


@dataclass(frozen=True)
class SnapshotSeries:
    snapshots: tuple[Snapshot, ...]

    def __post_init__(self) -> None:
        _validate_snapshot_sequence(self.snapshots)

    @property
    def times(self) -> np.ndarray:
        return np.asarray([snapshot.time for snapshot in self.snapshots])

    @property
    def provenance_snapshots(self) -> tuple[Snapshot, ...]:
        return self.snapshots

    def sample(self, events: object) -> np.ndarray:
        return _sample_snapshot_events(self.times, self.snapshots.__getitem__, events)

    def loading_document(self) -> dict[str, object]:
        return {
            "mode": "eager",
            "snapshot_count": len(self.snapshots),
            "resident_snapshot_count": len(self.snapshots),
            "peak_cached_snapshots": len(self.snapshots),
        }


@dataclass
class LazySnapshotSeries:
    """Metadata-resident series with a bounded least-recently-used data cache."""

    descriptors: tuple[SnapshotDescriptor, ...]
    loader: Callable[[int], Snapshot]
    cache_size: int = 2
    _cache: OrderedDict[int, Snapshot] = field(
        init=False, default_factory=OrderedDict, repr=False
    )
    _cache_hits: int = field(init=False, default=0, repr=False)
    _cache_misses: int = field(init=False, default=0, repr=False)
    _cache_evictions: int = field(init=False, default=0, repr=False)
    _peak_cached_snapshots: int = field(init=False, default=0, repr=False)

    def __post_init__(self) -> None:
        if (
            not isinstance(self.cache_size, int)
            or isinstance(self.cache_size, bool)
            or self.cache_size < 2
        ):
            raise ValueError("snapshot cache size must be an integer of at least two")
        _validate_snapshot_sequence(self.descriptors)

    @property
    def times(self) -> np.ndarray:
        return np.asarray([snapshot.time for snapshot in self.descriptors])

    @property
    def provenance_snapshots(self) -> tuple[SnapshotDescriptor, ...]:
        return self.descriptors

    def _snapshot(self, index: int) -> Snapshot:
        if index in self._cache:
            self._cache_hits += 1
            snapshot = self._cache.pop(index)
            self._cache[index] = snapshot
            return snapshot
        if len(self._cache) >= self.cache_size:
            self._cache.popitem(last=False)
            self._cache_evictions += 1
        snapshot = self.loader(index)
        descriptor = self.descriptors[index]
        if (
            snapshot.time != descriptor.time
            or snapshot.cycle != descriptor.cycle
            or snapshot.state_path != descriptor.state_path
            or snapshot.adm_path != descriptor.adm_path
            or snapshot.shape_xyz != descriptor.shape_xyz
            or snapshot.source_level != descriptor.source_level
            or snapshot.source_meshblock_count != descriptor.source_meshblock_count
            or snapshot.available_leaf_levels != descriptor.available_leaf_levels
            or snapshot.source_storage != descriptor.source_storage
            or not np.allclose(
                snapshot.lower, descriptor.lower, rtol=0.0, atol=2.0e-13
            )
            or not np.allclose(
                snapshot.spacing, descriptor.spacing, rtol=0.0, atol=2.0e-13
            )
        ):
            raise RuntimeError(
                f"snapshot {index} metadata changed after its initial scan"
            )
        self._cache_misses += 1
        self._cache[index] = snapshot
        self._peak_cached_snapshots = max(
            self._peak_cached_snapshots, len(self._cache)
        )
        return snapshot

    def sample(self, events: object) -> np.ndarray:
        return _sample_snapshot_events(self.times, self._snapshot, events)

    def loading_document(self) -> dict[str, object]:
        return {
            "mode": "lazy_lru",
            "snapshot_count": len(self.descriptors),
            "cache_capacity": self.cache_size,
            "resident_snapshot_count": len(self._cache),
            "cache_hits": self._cache_hits,
            "cache_misses": self._cache_misses,
            "cache_evictions": self._cache_evictions,
            "peak_cached_snapshots": self._peak_cached_snapshots,
        }


SourceSeries = SnapshotSeries | LazySnapshotSeries


@dataclass
class SamplingAudit:
    maximum_four_velocity_normalization_error: float = 0.0
    maximum_ideal_mhd_residual: float = 0.0
    minimum_local_lapse: float = math.inf
    minimum_local_spatial_metric_eigenvalue: float = math.inf
    minimum_jacobian_absolute_determinant: float = math.inf
    mapped_global_time_minimum: float = math.inf
    mapped_global_time_maximum: float = -math.inf
    event_count: int = 0

    def update(
        self,
        jacobian: np.ndarray,
        metric: np.ndarray,
        four_velocity: np.ndarray,
        faraday: np.ndarray,
        global_times: np.ndarray,
    ) -> None:
        gamma = metric[:, 1:, 1:]
        beta_lower = metric[:, 0, 1:]
        beta = np.linalg.solve(gamma, beta_lower[..., None])[..., 0]
        lapse2 = np.einsum("ni,ni->n", beta_lower, beta) - metric[:, 0, 0]
        if np.any(lapse2 <= 0.0):
            raise RuntimeError("pulled-back local metric has non-positive lapse squared")
        eigenvalues = np.linalg.eigvalsh(gamma)
        if np.any(eigenvalues <= 0.0):
            raise RuntimeError(
                "pulled-back local spatial metric is not positive definite"
            )
        normalization = np.einsum(
            "ni,nij,nj->n", four_velocity, metric, four_velocity
        )
        ideal = np.einsum("nij,nj->ni", faraday, four_velocity)
        faraday_scale = np.maximum(np.max(np.abs(faraday), axis=(1, 2)), 1.0)
        velocity_scale = np.maximum(np.max(np.abs(four_velocity), axis=1), 1.0)
        ideal_relative = np.max(np.abs(ideal), axis=1) / (
            faraday_scale * velocity_scale
        )
        determinants = np.abs(np.linalg.det(jacobian))
        self.maximum_four_velocity_normalization_error = max(
            self.maximum_four_velocity_normalization_error,
            float(np.max(np.abs(normalization + 1.0))),
        )
        self.maximum_ideal_mhd_residual = max(
            self.maximum_ideal_mhd_residual, float(np.max(ideal_relative))
        )
        self.minimum_local_lapse = min(
            self.minimum_local_lapse, math.sqrt(float(np.min(lapse2)))
        )
        self.minimum_local_spatial_metric_eigenvalue = min(
            self.minimum_local_spatial_metric_eigenvalue, float(np.min(eigenvalues))
        )
        self.minimum_jacobian_absolute_determinant = min(
            self.minimum_jacobian_absolute_determinant, float(np.min(determinants))
        )
        self.mapped_global_time_minimum = min(
            self.mapped_global_time_minimum, float(np.min(global_times))
        )
        self.mapped_global_time_maximum = max(
            self.mapped_global_time_maximum, float(np.max(global_times))
        )
        self.event_count += int(global_times.size)

    def document(self) -> dict[str, float | int]:
        return {
            "sampled_spacetime_event_count": self.event_count,
            "maximum_four_velocity_normalization_error": (
                self.maximum_four_velocity_normalization_error
            ),
            "maximum_ideal_mhd_residual": self.maximum_ideal_mhd_residual,
            "minimum_local_lapse": self.minimum_local_lapse,
            "minimum_local_spatial_metric_eigenvalue": (
                self.minimum_local_spatial_metric_eigenvalue
            ),
            "minimum_jacobian_absolute_determinant": (
                self.minimum_jacobian_absolute_determinant
            ),
            "mapped_global_time_range": [
                self.mapped_global_time_minimum,
                self.mapped_global_time_maximum,
            ],
        }


def _fixed_level_snapshot(
    state: dict[str, object],
    adm: dict[str, object],
    adm_alignment: list[int],
    state_path: Path,
    adm_path: Path,
    source_level: int,
) -> FixedLevelSnapshot:
    logical_all = np.asarray(state["mb_logical"], dtype=np.int64)
    selected = np.flatnonzero(logical_all[:, 3] == source_level)
    if selected.size == 0:
        available = sorted(set(int(value) for value in logical_all[:, 3]))
        raise RuntimeError(
            f"source_level={source_level} has no leaf MeshBlocks; "
            f"available leaf levels are {available}"
        )
    variable_blocks = []
    try:
        for name in STATE_VARIABLES:
            variable_blocks.append(
                np.stack([state["mb_data"][name][index] for index in selected])
            )
        for name in ADM_VARIABLES:
            variable_blocks.append(
                np.stack(
                    [
                        adm["mb_data"][name][adm_alignment[index]]
                        for index in selected
                    ]
                )
            )
        values = np.stack(variable_blocks, axis=1)
    except ValueError as error:
        raise RuntimeError(
            "selected source-level MeshBlocks do not have a common output shape"
        ) from error

    block_shape_xyz = (
        int(values.shape[4]),
        int(values.shape[3]),
        int(values.shape[2]),
    )
    lower = np.asarray(
        (state["x1min"], state["x2min"], state["x3min"]), dtype=np.float64
    )
    upper = np.asarray(
        (state["x1max"], state["x2max"], state["x3max"]), dtype=np.float64
    )
    root_shape = np.asarray(
        (state["Nx1"], state["Nx2"], state["Nx3"]), dtype=np.int64
    )
    fine_shape = root_shape * (1 << source_level)
    spacing = (upper - lower) / fine_shape
    logical = logical_all[selected, :3]
    geometries = np.asarray(state["mb_geometry"], dtype=np.float64)[selected]
    for row, geometry in zip(logical, geometries, strict=True):
        cell_begin = row * np.asarray(block_shape_xyz)
        cell_end = cell_begin + np.asarray(block_shape_xyz)
        if np.any(cell_begin < 0) or np.any(cell_end > fine_shape):
            raise RuntimeError(
                "source-level logical MeshBlock lies outside its declared root grid"
            )
        expected_lower = lower + cell_begin * spacing
        expected_upper = lower + cell_end * spacing
        stored_lower = geometry[[0, 2, 4]]
        stored_upper = geometry[[1, 3, 5]]
        scale = max(1.0, float(np.max(np.abs(geometry))))
        if not np.allclose(
            stored_lower,
            expected_lower,
            rtol=0.0,
            atol=2.0e-12 * scale,
        ) or not np.allclose(
            stored_upper,
            expected_upper,
            rtol=0.0,
            atol=2.0e-12 * scale,
        ):
            raise RuntimeError(
                "source-level logical and physical MeshBlock geometries disagree"
            )
    return FixedLevelSnapshot(
        time=float(state["time"]),
        cycle=int(state["cycle"]),
        lower=lower,
        spacing=spacing,
        shape_xyz=tuple(int(value) for value in fine_shape),
        block_shape_xyz=block_shape_xyz,
        block_logical=logical,
        values=values,
        state_path=state_path,
        adm_path=adm_path,
        source_level=source_level,
        available_leaf_levels=tuple(
            sorted(set(int(value) for value in logical_all[:, 3]))
        ),
    )


def _read_snapshot_pair(
    entry: dict[str, object],
    directory: Path,
    load_variables: bool,
) -> tuple[Path, Path, dict[str, object], dict[str, object], list[int]]:
    state_path = _resolve_path(entry.get("state"), directory, "snapshot state")
    adm_path = _resolve_path(entry.get("adm"), directory, "snapshot ADM")
    state_variables = STATE_VARIABLES if load_variables else ()
    adm_variables = ADM_VARIABLES if load_variables else ()
    state = bin_convert.read_binary(str(state_path), variables=state_variables)
    adm = bin_convert.read_binary(str(adm_path), variables=adm_variables)
    missing_state = sorted(set(STATE_VARIABLES).difference(state["var_names"]))
    missing_adm = sorted(set(ADM_VARIABLES).difference(adm["var_names"]))
    if missing_state:
        raise RuntimeError(
            f"state dump is missing fields: {', '.join(missing_state)}"
        )
    if missing_adm:
        raise RuntimeError(f"ADM dump is missing fields: {', '.join(missing_adm)}")
    adm_alignment = static._aligned_adm_blocks(state, adm)
    time_scale = max(abs(float(state["time"])), abs(float(adm["time"])), 1.0)
    if int(state["cycle"]) != int(adm["cycle"]) or not math.isclose(
        float(state["time"]),
        float(adm["time"]),
        rel_tol=0.0,
        abs_tol=64.0 * np.finfo(float).eps * time_scale,
    ):
        raise RuntimeError("state and ADM snapshots are not co-temporal")
    if "time" in entry and not math.isclose(
        float(entry["time"]),
        float(state["time"]),
        rel_tol=0.0,
        abs_tol=64.0 * np.finfo(float).eps * time_scale,
    ):
        raise RuntimeError("manifest snapshot time does not match its binary dump")
    grid_labels = (
        "Nx1",
        "Nx2",
        "Nx3",
        "x1min",
        "x1max",
        "x2min",
        "x2max",
        "x3min",
        "x3max",
    )
    for label in grid_labels:
        if not math.isclose(
            float(state[label]), float(adm[label]), rel_tol=0.0, abs_tol=2.0e-13
        ):
            raise RuntimeError(f"state and ADM source grids differ in {label}")
    return state_path, adm_path, state, adm, adm_alignment


def _selected_source_level(
    state: dict[str, object], source_level: int | None
) -> tuple[int, tuple[int, ...]]:
    state_levels = np.asarray(state["mb_logical"], dtype=np.int64)[:, 3]
    available_levels = tuple(sorted(set(int(value) for value in state_levels)))
    if source_level is None and len(available_levels) > 1:
        raise RuntimeError(
            "AMR/SMR snapshot requires an explicit manifest source_level; "
            f"available leaf levels are {list(available_levels)}"
        )
    selected_level = available_levels[0] if source_level is None else source_level
    if selected_level not in available_levels:
        raise RuntimeError(
            f"source_level={selected_level} has no leaf MeshBlocks; "
            f"available leaf levels are {list(available_levels)}"
        )
    return selected_level, available_levels


def _scan_snapshot(
    entry: dict[str, object], directory: Path, source_level: int | None = None
) -> SnapshotDescriptor:
    state_path, adm_path, state, _, _ = _read_snapshot_pair(
        entry, directory, load_variables=False
    )
    selected_level, available_levels = _selected_source_level(state, source_level)
    root_shape = np.asarray(
        (state["Nx1"], state["Nx2"], state["Nx3"]), dtype=np.int64
    )
    if np.any(root_shape < 2):
        raise RuntimeError("worldtube interpolation requires at least two cells per axis")
    shape = root_shape * (1 << selected_level)
    lower = np.asarray(
        (state["x1min"], state["x2min"], state["x3min"]), dtype=np.float64
    )
    upper = np.asarray(
        (state["x1max"], state["x2max"], state["x3max"]), dtype=np.float64
    )
    sparse = source_level is not None or selected_level != 0
    levels = np.asarray(state["mb_logical"], dtype=np.int64)[:, 3]
    selected_count = int(np.count_nonzero(levels == selected_level))
    return SnapshotDescriptor(
        time=float(state["time"]),
        cycle=int(state["cycle"]),
        lower=lower,
        spacing=(upper - lower) / shape,
        shape_xyz=tuple(int(value) for value in shape),
        state_path=state_path,
        adm_path=adm_path,
        source_level=selected_level,
        source_meshblock_count=(selected_count if sparse else int(state["n_mbs"])),
        available_leaf_levels=available_levels,
        source_storage=("sparse_fixed_leaf_level" if sparse else "dense_uniform"),
    )


def _load_snapshot(
    entry: dict[str, object], directory: Path, source_level: int | None = None
) -> Snapshot:
    state_path, adm_path, state, adm, adm_alignment = _read_snapshot_pair(
        entry, directory, load_variables=True
    )
    static._check_dump(state, STATE_VARIABLES, "state")
    static._check_dump(adm, ADM_VARIABLES, "ADM")
    selected_level, available_levels = _selected_source_level(state, source_level)
    if source_level is not None or selected_level != 0:
        return _fixed_level_snapshot(
            state,
            adm,
            adm_alignment,
            state_path,
            adm_path,
            selected_level,
        )
    values = np.stack(
        [
            closure.assemble_uniform_grid(state, name)
            for name in STATE_VARIABLES
        ]
        + [
            closure.assemble_uniform_grid(adm, name)
            for name in ADM_VARIABLES
        ]
    )
    lower = np.asarray((state["x1min"], state["x2min"], state["x3min"]), dtype=float)
    upper = np.asarray((state["x1max"], state["x2max"], state["x3max"]), dtype=float)
    shape_xyz = (int(state["Nx1"]), int(state["Nx2"]), int(state["Nx3"]))
    if min(shape_xyz) < 2:
        raise RuntimeError("worldtube interpolation requires at least two cells per axis")
    spacing = (upper - lower) / np.asarray(shape_xyz)
    return UniformSnapshot(
        time=float(state["time"]),
        cycle=int(state["cycle"]),
        lower=lower,
        spacing=spacing,
        shape_xyz=shape_xyz,
        values=values,
        state_path=state_path,
        adm_path=adm_path,
        source_level=selected_level,
        source_meshblock_count=int(state["n_mbs"]),
        available_leaf_levels=available_levels,
    )


def _global_geometry(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    offset = len(STATE_VARIABLES)
    gamma = np.asarray(
        (
            (values[offset], values[offset + 1], values[offset + 2]),
            (values[offset + 1], values[offset + 3], values[offset + 4]),
            (values[offset + 2], values[offset + 4], values[offset + 5]),
        )
    )
    alpha = values[offset + 6]
    beta = values[offset + 7 : offset + 10]
    return gamma, alpha, beta


def _global_tensors(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    gamma, alpha, beta = _global_geometry(values)
    determinant = float(np.linalg.det(gamma))
    if (
        not alpha > 0.0
        or not determinant > 0.0
        or np.min(np.linalg.eigvalsh(gamma)) <= 0.0
    ):
        raise RuntimeError("interpolated global ADM geometry is invalid")
    primitive_u = values[1:4]
    lorentz = math.sqrt(1.0 + float(primitive_u @ gamma @ primitive_u))
    four_velocity = np.empty(4, dtype=np.float64)
    four_velocity[0] = lorentz / alpha
    four_velocity[1:] = primitive_u - beta * four_velocity[0]
    bcc = values[5:8]
    transport_velocity = four_velocity[1:] / four_velocity[0]
    electric = -np.cross(transport_velocity, bcc)
    faraday = frame.faraday_tensor(electric, bcc)
    metric = static.spacetime_metric_from_adm(gamma, float(alpha), beta)
    return metric, four_velocity, faraday


def transform_grmhd_sample(
    values: object,
    jacobian: object,
    global_length_in_local_units: float = 1.0,
    density_renormalization: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Transform one interpolated global primitive into local GRMHD variables.

    Returns ``(state, local_metric, local_four_velocity, local_faraday)``.  The
    state order is ``rho,u1,u2,u3,pgas,sqrt(gamma)B1..3``.
    """

    source = _finite_array(values, (len(ALL_VARIABLES),), "global GRMHD sample")
    mapping = _finite_array(jacobian, (4, 4), "spacetime jacobian")
    length_scale = _finite_positive(
        global_length_in_local_units, "global_length_in_local_units"
    )
    density_scale = _finite_positive(density_renormalization, "density_renormalization")
    if source[0] <= 0.0 or source[4] <= 0.0:
        raise RuntimeError("interpolated density and pressure must be positive")
    global_metric, global_four_velocity, global_faraday = _global_tensors(source)
    local_metric = length_scale**2 * (mapping.T @ global_metric @ mapping)
    local_four_velocity = np.linalg.solve(mapping, global_four_velocity) / length_scale
    local_faraday = (
        length_scale
        * math.sqrt(density_scale)
        * (mapping.T @ global_faraday @ mapping)
    )
    gamma = local_metric[1:, 1:]
    beta_lower = local_metric[0, 1:]
    try:
        beta = np.linalg.solve(gamma, beta_lower)
    except np.linalg.LinAlgError as error:
        raise RuntimeError("local spatial metric is singular") from error
    lapse2 = float(beta_lower @ beta - local_metric[0, 0])
    if lapse2 <= 0.0 or np.min(np.linalg.eigvalsh(gamma)) <= 0.0:
        raise RuntimeError("transformed local ADM geometry is invalid")
    lapse = math.sqrt(lapse2)
    if lapse * local_four_velocity[0] <= 0.0:
        raise RuntimeError("transformed fluid four-velocity is not future directed")
    primitive_u = local_four_velocity[1:] + beta * local_four_velocity[0]
    _, local_bcc = frame.electric_magnetic_from_faraday(local_faraday)
    state = np.asarray(
        (
            source[0] * density_scale / length_scale**2,
            primitive_u[0],
            primitive_u[1],
            primitive_u[2],
            source[4] * density_scale / length_scale**2,
            local_bcc[0],
            local_bcc[1],
            local_bcc[2],
        )
    )
    return state, local_metric, local_four_velocity, local_faraday


class GlobalWorldtubeSampler:
    def __init__(
        self,
        snapshots: SourceSeries,
        frames: AffineFrameSeries,
        global_length_in_local_units: float,
        density_renormalization: float,
    ):
        self.snapshots = snapshots
        self.frames = frames
        self.length_scale = _finite_positive(
            global_length_in_local_units, "global_length_in_local_units"
        )
        self.density_scale = _finite_positive(
            density_renormalization, "density_renormalization"
        )
        self.audit = SamplingAudit()

    def sample(
        self, local_time: float, local_positions: object
    ) -> tuple[np.ndarray, np.ndarray]:
        positions = np.asarray(local_positions, dtype=np.float64)
        if (
            positions.ndim != 2
            or positions.shape[1] != 3
            or not np.isfinite(positions).all()
        ):
            raise ValueError("local sample positions must have finite shape (N,3)")
        worldline, instantaneous = self.frames.evaluate(local_time)
        events = worldline[None, :] + positions @ instantaneous.spatial_legs.T
        values = self.snapshots.sample(events)
        count = positions.shape[0]
        states = np.empty((count, 8), dtype=np.float64)
        metrics = np.empty((count, 4, 4), dtype=np.float64)
        four_velocities = np.empty((count, 4), dtype=np.float64)
        faradays = np.empty((count, 4, 4), dtype=np.float64)
        jacobians = np.empty((count, 4, 4), dtype=np.float64)
        for sample in range(count):
            jacobian = instantaneous.jacobian(positions[sample])
            state, metric, four_velocity, faraday = transform_grmhd_sample(
                values[sample], jacobian, self.length_scale, self.density_scale
            )
            states[sample] = state
            metrics[sample] = metric
            four_velocities[sample] = four_velocity
            faradays[sample] = faraday
            jacobians[sample] = jacobian
        self.audit.update(
            jacobians, metrics, four_velocities, faradays, events[:, 0]
        )
        return states, faradays


@dataclass(frozen=True)
class CubeGeometry:
    center: np.ndarray
    half_width: float
    cells: int

    def __post_init__(self) -> None:
        center = _finite_array(self.center, (3,), "worldtube center")
        half_width = _finite_positive(self.half_width, "worldtube half_width")
        if not isinstance(self.cells, int) or self.cells < 1:
            raise ValueError("worldtube cells_per_edge must be a positive integer")
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "half_width", half_width)

    @property
    def spacing(self) -> float:
        return 2.0 * self.half_width / self.cells

    def face_positions(
        self,
        name: str,
        u_coordinates: np.ndarray,
        v_coordinates: np.ndarray,
        exterior_offset: float = 0.0,
    ) -> np.ndarray:
        orientation = worldtube.ORIENTATIONS[name]
        u_values, v_values = np.broadcast_arrays(u_coordinates, v_coordinates)
        result = np.broadcast_to(self.center, u_values.shape + (3,)).copy()
        result[..., orientation.normal_axis] += orientation.normal_sign * (
            self.half_width + exterior_offset
        )
        result[..., orientation.u_axis] += orientation.u_sign * u_values
        result[..., orientation.v_axis] += orientation.v_sign * v_values
        return result


def _quadrature(order: int) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(order, int) or order < 1 or order > 8:
        raise ValueError("quadrature_order must be an integer from one through eight")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return np.asarray(nodes), np.asarray(weights)


def _sample_endpoint_state(
    sampler: GlobalWorldtubeSampler,
    geometry: CubeGeometry,
    name: str,
    local_time: float,
) -> np.ndarray:
    centers = -geometry.half_width + (np.arange(geometry.cells) + 0.5) * geometry.spacing
    u, v = np.meshgrid(centers, centers, indexing="xy")
    positions = geometry.face_positions(
        name, u, v, exterior_offset=0.5 * geometry.spacing
    )
    state, _ = sampler.sample(local_time, positions.reshape(-1, 3))
    return state.reshape(geometry.cells, geometry.cells, 8).transpose(2, 0, 1)


def _sample_endpoint_flux(
    sampler: GlobalWorldtubeSampler,
    geometry: CubeGeometry,
    name: str,
    local_time: float,
    quadrature_order: int,
) -> np.ndarray:
    orientation = worldtube.ORIENTATIONS[name]
    nodes, weights = _quadrature(quadrature_order)
    centers = -geometry.half_width + (np.arange(geometry.cells) + 0.5) * geometry.spacing
    u = (
        centers[None, :, None, None]
        + 0.5 * geometry.spacing * nodes[None, None, None, :]
    )
    v = (
        centers[:, None, None, None]
        + 0.5 * geometry.spacing * nodes[None, None, :, None]
    )
    u, v = np.broadcast_arrays(u, v)
    positions = geometry.face_positions(name, u, v)
    _, faraday = sampler.sample(local_time, positions.reshape(-1, 3))
    magnetic = np.column_stack(
        (faraday[:, 2, 3], faraday[:, 3, 1], faraday[:, 1, 2])
    )
    integrand = orientation.normal_sign * magnetic[:, orientation.normal_axis]
    integrand = integrand.reshape(
        geometry.cells, geometry.cells, quadrature_order, quadrature_order
    )
    quadrature_weight = (
        (0.5 * geometry.spacing) ** 2
        * weights[:, None]
        * weights[None, :]
    )
    return np.sum(integrand * quadrature_weight[None, None, :, :], axis=(2, 3))


def _edge_positions(
    geometry: CubeGeometry,
    name: str,
    along_u: bool,
    quadrature_order: int,
) -> tuple[np.ndarray, int, int]:
    nodes, _ = _quadrature(quadrature_order)
    centers = -geometry.half_width + (np.arange(geometry.cells) + 0.5) * geometry.spacing
    vertices = -geometry.half_width + np.arange(geometry.cells + 1) * geometry.spacing
    if along_u:
        u = centers[None, :, None] + 0.5 * geometry.spacing * nodes[None, None, :]
        v = vertices[:, None, None]
        shape = (geometry.cells + 1, geometry.cells)
    else:
        u = vertices[None, :, None]
        v = centers[:, None, None] + 0.5 * geometry.spacing * nodes[None, None, :]
        shape = (geometry.cells, geometry.cells + 1)
    u, v = np.broadcast_arrays(u, v)
    return geometry.face_positions(name, u, v), shape[0], shape[1]


def _sample_interval_emf(
    sampler: GlobalWorldtubeSampler,
    geometry: CubeGeometry,
    name: str,
    left_time: float,
    right_time: float,
    quadrature_order: int,
    along_u: bool,
) -> np.ndarray:
    orientation = worldtube.ORIENTATIONS[name]
    nodes, weights = _quadrature(quadrature_order)
    positions, rows, columns = _edge_positions(
        geometry, name, along_u, quadrature_order
    )
    component = orientation.u_axis if along_u else orientation.v_axis
    direction_sign = orientation.u_sign if along_u else orientation.v_sign
    result = np.zeros((rows, columns), dtype=np.float64)
    midpoint = 0.5 * (left_time + right_time)
    half_dt = 0.5 * (right_time - left_time)
    spatial_weight = 0.5 * geometry.spacing * weights
    for time_node, time_weight in zip(nodes, weights, strict=True):
        local_time = midpoint + half_dt * time_node
        _, faraday = sampler.sample(local_time, positions.reshape(-1, 3))
        electric = -faraday[:, 0, 1:]
        line_values = direction_sign * electric[:, component]
        line_values = line_values.reshape(rows, columns, quadrature_order)
        line_integral = np.sum(
            line_values * spatial_weight[None, None, :], axis=2
        )
        result += 0.5 * time_weight * line_integral
    return result


def sample_worldtube(
    sampler: GlobalWorldtubeSampler,
    geometry: CubeGeometry,
    times: Iterable[float],
    quadrature_order: int = 2,
) -> tuple[dict[str, worldtube.FaceData], dict[str, object]]:
    times_array = worldtube.validate_times(times)
    _quadrature(quadrature_order)
    faces = {
        name: worldtube.FaceData(
            np.empty((times_array.size, 8, geometry.cells, geometry.cells)),
            np.empty((times_array.size, geometry.cells, geometry.cells)),
            np.empty(
                (times_array.size - 1, geometry.cells + 1, geometry.cells)
            ),
            np.empty(
                (times_array.size - 1, geometry.cells, geometry.cells + 1)
            ),
        )
        for name in worldtube.FACE_NAMES
    }
    for interval in range(times_array.size - 1):
        local_time = float(times_array[interval])
        for name in worldtube.FACE_NAMES:
            faces[name].cell_state[interval] = _sample_endpoint_state(
                sampler, geometry, name, local_time
            )
            faces[name].normal_flux[interval] = _sample_endpoint_flux(
                sampler, geometry, name, local_time, quadrature_order
            )
        for name in worldtube.FACE_NAMES:
            faces[name].emf_u[interval] = _sample_interval_emf(
                sampler,
                geometry,
                name,
                local_time,
                float(times_array[interval + 1]),
                quadrature_order,
                True,
            )
            faces[name].emf_v[interval] = _sample_interval_emf(
                sampler,
                geometry,
                name,
                local_time,
                float(times_array[interval + 1]),
                quadrature_order,
                False,
            )

    final = times_array.size - 1
    final_time = float(times_array[final])
    for name in worldtube.FACE_NAMES:
        faces[name].cell_state[final] = _sample_endpoint_state(
            sampler, geometry, name, final_time
        )
        faces[name].normal_flux[final] = _sample_endpoint_flux(
            sampler, geometry, name, final_time, quadrature_order
        )

    flux_projected, flux_diagnostics = frame.project_closed_surface_fluxes(
        times_array, faces
    )
    projected, emf_diagnostics = frame.project_moving_samples(
        times_array, flux_projected
    )
    return projected, {
        "classification": OUTPUT_CLASSIFICATION,
        "quadrature_order": quadrature_order,
        "raw_sampling": sampler.audit.document(),
        "snapshot_loading": sampler.snapshots.loading_document(),
        "closed_flux_projection": flux_diagnostics,
        "edge_emf_projection": emf_diagnostics,
    }


def _frame_document(value: object, directory: Path) -> dict[str, object]:
    if isinstance(value, dict):
        return value
    path = _resolve_path(value, directory, "frame series")
    document = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(document, dict):
        raise ValueError("frame-series file must contain a JSON object")
    return document


def extract_manifest(
    manifest_path: Path,
) -> tuple[
    np.ndarray,
    dict[str, worldtube.FaceData],
    dict[str, object],
    dict[str, object],
]:
    path = manifest_path.expanduser().resolve(strict=True)
    document = json.loads(path.read_text(encoding="utf-8"))
    if (
        not isinstance(document, dict)
        or document.get("classification") != INPUT_CLASSIFICATION
    ):
        raise ValueError("global-worldtube manifest classification is unsupported")
    hash_source_files = document.get("hash_source_files", False)
    if not isinstance(hash_source_files, bool):
        raise ValueError("hash_source_files must be true or false")
    entries = document.get("snapshots")
    if not isinstance(entries, list) or len(entries) < 2 or not all(
        isinstance(entry, dict) for entry in entries
    ):
        raise ValueError("manifest requires at least two snapshot objects")
    source_level_value = document.get("source_level")
    if source_level_value is None:
        source_level = None
    elif (
        not isinstance(source_level_value, int)
        or isinstance(source_level_value, bool)
        or source_level_value < 0
    ):
        raise ValueError("source_level must be a nonnegative integer when provided")
    else:
        source_level = source_level_value
    snapshot_cache_size = document.get("snapshot_cache_size", 2)
    if (
        not isinstance(snapshot_cache_size, int)
        or isinstance(snapshot_cache_size, bool)
        or snapshot_cache_size < 2
    ):
        raise ValueError("snapshot_cache_size must be an integer of at least two")
    snapshot_entries = tuple(entries)
    descriptors = tuple(
        _scan_snapshot(entry, path.parent, source_level)
        for entry in snapshot_entries
    )

    def load_snapshot(index: int) -> Snapshot:
        return _load_snapshot(snapshot_entries[index], path.parent, source_level)

    snapshots = LazySnapshotSeries(
        descriptors, load_snapshot, cache_size=snapshot_cache_size
    )
    frame_document = _frame_document(document.get("frame"), path.parent)
    frames = AffineFrameSeries.from_document(frame_document)
    sample_times = document.get("sample_times", frames.times)
    times = worldtube.validate_times(sample_times)
    target = document.get("worldtube")
    if not isinstance(target, dict):
        raise ValueError("manifest worldtube must be an object")
    cells = target.get("cells_per_edge")
    if not isinstance(cells, int):
        raise ValueError("worldtube cells_per_edge must be an integer")
    geometry = CubeGeometry(
        center=np.asarray(target.get("center", (0.0, 0.0, 0.0))),
        half_width=float(target.get("half_width")),
        cells=cells,
    )
    quadrature_order = document.get("quadrature_order", 2)
    if not isinstance(quadrature_order, int):
        raise ValueError("quadrature_order must be an integer")
    length_scale = _finite_positive(
        document.get("global_length_in_local_units", 1.0),
        "global_length_in_local_units",
    )
    density_scale = _finite_positive(
        document.get("density_renormalization", 1.0), "density_renormalization"
    )
    sampler = GlobalWorldtubeSampler(
        snapshots, frames, length_scale, density_scale
    )
    faces, diagnostics = sample_worldtube(
        sampler, geometry, times, quadrature_order
    )
    source_snapshots = []
    for snapshot in snapshots.provenance_snapshots:
        provenance = {
            "time": snapshot.time,
            "cycle": snapshot.cycle,
            "state": str(snapshot.state_path),
            "adm": str(snapshot.adm_path),
            "state_size_bytes": snapshot.state_path.stat().st_size,
            "adm_size_bytes": snapshot.adm_path.stat().st_size,
            "shape_xyz": list(snapshot.shape_xyz),
            "cell_spacing": snapshot.spacing.tolist(),
            "source_level": snapshot.source_level,
            "source_storage": snapshot.source_storage,
            "selected_leaf_meshblocks": snapshot.source_meshblock_count,
            "available_leaf_levels": list(snapshot.available_leaf_levels),
        }
        if hash_source_files:
            provenance["state_sha256"] = static.file_sha256(snapshot.state_path)
            provenance["adm_sha256"] = static.file_sha256(snapshot.adm_path)
        source_snapshots.append(provenance)
    metadata = {
        "producer_classification": OUTPUT_CLASSIFICATION,
        "source_manifest": str(path),
        "source_manifest_sha256": static.file_sha256(path),
        "frame_series_contract": frame_document,
        "source_file_hashes_recorded": hash_source_files,
        "source_snapshots": source_snapshots,
        "state_variables": ["rho", "u1", "u2", "u3", "pgas", "bcc1", "bcc2", "bcc3"],
        "dynamical_grmhd_state": True,
        "fluid_state_frame": "inner_coordinate",
        "velocity_representation": "inner Eulerian spatial four-velocity u^i",
        "magnetic_representation": (
            "inner sqrt(gamma) B^i from pulled Faraday two-form"
        ),
        "center": geometry.center.tolist(),
        "half_width": geometry.half_width,
        "grid_spacing": geometry.spacing,
        "state_sampling": (
            "outer-adjacent local cell center with spacetime interpolation"
        ),
        "flux_sampling": "Gauss integral of pulled-back Faraday spatial two-form",
        "emf_sampling": (
            "Gauss spacetime integral of pulled-back Faraday time-edge one-form"
        ),
        "global_length_in_local_units": length_scale,
        "density_renormalization": density_scale,
        "source_level": snapshots.provenance_snapshots[0].source_level,
        "sampling_diagnostics": diagnostics,
        "limitations": [
            "one-way coupling; the inner accretor does not modify the global disk",
            (
                "each trilinear stencil must remain on the selected fixed leaf level; "
                "coarse-fine interpolation is deliberately rejected"
            ),
            (
                "linear temporal and trilinear spatial interpolation of primitive "
                "plus ADM fields"
            ),
            (
                "inner total metric should approach the transformed global metric "
                "at the replay boundary"
            ),
            (
                "snapshot data use a bounded lazy cache, but the current binary "
                "reader still has a one-snapshot-pair parsing memory peak"
            ),
        ],
    }
    return times, faces, metadata, diagnostics


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--diagnostics", type=Path)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    times, faces, metadata, diagnostics = extract_manifest(arguments.manifest)
    validation = worldtube.write_worldtube(arguments.output, times, faces, metadata)
    report: dict[str, Any] = dict(diagnostics)
    report["validation"] = validation
    encoded = json.dumps(report, indent=2, sort_keys=True)
    if arguments.diagnostics is not None:
        arguments.diagnostics.parent.mkdir(parents=True, exist_ok=True)
        arguments.diagnostics.write_text(encoded + "\n", encoding="utf-8")
    print(encoded)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
