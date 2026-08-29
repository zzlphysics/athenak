#!/usr/bin/env python3
"""Audit global-snapshot worldtube coverage without loading fluid variables.

The preflight traverses the exact endpoint-state, face-flux, and spacetime-edge
quadrature loci used by ``extract_global_worldtube.py``.  It checks mapped time
coverage, the source cell-center envelope, every trilinear stencil on both time
endpoints, fixed-level AMR/SMR MeshBlock coverage, and the affine Jacobian.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from itertools import product
import json
import math
from pathlib import Path

import numpy as np

import extract_global_worldtube as extract
import worldtube_flux_emf as worldtube


CLASSIFICATION = "athenak-emri-global-worldtube-preflight-v1"


@dataclass
class SnapshotCoverage:
    descriptor: extract.SnapshotDescriptor
    queried_position_count: int = 0
    minimum_cell_center_envelope_margin_cells: float = math.inf
    minimum_additional_stencil_halo_cells: int = 2**31 - 1
    _block_lookup: set[tuple[int, int, int]] = field(
        init=False, default_factory=set, repr=False
    )
    _halo_cache: dict[tuple[int, int, int], int] = field(
        init=False, default_factory=dict, repr=False
    )

    def __post_init__(self) -> None:
        self._block_lookup = {
            tuple(int(value) for value in logical)
            for logical in self.descriptor.block_logical
        }

    def _covered_cell(self, cell: np.ndarray) -> bool:
        if np.any(cell < 0) or np.any(cell >= self.descriptor.shape_xyz):
            return False
        if self.descriptor.source_storage == "dense_uniform":
            return True
        block = cell // np.asarray(self.descriptor.block_shape_xyz)
        return tuple(int(value) for value in block) in self._block_lookup

    def _additional_halo(self, lower: np.ndarray) -> int:
        if self.descriptor.source_storage == "dense_uniform":
            upper_margin = np.asarray(self.descriptor.shape_xyz) - lower - 2
            return int(min(np.min(lower), np.min(upper_margin)))
        key = tuple(int(value) for value in lower)
        if key in self._halo_cache:
            return self._halo_cache[key]
        halo = 0
        while True:
            candidate_halo = halo + 1
            lower_bound = lower - candidate_halo
            upper_bound = lower + 1 + candidate_halo
            ranges = (
                range(int(lower_bound[axis]), int(upper_bound[axis]) + 1)
                for axis in range(3)
            )
            layer_covered = True
            for values in product(*ranges):
                cell = np.asarray(values)
                on_boundary = np.any(cell == lower_bound) or np.any(
                    cell == upper_bound
                )
                if on_boundary and not self._covered_cell(cell):
                    layer_covered = False
                    break
            if not layer_covered:
                self._halo_cache[key] = halo
                return halo
            halo = candidate_halo

    def audit(self, positions: np.ndarray) -> None:
        points = np.asarray(positions, dtype=np.float64)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("preflight positions must have finite shape (N,3)")
        shape = np.asarray(self.descriptor.shape_xyz)
        coordinate = (
            (points - self.descriptor.lower) / self.descriptor.spacing - 0.5
        )
        tolerance = 128.0 * np.finfo(float).eps * np.maximum(
            1.0, np.abs(coordinate)
        )
        outside = (coordinate < -tolerance) | (
            coordinate > shape[None, :] - 1 + tolerance
        )
        if np.any(outside):
            sample, axis = np.argwhere(outside)[0]
            raise ValueError(
                "worldtube preflight leaves the source cell-center envelope: "
                f"snapshot_time={self.descriptor.time:.17g}, "
                f"point={points[sample].tolist()}, axis={axis + 1}, "
                f"coordinate_index={coordinate[sample, axis]:.17g}"
            )
        coordinate = np.clip(coordinate, 0.0, shape - 1.0)
        lower = np.floor(coordinate).astype(np.int64)
        lower = np.minimum(lower, shape - 2)
        envelope_margin = np.minimum(coordinate, shape - 1 - coordinate)
        self.minimum_cell_center_envelope_margin_cells = min(
            self.minimum_cell_center_envelope_margin_cells,
            float(np.min(envelope_margin)),
        )
        for stencil in np.unique(lower, axis=0):
            for dx in (0, 1):
                for dy in (0, 1):
                    for dz in (0, 1):
                        cell = stencil + np.asarray((dx, dy, dz))
                        if not self._covered_cell(cell):
                            block = cell // np.asarray(
                                self.descriptor.block_shape_xyz
                            )
                            raise ValueError(
                                "worldtube preflight stencil crosses an unavailable "
                                f"level-{self.descriptor.source_level} leaf block: "
                                f"snapshot_time={self.descriptor.time:.17g}, "
                                f"cell={cell.tolist()}, "
                                f"logical_block={tuple(int(x) for x in block)}"
                            )
            self.minimum_additional_stencil_halo_cells = min(
                self.minimum_additional_stencil_halo_cells,
                self._additional_halo(stencil),
            )
        self.queried_position_count += points.shape[0]

    def document(self) -> dict[str, object]:
        queried = self.queried_position_count > 0
        return {
            "time": self.descriptor.time,
            "cycle": self.descriptor.cycle,
            "source_level": self.descriptor.source_level,
            "source_storage": self.descriptor.source_storage,
            "selected_leaf_meshblocks": self.descriptor.source_meshblock_count,
            "queried_position_count": self.queried_position_count,
            "minimum_cell_center_envelope_margin_cells": (
                self.minimum_cell_center_envelope_margin_cells if queried else None
            ),
            "minimum_additional_stencil_halo_cells": (
                self.minimum_additional_stencil_halo_cells if queried else None
            ),
        }


@dataclass
class CoverageSeries:
    descriptors: tuple[extract.SnapshotDescriptor, ...]
    snapshots: tuple[SnapshotCoverage, ...] = field(init=False)
    mapped_global_time_minimum: float = math.inf
    mapped_global_time_maximum: float = -math.inf

    def __post_init__(self) -> None:
        extract._validate_snapshot_sequence(self.descriptors)
        self.snapshots = tuple(
            SnapshotCoverage(descriptor) for descriptor in self.descriptors
        )

    @property
    def times(self) -> np.ndarray:
        return np.asarray([descriptor.time for descriptor in self.descriptors])

    def audit(self, events: np.ndarray) -> None:
        points = np.asarray(events, dtype=np.float64)
        if points.ndim != 2 or points.shape[1] != 4 or not np.isfinite(points).all():
            raise ValueError("preflight global events must have finite shape (N,4)")
        times = self.times
        tolerance = 128.0 * np.finfo(float).eps * np.maximum(
            1.0, np.abs(points[:, 0])
        )
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
        for interval in np.unique(intervals):
            selected = intervals == interval
            positions = points[selected, 1:]
            self.snapshots[int(interval)].audit(positions)
            self.snapshots[int(interval) + 1].audit(positions)
        self.mapped_global_time_minimum = min(
            self.mapped_global_time_minimum, float(np.min(points[:, 0]))
        )
        self.mapped_global_time_maximum = max(
            self.mapped_global_time_maximum, float(np.max(points[:, 0]))
        )


@dataclass
class CoverageSampler:
    snapshots: CoverageSeries
    frames: extract.AffineFrameSeries
    local_event_count: int = 0
    minimum_jacobian_absolute_determinant: float = math.inf
    maximum_jacobian_condition_number: float = 0.0

    def sample(
        self, local_time: float, local_positions: object
    ) -> tuple[np.ndarray, np.ndarray]:
        positions = np.asarray(local_positions, dtype=np.float64)
        worldline, instantaneous = self.frames.evaluate(local_time)
        events = worldline[None, :] + positions @ instantaneous.spatial_legs.T
        self.snapshots.audit(events)
        for position in positions:
            jacobian = instantaneous.jacobian(position)
            self.minimum_jacobian_absolute_determinant = min(
                self.minimum_jacobian_absolute_determinant,
                abs(float(np.linalg.det(jacobian))),
            )
            self.maximum_jacobian_condition_number = max(
                self.maximum_jacobian_condition_number,
                float(np.linalg.cond(jacobian)),
            )
        self.local_event_count += positions.shape[0]
        return (
            np.zeros((positions.shape[0], 8)),
            np.zeros((positions.shape[0], 4, 4)),
        )


def _sampling_contract(
    scan: extract.SnapshotManifestScan,
) -> tuple[
    extract.AffineFrameSeries,
    np.ndarray,
    extract.CubeGeometry,
    int,
    dict[str, object],
]:
    document = scan.document
    frame_document = extract._frame_document(document.get("frame"), scan.path.parent)
    frames = extract.AffineFrameSeries.from_document(frame_document)
    times = worldtube.validate_times(document.get("sample_times", frames.times))
    target = document.get("worldtube")
    if not isinstance(target, dict):
        raise ValueError("manifest worldtube must be an object")
    cells = target.get("cells_per_edge")
    if not isinstance(cells, int) or isinstance(cells, bool):
        raise ValueError("worldtube cells_per_edge must be an integer")
    geometry = extract.CubeGeometry(
        center=np.asarray(target.get("center", (0.0, 0.0, 0.0))),
        half_width=float(target.get("half_width")),
        cells=cells,
    )
    quadrature_order = document.get("quadrature_order", 2)
    if not isinstance(quadrature_order, int) or isinstance(quadrature_order, bool):
        raise ValueError("quadrature_order must be an integer")
    extract._quadrature(quadrature_order)
    return frames, times, geometry, quadrature_order, frame_document


def audit_manifest(manifest_path: Path) -> dict[str, object]:
    scan = extract.scan_snapshot_manifest(manifest_path)
    frames, times, geometry, quadrature_order, frame_document = _sampling_contract(
        scan
    )
    coverage = CoverageSeries(scan.descriptors)
    sampler = CoverageSampler(coverage, frames)
    for interval in range(times.size - 1):
        local_time = float(times[interval])
        for name in worldtube.FACE_NAMES:
            extract._sample_endpoint_state(sampler, geometry, name, local_time)
            extract._sample_endpoint_flux(
                sampler, geometry, name, local_time, quadrature_order
            )
        for name in worldtube.FACE_NAMES:
            extract._sample_interval_emf(
                sampler,
                geometry,
                name,
                local_time,
                float(times[interval + 1]),
                quadrature_order,
                True,
            )
            extract._sample_interval_emf(
                sampler,
                geometry,
                name,
                local_time,
                float(times[interval + 1]),
                quadrature_order,
                False,
            )
    final_time = float(times[-1])
    for name in worldtube.FACE_NAMES:
        extract._sample_endpoint_state(sampler, geometry, name, final_time)
        extract._sample_endpoint_flux(
            sampler, geometry, name, final_time, quadrature_order
        )
    queried = [
        snapshot
        for snapshot in coverage.snapshots
        if snapshot.queried_position_count > 0
    ]
    minimum_envelope = min(
        snapshot.minimum_cell_center_envelope_margin_cells
        for snapshot in queried
    )
    minimum_halo = min(
        snapshot.minimum_additional_stencil_halo_cells for snapshot in queried
    )
    length_scale = float(
        scan.document.get("global_length_in_local_units", 1.0)
    )
    if not math.isfinite(length_scale) or length_scale <= 0.0:
        raise ValueError("global_length_in_local_units must be finite and positive")
    warnings = []
    if minimum_halo < 1:
        warnings.append(
            "some interpolation stencils have no additional selected-level cell halo"
        )
    if minimum_envelope < 1.0:
        warnings.append(
            "some samples lie less than one source cell from the center envelope"
        )
    return {
        "classification": CLASSIFICATION,
        "passed": True,
        "source_manifest": str(scan.path),
        "source_manifest_sha256": extract.static.file_sha256(scan.path),
        "source_snapshot_count": len(scan.descriptors),
        "selected_source_level": scan.descriptors[0].source_level,
        "sample_times": times.tolist(),
        "worldtube": {
            "center": geometry.center.tolist(),
            "half_width": geometry.half_width,
            "cells_per_edge": geometry.cells,
            "grid_spacing": geometry.spacing,
        },
        "quadrature_order": quadrature_order,
        "frame_series_contract": frame_document,
        "sampled_local_event_count": sampler.local_event_count,
        "mapped_global_time_range": [
            coverage.mapped_global_time_minimum,
            coverage.mapped_global_time_maximum,
        ],
        "minimum_jacobian_absolute_determinant": (
            sampler.minimum_jacobian_absolute_determinant
        ),
        "minimum_scaled_jacobian_absolute_determinant": (
            sampler.minimum_jacobian_absolute_determinant * length_scale**4
        ),
        "maximum_jacobian_condition_number": (
            sampler.maximum_jacobian_condition_number
        ),
        "minimum_cell_center_envelope_margin_cells": minimum_envelope,
        "minimum_additional_stencil_halo_cells": minimum_halo,
        "warnings": warnings,
        "snapshots": [snapshot.document() for snapshot in coverage.snapshots],
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    try:
        report = audit_manifest(arguments.manifest)
    except (OSError, RuntimeError, ValueError) as error:
        report = {
            "classification": CLASSIFICATION,
            "passed": False,
            "source_manifest": str(arguments.manifest.expanduser().resolve()),
            "error": str(error),
        }
    encoded = json.dumps(report, indent=2, sort_keys=True)
    if arguments.output is not None:
        output = arguments.output.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded + "\n", encoding="utf-8")
    print(encoded)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
