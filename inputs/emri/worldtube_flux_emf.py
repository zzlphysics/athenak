#!/usr/bin/env python3
"""Validate and conservatively regrid cubical CT worldtube data.

The format stores outward normal magnetic *fluxes* on face cells and time-averaged,
line-integrated electric fields on face edges.  Its tensor-product transfer operators
commute with the face-local discrete Faraday curl.  This module is the host-side format
and validation layer; AthenaK outer extraction and inner RK-stage replay use the same
topological contract but are implemented separately.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import struct
from typing import Iterable
import zlib

import numpy as np


CLASSIFICATION = "athenak-emri-cubical-flux-emf-worldtube-v1"
OUTER_STREAM_CLASSIFICATION = "athenak-emri-outer-worldtube-stream-v1"
INNER_BINARY_CLASSIFICATION = "athenak-emri-inner-worldtube-binary-v1"
INNER_BINARY_MAGIC = b"AEMRIWTBIN0001\x00\x00"
INNER_BINARY_VERSION = 1
INNER_BINARY_HEADER = struct.Struct("<16sIIII4dQ")
FACE_NAMES = ("x1m", "x1p", "x2m", "x2p", "x3m", "x3p")


@dataclass(frozen=True)
class FaceOrientation:
    """Right-handed local face axes with ``u cross v = outward normal``."""

    normal_axis: int
    normal_sign: int
    u_axis: int
    u_sign: int
    v_axis: int
    v_sign: int


ORIENTATIONS = {
    "x1p": FaceOrientation(0, 1, 1, 1, 2, 1),
    "x1m": FaceOrientation(0, -1, 1, 1, 2, -1),
    "x2p": FaceOrientation(1, 1, 2, 1, 0, 1),
    "x2m": FaceOrientation(1, -1, 2, 1, 0, -1),
    "x3p": FaceOrientation(2, 1, 0, 1, 1, 1),
    "x3m": FaceOrientation(2, -1, 0, 1, 1, -1),
}


@dataclass
class FaceData:
    """One face sampled at state times and over the intervening time intervals.

    ``cell_state`` has shape ``(nt, nvar, nv, nu)`` and contains intensive face-cell
    values. ``normal_flux`` has shape ``(nt, nv, nu)`` and is outward integrated flux.
    ``emf_u`` has shape ``(nt-1, nv+1, nu)``; ``emf_v`` has shape
    ``(nt-1, nv, nu+1)``.  EMFs are line integrals oriented along local +u and +v.
    """

    cell_state: np.ndarray
    normal_flux: np.ndarray
    emf_u: np.ndarray
    emf_v: np.ndarray


def _finite_array(values: np.ndarray, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if not np.isfinite(array).all():
        raise ValueError(f"{name} contains a non-finite value")
    return array


def validate_times(times: Iterable[float]) -> np.ndarray:
    result = _finite_array(np.asarray(tuple(times), dtype=np.float64), "times")
    if result.ndim != 1 or result.size < 2:
        raise ValueError("worldtube requires at least two one-dimensional times")
    if not np.all(np.diff(result) > 0.0):
        raise ValueError("worldtube times must increase strictly")
    return result


def validate_face(face: FaceData, times: np.ndarray, name: str) -> FaceData:
    state = _finite_array(face.cell_state, f"{name} cell_state")
    flux = _finite_array(face.normal_flux, f"{name} normal_flux")
    emf_u = _finite_array(face.emf_u, f"{name} emf_u")
    emf_v = _finite_array(face.emf_v, f"{name} emf_v")
    if state.ndim != 4 or state.shape[0] != times.size or state.shape[1] < 1:
        raise ValueError(f"{name} cell_state must have shape (nt,nvar,nv,nu)")
    nt, _, nv, nu = state.shape
    expected = {
        "normal_flux": (nt, nv, nu),
        "emf_u": (nt - 1, nv + 1, nu),
        "emf_v": (nt - 1, nv, nu + 1),
    }
    for field_name, array in (
        ("normal_flux", flux),
        ("emf_u", emf_u),
        ("emf_v", emf_v),
    ):
        if array.shape != expected[field_name]:
            raise ValueError(
                f"{name} {field_name} has shape {array.shape}, expected "
                f"{expected[field_name]}"
            )
    return FaceData(state, flux, emf_u, emf_v)


def faraday_update(
    normal_flux: np.ndarray,
    emf_u: np.ndarray,
    emf_v: np.ndarray,
    dt: float,
) -> np.ndarray:
    """Advance one face's integrated normal flux by its oriented edge circulation."""

    normal_flux = np.asarray(normal_flux, dtype=np.float64)
    emf_u = np.asarray(emf_u, dtype=np.float64)
    emf_v = np.asarray(emf_v, dtype=np.float64)
    nv, nu = normal_flux.shape
    if emf_u.shape != (nv + 1, nu) or emf_v.shape != (nv, nu + 1):
        raise ValueError("EMF shapes are incompatible with normal_flux")
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("Faraday timestep must be finite and positive")
    circulation = (
        emf_u[:-1, :]
        + emf_v[:, 1:]
        - emf_u[1:, :]
        - emf_v[:, :-1]
    )
    return normal_flux - dt * circulation


def faraday_residuals(face: FaceData, times: np.ndarray) -> np.ndarray:
    residuals = []
    for index, dt in enumerate(np.diff(times)):
        predicted = faraday_update(
            face.normal_flux[index], face.emf_u[index], face.emf_v[index], float(dt)
        )
        residuals.append(face.normal_flux[index + 1] - predicted)
    return np.asarray(residuals)


def _overlap_lengths(source_cells: int, target_cells: int) -> np.ndarray:
    if source_cells < 1 or target_cells < 1:
        raise ValueError("source and target grid sizes must be positive")
    source_edges = np.linspace(0.0, 1.0, source_cells + 1)
    target_edges = np.linspace(0.0, 1.0, target_cells + 1)
    result = np.zeros((target_cells, source_cells), dtype=np.float64)
    for target in range(target_cells):
        for source in range(source_cells):
            result[target, source] = max(
                0.0,
                min(target_edges[target + 1], source_edges[source + 1])
                - max(target_edges[target], source_edges[source]),
            )
    return result


def integral_transfer_matrix(source_cells: int, target_cells: int) -> np.ndarray:
    """Map source cell integrals into target cell integrals, preserving their sum."""

    return _overlap_lengths(source_cells, target_cells) * source_cells


def average_transfer_matrix(source_cells: int, target_cells: int) -> np.ndarray:
    """Map piecewise-constant source averages into target cell averages."""

    return _overlap_lengths(source_cells, target_cells) * target_cells


def nodal_linear_matrix(source_cells: int, target_cells: int) -> np.ndarray:
    """Interpolate source nodes to target nodes on the same unit interval."""

    source_nodes = np.linspace(0.0, 1.0, source_cells + 1)
    target_nodes = np.linspace(0.0, 1.0, target_cells + 1)
    result = np.zeros((target_cells + 1, source_cells + 1), dtype=np.float64)
    for target, coordinate in enumerate(target_nodes):
        if coordinate >= 1.0:
            result[target, -1] = 1.0
            continue
        lower = min(int(math.floor(coordinate * source_cells)), source_cells - 1)
        fraction = (coordinate - source_nodes[lower]) * source_cells
        result[target, lower] = 1.0 - fraction
        result[target, lower + 1] = fraction
    return result


def resample_face(face: FaceData, target_nv: int, target_nu: int) -> FaceData:
    """Mimetically regrid one face while retaining the discrete Faraday relation."""

    _, _, source_nv, source_nu = face.cell_state.shape
    flux_u = integral_transfer_matrix(source_nu, target_nu)
    flux_v = integral_transfer_matrix(source_nv, target_nv)
    average_u = average_transfer_matrix(source_nu, target_nu)
    average_v = average_transfer_matrix(source_nv, target_nv)
    node_u = nodal_linear_matrix(source_nu, target_nu)
    node_v = nodal_linear_matrix(source_nv, target_nv)

    state = np.einsum(
        "av,tnvu,bu->tnab", average_v, face.cell_state, average_u, optimize=True
    )
    normal_flux = np.einsum(
        "av,tvu,bu->tab", flux_v, face.normal_flux, flux_u, optimize=True
    )
    # u-directed line integrals are conservatively transferred along u and linearly
    # interpolated across v nodes.  The v-directed form is the tensor-product dual.
    emf_u = np.einsum(
        "av,tvu,bu->tab", node_v, face.emf_u, flux_u, optimize=True
    )
    emf_v = np.einsum(
        "av,tvu,bu->tab", flux_v, face.emf_v, node_u, optimize=True
    )
    return FaceData(state, normal_flux, emf_u, emf_v)


def _canonical_edge(
    face_name: str, side: str, face: FaceData
) -> tuple[tuple[tuple[int, int], tuple[int, int]], np.ndarray]:
    orientation = ORIENTATIONS[face_name]
    if side in ("vmin", "vmax"):
        tangent_axis = orientation.u_axis
        tangent_sign = orientation.u_sign
        v_local_sign = -1 if side == "vmin" else 1
        fixed = (
            (orientation.normal_axis, orientation.normal_sign),
            (orientation.v_axis, orientation.v_sign * v_local_sign),
        )
        values = face.emf_u[:, 0 if side == "vmin" else -1, :]
    elif side in ("umin", "umax"):
        tangent_axis = orientation.v_axis
        tangent_sign = orientation.v_sign
        u_local_sign = -1 if side == "umin" else 1
        fixed = (
            (orientation.normal_axis, orientation.normal_sign),
            (orientation.u_axis, orientation.u_sign * u_local_sign),
        )
        values = face.emf_v[:, :, 0 if side == "umin" else -1]
    else:
        raise ValueError(f"unknown face edge side {side}")
    key = tuple(sorted(fixed))
    free_axis = ({0, 1, 2} - {key[0][0], key[1][0]}).pop()
    if free_axis != tangent_axis:
        raise RuntimeError("invalid cubical face orientation table")
    if tangent_sign < 0:
        values = -values[:, ::-1]
    return key, values


def edge_consistency_residuals(faces: dict[str, FaceData]) -> dict[str, float]:
    occurrences: dict[tuple[tuple[int, int], tuple[int, int]], list[np.ndarray]] = {}
    for face_name in FACE_NAMES:
        for side in ("umin", "umax", "vmin", "vmax"):
            key, values = _canonical_edge(face_name, side, faces[face_name])
            occurrences.setdefault(key, []).append(values)
    if len(occurrences) != 12 or any(len(values) != 2 for values in occurrences.values()):
        raise RuntimeError("cubical edge topology did not produce twelve paired edges")
    result: dict[str, float] = {}
    for key, values in occurrences.items():
        if values[0].shape != values[1].shape:
            raise ValueError(f"worldtube faces use incompatible resolution at edge {key}")
        label = ",".join(f"x{axis + 1}{'p' if sign > 0 else 'm'}" for axis, sign in key)
        result[label] = float(np.max(np.abs(values[0] - values[1])))
    return result


def validate_worldtube(
    times: Iterable[float],
    faces: dict[str, FaceData],
    absolute_tolerance: float = 1.0e-12,
    relative_tolerance: float = 1.0e-10,
) -> dict[str, object]:
    times_array = validate_times(times)
    if set(faces) != set(FACE_NAMES):
        missing = sorted(set(FACE_NAMES).difference(faces))
        extra = sorted(set(faces).difference(FACE_NAMES))
        raise ValueError(f"worldtube face mismatch: missing={missing}, extra={extra}")
    checked = {
        name: validate_face(faces[name], times_array, name) for name in FACE_NAMES
    }
    faraday_max: dict[str, float] = {}
    faraday_scale: dict[str, float] = {}
    for name, face in checked.items():
        residual = faraday_residuals(face, times_array)
        faraday_max[name] = float(np.max(np.abs(residual)))
        faraday_scale[name] = max(float(np.max(np.abs(face.normal_flux))), 1.0)
    edge_max = edge_consistency_residuals(checked)
    edge_scale = max(
        1.0,
        *(float(np.max(np.abs(face.emf_u))) for face in checked.values()),
        *(float(np.max(np.abs(face.emf_v))) for face in checked.values()),
    )
    net_flux = np.asarray(
        [
            sum(float(np.sum(checked[name].normal_flux[index])) for name in FACE_NAMES)
            for index in range(times_array.size)
        ]
    )
    flux_scale = max(
        1.0,
        max(
            sum(float(np.sum(np.abs(checked[name].normal_flux[index])))
                for name in FACE_NAMES)
            for index in range(times_array.size)
        ),
    )
    failures = []
    for name in FACE_NAMES:
        limit = absolute_tolerance + relative_tolerance * faraday_scale[name]
        if faraday_max[name] > limit:
            failures.append(f"{name} Faraday residual {faraday_max[name]:.6g}>{limit:.6g}")
    edge_limit = absolute_tolerance + relative_tolerance * edge_scale
    for name, residual in edge_max.items():
        if residual > edge_limit:
            failures.append(f"edge {name} residual {residual:.6g}>{edge_limit:.6g}")
    flux_limit = absolute_tolerance + relative_tolerance * flux_scale
    if float(np.max(np.abs(net_flux))) > flux_limit:
        failures.append(
            f"closed-surface flux {float(np.max(np.abs(net_flux))):.6g}>{flux_limit:.6g}"
        )
    if failures:
        raise ValueError("invalid flux/EMF worldtube: " + "; ".join(failures))
    return {
        "classification": CLASSIFICATION,
        "sample_count": int(times_array.size),
        "time_range": [float(times_array[0]), float(times_array[-1])],
        "maximum_faraday_residual_by_face": faraday_max,
        "maximum_shared_edge_emf_residual": max(edge_max.values()),
        "maximum_closed_surface_flux": float(np.max(np.abs(net_flux))),
        "absolute_tolerance": absolute_tolerance,
        "relative_tolerance": relative_tolerance,
    }


def write_worldtube(
    path: Path,
    times: Iterable[float],
    faces: dict[str, FaceData],
    metadata: dict[str, object] | None = None,
) -> dict[str, object]:
    times_array = validate_times(times)
    diagnostics = validate_worldtube(times_array, faces)
    document = dict(metadata or {})
    document["classification"] = CLASSIFICATION
    document["face_orientations"] = {
        name: ORIENTATIONS[name].__dict__ for name in FACE_NAMES
    }
    arrays: dict[str, np.ndarray] = {
        "metadata_json": np.asarray(json.dumps(document, sort_keys=True)),
        "times": times_array,
    }
    for name in FACE_NAMES:
        face = validate_face(faces[name], times_array, name)
        arrays[f"{name}_cell_state"] = face.cell_state
        arrays[f"{name}_normal_flux"] = face.normal_flux
        arrays[f"{name}_emf_u"] = face.emf_u
        arrays[f"{name}_emf_v"] = face.emf_v
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as stream:
        np.savez(stream, **arrays)
    return diagnostics


def read_worldtube(path: Path) -> tuple[np.ndarray, dict[str, FaceData], dict[str, object]]:
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata_json"]))
        if metadata.get("classification") != CLASSIFICATION:
            raise ValueError("worldtube classification is missing or unsupported")
        recorded = metadata.get("face_orientations")
        expected = {name: ORIENTATIONS[name].__dict__ for name in FACE_NAMES}
        if recorded != expected:
            raise ValueError("worldtube face orientation metadata is incompatible")
        times = np.asarray(archive["times"], dtype=np.float64)
        faces = {
            name: FaceData(
                np.asarray(archive[f"{name}_cell_state"], dtype=np.float64),
                np.asarray(archive[f"{name}_normal_flux"], dtype=np.float64),
                np.asarray(archive[f"{name}_emf_u"], dtype=np.float64),
                np.asarray(archive[f"{name}_emf_v"], dtype=np.float64),
            )
            for name in FACE_NAMES
        }
    return validate_times(times), faces, metadata


def resample_worldtube(
    times: np.ndarray,
    faces: dict[str, FaceData],
    target_cells: int,
) -> dict[str, FaceData]:
    if target_cells < 1:
        raise ValueError("target_cells must be positive")
    result = {
        name: resample_face(
            validate_face(faces[name], times, name), target_cells, target_cells
        )
        for name in FACE_NAMES
    }
    validate_worldtube(times, result)
    return result


def _read_exact_binary(path: Path, dtype: str, shape: tuple[int, ...]) -> np.ndarray:
    expected = math.prod(shape)
    values = np.fromfile(path, dtype=np.dtype(dtype))
    if values.size != expected:
        raise ValueError(
            f"{path} contains {values.size} values, expected exactly {expected}"
        )
    result = np.asarray(values.reshape(shape), dtype=np.float64)
    if not np.isfinite(result).all():
        raise ValueError(f"{path} contains a non-finite value")
    return result


def read_outer_stream(
    manifest_path: Path,
) -> tuple[np.ndarray, dict[str, FaceData], dict[str, object]]:
    """Load the strict streaming files written by AthenaK's outer recorder."""

    with manifest_path.open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    if manifest.get("classification") != OUTER_STREAM_CLASSIFICATION:
        raise ValueError("outer stream classification is missing or unsupported")
    if manifest.get("target_classification") != CLASSIFICATION:
        raise ValueError("outer stream targets an incompatible worldtube schema")
    if manifest.get("complete") is not True:
        raise ValueError("outer stream is incomplete; finalize the AthenaK run first")
    if manifest.get("binary_dtype") != "<f8":
        raise ValueError("outer stream must use little-endian float64 binary data")
    nt = manifest.get("nt")
    ninterval = manifest.get("ninterval")
    nvar = manifest.get("nvar")
    cells = manifest.get("cells_per_face_axis")
    if (
        not isinstance(nt, int)
        or nt < 2
        or not isinstance(ninterval, int)
        or ninterval != nt - 1
        or not isinstance(nvar, int)
        or nvar < 1
        or not isinstance(cells, int)
        or cells < 1
    ):
        raise ValueError("outer stream dimensions or interval count are invalid")
    face_files = manifest.get("faces")
    if not isinstance(face_files, dict) or set(face_files) != set(FACE_NAMES):
        raise ValueError("outer stream must name exactly the six cubical faces")

    directory = manifest_path.resolve().parent

    def resolve_file(filename: object) -> Path:
        if not isinstance(filename, str) or not filename:
            raise ValueError("outer stream contains an invalid binary filename")
        candidate = (directory / filename).resolve()
        if candidate.parent != directory:
            raise ValueError("outer stream binary filenames must stay beside the manifest")
        return candidate

    times = _read_exact_binary(
        resolve_file(manifest.get("times_file")), "<f8", (nt,)
    )
    times = validate_times(times)
    faces: dict[str, FaceData] = {}
    expected_keys = {"cell_state", "normal_flux", "emf_u", "emf_v"}
    for name in FACE_NAMES:
        files = face_files[name]
        if not isinstance(files, dict) or set(files) != expected_keys:
            raise ValueError(f"outer stream file table for {name} is invalid")
        faces[name] = FaceData(
            cell_state=_read_exact_binary(
                resolve_file(files["cell_state"]), "<f8", (nt, nvar, cells, cells)
            ),
            normal_flux=_read_exact_binary(
                resolve_file(files["normal_flux"]), "<f8", (nt, cells, cells)
            ),
            emf_u=_read_exact_binary(
                resolve_file(files["emf_u"]),
                "<f8",
                (ninterval, cells + 1, cells),
            ),
            emf_v=_read_exact_binary(
                resolve_file(files["emf_v"]),
                "<f8",
                (ninterval, cells, cells + 1),
            ),
        )
    validate_worldtube(times, faces)
    metadata = {
        "outer_stream_manifest": manifest,
        "state_variables": manifest.get("state_variables", []),
        "center": manifest.get("center"),
        "half_width": manifest.get("half_width"),
        "grid_spacing": manifest.get("grid_spacing"),
    }
    return times, faces, metadata


def _crc32_update(value: int, payload: bytes) -> int:
    return zlib.crc32(payload, value) & 0xFFFFFFFF


def write_inner_binary(
    path: Path,
    times: Iterable[float],
    faces: dict[str, FaceData],
    metadata: dict[str, object] | None = None,
) -> dict[str, object]:
    """Write a strict little-endian stream for bounded-slab C++ device replay."""

    times_array = validate_times(times)
    diagnostics = validate_worldtube(times_array, faces)
    checked = {
        name: validate_face(faces[name], times_array, name) for name in FACE_NAMES
    }
    reference_shape = checked[FACE_NAMES[0]].cell_state.shape
    _, nvar, cells_v, cells_u = reference_shape
    if cells_u != cells_v:
        raise ValueError("inner replay currently requires square face grids")
    for name in FACE_NAMES:
        if checked[name].cell_state.shape != reference_shape:
            raise ValueError("inner replay requires one shared face resolution and nvar")
    document = dict(metadata or {})
    center = document.get("center", (0.0, 0.0, 0.0))
    half_width = document.get("half_width", 1.0)
    if (
        not isinstance(center, (list, tuple))
        or len(center) != 3
        or not all(math.isfinite(float(value)) for value in center)
        or not math.isfinite(float(half_width))
        or float(half_width) <= 0.0
    ):
        raise ValueError("inner replay metadata requires finite center and half_width")

    path.parent.mkdir(parents=True, exist_ok=True)
    checksum = 0
    with path.open("w+b") as stream:
        stream.write(
            INNER_BINARY_HEADER.pack(
                INNER_BINARY_MAGIC,
                INNER_BINARY_VERSION,
                cells_u,
                nvar,
                times_array.size,
                float(center[0]),
                float(center[1]),
                float(center[2]),
                float(half_width),
                0,
            )
        )

        def append(array: np.ndarray) -> None:
            nonlocal checksum
            payload = np.ascontiguousarray(array, dtype="<f8").tobytes(order="C")
            checksum = _crc32_update(checksum, payload)
            stream.write(payload)

        append(times_array)
        for name in FACE_NAMES:
            append(checked[name].cell_state)
            append(checked[name].normal_flux)
            append(checked[name].emf_u)
            append(checked[name].emf_v)
        stream.seek(0)
        stream.write(
            INNER_BINARY_HEADER.pack(
                INNER_BINARY_MAGIC,
                INNER_BINARY_VERSION,
                cells_u,
                nvar,
                times_array.size,
                float(center[0]),
                float(center[1]),
                float(center[2]),
                float(half_width),
                checksum,
            )
        )

    sidecar = {
        "classification": INNER_BINARY_CLASSIFICATION,
        "binary_file": path.name,
        "format": "<16sIIII4dQ followed by little-endian float64 arrays",
        "crc32_payload": f"{checksum:08x}",
        "face_order": list(FACE_NAMES),
        "array_order_per_face": ["cell_state", "normal_flux", "emf_u", "emf_v"],
        "nt": int(times_array.size),
        "nvar": int(nvar),
        "cells_per_face_axis": int(cells_u),
        "center": [float(value) for value in center],
        "half_width": float(half_width),
        "state_variables": document.get("state_variables", []),
        "source_metadata": document,
    }
    sidecar_path = path.with_suffix(path.suffix + ".json")
    sidecar_path.write_text(json.dumps(sidecar, indent=2, sort_keys=True), encoding="utf-8")
    diagnostics["binary_file"] = str(path)
    diagnostics["payload_crc32"] = f"{checksum:08x}"
    return diagnostics


def read_inner_binary(
    path: Path,
) -> tuple[np.ndarray, dict[str, FaceData], dict[str, object]]:
    """Reference reader used to verify files before the C++ replay consumes them."""

    payload = path.read_bytes()
    if len(payload) < INNER_BINARY_HEADER.size:
        raise ValueError("inner replay binary is shorter than its header")
    header = INNER_BINARY_HEADER.unpack_from(payload)
    magic, version, cells, nvar, nt, *remaining = header
    center = remaining[:3]
    half_width = remaining[3]
    expected_checksum = remaining[4]
    if magic != INNER_BINARY_MAGIC or version != INNER_BINARY_VERSION:
        raise ValueError("inner replay magic or version is unsupported")
    if cells < 1 or nvar < 1 or nt < 2:
        raise ValueError("inner replay dimensions are invalid")
    binary_payload = payload[INNER_BINARY_HEADER.size :]
    checksum = _crc32_update(0, binary_payload)
    if checksum != expected_checksum:
        raise ValueError("inner replay payload checksum mismatch")
    values = np.frombuffer(binary_payload, dtype="<f8")
    cursor = 0

    def take(shape: tuple[int, ...]) -> np.ndarray:
        nonlocal cursor
        count = math.prod(shape)
        if cursor + count > values.size:
            raise ValueError("inner replay binary ended inside an array")
        result = np.asarray(values[cursor : cursor + count].reshape(shape), dtype=np.float64)
        cursor += count
        if not np.isfinite(result).all():
            raise ValueError("inner replay binary contains a non-finite value")
        return result

    times = validate_times(take((nt,)))
    faces = {}
    for name in FACE_NAMES:
        faces[name] = FaceData(
            cell_state=take((nt, nvar, cells, cells)),
            normal_flux=take((nt, cells, cells)),
            emf_u=take((nt - 1, cells + 1, cells)),
            emf_v=take((nt - 1, cells, cells + 1)),
        )
    if cursor != values.size:
        raise ValueError("inner replay binary has unclassified trailing data")
    diagnostics = validate_worldtube(times, faces)
    metadata = {
        "classification": INNER_BINARY_CLASSIFICATION,
        "center": list(center),
        "half_width": half_width,
        "payload_crc32": f"{checksum:08x}",
        "diagnostics": diagnostics,
    }
    sidecar_path = path.with_suffix(path.suffix + ".json")
    if sidecar_path.exists():
        with sidecar_path.open(encoding="utf-8") as stream:
            sidecar = json.load(stream)
        sidecar_checksum = sidecar.get(
            "crc32_payload", sidecar.get("payload_crc32")
        )
        if (
            sidecar.get("classification") != INNER_BINARY_CLASSIFICATION
            or sidecar_checksum != f"{checksum:08x}"
        ):
            raise ValueError("inner replay sidecar does not match its binary payload")
        metadata["state_variables"] = sidecar.get("state_variables", [])
        metadata["source_metadata"] = sidecar.get("source_metadata", {})
    return times, faces, metadata


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("input", type=Path)
    resample_parser = subparsers.add_parser("resample")
    resample_parser.add_argument("input", type=Path)
    resample_parser.add_argument("output", type=Path)
    resample_parser.add_argument("--cells-per-face", type=int, required=True)
    pack_parser = subparsers.add_parser("pack-outer")
    pack_parser.add_argument("manifest", type=Path)
    pack_parser.add_argument("output", type=Path)
    prepare_parser = subparsers.add_parser("prepare-inner")
    prepare_parser.add_argument("input", type=Path)
    prepare_parser.add_argument("output", type=Path)
    inspect_parser = subparsers.add_parser("inspect-inner")
    inspect_parser.add_argument("input", type=Path)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    if arguments.command == "pack-outer":
        times, faces, metadata = read_outer_stream(arguments.manifest)
        diagnostics = write_worldtube(arguments.output, times, faces, metadata)
        print(json.dumps(diagnostics, indent=2))
        return
    if arguments.command == "inspect-inner":
        times, faces, metadata = read_inner_binary(arguments.input)
        diagnostics = validate_worldtube(times, faces)
        diagnostics["metadata"] = metadata
        print(json.dumps(diagnostics, indent=2))
        return
    times, faces, metadata = read_worldtube(arguments.input)
    if arguments.command == "prepare-inner":
        diagnostics = write_inner_binary(arguments.output, times, faces, metadata)
    elif arguments.command == "validate":
        diagnostics = validate_worldtube(times, faces)
    else:
        faces = resample_worldtube(times, faces, arguments.cells_per_face)
        diagnostics = write_worldtube(arguments.output, times, faces, metadata)
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
