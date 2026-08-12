#!/usr/bin/env python3
"""Stream-compare physical fields in two closed AthenaK restart files.

The comparator is intentionally narrower than the opaque restart-payload
comparator.  It understands MHD restarts with an optional prescribed ADM
metric and reports numerical differences separately for

* active and ghost cells of ``mhd::MHD::u0``;
* active faces of each constrained-transport magnetic-field component; and
* active and ghost cells of the 17 prescribed-ADM variables.

Only one MeshBlock from each file is resident at a time.  ParameterInput and
tree metadata are parsed by ``read_athenak_restart_metadata.py``; array shapes,
the MHD scalar count, and the per-MeshBlock ``data_size`` record are then
derived and validated independently.  Files containing Hydro, radiation,
Z4c, or the turbulence driver are rejected rather than silently interpreting
an unsupported slot layout.

Inputs must be closed single-file restarts.  File identity, size, mtime, and
ctime are checked before and after streaming so a changing checkpoint fails
closed.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import struct
import sys
from typing import BinaryIO, Callable, Sequence

import numpy as np


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
if str(SCRIPT_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

from read_athenak_restart_metadata import (  # noqa: E402
    RestartMetadata,
    read_restart_metadata,
)


ADM_NAMES = (
    "adm_gxx", "adm_gxy", "adm_gxz", "adm_gyy", "adm_gyz", "adm_gzz",
    "adm_Kxx", "adm_Kxy", "adm_Kxz", "adm_Kyy", "adm_Kyz", "adm_Kzz",
    "adm_psi4", "adm_alpha", "adm_betax", "adm_betay", "adm_betaz",
)
UNSUPPORTED_FIELD_BLOCKS = ("hydro", "radiation", "z4c", "turb_driving")


@dataclass(frozen=True)
class FileSignature:
    device: int
    inode: int
    size: int
    mtime_ns: int
    ctime_ns: int


@dataclass(frozen=True)
class RestartFieldLayout:
    real_dtype: np.dtype
    unsigned_dtype: np.dtype
    real_bytes: int
    nghost: int
    block_shape: tuple[int, int, int]
    stored_shape: tuple[int, int, int]
    active_starts: tuple[int, int, int]
    active_stops: tuple[int, int, int]
    nmhd_base: int
    nscalars: int
    nmhd: int
    nadm: int
    mhd_elements: int
    face_elements: tuple[int, int, int]
    adm_elements: int
    mhd_offset: int
    face_offsets: tuple[int, int, int]
    adm_offset: int
    data_size: int


def _signature_from_stat(stat_result: os.stat_result) -> FileSignature:
    return FileSignature(
        device=stat_result.st_dev,
        inode=stat_result.st_ino,
        size=stat_result.st_size,
        mtime_ns=stat_result.st_mtime_ns,
        ctime_ns=stat_result.st_ctime_ns,
    )


def _file_signature(path: Path) -> FileSignature:
    return _signature_from_stat(path.stat())


def _assert_path_unchanged(
    path: Path, expected: FileSignature, checkpoint: str
) -> None:
    if _file_signature(path) != expected:
        raise RuntimeError(
            f"restart changed {checkpoint}: {path}; compare only closed files"
        )


def _assert_stream_unchanged(
    stream: BinaryIO, path: Path, expected: FileSignature, checkpoint: str
) -> None:
    if _signature_from_stat(os.fstat(stream.fileno())) != expected:
        raise RuntimeError(
            f"open restart changed {checkpoint}: {path}; compare only closed files"
        )


def _parameter_int(metadata: RestartMetadata, block: str, key: str) -> int:
    try:
        return int(metadata.parameters[block][key])
    except (KeyError, ValueError) as exc:
        raise ValueError(
            f"restart parameter <{block}>/{key} is missing or invalid"
        ) from exc


def _parameter_float(metadata: RestartMetadata, block: str, key: str) -> float:
    try:
        return float(metadata.parameters[block][key])
    except (KeyError, ValueError) as exc:
        raise ValueError(
            f"restart parameter <{block}>/{key} is missing or invalid"
        ) from exc


def _derive_layout(metadata: RestartMetadata) -> RestartFieldLayout:
    parameters = metadata.parameters
    unsupported = [name for name in UNSUPPORTED_FIELD_BLOCKS if name in parameters]
    if unsupported:
        raise ValueError(
            "restart field comparator does not support physics block(s): "
            + ", ".join(f"<{name}>" for name in unsupported)
        )
    if "mhd" not in parameters:
        raise ValueError("restart has no <mhd> block")

    nghost = _parameter_int(metadata, "mesh", "nghost")
    if nghost < 0:
        raise ValueError("restart has a negative ghost-zone count")
    block_shape = tuple(
        _parameter_int(metadata, "meshblock", f"nx{axis}")
        for axis in range(1, 4)
    )
    if any(size < 1 for size in block_shape):
        raise ValueError("restart has a non-positive MeshBlock dimension")

    stored_shape = tuple(
        size + 2 * nghost if size > 1 else 1 for size in block_shape
    )
    active_starts = tuple(nghost if size > 1 else 0 for size in block_shape)
    active_stops = tuple(
        start + size for start, size in zip(active_starts, block_shape)
    )

    eos = parameters["mhd"].get("eos")
    if eos == "ideal":
        nmhd_base = 5
    elif eos == "isothermal":
        nmhd_base = 4
    else:
        raise ValueError(f"unsupported or missing <mhd>/eos: {eos!r}")
    try:
        nscalars = int(parameters["mhd"].get("nscalars", "0"))
    except ValueError as exc:
        raise ValueError("restart parameter <mhd>/nscalars is invalid") from exc
    if nscalars < 0:
        raise ValueError("restart has a negative MHD scalar count")
    nmhd = nmhd_base + nscalars
    nadm = len(ADM_NAMES) if "adm" in parameters else 0

    nout1, nout2, nout3 = stored_shape
    cell_elements = nout1 * nout2 * nout3
    mhd_elements = nmhd * cell_elements
    face_elements = (
        (nout1 + 1) * nout2 * nout3,
        nout1 * (nout2 + 1) * nout3,
        nout1 * nout2 * (nout3 + 1),
    )
    adm_elements = nadm * cell_elements
    real_bytes = metadata.real_bytes
    if real_bytes not in (4, 8):
        raise ValueError(f"unsupported Real size: {real_bytes}")
    byte_prefix = "<" if metadata.byte_order == "little" else ">"
    real_dtype = np.dtype(f"{byte_prefix}f{real_bytes}")
    unsigned_dtype = np.dtype(f"{byte_prefix}u{real_bytes}")

    mhd_offset = 0
    face1_offset = mhd_elements * real_bytes
    face2_offset = face1_offset + face_elements[0] * real_bytes
    face3_offset = face2_offset + face_elements[1] * real_bytes
    adm_offset = face3_offset + face_elements[2] * real_bytes
    data_size = adm_offset + adm_elements * real_bytes
    return RestartFieldLayout(
        real_dtype=real_dtype,
        unsigned_dtype=unsigned_dtype,
        real_bytes=real_bytes,
        nghost=nghost,
        block_shape=block_shape,
        stored_shape=stored_shape,
        active_starts=active_starts,
        active_stops=active_stops,
        nmhd_base=nmhd_base,
        nscalars=nscalars,
        nmhd=nmhd,
        nadm=nadm,
        mhd_elements=mhd_elements,
        face_elements=face_elements,
        adm_elements=adm_elements,
        mhd_offset=mhd_offset,
        face_offsets=(face1_offset, face2_offset, face3_offset),
        adm_offset=adm_offset,
        data_size=data_size,
    )


def _topology_signature(metadata: RestartMetadata) -> tuple[object, ...]:
    return (
        metadata.byte_order,
        metadata.real_bytes,
        metadata.nmb_total,
        metadata.root_level,
        _parameter_int(metadata, "mesh", "nghost"),
        *(
            _parameter_int(metadata, block, f"nx{axis}")
            for block in ("mesh", "meshblock")
            for axis in range(1, 4)
        ),
        *(
            _parameter_float(metadata, "mesh", f"x{axis}{bound}")
            for axis in range(1, 4)
            for bound in ("min", "max")
        ),
        metadata.locations,
    )


def _mhd_names(layout: RestartFieldLayout) -> tuple[str, ...]:
    base = ("IDN", "IM1", "IM2", "IM3")
    if layout.nmhd_base == 5:
        base += ("IEN",)
    return base + tuple(f"scalar_{index}" for index in range(layout.nscalars))


def _empty_accumulator() -> dict[str, object]:
    return {
        "elements": 0,
        "finite_pairs": 0,
        "exact_equal": True,
        "within_tolerance": True,
        "max_abs": 0.0,
        "max_abs_location": None,
        "l1_sum": 0.0,
        "reference_l1_sum": 0.0,
        "left_nonfinite": {
            "nan": 0,
            "positive_infinity": 0,
            "negative_infinity": 0,
        },
        "right_nonfinite": {
            "nan": 0,
            "positive_infinity": 0,
            "negative_infinity": 0,
        },
        "nonfinite_class_mismatches": 0,
        "first_difference": None,
    }


def _nonfinite_counts(values: np.ndarray) -> dict[str, int]:
    return {
        "nan": int(np.count_nonzero(np.isnan(values))),
        "positive_infinity": int(np.count_nonzero(np.isposinf(values))),
        "negative_infinity": int(np.count_nonzero(np.isneginf(values))),
    }


def _json_value(value: np.generic | float) -> float | str:
    converted = float(value)
    if math.isnan(converted):
        return "NaN"
    if math.isinf(converted):
        return "+Infinity" if converted > 0 else "-Infinity"
    return converted


def _accumulate_region(
    accumulator: dict[str, object],
    left: np.ndarray,
    right: np.ndarray,
    selector: np.ndarray,
    describe_index: Callable[[tuple[int, ...]], dict[str, object]],
    *,
    unsigned_dtype: np.dtype,
    rtol: float,
    atol: float,
) -> None:
    if left.shape != right.shape:
        raise ValueError(f"array shape mismatch: {left.shape} != {right.shape}")
    selected = np.broadcast_to(selector, left.shape)
    selected_count = int(np.count_nonzero(selected))
    accumulator["elements"] = int(accumulator["elements"]) + selected_count
    if selected_count == 0:
        return

    left_selected = left[selected]
    right_selected = right[selected]
    for side, values in (("left", left_selected), ("right", right_selected)):
        counts = _nonfinite_counts(values)
        stored = accumulator[f"{side}_nonfinite"]
        assert isinstance(stored, dict)
        for name, count in counts.items():
            stored[name] = int(stored[name]) + count

    nonfinite_mismatches = (
        np.count_nonzero(np.isnan(left_selected) != np.isnan(right_selected))
        + np.count_nonzero(np.isposinf(left_selected) != np.isposinf(right_selected))
        + np.count_nonzero(np.isneginf(left_selected) != np.isneginf(right_selected))
    )
    accumulator["nonfinite_class_mismatches"] = (
        int(accumulator["nonfinite_class_mismatches"])
        + int(nonfinite_mismatches)
    )

    raw_mismatch = left.view(unsigned_dtype) != right.view(unsigned_dtype)
    selected_mismatch = raw_mismatch & selected
    if np.any(selected_mismatch):
        accumulator["exact_equal"] = False
        if accumulator["first_difference"] is None:
            flat_index = int(np.flatnonzero(selected_mismatch)[0])
            index = tuple(int(value) for value in np.unravel_index(flat_index, left.shape))
            report = describe_index(index)
            report["left"] = _json_value(left[index])
            report["right"] = _json_value(right[index])
            accumulator["first_difference"] = report

    left_finite = np.isfinite(left)
    right_finite = np.isfinite(right)
    finite_pair = selected & left_finite & right_finite
    finite_count = int(np.count_nonzero(finite_pair))
    accumulator["finite_pairs"] = int(accumulator["finite_pairs"]) + finite_count
    if finite_count != selected_count:
        accumulator["within_tolerance"] = False
    if finite_count == 0:
        return

    left_values = np.asarray(left[finite_pair], dtype=np.float64)
    right_values = np.asarray(right[finite_pair], dtype=np.float64)
    with np.errstate(over="ignore", invalid="ignore"):
        difference = np.abs(right_values - left_values)
        tolerance = atol + rtol * np.abs(left_values)
    if np.any(difference > tolerance):
        accumulator["within_tolerance"] = False

    block_l1 = float(np.sum(difference, dtype=np.longdouble))
    reference_l1 = float(np.sum(np.abs(left_values), dtype=np.longdouble))
    accumulator["l1_sum"] = float(accumulator["l1_sum"]) + block_l1
    accumulator["reference_l1_sum"] = (
        float(accumulator["reference_l1_sum"]) + reference_l1
    )
    local_max_position = int(np.argmax(difference))
    local_max = float(difference[local_max_position])
    if local_max > float(accumulator["max_abs"]) or (
        accumulator["max_abs_location"] is None and local_max == 0.0
    ):
        finite_flat_indices = np.flatnonzero(finite_pair)
        full_flat_index = int(finite_flat_indices[local_max_position])
        index = tuple(
            int(value) for value in np.unravel_index(full_flat_index, left.shape)
        )
        accumulator["max_abs"] = local_max
        accumulator["max_abs_location"] = describe_index(index)


def _finalize(accumulator: dict[str, object]) -> dict[str, object]:
    reference = float(accumulator["reference_l1_sum"])
    difference = float(accumulator["l1_sum"])
    if reference > 0.0:
        relative = difference / reference
    elif difference == 0.0:
        relative = 0.0
    else:
        relative = math.inf
    result = {
        **accumulator,
        "relative_l1": _json_value(relative),
        "match": bool(accumulator["within_tolerance"]),
    }
    # Subtracting two finite values can overflow, and an L1 sum can exceed the
    # representable Real range.  Keep CLI JSON standards-compliant in that case.
    for name in ("max_abs", "l1_sum", "reference_l1_sum"):
        result[name] = _json_value(float(result[name]))
    return result


def _logical_report(metadata: RestartMetadata, gid: int) -> dict[str, int]:
    location = metadata.locations[gid]
    return {
        "level": location.level,
        "lx1": location.lx1,
        "lx2": location.lx2,
        "lx3": location.lx3,
    }


def _active_spatial_mask(layout: RestartFieldLayout) -> np.ndarray:
    nout1, nout2, nout3 = layout.stored_shape
    starts1, starts2, starts3 = layout.active_starts
    stops1, stops2, stops3 = layout.active_stops
    mask = np.zeros((nout3, nout2, nout1), dtype=bool)
    mask[starts3:stops3, starts2:stops2, starts1:stops1] = True
    return mask


def _active_face_mask(
    layout: RestartFieldLayout, component: int
) -> np.ndarray:
    nout1, nout2, nout3 = layout.stored_shape
    starts1, starts2, starts3 = layout.active_starts
    stops1, stops2, stops3 = layout.active_stops
    shape = [nout3, nout2, nout1]
    shape[2 - component] += 1
    mask = np.zeros(tuple(shape), dtype=bool)
    slices = [
        slice(starts3, stops3),
        slice(starts2, stops2),
        slice(starts1, stops1),
    ]
    slices[2 - component] = slice(
        (starts3, starts2, starts1)[2 - component],
        (stops3, stops2, stops1)[2 - component] + 1,
    )
    mask[tuple(slices)] = True
    return mask


def _read_exact(stream: BinaryIO, size: int, description: str) -> bytes:
    pieces: list[bytes] = []
    remaining = size
    while remaining:
        piece = stream.read(remaining)
        if not piece:
            break
        pieces.append(piece)
        remaining -= len(piece)
    if remaining:
        raise RuntimeError(
            f"restart ended early while reading {description}: expected {size} bytes, "
            f"got {size - remaining}"
        )
    return pieces[0] if len(pieces) == 1 else b"".join(pieces)


def _validate_file_layout(
    stream: BinaryIO,
    path: Path,
    signature: FileSignature,
    metadata: RestartMetadata,
    layout: RestartFieldLayout,
) -> int:
    # The supported MHD/optional-ADM layout has no Step-3 object state.  Therefore the
    # data-size record immediately follows metadata and the field slots consume the EOF.
    size_record_offset = metadata.metadata_end
    field_start = size_record_offset + struct.calcsize("Q")
    expected_size = field_start + metadata.nmb_total * layout.data_size
    if signature.size != expected_size:
        raise ValueError(
            f"restart size disagrees with derived field layout for {path}: "
            f"expected {expected_size}, found {signature.size}"
        )
    stream.seek(size_record_offset)
    byte_prefix = "<" if metadata.byte_order == "little" else ">"
    stored_data_size = struct.unpack(
        f"{byte_prefix}Q",
        _read_exact(stream, struct.calcsize("Q"), "per-MeshBlock data_size"),
    )[0]
    if stored_data_size != layout.data_size:
        raise ValueError(
            f"restart data_size mismatch for {path}: stored {stored_data_size}, "
            f"derived {layout.data_size}"
        )
    return field_start


def compare_restart_fields(
    left: Path | str,
    right: Path | str,
    *,
    rtol: float = 0.0,
    atol: float = 0.0,
) -> dict[str, object]:
    """Stream and compare supported physical fields from two stable restarts."""

    if not (math.isfinite(rtol) and math.isfinite(atol) and rtol >= 0 and atol >= 0):
        raise ValueError("rtol and atol must be finite and non-negative")
    left_path = Path(left).resolve(strict=True)
    right_path = Path(right).resolve(strict=True)
    left_signature = _file_signature(left_path)
    right_signature = _file_signature(right_path)
    left_metadata = read_restart_metadata(left_path)
    right_metadata = read_restart_metadata(right_path)
    _assert_path_unchanged(left_path, left_signature, "while metadata was read")
    _assert_path_unchanged(right_path, right_signature, "while metadata was read")

    if _topology_signature(left_metadata) != _topology_signature(right_metadata):
        raise ValueError("restart topology or field ABI differs")
    left_layout = _derive_layout(left_metadata)
    right_layout = _derive_layout(right_metadata)
    if left_layout != right_layout:
        raise ValueError("restart-derived physical field layouts differ")
    layout = left_layout

    mhd_active = _empty_accumulator()
    mhd_ghost = _empty_accumulator()
    face_active = _empty_accumulator()
    face_components = [_empty_accumulator() for _ in range(3)]
    adm_active = _empty_accumulator()
    adm_ghost = _empty_accumulator()
    active_mask = _active_spatial_mask(layout)
    ghost_mask = ~active_mask
    face_masks = tuple(_active_face_mask(layout, component) for component in range(3))
    mhd_names = _mhd_names(layout)

    nout1, nout2, nout3 = layout.stored_shape
    with left_path.open("rb") as left_stream, right_path.open("rb") as right_stream:
        _assert_stream_unchanged(
            left_stream, left_path, left_signature, "before field streaming"
        )
        _assert_stream_unchanged(
            right_stream, right_path, right_signature, "before field streaming"
        )
        left_field_start = _validate_file_layout(
            left_stream, left_path, left_signature, left_metadata, layout
        )
        right_field_start = _validate_file_layout(
            right_stream, right_path, right_signature, right_metadata, layout
        )
        left_stream.seek(left_field_start)
        right_stream.seek(right_field_start)

        for gid in range(left_metadata.nmb_total):
            left_block = _read_exact(
                left_stream, layout.data_size, f"left MeshBlock GID {gid}"
            )
            right_block = _read_exact(
                right_stream, layout.data_size, f"right MeshBlock GID {gid}"
            )
            logical = _logical_report(left_metadata, gid)

            left_mhd = np.frombuffer(
                left_block,
                dtype=layout.real_dtype,
                count=layout.mhd_elements,
                offset=layout.mhd_offset,
            ).reshape((layout.nmhd, nout3, nout2, nout1))
            right_mhd = np.frombuffer(
                right_block,
                dtype=layout.real_dtype,
                count=layout.mhd_elements,
                offset=layout.mhd_offset,
            ).reshape(left_mhd.shape)

            def describe_mhd(index: tuple[int, ...]) -> dict[str, object]:
                variable, k, j, i = index
                return {
                    "gid": gid,
                    "logical": logical,
                    "variable": variable,
                    "variable_name": mhd_names[variable],
                    "k": k,
                    "j": j,
                    "i": i,
                }

            for accumulator, mask in (
                (mhd_active, active_mask[None, ...]),
                (mhd_ghost, ghost_mask[None, ...]),
            ):
                _accumulate_region(
                    accumulator,
                    left_mhd,
                    right_mhd,
                    mask,
                    describe_mhd,
                    unsigned_dtype=layout.unsigned_dtype,
                    rtol=rtol,
                    atol=atol,
                )

            face_shapes = (
                (nout3, nout2, nout1 + 1),
                (nout3, nout2 + 1, nout1),
                (nout3 + 1, nout2, nout1),
            )
            for component, (offset, elements, shape, mask) in enumerate(zip(
                layout.face_offsets, layout.face_elements, face_shapes, face_masks
            )):
                left_face = np.frombuffer(
                    left_block,
                    dtype=layout.real_dtype,
                    count=elements,
                    offset=offset,
                ).reshape(shape)
                right_face = np.frombuffer(
                    right_block,
                    dtype=layout.real_dtype,
                    count=elements,
                    offset=offset,
                ).reshape(shape)

                def describe_face(
                    index: tuple[int, ...], component: int = component
                ) -> dict[str, object]:
                    k, j, i = index
                    return {
                        "gid": gid,
                        "logical": logical,
                        "component": component + 1,
                        "component_name": f"x{component + 1}f",
                        "k": k,
                        "j": j,
                        "i": i,
                    }

                for accumulator in (face_active, face_components[component]):
                    _accumulate_region(
                        accumulator,
                        left_face,
                        right_face,
                        mask,
                        describe_face,
                        unsigned_dtype=layout.unsigned_dtype,
                        rtol=rtol,
                        atol=atol,
                    )

            if layout.nadm:
                left_adm = np.frombuffer(
                    left_block,
                    dtype=layout.real_dtype,
                    count=layout.adm_elements,
                    offset=layout.adm_offset,
                ).reshape((layout.nadm, nout3, nout2, nout1))
                right_adm = np.frombuffer(
                    right_block,
                    dtype=layout.real_dtype,
                    count=layout.adm_elements,
                    offset=layout.adm_offset,
                ).reshape(left_adm.shape)

                def describe_adm(index: tuple[int, ...]) -> dict[str, object]:
                    variable, k, j, i = index
                    return {
                        "gid": gid,
                        "logical": logical,
                        "variable": variable,
                        "variable_name": ADM_NAMES[variable],
                        "k": k,
                        "j": j,
                        "i": i,
                    }

                for accumulator, mask in (
                    (adm_active, active_mask[None, ...]),
                    (adm_ghost, ghost_mask[None, ...]),
                ):
                    _accumulate_region(
                        accumulator,
                        left_adm,
                        right_adm,
                        mask,
                        describe_adm,
                        unsigned_dtype=layout.unsigned_dtype,
                        rtol=rtol,
                        atol=atol,
                    )

        _assert_stream_unchanged(
            left_stream, left_path, left_signature, "after field streaming"
        )
        _assert_stream_unchanged(
            right_stream, right_path, right_signature, "after field streaming"
        )

    _assert_path_unchanged(left_path, left_signature, "during comparison")
    _assert_path_unchanged(right_path, right_signature, "during comparison")

    mhd_active_result = _finalize(mhd_active)
    mhd_ghost_result = _finalize(mhd_ghost)
    face_result = _finalize(face_active)
    face_component_results = {
        f"x{component + 1}f": _finalize(accumulator)
        for component, accumulator in enumerate(face_components)
    }
    same_endpoint = bool(
        left_metadata.time == right_metadata.time
        and left_metadata.cycle == right_metadata.cycle
    )
    included_results = [mhd_active_result, mhd_ghost_result, face_result]
    adm_result: dict[str, object]
    if layout.nadm:
        adm_active_result = _finalize(adm_active)
        adm_ghost_result = _finalize(adm_ghost)
        included_results.extend((adm_active_result, adm_ghost_result))
        adm_result = {
            "present": True,
            "active": adm_active_result,
            "ghost": adm_ghost_result,
        }
    else:
        adm_result = {"present": False}

    authoritative_fields_match = bool(
        mhd_active_result["match"] and face_result["match"]
    )
    fields_match = all(bool(result["match"]) for result in included_results)
    return {
        "kind": "athenak_restart_field_comparison",
        "left": str(left_path),
        "right": str(right_path),
        "match": bool(same_endpoint and fields_match),
        # Active conserved cells and CT faces are the authoritative restart state.
        # Ghost cells are reconstructed during initialization, while a prescribed
        # dynamic ADM metric is deterministically rebuilt from time.  Keep reporting
        # those caches without conflating them with an evolved-state mismatch.
        "authoritative_match": bool(same_endpoint and authoritative_fields_match),
        "all_stored_fields_match": bool(same_endpoint and fields_match),
        "same_endpoint": same_endpoint,
        "endpoint": {
            "left": {
                "time": left_metadata.time,
                "cycle": left_metadata.cycle,
                "last_dt": left_metadata.last_dt,
            },
            "right": {
                "time": right_metadata.time,
                "cycle": right_metadata.cycle,
                "last_dt": right_metadata.last_dt,
            },
        },
        "topology": {
            "match": True,
            "nmb_total": left_metadata.nmb_total,
            "root_level": left_metadata.root_level,
        },
        "layout": {
            "real_bytes": layout.real_bytes,
            "byte_order": left_metadata.byte_order,
            "nghost": layout.nghost,
            "meshblock_shape": list(layout.block_shape),
            "stored_shape": list(layout.stored_shape),
            "nmhd_base": layout.nmhd_base,
            "nscalars": layout.nscalars,
            "nmhd_total": layout.nmhd,
            "nadm": layout.nadm,
            "data_size": layout.data_size,
            "streaming_unit": "one MeshBlock per file",
            "raw_block_pair_bytes": 2 * layout.data_size,
        },
        "tolerance": {"rtol": rtol, "atol": atol},
        "mhd_u0": {
            "active": mhd_active_result,
            "ghost": mhd_ghost_result,
        },
        "face_b": {
            "active_faces": face_result,
            "components": face_component_results,
            "ghost_faces_compared": False,
        },
        "adm": adm_result,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument("--rtol", type=float, default=0.0)
    parser.add_argument("--atol", type=float, default=0.0)
    parser.add_argument(
        "--compact", action="store_true", help="emit compact rather than indented JSON"
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        result = compare_restart_fields(
            args.left, args.right, rtol=args.rtol, atol=args.atol
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            result,
            indent=None if args.compact else 2,
            sort_keys=True,
            allow_nan=False,
        )
    )
    return 0 if result["match"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
