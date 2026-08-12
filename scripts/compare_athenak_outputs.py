#!/usr/bin/env python3
"""Compare two *closed* AthenaK binary outputs or history endpoints.

Binary fields are loaded one variable at a time through
``vis/python/bin_convert.py``.  MeshBlocks are paired by their canonical
logical location ``(level, lx3, lx2, lx1)`` rather than by file/GID order, so
MPI load-balancing differences do not produce false mismatches.

For each field, ``L1`` is ``sum(abs(right-left))`` over finite pairs and
``relative_L1`` divides that value by ``sum(abs(left))``.  ``max_relative``
uses the left value as its pointwise denominator and is infinite for a
nonzero difference from a zero reference.  Non-finite values are counted
separately and make the tolerance-based verdict fail.

This tool deliberately does not try to lock or copy its inputs.  Only invoke
it after both producers have closed the files.  It checks size and mtime before
and after the comparison and rejects files that changed while being read.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable, Sequence

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
VIS_PYTHON = REPOSITORY_ROOT / "vis" / "python"
if str(VIS_PYTHON) not in sys.path:
    sys.path.insert(0, str(VIS_PYTHON))

import athena_read  # noqa: E402
import bin_convert  # noqa: E402


TOPOLOGY_SCALARS = (
    "Nx1",
    "Nx2",
    "Nx3",
    "nx1_mb",
    "nx2_mb",
    "nx3_mb",
    "nx1_out_mb",
    "nx2_out_mb",
    "nx3_out_mb",
    "n_mbs",
)


def _file_signature(path: Path) -> tuple[int, int]:
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns


def _assert_unchanged(path: Path, signature: tuple[int, int]) -> None:
    if _file_signature(path) != signature:
        raise RuntimeError(
            f"{path} changed while it was being read; compare only closed files"
        )


def canonical_order(dump: dict[str, Any]) -> np.ndarray:
    """Return MeshBlock indices sorted by (level, lx3, lx2, lx1)."""

    logical = np.asarray(dump["mb_logical"])
    if logical.shape != (int(dump["n_mbs"]), 4):
        raise ValueError(
            f"invalid mb_logical shape {logical.shape}; expected "
            f"({dump['n_mbs']}, 4)"
        )
    if len({tuple(int(value) for value in row) for row in logical}) != len(logical):
        raise ValueError("duplicate MeshBlock logical locations in binary output")
    # np.lexsort uses the last key as primary.
    return np.lexsort(
        (logical[:, 0], logical[:, 1], logical[:, 2], logical[:, 3])
    )


def _nonfinite_counts(values: np.ndarray) -> dict[str, int]:
    return {
        "nan": int(np.count_nonzero(np.isnan(values))),
        "positive_infinity": int(np.count_nonzero(np.isposinf(values))),
        "negative_infinity": int(np.count_nonzero(np.isneginf(values))),
    }


def _empty_accumulator() -> dict[str, Any]:
    return {
        "elements": 0,
        "finite_pairs": 0,
        "array_equal": True,
        "allclose": True,
        "max_abs": 0.0,
        "max_abs_location": None,
        "l1_sum": 0.0,
        "reference_l1_sum": 0.0,
        "max_relative": 0.0,
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
    }


def _accumulate_arrays(
    accumulator: dict[str, Any],
    left: np.ndarray,
    right: np.ndarray,
    *,
    rtol: float,
    atol: float,
    logical_location: Sequence[int] | None = None,
    block_geometry: Sequence[float] | None = None,
    block_indices: Sequence[int] | None = None,
    configured_block_shape: Sequence[int] | None = None,
) -> None:
    if left.shape != right.shape:
        raise ValueError(f"array shape mismatch: {left.shape} != {right.shape}")

    accumulator["elements"] += int(left.size)
    accumulator["array_equal"] &= bool(np.array_equal(left, right))
    accumulator["allclose"] &= bool(
        np.allclose(left, right, rtol=rtol, atol=atol, equal_nan=False)
    )

    for side, values in (("left", left), ("right", right)):
        counts = _nonfinite_counts(values)
        for name, count in counts.items():
            accumulator[f"{side}_nonfinite"][name] += count

    left_finite = np.isfinite(left)
    right_finite = np.isfinite(right)
    accumulator["nonfinite_class_mismatches"] += int(
        np.count_nonzero(
            np.isnan(left) != np.isnan(right)
        )
        + np.count_nonzero(np.isposinf(left) != np.isposinf(right))
        + np.count_nonzero(np.isneginf(left) != np.isneginf(right))
    )
    finite_pair = left_finite & right_finite
    finite_count = int(np.count_nonzero(finite_pair))
    accumulator["finite_pairs"] += finite_count
    if finite_count == 0:
        return

    # Boolean indexing bounds all temporaries to one MeshBlock in binary mode.
    left_finite_values = np.asarray(left[finite_pair], dtype=np.float64)
    right_finite_values = np.asarray(right[finite_pair], dtype=np.float64)
    difference = np.abs(right_finite_values - left_finite_values)
    reference_abs = np.abs(left_finite_values)
    block_max_flat_index = int(np.argmax(difference))
    block_max_abs = float(difference[block_max_flat_index])
    # Retain the first canonical occurrence when the maximum is tied.  This
    # makes the reported location deterministic across MPI/GID orderings.
    if block_max_abs > accumulator["max_abs"] or (
        accumulator["max_abs_location"] is None and block_max_abs == 0.0
    ):
        finite_flat_indices = np.flatnonzero(finite_pair)
        full_flat_index = int(finite_flat_indices[block_max_flat_index])
        local_index = np.unravel_index(full_flat_index, left.shape)
        accumulator["max_abs"] = block_max_abs
        location_report = {
            "logical": (
                [int(value) for value in logical_location]
                if logical_location is not None
                else None
            ),
            "local_kji": [int(value) for value in local_index],
            "left": float(left[local_index]),
            "right": float(right[local_index]),
            "signed_difference": float(right[local_index] - left[local_index]),
        }
        if (
            block_geometry is not None
            and block_indices is not None
            and configured_block_shape is not None
            and len(local_index) == 3
        ):
            k_index, j_index, i_index = local_index
            x1min, x1max, x2min, x2max, x3min, x3max = (
                float(value) for value in block_geometry
            )
            nx1, nx2, nx3 = (int(value) for value in configured_block_shape)
            i_index += int(block_indices[0])
            j_index += int(block_indices[2])
            k_index += int(block_indices[4])
            location_report["cell_center_xyz"] = [
                x1min + (i_index + 0.5) * (x1max - x1min) / nx1,
                x2min + (j_index + 0.5) * (x2max - x2min) / nx2,
                x3min + (k_index + 0.5) * (x3max - x3min) / nx3,
            ]
        accumulator["max_abs_location"] = location_report
    accumulator["l1_sum"] += float(np.sum(difference, dtype=np.float64))
    accumulator["reference_l1_sum"] += float(
        np.sum(reference_abs, dtype=np.float64)
    )

    nonzero_reference = reference_abs != 0.0
    if np.any(nonzero_reference):
        accumulator["max_relative"] = max(
            accumulator["max_relative"],
            float(
                np.max(
                    difference[nonzero_reference] / reference_abs[nonzero_reference]
                )
            ),
        )
    if np.any((~nonzero_reference) & (difference != 0.0)):
        accumulator["max_relative"] = math.inf


def _finish_accumulator(accumulator: dict[str, Any]) -> dict[str, Any]:
    result = dict(accumulator)
    count = result["finite_pairs"]
    result["l1_mean"] = result["l1_sum"] / count if count else None
    denominator = result.pop("reference_l1_sum")
    if denominator != 0.0:
        result["relative_l1"] = result["l1_sum"] / denominator
    elif result["l1_sum"] == 0.0:
        result["relative_l1"] = 0.0
    else:
        result["relative_l1"] = math.inf
    result["nonfinite_total"] = sum(result["left_nonfinite"].values()) + sum(
        result["right_nonfinite"].values()
    )
    result["match"] = bool(result["allclose"] and result["nonfinite_total"] == 0)
    return result


def _array_summary(
    left: np.ndarray, right: np.ndarray, *, rtol: float, atol: float
) -> dict[str, Any]:
    if left.shape != right.shape:
        return {
            "shape_equal": False,
            "left_shape": list(left.shape),
            "right_shape": list(right.shape),
            "array_equal": False,
            "allclose": False,
            "match": False,
        }
    accumulator = _empty_accumulator()
    _accumulate_arrays(accumulator, left, right, rtol=rtol, atol=atol)
    result = _finish_accumulator(accumulator)
    result["shape_equal"] = True
    result["left_shape"] = list(left.shape)
    result["right_shape"] = list(right.shape)
    return result


def _compare_topology(
    left: dict[str, Any],
    right: dict[str, Any],
    left_order: np.ndarray,
    right_order: np.ndarray,
) -> dict[str, Any]:
    scalars = {
        name: {
            "left": int(left[name]),
            "right": int(right[name]),
            "equal": int(left[name]) == int(right[name]),
        }
        for name in TOPOLOGY_SCALARS
    }
    logical = _array_summary(
        np.asarray(left["mb_logical"])[left_order],
        np.asarray(right["mb_logical"])[right_order],
        rtol=0.0,
        atol=0.0,
    )
    indices = _array_summary(
        np.asarray(left["mb_index"])[left_order],
        np.asarray(right["mb_index"])[right_order],
        rtol=0.0,
        atol=0.0,
    )
    geometry = _array_summary(
        np.asarray(left["mb_geometry"])[left_order],
        np.asarray(right["mb_geometry"])[right_order],
        rtol=0.0,
        atol=0.0,
    )
    return {
        "scalars": scalars,
        "logical_locations": logical,
        "output_indices": indices,
        "geometry": geometry,
        "match": bool(
            all(item["equal"] for item in scalars.values())
            and logical["match"]
            and indices["match"]
            and geometry["match"]
        ),
    }


def _compare_loaded_field(
    left: dict[str, Any],
    right: dict[str, Any],
    variable: str,
    left_order: np.ndarray,
    right_order: np.ndarray,
    *,
    rtol: float,
    atol: float,
) -> dict[str, Any]:
    left_blocks = left["mb_data"][variable]
    right_blocks = right["mb_data"][variable]
    if len(left_blocks) != len(right_blocks):
        return {
            "match": False,
            "error": (
                f"MeshBlock count mismatch: {len(left_blocks)} != "
                f"{len(right_blocks)}"
            ),
        }

    accumulator = _empty_accumulator()
    for left_gid, right_gid in zip(left_order, right_order):
        left_values = np.asarray(left_blocks[int(left_gid)])
        right_values = np.asarray(right_blocks[int(right_gid)])
        logical_location = np.asarray(left["mb_logical"])[int(left_gid)]
        block_geometry = np.asarray(left["mb_geometry"])[int(left_gid)]
        block_indices = np.asarray(left["mb_index"])[int(left_gid)]
        try:
            _accumulate_arrays(
                accumulator,
                left_values,
                right_values,
                rtol=rtol,
                atol=atol,
                logical_location=logical_location,
                block_geometry=block_geometry,
                block_indices=block_indices,
                configured_block_shape=(
                    left["nx1_mb"],
                    left["nx2_mb"],
                    left["nx3_mb"],
                ),
            )
        except ValueError as exc:
            return {
                "match": False,
                "error": str(exc),
                "left_logical": np.asarray(left["mb_logical"])[int(left_gid)].tolist(),
                "right_logical": np.asarray(right["mb_logical"])[int(right_gid)].tolist(),
            }
    return _finish_accumulator(accumulator)


def _validate_requested_names(names: Sequence[str] | None) -> list[str] | None:
    if names is None:
        return None
    result: list[str] = []
    for name in names:
        for item in name.split(","):
            item = item.strip()
            if item and item not in result:
                result.append(item)
    if not result:
        raise ValueError("at least one non-empty variable/column name is required")
    return result


def compare_binary_files(
    left_path: Path,
    right_path: Path,
    *,
    variables: Sequence[str] | None = None,
    rtol: float = 0.0,
    atol: float = 0.0,
    time_atol: float = 0.0,
) -> dict[str, Any]:
    """Compare two AthenaK binary files while holding one field at a time."""

    left_path = left_path.resolve(strict=True)
    right_path = right_path.resolve(strict=True)
    left_signature = _file_signature(left_path)
    right_signature = _file_signature(right_path)
    requested = _validate_requested_names(variables)

    left_metadata = bin_convert.read_binary(str(left_path), variables=())
    right_metadata = bin_convert.read_binary(str(right_path), variables=())
    left_order = canonical_order(left_metadata)
    right_order = canonical_order(right_metadata)
    topology = _compare_topology(
        left_metadata, right_metadata, left_order, right_order
    )

    left_names = list(left_metadata["var_names"])
    right_names = list(right_metadata["var_names"])
    names_equal = left_names == right_names
    if requested is None:
        selected = [name for name in left_names if name in right_names]
        required_names_present = set(left_names) == set(right_names)
    else:
        selected = [
            name for name in requested if name in left_names and name in right_names
        ]
        required_names_present = len(selected) == len(requested)

    fields: dict[str, Any] = {}
    if topology["match"]:
        for variable in selected:
            left_field = bin_convert.read_binary(
                str(left_path), variables=(variable,)
            )
            right_field = bin_convert.read_binary(
                str(right_path), variables=(variable,)
            )
            fields[variable] = _compare_loaded_field(
                left_field,
                right_field,
                variable,
                left_order,
                right_order,
                rtol=rtol,
                atol=atol,
            )
            del left_field, right_field
            gc.collect()

    _assert_unchanged(left_path, left_signature)
    _assert_unchanged(right_path, right_signature)

    cycle_equal = int(left_metadata["cycle"]) == int(right_metadata["cycle"])
    time_equal = math.isclose(
        float(left_metadata["time"]),
        float(right_metadata["time"]),
        rel_tol=0.0,
        abs_tol=time_atol,
    )
    fields_match = bool(selected) and all(
        item.get("match", False) for item in fields.values()
    )
    result = {
        "kind": "athenak_binary_comparison",
        "left": str(left_path),
        "right": str(right_path),
        "cycle": {
            "left": int(left_metadata["cycle"]),
            "right": int(right_metadata["cycle"]),
            "equal": cycle_equal,
        },
        "time": {
            "left": float(left_metadata["time"]),
            "right": float(right_metadata["time"]),
            "absolute_difference": abs(
                float(right_metadata["time"]) - float(left_metadata["time"])
            ),
            "atol": time_atol,
            "equal": time_equal,
        },
        "variables": {
            "left": left_names,
            "right": right_names,
            "same_ordered_list": names_equal,
            "missing_from_left": sorted(set(right_names) - set(left_names)),
            "missing_from_right": sorted(set(left_names) - set(right_names)),
            "requested": requested,
            "compared": selected,
            "required_names_present": required_names_present,
        },
        "topology": topology,
        "fields": fields,
    }
    result["match"] = bool(
        cycle_equal
        and time_equal
        and topology["match"]
        and required_names_present
        and (names_equal if requested is None else True)
        and fields_match
    )
    result["array_equal"] = bool(
        cycle_equal
        and float(left_metadata["time"]) == float(right_metadata["time"])
        and topology["match"]
        and required_names_present
        and (names_equal if requested is None else True)
        and bool(selected)
        and all(item.get("array_equal", False) for item in fields.values())
    )
    return result


def _latest_overlapping_rows(
    left_time: np.ndarray, right_time: np.ndarray, time_atol: float
) -> tuple[int, int] | None:
    """Find the latest time represented in both monotonic histories."""

    left_index = len(left_time) - 1
    right_index = len(right_time) - 1
    while left_index >= 0 and right_index >= 0:
        left_value = float(left_time[left_index])
        right_value = float(right_time[right_index])
        if math.isclose(left_value, right_value, rel_tol=0.0, abs_tol=time_atol):
            return left_index, right_index
        if left_value > right_value:
            left_index -= 1
        else:
            right_index -= 1
    return None


def compare_history_endpoints(
    left_path: Path,
    right_path: Path,
    *,
    columns: Sequence[str] | None = None,
    rtol: float = 0.0,
    atol: float = 0.0,
    time_atol: float = 1.0e-12,
) -> dict[str, Any]:
    """Compare the latest row whose time occurs in both history segments."""

    left_path = left_path.resolve(strict=True)
    right_path = right_path.resolve(strict=True)
    left_signature = _file_signature(left_path)
    right_signature = _file_signature(right_path)
    requested = _validate_requested_names(columns)
    left = athena_read.hst(str(left_path), raw=False)
    right = athena_read.hst(str(right_path), raw=False)
    if "time" not in left or "time" not in right:
        raise ValueError("both histories must contain a time column")
    if len(left["time"]) == 0 or len(right["time"]) == 0:
        raise ValueError("both histories must contain at least one data row")

    left_names = list(left)
    right_names = list(right)
    if requested is None:
        selected = [name for name in left_names if name in right_names]
        required_names_present = set(left_names) == set(right_names)
    else:
        selected = [
            name for name in requested if name in left_names and name in right_names
        ]
        required_names_present = len(selected) == len(requested)

    pair = _latest_overlapping_rows(
        np.asarray(left["time"]), np.asarray(right["time"]), time_atol
    )
    column_results: dict[str, Any] = {}
    overlap: dict[str, Any]
    if pair is None:
        overlap = {"found": False, "time_atol": time_atol}
    else:
        left_index, right_index = pair
        overlap = {
            "found": True,
            "left_index": left_index,
            "right_index": right_index,
            "left_time": float(left["time"][left_index]),
            "right_time": float(right["time"][right_index]),
            "left_is_terminal": left_index == len(left["time"]) - 1,
            "right_is_terminal": right_index == len(right["time"]) - 1,
            "time_atol": time_atol,
        }
        for name in selected:
            left_value = np.asarray([left[name][left_index]], dtype=np.float64)
            right_value = np.asarray([right[name][right_index]], dtype=np.float64)
            column_atol = max(atol, time_atol) if name == "time" else atol
            summary = _array_summary(
                left_value, right_value, rtol=rtol, atol=column_atol
            )
            summary["left"] = float(left_value[0])
            summary["right"] = float(right_value[0])
            column_results[name] = summary

    _assert_unchanged(left_path, left_signature)
    _assert_unchanged(right_path, right_signature)
    columns_match = bool(selected) and all(
        item.get("match", False) for item in column_results.values()
    )
    result = {
        "kind": "athenak_history_endpoint_comparison",
        "left": str(left_path),
        "right": str(right_path),
        "left_rows": len(left["time"]),
        "right_rows": len(right["time"]),
        "columns": {
            "left": left_names,
            "right": right_names,
            "same_ordered_list": left_names == right_names,
            "missing_from_left": sorted(set(right_names) - set(left_names)),
            "missing_from_right": sorted(set(left_names) - set(right_names)),
            "requested": requested,
            "compared": selected,
            "required_names_present": required_names_present,
        },
        "overlap": overlap,
        "endpoint_columns": column_results,
    }
    result["match"] = bool(
        overlap["found"]
        and required_names_present
        and (left_names == right_names if requested is None else True)
        and columns_match
    )
    result["array_equal"] = bool(
        overlap["found"]
        and required_names_present
        and (left_names == right_names if requested is None else True)
        and bool(selected)
        and all(
            item.get("array_equal", False) for item in column_results.values()
        )
    )
    return result


def _format_number(value: Any) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        if math.isinf(value):
            return "inf" if value > 0.0 else "-inf"
        return f"{value:.17g}"
    return str(value)


def _format_metric(name: str, result: dict[str, Any]) -> str:
    if "error" in result:
        return f"  {name}: ERROR {result['error']}"
    line = (
        f"  {name}: array_equal={result['array_equal']} "
        f"allclose={result['allclose']} "
        f"max_abs={_format_number(result.get('max_abs'))} "
        f"L1={_format_number(result.get('l1_sum'))} "
        f"relative_L1={_format_number(result.get('relative_l1'))} "
        f"max_relative={_format_number(result.get('max_relative'))} "
        f"nonfinite=({result.get('left_nonfinite')}, "
        f"{result.get('right_nonfinite')})"
    )
    location = result.get("max_abs_location")
    if location is not None:
        line += (
            " max_abs_at="
            f"logical{location['logical']}/kji{location['local_kji']} "
            f"left={_format_number(location['left'])} "
            f"right={_format_number(location['right'])} "
            f"right-left={_format_number(location['signed_difference'])}"
        )
        if "cell_center_xyz" in location:
            line += f" xyz={location['cell_center_xyz']}"
    return line


def format_text(result: dict[str, Any]) -> str:
    verdict = "MATCH" if result["match"] else "DIFFER"
    lines = [
        f"{result['kind']}: {verdict}",
        f"left:  {result['left']}",
        f"right: {result['right']}",
    ]
    if result["kind"] == "athenak_binary_comparison":
        cycle = result["cycle"]
        time = result["time"]
        lines.extend(
            (
                f"cycle: {cycle['left']} vs {cycle['right']} equal={cycle['equal']}",
                "time: "
                f"{_format_number(time['left'])} vs {_format_number(time['right'])} "
                f"abs_diff={_format_number(time['absolute_difference'])} "
                f"equal={time['equal']}",
                f"topology: match={result['topology']['match']}",
                "variable lists: "
                f"same_order={result['variables']['same_ordered_list']} "
                f"missing_left={result['variables']['missing_from_left']} "
                f"missing_right={result['variables']['missing_from_right']}",
                "variables compared: "
                + (", ".join(result["variables"]["compared"]) or "(none)"),
            )
        )
        topology = result["topology"]
        unequal_scalars = [
            f"{name}={item['left']}/{item['right']}"
            for name, item in topology["scalars"].items()
            if not item["equal"]
        ]
        if unequal_scalars:
            lines.append("  topology scalar differences: " + ", ".join(unequal_scalars))
        for name, key in (
            ("logical_locations", "logical locations"),
            ("output_indices", "output indices"),
            ("geometry", "geometry"),
        ):
            summary = topology[name]
            if not summary["match"]:
                lines.append(_format_metric(key, summary))
        for name, summary in result["fields"].items():
            lines.append(_format_metric(name, summary))
    else:
        overlap = result["overlap"]
        if overlap["found"]:
            lines.append(
                "overlap: "
                f"left[{overlap['left_index']}]={_format_number(overlap['left_time'])}, "
                f"right[{overlap['right_index']}]={_format_number(overlap['right_time'])}"
            )
        else:
            lines.append("overlap: none")
        lines.append(
            "columns compared: "
            + (", ".join(result["columns"]["compared"]) or "(none)")
        )
        if not result["columns"]["same_ordered_list"]:
            lines.append(
                "column lists differ: "
                f"missing_left={result['columns']['missing_from_left']} "
                f"missing_right={result['columns']['missing_from_right']}"
            )
        for name, summary in result["endpoint_columns"].items():
            lines.append(_format_metric(name, summary))
    lines.append(f"array_equal: {result['array_equal']}")
    return "\n".join(lines)


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        if math.isnan(value):
            return "nan"
        return "inf" if value > 0.0 else "-inf"
    return value


def _nonnegative_float(text: str) -> float:
    value = float(text)
    if not math.isfinite(value) or value < 0.0:
        raise argparse.ArgumentTypeError("must be a finite non-negative number")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    binary = subparsers.add_parser(
        "bin", help="compare two closed AthenaK .bin outputs"
    )
    binary.add_argument("left", type=Path)
    binary.add_argument("right", type=Path)
    binary.add_argument(
        "--variables",
        nargs="+",
        help="variables to compare (space- or comma-separated; default: all)",
    )
    binary.add_argument("--rtol", type=_nonnegative_float, default=0.0)
    binary.add_argument("--atol", type=_nonnegative_float, default=0.0)
    binary.add_argument("--time-atol", type=_nonnegative_float, default=0.0)
    binary.add_argument("--json", action="store_true", help="emit JSON to stdout")

    history = subparsers.add_parser(
        "hst", help="compare the latest overlapping row of two closed .hst files"
    )
    history.add_argument("left", type=Path)
    history.add_argument("right", type=Path)
    history.add_argument(
        "--columns",
        nargs="+",
        help="columns to compare (space- or comma-separated; default: all)",
    )
    history.add_argument("--rtol", type=_nonnegative_float, default=0.0)
    history.add_argument("--atol", type=_nonnegative_float, default=0.0)
    history.add_argument("--time-atol", type=_nonnegative_float, default=1.0e-12)
    history.add_argument("--json", action="store_true", help="emit JSON to stdout")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "bin":
            result = compare_binary_files(
                args.left,
                args.right,
                variables=args.variables,
                rtol=args.rtol,
                atol=args.atol,
                time_atol=args.time_atol,
            )
        else:
            result = compare_history_endpoints(
                args.left,
                args.right,
                columns=args.columns,
                rtol=args.rtol,
                atol=args.atol,
                time_atol=args.time_atol,
            )
    except (OSError, RuntimeError, TypeError, ValueError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(_json_safe(result), indent=2, sort_keys=True, allow_nan=False))
    else:
        print(format_text(result))
    return 0 if result["match"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
