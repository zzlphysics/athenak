#!/usr/bin/env python3
"""Strictly validate and summarize one AthenaK FOFC spatial event-log cycle.

The ``# fofc_spatial_v1`` records are diagnostic comments embedded next to the
traditional integer event row.  This tool validates the fixed v1 schema and the
whole log as an ordered sequence of closed ``summary -> bins -> event`` groups.
It then reports the selected cycle as one conservative object: exactly one
summary and traditional row; ``summary.count == summary.nfofc == event.fofc``;
unique in-range histogram keys; and a histogram and all six marginals that sum
to the authoritative FOFC count.  A restart/IC prefix may be reported as
``unattributed`` only in the first group and only in the canonical
overflow/other/unknown bin.  Use ``--require-unattributed-zero`` when a clean
selected-cycle attribution is a qualification policy.

The input is opened read-only with ``O_NOFOLLOW``.  Device, inode, size, mtime,
and ctime must remain unchanged for the complete streaming read.  The tool does
not create, lock, truncate, or otherwise modify the event log.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import re
import stat
import sys
from typing import Any, Iterable, Iterator, Mapping, Sequence


EVENT_COLUMNS = (
    "cycle",
    "eos_dfloor",
    "eos_efloor",
    "eos_tfloor",
    "eos_vceil",
    "eos_fail",
    "c2p_it",
    "fofc",
    "cons_adjust",
    "mag_adjust",
    "c2p_calls",
    "fofc_tests",
)
LEVEL_BINS = 33
STAGE_LABELS = ("other", "1", "2", "3")
REASON_LABELS = (
    "unknown",
    "dmp_preflag",
    "scalar",
    "cons_density_floor",
    "cons_energy_floor",
    "prim_density_floor",
    "prim_temperature_floor",
    "rho_too_big",
    "rho_too_small",
    "nans_in_cons",
    "mag_too_big",
    "bracketing_failed",
    "no_solution",
    "invalid_geometry",
    "other_c2p",
)
REASON_INDEX = {reason: index for index, reason in enumerate(REASON_LABELS)}
R_CYL_EDGES = (2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
ABS_Z_EDGES = (0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
LAPSE_EDGES = (0.2, 0.4, 0.6, 0.8, 1.0)
CANONICAL_UNATTRIBUTED_KEY = (
    LEVEL_BINS - 1,
    0,
    REASON_INDEX["unknown"],
    len(R_CYL_EDGES),
    len(ABS_Z_EDGES),
    len(LAPSE_EDGES),
)
UINT64_MAX = (1 << 64) - 1
INT_MAX = (1 << 31) - 1
KEY_PATTERN = re.compile(r"[a-z][a-z0-9_]*\Z")
UNSIGNED_PATTERN = re.compile(r"[0-9]+\Z")

SCHEMA_KEYS = frozenset(
    {
        "kind",
        "level_bins",
        "stage_bins",
        "reason_bins",
        "r_cyl_edges",
        "abs_z_edges",
        "lapse_edges",
        "center1",
        "center2",
        "center3",
    }
)
SUMMARY_KEYS = frozenset({"kind", "cycle", "count", "nfofc", "unattributed"})
BIN_KEYS = frozenset(
    {
        "kind",
        "cycle",
        "level_bin",
        "stage_bin",
        "reason",
        "r_cyl_bin",
        "abs_z_bin",
        "lapse_bin",
        "count",
    }
)


class FofcTelemetryError(ValueError):
    """The event log does not satisfy the strict FOFC v1 contract."""


def _failure(source: str, line_number: int | None, message: str) -> None:
    location = source if line_number is None else f"{source}:{line_number}"
    raise FofcTelemetryError(f"{location}: {message}")


def _require_keys(
    record: Mapping[str, str], expected: frozenset[str], source: str,
    line_number: int,
) -> None:
    actual = set(record)
    if actual == expected:
        return
    missing = sorted(expected - actual)
    extra = sorted(actual - expected)
    _failure(
        source,
        line_number,
        f"record keys are not exact; missing={missing}, extra={extra}",
    )


def _unsigned(
    text: str, name: str, source: str, line_number: int,
    maximum: int = UINT64_MAX,
) -> int:
    if UNSIGNED_PATTERN.fullmatch(text) is None:
        _failure(source, line_number, f"{name} is not an unsigned decimal integer")
    canonical = text.lstrip("0") or "0"
    limit = str(maximum)
    if len(canonical) > len(limit) or (
        len(canonical) == len(limit) and canonical > limit
    ):
        _failure(source, line_number, f"{name} exceeds {maximum}")
    return int(canonical)


def _finite_float(text: str, name: str, source: str, line_number: int) -> float:
    try:
        value = float(text)
    except ValueError:
        _failure(source, line_number, f"{name} is not a floating-point number")
    if not math.isfinite(value):
        _failure(source, line_number, f"{name} is not finite")
    return value


def _parse_record(line: str, source: str, line_number: int) -> dict[str, str]:
    fields = line.split()
    if fields[:2] != ["#", "fofc_spatial_v1"]:
        _failure(source, line_number, "invalid FOFC telemetry prefix")
    record: dict[str, str] = {}
    for field in fields[2:]:
        if "=" not in field:
            _failure(source, line_number, f"telemetry token lacks '=': {field!r}")
        key, value = field.split("=", 1)
        if KEY_PATTERN.fullmatch(key) is None or not value:
            _failure(source, line_number, f"invalid telemetry token: {field!r}")
        if key in record:
            _failure(source, line_number, f"duplicate telemetry key {key!r}")
        record[key] = value
    if "kind" not in record:
        _failure(source, line_number, "telemetry record lacks kind")
    return record


def _parse_csv_floats(
    text: str, name: str, expected: Sequence[float], source: str,
    line_number: int,
) -> tuple[float, ...]:
    fields = text.split(",")
    values = tuple(
        _finite_float(field, name, source, line_number) for field in fields
    )
    if values != tuple(expected):
        _failure(
            source,
            line_number,
            f"{name} changed: expected {tuple(expected)}, found {values}",
        )
    return values


def _parse_schema(
    record: Mapping[str, str], source: str, line_number: int,
) -> dict[str, Any]:
    _require_keys(record, SCHEMA_KEYS, source, line_number)
    expected_literals = {
        "level_bins": "0..31,overflow",
        "stage_bins": ",".join(STAGE_LABELS),
        "reason_bins": ",".join(REASON_LABELS),
    }
    for name, expected in expected_literals.items():
        if record[name] != expected:
            _failure(
                source,
                line_number,
                f"{name} changed: expected {expected!r}, found {record[name]!r}",
            )
    r_edges = _parse_csv_floats(
        record["r_cyl_edges"], "r_cyl_edges", R_CYL_EDGES, source, line_number
    )
    z_edges = _parse_csv_floats(
        record["abs_z_edges"], "abs_z_edges", ABS_Z_EDGES, source, line_number
    )
    lapse_edges = _parse_csv_floats(
        record["lapse_edges"], "lapse_edges", LAPSE_EDGES, source, line_number
    )
    center = [
        _finite_float(record[f"center{axis}"], f"center{axis}", source, line_number)
        for axis in (1, 2, 3)
    ]
    return {
        "version": 1,
        "level_bins": [str(value) for value in range(32)] + ["overflow"],
        "stage_bins": list(STAGE_LABELS),
        "reason_bins": list(REASON_LABELS),
        "r_cyl_edges": list(r_edges),
        "abs_z_edges": list(z_edges),
        "lapse_edges": list(lapse_edges),
        "center": center,
    }


def _parse_summary(
    record: Mapping[str, str], source: str, line_number: int,
) -> dict[str, int]:
    _require_keys(record, SUMMARY_KEYS, source, line_number)
    return {
        "cycle": _unsigned(record["cycle"], "cycle", source, line_number, INT_MAX),
        "count": _unsigned(record["count"], "count", source, line_number),
        "nfofc": _unsigned(record["nfofc"], "nfofc", source, line_number),
        "unattributed": _unsigned(
            record["unattributed"], "unattributed", source, line_number
        ),
        "line": line_number,
    }


def _bounded_bin(
    record: Mapping[str, str], name: str, bins: int, source: str,
    line_number: int,
) -> int:
    value = _unsigned(record[name], name, source, line_number, INT_MAX)
    if value >= bins:
        _failure(source, line_number, f"{name}={value} is outside 0..{bins - 1}")
    return value


def _parse_bin(
    record: Mapping[str, str], source: str, line_number: int,
) -> dict[str, Any]:
    _require_keys(record, BIN_KEYS, source, line_number)
    reason = record["reason"]
    if reason not in REASON_INDEX:
        _failure(source, line_number, f"unknown reason bin {reason!r}")
    count = _unsigned(record["count"], "count", source, line_number)
    if count == 0:
        _failure(source, line_number, "zero-count bin must not be emitted")
    return {
        "cycle": _unsigned(record["cycle"], "cycle", source, line_number, INT_MAX),
        "level_bin": _bounded_bin(
            record, "level_bin", LEVEL_BINS, source, line_number
        ),
        "stage_bin": _bounded_bin(
            record, "stage_bin", len(STAGE_LABELS), source, line_number
        ),
        "reason": reason,
        "r_cyl_bin": _bounded_bin(
            record, "r_cyl_bin", len(R_CYL_EDGES) + 1, source, line_number
        ),
        "abs_z_bin": _bounded_bin(
            record, "abs_z_bin", len(ABS_Z_EDGES) + 1, source, line_number
        ),
        "lapse_bin": _bounded_bin(
            record, "lapse_bin", len(LAPSE_EDGES) + 1, source, line_number
        ),
        "count": count,
        "line": line_number,
    }


def _parse_event_row(
    fields: Sequence[str], source: str, line_number: int,
) -> dict[str, int]:
    if len(fields) != len(EVENT_COLUMNS):
        _failure(
            source,
            line_number,
            f"event row has {len(fields)} columns, expected {len(EVENT_COLUMNS)}",
        )
    values = []
    for index, (name, field) in enumerate(zip(EVENT_COLUMNS, fields)):
        maximum = INT_MAX if name in ("cycle", "c2p_it") else UINT64_MAX
        values.append(_unsigned(field, name, source, line_number, maximum))
        if index == 0 and values[-1] < 0:
            _failure(source, line_number, "event cycle is negative")
    return dict(zip(EVENT_COLUMNS, values))


def _signature(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _checked_path(path: Path) -> Path:
    return Path(os.path.abspath(os.path.expanduser(os.fspath(path))))


def _read_lines_stably(
    path: Path,
) -> tuple[Iterator[tuple[int, str]], dict[str, Any], Path]:
    """Return a streaming line iterator plus state finalized by its consumer."""

    checked = _checked_path(path)
    initial = checked.lstat()
    if stat.S_ISLNK(initial.st_mode) or not stat.S_ISREG(initial.st_mode):
        raise FofcTelemetryError(
            f"{checked}: input must be a regular, non-symlink file"
        )
    expected = _signature(initial)
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    nofollow = getattr(os, "O_NOFOLLOW", None)
    if nofollow is None:
        raise FofcTelemetryError("O_NOFOLLOW is required for strict read-only input")
    descriptor = os.open(checked, flags | nofollow)
    stream = os.fdopen(descriptor, "rb", closefd=True)
    opened = os.fstat(stream.fileno())
    if _signature(opened) != expected or not stat.S_ISREG(opened.st_mode):
        stream.close()
        raise FofcTelemetryError(f"{checked}: input changed between lstat and open")

    state: dict[str, Any] = {
        "path": checked,
        "expected": expected,
        "stream": stream,
        "digest": hashlib.sha256(),
        "bytes_read": 0,
        "ends_with_newline": False,
        "complete": False,
    }

    def lines() -> Iterator[tuple[int, str]]:
        try:
            for line_number, raw_line in enumerate(stream, start=1):
                state["digest"].update(raw_line)
                state["bytes_read"] += len(raw_line)
                state["ends_with_newline"] = raw_line.endswith(b"\n")
                try:
                    line = raw_line.decode("utf-8")
                except UnicodeDecodeError as exc:
                    _failure(str(checked), line_number, "event log is not UTF-8")
                    raise AssertionError("unreachable") from exc
                yield line_number, line.rstrip("\r\n")
            if state["bytes_read"] == 0 or not state["ends_with_newline"]:
                raise FofcTelemetryError(
                    f"{checked}: event log is not terminated by a newline"
                )
            descriptor_final = os.fstat(stream.fileno())
            if (
                not stat.S_ISREG(descriptor_final.st_mode)
                or _signature(descriptor_final) != expected
            ):
                raise FofcTelemetryError(
                    f"{checked}: opened input changed while being read"
                )
            state["descriptor_stable"] = True
            state["complete"] = True
        finally:
            if not stream.closed:
                stream.close()

    return lines(), state, checked


def _finish_stable_read(state: Mapping[str, Any]) -> dict[str, Any]:
    path = state["path"]
    if not state["complete"]:
        raise FofcTelemetryError(f"{path}: event log read did not reach EOF")
    if not state.get("descriptor_stable"):
        raise FofcTelemetryError(f"{path}: input descriptor was not stable")
    expected = state["expected"]
    try:
        final = path.lstat()
    except FileNotFoundError as exc:
        raise FofcTelemetryError(f"{path}: input disappeared while being read") from exc
    if stat.S_ISLNK(final.st_mode) or not stat.S_ISREG(final.st_mode):
        raise FofcTelemetryError(f"{path}: input type changed while being read")
    if _signature(final) != expected:
        raise FofcTelemetryError(f"{path}: input changed while being read")
    device, inode, size, mtime_ns, ctime_ns = expected
    return {
        "path": str(path),
        "device": device,
        "inode": inode,
        "size": size,
        "mtime_ns": mtime_ns,
        "ctime_ns": ctime_ns,
        "sha256": state["digest"].hexdigest(),
        "stable_during_read": True,
        "open_mode": "read_only_nofollow",
    }


def _bin_key(row: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        row["level_bin"],
        row["stage_bin"],
        REASON_INDEX[row["reason"]],
        row["r_cyl_bin"],
        row["abs_z_bin"],
        row["lapse_bin"],
    )


def _format_edge(value: float) -> str:
    return format(value, "g")


def _interval_labels(edges: Sequence[float], first_floor: float = 0.0) -> list[str]:
    labels = [f"[{_format_edge(first_floor)},{_format_edge(edges[0])})"]
    labels.extend(
        f"[{_format_edge(left)},{_format_edge(right)})"
        for left, right in zip(edges, edges[1:])
    )
    labels.append(f"[{_format_edge(edges[-1])},inf)")
    return labels


def _marginal(
    bins: Sequence[Mapping[str, Any]], field: str, labels: Sequence[str],
    total: int,
) -> dict[str, Any]:
    counts = Counter(row[field] for row in bins)
    weighted = Counter()
    for row in bins:
        weighted[row[field]] += row["count"]
    rows = []
    for index, label in enumerate(labels):
        count = weighted[index]
        if count == 0:
            continue
        rows.append(
            {
                "bin": index,
                "label": label,
                "count": count,
                "fraction": count / total if total else 0.0,
                "records": counts[index],
            }
        )
    marginal_sum = sum(weighted.values())
    if marginal_sum != total:
        raise AssertionError(f"internal {field} marginal conservation failure")
    return {"sum": marginal_sum, "conservative": True, "bins": rows}


def _reason_marginal(
    bins: Sequence[Mapping[str, Any]], total: int,
) -> dict[str, Any]:
    weighted = Counter()
    records = Counter()
    for row in bins:
        weighted[row["reason"]] += row["count"]
        records[row["reason"]] += 1
    rows = []
    for index, reason in enumerate(REASON_LABELS):
        count = weighted[reason]
        if count == 0:
            continue
        rows.append(
            {
                "bin": index,
                "label": reason,
                "count": count,
                "fraction": count / total if total else 0.0,
                "records": records[reason],
            }
        )
    marginal_sum = sum(weighted.values())
    if marginal_sum != total:
        raise AssertionError("internal reason marginal conservation failure")
    return {"sum": marginal_sum, "conservative": True, "bins": rows}


def _decorated_bin(row: Mapping[str, Any], total: int) -> dict[str, Any]:
    radius_labels = _interval_labels(R_CYL_EDGES)
    z_labels = _interval_labels(ABS_Z_EDGES)
    lapse_labels = _interval_labels(LAPSE_EDGES)
    radius_labels[-1] = "other(r_cyl>=64_or_nonfinite)"
    z_labels[-1] = "other(abs_z>=16_or_nonfinite)"
    lapse_labels[-1] = "other(lapse<0,lapse>=1,or_nonfinite)"
    level_bin = row["level_bin"]
    stage_bin = row["stage_bin"]
    return {
        "level_bin": level_bin,
        "level": str(level_bin) if level_bin < 32 else "overflow",
        "stage_bin": stage_bin,
        "stage": STAGE_LABELS[stage_bin],
        "reason": row["reason"],
        "r_cyl_bin": row["r_cyl_bin"],
        "r_cyl": radius_labels[row["r_cyl_bin"]],
        "abs_z_bin": row["abs_z_bin"],
        "abs_z": z_labels[row["abs_z_bin"]],
        "lapse_bin": row["lapse_bin"],
        "lapse": lapse_labels[row["lapse_bin"]],
        "count": row["count"],
        "fraction": row["count"] / total if total else 0.0,
    }


def _global_ratio(
    event: Mapping[str, int], numerator: str, denominator: str,
) -> dict[str, Any]:
    numerator_value = event[numerator]
    denominator_value = event[denominator]
    defined = denominator_value != 0
    return {
        "numerator": numerator,
        "numerator_value": numerator_value,
        "denominator": denominator,
        "denominator_value": denominator_value,
        "defined": defined,
        "value": numerator_value / denominator_value if defined else None,
    }


def _validate_closed_cycle(
    summary: Mapping[str, int], event: Mapping[str, int],
    bins: Sequence[Mapping[str, Any]], source: str, cycle_index: int,
) -> tuple[int, int]:
    """Validate one complete ordered group and return its bin/prefix-bin sums."""

    cycle = summary["cycle"]
    if summary["count"] != summary["nfofc"]:
        _failure(
            source,
            summary["line"],
            f"cycle {cycle}: summary.count={summary['count']} differs from "
            f"summary.nfofc={summary['nfofc']}",
        )
    if summary["nfofc"] != event["fofc"]:
        _failure(
            source,
            summary["line"],
            f"cycle {cycle}: summary.nfofc={summary['nfofc']} differs from "
            f"event.fofc={event['fofc']}",
        )
    bin_sum = sum(row["count"] for row in bins)
    if bin_sum != event["fofc"]:
        _failure(
            source,
            summary["line"],
            f"cycle {cycle}: histogram sum={bin_sum} differs from "
            f"event.fofc={event['fofc']}",
        )

    unattributed = summary["unattributed"]
    canonical_count = sum(
        row["count"] for row in bins
        if _bin_key(row) == CANONICAL_UNATTRIBUTED_KEY
    )
    if canonical_count != unattributed:
        _failure(
            source,
            summary["line"],
            f"cycle {cycle}: canonical overflow/other/unknown bin count="
            f"{canonical_count} differs from unattributed={unattributed} "
            "(level_bin=32 stage_bin=0 reason=unknown r_cyl_bin=6 "
            "abs_z_bin=6 lapse_bin=5)",
        )
    for row in bins:
        if row["stage_bin"] == 0 and _bin_key(row) != CANONICAL_UNATTRIBUTED_KEY:
            _failure(
                source,
                row["line"],
                f"cycle {cycle}: stage_bin=0 is reserved for the canonical "
                "unattributed prefix bin",
            )
    if unattributed > 0:
        if cycle_index != 0:
            _failure(
                source,
                summary["line"],
                f"cycle {cycle}: unattributed={unattributed} is permitted only "
                "in the first telemetry cycle",
            )
    return bin_sum, canonical_count


def analyze_event_log(
    path: Path, cycle: int, *, top_bins: int = 20,
    require_unattributed_zero: bool = False,
) -> dict[str, Any]:
    """Validate *cycle* and return its conservative FOFC spatial summary."""

    if not isinstance(cycle, int) or isinstance(cycle, bool) or cycle < 0:
        raise ValueError("cycle must be a non-negative integer")
    if cycle > INT_MAX:
        raise ValueError(f"cycle exceeds {INT_MAX}")
    if not isinstance(top_bins, int) or isinstance(top_bins, bool) or top_bins < 0:
        raise ValueError("top_bins must be a non-negative integer")
    if not isinstance(require_unattributed_zero, bool):
        raise ValueError("require_unattributed_zero must be a boolean")

    line_iterator, read_state, checked = _read_lines_stably(path)
    source = str(checked)
    headers: list[tuple[str, ...]] = []
    schema: dict[str, Any] | None = None
    active_summary: dict[str, int] | None = None
    active_bins: list[dict[str, Any]] = []
    active_keys: set[tuple[Any, ...]] = set()
    selected_group: tuple[
        dict[str, int], dict[str, int], list[dict[str, Any]], int, int, int
    ] | None = None
    last_completed_cycle: int | None = None
    telemetry_records = 0
    event_rows = 0
    completed_groups = 0

    for line_number, line in line_iterator:
        fields = line.split()
        if not fields:
            continue
        if fields[:2] == ["#", "fofc_spatial_v1"]:
            telemetry_records += 1
            record = _parse_record(line, source, line_number)
            kind = record["kind"]
            if kind == "schema":
                if headers != [EVENT_COLUMNS]:
                    _failure(
                        source,
                        line_number,
                        "schema record does not follow one exact event header",
                    )
                if schema is not None:
                    _failure(source, line_number, "duplicate schema record")
                if active_summary is not None or event_rows != 0:
                    _failure(source, line_number, "schema record is out of order")
                schema = _parse_schema(record, source, line_number)
            elif kind == "summary":
                if schema is None:
                    _failure(source, line_number, "summary record precedes schema")
                summary = _parse_summary(record, source, line_number)
                if active_summary is not None:
                    _failure(
                        source,
                        line_number,
                        f"cycle {active_summary['cycle']}: missing event row before "
                        f"summary for cycle {summary['cycle']}",
                    )
                if (
                    last_completed_cycle is not None
                    and summary["cycle"] <= last_completed_cycle
                ):
                    _failure(
                        source,
                        line_number,
                        "telemetry cycles are not strictly increasing: "
                        f"previous={last_completed_cycle}, "
                        f"current={summary['cycle']}",
                    )
                if (
                    last_completed_cycle is not None
                    and summary["cycle"] != last_completed_cycle + 1
                ):
                    _failure(
                        source,
                        line_number,
                        "telemetry cycles are not consecutive under the dcycle=1 "
                        f"contract: expected={last_completed_cycle + 1}, "
                        f"current={summary['cycle']}",
                    )
                active_summary = summary
                active_bins = []
                active_keys = set()
            elif kind == "bin":
                if schema is None:
                    _failure(source, line_number, "bin record precedes schema")
                histogram_bin = _parse_bin(record, source, line_number)
                if active_summary is None:
                    _failure(
                        source,
                        line_number,
                        f"bin for cycle {histogram_bin['cycle']} has no active summary",
                    )
                if histogram_bin["cycle"] != active_summary["cycle"]:
                    _failure(
                        source,
                        line_number,
                        f"bin cycle {histogram_bin['cycle']} differs from active "
                        f"summary cycle {active_summary['cycle']}",
                    )
                key = _bin_key(histogram_bin)
                if key in active_keys:
                    _failure(
                        source,
                        line_number,
                        f"duplicate histogram key for cycle "
                        f"{active_summary['cycle']}: {key}",
                    )
                active_keys.add(key)
                active_bins.append(histogram_bin)
            else:
                _failure(source, line_number, f"unknown telemetry kind {kind!r}")
            continue
        if fields[0] == "#" and len(fields) > 1 and fields[1] == "cycle":
            header = tuple(fields[1:])
            if header != EVENT_COLUMNS:
                _failure(source, line_number, "event header does not match schema")
            if schema is not None or active_summary is not None or event_rows != 0:
                _failure(source, line_number, "event header is out of order")
            if headers:
                _failure(source, line_number, "duplicate event header")
            headers.append(header)
            continue
        if fields[0].startswith("#"):
            continue
        if headers != [EVENT_COLUMNS]:
            _failure(source, line_number, "event data does not follow one unique header")
        if schema is None:
            _failure(source, line_number, "event row precedes schema")
        event = _parse_event_row(fields, source, line_number)
        event_rows += 1
        if active_summary is None:
            _failure(
                source,
                line_number,
                f"event row for cycle {event['cycle']} has no active summary",
            )
        if event["cycle"] != active_summary["cycle"]:
            _failure(
                source,
                line_number,
                f"event cycle {event['cycle']} differs from active summary cycle "
                f"{active_summary['cycle']}",
            )
        bin_sum, canonical_count = _validate_closed_cycle(
            active_summary, event, active_bins, source, completed_groups
        )
        completed_cycle = active_summary["cycle"]
        if completed_cycle == cycle:
            if selected_group is not None:
                _failure(source, line_number, f"duplicate telemetry cycle {cycle}")
            selected_group = (
                active_summary,
                event,
                list(active_bins),
                bin_sum,
                canonical_count,
                completed_groups,
            )
        completed_groups += 1
        last_completed_cycle = completed_cycle
        active_summary = None
        active_bins = []
        active_keys = set()

    input_record = _finish_stable_read(read_state)
    if headers != [EVENT_COLUMNS]:
        _failure(
            source,
            None,
            f"event header must be exactly one copy of {EVENT_COLUMNS}",
        )
    if schema is None:
        _failure(source, None, "expected one schema record, found 0")
    if active_summary is not None:
        _failure(
            source,
            active_summary["line"],
            f"cycle {active_summary['cycle']}: incomplete telemetry group at EOF; "
            "missing event row",
        )
    if selected_group is None:
        _failure(
            source,
            None,
            f"cycle {cycle}: expected one complete summary/bin/event group, found 0",
        )

    summary, event, selected_bins, bin_sum, canonical_count, cycle_index = selected_group
    if require_unattributed_zero and summary["unattributed"] != 0:
        _failure(
            source,
            summary["line"],
            f"cycle {cycle}: unattributed={summary['unattributed']} violates "
            "require_unattributed_zero policy",
        )

    total = event["fofc"]
    level_labels = [str(value) for value in range(32)] + ["overflow"]
    radius_labels = _interval_labels(R_CYL_EDGES)
    z_labels = _interval_labels(ABS_Z_EDGES)
    lapse_labels = _interval_labels(LAPSE_EDGES)
    radius_labels[-1] = "other(r_cyl>=64_or_nonfinite)"
    z_labels[-1] = "other(abs_z>=16_or_nonfinite)"
    lapse_labels[-1] = "other(lapse<0,lapse>=1,or_nonfinite)"
    marginals = {
        "level": _marginal(selected_bins, "level_bin", level_labels, total),
        "stage": _marginal(
            selected_bins, "stage_bin", STAGE_LABELS, total
        ),
        "reason": _reason_marginal(selected_bins, total),
        "r_cyl": _marginal(
            selected_bins, "r_cyl_bin", radius_labels, total
        ),
        "abs_z": _marginal(selected_bins, "abs_z_bin", z_labels, total),
        "lapse": _marginal(
            selected_bins, "lapse_bin", lapse_labels, total
        ),
    }
    marginal_sums = {name: value["sum"] for name, value in marginals.items()}
    if any(value != total for value in marginal_sums.values()):
        raise AssertionError("internal marginal conservation failure")

    ordered_bins = sorted(
        selected_bins,
        key=lambda row: (-row["count"], _bin_key(row)),
    )
    clean_summary = {name: summary[name] for name in ("count", "nfofc", "unattributed")}
    return {
        "kind": "athenak_fofc_spatial_v1_cycle",
        "valid": True,
        "input": input_record,
        "cycle": cycle,
        "policy": {
            "require_unattributed_zero": require_unattributed_zero,
        },
        "schema": schema,
        "event": event,
        "summary": clean_summary,
        "unattributed": {
            "count": summary["unattributed"],
            "present": summary["unattributed"] > 0,
            "selected_cycle_is_first_telemetry_cycle": cycle_index == 0,
            "permitted_only_on_first_telemetry_cycle": True,
            "canonical_bin": {
                "level_bin": CANONICAL_UNATTRIBUTED_KEY[0],
                "stage_bin": CANONICAL_UNATTRIBUTED_KEY[1],
                "reason": "unknown",
                "r_cyl_bin": CANONICAL_UNATTRIBUTED_KEY[3],
                "abs_z_bin": CANONICAL_UNATTRIBUTED_KEY[4],
                "lapse_bin": CANONICAL_UNATTRIBUTED_KEY[5],
            },
            "canonical_bin_count": canonical_count,
            "canonical_bin_equals_unattributed": (
                canonical_count == summary["unattributed"]
            ),
            "policy_satisfied": (
                not require_unattributed_zero or summary["unattributed"] == 0
            ),
        },
        "global_ratios": {
            "fofc_per_test": _global_ratio(event, "fofc", "fofc_tests"),
            "cons_adjust_per_c2p_call": _global_ratio(
                event, "cons_adjust", "c2p_calls"
            ),
            "mag_adjust_per_c2p_call": _global_ratio(
                event, "mag_adjust", "c2p_calls"
            ),
        },
        "checks": {
            "unique_schema": True,
            "unique_summary": True,
            "unique_event_row": True,
            "all_cycles_ordered_and_closed": True,
            "cycles_strictly_increasing": True,
            "cycles_consecutive_under_dcycle_one": True,
            "all_cycles_conservative": True,
            "summary_count_equals_nfofc": True,
            "summary_nfofc_equals_event_fofc": True,
            "unattributed_zero": summary["unattributed"] == 0,
            "unattributed_policy_satisfied": True,
            "unattributed_canonical_bin_exact": (
                canonical_count == summary["unattributed"]
            ),
            "bin_keys_unique": True,
            "bin_ranges_valid": True,
            "bin_sum": bin_sum,
            "bin_sum_equals_event_fofc": True,
            "marginal_sums": marginal_sums,
            "marginals_conservative": True,
        },
        "records": {
            "telemetry_total": telemetry_records,
            "event_rows_total": event_rows,
            "complete_cycle_groups": completed_groups,
            "selected_cycle_group_index": cycle_index,
            "selected_nonzero_bins": len(selected_bins),
        },
        "marginals": marginals,
        "top_bins": [
            _decorated_bin(row, total) for row in ordered_bins[:top_bins]
        ],
    }


def _nonnegative_integer(text: str) -> int:
    if UNSIGNED_PATTERN.fullmatch(text) is None:
        raise argparse.ArgumentTypeError("must be a non-negative integer")
    canonical = text.lstrip("0") or "0"
    limit = str(INT_MAX)
    if len(canonical) > len(limit) or (
        len(canonical) == len(limit) and canonical > limit
    ):
        raise argparse.ArgumentTypeError(f"must not exceed {INT_MAX}")
    return int(canonical)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("event_log", type=Path)
    parser.add_argument("--cycle", required=True, type=_nonnegative_integer)
    parser.add_argument("--top-bins", type=_nonnegative_integer, default=20)
    parser.add_argument(
        "--require-unattributed-zero",
        action="store_true",
        help=(
            "reject a selected cycle with a restart/IC prefix that cannot be "
            "spatially attributed"
        ),
    )
    style = parser.add_mutually_exclusive_group()
    style.add_argument(
        "--compact",
        action="store_const",
        dest="json_style",
        const="compact",
        help="emit one-line compact JSON",
    )
    style.add_argument(
        "--pretty",
        action="store_const",
        dest="json_style",
        const="pretty",
        help="emit indented JSON (default)",
    )
    parser.set_defaults(json_style="pretty")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        result = analyze_event_log(
            args.event_log,
            args.cycle,
            top_bins=args.top_bins,
            require_unattributed_zero=args.require_unattributed_zero,
        )
    except (FofcTelemetryError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if args.json_style == "compact":
        print(json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False))
    else:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
