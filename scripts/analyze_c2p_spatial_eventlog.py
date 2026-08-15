#!/usr/bin/env python3
"""Strictly validate and summarize AthenaK C2P spatial telemetry.

``# c2p_spatial_v1`` records are diagnostic comments adjacent to the traditional
event row.  This reader treats each cycle as one conservative object: both
intervention summaries must be present, joint/stage/geometry histograms must each
sum to the matching traditional counter, keys must be unique and in range, and
restart-carried unattributed counts are permitted only in the first telemetry cycle.

The input is streamed from a regular, non-symlink file opened with ``O_NOFOLLOW``.
Its device, inode, size, mtime, and ctime must remain unchanged during the read.
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
from typing import Any, Iterable, Mapping, Sequence


EVENT_COLUMNS = (
    "cycle", "eos_dfloor", "eos_efloor", "eos_tfloor", "eos_vceil",
    "eos_fail", "c2p_it", "fofc", "cons_adjust", "mag_adjust",
    "c2p_calls", "fofc_tests",
)
INTERVENTIONS = ("cons_adjust", "mag_adjust")
LEVEL_BINS = 33
R_CYL_EDGES = (2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
ABS_Z_EDGES = (0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
LAPSE_EDGES = (0.2, 0.4, 0.6, 0.8, 1.0)
DENSITY_RATIO_EDGES = (1.0, 2.0, 4.0, 16.0, 64.0, 256.0)
MAG_RATIO_EDGES = (0.01, 0.1, 0.5, 1.0, 2.0, 10.0)
QUANTITY_BINS = 8
QUANTITY_INVALID_BIN = 7
STAGE_BINS = 4
GEOMETRY_BINS = 2
UINT64_MAX = (1 << 64) - 1
INT_MAX = (1 << 31) - 1
KEY_PATTERN = re.compile(r"[a-z][a-z0-9_]*\Z")
UNSIGNED_PATTERN = re.compile(r"[0-9]+\Z")

SCHEMA_KEYS = frozenset({
    "kind", "intervention_bins", "level_bins", "r_cyl_edges",
    "abs_z_edges", "lapse_edges", "density_floor_ratio_edges",
    "magnetization_limit_ratio_edges", "quantity_invalid_bin",
    "stage_bins", "geometry_bins", "center1", "center2", "center3",
})
SUMMARY_KEYS = frozenset({
    "kind", "cycle", "intervention", "count", "authoritative", "unattributed",
})
BIN_KEYS = frozenset({
    "kind", "cycle", "intervention", "level_bin", "r_cyl_bin", "abs_z_bin",
    "lapse_bin", "density_floor_ratio_bin",
    "magnetization_limit_ratio_bin", "count",
})
STAGE_KEYS = frozenset({
    "kind", "cycle", "intervention", "stage_bin", "count",
})
GEOMETRY_KEYS = frozenset({
    "kind", "cycle", "intervention", "valid", "count",
})


class C2PTelemetryError(ValueError):
    """The event log does not satisfy the C2P spatial-v1 contract."""


def _fail(source: str, line: int | None, message: str) -> None:
    location = source if line is None else f"{source}:{line}"
    raise C2PTelemetryError(f"{location}: {message}")


def _exact_keys(
    record: Mapping[str, str], expected: frozenset[str], source: str, line: int,
) -> None:
    actual = set(record)
    if actual != expected:
        _fail(
            source, line,
            f"record keys are not exact; missing={sorted(expected-actual)}, "
            f"extra={sorted(actual-expected)}",
        )


def _unsigned(
    text: str, name: str, source: str, line: int, maximum: int = UINT64_MAX,
) -> int:
    if UNSIGNED_PATTERN.fullmatch(text) is None:
        _fail(source, line, f"{name} is not an unsigned decimal integer")
    canonical = text.lstrip("0") or "0"
    limit = str(maximum)
    if len(canonical) > len(limit) or (
        len(canonical) == len(limit) and canonical > limit
    ):
        _fail(source, line, f"{name} exceeds {maximum}")
    return int(canonical)


def _bounded(
    record: Mapping[str, str], name: str, bins: int, source: str, line: int,
) -> int:
    value = _unsigned(record[name], name, source, line, INT_MAX)
    if value >= bins:
        _fail(source, line, f"{name}={value} is outside [0,{bins})")
    return value


def _record(line: str, source: str, line_number: int) -> dict[str, str]:
    fields = line.split()
    if fields[:2] != ["#", "c2p_spatial_v1"]:
        _fail(source, line_number, "invalid C2P telemetry prefix")
    parsed: dict[str, str] = {}
    for field in fields[2:]:
        if "=" not in field:
            _fail(source, line_number, f"telemetry token lacks '=': {field!r}")
        key, value = field.split("=", 1)
        if KEY_PATTERN.fullmatch(key) is None or not value:
            _fail(source, line_number, f"invalid telemetry token: {field!r}")
        if key in parsed:
            _fail(source, line_number, f"duplicate telemetry key {key!r}")
        parsed[key] = value
    if "kind" not in parsed:
        _fail(source, line_number, "telemetry record lacks kind")
    return parsed


def _finite(value: str, name: str, source: str, line: int) -> float:
    try:
        result = float(value)
    except ValueError:
        _fail(source, line, f"{name} is not a floating-point number")
    if not math.isfinite(result):
        _fail(source, line, f"{name} is not finite")
    return result


def _csv_edges(
    value: str, name: str, expected: Sequence[float], source: str, line: int,
) -> list[float]:
    actual = [_finite(item, name, source, line) for item in value.split(",")]
    if actual != list(expected):
        _fail(source, line, f"{name} changed: expected {list(expected)}, got {actual}")
    return actual


def _parse_schema(
    record: Mapping[str, str], source: str, line: int,
) -> dict[str, Any]:
    _exact_keys(record, SCHEMA_KEYS, source, line)
    expected_literals = {
        "intervention_bins": ",".join(INTERVENTIONS),
        "level_bins": "0..31,overflow",
        "quantity_invalid_bin": str(QUANTITY_INVALID_BIN),
        "stage_bins": "other,1,2,3",
        "geometry_bins": "invalid,valid",
    }
    for name, expected in expected_literals.items():
        if record[name] != expected:
            _fail(
                source, line,
                f"{name} changed: expected {expected!r}, got {record[name]!r}",
            )
    return {
        "version": 1,
        "intervention_bins": list(INTERVENTIONS),
        "level_bins": [str(value) for value in range(32)] + ["overflow"],
        "r_cyl_edges": _csv_edges(
            record["r_cyl_edges"], "r_cyl_edges", R_CYL_EDGES, source, line),
        "abs_z_edges": _csv_edges(
            record["abs_z_edges"], "abs_z_edges", ABS_Z_EDGES, source, line),
        "lapse_edges": _csv_edges(
            record["lapse_edges"], "lapse_edges", LAPSE_EDGES, source, line),
        "density_floor_ratio_edges": _csv_edges(
            record["density_floor_ratio_edges"], "density_floor_ratio_edges",
            DENSITY_RATIO_EDGES, source, line),
        "magnetization_limit_ratio_edges": _csv_edges(
            record["magnetization_limit_ratio_edges"],
            "magnetization_limit_ratio_edges", MAG_RATIO_EDGES, source, line),
        "quantity_invalid_bin": QUANTITY_INVALID_BIN,
        "stage_bins": ["other", "1", "2", "3"],
        "geometry_bins": ["invalid", "valid"],
        "center": [
            _finite(record[f"center{axis}"], f"center{axis}", source, line)
            for axis in (1, 2, 3)
        ],
    }


def _new_group(summary: dict[str, Any]) -> dict[str, Any]:
    return {
        "summary": summary,
        "bins": Counter(),
        "stages": Counter(),
        "geometry": Counter(),
    }


def _counter_marginal(
    bins: Mapping[tuple[int, ...], int], dimension: int,
) -> dict[str, int]:
    result: Counter[int] = Counter()
    for key, count in bins.items():
        result[key[dimension]] += count
    return {str(key): result[key] for key in sorted(result)}


def _fraction(count: int, total: int) -> float:
    return count/total if total else 0.0


def analyze_lines(
    lines: Iterable[str], *, source: str = "<memory>", cycle: int | None = None,
    require_unattributed_zero: bool = False,
) -> dict[str, Any]:
    schema: dict[str, Any] | None = None
    groups: dict[tuple[int, str], dict[str, Any]] = {}
    events: dict[int, dict[str, int]] = {}
    event_header: tuple[str, ...] | None = None

    for line_number, raw_line in enumerate(lines, 1):
        line = raw_line.rstrip("\r\n")
        if line.startswith("#  cycle"):
            names = tuple(line[1:].split())
            if names != EVENT_COLUMNS:
                _fail(source, line_number, f"traditional event schema changed: {names}")
            if event_header is not None:
                _fail(source, line_number, "duplicate traditional event header")
            event_header = names
            continue
        if line.startswith("# c2p_spatial_v1 "):
            record = _record(line, source, line_number)
            kind = record["kind"]
            if kind == "schema":
                if schema is not None:
                    _fail(source, line_number, "duplicate C2P telemetry schema")
                schema = _parse_schema(record, source, line_number)
                continue
            if schema is None:
                _fail(source, line_number, "C2P record precedes schema")
            expected = {
                "summary": SUMMARY_KEYS, "bin": BIN_KEYS,
                "stage": STAGE_KEYS, "geometry": GEOMETRY_KEYS,
            }.get(kind)
            if expected is None:
                _fail(source, line_number, f"unknown C2P record kind {kind!r}")
            _exact_keys(record, expected, source, line_number)
            record_cycle = _unsigned(
                record["cycle"], "cycle", source, line_number, INT_MAX
            )
            intervention = record["intervention"]
            if intervention not in INTERVENTIONS:
                _fail(source, line_number, f"unknown intervention {intervention!r}")
            group_key = (record_cycle, intervention)
            if kind == "summary":
                if group_key in groups:
                    _fail(source, line_number, f"duplicate summary for {group_key}")
                groups[group_key] = _new_group({
                    "cycle": record_cycle,
                    "intervention": intervention,
                    "count": _unsigned(record["count"], "count", source, line_number),
                    "authoritative": _unsigned(
                        record["authoritative"], "authoritative", source, line_number),
                    "unattributed": _unsigned(
                        record["unattributed"], "unattributed", source, line_number),
                    "line": line_number,
                })
                continue
            if group_key not in groups:
                _fail(source, line_number, f"{kind} precedes summary for {group_key}")
            count = _unsigned(record["count"], "count", source, line_number)
            if count == 0:
                _fail(source, line_number, f"zero-count {kind} record is not canonical")
            group = groups[group_key]
            if kind == "bin":
                key = (
                    _bounded(record, "level_bin", LEVEL_BINS, source, line_number),
                    _bounded(
                        record, "r_cyl_bin", len(R_CYL_EDGES)+1,
                        source, line_number,
                    ),
                    _bounded(
                        record, "abs_z_bin", len(ABS_Z_EDGES)+1,
                        source, line_number,
                    ),
                    _bounded(
                        record, "lapse_bin", len(LAPSE_EDGES)+1,
                        source, line_number,
                    ),
                    _bounded(record, "density_floor_ratio_bin", QUANTITY_BINS,
                             source, line_number),
                    _bounded(record, "magnetization_limit_ratio_bin", QUANTITY_BINS,
                             source, line_number),
                )
                target = group["bins"]
            elif kind == "stage":
                key = (_bounded(record, "stage_bin", STAGE_BINS, source, line_number),)
                target = group["stages"]
            else:
                key = (_bounded(record, "valid", GEOMETRY_BINS, source, line_number),)
                target = group["geometry"]
            if key in target:
                _fail(source, line_number, f"duplicate {kind} key {key} for {group_key}")
            target[key] = count
            continue
        if line and not line.startswith("#"):
            if event_header is None:
                _fail(source, line_number, "traditional event row precedes its header")
            fields = line.split()
            if len(fields) != len(EVENT_COLUMNS):
                _fail(source, line_number, "traditional event row has wrong field count")
            values = [
                _unsigned(value, name, source, line_number,
                          INT_MAX if name in ("cycle", "c2p_it") else UINT64_MAX)
                for name, value in zip(EVENT_COLUMNS, fields)
            ]
            row = dict(zip(EVENT_COLUMNS, values))
            if row["cycle"] in events:
                _fail(source, line_number, f"duplicate traditional cycle {row['cycle']}")
            events[row["cycle"]] = row

    if schema is None:
        _fail(source, None, "no C2P spatial-v1 schema found")
    telemetry_cycles = sorted({key[0] for key in groups})
    if not telemetry_cycles:
        _fail(source, None, "no C2P telemetry summaries found")
    for previous, current in zip(telemetry_cycles, telemetry_cycles[1:]):
        if current != previous + 1:
            _fail(
                source, None,
                f"telemetry cycles are not consecutive: {previous}->{current}",
            )

    first_cycle = telemetry_cycles[0]
    for telemetry_cycle in telemetry_cycles:
        if telemetry_cycle not in events:
            _fail(source, None, f"cycle {telemetry_cycle}: missing traditional event row")
        row = events[telemetry_cycle]
        for intervention in INTERVENTIONS:
            key = (telemetry_cycle, intervention)
            if key not in groups:
                _fail(
                    source, None,
                    f"cycle {telemetry_cycle}: missing {intervention} summary",
                )
            group = groups[key]
            summary = group["summary"]
            expected = row[intervention]
            joint_total = sum(group["bins"].values())
            stage_total = sum(group["stages"].values())
            geometry_total = sum(group["geometry"].values())
            if not (summary["count"] == summary["authoritative"] == expected):
                _fail(
                    source, summary["line"],
                    f"{key}: summary count/authoritative/event mismatch: "
                    f"{summary['count']}/{summary['authoritative']}/{expected}",
                )
            if (joint_total, stage_total, geometry_total) != (
                expected, expected, expected
            ):
                _fail(
                    source, summary["line"],
                    f"{key}: joint/stage/geometry totals "
                    f"{joint_total}/{stage_total}/{geometry_total} != {expected}",
                )
            unattributed = summary["unattributed"]
            if unattributed and telemetry_cycle != first_cycle:
                _fail(
                    source, summary["line"],
                    "unattributed count is allowed only in first cycle",
                )
            if unattributed:
                canonical_joint = (
                    LEVEL_BINS-1, len(R_CYL_EDGES), len(ABS_Z_EDGES),
                    len(LAPSE_EDGES), QUANTITY_INVALID_BIN, QUANTITY_INVALID_BIN,
                )
                if (group["bins"].get(canonical_joint, 0) < unattributed or
                    group["stages"].get((0,), 0) < unattributed or
                    group["geometry"].get((0,), 0) < unattributed):
                    _fail(
                        source, summary["line"],
                        "unattributed prefix is not in canonical bins",
                    )

    selected_cycle = telemetry_cycles[-1] if cycle is None else cycle
    if selected_cycle not in telemetry_cycles:
        _fail(source, None, f"selected cycle {selected_cycle} has no C2P telemetry")
    interventions: dict[str, Any] = {}
    for intervention in INTERVENTIONS:
        group = groups[(selected_cycle, intervention)]
        total = group["summary"]["count"]
        marginals = {
            "level_bin": _counter_marginal(group["bins"], 0),
            "r_cyl_bin": _counter_marginal(group["bins"], 1),
            "abs_z_bin": _counter_marginal(group["bins"], 2),
            "lapse_bin": _counter_marginal(group["bins"], 3),
            "density_floor_ratio_bin": _counter_marginal(group["bins"], 4),
            "magnetization_limit_ratio_bin": _counter_marginal(group["bins"], 5),
            "stage_bin": {
                str(key[0]): value
                for key, value in sorted(group["stages"].items())
            },
            "geometry_valid": {
                str(key[0]): value for key, value in sorted(group["geometry"].items())
            },
        }
        density = {
            int(key): value
            for key, value in marginals["density_floor_ratio_bin"].items()
        }
        magnetization = {
            int(key): value
            for key, value in marginals["magnetization_limit_ratio_bin"].items()
        }
        radius = {int(key): value for key, value in marginals["r_cyl_bin"].items()}
        geometry = {int(key): value for key, value in marginals["geometry_valid"].items()}
        invalid_quantity_count = sum(
            count for key, count in group["bins"].items()
            if key[4] == QUANTITY_INVALID_BIN or key[5] == QUANTITY_INVALID_BIN
        )
        interventions[intervention] = {
            "summary": {
                name: value for name, value in group["summary"].items()
                if name != "line"
            },
            "nonzero_joint_bins": len(group["bins"]),
            "marginals": marginals,
            "derived_fractions": {
                "r_cyl_lt_16": _fraction(
                    sum(radius.get(index, 0) for index in range(4)), total
                ),
                "density_lt_16_times_floor": _fraction(
                    sum(density.get(index, 0) for index in range(4)), total),
                "density_ge_256_times_floor": _fraction(density.get(6, 0), total),
                "magnetization_at_or_above_limit": _fraction(
                    sum(magnetization.get(index, 0) for index in (4, 5, 6)), total),
                "quantity_invalid": _fraction(invalid_quantity_count, total),
                "geometry_invalid": _fraction(geometry.get(0, 0), total),
            },
        }
    if require_unattributed_zero and any(
        interventions[name]["summary"]["unattributed"] != 0 for name in INTERVENTIONS
    ):
        _fail(source, None, f"cycle {selected_cycle}: unattributed telemetry is nonzero")
    return {
        "source": source,
        "schema": schema,
        "telemetry_cycles": telemetry_cycles,
        "selected_cycle": selected_cycle,
        "selected_cycle_is_first_telemetry_cycle": selected_cycle == first_cycle,
        "event": events[selected_cycle],
        "interventions": interventions,
    }


def analyze_path(
    path: Path, *, cycle: int | None = None,
    require_unattributed_zero: bool = False,
) -> dict[str, Any]:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise C2PTelemetryError(f"{path}: cannot open read-only: {exc}") from exc
    digest = hashlib.sha256()
    try:
        with os.fdopen(descriptor, "rb", closefd=True) as stream:
            before = os.fstat(stream.fileno())
            if not stat.S_ISREG(before.st_mode):
                _fail(str(path), None, "input is not a regular file")

            def decoded_lines() -> Iterable[str]:
                for line_number, raw in enumerate(stream, 1):
                    digest.update(raw)
                    try:
                        yield raw.decode("utf-8")
                    except UnicodeDecodeError as exc:
                        _fail(str(path), line_number, f"input is not UTF-8: {exc}")

            report = analyze_lines(
                decoded_lines(), source=str(path), cycle=cycle,
                require_unattributed_zero=require_unattributed_zero,
            )
            after = os.fstat(stream.fileno())
    except C2PTelemetryError:
        raise
    def identity(value: os.stat_result) -> tuple[int, int, int, int, int]:
        return (
            value.st_dev, value.st_ino, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns,
        )

    if identity(before) != identity(after):
        _fail(str(path), None, "input changed during streaming read")
    report["source_identity"] = {
        "device": before.st_dev, "inode": before.st_ino, "size": before.st_size,
        "mtime_ns": before.st_mtime_ns, "ctime_ns": before.st_ctime_ns,
        "sha256": digest.hexdigest(),
    }
    return report


def _text(report: Mapping[str, Any]) -> str:
    lines = [
        f"C2P spatial telemetry cycle {report['selected_cycle']}",
        f"source: {report['source']}",
    ]
    for name in INTERVENTIONS:
        item = report["interventions"][name]
        summary = item["summary"]
        fractions = item["derived_fractions"]
        lines.append(
            f"{name}: count={summary['count']} unattributed={summary['unattributed']} "
            f"R_cyl<16={fractions['r_cyl_lt_16']:.6%} "
            f"D<16*floor={fractions['density_lt_16_times_floor']:.6%} "
            f"D>=256*floor={fractions['density_ge_256_times_floor']:.6%} "
            f"mag>=limit={fractions['magnetization_at_or_above_limit']:.6%} "
            f"invalid_geometry={fractions['geometry_invalid']:.6%}"
        )
        lines.append(
            "  level bins: " + json.dumps(item["marginals"]["level_bin"], sort_keys=True)
        )
        lines.append(
            "  lapse bins: " + json.dumps(item["marginals"]["lapse_bin"], sort_keys=True)
        )
        lines.append(
            "  stage bins: " + json.dumps(item["marginals"]["stage_bin"], sort_keys=True)
        )
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("event_log", type=Path)
    parser.add_argument("--cycle", type=int)
    parser.add_argument("--require-unattributed-zero", action="store_true")
    parser.add_argument("--json", action="store_true", dest="as_json")
    args = parser.parse_args(argv)
    try:
        report = analyze_path(
            args.event_log, cycle=args.cycle,
            require_unattributed_zero=args.require_unattributed_zero,
        )
    except C2PTelemetryError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    if args.as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(_text(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
