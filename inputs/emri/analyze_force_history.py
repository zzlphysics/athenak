#!/usr/bin/env python3
"""Time-average EMRI wind-tunnel force histories and form controlled contrasts."""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
import re
import statistics
import sys
from pathlib import Path


HEADER_PATTERN = re.compile(r"\[\d+\]=([^\s]+)")
AXES = ("x", "y", "z")


class ForceHistory:
    def __init__(self, label: str, path: Path) -> None:
        self.label = label
        self.path = path
        self.times: list[float] = []
        self.values: dict[str, list[float]] = {}
        self._read()

    def _read(self) -> None:
        columns: list[str] = []
        rows_by_time: dict[float, dict[str, float]] = {}
        with self.path.open(encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, start=1):
                if line.startswith("#"):
                    names = HEADER_PATTERN.findall(line)
                    if names:
                        columns = names
                    continue
                if not line.strip():
                    continue
                if not columns:
                    raise ValueError(f"{self.path}: data precedes the column header")
                fields = line.split()
                if len(fields) != len(columns):
                    raise ValueError(
                        f"{self.path}:{line_number}: expected {len(columns)} values, "
                        f"found {len(fields)}"
                    )
                values = [float(field) for field in fields]
                if not all(math.isfinite(value) for value in values):
                    raise ValueError(f"{self.path}:{line_number}: non-finite value")
                row = dict(zip(columns, values))
                rows_by_time[row["time"]] = row  # retain the last copy after a restart

        required = {"time", "mdot", "geo_resid"}
        for axis in AXES:
            required.update({f"Fmom_{axis}", f"Fnewt_{axis}"})
            for radius in (1, 2, 3):
                required.add(f"Frel{radius}_{axis}")
        missing = required.difference(columns)
        if missing:
            raise ValueError(f"{self.path}: missing columns: {', '.join(sorted(missing))}")
        if len(rows_by_time) < 2:
            raise ValueError(f"{self.path}: at least two distinct history times are required")

        self.times = sorted(rows_by_time)
        rows = [rows_by_time[time] for time in self.times]
        for name in columns:
            self.values[name] = [row[name] for row in rows]
        for axis in AXES:
            momentum = self.values[f"Fmom_{axis}"]
            newtonian = self.values[f"Fnewt_{axis}"]
            self.values[f"Ftotal_newt_{axis}"] = [
                -mom + grav for mom, grav in zip(momentum, newtonian)
            ]
            for radius in (1, 2, 3):
                relativistic = self.values[f"Frel{radius}_{axis}"]
                self.values[f"Ftotal{radius}_{axis}"] = [
                    -mom + grav for mom, grav in zip(momentum, relativistic)
                ]
            self.values[f"dFrel21_{axis}"] = [
                outer - inner
                for outer, inner in zip(
                    self.values[f"Frel2_{axis}"], self.values[f"Frel1_{axis}"]
                )
            ]
            self.values[f"dFrel32_{axis}"] = [
                outer - inner
                for outer, inner in zip(
                    self.values[f"Frel3_{axis}"], self.values[f"Frel2_{axis}"]
                )
            ]

    def value_at(self, quantity: str, time: float) -> float:
        values = self.values[quantity]
        if time <= self.times[0]:
            return values[0]
        if time >= self.times[-1]:
            return values[-1]
        upper = bisect.bisect_right(self.times, time)
        lower = upper - 1
        fraction = ((time-self.times[lower])
                    /(self.times[upper]-self.times[lower]))
        return values[lower] + fraction*(values[upper]-values[lower])

    def average(self, quantity: str, start: float, stop: float) -> float:
        knots = [start]
        knots.extend(time for time in self.times if start < time < stop)
        knots.append(stop)
        integral = 0.0
        previous_time = knots[0]
        previous_value = self.value_at(quantity, previous_time)
        for time in knots[1:]:
            value = self.value_at(quantity, time)
            integral += 0.5*(previous_value+value)*(time-previous_time)
            previous_time = time
            previous_value = value
        return integral/(stop-start)


def default_quantities() -> list[str]:
    names = ["mdot"]
    for prefix in ("Ftotal1", "Ftotal2", "Ftotal3", "dFrel21", "dFrel32"):
        names.extend(f"{prefix}_{axis}" for axis in AXES)
    return names


def summarize(history: ForceHistory, quantities: list[str], start: float, stop: float,
              requested_blocks: int, kind: str = "case") -> list[dict[str, object]]:
    sampled_intervals = sum(
        left < stop and right > start
        for left, right in zip(history.times[:-1], history.times[1:])
    )
    blocks = min(requested_blocks, sampled_intervals)
    results: list[dict[str, object]] = []
    for quantity in quantities:
        if quantity not in history.values:
            raise ValueError(f"{history.label}: unknown quantity {quantity}")
        mean = history.average(quantity, start, stop)
        block_means = []
        for block in range(blocks):
            block_start = start+(stop-start)*block/blocks
            block_stop = start+(stop-start)*(block+1)/blocks
            block_means.append(history.average(quantity, block_start, block_stop))
        standard_error = (statistics.stdev(block_means)/math.sqrt(blocks)
                          if blocks > 1 else math.nan)
        results.append({
            "kind": kind,
            "label": history.label,
            "quantity": quantity,
            "tmin": start,
            "tmax": stop,
            "mean": mean,
            "block_se": standard_error,
            "nblocks": blocks,
        })
    return results


def difference_history(label: str, left: ForceHistory, right: ForceHistory,
                       quantities: list[str], start: float, stop: float) -> ForceHistory:
    overlap_times = {start, stop}
    overlap_times.update(time for time in left.times if start < time < stop)
    overlap_times.update(time for time in right.times if start < time < stop)
    result = ForceHistory.__new__(ForceHistory)
    result.label = label
    result.path = Path(f"{left.path}-{right.path}")
    result.times = sorted(overlap_times)
    result.values = {
        quantity: [left.value_at(quantity, time)-right.value_at(quantity, time)
                   for time in result.times]
        for quantity in quantities
    }
    return result


def parse_labeled_path(specification: str) -> tuple[str, Path]:
    if "=" not in specification:
        path = Path(specification)
        return path.stem, path
    label, path_text = specification.split("=", 1)
    if not label or not path_text:
        raise argparse.ArgumentTypeError("case must be LABEL=FILE or FILE")
    return label, Path(path_text)


def parse_difference(specification: str) -> tuple[str, str, str]:
    if "=" not in specification or "-" not in specification.split("=", 1)[1]:
        raise argparse.ArgumentTypeError("difference must be NAME=LEFT-RIGHT")
    name, operands = specification.split("=", 1)
    left, right = operands.split("-", 1)
    if not name or not left or not right:
        raise argparse.ArgumentTypeError("difference must be NAME=LEFT-RIGHT")
    return name, left, right


def print_results(results: list[dict[str, object]], output_format: str) -> None:
    fields = ("kind", "label", "quantity", "tmin", "tmax", "mean", "block_se",
              "nblocks")
    if output_format == "json":
        json_results = [
            {key: (None if isinstance(value, float) and not math.isfinite(value)
                   else value)
             for key, value in row.items()}
            for row in results
        ]
        print(json.dumps(json_results, indent=2, allow_nan=False))
        return
    if output_format == "csv":
        writer = csv.DictWriter(sys.stdout, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
        return
    print(f"{'type':8s} {'label':18s} {'quantity':16s} {'mean':>15s} "
          f"{'block SE':>15s} {'blocks':>6s}")
    for row in results:
        print(f"{row['kind']:8s} {row['label']:18s} {row['quantity']:16s} "
              f"{row['mean']:15.7e} {row['block_se']:15.7e} {row['nblocks']:6d}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="+", type=parse_labeled_path,
                        metavar="[LABEL=]FILE")
    parser.add_argument("--tmin", type=float, help="start of common averaging interval")
    parser.add_argument("--tmax", type=float, help="end of common averaging interval")
    parser.add_argument("--blocks", type=int, default=8,
                        help="number of equal-duration blocks for the standard error")
    parser.add_argument("--quantity", action="append", dest="quantities",
                        help="quantity to report; may be repeated")
    parser.add_argument("--difference", action="append", type=parse_difference,
                        default=[], metavar="NAME=LEFT-RIGHT")
    parser.add_argument("--no-standard-contrasts", action="store_true",
                        help="do not add full-frame_only, frame_only-isolated contrasts")
    parser.add_argument("--format", choices=("table", "csv", "json"), default="table")
    args = parser.parse_args()
    if args.blocks < 1:
        parser.error("--blocks must be positive")

    try:
        histories = {label: ForceHistory(label, path) for label, path in args.cases}
        if len(histories) != len(args.cases):
            raise ValueError("case labels must be unique")
        start = (args.tmin if args.tmin is not None
                 else max(history.times[0] for history in histories.values()))
        stop = (args.tmax if args.tmax is not None
                else min(history.times[-1] for history in histories.values()))
        if not math.isfinite(start) or not math.isfinite(stop) or not start < stop:
            raise ValueError("the cases have no positive common averaging interval")
        for history in histories.values():
            if start < history.times[0] or stop > history.times[-1]:
                raise ValueError(
                    f"averaging interval lies outside {history.label}'s time range"
                )
        quantities = args.quantities or default_quantities()
        results: list[dict[str, object]] = []
        for history in histories.values():
            results.extend(summarize(
                history, quantities, start, stop, args.blocks
            ))

        differences = list(args.difference)
        if not args.no_standard_contrasts:
            if {"full", "frame_only"}.issubset(histories):
                differences.append(("primary_metric_contrast", "full", "frame_only"))
            if {"frame_only", "isolated"}.issubset(histories):
                differences.append(("orbital_frame_contrast", "frame_only", "isolated"))
            if {"full", "isolated"}.issubset(histories):
                differences.append(("combined_contrast", "full", "isolated"))
        seen_differences: set[str] = set()
        for name, left_label, right_label in differences:
            if name in seen_differences:
                raise ValueError(f"duplicate difference label: {name}")
            seen_differences.add(name)
            if left_label not in histories or right_label not in histories:
                raise ValueError(f"unknown case in difference {name}")
            contrast = difference_history(
                name, histories[left_label], histories[right_label], quantities,
                start, stop
            )
            results.extend(summarize(
                contrast, quantities, start, stop, args.blocks, kind="contrast"
            ))
        print_results(results, args.format)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
