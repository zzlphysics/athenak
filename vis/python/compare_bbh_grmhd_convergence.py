#!/usr/bin/env python3
"""Time-align L/M/H BBH GRMHD histories and quantify tier differences.

The reported Richardson-like order is an empirical inner-hierarchy indicator.  The
campaign keeps the bulk-disk floor fixed at L4, so it is not a formal convergence order
for MRI turbulence or any diagnostic dominated by that unchanged disk resolution.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np


def read_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or "time" not in reader.fieldnames:
            raise RuntimeError(f"{path}: merged history CSV has no time column")
        values = {name: [] for name in reader.fieldnames}
        for row in reader:
            for name in reader.fieldnames:
                values[name].append(float(row[name]))
    result = {name: np.asarray(column) for name, column in values.items()}
    if result["time"].size < 2 or np.any(np.diff(result["time"]) <= 0.0):
        raise RuntimeError(f"{path}: time must be strictly increasing")
    return result


def derive(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    diagnostics: dict[str, np.ndarray] = {}
    if "mass" in data:
        diagnostics["mass_rel"] = data["mass"] / data["mass"][0] - 1.0
    if "tot-E" in data:
        diagnostics["total_energy_rel"] = data["tot-E"] / data["tot-E"][0] - 1.0
    magnetic = [data[name] for name in ("1-ME", "2-ME", "3-ME") if name in data]
    if magnetic:
        total = sum(magnetic)
        diagnostics["magnetic_energy_rel"] = total / total[0] - 1.0
    kinetic = [data[name] for name in ("1-KE", "2-KE", "3-KE") if name in data]
    if kinetic:
        total = sum(kinetic)
        diagnostics["kinetic_energy_rel"] = total / total[0] - 1.0
    if "dt" in data:
        diagnostics["root_dt"] = data["dt"]
    return diagnostics


def parse_series(values: list[str]) -> list[tuple[str, Path]]:
    series = []
    labels = set()
    for value in values:
        if "=" not in value:
            raise ValueError(f"series must be LABEL=merged-history.csv: {value}")
        label, path_text = value.split("=", 1)
        if not label or label in labels:
            raise ValueError(f"invalid or duplicate series label: {label!r}")
        labels.add(label)
        series.append((label, Path(path_text).expanduser().resolve(strict=True)))
    if len(series) < 2:
        raise ValueError("at least two series are required")
    return series


def l2(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def analyze(
    series: list[tuple[str, dict[str, np.ndarray]]],
    diagnostic_names: list[str],
    refinement_ratio: float,
) -> tuple[np.ndarray, dict[str, dict[str, object]], dict[str, dict[str, np.ndarray]]]:
    overlap_min = max(float(data["time"][0]) for _, data in series)
    overlap_max = min(float(data["time"][-1]) for _, data in series)
    if not overlap_max > overlap_min:
        raise RuntimeError("history series have no overlapping time interval")
    reference_time = series[-1][1]["time"]
    common_time = reference_time[
        (reference_time >= overlap_min) & (reference_time <= overlap_max)
    ]
    if common_time.size < 2:
        raise RuntimeError("fewer than two reference samples fall in the overlap")

    summaries: dict[str, dict[str, object]] = {}
    aligned: dict[str, dict[str, np.ndarray]] = {}
    for diagnostic in diagnostic_names:
        aligned[diagnostic] = {}
        for label, data in series:
            derived = derive(data)
            if diagnostic not in derived:
                raise RuntimeError(f"{label}: diagnostic unavailable: {diagnostic}")
            aligned[diagnostic][label] = np.interp(
                common_time, data["time"], derived[diagnostic]
            )
        pair_differences = []
        pair_norms = {}
        for (left_label, _), (right_label, _) in zip(series[:-1], series[1:]):
            difference = (
                aligned[diagnostic][left_label] - aligned[diagnostic][right_label]
            )
            pair_differences.append(difference)
            pair_norms[f"{left_label}-{right_label}"] = l2(difference)
        summary: dict[str, object] = {
            "pair_l2_absolute": pair_norms,
            "reference_l2": l2(aligned[diagnostic][series[-1][0]]),
        }
        if len(pair_differences) == 2:
            numerator = l2(pair_differences[0])
            denominator = l2(pair_differences[1])
            summary["empirical_order"] = (
                math.log(numerator / denominator) / math.log(refinement_ratio)
                if numerator > 0.0 and denominator > 0.0
                else None
            )
            summary["difference_ratio"] = (
                numerator / denominator if denominator > 0.0 else None
            )
        summaries[diagnostic] = summary
    return common_time, summaries, aligned


def render(
    output_path: Path,
    common_time: np.ndarray,
    series: list[tuple[str, dict[str, np.ndarray]]],
    diagnostic_names: list[str],
    aligned: dict[str, dict[str, np.ndarray]],
    dpi: int,
) -> None:
    cache = Path(tempfile.gettempdir()) / "athenak-matplotlib-cache"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib

    matplotlib.use("agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(
        len(diagnostic_names), 2, figsize=(12, 3.6 * len(diagnostic_names)), squeeze=False
    )
    for row, diagnostic in enumerate(diagnostic_names):
        positive_difference = False
        for label, _ in series:
            axes[row, 0].plot(common_time, aligned[diagnostic][label], label=label)
        for (left_label, _), (right_label, _) in zip(series[:-1], series[1:]):
            difference = np.abs(
                aligned[diagnostic][left_label] - aligned[diagnostic][right_label]
            )
            positive_difference = positive_difference or bool(np.any(difference > 0.0))
            axes[row, 1].plot(
                common_time, difference, label=f"|{left_label}-{right_label}|"
            )
        axes[row, 0].set_ylabel(diagnostic)
        axes[row, 1].set_ylabel("absolute tier difference")
        if positive_difference:
            axes[row, 1].set_yscale("log", nonpositive="clip")
        for axis in axes[row]:
            axis.set_xlabel("t / M")
            axis.grid(alpha=0.25)
            axis.legend(fontsize=8)
    figure.suptitle("BBH GRMHD time-aligned tier comparison")
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.98))
    figure.savefig(output_path, dpi=dpi)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("series", nargs="+", help="ordered LABEL=merged-history.csv")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--diagnostics",
        default="mass_rel,total_energy_rel,magnetic_energy_rel",
    )
    parser.add_argument("--refinement-ratio", type=float, default=2.0)
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()
    if args.refinement_ratio <= 1.0:
        raise SystemExit("--refinement-ratio must exceed 1")
    parsed = parse_series(args.series)
    loaded = [(label, read_csv(path)) for label, path in parsed]
    diagnostics = [name.strip() for name in args.diagnostics.split(",") if name.strip()]
    if not diagnostics:
        raise SystemExit("at least one diagnostic is required")
    common_time, summaries, aligned = analyze(
        loaded, diagnostics, args.refinement_ratio
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    aligned_path = args.output_dir / "time-aligned-diagnostics.csv"
    with aligned_path.open("w", newline="", encoding="utf-8") as stream:
        fieldnames = ["time"] + [
            f"{diagnostic}:{label}"
            for diagnostic in diagnostics
            for label, _ in loaded
        ]
        writer = csv.writer(stream)
        writer.writerow(fieldnames)
        for index, time_value in enumerate(common_time):
            writer.writerow(
                [time_value]
                + [
                    aligned[diagnostic][label][index]
                    for diagnostic in diagnostics
                    for label, _ in loaded
                ]
            )
    figure_path = args.output_dir / "convergence-history.png"
    render(figure_path, common_time, loaded, diagnostics, aligned, args.dpi)
    report = {
        "schema": 1,
        "classification": "athenak-bbh-grmhd-history-tier-comparison",
        "series": [
            {"label": label, "path": str(path)} for label, path in parsed
        ],
        "overlap_time_M": [float(common_time[0]), float(common_time[-1])],
        "samples": int(common_time.size),
        "refinement_ratio": args.refinement_ratio,
        "diagnostics": summaries,
        "aligned_csv": str(aligned_path),
        "figure": str(figure_path),
        "interpretation_limit": (
            "The L/M/H inner moving-hole spacing changes by 2:1, but the bulk disk "
            "floor remains L4. Empirical orders are not MRI-turbulence convergence "
            "orders."
        ),
    }
    report_path = args.output_dir / "convergence-report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(report_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
