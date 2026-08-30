#!/usr/bin/env python3
"""Plot and summarize the short r56 production-grid EMRI qualification run."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

from analyze_force_history import AXES, ForceHistory  # noqa: E402


SHELL_LABELS = (r"$0.5\,r_a$", r"$1\,r_a$", r"$2\,r_a$")
SHELL_COLORS = ("#3b82f6", "#f59e0b", "#dc2626")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history", type=Path, required=True)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--restart-summary", type=Path, required=True)
    parser.add_argument("--io-summary", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    return parser.parse_args()


def vector(history: ForceHistory, prefix: str) -> np.ndarray:
    return np.column_stack([
        history.values[history.force_name(prefix, axis)] for axis in AXES
    ])


def save_figure(figure: plt.Figure, base: Path) -> None:
    figure.savefig(base.with_suffix(".png"), dpi=220, bbox_inches="tight")
    figure.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    arguments = parse_arguments()
    output = arguments.output_directory.resolve()
    output.mkdir(parents=True, exist_ok=True)

    history = ForceHistory("r56-full", arguments.history.resolve(strict=True))
    pilot = json.loads(arguments.pilot.resolve(strict=True).read_text())
    restart = json.loads(arguments.restart_summary.resolve(strict=True).read_text())
    io_summary = json.loads(arguments.io_summary.resolve(strict=True).read_text())
    upstream = pilot["bhl_plan"]["upstream_state"]
    scales = pilot["bhl_plan"]["scales"]

    time = np.asarray(history.times)
    crossing_time = float(scales["capture_crossing_time"])
    capture_radius = float(scales["capture_radius_in_secondary_masses"])
    causal_speed = float(upstream["causal_crossing_speed_proxy"])
    causal_reach = causal_speed * time
    rho = float(upstream["rho"])
    effective_speed_squared = float(upstream["capture_effective_speed_squared"])
    mdot_proxy = 4.0 * math.pi * rho / effective_speed_squared**1.5
    mdot = np.asarray(history.values[history.mdot_name])

    wind = np.asarray(upstream["spatial_four_velocity"], dtype=float)
    wind_hat = wind / np.linalg.norm(wind)
    total = [vector(history, f"Ftotal{index}") for index in (1, 2, 3)]
    frel = [vector(history, f"Frel{index}") for index in (1, 2, 3)]
    delta_total = [values - values[0] for values in total]
    delta_parallel = [values @ wind_hat for values in delta_total]
    d21 = frel[1] - frel[0]
    d32 = frel[2] - frel[1]
    norm = lambda values: np.linalg.norm(values, axis=1)
    closure21 = norm(d21) / norm(total[1])
    closure32 = norm(d32) / norm(total[2])

    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "legend.fontsize": 8.5,
        "figure.titlesize": 13,
    })

    figure, axes = plt.subplots(2, 2, figsize=(11.2, 7.7), sharex=True)
    ax = axes[0, 0]
    ax.plot(time, 100.0 * mdot / mdot_proxy, "o-", color="#2563eb", lw=2)
    ax.set_ylabel(r"$100\,\dot M/\dot M_{\rm BHL,proxy}$  [%]")
    ax.set_title("(a) AMR startup + early accretion response")
    ax.grid(alpha=0.25)
    ax.annotate(
        f"final = {100.0*mdot[-1]/mdot_proxy:.2f}%\n"
        f"raw growth = {mdot[-1]/mdot[0]:.1f}x",
        xy=(time[-1], 100.0*mdot[-1]/mdot_proxy),
        xytext=(-92, -5), textcoords="offset points",
        arrowprops={"arrowstyle": "->", "color": "0.25"},
        bbox={"boxstyle": "round,pad=0.25", "fc": "white", "alpha": 0.85},
    )
    top = ax.secondary_xaxis(
        "top",
        functions=(lambda value: 100.0*value/crossing_time,
                   lambda value: value*crossing_time/100.0),
    )
    top.set_xlabel("capture-crossing fraction [%]")

    ax = axes[0, 1]
    for column, label, color in zip(range(3), ("radial x", "prograde y", "vertical z"),
                                    ("#2563eb", "#dc2626", "#059669")):
        ax.plot(time, total[2][:, column], "o-", color=color, label=label)
    ax.axhline(0.0, color="0.35", lw=0.8)
    ax.set_ylabel(r"$F_{\rm total}(2r_a)$ [code units]")
    ax.set_title("(b) Raw force estimator during AMR startup")
    ax.legend(ncol=1)
    ax.grid(alpha=0.25)

    ax = axes[1, 0]
    for values, label, color in zip(delta_parallel, SHELL_LABELS, SHELL_COLORS):
        ax.plot(time, values, "o-", color=color, label=label)
    ax.axhline(0.0, color="0.35", lw=0.8)
    ax.set_xlabel(r"$t/m_2$")
    ax.set_ylabel(r"$\Delta\mathbf{F}_{\rm total}\cdot\hat{u}_{\infty}$")
    ax.set_title("(c) Baseline difference (AMR + transient)")
    ax.legend()
    ax.grid(alpha=0.25)

    ax = axes[1, 1]
    ax.plot(time, 100.0*closure21, "o-", color=SHELL_COLORS[1],
            label=r"$|F_{\rm rel}(r_a)-F_{\rm rel}(0.5r_a)|/|F_{\rm total}(r_a)|$")
    ax.plot(time, 100.0*closure32, "o-", color=SHELL_COLORS[2],
            label=r"$|F_{\rm rel}(2r_a)-F_{\rm rel}(r_a)|/|F_{\rm total}(2r_a)|$")
    ax.axhline(15.0, color="black", ls="--", lw=1.1, label="provisional 15% gate")
    ax.set_xlabel(r"$t/m_2$")
    ax.set_ylabel("outer-shell increment [%]")
    ax.set_title("(d) Force-shell closure is not yet assessable")
    ax.legend(loc="best")
    ax.grid(alpha=0.25)
    figure.suptitle(
        "EMRI r56 production-grid qualification: early transient only\n"
        f"q=1e-5, r/M=56, t_end/t_cross={time[-1]/crossing_time:.4f}, "
        f"M_fast,proxy={upstream['fast_mach_proxy']:.3f}"
    )
    figure.tight_layout()
    save_figure(figure, output / "early_transient_physics")

    momentum = vector(history, "Fmom")
    newtonian = vector(history, "Fnewt")
    raw_terms = {
        r"$-F_{\rm mom}$": -momentum @ wind_hat,
        r"$F_{\rm Newt}$": newtonian @ wind_hat,
        r"$F_{\rm rel}(2r_a)$": frel[2] @ wind_hat,
        r"$F_{\rm total}(2r_a)$": total[2] @ wind_hat,
    }
    figure, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), sharex=True)
    for label, values in raw_terms.items():
        axes[0].plot(time, values, "o-", lw=1.8, label=label)
    axes[0].axhline(0.0, color="0.35", lw=0.8)
    axes[0].set_xlabel(r"$t/m_2$")
    axes[0].set_ylabel(r"$\mathbf{F}\cdot\hat{u}_{\infty}$ [code units]")
    axes[0].set_title("(a) Raw force decomposition")
    axes[0].legend()
    axes[0].grid(alpha=0.25)

    delta_terms = {
        r"$-\Delta F_{\rm mom}$": (-momentum + momentum[0]) @ wind_hat,
        r"$\Delta F_{\rm rel}(2r_a)$": (frel[2] - frel[2][0]) @ wind_hat,
        r"$\Delta F_{\rm total}(2r_a)$": delta_total[2] @ wind_hat,
    }
    for label, values in delta_terms.items():
        axes[1].plot(time, values, "o-", lw=1.8, label=label)
    axes[1].axhline(0.0, color="0.35", lw=0.8)
    axes[1].set_xlabel(r"$t/m_2$")
    axes[1].set_title("(b) Difference from the cycle-0 baseline")
    axes[1].legend()
    axes[1].grid(alpha=0.25)
    figure.suptitle(
        "Force terms projected on the recorded upstream velocity\n"
        "(baseline differences still contain AMR-quadrature changes)"
    )
    figure.tight_layout()
    save_figure(figure, output / "force_decomposition")

    scale_points = [
        ("momentum surface", 3.0, "#111827"),
        ("causal reach at t_end", float(causal_reach[-1]), "#7c3aed"),
        ("matching radius", float(scales["matching_radius"]), "#0f766e"),
        ("0.5 r_a", 0.5*capture_radius, SHELL_COLORS[0]),
        ("r_a", capture_radius, SHELL_COLORS[1]),
        ("2 r_a", 2.0*capture_radius, SHELL_COLORS[2]),
        ("magnetic coherence", float(scales["limiting_environment_scale"]), "#64748b"),
        ("Hill radius", float(scales["hill_radius"]), "#334155"),
    ]
    figure, (ax, text_ax) = plt.subplots(
        2, 1, figsize=(11.2, 5.4), gridspec_kw={"height_ratios": [1.25, 1.0]}
    )
    ax.set_xscale("log")
    for row, (label, value, color) in enumerate(scale_points):
        y = row % 2
        ax.scatter(value, y, s=65, color=color, zorder=3)
        ax.vlines(value, -0.2, y, color=color, lw=1, alpha=0.55)
        ax.annotate(f"{label}\n{value:.3g} m2", (value, y),
                    xytext=(4, 8 if y == 0 else -28), textcoords="offset points",
                    color=color, fontsize=8.5)
    ax.axvspan(3.0, float(causal_reach[-1]), color="#7c3aed", alpha=0.10,
               label="region causally reached by end of run")
    ax.set_ylim(-0.35, 1.35)
    ax.set_yticks([])
    ax.set_xlabel(r"secondary-frame length scale [$m_2$]")
    ax.set_title("Spatial-scale separation at the final saved time")
    ax.grid(axis="x", which="both", alpha=0.2)
    ax.legend(loc="upper right")

    endpoint = restart["endpoint_comparison"]
    cache = restart["cold_restart"]
    lines = [
        (
            f"Upstream: beta={upstream['plasma_beta']:.1f}, "
            f"sigma={upstream['magnetization_sigma']:.2e}, "
            f"v={upstream['three_speed']:.4f}, "
            f"M_fast,proxy={upstream['fast_mach_proxy']:.3f}"
        ),
        (
            f"Time coverage: t_end={time[-1]:.3f} m2 = "
            f"{100.0*time[-1]/crossing_time:.3f}% of one capture crossing; "
            f"causal reach={causal_reach[-1]:.3f} m2"
        ),
        (
            "Restart identity: MHD/CT/ADM endpoint max differences = "
            f"{endpoint['mhd_active_max_abs']:.0f}/"
            f"{endpoint['face_active_max_abs']:.0f}/"
            f"{endpoint['adm_active_max_abs']:.0f}; history max difference = "
            f"{restart['history_comparison']['maximum_absolute_difference']:.0f}"
        ),
        (
            f"Numerics: max|divB|={restart['constraint']['continuous_maximum_divB']:.2e}, "
            f"cache proposal/tolerance="
            f"{cache['maximum_mixed_scale_proposed_change']/cache['mixed_scale_acceptance_tolerance']:.2e}, "
            f"peak A100 memory={io_summary['peak_gpu_memory_MiB']:.0f} MiB"
        ),
    ]
    text_ax.axis("off")
    text_ax.text(0.01, 0.95, "\n\n".join(lines), va="top", ha="left",
                 family="monospace", fontsize=9.2,
                 bbox={"boxstyle": "round,pad=0.55", "fc": "#f8fafc", "ec": "#cbd5e1"})
    figure.tight_layout()
    save_figure(figure, output / "scale_and_numerical_health")

    cycle_rows = restart["phases"]["continuous_fresh"]["cycles"]
    cycles = np.asarray([row["cycle"] for row in cycle_rows[1:]], dtype=int)
    elapsed = np.asarray([row["elapsed"] for row in cycle_rows], dtype=float)
    cycle_wall = np.diff(elapsed)
    figure, ax = plt.subplots(figsize=(8.4, 4.5))
    ax.semilogy(cycles, cycle_wall, "o-", color="#7c3aed", lw=2)
    ax.axhspan(108.0, 121.0, color="#f59e0b", alpha=0.16,
               label="full-topology cycle range")
    ax.set_xlabel("completed root cycle")
    ax.set_ylabel("wall time per root cycle [s]")
    ax.set_title("AMR startup cost converges to the production topology")
    ax.grid(which="both", alpha=0.25)
    ax.legend()
    figure.tight_layout()
    save_figure(figure, output / "runtime_amr_ramp")

    csv_columns = ["time", "crossing_fraction", "causal_reach", "mdot_hat",
                   "mdot_proxy_fraction"]
    for shell in (1, 2, 3):
        csv_columns.extend(f"Ftotal{shell}H_{axis}" for axis in AXES)
        csv_columns.append(f"delta_Ftotal{shell}_parallel")
    csv_columns.extend(("closure21", "closure32"))
    with (output / "derived_history.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=csv_columns)
        writer.writeheader()
        for index, value in enumerate(time):
            row: dict[str, float] = {
                "time": float(value),
                "crossing_fraction": float(value/crossing_time),
                "causal_reach": float(causal_reach[index]),
                "mdot_hat": float(mdot[index]),
                "mdot_proxy_fraction": float(mdot[index]/mdot_proxy),
                "closure21": float(closure21[index]),
                "closure32": float(closure32[index]),
            }
            for shell in (1, 2, 3):
                for axis_index, axis in enumerate(AXES):
                    row[f"Ftotal{shell}H_{axis}"] = float(
                        total[shell-1][index, axis_index]
                    )
                row[f"delta_Ftotal{shell}_parallel"] = float(
                    delta_parallel[shell-1][index]
                )
            writer.writerow(row)

    delta32 = d32[-1] - d32[0]
    summary = {
        "classification": "athenak-emri-r56-early-transient-analysis-v1",
        "case_id": pilot["case_id"],
        "scope": "early numerical and causal transient; not a stationary BHL result",
        "samples": len(time),
        "time_end_in_secondary_masses": float(time[-1]),
        "capture_crossing_time_in_secondary_masses": crossing_time,
        "capture_crossing_fraction": float(time[-1]/crossing_time),
        "causal_reach_in_secondary_masses": float(causal_reach[-1]),
        "upstream": {
            "rho": rho,
            "plasma_beta": upstream["plasma_beta"],
            "magnetization_sigma": upstream["magnetization_sigma"],
            "three_speed": upstream["three_speed"],
            "fast_mach_proxy": upstream["fast_mach_proxy"],
            "wind_unit_vector": wind_hat.tolist(),
        },
        "accretion": {
            "mdot_initial": float(mdot[0]),
            "mdot_final": float(mdot[-1]),
            "growth_factor": float(mdot[-1]/mdot[0]),
            "bhl_proxy": mdot_proxy,
            "final_to_bhl_proxy": float(mdot[-1]/mdot_proxy),
            "fractional_change_last_saved_interval": float(
                mdot[-1]/mdot[-2] - 1.0
            ),
        },
        "force": {
            "Ftotal3_final": total[2][-1].tolist(),
            "Ftotal3_vertical_fraction": float(
                abs(total[2][-1, 2])/np.linalg.norm(total[2][-1])
            ),
            "outer_increment_fraction_at_end": float(closure32[-1]),
            "dFrel32_change_relative_to_initial": float(
                np.linalg.norm(delta32)/np.linalg.norm(d32[0])
            ),
            "delta_Ftotal_final_by_shell": [values[-1].tolist()
                                             for values in delta_total],
            "Ftotal3_norm_fractional_change_last_saved_interval": float(
                np.linalg.norm(total[2][-1])/np.linalg.norm(total[2][-2]) - 1.0
            ),
        },
        "numerical_health": {
            "endpoint_all_stored_fields_match": endpoint["all_stored_fields_match"],
            "history_exact_after_time_deduplication": restart[
                "history_comparison"
            ]["exact_after_time_deduplication"],
            "maximum_divB": restart["constraint"]["continuous_maximum_divB"],
            "cache_proposal_to_tolerance": (
                cache["maximum_mixed_scale_proposed_change"]
                / cache["mixed_scale_acceptance_tolerance"]
            ),
            "peak_gpu_memory_MiB": io_summary["peak_gpu_memory_MiB"],
        },
    }
    (output / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
