#!/usr/bin/env python3
"""Plot the EMRI delta-T force zero test and the r56 horizon resolution audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DOMAIN_LOWER = np.array(
    [-2890.2711304431523, -2645.0832268188583, -1318.1407462145046]
)
DOMAIN_UPPER = np.array(
    [1845.119124744866, 1845.1185514675913, 1320.1180539362786]
)
ROOT_DIMENSIONS = np.array([128, 112, 80])
PHYSICAL_REFINEMENT_LEVELS = 9
SECONDARY_HORIZON_RADIUS = 2.0  # chi_2=0, r_+=2m_2
ATHENAK_PAPER_FINE_DX = np.array([0.125, 0.08333333333333333, 0.0625])
ATHENAK_PAPER_SPIN = 0.9375
ATHENAK_PAPER_HORIZON_RADIUS = 1.0 + np.sqrt(1.0 - ATHENAK_PAPER_SPIN**2)


def _last_row(path: Path) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(path))[-1]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subtracted-history", type=Path, required=True)
    parser.add_argument("--raw-history", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()

    subtracted = _last_row(args.subtracted_history)
    raw = _last_row(args.raw_history)
    frel_slice = slice(12, 21)
    fmom_slice = slice(6, 9)
    labels = [f"R{k}{axis}" for k in (1, 2, 3) for axis in "xyz"]

    root_dx = (DOMAIN_UPPER - DOMAIN_LOWER) / ROOT_DIMENSIONS
    finest_dx = root_dx / 2**PHYSICAL_REFINEMENT_LEVELS
    cells_per_horizon_radius = SECONDARY_HORIZON_RADIUS / finest_dx

    summary = {
        "classification": "athenak-emri-force-background-and-resolution-audit-v1",
        "force_zero_test": {
            "scope": "local CPU t=0 analytic Taylor-gradient smoke test",
            "raw_Frel_vectors": raw[frel_slice].reshape(3, 3).tolist(),
            "deltaT_Frel_vectors": subtracted[frel_slice].reshape(3, 3).tolist(),
            "raw_Frel_norm": float(np.linalg.norm(raw[frel_slice])),
            "deltaT_Frel_norm": float(np.linalg.norm(subtracted[frel_slice])),
            "maximum_absolute_deltaT_Frel_component": float(
                np.max(np.abs(subtracted[frel_slice]))
            ),
            "raw_Fmom": raw[fmom_slice].tolist(),
            "analytic_subtracted_Fmom": subtracted[fmom_slice].tolist(),
            "note": (
                "The remaining Fmom is interpolation/quadrature bias on the extraction "
                "sphere; measure it with the fixed-topology restart control."
            ),
        },
        "r56_resolution": {
            "root_dx_in_secondary_masses": root_dx.tolist(),
            "finest_dx_in_secondary_masses": finest_dx.tolist(),
            "physical_refinement_levels": PHYSICAL_REFINEMENT_LEVELS,
            "secondary_spin": 0.0,
            "horizon_radius_in_secondary_masses": SECONDARY_HORIZON_RADIUS,
            "cells_per_horizon_radius_by_axis": cells_per_horizon_radius.tolist(),
            "cells_per_horizon_diameter_by_axis": (
                2.0 * cells_per_horizon_radius
            ).tolist(),
            "cell_diagonal_in_secondary_masses": float(
                np.linalg.norm(finest_dx)
            ),
            "athenak_paper_fine_dx_in_black_hole_masses":
                ATHENAK_PAPER_FINE_DX.tolist(),
            "athenak_paper_spin": ATHENAK_PAPER_SPIN,
            "athenak_paper_cells_per_horizon_radius": (
                ATHENAK_PAPER_HORIZON_RADIUS / ATHENAK_PAPER_FINE_DX
            ).tolist(),
            "assessment": (
                "Current finest spacing is comparable to the paper's 96^3--128^3 "
                "near-horizon range, but magnetic flux/force convergence is not "
                "established without a one-level-finer local rerun."
            ),
        },
    }

    output = args.output_directory
    output.mkdir(parents=True, exist_ok=True)
    (output / "validation_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )

    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.5), constrained_layout=True)

    x = np.arange(len(labels))
    width = 0.38
    axes[0].bar(x - width / 2, np.abs(raw[frel_slice]), width, label="raw $T^{\\mu\\nu}$")
    axes[0].bar(
        x + width / 2,
        np.abs(subtracted[frel_slice]),
        width,
        label="$\\delta T^{\\mu\\nu}$",
    )
    axes[0].set_yscale("log")
    axes[0].set_xticks(x, labels, rotation=45)
    axes[0].set_ylabel("absolute force component")
    axes[0].set_title("Analytic-wind force zero test")
    axes[0].legend(frameon=False)

    paper_labels = ["AthenaK 64³", "AthenaK 96³", "AthenaK 128³"]
    current_labels = ["current x", "current y", "current z"]
    paper_cells_per_horizon = ATHENAK_PAPER_HORIZON_RADIUS / ATHENAK_PAPER_FINE_DX
    all_values = np.concatenate([paper_cells_per_horizon, cells_per_horizon_radius])
    colors = ["#999999"] * 3 + ["#2474b5"] * 3
    axes[1].barh(paper_labels + current_labels, all_values, color=colors)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("cells per Kerr horizon radius $r_+$")
    axes[1].set_title("Near-horizon spatial resolution")
    for index, value in enumerate(all_values):
        suffix = (
            f"  ($\\Delta x={ATHENAK_PAPER_FINE_DX[index]:.4f}M$)"
            if index < 3
            else f"  ($\\Delta x={finest_dx[index-3]:.4f}m_2$)"
        )
        axes[1].text(value + 0.4, index, f"{value:.1f}{suffix}", va="center", fontsize=9)
    axes[1].set_xlim(0.0, 40.0)

    fig.savefig(output / "force_background_and_resolution.png", dpi=180)
    fig.savefig(output / "force_background_and_resolution.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
