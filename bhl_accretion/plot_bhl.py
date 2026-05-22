#!/usr/bin/env python3
"""Quick-look plots for AthenaK GR BHL runs.

This script is intentionally a thin wrapper around AthenaK's existing
``vis/python/plot_slice.py`` for AMR-aware binary slices, plus a small history
plotter for the user history written by ``src/pgen/gr_bhl.cpp``.
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path
import subprocess
import sys
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
VIS_DIR = REPO_ROOT / "vis" / "python"


def import_athena_read():
    sys.path.insert(0, str(VIS_DIR))
    import athena_read  # pylint: disable=import-error,import-outside-toplevel

    return athena_read


def find_history(run_dir: Path, basename: str | None) -> Path | None:
    if basename:
        candidates = sorted(run_dir.glob(f"{basename}*.user.hst"))
        candidates += sorted(run_dir.glob(f"{basename}*.hst"))
    else:
        candidates = sorted(run_dir.glob("*.user.hst"))
        candidates += sorted(run_dir.glob("*.hst"))
    return candidates[-1] if candidates else None


def find_dumps(run_dir: Path, basename: str | None) -> list[Path]:
    bin_dir = run_dir / "bin"
    if basename:
        patterns = [
            bin_dir / f"{basename}.*.bin",
            run_dir / f"{basename}.*.bin",
        ]
    else:
        patterns = [bin_dir / "*.bin", run_dir / "*.bin"]
    dumps: list[Path] = []
    for pattern in patterns:
        dumps.extend(Path(p) for p in glob.glob(str(pattern)))
    return sorted(set(dumps))


def require_numpy():
    try:
        import numpy as np  # pylint: disable=import-outside-toplevel
    except ImportError as exc:
        raise SystemExit(
            "numpy is required for history plotting. Install it in the active Python "
            "environment, e.g. add py-numpy to the grmhd Spack environment."
        ) from exc

    return np


def require_pyplot():
    try:
        import matplotlib  # pylint: disable=import-outside-toplevel

        matplotlib.use("agg")
        import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required for plotting. Install it in the active Python "
            "environment, e.g. add py-matplotlib to the grmhd Spack environment."
        ) from exc

    return plt


def moving_average(values, width: int):
    np = require_numpy()
    if width <= 1 or len(values) < width:
        return values
    kernel = np.ones(width) / width
    return np.convolve(values, kernel, mode="same")


def first_existing_key(data: dict, prefixes: Iterable[str]) -> str | None:
    for prefix in prefixes:
        for key in data:
            if key.startswith(prefix):
                return key
    return None


def keys_with_prefix(data: dict, prefixes: Iterable[str]) -> list[str]:
    return [key for key in data if any(key.startswith(prefix) for prefix in prefixes)]


def plot_history(
    hst_file: Path,
    output_dir: Path,
    r_acc: float,
    v_inf: float,
    smooth: int,
) -> None:
    np = require_numpy()
    plt = require_pyplot()
    athena_read = import_athena_read()
    data = athena_read.hst(str(hst_file))
    time = data["time"]
    tau_a = r_acc / v_inf
    x = time / tau_a

    mdot_key = first_existing_key(data, ("mdotHL_",))
    mdot_raw_key = first_existing_key(data, ("mdot_",))
    varphi_key = first_existing_key(data, ("varphi_",))
    phi_key = first_existing_key(data, ("phi_",))
    edot_key = first_existing_key(data, ("edot_",))
    jdot_key = first_existing_key(data, ("jdot_",))
    drag_keys = keys_with_prefix(data, ("fd_x_", "fd_y_", "fd_z_"))

    panel_count = 1
    panel_count += 1 if (varphi_key or phi_key) else 0
    panel_count += 1 if (edot_key or jdot_key) else 0
    panel_count += 1 if drag_keys else 0
    fig, axes = plt.subplots(panel_count, 1, figsize=(8.0, 2.2 + 1.6*panel_count),
                             sharex=True)
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    axis_iter = iter(axes)

    ax = next(axis_iter)
    if mdot_key:
        y = data[mdot_key]
        ax.plot(x, y, color="0.65", linewidth=0.8, label=mdot_key)
        ax.plot(x, moving_average(y, smooth), color="tab:blue", linewidth=1.5,
                label=f"{mdot_key}, smoothed")
        ax.set_ylabel(r"$\dot{M}/\dot{M}_{HL}$")
    elif mdot_raw_key:
        y = data[mdot_raw_key]
        ax.plot(x, y, color="tab:blue", linewidth=1.2, label=mdot_raw_key)
        ax.set_ylabel(r"$\dot{M}$")
    else:
        ax.text(0.5, 0.5, "No mdot history column found", ha="center", va="center")
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8)

    if varphi_key or phi_key:
        ax = next(axis_iter)
        key = varphi_key or phi_key
        y = data[key]
        ax.plot(x, y, color="0.65", linewidth=0.8, label=key)
        ax.plot(x, moving_average(y, smooth), color="tab:orange", linewidth=1.5,
                label=f"{key}, smoothed")
        ax.set_ylabel(r"$\varphi$" if varphi_key else r"$\Phi$")
        ax.grid(alpha=0.25)
        ax.legend(loc="best", fontsize=8)

    if edot_key or jdot_key:
        ax = next(axis_iter)
        if edot_key:
            ax.plot(x, moving_average(data[edot_key], smooth), color="tab:green",
                    linewidth=1.3, label=edot_key)
        if jdot_key:
            ax.plot(x, moving_average(data[jdot_key], smooth), color="tab:red",
                    linewidth=1.3, label=jdot_key)
        ax.set_ylabel("flux")
        ax.grid(alpha=0.25)
        ax.legend(loc="best", fontsize=8)

    if drag_keys:
        ax = next(axis_iter)
        colors = {"x": "tab:blue", "y": "tab:orange", "z": "tab:green"}
        for key in drag_keys:
            component = key.split("_")[1]
            ax.plot(x, moving_average(data[key], smooth), linewidth=1.2,
                    color=colors.get(component), label=key)
        ax.axhline(0.0, color="0.2", linewidth=0.7, alpha=0.5)
        ax.set_ylabel(r"$F^i_{drag}$")
        ax.grid(alpha=0.25)
        ax.legend(loc="best", fontsize=8)

    axes[-1].set_xlabel(r"$t/\tau_a$")
    fig.suptitle(hst_file.name)
    fig.tight_layout()
    fig.savefig(output_dir / "history_bhl.png", dpi=180)
    plt.close(fig)


def run_slice_plot(
    dump: Path,
    variable: str,
    output_file: Path,
    dimension: str,
    r_max: float,
    norm: str | None,
    cmap: str,
    dpi: int,
    extra: list[str],
) -> None:
    cmd = [
        sys.executable,
        str(VIS_DIR / "plot_slice.py"),
        str(dump),
        variable,
        str(output_file),
        "-d",
        dimension,
        "-l",
        "0.0",
        "--r_max",
        str(r_max),
        "-c",
        cmap,
        "--horizon",
        "--horizon_mask",
        "--notex",
        "--dpi",
        str(dpi),
    ]
    if norm:
        cmd += ["-n", norm]
    cmd += extra
    subprocess.run(cmd, check=True)


def make_slice_gallery(
    dumps: list[Path],
    output_dir: Path,
    variables: list[str],
    r_max: float,
    latest_only: bool,
    stride: int,
    dpi: int,
) -> None:
    if latest_only and dumps:
        dumps = [dumps[-1]]
    elif stride > 1:
        dumps = dumps[::stride]

    for dump in dumps:
        stem = dump.stem
        for variable in variables:
            safe_var = variable.replace("derived:", "derived_").replace("/", "_")
            norm = "log" if variable in ("dens", "derived:pgas", "derived:sigma_rel") else None
            cmap = "magma" if variable == "dens" else "viridis"

            for dimension, plane in (("y", "xz"), ("z", "xy")):
                output_file = output_dir / f"{stem}.{plane}.{safe_var}.png"
                try:
                    run_slice_plot(
                        dump=dump,
                        variable=variable,
                        output_file=output_file,
                        dimension=dimension,
                        r_max=r_max,
                        norm=norm,
                        cmap=cmap,
                        dpi=dpi,
                        extra=[],
                    )
                except subprocess.CalledProcessError as exc:
                    print(f"warning: failed plotting {variable} from {dump.name}: {exc}",
                          file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create quick-look BHL history and slice plots from AthenaK outputs."
    )
    parser.add_argument("--run-dir", default=".", help="directory containing hst and bin/")
    parser.add_argument("--basename", help="AthenaK job basename, e.g. bhl_b10nr53")
    parser.add_argument("--output-dir", default="bhl_plots", help="plot output directory")
    parser.add_argument("--r-acc", type=float, default=200.0,
                        help="BHL accretion radius used for t/tau_a")
    parser.add_argument("--v-inf", type=float, default=0.1,
                        help="wind speed used for t/tau_a")
    parser.add_argument("--r-max", type=float, default=300.0,
                        help="half-width for slice plots")
    parser.add_argument("--smooth", type=int, default=25,
                        help="history moving-average window in samples")
    parser.add_argument("--variables", default="dens,derived:pgas,derived:beta_inv_rel",
                        help="comma-separated variables for plot_slice.py")
    parser.add_argument("--latest-only", action="store_true",
                        help="plot only the latest binary dump")
    parser.add_argument("--stride", type=int, default=1,
                        help="plot every Nth dump when not using --latest-only")
    parser.add_argument("--dpi", type=int, default=180, help="image DPI")
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = run_dir / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    hst_file = find_history(run_dir, args.basename)
    if hst_file:
        plot_history(hst_file, output_dir, args.r_acc, args.v_inf, args.smooth)
        print(f"wrote {output_dir / 'history_bhl.png'}")
    else:
        print("warning: no history file found", file=sys.stderr)

    dumps = find_dumps(run_dir, args.basename)
    if dumps:
        variables = [item.strip() for item in args.variables.split(",") if item.strip()]
        make_slice_gallery(dumps, output_dir, variables, args.r_max, args.latest_only,
                           args.stride, args.dpi)
        print(f"wrote slice plots to {output_dir}")
    else:
        print("warning: no binary dumps found", file=sys.stderr)


if __name__ == "__main__":
    main()
