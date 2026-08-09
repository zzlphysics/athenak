#!/usr/bin/env python3
"""Create AMR-aware BBH GRMHD dashboards from AthenaK binary output.

The stock ``plot_slice.py`` remains the general single-panel plotting tool.  This
script follows its binary-slice selection logic, but reads every requested field in
one pass so that multi-panel figures remain practical for gigabyte-scale AMR dumps.

The current ``mhd_w_bcc`` output does not contain the dynamical ADM metric.  Density,
pressure, temperature, magnetic components, and primitive velocity components are
therefore plotted exactly as stored.  Quantities involving contractions of magnetic
or velocity components are explicitly labelled ``proxy`` and must not be interpreted
as covariant GR diagnostics.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import glob
import json
import math
import os
from pathlib import Path
import struct
import sys
import tempfile
from typing import Callable

import numpy as np


PROXY_WARNING = (
    "Coordinate-component proxy: the mhd_w_bcc dump does not contain the "
    "dynamical ADM metric needed for a covariant GR contraction."
)


@dataclass(frozen=True)
class BinaryHeader:
    time: float
    cycle: int
    location_size: int
    variable_size: int
    variables: tuple[str, ...]
    parameters: dict[str, dict[str, str]]
    data_offset: int


@dataclass
class SliceBlock:
    extent: tuple[float, float, float, float]
    level: int
    logical_location: tuple[int, int, int]
    slice_shape: tuple[int, int]
    fields: dict[str, np.ndarray]


@dataclass
class SliceData:
    header: BinaryHeader
    plane: str
    location: float
    blocks: list[SliceBlock]
    level_counts: Counter
    selected_level_counts: Counter


@dataclass(frozen=True)
class Panel:
    name: str
    label: str
    dependencies: tuple[str, ...]
    calculate: Callable[[dict[str, np.ndarray], int], np.ndarray]
    cmap: str
    scale: str
    proxy: bool = False


def _identity(name: str) -> Callable[[dict[str, np.ndarray], int], np.ndarray]:
    return lambda fields, level: fields[name]


def _b2(fields: dict[str, np.ndarray]) -> np.ndarray:
    return fields["bcc1"] ** 2 + fields["bcc2"] ** 2 + fields["bcc3"] ** 2


PANELS = {
    "dens": Panel("dens", r"$\rho$", ("dens",), _identity("dens"), "magma", "log"),
    "press": Panel(
        "press", r"$p_\mathrm{gas}$", ("press",), _identity("press"), "inferno", "log"
    ),
    "temperature": Panel(
        "temperature",
        r"$T$ (code output)",
        ("temperature",),
        _identity("temperature"),
        "plasma",
        "log",
    ),
    "bmag": Panel(
        "bmag",
        r"$\sqrt{(B^x)^2+(B^y)^2+(B^z)^2}$ (proxy)",
        ("dens", "bcc1", "bcc2", "bcc3"),
        lambda fields, level: np.sqrt(_b2(fields)),
        "viridis",
        "log",
        True,
    ),
    "beta_inv_proxy": Panel(
        "beta_inv_proxy",
        r"$[(B^x)^2+(B^y)^2+(B^z)^2]/(2p)$ (proxy)",
        ("dens", "press", "bcc1", "bcc2", "bcc3"),
        lambda fields, level: np.divide(
            0.5 * _b2(fields),
            fields["press"],
            out=np.full_like(fields["press"], np.nan),
            where=fields["press"] > 0.0,
        ),
        "cividis",
        "log",
        True,
    ),
    "sigma_proxy": Panel(
        "sigma_proxy",
        r"$[(B^x)^2+(B^y)^2+(B^z)^2]/\rho$ (proxy)",
        ("dens", "bcc1", "bcc2", "bcc3"),
        lambda fields, level: np.divide(
            _b2(fields),
            fields["dens"],
            out=np.full_like(fields["dens"], np.nan),
            where=fields["dens"] > 0.0,
        ),
        "cubehelix",
        "log",
        True,
    ),
    "velmag_proxy": Panel(
        "velmag_proxy",
        r"$\sqrt{(\tilde u^x)^2+(\tilde u^y)^2+(\tilde u^z)^2}$ (proxy)",
        ("dens", "velx", "vely", "velz"),
        lambda fields, level: np.sqrt(
            fields["velx"] ** 2 + fields["vely"] ** 2 + fields["velz"] ** 2
        ),
        "magma",
        "linear",
        True,
    ),
    "divb": Panel(
        "divb",
        r"$\nabla\!\cdot\!B$",
        ("divb",),
        _identity("divb"),
        "coolwarm",
        "symlog",
    ),
    "level": Panel(
        "level",
        "physical AMR level",
        (),
        lambda fields, level: np.full(next(iter(fields.values())).shape, level)
        if fields
        else np.asarray(level),
        "turbo",
        "level",
    ),
}


def parse_parameter_header(text: str) -> dict[str, dict[str, str]]:
    """Parse the athinput copy embedded in an AthenaK binary header."""

    parameters: dict[str, dict[str, str]] = {}
    section: str | None = None
    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            section = line[1:-1]
            parameters.setdefault(section, {})
            continue
        if section is not None and "=" in line:
            key, value = line.split("=", 1)
            parameters[section][key.strip()] = value.strip()
    return parameters


def read_binary_header(path: Path) -> BinaryHeader:
    """Read the AthenaK binary preheader and embedded input parameters."""

    with path.open("rb") as stream:
        signature = stream.readline().decode("ascii")
        if signature != "Athena binary output version=1.1\n":
            raise RuntimeError(f"{path}: unsupported Athena binary signature")
        preheader_count = int(stream.readline().split(b"=", 1)[1])
        preheader = {}
        for _ in range(preheader_count - 1):
            key, value = stream.readline().decode("ascii").split("=", 1)
            preheader[key.strip()] = value.strip()
        number_variables = int(stream.readline().split(b"=", 1)[1])
        variable_line = stream.readline().decode("ascii").split(":", 1)[1]
        variables = tuple(variable_line.split())
        if len(variables) != number_variables:
            raise RuntimeError(
                f"{path}: header declares {number_variables} variables but lists "
                f"{len(variables)}"
            )
        header_size = int(stream.readline().split(b"=", 1)[1])
        parameter_text = stream.read(header_size).decode("ascii", "replace")
        return BinaryHeader(
            time=float(preheader["time"]),
            cycle=int(preheader["cycle"]),
            location_size=int(preheader["size of location"]),
            variable_size=int(preheader["size of variable"]),
            variables=variables,
            parameters=parse_parameter_header(parameter_text),
            data_offset=stream.tell(),
        )


def _get_int(header: BinaryHeader, section: str, key: str) -> int:
    try:
        return int(header.parameters[section][key])
    except KeyError as exc:
        raise RuntimeError(f"binary header has no <{section}>/{key}") from exc


def _get_float(header: BinaryHeader, section: str, key: str) -> float:
    try:
        return float(header.parameters[section][key])
    except KeyError as exc:
        raise RuntimeError(f"binary header has no <{section}>/{key}") from exc


def _target_block_and_index(
    location: float,
    domain_min: float,
    domain_max: float,
    root_blocks: int,
    cells_per_block: int,
    level: int,
) -> tuple[int, int]:
    if location <= domain_min:
        return 0, 0
    blocks_at_level = root_blocks * 2**level
    if location >= domain_max:
        return blocks_at_level - 1, cells_per_block - 1
    normalized = (location - domain_min) / (domain_max - domain_min)
    number_cells = cells_per_block * blocks_at_level
    cell = min(int(normalized * number_cells), number_cells - 1)
    return cell // cells_per_block, cell % cells_per_block


def _overlaps_extent(
    extent: tuple[float, float, float, float], half_width: float | None
) -> bool:
    if half_width is None:
        return True
    return not (
        extent[1] <= -half_width
        or extent[0] >= half_width
        or extent[3] <= -half_width
        or extent[2] >= half_width
    )


def read_slice(
    path: Path,
    panel_names: list[str],
    plane: str,
    location: float,
    half_width: float | None,
) -> SliceData:
    """Read all requested fields on one AMR slice in a single file scan."""

    header = read_binary_header(path)
    if header.location_size not in (4, 8) or header.variable_size not in (4, 8):
        raise RuntimeError(
            "only 4- and 8-byte AthenaK location/variable data are supported"
        )

    dependencies = {
        dependency
        for panel_name in panel_names
        for dependency in PANELS[panel_name].dependencies
    }
    missing = sorted(dependencies.difference(header.variables))
    if missing:
        raise RuntimeError(
            f"{path}: panels require missing variables: {', '.join(missing)}; "
            f"available: {', '.join(header.variables)}"
        )
    variable_indices = {name: header.variables.index(name) for name in dependencies}
    nghost = _get_int(header, "mesh", "nghost")
    axis = {"x": 0, "y": 1, "z": 2}[plane]
    mesh_cells = tuple(_get_int(header, "mesh", f"nx{i}") for i in (1, 2, 3))
    block_cells_config = tuple(
        _get_int(header, "meshblock", f"nx{i}") for i in (1, 2, 3)
    )
    domain_min = _get_float(header, "mesh", f"x{axis + 1}min")
    domain_max = _get_float(header, "mesh", f"x{axis + 1}max")
    root_blocks = mesh_cells[axis] // block_cells_config[axis]

    location_dtype = np.dtype("<f4" if header.location_size == 4 else "<f8")
    variable_dtype = np.dtype("<f4" if header.variable_size == 4 else "<f8")
    level_counts: Counter = Counter()
    selected_level_counts: Counter = Counter()
    blocks: list[SliceBlock] = []

    with path.open("rb") as stream:
        stream.seek(header.data_offset)
        file_size = path.stat().st_size
        first_record = True
        while stream.tell() < file_size:
            index_bytes = stream.read(24)
            if not index_bytes:
                break
            if len(index_bytes) != 24:
                raise RuntimeError(f"{path}: truncated MeshBlock index record")
            indices = np.asarray(struct.unpack("=6i", index_bytes)) - nghost
            logical = struct.unpack("=4i", stream.read(16))
            block_i, block_j, block_k, level = logical
            level_counts[level] += 1
            block_cells = (
                int(indices[1] - indices[0] + 1),
                int(indices[3] - indices[2] + 1),
                int(indices[5] - indices[4] + 1),
            )
            if first_record:
                if block_cells != block_cells_config:
                    raise RuntimeError(
                        f"{path}: output block shape {block_cells} differs from "
                        f"configured shape {block_cells_config}"
                    )
                first_record = False
            number_cells = math.prod(block_cells)
            variable_bytes = number_cells * header.variable_size
            record_data_bytes = len(header.variables) * variable_bytes
            target_block, slice_index = _target_block_and_index(
                location,
                domain_min,
                domain_max,
                root_blocks,
                block_cells[axis],
                level,
            )
            block_location = (block_i, block_j, block_k)[axis]
            if block_location != target_block:
                stream.seek(6 * header.location_size + record_data_bytes, os.SEEK_CUR)
                continue

            limits_buffer = stream.read(6 * header.location_size)
            limits = np.frombuffer(limits_buffer, dtype=location_dtype)
            if limits.size != 6:
                raise RuntimeError(f"{path}: truncated MeshBlock coordinate record")
            if plane == "x":
                extent = (limits[2], limits[3], limits[4], limits[5])
            elif plane == "y":
                extent = (limits[0], limits[1], limits[4], limits[5])
            else:
                extent = (limits[0], limits[1], limits[2], limits[3])
            extent_tuple = tuple(float(value) for value in extent)
            if not _overlaps_extent(extent_tuple, half_width):
                stream.seek(record_data_bytes, os.SEEK_CUR)
                continue

            data_start = stream.tell()
            fields = {}
            for name, variable_index in sorted(
                variable_indices.items(), key=lambda item: item[1]
            ):
                stream.seek(data_start + variable_index * variable_bytes)
                values = np.fromfile(stream, dtype=variable_dtype, count=number_cells)
                if values.size != number_cells:
                    raise RuntimeError(f"{path}: truncated cell data for {name}")
                values = values.reshape(block_cells[2], block_cells[1], block_cells[0])
                if plane == "x":
                    fields[name] = values[:, :, slice_index]
                elif plane == "y":
                    fields[name] = values[:, slice_index, :]
                else:
                    fields[name] = values[slice_index, :, :]
            stream.seek(data_start + record_data_bytes)
            selected_level_counts[level] += 1
            blocks.append(
                SliceBlock(
                    extent=extent_tuple,
                    level=level,
                    logical_location=(block_i, block_j, block_k),
                    slice_shape=(
                        (block_cells[2], block_cells[1])
                        if plane == "x"
                        else (block_cells[2], block_cells[0])
                        if plane == "y"
                        else (block_cells[1], block_cells[0])
                    ),
                    fields=fields,
                )
            )

    if not blocks:
        raise RuntimeError(f"{path}: slice selected no MeshBlocks")
    return SliceData(
        header=header,
        plane=plane,
        location=location,
        blocks=blocks,
        level_counts=level_counts,
        selected_level_counts=selected_level_counts,
    )


def panel_values(
    slice_data: SliceData, panel: Panel, density_threshold: float | None
) -> list[np.ndarray]:
    values = []
    for block in slice_data.blocks:
        if panel.name == "level":
            value = np.full(block.slice_shape, block.level, dtype=float)
        else:
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                value = panel.calculate(block.fields, block.level)
            if panel.proxy and density_threshold is not None:
                value = np.where(block.fields["dens"] >= density_threshold, value, np.nan)
        values.append(np.asarray(value))
    return values


def calculate_limits(values: list[np.ndarray], scale: str) -> tuple[float, float]:
    finite_parts = [value[np.isfinite(value)] for value in values]
    finite_parts = [value for value in finite_parts if value.size]
    if not finite_parts:
        raise RuntimeError("panel contains no finite values")
    finite = np.concatenate(finite_parts)
    if scale == "level":
        return float(np.min(finite)), float(np.max(finite))
    if scale == "log":
        finite = finite[finite > 0.0]
        if not finite.size:
            raise RuntimeError("logarithmic panel contains no positive finite values")
        low, high = np.nanpercentile(finite, (1.0, 99.5))
    elif scale == "symlog":
        bound = float(np.nanpercentile(np.abs(finite), 99.5))
        if bound == 0.0:
            bound = float(np.finfo(finite.dtype).tiny)
        low, high = -bound, bound
    else:
        low, high = np.nanpercentile(finite, (0.5, 99.5))
        if low < 0.0 < high:
            bound = max(abs(low), abs(high))
            low, high = -bound, bound
    if not high > low:
        delta = max(abs(float(low)) * 1.0e-6, 1.0e-30)
        low, high = float(low) - delta, float(high) + delta
    return float(low), float(high)


def load_trajectory_position(
    trajectory_path: Path, simulation_time: float, time_offset: float
) -> dict[str, object]:
    """Interpolate the two coordinate centers from a dynbbh 21-column table."""

    columns = np.loadtxt(trajectory_path, comments="#", usecols=range(7), ndmin=2)
    if columns.shape[0] < 2:
        raise RuntimeError(f"{trajectory_path}: trajectory needs at least two rows")
    table_time = simulation_time + time_offset
    if table_time < columns[0, 0] or table_time > columns[-1, 0]:
        raise RuntimeError(
            f"trajectory time {table_time} is outside "
            f"[{columns[0, 0]}, {columns[-1, 0]}]"
        )
    centers = []
    for start in (1, 4):
        centers.append(
            [
                float(np.interp(table_time, columns[:, 0], columns[:, component]))
                for component in range(start, start + 3)
            ]
        )
    return {"table_time": table_time, "centers": centers}


def _project_center(center: list[float], plane: str) -> tuple[float, float, float]:
    if plane == "x":
        return center[1], center[2], center[0]
    if plane == "y":
        return center[0], center[2], center[1]
    return center[0], center[1], center[2]


def create_dashboard(
    slice_data: SliceData,
    panel_names: list[str],
    output_path: Path,
    half_width: float | None,
    dpi: int,
    show_grid: bool,
    trajectory: dict[str, object] | None,
    rho_mask_fraction: float,
) -> dict[str, object]:
    """Render a multi-panel AMR dashboard and return its numerical summary."""

    cache = Path(tempfile.gettempdir()) / "athenak-matplotlib-cache"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib

    matplotlib.use("agg")
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.patches import Rectangle

    number_panels = len(panel_names)
    columns = min(3, number_panels)
    rows = math.ceil(number_panels / columns)
    figure, axes = plt.subplots(
        rows, columns, figsize=(5.4 * columns, 4.8 * rows), squeeze=False
    )
    axes_flat = axes.ravel()
    if any(PANELS[name].proxy for name in panel_names):
        density_max = max(
            float(np.nanmax(block.fields["dens"])) for block in slice_data.blocks
        )
        density_threshold = (
            density_max * rho_mask_fraction if rho_mask_fraction > 0.0 else None
        )
    else:
        density_threshold = None
    summary: dict[str, object] = {
        "panels": {},
        "proxy_density_mask_fraction": rho_mask_fraction,
        "proxy_density_threshold": density_threshold,
    }
    sorted_blocks = sorted(slice_data.blocks, key=lambda block: block.level)
    xlabel, ylabel = {
        "x": ("y / M", "z / M"),
        "y": ("x / M", "z / M"),
        "z": ("x / M", "y / M"),
    }[slice_data.plane]

    for axis_plot, panel_name in zip(axes_flat, panel_names):
        panel = PANELS[panel_name]
        values = panel_values(slice_data, panel, density_threshold)
        low, high = calculate_limits(values, panel.scale)
        symlog_linthresh = None
        if panel.scale == "log":
            norm = colors.LogNorm(vmin=low, vmax=high, clip=True)
        elif panel.scale == "symlog":
            symlog_linthresh = max(high * 1.0e-6, np.finfo(float).tiny)
            norm = colors.SymLogNorm(
                linthresh=symlog_linthresh,
                vmin=low,
                vmax=high,
                clip=True,
            )
        elif panel.scale == "level":
            boundaries = np.arange(math.floor(low) - 0.5, math.ceil(high) + 1.5)
            norm = colors.BoundaryNorm(boundaries, matplotlib.colormaps[panel.cmap].N)
        else:
            norm = colors.Normalize(vmin=low, vmax=high, clip=True)

        values_by_id = {
            id(block): value for block, value in zip(slice_data.blocks, values)
        }
        for block in sorted_blocks:
            x_min, x_max, y_min, y_max = block.extent
            value = values_by_id[id(block)]
            x_edges = np.linspace(x_min, x_max, value.shape[1] + 1)
            y_edges = np.linspace(y_min, y_max, value.shape[0] + 1)
            axis_plot.pcolormesh(
                x_edges,
                y_edges,
                value,
                cmap=panel.cmap,
                norm=norm,
                shading="flat",
                rasterized=True,
            )
            if show_grid:
                axis_plot.add_patch(
                    Rectangle(
                        (x_min, y_min),
                        x_max - x_min,
                        y_max - y_min,
                        fill=False,
                        edgecolor="white",
                        linewidth=0.18,
                        alpha=0.22,
                    )
                )

        if trajectory is not None:
            marker_styles = (("o", "white"), ("s", "cyan"))
            for index, center in enumerate(trajectory["centers"]):
                x_bh, y_bh, normal_bh = _project_center(center, slice_data.plane)
                distance = abs(normal_bh - slice_data.location)
                marker, color = marker_styles[index]
                axis_plot.scatter(
                    [x_bh],
                    [y_bh],
                    marker=marker,
                    s=42,
                    facecolors="none",
                    edgecolors=color,
                    linewidths=1.1,
                    alpha=1.0 if distance < 1.0 else 0.45,
                    zorder=20,
                )

        if half_width is not None:
            axis_plot.set_xlim(-half_width, half_width)
            axis_plot.set_ylim(-half_width, half_width)
        else:
            all_extents = np.asarray([block.extent for block in slice_data.blocks])
            axis_plot.set_xlim(np.min(all_extents[:, 0]), np.max(all_extents[:, 1]))
            axis_plot.set_ylim(np.min(all_extents[:, 2]), np.max(all_extents[:, 3]))
        axis_plot.set_aspect("equal")
        axis_plot.set_xlabel(xlabel)
        axis_plot.set_ylabel(ylabel)
        axis_plot.set_title(panel.label, fontsize=10)
        colorbar = figure.colorbar(
            ScalarMappable(norm=norm, cmap=panel.cmap), ax=axis_plot, pad=0.02, shrink=0.9
        )
        if panel.scale == "level":
            colorbar.set_ticks(np.arange(math.floor(low), math.ceil(high) + 1))
        elif panel.scale == "symlog":
            inner_tick = math.sqrt(high * symlog_linthresh)
            ticks = (-high, -inner_tick, 0.0, inner_tick, high)
            colorbar.set_ticks(ticks)
            colorbar.set_ticklabels([f"{tick:.1e}" for tick in ticks])
        summary["panels"][panel_name] = {
            "display_min": low,
            "display_max": high,
            "scale": panel.scale,
            "proxy": panel.proxy,
        }

    for unused_axis in axes_flat[number_panels:]:
        unused_axis.set_visible(False)
    title = (
        f"AthenaK BBH GRMHD | t={slice_data.header.time:.6g} M, "
        f"cycle={slice_data.header.cycle} | {slice_data.plane}="
        f"{slice_data.location:g} M"
    )
    figure.suptitle(title, fontsize=13, y=0.985)
    if any(PANELS[name].proxy for name in panel_names):
        figure.text(
            0.5,
            0.012,
            "Panels marked proxy use coordinate components; covariant ADM contractions "
            "require metric output. Low-density atmosphere is masked.",
            ha="center",
            va="bottom",
            fontsize=8,
            color="darkred",
        )
    figure.tight_layout(rect=(0.0, 0.04, 1.0, 0.94), h_pad=3.2, w_pad=2.0)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=dpi)
    plt.close(figure)
    return summary


def create_history_plot(
    history_path: Path, output_path: Path, dpi: int
) -> dict[str, object]:
    """Plot conservation and energy history using AthenaK's bundled hst reader."""

    cache = Path(tempfile.gettempdir()) / "athenak-matplotlib-cache"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib

    matplotlib.use("agg")
    import matplotlib.pyplot as plt

    try:
        from athena_read import hst
    except ImportError:
        from .athena_read import hst

    data = hst(str(history_path))
    time = data["time"]
    figure, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    axes = axes.ravel()
    if "mass" in data:
        mass_reference = data["mass"][0]
        axes[0].plot(time, data["mass"] / mass_reference - 1.0)
        axes[0].set_ylabel(r"$M/M_0-1$")
    momentum_names = [name for name in ("1-mom", "2-mom", "3-mom") if name in data]
    for name in momentum_names:
        axes[1].plot(time, data[name], label=name)
    if momentum_names:
        axes[1].legend(fontsize=8)
    axes[1].set_ylabel("domain momentum")
    kinetic_names = [name for name in ("1-KE", "2-KE", "3-KE") if name in data]
    magnetic_names = [name for name in ("1-ME", "2-ME", "3-ME") if name in data]
    if kinetic_names:
        kinetic = sum(data[name] for name in kinetic_names)
        axes[2].plot(time, kinetic, label="kinetic")
    if magnetic_names:
        magnetic = sum(data[name] for name in magnetic_names)
        axes[2].plot(time, magnetic, label="magnetic")
    if kinetic_names or magnetic_names:
        axes[2].legend(fontsize=8)
        axes[2].set_yscale("log")
    axes[2].set_ylabel("integrated energy")
    axes[3].plot(time, data["dt"])
    axes[3].set_ylabel(r"$\Delta t_\mathrm{root}$")
    for axis in axes:
        axis.set_xlabel(r"$t/M$")
        axis.grid(alpha=0.25)
    figure.suptitle(f"AthenaK history: {history_path.name}")
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=dpi)
    plt.close(figure)
    diagnostics: dict[str, float] = {}
    if "mass" in data:
        diagnostics["mass_relative_change"] = float(
            data["mass"][-1] / data["mass"][0] - 1.0
        )
    if "tot-E" in data:
        diagnostics["total_energy_relative_change"] = float(
            data["tot-E"][-1] / data["tot-E"][0] - 1.0
        )
    if kinetic_names:
        kinetic = sum(data[name] for name in kinetic_names)
        diagnostics["kinetic_energy_initial"] = float(kinetic[0])
        diagnostics["kinetic_energy_final"] = float(kinetic[-1])
    if magnetic_names:
        magnetic = sum(data[name] for name in magnetic_names)
        diagnostics["magnetic_energy_initial"] = float(magnetic[0])
        diagnostics["magnetic_energy_final"] = float(magnetic[-1])
        diagnostics["magnetic_energy_growth_factor"] = float(
            magnetic[-1] / magnetic[0]
        )
    return {
        "input": str(history_path.resolve()),
        "output": str(output_path.resolve()),
        "samples": int(time.size),
        "time_min": float(time.min()),
        "time_max": float(time.max()),
        "columns": sorted(data),
        "diagnostics": diagnostics,
        "interpretation_notes": [
            "A prescribed time-dependent BBH metric can exchange coordinate energy "
            "and momentum with the fluid, so tot-E and domain momentum are not "
            "strict conservation invariants.",
            "Domain mass can change through excision, atmosphere recovery, accretion, "
            "and boundary flux. Publication analysis requires moving-surface fluxes."
        ],
    }


def expand_inputs(patterns: list[str]) -> list[Path]:
    paths = []
    for pattern in patterns:
        path = Path(pattern)
        if path.is_dir():
            paths.extend(sorted(path.glob("*.bin")))
            continue
        matches = [Path(match) for match in glob.glob(pattern)]
        if matches:
            paths.extend(match for match in matches if match.is_file())
        elif path.exists() and path.is_file():
            paths.append(path)
        else:
            raise FileNotFoundError(pattern)
    unique = sorted({path.resolve() for path in paths})
    if not unique:
        raise RuntimeError("no AthenaK binary inputs selected")
    return unique


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "inputs", nargs="+", help="AthenaK .bin files, directories, or globs"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True, help="directory for figures and JSON"
    )
    parser.add_argument("--plane", choices=("x", "y", "z"), default="z")
    parser.add_argument("--location", type=float, default=0.0, help="slice coordinate")
    parser.add_argument(
        "--extent",
        type=float,
        default=80.0,
        help="plot half-width in M; use 0 for the full selected plane",
    )
    parser.add_argument(
        "--panels",
        default="dens,press,temperature,bmag,beta_inv_proxy,level",
        help=f"comma-separated panels; choices: {','.join(PANELS)}",
    )
    parser.add_argument("--grid", action="store_true", help="outline AMR MeshBlocks")
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument(
        "--rho-mask-fraction",
        type=float,
        default=1.0e-8,
        help="mask proxy panels below this fraction of slice maximum density; 0 disables",
    )
    parser.add_argument(
        "--trajectory", type=Path, help="dynbbh 21-column trajectory for BH markers"
    )
    parser.add_argument(
        "--trajectory-time-offset",
        type=float,
        help="override <problem>/trajectory_time_offset from the binary header",
    )
    parser.add_argument("--history", type=Path, help="optional AthenaK .hst file")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    panel_names = [name.strip() for name in args.panels.split(",") if name.strip()]
    invalid_panels = sorted(set(panel_names).difference(PANELS))
    if invalid_panels:
        raise SystemExit(
            f"invalid panels: {', '.join(invalid_panels)}; choices: {', '.join(PANELS)}"
        )
    if not panel_names:
        raise SystemExit("at least one panel is required")
    half_width = None if args.extent == 0.0 else args.extent
    if half_width is not None and half_width <= 0.0:
        raise SystemExit("--extent must be positive, or 0 for the full plane")
    if args.rho_mask_fraction < 0.0:
        raise SystemExit("--rho-mask-fraction must be non-negative")
    input_paths = expand_inputs(args.inputs)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    run_manifest = {
        "schema": 1,
        "classification": "athenak-bbh-grmhd-amr-postprocessing",
        "plane": args.plane,
        "location_M": args.location,
        "extent_M": half_width,
        "panels": panel_names,
        "proxy_density_mask_fraction": args.rho_mask_fraction,
        "proxy_warning": PROXY_WARNING
        if any(PANELS[name].proxy for name in panel_names)
        else None,
        "frames": [],
    }

    for input_path in input_paths:
        slice_data = read_slice(
            input_path, panel_names, args.plane, args.location, half_width
        )
        trajectory = None
        if args.trajectory is not None:
            if args.trajectory_time_offset is not None:
                time_offset = args.trajectory_time_offset
            else:
                time_offset = float(
                    slice_data.header.parameters.get("problem", {}).get(
                        "trajectory_time_offset", "0"
                    )
                )
            trajectory = load_trajectory_position(
                args.trajectory, slice_data.header.time, time_offset
            )
        stem = (
            input_path.name[:-4]
            if input_path.name.endswith(".bin")
            else input_path.name
        )
        figure_path = args.output_dir / (
            f"{stem}.dashboard_{args.plane}{args.location:g}.png"
        )
        numerical_summary = create_dashboard(
            slice_data,
            panel_names,
            figure_path,
            half_width,
            args.dpi,
            args.grid,
            trajectory,
            args.rho_mask_fraction,
        )
        frame_summary = {
            "input": str(input_path),
            "input_bytes": input_path.stat().st_size,
            "output": str(figure_path.resolve()),
            "time_M": slice_data.header.time,
            "cycle": slice_data.header.cycle,
            "variables": list(slice_data.header.variables),
            "meshblocks": int(sum(slice_data.level_counts.values())),
            "meshblocks_by_level": {
                str(level): count
                for level, count in sorted(slice_data.level_counts.items())
            },
            "slice_meshblocks": len(slice_data.blocks),
            "slice_meshblocks_by_level": {
                str(level): count
                for level, count in sorted(slice_data.selected_level_counts.items())
            },
            "trajectory": trajectory,
            **numerical_summary,
        }
        frame_manifest_path = figure_path.with_suffix(".json")
        frame_manifest_path.write_text(
            json.dumps(frame_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        run_manifest["frames"].append(frame_summary)
        print(figure_path)

    if args.history is not None:
        history_output = args.output_dir / f"{args.history.stem}.history.png"
        run_manifest["history"] = create_history_plot(
            args.history, history_output, args.dpi
        )
        print(history_output)

    manifest_path = args.output_dir / "postprocess-manifest.json"
    manifest_path.write_text(
        json.dumps(run_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(manifest_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.stderr.close()
        raise SystemExit(1)
