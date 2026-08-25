#!/usr/bin/env python3
"""Create AMR-aware BBH GRMHD dashboards from AthenaK binary output.

The stock ``plot_slice.py`` remains the general single-panel plotting tool.  This
script follows its binary-slice selection logic, but reads every requested field in
one pass so that multi-panel figures remain practical for gigabyte-scale AMR dumps.

The ordinary ``mhd_w_bcc`` output does not contain the dynamical ADM metric.  Density,
pressure, temperature, densitized magnetic components, and primitive velocity
components are therefore plotted exactly as stored.  Quantities involving contractions
of those components are explicitly labelled ``proxy`` and must not be interpreted as
covariant GR diagnostics.  The separate ``mhd_gr_diagnostics`` output contains native
metric-aware magnetic invariants and the Lorentz factor.  New files also contain
``gr_excision_mask``; physical GR panels automatically hide excised cells.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass, field
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
    "Stored-component proxy: DynGRMHD bcc is sqrt(det(gamma_ij))*B^i, and the "
    "mhd_w_bcc dump does not contain the metric needed for a covariant contraction."
)

EXCISION_MASK_FIELD = "gr_excision_mask"
MASKED_GR_PANELS = {"gr_bsq", "gr_lorentz", "gr_sigma", "gr_beta_inv"}


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
    presliced: bool
    sampling_metadata: dict[str, object] = field(default_factory=dict)


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
        r"$\sqrt{(\tilde B^x)^2+(\tilde B^y)^2+(\tilde B^z)^2}$ (stored-bcc proxy)",
        ("dens", "bcc1", "bcc2", "bcc3"),
        lambda fields, level: np.sqrt(_b2(fields)),
        "viridis",
        "log",
        True,
    ),
    "beta_inv_proxy": Panel(
        "beta_inv_proxy",
        r"$[(\tilde B^x)^2+(\tilde B^y)^2+(\tilde B^z)^2]/(2p)$ (proxy)",
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
        r"$[(\tilde B^x)^2+(\tilde B^y)^2+(\tilde B^z)^2]/\rho$ (proxy)",
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
    "gr_bsq": Panel(
        "gr_bsq",
        r"$b^2=b^\mu b_\mu$",
        ("gr_bsq",),
        _identity("gr_bsq"),
        "viridis",
        "log",
    ),
    "gr_lorentz": Panel(
        "gr_lorentz",
        r"$W=\alpha u^0$",
        ("gr_lorentz",),
        _identity("gr_lorentz"),
        "magma",
        "linear",
    ),
    "gr_sigma": Panel(
        "gr_sigma",
        r"$\sigma=b^2/\rho$",
        ("gr_sigma",),
        _identity("gr_sigma"),
        "cubehelix",
        "log",
    ),
    "gr_beta_inv": Panel(
        "gr_beta_inv",
        r"$\beta_\mathrm{mag}^{-1}=b^2/(2p)$",
        ("gr_beta_inv",),
        _identity("gr_beta_inv"),
        "cividis",
        "log",
    ),
    "excision_mask": Panel(
        "excision_mask",
        "excision mask (1 = excluded)",
        (EXCISION_MASK_FIELD,),
        _identity(EXCISION_MASK_FIELD),
        "gray_r",
        "linear",
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
    extent: tuple[float, float, float, float],
    half_width: float | None,
    center: tuple[float, float] = (0.0, 0.0),
) -> bool:
    if half_width is None:
        return True
    return not (
        extent[1] <= center[0] - half_width
        or extent[0] >= center[0] + half_width
        or extent[3] <= center[1] - half_width
        or extent[2] >= center[1] + half_width
    )


def _bracketing_cell_centers(
    location: float,
    domain_min: float,
    domain_max: float,
    root_cells: int,
    cells_per_block: int,
    level: int,
) -> tuple[tuple[int, int], tuple[int, int], float]:
    """Return block/cell indices and the upper-center interpolation weight."""

    number_cells = root_cells * 2**level
    spacing = (domain_max - domain_min) / number_cells
    fractional_center_index = (location - domain_min) / spacing - 0.5
    lower_index = math.floor(fractional_center_index)
    upper_index = lower_index + 1
    if lower_index < 0:
        lower_index = upper_index = 0
        upper_weight = 0.0
    elif upper_index >= number_cells:
        lower_index = upper_index = number_cells - 1
        upper_weight = 0.0
    else:
        lower_center = domain_min + (lower_index + 0.5) * spacing
        upper_weight = (location - lower_center) / spacing
        upper_weight = min(max(upper_weight, 0.0), 1.0)
    lower = lower_index // cells_per_block, lower_index % cells_per_block
    upper = upper_index // cells_per_block, upper_index % cells_per_block
    return lower, upper, upper_weight


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
    # New diagnostics carry an exact cell mask from the evolution.  Continue to read
    # legacy four-field files, but automatically load and apply the mask when present.
    if (
        EXCISION_MASK_FIELD in header.variables
        and any(panel_name in MASKED_GR_PANELS for panel_name in panel_names)
    ):
        dependencies.add(EXCISION_MASK_FIELD)
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
    file_presliced: bool | None = None

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
            presliced = (
                block_cells[axis] == 1 and block_cells_config[axis] > 1
            )
            if first_record:
                shape_is_supported = all(
                    actual == configured or
                    (dimension == axis and actual == 1 and configured > 1)
                    for dimension, (actual, configured) in enumerate(
                        zip(block_cells, block_cells_config)
                    )
                )
                if not shape_is_supported:
                    raise RuntimeError(
                        f"{path}: output block shape {block_cells} differs from "
                        f"configured shape {block_cells_config}; the file may be sliced "
                        f"along a different axis"
                    )
                file_presliced = presliced
                first_record = False
            elif presliced != file_presliced:
                raise RuntimeError(f"{path}: inconsistent MeshBlock slice shapes")
            number_cells = math.prod(block_cells)
            variable_bytes = number_cells * header.variable_size
            record_data_bytes = len(header.variables) * variable_bytes
            if presliced:
                slice_index = 0
            else:
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
                    stream.seek(
                        6 * header.location_size + record_data_bytes, os.SEEK_CUR
                    )
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
        presliced=bool(file_presliced),
    )


def read_linearly_interpolated_slice(
    path: Path,
    panel_names: list[str],
    plane: str,
    location: float,
    half_width: float | None,
    in_plane_center: tuple[float, float] = (0.0, 0.0),
    raster_resolution: int | None = None,
) -> SliceData:
    """Interpolate an arbitrary AMR plane between bracketing cell centers.

    Both normal-direction samples must have identical projected leaf topology.  The
    routine aborts instead of silently mixing refinement levels when a moving plane
    crosses an asymmetric AMR boundary.
    """

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
    if (
        EXCISION_MASK_FIELD in header.variables
        and any(panel_name in MASKED_GR_PANELS for panel_name in panel_names)
    ):
        dependencies.add(EXCISION_MASK_FIELD)
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
    if not domain_min <= location <= domain_max:
        raise RuntimeError(
            f"{path}: slice location {location} lies outside "
            f"[{domain_min}, {domain_max}]"
        )
    if raster_resolution is not None:
        if raster_resolution < 32:
            raise ValueError("raster_resolution must be at least 32")
        if half_width is None:
            raise ValueError("rasterized interpolation requires a finite extent")
        if "level" in panel_names:
            raise ValueError("rasterized interpolation does not support the level panel")

    location_dtype = np.dtype("<f4" if header.location_size == 4 else "<f8")
    variable_dtype = np.dtype("<f4" if header.variable_size == 4 else "<f8")
    level_counts: Counter = Counter()
    brackets: dict[int, tuple[tuple[int, int], tuple[int, int], float]] = {}
    side_blocks: dict[str, dict[tuple[int, int, int], SliceBlock]] = {
        "lower": {},
        "upper": {},
    }
    normal_coordinates: dict[str, dict[tuple[int, int, int], float]] = {
        "lower": {},
        "upper": {},
    }

    def projected_key(
        level: int, logical_location: tuple[int, int, int]
    ) -> tuple[int, int, int]:
        coordinates = list(logical_location)
        del coordinates[axis]
        return level, coordinates[0], coordinates[1]

    with path.open("rb") as stream:
        stream.seek(header.data_offset)
        file_size = path.stat().st_size
        while stream.tell() < file_size:
            index_bytes = stream.read(24)
            if not index_bytes:
                break
            if len(index_bytes) != 24:
                raise RuntimeError(f"{path}: truncated MeshBlock index record")
            indices = np.asarray(struct.unpack("=6i", index_bytes)) - nghost
            logical = struct.unpack("=4i", stream.read(16))
            block_i, block_j, block_k, level = logical
            logical_location = (block_i, block_j, block_k)
            level_counts[level] += 1
            block_cells = (
                int(indices[1] - indices[0] + 1),
                int(indices[3] - indices[2] + 1),
                int(indices[5] - indices[4] + 1),
            )
            if block_cells[axis] == 1 and block_cells_config[axis] > 1:
                raise RuntimeError(
                    f"{path}: the output already stores one cell along {plane}; "
                    "linear spatial interpolation requires a full 3-D dump"
                )
            if block_cells != block_cells_config:
                raise RuntimeError(
                    f"{path}: output block shape {block_cells} differs from "
                    f"configured shape {block_cells_config}"
                )
            number_cells = math.prod(block_cells)
            variable_bytes = number_cells * header.variable_size
            record_data_bytes = len(header.variables) * variable_bytes
            bracket = brackets.setdefault(
                level,
                _bracketing_cell_centers(
                    location,
                    domain_min,
                    domain_max,
                    mesh_cells[axis],
                    block_cells[axis],
                    level,
                ),
            )
            lower_target, upper_target, _ = bracket
            block_location = logical_location[axis]
            targets: dict[str, int] = {}
            if block_location == lower_target[0]:
                targets["lower"] = lower_target[1]
            if block_location == upper_target[0]:
                targets["upper"] = upper_target[1]
            if not targets:
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
            if not _overlaps_extent(extent_tuple, half_width, in_plane_center):
                stream.seek(record_data_bytes, os.SEEK_CUR)
                continue

            data_start = stream.tell()
            fields_by_side: dict[str, dict[str, np.ndarray]] = {
                side: {} for side in targets
            }
            for name, variable_index in sorted(
                variable_indices.items(), key=lambda item: item[1]
            ):
                stream.seek(data_start + variable_index * variable_bytes)
                values = np.fromfile(stream, dtype=variable_dtype, count=number_cells)
                if values.size != number_cells:
                    raise RuntimeError(f"{path}: truncated cell data for {name}")
                values = values.reshape(block_cells[2], block_cells[1], block_cells[0])
                for side, slice_index in targets.items():
                    if plane == "x":
                        fields_by_side[side][name] = values[:, :, slice_index]
                    elif plane == "y":
                        fields_by_side[side][name] = values[:, slice_index, :]
                    else:
                        fields_by_side[side][name] = values[slice_index, :, :]
            stream.seek(data_start + record_data_bytes)
            slice_shape = (
                (block_cells[2], block_cells[1])
                if plane == "x"
                else (block_cells[2], block_cells[0])
                if plane == "y"
                else (block_cells[1], block_cells[0])
            )
            key = projected_key(level, logical_location)
            for side, fields in fields_by_side.items():
                if key in side_blocks[side]:
                    raise RuntimeError(
                        f"{path}: duplicate projected {side} block for {key}"
                    )
                side_blocks[side][key] = SliceBlock(
                    extent=extent_tuple,
                    level=level,
                    logical_location=logical_location,
                    slice_shape=slice_shape,
                    fields=fields,
                )
                target = lower_target if side == "lower" else upper_target
                global_cell = target[0] * block_cells[axis] + target[1]
                spacing = (domain_max - domain_min) / (mesh_cells[axis] * 2**level)
                normal_coordinates[side][key] = (
                    domain_min + (global_cell + 0.5) * spacing
                )

    lower_blocks = side_blocks["lower"]
    upper_blocks = side_blocks["upper"]
    if raster_resolution is not None:
        return _rasterize_interpolated_slice(
            header,
            plane,
            location,
            float(half_width),
            in_plane_center,
            raster_resolution,
            lower_blocks,
            upper_blocks,
            normal_coordinates,
            level_counts,
        )
    if lower_blocks.keys() != upper_blocks.keys():
        missing_lower = sorted(upper_blocks.keys() - lower_blocks.keys())
        missing_upper = sorted(lower_blocks.keys() - upper_blocks.keys())
        raise RuntimeError(
            f"{path}: AMR topology differs across moving {plane}={location:g}; "
            f"missing below={missing_lower[:4]}, missing above={missing_upper[:4]}"
        )
    if not lower_blocks:
        raise RuntimeError(f"{path}: interpolated slice selected no MeshBlocks")

    blocks: list[SliceBlock] = []
    selected_level_counts: Counter = Counter()
    for key in sorted(lower_blocks):
        low_block = lower_blocks[key]
        high_block = upper_blocks[key]
        if (
            low_block.slice_shape != high_block.slice_shape
            or not np.allclose(low_block.extent, high_block.extent, rtol=0.0, atol=0.0)
            or low_block.fields.keys() != high_block.fields.keys()
        ):
            raise RuntimeError(f"{path}: incompatible moving-plane blocks: {key}")
        upper_weight = brackets[key[0]][2]
        lower_weight = 1.0 - upper_weight
        fields = {
            name: lower_weight * low_block.fields[name]
            + upper_weight * high_block.fields[name]
            for name in low_block.fields
        }
        blocks.append(
            SliceBlock(
                extent=low_block.extent,
                level=low_block.level,
                logical_location=low_block.logical_location,
                slice_shape=low_block.slice_shape,
                fields=fields,
            )
        )
        selected_level_counts[low_block.level] += 1
    return SliceData(
        header=header,
        plane=plane,
        location=location,
        blocks=blocks,
        level_counts=level_counts,
        selected_level_counts=selected_level_counts,
        presliced=False,
    )


def _rasterize_interpolated_slice(
    header: BinaryHeader,
    plane: str,
    location: float,
    half_width: float,
    in_plane_center: tuple[float, float],
    resolution: int,
    lower_blocks: dict[tuple[int, int, int], SliceBlock],
    upper_blocks: dict[tuple[int, int, int], SliceBlock],
    normal_coordinates: dict[str, dict[tuple[int, int, int], float]],
    level_counts: Counter,
) -> SliceData:
    """Resample different AMR topologies before normal-direction interpolation."""

    from scipy.ndimage import map_coordinates

    if not lower_blocks or not upper_blocks:
        raise RuntimeError("rasterized interpolation has an empty bracketing side")
    field_names = set(next(iter(lower_blocks.values())).fields)
    if any(set(block.fields) != field_names for block in lower_blocks.values()):
        raise RuntimeError("lower-side raster blocks have inconsistent fields")
    if any(set(block.fields) != field_names for block in upper_blocks.values()):
        raise RuntimeError("upper-side raster blocks have inconsistent fields")

    horizontal_min = in_plane_center[0] - half_width
    horizontal_max = in_plane_center[0] + half_width
    vertical_min = in_plane_center[1] - half_width
    vertical_max = in_plane_center[1] + half_width
    horizontal = np.linspace(
        horizontal_min,
        horizontal_max,
        resolution,
        endpoint=False,
        dtype=float,
    )
    vertical = np.linspace(
        vertical_min,
        vertical_max,
        resolution,
        endpoint=False,
        dtype=float,
    )
    pixel_size = 2.0 * half_width / resolution
    horizontal += 0.5 * pixel_size
    vertical += 0.5 * pixel_size

    def rasterize_side(
        blocks: dict[tuple[int, int, int], SliceBlock],
        coordinates: dict[tuple[int, int, int], float],
    ) -> tuple[dict[str, np.ndarray], np.ndarray]:
        raster_fields = {
            name: np.full((resolution, resolution), np.nan, dtype=float)
            for name in field_names
        }
        normal_raster = np.full((resolution, resolution), np.nan, dtype=float)
        for key, block in sorted(blocks.items(), key=lambda item: item[0][0]):
            x_min, x_max, y_min, y_max = block.extent
            x_indices = np.flatnonzero(
                (horizontal >= x_min) & (horizontal < x_max)
            )
            y_indices = np.flatnonzero((vertical >= y_min) & (vertical < y_max))
            if not x_indices.size or not y_indices.size:
                continue
            sample = next(iter(block.fields.values()))
            cell_width = (x_max - x_min) / sample.shape[1]
            cell_height = (y_max - y_min) / sample.shape[0]
            x_coordinates = (horizontal[x_indices] - x_min) / cell_width - 0.5
            y_coordinates = (vertical[y_indices] - y_min) / cell_height - 0.5
            yy, xx = np.meshgrid(y_coordinates, x_coordinates, indexing="ij")
            target = np.ix_(y_indices, x_indices)
            for name, values in block.fields.items():
                order = 0 if name == EXCISION_MASK_FIELD else 1
                raster_fields[name][target] = map_coordinates(
                    values,
                    (yy, xx),
                    order=order,
                    mode="nearest",
                    prefilter=False,
                )
            normal_raster[target] = coordinates[key]
        return raster_fields, normal_raster

    lower_fields, lower_coordinate = rasterize_side(
        lower_blocks, normal_coordinates["lower"]
    )
    upper_fields, upper_coordinate = rasterize_side(
        upper_blocks, normal_coordinates["upper"]
    )
    lower_valid = np.isfinite(lower_coordinate)
    upper_valid = np.isfinite(upper_coordinate)
    uncovered = ~(lower_valid | upper_valid)
    if np.any(uncovered):
        missing = int(np.count_nonzero(uncovered))
        raise RuntimeError(
            f"rasterized moving plane has {missing} cells missing on both sides"
        )
    paired = lower_valid & upper_valid
    separation = upper_coordinate[paired] - lower_coordinate[paired]
    if np.any(separation <= 0.0):
        raise RuntimeError("rasterized bracketing cell centers are not ordered")
    upper_weight = np.zeros((resolution, resolution), dtype=float)
    upper_weight[paired] = np.clip(
        (location - lower_coordinate[paired]) / separation, 0.0, 1.0
    )
    fields: dict[str, np.ndarray] = {}
    for name in field_names:
        value = np.where(lower_valid, lower_fields[name], upper_fields[name])
        if name == EXCISION_MASK_FIELD:
            value[paired] = np.maximum(
                lower_fields[name][paired], upper_fields[name][paired]
            )
        else:
            value[paired] = (
                (1.0 - upper_weight[paired]) * lower_fields[name][paired]
                + upper_weight[paired] * upper_fields[name][paired]
            )
        fields[name] = value
    maximum_level = max(
        max(key[0] for key in lower_blocks), max(key[0] for key in upper_blocks)
    )
    block = SliceBlock(
        extent=(horizontal_min, horizontal_max, vertical_min, vertical_max),
        level=maximum_level,
        logical_location=(0, 0, 0),
        slice_shape=(resolution, resolution),
        fields=fields,
    )
    return SliceData(
        header=header,
        plane=plane,
        location=location,
        blocks=[block],
        level_counts=level_counts,
        selected_level_counts=Counter({maximum_level: 1}),
        presliced=False,
        sampling_metadata={
            "raster_resolution": resolution,
            "two_sided_fraction": float(np.count_nonzero(paired) / paired.size),
            "one_sided_amr_fallback_fraction": float(
                np.count_nonzero(lower_valid ^ upper_valid) / paired.size
            ),
            "one_sided_amr_fallback_policy": (
                "use the available leaf-side sample without normal interpolation"
            ),
        },
    )


def read_cell_face_interpolated_slice(
    path: Path,
    panel_names: list[str],
    plane: str,
    location: float,
    half_width: float | None,
) -> SliceData:
    """Average the two cell-center planes bracketing an exact cell face.

    This is intended for equatorial slices such as ``z=0`` in a symmetric even-cell
    domain.  It requires a full 3-D output.  A file that AthenaK already reduced to one
    stored cell along the slice axis cannot be reconstructed after the fact.
    """

    header = read_binary_header(path)
    axis = {"x": 0, "y": 1, "z": 2}[plane]
    domain_min = _get_float(header, "mesh", f"x{axis + 1}min")
    domain_max = _get_float(header, "mesh", f"x{axis + 1}max")
    domain_width = domain_max - domain_min
    if not domain_width > 0.0:
        raise RuntimeError(f"{path}: invalid slice-axis domain")
    epsilon = domain_width * 1.0e-12
    lower = read_slice(
        path, panel_names, plane, location - epsilon, half_width
    )
    upper = read_slice(
        path, panel_names, plane, location + epsilon, half_width
    )
    if lower.presliced or upper.presliced:
        raise RuntimeError(
            f"{path}: the output already stores one cell along {plane}; "
            "two-sided spatial interpolation requires a full 3-D dump"
        )

    mesh_cells = _get_int(header, "mesh", f"nx{axis + 1}")
    for level in sorted(lower.selected_level_counts):
        spacing = domain_width / (mesh_cells * 2**level)
        face_coordinate = (location - domain_min) / spacing
        if not math.isclose(
            face_coordinate,
            round(face_coordinate),
            rel_tol=0.0,
            abs_tol=1.0e-9,
        ):
            raise RuntimeError(
                f"{path}: {plane}={location:g} is not a cell face at physical "
                f"level {level}; general off-face interpolation is not implemented"
            )

    def projected_key(block: SliceBlock) -> tuple[int, int, int]:
        coordinates = list(block.logical_location)
        del coordinates[axis]
        return block.level, coordinates[0], coordinates[1]

    lower_blocks = {projected_key(block): block for block in lower.blocks}
    upper_blocks = {projected_key(block): block for block in upper.blocks}
    if lower_blocks.keys() != upper_blocks.keys():
        missing_lower = sorted(upper_blocks.keys() - lower_blocks.keys())
        missing_upper = sorted(lower_blocks.keys() - upper_blocks.keys())
        raise RuntimeError(
            f"{path}: AMR topology differs across {plane}={location:g}; "
            f"missing below={missing_lower[:4]}, missing above={missing_upper[:4]}"
        )

    blocks: list[SliceBlock] = []
    for key in sorted(lower_blocks):
        low_block = lower_blocks[key]
        high_block = upper_blocks[key]
        if (
            low_block.slice_shape != high_block.slice_shape
            or not np.allclose(low_block.extent, high_block.extent, rtol=0.0, atol=0.0)
            or low_block.fields.keys() != high_block.fields.keys()
        ):
            raise RuntimeError(f"{path}: incompatible blocks across slice face: {key}")
        fields = {
            name: 0.5 * (low_block.fields[name] + high_block.fields[name])
            for name in low_block.fields
        }
        blocks.append(
            SliceBlock(
                extent=low_block.extent,
                level=low_block.level,
                logical_location=low_block.logical_location,
                slice_shape=low_block.slice_shape,
                fields=fields,
            )
        )

    return SliceData(
        header=header,
        plane=plane,
        location=location,
        blocks=blocks,
        level_counts=lower.level_counts,
        selected_level_counts=lower.selected_level_counts,
        presliced=False,
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
            if panel.name in MASKED_GR_PANELS and EXCISION_MASK_FIELD in block.fields:
                value = np.where(
                    block.fields[EXCISION_MASK_FIELD] < 0.5, value, np.nan
                )
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
    fixed_limits: dict[str, tuple[float, float]] | None = None,
    coordinate_origin: tuple[float, float] | None = None,
    coordinate_origin_label: str | None = None,
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
        "coordinate_origin": coordinate_origin,
        "coordinate_origin_label": coordinate_origin_label,
        "sampling_metadata": slice_data.sampling_metadata,
    }
    sorted_blocks = sorted(slice_data.blocks, key=lambda block: block.level)
    horizontal_name, vertical_name = {
        "x": ("y", "z"),
        "y": ("x", "z"),
        "z": ("x", "y"),
    }[slice_data.plane]
    if coordinate_origin is None:
        xlabel = f"{horizontal_name} / M"
        ylabel = f"{vertical_name} / M"
        horizontal_offset = vertical_offset = 0.0
    else:
        label = coordinate_origin_label or "origin"
        xlabel = f"({horizontal_name}-{horizontal_name}_{label}) / M"
        ylabel = f"({vertical_name}-{vertical_name}_{label}) / M"
        horizontal_offset, vertical_offset = coordinate_origin

    for axis_plot, panel_name in zip(axes_flat, panel_names):
        panel = PANELS[panel_name]
        values = panel_values(slice_data, panel, density_threshold)
        if fixed_limits is not None and panel_name in fixed_limits:
            low, high = fixed_limits[panel_name]
            if not (math.isfinite(low) and math.isfinite(high) and high > low):
                raise ValueError(
                    f"invalid fixed limits for {panel_name}: {low}, {high}"
                )
            if panel.scale == "log" and low <= 0.0:
                raise ValueError(
                    f"logarithmic fixed limit for {panel_name} must be positive"
                )
            if panel.scale == "symlog" and not math.isclose(
                low, -high, rel_tol=1.0e-12, abs_tol=0.0
            ):
                raise ValueError(
                    f"symmetric-log fixed limits for {panel_name} must be symmetric"
                )
            limit_policy = "fixed-across-series"
        else:
            low, high = calculate_limits(values, panel.scale)
            limit_policy = "per-frame-robust-percentile"
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
            x_edges = (
                np.linspace(x_min, x_max, value.shape[1] + 1) - horizontal_offset
            )
            y_edges = (
                np.linspace(y_min, y_max, value.shape[0] + 1) - vertical_offset
            )
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
                x_bh -= horizontal_offset
                y_bh -= vertical_offset
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

        if coordinate_origin is not None:
            if half_width is None:
                raise ValueError("a shifted coordinate origin requires a finite extent")
            axis_plot.set_xlim(-half_width, half_width)
            axis_plot.set_ylim(-half_width, half_width)
        elif half_width is not None:
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
            "display_limit_policy": limit_policy,
            "scale": panel.scale,
            "proxy": panel.proxy,
        }

    for unused_axis in axes_flat[number_panels:]:
        unused_axis.set_visible(False)
    if coordinate_origin_label is None:
        title = (
            f"AthenaK BBH GRMHD | t={slice_data.header.time:.6g} M, "
            f"cycle={slice_data.header.cycle} | {slice_data.plane}="
            f"{slice_data.location:g} M"
        )
    else:
        title = (
            f"AthenaK BBH GRMHD | t={slice_data.header.time:.6g} M, "
            f"cycle={slice_data.header.cycle}\n"
            f"{coordinate_origin_label}-following plane | {slice_data.plane}="
            f"{slice_data.location:g} M"
        )
    if any(PANELS[name].proxy for name in panel_names):
        figure.text(
            0.5,
            0.012,
            "Panels marked proxy use stored components; DynGRMHD bcc is densitized by "
            "sqrt(det gamma). Covariant contractions require mhd_gr_diagnostics. "
            "Low-density atmosphere is masked.",
            ha="center",
            va="bottom",
            fontsize=8,
            color="darkred",
        )
    # Reserve the title band explicitly after laying out math-heavy panel labels.
    # Otherwise tight_layout can move native-GR axes over the figure suptitle.
    figure.tight_layout(rect=(0.0, 0.04, 1.0, 0.90), h_pad=3.2, w_pad=2.0)
    figure.suptitle(title, fontsize=13, y=0.975)
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
    bbh_history = "baryon_m" in data
    momentum_names: list[str] = []
    kinetic_names: list[str] = []
    magnetic_names: list[str] = []
    if bbh_history:
        baryon_reference = data["baryon_m"][0]
        axes[0].plot(time, data["baryon_m"] / baryon_reference - 1.0)
        axes[0].set_ylabel(r"$M_b/M_{b,0}-1$")
        axes[1].plot(time, data["bh_sep"], label="separation")
        axes[1].set_ylabel(r"$d/M$")
        omega_axis = axes[1].twinx()
        omega_axis.plot(time, data["orb_omega"], color="tab:orange", label="omega")
        omega_axis.set_ylabel(r"$M\Omega$")
        for name, label in (
            ("pgas_prp", r"$\int p\sqrt{\gamma}\,d^3x$"),
            ("emag_prp", r"$\int b^2\sqrt{\gamma}\,d^3x/2$"),
        ):
            axes[2].plot(time, data[name], label=label)
        axes[2].set_yscale("log")
        axes[2].set_ylabel("proper-volume integrals")
        axes[2].legend(fontsize=8)
        axes[3].plot(time, data["rho_max"], label=r"$\rho_\mathrm{max}$")
        axes[3].plot(time, data["sigma_max"], label=r"$\sigma_\mathrm{max}$")
        axes[3].set_yscale("log")
        axes[3].set_ylabel("outside-excision maxima")
        axes[3].legend(fontsize=8)
    else:
        if "mass" in data:
            mass_reference = data["mass"][0]
            axes[0].plot(time, data["mass"] / mass_reference - 1.0)
            axes[0].set_ylabel(r"$M/M_0-1$")
        momentum_names = [
            name for name in ("1-mom", "2-mom", "3-mom") if name in data
        ]
        for name in momentum_names:
            axes[1].plot(time, data[name], label=name)
        if momentum_names:
            axes[1].legend(fontsize=8)
        axes[1].set_ylabel("domain momentum")
        kinetic_names = [
            name for name in ("1-KE", "2-KE", "3-KE") if name in data
        ]
        magnetic_names = [
            name for name in ("1-ME", "2-ME", "3-ME") if name in data
        ]
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
    if bbh_history:
        diagnostics.update(
            {
                "baryon_mass_relative_change": float(
                    data["baryon_m"][-1] / data["baryon_m"][0] - 1.0
                ),
                "inner_mass_fraction_initial": float(
                    data["inner_D"][0] / data["baryon_m"][0]
                ),
                "inner_mass_fraction_final": float(
                    data["inner_D"][-1] / data["baryon_m"][-1]
                ),
                "separation_final": float(data["bh_sep"][-1]),
                "orbital_omega_final": float(data["orb_omega"][-1]),
                "D_weighted_lorentz_final": float(
                    data["lor_D"][-1] / data["baryon_m"][-1]
                ),
                "D_weighted_sigma_final": float(
                    data["sigma_D"][-1] / data["baryon_m"][-1]
                ),
                "integrated_beta_inverse_final": float(
                    data["emag_prp"][-1] / data["pgas_prp"][-1]
                ),
                "rho_max_final": float(data["rho_max"][-1]),
                "sigma_max_final": float(data["sigma_max"][-1]),
            }
        )
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
        "interpretation_notes": (
            [
                "BBH user integrals omit the excision-floor mask. baryon_m and inner_D "
                "use the densitized conserved rest mass; rho_prp, pgas_prp, and "
                "emag_prp use proper spatial volume.",
                "angmom_z is a coordinate angular-momentum proxy about the instantaneous "
                "binary mass center. Publication accretion claims still require "
                "moving-surface mass and magnetic-flux diagnostics.",
            ]
            if bbh_history
            else [
                "A prescribed time-dependent BBH metric can exchange coordinate energy "
                "and momentum with the fluid, so tot-E and domain momentum are not "
                "strict conservation invariants.",
                "Domain mass can change through excision, atmosphere recovery, accretion, "
                "and boundary flux. Publication analysis requires moving-surface fluxes.",
            ]
        ),
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
    parser.add_argument(
        "--interpolate-plane",
        action="store_true",
        help="average the two full-3D cell-center planes bracketing a cell face",
    )
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
        "plane_sampling": (
            "two-sided-cell-center-linear-interpolation"
            if args.interpolate_plane
            else "stored-or-nearest-cell"
        ),
        "frames": [],
    }

    for input_path in input_paths:
        reader = (
            read_cell_face_interpolated_slice
            if args.interpolate_plane
            else read_slice
        )
        slice_data = reader(input_path, panel_names, args.plane, args.location, half_width)
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
            "excision_mask_available": EXCISION_MASK_FIELD
            in slice_data.header.variables,
            "excision_mask_applied": EXCISION_MASK_FIELD
            in slice_data.header.variables
            and any(name in MASKED_GR_PANELS for name in panel_names),
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
