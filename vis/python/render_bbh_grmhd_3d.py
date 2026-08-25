#!/usr/bin/env python3
"""Render BBH GRMHD sigma and density isosurfaces from full 3-D AMR dumps."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys

import numpy as np

from make_bbh_grmhd_movie import load_trajectory_table, trajectory_position
from plot_bbh_grmhd import read_binary_header


def _parameter_int(header, section: str, key: str) -> int:
    return int(header.parameters[section][key])


def read_uniform_volume(
    path: Path,
    fields: list[str],
    center: tuple[float, float, float],
    half_width: float,
    resolution: int,
) -> dict[str, np.ndarray]:
    """Trilinearly resample leaf MeshBlocks onto a uniform cell-centered cube."""

    from scipy.ndimage import map_coordinates

    header = read_binary_header(path)
    missing = sorted(set(fields).difference(header.variables))
    if missing:
        raise RuntimeError(f"{path}: missing fields {', '.join(missing)}")
    if header.location_size not in (4, 8) or header.variable_size not in (4, 8):
        raise RuntimeError("unsupported AthenaK binary scalar size")
    variable_indices = {name: header.variables.index(name) for name in fields}
    nghost = _parameter_int(header, "mesh", "nghost")
    location_dtype = np.dtype("<f4" if header.location_size == 4 else "<f8")
    variable_dtype = np.dtype("<f4" if header.variable_size == 4 else "<f8")
    axes = []
    spacing = 2.0 * half_width / resolution
    for origin in center:
        coordinates = np.linspace(
            origin - half_width,
            origin + half_width,
            resolution,
            endpoint=False,
            dtype=float,
        )
        axes.append(coordinates + 0.5 * spacing)
    x_axis, y_axis, z_axis = axes
    output = {
        name: np.full((resolution, resolution, resolution), np.nan, dtype=np.float32)
        for name in fields
    }

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
            logical_bytes = stream.read(16)
            if len(logical_bytes) != 16:
                raise RuntimeError(f"{path}: truncated MeshBlock logical record")
            block_cells = (
                int(indices[1] - indices[0] + 1),
                int(indices[3] - indices[2] + 1),
                int(indices[5] - indices[4] + 1),
            )
            number_cells = math.prod(block_cells)
            variable_bytes = number_cells * header.variable_size
            record_data_bytes = len(header.variables) * variable_bytes
            limits_buffer = stream.read(6 * header.location_size)
            limits = np.frombuffer(limits_buffer, dtype=location_dtype)
            if limits.size != 6:
                raise RuntimeError(f"{path}: truncated MeshBlock coordinate record")
            x_min, x_max, y_min, y_max, z_min, z_max = map(float, limits)
            x_indices = np.flatnonzero((x_axis >= x_min) & (x_axis < x_max))
            y_indices = np.flatnonzero((y_axis >= y_min) & (y_axis < y_max))
            z_indices = np.flatnonzero((z_axis >= z_min) & (z_axis < z_max))
            if not x_indices.size or not y_indices.size or not z_indices.size:
                stream.seek(record_data_bytes, os.SEEK_CUR)
                continue
            x_cell = (x_max - x_min) / block_cells[0]
            y_cell = (y_max - y_min) / block_cells[1]
            z_cell = (z_max - z_min) / block_cells[2]
            sample_x = (x_axis[x_indices] - x_min) / x_cell - 0.5
            sample_y = (y_axis[y_indices] - y_min) / y_cell - 0.5
            sample_z = (z_axis[z_indices] - z_min) / z_cell - 0.5
            zz, yy, xx = np.meshgrid(sample_z, sample_y, sample_x, indexing="ij")
            target = np.ix_(z_indices, y_indices, x_indices)
            data_start = stream.tell()
            for name, variable_index in sorted(
                variable_indices.items(), key=lambda item: item[1]
            ):
                stream.seek(data_start + variable_index * variable_bytes)
                values = np.fromfile(stream, dtype=variable_dtype, count=number_cells)
                if values.size != number_cells:
                    raise RuntimeError(f"{path}: truncated cell data for {name}")
                values = values.reshape(block_cells[2], block_cells[1], block_cells[0])
                order = 0 if name == "gr_excision_mask" else 1
                output[name][target] = map_coordinates(
                    values,
                    (zz, yy, xx),
                    order=order,
                    mode="nearest",
                    prefilter=False,
                ).astype(np.float32, copy=False)
            stream.seek(data_start + record_data_bytes)
    for name, values in output.items():
        if not np.isfinite(values).all():
            missing_cells = int(np.count_nonzero(~np.isfinite(values)))
            raise RuntimeError(f"{path}: {name} has {missing_cells} uncovered voxels")
    return output


def extract_surface(
    values: np.ndarray,
    level: float,
    spacing: float,
    half_width: float,
    step_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    from skimage.measure import marching_cubes

    vertices, faces, _, _ = marching_cubes(
        values,
        level=level,
        spacing=(spacing, spacing, spacing),
        step_size=step_size,
        allow_degenerate=False,
    )
    vertices += -half_width + 0.5 * spacing
    return vertices[:, [2, 1, 0]], faces


def add_surface(axis, vertices, faces, color: str, alpha: float):
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    surface = Poly3DCollection(
        vertices[faces],
        facecolor=color,
        edgecolor="none",
        linewidth=0.0,
        alpha=alpha,
    )
    axis.add_collection3d(surface)
    return surface


def configure_axis(axis, half_width: float, title: str) -> None:
    axis.set_xlim(-half_width, half_width)
    axis.set_ylim(-half_width, half_width)
    axis.set_zlim(-half_width, half_width)
    axis.set_box_aspect((1, 1, 1))
    axis.set_xlabel(r"$(x-x_{\rm BH1})/M$")
    axis.set_ylabel(r"$(y-y_{\rm BH1})/M$")
    axis.set_zlabel(r"$(z-z_{\rm BH1})/M$")
    axis.set_title(title, fontsize=11)
    axis.view_init(elev=22.0, azim=-55.0)
    axis.grid(False)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sigma-file", type=Path, required=True)
    parser.add_argument("--density-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--trajectory-time-offset", type=float, default=0.0)
    parser.add_argument("--follow-bh", choices=(1, 2), type=int, default=1)
    parser.add_argument("--extent", type=float, default=40.0)
    parser.add_argument("--resolution", type=int, default=192)
    parser.add_argument("--marching-step", type=int, default=2)
    parser.add_argument("--density-level", type=float, action="append", default=[])
    parser.add_argument("--rotation-frames", type=int, default=0)
    parser.add_argument("--rotation-fps", type=int, default=30)
    parser.add_argument("--dpi", type=int, default=140)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.extent <= 0.0 or args.resolution < 32 or args.marching_step <= 0:
        raise SystemExit("extent and marching-step must be positive; resolution >= 32")
    if args.rotation_frames < 0 or args.rotation_fps <= 0 or args.dpi <= 0:
        raise SystemExit("invalid rotation or dpi settings")
    sigma_path = args.sigma_file.expanduser().resolve(strict=True)
    density_path = args.density_file.expanduser().resolve(strict=True)
    trajectory_path = args.trajectory.expanduser().resolve(strict=True)
    sigma_header = read_binary_header(sigma_path)
    density_header = read_binary_header(density_path)
    if sigma_header.cycle != density_header.cycle or not math.isclose(
        sigma_header.time, density_header.time, rel_tol=0.0, abs_tol=1.0e-9
    ):
        raise RuntimeError("sigma and density dumps do not have the same time/cycle")
    trajectory_table = load_trajectory_table(trajectory_path)
    trajectory = trajectory_position(
        trajectory_table, sigma_header.time, args.trajectory_time_offset
    )
    followed = trajectory["centers"][args.follow_bh - 1]
    center = tuple(float(value) for value in followed)

    dependency_path = (
        Path(__file__).resolve().parents[2] / "output/postprocess/.jet-vis-deps"
    )
    if dependency_path.exists():
        sys.path.insert(0, str(dependency_path))
    sigma_fields = read_uniform_volume(
        sigma_path,
        ["gr_sigma", "gr_excision_mask"],
        center,
        args.extent,
        args.resolution,
    )
    density_fields = read_uniform_volume(
        density_path,
        ["dens"],
        center,
        args.extent,
        args.resolution,
    )
    sigma = sigma_fields["gr_sigma"].astype(float)
    sigma[sigma_fields["gr_excision_mask"] >= 0.5] = 1.0e-12
    sigma = np.maximum(sigma, 1.0e-12)
    density = np.maximum(density_fields["dens"].astype(float), 1.0e-30)
    spacing = 2.0 * args.extent / args.resolution
    sigma_vertices, sigma_faces = extract_surface(
        np.log10(sigma), 0.0, spacing, args.extent, args.marching_step
    )
    density_levels = args.density_level or [1.0e-3, 1.0e-2, 1.0e-1]
    density_surfaces = [
        (
            level,
            *extract_surface(
                np.log10(density),
                math.log10(level),
                spacing,
                args.extent,
                args.marching_step,
            ),
        )
        for level in density_levels
    ]
    from skimage import __version__ as skimage_version

    import matplotlib

    matplotlib.use("agg")
    import matplotlib.pyplot as plt

    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(figsize=(12.8, 6.4), facecolor="#08111f")
    sigma_axis = figure.add_subplot(1, 2, 1, projection="3d", facecolor="#08111f")
    density_axis = figure.add_subplot(1, 2, 2, projection="3d", facecolor="#08111f")
    for axis in (sigma_axis, density_axis):
        axis.tick_params(colors="white", labelsize=8)
        axis.xaxis.label.set_color("white")
        axis.yaxis.label.set_color("white")
        axis.zaxis.label.set_color("white")
        axis.title.set_color("white")
        for pane in (axis.xaxis.pane, axis.yaxis.pane, axis.zaxis.pane):
            pane.set_facecolor((0.03, 0.07, 0.12, 1.0))
            pane.set_edgecolor((0.35, 0.45, 0.6, 0.3))
    add_surface(sigma_axis, sigma_vertices, sigma_faces, "#7dd3fc", 0.48)
    density_colors = ("#2563eb", "#f59e0b", "#ef4444")
    density_alphas = (0.08, 0.16, 0.32)
    for index, (level, vertices, faces) in enumerate(density_surfaces):
        add_surface(
            density_axis,
            vertices,
            faces,
            density_colors[min(index, len(density_colors) - 1)],
            density_alphas[min(index, len(density_alphas) - 1)],
        )
    relative_centers = [
        np.asarray(center_value, dtype=float) - np.asarray(center, dtype=float)
        for center_value in trajectory["centers"]
    ]
    for axis in (sigma_axis, density_axis):
        for index, relative in enumerate(relative_centers):
            axis.scatter(
                [relative[0]],
                [relative[1]],
                [relative[2]],
                s=34,
                c="white" if index == 0 else "#22d3ee",
                edgecolors="black",
                linewidths=0.5,
                depthshade=False,
            )
    configure_axis(sigma_axis, args.extent, r"jet boundary: $\sigma=1$")
    configure_axis(
        density_axis,
        args.extent,
        "density isosurfaces: " + ", ".join(f"{level:g}" for level in density_levels),
    )
    figure.suptitle(
        f"AthenaK BBH GRMHD | t={sigma_header.time:.6g} M, "
        f"cycle={sigma_header.cycle} | BH{args.follow_bh}-centered",
        color="white",
        fontsize=13,
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    image_path = args.output_dir / "sigma1-density-3d.png"
    figure.savefig(image_path, dpi=args.dpi, facecolor=figure.get_facecolor())

    movie_path = None
    if args.rotation_frames:
        frames_dir = args.output_dir / "rotation-frames"
        frames_dir.mkdir(parents=True, exist_ok=True)
        if any(frames_dir.iterdir()):
            raise RuntimeError(f"refusing non-empty frame directory: {frames_dir}")
        for index in range(args.rotation_frames):
            azimuth = -55.0 + 360.0 * index / args.rotation_frames
            sigma_axis.view_init(elev=22.0, azim=azimuth)
            density_axis.view_init(elev=22.0, azim=azimuth)
            figure.savefig(
                frames_dir / f"frame-{index:06d}.png",
                dpi=args.dpi,
                facecolor=figure.get_facecolor(),
            )
            if index == 0 or (index + 1) % 20 == 0 or index + 1 == args.rotation_frames:
                print(f"rotation {index + 1}/{args.rotation_frames}", flush=True)
        movie_path = args.output_dir / "sigma1-density-3d-rotation.mp4"
        subprocess.run(
            [
                "ffmpeg",
                "-hide_banner",
                "-loglevel",
                "error",
                "-y",
                "-framerate",
                str(args.rotation_fps),
                "-i",
                str(frames_dir / "frame-%06d.png"),
                "-c:v",
                "libx264",
                "-preset",
                "slow",
                "-crf",
                "18",
                "-pix_fmt",
                "yuv420p",
                "-movflags",
                "+faststart",
                str(movie_path),
            ],
            check=True,
        )
    plt.close(figure)

    manifest = {
        "schema": 1,
        "classification": "athenak-bbh-grmhd-3d-isosurface",
        "sigma_input": str(sigma_path),
        "density_input": str(density_path),
        "time_M": sigma_header.time,
        "cycle": sigma_header.cycle,
        "follow_bh": args.follow_bh,
        "coordinate_origin_M": center,
        "extent_M": args.extent,
        "uniform_resolution": args.resolution,
        "uniform_spacing_M": spacing,
        "amr_resampling": "trilinear leaf-block to uniform cell-centered cube",
        "sigma_isosurface": 1.0,
        "density_isosurfaces": density_levels,
        "marching_cubes_step_size": args.marching_step,
        "scikit_image_version": skimage_version,
        "image": str(image_path),
        "image_sha256": file_sha256(image_path),
        "rotation_movie": str(movie_path) if movie_path else None,
        "rotation_movie_sha256": file_sha256(movie_path) if movie_path else None,
        "trajectory": str(trajectory_path),
        "trajectory_time_offset_M": args.trajectory_time_offset,
    }
    manifest_path = args.output_dir / "isosurface-manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(image_path)
    if movie_path:
        print(movie_path)
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
