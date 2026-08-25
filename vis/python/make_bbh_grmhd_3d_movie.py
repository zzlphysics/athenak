#!/usr/bin/env python3
"""Render a time-uniform BH-centered movie of the three-dimensional sigma=1 surface."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
from pathlib import Path, PurePosixPath
import subprocess
import sys

import numpy as np

from make_bbh_grmhd_movie import (
    MovieFrame,
    file_sha256,
    load_trajectory_table,
    select_frames,
    trajectory_position,
)
from render_bbh_grmhd_3d import (
    add_surface,
    configure_axis,
    extract_surface,
    read_uniform_volume,
)


def select_uniform_time_frames(
    frames: list[MovieFrame], requested: int
) -> list[MovieFrame]:
    if requested >= len(frames):
        return frames
    times = np.asarray([frame.time for frame in frames])
    targets = np.linspace(times[0], times[-1], requested)
    selected_indices = []
    for target in targets:
        index = int(np.argmin(np.abs(times - target)))
        if not selected_indices or index != selected_indices[-1]:
            selected_indices.append(index)
    return [frames[index] for index in selected_indices]


def _render_3d_frame(task: dict[str, object]) -> dict[str, object]:
    dependency_path = (
        Path(__file__).resolve().parents[2] / "output/postprocess/.jet-vis-deps"
    )
    if dependency_path.exists():
        sys.path.insert(0, str(dependency_path))
    frame = MovieFrame(
        path=Path(str(task["path"])),
        relative=PurePosixPath(str(task["relative"])),
        size=int(task["size"]),
        sha256=str(task["sha256"]),
        time=float(task["time"]),
        cycle=int(task["cycle"]),
    )
    center_values = task["center"]
    center = tuple(float(value) for value in center_values)
    extent = float(task["extent"])
    resolution = int(task["resolution"])
    spacing = 2.0 * extent / resolution
    fields = read_uniform_volume(
        frame.path,
        ["gr_sigma", "gr_excision_mask"],
        center,
        extent,
        resolution,
    )
    sigma = np.maximum(fields["gr_sigma"].astype(float), 1.0e-12)
    sigma[fields["gr_excision_mask"] >= 0.5] = 1.0e-12
    log_sigma = np.log10(sigma)
    surface_exists = float(np.min(log_sigma)) < 0.0 < float(np.max(log_sigma))
    if surface_exists:
        vertices, faces = extract_surface(
            log_sigma,
            0.0,
            spacing,
            extent,
            int(task["marching_step"]),
        )
    else:
        vertices = np.empty((0, 3))
        faces = np.empty((0, 3), dtype=int)

    import matplotlib

    matplotlib.use("agg")
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(7.0, 6.6), facecolor="#08111f")
    axis = figure.add_subplot(1, 1, 1, projection="3d", facecolor="#08111f")
    axis.tick_params(colors="white", labelsize=8)
    axis.xaxis.label.set_color("white")
    axis.yaxis.label.set_color("white")
    axis.zaxis.label.set_color("white")
    axis.title.set_color("white")
    for pane in (axis.xaxis.pane, axis.yaxis.pane, axis.zaxis.pane):
        pane.set_facecolor((0.03, 0.07, 0.12, 1.0))
        pane.set_edgecolor((0.35, 0.45, 0.6, 0.3))
    if surface_exists:
        add_surface(axis, vertices, faces, "#7dd3fc", 0.52)
    else:
        axis.text2D(
            0.5,
            0.5,
            r"no $\sigma=1$ surface",
            transform=axis.transAxes,
            color="white",
            ha="center",
        )
    centers = task["trajectory_centers"]
    for index, bh_center in enumerate(centers):
        relative = np.asarray(bh_center, dtype=float) - np.asarray(center, dtype=float)
        axis.scatter(
            [relative[0]],
            [relative[1]],
            [relative[2]],
            s=38,
            c="white" if index == int(task["follow_bh"]) - 1 else "#22d3ee",
            edgecolors="black",
            linewidths=0.5,
            depthshade=False,
        )
    configure_axis(axis, extent, r"jet boundary: $\sigma=1$")
    axis.view_init(elev=float(task["elevation"]), azim=float(task["azimuth"]))
    figure.suptitle(
        f"AthenaK BBH GRMHD | t={frame.time:.6g} M, cycle={frame.cycle} | "
        f"BH{task['follow_bh']}-centered",
        color="white",
        fontsize=12,
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    output_path = Path(str(task["output_path"]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=int(task["dpi"]), facecolor=figure.get_facecolor())
    plt.close(figure)
    record = {
        "input": str(frame.path),
        "relative_input": frame.relative.as_posix(),
        "input_bytes": frame.size,
        "input_sha256": frame.sha256,
        "time_M": frame.time,
        "cycle": frame.cycle,
        "coordinate_origin_M": center,
        "surface_exists": surface_exists,
        "sigma_min": float(np.min(sigma)),
        "sigma_max": float(np.max(sigma)),
        "vertices": int(vertices.shape[0]),
        "triangles": int(faces.shape[0]),
        "output": str(output_path),
    }
    output_path.with_suffix(".json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--path-prefix", action="append", default=[])
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--trajectory-time-offset", type=float, default=0.0)
    parser.add_argument("--follow-bh", choices=(1, 2), type=int, default=1)
    parser.add_argument("--time-frames", type=int, default=32)
    parser.add_argument("--extent", type=float, default=40.0)
    parser.add_argument("--resolution", type=int, default=128)
    parser.add_argument("--marching-step", type=int, default=2)
    parser.add_argument("--elevation", type=float, default=22.0)
    parser.add_argument("--azimuth", type=float, default=-55.0)
    parser.add_argument("--fps", type=int, default=8)
    parser.add_argument("--dpi", type=int, default=120)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--verify-sha256", action="store_true")
    parser.add_argument("--movie-name", default="bbh-sigma1-3d-evolution.mp4")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if (
        args.time_frames < 2
        or args.extent <= 0.0
        or args.resolution < 32
        or args.marching_step <= 0
        or args.fps <= 0
        or args.dpi <= 0
        or args.workers <= 0
    ):
        raise SystemExit("invalid positive movie settings")
    dataset = args.dataset.expanduser().resolve(strict=True)
    trajectory_path = args.trajectory.expanduser().resolve(strict=True)
    dependency_path = (
        Path(__file__).resolve().parents[2] / "output/postprocess/.jet-vis-deps"
    )
    if dependency_path.exists():
        sys.path.insert(0, str(dependency_path))
    from skimage import __version__ as skimage_version

    all_frames = select_frames(
        dataset,
        "mhd_gr_diagnostics",
        tuple(args.path_prefix),
        args.verify_sha256,
    )
    frames = select_uniform_time_frames(all_frames, args.time_frames)
    trajectory_table = load_trajectory_table(trajectory_path)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frames_dir = args.output_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    if any(frames_dir.iterdir()):
        raise RuntimeError(f"refusing non-empty frame directory: {frames_dir}")

    tasks = []
    for index, frame in enumerate(frames):
        trajectory = trajectory_position(
            trajectory_table, frame.time, args.trajectory_time_offset
        )
        center = trajectory["centers"][args.follow_bh - 1]
        tasks.append(
            {
                "path": str(frame.path),
                "relative": frame.relative.as_posix(),
                "size": frame.size,
                "sha256": frame.sha256,
                "time": frame.time,
                "cycle": frame.cycle,
                "center": center,
                "trajectory_centers": trajectory["centers"],
                "follow_bh": args.follow_bh,
                "extent": args.extent,
                "resolution": args.resolution,
                "marching_step": args.marching_step,
                "elevation": args.elevation,
                "azimuth": args.azimuth,
                "dpi": args.dpi,
                "output_path": str(frames_dir / f"frame-{index:06d}.png"),
            }
        )
    by_cycle: dict[int, dict[str, object]] = {}
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(_render_3d_frame, task): task for task in tasks}
        for index, future in enumerate(as_completed(futures), start=1):
            record = future.result()
            by_cycle[int(record["cycle"])] = record
            if index == 1 or index % 8 == 0 or index == len(tasks):
                print(f"rendered {index}/{len(tasks)}", flush=True)
    rendered = [by_cycle[frame.cycle] for frame in frames]

    movie_path = args.output_dir / args.movie_name
    subprocess.run(
        [
            "ffmpeg",
            "-hide_banner",
            "-loglevel",
            "error",
            "-y",
            "-framerate",
            str(args.fps),
            "-i",
            str(frames_dir / "frame-%06d.png"),
            "-vf",
            "pad=ceil(iw/2)*2:ceil(ih/2)*2",
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
    manifest = {
        "schema": 1,
        "classification": "athenak-bbh-grmhd-sigma1-3d-time-movie",
        "dataset": str(dataset),
        "path_prefixes": args.path_prefix,
        "input_sha256_verified": args.verify_sha256,
        "selection": "nearest stored dump to uniformly spaced physical times",
        "requested_time_frames": args.time_frames,
        "frames": rendered,
        "follow_bh": args.follow_bh,
        "extent_M": args.extent,
        "uniform_resolution": args.resolution,
        "uniform_spacing_M": 2.0 * args.extent / args.resolution,
        "marching_cubes_step_size": args.marching_step,
        "scikit_image_version": skimage_version,
        "sigma_isosurface": 1.0,
        "camera": {"elevation_deg": args.elevation, "azimuth_deg": args.azimuth},
        "fps": args.fps,
        "duration_seconds": len(rendered) / args.fps,
        "movie": str(movie_path),
        "movie_bytes": movie_path.stat().st_size,
        "movie_sha256": file_sha256(movie_path),
        "trajectory": str(trajectory_path),
        "trajectory_time_offset_M": args.trajectory_time_offset,
    }
    manifest_path = args.output_dir / "movie-manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(movie_path)
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
