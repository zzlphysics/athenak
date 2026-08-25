#!/usr/bin/env python3
"""Render a fixed-color AMR slice that follows one BBH trajectory center."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
from pathlib import Path, PurePosixPath

from make_bbh_grmhd_movie import (
    MovieFrame,
    encode_movie,
    file_sha256,
    load_trajectory_table,
    select_frames,
    trajectory_position,
)
from plot_bbh_grmhd import (
    PANELS,
    create_dashboard,
    read_binary_header,
    read_linearly_interpolated_slice,
)


def parse_panel_limits(
    specifications: list[str], panel_names: list[str]
) -> dict[str, tuple[float, float]]:
    limits: dict[str, tuple[float, float]] = {}
    for specification in specifications:
        try:
            name, bounds = specification.split("=", 1)
            lower_text, upper_text = bounds.split(",", 1)
            lower, upper = float(lower_text), float(upper_text)
        except ValueError as exc:
            raise ValueError(
                f"invalid --panel-limit {specification!r}; expected name=min,max"
            ) from exc
        if name not in panel_names:
            raise ValueError(f"limit supplied for an unselected panel: {name}")
        if name in limits:
            raise ValueError(f"duplicate panel limit: {name}")
        if not (math.isfinite(lower) and math.isfinite(upper) and upper > lower):
            raise ValueError(f"invalid panel limit for {name}: {lower}, {upper}")
        if PANELS[name].scale == "log" and lower <= 0.0:
            raise ValueError(f"log panel {name} requires a positive lower limit")
        limits[name] = (lower, upper)
    missing = sorted(set(panel_names).difference(limits))
    if missing:
        raise ValueError(
            "a following-plane movie requires explicit fixed limits for every panel; "
            f"missing: {', '.join(missing)}"
        )
    return limits


def plane_geometry(
    center: list[float], plane: str
) -> tuple[float, tuple[float, float]]:
    if plane == "x":
        return center[0], (center[1], center[2])
    if plane == "y":
        return center[1], (center[0], center[2])
    return center[2], (center[0], center[1])


def _render_following_frame(task: dict[str, object]) -> dict[str, object]:
    frame = MovieFrame(
        path=Path(str(task["path"])),
        relative=PurePosixPath(str(task["relative"])),
        size=int(task["size"]),
        sha256=str(task["sha256"]),
        time=float(task["time"]),
        cycle=int(task["cycle"]),
    )
    panel_names = [str(name) for name in task["panel_names"]]
    fixed_limits = {
        str(name): (float(bounds[0]), float(bounds[1]))
        for name, bounds in dict(task["fixed_limits"]).items()
    }
    origin_values = task["coordinate_origin"]
    coordinate_origin = (float(origin_values[0]), float(origin_values[1]))
    slice_data = read_linearly_interpolated_slice(
        frame.path,
        panel_names,
        str(task["plane"]),
        float(task["location"]),
        float(task["extent"]),
        coordinate_origin,
        int(task["raster_resolution"]),
    )
    output_path = Path(str(task["output_path"]))
    summary = create_dashboard(
        slice_data,
        panel_names,
        output_path,
        float(task["extent"]),
        int(task["dpi"]),
        bool(task["show_grid"]),
        dict(task["trajectory"]),
        float(task["rho_mask_fraction"]),
        fixed_limits=fixed_limits,
        coordinate_origin=coordinate_origin,
        coordinate_origin_label=str(task["coordinate_origin_label"]),
    )
    record = {
        "input": str(frame.path),
        "relative_input": frame.relative.as_posix(),
        "input_bytes": frame.size,
        "input_sha256": frame.sha256,
        "output": str(output_path),
        "time_M": frame.time,
        "cycle": frame.cycle,
        "slice_location_M": float(task["location"]),
        "coordinate_origin_M": coordinate_origin,
        "plane_sampling": "amr-rasterized-two-cell-center-linear-interpolation",
        "raster_resolution": int(task["raster_resolution"]),
        **summary,
    }
    output_path.with_suffix(".json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--stream",
        choices=("mhd_w_bcc", "mhd_gr_diagnostics"),
        required=True,
    )
    parser.add_argument("--path-prefix", action="append", default=[])
    parser.add_argument("--panels", required=True)
    parser.add_argument(
        "--panel-limit",
        action="append",
        default=[],
        help="fixed display range as panel=min,max; repeat for every panel",
    )
    parser.add_argument("--plane", choices=("x", "y", "z"), default="y")
    parser.add_argument("--follow-bh", choices=(1, 2), type=int, default=1)
    parser.add_argument("--extent", type=float, default=40.0)
    parser.add_argument("--raster-resolution", type=int, default=512)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--trajectory-time-offset", type=float)
    parser.add_argument("--rho-mask-fraction", type=float, default=1.0e-8)
    parser.add_argument("--frame-step", type=int, default=1)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--dpi", type=int, default=120)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--grid", action="store_true")
    parser.add_argument("--verify-sha256", action="store_true")
    parser.add_argument("--movie-name", default="bbh-following-plane.mp4")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    panel_names = [name.strip() for name in args.panels.split(",") if name.strip()]
    invalid = sorted(set(panel_names).difference(PANELS))
    if invalid or not panel_names:
        raise SystemExit(f"invalid or empty panels: {', '.join(invalid)}")
    try:
        fixed_limits = parse_panel_limits(args.panel_limit, panel_names)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    if (
        args.extent <= 0.0
        or args.raster_resolution < 32
        or args.frame_step <= 0
        or args.fps <= 0
        or args.dpi <= 0
        or args.workers <= 0
    ):
        raise SystemExit(
            "extent, frame-step, fps, dpi, and workers must be positive; "
            "raster-resolution must be at least 32"
        )
    if args.rho_mask_fraction < 0.0:
        raise SystemExit("rho-mask-fraction must be non-negative")
    if Path(args.movie_name).name != args.movie_name or not args.movie_name.endswith(
        ".mp4"
    ):
        raise SystemExit("--movie-name must be a plain .mp4 filename")

    dataset = args.dataset.expanduser().resolve(strict=True)
    trajectory_path = args.trajectory.expanduser().resolve(strict=True)
    frames = select_frames(
        dataset,
        args.stream,
        tuple(args.path_prefix),
        args.verify_sha256,
    )[:: args.frame_step]
    trajectory_columns = load_trajectory_table(trajectory_path)
    first_header = read_binary_header(frames[0].path)
    if args.trajectory_time_offset is None:
        time_offset = float(
            first_header.parameters.get("problem", {}).get(
                "trajectory_time_offset", "0"
            )
        )
    else:
        time_offset = args.trajectory_time_offset

    args.output_dir.mkdir(parents=True, exist_ok=True)
    frames_dir = args.output_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    if any(frames_dir.iterdir()):
        raise RuntimeError(f"refusing non-empty frame directory: {frames_dir}")
    limits_payload = {
        "schema": 1,
        "policy": "user-specified fixed range for every panel and frame",
        "stream": args.stream,
        "panels": {
            name: {
                "display_min": bounds[0],
                "display_max": bounds[1],
                "scale": PANELS[name].scale,
            }
            for name, bounds in fixed_limits.items()
        },
        "frames": len(frames),
        "cycle_min": frames[0].cycle,
        "cycle_max": frames[-1].cycle,
        "time_min_M": frames[0].time,
        "time_max_M": frames[-1].time,
        "plane_sampling": (
            "trajectory-following-amr-rasterized-two-cell-center-linear-interpolation"
        ),
        "follow_bh": args.follow_bh,
        "extent_M": args.extent,
        "raster_resolution": args.raster_resolution,
    }
    limits_path = args.output_dir / "fixed-color-limits.json"
    limits_path.write_text(
        json.dumps(limits_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    tasks: list[dict[str, object]] = []
    for index, frame in enumerate(frames):
        trajectory = trajectory_position(trajectory_columns, frame.time, time_offset)
        center = trajectory["centers"][args.follow_bh - 1]
        location, coordinate_origin = plane_geometry(center, args.plane)
        tasks.append(
            {
                "path": str(frame.path),
                "relative": frame.relative.as_posix(),
                "size": frame.size,
                "sha256": frame.sha256,
                "time": frame.time,
                "cycle": frame.cycle,
                "panel_names": panel_names,
                "fixed_limits": fixed_limits,
                "plane": args.plane,
                "location": location,
                "extent": args.extent,
                "raster_resolution": args.raster_resolution,
                "dpi": args.dpi,
                "show_grid": args.grid,
                "trajectory": trajectory,
                "coordinate_origin": coordinate_origin,
                "coordinate_origin_label": f"BH{args.follow_bh}",
                "rho_mask_fraction": args.rho_mask_fraction,
                "output_path": str(frames_dir / f"frame-{index:06d}.png"),
            }
        )

    rendered: list[dict[str, object]] = []
    if args.workers == 1:
        for index, task in enumerate(tasks, start=1):
            rendered.append(_render_following_frame(task))
            if index == 1 or index % 25 == 0 or index == len(tasks):
                print(f"rendered {index}/{len(tasks)}", flush=True)
    else:
        by_cycle: dict[int, dict[str, object]] = {}
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(_render_following_frame, task): task for task in tasks
            }
            for index, future in enumerate(as_completed(futures), start=1):
                record = future.result()
                by_cycle[int(record["cycle"])] = record
                if index == 1 or index % 25 == 0 or index == len(tasks):
                    print(f"rendered {index}/{len(tasks)}", flush=True)
        rendered = [by_cycle[frame.cycle] for frame in frames]

    movie_path = args.output_dir / args.movie_name
    encode_movie(frames_dir, args.fps, movie_path)
    manifest = {
        "schema": 1,
        "classification": "athenak-bbh-grmhd-following-plane-fixed-color-movie",
        "dataset": str(dataset),
        "stream": args.stream,
        "path_prefixes": args.path_prefix,
        "input_sha256_verified": args.verify_sha256,
        "panels": panel_names,
        "fixed_color_limits": str(limits_path),
        "plane": args.plane,
        "plane_sampling": limits_payload["plane_sampling"],
        "follow_bh": args.follow_bh,
        "coordinate_system": f"camera centered on BH{args.follow_bh}",
        "extent_M": args.extent,
        "raster_resolution": args.raster_resolution,
        "frame_step": args.frame_step,
        "frames": rendered,
        "fps": args.fps,
        "duration_seconds": len(rendered) / args.fps,
        "movie": str(movie_path),
        "movie_bytes": movie_path.stat().st_size,
        "movie_sha256": file_sha256(movie_path),
        "trajectory": str(trajectory_path),
        "trajectory_time_offset_M": time_offset,
        "ffmpeg_codec": "libx264 crf=18 preset=slow pix_fmt=yuv420p",
    }
    manifest_path = args.output_dir / "movie-manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(movie_path)
    print(manifest_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        raise SystemExit(1)
