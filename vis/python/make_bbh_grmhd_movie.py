#!/usr/bin/env python3
"""Render a checksum-bound BBH GRMHD movie with one color scale for all frames."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess

import numpy as np

from analyze_bbh_grmhd_campaign import (
    classify_binary,
    is_athenak_binary_output,
    load_verified_files,
    verify_file,
)
from plot_bbh_grmhd import (
    PANELS,
    calculate_limits,
    create_dashboard,
    panel_values,
    read_binary_header,
    read_cell_face_interpolated_slice,
    read_slice,
)


@dataclass(frozen=True)
class MovieFrame:
    path: Path
    relative: PurePosixPath
    size: int
    sha256: str
    time: float
    cycle: int


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def select_frames(
    dataset: Path,
    stream_name: str,
    path_prefixes: tuple[str, ...],
    verify_sha256: bool,
) -> list[MovieFrame]:
    records = load_verified_files(dataset)
    token = f".{stream_name}."
    candidates = [
        record
        for record in records
        if record.path.suffix == ".bin"
        and token in record.path.name
        and (
            not path_prefixes
            or any(
                record.relative.as_posix().startswith(prefix)
                for prefix in path_prefixes
            )
        )
    ]
    if not candidates:
        raise RuntimeError(
            f"no ACK-bound {stream_name} files match the requested paths"
        )

    frames_by_cycle: dict[int, MovieFrame] = {}
    for record in candidates:
        verify_file(record, verify_sha256)
        if not is_athenak_binary_output(record.path):
            raise RuntimeError(f"not an AthenaK binary output: {record.path}")
        header = read_binary_header(record.path)
        classified = classify_binary(record.path.name, header.variables)
        if classified != stream_name:
            raise RuntimeError(
                f"{record.path}: expected {stream_name}, classified as {classified}"
            )
        frame = MovieFrame(
            record.path,
            record.relative,
            record.size,
            record.sha256,
            header.time,
            header.cycle,
        )
        previous = frames_by_cycle.get(frame.cycle)
        if previous is not None:
            raise RuntimeError(
                f"duplicate cycle {frame.cycle}: {previous.relative} "
                f"and {frame.relative}"
            )
        frames_by_cycle[frame.cycle] = frame

    frames = sorted(
        frames_by_cycle.values(), key=lambda frame: (frame.cycle, frame.time)
    )
    for previous, current in zip(frames, frames[1:]):
        if current.cycle <= previous.cycle or current.time <= previous.time:
            raise RuntimeError(
                "movie frames are not strictly increasing in cycle and time"
            )
        if (
            stream_name.startswith("bbh_local_")
            and current.cycle != previous.cycle + 1
        ):
            raise RuntimeError(
                "movie timeline has a root-cycle gap: "
                f"{previous.cycle} -> {current.cycle}"
            )
    return frames


def _sample_finite(values: list[np.ndarray], scale: str, limit: int) -> np.ndarray:
    finite_parts = [value[np.isfinite(value)].ravel() for value in values]
    finite_parts = [value for value in finite_parts if value.size]
    if not finite_parts:
        return np.empty(0, dtype=float)
    finite = np.concatenate(finite_parts)
    if scale == "log":
        finite = finite[finite > 0.0]
    if finite.size > limit:
        indices = np.linspace(0, finite.size - 1, limit, dtype=np.int64)
        finite = finite[indices]
    return finite.astype(float, copy=False)


def fixed_color_limits(
    frames: list[MovieFrame],
    panel_names: list[str],
    plane: str,
    location: float,
    half_width: float | None,
    rho_mask_fraction: float,
    samples_per_frame: int,
    interpolate_plane: bool,
) -> dict[str, tuple[float, float]]:
    samples: dict[str, list[np.ndarray]] = {name: [] for name in panel_names}
    for index, frame in enumerate(frames, start=1):
        reader = read_cell_face_interpolated_slice if interpolate_plane else read_slice
        slice_data = reader(frame.path, panel_names, plane, location, half_width)
        if any(PANELS[name].proxy for name in panel_names):
            density_max = max(
                float(np.nanmax(block.fields["dens"])) for block in slice_data.blocks
            )
            density_threshold = (
                density_max * rho_mask_fraction if rho_mask_fraction > 0.0 else None
            )
        else:
            density_threshold = None
        for panel_name in panel_names:
            panel = PANELS[panel_name]
            sampled = _sample_finite(
                panel_values(slice_data, panel, density_threshold),
                panel.scale,
                samples_per_frame,
            )
            if sampled.size:
                samples[panel_name].append(sampled)
        if index == 1 or index % 50 == 0 or index == len(frames):
            print(f"color-limit pass {index}/{len(frames)}", flush=True)

    limits: dict[str, tuple[float, float]] = {}
    for panel_name in panel_names:
        if not samples[panel_name]:
            raise RuntimeError(f"no finite samples for panel {panel_name}")
        limits[panel_name] = calculate_limits(
            samples[panel_name], PANELS[panel_name].scale
        )
    return limits


def load_trajectory_table(path: Path | None) -> np.ndarray | None:
    if path is None:
        return None
    columns = np.loadtxt(path, comments="#", usecols=range(7), ndmin=2)
    if columns.shape[0] < 2 or np.any(np.diff(columns[:, 0]) <= 0.0):
        raise RuntimeError(f"invalid trajectory table: {path}")
    return columns


def trajectory_position(
    columns: np.ndarray | None, simulation_time: float, time_offset: float
) -> dict[str, object] | None:
    if columns is None:
        return None
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


def _render_frame(task: dict[str, object]) -> dict[str, object]:
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
    half_width_value = task["half_width"]
    half_width = None if half_width_value is None else float(half_width_value)
    reader = (
        read_cell_face_interpolated_slice
        if bool(task["interpolate_plane"])
        else read_slice
    )
    slice_data = reader(
        frame.path,
        panel_names,
        str(task["plane"]),
        float(task["location"]),
        half_width,
    )
    output_path = Path(str(task["output_path"]))
    summary = create_dashboard(
        slice_data,
        panel_names,
        output_path,
        half_width,
        int(task["dpi"]),
        bool(task["show_grid"]),
        task["trajectory"],
        float(task["rho_mask_fraction"]),
        fixed_limits=fixed_limits,
    )
    record = {
        "input": str(frame.path),
        "relative_input": frame.relative.as_posix(),
        "input_bytes": frame.size,
        "input_sha256": frame.sha256,
        "output": str(output_path),
        "time_M": frame.time,
        "cycle": frame.cycle,
        **summary,
    }
    output_path.with_suffix(".json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return record


def encode_movie(frames_dir: Path, fps: int, output_path: Path) -> None:
    command = [
        "ffmpeg",
        "-hide_banner",
        "-loglevel",
        "error",
        "-y",
        "-framerate",
        str(fps),
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
        str(output_path),
    ]
    subprocess.run(command, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--stream",
        choices=(
            "bbh_local_w",
            "bbh_local_gr",
            "mhd_w_bcc",
            "mhd_gr_diagnostics",
        ),
        required=True,
    )
    parser.add_argument(
        "--path-prefix",
        action="append",
        default=[],
        help="ACK-relative production prefix; repeat to form a continuous chain",
    )
    parser.add_argument("--panels", required=True)
    parser.add_argument("--plane", choices=("x", "y", "z"), default="z")
    parser.add_argument("--location", type=float, default=0.0)
    parser.add_argument(
        "--extent",
        type=float,
        default=0.0,
        help="plot half-width in M; 0 uses the complete stored local region",
    )
    parser.add_argument("--trajectory", type=Path)
    parser.add_argument("--trajectory-time-offset", type=float)
    parser.add_argument("--rho-mask-fraction", type=float, default=1.0e-8)
    parser.add_argument("--samples-per-frame", type=int, default=4096)
    parser.add_argument("--frame-step", type=int, default=1)
    parser.add_argument("--fps", type=int, default=24)
    parser.add_argument("--dpi", type=int, default=120)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--grid", action="store_true")
    parser.add_argument(
        "--interpolate-plane",
        action="store_true",
        help="average the two full-3D cell-center planes bracketing a cell face",
    )
    parser.add_argument("--verify-sha256", action="store_true")
    parser.add_argument("--movie-name", default="bbh-grmhd.mp4")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    panel_names = [name.strip() for name in args.panels.split(",") if name.strip()]
    invalid = sorted(set(panel_names).difference(PANELS))
    if invalid or not panel_names:
        raise SystemExit(f"invalid or empty panels: {', '.join(invalid)}")
    if args.frame_step <= 0 or args.fps <= 0 or args.dpi <= 0 or args.workers <= 0:
        raise SystemExit("frame-step, fps, dpi, and workers must be positive")
    if args.samples_per_frame < 128:
        raise SystemExit("--samples-per-frame must be at least 128")
    if args.extent < 0.0 or args.rho_mask_fraction < 0.0:
        raise SystemExit("extent and rho-mask-fraction must be non-negative")
    if Path(args.movie_name).name != args.movie_name or not args.movie_name.endswith(
        ".mp4"
    ):
        raise SystemExit("--movie-name must be a plain .mp4 filename")

    dataset = args.dataset.expanduser().resolve(strict=True)
    frames = select_frames(
        dataset,
        args.stream,
        tuple(args.path_prefix),
        args.verify_sha256,
    )
    frames = frames[:: args.frame_step]
    if not frames:
        raise RuntimeError("frame selection is empty")
    half_width = None if args.extent == 0.0 else args.extent

    args.output_dir.mkdir(parents=True, exist_ok=True)
    frames_dir = args.output_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    if any(frames_dir.iterdir()):
        raise RuntimeError(f"refusing non-empty frame directory: {frames_dir}")

    limits = fixed_color_limits(
        frames,
        panel_names,
        args.plane,
        args.location,
        half_width,
        args.rho_mask_fraction,
        args.samples_per_frame,
        args.interpolate_plane,
    )
    limits_payload = {
        "schema": 1,
        "policy": (
            "one robust sampled percentile range per panel for the complete movie"
        ),
        "stream": args.stream,
        "panels": {
            name: {
                "display_min": bounds[0],
                "display_max": bounds[1],
                "scale": PANELS[name].scale,
            }
            for name, bounds in limits.items()
        },
        "frames": len(frames),
        "cycle_min": frames[0].cycle,
        "cycle_max": frames[-1].cycle,
        "time_min_M": frames[0].time,
        "time_max_M": frames[-1].time,
        "samples_per_frame_per_panel": args.samples_per_frame,
        "rho_mask_fraction": args.rho_mask_fraction,
        "plane_sampling": (
            "two-sided-cell-center-linear-interpolation"
            if args.interpolate_plane
            else "stored-or-nearest-cell"
        ),
    }
    limits_path = args.output_dir / "fixed-color-limits.json"
    limits_path.write_text(
        json.dumps(limits_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    trajectory_columns = load_trajectory_table(
        args.trajectory.expanduser().resolve(strict=True) if args.trajectory else None
    )
    first_header = read_binary_header(frames[0].path)
    if args.trajectory_time_offset is None:
        time_offset = float(
            first_header.parameters.get("problem", {}).get(
                "trajectory_time_offset", "0"
            )
        )
    else:
        time_offset = args.trajectory_time_offset

    tasks = []
    for index, frame in enumerate(frames):
        tasks.append(
            {
                "path": str(frame.path),
                "relative": frame.relative.as_posix(),
                "size": frame.size,
                "sha256": frame.sha256,
                "time": frame.time,
                "cycle": frame.cycle,
                "panel_names": panel_names,
                "fixed_limits": limits,
                "plane": args.plane,
                "location": args.location,
                "half_width": half_width,
                "dpi": args.dpi,
                "show_grid": args.grid,
                "interpolate_plane": args.interpolate_plane,
                "trajectory": trajectory_position(
                    trajectory_columns, frame.time, time_offset
                ),
                "rho_mask_fraction": args.rho_mask_fraction,
                "output_path": str(frames_dir / f"frame-{index:06d}.png"),
            }
        )

    rendered: list[dict[str, object]] = []
    if args.workers == 1:
        for index, task in enumerate(tasks, start=1):
            rendered.append(_render_frame(task))
            if index == 1 or index % 25 == 0 or index == len(tasks):
                print(f"rendered {index}/{len(tasks)}", flush=True)
    else:
        by_cycle: dict[int, dict[str, object]] = {}
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(_render_frame, task): task for task in tasks}
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
        "classification": "athenak-bbh-grmhd-fixed-color-movie",
        "dataset": str(dataset),
        "stream": args.stream,
        "path_prefixes": args.path_prefix,
        "input_sha256_verified": args.verify_sha256,
        "panels": panel_names,
        "fixed_color_limits": str(limits_path),
        "plane_sampling": limits_payload["plane_sampling"],
        "frame_step": args.frame_step,
        "frames": rendered,
        "fps": args.fps,
        "duration_seconds": len(rendered) / args.fps,
        "movie": str(movie_path),
        "movie_bytes": movie_path.stat().st_size,
        "movie_sha256": file_sha256(movie_path),
        "trajectory": str(args.trajectory.resolve()) if args.trajectory else None,
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
