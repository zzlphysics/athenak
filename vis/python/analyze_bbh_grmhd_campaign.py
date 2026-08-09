#!/usr/bin/env python3
"""Analyze only checksum-verified, closed BBH GRMHD campaign outputs.

Each input is a completed segment directory produced by ``pull_ready_outputs.py``.
The analysis inventory is built exclusively from file records in local ACKs.  Size
and SHA256 are rechecked before a binary or history file is read, preventing a
partially transferred or subsequently changed output from entering a science plot.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import sys
import tempfile

import numpy as np

from athena_read import hst
from plot_bbh_grmhd import (
    PROXY_WARNING,
    create_dashboard,
    create_history_plot,
    load_trajectory_position,
    read_binary_header,
    read_slice,
)


@dataclass(frozen=True)
class VerifiedFile:
    segment: Path
    relative: PurePosixPath
    path: Path
    size: int
    sha256: str
    ack: Path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_relative(value: object) -> PurePosixPath:
    if not isinstance(value, str):
        raise ValueError(f"manifest path is not a string: {value!r}")
    relative = PurePosixPath(value)
    if (
        relative.is_absolute()
        or not relative.parts
        or any(part in ("", ".", "..") for part in relative.parts)
    ):
        raise ValueError(f"unsafe ACK path: {value!r}")
    return relative


def normalized_file_records(
    payload: dict[str, object], source: Path
) -> dict[str, tuple[int, str]]:
    files = payload.get("files")
    if not isinstance(files, list):
        raise RuntimeError(f"{source}: no file records")
    records: dict[str, tuple[int, str]] = {}
    for record in files:
        if not isinstance(record, dict):
            raise RuntimeError(f"{source}: malformed file record")
        relative = safe_relative(record.get("path"))
        size = int(record.get("size", -1))
        digest = str(record.get("sha256", ""))
        if (
            size < 0
            or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
        ):
            raise RuntimeError(f"{source}: invalid record for {relative}")
        key = relative.as_posix()
        if key in records:
            raise RuntimeError(f"{source}: duplicate record for {relative}")
        records[key] = (size, digest)
    return records


def load_verified_files(segment: Path) -> list[VerifiedFile]:
    """Load immutable file records from every local ACK in one segment."""

    segment = segment.expanduser().resolve(strict=True)
    ack_paths = sorted((segment / ".acks").glob("*.ack"))
    if not ack_paths:
        raise RuntimeError(f"{segment}: no local transfer ACKs found")
    records: dict[str, VerifiedFile] = {}
    for ack_path in ack_paths:
        payload = json.loads(ack_path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise RuntimeError(f"{ack_path}: ACK is not an object")
        if payload.get("schema") != 1:
            raise RuntimeError(f"{ack_path}: unsupported ACK schema")
        manifest_name = payload.get("manifest")
        if (
            not isinstance(manifest_name, str)
            or PurePosixPath(manifest_name).name != manifest_name
            or not manifest_name.endswith(".manifest.ready")
        ):
            raise RuntimeError(f"{ack_path}: unsafe ready-manifest name")
        manifest_path = segment / ".manifests" / manifest_name
        manifest_bytes = manifest_path.read_bytes()
        if hashlib.sha256(manifest_bytes).hexdigest() != payload.get(
            "manifest_sha256"
        ):
            raise RuntimeError(f"{ack_path}: local ready-manifest hash differs")
        manifest_payload = json.loads(manifest_bytes)
        if not isinstance(manifest_payload, dict):
            raise RuntimeError(f"{manifest_path}: manifest is not an object")
        if manifest_payload.get("schema") != 1:
            raise RuntimeError(f"{manifest_path}: unsupported manifest schema")
        if manifest_payload.get("segment") != payload.get("segment"):
            raise RuntimeError(f"{ack_path}: ACK/manifest segment differs")
        ack_records = normalized_file_records(payload, ack_path)
        manifest_records = normalized_file_records(manifest_payload, manifest_path)
        if ack_records != manifest_records:
            raise RuntimeError(f"{ack_path}: ACK/manifest file records differ")
        for relative_text, (size, digest) in ack_records.items():
            relative = PurePosixPath(relative_text)
            path = segment.joinpath(*relative.parts).resolve()
            try:
                path.relative_to(segment)
            except ValueError as exc:
                raise RuntimeError(
                    f"{ack_path}: path escapes segment: {relative}"
                ) from exc
            candidate = VerifiedFile(segment, relative, path, size, digest, ack_path)
            previous = records.get(relative.as_posix())
            if previous is not None and (
                previous.size != candidate.size or previous.sha256 != candidate.sha256
            ):
                raise RuntimeError(
                    f"conflicting ACK records for {relative}: "
                    f"{previous.ack} and {ack_path}"
                )
            records[relative.as_posix()] = candidate
    return list(records.values())


def verify_file(record: VerifiedFile, verify_sha256: bool) -> None:
    if not record.path.is_file() or record.path.is_symlink():
        raise RuntimeError(f"verified output is absent or not regular: {record.path}")
    if record.path.stat().st_size != record.size:
        raise RuntimeError(f"verified output size changed: {record.path}")
    if verify_sha256 and file_sha256(record.path) != record.sha256:
        raise RuntimeError(f"verified output SHA256 changed: {record.path}")


def classify_binary(name: str, variables: tuple[str, ...]) -> str:
    if ".mhd_w_bcc." in name:
        return "mhd_w_bcc"
    if ".mhd_divb." in name or variables == ("divb",):
        return "mhd_divb"
    return "other"


def configured_output_dt(header, stream_name: str) -> float | None:
    for section in header.parameters.values():
        if section.get("file_type") != "bin":
            continue
        if section.get("variable") == stream_name or section.get("id") == stream_name:
            try:
                value = float(section["dt"])
            except (KeyError, ValueError):
                return None
            return value if value > 0.0 and math.isfinite(value) else None
    return None


def meshblock_count(path: Path, header) -> int:
    try:
        cells = math.prod(
            int(header.parameters["meshblock"][f"nx{axis}"]) for axis in (1, 2, 3)
        )
    except KeyError as exc:
        raise RuntimeError(f"{path}: binary header lacks MeshBlock dimensions") from exc
    record_bytes = (
        24
        + 16
        + 6 * header.location_size
        + len(header.variables) * cells * header.variable_size
    )
    payload_bytes = path.stat().st_size - header.data_offset
    count, remainder = divmod(payload_bytes, record_bytes)
    if remainder:
        raise RuntimeError(
            f"{path}: {remainder} bytes remain after fixed MeshBlock records; "
            "file may be incomplete"
        )
    return count


def cadence_summary(frames: list[dict[str, object]]) -> dict[str, object]:
    times = np.asarray(sorted(float(frame["time_M"]) for frame in frames))
    gaps = np.diff(times)
    positive_gaps = gaps[gaps > 0.0]
    result: dict[str, object] = {
        "frames": int(times.size),
        "time_min_M": float(times.min()),
        "time_max_M": float(times.max()),
    }
    if positive_gaps.size:
        result.update(
            {
                "gap_min_M": float(positive_gaps.min()),
                "gap_median_M": float(np.median(positive_gaps)),
                "gap_max_M": float(positive_gaps.max()),
            }
        )
    if np.any(gaps == 0.0):
        result["duplicate_time_records"] = int(np.count_nonzero(gaps == 0.0))
    return result


def storage_projection(
    streams: dict[str, list[dict[str, object]]], target_time: float
) -> dict[str, object]:
    projection: dict[str, object] = {"target_time_M": target_time, "streams": {}}
    total = 0.0
    for name, frames in streams.items():
        if not frames:
            continue
        ordered = sorted(frames, key=lambda frame: float(frame["time_M"]))
        times = np.asarray([float(frame["time_M"]) for frame in ordered])
        sizes = np.asarray([int(frame["bytes"]) for frame in ordered], dtype=float)
        positive_gaps = np.diff(times)
        positive_gaps = positive_gaps[positive_gaps > 0.0]
        configured = [
            float(frame["configured_dt_M"])
            for frame in ordered
            if frame.get("configured_dt_M") is not None
        ]
        if not positive_gaps.size and not configured:
            stream_projection = {
                "status": "needs_at_least_two_frames",
                "median_frame_bytes": float(np.median(sizes)),
            }
        else:
            cadence = (
                float(np.median(positive_gaps))
                if positive_gaps.size
                else float(np.median(configured))
            )
            remaining_frames = max(0, math.ceil((target_time - times.max()) / cadence))
            remaining_bytes = remaining_frames * float(np.median(sizes))
            stream_projection = {
                "status": (
                    "estimated_from_observed_cadence"
                    if positive_gaps.size
                    else "estimated_from_configured_cadence"
                ),
                "cadence_M": cadence,
                "median_frame_bytes": float(np.median(sizes)),
                "remaining_frames": remaining_frames,
                "remaining_bytes": remaining_bytes,
            }
            total += remaining_bytes
        projection["streams"][name] = stream_projection
    projection["remaining_binary_bytes"] = total
    projection["remaining_binary_TiB"] = total / 2**40
    return projection


def merge_histories(paths: list[Path]) -> tuple[list[str], dict[str, np.ndarray]]:
    """Merge restart segments by time; later segment arguments replace overlaps."""

    if not paths:
        raise RuntimeError("no verified history files")
    columns: list[str] | None = None
    rows: dict[float, tuple[float, ...]] = {}
    for path in paths:
        data = hst(str(path))
        current_columns = list(data)
        if columns is None:
            columns = current_columns
        elif current_columns != columns:
            raise RuntimeError(
                f"history columns differ in {path}: {current_columns} != {columns}"
            )
        if not current_columns or current_columns[0] != "time":
            raise RuntimeError(f"{path}: first history column is not time")
        for index, time_value in enumerate(data["time"]):
            rows[float(time_value)] = tuple(
                float(data[column][index]) for column in current_columns
            )
    assert columns is not None
    ordered_rows = [rows[time_value] for time_value in sorted(rows)]
    merged = {
        column: np.asarray([row[index] for row in ordered_rows])
        for index, column in enumerate(columns)
    }
    return columns, merged


def write_history_products(
    output_dir: Path, columns: list[str], data: dict[str, np.ndarray]
) -> tuple[Path, Path]:
    csv_path = output_dir / "merged-history.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*(data[column] for column in columns)))

    hst_path = output_dir / "merged-history.hst"
    with hst_path.open("w", encoding="utf-8") as stream:
        stream.write("# Athena++ history data\n")
        stream.write("# " + " ".join(
            f"[{index}]={column}" for index, column in enumerate(columns)
        ) + "\n")
        for row in zip(*(data[column] for column in columns)):
            stream.write(" ".join(f"{value:.17e}" for value in row) + "\n")
    return csv_path, hst_path


def render_timeline(
    streams: dict[str, list[dict[str, object]]], output_path: Path, dpi: int
) -> None:
    cache = Path(tempfile.gettempdir()) / "athenak-matplotlib-cache"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib

    matplotlib.use("agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(11, 3.8))
    for row, (name, frames) in enumerate(sorted(streams.items())):
        if frames:
            times = [float(frame["time_M"]) for frame in frames]
            axis.scatter(times, [row] * len(times), s=28, label=name)
    axis.set_xlabel("t / M")
    axis.set_ylabel("verified output stream")
    axis.set_yticks(range(len(streams)), sorted(streams))
    axis.grid(axis="x", alpha=0.3)
    axis.set_title("Checksum-verified closed-file output cadence")
    figure.tight_layout()
    figure.savefig(output_path, dpi=dpi)
    plt.close(figure)


def selected_frames(
    frames: list[dict[str, object]], every: int
) -> list[dict[str, object]]:
    ordered = sorted(frames, key=lambda frame: float(frame["time_M"]))
    if every <= 0 or not ordered:
        return []
    selected = ordered[::every]
    if selected[-1] is not ordered[-1]:
        selected.append(ordered[-1])
    return selected


def render_dashboards(
    streams: dict[str, list[dict[str, object]]],
    output_dir: Path,
    every: int,
    plane: str,
    location: float,
    extent: float,
    dpi: int,
    trajectory_path: Path | None,
) -> list[dict[str, object]]:
    rendered: list[dict[str, object]] = []
    choices = {
        "mhd_w_bcc": ["dens", "press", "temperature", "bmag", "beta_inv_proxy", "level"],
        "mhd_divb": ["divb", "level"],
    }
    for stream_name, panel_names in choices.items():
        for frame in selected_frames(streams.get(stream_name, []), every):
            input_path = Path(str(frame["path"]))
            slice_data = read_slice(input_path, panel_names, plane, location, extent)
            trajectory = None
            if trajectory_path is not None:
                offset = float(
                    slice_data.header.parameters.get("problem", {}).get(
                        "trajectory_time_offset", "0"
                    )
                )
                trajectory = load_trajectory_position(
                    trajectory_path, slice_data.header.time, offset
                )
            frame_dir = output_dir / "frames" / stream_name
            figure_path = frame_dir / (
                f"t{slice_data.header.time:012.6f}_c{slice_data.header.cycle:08d}_"
                f"{plane}{location:g}.png"
            )
            numerical = create_dashboard(
                slice_data,
                panel_names,
                figure_path,
                extent,
                dpi,
                False,
                trajectory,
                1.0e-8,
            )
            record = {
                "input": str(input_path),
                "output": str(figure_path),
                "stream": stream_name,
                "time_M": slice_data.header.time,
                "cycle": slice_data.header.cycle,
                "panels": panel_names,
                **numerical,
            }
            figure_path.with_suffix(".json").write_text(
                json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            rendered.append(record)
            print(figure_path)
    return rendered


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("segments", nargs="+", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--no-sha256",
        action="store_true",
        help="trust ACK size only (unsafe and intended only for quick local debugging)",
    )
    parser.add_argument(
        "--render-every",
        type=int,
        default=1,
        help="render every Nth frame plus the last; 0 builds inventory/history only",
    )
    parser.add_argument("--target-time", type=float, default=3500.0)
    parser.add_argument("--plane", choices=("x", "y", "z"), default="z")
    parser.add_argument("--location", type=float, default=0.0)
    parser.add_argument("--extent", type=float, default=80.0)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--trajectory", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.render_every < 0:
        raise SystemExit("--render-every must be non-negative")
    if args.extent <= 0.0:
        raise SystemExit("--extent must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_records: list[VerifiedFile] = []
    for segment in args.segments:
        all_records.extend(load_verified_files(segment))
    science_records = [
        record
        for record in all_records
        if record.path.suffix == ".bin" or record.path.suffix == ".hst"
    ]
    if not science_records:
        raise RuntimeError("ACKs contain no binary or history science output")
    for record in science_records:
        verify_file(record, not args.no_sha256)
        print(f"verified {record.relative} ({record.size} bytes)", file=sys.stderr)

    streams: dict[str, list[dict[str, object]]] = {
        "mhd_w_bcc": [],
        "mhd_divb": [],
        "other": [],
    }
    for record in science_records:
        if record.path.suffix != ".bin":
            continue
        header = read_binary_header(record.path)
        stream_name = classify_binary(record.path.name, header.variables)
        streams[stream_name].append(
            {
                "path": str(record.path),
                "segment": str(record.segment),
                "ack": str(record.ack),
                "bytes": record.size,
                "sha256": record.sha256,
                "time_M": header.time,
                "cycle": header.cycle,
                "variables": list(header.variables),
                "meshblocks": meshblock_count(record.path, header),
                "configured_dt_M": configured_output_dt(header, stream_name),
            }
        )
    for frames in streams.values():
        frames.sort(key=lambda frame: (float(frame["time_M"]), str(frame["path"])))
    streams = {name: frames for name, frames in streams.items() if frames}

    history_records = [
        record for record in science_records if record.path.suffix == ".hst"
    ]
    # Preserve command-line segment order so later restart segments replace overlap.
    history_paths = [record.path for record in history_records]
    history_summary = None
    if history_paths:
        columns, history_data = merge_histories(history_paths)
        csv_path, hst_path = write_history_products(
            args.output_dir, columns, history_data
        )
        history_summary = create_history_plot(
            hst_path, args.output_dir / "merged-history.png", args.dpi
        )
        history_summary["csv"] = str(csv_path)
        history_summary["verified_sources"] = [str(path) for path in history_paths]

    timeline_path = args.output_dir / "verified-output-timeline.png"
    render_timeline(streams, timeline_path, args.dpi)
    rendered = render_dashboards(
        streams,
        args.output_dir,
        args.render_every,
        args.plane,
        args.location,
        args.extent,
        args.dpi,
        args.trajectory.expanduser().resolve(strict=True) if args.trajectory else None,
    )
    report = {
        "schema": 1,
        "classification": "athenak-bbh-grmhd-verified-campaign-analysis",
        "segments": [str(path.expanduser().resolve()) for path in args.segments],
        "input_policy": {
            "local_ack_required": True,
            "size_verified": True,
            "sha256_verified": not args.no_sha256,
            "active_or_partial_files_accepted": False,
        },
        "streams": streams,
        "cadence": {
            name: cadence_summary(frames) for name, frames in streams.items()
        },
        "storage_projection": storage_projection(streams, args.target_time),
        "history": history_summary,
        "timeline": str(timeline_path),
        "rendered_frames": rendered,
        "publication_readiness": {
            "exact_outputs": [
                "density",
                "pressure",
                "temperature",
                "coordinate B",
                "primitive velocity",
                "divB",
            ],
            "proxy_warning": PROXY_WARNING,
            "not_available_from_current_dump": [
                "covariant magnetization and plasma beta",
                "MRI quality factors Q_theta and Q_phi",
                "moving-horizon mass and magnetic fluxes",
                "covariant Bernoulli parameter and Lorentz factor",
            ],
        },
    }
    report_path = args.output_dir / "campaign-analysis.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(report_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.stderr.close()
        raise SystemExit(1)
