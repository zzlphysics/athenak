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
import struct
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
    gr_diagnostics = (
        "gr_bsq",
        "gr_lorentz",
        "gr_sigma",
        "gr_beta_inv",
    )
    if ".bbh_local_w." in name:
        return "bbh_local_w"
    if ".bbh_local_gr." in name:
        return "bbh_local_gr"
    if ".mhd_w_bcc." in name:
        return "mhd_w_bcc"
    if (
        ".mhd_gr_diagnostics." in name
        or variables == gr_diagnostics
        or variables == gr_diagnostics + ("gr_excision_mask",)
    ):
        return "mhd_gr_diagnostics"
    if ".mhd_divb." in name or variables == ("divb",):
        return "mhd_divb"
    return "other"


def classify_history(name: str) -> str:
    """Return the AthenaK physics suffix of a history filename."""

    for physics in ("user", "mhd", "hydro", "z4c"):
        if name.endswith(f".{physics}.hst"):
            return physics
    return "other"


def is_event_log(path: Path) -> bool:
    with path.open("r", encoding="utf-8") as stream:
        return stream.readline() == "# Athena event counter data\n"


def configured_output_dt(header, stream_name: str) -> float | None:
    # Prefer the explicit id because local and global streams may contain the same
    # variable group at different cadences.
    sections = [
        section
        for section in header.parameters.values()
        if section.get("file_type") == "bin"
    ]
    matches = [section for section in sections if section.get("id") == stream_name]
    if not matches:
        matches = [
            section for section in sections if section.get("variable") == stream_name
        ]
    if len(matches) != 1:
        return None
    try:
        value = float(matches[0]["dt"])
    except (KeyError, ValueError):
        return None
    return value if value > 0.0 and math.isfinite(value) else None


def configured_output_dcycle(header, stream_name: str) -> int | None:
    sections = [
        section
        for section in header.parameters.values()
        if section.get("file_type") == "bin"
    ]
    matches = [section for section in sections if section.get("id") == stream_name]
    if not matches:
        matches = [
            section for section in sections if section.get("variable") == stream_name
        ]
    if len(matches) != 1:
        return None
    try:
        value = int(matches[0]["dcycle"])
    except (KeyError, ValueError):
        return None
    return value if value > 0 else None


def configured_output_region(header, stream_name: str) -> dict[str, object] | None:
    for section in header.parameters.values():
        if section.get("file_type") != "bin":
            continue
        if section.get("id") != stream_name or "region_center" not in section:
            continue
        try:
            default_width = float(section.get("region_half_width", "0"))
            widths = [
                float(section.get(f"region_half_width{axis}", default_width))
                for axis in (1, 2, 3)
            ]
            slice_axis = int(section.get("region_slice_axis", "0"))
            slice_offset = float(section.get("region_slice_offset", "0"))
        except ValueError:
            return None
        return {
            "center": section["region_center"],
            "half_width_M": widths,
            "slice_axis": slice_axis,
            "slice_offset_M": slice_offset,
        }
    return None


def meshblock_level_counts(path: Path, header) -> dict[int, int]:
    levels: dict[int, int] = {}
    with path.open("rb") as stream:
        stream.seek(header.data_offset)
        file_size = path.stat().st_size
        while stream.tell() < file_size:
            indices_buffer = stream.read(24)
            if len(indices_buffer) != 24:
                raise RuntimeError(f"{path}: truncated MeshBlock index record")
            indices = struct.unpack("=6i", indices_buffer)
            shape = (
                indices[1]-indices[0]+1,
                indices[3]-indices[2]+1,
                indices[5]-indices[4]+1,
            )
            if any(cells <= 0 for cells in shape):
                raise RuntimeError(f"{path}: invalid MeshBlock output shape {shape}")
            logical = stream.read(16)
            if len(logical) != 16:
                raise RuntimeError(f"{path}: truncated MeshBlock logical record")
            level = struct.unpack("=4i", logical)[3]
            levels[level] = levels.get(level, 0) + 1
            data_bytes = (
                len(header.variables) * math.prod(shape) * header.variable_size
            )
            stream.seek(6 * header.location_size + data_bytes, os.SEEK_CUR)
        if stream.tell() != file_size:
            raise RuntimeError(f"{path}: MeshBlock scan did not reach end of file")
    return levels


def subcycling_work_model(level_counts: dict[int, int]) -> dict[str, object]:
    if not level_counts:
        raise RuntimeError("cannot model subcycling work for an empty hierarchy")
    finest_level = max(level_counts)
    subcycled = sum(count * 2**level for level, count in level_counts.items())
    global_small_step = sum(level_counts.values()) * 2**finest_level
    return {
        "finest_physical_level": finest_level,
        "subcycled_meshblock_updates_per_root_step": subcycled,
        "global_finest_dt_meshblock_updates_per_root_step": global_small_step,
        "meshblock_update_reduction_fraction": 1.0 - subcycled / global_small_step,
        "global_to_subcycled_update_ratio": global_small_step / subcycled,
        "scope": "MeshBlock-update count only; excludes communication and output I/O",
    }


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


def configured_restart_dt(header) -> float | None:
    for section in header.parameters.values():
        if section.get("file_type") == "rst":
            try:
                value = float(section["dt"])
            except (KeyError, ValueError):
                return None
            return value if value > 0.0 and math.isfinite(value) else None
    return None


def storage_projection(
    streams: dict[str, list[dict[str, object]]],
    target_time: float,
    restart_sizes: list[int],
    restart_dt: float | None,
    segment_span: float | None,
    root_dt: float,
    root_step_seconds: float | None,
    drain_mib_s: float,
) -> dict[str, object]:
    projection: dict[str, object] = {"target_time_M": target_time, "streams": {}}
    binary_total = 0.0
    bytes_per_simulation_M = 0.0
    latest_time = max(
        float(frame["time_M"])
        for frames in streams.values()
        for frame in frames
    )
    segment_binary_bytes = 0.0
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
        configured_cycles = [
            int(frame["configured_root_dcycle"])
            for frame in ordered
            if frame.get("configured_root_dcycle") is not None
        ]
        if not configured and configured_cycles:
            configured = [float(np.median(configured_cycles)) * root_dt]
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
            binary_total += remaining_bytes
            median_size = float(np.median(sizes))
            bytes_per_simulation_M += median_size / cadence
            if segment_span is not None:
                segment_binary_bytes += math.ceil(segment_span / cadence) * median_size
        projection["streams"][name] = stream_projection
    projection["remaining_binary_bytes"] = binary_total
    projection["remaining_binary_TiB"] = binary_total / 2**40

    restart_projection: dict[str, object]
    restart_archive_bytes = 0.0
    remote_restart_working_set = 0.0
    if restart_sizes and (restart_dt is not None or segment_span is not None):
        median_restart = float(np.median(restart_sizes))
        effective_dt = (
            min(value for value in (restart_dt, segment_span) if value is not None)
        )
        remaining_restarts = max(
            0, math.ceil((target_time - latest_time) / effective_dt)
        )
        restart_archive_bytes = remaining_restarts * median_restart
        remote_restart_working_set = 2.0 * median_restart
        bytes_per_simulation_M += median_restart / effective_dt
        restart_projection = {
            "configured_dt_M": restart_dt,
            "forced_segment_span_M": segment_span,
            "effective_archive_cadence_M": effective_dt,
            "median_restart_bytes": median_restart,
            "remaining_restarts": remaining_restarts,
            "remaining_archive_bytes": restart_archive_bytes,
            "remote_retention_generations": 2,
            "remote_retained_restart_bytes": remote_restart_working_set,
        }
    else:
        restart_projection = {
            "status": "needs_a_verified_restart_and_cadence",
            "configured_dt_M": restart_dt,
            "forced_segment_span_M": segment_span,
        }
    projection["restarts"] = restart_projection
    archive_total = binary_total + restart_archive_bytes
    projection["remaining_archive_bytes"] = archive_total
    projection["remaining_archive_TiB"] = archive_total / 2**40
    if segment_span is not None:
        projection["remote_segment_working_set_bytes"] = (
            segment_binary_bytes + remote_restart_working_set
        )
    if root_step_seconds is not None:
        simulation_M_per_second = root_dt / root_step_seconds
        generation_bytes_s = bytes_per_simulation_M * simulation_M_per_second
        drain_bytes_s = drain_mib_s * 2**20
        projection["transfer_budget"] = {
            "root_dt_M": root_dt,
            "root_step_seconds": root_step_seconds,
            "average_generation_MiB_s": generation_bytes_s / 2**20,
            "assumed_sustained_drain_MiB_s": drain_mib_s,
            "drain_to_generation_ratio": (
                drain_bytes_s / generation_bytes_s
                if generation_bytes_s > 0.0
                else None
            ),
            "average_drain_has_headroom": drain_bytes_s > generation_bytes_s,
        }
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


def summarize_event_logs(paths: list[Path]) -> dict[str, object]:
    """Merge overlapping restart event logs and summarize interval counters."""

    base_columns = (
        "cycle",
        "eos_dfloor",
        "eos_efloor",
        "eos_tfloor",
        "eos_vceil",
        "eos_fail",
        "c2p_it",
        "fofc",
    )
    appended_columns = (
        "cons_adjust",
        "mag_adjust",
        "c2p_calls",
        "fofc_tests",
    )
    all_columns = base_columns + appended_columns
    supported_schemas = {frozenset(base_columns), frozenset(all_columns)}
    rows: dict[int, dict[str, int]] = {}
    for path in paths:
        header: tuple[str, ...] | None = None
        with path.open("r", encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, start=1):
                fields = line.split()
                if not fields:
                    continue
                if fields[0] == "#" and len(fields) > 1 and fields[1] == "cycle":
                    header = tuple(fields[1:])
                    if (
                        len(header) != len(set(header))
                        or frozenset(header) not in supported_schemas
                    ):
                        raise RuntimeError(
                            f"{path}:{line_number}: unsupported event-log header"
                        )
                    continue
                if fields[0].startswith("#"):
                    continue
                if header is None or len(fields) != len(header):
                    raise RuntimeError(f"{path}:{line_number}: malformed event log")
                try:
                    row = {
                        name: int(value) for name, value in zip(header, fields)
                    }
                except ValueError as error:
                    raise RuntimeError(
                        f"{path}:{line_number}: malformed event log"
                    ) from error
                for name in appended_columns:
                    row.setdefault(name, 0)
                rows[row["cycle"]] = row
        if header is None:
            raise RuntimeError(f"{path}: event-log header is missing or unsupported")

    ordered = [rows[cycle] for cycle in sorted(rows)]
    corrective_columns = tuple(
        name
        for name in all_columns[1:]
        if name not in ("c2p_it", "c2p_calls", "fofc_tests")
    )
    records_with_events = sum(
        any(row[name] != 0 for name in corrective_columns) for row in ordered
    )
    result: dict[str, object] = {
        "verified_sources": [str(path) for path in paths],
        "records": len(ordered),
        "records_with_events": records_with_events,
        "records_without_corrective_events": len(ordered) - records_with_events,
        "totals": {
            name: sum(row[name] for row in ordered)
            for name in all_columns[1:]
            if name != "c2p_it"
        },
        "c2p_iteration_max": max((row["c2p_it"] for row in ordered), default=0),
    }
    if ordered:
        result["cycle_min"] = ordered[0]["cycle"]
        result["cycle_max"] = ordered[-1]["cycle"]
    return result


def write_history_products(
    output_dir: Path,
    columns: list[str],
    data: dict[str, np.ndarray],
    stem: str = "merged-history",
) -> tuple[Path, Path]:
    csv_path = output_dir / f"{stem}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(zip(*(data[column] for column in columns)))

    hst_path = output_dir / f"{stem}.hst"
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
        "mhd_w_bcc": [
            "dens",
            "press",
            "temperature",
            "bmag",
            "beta_inv_proxy",
            "level",
        ],
        "mhd_gr_diagnostics": [
            "gr_bsq",
            "gr_lorentz",
            "gr_sigma",
            "gr_beta_inv",
            "level",
        ],
        "mhd_divb": ["divb", "level"],
        "bbh_local_w": [
            "dens",
            "press",
            "temperature",
            "bmag",
            "beta_inv_proxy",
            "level",
        ],
        "bbh_local_gr": [
            "gr_bsq",
            "gr_lorentz",
            "gr_sigma",
            "gr_beta_inv",
            "level",
        ],
    }
    for stream_name, panel_names in choices.items():
        for frame in selected_frames(streams.get(stream_name, []), every):
            input_path = Path(str(frame["path"]))
            local_region = frame.get("region") is not None
            slice_data = read_slice(
                input_path,
                panel_names,
                plane,
                location,
                None if local_region else extent,
            )
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
                None if local_region else extent,
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


def write_markdown_report(report: dict[str, object], path: Path) -> None:
    lines = [
        "# Verified BBH GRMHD campaign analysis",
        "",
        "## Integrity",
        "",
        "Only files chained through an immutable ready manifest and local ACK are "
        "included. Science-file size and SHA256 checks are recorded in "
        "`campaign-analysis.json`.",
        "",
        "## Output streams",
        "",
        "| stream | frames | time range (M) | median gap (M) | latest blocks | "
        "subcycling update reduction |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    streams = report["streams"]
    cadence = report["cadence"]
    for name in sorted(streams):
        frames = streams[name]
        summary = cadence[name]
        latest = max(frames, key=lambda frame: float(frame["time_M"]))
        work = latest["subcycling_work_model"]
        lines.append(
            f"| {name} | {summary['frames']} | "
            f"{summary['time_min_M']:.6g}–{summary['time_max_M']:.6g} | "
            f"{summary.get('gap_median_M', 'n/a')} | {latest['meshblocks']} | "
            f"{100.0 * work['meshblock_update_reduction_fraction']:.2f}% |"
        )

    storage = report["storage_projection"]
    lines.extend(
        [
            "",
            "## Storage and drain budget",
            "",
            f"- Target time: {storage['target_time_M']:.6g}M",
            f"- Remaining binary output: {storage['remaining_binary_TiB']:.3f} TiB",
            f"- Remaining total archive: {storage['remaining_archive_TiB']:.3f} TiB",
            f"- Remote segment working set: "
            f"{storage.get('remote_segment_working_set_bytes', 0) / 2**30:.1f} GiB",
        ]
    )
    transfer = storage.get("transfer_budget")
    if transfer is not None:
        lines.extend(
            [
                f"- Average generation: {transfer['average_generation_MiB_s']:.3f} "
                "MiB/s",
                f"- Assumed sustained drain: "
                f"{transfer['assumed_sustained_drain_MiB_s']:.3f} MiB/s",
                f"- Drain/generation headroom: "
                f"{transfer['drain_to_generation_ratio']:.2f}x",
            ]
        )

    history = report.get("history")
    lines.extend(["", "## History diagnostics", ""])
    if not history:
        lines.append("No verified history file was available.")
    else:
        for physics, summary in sorted(history.items()):
            lines.append(f"### {physics} history")
            lines.append("")
            for name, value in sorted(summary["diagnostics"].items()):
                lines.append(f"- `{name}`: {value:.8g}")
            lines.append("")
            for note in summary.get("interpretation_notes", []):
                lines.append(f"> {note}")
                lines.append("")

    events = report.get("event_counters")
    lines.extend(["## Numerical event counters", ""])
    if events is None:
        lines.append("No verified Athena event log was available.")
    else:
        lines.append(f"- Logged root endpoints: {events.get('records', 0)}")
        lines.append(
            "- Endpoints containing corrective/failure events: "
            f"{events['records_with_events']}"
        )
        lines.append(f"- Maximum C2P iterations: {events['c2p_iteration_max']}")
        for name, value in sorted(events["totals"].items()):
            lines.append(f"- `{name}` total: {value}")
        lines.append("")

    readiness = report["publication_readiness"]
    lines.extend(["## Publication gates still open", ""])
    for item in readiness["not_available_from_current_dump"]:
        lines.append(f"- {item}")
    for item in readiness["campaign_systematics_still_required"]:
        lines.append(f"- {item}")
    lines.extend(
        [
            "",
            "Coordinate-component magnetic and velocity panels remain explicitly "
            "labelled as proxies until synchronized metric data are present.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


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
    parser.add_argument(
        "--segment-span",
        type=float,
        default=100.0,
        help="simulation time per cloud segment; 0 disables forced-restart budgeting",
    )
    parser.add_argument("--root-dt", type=float, default=4.8)
    parser.add_argument(
        "--root-step-seconds",
        type=float,
        help="measured wall seconds per root step, enabling network-rate budgeting",
    )
    parser.add_argument("--drain-mib-s", type=float, default=8.0)
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
    if args.segment_span < 0.0 or args.root_dt <= 0.0 or args.drain_mib_s <= 0.0:
        raise SystemExit(
            "segment span must be non-negative; root dt and drain rate must be positive"
        )
    if args.root_step_seconds is not None and args.root_step_seconds <= 0.0:
        raise SystemExit("--root-step-seconds must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_records: list[VerifiedFile] = []
    for segment in args.segments:
        all_records.extend(load_verified_files(segment))
    science_records = [
        record
        for record in all_records
        if record.path.suffix in (".bin", ".hst", ".log")
    ]
    if not science_records:
        raise RuntimeError("ACKs contain no binary, history, or event-log science output")
    for record in science_records:
        verify_file(record, not args.no_sha256)
        print(f"verified {record.relative} ({record.size} bytes)", file=sys.stderr)

    streams: dict[str, list[dict[str, object]]] = {
        "mhd_w_bcc": [],
        "mhd_gr_diagnostics": [],
        "mhd_divb": [],
        "bbh_local_w": [],
        "bbh_local_gr": [],
        "other": [],
    }
    for record in science_records:
        if record.path.suffix != ".bin":
            continue
        header = read_binary_header(record.path)
        stream_name = classify_binary(record.path.name, header.variables)
        level_counts = meshblock_level_counts(record.path, header)
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
                "meshblocks": sum(level_counts.values()),
                "meshblocks_by_physical_level": {
                    str(level): count for level, count in sorted(level_counts.items())
                },
                "subcycling_work_model": subcycling_work_model(level_counts),
                "configured_dt_M": configured_output_dt(header, stream_name),
                "configured_root_dcycle": configured_output_dcycle(
                    header, stream_name
                ),
                "region": configured_output_region(header, stream_name),
            }
        )
    for frames in streams.values():
        frames.sort(key=lambda frame: (float(frame["time_M"]), str(frame["path"])))
    streams = {name: frames for name, frames in streams.items() if frames}

    history_records = [record for record in science_records if record.path.suffix == ".hst"]
    event_records = [
        record
        for record in science_records
        if record.path.suffix == ".log" and is_event_log(record.path)
    ]
    restart_records = [record for record in all_records if record.path.suffix == ".rst"]
    for record in restart_records:
        verify_file(record, False)
    # Preserve command-line segment order so later restart segments replace overlap.
    history_groups: dict[str, list[Path]] = {}
    for record in history_records:
        history_groups.setdefault(classify_history(record.path.name), []).append(record.path)
    history_summary: dict[str, dict[str, object]] = {}
    for physics, history_paths in sorted(history_groups.items()):
        columns, history_data = merge_histories(history_paths)
        stem = "merged-history" if physics == "mhd" else f"merged-{physics}-history"
        csv_path, hst_path = write_history_products(
            args.output_dir, columns, history_data, stem=stem
        )
        summary = create_history_plot(
            hst_path, args.output_dir / f"{stem}.png", args.dpi
        )
        summary["csv"] = str(csv_path)
        summary["verified_sources"] = [str(path) for path in history_paths]
        history_summary[physics] = summary

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
    native_gr_diagnostics = bool(
        streams.get("mhd_gr_diagnostics") or streams.get("bbh_local_gr")
    )
    native_bbh_history = "user" in history_summary
    exact_outputs = [
        "density",
        "pressure",
        "temperature",
        "stored bcc components (densitized in DynGRMHD)",
        "primitive velocity",
        "divB",
    ]
    unavailable = [
        "MRI quality factors Q_theta and Q_phi",
        "moving-horizon mass and magnetic fluxes",
        "covariant Bernoulli parameter",
    ]
    if native_gr_diagnostics:
        exact_outputs.extend(
            [
                "comoving magnetic invariant b^2",
                "Lorentz factor W",
                "magnetization b^2/rho",
                "inverse magnetic beta b^2/(2p)",
            ]
        )
    else:
        unavailable.insert(0, "covariant magnetization, plasma beta, and Lorentz factor")
    if native_bbh_history:
        exact_outputs.extend(
            [
                "outside-excision global baryon mass and GRMHD proper-volume integrals",
                "BBH positions, separation, orbital angular frequency, and term masses",
                "outside-excision rho and magnetization maxima",
            ]
        )

    report = {
        "schema": 2,
        "classification": "athenak-bbh-grmhd-verified-campaign-analysis",
        "segments": [str(path.expanduser().resolve()) for path in args.segments],
        "input_policy": {
            "local_ack_required": True,
            "science_file_size_verified": True,
            "science_file_sha256_verified": not args.no_sha256,
            "restart_size_verified_for_budget": True,
            "restart_sha256_trusted_from_transfer_ack": True,
            "active_or_partial_files_accepted": False,
        },
        "streams": streams,
        "cadence": {
            name: cadence_summary(frames) for name, frames in streams.items()
        },
        "storage_projection": storage_projection(
            streams,
            args.target_time,
            [record.size for record in restart_records],
            configured_restart_dt(
                max(
                    (
                        read_binary_header(Path(str(frame["path"])))
                        for frames in streams.values()
                        for frame in frames
                    ),
                    key=lambda header: header.time,
                )
            ),
            args.segment_span if args.segment_span > 0.0 else None,
            args.root_dt,
            args.root_step_seconds,
            args.drain_mib_s,
        ),
        "history": history_summary,
        "event_counters": (
            summarize_event_logs([record.path for record in event_records])
            if event_records
            else None
        ),
        "timeline": str(timeline_path),
        "rendered_frames": rendered,
        "publication_readiness": {
            "native_gr_diagnostics_available": native_gr_diagnostics,
            "native_bbh_history_available": native_bbh_history,
            "exact_outputs": exact_outputs,
            "proxy_warning": PROXY_WARNING,
            "not_available_from_current_dump": unavailable,
            "campaign_systematics_still_required": [
                "bulk-disk resolution sequence beyond the fixed physical-L4 floor",
                "horizon-following magnetic leakage/forcing validation",
                "half-CFL and metric finite-difference-step sensitivity",
                "floor, excision, and outer-boundary sensitivity",
            ],
        },
    }
    report_path = args.output_dir / "campaign-analysis.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    markdown_path = args.output_dir / "campaign-analysis.md"
    write_markdown_report(report, markdown_path)
    print(markdown_path)
    print(report_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.stderr.close()
        raise SystemExit(1)
