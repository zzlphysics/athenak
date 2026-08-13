#!/usr/bin/env python3
"""Qualify a closed, scheduled prefix of an interrupted AthenaK segment.

This is an offline, fail-closed recovery tool.  It never resumes a process, never
modifies the original state/evidence trees, and never weakens the complete-segment
checker.  A recovery pass is a distinct qualification mode that may only end at an
original-plan ``scheduled`` restart whose Execute cadence state is replayed exactly.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import socket
import stat
import struct
import sys
import tempfile
import time
from typing import Any, Iterable


sys.dont_write_bytecode = True
SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIRECTORY = SCRIPT_PATH.parent
if str(SCRIPT_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

import audit_athenak_restart as RESTART_AUDITOR
import check_athenak_segment as CHECKER
from audit_athenak_restart import audit_restart
from output_integrity import (
    HASH_CHUNK_BYTES,
    _assert_closed,
    stable_sha256,
    strict_json_loads,
)
from read_athenak_restart_metadata import read_restart_metadata


SCHEMA = 1
QUALIFICATION_MODE = "scheduled_prefix_recovery_v1"
RECOVERY_POLICY = {
    "kind": QUALIFICATION_MODE,
    "candidate_kind": "original_plan_scheduled_restart_only",
    "selection": "highest_complete_candidate",
    "later_complete_scientific_failure": "fatal_no_fallback",
    "cadence_replay": "execute_only_binary64_v1",
    "original_trees": "read_only_no_delete",
    "derived_text": "exact_original_byte_prefix",
    "unknown_state_nodes": "fatal",
    "lifecycle": "nonzero_completion_or_same_boot_closed_processes",
}
SAME_BOOT_QUIESCENCE_TIMEOUT_SECONDS = 10.0
SAME_BOOT_QUIESCENCE_POLL_SECONDS = 0.1
BOOT_ID_PATH = Path("/proc/sys/kernel/random/boot_id")
PROC_STAT_PATH = Path("/proc/stat")


class RecoveryFailure(RuntimeError):
    """A permanent failure to prove a safe scheduled-prefix recovery."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RecoveryFailure(message)


@dataclass
class BoundDirectory:
    """Directory identity pinned by O_DIRECTORY/O_NOFOLLOW for this recovery."""

    path: Path
    fd: int
    device: int
    inode: int
    owner_uid: int
    mode: int

    @classmethod
    def open(cls, path: Path, label: str, *, required_mode: int | None = None,
             require_empty: bool = False) -> "BoundDirectory":
        absolute = path.expanduser().absolute()
        flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_CLOEXEC", 0) | \
            getattr(os, "O_NOFOLLOW", 0)
        try:
            descriptor = os.open(absolute, flags)
            info = os.fstat(descriptor)
            path_info = absolute.lstat()
        except OSError as exc:
            raise RecoveryFailure(f"cannot bind {label}: {exc}") from exc
        signature = (info.st_dev, info.st_ino)
        if (not stat.S_ISDIR(info.st_mode) or stat.S_ISLNK(path_info.st_mode) or
                signature != (path_info.st_dev, path_info.st_ino) or
                info.st_uid != os.geteuid() or
                (required_mode is not None and
                 stat.S_IMODE(info.st_mode) != required_mode)):
            os.close(descriptor)
            raise RecoveryFailure(
                f"{label} is not a stable owner-owned directory with required mode")
        if require_empty and os.listdir(descriptor):
            os.close(descriptor)
            raise RecoveryFailure(f"{label} must be empty")
        return cls(absolute.resolve(strict=True), descriptor, info.st_dev,
                   info.st_ino, info.st_uid, stat.S_IMODE(info.st_mode))

    def audit(self, label: str, *, require_empty: bool = False) -> dict[str, Any]:
        try:
            info = os.fstat(self.fd)
            path_info = self.path.lstat()
        except OSError as exc:
            raise RecoveryFailure(f"cannot re-audit {label}: {exc}") from exc
        _require(stat.S_ISDIR(info.st_mode) and
                 (info.st_dev, info.st_ino, info.st_uid,
                  stat.S_IMODE(info.st_mode)) ==
                 (self.device, self.inode, self.owner_uid, self.mode) and
                 (path_info.st_dev, path_info.st_ino, path_info.st_uid,
                  stat.S_IMODE(path_info.st_mode)) ==
                 (self.device, self.inode, self.owner_uid, self.mode) and
                 not stat.S_ISLNK(path_info.st_mode),
                 f"{label} path or descriptor identity changed")
        if require_empty:
            _require(not os.listdir(self.fd), f"{label} is no longer empty")
        return {
            "path": str(self.path), "device": self.device, "inode": self.inode,
            "owner_uid": self.owner_uid, "mode": f"{self.mode:04o}",
            "descriptor_bound": True,
        }

    def close(self) -> None:
        if self.fd >= 0:
            os.close(self.fd)
            self.fd = -1


def _canonical_json_sha256(value: Any) -> str:
    raw = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _binding(path: Path) -> dict[str, Any]:
    try:
        audited = stable_sha256(path)
    except (OSError, RuntimeError, ValueError) as exc:
        raise RecoveryFailure(f"cannot bind closed file {path}: {exc}") from exc
    return {
        "path": audited["path"],
        **{name: audited[name] for name in
           ("device", "inode", "size", "mtime_ns", "ctime_ns")},
        "sha256": audited["sha256"],
        "closure_check": audited["closure_check"],
    }


def _binding_from_open_fd(path: Path, descriptor: int) -> dict[str, Any]:
    """Hash the exact inode already pinned by a no-follow directory walk."""

    before = os.fstat(descriptor)
    _require(stat.S_ISREG(before.st_mode), f"pinned artifact is not regular: {path}")
    signature = {
        "device": before.st_dev, "inode": before.st_ino,
        "size": before.st_size, "mtime_ns": before.st_mtime_ns,
        "ctime_ns": before.st_ctime_ns,
    }
    exempt = {(os.getpid(), descriptor)}
    try:
        _assert_closed(path, signature, exempt)
        os.lseek(descriptor, 0, os.SEEK_SET)
        digest = hashlib.sha256()
        while True:
            piece = os.read(descriptor, HASH_CHUNK_BYTES)
            if not piece:
                break
            digest.update(piece)
        after = os.fstat(descriptor)
        _require((after.st_dev, after.st_ino, after.st_size,
                  after.st_mtime_ns, after.st_ctime_ns) ==
                 (before.st_dev, before.st_ino, before.st_size,
                  before.st_mtime_ns, before.st_ctime_ns),
                 f"pinned artifact changed while hashing: {path}")
        current = path.lstat()
        _require((current.st_dev, current.st_ino, current.st_size,
                  current.st_mtime_ns, current.st_ctime_ns) ==
                 (before.st_dev, before.st_ino, before.st_size,
                  before.st_mtime_ns, before.st_ctime_ns),
                 f"pinned artifact path changed while hashing: {path}")
        _assert_closed(path, signature, exempt)
    except (OSError, RuntimeError, ValueError) as exc:
        raise RecoveryFailure(f"cannot bind pinned file {path}: {exc}") from exc
    return {
        "path": str(path.resolve(strict=True)), **signature,
        "sha256": digest.hexdigest(), "closure_check": "linux_proc_fd",
    }


def _same_binding(left: Any, right: Any) -> bool:
    return isinstance(left, dict) and isinstance(right, dict) and all(
        left.get(name) == right.get(name) for name in
        ("path", "device", "inode", "size", "mtime_ns", "ctime_ns", "sha256")
    )


def _require_settled(binding: dict[str, Any], minimum_age_seconds: float,
                     label: str) -> None:
    age = time.time() - float(binding["mtime_ns"]) / 1.0e9
    _require(math.isfinite(age) and age >= minimum_age_seconds,
             f"{label} is not settled for {minimum_age_seconds:g} seconds")


def _immutable_json(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    absolute = path.expanduser().absolute()
    try:
        info = absolute.lstat()
    except OSError as exc:
        raise RecoveryFailure(f"cannot stat immutable JSON {absolute}: {exc}") from exc
    _require(stat.S_ISREG(info.st_mode) and not stat.S_ISLNK(info.st_mode) and
             not (stat.S_IMODE(info.st_mode) & 0o222),
             f"immutable JSON must be a read-only regular file: {absolute}")
    binding = _binding(absolute)
    try:
        raw = absolute.read_bytes()
        value = strict_json_loads(raw, str(absolute))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise RecoveryFailure(f"invalid immutable JSON {absolute}: {exc}") from exc
    _require(isinstance(value, dict), f"JSON root is not an object: {absolute}")
    _require(hashlib.sha256(raw).hexdigest() == binding["sha256"],
             f"JSON changed while parsing: {absolute}")
    return value, binding


def _planned_file(plan: dict[str, Any], name: str, actual: Path) -> dict[str, Any]:
    record = plan.get("tools", {}).get(name)
    _require(isinstance(record, dict), f"plan does not bind tools.{name}")
    canonical = actual.resolve(strict=True)
    _require(record.get("path") == str(canonical),
             f"current {name} path differs from original plan")
    current = _binding(canonical)
    _require(current["size"] == record.get("size") and
             current["sha256"] == record.get("sha256"),
             f"current {name} bytes differ from original plan")
    return current


def _private_empty_directory(path: Path, label: str) -> Path:
    absolute = path.expanduser().absolute()
    try:
        info = absolute.lstat()
    except OSError as exc:
        raise RecoveryFailure(f"cannot stat {label}: {exc}") from exc
    _require(stat.S_ISDIR(info.st_mode) and not stat.S_ISLNK(info.st_mode),
             f"{label} must be a real directory")
    resolved = absolute.resolve(strict=True)
    _require(absolute == resolved and info.st_uid == os.geteuid() and
             stat.S_IMODE(info.st_mode) == 0o700,
             f"{label} must be canonical, owner-owned, and mode 0700")
    _require(not any(resolved.iterdir()), f"{label} must be empty")
    return resolved


def _separate_directories(paths: Iterable[Path]) -> None:
    values = list(paths)
    for index, left in enumerate(values):
        for right in values[index + 1:]:
            _require(left != right, "original and recovery directories must differ")
            for child, parent in ((left, right), (right, left)):
                try:
                    child.relative_to(parent)
                except ValueError:
                    continue
                raise RecoveryFailure(
                    "original and recovery directories may not contain one another")


def _float32(value: float) -> float:
    try:
        return struct.unpack("=f", struct.pack("=f", value))[0]
    except (OverflowError, struct.error) as exc:
        raise RecoveryFailure(f"output time is not float32-representable: {value}") \
            from exc


def replay_execute_prefix(plan: dict[str, Any], final_cycle: int) \
        -> tuple[dict[str, list[dict[str, Any]]], dict[str, dict[str, Any]], float]:
    """Replay only ``Outputs::Execute`` through ``final_cycle``."""

    expected = plan["expected"]
    source_cycle = int(expected["source_cycle"])
    _require(source_cycle < final_cycle <= int(expected["final_cycle"]),
             "recovery cycle is outside the original segment")
    states: dict[str, dict[str, Any]] = {}
    for output in plan["outputs"]:
        block = output["block"]
        parameters = output["parameters"]
        try:
            states[block] = {
                "file_number": int(parameters.get("file_number", "0")),
                "last_time": float(parameters.get("last_time", "-1")),
                "last_write_cycle": int(parameters.get("last_write_cycle", "-1")),
                "dt": (float(parameters["dt"])
                       if output["cadence_mode"] == "dt" else 0.0),
                "dcycle": (int(parameters["dcycle"])
                           if output["cadence_mode"] == "dcycle" else 0),
                "writes": [],
            }
        except (KeyError, TypeError, ValueError) as exc:
            raise RecoveryFailure(f"invalid cadence state for {block}") from exc

    def due(state: dict[str, Any], cycle: int, output_time: float) -> bool:
        return ((state["dt"] > 0.0 and
                 _float32(output_time) >=
                 _float32(state["last_time"] + state["dt"]) and
                 _float32(output_time) < _float32(expected["tlim"])) or
                (state["dcycle"] > 0 and cycle % state["dcycle"] == 0))

    current_time = float(expected["source_time"])
    root_dt = float(expected["root_dt"])
    for cycle in range(source_cycle + 1, final_cycle + 1):
        current_time += root_dt
        for output in plan["outputs"]:
            state = states[output["block"]]
            if not due(state, cycle, current_time):
                continue
            write: dict[str, Any] = {
                "cycle": cycle, "time": current_time, "kind": "scheduled",
            }
            if output["numbered"]:
                write["file_number"] = state["file_number"]
                state["file_number"] += 1
            state["writes"].append(write)
            state["last_time"] = (current_time if state["last_time"] < 0.0
                                  else state["last_time"] + state["dt"])
            state["last_write_cycle"] = cycle
    return ({block: value["writes"] for block, value in states.items()},
            states, current_time)


def _restart_output(plan: dict[str, Any]) -> dict[str, Any]:
    rows = [row for row in plan["outputs"] if row.get("file_type") == "rst"]
    _require(len(rows) == 1, "plan must contain one restart stream")
    return rows[0]


def _relative_for(output: dict[str, Any], write: dict[str, Any]) -> str:
    template = output.get("relative_path_template")
    _require(isinstance(template, str), f"{output.get('block')} has no path template")
    return template.format(file_number=write["file_number"])


def _safe_original_path(root: Path, relative: str, *, must_exist: bool = True) -> Path:
    rel = Path(relative)
    _require(not rel.is_absolute() and relative == rel.as_posix() and
             all(part not in ("", ".", "..") for part in rel.parts),
             f"noncanonical planned relative path: {relative!r}")
    candidate = root / rel
    current = candidate.absolute()
    while current != root.parent:
        _require(not current.is_symlink(), f"original state path traverses symlink: {current}")
        if current == root:
            break
        current = current.parent
    if must_exist:
        resolved = candidate.resolve(strict=True)
        try:
            resolved.relative_to(root)
        except ValueError as exc:
            raise RecoveryFailure(f"original state path escapes tree: {candidate}") from exc
        return resolved
    return candidate


def _inventory_original_tree(state_dir: Path, known: set[str],
                             bound: BoundDirectory | None = None) \
        -> dict[str, Path]:
    actual: dict[str, Path] = {}
    known_directories = {"."}
    for relative in known:
        parent = Path(relative).parent
        while str(parent) != ".":
            known_directories.add(parent.as_posix())
            parent = parent.parent
    if bound is not None:
        bound.audit("original state directory")
    pending: list[tuple[Path, int | None]] = [
        (state_dir, os.dup(bound.fd) if bound is not None else None)]
    while pending:
        directory, directory_fd = pending.pop()
        try:
            entries = list(os.scandir(directory_fd if directory_fd is not None
                                      else directory))
        except OSError as exc:
            if directory_fd is not None:
                os.close(directory_fd)
            raise RecoveryFailure(f"cannot inventory original state: {exc}") from exc
        for entry in entries:
            node = entry.stat(follow_symlinks=False)
            path = directory / entry.name
            _require(not stat.S_ISLNK(node.st_mode),
                     f"original state contains symlink: {path}")
            if stat.S_ISDIR(node.st_mode):
                relative_directory = path.relative_to(state_dir).as_posix()
                _require(relative_directory in known_directories,
                         "original state contains unknown extra directory: "
                         f"{relative_directory}")
                child_fd = None
                if directory_fd is not None:
                    child_fd = os.open(
                        entry.name, os.O_RDONLY | os.O_DIRECTORY |
                        getattr(os, "O_CLOEXEC", 0) |
                        getattr(os, "O_NOFOLLOW", 0), dir_fd=directory_fd)
                    child_info = os.fstat(child_fd)
                    _require((child_info.st_dev, child_info.st_ino) ==
                             (node.st_dev, node.st_ino),
                             f"original directory changed while opening: {path}")
                pending.append((path, child_fd))
                continue
            _require(stat.S_ISREG(node.st_mode),
                     f"original state contains non-regular node: {path}")
            relative = path.relative_to(state_dir).as_posix()
            _require(relative in known,
                     f"original state contains unknown extra artifact: {relative}")
            if directory_fd is not None:
                file_fd = os.open(
                    entry.name, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                    getattr(os, "O_NOFOLLOW", 0), dir_fd=directory_fd)
                try:
                    file_info = os.fstat(file_fd)
                    _require(stat.S_ISREG(file_info.st_mode) and
                             (file_info.st_dev, file_info.st_ino) ==
                             (node.st_dev, node.st_ino),
                             f"original file changed while opening: {path}")
                    live = _binding_from_open_fd(path, file_fd)
                    _require((live["device"], live["inode"]) ==
                             (file_info.st_dev, file_info.st_ino),
                             f"original file path changed while binding: {path}")
                finally:
                    os.close(file_fd)
            actual[relative] = path.resolve(strict=True)
        if directory_fd is not None:
            os.close(directory_fd)
    if bound is not None:
        bound.audit("original state directory")
    return actual


def _known_paths(plan: dict[str, Any]) -> tuple[set[str], dict[str, dict[str, Any]]]:
    known = set(plan["required_files"])
    numbered: dict[str, dict[str, Any]] = {}
    for output in plan["outputs"]:
        if not output["numbered"]:
            continue
        for write in output["expected_writes"]:
            relative = _relative_for(output, write)
            _require(relative not in numbered,
                     f"plan maps two logical writes to {relative}")
            numbered[relative] = {"output": output, "write": write}
            known.add(relative)
    return known, numbered


def _expected_restart_size(path: Path) -> tuple[Any, int]:
    metadata = read_restart_metadata(path)
    layout = RESTART_AUDITOR._derive_layout(metadata)
    expected_size = metadata.metadata_end + struct.calcsize("Q") + \
        metadata.nmb_total * layout.data_size
    return metadata, expected_size


def classify_restart(path: Path, minimum_complete_header_bytes: int = 0) \
        -> dict[str, Any]:
    """Distinguish an obvious truncated write from a complete invalid restart."""
    try:
        return CHECKER.classify_recovery_restart(path)
    except Exception as exc:
        raise RecoveryFailure(f"cannot classify restart {path}: {exc}") from exc


def select_highest_complete_candidate(
        candidates: list[dict[str, Any]]) -> dict[str, Any]:
    """Select once; callers must never catch qualification failure and fall back."""

    ordered = sorted(candidates, key=lambda row: row["expected_write"]["cycle"],
                     reverse=True)
    invalid = next((row for row in ordered
                    if row["classification"] == "complete_invalid"), None)
    complete = next((row for row in ordered
                     if row["classification"] == "complete"), None)
    if invalid is not None and (complete is None or
                                invalid["expected_write"]["cycle"] >
                                complete["expected_write"]["cycle"]):
        raise RecoveryFailure(
            "a later complete restart is invalid; fallback to an older prefix is "
            f"forbidden: {invalid['relative_path']}: {invalid.get('failure')}")
    _require(complete is not None, "no complete scheduled restart is recoverable")
    return complete


def _assert_execute_endpoint(plan: dict[str, Any], launch: dict[str, Any],
                             audit: dict[str, Any], final_cycle: int,
                             final_time: float,
                             states: dict[str, dict[str, Any]]) -> Any:
    metadata = read_restart_metadata(Path(audit["path"]))
    _require(metadata.cycle == final_cycle and metadata.time == final_time,
             "scheduled restart metadata does not match its planned cycle/time")
    for output in plan["outputs"]:
        block = output["block"]
        endpoint = metadata.parameters.get(block)
        _require(isinstance(endpoint, dict), f"scheduled restart lacks <{block}>")
        try:
            observed = {
                "file_number": int(endpoint.get("file_number", "0")),
                "last_time": float(endpoint["last_time"]),
                "last_write_cycle": int(endpoint["last_write_cycle"]),
            }
        except (KeyError, ValueError) as exc:
            raise RecoveryFailure(f"scheduled restart has invalid {block} cadence") \
                from exc
        expected_state = {name: states[block][name] for name in
                          ("file_number", "last_time", "last_write_cycle")}
        _require(observed == expected_state,
                 f"{block} cadence differs from Execute-only replay: "
                 f"{observed} != {expected_state}")

    partial_expected = dict(plan["expected"])
    partial_expected.update({
        "final_cycle": final_cycle,
        "final_time": final_time,
        "root_steps": final_cycle - int(plan["expected"]["source_cycle"]),
    })
    CHECKER.audit_restart_contract(metadata, "recovered endpoint", plan,
                                   partial_expected)
    source = read_restart_metadata(Path(plan["inputs"]["source_restart"]["path"]))
    trajectory_tokens = [token for token in launch.get("athena_argv", [])
                         if isinstance(token, str) and
                         token.startswith("problem/trajectory_file=")]
    _require(len(trajectory_tokens) == 1,
             "launch record lacks exact runtime trajectory rebinding")
    runtime_trajectory = trajectory_tokens[0].split("=", 1)[1]
    source_trajectory = source.parameters.get("problem", {}).get("trajectory_file")
    exact_rebindings = {
            "output3/dt": {
                "source": [source.parameters["output3"]["dt"]],
                "endpoint": repr(plan["expected"]["root_dt"]),
            },
            "problem/trajectory_file": {
                "source": [source_trajectory], "endpoint": runtime_trajectory,
            },
        }
    transition = plan.get("capacity_transition", {})
    if transition.get("kind") == "increase_v1":
        exact_rebindings["mesh_refinement/max_nmb_per_rank"] = {
            "source": [str(transition["source_max_nmb_per_rank"])],
            "endpoint": str(transition["target_max_nmb_per_rank"]),
        }
    CHECKER.compare_parameters(
        source, metadata, plan["policy"]["mutable_parameters"],
        exact_rebindings)
    return metadata


def _exact_text_prefix(raw: bytes, data_rows: int, label: str) \
        -> tuple[bytes, bytes]:
    _require(data_rows > 0, f"{label} prefix must contain data")
    rows = 0
    offset = 0
    for line in raw.splitlines(keepends=True):
        offset += len(line)
        stripped = line.strip()
        if stripped and not stripped.startswith(b"#"):
            rows += 1
            if rows == data_rows:
                _require(line.endswith(b"\n"),
                         f"{label} target row is not a closed newline-terminated row")
                return raw[:offset], raw[offset:]
    raise RecoveryFailure(f"{label} is truncated before logical prefix row {data_rows}")


def _audit_text_bytes(raw: bytes, callback: Any) -> Any:
    descriptor, name = tempfile.mkstemp(prefix="athenak-prefix-audit-", suffix=".txt")
    path = Path(name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(path, 0o400)
        return callback(path)
    finally:
        try:
            path.unlink()
        except FileNotFoundError:
            pass


def _install_prefix_bytes(root: Path, relative: str, raw: bytes,
                          root_binding: BoundDirectory) -> dict[str, Any]:
    root_binding.audit("recovery state directory")
    destination = root / Path(relative)
    parts = Path(relative).parts
    _require(parts and all(part not in ("", ".", "..") for part in parts),
             "derived prefix has a noncanonical relative path")
    parent_fd = os.dup(root_binding.fd)
    try:
        for component in parts[:-1]:
            try:
                os.mkdir(component, mode=0o700, dir_fd=parent_fd)
                os.fsync(parent_fd)
            except FileExistsError:
                pass
            child_fd = os.open(
                component, os.O_RDONLY | os.O_DIRECTORY |
                getattr(os, "O_CLOEXEC", 0) |
                getattr(os, "O_NOFOLLOW", 0), dir_fd=parent_fd)
            child_info = os.fstat(child_fd)
            _require(stat.S_ISDIR(child_info.st_mode) and
                     child_info.st_uid == os.geteuid() and
                     stat.S_IMODE(child_info.st_mode) == 0o700,
                     "derived prefix parent is not an owner-only directory")
            os.close(parent_fd)
            parent_fd = child_fd
        try:
            os.stat(parts[-1], dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise RecoveryFailure(f"refusing to overwrite derived prefix: {destination}")
        temporary_name = f".{parts[-1]}.{os.getpid()}.{time.time_ns()}.tmp"
        descriptor: int | None = None
        try:
            descriptor = os.open(
                temporary_name, os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
                0o600, dir_fd=parent_fd)
            with os.fdopen(descriptor, "wb", closefd=True) as stream:
                descriptor = None
                stream.write(raw)
                stream.flush()
                os.fsync(stream.fileno())
                os.fchmod(stream.fileno(), 0o400)
            try:
                os.link(temporary_name, parts[-1], src_dir_fd=parent_fd,
                        dst_dir_fd=parent_fd, follow_symlinks=False)
            except FileExistsError as exc:
                raise RecoveryFailure(
                    f"refusing to replace derived prefix: {destination}") from exc
            os.unlink(temporary_name, dir_fd=parent_fd)
            os.fsync(parent_fd)
        except BaseException:
            if descriptor is not None:
                os.close(descriptor)
            try:
                os.unlink(temporary_name, dir_fd=parent_fd)
            except FileNotFoundError:
                pass
            raise
    finally:
        os.close(parent_fd)
    root_binding.audit("recovery state directory")
    return _binding(destination)


def _install_bound_json(root_binding: BoundDirectory, basename: str,
                        payload: Any) -> dict[str, Any]:
    """Install immutable JSON through a held directory descriptor."""

    _require(Path(basename).name == basename and basename not in ("", ".", ".."),
             "bound JSON destination must be one canonical basename")
    encoded = (json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) +
               "\n").encode("utf-8")
    return _install_prefix_bytes(
        root_binding.path, basename, encoded, root_binding)


def _parse_utc(value: Any) -> float:
    _require(isinstance(value, str) and value,
             "launch record created_utc is missing")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise RecoveryFailure("launch created_utc is invalid") from exc
    _require(parsed.tzinfo is not None, "launch created_utc lacks timezone")
    return parsed.timestamp()


def _boot_epoch() -> int:
    try:
        for line in PROC_STAT_PATH.read_text(encoding="ascii").splitlines():
            if line.startswith("btime "):
                return int(line.split()[1])
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        raise RecoveryFailure(f"cannot read Linux boot epoch: {exc}") from exc
    raise RecoveryFailure("Linux /proc/stat has no btime")


def _process_identity_observation(pid: Any, start_ticks: Any) -> dict[str, Any]:
    _require(isinstance(pid, int) and pid > 0 and
             isinstance(start_ticks, int) and start_ticks >= 0,
             "launch process identity is invalid")
    stat_path = Path("/proc") / str(pid) / "stat"
    try:
        text = stat_path.read_text(encoding="ascii")
    except FileNotFoundError:
        return {"pid": pid, "recorded_start_time_ticks": start_ticks,
                "state": "disappeared", "original_identity_gone": True}
    except OSError as exc:
        raise RecoveryFailure(f"cannot inspect {stat_path}: {exc}") from exc
    try:
        observed = int(text[text.rfind(")") + 2:].split()[19])
    except (ValueError, IndexError) as exc:
        raise RecoveryFailure(f"cannot parse {stat_path}") from exc
    if observed == start_ticks:
        return {"pid": pid, "recorded_start_time_ticks": start_ticks,
                "observed_start_time_ticks": observed, "state": "still_live",
                "original_identity_gone": False}
    return {"pid": pid, "recorded_start_time_ticks": start_ticks,
            "observed_start_time_ticks": observed, "state": "pid_reused",
            "original_identity_gone": True}


def _process_identity_gone(pid: Any, start_ticks: Any) -> dict[str, Any]:
    observation = _process_identity_observation(pid, start_ticks)
    _require(observation["original_identity_gone"] is True,
             f"original segment process identity is still live: pid {pid}")
    return observation


def _same_boot_closure(launch: dict[str, Any],
                       plan: dict[str, Any] | None = None) -> dict[str, Any]:
    boot_id = BOOT_ID_PATH.read_text(encoding="ascii").strip()
    _require(re.fullmatch(r"[0-9a-fA-F-]{36}", boot_id) is not None,
             "current Linux boot id is invalid")
    boot_epoch = _boot_epoch()
    launch_epoch = _parse_utc(launch.get("created_utc"))
    now = time.time()
    _require(boot_epoch <= launch_epoch <= now + 1.0,
             "boot changed after the recorded launch; same-boot closure is impossible")
    _require(launch.get("hostname") == socket.gethostname(),
             "launch hostname differs from the recovery host")
    holder = launch.get("input_transport")
    _require(isinstance(holder, dict),
             "launch record lacks launcher-holder identity")
    ranks = launch.get("ranks")
    _require(isinstance(ranks, list), "launch record lacks rank identities")
    identity_specs = [{
        "role": "mpirun", "pid": launch.get("mpirun_pid"),
        "recorded_start_time_ticks": launch.get("mpirun_start_time_ticks"),
    }, {
        "role": "launcher_holder", "pid": holder.get("holder_pid"),
        "recorded_start_time_ticks": holder.get("holder_start_time_ticks"),
    }]
    for rank in ranks:
        _require(isinstance(rank, dict), "launch rank identity is invalid")
        identity_specs.append({
            "role": f"rank:{rank.get('global_rank')}", "pid": rank.get("pid"),
            "recorded_start_time_ticks": rank.get("start_time_ticks"),
        })
    for identity in identity_specs:
        _require(isinstance(identity["pid"], int) and identity["pid"] > 0 and
                 isinstance(identity["recorded_start_time_ticks"], int) and
                 identity["recorded_start_time_ticks"] >= 0,
                 "launch process identity is invalid")
    managed_group = launch.get("managed_process_group")
    _require(isinstance(managed_group, dict) and
             managed_group.get("pgid") == launch.get("mpirun_pid"),
             "launch record lacks exact managed process group")
    deadline = time.monotonic() + SAME_BOOT_QUIESCENCE_TIMEOUT_SECONDS
    attempts = 0
    identities: list[dict[str, Any]] = []
    live_gpu: dict[str, Any] | None = None
    gpu_failure: str | None = None
    managed_group_gone = False
    while True:
        attempts += 1
        identities = []
        for identity in identity_specs:
            observed = _process_identity_observation(
                identity["pid"], identity["recorded_start_time_ticks"])
            observed["role"] = identity["role"]
            identities.append(observed)
        try:
            os.killpg(managed_group["pgid"], 0)
        except ProcessLookupError:
            managed_group_gone = True
        except PermissionError as exc:
            raise RecoveryFailure(
                "cannot prove original managed process group is gone") from exc
        else:
            managed_group_gone = False
        live_gpu = None
        gpu_failure = None
        if plan is not None:
            try:
                live_gpu = CHECKER.audit_live_gpu_quiescence(
                    plan, Path(plan["launch_contract"]["evidence"]["gpu_before"]))
            except Exception as exc:
                gpu_failure = str(exc)
        identities_gone = all(
            row["original_identity_gone"] is True for row in identities)
        if (identities_gone and managed_group_gone and
                (plan is None or live_gpu is not None)):
            break
        if time.monotonic() >= deadline:
            live_roles = [row["role"] for row in identities
                          if not row["original_identity_gone"]]
            detail = (f"live identities={live_roles!r}, "
                      f"managed_process_group_gone={managed_group_gone}, "
                      f"gpu={gpu_failure or 'quiet'}")
            raise RecoveryFailure(
                "same-boot closure did not become quiescent within "
                f"{SAME_BOOT_QUIESCENCE_TIMEOUT_SECONDS:g}s: {detail}")
        time.sleep(SAME_BOOT_QUIESCENCE_POLL_SECONDS)
    closed_artifacts: dict[str, Any] = {}
    if plan is not None:
        evidence = plan["launch_contract"]["evidence"]
        # A live launcher keeps run.log open.  Requiring a closed run log prevents
        # recovery from racing its normal completion publication.
        for name in ("run_log", "gpu_before"):
            closed_artifacts[name] = _binding(Path(evidence[name]))
        for name in ("exit_status", "gpu_after"):
            path = Path(evidence[name])
            if os.path.lexists(path):
                closed_artifacts[name] = _binding(path)
    return {
        "kind": "same_boot_closed_processes_v1",
        "boot_id": boot_id.lower(),
        "boot_epoch": boot_epoch,
        "launch_epoch": launch_epoch,
        "hostname": socket.gethostname(),
        "process_identities": identities,
        "all_original_identities_gone": True,
        "managed_process_group": managed_group,
        "managed_process_group_gone": managed_group_gone,
        "bounded_quiescence": {
            "kind": "poll_all_identities_group_and_gpu_v1",
            "timeout_seconds": SAME_BOOT_QUIESCENCE_TIMEOUT_SECONDS,
            "poll_seconds": SAME_BOOT_QUIESCENCE_POLL_SECONDS,
            "attempts": attempts,
            "all_quiet": True,
        },
        "closed_evidence_artifacts": closed_artifacts,
        "live_gpu_quiescence": live_gpu,
    }


def _nonzero_completion(path: Path, plan: dict[str, Any],
                        launch: dict[str, Any]) -> dict[str, Any]:
    record, binding = _immutable_json(path)
    _require(record.get("schema") == SCHEMA and
             record.get("kind") == "athenak_segment_completion" and
             record.get("status") == "ready" and
             isinstance(record.get("return_code"), int) and
             not isinstance(record.get("return_code"), bool) and
             record["return_code"] != 0,
             "recovery requires a nonzero completion record")
    _require(record.get("state_dir") == plan["launch_contract"]["state_dir"] and
             record.get("world_size") == plan["policy"]["ranks"],
             "completion state/world binding differs from original plan")
    expected_path = plan["launch_contract"]["evidence"]["completion_record"]
    _require(str(path.resolve(strict=True)) == str(Path(expected_path).resolve(strict=True)),
             "completion record differs from original plan")
    artifacts = record.get("artifacts")
    _require(isinstance(artifacts, dict) and set(artifacts) == {
        "plan", "launch_record", "run_log", "exit_status", "gpu_before", "gpu_after",
    }, "nonzero completion artifact set is not exact")
    evidence = plan["launch_contract"]["evidence"]
    expected_paths = {
        "plan": plan["launch_contract"]["plan_path"],
        "launch_record": evidence["launch_record"],
        "run_log": evidence["run_log"],
        "exit_status": evidence["exit_status"],
        "gpu_before": evidence["gpu_before"],
        "gpu_after": evidence["gpu_after"],
    }
    checked: dict[str, Any] = {}
    for name, planned in artifacts.items():
        _require(isinstance(planned, dict) and isinstance(planned.get("path"), str),
                 f"completion artifact {name} binding is invalid")
        _require(Path(planned["path"]).resolve(strict=True) ==
                 Path(expected_paths[name]).resolve(strict=True),
                 f"completion artifact {name} differs from original plan")
        fresh = _binding(Path(planned["path"]))
        _require(_same_binding(planned, fresh),
                 f"completion artifact {name} changed")
        checked[name] = fresh
    try:
        exit_value = int(Path(checked["exit_status"]["path"]).read_text(
            encoding="ascii").strip())
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        raise RecoveryFailure("completion exit status is invalid") from exc
    _require(exit_value == record["return_code"],
             "completion return code differs from exit.status")
    try:
        canonical = CHECKER.audit_completion_record(
            path, plan, Path(plan["launch_contract"]["state_dir"]), launch,
            checked, expected_return_code=exit_value)
        gpu_audit = CHECKER.audit_gpus(
            Path(checked["gpu_before"]["path"]),
            Path(checked["gpu_after"]["path"]),
            plan["policy"]["ranks"],
            plan["policy"]["gpu_exit_memory_mib_max"],
            plan["policy"]["gpu_ecc"])
        live_gpu = CHECKER.audit_live_gpu_quiescence(
            plan, Path(checked["gpu_before"]["path"]))
    except Exception as exc:
        raise RecoveryFailure(
            f"nonzero completion fails canonical lifecycle audit: {exc}") from exc
    mpirun_identity = _process_identity_gone(
        launch.get("mpirun_pid"), launch.get("mpirun_start_time_ticks"))
    mpirun_identity["role"] = "mpirun"
    live_identities = [mpirun_identity]
    holder = launch.get("input_transport")
    _require(isinstance(holder, dict),
             "launch record lacks launcher-holder identity")
    holder_identity = _process_identity_gone(
        holder.get("holder_pid"), holder.get("holder_start_time_ticks"))
    holder_identity["role"] = "launcher_holder"
    live_identities.append(holder_identity)
    for rank in launch.get("ranks", []):
        identity = _process_identity_gone(
            rank.get("pid"), rank.get("start_time_ticks"))
        identity["role"] = f"rank:{rank.get('global_rank')}"
        live_identities.append(identity)
    return {"kind": "nonzero_completion_v1", "return_code": exit_value,
            "completion_record": binding, "artifacts": checked,
            "quiescence": canonical["quiescence"],
            "canonical_completion_audit": canonical,
            "gpu_audit": gpu_audit,
            "live_gpu_quiescence": live_gpu,
            "live_identity_recheck": live_identities}


def _lifecycle(plan: dict[str, Any], launch: dict[str, Any],
               completion_path: Path) -> dict[str, Any]:
    expected = Path(
        plan["launch_contract"]["evidence"]["completion_record"]).absolute()
    _require(completion_path.expanduser().absolute() == expected,
             "completion path differs from the original plan")
    if os.path.lexists(completion_path):
        return _nonzero_completion(completion_path, plan, launch)
    return _same_boot_closure(launch, plan)


def _audit_run_log_prefix(path: Path, expected: dict[str, Any]) -> dict[str, Any]:
    try:
        lines = CHECKER._read_stable_bytes(path).decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise RecoveryFailure("run log is not UTF-8") from exc
    limits = [match for line in lines
              if (match := CHECKER.LIMIT_STATE_RE.fullmatch(line.strip()))]
    _require(len(limits) == 1 and int(limits[0].group(2)) == expected["final_cycle"] and
             float(limits[0].group(1)) == float(f"{expected['tlim']:.6e}"),
             "run log does not bind the original plan limits")
    diagnostics = [match for line in lines
                   if (match := CHECKER.DIAGNOSTIC_RE.fullmatch(line.strip()))]
    observed_cycles = [int(match.group(2)) for match in diagnostics]
    _require(observed_cycles == sorted(set(observed_cycles)),
             "run-log diagnostic cycles are duplicated or unordered")
    final_cycle = expected["recovery_final_cycle"]
    prefix = [match for match in diagnostics if int(match.group(2)) <= final_cycle]
    prefix_cycles = [int(match.group(2)) for match in prefix]
    wanted_cycles = list(range(expected["source_cycle"], final_cycle + 1))
    _require(prefix_cycles == wanted_cycles,
             "run log lacks an exact diagnostic prefix through the recovered restart")
    current = float(expected["source_time"])
    times = [current]
    for _ in range(final_cycle - expected["source_cycle"]):
        current += float(expected["root_dt"])
        times.append(current)
    _require([float(match.group(3)) for match in prefix] ==
             [float(f"{value:.6e}") for value in times] and
             [float(match.group(4)) for match in prefix] ==
             [float(f"{expected['root_dt']:.6e}")] * len(prefix),
             "run-log prefix does not follow fixed sequential root_dt")
    elapsed = [float(match.group(1)) for match in prefix]
    _require(all(math.isfinite(value) and value >= 0.0 for value in elapsed) and
             all(right >= left for left, right in zip(elapsed, elapsed[1:])),
             "run-log prefix elapsed values are invalid")
    caches = [match for line in lines
              if (match := CHECKER.CACHE_RE.fullmatch(line.strip()))]
    _require(len(caches) == 1,
             "run log must contain one restart-cache qualification")
    solver_failures, nonfinite = int(caches[0].group(1)), int(caches[0].group(2))
    raw_relative, absolute, mixed, tolerance = map(float, caches[0].groups()[2:])
    _require(solver_failures == 0 and nonfinite == 0 and
             all(math.isfinite(value) and value >= 0.0 for value in
                 (raw_relative, absolute, mixed, tolerance)) and
             tolerance > 0.0 and mixed <= tolerance,
             "run-log restart-cache qualification failed")
    return {
        "original_nlim": int(limits[0].group(2)),
        "original_tlim": float(limits[0].group(1)),
        "root_step_diagnostics": {
            "rows": len(prefix), "cycle_min": prefix_cycles[0],
            "cycle_max": prefix_cycles[-1], "fixed_dt": expected["root_dt"],
            "all_recovered_prefix_cycles_present": True,
        },
        "cache": {
            "solver_failures": solver_failures,
            "nonfinite_proposed_values": nonfinite,
            "max_raw_component_relative_change": raw_relative,
            "max_absolute_change": absolute,
            "max_mixed_scale_change": mixed,
            "mixed_scale_tolerance": tolerance,
        },
    }


def _prefix_text_audits(plan: dict[str, Any], actual: dict[str, Path],
                        root_steps: int, source_cycle: int, final_cycle: int,
                        source_time: float, final_time: float,
                        source_mass: float) \
        -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any],
                 dict[str, Any]]:
    derived: list[dict[str, Any]] = []
    histories: list[dict[str, Any]] = []
    event_result: dict[str, Any] | None = None
    baryon_result: dict[str, Any] | None = None
    for output in plan["outputs"]:
        if output["numbered"]:
            continue
        for relative in output["required_unnumbered_paths"]:
            _require(relative in actual, f"required unnumbered output is missing: {relative}")
            source_binding = _binding(actual[relative])
            raw = CHECKER._read_stable_bytes(actual[relative])
            _require(hashlib.sha256(raw).hexdigest() == source_binding["sha256"],
                     f"unnumbered output changed while deriving prefix: {relative}")
            prefix, suffix = _exact_text_prefix(raw, root_steps, relative)
            forensic = {
                "relative_path": relative,
                "source": source_binding,
                "prefix_size": len(prefix),
                "prefix_sha256": hashlib.sha256(prefix).hexdigest(),
                "suffix_size": len(suffix),
                "suffix_sha256": hashlib.sha256(suffix).hexdigest(),
            }
            if output["file_type"] == "hst":
                audited = _audit_text_bytes(prefix, lambda path: CHECKER.audit_history(
                    path, source_cycle, final_cycle, source_time,
                    plan["expected"]["root_dt"], final_time,
                    plan["policy"]["endpoint_time_ulp_tolerance"]))
                histories.append({"relative_path": relative, "history": audited,
                                  "prefix": prefix, "forensic": forensic})
            elif output["file_type"] == "log":
                _require(event_result is None,
                         "plan has more than one authoritative event log")
                event_result = _audit_text_bytes(prefix, lambda path:
                    CHECKER.audit_event_log(
                        path, source_cycle, final_cycle,
                        plan["policy"]["event_thresholds"],
                        plan["policy"]["event_absolute_thresholds"]))
                derived.append({"relative_path": relative, "prefix": prefix,
                                "forensic": forensic, "event": event_result,
                                "file_type": "log", "block": output["block"]})
            else:
                raise RecoveryFailure(
                    f"unsupported unnumbered recovery output type: {output['file_type']}")
    baryon = [row for row in histories
              if "baryon_m" in row["history"].get("columns", [])]
    _require(len(baryon) == 1, "recovery prefix must have one baryon history")
    baryon_result = CHECKER.audit_baryon_mass(
        baryon[0]["history"],
        plan["policy"]["baryon_mass_fractional_loss"]["hard_per_root_step"],
        root_steps, source_mass)
    _require(event_result is not None, "recovery prefix lacks an event log")
    for row in histories:
        derived.append({"relative_path": row["relative_path"],
                        "prefix": row["prefix"], "forensic": row["forensic"],
                        "history": row["history"], "file_type": "hst"})
    return derived, histories, event_result, baryon_result


def _audit_source_chain(plan: dict[str, Any], source_binding: dict[str, Any],
                        source_cycle: int, source_time: float,
                        source_mass: float) \
        -> tuple[dict[str, Any], dict[str, Any] | None]:
    qualification = plan["source_qualification"]
    mass_evidence = CHECKER._verify_planned_file(
        qualification["source_baryon_mass"]["evidence"],
        "source baryon-mass evidence")
    if qualification["mode"] == "anchor_full_audit":
        return ({"mode": "anchor_full_audit",
                 "source_baryon_mass_evidence": mass_evidence}, None)
    parent_record = qualification["parent_segment_pass"]
    parent_binding = CHECKER._verify_planned_file(
        parent_record, "source parent segment pass")
    parent_path = Path(parent_binding["path"])
    CHECKER._require_immutable_ready(
        parent_path, ".pass.ready", "source parent segment pass")
    parent, _ = CHECKER._load_json(parent_path)
    endpoint = parent.get("bindings", {}).get("endpoint_restart")
    trajectory = parent.get("bindings", {}).get("trajectory")
    _require(isinstance(endpoint, dict) and all(
        endpoint.get(name) == source_binding.get(name)
        for name in ("path", "size", "sha256")),
        "source parent pass does not bind the live source restart")
    _require(isinstance(trajectory, dict) and
             trajectory.get("size") == plan["inputs"]["trajectory"]["size"] and
             trajectory.get("sha256") == plan["inputs"]["trajectory"]["sha256"],
             "source parent pass trajectory differs from the plan")
    histories = [row for row in parent.get("output_inventory", [])
                 if isinstance(row, dict) and
                 isinstance(row.get("history"), dict) and
                 "baryon_m" in row["history"].get("columns", [])]
    _require(len(histories) == 1 and all(
        histories[0].get(name) == mass_evidence.get(name)
        for name in ("path", "size", "sha256")),
        "source parent pass does not bind the source baryon history")
    try:
        parent_mass = float(
            parent["scientific_threshold_audit"]["baryon_mass"]["last"])
    except (KeyError, TypeError, ValueError) as exc:
        raise RecoveryFailure("source parent pass lacks endpoint baryon mass") from exc
    _require(math.isfinite(parent_mass) and parent_mass == source_mass,
             "source parent baryon mass differs from live evidence")
    provenance = CHECKER.audit_parent_qualification_provenance(
        parent, qualification, source_binding, mass_evidence,
        source_cycle, source_time, parent_path)
    return ({
        "mode": "parent_segment_pass", "parent_pass": parent_binding,
        "source_baryon_mass_evidence": mass_evidence,
        "parent_qualification": provenance,
    }, parent)


def _audit_numbered_prefix(plan: dict[str, Any], actual: dict[str, Path],
                           schedules: dict[str, list[dict[str, Any]]],
                           source_metadata: Any, endpoint_metadata: Any,
                           final_cycle: int, final_time: float,
                           history_audits: list[dict[str, Any]]) \
        -> tuple[list[dict[str, Any]], list[dict[str, Any]], float,
                 list[dict[str, Any]]]:
    inventory: list[dict[str, Any]] = []
    binary_by_block: dict[str, list[dict[str, Any]]] = {}
    divb: list[dict[str, Any]] = []
    for output in plan["outputs"]:
        if not output["numbered"]:
            continue
        for write in schedules[output["block"]]:
            relative = _relative_for(output, write)
            _require(relative in actual,
                     f"logical prefix is missing planned output: {relative}")
            path = actual[relative]
            if output["file_type"] == "rst":
                audited_restart = audit_restart(path)
                _require(audited_restart["metadata"]["cycle"] == write["cycle"] and
                         audited_restart["metadata"]["time"] == write["time"],
                         f"restart {relative} differs from Execute replay")
                audited: dict[str, Any] = {
                    **_binding(path), "restart_audit": audited_restart,
                }
            elif output["file_type"] == "bin":
                audited = CHECKER.audit_binary(
                    path, plan["expected"]["source_cycle"], final_cycle,
                    plan["expected"]["source_time"], final_time,
                    output["parameters"], endpoint_metadata,
                    source_metadata.parameters, output["block"],
                    output["expected_binary_variables"])
                _require(audited["cycle"] == write["cycle"] and
                         audited["time"] == write["time"],
                         f"binary {relative} differs from Execute replay")
                binary_by_block.setdefault(output["block"], []).append(audited)
                for variable, maximum in audited["variable_max_abs"].items():
                    if variable.lower() == "divb":
                        divb.append({"path": str(path), "cycle": audited["cycle"],
                                     "time": audited["time"], "max_abs": maximum})
            else:
                audited = _binding(path)
            inventory.append({
                "block": output["block"], "file_type": output["file_type"],
                "file_number": write["file_number"], "expected_write": write,
                **audited,
            })
    CHECKER.audit_binary_pairs(plan, binary_by_block)
    bbh = [row["history"] for row in history_audits
           if {"bh1_x", "bh2_x", "bh1_mass", "bh2_mass"}.issubset(
               row["history"].get("columns", []))]
    _require(len(bbh) == 1,
             "recovery prefix must have one BBH-center history")
    centers = CHECKER._bbh_history_centers(bbh[0])
    full_by_cycle: dict[int, dict[str, Any]] = {}
    for rows in binary_by_block.values():
        for audited in rows:
            if audited["topology"]["scope"] != "full_domain":
                continue
            previous = full_by_cycle.get(audited["cycle"])
            if previous is not None:
                _require(previous["time"] == audited["time"] and
                         previous["_topology_records"] ==
                         audited["_topology_records"],
                         "same-cycle full-domain binary topologies disagree")
            else:
                full_by_cycle[audited["cycle"]] = audited
    output_by_block = {row["block"]: row for row in plan["outputs"]}
    coverage: list[dict[str, Any]] = []
    for block, rows in binary_by_block.items():
        parameters = output_by_block[block]["parameters"]
        for audited in rows:
            result = CHECKER.audit_selected_binary_coverage(
                audited, parameters, centers.get(audited["cycle"]),
                full_by_cycle.get(audited["cycle"]))
            audited["selection_coverage"] = result
            coverage.append({"path": audited["path"], "cycle": audited["cycle"],
                             **result})
    _require(divb, "recovery prefix contains no audited divB output")
    worst = max(row["max_abs"] for row in divb)
    hard = float(plan["policy"]["divb_max_abs"]["hard"])
    _require(worst <= hard, f"prefix divB maximum {worst} exceeds {hard}")
    coverage_by_path = {row["path"]: row for row in coverage}
    for row in inventory:
        if row.get("path") in coverage_by_path:
            row["selection_coverage"] = {
                key: value for key, value in coverage_by_path[row["path"]].items()
                if key != "path"
            }
        row.pop("_logical_locations", None)
        row.pop("_topology_records", None)
    return inventory, divb, worst, coverage


def recover(args: argparse.Namespace) -> dict[str, Any]:
    plan, plan_binding = _immutable_json(args.plan)
    try:
        expected = CHECKER.validate_plan(plan)
    except Exception as exc:
        raise RecoveryFailure(f"original plan is not a current strict plan: {exc}") \
            from exc
    _require(plan.get("policy", {}).get("scheduled_prefix_recovery") ==
             RECOVERY_POLICY,
             "original plan lacks the exact scheduled-prefix recovery policy")
    recovery_tool = _planned_file(plan, "prefix_recovery", SCRIPT_PATH)
    checker_tool = _planned_file(plan, "segment_checker",
                                 SCRIPT_DIRECTORY / "check_athenak_segment.py")
    integrity_tool = _planned_file(plan, "output_integrity",
                                   SCRIPT_DIRECTORY / "output_integrity.py")
    auditor_tool = _planned_file(plan, "restart_auditor",
                                 SCRIPT_DIRECTORY / "audit_athenak_restart.py")
    reader_tool = _planned_file(plan, "restart_metadata_reader",
                                SCRIPT_DIRECTORY / "read_athenak_restart_metadata.py")
    planned_tools = {
        name: CHECKER._verify_planned_file(plan["tools"][name],
                                           f"original plan tool {name}")
        for name in CHECKER.PLANNED_TOOL_NAMES
    }
    repository = CHECKER._verify_planned_repository(
        plan["inputs"]["repo"], plan["tools"]["git"])
    binary_binding = CHECKER._verify_planned_file(
        plan["inputs"]["binary"], "AthenaK binary")
    trajectory_binding = CHECKER._verify_planned_file(
        plan["inputs"]["trajectory"], "trajectory")

    original_state = args.state_dir.expanduser().absolute().resolve(strict=True)
    original_evidence = Path(plan["launch_contract"]["evidence_dir"]).resolve(
        strict=True)
    _require(str(original_state) == plan["launch_contract"]["state_dir"],
             "--state-dir differs from original plan")
    original_state_binding = BoundDirectory.open(
        original_state, "original state directory")
    original_evidence_binding = BoundDirectory.open(
        original_evidence, "original evidence directory")
    recovery_state = _private_empty_directory(
        args.recovery_state_dir, "recovery state directory")
    recovery_evidence = _private_empty_directory(
        args.recovery_evidence_dir, "recovery evidence directory")
    recovery_state_binding = BoundDirectory.open(
        recovery_state, "recovery state directory", required_mode=0o700,
        require_empty=True)
    recovery_evidence_binding = BoundDirectory.open(
        recovery_evidence, "recovery evidence directory", required_mode=0o700,
        require_empty=True)
    _separate_directories(
        (original_state, original_evidence, recovery_state, recovery_evidence))
    output = args.output.expanduser().absolute()
    record_path = recovery_evidence / "segment.prefix.recovery.ready"
    _require(output == recovery_evidence / "segment.prefix.pass.ready",
             "--output must be recovery evidence/segment.prefix.pass.ready")
    _require(args.plan.resolve(strict=True) ==
             Path(plan["launch_contract"]["plan_path"]).resolve(strict=True),
             "--plan differs from original plan-bound path")
    planned_evidence_paths = {
        "plan": Path(plan["launch_contract"]["plan_path"]).resolve(strict=True),
        **{
            name: Path(value).resolve(strict=True)
            for name, value in plan["launch_contract"]["evidence"].items()
            if os.path.lexists(value)
        },
    }
    for name, path in planned_evidence_paths.items():
        _require(path.parent == original_evidence,
                 f"planned original evidence {name} is outside evidence_dir")
    all_evidence_names = {
        Path(plan["launch_contract"]["plan_path"]).name,
        *(Path(value).name
          for value in plan["launch_contract"]["evidence"].values()),
    }
    _require(len(all_evidence_names) ==
             1 + len(plan["launch_contract"]["evidence"]),
             "original evidence plan maps duplicate basenames")
    initial_evidence = _inventory_original_tree(
        original_evidence, all_evidence_names, original_evidence_binding)
    required_evidence_names = {
        Path(plan["launch_contract"]["plan_path"]).name,
        *(Path(plan["launch_contract"]["evidence"][name]).name
          for name in ("launch_record", "run_log", "gpu_before")),
    }
    _require(required_evidence_names.issubset(initial_evidence),
             "original evidence lacks plan/launch/run-log/GPU-before")
    initial_evidence_bindings = {
        relative: _binding(path) for relative, path in initial_evidence.items()
    }
    launch_path = Path(plan["launch_contract"]["evidence"]["launch_record"])
    _require(args.launch_record.resolve(strict=True) == launch_path.resolve(strict=True),
             "--launch-record differs from original plan")
    launch, launch_binding = _immutable_json(args.launch_record)
    _require(_same_binding(
        initial_evidence_bindings[args.plan.name], plan_binding) and
             _same_binding(
                 initial_evidence_bindings[args.launch_record.name],
                 launch_binding),
             "plan/launch binding differs from bound original evidence tree")
    gpu_before = _binding(Path(plan["launch_contract"]["evidence"]["gpu_before"]))
    try:
        CHECKER.audit_launch_record(args.launch_record, plan, args.plan,
                                    original_state, gpu_before)
    except Exception as exc:
        raise RecoveryFailure(f"original launch record is invalid: {exc}") from exc
    completion = (args.completion_record if args.completion_record is not None else
                  Path(plan["launch_contract"]["evidence"]["completion_record"]))
    lifecycle = _lifecycle(plan, launch, completion)
    lifecycle_evidence = (lifecycle.get("artifacts", {})
                          if lifecycle["kind"] == "nonzero_completion_v1"
                          else lifecycle.get("closed_evidence_artifacts", {}))
    for name, artifact in lifecycle_evidence.items():
        basename = Path(artifact["path"]).name
        _require(basename in initial_evidence_bindings and
                 _same_binding(artifact, initial_evidence_bindings[basename]),
                 f"lifecycle {name} binding differs from original evidence tree")

    known, numbered = _known_paths(plan)
    actual = _inventory_original_tree(
        original_state, known, original_state_binding)
    initial_original_bindings = {
        relative: _binding(path) for relative, path in actual.items()
    }
    minimum_age = float(
        plan["policy"]["output_integrity"]["minimum_closed_file_age_seconds"])
    _require(math.isfinite(minimum_age) and minimum_age >= 0.0,
             "original plan has invalid closed-file settling interval")
    for relative, artifact in initial_original_bindings.items():
        _require_settled(artifact, minimum_age,
                         f"original state artifact {relative}")
    if lifecycle["kind"] == "nonzero_completion_v1":
        _require_settled(lifecycle["completion_record"], minimum_age,
                         "nonzero completion record")
        for name, artifact in lifecycle["artifacts"].items():
            _require_settled(artifact, minimum_age,
                             f"completion artifact {name}")
    else:
        for name, artifact in lifecycle["closed_evidence_artifacts"].items():
            _require_settled(artifact, minimum_age,
                             f"same-boot evidence artifact {name}")
    source_path = Path(plan["inputs"]["source_restart"]["path"])
    source_metadata = read_restart_metadata(source_path)
    minimum_restart_header = source_metadata.metadata_end + struct.calcsize("Q")
    restart_output = _restart_output(plan)
    scheduled_restarts = [row for row in restart_output["expected_writes"]
                          if row.get("kind") == "scheduled"]
    _require(scheduled_restarts, "original plan contains no scheduled restart")
    candidates: list[dict[str, Any]] = []
    for write in sorted(scheduled_restarts, key=lambda row: row["cycle"], reverse=True):
        relative = _relative_for(restart_output, write)
        if relative not in actual:
            candidates.append({"relative_path": relative, "expected_write": write,
                               "classification": "absent"})
            continue
        classified = classify_restart(
            actual[relative], minimum_complete_header_bytes=minimum_restart_header)
        candidates.append({"relative_path": relative, "expected_write": write,
                           **classified})
    complete = select_highest_complete_candidate(candidates)
    selected_write = complete["expected_write"]
    final_cycle = int(selected_write["cycle"])
    schedules, states, final_time = replay_execute_prefix(plan, final_cycle)
    _require(selected_write["time"] == final_time,
             "selected restart time differs from Execute replay")
    endpoint_metadata = _assert_execute_endpoint(
        plan, launch, complete["audit"], final_cycle, final_time, states)

    source_audit = audit_restart(source_path)
    _require(source_audit["sha256"] == plan["inputs"]["source_restart"]["sha256"],
             "live source restart differs from original plan")
    source_binding = {
        "path": source_audit["path"], **source_audit["signature"],
        "sha256": source_audit["sha256"],
        "closure_check": source_audit["closure_check"],
    }
    partial_expected = dict(expected)
    partial_expected.update({"final_cycle": final_cycle, "final_time": final_time,
                             "root_steps": final_cycle - expected["source_cycle"]})
    run_log_expected = dict(expected)
    run_log_expected["recovery_final_cycle"] = final_cycle
    run_log_prefix_audit = _audit_run_log_prefix(
        Path(plan["launch_contract"]["evidence"]["run_log"]),
        run_log_expected)
    run_log_name = Path(plan["launch_contract"]["evidence"]["run_log"]).name
    _require(_same_binding(
        _binding(Path(plan["launch_contract"]["evidence"]["run_log"])),
        initial_evidence_bindings[run_log_name]),
        "run log changed between lifecycle binding and prefix audit")
    CHECKER.audit_restart_contract(source_metadata, "source", plan, partial_expected)
    source_mass = CHECKER.audit_source_baryon_evidence(
        plan["source_qualification"], expected["source_time"],
        ({} if plan["source_qualification"]["mode"] == "parent_segment_pass"
         else None))
    source_chain, source_parent = _audit_source_chain(
        plan, source_binding, expected["source_cycle"],
        expected["source_time"], source_mass)

    derived, histories, events, baryon = _prefix_text_audits(
        plan, actual, partial_expected["root_steps"], expected["source_cycle"],
        final_cycle, expected["source_time"], final_time, source_mass)
    numbered_inventory, divb_files, worst_divb, selection_coverage = \
        _audit_numbered_prefix(
        plan, actual, schedules, source_metadata, endpoint_metadata,
        final_cycle, final_time, histories)
    baryon_histories = [row["history"] for row in histories
                        if "baryon_m" in row["history"].get("columns", [])]
    _require(len(baryon_histories) == 1,
             "recovery prefix lacks one baryon advisory history")
    baryon_policy = plan["policy"]["baryon_mass_fractional_loss"]
    event_advisories = CHECKER.audit_event_ratio_advisories(
        events, plan["policy"]["yellow_event_thresholds"])
    baryon_advisory = CHECKER.audit_baryon_mass_advisory(
        baryon_histories[0], source_mass, expected["source_cycle"],
        baryon_policy["yellow_per_root_step"],
        baryon_policy["yellow_per_48M"],
        baryon_policy["rolling_window_root_steps"])
    floor_trends = CHECKER.audit_floor_rate_trends(events, source_parent)
    divb_yellow = float(plan["policy"]["divb_max_abs"]["yellow"])
    worst_divb_file = max(divb_files, key=lambda row: row["max_abs"])
    divb_advisory = {
        "severity": "yellow" if worst_divb > divb_yellow else "green",
        "yellow_max_abs_exclusive_min": divb_yellow,
        "observed_max_abs": worst_divb,
        "path": worst_divb_file["path"],
        "cycle": worst_divb_file["cycle"],
        "time": worst_divb_file["time"],
    }
    scientific_advisories = {
        "schema": "athenak_scientific_advisories_v1",
        "severity": (
            "yellow" if event_advisories["severity"] == "yellow" or
            baryon_advisory["severity"] == "yellow" or
            divb_advisory["severity"] == "yellow" else "green"),
        "event_ratios": event_advisories,
        "baryon_mass": baryon_advisory,
        "divb": divb_advisory,
        "floor_rates": floor_trends,
        "pass_fail_effect": "none_yellow_advisories_are_nonfatal",
    }
    events.pop("_rows", None)

    prefix_paths = {
        _relative_for(output, write)
        for output in plan["outputs"] if output["numbered"]
        for write in schedules[output["block"]]
    } | set(plan["required_files"])
    suffix_forensics: list[dict[str, Any]] = []
    for relative, path in sorted(actual.items()):
        if relative in prefix_paths:
            continue
        logical = numbered.get(relative)
        binding = _binding(path)
        suffix_forensics.append({
            "relative_path": relative,
            "classification": "planned_post_prefix_artifact",
            "logical_write": logical["write"] if logical else None,
            "binding": binding,
        })

    derived_inventory: list[dict[str, Any]] = []
    derived_forensics: list[dict[str, Any]] = []
    for row in derived:
        installed = _install_prefix_bytes(
            recovery_state, row["relative_path"], row["prefix"],
            recovery_state_binding)
        forensic = dict(row["forensic"])
        forensic["derived"] = installed
        derived_forensics.append(forensic)
        inventory_row: dict[str, Any] = {
            "block": row.get("block") or next(
                output["block"] for output in plan["outputs"]
                if row["relative_path"] in output["required_unnumbered_paths"]),
            "file_type": row["file_type"], "file_number": None,
            "expected_write": None, **installed,
        }
        if "history" in row:
            history = dict(row["history"])
            history.pop("_row_values", None)
            inventory_row["history"] = history
        if "event" in row:
            inventory_row["event_log"] = row["event"]
        derived_inventory.append(inventory_row)

    endpoint_binding = {
        "path": complete["audit"]["path"],
        **complete["audit"]["signature"],
        "sha256": complete["audit"]["sha256"],
        "closure_check": complete["audit"]["closure_check"],
    }
    logical_prefix = {
        "source_cycle": expected["source_cycle"],
        "source_time": expected["source_time"],
        "final_cycle": final_cycle,
        "final_time": final_time,
        "root_steps": partial_expected["root_steps"],
        "execute_only": True,
        "expected_numbered_paths": sorted(prefix_paths - set(plan["required_files"])),
        "required_unnumbered_paths": sorted(plan["required_files"]),
        "all_expected_prefix_artifacts_present": True,
    }
    record = {
        "schema": SCHEMA,
        "kind": "athenak_segment_prefix_recovery",
        "status": "prepared",
        "publication_transaction": {
            "kind": "prepared_record_then_pass_commit_v1",
            "commit_filename": "segment.prefix.pass.ready",
            "prepared_record_alone_is_not_consumable": True,
        },
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "qualification_mode": QUALIFICATION_MODE,
        "policy": RECOVERY_POLICY,
        "original_plan": plan_binding,
        "original_launch_record": launch_binding,
        "lifecycle": lifecycle,
        "selected_scheduled_restart": {
            "block": restart_output["block"],
            "relative_path": complete["relative_path"],
            "expected_write": selected_write,
            "binding": endpoint_binding,
        },
        "candidate_inventory": candidates,
        "logical_prefix": logical_prefix,
        "run_log_prefix_audit": run_log_prefix_audit,
        "derived_text_prefixes": derived_forensics,
        "suffix_forensics": suffix_forensics,
        "original_state_tree": {
            "path": str(original_state), "read_only_access": True,
            "regular_files": len(actual), "unknown_nodes": 0,
            "directory_binding": original_state_binding.audit(
                "original state directory"),
        },
        "original_evidence_directory_binding": original_evidence_binding.audit(
            "original evidence directory"),
        "recovery_state_directory_binding": recovery_state_binding.audit(
            "recovery state directory"),
        "recovery_evidence_directory_binding": recovery_evidence_binding.audit(
            "recovery evidence directory"),
        "recovery_state_dir": str(recovery_state),
        "recovery_evidence_dir": str(recovery_evidence),
        "tools": {
            **planned_tools,
            "prefix_recovery": recovery_tool,
            "segment_checker": checker_tool,
            "output_integrity": integrity_tool,
            "restart_auditor": auditor_tool,
            "restart_metadata_reader": reader_tool,
        },
        "repository": repository,
        "binary": binary_binding,
        "trajectory": trajectory_binding,
        "source_chain": source_chain,
    }
    # Complete every fallible source/state/evidence proof before publishing the
    # prepared transaction record.  A crash after this point can leave only an
    # explicitly non-consumable ``prepared`` record, never a pass-shaped orphan.
    recovery_actual = _inventory_original_tree(
        recovery_state, set(plan["required_files"]), recovery_state_binding)
    recovery_by_relative = {
        row["relative_path"]: row["derived"] for row in derived_forensics
    }
    _require(set(recovery_actual) == set(recovery_by_relative),
             "recovery state tree differs from exact derived-prefix inventory")
    for relative, path in recovery_actual.items():
        _require(_same_binding(_binding(path), recovery_by_relative[relative]),
                 f"derived prefix changed before publication: {relative}")
    _require(_same_binding(_binding(args.plan), plan_binding),
             "original plan changed during recovery")
    _require(_same_binding(_binding(args.launch_record), launch_binding),
             "original launch record changed during recovery")
    _require(_same_binding(_binding(Path(endpoint_binding["path"])), endpoint_binding),
             "selected restart changed during recovery")
    final_original = _inventory_original_tree(
        original_state, known, original_state_binding)
    _require(set(final_original) == set(actual),
             "original state tree changed during recovery")
    for relative, path in final_original.items():
        _require(_same_binding(_binding(path), initial_original_bindings[relative]),
                 f"original state artifact changed during recovery: {relative}")
    final_evidence = _inventory_original_tree(
        original_evidence, all_evidence_names, original_evidence_binding)
    _require(set(final_evidence) == set(initial_evidence),
             "original evidence tree changed during recovery")
    for relative, path in final_evidence.items():
        _require(_same_binding(_binding(path), initial_evidence_bindings[relative]),
                 f"original evidence changed during recovery: {relative}")
    original_state_binding.audit("original state directory")
    original_evidence_binding.audit("original evidence directory")
    recovery_state_binding.audit("recovery state directory")
    recovery_evidence_binding.audit("recovery evidence directory", require_empty=True)
    planned_completion = Path(
        plan["launch_contract"]["evidence"]["completion_record"])
    if lifecycle["kind"] == "same_boot_closed_processes_v1":
        _require(not os.path.lexists(planned_completion),
                 "normal completion appeared during same-boot recovery; rerun "
                 "qualification from that immutable lifecycle")
    recovery_evidence_binding.audit("recovery evidence directory")
    try:
        record_binding = _install_bound_json(
            recovery_evidence_binding, record_path.name, record)
    except (OSError, RuntimeError, ValueError) as exc:
        raise RecoveryFailure(f"cannot publish recovery record: {exc}") from exc

    baryon_history_rows = [row for row in derived_inventory
                           if isinstance(row.get("history"), dict) and
                           "baryon_m" in row["history"].get("columns", [])]
    _require(len(baryon_history_rows) == 1,
             "derived inventory lacks unique baryon history")
    report = {
        "schema": SCHEMA,
        "kind": "athenak_segment_pass",
        "status": "pass",
        "qualification_mode": QUALIFICATION_MODE,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "expected": partial_expected,
        "bindings": {
            "original_plan": plan_binding,
            "original_launch_record": launch_binding,
            "source_restart": source_binding,
            "endpoint_restart": endpoint_binding,
            "trajectory": trajectory_binding,
            "recovery_record": record_binding,
            "tools": record["tools"],
        },
        "source_restart_audit": source_audit,
        "endpoint_restart_audit": complete["audit"],
        "restart_contract_audit": {
            "source": CHECKER.audit_restart_contract(
                source_metadata, "source", plan, partial_expected),
            "endpoint": CHECKER.audit_restart_contract(
                endpoint_metadata, "recovered endpoint", plan, partial_expected),
        },
        "event_log_audit": events,
        "run_log_prefix_audit": run_log_prefix_audit,
        "scientific_threshold_audit": {
            "nonfinite_count": 0,
            "divb": {
                "hard_max_abs": plan["policy"]["divb_max_abs"]["hard"],
                "observed_max_abs": worst_divb, "files": divb_files,
            },
            "baryon_mass": {
                "path": baryon_history_rows[0]["path"], **baryon,
            },
        },
        "scientific_advisories": scientific_advisories,
        "binary_selection_coverage_audit": selection_coverage,
        "output_inventory": [*numbered_inventory, *derived_inventory],
        "recovery_provenance": {
            "kind": QUALIFICATION_MODE,
            "policy": RECOVERY_POLICY,
            "record": record_binding,
            "original_plan": plan_binding,
            "original_launch_record": launch_binding,
            "selected_scheduled_restart": record["selected_scheduled_restart"],
            "logical_prefix": logical_prefix,
            "run_log_prefix_audit": run_log_prefix_audit,
            "lifecycle": lifecycle,
            "derived_text_prefixes": derived_forensics,
            "suffix_forensics": suffix_forensics,
            "original_trees_unchanged": True,
        },
    }
    original_state_binding.close()
    original_evidence_binding.close()
    recovery_state_binding.close()
    recovery_evidence_binding.close()
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--state-dir", required=True, type=Path)
    parser.add_argument("--launch-record", required=True, type=Path)
    parser.add_argument("--completion-record", type=Path)
    parser.add_argument("--recovery-state-dir", required=True, type=Path)
    parser.add_argument("--recovery-evidence-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        report = recover(args)
        evidence_root = args.recovery_evidence_dir.expanduser().absolute()
        bound = BoundDirectory.open(
            evidence_root, "recovery evidence directory", required_mode=0o700)
        try:
            bound.audit("recovery evidence directory")
            _require(args.output.expanduser().absolute().parent == bound.path,
                     "report output is not a direct recovery-evidence child")
            _require(set(os.listdir(bound.fd)) == {
                "segment.prefix.recovery.ready"},
                "recovery evidence directory differs from pre-publication allowlist")
            if report["recovery_provenance"]["lifecycle"]["kind"] == \
                    "same_boot_closed_processes_v1":
                plan_path = Path(report["bindings"]["original_plan"]["path"])
                original_plan, _ = _immutable_json(plan_path)
                completion_path = Path(
                    original_plan["launch_contract"]["evidence"]
                    ["completion_record"])
                _require(not os.path.lexists(completion_path),
                         "normal completion appeared before recovery pass commit")
            _install_bound_json(bound, args.output.name, report)
            bound.audit("recovery evidence directory")
            _require(set(os.listdir(bound.fd)) == {
                "segment.prefix.recovery.ready", "segment.prefix.pass.ready"},
                "recovery evidence directory differs from final allowlist")
        finally:
            bound.close()
    except (OSError, RuntimeError, ValueError, RecoveryFailure) as exc:
        print(json.dumps({
            "schema": SCHEMA, "kind": "athenak_segment_prefix_recovery",
            "status": "fail", "message": str(exc),
        }, sort_keys=True), file=sys.stderr)
        return 1
    print(json.dumps({
        "schema": SCHEMA, "kind": "athenak_segment_prefix_recovery",
        "status": "pass", "report": str(args.output.resolve()),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
