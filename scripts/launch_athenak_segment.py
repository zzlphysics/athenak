#!/usr/bin/env python3
"""Launch exactly one immutable, cycle-limited AthenaK GPU segment.

This is deliberately a launcher, not a monitor.  It proves the launch once,
publishes a closed ``segment.launch.ready`` record, waits for the MPI job it
started, and publishes the exit/GPU evidence.  It never chains another segment
and makes no runtime health or shutdown decision after launch qualification.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import socket
import stat
import struct
import subprocess
import sys
import time
from typing import Any, Callable, Mapping, Protocol, Sequence


sys.dont_write_bytecode = True
SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIRECTORY = Path(__file__).resolve().parent
OUTPUT_INTEGRITY_PATH = (SCRIPT_DIRECTORY / "output_integrity.py").resolve()
SEGMENT_CHECKER_PATH = (SCRIPT_DIRECTORY / "check_athenak_segment.py").resolve()
RESTART_AUDITOR_PATH = (SCRIPT_DIRECTORY / "audit_athenak_restart.py").resolve()
RESTART_READER_PATH = (
    SCRIPT_DIRECTORY / "read_athenak_restart_metadata.py"
).resolve()
if str(SCRIPT_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

from output_integrity import (
    HASH_CHUNK_BYTES,
    _assert_closed,
    _assert_path_signature,
    _assert_stream_signature,
    _open_regular_nofollow,
    install_immutable_json,
    stable_sha256,
)
from audit_athenak_restart import audit_restart
from check_athenak_segment import (
    CheckFailure as PlanCheckFailure,
    validate_plan as validate_strict_plan,
)
from read_athenak_restart_metadata import (
    metadata_dict as restart_metadata_dict,
    read_restart_metadata,
)


SCHEMA = 1
MAX_NMB_PER_RANK = 16384
CADENCE_MULTIPLE_MAX_ULPS = 2
CAPACITY_MEMORY_MIB_PER_SLOT_NUMERATOR = 1433
CAPACITY_MEMORY_MIB_PER_SLOT_DENOMINATOR = 100
CAPACITY_GPU_USABLE_FRACTION_NUMERATOR = 4
CAPACITY_GPU_USABLE_FRACTION_DENOMINATOR = 5
GIB_BYTES = 1 << 30
DISK_PREFLIGHT_KIND = "statvfs_unique_filesystem_budget_v1"
DISK_PREFLIGHT_ACCOUNTING = "group_roles_by_st_dev_once_v1"
DISK_PREFLIGHT_FORMULA = (
    "per_filesystem_required_free_bytes=max(additional_hard_minimum_free_bytes,"
    "max(minimum_reserve_bytes,minimum_reserve_restart_multiples*"
    "source_restart_size_bytes)+sum(role_contribution_bytes))"
)
DEFAULT_PROOF_TIMEOUT_SECONDS = 120.0
DEFAULT_PROOF_POLL_SECONDS = 0.25
DEFAULT_QUIESCENCE_TIMEOUT_SECONDS = 10.0
DEFAULT_QUIESCENCE_POLL_SECONDS = 0.1
GPU_QUERY = (
    "index,uuid,pci.bus_id,ecc.errors.uncorrected.volatile.total,"
    "ecc.errors.corrected.volatile.total,memory.total,memory.used"
)
GPU_APP_QUERY = "pid,gpu_uuid"
GPU_VISIBILITY_ENVIRONMENT_KEYS = ("CUDA_VISIBLE_DEVICES", "KOKKOS_VISIBLE_DEVICES")
KOKKOS_GPU_MAP_TOKEN = "--kokkos-map-device-id-by=mpi_rank"
INPUT_TRANSPORT_KIND = "linux_proc_holder_fd_v1"
SOURCE_RESTART_FD = 200
TRAJECTORY_FD = 201
STATE_DIRECTORY_FD = 202
EVIDENCE_DIRECTORY_FD = 203
MPI_LAUNCHER_FD = 204
BINARY_EXECUTABLE_FD = 205
STAGING_DIRECTORY_PREFLIGHT_FD = 206
HOLDER_PID_TOKEN = "{holder_pid}"
DIRECTORY_TRANSPORT_KIND = "linux_proc_holder_dirfd_v1"
EXECUTABLE_TRANSPORT_KIND = "linux_proc_holder_execfd_v1"
LAUNCH_ENVIRONMENT_KIND = "explicit_values_with_rank_projection_v3"
RANK_ENVIRONMENT_PROJECTION_KIND = \
    "prrte_openmpi_pmix_single_node_projection_v2"
MCA_CONFIGURATION_KIND = "openmpi_prrte_pmix_default_files_v1"
LAUNCH_ENVIRONMENT_KEYS = (
    "HOME", "LANG", "LC_ALL", "CUDA_DEVICE_ORDER",
    "PRTE_MCA_schizo_proxy",
)
RANK_INHERITED_LAUNCH_ENVIRONMENT_KEYS = (
    "HOME", "LANG", "LC_ALL", "CUDA_DEVICE_ORDER",
)
RANK_CONSUMED_LAUNCH_ENVIRONMENT_KEYS = ("PRTE_MCA_schizo_proxy",)
RANK_ENVIRONMENT_KEYS = (
    "HOME", "LANG", "LC_ALL", "CUDA_DEVICE_ORDER",
    "OMPI_ARGV", "OMPI_COMMAND",
    "OMPI_COMM_WORLD_LOCAL_RANK", "OMPI_COMM_WORLD_LOCAL_SIZE",
    "OMPI_COMM_WORLD_NODE_RANK", "OMPI_COMM_WORLD_RANK",
    "OMPI_COMM_WORLD_SIZE", "OMPI_FILE_LOCATION",
    "OMPI_MCA_cpu_type", "OMPI_MCA_initial_wdir", "OMPI_MCA_num_procs",
    "OMPI_NUM_APP_CTX", "OMPI_UNIVERSE_SIZE", "OMPI_WORLD_LOCAL_SIZE",
    "OMPI_WORLD_SIZE", "PMIX_BFROP_BUFFER_TYPE", "PMIX_GDS_MODULE",
    "PMIX_HOSTNAME", "PMIX_NAMESPACE", "PMIX_PARAM_FILE_PASSED",
    "PMIX_RANK", "PMIX_SECURITY_MODE", "PMIX_SERVER_TMPDIR",
    "PMIX_SERVER_URI21", "PMIX_SERVER_URI2", "PMIX_SERVER_URI3",
    "PMIX_SERVER_URI41", "PMIX_SERVER_URI4", "PMIX_SYSTEM_TMPDIR",
    "PMIX_VERSION", "PRTE_LAUNCHED", "PRTE_SHARED_FS",
    "OPAL_USER_PARAMS_GIVEN", "PWD", "ZES_ENABLE_SYSMAN",
)
RANK_FIXED_ENVIRONMENT_VALUES = {
    "OMPI_MCA_cpu_type": "x86_64", "OMPI_NUM_APP_CTX": "1",
    "OMPI_UNIVERSE_SIZE": "32",
    "PMIX_BFROP_BUFFER_TYPE": "PMIX_BFROP_BUFFER_NON_DESC",
    "PMIX_GDS_MODULE": "shmem2,hash", "PMIX_PARAM_FILE_PASSED": "1",
    "PMIX_SECURITY_MODE": "native", "PMIX_SYSTEM_TMPDIR": "/tmp",
    "PMIX_VERSION": "5.0.9a1", "PRTE_LAUNCHED": "1",
    "PRTE_SHARED_FS": "FALSE", "OPAL_USER_PARAMS_GIVEN": "1",
    "ZES_ENABLE_SYSMAN": "1",
}
RANK_DERIVED_ENVIRONMENT_VALUES = {
    "OMPI_COMMAND": "athena_argv[0]",
    "OMPI_ARGV": "space_join(athena_argv[1:])",
    "OMPI_COMM_WORLD_RANK": "global_rank",
    "OMPI_COMM_WORLD_SIZE": "world_size",
    "OMPI_COMM_WORLD_LOCAL_RANK": "local_rank",
    "OMPI_COMM_WORLD_LOCAL_SIZE": "world_size_single_node",
    "OMPI_COMM_WORLD_NODE_RANK": "local_rank_single_node",
    "OMPI_FILE_LOCATION": "/tmp/ompi.<launcher_pid>/1/<global_rank>",
    "OMPI_MCA_initial_wdir": "state_dir", "OMPI_MCA_num_procs": "world_size",
    "OMPI_WORLD_LOCAL_SIZE": "world_size_single_node",
    "OMPI_WORLD_SIZE": "world_size", "PMIX_HOSTNAME": "launcher_hostname",
    "PMIX_NAMESPACE": "prterun-<hostname>-<launcher_pid>@1",
    "PMIX_RANK": "global_rank",
    "PMIX_SERVER_TMPDIR": "/tmp/ompi.<launcher_pid>",
    "PMIX_SERVER_URI21": "shared_namespace_tcp4_loopback_uri",
    "PMIX_SERVER_URI2": "shared_namespace_tcp4_loopback_uri",
    "PMIX_SERVER_URI3": "shared_namespace_tcp4_loopback_uri",
    "PMIX_SERVER_URI41": "shared_namespace_tcp4_loopback_uri",
    "PMIX_SERVER_URI4": "shared_namespace_tcp4_loopback_uri",
    "PWD": "state_dir",
}
RANK_URI_ENVIRONMENT_KEYS = (
    "PMIX_SERVER_URI21", "PMIX_SERVER_URI2", "PMIX_SERVER_URI3",
    "PMIX_SERVER_URI41", "PMIX_SERVER_URI4",
)
MCA_CONFIGURATION_LAYOUT = (
    ("home", "openmpi", ".openmpi/mca-params.conf"),
    ("home", "prte", ".prte/mca-params.conf"),
    ("home", "pmix", ".pmix/mca-params.conf"),
    ("prefix", "openmpi", "etc/openmpi-mca-params.conf"),
    ("prefix", "prte", "etc/prte-mca-params.conf"),
    ("prefix", "pmix", "etc/pmix-mca-params.conf"),
)
EVIDENCE_NAMES = (
    "launch_record", "completion_record", "run_log", "exit_status",
    "gpu_before", "gpu_after",
)
MANAGED_TERMINATION_SIGNALS = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
FORBIDDEN_OVERRIDE = re.compile(
    r"^(?:output\d+/(?:dt|dcycle|file_number|last_time|last_write_cycle)|"
    r"job/basename|problem/trajectory_file|mesh(?:block|_refinement)?/[^=]+)="
)


class LaunchFailure(RuntimeError):
    """A condition for which launch provenance cannot be proven."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise LaunchFailure(message)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _load_strict_json(raw: bytes, path: Path) -> dict[str, Any]:
    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise LaunchFailure(f"{path}: duplicate JSON key {key!r}")
            result[key] = value
        return result

    def no_constants(value: str) -> None:
        raise LaunchFailure(f"{path}: non-finite JSON number {value}")

    try:
        value = json.loads(
            raw.decode("utf-8"), object_pairs_hook=no_duplicates,
            parse_constant=no_constants,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise LaunchFailure(f"{path}: invalid JSON: {exc}") from exc
    _require(isinstance(value, dict), f"{path}: JSON root must be an object")
    return value


def _read_immutable_plan(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    absolute = path.expanduser().absolute()
    try:
        before = absolute.lstat()
        resolved = absolute.resolve(strict=True)
    except OSError as exc:
        raise LaunchFailure(f"cannot stat plan {absolute}: {exc}") from exc
    _require(absolute == resolved, "plan path must not traverse symlinks")
    _require(stat.S_ISREG(before.st_mode), "plan must be a regular non-symlink file")
    _require(not stat.S_ISLNK(before.st_mode), "plan must not be a symlink")
    _require(not (stat.S_IMODE(before.st_mode) & 0o222),
             "plan must be immutable (no write bits)")
    try:
        binding = stable_sha256(absolute)
        raw = absolute.read_bytes()
        after = absolute.lstat()
    except (OSError, ValueError, RuntimeError) as exc:
        raise LaunchFailure(f"cannot read immutable plan {absolute}: {exc}") from exc
    signature = ("st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns")
    _require(all(getattr(before, name) == getattr(after, name) for name in signature),
             "plan changed while being read")
    _require(hashlib.sha256(raw).hexdigest() == binding["sha256"],
             "plan changed between stable hash and parse")
    return _load_strict_json(raw, absolute), binding


def _validate_file_record(value: Any, label: str) -> dict[str, Any]:
    _require(isinstance(value, dict), f"{label} must be an object")
    path = value.get("path")
    size = value.get("size")
    digest = value.get("sha256")
    _require(isinstance(path, str) and Path(path).is_absolute(),
             f"{label}.path must be absolute")
    _require(isinstance(size, int) and not isinstance(size, bool) and size >= 0,
             f"{label}.size must be a nonnegative integer")
    _require(isinstance(digest, str) and
             re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
             f"{label}.sha256 must be lowercase SHA-256")
    return value


def _audit_planned_file(value: Any, label: str,
                        *, executable: bool = False) -> dict[str, Any]:
    planned = _validate_file_record(value, label)
    try:
        resolved = Path(planned["path"]).expanduser().resolve(strict=True)
    except OSError as exc:
        raise LaunchFailure(f"cannot resolve {label}: {exc}") from exc
    _require(str(resolved) == planned["path"], f"{label}.path is not canonical")
    _require(resolved.is_file(), f"{label} is not a regular file")
    if executable:
        _require(os.access(resolved, os.X_OK), f"{label} is not executable")
    try:
        audited = stable_sha256(resolved)
    except (OSError, ValueError, RuntimeError) as exc:
        raise LaunchFailure(f"cannot audit {label}: {exc}") from exc
    _require(audited["size"] == planned["size"] and
             audited["sha256"] == planned["sha256"],
             f"{label} differs from immutable plan")
    return audited


def _run_text(command: Sequence[str], runtime: "Runtime",
              *, environment: Mapping[str, str] | None = None) -> str:
    try:
        result = runtime.run(
            list(command), check=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True,
            **({"env": dict(environment)} if environment is not None else {}),
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        detail = getattr(exc, "stderr", "") or str(exc)
        raise LaunchFailure(f"command failed: {list(command)!r}: {detail.strip()}") \
            from exc
    return result.stdout.strip()


def _audit_repository(value: Any, git_plan: Any,
                      runtime: "Runtime") -> dict[str, Any]:
    _require(isinstance(value, dict), "inputs.repo must be an object")
    path = value.get("path")
    commit = value.get("commit")
    _require(value.get("clean") is True, "inputs.repo.clean must be true")
    _require(isinstance(path, str) and Path(path).is_absolute(),
             "inputs.repo.path must be absolute")
    _require(isinstance(commit, str) and
             re.fullmatch(r"[0-9a-f]{40,64}", commit) is not None,
             "inputs.repo.commit is invalid")
    try:
        repo = Path(path).resolve(strict=True)
    except OSError as exc:
        raise LaunchFailure(f"cannot resolve planned repository: {exc}") from exc
    _require(repo.is_dir() and str(repo) == path,
             "planned repository is not a canonical directory")
    git_tool = _audit_planned_file(
        git_plan, "tools.git", executable=True)
    git_path = str(git_tool["path"])
    git_environment = {
        "LANG": "C", "LC_ALL": "C",
        "GIT_CONFIG_NOSYSTEM": "1", "GIT_CONFIG_GLOBAL": "/dev/null",
    }
    git_configuration = [
        "-c", "core.fsmonitor=false",
        "-c", "core.untrackedCache=false",
        "-c", "core.hooksPath=/dev/null",
    ]
    top = Path(_run_text(
        [git_path, "-C", str(repo), *git_configuration,
         "rev-parse", "--show-toplevel"], runtime,
        environment=git_environment,
    )).resolve(strict=True)
    _require(top == repo, "planned repository is not its worktree root")
    actual_commit = _run_text(
        [git_path, "-C", str(repo), *git_configuration,
         "rev-parse", "--verify", "HEAD"], runtime,
        environment=git_environment)
    _require(actual_commit == commit, "planned repository commit changed")
    status = _run_text(
        [git_path, "-C", str(repo), *git_configuration,
         "status", "--porcelain=v1",
         "--untracked-files=all"], runtime, environment=git_environment)
    _require(not status, "planned repository is no longer clean")
    git_tool_after = _audit_planned_file(
        git_plan, "tools.git after repository audit", executable=True)
    _require(_same_artifact_binding(git_tool, git_tool_after),
             "plan-bound Git executable changed during repository audit")
    stable_git_tool = {
        key: item for key, item in git_tool_after.items()
        if key != "age_seconds"
    }
    return {
        "path": str(repo), "commit": actual_commit, "clean": True,
        "git_tool": stable_git_tool,
        "git_environment_policy": "explicit_clean_environment_v1",
        "git_environment": git_environment,
        "git_configuration": git_configuration,
    }


def _integer(value: Any, label: str) -> int:
    _require(isinstance(value, int) and not isinstance(value, bool),
             f"{label} must be an integer")
    return value


def _ceil_ratio(numerator: int, denominator: int) -> int:
    return (numerator + denominator - 1) // denominator


def _capacity_transition_contract(plan: Mapping[str, Any]) -> dict[str, Any]:
    value = plan.get("capacity_transition")
    _require(isinstance(value, dict) and set(value) == {
        "kind", "parameter", "source_max_nmb_per_rank",
        "target_max_nmb_per_rank", "maximum_target_max_nmb_per_rank",
        "runtime_override", "gpu_memory_model",
    }, "capacity_transition must be the exact strict object")
    source = _integer(value.get("source_max_nmb_per_rank"),
                      "capacity source")
    target = _integer(value.get("target_max_nmb_per_rank"),
                      "capacity target")
    _require(1 <= source <= target <= MAX_NMB_PER_RANK and
             value.get("maximum_target_max_nmb_per_rank") == MAX_NMB_PER_RANK and
             value.get("parameter") ==
             "mesh_refinement/max_nmb_per_rank",
             "capacity transition bounds are invalid")
    runtime_override = value.get("runtime_override")
    if source == target:
        _require(value.get("kind") == "unchanged_v1" and
                 runtime_override is None,
                 "unchanged capacity transition may not carry an override")
    else:
        _require(value.get("kind") == "increase_v1" and
                 runtime_override ==
                 f"mesh_refinement/max_nmb_per_rank={target}",
                 "capacity increase lacks the canonical override")
    required = _ceil_ratio(
        target * CAPACITY_MEMORY_MIB_PER_SLOT_NUMERATOR,
        CAPACITY_MEMORY_MIB_PER_SLOT_DENOMINATOR)
    minimum_total = _ceil_ratio(
        target * CAPACITY_MEMORY_MIB_PER_SLOT_NUMERATOR *
        CAPACITY_GPU_USABLE_FRACTION_DENOMINATOR,
        CAPACITY_MEMORY_MIB_PER_SLOT_DENOMINATOR *
        CAPACITY_GPU_USABLE_FRACTION_NUMERATOR)
    _require(value.get("gpu_memory_model") == {
        "kind": "fixed_conservative_per_meshblock_slot_v1",
        "mib_per_slot_numerator":
            CAPACITY_MEMORY_MIB_PER_SLOT_NUMERATOR,
        "mib_per_slot_denominator":
            CAPACITY_MEMORY_MIB_PER_SLOT_DENOMINATOR,
        "usable_fraction_numerator":
            CAPACITY_GPU_USABLE_FRACTION_NUMERATOR,
        "usable_fraction_denominator":
            CAPACITY_GPU_USABLE_FRACTION_DENOMINATOR,
        "required_per_rank_memory_mib_ceiling": required,
        "minimum_gpu_memory_total_mib": minimum_total,
    }, "capacity GPU-memory model is not canonical")
    source_record = plan.get("source")
    _require(isinstance(source_record, Mapping),
             "plan source must be an object for capacity transition")
    source_parameters = source_record.get("parameters")
    _require(isinstance(source_parameters, dict) and
             source_parameters.get("mesh_refinement", {}).get(
                 "max_nmb_per_rank") == str(source),
             "source restart capacity differs from capacity transition")
    return value


def _finite(value: Any, label: str) -> float:
    _require(isinstance(value, (int, float)) and not isinstance(value, bool),
             f"{label} must be a number")
    number = float(value)
    _require(math.isfinite(number), f"{label} must be finite")
    return number


def _restart_cadence_transition_contract(
        plan: Mapping[str, Any], root_dt: float) -> dict[str, Any]:
    value = plan.get("restart_cadence_transition")
    _require(isinstance(value, dict) and set(value) == {
        "kind", "block", "parameter", "source_dt", "target_dt", "root_dt",
        "source_root_step_multiple", "target_root_step_multiple", "phase",
        "runtime_override",
    }, "restart_cadence_transition must be the exact strict object")
    source_record = plan.get("source")
    _require(isinstance(source_record, Mapping),
             "source must bind restart cadence")
    source_parameters = source_record.get("parameters")
    _require(isinstance(source_parameters, dict),
             "source parameters must bind restart cadence")
    restart_blocks = sorted(
        name for name, block in source_parameters.items()
        if isinstance(name, str) and re.fullmatch(r"output\d+", name) and
        isinstance(block, dict) and
        block.get("file_type", "").strip() == "rst")
    _require(len(restart_blocks) == 1,
             "source must contain exactly one restart output")
    block_name = restart_blocks[0]
    block = source_parameters[block_name]
    _require("dt" in block and "dcycle" not in block,
             f"source {block_name} must use dt without dcycle")
    try:
        source_dt = float(block["dt"])
        phase = {
            "file_number": int(block.get("file_number", "0")),
            "last_time": float(block.get("last_time", "-1")),
            "last_write_cycle": int(block.get("last_write_cycle", "-1")),
        }
    except (TypeError, ValueError) as exc:
        raise LaunchFailure("source restart cadence state is invalid") from exc
    planned_phase = value.get("phase")
    _require(type(value.get("source_dt")) is float and
             type(value.get("target_dt")) is float and
             type(value.get("root_dt")) is float and
             type(value.get("source_root_step_multiple")) is int and
             type(value.get("target_root_step_multiple")) is int and
             isinstance(planned_phase, dict) and
             set(planned_phase) == {
                 "file_number", "last_time", "last_write_cycle"} and
             type(planned_phase.get("file_number")) is int and
             type(planned_phase.get("last_time")) is float and
             type(planned_phase.get("last_write_cycle")) is int,
             "restart cadence transition numeric types are not canonical")
    target_dt = _finite(value.get("target_dt"), "restart target dt")
    planned_source_dt = _finite(value.get("source_dt"), "restart source dt")
    planned_root_dt = _finite(value.get("root_dt"), "restart root dt")
    source_ratio = planned_source_dt / root_dt
    target_ratio = target_dt / root_dt
    source_nearest = round(source_ratio) if math.isfinite(source_ratio) else 0
    target_nearest = round(target_ratio) if math.isfinite(target_ratio) else 0
    source_tolerance = (CADENCE_MULTIPLE_MAX_ULPS * max(
        math.ulp(source_ratio), math.ulp(float(source_nearest)))) \
        if math.isfinite(source_ratio) else 0.0
    target_tolerance = (CADENCE_MULTIPLE_MAX_ULPS * max(
        math.ulp(target_ratio), math.ulp(float(target_nearest)))) \
        if math.isfinite(target_ratio) else 0.0
    _require(source_dt == planned_source_dt and planned_source_dt > 0.0 and
             target_dt > 0.0 and target_dt <= planned_source_dt and
             planned_root_dt == root_dt and math.isfinite(source_ratio) and
             math.isfinite(target_ratio) and source_nearest >= 1 and
             target_nearest >= 1 and
             abs(source_ratio - source_nearest) <= source_tolerance and
             abs(target_ratio - target_nearest) <= target_tolerance and
             phase["file_number"] >= 0 and
             math.isfinite(phase["last_time"]) and
             value.get("block") == block_name and
             value.get("parameter") == f"{block_name}/dt" and
             value.get("source_root_step_multiple") == source_nearest and
             value.get("target_root_step_multiple") == target_nearest and
             planned_phase == phase,
             "restart cadence transition source/target/phase is not canonical")
    if target_dt == planned_source_dt:
        _require(value.get("kind") == "unchanged_v1" and
                 value.get("runtime_override") is None,
                 "unchanged restart cadence may not carry an override")
    else:
        _require(value.get("kind") == "tighten_v1" and
                 value.get("runtime_override") ==
                 f"{block_name}/dt={target_dt!r}",
                 "restart cadence tightening lacks its canonical override")
    return value


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _validate_disk_preflight_contract(value: Any,
                                      source_restart_size: int,
                                      trajectory_size: int) -> dict[str, Any]:
    _require(isinstance(value, dict),
             "launch_contract.disk_preflight must be an object")
    peak_gib = _integer(
        value.get("planned_peak_output_gib"),
        "disk_preflight.planned_peak_output_gib")
    _require(peak_gib > 0, "planned peak output GiB must be positive")
    peak_bytes = peak_gib * GIB_BYTES
    expected = {
        "kind": DISK_PREFLIGHT_KIND,
        "accounting": DISK_PREFLIGHT_ACCOUNTING,
        "formula": DISK_PREFLIGHT_FORMULA,
        "used_percent_exclusive_max": 75,
        "minimum_reserve_bytes": 50 * GIB_BYTES,
        "minimum_reserve_restart_multiples": 2,
        "additional_hard_minimum_free_bytes": 180 * GIB_BYTES,
        "source_restart_size_bytes": source_restart_size,
        "source_restart_staging_bytes": source_restart_size,
        "trajectory_staging_bytes": trajectory_size,
        "staging_copy_bytes": source_restart_size + trajectory_size,
        "planned_peak_output_gib": peak_gib,
        "planned_peak_output_bytes": peak_bytes,
        "role_contributions_bytes": {
            "state_dir": {"planned_peak_output_bytes": peak_bytes},
            "staging_dir": {
                "source_restart_staging_bytes": source_restart_size,
                "trajectory_staging_bytes": trajectory_size,
            },
        },
        "bound_directory_fds": {
            "state_dir": STATE_DIRECTORY_FD,
            "staging_dir": STAGING_DIRECTORY_PREFLIGHT_FD,
        },
    }
    _require(value == expected,
             "launch disk-preflight contract is not the exact capacity policy")
    return value


def _disk_preflight_snapshot(
        contract: Mapping[str, Any], phase: str,
        role_bindings: Mapping[str, Mapping[str, Any]],
        statvfs_reader: Callable[[Any], Any]) -> dict[str, Any]:
    """Measure and enforce one conservative budget per unique ``st_dev``."""

    _require(phase in ("before_staging", "before_spawn"),
             "disk preflight phase is invalid")
    _require(set(role_bindings) == {"state_dir", "staging_dir"},
             "disk preflight must bind state and staging directories")
    grouped: dict[int, list[str]] = {}
    for role, binding in role_bindings.items():
        device = binding.get("device")
        _require(isinstance(device, int) and not isinstance(device, bool) and device >= 0,
                 f"disk preflight {role} device is invalid")
        grouped.setdefault(device, []).append(role)

    reserve = max(
        int(contract["minimum_reserve_bytes"]),
        int(contract["minimum_reserve_restart_multiples"]) *
        int(contract["source_restart_size_bytes"]),
    )
    role_contributions = {
        "state_dir": int(contract["planned_peak_output_bytes"]),
        "staging_dir": int(contract["staging_copy_bytes"]),
    }
    filesystems: list[dict[str, Any]] = []
    for device in sorted(grouped):
        roles = sorted(grouped[device])
        contributions = {role: role_contributions[role] for role in roles}
        contribution_total = sum(contributions.values())
        required_free = max(
            int(contract["additional_hard_minimum_free_bytes"]),
            reserve + contribution_total,
        )
        observations: dict[str, dict[str, Any]] = {}
        for role in roles:
            binding = role_bindings[role]
            try:
                fs = statvfs_reader(binding["target"])
            except OSError as exc:
                raise LaunchFailure(
                    f"cannot read {phase} statvfs for {role} on device "
                    f"{device}: {exc}") from exc
            fragment = int(getattr(fs, "f_frsize", 0))
            total_blocks = int(getattr(fs, "f_blocks", -1))
            free_blocks = int(getattr(fs, "f_bfree", -1))
            available_blocks = int(getattr(fs, "f_bavail", -1))
            _require(fragment > 0 and total_blocks > 0 and
                     0 <= available_blocks <= free_blocks <= total_blocks,
                     f"disk preflight {phase} returned invalid statvfs counts")
            free_bytes = free_blocks * fragment
            available_bytes = available_blocks * fragment
            effective_free = min(free_bytes, available_bytes)
            used_blocks = total_blocks - free_blocks
            used_numerator = used_blocks * 100
            used_denominator = used_blocks + available_blocks
            _require(used_denominator > 0,
                     f"disk preflight {phase} has no usable filesystem blocks")
            used_pass = used_numerator < (
                used_denominator * int(contract["used_percent_exclusive_max"]))
            free_pass = effective_free >= required_free
            _require(used_pass,
                     f"disk preflight {phase} rejects {role} on device {device}: "
                     f"used space is at or above "
                     f"{contract['used_percent_exclusive_max']}%")
            _require(free_pass,
                     f"disk preflight {phase} rejects {role} on device {device}: "
                     f"effective free {effective_free} bytes is below required "
                     f"{required_free} bytes")
            observations[role] = {
                "access": dict(binding["access"]),
                "statvfs": {
                    "fragment_size_bytes": fragment,
                    "total_blocks": total_blocks,
                    "free_blocks": free_blocks,
                    "available_blocks": available_blocks,
                    "total_bytes": total_blocks * fragment,
                    "free_bytes": free_bytes,
                    "available_bytes": available_bytes,
                    "effective_free_bytes": effective_free,
                    "used_bytes": used_blocks * fragment,
                    "used_percent_numerator": used_numerator,
                    "used_percent_denominator": used_denominator,
                },
                "used_percent_pass": True,
                "free_bytes_pass": True,
            }
        filesystems.append({
            "device": device,
            "roles": roles,
            "observations": observations,
            "role_contribution_bytes": contributions,
            "contribution_bytes_total": contribution_total,
            "reserve_bytes": reserve,
            "required_free_bytes": required_free,
            "status": "pass",
        })
    return {
        "kind": DISK_PREFLIGHT_KIND,
        "phase": phase,
        "created_utc": _utc_now(),
        "contract_sha256": _canonical_sha256(contract),
        "filesystems": filesystems,
        "status": "pass",
    }


def _athena_float32(value: float) -> float:
    try:
        return struct.unpack("=f", struct.pack("=f", value))[0]
    except (OverflowError, struct.error) as exc:
        raise LaunchFailure(
            f"output cadence value is not representable as float32: {value}") from exc


def _production_validate_plan(plan: dict[str, Any]) -> dict[str, Any]:
    try:
        validated = validate_strict_plan(plan)
    except (PlanCheckFailure, OSError, ValueError, RuntimeError) as exc:
        raise LaunchFailure(f"strict plan validation failed: {exc}") from exc
    _require(isinstance(validated, dict),
             "strict plan validator did not return a validation record")
    return validated


def _production_validate_source(path: Path,
                                plan: Mapping[str, Any]) -> dict[str, Any]:
    """Independently replay the source contract before any GPU process starts."""

    try:
        actual_audit = audit_restart(path)
        metadata = read_restart_metadata(path)
    except (OSError, ValueError, RuntimeError) as exc:
        raise LaunchFailure(f"live source restart audit failed: {exc}") from exc
    inputs = plan["inputs"]
    planned_audit = plan["source_qualification"]["audit"]
    _require(actual_audit.get("valid") is True and
             actual_audit.get("sha256") == inputs["source_restart"]["sha256"] and
             actual_audit.get("sha256") == planned_audit.get("sha256") and
             all(actual_audit.get(name) == planned_audit.get(name) for name in (
                 "metadata", "layout", "topology", "stored_reals",
             )), "live source restart differs from its planned full audit")
    output3 = metadata.parameters.get("output3")
    _require(isinstance(output3, dict) and "dcycle" not in output3 and
             "dt" in output3,
             "live source output3 must use dt and must not serialize dcycle")
    try:
        source_dt = float(output3["dt"])
        last_time = float(output3.get("last_time", "-1"))
        file_number = int(output3.get("file_number", "0"))
    except (TypeError, ValueError) as exc:
        raise LaunchFailure("live source output3 cadence state is invalid") from exc
    expected = plan["expected"]
    root_steps = int(expected["root_steps"])
    _require(file_number + root_steps < 100000,
             "planned output3 file_number reaches the C++ five-digit limit")
    try:
        capacity = int(metadata.parameters["mesh_refinement"]["max_nmb_per_rank"])
        ranks = int(plan["policy"]["ranks"])
    except (KeyError, TypeError, ValueError) as exc:
        raise LaunchFailure(
            "live source restart lacks canonical rank-capacity parameters") from exc
    source_summary = restart_metadata_dict(metadata, ranks, capacity)
    source_summary["parameters"] = metadata.parameters
    parameters_sha256 = _canonical_sha256(metadata.parameters)
    source_summary["parameters_sha256"] = parameters_sha256
    _require(source_summary == plan.get("source"),
             "live source summary differs from the immutable plan")
    _require(plan.get("parameter_snapshots", {}).get("source_sha256") ==
             parameters_sha256,
             "live source parameter hash differs from the planned snapshot")

    root_dt = float(expected["root_dt"])
    tlim = float(expected["tlim"])
    _require(math.isfinite(source_dt) and source_dt >= 0.0 and
             math.isfinite(last_time) and file_number >= 0,
             "live source output3 cadence state is non-finite or negative")
    if plan["source_qualification"]["mode"] == "parent_segment_pass":
        _require(source_dt == root_dt,
                 "parent-chain source output3/dt differs from root_dt")
    current_time = float(expected["source_time"])
    writes: list[dict[str, Any]] = []
    for cycle in range(int(expected["source_cycle"]) + 1,
                       int(expected["final_cycle"]) + 1):
        current_time += root_dt
        due = (_athena_float32(current_time) >=
               _athena_float32(last_time + root_dt) and
               _athena_float32(current_time) < _athena_float32(tlim))
        _require(due,
                 f"live source output3 cadence is not due at root cycle {cycle}")
        writes.append({
            "cycle": cycle, "time": current_time,
            "kind": "scheduled", "file_number": file_number,
        })
        file_number += 1
        last_time = current_time if last_time < 0.0 else last_time + root_dt
    planned_output3 = next(
        (row for row in plan["outputs"] if row.get("block") == "output3"), None)
    runtime_output_snapshots = {
        row["block"]: row["parameters"] for row in plan["outputs"]
    }
    _require(plan.get("parameter_snapshots", {}).get("output_blocks") ==
             runtime_output_snapshots,
             "planned output parameter snapshots are not exact")
    _require(isinstance(planned_output3, dict) and
             planned_output3.get("expected_writes") == writes and
             planned_output3.get("expected_endpoint_state") == {
                 "file_number": file_number,
                 "last_time": last_time,
                 "last_write_cycle": int(expected["final_cycle"]),
             }, "planned output3 cadence differs from live-source float32 replay")
    return {
        "source_restart_sha256": actual_audit["sha256"],
        "source_parameters_sha256": parameters_sha256,
        "full_audit": actual_audit,
        "output3_root_cycle_replay": {
            "first_cycle": writes[0]["cycle"],
            "last_cycle": writes[-1]["cycle"],
            "write_count": len(writes),
            "all_scheduled": True,
            "cadence": root_dt,
        },
    }


def _wall_time_token(seconds: int) -> str:
    _require(seconds > 0, "launch wall time must be positive")
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds_part = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds_part:02d}"


def _safe_directory(path: Path, label: str) -> Path:
    absolute = path.expanduser().absolute()
    try:
        resolved = absolute.resolve(strict=True)
        info = absolute.lstat()
    except OSError as exc:
        raise LaunchFailure(f"cannot resolve {label} {absolute}: {exc}") from exc
    _require(absolute == resolved, f"{label} must not traverse symlinks")
    _require(stat.S_ISDIR(info.st_mode), f"{label} must be a directory")
    return resolved


def _validate_launch_environment(value: Any) -> dict[str, str]:
    """Return the exact, immutable environment permitted for MPI/AthenaK."""

    _require(isinstance(value, dict) and
             set(value) == {"kind", "values", "sha256", "rank_projection"} and
             value.get("kind") == LAUNCH_ENVIRONMENT_KIND,
             f"launch environment kind must be {LAUNCH_ENVIRONMENT_KIND}")
    values = value.get("values")
    _require(isinstance(values, dict) and set(values) == set(LAUNCH_ENVIRONMENT_KEYS),
             "launch environment must contain exactly HOME/LANG/LC_ALL/"
             "CUDA_DEVICE_ORDER/PRTE_MCA_schizo_proxy")
    _require(all(isinstance(key, str) and isinstance(item, str) and
                 key and "=" not in key and "\x00" not in key and
                 "\x00" not in item for key, item in values.items()),
             "launch environment contains an invalid name or value")
    _require(values["LANG"] == "C" and values["LC_ALL"] == "C" and
             values["CUDA_DEVICE_ORDER"] == "PCI_BUS_ID" and
             values["PRTE_MCA_schizo_proxy"] == "ompi",
             "launch environment locale/CUDA/PRRTE personality is not canonical")
    home = _safe_directory(Path(values["HOME"]), "launch HOME")
    home_info = home.stat()
    _require(str(home) == values["HOME"] and home_info.st_uid == os.geteuid() and
             not (stat.S_IMODE(home_info.st_mode) & 0o022),
             "launch HOME must be canonical, launcher-owned, and not group/other "
             "writable")
    canonical = {key: values[key] for key in LAUNCH_ENVIRONMENT_KEYS}
    _require(value.get("sha256") == _canonical_sha256(canonical),
             "launch environment SHA-256 is not canonical")
    projection = value.get("rank_projection")
    _require(isinstance(projection, dict) and
             set(projection) == {
                 "kind", "inherited_values", "consumed_absent", "exact_keys",
                 "fixed_values", "derived_values", "sha256",
             } and
             projection.get("kind") == RANK_ENVIRONMENT_PROJECTION_KIND,
             "rank environment projection kind is not canonical")
    projection_payload = {
        "inherited_values": {
            key: canonical[key] for key in RANK_INHERITED_LAUNCH_ENVIRONMENT_KEYS
        },
        "consumed_absent": list(RANK_CONSUMED_LAUNCH_ENVIRONMENT_KEYS),
        "exact_keys": list(RANK_ENVIRONMENT_KEYS),
        "fixed_values": dict(RANK_FIXED_ENVIRONMENT_VALUES),
        "derived_values": dict(RANK_DERIVED_ENVIRONMENT_VALUES),
    }
    _require(projection.get("inherited_values") ==
             projection_payload["inherited_values"] and
             projection.get("consumed_absent") ==
             projection_payload["consumed_absent"] and
             projection.get("sha256") == _canonical_sha256(projection_payload),
             "rank environment projection differs from the exact PRRTE contract")
    return canonical


def _snapshot_mca_configuration_file(scope: str, project: str,
                                     path: Path) -> dict[str, Any]:
    """Independently snapshot one default MCA path, including absence."""

    absolute = path.absolute()
    try:
        info = absolute.lstat()
    except FileNotFoundError:
        _require(not os.path.lexists(absolute),
                 f"MCA configuration path is not safely absent: {absolute}")
        return {
            "scope": scope, "project": project,
            "path": str(absolute), "state": "absent",
        }
    except OSError as exc:
        raise LaunchFailure(f"cannot inspect MCA configuration {absolute}: {exc}") \
            from exc
    try:
        resolved = absolute.resolve(strict=True)
    except OSError as exc:
        raise LaunchFailure(f"cannot resolve MCA configuration {absolute}: {exc}") \
            from exc
    mode = stat.S_IMODE(info.st_mode)
    _require(resolved == absolute and not stat.S_ISLNK(info.st_mode) and
             stat.S_ISREG(info.st_mode) and info.st_uid in (0, os.geteuid()) and
             not (mode & 0o022),
             f"MCA configuration is not canonical, regular, safely owned/mode: "
             f"{absolute}")
    try:
        binding = stable_sha256(absolute)
    except (OSError, ValueError, RuntimeError) as exc:
        raise LaunchFailure(f"cannot audit MCA configuration {absolute}: {exc}") \
            from exc
    return {
        "scope": scope, "project": project, "path": str(absolute),
        "state": "present", "device": binding["device"],
        "inode": binding["inode"], "owner_uid": info.st_uid,
        "mode": f"{mode:04o}", "size": binding["size"],
        "mtime_ns": binding["mtime_ns"], "ctime_ns": binding["ctime_ns"],
        "sha256": binding["sha256"],
        "closure_check": binding["closure_check"],
    }


def _snapshot_mca_prefix_directory(path: Path) -> dict[str, Any]:
    absolute = path.absolute()
    try:
        resolved = absolute.resolve(strict=True)
        info = absolute.lstat()
    except OSError as exc:
        raise LaunchFailure(f"cannot bind MCA prefix {absolute}: {exc}") from exc
    mode = stat.S_IMODE(info.st_mode)
    _require(resolved == absolute and not stat.S_ISLNK(info.st_mode) and
             stat.S_ISDIR(info.st_mode) and info.st_uid in (0, os.geteuid()) and
             not (mode & 0o022),
             "MCA prefix must be canonical, non-symlink, safely owned/mode")
    return {
        "path": str(absolute), "device": info.st_dev, "inode": info.st_ino,
        "owner_uid": info.st_uid, "mode": f"{mode:04o}",
    }


def _audit_mca_configuration(value: Any, home: str) -> dict[str, Any]:
    """Validate the plan binding and return a fresh six-path snapshot."""

    _require(isinstance(value, dict) and
             set(value) == {"kind", "home", "prefix", "prefix_directory",
                            "files", "sha256"} and
             value.get("kind") == MCA_CONFIGURATION_KIND,
             f"MCA configuration contract must be {MCA_CONFIGURATION_KIND}")
    prefix_value = value.get("prefix")
    _require(isinstance(prefix_value, str) and Path(prefix_value).is_absolute(),
             "MCA configuration prefix must be absolute")
    prefix_directory = _snapshot_mca_prefix_directory(Path(prefix_value))
    prefix = Path(prefix_directory["path"])
    _require(value.get("home") == home and
             value.get("prefix_directory") == prefix_directory,
             "MCA configuration HOME/prefix directory differs from plan")
    planned_files = value.get("files")
    _require(isinstance(planned_files, list) and
             len(planned_files) == len(MCA_CONFIGURATION_LAYOUT),
             "MCA configuration file set is not exact")
    current_files: list[dict[str, Any]] = []
    for index, (scope, project, relative) in enumerate(MCA_CONFIGURATION_LAYOUT):
        expected_path = (Path(home) if scope == "home" else prefix) / relative
        current = _snapshot_mca_configuration_file(
            scope, project, expected_path)
        _require(planned_files[index] == current,
                 f"MCA configuration {expected_path} differs from immutable plan")
        current_files.append(current)
    payload = {
        "kind": MCA_CONFIGURATION_KIND,
        "home": home, "prefix": str(prefix),
        "prefix_directory": prefix_directory, "files": current_files,
    }
    _require(value.get("sha256") == _canonical_sha256(payload),
             "MCA configuration contract SHA-256 is not canonical")
    return {**payload, "sha256": value["sha256"]}


def _directory_signature(info: os.stat_result) -> tuple[int, int, int, int]:
    return (info.st_dev, info.st_ino, info.st_uid, stat.S_IMODE(info.st_mode))


def _open_bound_directory(path: Path, label: str) -> tuple[Path, int, dict[str, Any]]:
    resolved = _safe_directory(path, label)
    before = resolved.lstat()
    _require(before.st_uid == os.geteuid(), f"{label} must be owned by launcher")
    _require(not (stat.S_IMODE(before.st_mode) & 0o022),
             f"{label} must not be group/other writable")
    try:
        descriptor = os.open(
            resolved, os.O_RDONLY | os.O_DIRECTORY |
            getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
    except OSError as exc:
        raise LaunchFailure(f"cannot bind {label}: {exc}") from exc
    try:
        bound = os.fstat(descriptor)
        _require(_directory_signature(bound) == _directory_signature(before),
                 f"{label} changed while being bound")
    except BaseException:
        os.close(descriptor)
        raise
    return resolved, descriptor, {
        "path": str(resolved), "device": bound.st_dev, "inode": bound.st_ino,
        "owner_uid": bound.st_uid, "mode": f"{stat.S_IMODE(bound.st_mode):04o}",
    }


def _validate_new_target(path: Path, label: str) -> Path:
    absolute = path.expanduser().absolute()
    bound_evidence_parent = Path(
        f"/proc/{os.getpid()}/fd/{EVIDENCE_DIRECTORY_FD}")
    if absolute.parent == bound_evidence_parent:
        _require(_fd_is_open(EVIDENCE_DIRECTORY_FD) and
                 stat.S_ISDIR(os.fstat(EVIDENCE_DIRECTORY_FD).st_mode),
                 "bound evidence directory descriptor is unavailable")
    else:
        _safe_directory(absolute.parent, f"{label} parent")
    _require(not os.path.lexists(absolute), f"refusing to overwrite {label}: {absolute}")
    return absolute


def _paths_separate(left: Path, right: Path, message: str) -> None:
    for first, second in ((left, right), (right, left)):
        try:
            first.relative_to(second)
        except ValueError:
            continue
        raise LaunchFailure(message)


def _audit_open_descriptor(descriptor: int, planned: Mapping[str, Any],
                           label: str) -> dict[str, Any]:
    """Hash a fixed read-only descriptor without consulting its pathname."""

    before = os.fstat(descriptor)
    _require(stat.S_ISREG(before.st_mode), f"{label} descriptor is not regular")
    digest = hashlib.sha256()
    offset = 0
    while offset < before.st_size:
        chunk = os.pread(descriptor, min(HASH_CHUNK_BYTES, before.st_size - offset), offset)
        _require(bool(chunk), f"short read while hashing {label}")
        digest.update(chunk)
        offset += len(chunk)
    after = os.fstat(descriptor)
    signature = (before.st_dev, before.st_ino, before.st_size,
                 before.st_mtime_ns, before.st_ctime_ns)
    _require(signature == (
        after.st_dev, after.st_ino, after.st_size,
        after.st_mtime_ns, after.st_ctime_ns,
    ), f"{label} descriptor changed while hashing")
    _require(before.st_size == planned["size"] and
             digest.hexdigest() == planned["sha256"],
             f"{label} descriptor differs from immutable plan")
    return {
        "path": planned["path"],
        "device": before.st_dev,
        "inode": before.st_ino,
        "size": before.st_size,
        "mtime_ns": before.st_mtime_ns,
        "ctime_ns": before.st_ctime_ns,
        "sha256": digest.hexdigest(),
        "closure_check": "fixed_read_only_descriptor",
    }


def _staging_directory(path: Path, state_dir: Path) -> tuple[Path, int]:
    staging = _safe_directory(path, "staging directory")
    _paths_separate(staging, state_dir,
                    "staging and state directories must not contain each other")
    info = staging.stat()
    _require(info.st_uid == os.geteuid(), "staging directory must be owned by launcher")
    _require(not (stat.S_IMODE(info.st_mode) & 0o022),
             "staging directory must not be group/other writable")
    _require(stat.S_IMODE(info.st_mode) & 0o200,
             "empty staging directory must initially be owner-writable")
    _require(not any(staging.iterdir()), "staging directory must initially be empty")
    _require(not _fd_is_open(STAGING_DIRECTORY_PREFLIGHT_FD),
             f"refusing to replace existing process fd "
             f"{STAGING_DIRECTORY_PREFLIGHT_FD}")
    descriptor: int | None = None
    try:
        descriptor = os.open(
            staging,
            os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NOFOLLOW", 0),
        )
    except OSError as exc:
        raise LaunchFailure(f"cannot bind staging directory descriptor: {exc}") from exc
    try:
        bound = os.fstat(descriptor)
        _require(bound.st_dev == info.st_dev and bound.st_ino == info.st_ino,
                 "staging directory changed while being opened")
        if descriptor != STAGING_DIRECTORY_PREFLIGHT_FD:
            os.dup2(
                descriptor, STAGING_DIRECTORY_PREFLIGHT_FD, inheritable=False)
            os.close(descriptor)
            descriptor = None
        else:
            os.set_inheritable(descriptor, False)
        fixed = os.fstat(STAGING_DIRECTORY_PREFLIGHT_FD)
        _require(fixed.st_dev == info.st_dev and fixed.st_ino == info.st_ino,
                 "fixed staging directory descriptor changed while being bound")
    except BaseException:
        if descriptor is not None:
            os.close(descriptor)
        if _fd_is_open(STAGING_DIRECTORY_PREFLIGHT_FD):
            os.close(STAGING_DIRECTORY_PREFLIGHT_FD)
        raise
    return staging, STAGING_DIRECTORY_PREFLIGHT_FD


def _stage_from_audited_descriptor(source: Mapping[str, Any],
                                   staged_plan: Any,
                                   staging_dir: Path,
                                   staging_directory_fd: int,
                                   label: str) -> tuple[dict[str, Any], int]:
    """Copy the exact descriptor being hashed into a no-replace staged file."""

    planned = _validate_file_record(staged_plan, f"launch_contract.{label}")
    destination = Path(planned["path"]).expanduser().absolute()
    _require(destination.parent == staging_dir and destination.name not in ("", ".", ".."),
             f"{label} must be a direct child of staging directory")
    _require(planned["size"] == source["size"] and
             planned["sha256"] == source["sha256"],
             f"{label} content identity must equal its original planned input")
    try:
        os.stat(destination.name, dir_fd=staging_directory_fd,
                follow_symlinks=False)
    except FileNotFoundError:
        pass
    else:
        raise LaunchFailure(f"refusing pre-existing staged target: {destination}")

    temporary = staging_dir / (
        f".{destination.name}.{os.getpid()}.{time.time_ns()}.stage"
    )
    descriptor: int | None = None
    staged_fd: int | None = None
    source_stream: Any | None = None
    destination_signature: tuple[int, int, int, int] | None = None
    destination_installed = False
    committed = False
    digest = hashlib.sha256()
    try:
        try:
            try:
                checked_path, source_stream, source_signature = \
                    _open_regular_nofollow(source["path"])
                _require(source_signature["size"] == source["size"],
                         f"original {label} size differs from plan")
                exempt = {(os.getpid(), source_stream.fileno())}
                _assert_closed(checked_path, source_signature, exempt)
                descriptor = os.open(
                    temporary.name, os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                    getattr(os, "O_CLOEXEC", 0) |
                    getattr(os, "O_NOFOLLOW", 0), 0o600,
                    dir_fd=staging_directory_fd,
                )
                with os.fdopen(descriptor, "wb", closefd=True) as destination_stream:
                    descriptor = None
                    for chunk in iter(
                            lambda: source_stream.read(HASH_CHUNK_BYTES), b""):
                        digest.update(chunk)
                        destination_stream.write(chunk)
                    destination_stream.flush()
                    os.fsync(destination_stream.fileno())
                    os.fchmod(destination_stream.fileno(), 0o444)
                    os.fsync(destination_stream.fileno())
                _assert_stream_signature(
                    source_stream, checked_path, source_signature, "while staging")
                _assert_path_signature(
                    checked_path, source_signature, "while staging")
                _assert_closed(checked_path, source_signature, exempt)
                source_stream.close()
                source_stream = None
                _assert_path_signature(checked_path, source_signature, "after staging")
                _assert_closed(checked_path, source_signature)
                _require(digest.hexdigest() == source["sha256"],
                         f"original {label} bytes differ from immutable plan")
                temporary_info = os.stat(
                    temporary.name, dir_fd=staging_directory_fd,
                    follow_symlinks=False)
                _require(stat.S_ISREG(temporary_info.st_mode),
                         f"temporary staged {label} is not regular")
                destination_signature = (
                    temporary_info.st_dev, temporary_info.st_ino,
                    temporary_info.st_size, temporary_info.st_mtime_ns,
                )
                try:
                    os.link(
                        temporary.name, destination.name,
                        src_dir_fd=staging_directory_fd,
                        dst_dir_fd=staging_directory_fd,
                        follow_symlinks=False,
                    )
                except FileExistsError as exc:
                    raise LaunchFailure(
                        f"refusing pre-existing staged target: {destination}") from exc
                destination_installed = True
                linked = os.stat(
                    destination.name, dir_fd=staging_directory_fd,
                    follow_symlinks=False)
                _require(stat.S_ISREG(linked.st_mode) and
                         (linked.st_dev, linked.st_ino) ==
                         (temporary_info.st_dev, temporary_info.st_ino),
                         f"new staged {label} link does not bind the copied inode")
                staged_fd = os.open(
                    destination.name,
                    os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                    getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=staging_directory_fd,
                )
            finally:
                if source_stream is not None:
                    source_stream.close()
                if descriptor is not None:
                    os.close(descriptor)
                try:
                    os.unlink(temporary.name, dir_fd=staging_directory_fd)
                except FileNotFoundError:
                    pass
                os.fsync(staging_directory_fd)
        except (OSError, ValueError, RuntimeError) as exc:
            if isinstance(exc, LaunchFailure):
                raise
            raise LaunchFailure(f"cannot stage {label}: {exc}") from exc
        _require(staged_fd is not None, f"failed to bind staged {label} descriptor")
        audited = _audit_open_descriptor(staged_fd, planned, f"staged {label}")
        result_fd = staged_fd
        staged_fd = None
        committed = True
        return audited, result_fd
    finally:
        if staged_fd is not None:
            try:
                os.close(staged_fd)
            except OSError:
                pass
        if destination_installed and destination_signature is not None and not committed:
            try:
                current = os.stat(
                    destination.name, dir_fd=staging_directory_fd,
                    follow_symlinks=False)
            except FileNotFoundError:
                pass
            else:
                current_signature = (
                    current.st_dev, current.st_ino, current.st_size,
                    current.st_mtime_ns,
                )
                if (stat.S_ISREG(current.st_mode) and
                        current_signature == destination_signature):
                    os.unlink(destination.name, dir_fd=staging_directory_fd)
                    os.fsync(staging_directory_fd)
                else:
                    raise LaunchFailure(
                        f"refusing cleanup of changed staged {label} inode")


def _seal_staging_directory(staging_dir: Path, descriptor: int) -> dict[str, Any]:
    os.fchmod(descriptor, 0o555)
    os.fsync(descriptor)
    info = os.fstat(descriptor)
    _require(stat.S_IMODE(info.st_mode) == 0o555,
             "staging directory did not become read-only")
    return {
        "path": str(staging_dir),
        "device": info.st_dev,
        "inode": info.st_ino,
        "mode": "0555",
    }


def _cleanup_bound_staged_inputs(staging_fd: int,
                                 staging_record: Mapping[str, Any],
                                 staged_records: Sequence[Mapping[str, Any]]) -> None:
    """Remove only this launch's inode-bound direct children after a pre-spawn fail."""

    info = os.fstat(staging_fd)
    _require(stat.S_ISDIR(info.st_mode) and
             info.st_dev == staging_record.get("device") and
             info.st_ino == staging_record.get("inode"),
             "refusing cleanup through a changed staging directory descriptor")
    os.fchmod(staging_fd, 0o700)
    for record in staged_records:
        path = Path(str(record.get("path", "")))
        _require(path.name not in ("", ".", "..") and
                 str(path.parent) == str(staging_record.get("path")),
                 "refusing cleanup of a non-direct staged child")
        entry = os.stat(path.name, dir_fd=staging_fd, follow_symlinks=False)
        expected = tuple(record.get(name) for name in (
            "device", "inode", "size", "mtime_ns", "ctime_ns",
        ))
        actual = (entry.st_dev, entry.st_ino, entry.st_size,
                  entry.st_mtime_ns, entry.st_ctime_ns)
        _require(stat.S_ISREG(entry.st_mode) and actual == expected,
                 f"refusing cleanup of changed staged inode {path.name}")
        os.unlink(path.name, dir_fd=staging_fd)
    os.fsync(staging_fd)
    _require(not os.listdir(staging_fd),
             "staging directory is not empty after bounded pre-spawn cleanup")


@dataclass
class HolderInputs:
    pid: int
    start_time_ticks: int
    roles: dict[str, dict[str, Any]]
    start_time_probe: Callable[[int], int | None] = field(repr=False)
    closed: bool = False

    def audit_metadata(self) -> dict[str, Any]:
        _require(not self.closed, "input holder descriptors are already closed")
        current_start = self.start_time_probe(self.pid)
        _require(current_start == self.start_time_ticks,
                 "input holder process identity changed")
        roles: dict[str, dict[str, Any]] = {}
        for role, record in self.roles.items():
            fd = record["fd"]
            try:
                flags = fcntl.fcntl(fd, fcntl.F_GETFL)
                _require(flags & os.O_ACCMODE == os.O_RDONLY,
                         f"holder {role} descriptor is not read-only")
                info = os.fstat(fd)
            except (OSError, ValueError, RuntimeError) as exc:
                if isinstance(exc, LaunchFailure):
                    raise
                raise LaunchFailure(
                    f"cannot audit holder {role} descriptor: {exc}") from exc
            _require(_file_signature(record) == (
                info.st_dev, info.st_ino, info.st_size,
                info.st_mtime_ns, info.st_ctime_ns,
            ), f"holder {role} descriptor identity changed")
            roles[role] = {
                "fd": fd, "proc_path": record["proc_path"], "role": role,
                "access_mode": "read_only",
                "device": info.st_dev, "inode": info.st_ino,
                "size": info.st_size, "mtime_ns": info.st_mtime_ns,
                "ctime_ns": info.st_ctime_ns,
            }
        return {
            "kind": INPUT_TRANSPORT_KIND,
            "holder_pid": self.pid,
            "holder_start_time_ticks": self.start_time_ticks,
            "roles": roles,
        }

    def audit(self) -> dict[str, Any]:
        metadata = self.audit_metadata()
        audited_roles: dict[str, dict[str, Any]] = {}
        for role, record in self.roles.items():
            fd = record["fd"]
            audited = _audit_open_descriptor(fd, record, f"holder {role}")
            _require(all(audited[name] == record[name] for name in (
                "device", "inode", "size", "mtime_ns", "ctime_ns", "sha256",
            )), f"holder {role} descriptor identity changed")
            audited_roles[role] = {
                **audited,
                **metadata["roles"][role],
            }
        return {
            **{name: metadata[name] for name in (
                "kind", "holder_pid", "holder_start_time_ticks")},
            "roles": audited_roles,
        }

    def close(self) -> None:
        if self.closed:
            return
        for record in self.roles.values():
            try:
                os.close(record["fd"])
            except OSError:
                pass
        self.closed = True

    def prove_closed(self) -> dict[str, Any]:
        _require(self.closed, "holder descriptors must be closed before closure proof")
        states: dict[str, dict[str, Any]] = {}
        for role, record in self.roles.items():
            try:
                os.fstat(record["fd"])
            except OSError:
                gone = True
            else:
                gone = False
            _require(gone, f"holder {role} descriptor remained open")
            states[role] = {"fd": record["fd"], "closed": True}
        return {"all_holder_fds_closed": True, "roles": states}

    def __del__(self) -> None:  # pragma: no cover - defensive leak prevention
        self.close()


@dataclass
class HolderDirectories:
    """Stable state/evidence directory identities held by fixed procfs FDs."""

    pid: int
    start_time_ticks: int
    roles: dict[str, dict[str, Any]]
    start_time_probe: Callable[[int], int | None] = field(repr=False)
    closed: bool = False

    def audit(self, *, require_paths: bool = True) -> dict[str, Any]:
        _require(not self.closed, "directory holder descriptors are already closed")
        _require(self.start_time_probe(self.pid) == self.start_time_ticks,
                 "directory holder process identity changed")
        audited: dict[str, dict[str, Any]] = {}
        for role, record in self.roles.items():
            try:
                info = os.fstat(record["fd"])
            except OSError as exc:
                raise LaunchFailure(
                    f"cannot audit holder {role} directory descriptor: {exc}") from exc
            expected = (
                record["device"], record["inode"], record["owner_uid"],
                int(record["mode"], 8),
            )
            _require(stat.S_ISDIR(info.st_mode) and
                     _directory_signature(info) == expected,
                     f"holder {role} directory identity changed")
            if require_paths:
                try:
                    path_info = Path(record["path"]).lstat()
                except OSError as exc:
                    raise LaunchFailure(
                        f"planned {role} directory path disappeared: {exc}") from exc
                _require(_directory_signature(path_info) == expected,
                         f"planned {role} directory path was replaced")
            audited[role] = {
                **record, "access_mode": "read_only_directory_descriptor",
            }
        return {
            "kind": DIRECTORY_TRANSPORT_KIND,
            "holder_pid": self.pid,
            "holder_start_time_ticks": self.start_time_ticks,
            "roles": audited,
        }

    def close(self) -> None:
        if self.closed:
            return
        for record in self.roles.values():
            try:
                os.close(record["fd"])
            except OSError:
                pass
        self.closed = True

    def prove_closed(self) -> dict[str, Any]:
        _require(self.closed, "directory holder must be closed before proof")
        states: dict[str, dict[str, Any]] = {}
        for role, record in self.roles.items():
            try:
                os.fstat(record["fd"])
            except OSError:
                gone = True
            else:
                gone = False
            _require(gone, f"holder {role} directory descriptor remained open")
            states[role] = {"fd": record["fd"], "closed": True}
        return {"all_directory_fds_closed": True, "roles": states}

    def __del__(self) -> None:  # pragma: no cover - defensive leak prevention
        self.close()


@dataclass
class HolderExecutables:
    """Read-only executable descriptors immune to pathname replacement."""

    pid: int
    start_time_ticks: int
    roles: dict[str, dict[str, Any]]
    start_time_probe: Callable[[int], int | None] = field(repr=False)
    closed: bool = False

    def audit(self, *, hash_content: bool = True) -> dict[str, Any]:
        _require(not self.closed, "executable holder descriptors are already closed")
        _require(self.start_time_probe(self.pid) == self.start_time_ticks,
                 "executable holder process identity changed")
        audited: dict[str, dict[str, Any]] = {}
        for role, record in self.roles.items():
            if hash_content:
                binding = _audit_open_descriptor(
                    record["fd"], record, f"holder executable {role}")
            else:
                info = os.fstat(record["fd"])
                _require(_file_signature(record) == (
                    info.st_dev, info.st_ino, info.st_size,
                    info.st_mtime_ns, info.st_ctime_ns,
                ), f"holder executable {role} identity changed")
                binding = {name: record[name] for name in (
                    "path", "device", "inode", "size", "mtime_ns", "ctime_ns",
                    "sha256", "closure_check",
                )}
            audited[role] = {
                **binding, "role": role, "fd": record["fd"],
                "proc_path": record["proc_path"], "access_mode": "read_only",
            }
        return {
            "kind": EXECUTABLE_TRANSPORT_KIND,
            "holder_pid": self.pid,
            "holder_start_time_ticks": self.start_time_ticks,
            "roles": audited,
        }

    def close(self) -> None:
        if self.closed:
            return
        for record in self.roles.values():
            try:
                os.close(record["fd"])
            except OSError:
                pass
        self.closed = True

    def __del__(self) -> None:  # pragma: no cover
        self.close()


def _fd_is_open(descriptor: int) -> bool:
    try:
        os.fstat(descriptor)
    except OSError:
        return False
    return True


def _install_holder_inputs(source_fd: int, trajectory_fd: int,
                           source_record: Mapping[str, Any],
                           trajectory_record: Mapping[str, Any],
                           inspector: "Inspector") -> HolderInputs:
    for descriptor in (SOURCE_RESTART_FD, TRAJECTORY_FD):
        _require(not _fd_is_open(descriptor),
                 f"refusing to replace existing process fd {descriptor}")
    installed: list[int] = []
    try:
        os.dup2(source_fd, SOURCE_RESTART_FD, inheritable=False)
        installed.append(SOURCE_RESTART_FD)
        os.dup2(trajectory_fd, TRAJECTORY_FD, inheritable=False)
        installed.append(TRAJECTORY_FD)
        holder_pid = os.getpid()
        holder_start = inspector.start_time_ticks(holder_pid)
        _require(isinstance(holder_start, int),
                 "cannot bind launcher holder process starttime")
        roles: dict[str, dict[str, Any]] = {}
        for role, descriptor, planned in (
            ("source_restart", SOURCE_RESTART_FD, source_record),
            ("trajectory", TRAJECTORY_FD, trajectory_record),
        ):
            info = os.fstat(descriptor)
            _require(_file_signature(planned) == (
                info.st_dev, info.st_ino, info.st_size,
                info.st_mtime_ns, info.st_ctime_ns,
            ), f"holder {role} descriptor does not bind the staged inode")
            roles[role] = {
                **dict(planned),
                "fd": descriptor,
                "proc_path": f"/proc/{holder_pid}/fd/{descriptor}",
                "role": role,
            }
        holder = HolderInputs(
            holder_pid, holder_start, roles, inspector.start_time_ticks)
        holder.audit_metadata()
        return holder
    except BaseException:
        for descriptor in installed:
            try:
                os.close(descriptor)
            except OSError:
                pass
        raise


def _install_holder_directories(
        state_fd: int, evidence_fd: int,
        state_record: Mapping[str, Any], evidence_record: Mapping[str, Any],
        inspector: "Inspector") -> HolderDirectories:
    for descriptor in (STATE_DIRECTORY_FD, EVIDENCE_DIRECTORY_FD):
        _require(not _fd_is_open(descriptor),
                 f"refusing to replace existing process fd {descriptor}")
    installed: list[int] = []
    try:
        os.dup2(state_fd, STATE_DIRECTORY_FD, inheritable=False)
        installed.append(STATE_DIRECTORY_FD)
        os.dup2(evidence_fd, EVIDENCE_DIRECTORY_FD, inheritable=False)
        installed.append(EVIDENCE_DIRECTORY_FD)
        holder_pid = os.getpid()
        holder_start = inspector.start_time_ticks(holder_pid)
        _require(isinstance(holder_start, int),
                 "cannot bind directory holder process starttime")
        roles: dict[str, dict[str, Any]] = {}
        for role, descriptor, planned in (
            ("state_dir", STATE_DIRECTORY_FD, state_record),
            ("evidence_dir", EVIDENCE_DIRECTORY_FD, evidence_record),
        ):
            info = os.fstat(descriptor)
            expected = (
                planned["device"], planned["inode"], planned["owner_uid"],
                int(planned["mode"], 8),
            )
            _require(_directory_signature(info) == expected,
                     f"holder {role} descriptor does not bind planned directory")
            roles[role] = {
                **dict(planned), "fd": descriptor, "role": role,
                "proc_path": f"/proc/{holder_pid}/fd/{descriptor}",
            }
        holder = HolderDirectories(
            holder_pid, holder_start, roles, inspector.start_time_ticks)
        holder.audit()
        return holder
    except BaseException:
        for descriptor in installed:
            try:
                os.close(descriptor)
            except OSError:
                pass
        raise


def _open_executable_descriptor(record: Mapping[str, Any], label: str) -> int:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(record["path"], flags)
    except OSError as exc:
        raise LaunchFailure(f"cannot bind {label} executable descriptor: {exc}") from exc
    try:
        _audit_open_descriptor(descriptor, record, label)
    except BaseException:
        os.close(descriptor)
        raise
    return descriptor


def _install_holder_executables(
        launcher_fd: int, binary_fd: int,
        launcher_record: Mapping[str, Any], binary_record: Mapping[str, Any],
        inspector: "Inspector") -> HolderExecutables:
    for descriptor in (MPI_LAUNCHER_FD, BINARY_EXECUTABLE_FD):
        _require(not _fd_is_open(descriptor),
                 f"refusing to replace existing process fd {descriptor}")
    installed: list[int] = []
    try:
        os.dup2(launcher_fd, MPI_LAUNCHER_FD, inheritable=False)
        installed.append(MPI_LAUNCHER_FD)
        os.dup2(binary_fd, BINARY_EXECUTABLE_FD, inheritable=False)
        installed.append(BINARY_EXECUTABLE_FD)
        holder_pid = os.getpid()
        holder_start = inspector.start_time_ticks(holder_pid)
        _require(isinstance(holder_start, int),
                 "cannot bind executable holder process starttime")
        roles: dict[str, dict[str, Any]] = {}
        for role, descriptor, planned in (
            ("launcher", MPI_LAUNCHER_FD, launcher_record),
            ("binary", BINARY_EXECUTABLE_FD, binary_record),
        ):
            info = os.fstat(descriptor)
            _require(_file_signature(planned) == (
                info.st_dev, info.st_ino, info.st_size,
                info.st_mtime_ns, info.st_ctime_ns,
            ), f"holder executable {role} does not bind planned inode")
            roles[role] = {
                **dict(planned), "role": role, "fd": descriptor,
                "proc_path": f"/proc/{holder_pid}/fd/{descriptor}",
            }
        holder = HolderExecutables(
            holder_pid, holder_start, roles, inspector.start_time_ticks)
        holder.audit(hash_content=False)
        return holder
    except BaseException:
        for descriptor in installed:
            try:
                os.close(descriptor)
            except OSError:
                pass
        raise


def _probe_proc_holder_access(
        input_holder: HolderInputs, directory_holder: HolderDirectories,
        executable_holder: HolderExecutables) -> dict[str, Any]:
    """Fork a real peer which reopens every holder FD through parent procfs."""

    _require(sys.platform == "linux" and hasattr(os, "fork"),
             "proc-holder access probe requires Linux fork/procfs")
    role_sets = (
        ("input", input_holder.roles),
        ("directory", directory_holder.roles),
        ("executable", executable_holder.roles),
    )
    read_fd, write_fd = os.pipe2(getattr(os, "O_CLOEXEC", 0))
    child = os.fork()
    if child == 0:  # pragma: no cover - assertions are evaluated in parent
        try:
            os.close(read_fd)
            samples: dict[tuple[str, str], tuple[bytes, bytes]] = {}
            descriptors: list[int] = []
            for family, roles in role_sets:
                for role, record in roles.items():
                    descriptor = int(record["fd"])
                    descriptors.append(descriptor)
                    if family != "directory":
                        info = os.fstat(descriptor)
                        first = os.pread(descriptor, min(64, info.st_size), 0)
                        tail_offset = max(0, info.st_size - 64)
                        last = os.pread(descriptor, min(64, info.st_size), tail_offset)
                        samples[(family, role)] = (first, last)
            for descriptor in descriptors:
                os.close(descriptor)
            for family, roles in role_sets:
                for role, record in roles.items():
                    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
                    if family == "directory":
                        flags |= os.O_DIRECTORY
                    if family == "executable":
                        try:
                            executable_access = os.access(
                                record["proc_path"], os.X_OK, effective_ids=True)
                        except TypeError:  # pragma: no cover - older Python fallback
                            executable_access = os.access(record["proc_path"], os.X_OK)
                        if not executable_access:
                            raise RuntimeError(
                                f"{family}:{role} proc path is not executable")
                    reopened = os.open(record["proc_path"], flags)
                    try:
                        info = os.fstat(reopened)
                        if family == "directory":
                            if not stat.S_ISDIR(info.st_mode):
                                raise RuntimeError(f"{family}:{role} is not a directory")
                        else:
                            expected = _file_signature(record)
                            actual = (
                                info.st_dev, info.st_ino, info.st_size,
                                info.st_mtime_ns, info.st_ctime_ns,
                            )
                            if actual != expected:
                                raise RuntimeError(
                                    f"{family}:{role} descriptor identity changed")
                            first, last = samples[(family, role)]
                            if (os.pread(reopened, len(first), 0) != first or
                                    os.pread(reopened, len(last),
                                             max(0, info.st_size - 64)) != last):
                                raise RuntimeError(
                                    f"{family}:{role} sample bytes changed")
                    finally:
                        os.close(reopened)
            os.write(write_fd, b"ok")
            os._exit(0)
        except BaseException as exc:
            try:
                os.write(write_fd, (f"error:{type(exc).__name__}:{exc}").encode(
                    "utf-8", errors="replace")[:4096])
            finally:
                os._exit(1)
    os.close(write_fd)
    try:
        chunks: list[bytes] = []
        while True:
            chunk = os.read(read_fd, 4096)
            if not chunk:
                break
            chunks.append(chunk)
    finally:
        os.close(read_fd)
    _, status = os.waitpid(child, 0)
    detail = b"".join(chunks).decode("utf-8", errors="replace")
    _require(os.WIFEXITED(status) and os.WEXITSTATUS(status) == 0 and detail == "ok",
             f"peer proc-holder access probe failed: {detail or status}")
    return {
        "kind": "fork_peer_procfs_reopen_v1",
        "peer_pid": child,
        "holder_pid": input_holder.pid,
        "families": {
            family: sorted(roles) for family, roles in role_sets
        },
        "all_reopened_and_sampled": True,
        "executable_access": {
            role: "effective_uid_x_ok"
            for role in sorted(executable_holder.roles)
        },
    }


def _file_signature(record: Mapping[str, Any]) -> tuple[int, int, int, int, int]:
    return tuple(int(record[name]) for name in
                 ("device", "inode", "size", "mtime_ns", "ctime_ns"))


def _assert_path_binding(record: Mapping[str, Any], label: str) -> None:
    path = Path(str(record["path"]))
    try:
        current = path.stat()
    except OSError as exc:
        raise LaunchFailure(f"cannot restat {label}: {exc}") from exc
    actual = (current.st_dev, current.st_ino, current.st_size,
              current.st_mtime_ns, current.st_ctime_ns)
    _require(actual == _file_signature(record), f"{label} changed after preflight")


@dataclass(frozen=True)
class GPURecord:
    index: int
    uuid: str
    pci_bus_id: str
    cuda_ordinal: int
    uncorrected_ecc: int
    corrected_ecc: int
    memory_total_mib: int
    memory_used_mib: int

    def csv_line(self) -> str:
        return (f"{self.index},{self.uuid},{self.pci_bus_id},{self.cuda_ordinal},"
                f"{self.uncorrected_ecc},"
                f"{self.corrected_ecc},{self.memory_total_mib},"
                f"{self.memory_used_mib}\n")

    def as_dict(self) -> dict[str, Any]:
        return {
            "index": self.index,
            "uuid": self.uuid,
            "pci_bus_id": self.pci_bus_id,
            "cuda_ordinal": self.cuda_ordinal,
            "uncorrected_ecc": self.uncorrected_ecc,
            "corrected_ecc": self.corrected_ecc,
            "memory_total_mib": self.memory_total_mib,
            "memory_used_mib": self.memory_used_mib,
        }


@dataclass(frozen=True)
class GPUApplication:
    pid: int
    uuid: str


def _nvidia_query(query: str, runtime: "Runtime") -> list[list[str]]:
    _require(isinstance(runtime.nvidia_smi, str) and
             Path(runtime.nvidia_smi).is_absolute(),
             "nvidia-smi must be an absolute plan-bound executable")
    output = _run_text([
        runtime.nvidia_smi, f"--query-gpu={query}",
        "--format=csv,noheader,nounits",
    ], runtime, environment=runtime.nvidia_environment)
    if not output:
        return []
    return [[field.strip() for field in line.split(",")]
            for line in output.splitlines() if line.strip()]


def _pci_bus_key(value: str) -> tuple[int, int, int, int]:
    match = re.fullmatch(
        r"(?:0[xX])?([0-9A-Fa-f]{4,8}):([0-9A-Fa-f]{2}):"
        r"([0-9A-Fa-f]{2})\.([0-7])", value,
    )
    _require(match is not None, f"invalid NVIDIA PCI bus id: {value!r}")
    assert match is not None
    return tuple(int(component, 16) for component in match.groups())


def query_gpu_inventory(runtime: "Runtime") -> list[GPURecord]:
    rows = _nvidia_query(GPU_QUERY, runtime)
    parsed: list[tuple[int, str, str, int, int, int, int]] = []
    for row in rows:
        _require(len(row) == 7, f"unexpected nvidia-smi GPU row: {row!r}")
        try:
            values = (
                int(row[0]), row[1], row[2], int(row[3]), int(row[4]),
                int(row[5]), int(row[6])
            )
        except ValueError as exc:
            raise LaunchFailure(f"non-integer nvidia-smi GPU row: {row!r}") from exc
        (index, uuid, pci_bus_id, uncorrected, corrected,
         memory_total, memory) = values
        _pci_bus_key(pci_bus_id)
        _require(index >= 0 and uuid.startswith("GPU-") and
                 uncorrected >= 0 and corrected >= 0 and
                 memory_total > 0 and 0 <= memory <= memory_total,
                 f"invalid nvidia-smi GPU row: {row!r}")
        parsed.append(values)
    parsed.sort(key=lambda item: _pci_bus_key(item[2]))
    result = [
        GPURecord(index, uuid, pci_bus_id, ordinal, uncorrected, corrected,
                  memory_total, memory)
        for ordinal, (index, uuid, pci_bus_id, uncorrected, corrected,
                      memory_total, memory)
        in enumerate(parsed)
    ]
    _require(len({item.index for item in result}) == len(result),
             "nvidia-smi returned duplicate GPU indices")
    _require(len({item.uuid for item in result}) == len(result),
             "nvidia-smi returned duplicate GPU UUIDs")
    _require(len({item.pci_bus_id.lower() for item in result}) == len(result),
             "nvidia-smi returned duplicate GPU PCI bus ids")
    return result


def _gpu_capacity_preflight(
        gpus: Sequence[GPURecord],
        capacity_transition: Mapping[str, Any]) -> dict[str, Any]:
    model = capacity_transition["gpu_memory_model"]
    minimum_total = int(model["minimum_gpu_memory_total_mib"])
    required = int(model["required_per_rank_memory_mib_ceiling"])
    _require(bool(gpus), "GPU capacity preflight requires at least one device")
    _require(all(gpu.memory_total_mib >= minimum_total for gpu in gpus),
             "preflight GPU total memory is below the plan-bound 80% capacity "
             "model")
    _require(all(gpu.memory_used_mib + required <= gpu.memory_total_mib
                 for gpu in gpus),
             "preflight GPU available memory is below the plan-bound capacity "
             "requirement")
    return {
        "kind": "per_rank_memory_total_and_available_gate_v1",
        "required_per_rank_memory_mib_ceiling": required,
        "minimum_gpu_memory_total_mib": minimum_total,
        "minimum_observed_gpu_memory_total_mib":
            min(gpu.memory_total_mib for gpu in gpus),
        "maximum_observed_gpu_memory_used_mib":
            max(gpu.memory_used_mib for gpu in gpus),
        "minimum_observed_gpu_memory_available_mib":
            min(gpu.memory_total_mib - gpu.memory_used_mib for gpu in gpus),
        "all_ranks_pass": True,
    }


def query_gpu_applications(runtime: "Runtime") -> list[GPUApplication]:
    _require(isinstance(runtime.nvidia_smi, str) and
             Path(runtime.nvidia_smi).is_absolute(),
             "nvidia-smi must be an absolute plan-bound executable")
    try:
        result = runtime.run([
            runtime.nvidia_smi, f"--query-compute-apps={GPU_APP_QUERY}",
            "--format=csv,noheader,nounits",
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
           **({"env": dict(runtime.nvidia_environment)}
              if runtime.nvidia_environment is not None else {}))
    except (OSError, subprocess.CalledProcessError) as exc:
        detail = getattr(exc, "stderr", "") or str(exc)
        raise LaunchFailure(f"cannot query GPU applications: {detail.strip()}") from exc
    applications: list[GPUApplication] = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        row = [field.strip() for field in line.split(",")]
        _require(len(row) == 2, f"unexpected GPU application row: {row!r}")
        try:
            pid = int(row[0])
        except ValueError as exc:
            raise LaunchFailure(f"invalid GPU application PID: {row!r}") from exc
        _require(pid > 0 and row[1].startswith("GPU-"),
                 f"invalid GPU application row: {row!r}")
        applications.append(GPUApplication(pid, row[1]))
    return applications


class Inspector(Protocol):
    def descendants(self, root_pid: int) -> set[int]: ...
    def cmdline(self, pid: int) -> list[str]: ...
    def environment(self, pid: int) -> dict[str, str]: ...
    def executable(self, pid: int) -> dict[str, Any]: ...
    def start_time_ticks(self, pid: int) -> int | None: ...


class LinuxProcInspector:
    """Read process identity directly from Linux procfs."""

    def __init__(self, proc_root: Path = Path("/proc")) -> None:
        self.proc_root = proc_root

    def _read_bytes(self, pid: int, name: str) -> bytes:
        try:
            return (self.proc_root / str(pid) / name).read_bytes()
        except (OSError, ProcessLookupError) as exc:
            raise LaunchFailure(f"cannot read /proc/{pid}/{name}: {exc}") from exc

    def cmdline(self, pid: int) -> list[str]:
        raw = self._read_bytes(pid, "cmdline")
        try:
            return [token.decode("utf-8") for token in raw.rstrip(b"\0").split(b"\0")]
        except UnicodeDecodeError as exc:
            raise LaunchFailure(f"/proc/{pid}/cmdline is not UTF-8") from exc

    def environment(self, pid: int) -> dict[str, str]:
        raw = self._read_bytes(pid, "environ")
        result: dict[str, str] = {}
        try:
            tokens = raw.rstrip(b"\0").split(b"\0")
            for token in tokens:
                name, separator, value = token.partition(b"=")
                if separator:
                    result[name.decode("utf-8")] = value.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise LaunchFailure(f"/proc/{pid}/environ is not UTF-8") from exc
        return result

    def descendants(self, root_pid: int) -> set[int]:
        parents: dict[int, int] = {}
        try:
            processes = list(self.proc_root.iterdir())
        except OSError as exc:
            raise LaunchFailure(f"cannot enumerate procfs: {exc}") from exc
        for process in processes:
            if not process.name.isdigit():
                continue
            try:
                # Field 2 is parenthesized and may contain spaces or ')'.
                text = (process / "stat").read_text(encoding="ascii")
                suffix = text[text.rfind(")") + 2:].split()
                parents[int(process.name)] = int(suffix[1])
            except (OSError, ValueError, IndexError, UnicodeDecodeError):
                continue
        descendants: set[int] = set()
        changed = True
        while changed:
            changed = False
            for pid, parent in parents.items():
                if pid not in descendants and (parent == root_pid or parent in descendants):
                    descendants.add(pid)
                    changed = True
        return descendants

    def start_time_ticks(self, pid: int) -> int | None:
        try:
            text = (self.proc_root / str(pid) / "stat").read_text(encoding="ascii")
        except (FileNotFoundError, ProcessLookupError):
            return None
        except OSError as exc:
            raise LaunchFailure(f"cannot read /proc/{pid}/stat: {exc}") from exc
        try:
            suffix = text[text.rfind(")") + 2:].split()
            # proc(5) field 22; suffix starts at field 3 (state).
            start_time = int(suffix[19])
        except (ValueError, IndexError) as exc:
            raise LaunchFailure(f"cannot parse /proc/{pid}/stat starttime") from exc
        _require(start_time >= 0, f"/proc/{pid}/stat has negative starttime")
        return start_time

    def executable(self, pid: int) -> dict[str, Any]:
        path = self.proc_root / str(pid) / "exe"
        try:
            descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_CLOEXEC", 0))
        except OSError as exc:
            raise LaunchFailure(f"cannot open /proc/{pid}/exe: {exc}") from exc
        digest = hashlib.sha256()
        try:
            before = os.fstat(descriptor)
            with os.fdopen(descriptor, "rb", closefd=True) as stream:
                for chunk in iter(lambda: stream.read(HASH_CHUNK_BYTES), b""):
                    digest.update(chunk)
                after = os.fstat(stream.fileno())
        except BaseException:
            try:
                os.close(descriptor)
            except OSError:
                pass
            raise
        signature = (before.st_dev, before.st_ino, before.st_size,
                     before.st_mtime_ns, before.st_ctime_ns)
        _require(signature == (
            after.st_dev, after.st_ino, after.st_size,
            after.st_mtime_ns, after.st_ctime_ns,
        ), f"/proc/{pid}/exe changed while hashing")
        return {
            "device": before.st_dev,
            "inode": before.st_ino,
            "size": before.st_size,
            "mtime_ns": before.st_mtime_ns,
            "ctime_ns": before.st_ctime_ns,
            "sha256": digest.hexdigest(),
        }


@dataclass
class Runtime:
    run: Callable[..., Any] = subprocess.run
    popen: Callable[..., Any] = subprocess.Popen
    sleep: Callable[[float], None] = time.sleep
    monotonic: Callable[[], float] = time.monotonic
    hostname: Callable[[], str] = socket.gethostname
    environment: Callable[[], Mapping[str, str]] = lambda: os.environ
    inspector: Inspector = field(default_factory=LinuxProcInspector)
    nvidia_smi: str | None = None
    nvidia_environment: Mapping[str, str] | None = None
    getpgid: Callable[[int], int] = os.getpgid
    killpg: Callable[[int, int], None] = os.killpg
    get_signal_handler: Callable[[int], Any] = signal.getsignal
    set_signal_handler: Callable[[int, Any], Any] = signal.signal
    statvfs: Callable[[Any], Any] = os.statvfs
    fstatvfs: Callable[[int], Any] = os.fstatvfs
    plan_validator: Callable[[dict[str, Any]], dict[str, Any]] = \
        _production_validate_plan
    source_validator: Callable[[Path, Mapping[str, Any]], dict[str, Any]] = \
        _production_validate_source


@dataclass(frozen=True)
class PreparedLaunch:
    plan: dict[str, Any]
    plan_binding: dict[str, Any]
    state_dir: Path
    world_size: int
    wall_time_seconds: int
    plan_validation: dict[str, Any]
    source_validation: dict[str, Any]
    proc_access_probe: dict[str, Any]
    launch_environment: dict[str, str]
    launch_environment_sha256: str
    mca_configuration_preflight: dict[str, Any]
    execution_tools: dict[str, dict[str, Any]]
    repository: dict[str, Any]
    binary: dict[str, Any]
    source_restart: dict[str, Any]
    trajectory: dict[str, Any]
    staging_directory: dict[str, Any]
    staged_source_restart: dict[str, Any]
    staged_trajectory: dict[str, Any]
    staging_directory_fd: int = field(repr=False)
    disk_preflight_before_staging: dict[str, Any]
    input_holder: HolderInputs
    directory_holder: HolderDirectories
    executable_holder: HolderExecutables
    launcher: dict[str, Any]
    mpi_argv: tuple[str, ...]
    athena_argv: tuple[str, ...]
    launch_argv: tuple[str, ...]
    gpus: tuple[GPURecord, ...]
    gpu_capacity_preflight: dict[str, Any]
    gpu_visibility_environment: dict[str, str | None]

    def close(self) -> None:
        self.input_holder.close()
        self.directory_holder.close()
        self.executable_holder.close()
        descriptor = self.staging_directory_fd
        if descriptor >= 0:
            object.__setattr__(self, "staging_directory_fd", -1)
            if _fd_is_open(descriptor):
                os.close(descriptor)

    def __del__(self) -> None:  # pragma: no cover - defensive library-use fallback
        try:
            self.close()
        except (AttributeError, OSError):
            pass


def prepare_launch(plan_path: Path, state_dir: Path,
                   runtime: Runtime | None = None) -> PreparedLaunch:
    runtime = runtime or Runtime()
    plan, plan_binding = _read_immutable_plan(plan_path)
    _require(plan.get("schema") == SCHEMA, f"plan schema must equal {SCHEMA}")
    inputs = plan.get("inputs")
    expected = plan.get("expected")
    policy = plan.get("policy")
    tools = plan.get("tools")
    contract = plan.get("launch_contract")
    _require(isinstance(inputs, dict), "plan inputs must be an object")
    _require(isinstance(expected, dict), "plan expected must be an object")
    _require(isinstance(policy, dict), "plan policy must be an object")
    _require(isinstance(tools, dict) and tools.get("hash_algorithm") == "sha256",
             "plan tools must bind SHA-256 tool identities")
    _require(isinstance(contract, dict), "plan launch_contract must be an object")

    execution_tool_paths = {
        "segment_launcher": SCRIPT_PATH,
        "segment_checker": SEGMENT_CHECKER_PATH,
        "output_integrity": OUTPUT_INTEGRITY_PATH,
        "restart_auditor": RESTART_AUDITOR_PATH,
        "restart_metadata_reader": RESTART_READER_PATH,
        "nvidia_smi": Path(str(tools.get("nvidia_smi", {}).get("path", ""))),
    }
    execution_tools: dict[str, dict[str, Any]] = {}
    for name, actual_path in execution_tool_paths.items():
        binding = _audit_planned_file(
            tools.get(name), f"tools.{name}", executable=True)
        _require(binding["path"] == str(actual_path),
                 f"current {name} differs from plan-bound path")
        execution_tools[name] = binding
    planned_nvidia_smi = execution_tools["nvidia_smi"]["path"]
    _require(runtime.nvidia_smi in (None, planned_nvidia_smi),
             "runtime nvidia-smi differs from immutable plan")
    runtime.nvidia_smi = planned_nvidia_smi

    plan_validation = runtime.plan_validator(plan)

    launch_environment = _validate_launch_environment(contract.get("environment"))
    launch_environment_sha256 = _canonical_sha256(launch_environment)
    nvidia_environment = {
        key: launch_environment[key] for key in ("HOME", "LANG", "LC_ALL")
    }
    _require(runtime.nvidia_environment in (None, nvidia_environment),
             "runtime nvidia-smi environment differs from immutable plan")
    runtime.nvidia_environment = nvidia_environment

    resolved_state = _safe_directory(state_dir, "state directory")
    try:
        state_entries = list(resolved_state.iterdir())
    except OSError as exc:
        raise LaunchFailure(f"cannot inventory state directory: {exc}") from exc
    _require(not state_entries, "state directory must be empty before segment launch")
    world_size = _integer(policy.get("ranks"), "policy.ranks")
    _require(world_size > 0, "policy.ranks must be positive")
    wall_time = _integer(contract.get("wall_time_seconds"),
                         "launch_contract.wall_time_seconds")
    _require(wall_time > 0, "launch wall time must be positive")
    _require(contract.get("state_dir") == str(resolved_state),
             "state directory differs from immutable launch contract")
    _require(contract.get("world_size") == world_size and
             contract.get("gpu_count") == world_size,
             "launch contract rank/GPU counts differ from policy.ranks")

    evidence_value = contract.get("evidence")
    evidence_dir_value = contract.get("evidence_dir")
    _require(isinstance(evidence_value, dict) and
             set(evidence_value) == set(EVIDENCE_NAMES),
             "launch contract evidence set is not exact")
    _require(isinstance(evidence_dir_value, str) and
             Path(evidence_dir_value).is_absolute(),
             "launch contract evidence_dir must be absolute")
    resolved_evidence = _safe_directory(
        Path(evidence_dir_value), "evidence directory")
    _require(str(resolved_evidence) == evidence_dir_value,
             "launch contract evidence_dir is not canonical")
    for name in EVIDENCE_NAMES:
        value = evidence_value[name]
        _require(isinstance(value, str) and Path(value).is_absolute() and
                 Path(value).parent == resolved_evidence and
                 Path(value).name not in ("", ".", ".."),
                 f"launch contract evidence.{name} is not a direct canonical child")

    directory_transport = contract.get("directory_transport")
    _require(isinstance(directory_transport, dict) and
             directory_transport.get("kind") == DIRECTORY_TRANSPORT_KIND and
             directory_transport.get("holder_pid_token") == HOLDER_PID_TOKEN,
             f"directory transport must be {DIRECTORY_TRANSPORT_KIND}")
    directory_roles = directory_transport.get("roles")
    _require(isinstance(directory_roles, dict) and
             set(directory_roles) == {"state_dir", "evidence_dir"},
             "directory transport roles must be exactly state_dir/evidence_dir")
    for role, descriptor, planned_path in (
        ("state_dir", STATE_DIRECTORY_FD, resolved_state),
        ("evidence_dir", EVIDENCE_DIRECTORY_FD, resolved_evidence),
    ):
        role_value = directory_roles.get(role)
        _require(isinstance(role_value, dict) and
                 role_value == {
                     "role": role, "fd": descriptor,
                     "planned_path": str(planned_path),
                     "proc_path_template":
                         f"/proc/{HOLDER_PID_TOKEN}/fd/{descriptor}",
                 }, f"directory transport role {role} is not canonical")

    executable_transport = contract.get("executable_transport")
    _require(isinstance(executable_transport, dict) and
             executable_transport.get("kind") == EXECUTABLE_TRANSPORT_KIND and
             executable_transport.get("holder_pid_token") == HOLDER_PID_TOKEN,
             f"executable transport must be {EXECUTABLE_TRANSPORT_KIND}")
    executable_roles = executable_transport.get("roles")
    _require(isinstance(executable_roles, dict) and
             set(executable_roles) == {"launcher", "binary"},
             "executable transport roles must be exactly launcher/binary")
    for role, descriptor, parent in (
        ("launcher", MPI_LAUNCHER_FD, contract.get("launcher")),
        ("binary", BINARY_EXECUTABLE_FD, inputs.get("binary")),
    ):
        role_value = executable_roles.get(role)
        _require(isinstance(parent, dict) and isinstance(role_value, dict) and
                 role_value == {
                     "role": role, "fd": descriptor,
                     "parent_path": parent.get("path"),
                     "proc_path_template":
                         f"/proc/{HOLDER_PID_TOKEN}/fd/{descriptor}",
                 }, f"executable transport role {role} is not canonical")
    repository = _audit_repository(inputs.get("repo"), tools.get("git"), runtime)
    binary = _audit_planned_file(inputs.get("binary"), "inputs.binary",
                                 executable=True)
    source = _audit_planned_file(inputs.get("source_restart"),
                                 "inputs.source_restart")
    source_validation = runtime.source_validator(
        Path(source["path"]), plan)
    trajectory = _audit_planned_file(inputs.get("trajectory"), "inputs.trajectory")
    launcher_plan = contract.get("launcher")
    launcher = _audit_planned_file(launcher_plan, "launch_contract.launcher",
                                   executable=True)
    mca_configuration_preflight = _audit_mca_configuration(
        contract.get("mca_configuration"), launch_environment["HOME"])
    final_cycle = _integer(expected.get("final_cycle"), "expected.final_cycle")
    tlim = _finite(expected.get("tlim"), "expected.tlim")
    root_dt = _finite(expected.get("root_dt"), "expected.root_dt")
    _require(final_cycle >= 0 and tlim > 0.0,
             "planned cycle/time limits must be nonnegative/positive")
    _require(root_dt > 0.0, "planned root_dt must be positive")
    capacity_transition = _capacity_transition_contract(plan)
    restart_cadence_transition = _restart_cadence_transition_contract(
        plan, root_dt)
    expected_overrides = {
        "time/nlim": final_cycle,
        "time/tlim": tlim,
        "output3/dt": root_dt,
    }
    if restart_cadence_transition["kind"] == "tighten_v1":
        expected_overrides[restart_cadence_transition["parameter"]] = \
            restart_cadence_transition["target_dt"]
    if capacity_transition["kind"] == "increase_v1":
        expected_overrides["mesh_refinement/max_nmb_per_rank"] = \
            capacity_transition["target_max_nmb_per_rank"]
    overrides = plan.get("command_overrides")
    _require(overrides == expected_overrides,
             "plan command_overrides must contain only the exact limits, divB "
             "cadence, optional restart tightening, and optional capacity increase")

    transport = contract.get("input_transport")
    _require(isinstance(transport, dict) and
             transport.get("kind") == INPUT_TRANSPORT_KIND and
             transport.get("holder_pid_token") == HOLDER_PID_TOKEN,
             f"launch input transport must be {INPUT_TRANSPORT_KIND}")
    roles = transport.get("roles")
    _require(isinstance(roles, dict) and
             set(roles) == {"source_restart", "trajectory"},
             "input transport roles must be exactly source_restart/trajectory")
    expected_role_contracts = {
        "source_restart": (SOURCE_RESTART_FD, inputs["source_restart"]),
        "trajectory": (TRAJECTORY_FD, inputs["trajectory"]),
    }
    for role, (descriptor, parent) in expected_role_contracts.items():
        value = roles.get(role)
        _require(isinstance(value, dict) and value.get("role") == role and
                 value.get("fd") == descriptor and
                 value.get("proc_path_template") ==
                 f"/proc/{HOLDER_PID_TOKEN}/fd/{descriptor}",
                 f"input transport role {role} is not canonical")
        staged_file = _validate_file_record(
            value.get("staged_file"), f"input_transport.roles.{role}.staged_file")
        _require(staged_file["size"] == parent["size"] and
                 staged_file["sha256"] == parent["sha256"],
                 f"input transport role {role} does not bind parent content identity")

    source_parameters = plan.get("source", {}).get("parameters", {})
    serialized_trajectory = source_parameters.get("problem", {}).get(
        "trajectory_file") if isinstance(source_parameters, dict) else None
    _require(isinstance(serialized_trajectory, str) and serialized_trajectory,
             "plan source lacks serialized problem/trajectory_file identity")
    parent_identity = transport.get("parent_content_identity")
    _require(parent_identity == {
        "source_restart_sha256": inputs["source_restart"]["sha256"],
        "trajectory_sha256": inputs["trajectory"]["sha256"],
        "source_serialized_trajectory_path": serialized_trajectory,
    }, "input transport parent content identity is not exact")
    trajectory_rebinding = transport.get("trajectory_rebinding")
    trajectory_template = f"/proc/{HOLDER_PID_TOKEN}/fd/{TRAJECTORY_FD}"
    _require(trajectory_rebinding == {
        "parameter": "problem/trajectory_file",
        "parent_sha256": inputs["trajectory"]["sha256"],
        "runtime_value_template": trajectory_template,
    }, "trajectory rebinding authorization is not exact")
    staging_value = transport.get("staging_dir")
    _require(isinstance(staging_value, str) and Path(staging_value).is_absolute(),
             "input_transport.staging_dir must be an absolute path")
    disk_contract = _validate_disk_preflight_contract(
        contract.get("disk_preflight"), inputs["source_restart"]["size"],
        inputs["trajectory"]["size"])
    resolved_staging = _safe_directory(
        Path(staging_value), "staging directory")
    _require(str(resolved_staging) == staging_value,
             "input_transport.staging_dir must be canonical")

    source_template = f"/proc/{HOLDER_PID_TOKEN}/fd/{SOURCE_RESTART_FD}"
    state_template = f"/proc/{HOLDER_PID_TOKEN}/fd/{STATE_DIRECTORY_FD}"
    binary_template = f"/proc/{HOLDER_PID_TOKEN}/fd/{BINARY_EXECUTABLE_FD}"
    divb_cadence_override = f"output3/dt={root_dt!r}"
    canonical_athena_template = (
        binary_template, KOKKOS_GPU_MAP_TOKEN,
        "-r", source_template,
        "-d", state_template, "-t", _wall_time_token(wall_time),
        f"time/nlim={final_cycle}", f"time/tlim={tlim!r}",
        f"problem/trajectory_file={trajectory_template}",
        divb_cadence_override,
        *(() if restart_cadence_transition["kind"] == "unchanged_v1" else
          (restart_cadence_transition["runtime_override"],)),
        *(() if capacity_transition["kind"] == "unchanged_v1" else
          (capacity_transition["runtime_override"],)),
    )
    canonical_mpi = (
        launcher_plan["path"], "--allow-run-as-root", "--bind-to", "none",
        "-np", str(world_size),
    )
    _require(contract.get("athena_argv_template") == list(canonical_athena_template),
             "launch_contract.athena_argv_template is not canonical")
    _require("athena_argv" not in contract,
             "plan must bind a holder-PID template, not pathname-only Athena argv")
    _require(contract.get("mpi_argv") == list(canonical_mpi),
             "launch_contract.mpi_argv is not canonical")
    permitted_trajectory = f"problem/trajectory_file={trajectory_template}"
    permitted_overrides = {permitted_trajectory, divb_cadence_override}
    if restart_cadence_transition["kind"] == "tighten_v1":
        permitted_overrides.add(restart_cadence_transition["runtime_override"])
    if capacity_transition["kind"] == "increase_v1":
        permitted_overrides.add(capacity_transition["runtime_override"])
    _require(not any(FORBIDDEN_OVERRIDE.match(token) and token not in permitted_overrides
                     for token in canonical_athena_template),
             "canonical Athena argv contains a forbidden output/provenance override")
    _require(canonical_athena_template.count(permitted_trajectory) == 1 and
             canonical_athena_template.count(divb_cadence_override) == 1,
             "canonical Athena argv must contain exact trajectory/divB overrides")
    mesh_tokens = [token for token in canonical_athena_template
                   if token.startswith(("mesh/", "meshblock/",
                                        "mesh_refinement/"))]
    _require(mesh_tokens == ([] if capacity_transition["kind"] == "unchanged_v1"
                             else [capacity_transition["runtime_override"]]),
             "canonical Athena argv contains an extra Mesh parameter")

    gpus = query_gpu_inventory(runtime)
    _require(len(gpus) == world_size,
             f"planned {world_size} ranks require exactly {world_size} GPUs; "
             f"found {len(gpus)}")
    _require([gpu.cuda_ordinal for gpu in gpus] == list(range(world_size)),
             "PCI-ordered CUDA ordinals must be contiguous from zero")
    _require(all(gpu.uncorrected_ecc == 0 and gpu.corrected_ecc == 0 for gpu in gpus),
             "preflight GPU ECC counters must all be zero")
    gpu_capacity_preflight = _gpu_capacity_preflight(
        gpus, capacity_transition)
    active = query_gpu_applications(runtime)
    _require(not active, f"GPU compute applications already active: {active!r}")
    _paths_separate(
        resolved_state, resolved_evidence,
        "state and launch evidence directories must not contain each other")
    evidence_entries = list(resolved_evidence.iterdir())
    allowed_evidence_entries = {
        Path(plan_binding["path"])
    } if Path(plan_binding["path"]).parent == resolved_evidence else set()
    _require(set(evidence_entries) <= allowed_evidence_entries,
             "evidence directory must initially contain only the immutable plan")

    state_path, state_fd, state_record = _open_bound_directory(
        resolved_state, "state directory")
    staging_fd: int | None = None
    try:
        staging_dir, staging_fd = _staging_directory(
            resolved_staging, state_path)
        state_info = os.fstat(state_fd)
        staging_info = os.fstat(staging_fd)
        disk_preflight_before_staging = _disk_preflight_snapshot(
            disk_contract,
            "before_staging",
            {
                "state_dir": {
                    "device": state_info.st_dev, "target": state_fd,
                    "access": {
                        "method": "bound_directory_fstatvfs_v1",
                        "fd": state_fd, "planned_path": str(state_path),
                    },
                },
                "staging_dir": {
                    "device": staging_info.st_dev, "target": staging_fd,
                    "access": {
                        "method": "bound_directory_fstatvfs_v1",
                        "fd": staging_fd, "planned_path": str(staging_dir),
                    },
                },
            },
            runtime.fstatvfs,
        )
    except BaseException:
        os.close(state_fd)
        if staging_fd is not None and _fd_is_open(staging_fd):
            os.close(staging_fd)
        raise
    evidence_path: Path | None = None
    evidence_fd: int | None = None
    directory_holder: HolderDirectories | None = None
    launcher_fd: int | None = None
    binary_fd: int | None = None
    executable_holder: HolderExecutables | None = None
    staged_source_fd: int | None = None
    staged_trajectory_fd: int | None = None
    staged_source: dict[str, Any] | None = None
    staged_trajectory: dict[str, Any] | None = None
    staging_record: dict[str, Any] | None = None
    holder: HolderInputs | None = None
    try:
        evidence_path, evidence_fd, evidence_record = _open_bound_directory(
            resolved_evidence, "evidence directory")
        directory_holder = _install_holder_directories(
            state_fd, evidence_fd, state_record, evidence_record,
            runtime.inspector)
        launcher_fd = _open_executable_descriptor(launcher, "MPI launcher")
        binary_fd = _open_executable_descriptor(binary, "Athena binary")
        executable_holder = _install_holder_executables(
            launcher_fd, binary_fd, launcher, binary, runtime.inspector)
        _paths_separate(
            evidence_path, staging_dir,
            "staging and launch evidence directories must not contain each other")
        staged_source, staged_source_fd = _stage_from_audited_descriptor(
            inputs["source_restart"], roles["source_restart"]["staged_file"],
            staging_dir, staging_fd, "source_restart")
        staged_trajectory, staged_trajectory_fd = _stage_from_audited_descriptor(
            inputs["trajectory"], roles["trajectory"]["staged_file"],
            staging_dir, staging_fd, "trajectory")
        staging_record = _seal_staging_directory(staging_dir, staging_fd)
        holder = _install_holder_inputs(
            staged_source_fd, staged_trajectory_fd,
            staged_source, staged_trajectory, runtime.inspector)
    except BaseException:
        if holder is not None:
            holder.close()
        if directory_holder is not None:
            directory_holder.close()
        if executable_holder is not None:
            executable_holder.close()
        if staging_fd is not None:
            created = [record for record in (staged_source, staged_trajectory)
                       if record is not None]
            if created:
                bound = os.fstat(staging_fd)
                _cleanup_bound_staged_inputs(
                    staging_fd,
                    {"path": str(resolved_staging), "device": bound.st_dev,
                     "inode": bound.st_ino},
                    created,
                )
            os.close(staging_fd)
            staging_fd = None
        raise
    finally:
        for descriptor in (
                staged_source_fd, staged_trajectory_fd,
                state_fd, evidence_fd, launcher_fd, binary_fd):
            if descriptor is not None:
                try:
                    os.close(descriptor)
                except OSError:
                    pass
    _require(holder is not None and directory_holder is not None and
             executable_holder is not None and staging_fd is not None and
             staging_record is not None and staged_source is not None and
             staged_trajectory is not None,
             "failed to install holder descriptors")
    try:
        _require(holder.pid == directory_holder.pid,
                 "input and directory holders must share one process identity")
        _require(holder.pid == executable_holder.pid,
                 "all holder transports must share one process identity")
        proc_access_probe = _probe_proc_holder_access(
            holder, directory_holder, executable_holder)
        canonical_athena = tuple(
            token.replace(HOLDER_PID_TOKEN, str(holder.pid))
            for token in canonical_athena_template
        )
        return PreparedLaunch(
            plan=plan, plan_binding=plan_binding, state_dir=resolved_state,
            world_size=world_size, wall_time_seconds=wall_time,
            plan_validation=plan_validation,
            source_validation=source_validation,
            proc_access_probe=proc_access_probe,
            launch_environment=launch_environment,
            launch_environment_sha256=launch_environment_sha256,
            mca_configuration_preflight=mca_configuration_preflight,
            execution_tools=execution_tools,
            repository=repository, binary=binary, source_restart=source,
            trajectory=trajectory, staging_directory=staging_record,
            staged_source_restart=staged_source, staged_trajectory=staged_trajectory,
            staging_directory_fd=staging_fd,
            disk_preflight_before_staging=disk_preflight_before_staging,
            input_holder=holder, directory_holder=directory_holder,
            executable_holder=executable_holder,
            launcher=launcher,
            mpi_argv=canonical_mpi, athena_argv=canonical_athena,
            launch_argv=canonical_mpi + canonical_athena, gpus=tuple(gpus),
            gpu_capacity_preflight=gpu_capacity_preflight,
            gpu_visibility_environment={
                "CUDA_VISIBLE_DEVICES": None,
                "KOKKOS_VISIBLE_DEVICES": None,
                "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
            },
        )
    except BaseException:
        holder.close()
        directory_holder.close()
        executable_holder.close()
        if _fd_is_open(staging_fd):
            os.close(staging_fd)
        raise


def _same_executable(actual: Mapping[str, Any],
                     expected: Mapping[str, Any]) -> bool:
    return all(actual.get(name) == expected.get(name) for name in
               ("device", "inode", "size", "mtime_ns", "ctime_ns", "sha256"))


def _parse_rank_environment(environment: Mapping[str, str], world_size: int,
                            pid: int, launch_environment: Mapping[str, str],
                            *, launcher_pid: int, hostname: str,
                            state_dir: Path,
                            athena_argv: Sequence[str]) \
        -> tuple[int, int, dict[str, str]]:
    _require(all(isinstance(key, str) and isinstance(item, str)
                 for key, item in environment.items()),
             f"rank pid {pid} environment is not a string mapping")
    actual_keys = set(environment)
    expected_keys = set(RANK_ENVIRONMENT_KEYS)
    _require(actual_keys == expected_keys,
             f"rank pid {pid} runtime environment key closure differs: "
             f"missing={sorted(expected_keys - actual_keys)!r}, "
             f"unexpected={sorted(actual_keys - expected_keys)!r}")
    selected = {key: environment[key] for key in RANK_ENVIRONMENT_KEYS}
    try:
        global_rank = int(selected["OMPI_COMM_WORLD_RANK"])
        global_size = int(selected["OMPI_COMM_WORLD_SIZE"])
        local_rank = int(selected["OMPI_COMM_WORLD_LOCAL_RANK"])
        local_size = int(selected["OMPI_COMM_WORLD_LOCAL_SIZE"])
    except ValueError as exc:
        raise LaunchFailure(f"rank pid {pid} has non-integer MPI environment") from exc
    _require(global_size == world_size and local_size == world_size,
             f"rank pid {pid} is not in the planned single-node MPI world")
    _require(0 <= global_rank < world_size and 0 <= local_rank < world_size,
             f"rank pid {pid} has out-of-range MPI rank")
    namespace_family = f"prterun-{hostname}-{launcher_pid}"
    namespace = f"{namespace_family}@1"
    expected = {
        **{key: launch_environment[key]
           for key in RANK_INHERITED_LAUNCH_ENVIRONMENT_KEYS},
        **RANK_FIXED_ENVIRONMENT_VALUES,
        "OMPI_COMMAND": str(athena_argv[0]),
        "OMPI_ARGV": " ".join(athena_argv[1:]),
        "OMPI_COMM_WORLD_RANK": str(global_rank),
        "OMPI_COMM_WORLD_SIZE": str(world_size),
        "OMPI_COMM_WORLD_LOCAL_RANK": str(local_rank),
        "OMPI_COMM_WORLD_LOCAL_SIZE": str(world_size),
        "OMPI_COMM_WORLD_NODE_RANK": str(local_rank),
        "OMPI_FILE_LOCATION": f"/tmp/ompi.{launcher_pid}/1/{global_rank}",
        "OMPI_MCA_initial_wdir": str(state_dir),
        "OMPI_MCA_num_procs": str(world_size),
        "OMPI_WORLD_LOCAL_SIZE": str(world_size),
        "OMPI_WORLD_SIZE": str(world_size),
        "PMIX_HOSTNAME": hostname,
        "PMIX_NAMESPACE": namespace,
        "PMIX_RANK": str(global_rank),
        "PMIX_SERVER_TMPDIR": f"/tmp/ompi.{launcher_pid}",
        "PWD": str(state_dir),
    }
    for key, expected_value in expected.items():
        _require(selected.get(key) == expected_value,
                 f"rank pid {pid} environment {key} differs from derived contract")
    uri_values = {selected[key] for key in RANK_URI_ENVIRONMENT_KEYS}
    _require(len(uri_values) == 1,
             f"rank pid {pid} PMIx server URI aliases differ")
    uri = next(iter(uri_values))
    uri_match = re.fullmatch(
        rf"{re.escape(namespace_family)}@0\.0;tcp4://127\.0\.0\.1:(\d+)",
        uri)
    _require(uri_match is not None and 1 <= int(uri_match.group(1)) <= 65535,
             f"rank pid {pid} PMIx server URI is not the derived loopback URI")
    return global_rank, local_rank, selected


def prove_running_launch(process: Any, prepared: PreparedLaunch,
                         runtime: Runtime | None = None,
                         timeout_seconds: float = DEFAULT_PROOF_TIMEOUT_SECONDS,
                         poll_seconds: float = DEFAULT_PROOF_POLL_SECONDS,
                         plan_path: Path | None = None) -> dict[str, Any]:
    runtime = runtime or Runtime()
    _require(math.isfinite(timeout_seconds) and timeout_seconds > 0,
             "proof timeout must be finite and positive")
    _require(math.isfinite(poll_seconds) and poll_seconds >= 0,
             "proof polling interval must be finite and nonnegative")
    deadline = runtime.monotonic() + timeout_seconds
    last_error = "rank processes have not appeared"
    while True:
        return_code = process.poll()
        if return_code is not None:
            raise LaunchFailure(
                f"MPI launcher exited with status {return_code} before proof: {last_error}")
        try:
            mpirun_cmdline = runtime.inspector.cmdline(process.pid)
            _require(mpirun_cmdline == list(prepared.launch_argv),
                     "live mpirun argv differs from canonical launch argv")
            mpirun_exe = runtime.inspector.executable(process.pid)
            _require(_same_executable(mpirun_exe, prepared.launcher),
                     "live mpirun executable differs from planned launcher")
            mpirun_start_time = runtime.inspector.start_time_ticks(process.pid)
            _require(isinstance(mpirun_start_time, int),
                     "live mpirun process lacks a stable starttime identity")

            candidates: list[tuple[int, dict[str, Any]]] = []
            for pid in sorted(runtime.inspector.descendants(process.pid)):
                executable = runtime.inspector.executable(pid)
                if _same_executable(executable, prepared.binary):
                    candidates.append((pid, executable))
            _require(len(candidates) == prepared.world_size,
                     f"expected exactly {prepared.world_size} live Athena ranks, "
                     f"found {len(candidates)}")

            applications = query_gpu_applications(runtime)
            application_map: dict[int, list[str]] = {}
            for application in applications:
                application_map.setdefault(application.pid, []).append(application.uuid)
            candidate_pids = {pid for pid, _ in candidates}
            _require(set(application_map) == candidate_pids,
                     "GPU compute contexts are not owned exactly by the Athena ranks")
            gpu_by_uuid = {gpu.uuid: gpu for gpu in prepared.gpus}
            hostname = runtime.hostname()
            ranks: list[dict[str, Any]] = []
            for pid, executable in candidates:
                cmdline = runtime.inspector.cmdline(pid)
                _require(cmdline == list(prepared.athena_argv),
                         f"rank pid {pid} argv differs from canonical Athena argv")
                environment = runtime.inspector.environment(pid)
                start_time = runtime.inspector.start_time_ticks(pid)
                _require(isinstance(start_time, int),
                         f"rank pid {pid} lacks a stable starttime identity")
                global_rank, local_rank, selected_environment = \
                    _parse_rank_environment(
                        environment, prepared.world_size, pid,
                        prepared.launch_environment,
                        launcher_pid=process.pid, hostname=hostname,
                        state_dir=prepared.state_dir,
                        athena_argv=prepared.athena_argv)
                _require(all(not selected_environment.get(name)
                             for name in GPU_VISIBILITY_ENVIRONMENT_KEYS),
                         f"rank pid {pid} has a forbidden GPU visibility override")
                _require(selected_environment.get("CUDA_DEVICE_ORDER") == "PCI_BUS_ID",
                         f"rank pid {pid} lacks canonical CUDA PCI ordering")
                uuids = application_map.get(pid, [])
                _require(len(uuids) == 1,
                         f"rank pid {pid} must own exactly one GPU context; got {uuids!r}")
                uuid = uuids[0]
                _require(uuid in gpu_by_uuid,
                         f"rank pid {pid} uses an unplanned GPU UUID {uuid}")
                gpu = gpu_by_uuid[uuid]
                _require(gpu.cuda_ordinal == local_rank,
                         f"rank pid {pid} GPU UUID does not match its PCI-ordered "
                         "CUDA ordinal")
                ranks.append({
                    "global_rank": global_rank,
                    "local_rank": local_rank,
                    "pid": pid,
                    "start_time_ticks": start_time,
                    "hostname": hostname,
                    "gpu_index": gpu.index,
                    "gpu_uuid": gpu.uuid,
                    "gpu_pci_bus_id": gpu.pci_bus_id,
                    "gpu_cuda_ordinal": gpu.cuda_ordinal,
                    "executable": executable,
                    "cmdline": cmdline,
                    "mpi_environment": selected_environment,
                })
            ranks.sort(key=lambda row: row["global_rank"])
            _require([row["global_rank"] for row in ranks] ==
                     list(range(prepared.world_size)),
                     "global MPI ranks are not exactly 0..world_size-1")
            _require(sorted(row["local_rank"] for row in ranks) ==
                     list(range(prepared.world_size)),
                     "local MPI ranks are not exactly 0..world_size-1")
            _require(len({row["pid"] for row in ranks}) == prepared.world_size,
                     "Athena rank PIDs are not unique")
            _require(len({row["gpu_uuid"] for row in ranks}) == prepared.world_size,
                     "Athena ranks are not mapped one-to-one onto GPUs")
            _require(len({row["mpi_environment"]["PMIX_SERVER_URI21"]
                          for row in ranks}) == 1,
                     "Athena ranks do not share one PMIx server URI")
            input_transport = prepared.input_holder.audit()
            directory_transport = prepared.directory_holder.audit()
            executable_transport = prepared.executable_holder.audit()
            repository = _audit_repository(
                prepared.plan["inputs"]["repo"],
                prepared.plan["tools"]["git"], runtime)
            execution_tools_at_launch = {
                name: _audit_planned_file(
                    prepared.plan["tools"][name],
                    f"tools.{name} at launch", executable=True)
                for name in prepared.execution_tools
            }
            _require(all(_same_artifact_binding(
                         prepared.execution_tools[name], binding)
                         for name, binding in execution_tools_at_launch.items()),
                     "plan-bound launcher execution tools changed after preflight")
            mca_configuration_at_rank_proof = _audit_mca_configuration(
                prepared.plan["launch_contract"].get("mca_configuration"),
                prepared.launch_environment["HOME"])
            if plan_path is not None:
                current_plan, current_binding = _read_immutable_plan(plan_path)
                _require(current_plan == prepared.plan and
                         current_binding["sha256"] == prepared.plan_binding["sha256"],
                         "immutable plan changed after preflight")
            _require(process.poll() is None,
                     "MPI launcher exited while launch proof was being closed")
            _require(runtime.inspector.start_time_ticks(process.pid) ==
                     mpirun_start_time and
                     runtime.inspector.cmdline(process.pid) == mpirun_cmdline and
                     _same_executable(runtime.inspector.executable(process.pid),
                                      prepared.launcher),
                     "mpirun identity changed while launch proof was being closed")
            for rank in ranks:
                _require(
                    runtime.inspector.start_time_ticks(rank["pid"]) ==
                    rank["start_time_ticks"] and
                    runtime.inspector.cmdline(rank["pid"]) == rank["cmdline"] and
                    _same_executable(
                        runtime.inspector.executable(rank["pid"]), prepared.binary),
                    f"rank {rank['global_rank']} identity changed while launch proof "
                    "was being closed",
                )
            closing_applications = query_gpu_applications(runtime)
            _require(sorted((item.pid, item.uuid) for item in closing_applications) ==
                     sorted((rank["pid"], rank["gpu_uuid"]) for rank in ranks),
                     "GPU contexts changed while launch proof was being closed")
            prepared.input_holder.audit_metadata()
            return {
                "mpirun_pid": process.pid,
                "mpirun_start_time_ticks": mpirun_start_time,
                "mpirun_executable": mpirun_exe,
                "mpirun_cmdline": mpirun_cmdline,
                "hostname": hostname,
                "ranks": ranks,
                "input_transport": input_transport,
                "directory_transport": directory_transport,
                "executable_transport": executable_transport,
                "repository_at_launch": repository,
                "execution_tools_at_launch": execution_tools_at_launch,
                "mca_configuration_at_rank_proof":
                    mca_configuration_at_rank_proof,
                "gpu_mapping_basis": (
                    "kokkos_mpi_rank_token_plus_ompi_local_rank_plus_"
                    "nvidia_compute_context_uuid"
                ),
            }
        except LaunchFailure as exc:
            last_error = str(exc)
        if runtime.monotonic() >= deadline:
            raise LaunchFailure(f"launch proof timed out: {last_error}")
        runtime.sleep(poll_seconds)


def _install_immutable_bytes(path: Path, content: bytes) -> dict[str, Any]:
    """Install a small evidence file atomically, read-only, without replacement."""

    destination = _validate_new_target(path, "evidence file")
    parent = destination.parent
    temporary = parent / f".{destination.name}.{os.getpid()}.{time.time_ns()}.tmp"
    descriptor: int | None = None
    try:
        descriptor = os.open(
            temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL |
            getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0), 0o600,
        )
        with os.fdopen(descriptor, "wb", closefd=True) as stream:
            descriptor = None
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
            os.fchmod(stream.fileno(), 0o444)
            os.fsync(stream.fileno())
        try:
            os.link(temporary, destination, follow_symlinks=False)
        except FileExistsError as exc:
            raise LaunchFailure(f"refusing to overwrite evidence: {destination}") from exc
    finally:
        if descriptor is not None:
            os.close(descriptor)
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        directory = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    return stable_sha256(destination)


def _gpu_csv(gpus: Sequence[GPURecord]) -> bytes:
    return "".join(gpu.csv_line() for gpu in gpus).encode("ascii")


def _gpu_identity(
        gpus: Sequence[GPURecord]) -> list[tuple[int, str, str, int, int]]:
    return [
        (gpu.index, gpu.uuid, gpu.pci_bus_id.lower(), gpu.cuda_ordinal,
         gpu.memory_total_mib)
        for gpu in gpus
    ]


def prove_process_quiescence(proof: Mapping[str, Any],
                             managed_process_group: Mapping[str, Any],
                             runtime: Runtime,
                             timeout_seconds: float =
                             DEFAULT_QUIESCENCE_TIMEOUT_SECONDS,
                             poll_seconds: float =
                             DEFAULT_QUIESCENCE_POLL_SECONDS) -> dict[str, Any]:
    """Boundedly prove process/group/GPU quiescence after launcher exit."""

    _require(math.isfinite(timeout_seconds) and timeout_seconds > 0.0,
             "quiescence timeout must be finite and positive")
    _require(math.isfinite(poll_seconds) and poll_seconds > 0.0,
             "quiescence poll interval must be finite and positive")

    identities = [{
        "role": "mpirun",
        "pid": proof["mpirun_pid"],
        "recorded_start_time_ticks": proof["mpirun_start_time_ticks"],
    }]
    identities.extend({
        "role": f"rank:{rank['global_rank']}",
        "pid": rank["pid"],
        "recorded_start_time_ticks": rank["start_time_ticks"],
    } for rank in proof["ranks"])
    managed_group = managed_process_group
    _require(isinstance(managed_group, dict) and
             managed_group.get("pgid") == proof["mpirun_pid"],
             "launch proof lacks exact managed process group")
    deadline = runtime.monotonic() + timeout_seconds
    observations: list[dict[str, Any]] = []
    applications: list[dict[str, Any]] = []
    group_exists = True
    while True:
        observations = []
        identities_gone = True
        for identity in identities:
            current = runtime.inspector.start_time_ticks(identity["pid"])
            gone = current != identity["recorded_start_time_ticks"]
            identities_gone = identities_gone and gone
            observations.append({
                **identity,
                "state": (
                    "disappeared" if current is None else
                    "pid_reused" if gone else "still_live"
                ),
                "observed_start_time_ticks": current,
                "original_identity_gone": gone,
            })
        applications = query_gpu_applications(runtime)
        group_exists = _process_group_exists(managed_group["pgid"], runtime)
        if identities_gone and not applications and not group_exists:
            break
        if runtime.monotonic() >= deadline:
            live_roles = [
                row["role"] for row in observations
                if not row["original_identity_gone"]
            ]
            if live_roles:
                raise LaunchFailure(
                    "recorded process identity is still live after bounded "
                    f"quiescence wait: {live_roles!r}")
            if applications:
                raise LaunchFailure(
                    "GPU compute contexts remain after bounded quiescence "
                    f"wait: {applications!r}")
            raise LaunchFailure(
                "managed MPI process group remains after bounded quiescence wait")
        runtime.sleep(min(poll_seconds, max(0.0, deadline - runtime.monotonic())))
    return {
        "gpu_compute_contexts_empty": True,
        "process_identities": observations,
        "all_original_identities_gone": True,
        "managed_process_group": managed_group,
        "managed_process_group_gone": True,
    }


def _open_new_log(path: Path):
    destination = _validate_new_target(path, "run log")
    descriptor = os.open(
        destination, os.O_WRONLY | os.O_CREAT | os.O_EXCL |
        getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0), 0o644,
    )
    return destination, os.fdopen(descriptor, "wb", buffering=0, closefd=True)


def _process_group_exists(pgid: int, runtime: Runtime) -> bool:
    try:
        runtime.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    except PermissionError as exc:
        raise LaunchFailure(
            f"cannot inspect managed process group {pgid}: {exc}") from exc
    return True


def _install_managed_signal_handlers(runtime: Runtime) -> dict[int, Any]:
    """Convert launcher termination signals into the managed cleanup path."""

    previous: dict[int, Any] = {}

    def interrupt(signum: int, _frame: Any) -> None:
        try:
            name = signal.Signals(signum).name
        except ValueError:  # pragma: no cover - only registered signals arrive
            name = str(signum)
        raise LaunchFailure(f"launcher received managed termination signal {name}")

    try:
        for signum in MANAGED_TERMINATION_SIGNALS:
            previous[signum] = runtime.get_signal_handler(signum)
            runtime.set_signal_handler(signum, interrupt)
    except (OSError, ValueError) as exc:
        for signum, handler in previous.items():
            try:
                runtime.set_signal_handler(signum, handler)
            except (OSError, ValueError):
                pass
        raise LaunchFailure(
            f"cannot install managed termination-signal handlers: {exc}") from exc
    return previous


def _restore_signal_handlers(runtime: Runtime,
                             previous: Mapping[int, Any]) -> None:
    failures: list[str] = []
    for signum, handler in previous.items():
        try:
            runtime.set_signal_handler(signum, handler)
        except (OSError, ValueError) as exc:
            failures.append(f"{signum}:{exc}")
    _require(not failures,
             f"cannot restore launcher signal handlers: {', '.join(failures)}")


def _cleanup_managed_group(process: Any, pgid: int, runtime: Runtime,
                           *, grace_seconds: float = 10.0) -> dict[str, Any]:
    """Terminate the exact new session and prove no GPU/process work remains."""

    _require(isinstance(pgid, int) and pgid == process.pid and pgid > 1,
             "refusing cleanup without an exact launcher-owned process group")
    _require(math.isfinite(grace_seconds) and 0.0 < grace_seconds <= 60.0,
             "cleanup grace period must be in (0, 60] seconds")
    signals_sent: list[str] = []
    for number, label in ((signal.SIGTERM, "SIGTERM"),
                          (signal.SIGKILL, "SIGKILL")):
        if _process_group_exists(pgid, runtime):
            try:
                runtime.killpg(pgid, number)
            except ProcessLookupError:
                pass
            except OSError as exc:
                raise LaunchFailure(
                    f"cannot send {label} to managed process group {pgid}: {exc}") \
                    from exc
            else:
                signals_sent.append(label)
        deadline = time.monotonic() + grace_seconds
        while True:
            # Reap an exited group leader opportunistically.  A zombie leader
            # keeps killpg(..., 0) reporting the group as present until waitpid.
            try:
                process.poll()
            except (OSError, ProcessLookupError):
                pass
            group_live = _process_group_exists(pgid, runtime)
            applications = query_gpu_applications(runtime)
            if not group_live and not applications:
                try:
                    process.wait(timeout=0)
                except (OSError, subprocess.TimeoutExpired):
                    pass
                return {
                    "managed_pgid": pgid, "signals_sent": signals_sent,
                    "process_group_gone": True,
                    "gpu_compute_contexts_empty": True,
                }
            if time.monotonic() >= deadline:
                break
            runtime.sleep(min(0.05, max(0.0, deadline - time.monotonic())))
    applications = query_gpu_applications(runtime)
    raise LaunchFailure(
        f"managed process group {pgid} did not quiesce after SIGKILL; "
        f"GPU applications={applications!r}")


@dataclass(frozen=True)
class OutputPaths:
    launch_record: Path
    completion_record: Path
    run_log: Path
    exit_status: Path
    gpu_before: Path
    gpu_after: Path


def _output_mapping(outputs: OutputPaths) -> dict[str, Path]:
    return {name: Path(getattr(outputs, name)) for name in EVIDENCE_NAMES}


def _evidence_io_path(prepared: PreparedLaunch, path: Path) -> Path:
    role = prepared.directory_holder.roles["evidence_dir"]
    human = path.expanduser().absolute()
    _require(human.parent == Path(role["path"]) and
             human.name not in ("", ".", ".."),
             f"evidence target is outside bound evidence directory: {human}")
    return Path(role["proc_path"]) / human.name


def _humanize_binding(binding: Mapping[str, Any], human: Path) -> dict[str, Any]:
    return {**dict(binding), "path": str(human.expanduser().absolute())}


def _install_evidence_bytes(prepared: PreparedLaunch, human: Path,
                            content: bytes) -> dict[str, Any]:
    binding = _install_immutable_bytes(_evidence_io_path(prepared, human), content)
    return _humanize_binding(binding, human)


def _install_evidence_json(prepared: PreparedLaunch, human: Path,
                           payload: Any) -> dict[str, Any]:
    try:
        encoded = (
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
        binding = _install_immutable_bytes(
            _evidence_io_path(prepared, human), encoded)
    except (OSError, ValueError, RuntimeError) as exc:
        raise LaunchFailure(f"cannot publish evidence JSON: {exc}") from exc
    return _humanize_binding(binding, human)


def _fresh_evidence_binding(prepared: PreparedLaunch, human: Path,
                            label: str) -> dict[str, Any]:
    binding = _fresh_artifact_binding(_evidence_io_path(prepared, human), label)
    return _humanize_binding(binding, human)


def _fresh_artifact_binding(path: Path, label: str) -> dict[str, Any]:
    try:
        return stable_sha256(path)
    except (OSError, ValueError, RuntimeError) as exc:
        raise LaunchFailure(f"cannot bind completed {label}: {exc}") from exc


def _same_artifact_binding(left: Mapping[str, Any],
                           right: Mapping[str, Any]) -> bool:
    return all(left.get(name) == right.get(name) for name in (
        "path", "device", "inode", "size", "mtime_ns", "ctime_ns", "sha256",
    ))


def publish_completion_record(
    destination: Path,
    *,
    return_code: int,
    prepared: PreparedLaunch,
    plan_path: Path,
    launch_record: Path,
    run_log: Path,
    exit_status: Path,
    gpu_before: Path,
    gpu_after: Path,
    quiescence: Mapping[str, Any],
    input_transport: Mapping[str, Any],
    expected_artifacts: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Publish the unique lifecycle-complete marker after all evidence is closed."""

    artifacts = {
        "plan": _fresh_artifact_binding(plan_path, "plan"),
        "launch_record": _fresh_evidence_binding(
            prepared, launch_record, "launch record"),
        "run_log": _fresh_evidence_binding(prepared, run_log, "run log"),
        "exit_status": _fresh_evidence_binding(
            prepared, exit_status, "exit status"),
        "gpu_before": _fresh_evidence_binding(
            prepared, gpu_before, "GPU-before snapshot"),
        "gpu_after": _fresh_evidence_binding(
            prepared, gpu_after, "GPU-after snapshot"),
    }
    _require(set(expected_artifacts) == set(artifacts),
             "completion expected-artifact set is incomplete")
    for name, binding in artifacts.items():
        _require(_same_artifact_binding(binding, expected_artifacts[name]),
                 f"{name} changed before completion publication")
    payload = {
        "schema": SCHEMA,
        "kind": "athenak_segment_completion",
        "status": "ready",
        "created_utc": _utc_now(),
        "return_code": return_code,
        "state_dir": str(prepared.state_dir),
        "world_size": prepared.world_size,
        "artifacts": artifacts,
        "quiescence": dict(quiescence),
        "input_transport": dict(input_transport),
        "directory_transport": prepared.directory_holder.audit(),
        "executable_transport": prepared.executable_holder.audit(),
        "publication_contract": {
            "unique_lifecycle_complete_marker": True,
            "published_after_closed_artifacts": [
                "launch_record", "run_log", "exit_status",
                "gpu_before", "gpu_after",
            ],
        },
    }
    return _install_evidence_json(prepared, destination, payload)


def _run_segment_with_holder(
        prepared: PreparedLaunch, plan_path: Path, outputs: OutputPaths,
        runtime: Runtime,
        proof_timeout_seconds: float) -> int:
    _require(math.isfinite(proof_timeout_seconds) and proof_timeout_seconds > 0.0,
             "proof timeout must be finite and positive")
    current_plan, current_plan_binding = _read_immutable_plan(plan_path)
    _require(current_plan == prepared.plan and
             _same_artifact_binding(current_plan_binding, prepared.plan_binding),
             "run plan path/identity differs from prepared immutable plan")
    contract = prepared.plan["launch_contract"]
    _require(contract.get("plan_path") == prepared.plan_binding["path"],
             "launch contract plan_path differs from immutable plan")
    normalized = OutputPaths(*(
        _validate_new_target(value, name) for value, name in (
            (outputs.launch_record, "launch record"),
            (outputs.completion_record, "completion record"),
            (outputs.run_log, "run log"),
            (outputs.exit_status, "exit status"),
            (outputs.gpu_before, "GPU-before snapshot"),
            (outputs.gpu_after, "GPU-after snapshot"),
        )
    ))
    _require(len(set(normalized.__dict__.values())) == 6,
             "all output evidence paths must be distinct")
    normalized_mapping = {
        name: str(path) for name, path in _output_mapping(normalized).items()
    }
    _require(normalized_mapping == contract.get("evidence"),
             "CLI evidence paths differ from immutable launch contract")
    evidence_parents = {path.parent for path in normalized.__dict__.values()}
    _require(len(evidence_parents) == 1,
             "all launch evidence must use one dedicated directory")
    evidence_directory = next(iter(evidence_parents))
    _paths_separate(
        evidence_directory, prepared.state_dir,
        "state and launch evidence directories must not contain each other")
    staging_path = Path(prepared.staging_directory["path"])
    _paths_separate(
        evidence_directory, staging_path,
        "staging and launch evidence directories must not contain each other")
    prepared.directory_holder.audit()
    prepared.executable_holder.audit()
    _require(not os.listdir(prepared.directory_holder.roles["state_dir"]["fd"]),
             "state directory ceased to be empty after preflight")
    for name, initial in prepared.execution_tools.items():
        current = _audit_planned_file(
            prepared.plan["tools"][name], f"tools.{name} before spawn",
            executable=True)
        _require(_same_artifact_binding(initial, current),
                 f"plan-bound tool {name} changed before spawn")
    _require(runtime.nvidia_smi == prepared.execution_tools["nvidia_smi"]["path"],
             "runtime nvidia-smi is no longer plan-bound")
    before = query_gpu_inventory(runtime)
    _require(_gpu_identity(before) == _gpu_identity(prepared.gpus),
             "GPU PCI/UUID inventory changed after preflight")
    _require(all(gpu.uncorrected_ecc == 0 and gpu.corrected_ecc == 0
                 for gpu in before),
             "GPU ECC counters changed after preflight")
    capacity_transition = _capacity_transition_contract(prepared.plan)
    gpu_capacity_preflight = _gpu_capacity_preflight(
        before, capacity_transition)
    _require(not query_gpu_applications(runtime),
             "a GPU compute application appeared after preflight")
    prepared.input_holder.audit()
    prepared.directory_holder.audit()
    disk_contract = _validate_disk_preflight_contract(
        contract.get("disk_preflight"),
        prepared.plan["inputs"]["source_restart"]["size"],
        prepared.plan["inputs"]["trajectory"]["size"],
    )
    state_preflight_fd = prepared.directory_holder.roles["state_dir"]["fd"]
    staging_preflight_fd = prepared.staging_directory_fd
    _require(state_preflight_fd ==
             disk_contract["bound_directory_fds"]["state_dir"] and
             staging_preflight_fd ==
             disk_contract["bound_directory_fds"]["staging_dir"],
             "bound disk-preflight descriptors differ from the plan")
    state_fs_info = os.fstat(state_preflight_fd)
    staging_fs_info = os.fstat(staging_preflight_fd)
    try:
        disk_preflight_before_spawn = _disk_preflight_snapshot(
            disk_contract,
            "before_spawn",
            {
                "state_dir": {
                    "device": state_fs_info.st_dev,
                    "target": state_preflight_fd,
                    "access": {
                        "method": "bound_fd_fstatvfs_v1",
                        "fd": state_preflight_fd,
                    },
                },
                "staging_dir": {
                    "device": staging_fs_info.st_dev,
                    "target": staging_preflight_fd,
                    "access": {
                        "method": "bound_fd_fstatvfs_v1",
                        "fd": staging_preflight_fd,
                    },
                },
            },
            runtime.fstatvfs,
        )
    except BaseException:
        _cleanup_bound_staged_inputs(
            staging_preflight_fd,
            prepared.staging_directory,
            (prepared.staged_source_restart, prepared.staged_trajectory),
        )
        raise
    before_binding = _install_evidence_bytes(
        prepared,
        normalized.gpu_before, _gpu_csv(before))
    _, log_stream = _open_new_log(_evidence_io_path(prepared, normalized.run_log))
    process: Any | None = None
    managed_pgid: int | None = None
    previous_signal_handlers: dict[int, Any] = {}
    try:
        previous_signal_handlers = _install_managed_signal_handlers(runtime)
        child_environment = dict(prepared.launch_environment)
        _require(_canonical_sha256(child_environment) ==
                 prepared.launch_environment_sha256,
                 "prepared launch environment identity changed")
        mca_configuration_before_spawn = _audit_mca_configuration(
            prepared.plan["launch_contract"].get("mca_configuration"),
            child_environment["HOME"])
        prepared.directory_holder.audit()
        prepared.executable_holder.audit()
        process = runtime.popen(
            list(prepared.launch_argv),
            executable=prepared.executable_holder.roles["launcher"]["proc_path"],
            cwd=prepared.directory_holder.roles["state_dir"]["proc_path"],
            stdin=subprocess.DEVNULL, stdout=log_stream, stderr=subprocess.STDOUT,
            start_new_session=True, close_fds=True, shell=False,
            env=child_environment,
        )
        # start_new_session=True makes the child both session leader and process
        # group leader. Bind that cleanup identity immediately: even if the
        # subsequent getpgid observation fails, every exceptional path can still
        # terminate the whole group before any holder descriptor is closed.
        managed_pgid = int(process.pid)
        _require(managed_pgid > 1,
                 "MPI launcher returned an unsafe managed process-group identity")
        try:
            observed_pgid = runtime.getpgid(process.pid)
        except (OSError, ProcessLookupError) as exc:
            raise LaunchFailure(
                f"cannot bind new MPI process group identity: {exc}") from exc
        _require(observed_pgid == managed_pgid,
                 "MPI launcher did not create its own managed session/process group")
        proof = prove_running_launch(
            process, prepared, runtime, timeout_seconds=proof_timeout_seconds,
            plan_path=plan_path,
        )
        record = {
            "schema": SCHEMA,
            "kind": "athenak_segment_launch",
            "status": "ready",
            "created_utc": _utc_now(),
            "plan_path": str(plan_path.expanduser().resolve(strict=True)),
            "plan_sha256": prepared.plan_binding["sha256"],
            "state_dir": str(prepared.state_dir),
            "world_size": prepared.world_size,
            "gpu_count": len(prepared.gpus),
            "binary_sha256": prepared.plan["inputs"]["binary"]["sha256"],
            "source_restart_sha256": (
                prepared.plan["inputs"]["source_restart"]["sha256"]
            ),
            "trajectory_sha256": prepared.plan["inputs"]["trajectory"]["sha256"],
            "original_inputs": {
                "source_restart": prepared.source_restart,
                "trajectory": prepared.trajectory,
            },
            "staging_directory": prepared.staging_directory,
            "staged_inputs": {
                "source_restart": prepared.staged_source_restart,
                "trajectory": prepared.staged_trajectory,
            },
            "disk_preflight": {
                "contract": disk_contract,
                "before_staging": prepared.disk_preflight_before_staging,
                "before_spawn": disk_preflight_before_spawn,
            },
            "input_transport_contract": (
                prepared.plan["launch_contract"]["input_transport"]
            ),
            "directory_transport_contract": (
                prepared.plan["launch_contract"]["directory_transport"]
            ),
            "executable_transport_contract": (
                prepared.plan["launch_contract"]["executable_transport"]
            ),
            "repository_preflight": prepared.repository,
            "plan_validation": prepared.plan_validation,
            "source_validation": prepared.source_validation,
            "proc_access_probe": prepared.proc_access_probe,
            "launcher": prepared.plan["launch_contract"]["launcher"],
            "mpi_argv": list(prepared.mpi_argv),
            "athena_argv": list(prepared.athena_argv),
            "launch_argv": list(prepared.launch_argv),
            "launch_environment":
                prepared.plan["launch_contract"]["environment"],
            "mca_configuration_contract":
                prepared.plan["launch_contract"]["mca_configuration"],
            "mca_configuration": {
                "preflight": prepared.mca_configuration_preflight,
                "before_spawn": mca_configuration_before_spawn,
                "at_rank_proof": proof["mca_configuration_at_rank_proof"],
            },
            "gpu_visibility_environment": {
                name: child_environment.get(name)
                for name in (*GPU_VISIBILITY_ENVIRONMENT_KEYS, "CUDA_DEVICE_ORDER")
            },
            "gpu_before": {
                "path": str(normalized.gpu_before),
                "sha256": before_binding["sha256"],
                "devices": [gpu.as_dict() for gpu in before],
            },
            "capacity_transition": capacity_transition,
            "gpu_capacity_preflight": gpu_capacity_preflight,
            "launcher_tool": proof["execution_tools_at_launch"][
                "segment_launcher"],
            "nvidia_smi": proof["execution_tools_at_launch"]["nvidia_smi"],
            "managed_process_group": {
                "pgid": managed_pgid,
                "new_session": True,
                "failure_cleanup": "SIGTERM_then_SIGKILL_with_quiescence_proof",
            },
            **proof,
        }
        launch_binding = _install_evidence_json(
            prepared, normalized.launch_record, record)
        return_code = int(process.wait())
        try:
            log_stream.flush()
            os.fsync(log_stream.fileno())
            os.fchmod(log_stream.fileno(), 0o444)
            os.fsync(log_stream.fileno())
        finally:
            log_stream.close()
        quiescence = prove_process_quiescence(
            proof, record["managed_process_group"], runtime)
        after = query_gpu_inventory(runtime)
        _require(len(after) == prepared.world_size,
                 "GPU inventory changed before exit evidence capture")
        _require(_gpu_identity(after) == _gpu_identity(prepared.gpus),
                 "GPU PCI/UUID inventory changed during segment")
        prepared.directory_holder.audit()
        after_binding = _install_evidence_bytes(
            prepared, normalized.gpu_after, _gpu_csv(after))
        exit_binding = _install_evidence_bytes(
            prepared, normalized.exit_status, f"{return_code}\n".encode("ascii"))
        log_binding = _fresh_evidence_binding(
            prepared, normalized.run_log, "run log")
        holder_at_completion = prepared.input_holder.audit()
        prepared.input_holder.close()
        holder_closure = prepared.input_holder.prove_closed()
        completion_transport = {
            "kind": INPUT_TRANSPORT_KIND,
            "at_launch": proof["input_transport"],
            "at_completion": holder_at_completion,
            "closure": holder_closure,
        }
        # This is deliberately the last creation and sole lifecycle marker.
        publish_completion_record(
            normalized.completion_record,
            return_code=return_code,
            prepared=prepared,
            plan_path=plan_path,
            launch_record=normalized.launch_record,
            run_log=normalized.run_log,
            exit_status=normalized.exit_status,
            gpu_before=normalized.gpu_before,
            gpu_after=normalized.gpu_after,
            quiescence=quiescence,
            input_transport=completion_transport,
            expected_artifacts={
                "plan": prepared.plan_binding,
                "launch_record": launch_binding,
                "run_log": log_binding,
                "exit_status": exit_binding,
                "gpu_before": before_binding,
                "gpu_after": after_binding,
            },
        )
        return return_code
    except BaseException as primary:
        secondary_errors: list[str] = []
        if previous_signal_handlers:
            try:
                _restore_signal_handlers(runtime, previous_signal_handlers)
            except BaseException as restore_error:
                secondary_errors.append(f"signal-handler restore failed: {restore_error}")
            else:
                previous_signal_handlers = {}
        try:
            if process is not None and managed_pgid is not None:
                _cleanup_managed_group(process, managed_pgid, runtime)
        except BaseException as cleanup_error:
            secondary_errors.append(f"managed process cleanup failed: {cleanup_error}")
        if secondary_errors:
            raise LaunchFailure(
                f"segment failed ({primary}); " + "; ".join(secondary_errors)
            ) from primary
        raise
    finally:
        if previous_signal_handlers:
            _restore_signal_handlers(runtime, previous_signal_handlers)
        if not log_stream.closed:
            try:
                log_stream.flush()
                os.fsync(log_stream.fileno())
                os.fchmod(log_stream.fileno(), 0o444)
                os.fsync(log_stream.fileno())
            finally:
                log_stream.close()


def run_segment(prepared: PreparedLaunch, plan_path: Path, outputs: OutputPaths,
                runtime: Runtime | None = None,
                proof_timeout_seconds: float = DEFAULT_PROOF_TIMEOUT_SECONDS) -> int:
    """Run one segment and close its holder descriptors on every return path."""

    runtime = runtime or Runtime()
    try:
        return _run_segment_with_holder(
            prepared, plan_path, outputs, runtime, proof_timeout_seconds)
    finally:
        prepared.close()


def _prepared_summary(prepared: PreparedLaunch) -> dict[str, Any]:
    capacity_transition = _capacity_transition_contract(prepared.plan)
    return {
        "status": "prepared",
        "plan_path": prepared.plan_binding["path"],
        "plan_sha256": prepared.plan_binding["sha256"],
        "state_dir": str(prepared.state_dir),
        "world_size": prepared.world_size,
        "wall_time_seconds": prepared.wall_time_seconds,
        "mpi_argv": list(prepared.mpi_argv),
        "athena_argv": list(prepared.athena_argv),
        "launch_argv": list(prepared.launch_argv),
        "staging_directory": prepared.staging_directory,
        "staged_inputs": {
            "source_restart": prepared.staged_source_restart,
            "trajectory": prepared.staged_trajectory,
        },
        "disk_preflight_before_staging": prepared.disk_preflight_before_staging,
        "gpu_visibility_environment": prepared.gpu_visibility_environment,
        "launch_environment": prepared.plan["launch_contract"]["environment"],
        "directory_transport": prepared.directory_holder.audit(),
        "executable_transport": prepared.executable_holder.audit(),
        "proc_access_probe": prepared.proc_access_probe,
        "gpus": [gpu.as_dict() for gpu in prepared.gpus],
        "capacity_transition": capacity_transition,
        "gpu_capacity_preflight": prepared.gpu_capacity_preflight,
    }


def main(argv: Sequence[str] | None = None, runtime: Runtime | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--state-dir", type=Path, required=True)
    parser.add_argument("--launch-record", type=Path)
    parser.add_argument("--completion-record", type=Path)
    parser.add_argument("--run-log", type=Path)
    parser.add_argument("--exit-status", type=Path)
    parser.add_argument("--gpu-before", type=Path)
    parser.add_argument("--gpu-after", type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--proof-timeout-seconds", type=float,
                        default=DEFAULT_PROOF_TIMEOUT_SECONDS)
    args = parser.parse_args(argv)
    runtime = runtime or Runtime()
    names = EVIDENCE_NAMES
    if not args.prepare_only:
        missing = [f"--{name.replace('_', '-')}" for name in names
                   if getattr(args, name) is None]
        if missing:
            parser.error(f"non-prepare launch requires {', '.join(missing)}")
    if (not math.isfinite(args.proof_timeout_seconds) or
            args.proof_timeout_seconds <= 0.0):
        parser.error("proof timeout must be finite and positive")
    prepared: PreparedLaunch | None = None
    try:
        prepared = prepare_launch(args.plan, args.state_dir, runtime)
        if args.prepare_only:
            try:
                print(json.dumps(_prepared_summary(prepared), sort_keys=True,
                                 allow_nan=False))
            finally:
                prepared.close()
            return 0
        return run_segment(
            prepared, args.plan,
            OutputPaths(*(getattr(args, name) for name in names)),
            runtime, args.proof_timeout_seconds,
        )
    except (OSError, LaunchFailure, ValueError) as exc:
        if prepared is not None:
            prepared.close()
        parser.error(str(exc))
    return 2  # pragma: no cover


if __name__ == "__main__":
    raise SystemExit(main())
