#!/usr/bin/env python3
"""Fail-closed qualification of one completed AthenaK campaign segment.

The checker deliberately publishes no partial or negative result.  A successful
qualification is atomically created as an immutable ``*.pass.ready`` JSON report;
recent files return a retry status so a closed-file hand-off cannot race final I/O.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import hashlib
import json
import math
import os
from pathlib import Path
import pwd
import re
import stat
import struct
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

sys.dont_write_bytecode = True
SCRIPT_DIRECTORY = Path(__file__).resolve().parent
if str(SCRIPT_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIRECTORY))

from audit_athenak_restart import audit_restart
from output_integrity import (
    HASH_CHUNK_BYTES,
    _assert_closed,
    _assert_path_signature,
    _assert_stream_signature,
    _open_regular_nofollow,
    install_immutable_json,
    stable_sha256,
    strict_json_loads,
)
from read_athenak_restart_metadata import (
    RestartMetadata,
    RestartTruncationError,
    load_balance,
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
PREFIX_RECOVERY_POLICY = {
    "kind": "scheduled_prefix_recovery_v1",
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
MIN_READY_AGE_SECONDS = 120.0
RETRY_EXIT_STATUS = 75
EVENT_COLUMNS = (
    "cycle",
    "eos_dfloor",
    "eos_efloor",
    "eos_tfloor",
    "eos_vceil",
    "eos_fail",
    "c2p_it",
    "fofc",
    "cons_adjust",
    "mag_adjust",
    "c2p_calls",
    "fofc_tests",
)
GPU_HEADER = (
    "index", "uuid", "pci_bus_id", "cuda_ordinal", "uncorr", "corr",
    "memory_total", "memory_used",
)
PCI_BUS_ID_RE = re.compile(
    r"^(?:[0-9A-Fa-f]{4}|[0-9A-Fa-f]{8}):"
    r"[0-9A-Fa-f]{2}:[0-9A-Fa-f]{2}\.[0-7]$"
)


def _pci_bus_sort_key(value: str) -> tuple[int, int, int, int]:
    match = PCI_BUS_ID_RE.fullmatch(value)
    if match is None:
        raise CheckFailure(f"invalid GPU PCI bus identity: {value!r}")
    domain_text, bus_text, device_text, function_text = re.split(r"[:.]", value)
    return tuple(int(part, 16) for part in (
        domain_text, bus_text, device_text, function_text))
PLANNED_TOOL_NAMES = (
    "planner",
    "segment_checker",
    "segment_launcher",
    "restart_metadata_reader",
    "output_integrity",
    "restart_auditor",
    "restart_layout",
    "prefix_recovery",
    "git",
    "nvidia_smi",
)
CANONICAL_MUTABLE_PARAMETERS = [
    "time/nlim",
    "time/tlim",
    "time/restart_dt_growth",
    "output*/file_number",
    "output*/last_time",
    "output*/last_write_cycle",
]
INPUT_TRANSPORT_KIND = "linux_proc_holder_fd_v1"
HOLDER_PID_TOKEN = "{holder_pid}"
SOURCE_RESTART_FD = 200
TRAJECTORY_FD = 201
STATE_DIRECTORY_FD = 202
EVIDENCE_DIRECTORY_FD = 203
LAUNCHER_EXECUTABLE_FD = 204
BINARY_EXECUTABLE_FD = 205
DIRECTORY_TRANSPORT_KIND = "linux_proc_holder_dirfd_v1"
EXECUTABLE_TRANSPORT_KIND = "linux_proc_holder_execfd_v1"
LAUNCH_ENVIRONMENT_KIND = "explicit_values_with_rank_projection_v3"
RANK_ENVIRONMENT_PROJECTION_KIND = \
    "prrte_openmpi_pmix_single_node_projection_v2"
MCA_CONFIGURATION_KIND = "openmpi_prrte_pmix_default_files_v1"
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
GIT_ENVIRONMENT = {
    "LANG": "C",
    "LC_ALL": "C",
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_CONFIG_GLOBAL": "/dev/null",
}
GIT_CONFIGURATION = [
    "-c", "core.fsmonitor=false",
    "-c", "core.untrackedCache=false",
    "-c", "core.hooksPath=/dev/null",
]
FLOAT_TOKEN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
CACHE_RE = re.compile(
    rf"^Strict-subcycling restart cache reconstruction passed: "
    rf"solver failures=(\d+), non-finite proposed values=(\d+), "
    rf"maximum raw component-relative proposed change=({FLOAT_TOKEN}), "
    rf"maximum absolute proposed change=({FLOAT_TOKEN}), "
    rf"maximum mixed-scale proposed change=({FLOAT_TOKEN}), "
    rf"mixed-scale acceptance tolerance=({FLOAT_TOKEN})\.$"
)
FINAL_STATE_RE = re.compile(
    rf"^time=({FLOAT_TOKEN})\s+cycle=(-?\d+)\s*$"
)
LIMIT_STATE_RE = re.compile(
    rf"^tlim=({FLOAT_TOKEN})\s+nlim=(-?\d+)\s*$"
)
TERMINATION_RE = re.compile(r"^Terminating on (cycle|time|wall clock) limit\s*$")
DIAGNOSTIC_RE = re.compile(
    rf"^elapsed=({FLOAT_TOKEN})\s+cycle=(-?\d+)\s+"
    rf"time=({FLOAT_TOKEN})\s+dt=({FLOAT_TOKEN})\s*$"
)


class CheckFailure(RuntimeError):
    """A permanent, structured qualification failure."""


class NotReady(RuntimeError):
    """A transient closed-file readiness failure."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise CheckFailure(message)


def _number(value: Any, name: str) -> float:
    _require(isinstance(value, (int, float)) and not isinstance(value, bool),
             f"{name} must be a number")
    result = float(value)
    _require(math.isfinite(result), f"{name} must be finite")
    return result


def _integer(value: Any, name: str) -> int:
    _require(isinstance(value, int) and not isinstance(value, bool),
             f"{name} must be an integer")
    return value


def _ceil_ratio(numerator: int, denominator: int) -> int:
    return (numerator + denominator - 1) // denominator


def _validate_capacity_transition(plan: dict[str, Any]) -> dict[str, Any]:
    """Validate the exact, single-parameter capacity-transition contract."""

    value = plan.get("capacity_transition")
    _require(isinstance(value, dict) and set(value) == {
        "kind", "parameter", "source_max_nmb_per_rank",
        "target_max_nmb_per_rank", "maximum_target_max_nmb_per_rank",
        "runtime_override", "gpu_memory_model",
    }, "plan capacity_transition must be the exact strict object")
    source = _integer(
        value.get("source_max_nmb_per_rank"),
        "capacity_transition.source_max_nmb_per_rank")
    target = _integer(
        value.get("target_max_nmb_per_rank"),
        "capacity_transition.target_max_nmb_per_rank")
    _require(1 <= source <= target <= MAX_NMB_PER_RANK and
             value.get("maximum_target_max_nmb_per_rank") == MAX_NMB_PER_RANK and
             value.get("parameter") ==
             "mesh_refinement/max_nmb_per_rank",
             "capacity transition source/target/bound is invalid")
    if source == target:
        _require(value.get("kind") == "unchanged_v1" and
                 value.get("runtime_override") is None,
                 "unchanged capacity transition may not carry an override")
    else:
        _require(value.get("kind") == "increase_v1" and
                 value.get("runtime_override") ==
                 f"mesh_refinement/max_nmb_per_rank={target}",
                 "capacity increase lacks its sole canonical runtime override")

    model = value.get("gpu_memory_model")
    required_mib = _ceil_ratio(
        target * CAPACITY_MEMORY_MIB_PER_SLOT_NUMERATOR,
        CAPACITY_MEMORY_MIB_PER_SLOT_DENOMINATOR)
    minimum_total_mib = _ceil_ratio(
        target * CAPACITY_MEMORY_MIB_PER_SLOT_NUMERATOR *
        CAPACITY_GPU_USABLE_FRACTION_DENOMINATOR,
        CAPACITY_MEMORY_MIB_PER_SLOT_DENOMINATOR *
        CAPACITY_GPU_USABLE_FRACTION_NUMERATOR)
    _require(model == {
        "kind": "fixed_conservative_per_meshblock_slot_v1",
        "mib_per_slot_numerator":
            CAPACITY_MEMORY_MIB_PER_SLOT_NUMERATOR,
        "mib_per_slot_denominator":
            CAPACITY_MEMORY_MIB_PER_SLOT_DENOMINATOR,
        "usable_fraction_numerator":
            CAPACITY_GPU_USABLE_FRACTION_NUMERATOR,
        "usable_fraction_denominator":
            CAPACITY_GPU_USABLE_FRACTION_DENOMINATOR,
        "required_per_rank_memory_mib_ceiling": required_mib,
        "minimum_gpu_memory_total_mib": minimum_total_mib,
    }, "capacity transition GPU-memory model is not canonical")
    source_record = plan.get("source")
    _require(isinstance(source_record, dict),
             "plan source must be an object for capacity transition")
    source_parameters = source_record.get("parameters")
    _require(isinstance(source_parameters, dict),
             "plan source parameters must be an object")
    try:
        serialized_source = source_parameters[
            "mesh_refinement"]["max_nmb_per_rank"]
    except (KeyError, TypeError) as exc:
        raise CheckFailure(
            "plan source lacks mesh_refinement/max_nmb_per_rank") from exc
    _require(serialized_source == str(source),
             "plan source capacity is not the canonical transition source")
    return value


def _cadence_root_step_multiple(value: float, root_dt: float,
                                label: str) -> int:
    ratio = value / root_dt
    nearest = round(ratio) if math.isfinite(ratio) else 0
    tolerance = (CADENCE_MULTIPLE_MAX_ULPS *
                 max(math.ulp(ratio), math.ulp(float(nearest)))) \
        if math.isfinite(ratio) else 0.0
    _require(value > 0.0 and math.isfinite(ratio) and nearest >= 1 and
             abs(ratio - nearest) <= tolerance,
             f"{label} must be a positive integer root-dt multiple")
    return nearest


def _validate_restart_cadence_transition(
        plan: dict[str, Any], root_dt: float) -> dict[str, Any]:
    """Validate the exact optional tightening of the sole restart stream."""

    value = plan.get("restart_cadence_transition")
    _require(isinstance(value, dict) and set(value) == {
        "kind", "block", "parameter", "source_dt", "target_dt", "root_dt",
        "source_root_step_multiple", "target_root_step_multiple", "phase",
        "runtime_override",
    }, "plan restart_cadence_transition must be the exact strict object")
    source_record = plan.get("source")
    _require(isinstance(source_record, dict),
             "plan source must bind restart cadence")
    source_parameters = source_record.get("parameters")
    _require(isinstance(source_parameters, dict),
             "plan source parameters must bind restart cadence")
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
             f"source {block_name} restart output must use dt without dcycle")
    try:
        serialized_source_dt = float(block["dt"])
        serialized_phase = {
            "file_number": int(block.get("file_number", "0")),
            "last_time": float(block.get("last_time", "-1")),
            "last_write_cycle": int(block.get("last_write_cycle", "-1")),
        }
    except (TypeError, ValueError) as exc:
        raise CheckFailure(
            f"source {block_name} restart cadence state is invalid") from exc
    _require(serialized_phase["file_number"] >= 0 and
             math.isfinite(serialized_phase["last_time"]),
             f"source {block_name} restart cadence phase is invalid")
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
    source_dt = _number(value.get("source_dt"),
                        "restart_cadence_transition.source_dt")
    target_dt = _number(value.get("target_dt"),
                        "restart_cadence_transition.target_dt")
    transition_root_dt = _number(
        value.get("root_dt"), "restart_cadence_transition.root_dt")
    source_multiple = _cadence_root_step_multiple(
        source_dt, root_dt, "restart cadence source_dt")
    target_multiple = _cadence_root_step_multiple(
        target_dt, root_dt, "restart cadence target_dt")
    _require(value.get("block") == block_name and
             value.get("parameter") == f"{block_name}/dt" and
             transition_root_dt == root_dt and
             serialized_source_dt == source_dt and target_dt <= source_dt and
             value.get("source_root_step_multiple") == source_multiple and
             value.get("target_root_step_multiple") == target_multiple and
             planned_phase == serialized_phase,
             "restart cadence transition source/target/phase is not canonical")
    if target_dt == source_dt:
        _require(value.get("kind") == "unchanged_v1" and
                 value.get("runtime_override") is None,
                 "unchanged restart cadence may not carry an override")
    else:
        _require(target_dt < source_dt and value.get("kind") == "tighten_v1" and
                 value.get("runtime_override") ==
                 f"{block_name}/dt={target_dt!r}",
                 "restart cadence tightening lacks its canonical override")
    return value


def _gpu_capacity_preflight_evidence(
        devices: list[dict[str, Any]],
        transition: dict[str, Any]) -> dict[str, Any]:
    _require(bool(devices), "GPU capacity preflight requires device evidence")
    model = transition["gpu_memory_model"]
    required = int(model["required_per_rank_memory_mib_ceiling"])
    minimum_total = int(model["minimum_gpu_memory_total_mib"])
    _require(all(device["memory_total_mib"] >= minimum_total
                 for device in devices),
             "GPU total memory is below the capacity-transition 80% gate")
    _require(all(device["memory_used_mib"] + required <=
                 device["memory_total_mib"] for device in devices),
             "GPU available memory is below the capacity-transition gate")
    return {
        "kind": "per_rank_memory_total_and_available_gate_v1",
        "required_per_rank_memory_mib_ceiling": required,
        "minimum_gpu_memory_total_mib": minimum_total,
        "minimum_observed_gpu_memory_total_mib":
            min(device["memory_total_mib"] for device in devices),
        "maximum_observed_gpu_memory_used_mib":
            max(device["memory_used_mib"] for device in devices),
        "minimum_observed_gpu_memory_available_mib":
            min(device["memory_total_mib"] - device["memory_used_mib"]
                for device in devices),
        "all_ranks_pass": True,
    }


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
            "staging_dir": 206,
        },
    }
    _require(value == expected,
             "launch disk-preflight contract is not the exact capacity policy")
    return value


def _validate_disk_preflight_snapshot(value: Any, contract: dict[str, Any],
                                      phase: str,
                                      state_dir: str,
                                      staging_dir: str,
                                      expected_devices: dict[str, int]) \
        -> dict[str, Any]:
    """Validate immutable arithmetic evidence without trusting its conclusions."""

    _require(isinstance(value, dict) and set(value) == {
        "kind", "phase", "created_utc", "contract_sha256", "filesystems",
        "status",
    }, f"disk preflight {phase} evidence shape is invalid")
    _require(value.get("kind") == DISK_PREFLIGHT_KIND and
             value.get("phase") == phase and value.get("status") == "pass" and
             isinstance(value.get("created_utc"), str) and value["created_utc"] and
             value.get("contract_sha256") == _canonical_json_sha256(contract),
             f"disk preflight {phase} identity is invalid")
    filesystems = value.get("filesystems")
    _require(isinstance(filesystems, list) and 1 <= len(filesystems) <= 2,
             f"disk preflight {phase} must cover one or two filesystems")
    expected_access = ({
        "state_dir": {"method": "bound_fd_fstatvfs_v1",
                      "fd": contract["bound_directory_fds"]["state_dir"]},
        "staging_dir": {"method": "bound_fd_fstatvfs_v1",
                        "fd": contract["bound_directory_fds"]["staging_dir"]},
    } if phase == "before_spawn" else {})
    expected_contributions = {
        "state_dir": contract["planned_peak_output_bytes"],
        "staging_dir": contract["staging_copy_bytes"],
    }
    seen_roles: set[str] = set()
    seen_devices: set[int] = set()
    reserve = max(
        contract["minimum_reserve_bytes"],
        contract["minimum_reserve_restart_multiples"] *
        contract["source_restart_size_bytes"],
    )
    for row in filesystems:
        _require(isinstance(row, dict) and set(row) == {
            "device", "roles", "observations", "role_contribution_bytes",
            "contribution_bytes_total", "reserve_bytes", "required_free_bytes",
            "status",
        }, f"disk preflight {phase} filesystem evidence shape is invalid")
        device = _integer(row.get("device"), f"disk preflight {phase} device")
        roles = row.get("roles")
        _require(device >= 0 and device not in seen_devices and
                 isinstance(roles, list) and roles == sorted(roles) and roles and
                 set(roles) <= {"state_dir", "staging_dir"} and
                 not (set(roles) & seen_roles),
                 f"disk preflight {phase} filesystem roles/devices are invalid")
        seen_devices.add(device)
        seen_roles.update(roles)
        _require(all(expected_devices.get(role) == device for role in roles),
                 f"disk preflight {phase} device binding is invalid")
        contributions = {role: expected_contributions[role] for role in roles}
        contribution_total = sum(contributions.values())
        required_free = max(
            contract["additional_hard_minimum_free_bytes"],
            reserve + contribution_total,
        )
        _require(row.get("role_contribution_bytes") == contributions and
                 row.get("contribution_bytes_total") == contribution_total and
                 row.get("reserve_bytes") == reserve and
                 row.get("required_free_bytes") == required_free,
                 f"disk preflight {phase} budget arithmetic is invalid")
        observations = row.get("observations")
        _require(isinstance(observations, dict) and set(observations) == set(roles),
                 f"disk preflight {phase} does not measure every role")
        for role in roles:
            observation = observations[role]
            access = observation.get("access") if isinstance(observation, dict) \
                else None
            expected_path = state_dir if role == "state_dir" else staging_dir
            access_valid = (access == expected_access[role]
                            if phase == "before_spawn" else
                            isinstance(access, dict) and
                            set(access) == {"method", "fd", "planned_path"} and
                            access.get("method") ==
                            "bound_directory_fstatvfs_v1" and
                            isinstance(access.get("fd"), int) and
                            not isinstance(access.get("fd"), bool) and
                            access["fd"] >= 0 and
                            access.get("planned_path") == expected_path)
            _require(isinstance(observation, dict) and set(observation) == {
                "access", "statvfs", "used_percent_pass", "free_bytes_pass",
            } and access_valid,
                     f"disk preflight {phase} {role} access proof is invalid")
            fs = observation.get("statvfs")
            _require(isinstance(fs, dict) and set(fs) == {
                "fragment_size_bytes", "total_blocks", "free_blocks",
                "available_blocks", "total_bytes", "free_bytes",
                "available_bytes", "effective_free_bytes", "used_bytes",
                "used_percent_numerator", "used_percent_denominator",
            }, f"disk preflight {phase} statvfs evidence shape is invalid")
            fragment = _integer(
                fs.get("fragment_size_bytes"), "statvfs fragment size")
            total_blocks = _integer(fs.get("total_blocks"), "statvfs total blocks")
            free_blocks = _integer(fs.get("free_blocks"), "statvfs free blocks")
            available_blocks = _integer(
                fs.get("available_blocks"), "statvfs available blocks")
            _require(fragment > 0 and total_blocks > 0 and
                     0 <= available_blocks <= free_blocks <= total_blocks,
                     f"disk preflight {phase} statvfs block counts are invalid")
            free_bytes = free_blocks * fragment
            available_bytes = available_blocks * fragment
            used_bytes = (total_blocks - free_blocks) * fragment
            used_numerator = (total_blocks - free_blocks) * 100
            used_denominator = total_blocks - free_blocks + available_blocks
            _require(used_denominator > 0,
                     f"disk preflight {phase} has no usable filesystem blocks")
            _require(fs == {
                "fragment_size_bytes": fragment,
                "total_blocks": total_blocks,
                "free_blocks": free_blocks,
                "available_blocks": available_blocks,
                "total_bytes": total_blocks * fragment,
                "free_bytes": free_bytes,
                "available_bytes": available_bytes,
                "effective_free_bytes": min(free_bytes, available_bytes),
                "used_bytes": used_bytes,
                "used_percent_numerator": used_numerator,
                "used_percent_denominator": used_denominator,
            }, f"disk preflight {phase} statvfs byte arithmetic is invalid")
            used_pass = used_numerator < (
                used_denominator * contract["used_percent_exclusive_max"])
            free_pass = min(free_bytes, available_bytes) >= required_free
            _require(used_pass and free_pass and
                     observation.get("used_percent_pass") is True and
                     observation.get("free_bytes_pass") is True,
                     f"disk preflight {phase} does not prove {role} capacity")
        _require(row.get("status") == "pass",
                 f"disk preflight {phase} filesystem status is invalid")
    _require(seen_roles == {"state_dir", "staging_dir"},
             f"disk preflight {phase} does not cover both directory roles")
    if phase == "before_staging":
        measured_fds = {
            observation["access"]["fd"]
            for row in filesystems
            for observation in row["observations"].values()
        }
        _require(len(measured_fds) == 2,
                 "before-staging disk preflight must bind two distinct directories")
    return value


def _load_json(path: Path) -> tuple[dict[str, Any], bytes]:
    raw = _read_stable_bytes(path)
    try:
        value = strict_json_loads(raw, str(path))
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise CheckFailure(f"{path}: invalid JSON: {exc}") from exc
    _require(isinstance(value, dict), f"{path}: JSON root must be an object")
    return value, raw


def _signature(path: Path) -> tuple[int, int, int, int, int]:
    try:
        value = path.lstat()
    except OSError as exc:
        raise CheckFailure(f"cannot stat {path}: {exc}") from exc
    _require(stat.S_ISREG(value.st_mode), f"{path}: must be a regular, non-symlink file")
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns,
            value.st_ctime_ns)


def _reject_output_ancestors(path: Path) -> None:
    """Reject a symlink at the report path or any already-existing parent."""

    current = Path(path).absolute()
    while True:
        if current.is_symlink():
            raise CheckFailure(f"report path traverses a symlink: {current}")
        if current.parent == current:
            break
        current = current.parent


def _read_stable_bytes(path: Path) -> bytes:
    checked_path, stream, expected = _open_regular_nofollow(path)
    try:
        exempt = {(os.getpid(), stream.fileno())}
        _assert_closed(checked_path, expected, exempt)
        pieces: list[bytes] = []
        for chunk in iter(lambda: stream.read(HASH_CHUNK_BYTES), b""):
            pieces.append(chunk)
        _assert_stream_signature(stream, checked_path, expected, "while reading")
        _assert_path_signature(checked_path, expected, "while reading")
        _assert_closed(checked_path, expected, exempt)
    finally:
        stream.close()
    _assert_path_signature(checked_path, expected, "after reading")
    _assert_closed(checked_path, expected)
    raw = b"".join(pieces)
    return raw


def _hash_record(path: Path) -> dict[str, Any]:
    audited = stable_sha256(path)
    return {
        "path": audited["path"],
        **{name: audited[name] for name in
           ("device", "inode", "size", "mtime_ns", "ctime_ns")},
        "sha256": audited["sha256"],
        "closure_check": audited["closure_check"],
    }


def _validate_planned_file_record(value: Any, name: str) -> dict[str, Any]:
    _require(isinstance(value, dict), f"{name} must be an object")
    path = value.get("path")
    size = value.get("size")
    digest = value.get("sha256")
    _require(isinstance(path, str) and bool(path) and Path(path).is_absolute(),
             f"{name}.path must be an absolute nonempty path")
    _require(isinstance(size, int) and not isinstance(size, bool) and size >= 0,
             f"{name}.size must be a nonnegative integer")
    _require(isinstance(digest, str) and
             re.fullmatch(r"[0-9a-f]{64}", digest) is not None,
             f"{name}.sha256 must be lowercase SHA-256")
    return value


def _verify_planned_file(value: dict[str, Any], name: str) -> dict[str, Any]:
    planned_path = Path(value["path"])
    canonical = planned_path.resolve(strict=True)
    _require(value["path"] == str(canonical),
             f"{name} path is not canonical absolute")
    record = _hash_record(canonical)
    _require(record["size"] == value["size"] and
             record["sha256"] == value["sha256"],
             f"{name} differs from the immutable plan")
    return record


def _verify_current_tool(planned: dict[str, Any], path: Path,
                         name: str) -> dict[str, Any]:
    canonical = path.resolve(strict=True)
    _require(planned.get("path") == str(canonical),
             f"current {name} path differs from the plan-bound tool")
    return _verify_planned_file(planned, f"current {name}")


def _git_output(repo: Path, git_path: str, *arguments: str) -> str:
    try:
        result = subprocess.run(
            [git_path, "-C", str(repo), *GIT_CONFIGURATION, *arguments], check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            env=GIT_ENVIRONMENT,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        detail = getattr(exc, "stderr", "") or str(exc)
        raise CheckFailure(f"cannot verify planned repository: {detail.strip()}") from exc
    return result.stdout.strip()


def _verify_planned_repository(value: dict[str, Any],
                               git_record: dict[str, Any]) -> dict[str, Any]:
    _verify_planned_file(git_record, "planned git executable")
    git_path = git_record["path"]
    repo = Path(value["path"]).resolve(strict=True)
    _require(repo.is_dir(), "planned repository is not a directory")
    top = Path(_git_output(
        repo, git_path, "rev-parse", "--show-toplevel")).resolve(strict=True)
    _require(top == repo, "planned repository path is not its worktree root")
    commit = _git_output(repo, git_path, "rev-parse", "--verify", "HEAD")
    _require(commit == value["commit"], "planned repository commit changed")
    status = _git_output(
        repo, git_path, "status", "--porcelain=v1", "--untracked-files=all")
    _require(not status, "planned repository is no longer clean")
    return {
        "path": str(repo), "commit": commit, "clean": True,
        "git_environment_policy": "explicit_clean_environment_v1",
        "git_environment": dict(GIT_ENVIRONMENT),
        "git_configuration": list(GIT_CONFIGURATION),
    }


def _assert_binding_unchanged(path: Path, binding: dict[str, Any]) -> None:
    current = _signature(path)
    expected = tuple(binding[name] for name in
                     ("device", "inode", "size", "mtime_ns", "ctime_ns"))
    _require(current == expected, f"bound evidence changed after inspection: {path}")


def _check_ready(paths: Iterable[Path], now: float | None = None) -> None:
    current = time.time() if now is None else now
    recent: list[dict[str, Any]] = []
    for path in dict.fromkeys(paths):
        signature = _signature(path)
        age = current - signature[3] / 1.0e9
        if age < MIN_READY_AGE_SECONDS:
            recent.append({"path": str(path), "age_seconds": max(age, 0.0)})
    if recent:
        raise NotReady(json.dumps({
            "minimum_age_seconds": MIN_READY_AGE_SECONDS,
            "recent": recent,
        }, sort_keys=True))


def _close_in_ulps(left: float, right: float, ulps: int) -> bool:
    if not (math.isfinite(left) and math.isfinite(right)):
        return False
    if left == right:
        return True
    tolerance = max(math.ulp(left), math.ulp(right)) * ulps
    return abs(left - right) <= tolerance


def _sequential_add(initial: float, increment: float, count: int) -> float:
    """Reproduce C++ binary64 ``time = time + dt`` without reassociation."""

    result = float(initial)
    for _ in range(count):
        result += increment
    return result


def _canonical_json_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _canonical_launch_environment(value: Any) -> dict[str, str]:
    _require(isinstance(value, dict) and
             set(value) == {"kind", "values", "sha256", "rank_projection"} and
             value.get("kind") == LAUNCH_ENVIRONMENT_KIND,
             f"launch environment kind must be {LAUNCH_ENVIRONMENT_KIND}")
    values = value.get("values")
    expected_keys = {
        "HOME", "LANG", "LC_ALL", "CUDA_DEVICE_ORDER",
        "PRTE_MCA_schizo_proxy",
    }
    _require(isinstance(values, dict) and set(values) == expected_keys and
             all(isinstance(item, str) and "\x00" not in item
                 for item in values.values()),
             "launch environment values are not the exact safe string set")
    try:
        passwd_home = pwd.getpwuid(os.geteuid()).pw_dir
    except KeyError as exc:
        raise CheckFailure("effective uid has no passwd database entry") from exc
    _require(isinstance(passwd_home, str) and passwd_home and
             Path(passwd_home).is_absolute(),
             "effective uid passwd HOME is invalid")
    home_absolute = Path(values["HOME"]).absolute()
    try:
        home = home_absolute.resolve(strict=True)
        home_info = home_absolute.lstat()
    except OSError as exc:
        raise CheckFailure(f"cannot audit launch HOME: {exc}") from exc
    _require(values == {
        "HOME": passwd_home,
        "LANG": "C",
        "LC_ALL": "C",
        "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
        "PRTE_MCA_schizo_proxy": "ompi",
    } and str(home) == passwd_home and home_absolute == home and
             stat.S_ISDIR(home_info.st_mode) and
             not stat.S_ISLNK(home_info.st_mode) and
             home_info.st_uid == os.geteuid() and
             not (stat.S_IMODE(home_info.st_mode) & 0o022),
             "launch HOME/locale/CUDA values are not canonical and safe")
    canonical = {
        "HOME": values["HOME"], "LANG": "C", "LC_ALL": "C",
        "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
        "PRTE_MCA_schizo_proxy": "ompi",
    }
    _require(value.get("sha256") == _canonical_json_sha256(canonical),
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
             projection.get("sha256") ==
             _canonical_json_sha256(projection_payload),
             "rank environment projection differs from the exact PRRTE contract")
    return canonical


def _snapshot_mca_configuration_file(scope: str, project: str,
                                     path: Path) -> dict[str, Any]:
    """Checker-owned re-observation of one default MCA parameter path."""

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
        raise CheckFailure(f"cannot inspect MCA configuration {absolute}: {exc}") \
            from exc
    try:
        resolved = absolute.resolve(strict=True)
    except OSError as exc:
        raise CheckFailure(f"cannot resolve MCA configuration {absolute}: {exc}") \
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
        raise CheckFailure(f"cannot audit MCA configuration {absolute}: {exc}") \
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
        raise CheckFailure(f"cannot bind MCA prefix {absolute}: {exc}") from exc
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
    _validate_mca_configuration_contract(value, home)
    _require(isinstance(value, dict) and
             set(value) == {"kind", "home", "prefix", "prefix_directory",
                            "files", "sha256"} and
             value.get("kind") == MCA_CONFIGURATION_KIND,
             f"MCA configuration contract must be {MCA_CONFIGURATION_KIND}")
    prefix_directory = _snapshot_mca_prefix_directory(Path(value["prefix"]))
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
        current = _snapshot_mca_configuration_file(scope, project, expected_path)
        _require(planned_files[index] == current,
                 f"MCA configuration {expected_path} differs from immutable plan")
        current_files.append(current)
    payload = {
        "kind": MCA_CONFIGURATION_KIND,
        "home": home, "prefix": str(prefix),
        "prefix_directory": prefix_directory, "files": current_files,
    }
    _require(value.get("sha256") == _canonical_json_sha256(payload),
             "MCA configuration contract SHA-256 is not canonical")
    return {**payload, "sha256": value["sha256"]}


def _validate_mca_configuration_contract(value: Any,
                                         home: str) -> dict[str, Any]:
    """Validate the immutable six-path contract without trusting live files."""

    _require(isinstance(value, dict) and
             set(value) == {"kind", "home", "prefix", "prefix_directory",
                            "files", "sha256"} and
             value.get("kind") == MCA_CONFIGURATION_KIND,
             f"MCA configuration contract must be {MCA_CONFIGURATION_KIND}")
    prefix_value = value.get("prefix")
    prefix_directory = value.get("prefix_directory")
    _require(isinstance(prefix_value, str) and Path(prefix_value).is_absolute() and
             isinstance(prefix_directory, dict) and
             set(prefix_directory) == {
                 "path", "device", "inode", "owner_uid", "mode",
             } and prefix_directory.get("path") == prefix_value and
             all(isinstance(prefix_directory.get(name), int) and
                 prefix_directory[name] >= 0
                 for name in ("device", "inode", "owner_uid")) and
             prefix_directory["owner_uid"] in (0, os.geteuid()) and
             isinstance(prefix_directory.get("mode"), str) and
             re.fullmatch(r"0[0-7]{3}", prefix_directory["mode"]) is not None and
             not (int(prefix_directory["mode"], 8) & 0o022) and
             value.get("home") == home,
             "MCA configuration HOME/prefix directory contract is invalid")
    prefix = Path(prefix_value)
    files = value.get("files")
    _require(isinstance(files, list) and len(files) == len(MCA_CONFIGURATION_LAYOUT),
             "MCA configuration file set is not exact")
    for index, (scope, project, relative) in enumerate(MCA_CONFIGURATION_LAYOUT):
        record = files[index]
        expected_path = (Path(home) if scope == "home" else prefix) / relative
        _require(isinstance(record, dict) and
                 record.get("scope") == scope and
                 record.get("project") == project and
                 record.get("path") == str(expected_path) and
                 record.get("state") in ("absent", "present"),
                 f"MCA configuration record {index} is not canonical")
        if record.get("state") == "absent":
            _require(set(record) == {"scope", "project", "path", "state"},
                     f"absent MCA configuration record {index} is not exact")
        else:
            _require(set(record) == {
                "scope", "project", "path", "state", "device", "inode",
                "owner_uid", "mode", "size", "mtime_ns", "ctime_ns", "sha256",
                "closure_check",
            } and all(isinstance(record[name], int) and record[name] >= 0
                      for name in ("device", "inode", "owner_uid", "size",
                                   "mtime_ns", "ctime_ns")) and
                     record["owner_uid"] in (0, os.geteuid()) and
                     isinstance(record["mode"], str) and
                     re.fullmatch(r"0[0-7]{3}", record["mode"]) is not None and
                     not (int(record["mode"], 8) & 0o022) and
                     isinstance(record["sha256"], str) and
                     re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) is not None and
                     record["closure_check"] == "linux_proc_fd",
                     f"present MCA configuration record {index} is not exact/safe")
    payload = {
        "kind": MCA_CONFIGURATION_KIND, "home": home,
        "prefix": str(prefix), "prefix_directory": prefix_directory,
        "files": files,
    }
    _require(value.get("sha256") == _canonical_json_sha256(payload),
             "MCA configuration contract SHA-256 is not canonical")
    return value


def _wall_time_token(seconds: int) -> str:
    _require(seconds > 0, "launch wall time must be positive")
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds_part = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds_part:02d}"


def _validate_rank_environment(value: Any, *, global_rank: int,
                               local_rank: int, world_size: int,
                               launcher_pid: int, hostname: str,
                               state_dir: str,
                               athena_argv: Sequence[str]) -> dict[str, str]:
    """Independently re-derive the complete observed 5.0.9 rank environment."""

    _require(isinstance(value, dict) and
             all(isinstance(key, str) and isinstance(item, str)
                 for key, item in value.items()),
             f"launch rank {global_rank} environment is not a string mapping")
    actual_keys = set(value)
    expected_keys = set(RANK_ENVIRONMENT_KEYS)
    _require(actual_keys == expected_keys,
             f"launch rank {global_rank} environment key closure differs: "
             f"missing={sorted(expected_keys - actual_keys)!r}, "
             f"unexpected={sorted(actual_keys - expected_keys)!r}")
    namespace_family = f"prterun-{hostname}-{launcher_pid}"
    namespace = f"{namespace_family}@1"
    inherited = {
        "HOME": str(Path(pwd.getpwuid(os.geteuid()).pw_dir).resolve(strict=True)),
        "LANG": "C", "LC_ALL": "C", "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
    }
    expected = {
        **inherited, **RANK_FIXED_ENVIRONMENT_VALUES,
        "OMPI_COMMAND": str(athena_argv[0]),
        "OMPI_ARGV": " ".join(athena_argv[1:]),
        "OMPI_COMM_WORLD_RANK": str(global_rank),
        "OMPI_COMM_WORLD_SIZE": str(world_size),
        "OMPI_COMM_WORLD_LOCAL_RANK": str(local_rank),
        "OMPI_COMM_WORLD_LOCAL_SIZE": str(world_size),
        "OMPI_COMM_WORLD_NODE_RANK": str(local_rank),
        "OMPI_FILE_LOCATION": f"/tmp/ompi.{launcher_pid}/1/{global_rank}",
        "OMPI_MCA_initial_wdir": state_dir,
        "OMPI_MCA_num_procs": str(world_size),
        "OMPI_WORLD_LOCAL_SIZE": str(world_size),
        "OMPI_WORLD_SIZE": str(world_size), "PMIX_HOSTNAME": hostname,
        "PMIX_NAMESPACE": namespace, "PMIX_RANK": str(global_rank),
        "PMIX_SERVER_TMPDIR": f"/tmp/ompi.{launcher_pid}", "PWD": state_dir,
    }
    for key, expected_value in expected.items():
        _require(value.get(key) == expected_value,
                 f"launch rank {global_rank} environment {key} differs from "
                 "derived contract")
    uri_values = {value[key] for key in RANK_URI_ENVIRONMENT_KEYS}
    _require(len(uri_values) == 1,
             f"launch rank {global_rank} PMIx server URI aliases differ")
    uri = next(iter(uri_values))
    match = re.fullmatch(
        rf"{re.escape(namespace_family)}@0\.0;tcp4://127\.0\.0\.1:(\d+)",
        uri)
    _require(match is not None and 1 <= int(match.group(1)) <= 65535,
             f"launch rank {global_rank} PMIx server URI is invalid")
    return {key: value[key] for key in RANK_ENVIRONMENT_KEYS}


def validate_plan(plan: dict[str, Any]) -> dict[str, Any]:
    _require(plan.get("schema") == SCHEMA, f"plan schema must equal {SCHEMA}")
    expected = plan.get("expected")
    policy = plan.get("policy")
    inputs = plan.get("inputs")
    outputs = plan.get("outputs")
    overrides = plan.get("command_overrides")
    _require(isinstance(expected, dict), "plan expected must be an object")
    _require(isinstance(policy, dict), "plan policy must be an object")
    _require(isinstance(inputs, dict), "plan inputs must be an object")
    _require(isinstance(outputs, list), "plan outputs must be an array")
    _require(isinstance(overrides, dict), "plan command_overrides must be an object")

    source_cycle = _integer(expected.get("source_cycle"), "expected.source_cycle")
    final_cycle = _integer(expected.get("final_cycle"), "expected.final_cycle")
    root_steps = _integer(expected.get("root_steps"), "expected.root_steps")
    guard_steps = _integer(
        expected.get("tlim_guard_steps"), "expected.tlim_guard_steps")
    source_time = _number(expected.get("source_time"), "expected.source_time")
    final_time = _number(expected.get("final_time"), "expected.final_time")
    root_dt = _number(expected.get("root_dt"), "expected.root_dt")
    tlim = _number(expected.get("tlim"), "expected.tlim")
    _require(root_steps > 0 and final_cycle == source_cycle + root_steps,
             "expected final_cycle must equal source_cycle + positive root_steps")
    _require(guard_steps > 0, "expected tlim_guard_steps must be positive")
    _require(root_dt > 0.0 and final_time > source_time,
             "expected root_dt and segment time span must be positive")
    recomputed_final = _sequential_add(source_time, root_dt, root_steps)
    recomputed_tlim = _sequential_add(recomputed_final, root_dt, guard_steps)
    _require(final_time == recomputed_final,
             "expected final_time is not the sequential binary64 root-step endpoint")
    _require(tlim == recomputed_tlim,
             "expected tlim is not the sequential binary64 guard endpoint")
    _require(tlim > final_time, "expected tlim must be a nonbinding upper guard")
    capacity_transition = _validate_capacity_transition(plan)
    restart_cadence_transition = _validate_restart_cadence_transition(
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
    _require(overrides == expected_overrides and
             overrides.get("time/nlim") == final_cycle and
             _number(overrides.get("time/tlim"),
                     "command_overrides.time/tlim") == tlim and
             _number(overrides.get("output3/dt"),
                     "command_overrides.output3/dt") == root_dt,
             "command_overrides must be exactly the canonical limits, divB "
             "cadence, optional restart tightening, and optional capacity increase")

    for name in ("binary", "trajectory", "source_restart"):
        _validate_planned_file_record(inputs.get(name), f"inputs.{name}")
    source_input = inputs["source_restart"]
    repo_input = inputs.get("repo")
    _require(isinstance(repo_input, dict) and
             isinstance(repo_input.get("path"), str) and
             Path(repo_input["path"]).is_absolute() and
             re.fullmatch(r"[0-9a-f]{40,64}", str(repo_input.get("commit", "")))
             is not None and repo_input.get("clean") is True,
             "inputs.repo must bind an absolute clean repository and commit")
    tools = plan.get("tools")
    _require(isinstance(tools, dict) and tools.get("hash_algorithm") == "sha256",
             "plan tools must specify SHA-256")
    for name in PLANNED_TOOL_NAMES:
        _validate_planned_file_record(tools.get(name), f"tools.{name}")
    source_qualification = plan.get("source_qualification")
    _require(isinstance(source_qualification, dict) and
             source_qualification.get("mode") in
             ("anchor_full_audit", "parent_segment_pass"),
             "plan must contain exactly one supported source qualification mode")
    source_audit = source_qualification.get("audit")
    _require(isinstance(source_audit, dict) and
             source_audit.get("valid") is True and
             source_audit.get("sha256") == source_input["sha256"] and
             source_audit.get("stored_reals", {}).get("nonfinite_count") == 0 and
             source_audit.get("stored_reals", {}).get("finite_count") ==
             source_audit.get("stored_reals", {}).get("count") and
             source_audit.get("layout", {}).get("expected_file_size") ==
             source_input["size"] and
             source_audit.get("topology", {}).get(
                 "complete_leaf_coverage") is True,
             "plan source full audit is invalid or incomplete")
    if source_qualification["mode"] == "parent_segment_pass":
        parent = source_qualification.get("parent_segment_pass")
        _validate_planned_file_record(parent, "source_qualification.parent_segment_pass")
        _require(source_qualification.get("parent_qualification_mode") in
                 ("complete_segment_v1", "legacy_complete_segment_v1",
                  "scheduled_prefix_recovery_v1"),
                 "parent source qualification lacks an explicit supported mode")
    else:
        _require("parent_segment_pass" not in source_qualification,
                 "anchor source qualification may not contain a parent pass")
        _require("parent_qualification_mode" not in source_qualification,
                 "anchor source qualification may not contain a parent mode")
    source_baryon = source_qualification.get("source_baryon_mass")
    _require(isinstance(source_baryon, dict) and
             _number(source_baryon.get("value"), "source_baryon_mass.value") > 0.0 and
             _number(source_baryon.get("time"), "source_baryon_mass.time") == source_time,
             "source qualification must bind positive endpoint baryon mass")
    _validate_planned_file_record(
        source_baryon.get("evidence"), "source_baryon_mass.evidence")

    termination_policy = policy.get("segment_termination")
    _require(isinstance(termination_policy, dict),
             "policy.segment_termination must be an object")
    ranks = _integer(policy.get("ranks"), "policy.ranks")
    capacity_policy = policy.get("capacity")
    _require(isinstance(capacity_policy, dict) and
             capacity_policy.get("ranks") == ranks,
             "policy.capacity.ranks must equal policy.ranks")
    memory_limit = _number(policy.get("gpu_exit_memory_mib_max"),
                           "policy.gpu_exit_memory_mib_max")
    ulps = _integer(policy.get("endpoint_time_ulp_tolerance"),
                    "policy.endpoint_time_ulp_tolerance")
    _require(termination_policy.get("endpoint_time_max_ulps") == ulps,
             "top-level and termination endpoint ULP policies disagree")
    _require(ulps == 0,
             "endpoint time ULP tolerance must be exactly zero")
    _require(termination_policy == {
        "primary": "cycle_limit",
        "endpoint_time_max_ulps": 0,
        "tlim_role": "nonbinding_guard",
        "require_exact_final_cycle": True,
        "require_exact_binary64_endpoint": True,
    }, "segment termination policy is not the exact strict campaign policy")
    fixed_dt = policy.get("fixed_root_dt_mode")
    _require(fixed_dt == {
        "enabled": True,
        "root_dt": root_dt,
        "root_dt_parameter": "time/root_dt_max",
        "require_source_last_dt_exact": True,
        "sequential_binary64_addition": True,
    }, "fixed-root-dt policy is not the exact campaign policy")
    _require(ranks > 0, "policy.ranks must be positive")
    _require(memory_limit >= 0.0, "GPU exit-memory limit must be nonnegative")
    _require(0 <= ulps <= 64, "endpoint ULP tolerance must be in [0,64]")
    launch_contract = plan.get("launch_contract")
    _require(isinstance(launch_contract, dict),
             "plan launch_contract must be an object")
    launcher = _validate_planned_file_record(
        launch_contract.get("launcher"), "launch_contract.launcher")
    state_dir = launch_contract.get("state_dir")
    wall_seconds = launch_contract.get("wall_time_seconds")
    _require(isinstance(state_dir, str) and Path(state_dir).is_absolute(),
             "launch_contract.state_dir must be absolute")
    evidence_dir = launch_contract.get("evidence_dir")
    plan_path = launch_contract.get("plan_path")
    _require(isinstance(evidence_dir, str) and Path(evidence_dir).is_absolute() and
             isinstance(plan_path, str) and
             plan_path == str(Path(evidence_dir) / "segment.plan.json"),
             "launch_contract.plan_path must be the canonical direct evidence plan")
    launch_environment = _canonical_launch_environment(
        launch_contract.get("environment"))
    _validate_mca_configuration_contract(
        launch_contract.get("mca_configuration"), launch_environment["HOME"])
    directory_transport = launch_contract.get("directory_transport")
    _require(isinstance(directory_transport, dict) and
             set(directory_transport) == {"kind", "holder_pid_token", "roles"} and
             directory_transport.get("kind") == DIRECTORY_TRANSPORT_KIND and
             directory_transport.get("holder_pid_token") == HOLDER_PID_TOKEN,
             "launch contract lacks canonical fixed directory descriptors")
    directory_roles = directory_transport.get("roles")
    _require(isinstance(directory_roles, dict) and
             set(directory_roles) == {"state_dir", "evidence_dir"},
             "directory transport roles must be exactly state_dir/evidence_dir")
    for role, descriptor, planned_path in (
        ("state_dir", STATE_DIRECTORY_FD, state_dir),
        ("evidence_dir", EVIDENCE_DIRECTORY_FD, evidence_dir),
    ):
        _require(directory_roles.get(role) == {
            "role": role,
            "fd": descriptor,
            "planned_path": planned_path,
            "proc_path_template": f"/proc/{HOLDER_PID_TOKEN}/fd/{descriptor}",
        }, f"directory transport role {role} is not canonical")
    executable_transport = launch_contract.get("executable_transport")
    _require(isinstance(executable_transport, dict) and
             set(executable_transport) == {"kind", "holder_pid_token", "roles"} and
             executable_transport.get("kind") == EXECUTABLE_TRANSPORT_KIND and
             executable_transport.get("holder_pid_token") == HOLDER_PID_TOKEN,
             "launch contract lacks canonical executable descriptors")
    executable_roles = executable_transport.get("roles")
    _require(isinstance(executable_roles, dict) and
             set(executable_roles) == {"launcher", "binary"},
             "executable transport roles must be exactly launcher/binary")
    for role, descriptor, parent_path in (
        ("launcher", LAUNCHER_EXECUTABLE_FD,
         launch_contract.get("launcher", {}).get("path")),
        ("binary", BINARY_EXECUTABLE_FD, inputs["binary"]["path"]),
    ):
        _require(executable_roles.get(role) == {
            "role": role, "fd": descriptor, "parent_path": parent_path,
            "proc_path_template": f"/proc/{HOLDER_PID_TOKEN}/fd/{descriptor}",
        }, f"executable transport role {role} is not canonical")
    transport = launch_contract.get("input_transport")
    _require(isinstance(transport, dict) and
             transport.get("kind") == INPUT_TRANSPORT_KIND and
             transport.get("holder_pid_token") == HOLDER_PID_TOKEN,
             "launch contract lacks the canonical fixed-descriptor input transport")
    staging_dir = transport.get("staging_dir")
    _require(isinstance(staging_dir, str) and Path(staging_dir).is_absolute(),
             "input_transport.staging_dir must be absolute")
    _validate_disk_preflight_contract(
        launch_contract.get("disk_preflight"), source_input["size"],
        inputs["trajectory"]["size"])
    roles = transport.get("roles")
    _require(isinstance(roles, dict) and
             set(roles) == {"source_restart", "trajectory"},
             "input transport roles must be exactly source_restart/trajectory")
    source_template = f"/proc/{HOLDER_PID_TOKEN}/fd/{SOURCE_RESTART_FD}"
    trajectory_template = f"/proc/{HOLDER_PID_TOKEN}/fd/{TRAJECTORY_FD}"
    expected_roles = {
        "source_restart": (SOURCE_RESTART_FD, source_template,
                           inputs["source_restart"]),
        "trajectory": (TRAJECTORY_FD, trajectory_template,
                       inputs["trajectory"]),
    }
    staged: dict[str, dict[str, Any]] = {}
    for role, (descriptor, proc_template, parent_input) in expected_roles.items():
        role_record = roles.get(role)
        _require(isinstance(role_record, dict) and
                 role_record.get("role") == role and
                 role_record.get("fd") == descriptor and
                 role_record.get("proc_path_template") == proc_template,
                 f"input transport role {role} is not canonical")
        staged[role] = _validate_planned_file_record(
            role_record.get("staged_file"),
            f"input_transport.roles.{role}.staged_file")
        _require(Path(staged[role]["path"]).parent == Path(staging_dir) and
                 staged[role]["size"] == parent_input["size"] and
                 staged[role]["sha256"] == parent_input["sha256"],
                 f"staged {role} target does not bind its parent input bytes")
    serialized_trajectory = plan.get("source", {}).get(
        "parameters", {}).get("problem", {}).get("trajectory_file")
    _require(isinstance(serialized_trajectory, str) and serialized_trajectory,
             "plan source lacks serialized trajectory provenance")
    _require(transport.get("parent_content_identity") == {
        "source_restart_sha256": inputs["source_restart"]["sha256"],
        "trajectory_sha256": inputs["trajectory"]["sha256"],
        "source_serialized_trajectory_path": serialized_trajectory,
    }, "input transport parent content identity is not exact")
    _require(transport.get("trajectory_rebinding") == {
        "parameter": "problem/trajectory_file",
        "parent_sha256": inputs["trajectory"]["sha256"],
        "runtime_value_template": trajectory_template,
    }, "input transport trajectory rebinding is not exact")
    _require(isinstance(wall_seconds, int) and not isinstance(wall_seconds, bool) and
             wall_seconds > 0,
             "launch_contract.wall_time_seconds must be a positive integer")
    expected_mpi = [launcher["path"], "--allow-run-as-root", "--bind-to", "none",
                    "-np", str(ranks)]
    expected_athena = [
        f"/proc/{HOLDER_PID_TOKEN}/fd/{BINARY_EXECUTABLE_FD}",
        "--kokkos-map-device-id-by=mpi_rank",
        "-r", source_template,
        "-d", f"/proc/{HOLDER_PID_TOKEN}/fd/{STATE_DIRECTORY_FD}",
        "-t", _wall_time_token(wall_seconds),
        f"time/nlim={final_cycle}", f"time/tlim={tlim!r}",
        f"problem/trajectory_file={trajectory_template}",
        f"output3/dt={root_dt!r}",
    ]
    if restart_cadence_transition["kind"] == "tighten_v1":
        expected_athena.append(restart_cadence_transition["runtime_override"])
    if capacity_transition["kind"] == "increase_v1":
        expected_athena.append(capacity_transition["runtime_override"])
    _require(launch_contract.get("world_size") == ranks and
             launch_contract.get("gpu_count") == ranks and
             launch_contract.get("mpi_argv") == expected_mpi and
             launch_contract.get("athena_argv_template") == expected_athena and
             "athena_argv" not in launch_contract,
             "plan launch_contract is not the canonical single-node MPI launch")
    evidence = launch_contract.get("evidence")
    _require(isinstance(evidence_dir, str) and Path(evidence_dir).is_absolute() and
             isinstance(evidence, dict),
             "launch contract must bind a dedicated evidence directory")
    expected_evidence = {
        "launch_record": str(Path(evidence_dir) / "segment.launch.ready"),
        "completion_record": str(Path(evidence_dir) / "segment.completion.ready"),
        "run_log": str(Path(evidence_dir) / "run.log"),
        "exit_status": str(Path(evidence_dir) / "exit.status"),
        "gpu_before": str(Path(evidence_dir) / "gpu-before.csv"),
        "gpu_after": str(Path(evidence_dir) / "gpu-after.csv"),
    }
    _require(evidence == expected_evidence,
             "launch contract evidence paths are not canonical")
    _require(launch_environment == launch_contract["environment"]["values"],
             "launch environment canonicalization is internally inconsistent")
    mutable = policy.get("mutable_parameters")
    _require(mutable == CANONICAL_MUTABLE_PARAMETERS,
             "policy.mutable_parameters differs from the strict whitelist")
    _require(capacity_policy == {
        "ranks": ranks,
        "minimum_headroom_blocks_hard": 64,
        "minimum_headroom_blocks_yellow": 128,
    }, "policy.capacity differs from the strict campaign thresholds")
    normalized_events = _normalize_event_thresholds(policy.get("event_thresholds"))
    _require(normalized_events == [
        {"name": "fofc_per_test", "numerator": "fofc",
         "denominator": "fofc_tests", "max_ratio": 0.005},
        {"name": "cons_adjust_per_c2p_call", "numerator": "cons_adjust",
         "denominator": "c2p_calls",
         "max_ratio": 0.005},
        {"name": "mag_adjust_per_c2p_call", "numerator": "mag_adjust",
         "denominator": "c2p_calls", "max_ratio": 0.005},
    ], "policy event ratios differ from strict campaign thresholds")
    absolute_events = policy.get("event_absolute_thresholds")
    _require(absolute_events == {
        "hard_equal_zero": ["eos_fail", "eos_vceil"],
        "c2p_iterations_exclusive_max": 25,
    }, "absolute event policy differs from the exact strict thresholds")
    _require(policy.get("nonfinite_count_max") == 0,
             "policy.nonfinite_count_max must be exactly zero")
    _require(policy.get("scheduled_prefix_recovery") == PREFIX_RECOVERY_POLICY,
             "scheduled-prefix recovery policy differs from the exact contract")
    divb_policy = policy.get("divb_max_abs")
    mass_policy = policy.get("baryon_mass_fractional_loss")
    _require(isinstance(divb_policy, dict), "policy.divb_max_abs must be an object")
    _require(isinstance(mass_policy, dict),
             "policy.baryon_mass_fractional_loss must be an object")
    _require(divb_policy == {"hard": 1.0e-8, "yellow": 1.0e-11},
             "divB policy differs from strict campaign thresholds")
    _require(mass_policy == {
        "hard_per_root_step": 0.005,
        "yellow_per_root_step": 0.0025,
        "yellow_per_48M": 0.02,
        "rolling_window_root_steps": 10,
    }, "baryon policy differs from strict campaign thresholds")
    _require(memory_limit == 100.0,
             "GPU exit-memory limit must equal the strict value 100 MiB")
    _require(policy.get("gpu_ecc") == {
        "corrected_before_max": 0, "corrected_after_max": 0,
        "uncorrected_before_max": 0, "uncorrected_after_max": 0,
    }, "GPU ECC policy differs from strict campaign thresholds")
    _require(policy.get("restart_contract") == {
        "real_bytes": 8,
        "subcycling": "level",
        "allow_legacy_subcycling_restart": False,
        "allow_legacy_ghost_event_counters": False,
        "amr_counter_version": 1,
        "event_counter_version": 2,
        "restart_cache_contract_version": 1,
        "pending_event_counters_all_zero": True,
        "level_subcycling_costs_exact": True,
    }, "restart contract differs from the canonical strict contract")
    _require(policy.get("output_integrity") == {
        "minimum_closed_file_age_seconds": MIN_READY_AGE_SECONDS,
        "require_no_open_descriptors": True,
        "require_stable_size_mtime_while_hashing": True,
        "require_sha256": True,
        "refuse_manifest_overwrite": True,
    }, "output-integrity policy differs from the checker implementation")
    _require(policy.get("remote_disk") == {
        "warn_percent": 65,
        "do_not_start_percent": 75,
        "synchronized_stop_percent": 80,
        "emergency_stop_percent": 85,
        "minimum_reserve_gib": 50,
        "minimum_reserve_restart_multiples": 2,
    }, "remote-disk policy differs from the exact strict thresholds")
    _require(policy.get("yellow_event_thresholds") == [
        {"name": "fofc_per_test", "numerator": "fofc",
         "denominator": "fofc_tests", "max_ratio": 0.001,
         "consecutive_rows": 3},
        {"name": "cons_adjust_per_c2p_call", "numerator": "cons_adjust",
         "denominator": "c2p_calls", "max_ratio": 0.001,
         "consecutive_rows": 3},
        {"name": "mag_adjust_per_c2p_call", "numerator": "mag_adjust",
         "denominator": "c2p_calls", "max_ratio": 0.001,
         "consecutive_rows": 3},
    ], "yellow event policy differs from strict campaign thresholds")

    seen_blocks: set[str] = set()
    source_plan_parameters = plan.get("source", {}).get("parameters")
    _require(isinstance(source_plan_parameters, dict),
             "plan source parameter snapshot is required")
    declared_unnumbered: list[str] = []
    restart_outputs = 0
    canonical_topology_output: dict[str, Any] | None = None
    for index, output in enumerate(outputs):
        _require(isinstance(output, dict), f"outputs[{index}] must be an object")
        block = output.get("block")
        _require(isinstance(block, str) and re.fullmatch(r"output\d+", block) is not None,
                 f"outputs[{index}].block must be outputN")
        _require(block not in seen_blocks, f"duplicate output block {block}")
        seen_blocks.add(block)
        block_index = int(block[6:])
        _require(output.get("index") == block_index,
                 f"outputs[{index}].index differs from {block}")
        file_type = output.get("file_type")
        _require(isinstance(file_type, str) and bool(file_type),
                 f"outputs[{index}].file_type must be nonempty")
        if file_type == "rst":
            restart_outputs += 1
        _require(output.get("enabled") is True,
                 f"outputs[{index}] must be enabled")
        cadence_mode = output.get("cadence_mode")
        cadence = _number(output.get("cadence"), f"outputs[{index}].cadence")
        _require(cadence_mode in ("dt", "dcycle") and cadence > 0.0,
                 f"outputs[{index}] has invalid cadence")
        if cadence_mode == "dcycle":
            _require(cadence.is_integer(),
                     f"outputs[{index}] dcycle must be an integer")
        _require(isinstance(output.get("numbered"), bool),
                 f"outputs[{index}].numbered must be boolean")
        template = output.get("relative_path_template")
        required_unnumbered = output.get("required_unnumbered_paths")
        if output["numbered"]:
            _require(isinstance(template, str) and template,
                     f"outputs[{index}].relative_path_template is required")
            _validate_relative_template(template, True)
            _require(required_unnumbered in (None, []),
                     f"numbered {block} may not declare unnumbered paths")
        else:
            _require(template is None or isinstance(template, str),
                     f"unnumbered {block} has an invalid path template")
            _require(isinstance(required_unnumbered, list) and
                     bool(required_unnumbered),
                     f"unnumbered {block} must declare required_unnumbered_paths")
            for relative in required_unnumbered:
                _require(isinstance(relative, str),
                         f"{block} has a non-string unnumbered path")
                _validate_relative_template(relative, False)
                declared_unnumbered.append(relative)
            if template is not None:
                _require(template in required_unnumbered,
                         f"{block} template is absent from required unnumbered paths")
        _require(isinstance(output.get("parameters"), dict),
                 f"outputs[{index}].parameters must be an object")
        parameters = output["parameters"]
        source_output_parameters = source_plan_parameters.get(block)
        _require(isinstance(source_output_parameters, dict),
                 f"source parameter snapshot lacks {block}")
        expected_runtime_parameters = dict(source_output_parameters)
        if block == "output3":
            try:
                output3_gid = int(source_output_parameters.get("gid", "-1"))
            except ValueError as exc:
                raise CheckFailure("canonical output3 gid is invalid") from exc
            _require(
                source_output_parameters.get("file_type", "").strip() == "bin" and
                source_output_parameters.get("variable", "").strip() ==
                "mhd_divb" and output3_gid < 0 and
                "region_center" not in source_output_parameters and
                not any(f"slice_x{axis}" in source_output_parameters
                        for axis in (1, 2, 3)) and
                source_output_parameters.get(
                    "ghost_zones", "false").strip().lower() in ("false", "0") and
                "dcycle" not in source_output_parameters and
                "dt" in source_output_parameters,
                "canonical output3 is not a ghost-free full-domain mhd_divb stream")
            try:
                source_output3_dt = float(source_output_parameters["dt"])
            except ValueError as exc:
                raise CheckFailure("canonical source output3/dt is invalid") from exc
            _require(
                math.isfinite(source_output3_dt) and source_output3_dt >= 0.0 and
                (source_qualification["mode"] == "anchor_full_audit" or
                 source_output3_dt == root_dt),
                "source output3/dt must be finite/nonnegative and a chain source "
                "must already equal root_dt")
            expected_runtime_parameters["dt"] = repr(root_dt)
        if (block == restart_cadence_transition["block"] and
                restart_cadence_transition["kind"] == "tighten_v1"):
            expected_runtime_parameters["dt"] = repr(
                restart_cadence_transition["target_dt"])
        _require(parameters == expected_runtime_parameters,
                 f"{block}: planned runtime parameter snapshot differs from the "
                 "source plus exact canonical override")
        _require(parameters.get("file_type", "").strip() == file_type,
                 f"{block} file_type differs from its parameter snapshot")
        if cadence_mode == "dt":
            _require("dt" in parameters and "dcycle" not in parameters and
                     float(parameters["dt"]) == cadence,
                     f"{block} dt cadence differs from its parameter snapshot")
        else:
            _require("dcycle" in parameters and
                     int(parameters["dcycle"]) == int(cadence),
                     f"{block} dcycle differs from its parameter snapshot")
            _require("dt" not in parameters or block == "output3",
                     f"{block}: only canonical output3 may serialize dt+dcycle")
        if file_type in ("hst", "log"):
            _require(cadence_mode == "dcycle" and cadence == 1.0,
                     f"{block} {file_type} must output every root cycle")
        if output["numbered"]:
            _require(file_type in ("bin", "rst"),
                     f"{block}: strict qualifier only supports numbered bin/rst")
        if file_type == "bin":
            variable = parameters.get("variable", "").strip()
            expected_variables = _expected_binary_variables(
                variable, source_plan_parameters)
            _require(output.get("expected_binary_variables") == expected_variables,
                     f"{block}: planned binary variable labels are not canonical")
        else:
            _require("expected_binary_variables" not in output,
                     f"{block}: non-binary output declares binary variable labels")
        expected_writes = output.get("expected_writes")
        _require(isinstance(expected_writes, list),
                 f"{block}.expected_writes must be an array")
        for write_index, write in enumerate(expected_writes):
            _require(isinstance(write, dict),
                     f"{block}.expected_writes[{write_index}] must be an object")
            _integer(write.get("cycle"),
                     f"{block}.expected_writes[{write_index}].cycle")
            _number(write.get("time"),
                    f"{block}.expected_writes[{write_index}].time")
            _require(write.get("kind") in
                     ("scheduled", "forced_final", "restart_final_rewrite"),
                     f"{block}.expected_writes[{write_index}] has invalid kind")
            if output["numbered"]:
                _integer(write.get("file_number"),
                         f"{block}.expected_writes[{write_index}].file_number")
            else:
                _require("file_number" not in write,
                         f"unnumbered {block} expected write has a file number")
        expected_endpoint_state = output.get("expected_endpoint_state")
        _require(isinstance(expected_endpoint_state, dict),
                 f"{block}.expected_endpoint_state must be an object")
        _integer(expected_endpoint_state.get("file_number"),
                 f"{block}.expected_endpoint_state.file_number")
        _number(expected_endpoint_state.get("last_time"),
                f"{block}.expected_endpoint_state.last_time")
        _integer(expected_endpoint_state.get("last_write_cycle"),
                 f"{block}.expected_endpoint_state.last_write_cycle")
        if block == "output3":
            canonical_topology_output = output
    _require(restart_outputs == 1,
             "strict segment plan must declare exactly one restart output")
    replayed_schedules, replayed_endpoint_states = _replay_output_plan(
        outputs, source_plan_parameters, expected)
    for output in outputs:
        block = output["block"]
        _require(output["expected_writes"] == replayed_schedules[block],
                 f"{block}: plan expected_writes differs from independent "
                 "source-phase cadence replay")
        _require(output["expected_endpoint_state"] ==
                 replayed_endpoint_states[block],
                 f"{block} endpoint cadence state: plan expected_endpoint_state "
                 "differs from independent source-phase cadence replay")
    _require("output3" in seen_blocks,
             "strict segment plan lacks canonical output3 topology reference")
    _require(canonical_topology_output is not None and
             canonical_topology_output["file_type"] == "bin" and
             canonical_topology_output["cadence_mode"] == "dt" and
             canonical_topology_output["cadence"] == root_dt and
             canonical_topology_output["parameters"].get("variable") ==
             "mhd_divb" and
             [write["cycle"] for write in
              canonical_topology_output["expected_writes"]] ==
             list(range(source_cycle + 1, final_cycle + 1)) and
             all(write["kind"] == "scheduled" for write in
                 canonical_topology_output["expected_writes"]),
             "canonical output3 must provide exactly one scheduled full-domain "
             "divB topology reference per root cycle")
    source_output3 = source_plan_parameters["output3"]
    try:
        source_file_number = int(source_output3.get("file_number", "0"))
        source_last_time = float(source_output3.get("last_time", "-1"))
    except ValueError as exc:
        raise CheckFailure("source output3 cadence state is invalid") from exc
    replayed_last_time = source_last_time
    for write in canonical_topology_output["expected_writes"]:
        replayed_last_time = (write["time"] if replayed_last_time < 0.0
                              else replayed_last_time + root_dt)
    _require(canonical_topology_output["expected_endpoint_state"] == {
        "file_number": source_file_number + root_steps,
        "last_time": replayed_last_time,
        "last_write_cycle": final_cycle,
    }, "canonical output3 endpoint cadence state differs from exact root-step replay")
    required_files = plan.get("required_files")
    _require(isinstance(required_files, list) and all(
        isinstance(item, str) for item in required_files
    ), "plan.required_files must be an array of relative paths")
    for relative in required_files:
        _validate_relative_template(relative, False)
    _require(len(required_files) == len(set(required_files)),
             "plan.required_files contains duplicates")
    _require(sorted(required_files) == sorted(declared_unnumbered),
             "plan.required_files must equal all declared unnumbered output paths")
    return {
        "source_cycle": source_cycle,
        "final_cycle": final_cycle,
        "source_time": source_time,
        "final_time": final_time,
        "root_steps": root_steps,
        "root_dt": root_dt,
        "tlim": tlim,
        "ranks": ranks,
        "gpu_exit_memory_mib_max": memory_limit,
        "endpoint_time_ulp_tolerance": ulps,
    }


def _validate_relative_template(template: str, numbered: bool) -> None:
    _require("\x00" not in template, "output path template contains NUL")
    fields = re.findall(r"\{([^{}]+)\}", template)
    if numbered:
        _require(fields == ["file_number:05d"],
                 "numbered output template must contain exactly {file_number:05d}")
    else:
        _require(not fields and "{" not in template and "}" not in template,
                 "unnumbered output template may not contain replacement fields")
    probe = template.format(file_number=0) if numbered else template
    path = Path(probe)
    _require(not path.is_absolute() and all(part not in ("", ".", "..")
                                            for part in path.parts),
             f"unsafe output relative path template: {template}")


def _expected_binary_variables(variable: str,
                               parameters: dict[str, dict[str, str]]) -> list[str]:
    if variable == "mhd_divb":
        return ["divb"]
    if variable == "mhd_gr_diagnostics":
        return [
            "gr_bsq", "gr_lorentz", "gr_sigma", "gr_beta_inv",
            "gr_excision_mask",
        ]
    if variable == "mhd_w_bcc":
        mhd = parameters.get("mhd", {})
        adm = parameters.get("adm", {})
        try:
            nscalars = int(mhd.get("nscalars", "0"))
        except ValueError as exc:
            raise CheckFailure("strict mhd_w_bcc output has invalid mhd/nscalars") \
                from exc
        dynamic = adm.get("dynamic", "false").strip().lower() in ("true", "1")
        _require(mhd.get("eos", "").strip() == "ideal" and
                 nscalars == 0 and dynamic,
                 "strict BBH mhd_w_bcc requires ideal EOS, zero scalars, dynamic GR")
        return [
            "dens", "velx", "vely", "velz", "press",
            "bcc1", "bcc2", "bcc3", "temperature",
        ]
    raise CheckFailure(
        f"unsupported strict BBH binary variable group: {variable!r}")


def _parse_exit_status(path: Path) -> int:
    raw = _read_stable_bytes(path)
    try:
        text = raw.decode("ascii").strip()
    except UnicodeDecodeError as exc:
        raise CheckFailure(f"{path}: exit status is not ASCII") from exc
    _require(re.fullmatch(r"[+-]?\d+", text) is not None,
             f"{path}: exit status must contain exactly one integer")
    return int(text)


def audit_run_log(path: Path, expected: dict[str, Any]) -> dict[str, Any]:
    try:
        lines = _read_stable_bytes(path).decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise CheckFailure(f"{path}: run log is not UTF-8") from exc
    terminations = [TERMINATION_RE.fullmatch(line.strip()) for line in lines]
    terminations = [match.group(1) for match in terminations if match is not None]
    _require(terminations == ["cycle"],
             f"run log must contain exactly one cycle-limit termination; got {terminations}")

    finals = [FINAL_STATE_RE.fullmatch(line.strip()) for line in lines]
    finals = [match for match in finals if match is not None]
    _require(len(finals) == 1, "run log must contain exactly one final time/cycle line")
    final_time, final_cycle = float(finals[0].group(1)), int(finals[0].group(2))
    _require(final_cycle == expected["final_cycle"],
             f"run-log final cycle {final_cycle} != {expected['final_cycle']}")
    _require(final_time == float(f"{expected['final_time']:.6e}"),
             f"run-log final time {final_time!r} != expected {expected['final_time']!r}")

    limits = [LIMIT_STATE_RE.fullmatch(line.strip()) for line in lines]
    limits = [match for match in limits if match is not None]
    _require(len(limits) == 1, "run log must contain exactly one tlim/nlim line")
    logged_tlim, logged_nlim = float(limits[0].group(1)), int(limits[0].group(2))
    _require(logged_nlim == expected["final_cycle"],
             "run-log nlim does not equal the planned final cycle")
    _require(logged_tlim == float(f"{expected['tlim']:.6e}"),
             "run-log tlim differs from the plan's diagnostic representation")

    diagnostics = [DIAGNOSTIC_RE.fullmatch(line.strip()) for line in lines]
    diagnostics = [match for match in diagnostics if match is not None]
    expected_cycles = list(range(
        expected["source_cycle"], expected["final_cycle"] + 1))
    observed_cycles = [int(match.group(2)) for match in diagnostics]
    _require(observed_cycles == expected_cycles,
             "run-log diagnostics must contain every root cycle from source through "
             "the Finalize endpoint exactly once and in order")
    elapsed_values = [float(match.group(1)) for match in diagnostics]
    _require(all(math.isfinite(value) and value >= 0.0 for value in elapsed_values) and
             all(right >= left for left, right in
                 zip(elapsed_values, elapsed_values[1:])),
             "run-log elapsed times must be finite, nonnegative, and monotone")
    expected_times: list[float] = []
    current_time = expected["source_time"]
    expected_times.append(current_time)
    for _ in range(expected["root_steps"]):
        current_time += expected["root_dt"]
        expected_times.append(current_time)
    # OutputCycleDiagnostics intentionally uses scientific setprecision(6).  Compare to
    # that exact public representation instead of pretending the text retains all Real8
    # bits; the endpoint restart and plan independently bind the full-precision value.
    printed_times = [float(f"{value:.6e}") for value in expected_times]
    observed_times = [float(match.group(3)) for match in diagnostics]
    observed_dts = [float(match.group(4)) for match in diagnostics]
    printed_dt = float(f"{expected['root_dt']:.6e}")
    _require(observed_times == printed_times,
             "run-log diagnostic times do not follow sequential fixed root_dt steps")
    _require(observed_dts == [printed_dt] * len(diagnostics),
             "run-log diagnostics contain a non-fixed root dt")

    caches = [CACHE_RE.fullmatch(line.strip()) for line in lines]
    caches = [match for match in caches if match is not None]
    _require(len(caches) == 1,
             "run log must contain exactly one anchored cache-qualification line")
    solver_failures, nonfinite = int(caches[0].group(1)), int(caches[0].group(2))
    raw_relative, absolute, mixed, tolerance = map(float, caches[0].groups()[2:])
    _require(solver_failures == 0 and nonfinite == 0,
             "restart cache qualification reported integer failures")
    _require(all(math.isfinite(value) and value >= 0.0 for value in
                 (raw_relative, absolute, mixed, tolerance)),
             "restart cache qualification contains invalid floating-point values")
    _require(tolerance > 0.0 and mixed <= tolerance,
             f"restart cache mixed change {mixed} exceeds tolerance {tolerance}")
    return {
        "termination": "cycle_limit",
        "time": final_time,
        "cycle": final_cycle,
        "tlim": logged_tlim,
        "nlim": logged_nlim,
        "root_step_diagnostics": {
            "rows": len(diagnostics),
            "cycle_min": observed_cycles[0],
            "cycle_max": observed_cycles[-1],
            "fixed_dt": expected["root_dt"],
            "all_root_cycles_present": True,
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


def _same_planned_file(actual: Any, planned: Any) -> bool:
    return isinstance(actual, dict) and isinstance(planned, dict) and all(
        actual.get(name) == planned.get(name) for name in ("path", "size", "sha256")
    )


def _same_file_content(actual: Any, planned: Any) -> bool:
    return isinstance(actual, dict) and isinstance(planned, dict) and all(
        actual.get(name) == planned.get(name) for name in ("size", "sha256")
    )


def _audit_directory_transport_record(
        actual: Any, planned: dict[str, Any], holder_pid: int,
        holder_start_time_ticks: int) -> dict[str, Any]:
    _require(isinstance(actual, dict) and
             set(actual) == {"kind", "holder_pid", "holder_start_time_ticks",
                             "roles"} and
             actual.get("kind") == DIRECTORY_TRANSPORT_KIND and
             actual.get("holder_pid") == holder_pid and
             actual.get("holder_start_time_ticks") == holder_start_time_ticks,
             "directory-holder process identity differs from input holder")
    roles = actual.get("roles")
    planned_roles = planned.get("roles")
    _require(isinstance(roles, dict) and isinstance(planned_roles, dict) and
             set(roles) == {"state_dir", "evidence_dir"},
             "directory-holder roles are not exact")
    for role, descriptor in (("state_dir", STATE_DIRECTORY_FD),
                             ("evidence_dir", EVIDENCE_DIRECTORY_FD)):
        expected = planned_roles[role]
        row = roles[role]
        _require(isinstance(row, dict) and
                 set(row) == {"path", "device", "inode", "owner_uid", "mode",
                              "fd", "proc_path", "role", "access_mode"} and
                 row.get("role") == role and row.get("fd") == descriptor and
                 row.get("path") == expected["planned_path"] and
                 row.get("proc_path") == f"/proc/{holder_pid}/fd/{descriptor}" and
                 row.get("access_mode") == "read_only_directory_descriptor" and
                 isinstance(row.get("device"), int) and
                 isinstance(row.get("inode"), int) and
                 row.get("owner_uid") == os.geteuid() and
                 isinstance(row.get("mode"), str) and
                 re.fullmatch(r"0[0-7]{3}", row["mode"]) is not None and
                 not (int(row["mode"], 8) & 0o022),
                 f"directory-holder role {role} is invalid")
        try:
            path_info = Path(row["path"]).lstat()
        except OSError as exc:
            raise CheckFailure(
                f"cannot audit held {role} directory path: {exc}") from exc
        _require(stat.S_ISDIR(path_info.st_mode) and
                 not stat.S_ISLNK(path_info.st_mode) and
                 (path_info.st_dev, path_info.st_ino, path_info.st_uid,
                  stat.S_IMODE(path_info.st_mode)) ==
                 (row["device"], row["inode"], row["owner_uid"],
                  int(row["mode"], 8)),
                 f"held {role} directory path was replaced")
    return actual


def _audit_executable_transport_record(
        actual: Any, planned: dict[str, Any], holder_pid: int,
        holder_start_time_ticks: int,
        plan: dict[str, Any]) -> dict[str, Any]:
    _require(isinstance(actual, dict) and
             set(actual) == {"kind", "holder_pid", "holder_start_time_ticks",
                             "roles"} and
             actual.get("kind") == EXECUTABLE_TRANSPORT_KIND and
             actual.get("holder_pid") == holder_pid and
             actual.get("holder_start_time_ticks") == holder_start_time_ticks,
             "executable-holder process identity differs from input holder")
    roles = actual.get("roles")
    planned_roles = planned.get("roles")
    _require(isinstance(roles, dict) and isinstance(planned_roles, dict) and
             set(roles) == {"launcher", "binary"},
             "executable-holder roles are not exact")
    expected_parents = {
        "launcher": plan["launch_contract"]["launcher"],
        "binary": plan["inputs"]["binary"],
    }
    for role, descriptor in (("launcher", LAUNCHER_EXECUTABLE_FD),
                             ("binary", BINARY_EXECUTABLE_FD)):
        row = roles[role]
        expected_contract = planned_roles[role]
        expected_parent = expected_parents[role]
        _require(isinstance(row, dict) and
                 set(row) == {"path", "device", "inode", "size", "mtime_ns",
                              "ctime_ns", "sha256", "closure_check", "role", "fd",
                              "proc_path", "access_mode"} and
                 row.get("role") == role and row.get("fd") == descriptor and
                 row.get("path") == expected_parent["path"] and
                 row.get("proc_path") == f"/proc/{holder_pid}/fd/{descriptor}" and
                 row.get("access_mode") == "read_only" and
                 row.get("closure_check") == "fixed_read_only_descriptor" and
                 _same_file_content(row, expected_parent) and
                 row.get("path") == expected_contract["parent_path"] and
                 all(isinstance(row.get(name), int) for name in (
                     "device", "inode", "size", "mtime_ns", "ctime_ns")),
                 f"executable-holder role {role} differs from planned bytes")
    return actual


def audit_launch_record(path: Path, plan: dict[str, Any], plan_path: Path,
                        state_dir: Path,
                        gpu_before_binding: dict[str, Any] | None = None
                        ) -> dict[str, Any]:
    """Bind immutable, tokenized launch evidence to every execution-critical input."""

    _require_immutable_ready(path, ".launch.ready", "launch record")
    record, raw = _load_json(path)
    _require(record.get("schema") == SCHEMA and
             record.get("kind") == "athenak_segment_launch" and
             record.get("status") == "ready",
             "launch record has an unsupported schema/kind/status")
    athena_argv = record.get("athena_argv")
    mpi_argv = record.get("mpi_argv")
    _require(isinstance(athena_argv, list) and athena_argv and
             all(isinstance(token, str) and token and "\x00" not in token
                 for token in athena_argv),
             "launch athena_argv must be a nonempty token array")
    _require(isinstance(mpi_argv, list) and len(mpi_argv) == 6 and
             all(isinstance(token, str) and token and "\x00" not in token
                 for token in mpi_argv),
             "launch mpi_argv must be the canonical six-token MPI invocation")
    launcher = record.get("launcher")
    _require(isinstance(launcher, dict), "launch record must bind the MPI launcher")
    _validate_planned_file_record(launcher, "launch.launcher")
    launcher_binding = _verify_planned_file(launcher, "MPI launcher")
    _require(mpi_argv == [launcher["path"], "--allow-run-as-root", "--bind-to",
                          "none", "-np", str(plan["policy"]["ranks"])],
             "launch mpi_argv does not encode exactly the planned single application")
    wall_limit = plan.get("launch_contract", {}).get("wall_time_seconds")
    _require(isinstance(wall_limit, int) and not isinstance(wall_limit, bool) and
             wall_limit > 0,
             "plan launch_contract must bind a positive integer wall_time_seconds")
    inputs = plan["inputs"]
    launch_contract = plan.get("launch_contract")
    transport_contract = launch_contract.get("input_transport", {})
    _require(record.get("input_transport_contract") == transport_contract,
             "launch record input-transport contract differs from the plan")
    launch_transport = record.get("input_transport")
    _require(isinstance(launch_transport, dict) and
             launch_transport.get("kind") == INPUT_TRANSPORT_KIND and
             isinstance(launch_transport.get("holder_pid"), int) and
             launch_transport["holder_pid"] > 0 and
             isinstance(launch_transport.get("holder_start_time_ticks"), int) and
             launch_transport["holder_start_time_ticks"] > 0,
             "launch record lacks a fixed-descriptor input-holder identity")
    holder_pid = launch_transport["holder_pid"]
    holder_start_time_ticks = launch_transport["holder_start_time_ticks"]
    directory_transport = _audit_directory_transport_record(
        record.get("directory_transport"),
        launch_contract.get("directory_transport", {}),
        holder_pid, holder_start_time_ticks)
    _require(record.get("directory_transport_contract") ==
             launch_contract.get("directory_transport"),
             "launch directory-transport contract differs from the plan")
    executable_transport = _audit_executable_transport_record(
        record.get("executable_transport"),
        launch_contract.get("executable_transport", {}),
        holder_pid, holder_start_time_ticks, plan)
    _require(record.get("executable_transport_contract") ==
             launch_contract.get("executable_transport"),
             "launch executable-transport contract differs from the plan")
    expected_athena = [
        token.replace(HOLDER_PID_TOKEN, str(holder_pid))
        for token in launch_contract.get("athena_argv_template", [])
    ]
    _require(isinstance(launch_contract, dict) and
             launch_contract.get("mpi_argv") == mpi_argv and
             _same_planned_file(launcher, launch_contract.get("launcher", {})) and
             launch_contract.get("state_dir") == str(state_dir.resolve(strict=True)) and
             launch_contract.get("world_size") == plan["policy"]["ranks"] and
             launch_contract.get("gpu_count") == plan["policy"]["ranks"],
             "plan launch_contract does not contain the canonical MPI/Athena argv")
    _require(athena_argv == expected_athena,
             "launch record Athena argv differs from the exact planned argv")
    _require(record.get("world_size") == plan["policy"]["ranks"],
             "launch record MPI world size differs from planned ranks")
    _require(record.get("gpu_count") == plan["policy"]["ranks"],
             "launch record GPU count differs from planned ranks")
    _require(record.get("plan_path") == str(plan_path.resolve(strict=True)) and
             launch_contract.get("plan_path") == str(plan_path.resolve(strict=True)),
             "launch record plan path differs from the checked plan")
    _require(record.get("plan_sha256") == hashlib.sha256(
        _read_stable_bytes(plan_path)).hexdigest(),
        "launch record plan hash differs from the checked plan")
    _require(record.get("state_dir") == str(state_dir.resolve(strict=True)),
             "launch record state directory differs from the checked state directory")
    _require(record.get("binary_sha256") == inputs["binary"]["sha256"] and
             record.get("source_restart_sha256") ==
             inputs["source_restart"]["sha256"] and
             record.get("trajectory_sha256") == inputs["trajectory"]["sha256"],
             "launch record input hashes differ from the plan")
    launch_argv = record.get("launch_argv")
    _require(launch_argv == [*mpi_argv, *athena_argv],
             "launch record does not bind the exact shell-free launch argv")
    mpirun_pid = record.get("mpirun_pid")
    _require(isinstance(mpirun_pid, int) and not isinstance(mpirun_pid, bool) and
             mpirun_pid > 0 and
             isinstance(record.get("mpirun_start_time_ticks"), int) and
             record["mpirun_start_time_ticks"] > 0 and
             record.get("mpirun_cmdline") == launch_argv and
             _same_file_content(record.get("mpirun_executable"), launcher),
             "launch record lacks canonical live mpirun process proof")
    _require(holder_pid != mpirun_pid,
             "input holder and mpirun must be distinct process identities")
    managed_process_group = record.get("managed_process_group")
    _require(managed_process_group == {
        "pgid": mpirun_pid,
        "new_session": True,
        "failure_cleanup": "SIGTERM_then_SIGKILL_with_quiescence_proof",
    }, "launch record lacks the exact managed process-group contract")
    proc_access_probe = record.get("proc_access_probe")
    _require(isinstance(proc_access_probe, dict) and
             set(proc_access_probe) == {
                 "kind", "peer_pid", "holder_pid", "families",
                 "all_reopened_and_sampled", "executable_access",
             } and
             proc_access_probe.get("kind") == "fork_peer_procfs_reopen_v1" and
             isinstance(proc_access_probe.get("peer_pid"), int) and
             not isinstance(proc_access_probe.get("peer_pid"), bool) and
             proc_access_probe["peer_pid"] > 0 and
             proc_access_probe.get("holder_pid") == holder_pid and
             proc_access_probe.get("families") == {
                 "input": ["source_restart", "trajectory"],
                 "directory": ["evidence_dir", "state_dir"],
                 "executable": ["binary", "launcher"],
             } and proc_access_probe.get("executable_access") == {
                 "binary": "effective_uid_x_ok",
                 "launcher": "effective_uid_x_ok",
             } and proc_access_probe.get("all_reopened_and_sampled") is True,
             "launch record lacks exact peer proc-holder access proof")
    original_inputs = record.get("original_inputs")
    staged_inputs = record.get("staged_inputs")
    staging_directory = record.get("staging_directory")
    _require(isinstance(original_inputs, dict) and
             set(original_inputs) == {"source_restart", "trajectory"} and
             all(_same_planned_file(original_inputs[name], inputs[name])
                 for name in original_inputs),
             "launch record original-input bindings differ from the plan")
    _require(isinstance(staged_inputs, dict) and
             set(staged_inputs) == {"source_restart", "trajectory"} and
             isinstance(staging_directory, dict) and
             staging_directory.get("path") == transport_contract.get("staging_dir") and
             staging_directory.get("mode") == "0555" and
             all(isinstance(staging_directory.get(name), int)
                 for name in ("device", "inode")),
             "launch record staging bindings are invalid")
    disk_contract = _validate_disk_preflight_contract(
        launch_contract.get("disk_preflight"), inputs["source_restart"]["size"],
        inputs["trajectory"]["size"])
    disk_evidence = record.get("disk_preflight")
    _require(isinstance(disk_evidence, dict) and set(disk_evidence) == {
        "contract", "before_staging", "before_spawn",
    } and disk_evidence.get("contract") == disk_contract,
             "launch record disk-preflight contract differs from the plan")
    _validate_disk_preflight_snapshot(
        disk_evidence["before_staging"], disk_contract, "before_staging",
        state_dir=str(state_dir.resolve(strict=True)),
        staging_dir=str(Path(transport_contract["staging_dir"])),
        expected_devices={
            "state_dir": directory_transport["roles"]["state_dir"]["device"],
            "staging_dir": staging_directory["device"],
        },
    )
    _validate_disk_preflight_snapshot(
        disk_evidence["before_spawn"], disk_contract, "before_spawn",
        state_dir=str(state_dir.resolve(strict=True)),
        staging_dir=str(Path(transport_contract["staging_dir"])),
        expected_devices={
            "state_dir": directory_transport["roles"]["state_dir"]["device"],
            "staging_dir": staging_directory["device"],
        },
    )
    planned_roles = transport_contract.get("roles", {})
    _require(isinstance(planned_roles, dict) and
             set(planned_roles) == {"source_restart", "trajectory"},
             "planned input-holder roles are not exact")
    held_roles = launch_transport.get("roles")
    _require(isinstance(held_roles, dict) and
             set(held_roles) == {"source_restart", "trajectory"},
             "launch input-holder roles are not exact")
    for role, descriptor in (("source_restart", SOURCE_RESTART_FD),
                             ("trajectory", TRAJECTORY_FD)):
        planned_staged = planned_roles[role]["staged_file"]
        staged = staged_inputs[role]
        held = held_roles[role]
        _require(_same_planned_file(staged, planned_staged) and
                 _same_planned_file(held, planned_staged) and
                 held.get("role") == role and held.get("fd") == descriptor and
                 held.get("proc_path") == f"/proc/{holder_pid}/fd/{descriptor}" and
                 held.get("access_mode") == "read_only" and
                 all(held.get(name) == staged.get(name) for name in (
                     "device", "inode", "size", "mtime_ns", "ctime_ns", "sha256",
                 )),
                 f"launch held descriptor for {role} differs from the staged bytes")
    hostname = record.get("hostname")
    _require(isinstance(hostname, str) and hostname,
             "launch record hostname must be nonempty")
    ranks = record.get("ranks")
    world_size = plan["policy"]["ranks"]
    _require(isinstance(ranks, list) and len(ranks) == world_size,
             "launch record must prove exactly one live process per MPI rank")
    gpu_proof = record.get("gpu_before")
    _require(isinstance(gpu_proof, dict) and
             isinstance(gpu_proof.get("path"), str) and
             isinstance(gpu_proof.get("sha256"), str) and
             isinstance(gpu_proof.get("devices"), list),
             "launch record lacks GPU-before inventory proof")
    if gpu_before_binding is not None:
        _require(gpu_proof["path"] == gpu_before_binding["path"] and
                 gpu_proof["sha256"] == gpu_before_binding["sha256"],
                 "launch record GPU-before binding differs from checked evidence")
    devices = gpu_proof["devices"]
    _require(len(devices) == world_size,
             "launch record GPU-before inventory has the wrong size")
    device_by_ordinal: dict[int, dict[str, Any]] = {}
    for ordinal, device in enumerate(devices):
        _require(isinstance(device, dict) and
                 isinstance(device.get("index"), int) and
                 device.get("cuda_ordinal") == ordinal and
                 isinstance(device.get("uuid"), str) and
                 device["uuid"].startswith("GPU-") and
                 isinstance(device.get("pci_bus_id"), str) and
                 PCI_BUS_ID_RE.fullmatch(device["pci_bus_id"]) is not None and
                 device.get("uncorrected_ecc") == 0 and
                 device.get("corrected_ecc") == 0 and
                 isinstance(device.get("memory_total_mib"), int) and
                 device["memory_total_mib"] > 0 and
                 isinstance(device.get("memory_used_mib"), int) and
                 0 <= device["memory_used_mib"] <= device["memory_total_mib"],
                 "launch record GPU-before device proof is invalid")
        device_by_ordinal[ordinal] = device
    _require(len({device["index"] for device in devices}) == world_size and
             len({device["uuid"] for device in devices}) == world_size and
             len({device["pci_bus_id"].lower() for device in devices}) == world_size,
             "launch GPU inventory contains duplicate physical identities")
    _require([_pci_bus_sort_key(device["pci_bus_id"]) for device in devices] ==
             sorted(_pci_bus_sort_key(device["pci_bus_id"])
                    for device in devices),
             "launch GPU inventory is not PCI ordered")
    capacity_transition = _validate_capacity_transition(plan)
    capacity_preflight = _gpu_capacity_preflight_evidence(
        devices, capacity_transition)
    _require(record.get("capacity_transition") == capacity_transition and
             record.get("gpu_capacity_preflight") == capacity_preflight,
             "launch record does not prove the plan-bound GPU capacity gate")
    seen_pids: set[int] = set()
    seen_uuids: set[str] = set()
    seen_local_ranks: set[int] = set()
    for global_rank, rank in enumerate(ranks):
        _require(isinstance(rank, dict), "launch rank proof must be an object")
        local_rank = rank.get("local_rank")
        pid = rank.get("pid")
        environment = rank.get("mpi_environment")
        _require(rank.get("global_rank") == global_rank and
                 isinstance(local_rank, int) and not isinstance(local_rank, bool) and
                 0 <= local_rank < world_size and local_rank not in seen_local_ranks and
                 isinstance(pid, int) and not isinstance(pid, bool) and pid > 0 and
                 isinstance(rank.get("start_time_ticks"), int) and
                 rank["start_time_ticks"] > 0 and
                 pid != mpirun_pid and pid != holder_pid and pid not in seen_pids and
                 rank.get("hostname") == hostname and
                 rank.get("gpu_cuda_ordinal") == local_rank and
                 isinstance(rank.get("gpu_index"), int) and
                 rank.get("gpu_uuid") ==
                 device_by_ordinal[local_rank]["uuid"] and
                 rank.get("gpu_index") == device_by_ordinal[local_rank]["index"] and
                 rank.get("gpu_pci_bus_id") ==
                 device_by_ordinal[local_rank]["pci_bus_id"] and
                 rank.get("gpu_uuid") not in seen_uuids and
                 rank.get("cmdline") == athena_argv and
                 _same_file_content(rank.get("executable"), inputs["binary"]),
                 f"launch record rank {global_rank} process/GPU proof is invalid")
        _validate_rank_environment(
            environment, global_rank=global_rank, local_rank=local_rank,
            world_size=world_size, launcher_pid=mpirun_pid, hostname=hostname,
            state_dir=str(state_dir.resolve(strict=True)),
            athena_argv=athena_argv)
        seen_pids.add(pid)
        seen_uuids.add(rank["gpu_uuid"])
        seen_local_ranks.add(local_rank)
    _require(seen_local_ranks == set(range(world_size)),
             "launch record local MPI ranks are not exactly 0..world_size-1")
    _require(seen_uuids == {device["uuid"] for device in devices},
             "launch record rank GPU contexts do not cover the planned inventory")
    _require(len({rank["mpi_environment"]["PMIX_SERVER_URI21"]
                  for rank in ranks}) == 1,
             "launch record ranks do not share one PMIx server URI")
    _require(record.get("gpu_mapping_basis") ==
             "kokkos_mpi_rank_token_plus_ompi_local_rank_plus_"
             "nvidia_compute_context_uuid",
             "launch record lacks the required observed GPU mapping proof")
    launcher_tool = record.get("launcher_tool")
    planned_launcher_tool = plan.get("tools", {}).get("segment_launcher")
    _validate_planned_file_record(launcher_tool, "launch.launcher_tool")
    _require(_same_planned_file(launcher_tool, planned_launcher_tool),
             "launch record launcher tool differs from the plan-bound tool")
    launcher_tool_binding = _verify_planned_file(
        planned_launcher_tool, "segment launcher tool")
    nvidia_smi = record.get("nvidia_smi")
    _validate_planned_file_record(nvidia_smi, "launch.nvidia_smi")
    _require(_same_planned_file(nvidia_smi, plan["tools"]["nvidia_smi"]),
             "launch nvidia-smi tool differs from the plan-bound executable")
    nvidia_smi_binding = _verify_planned_file(
        plan["tools"]["nvidia_smi"], "plan-bound nvidia-smi executable")
    _require(os.access(Path(nvidia_smi_binding["path"]), os.X_OK),
             "plan-bound nvidia-smi is not executable at audit time")
    launch_environment = record.get("launch_environment")
    _require(launch_environment == launch_contract.get("environment") and
             _canonical_launch_environment(launch_environment) ==
             launch_contract["environment"]["values"],
             "launch environment evidence differs from the exact plan contract")
    mca_contract = launch_contract.get("mca_configuration")
    mca_current = _audit_mca_configuration(
        mca_contract, launch_contract["environment"]["values"]["HOME"])
    _require(record.get("mca_configuration_contract") == mca_contract,
             "launch MCA configuration contract differs from the plan")
    mca_evidence = record.get("mca_configuration")
    _require(isinstance(mca_evidence, dict) and
             set(mca_evidence) == {"preflight", "before_spawn", "at_rank_proof"}
             and all(mca_evidence[stage] == mca_contract
                     for stage in ("preflight", "before_spawn", "at_rank_proof"))
             and mca_current == mca_contract,
             "launch MCA configuration snapshots do not close the plan binding")
    execution_tools_at_launch = record.get("execution_tools_at_launch")
    _require(isinstance(execution_tools_at_launch, dict) and
             set(execution_tools_at_launch) == {
                 "segment_launcher", "segment_checker", "output_integrity",
                 "restart_auditor", "restart_metadata_reader", "nvidia_smi",
             } and all(
                 _same_planned_file(execution_tools_at_launch[name],
                                    plan["tools"][name])
                 for name in execution_tools_at_launch),
             "launch execution-tool proof set differs from plan-bound executables")
    expected_repository = plan["inputs"]["repo"]
    repository_preflight = record.get("repository_preflight")
    repository_at_launch = record.get("repository_at_launch")
    _require(isinstance(repository_preflight, dict) and
             repository_at_launch == repository_preflight and
             repository_preflight.get("path") == expected_repository["path"] and
             repository_preflight.get("commit") == expected_repository["commit"] and
             repository_preflight.get("clean") is True and
             _same_planned_file(repository_preflight.get("git_tool"),
                                plan["tools"]["git"]) and
             repository_preflight.get("git_environment_policy") ==
             "explicit_clean_environment_v1" and
             repository_preflight.get("git_environment") == GIT_ENVIRONMENT and
             repository_preflight.get("git_configuration") == GIT_CONFIGURATION,
             "launch repository proofs differ from the plan-bound clean repository")
    _require(record.get("gpu_visibility_environment") == {
        "CUDA_VISIBLE_DEVICES": None,
        "KOKKOS_VISIBLE_DEVICES": None,
        "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
    }, "launch record GPU visibility environment is not canonical")
    forbidden = re.compile(
        r"^(?:output\d+/(?:dt|dcycle|file_number|last_time|last_write_cycle)|"
        r"job/basename|problem/trajectory_file|mesh(?:block|_refinement)?/[^=]+)="
    )
    permitted = {
        f"problem/trajectory_file=/proc/{holder_pid}/fd/{TRAJECTORY_FD}",
        f"output3/dt={plan['expected']['root_dt']!r}",
    }
    restart_cadence_transition = _validate_restart_cadence_transition(
        plan, float(plan["expected"]["root_dt"]))
    if restart_cadence_transition["kind"] == "tighten_v1":
        permitted.add(restart_cadence_transition["runtime_override"])
    capacity_transition = _validate_capacity_transition(plan)
    if capacity_transition["kind"] == "increase_v1":
        permitted.add(capacity_transition["runtime_override"])
    _require(not any(forbidden.match(token) and token not in permitted
                     for token in [*mpi_argv, *athena_argv]),
             "launch record contains a forbidden output/provenance override")
    mesh_tokens = [token for token in athena_argv
                   if token.startswith(("mesh/", "meshblock/",
                                        "mesh_refinement/"))]
    _require(mesh_tokens == ([] if capacity_transition["kind"] == "unchanged_v1"
                             else [capacity_transition["runtime_override"]]),
             "launch record contains an extra Mesh parameter")
    return {
        "sha256": hashlib.sha256(raw).hexdigest(),
        "athena_argv": athena_argv,
        "mpi_argv": mpi_argv,
        "launcher": launcher_binding,
        "world_size": record["world_size"],
        "gpu_count": record["gpu_count"],
        "state_dir": record["state_dir"],
        "plan_sha256": record["plan_sha256"],
        "launch_argv": launch_argv,
        "mpirun_pid": mpirun_pid,
        "mpirun_start_time_ticks": record["mpirun_start_time_ticks"],
        "managed_process_group": managed_process_group,
        "proc_access_probe": proc_access_probe,
        "hostname": hostname,
        "ranks": ranks,
        "gpu_before": gpu_proof,
        "capacity_transition": capacity_transition,
        "gpu_capacity_preflight": record["gpu_capacity_preflight"],
        "launcher_tool": launcher_tool_binding,
        "input_transport": launch_transport,
        "original_inputs": original_inputs,
        "staged_inputs": staged_inputs,
        "staging_directory": staging_directory,
        "disk_preflight": disk_evidence,
        "repository_preflight": repository_preflight,
        "repository_at_launch": repository_at_launch,
        "nvidia_smi": nvidia_smi_binding,
        "launch_environment": launch_environment,
        "directory_transport": directory_transport,
        "executable_transport": executable_transport,
        "execution_tools_at_launch": execution_tools_at_launch,
    }


def _normalize_event_thresholds(value: Any) -> list[dict[str, Any]]:
    if isinstance(value, dict):
        _require(value.get("hard_equal_zero") == ["eos_fail", "eos_vceil"],
                 "event hard-zero policy must be exactly eos_fail/eos_vceil")
        c2p_max = _integer(value.get("c2p_iterations_exclusive_max"),
                           "event_thresholds.c2p_iterations_exclusive_max")
        _require(c2p_max > 0, "C2P iteration maximum must be positive")
        value = value.get("hard_ratios")
    _require(isinstance(value, list) and value,
             "event hard-ratio thresholds must be a nonempty array")
    result: list[dict[str, Any]] = []
    for index, item in enumerate(value):
        _require(isinstance(item, dict), f"event threshold {index} must be an object")
        numerator = item.get("numerator")
        denominator = item.get("denominator")
        maximum = _number(item.get("max_ratio"),
                          f"event_thresholds[{index}].max_ratio")
        _require(numerator in EVENT_COLUMNS[1:] and denominator in EVENT_COLUMNS[1:],
                 f"event threshold {index} uses an unsupported event column")
        _require(maximum >= 0.0, f"event threshold {index} must be nonnegative")
        result.append({
            "name": item.get("name"),
            "numerator": numerator,
            "denominator": denominator,
            "max_ratio": maximum,
        })
    return result


def _normalize_yellow_event_thresholds(value: Any) -> list[dict[str, Any]]:
    """Validate non-fatal ratio warnings without weakening the hard rules."""

    _require(isinstance(value, list) and value,
             "yellow event thresholds must be a nonempty array")
    normalized = _normalize_event_thresholds(value)
    result: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, (raw, rule) in enumerate(zip(value, normalized)):
        name = rule.get("name")
        consecutive = _integer(raw.get("consecutive_rows"),
                               f"yellow_event_thresholds[{index}].consecutive_rows")
        _require(isinstance(name, str) and bool(name) and name not in seen,
                 f"yellow event threshold {index} has an invalid or duplicate name")
        _require(consecutive > 0,
                 f"yellow event threshold {index} consecutive_rows must be positive")
        seen.add(name)
        result.append({**rule, "consecutive_rows": consecutive})
    return result


def audit_event_log(path: Path, source_cycle: int, final_cycle: int,
                    thresholds: Any,
                    absolute_thresholds: dict[str, Any] | None = None) -> dict[str, Any]:
    rules = _normalize_event_thresholds(thresholds)
    if absolute_thresholds is None and isinstance(thresholds, dict):
        absolute_thresholds = thresholds
    _require(isinstance(absolute_thresholds, dict),
             "absolute event thresholds are required")
    _require(absolute_thresholds.get("hard_equal_zero") ==
             ["eos_fail", "eos_vceil"],
             "absolute event hard-zero schema changed")
    c2p_exclusive_max = _integer(
        absolute_thresholds.get("c2p_iterations_exclusive_max"),
        "event_absolute_thresholds.c2p_iterations_exclusive_max")
    try:
        lines = _read_stable_bytes(path).decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise CheckFailure(f"{path}: event log is not UTF-8") from exc
    headers: list[tuple[str, ...]] = []
    rows: list[dict[str, int]] = []
    for line_number, line in enumerate(lines, start=1):
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "#" and len(fields) > 1 and fields[1] == "cycle":
            headers.append(tuple(fields[1:]))
            continue
        if fields[0].startswith("#"):
            continue
        _require(len(headers) == 1, f"{path}:{line_number}: data precedes unique header")
        _require(len(fields) == len(EVENT_COLUMNS),
                 f"{path}:{line_number}: event row has wrong column count")
        try:
            values = [int(field) for field in fields]
        except ValueError as exc:
            raise CheckFailure(f"{path}:{line_number}: event value is not an integer") from exc
        row = dict(zip(EVENT_COLUMNS, values))
        _require(all(value >= 0 for value in values),
                 f"{path}:{line_number}: negative event value")
        rows.append(row)
    _require(headers == [EVENT_COLUMNS],
             f"event log must have exactly the strict schema {EVENT_COLUMNS}")
    cycles = [row["cycle"] for row in rows]
    expected_cycles = list(range(source_cycle + 1, final_cycle + 1))
    _require(cycles == expected_cycles,
             "event cycles must be ordered, unique, and exactly source+1..endpoint")
    ratio_observations: list[dict[str, Any]] = []
    for row in rows:
        _require(row["eos_fail"] == 0 and row["eos_vceil"] == 0,
                 f"cycle {row['cycle']}: eos_fail and eos_vceil must be zero")
        if c2p_exclusive_max is not None:
            _require(row["c2p_it"] < c2p_exclusive_max,
                     f"cycle {row['cycle']}: c2p_it {row['c2p_it']} is not below "
                     f"{c2p_exclusive_max}")
        for rule in rules:
            denominator = row[rule["denominator"]]
            _require(denominator > 0,
                     f"cycle {row['cycle']}: {rule['denominator']} denominator is zero")
            ratio = row[rule["numerator"]] / denominator
            _require(math.isfinite(ratio) and ratio <= rule["max_ratio"],
                     f"cycle {row['cycle']}: {rule['numerator']}/{rule['denominator']} "
                     f"={ratio} exceeds {rule['max_ratio']}")
    for rule in rules:
        ratios = [row[rule["numerator"]] / row[rule["denominator"]]
                  for row in rows]
        maximum_index = max(range(len(ratios)), key=ratios.__getitem__)
        ratio_observations.append({
            **rule,
            "rows_checked": len(rows),
            "maximum_ratio": ratios[maximum_index],
            "maximum_cycle": rows[maximum_index]["cycle"],
        })
    return {
        "schema": list(EVENT_COLUMNS),
        "rows": len(rows),
        "cycle_min": cycles[0] if cycles else None,
        "cycle_max": cycles[-1] if cycles else None,
        "thresholds": rules,
        "c2p_iterations_exclusive_max": c2p_exclusive_max,
        "totals": {
            name: sum(row[name] for row in rows) for name in EVENT_COLUMNS[1:]
        },
        "hard_ratio_observations": ratio_observations,
        "_rows": rows,
    }


def audit_event_ratio_advisories(event_audit: dict[str, Any],
                                 yellow_thresholds: Any) -> dict[str, Any]:
    """Report sustained yellow ratios while leaving hard pass/fail unchanged."""

    rules = _normalize_yellow_event_thresholds(yellow_thresholds)
    rows = event_audit.get("_rows")
    _require(isinstance(rows, list) and rows,
             "event audit lacks rows needed for scientific advisories")
    audits: list[dict[str, Any]] = []
    for rule in rules:
        observations: list[tuple[int, float]] = []
        for row in rows:
            _require(isinstance(row, dict), "event advisory row must be an object")
            denominator = _integer(row.get(rule["denominator"]),
                                   f"cycle event {rule['denominator']}")
            numerator = _integer(row.get(rule["numerator"]),
                                 f"cycle event {rule['numerator']}")
            _require(denominator > 0,
                     f"event advisory denominator {rule['denominator']} is zero")
            observations.append((_integer(row.get("cycle"), "event cycle"),
                                 numerator / denominator))

        runs: list[list[tuple[int, float]]] = []
        active: list[tuple[int, float]] = []
        for observation in observations:
            if observation[1] > rule["max_ratio"]:
                active.append(observation)
            elif active:
                runs.append(active)
                active = []
        if active:
            runs.append(active)
        qualifying = [run for run in runs if len(run) >= rule["consecutive_rows"]]
        maximum_cycle, maximum_ratio = max(observations, key=lambda item: item[1])
        triggered_runs: list[dict[str, Any]] = []
        for run in qualifying:
            peak_cycle, peak_ratio = max(run, key=lambda item: item[1])
            triggered_runs.append({
                "cycle_start": run[0][0],
                "cycle_end": run[-1][0],
                "rows": len(run),
                "maximum_ratio": peak_ratio,
                "maximum_cycle": peak_cycle,
            })
        audits.append({
            "name": rule["name"],
            "severity": "yellow" if triggered_runs else "green",
            "numerator": rule["numerator"],
            "denominator": rule["denominator"],
            "yellow_ratio_exclusive_min": rule["max_ratio"],
            "required_consecutive_rows": rule["consecutive_rows"],
            "rows_checked": len(observations),
            "cycle_min": observations[0][0],
            "cycle_max": observations[-1][0],
            "maximum_ratio": maximum_ratio,
            "maximum_cycle": maximum_cycle,
            "maximum_consecutive_exceedances": max(
                (len(run) for run in runs), default=0),
            "triggered_runs": triggered_runs,
        })
    return {
        "severity": ("yellow" if any(row["severity"] == "yellow" for row in audits)
                     else "green"),
        "ratios": audits,
    }


def audit_floor_rate_trends(event_audit: dict[str, Any],
                            parent: dict[str, Any] | None) -> dict[str, Any]:
    """Record floor rates normalized by C2P calls; never turn a trend into failure."""

    floor_columns = ("eos_dfloor", "eos_efloor", "eos_tfloor")

    def summary(audit: dict[str, Any], label: str) -> dict[str, Any]:
        totals = audit.get("totals")
        _require(isinstance(totals, dict), f"{label} event audit lacks totals")
        calls = _integer(totals.get("c2p_calls"), f"{label} c2p_calls")
        _require(calls > 0, f"{label} c2p_calls must be positive")
        cycles = (_integer(audit.get("cycle_min"), f"{label} cycle_min"),
                  _integer(audit.get("cycle_max"), f"{label} cycle_max"))
        counts = {name: _integer(totals.get(name), f"{label} {name}")
                  for name in floor_columns}
        _require(all(value >= 0 for value in counts.values()),
                 f"{label} floor counts must be nonnegative")
        return {
            "cycle_min": cycles[0], "cycle_max": cycles[1],
            "c2p_calls": calls, "counts": counts,
            "rates_per_c2p_call": {
                name: counts[name] / calls for name in floor_columns
            },
        }

    current = summary(event_audit, "current")
    if parent is None:
        comparison: dict[str, Any] = {
            "status": "unavailable_anchor",
            "reason": "no parent segment pass",
        }
    else:
        parent_summary = summary(parent.get("event_log_audit", {}), "parent")
        rates: dict[str, Any] = {}
        for name in floor_columns:
            current_rate = current["rates_per_c2p_call"][name]
            parent_rate = parent_summary["rates_per_c2p_call"][name]
            rates[name] = {
                "current_rate": current_rate,
                "parent_rate": parent_rate,
                "absolute_change": current_rate - parent_rate,
                "current_over_parent": (
                    current_rate / parent_rate if parent_rate > 0.0 else None),
                "trend": ("increased" if current_rate > parent_rate else
                          "decreased" if current_rate < parent_rate else "unchanged"),
            }
        comparison = {
            "status": "available",
            "parent": parent_summary,
            "rates": rates,
        }
    return {
        "severity": "green",
        "classification": "trend_only_no_pass_fail_threshold",
        "current": current,
        "parent_comparison": comparison,
    }


def parse_gpu_csv(path: Path) -> list[dict[str, Any]]:
    try:
        text = _read_stable_bytes(path).decode("utf-8")
    except UnicodeDecodeError as exc:
        raise CheckFailure(f"{path}: GPU CSV is not UTF-8") from exc
    parsed = list(csv.reader(text.splitlines(), skipinitialspace=True))
    parsed = [[field.strip() for field in row] for row in parsed if row]
    _require(parsed, f"{path}: GPU CSV has no device rows")
    result: list[dict[str, Any]] = []
    for line_number, row in enumerate(parsed, start=1):
        _require(len(row) == len(GPU_HEADER),
                 f"{path}: GPU row {line_number} must have exactly "
                 f"{len(GPU_HEADER)} columns")
        try:
            (index, cuda_ordinal, uncorrected, corrected,
             memory_total, memory_used) = (
                int(row[0]), int(row[3]), int(row[4]), int(row[5]),
                int(row[6]), int(row[7])
            )
        except ValueError as exc:
            raise CheckFailure(f"{path}: GPU row {line_number} has a non-integer field") from exc
        _require(index >= 0 and cuda_ordinal >= 0 and uncorrected >= 0 and
                 corrected >= 0 and memory_total > 0 and
                 0 <= memory_used <= memory_total,
                 f"{path}: GPU row {line_number} contains a negative value")
        _require(row[1].startswith("GPU-") and
                 PCI_BUS_ID_RE.fullmatch(row[2]) is not None,
                 f"{path}: GPU row {line_number} has an invalid UUID/PCI identity")
        result.append({
            "index": index,
            "uuid": row[1],
            "pci_bus_id": row[2].upper(),
            "cuda_ordinal": cuda_ordinal,
            "uncorrected_ecc": uncorrected,
            "corrected_ecc": corrected,
            "memory_total_mib": memory_total,
            "memory_used_mib": memory_used,
        })
    _require(len({row["index"] for row in result}) == len(result),
             f"{path}: duplicate GPU index")
    _require(len({row["uuid"] for row in result}) == len(result),
             f"{path}: duplicate GPU UUID")
    _require(len({row["pci_bus_id"] for row in result}) == len(result),
             f"{path}: duplicate GPU PCI bus identity")
    _require(len({row["cuda_ordinal"] for row in result}) == len(result),
             f"{path}: duplicate CUDA ordinal")
    _require([row["cuda_ordinal"] for row in result] == list(range(len(result))),
             f"{path}: rows are not in contiguous PCI-ordered CUDA ordinal order")
    _require([_pci_bus_sort_key(row["pci_bus_id"]) for row in result] == sorted(
        _pci_bus_sort_key(row["pci_bus_id"]) for row in result),
        f"{path}: GPU rows are not sorted by PCI bus identity")
    return result


def audit_gpus(before_path: Path, after_path: Path, ranks: int,
               exit_memory_mib_max: float,
               ecc_policy: dict[str, Any] | None = None) -> dict[str, Any]:
    before = parse_gpu_csv(before_path)
    after = parse_gpu_csv(after_path)
    _require(len(before) == ranks and len(after) == ranks,
             f"GPU snapshots must each contain exactly {ranks} devices")
    _require([row["cuda_ordinal"] for row in before] == list(range(ranks)) and
             [row["cuda_ordinal"] for row in after] == list(range(ranks)),
             "CUDA ordinals must be contiguous 0..ranks-1")
    _require([(row["uuid"], row["pci_bus_id"], row["cuda_ordinal"],
               row["memory_total_mib"])
              for row in before] ==
             [(row["uuid"], row["pci_bus_id"], row["cuda_ordinal"],
               row["memory_total_mib"])
              for row in after],
             "GPU UUID/PCI/CUDA/total-memory identities changed between "
             "preflight and exit")
    if ecc_policy is None:
        ecc_policy = {
            "uncorrected_before_max": 0, "corrected_before_max": 0,
            "uncorrected_after_max": 0, "corrected_after_max": 0,
        }
    expected_ecc_policy = {
        "uncorrected_before_max": 0, "corrected_before_max": 0,
        "uncorrected_after_max": 0, "corrected_after_max": 0,
    }
    _require(ecc_policy == expected_ecc_policy,
             "GPU ECC policy must require all four snapshots to equal zero")
    _require(all(row["uncorrected_ecc"] <= ecc_policy["uncorrected_before_max"] and
                 row["corrected_ecc"] <= ecc_policy["corrected_before_max"]
                 for row in before), "preflight GPU ECC count exceeds policy")
    _require(all(row["uncorrected_ecc"] <= ecc_policy["uncorrected_after_max"] and
                 row["corrected_ecc"] <= ecc_policy["corrected_after_max"]
                 for row in after), "exit GPU ECC count exceeds policy")
    _require(all(row["memory_used_mib"] <= exit_memory_mib_max for row in after),
             "exit GPU memory exceeds the planned quiescent threshold")
    return {"before": before, "after": after, "ecc_policy": ecc_policy,
            "exit_memory_mib_max": exit_memory_mib_max}


def audit_live_gpu_quiescence(plan: dict[str, Any], before_path: Path) \
        -> dict[str, Any]:
    """Query the plan-bound NVIDIA tool and prove identities/ECC/apps are quiet."""

    planned = plan["tools"]["nvidia_smi"]
    live_tool = _verify_planned_file(planned, "plan-bound nvidia-smi")
    executable = live_tool["path"]
    _require(os.access(executable, os.X_OK),
             "plan-bound nvidia-smi is not executable")
    values = plan["launch_contract"]["environment"]["values"]
    environment = {name: values[name] for name in ("HOME", "LANG", "LC_ALL")}

    def query(argument: str) -> list[list[str]]:
        try:
            result = subprocess.run(
                [executable, argument, "--format=csv,noheader,nounits"],
                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                text=True, env=environment)
        except (OSError, subprocess.CalledProcessError) as exc:
            detail = getattr(exc, "stderr", "") or str(exc)
            raise CheckFailure(
                f"cannot query live GPU quiescence: {detail.strip()}") from exc
        return [[field.strip() for field in line.split(",")]
                for line in result.stdout.splitlines() if line.strip()]

    gpu_rows = query(
        "--query-gpu=index,uuid,pci.bus_id,"
        "ecc.errors.uncorrected.volatile.total,"
        "ecc.errors.corrected.volatile.total,memory.total,memory.used")
    parsed: list[dict[str, Any]] = []
    for row in gpu_rows:
        _require(len(row) == 7, "live GPU query returned a malformed row")
        try:
            index, uncorr, corr, total, used = map(
                int, (row[0], row[3], row[4], row[5], row[6]))
        except ValueError as exc:
            raise CheckFailure("live GPU query returned a non-integer field") from exc
        pci = row[2].upper()
        _require(index >= 0 and row[1].startswith("GPU-") and
                 PCI_BUS_ID_RE.fullmatch(pci) is not None and
                 uncorr == 0 and corr == 0 and total > 0 and
                 0 <= used <= total,
                 "live GPU identity/ECC/memory row is invalid")
        parsed.append({
            "index": index, "uuid": row[1], "pci_bus_id": pci,
            "uncorrected_ecc": uncorr, "corrected_ecc": corr,
            "memory_total_mib": total, "memory_used_mib": used,
        })
    parsed.sort(key=lambda row: _pci_bus_sort_key(row["pci_bus_id"]))
    for ordinal, row in enumerate(parsed):
        row["cuda_ordinal"] = ordinal
    before = parse_gpu_csv(before_path)
    _require(len(parsed) == plan["policy"]["ranks"] and
             [(row["index"], row["uuid"], row["pci_bus_id"],
               row["cuda_ordinal"], row["memory_total_mib"])
              for row in parsed] ==
             [(row["index"], row["uuid"], row["pci_bus_id"],
               row["cuda_ordinal"], row["memory_total_mib"])
              for row in before],
             "live GPU inventory differs from launch preflight")
    _require(all(row["memory_used_mib"] <=
                 plan["policy"]["gpu_exit_memory_mib_max"] for row in parsed),
             "live GPU memory exceeds the quiescent threshold")
    applications = query("--query-compute-apps=pid,gpu_uuid")
    _require(not applications,
             "live GPU compute contexts are not empty after segment closure")
    return {
        "nvidia_smi": live_tool,
        "inventory": parsed,
        "launch_inventory": before,
        "compute_contexts": [],
        "compute_contexts_empty": True,
        "ecc_all_zero": True,
    }


def _pending_events_zero(metadata: RestartMetadata, label: str) -> None:
    counters = metadata.event_counters
    _require(counters is not None, f"{label} restart lacks pending event counters")
    values = counters.as_dict()
    _require(all(value == 0 for value in values.values()),
             f"{label} restart has nonzero pending event counters: {values}")


def audit_restart_contract(metadata: RestartMetadata, label: str,
                           plan: dict[str, Any],
                           expected: dict[str, Any]) -> dict[str, Any]:
    """Re-prove the strict restart and load-balance contract at either endpoint."""

    policy = plan["policy"]
    restart_policy = policy.get("restart_contract")
    capacity_policy = policy.get("capacity")
    _require(isinstance(restart_policy, dict),
             "plan restart_contract policy is missing")
    _require(isinstance(capacity_policy, dict),
             "plan capacity policy is missing")
    required_real_bytes = _integer(
        restart_policy.get("real_bytes"), "restart_contract.real_bytes")
    _require(required_real_bytes == 8 and metadata.real_bytes == required_real_bytes,
             f"{label} restart must use Real8")
    _require(metadata.last_dt == expected["root_dt"],
             f"{label} restart last_dt {metadata.last_dt!r} differs from fixed "
             f"root dt {expected['root_dt']!r}")
    time_parameters = metadata.parameters.get("time", {})
    _require(restart_policy.get("subcycling") == "level" and
             time_parameters.get("subcycling", "").strip() == "level",
             f"{label} restart does not use strict level subcycling")
    for key in ("allow_legacy_subcycling_restart",
                "allow_legacy_ghost_event_counters"):
        value = time_parameters.get(key, "").strip().lower()
        _require(restart_policy.get(key) is False and value in ("false", "0"),
                 f"{label} restart permits legacy mode {key}")

    required_amr = _integer(restart_policy.get("amr_counter_version"),
                            "restart_contract.amr_counter_version")
    required_event = _integer(restart_policy.get("event_counter_version"),
                              "restart_contract.event_counter_version")
    required_cache = _integer(
        restart_policy.get("restart_cache_contract_version"),
        "restart_contract.restart_cache_contract_version")
    _require(required_amr == 1 and metadata.amr_cycle_counters is not None,
             f"{label} restart lacks the required v1 AMR counters")
    _require(len(metadata.amr_cycle_counters) == metadata.nmb_total and
             all(isinstance(value, int) and value >= 0
                 for value in metadata.amr_cycle_counters),
             f"{label} restart has invalid or negative AMR counters")
    _require(metadata.event_counter_version == required_event,
             f"{label} restart event-counter version is not {required_event}")
    _require(metadata.restart_cache_contract_version == required_cache,
             f"{label} restart cache-contract version is not {required_cache}")
    _pending_events_zero(metadata, label)

    expected_costs = tuple(
        math.ldexp(1.0, location.level - metadata.root_level)
        for location in metadata.locations
    )
    _require(restart_policy.get("level_subcycling_costs_exact") is True and
             metadata.costs == expected_costs,
             f"{label} restart MeshBlock costs are not exact level costs")
    try:
        serialized_capacity = metadata.parameters[
            "mesh_refinement"]["max_nmb_per_rank"]
        capacity = int(serialized_capacity)
    except (KeyError, ValueError) as exc:
        raise CheckFailure(
            f"{label} restart has invalid max_nmb_per_rank") from exc
    _require(serialized_capacity == str(capacity),
             f"{label} restart max_nmb_per_rank is not canonical decimal")
    capacity_transition = _validate_capacity_transition(plan)
    is_source = metadata.cycle == expected["source_cycle"]
    required_capacity = capacity_transition[
        "source_max_nmb_per_rank" if is_source else
        "target_max_nmb_per_rank"]
    _require(capacity == required_capacity,
             f"{label} restart max_nmb_per_rank {capacity} differs from exact "
             f"planned {'source' if is_source else 'target'} capacity "
             f"{required_capacity}")
    ranks = _integer(policy.get("ranks"), "policy.ranks")
    hard_headroom = _integer(
        capacity_policy.get("minimum_headroom_blocks_hard"),
        "capacity.minimum_headroom_blocks_hard")
    _require(capacity > 0 and hard_headroom >= 0,
             "capacity and hard headroom must be nonnegative")
    partition_capacity = capacity_transition["target_max_nmb_per_rank"]
    try:
        serialized_partitions = load_balance(metadata.costs, ranks, capacity)
        target_partitions = load_balance(
            metadata.costs, ranks, partition_capacity)
    except ValueError as exc:
        raise CheckFailure(f"{label} restart cannot be partitioned: {exc}") from exc
    serialized_headroom = min(
        capacity - row.blocks for row in serialized_partitions)
    target_headroom = min(
        partition_capacity - row.blocks for row in target_partitions)
    _require(target_headroom >= hard_headroom,
             f"{label} restart target partition headroom {target_headroom} is below hard "
             f"minimum {hard_headroom}")
    return {
        "real_bytes": metadata.real_bytes,
        "last_dt": metadata.last_dt,
        "subcycling": "level",
        "amr_counter_version": required_amr,
        "event_counter_version": metadata.event_counter_version,
        "restart_cache_contract_version": metadata.restart_cache_contract_version,
        "pending_event_counters_all_zero": True,
        "level_subcycling_costs_exact": True,
        "ranks": ranks,
        "max_nmb_per_rank": capacity,
        "serialized_capacity_headroom": serialized_headroom,
        "target_partition_capacity": partition_capacity,
        "minimum_capacity_headroom": target_headroom,
        "hard_capacity_headroom": hard_headroom,
        "capacity_transition_role": "source" if is_source else "target",
        "capacity_transition": capacity_transition,
    }


def _flatten_parameters(parameters: dict[str, dict[str, str]]) -> dict[str, str]:
    return {
        f"{block}/{key}": value
        for block, values in parameters.items()
        for key, value in values.items()
    }


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def compare_parameters(source: RestartMetadata, endpoint: RestartMetadata,
                       mutable_patterns: list[str],
                       exact_rebindings: dict[str, dict[str, Any]] | None = None
                       ) -> dict[str, Any]:
    source_flat = _flatten_parameters(source.parameters)
    endpoint_flat = _flatten_parameters(endpoint.parameters)
    exact_rebindings = exact_rebindings or {}
    differences: list[dict[str, Any]] = []
    mutable_differences: list[str] = []
    observed_exact_rebindings: list[str] = []
    for key in sorted(source_flat.keys() | endpoint_flat.keys()):
        left, right = source_flat.get(key), endpoint_flat.get(key)
        if left == right:
            if key in exact_rebindings:
                rule = exact_rebindings[key]
                _require(left in rule.get("source", []) and
                         right == rule.get("endpoint"),
                         f"exact parameter rebinding {key} has invalid unchanged "
                         "source/endpoint values")
            continue
        if key in exact_rebindings:
            rule = exact_rebindings[key]
            _require(left in rule.get("source", []) and
                     right == rule.get("endpoint"),
                     f"exact parameter rebinding {key} differs from its bound "
                     "source/endpoint values")
            observed_exact_rebindings.append(key)
            continue
        if any(fnmatch.fnmatchcase(key, pattern) for pattern in mutable_patterns):
            mutable_differences.append(key)
        else:
            differences.append({"parameter": key, "source": left, "endpoint": right})
    _require(not differences, f"immutable restart parameters changed: {differences[:8]}")
    for key, rule in exact_rebindings.items():
        _require(source_flat.get(key) in rule.get("source", []) and
                 endpoint_flat.get(key) == rule.get("endpoint"),
                 f"required exact parameter rebinding {key} is absent or invalid")
    return {
        "mutable_patterns": mutable_patterns,
        "observed_mutable_differences": mutable_differences,
        "exact_rebindings": exact_rebindings,
        "observed_exact_rebindings": observed_exact_rebindings,
    }


def _safe_output_path(state_dir: Path, relative: str) -> Path:
    root = state_dir.resolve(strict=True)
    candidate = state_dir.joinpath(relative)
    current = candidate.absolute()
    while current != root.parent:
        _require(not current.is_symlink(),
                 f"output path traverses a symlink: {current}")
        if current == root:
            break
        current = current.parent
    resolved = candidate.resolve(strict=True)
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise CheckFailure(f"output escapes state directory: {candidate}") from exc
    return resolved


def _float32(value: float) -> float:
    try:
        return struct.unpack("=f", struct.pack("=f", value))[0]
    except (OverflowError, struct.error) as exc:
        raise CheckFailure(f"value cannot be represented as Athena output float32: {value}") \
            from exc


def _replay_output_plan(
        outputs: list[dict[str, Any]],
        source_parameters: dict[str, dict[str, str]],
        expected: dict[str, Any],
        ) -> tuple[dict[str, list[dict[str, Any]]],
                   dict[str, dict[str, Any]]]:
    """Replay Execute/Finalize from immutable source parameters alone."""

    states: dict[str, dict[str, Any]] = {}
    for output in outputs:
        block = output["block"]
        source_block = source_parameters.get(block)
        _require(source_block is not None, f"source restart parameter dump lacks {block}")
        try:
            file_number = int(source_block.get("file_number", "0"))
            last_time = float(source_block.get("last_time", "-1"))
            last_write_cycle = int(source_block.get("last_write_cycle", "-1"))
        except ValueError as exc:
            raise CheckFailure(f"{block}: invalid serialized cadence state") from exc
        _require(file_number >= 0 and math.isfinite(last_time),
                 f"{block}: invalid source cadence state")
        cadence_mode = output.get("cadence_mode")
        cadence_value = _number(output.get("cadence"), f"{block}.cadence")
        _require(cadence_mode in ("dt", "dcycle") and cadence_value > 0.0,
                 f"{block}: unsupported or nonpositive cadence")
        if cadence_mode == "dcycle":
            _require(cadence_value.is_integer(),
                     f"{block}: dcycle cadence must be an integer")
            dt, dcycle = 0.0, int(cadence_value)
        else:
            dt, dcycle = cadence_value, 0
        states[block] = {
            "file_number": file_number,
            "last_time": last_time,
            "last_write_cycle": last_write_cycle,
            "dt": dt,
            "dcycle": dcycle,
            "wrote_this_run": False,
            "writes": [],
        }

    def due(state: dict[str, Any], cycle: int, output_time: float,
            enforce_time_limit: bool) -> bool:
        time_due = (state["dt"] > 0.0 and
                    _float32(output_time) >=
                    _float32(state["last_time"] + state["dt"]) and
                    (not enforce_time_limit or
                     _float32(output_time) < _float32(expected["tlim"])))
        cycle_due = (state["dcycle"] > 0 and
                     cycle % state["dcycle"] == 0)
        return time_due or cycle_due

    def write(output: dict[str, Any], state: dict[str, Any], cycle: int,
              output_time: float, kind: str, advance_cadence: bool) -> None:
        row: dict[str, Any] = {
            "cycle": cycle, "time": output_time, "kind": kind,
        }
        if output["numbered"]:
            row["file_number"] = state["file_number"]
            state["file_number"] += 1
        state["writes"].append(row)
        if advance_cadence:
            if state["last_time"] < 0.0:
                state["last_time"] = output_time
            else:
                state["last_time"] += state["dt"]
        state["last_write_cycle"] = cycle
        state["wrote_this_run"] = True

    current_time = expected["source_time"]
    for cycle in range(expected["source_cycle"] + 1, expected["final_cycle"] + 1):
        current_time += expected["root_dt"]
        for output in outputs:
            state = states[output["block"]]
            if due(state, cycle, current_time, True):
                write(output, state, cycle, current_time, "scheduled", True)

    # Outputs stores all non-restart objects by push-front and the sole restart at the
    # tail.  Therefore Finalize visits non-restarts in reverse block order, then restart.
    non_restarts = sorted(
        (output for output in outputs if output["file_type"] != "rst"),
        key=lambda output: output["index"], reverse=True)
    restarts = [output for output in outputs
                if output["file_type"] == "rst"]
    _require(len(restarts) == 1,
             "strict segment qualification requires exactly one restart output")
    final_parameter_state_changed = False
    for output in [*non_restarts, *restarts]:
        state = states[output["block"]]
        is_restart = output["file_type"] == "rst"
        wrote_current = (state["wrote_this_run"] and
                         state["last_write_cycle"] == expected["final_cycle"])
        if wrote_current and (not is_restart or not final_parameter_state_changed):
            continue
        advance = (state["last_write_cycle"] != expected["final_cycle"] and
                   due(state, expected["final_cycle"], expected["final_time"], False))
        kind = "restart_final_rewrite" if is_restart and wrote_current else "forced_final"
        write(output, state, expected["final_cycle"], expected["final_time"],
              kind, advance)
        if not is_restart:
            final_parameter_state_changed = True

    schedules = {block: state["writes"] for block, state in states.items()}
    endpoint_states = {
        block: {
            "file_number": state["file_number"],
            "last_time": state["last_time"],
            "last_write_cycle": state["last_write_cycle"],
        }
        for block, state in states.items()
    }
    return schedules, endpoint_states


def _replay_all_output_writes(plan: dict[str, Any], source: RestartMetadata,
                              endpoint: RestartMetadata,
                              expected: dict[str, Any]
                              ) -> dict[str, list[dict[str, Any]]]:
    """Replay Execute/Finalize and re-prove the endpoint restart state."""

    schedules, endpoint_states = _replay_output_plan(
        plan["outputs"], source.parameters, expected)
    for block, replayed in endpoint_states.items():
        endpoint_block = endpoint.parameters.get(block)
        _require(endpoint_block is not None,
                 f"endpoint restart parameter dump lacks {block}")
        try:
            endpoint_state = (
                int(endpoint_block.get("file_number", "0")),
                float(endpoint_block["last_time"]),
                int(endpoint_block["last_write_cycle"]),
            )
        except (KeyError, ValueError) as exc:
            raise CheckFailure(f"{block}: endpoint lacks valid cadence state") from exc
        replayed_state = (
            replayed["file_number"], replayed["last_time"],
            replayed["last_write_cycle"])
        _require(endpoint_state == replayed_state,
                 f"{block}: endpoint cadence state does not match independent "
                 f"Execute/Finalize replay: {endpoint_state} != {replayed_state}")
    return schedules


def _assert_exact_numbered_output_set(state_dir: Path, template: str,
                                      writes: list[dict[str, Any]],
                                      block: str) -> None:
    """Reject every missing, extra, or noncanonically padded numbered artifact."""

    marker = "{file_number:05d}"
    prefix_text, suffix_text = template.split(marker)
    parent_relative = str(Path(prefix_text).parent)
    parent = _safe_output_path(
        state_dir, parent_relative if parent_relative != "." else "")
    name_prefix = Path(prefix_text).name
    actual_numbered: set[str] = set()
    for child in parent.iterdir():
        name = child.name
        if not name.startswith(name_prefix) or not name.endswith(suffix_text):
            continue
        middle_end = len(name) - len(suffix_text) if suffix_text else len(name)
        middle = name[len(name_prefix):middle_end]
        if re.fullmatch(r"\d+", middle):
            actual_numbered.add(name)
    expected_names = {
        Path(template.format(file_number=write["file_number"])).name
        for write in writes
    }
    _require(actual_numbered == expected_names,
             f"{block}: numbered output directory set differs from cadence replay: "
             f"actual={sorted(actual_numbered)}, expected={sorted(expected_names)}")


def expected_output_paths(plan: dict[str, Any], source: RestartMetadata,
                          endpoint: RestartMetadata, state_dir: Path,
                          expected: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    if expected is None:
        expected = validate_plan(plan)
    inventory: list[dict[str, Any]] = []
    seen: set[Path] = set()
    schedules = _replay_all_output_writes(plan, source, endpoint, expected)
    restart_cadence_transition = _validate_restart_cadence_transition(
        plan, float(expected["root_dt"]))
    for output in plan["outputs"]:
        block = output["block"]
        source_block = source.parameters.get(block)
        endpoint_block = endpoint.parameters.get(block)
        _require(source_block is not None and endpoint_block is not None,
                 f"restart parameter dump lacks {block}")
        expected_runtime_block = dict(source_block)
        if block == "output3":
            expected_runtime_block["dt"] = repr(expected["root_dt"])
        if (block == restart_cadence_transition["block"] and
                restart_cadence_transition["kind"] == "tighten_v1"):
            expected_runtime_block["dt"] = repr(
                restart_cadence_transition["target_dt"])
        _require(output["parameters"] == expected_runtime_block,
                 f"plan snapshot for {block} differs from the source restart plus "
                 "the exact canonical runtime override")
        snapshot_sha = output.get("parameters_sha256")
        _require(isinstance(snapshot_sha, str) and
                 snapshot_sha == _canonical_sha256(expected_runtime_block),
                 f"plan snapshot hash for {block} is invalid")
        template = output["relative_path_template"]
        writes = schedules[block]
        if output["numbered"]:
            _require(writes, f"{block}: numbered output produced no expected files")
            rows = [(write, template.format(file_number=write["file_number"]))
                    for write in writes]
            _assert_exact_numbered_output_set(state_dir, template, writes, block)
        else:
            rows = [(None, relative)
                    for relative in output["required_unnumbered_paths"]]
        for write, relative in rows:
            path = _safe_output_path(state_dir, relative)
            _require(path not in seen, f"duplicate planned output path: {path}")
            seen.add(path)
            inventory.append({
                "block": block,
                "file_type": output.get("file_type"),
                "file_number": (write["file_number"]
                                if isinstance(write, dict) else None),
                "path": path,
                "inspect_binary": bool(output.get("inspect_binary", False)),
                "parameters": output["parameters"],
                "expected_binary_variables": output.get(
                    "expected_binary_variables"),
                "expected_write": write,
            })
    return inventory


def _assert_exact_state_tree(state_dir: Path,
                             expected_paths: Iterable[Path]) -> dict[str, Any]:
    """Reject every unplanned file and every non-regular node below state_dir."""

    root_absolute = state_dir.absolute()
    root_stat = root_absolute.lstat()
    _require(stat.S_ISDIR(root_stat.st_mode) and
             not stat.S_ISLNK(root_stat.st_mode),
             "state directory must be a real, non-symlink directory")
    root = root_absolute.resolve(strict=True)
    expected = {path.resolve(strict=True) for path in expected_paths}
    expected_directories: set[Path] = {root}
    for path in expected:
        parent = path.parent
        while True:
            expected_directories.add(parent)
            if parent == root:
                break
            try:
                parent.relative_to(root)
            except ValueError as exc:
                raise CheckFailure(
                    f"expected state output escapes state directory: {path}") from exc
            parent = parent.parent
    actual: set[Path] = set()
    actual_directories: set[Path] = set()
    pending = [root]
    while pending:
        directory = pending.pop()
        actual_directories.add(directory.resolve(strict=True))
        try:
            entries = list(os.scandir(directory))
        except OSError as exc:
            raise CheckFailure(f"cannot inventory state directory {directory}: {exc}") \
                from exc
        for entry in entries:
            path = Path(entry.path)
            try:
                node = entry.stat(follow_symlinks=False)
            except OSError as exc:
                raise CheckFailure(f"cannot stat state node {path}: {exc}") from exc
            _require(not stat.S_ISLNK(node.st_mode),
                     f"state directory contains a symlink: {path}")
            if stat.S_ISDIR(node.st_mode):
                pending.append(path)
            else:
                _require(stat.S_ISREG(node.st_mode),
                         f"state directory contains a non-regular node: {path}")
                actual.add(path.resolve(strict=True))
    _require(actual == expected,
             "state directory regular-file set differs from the exact replayed "
             f"allowlist: missing={sorted(str(path) for path in expected - actual)}, "
             f"extra={sorted(str(path) for path in actual - expected)}")
    _require(actual_directories == expected_directories,
             "state directory set differs from exact replayed allowlist: "
             f"missing={sorted(str(path) for path in expected_directories - actual_directories)}, "
             f"extra={sorted(str(path) for path in actual_directories - expected_directories)}")
    return {
        "regular_files": len(actual),
        "directories": len(actual_directories),
        "exact_replayed_allowlist": True,
        "nonregular_nodes": 0,
    }


def audit_history(path: Path, source_cycle: int, final_cycle: int,
                  source_time: float, root_dt: float, endpoint_time: float,
                  endpoint_ulps: int) -> dict[str, Any]:
    """Validate a fresh Athena history segment and return column endpoints."""

    try:
        lines = _read_stable_bytes(path).decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise CheckFailure(f"{path}: history file is not UTF-8") from exc
    headers: list[list[str]] = []
    rows: list[list[float]] = []
    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            names = re.findall(r"\[\d+\]=([^\s]+)", stripped)
            if names:
                headers.append(names)
            continue
        _require(len(headers) == 1,
                 f"{path}:{line_number}: history data precedes a unique header")
        fields = stripped.split()
        _require(len(fields) == len(headers[0]),
                 f"{path}:{line_number}: history row has wrong column count")
        try:
            values = [float(value) for value in fields]
        except ValueError as exc:
            raise CheckFailure(f"{path}:{line_number}: malformed history value") from exc
        _require(all(math.isfinite(value) for value in values),
                 f"{path}:{line_number}: non-finite history value")
        rows.append(values)
    _require(len(headers) == 1 and len(headers[0]) == len(set(headers[0])),
             f"{path}: history must contain one unique-column header")
    _require(rows, f"{path}: history contains no data rows")
    columns = {name: [row[index] for row in rows]
               for index, name in enumerate(headers[0])}
    row_dicts = [dict(zip(headers[0], row)) for row in rows]
    _require("time" in columns, f"{path}: history has no time column")
    _require("dt" in columns, f"{path}: history has no dt column")
    times = columns["time"]
    root_steps = final_cycle - source_cycle
    _require(root_steps > 0 and len(times) == root_steps,
             f"{path}: history must contain exactly one row per completed root step")
    expected_times: list[float] = []
    current_time = source_time
    for _ in range(root_steps):
        current_time += root_dt
        expected_times.append(current_time)
    _require(all(_close_in_ulps(actual, planned, endpoint_ulps)
                 for actual, planned in zip(times, expected_times)),
             f"{path}: history times do not follow every sequential fixed root step")
    _require(all(_close_in_ulps(value, root_dt, endpoint_ulps)
                 for value in columns["dt"]),
             f"{path}: history contains a non-fixed root dt")
    _require(_close_in_ulps(times[-1], endpoint_time, endpoint_ulps),
             f"{path}: history does not reach the planned endpoint")
    baryon_step_loss: dict[str, Any] | None = None
    if "baryon_m" in columns:
        masses = columns["baryon_m"]
        _require(all(value > 0.0 for value in masses),
                 f"{path}: baryon_m must remain strictly positive")
        losses = [max(0.0, (left - right) / left)
                  for left, right in zip(masses, masses[1:])]
        maximum = max(losses, default=0.0)
        maximum_index = losses.index(maximum) if losses else None
        baryon_step_loss = {
            "intervals_checked": len(losses),
            "maximum_fractional_loss": maximum,
            "from_time": (times[maximum_index]
                          if maximum_index is not None else None),
            "to_time": (times[maximum_index + 1]
                        if maximum_index is not None else None),
        }
    return {
        "columns": headers[0],
        "rows": len(rows),
        "time_min": times[0],
        "time_max": times[-1],
        "implicit_cycle_min": source_cycle + 1,
        "implicit_cycle_max": final_cycle,
        "all_root_cycles_present": True,
        "fixed_root_dt": root_dt,
        "endpoint_time_matches": _close_in_ulps(
            times[-1], endpoint_time, endpoint_ulps),
        "column_first": {name: values[0] for name, values in columns.items()},
        "column_last": {name: values[-1] for name, values in columns.items()},
        "all_values_finite": True,
        "baryon_step_loss": baryon_step_loss,
        "_row_values": row_dicts,
    }


def audit_baryon_mass(history: dict[str, Any], hard_per_root_step: Any,
                      root_steps: int, source_mass: Any) -> dict[str, Any]:
    """Apply both adjacent-step and cumulative baryon-loss hard limits."""

    _require("baryon_m" in history.get("columns", []),
             "history does not contain baryon_m")
    source = _number(source_mass, "source baryon mass")
    first = _number(history["column_first"].get("baryon_m"),
                    "first evolved baryon mass")
    last = _number(history["column_last"].get("baryon_m"),
                   "last baryon mass")
    _require(source > 0.0 and first > 0.0 and last > 0.0,
             "baryon-mass history contains a nonphysical endpoint")
    hard = _number(hard_per_root_step,
                   "baryon_mass_fractional_loss.hard_per_root_step")
    _require(hard >= 0.0, "baryon-mass loss threshold must be nonnegative")
    step_loss = history.get("baryon_step_loss")
    _require(isinstance(step_loss, dict),
             "baryon history did not produce adjacent-step loss diagnostics")
    maximum_step_loss = _number(
        step_loss.get("maximum_fractional_loss"),
        "observed baryon-mass maximum adjacent-step loss")
    source_first_loss = max(0.0, (source - first) / source)
    maximum_step_loss = max(maximum_step_loss, source_first_loss)
    _require(maximum_step_loss <= hard,
             f"baryon-mass adjacent-step fractional loss {maximum_step_loss} "
             f"exceeds {hard}")
    total_loss = max(0.0, (source - last) / source)
    total_limit = hard * root_steps
    _require(math.isfinite(total_limit) and total_loss <= total_limit,
             f"baryon-mass fractional loss {total_loss} exceeds segment "
             f"limit {total_limit}")
    return {
        "source": source,
        "first": first,
        "last": last,
        "fractional_loss": total_loss,
        "hard_fractional_loss": total_limit,
        "adjacent_step": {
            **step_loss,
            "intervals_checked": root_steps,
            "source_to_first_fractional_loss": source_first_loss,
            "maximum_fractional_loss": maximum_step_loss,
            "hard_fractional_loss": hard,
        },
    }


def audit_baryon_mass_advisory(history: dict[str, Any], source_mass: Any,
                               source_cycle: int, yellow_per_root_step: Any,
                               yellow_per_window: Any,
                               window_root_steps: Any) -> dict[str, Any]:
    """Compute per-step and sliding-window yellow baryon-loss diagnostics."""

    rows = history.get("_row_values")
    _require(isinstance(rows, list) and rows,
             "baryon history lacks rows needed for scientific advisories")
    source = _number(source_mass, "source baryon mass")
    _require(source > 0.0, "source baryon mass must be positive")
    step_threshold = _number(yellow_per_root_step,
                             "baryon yellow_per_root_step")
    window_threshold = _number(yellow_per_window, "baryon yellow_per_48M")
    window_steps = _integer(window_root_steps, "baryon rolling_window_root_steps")
    _require(step_threshold >= 0.0 and window_threshold >= 0.0 and window_steps > 0,
             "baryon yellow thresholds and window must be nonnegative/positive")
    root_dt = _number(history.get("fixed_root_dt"), "baryon fixed_root_dt")
    first_cycle = _integer(history.get("implicit_cycle_min"),
                           "baryon implicit_cycle_min")
    _require(first_cycle == source_cycle + 1,
             "baryon advisory source cycle differs from history")
    masses = [source]
    times = [_number(history.get("time_min"), "baryon time_min") - root_dt]
    for index, row in enumerate(rows):
        _require(isinstance(row, dict), f"baryon row {index} must be an object")
        mass = _number(row.get("baryon_m"), f"baryon row {index} mass")
        row_time = _number(row.get("time"), f"baryon row {index} time")
        _require(mass > 0.0, f"baryon row {index} mass must be positive")
        masses.append(mass)
        times.append(row_time)
    cycles = list(range(source_cycle, source_cycle + len(masses)))

    adjacent = [max(0.0, (left - right) / left)
                for left, right in zip(masses, masses[1:])]
    maximum_step_index = max(range(len(adjacent)), key=adjacent.__getitem__)
    maximum_step = adjacent[maximum_step_index]
    step_audit = {
        "severity": "yellow" if maximum_step > step_threshold else "green",
        "yellow_fractional_loss_exclusive_min": step_threshold,
        "intervals_checked": len(adjacent),
        "maximum_fractional_loss": maximum_step,
        "cycle_from": cycles[maximum_step_index],
        "cycle_to": cycles[maximum_step_index + 1],
        "time_from": times[maximum_step_index],
        "time_to": times[maximum_step_index + 1],
    }

    rolling: list[tuple[int, float]] = []
    for index in range(len(masses) - window_steps):
        rolling.append((index, max(
            0.0, (masses[index] - masses[index + window_steps]) / masses[index])))
    if rolling:
        maximum_window_index, maximum_window = max(rolling, key=lambda item: item[1])
        window_audit: dict[str, Any] = {
            "severity": ("yellow" if maximum_window > window_threshold else "green"),
            "status": "evaluated",
            "yellow_fractional_loss_exclusive_min": window_threshold,
            "window_root_steps": window_steps,
            "window_duration": root_dt * window_steps,
            "windows_checked": len(rolling),
            "maximum_fractional_loss": maximum_window,
            "cycle_from": cycles[maximum_window_index],
            "cycle_to": cycles[maximum_window_index + window_steps],
            "time_from": times[maximum_window_index],
            "time_to": times[maximum_window_index + window_steps],
        }
    else:
        window_audit = {
            "severity": "green",
            "status": "insufficient_segment_span",
            "yellow_fractional_loss_exclusive_min": window_threshold,
            "window_root_steps": window_steps,
            "window_duration": root_dt * window_steps,
            "windows_checked": 0,
            "maximum_fractional_loss": None,
            "cycle_from": None, "cycle_to": None,
            "time_from": None, "time_to": None,
        }
    return {
        "severity": ("yellow" if step_audit["severity"] == "yellow" or
                     window_audit["severity"] == "yellow" else "green"),
        "per_root_step": step_audit,
        "rolling_window": window_audit,
    }


def audit_source_baryon_evidence(source_qualification: dict[str, Any],
                                 source_time: float,
                                 parent: dict[str, Any] | None = None) -> float:
    """Independently recover source mass from the bound history bytes."""

    planned = source_qualification["source_baryon_mass"]
    planned_mass = _number(planned.get("value"), "source_baryon_mass.value")
    evidence = planned["evidence"]
    if source_qualification["mode"] == "parent_segment_pass":
        _require(parent is not None,
                 "parent source qualification lacks its checked pass report")
    path = Path(evidence["path"])
    try:
        lines = _read_stable_bytes(path).decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise CheckFailure(f"{path}: source baryon history is not UTF-8") from exc
    headers = [re.findall(r"\[\d+\]=([^\s]+)", line)
               for line in lines if line.lstrip().startswith("#")]
    headers = [header for header in headers if header]
    _require(len(headers) == 1 and "time" in headers[0] and
             "baryon_m" in headers[0] and
             len(headers[0]) == len(set(headers[0])),
             "source baryon evidence must have one time/baryon_m header")
    last: list[float] | None = None
    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        try:
            row = [float(token) for token in stripped.split()]
        except ValueError as exc:
            raise CheckFailure(
                f"{path}:{line_number}: malformed source baryon history row"
            ) from exc
        _require(len(row) == len(headers[0]) and
                 all(math.isfinite(value) for value in row),
                 f"{path}:{line_number}: invalid source baryon history row")
        last = row
    _require(last is not None, "source baryon evidence contains no data row")
    columns = {name: index for index, name in enumerate(headers[0])}
    _require(last[columns["time"]] == source_time,
             "source baryon evidence does not end at the source time")
    observed = last[columns["baryon_m"]]
    _require(math.isfinite(observed) and observed > 0.0 and observed == planned_mass,
             "planned source baryon mass differs from independently checked evidence")
    return observed


def classify_recovery_restart(path: Path) -> dict[str, Any]:
    """Independently distinguish a truncated write from a complete invalid one."""

    if not os.path.lexists(path):
        return {"classification": "absent"}
    binding = _hash_record(path)
    try:
        metadata = read_restart_metadata(path)
    except RestartTruncationError as exc:
        return {"classification": "incomplete_truncated", "binding": binding,
                "failure": str(exc)}
    except (OSError, RuntimeError, ValueError) as exc:
        return {"classification": "complete_invalid", "binding": binding,
                "failure": str(exc)}
    try:
        auditor_module = sys.modules[audit_restart.__module__]
        layout = auditor_module._derive_layout(metadata)
        expected_size = metadata.metadata_end + struct.calcsize("Q") + \
            metadata.nmb_total * layout.data_size
    except (AttributeError, OSError, RuntimeError, ValueError):
        return {"classification": "complete_invalid", "binding": binding,
                "failure": "restart layout cannot be derived"}
    if binding["size"] < expected_size:
        return {"classification": "incomplete_truncated", "binding": binding,
                "expected_size": expected_size}
    if binding["size"] > expected_size:
        return {"classification": "complete_invalid", "binding": binding,
                "expected_size": expected_size,
                "failure": "restart contains trailing bytes"}
    try:
        audited = audit_restart(path)
    except (OSError, RuntimeError, ValueError) as exc:
        return {"classification": "complete_invalid", "binding": binding,
                "failure": str(exc)}
    return {"classification": "complete", "binding": binding,
            "audit": audited}


def _audit_recovery_candidates(original_plan: dict[str, Any], record: dict[str, Any],
                               selected: dict[str, Any]) -> dict[str, Any]:
    """Rebuild the complete scheduled-restart candidate set from original plan."""

    validate_plan(original_plan)
    restarts = [row for row in original_plan["outputs"]
                if row.get("file_type") == "rst"]
    _require(len(restarts) == 1, "recovery original plan lacks one restart stream")
    restart = restarts[0]
    scheduled = [row for row in restart["expected_writes"]
                 if row.get("kind") == "scheduled"]
    _require(scheduled, "recovery original plan has no scheduled restart")
    expected_rows: dict[str, dict[str, Any]] = {}
    for write in scheduled:
        relative = restart["relative_path_template"].format(
            file_number=write["file_number"])
        _require(relative not in expected_rows,
                 "recovery original plan has duplicate scheduled restart path")
        expected_rows[relative] = write
    candidates = record.get("candidate_inventory")
    _require(isinstance(candidates, list) and len(candidates) == len(expected_rows),
             "recovery candidate inventory is not the complete scheduled set")
    actual_rows: dict[str, dict[str, Any]] = {}
    for row in candidates:
        _require(isinstance(row, dict) and
                 isinstance(row.get("relative_path"), str) and
                 row["relative_path"] not in actual_rows and
                 row["relative_path"] in expected_rows and
                 row.get("expected_write") == expected_rows[row["relative_path"]],
                 "recovery candidate key/write differs from original plan")
        actual_rows[row["relative_path"]] = row
    _require(set(actual_rows) == set(expected_rows),
             "recovery candidate inventory omits a scheduled restart")
    state_dir = Path(original_plan["launch_contract"]["state_dir"]).resolve(
        strict=True)
    live_rows: list[dict[str, Any]] = []
    for relative, write in expected_rows.items():
        rel = Path(relative)
        _require(not rel.is_absolute() and relative == rel.as_posix() and
                 all(part not in ("", ".", "..") for part in rel.parts),
                 "recovery candidate relative path is not canonical")
        candidate_path = state_dir / rel
        current = candidate_path
        while current != state_dir:
            _require(not current.is_symlink(),
                     f"recovery candidate path traverses symlink: {current}")
            parent = current.parent
            _require(parent != current,
                     "recovery candidate path does not descend from state root")
            current = parent
        classified = classify_recovery_restart(candidate_path)
        classification = classified["classification"]
        binding = classified.get("binding")
        reported = actual_rows[relative]
        _require(reported.get("classification") == classification,
                 f"recovery candidate live classification changed: {relative}")
        if binding is None:
            _require("binding" not in reported,
                     f"absent recovery candidate has a binding: {relative}")
        else:
            _require(_same_full_binding(reported.get("binding"), binding),
                     f"recovery candidate binding changed: {relative}")
        live_rows.append({"relative_path": relative, "expected_write": write,
                          "classification": classification,
                          "binding": binding})
    ordered = sorted(live_rows,
                     key=lambda row: row["expected_write"]["cycle"], reverse=True)
    highest_complete = next((row for row in ordered
                             if row["classification"] == "complete"), None)
    _require(highest_complete is not None,
             "recovery no longer has a complete scheduled restart")
    higher = ordered[:ordered.index(highest_complete)]
    _require(all(row["classification"] in ("absent", "incomplete_truncated")
                 for row in higher),
             "a later complete-invalid restart forbids recovery fallback")
    _require(selected.get("block") == restart["block"] and
             selected.get("relative_path") == highest_complete["relative_path"] and
             selected.get("expected_write") == highest_complete["expected_write"] and
             _same_full_binding(selected.get("binding"),
                                highest_complete["binding"]),
             "recovery did not select the live highest complete scheduled restart")
    return {"scheduled_candidates": len(live_rows),
            "selected_relative_path": highest_complete["relative_path"],
            "live_reclassification": live_rows}


def _replay_recovery_execute_prefix(
        plan: dict[str, Any], final_cycle: int
        ) -> tuple[dict[str, list[dict[str, Any]]],
                   dict[str, dict[str, Any]], float]:
    """Replay only AthenaK's Execute output loop through a scheduled restart."""

    expected = plan["expected"]
    _require(expected["source_cycle"] < final_cycle <= expected["final_cycle"],
             "recovery endpoint cycle is outside the original plan")
    states: dict[str, dict[str, Any]] = {}
    schedules: dict[str, list[dict[str, Any]]] = {}
    for output in plan["outputs"]:
        block = output["block"]
        parameters = output["parameters"]
        try:
            file_number = int(parameters.get("file_number", "0"))
            last_time = float(parameters.get("last_time", "-1"))
            last_write_cycle = int(parameters.get("last_write_cycle", "-1"))
        except ValueError as exc:
            raise CheckFailure(
                f"{block}: invalid recovery source cadence state") from exc
        cadence = _number(output.get("cadence"), f"{block}.cadence")
        if output.get("cadence_mode") == "dt":
            dt, dcycle = cadence, 0
        else:
            _require(cadence.is_integer(),
                     f"{block}: recovery dcycle is not integral")
            dt, dcycle = 0.0, int(cadence)
        states[block] = {
            "file_number": file_number, "last_time": last_time,
            "last_write_cycle": last_write_cycle, "dt": dt,
            "dcycle": dcycle,
        }
        schedules[block] = []

    current_time = expected["source_time"]
    for cycle in range(expected["source_cycle"] + 1, final_cycle + 1):
        current_time += expected["root_dt"]
        for output in plan["outputs"]:
            state = states[output["block"]]
            due = ((state["dt"] > 0.0 and
                    _float32(current_time) >=
                    _float32(state["last_time"] + state["dt"]) and
                    _float32(current_time) < _float32(expected["tlim"])) or
                   (state["dcycle"] > 0 and cycle % state["dcycle"] == 0))
            if not due:
                continue
            write: dict[str, Any] = {
                "cycle": cycle, "time": current_time, "kind": "scheduled",
            }
            if output["numbered"]:
                write["file_number"] = state["file_number"]
                state["file_number"] += 1
            schedules[output["block"]].append(write)
            state["last_time"] = (current_time if state["last_time"] < 0.0
                                  else state["last_time"] + state["dt"])
            state["last_write_cycle"] = cycle
    return schedules, states, current_time


def _audit_recovery_derived_prefixes(
        original_plan: dict[str, Any], record: dict[str, Any],
        parent: dict[str, Any]) -> dict[str, Any]:
    """Reopen every original/derived text and prove the exact byte split."""

    selected = record.get("selected_scheduled_restart")
    _require(isinstance(selected, dict) and
             isinstance(selected.get("expected_write"), dict),
             "recovery record lacks its selected scheduled restart")
    selected_write = selected["expected_write"]
    final_cycle = _integer(selected_write.get("cycle"), "recovery final cycle")
    schedules, _, final_time = _replay_recovery_execute_prefix(
        original_plan, final_cycle)
    source_metadata = read_restart_metadata(
        Path(original_plan["inputs"]["source_restart"]["path"]))
    endpoint_metadata = read_restart_metadata(
        Path(selected["binding"]["path"]))
    root_steps = final_cycle - original_plan["expected"]["source_cycle"]
    expected_relatives = set(original_plan.get("required_files", []))
    expected_numbered_rows: dict[str, dict[str, Any]] = {}
    output_for_required: dict[str, dict[str, Any]] = {}
    for output in original_plan["outputs"]:
        if output["numbered"]:
            for write in schedules[output["block"]]:
                relative = output["relative_path_template"].format(
                    file_number=write["file_number"])
                _require(relative not in expected_numbered_rows,
                         "recovery Execute replay maps duplicate numbered path")
                expected_numbered_rows[relative] = {
                    "block": output["block"], "file_type": output["file_type"],
                    "file_number": write["file_number"],
                    "expected_write": write,
                }
        else:
            for relative in output["required_unnumbered_paths"]:
                _require(relative not in output_for_required,
                         "original plan maps duplicate required output")
                output_for_required[relative] = output
    _require(set(output_for_required) == expected_relatives,
             "original required-files/output mapping is not exact")
    expected_logical = {
        "source_cycle": original_plan["expected"]["source_cycle"],
        "source_time": original_plan["expected"]["source_time"],
        "final_cycle": final_cycle, "final_time": final_time,
        "root_steps": root_steps, "execute_only": True,
        "expected_numbered_paths": sorted(expected_numbered_rows),
        "required_unnumbered_paths": sorted(expected_relatives),
        "all_expected_prefix_artifacts_present": True,
    }
    _require(record.get("logical_prefix") == expected_logical,
             "recovery logical prefix differs from independent Execute replay")
    _require(selected_write.get("time") == final_time and
             selected.get("block") in schedules and
             selected_write in schedules[selected["block"]],
             "recovery selected restart is not an Execute-scheduled write")

    rows = record.get("derived_text_prefixes")
    _require(isinstance(rows, list) and len(rows) == len(expected_relatives),
             "recovery derived-prefix inventory is not exact")
    by_relative: dict[str, dict[str, Any]] = {}
    original_state = Path(original_plan["launch_contract"]["state_dir"]).resolve(
        strict=True)
    recovery_state = Path(record["recovery_state_dir"]).resolve(strict=True)
    for row in rows:
        _require(isinstance(row, dict) and
                 isinstance(row.get("relative_path"), str) and
                 row["relative_path"] in expected_relatives and
                 row["relative_path"] not in by_relative,
                 "recovery derived-prefix relative path is invalid")
        relative = row["relative_path"]
        source_path = _safe_output_path(original_state, relative)
        derived_path = _safe_output_path(recovery_state, relative)
        source_binding = _verify_planned_file(
            row.get("source"), f"recovery original text {relative}")
        derived_binding = _verify_planned_file(
            row.get("derived"), f"recovery derived text {relative}")
        _require(source_binding["path"] == str(source_path) and
                 derived_binding["path"] == str(derived_path),
                 "recovery text prefix binding uses the wrong tree/path")
        source_raw = _read_stable_bytes(source_path)
        derived_raw = _read_stable_bytes(derived_path)
        prefix_size = _integer(row.get("prefix_size"), "recovery prefix_size")
        suffix_size = _integer(row.get("suffix_size"), "recovery suffix_size")
        _require(prefix_size == len(derived_raw) and
                 prefix_size + suffix_size == len(source_raw) and
                 source_raw[:prefix_size] == derived_raw and
                 row.get("prefix_sha256") == hashlib.sha256(derived_raw).hexdigest() and
                 row.get("suffix_sha256") ==
                 hashlib.sha256(source_raw[prefix_size:]).hexdigest() and
                 derived_raw.endswith(b"\n"),
                 "recovery derived text is not the exact newline-closed source prefix")
        output = output_for_required[relative]
        if output["file_type"] == "hst":
            semantic = audit_history(
                derived_path, original_plan["expected"]["source_cycle"],
                final_cycle, original_plan["expected"]["source_time"],
                original_plan["expected"]["root_dt"], final_time,
                original_plan["policy"]["endpoint_time_ulp_tolerance"])
            persisted_semantic = dict(semantic)
            persisted_semantic.pop("_row_values", None)
        elif output["file_type"] == "log":
            semantic = audit_event_log(
                derived_path, original_plan["expected"]["source_cycle"],
                final_cycle, original_plan["policy"]["event_thresholds"],
                original_plan["policy"]["event_absolute_thresholds"])
            persisted_semantic = dict(semantic)
            persisted_semantic.pop("_rows", None)
        else:
            raise CheckFailure(
                f"unsupported recovery unnumbered type: {output['file_type']}")
        by_relative[relative] = {"source": source_binding,
                                 "derived": derived_binding,
                                 "prefix_size": prefix_size,
                                 "suffix_size": suffix_size,
                                 "block": output["block"],
                                 "file_type": output["file_type"],
                                 "semantic": semantic,
                                 "persisted_semantic": persisted_semantic}
    _require(set(by_relative) == expected_relatives,
             "recovery derived-prefix inventory omits required text")
    unnumbered = [row for row in parent.get("output_inventory", [])
                  if isinstance(row, dict) and row.get("file_number") is None]
    inventory_by_path = {row.get("path"): row for row in unnumbered}
    _require(len(inventory_by_path) == len(unnumbered) == len(by_relative),
             "recovery output inventory has duplicate/extra unnumbered entries")
    for evidence in by_relative.values():
        derived = evidence["derived"]
        inventory = inventory_by_path.get(derived["path"])
        _require(isinstance(inventory, dict) and
                 _same_full_binding(inventory, derived) and
                 inventory.get("block") == evidence["block"] and
                 inventory.get("file_type") == evidence["file_type"] and
                 inventory.get("file_number") is None and
                 inventory.get("expected_write") is None and
                 ((evidence["file_type"] == "hst" and
                   inventory.get("history") == evidence["persisted_semantic"]) or
                  (evidence["file_type"] == "log" and
                   inventory.get("event_log") ==
                   evidence["persisted_semantic"])),
                 "recovery output inventory does not bind every derived prefix")
    expected_numbered = set(expected_numbered_rows)
    numbered_rows = [row for row in parent.get("output_inventory", [])
                     if isinstance(row, dict) and row.get("file_number") is not None]
    state_rows: dict[str, dict[str, Any]] = {}
    binary_by_block: dict[str, list[dict[str, Any]]] = {}
    divb_rows: list[dict[str, Any]] = []
    full_by_cycle: dict[int, dict[str, Any]] = {}
    for row in numbered_rows:
        path_value = row.get("path")
        _require(isinstance(path_value, str),
                 "recovery numbered output lacks an absolute path")
        path = Path(path_value).resolve(strict=True)
        try:
            relative = path.relative_to(original_state).as_posix()
        except ValueError as exc:
            raise CheckFailure(
                "recovery numbered output is outside original state") from exc
        _require(relative in expected_numbered and relative not in state_rows,
                 "recovery numbered output is duplicate or outside logical prefix")
        live = _verify_planned_file(row, f"recovery numbered output {relative}")
        _require(_same_full_binding(row, live),
                 f"recovery numbered output binding changed: {relative}")
        expected_row = expected_numbered_rows[relative]
        _require(all(row.get(name) == expected_row[name] for name in (
                     "block", "file_type", "file_number", "expected_write")),
                 f"recovery numbered output metadata differs from replay: {relative}")
        output = next(item for item in original_plan["outputs"]
                      if item["block"] == expected_row["block"])
        if output["file_type"] == "rst":
            audited_restart = audit_restart(path)
            _require(audited_restart.get("valid") is True and
                     audited_restart.get("metadata", {}).get("cycle") ==
                     expected_row["expected_write"]["cycle"] and
                     audited_restart.get("metadata", {}).get("time") ==
                     expected_row["expected_write"]["time"] and
                     row.get("restart_audit") == audited_restart,
                     f"recovery prefix restart is scientifically invalid: {relative}")
        elif output["file_type"] == "bin":
            audited_binary = audit_binary(
                path, original_plan["expected"]["source_cycle"], final_cycle,
                original_plan["expected"]["source_time"], final_time,
                output["parameters"], endpoint_metadata,
                source_metadata.parameters, output["block"],
                output["expected_binary_variables"])
            _require(audited_binary["cycle"] ==
                     expected_row["expected_write"]["cycle"] and
                     audited_binary["time"] ==
                     expected_row["expected_write"]["time"],
                     f"recovery prefix binary differs from replay: {relative}")
            binary_by_block.setdefault(output["block"], []).append(audited_binary)
            if audited_binary["topology"]["scope"] == "full_domain":
                previous = full_by_cycle.get(audited_binary["cycle"])
                _require(previous is None or
                         (previous["time"] == audited_binary["time"] and
                          previous["_topology_records"] ==
                          audited_binary["_topology_records"]),
                         "recovery full-domain binary topologies disagree")
                full_by_cycle.setdefault(audited_binary["cycle"], audited_binary)
            for variable, maximum in audited_binary["variable_max_abs"].items():
                if variable.lower() == "divb":
                    divb_rows.append({
                        "path": str(path), "cycle": audited_binary["cycle"],
                        "time": audited_binary["time"], "max_abs": maximum,
                    })
        state_rows[relative] = live
    _require(set(state_rows) == expected_numbered,
             "recovery output inventory does not exactly cover numbered prefix")
    audit_binary_pairs(original_plan, binary_by_block)
    bbh_histories = [evidence["semantic"] for evidence in by_relative.values()
                     if evidence["file_type"] == "hst" and
                     {"bh1_x", "bh2_x", "bh1_mass", "bh2_mass"}.issubset(
                         evidence["semantic"].get("columns", []))]
    _require(len(bbh_histories) == 1,
             "recovery prefix lacks exactly one BBH center history")
    centers = _bbh_history_centers(bbh_histories[0])
    reported_coverage = parent.get("binary_selection_coverage_audit")
    _require(isinstance(reported_coverage, list),
             "recovery pass lacks binary selection coverage audit")
    coverage: list[dict[str, Any]] = []
    for block, audited_rows in binary_by_block.items():
        parameters = next(output["parameters"] for output in original_plan["outputs"]
                          if output["block"] == block)
        for audited in audited_rows:
            result = audit_selected_binary_coverage(
                audited, parameters, centers.get(audited["cycle"]),
                full_by_cycle.get(audited["cycle"]))
            coverage.append({"path": audited["path"], "cycle": audited["cycle"],
                             **result})
    _require(reported_coverage == coverage,
             "recovery pass binary selection coverage differs from live audit")
    _require(divb_rows and
             max(row["max_abs"] for row in divb_rows) <=
             original_plan["policy"]["divb_max_abs"]["hard"],
             "recovery prefix divB exceeds hard threshold")
    reported_divb = parent.get("scientific_threshold_audit", {}).get("divb")
    _require(isinstance(reported_divb, dict) and
             reported_divb.get("files") == divb_rows and
             reported_divb.get("observed_max_abs") ==
             max(row["max_abs"] for row in divb_rows) and
             reported_divb.get("hard_max_abs") ==
             original_plan["policy"]["divb_max_abs"]["hard"],
             "recovery pass divB audit differs from live binary audit")

    known_numbered: set[str] = set()
    for output in original_plan["outputs"]:
        if not output.get("numbered"):
            continue
        for write in output["expected_writes"]:
            known_numbered.add(output["relative_path_template"].format(
                file_number=write["file_number"]))
    suffix_rows = record.get("suffix_forensics")
    _require(isinstance(suffix_rows, list),
             "recovery record lacks suffix forensics")
    suffix_by_relative: dict[str, dict[str, Any]] = {}
    for row in suffix_rows:
        _require(isinstance(row, dict) and
                 isinstance(row.get("relative_path"), str) and
                 row["relative_path"] in known_numbered - expected_numbered and
                 row["relative_path"] not in suffix_by_relative and
                 row.get("classification") == "planned_post_prefix_artifact",
                 "recovery suffix forensic row is invalid")
        relative = row["relative_path"]
        expected_path = _safe_output_path(original_state, relative)
        binding = _verify_planned_file(
            row.get("binding"), f"recovery suffix artifact {relative}")
        _require(binding["path"] == str(expected_path) and
                 _same_full_binding(row.get("binding"), binding),
                 f"recovery suffix binding changed: {relative}")
        suffix_by_relative[relative] = binding

    actual_original: set[str] = set()
    actual_directories: set[str] = {"."}
    for directory, directories, files in os.walk(original_state, followlinks=False):
        directory_path = Path(directory)
        for name in directories:
            child = directory_path / name
            _require(not child.is_symlink(),
                     "recovery original state contains a symlink directory")
            actual_directories.add(child.relative_to(original_state).as_posix())
        for name in files:
            path = directory_path / name
            info = path.lstat()
            _require(stat.S_ISREG(info.st_mode) and not stat.S_ISLNK(info.st_mode),
                     "recovery original state contains a non-regular node")
            actual_original.add(path.relative_to(original_state).as_posix())
    expected_original = expected_numbered | set(suffix_by_relative) | \
        expected_relatives
    _require(actual_original == expected_original,
             "recovery original state tree differs from prefix/suffix provenance")
    expected_directories = {"."}
    for relative in expected_original:
        parent_path = Path(relative).parent
        while parent_path.as_posix() != ".":
            expected_directories.add(parent_path.as_posix())
            parent_path = parent_path.parent
    _require(actual_directories == expected_directories,
             "recovery original state contains an unknown directory")
    return {
        "derived_prefixes": len(by_relative),
        "relative_paths": sorted(by_relative),
        "numbered_prefixes": len(state_rows),
        "suffix_artifacts": len(suffix_by_relative),
        "original_state_exact": True,
    }


def _audit_recovery_original_launch(
        parent: dict[str, Any], provenance: dict[str, Any],
        record: dict[str, Any], original_plan: dict[str, Any],
        original_plan_path: Path) -> dict[str, Any]:
    """Revalidate the complete original plan/launch/tool trust anchor."""

    validate_plan(original_plan)
    _require(original_plan.get("policy", {}).get(
                 "scheduled_prefix_recovery") == PREFIX_RECOVERY_POLICY,
             "recovery original plan lacks the exact recovery policy")
    _require(original_plan_path.resolve(strict=True) == Path(
                 original_plan["launch_contract"]["plan_path"]).resolve(strict=True),
             "recovery original plan path differs from its launch contract")
    state_dir = Path(original_plan["launch_contract"]["state_dir"]).resolve(
        strict=True)
    original_launch = provenance.get("original_launch_record")
    _validate_planned_file_record(
        original_launch, "recovery_provenance.original_launch_record")
    _require(_same_full_binding(
        original_launch, parent.get("bindings", {}).get("original_launch_record")) and
             _same_full_binding(original_launch,
                                record.get("original_launch_record")),
             "recovery original-launch bindings disagree")
    launch_binding = _verify_planned_file(
        original_launch, "recovery original launch record")
    launch_path = Path(launch_binding["path"])
    _require_immutable_ready(
        launch_path, ".launch.ready", "recovery original launch record")
    planned_launch_path = Path(
        original_plan["launch_contract"]["evidence"]["launch_record"])
    _require(launch_path.resolve(strict=True) == planned_launch_path.resolve(strict=True),
             "recovery launch record path differs from original plan")
    gpu_before_path = Path(
        original_plan["launch_contract"]["evidence"]["gpu_before"])
    gpu_before = _hash_record(gpu_before_path)
    launch_audit = audit_launch_record(
        launch_path, original_plan, original_plan_path, state_dir, gpu_before)
    _require(launch_audit["sha256"] == launch_binding["sha256"],
             "recovery launch record changed during canonical audit")

    record_tools = record.get("tools")
    parent_tools = parent.get("bindings", {}).get("tools")
    _require(isinstance(record_tools, dict) and isinstance(parent_tools, dict),
             "recovery record/pass lacks tool bindings")
    live_tools = {
        name: _verify_planned_file(original_plan["tools"][name],
                                   f"recovery original-plan tool {name}")
        for name in PLANNED_TOOL_NAMES
    }
    for name, live in live_tools.items():
        _require(_same_planned_file(record_tools.get(name),
                                    original_plan["tools"][name]) and
                 _same_full_binding(record_tools.get(name), live) and
                 _same_full_binding(parent_tools.get(name), live),
                 f"recovery {name} tool does not match original plan/live bytes")
    _verify_planned_repository(
        original_plan["inputs"]["repo"], original_plan["tools"]["git"])
    launch_payload, _ = _load_json(launch_path)
    return {
        "original_plan_validated": True,
        "original_launch_record": launch_binding,
        "launch_audit": launch_audit,
        "launch_record_payload": launch_payload,
        "plan_bound_tools": live_tools,
        "state_dir": str(state_dir),
    }


def _audit_recovery_directories(
        record: dict[str, Any], original_plan: dict[str, Any],
        record_path: Path, parent_path: Path | None) -> dict[str, Any]:
    """Rebind all four recovery directories and their canonical artifact paths."""

    def live(path: Path, label: str, *, owner_only: bool) -> dict[str, Any]:
        absolute = path.expanduser().absolute()
        info = absolute.lstat()
        _require(stat.S_ISDIR(info.st_mode) and not stat.S_ISLNK(info.st_mode) and
                 absolute == absolute.resolve(strict=True) and
                 info.st_uid == os.geteuid() and
                 (not owner_only or stat.S_IMODE(info.st_mode) == 0o700),
                 f"{label} is not a canonical owner-owned directory")
        return {
            "path": str(absolute), "device": info.st_dev, "inode": info.st_ino,
            "owner_uid": info.st_uid, "mode": f"{stat.S_IMODE(info.st_mode):04o}",
            "descriptor_bound": True,
        }

    original_state = Path(original_plan["launch_contract"]["state_dir"])
    original_evidence = Path(original_plan["launch_contract"]["evidence_dir"])
    recovery_state_value = record.get("recovery_state_dir")
    recovery_evidence_value = record.get("recovery_evidence_dir")
    _require(isinstance(recovery_state_value, str) and
             isinstance(recovery_evidence_value, str),
             "recovery record lacks canonical recovery directories")
    recovery_state = Path(recovery_state_value)
    recovery_evidence = Path(recovery_evidence_value)
    checks = {
        "original_state": (
            live(original_state, "recovery original state", owner_only=False),
            record.get("original_state_tree", {}).get("directory_binding")),
        "original_evidence": (
            live(original_evidence, "recovery original evidence", owner_only=False),
            record.get("original_evidence_directory_binding")),
        "recovery_state": (
            live(recovery_state, "recovery derived state", owner_only=True),
            record.get("recovery_state_directory_binding")),
        "recovery_evidence": (
            live(recovery_evidence, "recovery evidence", owner_only=True),
            record.get("recovery_evidence_directory_binding")),
    }
    for label, (observed, planned) in checks.items():
        _require(observed == planned,
                 f"{label} directory identity differs from recovery record")
    _require(record_path.resolve(strict=True) ==
             (recovery_evidence / "segment.prefix.recovery.ready").resolve(
                 strict=True),
             "recovery record is outside its exact recovery evidence path")
    _require(parent_path is not None and parent_path.resolve(strict=True) ==
             (recovery_evidence / "segment.prefix.pass.ready").resolve(strict=True),
             "recovery parent pass is outside its exact recovery evidence path")
    actual_evidence = {entry.name for entry in os.scandir(recovery_evidence)}
    _require(actual_evidence == {
        "segment.prefix.recovery.ready", "segment.prefix.pass.ready"},
        "recovery evidence directory contains an unknown node")
    return {name: observed for name, (observed, _) in checks.items()}


def _audit_stored_live_gpu(value: Any, plan: dict[str, Any],
                           before_path: Path) -> dict[str, Any]:
    tool = _verify_planned_file(
        plan["tools"]["nvidia_smi"], "recovery plan-bound nvidia-smi")
    before = parse_gpu_csv(before_path)
    _require(isinstance(value, dict) and
             value.get("compute_contexts") == [] and
             value.get("compute_contexts_empty") is True and
             value.get("ecc_all_zero") is True and
             _same_full_binding(value.get("nvidia_smi"),
                                tool) and
             value.get("launch_inventory") == before,
             "stored recovery live-GPU quiescence proof is invalid")
    inventory = value.get("inventory")
    _require(isinstance(inventory, list) and
             len(inventory) == plan["policy"]["ranks"] and
             [(row.get("index"), row.get("uuid"), row.get("pci_bus_id"),
              row.get("cuda_ordinal"), row.get("memory_total_mib"))
              for row in inventory] ==
             [(row["index"], row["uuid"], row["pci_bus_id"],
               row["cuda_ordinal"], row["memory_total_mib"])
              for row in before] and
             all(row.get("uncorrected_ecc") == 0 and
                 row.get("corrected_ecc") == 0 and
                 isinstance(row.get("memory_used_mib"), int) and
                 not isinstance(row.get("memory_used_mib"), bool) and
                 0 <= row["memory_used_mib"] <= row.get("memory_total_mib", -1) and
                 row["memory_used_mib"] <=
                 plan["policy"]["gpu_exit_memory_mib_max"]
                 for row in inventory),
             "stored recovery GPU identity/ECC/memory proof is invalid")
    return {"nvidia_smi": tool, "launch_inventory": before,
            "producer_inventory": inventory,
            "compute_contexts_empty": True, "ecc_all_zero": True}


def _audit_recovery_lifecycle(
        lifecycle: dict[str, Any], record: dict[str, Any],
        plan: dict[str, Any], anchor: dict[str, Any]) -> dict[str, Any]:
    """Re-prove lifecycle closure from the canonically audited original launch."""

    launch = anchor["launch_record_payload"]
    launch_audit = anchor["launch_audit"]
    evidence = plan["launch_contract"]["evidence"]
    expected_identities = [{
        "role": "mpirun", "pid": launch_audit["mpirun_pid"],
        "recorded_start_time_ticks": launch_audit["mpirun_start_time_ticks"],
    }, {
        "role": "launcher_holder",
        "pid": launch_audit["input_transport"]["holder_pid"],
        "recorded_start_time_ticks":
            launch_audit["input_transport"]["holder_start_time_ticks"],
    }]
    expected_identities.extend({
        "role": f"rank:{rank['global_rank']}", "pid": rank["pid"],
        "recorded_start_time_ticks": rank["start_time_ticks"],
    } for rank in launch_audit["ranks"])

    reported_identities = (lifecycle.get("process_identities")
                           if lifecycle["kind"] ==
                           "same_boot_closed_processes_v1"
                           else lifecycle.get("live_identity_recheck"))
    _require(isinstance(reported_identities, list) and
             len(reported_identities) == len(expected_identities),
             "recovery lifecycle process identity set is incomplete")
    for reported, expected in zip(reported_identities, expected_identities):
        _require(isinstance(reported, dict) and
                 all(reported.get(name) == expected[name]
                     for name in ("role", "pid", "recorded_start_time_ticks")) and
                 reported.get("original_identity_gone") is True and
                 reported.get("state") in ("disappeared", "pid_reused") and
                 ((reported.get("state") == "disappeared" and
                   reported.get("observed_start_time_ticks") is None) or
                  (reported.get("state") == "pid_reused" and
                   isinstance(reported.get("observed_start_time_ticks"), int) and
                   reported["observed_start_time_ticks"] !=
                   reported["recorded_start_time_ticks"])),
                 "recovery lifecycle identity differs from original launch")

    gpu_before = Path(evidence["gpu_before"])
    producer_gpu = _audit_stored_live_gpu(
        lifecycle.get("live_gpu_quiescence"), plan, gpu_before)

    evidence_bindings: dict[str, dict[str, Any]] = {
        "plan": record["original_plan"],
        "launch_record": record["original_launch_record"],
    }
    expected_paths = {
        "plan": plan["launch_contract"]["plan_path"],
        **evidence,
    }
    if lifecycle["kind"] == "nonzero_completion_v1":
        return_code = lifecycle.get("return_code")
        _require(isinstance(return_code, int) and
                 not isinstance(return_code, bool) and return_code != 0,
                 "recovery nonzero return code is invalid")
        artifacts = lifecycle.get("artifacts")
        _require(isinstance(artifacts, dict) and set(artifacts) == {
            "plan", "launch_record", "run_log", "exit_status",
            "gpu_before", "gpu_after",
        }, "recovery nonzero lifecycle artifact set is not exact")
        _require(_same_full_binding(artifacts.get("plan"),
                                    record["original_plan"]) and
                 _same_full_binding(artifacts.get("launch_record"),
                                    record["original_launch_record"]),
                 "recovery completion overwrites trusted plan/launch binding")
        evidence_bindings.update(artifacts)
        completion = lifecycle.get("completion_record")
        _validate_planned_file_record(
            completion, "recovery lifecycle completion_record")
        completion_live = _verify_planned_file(
            completion, "recovery lifecycle completion_record")
        completion_path = Path(completion_live["path"])
        _require(completion_path.resolve(strict=True) ==
                 Path(evidence["completion_record"]).resolve(strict=True),
                 "recovery completion path differs from original plan")
        try:
            exit_value = int(Path(artifacts["exit_status"]["path"]).read_text(
                encoding="ascii").strip())
        except (OSError, UnicodeDecodeError, ValueError) as exc:
            raise CheckFailure("recovery exit.status is invalid") from exc
        _require(exit_value == return_code,
                 "recovery completion return code differs from exit.status")
        canonical = audit_completion_record(
            completion_path, plan,
            Path(plan["launch_contract"]["state_dir"]), launch_audit,
            artifacts, expected_return_code=return_code)
        gpu_audit = audit_gpus(
            Path(artifacts["gpu_before"]["path"]),
            Path(artifacts["gpu_after"]["path"]), plan["policy"]["ranks"],
            plan["policy"]["gpu_exit_memory_mib_max"],
            plan["policy"]["gpu_ecc"])
        _require(lifecycle.get("canonical_completion_audit") == canonical and
                 lifecycle.get("gpu_audit") == gpu_audit,
                 "recovery canonical completion/GPU audit differs from live replay")
        evidence_bindings["completion_record"] = completion_live
    else:
        bounded = lifecycle.get("bounded_quiescence")
        _require(lifecycle.get("hostname") == launch.get("hostname") and
                 lifecycle.get("all_original_identities_gone") is True and
                 lifecycle.get("managed_process_group") ==
                 launch_audit.get("managed_process_group") and
                 lifecycle.get("managed_process_group_gone") is True and
                 isinstance(bounded, dict) and
                 bounded.get("kind") ==
                 "poll_all_identities_group_and_gpu_v1" and
                 bounded.get("timeout_seconds") ==
                 SAME_BOOT_QUIESCENCE_TIMEOUT_SECONDS and
                 bounded.get("poll_seconds") ==
                 SAME_BOOT_QUIESCENCE_POLL_SECONDS and
                 isinstance(bounded.get("attempts"), int) and
                 not isinstance(bounded.get("attempts"), bool) and
                 bounded["attempts"] > 0 and
                 bounded.get("all_quiet") is True and
                 not os.path.lexists(evidence["completion_record"]),
                 "recovery same-boot lifecycle host/completion proof is invalid")
        try:
            launch_epoch = __import__("datetime").datetime.fromisoformat(
                launch["created_utc"].replace("Z", "+00:00")).timestamp()
        except (KeyError, TypeError, ValueError) as exc:
            raise CheckFailure("cannot parse original launch timestamp") from exc
        _require(isinstance(lifecycle.get("boot_id"), str) and
                 re.fullmatch(r"[0-9a-f-]{36}", lifecycle["boot_id"]) is not None and
                 isinstance(lifecycle.get("boot_epoch"), int) and
                 lifecycle.get("launch_epoch") == launch_epoch and
                 isinstance(lifecycle.get("launch_epoch"), (int, float)) and
                 lifecycle["boot_epoch"] <= lifecycle["launch_epoch"],
                 "recovery same-boot epoch/ID proof is invalid")
        closed = lifecycle.get("closed_evidence_artifacts")
        _require(isinstance(closed, dict) and
                 {"run_log", "gpu_before"}.issubset(closed) and
                 set(closed).issubset({
                     "run_log", "gpu_before", "exit_status", "gpu_after"}),
                 "recovery same-boot closed-evidence set is invalid")
        evidence_bindings.update(closed)

    expected_by_name: dict[str, dict[str, Any]] = {}
    for role, binding in evidence_bindings.items():
        _validate_planned_file_record(binding, "recovery original evidence")
        _require(role in expected_paths and
                 Path(binding["path"]).resolve(strict=True) ==
                 Path(expected_paths[role]).resolve(strict=True),
                 f"recovery {role} path differs from original plan")
        name = Path(binding["path"]).name
        _require(name not in expected_by_name or
                 _same_full_binding(expected_by_name[name], binding),
                 "recovery original evidence basename collision")
        expected_by_name[name] = binding
    evidence_dir = Path(plan["launch_contract"]["evidence_dir"]).resolve(strict=True)
    _require(all(Path(binding["path"]).resolve(strict=True).parent == evidence_dir
                 for binding in evidence_bindings.values()),
             "recovery lifecycle evidence is outside original evidence_dir")
    actual_names: set[str] = set()
    for entry in os.scandir(evidence_dir):
        info = entry.stat(follow_symlinks=False)
        _require(stat.S_ISREG(info.st_mode) and not stat.S_ISLNK(info.st_mode),
                 "recovery original evidence contains a non-regular node")
        actual_names.add(entry.name)
    _require(actual_names == set(expected_by_name),
             "recovery original evidence directory set is not exact")
    for name, binding in expected_by_name.items():
        live = _verify_planned_file(binding, f"recovery evidence {name}")
        _require(_same_full_binding(live, binding),
                 f"recovery original evidence binding changed: {name}")
    return {
        "kind": lifecycle["kind"],
        "launch_identity_set_exact": True,
        "producer_gpu_quiescence": producer_gpu,
        "portable_static_revalidation": True,
        "original_evidence_exact": True,
    }


def _audit_recovery_source_chain(
        record: dict[str, Any], plan: dict[str, Any],
        seen: set[tuple[str, str]], depth: int) -> dict[str, Any]:
    """Rebuild and recursively qualify the original segment's source chain."""

    _require(depth <= 64, "recovery source chain exceeds maximum depth")
    qualification = plan["source_qualification"]
    source_path = Path(plan["inputs"]["source_restart"]["path"])
    source_audit = _verify_endpoint_audit(audit_restart(source_path), source_path)
    source_binding = {
        "path": source_audit["path"], **source_audit["signature"],
        "sha256": source_audit["sha256"],
        "closure_check": source_audit["closure_check"],
    }
    _require(_same_planned_file(source_binding,
                                plan["inputs"]["source_restart"]),
             "recovery original source differs from original plan")
    audit_restart_contract(
        read_restart_metadata(source_path), "recovery original source",
        plan, plan["expected"])
    _require(source_audit["layout"] == qualification["audit"]["layout"] and
             source_audit["stored_reals"] ==
             qualification["audit"]["stored_reals"] and
             source_audit["topology"] == qualification["audit"]["topology"],
             "recovery original source audit differs from original plan")
    mass = _verify_planned_file(
        qualification["source_baryon_mass"]["evidence"],
        "recovery original source baryon evidence")
    if qualification["mode"] == "anchor_full_audit":
        audit_source_baryon_evidence(
            qualification, plan["expected"]["source_time"])
        expected = {
            "mode": "anchor_full_audit",
            "source_baryon_mass_evidence": mass,
        }
        _require(record.get("source_chain") == expected,
                 "recovery anchor source_chain differs from live evidence")
        return expected

    parent_binding = _verify_planned_file(
        qualification["parent_segment_pass"],
        "recovery original source parent pass")
    identity = (parent_binding["path"], parent_binding["sha256"])
    _require(identity not in seen,
             "recovery source qualification chain contains a cycle")
    parent_path = Path(parent_binding["path"])
    _require_immutable_ready(parent_path, ".pass.ready",
                             "recovery source parent pass")
    parent, _ = _load_json(parent_path)
    audit_source_baryon_evidence(
        qualification, plan["expected"]["source_time"], parent)
    endpoint = parent.get("bindings", {}).get("endpoint_restart")
    trajectory = parent.get("bindings", {}).get("trajectory")
    _require(_same_full_binding(endpoint, source_binding) and
             _same_planned_file(trajectory, plan["inputs"]["trajectory"]),
             "recovery source parent does not bind original source/trajectory")
    histories = [row for row in parent.get("output_inventory", [])
                 if isinstance(row, dict) and
                 isinstance(row.get("history"), dict) and
                 "baryon_m" in row["history"].get("columns", [])]
    _require(len(histories) == 1 and _same_full_binding(histories[0], mass),
             "recovery source parent does not bind baryon evidence")
    provenance = audit_parent_qualification_provenance(
        parent, qualification, source_binding, mass,
        plan["expected"]["source_cycle"], plan["expected"]["source_time"],
        parent_path, _chain_seen=seen | {identity}, _chain_depth=depth + 1)
    expected = {
        "mode": "parent_segment_pass", "parent_pass": parent_binding,
        "source_baryon_mass_evidence": mass,
        "parent_qualification": provenance,
    }
    _require(record.get("source_chain") == expected,
             "recovery source_chain differs from recursive live qualification")
    return expected


def audit_parent_qualification_provenance(
        parent: dict[str, Any], source_qualification: dict[str, Any],
        source_binding: dict[str, Any], source_mass_evidence: dict[str, Any],
        source_cycle: int, source_time: float,
        parent_path: Path | None = None, *,
        _chain_seen: set[tuple[str, str]] | None = None,
        _chain_depth: int = 0) -> dict[str, Any]:
    """Re-prove the explicit complete/recovery mode of an immutable parent pass."""

    mode = parent.get("qualification_mode")
    if mode is None:
        _require("recovery_provenance" not in parent and
                 parent.get("completion_record_audit", {}).get("return_code") == 0 and
                 parent.get("run_log_audit", {}).get("termination") == "cycle_limit",
                 "legacy parent without qualification_mode lacks a complete "
                 "zero-exit cycle-limit lifecycle")
        mode = "legacy_complete_segment_v1"
    _require(mode == source_qualification.get("parent_qualification_mode") and
             mode in ("complete_segment_v1", "legacy_complete_segment_v1",
                      "scheduled_prefix_recovery_v1"),
             "parent pass qualification mode differs from the source plan")
    if mode in ("complete_segment_v1", "legacy_complete_segment_v1"):
        _require("recovery_provenance" not in parent,
                 "complete parent pass may not contain recovery provenance")
        return {"qualification_mode": mode}

    provenance = parent.get("recovery_provenance")
    advisories = parent.get("scientific_advisories")
    _require(isinstance(provenance, dict) and
             provenance.get("kind") == "scheduled_prefix_recovery_v1" and
             provenance.get("policy") == PREFIX_RECOVERY_POLICY and
             provenance.get("original_trees_unchanged") is True,
             "recovery parent policy/provenance is invalid")
    _require(isinstance(advisories, dict) and
             advisories.get("schema") == "athenak_scientific_advisories_v1" and
             advisories.get("severity") in ("green", "yellow") and
             advisories.get("pass_fail_effect") ==
             "none_yellow_advisories_are_nonfatal" and
             isinstance(advisories.get("floor_rates"), dict),
             "recovery parent lacks explicit scientific advisories")
    selected = provenance.get("selected_scheduled_restart")
    logical = provenance.get("logical_prefix")
    lifecycle = provenance.get("lifecycle")
    run_log_prefix = provenance.get("run_log_prefix_audit")
    _require(isinstance(selected, dict) and
             selected.get("expected_write", {}).get("kind") == "scheduled" and
             selected.get("expected_write", {}).get("cycle") == source_cycle and
             selected.get("expected_write", {}).get("time") == source_time and
             isinstance(logical, dict) and
             logical.get("final_cycle") == source_cycle and
             logical.get("final_time") == source_time and
             logical.get("execute_only") is True and
             logical.get("all_expected_prefix_artifacts_present") is True,
             "recovery parent scheduled endpoint/logical-prefix proof is invalid")
    _require(isinstance(run_log_prefix, dict) and
             run_log_prefix.get("root_step_diagnostics", {}).get("cycle_max") ==
             source_cycle and
             run_log_prefix.get("root_step_diagnostics", {}).get(
                 "all_recovered_prefix_cycles_present") is True and
             run_log_prefix.get("cache", {}).get("solver_failures") == 0 and
             run_log_prefix.get("cache", {}).get(
                 "nonfinite_proposed_values") == 0,
             "recovery parent run-log/cache prefix proof is invalid")
    _require(_same_full_binding(selected.get("binding"), source_binding),
             "recovery selected restart differs from the live source restart")
    _require(isinstance(lifecycle, dict) and lifecycle.get("kind") in
             ("nonzero_completion_v1", "same_boot_closed_processes_v1"),
             "recovery parent lifecycle proof is invalid")
    if lifecycle["kind"] == "nonzero_completion_v1":
        _require(isinstance(lifecycle.get("return_code"), int) and
                 not isinstance(lifecycle.get("return_code"), bool) and
                 lifecycle["return_code"] != 0 and
                 isinstance(lifecycle.get("completion_record"), dict),
                 "recovery parent nonzero-completion proof is invalid")
    else:
        _require(lifecycle.get("all_original_identities_gone") is True and
                 isinstance(lifecycle.get("boot_id"), str) and
                 isinstance(lifecycle.get("process_identities"), list) and
                 all(isinstance(row, dict) and
                     row.get("original_identity_gone") is True
                     for row in lifecycle["process_identities"]),
                 "recovery parent same-boot closure proof is invalid")

    recovery_record = provenance.get("record")
    _validate_planned_file_record(recovery_record, "recovery_provenance.record")
    _require(_same_full_binding(
        recovery_record, parent.get("bindings", {}).get("recovery_record")),
        "recovery record bindings disagree inside parent pass")
    record_binding = _verify_planned_file(recovery_record, "recovery record")
    record_path = Path(record_binding["path"])
    _require_immutable_ready(record_path, ".recovery.ready", "recovery record")
    record, _ = _load_json(record_path)
    original_plan = provenance.get("original_plan")
    _validate_planned_file_record(
        original_plan, "recovery_provenance.original_plan")
    _require(_same_full_binding(
        original_plan, parent.get("bindings", {}).get("original_plan")),
        "recovery original-plan bindings disagree inside parent pass")
    original_plan_binding = _verify_planned_file(
        original_plan, "recovery original plan")
    original_plan_path = Path(original_plan_binding["path"])
    _require_immutable_ready(
        original_plan_path, "segment.plan.json", "recovery original plan")
    original_plan_payload, _ = _load_json(original_plan_path)
    _require(record.get("schema") == SCHEMA and
             record.get("kind") == "athenak_segment_prefix_recovery" and
             record.get("status") == "prepared" and
             record.get("publication_transaction") == {
                 "kind": "prepared_record_then_pass_commit_v1",
                 "commit_filename": "segment.prefix.pass.ready",
                 "prepared_record_alone_is_not_consumable": True,
             } and
             record.get("qualification_mode") == mode and
             record.get("policy") == PREFIX_RECOVERY_POLICY and
             record.get("original_plan") == original_plan and
             record.get("original_launch_record") ==
             provenance.get("original_launch_record") and
             record.get("selected_scheduled_restart") == selected and
             record.get("logical_prefix") == logical and
             record.get("run_log_prefix_audit") == run_log_prefix and
             record.get("lifecycle") == lifecycle and
             record.get("derived_text_prefixes") ==
             provenance.get("derived_text_prefixes") and
             record.get("suffix_forensics") ==
             provenance.get("suffix_forensics"),
             "immutable recovery record differs from parent provenance")
    original_anchor_audit = _audit_recovery_original_launch(
        parent, provenance, record, original_plan_payload, original_plan_path)
    directory_audit = _audit_recovery_directories(
        record, original_plan_payload, record_path, parent_path)
    lifecycle_audit = _audit_recovery_lifecycle(
        lifecycle, record, original_plan_payload, original_anchor_audit)
    source_chain_audit = _audit_recovery_source_chain(
        record, original_plan_payload, _chain_seen or set(), _chain_depth)
    original_anchor_audit = {
        name: value for name, value in original_anchor_audit.items()
        if name != "launch_record_payload"
    }
    candidate_audit = _audit_recovery_candidates(
        original_plan_payload, record, selected)
    prefix_audit = _audit_recovery_derived_prefixes(
        original_plan_payload, record, parent)
    derived = provenance.get("derived_text_prefixes")
    _require(isinstance(derived, list) and any(
        isinstance(row, dict) and
        _same_full_binding(row.get("derived"), source_mass_evidence)
        for row in derived),
        "source baryon history is not a recovery-derived byte prefix")
    history_path = Path(source_mass_evidence["path"])
    history_info = history_path.lstat()
    recovery_state_value = record.get("recovery_state_dir")
    _require(isinstance(recovery_state_value, str),
             "recovery record lacks recovery_state_dir")
    recovery_state = Path(recovery_state_value).resolve(strict=True)
    recovery_state_info = recovery_state.lstat()
    try:
        history_path.resolve(strict=True).relative_to(recovery_state)
    except ValueError as exc:
        raise CheckFailure(
            "recovery baryon prefix is outside recovery state") from exc
    _require(stat.S_ISREG(history_info.st_mode) and
             not stat.S_ISLNK(history_info.st_mode) and
             not (stat.S_IMODE(history_info.st_mode) & 0o222) and
             stat.S_ISDIR(recovery_state_info.st_mode) and
             not stat.S_ISLNK(recovery_state_info.st_mode) and
             recovery_state_info.st_uid == os.geteuid() and
             stat.S_IMODE(recovery_state_info.st_mode) == 0o700,
             "recovery state/history must remain owner-only and immutable")
    return {
        "qualification_mode": mode,
        "recovery_record": record_binding,
        "selected_scheduled_restart": selected,
        "logical_prefix": logical,
        "lifecycle": lifecycle,
        "original_anchor_audit": original_anchor_audit,
        "directory_audit": directory_audit,
        "lifecycle_audit": lifecycle_audit,
        "source_chain_audit": source_chain_audit,
        "candidate_audit": candidate_audit,
        "prefix_audit": prefix_audit,
    }


def _read_exact(stream: Any, size: int, path: Path, description: str) -> bytes:
    data = stream.read(size)
    if len(data) != size:
        raise CheckFailure(f"{path}: truncated {description}")
    return data


def _finite_packed(data: bytes, code: str) -> tuple[int, float]:
    count = 0
    maximum = 0.0
    for (value,) in struct.iter_unpack("=" + code, data):
        if not math.isfinite(value):
            raise CheckFailure("binary output contains a non-finite stored field")
        maximum = max(maximum, abs(value))
        count += 1
    return count, maximum


def _parse_binary_parameters(raw: bytes, path: Path) -> dict[str, dict[str, str]]:
    try:
        lines = raw.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise CheckFailure(f"{path}: binary parameter header is not UTF-8") from exc
    parameters: dict[str, dict[str, str]] = {}
    block: str | None = None
    saw_end = False
    for raw_line in lines:
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            name = line[1:-1]
            if name == "par_end":
                saw_end = True
                break
            _require(bool(name), f"{path}: empty binary parameter block")
            _require(name not in parameters,
                     f"{path}: duplicate binary parameter block <{name}>")
            block = name
            parameters[block] = {}
            continue
        _require(block is not None and "=" in line,
                 f"{path}: malformed binary parameter header")
        key, value = line.split("=", 1)
        key = key.strip()
        _require(bool(key), f"{path}: empty binary parameter name")
        _require(key not in parameters[block],
                 f"{path}: duplicate binary parameter {block}/{key}")
        parameters[block][key] = value.strip()
    _require(saw_end, f"{path}: binary parameter header lacks <par_end>")
    return parameters


def _binary_output_is_full(parameters: dict[str, str]) -> bool:
    try:
        gid = int(parameters.get("gid", "-1"))
    except ValueError as exc:
        raise CheckFailure("binary output has an invalid gid parameter") from exc
    return (gid < 0 and "region_center" not in parameters and
            not any(f"slice_x{axis}" in parameters for axis in (1, 2, 3)))


def _float_width_value(value: float, width: int) -> float:
    code = "f" if width == 4 else "d"
    return struct.unpack("=" + code, struct.pack("=" + code, value))[0]


def _geometry_matches(actual: float, expected: float, width: int) -> bool:
    rounded = _float_width_value(expected, width)
    if actual == rounded:
        return True
    # Coordinate construction can round once in RegionSize arithmetic before it is
    # serialized.  A small width-aware ULP allowance accepts that operation but not a
    # displaced or fabricated MeshBlock.
    tolerance = max(math.ulp(actual), math.ulp(rounded)) * 8
    if width == 4:
        tolerance = max(tolerance, max(abs(rounded), 2.0**-126) * 2.0**-23 * 8)
    return abs(actual - rounded) <= tolerance


def _validate_binary_topology(
        path: Path, parameters: dict[str, dict[str, str]],
        output_parameters: dict[str, str], records: list[dict[str, Any]],
        location_size: int) -> dict[str, Any]:
    try:
        mesh = parameters["mesh"]
        meshblock = parameters["meshblock"]
        num_levels = int(parameters["mesh_refinement"]["num_levels"])
        mesh_cells = tuple(int(mesh[f"nx{axis}"]) for axis in (1, 2, 3))
        block_cells = tuple(int(meshblock[f"nx{axis}"]) for axis in (1, 2, 3))
        minima = tuple(float(mesh[f"x{axis}min"]) for axis in (1, 2, 3))
        maxima = tuple(float(mesh[f"x{axis}max"]) for axis in (1, 2, 3))
        nghost = int(mesh["nghost"])
    except (KeyError, ValueError) as exc:
        raise CheckFailure(f"{path}: invalid mesh metadata in binary header") from exc
    _require(all(value > 0 for value in mesh_cells + block_cells) and nghost >= 0 and
             num_levels > 0,
             f"{path}: nonpositive mesh dimensions or negative nghost")
    _require(all(total % block == 0
                 for total, block in zip(mesh_cells, block_cells)),
             f"{path}: mesh dimensions are not integral MeshBlocks")
    root_blocks = tuple(total // block
                        for total, block in zip(mesh_cells, block_cells))
    _require(all(math.isfinite(left) and math.isfinite(right) and right > left
                 for left, right in zip(minima, maxima)),
             f"{path}: invalid root-domain geometry")

    try:
        include_ghost = output_parameters.get("ghost_zones", "false").strip().lower()
        _require(include_ghost in ("false", "0", "true", "1"),
                 f"{path}: invalid ghost_zones parameter")
        use_ghost = include_ghost in ("true", "1")
        region_slice_axis = int(output_parameters.get("region_slice_axis", "0"))
    except ValueError as exc:
        raise CheckFailure(f"{path}: invalid binary slice parameters") from exc
    fixed_slices = {axis for axis in (1, 2, 3)
                    if f"slice_x{axis}" in output_parameters}
    fixed_slice_coordinates: dict[int, float] = {}
    for axis in fixed_slices:
        try:
            coordinate = float(output_parameters[f"slice_x{axis}"])
        except ValueError as exc:
            raise CheckFailure(
                f"{path}: invalid fixed slice_x{axis} coordinate") from exc
        _require(math.isfinite(coordinate),
                 f"{path}: non-finite fixed slice_x{axis} coordinate")
        fixed_slice_coordinates[axis] = coordinate
    if region_slice_axis:
        _require(region_slice_axis in (1, 2, 3) and not fixed_slices,
                 f"{path}: invalid moving-region slice axis")
        fixed_slices.add(region_slice_axis)

    topology: set[tuple[int, int, int, int]] = set()
    active_dimensions = tuple(total > 1 for total in mesh_cells)
    max_level = 0
    for block_number, record in enumerate(records):
        indices = record["indices"]
        logical = record["logical"]
        geometry = record["geometry"]
        lx = logical[:3]
        level = logical[3]
        _require(0 <= level < num_levels and all(value >= 0 for value in lx),
                 f"{path}: MeshBlock {block_number} has an invalid logical level/location")
        for axis, value in enumerate(lx):
            upper = (root_blocks[axis] << level
                     if active_dimensions[axis] else root_blocks[axis])
            _require(value < upper,
                     f"{path}: MeshBlock {block_number} logical coordinate is out of range")
        identity = (level, lx[0], lx[1], lx[2])
        _require(identity not in topology,
                 f"{path}: duplicate MeshBlock logical location {identity}")
        topology.add(identity)
        max_level = max(max_level, level)

        extents = tuple(indices[2 * axis + 1] - indices[2 * axis] + 1
                        for axis in range(3))
        expected_extents = tuple(
            1 if (not active_dimensions[axis] or axis + 1 in fixed_slices) else
            block_cells[axis] + (2 * nghost if use_ghost else 0)
            for axis in range(3)
        )
        _require(extents == expected_extents,
                 f"{path}: MeshBlock {block_number} index extents {extents} differ "
                 f"from header/output contract {expected_extents}")
        for axis in range(3):
            subdivisions = (root_blocks[axis] << level
                            if active_dimensions[axis] else root_blocks[axis])
            width = (maxima[axis] - minima[axis]) / subdivisions
            expected_left = minima[axis] + lx[axis] * width
            expected_right = minima[axis] + (lx[axis] + 1) * width
            _require(_geometry_matches(geometry[2 * axis], expected_left, location_size) and
                     _geometry_matches(geometry[2 * axis + 1], expected_right,
                                       location_size),
                     f"{path}: MeshBlock {block_number} geometry disagrees with its "
                     "logical location and root domain")
            start = nghost if active_dimensions[axis] else 0
            if axis + 1 in fixed_slice_coordinates:
                coordinate = fixed_slice_coordinates[axis + 1]
                _require(expected_left <= coordinate < expected_right,
                         f"{path}: MeshBlock {block_number} does not contain its "
                         f"fixed x{axis + 1} slice")
                relative = ((coordinate - expected_left) /
                            (expected_right - expected_left))
                selected = int(relative * block_cells[axis]) + start
                expected_pair = (selected, selected)
            elif axis + 1 == region_slice_axis:
                # The moving coordinate is independently recovered from the history at
                # this output cycle, then its exact CellCenterIndex is checked below.
                expected_pair = None
            else:
                extent = (block_cells[axis] + 2 * nghost
                          if use_ghost and active_dimensions[axis]
                          else block_cells[axis] if active_dimensions[axis] else 1)
                expected_pair = ((0 if use_ghost and active_dimensions[axis]
                                  else start),
                                 (0 if use_ghost and active_dimensions[axis]
                                  else start) + extent - 1)
            if expected_pair is not None:
                actual_pair = (indices[2 * axis], indices[2 * axis + 1])
                _require(actual_pair == expected_pair,
                         f"{path}: MeshBlock {block_number} x{axis + 1} indices "
                         f"{actual_pair} differ from exact output indices "
                         f"{expected_pair}")

    # No leaf may coexist with one of its ancestors.  Together with unique logical
    # locations this proves pairwise non-overlap without expanding the finest lattice.
    for level, lx1, lx2, lx3 in topology:
        for ancestor_level in range(level):
            shift = level - ancestor_level
            ancestor_coordinates = tuple(
                coordinate >> shift if active_dimensions[axis] else coordinate
                for axis, coordinate in enumerate((lx1, lx2, lx3)))
            ancestor = (ancestor_level, *ancestor_coordinates)
            _require(ancestor not in topology,
                     f"{path}: overlapping ancestor/descendant MeshBlocks")

    full_domain = _binary_output_is_full(output_parameters)
    coverage_units: int | None = None
    domain_units: int | None = None
    if full_domain:
        dimensions = sum(active_dimensions)
        coverage_units = sum(1 << (dimensions * (max_level - level))
                             for level, *_ in topology)
        domain_units = math.prod(root_blocks) * (1 << (dimensions * max_level))
        _require(coverage_units == domain_units,
                 f"{path}: full-domain AMR leaf coverage is incomplete "
                 f"({coverage_units} != {domain_units})")
    return {
        "scope": "full_domain" if full_domain else "selected_region_or_slice",
        "root_meshblocks": list(root_blocks),
        "root_domain": [[left, right] for left, right in zip(minima, maxima)],
        "meshblock_cells": list(block_cells),
        "active_dimensions": list(active_dimensions),
        "nghost": nghost,
        "ghost_zones": use_ghost,
        "maximum_physical_level": max_level,
        "configured_physical_levels": num_levels,
        "logical_locations_unique": True,
        "ancestor_descendant_overlap": False,
        "coverage_units": coverage_units,
        "domain_units": domain_units,
        "complete_leaf_coverage": full_domain,
    }


def audit_binary(path: Path, source_cycle: int, endpoint_cycle: int,
                 source_time: float, endpoint_time: float,
                 output_parameters: dict[str, str] | None = None,
                 endpoint_metadata: RestartMetadata | None = None,
                 expected_header_parameters: dict[str, dict[str, str]] | None = None,
                 output_block: str | None = None,
                 expected_variables: list[str] | None = None) -> dict[str, Any]:
    checked_path, raw, closed_signature = _open_regular_nofollow(path)
    digest = hashlib.sha256()
    field_count = 0
    max_abs = 0.0
    variable_max_abs: dict[str, float] = {}
    block_count = 0
    records: list[dict[str, Any]] = []
    try:
        with raw:
            exempt = {(os.getpid(), raw.fileno())}
            _assert_closed(checked_path, closed_signature, exempt)
            class HashReader:
                def read(self, size: int = -1) -> bytes:
                    data = raw.read(size)
                    digest.update(data)
                    return data

                def readline(self, size: int = -1) -> bytes:
                    data = raw.readline(size)
                    digest.update(data)
                    return data

                def tell(self) -> int:
                    return raw.tell()

            stream = HashReader()
            signature = stream.readline().decode("ascii")
            _require(signature == "Athena binary output version=1.1\n",
                     f"{path}: unsupported Athena binary signature")
            preheader_line = stream.readline()
            _require(preheader_line.startswith(b"  size of preheader="),
                     f"{path}: malformed preheader count")
            preheader_count = int(preheader_line.split(b"=", 1)[1])
            _require(preheader_count == 5, f"{path}: unsupported preheader count")
            preheader: dict[str, str] = {}
            for _ in range(preheader_count - 1):
                line = stream.readline().decode("ascii")
                _require("=" in line, f"{path}: malformed binary preheader")
                key, value = line.split("=", 1)
                preheader[key.strip()] = value.strip()
            number_line = stream.readline()
            _require(number_line.startswith(b"  number of variables="),
                     f"{path}: malformed variable count")
            variable_count = int(number_line.split(b"=", 1)[1])
            variable_line = stream.readline().decode("ascii")
            _require(":" in variable_line, f"{path}: malformed variable list")
            variables = variable_line.split(":", 1)[1].split()
            _require(variable_count > 0 and len(variables) == variable_count and
                     len(set(variables)) == variable_count,
                     f"{path}: invalid binary variable list")
            header_line = stream.readline()
            _require(header_line.startswith(b"  header offset="),
                     f"{path}: malformed parameter-header size")
            header_size = int(header_line.split(b"=", 1)[1])
            _require(0 <= header_size <= 16 * 1024 * 1024,
                     f"{path}: unreasonable parameter-header size")
            parameter_raw = _read_exact(
                stream, header_size, path, "binary parameter header")
            parameters = _parse_binary_parameters(parameter_raw, path)
            if expected_variables is None and output_parameters is not None:
                variable_group = output_parameters.get("variable", "").strip()
                if variable_group in ("mhd_divb", "mhd_gr_diagnostics"):
                    expected_variables = _expected_binary_variables(
                        variable_group, parameters)
                elif variable_group == "mhd_w_bcc" and expected_header_parameters:
                    expected_variables = _expected_binary_variables(
                        variable_group, expected_header_parameters)
            _require(expected_variables is not None and
                     variables == expected_variables,
                     f"{path}: binary variable list {variables} differs from exact "
                     f"planned semantics {expected_variables}")
            if expected_header_parameters is not None:
                for block in ("mesh", "meshblock", "mesh_refinement"):
                    _require(parameters.get(block) == expected_header_parameters.get(block),
                             f"{path}: binary {block} header differs from source restart")
                _require(output_block is not None and output_parameters is not None and
                         output_block in parameters,
                         f"{path}: planned output block is absent from binary header")
                mutable_output_keys = {
                    "file_number", "last_time", "last_write_cycle",
                }
                immutable_runtime = {
                    key: value for key, value in parameters[output_block].items()
                    if key not in mutable_output_keys
                }
                immutable_planned = {
                    key: value for key, value in output_parameters.items()
                    if key not in mutable_output_keys
                }
                _require(immutable_runtime == immutable_planned,
                         f"{path}: binary output header differs from planned parameters")
            cycle = int(preheader["cycle"])
            output_time = float(preheader["time"])
            location_size = int(preheader["size of location"])
            variable_size = int(preheader["size of variable"])
            _require(math.isfinite(output_time), f"{path}: non-finite output time")
            _require(source_cycle <= cycle <= endpoint_cycle,
                     f"{path}: cycle {cycle} lies outside segment")
            _require(source_time <= output_time <= endpoint_time,
                     f"{path}: time {output_time} lies outside segment")
            _require(location_size in (4, 8) and variable_size in (4, 8),
                     f"{path}: unsupported floating-point widths")
            location_code = "f" if location_size == 4 else "d"
            variable_code = "f" if variable_size == 4 else "d"
            file_size = closed_signature["size"]
            while stream.tell() < file_size:
                indices_raw = _read_exact(stream, 24, path, "MeshBlock index record")
                logical_raw = _read_exact(
                    stream, 16, path, "MeshBlock logical location")
                indices = struct.unpack("=6i", indices_raw)
                logical = struct.unpack("=4i", logical_raw)
                dimensions = (indices[1] - indices[0] + 1,
                              indices[3] - indices[2] + 1,
                              indices[5] - indices[4] + 1)
                _require(all(0 < extent <= 1_000_000 for extent in dimensions),
                         f"{path}: invalid MeshBlock dimensions {dimensions}")
                cells = math.prod(dimensions)
                coordinates = _read_exact(stream, 6 * location_size, path,
                                          "MeshBlock coordinate record")
                _finite_packed(coordinates, location_code)
                geometry = struct.unpack("=6" + location_code, coordinates)
                records.append({
                    "indices": indices,
                    "logical": logical,
                    "geometry": geometry,
                    "topology_bytes": indices_raw + logical_raw + coordinates,
                })
                bytes_per_variable = cells * variable_size
                data_size = variable_count * bytes_per_variable
                _require(data_size <= file_size - stream.tell(),
                         f"{path}: MeshBlock field payload exceeds file")
                for variable in variables:
                    remaining = bytes_per_variable
                    current_variable_max = variable_max_abs.get(variable, 0.0)
                    while remaining:
                        amount = min(8 * 1024 * 1024, remaining)
                        amount -= amount % variable_size
                        data = _read_exact(
                            stream, amount, path,
                            f"MeshBlock {variable} field payload")
                        count, current_max = _finite_packed(data, variable_code)
                        field_count += count
                        max_abs = max(max_abs, current_max)
                        current_variable_max = max(current_variable_max, current_max)
                        remaining -= amount
                    variable_max_abs[variable] = current_variable_max
                block_count += 1
            _require(stream.tell() == file_size and block_count > 0 and field_count > 0,
                     f"{path}: invalid or empty binary payload")
            _assert_stream_signature(raw, checked_path, closed_signature,
                                     "during binary audit")
            _assert_path_signature(checked_path, closed_signature,
                                   "during binary audit")
            _assert_closed(checked_path, closed_signature, exempt)
    except (KeyError, ValueError, UnicodeDecodeError, struct.error) as exc:
        if isinstance(exc, CheckFailure):
            raise
        raise CheckFailure(f"{path}: malformed Athena binary output: {exc}") from exc
    _assert_path_signature(checked_path, closed_signature, "after binary audit")
    _assert_closed(checked_path, closed_signature)
    selected_parameters = output_parameters
    if selected_parameters is None:
        binary_blocks = [value for key, value in parameters.items()
                         if re.fullmatch(r"output\d+", key) and
                         value.get("file_type", "").strip() == "bin"]
        _require(len(binary_blocks) == 1,
                 f"{path}: output_parameters are required for a multi-output header")
        selected_parameters = binary_blocks[0]
    topology = _validate_binary_topology(
        path, parameters, selected_parameters, records, location_size)
    topology_digest = hashlib.sha256()
    for record in records:
        topology_digest.update(record["topology_bytes"])
    logical_locations = tuple(record["logical"] for record in records)
    endpoint_topology_matches: bool | None = None
    if (endpoint_metadata is not None and topology["scope"] == "full_domain" and
            _close_in_ulps(output_time, endpoint_metadata.time, 0) and
            cycle == endpoint_metadata.cycle):
        expected_locations = tuple(
            (location.lx1, location.lx2, location.lx3,
             location.level - endpoint_metadata.root_level)
            for location in endpoint_metadata.locations
        )
        endpoint_topology_matches = logical_locations == expected_locations
        _require(block_count == endpoint_metadata.nmb_total and
                 endpoint_topology_matches,
                 f"{path}: endpoint binary topology/count disagrees with restart")
    return {
        "path": str(path),
        "sha256": digest.hexdigest(),
        **closed_signature,
        "closure_check": "linux_proc_fd",
        "time": output_time,
        "cycle": cycle,
        "variables": variables,
        "meshblocks": block_count,
        "topology_sha256": topology_digest.hexdigest(),
        "topology": topology,
        "endpoint_topology_matches_restart": endpoint_topology_matches,
        "_logical_locations": logical_locations,
        "_topology_records": tuple(
            (record["indices"], record["logical"], record["geometry"])
            for record in records
        ),
        "stored_field_count": field_count,
        "stored_fields_finite": True,
        "max_abs_stored_field": max_abs,
        "variable_max_abs": variable_max_abs,
    }


def _binary_pair_signature(parameters: dict[str, str]) -> str:
    ignored = {"variable", "id", "file_number", "last_time", "last_write_cycle"}
    return _canonical_sha256({
        key: value for key, value in parameters.items() if key not in ignored
    })


def audit_binary_pairs(plan: dict[str, Any],
                       audited_by_block: dict[str, list[dict[str, Any]]]
                       ) -> list[dict[str, Any]]:
    """Require each primitive-state stream to have an exact GR topology twin."""

    paired_variables = {"mhd_w_bcc", "mhd_gr_diagnostics"}
    groups: dict[str, dict[str, str]] = {}
    for output in plan["outputs"]:
        if not output.get("inspect_binary"):
            continue
        variable = output["parameters"].get("variable", "").strip()
        if variable not in paired_variables:
            continue
        signature = _binary_pair_signature(output["parameters"])
        members = groups.setdefault(signature, {})
        _require(variable not in members,
                 f"duplicate paired binary stream for {variable}")
        members[variable] = output["block"]
    _require(groups, "no W/GR binary topology pair is declared")
    summaries: list[dict[str, Any]] = []
    for signature, members in sorted(groups.items()):
        _require(set(members) == paired_variables,
                 "each W/GR scope and cadence must declare both binary streams")
        left_block = members["mhd_w_bcc"]
        right_block = members["mhd_gr_diagnostics"]
        left = audited_by_block.get(left_block, [])
        right = audited_by_block.get(right_block, [])
        _require(len(left) == len(right) and bool(left),
                 f"paired streams {left_block}/{right_block} have unequal file counts")
        paired_rows: list[dict[str, Any]] = []
        for left_row, right_row in zip(left, right):
            _require(left_row["cycle"] == right_row["cycle"] and
                     left_row["time"] == right_row["time"],
                     f"paired streams {left_block}/{right_block} have mismatched "
                     "cycle/time")
            _require(left_row["meshblocks"] == right_row["meshblocks"] and
                     left_row["topology_sha256"] == right_row["topology_sha256"],
                     f"paired streams {left_block}/{right_block} have mismatched "
                     "MeshBlock indices/logical locations/geometry")
            paired_rows.append({
                "cycle": left_row["cycle"], "time": left_row["time"],
                "meshblocks": left_row["meshblocks"],
                "topology_sha256": left_row["topology_sha256"],
            })
        summaries.append({
            "signature": signature,
            "w_block": left_block,
            "gr_block": right_block,
            "files": paired_rows,
        })
    return summaries


def _bbh_history_centers(history: dict[str, Any]) -> dict[int, dict[str, tuple[float, ...]]]:
    required = {
        "bh1_x", "bh1_y", "bh1_z", "bh2_x", "bh2_y", "bh2_z",
        "bh1_mass", "bh2_mass",
    }
    _require(required.issubset(history.get("columns", [])),
             "moving BBH output requires a history file with hole positions and masses")
    rows = history.get("_row_values")
    _require(isinstance(rows, list), "BBH history row evidence is unavailable")
    first_cycle = _integer(history.get("implicit_cycle_min"),
                           "history.implicit_cycle_min")
    centers: dict[int, dict[str, tuple[float, ...]]] = {}
    for offset, row in enumerate(rows):
        _require(isinstance(row, dict), "invalid BBH history row")
        bh1 = tuple(_number(row[f"bh1_{axis}"], f"bh1_{axis}")
                    for axis in ("x", "y", "z"))
        bh2 = tuple(_number(row[f"bh2_{axis}"], f"bh2_{axis}")
                    for axis in ("x", "y", "z"))
        mass1 = _number(row["bh1_mass"], "bh1_mass")
        mass2 = _number(row["bh2_mass"], "bh2_mass")
        total = mass1 + mass2
        _require(mass1 >= 0.0 and mass2 >= 0.0 and total > 0.0,
                 "BBH history contains invalid component masses")
        com = tuple((mass1 * bh1[axis] + mass2 * bh2[axis]) / total
                    for axis in range(3))
        centers[first_cycle + offset] = {"bh1": bh1, "bh2": bh2, "bbh_com": com}
    return centers


def audit_selected_binary_coverage(
        audited: dict[str, Any], output_parameters: dict[str, str],
        center: dict[str, tuple[float, ...]] | None,
        full_reference: dict[str, Any] | None = None) -> dict[str, Any]:
    """Prove selection by replaying it over a same-cycle full leaf topology."""

    if audited["topology"]["scope"] == "full_domain":
        return {"scope": "full_domain", "complete": True}
    try:
        gid = int(output_parameters.get("gid", "-1"))
    except ValueError as exc:
        raise CheckFailure("selected binary output has an invalid gid") from exc
    _require(gid < 0,
             "gid-selected binary topology cannot be proven complete by this qualifier")
    records = audited.get("_topology_records")
    _require(isinstance(records, tuple) and bool(records),
             "selected binary topology records are unavailable")
    _require(isinstance(full_reference, dict) and
             full_reference.get("topology", {}).get("scope") == "full_domain" and
             full_reference.get("topology", {}).get(
                 "complete_leaf_coverage") is True and
             full_reference.get("cycle") == audited.get("cycle") and
             full_reference.get("time") == audited.get("time"),
             "selected binary lacks a complete same-cycle full-domain topology reference")
    full_records = full_reference.get("_topology_records")
    _require(isinstance(full_records, tuple) and bool(full_records),
             "same-cycle full-domain topology records are unavailable")

    region_name = output_parameters.get("region_center")
    root_domain = audited["topology"].get("root_domain")
    _require(isinstance(root_domain, list) and len(root_domain) == 3,
             "selected binary root-domain evidence is unavailable")
    requested_bounds = tuple(tuple(float(value) for value in pair)
                             for pair in root_domain)
    region_center: tuple[float, ...] | None = None
    if region_name is not None:
        _require(center is not None and region_name in center,
                 "selected binary coverage requires a supported history-bound "
                 "moving region")
        region_center = center[region_name]
        try:
            default_half = float(output_parameters.get("region_half_width", "0"))
            half_widths = tuple(float(output_parameters.get(
                f"region_half_width{axis}", str(default_half)))
                for axis in (1, 2, 3))
        except ValueError as exc:
            raise CheckFailure(
                "moving output region has invalid half widths") from exc
        _require(all(math.isfinite(value) and value > 0.0
                     for value in half_widths),
                 "moving output region half widths must be positive and finite")
        requested_bounds = tuple(
            (region_center[axis] - half_widths[axis],
             region_center[axis] + half_widths[axis])
            for axis in range(3))

    slice_coordinates: dict[int, float] = {}
    for axis in range(3):
        key = f"slice_x{axis + 1}"
        if key in output_parameters:
            try:
                slice_coordinates[axis] = float(output_parameters[key])
            except ValueError as exc:
                raise CheckFailure(f"invalid {key} coordinate") from exc
    try:
        moving_axis = int(output_parameters.get("region_slice_axis", "0"))
        moving_offset = float(output_parameters.get("region_slice_offset", "0"))
    except ValueError as exc:
        raise CheckFailure("moving output slice has invalid parameters") from exc
    if moving_axis:
        _require(region_center is not None and moving_axis in (1, 2, 3) and
                 not slice_coordinates and math.isfinite(moving_offset),
                 "moving output slice contract is invalid")
        slice_coordinates[moving_axis - 1] = (
            region_center[moving_axis - 1] + moving_offset)
    _require(all(math.isfinite(value) for value in slice_coordinates.values()),
             "selected output contains a non-finite slice coordinate")

    def selected_by_contract(record: tuple[Any, Any, Any]) -> bool:
        geometry = record[2]
        if region_name is not None and not all(
                geometry[2 * axis + 1] > requested_bounds[axis][0] and
                geometry[2 * axis] < requested_bounds[axis][1]
                for axis in range(3)):
            return False
        return all(geometry[2 * axis] <= coordinate <
                   geometry[2 * axis + 1]
                   for axis, coordinate in slice_coordinates.items())

    expected_records = tuple(record for record in full_records
                             if selected_by_contract(record))
    _require(expected_records,
             "selected output contract selects no full-domain MeshBlocks")
    expected_logical = tuple(record[1] for record in expected_records)
    actual_logical = tuple(record[1] for record in records)
    _require(actual_logical == expected_logical,
             "selected binary logical locations differ from the exact same-cycle "
             "full-domain selection")
    expected_by_logical = {record[1]: record for record in expected_records}
    topology = audited["topology"]
    block_cells = topology.get("meshblock_cells")
    active = topology.get("active_dimensions")
    nghost = topology.get("nghost")
    _require(isinstance(block_cells, list) and len(block_cells) == 3 and
             isinstance(active, list) and len(active) == 3 and
             isinstance(nghost, int) and nghost >= 0,
             "selected binary lacks index-replay metadata")
    for indices, logical, geometry in records:
        reference_geometry = expected_by_logical[logical][2]
        _require(geometry == reference_geometry,
                 "selected binary geometry differs from its full-domain reference")
        for axis, coordinate in slice_coordinates.items():
            relative = ((coordinate - geometry[2 * axis]) /
                        (geometry[2 * axis + 1] - geometry[2 * axis]))
            selected_index = int(relative * block_cells[axis]) + (
                nghost if active[axis] else 0)
            actual_pair = (indices[2 * axis], indices[2 * axis + 1])
            _require(actual_pair == (selected_index, selected_index),
                     f"selected binary x{axis + 1} indices {actual_pair} differ "
                     f"from exact CellCenterIndex {selected_index}")
    return {
        "scope": "selected_region_or_slice",
        "complete": True,
        "center_name": region_name,
        "center": list(region_center) if region_center is not None else None,
        "slice_coordinates": {
            str(axis + 1): value for axis, value in slice_coordinates.items()
        },
        "selected_meshblocks": len(records),
        "reference_path": full_reference.get("path"),
        "reference_sha256": full_reference.get("sha256"),
        "reference_meshblocks": len(full_records),
        "selection_replayed_exactly": True,
    }


def _verify_endpoint_audit(result: Any, endpoint: Path) -> dict[str, Any]:
    _require(isinstance(result, dict) and result.get("valid") is True,
             f"complete restart audit did not validate {endpoint}")
    stored = result.get("stored_reals")
    _require(isinstance(stored, dict) and stored.get("nonfinite_count") == 0 and
             stored.get("finite_count") == stored.get("count"),
             "endpoint restart complete audit did not prove all stored Reals finite")
    layout = result.get("layout")
    _require(isinstance(layout, dict) and
             layout.get("expected_file_size") == _signature(endpoint)[2],
             "endpoint restart complete audit did not prove exact file layout")
    _require(result.get("topology", {}).get("complete_leaf_coverage") is True,
             "restart complete audit did not prove complete AMR leaf topology")
    return result


def _planned_endpoint_restart_path(plan: dict[str, Any], state_dir: Path) -> Path:
    restarts = [output for output in plan["outputs"]
                if output["file_type"] == "rst"]
    _require(len(restarts) == 1, "plan must have exactly one restart stream")
    output = restarts[0]
    writes = output.get("expected_writes")
    _require(isinstance(writes, list) and writes,
             "restart stream must have at least one planned write")
    final_write = writes[-1]
    _require(final_write.get("cycle") == plan["expected"]["final_cycle"] and
             final_write.get("time") == plan["expected"]["final_time"],
             "last planned restart is not the exact segment endpoint")
    relative = output["relative_path_template"].format(
        file_number=final_write["file_number"])
    return _safe_output_path(state_dir, relative)


def _planned_event_log_path(plan: dict[str, Any], state_dir: Path) -> Path:
    logs = [output for output in plan["outputs"] if output["file_type"] == "log"]
    _require(len(logs) == 1, "plan must have exactly one authoritative event log")
    paths = logs[0].get("required_unnumbered_paths")
    _require(isinstance(paths, list) and len(paths) == 1,
             "event log output must bind exactly one required path")
    return _safe_output_path(state_dir, paths[0])


def _plan_source_sha(plan: dict[str, Any]) -> str:
    return str(plan["inputs"]["source_restart"]["sha256"])


def _require_immutable_ready(path: Path, suffix: str, label: str) -> None:
    absolute = path.absolute()
    info = absolute.lstat()
    _require(absolute.name.endswith(suffix) and stat.S_ISREG(info.st_mode) and
             not stat.S_ISLNK(info.st_mode) and
             not (stat.S_IMODE(info.st_mode) & 0o222),
             f"{label} must be an immutable regular {suffix} file")


def _same_full_binding(actual: Any, expected: dict[str, Any]) -> bool:
    return isinstance(actual, dict) and all(
        actual.get(name) == expected.get(name) for name in (
            "path", "device", "inode", "size", "mtime_ns", "ctime_ns", "sha256",
        )
    )


def audit_completion_record(
        path: Path, plan: dict[str, Any], state_dir: Path,
        launch: dict[str, Any], artifact_bindings: dict[str, dict[str, Any]],
        expected_return_code: int = 0) \
        -> dict[str, Any]:
    """Consume the sole marker that proves the launcher lifecycle is complete."""

    _require_immutable_ready(path, ".completion.ready", "completion record")
    record, raw = _load_json(path)
    _require(set(artifact_bindings) == {
        "plan", "launch_record", "run_log", "exit_status",
        "gpu_before", "gpu_after",
    }, "checker completion artifact binding set is incomplete")
    _require(isinstance(expected_return_code, int) and
             not isinstance(expected_return_code, bool),
             "expected completion return code must be an integer")
    _require(record.get("schema") == SCHEMA and
             record.get("kind") == "athenak_segment_completion" and
             record.get("status") == "ready" and
             record.get("return_code") == expected_return_code and
             record.get("state_dir") == str(state_dir.resolve(strict=True)) and
             record.get("world_size") == plan["policy"]["ranks"],
             "completion record schema/status/exit/state/world binding is invalid")
    artifacts = record.get("artifacts")
    _require(isinstance(artifacts, dict) and
             set(artifacts) == set(artifact_bindings),
             "completion record artifact set is not exact")
    for name, binding in artifact_bindings.items():
        fresh = _hash_record(Path(binding["path"]))
        _require(_same_full_binding(binding, fresh),
                 f"independent {name} binding differs from its current closed bytes")
        _require(_same_full_binding(artifacts[name], binding),
                 f"completion record {name} binding differs from independently "
                 "audited evidence")
    publication = record.get("publication_contract")
    _require(publication == {
        "unique_lifecycle_complete_marker": True,
        "published_after_closed_artifacts": [
            "launch_record", "run_log", "exit_status", "gpu_before", "gpu_after",
        ],
    }, "completion publication contract is not exact")
    quiescence = record.get("quiescence")
    _require(isinstance(quiescence, dict) and
             quiescence.get("gpu_compute_contexts_empty") is True and
             quiescence.get("all_original_identities_gone") is True and
             quiescence.get("managed_process_group") ==
             launch["managed_process_group"] and
             quiescence.get("managed_process_group_gone") is True,
             "completion record does not prove process/GPU quiescence")
    identities = quiescence.get("process_identities")
    _require(isinstance(identities, list),
             "completion quiescence lacks process identities")
    expected_identities = [{
        "role": "mpirun",
        "pid": launch["mpirun_pid"],
        "recorded_start_time_ticks": launch["mpirun_start_time_ticks"],
    }]
    expected_identities.extend({
        "role": f"rank:{rank['global_rank']}",
        "pid": rank["pid"],
        "recorded_start_time_ticks": rank["start_time_ticks"],
    } for rank in launch["ranks"])
    _require(len(identities) == len(expected_identities),
             "completion quiescence process count differs from launch proof")
    for actual, expected_identity in zip(identities, expected_identities):
        _require(isinstance(actual, dict) and
                 all(actual.get(name) == value
                     for name, value in expected_identity.items()) and
                 actual.get("state") in ("disappeared", "pid_reused") and
                 actual.get("original_identity_gone") is True and
                 ((actual.get("state") == "disappeared" and
                   actual.get("observed_start_time_ticks") is None) or
                  (actual.get("state") == "pid_reused" and
                   isinstance(actual.get("observed_start_time_ticks"), int) and
                   actual["observed_start_time_ticks"] !=
                   actual["recorded_start_time_ticks"])),
                 "completion process quiescence differs from launch identities")
    transport = record.get("input_transport")
    launch_transport = launch.get("input_transport")
    _require(isinstance(transport, dict) and
             transport.get("kind") == INPUT_TRANSPORT_KIND and
             transport.get("at_launch") == launch_transport,
             "completion input transport does not bind the launch-time held FDs")
    at_completion = transport.get("at_completion")
    closure = transport.get("closure")
    _require(isinstance(at_completion, dict) and
             at_completion.get("kind") == INPUT_TRANSPORT_KIND and
             at_completion.get("holder_pid") ==
             launch_transport.get("holder_pid") and
             at_completion.get("holder_start_time_ticks") ==
             launch_transport.get("holder_start_time_ticks") and
             at_completion.get("roles") == launch_transport.get("roles"),
             "held input descriptors changed between launch and completion")
    _require(isinstance(closure, dict) and
             closure.get("all_holder_fds_closed") is True and
             set(closure.get("roles", {})) == {"source_restart", "trajectory"},
             "completion does not prove both input holder descriptors closed")
    for role, descriptor in (("source_restart", SOURCE_RESTART_FD),
                             ("trajectory", TRAJECTORY_FD)):
        _require(closure["roles"][role] == {"fd": descriptor, "closed": True},
                 f"completion holder closure for {role} is invalid")
    completion_directories = _audit_directory_transport_record(
        record.get("directory_transport"),
        plan.get("launch_contract", {}).get("directory_transport", {}),
        launch_transport["holder_pid"],
        launch_transport["holder_start_time_ticks"])
    _require(completion_directories == launch.get("directory_transport"),
             "completion directory-holder proof changed since launch")
    completion_executables = _audit_executable_transport_record(
        record.get("executable_transport"),
        plan.get("launch_contract", {}).get("executable_transport", {}),
        launch_transport["holder_pid"],
        launch_transport["holder_start_time_ticks"], plan)
    _require(completion_executables == launch.get("executable_transport"),
             "completion executable-holder proof changed since launch")
    marker = _signature(path)
    _require(all(marker[3] >= binding["mtime_ns"] and
                 marker[4] >= binding["ctime_ns"]
                 for binding in artifact_bindings.values()),
             "completion record is not newer than every closed lifecycle artifact")
    evidence_dir = path.parent.resolve(strict=True)
    expected_nodes = {
        path.resolve(strict=True),
        *(Path(binding["path"]).resolve(strict=True)
          for binding in artifact_bindings.values()
          if Path(binding["path"]).resolve(strict=True).parent == evidence_dir),
    }
    actual_nodes: set[Path] = set()
    for entry in os.scandir(evidence_dir):
        node = entry.stat(follow_symlinks=False)
        _require(stat.S_ISREG(node.st_mode) and not stat.S_ISLNK(node.st_mode),
                 f"evidence directory contains a non-regular node: {entry.path}")
        actual_nodes.add(Path(entry.path).resolve(strict=True))
    _require(actual_nodes == expected_nodes,
             "evidence directory does not contain exactly the lifecycle artifacts "
             "plus the sole completion marker")
    return {
        "sha256": hashlib.sha256(raw).hexdigest(),
        "artifacts": artifacts,
        "quiescence": quiescence,
        "input_transport": transport,
        "directory_transport": completion_directories,
        "executable_transport": completion_executables,
        "publication_contract": publication,
        "return_code": expected_return_code,
    }


def check_segment(args: argparse.Namespace) -> dict[str, Any]:
    output = args.output
    _reject_output_ancestors(output)
    _require(output.name.endswith(".pass.ready"),
             "--output must end in .pass.ready")
    _require(not output.exists() and not output.is_symlink(),
             f"refusing to overwrite report: {output}")
    plan, plan_raw = _load_json(args.plan)
    expected = validate_plan(plan)
    _require(args.plan.resolve(strict=True) ==
             Path(plan["launch_contract"]["plan_path"]).resolve(strict=True),
             "--plan differs from the exact plan-bound evidence path")
    state_absolute = args.state_dir.absolute()
    state_stat = state_absolute.lstat()
    _require(stat.S_ISDIR(state_stat.st_mode) and
             not stat.S_ISLNK(state_stat.st_mode),
             "--state-dir must be a real, non-symlink directory")
    state_resolved = state_absolute.resolve(strict=True)
    _require(plan.get("launch_contract", {}).get("state_dir") ==
             str(state_resolved),
             "--state-dir differs from the immutable plan")
    planned_evidence = plan.get("launch_contract", {}).get("evidence", {})
    for argument_name, path in (
        ("launch_record", args.launch_record),
        ("completion_record", args.completion_record),
        ("run_log", args.run_log),
        ("exit_status", args.exit_status),
        ("gpu_before", args.gpu_before),
        ("gpu_after", args.gpu_after),
    ):
        _require(isinstance(planned_evidence.get(argument_name), str) and
                 path.resolve(strict=True) ==
                 Path(planned_evidence[argument_name]).resolve(strict=True),
                 f"--{argument_name.replace('_', '-')} differs from the exact "
                 "plan-bound evidence path")
    planned_endpoint = _planned_endpoint_restart_path(plan, state_resolved)
    planned_event_log = _planned_event_log_path(plan, state_resolved)
    _require(args.endpoint_restart.resolve(strict=True) == planned_endpoint,
             "--endpoint-restart is not the independently replayed final restart")
    _require(args.event_log.resolve(strict=True) == planned_event_log,
             "--event-log is not the unique plan-required event log")
    source_path = Path(plan["inputs"]["source_restart"]["path"])
    initial_ready = [args.plan, args.launch_record, source_path,
                     args.endpoint_restart, args.run_log,
                     args.event_log, args.exit_status, args.gpu_before, args.gpu_after,
                     args.completion_record]
    _check_ready(initial_ready)
    plan_binding = _hash_record(args.plan)
    _require(plan_binding["sha256"] == hashlib.sha256(plan_raw).hexdigest(),
             "plan changed between parsing and closed-file binding")
    exit_binding = _hash_record(args.exit_status)
    run_binding = _hash_record(args.run_log)
    event_binding = _hash_record(args.event_log)
    gpu_before_binding = _hash_record(args.gpu_before)
    gpu_after_binding = _hash_record(args.gpu_after)
    launch_binding = _hash_record(args.launch_record)
    completion_binding = _hash_record(args.completion_record)
    _require(_parse_exit_status(args.exit_status) == 0,
             "process exit status must be zero")
    launch = audit_launch_record(
        args.launch_record, plan, args.plan, args.state_dir, gpu_before_binding)

    repository_binding = _verify_planned_repository(
        plan["inputs"]["repo"], plan["tools"]["git"])
    binary_binding = _verify_planned_file(plan["inputs"]["binary"], "AthenaK binary")
    trajectory_binding = _verify_planned_file(
        plan["inputs"]["trajectory"], "trajectory")
    planned_tool_bindings = {
        name: _verify_planned_file(plan["tools"][name], f"planning tool {name}")
        for name in PLANNED_TOOL_NAMES
    }
    current_tool_bindings = {
        "segment_checker": _verify_current_tool(
            plan["tools"]["segment_checker"], Path(__file__), "segment checker"),
        "restart_auditor": _verify_current_tool(
            plan["tools"]["restart_auditor"],
            Path(sys.modules[audit_restart.__module__].__file__), "restart auditor"),
        "output_integrity": _verify_current_tool(
            plan["tools"]["output_integrity"],
            Path(sys.modules[stable_sha256.__module__].__file__), "output integrity"),
        "restart_metadata_reader": _verify_current_tool(
            plan["tools"]["restart_metadata_reader"],
            Path(sys.modules[read_restart_metadata.__module__].__file__),
            "restart metadata reader"),
        "restart_layout": _verify_current_tool(
            plan["tools"]["restart_layout"],
            SCRIPT_DIRECTORY / "compare_athenak_restart_fields.py",
            "restart layout tool"),
    }

    source_audit = _verify_endpoint_audit(audit_restart(source_path), source_path)
    source_binding = {
        "path": source_audit["path"],
        **source_audit["signature"],
        "sha256": source_audit["sha256"],
        "closure_check": source_audit["closure_check"],
    }
    source_metadata = read_restart_metadata(source_path)
    endpoint_audit = _verify_endpoint_audit(
        audit_restart(args.endpoint_restart), args.endpoint_restart)
    endpoint_binding = {
        "path": endpoint_audit["path"],
        **endpoint_audit["signature"],
        "sha256": endpoint_audit["sha256"],
        "closure_check": endpoint_audit["closure_check"],
    }
    endpoint_metadata = read_restart_metadata(args.endpoint_restart)
    source_sha = str(source_binding["sha256"])
    _require(source_sha == _plan_source_sha(plan),
             "source restart SHA-256 differs from the signed plan")
    _require(source_binding["size"] == plan["inputs"]["source_restart"]["size"],
             "source restart size differs from the signed plan")
    source_qualification = plan["source_qualification"]
    planned_source_audit = source_qualification["audit"]
    _require(planned_source_audit["sha256"] == source_audit["sha256"] and
             planned_source_audit["layout"] == source_audit["layout"] and
             planned_source_audit["stored_reals"] == source_audit["stored_reals"] and
             planned_source_audit["topology"] == source_audit["topology"],
             "live source audit differs from the plan's full source audit")
    source_mass_evidence = _verify_planned_file(
        source_qualification["source_baryon_mass"]["evidence"],
        "source baryon mass evidence")
    parent: dict[str, Any] | None = None
    parent_binding: dict[str, Any] | None = None
    parent_qualification_audit: dict[str, Any] | None = None
    if source_qualification["mode"] == "parent_segment_pass":
        parent_record = source_qualification["parent_segment_pass"]
        parent_path = Path(parent_record["path"])
        parent_stat = parent_path.lstat()
        _require(parent_path.name.endswith(".pass.ready") and
                 stat.S_ISREG(parent_stat.st_mode) and
                 not stat.S_ISLNK(parent_stat.st_mode) and
                 not (stat.S_IMODE(parent_stat.st_mode) & 0o222),
                 "parent segment pass must be immutable .pass.ready evidence")
        parent_binding = _verify_planned_file(
            parent_record, "parent segment pass")
        parent, _ = _load_json(parent_path)
        _require(parent.get("kind") == "athenak_segment_pass" and
                 parent.get("status") == "pass" and
                 parent.get("bindings", {}).get("endpoint_restart", {}).get("sha256") ==
                 source_binding["sha256"] and
                 parent.get("bindings", {}).get("endpoint_restart", {}).get("path") ==
                 str(source_path.absolute()) and
                 parent.get("bindings", {}).get("endpoint_restart", {}).get("size") ==
                 source_binding["size"] and
                 parent.get("bindings", {}).get("trajectory", {}).get("size") ==
                 plan["inputs"]["trajectory"]["size"] and
                 parent.get("bindings", {}).get("trajectory", {}).get("sha256") ==
                 plan["inputs"]["trajectory"]["sha256"] and
                 parent.get("endpoint_restart_audit", {}).get("stored_reals", {}).get(
                     "nonfinite_count") == 0 and
                 parent.get("endpoint_restart_audit", {}).get("layout") ==
                 source_audit["layout"] and
                 parent.get("endpoint_restart_audit", {}).get("stored_reals") ==
                 source_audit["stored_reals"] and
                 parent.get("expected", {}).get("final_cycle") ==
                 expected["source_cycle"] and
                 parent.get("expected", {}).get("final_time") ==
                 expected["source_time"],
                 "parent pass does not bind this fully audited source restart")
        parent_inventory = parent.get("output_inventory")
        _require(isinstance(parent_inventory, list),
                 "parent pass lacks an output inventory")
        parent_histories = [
            row for row in parent_inventory
            if isinstance(row, dict) and
            isinstance(row.get("history"), dict) and
            "baryon_m" in row["history"].get("columns", [])
        ]
        _require(len(parent_histories) == 1 and all(
            parent_histories[0].get(name) == source_mass_evidence.get(name)
            for name in ("path", "size", "sha256")),
            "source history evidence is not the unique parent-audited baryon history")
        parent_qualification_audit = audit_parent_qualification_provenance(
            parent, source_qualification, source_binding, source_mass_evidence,
            expected["source_cycle"], expected["source_time"], parent_path)
    source_baryon_mass = audit_source_baryon_evidence(
        source_qualification, expected["source_time"], parent)
    if parent is not None:
        _require(_number(
            parent_histories[0]["history"].get("column_last", {}).get("baryon_m"),
            "parent output-inventory endpoint baryon mass") ==
            source_baryon_mass and
            _number(parent.get("scientific_threshold_audit", {}).get(
                "baryon_mass", {}).get("last"),
                "parent scientific endpoint baryon mass") == source_baryon_mass,
            "parent audited baryon mass differs from independently parsed source "
            "history")
    _require(source_metadata.cycle == expected["source_cycle"],
             "source restart cycle differs from the plan")
    _require(_close_in_ulps(source_metadata.time, expected["source_time"],
                            expected["endpoint_time_ulp_tolerance"]),
             "source restart time differs from the plan")
    source_summary = plan.get("source")
    _require(isinstance(source_summary, dict), "plan source summary is missing")
    for key, actual in (
        ("cycle", source_metadata.cycle), ("time", source_metadata.time),
        ("real_bytes", source_metadata.real_bytes),
        ("nmb_total", source_metadata.nmb_total),
        ("event_counter_version", source_metadata.event_counter_version),
        ("restart_cache_contract_version",
         source_metadata.restart_cache_contract_version),
    ):
        _require(source_summary.get(key) == actual,
                 f"plan source summary {key} differs from source restart")
    _require(source_summary.get("parameters") == source_metadata.parameters,
             "plan source parameter snapshot differs from source restart")
    source_parameter_sha = _canonical_sha256(source_metadata.parameters)
    _require(source_summary.get("parameters_sha256") == source_parameter_sha,
             "plan source parameter hash is invalid")
    parameter_snapshots = plan.get("parameter_snapshots")
    _require(isinstance(parameter_snapshots, dict) and
             parameter_snapshots.get("source_sha256") == source_parameter_sha,
             "plan parameter_snapshots source hash is invalid")
    _require(endpoint_metadata.cycle == expected["final_cycle"],
             "endpoint restart cycle differs from the plan")
    _require(_close_in_ulps(endpoint_metadata.time, expected["final_time"],
                            expected["endpoint_time_ulp_tolerance"]),
             "endpoint restart time differs from the plan")
    source_contract = audit_restart_contract(
        source_metadata, "source", plan, expected)
    endpoint_contract = audit_restart_contract(
        endpoint_metadata, "endpoint", plan, expected)
    trajectory_tokens = [
        token for token in launch["athena_argv"]
        if token.startswith("problem/trajectory_file=")
    ]
    _require(len(trajectory_tokens) == 1,
             "launch evidence lacks one exact trajectory provenance override")
    runtime_trajectory = trajectory_tokens[0].split("=", 1)[1]
    source_trajectory = source_metadata.parameters.get(
        "problem", {}).get("trajectory_file")
    _require(isinstance(source_trajectory, str) and bool(source_trajectory),
             "source restart lacks serialized trajectory provenance")
    exact_rebindings = {
            "output3/dt": {
                "source": [source_metadata.parameters["output3"]["dt"]],
                "endpoint": repr(expected["root_dt"]),
            },
            "problem/trajectory_file": {
                "source": [source_trajectory], "endpoint": runtime_trajectory,
            },
        }
    restart_cadence_transition = _validate_restart_cadence_transition(
        plan, float(expected["root_dt"]))
    if restart_cadence_transition["kind"] == "tighten_v1":
        restart_block = restart_cadence_transition["block"]
        exact_rebindings[restart_cadence_transition["parameter"]] = {
            "source": [source_metadata.parameters[restart_block]["dt"]],
            "endpoint": repr(restart_cadence_transition["target_dt"]),
        }
    capacity_transition = _validate_capacity_transition(plan)
    if capacity_transition["kind"] == "increase_v1":
        exact_rebindings["mesh_refinement/max_nmb_per_rank"] = {
            "source": [str(capacity_transition["source_max_nmb_per_rank"])],
            "endpoint": str(capacity_transition["target_max_nmb_per_rank"]),
        }
    parameters = compare_parameters(
        source_metadata, endpoint_metadata, plan["policy"]["mutable_parameters"],
        exact_rebindings)

    run = audit_run_log(args.run_log, expected)
    events = audit_event_log(args.event_log, expected["source_cycle"],
                             expected["final_cycle"],
                             plan["policy"]["event_thresholds"],
                             plan["policy"]["event_absolute_thresholds"])
    gpus = audit_gpus(args.gpu_before, args.gpu_after, expected["ranks"],
                      expected["gpu_exit_memory_mib_max"],
                      plan["policy"].get("gpu_ecc"))
    capacity_transition = _validate_capacity_transition(plan)
    capacity_preflight = _gpu_capacity_preflight_evidence(
        gpus["before"], capacity_transition)
    gpus["capacity_preflight"] = {
        "capacity_transition": capacity_transition,
        "evidence": capacity_preflight,
    }
    completion = audit_completion_record(
        args.completion_record, plan, state_resolved, launch, {
            "plan": plan_binding,
            "launch_record": launch_binding,
            "run_log": run_binding,
            "exit_status": exit_binding,
            "gpu_before": gpu_before_binding,
            "gpu_after": gpu_after_binding,
        })
    _require(completion["sha256"] == completion_binding["sha256"],
             "completion marker changed between binding and semantic audit")

    inventory = expected_output_paths(
        plan, source_metadata, endpoint_metadata, args.state_dir, expected)
    replayed_restarts = [row for row in inventory if row["file_type"] == "rst"]
    _require(replayed_restarts and
             replayed_restarts[-1]["expected_write"]["cycle"] ==
             expected["final_cycle"] and
             replayed_restarts[-1]["expected_write"]["time"] ==
             expected["final_time"] and
             args.endpoint_restart.resolve(strict=True) ==
             replayed_restarts[-1]["path"],
             "--endpoint-restart is not the independently replayed final restart")
    _check_ready(row["path"] for row in inventory)
    state_tree = _assert_exact_state_tree(
        state_resolved, (row["path"] for row in inventory))
    replayed_schedules = _replay_all_output_writes(
        plan, source_metadata, endpoint_metadata, expected)
    for output in plan["outputs"]:
        _require(output["expected_writes"] ==
                 replayed_schedules[output["block"]],
                 f"{output['block']}: plan expected_writes differs from independent "
                 "cadence replay")
        endpoint_block = endpoint_metadata.parameters[output["block"]]
        actual_endpoint_state = {
            "file_number": int(endpoint_block.get("file_number", "0")),
            "last_time": float(endpoint_block["last_time"]),
            "last_write_cycle": int(endpoint_block["last_write_cycle"]),
        }
        _require(output["expected_endpoint_state"] == actual_endpoint_state,
                 f"{output['block']}: planned endpoint cadence state differs from restart")
    output_rows: list[dict[str, Any]] = []
    binary_by_block: dict[str, list[dict[str, Any]]] = {}
    divb_maxima: list[dict[str, Any]] = []
    baryon_histories: list[dict[str, Any]] = []
    bbh_histories: list[dict[str, Any]] = []
    for row in inventory:
        if row["inspect_binary"]:
            audited = audit_binary(
                row["path"], expected["source_cycle"], expected["final_cycle"],
                expected["source_time"], endpoint_metadata.time,
                row["parameters"], endpoint_metadata, source_metadata.parameters,
                row["block"], row["expected_binary_variables"])
            binary_by_block.setdefault(row["block"], []).append(audited)
            expected_write = row["expected_write"]
            _require(audited["cycle"] == expected_write["cycle"] and
                     audited["time"] == expected_write["time"],
                     f"{row['path']}: binary cycle/time differs from cadence replay")
            for variable, maximum in audited["variable_max_abs"].items():
                if variable.lower() == "divb":
                    divb_maxima.append({
                        "path": str(row["path"]),
                        "time": audited["time"],
                        "cycle": audited["cycle"],
                        "max_abs": maximum,
                    })
        elif row["file_type"] == "hst" or row["path"].suffix == ".hst":
            history_binding = _hash_record(row["path"])
            history = audit_history(
                row["path"], expected["source_cycle"], expected["final_cycle"],
                expected["source_time"], expected["root_dt"],
                endpoint_metadata.time,
                expected["endpoint_time_ulp_tolerance"])
            _assert_binding_unchanged(row["path"], history_binding)
            audited = {**history_binding, "history": history}
            if "baryon_m" in history["columns"]:
                baryon_histories.append({"path": str(row["path"]), **history})
            if {"bh1_x", "bh2_x", "bh1_mass", "bh2_mass"}.issubset(
                    history["columns"]):
                bbh_histories.append({"path": str(row["path"]), **history})
        elif row["file_type"] == "rst":
            restart_audit = _verify_endpoint_audit(
                audit_restart(row["path"]), row["path"])
            expected_write = row["expected_write"]
            metadata_record = restart_audit.get("metadata", {})
            _require(isinstance(expected_write, dict) and
                     metadata_record.get("cycle") == expected_write["cycle"] and
                     metadata_record.get("time") == expected_write["time"],
                     f"{row['path']}: restart cycle/time differs from cadence replay")
            audited = {
                "path": restart_audit["path"],
                **restart_audit["signature"],
                "sha256": restart_audit["sha256"],
                "closure_check": restart_audit["closure_check"],
                "restart_audit": restart_audit,
            }
        else:
            audited = _hash_record(row["path"])
        output_rows.append({
            "block": row["block"],
            "file_number": row["file_number"],
            "expected_write": row["expected_write"],
            **audited,
        })

    binary_pairs = audit_binary_pairs(plan, binary_by_block)
    _require(len(bbh_histories) == 1,
             "exactly one required history file must carry BBH trajectory columns")
    centers = _bbh_history_centers(bbh_histories[0])
    full_reference_by_cycle: dict[int, dict[str, Any]] = {}
    for rows in binary_by_block.values():
        for audited in rows:
            if audited["topology"]["scope"] != "full_domain":
                continue
            previous = full_reference_by_cycle.get(audited["cycle"])
            if previous is not None:
                previous_geometry = tuple(
                    (record[1], record[2])
                    for record in previous["_topology_records"])
                current_geometry = tuple(
                    (record[1], record[2])
                    for record in audited["_topology_records"])
                _require(previous["time"] == audited["time"] and
                         previous_geometry == current_geometry,
                         "same-cycle full-domain binary streams disagree on the AMR "
                         "leaf topology")
            else:
                full_reference_by_cycle[audited["cycle"]] = audited
    selection_coverage: list[dict[str, Any]] = []
    output_plan_by_block = {output["block"]: output for output in plan["outputs"]}
    for block, rows in binary_by_block.items():
        output_parameters = output_plan_by_block[block]["parameters"]
        for audited in rows:
            cycle = audited["cycle"]
            center = centers.get(cycle)
            coverage = audit_selected_binary_coverage(
                audited, output_parameters, center,
                full_reference_by_cycle.get(cycle))
            audited["selection_coverage"] = coverage
            selection_coverage.append({
                "path": audited["path"], "cycle": cycle, **coverage,
            })
    coverage_by_path = {row["path"]: row for row in selection_coverage}
    for audited in output_rows:
        if audited.get("path") in coverage_by_path:
            audited["selection_coverage"] = {
                key: value for key, value in
                coverage_by_path[audited["path"]].items() if key != "path"
            }

    divb_hard = _number(plan["policy"]["divb_max_abs"]["hard"],
                         "divb_max_abs.hard")
    _require(divb_maxima, "no inspected binary output contained the divB field")
    worst_divb = max(item["max_abs"] for item in divb_maxima)
    _require(worst_divb <= divb_hard,
             f"segment divB maximum {worst_divb} exceeds {divb_hard}")
    _require(len(baryon_histories) == 1,
             "exactly one required history file must contain baryon_m")
    baryon_history = baryon_histories[0]
    _require(baryon_history["endpoint_time_matches"],
             "baryon-mass history does not reach the exact segment endpoint")
    baryon_audit = audit_baryon_mass(
        baryon_history,
        plan["policy"]["baryon_mass_fractional_loss"]["hard_per_root_step"],
        expected["root_steps"],
        source_baryon_mass)
    baryon_policy = plan["policy"]["baryon_mass_fractional_loss"]
    baryon_advisory = audit_baryon_mass_advisory(
        baryon_history, source_baryon_mass, expected["source_cycle"],
        baryon_policy["yellow_per_root_step"],
        baryon_policy["yellow_per_48M"],
        baryon_policy["rolling_window_root_steps"])
    event_advisories = audit_event_ratio_advisories(
        events, plan["policy"]["yellow_event_thresholds"])
    floor_trends = audit_floor_rate_trends(events, parent)
    divb_yellow = _number(plan["policy"]["divb_max_abs"]["yellow"],
                           "divb_max_abs.yellow")
    worst_divb_file = max(divb_maxima, key=lambda item: item["max_abs"])
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
        "severity": ("yellow" if event_advisories["severity"] == "yellow" or
                     baryon_advisory["severity"] == "yellow" or
                     divb_advisory["severity"] == "yellow" else "green"),
        "event_ratios": event_advisories,
        "baryon_mass": baryon_advisory,
        "divb": divb_advisory,
        "floor_rates": floor_trends,
        "pass_fail_effect": "none_yellow_advisories_are_nonfatal",
    }
    scientific = {
        "nonfinite_count": 0,
        "divb": {"hard_max_abs": divb_hard, "observed_max_abs": worst_divb,
                 "files": divb_maxima},
        "baryon_mass": {
            "path": baryon_history["path"],
            **baryon_audit,
        },
    }

    # Row-level scratch arrays are retained until every hard and advisory proof has
    # finished; only compact, concrete cycle/rate summaries enter the pass report.
    events.pop("_rows", None)
    for rows in binary_by_block.values():
        for audited in rows:
            audited.pop("_logical_locations", None)
            audited.pop("_topology_records", None)
    for audited in output_rows:
        audited.pop("_logical_locations", None)
        audited.pop("_topology_records", None)
        history = audited.get("history")
        if isinstance(history, dict):
            history.pop("_row_values", None)
    for history in (*baryon_histories, *bbh_histories):
        history.pop("_row_values", None)

    script_path = Path(__file__).resolve()
    audit_tool = Path(sys.modules[audit_restart.__module__].__file__).resolve()
    integrity_tool = Path(sys.modules[stable_sha256.__module__].__file__).resolve()
    metadata_tool = SCRIPT_DIRECTORY / "read_athenak_restart_metadata.py"
    layout_tool = SCRIPT_DIRECTORY / "compare_athenak_restart_fields.py"
    for path, binding in (
        (args.plan, plan_binding), (source_path, source_binding),
        (args.launch_record, launch_binding),
        (args.completion_record, completion_binding),
        (args.endpoint_restart, endpoint_binding),
        (args.run_log, run_binding), (args.event_log, event_binding),
        (args.exit_status, exit_binding), (args.gpu_before, gpu_before_binding),
        (args.gpu_after, gpu_after_binding),
        (Path(plan["inputs"]["binary"]["path"]), binary_binding),
        (Path(plan["inputs"]["trajectory"]["path"]), trajectory_binding),
        (Path(source_qualification["source_baryon_mass"]["evidence"]["path"]),
         source_mass_evidence),
    ):
        _assert_binding_unchanged(path, binding)
    for name, binding in planned_tool_bindings.items():
        _assert_binding_unchanged(Path(plan["tools"][name]["path"]), binding)
    for name, binding in current_tool_bindings.items():
        _assert_binding_unchanged(Path(plan["tools"][name]["path"]), binding)
    _require(_verify_planned_repository(
        plan["inputs"]["repo"], plan["tools"]["git"]) ==
             repository_binding,
             "planned repository changed during segment qualification")
    for row in output_rows:
        _assert_binding_unchanged(Path(row["path"]), row)
    if parent_binding is not None:
        _assert_binding_unchanged(
            Path(source_qualification["parent_segment_pass"]["path"]),
            parent_binding)
    final_state_tree = _assert_exact_state_tree(
        state_resolved, (Path(row["path"]) for row in output_rows))
    _require(final_state_tree == state_tree,
             "state directory tree changed during segment qualification")
    bindings = {
        "plan": plan_binding,
        "launch_record": launch_binding,
        "completion_record": completion_binding,
        "source_restart": source_binding,
        "endpoint_restart": endpoint_binding,
        "repository": repository_binding,
        "binary": binary_binding,
        "trajectory": trajectory_binding,
        "source_baryon_mass_evidence": source_mass_evidence,
        "parent_segment_pass": parent_binding,
        "run_log": run_binding,
        "event_log": event_binding,
        "exit_status": exit_binding,
        "gpu_before": gpu_before_binding,
        "gpu_after": gpu_after_binding,
        "tools": [_hash_record(script_path), _hash_record(audit_tool),
                  _hash_record(integrity_tool), _hash_record(metadata_tool),
                  _hash_record(layout_tool)],
        "planning_tools": plan.get("tools"),
        "verified_planning_tools": planned_tool_bindings,
        "verified_current_tools": current_tool_bindings,
    }
    return {
        "schema": SCHEMA,
        "kind": "athenak_segment_pass",
        "status": "pass",
        "qualification_mode": "complete_segment_v1",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "expected": expected,
        "bindings": bindings,
        "source_restart_audit": source_audit,
        "parent_qualification_audit": parent_qualification_audit,
        "endpoint_restart_audit": endpoint_audit,
        "restart_contract_audit": {
            "source": source_contract,
            "endpoint": endpoint_contract,
        },
        "run_log_audit": run,
        "launch_record_audit": launch,
        "completion_record_audit": completion,
        "event_log_audit": events,
        "gpu_audit": gpus,
        "parameter_audit": parameters,
        "scientific_threshold_audit": scientific,
        "scientific_advisories": scientific_advisories,
        "binary_topology_pair_audit": binary_pairs,
        "binary_selection_coverage_audit": selection_coverage,
        "state_directory_audit": state_tree,
        "output_inventory": output_rows,
    }


def publish_report(path: Path, report: dict[str, Any]) -> None:
    _require(path.name.endswith(".pass.ready"),
             "pass report filename must end in .pass.ready")
    _reject_output_ancestors(path)
    try:
        install_immutable_json(path, report)
    except FileExistsError as exc:
        raise CheckFailure(f"refusing to overwrite report: {path}") from exc


def _emit_stderr(status: str, message: str) -> None:
    print(json.dumps({
        "schema": SCHEMA,
        "kind": "athenak_segment_check",
        "status": status,
        "message": message,
    }, sort_keys=True), file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--launch-record", type=Path, required=True)
    parser.add_argument("--completion-record", type=Path, required=True)
    parser.add_argument("--endpoint-restart", type=Path, required=True)
    parser.add_argument("--state-dir", type=Path, required=True)
    parser.add_argument("--run-log", type=Path, required=True)
    parser.add_argument("--event-log", type=Path, required=True)
    parser.add_argument("--exit-status", type=Path, required=True)
    parser.add_argument("--gpu-before", type=Path, required=True)
    parser.add_argument("--gpu-after", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        report = check_segment(args)
        publish_report(args.output, report)
    except NotReady as exc:
        _emit_stderr("retry", str(exc))
        return RETRY_EXIT_STATUS
    except (CheckFailure, OSError, RuntimeError, ValueError) as exc:
        _emit_stderr("fail", str(exc))
        return 1
    print(json.dumps({
        "schema": SCHEMA,
        "kind": "athenak_segment_check",
        "status": "pass",
        "report": str(args.output.resolve()),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
