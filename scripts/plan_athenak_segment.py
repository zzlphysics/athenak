#!/usr/bin/env python3
"""Create a fail-closed, immutable execution plan for one AthenaK segment."""

from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path, PurePosixPath
import pwd
import re
import shutil
import stat
import struct
import subprocess
import sys
from typing import Any


# Importing a repository tool must not make an otherwise-clean production clone dirty.
sys.dont_write_bytecode = True
SCRIPT_PATH = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT_PATH.parent
READER_PATH = SCRIPT_DIR / "read_athenak_restart_metadata.py"
READER_SPEC = importlib.util.spec_from_file_location(
    "athenak_segment_restart_metadata", READER_PATH
)
if READER_SPEC is None or READER_SPEC.loader is None:  # pragma: no cover - installation
    raise RuntimeError(f"cannot load restart metadata reader: {READER_PATH}")
READER = importlib.util.module_from_spec(READER_SPEC)
sys.modules[READER_SPEC.name] = READER
READER_SPEC.loader.exec_module(READER)
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
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
from audit_athenak_restart import audit_restart


SCHEMA_VERSION = 1
GIT_PATH = Path(shutil.which("git") or "").resolve(strict=True)
MAX_ATHENA_CYCLE = 2**31 - 1
MIN_CAPACITY_HEADROOM = 64
YELLOW_CAPACITY_HEADROOM = 128
GIB_BYTES = 1 << 30
DEFAULT_PLANNED_PEAK_OUTPUT_GIB = 200
DISK_PREFLIGHT_KIND = "statvfs_unique_filesystem_budget_v1"
DISK_PREFLIGHT_ACCOUNTING = "group_roles_by_st_dev_once_v1"
DISK_PREFLIGHT_FORMULA = (
    "per_filesystem_required_free_bytes=max(additional_hard_minimum_free_bytes,"
    "max(minimum_reserve_bytes,minimum_reserve_restart_multiples*"
    "source_restart_size_bytes)+sum(role_contribution_bytes))"
)
GIT_ENVIRONMENT = {
    "LANG": "C",
    "LC_ALL": "C",
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_CONFIG_GLOBAL": "/dev/null",
}
GIT_CONFIGURATION = (
    "-c", "core.fsmonitor=false",
    "-c", "core.untrackedCache=false",
    "-c", "core.hooksPath=/dev/null",
)


def _file_record(path: Path, *, executable: bool = False) -> dict[str, Any]:
    resolved = path.expanduser().resolve(strict=True)
    if not resolved.is_file():
        raise ValueError(f"not a regular file: {resolved}")
    if executable and not os.access(resolved, os.X_OK):
        raise ValueError(f"AthenaK binary is not executable: {resolved}")
    audited = stable_sha256(resolved)
    return {
        "path": str(resolved),
        "size": audited["size"],
        "sha256": audited["sha256"],
        "closure_check": audited["closure_check"],
    }


def _canonical_executable(name: str) -> Path:
    discovered = shutil.which(name)
    if not discovered:
        raise ValueError(f"required executable is unavailable: {name}")
    try:
        resolved = Path(discovered).resolve(strict=True)
    except OSError as exc:
        raise ValueError(f"cannot resolve required executable {name}: {exc}") from exc
    if not resolved.is_file() or not os.access(resolved, os.X_OK):
        raise ValueError(f"required executable is not a regular executable: {resolved}")
    return resolved


def _run_git(repo: Path, *arguments: str) -> str:
    command = [str(GIT_PATH), "-C", str(repo), *GIT_CONFIGURATION, *arguments]
    try:
        result = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=GIT_ENVIRONMENT,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        detail = getattr(exc, "stderr", "") or str(exc)
        raise ValueError(f"git command failed: {' '.join(command)}: {detail.strip()}") \
            from exc
    return result.stdout.strip()


def _repository_record(path: Path) -> dict[str, Any]:
    repo = path.expanduser().resolve(strict=True)
    if not repo.is_dir():
        raise ValueError(f"repository is not a directory: {repo}")
    top = Path(_run_git(repo, "rev-parse", "--show-toplevel")).resolve(strict=True)
    if top != repo:
        raise ValueError(f"--repo must name the Git worktree root: {repo} != {top}")
    commit = _run_git(repo, "rev-parse", "--verify", "HEAD")
    if re.fullmatch(r"[0-9a-f]{40,64}", commit) is None:
        raise ValueError(f"invalid Git commit identity: {commit!r}")
    status = _run_git(repo, "status", "--porcelain=v1", "--untracked-files=all")
    if status:
        first = status.splitlines()[0]
        raise ValueError(f"repository is not clean (first entry: {first})")
    return {"path": str(repo), "commit": commit, "clean": True}


def _launch_environment() -> dict[str, Any]:
    """Bind the sole explicit environment permitted for MPI/AthenaK."""

    try:
        passwd_home = pwd.getpwuid(os.geteuid()).pw_dir
    except KeyError as exc:
        raise ValueError("effective uid has no passwd database entry") from exc
    if not passwd_home or not Path(passwd_home).is_absolute():
        raise ValueError("passwd HOME must be a nonempty absolute path")
    absolute = Path(passwd_home).absolute()
    try:
        home = absolute.resolve(strict=True)
        info = absolute.lstat()
    except OSError as exc:
        raise ValueError(f"cannot bind passwd HOME {absolute}: {exc}") from exc
    if (absolute != home or str(home) != passwd_home or
            not stat.S_ISDIR(info.st_mode) or stat.S_ISLNK(info.st_mode) or
            info.st_uid != os.geteuid() or stat.S_IMODE(info.st_mode) & 0o022):
        raise ValueError(
            "passwd HOME must be canonical, launcher-owned, non-symlink, and not "
            "group/other writable")
    values = {
        "HOME": str(home),
        "LANG": "C",
        "LC_ALL": "C",
        "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
    }
    return {
        "kind": "explicit_values_v1",
        "values": values,
        "sha256": _canonical_sha256(values),
    }


def _strict_false(parameters: dict[str, dict[str, str]], block: str, key: str) -> None:
    try:
        value = parameters[block][key].strip().lower()
    except KeyError as exc:
        raise ValueError(f"restart parameter <{block}>/{key} is required") from exc
    if value not in ("false", "0"):
        raise ValueError(f"restart parameter <{block}>/{key} must be false")


def _positive_finite(value: float, label: str) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{label} must be finite and positive")
    return value


def _positive_integer(value: int, label: str) -> int:
    if value < 1:
        raise ValueError(f"{label} must be positive")
    return value


def _disk_preflight_contract(source_restart_size: int, trajectory_size: int,
                             planned_peak_output_gib: int) -> dict[str, Any]:
    """Bind the conservative, per-filesystem launch capacity budget."""

    peak_gib = _positive_integer(
        planned_peak_output_gib, "--planned-peak-output-gib")
    peak_bytes = peak_gib * GIB_BYTES
    return {
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
        "bound_directory_fds": {"state_dir": 202, "staging_dir": 206},
    }


def _sequential_add(initial: float, increment: float, count: int) -> float:
    """Match C++ binary64 ``time += root_dt`` without algebraic reassociation."""

    result = float(initial)
    for _ in range(count):
        result += increment
    return result


def _wall_time_token(seconds: int) -> str:
    """Format the wall limit exactly as AthenaK's ``%d:%d:%d`` parser expects."""

    _positive_integer(seconds, "--wall-time-seconds")
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds_part = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds_part:02d}"


def _read_stable_bytes(path: Path) -> tuple[bytes, dict[str, Any]]:
    """Read one closed regular file through the descriptor whose identity is bound."""

    checked_path, stream, signature = _open_regular_nofollow(path)
    digest = hashlib.sha256()
    pieces: list[bytes] = []
    try:
        exempt = {(os.getpid(), stream.fileno())}
        _assert_closed(checked_path, signature, exempt)
        for piece in iter(lambda: stream.read(HASH_CHUNK_BYTES), b""):
            pieces.append(piece)
            digest.update(piece)
        _assert_stream_signature(stream, checked_path, signature, "while reading")
        _assert_path_signature(checked_path, signature, "while reading")
        _assert_closed(checked_path, signature, exempt)
    finally:
        stream.close()
    _assert_path_signature(checked_path, signature, "after reading")
    _assert_closed(checked_path, signature)
    return b"".join(pieces), {
        "path": str(checked_path), **signature, "sha256": digest.hexdigest(),
        "closure_check": "linux_proc_fd",
    }


def _load_parent_segment_pass(path: Path, source_restart: dict[str, Any],
                              source_audit: dict[str, Any],
                              source_time: float, source_cycle: int,
                              source_baryon: dict[str, Any],
                              trajectory: dict[str, Any]
                              ) -> tuple[dict[str, Any], dict[str, Any]]:
    """Qualify the immutable predecessor and bind its audited history output."""

    parent_path = path.expanduser().absolute()
    if not parent_path.name.endswith(".pass.ready"):
        raise ValueError("--parent-segment-pass must end in .pass.ready")
    parent_stat = parent_path.lstat()
    if not stat.S_ISREG(parent_stat.st_mode) or stat.S_ISLNK(parent_stat.st_mode):
        raise ValueError("parent segment pass must be a regular non-symlink file")
    if stat.S_IMODE(parent_stat.st_mode) & 0o222:
        raise ValueError("parent segment pass must be immutable (no write bits)")
    raw, binding = _read_stable_bytes(parent_path)
    try:
        parent = strict_json_loads(raw, str(parent_path))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"parent segment pass is invalid JSON: {exc}") from exc
    if not isinstance(parent, dict):
        raise ValueError("parent segment pass JSON root must be an object")
    if (parent.get("schema") != SCHEMA_VERSION or
            parent.get("kind") != "athenak_segment_pass" or
            parent.get("status") != "pass"):
        raise ValueError("parent segment pass has unsupported schema/kind/status")
    endpoint_binding = parent.get("bindings", {}).get("endpoint_restart")
    endpoint_audit = parent.get("endpoint_restart_audit")
    expected = parent.get("expected")
    if not isinstance(endpoint_binding, dict) or not isinstance(endpoint_audit, dict):
        raise ValueError("parent segment pass lacks endpoint restart evidence")
    if not isinstance(expected, dict):
        raise ValueError("parent segment pass lacks expected endpoint evidence")
    if (endpoint_binding.get("path") != source_restart["path"] or
            endpoint_binding.get("size") != source_restart["size"] or
            endpoint_binding.get("sha256") != source_restart["sha256"]):
        raise ValueError("parent endpoint binding does not identify --source-restart")
    parent_trajectory = parent.get("bindings", {}).get("trajectory")
    if (not isinstance(parent_trajectory, dict) or
            parent_trajectory.get("size") != trajectory["size"] or
            parent_trajectory.get("sha256") != trajectory["sha256"]):
        raise ValueError(
            "parent pass trajectory content does not identify --trajectory")
    audit_metadata = endpoint_audit.get("metadata")
    audit_layout = endpoint_audit.get("layout")
    stored = endpoint_audit.get("stored_reals")
    if (endpoint_audit.get("valid") is not True or
            endpoint_audit.get("sha256") != source_restart["sha256"] or
            endpoint_audit.get("signature", {}).get("size") !=
            source_restart["size"] or
            not isinstance(audit_metadata, dict) or
            audit_metadata.get("cycle") != source_cycle or
            audit_metadata.get("time") != source_time or
            not isinstance(audit_layout, dict) or
            audit_layout.get("expected_file_size") != source_restart["size"] or
            not isinstance(stored, dict) or
            stored.get("nonfinite_count") != 0 or
            stored.get("finite_count") != stored.get("count") or
            endpoint_audit.get("topology", {}).get(
                "complete_leaf_coverage") is not True):
        raise ValueError("parent endpoint audit is invalid or incomplete")
    if (expected.get("final_cycle") != source_cycle or
            expected.get("final_time") != source_time):
        raise ValueError("parent planned endpoint does not match --source-restart")
    # Reconcile the parent's embedded audit with the new planner's independent live
    # audit, so a structurally plausible hand-written pass cannot qualify the source.
    if (endpoint_audit.get("layout") != source_audit.get("layout") or
            endpoint_audit.get("stored_reals") != source_audit.get("stored_reals") or
            endpoint_audit.get("topology") != source_audit.get("topology")):
        raise ValueError("parent endpoint audit differs from the live source audit")
    inventory = parent.get("output_inventory")
    if not isinstance(inventory, list):
        raise ValueError("parent segment pass lacks output inventory")
    baryon_histories = [
        row for row in inventory
        if isinstance(row, dict) and
        isinstance(row.get("history"), dict) and
        "baryon_m" in row["history"].get("columns", [])
    ]
    if len(baryon_histories) != 1:
        raise ValueError("parent pass must bind exactly one baryon history output")
    history = baryon_histories[0]
    source_history = source_baryon["evidence"]
    if not all(history.get(name) == source_history.get(name)
               for name in ("path", "size", "sha256")):
        raise ValueError(
            "--source-history is not the parent pass audited baryon history output")
    try:
        inventory_last = float(history["history"]["column_last"]["baryon_m"])
        scientific_last = float(
            parent["scientific_threshold_audit"]["baryon_mass"]["last"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("parent pass lacks audited endpoint baryon mass") from exc
    if (not math.isfinite(inventory_last) or
            not math.isfinite(scientific_last) or
            inventory_last != source_baryon["value"] or
            scientific_last != source_baryon["value"]):
        raise ValueError(
            "parent audited baryon mass differs from independently parsed "
            "--source-history")
    return parent, {
        "path": binding["path"], "size": binding["size"],
        "sha256": binding["sha256"],
        "closure_check": binding["closure_check"],
    }


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _read_source_baryon_mass(path: Path, expected_time: float) -> dict[str, Any]:
    raw, audited = _read_stable_bytes(path)
    try:
        lines = raw.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise ValueError(f"source history is not UTF-8: {path}") from exc
    headers = [re.findall(r"\[\d+\]=([^\s]+)", line)
               for line in lines if line.lstrip().startswith("#")]
    headers = [row for row in headers if row]
    if len(headers) != 1 or "time" not in headers[0] or "baryon_m" not in headers[0]:
        raise ValueError("source history must have one header with time and baryon_m")
    rows: list[list[float]] = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        try:
            values = [float(value) for value in stripped.split()]
        except ValueError as exc:
            raise ValueError("source history contains a malformed value") from exc
        if len(values) != len(headers[0]) or not all(math.isfinite(v) for v in values):
            raise ValueError("source history row is non-finite or has wrong width")
        rows.append(values)
    if not rows:
        raise ValueError("source history has no data rows")
    columns = {name: index for index, name in enumerate(headers[0])}
    last = rows[-1]
    if last[columns["time"]] != expected_time:
        raise ValueError("source history does not end at the source restart time")
    mass = last[columns["baryon_m"]]
    if not math.isfinite(mass) or mass <= 0.0:
        raise ValueError("source history endpoint baryon mass is not positive and finite")
    return {
        "value": mass,
        "time": last[columns["time"]],
        "evidence": {
            "path": audited["path"], "size": audited["size"],
            "sha256": audited["sha256"],
            "closure_check": audited["closure_check"],
        },
    }


def _bool_parameter(block: dict[str, str], key: str, default: bool = False) -> bool:
    text = block.get(key, "true" if default else "false").strip().lower()
    if text in ("true", "1"):
        return True
    if text in ("false", "0"):
        return False
    raise ValueError(f"invalid boolean output parameter {key}={block.get(key)!r}")


def _output_template(
        basename: str, block: dict[str, str]) -> tuple[bool, str | None, bool]:
    """Return numbered/template/binary-inspection metadata for supported outputs."""

    file_type = block.get("file_type", "").strip()
    if not file_type:
        raise ValueError("output block has no file_type")
    variable = block.get("variable", "").strip()
    file_id = block.get("id", variable).strip()
    per_rank = _bool_parameter(block, "single_file_per_rank", False)
    if per_rank:
        raise ValueError("per-rank segment outputs are outside this campaign contract")
    if file_type == "bin":
        if not file_id:
            raise ValueError("binary output block has no variable/id")
        return True, f"bin/{basename}.{file_id}.{{file_number:05d}}.bin", True
    if file_type == "rst":
        return True, f"rst/{basename}.{{file_number:05d}}.rst", False
    if file_type == "log":
        return False, f"{basename}.log", False
    if file_type == "hst":
        # One hst block may emit several physics-specific files.  Their exact paths are
        # deliberately supplied with --required-unnumbered and stored at plan top level.
        return False, None, False
    if file_type == "tab":
        if not file_id:
            raise ValueError("tab output block has no variable/id")
        return True, f"tab/{basename}.{file_id}.{{file_number:05d}}.tab", False
    if file_type == "vtk":
        if not file_id:
            raise ValueError("vtk output block has no variable/id")
        gid = block.get("gid")
        gid_component = f"{gid}." if gid is not None and int(gid) >= 0 else ""
        return True, (
            f"vtk/{basename}.{file_id}.{gid_component}{{file_number:05d}}.vtk"
        ), False
    if file_type == "cart":
        if not file_id:
            raise ValueError("cart output block has no variable/id")
        return True, f"cart/{basename}.{file_id}.{{file_number:05d}}.bin", False
    raise ValueError(f"unsupported output file_type in strict segment plan: {file_type}")


def _output_blocks(
        parameters: dict[str, dict[str, str]],
        required_unnumbered: list[str]) -> list[dict[str, Any]]:
    try:
        basename = parameters["job"]["basename"].strip()
    except KeyError as exc:
        raise ValueError("restart parameter <job>/basename is required") from exc
    if not basename or "/" in basename:
        raise ValueError("restart <job>/basename is empty or contains a slash")
    required_set = set(required_unnumbered)
    outputs: list[dict[str, Any]] = []
    names = sorted(
        (name for name in parameters if re.fullmatch(r"output[0-9]+", name)),
        key=lambda name: int(name[6:]),
    )
    for name in names:
        snapshot = copy.deepcopy(parameters[name])
        has_dt = "dt" in snapshot
        has_dcycle = "dcycle" in snapshot
        if has_dcycle:
            # Outputs::Outputs deliberately gives dcycle precedence when a restart
            # serializes both keys.  Preserve that general AthenaK behavior here;
            # the canonical output3 stream is validated separately and may not use
            # dcycle at all.
            try:
                cadence = _positive_integer(
                    int(snapshot["dcycle"]), f"<{name}>/dcycle"
                )
            except ValueError as exc:
                raise ValueError(f"<{name}>/dcycle must be a positive integer") from exc
            cadence_mode = "dcycle"
        elif has_dt:
            cadence = _positive_finite(float(snapshot["dt"]), f"<{name}>/dt")
            cadence_mode = "dt"
        else:
            raise ValueError(f"<{name}> must specify dt or dcycle")
        numbered, template, inspect_binary = _output_template(basename, snapshot)
        file_type = snapshot["file_type"].strip()
        matching_required: list[str] = []
        if not numbered:
            if template is not None and template in required_set:
                matching_required.append(template)
            elif file_type == "hst":
                prefix = f"{basename}."
                matching_required = sorted(
                    value for value in required_set
                    if value.startswith(prefix) and value.endswith(".hst")
                )
                if not matching_required:
                    raise ValueError(
                        f"<{name}> requires at least one matching "
                        "--required-unnumbered history path"
                    )
                # The checker consumes one template per output block and separately
                # inventories every required unnumbered file.  Bind this block to the
                # first explicit history artifact rather than guessing physics modules.
                template = matching_required[0]
        outputs.append({
            "block": name,
            "index": int(name[6:]),
            "file_type": file_type,
            "id": snapshot.get("id", snapshot.get("variable")),
            "enabled": True,
            "cadence_mode": cadence_mode,
            "cadence": cadence,
            "numbered": numbered,
            "relative_path_template": template,
            "required_unnumbered_paths": matching_required,
            "inspect_binary": inspect_binary,
            "parameters": snapshot,
            "parameters_sha256": _canonical_sha256(snapshot),
        })
    if not outputs:
        raise ValueError("restart contains no enabled output blocks")
    return outputs


def _runtime_parameters_with_topology_reference(
        source: dict[str, dict[str, str]], root_dt: float,
        *, anchor: bool) -> dict[str, dict[str, str]]:
    """Apply the sole exact scientific-output override used by strict segments."""

    runtime = copy.deepcopy(source)
    block = runtime.get("output3")
    if not isinstance(block, dict):
        raise ValueError("strict segment requires canonical <output3> divB stream")
    try:
        gid = int(block.get("gid", "-1"))
    except ValueError as exc:
        raise ValueError("canonical <output3>/gid is invalid") from exc
    if (block.get("file_type", "").strip() != "bin" or
            block.get("variable", "").strip() != "mhd_divb" or
            gid >= 0 or "region_center" in block or
            any(f"slice_x{axis}" in block for axis in (1, 2, 3)) or
            block.get("ghost_zones", "false").strip().lower() not in
            ("false", "0")):
        raise ValueError(
            "canonical <output3> must be a ghost-free full-domain mhd_divb binary "
            "stream")
    if "dcycle" in block:
        raise ValueError("canonical <output3> may not contain dcycle")
    try:
        source_dt = float(block["dt"])
    except (KeyError, ValueError) as exc:
        raise ValueError("canonical <output3> must contain dt") from exc
    if not math.isfinite(source_dt) or source_dt < 0.0:
        raise ValueError("source <output3>/dt must be finite and nonnegative")
    if not anchor and source_dt != root_dt:
        raise ValueError(
            "parent-chain source restart must serialize output3/dt=root_dt")
    block["dt"] = repr(root_dt)
    return runtime


def _expected_binary_variables(variable: str,
                               parameters: dict[str, dict[str, str]]) -> list[str]:
    """Derive the exact strict-BBH binary labels from BaseTypeOutput semantics."""

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
            raise ValueError("strict mhd_w_bcc output has invalid mhd/nscalars") from exc
        dynamic = adm.get("dynamic", "false").strip().lower() in ("true", "1")
        if mhd.get("eos", "").strip() != "ideal" or nscalars != 0 or not dynamic:
            raise ValueError(
                "strict BBH mhd_w_bcc contract requires ideal EOS, zero scalars, "
                "and dynamic GR")
        return [
            "dens", "velx", "vely", "velz", "press",
            "bcc1", "bcc2", "bcc3", "temperature",
        ]
    raise ValueError(f"unsupported strict BBH binary variable group: {variable!r}")


def _float32(value: float) -> float:
    try:
        return struct.unpack("=f", struct.pack("=f", value))[0]
    except (OverflowError, struct.error) as exc:
        raise ValueError(f"output time is not representable as float32: {value}") from exc


def _replay_output_contract(outputs: list[dict[str, Any]], source_time: float,
                            source_cycle: int, final_time: float, final_cycle: int,
                            root_dt: float, tlim: float) -> None:
    """Seal the exact Execute/Finalize output schedule into each output snapshot."""

    states: dict[str, dict[str, Any]] = {}
    for output in outputs:
        block = output["block"]
        parameters = output["parameters"]
        try:
            state = {
                "file_number": int(parameters.get("file_number", "0")),
                "last_time": float(parameters.get("last_time", "-1")),
                "last_write_cycle": int(parameters.get("last_write_cycle", "-1")),
                "dt": float(parameters["dt"])
                      if output["cadence_mode"] == "dt" else 0.0,
                "dcycle": int(parameters["dcycle"])
                          if output["cadence_mode"] == "dcycle" else 0,
                "wrote_this_run": False,
                "writes": [],
            }
        except (KeyError, ValueError) as exc:
            raise ValueError(f"<{block}> has invalid serialized cadence state") from exc
        if (state["file_number"] < 0 or not math.isfinite(state["last_time"]) or
                state["dt"] < 0.0 or state["dcycle"] < 0):
            raise ValueError(f"<{block}> has invalid serialized cadence state")
        states[block] = state

    def due(state: dict[str, Any], cycle: int, output_time: float,
            enforce_time_limit: bool) -> bool:
        time_due = (state["dt"] > 0.0 and
                    _float32(output_time) >=
                    _float32(state["last_time"] + state["dt"]) and
                    (not enforce_time_limit or
                     _float32(output_time) < _float32(tlim)))
        cycle_due = state["dcycle"] > 0 and cycle % state["dcycle"] == 0
        return time_due or cycle_due

    def write(output: dict[str, Any], state: dict[str, Any], cycle: int,
              output_time: float, kind: str, advance: bool) -> None:
        row: dict[str, Any] = {"cycle": cycle, "time": output_time, "kind": kind}
        if output["numbered"]:
            row["file_number"] = state["file_number"]
            state["file_number"] += 1
        state["writes"].append(row)
        if advance:
            if state["last_time"] < 0.0:
                state["last_time"] = output_time
            else:
                state["last_time"] += state["dt"]
        state["last_write_cycle"] = cycle
        state["wrote_this_run"] = True

    output_time = source_time
    for cycle in range(source_cycle + 1, final_cycle + 1):
        output_time += root_dt
        for output in outputs:
            state = states[output["block"]]
            if due(state, cycle, output_time, True):
                write(output, state, cycle, output_time, "scheduled", True)

    non_restarts = sorted(
        (output for output in outputs if output["file_type"] != "rst"),
        key=lambda output: output["index"], reverse=True)
    restarts = [output for output in outputs if output["file_type"] == "rst"]
    if len(restarts) != 1:
        raise ValueError("strict segment plan requires exactly one restart output")
    final_parameter_state_changed = False
    for output in [*non_restarts, *restarts]:
        state = states[output["block"]]
        is_restart = output["file_type"] == "rst"
        wrote_current = (state["wrote_this_run"] and
                         state["last_write_cycle"] == final_cycle)
        if wrote_current and (not is_restart or not final_parameter_state_changed):
            continue
        advance = (state["last_write_cycle"] != final_cycle and
                   due(state, final_cycle, final_time, False))
        kind = "restart_final_rewrite" if is_restart and wrote_current else "forced_final"
        write(output, state, final_cycle, final_time, kind, advance)
        if not is_restart:
            final_parameter_state_changed = True

    for output in outputs:
        state = states[output["block"]]
        output["expected_writes"] = state["writes"]
        output["expected_endpoint_state"] = {
            "file_number": state["file_number"],
            "last_time": state["last_time"],
            "last_write_cycle": state["last_write_cycle"],
        }


def _required_paths(values: list[str]) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        path = PurePosixPath(value)
        if (not value or path.is_absolute() or value != path.as_posix()
                or any(part in ("", ".", "..") for part in path.parts)):
            raise ValueError(
                f"--required-unnumbered must be a normalized relative path: {value!r}"
            )
        if value in seen:
            raise ValueError(f"duplicate --required-unnumbered path: {value}")
        seen.add(value)
        result.append(value)
    return sorted(result)


def _policy(root_dt: float, ranks: int) -> dict[str, Any]:
    hard_ratios = [
        {
            "name": "fofc_per_test",
            "numerator": "fofc",
            "denominator": "fofc_tests",
            "max_ratio": 0.005,
        },
        {
            "name": "cons_adjust_per_c2p_call",
            "numerator": "cons_adjust",
            "denominator": "c2p_calls",
            "max_ratio": 0.005,
        },
        {
            "name": "mag_adjust_per_c2p_call",
            "numerator": "mag_adjust",
            "denominator": "c2p_calls",
            "max_ratio": 0.005,
        },
    ]
    yellow_ratios = [
        {**item, "max_ratio": 0.001, "consecutive_rows": 3}
        for item in hard_ratios
    ]
    return {
        "ranks": ranks,
        "endpoint_time_ulp_tolerance": 0,
        "segment_termination": {
            "primary": "cycle_limit",
            "endpoint_time_max_ulps": 0,
            "tlim_role": "nonbinding_guard",
            "require_exact_final_cycle": True,
            "require_exact_binary64_endpoint": True,
        },
        "fixed_root_dt_mode": {
            "enabled": True,
            "root_dt": root_dt,
            "root_dt_parameter": "time/root_dt_max",
            "require_source_last_dt_exact": True,
            "sequential_binary64_addition": True,
        },
        "restart_contract": {
            "real_bytes": 8,
            "subcycling": "level",
            "allow_legacy_subcycling_restart": False,
            "allow_legacy_ghost_event_counters": False,
            "amr_counter_version": READER.AMR_COUNTER_VERSION,
            "event_counter_version": READER.EVENT_COUNTER_VERSION,
            "restart_cache_contract_version": READER.RESTART_CACHE_CONTRACT_VERSION,
            "pending_event_counters_all_zero": True,
            "level_subcycling_costs_exact": True,
        },
        "capacity": {
            "ranks": ranks,
            "minimum_headroom_blocks_hard": MIN_CAPACITY_HEADROOM,
            "minimum_headroom_blocks_yellow": YELLOW_CAPACITY_HEADROOM,
        },
        "output_integrity": {
            "minimum_closed_file_age_seconds": 120.0,
            "require_no_open_descriptors": True,
            "require_stable_size_mtime_while_hashing": True,
            "require_sha256": True,
            "refuse_manifest_overwrite": True,
        },
        "event_thresholds": hard_ratios,
        "event_absolute_thresholds": {
            "hard_equal_zero": ["eos_fail", "eos_vceil"],
            "c2p_iterations_exclusive_max": 25,
        },
        "yellow_event_thresholds": yellow_ratios,
        "nonfinite_count_max": 0,
        "divb_max_abs": {"hard": 1.0e-8, "yellow": 1.0e-11},
        "baryon_mass_fractional_loss": {
            "hard_per_root_step": 0.005,
            "yellow_per_root_step": 0.0025,
            "yellow_per_48M": 0.02,
            "rolling_window_root_steps": 10,
        },
        "gpu_exit_memory_mib_max": 100.0,
        "gpu_ecc": {
            "corrected_before_max": 0,
            "corrected_after_max": 0,
            "uncorrected_before_max": 0,
            "uncorrected_after_max": 0,
        },
        "remote_disk": {
            "warn_percent": 65,
            "do_not_start_percent": 75,
            "synchronized_stop_percent": 80,
            "emergency_stop_percent": 85,
            "minimum_reserve_gib": 50,
            "minimum_reserve_restart_multiples": 2,
        },
        "mutable_parameters": [
            "time/nlim",
            "time/tlim",
            "time/restart_dt_growth",
            "output*/file_number",
            "output*/last_time",
            "output*/last_write_cycle",
        ],
    }


def build_plan(args: argparse.Namespace) -> dict[str, Any]:
    root_steps = _positive_integer(args.root_steps, "--root-steps")
    guard_steps = _positive_integer(args.tlim_guard_steps, "--tlim-guard-steps")
    ranks = _positive_integer(args.ranks, "--ranks")
    root_dt = _positive_finite(args.root_dt, "--root-dt")
    required = _required_paths(args.required_unnumbered)

    repo = _repository_record(args.repo)
    binary = _file_record(args.binary, executable=True)
    trajectory = _file_record(args.trajectory)
    source_restart = _file_record(args.source_restart)
    launcher = _file_record(args.launcher, executable=True)
    state_dir = args.state_dir.expanduser().resolve(strict=True)
    if not state_dir.is_dir():
        raise ValueError(f"--state-dir is not a directory: {state_dir}")
    if any(state_dir.iterdir()):
        raise ValueError("--state-dir must be a new empty segment directory")
    staging_dir = args.staging_dir.expanduser().resolve(strict=True)
    evidence_dir = args.evidence_dir.expanduser().resolve(strict=True)
    for directory, label in ((staging_dir, "--staging-dir"),
                             (evidence_dir, "--evidence-dir")):
        if not directory.is_dir() or directory.is_symlink():
            raise ValueError(f"{label} must be a real directory: {directory}")
    plan_path = evidence_dir / "segment.plan.json"
    output_path = args.output.expanduser().absolute()
    if (not output_path.is_absolute() or output_path != plan_path or
            output_path.parent != evidence_dir):
        raise ValueError(
            f"--output must be the canonical direct evidence child {plan_path}")
    if os.path.lexists(output_path):
        raise ValueError(f"refusing to overwrite existing plan: {output_path}")
    for directory, label in ((staging_dir, "--staging-dir"),
                             (evidence_dir, "--evidence-dir")):
        if any(directory.iterdir()):
            raise ValueError(f"{label} must be a new empty directory")
    for left, right in ((state_dir, staging_dir), (state_dir, evidence_dir),
                        (staging_dir, evidence_dir)):
        for first, second in ((left, right), (right, left)):
            try:
                first.relative_to(second)
            except ValueError:
                continue
            raise ValueError("state, staging, and evidence directories must be separate")
    wall_time_seconds = _positive_integer(
        args.wall_time_seconds, "--wall-time-seconds")
    source_full_audit = audit_restart(Path(source_restart["path"]))
    if (source_full_audit.get("valid") is not True or
            source_full_audit.get("sha256") != source_restart["sha256"] or
            source_full_audit.get("signature", {}).get("size") !=
            source_restart["size"] or
            source_full_audit.get("layout", {}).get("expected_file_size") !=
            source_restart["size"] or
            source_full_audit.get("stored_reals", {}).get("nonfinite_count") != 0 or
            source_full_audit.get("stored_reals", {}).get("finite_count") !=
            source_full_audit.get("stored_reals", {}).get("count") or
            source_full_audit.get("topology", {}).get(
                "complete_leaf_coverage") is not True):
        raise ValueError(
            "source restart complete audit did not prove finite exact layout and "
            "complete AMR leaf coverage")
    metadata = READER.read_restart_metadata(Path(source_restart["path"]))
    parameters = metadata.parameters

    if metadata.real_bytes != 8:
        raise ValueError(f"source restart uses Real{metadata.real_bytes}, not Real8")
    if parameters.get("time", {}).get("subcycling", "").strip() != "level":
        raise ValueError("source restart does not use strict level subcycling")
    _strict_false(parameters, "time", "allow_legacy_subcycling_restart")
    _strict_false(parameters, "time", "allow_legacy_ghost_event_counters")
    if metadata.amr_cycle_counters is None:
        raise ValueError("source restart has no AMR counter extension")
    if metadata.event_counter_version != READER.EVENT_COUNTER_VERSION:
        raise ValueError(
            "source restart does not carry the required v2 event-counter extension"
        )
    if metadata.restart_cache_contract_version != \
            READER.RESTART_CACHE_CONTRACT_VERSION:
        raise ValueError("source restart has no v1 restart-cache contract")
    if metadata.event_counters is None:
        raise ValueError("source restart has no pending event counters")
    pending_events = metadata.event_counters.as_dict()
    if any(value != 0 for value in pending_events.values()):
        raise ValueError(f"source restart has nonzero pending event counters: {pending_events}")

    expected_costs = tuple(
        math.ldexp(1.0, location.level - metadata.root_level)
        for location in metadata.locations
    )
    if metadata.costs != expected_costs:
        raise ValueError("source restart MeshBlock costs are not exact level costs")
    try:
        capacity = int(parameters["mesh_refinement"]["max_nmb_per_rank"])
    except (KeyError, ValueError) as exc:
        raise ValueError("source restart has invalid max_nmb_per_rank") from exc
    if capacity < 1:
        raise ValueError("source restart max_nmb_per_rank must be positive")
    partitions = READER.load_balance(metadata.costs, ranks, capacity)
    minimum_headroom = min(capacity - row.blocks for row in partitions)
    if minimum_headroom < MIN_CAPACITY_HEADROOM:
        raise ValueError(
            f"rank capacity headroom {minimum_headroom} is below "
            f"hard minimum {MIN_CAPACITY_HEADROOM}"
        )

    try:
        serialized_root_dt = float(parameters["time"]["root_dt_max"])
    except (KeyError, ValueError) as exc:
        raise ValueError("source restart has invalid <time>/root_dt_max") from exc
    if not math.isfinite(serialized_root_dt) or serialized_root_dt != root_dt:
        raise ValueError(
            f"source root_dt_max {serialized_root_dt!r} is not exactly {root_dt!r}"
        )
    if metadata.last_dt != root_dt:
        raise ValueError(
            f"source last_dt {metadata.last_dt!r} is not fixed root dt {root_dt!r}"
        )
    source_trajectory = parameters.get("problem", {}).get("trajectory_file")
    if source_trajectory is None:
        raise ValueError("source restart has no <problem>/trajectory_file")

    source_history_path = args.source_history.expanduser().resolve(strict=True)
    source_baryon_mass = _read_source_baryon_mass(
        source_history_path, metadata.time)
    parent_segment_pass = getattr(args, "parent_segment_pass", None)
    anchor = bool(getattr(args, "anchor", False))
    if anchor == (parent_segment_pass is not None):
        raise ValueError("specify exactly one of --anchor and --parent-segment-pass")
    if anchor:
        try:
            serialized_trajectory = Path(source_trajectory).expanduser().resolve(
                strict=True)
        except OSError as exc:
            raise ValueError(
                f"anchor serialized trajectory path is unavailable: "
                f"{source_trajectory}") from exc
        if serialized_trajectory != Path(trajectory["path"]):
            raise ValueError(
                "--trajectory does not resolve to the anchor source restart "
                "trajectory_file")
        source_qualification: dict[str, Any] = {
            "mode": "anchor_full_audit",
            "audit": source_full_audit,
            "source_baryon_mass": source_baryon_mass,
        }
    else:
        _, parent_record = _load_parent_segment_pass(
            parent_segment_pass, source_restart, source_full_audit,
            metadata.time, metadata.cycle,
            source_baryon_mass, trajectory)
        source_qualification = {
            "mode": "parent_segment_pass",
            "audit": source_full_audit,
            "parent_segment_pass": parent_record,
            "source_baryon_mass": source_baryon_mass,
        }

    final_cycle = metadata.cycle + root_steps
    if final_cycle > MAX_ATHENA_CYCLE:
        raise ValueError(f"planned final cycle exceeds AthenaK int range: {final_cycle}")
    final_time = _sequential_add(metadata.time, root_dt, root_steps)
    tlim = _sequential_add(final_time, root_dt, guard_steps)
    if not math.isfinite(final_time) or not math.isfinite(tlim) or tlim <= final_time:
        raise ValueError("planned endpoint/tlim is non-finite or non-increasing")

    runtime_parameters = _runtime_parameters_with_topology_reference(
        parameters, root_dt, anchor=anchor)
    outputs = _output_blocks(runtime_parameters, required)
    for output in outputs:
        if output["inspect_binary"]:
            variable = output["parameters"].get("variable", "").strip()
            output["expected_binary_variables"] = _expected_binary_variables(
                variable, parameters)
    _replay_output_contract(
        outputs, metadata.time, metadata.cycle, final_time, final_cycle, root_dt, tlim)
    output3 = next((output for output in outputs
                    if output["block"] == "output3"), None)
    if (output3 is None or output3["cadence_mode"] != "dt" or
            output3["cadence"] != root_dt or
            [write["cycle"] for write in output3["expected_writes"]] !=
            list(range(metadata.cycle + 1, final_cycle + 1)) or
            any(write["kind"] != "scheduled"
                for write in output3["expected_writes"])):
        raise ValueError(
            "canonical output3/dt=root_dt does not produce exactly one scheduled "
            "full-domain topology reference per root cycle")
    parameter_snapshot = copy.deepcopy(parameters)
    parameter_sha = _canonical_sha256(parameter_snapshot)
    source_summary = READER.metadata_dict(metadata, ranks, capacity)
    source_summary["parameters"] = parameter_snapshot
    source_summary["parameters_sha256"] = parameter_sha
    mpi_argv = [launcher["path"], "--allow-run-as-root", "--bind-to", "none",
                "-np", str(ranks)]
    staged_source_restart = {
        "path": str(staging_dir / "source.rst"),
        "size": source_restart["size"], "sha256": source_restart["sha256"],
    }
    staged_trajectory = {
        "path": str(staging_dir / "trajectory.dat"),
        "size": trajectory["size"], "sha256": trajectory["sha256"],
    }
    holder_pid_token = "{holder_pid}"
    source_proc_template = f"/proc/{holder_pid_token}/fd/200"
    trajectory_proc_template = f"/proc/{holder_pid_token}/fd/201"
    state_proc_template = f"/proc/{holder_pid_token}/fd/202"
    evidence_proc_template = f"/proc/{holder_pid_token}/fd/203"
    launcher_proc_template = f"/proc/{holder_pid_token}/fd/204"
    binary_proc_template = f"/proc/{holder_pid_token}/fd/205"
    athena_argv_template = [
        binary_proc_template, "--kokkos-map-device-id-by=mpi_rank",
        "-r", source_proc_template, "-d", state_proc_template,
        "-t", _wall_time_token(wall_time_seconds), f"time/nlim={final_cycle}",
        f"time/tlim={tlim!r}",
        f"problem/trajectory_file={trajectory_proc_template}",
        f"output3/dt={root_dt!r}",
    ]
    evidence = {
        "launch_record": str(evidence_dir / "segment.launch.ready"),
        "completion_record": str(evidence_dir / "segment.completion.ready"),
        "run_log": str(evidence_dir / "run.log"),
        "exit_status": str(evidence_dir / "exit.status"),
        "gpu_before": str(evidence_dir / "gpu-before.csv"),
        "gpu_after": str(evidence_dir / "gpu-after.csv"),
    }

    return {
        "schema": SCHEMA_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "policy": _policy(root_dt, ranks),
        "tools": {
            "planner": _file_record(SCRIPT_PATH),
            "segment_checker": _file_record(
                SCRIPT_DIR / "check_athenak_segment.py"),
            "segment_launcher": _file_record(
                SCRIPT_DIR / "launch_athenak_segment.py"),
            "restart_metadata_reader": _file_record(READER_PATH),
            "output_integrity": _file_record(SCRIPT_DIR / "output_integrity.py"),
            "restart_auditor": _file_record(SCRIPT_DIR / "audit_athenak_restart.py"),
            "restart_layout": _file_record(
                SCRIPT_DIR / "compare_athenak_restart_fields.py"),
            "git": _file_record(GIT_PATH, executable=True),
            "nvidia_smi": _file_record(
                _canonical_executable("nvidia-smi"), executable=True),
            "hash_algorithm": "sha256",
        },
        "inputs": {
            "repo": repo,
            "binary": binary,
            "trajectory": trajectory,
            "source_restart": source_restart,
        },
        "source_qualification": source_qualification,
        "source": source_summary,
        "expected": {
            "source_cycle": metadata.cycle,
            "source_time": metadata.time,
            "root_steps": root_steps,
            "root_dt": root_dt,
            "tlim_guard_steps": guard_steps,
            "final_cycle": final_cycle,
            "final_time": final_time,
            "tlim": tlim,
        },
        "command_overrides": {
            "time/nlim": final_cycle,
            "time/tlim": tlim,
            "output3/dt": root_dt,
        },
        "launch_contract": {
            "launcher": launcher,
            "mpi_argv": mpi_argv,
            "athena_argv_template": athena_argv_template,
            "state_dir": str(state_dir),
            "plan_path": str(plan_path),
            "environment": _launch_environment(),
            "directory_transport": {
                "kind": "linux_proc_holder_dirfd_v1",
                "holder_pid_token": holder_pid_token,
                "roles": {
                    "state_dir": {
                        "role": "state_dir", "fd": 202,
                        "planned_path": str(state_dir),
                        "proc_path_template": state_proc_template,
                    },
                    "evidence_dir": {
                        "role": "evidence_dir", "fd": 203,
                        "planned_path": str(evidence_dir),
                        "proc_path_template": evidence_proc_template,
                    },
                },
            },
            "executable_transport": {
                "kind": "linux_proc_holder_execfd_v1",
                "holder_pid_token": holder_pid_token,
                "roles": {
                    "launcher": {
                        "role": "launcher", "fd": 204,
                        "parent_path": launcher["path"],
                        "proc_path_template": launcher_proc_template,
                    },
                    "binary": {
                        "role": "binary", "fd": 205,
                        "parent_path": binary["path"],
                        "proc_path_template": binary_proc_template,
                    },
                },
            },
            "input_transport": {
                "kind": "linux_proc_holder_fd_v1",
                "holder_pid_token": holder_pid_token,
                "staging_dir": str(staging_dir),
                "roles": {
                    "source_restart": {
                        "role": "source_restart", "fd": 200,
                        "proc_path_template": source_proc_template,
                        "staged_file": staged_source_restart,
                    },
                    "trajectory": {
                        "role": "trajectory", "fd": 201,
                        "proc_path_template": trajectory_proc_template,
                        "staged_file": staged_trajectory,
                    },
                },
                "parent_content_identity": {
                    "source_restart_sha256": source_restart["sha256"],
                    "trajectory_sha256": trajectory["sha256"],
                    "source_serialized_trajectory_path": source_trajectory,
                },
                "trajectory_rebinding": {
                    "parameter": "problem/trajectory_file",
                    "parent_sha256": trajectory["sha256"],
                    "runtime_value_template": trajectory_proc_template,
                },
            },
            "disk_preflight": _disk_preflight_contract(
                source_restart["size"], trajectory["size"],
                args.planned_peak_output_gib),
            "evidence_dir": str(evidence_dir),
            "evidence": evidence,
            "wall_time_seconds": wall_time_seconds,
            "world_size": ranks,
            "gpu_count": ranks,
        },
        "parameter_snapshots": {
            "source_sha256": parameter_sha,
            "output_blocks": {
                output["block"]: output["parameters"] for output in outputs
            },
        },
        "outputs": outputs,
        "required_files": required,
    }


def _write_new_json(path: Path, payload: dict[str, Any]) -> None:
    output = path.expanduser().absolute()
    if os.path.lexists(output):
        raise ValueError(f"refusing to overwrite existing plan: {output}")
    try:
        install_immutable_json(output, payload)
    except FileExistsError as exc:
        raise ValueError(f"refusing to overwrite existing plan: {output}") from exc


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--source-restart", required=True, type=Path)
    parser.add_argument("--source-history", required=True, type=Path)
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument("--anchor", action="store_true")
    source_group.add_argument(
        "--parent-segment-pass", "--parent-pass",
        dest="parent_segment_pass", type=Path,
        help="immutable .pass.ready report whose endpoint is --source-restart",
    )
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--trajectory", required=True, type=Path)
    parser.add_argument("--root-steps", required=True, type=int)
    parser.add_argument("--root-dt", required=True, type=float)
    parser.add_argument("--tlim-guard-steps", required=True, type=int)
    parser.add_argument("--ranks", required=True, type=int)
    parser.add_argument("--launcher", required=True, type=Path)
    parser.add_argument("--state-dir", required=True, type=Path)
    parser.add_argument("--staging-dir", required=True, type=Path)
    parser.add_argument("--evidence-dir", required=True, type=Path)
    parser.add_argument("--wall-time-seconds", required=True, type=int)
    parser.add_argument(
        "--planned-peak-output-gib", type=int,
        default=DEFAULT_PLANNED_PEAK_OUTPUT_GIB,
        help=("conservative peak GiB of new state output only; staging copy and "
              "reserve are added separately (default: 200)"),
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--required-unnumbered",
        action="append",
        default=[],
        help="relative state-directory file that must exist after this segment",
    )
    args = parser.parse_args()
    try:
        plan = build_plan(args)
        _write_new_json(args.output, plan)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(args.output.expanduser().absolute())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
