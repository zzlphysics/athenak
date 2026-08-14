"""Adversarial unit tests for scheduled AthenaK prefix recovery."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import socket
import struct
from types import SimpleNamespace
from typing import Any

import pytest


ROOT = Path(__file__).resolve().parents[2]
import sys
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tst"))
import recover_athenak_segment_prefix as RECOVERY
import plan_athenak_segment as PLANNER


def _write_closed(path: Path, content: bytes) -> None:
    path.write_bytes(content)
    path.chmod(0o400)


def test_recovery_policy_is_identical_in_all_producers_and_consumers() -> None:
    assert RECOVERY.RECOVERY_POLICY == PLANNER.PREFIX_RECOVERY_POLICY
    assert RECOVERY.RECOVERY_POLICY == RECOVERY.CHECKER.PREFIX_RECOVERY_POLICY


def test_half_written_latest_falls_back_to_highest_complete() -> None:
    older = {
        "relative_path": "rst/run.00003.rst",
        "expected_write": {"cycle": 30},
        "classification": "complete",
        "audit": {"valid": True},
    }
    latest = {
        "relative_path": "rst/run.00004.rst",
        "expected_write": {"cycle": 40},
        "classification": "incomplete_truncated",
    }

    assert RECOVERY.select_highest_complete_candidate([older, latest]) is older


def test_later_complete_invalid_forbids_fallback() -> None:
    older = {
        "relative_path": "rst/run.00003.rst",
        "expected_write": {"cycle": 30},
        "classification": "complete",
    }
    latest = {
        "relative_path": "rst/run.00004.rst",
        "expected_write": {"cycle": 40},
        "classification": "complete_invalid",
        "failure": "non-finite stored field",
    }

    with pytest.raises(RECOVERY.RecoveryFailure, match="fallback.*forbidden"):
        RECOVERY.select_highest_complete_candidate([older, latest])


def test_highest_complete_is_selected_once_for_scientific_qualification() -> None:
    older = {
        "relative_path": "rst/run.00003.rst",
        "expected_write": {"cycle": 30},
        "classification": "complete",
    }
    latest = {
        "relative_path": "rst/run.00004.rst",
        "expected_write": {"cycle": 40},
        "classification": "complete",
    }

    selected = RECOVERY.select_highest_complete_candidate([older, latest])
    assert selected is latest
    # The production path validates only this returned candidate and has no retry loop.
    with pytest.raises(RECOVERY.RecoveryFailure, match="scientific"):
        raise RECOVERY.RecoveryFailure("scientific gate failed")


def test_obvious_short_restart_is_incomplete_not_invalid(tmp_path: Path) -> None:
    path = tmp_path / "latest.rst"
    _write_closed(path, b"half header")

    with pytest.raises(RECOVERY.CHECKER.RestartTruncationError):
        RECOVERY.CHECKER.read_restart_metadata(path)

    result = RECOVERY.classify_restart(path, minimum_complete_header_bytes=100)

    assert result["classification"] == "incomplete_truncated"


def test_complete_nonfinite_restart_is_invalid_not_incomplete(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "latest.rst"
    _write_closed(path, b"complete-size")
    monkeypatch.setattr(
        RECOVERY.CHECKER, "read_restart_metadata",
        lambda unused: SimpleNamespace(metadata_end=0, nmb_total=1))
    auditor_module_name = RECOVERY.CHECKER.audit_restart.__module__
    auditor_module = sys.modules[auditor_module_name]
    monkeypatch.setattr(
        auditor_module, "_derive_layout",
        lambda unused: SimpleNamespace(data_size=len(b"complete-size") - 8))

    def fail_nonfinite(unused):
        raise ValueError("non-finite stored Real")

    fail_nonfinite.__module__ = auditor_module_name
    monkeypatch.setattr(RECOVERY.CHECKER, "audit_restart", fail_nonfinite)

    result = RECOVERY.classify_restart(path, minimum_complete_header_bytes=1)

    assert result["classification"] == "complete_invalid"
    assert "non-finite" in result["failure"]


def test_execute_only_replay_never_adds_finalize_write() -> None:
    plan = {
        "expected": {
            "source_cycle": 0, "source_time": 0.0, "final_cycle": 6,
            "root_dt": 1.0, "tlim": 8.0,
        },
        "outputs": [
            {
                "block": "output3", "numbered": True, "cadence_mode": "dt",
                "parameters": {"dt": "1", "file_number": "0",
                               "last_time": "0", "last_write_cycle": "0"},
            },
            {
                "block": "output4", "numbered": True, "cadence_mode": "dt",
                "parameters": {"dt": "2", "file_number": "7",
                               "last_time": "0", "last_write_cycle": "0"},
            },
        ],
    }

    schedules, states, endpoint = RECOVERY.replay_execute_prefix(plan, 4)

    assert endpoint == 4.0
    assert [row["cycle"] for row in schedules["output3"]] == [1, 2, 3, 4]
    assert [row["cycle"] for row in schedules["output4"]] == [2, 4]
    assert all(row["kind"] == "scheduled"
               for rows in schedules.values() for row in rows)
    assert states["output4"]["file_number"] == 9


def test_history_truncation_before_target_is_fatal() -> None:
    raw = b"# [1]=time [2]=baryon_m\n1 2\n2 2\n"

    with pytest.raises(RECOVERY.RecoveryFailure, match="truncated"):
        RECOVERY._exact_text_prefix(raw, 3, "history")


def test_prefix_is_original_bytes_and_suffix_is_forensic() -> None:
    raw = b"# header\n1 2\n2 3\n3 4\npartial"

    prefix, suffix = RECOVERY._exact_text_prefix(raw, 2, "history")

    assert prefix == b"# header\n1 2\n2 3\n"
    assert suffix == b"3 4\npartial"
    assert prefix + suffix == raw


def test_event_prefix_keeps_target_telemetry_but_not_next_cycle_comments() -> None:
    raw = (
        b"# event header\n"
        b"# fofc_spatial_v1 kind=summary cycle=322 count=1 nfofc=1 "
        b"unattributed=0\n"
        b"# fofc_spatial_v1 kind=bin cycle=322 level_bin=8 count=1\n"
        b"322 1\n"
        b"# fofc_spatial_v1 kind=summary cycle=323 count=2 nfofc=2 "
        b"unattributed=0\n"
        b"# fofc_spatial_v1 kind=bin cycle=323 level_bin=9 count=2\n"
        b"323 2\n"
    )

    prefix, suffix = RECOVERY._exact_text_prefix(raw, 1, "event")

    assert b"kind=summary cycle=322" in prefix
    assert b"kind=bin cycle=322" in prefix
    assert prefix.endswith(b"322 1\n")
    assert b"cycle=323" not in prefix
    assert suffix.startswith(b"# fofc_spatial_v1 kind=summary cycle=323")
    assert prefix + suffix == raw


def test_run_log_prefix_replays_cache_and_fixed_root_steps(tmp_path: Path) -> None:
    path = tmp_path / "run.log"
    path.write_text(
        "tlim=2.400000e+02 nlim=40\n"
        "Strict-subcycling restart cache reconstruction passed: solver failures=0, "
        "non-finite proposed values=0, maximum raw component-relative proposed "
        "change=1.000000e-12, maximum absolute proposed change=1.000000e-12, "
        "maximum mixed-scale proposed change=1.000000e-12, mixed-scale acceptance "
        "tolerance=1.000000e-10.\n"
        "elapsed=0 cycle=30 time=1.440000e+02 dt=4.800000e+00\n"
        "elapsed=1 cycle=31 time=1.488000e+02 dt=4.800000e+00\n"
        "elapsed=2 cycle=32 time=1.536000e+02 dt=4.800000e+00\n",
        encoding="utf-8")
    path.chmod(0o400)
    expected = {
        "source_cycle": 30, "source_time": 144.0, "root_dt": 4.8,
        "final_cycle": 40, "tlim": 240.0, "recovery_final_cycle": 32,
    }

    result = RECOVERY._audit_run_log_prefix(path, expected)

    assert result["root_step_diagnostics"]["cycle_max"] == 32
    assert result["cache"]["solver_failures"] == 0
    assert result["original_limit_evidence"] == {
        "kind": "driver_limit_state_v1",
        "driver_limit_state_rows": 1,
        "nlim": 40,
        "tlim": 240.0,
    }


def test_unknown_extra_file_and_directory_are_rejected(tmp_path: Path) -> None:
    state = tmp_path / "state"
    state.mkdir()
    (state / "rst").mkdir()
    _write_closed(state / "rst" / "run.00001.rst", b"x")
    _write_closed(state / "surprise", b"x")

    with pytest.raises(RECOVERY.RecoveryFailure, match="unknown extra artifact"):
        RECOVERY._inventory_original_tree(state.resolve(), {"rst/run.00001.rst"})

    (state / "surprise").unlink()
    (state / "unexpected-directory").mkdir()
    with pytest.raises(RECOVERY.RecoveryFailure, match="unknown extra directory"):
        RECOVERY._inventory_original_tree(state.resolve(), {"rst/run.00001.rst"})


def test_boot_change_rejects_same_boot_closure(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    boot_id = tmp_path / "boot_id"
    boot_id.write_text("01234567-89ab-cdef-0123-456789abcdef\n", encoding="ascii")
    proc_stat = tmp_path / "stat"
    proc_stat.write_text("btime 200\n", encoding="ascii")
    monkeypatch.setattr(RECOVERY, "BOOT_ID_PATH", boot_id)
    monkeypatch.setattr(RECOVERY, "PROC_STAT_PATH", proc_stat)
    launch = {
        "created_utc": datetime.fromtimestamp(100, timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "mpirun_pid": 99999999, "mpirun_start_time_ticks": 1, "ranks": [],
    }

    with pytest.raises(RECOVERY.RecoveryFailure, match="boot changed"):
        RECOVERY._same_boot_closure(launch)


def test_same_boot_closure_requires_launcher_holder_identity(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    boot_id = tmp_path / "boot_id"
    boot_id.write_text("01234567-89ab-cdef-0123-456789abcdef\n", encoding="ascii")
    proc_stat = tmp_path / "stat"
    proc_stat.write_text("btime 1\n", encoding="ascii")
    monkeypatch.setattr(RECOVERY, "BOOT_ID_PATH", boot_id)
    monkeypatch.setattr(RECOVERY, "PROC_STAT_PATH", proc_stat)
    monkeypatch.setattr(
        RECOVERY, "_process_identity_observation",
        lambda pid, ticks: {"pid": pid, "recorded_start_time_ticks": ticks,
                            "observed_start_time_ticks": None,
                            "state": "disappeared",
                            "original_identity_gone": True})
    launch = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "mpirun_pid": 10, "mpirun_start_time_ticks": 100,
        "managed_process_group": {"pgid": 10, "new_session": True,
                                  "failure_cleanup":
                                  "SIGTERM_then_SIGKILL_with_quiescence_proof"},
        "input_transport": {"holder_pid": 11,
                            "holder_start_time_ticks": 101},
        "ranks": [{"global_rank": 0, "pid": 12,
                   "start_time_ticks": 102}],
    }
    monkeypatch.setattr(
        RECOVERY.os, "killpg",
        lambda pgid, number: (_ for _ in ()).throw(ProcessLookupError(pgid)))

    result = RECOVERY._same_boot_closure(launch)

    assert [row["pid"] for row in result["process_identities"]] == [10, 11, 12]
    assert result["process_identities"][1]["role"] == "launcher_holder"
    assert result["bounded_quiescence"]["attempts"] == 1


def test_same_boot_closure_boundedly_rechecks_all_quiescence(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    boot_id = tmp_path / "boot_id"
    boot_id.write_text("01234567-89ab-cdef-0123-456789abcdef\n", encoding="ascii")
    proc_stat = tmp_path / "stat"
    proc_stat.write_text("btime 1\n", encoding="ascii")
    monkeypatch.setattr(RECOVERY, "BOOT_ID_PATH", boot_id)
    monkeypatch.setattr(RECOVERY, "PROC_STAT_PATH", proc_stat)
    calls: dict[int, int] = {}

    def observe(pid: int, ticks: int) -> dict[str, Any]:
        calls[pid] = calls.get(pid, 0) + 1
        gone = calls[pid] > 1
        return {
            "pid": pid, "recorded_start_time_ticks": ticks,
            "observed_start_time_ticks": None if gone else ticks,
            "state": "disappeared" if gone else "still_live",
            "original_identity_gone": gone,
        }

    monkeypatch.setattr(RECOVERY, "_process_identity_observation", observe)
    group_calls = [0]

    def killpg(pgid: int, number: int) -> None:
        group_calls[0] += 1
        if group_calls[0] > 1:
            raise ProcessLookupError(pgid)

    monkeypatch.setattr(RECOVERY.os, "killpg", killpg)
    clock = [0.0]
    monkeypatch.setattr(RECOVERY.time, "monotonic", lambda: clock[0])
    monkeypatch.setattr(
        RECOVERY.time, "sleep", lambda seconds: clock.__setitem__(0, clock[0] + seconds))
    launch = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "hostname": socket.gethostname(),
        "mpirun_pid": 10, "mpirun_start_time_ticks": 100,
        "managed_process_group": {"pgid": 10, "new_session": True,
                                  "failure_cleanup":
                                  "SIGTERM_then_SIGKILL_with_quiescence_proof"},
        "input_transport": {"holder_pid": 11,
                            "holder_start_time_ticks": 101},
        "ranks": [{"global_rank": 0, "pid": 12,
                   "start_time_ticks": 102}],
    }

    result = RECOVERY._same_boot_closure(launch)

    assert result["bounded_quiescence"]["attempts"] == 2
    assert group_calls[0] == 2
    assert calls == {10: 2, 11: 2, 12: 2}


def test_tool_tamper_is_rejected(tmp_path: Path) -> None:
    tool = tmp_path / "recover"
    _write_closed(tool, b"version one")
    planned = {
        "path": str(tool.resolve()), "size": tool.stat().st_size,
        "sha256": hashlib.sha256(b"version one").hexdigest(),
    }
    plan = {"tools": {"prefix_recovery": planned}}
    tool.chmod(0o600)
    tool.write_bytes(b"version two")
    tool.chmod(0o400)

    with pytest.raises(RECOVERY.RecoveryFailure, match="bytes differ"):
        RECOVERY._planned_file(plan, "prefix_recovery", tool)


def test_recovery_directories_must_be_owner_only_and_empty(tmp_path: Path) -> None:
    directory = tmp_path / "recovery"
    directory.mkdir(mode=0o755)
    directory.chmod(0o755)

    with pytest.raises(RECOVERY.RecoveryFailure, match="mode 0700"):
        RECOVERY._private_empty_directory(directory, "recovery")

    directory.chmod(0o700)
    assert RECOVERY._private_empty_directory(directory, "recovery") == \
        directory.resolve()
    (directory / "old").write_text("x", encoding="ascii")
    with pytest.raises(RECOVERY.RecoveryFailure, match="must be empty"):
        RECOVERY._private_empty_directory(directory, "recovery")


def test_zero_completion_cannot_fall_back_to_same_boot_recovery(
        tmp_path: Path) -> None:
    completion = tmp_path / "segment.completion.ready"
    completion.write_text(json.dumps({
        "schema": RECOVERY.SCHEMA,
        "kind": "athenak_segment_completion",
        "status": "ready",
        "return_code": 0,
    }), encoding="utf-8")
    completion.chmod(0o400)
    plan = {"launch_contract": {"evidence": {
        "completion_record": str(completion.absolute()),
    }}}

    with pytest.raises(RECOVERY.RecoveryFailure,
                       match="requires a nonzero completion"):
        RECOVERY._lifecycle(plan, {}, completion)


def _write_live_nvidia_smi(path: Path, ranks: int) -> None:
    gpu_lines = "\n".join(
        f"{rank}, GPU-{rank}, 00000000:{16 + rank:02X}:00.0, 0, 0, "
        f"32768, {5 + rank}"
        for rank in range(ranks))
    path.write_text(
        "#!/bin/sh\n"
        "case \"$1\" in\n"
        f"  --query-gpu=*) printf '%s\\n' '{gpu_lines}' ;;\n"
        "  --query-compute-apps=*) exit 0 ;;\n"
        "  *) exit 2 ;;\n"
        "esac\n",
        encoding="ascii")
    path.chmod(0o755)


def _write_planned_binary(path: Path, metadata: Any,
                          output: dict[str, Any], write: dict[str, Any]) -> None:
    blocks = ("mesh", "meshblock", "mesh_refinement", output["block"])
    header_lines: list[str] = []
    for block in blocks:
        header_lines.append(f"<{block}>")
        parameters = (output["parameters"] if block == output["block"] else
                      metadata.parameters[block])
        header_lines.extend(f"{name} = {value}"
                            for name, value in parameters.items())
    header_lines.append("<par_end>")
    parameter_header = ("\n".join(header_lines) + "\n").encode("utf-8")
    variables = output["expected_binary_variables"]
    content = bytearray(
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        + f"  time={write['time']!r}\n".encode("ascii")
        + f"  cycle={write['cycle']}\n".encode("ascii")
        + b"  size of location=8\n"
        + b"  size of variable=4\n"
        + f"  number of variables={len(variables)}\n".encode("ascii")
        + f"  variables:  {'  '.join(variables)}  \n".encode("ascii")
        + f"  header offset={len(parameter_header)}\n".encode("ascii")
        + parameter_header)
    mesh = metadata.parameters["mesh"]
    meshblock = metadata.parameters["meshblock"]
    mesh_cells = tuple(int(mesh[f"nx{axis}"]) for axis in range(1, 4))
    block_cells = tuple(int(meshblock[f"nx{axis}"]) for axis in range(1, 4))
    root_blocks = tuple(mesh_cells[index] // block_cells[index]
                        for index in range(3))
    minima = tuple(float(mesh[f"x{axis}min"]) for axis in range(1, 4))
    maxima = tuple(float(mesh[f"x{axis}max"]) for axis in range(1, 4))
    for location in metadata.locations:
        logical = (location.lx1, location.lx2, location.lx3,
                   location.level - metadata.root_level)
        nghost = int(mesh["nghost"])
        content.extend(struct.pack(
            "=6i", *(value for extent in block_cells
                     for value in (nghost, nghost + extent - 1))))
        content.extend(struct.pack("=4i", *logical))
        geometry: list[float] = []
        for axis, coordinate in enumerate(logical[:3]):
            subdivisions = root_blocks[axis] << logical[3]
            width = maxima[axis] - minima[axis]
            geometry.extend((
                minima[axis] + width * coordinate / subdivisions,
                minima[axis] + width * (coordinate + 1) / subdivisions,
            ))
        content.extend(struct.pack("=6d", *geometry))
        content.extend(b"\0" * (
            math.prod(block_cells) * 4 * len(variables)))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    path.chmod(0o400)


def _run_log_bytes(expected: dict[str, Any], *,
                   include_driver_limit_state: bool = True,
                   limit_state: str | None = None) -> bytes:
    current = float(expected["source_time"])
    rows = []
    for cycle in range(expected["source_cycle"], expected["final_cycle"] + 1):
        if cycle != expected["source_cycle"]:
            current += float(expected["root_dt"])
        rows.append(
            f"elapsed={float(cycle - expected['source_cycle']):.6e} "
            f"cycle={cycle} time={current:.6e} "
            f"dt={expected['root_dt']:.6e}")
    cache = (
        "Strict-subcycling restart cache reconstruction passed: solver failures=0, "
        "non-finite proposed values=0, maximum raw component-relative proposed "
        "change=1.0e-13, maximum absolute proposed change=1.0e-16, maximum "
        "mixed-scale proposed change=1.0e-14, mixed-scale acceptance "
        "tolerance=1.0e-12.")
    limit = (limit_state if limit_state is not None else
             f"tlim={expected['tlim']:.6e} nlim={expected['final_cycle']}")
    prefix = f"{limit}\n" if include_driver_limit_state else ""
    return (prefix + f"{cache}\n" + "\n".join(rows) + "\n").encode("ascii")


def _interrupted_limit_context(expected: dict[str, Any], lifecycle_kind: str,
                               *, nlim: int | None = None,
                               tlim: float | None = None) \
        -> tuple[dict[str, Any], dict[str, Any]]:
    launch_binding = {
        "path": "/evidence/segment.launch.ready",
        "device": 1, "inode": 2, "size": 3, "mtime_ns": 4, "ctime_ns": 5,
        "sha256": "a" * 64, "closure_check": "linux_proc_fd",
    }
    launch_audit = {
        "sha256": launch_binding["sha256"],
        "athena_argv": [
            "/proc/123/fd/205", "-r", "/proc/123/fd/200",
            f"time/nlim={expected['final_cycle'] if nlim is None else nlim}",
            f"time/tlim={expected['tlim'] if tlim is None else tlim!r}",
        ],
    }
    if lifecycle_kind == "same_boot_closed_processes_v1":
        lifecycle = {
            "kind": lifecycle_kind,
            "all_original_identities_gone": True,
            "managed_process_group_gone": True,
            "bounded_quiescence": {"all_quiet": True},
        }
    else:
        lifecycle = {
            "kind": "nonzero_completion_v1", "return_code": 9,
            "completion_record": {"path": "/evidence/segment.completion.ready"},
            "canonical_completion_audit": {"return_code": 9},
            "artifacts": {"launch_record": launch_binding},
        }
    return lifecycle, {"launch_audit": launch_audit,
                       "launch_binding": launch_binding}


@pytest.mark.parametrize("lifecycle_kind", [
    "same_boot_closed_processes_v1", "nonzero_completion_v1",
])
def test_interrupted_run_log_without_driver_limit_uses_exact_launch_argv(
        tmp_path: Path, lifecycle_kind: str) -> None:
    expected = {
        "source_cycle": 30, "source_time": 144.0, "root_dt": 4.8,
        "final_cycle": 40, "tlim": 240.0, "recovery_final_cycle": 32,
    }
    path = tmp_path / f"{lifecycle_kind}.log"
    _write_closed(path, _run_log_bytes(
        expected, include_driver_limit_state=False))
    lifecycle, launch = _interrupted_limit_context(expected, lifecycle_kind)

    result = RECOVERY._audit_run_log_prefix(
        path, expected, lifecycle=lifecycle, **launch)

    evidence = result["original_limit_evidence"]
    assert evidence["kind"] == "immutable_plan_launch_argv_v1"
    assert evidence["lifecycle_kind"] == lifecycle_kind
    assert evidence["plan_expected"] == {"final_cycle": 40, "tlim": 240.0}
    assert evidence["exact_plan_limit_token_match"] is True


def test_missing_driver_limit_rejects_launch_argv_plan_mismatch(
        tmp_path: Path) -> None:
    expected = {
        "source_cycle": 30, "source_time": 144.0, "root_dt": 4.8,
        "final_cycle": 40, "tlim": 240.0, "recovery_final_cycle": 32,
    }
    path = tmp_path / "run.log"
    _write_closed(path, _run_log_bytes(
        expected, include_driver_limit_state=False))
    lifecycle, launch = _interrupted_limit_context(
        expected, "same_boot_closed_processes_v1", nlim=41)

    with pytest.raises(RECOVERY.RecoveryFailure,
                       match="actual launch argv does not exactly bind"):
        RECOVERY._audit_run_log_prefix(
            path, expected, lifecycle=lifecycle, **launch)


@pytest.mark.parametrize("limit_state, message", [
    ("tlim=2.400000e+02 nlim=41", "contradicts"),
    ("tlim=2.400000e+02 nlim=40\ntlim=2.400000e+02 nlim=40", "multiple"),
])
def test_present_wrong_or_duplicate_driver_limit_never_uses_fallback(
        tmp_path: Path, limit_state: str, message: str) -> None:
    expected = {
        "source_cycle": 30, "source_time": 144.0, "root_dt": 4.8,
        "final_cycle": 40, "tlim": 240.0, "recovery_final_cycle": 32,
    }
    path = tmp_path / f"{message}.log"
    _write_closed(path, _run_log_bytes(expected, limit_state=limit_state))
    lifecycle, launch = _interrupted_limit_context(
        expected, "same_boot_closed_processes_v1")

    with pytest.raises(RECOVERY.RecoveryFailure, match=message):
        RECOVERY._audit_run_log_prefix(
            path, expected, lifecycle=lifecycle, **launch)


@pytest.fixture(scope="module", params=(
    "nonzero_completion", "same_boot_sigterm",
))
def strict_recovery_chain(tmp_path_factory: pytest.TempPathFactory,
                          request: pytest.FixtureRequest):
    import test_suite.test_bbh_segment_launch as launch_tests
    import test_suite.test_bbh_segment_plan as plan_tests

    lifecycle_mode = str(request.param)
    tmp_path = tmp_path_factory.mktemp(f"strict-prefix-recovery-{lifecycle_mode}")
    campaign = plan_tests._campaign(tmp_path)
    _write_live_nvidia_smi(campaign["nvidia_smi"], 8)
    planned = plan_tests._run(campaign, root_steps="8")
    assert planned.returncode == 0, planned.stderr
    plan_path = campaign["output"]
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    expected = RECOVERY.CHECKER.validate_plan(plan)

    fake_run = launch_tests.FakeRun(8)
    prepare_runtime = launch_tests.LAUNCHER.Runtime(
        run=fake_run,
        nvidia_smi=plan["tools"]["nvidia_smi"]["path"],
        statvfs=launch_tests._ample_statvfs,
        fstatvfs=launch_tests._ample_statvfs)
    prepared = launch_tests.LAUNCHER.prepare_launch(
        plan_path, campaign["state_dir"], prepare_runtime)
    runtime = launch_tests._proof_runtime({"run": fake_run}, prepared)
    process = launch_tests.FakeProcess(fake_run, return_code=9)

    signal_handlers: dict[int, Any] | None = None
    if lifecycle_mode == "same_boot_sigterm":
        launch_hostname = socket.gethostname()
        runtime.hostname = lambda: launch_hostname
        runtime.inspector.environment = lambda pid: launch_tests._rank_environment(
            prepared, fake_run.rank_pids.index(pid), hostname=launch_hostname)
        signal_handlers = {
            number: object()
            for number in launch_tests.LAUNCHER.MANAGED_TERMINATION_SIGNALS
        }
        runtime.get_signal_handler = lambda number: signal_handlers[number]

        def set_signal_handler(number: int, handler: Any) -> Any:
            old = signal_handlers[number]
            signal_handlers[number] = handler
            return old

        runtime.set_signal_handler = set_signal_handler
        original_wait = process.wait
        injected = False

        def sigterm_wait(timeout=None):
            nonlocal injected
            if not injected and timeout is None:
                injected = True
                signal_handlers[signal.SIGTERM](signal.SIGTERM, None)
            return original_wait(timeout=timeout)

        process.wait = sigterm_wait

    def fake_popen(argv, **kwargs):
        fake_run.launched = True
        kwargs["stdout"].write(_run_log_bytes(
            expected, include_driver_limit_state=False))
        return process

    runtime.popen = fake_popen
    evidence = {name: Path(value) for name, value in
                plan["launch_contract"]["evidence"].items()}
    launch_paths = launch_tests.LAUNCHER.OutputPaths(*(
        evidence[name] for name in launch_tests.LAUNCHER.EVIDENCE_NAMES))
    if lifecycle_mode == "nonzero_completion":
        assert launch_tests.LAUNCHER.run_segment(
            prepared, plan_path, launch_paths, runtime,
            proof_timeout_seconds=1) == 9
    else:
        with pytest.raises(launch_tests.LAUNCHER.LaunchFailure,
                           match="managed termination signal SIGTERM"):
            launch_tests.LAUNCHER.run_segment(
                prepared, plan_path, launch_paths, runtime,
                proof_timeout_seconds=1)
        assert signal_handlers is not None
        assert all(not callable(handler) for handler in signal_handlers.values())
        assert not evidence["completion_record"].exists()

    recovery_cycle = expected["source_cycle"] + 4
    schedules, states, recovery_time = RECOVERY.replay_execute_prefix(
        plan, recovery_cycle)
    runtime_trajectory = next(
        token.split("=", 1)[1] for token in prepared.athena_argv
        if token.startswith("problem/trajectory_file="))
    restart_output = next(output for output in plan["outputs"]
                          if output["file_type"] == "rst")
    selected_write = schedules[restart_output["block"]][-1]
    selected_relative = restart_output["relative_path_template"].format(
        file_number=selected_write["file_number"])
    selected_path = campaign["state_dir"] / selected_relative
    selected_path.parent.mkdir(parents=True)
    plan_tests._write_restart(
        selected_path, Path(runtime_trajectory), cycle=recovery_cycle,
        time_value=recovery_time, output_states=states,
        nlim=expected["final_cycle"], tlim=expected["tlim"])
    selected_path.chmod(0o400)

    later_write = next(
        write for write in restart_output["expected_writes"]
        if write["kind"] == "scheduled" and write["cycle"] > recovery_cycle)
    later_relative = restart_output["relative_path_template"].format(
        file_number=later_write["file_number"])
    later_path = campaign["state_dir"] / later_relative
    _write_closed(later_path, b"interrupted restart header")

    source_metadata = RECOVERY.read_restart_metadata(campaign["restart"])
    for output in plan["outputs"]:
        if not output["numbered"] or output["file_type"] == "rst":
            continue
        for write in schedules[output["block"]]:
            relative = output["relative_path_template"].format(
                file_number=write["file_number"])
            _write_planned_binary(
                campaign["state_dir"] / relative,
                source_metadata, output, write)

    times: list[tuple[int, float]] = []
    current = float(expected["source_time"])
    for cycle in range(expected["source_cycle"] + 1,
                       expected["final_cycle"] + 1):
        current += float(expected["root_dt"])
        times.append((cycle, current))
    for output in plan["outputs"]:
        if output["numbered"]:
            continue
        for relative in output["required_unnumbered_paths"]:
            path = campaign["state_dir"] / relative
            if output["file_type"] == "hst":
                if "user" in path.name:
                    header = ("# Athena++ history data\n"
                              "# [1]=time [2]=dt [3]=baryon_m "
                              "[4]=bh1_x [5]=bh1_y [6]=bh1_z "
                              "[7]=bh2_x [8]=bh2_y [9]=bh2_z "
                              "[10]=bh1_mass [11]=bh2_mass\n")
                    rows = [f"{time_value!r} {expected['root_dt']!r} 100.0 "
                            "-0.25 0 0 0.25 0 0 1 1"
                            for _, time_value in times]
                else:
                    header = "# Athena++ history data\n# [1]=time [2]=dt\n"
                    rows = [f"{time_value!r} {expected['root_dt']!r}"
                            for _, time_value in times]
                path.write_text(header + "\n".join(rows) + "\n", encoding="utf-8")
            else:
                header = "# " + " ".join(RECOVERY.CHECKER.EVENT_COLUMNS) + "\n"
                rows = [f"{cycle} 0 0 0 0 0 1 0 0 0 1 1"
                        for cycle, _ in times]
                path.write_text(header + "\n".join(rows) + "\n", encoding="utf-8")
            path.chmod(0o400)

    recovery_state = tmp_path / "recovery-state"
    recovery_evidence = tmp_path / "recovery-evidence"
    recovery_state.mkdir(mode=0o700)
    recovery_evidence.mkdir(mode=0o700)
    recovery_state.chmod(0o700)
    recovery_evidence.chmod(0o700)
    pass_path = recovery_evidence / "segment.prefix.pass.ready"
    original_time = RECOVERY.time.time
    RECOVERY.time.time = lambda: original_time() + 1000.0
    original_identity = RECOVERY._process_identity_gone
    original_observation = RECOVERY._process_identity_observation
    original_killpg = RECOVERY.os.killpg
    RECOVERY._process_identity_gone = lambda pid, ticks: {
        "pid": pid, "recorded_start_time_ticks": ticks,
        "observed_start_time_ticks": None,
        "state": "disappeared", "original_identity_gone": True,
    }
    if lifecycle_mode == "same_boot_sigterm":
        RECOVERY._process_identity_observation = RECOVERY._process_identity_gone
        RECOVERY.os.killpg = lambda pgid, number: (
            (_ for _ in ()).throw(ProcessLookupError(pgid)))
    try:
        report = RECOVERY.recover(SimpleNamespace(
            plan=plan_path, state_dir=campaign["state_dir"],
            launch_record=evidence["launch_record"],
            completion_record=evidence["completion_record"],
            recovery_state_dir=recovery_state,
            recovery_evidence_dir=recovery_evidence, output=pass_path))
        bound = RECOVERY.BoundDirectory.open(
            recovery_evidence, "recovery evidence directory",
            required_mode=0o700)
        try:
            RECOVERY._install_bound_json(bound, pass_path.name, report)
        finally:
            bound.close()
        status = 0
    finally:
        RECOVERY.os.killpg = original_killpg
        RECOVERY._process_identity_observation = original_observation
        RECOVERY._process_identity_gone = original_identity
        RECOVERY.time.time = original_time
    assert status == 0
    parent = json.loads(pass_path.read_text(encoding="utf-8"))
    record_path = recovery_evidence / "segment.prefix.recovery.ready"
    record = json.loads(record_path.read_text(encoding="utf-8"))
    history = next(
        row for row in parent["output_inventory"]
        if isinstance(row.get("history"), dict) and
        "baryon_m" in row["history"]["columns"])
    return {
        "campaign": campaign, "plan": plan, "plan_path": plan_path,
        "parent": parent, "pass_path": pass_path, "record": record,
        "record_path": record_path, "selected_path": selected_path,
        "later_path": later_path, "recovery_cycle": recovery_cycle,
        "recovery_time": recovery_time, "history": history,
        "recovery_state": recovery_state, "lifecycle_mode": lifecycle_mode,
    }


def test_full_interrupted_recovery_without_finalize_limit_state(
        strict_recovery_chain) -> None:
    chain = strict_recovery_chain
    evidence = chain["record"]["run_log_prefix_audit"][
        "original_limit_evidence"]

    assert evidence["kind"] == "immutable_plan_launch_argv_v1"
    assert evidence["driver_limit_state_rows"] == 0
    assert evidence["lifecycle_kind"] == (
        "nonzero_completion_v1"
        if chain["lifecycle_mode"] == "nonzero_completion"
        else "same_boot_closed_processes_v1")


def test_real_recovery_producer_commit_checker_and_planner_consumer(
        strict_recovery_chain, tmp_path: Path) -> None:
    import test_suite.test_bbh_segment_plan as plan_tests

    chain = strict_recovery_chain
    parent = chain["parent"]
    source = RECOVERY.CHECKER._hash_record(chain["selected_path"])
    history = RECOVERY.CHECKER._hash_record(Path(chain["history"]["path"]))
    result = RECOVERY.CHECKER.audit_parent_qualification_provenance(
        parent,
        {"parent_qualification_mode": RECOVERY.QUALIFICATION_MODE},
        source, history, chain["recovery_cycle"], chain["recovery_time"],
        chain["pass_path"])
    assert result["qualification_mode"] == RECOVERY.QUALIFICATION_MODE
    assert chain["record"]["status"] == "prepared"
    assert chain["pass_path"].stat().st_mode & 0o222 == 0

    campaign = dict(chain["campaign"])
    campaign["restart"] = chain["selected_path"]
    campaign["source_history"] = Path(chain["history"]["path"])
    campaign["state_dir"] = tmp_path / "next-state"
    campaign["staging_dir"] = tmp_path / "next-staging"
    campaign["evidence_dir"] = tmp_path / "next-evidence"
    for key in ("state_dir", "staging_dir", "evidence_dir"):
        campaign[key].mkdir()
    campaign["output"] = campaign["evidence_dir"] / "segment.plan.json"
    planned = plan_tests._run_with_parent(campaign, chain["pass_path"])
    assert planned.returncode == 0, planned.stderr
    child = json.loads(campaign["output"].read_text(encoding="utf-8"))
    assert child["source_qualification"]["parent_qualification_mode"] == \
        RECOVERY.QUALIFICATION_MODE


def test_real_plan_rejects_a_later_complete_invalid_candidate(
        strict_recovery_chain, monkeypatch: pytest.MonkeyPatch) -> None:
    chain = strict_recovery_chain
    record = json.loads(json.dumps(chain["record"]))
    later_relative = chain["later_path"].relative_to(
        Path(chain["plan"]["launch_contract"]["state_dir"])).as_posix()
    later_row = next(row for row in record["candidate_inventory"]
                     if row["relative_path"] == later_relative)
    later_row["classification"] = "complete_invalid"
    original_classifier = RECOVERY.CHECKER.classify_recovery_restart

    def classify(path: Path) -> dict[str, Any]:
        if path == chain["later_path"]:
            return {
                "classification": "complete_invalid",
                "binding": RECOVERY.CHECKER._hash_record(path),
                "failure": "non-finite stored Real",
            }
        return original_classifier(path)

    monkeypatch.setattr(RECOVERY.CHECKER, "classify_recovery_restart", classify)
    with pytest.raises(
            RECOVERY.CHECKER.CheckFailure,
            match="later complete-invalid restart forbids recovery fallback"):
        RECOVERY.CHECKER._audit_recovery_candidates(
            chain["plan"], record,
            chain["parent"]["recovery_provenance"]
            ["selected_scheduled_restart"])
