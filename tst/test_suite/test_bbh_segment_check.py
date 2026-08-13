"""Fail-closed regressions for the BBH campaign segment qualifier."""

from __future__ import annotations

from dataclasses import replace
import copy
import importlib.util
import json
import math
import os
from pathlib import Path
import pwd
import stat
import struct
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
from read_athenak_restart_metadata import (  # noqa: E402
    LogicalLocation,
    RestartEventCounters,
)
CHECKER_PATH = SCRIPTS / "check_athenak_segment.py"
SPEC = importlib.util.spec_from_file_location("athenak_segment_check", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
CHECKER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CHECKER
SPEC.loader.exec_module(CHECKER)


EXPECTED = {
    "source_cycle": 10,
    "final_cycle": 13,
    "source_time": 48.0,
    "final_time": 62.39999999999999,
    "root_steps": 3,
    "root_dt": 4.8,
    "tlim": 67.19999999999999,
    "ranks": 2,
    "gpu_exit_memory_mib_max": 100.0,
    "endpoint_time_ulp_tolerance": 0,
}
CACHE_LINE = (
    "Strict-subcycling restart cache reconstruction passed: solver failures=0, "
    "non-finite proposed values=0, maximum raw component-relative proposed "
    "change=2.0e-13, maximum absolute proposed change=3.0e-16, maximum "
    "mixed-scale proposed change=1.0e-14, mixed-scale acceptance "
    "tolerance=9.0e-13."
)
EVENT_HEADER = "#  " + " ".join(CHECKER.EVENT_COLUMNS)
THRESHOLDS = [{
    "numerator": "fofc",
    "denominator": "fofc_tests",
    "max_ratio": 1.0,
}]
ABSOLUTE_THRESHOLDS = {
    "hard_equal_zero": ["eos_fail", "eos_vceil"],
    "c2p_iterations_exclusive_max": 25,
}


def _environment_contract() -> dict[str, object]:
    values = {
        "HOME": str(Path(pwd.getpwuid(os.geteuid()).pw_dir).resolve(strict=True)),
        "LANG": "C", "LC_ALL": "C", "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
    }
    return {
        "kind": CHECKER.LAUNCH_ENVIRONMENT_KIND,
        "values": values,
        "sha256": CHECKER._canonical_json_sha256(values),
    }


def _plan() -> dict[str, object]:
    return {
        "schema": 1,
        "expected": {
            "source_cycle": 10, "source_time": 48.0, "root_steps": 3,
            "root_dt": 4.8, "tlim_guard_steps": 1, "final_cycle": 13,
            "final_time": 62.39999999999999,
            "tlim": 67.19999999999999,
        },
        "command_overrides": {
            "time/nlim": 13, "time/tlim": 67.19999999999999,
            "output3/dt": 4.8,
        },
        "inputs": {
            "repo": {"path": "/campaign/repo", "commit": "a" * 40,
                     "clean": True},
            "binary": {"path": "/campaign/athena", "size": 1,
                       "sha256": "1" * 64},
            "trajectory": {"path": "/campaign/trajectory.dat", "size": 1,
                           "sha256": "2" * 64},
            "source_restart": {
                "path": "/campaign/source.rst", "size": 1,
                "sha256": "0" * 64,
            },
        },
        "tools": {
            "hash_algorithm": "sha256",
            "planner": {"path": "/campaign/plan.py", "size": 1,
                        "sha256": "3" * 64},
            "segment_checker": {"path": "/campaign/check.py", "size": 1,
                                "sha256": "7" * 64},
            "segment_launcher": {"path": "/campaign/launch.py", "size": 1,
                                 "sha256": "8" * 64},
            "restart_metadata_reader": {
                "path": "/campaign/read.py", "size": 1, "sha256": "4" * 64,
            },
            "output_integrity": {
                "path": "/campaign/integrity.py", "size": 1,
                "sha256": "5" * 64,
            },
            "restart_auditor": {
                "path": "/campaign/audit.py", "size": 1, "sha256": "9" * 64,
            },
            "restart_layout": {
                "path": "/campaign/layout.py", "size": 1, "sha256": "a" * 64,
            },
            "git": {"path": "/usr/bin/git", "size": 1,
                    "sha256": "c" * 64},
            "nvidia_smi": {"path": "/usr/bin/nvidia-smi", "size": 1,
                           "sha256": "d" * 64},
        },
        "source_qualification": {
            "mode": "anchor_full_audit",
            "audit": {
                "valid": True, "sha256": "0" * 64,
                "stored_reals": {"count": 1, "finite_count": 1,
                                  "nonfinite_count": 0},
                "layout": {"expected_file_size": 1},
                "topology": {"complete_leaf_coverage": True},
            },
            "source_baryon_mass": {
                "value": 100.0, "time": 48.0,
                "evidence": {"path": "/campaign/source.user.hst", "size": 1,
                             "sha256": "6" * 64},
            },
        },
        "policy": {
            "ranks": 2, "endpoint_time_ulp_tolerance": 0,
            "segment_termination": {
                "primary": "cycle_limit", "endpoint_time_max_ulps": 0,
                "tlim_role": "nonbinding_guard",
                "require_exact_final_cycle": True,
                "require_exact_binary64_endpoint": True,
            },
            "fixed_root_dt_mode": {
                "enabled": True, "root_dt": 4.8,
                "root_dt_parameter": "time/root_dt_max",
                "require_source_last_dt_exact": True,
                "sequential_binary64_addition": True,
            },
            "capacity": {"ranks": 2,
                         "minimum_headroom_blocks_hard": 64,
                         "minimum_headroom_blocks_yellow": 128},
            "gpu_exit_memory_mib_max": 100.0,
            "gpu_ecc": {
                "corrected_before_max": 0, "corrected_after_max": 0,
                "uncorrected_before_max": 0, "uncorrected_after_max": 0,
            },
            "mutable_parameters": CHECKER.CANONICAL_MUTABLE_PARAMETERS,
            "event_thresholds": [{
                "name": "fofc_per_test", "numerator": "fofc",
                "denominator": "fofc_tests", "max_ratio": 0.005,
            }, {
                "name": "cons_adjust_per_c2p_call",
                "numerator": "cons_adjust", "denominator": "c2p_calls",
                "max_ratio": 0.005,
            }],
            "event_absolute_thresholds": ABSOLUTE_THRESHOLDS,
            "yellow_event_thresholds": [{
                "name": "fofc_per_test", "numerator": "fofc",
                "denominator": "fofc_tests", "max_ratio": 0.001,
                "consecutive_rows": 3,
            }, {
                "name": "cons_adjust_per_c2p_call",
                "numerator": "cons_adjust", "denominator": "c2p_calls",
                "max_ratio": 0.001, "consecutive_rows": 3,
            }],
            "nonfinite_count_max": 0,
            "divb_max_abs": {"hard": 1.0e-8, "yellow": 1.0e-11},
            "baryon_mass_fractional_loss": {
                "hard_per_root_step": 0.005,
                "yellow_per_root_step": 0.0025,
                "yellow_per_48M": 0.02,
            },
            "restart_contract": {
                "real_bytes": 8, "subcycling": "level",
                "allow_legacy_subcycling_restart": False,
                "allow_legacy_ghost_event_counters": False,
                "amr_counter_version": 1, "event_counter_version": 2,
                "restart_cache_contract_version": 1,
                "pending_event_counters_all_zero": True,
                "level_subcycling_costs_exact": True,
            },
            "output_integrity": {
                "minimum_closed_file_age_seconds": 120.0,
                "require_no_open_descriptors": True,
                "require_stable_size_mtime_while_hashing": True,
                "require_sha256": True,
                "refuse_manifest_overwrite": True,
            },
            "remote_disk": {
                "warn_percent": 65,
                "do_not_start_percent": 75,
                "synchronized_stop_percent": 80,
                "emergency_stop_percent": 85,
                "minimum_reserve_gib": 50,
                "minimum_reserve_restart_multiples": 2,
            },
        },
        "launch_contract": {
            "launcher": {"path": "/campaign/mpirun", "size": 1,
                         "sha256": "b" * 64},
            "mpi_argv": ["/campaign/mpirun", "--allow-run-as-root", "--bind-to",
                         "none", "-np", "2"],
            "athena_argv_template": [
                "/proc/{holder_pid}/fd/205",
                "--kokkos-map-device-id-by=mpi_rank",
                "-r", "/proc/{holder_pid}/fd/200", "-d",
                "/proc/{holder_pid}/fd/202",
                "-t", "01:00:00", "time/nlim=13",
                "time/tlim=67.19999999999999",
                "problem/trajectory_file=/proc/{holder_pid}/fd/201",
                "output3/dt=4.8",
            ],
            "state_dir": "/campaign/state", "wall_time_seconds": 3600,
            "plan_path": "/campaign/evidence/segment.plan.json",
            "environment": _environment_contract(),
            "directory_transport": {
                "kind": CHECKER.DIRECTORY_TRANSPORT_KIND,
                "holder_pid_token": CHECKER.HOLDER_PID_TOKEN,
                "roles": {
                    "state_dir": {
                        "role": "state_dir", "fd": 202,
                        "planned_path": "/campaign/state",
                        "proc_path_template": "/proc/{holder_pid}/fd/202",
                    },
                    "evidence_dir": {
                        "role": "evidence_dir", "fd": 203,
                        "planned_path": "/campaign/evidence",
                        "proc_path_template": "/proc/{holder_pid}/fd/203",
                    },
                },
            },
            "executable_transport": {
                "kind": CHECKER.EXECUTABLE_TRANSPORT_KIND,
                "holder_pid_token": CHECKER.HOLDER_PID_TOKEN,
                "roles": {
                    "launcher": {
                        "role": "launcher", "fd": 204,
                        "parent_path": "/campaign/mpirun",
                        "proc_path_template": "/proc/{holder_pid}/fd/204",
                    },
                    "binary": {
                        "role": "binary", "fd": 205,
                        "parent_path": "/campaign/athena",
                        "proc_path_template": "/proc/{holder_pid}/fd/205",
                    },
                },
            },
            "input_transport": {
                "kind": CHECKER.INPUT_TRANSPORT_KIND,
                "holder_pid_token": CHECKER.HOLDER_PID_TOKEN,
                "staging_dir": "/campaign/staging",
                "roles": {
                    "source_restart": {
                        "role": "source_restart", "fd": 200,
                        "proc_path_template": "/proc/{holder_pid}/fd/200",
                        "staged_file": {
                            "path": "/campaign/staging/source.rst", "size": 1,
                            "sha256": "0" * 64,
                        },
                    },
                    "trajectory": {
                        "role": "trajectory", "fd": 201,
                        "proc_path_template": "/proc/{holder_pid}/fd/201",
                        "staged_file": {
                            "path": "/campaign/staging/trajectory.dat", "size": 1,
                            "sha256": "2" * 64,
                        },
                    },
                },
                "parent_content_identity": {
                    "source_restart_sha256": "0" * 64,
                    "trajectory_sha256": "2" * 64,
                    "source_serialized_trajectory_path":
                        "/campaign/trajectory.dat",
                },
                "trajectory_rebinding": {
                    "parameter": "problem/trajectory_file",
                    "parent_sha256": "2" * 64,
                    "runtime_value_template": "/proc/{holder_pid}/fd/201",
                },
            },
            "evidence_dir": "/campaign/evidence",
            "evidence": {
                "launch_record": "/campaign/evidence/segment.launch.ready",
                "completion_record": "/campaign/evidence/segment.completion.ready",
                "run_log": "/campaign/evidence/run.log",
                "exit_status": "/campaign/evidence/exit.status",
                "gpu_before": "/campaign/evidence/gpu-before.csv",
                "gpu_after": "/campaign/evidence/gpu-after.csv",
            },
            "world_size": 2, "gpu_count": 2,
        },
        "source": {"parameters": {
            "job": {"basename": "run"},
            "mesh": {"nghost": "0", "nx1": "4", "x1min": "-1",
                     "x1max": "1", "nx2": "4", "x2min": "-1",
                     "x2max": "1", "nx3": "4", "x3min": "-1",
                     "x3max": "1"},
            "meshblock": {"nx1": "2", "nx2": "2", "nx3": "2"},
            "mesh_refinement": {"num_levels": "2"},
            "mhd": {"eos": "ideal", "nscalars": "0"},
            "adm": {"dynamic": "true"},
            "problem": {"trajectory_file": "/campaign/trajectory.dat"},
            "output1": {"file_type": "hst", "dcycle": "1"},
            "output2": {"file_type": "bin", "variable": "mhd_divb",
                        "dt": "4.8", "ghost_zones": "false"},
            "output3": {"file_type": "bin", "variable": "mhd_divb",
                        "dt": "4.8", "ghost_zones": "false"},
            "output4": {"file_type": "rst", "dt": "4.8"},
            "output5": {"file_type": "log", "dcycle": "1"},
        }},
        "outputs": [{
            "block": "output1", "index": 1, "file_type": "hst",
            "enabled": True, "cadence_mode": "dcycle", "cadence": 1,
            "numbered": False,
            "relative_path_template": "run.user.hst",
            "required_unnumbered_paths": ["run.user.hst"],
            "inspect_binary": False,
            "parameters": {"file_type": "hst", "dcycle": "1"},
            "expected_writes": [
                {"cycle": 11, "time": 52.8, "kind": "scheduled"},
                {"cycle": 12, "time": 57.599999999999994,
                 "kind": "scheduled"},
                {"cycle": 13, "time": 62.39999999999999,
                 "kind": "scheduled"},
            ],
            "expected_endpoint_state": {
                "file_number": 0, "last_time": 52.8, "last_write_cycle": 13,
            },
        }, {
            "block": "output2", "index": 2, "file_type": "bin",
            "enabled": True, "cadence_mode": "dt", "cadence": 4.8,
            "numbered": True,
            "relative_path_template": "bin/run.divB.{file_number:05d}.bin",
            "required_unnumbered_paths": [], "inspect_binary": True,
            "parameters": {"file_type": "bin", "variable": "mhd_divb",
                           "dt": "4.8", "ghost_zones": "false"},
            "expected_binary_variables": ["divb"],
            "expected_writes": [
                {"cycle": 11, "time": 52.8, "kind": "scheduled",
                 "file_number": 0},
                {"cycle": 12, "time": 57.599999999999994,
                 "kind": "scheduled", "file_number": 1},
                {"cycle": 13, "time": 62.39999999999999,
                 "kind": "scheduled", "file_number": 2},
            ],
            "expected_endpoint_state": {
                "file_number": 3, "last_time": 62.39999999999999,
                "last_write_cycle": 13,
            },
        }, {
            "block": "output3", "index": 3, "file_type": "bin",
            "enabled": True, "cadence_mode": "dt", "cadence": 4.8,
            "numbered": True,
            "relative_path_template": "bin/run.topology.{file_number:05d}.bin",
            "required_unnumbered_paths": [], "inspect_binary": True,
            "parameters": {"file_type": "bin", "variable": "mhd_divb",
                           "dt": "4.8", "ghost_zones": "false"},
            "expected_binary_variables": ["divb"],
            "expected_writes": [
                {"cycle": 11, "time": 52.8, "kind": "scheduled",
                 "file_number": 0},
                {"cycle": 12, "time": 57.599999999999994,
                 "kind": "scheduled", "file_number": 1},
                {"cycle": 13, "time": 62.39999999999999,
                 "kind": "scheduled", "file_number": 2},
            ],
            "expected_endpoint_state": {
                "file_number": 3, "last_time": 62.39999999999999,
                "last_write_cycle": 13,
            },
        }, {
            "block": "output4", "index": 4, "file_type": "rst",
            "enabled": True, "cadence_mode": "dt", "cadence": 4.8,
            "numbered": True,
            "relative_path_template": "rst/run.{file_number:05d}.rst",
            "required_unnumbered_paths": [], "inspect_binary": False,
            "parameters": {"file_type": "rst", "dt": "4.8"},
            "expected_writes": [
                {"cycle": 11, "time": 52.8, "kind": "scheduled",
                 "file_number": 0},
                {"cycle": 12, "time": 57.599999999999994,
                 "kind": "scheduled", "file_number": 1},
                {"cycle": 13, "time": 62.39999999999999,
                 "kind": "scheduled", "file_number": 2},
            ],
            "expected_endpoint_state": {
                "file_number": 3, "last_time": 62.39999999999999,
                "last_write_cycle": 13,
            },
        }, {
            "block": "output5", "index": 5, "file_type": "log",
            "enabled": True, "cadence_mode": "dcycle", "cadence": 1,
            "numbered": False,
            "relative_path_template": "events.log",
            "required_unnumbered_paths": ["events.log"],
            "inspect_binary": False,
            "parameters": {"file_type": "log", "dcycle": "1"},
            "expected_writes": [
                {"cycle": 11, "time": 52.8, "kind": "scheduled"},
                {"cycle": 12, "time": 57.599999999999994,
                 "kind": "scheduled"},
                {"cycle": 13, "time": 62.39999999999999,
                 "kind": "scheduled"},
            ],
            "expected_endpoint_state": {
                "file_number": 0, "last_time": 52.8,
                "last_write_cycle": 13,
            },
        }],
        "required_files": ["run.user.hst", "events.log"],
    }


def _run_log(path: Path, termination: str = "cycle") -> None:
    diagnostics = "\n".join(
        f"elapsed={float(cycle - 9):.6e} cycle={cycle} "
        f"time={time_value:.6e} dt=4.800000e+00"
        for cycle, time_value in zip((10, 11, 12, 13),
                                     (48.0, 52.8, 57.6, 62.4))
    )
    path.write_text(
        f"{CACHE_LINE}\n"
        f"{diagnostics}\n"
        f"Terminating on {termination} limit\n"
        "time=62.399999999999999 cycle=13\n"
        "tlim=67.200000000000003 nlim=13\n",
        encoding="utf-8",
    )


def test_checker_consumes_planner_contract_schema() -> None:
    assert CHECKER.validate_plan(_plan()) == EXPECTED


@pytest.mark.parametrize(
    ("token", "message"),
    (
        ("3600", "canonical"),
        ("1:00:00", "canonical"),
        ("01:60:00", "canonical"),
    ),
)
def test_validate_plan_rejects_noncanonical_wall_time_tokens(
        token: str, message: str) -> None:
    plan = _plan()
    plan["launch_contract"]["athena_argv_template"][7] = token
    with pytest.raises(CHECKER.CheckFailure, match=message):
        CHECKER.validate_plan(plan)


def test_validate_plan_requires_exactly_one_kokkos_rank_mapping_token() -> None:
    plan = _plan()
    plan["launch_contract"]["athena_argv_template"].remove(
        "--kokkos-map-device-id-by=mpi_rank")
    with pytest.raises(CHECKER.CheckFailure, match="canonical"):
        CHECKER.validate_plan(plan)


def test_validate_plan_requires_every_execution_tool_binding() -> None:
    plan = _plan()
    del plan["tools"]["restart_auditor"]
    with pytest.raises(CHECKER.CheckFailure, match="restart_auditor"):
        CHECKER.validate_plan(plan)


@pytest.mark.parametrize(
    ("mutate", "message"),
    (
        (lambda plan: plan["tools"].pop("nvidia_smi"), "nvidia_smi"),
        (lambda plan: plan["launch_contract"].__setitem__(
            "plan_path", "/campaign/evidence/other.json"), "plan_path"),
        (lambda plan: plan["launch_contract"]["environment"]["values"].__setitem__(
            "LANG", "en_US.UTF-8"), "locale"),
        (lambda plan: plan["launch_contract"]["environment"].__setitem__(
            "sha256", "f" * 64), "SHA-256"),
        (lambda plan: plan["launch_contract"]["environment"].__setitem__(
            "unreviewed", True), "environment kind"),
        (lambda plan: plan["launch_contract"]["directory_transport"]["roles"][
            "state_dir"].__setitem__("fd", 999), "directory transport"),
        (lambda plan: plan["launch_contract"]["directory_transport"].__setitem__(
            "unreviewed", True), "fixed directory descriptors"),
        (lambda plan: plan["launch_contract"]["executable_transport"]["roles"][
            "binary"].__setitem__("fd", 204), "executable transport"),
        (lambda plan: plan["launch_contract"]["executable_transport"].__setitem__(
            "unreviewed", True), "executable descriptors"),
    ),
)
def test_validate_plan_rejects_tampered_execution_provenance(
        mutate, message: str) -> None:
    plan = _plan()
    mutate(plan)
    with pytest.raises(CHECKER.CheckFailure, match=message):
        CHECKER.validate_plan(plan)


@pytest.mark.parametrize("field", ("final_time", "tlim"))
def test_validate_plan_recomputes_binary64_endpoints(field: str) -> None:
    plan = _plan()
    # These human-friendly decimal endpoints differ from repeated binary64 addition.
    plan["expected"][field] = 62.4 if field == "final_time" else 67.2
    if field == "tlim":
        plan["command_overrides"]["time/tlim"] = 67.2
    with pytest.raises(CHECKER.CheckFailure, match="sequential binary64"):
        CHECKER.validate_plan(plan)


def _restart_metadata(**changes):
    metadata = CHECKER.RestartMetadata(
        source="fixture.rst", file_size=4096, parameter_end=512,
        metadata_end=1024, max_read_offset=1024, byte_order="little",
        real_bytes=8, nmb_total=2, root_level=0, time=48.0, last_dt=4.8,
        cycle=10,
        locations=(LogicalLocation(0, 0, 0, 0), LogicalLocation(1, 0, 0, 0)),
        costs=(1.0, 1.0), amr_cycle_counters=(0, 1),
        event_counters=RestartEventCounters(*(0 for _ in range(11))),
        event_counter_version=2, restart_cache_contract_version=1,
        parameters={
            "time": {
                "subcycling": "level",
                "allow_legacy_subcycling_restart": "false",
                "allow_legacy_ghost_event_counters": "false",
            },
            "mesh_refinement": {"max_nmb_per_rank": "65"},
        },
    )
    return replace(metadata, **changes)


def _contract_plan() -> dict[str, object]:
    plan = _plan()
    plan["policy"]["capacity"] = {
        "ranks": 2, "minimum_headroom_blocks_hard": 64,
    }
    return plan


def test_restart_contract_reproves_fixed_dt_real8_extensions_and_capacity() -> None:
    result = CHECKER.audit_restart_contract(
        _restart_metadata(), "endpoint", _contract_plan(), EXPECTED)
    assert result["last_dt"] == 4.8
    assert result["minimum_capacity_headroom"] == 64
    assert result["pending_event_counters_all_zero"] is True


@pytest.mark.parametrize(
    ("changes", "pattern"),
    (
        ({"last_dt": 2.4}, "last_dt"),
        ({"real_bytes": 4}, "Real8"),
        ({"event_counter_version": 1}, "event-counter version"),
        ({"restart_cache_contract_version": None}, "cache-contract version"),
        ({"amr_cycle_counters": (0, -1)}, "negative AMR counters"),
        ({"costs": (1.0, 2.0)}, "exact level costs"),
    ),
)
def test_restart_contract_rejects_weakened_endpoint(changes, pattern) -> None:
    with pytest.raises(CHECKER.CheckFailure, match=pattern):
        CHECKER.audit_restart_contract(
            _restart_metadata(**changes), "endpoint", _contract_plan(), EXPECTED)


def test_restart_contract_rejects_insufficient_repartitioned_headroom() -> None:
    metadata = _restart_metadata(parameters={
        **_restart_metadata().parameters,
        "mesh_refinement": {"max_nmb_per_rank": "64"},
    })
    with pytest.raises(CHECKER.CheckFailure, match="headroom 63"):
        CHECKER.audit_restart_contract(
            metadata, "endpoint", _contract_plan(), EXPECTED)


def _event_row(cycle: int) -> str:
    # All values are nonnegative, neither hard failure counter fires, and both
    # scientifically meaningful denominators are positive.
    return f"{cycle} 1 2 3 0 0 11 4 5 6 100 20"


def _event_log(path: Path, cycles: list[int]) -> None:
    path.write_text(
        EVENT_HEADER + "\n" + "\n".join(_event_row(cycle) for cycle in cycles) + "\n",
        encoding="utf-8",
    )


def _gpu(path: Path, *, uuid1: str = "GPU-b", ecc: int = 0,
         four_columns: bool = False) -> None:
    rows = [
        f"0,GPU-a,0000:01:00.0,0,0,{ecc},5",
        f"1,{uuid1},0000:02:00.0,1,0,0,7",
    ]
    if four_columns:
        rows[0] = "0,GPU-a,0,5"
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def test_exit_zero_cannot_turn_wall_limit_into_pass(tmp_path: Path) -> None:
    status = tmp_path / "exit.status"
    status.write_text("0\n", encoding="ascii")
    run_log = tmp_path / "run.log"
    _run_log(run_log, termination="wall clock")

    assert CHECKER._parse_exit_status(status) == 0
    with pytest.raises(CHECKER.CheckFailure, match="cycle-limit termination"):
        CHECKER.audit_run_log(run_log, EXPECTED)


def test_run_log_proves_every_fixed_root_step(tmp_path: Path) -> None:
    run_log = tmp_path / "run.log"
    _run_log(run_log)
    result = CHECKER.audit_run_log(run_log, EXPECTED)
    assert result["root_step_diagnostics"] == {
        "rows": 4,
        "cycle_min": 10,
        "cycle_max": 13,
        "fixed_dt": 4.8,
        "all_root_cycles_present": True,
    }


@pytest.mark.parametrize(
    ("old", "new", "message"),
    (
        ("cycle=12 time=5.760000e+01", "cycle=12 time=5.760001e+01",
         "times do not follow"),
        ("cycle=12 time=5.760000e+01 dt=4.800000e+00",
         "cycle=12 time=5.760000e+01 dt=2.400000e+00", "non-fixed"),
        ("elapsed=3.000000e+00 cycle=12 time=5.760000e+01 dt=4.800000e+00\n",
         "", "every root cycle"),
    ),
)
def test_run_log_rejects_changed_or_missing_root_step(
        tmp_path: Path, old: str, new: str, message: str) -> None:
    run_log = tmp_path / "run.log"
    _run_log(run_log)
    run_log.write_text(run_log.read_text(encoding="utf-8").replace(old, new),
                       encoding="utf-8")
    with pytest.raises(CHECKER.CheckFailure, match=message):
        CHECKER.audit_run_log(run_log, EXPECTED)


def test_one_step_run_log_has_source_and_finalize_endpoint_diagnostics(
        tmp_path: Path) -> None:
    expected = {
        **EXPECTED, "final_cycle": 11, "root_steps": 1,
        "final_time": 52.8, "tlim": 57.6,
    }
    run_log = tmp_path / "one-step.log"
    run_log.write_text(
        f"{CACHE_LINE}\n"
        "elapsed=1.000000e+00 cycle=10 time=4.800000e+01 dt=4.800000e+00\n"
        "elapsed=2.000000e+00 cycle=11 time=5.280000e+01 dt=4.800000e+00\n"
        "Terminating on cycle limit\n"
        "time=5.280000e+01 cycle=11\n"
        "tlim=5.760000e+01 nlim=11\n",
        encoding="utf-8",
    )
    result = CHECKER.audit_run_log(run_log, expected)
    assert result["root_step_diagnostics"]["rows"] == 2
    run_log.write_text(run_log.read_text(encoding="utf-8").replace(
        "elapsed=2.000000e+00 cycle=11 time=5.280000e+01 dt=4.800000e+00\n",
        ""), encoding="utf-8")
    with pytest.raises(CHECKER.CheckFailure, match="Finalize endpoint"):
        CHECKER.audit_run_log(run_log, expected)


@pytest.mark.parametrize(
    "cycles",
    (
        [11, 13],       # missing
        [11, 12, 12, 13],  # duplicate
        [11, 13, 12],   # out of order
    ),
    ids=("missing", "duplicate", "out-of-order"),
)
def test_event_log_requires_every_cycle_exactly_once(
        tmp_path: Path, cycles: list[int]) -> None:
    event_log = tmp_path / "events.log"
    _event_log(event_log, cycles)
    with pytest.raises(CHECKER.CheckFailure, match="ordered, unique"):
        CHECKER.audit_event_log(
            event_log, 10, 13, THRESHOLDS, ABSOLUTE_THRESHOLDS)


def test_four_column_gpu_snapshot_is_rejected(tmp_path: Path) -> None:
    snapshot = tmp_path / "gpu.csv"
    _gpu(snapshot, four_columns=True)
    with pytest.raises(CHECKER.CheckFailure, match="exactly 7 columns"):
        CHECKER.parse_gpu_csv(snapshot)


def test_nonzero_gpu_ecc_is_rejected(tmp_path: Path) -> None:
    before, after = tmp_path / "before.csv", tmp_path / "after.csv"
    _gpu(before)
    _gpu(after, ecc=1)
    with pytest.raises(CHECKER.CheckFailure, match="ECC"):
        CHECKER.audit_gpus(before, after, 2, 100.0)


def test_changed_gpu_uuid_is_rejected(tmp_path: Path) -> None:
    before, after = tmp_path / "before.csv", tmp_path / "after.csv"
    _gpu(before)
    _gpu(after, uuid1="GPU-replaced")
    with pytest.raises(CHECKER.CheckFailure, match="UUID"):
        CHECKER.audit_gpus(before, after, 2, 100.0)


def _launch_record_fixture(tmp_path: Path):
    plan = _plan()
    state_dir = tmp_path / "state"
    state_dir.mkdir()
    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    staging_dir = tmp_path / "staging"
    staging_dir.mkdir()
    plan_path = evidence_dir / "segment.plan.json"
    launcher = tmp_path / "mpirun"
    launcher.write_bytes(b"launcher")
    launcher.chmod(0o755)
    binary = tmp_path / "athena"
    binary.write_bytes(b"athena")
    binary.chmod(0o755)
    source = tmp_path / "source.rst"
    source.write_bytes(b"restart")
    trajectory = tmp_path / "trajectory.dat"
    trajectory.write_bytes(b"trajectory")

    def planned(path: Path) -> dict[str, object]:
        binding = CHECKER.stable_sha256(path)
        return {name: binding[name] for name in ("path", "size", "sha256")}

    launcher_record = planned(launcher)
    plan["inputs"]["binary"] = planned(binary)
    plan["inputs"]["source_restart"] = planned(source)
    plan["inputs"]["trajectory"] = planned(trajectory)
    plan["source"]["parameters"]["problem"]["trajectory_file"] = str(
        trajectory.resolve())
    plan["source_qualification"]["audit"]["sha256"] = (
        plan["inputs"]["source_restart"]["sha256"])
    plan["source_qualification"]["audit"]["layout"]["expected_file_size"] = (
        plan["inputs"]["source_restart"]["size"])

    staged_source_path = staging_dir / "source.rst"
    staged_source_path.write_bytes(source.read_bytes())
    staged_trajectory_path = staging_dir / "trajectory.dat"
    staged_trajectory_path.write_bytes(trajectory.read_bytes())
    staged_source_path.chmod(0o444)
    staged_trajectory_path.chmod(0o444)
    staged_source = CHECKER.stable_sha256(staged_source_path)
    staged_trajectory = CHECKER.stable_sha256(staged_trajectory_path)
    staging_stat = staging_dir.stat()
    staging_dir.chmod(0o555)
    staging_stat = staging_dir.stat()

    holder_pid = 400
    source_proc = f"/proc/{holder_pid}/fd/200"
    trajectory_proc = f"/proc/{holder_pid}/fd/201"
    state_proc = f"/proc/{holder_pid}/fd/202"
    mpi_argv = [str(launcher.resolve()), "--allow-run-as-root", "--bind-to",
                "none", "-np", "2"]
    athena_argv = [
        f"/proc/{holder_pid}/fd/205", "--kokkos-map-device-id-by=mpi_rank", "-r",
        source_proc, "-d", state_proc, "-t", "01:00:00", "time/nlim=13",
        "time/tlim=67.19999999999999",
        f"problem/trajectory_file={trajectory_proc}", "output3/dt=4.8",
    ]
    launcher_tool = CHECKER.stable_sha256(
        SCRIPTS / "launch_athenak_segment.py")
    plan["tools"]["segment_launcher"] = launcher_tool
    nvidia_smi = tmp_path / "nvidia-smi"
    nvidia_smi.write_bytes(b"fake nvidia-smi\n")
    nvidia_smi.chmod(0o755)
    plan["tools"]["nvidia_smi"] = CHECKER.stable_sha256(nvidia_smi)
    plan["launch_contract"]["input_transport"] = {
        "kind": CHECKER.INPUT_TRANSPORT_KIND,
        "holder_pid_token": CHECKER.HOLDER_PID_TOKEN,
        "staging_dir": str(staging_dir.resolve()),
        "roles": {
            "source_restart": {
                "role": "source_restart", "fd": 200,
                "proc_path_template": "/proc/{holder_pid}/fd/200",
                "staged_file": planned(staged_source_path),
            },
            "trajectory": {
                "role": "trajectory", "fd": 201,
                "proc_path_template": "/proc/{holder_pid}/fd/201",
                "staged_file": planned(staged_trajectory_path),
            },
        },
        "parent_content_identity": {
            "source_restart_sha256": plan["inputs"]["source_restart"]["sha256"],
            "trajectory_sha256": plan["inputs"]["trajectory"]["sha256"],
            "source_serialized_trajectory_path": str(trajectory.resolve()),
        },
        "trajectory_rebinding": {
            "parameter": "problem/trajectory_file",
            "parent_sha256": plan["inputs"]["trajectory"]["sha256"],
            "runtime_value_template": "/proc/{holder_pid}/fd/201",
        },
    }
    directory_contract = {
        "kind": CHECKER.DIRECTORY_TRANSPORT_KIND,
        "holder_pid_token": CHECKER.HOLDER_PID_TOKEN,
        "roles": {
            "state_dir": {
                "role": "state_dir", "fd": 202,
                "planned_path": str(state_dir.resolve()),
                "proc_path_template": "/proc/{holder_pid}/fd/202",
            },
            "evidence_dir": {
                "role": "evidence_dir", "fd": 203,
                "planned_path": str(evidence_dir.resolve()),
                "proc_path_template": "/proc/{holder_pid}/fd/203",
            },
        },
    }
    executable_contract = {
        "kind": CHECKER.EXECUTABLE_TRANSPORT_KIND,
        "holder_pid_token": CHECKER.HOLDER_PID_TOKEN,
        "roles": {
            "launcher": {
                "role": "launcher", "fd": 204,
                "parent_path": str(launcher.resolve()),
                "proc_path_template": "/proc/{holder_pid}/fd/204",
            },
            "binary": {
                "role": "binary", "fd": 205,
                "parent_path": str(binary.resolve()),
                "proc_path_template": "/proc/{holder_pid}/fd/205",
            },
        },
    }
    plan["launch_contract"] = {
        "wall_time_seconds": 3600, "launcher": launcher_record,
        "mpi_argv": mpi_argv,
        "athena_argv_template": [
            token.replace(str(holder_pid), "{holder_pid}") for token in athena_argv
        ],
        "state_dir": str(state_dir.resolve()), "world_size": 2, "gpu_count": 2,
        "input_transport": plan["launch_contract"]["input_transport"],
        "evidence_dir": str(evidence_dir.resolve()),
        "plan_path": str(plan_path.resolve()),
        "environment": _environment_contract(),
        "directory_transport": directory_contract,
        "executable_transport": executable_contract,
    }
    plan_path.write_text(json.dumps(plan), encoding="utf-8")
    transport = {
        "kind": CHECKER.INPUT_TRANSPORT_KIND,
        "holder_pid": holder_pid,
        "holder_start_time_ticks": 450,
        "roles": {
            "source_restart": {
                **staged_source, "role": "source_restart", "fd": 200,
                "proc_path": source_proc, "access_mode": "read_only",
            },
            "trajectory": {
                **staged_trajectory, "role": "trajectory", "fd": 201,
                "proc_path": trajectory_proc, "access_mode": "read_only",
            },
        },
    }
    def held_directory(role: str, path: Path, fd: int) -> dict[str, object]:
        info = path.stat()
        return {
            "path": str(path.resolve()), "device": info.st_dev,
            "inode": info.st_ino, "owner_uid": info.st_uid,
            "mode": f"{stat.S_IMODE(info.st_mode):04o}",
            "fd": fd, "role": role,
            "proc_path": f"/proc/{holder_pid}/fd/{fd}",
            "access_mode": "read_only_directory_descriptor",
        }

    directory_transport = {
        "kind": CHECKER.DIRECTORY_TRANSPORT_KIND,
        "holder_pid": holder_pid, "holder_start_time_ticks": 450,
        "roles": {
            "state_dir": held_directory("state_dir", state_dir, 202),
            "evidence_dir": held_directory("evidence_dir", evidence_dir, 203),
        },
    }
    def held_executable(role: str, binding: dict[str, object], fd: int) \
            -> dict[str, object]:
        canonical = {name: binding[name] for name in (
            "path", "device", "inode", "size", "mtime_ns", "ctime_ns", "sha256",
        )}
        return {
            **canonical, "closure_check": "fixed_read_only_descriptor",
            "role": role, "fd": fd,
            "proc_path": f"/proc/{holder_pid}/fd/{fd}",
            "access_mode": "read_only",
        }

    executable_transport = {
        "kind": CHECKER.EXECUTABLE_TRANSPORT_KIND,
        "holder_pid": holder_pid, "holder_start_time_ticks": 450,
        "roles": {
            "launcher": held_executable(
                "launcher", CHECKER.stable_sha256(launcher), 204),
            "binary": held_executable(
                "binary", CHECKER.stable_sha256(binary), 205),
        },
    }
    repository_proof = {
        "path": plan["inputs"]["repo"]["path"],
        "commit": plan["inputs"]["repo"]["commit"], "clean": True,
        "git_tool": plan["tools"]["git"],
        "git_environment_policy": "explicit_clean_environment_v1",
        "git_environment": dict(CHECKER.GIT_ENVIRONMENT),
        "git_configuration": list(CHECKER.GIT_CONFIGURATION),
    }
    launch = {
        "schema": 1, "kind": "athenak_segment_launch", "status": "ready",
        "mpi_argv": mpi_argv, "athena_argv": athena_argv,
        "launcher": launcher_record, "world_size": 2, "gpu_count": 2,
        "plan_path": str(plan_path.resolve()),
        "plan_sha256": CHECKER.stable_sha256(plan_path)["sha256"],
        "state_dir": str(state_dir.resolve()),
        "binary_sha256": plan["inputs"]["binary"]["sha256"],
        "source_restart_sha256": plan["inputs"]["source_restart"]["sha256"],
        "trajectory_sha256": plan["inputs"]["trajectory"]["sha256"],
        "original_inputs": {
            "source_restart": CHECKER.stable_sha256(source),
            "trajectory": CHECKER.stable_sha256(trajectory),
        },
        "staging_directory": {
            "path": str(staging_dir.resolve()), "device": staging_stat.st_dev,
            "inode": staging_stat.st_ino, "mode": "0555",
        },
        "staged_inputs": {
            "source_restart": staged_source, "trajectory": staged_trajectory,
        },
        "input_transport_contract": plan["launch_contract"]["input_transport"],
        "input_transport": transport,
        "directory_transport_contract": copy.deepcopy(directory_contract),
        "directory_transport": directory_transport,
        "executable_transport": executable_transport,
        "executable_transport_contract": copy.deepcopy(executable_contract),
        "repository_preflight": repository_proof,
        "repository_at_launch": repository_proof,
        "launch_argv": [*mpi_argv, *athena_argv],
        "mpirun_pid": 500, "mpirun_start_time_ticks": 550,
        "managed_process_group": {
            "pgid": 500, "new_session": True,
            "failure_cleanup": "SIGTERM_then_SIGKILL_with_quiescence_proof",
        },
        "proc_access_probe": {
            "kind": "fork_peer_procfs_reopen_v1", "peer_pid": 399,
            "holder_pid": holder_pid,
            "families": {
                "input": ["source_restart", "trajectory"],
                "directory": ["evidence_dir", "state_dir"],
                "executable": ["binary", "launcher"],
            },
            "executable_access": {
                "binary": "effective_uid_x_ok",
                "launcher": "effective_uid_x_ok",
            },
            "all_reopened_and_sampled": True,
        },
        "mpirun_executable": CHECKER.stable_sha256(launcher),
        "mpirun_cmdline": [*mpi_argv, *athena_argv],
        "hostname": "gpu-node",
        "ranks": [{
            "global_rank": rank, "local_rank": rank, "pid": 600 + rank,
            "start_time_ticks": 650 + rank,
            "hostname": "gpu-node", "gpu_index": rank,
            "gpu_uuid": f"GPU-{rank}",
            "gpu_pci_bus_id": f"00000000:{16 + rank:02x}:00.0",
            "gpu_cuda_ordinal": rank,
            "executable": CHECKER.stable_sha256(binary),
            "cmdline": athena_argv,
            "mpi_environment": {
                "OMPI_COMM_WORLD_RANK": str(rank),
                "OMPI_COMM_WORLD_SIZE": "2",
                "OMPI_COMM_WORLD_LOCAL_RANK": str(rank),
                "OMPI_COMM_WORLD_LOCAL_SIZE": "2",
                **_environment_contract()["values"],
            },
        } for rank in range(2)],
        "gpu_before": {
            "path": str((tmp_path / "gpu-before.csv").resolve()),
            "sha256": "c" * 64,
            "devices": [{
                "index": rank, "uuid": f"GPU-{rank}",
                "pci_bus_id": f"00000000:{16 + rank:02x}:00.0",
                "cuda_ordinal": rank,
                "uncorrected_ecc": 0, "corrected_ecc": 0, "memory_mib": 5,
            } for rank in range(2)],
        },
        "gpu_mapping_basis": (
            "kokkos_mpi_rank_token_plus_ompi_local_rank_plus_"
            "nvidia_compute_context_uuid"
        ),
        "launcher_tool": launcher_tool,
        "launch_environment": _environment_contract(),
        "nvidia_smi": plan["tools"]["nvidia_smi"],
        "execution_tools_at_launch": {
            name: plan["tools"][name] for name in (
                "segment_launcher", "segment_checker", "output_integrity",
                "restart_auditor", "restart_metadata_reader", "nvidia_smi",
            )
        },
        "gpu_visibility_environment": {
            "CUDA_VISIBLE_DEVICES": None,
            "KOKKOS_VISIBLE_DEVICES": None,
            "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
        },
    }
    launch_path = tmp_path / "segment.launch.ready"
    launch_path.write_text(json.dumps(launch), encoding="utf-8")
    launch_path.chmod(0o444)
    return plan, plan_path, state_dir, launch_path, launch


def test_launch_record_binds_exact_mpi_and_athena_argv(tmp_path: Path) -> None:
    plan, plan_path, state_dir, launch_path, launch = _launch_record_fixture(tmp_path)
    assert CHECKER.audit_launch_record(
        launch_path, plan, plan_path, state_dir)["world_size"] == 2
    launch["athena_argv"] = [*launch["athena_argv"], "output2/dt=99"]
    launch_path.chmod(0o644)
    launch_path.write_text(json.dumps(launch), encoding="utf-8")
    launch_path.chmod(0o444)
    with pytest.raises(CHECKER.CheckFailure, match="exact planned argv"):
        CHECKER.audit_launch_record(launch_path, plan, plan_path, state_dir)


def test_launch_record_requires_bound_nvidia_smi_to_remain_executable(
        tmp_path: Path) -> None:
    plan, plan_path, state_dir, launch_path, _ = _launch_record_fixture(tmp_path)
    Path(plan["tools"]["nvidia_smi"]["path"]).chmod(0o644)

    with pytest.raises(CHECKER.CheckFailure, match="nvidia-smi is not executable"):
        CHECKER.audit_launch_record(launch_path, plan, plan_path, state_dir)


@pytest.mark.parametrize(
    ("mutate", "message"),
    (
        (lambda launch: launch["launch_environment"].__setitem__(
            "sha256", "f" * 64), "environment"),
        (lambda launch: launch["nvidia_smi"].__setitem__(
            "sha256", "f" * 64), "nvidia-smi"),
        (lambda launch: launch["directory_transport"]["roles"][
            "state_dir"].__setitem__("inode", -1), "state_dir"),
        (lambda launch: launch["directory_transport"]["roles"][
            "state_dir"].__setitem__("unreviewed", True), "state_dir"),
        (lambda launch: launch["directory_transport_contract"]["roles"][
            "state_dir"].__setitem__("fd", 999), "directory-transport contract"),
        (lambda launch: launch["executable_transport"]["roles"][
            "binary"].__setitem__("sha256", "f" * 64), "binary"),
        (lambda launch: launch["executable_transport"]["roles"][
            "binary"].__setitem__("unreviewed", True), "binary"),
        (lambda launch: launch["ranks"][0]["mpi_environment"].__setitem__(
            "LANG", "C.UTF-8"), "MPI environment"),
        (lambda launch: launch["managed_process_group"].__setitem__(
            "pgid", 999), "managed process-group"),
        (lambda launch: launch["proc_access_probe"]["families"][
            "executable"].remove("binary"), "proc-holder access proof"),
    ),
)
def test_launch_record_rejects_tampered_environment_tools_and_holders(
        tmp_path: Path, mutate, message: str) -> None:
    plan, plan_path, state_dir, launch_path, launch = _launch_record_fixture(tmp_path)
    mutate(launch)
    launch_path.chmod(0o644)
    launch_path.write_text(json.dumps(launch), encoding="utf-8")
    launch_path.chmod(0o444)
    with pytest.raises(CHECKER.CheckFailure, match=message):
        CHECKER.audit_launch_record(launch_path, plan, plan_path, state_dir)


def test_launch_record_rejects_hand_written_rank_or_launcher_proof(
        tmp_path: Path) -> None:
    plan, plan_path, state_dir, launch_path, record = _launch_record_fixture(tmp_path)
    del record["mpirun_start_time_ticks"]
    launch_path.chmod(0o644)
    launch_path.write_text(json.dumps(record), encoding="utf-8")
    launch_path.chmod(0o444)
    with pytest.raises(CHECKER.CheckFailure, match="mpirun process proof"):
        CHECKER.audit_launch_record(launch_path, plan, plan_path, state_dir)


def test_nan_event_threshold_is_rejected() -> None:
    with pytest.raises(CHECKER.CheckFailure, match="finite"):
        CHECKER._normalize_event_thresholds([{
            "numerator": "fofc",
            "denominator": "fofc_tests",
            "max_ratio": float("nan"),
        }])


def test_pass_report_is_read_only_and_never_overwritten(tmp_path: Path) -> None:
    report = tmp_path / "segment.pass.ready"
    first = {"schema": 1, "kind": "athenak_segment_pass", "status": "pass"}
    CHECKER.publish_report(report, first)
    original = report.read_bytes()
    assert stat.S_IMODE(report.stat().st_mode) & 0o222 == 0

    with pytest.raises(CHECKER.CheckFailure, match="overwrite"):
        CHECKER.publish_report(report, {**first, "unexpected": True})
    assert report.read_bytes() == original
    assert json.loads(original)["status"] == "pass"


def test_recent_evidence_is_retriable(tmp_path: Path) -> None:
    recent = tmp_path / "endpoint.rst"
    recent.write_bytes(b"not inspected yet")
    with pytest.raises(CHECKER.NotReady, match="minimum_age_seconds"):
        CHECKER._check_ready([recent])


def test_baryon_hard_limit_is_applied_to_each_adjacent_step(tmp_path: Path) -> None:
    history_path = tmp_path / "run.user.hst"
    history_path.write_text(
        "# Athena++ history data\n"
        "# [1]=time [2]=dt [3]=baryon_m\n"
        "5.28e1 4.8 1.0000000000000000e2\n"
        "5.76e1 4.8 9.4000000000000000e1\n"
        "6.24e1 4.8 9.3900000000000000e1\n",
        encoding="utf-8",
    )
    history = CHECKER.audit_history(
        history_path, 10, 13, 48.0, 4.8, 62.39999999999999, 8)

    assert history["baryon_step_loss"]["maximum_fractional_loss"] == pytest.approx(0.06)
    # The total 6.1% loss is below the old N*5% segment bound (15%), but the
    # first adjacent step violates the actual 0.5% per-root-step hard limit.
    with pytest.raises(CHECKER.CheckFailure, match="adjacent-step"):
        CHECKER.audit_baryon_mass(history, 0.005, 3, 100.0)


def test_baryon_hard_limit_includes_source_to_first_step(tmp_path: Path) -> None:
    history_path = tmp_path / "run.user.hst"
    history_path.write_text(
        "# Athena++ history data\n"
        "# [1]=time [2]=dt [3]=baryon_m\n"
        "5.28e1 4.8 9.4e1\n5.76e1 4.8 9.4e1\n6.24e1 4.8 9.4e1\n",
        encoding="utf-8")
    history = CHECKER.audit_history(
        history_path, 10, 13, 48.0, 4.8, 62.39999999999999, 8)
    assert history["baryon_step_loss"]["maximum_fractional_loss"] == 0.0
    with pytest.raises(CHECKER.CheckFailure, match="adjacent-step"):
        CHECKER.audit_baryon_mass(history, 0.005, 3, 100.0)


def test_anchor_source_baryon_is_reparsed_not_trusted(tmp_path: Path) -> None:
    history = tmp_path / "source.user.hst"
    history.write_text(
        "# [1]=time [2]=dt [3]=baryon_m\n48 4.8 100\n",
        encoding="utf-8")
    audited = CHECKER.stable_sha256(history)
    qualification = {
        "mode": "anchor_full_audit",
        "source_baryon_mass": {
            "value": 101.0, "time": 48.0,
            "evidence": {key: audited[key] for key in
                         ("path", "size", "sha256", "closure_check")},
        },
    }
    with pytest.raises(CHECKER.CheckFailure, match="differs"):
        CHECKER.audit_source_baryon_evidence(qualification, 48.0)


def test_anchor_source_baryon_rejects_evidence_past_source_time(
        tmp_path: Path) -> None:
    history = tmp_path / "source.user.hst"
    history.write_text(
        "# [1]=time [2]=baryon_m\n48 100\n52.8 99\n", encoding="utf-8")
    audited = CHECKER.stable_sha256(history)
    qualification = {
        "mode": "anchor_full_audit",
        "source_baryon_mass": {
            "value": 100.0, "time": 48.0,
            "evidence": {key: audited[key] for key in
                         ("path", "size", "sha256", "closure_check")},
        },
    }
    with pytest.raises(CHECKER.CheckFailure, match="does not end"):
        CHECKER.audit_source_baryon_evidence(qualification, 48.0)


def test_parent_source_baryon_is_reloaded_not_trusted(tmp_path: Path) -> None:
    history_path = tmp_path / "parent.user.hst"
    history_path.write_text(
        "# [1]=time [2]=baryon_m\n48 99\n", encoding="utf-8")
    audited = CHECKER.stable_sha256(history_path)
    record = {key: audited[key] for key in
              ("path", "size", "sha256", "closure_check")}
    qualification = {
        "mode": "parent_segment_pass", "parent_segment_pass": record,
        "source_baryon_mass": {"value": 100.0, "time": 48.0,
                               "evidence": record},
    }
    parent = {"scientific_threshold_audit": {"baryon_mass": {"last": 100.0}}}
    with pytest.raises(CHECKER.CheckFailure, match="differs"):
        CHECKER.audit_source_baryon_evidence(qualification, 48.0, parent)


@pytest.mark.parametrize(
    ("rows", "message"),
    (
        (["5.28e1 4.8 100", "6.24e1 4.8 99"], "one row per"),
        (["5.28e1 4.8 100", "5.76e1 2.4 99", "6.24e1 4.8 98"],
         "non-fixed"),
        (["5.28e1 4.8 100", "5.77e1 4.8 99", "6.24e1 4.8 98"],
         "sequential fixed"),
    ),
)
def test_history_rejects_missing_or_variable_root_steps(
        tmp_path: Path, rows: list[str], message: str) -> None:
    history_path = tmp_path / "run.user.hst"
    history_path.write_text(
        "# Athena++ history data\n"
        "# [1]=time [2]=dt [3]=baryon_m\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    with pytest.raises(CHECKER.CheckFailure, match=message):
        CHECKER.audit_history(
            history_path, 10, 13, 48.0, 4.8, 62.39999999999999, 8)


def _write_binary(path: Path, logical_locations: list[tuple[int, int, int, int]],
                  *, variable: str | list[str] = "divb",
                  output_parameters: dict[str, str] | None = None,
                  mesh_cells: tuple[int, int, int] = (4, 4, 4),
                  block_cells: tuple[int, int, int] = (2, 2, 2),
                  nghost: int = 0) -> None:
    output_parameters = output_parameters or {
        "file_type": "bin", "variable": "mhd_divb", "dt": "4.8",
        "ghost_zones": "false",
    }
    output_lines = "\n".join(
        f"{key} = {value}" for key, value in output_parameters.items())
    parameter_header = (
        f"<mesh>\nnghost = {nghost}\nnx1 = {mesh_cells[0]}\n"
        "x1min = -1\nx1max = 1\n"
        f"nx2 = {mesh_cells[1]}\nx2min = -1\nx2max = 1\n"
        f"nx3 = {mesh_cells[2]}\nx3min = -1\nx3max = 1\n"
        f"<meshblock>\nnx1 = {block_cells[0]}\nnx2 = {block_cells[1]}\n"
        f"nx3 = {block_cells[2]}\n"
        "<mesh_refinement>\nnum_levels = 2\n"
        f"<output2>\n{output_lines}\n<par_end>\n"
    ).encode("utf-8")
    variables = [variable] if isinstance(variable, str) else variable
    content = bytearray(
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        b"  time=5.2800000000000000e+01\n"
        b"  cycle=11\n"
        b"  size of location=8\n"
        b"  size of variable=4\n"
        + f"  number of variables={len(variables)}\n".encode("ascii")
        + f"  variables:  {'  '.join(variables)}  \n".encode("ascii")
        + f"  header offset={len(parameter_header)}\n".encode("ascii")
        + parameter_header
    )
    for lx1, lx2, lx3, level in logical_locations:
        active = tuple(value > 1 for value in mesh_cells)
        extents = tuple((block_cells[axis] + 2 * nghost)
                        if active[axis] else 1 for axis in range(3))
        content.extend(struct.pack(
            "=6i", *(value for extent in extents for value in (0, extent - 1))))
        content.extend(struct.pack("=4i", lx1, lx2, lx3, level))
        root_blocks = tuple(mesh_cells[axis] // block_cells[axis]
                            for axis in range(3))
        geometry: list[float] = []
        for axis, logical in enumerate((lx1, lx2, lx3)):
            subdivisions = (root_blocks[axis] << level
                            if active[axis] else root_blocks[axis])
            geometry.extend((
                -1.0 + 2.0 * logical / subdivisions,
                -1.0 + 2.0 * (logical + 1) / subdivisions,
            ))
        content.extend(struct.pack("=6d", *geometry))
        cells = math.prod(extents)
        for _ in variables:
            content.extend(struct.pack(
                f"={cells}f", *(float(index) for index in range(cells))))
    path.write_bytes(content)


def _root_binary_locations() -> list[tuple[int, int, int, int]]:
    return [(lx1, lx2, lx3, 0)
            for lx3 in range(2) for lx2 in range(2) for lx1 in range(2)]


def test_binary_audit_proves_complete_amr_leaf_coverage(tmp_path: Path) -> None:
    path = tmp_path / "divb.bin"
    parameters = {
        "file_type": "bin", "variable": "mhd_divb", "dt": "4.8",
        "ghost_zones": "false",
    }
    _write_binary(path, _root_binary_locations(), output_parameters=parameters)
    result = CHECKER.audit_binary(
        path, 10, 13, 48.0, 62.39999999999999, parameters)
    assert result["meshblocks"] == 8
    assert result["topology"]["complete_leaf_coverage"] is True
    assert result["topology"]["coverage_units"] == result["topology"]["domain_units"]


def test_binary_audit_rejects_one_whole_meshblock_removed(tmp_path: Path) -> None:
    path = tmp_path / "truncated-on-record-boundary.bin"
    parameters = {
        "file_type": "bin", "variable": "mhd_divb", "dt": "4.8",
        "ghost_zones": "false",
    }
    _write_binary(path, _root_binary_locations()[:-1], output_parameters=parameters)
    with pytest.raises(CHECKER.CheckFailure, match="leaf coverage is incomplete"):
        CHECKER.audit_binary(
            path, 10, 13, 48.0, 62.39999999999999, parameters)


def test_binary_header_binds_mesh_refinement_parameters(tmp_path: Path) -> None:
    path = tmp_path / "divb.bin"
    parameters = {
        "file_type": "bin", "variable": "mhd_divb", "dt": "4.8",
        "ghost_zones": "false",
    }
    _write_binary(path, _root_binary_locations(), output_parameters=parameters)
    source_parameters = {
        "mesh": {
            "nghost": "0", "nx1": "4", "x1min": "-1", "x1max": "1",
            "nx2": "4", "x2min": "-1", "x2max": "1",
            "nx3": "4", "x3min": "-1", "x3max": "1",
        },
        "meshblock": {"nx1": "2", "nx2": "2", "nx3": "2"},
        "mesh_refinement": {"num_levels": "3"},
        "output2": parameters,
    }
    with pytest.raises(CHECKER.CheckFailure, match="mesh_refinement"):
        CHECKER.audit_binary(
            path, 10, 13, 48.0, 62.39999999999999, parameters,
            expected_header_parameters=source_parameters, output_block="output2")


@pytest.mark.parametrize(
    ("mesh_cells", "block_cells", "locations", "expected"),
    (
        ((4, 1, 1), (2, 1, 1), [(0, 0, 0, 0), (1, 0, 0, 0)], 2),
        ((4, 4, 1), (2, 2, 1),
         [(x, y, 0, 0) for y in range(2) for x in range(2)], 4),
    ),
    ids=("one-dimensional", "two-dimensional"),
)
def test_binary_topology_accepts_inactive_dimensions_with_ghosts(
        tmp_path: Path, mesh_cells: tuple[int, int, int],
        block_cells: tuple[int, int, int],
        locations: list[tuple[int, int, int, int]], expected: int) -> None:
    path = tmp_path / "lower-dimensional.bin"
    parameters = {
        "file_type": "bin", "variable": "mhd_divb", "dt": "4.8",
        "ghost_zones": "true",
    }
    _write_binary(path, locations, output_parameters=parameters,
                  mesh_cells=mesh_cells, block_cells=block_cells, nghost=2)
    audited = CHECKER.audit_binary(
        path, 10, 13, 48.0, 62.39999999999999, parameters)
    assert audited["meshblocks"] == expected
    assert audited["topology"]["complete_leaf_coverage"] is True


def test_binary_variables_must_exactly_match_planned_semantics(
        tmp_path: Path) -> None:
    path = tmp_path / "gr.bin"
    parameters = {
        "file_type": "bin", "variable": "mhd_gr_diagnostics", "dt": "4.8",
        "ghost_zones": "false",
    }
    wrong = ["gr_bsq", "gr_lorentz", "gr_sigma", "gr_beta_inv"]
    _write_binary(path, _root_binary_locations(), variable=wrong,
                  output_parameters=parameters)
    with pytest.raises(CHECKER.CheckFailure, match="exact planned semantics"):
        CHECKER.audit_binary(
            path, 10, 13, 48.0, 62.39999999999999, parameters,
            expected_variables=CHECKER._expected_binary_variables(
                "mhd_gr_diagnostics", {}))


def test_binary_topology_rejects_shifted_non_slice_indices() -> None:
    parameters = {
        "mesh": {
            "nghost": "4", "nx1": "16", "x1min": "-1", "x1max": "1",
            "nx2": "16", "x2min": "-1", "x2max": "1",
            "nx3": "16", "x3min": "-1", "x3max": "1",
        },
        "meshblock": {"nx1": "16", "nx2": "16", "nx3": "16"},
        "mesh_refinement": {"num_levels": "1"},
    }
    output = {
        "file_type": "bin", "variable": "mhd_divb", "ghost_zones": "false",
        "region_center": "bbh_com", "region_half_width": "1",
        "region_slice_axis": "3",
    }
    records = [{
        "indices": (999, 1014, 4, 19, -500, -500),
        "logical": (0, 0, 0, 0),
        "geometry": (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0),
    }]
    with pytest.raises(CHECKER.CheckFailure, match="exact output indices"):
        CHECKER._validate_binary_topology(
            Path("shifted.bin"), parameters, output, records, 8)


def test_selected_binary_requires_full_reference_and_exact_cell_index() -> None:
    topology = {
        "scope": "selected_region_or_slice", "root_domain": [[-1, 1]] * 3,
        "meshblock_cells": [16, 16, 16], "active_dimensions": [True] * 3,
        "nghost": 4,
    }
    selected = {
        "path": "/selected.bin", "sha256": "a" * 64,
        "cycle": 11, "time": 52.8, "topology": topology,
        "_topology_records": (((4, 19, 4, 19, -500, -500),
                                (0, 0, 0, 0),
                                (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)),),
    }
    output = {
        "region_center": "bbh_com", "region_half_width": "1",
        "region_slice_axis": "3", "region_slice_offset": "0", "gid": "-1",
    }
    center = {"bbh_com": (0.0, 0.0, 0.0)}
    with pytest.raises(CHECKER.CheckFailure, match="full-domain topology"):
        CHECKER.audit_selected_binary_coverage(selected, output, center)
    full = {
        "path": "/full.bin", "sha256": "b" * 64,
        "cycle": 11, "time": 52.8,
        "topology": {"scope": "full_domain", "complete_leaf_coverage": True},
        "_topology_records": (((4, 19, 4, 19, 4, 19),
                                (0, 0, 0, 0),
                                (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)),),
    }
    with pytest.raises(CHECKER.CheckFailure, match="CellCenterIndex"):
        CHECKER.audit_selected_binary_coverage(selected, output, center, full)


def test_state_tree_rejects_extra_and_nonregular_nodes(tmp_path: Path) -> None:
    state = tmp_path / "state"
    state.mkdir()
    planned = state / "planned.dat"
    planned.write_bytes(b"planned")
    assert CHECKER._assert_exact_state_tree(state, [planned])[
        "exact_replayed_allowlist"] is True
    empty_extra = state / "unexpected-empty-directory"
    empty_extra.mkdir()
    with pytest.raises(CHECKER.CheckFailure, match="directory set"):
        CHECKER._assert_exact_state_tree(state, [planned])
    empty_extra.rmdir()
    extra = state / "extra.dat"
    extra.write_bytes(b"extra")
    with pytest.raises(CHECKER.CheckFailure, match="extra="):
        CHECKER._assert_exact_state_tree(state, [planned])
    extra.unlink()
    fifo = state / "fifo"
    import os
    os.mkfifo(fifo)
    with pytest.raises(CHECKER.CheckFailure, match="non-regular"):
        CHECKER._assert_exact_state_tree(state, [planned])


def test_binary_parameter_header_rejects_duplicate_blocks_and_keys() -> None:
    with pytest.raises(CHECKER.CheckFailure, match="duplicate binary parameter"):
        CHECKER._parse_binary_parameters(
            b"<mesh>\nnx1=4\nnx1=8\n<par_end>\n", Path("duplicate-key.bin"))
    with pytest.raises(CHECKER.CheckFailure, match="duplicate binary parameter block"):
        CHECKER._parse_binary_parameters(
            b"<mesh>\nnx1=4\n<mesh>\nnx2=4\n<par_end>\n",
            Path("duplicate-block.bin"))


def test_completion_record_binds_all_artifacts_quiescence_and_holder_closure(
        tmp_path: Path) -> None:
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    plan_path = tmp_path / "plan.json"
    plan_path.write_text("{}\n", encoding="utf-8")
    paths = {
        "plan": plan_path,
        "launch_record": evidence / "segment.launch.ready",
        "run_log": evidence / "run.log",
        "exit_status": evidence / "exit.status",
        "gpu_before": evidence / "gpu-before.csv",
        "gpu_after": evidence / "gpu-after.csv",
    }
    for name, path in paths.items():
        if name != "plan":
            path.write_text(name + "\n", encoding="utf-8")
        path.chmod(0o444)
    bindings = {name: CHECKER._hash_record(path)
                for name, path in paths.items()}
    holder = {
        "kind": CHECKER.INPUT_TRANSPORT_KIND,
        "holder_pid": 123, "holder_start_time_ticks": 456,
        "roles": {
            "source_restart": {"fd": 200, "sha256": "a" * 64},
            "trajectory": {"fd": 201, "sha256": "b" * 64},
        },
    }
    def directory_role(role: str, path: Path, fd: int) -> dict[str, object]:
        info = path.stat()
        return {
            "path": str(path.resolve()), "device": info.st_dev,
            "inode": info.st_ino, "owner_uid": info.st_uid,
            "mode": f"{stat.S_IMODE(info.st_mode):04o}",
            "fd": fd, "role": role,
            "proc_path": f"/proc/123/fd/{fd}",
            "access_mode": "read_only_directory_descriptor",
        }
    directories = {
        "kind": CHECKER.DIRECTORY_TRANSPORT_KIND,
        "holder_pid": 123, "holder_start_time_ticks": 456,
        "roles": {
            "state_dir": directory_role("state_dir", tmp_path, 202),
            "evidence_dir": directory_role("evidence_dir", evidence, 203),
        },
    }
    directory_contract = {
        "kind": CHECKER.DIRECTORY_TRANSPORT_KIND,
        "holder_pid_token": CHECKER.HOLDER_PID_TOKEN,
        "roles": {
            "state_dir": {"role": "state_dir", "fd": 202,
                          "planned_path": str(tmp_path.resolve()),
                          "proc_path_template": "/proc/{holder_pid}/fd/202"},
            "evidence_dir": {"role": "evidence_dir", "fd": 203,
                             "planned_path": str(evidence.resolve()),
                             "proc_path_template": "/proc/{holder_pid}/fd/203"},
        },
    }
    launch = {
        "mpirun_pid": 500, "mpirun_start_time_ticks": 600,
        "ranks": [{"global_rank": 0, "pid": 700,
                   "start_time_ticks": 800}],
        "input_transport": holder,
        "directory_transport": directories,
        "executable_transport": {
            "kind": CHECKER.EXECUTABLE_TRANSPORT_KIND,
            "holder_pid": 123, "holder_start_time_ticks": 456,
            "roles": {},
        },
    }
    completion = evidence / "segment.completion.ready"
    # Completion's executable proof is exercised by the launch fixture; this focused
    # unit still needs exact real-file roles for the common audit helper.
    launcher_file = tmp_path / "mpirun"
    binary_file = tmp_path / "athena"
    launcher_file.write_bytes(b"launcher")
    binary_file.write_bytes(b"binary")
    for executable in (launcher_file, binary_file):
        executable.chmod(0o755)
    executable_contract = {
        "kind": CHECKER.EXECUTABLE_TRANSPORT_KIND,
        "holder_pid_token": CHECKER.HOLDER_PID_TOKEN,
        "roles": {
            "launcher": {"role": "launcher", "fd": 204,
                         "parent_path": str(launcher_file.resolve()),
                         "proc_path_template": "/proc/{holder_pid}/fd/204"},
            "binary": {"role": "binary", "fd": 205,
                       "parent_path": str(binary_file.resolve()),
                       "proc_path_template": "/proc/{holder_pid}/fd/205"},
        },
    }
    def executable_role(role: str, path: Path, fd: int) -> dict[str, object]:
        binding = CHECKER.stable_sha256(path)
        canonical = {name: binding[name] for name in (
            "path", "device", "inode", "size", "mtime_ns", "ctime_ns", "sha256",
        )}
        return {
            **canonical, "closure_check": "fixed_read_only_descriptor",
            "role": role, "fd": fd, "proc_path": f"/proc/123/fd/{fd}",
            "access_mode": "read_only",
        }
    executables = {
        "kind": CHECKER.EXECUTABLE_TRANSPORT_KIND,
        "holder_pid": 123, "holder_start_time_ticks": 456,
        "roles": {
            "launcher": executable_role("launcher", launcher_file, 204),
            "binary": executable_role("binary", binary_file, 205),
        },
    }
    launch["executable_transport"] = executables
    payload = {
        "schema": 1, "kind": "athenak_segment_completion", "status": "ready",
        "return_code": 0, "state_dir": str(tmp_path.resolve()), "world_size": 1,
        "artifacts": copy.deepcopy(bindings),
        "quiescence": {
            "gpu_compute_contexts_empty": True,
            "all_original_identities_gone": True,
            "process_identities": [{
                "role": "mpirun", "pid": 500,
                "recorded_start_time_ticks": 600, "state": "disappeared",
                "observed_start_time_ticks": None, "original_identity_gone": True,
            }, {
                "role": "rank:0", "pid": 700,
                "recorded_start_time_ticks": 800, "state": "pid_reused",
                "observed_start_time_ticks": 801,
                "original_identity_gone": True,
            }],
        },
        "input_transport": {
            "kind": CHECKER.INPUT_TRANSPORT_KIND,
            "at_launch": holder, "at_completion": holder,
            "closure": {"all_holder_fds_closed": True, "roles": {
                "source_restart": {"fd": 200, "closed": True},
                "trajectory": {"fd": 201, "closed": True},
            }},
        },
        "directory_transport": directories,
        "executable_transport": executables,
        "publication_contract": {
            "unique_lifecycle_complete_marker": True,
            "published_after_closed_artifacts": [
                "launch_record", "run_log", "exit_status",
                "gpu_before", "gpu_after",
            ],
        },
    }
    completion.write_text(json.dumps(payload), encoding="utf-8")
    completion.chmod(0o444)
    plan = {
        "policy": {"ranks": 1},
        "inputs": {"binary": CHECKER.stable_sha256(binary_file)},
        "launch_contract": {
        "launcher": CHECKER.stable_sha256(launcher_file),
        "directory_transport": directory_contract,
        "executable_transport": executable_contract,
    }}
    audited = CHECKER.audit_completion_record(
        completion, plan, tmp_path, launch, bindings)
    assert audited["return_code"] == 0
    payload["artifacts"]["run_log"]["sha256"] = "f" * 64
    completion.chmod(0o644)
    completion.write_text(json.dumps(payload), encoding="utf-8")
    completion.chmod(0o444)
    with pytest.raises(CHECKER.CheckFailure, match="run_log binding"):
        CHECKER.audit_completion_record(completion, plan, tmp_path, launch, bindings)
    payload["artifacts"]["run_log"] = copy.deepcopy(bindings["run_log"])
    payload["directory_transport"]["roles"]["state_dir"]["inode"] = -1
    completion.chmod(0o644)
    completion.write_text(json.dumps(payload), encoding="utf-8")
    completion.chmod(0o444)
    with pytest.raises(CHECKER.CheckFailure, match="state_dir"):
        CHECKER.audit_completion_record(completion, plan, tmp_path, launch, bindings)
    payload["directory_transport"]["roles"]["state_dir"]["inode"] = (
        tmp_path.stat().st_ino)
    payload["executable_transport"]["roles"]["binary"]["sha256"] = "f" * 64
    completion.chmod(0o644)
    completion.write_text(json.dumps(payload), encoding="utf-8")
    completion.chmod(0o444)
    with pytest.raises(CHECKER.CheckFailure, match="binary"):
        CHECKER.audit_completion_record(completion, plan, tmp_path, launch, bindings)


def test_w_gr_pair_rejects_same_cycle_with_different_topology() -> None:
    shared = {"file_type": "bin", "dt": "4.8", "ghost_zones": "false"}
    plan = {"outputs": [
        {"block": "output2", "inspect_binary": True,
         "parameters": {**shared, "variable": "mhd_w_bcc"}},
        {"block": "output5", "inspect_binary": True,
         "parameters": {**shared, "variable": "mhd_gr_diagnostics"}},
    ]}
    left = {"cycle": 11, "time": 52.8, "meshblocks": 8,
            "topology_sha256": "a" * 64}
    right = {**left, "topology_sha256": "b" * 64}
    with pytest.raises(CHECKER.CheckFailure, match="indices/logical"):
        CHECKER.audit_binary_pairs(
            plan, {"output2": [left], "output5": [right]})


@pytest.mark.parametrize("mutation", ("missing", "extra", "short-padding"))
def test_numbered_output_set_rejects_missing_extra_or_noncanonical_names(
        tmp_path: Path, mutation: str) -> None:
    directory = tmp_path / "bin"
    directory.mkdir()
    template = "bin/run.divB.{file_number:05d}.bin"
    writes = [{"file_number": number} for number in (4, 5, 6)]
    for write in writes:
        (directory / f"run.divB.{write['file_number']:05d}.bin").write_bytes(b"x")
    if mutation == "missing":
        (directory / "run.divB.00005.bin").unlink()
    elif mutation == "extra":
        (directory / "run.divB.00007.bin").write_bytes(b"x")
    else:
        (directory / "run.divB.4.bin").write_bytes(b"x")
    with pytest.raises(CHECKER.CheckFailure, match="directory set"):
        CHECKER._assert_exact_numbered_output_set(
            tmp_path, template, writes, "output2")


def test_retry_exit_never_creates_pass_report(tmp_path: Path, capsys) -> None:
    plan = _plan()
    state_dir = tmp_path / "state"
    state_dir.mkdir()
    restart_dir = state_dir / "rst"
    restart_dir.mkdir()
    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    evidence = {
        "source.rst": tmp_path / "source.rst",
        "events.log": state_dir / "events.log",
        "launch": evidence_dir / "segment.launch.ready",
        "completion": evidence_dir / "segment.completion.ready",
        "run.log": evidence_dir / "run.log",
        "exit.status": evidence_dir / "exit.status",
        "gpu-before.csv": evidence_dir / "gpu-before.csv",
        "gpu-after.csv": evidence_dir / "gpu-after.csv",
    }
    evidence["endpoint.rst"] = restart_dir / "run.00002.rst"
    for path in evidence.values():
        path.write_bytes(b"recent\n")
    plan["inputs"]["source_restart"]["path"] = str(evidence["source.rst"])
    plan["launch_contract"]["state_dir"] = str(state_dir)
    plan["launch_contract"]["evidence_dir"] = str(evidence_dir)
    plan["launch_contract"]["plan_path"] = str(evidence_dir / "segment.plan.json")
    plan["launch_contract"]["directory_transport"]["roles"][
        "state_dir"]["planned_path"] = str(state_dir)
    plan["launch_contract"]["directory_transport"]["roles"][
        "evidence_dir"]["planned_path"] = str(evidence_dir)
    plan["launch_contract"]["evidence"] = {
        "launch_record": str(evidence["launch"]),
        "completion_record": str(evidence["completion"]),
        "run_log": str(evidence["run.log"]),
        "exit_status": str(evidence["exit.status"]),
        "gpu_before": str(evidence["gpu-before.csv"]),
        "gpu_after": str(evidence["gpu-after.csv"]),
    }
    plan_path = evidence_dir / "segment.plan.json"
    plan_path.write_text(json.dumps(plan), encoding="utf-8")
    report = tmp_path / "segment.pass.ready"

    status = CHECKER.main([
        "--plan", str(plan_path),
        "--launch-record", str(evidence["launch"]),
        "--completion-record", str(evidence["completion"]),
        "--endpoint-restart", str(evidence["endpoint.rst"]),
        "--state-dir", str(state_dir),
        "--run-log", str(evidence["run.log"]),
        "--event-log", str(evidence["events.log"]),
        "--exit-status", str(evidence["exit.status"]),
        "--gpu-before", str(evidence["gpu-before.csv"]),
        "--gpu-after", str(evidence["gpu-after.csv"]),
        "--output", str(report),
    ])

    captured = capsys.readouterr()
    assert status == CHECKER.RETRY_EXIT_STATUS
    assert json.loads(captured.err)["status"] == "retry"
    assert not report.exists()
