"""Regressions for the strict AthenaK segment launcher.

No test in this module starts MPI or an AthenaK executable.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import pwd
import signal
import shutil
import subprocess
import sys
from types import SimpleNamespace
from typing import Any

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
LAUNCHER_PATH = SCRIPTS / "launch_athenak_segment.py"
SPEC = importlib.util.spec_from_file_location("athenak_segment_launch", LAUNCHER_PATH)
assert SPEC is not None and SPEC.loader is not None
LAUNCHER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = LAUNCHER
SPEC.loader.exec_module(LAUNCHER)


@pytest.fixture(autouse=True)
def _close_fixed_holder_descriptors_after_test():
    """Keep independent tests from retaining the launcher's fixed process FDs."""

    yield
    for descriptor in (
            LAUNCHER.SOURCE_RESTART_FD, LAUNCHER.TRAJECTORY_FD,
            LAUNCHER.STATE_DIRECTORY_FD, LAUNCHER.EVIDENCE_DIRECTORY_FD,
            LAUNCHER.MPI_LAUNCHER_FD, LAUNCHER.BINARY_EXECUTABLE_FD,
            LAUNCHER.STAGING_DIRECTORY_PREFLIGHT_FD):
        try:
            os.close(descriptor)
        except OSError:
            pass


def _git(repo: Path, *arguments: str) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo), *arguments], check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    return result.stdout.strip()


def _record(path: Path) -> dict[str, Any]:
    content = path.read_bytes()
    return {
        "path": str(path.resolve()),
        "size": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
        "closure_check": "fixture",
    }


def _mca_contract(home: str, prefix: Path) -> dict[str, Any]:
    prefix_directory = LAUNCHER._snapshot_mca_prefix_directory(prefix)
    prefix = Path(prefix_directory["path"])
    files = [
        LAUNCHER._snapshot_mca_configuration_file(
            scope, project,
            (Path(home) if scope == "home" else prefix) / relative)
        for scope, project, relative in LAUNCHER.MCA_CONFIGURATION_LAYOUT
    ]
    payload = {
        "kind": LAUNCHER.MCA_CONFIGURATION_KIND,
        "home": home, "prefix": str(prefix),
        "prefix_directory": prefix_directory, "files": files,
    }
    return {**payload, "sha256": LAUNCHER._canonical_sha256(payload)}


def _ample_statvfs(_target: Any) -> Any:
    return SimpleNamespace(
        f_frsize=1 << 30, f_blocks=1000, f_bfree=900, f_bavail=900)


def _exactly_75_percent_used_statvfs(_target: Any) -> Any:
    return SimpleNamespace(
        f_frsize=1 << 30, f_blocks=1000, f_bfree=250, f_bavail=250)


def _disk_contract(source_size: int, trajectory_size: int = 0,
                   peak_gib: int = 200) -> dict[str, Any]:
    peak_bytes = peak_gib * LAUNCHER.GIB_BYTES
    return {
        "kind": LAUNCHER.DISK_PREFLIGHT_KIND,
        "accounting": LAUNCHER.DISK_PREFLIGHT_ACCOUNTING,
        "formula": LAUNCHER.DISK_PREFLIGHT_FORMULA,
        "used_percent_exclusive_max": 75,
        "minimum_reserve_bytes": 50 * LAUNCHER.GIB_BYTES,
        "minimum_reserve_restart_multiples": 2,
        "additional_hard_minimum_free_bytes": 180 * LAUNCHER.GIB_BYTES,
        "source_restart_size_bytes": source_size,
        "source_restart_staging_bytes": source_size,
        "trajectory_staging_bytes": trajectory_size,
        "staging_copy_bytes": source_size + trajectory_size,
        "planned_peak_output_gib": peak_gib,
        "planned_peak_output_bytes": peak_bytes,
        "role_contributions_bytes": {
            "state_dir": {"planned_peak_output_bytes": peak_bytes},
            "staging_dir": {
                "source_restart_staging_bytes": source_size,
                "trajectory_staging_bytes": trajectory_size,
            },
        },
        "bound_directory_fds": {
            "state_dir": LAUNCHER.STATE_DIRECTORY_FD,
            "staging_dir": LAUNCHER.STAGING_DIRECTORY_PREFLIGHT_FD,
        },
    }


class FakeRun:
    def __init__(self, world_size: int) -> None:
        self.world_size = world_size
        self.launched = False
        self.completed = False
        self.extra_application = False
        self.permute_applications = False
        self.reverse_nvidia_indices = False
        self.retain_applications_after_wait = False
        self.retain_process_identity_after_wait = False
        self.reuse_process_ids_after_wait = False
        self.cleanup_signals: list[int] = []
        self.rank_pids = [600 + index for index in range(world_size)]
        self.memory_total_mib = 32768
        self.memory_used_base_mib = 5

    def __call__(self, command, **kwargs):
        if Path(command[0]).name != "nvidia-smi":
            return subprocess.run(command, **kwargs)
        if command[1].startswith("--query-gpu="):
            output = "".join(
                f"{self.world_size - 1 - index if self.reverse_nvidia_indices else index}, "
                f"GPU-{index}, 00000000:{16 + index:02X}:00.0, "
                f"0, 0, {self.memory_total_mib}, "
                f"{self.memory_used_base_mib + index}\n"
                for index in range(self.world_size)
            )
        elif command[1].startswith("--query-compute-apps="):
            applications_active = self.launched or (
                self.completed and self.retain_applications_after_wait
            )
            output = "" if not applications_active else "".join(
                f"{pid}, GPU-{self.world_size - 1 - index if self.permute_applications else index}\n"
                for index, pid in enumerate(self.rank_pids)
            )
            if self.launched and self.extra_application:
                output += "999,GPU-0\n"
        else:  # pragma: no cover - catches an accidental new external dependency
            raise AssertionError(command)
        return subprocess.CompletedProcess(command, 0, output, "")


def _campaign(tmp_path: Path, world_size: int = 2) -> dict[str, Any]:
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "config", "user.email", "launch-test@example.invalid")
    _git(repo, "config", "user.name", "Launch Test")
    (repo / "tracked").write_text("clean\n", encoding="ascii")
    _git(repo, "add", "tracked")
    _git(repo, "commit", "-q", "-m", "fixture")

    inputs = tmp_path / "inputs"
    inputs.mkdir()
    binary = inputs / "athena"
    binary.write_bytes(b"fake-athena-binary\n")
    binary.chmod(0o755)
    source = inputs / "source.rst"
    source.write_bytes(b"closed-restart\n")
    trajectory = inputs / "trajectory.dat"
    trajectory.write_text("0 0 0\n", encoding="ascii")
    mpirun = inputs / "mpirun"
    mpirun.write_bytes(b"fake-openmpi-launcher\n")
    mpirun.chmod(0o755)
    nvidia_smi = inputs / "nvidia-smi"
    nvidia_smi.write_bytes(b"fake-nvidia-smi\n")
    nvidia_smi.chmod(0o755)
    state = tmp_path / "state"
    state.mkdir()
    staging = tmp_path / "staging"
    staging.mkdir(mode=0o700)
    evidence = tmp_path / "evidence"
    evidence.mkdir(mode=0o700)
    mca_prefix = tmp_path / "openmpi-prefix"
    mca_prefix.mkdir()
    staged_source = staging / "source.rst"
    staged_trajectory = staging / "trajectory.dat"
    source_record = _record(source)
    trajectory_record = _record(trajectory)
    staged_source_record = {**source_record, "path": str(staged_source.resolve())}
    staged_trajectory_record = {
        **trajectory_record, "path": str(staged_trajectory.resolve())
    }

    final_cycle = 13
    tlim = 67.19999999999999
    root_dt = 4.8
    wall = 28800
    holder = "{holder_pid}"
    athena_argv_template = [
        f"/proc/{holder}/fd/205", "--kokkos-map-device-id-by=mpi_rank",
        "-r", f"/proc/{holder}/fd/200",
        "-d", f"/proc/{holder}/fd/202", "-t", "08:00:00",
        f"time/nlim={final_cycle}", f"time/tlim={tlim!r}",
        f"problem/trajectory_file=/proc/{holder}/fd/201",
        f"output3/dt={root_dt!r}",
    ]
    mpi_argv = [
        str(mpirun.resolve()), "--allow-run-as-root", "--bind-to", "none",
        "-np", str(world_size),
    ]
    evidence_paths = {
        "launch_record": str((evidence / "segment.launch.ready").resolve()),
        "completion_record": str((evidence / "segment.completion.ready").resolve()),
        "run_log": str((evidence / "run.log").resolve()),
        "exit_status": str((evidence / "exit.status").resolve()),
        "gpu_before": str((evidence / "gpu-before.csv").resolve()),
        "gpu_after": str((evidence / "gpu-after.csv").resolve()),
    }
    launch_environment = {
        "HOME": str(Path(pwd.getpwuid(os.geteuid()).pw_dir).resolve(strict=True)),
        "LANG": "C", "LC_ALL": "C", "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
        "PRTE_MCA_schizo_proxy": "ompi",
    }
    rank_projection_payload = {
        "inherited_values": {
            key: launch_environment[key]
            for key in ("HOME", "LANG", "LC_ALL", "CUDA_DEVICE_ORDER")
        },
        "consumed_absent": ["PRTE_MCA_schizo_proxy"],
        "exact_keys": list(LAUNCHER.RANK_ENVIRONMENT_KEYS),
        "fixed_values": dict(LAUNCHER.RANK_FIXED_ENVIRONMENT_VALUES),
        "derived_values": dict(LAUNCHER.RANK_DERIVED_ENVIRONMENT_VALUES),
    }
    plan = {
        "schema": 1,
        "policy": {"ranks": world_size},
        "expected": {
            "final_cycle": final_cycle, "tlim": tlim, "root_dt": root_dt,
        },
        "command_overrides": {
            "time/nlim": final_cycle, "time/tlim": tlim,
            "output3/dt": root_dt,
        },
        "capacity_transition": {
            "kind": "unchanged_v1",
            "parameter": "mesh_refinement/max_nmb_per_rank",
            "source_max_nmb_per_rank": 1024,
            "target_max_nmb_per_rank": 1024,
            "maximum_target_max_nmb_per_rank": 16384,
            "runtime_override": None,
            "gpu_memory_model": {
                "kind": "fixed_conservative_per_meshblock_slot_v1",
                "mib_per_slot_numerator": 1433,
                "mib_per_slot_denominator": 100,
                "usable_fraction_numerator": 4,
                "usable_fraction_denominator": 5,
                "required_per_rank_memory_mib_ceiling": 14674,
                "minimum_gpu_memory_total_mib": 18343,
            },
        },
        "restart_cadence_transition": {
            "kind": "unchanged_v1", "block": "output4",
            "parameter": "output4/dt", "source_dt": 19.2,
            "target_dt": 19.2, "root_dt": 4.8,
            "source_root_step_multiple": 4,
            "target_root_step_multiple": 4,
            "phase": {"file_number": 9, "last_time": 48.0,
                      "last_write_cycle": 10},
            "runtime_override": None,
        },
        "source": {
            "parameters": {
                "problem": {"trajectory_file": str(trajectory.resolve())},
                "mesh_refinement": {"max_nmb_per_rank": "1024"},
                "output4": {
                    "file_type": "rst", "dt": "19.2", "file_number": "9",
                    "last_time": "48.0", "last_write_cycle": "10",
                },
            },
        },
        "inputs": {
            "repo": {
                "path": str(repo.resolve()),
                "commit": _git(repo, "rev-parse", "HEAD"),
                "clean": True,
            },
            "binary": _record(binary),
            "source_restart": source_record,
            "trajectory": trajectory_record,
        },
        "tools": {
            "git": _record(Path(shutil.which("git") or "").resolve(strict=True)),
            "segment_launcher": _record(LAUNCHER_PATH),
            "segment_checker": _record(SCRIPTS / "check_athenak_segment.py"),
            "output_integrity": _record(SCRIPTS / "output_integrity.py"),
            "restart_auditor": _record(SCRIPTS / "audit_athenak_restart.py"),
            "restart_metadata_reader": _record(
                SCRIPTS / "read_athenak_restart_metadata.py"),
            "nvidia_smi": _record(nvidia_smi),
            "hash_algorithm": "sha256",
        },
        "launch_contract": {
            "state_dir": str(state.resolve()),
            "wall_time_seconds": wall,
            "world_size": world_size,
            "gpu_count": world_size,
            "plan_path": str((evidence / "segment.plan.json").resolve()),
            "evidence_dir": str(evidence.resolve()),
            "evidence": evidence_paths,
            "environment": {
                "kind": LAUNCHER.LAUNCH_ENVIRONMENT_KIND,
                "values": launch_environment,
                "sha256": LAUNCHER._canonical_sha256(launch_environment),
                "rank_projection": {
                    "kind": LAUNCHER.RANK_ENVIRONMENT_PROJECTION_KIND,
                    **rank_projection_payload,
                    "sha256": LAUNCHER._canonical_sha256(
                        rank_projection_payload),
                },
            },
            "mca_configuration": _mca_contract(
                launch_environment["HOME"], mca_prefix),
            "directory_transport": {
                "kind": "linux_proc_holder_dirfd_v1",
                "holder_pid_token": holder,
                "roles": {
                    "state_dir": {
                        "role": "state_dir", "fd": 202,
                        "planned_path": str(state.resolve()),
                        "proc_path_template": f"/proc/{holder}/fd/202",
                    },
                    "evidence_dir": {
                        "role": "evidence_dir", "fd": 203,
                        "planned_path": str(evidence.resolve()),
                        "proc_path_template": f"/proc/{holder}/fd/203",
                    },
                },
            },
            "executable_transport": {
                "kind": "linux_proc_holder_execfd_v1",
                "holder_pid_token": holder,
                "roles": {
                    "launcher": {
                        "role": "launcher", "fd": 204,
                        "parent_path": str(mpirun.resolve()),
                        "proc_path_template": f"/proc/{holder}/fd/204",
                    },
                    "binary": {
                        "role": "binary", "fd": 205,
                        "parent_path": str(binary.resolve()),
                        "proc_path_template": f"/proc/{holder}/fd/205",
                    },
                },
            },
            "input_transport": {
                "kind": "linux_proc_holder_fd_v1",
                "holder_pid_token": holder,
                "staging_dir": str(staging.resolve()),
                "roles": {
                    "source_restart": {
                        "role": "source_restart", "fd": 200,
                        "proc_path_template": f"/proc/{holder}/fd/200",
                        "staged_file": staged_source_record,
                    },
                    "trajectory": {
                        "role": "trajectory", "fd": 201,
                        "proc_path_template": f"/proc/{holder}/fd/201",
                        "staged_file": staged_trajectory_record,
                    },
                },
                "parent_content_identity": {
                    "source_restart_sha256": source_record["sha256"],
                    "trajectory_sha256": trajectory_record["sha256"],
                    "source_serialized_trajectory_path": str(trajectory.resolve()),
                },
                "trajectory_rebinding": {
                    "parameter": "problem/trajectory_file",
                    "parent_sha256": trajectory_record["sha256"],
                    "runtime_value_template": f"/proc/{holder}/fd/201",
                },
            },
            "disk_preflight": _disk_contract(
                source_record["size"], trajectory_record["size"]),
            "launcher": _record(mpirun),
            "mpi_argv": mpi_argv,
            "athena_argv_template": athena_argv_template,
        },
    }
    plan_path = evidence / "segment.plan.json"
    plan_path.write_text(
        json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    plan_path.chmod(0o444)
    fake_run = FakeRun(world_size)
    runtime = LAUNCHER.Runtime(
        run=fake_run, statvfs=_ample_statvfs, fstatvfs=_ample_statvfs)
    runtime.plan_validator = lambda plan: {"fixture": "plan-validated"}
    runtime.source_validator = lambda path, plan: {
        "fixture": "source-validated",
    }
    return {
        "repo": repo, "binary": binary, "source": source,
        "trajectory": trajectory, "mpirun": mpirun, "state": state,
        "staging": staging, "staged_source": staged_source,
        "staged_trajectory": staged_trajectory,
        "evidence": evidence,
        "mca_prefix": mca_prefix,
        "plan": plan, "plan_path": plan_path, "run": fake_run,
        "runtime": runtime,
    }


def _rewrite_plan(campaign: dict[str, Any], mutate) -> None:
    path = campaign["plan_path"]
    path.chmod(0o644)
    plan = json.loads(path.read_text(encoding="utf-8"))
    mutate(plan)
    path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
    path.chmod(0o444)


def _set_capacity_increase(plan: dict[str, Any], target: int = 1280) -> None:
    required = (target * 1433 + 99) // 100
    minimum_total = (target * 1433 * 5 + 399) // 400
    token = f"mesh_refinement/max_nmb_per_rank={target}"
    plan["capacity_transition"] = {
        "kind": "increase_v1",
        "parameter": "mesh_refinement/max_nmb_per_rank",
        "source_max_nmb_per_rank": 1024,
        "target_max_nmb_per_rank": target,
        "maximum_target_max_nmb_per_rank": 16384,
        "runtime_override": token,
        "gpu_memory_model": {
            "kind": "fixed_conservative_per_meshblock_slot_v1",
            "mib_per_slot_numerator": 1433,
            "mib_per_slot_denominator": 100,
            "usable_fraction_numerator": 4,
            "usable_fraction_denominator": 5,
            "required_per_rank_memory_mib_ceiling": required,
            "minimum_gpu_memory_total_mib": minimum_total,
        },
    }
    plan["command_overrides"]["mesh_refinement/max_nmb_per_rank"] = target
    plan["launch_contract"]["athena_argv_template"].append(token)


def _set_restart_cadence_tightening(plan: dict[str, Any]) -> None:
    plan["source"]["parameters"]["output4"]["dt"] = "48.0"
    plan["restart_cadence_transition"] = {
        "kind": "tighten_v1", "block": "output4",
        "parameter": "output4/dt", "source_dt": 48.0,
        "target_dt": 19.2, "root_dt": 4.8,
        "source_root_step_multiple": 10,
        "target_root_step_multiple": 4,
        "phase": {"file_number": 9, "last_time": 48.0,
                  "last_write_cycle": 10},
        "runtime_override": "output4/dt=19.2",
    }
    plan["command_overrides"]["output4/dt"] = 19.2
    plan["launch_contract"]["athena_argv_template"].append(
        "output4/dt=19.2")


def _set_global_cadence_tightening(plan: dict[str, Any]) -> None:
    expected = plan["expected"]
    expected.update({
        "source_cycle": 10, "source_time": 48.0,
        "final_time": 62.39999999999999,
    })
    phase = {"file_number": 4, "last_time": 48.0,
             "last_write_cycle": 10}
    source_schedule = [{
        "cycle": 13, "time": 62.39999999999999,
        "kind": "forced_final", "file_number": 4,
    }]
    source_endpoint = {
        "file_number": 5, "last_time": 48.0, "last_write_cycle": 13,
    }
    target_schedule = [{
        "cycle": 13, "time": 62.39999999999999,
        "kind": "scheduled", "file_number": 4,
    }]
    target_endpoint = {
        "file_number": 5, "last_time": 58.0, "last_write_cycle": 13,
    }
    plan["source"]["parameters"]["job"] = {"basename": "run"}
    streams = []
    plan["outputs"] = []
    for block, variable in LAUNCHER.GLOBAL_CADENCE_STREAMS:
        source_parameters = {
            "file_type": "bin", "variable": variable, "dt": "48.0",
            "file_number": "4", "last_time": "48.0",
            "last_write_cycle": "10", "ghost_zones": "false",
        }
        plan["source"]["parameters"][block] = source_parameters
        runtime_parameters = {**source_parameters, "dt": "10.0"}
        streams.append({
            "block": block, "parameter": f"{block}/dt",
            "file_type": "bin", "variable": variable,
            "source_dt": 48.0, "target_dt": 10.0,
            "phase": dict(phase),
            "runtime_override": f"{block}/dt=10.0",
            "source_schedule": [dict(row) for row in source_schedule],
            "target_schedule": [dict(row) for row in target_schedule],
            "source_endpoint_state": dict(source_endpoint),
            "target_endpoint_state": dict(target_endpoint),
        })
        plan["outputs"].append({
            "block": block, "file_type": "bin", "cadence_mode": "dt",
            "cadence": 10.0, "numbered": True,
            "relative_path_template":
                f"bin/run.{variable}.{{file_number:05d}}.bin",
            "parameters": runtime_parameters,
            "parameters_sha256": LAUNCHER._canonical_sha256(
                runtime_parameters),
            "expected_writes": [dict(row) for row in target_schedule],
            "expected_endpoint_state": dict(target_endpoint),
        })
    plan["global_cadence_transition"] = {
        "kind": "tighten_pair_v1", "target_dt": 10.0,
        "streams": streams,
        "runtime_overrides": ["output2/dt=10.0", "output5/dt=10.0"],
    }
    plan["command_overrides"].update({
        "output2/dt": 10.0, "output5/dt": 10.0,
    })
    plan["launch_contract"]["athena_argv_template"].extend(
        ["output2/dt=10.0", "output5/dt=10.0"])


def test_prepare_derives_only_the_canonical_token_arrays(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)

    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])

    assert prepared.mpi_argv == tuple(
        campaign["plan"]["launch_contract"]["mpi_argv"])
    assert prepared.athena_argv == tuple(
        token.replace("{holder_pid}", str(prepared.input_holder.pid))
        for token in campaign["plan"]["launch_contract"]["athena_argv_template"]
    )
    assert prepared.launch_argv == prepared.mpi_argv + prepared.athena_argv
    assert [gpu.uuid for gpu in prepared.gpus] == ["GPU-0", "GPU-1"]
    assert prepared.plan["launch_contract"]["mca_configuration"]["prefix"] == \
        str(campaign["mca_prefix"].resolve())
    assert campaign["mpirun"].parent != campaign["mca_prefix"] / "bin"
    assert prepared.staged_source_restart["path"] == str(
        campaign["staged_source"].resolve())
    assert campaign["staged_source"].read_bytes() == campaign["source"].read_bytes()
    assert campaign["staged_trajectory"].read_bytes() == \
        campaign["trajectory"].read_bytes()
    assert campaign["staging"].stat().st_mode & 0o777 == 0o555
    assert campaign["staged_source"].stat().st_mode & 0o777 == 0o444
    assert prepared.athena_argv.count(
        f"problem/trajectory_file=/proc/{prepared.input_holder.pid}/fd/201") == 1
    assert prepared.athena_argv.count(
        f"/proc/{prepared.input_holder.pid}/fd/200") == 1
    assert prepared.athena_argv[0] == \
        f"/proc/{prepared.input_holder.pid}/fd/205"
    assert prepared.athena_argv[5] == \
        f"/proc/{prepared.input_holder.pid}/fd/202"
    assert prepared.proc_access_probe["all_reopened_and_sampled"] is True
    holder = prepared.input_holder.audit()
    assert holder["kind"] == "linux_proc_holder_fd_v1"
    assert holder["roles"]["source_restart"]["sha256"] == \
        prepared.source_restart["sha256"]


def test_disk_budget_counts_one_reserve_on_a_shared_filesystem() -> None:
    contract = _disk_contract(
        10 * LAUNCHER.GIB_BYTES, 0, peak_gib=200)
    calls: list[object] = []

    def observe(target: object) -> Any:
        calls.append(target)
        return _ample_statvfs(target)

    snapshot = LAUNCHER._disk_preflight_snapshot(
        contract, "before_spawn", {
            "state_dir": {
                "device": 7, "target": LAUNCHER.STATE_DIRECTORY_FD,
                "access": {"method": "bound_fd_fstatvfs_v1",
                           "fd": LAUNCHER.STATE_DIRECTORY_FD},
            },
            "staging_dir": {
                "device": 7,
                "target": LAUNCHER.STAGING_DIRECTORY_PREFLIGHT_FD,
                "access": {"method": "bound_fd_fstatvfs_v1",
                           "fd": LAUNCHER.STAGING_DIRECTORY_PREFLIGHT_FD},
            },
        }, observe)

    # Both paths are really measured, while the shared-device budget below
    # contains one reserve and one copy of each role contribution.
    assert len(calls) == 2
    assert set(calls) == {
        LAUNCHER.STATE_DIRECTORY_FD,
        LAUNCHER.STAGING_DIRECTORY_PREFLIGHT_FD,
    }
    filesystem = snapshot["filesystems"][0]
    assert filesystem["reserve_bytes"] == 50 * LAUNCHER.GIB_BYTES
    assert filesystem["required_free_bytes"] == 260 * LAUNCHER.GIB_BYTES


def test_disk_budget_keeps_independent_reserves_on_split_filesystems() -> None:
    contract = _disk_contract(
        10 * LAUNCHER.GIB_BYTES, 0, peak_gib=200)
    snapshot = LAUNCHER._disk_preflight_snapshot(
        contract, "before_staging", {
            "state_dir": {
                "device": 7, "target": "/state",
                "access": {"method": "path_statvfs_v1", "path": "/state"},
            },
            "staging_dir": {
                "device": 8, "target": "/staging",
                "access": {"method": "path_statvfs_v1", "path": "/staging"},
            },
        }, _ample_statvfs)

    by_role = {
        tuple(row["roles"]): row for row in snapshot["filesystems"]
    }
    assert by_role[("state_dir",)]["required_free_bytes"] == \
        250 * LAUNCHER.GIB_BYTES
    # The staging filesystem gets its own reserve and hard 180-GiB floor;
    # neither reserve is merged into or omitted from the other filesystem.
    assert by_role[("staging_dir",)]["reserve_bytes"] == \
        50 * LAUNCHER.GIB_BYTES
    assert by_role[("staging_dir",)]["required_free_bytes"] == \
        180 * LAUNCHER.GIB_BYTES


def test_first_disk_failure_precedes_staging_and_never_calls_popen(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    popen_called = False

    def forbidden_popen(*args, **kwargs):
        nonlocal popen_called
        popen_called = True
        raise AssertionError("disk-rejected launch called Popen")

    campaign["runtime"].fstatvfs = _exactly_75_percent_used_statvfs
    campaign["runtime"].popen = forbidden_popen
    with pytest.raises(LAUNCHER.LaunchFailure, match="at or above 75%"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])

    assert popen_called is False
    assert not any(campaign["staging"].iterdir())
    assert campaign["staging"].stat().st_mode & 0o777 == 0o700


def test_source_path_swap_during_descriptor_copy_is_rejected(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    campaign = _campaign(tmp_path)
    original_open = LAUNCHER._open_regular_nofollow
    backup = campaign["source"].with_suffix(".original")
    swapped = False

    def swapping_open(path):
        nonlocal swapped
        result = original_open(path)
        if Path(path) == campaign["source"] and not swapped:
            campaign["source"].rename(backup)
            campaign["source"].write_bytes(b"same-size-evil!\n")
            swapped = True
        return result

    monkeypatch.setattr(LAUNCHER, "_open_regular_nofollow", swapping_open)

    with pytest.raises(LAUNCHER.LaunchFailure, match="cannot stage.*file changed"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])

    assert not campaign["staged_source"].exists()
    assert not any(campaign["staging"].iterdir())


def test_stage_link_is_removed_if_post_link_descriptor_audit_fails(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    campaign = _campaign(tmp_path)
    original_audit = LAUNCHER._audit_open_descriptor

    def failing_audit(descriptor, planned, label):
        if label == "staged source_restart":
            raise LAUNCHER.LaunchFailure("injected post-link audit failure")
        return original_audit(descriptor, planned, label)

    monkeypatch.setattr(LAUNCHER, "_audit_open_descriptor", failing_audit)
    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="injected post-link audit failure"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])

    assert not campaign["staged_source"].exists()
    assert not any(campaign["staging"].iterdir())


def test_first_disk_gate_measures_bound_descriptors_not_paths(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    targets: list[object] = []

    def forbidden_path_statvfs(_target: object) -> Any:
        raise AssertionError("first disk gate used a path statvfs")

    def observe_fd(descriptor: object) -> Any:
        assert isinstance(descriptor, int)
        os.fstat(descriptor)
        targets.append(descriptor)
        return _ample_statvfs(descriptor)

    campaign["runtime"].statvfs = forbidden_path_statvfs
    campaign["runtime"].fstatvfs = observe_fd
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        assert len(targets) == 2
        observations = prepared.disk_preflight_before_staging["filesystems"][0][
            "observations"]
        assert observations["state_dir"]["access"]["method"] == \
            "bound_directory_fstatvfs_v1"
        assert observations["staging_dir"]["access"]["fd"] == \
            LAUNCHER.STAGING_DIRECTORY_PREFLIGHT_FD
    finally:
        prepared.close()


def test_prepare_accepts_exact_1024_to_1280_capacity_transition(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(campaign, _set_capacity_increase)

    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        assert prepared.athena_argv[-1] == \
            "mesh_refinement/max_nmb_per_rank=1280"
        assert prepared.athena_argv.count(
            "mesh_refinement/max_nmb_per_rank=1280") == 1
    finally:
        prepared.close()


def test_prepare_accepts_exact_restart_cadence_tightening(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(campaign, _set_restart_cadence_tightening)

    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        assert prepared.athena_argv[-1] == "output4/dt=19.2"
        assert prepared.athena_argv.count("output4/dt=19.2") == 1
    finally:
        prepared.close()


def test_prepare_accepts_exact_paired_global_cadence_tightening(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(campaign, _set_global_cadence_tightening)

    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        assert prepared.athena_argv[-2:] == (
            "output2/dt=10.0", "output5/dt=10.0")
        assert prepared.athena_argv.count("output2/dt=10.0") == 1
        assert prepared.athena_argv.count("output5/dt=10.0") == 1
    finally:
        prepared.close()


def test_launcher_prelaunch_rejects_duplicate_numbered_filename_templates(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)

    def mutate(plan: dict[str, Any]) -> None:
        _set_global_cadence_tightening(plan)
        source = plan["source"]["parameters"]["output5"]
        source["id"] = "mhd_w_bcc"
        output = next(row for row in plan["outputs"]
                      if row["block"] == "output5")
        output["parameters"]["id"] = "mhd_w_bcc"
        output["parameters_sha256"] = LAUNCHER._canonical_sha256(
            output["parameters"])
        output["relative_path_template"] = \
            "bin/run.mhd_w_bcc.{file_number:05d}.bin"

    _rewrite_plan(campaign, mutate)

    with pytest.raises(LAUNCHER.LaunchFailure, match="template collides"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


@pytest.mark.parametrize("file_id", (
    "unique-but-noncanonical", " mhd_w_bcc "))
def test_launcher_rejects_noncanonical_global_pair_id(
        tmp_path: Path, file_id: str) -> None:
    campaign = _campaign(tmp_path)

    def mutate(plan: dict[str, Any]) -> None:
        _set_global_cadence_tightening(plan)
        source = plan["source"]["parameters"]["output2"]
        source["id"] = file_id
        output = next(row for row in plan["outputs"]
                      if row["block"] == "output2")
        output["parameters"]["id"] = file_id
        output["parameters_sha256"] = LAUNCHER._canonical_sha256(
            output["parameters"])
        output["relative_path_template"] = (
            f"bin/run.{file_id.strip()}.{{file_number:05d}}.bin")

    _rewrite_plan(campaign, mutate)

    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="global cadence source/target/phase"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


@pytest.mark.parametrize(
    ("write_number", "endpoint_number", "message"),
    (
        (100000, 100001, "expected file_number 100000"),
        (99999, 100000, "endpoint next file_number 100000"),
    ),
)
def test_launcher_prelaunch_rejects_numbered_filename_counter_exhaustion(
        tmp_path: Path, write_number: int, endpoint_number: int,
        message: str) -> None:
    campaign = _campaign(tmp_path)

    def mutate(plan: dict[str, Any]) -> None:
        _set_global_cadence_tightening(plan)
        output = next(row for row in plan["outputs"]
                      if row["block"] == "output2")
        output["expected_writes"][0]["file_number"] = write_number
        output["expected_endpoint_state"]["file_number"] = endpoint_number

    _rewrite_plan(campaign, mutate)

    with pytest.raises(LAUNCHER.LaunchFailure, match=message):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


@pytest.mark.parametrize(
    "tamper",
    ("command", "argv", "phase", "source_schedule", "target_schedule",
     "endpoint", "pair", "one_override", "output_schedule",
     "output_endpoint", "output_parameters"),
)
def test_prepare_rejects_paired_global_cadence_transition_tamper(
        tmp_path: Path, tamper: str) -> None:
    campaign = _campaign(tmp_path)

    def mutate(plan: dict[str, Any]) -> None:
        _set_global_cadence_tightening(plan)
        transition = plan["global_cadence_transition"]
        if tamper == "command":
            plan["command_overrides"]["output2/dt"] = 9.6
        elif tamper == "argv":
            plan["launch_contract"]["athena_argv_template"][-2] = \
                "output2/dt=9.6"
        elif tamper == "phase":
            transition["streams"][0]["phase"]["last_time"] = 43.2
        elif tamper == "source_schedule":
            transition["streams"][0]["source_schedule"][0]["kind"] = \
                "scheduled"
        elif tamper == "target_schedule":
            transition["streams"][1]["target_schedule"][0]["cycle"] = 12
        elif tamper == "endpoint":
            transition["streams"][0]["target_endpoint_state"][
                "file_number"] += 1
        elif tamper == "pair":
            plan["source"]["parameters"]["output5"]["last_time"] = "43.2"
        elif tamper == "one_override":
            transition["runtime_overrides"].pop()
        elif tamper == "output_schedule":
            plan["outputs"][0]["expected_writes"][0]["cycle"] = 12
        elif tamper == "output_endpoint":
            plan["outputs"][1]["expected_endpoint_state"]["file_number"] += 1
        else:
            plan["outputs"][0]["parameters"]["dt"] = "9.6"

    _rewrite_plan(campaign, mutate)

    with pytest.raises(LAUNCHER.LaunchFailure):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_launcher_rejects_present_null_global_cadence_extension(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(
        campaign,
        lambda plan: plan.__setitem__("global_cadence_transition", None))

    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="global_cadence_transition"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


@pytest.mark.parametrize(
    "mutation",
    (
        lambda transition: transition.__setitem__("target_dt", 10),
        lambda transition: transition["streams"][0].__setitem__(
            "source_dt", 48),
        lambda transition: transition["streams"][1].__setitem__(
            "target_dt", 10),
        lambda transition: transition["streams"][0]["phase"].__setitem__(
            "file_number", 4.0),
        lambda transition: transition["streams"][0]["phase"].__setitem__(
            "last_time", 48),
        lambda transition: transition["streams"][0]["phase"].__setitem__(
            "last_write_cycle", 10.0),
        lambda transition: transition["streams"][0]["source_schedule"][0].
            __setitem__("file_number", 4.0),
        lambda transition: transition["streams"][1]["target_schedule"][0].
            __setitem__("time", 62),
        lambda transition: transition["streams"][0]["source_endpoint_state"].
            __setitem__("file_number", 5.0),
        lambda transition: transition["streams"][1]["target_endpoint_state"].
            __setitem__("last_time", 58),
    ),
    ids=("integer-target", "integer-source", "integer-stream-target",
         "float-file-number", "integer-last-time", "float-last-cycle",
         "float-write-number", "integer-write-time",
         "float-endpoint-number", "integer-endpoint-time"),
)
def test_launcher_rejects_noncanonical_global_cadence_numeric_types(
        tmp_path: Path, mutation) -> None:
    campaign = _campaign(tmp_path)

    def mutate(plan: dict[str, Any]) -> None:
        _set_global_cadence_tightening(plan)
        mutation(plan["global_cadence_transition"])

    _rewrite_plan(campaign, mutate)

    with pytest.raises(LAUNCHER.LaunchFailure, match="numeric types|numeric type"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_prepare_accepts_decimal_restart_multiple_with_division_roundoff(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)

    def mutate(plan: dict[str, Any]) -> None:
        plan["source"]["parameters"]["output4"]["dt"] = "33.6"
        plan["restart_cadence_transition"].update({
            "source_dt": 33.6, "target_dt": 33.6,
            "source_root_step_multiple": 7,
            "target_root_step_multiple": 7,
        })

    _rewrite_plan(campaign, mutate)

    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        assert "output4/dt=33.6" not in prepared.athena_argv
    finally:
        prepared.close()


@pytest.mark.parametrize("tamper", ("command", "argv", "phase"))
def test_prepare_rejects_restart_cadence_transition_tamper(
        tmp_path: Path, tamper: str) -> None:
    campaign = _campaign(tmp_path)

    def mutate(plan: dict[str, Any]) -> None:
        _set_restart_cadence_tightening(plan)
        if tamper == "command":
            plan["command_overrides"]["output4/dt"] = 9.6
        elif tamper == "argv":
            plan["launch_contract"]["athena_argv_template"][-1] = \
                "output4/dt=9.6"
        else:
            plan["restart_cadence_transition"]["phase"]["last_time"] = 52.8

    _rewrite_plan(campaign, mutate)

    with pytest.raises(LAUNCHER.LaunchFailure):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


@pytest.mark.parametrize(
    "mutation",
    (
        lambda transition: transition.__setitem__("source_dt", 19),
        lambda transition: transition.__setitem__(
            "source_root_step_multiple", 4.0),
        lambda transition: transition["phase"].__setitem__(
            "file_number", 9.0),
        lambda transition: transition["phase"].__setitem__(
            "last_time", 48),
        lambda transition: transition["phase"].__setitem__(
            "last_write_cycle", 10.0),
    ),
    ids=("integer-dt", "float-multiple", "float-file-number",
         "integer-last-time", "float-last-write-cycle"),
)
def test_launcher_rejects_noncanonical_restart_cadence_numeric_types(
        tmp_path: Path, mutation) -> None:
    campaign = _campaign(tmp_path)

    def mutate(plan: dict[str, Any]) -> None:
        mutation(plan["restart_cadence_transition"])

    _rewrite_plan(campaign, mutate)

    with pytest.raises(LAUNCHER.LaunchFailure, match="numeric types"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_prepare_rejects_low_total_memory_before_popen(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(campaign, _set_capacity_increase)
    campaign["run"].memory_total_mib = 22000

    with pytest.raises(LAUNCHER.LaunchFailure, match="GPU total memory"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_prepare_rejects_high_used_memory_without_compute_app(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(campaign, _set_capacity_increase)
    campaign["run"].memory_used_base_mib = 15000

    with pytest.raises(LAUNCHER.LaunchFailure, match="available memory"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])
    assert campaign["run"].launched is False


@pytest.mark.parametrize("token", (
    "mesh_refinement/max_nmb_per_rank=1280",
    "mesh_refinement/refinement=static",
    "mesh/nx1=32",
    "meshblock/nx1=8",
))
def test_prepare_rejects_duplicate_or_extra_mesh_override(
        tmp_path: Path, token: str) -> None:
    campaign = _campaign(tmp_path)

    def mutate(plan: dict[str, Any]) -> None:
        _set_capacity_increase(plan)
        plan["launch_contract"]["athena_argv_template"].append(token)

    _rewrite_plan(campaign, mutate)
    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="athena_argv_template is not canonical"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_preexisting_staged_target_is_never_overwritten(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    campaign["staged_source"].write_bytes(b"sentinel\n")

    with pytest.raises(LAUNCHER.LaunchFailure, match="initially be empty"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])

    assert campaign["staged_source"].read_bytes() == b"sentinel\n"


def test_prepare_only_never_calls_popen(tmp_path: Path, capsys) -> None:
    campaign = _campaign(tmp_path)

    def forbidden_popen(*args, **kwargs):  # pragma: no cover - must never run
        raise AssertionError("prepare-only attempted to launch MPI")

    campaign["runtime"].popen = forbidden_popen
    status = LAUNCHER.main([
        "--plan", str(campaign["plan_path"]),
        "--state-dir", str(campaign["state"]),
        "--prepare-only",
    ], campaign["runtime"])

    assert status == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["status"] == "prepared"
    assert summary["launch_argv"] == (
        campaign["plan"]["launch_contract"]["mpi_argv"] +
        [token.replace("{holder_pid}", str(os.getpid())) for token in
         campaign["plan"]["launch_contract"]["athena_argv_template"]]
    )
    assert not LAUNCHER._fd_is_open(LAUNCHER.SOURCE_RESTART_FD)
    assert not LAUNCHER._fd_is_open(LAUNCHER.TRAJECTORY_FD)


def test_cuda_ordinals_follow_pci_order_not_nvidia_index(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    campaign["run"].reverse_nvidia_indices = True

    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])

    assert [gpu.index for gpu in prepared.gpus] == [1, 0]
    assert [gpu.cuda_ordinal for gpu in prepared.gpus] == [0, 1]
    assert [gpu.uuid for gpu in prepared.gpus] == ["GPU-0", "GPU-1"]


def test_writable_plan_is_rejected(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    campaign["plan_path"].chmod(0o644)

    with pytest.raises(LAUNCHER.LaunchFailure, match="immutable"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_changed_source_hash_is_rejected(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    campaign["source"].write_bytes(b"tampered\n")

    with pytest.raises(LAUNCHER.LaunchFailure, match="source_restart.*differs"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_dirty_planned_repository_is_rejected(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    (campaign["repo"] / "untracked").write_text("dirty", encoding="ascii")

    with pytest.raises(LAUNCHER.LaunchFailure, match="no longer clean"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_noncanonical_or_extra_athena_override_is_rejected(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(
        campaign,
        lambda plan: plan["launch_contract"]["athena_argv_template"].append(
            "output2/dt=1.0"))

    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="athena_argv_template is not canonical"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_divb_cadence_must_be_exact_root_dt_override(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(
        campaign,
        lambda plan: plan["launch_contract"]["athena_argv_template"].__setitem__(
            11, "output3/dt=9.6"))

    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="athena_argv_template is not canonical"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_divb_dcycle_override_is_never_accepted(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(
        campaign,
        lambda plan: plan["launch_contract"]["athena_argv_template"].__setitem__(
            11, "output3/dcycle=1"))

    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="athena_argv_template is not canonical"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_wall_time_seconds_token_is_rejected_in_favor_of_hhmmss(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(
        campaign,
        lambda plan: plan["launch_contract"]["athena_argv_template"].__setitem__(
            7, "28800"))

    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="athena_argv_template is not canonical"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_noncanonical_mpi_options_are_rejected(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(campaign, lambda plan: plan["launch_contract"]["mpi_argv"].remove(
        "--bind-to"))

    with pytest.raises(LAUNCHER.LaunchFailure, match="mpi_argv is not canonical"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_missing_kokkos_mpi_rank_mapping_token_is_rejected(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(
        campaign,
        lambda plan: plan["launch_contract"]["athena_argv_template"].remove(
            "--kokkos-map-device-id-by=mpi_rank"))

    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="athena_argv_template is not canonical"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_command_overrides_must_be_the_exact_planned_limits(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    _rewrite_plan(campaign, lambda plan: plan["command_overrides"].update({
        "output2/dt": 1.0,
    }))

    with pytest.raises(LAUNCHER.LaunchFailure, match="command_overrides"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_existing_gpu_application_is_rejected_before_launch(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    campaign["run"].launched = True

    with pytest.raises(LAUNCHER.LaunchFailure, match="already active"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_inherited_gpu_visibility_override_is_not_in_launch_environment(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    campaign["runtime"].environment = lambda: {"CUDA_VISIBLE_DEVICES": "1,0"}

    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        assert "CUDA_VISIBLE_DEVICES" not in prepared.launch_environment
        assert "PATH" not in prepared.launch_environment
        assert "LD_PRELOAD" not in prepared.launch_environment
    finally:
        prepared.input_holder.close()
        prepared.directory_holder.close()
        prepared.executable_holder.close()


@pytest.mark.parametrize("name", [
    "LD_PRELOAD", "LD_LIBRARY_PATH", "PATH", "PYTHONPATH",
    "OMPI_MCA_pml", "PMIX_MCA_gds", "CUDA_VISIBLE_DEVICES",
    "KOKKOS_VISIBLE_DEVICES",
])
def test_dangerous_parent_environment_is_never_forwarded(
        tmp_path: Path, name: str) -> None:
    campaign = _campaign(tmp_path)
    campaign["runtime"].environment = lambda: {name: "/attacker/injection"}

    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        assert prepared.launch_environment == \
            campaign["plan"]["launch_contract"]["environment"]["values"]
        assert name not in prepared.launch_environment
    finally:
        prepared.input_holder.close()
        prepared.directory_holder.close()
        prepared.executable_holder.close()


def test_parent_cannot_override_bound_prrte_personality(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    campaign["runtime"].environment = lambda: {
        "PRTE_MCA_schizo_proxy": "prte",
    }

    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        assert prepared.launch_environment["PRTE_MCA_schizo_proxy"] == "ompi"
        assert prepared.launch_environment == \
            campaign["plan"]["launch_contract"]["environment"]["values"]
    finally:
        prepared.input_holder.close()
        prepared.directory_holder.close()
        prepared.executable_holder.close()


@pytest.mark.parametrize("mode", ["retained", "missing", "wrong_hash"])
def test_prrte_rank_environment_projection_is_exactly_plan_bound(
        tmp_path: Path, mode: str) -> None:
    campaign = _campaign(tmp_path)

    def mutate(plan: dict[str, Any]) -> None:
        projection = plan["launch_contract"]["environment"]["rank_projection"]
        if mode == "retained":
            projection["consumed_absent"] = []
        elif mode == "missing":
            del projection["inherited_values"]["LANG"]
        else:
            projection["sha256"] = "f" * 64

    _rewrite_plan(campaign, mutate)
    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="rank environment projection"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_runtime_cannot_override_plan_bound_nvidia_smi(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    campaign["runtime"].nvidia_smi = "/attacker/nvidia-smi"

    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="runtime nvidia-smi differs"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_proc_holder_probe_rejects_inaccessible_executable_before_popen(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    campaign = _campaign(tmp_path)
    real_access = LAUNCHER.os.access

    def deny_binary_proc_exec(path, mode, *args, **kwargs):
        text = os.fspath(path)
        if text.endswith(f"/fd/{LAUNCHER.BINARY_EXECUTABLE_FD}") and mode == os.X_OK:
            return False
        return real_access(path, mode, *args, **kwargs)

    monkeypatch.setattr(LAUNCHER.os, "access", deny_binary_proc_exec)
    popen_called = False

    def forbidden_popen(*args, **kwargs):
        nonlocal popen_called
        popen_called = True
        raise AssertionError("inaccessible exec holder reached Popen")

    campaign["runtime"].popen = forbidden_popen
    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="proc-holder access probe failed.*not executable"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])
    assert popen_called is False
    assert all(not LAUNCHER._fd_is_open(fd) for fd in (
        LAUNCHER.SOURCE_RESTART_FD, LAUNCHER.TRAJECTORY_FD,
        LAUNCHER.STATE_DIRECTORY_FD, LAUNCHER.EVIDENCE_DIRECTORY_FD,
        LAUNCHER.MPI_LAUNCHER_FD, LAUNCHER.BINARY_EXECUTABLE_FD,
    ))


def test_repository_audit_uses_plan_bound_git_and_removes_git_environment(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    observed_git_environments: list[dict[str, str]] = []
    original_run = campaign["run"]

    def recording_run(command, **kwargs):
        if command[0] == campaign["plan"]["tools"]["git"]["path"]:
            observed_git_environments.append(dict(kwargs["env"]))
        return original_run(command, **kwargs)

    campaign["runtime"].run = recording_run
    campaign["runtime"].environment = lambda: {
        "PATH": os.environ.get("PATH", ""),
        "GIT_DIR": "/attacker/repository",
        "GIT_WORK_TREE": "/attacker/worktree",
        "GIT_CONFIG_GLOBAL": "/attacker/config",
    }
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        assert len(observed_git_environments) == 3
        assert all(environment == {
            "LANG": "C", "LC_ALL": "C",
            "GIT_CONFIG_NOSYSTEM": "1", "GIT_CONFIG_GLOBAL": "/dev/null",
        } for environment in observed_git_environments)
        assert prepared.repository["git_tool"]["sha256"] == \
            campaign["plan"]["tools"]["git"]["sha256"]
        assert prepared.repository["git_environment_policy"] == \
            "explicit_clean_environment_v1"
    finally:
        prepared.input_holder.close()


class FakeProcess:
    def __init__(self, fake_run: FakeRun, return_code: int = 0) -> None:
        self.pid = 500
        self.fake_run = fake_run
        self.return_code = return_code
        self.exited = False
        self.terminated = False
        self.killed = False

    def poll(self):
        return self.return_code if self.exited else None

    def wait(self, timeout=None):
        self.exited = True
        self.fake_run.launched = False
        self.fake_run.completed = True
        return self.return_code

    def terminate(self):
        self.terminated = True
        self.exited = True
        self.fake_run.launched = False

    def kill(self):
        self.killed = True
        self.exited = True
        self.fake_run.launched = False


def _rank_environment(prepared: Any, rank: int, *, launcher_pid: int = 500,
                      hostname: str = "v100-node", port: int = 43210) \
        -> dict[str, str]:
    namespace_family = f"prterun-{hostname}-{launcher_pid}"
    namespace = f"{namespace_family}@1"
    uri = f"{namespace_family}@0.0;tcp4://127.0.0.1:{port}"
    state = str(prepared.state_dir)
    return {
        "HOME": prepared.launch_environment["HOME"],
        "LANG": "C", "LC_ALL": "C", "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
        "OMPI_ARGV": " ".join(prepared.athena_argv[1:]),
        "OMPI_COMMAND": prepared.athena_argv[0],
        "OMPI_COMM_WORLD_LOCAL_RANK": str(rank),
        "OMPI_COMM_WORLD_LOCAL_SIZE": str(prepared.world_size),
        "OMPI_COMM_WORLD_NODE_RANK": str(rank),
        "OMPI_COMM_WORLD_RANK": str(rank),
        "OMPI_COMM_WORLD_SIZE": str(prepared.world_size),
        "OMPI_FILE_LOCATION": f"/tmp/ompi.{launcher_pid}/1/{rank}",
        "OMPI_MCA_cpu_type": "x86_64", "OMPI_MCA_initial_wdir": state,
        "OMPI_MCA_num_procs": str(prepared.world_size),
        "OMPI_NUM_APP_CTX": "1", "OMPI_UNIVERSE_SIZE": "32",
        "OMPI_WORLD_LOCAL_SIZE": str(prepared.world_size),
        "OMPI_WORLD_SIZE": str(prepared.world_size),
        "PMIX_BFROP_BUFFER_TYPE": "PMIX_BFROP_BUFFER_NON_DESC",
        "PMIX_GDS_MODULE": "shmem2,hash", "PMIX_HOSTNAME": hostname,
        "PMIX_NAMESPACE": namespace, "PMIX_PARAM_FILE_PASSED": "1",
        "PMIX_RANK": str(rank), "PMIX_SECURITY_MODE": "native",
        "PMIX_SERVER_TMPDIR": f"/tmp/ompi.{launcher_pid}",
        "PMIX_SERVER_URI21": uri, "PMIX_SERVER_URI2": uri,
        "PMIX_SERVER_URI3": uri, "PMIX_SERVER_URI41": uri,
        "PMIX_SERVER_URI4": uri, "PMIX_SYSTEM_TMPDIR": "/tmp",
        "PMIX_VERSION": "5.0.9a1", "PRTE_LAUNCHED": "1",
        "PRTE_SHARED_FS": "FALSE", "OPAL_USER_PARAMS_GIVEN": "1",
        "PWD": state, "ZES_ENABLE_SYSMAN": "1",
    }


class FakeInspector:
    def __init__(self, prepared, fake_run: FakeRun,
                 *, bad_cmdline: bool = False,
                 gpu_visibility: str | None = None) -> None:
        self.prepared = prepared
        self.fake_run = fake_run
        self.bad_cmdline = bad_cmdline
        self.gpu_visibility = gpu_visibility

    def descendants(self, root_pid: int) -> set[int]:
        assert root_pid == 500
        return set(self.fake_run.rank_pids) | {550}

    def cmdline(self, pid: int) -> list[str]:
        if pid == 500:
            return list(self.prepared.launch_argv)
        if pid in self.fake_run.rank_pids:
            result = list(self.prepared.athena_argv)
            if self.bad_cmdline and pid == self.fake_run.rank_pids[-1]:
                result.append("output2/dt=1")
            return result
        return ["orted"]

    def environment(self, pid: int) -> dict[str, str]:
        index = self.fake_run.rank_pids.index(pid)
        result = _rank_environment(self.prepared, index)
        if self.gpu_visibility is not None:
            result["CUDA_VISIBLE_DEVICES"] = self.gpu_visibility
        return result

    def executable(self, pid: int) -> dict[str, Any]:
        if pid == 500:
            return dict(self.prepared.launcher)
        if pid in self.fake_run.rank_pids:
            return dict(self.prepared.binary)
        return {
            "device": 999, "inode": 999, "size": 1, "mtime_ns": 1,
            "ctime_ns": 1, "sha256": "f" * 64,
        }

    def start_time_ticks(self, pid: int) -> int | None:
        starts = {500: 1000, **{
            rank_pid: 2000 + rank
            for rank, rank_pid in enumerate(self.fake_run.rank_pids)
        }}
        if pid not in starts:
            return None
        if self.fake_run.launched or self.fake_run.retain_process_identity_after_wait:
            return starts[pid]
        if self.fake_run.reuse_process_ids_after_wait:
            return starts[pid] + 10000
        return None


def _proof_runtime(campaign: dict[str, Any], prepared,
                   *, bad_cmdline: bool = False,
                   gpu_visibility: str | None = None) -> LAUNCHER.Runtime:
    clock = [0.0]

    def fake_sleep(seconds: float) -> None:
        clock[0] += seconds

    def fake_killpg(pgid: int, number: int) -> None:
        assert pgid == 500
        live = (campaign["run"].launched or
                campaign["run"].retain_applications_after_wait or
                campaign["run"].retain_process_identity_after_wait)
        if number == 0:
            if live:
                return
            raise ProcessLookupError(pgid)
        campaign["run"].cleanup_signals.append(number)
        campaign["run"].launched = False
        campaign["run"].completed = True
        campaign["run"].retain_applications_after_wait = False
        campaign["run"].retain_process_identity_after_wait = False

    return LAUNCHER.Runtime(
        run=campaign["run"], sleep=fake_sleep,
        monotonic=lambda: clock[0],
        hostname=lambda: "v100-node",
        inspector=FakeInspector(prepared, campaign["run"],
                                bad_cmdline=bad_cmdline,
                                gpu_visibility=gpu_visibility),
        nvidia_smi=prepared.execution_tools["nvidia_smi"]["path"],
        getpgid=lambda pid: pid,
        killpg=fake_killpg,
        statvfs=_ample_statvfs, fstatvfs=_ample_statvfs,
    )


def _lifecycle_fixture(tmp_path: Path):
    campaign = _campaign(tmp_path)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    runtime = _proof_runtime(campaign, prepared)
    process = FakeProcess(campaign["run"])

    def fake_popen(argv, **kwargs):
        campaign["run"].launched = True
        return process

    runtime.popen = fake_popen
    evidence = campaign["evidence"]
    paths = LAUNCHER.OutputPaths(
        evidence / "segment.launch.ready", evidence / "segment.completion.ready",
        evidence / "run.log", evidence / "exit.status",
        evidence / "gpu-before.csv", evidence / "gpu-after.csv",
    )
    return campaign, prepared, runtime, process, paths


def test_bound_second_disk_failure_cleans_only_staged_inputs_without_popen(
        tmp_path: Path) -> None:
    campaign, prepared, runtime, _, paths = _lifecycle_fixture(tmp_path)
    popen_called = False

    def forbidden_popen(*args, **kwargs):
        nonlocal popen_called
        popen_called = True
        raise AssertionError("second disk rejection called Popen")

    runtime.fstatvfs = _exactly_75_percent_used_statvfs
    runtime.popen = forbidden_popen
    with pytest.raises(LAUNCHER.LaunchFailure, match="before_spawn.*75%"):
        LAUNCHER.run_segment(
            prepared, campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)

    assert popen_called is False
    assert not any(campaign["staging"].iterdir())
    assert campaign["staging"].stat().st_mode & 0o777 == 0o700
    assert not paths.launch_record.exists()
    assert not paths.gpu_before.exists()


def test_live_proof_binds_every_rank_pid_argv_exe_and_gpu(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    runtime = _proof_runtime(campaign, prepared)
    campaign["run"].launched = True
    process = FakeProcess(campaign["run"])

    proof = LAUNCHER.prove_running_launch(
        process, prepared, runtime, timeout_seconds=1,
        plan_path=campaign["plan_path"])

    assert proof["mpirun_pid"] == 500
    assert [row["global_rank"] for row in proof["ranks"]] == [0, 1]
    assert [row["gpu_uuid"] for row in proof["ranks"]] == ["GPU-0", "GPU-1"]
    assert all(row["cmdline"] == list(prepared.athena_argv)
               for row in proof["ranks"])


def test_holder_survives_original_input_path_replacement(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    runtime = _proof_runtime(campaign, prepared)
    source_original = campaign["source"].with_suffix(".bound")
    trajectory_original = campaign["trajectory"].with_suffix(".bound")
    campaign["source"].rename(source_original)
    campaign["source"].write_bytes(b"untrusted replacement\n")
    campaign["trajectory"].rename(trajectory_original)
    campaign["trajectory"].write_bytes(b"untrusted replacement\n")
    campaign["run"].launched = True

    proof = LAUNCHER.prove_running_launch(
        FakeProcess(campaign["run"]), prepared, runtime,
        timeout_seconds=1, plan_path=campaign["plan_path"])

    assert proof["input_transport"]["roles"]["source_restart"]["sha256"] == \
        prepared.source_restart["sha256"]
    assert proof["input_transport"]["roles"]["trajectory"]["sha256"] == \
        prepared.trajectory["sha256"]


def test_run_rejects_early_holder_fd_close_and_closes_remaining_fd(
        tmp_path: Path) -> None:
    campaign, prepared, runtime, _, paths = _lifecycle_fixture(tmp_path)
    os.close(LAUNCHER.SOURCE_RESTART_FD)

    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="cannot audit holder source_restart descriptor"):
        LAUNCHER.run_segment(
            prepared, campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)

    assert not LAUNCHER._fd_is_open(LAUNCHER.SOURCE_RESTART_FD)
    assert not LAUNCHER._fd_is_open(LAUNCHER.TRAJECTORY_FD)
    assert not paths.launch_record.exists()


def test_live_proof_rejects_one_rank_with_noncanonical_argv(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    runtime = _proof_runtime(campaign, prepared, bad_cmdline=True)
    campaign["run"].launched = True
    ticks = iter((0.0, 2.0))
    runtime.monotonic = lambda: next(ticks)

    with pytest.raises(LAUNCHER.LaunchFailure, match="timed out.*argv differs"):
        LAUNCHER.prove_running_launch(
            FakeProcess(campaign["run"]), prepared, runtime,
            timeout_seconds=1, poll_seconds=0)


def test_live_proof_rejects_an_unrelated_gpu_process(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    runtime = _proof_runtime(campaign, prepared)
    campaign["run"].launched = True
    campaign["run"].extra_application = True
    ticks = iter((0.0, 2.0))
    runtime.monotonic = lambda: next(ticks)

    with pytest.raises(LAUNCHER.LaunchFailure, match="timed out.*exactly"):
        LAUNCHER.prove_running_launch(
            FakeProcess(campaign["run"]), prepared, runtime,
            timeout_seconds=1, poll_seconds=0)


def test_missing_openmpi_rank_environment_is_rejected() -> None:
    with pytest.raises(LAUNCHER.LaunchFailure, match="key closure differs"):
        LAUNCHER._parse_rank_environment({}, 2, 600, {
            "HOME": "/root", "LANG": "C", "LC_ALL": "C",
            "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
            "PRTE_MCA_schizo_proxy": "ompi",
        }, launcher_pid=500, hostname="v100-node", state_dir=Path("/state"),
            athena_argv=("/proc/1/fd/205", "arg"))


def test_rank_projection_requires_prrte_personality_to_be_consumed() -> None:
    environment = {
        "OMPI_COMM_WORLD_RANK": "0",
        "OMPI_COMM_WORLD_SIZE": "1",
        "OMPI_COMM_WORLD_LOCAL_RANK": "0",
        "OMPI_COMM_WORLD_LOCAL_SIZE": "1",
        "HOME": "/root", "LANG": "C", "LC_ALL": "C",
        "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
        "PRTE_MCA_schizo_proxy": "ompi",
    }
    with pytest.raises(LAUNCHER.LaunchFailure, match="key closure differs"):
        LAUNCHER._parse_rank_environment(environment, 1, 600, {
            "HOME": "/root", "LANG": "C", "LC_ALL": "C",
            "CUDA_DEVICE_ORDER": "PCI_BUS_ID",
            "PRTE_MCA_schizo_proxy": "ompi",
        }, launcher_pid=500, hostname="v100-node", state_dir=Path("/state"),
            athena_argv=("/proc/1/fd/205", "arg"))


@pytest.mark.parametrize("name", [
    "OMPI_MCA_pml", "PMIX_MCA_gds", "PRTE_MCA_unreviewed",
])
def test_rank_environment_rejects_every_unknown_runtime_prefix_key(
        tmp_path: Path, name: str) -> None:
    campaign = _campaign(tmp_path, world_size=1)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        environment = _rank_environment(prepared, 0)
        environment[name] = "attacker"
        with pytest.raises(LAUNCHER.LaunchFailure, match="key closure differs"):
            LAUNCHER._parse_rank_environment(
                environment, 1, 600, prepared.launch_environment,
                launcher_pid=500, hostname="v100-node",
                state_dir=prepared.state_dir,
                athena_argv=prepared.athena_argv)
    finally:
        prepared.close()


@pytest.mark.parametrize(("key", "value", "message"), [
    ("OMPI_ARGV", "attacker", "OMPI_ARGV"),
    ("PMIX_NAMESPACE", "attacker", "PMIX_NAMESPACE"),
    ("PMIX_SERVER_URI4", "bad", "URI aliases"),
])
def test_rank_environment_rejects_tampered_derived_values(
        tmp_path: Path, key: str, value: str, message: str) -> None:
    campaign = _campaign(tmp_path, world_size=1)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        environment = _rank_environment(prepared, 0)
        environment[key] = value
        with pytest.raises(LAUNCHER.LaunchFailure, match=message):
            LAUNCHER._parse_rank_environment(
                environment, 1, 600, prepared.launch_environment,
                launcher_pid=500, hostname="v100-node",
                state_dir=prepared.state_dir,
                athena_argv=prepared.athena_argv)
    finally:
        prepared.close()


def test_rank_environment_accepts_namespace_family_server_rank_uri(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path, world_size=1)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        environment = _rank_environment(prepared, 0)
        assert environment["PMIX_NAMESPACE"] == "prterun-v100-node-500@1"
        assert environment["PMIX_SERVER_URI21"] == (
            "prterun-v100-node-500@0.0;tcp4://127.0.0.1:43210")
        _, _, selected = LAUNCHER._parse_rank_environment(
            environment, 1, 600, prepared.launch_environment,
            launcher_pid=500, hostname="v100-node",
            state_dir=prepared.state_dir,
            athena_argv=prepared.athena_argv)
        assert selected == environment
    finally:
        prepared.close()


def test_rank_environment_rejects_namespace_identifier_as_server_uri_family(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path, world_size=1)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        environment = _rank_environment(prepared, 0)
        invalid_uri = (
            f'{environment["PMIX_NAMESPACE"]}@0.0;tcp4://127.0.0.1:43210')
        for key in LAUNCHER.RANK_URI_ENVIRONMENT_KEYS:
            environment[key] = invalid_uri
        with pytest.raises(LAUNCHER.LaunchFailure,
                           match="not the derived loopback URI"):
            LAUNCHER._parse_rank_environment(
                environment, 1, 600, prepared.launch_environment,
                launcher_pid=500, hostname="v100-node",
                state_dir=prepared.state_dir,
                athena_argv=prepared.athena_argv)
    finally:
        prepared.close()


def test_prepare_rejects_mca_file_created_after_plan(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    config = campaign["mca_prefix"] / "etc" / "openmpi-mca-params.conf"
    config.parent.mkdir()
    config.write_text("pml = attacker\n", encoding="ascii")
    config.chmod(0o644)

    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="MCA configuration.*differs"):
        LAUNCHER.prepare_launch(
            campaign["plan_path"], campaign["state"], campaign["runtime"])


def test_live_proof_rejects_rank_specific_cuda_visibility(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    runtime = _proof_runtime(campaign, prepared, gpu_visibility="0")
    campaign["run"].launched = True
    ticks = iter((0.0, 2.0))
    runtime.monotonic = lambda: next(ticks)

    with pytest.raises(LAUNCHER.LaunchFailure, match="timed out.*key closure"):
        LAUNCHER.prove_running_launch(
            FakeProcess(campaign["run"]), prepared, runtime,
            timeout_seconds=1, poll_seconds=0)


def test_live_proof_rejects_rank_to_gpu_permutation(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    runtime = _proof_runtime(campaign, prepared)
    campaign["run"].launched = True
    campaign["run"].permute_applications = True
    ticks = iter((0.0, 2.0))
    runtime.monotonic = lambda: next(ticks)

    with pytest.raises(LAUNCHER.LaunchFailure, match="timed out.*PCI-ordered"):
        LAUNCHER.prove_running_launch(
            FakeProcess(campaign["run"]), prepared, runtime,
            timeout_seconds=1, poll_seconds=0)


def test_full_lifecycle_uses_no_shell_and_publishes_closed_evidence(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    campaign = _campaign(tmp_path)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    runtime = _proof_runtime(campaign, prepared)
    captured: dict[str, Any] = {}
    process = FakeProcess(campaign["run"])

    def fake_popen(argv, **kwargs):
        captured["argv"] = argv
        captured["kwargs"] = kwargs
        campaign["run"].launched = True
        return process

    runtime.popen = fake_popen
    evidence = campaign["evidence"]
    paths = LAUNCHER.OutputPaths(
        evidence / "segment.launch.ready", evidence / "segment.completion.ready",
        evidence / "run.log",
        evidence / "exit.status", evidence / "gpu-before.csv",
        evidence / "gpu-after.csv",
    )
    installs: list[tuple[str, str]] = []
    original_install = LAUNCHER._install_immutable_bytes

    def counted_install(path: Path, content: bytes):
        installs.append(("bytes", path.name))
        return original_install(path, content)

    monkeypatch.setattr(LAUNCHER, "_install_immutable_bytes", counted_install)

    status = LAUNCHER.run_segment(
        prepared, campaign["plan_path"], paths, runtime,
        proof_timeout_seconds=1)

    assert status == 0
    assert captured["argv"] == list(prepared.launch_argv)
    assert captured["kwargs"]["shell"] is False
    assert captured["kwargs"]["executable"] == \
        f"/proc/{prepared.input_holder.pid}/fd/204"
    assert captured["kwargs"]["cwd"] == \
        f"/proc/{prepared.input_holder.pid}/fd/202"
    assert captured["kwargs"]["env"] == prepared.launch_environment
    assert "CUDA_VISIBLE_DEVICES" not in captured["kwargs"]["env"]
    assert "KOKKOS_VISIBLE_DEVICES" not in captured["kwargs"]["env"]
    assert captured["kwargs"]["env"]["CUDA_DEVICE_ORDER"] == "PCI_BUS_ID"
    record = json.loads(paths.launch_record.read_text(encoding="utf-8"))
    assert record["status"] == "ready"
    assert record["launch_argv"] == list(prepared.launch_argv)
    assert record["launch_environment"] == \
        campaign["plan"]["launch_contract"]["environment"]
    assert record["nvidia_smi"]["sha256"] == \
        campaign["plan"]["tools"]["nvidia_smi"]["sha256"]
    assert record["executable_transport"]["roles"]["launcher"]["fd"] == 204
    assert record["executable_transport"]["roles"]["binary"]["fd"] == 205
    assert len(record["ranks"]) == prepared.world_size
    assert paths.exit_status.read_text(encoding="ascii") == "0\n"
    assert paths.gpu_before.read_text(encoding="ascii").splitlines()[0] == (
        "0,GPU-0,00000000:10:00.0,0,0,0,32768,5")
    assert paths.gpu_before.read_text(encoding="ascii").splitlines()[1] == (
        "1,GPU-1,00000000:11:00.0,1,0,0,32768,6")
    completion = json.loads(paths.completion_record.read_text(encoding="utf-8"))
    assert completion["kind"] == "athenak_segment_completion"
    assert completion["status"] == "ready"
    assert completion["return_code"] == 0
    assert completion["quiescence"]["all_original_identities_gone"] is True
    assert completion["input_transport"]["kind"] == \
        "linux_proc_holder_fd_v1"
    assert completion["input_transport"]["closure"][
        "all_holder_fds_closed"] is True
    assert set(completion["input_transport"]["closure"]["roles"]) == {
        "source_restart", "trajectory",
    }
    assert not LAUNCHER._fd_is_open(LAUNCHER.SOURCE_RESTART_FD)
    assert not LAUNCHER._fd_is_open(LAUNCHER.TRAJECTORY_FD)
    assert set(completion["artifacts"]) == {
        "plan", "launch_record", "run_log", "exit_status",
        "gpu_before", "gpu_after",
    }
    assert installs.count(("bytes", paths.gpu_before.name)) == 1
    assert installs.count(("bytes", paths.gpu_after.name)) == 1
    assert installs.count(("bytes", paths.exit_status.name)) == 1
    assert installs.count(("bytes", paths.launch_record.name)) == 1
    assert installs.count(("bytes", paths.completion_record.name)) == 1
    assert installs[-1] == ("bytes", paths.completion_record.name)
    for path in paths.__dict__.values():
        assert path.stat().st_mode & 0o222 == 0


def test_completion_is_withheld_when_gpu_context_remains(tmp_path: Path) -> None:
    campaign, prepared, runtime, process, paths = _lifecycle_fixture(tmp_path)
    campaign["run"].retain_applications_after_wait = True

    with pytest.raises(LAUNCHER.LaunchFailure, match="GPU compute contexts remain"):
        LAUNCHER.run_segment(
            prepared, campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)

    assert process.exited is True
    assert campaign["run"].cleanup_signals == [signal.SIGTERM]
    assert not paths.completion_record.exists()
    assert not LAUNCHER._fd_is_open(LAUNCHER.SOURCE_RESTART_FD)
    assert not LAUNCHER._fd_is_open(LAUNCHER.TRAJECTORY_FD)


def test_completion_is_withheld_when_recorded_pid_identity_remains(
        tmp_path: Path) -> None:
    campaign, prepared, runtime, process, paths = _lifecycle_fixture(tmp_path)
    campaign["run"].retain_process_identity_after_wait = True

    with pytest.raises(LAUNCHER.LaunchFailure, match="identity is still live"):
        LAUNCHER.run_segment(
            prepared, campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)

    assert process.exited is True
    assert campaign["run"].cleanup_signals == [signal.SIGTERM]
    assert not paths.completion_record.exists()


def test_wait_interruption_terminates_whole_managed_group_before_holder_close(
        tmp_path: Path) -> None:
    campaign, prepared, runtime, process, paths = _lifecycle_fixture(tmp_path)
    original_wait = process.wait
    first = True

    def interrupted_wait(timeout=None):
        nonlocal first
        if first and timeout is None:
            first = False
            raise KeyboardInterrupt("injected SIGINT")
        return original_wait(timeout=timeout)

    process.wait = interrupted_wait
    with pytest.raises(KeyboardInterrupt, match="injected SIGINT"):
        LAUNCHER.run_segment(
            prepared, campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)

    assert campaign["run"].cleanup_signals == [signal.SIGTERM]
    assert paths.launch_record.exists()
    assert not paths.completion_record.exists()
    assert not LAUNCHER._fd_is_open(LAUNCHER.SOURCE_RESTART_FD)
    assert not LAUNCHER._fd_is_open(LAUNCHER.TRAJECTORY_FD)
    assert not LAUNCHER._fd_is_open(LAUNCHER.MPI_LAUNCHER_FD)
    assert not LAUNCHER._fd_is_open(LAUNCHER.BINARY_EXECUTABLE_FD)


def test_sigterm_uses_managed_group_cleanup_and_restores_signal_policy(
        tmp_path: Path) -> None:
    campaign, prepared, runtime, process, paths = _lifecycle_fixture(tmp_path)
    sentinels = {
        signum: object() for signum in LAUNCHER.MANAGED_TERMINATION_SIGNALS
    }
    handlers: dict[int, Any] = dict(sentinels)
    runtime.get_signal_handler = lambda signum: handlers[signum]

    def set_handler(signum: int, handler: Any) -> Any:
        old = handlers[signum]
        handlers[signum] = handler
        return old

    runtime.set_signal_handler = set_handler
    original_wait = process.wait
    injected = False

    def sigterm_wait(timeout=None):
        nonlocal injected
        if not injected and timeout is None:
            injected = True
            handlers[signal.SIGTERM](signal.SIGTERM, None)
        return original_wait(timeout=timeout)

    process.wait = sigterm_wait
    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="managed termination signal SIGTERM"):
        LAUNCHER.run_segment(
            prepared, campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)

    assert campaign["run"].cleanup_signals == [signal.SIGTERM]
    assert handlers == sentinels
    assert paths.launch_record.exists()
    assert not paths.completion_record.exists()
    assert all(not LAUNCHER._fd_is_open(fd) for fd in (
        LAUNCHER.SOURCE_RESTART_FD, LAUNCHER.TRAJECTORY_FD,
        LAUNCHER.STATE_DIRECTORY_FD, LAUNCHER.EVIDENCE_DIRECTORY_FD,
        LAUNCHER.MPI_LAUNCHER_FD, LAUNCHER.BINARY_EXECUTABLE_FD,
    ))


def test_completion_accepts_pid_reuse_but_records_it(tmp_path: Path) -> None:
    campaign, prepared, runtime, _, paths = _lifecycle_fixture(tmp_path)
    campaign["run"].reuse_process_ids_after_wait = True

    assert LAUNCHER.run_segment(
        prepared, campaign["plan_path"], paths, runtime,
        proof_timeout_seconds=1) == 0

    completion = json.loads(paths.completion_record.read_text(encoding="utf-8"))
    observations = completion["quiescence"]["process_identities"]
    assert observations
    assert {row["state"] for row in observations} == {"pid_reused"}
    assert all(row["original_identity_gone"] is True for row in observations)


def test_completion_boundedly_polls_until_managed_group_disappears(
        tmp_path: Path) -> None:
    campaign, prepared, runtime, _, paths = _lifecycle_fixture(tmp_path)
    original_killpg = runtime.killpg
    post_exit_polls = [0]

    def transient_group(pgid: int, number: int) -> None:
        if number == 0 and campaign["run"].completed:
            post_exit_polls[0] += 1
            if post_exit_polls[0] < 3:
                return
        original_killpg(pgid, number)

    runtime.killpg = transient_group
    assert LAUNCHER.run_segment(
        prepared, campaign["plan_path"], paths, runtime,
        proof_timeout_seconds=1) == 0

    completion = json.loads(paths.completion_record.read_text(encoding="utf-8"))
    assert post_exit_polls[0] == 3
    assert completion["quiescence"]["managed_process_group_gone"] is True


def test_completion_is_withheld_when_managed_group_never_disappears(
        tmp_path: Path) -> None:
    campaign, prepared, runtime, process, paths = _lifecycle_fixture(tmp_path)
    original_killpg = runtime.killpg

    def persistent_group(pgid: int, number: int) -> None:
        if number == 0 and campaign["run"].completed:
            return
        original_killpg(pgid, number)

    runtime.killpg = persistent_group
    with pytest.raises(
            LAUNCHER.LaunchFailure,
            match="managed MPI process group remains after bounded"):
        LAUNCHER.run_segment(
            prepared, campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)

    assert process.exited is True
    assert campaign["run"].cleanup_signals == [signal.SIGTERM, signal.SIGKILL]
    assert not paths.completion_record.exists()


def test_completion_is_withheld_when_gpu_after_install_fails(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    campaign, prepared, runtime, _, paths = _lifecycle_fixture(tmp_path)
    original = LAUNCHER._install_immutable_bytes

    def fail_gpu_after(path: Path, content: bytes):
        if path.name == paths.gpu_after.name:
            raise LAUNCHER.LaunchFailure("injected GPU-after failure")
        return original(path, content)

    monkeypatch.setattr(LAUNCHER, "_install_immutable_bytes", fail_gpu_after)
    with pytest.raises(LAUNCHER.LaunchFailure, match="injected GPU-after"):
        LAUNCHER.run_segment(
            prepared, campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)

    assert not paths.completion_record.exists()


def test_completion_rejects_tampered_launch_record_before_last_marker(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    campaign, prepared, runtime, _, paths = _lifecycle_fixture(tmp_path)
    original = LAUNCHER._install_immutable_bytes

    def tamper_after_exit(path: Path, content: bytes):
        binding = original(path, content)
        if path.name == paths.exit_status.name:
            paths.launch_record.chmod(0o644)
            paths.launch_record.write_text("{}\n", encoding="utf-8")
            paths.launch_record.chmod(0o444)
        return binding

    monkeypatch.setattr(LAUNCHER, "_install_immutable_bytes", tamper_after_exit)
    with pytest.raises(LAUNCHER.LaunchFailure, match="launch_record changed"):
        LAUNCHER.run_segment(
            prepared, campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)

    assert not paths.completion_record.exists()


def test_launch_rejects_staged_source_tamper_before_spawn(tmp_path: Path) -> None:
    campaign, prepared, runtime, _, paths = _lifecycle_fixture(tmp_path)
    campaign["staging"].chmod(0o755)
    campaign["staged_source"].chmod(0o644)
    campaign["staged_source"].write_bytes(b"tampered restart\n")
    campaign["staged_source"].chmod(0o444)
    campaign["staging"].chmod(0o555)

    with pytest.raises(
            LAUNCHER.LaunchFailure,
            match="holder source_restart descriptor identity changed"):
        LAUNCHER.run_segment(
            prepared, campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)

    assert not paths.launch_record.exists()
    assert not paths.completion_record.exists()


def test_executable_path_replacement_cannot_change_spawned_inodes(
        tmp_path: Path) -> None:
    campaign, prepared, runtime, _, paths = _lifecycle_fixture(tmp_path)
    captured: dict[str, Any] = {}
    process = FakeProcess(campaign["run"])

    def fake_popen(argv, **kwargs):
        captured["argv"] = list(argv)
        captured["executable"] = kwargs["executable"]
        campaign["run"].launched = True
        return process

    runtime.popen = fake_popen
    for path in (campaign["binary"], campaign["mpirun"]):
        bound = path.with_suffix(path.suffix + ".planned-inode")
        path.rename(bound)
        path.write_bytes(b"attacker replacement executable\n")
        path.chmod(0o755)

    assert LAUNCHER.run_segment(
        prepared, campaign["plan_path"], paths, runtime,
        proof_timeout_seconds=1) == 0
    assert captured["executable"] == \
        f"/proc/{prepared.input_holder.pid}/fd/204"
    assert captured["argv"][len(prepared.mpi_argv)] == \
        f"/proc/{prepared.input_holder.pid}/fd/205"


@pytest.mark.parametrize("role", ["state_dir", "evidence_dir"])
def test_directory_path_replacement_is_rejected_before_spawn(
        tmp_path: Path, role: str) -> None:
    campaign, prepared, runtime, _, paths = _lifecycle_fixture(tmp_path)
    popen_called = False

    def forbidden_popen(*args, **kwargs):
        nonlocal popen_called
        popen_called = True
        raise AssertionError("replaced directory reached Popen")

    runtime.popen = forbidden_popen
    path = campaign["state"] if role == "state_dir" else campaign["evidence"]
    bound = path.with_name(path.name + ".planned-inode")
    path.rename(bound)
    path.mkdir(mode=0o700)
    if role == "evidence_dir":
        # Preserve the immutable plan pathname so plan-identity checking reaches
        # the directory-holder identity proof rather than failing on ENOENT.
        os.link(bound / "segment.plan.json", path / "segment.plan.json")

    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="run plan path/identity|directory path was replaced"):
        LAUNCHER.run_segment(
            prepared, path / "segment.plan.json" if role == "evidence_dir"
            else campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)
    assert popen_called is False


def test_eight_rank_proc_fd_argv_and_gpu_mapping_are_exact(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path, world_size=8)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    try:
        assert prepared.world_size == 8
        assert prepared.mpi_argv[-1] == "8"
        assert prepared.athena_argv[0] == \
            f"/proc/{prepared.input_holder.pid}/fd/205"
        assert prepared.athena_argv[5] == \
            f"/proc/{prepared.input_holder.pid}/fd/202"
        runtime = _proof_runtime(campaign, prepared)
        campaign["run"].launched = True
        proof = LAUNCHER.prove_running_launch(
            FakeProcess(campaign["run"]), prepared, runtime,
            timeout_seconds=1, plan_path=campaign["plan_path"])
        assert [row["global_rank"] for row in proof["ranks"]] == list(range(8))
        assert [row["gpu_cuda_ordinal"] for row in proof["ranks"]] == list(range(8))
    finally:
        prepared.input_holder.close()
        prepared.directory_holder.close()
        prepared.executable_holder.close()


def test_nonprepare_cli_requires_explicit_completion_path(
        tmp_path: Path, capsys) -> None:
    campaign = _campaign(tmp_path)
    evidence = campaign["evidence"]

    with pytest.raises(SystemExit) as error:
        LAUNCHER.main([
            "--plan", str(campaign["plan_path"]),
            "--state-dir", str(campaign["state"]),
            "--launch-record", str(evidence / "segment.launch.ready"),
            "--run-log", str(evidence / "run.log"),
            "--exit-status", str(evidence / "exit.status"),
            "--gpu-before", str(evidence / "gpu-before.csv"),
            "--gpu-after", str(evidence / "gpu-after.csv"),
        ], campaign["runtime"])

    assert error.value.code == 2
    assert "--completion-record" in capsys.readouterr().err


def test_launch_evidence_must_be_outside_state_directory(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    evidence = campaign["state"] / "evidence"
    evidence.mkdir(exist_ok=True)
    paths = LAUNCHER.OutputPaths(
        evidence / "segment.launch.ready", evidence / "segment.completion.ready",
        evidence / "run.log",
        evidence / "exit.status", evidence / "gpu-before.csv",
        evidence / "gpu-after.csv",
    )

    with pytest.raises(LAUNCHER.LaunchFailure, match="CLI evidence paths differ"):
        LAUNCHER.run_segment(prepared, campaign["plan_path"], paths,
                             campaign["runtime"], proof_timeout_seconds=1)


def test_launch_evidence_must_share_one_dedicated_directory(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    first, second = tmp_path / "evidence-a", tmp_path / "evidence-b"
    first.mkdir()
    second.mkdir()
    paths = LAUNCHER.OutputPaths(
        first / "segment.launch.ready", first / "segment.completion.ready",
        first / "run.log",
        first / "exit.status", first / "gpu-before.csv",
        second / "gpu-after.csv",
    )

    with pytest.raises(LAUNCHER.LaunchFailure, match="CLI evidence paths differ"):
        LAUNCHER.run_segment(prepared, campaign["plan_path"], paths,
                             campaign["runtime"], proof_timeout_seconds=1)


def test_state_directory_must_not_be_nested_in_evidence_directory(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    paths = LAUNCHER.OutputPaths(
        tmp_path / "segment.launch.ready", tmp_path / "segment.completion.ready",
        tmp_path / "run.log",
        tmp_path / "exit.status", tmp_path / "gpu-before.csv",
        tmp_path / "gpu-after.csv",
    )

    with pytest.raises(LAUNCHER.LaunchFailure, match="CLI evidence paths differ"):
        LAUNCHER.run_segment(prepared, campaign["plan_path"], paths,
                             campaign["runtime"], proof_timeout_seconds=1)


def test_launch_proof_failure_cleans_up_only_unqualified_job(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    prepared = LAUNCHER.prepare_launch(
        campaign["plan_path"], campaign["state"], campaign["runtime"])
    runtime = _proof_runtime(campaign, prepared, bad_cmdline=True)
    ticks = iter((0.0, 2.0))
    runtime.monotonic = lambda: next(ticks)
    process = FakeProcess(campaign["run"])
    def fake_popen(argv, **kwargs):
        campaign["run"].launched = True
        return process

    runtime.popen = fake_popen
    evidence = campaign["evidence"]
    paths = LAUNCHER.OutputPaths(
        evidence / "segment.launch.ready", evidence / "segment.completion.ready",
        evidence / "run.log",
        evidence / "exit.status", evidence / "gpu-before.csv",
        evidence / "gpu-after.csv",
    )

    with pytest.raises(LAUNCHER.LaunchFailure, match="timed out"):
        LAUNCHER.run_segment(
            prepared, campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)

    assert campaign["run"].cleanup_signals == [signal.SIGTERM]
    assert not paths.launch_record.exists()
    assert not paths.exit_status.exists()
    assert not LAUNCHER._fd_is_open(LAUNCHER.SOURCE_RESTART_FD)
    assert not LAUNCHER._fd_is_open(LAUNCHER.TRAJECTORY_FD)


def test_getpgid_failure_cleans_whole_group_before_any_holder_closes(
        tmp_path: Path) -> None:
    campaign, prepared, runtime, process, paths = _lifecycle_fixture(tmp_path)
    runtime.getpgid = lambda pid: (_ for _ in ()).throw(
        OSError("injected getpgid failure"))
    original_killpg = runtime.killpg
    open_during_cleanup: list[int] = []
    all_holder_fds = (
        LAUNCHER.SOURCE_RESTART_FD, LAUNCHER.TRAJECTORY_FD,
        LAUNCHER.STATE_DIRECTORY_FD, LAUNCHER.EVIDENCE_DIRECTORY_FD,
        LAUNCHER.MPI_LAUNCHER_FD, LAUNCHER.BINARY_EXECUTABLE_FD,
    )

    def asserting_killpg(pgid: int, number: int) -> None:
        if number in (signal.SIGTERM, signal.SIGKILL):
            assert all(LAUNCHER._fd_is_open(fd) for fd in all_holder_fds)
            open_during_cleanup.extend(all_holder_fds)
        original_killpg(pgid, number)

    runtime.killpg = asserting_killpg
    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="cannot bind new MPI process group"):
        LAUNCHER.run_segment(
            prepared, campaign["plan_path"], paths, runtime,
            proof_timeout_seconds=1)

    assert process.terminated is False
    assert campaign["run"].cleanup_signals == [signal.SIGTERM]
    assert set(open_during_cleanup) == set(all_holder_fds)
    assert all(not LAUNCHER._fd_is_open(fd) for fd in all_holder_fds)
    assert not paths.launch_record.exists()
    assert not paths.completion_record.exists()


def test_real_planner_artifact_is_accepted_by_launcher_prepare(
        tmp_path: Path) -> None:
    """Exercise the planner/launcher schema boundary without starting MPI."""

    import test_suite.test_bbh_segment_plan as plan_tests

    campaign = plan_tests._campaign(tmp_path)
    result = plan_tests._run(campaign)
    assert result.returncode == 0, result.stderr
    plan = json.loads(campaign["output"].read_text(encoding="utf-8"))
    ranks = plan["policy"]["ranks"]
    fake_run = FakeRun(ranks)
    runtime = LAUNCHER.Runtime(
        run=fake_run, statvfs=_ample_statvfs, fstatvfs=_ample_statvfs)
    prepared = LAUNCHER.prepare_launch(
        campaign["output"], campaign["state_dir"], runtime)
    try:
        assert prepared.mpi_argv == tuple(plan["launch_contract"]["mpi_argv"])
        assert prepared.athena_argv == tuple(
            token.replace("{holder_pid}", str(prepared.input_holder.pid))
            for token in plan["launch_contract"]["athena_argv_template"]
        )
        assert prepared.athena_argv[-2].startswith(
            "problem/trajectory_file=/proc/")
        assert prepared.athena_argv[-1] == "output3/dt=4.8"
    finally:
        prepared.input_holder.close()


def test_real_planner_global_cadence_transition_crosses_checker_and_launcher(
        tmp_path: Path) -> None:
    """Exercise the real 48M -> 10M paired contract across all three tools."""

    import check_athenak_segment as checker
    import test_suite.test_bbh_segment_plan as plan_tests

    campaign = plan_tests._campaign(tmp_path)
    result = plan_tests._run(
        campaign, root_steps="20", target_global_dt="10")
    assert result.returncode == 0, result.stderr
    plan = json.loads(campaign["output"].read_text(encoding="utf-8"))
    assert checker.validate_plan(plan)["final_cycle"] == 51

    fake_run = FakeRun(plan["policy"]["ranks"])
    runtime = LAUNCHER.Runtime(
        run=fake_run, statvfs=_ample_statvfs, fstatvfs=_ample_statvfs)
    prepared = LAUNCHER.prepare_launch(
        campaign["output"], campaign["state_dir"], runtime)
    try:
        assert prepared.athena_argv[-2:] == (
            "output2/dt=10.0", "output5/dt=10.0")
        transition = prepared.plan["global_cadence_transition"]
        assert transition["kind"] == "tighten_pair_v1"
        assert transition["streams"][0]["target_schedule"] == \
            transition["streams"][1]["target_schedule"]
    finally:
        prepared.close()


def test_real_planner_mocked_eight_rank_lifecycle_passes_checker_consumers(
        tmp_path: Path) -> None:
    """Exercise planner -> launcher -> launch/completion checker as one chain."""

    import check_athenak_segment as checker
    import test_suite.test_bbh_segment_plan as plan_tests

    campaign = plan_tests._campaign(tmp_path)
    result = plan_tests._run(campaign, ranks="8")
    assert result.returncode == 0, result.stderr
    plan_path = campaign["output"]
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    assert checker.validate_plan(plan)["ranks"] == 8

    fake_run = FakeRun(8)
    prepare_runtime = LAUNCHER.Runtime(
        run=fake_run, statvfs=_ample_statvfs, fstatvfs=_ample_statvfs)
    prepared = LAUNCHER.prepare_launch(
        plan_path, campaign["state_dir"], prepare_runtime)
    runtime = _proof_runtime({"run": fake_run}, prepared)
    process = FakeProcess(fake_run)

    def fake_popen(argv, **kwargs):
        fake_run.launched = True
        return process

    runtime.popen = fake_popen
    evidence = {
        name: Path(value) for name, value in
        plan["launch_contract"]["evidence"].items()
    }
    outputs = LAUNCHER.OutputPaths(*(
        evidence[name] for name in LAUNCHER.EVIDENCE_NAMES
    ))
    assert LAUNCHER.run_segment(
        prepared, plan_path, outputs, runtime, proof_timeout_seconds=1) == 0

    gpu_before_binding = checker.stable_sha256(evidence["gpu_before"])
    launch = checker.audit_launch_record(
        evidence["launch_record"], plan, plan_path, campaign["state_dir"],
        gpu_before_binding)
    assert launch["world_size"] == 8
    bindings = {
        "plan": checker.stable_sha256(plan_path),
        "launch_record": checker.stable_sha256(evidence["launch_record"]),
        "run_log": checker.stable_sha256(evidence["run_log"]),
        "exit_status": checker.stable_sha256(evidence["exit_status"]),
        "gpu_before": gpu_before_binding,
        "gpu_after": checker.stable_sha256(evidence["gpu_after"]),
    }
    completion = checker.audit_completion_record(
        evidence["completion_record"], plan, campaign["state_dir"],
        launch, bindings)
    assert completion["return_code"] == 0
    assert completion["quiescence"]["all_original_identities_gone"] is True


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (
            lambda plan: plan["outputs"][2]["expected_writes"].pop(),
            "strict plan validation failed.*output3",
        ),
        (
            lambda plan: plan["outputs"][2]["expected_endpoint_state"].update(
                {"file_number": 999999}),
            "strict plan validation failed.*output3.*file_number",
        ),
        (
            lambda plan: plan["parameter_snapshots"].update(
                {"source_sha256": "0" * 64}),
            "live source parameter hash",
        ),
        (
            lambda plan: plan["parameter_snapshots"]["output_blocks"].update(
                {"output3": {}}),
            "planned output parameter snapshots",
        ),
    ],
)
def test_real_planner_tamper_is_rejected_before_popen(
        tmp_path: Path, mutate, message: str) -> None:
    import test_suite.test_bbh_segment_plan as plan_tests

    campaign = plan_tests._campaign(tmp_path)
    result = plan_tests._run(campaign)
    assert result.returncode == 0, result.stderr
    plan_path = campaign["output"]
    plan_path.chmod(0o644)
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    mutate(plan)
    plan_path.write_text(
        json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    plan_path.chmod(0o444)
    fake_run = FakeRun(plan["policy"]["ranks"])
    popen_called = False

    def forbidden_popen(*args, **kwargs):
        nonlocal popen_called
        popen_called = True
        raise AssertionError("tampered plan reached Popen")

    runtime = LAUNCHER.Runtime(
        run=fake_run, popen=forbidden_popen,
        statvfs=_ample_statvfs, fstatvfs=_ample_statvfs)
    with pytest.raises(LAUNCHER.LaunchFailure, match=message):
        LAUNCHER.prepare_launch(plan_path, campaign["state_dir"], runtime)
    assert popen_called is False


@pytest.mark.parametrize("tamper", ("c329_write", "endpoint_state"))
def test_real_c322_restart_transition_tamper_is_rejected_before_popen(
        tmp_path: Path, tamper: str) -> None:
    import test_suite.test_bbh_segment_plan as plan_tests

    campaign = plan_tests._campaign(
        tmp_path, cycle=322, time_value=1545.599999999991,
        output3_last_time="1545.599999999991",
        output3_last_write_cycle=322,
        output4_dt="48.0", output4_file_number=45,
        output4_last_time="1540.8000000000002",
        output4_last_write_cycle=322)
    result = plan_tests._run(
        campaign, root_steps="50", target_restart_dt="19.2")
    assert result.returncode == 0, result.stderr
    plan_path = campaign["output"]
    plan_path.chmod(0o644)
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    output4 = next(
        output for output in plan["outputs"] if output["block"] == "output4")
    if tamper == "c329_write":
        output4["expected_writes"] = [
            write for write in output4["expected_writes"]
            if write["cycle"] != 329]
    else:
        output4["expected_endpoint_state"]["file_number"] += 1
    plan_path.write_text(
        json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    plan_path.chmod(0o444)
    popen_called = False

    def forbidden_popen(*args, **kwargs):
        nonlocal popen_called
        popen_called = True
        raise AssertionError("tampered restart schedule reached Popen")

    runtime = LAUNCHER.Runtime(
        run=FakeRun(plan["policy"]["ranks"]), popen=forbidden_popen,
        statvfs=_ample_statvfs, fstatvfs=_ample_statvfs)
    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="strict plan validation failed.*output4"):
        LAUNCHER.prepare_launch(plan_path, campaign["state_dir"], runtime)
    assert popen_called is False


def test_live_source_dcycle_is_rejected_before_popen(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import test_suite.test_bbh_segment_plan as plan_tests

    campaign = plan_tests._campaign(tmp_path)
    result = plan_tests._run(campaign)
    assert result.returncode == 0, result.stderr
    plan = json.loads(campaign["output"].read_text(encoding="utf-8"))
    real_reader = LAUNCHER.read_restart_metadata

    def reader_with_dcycle(path):
        metadata = real_reader(path)
        metadata.parameters["output3"]["dcycle"] = "1"
        return metadata

    monkeypatch.setattr(LAUNCHER, "read_restart_metadata", reader_with_dcycle)
    fake_run = FakeRun(plan["policy"]["ranks"])
    runtime = LAUNCHER.Runtime(
        run=fake_run, statvfs=_ample_statvfs, fstatvfs=_ample_statvfs)
    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="output3 must use dt.*must not serialize dcycle"):
        LAUNCHER.prepare_launch(
            campaign["output"], campaign["state_dir"], runtime)


def test_output3_file_number_limit_is_rejected_before_popen(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import test_suite.test_bbh_segment_plan as plan_tests

    campaign = plan_tests._campaign(tmp_path)
    result = plan_tests._run(campaign)
    assert result.returncode == 0, result.stderr
    plan = json.loads(campaign["output"].read_text(encoding="utf-8"))
    real_reader = LAUNCHER.read_restart_metadata

    def reader_near_limit(path):
        metadata = real_reader(path)
        metadata.parameters["output3"]["file_number"] = "99950"
        plan["source"]["parameters"]["output3"]["file_number"] = "99950"
        return metadata

    monkeypatch.setattr(LAUNCHER, "read_restart_metadata", reader_near_limit)
    fake_run = FakeRun(plan["policy"]["ranks"])
    runtime = LAUNCHER.Runtime(
        run=fake_run, statvfs=_ample_statvfs, fstatvfs=_ample_statvfs)
    runtime.plan_validator = lambda value: {"fixture": "plan-validated"}
    with pytest.raises(LAUNCHER.LaunchFailure,
                       match="file_number reaches.*five-digit limit"):
        runtime.source_validator(campaign["restart"], plan)
