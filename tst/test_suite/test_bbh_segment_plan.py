"""Tests for immutable, cycle-limited AthenaK segment plans."""

from __future__ import annotations

import json
import hashlib
import os
from pathlib import Path
import struct
import subprocess
import sys

import numpy as np

import pytest


ROOT = Path(__file__).resolve().parents[2]
PLANNER = ROOT / "scripts" / "plan_athenak_segment.py"
AMR_MAGIC = 0x41544B414D524331
EVENT_MAGIC = 0x41544B4556543031
CACHE_MAGIC = 0x41544B5343433031


def _git(repo: Path, *arguments: str) -> None:
    subprocess.run(
        ["git", "-C", str(repo), *arguments],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )


def _write_restart(
        path: Path, trajectory: Path, *, real_bytes: int = 8,
        subcycling: str = "level", legacy_subcycling: str = "false",
        legacy_events: str = "false", amr_extension: bool = True,
        event_version: int | None = 2, cache_version: int | None = 1,
        pending_events: tuple[int, ...] | None = None,
        exact_costs: bool = True, capacity: int = 1024,
        root_dt_parameter: str = "4.8", last_dt: float = 4.8,
        output3_dt: str = "4.8", output3_last_time: str = "148.8",
        output3_last_write_cycle: int = 31, cycle: int = 31,
        time_value: float = 148.8,
        output_states: dict[str, dict[str, object]] | None = None,
        nlim: int = -1, tlim: float = 628.8) -> None:
    real_code = "d" if real_bytes == 8 else "f"
    states: dict[str, dict[str, object]] = {
        "output1": {"file_number": 0, "last_time": -1,
                    "last_write_cycle": 31},
        "output2": {"file_number": 4, "last_time": 120,
                    "last_write_cycle": 25},
        "output3": {"file_number": 9, "last_time": output3_last_time,
                    "last_write_cycle": output3_last_write_cycle},
        "output4": {"file_number": 9, "last_time": 148.8,
                    "last_write_cycle": 31},
        "output5": {"file_number": 4, "last_time": 120,
                    "last_write_cycle": 25},
        "output6": {"file_number": 0, "last_time": -1,
                    "last_write_cycle": 31},
    }
    if output_states is not None:
        states = output_states

    def output_value(block: str, name: str) -> object:
        return states[block][name]

    parameters = f"""<job>
basename = effective_bbh_4pn_V100Q
<mesh>
nghost = 4
nx1 = 16
x1min = -1
x1max = 1
nx2 = 16
x2min = -1
x2max = 1
nx3 = 16
x3min = -1
x3max = 1
<meshblock>
nx1 = 16
nx2 = 16
nx3 = 16
<mesh_refinement>
refinement = adaptive
num_levels = 2
max_nmb_per_rank = {capacity}
<mhd>
eos = ideal
nscalars = 0
<adm>
dynamic = true
<time>
subcycling = {subcycling}
allow_legacy_subcycling_restart = {legacy_subcycling}
allow_legacy_ghost_event_counters = {legacy_events}
root_dt_max = {root_dt_parameter}
nlim = {nlim}
tlim = {tlim!r}
restart_dt_growth = 4.8
<problem>
trajectory_file = {trajectory}
<output1>
file_type = hst
dcycle = 1
file_number = {output_value('output1', 'file_number')}
last_time = {output_value('output1', 'last_time')}
last_write_cycle = {output_value('output1', 'last_write_cycle')}
<output2>
file_type = bin
variable = mhd_w_bcc
dt = 48.0
file_number = {output_value('output2', 'file_number')}
last_time = {output_value('output2', 'last_time')}
last_write_cycle = {output_value('output2', 'last_write_cycle')}
ghost_zones = false
<output3>
file_type = bin
variable = mhd_divb
dt = {output3_dt}
file_number = {output_value('output3', 'file_number')}
last_time = {output_value('output3', 'last_time')}
last_write_cycle = {output_value('output3', 'last_write_cycle')}
ghost_zones = false
<output4>
file_type = rst
dt = 19.2
file_number = {output_value('output4', 'file_number')}
last_time = {output_value('output4', 'last_time')}
last_write_cycle = {output_value('output4', 'last_write_cycle')}
single_file_per_rank = false
<output5>
file_type = bin
variable = mhd_gr_diagnostics
dt = 48.0
file_number = {output_value('output5', 'file_number')}
last_time = {output_value('output5', 'last_time')}
last_write_cycle = {output_value('output5', 'last_write_cycle')}
ghost_zones = false
<output6>
file_type = log
dcycle = 1
file_number = {output_value('output6', 'file_number')}
last_time = {output_value('output6', 'last_time')}
last_write_cycle = {output_value('output6', 'last_write_cycle')}
<par_end>
""".encode("utf-8")
    region_size = (-1.0, -1.0, -1.0, 1.0, 1.0, 1.0,
                   0.125, 0.125, 0.125)
    indices = (4, 16, 16, 16, 4, 19, 4, 19, 4, 19,
               8, 8, 8, 4, 11, 4, 11, 4, 11)
    locations = [
        (lx1, lx2, lx3, 1)
        for lx3 in range(2)
        for lx2 in range(2)
        for lx1 in range(2)
    ]
    costs = [2.0] * 8
    if not exact_costs:
        costs[-1] = 1.0
    events = pending_events if pending_events is not None else (0,) * 10
    if len(events) != 10:
        raise ValueError("fixture pending event tuple must contain ten sums")

    content = bytearray(parameters)
    content.extend(struct.pack("<ii", 8, 0))
    content.extend(struct.pack(f"<9{real_code}", *region_size))
    content.extend(struct.pack("<19i", *indices))
    content.extend(struct.pack("<19i", *indices))
    content.extend(struct.pack(
        f"<{real_code}{real_code}i", time_value, last_dt, cycle))
    for location in locations:
        content.extend(struct.pack("<4i", *location))
    for cost in costs:
        content.extend(struct.pack("<f", cost))
    if amr_extension:
        content.extend(struct.pack("<Qii", AMR_MAGIC, 1, 8))
        content.extend(struct.pack("<8i", 0, 0, 1, 1, 2, 2, 3, 3))
    if event_version is not None:
        content.extend(struct.pack("<Qii", EVENT_MAGIC, event_version, 10))
        content.extend(struct.pack("<10Q", *events))
        content.extend(struct.pack("<i", 0))
    if cache_version is not None:
        content.extend(struct.pack("<Qi", CACHE_MAGIC, cache_version))
    stored_extent = 16 + 2 * 4
    cells = stored_extent ** 3
    values_per_block = (
        5 * cells
        + (stored_extent + 1) * stored_extent * stored_extent
        + stored_extent * (stored_extent + 1) * stored_extent
        + stored_extent * stored_extent * (stored_extent + 1)
        + 17 * cells
    )
    field_bytes = np.zeros(values_per_block * len(locations), dtype="<f8").tobytes()
    content.extend(struct.pack("<Q", values_per_block * real_bytes))
    if real_bytes == 8:
        content.extend(field_bytes)
    else:
        content.extend(np.zeros(values_per_block * len(locations),
                                dtype="<f4").tobytes())
    path.write_bytes(content)


def _campaign(tmp_path: Path, **restart_options: object) -> dict[str, Path]:
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "config", "user.email", "segment-test@example.invalid")
    _git(repo, "config", "user.name", "Segment Test")
    binary = repo / "athena"
    binary.write_bytes(b"test-athenak-binary")
    binary.chmod(0o755)
    launcher = repo / "mpirun"
    launcher.write_bytes(b"test-mpirun")
    launcher.chmod(0o755)
    nvidia_smi = tmp_path / "nvidia-smi"
    nvidia_smi.write_bytes(b"test-nvidia-smi")
    nvidia_smi.chmod(0o755)
    trajectory = repo / "trajectory.dat"
    trajectory.write_text("0 0 0\n", encoding="ascii")
    restart = repo / "source.rst"
    _write_restart(restart, trajectory.resolve(), **restart_options)
    source_history = repo / "source.user.hst"
    source_history.write_text(
        "# Athena++ history data\n"
        "# [1]=time [2]=dt [3]=baryon_m\n"
        "1.4880000000000000e2 4.8 1.0e2\n", encoding="utf-8")
    _git(repo, "add", "athena", "mpirun", "trajectory.dat", "source.rst",
         "source.user.hst")
    _git(repo, "commit", "-q", "-m", "fixture")
    state_dir = tmp_path / "state"
    state_dir.mkdir()
    staging_dir = tmp_path / "staging"
    staging_dir.mkdir()
    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    return {
        "repo": repo,
        "binary": binary,
        "trajectory": trajectory,
        "launcher": launcher,
        "nvidia_smi": nvidia_smi,
        "restart": restart,
        "source_history": source_history,
        "state_dir": state_dir,
        "staging_dir": staging_dir,
        "evidence_dir": evidence_dir,
        "output": evidence_dir / "segment.plan.json",
    }


def _command(campaign: dict[str, Path], **overrides: str) -> list[str]:
    values = {
        "root_steps": "100",
        "root_dt": "4.8",
        "guard_steps": "2",
        "ranks": "8",
        "wall": "28800",
        "planned_peak_output_gib": "200",
        **overrides,
    }
    command = [
        sys.executable,
        str(PLANNER),
        "--repo", str(campaign["repo"]),
        "--source-restart", str(campaign["restart"]),
        "--source-history", str(campaign["source_history"]),
        "--anchor",
        "--binary", str(campaign["binary"]),
        "--trajectory", str(campaign["trajectory"]),
        "--root-steps", values["root_steps"],
        "--root-dt", values["root_dt"],
        "--tlim-guard-steps", values["guard_steps"],
        "--ranks", values["ranks"],
        "--launcher", str(campaign["launcher"]),
        "--state-dir", str(campaign["state_dir"]),
        "--staging-dir", str(campaign["staging_dir"]),
        "--evidence-dir", str(campaign["evidence_dir"]),
        "--wall-time-seconds", values["wall"],
        "--planned-peak-output-gib", values["planned_peak_output_gib"],
        "--required-unnumbered", "effective_bbh_4pn_V100Q.mhd.hst",
        "--required-unnumbered", "effective_bbh_4pn_V100Q.user.hst",
        "--required-unnumbered", "effective_bbh_4pn_V100Q.log",
        "--output", str(campaign["output"]),
    ]
    if "target_capacity" in overrides:
        command.extend([
            "--target-max-nmb-per-rank", values["target_capacity"],
        ])
    return command


def _run(campaign: dict[str, Path], **overrides: str) -> subprocess.CompletedProcess[str]:
    environment = dict(os.environ)
    environment["PATH"] = (
        f"{campaign['nvidia_smi'].parent}{os.pathsep}"
        f"{environment.get('PATH', '')}"
    )
    return subprocess.run(
        _command(campaign, **overrides),
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=environment,
    )


def test_binary64_endpoint_cycle_guard_and_snapshots(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)

    result = _run(campaign)

    assert result.returncode == 0, result.stderr
    plan = json.loads(campaign["output"].read_text(encoding="utf-8"))
    assert plan["expected"] == {
        "source_cycle": 31,
        "source_time": 148.8,
        "root_steps": 100,
        "root_dt": 4.8,
        "tlim_guard_steps": 2,
        "final_cycle": 131,
        "final_time": 628.7999999999997,
        "tlim": 638.3999999999996,
    }
    assert plan["command_overrides"] == {
        "time/nlim": 131,
        "time/tlim": 638.3999999999996,
        "output3/dt": 4.8,
    }
    assert plan["capacity_transition"]["kind"] == "unchanged_v1"
    assert plan["capacity_transition"]["source_max_nmb_per_rank"] == 1024
    assert plan["capacity_transition"]["target_max_nmb_per_rank"] == 1024
    assert plan["capacity_transition"]["runtime_override"] is None
    assert plan["source_qualification"]["mode"] == "anchor_full_audit"
    assert plan["source_qualification"]["audit"]["valid"] is True
    assert plan["source_qualification"]["audit"]["stored_reals"][
        "nonfinite_count"] == 0
    assert plan["source_qualification"]["source_baryon_mass"]["value"] == 100.0
    assert plan["launch_contract"]["mpi_argv"] == [
        str(campaign["launcher"].resolve()), "--allow-run-as-root", "--bind-to",
        "none", "-np", "8",
    ]
    restart_size = plan["inputs"]["source_restart"]["size"]
    trajectory_size = plan["inputs"]["trajectory"]["size"]
    assert plan["launch_contract"]["disk_preflight"] == {
        "kind": "statvfs_unique_filesystem_budget_v1",
        "accounting": "group_roles_by_st_dev_once_v1",
        "formula": (
            "per_filesystem_required_free_bytes=max("
            "additional_hard_minimum_free_bytes,max(minimum_reserve_bytes,"
            "minimum_reserve_restart_multiples*source_restart_size_bytes)+"
            "sum(role_contribution_bytes))"
        ),
        "used_percent_exclusive_max": 75,
        "minimum_reserve_bytes": 50 * (1 << 30),
        "minimum_reserve_restart_multiples": 2,
        "additional_hard_minimum_free_bytes": 180 * (1 << 30),
        "source_restart_size_bytes": restart_size,
        "source_restart_staging_bytes": restart_size,
        "trajectory_staging_bytes": trajectory_size,
        "staging_copy_bytes": restart_size + trajectory_size,
        "planned_peak_output_gib": 200,
        "planned_peak_output_bytes": 200 * (1 << 30),
        "role_contributions_bytes": {
            "state_dir": {"planned_peak_output_bytes": 200 * (1 << 30)},
            "staging_dir": {
                "source_restart_staging_bytes": restart_size,
                "trajectory_staging_bytes": trajectory_size,
            },
        },
        "bound_directory_fds": {"state_dir": 202, "staging_dir": 206},
    }
    assert plan["launch_contract"]["athena_argv_template"] == [
        "/proc/{holder_pid}/fd/205", "--kokkos-map-device-id-by=mpi_rank",
        "-r", "/proc/{holder_pid}/fd/200",
        "-d", "/proc/{holder_pid}/fd/202", "-t", "08:00:00",
        "time/nlim=131", "time/tlim=638.3999999999996",
        "problem/trajectory_file=/proc/{holder_pid}/fd/201",
        "output3/dt=4.8",
    ]
    assert plan["policy"]["segment_termination"]["primary"] == "cycle_limit"
    assert plan["launch_contract"]["plan_path"] == str(campaign["output"])
    assert plan["launch_contract"]["directory_transport"] == {
        "kind": "linux_proc_holder_dirfd_v1",
        "holder_pid_token": "{holder_pid}",
        "roles": {
            "state_dir": {
                "role": "state_dir", "fd": 202,
                "planned_path": str(campaign["state_dir"]),
                "proc_path_template": "/proc/{holder_pid}/fd/202",
            },
            "evidence_dir": {
                "role": "evidence_dir", "fd": 203,
                "planned_path": str(campaign["evidence_dir"]),
                "proc_path_template": "/proc/{holder_pid}/fd/203",
            },
        },
    }


    assert plan["launch_contract"]["executable_transport"] == {
        "kind": "linux_proc_holder_execfd_v1",
        "holder_pid_token": "{holder_pid}",
        "roles": {
            "launcher": {
                "role": "launcher", "fd": 204,
                "parent_path": str(campaign["launcher"].resolve()),
                "proc_path_template": "/proc/{holder_pid}/fd/204",
            },
            "binary": {
                "role": "binary", "fd": 205,
                "parent_path": str(campaign["binary"].resolve()),
                "proc_path_template": "/proc/{holder_pid}/fd/205",
            },
        },
    }
    environment = plan["launch_contract"]["environment"]
    assert environment["kind"] == "explicit_values_with_rank_projection_v2"
    assert environment["values"]["LANG"] == environment["values"]["LC_ALL"] == "C"
    assert environment["values"]["CUDA_DEVICE_ORDER"] == "PCI_BUS_ID"
    assert environment["values"]["PRTE_MCA_schizo_proxy"] == "ompi"
    assert environment["sha256"] == hashlib.sha256(json.dumps(
        environment["values"], sort_keys=True, separators=(",", ":"),
        ensure_ascii=True, allow_nan=False).encode("utf-8")).hexdigest()
    projection = environment["rank_projection"]
    assert projection["kind"] == "prrte_consumed_projection_v1"
    assert projection["inherited_values"] == {
        key: environment["values"][key]
        for key in ("HOME", "LANG", "LC_ALL", "CUDA_DEVICE_ORDER")
    }
    assert projection["consumed_absent"] == ["PRTE_MCA_schizo_proxy"]
    assert projection["sha256"] == hashlib.sha256(json.dumps({
        "inherited_values": projection["inherited_values"],
        "consumed_absent": projection["consumed_absent"],
    }, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False).encode("utf-8")).hexdigest()
    assert plan["tools"]["nvidia_smi"]["path"] == str(
        campaign["nvidia_smi"].resolve())
    assert plan["policy"]["endpoint_time_ulp_tolerance"] == 0
    assert plan["policy"]["fixed_root_dt_mode"]["enabled"] is True
    assert plan["source"]["real_bytes"] == 8
    assert plan["source"]["event_counter_version"] == 2
    assert plan["source"]["restart_cache_contract_version"] == 1
    assert plan["source"]["level_subcycling_costs_exact"] is True
    assert plan["inputs"]["repo"]["clean"] is True
    assert len(plan["inputs"]["repo"]["commit"]) == 40
    assert all(len(plan["inputs"][name]["sha256"]) == 64 for name in (
        "binary", "trajectory", "source_restart"
    ))
    outputs = {row["block"]: row for row in plan["outputs"]}
    assert len(outputs["output1"]["expected_writes"]) == 100
    assert outputs["output2"]["expected_writes"][-1] == {
        "cycle": 131, "time": 628.7999999999997,
        "kind": "forced_final", "file_number": 14,
    }
    assert outputs["output4"]["expected_writes"][-1] == {
        "cycle": 131, "time": 628.7999999999997,
        "kind": "restart_final_rewrite", "file_number": 34,
    }
    assert outputs["output4"]["expected_endpoint_state"] == {
        "file_number": 35, "last_time": 628.8000000000001,
        "last_write_cycle": 131,
    }
    assert outputs["output2"]["relative_path_template"] == (
        "bin/effective_bbh_4pn_V100Q.mhd_w_bcc.{file_number:05d}.bin"
    )
    assert outputs["output4"]["relative_path_template"] == (
        "rst/effective_bbh_4pn_V100Q.{file_number:05d}.rst"
    )
    assert set(outputs["output1"]["required_unnumbered_paths"]) == {
        "effective_bbh_4pn_V100Q.mhd.hst",
        "effective_bbh_4pn_V100Q.user.hst",
    }
    assert plan["required_files"] == [
        "effective_bbh_4pn_V100Q.log",
        "effective_bbh_4pn_V100Q.mhd.hst",
        "effective_bbh_4pn_V100Q.user.hst",
    ]
    assert plan["source"]["parameters"]["job"]["basename"] == (
        "effective_bbh_4pn_V100Q"
    )
    assert campaign["output"].stat().st_mode & 0o222 == 0


def test_plan_binds_only_a_strict_capacity_increase(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)

    result = _run(campaign, target_capacity="1280")

    assert result.returncode == 0, result.stderr
    plan = json.loads(campaign["output"].read_text(encoding="utf-8"))
    transition = plan["capacity_transition"]
    assert transition["kind"] == "increase_v1"
    assert transition["source_max_nmb_per_rank"] == 1024
    assert transition["target_max_nmb_per_rank"] == 1280
    assert transition["runtime_override"] == \
        "mesh_refinement/max_nmb_per_rank=1280"
    assert transition["gpu_memory_model"] == {
        "kind": "fixed_conservative_per_meshblock_slot_v1",
        "mib_per_slot_numerator": 1433,
        "mib_per_slot_denominator": 100,
        "usable_fraction_numerator": 4,
        "usable_fraction_denominator": 5,
        "required_per_rank_memory_mib_ceiling": 18343,
        "minimum_gpu_memory_total_mib": 22928,
    }
    assert plan["command_overrides"][
        "mesh_refinement/max_nmb_per_rank"] == 1280
    assert plan["launch_contract"]["athena_argv_template"][-1] == \
        "mesh_refinement/max_nmb_per_rank=1280"
    assert plan["policy"]["mutable_parameters"] == [
        "time/nlim", "time/tlim", "time/restart_dt_growth",
        "output*/file_number", "output*/last_time",
        "output*/last_write_cycle",
    ]


@pytest.mark.parametrize("target", ("1024", "1000", "16385"))
def test_plan_rejects_equal_decreasing_or_excess_capacity_target(
        tmp_path: Path, target: str) -> None:
    campaign = _campaign(tmp_path)

    result = _run(campaign, target_capacity=target)

    assert result.returncode != 0
    assert "target-max-nmb-per-rank" in result.stderr


def test_plan_rejects_duplicate_capacity_target(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    command = _command(campaign, target_capacity="1280")
    command.extend(["--target-max-nmb-per-rank", "1536"])
    environment = dict(os.environ)
    environment["PATH"] = (
        f"{campaign['nvidia_smi'].parent}{os.pathsep}"
        f"{environment.get('PATH', '')}"
    )

    result = subprocess.run(
        command, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, env=environment)

    assert result.returncode != 0
    assert "may be specified only once" in result.stderr


def test_plan_capacity_increase_rescues_near_full_source(tmp_path: Path) -> None:
    # Eight blocks over eight ranks exactly fill a serialized capacity of one.
    campaign = _campaign(tmp_path, capacity=1)

    result = _run(campaign, target_capacity="65")

    assert result.returncode == 0, result.stderr
    plan = json.loads(campaign["output"].read_text(encoding="utf-8"))
    assert plan["source"]["partition"]["minimum_capacity_headroom"] == 0
    assert plan["capacity_transition"]["target_max_nmb_per_rank"] == 65


def test_plan_output_must_be_canonical_direct_evidence_child(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    campaign["output"] = tmp_path / "wrong.plan.json"
    result = _run(campaign)
    assert result.returncode != 0
    assert "canonical direct evidence child" in result.stderr
    assert not campaign["output"].exists()


@pytest.mark.parametrize(
    ("restart_options", "message"),
    [
        ({"real_bytes": 4}, "not Real8"),
        ({"subcycling": "none"}, "strict level subcycling"),
        ({"legacy_subcycling": "true"}, "allow_legacy_subcycling_restart"),
        ({"legacy_events": "true"}, "allow_legacy_ghost_event_counters"),
        ({"amr_extension": False}, "AMR counter extension"),
        ({"event_version": 1}, "required v2 event-counter"),
        ({"event_version": None}, "required v2 event-counter"),
        ({"cache_version": None}, "v1 restart-cache contract"),
        ({"pending_events": (0, 0, 0, 0, 1, 0, 0, 0, 0, 0)},
         "nonzero pending event counters"),
        ({"exact_costs": False}, "exact level costs"),
        ({"capacity": 64}, "capacity headroom"),
        ({"root_dt_parameter": "4.800000000000001"}, "root_dt_max"),
        ({"last_dt": 2.4}, "source last_dt"),
    ],
)
def test_metadata_contract_rejections(
        tmp_path: Path, restart_options: dict[str, object], message: str) -> None:
    campaign = _campaign(tmp_path, **restart_options)

    result = _run(campaign)

    assert result.returncode != 0
    assert message in result.stderr
    assert not campaign["output"].exists()


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"root_steps": "0"}, "--root-steps must be positive"),
        ({"root_dt": "nan"}, "--root-dt must be finite and positive"),
        ({"root_dt": "inf"}, "--root-dt must be finite and positive"),
        ({"guard_steps": "0"}, "--tlim-guard-steps must be positive"),
        ({"ranks": "0"}, "--ranks must be positive"),
        ({"wall": "0"}, "--wall-time-seconds must be positive"),
        ({"planned_peak_output_gib": "0"},
         "--planned-peak-output-gib must be positive"),
    ],
)
def test_invalid_plan_parameters_are_rejected(
        tmp_path: Path, overrides: dict[str, str], message: str) -> None:
    campaign = _campaign(tmp_path)

    result = _run(campaign, **overrides)

    assert result.returncode != 0
    assert message in result.stderr
    assert not campaign["output"].exists()


def test_explicit_planned_peak_output_is_bound_as_integer_bytes(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)

    result = _run(campaign, planned_peak_output_gib="110")

    assert result.returncode == 0, result.stderr
    contract = json.loads(campaign["output"].read_text(
        encoding="utf-8"))["launch_contract"]["disk_preflight"]
    assert contract["planned_peak_output_gib"] == 110
    assert contract["planned_peak_output_bytes"] == 110 * (1 << 30)


def test_existing_plan_is_never_overwritten(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    campaign["output"].write_text("sentinel\n", encoding="ascii")

    result = _run(campaign)

    assert result.returncode != 0
    assert "refusing to overwrite existing plan" in result.stderr
    assert campaign["output"].read_text(encoding="ascii") == "sentinel\n"


def test_dirty_repository_is_rejected(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    (campaign["repo"] / "untracked").write_text("dirty", encoding="ascii")

    result = _run(campaign)

    assert result.returncode != 0
    assert "repository is not clean" in result.stderr
    assert not campaign["output"].exists()


@pytest.mark.parametrize(
    ("seconds", "token"),
    (("1", "00:00:01"), ("360001", "100:00:01")),
)
def test_wall_limit_is_serialized_as_strict_hhmmss(
        tmp_path: Path, seconds: str, token: str) -> None:
    campaign = _campaign(tmp_path)
    result = _run(campaign, wall=seconds)
    assert result.returncode == 0, result.stderr
    plan = json.loads(campaign["output"].read_text(encoding="utf-8"))
    argv = plan["launch_contract"]["athena_argv_template"]
    assert argv[argv.index("-t") + 1] == token


def test_output3_replay_accepts_two_root_steps_of_serialized_phase_lag(
        tmp_path: Path) -> None:
    """The cadence phase is restart state, not a proxy for physical time."""

    campaign = _campaign(
        tmp_path, output3_dt="19.2", output3_last_time="139.2",
        output3_last_write_cycle=28)
    result = _run(campaign)
    assert result.returncode == 0, result.stderr
    plan = json.loads(campaign["output"].read_text(encoding="utf-8"))
    output3 = next(row for row in plan["outputs"] if row["block"] == "output3")
    assert [row["cycle"] for row in output3["expected_writes"]] == list(
        range(32, 132))
    assert output3["expected_endpoint_state"] == {
        "file_number": 109,
        "last_time": 619.1999999999998,
        "last_write_cycle": 131,
    }


def test_output3_replay_rejects_phase_that_would_skip_a_root_cycle(
        tmp_path: Path) -> None:
    campaign = _campaign(
        tmp_path, output3_dt="19.2", output3_last_time="151.2",
        output3_last_write_cycle=31)
    result = _run(campaign)
    assert result.returncode != 0
    assert "exactly one scheduled" in result.stderr
    assert not campaign["output"].exists()


def test_anchor_can_rebind_disabled_output3_dt_to_root_cadence(
        tmp_path: Path) -> None:
    campaign = _campaign(
        tmp_path, output3_dt="0", output3_last_time="139.2",
        output3_last_write_cycle=28)
    result = _run(campaign)
    assert result.returncode == 0, result.stderr
    plan = json.loads(campaign["output"].read_text(encoding="utf-8"))
    output3 = next(row for row in plan["outputs"] if row["block"] == "output3")
    assert output3["parameters"]["dt"] == "4.8"


def _make_parent_pass(campaign: dict[str, Path], tmp_path: Path,
                      mutation=None, *, immutable: bool = True) -> Path:
    anchor_result = _run(campaign)
    assert anchor_result.returncode == 0, anchor_result.stderr
    anchor = json.loads(campaign["output"].read_text(encoding="utf-8"))
    parent = {
        "schema": 1,
        "kind": "athenak_segment_pass",
        "status": "pass",
        "qualification_mode": "complete_segment_v1",
        "expected": {"final_cycle": 31, "final_time": 148.8},
        "bindings": {
            "endpoint_restart": anchor["inputs"]["source_restart"],
            "trajectory": anchor["inputs"]["trajectory"],
        },
        "endpoint_restart_audit": anchor["source_qualification"]["audit"],
        "scientific_threshold_audit": {"baryon_mass": {"last": 100.0}},
        "scientific_advisories": {
            "schema": "athenak_scientific_advisories_v1",
            "severity": "green",
            "floor_rates": {},
            "pass_fail_effect": "none_yellow_advisories_are_nonfatal",
        },
        "output_inventory": [{
            **anchor["source_qualification"]["source_baryon_mass"]["evidence"],
            "history": {
                "columns": ["time", "dt", "baryon_m"],
                "column_last": {"time": 148.8, "dt": 4.8,
                                "baryon_m": 100.0},
            },
        }],
    }
    if mutation is not None:
        mutation(parent)
    parent_path = tmp_path / "previous.pass.ready"
    parent_path.write_text(json.dumps(parent), encoding="utf-8")
    if immutable:
        parent_path.chmod(0o444)
    next_state = tmp_path / "next-state"
    next_staging = tmp_path / "next-staging"
    next_evidence = tmp_path / "next-evidence"
    for directory in (next_state, next_staging, next_evidence):
        directory.mkdir()
    campaign["state_dir"] = next_state
    campaign["staging_dir"] = next_staging
    campaign["evidence_dir"] = next_evidence
    campaign["output"] = next_evidence / "segment.plan.json"
    return parent_path


def _run_with_parent(campaign: dict[str, Path], parent: Path
                     ) -> subprocess.CompletedProcess[str]:
    command = _command(campaign)
    anchor_index = command.index("--anchor")
    command[anchor_index:anchor_index + 1] = [
        "--parent-segment-pass", str(parent)]
    environment = dict(os.environ)
    environment["PATH"] = (
        f"{campaign['nvidia_smi'].parent}{os.pathsep}"
        f"{environment.get('PATH', '')}"
    )
    return subprocess.run(
        command, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, env=environment)


def _make_recovery_parent_pass(
        campaign: dict[str, Path], tmp_path: Path, mutation=None) -> Path:
    anchor_result = _run(campaign)
    assert anchor_result.returncode == 0, anchor_result.stderr
    anchor = json.loads(campaign["output"].read_text(encoding="utf-8"))
    recovery_state = tmp_path / "recovered-state"
    recovery_evidence = tmp_path / "recovered-evidence"
    recovery_state.mkdir(mode=0o700)
    recovery_evidence.mkdir(mode=0o700)
    recovery_state.chmod(0o700)
    recovery_evidence.chmod(0o700)
    derived_history = recovery_state / "effective_bbh_4pn_V100Q.user.hst"
    derived_history.write_bytes(campaign["source_history"].read_bytes())
    derived_history.chmod(0o444)

    def record(path: Path) -> dict[str, object]:
        return {
            "path": str(path.resolve()), "size": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "closure_check": "linux_proc_fd",
        }

    endpoint = anchor["inputs"]["source_restart"]
    selected = {
        "block": "output4", "relative_path": "rst/run.00009.rst",
        "expected_write": {"cycle": 31, "time": 148.8,
                           "kind": "scheduled", "file_number": 9},
        "binding": endpoint,
    }
    logical = {
        "source_cycle": 1, "source_time": 4.8,
        "final_cycle": 31, "final_time": 148.8, "root_steps": 30,
        "execute_only": True, "expected_numbered_paths": [],
        "required_unnumbered_paths": [],
        "all_expected_prefix_artifacts_present": True,
    }
    lifecycle = {
        "kind": "same_boot_closed_processes_v1",
        "all_original_identities_gone": True,
    }
    run_log_prefix = {
        "root_step_diagnostics": {
            "cycle_max": 31,
            "all_recovered_prefix_cycles_present": True,
        },
        "cache": {"solver_failures": 0, "nonfinite_proposed_values": 0},
    }
    history_binding = record(derived_history)
    candidates = [{
        "relative_path": "rst/run.00010.rst",
        "expected_write": {"cycle": 40},
        "classification": "incomplete_truncated",
    }, {
        "relative_path": selected["relative_path"],
        "expected_write": selected["expected_write"],
        "classification": "complete",
    }]
    original_plan = recovery_evidence / "segment.plan.json"
    original_plan.write_text(json.dumps({
        "policy": {"scheduled_prefix_recovery": {
            "kind": "scheduled_prefix_recovery_v1",
            "candidate_kind": "original_plan_scheduled_restart_only",
            "selection": "highest_complete_candidate",
            "later_complete_scientific_failure": "fatal_no_fallback",
            "cadence_replay": "execute_only_binary64_v1",
            "original_trees": "read_only_no_delete",
            "derived_text": "exact_original_byte_prefix",
            "unknown_state_nodes": "fatal",
            "lifecycle": "nonzero_completion_or_same_boot_closed_processes",
        }},
        "outputs": [{
            "block": "output4", "file_type": "rst",
            "relative_path_template": "rst/run.{file_number:05d}.rst",
            "expected_writes": [selected["expected_write"]],
        }],
    }), encoding="utf-8")
    original_plan.chmod(0o444)
    original_plan_binding = record(original_plan)
    recovery_record_payload = {
        "schema": 1, "kind": "athenak_segment_prefix_recovery",
        "status": "pass", "qualification_mode":
            "scheduled_prefix_recovery_v1",
        "policy": {
            "kind": "scheduled_prefix_recovery_v1",
            "candidate_kind": "original_plan_scheduled_restart_only",
            "selection": "highest_complete_candidate",
            "later_complete_scientific_failure": "fatal_no_fallback",
            "cadence_replay": "execute_only_binary64_v1",
            "original_trees": "read_only_no_delete",
            "derived_text": "exact_original_byte_prefix",
            "unknown_state_nodes": "fatal",
            "lifecycle": "nonzero_completion_or_same_boot_closed_processes",
        },
        "original_plan": original_plan_binding,
        "selected_scheduled_restart": selected,
        "logical_prefix": logical,
        "lifecycle": lifecycle,
        "run_log_prefix_audit": run_log_prefix,
        "candidate_inventory": candidates,
        "recovery_state_dir": str(recovery_state.resolve()),
    }
    recovery_record = recovery_evidence / "segment.prefix.recovery.ready"
    recovery_record.write_text(json.dumps(recovery_record_payload), encoding="utf-8")
    recovery_record.chmod(0o444)
    recovery_record_binding = record(recovery_record)
    parent = {
        "schema": 1, "kind": "athenak_segment_pass", "status": "pass",
        "qualification_mode": "scheduled_prefix_recovery_v1",
        "expected": {"final_cycle": 31, "final_time": 148.8},
        "bindings": {
            "endpoint_restart": endpoint,
            "trajectory": anchor["inputs"]["trajectory"],
            "recovery_record": recovery_record_binding,
            "original_plan": original_plan_binding,
        },
        "endpoint_restart_audit": anchor["source_qualification"]["audit"],
        "scientific_threshold_audit": {"baryon_mass": {"last": 100.0}},
        "scientific_advisories": {
            "schema": "athenak_scientific_advisories_v1",
            "severity": "green", "floor_rates": {},
            "pass_fail_effect": "none_yellow_advisories_are_nonfatal",
        },
        "output_inventory": [{
            **history_binding,
            "history": {
                "columns": ["time", "dt", "baryon_m"],
                "column_last": {"time": 148.8, "dt": 4.8,
                                "baryon_m": 100.0},
            },
        }],
        "recovery_provenance": {
            "kind": "scheduled_prefix_recovery_v1",
            "policy": recovery_record_payload["policy"],
            "record": recovery_record_binding,
            "original_plan": original_plan_binding,
            "selected_scheduled_restart": selected,
            "logical_prefix": logical,
            "lifecycle": lifecycle,
            "run_log_prefix_audit": run_log_prefix,
            "derived_text_prefixes": [{"derived": history_binding}],
            "original_trees_unchanged": True,
        },
    }
    if mutation is not None:
        mutation(parent)
    parent_path = recovery_evidence / "segment.prefix.pass.ready"
    parent_path.write_text(json.dumps(parent), encoding="utf-8")
    parent_path.chmod(0o444)
    campaign["source_history"] = derived_history
    next_state = tmp_path / "next-recovery-state"
    next_staging = tmp_path / "next-recovery-staging"
    next_evidence = tmp_path / "next-recovery-evidence"
    for directory in (next_state, next_staging, next_evidence):
        directory.mkdir()
    campaign["state_dir"] = next_state
    campaign["staging_dir"] = next_staging
    campaign["evidence_dir"] = next_evidence
    campaign["output"] = next_evidence / "segment.plan.json"
    return parent_path


def test_planner_rejects_handwritten_legacy_recovery_parent(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    parent = _make_recovery_parent_pass(campaign, tmp_path)

    result = _run_with_parent(campaign, parent)

    assert result.returncode != 0
    assert "recovery parent strict provenance failed" in result.stderr


def test_planner_rejects_recovery_parent_provenance_tamper(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    parent = _make_recovery_parent_pass(
        campaign, tmp_path,
        lambda report: report["recovery_provenance"][
            "selected_scheduled_restart"]["expected_write"].__setitem__(
                "kind", "forced_final"))

    result = _run_with_parent(campaign, parent)

    assert result.returncode != 0
    assert "recovery parent" in result.stderr


def test_parent_segment_pass_chains_endpoint_and_baryon_source(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    parent = _make_parent_pass(campaign, tmp_path)
    result = _run_with_parent(campaign, parent)
    assert result.returncode == 0, result.stderr
    plan = json.loads(campaign["output"].read_text(encoding="utf-8"))
    qualification = plan["source_qualification"]
    assert qualification["mode"] == "parent_segment_pass"
    assert qualification["audit"]["valid"] is True
    assert qualification["source_baryon_mass"]["evidence"]["path"] == str(
        campaign["source_history"].resolve())
    assert qualification["source_baryon_mass"]["value"] == 100.0


def test_legacy_complete_parent_without_mode_remains_chainable(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    parent = _make_parent_pass(
        campaign, tmp_path,
        lambda report: (
            report.pop("qualification_mode"),
            report.__setitem__("completion_record_audit", {"return_code": 0}),
            report.__setitem__("run_log_audit", {"termination": "cycle_limit"}),
        ))

    result = _run_with_parent(campaign, parent)

    assert result.returncode == 0, result.stderr
    plan = json.loads(campaign["output"].read_text(encoding="utf-8"))
    assert plan["source_qualification"]["parent_qualification_mode"] == \
        "legacy_complete_segment_v1"


def test_missing_mode_cannot_disguise_recovery_as_legacy_complete(
        tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    parent = _make_parent_pass(
        campaign, tmp_path,
        lambda report: (
            report.pop("qualification_mode"),
            report.__setitem__("completion_record_audit", {"return_code": 0}),
            report.__setitem__("run_log_audit", {"termination": "cycle_limit"}),
            report.__setitem__("recovery_provenance", {}),
        ))

    result = _run_with_parent(campaign, parent)

    assert result.returncode != 0
    assert "legacy parent" in result.stderr


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        (lambda parent: parent["bindings"]["endpoint_restart"].__setitem__(
            "sha256", "f" * 64), "does not identify"),
        (lambda parent: parent["expected"].__setitem__("final_time", 149.0),
         "planned endpoint"),
        (lambda parent: parent["scientific_threshold_audit"][
            "baryon_mass"].__setitem__("last", 0.0), "baryon mass"),
    ),
)
def test_parent_segment_pass_rejects_broken_chain(
        tmp_path: Path, mutation, message: str) -> None:
    campaign = _campaign(tmp_path)
    parent = _make_parent_pass(campaign, tmp_path, mutation)
    result = _run_with_parent(campaign, parent)
    assert result.returncode != 0
    assert message in result.stderr
    assert not campaign["output"].exists()


def test_parent_segment_pass_must_be_immutable(tmp_path: Path) -> None:
    campaign = _campaign(tmp_path)
    parent = _make_parent_pass(campaign, tmp_path, immutable=False)
    result = _run_with_parent(campaign, parent)
    assert result.returncode != 0
    assert "immutable" in result.stderr
