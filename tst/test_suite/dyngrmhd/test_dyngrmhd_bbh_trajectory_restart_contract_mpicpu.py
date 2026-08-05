"""Restart-contract regression for effective-BBH table trajectories."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
INPUT_FILE = ROOT / "inputs" / "dyngr" / "effective_bbh.athinput"
BASENAME = "bbh_trajectory_restart_contract"
MPI_RANKS = 2
TIMEOUT_SECONDS = 45
CONTRACT_FAILURE = (
    r"dynbbh trajectory or metric parameters differ from the restart contract"
)


def _write_table(path: Path, *, changed_middle_position=False) -> None:
    """Write a safe stationary, equal-mass table spanning every test stencil."""
    rows = []
    for row_index, time in enumerate((-0.5, 0.0, 0.5)):
        x1 = 10.000001 if changed_middle_position and row_index == 1 else 10.0
        state = (
            x1, 0.0, 0.0,
            -10.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
            0.5, 0.5,
        )
        rows.append(" ".join(f"{value:.17e}" for value in (time, *state)))
    path.write_text(
        "# t x1 y1 z1 x2 y2 z2 vx1 vy1 vz1 vx2 vy2 vz2 "
        "a1x a1y a1z a2x a2y a2z m1_full m2_full\n"
        + "\n".join(rows)
        + "\n",
        encoding="utf-8",
    )


def _base_command(run_dir: Path) -> list[str]:
    timeout = shutil.which("timeout")
    assert timeout is not None, "GNU timeout is required by this restart test"
    mpi_launcher = shutil.which("mpirun") or shutil.which("mpiexec")
    if mpi_launcher is None:
        # Some module/spack builds expose the MPI compiler to CMake without adding
        # the matching launcher to PATH.  Use the compiler's sibling launcher so
        # this remains the same MPI implementation that linked ./athena.
        cache_path = Path("..").resolve() / "CMakeCache.txt"
        if cache_path.is_file():
            for line in cache_path.read_text(encoding="utf-8").splitlines():
                if line.startswith("MPI_CXX_COMPILER:FILEPATH="):
                    compiler = Path(line.split("=", 1)[1])
                    for name in ("mpirun", "mpiexec"):
                        candidate = compiler.with_name(name)
                        if candidate.is_file():
                            mpi_launcher = str(candidate)
                            break
                    break
    assert mpi_launcher is not None, "an MPI launcher is required by this restart test"
    return [
        timeout,
        "--signal=TERM",
        "--kill-after=5s",
        f"{TIMEOUT_SECONDS}s",
        mpi_launcher,
        "-np",
        str(MPI_RANKS),
        "./athena",
        "-d",
        str(run_dir),
    ]


def _run(command: list[str], label: str) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    with Path(testutils.LOG_FILE_PATH).open("a", encoding="utf-8") as log_file:
        log_file.write(f"\n{label}: {' '.join(command)}\n")
        log_file.write(result.stdout)
    assert result.returncode not in (124, 137), (
        f"{label} timed out after {TIMEOUT_SECONDS} seconds:\n{result.stdout[-8000:]}"
    )
    return result


def _initial_command(run_dir: Path, trajectory: Path) -> list[str]:
    return _base_command(run_dir) + [
        "-i",
        str(INPUT_FILE),
        f"job/basename={BASENAME}",
        "mesh/nghost=4",
        "mesh/nx1=8",
        "mesh/nx2=8",
        "mesh/nx3=8",
        "meshblock/nx1=4",
        "meshblock/nx2=8",
        "meshblock/nx3=8",
        "mesh_refinement/refinement=adaptive",
        "time/nlim=1",
        "time/tlim=0.01",
        "problem/trajectory_mode=table",
        f"problem/trajectory_file={trajectory.resolve()}",
        "problem/trajectory_time_offset=0.0",
        "output1/dt=0.0",
        "output2/dt=0.0",
    ]


def _restart_command(
    run_dir: Path, checkpoint: Path, label: str, overrides=()
) -> list[str]:
    return _base_command(run_dir) + [
        "-r",
        str(checkpoint.resolve()),
        f"job/basename={BASENAME}_{label}",
        "time/nlim=2",
        "time/tlim=0.02",
        "output3/dt=0.0",
        *overrides,
    ]


def _require_success(result: subprocess.CompletedProcess[str], label: str) -> None:
    assert result.returncode == 0, f"{label} failed:\n{result.stdout[-8000:]}"


def _load_balance_efficiency(
    result: subprocess.CompletedProcess[str], label: str
) -> float:
    matches = re.findall(
        r"load balancing efficiency = ([0-9.eE+-]+) \(this run segment\)",
        result.stdout,
    )
    assert matches, f"{label} did not report segment-local load balance efficiency"
    return float(matches[-1])


def _require_contract_failure(
    result: subprocess.CompletedProcess[str], label: str
) -> None:
    assert result.returncode != 0, f"{label} unexpectedly succeeded"
    assert re.search(CONTRACT_FAILURE, result.stdout, flags=re.IGNORECASE), (
        f"{label} did not fail through the trajectory contract:\n"
        f"{result.stdout[-8000:]}"
    )


def test_table_trajectory_restart_contract_is_path_independent_and_fail_closed(
    tmp_path,
):
    """Content identity is path-free; content and metric-time changes are fatal."""
    original_table = tmp_path / "original" / "trajectory.dat"
    relocated_table = tmp_path / "relocated" / "renamed-trajectory.dat"
    original_table.parent.mkdir()
    relocated_table.parent.mkdir()
    _write_table(original_table)
    shutil.copyfile(original_table, relocated_table)
    assert original_table.read_bytes() == relocated_table.read_bytes()

    initial_dir = tmp_path / "initial"
    initial = _run(
        _initial_command(initial_dir, original_table), "initial table run"
    )
    _require_success(initial, "initial table run")
    assert _load_balance_efficiency(initial, "initial table run") == 1.0
    fingerprints = re.findall(
        r"fingerprint=(content128-v1:[0-9a-f]{32}:[0-9]+)", initial.stdout
    )
    assert len(fingerprints) == 1, initial.stdout[-8000:]
    checkpoints = sorted((initial_dir / "rst").glob(f"{BASENAME}.*.rst"))
    assert checkpoints, "initial run did not write a checkpoint"
    checkpoint = checkpoints[-1]

    moved = _run(
        _restart_command(
            tmp_path / "same-content-new-path",
            checkpoint,
            "moved",
            (f"problem/trajectory_file={relocated_table.resolve()}",),
        ),
        "same-content relocated restart",
    )
    _require_success(moved, "same-content relocated restart")
    assert _load_balance_efficiency(moved, "same-content relocated restart") == 1.0
    assert f"file={relocated_table.resolve()}" in moved.stdout
    assert f"fingerprint={fingerprints[0]}" in moved.stdout

    # Mutate a numeric table entry after the checkpoint.  The altered table is still
    # well formed, subluminal, covered in time, and metric-safe, so only provenance
    # identity should reject this restart.
    _write_table(original_table, changed_middle_position=True)
    changed_content = _run(
        _restart_command(
            tmp_path / "changed-content", checkpoint, "changed_content"
        ),
        "changed-content restart",
    )
    _require_contract_failure(changed_content, "changed-content restart")

    changed_fd_step = _run(
        _restart_command(
            tmp_path / "changed-fd-step",
            checkpoint,
            "changed_fd_step",
            (
                f"problem/trajectory_file={relocated_table.resolve()}",
                "problem/metric_fd_step=1.0e-4",
            ),
        ),
        "changed metric_fd_step restart",
    )
    _require_contract_failure(changed_fd_step, "changed metric_fd_step restart")

    changed_time_offset = _run(
        _restart_command(
            tmp_path / "changed-time-offset",
            checkpoint,
            "changed_time_offset",
            (
                f"problem/trajectory_file={relocated_table.resolve()}",
                "problem/trajectory_time_offset=0.125",
            ),
        ),
        "changed trajectory_time_offset restart",
    )
    _require_contract_failure(
        changed_time_offset, "changed trajectory_time_offset restart"
    )
