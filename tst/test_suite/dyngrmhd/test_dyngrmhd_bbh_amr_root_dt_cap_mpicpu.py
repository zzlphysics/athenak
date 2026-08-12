"""Root-step-cap regression for moving-BBH AMR lookahead under level subcycling."""

from __future__ import annotations

import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
import read_athenak_restart_metadata as restart_metadata  # noqa: E402


INPUT_FILE = ROOT / "inputs" / "dyngr" / "effective_bbh_amr_subcycle_smoke.athinput"
ROOT_DT_MAX = 1.0e-3
TRANSVERSE_MIN = -4.9996
TRANSVERSE_MAX = 3.0004
TIMEOUT_SECONDS = 60


def _write_capped_input(path: Path) -> None:
    text = INPUT_FILE.read_text(encoding="utf-8")
    anchor = "subcycling = level\n"
    assert text.count(anchor) == 1
    path.write_text(
        text.replace(anchor, f"{anchor}root_dt_max = {ROOT_DT_MAX:.17g}\n"),
        encoding="utf-8",
    )


def _write_stationary_table(path: Path) -> None:
    """Write immutable coverage well beyond either segment boundary."""
    state = (
        10.0,
        0.0,
        0.0,
        -10.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.5,
        0.5,
    )
    rows = [
        " ".join(f"{value:.17e}" for value in (time, *state))
        for time in (-0.5, 0.0, 0.5)
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _mpi_launcher() -> str:
    launcher = shutil.which("mpirun") or shutil.which("mpiexec")
    if launcher is not None:
        return launcher
    cache = Path("..").resolve() / "CMakeCache.txt"
    if cache.is_file():
        for line in cache.read_text(encoding="utf-8").splitlines():
            if line.startswith("MPI_CXX_COMPILER:FILEPATH="):
                compiler = Path(line.split("=", 1)[1])
                for name in ("mpirun", "mpiexec"):
                    candidate = compiler.with_name(name)
                    if candidate.is_file():
                        return str(candidate)
    raise AssertionError("an MPI launcher matching the test binary is required")


def _run(
    input_path: Path,
    run_dir: Path,
    *,
    ranks: int,
    tlim: float = ROOT_DT_MAX,
    overrides: tuple[str, ...] = (),
) -> str:
    timeout = shutil.which("timeout")
    assert timeout is not None, "GNU timeout is required by this MPI regression"
    command = [
        timeout,
        "--signal=TERM",
        "--kill-after=5s",
        f"{TIMEOUT_SECONDS}s",
        _mpi_launcher(),
        "-np",
        str(ranks),
        "./athena",
        "-i",
        str(input_path),
        "-d",
        str(run_dir),
        "time/nlim=1",
        f"time/tlim={tlim:.17g}",
        "output1/dcycle=0",
        "output2/dcycle=0",
        "output3/dcycle=0",
        "mesh/x1min=-12.0",
        "mesh/x1max=12.0",
        f"mesh/x2min={TRANSVERSE_MIN:.17g}",
        f"mesh/x2max={TRANSVERSE_MAX:.17g}",
        f"mesh/x3min={TRANSVERSE_MIN:.17g}",
        f"mesh/x3max={TRANSVERSE_MAX:.17g}",
        "mesh_refinement/num_levels=2",
        "problem/omega=0.0",
        "problem/refinement_radius=1.0",
        "problem/refinement_radius_level_1=1.0",
        "problem/refinement_horizon_factor=1.0",
        *overrides,
    ]
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=TIMEOUT_SECONDS + 10,
    )
    with Path(testutils.LOG_FILE_PATH).open("a", encoding="utf-8") as log_file:
        log_file.write(f"\nBBH AMR root-cap command: {' '.join(command)}\n")
        log_file.write(result.stdout)
    assert result.returncode == 0, result.stdout[-8000:]
    return result.stdout


def _run_restart(restart: Path, run_dir: Path) -> str:
    timeout = shutil.which("timeout")
    assert timeout is not None, "GNU timeout is required by this MPI regression"
    command = [
        timeout,
        "--signal=TERM",
        "--kill-after=5s",
        f"{TIMEOUT_SECONDS}s",
        _mpi_launcher(),
        "-np",
        "1",
        "./athena",
        "-r",
        str(restart),
        "-d",
        str(run_dir),
        "time/nlim=2",
        "time/tlim=0.1",
    ]
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=TIMEOUT_SECONDS + 10,
    )
    with Path(testutils.LOG_FILE_PATH).open("a", encoding="utf-8") as log_file:
        log_file.write(f"\nBBH AMR clipped-restart command: {' '.join(command)}\n")
        log_file.write(result.stdout)
    assert result.returncode == 0, result.stdout[-8000:]
    return result.stdout


def _initial_dt(log: str) -> float:
    match = re.search(
        r"(?m)^elapsed=\S+ cycle=0 time=0\.000000e\+00 dt=(\S+)\s*$", log
    )
    assert match is not None, log[-4000:]
    return float(match.group(1))


def _terminal_contract(log: str) -> tuple[float, int, int, int]:
    terminal = re.findall(r"(?m)^time=(\S+) cycle=(\d+)\s*$", log)
    blocks = re.findall(r"Current number of MeshBlocks\s*=\s*(\d+)", log)
    amr = re.findall(r"(\d+) MeshBlocks created, (\d+) deleted by AMR", log)
    assert terminal and blocks and amr, log[-4000:]
    time, cycle = terminal[-1]
    created, deleted = amr[-1]
    assert int(cycle) == 1
    assert int(deleted) == 0
    return float(time), int(blocks[-1]), int(created), int(deleted)


def test_track_amr_lookahead_honors_enforced_root_dt_max(tmp_path):
    """The AMR halo must not target a root step the Driver cannot take."""
    capped_input = tmp_path / "effective_bbh_amr_cap_smoke.athinput"
    _write_capped_input(capped_input)

    uncapped_log = _run(INPUT_FILE, tmp_path / "uncapped", ranks=1)
    capped_log = _run(capped_input, tmp_path / "capped", ranks=1)
    capped_mpi_log = _run(capped_input, tmp_path / "capped_mpi", ranks=2)

    uncapped_dt = _initial_dt(uncapped_log)
    capped_dt = _initial_dt(capped_log)
    assert math.isclose(capped_dt, uncapped_dt, rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(uncapped_dt, ROOT_DT_MAX, rel_tol=0.0, abs_tol=1.0e-12)

    uncapped_time, uncapped_nmb, uncapped_created, _ = _terminal_contract(
        uncapped_log
    )
    capped_time, capped_nmb, capped_created, _ = _terminal_contract(capped_log)
    assert math.isclose(capped_time, uncapped_time, rel_tol=0.0, abs_tol=1.0e-12)

    # Both runs take the same physical root step.  Looking one reachable step through
    # the checkpoint, the uncapped factor-two growth bound gives a 0.0005M half-sample
    # padding, whereas the capped bound gives 0.00025M.  Deliberately offset root-block
    # faces make this distinguish four parent blocks without relying on roundoff or
    # evolved state.
    assert (uncapped_nmb, uncapped_created) == (134, 70)
    assert (capped_nmb, capped_created) == (106, 42)
    capped_mpi_time, capped_mpi_nmb, capped_mpi_created, _ = _terminal_contract(
        capped_mpi_log
    )
    assert math.isclose(
        capped_mpi_time, capped_time, rel_tol=0.0, abs_tol=1.0e-12
    )
    assert (capped_mpi_nmb, capped_mpi_created) == (106, 42)
    assert "root-step lookahead cap=" not in uncapped_log
    assert "root-step lookahead cap=0.001" in capped_log
    assert "root-step lookahead cap=0.001" in capped_mpi_log
    assert "Number of parallel ranks = 2" in capped_mpi_log


def test_track_amr_checkpoint_tlim_does_not_change_future_topology(tmp_path):
    """A segment boundary must not masquerade as the end of the trajectory."""
    capped_input = tmp_path / "effective_bbh_amr_cap_smoke.athinput"
    _write_capped_input(capped_input)
    trajectory = tmp_path / "stationary_trajectory.dat"
    _write_stationary_table(trajectory)
    modes = {
        "circular": (),
        "table": (
            "problem/trajectory_mode=table",
            f"problem/trajectory_file={trajectory.resolve()}",
            "problem/trajectory_time_offset=0.0",
        ),
    }

    for label, overrides in modes.items():
        endpoint_log = _run(
            capped_input,
            tmp_path / f"{label}_checkpoint_endpoint",
            ranks=1,
            tlim=ROOT_DT_MAX,
            overrides=overrides,
        )
        continuing_log = _run(
            capped_input,
            tmp_path / f"{label}_continuing_segment",
            ranks=1,
            tlim=0.1,
            overrides=overrides,
        )
        endpoint = _terminal_contract(endpoint_log)
        continuing = _terminal_contract(continuing_log)

        # Both executions complete exactly one capped root step from identical data.
        # The first stops at tlim while the second can continue.  Their post-AMR
        # topology must be identical so the checkpoint is already safe for the next
        # segment.  The former implementation added the entire unsampled lookahead at
        # tlim and created 28 extra children here.
        assert math.isclose(
            endpoint[0], continuing[0], rel_tol=0.0, abs_tol=1.0e-12
        )
        assert endpoint[1:] == continuing[1:] == (106, 42, 0)


def test_short_tlim_clip_preserves_restart_growth_lookahead(tmp_path):
    """A clipped endpoint must already cover the full step restored on restart."""
    capped_input = tmp_path / "effective_bbh_amr_cap_smoke.athinput"
    _write_capped_input(capped_input)
    stationary = ("problem/omega=0.0",)

    clipped_log = _run(
        capped_input,
        tmp_path / "short_clipped_endpoint",
        ranks=1,
        tlim=ROOT_DT_MAX / 4.0,
        overrides=stationary,
    )
    continuing_log = _run(
        capped_input,
        tmp_path / "full_reachable_step",
        ranks=1,
        tlim=0.1,
        overrides=stationary,
    )
    clipped = _terminal_contract(clipped_log)
    continuing = _terminal_contract(continuing_log)

    assert math.isclose(
        clipped[0], ROOT_DT_MAX / 4.0, rel_tol=0.0, abs_tol=1.0e-12
    )
    assert math.isclose(
        continuing[0], ROOT_DT_MAX, rel_tol=0.0, abs_tol=1.0e-12
    )
    # The targets are stationary.  Although tlim clips the first evolution step, the
    # checkpoint retains ROOT_DT_MAX as its restart growth reference.  Its post-AMR
    # hierarchy must therefore cover the same reachable next step as the continuing
    # execution.  The old pmesh->dt-based lookahead covered only half this interval.
    assert clipped[1:] == continuing[1:] == (106, 42, 0)


def test_short_tlim_checkpoint_restores_full_step_and_tree(tmp_path):
    """Exercise serialization, BuildTreeFromRestart, and restored growth together."""
    capped_input = tmp_path / "effective_bbh_amr_cap_smoke.athinput"
    _write_capped_input(capped_input)
    clipped_dir = tmp_path / "clipped_checkpoint"
    _run(
        capped_input,
        clipped_dir,
        ranks=1,
        tlim=ROOT_DT_MAX / 4.0,
        overrides=("problem/omega=0.0", "output3/dcycle=1"),
    )
    clipped_restarts = sorted((clipped_dir / "rst").glob("*.rst"))
    assert [path.name.split(".")[-2] for path in clipped_restarts] == [
        "00000",
        "00001",
    ]
    checkpoint = clipped_restarts[-1]
    source = restart_metadata.read_restart_metadata(checkpoint)
    assert source.cycle == 1
    assert math.isclose(
        source.time, ROOT_DT_MAX / 4.0, rel_tol=0.0, abs_tol=1.0e-12
    )
    assert math.isclose(
        source.last_dt, ROOT_DT_MAX / 4.0, rel_tol=0.0, abs_tol=1.0e-12
    )
    assert math.isclose(
        float(source.parameters["time"]["restart_dt_growth"]),
        ROOT_DT_MAX,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    assert source.nmb_total == 106
    assert source.amr_cycle_counters is not None

    resumed_dir = tmp_path / "resumed_full_step"
    restart_log = _run_restart(checkpoint, resumed_dir)
    first = re.search(
        r"(?m)^elapsed=\S+ cycle=(\d+) time=(\S+) dt=(\S+)\s*$",
        restart_log,
    )
    assert first is not None, restart_log[-4000:]
    assert int(first.group(1)) == 1
    assert math.isclose(
        float(first.group(2)), ROOT_DT_MAX / 4.0, rel_tol=0.0, abs_tol=1.0e-12
    )
    assert math.isclose(
        float(first.group(3)), ROOT_DT_MAX, rel_tol=0.0, abs_tol=1.0e-12
    )

    resumed_restarts = sorted((resumed_dir / "rst").glob("*.rst"))
    assert [path.name.split(".")[-2] for path in resumed_restarts] == ["00002"]
    resumed = restart_metadata.read_restart_metadata(resumed_restarts[0])
    assert resumed.cycle == 2
    assert math.isclose(
        resumed.time,
        ROOT_DT_MAX / 4.0 + ROOT_DT_MAX,
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    assert math.isclose(
        resumed.last_dt, ROOT_DT_MAX, rel_tol=0.0, abs_tol=1.0e-12
    )
    assert source.locations == resumed.locations
    assert source.costs == resumed.costs
    assert resumed.amr_cycle_counters is not None
    assert all(
        resumed_counter == source_counter + 1
        for source_counter, resumed_counter in zip(
            source.amr_cycle_counters, resumed.amr_cycle_counters
        )
    )
