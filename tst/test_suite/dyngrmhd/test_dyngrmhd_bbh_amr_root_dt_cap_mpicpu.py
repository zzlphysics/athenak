"""Root-step-cap regression for moving-BBH AMR lookahead under level subcycling."""

from __future__ import annotations

import math
import re
import shutil
import subprocess
from pathlib import Path

import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
INPUT_FILE = ROOT / "inputs" / "dyngr" / "effective_bbh_amr_subcycle_smoke.athinput"
ROOT_DT_MAX = 1.0e-3
TIMEOUT_SECONDS = 60


def _write_capped_input(path: Path) -> None:
    text = INPUT_FILE.read_text(encoding="utf-8")
    anchor = "subcycling = level\n"
    assert text.count(anchor) == 1
    path.write_text(
        text.replace(anchor, f"{anchor}root_dt_max = {ROOT_DT_MAX:.17g}\n"),
        encoding="utf-8",
    )


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


def _run(input_path: Path, run_dir: Path, *, ranks: int) -> str:
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
        "time/tlim=1.0e-3",
        "output1/dcycle=0",
        "output2/dcycle=0",
        "output3/dcycle=0",
        "mesh/x1min=-12.0",
        "mesh/x1max=12.0",
        "mesh/x2min=-4.9985",
        "mesh/x2max=3.0015",
        "mesh/x3min=-4.9985",
        "mesh/x3max=3.0015",
        "mesh_refinement/num_levels=2",
        "problem/omega=0.0",
        "problem/refinement_radius=1.0",
        "problem/refinement_radius_level_1=1.0",
        "problem/refinement_horizon_factor=1.0",
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

    # Both runs take the same physical root step to tlim.  At that endpoint the uncapped
    # factor-two growth bound gives a 0.002M padding, whereas the reachable capped step
    # gives 0.001M.  The deliberately offset root-block faces make this distinguish four
    # parent blocks without relying on roundoff or evolved state.
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
