"""Regression tests for scheduled output cadence across finalization and restart."""

from __future__ import annotations

import math
import re
import subprocess
from pathlib import Path

import test_suite.testutils as testutils


BASENAME = "output_finalize_cadence"
OUTPUT_DT = 0.1
TIME_TOLERANCE = 2.0e-14


def _write_input(
    path: Path,
    output_dt: float = OUTPUT_DT,
    *,
    restart_dcycle: int | None = None,
) -> None:
    restart_cadence = (
        f"dcycle = {restart_dcycle}"
        if restart_dcycle is not None
        else f"dt = {output_dt:.17g}"
    )
    path.write_text(
        f"""
<job>
basename = {BASENAME}

<mesh>
nghost = 2
nx1 = 8
x1min = 0.0
x1max = 1.0
ix1_bc = periodic
ox1_bc = periodic
nx2 = 1
x2min = 0.0
x2max = 1.0
ix2_bc = periodic
ox2_bc = periodic
nx3 = 1
x3min = 0.0
x3max = 1.0
ix3_bc = periodic
ox3_bc = periodic

<meshblock>
nx1 = 8
nx2 = 1
nx3 = 1

<mesh_refinement>
refinement = none

<time>
evolution = dynamic
integrator = rk2
cfl_number = 0.4
nlim = -1
tlim = 1.0
ndiag = 1

<hydro>
eos = ideal
reconstruct = plm
rsolver = llf
gamma = 1.6666666666666667

<problem>
pgen_name = linear_wave
wave_flag = 0
amp = 1.0e-6
dens = 1.0
pgas = 0.6
vx0 = 0.0
along_x1 = true
along_x2 = false
along_x3 = false

<output1>
file_type = bin
variable = hydro_w
dt = {output_dt:.17g}
ghost_zones = false

<output2>
file_type = rst
{restart_cadence}

<output3>
file_type = hst
data_format = %.17e
dt = {output_dt:.17g}
""".lstrip(),
        encoding="utf-8",
    )


def _run(
    run_dir: Path,
    input_path: Path | None,
    restart: Path | None,
    *,
    nlim: int,
    tlim: float,
) -> str:
    assert (input_path is None) != (restart is None)
    command = ["./athena"]
    if input_path is not None:
        command += ["-i", str(input_path.resolve())]
    else:
        command += ["-r", str(restart.resolve())]
    command += [
        "-d",
        str(run_dir.resolve()),
        f"time/nlim={nlim}",
        f"time/tlim={tlim:.17g}",
    ]
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=20,
    )
    with Path(testutils.LOG_FILE_PATH).open("a", encoding="utf-8") as log_file:
        log_file.write(f"\nOutput-cadence command: {' '.join(command)}\n")
        log_file.write(result.stdout)
    assert result.returncode == 0, result.stdout[-8000:]
    return result.stdout


def _numbered_paths(directory: Path, kind: str, suffix: str) -> list[Path]:
    if kind == "bin":
        pattern = f"{BASENAME}.{suffix}.*.bin"
    else:
        pattern = f"{BASENAME}.*.rst"
    return sorted((directory / kind).glob(pattern))


def _numbers(paths: list[Path]) -> list[int]:
    return [int(path.name.split(".")[-2]) for path in paths]


def _binary_time(path: Path) -> float:
    header = path.read_bytes()[:4096].decode("ascii", errors="ignore")
    match = re.search(r"(?m)^\s*time=(\S+)", header)
    assert match is not None, f"No binary time in {path}"
    return float(match.group(1))


def _history_times(directory: Path) -> list[float]:
    path = directory / f"{BASENAME}.hydro.hst"
    assert path.is_file()
    return [
        float(line.split()[0])
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]


def _terminal_cycle(stdout: str) -> int:
    matches = re.findall(r"(?m)^time=\S+\s+cycle=(\d+)\s*$", stdout)
    assert matches, stdout[-4000:]
    return int(matches[-1])


def _restart_parameter_text(path: Path) -> str:
    raw = path.read_bytes()
    marker = b"<par_end>\n"
    offset = raw.find(marker)
    assert offset >= 0, f"No ParameterInput terminator in {path}"
    return raw[: offset + len(marker)].decode("utf-8")


def _parameter(text: str, block: str, name: str) -> str:
    section = re.search(
        rf"(?ms)^<{re.escape(block)}>\s*$\n(.*?)(?=^<[^>]+>\s*$|\Z)", text
    )
    assert section is not None, f"No <{block}> block"
    values = re.findall(
        rf"(?m)^\s*{re.escape(name)}\s*=\s*(\S+)", section.group(1)
    )
    assert len(values) == 1, f"Unexpected {block}/{name}: {values}"
    return values[0]


def _assert_checkpoint_contract(
    path: Path,
    *,
    last_time: float,
    next_file_number: int,
    last_write_cycle: int | None = None,
) -> None:
    text = _restart_parameter_text(path)
    for block in ("output1", "output2", "output3"):
        actual = float(_parameter(text, block, "last_time"))
        assert math.isclose(actual, last_time, rel_tol=0.0, abs_tol=TIME_TOLERANCE), (
            f"{path} {block}/last_time={actual}, expected {last_time}"
        )
    for block in ("output1", "output2"):
        actual = int(_parameter(text, block, "file_number"))
        assert actual == next_file_number, (
            f"{path} {block}/file_number={actual}, expected {next_file_number}"
        )
    if last_write_cycle is not None:
        for block in ("output1", "output2", "output3"):
            actual = int(_parameter(text, block, "last_write_cycle"))
            assert actual == last_write_cycle, (
                f"{path} {block}/last_write_cycle={actual}, "
                f"expected {last_write_cycle}"
            )


def _assert_times(actual: list[float], expected: list[float]) -> None:
    assert len(actual) == len(expected), (actual, expected)
    for value, target in zip(actual, expected):
        assert math.isclose(
            value, target, rel_tol=0.0, abs_tol=TIME_TOLERANCE
        ), (actual, expected)


def test_early_finalize_preserves_scheduled_phase_across_restart(tmp_path):
    input_path = tmp_path / "output_cadence.athinput"
    _write_input(input_path)

    fresh = tmp_path / "fresh"
    _run(fresh, input_path, None, nlim=1, tlim=1.0)
    fresh_bins = _numbered_paths(fresh, "bin", "hydro_w")
    fresh_restarts = _numbered_paths(fresh, "rst", "rst")
    assert _numbers(fresh_bins) == [0, 1]
    assert _numbers(fresh_restarts) == [0, 1]
    fresh_times = [_binary_time(path) for path in fresh_bins]
    assert fresh_times[0] == 0.0
    assert 0.0 < fresh_times[1] < OUTPUT_DT
    _assert_times(_history_times(fresh), fresh_times)
    _assert_checkpoint_contract(
        fresh_restarts[-1], last_time=0.0, next_file_number=2
    )

    split = tmp_path / "split_restart"
    _run(split, None, fresh_restarts[-1], nlim=-1, tlim=0.15)
    split_bins = _numbered_paths(split, "bin", "hydro_w")
    split_restarts = _numbered_paths(split, "rst", "rst")
    assert _numbers(split_bins) == [2, 3]
    assert _numbers(split_restarts) == [2, 3]
    split_times = [_binary_time(path) for path in split_bins]
    _assert_times(_history_times(split), split_times)
    _assert_checkpoint_contract(
        split_restarts[-1], last_time=0.1, next_file_number=4
    )

    continuous = tmp_path / "continuous"
    _run(continuous, input_path, None, nlim=-1, tlim=0.15)
    continuous_bins = _numbered_paths(continuous, "bin", "hydro_w")
    continuous_restarts = _numbered_paths(continuous, "rst", "rst")
    assert _numbers(continuous_bins) == [0, 1, 2]
    assert _numbers(continuous_restarts) == [0, 1, 2]
    continuous_times = [_binary_time(path) for path in continuous_bins]
    _assert_times(_history_times(continuous), continuous_times)
    # The perturbed wave makes its CFL step slightly smaller than 0.05.  Compare the
    # actual first phase crossing: split and continuous must select the same endpoint.
    _assert_times(split_times, continuous_times[1:])
    _assert_times([continuous_times[0], continuous_times[-1]], [0.0, 0.15])
    _assert_checkpoint_contract(
        continuous_restarts[-1], last_time=0.1, next_file_number=3
    )


def test_finalize_at_due_tlim_consumes_exactly_one_phase(tmp_path):
    input_path = tmp_path / "output_cadence.athinput"
    _write_input(input_path)

    due = tmp_path / "due_at_tlim"
    _run(due, input_path, None, nlim=-1, tlim=0.1)
    due_bins = _numbered_paths(due, "bin", "hydro_w")
    due_restarts = _numbered_paths(due, "rst", "rst")
    assert _numbers(due_bins) == [0, 1]
    assert _numbers(due_restarts) == [0, 1]
    _assert_times([_binary_time(path) for path in due_bins], [0.0, 0.1])
    _assert_times(_history_times(due), [0.0, 0.1])
    _assert_checkpoint_contract(
        due_restarts[-1], last_time=0.1, next_file_number=2
    )

    resumed = tmp_path / "resumed_after_due_tlim"
    _run(resumed, None, due_restarts[-1], nlim=-1, tlim=0.15)
    resumed_bins = _numbered_paths(resumed, "bin", "hydro_w")
    resumed_restarts = _numbered_paths(resumed, "rst", "rst")
    # No overdue 0.1 output is repeated after restart.  The only new file is the forced
    # 0.15 endpoint, whose unique restored file number proves no prior file was replaced.
    assert _numbers(resumed_bins) == [2]
    assert _numbers(resumed_restarts) == [2]
    _assert_times([_binary_time(path) for path in resumed_bins], [0.15])
    _assert_times(_history_times(resumed), [0.15])
    _assert_checkpoint_contract(
        resumed_restarts[-1], last_time=0.1, next_file_number=3
    )


def test_finalize_does_not_double_advance_an_overdue_scheduled_output(tmp_path):
    input_path = tmp_path / "overdue_output_cadence.athinput"
    short_output_dt = 1.0e-4
    _write_input(input_path, short_output_dt)

    run_dir = tmp_path / "scheduled_then_finalize"
    _run(run_dir, input_path, None, nlim=1, tlim=1.0)
    binaries = _numbered_paths(run_dir, "bin", "hydro_w")
    restarts = _numbered_paths(run_dir, "rst", "rst")
    # Initialize, the cycle-1 scheduled write, and the duplicate final snapshot each keep
    # a unique file.  Only the scheduled write is allowed to advance the cadence.
    assert _numbers(binaries) == [0, 1, 2]
    assert _numbers(restarts) == [0, 1, 2]
    times = [_binary_time(path) for path in binaries]
    assert times[0] == 0.0
    assert times[1] == times[2] and times[1] > short_output_dt
    _assert_times(_history_times(run_dir), times)
    _assert_checkpoint_contract(
        restarts[-1],
        last_time=short_output_dt,
        next_file_number=3,
        last_write_cycle=1,
    )

    # Re-finalizing that checkpoint without taking another step must not consume another
    # overdue cadence phase.  The persisted write-cycle marker makes this idempotent even
    # across process boundaries.
    zero_step = tmp_path / "zero_step_restart"
    _run(zero_step, None, restarts[-1], nlim=1, tlim=1.0)
    zero_step_bins = _numbered_paths(zero_step, "bin", "hydro_w")
    zero_step_restarts = _numbered_paths(zero_step, "rst", "rst")
    assert _numbers(zero_step_bins) == [3]
    assert _numbers(zero_step_restarts) == [3]
    _assert_times([_binary_time(path) for path in zero_step_bins], [times[-1]])
    _assert_times(_history_times(zero_step), [times[-1]])
    _assert_checkpoint_contract(
        zero_step_restarts[-1],
        last_time=short_output_dt,
        next_file_number=4,
        last_write_cycle=1,
    )


def test_intermediate_restart_preserves_per_backend_write_cycle(tmp_path):
    input_path = tmp_path / "mixed_output_cadence.athinput"
    _write_input(input_path, restart_dcycle=1)

    exact = tmp_path / "exact_tlim"
    stdout = _run(exact, input_path, None, nlim=-1, tlim=OUTPUT_DT)
    cycle = _terminal_cycle(stdout)
    restarts = _numbered_paths(exact, "rst", "rst")
    assert len(restarts) >= 3

    # The penultimate checkpoint is emitted by dcycle at the exact tlim endpoint before
    # Finalize.  Time-cadence outputs are deliberately suppressed there, so its header
    # must retain different last-write cycles for the different backends.
    intermediate = restarts[-2]
    intermediate_text = _restart_parameter_text(intermediate)
    assert int(_parameter(intermediate_text, "output1", "last_write_cycle")) < cycle
    assert int(_parameter(intermediate_text, "output2", "last_write_cycle")) == cycle
    assert int(_parameter(intermediate_text, "output3", "last_write_cycle")) < cycle

    zero_step = tmp_path / "zero_step_intermediate_restart"
    _run(zero_step, None, intermediate, nlim=cycle, tlim=OUTPUT_DT)
    zero_step_bins = _numbered_paths(zero_step, "bin", "hydro_w")
    zero_step_restarts = _numbered_paths(zero_step, "rst", "rst")
    assert len(zero_step_bins) == 1
    assert len(zero_step_restarts) == 1
    _assert_times([_binary_time(zero_step_bins[0])], [OUTPUT_DT])
    _assert_times(_history_times(zero_step), [OUTPUT_DT])

    final_text = _restart_parameter_text(zero_step_restarts[0])
    assert math.isclose(
        float(_parameter(final_text, "output1", "last_time")),
        OUTPUT_DT,
        rel_tol=0.0,
        abs_tol=TIME_TOLERANCE,
    )
    assert float(_parameter(final_text, "output2", "last_time")) == 0.0
    assert math.isclose(
        float(_parameter(final_text, "output3", "last_time")),
        OUTPUT_DT,
        rel_tol=0.0,
        abs_tol=TIME_TOLERANCE,
    )
    for block in ("output1", "output2", "output3"):
        assert int(_parameter(final_text, block, "last_write_cycle")) == cycle
