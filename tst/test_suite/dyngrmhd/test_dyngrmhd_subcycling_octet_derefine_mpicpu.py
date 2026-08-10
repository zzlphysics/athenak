"""3D magnetic-AMR regression for strict 2:1 level subcycling."""

import re
import shutil
import struct
import subprocess
from pathlib import Path

import numpy as np
import pytest

import athena_read
import bin_convert
import test_suite.testutils as testutils


INPUT_FILE = "inputs/subcycling_octet_derefine_dyngrmhd.athinput"
BASENAME = "subcycling_octet_derefine_dyngrmhd"
STATE_VARIABLE = "mhd_u_bcc"
DIVB_VARIABLE = "mhd_divb"
MPI_TIMEOUT_SECONDS = 180
MAX_DIVB = 2.0e-13
EXPECTED_FULL_UPDATES = 61
EXPECTED_SPLIT_UPDATES = 23
EXPECTED_RESTART_UPDATES = 38
MAX_TEST_OUTPUT_BYTES = 64 * 1024 * 1024


def _run_mpi(
    run_dir: Path,
    ranks: int,
    restart=None,
    nlim=4,
    overrides=(),
    input_file=INPUT_FILE,
):
    """Run Athena under a process-tree timeout and return this run's log."""
    timeout_executable = shutil.which("timeout")
    assert timeout_executable is not None, "GNU timeout is required by this MPI test"

    command = [
        timeout_executable,
        "--signal=TERM",
        "--kill-after=10s",
        f"{MPI_TIMEOUT_SECONDS}s",
        "mpirun",
        "-np",
        str(ranks),
        "./athena",
    ]
    if restart is None:
        command += ["-i", str(Path(input_file).resolve())]
    else:
        command += ["-r", str(Path(restart).resolve())]
    command += [
        "-d",
        str(run_dir),
        f"job/basename={BASENAME}",
        f"time/nlim={nlim}",
    ]
    command += list(overrides)

    log_path = Path(testutils.LOG_FILE_PATH)
    log_offset = log_path.stat().st_size if log_path.exists() else 0
    success = testutils.run_command(command)
    with log_path.open("rb") as log_file:
        log_file.seek(log_offset)
        log_output = log_file.read().decode("utf-8", errors="replace")
    assert success, (
        f"{ranks}-rank Athena run failed or exceeded {MPI_TIMEOUT_SECONDS} s:\n"
        f"{log_output[-8000:]}"
    )
    return log_output


def _run_mpi_expect_failure(run_dir: Path, input_file: Path):
    """Require invalid root-step parameters to fail promptly and diagnostically."""
    timeout_executable = shutil.which("timeout")
    assert timeout_executable is not None, "GNU timeout is required by this MPI test"
    command = [
        timeout_executable,
        "--signal=TERM",
        "--kill-after=10s",
        "30s",
        "mpirun",
        "-np",
        "1",
        "./athena",
        "-i",
        str(input_file.resolve()),
        "-d",
        str(run_dir),
        f"job/basename={BASENAME}",
        "time/nlim=0",
    ]
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    assert result.returncode not in (0, 124, 137), (
        "Invalid <time>/root_dt_max did not fail promptly:\n"
        f"{result.stdout[-8000:]}"
    )
    assert re.search(
        r"<time>/root_dt_max must be finite and greater than zero", result.stdout
    ), result.stdout[-8000:]


def _input_with_root_dt_max(
    tmp_path: Path, name: str, value=None, subcycling="level"
):
    """Create a test-local input while keeping the checked-in default cap-free."""
    source = Path(INPUT_FILE).resolve()
    contents = source.read_text(encoding="utf-8")
    replacement = f"subcycling = {subcycling}"
    if value is not None:
        replacement += f"\nroot_dt_max = {value}"
    contents, count = re.subn(
        r"(?m)^subcycling\s*=\s*\S+\s*$", replacement, contents, count=1
    )
    assert count == 1, f"Unable to locate <time>/subcycling in {source}"
    output = tmp_path / f"{name}.athinput"
    output.write_text(contents, encoding="utf-8")
    return output


def _history_time_and_dt(run_dir: Path):
    history_path = run_dir / f"{BASENAME}.mhd.hst"
    assert history_path.is_file(), f"Missing history file {history_path}"
    history = athena_read.hst(str(history_path), raw=True)
    return np.asarray(history["time"]), np.asarray(history["dt"])


def _restart_parameter_text(path: Path):
    contents = path.read_bytes()
    marker = b"<par_end>\n"
    marker_offset = contents.find(marker)
    assert marker_offset >= 0, f"Missing ParameterInput terminator in {path}"
    return contents[:marker_offset].decode("utf-8")


def _restart_time_dt_cycle(path: Path):
    contents = path.read_bytes()
    marker = b"<par_end>\n"
    marker_offset = contents.find(marker)
    assert marker_offset >= 0, f"Missing ParameterInput terminator in {path}"
    offset = marker_offset + len(marker)
    mesh_prefix = struct.calcsize("=ii9d19i19i")
    return struct.unpack_from("=ddi", contents, offset + mesh_prefix)


def _without_restart_growth_parameter(source: Path, destination: Path):
    """Create a byte-valid legacy checkpoint lacking the new optional parameter."""
    contents = source.read_bytes()
    marker = b"<par_end>\n"
    marker_offset = contents.find(marker)
    assert marker_offset >= 0, f"Missing ParameterInput terminator in {source}"
    parameter_text = contents[:marker_offset].decode("utf-8")
    parameter_text, count = re.subn(
        r"(?m)^restart_dt_growth\s*=.*\n", "", parameter_text, count=1
    )
    assert count == 1, "Checkpoint did not contain restart_dt_growth"
    destination.write_bytes(
        parameter_text.encode("utf-8") + contents[marker_offset:]
    )


def _read_dumps(run_dir: Path, variable: str):
    paths = sorted((run_dir / "bin").glob(f"{BASENAME}.{variable}.*.bin"))
    assert paths, f"No {variable} dumps found in {run_dir}"
    return [(path, bin_convert.read_binary(str(path))) for path in paths]


def _latest_by_cycle(dumps):
    """Discard the duplicate endpoint dump emitted during finalization."""
    result = {}
    for path, dump in dumps:
        result[dump["cycle"]] = (path, dump)
    return {cycle: item[1] for cycle, item in result.items()}


def _canonical_order(dump):
    logical = np.asarray(dump["mb_logical"])
    return np.asarray(
        sorted(
            range(dump["n_mbs"]),
            key=lambda gid: (
                int(logical[gid, 3]),
                int(logical[gid, 2]),
                int(logical[gid, 1]),
                int(logical[gid, 0]),
            ),
        )
    )


def _canonical_field(dump, variable: str):
    return np.stack(dump["mb_data"][variable], axis=0)[_canonical_order(dump)]


def _relative_level_counts(dump):
    levels = np.asarray(dump["mb_logical"])[:, 3].astype(int)
    root_level = int(np.min(levels))
    unique, counts = np.unique(levels - root_level, return_counts=True)
    return dict(zip(unique.tolist(), counts.tolist()))


def _fine_locations(dump):
    logical = np.asarray(dump["mb_logical"])
    root_level = int(np.min(logical[:, 3]))
    return {
        tuple(location.tolist())
        for location in logical[logical[:, 3] == root_level + 1]
    }


def _expected_octet(x1_begin: int):
    return {
        (i, j, k, 1)
        for i in (x1_begin, x1_begin + 1)
        for j in (0, 1)
        for k in (0, 1)
    }


def _assert_scripted_topology(dumps_by_cycle):
    assert set(range(5)).issubset(dumps_by_cycle)
    assert _relative_level_counts(dumps_by_cycle[0]) == {0: 4}
    for cycle, x1_begin in ((1, 0), (2, 2), (3, 4)):
        assert _relative_level_counts(dumps_by_cycle[cycle]) == {0: 3, 1: 8}
        assert _fine_locations(dumps_by_cycle[cycle]) == _expected_octet(x1_begin)
    assert _relative_level_counts(dumps_by_cycle[4]) == {0: 4}

    # Each move must de-refine all eight old children, not retain an overlapping slab.
    for old_cycle, new_cycle in ((1, 2), (2, 3)):
        assert _fine_locations(dumps_by_cycle[old_cycle]).isdisjoint(
            _fine_locations(dumps_by_cycle[new_cycle])
        )


def _assert_same_topology(reference, candidate):
    for name in ("Nx1", "Nx2", "Nx3", "nx1_mb", "nx2_mb", "nx3_mb", "n_mbs"):
        assert candidate[name] == reference[name], f"Topology differs in {name}"

    reference_order = _canonical_order(reference)
    candidate_order = _canonical_order(candidate)
    np.testing.assert_array_equal(
        np.asarray(candidate["mb_logical"])[candidate_order],
        np.asarray(reference["mb_logical"])[reference_order],
    )
    np.testing.assert_array_equal(
        np.asarray(candidate["mb_index"])[candidate_order],
        np.asarray(reference["mb_index"])[reference_order],
    )
    np.testing.assert_array_equal(
        np.asarray(candidate["mb_geometry"])[candidate_order],
        np.asarray(reference["mb_geometry"])[reference_order],
    )


def _assert_state_close(reference, candidate, description):
    _assert_same_topology(reference, candidate)
    assert candidate["var_names"] == reference["var_names"]
    assert candidate["cycle"] == reference["cycle"]
    assert candidate["time"] == pytest.approx(reference["time"], rel=2.0e-14)
    for variable in reference["var_names"]:
        reference_values = _canonical_field(reference, variable)
        candidate_values = _canonical_field(candidate, variable)
        assert np.isfinite(reference_values).all()
        assert np.isfinite(candidate_values).all()
        np.testing.assert_allclose(
            candidate_values,
            reference_values,
            rtol=2.0e-13,
            atol=2.0e-14,
            err_msg=f"{description} changed {variable}",
        )


def _assert_divb_bounded(*run_dirs):
    maxima_by_run = []
    for run_dir in run_dirs:
        maxima_by_cycle = {}
        for path, dump in _read_dumps(run_dir, DIVB_VARIABLE):
            values = _canonical_field(dump, "divb")
            assert np.isfinite(values).all(), f"Non-finite div(B) in {path}"
            maximum = float(np.max(np.abs(values)))
            assert maximum < MAX_DIVB, (
                f"Magnetic constraint violation in {path}: "
                f"max|div(B)|={maximum:.16e}"
            )
            maxima_by_cycle[dump["cycle"]] = maximum
        maxima_by_run.append(maxima_by_cycle)
    return maxima_by_run


def _meshblock_updates(log_output):
    matches = re.findall(r"MeshBlock-cycles\s*=\s*(\d+)", log_output)
    assert matches, "Run log did not report MeshBlock-cycles"
    return int(matches[-1])


def _amr_counts(log_output):
    matches = re.findall(r"(\d+) MeshBlocks created, (\d+) deleted by AMR", log_output)
    assert matches, "Run log did not report AMR creation/deletion counts"
    return tuple(map(int, matches[-1]))


def _communicated_blocks(log_output):
    matches = re.findall(r"(\d+) communicated for load balancing", log_output)
    assert matches, "Run log did not report load-balancing communication"
    return int(matches[-1])


def test_subcycling_root_dt_max_is_opt_in_validated_and_restartable(tmp_path):
    cap = 1.0e-4

    for index, invalid in enumerate(("0", "-1", "nan", "inf")):
        invalid_input = _input_with_root_dt_max(
            tmp_path, f"invalid_{index}", invalid
        )
        _run_mpi_expect_failure(tmp_path / f"invalid_run_{index}", invalid_input)

    # Omitting the parameter retains the historical level-subcycling CFL step.
    default_dir = tmp_path / "level_default"
    _run_mpi(default_dir, ranks=1, nlim=1)
    default_times, default_dt = _history_time_and_dt(default_dir)
    assert default_times[0] == pytest.approx(0.0, abs=0.0)
    assert default_dt[0] > 10.0 * cap

    # The cap applies to both the initial root step and the subsequent growth phase.
    cap_input = _input_with_root_dt_max(tmp_path, "capped", f"{cap:.17e}")
    continuous_dir = tmp_path / "capped_continuous"
    _run_mpi(
        continuous_dir,
        ranks=1,
        nlim=2,
        input_file=cap_input,
        overrides=("output4/dcycle=1",),
    )
    continuous_times, continuous_dt = _history_time_and_dt(continuous_dir)
    np.testing.assert_allclose(
        np.unique(continuous_times),
        np.asarray((0.0, cap, 2.0 * cap)),
        rtol=2.0e-14,
        atol=0.0,
    )
    np.testing.assert_allclose(continuous_dt, cap, rtol=2.0e-14, atol=0.0)

    # An explicit value is accepted but has no timestep effect in legacy mode.
    legacy_default_input = _input_with_root_dt_max(
        tmp_path, "legacy_default", subcycling="none"
    )
    legacy_capped_input = _input_with_root_dt_max(
        tmp_path, "legacy_with_ignored_cap", f"{cap:.17e}", subcycling="none"
    )
    legacy_default_dir = tmp_path / "legacy_default"
    legacy_capped_dir = tmp_path / "legacy_with_ignored_cap"
    _run_mpi(
        legacy_default_dir, ranks=1, nlim=2, input_file=legacy_default_input
    )
    _run_mpi(
        legacy_capped_dir, ranks=1, nlim=2, input_file=legacy_capped_input
    )
    legacy_default_time, legacy_default_dt = _history_time_and_dt(legacy_default_dir)
    legacy_capped_time, legacy_capped_dt = _history_time_and_dt(legacy_capped_dir)
    np.testing.assert_array_equal(legacy_capped_time, legacy_default_time)
    np.testing.assert_array_equal(legacy_capped_dt, legacy_default_dt)
    assert legacy_capped_dt[0] > 10.0 * cap

    # ParameterInput is serialized in the checkpoint and must re-arm the same cap
    # without a command-line override after restart.
    split_dir = tmp_path / "capped_split"
    _run_mpi(
        split_dir,
        ranks=1,
        nlim=1,
        input_file=cap_input,
        overrides=("output4/dcycle=1",),
    )
    split_restarts = sorted((split_dir / "rst").glob(f"{BASENAME}.*.rst"))
    assert split_restarts, "Expected capped split checkpoint"
    split_restart = split_restarts[-1]
    parameter_text = _restart_parameter_text(split_restart)
    assert re.search(
        rf"(?m)^\s*root_dt_max\s*=\s*{re.escape(f'{cap:.17e}')}\s*$",
        parameter_text,
    )

    restart_dir = tmp_path / "capped_restart"
    _run_mpi(
        restart_dir,
        ranks=1,
        restart=split_restart,
        nlim=2,
        overrides=("output4/dcycle=1",),
    )
    restart_times, restart_dt = _history_time_and_dt(restart_dir)
    np.testing.assert_allclose(restart_times, 2.0 * cap, rtol=2.0e-14, atol=0.0)
    np.testing.assert_allclose(restart_dt, cap, rtol=2.0e-14, atol=0.0)
    restart_restarts = sorted(
        (restart_dir / "rst").glob(f"{BASENAME}.*.rst")
    )
    assert restart_restarts, "Expected checkpoint after capped restart"
    assert "root_dt_max" in _restart_parameter_text(restart_restarts[-1])


def test_tlim_roundoff_tail_does_not_throttle_restart(tmp_path):
    cap = 1.0e-6
    endpoint = np.nextafter(cap, np.inf)
    tail = endpoint - cap
    assert 0.0 < tail < 64.0 * np.finfo(float).eps

    cap_input = _input_with_root_dt_max(
        tmp_path, "tail_capped", f"{cap:.17e}"
    )
    # linear_wave interprets tlim in wave periods and rescales it during initial-data
    # setup. Measure that deterministic factor so the evolved endpoint is exactly the
    # adjacent representable value above one root cap.
    probe_dir = tmp_path / "tail_scale_probe"
    _run_mpi(
        probe_dir,
        ranks=1,
        nlim=0,
        input_file=cap_input,
        overrides=("time/tlim=1.0", "output4/dcycle=1"),
    )
    probe_restarts = sorted(
        (probe_dir / "rst").glob(f"{BASENAME}.*.rst")
    )
    assert probe_restarts, "Expected a cycle-zero scale-probe checkpoint"
    probe_text = _restart_parameter_text(probe_restarts[-1])
    scale_match = re.search(r"(?m)^tlim\s*=\s*(\S+)", probe_text)
    assert scale_match, "Scale-probe checkpoint did not contain tlim"
    tlim_scale = float(scale_match.group(1))
    assert tlim_scale > 0.0

    tail_dir = tmp_path / "tail_endpoint"
    _run_mpi(
        tail_dir,
        ranks=1,
        nlim=-1,
        input_file=cap_input,
        overrides=(
            f"time/tlim={endpoint/tlim_scale:.17e}",
            "output4/dcycle=1",
        ),
    )
    tail_restarts = sorted(
        (tail_dir / "rst").glob(f"{BASENAME}.*.rst")
    )
    assert tail_restarts, "Expected a checkpoint after the roundoff tail"
    checkpoint = tail_restarts[-1]
    checkpoint_time, checkpoint_dt, checkpoint_cycle = _restart_time_dt_cycle(
        checkpoint
    )
    assert checkpoint_cycle == 2
    assert checkpoint_time == pytest.approx(endpoint, rel=0.0, abs=0.0)
    assert checkpoint_dt == pytest.approx(tail, rel=0.0, abs=0.0)

    parameter_text = _restart_parameter_text(checkpoint)
    match = re.search(r"(?m)^restart_dt_growth\s*=\s*(\S+)", parameter_text)
    assert match, "Checkpoint did not serialize restart_dt_growth"
    assert float(match.group(1)) == pytest.approx(cap, rel=2.0e-14, abs=0.0)

    # A new checkpoint resumes at the normal root cap, not 2*tail.
    resumed_dir = tmp_path / "tail_resumed"
    _run_mpi(
        resumed_dir,
        ranks=1,
        restart=checkpoint,
        nlim=3,
        overrides=("time/tlim=1.0", "output4/dcycle=1"),
    )
    resumed_time, resumed_dt = _history_time_and_dt(resumed_dir)
    np.testing.assert_allclose(resumed_time, endpoint + cap, rtol=2.0e-14, atol=0.0)
    np.testing.assert_allclose(resumed_dt, cap, rtol=2.0e-14, atol=0.0)

    # Existing checkpoints have no optional parameter.  A roundoff-scale legacy header
    # must remove only the growth constraint; the freshly evaluated CFL/root cap remains.
    legacy_checkpoint = tmp_path / "tail_legacy.rst"
    _without_restart_growth_parameter(checkpoint, legacy_checkpoint)
    legacy_dir = tmp_path / "tail_legacy_resumed"
    _run_mpi(
        legacy_dir,
        ranks=1,
        restart=legacy_checkpoint,
        nlim=3,
        overrides=("time/tlim=1.0", "output4/dcycle=1"),
    )
    legacy_time, legacy_dt = _history_time_and_dt(legacy_dir)
    np.testing.assert_allclose(legacy_time, endpoint + cap, rtol=2.0e-14, atol=0.0)
    np.testing.assert_allclose(legacy_dt, cap, rtol=2.0e-14, atol=0.0)


def test_subcycling_preserves_divb_on_anisotropic_2d_cells(tmp_path):
    """Cover the separate 2D Toth-Roe formula with dx1:dx2=1:2."""
    run_dir = tmp_path / "anisotropic_2d"
    log_output = _run_mpi(
        run_dir,
        ranks=1,
        nlim=2,
        overrides=("mesh/nx3=1", "meshblock/nx3=1"),
    )
    assert _meshblock_updates(log_output) == 15
    assert _amr_counts(log_output) == (6, 3)
    maxima = _assert_divb_bounded(run_dir)[0]
    assert maxima[1] < MAX_DIVB
    assert maxima[2] < MAX_DIVB


def test_subcycling_preserves_divb_through_3d_octet_amr_and_restart(tmp_path):
    one_dir = tmp_path / "one_rank"
    three_dir = tmp_path / "three_ranks"
    split_dir = tmp_path / "split_three_ranks"
    restart_dir = tmp_path / "restart_two_ranks"

    one_log = _run_mpi(one_dir, ranks=1)
    three_log = _run_mpi(three_dir, ranks=3)
    split_log = _run_mpi(
        split_dir,
        ranks=3,
        nlim=2,
        overrides=("output4/dcycle=2",),
    )
    split_restarts = sorted((split_dir / "rst").glob(f"{BASENAME}.*.rst"))
    assert split_restarts, "Expected a cycle-2 split checkpoint"
    restart_log = _run_mpi(
        restart_dir,
        ranks=2,
        restart=split_restarts[-1],
        nlim=4,
        overrides=("output4/dcycle=0",),
    )

    assert _meshblock_updates(one_log) == EXPECTED_FULL_UPDATES
    assert _meshblock_updates(three_log) == EXPECTED_FULL_UPDATES
    assert _meshblock_updates(split_log) == EXPECTED_SPLIT_UPDATES
    assert _meshblock_updates(restart_log) == EXPECTED_RESTART_UPDATES
    assert _amr_counts(one_log) == (21, 21)
    assert _amr_counts(three_log) == (21, 21)
    assert _amr_counts(split_log) == (14, 7)
    assert _amr_counts(restart_log) == (7, 14)
    assert _communicated_blocks(three_log) > 0
    assert _communicated_blocks(split_log) > 0
    assert _communicated_blocks(restart_log) > 0

    divb_maxima = _assert_divb_bounded(one_dir, three_dir, split_dir, restart_dir)
    # Cycle 1 catches missing face-area scaling on anisotropic cells; cycle 2 catches
    # stale coarse_b0 scratch consumed by complete-octet de-refinement.
    for maxima in divb_maxima[:3]:
        assert maxima[1] < MAX_DIVB
        assert maxima[2] < MAX_DIVB

    one_states = _latest_by_cycle(_read_dumps(one_dir, STATE_VARIABLE))
    three_states = _latest_by_cycle(_read_dumps(three_dir, STATE_VARIABLE))
    split_states = _latest_by_cycle(_read_dumps(split_dir, STATE_VARIABLE))
    restart_states = _latest_by_cycle(_read_dumps(restart_dir, STATE_VARIABLE))
    _assert_scripted_topology(one_states)
    _assert_scripted_topology(three_states)
    _assert_scripted_topology({**split_states, **restart_states})

    _assert_state_close(one_states[4], three_states[4], "MPI decomposition")
    _assert_state_close(three_states[2], split_states[2], "early termination")
    _assert_state_close(
        three_states[4], restart_states[4], "three-to-two-rank restart"
    )

    total_bytes = sum(
        path.stat().st_size for path in tmp_path.rglob("*") if path.is_file()
    )
    assert total_bytes < MAX_TEST_OUTPUT_BYTES, (
        f"Regression generated {total_bytes / 1024**2:.1f} MiB, "
        f"exceeding its storage budget"
    )
