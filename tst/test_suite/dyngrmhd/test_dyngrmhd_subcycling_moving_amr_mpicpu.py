"""Moving-AMR/restart regression for strict 2:1 level subcycling in dynGRMHD."""

import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pytest

import athena_read
import bin_convert
import test_suite.testutils as testutils


INPUT_FILE = "inputs/subcycling_moving_amr_dyngrmhd.athinput"
BASENAME = "subcycling_moving_amr_dyngrmhd"
STATE_VARIABLE = "mhd_u_bcc"
DIVB_VARIABLE = "mhd_divb"
MPI_TIMEOUT_SECONDS = 180
NEGATIVE_MPI_TIMEOUT_SECONDS = 30
EXPECTED_FULL_UPDATES = 104
EXPECTED_HALF_UPDATES = 52
# This bound is deliberately tight enough to catch stale coarse_b0 face data being
# copied into a newly de-refined parent.  The broken path produces 4.7--4.9e-13 in this
# problem, while a clean face-centered restriction stays below 9e-14 on both one and
# three MPI ranks.
MAX_DIVB = 2.0e-13


def _timed_mpi_command(
    run_dir: Path,
    ranks: int,
    restart=None,
    nlim=6,
    overrides=(),
    timeout_seconds=MPI_TIMEOUT_SECONDS,
):
    timeout_executable = shutil.which("timeout")
    assert timeout_executable is not None, "GNU timeout is required by this MPI test"

    athena_arguments = [
        "mpirun",
        "-np",
        str(ranks),
        "./athena",
    ]
    if restart is None:
        athena_arguments += ["-i", INPUT_FILE]
    else:
        athena_arguments += ["-r", str(Path(restart).resolve())]
        # Restart extensions must be detected from the file, not this mutable key.
        athena_arguments += ["mesh_refinement/restart_cycle_counters_version=0"]
    athena_arguments += [
        "-d",
        str(run_dir),
        f"job/basename={BASENAME}",
        f"time/nlim={nlim}",
    ]
    athena_arguments += list(overrides)
    command = [
        timeout_executable,
        "--signal=TERM",
        "--kill-after=10s",
        f"{timeout_seconds}s",
    ] + athena_arguments
    return command


def _run_timed_mpi(run_dir: Path, ranks: int, restart=None, nlim=6, overrides=()):
    """Run Athena under a process-tree timeout and return this run's log text."""
    log_path = Path(testutils.LOG_FILE_PATH)
    log_offset = log_path.stat().st_size if log_path.exists() else 0
    command = _timed_mpi_command(run_dir, ranks, restart, nlim, overrides)
    success = testutils.run_command(command)

    with log_path.open("rb") as log_file:
        log_file.seek(log_offset)
        log_output = log_file.read().decode("utf-8", errors="replace")
    assert success, (
        f"{ranks}-rank Athena run failed or exceeded {MPI_TIMEOUT_SECONDS} s:\n"
        f"{log_output[-8000:]}"
    )
    return log_output


def _run_timed_mpi_expect_failure(
    run_dir: Path, ranks: int, restart, expected_patterns, overrides=()
):
    """Require a prompt Athena failure, rather than a GNU-timeout termination."""
    command = _timed_mpi_command(
        run_dir,
        ranks,
        restart,
        overrides=overrides,
        timeout_seconds=NEGATIVE_MPI_TIMEOUT_SECONDS,
    )
    start = time.monotonic()
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    elapsed = time.monotonic() - start
    output = result.stdout

    with Path(testutils.LOG_FILE_PATH).open("a", encoding="utf-8") as log_file:
        log_file.write(f"\nExpected-failure command: {' '.join(command)}\n")
        log_file.write(output)

    assert result.returncode != 0, (
        "Invalid configuration or restart unexpectedly succeeded"
    )
    assert result.returncode not in (124, 137), (
        f"Invalid configuration or restart was killed by GNU timeout after "
        f"{elapsed:.1f} s:\n"
        f"{output[-8000:]}"
    )
    assert elapsed < NEGATIVE_MPI_TIMEOUT_SECONDS, (
        f"Invalid configuration or restart did not fail within "
        f"{NEGATIVE_MPI_TIMEOUT_SECONDS} s:\n{output[-8000:]}"
    )
    for pattern in expected_patterns:
        assert re.search(pattern, output, flags=re.IGNORECASE), (
            f"Invalid configuration or restart did not report /{pattern}/:\n"
            f"{output[-8000:]}"
        )
    return output


def test_subcycling_moving_amr_honors_block_capacity(tmp_path):
    """Repeated AMR repartitioning must stay within fixed MeshBlockPack capacity."""
    run_dir = tmp_path / "capacity_constrained"
    run_dir.mkdir()
    log_output = _run_timed_mpi(
        run_dir,
        ranks=3,
        overrides=("mesh_refinement/max_nmb_per_rank=5",),
    )

    assert "Current number of MeshBlocks = 14" in log_output
    assert "30 MeshBlocks created, 24 deleted by AMR" in log_output
    assert "exceeds <mesh_refinement>/max_nmb_per_rank" not in log_output


def _read_dumps(run_dir: Path, variable: str):
    pattern = f"{BASENAME}.{variable}.*.bin"
    paths = sorted((run_dir / "bin").glob(pattern))
    assert paths, f"No binary dumps matched {pattern}"
    return [(path, bin_convert.read_binary(str(path))) for path in paths]


def _restart_binary_payload_size(path: Path):
    contents = path.read_bytes()
    marker = b"<par_end>\n"
    marker_offset = contents.find(marker)
    assert marker_offset >= 0, f"Missing ParameterInput terminator in {path}"
    return len(contents) - marker_offset - len(marker)


def _latest_by_cycle(dumps):
    result = {}
    for path, dump in dumps:
        result[dump["cycle"]] = (path, dump)
    return {cycle: item[1] for cycle, item in result.items()}


def _final_dump(run_dir: Path, variable: str):
    return _read_dumps(run_dir, variable)[-1][1]


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
    values = np.stack(dump["mb_data"][variable], axis=0)
    return values[_canonical_order(dump)]


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
    np.testing.assert_allclose(
        np.asarray(candidate["mb_geometry"])[candidate_order],
        np.asarray(reference["mb_geometry"])[reference_order],
        rtol=0.0,
        atol=0.0,
    )


def _relative_level_counts(dump):
    levels = np.asarray(dump["mb_logical"])[:, 3].astype(int)
    root_level = int(np.min(levels))
    unique, counts = np.unique(levels - root_level, return_counts=True)
    return dict(zip(unique.tolist(), counts.tolist()))


def _fine_x1_center(dump):
    logical = np.asarray(dump["mb_logical"])
    root_level = int(np.min(logical[:, 3]))
    fine = logical[:, 3] == root_level + 1
    assert np.count_nonzero(fine) == 8
    geometry = np.asarray(dump["mb_geometry"])[fine]
    return float(np.mean(0.5 * (geometry[:, 0] + geometry[:, 1])))


def _fine_logical_locations(dump):
    logical = np.asarray(dump["mb_logical"])
    root_level = int(np.min(logical[:, 3]))
    return {
        tuple(location.tolist())
        for location in logical[logical[:, 3] == root_level + 1]
    }


def _assert_scripted_topology(dumps_by_cycle):
    expected_centers = {
        1: -0.375,
        2: -0.125,
        3: 0.125,
        5: 0.375,
        6: -0.375,
    }
    for cycle, center in expected_centers.items():
        assert _relative_level_counts(dumps_by_cycle[cycle]) == {0: 6, 1: 8}
        assert _fine_x1_center(dumps_by_cycle[cycle]) == pytest.approx(
            center, rel=0.0, abs=2.0e-14
        )
    assert _relative_level_counts(dumps_by_cycle[4]) == {0: 8}

    # Equal block counts but disjoint fine locations prove that each move both
    # derefined the old column and refined the new one at the same root endpoint.
    for old_cycle, new_cycle in ((1, 2), (2, 3), (5, 6)):
        assert _fine_logical_locations(dumps_by_cycle[old_cycle]).isdisjoint(
            _fine_logical_locations(dumps_by_cycle[new_cycle])
        )


def _subcycling_rank_map(dump, num_ranks: int):
    """Reproduce weighted Mesh::LoadBalance for a strict-subcycling mesh."""
    levels = np.asarray(dump["mb_logical"])[:, 3].astype(int)
    levels -= int(np.min(levels))
    block_costs = np.ldexp(np.ones(dump["n_mbs"], dtype=np.float64), levels)
    rank_map = np.empty(dump["n_mbs"], dtype=int)
    segment_end = dump["n_mbs"]
    remaining_cost = float(np.sum(block_costs, dtype=np.float64))
    for rank in range(num_ranks - 1, 0, -1):
        segment_begin = segment_end - 1
        segment_cost = float(block_costs[segment_begin])
        target_cost = remaining_cost / (rank + 1)
        while segment_begin > rank:
            candidate_cost = segment_cost + float(block_costs[segment_begin - 1])
            if abs(candidate_cost - target_cost) > abs(segment_cost - target_cost):
                break
            segment_begin -= 1
            segment_cost = candidate_cost
        rank_map[segment_begin:segment_end] = rank
        segment_end = segment_begin
        remaining_cost = max(0.0, remaining_cost - segment_cost)
    rank_map[:segment_end] = 0

    unique_ranks, counts = np.unique(rank_map, return_counts=True)
    np.testing.assert_array_equal(unique_ranks, np.arange(num_ranks))
    assert np.all(counts > 0)
    return rank_map


def _blocks_share_face(first_geometry, second_geometry, atol=1.0e-12):
    touches = []
    overlaps = []
    for direction in range(3):
        first_lo, first_hi = first_geometry[2 * direction:2 * direction + 2]
        second_lo, second_hi = second_geometry[2 * direction:2 * direction + 2]
        touches.append(
            abs(first_hi - second_lo) <= atol or abs(second_hi - first_lo) <= atol
        )
        overlaps.append(min(first_hi, second_hi) - max(first_lo, second_lo) > atol)
    if sum(touches) != 1:
        return False
    return all(overlaps[direction] for direction in range(3) if not touches[direction])


def _count_remote_coarse_fine_faces(dump, num_ranks: int):
    logical = np.asarray(dump["mb_logical"])
    geometry = np.asarray(dump["mb_geometry"])
    root_level = int(np.min(logical[:, 3]))
    rank_map = _subcycling_rank_map(dump, num_ranks)
    coarse_gids = np.flatnonzero(logical[:, 3] == root_level)
    fine_gids = np.flatnonzero(logical[:, 3] == root_level + 1)
    return sum(
        rank_map[coarse_gid] != rank_map[fine_gid]
        and _blocks_share_face(geometry[coarse_gid], geometry[fine_gid])
        for coarse_gid in coarse_gids
        for fine_gid in fine_gids
    )


def _assert_finite(dump):
    for variable in dump["var_names"]:
        assert np.isfinite(_canonical_field(dump, variable)).all(), (
            f"Non-finite values in {variable}"
        )


def _assert_state_close(
    reference, candidate, description, rtol=3.0e-6, atol=3.0e-7
):
    _assert_same_topology(reference, candidate)
    assert candidate["var_names"] == reference["var_names"]
    for variable in reference["var_names"]:
        np.testing.assert_allclose(
            _canonical_field(candidate, variable),
            _canonical_field(reference, variable),
            rtol=rtol,
            atol=atol,
            err_msg=f"{description} changed {variable}",
        )


def _assert_divb_close(reference, candidate, description):
    _assert_same_topology(reference, candidate)
    reference_divb = _canonical_field(reference, "divb")
    candidate_divb = _canonical_field(candidate, "divb")
    np.testing.assert_allclose(
        candidate_divb,
        reference_divb,
        rtol=0.0,
        atol=MAX_DIVB,
        err_msg=f"{description} changed div(B)",
    )
    assert np.max(np.abs(reference_divb)) < MAX_DIVB
    assert np.max(np.abs(candidate_divb)) < MAX_DIVB


def _assert_all_divb_bounded(*run_dirs):
    """Catch magnetic-constraint damage at the first dump, not only at the endpoint."""
    for run_dir in run_dirs:
        for path, dump in _read_dumps(run_dir, DIVB_VARIABLE):
            divb = _canonical_field(dump, "divb")
            assert np.isfinite(divb).all(), f"Non-finite div(B) in {path}"
            max_divb = float(np.max(np.abs(divb)))
            assert max_divb < MAX_DIVB, (
                f"Magnetic constraint violation in {path}: "
                f"max|div(B)|={max_divb:.16e}"
            )


def _meshblock_updates(log_output):
    matches = re.findall(r"MeshBlock-cycles\s*=\s*(\d+)", log_output)
    assert matches, "Run log did not report MeshBlock-cycles"
    return int(matches[-1])


def _assert_history_conservation(*run_dirs):
    histories = []
    for run_dir in run_dirs:
        history_path = run_dir / f"{BASENAME}.mhd.hst"
        assert history_path.is_file(), f"Missing history file {history_path}"
        histories.append(athena_read.hst(str(history_path), raw=True))

    for variable in ("mass", "1-mom", "2-mom", "3-mom", "tot-E"):
        values = np.concatenate([history[variable] for history in histories])
        assert np.isfinite(values).all(), f"Non-finite history values in {variable}"
        scale = max(float(np.max(np.abs(values))), 1.0e-30)
        relative_drift = float(np.ptp(values) / scale)
        assert relative_drift < 5.0e-9, (
            f"Conservation drift in {variable} is {relative_drift:.3e}"
        )


def test_subcycling_survives_moving_amr_mpi_and_restart(tmp_path):
    jcon_rejected_dir = tmp_path / "jcon_rejected"
    one_dir = tmp_path / "one_rank"
    three_dir = tmp_path / "three_ranks"
    split_dir = tmp_path / "split_three_ranks"
    restart_dir = tmp_path / "restart_two_ranks"
    counter_full_dir = tmp_path / "counter_full_three_ranks"
    counter_split_dir = tmp_path / "counter_split_three_ranks"
    counter_restart_dir = tmp_path / "counter_restart_two_ranks"
    rank_file_split_dir = tmp_path / "rank_file_split_three_ranks"
    rank_file_restart_dir = tmp_path / "rank_file_restart_three_ranks"
    rank_file_mismatch_dir = tmp_path / "rank_file_restart_two_ranks"
    rank_file_identity_dir = tmp_path / "rank_file_restart_mixed_identity"
    rank_file_truncated_dir = tmp_path / "rank_file_restart_truncated_rank"
    rank_file_missing_dir = tmp_path / "rank_file_restart_missing_rank"

    _run_timed_mpi_expect_failure(
        jcon_rejected_dir,
        ranks=1,
        restart=None,
        overrides=("output2/variable=mhd_jcon",),
        expected_patterns=(
            r"mhd_jcon output is not time-consistent with level subcycling",
        ),
    )
    one_log = _run_timed_mpi(one_dir, ranks=1)
    three_log = _run_timed_mpi(three_dir, ranks=3)
    split_log = _run_timed_mpi(split_dir, ranks=3, nlim=3)

    split_restarts = sorted((split_dir / "rst").glob(f"{BASENAME}.*.rst"))
    assert len(split_restarts) >= 2, "Expected initial and cycle-3 checkpoints"
    restart_log = _run_timed_mpi(
        restart_dir, ranks=2, restart=split_restarts[-1], nlim=6
    )

    # A two-cycle derefinement gate makes restart equivalence depend on restoring
    # ncyc_since_ref rather than merely reconstructing the same instantaneous tree.
    counter_overrides = ("mesh_refinement/refinement_interval=2",)
    counter_full_log = _run_timed_mpi(
        counter_full_dir, ranks=3, overrides=counter_overrides
    )
    counter_split_log = _run_timed_mpi(
        counter_split_dir, ranks=3, nlim=2, overrides=counter_overrides
    )
    counter_restarts = sorted(
        (counter_split_dir / "rst").glob(f"{BASENAME}.*.rst")
    )
    assert counter_restarts, "Expected a counter-gated split checkpoint"
    counter_restart_log = _run_timed_mpi(
        counter_restart_dir,
        ranks=2,
        restart=counter_restarts[-1],
        nlim=6,
        overrides=counter_overrides,
    )

    rank_file_overrides = (
        "output4/single_file_per_rank=true",
        "output4/dcycle=1",
    )
    rank_file_split_log = _run_timed_mpi(
        rank_file_split_dir, ranks=3, nlim=3, overrides=rank_file_overrides
    )
    rank_file_restarts = sorted(
        (rank_file_split_dir / "rst/rank_00000000").glob(
            f"{BASENAME}.*.rst"
        )
    )
    assert len(rank_file_restarts) >= 4, (
        "Expected one per-rank split checkpoint for each cycle from 0 through 3"
    )
    rank_file_restart_log = _run_timed_mpi(
        rank_file_restart_dir,
        ranks=3,
        restart=rank_file_restarts[-1],
        nlim=6,
        overrides=rank_file_overrides,
    )

    _run_timed_mpi_expect_failure(
        rank_file_mismatch_dir,
        ranks=2,
        restart=rank_file_restarts[-1],
        overrides=rank_file_overrides,
        expected_patterns=(
            r"Per-rank restart writer/current MPI rank layout is incompatible",
        ),
    )

    # Cycles 1 and 2 both have 14 MeshBlocks and the same three-rank writer
    # partition.  Retain both complete sets, then make an independent Frankenstein
    # set whose rank 1 field stream comes from the other checkpoint generation.
    identity_source = next(
        path for path in rank_file_restarts if path.name.endswith(".00001.rst")
    )
    identity_target = next(
        path for path in rank_file_restarts if path.name.endswith(".00002.rst")
    )
    mixed_rst_dir = tmp_path / "mixed_rank_checkpoint" / "rst"
    for rank in range(3):
        identity_source_path = (
            rank_file_split_dir
            / "rst"
            / f"rank_{rank:08d}"
            / identity_source.name
        )
        identity_target_path = (
            rank_file_split_dir
            / "rst"
            / f"rank_{rank:08d}"
            / identity_target.name
        )
        assert (
            _restart_binary_payload_size(identity_source_path)
            == _restart_binary_payload_size(identity_target_path)
        ), f"Cycles 1 and 2 have different writer partitions on rank {rank}"
        rank_dir = mixed_rst_dir / f"rank_{rank:08d}"
        rank_dir.mkdir(parents=True)
        shutil.copy2(identity_target_path, rank_dir / identity_target.name)
    shutil.copy2(
        rank_file_split_dir
        / "rst"
        / "rank_00000001"
        / identity_source.name,
        mixed_rst_dir / "rank_00000001" / identity_target.name,
    )
    _run_timed_mpi_expect_failure(
        rank_file_identity_dir,
        ranks=3,
        restart=mixed_rst_dir / "rank_00000000" / identity_target.name,
        overrides=rank_file_overrides,
        expected_patterns=(
            r"Per-rank checkpoint identity mismatch across rank files",
        ),
    )

    truncated_rank = 1
    truncated_rank_restart = (
        rank_file_split_dir
        / "rst"
        / f"rank_{truncated_rank:08d}"
        / rank_file_restarts[-1].name
    )
    original_size = truncated_rank_restart.stat().st_size
    assert original_size > 4096
    with truncated_rank_restart.open("r+b") as restart_file:
        restart_file.truncate(original_size - 4096)
    _run_timed_mpi_expect_failure(
        rank_file_truncated_dir,
        ranks=3,
        restart=rank_file_restarts[-1],
        overrides=rank_file_overrides,
        expected_patterns=(
            rf"Per-rank file read was truncated on rank {truncated_rank}:",
        ),
    )

    missing_rank = 2
    missing_rank_restart = (
        rank_file_split_dir
        / "rst"
        / f"rank_{missing_rank:08d}"
        / rank_file_restarts[-1].name
    )
    assert missing_rank_restart.is_file(), (
        f"Expected rank-{missing_rank} checkpoint {missing_rank_restart}"
    )
    missing_rank_restart.unlink()
    assert not missing_rank_restart.exists()
    _run_timed_mpi_expect_failure(
        rank_file_missing_dir,
        ranks=3,
        restart=rank_file_restarts[-1],
        overrides=rank_file_overrides,
        expected_patterns=(
            rf"Error: Unable to open restart file on rank {missing_rank}:",
            re.escape(str(missing_rank_restart.resolve())),
        ),
    )

    assert _meshblock_updates(one_log) == EXPECTED_FULL_UPDATES
    assert _meshblock_updates(three_log) == EXPECTED_FULL_UPDATES
    assert _meshblock_updates(split_log) == EXPECTED_HALF_UPDATES
    assert _meshblock_updates(restart_log) == EXPECTED_HALF_UPDATES
    assert _meshblock_updates(counter_full_log) == (
        _meshblock_updates(counter_split_log)
        + _meshblock_updates(counter_restart_log)
    )
    assert _meshblock_updates(three_log) == (
        _meshblock_updates(rank_file_split_log)
        + _meshblock_updates(rank_file_restart_log)
    )

    _assert_all_divb_bounded(
        one_dir,
        three_dir,
        split_dir,
        restart_dir,
        counter_full_dir,
        counter_split_dir,
        counter_restart_dir,
        rank_file_split_dir,
        rank_file_restart_dir,
    )

    # Cycle 2 is the first simultaneous move: the old fine column is de-refined while
    # a new column is refined.  This is the earliest point at which stale coarse_b0
    # boundary scratch used to be copied onto the new parent's active faces.
    three_divb_by_cycle = _latest_by_cycle(
        _read_dumps(three_dir, DIVB_VARIABLE)
    )
    assert 2 in three_divb_by_cycle
    cycle_two_divb = _canonical_field(three_divb_by_cycle[2], "divb")
    assert np.max(np.abs(cycle_two_divb)) < MAX_DIVB

    amr_counts = re.findall(
        r"(\d+) MeshBlocks created, (\d+) deleted by AMR", three_log
    )
    assert amr_counts and tuple(map(int, amr_counts[-1])) == (30, 24)
    communicated = re.findall(
        r"(\d+) communicated for load balancing", three_log
    )
    assert communicated and int(communicated[-1]) > 0

    one_states_by_cycle = _latest_by_cycle(_read_dumps(one_dir, STATE_VARIABLE))
    three_states_by_cycle = _latest_by_cycle(_read_dumps(three_dir, STATE_VARIABLE))
    split_states_by_cycle = _latest_by_cycle(_read_dumps(split_dir, STATE_VARIABLE))
    restart_states_by_cycle = _latest_by_cycle(_read_dumps(restart_dir, STATE_VARIABLE))
    assert set(range(7)).issubset(one_states_by_cycle)
    assert set(range(7)).issubset(three_states_by_cycle)
    assert {0, 1, 2, 3}.issubset(split_states_by_cycle)
    assert {4, 5, 6}.issubset(restart_states_by_cycle)
    _assert_scripted_topology(one_states_by_cycle)
    _assert_scripted_topology(three_states_by_cycle)
    _assert_scripted_topology(
        {**split_states_by_cycle, **restart_states_by_cycle}
    )

    remote_faces = [
        _count_remote_coarse_fine_faces(three_states_by_cycle[cycle], num_ranks=3)
        for cycle in (1, 2, 3, 5, 6)
    ]
    assert sum(count > 0 for count in remote_faces) >= 2, (
        f"Moving refinement did not exercise enough remote coarse/fine faces: "
        f"{remote_faces}"
    )

    one_state = _final_dump(one_dir, STATE_VARIABLE)
    three_state = _final_dump(three_dir, STATE_VARIABLE)
    split_state = split_states_by_cycle[3]
    restart_state = _final_dump(restart_dir, STATE_VARIABLE)
    one_divb = _final_dump(one_dir, DIVB_VARIABLE)
    three_divb = _final_dump(three_dir, DIVB_VARIABLE)
    restart_divb = _final_dump(restart_dir, DIVB_VARIABLE)
    for dump in (one_state, three_state, split_state, restart_state,
                 one_divb, three_divb, restart_divb):
        _assert_finite(dump)

    assert three_state["cycle"] == one_state["cycle"] == 6
    assert restart_state["cycle"] == 6
    assert three_state["time"] == pytest.approx(one_state["time"], rel=2.0e-14)
    assert restart_state["time"] == pytest.approx(one_state["time"], rel=2.0e-14)
    _assert_state_close(
        one_state,
        three_state,
        "MPI decomposition",
        rtol=8.0e-5,
        atol=1.0e-6,
    )
    _assert_state_close(
        three_states_by_cycle[3], split_state, "three-rank early termination"
    )
    _assert_state_close(three_state, restart_state, "3-to-2-rank restart")
    _assert_divb_close(one_divb, three_divb, "MPI decomposition")
    _assert_divb_close(three_divb, restart_divb, "3-to-2-rank restart")

    counter_full_state = _final_dump(counter_full_dir, STATE_VARIABLE)
    counter_restart_state = _final_dump(counter_restart_dir, STATE_VARIABLE)
    counter_full_divb = _final_dump(counter_full_dir, DIVB_VARIABLE)
    counter_restart_divb = _final_dump(counter_restart_dir, DIVB_VARIABLE)
    _assert_state_close(
        counter_full_state,
        counter_restart_state,
        "AMR counter-gated 3-to-2-rank restart",
        rtol=8.0e-5,
        atol=1.0e-6,
    )
    _assert_divb_close(
        counter_full_divb,
        counter_restart_divb,
        "AMR counter-gated 3-to-2-rank restart",
    )

    rank_file_restart_state = _final_dump(rank_file_restart_dir, STATE_VARIABLE)
    rank_file_restart_divb = _final_dump(rank_file_restart_dir, DIVB_VARIABLE)
    _assert_state_close(
        three_state,
        rank_file_restart_state,
        "single-file-per-rank restart",
    )
    _assert_divb_close(
        three_divb,
        rank_file_restart_divb,
        "single-file-per-rank restart",
    )

    _assert_history_conservation(one_dir)
    _assert_history_conservation(three_dir)
    _assert_history_conservation(split_dir, restart_dir)
    _assert_history_conservation(counter_full_dir)
    _assert_history_conservation(counter_split_dir, counter_restart_dir)
    _assert_history_conservation(rank_file_split_dir, rank_file_restart_dir)

    total_bytes = sum(
        path.stat().st_size for path in tmp_path.rglob("*") if path.is_file()
    )
    assert total_bytes < 64 * 1024 * 1024, (
        f"Regression generated {total_bytes / 1024**2:.1f} MiB, "
        f"exceeding its storage budget"
    )
