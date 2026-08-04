"""MPI decomposition regression for strict 2:1 level subcycling in dynGRMHD."""

import re
from pathlib import Path

import numpy as np
import pytest

import test_suite.testutils as testutils
import bin_convert


INPUT_FILE = "inputs/subcycling_smr_dyngrmhd.athinput"
BASENAME = "subcycling_smr_dyngrmhd"
STATE_VARIABLE = "mhd_w_bcc"
DIVB_VARIABLE = "mhd_divb"
EXPECTED_ROOT_CYCLE_UPDATES = 424


def _run_case(run_dir: Path, ranks: int):
    """Run one decomposition and return its final state and div(B) dumps."""
    log_path = Path(testutils.LOG_FILE_PATH)
    log_offset = log_path.stat().st_size if log_path.exists() else 0
    flags = ["-d", str(run_dir), f"job/basename={BASENAME}"]
    assert testutils.mpi_run(INPUT_FILE, flags, threads=ranks)
    with log_path.open("rb") as log_file:
        log_file.seek(log_offset)
        log_output = log_file.read().decode("utf-8", errors="replace")
    update_counts = re.findall(r"MeshBlock-cycles\s*=\s*(\d+)", log_output)
    assert update_counts, "Run log did not report MeshBlock-cycles"
    return (
        _read_final_dump(run_dir, STATE_VARIABLE),
        _read_final_dump(run_dir, DIVB_VARIABLE),
        int(update_counts[-1]),
    )


def _read_final_dump(run_dir: Path, variable: str):
    """Read the synchronized dump with the largest cycle number."""
    pattern = f"{BASENAME}.{variable}.*.bin"
    paths = sorted((run_dir / "bin").glob(pattern))
    assert len(paths) >= 2, f"Expected initial and final dumps matching {pattern}"
    dumps = [bin_convert.read_binary(str(path)) for path in paths]
    return max(dumps, key=lambda dump: (dump["cycle"], dump["time"]))


def _canonical_order(dump):
    """Order blocks by level and logical location, independent of output order."""
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
    """Check global mesh metadata and leaf-block topology exactly."""
    for name in ("Nx1", "Nx2", "Nx3", "nx1_mb", "nx2_mb", "nx3_mb", "n_mbs"):
        assert candidate[name] == reference[name], f"Topology differs in {name}"

    ref_order = _canonical_order(reference)
    candidate_order = _canonical_order(candidate)
    np.testing.assert_array_equal(
        np.asarray(candidate["mb_logical"])[candidate_order],
        np.asarray(reference["mb_logical"])[ref_order],
    )
    np.testing.assert_array_equal(
        np.asarray(candidate["mb_index"])[candidate_order],
        np.asarray(reference["mb_index"])[ref_order],
    )
    np.testing.assert_allclose(
        np.asarray(candidate["mb_geometry"])[candidate_order],
        np.asarray(reference["mb_geometry"])[ref_order],
        rtol=0.0,
        atol=0.0,
    )


def _subcycling_rank_map(dump, num_ranks: int):
    """Reproduce weighted Mesh::LoadBalance for a strict-subcycling mesh."""
    levels = np.asarray(dump["mb_logical"])[:, 3].astype(int)
    block_costs = np.ldexp(np.ones(dump["n_mbs"], dtype=np.float32), levels)
    num_blocks = dump["n_mbs"]
    rank_map = np.empty(num_blocks, dtype=int)
    rank = num_ranks - 1
    remaining_cost = np.float32(0.0)
    for cost in block_costs:
        remaining_cost = np.float32(remaining_cost + cost)
    target_cost = np.float32(remaining_cost / np.float32(num_ranks))
    rank_cost = np.float32(0.0)

    for gid in range(num_blocks - 1, -1, -1):
        rank_cost = np.float32(rank_cost + block_costs[gid])
        rank_map[gid] = rank
        if rank_cost >= target_cost and rank > 0:
            rank -= 1
            remaining_cost = np.float32(remaining_cost - rank_cost)
            rank_cost = np.float32(0.0)
            target_cost = np.float32(remaining_cost / np.float32(rank + 1))

    assert rank == 0
    return rank_map


def _blocks_share_face(coarse_geometry, fine_geometry, atol=1.0e-12):
    """Return true for a positive-area face shared by two leaf blocks."""
    touches = []
    overlaps = []
    for direction in range(3):
        clo, chi = coarse_geometry[2 * direction : 2 * direction + 2]
        flo, fhi = fine_geometry[2 * direction : 2 * direction + 2]
        touches.append(abs(chi - flo) <= atol or abs(fhi - clo) <= atol)
        overlaps.append(min(chi, fhi) - max(clo, flo) > atol)

    if sum(touches) != 1:
        return False
    return all(overlaps[direction] for direction in range(3) if not touches[direction])


def _count_remote_coarse_fine_faces(dump, num_ranks: int, coarse_level: int):
    """Infer remote faces for one adjacent level pair from global GID ordering."""
    logical = np.asarray(dump["mb_logical"])
    geometry = np.asarray(dump["mb_geometry"])
    rank_map = _subcycling_rank_map(dump, num_ranks)
    coarse_gids = np.flatnonzero(logical[:, 3] == coarse_level)
    fine_gids = np.flatnonzero(logical[:, 3] == coarse_level + 1)

    remote_faces = 0
    for coarse_gid in coarse_gids:
        for fine_gid in fine_gids:
            if rank_map[coarse_gid] == rank_map[fine_gid]:
                continue
            if _blocks_share_face(geometry[coarse_gid], geometry[fine_gid]):
                remote_faces += 1
    return remote_faces


def _assert_finite(dump):
    for variable in dump["var_names"]:
        values = _canonical_field(dump, variable)
        assert np.isfinite(values).all(), f"Non-finite values in {variable}"


def test_subcycling_is_mpi_decomposition_invariant(tmp_path):
    """Compare one and three ranks across deliberately remote coarse/fine faces."""
    one_state, one_divb, one_updates = _run_case(tmp_path / "one_rank", ranks=1)
    three_state, three_divb, three_updates = _run_case(
        tmp_path / "three_ranks", ranks=3
    )

    for dump in (one_state, one_divb, three_state, three_divb):
        assert dump["cycle"] == 1
        assert dump["time"] > 0.0
        _assert_finite(dump)

    assert three_state["time"] == pytest.approx(one_state["time"], rel=2.0e-14)
    assert three_divb["time"] == pytest.approx(one_divb["time"], rel=2.0e-14)

    _assert_same_topology(one_state, one_divb)
    _assert_same_topology(one_state, three_state)
    _assert_same_topology(one_state, three_divb)

    levels, counts = np.unique(one_state["mb_logical"][:, 3], return_counts=True)
    level_counts = dict(zip(levels.tolist(), counts.tolist()))
    assert level_counts == {0: 56, 1: 56, 2: 64}
    assert one_state["n_mbs"] == 176
    theoretical_updates = sum(count * (1 << level) for level, count in level_counts.items())
    assert theoretical_updates == EXPECTED_ROOT_CYCLE_UPDATES
    assert one_updates == EXPECTED_ROOT_CYCLE_UPDATES
    assert three_updates == EXPECTED_ROOT_CYCLE_UPDATES

    expected_remote_faces = {0: 6, 1: 7}
    for coarse_level, expected_faces in expected_remote_faces.items():
        remote_faces = _count_remote_coarse_fine_faces(
            one_state, num_ranks=3, coarse_level=coarse_level
        )
        assert remote_faces == expected_faces, (
            f"Unexpected remote level-{coarse_level}/level-{coarse_level + 1} "
            f"face count: {remote_faces}"
        )

    assert three_state["var_names"] == one_state["var_names"]
    for variable in one_state["var_names"]:
        np.testing.assert_allclose(
            _canonical_field(three_state, variable),
            _canonical_field(one_state, variable),
            rtol=3.0e-6,
            atol=3.0e-7,
            err_msg=f"MPI decomposition changed {variable}",
        )

    one_divb_values = _canonical_field(one_divb, "divb")
    three_divb_values = _canonical_field(three_divb, "divb")
    np.testing.assert_allclose(
        three_divb_values,
        one_divb_values,
        rtol=0.0,
        atol=1.0e-11,
        err_msg="MPI decomposition changed div(B)",
    )
    assert np.max(np.abs(one_divb_values)) < 1.0e-10
    assert np.max(np.abs(three_divb_values)) < 1.0e-10
