"""Moving BBH region output regression across trajectory phase and MPI layout."""

import subprocess
from pathlib import Path

import numpy as np

import bin_convert
import test_suite.testutils as testutils
from plot_bbh_grmhd import read_slice


INPUT_FILE = "../../../inputs/dyngr/effective_bbh_torus_smoke.athinput"
OMEGA = 0.01118033988749895
QUARTER_ORBIT_TIME = 0.5 * np.pi / OMEGA
HALF_WIDTH = 7.0


def _run_case(
    run_dir: Path, basename: str, ranks: int, trajectory_offset: float
) -> tuple[dict, Path]:
    flags = [
        "-d",
        str(run_dir),
        f"job/basename={basename}",
        "mesh/nx1=32",
        "mesh/nx2=32",
        "mesh/nx3=32",
        "meshblock/nx1=8",
        "meshblock/nx2=8",
        "meshblock/nx3=8",
        "time/nlim=0",
        "problem/history_inner_radius=20",
        f"problem/trajectory_time_offset={trajectory_offset:.17g}",
        "output6/region_center=bh1",
        f"output6/region_half_width={HALF_WIDTH:.17g}",
        "output3/dcycle=0",
        "output4/dcycle=0",
    ]
    assert testutils.mpi_run(INPUT_FILE, flags, threads=ranks)
    path = min((run_dir / "bin").glob(f"{basename}.bbh_local_w.*.bin"))
    return bin_convert.read_binary(str(path)), path


def _assert_region(data: dict, center: tuple[float, float, float]) -> None:
    assert data["n_mbs"] > 0
    assert data["n_mbs"] < 120
    assert (data["nx1_out_mb"], data["nx2_out_mb"], data["nx3_out_mb"]) == (
        8,
        8,
        1,
    )
    geometry = data["mb_geometry"]
    for axis in range(3):
        lower = geometry[:, 2 * axis]
        upper = geometry[:, 2 * axis + 1]
        assert np.all(upper > center[axis] - HALF_WIDTH)
        assert np.all(lower < center[axis] + HALF_WIDTH)
    assert np.all(geometry[:, 4] <= center[2])
    assert np.all(geometry[:, 5] > center[2])


def test_bbh_region_tracks_phase_and_is_mpi_decomposition_invariant(tmp_path):
    phase0, phase0_path = _run_case(
        tmp_path / "phase0", "moving_region_phase0", 1, 0.0
    )
    phase90, _ = _run_case(
        tmp_path / "phase90", "moving_region_phase90", 1, QUARTER_ORBIT_TIME
    )
    mpi2, _ = _run_case(tmp_path / "mpi2", "moving_region_mpi2", 2, 0.0)

    _assert_region(phase0, (10.0, 0.0, 0.0))
    _assert_region(phase90, (0.0, 10.0, 0.0))
    assert not np.array_equal(phase0["mb_logical"], phase90["mb_logical"])
    sliced = read_slice(phase0_path, ["dens", "bmag", "level"], "z", 0.0, None)
    assert len(sliced.blocks) == phase0["n_mbs"]
    assert {block.slice_shape for block in sliced.blocks} == {(8, 8)}

    assert phase0["n_mbs"] == mpi2["n_mbs"]
    np.testing.assert_array_equal(mpi2["mb_index"], phase0["mb_index"])
    np.testing.assert_array_equal(mpi2["mb_logical"], phase0["mb_logical"])
    np.testing.assert_array_equal(mpi2["mb_geometry"], phase0["mb_geometry"])
    for variable in phase0["var_names"]:
        np.testing.assert_array_equal(
            np.asarray(mpi2["mb_data"][variable]),
            np.asarray(phase0["mb_data"][variable]),
        )


def test_bbh_region_contract_and_selection_survive_restart(tmp_path):
    initial_dir = tmp_path / "initial"
    basename = "moving_region_restart"
    flags = [
        "-d",
        str(initial_dir),
        f"job/basename={basename}",
        "mesh/nx1=8",
        "mesh/nx2=8",
        "mesh/nx3=8",
        "meshblock/nx1=4",
        "meshblock/nx2=8",
        "meshblock/nx3=8",
        "time/nlim=0",
        "time/tlim=1",
        "problem/initial_data=uniform",
        "problem/history_inner_radius=20",
        "output2/dcycle=0",
        "output3/dcycle=0",
        "output4/dcycle=1",
        "output6/region_center=bh1",
        f"output6/region_half_width={HALF_WIDTH:.17g}",
    ]
    assert testutils.mpi_run(INPUT_FILE, flags, threads=1)

    checkpoint = max((initial_dir / "rst").glob(f"{basename}.*.rst"))
    before_path = max(
        (initial_dir / "bin").glob(f"{basename}.bbh_local_w.*.bin")
    )
    restart_dir = tmp_path / "restart"
    command = [
        "./athena",
        "-r",
        str(checkpoint.resolve()),
        "-d",
        str(restart_dir.resolve()),
        "time/nlim=0",
        "time/tlim=1",
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
        log_file.write(f"\nMoving-region restart command: {' '.join(command)}\n")
        log_file.write(result.stdout)
    assert result.returncode == 0, result.stdout[-8000:]

    after_path = max(
        (restart_dir / "bin").glob(f"{basename}.bbh_local_w.*.bin")
    )
    before = bin_convert.read_binary(str(before_path))
    after = bin_convert.read_binary(str(after_path))
    assert after["time"] == before["time"] == 0.0
    np.testing.assert_array_equal(after["mb_index"], before["mb_index"])
    np.testing.assert_array_equal(after["mb_logical"], before["mb_logical"])
    np.testing.assert_array_equal(after["mb_geometry"], before["mb_geometry"])
    for variable in before["var_names"]:
        np.testing.assert_array_equal(
            np.asarray(after["mb_data"][variable]),
            np.asarray(before["mb_data"][variable]),
        )


def test_multiple_derived_outputs_release_scratch_and_finalize_once(tmp_path):
    """Global and moving derived outputs must coexist without terminal duplicates."""
    run_dir = tmp_path / "derived_lifecycle"
    basename = "moving_region_derived_lifecycle"
    flags = [
        "-d",
        str(run_dir),
        f"job/basename={basename}",
        "mesh/nx1=8",
        "mesh/nx2=8",
        "mesh/nx3=8",
        "meshblock/nx1=4",
        "meshblock/nx2=8",
        "meshblock/nx3=8",
        "time/nlim=1",
        "time/tlim=1",
        "problem/initial_data=uniform",
        "problem/history_inner_radius=20",
        "output2/dcycle=0",
        "output3/variable=mhd_gr_diagnostics",
        "output3/dcycle=1",
        "output4/dcycle=0",
        "output6/variable=mhd_gr_diagnostics",
        "output6/id=bbh_local_gr",
        "output6/region_center=bbh_com",
        f"output6/region_half_width={HALF_WIDTH:.17g}",
        "output6/dcycle=1",
    ]
    assert testutils.mpi_run(INPUT_FILE, flags, threads=1)

    global_paths = sorted(
        (run_dir / "bin").glob(f"{basename}.mhd_gr_diagnostics.*.bin")
    )
    local_paths = sorted(
        (run_dir / "bin").glob(f"{basename}.bbh_local_gr.*.bin")
    )
    assert len(global_paths) == len(local_paths) == 2
    assert [bin_convert.read_binary(str(path))["cycle"] for path in global_paths] == [0, 1]
    assert [bin_convert.read_binary(str(path))["cycle"] for path in local_paths] == [0, 1]

    history_path = run_dir / f"{basename}.user.hst"
    rows = [
        line
        for line in history_path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert len(rows) == 2


def test_unknown_bbh_region_fails_before_creating_binary_file(tmp_path):
    run_dir = tmp_path / "unknown"
    basename = "moving_region_unknown"
    command = [
        "./athena",
        "-i",
        INPUT_FILE,
        "-d",
        str(run_dir.resolve()),
        f"job/basename={basename}",
        "mesh/nx1=8",
        "mesh/nx2=8",
        "mesh/nx3=8",
        "meshblock/nx1=4",
        "meshblock/nx2=8",
        "meshblock/nx3=8",
        "time/nlim=0",
        "problem/initial_data=uniform",
        "problem/history_inner_radius=20",
        "output2/dcycle=0",
        "output3/dcycle=0",
        "output4/dcycle=0",
        "output6/region_center=not_a_bbh_center",
        f"output6/region_half_width={HALF_WIDTH:.17g}",
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
        log_file.write(f"\nUnknown moving-region command: {' '.join(command)}\n")
        log_file.write(result.stdout)
    assert result.returncode != 0
    assert "rejected moving output region 'not_a_bbh_center'" in result.stdout
    assert not list((run_dir / "bin").glob("*.bin"))
