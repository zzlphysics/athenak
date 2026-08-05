"""MPI initialization regression for the effective-BBH reference FM torus."""

from pathlib import Path

import numpy as np
import pytest

import bin_convert
import test_suite.testutils as testutils


INPUT_FILE = "../../../inputs/dyngr/effective_bbh_torus_smoke.athinput"
STATE_VARIABLE = "mhd_w_bcc"
DIVB_VARIABLE = "mhd_divb"
COMMON_OVERRIDES = (
    "mesh/nx1=32",
    "mesh/nx2=32",
    "mesh/nx3=32",
    "meshblock/nx1=8",
    "meshblock/nx2=8",
    "meshblock/nx3=8",
    "time/nlim=0",
    "output4/dcycle=0",
)


def _run_case(run_dir: Path, basename: str, ranks: int, overrides=()):
    """Initialize one torus and return its state, div(B), and new log text."""
    log_path = Path(testutils.LOG_FILE_PATH)
    log_offset = log_path.stat().st_size if log_path.exists() else 0
    flags = ["-d", str(run_dir), f"job/basename={basename}"]
    flags.extend(COMMON_OVERRIDES)
    flags.extend(overrides)
    assert testutils.mpi_run(INPUT_FILE, flags, threads=ranks)
    with log_path.open("rb") as log_file:
        log_file.seek(log_offset)
        log_output = log_file.read().decode("utf-8", errors="replace")
    return (
        _read_first_dump(run_dir, basename, STATE_VARIABLE),
        _read_first_dump(run_dir, basename, DIVB_VARIABLE),
        log_output,
    )


def _read_first_dump(run_dir: Path, basename: str, variable: str):
    pattern = f"{basename}.{variable}.*.bin"
    paths = sorted((run_dir / "bin").glob(pattern))
    assert paths, f"No dumps matched {pattern}"
    dumps = [bin_convert.read_binary(str(path)) for path in paths]
    return min(dumps, key=lambda dump: (dump["cycle"], dump["time"]))


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


def _assert_finite(dump):
    for variable in dump["var_names"]:
        assert np.isfinite(_canonical_field(dump, variable)).all(), (
            f"Non-finite values in {variable}"
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
    np.testing.assert_allclose(
        np.asarray(candidate["mb_geometry"])[candidate_order],
        np.asarray(reference["mb_geometry"])[reference_order],
        rtol=0.0,
        atol=0.0,
    )


def test_bbh_torus_initial_data_is_mpi_decomposition_invariant(tmp_path):
    """The perturbed torus and CT field must not depend on the MPI partition."""
    one_state, one_divb, one_log = _run_case(
        tmp_path / "one_rank", "bbh_torus_one", ranks=1
    )
    two_state, two_divb, two_log = _run_case(
        tmp_path / "two_ranks", "bbh_torus_two", ranks=2
    )

    for dump in (one_state, one_divb, two_state, two_divb):
        assert dump["cycle"] == 0
        assert dump["time"] == pytest.approx(0.0, abs=0.0)
        _assert_finite(dump)
    _assert_same_topology(one_state, one_divb)
    _assert_same_topology(one_state, two_state)
    _assert_same_topology(one_state, two_divb)

    assert one_state["var_names"] == two_state["var_names"]
    for variable in one_state["var_names"]:
        np.testing.assert_array_equal(
            _canonical_field(two_state, variable),
            _canonical_field(one_state, variable),
            err_msg=f"MPI decomposition changed {variable}",
        )

    one_divb_values = _canonical_field(one_divb, "divb")
    two_divb_values = _canonical_field(two_divb, "divb")
    assert np.max(np.abs(one_divb_values)) < 1.0e-10
    assert np.max(np.abs(two_divb_values)) < 1.0e-10
    np.testing.assert_allclose(two_divb_values, one_divb_values, rtol=0.0, atol=1.0e-10)
    for log_output in (one_log, two_log):
        assert "magnetic norm=density_weighted_beta" in log_output
        assert "normalized target=100" in log_output


@pytest.mark.parametrize(
    ("case_name", "overrides", "expected_log", "expect_magnetic"),
    (
        (
            "peak",
            ("problem/torus_mag_norm=peak_beta",),
            "peak_beta unnormalized ratio=",
            True,
        ),
        (
            "pressure",
            ("problem/torus_mag_norm=integrated_pressure",),
            "integrated_pressure unnormalized ratio=",
            True,
        ),
        (
            "internal",
            ("problem/torus_mag_norm=integrated_internal_energy",),
            "integrated_internal_energy unnormalized ratio=",
            True,
        ),
        (
            "unmagnetized",
            ("problem/torus_magnetic_field=none",),
            "magnetic field=none",
            False,
        ),
    ),
)
def test_bbh_torus_magnetic_initialization_modes(
    tmp_path, case_name, overrides, expected_log, expect_magnetic
):
    """Every documented normalization and the unmagnetized branch must initialize."""
    state, divb, log_output = _run_case(
        tmp_path / case_name, f"bbh_torus_{case_name}", ranks=1,
        overrides=overrides
    )
    _assert_finite(state)
    _assert_finite(divb)
    magnetic_max = max(
        float(np.max(np.abs(_canonical_field(state, variable))))
        for variable in ("bcc1", "bcc2", "bcc3")
    )
    if expect_magnetic:
        assert magnetic_max > 0.0
        assert "normalized target=100" in log_output
    else:
        assert magnetic_max == 0.0
    assert expected_log in log_output
    assert np.max(np.abs(_canonical_field(divb, "divb"))) < 1.0e-10


def test_bbh_torus_general_reference_frame_initializes(tmp_path):
    """A scaled, translated, slowly moving reference potential remains finite."""
    overrides = (
        "problem/mass_scale1=0.5",
        "problem/mass_scale2=0.5",
        "problem/separation=10",
        "problem/omega=0.0223606797749979",
        "problem/torus_reference_mass=0.5",
        "problem/torus_reference_center1=5",
        "problem/torus_reference_center2=-3",
        "problem/torus_reference_center3=2",
        "problem/torus_reference_velocity1=0.01",
        "problem/torus_reference_velocity2=0.005",
        "problem/torus_reference_velocity3=0",
        "problem/torus_r_edge=9",
        "problem/torus_r_peak=14.5",
        "problem/torus_min_grid_peak_fraction=0.5",
        "problem/torus_min_magnetic_cells=16",
    )
    state, divb, log_output = _run_case(
        tmp_path / "general", "bbh_torus_general", ranks=1,
        overrides=overrides
    )
    _assert_finite(state)
    _assert_finite(divb)
    assert "M_ref=0.5 (trajectory total=0.5)" in log_output
    assert "center=(5,-3,2), velocity=(0.01,0.005,0)" in log_output
    assert np.max(np.abs(_canonical_field(divb, "divb"))) < 1.0e-10
