"""Metric-aware BBH user-history regression across MPI decompositions."""

from pathlib import Path

import numpy as np
import pytest

import athena_read
import bin_convert
import test_suite.testutils as testutils


INPUT_FILE = "../../../inputs/dyngr/effective_bbh_torus_smoke.athinput"
EXPECTED_COLUMNS = (
    "time",
    "dt",
    "bh1_x",
    "bh1_y",
    "bh1_z",
    "bh2_x",
    "bh2_y",
    "bh2_z",
    "bh_sep",
    "orb_omega",
    "bh1_mass",
    "bh2_mass",
    "baryon_m",
    "rho_prp",
    "pgas_prp",
    "emag_prp",
    "lor_D",
    "sigma_D",
    "angmom_z",
    "inner_D",
    "rho_max",
    "sigma_max",
)


def _run_case(
    run_dir: Path, basename: str, ranks: int
) -> tuple[dict[str, np.ndarray], dict, dict]:
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
        "problem/user_hist=true",
        "problem/history_inner_radius=20",
        "output3/variable=mhd_gr_diagnostics",
        "output4/dcycle=0",
    ]
    assert testutils.mpi_run(INPUT_FILE, flags, threads=ranks)
    path = run_dir / f"{basename}.user.hst"
    assert path.is_file()
    state_path = min((run_dir / "bin").glob(f"{basename}.mhd_w_bcc.*.bin"))
    diagnostic_path = min(
        (run_dir / "bin").glob(f"{basename}.mhd_gr_diagnostics.*.bin")
    )
    return (
        athena_read.hst(str(path)),
        bin_convert.read_binary(str(state_path)),
        bin_convert.read_binary(str(diagnostic_path)),
    )


def _last(history: dict[str, np.ndarray], name: str) -> float:
    return float(history[name][-1])


def test_bbh_user_history_is_physical_and_mpi_decomposition_invariant(tmp_path):
    """Trajectory metadata and GRMHD integrals must not be multiplied by rank count."""
    one, one_state, one_diagnostics = _run_case(
        tmp_path / "one", "bbh_history_one", ranks=1
    )
    two, _, _ = _run_case(tmp_path / "two", "bbh_history_two", ranks=2)

    assert tuple(one) == EXPECTED_COLUMNS
    assert tuple(two) == EXPECTED_COLUMNS
    assert one["time"].size == two["time"].size == 1
    for history in (one, two):
        for column in EXPECTED_COLUMNS:
            assert np.isfinite(history[column]).all(), column
        assert _last(history, "time") == pytest.approx(0.0, abs=0.0)
        assert _last(history, "bh1_x") == pytest.approx(10.0, abs=1.0e-14)
        assert _last(history, "bh2_x") == pytest.approx(-10.0, abs=1.0e-14)
        for column in ("bh1_y", "bh1_z", "bh2_y", "bh2_z"):
            assert _last(history, column) == pytest.approx(0.0, abs=1.0e-14)
        assert _last(history, "bh_sep") == pytest.approx(20.0, abs=1.0e-14)
        assert _last(history, "orb_omega") == pytest.approx(
            0.01118033988749895, rel=1.0e-14
        )
        assert _last(history, "bh1_mass") == pytest.approx(0.5, abs=1.0e-14)
        assert _last(history, "bh2_mass") == pytest.approx(0.5, abs=1.0e-14)
        assert _last(history, "baryon_m") > 0.0
        assert _last(history, "rho_prp") > 0.0
        assert _last(history, "pgas_prp") > 0.0
        assert _last(history, "emag_prp") > 0.0
        assert _last(history, "lor_D") >= _last(history, "baryon_m")
        assert _last(history, "sigma_D") >= 0.0
        assert 0.0 < _last(history, "inner_D") < _last(history, "baryon_m")
        assert _last(history, "rho_max") >= 0.9
        assert _last(history, "sigma_max") > 0.0

    density_max = max(
        float(np.max(values)) for values in one_state["mb_data"]["dens"]
    )
    sigma_max = max(
        float(np.max(values))
        for values in one_diagnostics["mb_data"]["gr_sigma"]
    )
    assert _last(one, "rho_max") == pytest.approx(density_max, rel=5.0e-7)
    assert _last(one, "sigma_max") == pytest.approx(sigma_max, rel=5.0e-7)

    for column in EXPECTED_COLUMNS:
        np.testing.assert_allclose(
            two[column], one[column], rtol=5.0e-13, atol=5.0e-14,
            err_msg=f"MPI decomposition changed {column}",
        )
