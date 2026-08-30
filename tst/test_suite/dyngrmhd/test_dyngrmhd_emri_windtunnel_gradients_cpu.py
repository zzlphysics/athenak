"""Analytic-gradient and constrained-transport EMRI wind-tunnel regression."""

from pathlib import Path

import numpy as np
import pytest

import bin_convert
import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
INPUT_FILE = str(ROOT / "inputs" / "emri" / "emri_windtunnel_smoke.athinput")
MAX_DIVB = 1.0e-12


def _run(run_dir: Path, basename: str, gradients: bool) -> tuple[dict, dict]:
    flags = [
        "-d",
        str(run_dir),
        f"job/basename={basename}",
        "time/nlim=0",
        "output1/variable=mhd_w_bcc",
        "output2/variable=mhd_divb",
        "output3/dt=0",
        "output4/dt=0",
    ]
    if gradients:
        flags += [
            "problem/dlnrho_dxh1=0.02",
            "problem/dlnpgas_dxh3=-0.01",
            "problem/du2_dxh1=-0.015",
            "problem/db1_dxh1=0.001",
            "problem/db1_dxh2=0.0003",
            "problem/db2_dxh2=-0.0004",
            "problem/db3_dxh3=-0.0006",
        ]
    assert testutils.run(INPUT_FILE, flags)
    state_path = min((run_dir / "bin").glob(f"{basename}.mhd_w_bcc.*.bin"))
    divb_path = min((run_dir / "bin").glob(f"{basename}.mhd_divb.*.bin"))
    return (
        bin_convert.read_binary(str(state_path)),
        bin_convert.read_binary(str(divb_path)),
    )


def _single_block_field(dump: dict, variable: str) -> np.ndarray:
    values = dump["mb_data"][variable]
    assert len(values) == 1
    return np.asarray(values[0])


def _maximum_divb(dump: dict) -> float:
    return max(float(np.max(np.abs(values))) for values in dump["mb_data"]["divb"])


def test_emri_analytic_gradients_preserve_uniform_limit_and_divb(tmp_path: Path) -> None:
    zero_state, zero_divb = _run(tmp_path / "zero", "emri_zero", gradients=False)
    gradient_state, gradient_divb = _run(
        tmp_path / "gradient", "emri_gradient", gradients=True
    )

    zero_density = _single_block_field(zero_state, "dens")
    zero_pressure = _single_block_field(zero_state, "press")
    assert float(zero_density[0, 0, 0]) == pytest.approx(1.0, abs=1.0e-7)
    assert float(zero_density[-1, -1, -1]) == pytest.approx(1.0, abs=1.0e-7)
    assert float(zero_pressure[0, 0, 0]) == pytest.approx(1.0e-2, rel=1.0e-6)
    assert float(zero_pressure[-1, -1, -1]) == pytest.approx(1.0e-2, rel=1.0e-6)

    density = _single_block_field(gradient_state, "dens")
    pressure = _single_block_field(gradient_state, "press")
    tangential_velocity = _single_block_field(gradient_state, "vely")
    assert float(density[0, 0, -1]) > float(density[0, 0, 0])
    assert float(pressure[-1, 0, 0]) < float(pressure[0, 0, 0])
    assert float(tangential_velocity[0, 0, -1]) < float(
        tangential_velocity[0, 0, 0]
    )

    assert _maximum_divb(zero_divb) < MAX_DIVB
    assert _maximum_divb(gradient_divb) < MAX_DIVB


def test_emri_user_boundary_is_active_level_safe(tmp_path: Path) -> None:
    flags = [
        "-d",
        str(tmp_path / "subcycle"),
        "job/basename=emri_user_bc_subcycle",
        "mesh/nx1=16",
        "mesh/nx2=16",
        "mesh/nx3=16",
        "meshblock/nx1=8",
        "meshblock/nx2=8",
        "meshblock/nx3=8",
        "mesh_refinement/refinement=adaptive",
        "mesh_refinement/num_levels=2",
        "mesh_refinement/max_nmb_per_rank=128",
        "time/subcycling=level",
        "time/nlim=2",
        "time/tlim=0.05",
        "adm/dynamic=true",
        "problem/user_hist=false",
        "problem/refinement_radius_level_1=4",
        "output1/dt=0",
        "output2/dt=0",
        "output3/dt=0",
        "output4/dt=0",
        "output5/dt=0",
    ]
    assert testutils.run(INPUT_FILE, flags)
