"""Outer-to-inner CT plus HLL fluid-worldtube replay regression."""

from pathlib import Path
import sys

import numpy as np

import bin_convert
import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import worldtube_flux_emf as worldtube  # noqa: E402


INPUT_FILE = str(EMRI_INPUTS / "emri_windtunnel_smoke.athinput")
MAX_DIVB = 1.0e-11


def _disable_ordinary_outputs() -> list[str]:
    return ["output1/dt=0", "output2/dt=0", "output3/dt=0", "output4/dt=0"]


def test_emri_outer_to_inner_hll_fluid_replay(tmp_path: Path) -> None:
    outer = tmp_path / "outer"
    stream_stem = outer / "tube"
    outer_flags = [
        "-d",
        str(outer),
        "job/basename=outer",
        "time/nlim=1",
        "time/tlim=0.02",
        "emri_worldtube/enabled=true",
        "emri_worldtube/mode=outer",
        "emri_worldtube/overwrite=true",
        f"emri_worldtube/file_basename={stream_stem}",
    ] + _disable_ordinary_outputs()
    assert testutils.run(INPUT_FILE, outer_flags)

    manifest = next(outer.glob("tube.cycle*.manifest.json"))
    times, faces, metadata = worldtube.read_outer_stream(manifest)
    assert metadata["state_variables"][-3:] == ["bcc1", "bcc2", "bcc3"]
    assert faces["x1m"].cell_state.shape[1] == 8
    packed = tmp_path / "tube.npz"
    binary = tmp_path / "tube.inner.bin"
    worldtube.write_worldtube(packed, times, faces, metadata)
    worldtube.write_inner_binary(binary, times, faces, metadata)

    inner = tmp_path / "inner"
    inner_flags = [
        "-d",
        str(inner),
        "job/basename=inner",
        "mesh/nx1=8",
        "mesh/nx2=8",
        "mesh/nx3=8",
        "mesh/x1min=-4",
        "mesh/x1max=4",
        "mesh/x2min=-4",
        "mesh/x2max=4",
        "mesh/x3min=-4",
        "mesh/x3max=4",
        "meshblock/nx1=4",
        "meshblock/nx2=4",
        "meshblock/nx3=4",
        "problem/force_surface_radius=2.5",
        "problem/force_outer_radius_1=2.8",
        "problem/force_outer_radius_2=3.2",
        "problem/force_outer_radius_3=3.7",
        "emri_worldtube/enabled=true",
        "emri_worldtube/mode=inner",
        f"emri_worldtube/file={binary}",
        "emri_worldtube/fluid_boundary=riemann",
        "time/integrator=rk3",
        "time/cfl_number=0.005",
        "time/nlim=20",
        "time/tlim=0.02",
        "output1/variable=mhd_w_bcc",
        "output1/dt=0.02",
        "output2/variable=mhd_divb",
        "output2/dt=0.02",
        "output3/dt=0",
        "output4/dt=0",
    ]
    assert testutils.run(INPUT_FILE, inner_flags)

    state_path = max((inner / "bin").glob("inner.mhd_w_bcc.*.bin"))
    divb_path = max((inner / "bin").glob("inner.mhd_divb.*.bin"))
    state = bin_convert.read_binary(str(state_path))["mb_data"]
    density = np.concatenate([np.asarray(values).ravel() for values in state["dens"]])
    assert np.isfinite(density).all()
    assert float(np.min(density)) > 0.0
    divb = bin_convert.read_binary(str(divb_path))["mb_data"]["divb"]
    maximum_divb = max(float(np.max(np.abs(values))) for values in divb)
    assert maximum_divb < MAX_DIVB
