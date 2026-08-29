"""Time-dependent Taylor-profile replay regression for the EMRI wind tunnel."""

from pathlib import Path
import math

import numpy as np

import bin_convert
import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
INPUT_FILE = str(ROOT / "inputs" / "emri" / "emri_windtunnel_smoke.athinput")
MAX_DIVB = 1.0e-12


def _write_profile(path: Path, upper_b3: float = 0.03) -> None:
    lower_gradient = [0.001, 0.0003, 0.0, 0.0, -0.0004, 0.0, 0.0, 0.0, -0.0006]
    upper_gradient = [0.002, 0.0006, 0.0, 0.0, -0.0008, 0.0, 0.0, 0.0, -0.0012]
    lower = (
        [1.0, 0.01, 0.5, 0.0, 0.0, 0.0, 0.0, 0.01]
        + [0.0] * 15
        + lower_gradient
    )
    upper = (
        [4.0, 0.04, 0.7, 0.0, 0.0, 0.0, 0.0, upper_b3]
        + [0.0] * 15
        + upper_gradient
    )
    assert len(lower) == len(upper) == 32
    primary_mass = 1000.0
    primary_spin = 500.0
    orbital_radius = 10000.0
    coordinate_radius = math.sqrt(orbital_radius**2 + primary_spin**2)
    orbital_omega = math.sqrt(primary_mass) / (
        orbital_radius**1.5 + primary_spin * math.sqrt(primary_mass)
    )
    lines = [
        "# athenak-emri-taylor-series-v2",
        f"# source_coordinate_radius_local_units: {coordinate_radius:.17g}",
        "# source_coordinate_angular_frequency_local_units: "
        f"{orbital_omega:.17g}",
        "# orbit_tolerance: 1e-12",
        "0 " + " ".join(str(value) for value in lower),
        "0.1 " + " ".join(str(value) for value in upper),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _run(
    run_dir: Path,
    profile: Path,
    flags: list[str],
    restart_output: bool = False,
) -> None:
    common = [
        "-d",
        str(run_dir),
        "adm/dynamic=true",
        f"problem/profile_file={profile}",
        "output1/variable=mhd_w_bcc",
        "output2/variable=mhd_divb",
        f"output3/dt={0.05 if restart_output else 0.0}",
        "output4/dt=0",
    ]
    assert testutils.run(INPUT_FILE, common + flags)


def _maximum_divb(path: Path) -> float:
    dump = bin_convert.read_binary(str(path))
    return max(
        float(np.max(np.abs(np.asarray(values))))
        for values in dump["mb_data"]["divb"]
    )


def test_emri_profile_interpolation_and_ct_replay(tmp_path: Path) -> None:
    profile = tmp_path / "profile.dat"
    _write_profile(profile)

    midpoint = tmp_path / "midpoint"
    _run(
        midpoint,
        profile,
        [
            "job/basename=midpoint",
            "time/nlim=0",
            "problem/profile_time_offset=0.05",
        ],
    )
    state_path = min((midpoint / "bin").glob("midpoint.mhd_w_bcc.*.bin"))
    divb_path = min((midpoint / "bin").glob("midpoint.mhd_divb.*.bin"))
    state = bin_convert.read_binary(str(state_path))["mb_data"]
    density = np.asarray(state["dens"][0])
    pressure = np.asarray(state["press"][0])
    non_excision = density > 0.5
    np.testing.assert_allclose(density[non_excision], 2.0, rtol=0.0, atol=1.0e-7)
    np.testing.assert_allclose(pressure[non_excision], 0.02, rtol=1.0e-6)
    assert _maximum_divb(divb_path) < MAX_DIVB

    evolution = tmp_path / "evolution"
    _run(
        evolution,
        profile,
        [
            "job/basename=evolution",
            "time/tlim=0.1",
            "time/nlim=2",
            "problem/profile_time_offset=0.0",
        ],
    )
    divb_outputs = sorted((evolution / "bin").glob("evolution.mhd_divb.*.bin"))
    assert len(divb_outputs) == 2
    assert max(_maximum_divb(path) for path in divb_outputs) < MAX_DIVB

    restart_source = tmp_path / "restart_source"
    _run(
        restart_source,
        profile,
        [
            "job/basename=profile_restart",
            "time/tlim=0.05",
            "time/nlim=1",
            "problem/profile_time_offset=0.0",
        ],
        restart_output=True,
    )
    checkpoint = max((restart_source / "rst").glob("profile_restart.*.rst"))
    resumed = tmp_path / "resumed"
    assert testutils.run_command(
        [
            "./athena",
            "-r",
            str(checkpoint),
            "-d",
            str(resumed),
            "time/tlim=0.1",
            "time/nlim=2",
            "output3/dt=0",
            "output4/dt=0",
        ]
    )

    _write_profile(profile, upper_b3=0.031)
    assert not testutils.run_command(
        [
            "./athena",
            "-r",
            str(checkpoint),
            "-d",
            str(tmp_path / "tampered"),
            "time/tlim=0.1",
            "time/nlim=2",
            "output3/dt=0",
            "output4/dt=0",
        ]
    )
