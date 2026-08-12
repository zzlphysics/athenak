"""Native metric-aware GRMHD output regression."""

from pathlib import Path

import numpy as np

import bin_convert
import test_suite.testutils as testutils


INPUT_FILE = "../../../inputs/dyngr/effective_bbh_amr_subcycle_smoke.athinput"
BASENAME = "effective_bbh_gr_diagnostics"


def _latest_dump(run_dir: Path, variable: str):
    paths = sorted((run_dir / "bin").glob(f"{BASENAME}.{variable}.*.bin"))
    assert paths, f"No dumps found for {variable}"
    dumps = [bin_convert.read_binary(str(path)) for path in paths]
    return max(dumps, key=lambda dump: (dump["cycle"], dump["time"]))


def _block_map(dump):
    return {tuple(int(value) for value in row): index
            for index, row in enumerate(dump["mb_logical"])}


def _field(dump, name: str, block: int) -> np.ndarray:
    return np.asarray(dump["mb_data"][name][block], dtype=np.float64)


def _reconstruct(fluid, metric, fluid_block: int, metric_block: int):
    gxx = _field(metric, "adm_gxx", metric_block)
    gxy = _field(metric, "adm_gxy", metric_block)
    gxz = _field(metric, "adm_gxz", metric_block)
    gyy = _field(metric, "adm_gyy", metric_block)
    gyz = _field(metric, "adm_gyz", metric_block)
    gzz = _field(metric, "adm_gzz", metric_block)
    determinant = (
        gxx * (gyy * gzz - gyz * gyz)
        - gxy * (gxy * gzz - gxz * gyz)
        + gxz * (gxy * gyz - gxz * gyy)
    )

    wx = _field(fluid, "velx", fluid_block)
    wy = _field(fluid, "vely", fluid_block)
    wz = _field(fluid, "velz", fluid_block)
    wx_lower = gxx * wx + gxy * wy + gxz * wz
    wy_lower = gxy * wx + gyy * wy + gyz * wz
    wz_lower = gxz * wx + gyz * wy + gzz * wz
    lorentz_squared = 1.0 + wx * wx_lower + wy * wy_lower + wz * wz_lower
    lorentz = np.sqrt(lorentz_squared)

    inverse_sqrt_determinant = 1.0 / np.sqrt(determinant)
    bx = _field(fluid, "bcc1", fluid_block) * inverse_sqrt_determinant
    by = _field(fluid, "bcc2", fluid_block) * inverse_sqrt_determinant
    bz = _field(fluid, "bcc3", fluid_block) * inverse_sqrt_determinant
    bx_lower = gxx * bx + gxy * by + gxz * bz
    by_lower = gxy * bx + gyy * by + gyz * bz
    bz_lower = gxz * bx + gyz * by + gzz * bz
    magnetic_squared = bx * bx_lower + by * by_lower + bz * bz_lower
    magnetic_velocity = (
        bx * wx_lower + by * wy_lower + bz * wz_lower
    ) / lorentz
    bsq = magnetic_velocity**2 + magnetic_squared / lorentz_squared
    density = _field(fluid, "dens", fluid_block)
    pressure = _field(fluid, "press", fluid_block)
    return {
        "gr_bsq": bsq,
        "gr_lorentz": lorentz,
        "gr_sigma": bsq / density,
        "gr_beta_inv": 0.5 * bsq / pressure,
    }


def test_native_grmhd_diagnostics_match_synchronized_metric_reconstruction(tmp_path):
    """The native output must undensitize bcc and use the ADM spatial metric."""
    flags = [
        "-d",
        str(tmp_path),
        f"job/basename={BASENAME}",
        "time/nlim=1",
        "output2/file_type=bin",
        "output2/variable=adm",
        "output2/dcycle=1",
        "output3/file_type=bin",
        "output3/variable=mhd_gr_diagnostics",
        "output3/dcycle=1",
    ]
    assert testutils.mpi_run(INPUT_FILE, flags, threads=1)

    fluid = _latest_dump(tmp_path, "mhd_w_bcc")
    metric = _latest_dump(tmp_path, "adm")
    diagnostics = _latest_dump(tmp_path, "mhd_gr_diagnostics")
    assert fluid["cycle"] == metric["cycle"] == diagnostics["cycle"] == 1
    assert fluid["time"] == metric["time"] == diagnostics["time"]
    assert diagnostics["var_names"] == [
        "gr_bsq",
        "gr_lorentz",
        "gr_sigma",
        "gr_beta_inv",
        "gr_excision_mask",
    ]

    fluid_blocks = _block_map(fluid)
    metric_blocks = _block_map(metric)
    diagnostic_blocks = _block_map(diagnostics)
    assert fluid_blocks.keys() == metric_blocks.keys() == diagnostic_blocks.keys()
    mask_values = set()
    for logical_location in fluid_blocks:
        expected = _reconstruct(
            fluid,
            metric,
            fluid_blocks[logical_location],
            metric_blocks[logical_location],
        )
        diagnostic_block = diagnostic_blocks[logical_location]
        for name, reference in expected.items():
            measured = _field(diagnostics, name, diagnostic_block)
            assert np.isfinite(measured).all()
            np.testing.assert_allclose(measured, reference, rtol=5.0e-6, atol=1.0e-20)
        assert np.min(_field(diagnostics, "gr_bsq", diagnostic_block)) >= 0.0
        assert np.min(_field(diagnostics, "gr_lorentz", diagnostic_block)) >= 1.0
        mask = _field(diagnostics, "gr_excision_mask", diagnostic_block)
        assert np.isfinite(mask).all()
        mask_values.update(float(value) for value in np.unique(mask))
    assert mask_values.issubset({0.0, 1.0})
    assert mask_values == {0.0, 1.0}
