"""Unit tests for the static global-to-local EMRI Taylor extractor."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import extract_static_taylor_worldtube as extractor  # noqa: E402


def _flat_adm(count: int) -> dict[str, np.ndarray]:
    zeros = np.zeros(count)
    ones = np.ones(count)
    return {
        "adm_gxx": ones,
        "adm_gxy": zeros,
        "adm_gxz": zeros,
        "adm_gyy": ones,
        "adm_gyz": zeros,
        "adm_gzz": ones,
        "adm_alpha": ones,
        "adm_betax": zeros,
        "adm_betay": zeros,
        "adm_betaz": zeros,
    }


def test_source_tetrad_is_orthonormal_and_dual() -> None:
    gamma = np.eye(3)
    metric = extractor.spacetime_metric_from_adm(gamma, 1.0, np.zeros(3))
    basis = np.eye(3)
    tetrad, coframe = extractor.build_source_tetrad(
        metric, np.asarray((0.0, 0.2, 0.0)), basis
    )
    minkowski = np.diag((-1.0, 1.0, 1.0, 1.0))
    np.testing.assert_allclose(tetrad @ metric @ tetrad.T, minkowski, atol=1.0e-14)
    np.testing.assert_allclose(coframe @ tetrad.T, np.eye(4), atol=1.0e-14)


def test_static_profile_recovers_manufactured_grmhd_state() -> None:
    rng = np.random.default_rng(5824)
    position = rng.uniform(-1.0, 1.0, size=(256, 3))
    position = position[np.sum(position**2, axis=1) < 1.0]
    count = position.shape[0]
    log_density_gradient = np.asarray((0.07, -0.03, 0.02))
    log_pressure_gradient = np.asarray((-0.01, 0.04, -0.02))
    velocity0 = np.asarray((0.2, -0.1, 0.05))
    velocity_gradient = np.asarray(
        ((0.01, 0.02, -0.01), (-0.03, 0.01, 0.02), (0.0, -0.02, 0.015))
    )
    magnetic0 = np.asarray((0.4, -0.2, 0.1))
    magnetic_gradient = np.asarray(
        ((0.03, -0.02, 0.01), (0.04, -0.01, -0.03), (0.02, 0.05, -0.02))
    )
    density0 = 1.3
    pressure0 = 0.08
    primitive_velocity = velocity0 + position @ velocity_gradient.T
    magnetic = magnetic0 + position @ magnetic_gradient.T
    primitive = {
        "dens": density0 * np.exp(position @ log_density_gradient),
        "press": pressure0 * np.exp(position @ log_pressure_gradient),
        "velx": primitive_velocity[:, 0],
        "vely": primitive_velocity[:, 1],
        "velz": primitive_velocity[:, 2],
        "bcc1": magnetic[:, 0],
        "bcc2": magnetic[:, 1],
        "bcc3": magnetic[:, 2],
    }
    cloud = extractor.SampleCloud(
        global_position=position.copy(),
        local_position=position.copy(),
        cell_volume=rng.uniform(0.4, 1.6, size=count),
        primitive=primitive,
        adm=_flat_adm(count),
    )
    parameters, diagnostics = extractor.fit_static_profile(cloud, np.eye(4), 1.0)

    np.testing.assert_allclose(parameters["rho0"], density0, atol=2.0e-14)
    np.testing.assert_allclose(parameters["pgas0"], pressure0, atol=2.0e-14)
    for direction in range(3):
        np.testing.assert_allclose(
            parameters[f"dlnrho_dxh{direction + 1}"],
            log_density_gradient[direction],
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            parameters[f"dlnpgas_dxh{direction + 1}"],
            log_pressure_gradient[direction],
            atol=2.0e-14,
        )
    for component in range(3):
        np.testing.assert_allclose(
            parameters[f"u{component + 1}"], velocity0[component], atol=2.0e-14
        )
        np.testing.assert_allclose(
            parameters[f"b{component + 1}"], magnetic0[component], atol=2.0e-14
        )
        for direction in range(3):
            np.testing.assert_allclose(
                parameters[f"du{component + 1}_dxh{direction + 1}"],
                velocity_gradient[component, direction],
                atol=2.0e-14,
            )
            np.testing.assert_allclose(
                parameters[f"db{component + 1}_dxh{direction + 1}"],
                magnetic_gradient[component, direction],
                atol=2.0e-14,
            )
    assert abs(diagnostics["magnetic_gradient_trace"]) < 1.0e-15
    assert diagnostics["density_log_weighted_rms"] < 1.0e-14
    assert diagnostics["velocity_weighted_rms"] < 1.0e-14
    assert diagnostics["magnetic_weighted_rms"] < 1.0e-14


def test_global_to_local_unit_rescaling() -> None:
    parameters = {name: 0.0 for name in extractor.PROFILE_PARAMETER_ORDER}
    parameters.update(
        rho0=8.0,
        pgas0=4.0,
        u1=0.25,
        b1=6.0,
        dlnrho_dxh1=10.0,
        dlnpgas_dxh2=8.0,
        du2_dxh3=12.0,
        db1_dxh2=14.0,
    )
    scaled = extractor.rescale_profile_parameters(parameters, 2.0, 9.0)
    np.testing.assert_allclose(scaled["rho0"], 18.0)
    np.testing.assert_allclose(scaled["pgas0"], 9.0)
    np.testing.assert_allclose(scaled["u1"], 0.25)
    np.testing.assert_allclose(scaled["b1"], 9.0)
    np.testing.assert_allclose(scaled["dlnrho_dxh1"], 5.0)
    np.testing.assert_allclose(scaled["dlnpgas_dxh2"], 4.0)
    np.testing.assert_allclose(scaled["du2_dxh3"], 6.0)
    np.testing.assert_allclose(scaled["db1_dxh2"], 10.5)
