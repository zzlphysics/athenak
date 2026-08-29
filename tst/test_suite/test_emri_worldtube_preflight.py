"""Tests for Kerr circular-frame generation and metadata-only coverage audits."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import build_kerr_circular_frame as circular  # noqa: E402
import extract_global_worldtube as extract  # noqa: E402
import extract_static_taylor_worldtube as static  # noqa: E402
import preflight_global_worldtube as preflight  # noqa: E402


def _descriptor(
    time: float,
    storage: str = "dense_uniform",
    logical: tuple[tuple[int, int, int], ...] = ((0, 0, 0),),
    block_shape: tuple[int, int, int] = (8, 8, 8),
) -> extract.SnapshotDescriptor:
    return extract.SnapshotDescriptor(
        time=time,
        cycle=round(10 * time),
        lower=np.asarray((-2.0, -2.0, -2.0)),
        spacing=np.asarray((0.5, 0.5, 0.5)),
        shape_xyz=(8, 8, 8),
        state_path=Path(f"state-{time}.bin"),
        adm_path=Path(f"adm-{time}.bin"),
        source_level=0 if storage == "dense_uniform" else 1,
        source_meshblock_count=len(logical),
        available_leaf_levels=(0,) if storage == "dense_uniform" else (0, 1),
        source_storage=storage,
        block_shape_xyz=block_shape,
        block_logical=np.asarray(logical),
    )


def test_circular_kerr_frequency_and_isco_match_schwarzschild_limits() -> None:
    omega = circular.circular_kerr_omega(1.0, 0.0, 10.0, 1)
    np.testing.assert_allclose(omega, 10.0**-1.5, rtol=0.0, atol=2.0e-17)
    np.testing.assert_allclose(circular.kerr_isco(1.0, 0.0, 1), 6.0, atol=0.0)
    np.testing.assert_allclose(circular.kerr_isco(1.0, 0.0, -1), 6.0, atol=0.0)


def test_generated_kerr_frame_is_orthonormal_at_knots() -> None:
    document = circular.build_frame_document(
        (0.0, 0.5, 1.0),
        primary_mass=1.0,
        dimensionless_spin=0.7,
        orbital_radius=8.0,
    )
    generator = document["generator"]
    mass = float(generator["primary_mass"])
    spin = float(generator["primary_spin_length"])
    expected = np.diag((-1.0, 1.0, 1.0, 1.0))
    for worldline, tangent, legs in zip(
        document["worldline"],
        document["worldline_tangent"],
        document["spatial_legs"],
        strict=True,
    ):
        position = np.asarray(worldline)[1:]
        metric = circular.kerr_schild_metric(position, mass, spin)
        tangent_array = np.asarray(tangent)
        norm2 = static.metric_inner(metric, tangent_array, tangent_array)
        tetrad = np.vstack(
            (tangent_array / np.sqrt(-norm2), np.asarray(legs).T)
        )
        np.testing.assert_allclose(
            tetrad @ metric @ tetrad.T, expected, rtol=0.0, atol=3.0e-14
        )
    diagnostics = document["diagnostics"]
    assert diagnostics["maximum_worldline_position_error_over_coordinate_radius"] < 1e-8
    assert diagnostics["maximum_interpolated_tetrad_gram_error"] < 1e-6
    assert diagnostics["minimum_origin_jacobian_absolute_determinant"] > 0.5


def test_circular_frame_rejects_an_orbit_inside_isco() -> None:
    try:
        circular.build_frame_document((0.0, 1.0), 1.0, 0.0, 5.9)
    except ValueError as error:
        assert "inside Kerr ISCO" in str(error)
    else:
        raise AssertionError("an orbit inside Schwarzschild ISCO was accepted")


def test_circular_frame_applies_extreme_mass_ratio_unit_scaling() -> None:
    scale = 1.0e5
    document = circular.build_frame_document(
        (2.0, 2.5, 3.0),
        primary_mass=1.0,
        dimensionless_spin=0.5,
        orbital_radius=8.0,
        global_length_in_local_units=scale,
    )
    np.testing.assert_allclose(document["times"], (2.0e5, 2.5e5, 3.0e5))
    np.testing.assert_allclose(
        np.asarray(document["worldline"])[:, 0], (2.0, 2.5, 3.0)
    )
    np.testing.assert_allclose(
        np.asarray(document["worldline_tangent"])[:, 0], 1.0 / scale
    )
    series = extract.AffineFrameSeries.from_document(document)
    _, instantaneous = series.evaluate(2.25e5)
    assert np.linalg.cond(instantaneous.jacobian()) < 10.0
    diagnostics = document["diagnostics"]
    assert diagnostics["minimum_scaled_origin_jacobian_absolute_determinant"] > 0.5


def test_sparse_preflight_accepts_same_block_and_measures_zero_halo() -> None:
    descriptor = _descriptor(
        0.0,
        storage="sparse_fixed_leaf_level",
        logical=((1, 1, 1),),
        block_shape=(2, 2, 2),
    )
    audit = preflight.SnapshotCoverage(descriptor)
    coordinate = np.asarray((2.2, 2.3, 2.4))
    point = descriptor.lower + (coordinate + 0.5) * descriptor.spacing
    audit.audit(point[None, :])
    assert audit.queried_position_count == 1
    assert audit.minimum_additional_stencil_halo_cells == 0


def test_sparse_preflight_rejects_a_coarse_fine_stencil() -> None:
    descriptor = _descriptor(
        0.0,
        storage="sparse_fixed_leaf_level",
        logical=((1, 1, 1),),
        block_shape=(2, 2, 2),
    )
    audit = preflight.SnapshotCoverage(descriptor)
    coordinate = np.asarray((3.2, 2.3, 2.4))
    point = descriptor.lower + (coordinate + 0.5) * descriptor.spacing
    try:
        audit.audit(point[None, :])
    except ValueError as error:
        assert "level-1 leaf block" in str(error)
        assert "logical_block=(2, 1, 1)" in str(error)
    else:
        raise AssertionError("a coarse-fine preflight stencil was accepted")


def test_coverage_series_checks_both_temporal_endpoints() -> None:
    descriptors = (_descriptor(0.0), _descriptor(1.0))
    series = preflight.CoverageSeries(descriptors)
    series.audit(np.asarray(((0.25, 0.0, 0.0, 0.0),)))
    assert series.snapshots[0].queried_position_count == 1
    assert series.snapshots[1].queried_position_count == 1
