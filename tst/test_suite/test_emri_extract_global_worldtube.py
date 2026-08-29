"""Tests for global GRMHD snapshot to moving local worldtube extraction."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import extract_global_worldtube as extract  # noqa: E402
import worldtube_flux_emf as worldtube  # noqa: E402


def _constant_global_values() -> np.ndarray:
    values = np.zeros(len(extract.ALL_VARIABLES))
    values[:8] = (2.0, 0.3, -0.1, 0.2, 0.4, 0.2, 0.05, -0.1)
    offset = len(extract.STATE_VARIABLES)
    values[offset : offset + 10] = (
        1.0,
        0.0,
        0.0,
        1.0,
        0.0,
        1.0,
        1.0,
        0.0,
        0.0,
        0.0,
    )
    return values


def _snapshot(time: float, values: np.ndarray | None = None) -> extract.UniformSnapshot:
    source = _constant_global_values() if values is None else values
    shape = (12, 12, 12)
    arrays = np.broadcast_to(source[:, None, None, None], (source.size, *shape)).copy()
    return extract.UniformSnapshot(
        time=time,
        cycle=round(10 * time),
        lower=np.asarray((-3.0, -3.0, -3.0)),
        spacing=np.asarray((0.5, 0.5, 0.5)),
        shape_xyz=shape,
        values=arrays,
        state_path=Path(f"state-{time}.bin"),
        adm_path=Path(f"adm-{time}.bin"),
    )


def _descriptor(snapshot: extract.UniformSnapshot) -> extract.SnapshotDescriptor:
    return extract.SnapshotDescriptor(
        time=snapshot.time,
        cycle=snapshot.cycle,
        lower=snapshot.lower,
        spacing=snapshot.spacing,
        shape_xyz=snapshot.shape_xyz,
        state_path=snapshot.state_path,
        adm_path=snapshot.adm_path,
        source_level=snapshot.source_level,
        source_meshblock_count=snapshot.source_meshblock_count,
        available_leaf_levels=snapshot.available_leaf_levels,
        source_storage=snapshot.source_storage,
    )


def _fixed_level_snapshot(
    time: float,
    include_right_block: bool = True,
    include_extra_block: bool = False,
) -> extract.FixedLevelSnapshot:
    logical = [(1, 1, 1)]
    if include_right_block:
        logical.append((2, 1, 1))
    if include_extra_block:
        logical.append((0, 0, 0))
    block_shape = (2, 2, 2)
    values = np.empty((len(logical), len(extract.ALL_VARIABLES), 2, 2, 2))
    background = _constant_global_values()
    for block, location in enumerate(logical):
        for z in range(2):
            for y in range(2):
                for x in range(2):
                    global_index = np.asarray(location) * 2 + (x, y, z)
                    values[block, :, z, y, x] = background
                    values[block, 0, z, y, x] = (
                        2.0
                        + 0.3 * global_index[0]
                        - 0.2 * global_index[1]
                        + 0.1 * global_index[2]
                    )
    return extract.FixedLevelSnapshot(
        time=time,
        cycle=round(10 * time),
        lower=np.asarray((-2.0, -2.0, -2.0)),
        spacing=np.asarray((0.5, 0.5, 0.5)),
        shape_xyz=(8, 8, 8),
        block_shape_xyz=block_shape,
        block_logical=np.asarray(logical),
        values=values,
        state_path=Path(f"state-fixed-{time}.bin"),
        adm_path=Path(f"adm-fixed-{time}.bin"),
        source_level=1,
        available_leaf_levels=(0, 1),
    )


def _frame_series(
    velocity: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> extract.AffineFrameSeries:
    times = np.asarray((0.0, 1.0))
    velocity_array = np.asarray(velocity)
    worldline = np.zeros((2, 4))
    worldline[:, 0] = times
    worldline[:, 1:] = times[:, None] * velocity_array
    tangent = np.broadcast_to(np.asarray((1.0, *velocity)), (2, 4)).copy()
    legs = np.zeros((2, 4, 3))
    legs[:, 1:, :] = np.eye(3)
    return extract.AffineFrameSeries(
        times=times,
        worldline=extract.HermiteSeries(times, worldline, tangent, "worldline"),
        spatial_legs=extract.HermiteSeries(
            times, legs, np.zeros_like(legs), "spatial legs"
        ),
    )


def test_uniform_snapshot_trilinear_interpolation_crosses_cell_boundaries() -> None:
    snapshot = _snapshot(0.0)
    z, y, x = np.meshgrid(
        np.arange(12), np.arange(12), np.arange(12), indexing="ij"
    )
    snapshot.values[0] = 1.2 + 0.3 * x - 0.2 * y + 0.1 * z
    points = np.asarray(
        ((-2.25, -1.75, -1.25), (0.1, 0.4, -0.2), (2.75, 2.75, 2.75))
    )
    sampled = snapshot.sample(points)[:, 0]
    cell_index = (points - snapshot.lower) / snapshot.spacing - 0.5
    expected = (
        1.2
        + 0.3 * cell_index[:, 0]
        - 0.2 * cell_index[:, 1]
        + 0.1 * cell_index[:, 2]
    )
    np.testing.assert_allclose(sampled, expected, rtol=0.0, atol=2.0e-14)


def test_fixed_level_interpolation_crosses_same_level_meshblock_boundary() -> None:
    snapshot = _fixed_level_snapshot(0.0)
    cell_coordinate = np.asarray((3.25, 2.3, 2.6))
    point = snapshot.lower + (cell_coordinate + 0.5) * snapshot.spacing
    sampled = snapshot.sample(point[None, :])[0, 0]
    expected = (
        2.0
        + 0.3 * cell_coordinate[0]
        - 0.2 * cell_coordinate[1]
        + 0.1 * cell_coordinate[2]
    )
    np.testing.assert_allclose(sampled, expected, rtol=0.0, atol=2.0e-14)


def test_fixed_level_interpolation_rejects_coarse_fine_interface() -> None:
    snapshot = _fixed_level_snapshot(0.0, include_right_block=False)
    cell_coordinate = np.asarray((3.25, 2.3, 2.6))
    point = snapshot.lower + (cell_coordinate + 0.5) * snapshot.spacing
    try:
        snapshot.sample(point[None, :])
    except ValueError as error:
        message = str(error)
        assert "level-1 leaf MeshBlock" in message
        assert "logical_block=(2, 1, 1)" in message
    else:
        raise AssertionError("coarse-fine interpolation stencil was accepted")


def test_snapshot_series_allows_fixed_level_topology_change_with_coverage() -> None:
    left = _fixed_level_snapshot(0.0)
    right = _fixed_level_snapshot(1.0, include_extra_block=True)
    series = extract.SnapshotSeries((left, right))
    cell_coordinate = np.asarray((3.25, 2.3, 2.6))
    point = left.lower + (cell_coordinate + 0.5) * left.spacing
    event = np.asarray(((0.4, *point),))
    expected = left.sample(point[None, :])
    np.testing.assert_allclose(series.sample(event), expected, atol=2.0e-14)


def test_unit_scaling_preserves_metric_and_rescales_grmhd_state() -> None:
    values = _constant_global_values()
    jacobian = 0.1 * np.eye(4)
    state, metric, four_velocity, faraday = extract.transform_grmhd_sample(
        values,
        jacobian,
        global_length_in_local_units=10.0,
        density_renormalization=4.0,
    )
    np.testing.assert_allclose(metric, np.diag((-1.0, 1.0, 1.0, 1.0)), atol=2.0e-15)
    np.testing.assert_allclose(state[:5], (0.08, 0.3, -0.1, 0.2, 0.016), atol=2.0e-15)
    np.testing.assert_allclose(state[5:], 0.2 * values[5:8], atol=2.0e-15)
    np.testing.assert_allclose(four_velocity @ metric @ four_velocity, -1.0, atol=2.0e-15)
    np.testing.assert_allclose(faraday @ four_velocity, 0.0, atol=2.0e-15)


def test_identity_map_recovers_curved_adm_primitives() -> None:
    values = _constant_global_values()
    offset = len(extract.STATE_VARIABLES)
    values[offset : offset + 10] = (
        1.3,
        0.08,
        -0.04,
        1.1,
        0.06,
        0.9,
        0.82,
        0.05,
        -0.03,
        0.02,
    )
    state, metric, four_velocity, faraday = extract.transform_grmhd_sample(
        values, np.eye(4)
    )
    gamma, alpha, beta = extract._global_geometry(values)
    expected_metric = extract.static.spacetime_metric_from_adm(
        gamma, alpha, beta
    )
    np.testing.assert_allclose(state, values[:8], rtol=0.0, atol=2.0e-15)
    np.testing.assert_allclose(metric, expected_metric, rtol=0.0, atol=2.0e-15)
    np.testing.assert_allclose(
        four_velocity @ metric @ four_velocity, -1.0, atol=2.0e-15
    )
    np.testing.assert_allclose(faraday @ four_velocity, 0.0, atol=2.0e-15)


def test_constant_global_state_builds_exact_static_ct_worldtube() -> None:
    snapshots = extract.SnapshotSeries((_snapshot(0.0), _snapshot(1.0)))
    sampler = extract.GlobalWorldtubeSampler(
        snapshots, _frame_series(), 1.0, 1.0
    )
    geometry = extract.CubeGeometry(np.zeros(3), 1.0, 4)
    times = np.asarray((0.0, 0.4, 1.0))
    faces, diagnostics = extract.sample_worldtube(
        sampler, geometry, times, quadrature_order=2
    )
    validation = worldtube.validate_worldtube(times, faces)
    source = _constant_global_values()
    lorentz = np.sqrt(1.0 + np.dot(source[1:4], source[1:4]))
    electric = -np.cross(source[1:4] / lorentz, source[5:8])
    area = geometry.spacing**2
    length = geometry.spacing
    for name in worldtube.FACE_NAMES:
        orientation = worldtube.ORIENTATIONS[name]
        expected_state = np.broadcast_to(
            source[None, :8, None, None], faces[name].cell_state.shape
        )
        np.testing.assert_allclose(
            faces[name].cell_state,
            expected_state,
            rtol=0.0,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            faces[name].normal_flux,
            orientation.normal_sign * source[5 + orientation.normal_axis] * area,
            rtol=0.0,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            faces[name].emf_u,
            orientation.u_sign * electric[orientation.u_axis] * length,
            rtol=0.0,
            atol=2.0e-14,
        )
        np.testing.assert_allclose(
            faces[name].emf_v,
            orientation.v_sign * electric[orientation.v_axis] * length,
            rtol=0.0,
            atol=2.0e-14,
        )
    assert validation["maximum_shared_edge_emf_residual"] == 0.0
    assert validation["maximum_closed_surface_flux"] < 2.0e-16
    assert diagnostics["raw_sampling"]["maximum_ideal_mhd_residual"] < 2.0e-16


def test_lazy_snapshot_series_matches_eager_and_loads_each_snapshot_once() -> None:
    source_snapshots = tuple(_snapshot(time) for time in (0.0, 0.5, 1.0))
    loader_calls = []

    def loader(index: int) -> extract.UniformSnapshot:
        loader_calls.append(index)
        return source_snapshots[index]

    lazy = extract.LazySnapshotSeries(
        tuple(_descriptor(snapshot) for snapshot in source_snapshots),
        loader,
        cache_size=2,
    )
    eager_sampler = extract.GlobalWorldtubeSampler(
        extract.SnapshotSeries(source_snapshots), _frame_series(), 1.0, 1.0
    )
    lazy_sampler = extract.GlobalWorldtubeSampler(
        lazy, _frame_series(), 1.0, 1.0
    )
    geometry = extract.CubeGeometry(np.zeros(3), 0.5, 2)
    times = np.asarray((0.0, 0.5, 1.0))
    eager_faces, _ = extract.sample_worldtube(
        eager_sampler, geometry, times, quadrature_order=1
    )
    lazy_faces, lazy_diagnostics = extract.sample_worldtube(
        lazy_sampler, geometry, times, quadrature_order=1
    )
    for name in worldtube.FACE_NAMES:
        np.testing.assert_array_equal(
            lazy_faces[name].cell_state, eager_faces[name].cell_state
        )
        np.testing.assert_array_equal(
            lazy_faces[name].normal_flux, eager_faces[name].normal_flux
        )
        np.testing.assert_array_equal(lazy_faces[name].emf_u, eager_faces[name].emf_u)
        np.testing.assert_array_equal(lazy_faces[name].emf_v, eager_faces[name].emf_v)
    loading = lazy_diagnostics["snapshot_loading"]
    assert loader_calls == [0, 1, 2]
    assert loading["cache_misses"] == 3
    assert loading["cache_evictions"] == 1
    assert loading["peak_cached_snapshots"] == 2
    assert loading["resident_snapshot_count"] == 2
    assert loading["cache_hits"] > 0


def test_moving_frame_includes_motional_emf_and_eulerian_velocity() -> None:
    source = _constant_global_values()
    frame_velocity = np.asarray((0.07, -0.04, 0.02))
    snapshots = extract.SnapshotSeries((_snapshot(0.0), _snapshot(1.0)))
    sampler = extract.GlobalWorldtubeSampler(
        snapshots, _frame_series(tuple(frame_velocity)), 1.0, 1.0
    )
    state, faraday = sampler.sample(0.5, np.asarray(((0.0, 0.0, 0.0),)))
    lorentz = np.sqrt(1.0 + np.dot(source[1:4], source[1:4]))
    global_electric = -np.cross(source[1:4] / lorentz, source[5:8])
    expected_electric = global_electric + np.cross(frame_velocity, source[5:8])
    local_electric = -faraday[0, 0, 1:]
    # The affine map preserves the t=constant slicing, so its shift absorbs the
    # coordinate motion while the Eulerian spatial four-velocity is unchanged.
    expected_u = source[1:4]
    np.testing.assert_allclose(local_electric, expected_electric, atol=2.0e-15)
    np.testing.assert_allclose(state[0, 1:4], expected_u, atol=2.0e-15)
    np.testing.assert_allclose(state[0, 5:8], source[5:8], atol=2.0e-15)
