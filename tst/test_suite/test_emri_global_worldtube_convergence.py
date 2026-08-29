"""Tests for offline-versus-online worldtube convergence analysis."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import analyze_global_worldtube_convergence as convergence  # noqa: E402
import run_global_worldtube_convergence as campaign  # noqa: E402
import worldtube_flux_emf as worldtube  # noqa: E402
import worldtube_frame  # noqa: E402


def _reference(cells: int, times: np.ndarray) -> dict[str, worldtube.FaceData]:
    generator = np.random.default_rng(6712)
    complex_ = worldtube_frame.CubeSurfaceComplex(cells)
    faces = {}
    for name in worldtube.FACE_NAMES:
        faces[name] = worldtube.FaceData(
            cell_state=generator.normal(size=(times.size, 3, cells, cells)),
            normal_flux=np.zeros((times.size, cells, cells)),
            emf_u=np.zeros((times.size - 1, cells + 1, cells)),
            emf_v=np.zeros((times.size - 1, cells, cells + 1)),
        )
    for interval, dt in enumerate(np.diff(times)):
        edge = generator.normal(size=complex_.edge_count)
        emf_u, emf_v = complex_.unpack_edges(edge)
        for name in worldtube.FACE_NAMES:
            faces[name].emf_u[interval] = emf_u[name]
            faces[name].emf_v[interval] = emf_v[name]
            faces[name].normal_flux[interval + 1] = worldtube.faraday_update(
                faces[name].normal_flux[interval],
                emf_u[name],
                emf_v[name],
                float(dt),
            )
    return faces


def test_matching_and_identical_coarsened_reference_have_zero_error() -> None:
    times = np.asarray((0.0, 0.1, 0.3, 0.55, 1.0))
    faces = _reference(3, times)
    selected_times, selected = worldtube.coarsen_worldtube_time(
        times, faces, (0, 2, 4)
    )
    report = convergence.compare_worldtubes(
        times,
        faces,
        selected_times,
        selected,
        state_variables=("a", "b", "c"),
    )
    assert report["reference_endpoint_indices"] == [0, 2, 4]
    assert report["state"]["absolute_linf"] == 0.0
    assert report["normal_flux"]["absolute_linf"] == 0.0
    assert report["emf"]["absolute_linf"] == 0.0


def test_comparison_reports_controlled_emf_error() -> None:
    times = np.asarray((0.0, 0.2, 0.5))
    faces = _reference(2, times)
    selected_times, selected = worldtube.coarsen_worldtube_time(
        times, faces, (0, 2)
    )
    complex_ = worldtube_frame.CubeSurfaceComplex(2)
    curl_free = np.zeros(complex_.edge_count)
    for key, edge_id in complex_._edge_ids.items():
        curl_free[edge_id] = 0.125 if key[0] == 0 else 0.0
    delta_u, delta_v = complex_.unpack_edges(curl_free)
    for name in worldtube.FACE_NAMES:
        selected[name].emf_u[0] += delta_u[name]
        selected[name].emf_v[0] += delta_v[name]
    report = convergence.compare_worldtubes(
        times, faces, selected_times, selected
    )
    assert report["emf"]["absolute_linf"] == 0.125
    assert report["state"]["absolute_linf"] == 0.0


def test_campaign_stride_and_identity_frame_contract() -> None:
    assert campaign.selected_indices(6, 1) == [0, 1, 2, 3, 4, 5]
    assert campaign.selected_indices(6, 2) == [0, 2, 4, 5]
    assert campaign.selected_indices(6, 4) == [0, 4, 5]
    document = campaign.identity_frame_document([2.0, 2.5, 3.0])
    frames = campaign.extract.AffineFrameSeries.from_document(document)
    worldline, instantaneous = frames.evaluate(2.25)
    np.testing.assert_allclose(worldline, (2.25, 0.0, 0.0, 0.0))
    np.testing.assert_allclose(instantaneous.jacobian(), np.eye(4))


def test_campaign_assessment_recovers_second_order_series() -> None:
    cases = {}
    for cells, spacing, flux, emf in (
        (16, 1.0, 1.0e-4, 8.0e-3),
        (32, 0.5, 2.5e-5, 2.0e-3),
    ):
        cases[f"n{cells}_s1_q3"] = {
            "source_spacing": spacing,
            "extraction_wall_seconds": 1.0,
            "comparison": {
                "state": {"relative_l2": 2.0e-8},
                "normal_flux": {"relative_l2": flux},
                "emf": {"relative_l2": emf},
            },
            "projection": {
                "maximum_relative_closed_flux_correction": 1.0e-16,
                "maximum_relative_edge_correction": 1.0e-3,
                "maximum_final_faraday_residual": 1.0e-16,
            },
        }
    report = {
        "outer_cells": [16, 32],
        "cadence_strides": [1],
        "quadrature_orders": [3],
        "cases": cases,
    }
    assessment = campaign.assess_campaign(report)
    series = assessment["spatial_convergence"]["s1_q3"]
    np.testing.assert_allclose(series["normal_flux"]["observed_orders"], (2.0,))
    np.testing.assert_allclose(series["emf"]["observed_orders"], (2.0,))
    assert assessment["smoke_gates"]["passed"]
