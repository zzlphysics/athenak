"""Tests for the numerical-ADM force convergence driver."""

from argparse import Namespace
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import run_adm_inner_replay_convergence as convergence  # noqa: E402


def test_case_matrix_has_unique_finest_cadence_reference() -> None:
    cases = convergence.convergence_cases([1, 2, 4], [1, 2])
    assert len(cases) == 6
    assert convergence.reference_case(cases).name == "mr4_s1"
    comparisons = {}
    for factor, error in ((1, 0.04), (2, 0.01), (4, 0.0)):
        comparisons[f"mr{factor}_s1"] = {
            "comparison": {
                "history": {
                    name: {"relative_l2": error}
                    for name in convergence.HISTORY_GROUPS
                },
                "adm": {
                    name: {"relative_l2": error}
                    for name in convergence.ADM_GROUPS
                },
            }
        }
    orders = convergence.observed_spatial_orders(comparisons, [1, 2, 4], 1)
    assert np.isclose(orders["adm"]["metric"][0], 2.0)
    assert np.isclose(orders["history"]["mdot"][0], 2.0)


def test_series_comparison_interpolates_and_normalizes_vector_error() -> None:
    reference_times = np.asarray((0.0, 1.0))
    reference = np.asarray(((1.0, 0.0), (1.0, 0.0)))
    candidate_times = np.asarray((0.0, 0.5, 1.0))
    candidate = np.asarray(((1.1, 0.0), (1.1, 0.0), (1.1, 0.0)))
    comparison = convergence.compare_series(
        reference_times, reference, candidate_times, candidate
    )
    assert np.isclose(comparison["absolute_rms"], 0.1)
    assert np.isclose(comparison["relative_l2"], 0.1)
    assert comparison["comparison_knots"] == 3

    zero_reference = np.zeros_like(reference)
    undefined = convergence.compare_series(
        reference_times, zero_reference, candidate_times, candidate
    )
    assert undefined["relative_l2"] is None


def test_assessment_separates_smoke_pass_from_science_readiness() -> None:
    history = {
        name: {"relative_l2": 0.01}
        for name in convergence.HISTORY_GROUPS
    }
    comparison = {
        "history": history,
        "adm": {
            "metric": {"relative_l2": 0.01},
            "extrinsic_curvature": {"relative_l2": 0.02},
        },
    }
    comparisons = {
        "reference": {"reference": True, "comparison": comparison},
        "candidate": {"reference": False, "comparison": comparison},
    }
    documents = {
        name: {
            "run_status": "completed",
            "structural_assessment": {"passed": True},
            "runtime_assessment": {"passed": True},
        }
        for name in comparisons
    }
    arguments = Namespace(
        maximum_metric_relative_l2=0.05,
        maximum_k_relative_l2=0.05,
        maximum_mdot_relative_l2=0.05,
        maximum_force_relative_l2=0.05,
    )
    assessment = convergence.assess_convergence(
        documents, comparisons, [1, 2], [1], arguments
    )
    assert assessment["passed"]
    assert not assessment["science_ready"]
    assert not assessment["production_matrix"]
