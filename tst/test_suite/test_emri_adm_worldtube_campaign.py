"""Tests for the ADM-driven production worldtube campaign."""

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EMRI_INPUTS = ROOT / "inputs" / "emri"
if str(EMRI_INPUTS) not in sys.path:
    sys.path.insert(0, str(EMRI_INPUTS))

import extract_global_worldtube as extract  # noqa: E402
import run_adm_worldtube_campaign as campaign  # noqa: E402


def _temporal_error(mapped: float, lower: float = 0.0, upper: float = 4.0) -> str:
    return (
        f"mapped global time {mapped:.17g} lies outside snapshot range "
        f"[{lower:.17g}, {upper:.17g}]"
    )


def test_campaign_matrix_selects_the_unique_finest_reference() -> None:
    cases = campaign.campaign_cases(
        [1.0, 0.5], [1, 2], [16, 32], [2, 3]
    )
    assert len(cases) == 16
    assert len({case.name for case in cases}) == 16
    reference = campaign.reference_case(cases)
    assert reference == campaign.CampaignCase(0.5, 1, 32, 3)
    assert reference.name == "fd0p5_s1_n32_q3"


def test_safe_time_trimming_removes_both_temporally_tilted_endpoints() -> None:
    calls = []

    def probe(times: np.ndarray) -> list[str]:
        calls.append(times.tolist())
        failures = []
        if times[0] < 1.0:
            failures.append(_temporal_error(-0.1))
        if times[-1] > 3.0:
            failures.append(_temporal_error(4.1))
        return failures

    safe, report = campaign.trim_safe_sample_times(
        np.asarray((0.0, 1.0, 2.0, 3.0, 4.0)), probe
    )
    np.testing.assert_array_equal(safe, (1.0, 2.0, 3.0))
    assert calls == [[0.0, 1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0]]
    assert report["dropped_leading_sample_times"] == [0.0]
    assert report["dropped_trailing_sample_times"] == [4.0]
    assert report["preflight_trimming_iterations"] == 1


def test_safe_time_trimming_rejects_spatial_preflight_failure() -> None:
    def probe(_: np.ndarray) -> list[str]:
        return ["worldtube stencil crosses an unavailable level-1 leaf block"]

    try:
        campaign.trim_safe_sample_times(np.asarray((0.0, 1.0, 2.0)), probe)
    except RuntimeError as error:
        assert "non-temporal reason" in str(error)
    else:
        raise AssertionError("safe-time trimming hid a spatial coverage failure")


def test_selected_indices_retain_the_final_safe_endpoint() -> None:
    assert campaign.selected_indices(6, 2) == [0, 2, 4, 5]
    assert campaign.selected_indices(5, 2) == [0, 2, 4]


def test_subsample_scan_thins_the_actual_source_snapshot_series() -> None:
    descriptors = tuple(
        extract.SnapshotDescriptor(
            time=float(index),
            cycle=index,
            lower=np.zeros(3),
            spacing=np.ones(3),
            shape_xyz=(4, 4, 4),
            state_path=Path(f"state-{index}.bin"),
            adm_path=Path(f"adm-{index}.bin"),
            source_level=0,
            source_meshblock_count=1,
            available_leaf_levels=(0,),
            source_storage="dense_uniform",
            block_shape_xyz=(4, 4, 4),
            block_logical=np.asarray(((0, 0, 0),)),
        )
        for index in range(6)
    )
    scan = extract.SnapshotManifestScan(
        path=Path("manifest.json"),
        document={},
        entries=tuple({"index": index} for index in range(6)),
        source_level=None,
        snapshot_cache_size=2,
        hash_source_files=False,
        descriptors=descriptors,
    )
    thinned, indices = campaign.subsample_scan(scan, 2)
    assert indices == [0, 2, 4, 5]
    assert [descriptor.time for descriptor in thinned.descriptors] == [
        0.0,
        2.0,
        4.0,
        5.0,
    ]
    assert [entry["index"] for entry in thinned.entries] == [0, 2, 4, 5]
