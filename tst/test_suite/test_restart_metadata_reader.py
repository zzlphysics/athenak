"""Tests for the field-free AthenaK restart metadata reader."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import struct
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
READER_PATH = ROOT / "scripts" / "read_athenak_restart_metadata.py"
SPEC = importlib.util.spec_from_file_location("restart_metadata_reader", READER_PATH)
assert SPEC is not None and SPEC.loader is not None
READER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = READER
SPEC.loader.exec_module(READER)


def write_fixture(
        path: Path, *, real_bytes: int = 8, refinement: str = "adaptive",
        nmb_total: int = 8, counter_marker: bool = True,
        counter_key: bool = True,
        event_marker: bool = True,
        event_version: int | None = None,
        event_count: int | None = None,
        event_maxit: int = 7,
        cache_contract_marker: bool = True,
        cache_contract_version: int | None = None,
        single_file_per_rank: str | None = None,
        sentinel: bytes = b"FIELD_DATA_MUST_NOT_BE_READ" * 128) -> int:
    """Write metadata for one refined root block followed by an opaque field payload."""

    if real_bytes not in (4, 8):
        raise ValueError("fixture Real must contain four or eight bytes")
    if refinement not in ("adaptive", "static"):
        raise ValueError("fixture refinement must be adaptive or static")
    if nmb_total not in (1, 8):
        raise ValueError("fixture supports one root or one refined octet")

    refinement_lines = [
        f"refinement = {refinement}",
        "max_nmb_per_rank = 4",
    ]
    if refinement == "adaptive":
        refinement_lines.append("num_levels = 2")
    if counter_key:
        refinement_lines.append("restart_cycle_counters_version = 1")
    refinement_parameters = ("\n".join(refinement_lines) + "\n").encode("ascii")
    output_parameters = b"file_type = rst\n"
    if single_file_per_rank is not None:
        output_parameters += (
            f"single_file_per_rank = {single_file_per_rank}\n".encode("ascii")
        )
    parameters = b"""<mesh>
nghost = 4
nx1 = 16
x1min = -1
x1max = 1
nx2 = 16
x2min = -1
x2max = 1
nx3 = 16
x3min = -1
x3max = 1
<meshblock>
nx1 = 16
nx2 = 16
nx3 = 16
<mesh_refinement>
""" + refinement_parameters + b"""<time>
subcycling = level
<output1>
""" + output_parameters + b"""<par_end>
"""
    region_size = (-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 0.125, 0.125, 0.125)
    mesh_indices = (4, 16, 16, 16, 4, 19, 4, 19, 4, 19,
                    8, 8, 8, 4, 11, 4, 11, 4, 11)
    block_indices = mesh_indices
    if nmb_total == 1:
        locations = [(0, 0, 0, 0)]
        costs = [1.0]
        counters = [7]
    else:
        locations = [
            (lx1, lx2, lx3, 1)
            for lx3 in range(2)
            for lx2 in range(2)
            for lx1 in range(2)
        ]
        costs = [2.0] * len(locations)
        counters = [0, 0, 1, 1, 2, 2, 3, 3]
    real_code = "f" if real_bytes == 4 else "d"
    content = bytearray(parameters)
    content.extend(struct.pack("<ii", len(locations), 0))
    content.extend(struct.pack(f"<9{real_code}", *region_size))
    content.extend(struct.pack("<19i", *mesh_indices))
    content.extend(struct.pack("<19i", *block_indices))
    content.extend(struct.pack(f"<{real_code}{real_code}i", 9.6, 4.8, 2))
    for location in locations:
        content.extend(struct.pack("<4i", *location))
    for cost in costs:
        content.extend(struct.pack("<f", cost))
    if counter_marker:
        content.extend(struct.pack(
            "<Qii", READER.AMR_COUNTER_MAGIC, READER.AMR_COUNTER_VERSION,
            len(locations)
        ))
        content.extend(struct.pack(f"<{len(counters)}i", *counters))
    if event_marker:
        event_sums = tuple(range(11, 11 + READER.EVENT_SUM_COUNTER_COUNT))
        content.extend(struct.pack(
            "<Qii",
            READER.EVENT_COUNTER_MAGIC,
            (READER.EVENT_COUNTER_VERSION
             if event_version is None else event_version),
            (READER.EVENT_SUM_COUNTER_COUNT
             if event_count is None else event_count),
        ))
        content.extend(struct.pack(
            f"<{READER.EVENT_SUM_COUNTER_COUNT}Q", *event_sums
        ))
        content.extend(struct.pack("<i", event_maxit))
    if cache_contract_marker:
        content.extend(struct.pack(
            "<Qi",
            READER.RESTART_CACHE_CONTRACT_MAGIC,
            (READER.RESTART_CACHE_CONTRACT_VERSION
             if cache_contract_version is None else cache_contract_version),
        ))
    metadata_end = len(content)
    content.extend(sentinel)
    path.write_bytes(content)
    return metadata_end


def test_restart_metadata_stops_before_field_data(tmp_path: Path) -> None:
    restart = tmp_path / "fixture.rst"
    expected_end = write_fixture(restart)

    metadata = READER.read_restart_metadata(restart)
    payload = READER.metadata_dict(metadata, nranks=2, max_blocks_per_rank=4)

    assert metadata.metadata_end == expected_end
    assert metadata.max_read_offset == expected_end
    assert metadata.max_read_offset < restart.stat().st_size
    assert metadata.nmb_total == 8
    assert metadata.root_level == 0
    assert metadata.time == 9.6
    assert metadata.cycle == 2
    assert metadata.physical_level_rows() == [{
        "physical_level": 1,
        "logical_level": 1,
        "blocks": 8,
        "cost": 16.0,
    }]
    assert payload["partition"]["ranks"] == [
        {"rank": 0, "gid_start": 0, "blocks": 4, "cost": 8.0},
        {"rank": 1, "gid_start": 4, "blocks": 4, "cost": 8.0},
    ]
    assert payload["event_counters"] == {
        "neos_dfloor": 11,
        "neos_efloor": 12,
        "neos_tfloor": 13,
        "neos_vceil": 14,
        "neos_fail": 15,
        "nfofc": 16,
        "ncons_adjust": 17,
        "nmag_adjust": 18,
        "nc2p_calls": 19,
        "nfofc_tests": 20,
        "maxit_c2p": 7,
    }
    assert metadata.event_counter_version == READER.EVENT_COUNTER_VERSION
    assert payload["event_counter_version"] == READER.EVENT_COUNTER_VERSION
    assert metadata.restart_cache_contract_version == 1
    assert payload["restart_cache_contract_version"] == 1


def test_capacity_partition_matches_cpp_regression() -> None:
    costs = [1.0] * 6 + [4.0, 4.0, 5.0, 5.0, 5.0]

    partitions = READER.load_balance(costs, nranks=2, max_blocks_per_rank=6)

    assert [(row.gid_start, row.blocks, row.cost) for row in partitions] == [
        (0, 6, 6.0),
        (6, 5, 23.0),
    ]


def test_static_refinement_does_not_require_num_levels(tmp_path: Path) -> None:
    restart = tmp_path / "static_fixture.rst"
    expected_end = write_fixture(
        restart,
        refinement="static",
        counter_marker=False,
        counter_key=False,
        event_marker=False,
        cache_contract_marker=False,
    )

    metadata = READER.read_restart_metadata(restart)

    assert metadata.metadata_end == expected_end
    assert metadata.max_read_offset == expected_end + 8
    assert metadata.nmb_total == 8
    assert metadata.amr_cycle_counters is None
    assert metadata.physical_level_rows()[0]["physical_level"] == 1


@pytest.mark.parametrize(
    ("counter_marker", "counter_key"),
    ((False, False), (False, True), (True, False), (True, True)),
)
def test_counter_extension_is_selected_only_by_serialized_magic(
        tmp_path: Path, counter_marker: bool, counter_key: bool) -> None:
    restart = tmp_path / f"counter_{int(counter_marker)}_{int(counter_key)}.rst"
    expected_end = write_fixture(
        restart, counter_marker=counter_marker, counter_key=counter_key
    )

    metadata = READER.read_restart_metadata(restart)

    assert metadata.metadata_end == expected_end
    assert (metadata.amr_cycle_counters is not None) == counter_marker
    if counter_marker:
        assert metadata.amr_cycle_counters == (0, 0, 1, 1, 2, 2, 3, 3)
        assert metadata.max_read_offset == expected_end
    else:
        assert metadata.max_read_offset == expected_end


def test_real4_single_block_does_not_probe_real8_or_sentinel(tmp_path: Path) -> None:
    restart = tmp_path / "real4_single_block.rst"
    expected_end = write_fixture(
        restart,
        real_bytes=4,
        nmb_total=1,
        counter_marker=True,
        counter_key=False,
        sentinel=b"REAL8_PROBE_MUST_NOT_TOUCH_THIS_SENTINEL",
    )

    metadata = READER.read_restart_metadata(restart)

    assert metadata.real_bytes == 4
    assert metadata.nmb_total == 1
    assert math.isclose(metadata.time, 9.6, rel_tol=2.0e-6, abs_tol=0.0)
    assert metadata.metadata_end == expected_end
    assert metadata.max_read_offset == expected_end
    assert metadata.amr_cycle_counters == (7,)


def test_legacy_restart_without_event_extension_restores_empty_state(
        tmp_path: Path) -> None:
    restart = tmp_path / "legacy_no_events.rst"
    expected_end = write_fixture(
        restart, event_marker=False, cache_contract_marker=False
    )

    metadata = READER.read_restart_metadata(restart)

    assert metadata.event_counters is None
    assert metadata.event_counter_version is None
    assert metadata.restart_cache_contract_version is None
    assert metadata.metadata_end == expected_end
    assert metadata.max_read_offset == expected_end + 8


def test_restart_cache_contract_can_follow_legacy_missing_event_extension(
        tmp_path: Path) -> None:
    restart = tmp_path / "cache_contract_without_events.rst"
    expected_end = write_fixture(restart, event_marker=False)

    metadata = READER.read_restart_metadata(restart)

    assert metadata.event_counters is None
    assert metadata.restart_cache_contract_version == 1
    assert metadata.metadata_end == expected_end
    assert metadata.max_read_offset == expected_end


def test_legacy_restart_without_cache_contract_rolls_back_before_fields(
        tmp_path: Path) -> None:
    restart = tmp_path / "legacy_no_cache_contract.rst"
    expected_end = write_fixture(restart, cache_contract_marker=False)

    metadata = READER.read_restart_metadata(restart)

    assert metadata.event_counters is not None
    assert metadata.restart_cache_contract_version is None
    assert metadata.metadata_end == expected_end
    assert metadata.max_read_offset == expected_end + 8


def test_invalid_restart_cache_contract_version_fails_closed(
        tmp_path: Path) -> None:
    restart = tmp_path / "invalid_cache_contract.rst"
    write_fixture(restart, cache_contract_version=2)

    with pytest.raises(ValueError, match="invalid restart cache-contract extension"):
        READER.read_restart_metadata(restart)


def test_legacy_ghost_event_counter_version_is_identified(tmp_path: Path) -> None:
    restart = tmp_path / "legacy_v1_events.rst"
    write_fixture(restart, event_version=READER.LEGACY_EVENT_COUNTER_VERSION)

    metadata = READER.read_restart_metadata(restart)

    assert metadata.event_counters is not None
    assert metadata.event_counter_version == READER.LEGACY_EVENT_COUNTER_VERSION


@pytest.mark.parametrize(
    ("event_version", "event_count", "pattern"),
    (
        (3, None, "invalid event-counter extension"),
        (None, 9, "invalid event-counter extension"),
    ),
)
def test_invalid_event_counter_header_fails_closed(
        tmp_path: Path, event_version: int | None, event_count: int | None,
        pattern: str) -> None:
    restart = tmp_path / "invalid_event_header.rst"
    write_fixture(
        restart, event_version=event_version, event_count=event_count
    )

    with pytest.raises(ValueError, match=pattern):
        READER.read_restart_metadata(restart)


def test_negative_event_iteration_maximum_fails_closed(tmp_path: Path) -> None:
    restart = tmp_path / "negative_event_maxit.rst"
    write_fixture(restart, event_maxit=-1)

    with pytest.raises(ValueError, match="negative maximum C2P iteration"):
        READER.read_restart_metadata(restart)


@pytest.mark.parametrize("value", ("true", "1"))
def test_per_rank_restart_boolean_true_is_rejected(
        tmp_path: Path, value: str) -> None:
    restart = tmp_path / f"per_rank_{value}.rst"
    write_fixture(restart, single_file_per_rank=value)

    with pytest.raises(ValueError, match="per-rank restart files are not supported"):
        READER.read_restart_metadata(restart)


@pytest.mark.parametrize("value", ("false", "0"))
def test_per_rank_restart_boolean_false_is_accepted(
        tmp_path: Path, value: str) -> None:
    restart = tmp_path / f"single_file_{value}.rst"
    write_fixture(restart, single_file_per_rank=value)

    assert READER.read_restart_metadata(restart).nmb_total == 8


def test_invalid_per_rank_restart_boolean_fails_closed(tmp_path: Path) -> None:
    restart = tmp_path / "invalid_boolean.rst"
    write_fixture(restart, single_file_per_rank="yes")

    with pytest.raises(ValueError, match="invalid boolean value"):
        READER.read_restart_metadata(restart)


def test_counter_array_is_included_in_metadata_capacity_limit(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    restart = tmp_path / "counter_capacity.rst"
    write_fixture(restart, counter_marker=True)
    assert READER._mesh_metadata_size(8, include_counters=False) == 160
    assert READER._mesh_metadata_size(8, include_counters=True) == 208
    assert READER._mesh_metadata_size(
        8, include_counters=True, include_event_counters=True
    ) == 308
    assert READER._mesh_metadata_size(
        8,
        include_counters=True,
        include_event_counters=True,
        include_restart_cache_contract=True,
    ) == 320
    monkeypatch.setattr(READER, "MAX_MESH_METADATA_BYTES", 310)

    with pytest.raises(
        ValueError,
        match=(
            "including AMR counters and event counters and restart cache contract"
        ),
    ):
        READER.read_restart_metadata(restart)


def test_partition_uses_cpp_ordered_binary64_sum() -> None:
    costs = [1.0, float(2**53), 1.0]

    partitions = READER.load_balance(
        costs, nranks=2, max_blocks_per_rank=3
    )

    assert [(row.gid_start, row.blocks) for row in partitions] == [
        (0, 2),
        (2, 1),
    ]
