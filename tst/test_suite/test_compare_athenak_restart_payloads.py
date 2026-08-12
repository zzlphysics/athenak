"""Tests for bounded, fail-closed AthenaK restart payload comparison."""

from __future__ import annotations

from dataclasses import replace
import importlib.util
from pathlib import Path
import struct
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
COMPARATOR_PATH = ROOT / "scripts" / "compare_athenak_restart_payloads.py"
SPEC = importlib.util.spec_from_file_location(
    "restart_payload_comparator", COMPARATOR_PATH
)
assert SPEC is not None and SPEC.loader is not None
COMPARATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = COMPARATOR
SPEC.loader.exec_module(COMPARATOR)


def write_fixture(
    path: Path,
    *,
    label: str = "alpha",
    payload: bytes = b"OBJECT_STATE\x00SIZE\x00" + struct.pack("<4d", 1, 2, 3, 4),
    first_cost: float = 2.0,
    first_event: int = 11,
) -> int:
    """Write valid refined metadata followed by an intentionally opaque payload."""

    if len(label) != 5:
        raise ValueError("test labels must have equal length")
    parameters = f"""<mesh>
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
refinement = adaptive
num_levels = 2
max_nmb_per_rank = 8
<problem>
label = {label}
<output1>
file_type = rst
single_file_per_rank = false
<par_end>
""".encode("ascii")
    region_size = (-1.0, -1.0, -1.0, 1.0, 1.0, 1.0,
                   0.125, 0.125, 0.125)
    indices = (4, 16, 16, 16, 4, 19, 4, 19, 4, 19,
               8, 8, 8, 4, 11, 4, 11, 4, 11)
    locations = [
        (lx1, lx2, lx3, 1)
        for lx3 in range(2)
        for lx2 in range(2)
        for lx1 in range(2)
    ]
    costs = [first_cost] + [2.0] * 7
    counters = tuple(range(8))

    content = bytearray(parameters)
    content.extend(struct.pack("<ii", 8, 0))
    content.extend(struct.pack("<9d", *region_size))
    content.extend(struct.pack("<19i", *indices))
    content.extend(struct.pack("<19i", *indices))
    content.extend(struct.pack("<ddi", 9.6, 4.8, 2))
    for location in locations:
        content.extend(struct.pack("<4i", *location))
    content.extend(struct.pack("<8f", *costs))
    content.extend(struct.pack(
        "<Qii",
        COMPARATOR.AMR_COUNTER_MAGIC,
        COMPARATOR.AMR_COUNTER_VERSION,
        8,
    ))
    content.extend(struct.pack("<8i", *counters))
    event_sums = [first_event] + list(
        range(12, 11 + COMPARATOR.EVENT_SUM_COUNTER_COUNT)
    )
    content.extend(struct.pack(
        "<Qii",
        COMPARATOR.EVENT_COUNTER_MAGIC,
        COMPARATOR.EVENT_COUNTER_VERSION,
        COMPARATOR.EVENT_SUM_COUNTER_COUNT,
    ))
    content.extend(struct.pack(
        f"<{COMPARATOR.EVENT_SUM_COUNTER_COUNT}Q", *event_sums
    ))
    content.extend(struct.pack("<i", 7))
    metadata_end = len(content)
    content.extend(payload)
    path.write_bytes(content)
    return metadata_end


def test_parameter_text_may_differ_when_payload_and_topology_match(
    tmp_path: Path,
) -> None:
    left = tmp_path / "left.rst"
    right = tmp_path / "right.rst"
    payload = b"opaque restart payload" * 13
    left_end = write_fixture(left, label="alpha", payload=payload)
    right_end = write_fixture(right, label="bravo", payload=payload)

    result = COMPARATOR.compare_restart_payloads(
        left, right, chunk_bytes=17
    )

    assert left_end == right_end
    assert result["match"] is True
    assert result["payload_equal"] is True
    assert result["payload_bytes"] == {
        "left": len(payload),
        "right": len(payload),
    }
    assert result["sha256"]["left"] == result["sha256"]["right"]
    assert result["first_difference_offset"] is None
    assert result["metadata"]["match"] is True
    assert all(result["metadata"]["fields"].values())
    assert result["field_data_read"] is True
    assert result["field_data_bytes_read"] == {
        "left": len(payload),
        "right": len(payload),
    }
    assert "parameters" not in result


def test_payload_difference_reports_relative_offset_and_exit_one(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    left = tmp_path / "left.rst"
    right = tmp_path / "right.rst"
    left_payload = bytearray(b"same-prefix-and-tail")
    right_payload = bytearray(left_payload)
    right_payload[7] ^= 0x01
    write_fixture(left, payload=bytes(left_payload))
    write_fixture(right, payload=bytes(right_payload))

    exit_code = COMPARATOR.main([
        str(left), str(right), "--chunk-bytes", "5", "--compact"
    ])
    output = capsys.readouterr()

    assert exit_code == 1
    assert '"payload_equal":false' in output.out
    result = COMPARATOR.compare_restart_payloads(left, right, chunk_bytes=5)
    assert result["match"] is False
    assert result["payload_equal"] is False
    assert result["first_difference_offset"] == 7
    assert result["sha256"]["left"] != result["sha256"]["right"]


def test_changed_file_is_rejected_fail_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    left = tmp_path / "left.rst"
    right = tmp_path / "right.rst"
    write_fixture(left)
    write_fixture(right)
    original_signature = COMPARATOR._file_signature
    calls = 0

    def changing_signature(path: Path):
        nonlocal calls
        calls += 1
        signature = original_signature(path)
        if calls == 5:
            return replace(signature, mtime_ns=signature.mtime_ns + 1)
        return signature

    monkeypatch.setattr(COMPARATOR, "_file_signature", changing_signature)

    with pytest.raises(RuntimeError, match="restart changed during comparison"):
        COMPARATOR.compare_restart_payloads(left, right, chunk_bytes=7)


def test_topology_difference_is_a_valid_mismatch(tmp_path: Path) -> None:
    left = tmp_path / "left.rst"
    right = tmp_path / "right.rst"
    payload = b"identical-field-state"
    write_fixture(left, payload=payload, first_cost=2.0)
    write_fixture(right, payload=payload, first_cost=4.0)

    result = COMPARATOR.compare_restart_payloads(left, right, chunk_bytes=8)

    assert result["match"] is False
    assert result["metadata"]["match"] is False
    assert result["metadata"]["fields"]["costs"] is False


def test_pending_event_difference_is_a_metadata_mismatch(tmp_path: Path) -> None:
    left = tmp_path / "left.rst"
    right = tmp_path / "right.rst"
    write_fixture(left, first_event=11)
    write_fixture(right, first_event=12)

    result = COMPARATOR.compare_restart_payloads(left, right, chunk_bytes=11)

    assert result["match"] is False
    assert result["payload_equal"] is True
    assert result["metadata"]["fields"]["event_counters"] is False
    assert result["payload_equal"] is True
    assert result["first_difference_offset"] is None
