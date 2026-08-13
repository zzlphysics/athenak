"""Tests for the streaming, fail-closed AthenaK restart audit."""

from __future__ import annotations

import hashlib
import importlib.util
import os
from pathlib import Path
import struct
import sys
from types import SimpleNamespace
from typing import Callable

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "audit_athenak_restart.py"
SPEC = importlib.util.spec_from_file_location("athenak_restart_auditor", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
AUDITOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = AUDITOR
SPEC.loader.exec_module(AUDITOR)


Arrays = dict[str, np.ndarray]


def _arrays() -> Arrays:
    return {
        "mhd": np.arange(6 * 4 * 4 * 4, dtype=np.float64).reshape(6, 4, 4, 4),
        "x1f": np.arange(4 * 4 * 5, dtype=np.float64).reshape(4, 4, 5) + 1000,
        "x2f": np.arange(4 * 5 * 4, dtype=np.float64).reshape(4, 5, 4) + 2000,
        "x3f": np.arange(5 * 4 * 4, dtype=np.float64).reshape(5, 4, 4) + 3000,
        "adm": (
            np.arange(17 * 4 * 4 * 4, dtype=np.float64).reshape(17, 4, 4, 4)
            + 4000
        ),
    }


def write_fixture(
    path: Path,
    *,
    mutate: Callable[[Arrays], None] | None = None,
    stored_data_size_delta: int = 0,
    trailing: bytes = b"",
) -> tuple[int, float]:
    arrays = _arrays()
    if mutate is not None:
        mutate(arrays)
    parameters = b"""<mesh>
nghost = 1
nx1 = 2
x1min = -1
x1max = 1
nx2 = 2
x2min = -1
x2max = 1
nx3 = 2
x3min = -1
x3max = 1
<meshblock>
nx1 = 2
nx2 = 2
nx3 = 2
<mesh_refinement>
refinement = none
num_levels = 1
<mhd>
eos = ideal
nscalars = 1
<adm>
dynamic = true
<time>
evolution = dynamic
<output1>
file_type = rst
single_file_per_rank = false
<par_end>
"""
    region_size = (-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    indices = (
        1, 2, 2, 2,
        1, 2, 1, 2, 1, 2,
        1, 1, 1,
        1, 1, 1, 1, 1, 1,
    )
    content = bytearray(parameters)
    content.extend(struct.pack("<ii", 1, 0))
    content.extend(struct.pack("<9d", *region_size))
    content.extend(struct.pack("<19i", *indices))
    content.extend(struct.pack("<19i", *indices))
    content.extend(struct.pack("<ddi", 9.6, 4.8, 2))
    content.extend(struct.pack("<4i", 0, 0, 0, 0))
    content.extend(struct.pack("<f", 1.0))
    field_bytes = b"".join(
        np.asarray(arrays[name], dtype="<f8").tobytes(order="C")
        for name in ("mhd", "x1f", "x2f", "x3f", "adm")
    )
    content.extend(struct.pack("<Q", len(field_bytes) + stored_data_size_delta))
    content.extend(field_bytes)
    content.extend(trailing)
    path.write_bytes(content)
    maximum = max(float(np.max(np.abs(array))) for array in arrays.values())
    return sum(array.size for array in arrays.values()), maximum


def test_valid_restart_streams_every_real_and_hashes_entire_file(
    tmp_path: Path,
) -> None:
    restart = tmp_path / "valid.rst"
    expected_count, expected_max = write_fixture(restart)

    result = AUDITOR.audit_restart(restart)

    assert result["valid"] is True
    assert result["kind"] == "athenak_restart_audit"
    assert result["sha256"] == hashlib.sha256(restart.read_bytes()).hexdigest()
    assert result["metadata"]["time"] == 9.6
    assert result["metadata"]["cycle"] == 2
    assert result["metadata"]["nmb_total"] == 1
    assert result["layout"]["expected_file_size"] == restart.stat().st_size
    assert result["stored_reals"]["count"] == expected_count
    assert result["stored_reals"]["finite_count"] == expected_count
    assert result["stored_reals"]["nonfinite_count"] == 0
    assert result["stored_reals"]["max_abs"] == expected_max


def test_nonfinite_real_fails_closed(tmp_path: Path) -> None:
    restart = tmp_path / "nonfinite.rst"

    def mutate(arrays: Arrays) -> None:
        arrays["adm"][3, 2, 1, 0] = np.nan

    write_fixture(restart, mutate=mutate)
    with pytest.raises(ValueError, match="non-finite stored Reals"):
        AUDITOR.audit_restart(restart)


@pytest.mark.parametrize(
    ("size_delta", "trailing", "pattern"),
    (
        (8, b"", "data_size mismatch"),
        (0, b"trailing", "size disagrees"),
    ),
)
def test_data_size_and_exact_eof_are_enforced(
    tmp_path: Path, size_delta: int, trailing: bytes, pattern: str
) -> None:
    restart = tmp_path / "bad-layout.rst"
    write_fixture(
        restart, stored_data_size_delta=size_delta, trailing=trailing
    )
    with pytest.raises(ValueError, match=pattern):
        AUDITOR.audit_restart(restart)


def test_truncated_meshblock_is_refused(tmp_path: Path) -> None:
    restart = tmp_path / "truncated.rst"
    write_fixture(restart)
    with restart.open("r+b") as stream:
        stream.truncate(restart.stat().st_size - 8)

    with pytest.raises(ValueError, match="size disagrees"):
        AUDITOR.audit_restart(restart)


def test_open_restart_and_symlink_are_refused(tmp_path: Path) -> None:
    restart = tmp_path / "valid.rst"
    write_fixture(restart)
    descriptor = os.open(restart, os.O_RDONLY)
    try:
        with pytest.raises(RuntimeError, match="refusing open output file"):
            AUDITOR.audit_restart(restart)
    finally:
        os.close(descriptor)

    link = tmp_path / "restart-link.rst"
    link.symlink_to(restart.name)
    with pytest.raises(ValueError, match="symbolic link"):
        AUDITOR.audit_restart(link)


def test_leaf_topology_rejects_partial_fine_octant() -> None:
    """Eight level-2 leaves cover only 1/8 of a 3-D one-block root domain."""

    metadata = SimpleNamespace(
        root_level=0,
        parameters={
            "mesh": {"nx1": "8", "nx2": "8", "nx3": "8"},
            "meshblock": {"nx1": "8", "nx2": "8", "nx3": "8"},
            "mesh_refinement": {"num_levels": "3"},
        },
        locations=tuple(
            SimpleNamespace(level=2, lx1=lx1, lx2=lx2, lx3=lx3)
            for lx3 in range(2) for lx2 in range(2) for lx1 in range(2)
        ),
    )
    with pytest.raises(ValueError, match="coverage is incomplete"):
        AUDITOR._audit_leaf_topology(metadata)
