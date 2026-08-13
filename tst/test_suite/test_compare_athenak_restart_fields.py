"""Focused tests for streaming AthenaK restart field comparison."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import struct
import sys
from typing import Callable

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "compare_athenak_restart_fields.py"
SPEC = importlib.util.spec_from_file_location("restart_field_comparator", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
FIELDS = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = FIELDS
SPEC.loader.exec_module(FIELDS)


Arrays = dict[str, np.ndarray]


def _arrays() -> Arrays:
    return {
        "mhd": np.arange(6 * 4 * 4 * 4, dtype=np.float64).reshape(6, 4, 4, 4),
        "x1f": np.arange(4 * 4 * 5, dtype=np.float64).reshape(4, 4, 5) + 1000,
        "x2f": np.arange(4 * 5 * 4, dtype=np.float64).reshape(4, 5, 4) + 2000,
        "x3f": np.arange(5 * 4 * 4, dtype=np.float64).reshape(5, 4, 4) + 3000,
        "adm": np.arange(17 * 4 * 4 * 4, dtype=np.float64).reshape(17, 4, 4, 4)
        + 4000,
    }


def write_fixture(
    path: Path,
    *,
    mutate: Callable[[Arrays], None] | None = None,
    stored_data_size_delta: int = 0,
    x1max: float = 1.0,
) -> None:
    arrays = _arrays()
    if mutate is not None:
        mutate(arrays)
    parameters = f"""<mesh>
nghost = 1
nx1 = 2
x1min = -1
x1max = {x1max}
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
""".encode("ascii")
    dx1 = (x1max + 1.0) / 2.0
    region_size = (-1.0, -1.0, -1.0, x1max, 1.0, 1.0, dx1, 1.0, 1.0)
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
    path.write_bytes(content)


def test_identical_fixture_reports_derived_layout_and_region_sizes(
    tmp_path: Path,
) -> None:
    left = tmp_path / "left.rst"
    right = tmp_path / "right.rst"
    write_fixture(left)
    write_fixture(right)

    result = FIELDS.compare_restart_fields(left, right)

    assert result["match"] is True
    assert result["authoritative_match"] is True
    assert result["all_stored_fields_match"] is True
    assert result["layout"]["nghost"] == 1
    assert result["layout"]["meshblock_shape"] == [2, 2, 2]
    assert result["layout"]["stored_shape"] == [4, 4, 4]
    assert result["layout"]["nmhd_base"] == 5
    assert result["layout"]["nscalars"] == 1
    assert result["layout"]["nadm"] == 17
    assert result["mhd_u0"]["active"]["elements"] == 6 * 2 * 2 * 2
    assert result["mhd_u0"]["ghost"]["elements"] == 6 * (4**3 - 2**3)
    assert result["face_b"]["active_faces"]["elements"] == 3 * 3 * 2 * 2
    assert result["face_b"]["ghost_faces"]["elements"] == (
        4 * 4 * 5 + 4 * 5 * 4 + 5 * 4 * 4 - 3 * 3 * 2 * 2
    )
    assert result["face_b"]["ghost_faces_compared"] is True
    assert result["adm"]["active"]["elements"] == 17 * 2 * 2 * 2
    assert result["adm"]["ghost"]["elements"] == 17 * (4**3 - 2**3)
    assert result["mhd_u0"]["active"]["first_difference"] is None


def test_differences_are_classified_by_physical_region(tmp_path: Path) -> None:
    left = tmp_path / "left.rst"
    right = tmp_path / "right.rst"
    write_fixture(left)

    def mutate(arrays: Arrays) -> None:
        arrays["mhd"][0, 1, 1, 1] += 1.5       # active cell
        arrays["mhd"][0, 1, 1, 3] += 2.5       # +x1 ghost cell
        arrays["x2f"][1, 2, 1] += 3.5           # active x2 face
        arrays["adm"][3, 0, 1, 1] += 4.5        # ADM ghost cell

    write_fixture(right, mutate=mutate)
    result = FIELDS.compare_restart_fields(left, right)

    assert result["match"] is False
    assert result["authoritative_match"] is False
    assert result["all_stored_fields_match"] is False
    assert result["mhd_u0"]["active"]["max_abs"] == pytest.approx(1.5)
    assert result["mhd_u0"]["active"]["l1_sum"] == pytest.approx(1.5)
    assert result["mhd_u0"]["active"]["first_difference"] == {
        "gid": 0,
        "logical": {"level": 0, "lx1": 0, "lx2": 0, "lx3": 0},
        "variable": 0,
        "variable_name": "IDN",
        "k": 1,
        "j": 1,
        "i": 1,
        "left": 21.0,
        "right": 22.5,
    }
    assert result["mhd_u0"]["ghost"]["max_abs"] == pytest.approx(2.5)
    assert result["mhd_u0"]["ghost"]["first_difference"]["i"] == 3
    assert result["face_b"]["active_faces"]["max_abs"] == pytest.approx(3.5)
    assert result["face_b"]["active_faces"]["first_difference"][
        "component_name"
    ] == "x2f"
    assert result["face_b"]["components"]["x1f"]["match"] is True
    assert result["face_b"]["components"]["x2f"]["match"] is False
    assert result["adm"]["active"]["match"] is True
    assert result["adm"]["ghost"]["max_abs"] == pytest.approx(4.5)


def test_nonfinite_is_reported_and_fails_even_without_a_finite_delta(
    tmp_path: Path,
) -> None:
    left = tmp_path / "left.rst"
    right = tmp_path / "right.rst"
    write_fixture(left)

    def mutate(arrays: Arrays) -> None:
        arrays["adm"][0, 1, 1, 1] = np.nan

    write_fixture(right, mutate=mutate)
    result = FIELDS.compare_restart_fields(left, right)

    active = result["adm"]["active"]
    assert result["match"] is False
    assert active["right_nonfinite"]["nan"] == 1
    assert active["left_nonfinite"]["nan"] == 0
    assert active["nonfinite_class_mismatches"] == 1
    assert active["within_tolerance"] is False
    assert active["first_difference"]["right"] == "NaN"


def test_derived_ghost_difference_is_distinct_from_authoritative_state(
    tmp_path: Path,
) -> None:
    left = tmp_path / "left.rst"
    right = tmp_path / "right.rst"
    write_fixture(left)

    def mutate(arrays: Arrays) -> None:
        arrays["mhd"][0, 0, 1, 1] += 0.5  # ghost cell only

    write_fixture(right, mutate=mutate)
    result = FIELDS.compare_restart_fields(left, right)

    assert result["match"] is False
    assert result["all_stored_fields_match"] is False
    assert result["authoritative_match"] is True
    assert result["mhd_u0"]["active"]["match"] is True
    assert result["mhd_u0"]["ghost"]["match"] is False


def test_face_ghost_difference_fails_all_stored_fields_only(tmp_path: Path) -> None:
    left = tmp_path / "left.rst"
    right = tmp_path / "right.rst"
    write_fixture(left)

    def mutate(arrays: Arrays) -> None:
        arrays["x1f"][0, 0, 0] += 0.75  # face ghost only

    write_fixture(right, mutate=mutate)
    result = FIELDS.compare_restart_fields(left, right)

    assert result["match"] is False
    assert result["all_stored_fields_match"] is False
    assert result["authoritative_match"] is True
    assert result["face_b"]["active_faces"]["match"] is True
    assert result["face_b"]["ghost_faces"]["match"] is False
    assert result["face_b"]["ghost_components"]["x1f"]["max_abs"] == pytest.approx(0.75)


def test_stored_data_size_and_topology_are_fail_closed(tmp_path: Path) -> None:
    left = tmp_path / "left.rst"
    bad_size = tmp_path / "bad-size.rst"
    different_topology = tmp_path / "different-topology.rst"
    write_fixture(left)
    write_fixture(bad_size, stored_data_size_delta=8)
    write_fixture(different_topology, x1max=2.0)

    with pytest.raises(ValueError, match="data_size mismatch"):
        FIELDS.compare_restart_fields(left, bad_size)
    with pytest.raises(ValueError, match="topology or field ABI differs"):
        FIELDS.compare_restart_fields(left, different_topology)


def test_cli_returns_one_for_a_physical_mismatch(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    left = tmp_path / "left.rst"
    right = tmp_path / "right.rst"
    write_fixture(left)

    def mutate(arrays: Arrays) -> None:
        arrays["mhd"][1, 1, 1, 1] += 0.25

    write_fixture(right, mutate=mutate)

    assert FIELDS.main([str(left), str(right), "--compact"]) == 1
    output = capsys.readouterr()
    assert '"kind": "athenak_restart_field_comparison"' in output.out
    assert '"match": false' in output.out
