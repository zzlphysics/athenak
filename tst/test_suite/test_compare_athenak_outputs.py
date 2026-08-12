"""Tests for the closed-output restart-equivalence comparison utility."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "vis" / "python"))

import bin_convert  # noqa: E402
from compare_athenak_outputs import (  # noqa: E402
    compare_binary_files,
    compare_history_endpoints,
)


VARIABLES = ("dens", "press")
BLOCK_SHAPE = (2, 2, 2)


def _field(logical: tuple[int, int, int, int], variable: str) -> np.ndarray:
    offset = 100.0 * logical[0] + (10.0 if variable == "press" else 0.0)
    return np.arange(8, dtype=np.float64).reshape(BLOCK_SHAPE) + offset


def _write_binary(
    path: Path,
    logical_order: list[tuple[int, int, int, int]],
    *,
    perturb: tuple[tuple[int, int, int, int], str, float] | None = None,
) -> None:
    parameter_dump = (
        "<mesh>\n"
        "nx1=4\n"
        "nx2=2\n"
        "nx3=2\n"
        "nghost=0\n"
        "x1min=-1.0\n"
        "x1max=1.0\n"
        "x2min=-1.0\n"
        "x2max=1.0\n"
        "x3min=-1.0\n"
        "x3max=1.0\n"
        "<meshblock>\n"
        "nx1=2\n"
        "nx2=2\n"
        "nx3=2\n"
    ).encode("ascii")
    preheader = (
        b"Athena binary output version=1.1\n"
        b"  size of preheader=5\n"
        b"  time=1.9200000000000000e+01\n"
        b"  cycle=4\n"
        b"  size of location=8\n"
        b"  size of variable=8\n"
        b"  number of variables=2\n"
        b"  variables: dens press\n"
        + f"  header offset={len(parameter_dump)}\n".encode("ascii")
        + parameter_dump
    )
    with path.open("wb") as stream:
        stream.write(preheader)
        for logical in logical_order:
            stream.write(np.asarray([0, 1, 0, 1, 0, 1], dtype=np.int32).tobytes())
            stream.write(np.asarray(logical, dtype=np.int32).tobytes())
            geometry = np.asarray(
                [float(logical[0]), float(logical[0] + 1), -1.0, 1.0, -1.0, 1.0],
                dtype=np.float64,
            )
            stream.write(geometry.tobytes())
            for variable in VARIABLES:
                values = _field(logical, variable)
                if perturb is not None and (logical, variable) == perturb[:2]:
                    values = values.copy()
                    values.flat[0] += perturb[2]
                stream.write(values.tobytes())


def _write_history(path: Path, rows: list[tuple[float, float, float]]) -> None:
    text = "# Athena++ history data\n# [0]=time [1]=dt [2]=mass\n"
    text += "".join(
        f"{time:.17e} {dt:.17e} {mass:.17e}\n" for time, dt, mass in rows
    )
    path.write_text(text, encoding="utf-8")


def test_binary_comparison_canonicalizes_meshblock_order(tmp_path: Path) -> None:
    logical = [(0, 0, 0, 0), (1, 0, 0, 0)]
    left = tmp_path / "left.bin"
    right = tmp_path / "right.bin"
    _write_binary(left, logical)
    _write_binary(right, list(reversed(logical)))

    result = compare_binary_files(left, right)

    assert result["match"]
    assert result["array_equal"]
    assert result["topology"]["match"]
    assert set(result["fields"]) == set(VARIABLES)
    assert all(summary["l1_sum"] == 0.0 for summary in result["fields"].values())


def test_binary_comparison_reports_norms_one_variable_at_a_time(
    tmp_path: Path,
) -> None:
    logical = [(0, 0, 0, 0), (1, 0, 0, 0)]
    left = tmp_path / "left.bin"
    right = tmp_path / "right.bin"
    _write_binary(left, logical)
    _write_binary(right, logical, perturb=(logical[1], "dens", 1.0))

    result = compare_binary_files(left, right, variables=["dens"])

    assert not result["match"]
    dens = result["fields"]["dens"]
    assert not dens["array_equal"]
    np.testing.assert_allclose(dens["max_abs"], 1.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(dens["l1_sum"], 1.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        dens["relative_l1"],
        1.0
        / sum(np.sum(np.abs(_field(location, "dens"))) for location in logical),
        rtol=1.0e-15,
        atol=0.0,
    )
    assert dens["nonfinite_total"] == 0
    assert dens["max_abs_location"] == {
        "logical": [1, 0, 0, 0],
        "local_kji": [0, 0, 0],
        "left": 100.0,
        "right": 101.0,
        "signed_difference": 1.0,
        "cell_center_xyz": [1.25, -0.5, -0.5],
    }


def test_bin_convert_selective_loading_preserves_file_metadata(tmp_path: Path) -> None:
    path = tmp_path / "state.bin"
    _write_binary(path, [(0, 0, 0, 0), (1, 0, 0, 0)])

    metadata = bin_convert.read_binary(str(path), variables=())
    pressure = bin_convert.read_binary(str(path), variables=("press",))

    assert metadata["var_names"] == list(VARIABLES)
    assert metadata["loaded_var_names"] == []
    assert metadata["mb_data"] == {}
    assert pressure["var_names"] == list(VARIABLES)
    assert pressure["loaded_var_names"] == ["press"]
    assert list(pressure["mb_data"]) == ["press"]
    np.testing.assert_array_equal(
        pressure["mb_data"]["press"][1], _field((1, 0, 0, 0), "press")
    )


def test_history_comparison_uses_latest_overlapping_endpoint(tmp_path: Path) -> None:
    left = tmp_path / "continuous.hst"
    right = tmp_path / "restarted.hst"
    _write_history(
        left,
        [(0.0, 0.0, 10.0), (9.6, 4.8, 9.9), (19.2, 4.8, 9.8)],
    )
    _write_history(
        right,
        [(9.6, 4.8, 9.9), (19.2 + 5.0e-14, 4.8, 9.8), (28.8, 4.8, 9.7)],
    )

    result = compare_history_endpoints(left, right, time_atol=1.0e-12)

    assert result["match"]
    assert not result["array_equal"]
    assert result["overlap"]["left_index"] == 2
    assert result["overlap"]["right_index"] == 1
    assert result["overlap"]["left_is_terminal"]
    assert not result["overlap"]["right_is_terminal"]
    assert result["endpoint_columns"]["mass"]["array_equal"]
