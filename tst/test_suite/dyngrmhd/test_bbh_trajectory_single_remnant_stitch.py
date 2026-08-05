"""Pure-Python tests for the canonical single-remnant trajectory stitcher."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
SCRIPT = REPOSITORY / "scripts" / "bbh_trajectory_stitch_single_remnant.py"
MODULE_SPEC = importlib.util.spec_from_file_location("single_remnant_stitch", SCRIPT)
assert MODULE_SPEC is not None and MODULE_SPEC.loader is not None
STITCHER = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(STITCHER)
FINAL_MASS = 0.95
FINAL_A = np.asarray([0.0, 0.0, 0.60])
KICK = np.asarray([0.002, -0.001, 0.0005])


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def _inspiral_table() -> np.ndarray:
    """Return a simple differentiable, normalized q=1 inspiral table."""
    time = np.linspace(0.0, 10.0, 1001)
    position1 = np.column_stack(
        (1.0 - 0.01 * time, 0.03 * time, np.zeros_like(time))
    )
    position2 = -position1
    velocity1 = np.tile(np.asarray([-0.01, 0.03, 0.0]), (time.size, 1))
    velocity2 = -velocity1
    zeros = np.zeros((time.size, 3))
    masses = np.full((time.size, 2), 0.5)
    return np.column_stack(
        (
            time,
            position1,
            position2,
            velocity1,
            velocity2,
            zeros,
            zeros,
            masses,
        )
    )


def _write_input(path: Path, table: np.ndarray | None = None) -> np.ndarray:
    if table is None:
        table = _inspiral_table()
    np.savetxt(path, table, fmt="%.17e")
    return table


def _command(
    input_path: Path,
    output_path: Path,
    *,
    use_separation: bool = False,
    source_labels: tuple[str, ...] = ("synthetic-linear-inspiral",),
    source_revision: str = "test-fixture-v1",
    source_provenance: Path | None = None,
) -> list[str]:
    command = [sys.executable, str(SCRIPT), str(input_path), str(output_path)]
    for label in source_labels:
        command.extend(["--source-model-label", label])
    command.extend(["--source-revision", source_revision])
    if source_provenance is not None:
        command.extend(["--source-provenance", str(source_provenance)])
    if use_separation:
        table = np.loadtxt(input_path)
        row = 200
        separation = np.linalg.norm(table[row, 1:4] - table[row, 4:7])
        command.extend(["--transition-start-separation", repr(float(separation))])
    else:
        command.extend(["--transition-start-time", "2.0"])
    command.extend(
        [
            "--transition-end-time",
            "8.0",
            "--final-mass",
            str(FINAL_MASS),
            "--final-a",
            *(str(value) for value in FINAL_A),
            "--kick",
            *(str(value) for value in KICK),
            "--postmerger-duration",
            "2.0",
            "--postmerger-dt",
            "0.01",
            "--max-relative-dxdt-v-error",
            "0.01",
        ]
    )
    return command


def _run(command: list[str], *, success: bool = True) -> subprocess.CompletedProcess:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if success and result.returncode != 0:
        pytest.fail(
            f"stitcher failed:\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    if not success and result.returncode == 0:
        pytest.fail("stitcher unexpectedly accepted invalid input")
    return result


def test_stitch_endpoints_derivative_hash_and_separation_selection(tmp_path):
    """The stitched worldlines have exact endpoints and analytic velocities."""
    input_path = tmp_path / "inspiral.dat"
    input_table = _write_input(input_path)
    source_provenance = tmp_path / "source-manifest.json"
    source_provenance.write_text('{"fixture": 1}\n', encoding="utf-8")
    output1 = tmp_path / "stitched-a.dat"
    output2 = tmp_path / "stitched-b.dat"
    _run(
        _command(
            input_path,
            output1,
            use_separation=True,
            source_provenance=source_provenance,
        )
    )
    _run(
        _command(
            input_path,
            output2,
            use_separation=True,
            source_provenance=source_provenance,
        )
    )

    stitched = np.loadtxt(output1)
    start_row = int(np.flatnonzero(stitched[:, 0] == 2.0)[0])
    end_row = int(np.flatnonzero(stitched[:, 0] == 8.0)[0])
    source_start = input_table[int(np.flatnonzero(input_table[:, 0] == 2.0)[0])]
    np.testing.assert_allclose(stitched[start_row], source_start, rtol=0.0, atol=2e-16)

    com_start = 0.5 * (source_start[1:4] + source_start[4:7])
    remnant_at_end = com_start + KICK * (8.0 - 2.0)
    np.testing.assert_allclose(stitched[end_row, 1:4], remnant_at_end, atol=2e-16)
    np.testing.assert_allclose(stitched[end_row, 4:7], remnant_at_end, atol=2e-16)
    np.testing.assert_allclose(stitched[end_row, 7:10], KICK, atol=2e-16)
    np.testing.assert_allclose(stitched[end_row, 10:13], KICK, atol=2e-16)
    np.testing.assert_allclose(stitched[end_row, 13:16], FINAL_A, atol=2e-16)
    np.testing.assert_allclose(stitched[end_row, 16:19], 0.0, atol=0.0)
    np.testing.assert_allclose(
        stitched[end_row, 19:21], [FINAL_MASS, 0.0], atol=0.0
    )
    np.testing.assert_allclose(
        stitched[-1, 1:4], com_start + KICK * 8.0, atol=2e-16
    )
    np.testing.assert_allclose(stitched[-1, 1:4], stitched[-1, 4:7], atol=0.0)

    for position_slice, velocity_slice in ((slice(1, 4), slice(7, 10)),
                                            (slice(4, 7), slice(10, 13))):
        numerical_velocity = np.column_stack(
            [
                np.gradient(stitched[:, component], stitched[:, 0], edge_order=2)
                for component in range(position_slice.start, position_slice.stop)
            ]
        )
        error = np.linalg.norm(
            numerical_velocity[1:-1] - stitched[1:-1, velocity_slice], axis=1
        )
        assert np.max(error) < 2.0e-4

    assert _sha256(output1) == _sha256(output2)
    provenance1 = json.loads(Path(str(output1) + ".provenance.json").read_text())
    provenance2 = json.loads(Path(str(output2) + ".provenance.json").read_text())
    assert provenance1["output"]["sha256"] == _sha256(output1)
    assert provenance2["output"]["sha256"] == _sha256(output2)
    assert provenance1["output"]["sha256"] == provenance2["output"]["sha256"]
    selection = provenance1["transition"]["start_selection"]
    assert selection["mode"] == "first-row-at-or-below-separation"
    assert selection["selected_row"] == 200
    assert provenance1["remnant"]["representation"] == "canonical-single-term-1"
    np.testing.assert_allclose(
        provenance1["remnant"]["chi"], FINAL_A / FINAL_MASS, atol=0.0
    )
    assert provenance1["remnant"]["a_to_chi_relation"] == (
        "a_length = Mf * chi; chi = a_length / Mf"
    )
    source_model = provenance1["trajectory_model"]
    assert source_model["declaration_source"] == "CLI"
    assert source_model["labels"] == ["synthetic-linear-inspiral"]
    assert "frozen-CBwaves" not in source_model["labels"]
    assert "local-instantaneous-4PN" not in source_model["labels"]
    assert source_model["source_revision"] == "test-fixture-v1"
    assert source_model["source_provenance"] == {
        "path": str(source_provenance.resolve()),
        "sha256": _sha256(source_provenance),
    }
    transition = provenance1["transition"]
    assert transition["sampling"]["row_count"] == 601
    assert transition["sampling"]["interval_count"] == 600
    assert transition["sampling"]["uniform_within_1e-12"]
    assert transition["sampling"]["median_dt"] == pytest.approx(0.01)
    assert transition["weight"]["formula_version"] == "cinf-bump-ratio-v1"
    diagnostic = provenance1["velocity_consistency"]
    assert diagnostic["row_count"] == stitched.shape[0]
    assert diagnostic["relative_l2"]["max"] < 0.01
    assert diagnostic["threshold"] == {
        "enforced": True,
        "limit": 0.01,
        "metric": "relative_l2.max",
        "passed": True,
    }

    original_hash = _sha256(output1)
    provenance_path = Path(str(output1) + ".provenance.json")
    original_provenance = provenance_path.read_bytes()
    refused = _run(
        _command(
            input_path,
            output1,
            use_separation=True,
            source_provenance=source_provenance,
        ),
        success=False,
    )
    assert "refusing to replace" in refused.stderr
    assert _sha256(output1) == original_hash

    forced = _command(
        input_path,
        output1,
        use_separation=True,
        source_provenance=source_provenance,
    )
    forced.append("--force")
    _run(forced)
    assert _sha256(output1) == original_hash
    assert provenance_path.read_bytes() == original_provenance


@pytest.mark.parametrize(
    ("old", "new", "message"),
    [
        (("--final-mass", "0.95"), ("--final-mass", "0.0"), "final-mass"),
        (
            ("--final-a", "0.0", "0.0", "0.6"),
            ("--final-a", "0.0", "0.0", "1.0"),
            "|a_f|",
        ),
        (
            ("--kick", "0.002", "-0.001", "0.0005"),
            ("--kick", "1.0", "0.0", "0.0"),
            "speed of light",
        ),
        (
            ("--transition-end-time", "8.0"),
            ("--transition-end-time", "1.0"),
            "transition end",
        ),
        (
            ("--postmerger-duration", "2.0"),
            ("--postmerger-duration", "-1.0"),
            "postmerger-duration",
        ),
        (
            ("--postmerger-dt", "0.01"),
            ("--postmerger-dt", "0.0"),
            "postmerger-dt",
        ),
        (
            ("--max-relative-dxdt-v-error", "0.01"),
            ("--max-relative-dxdt-v-error", "0.0"),
            "max-relative-dxdt-v-error",
        ),
    ],
)
def test_invalid_physical_parameters_are_rejected(tmp_path, old, new, message):
    input_path = tmp_path / "inspiral.dat"
    _write_input(input_path)
    command = _command(input_path, tmp_path / "invalid.dat")
    position = next(
        index
        for index in range(len(command))
        if tuple(command[index:index + len(old)]) == old
    )
    command[position:position + len(old)] = new
    result = _run(command, success=False)
    assert message in result.stderr


@pytest.mark.parametrize("defect", ["schema", "time", "nonfinite", "superextremal"])
def test_invalid_input_tables_are_rejected(tmp_path, defect):
    table = _inspiral_table()
    if defect == "schema":
        table = table[:, :-1]
        message = "exactly 21"
    elif defect == "time":
        table[20, 0] = table[19, 0]
        message = "strictly increasing"
    elif defect == "nonfinite":
        table[20, 1] = np.nan
        message = "non-finite"
    else:
        table[20, 13] = 0.6
        message = "superextremal"
    input_path = tmp_path / f"{defect}.dat"
    _write_input(input_path, table)
    result = _run(_command(input_path, tmp_path / "invalid.dat"), success=False)
    assert message in result.stderr


def test_measured_velocity_error_can_be_strictly_rejected(tmp_path):
    input_path = tmp_path / "inspiral.dat"
    _write_input(input_path)
    command = _command(input_path, tmp_path / "strict.dat")
    threshold = command.index("--max-relative-dxdt-v-error") + 1
    command[threshold] = "1.0e-4"
    result = _run(command, success=False)
    assert "full-table relative |dx/dt-v| error" in result.stderr


@pytest.mark.parametrize(
    ("missing_option", "message"),
    [
        ("--source-model-label", "--source-model-label"),
        ("--source-revision", "--source-revision"),
    ],
)
def test_source_declarations_are_required(tmp_path, missing_option, message):
    input_path = tmp_path / "inspiral.dat"
    _write_input(input_path)
    output_path = tmp_path / f"missing-{missing_option[2:]}.dat"
    command = _command(input_path, output_path)
    option = command.index(missing_option)
    del command[option:option + 2]
    result = _run(command, success=False)
    assert message in result.stderr
    assert not output_path.exists()
    assert not Path(str(output_path) + ".provenance.json").exists()


@pytest.mark.parametrize(
    ("source_labels", "source_revision", "message"),
    [
        (("   ",), "fixture-v1", "source-model-label"),
        (("duplicate", "duplicate"), "fixture-v1", "must be unique"),
        (("synthetic",), "   ", "source-revision"),
    ],
)
def test_source_declarations_must_be_auditable(
    tmp_path, source_labels, source_revision, message
):
    input_path = tmp_path / "inspiral.dat"
    _write_input(input_path)
    output_path = tmp_path / "invalid-source.dat"
    result = _run(
        _command(
            input_path,
            output_path,
            source_labels=source_labels,
            source_revision=source_revision,
        ),
        success=False,
    )
    assert message in result.stderr
    assert not output_path.exists()


def test_missing_source_provenance_is_rejected(tmp_path):
    input_path = tmp_path / "inspiral.dat"
    _write_input(input_path)
    output_path = tmp_path / "missing-source-provenance.dat"
    result = _run(
        _command(
            input_path,
            output_path,
            source_provenance=tmp_path / "does-not-exist.json",
        ),
        success=False,
    )
    assert "source-provenance is not a readable file" in result.stderr
    assert not output_path.exists()


def test_force_transaction_restores_both_old_files_on_second_replace_failure(
    tmp_path, monkeypatch
):
    """A failed provenance replace restores the byte-exact pre-call pair."""
    table_target = tmp_path / "table.dat"
    provenance_target = tmp_path / "table.dat.provenance.json"
    table_target.write_bytes(b"old table\n")
    provenance_target.write_bytes(b'{"old": "provenance"}\n')
    staged_table = tmp_path / ".new-table.tmp"
    staged_provenance = tmp_path / ".new-provenance.tmp"
    staged_table.write_bytes(b"new table\n")
    staged_provenance.write_bytes(b'{"new": "provenance"}\n')

    real_replace = STITCHER.os.replace

    def fail_provenance_replace(source, target):
        if Path(source) == staged_provenance:
            raise OSError("injected second replace failure")
        return real_replace(source, target)

    monkeypatch.setattr(STITCHER.os, "replace", fail_provenance_replace)
    with pytest.raises(OSError, match="injected second replace failure"):
        STITCHER._publish_transaction(
            (
                (staged_table, table_target),
                (staged_provenance, provenance_target),
            ),
            force=True,
        )

    assert table_target.read_bytes() == b"old table\n"
    assert provenance_target.read_bytes() == b'{"old": "provenance"}\n'
    assert not list(tmp_path.glob(".*.rollback"))


def test_no_clobber_transaction_removes_partial_pair_on_second_link_failure(
    tmp_path, monkeypatch
):
    """A no-clobber transaction never leaves only its first new artifact."""
    table_target = tmp_path / "table.dat"
    provenance_target = tmp_path / "table.dat.provenance.json"
    staged_table = tmp_path / ".new-table.tmp"
    staged_provenance = tmp_path / ".new-provenance.tmp"
    staged_table.write_bytes(b"new table\n")
    staged_provenance.write_bytes(b'{"new": "provenance"}\n')
    real_link = STITCHER.os.link

    def fail_provenance_link(source, target):
        if Path(source) == staged_provenance:
            raise OSError("injected second link failure")
        return real_link(source, target)

    monkeypatch.setattr(STITCHER.os, "link", fail_provenance_link)
    with pytest.raises(OSError, match="injected second link failure"):
        STITCHER._publish_transaction(
            (
                (staged_table, table_target),
                (staged_provenance, provenance_target),
            ),
            force=False,
        )

    assert not table_target.exists()
    assert not provenance_target.exists()
