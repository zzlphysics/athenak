"""Fail-closed unit tests for the effective-BBH 4PN campaign generator."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys

import pytest


REPOSITORY = Path(__file__).resolve().parents[3]
GENERATOR_PATH = (
    REPOSITORY / "inputs" / "dyngr" / "generate_effective_bbh_4pn_campaign.py"
)
MATRIX_PATH = REPOSITORY / "inputs" / "dyngr" / "effective_bbh_4pn_campaign.json"
MODULE_SPEC = importlib.util.spec_from_file_location(
    "effective_bbh_4pn_campaign_generator", GENERATOR_PATH
)
assert MODULE_SPEC is not None and MODULE_SPEC.loader is not None
GENERATOR = importlib.util.module_from_spec(MODULE_SPEC)
sys.modules[MODULE_SPEC.name] = GENERATOR
MODULE_SPEC.loader.exec_module(GENERATOR)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _state(
    position1: tuple[float, float, float],
    position2: tuple[float, float, float],
    *,
    mass1: float = 0.5,
    mass2: float = 0.5,
) -> tuple[float, ...]:
    zeros = (0.0, 0.0, 0.0)
    return position1 + position2 + zeros + zeros + zeros + zeros + (mass1, mass2)


def _write_trajectory(
    path: Path, rows: list[tuple[float, tuple[float, ...]]]
) -> None:
    text = "".join(
        " ".join(f"{value:.17e}" for value in (time, *state)) + "\n"
        for time, state in rows
    )
    path.write_text(text, encoding="ascii")


def test_interval_horizon_guard_captures_hidden_hermite_speed(tmp_path):
    """A 0.99 middle control inflates the certified guard despite zero row speeds."""
    trajectory = tmp_path / "hidden-speed.dat"
    displacement = 0.99 / 3.0
    rows = [
        (0.0, _state((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0))),
        (
            1.0,
            _state((-1.0 + displacement, 0.0, 0.0), (1.0, 0.0, 0.0)),
        ),
    ]
    _write_trajectory(trajectory, rows)

    factor = 2.5
    summary = GENERATOR.inspect_trajectory(
        trajectory, start_time=0.5, horizon_factor=factor, guard_limit=100.0
    )
    expected_row_guard = factor
    expected_interval_guard = factor / math.sqrt(1.0 - 0.99**2)
    assert summary.maximum_speed == pytest.approx(0.0, abs=0.0)
    assert summary.maximum_middle_control_speed == pytest.approx(0.99)
    assert summary.maximum_row_guard_radius == pytest.approx(expected_row_guard)
    assert summary.maximum_guard_radius == pytest.approx(expected_interval_guard)
    assert summary.maximum_guard_radius > 7.0 * summary.maximum_row_guard_radius


@pytest.mark.parametrize("middle_speed", [1.0, 1.01])
def test_non_subluminal_hermite_middle_control_is_rejected(
    tmp_path, middle_speed
):
    trajectory = tmp_path / f"superluminal-{middle_speed}.dat"
    # A three-unit interval makes the Hermite middle control exactly equal to
    # this displacement when both endpoint velocities vanish.  In particular,
    # the 1.0 case uses only exactly representable integer arithmetic.
    displacement = middle_speed
    _write_trajectory(
        trajectory,
        [
            (0.0, _state((-1.0, 0.0, 0.0), (1.0, 0.0, 0.0))),
            (
                3.0,
                _state(
                    (-1.0 + displacement, 0.0, 0.0),
                    (1.0, 0.0, 0.0),
                ),
            ),
        ],
    )
    with pytest.raises(
        GENERATOR.CampaignError,
        match="non-subluminal Hermite middle control",
    ):
        GENERATOR.inspect_trajectory(
            trajectory, start_time=1.5, horizon_factor=1.25, guard_limit=100.0
        )


def _write_synthetic_campaign_artifacts(tmp_path: Path) -> tuple[Path, ...]:
    """Create an internally hashed minimal campaign accepted by the CLI gates."""
    trajectory = tmp_path / "synthetic-campaign.dat"
    rows = [
        (1.0, _state((-10.0, 0.0, 0.0), (10.0, 0.0, 0.0))),
        (100.0, _state((-10.0, 0.0, 0.0), (10.0, 0.0, 0.0))),
        (
            201.0,
            _state(
                (0.0, 0.0, 0.0),
                (0.0, 0.0, 0.0),
                mass1=0.95,
                mass2=0.0,
            ),
        ),
    ]
    _write_trajectory(trajectory, rows)

    source_artifact = tmp_path / "synthetic-source.txt"
    source_artifact.write_text("synthetic source revision fixture\n", encoding="ascii")
    revision = "1" * 40
    matrix = json.loads(MATRIX_PATH.read_text(encoding="utf-8"))
    provenance = {
        "output": {
            "path": str(trajectory.resolve()),
            "rows": len(rows),
            "sha256": _sha256(trajectory),
            "time_bounds": [rows[0][0], rows[-1][0]],
        },
        "provenance_schema_version": 2,
        "remnant": {
            "a": [0.0, 0.0, 0.0],
            "a_to_chi_relation": "a_length = Mf * chi; chi = a_length / Mf",
            "chi": [0.0, 0.0, 0.0],
            "kick": [0.0, 0.0, 0.0],
            "mass": 0.95,
            "representation": "canonical-single-term-1",
        },
        "schema": {
            "columns": matrix["trajectory"]["columns"],
            "version": 1,
        },
        "trajectory_model": {
            "declaration_source": "CLI",
            "labels": ["frozen-CBwaves", "local-instantaneous-4PN"],
            "source_provenance": {
                "path": str(source_artifact.resolve()),
                "sha256": _sha256(source_artifact),
            },
            "source_revision": revision,
        },
    }
    provenance_path = tmp_path / "synthetic-campaign.provenance.json"
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    matrix["trajectory"]["expected_ascii_sha256"] = _sha256(trajectory)
    matrix["trajectory"]["expected_provenance_sha256"] = _sha256(
        provenance_path
    )
    matrix["trajectory"]["expected_source_revision"] = revision
    matrix["trajectory"]["expected_source_sha256"] = _sha256(source_artifact)
    custom_matrix = tmp_path / "synthetic-campaign-matrix.json"
    custom_matrix.write_text(
        json.dumps(matrix, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return trajectory, provenance_path, source_artifact, custom_matrix


def test_default_offset_covers_the_first_metric_fd_stencil(tmp_path):
    """The CLI default leaves a complete centered metric stencil at t=0."""
    trajectory, provenance, source_artifact, matrix_path = (
        _write_synthetic_campaign_artifacts(tmp_path)
    )
    command = [
        sys.executable,
        str(GENERATOR_PATH),
        "L",
        "--trajectory",
        str(trajectory),
        "--trajectory-provenance",
        str(provenance),
        "--source-artifact",
        str(source_artifact),
        "--matrix",
        str(matrix_path),
        "--gpus",
        "4",
        "--scratch-gib",
        "10000",
        "--validate-only",
    ]
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        pytest.fail(
            "synthetic campaign validation failed:\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    report = json.loads(result.stdout)
    matrix = json.loads(matrix_path.read_text(encoding="utf-8"))
    first_time = report["trajectory_time_range"][0]
    fd_step = matrix["common"]["time"]["metric_fd_step_M"]
    expected_offset = math.nextafter(first_time + fd_step, math.inf)
    actual_offset = report["trajectory_time_offset_M"]
    assert actual_offset == expected_offset
    assert first_time <= actual_offset - fd_step
    assert report["trajectory_time_range"][1] >= (
        actual_offset + report["tlim_M"] + fd_step
    )
