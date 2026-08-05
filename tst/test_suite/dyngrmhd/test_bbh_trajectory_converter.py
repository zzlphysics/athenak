"""Artifact-contract tests for the AnalyticalBBH trajectory converter."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import h5py
import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[3]
CONVERTER = ROOT / "scripts" / "bbh_trajectory_h5_to_ascii.py"
MODULE_SPEC = importlib.util.spec_from_file_location(
    "bbh_trajectory_converter", CONVERTER
)
assert MODULE_SPEC is not None and MODULE_SPEC.loader is not None
CONVERTER_MODULE = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(CONVERTER_MODULE)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_q1_trajectory(
    path: Path,
    *,
    bad_velocity=False,
    bad_remnant=False,
    unsafe_final_control=False,
):
    """Write a small, exactly differentiated q=1 nonspinning trajectory."""
    time = np.arange(5, dtype=np.float64)
    radial_coordinate = 1.0e-2 * (4.0 - time) ** 2
    radial_velocity = -2.0e-2 * (4.0 - time)
    zero = np.zeros_like(time)
    spin_z = 0.1 * time / time[-1]

    datasets = {
        "t": time,
        "x1": radial_coordinate,
        "y1": zero,
        "z1": zero,
        "x2": -radial_coordinate,
        "y2": zero,
        "z2": zero,
        "vx1": zero if bad_velocity else radial_velocity,
        "vy1": zero,
        "vz1": zero,
        "vx2": zero if bad_velocity else -radial_velocity,
        "vy2": zero,
        "vz2": zero,
        "a1x": zero,
        "a1y": zero,
        "a1z": spin_z,
        "a2x": zero,
        "a2y": zero,
        "a2z": spin_z.copy(),
        "m1_full": np.linspace(0.5, 0.475, time.size),
        "m2_full": np.linspace(0.5, 0.475, time.size),
    }
    if bad_remnant:
        datasets["a2z"][-1] += 0.05
    if unsafe_final_control:
        linear_position = 0.5 * time
        linear_velocity = np.full_like(time, 0.5)
        linear_velocity[-1] = -0.5
        datasets["x1"] = linear_position
        datasets["x2"] = -linear_position
        datasets["vx1"] = linear_velocity
        datasets["vx2"] = -linear_velocity

    with h5py.File(path, "w") as destination:
        for name, values in datasets.items():
            destination.create_dataset(name, data=values)


def _run(input_path: Path, output_path: Path, *options: str):
    return subprocess.run(
        [sys.executable, str(CONVERTER), str(input_path), str(output_path), *options],
        check=False,
        capture_output=True,
        text=True,
    )


def test_converter_writes_complete_manifest_and_passes_q1_profile(tmp_path):
    input_path = tmp_path / "trajectory.h5"
    output_path = tmp_path / "trajectory.dat"
    manifest_path = tmp_path / "contract.json"
    _write_q1_trajectory(input_path)

    result = _run(
        input_path,
        output_path,
        "--manifest",
        str(manifest_path),
        "--strict-velocity-mismatch",
        "--velocity-tolerance",
        "1e-12",
        "--validation-profile",
        "q1-nonspinning",
    )

    assert result.returncode == 0, result.stderr
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["schema_version"] == 1
    assert manifest["row_count"] == 5
    assert manifest["input_sha256"] == _sha256(input_path)
    assert manifest["output_sha256"] == _sha256(output_path)
    assert manifest["time_range"] == {
        "start": 0.0,
        "end": 4.0,
        "unit": "normalized mass time M",
    }
    assert manifest["mass_unit"]["input_value"] == 1.0
    assert manifest["initial_state"]["masses"] == [0.5, 0.5]
    assert manifest["initial_state"]["separation"] == 0.32
    assert manifest["initial_state"]["speeds"] == [0.08, 0.08]
    assert manifest["initial_state"]["spin_a"] == [[0.0, 0.0, 0.0]] * 2
    assert manifest["final_state"]["masses"] == [0.475, 0.475]
    assert manifest["final_state"]["separation"] == 0.0
    assert manifest["final_state"]["speeds"] == [0.0, 0.0]
    assert manifest["final_state"]["spin_a"] == [[0.0, 0.0, 0.1]] * 2
    assert manifest["velocity_mismatch"]["overall_max"]["row_index"] in (1, 2, 3)
    assert manifest["velocity_mismatch"]["overall_max"]["max_fractional"] < 1e-12
    certificate = manifest["hermite_velocity_certificate"]
    assert certificate["formula_version"] == (
        "dynbbh-quadratic-bezier-controls-v1"
    )
    assert certificate["certified_subluminal"] is True
    assert certificate["endpoint_max"]["speed"] == 0.08
    assert certificate["middle_control_max"]["interval_index"] == 0
    assert certificate["middle_control_max"]["speed"] == pytest.approx(0.07)
    assert manifest["validation"]["remnant_representation"] == (
        "duplicate-identical-terms"
    )


def test_converter_output_and_default_manifest_are_byte_stable(tmp_path):
    input_path = tmp_path / "trajectory.h5"
    _write_q1_trajectory(input_path)
    outputs = [tmp_path / "first.dat", tmp_path / "second.dat"]

    results = [_run(input_path, output) for output in outputs]

    assert all(result.returncode == 0 for result in results)
    manifests = [Path(f"{output}.manifest.json") for output in outputs]
    assert outputs[0].read_bytes() == outputs[1].read_bytes()
    assert manifests[0].read_bytes() == manifests[1].read_bytes()
    assert json.loads(manifests[0].read_text())["output_sha256"] == _sha256(outputs[0])


def test_strict_velocity_mismatch_rejects_without_artifacts(tmp_path):
    input_path = tmp_path / "bad-velocity.h5"
    output_path = tmp_path / "bad-velocity.dat"
    manifest_path = Path(f"{output_path}.manifest.json")
    _write_q1_trajectory(input_path, bad_velocity=True)

    result = _run(
        input_path,
        output_path,
        "--strict-velocity-mismatch",
        "--velocity-tolerance",
        "1e-6",
    )

    assert result.returncode != 0
    assert "exceeds tolerance" in result.stderr
    assert not output_path.exists()
    assert not manifest_path.exists()
    assert not list(tmp_path.glob(".*.tmp"))


@pytest.mark.parametrize("strict", [False, True])
def test_unsafe_final_hermite_control_is_always_rejected(tmp_path, strict):
    """Interior dx/dt checks cannot certify the final Hermite interval."""
    input_path = tmp_path / f"unsafe-control-{strict}.h5"
    output_path = tmp_path / f"unsafe-control-{strict}.dat"
    manifest_path = Path(f"{output_path}.manifest.json")
    _write_q1_trajectory(input_path, unsafe_final_control=True)
    options = (
        ("--strict-velocity-mismatch", "--velocity-tolerance", "1e-12")
        if strict
        else ()
    )

    result = _run(input_path, output_path, *options)

    assert result.returncode != 0
    assert "middle velocity control has |v|^2=2.25" in result.stderr
    assert not output_path.exists()
    assert not manifest_path.exists()
    assert not list(tmp_path.glob(".*.tmp"))


def test_q1_profile_rejects_ambiguous_remnant_without_artifacts(tmp_path):
    input_path = tmp_path / "bad-remnant.h5"
    output_path = tmp_path / "bad-remnant.dat"
    manifest_path = Path(f"{output_path}.manifest.json")
    _write_q1_trajectory(input_path, bad_remnant=True)

    result = _run(
        input_path,
        output_path,
        "--validation-profile",
        "q1-nonspinning",
    )

    assert result.returncode != 0
    assert "final remnant is neither" in result.stderr
    assert not output_path.exists()
    assert not manifest_path.exists()
    assert not list(tmp_path.glob(".*.tmp"))


def test_no_force_publish_preserves_concurrent_target_and_rolls_back_pair(
    tmp_path, monkeypatch
):
    """A target created after preflight is never overwritten or rolled back."""
    output_target = tmp_path / "trajectory.dat"
    manifest_target = tmp_path / "trajectory.dat.manifest.json"
    staged_output = tmp_path / ".trajectory.tmp"
    staged_manifest = tmp_path / ".manifest.tmp"
    staged_output.write_bytes(b"new trajectory\n")
    staged_manifest.write_bytes(b'{"new": "manifest"}\n')
    concurrent_manifest = b'{"user": "concurrent creation"}\n'
    real_link = CONVERTER_MODULE.os.link

    def create_manifest_before_link(source, target):
        if Path(target) == manifest_target:
            manifest_target.write_bytes(concurrent_manifest)
        return real_link(source, target)

    monkeypatch.setattr(CONVERTER_MODULE.os, "link", create_manifest_before_link)
    with pytest.raises(FileExistsError, match="refusing to replace"):
        CONVERTER_MODULE._publish_staged(
            [
                (staged_output, output_target),
                (staged_manifest, manifest_target),
            ],
            force=False,
        )

    assert not output_target.exists()
    assert manifest_target.read_bytes() == concurrent_manifest
    assert not staged_output.exists()
    assert not staged_manifest.exists()
    assert not list(tmp_path.glob(".*.rollback"))
