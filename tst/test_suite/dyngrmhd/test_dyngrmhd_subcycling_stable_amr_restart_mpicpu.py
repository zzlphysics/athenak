"""Stable-AMR restart equivalence for strict 2:1 dynGRMHD subcycling."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import bin_convert
import test_suite.testutils as testutils


ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import compare_athenak_restart_fields as restart_fields  # noqa: E402
import read_athenak_restart_metadata as restart_metadata  # noqa: E402


INPUT = ROOT / "inputs" / "dyngr" / "subcycling_stable_amr_restart_dyngrmhd.athinput"
BASENAME = "subcycling_stable_amr_restart_dyngrmhd"
RANKS = 2
SPLIT_CYCLE = 3
FINAL_CYCLE = 6
NGHOST = 4
TIMEOUT_SECONDS = 30
# A synchronized checkpoint is an exact continuation boundary.  Keep this regression
# bitwise: a tolerance would hide another uninitialized/stale restart cache.
RTOL = 0.0
ATOL = 0.0


def _mpi_launcher() -> str:
    launcher = shutil.which("mpirun") or shutil.which("mpiexec")
    if launcher is not None:
        return launcher
    cache = Path("..").resolve() / "CMakeCache.txt"
    if cache.is_file():
        for line in cache.read_text(encoding="utf-8").splitlines():
            if line.startswith("MPI_CXX_COMPILER:FILEPATH="):
                compiler = Path(line.split("=", 1)[1])
                for name in ("mpirun", "mpiexec"):
                    candidate = compiler.with_name(name)
                    if candidate.is_file():
                        return str(candidate)
    raise AssertionError("an MPI launcher matching the test binary is required")


def _run(
    run_dir: Path,
    *,
    nlim: int,
    restart: Path | None = None,
    overrides: tuple[str, ...] = (),
) -> str:
    command = [_mpi_launcher(), "-np", str(RANKS), "./athena"]
    if restart is None:
        command += ["-i", str(INPUT)]
    else:
        command += ["-r", str(restart.resolve())]
    command += [
        "-d",
        str(run_dir),
        f"job/basename={BASENAME}",
        f"time/nlim={nlim}",
        "time/tlim=0.1",
    ]
    command += list(overrides)

    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=TIMEOUT_SECONDS,
    )
    with Path(testutils.LOG_FILE_PATH).open("a", encoding="utf-8") as log_file:
        label = "restart" if restart is not None else "fresh"
        log_file.write(
            f"\nStable-AMR subcycling {label} command: {' '.join(command)}\n"
        )
        log_file.write(result.stdout)
    assert result.returncode == 0, result.stdout[-8000:]
    return result.stdout


def _run_expect_failure(
    run_dir: Path,
    *,
    restart: Path,
    expected: str,
    overrides: tuple[str, ...] = (),
) -> str:
    command = [
        _mpi_launcher(),
        "-np",
        str(RANKS),
        "./athena",
        "-r",
        str(restart.resolve()),
        "-d",
        str(run_dir),
        f"job/basename={BASENAME}",
        f"time/nlim={SPLIT_CYCLE}",
        "time/tlim=0.1",
        *overrides,
    ]
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=TIMEOUT_SECONDS,
    )
    with Path(testutils.LOG_FILE_PATH).open("a", encoding="utf-8") as log_file:
        log_file.write(
            f"\nExpected stable-AMR restart failure: {' '.join(command)}\n"
        )
        log_file.write(result.stdout)
    assert result.returncode != 0, result.stdout[-8000:]
    assert expected in result.stdout, result.stdout[-8000:]
    return result.stdout


def _latest_dump(run_dir: Path, variable: str) -> dict[str, object]:
    paths = sorted((run_dir / "bin").glob(f"{BASENAME}.{variable}.*.bin"))
    assert paths, f"no {variable} output found in {run_dir}"
    dumps = [bin_convert.read_binary(str(path)) for path in paths]
    return max(dumps, key=lambda dump: (dump["cycle"], dump["time"]))


def _latest_restart(run_dir: Path) -> tuple[Path, object]:
    paths = sorted((run_dir / "rst").glob(f"{BASENAME}.*.rst"))
    assert paths, f"no restart output found in {run_dir}"
    candidates = [(path, restart_metadata.read_restart_metadata(path)) for path in paths]
    return max(candidates, key=lambda item: (item[1].cycle, item[1].time))


def _canonical_order(dump: dict[str, object]) -> np.ndarray:
    logical = np.asarray(dump["mb_logical"])
    return np.lexsort((logical[:, 0], logical[:, 1], logical[:, 2], logical[:, 3]))


def _canonical_field(dump: dict[str, object], variable: str) -> np.ndarray:
    values = np.stack(dump["mb_data"][variable], axis=0)
    return values[_canonical_order(dump)]


def _assert_same_topology(
    reference: dict[str, object], candidate: dict[str, object]
) -> None:
    for name in ("Nx1", "Nx2", "Nx3", "nx1_mb", "nx2_mb", "nx3_mb", "n_mbs"):
        assert candidate[name] == reference[name], f"topology differs in {name}"

    reference_order = _canonical_order(reference)
    candidate_order = _canonical_order(candidate)
    for name in ("mb_logical", "mb_index", "mb_geometry"):
        np.testing.assert_array_equal(
            np.asarray(candidate[name])[candidate_order],
            np.asarray(reference[name])[reference_order],
            err_msg=f"topology differs in {name}",
        )


def _assert_active_and_ghost_close(
    reference: dict[str, object],
    candidate: dict[str, object],
    *,
    description: str,
    final_cycle: int,
) -> None:
    _assert_same_topology(reference, candidate)
    assert candidate["cycle"] == reference["cycle"] == final_cycle
    assert candidate["time"] == pytest.approx(reference["time"], rel=2.0e-14)
    assert candidate["var_names"] == reference["var_names"]

    for variable in reference["var_names"]:
        expected = _canonical_field(reference, variable)
        actual = _canonical_field(candidate, variable)
        assert expected.ndim == 4, f"unexpected {variable} field shape {expected.shape}"
        assert expected.shape[-3:] == (
            reference["nx3_mb"] + 2 * NGHOST,
            reference["nx2_mb"] + 2 * NGHOST,
            reference["nx1_mb"] + 2 * NGHOST,
        )
        assert np.isfinite(expected).all(), f"non-finite reference {variable}"
        assert np.isfinite(actual).all(), f"non-finite resumed {variable}"

        active = (
            slice(None),
            slice(NGHOST, -NGHOST),
            slice(NGHOST, -NGHOST),
            slice(NGHOST, -NGHOST),
        )
        np.testing.assert_allclose(
            actual[active],
            expected[active],
            rtol=RTOL,
            atol=ATOL,
            err_msg=f"{description} changed active {variable}",
        )

        ghost_mask = np.ones(expected.shape[-3:], dtype=bool)
        ghost_mask[
            NGHOST:-NGHOST, NGHOST:-NGHOST, NGHOST:-NGHOST
        ] = False
        np.testing.assert_allclose(
            actual[..., ghost_mask],
            expected[..., ghost_mask],
            rtol=RTOL,
            atol=ATOL,
            err_msg=f"{description} changed ghost {variable}",
        )


def _assert_equivalent_outputs(
    continuous_dir: Path,
    resumed_dir: Path,
    *,
    final_cycle: int,
    expected_level_counts: dict[int, int],
) -> None:
    continuous_u = _latest_dump(continuous_dir, "mhd_u_bcc")
    resumed_u = _latest_dump(resumed_dir, "mhd_u_bcc")
    continuous_w = _latest_dump(continuous_dir, "mhd_w_bcc")
    resumed_w = _latest_dump(resumed_dir, "mhd_w_bcc")
    continuous_divb = _latest_dump(continuous_dir, "mhd_divb")
    resumed_divb = _latest_dump(resumed_dir, "mhd_divb")

    levels = np.asarray(continuous_u["mb_logical"])[:, 3].astype(int)
    levels -= int(np.min(levels))
    unique, counts = np.unique(levels, return_counts=True)
    assert dict(zip(unique.tolist(), counts.tolist())) == expected_level_counts

    _assert_active_and_ghost_close(
        continuous_u,
        resumed_u,
        description="split restart conserved state",
        final_cycle=final_cycle,
    )
    _assert_active_and_ghost_close(
        continuous_w,
        resumed_w,
        description="split restart primitive state",
        final_cycle=final_cycle,
    )

    _assert_same_topology(continuous_divb, resumed_divb)
    continuous_divb_values = _canonical_field(continuous_divb, "divb")
    resumed_divb_values = _canonical_field(resumed_divb, "divb")
    np.testing.assert_allclose(
        resumed_divb_values,
        continuous_divb_values,
        rtol=0.0,
        atol=2.0e-13,
        err_msg="split restart changed div(B)",
    )
    assert float(np.max(np.abs(continuous_divb_values))) < 1.0e-11
    assert float(np.max(np.abs(resumed_divb_values))) < 1.0e-11

    continuous_restart, continuous_metadata = _latest_restart(continuous_dir)
    resumed_restart, resumed_metadata = _latest_restart(resumed_dir)
    assert continuous_metadata.cycle == resumed_metadata.cycle == final_cycle
    assert (
        continuous_metadata.restart_cache_contract_version
        == resumed_metadata.restart_cache_contract_version
        == restart_metadata.RESTART_CACHE_CONTRACT_VERSION
    )
    assert (
        continuous_metadata.event_counter_version
        == resumed_metadata.event_counter_version
        == restart_metadata.EVENT_COUNTER_VERSION
    )
    assert continuous_metadata.event_counters is not None
    assert resumed_metadata.event_counters == continuous_metadata.event_counters
    comparison = restart_fields.compare_restart_fields(
        continuous_restart, resumed_restart, rtol=RTOL, atol=ATOL
    )
    assert comparison["authoritative_match"], json.dumps(comparison, indent=2)
    assert comparison["all_stored_fields_match"], json.dumps(comparison, indent=2)


def test_stable_amr_continuous_and_split_restarts_are_equivalent(
    tmp_path: Path,
) -> None:
    """A synchronized restart must not change any state consumed by later steps."""
    continuous_dir = tmp_path / "continuous"
    split_dir = tmp_path / "split"
    resumed_dir = tmp_path / "resumed"

    continuous_log = _run(continuous_dir, nlim=FINAL_CYCLE)
    split_log = _run(split_dir, nlim=SPLIT_CYCLE)
    split_restart, split_metadata = _latest_restart(split_dir)
    assert split_metadata.cycle == SPLIT_CYCLE
    assert (
        split_metadata.restart_cache_contract_version
        == restart_metadata.RESTART_CACHE_CONTRACT_VERSION
    )
    resumed_log = _run(resumed_dir, nlim=FINAL_CYCLE, restart=split_restart)

    # Tracking creates one octet at cycle 1 and then leaves the hierarchy unchanged.
    # This distinguishes a restart-state failure from regrid motion.
    assert "7 MeshBlocks created, 0 deleted by AMR" in continuous_log
    assert "7 MeshBlocks created, 0 deleted by AMR" in split_log
    assert "0 MeshBlocks created, 0 deleted by AMR" in resumed_log

    _assert_equivalent_outputs(
        continuous_dir,
        resumed_dir,
        final_cycle=FINAL_CYCLE,
        expected_level_counts={0: 2, 1: 8},
    )


def test_three_physical_levels_are_bitwise_restart_equivalent(tmp_path: Path) -> None:
    """Exercise recursive subcycling at physical-boundary/coarse-fine corners."""
    continuous_dir = tmp_path / "continuous_three_levels"
    split_dir = tmp_path / "split_three_levels"
    resumed_dir = tmp_path / "resumed_three_levels"
    split_cycle = 5
    final_cycle = 10
    overrides = (
        "mesh_refinement/num_levels=3",
        "mesh_refinement/max_nmb_per_rank=64",
        "output4/dcycle=5",
    )

    continuous_log = _run(
        continuous_dir, nlim=final_cycle, overrides=overrides
    )
    split_log = _run(split_dir, nlim=split_cycle, overrides=overrides)
    split_restart, split_metadata = _latest_restart(split_dir)
    assert split_metadata.cycle == split_cycle
    assert (
        split_metadata.restart_cache_contract_version
        == restart_metadata.RESTART_CACHE_CONTRACT_VERSION
    )
    resumed_log = _run(
        resumed_dir,
        nlim=final_cycle,
        restart=split_restart,
        overrides=overrides,
    )

    assert "77 MeshBlocks created, 0 deleted by AMR" in continuous_log
    assert "77 MeshBlocks created, 0 deleted by AMR" in split_log
    assert "0 MeshBlocks created, 0 deleted by AMR" in resumed_log

    # The root block is fully covered by its children, so the leaf list contains the
    # two descendant levels even though the recursive scheduler spans three levels.
    _assert_equivalent_outputs(
        continuous_dir,
        resumed_dir,
        final_cycle=final_cycle,
        expected_level_counts={0: 16, 1: 64},
    )


def test_legacy_cache_contract_requires_one_time_audited_migration(
    tmp_path: Path,
) -> None:
    """A markerless checkpoint is never accepted silently by strict subcycling."""
    source_dir = tmp_path / "source"
    _run(source_dir, nlim=SPLIT_CYCLE)
    checkpoint, metadata = _latest_restart(source_dir)
    assert metadata.restart_cache_contract_version == 1

    contract_size = 8 + 4
    contract_start = metadata.metadata_end - contract_size
    contents = checkpoint.read_bytes()
    expected_magic = restart_metadata.RESTART_CACHE_CONTRACT_MAGIC.to_bytes(
        8, byteorder=sys.byteorder
    )
    assert contents[contract_start:contract_start + 8] == expected_magic
    legacy_checkpoint = tmp_path / "legacy_without_cache_contract.rst"
    legacy_checkpoint.write_bytes(
        contents[:contract_start] + contents[metadata.metadata_end:]
    )

    _run_expect_failure(
        tmp_path / "legacy_rejected",
        restart=legacy_checkpoint,
        expected="allow_legacy_subcycling_restart=true",
    )

    migrated_dir = tmp_path / "legacy_migrated"
    stdout = _run(
        migrated_dir,
        nlim=SPLIT_CYCLE + 1,
        restart=legacy_checkpoint,
        overrides=("time/allow_legacy_subcycling_restart=true",),
    )
    assert "Legacy strict-subcycling restart qualification passed" in stdout
    _, migrated_metadata = _latest_restart(migrated_dir)
    assert (
        migrated_metadata.restart_cache_contract_version
        == restart_metadata.RESTART_CACHE_CONTRACT_VERSION
    )
    assert migrated_metadata.parameters["time"][
        "allow_legacy_subcycling_restart"
    ] in {"0", "false"}
