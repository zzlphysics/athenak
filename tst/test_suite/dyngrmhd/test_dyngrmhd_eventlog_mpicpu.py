"""DynGRMHD event telemetry, MPI aggregation, and error-print cap regressions."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

import test_suite.testutils as testutils
import bin_convert


ROOT = Path(__file__).resolve().parents[3]
INPUT = ROOT / "tst" / "inputs" / "dyngrmhd_eventlog.athinput"
TIMEOUT_SECONDS = 30


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


def _run(run_dir: Path, *, ranks: int, overrides: tuple[str, ...] = ()) -> str:
    command = [] if ranks == 1 else [_mpi_launcher(), "-np", str(ranks)]
    command += [
        "./athena",
        "-i",
        str(INPUT),
        "-d",
        str(run_dir),
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
        log_file.write(f"\nDynGRMHD event-log command: {' '.join(command)}\n")
        log_file.write(result.stdout)
    assert result.returncode == 0, result.stdout[-8000:]
    return result.stdout


def _run_restart(
        restart: Path, run_dir: Path, *, nlim: int, ranks: int = 1,
        overrides: tuple[str, ...] = ()) -> str:
    command = [] if ranks == 1 else [_mpi_launcher(), "-np", str(ranks)]
    command += [
        "./athena",
        "-r",
        str(restart),
        "-d",
        str(run_dir),
        f"time/nlim={nlim}",
        "time/tlim=100",
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
        log_file.write(f"\nDynGRMHD event-log restart: {' '.join(command)}\n")
        log_file.write(result.stdout)
    assert result.returncode == 0, result.stdout[-8000:]
    return result.stdout


def _run_restart_expect_failure(
        restart: Path, run_dir: Path, *, ranks: int,
        expected: str, overrides: tuple[str, ...] = ()) -> None:
    command = [
        _mpi_launcher(),
        "-np",
        str(ranks),
        "./athena",
        "-r",
        str(restart),
        "-d",
        str(run_dir),
        "time/nlim=1",
        "time/tlim=100",
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
            f"\nExpected DynGRMHD event-log restart failure: "
            f"{' '.join(command)}\n"
        )
        log_file.write(result.stdout)
    assert result.returncode != 0, result.stdout[-8000:]
    assert expected in result.stdout, result.stdout[-8000:]


def _restart_metadata(restart: Path, *, ranks: int) -> dict[str, object]:
    command = [
        sys.executable,
        str(ROOT / "scripts" / "read_athenak_restart_metadata.py"),
        str(restart),
        "--ranks",
        str(ranks),
        "--json",
    ]
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
        timeout=TIMEOUT_SECONDS,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def _restart_parameter(restart: Path, block: str, name: str) -> str:
    raw = restart.read_bytes()
    marker = b"<par_end>\n"
    parameter_end = raw.find(marker)
    assert parameter_end >= 0, f"No ParameterInput terminator in {restart}"
    text = raw[:parameter_end].decode("utf-8")
    section = re.search(
        rf"(?ms)^<{re.escape(block)}>\s*$\n(.*?)(?=^<[^>]+>\s*$|\Z)", text
    )
    assert section is not None, f"No <{block}> block in {restart}"
    values = re.findall(
        rf"(?m)^\s*{re.escape(name)}\s*=\s*(\S+)", section.group(1)
    )
    assert len(values) == 1, f"Unexpected {block}/{name}: {values}"
    return values[0]


def _run_expect_failure(run_dir: Path, *overrides: str) -> str:
    command = [
        "./athena",
        "-i",
        str(INPUT),
        "-d",
        str(run_dir),
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
        log_file.write(f"\nDynGRMHD rejected-input command: {' '.join(command)}\n")
        log_file.write(result.stdout)
    assert result.returncode != 0, result.stdout[-8000:]
    return result.stdout


def _event_rows(run_dir: Path) -> list[dict[str, int]]:
    lines = (run_dir / "dyngrmhd_eventlog.log").read_text(
        encoding="utf-8"
    ).splitlines()
    headers = [line[1:].split() for line in lines if line.startswith("#  cycle")]
    assert len(headers) == 1
    names = headers[0]
    rows = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        values = [int(value) for value in line.split()]
        assert len(values) == len(names)
        rows.append(dict(zip(names, values)))
    return rows


def _state_dump_at_cycle(run_dir: Path, cycle: int) -> dict[str, object]:
    paths = sorted(
        (run_dir / "bin").glob("dyngrmhd_eventlog.mhd_u_bcc.*.bin")
    )
    assert paths, f"no mhd_u_bcc output found in {run_dir}"
    matches = []
    for path in paths:
        dump = bin_convert.read_binary(str(path))
        if dump["cycle"] == cycle:
            matches.append(dump)
    assert len(matches) == 1, (
        f"expected exactly one mhd_u_bcc output at cycle {cycle}, got "
        f"{len(matches)}"
    )
    return matches[0]


def _canonical_order(dump: dict[str, object]) -> np.ndarray:
    logical = np.asarray(dump["mb_logical"])
    return np.lexsort(
        (logical[:, 0], logical[:, 1], logical[:, 2], logical[:, 3])
    )


def _canonical_field(dump: dict[str, object], variable: str) -> np.ndarray:
    values = np.stack(dump["mb_data"][variable], axis=0)
    return values[_canonical_order(dump)]


def test_dense_event_rows_and_mpi_uint64_aggregation(tmp_path: Path) -> None:
    serial_dir = tmp_path / "serial"
    mpi_dir = tmp_path / "mpi"
    _run(serial_dir, ranks=1)
    _run(
        mpi_dir,
        ranks=2,
        overrides=("mesh/nx1=16", "meshblock/nx1=8"),
    )

    serial = _event_rows(serial_dir)
    mpi = _event_rows(mpi_dir)
    assert [row["cycle"] for row in serial] == [0, 1, 2]
    assert [row["cycle"] for row in mpi] == [0, 1, 2]
    for serial_row, mpi_row in zip(serial, mpi):
        for name in (
            "eos_dfloor",
            "eos_efloor",
            "eos_tfloor",
            "eos_vceil",
            "eos_fail",
            "fofc",
            "cons_adjust",
            "mag_adjust",
            "c2p_calls",
            "fofc_tests",
        ):
            assert mpi_row[name] == 2 * serial_row[name]
        assert mpi_row["c2p_it"] == serial_row["c2p_it"]
    # Telemetry counts physical active zones, not the six ghost cells in this 1D
    # problem.  This keeps the values invariant under ghost width and decomposition.
    assert [row["c2p_calls"] for row in serial] == [8, 16, 16]
    assert [row["fofc_tests"] for row in serial] == [0, 16, 16]


def test_active_zone_events_are_decomposition_invariant(tmp_path: Path) -> None:
    """The same physical mesh must not gain events merely from more ghost zones."""
    one_block_dir = tmp_path / "one_block"
    two_block_dir = tmp_path / "two_blocks"
    global_mesh = ("mesh/nx1=16",)
    _run(
        one_block_dir,
        ranks=1,
        overrides=(*global_mesh, "meshblock/nx1=16"),
    )
    _run(
        two_block_dir,
        ranks=2,
        overrides=(*global_mesh, "meshblock/nx1=8"),
    )

    assert _event_rows(two_block_dir) == _event_rows(one_block_dir)


def test_c2p_print_cap_does_not_truncate_failure_counter(tmp_path: Path) -> None:
    run_dir = tmp_path / "failure_cap"
    stdout = _run(
        run_dir,
        ranks=1,
        overrides=(
            "time/nlim=1",
            "mhd/c2p_iter=1",
            "mhd/c2perrs=2",
            "output2/dcycle=0",
        ),
    )
    assert stdout.count("An error occurred during the primitive solve:") == 2
    assert stdout.count("All future C2P errors") == 1

    rows = _event_rows(run_dir)
    assert [row["cycle"] for row in rows] == [0, 1]
    assert sum(row["eos_fail"] for row in rows) == 24
    assert sum(row["c2p_calls"] for row in rows) == 24
    assert sum(row["cons_adjust"] for row in rows) == 24


def test_restart_segments_have_one_dense_row_per_completed_cycle(tmp_path: Path) -> None:
    first_dir = tmp_path / "first"
    resumed_dir = tmp_path / "resumed"
    _run(first_dir, ranks=1)
    checkpoint = first_dir / "rst" / "dyngrmhd_eventlog.00001.rst"
    assert checkpoint.is_file()
    _run_restart(checkpoint, resumed_dir, nlim=2)

    first_rows = _event_rows(first_dir)
    resumed_rows = _event_rows(resumed_dir)
    # Restart initialization reconstructs primitives and is deliberately counted in the
    # next segment.  The cadence contract here is one row for each completed cycle, with
    # neither a replay of source cycle 1 nor a missing cycle 2.
    assert [row["cycle"] for row in first_rows[:2] + resumed_rows] == [0, 1, 2]
    assert [row["cycle"] for row in resumed_rows] == [2]
    assert resumed_rows[0]["c2p_calls"] > first_rows[2]["c2p_calls"]


def test_zero_step_restart_preserves_physical_main_state(tmp_path: Path) -> None:
    """Restart initialization must not alter checkpointed conserved MHD state."""
    source_dir = tmp_path / "source"
    resumed_dir = tmp_path / "zero_step_restart"
    mesh_overrides = ("mesh/nx1=16", "meshblock/nx1=8")
    _run(
        source_dir,
        ranks=2,
        overrides=(
            *mesh_overrides,
            "time/nlim=0",
            # Force the initial C2P pass through the conserved magnetization
            # adjustment.  The resulting checkpoint is already admissible, so a
            # second C2P pass during restart initialization must be idempotent.
            "problem/b3=1.0",
            "mhd/max_bsq=1.0",
            "output3/dcycle=1",
        ),
    )
    checkpoint = source_dir / "rst" / "dyngrmhd_eventlog.00000.rst"
    assert checkpoint.is_file()
    source_events = _event_rows(source_dir)
    assert [row["cycle"] for row in source_events] == [0]
    assert source_events[0]["mag_adjust"] > 0

    # nlim already equals the checkpoint cycle, so this executes initialization and
    # final output only.  Event/output metadata may legitimately advance; the physical
    # field comparison below deliberately excludes it.
    _run_restart(
        checkpoint,
        resumed_dir,
        nlim=0,
        ranks=2,
        overrides=("output3/dcycle=1",),
    )

    source = _state_dump_at_cycle(source_dir, cycle=0)
    resumed = _state_dump_at_cycle(resumed_dir, cycle=0)
    assert source["var_names"] == resumed["var_names"]
    np.testing.assert_array_equal(
        np.asarray(source["mb_logical"])[_canonical_order(source)],
        np.asarray(resumed["mb_logical"])[_canonical_order(resumed)],
    )
    for variable in source["var_names"]:
        np.testing.assert_array_equal(
            _canonical_field(resumed, variable),
            _canonical_field(source, variable),
            err_msg=f"zero-step restart changed physical field {variable}",
        )


def test_device_fofc_pending_is_drained_into_restart(tmp_path: Path) -> None:
    source_dir = tmp_path / "device_fofc_source"
    resumed_dir = tmp_path / "device_fofc_resumed"
    _run(
        source_dir,
        ranks=1,
        overrides=(
            "time/nlim=1",
            "mhd/c2p_iter=1",
            "output1/dcycle=10",
        ),
    )
    checkpoint = source_dir / "rst" / "dyngrmhd_eventlog.00001.rst"
    pending = _restart_metadata(checkpoint, ranks=1)["event_counters"]
    assert isinstance(pending, dict)
    assert pending["nfofc"] == 16

    _run_restart(checkpoint, resumed_dir, nlim=1)
    resumed_rows = _event_rows(resumed_dir)
    assert [row["cycle"] for row in resumed_rows] == [1]
    assert resumed_rows[0]["fofc"] == pending["nfofc"]


def test_pending_event_counters_survive_shared_and_per_rank_restarts(
        tmp_path: Path) -> None:
    mesh_overrides = ("mesh/nx1=16", "meshblock/nx1=8")

    shared_dir = tmp_path / "shared_source"
    shared_resumed_dir = tmp_path / "shared_resumed"
    _run(
        shared_dir,
        ranks=2,
        overrides=(
            *mesh_overrides,
            "time/nlim=3",
            "output1/dcycle=10",
        ),
    )
    shared_checkpoint = (
        shared_dir / "rst" / "dyngrmhd_eventlog.00001.rst"
    )
    metadata = _restart_metadata(shared_checkpoint, ranks=2)
    assert metadata["event_counter_version"] == 2
    pending = metadata["event_counters"]
    assert isinstance(pending, dict)
    assert pending["nfofc_tests"] == 32

    # A zero-step two-rank restart isolates restoration from further FOFC evolution.
    # The shared global value is placed on rank zero only; otherwise the event-log
    # Allreduce would incorrectly turn 32 into 64.
    _run_restart(
        shared_checkpoint,
        shared_resumed_dir,
        nlim=1,
        ranks=2,
    )
    shared_rows = _event_rows(shared_resumed_dir)
    assert [row["cycle"] for row in shared_rows] == [1]
    assert shared_rows[0]["fofc_tests"] == pending["nfofc_tests"]

    # Remove exactly the fixed-size event extension to emulate a legacy shared
    # checkpoint.  New code must rewind its marker probe and read the old object/data
    # payload at the original boundary.
    event_extension_bytes = 8 + 2 * 4 + 10 * 8 + 4
    metadata_end = int(metadata["metadata_end"])
    shared_bytes = shared_checkpoint.read_bytes()
    legacy_checkpoint = tmp_path / "legacy_without_event_extension.rst"
    legacy_checkpoint.write_bytes(
        shared_bytes[:metadata_end - event_extension_bytes]
        + shared_bytes[metadata_end:]
    )
    legacy_resumed_dir = tmp_path / "legacy_resumed"
    _run_restart(
        legacy_checkpoint,
        legacy_resumed_dir,
        nlim=1,
        ranks=2,
    )
    legacy_rows = _event_rows(legacy_resumed_dir)
    assert [row["cycle"] for row in legacy_rows] == [1]
    assert legacy_rows[0]["fofc_tests"] == 0

    # Version 1 used the same fields but included ghost-zone events.  It cannot be
    # converted to the active-zone v2 denominator.  Refuse it by default; an explicit
    # one-time migration discards pending totals and emits a clean v2 checkpoint.
    event_marker = (0x41544B4556543031).to_bytes(8, byteorder=sys.byteorder)
    event_start = shared_bytes.find(event_marker)
    assert event_start >= 0
    assert shared_bytes.find(event_marker, event_start + 1) == -1
    legacy_v1_checkpoint = tmp_path / "legacy_v1_events.rst"
    legacy_v1_bytes = bytearray(shared_bytes)
    legacy_v1_bytes[event_start + 8:event_start + 12] = (1).to_bytes(
        4, byteorder=sys.byteorder, signed=True
    )
    legacy_v1_checkpoint.write_bytes(legacy_v1_bytes)
    _run_restart_expect_failure(
        legacy_v1_checkpoint,
        tmp_path / "legacy_v1_rejected",
        ranks=2,
        expected="allow_legacy_ghost_event_counters=true",
    )
    migrated_v1_dir = tmp_path / "legacy_v1_migrated"
    migration_stdout = _run_restart(
        legacy_v1_checkpoint,
        migrated_v1_dir,
        nlim=1,
        ranks=2,
        overrides=("time/allow_legacy_ghost_event_counters=true",),
    )
    assert "discarding pending legacy v1 event counters" in migration_stdout
    migrated_rows = _event_rows(migrated_v1_dir)
    assert [row["cycle"] for row in migrated_rows] == [1]
    assert migrated_rows[0]["fofc_tests"] == 0
    migrated_checkpoint = (
        migrated_v1_dir / "rst" / "dyngrmhd_eventlog.00002.rst"
    )
    assert migrated_checkpoint.is_file()
    migrated_metadata = _restart_metadata(migrated_checkpoint, ranks=2)
    assert migrated_metadata["event_counter_version"] == 2
    assert migrated_metadata["event_counters"]["nfofc_tests"] == 0
    assert _restart_parameter(
        migrated_checkpoint, "time", "allow_legacy_ghost_event_counters"
    ) in ("false", "0")

    per_rank_dir = tmp_path / "per_rank_source"
    per_rank_resumed_dir = tmp_path / "per_rank_resumed"
    _run(
        per_rank_dir,
        ranks=2,
        overrides=(
            *mesh_overrides,
            "time/nlim=3",
            "output1/dcycle=10",
            "output2/single_file_per_rank=true",
        ),
    )
    rank_zero_checkpoint = (
        per_rank_dir
        / "rst"
        / "rank_00000000"
        / "dyngrmhd_eventlog.00001.rst"
    )
    assert rank_zero_checkpoint.is_file()
    assert (
        per_rank_dir
        / "rst"
        / "rank_00000001"
        / rank_zero_checkpoint.name
    ).is_file()
    _run_restart(
        rank_zero_checkpoint,
        per_rank_resumed_dir,
        nlim=1,
        ranks=2,
    )
    per_rank_rows = _event_rows(per_rank_resumed_dir)
    assert [row["cycle"] for row in per_rank_rows] == [1]
    assert per_rank_rows[0]["fofc_tests"] == pending["nfofc_tests"]

    legacy_per_rank_root = tmp_path / "legacy_per_rank" / "rst"
    for rank in range(2):
        source = (
            per_rank_dir
            / "rst"
            / f"rank_{rank:08d}"
            / rank_zero_checkpoint.name
        )
        source_bytes = source.read_bytes()
        event_start = source_bytes.find(event_marker)
        assert event_start >= 0
        assert source_bytes.find(event_marker, event_start + 1) == -1
        destination_dir = legacy_per_rank_root / f"rank_{rank:08d}"
        destination_dir.mkdir(parents=True)
        (destination_dir / source.name).write_bytes(
            source_bytes[:event_start]
            + source_bytes[event_start + event_extension_bytes:]
        )
    legacy_per_rank_resumed_dir = tmp_path / "legacy_per_rank_resumed"
    _run_restart(
        legacy_per_rank_root / "rank_00000000" / rank_zero_checkpoint.name,
        legacy_per_rank_resumed_dir,
        nlim=1,
        ranks=2,
    )
    legacy_per_rank_rows = _event_rows(legacy_per_rank_resumed_dir)
    assert [row["cycle"] for row in legacy_per_rank_rows] == [1]
    assert legacy_per_rank_rows[0]["fofc_tests"] == 0

    # A mixed-generation per-rank set must fail before field I/O.  Otherwise one rank
    # would restore a carried accumulator while another silently treated it as absent.
    mixed_root = tmp_path / "mixed_event_extension" / "rst"
    for rank in range(2):
        destination_dir = mixed_root / f"rank_{rank:08d}"
        destination_dir.mkdir(parents=True)
        source_root = per_rank_dir / "rst" if rank == 0 else legacy_per_rank_root
        shutil.copy2(
            source_root / f"rank_{rank:08d}" / rank_zero_checkpoint.name,
            destination_dir / rank_zero_checkpoint.name,
        )
    _run_restart_expect_failure(
        mixed_root / "rank_00000000" / rank_zero_checkpoint.name,
        tmp_path / "mixed_event_extension_run",
        ranks=2,
        expected="Per-rank restart files disagree on event-counter metadata",
    )

    # Presence alone is insufficient: a mixed v1/v2 per-rank set would otherwise
    # restore counters on one rank and discard them on another.
    mixed_version_root = tmp_path / "mixed_event_version" / "rst"
    for rank in range(2):
        source = (
            per_rank_dir / "rst" / f"rank_{rank:08d}" / rank_zero_checkpoint.name
        )
        source_bytes = bytearray(source.read_bytes())
        version_start = source_bytes.find(event_marker) + 8
        assert version_start >= 8
        if rank == 1:
            source_bytes[version_start:version_start + 4] = (1).to_bytes(
                4, byteorder=sys.byteorder, signed=True
            )
        destination_dir = mixed_version_root / f"rank_{rank:08d}"
        destination_dir.mkdir(parents=True)
        (destination_dir / source.name).write_bytes(source_bytes)
    _run_restart_expect_failure(
        mixed_version_root / "rank_00000000" / rank_zero_checkpoint.name,
        tmp_path / "mixed_event_version_run",
        ranks=2,
        expected="Per-rank restart files disagree on event-counter format version",
    )


def test_invalid_telemetry_and_magnetization_limits_are_rejected(tmp_path: Path) -> None:
    for index, value in enumerate(("0", "-1", "nan", "inf")):
        stdout = _run_expect_failure(
            tmp_path / f"max_bsq_{index}", f"mhd/max_bsq={value}"
        )
        assert "mhd/max_bsq must be finite and > 0" in stdout

    stdout = _run_expect_failure(tmp_path / "negative_c2perrs", "mhd/c2perrs=-1")
    assert "mhd/c2perrs must be nonnegative" in stdout

    stdout = _run_expect_failure(tmp_path / "zero_c2p_iter", "mhd/c2p_iter=0")
    assert "mhd/c2p_iter must be at least 1" in stdout
