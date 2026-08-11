"""Regression for weighted contiguous MeshBlock partitioning."""

import re
from pathlib import Path

import test_suite.testutils as testutils


INPUT_FILE = "inputs/load_balance_weighted.athinput"
CAPACITY_INPUT_FILE = "inputs/load_balance_capacity.athinput"


def test_weighted_partition_keeps_every_rank_nonempty():
    """A high-cost Z-order prefix must not leave lower MPI ranks empty."""
    log_path = Path(testutils.LOG_FILE_PATH)
    log_offset = log_path.stat().st_size if log_path.exists() else 0

    assert testutils.mpi_run(INPUT_FILE, flags=["-m"], threads=3)

    with log_path.open("rb") as log_file:
        log_file.seek(log_offset)
        log_output = log_file.read().decode("utf-8", errors="replace")

    assert "Total number of MeshBlocks = 4" in log_output
    rank_rows = re.findall(
        r"Rank = (\d+): (\d+) MeshBlocks, cost = (\d+)", log_output
    )
    assert rank_rows == [("0", "1", "4"), ("1", "1", "4"), ("2", "2", "3")]


def test_weighted_partition_honors_amr_block_capacity():
    """Cheap coarse blocks must not overflow a fixed-capacity rank partition."""
    log_path = Path(testutils.LOG_FILE_PATH)
    log_offset = log_path.stat().st_size if log_path.exists() else 0

    assert testutils.mpi_run(CAPACITY_INPUT_FILE, flags=["-m"], threads=2)

    with log_path.open("rb") as log_file:
        log_file.seek(log_offset)
        log_output = log_file.read().decode("utf-8", errors="replace")

    assert "Total number of MeshBlocks = 11" in log_output
    rank_rows = re.findall(
        r"Rank = (\d+): (\d+) MeshBlocks, cost = (\d+)", log_output
    )
    assert rank_rows == [("0", "6", "6"), ("1", "5", "23")]
