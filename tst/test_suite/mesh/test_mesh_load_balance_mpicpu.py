"""Regression for weighted contiguous MeshBlock partitioning."""

import re
from pathlib import Path

import test_suite.testutils as testutils


INPUT_FILE = "inputs/load_balance_weighted.athinput"


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
