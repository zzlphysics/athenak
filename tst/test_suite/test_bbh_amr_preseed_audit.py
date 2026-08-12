"""Tests for the read-only moving-BBH AMR preseed auditor."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
AUDITOR_PATH = ROOT / "scripts" / "audit_bbh_amr_preseed.py"
SPEC = importlib.util.spec_from_file_location("bbh_amr_preseed_auditor", AUDITOR_PATH)
assert SPEC is not None and SPEC.loader is not None
AUDITOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = AUDITOR
SPEC.loader.exec_module(AUDITOR)


# This preserves AthenaK's real header, logical-location, blank-line, and vertex-row
# layout.  As in production mesh_structure.dat, records need not occur in GID order.
REAL_MESH_STRUCTURE_FRAGMENT = """\
#MeshBlock 1 on rank=0 with cost=2
#  Logical level 3, location = (1 0 0)
-512 -1024 -1024
-256 -1024 -1024
-256 -768 -1024
-512 -768 -1024
-512 -1024 -1024
-512 -1024 -768
-256 -1024 -768
-256 -1024 -1024
-256 -1024 -768
-256 -768 -768
-256 -768 -1024
-256 -768 -768
-512 -768 -768
-512 -768 -1024
-512 -768 -768
-512 -1024 -768
-512 -1024 -1024


#MeshBlock 0 on rank=0 with cost=1
#  Logical level 2, location = (3 0 0)
512 -1024 -1024
1024 -1024 -1024
1024 -512 -1024
512 -512 -1024
512 -1024 -1024
512 -1024 -512
1024 -1024 -512
1024 -1024 -1024
1024 -1024 -512
1024 -512 -512
1024 -512 -1024
1024 -512 -512
512 -512 -512
512 -512 -1024
512 -512 -512
512 -1024 -512
512 -1024 -1024
"""


def test_parser_accepts_real_mesh_structure_layout_and_sorts_gids() -> None:
    records = AUDITOR.parse_mesh_structure_lines(
        REAL_MESH_STRUCTURE_FRAGMENT.splitlines(keepends=True)
    )

    assert [(record.gid, record.rank, record.cost) for record in records] == [
        (0, 0, 1.0),
        (1, 0, 2.0),
    ]
    assert records[0].location == AUDITOR.Location(2, 3, 0, 0)
    assert records[1].location == AUDITOR.Location(3, 1, 0, 0)


def test_parser_rejects_a_meshblock_without_logical_location() -> None:
    with pytest.raises(AUDITOR.AuditError, match="has no logical location"):
        AUDITOR.parse_mesh_structure_lines(
            [
                "#MeshBlock 0 on rank=0 with cost=1\n",
                "0 0 0\n",
                "#MeshBlock 1 on rank=0 with cost=1\n",
                "#  Logical level 0, location = (0 0 0)\n",
            ]
        )


def test_refined_region_model_reproduces_strict_neighbor_propagation() -> None:
    blocks = AUDITOR.parse_athinput_lines(
        """
<mesh>
nx1 = 32
x1min = -1
x1max = 1
nx2 = 32
x2min = -1
x2max = 1
nx3 = 32
x3min = -1
x3max = 1
<meshblock>
nx1 = 16
nx2 = 16
nx3 = 16
<mesh_refinement>
num_levels = 2
max_nmb_per_rank = 16
<refined_region1>
level = 1
x1min = -0.5
x1max = -0.5
x2min = -0.5
x2max = -0.5
x3min = -0.5
x3max = -0.5
""".splitlines()
    )
    geometry = AUDITOR.geometry_from_input(blocks)
    tree = AUDITOR.model_refined_regions(blocks, geometry)
    counts = {}
    for location in tree.leaf_locations():
        physical_level = location.level - geometry.root_level
        counts[physical_level] = counts.get(physical_level, 0) + 1

    # One of the eight root MeshBlocks becomes an octet.  All seven root neighbors
    # already exist, so strict neighbor propagation adds no further refinement.
    assert counts == {0: 7, 1: 8}
    assert AUDITOR.Location(geometry.root_level, 0, 0, 0) in tree.internal_locations()


def test_capacity_partition_replay_is_fail_closed() -> None:
    assert AUDITOR.replay_load_balance(
        [1.0] * 6 + [4.0, 4.0, 5.0, 5.0, 5.0], 2, 6
    ) == [(0, 6, 6.0), (6, 5, 23.0)]
    with pytest.raises(AUDITOR.AuditError, match="aggregate rank capacity"):
        AUDITOR.replay_load_balance([1.0] * 9, 2, 4)


def test_load_balance_uses_cpp_ordered_binary64_accumulation() -> None:
    costs = [1.0, 2.0**53, 1.0]

    # C++ accumulates left-to-right into one double.  Both unit increments are lost
    # separately at 2**53, whereas Python 3.12's compensated sum returns 2**53 + 2.
    assert AUDITOR.replay_load_balance(costs, 2, 2) == [
        (0, 2, 2.0**53),
        (2, 1, 1.0),
    ]


def test_invalid_optional_problem_numbers_raise_audit_error() -> None:
    invalid_real = AUDITOR.InputBlock("problem", {"mass_scale1": "not-a-real"})
    with pytest.raises(AUDITOR.AuditError, match="invalid real <problem>/mass_scale1"):
        AUDITOR.tracking_targets((0.0,) * 20, invalid_real)

    invalid_integer = AUDITOR.InputBlock(
        "problem", {"refinement_floor_level": "not-an-integer"}
    )
    geometry = AUDITOR.Geometry(
        (-1.0, -1.0, -1.0),
        (1.0, 1.0, 1.0),
        root_blocks=1,
        root_level=0,
        max_physical_level=1,
    )
    with pytest.raises(
        AUDITOR.AuditError,
        match="invalid integer <problem>/refinement_floor_level",
    ):
        AUDITOR.required_floor_parents([invalid_integer], geometry)


def test_current_v100_input_has_the_audited_pure_model_topology() -> None:
    input_path = (
        ROOT / "inputs" / "dyngr" / "effective_bbh_4pn_v100_qualification.athinput"
    )
    blocks = AUDITOR.parse_athinput(input_path)
    geometry = AUDITOR.geometry_from_input(blocks)
    tree = AUDITOR.model_refined_regions(blocks, geometry)
    locations = tree.ordered_leaf_locations()
    level_counts = [
        sum(
            location.level - geometry.root_level == physical_level
            for location in locations
        )
        for physical_level in range(geometry.max_physical_level + 1)
    ]
    costs = [
        float(1 << (location.level - geometry.root_level))
        for location in locations
    ]
    partitions = AUDITOR.replay_load_balance(costs, ranks=8, maximum_blocks=1024)

    assert len(locations) == 4320
    assert level_counts == [4, 420, 420, 420, 420, 420, 396, 588, 592, 640]
    assert int(sum(costs)) == 605884
    assert max(row[1] for row in partitions) == 724
