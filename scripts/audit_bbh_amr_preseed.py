#!/usr/bin/env python3
"""Reproduce and audit a moving-BBH cold-start AMR hierarchy.

The tool is intentionally read-only.  It reproduces AthenaK's ``<refined_region>``
tree construction, including the 26-neighbor strict-2:1 propagation performed by
``MeshBlockTree::Refine``.  It then compares every modeled leaf with a real
``mesh_structure.dat`` and checks that selected trajectory lookahead windows already
have every parent required by ``dynbbh::RefineTracker``.

Only cubic three-dimensional root grids with a power-of-two number of MeshBlocks per
axis are accepted.  This fail-closed restriction covers the effective-BBH campaign and
keeps the logical-location model exactly aligned with AthenaK's octree.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Iterable


POSITION = ((0, 1, 2), (3, 4, 5))
VELOCITY = ((6, 7, 8), (9, 10, 11))
SPIN = ((12, 13, 14), (15, 16, 17))
MASS = (18, 19)


class AuditError(ValueError):
    """An input or artifact cannot support a fail-closed audit."""


@dataclass(frozen=True, order=True)
class Location:
    level: int
    lx1: int
    lx2: int
    lx3: int

    @property
    def xyz(self) -> tuple[int, int, int]:
        return (self.lx1, self.lx2, self.lx3)

    def ancestor(self, level: int) -> "Location":
        if not 0 <= level <= self.level:
            raise AuditError("ancestor level lies outside the location path")
        shift = self.level - level
        return Location(
            level,
            self.lx1 >> shift,
            self.lx2 >> shift,
            self.lx3 >> shift,
        )


@dataclass(frozen=True)
class MeshRecord:
    gid: int
    rank: int
    cost: float
    location: Location


@dataclass(frozen=True)
class InputBlock:
    name: str
    values: dict[str, str]


@dataclass(frozen=True)
class Target:
    position: tuple[float, float, float]
    velocity: tuple[float, float, float]
    spin: tuple[float, float, float]
    mass: float


@dataclass(frozen=True)
class Geometry:
    lower: tuple[float, float, float]
    upper: tuple[float, float, float]
    root_blocks: int
    root_level: int
    max_physical_level: int

    def block_width(self, physical_level: int, axis: int) -> float:
        return (
            (self.upper[axis] - self.lower[axis])
            / self.root_blocks
            / (1 << physical_level)
        )

    def block_distance_squared(
        self,
        point: tuple[float, float, float],
        physical_level: int,
        xyz: tuple[int, int, int],
    ) -> float:
        total = 0.0
        for axis, (coordinate, index) in enumerate(zip(point, xyz)):
            width = self.block_width(physical_level, axis)
            lower = self.lower[axis] + index * width
            upper = lower + width
            delta = max(lower - coordinate, 0.0, coordinate - upper)
            total += delta * delta
        return total


class Node:
    """One node in the audit-only replica of AthenaK's MeshBlockTree."""

    def __init__(self, location: Location) -> None:
        self.location = location
        self.children: list[Node] | None = None


class Octree:
    """Minimal replica of AddNode/Refine for a complete cubic root grid."""

    def __init__(self, root_level: int) -> None:
        if root_level < 0:
            raise AuditError("root logical level cannot be negative")
        self.root = Node(Location(0, 0, 0, 0))
        self.root_level = root_level
        self._create_complete_root(self.root)

    @staticmethod
    def _make_children(node: Node) -> None:
        if node.children is not None:
            return
        location = node.location
        node.children = []
        for number in range(8):
            dx = number & 1
            dy = (number >> 1) & 1
            dz = (number >> 2) & 1
            node.children.append(
                Node(
                    Location(
                        location.level + 1,
                        2 * location.lx1 + dx,
                        2 * location.lx2 + dy,
                        2 * location.lx3 + dz,
                    )
                )
            )

    def _create_complete_root(self, node: Node) -> None:
        if node.location.level == self.root_level:
            return
        self._make_children(node)
        assert node.children is not None
        for child in node.children:
            self._create_complete_root(child)

    def refine(self, node: Node) -> None:
        """Match MeshBlockTree::Refine, including all 26 same-level neighbors."""

        if node.children is not None:
            return
        self._make_children(node)
        location = node.location
        limit = 1 << location.level
        for dz in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    if dx == dy == dz == 0:
                        continue
                    xyz = (
                        location.lx1 + dx,
                        location.lx2 + dy,
                        location.lx3 + dz,
                    )
                    if all(0 <= value < limit for value in xyz):
                        self.add_node(Location(location.level, *xyz))

    def add_node(self, target: Location) -> None:
        if target.level < 0 or any(
            value < 0 or value >= (1 << target.level) for value in target.xyz
        ):
            raise AuditError(f"logical target lies outside the octree: {target}")
        node = self.root
        while node.location.level < target.level:
            if node.children is None:
                self.refine(node)
            assert node.children is not None
            shift = target.level - node.location.level - 1
            dx = (target.lx1 >> shift) & 1
            dy = (target.lx2 >> shift) & 1
            dz = (target.lx3 >> shift) & 1
            node = node.children[dx + 2 * dy + 4 * dz]

    def _walk(self, node: Node, leaves: bool, result: set[Location]) -> None:
        if node.children is None:
            if leaves:
                result.add(node.location)
            return
        if not leaves:
            result.add(node.location)
        for child in node.children:
            self._walk(child, leaves, result)

    def leaf_locations(self) -> set[Location]:
        result: set[Location] = set()
        self._walk(self.root, True, result)
        return result

    def ordered_leaf_locations(self) -> list[Location]:
        """Return the child-0-through-7 Z-order AthenaK uses to assign GIDs."""

        result: list[Location] = []

        def append_leaves(node: Node) -> None:
            if node.children is None:
                result.append(node.location)
                return
            for child in node.children:
                append_leaves(child)

        append_leaves(self.root)
        return result

    def internal_locations(self) -> set[Location]:
        result: set[Location] = set()
        self._walk(self.root, False, result)
        return result


MESH_HEADER = re.compile(
    r"^#MeshBlock\s+(\d+)\s+on rank=(\d+)\s+with cost=([^\s]+)\s*$"
)
LOCATION_HEADER = re.compile(
    r"^#\s+Logical level\s+(\d+),\s+location = "
    r"\((\d+)\s+(\d+)\s+(\d+)\)\s*$"
)


def parse_mesh_structure_lines(lines: Iterable[str]) -> list[MeshRecord]:
    """Parse AthenaK mesh records while ignoring the following vertex geometry."""

    records: list[MeshRecord] = []
    pending: tuple[int, int, float] | None = None
    for line_number, line in enumerate(lines, 1):
        header = MESH_HEADER.match(line.rstrip("\n"))
        if header:
            if pending is not None:
                raise AuditError(
                    f"mesh record before line {line_number} has no logical location"
                )
            gid, rank = int(header.group(1)), int(header.group(2))
            try:
                cost = float(header.group(3))
            except ValueError as exc:
                raise AuditError(f"line {line_number}: invalid MeshBlock cost") from exc
            if not math.isfinite(cost) or not cost > 0.0:
                raise AuditError(f"line {line_number}: MeshBlock cost must be positive")
            pending = (gid, rank, cost)
            continue
        location = LOCATION_HEADER.match(line.rstrip("\n"))
        if location:
            if pending is None:
                raise AuditError(
                    f"line {line_number}: logical location has no MeshBlock header"
                )
            level, lx1, lx2, lx3 = map(int, location.groups())
            records.append(
                MeshRecord(*pending, Location(level, lx1, lx2, lx3))
            )
            pending = None
    if pending is not None:
        raise AuditError("final MeshBlock record has no logical location")
    if not records:
        raise AuditError("mesh_structure contains no MeshBlock records")
    records.sort(key=lambda record: record.gid)
    gids = [record.gid for record in records]
    if gids != list(range(len(records))):
        raise AuditError("MeshBlock GIDs must be unique and contiguous from zero")
    locations = [record.location for record in records]
    if len(set(locations)) != len(locations):
        raise AuditError("mesh_structure contains duplicate logical leaves")
    return records


def parse_mesh_structure(path: Path) -> list[MeshRecord]:
    try:
        with path.open(encoding="ascii") as stream:
            return parse_mesh_structure_lines(stream)
    except OSError as exc:
        raise AuditError(f"cannot read mesh structure {path}: {exc}") from exc


def parse_athinput_lines(lines: Iterable[str]) -> list[InputBlock]:
    blocks: list[InputBlock] = []
    current: InputBlock | None = None
    for line_number, raw_line in enumerate(lines, 1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            name = line[1:-1].strip()
            if not name:
                raise AuditError(f"line {line_number}: empty input block name")
            current = InputBlock(name, {})
            blocks.append(current)
            continue
        if current is None or "=" not in line:
            raise AuditError(f"line {line_number}: malformed parameter assignment")
        key, value = (part.strip() for part in line.split("=", 1))
        if not key or not value:
            raise AuditError(f"line {line_number}: empty parameter key or value")
        if key in current.values:
            raise AuditError(
                f"line {line_number}: duplicate parameter {current.name}/{key}"
            )
        current.values[key] = value
    if not blocks:
        raise AuditError("input contains no parameter blocks")
    return blocks


def parse_athinput(path: Path) -> list[InputBlock]:
    try:
        with path.open(encoding="utf-8") as stream:
            return parse_athinput_lines(stream)
    except OSError as exc:
        raise AuditError(f"cannot read input {path}: {exc}") from exc


def unique_block(blocks: list[InputBlock], name: str) -> InputBlock:
    matches = [block for block in blocks if block.name == name]
    if len(matches) != 1:
        raise AuditError(f"expected exactly one <{name}> block, found {len(matches)}")
    return matches[0]


def integer(block: InputBlock, name: str) -> int:
    try:
        return int(block.values[name])
    except KeyError as exc:
        raise AuditError(f"missing <{block.name}>/{name}") from exc
    except ValueError as exc:
        raise AuditError(f"invalid integer <{block.name}>/{name}") from exc


def real(block: InputBlock, name: str) -> float:
    try:
        value = float(block.values[name])
    except KeyError as exc:
        raise AuditError(f"missing <{block.name}>/{name}") from exc
    except ValueError as exc:
        raise AuditError(f"invalid real <{block.name}>/{name}") from exc
    if not math.isfinite(value):
        raise AuditError(f"non-finite <{block.name}>/{name}")
    return value


def optional_integer(block: InputBlock, name: str, default: int) -> int:
    """Read an optional integer through the audited fail-closed parser."""

    if name not in block.values:
        return default
    return integer(block, name)


def optional_real(block: InputBlock, name: str, default: float) -> float:
    """Read an optional finite real through the audited fail-closed parser."""

    if name not in block.values:
        return default
    return real(block, name)


def geometry_from_input(blocks: list[InputBlock]) -> Geometry:
    mesh = unique_block(blocks, "mesh")
    meshblock = unique_block(blocks, "meshblock")
    refinement = unique_block(blocks, "mesh_refinement")
    root_counts = []
    lower = []
    upper = []
    for axis in range(1, 4):
        cells = integer(mesh, f"nx{axis}")
        block_cells = integer(meshblock, f"nx{axis}")
        if cells <= 0 or block_cells <= 0 or cells % block_cells:
            raise AuditError(
                "mesh cells must be positive and divisible by MeshBlock cells"
            )
        root_counts.append(cells // block_cells)
        lower.append(real(mesh, f"x{axis}min"))
        upper.append(real(mesh, f"x{axis}max"))
    if len(set(root_counts)) != 1:
        raise AuditError("audit requires an equal root-MeshBlock count on all axes")
    root_blocks = root_counts[0]
    if root_blocks < 1 or root_blocks & (root_blocks - 1):
        raise AuditError("audit requires a power-of-two root-MeshBlock count")
    if any(not lo < hi for lo, hi in zip(lower, upper)):
        raise AuditError("mesh coordinate bounds must be strictly increasing")
    levels = integer(refinement, "num_levels")
    if levels < 2:
        raise AuditError("moving-BBH AMR audit requires at least two physical levels")
    return Geometry(
        tuple(lower),
        tuple(upper),
        root_blocks,
        root_blocks.bit_length() - 1,
        levels - 1,
    )


def _left_edge(index: int, count: int, lower: float, upper: float) -> float:
    return lower + (upper - lower) * index / count


def _region_index_bounds(
    target_level: int,
    region_lower: float,
    region_upper: float,
    domain_lower: float,
    domain_upper: float,
    root_blocks: int,
) -> tuple[int, int]:
    count = root_blocks << target_level
    minimum = 0
    while minimum < count:
        if _left_edge(minimum + 1, count, domain_lower, domain_upper) > region_lower:
            break
        minimum += 1
    maximum = minimum
    while maximum < count:
        if _left_edge(maximum + 1, count, domain_lower, domain_upper) >= region_upper:
            break
        maximum += 1
    if minimum & 1:
        minimum -= 1
    if maximum % 2 == 0:
        maximum += 1
    return minimum, maximum


def model_refined_regions(blocks: list[InputBlock], geometry: Geometry) -> Octree:
    tree = Octree(geometry.root_level)
    for block in blocks:
        if not block.name.startswith("refined_region"):
            continue
        target_level = integer(block, "level")
        if not 1 <= target_level <= geometry.max_physical_level:
            raise AuditError(f"<{block.name}> level lies outside the AMR hierarchy")
        bounds = []
        for axis in range(3):
            region_lower = real(block, f"x{axis + 1}min")
            region_upper = real(block, f"x{axis + 1}max")
            if region_lower > region_upper:
                raise AuditError(f"<{block.name}> has a reversed coordinate interval")
            if (
                region_lower < geometry.lower[axis]
                or region_upper > geometry.upper[axis]
            ):
                raise AuditError(f"<{block.name}> lies outside the root mesh")
            bounds.append(
                _region_index_bounds(
                    target_level,
                    region_lower,
                    region_upper,
                    geometry.lower[axis],
                    geometry.upper[axis],
                    geometry.root_blocks,
                )
            )
        logical_level = geometry.root_level + target_level
        for k in range(bounds[2][0], bounds[2][1], 2):
            for j in range(bounds[1][0], bounds[1][1], 2):
                for i in range(bounds[0][0], bounds[0][1], 2):
                    tree.add_node(Location(logical_level, i, j, k))
    return tree


def validate_leaf_partition(
    records: list[MeshRecord], geometry: Geometry
) -> list[str]:
    failures: list[str] = []
    leaves = {record.location for record in records}
    volume = Fraction(0, 1)
    for location in leaves:
        if location.level < geometry.root_level:
            failures.append(f"leaf below root logical level: {location}")
            continue
        limit = 1 << location.level
        if any(value < 0 or value >= limit for value in location.xyz):
            failures.append(f"leaf outside logical domain: {location}")
        volume += Fraction(1, 8**location.level)
        for level in range(geometry.root_level, location.level):
            if location.ancestor(level) in leaves:
                failures.append(f"overlapping ancestor and descendant at {location}")
                break
    if volume != 1:
        failures.append(f"logical leaves cover volume {volume}, expected 1")
    return failures


def read_trajectory(path: Path) -> list[tuple[float, tuple[float, ...]]]:
    rows = []
    try:
        with path.open(encoding="ascii") as stream:
            for line_number, raw_line in enumerate(stream, 1):
                payload = raw_line.split("#", 1)[0].strip()
                if not payload:
                    continue
                try:
                    values = tuple(float(token) for token in payload.split())
                except ValueError as exc:
                    raise AuditError(
                        f"trajectory line {line_number} contains a non-numeric field"
                    ) from exc
                if len(values) != 21 or not all(math.isfinite(value) for value in values):
                    raise AuditError(
                        f"trajectory line {line_number} must contain 21 finite fields"
                    )
                if rows and not values[0] > rows[-1][0]:
                    raise AuditError("trajectory times must be strictly increasing")
                rows.append((values[0], values[1:]))
    except OSError as exc:
        raise AuditError(f"cannot read trajectory {path}: {exc}") from exc
    if len(rows) < 2:
        raise AuditError("trajectory must contain at least two data rows")
    return rows


def interpolate_state(
    rows: list[tuple[float, tuple[float, ...]]], table_time: float
) -> tuple[float, ...]:
    if not rows[0][0] <= table_time <= rows[-1][0]:
        raise AuditError(f"sample time {table_time:.17g} lies outside the trajectory")
    upper_index = 0
    while upper_index < len(rows) and rows[upper_index][0] < table_time:
        upper_index += 1
    if upper_index == 0:
        return rows[0][1]
    if upper_index == len(rows):
        return rows[-1][1]
    lower = rows[upper_index - 1]
    upper = rows[upper_index]
    t0, state0 = lower
    t1, state1 = upper
    if table_time == t1:
        return state1
    interval = t1 - t0
    weight = (table_time - t0) / interval
    state = [a + weight * (b - a) for a, b in zip(state0, state1)]
    weight2 = weight * weight
    weight3 = weight2 * weight
    h00 = 2.0 * weight3 - 3.0 * weight2 + 1.0
    h10 = weight3 - 2.0 * weight2 + weight
    h01 = -2.0 * weight3 + 3.0 * weight2
    h11 = weight3 - weight2
    for hole in range(2):
        for position, velocity in zip(POSITION[hole], VELOCITY[hole]):
            p0, p1 = state0[position], state1[position]
            v0, v1 = state0[velocity], state1[velocity]
            state[position] = (
                h00 * p0 + h10 * interval * v0 + h01 * p1 + h11 * interval * v1
            )
            state[velocity] = (
                ((6.0 * weight2 - 6.0 * weight) / interval) * p0
                + (3.0 * weight2 - 4.0 * weight + 1.0) * v0
                + ((-6.0 * weight2 + 6.0 * weight) / interval) * p1
                + (3.0 * weight2 - 2.0 * weight) * v1
            )
    return tuple(state)


def _coincident(first: float, second: float) -> bool:
    scale = max(1.0, abs(first), abs(second))
    return abs(first - second) <= 128.0 * sys.float_info.epsilon * scale


def tracking_targets(state: tuple[float, ...], problem: InputBlock) -> tuple[Target, ...]:
    scales = (
        optional_real(problem, "mass_scale1", 1.0),
        optional_real(problem, "mass_scale2", 1.0),
    )
    targets = []
    for hole in range(2):
        scale = scales[hole]
        target = Target(
            tuple(state[index] for index in POSITION[hole]),
            tuple(state[index] for index in VELOCITY[hole]),
            tuple(state[index] * scale for index in SPIN[hole]),
            state[MASS[hole]] * scale,
        )
        if target.mass > 0.0:
            targets.append(target)
    if len(targets) != 2:
        return tuple(targets)
    first, second = targets
    coincident = all(
        _coincident(a, b)
        for first_values, second_values in (
            (first.position, second.position),
            (first.velocity, second.velocity),
            (first.spin, second.spin),
        )
        for a, b in zip(first_values, second_values)
    )
    if not coincident:
        return tuple(targets)
    return (
        Target(
            tuple((a + b) / 2.0 for a, b in zip(first.position, second.position)),
            tuple((a + b) / 2.0 for a, b in zip(first.velocity, second.velocity)),
            tuple((a + b) / 2.0 for a, b in zip(first.spin, second.spin)),
            first.mass + second.mass,
        ),
    )


def horizon_guard_radius(target: Target, factor: float) -> float:
    spin2 = sum(value * value for value in target.spin)
    horizon = target.mass + math.sqrt(max(target.mass * target.mass - spin2, 0.0))
    rest_enclosing = math.sqrt(spin2 + horizon * horizon)
    velocity2 = sum(value * value for value in target.velocity)
    if not velocity2 < 1.0:
        raise AuditError("trajectory sample contains a non-subluminal target")
    return factor * rest_enclosing / math.sqrt(1.0 - velocity2)


def intersecting_parents(
    geometry: Geometry,
    physical_level: int,
    point: tuple[float, float, float],
    radius: float,
) -> set[tuple[int, int, int]]:
    if not radius > 0.0 or not math.isfinite(radius):
        raise AuditError("tracker radius must be finite and positive")
    limit = geometry.root_blocks << physical_level
    ranges = []
    for axis, coordinate in enumerate(point):
        width = geometry.block_width(physical_level, axis)
        lower = max(
            0,
            math.floor((coordinate - radius - geometry.lower[axis]) / width) - 1,
        )
        upper = min(
            limit - 1,
            math.floor((coordinate + radius - geometry.lower[axis]) / width) + 1,
        )
        ranges.append(range(lower, upper + 1))
    result = set()
    for k in ranges[2]:
        for j in ranges[1]:
            for i in ranges[0]:
                xyz = (i, j, k)
                if geometry.block_distance_squared(
                    point, physical_level, xyz
                ) < radius * radius:
                    result.add(xyz)
    return result


def shell_radii(problem: InputBlock, maximum_level: int) -> dict[int, float]:
    explicit = any(
        f"refinement_radius_level_{level}" in problem.values
        for level in range(1, maximum_level + 1)
    )
    base = optional_real(problem, "refinement_radius", 6.0)
    ratio = optional_real(problem, "refinement_radius_ratio", 2.0)
    result = {}
    for level in range(1, maximum_level + 1):
        name = f"refinement_radius_level_{level}"
        if explicit and name not in problem.values:
            raise AuditError("explicit tracker shell table is incomplete")
        value = real(problem, name) if explicit else base * ratio ** (
            maximum_level - level
        )
        if not value > 0.0 or not math.isfinite(value):
            raise AuditError("tracker shell radii must be finite and positive")
        if level > 1 and value > result[level - 1]:
            raise AuditError("tracker shell radii must be non-increasing")
        result[level] = value
    return result


def required_parents_for_window(
    blocks: list[InputBlock],
    geometry: Geometry,
    trajectory: list[tuple[float, tuple[float, ...]]],
    start_time: float,
    horizon: float,
) -> dict[int, set[tuple[int, int, int]]]:
    problem = unique_block(blocks, "problem")
    offset = real(problem, "trajectory_time_offset")
    factor = optional_real(problem, "refinement_horizon_factor", 1.25)
    radii = shell_radii(problem, geometry.max_physical_level)
    padding = 0.25 * horizon
    required = {
        level: set() for level in range(1, geometry.max_physical_level + 1)
    }
    for simulation_time in (start_time, start_time + 0.5 * horizon, start_time + horizon):
        state = interpolate_state(trajectory, simulation_time + offset)
        for target in tracking_targets(state, problem):
            for child_level in range(1, geometry.max_physical_level + 1):
                radius = max(
                    radii[child_level], horizon_guard_radius(target, factor)
                ) + padding
                required[child_level].update(
                    intersecting_parents(
                        geometry,
                        child_level - 1,
                        target.position,
                        radius,
                    )
                )
    return required


def required_floor_parents(
    blocks: list[InputBlock], geometry: Geometry
) -> dict[int, set[tuple[int, int, int]]]:
    problem = unique_block(blocks, "problem")
    radius = optional_real(problem, "refinement_floor_radius", 0.0)
    floor_level = optional_integer(problem, "refinement_floor_level", 0)
    center = tuple(
        optional_real(problem, f"refinement_floor_center{axis}", 0.0)
        for axis in range(1, 4)
    )
    result = {}
    for child_level in range(1, floor_level + 1):
        result[child_level] = intersecting_parents(
            geometry, child_level - 1, center, radius
        )
    return result


def replay_load_balance(
    costs: list[float], ranks: int, maximum_blocks: int
) -> list[tuple[int, int, float]]:
    """Reproduce the capacity-constrained contiguous partition used by Mesh."""

    def ordered_binary64_sum(values: Iterable[float]) -> float:
        """Match C++ ``double total += static_cast<double>(float_cost)``."""

        total = 0.0
        for value in values:
            total += value
        return total

    blocks = len(costs)
    if ranks < 1 or blocks < ranks:
        raise AuditError("cannot assign at least one MeshBlock to every rank")
    if maximum_blocks < 1 or blocks > ranks * maximum_blocks:
        raise AuditError("MeshBlocks exceed aggregate rank capacity")
    starts = [0] * ranks
    counts = [0] * ranks
    segment_end = blocks
    remaining_cost = ordered_binary64_sum(costs)
    for rank in range(ranks - 1, 0, -1):
        minimum_begin = max(rank, segment_end - maximum_blocks)
        maximum_begin = min(segment_end - 1, rank * maximum_blocks)
        if minimum_begin > maximum_begin:
            raise AuditError(f"no feasible capacity boundary for rank {rank}")
        segment_begin = maximum_begin
        segment_cost = ordered_binary64_sum(costs[segment_begin:segment_end])
        target = remaining_cost / (rank + 1)
        while segment_begin > minimum_begin:
            candidate = segment_cost + costs[segment_begin - 1]
            if abs(candidate - target) > abs(segment_cost - target):
                break
            segment_begin -= 1
            segment_cost = candidate
        starts[rank] = segment_begin
        counts[rank] = segment_end - segment_begin
        segment_end = segment_begin
        remaining_cost = max(0.0, remaining_cost - segment_cost)
    counts[0] = segment_end
    return [
        (
            starts[rank],
            counts[rank],
            ordered_binary64_sum(
                costs[starts[rank]:starts[rank] + counts[rank]]
            ),
        )
        for rank in range(ranks)
    ]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit(args: argparse.Namespace) -> tuple[dict[str, object], list[str]]:
    blocks = parse_athinput(args.input)
    geometry = geometry_from_input(blocks)
    records = parse_mesh_structure(args.mesh_structure)
    trajectory = read_trajectory(args.trajectory)
    failures = validate_leaf_partition(records, geometry)

    observed = {record.location for record in records}
    modeled_tree = model_refined_regions(blocks, geometry)
    modeled = modeled_tree.leaf_locations()
    observed_only = observed - modeled
    modeled_only = modeled - observed
    if observed_only or modeled_only:
        failures.append(
            "modeled <refined_region> hierarchy differs from mesh_structure: "
            f"observed_only={len(observed_only)}, modeled_only={len(modeled_only)}"
        )

    level_counts = Counter(
        record.location.level - geometry.root_level for record in records
    )
    expected_costs = [
        float(1 << (record.location.level - geometry.root_level))
        for record in records
    ]
    costs_exact = all(
        record.cost == expected for record, expected in zip(records, expected_costs)
    )
    if not costs_exact:
        failures.append("mesh costs do not match strict level-2:1 subcycling")
    weighted_cost = int(sum(expected_costs))
    if args.expected_blocks is not None and len(records) != args.expected_blocks:
        failures.append(
            f"block count {len(records)} differs from expected {args.expected_blocks}"
        )
    if args.expected_weighted is not None and weighted_cost != args.expected_weighted:
        failures.append(
            "weighted cost "
            f"{weighted_cost} differs from expected {args.expected_weighted}"
        )

    internal = modeled_tree.internal_locations()
    window_rows = []
    for start in args.window_start:
        required = required_parents_for_window(
            blocks, geometry, trajectory, start, args.window_horizon
        )
        missing_by_level = {}
        requested_by_level = {}
        for child_level, parents in required.items():
            requested_by_level[str(child_level)] = len(parents)
            missing = [
                xyz
                for xyz in parents
                if Location(
                    geometry.root_level + child_level - 1, *xyz
                ) not in internal
            ]
            missing_by_level[str(child_level)] = len(missing)
            if missing:
                failures.append(
                    f"window {start:g}M lacks {len(missing)} L{child_level} parents"
                )
        window_rows.append(
            {
                "start_M": start,
                "sample_times_M": [
                    start,
                    start + 0.5 * args.window_horizon,
                    start + args.window_horizon,
                ],
                "requested_parents": requested_by_level,
                "missing_parents": missing_by_level,
                "direct_refine_flags": sum(missing_by_level.values()),
            }
        )

    floor_required = required_floor_parents(blocks, geometry)
    floor_missing = {}
    for child_level, parents in floor_required.items():
        missing = sum(
            Location(geometry.root_level + child_level - 1, *xyz) not in internal
            for xyz in parents
        )
        floor_missing[str(child_level)] = missing
        if missing:
            failures.append(f"persistent floor lacks {missing} L{child_level} parents")

    refinement = unique_block(blocks, "mesh_refinement")
    capacity = integer(refinement, "max_nmb_per_rank")
    partitions = replay_load_balance(expected_costs, args.ranks, capacity)
    partition_rows = [
        {
            "rank": rank,
            "gid_start": start,
            "blocks": count,
            "cost": cost,
        }
        for rank, (start, count, cost) in enumerate(partitions)
    ]

    report: dict[str, object] = {
        "status": "pass" if not failures else "fail",
        "input": str(args.input),
        "input_sha256": sha256(args.input),
        "trajectory": str(args.trajectory),
        "trajectory_sha256": sha256(args.trajectory),
        "mesh_structure": str(args.mesh_structure),
        "mesh_structure_sha256": sha256(args.mesh_structure),
        "meshblocks": len(records),
        "physical_level_counts": {
            str(level): level_counts[level]
            for level in range(geometry.max_physical_level + 1)
        },
        "weighted_meshblock_updates_per_root_step": weighted_cost,
        "strict_subcycling_costs_exact": costs_exact,
        "refined_region_model_exact": not observed_only and not modeled_only,
        "observed_only_leaves": len(observed_only),
        "modeled_only_leaves": len(modeled_only),
        "windows": window_rows,
        "floor_missing_parents": floor_missing,
        "partition": {
            "ranks": args.ranks,
            "max_nmb_per_rank": capacity,
            "minimum_blocks": min(row[1] for row in partitions),
            "maximum_blocks": max(row[1] for row in partitions),
            "minimum_cost": min(row[2] for row in partitions),
            "maximum_cost": max(row[2] for row in partitions),
            "minimum_capacity_headroom": min(capacity - row[1] for row in partitions),
            "rows": partition_rows,
        },
        "failures": failures,
    }
    return report, failures


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--trajectory", required=True, type=Path)
    parser.add_argument("--mesh-structure", required=True, type=Path)
    parser.add_argument(
        "--window-start",
        action="append",
        type=float,
        default=[],
        help="root-lookahead window start in M; repeat for multiple windows",
    )
    parser.add_argument(
        "--window-horizon",
        type=float,
        default=None,
        help="lookahead horizon in M (default: <time>/root_dt_max)",
    )
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--expected-blocks", type=int)
    parser.add_argument("--expected-weighted", type=int)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    blocks = parse_athinput(args.input)
    root_dt = real(unique_block(blocks, "time"), "root_dt_max")
    if args.window_horizon is None:
        args.window_horizon = root_dt
    if not args.window_horizon > 0.0 or not math.isfinite(args.window_horizon):
        parser.error("--window-horizon must be finite and positive")
    if not args.window_start:
        args.window_start = [0.0, args.window_horizon]
    if args.ranks < 1:
        parser.error("--ranks must be positive")
    return args


def main() -> int:
    try:
        args = parse_args()
        report, failures = audit(args)
    except (AuditError, OSError) as exc:
        print(f"AMR preseed audit error: {exc}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        print(f"status={report['status']}")
        print(f"meshblocks={report['meshblocks']}")
        print(
            "physical_levels="
            + ",".join(
                str(value) for value in report["physical_level_counts"].values()
            )
        )
        print(
            "weighted_updates="
            f"{report['weighted_meshblock_updates_per_root_step']}"
        )
        print(f"refined_region_model_exact={report['refined_region_model_exact']}")
        for row in report["windows"]:
            print(
                f"window_start_M={row['start_M']:g} "
                f"direct_refine_flags={row['direct_refine_flags']}"
            )
        partition = report["partition"]
        print(
            f"partition_blocks={partition['minimum_blocks']}.."
            f"{partition['maximum_blocks']} "
            f"partition_cost={partition['minimum_cost']:g}.."
            f"{partition['maximum_cost']:g} "
            f"minimum_headroom={partition['minimum_capacity_headroom']}"
        )
        for failure in failures:
            print(f"failure={failure}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
