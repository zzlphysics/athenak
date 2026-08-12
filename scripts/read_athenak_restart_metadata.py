#!/usr/bin/env python3
"""Inspect AthenaK single-file restart metadata without reading field arrays."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import json
import math
from pathlib import Path
import struct
from typing import BinaryIO


AMR_COUNTER_MAGIC = 0x41544B414D524331
AMR_COUNTER_VERSION = 1
EVENT_COUNTER_MAGIC = 0x41544B4556543031
EVENT_COUNTER_VERSION = 1
EVENT_SUM_COUNTER_COUNT = 10
MAX_PARAMETER_BYTES = 16 * 1024 * 1024
MAX_MESH_METADATA_BYTES = 512 * 1024 * 1024
REGION_INDEX_COUNT = 19


@dataclass(frozen=True)
class LogicalLocation:
    """Logical MeshBlock coordinates stored in a restart."""

    lx1: int
    lx2: int
    lx3: int
    level: int


@dataclass(frozen=True)
class RankPartition:
    """One contiguous GID range produced by AthenaK's LoadBalance."""

    rank: int
    gid_start: int
    blocks: int
    cost: float


@dataclass(frozen=True)
class RestartEventCounters:
    """Pending diagnostic accumulators persisted at a checkpoint boundary."""

    neos_dfloor: int
    neos_efloor: int
    neos_tfloor: int
    neos_vceil: int
    neos_fail: int
    nfofc: int
    ncons_adjust: int
    nmag_adjust: int
    nc2p_calls: int
    nfofc_tests: int
    maxit_c2p: int

    def as_dict(self) -> dict[str, int]:
        """Return fields in the stable restart serialization order."""

        return {
            "neos_dfloor": self.neos_dfloor,
            "neos_efloor": self.neos_efloor,
            "neos_tfloor": self.neos_tfloor,
            "neos_vceil": self.neos_vceil,
            "neos_fail": self.neos_fail,
            "nfofc": self.nfofc,
            "ncons_adjust": self.ncons_adjust,
            "nmag_adjust": self.nmag_adjust,
            "nc2p_calls": self.nc2p_calls,
            "nfofc_tests": self.nfofc_tests,
            "maxit_c2p": self.maxit_c2p,
        }


@dataclass(frozen=True)
class RestartMetadata:
    """Restart state read before the first field-data byte."""

    source: str
    file_size: int
    parameter_end: int
    metadata_end: int
    max_read_offset: int
    byte_order: str
    real_bytes: int
    nmb_total: int
    root_level: int
    time: float
    last_dt: float
    cycle: int
    locations: tuple[LogicalLocation, ...]
    costs: tuple[float, ...]
    amr_cycle_counters: tuple[int, ...] | None
    event_counters: RestartEventCounters | None
    parameters: dict[str, dict[str, str]]

    def physical_level_rows(self) -> list[dict[str, int | float]]:
        """Return MeshBlock counts and stored work cost for every populated level."""

        blocks: Counter[int] = Counter()
        costs: defaultdict[int, float] = defaultdict(float)
        for location, cost in zip(self.locations, self.costs):
            physical_level = location.level - self.root_level
            blocks[physical_level] += 1
            costs[physical_level] += cost
        return [
            {
                "physical_level": level,
                "logical_level": level + self.root_level,
                "blocks": blocks[level],
                "cost": costs[level],
            }
            for level in sorted(blocks)
        ]

    def counter_rows(self) -> list[dict[str, int]] | None:
        """Return AMR age-counter counts, split by physical level."""

        if self.amr_cycle_counters is None:
            return None
        counts: Counter[tuple[int, int]] = Counter()
        for location, counter in zip(self.locations, self.amr_cycle_counters):
            counts[(location.level - self.root_level, counter)] += 1
        return [
            {
                "physical_level": level,
                "cycles_since_refinement": counter,
                "blocks": counts[(level, counter)],
            }
            for level, counter in sorted(counts)
        ]


@dataclass(frozen=True)
class FixedHeader:
    """Validated native-ABI restart header preceding the MeshBlock lists."""

    nmb_total: int
    root_level: int
    root_blocks: tuple[int, int, int]
    time: float
    last_dt: float
    cycle: int


class AuditedReader:
    """Track the highest byte touched while retaining normal seekable-file behavior."""

    def __init__(self, stream: BinaryIO):
        self.stream = stream
        self.max_read_offset = stream.tell()

    def read(self, size: int) -> bytes:
        data = self.stream.read(size)
        self.max_read_offset = max(self.max_read_offset, self.stream.tell())
        return data

    def readline(self, size: int) -> bytes:
        data = self.stream.readline(size)
        self.max_read_offset = max(self.max_read_offset, self.stream.tell())
        return data

    def seek(self, offset: int) -> int:
        return self.stream.seek(offset)

    def tell(self) -> int:
        return self.stream.tell()


def _read_exact(reader: AuditedReader, size: int, description: str) -> bytes:
    pieces: list[bytes] = []
    remaining = size
    while remaining:
        piece = reader.read(remaining)
        if not piece:
            break
        pieces.append(piece)
        remaining -= len(piece)
    data = b"".join(pieces)
    if remaining:
        raise ValueError(
            f"truncated restart while reading {description}: expected {size} bytes, "
            f"got {len(data)}"
        )
    return data


def _read_parameter_dump(reader: AuditedReader) -> tuple[bytes, int]:
    pieces: list[bytes] = []
    total = 0
    while total < MAX_PARAMETER_BYTES:
        line = reader.readline(MAX_PARAMETER_BYTES - total + 1)
        if not line:
            raise ValueError("restart has no <par_end> parameter terminator")
        total += len(line)
        if total > MAX_PARAMETER_BYTES:
            raise ValueError(
                f"restart parameter dump exceeds {MAX_PARAMETER_BYTES} bytes"
            )
        pieces.append(line)
        if line.strip() == b"<par_end>":
            return b"".join(pieces), reader.tell()
    raise ValueError(f"restart parameter dump exceeds {MAX_PARAMETER_BYTES} bytes")


def _parse_parameters(raw: bytes) -> dict[str, dict[str, str]]:
    parameters: dict[str, dict[str, str]] = {}
    block: str | None = None
    for raw_line in raw.decode("utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            name = line[1:-1]
            if name == "par_end":
                break
            block = name
            parameters.setdefault(block, {})
            continue
        if block is not None and "=" in line:
            key, value = line.split("=", 1)
            parameters[block][key.strip()] = value.strip()
    return parameters


def _integer_parameter(
        parameters: dict[str, dict[str, str]], block: str, key: str) -> int:
    try:
        return int(parameters[block][key])
    except (KeyError, ValueError) as exc:
        raise ValueError(f"restart parameter <{block}>/{key} is missing or invalid") \
            from exc


def _float_parameter(
        parameters: dict[str, dict[str, str]], block: str, key: str) -> float:
    try:
        return float(parameters[block][key])
    except (KeyError, ValueError) as exc:
        raise ValueError(f"restart parameter <{block}>/{key} is missing or invalid") \
            from exc


def _expected_root_level(root_blocks: tuple[int, int, int]) -> int:
    maximum = max(root_blocks)
    level = 0
    while (1 << level) < maximum:
        level += 1
    return level


def _close(value: float, expected: float) -> bool:
    return math.isclose(value, expected, rel_tol=2.0e-6, abs_tol=1.0e-12)


def _mesh_metadata_size(
        nmb_total: int, include_counters: bool,
        include_event_counters: bool = False) -> int:
    """Return variable-size tree metadata bytes for the configured extension."""

    size = nmb_total * (4 * 4 + 4)
    if include_counters:
        size += struct.calcsize("Qii") + nmb_total * struct.calcsize("i")
    if include_event_counters:
        size += (
            struct.calcsize("Qii")
            + EVENT_SUM_COUNTER_COUNT * struct.calcsize("Q")
            + struct.calcsize("i")
        )
    return size


def _validate_metadata_capacity(
        nmb_total: int, include_counters: bool,
        include_event_counters: bool = False) -> None:
    metadata_bytes = _mesh_metadata_size(
        nmb_total, include_counters, include_event_counters
    )
    if metadata_bytes > MAX_MESH_METADATA_BYTES:
        extensions = []
        if include_counters:
            extensions.append("AMR counters")
        if include_event_counters:
            extensions.append("event counters")
        extension = f" including {' and '.join(extensions)}" if extensions else ""
        raise ValueError(
            f"MeshBlock metadata{extension} requires {metadata_bytes} bytes, "
            f"above safety limit"
        )


def _parse_fixed_header(
        reader: AuditedReader, file_size: int, parameter_end: int,
        parameters: dict[str, dict[str, str]], byte_order: str,
        real_bytes: int) -> FixedHeader:
    """Parse and validate only the fixed-size native-ABI header."""

    reader.seek(parameter_end)
    real_code = "f" if real_bytes == 4 else "d"
    int_struct = struct.Struct(f"{byte_order}ii")
    nmb_total, root_level = int_struct.unpack(
        _read_exact(reader, int_struct.size, "MeshBlock count and root level")
    )
    if nmb_total < 1:
        raise ValueError("non-positive MeshBlock count")
    _validate_metadata_capacity(nmb_total, include_counters=False)

    region_struct = struct.Struct(f"{byte_order}9{real_code}")
    mesh_size = region_struct.unpack(
        _read_exact(reader, region_struct.size, "mesh RegionSize")
    )
    indices_struct = struct.Struct(f"{byte_order}{REGION_INDEX_COUNT}i")
    mesh_indices = indices_struct.unpack(
        _read_exact(reader, indices_struct.size, "mesh RegionIndcs")
    )
    block_indices = indices_struct.unpack(
        _read_exact(reader, indices_struct.size, "MeshBlock RegionIndcs")
    )
    tail_struct = struct.Struct(f"{byte_order}{real_code}{real_code}i")
    time, last_dt, cycle = tail_struct.unpack(
        _read_exact(reader, tail_struct.size, "time, timestep, and cycle")
    )

    expected_mesh_cells = tuple(
        _integer_parameter(parameters, "mesh", f"nx{axis}")
        for axis in range(1, 4)
    )
    expected_block_cells = tuple(
        _integer_parameter(parameters, "meshblock", f"nx{axis}")
        for axis in range(1, 4)
    )
    expected_nghost = _integer_parameter(parameters, "mesh", "nghost")
    if mesh_indices[:4] != (expected_nghost, *expected_mesh_cells):
        raise ValueError("mesh RegionIndcs disagree with serialized parameters")
    if block_indices[:4] != (expected_nghost, *expected_block_cells):
        raise ValueError("MeshBlock RegionIndcs disagree with serialized parameters")
    if any(mesh % block for mesh, block in zip(
            expected_mesh_cells, expected_block_cells)):
        raise ValueError("root mesh dimensions are not divisible by MeshBlock dimensions")
    root_blocks = tuple(
        mesh // block for mesh, block in zip(expected_mesh_cells, expected_block_cells)
    )
    if root_level != _expected_root_level(root_blocks):
        raise ValueError("root logical level disagrees with root MeshBlock dimensions")

    expected_extents = (
        _float_parameter(parameters, "mesh", "x1min"),
        _float_parameter(parameters, "mesh", "x2min"),
        _float_parameter(parameters, "mesh", "x3min"),
        _float_parameter(parameters, "mesh", "x1max"),
        _float_parameter(parameters, "mesh", "x2max"),
        _float_parameter(parameters, "mesh", "x3max"),
    )
    if not all(_close(value, expected) for value, expected in zip(
            mesh_size[:6], expected_extents)):
        raise ValueError("mesh RegionSize disagrees with serialized parameters")
    expected_spacing = tuple(
        (expected_extents[axis + 3] - expected_extents[axis])
        / expected_mesh_cells[axis]
        for axis in range(3)
    )
    if not all(_close(value, expected) for value, expected in zip(
            mesh_size[6:], expected_spacing)):
        raise ValueError("mesh spacing in RegionSize is invalid")
    if not (math.isfinite(time) and math.isfinite(last_dt) and cycle >= 0):
        raise ValueError("restart time metadata is invalid")

    minimum_list_end = reader.tell() + _mesh_metadata_size(
        nmb_total, include_counters=False
    )
    if minimum_list_end > file_size:
        raise ValueError("restart is too small for its MeshBlock metadata")

    return FixedHeader(
        nmb_total=nmb_total,
        root_level=root_level,
        root_blocks=root_blocks,
        time=time,
        last_dt=last_dt,
        cycle=cycle,
    )


def _parse_layout(
        reader: AuditedReader, file_size: int, parameter_end: int,
        parameters: dict[str, dict[str, str]], byte_order: str,
        real_bytes: int) -> RestartMetadata:
    fixed = _parse_fixed_header(
        reader, file_size, parameter_end, parameters, byte_order, real_bytes
    )
    nmb_total = fixed.nmb_total
    root_level = fixed.root_level
    root_blocks = fixed.root_blocks
    time = fixed.time
    last_dt = fixed.last_dt
    cycle = fixed.cycle

    location_struct = struct.Struct(f"{byte_order}4i")
    location_bytes = _read_exact(
        reader, nmb_total * location_struct.size, "logical-location array"
    )
    locations = tuple(
        LogicalLocation(*values)
        for values in location_struct.iter_unpack(location_bytes)
    )
    cost_struct = struct.Struct(f"{byte_order}f")
    cost_bytes = _read_exact(
        reader, nmb_total * cost_struct.size, "MeshBlock cost array"
    )
    costs = tuple(values[0] for values in cost_struct.iter_unpack(cost_bytes))

    refinement = parameters.get("mesh_refinement", {}).get(
        "refinement", "none"
    )
    if refinement == "adaptive":
        maximum_physical_level = _integer_parameter(
            parameters, "mesh_refinement", "num_levels"
        ) - 1
        if maximum_physical_level < 0:
            raise ValueError("adaptive restart has fewer than one refinement level")
    elif refinement in ("static", "none"):
        # BuildTreeFromRestart initially permits every representable tree level and,
        # after reading the locations, sets max_level to the deepest stored level for
        # non-adaptive meshes.  In particular, valid static-SMR restarts do not carry
        # num_levels.
        maximum_physical_level = 31 - root_level
    else:
        raise ValueError(
            f"restart has invalid <mesh_refinement>/refinement: {refinement}"
        )
    seen: set[LogicalLocation] = set()
    for location in locations:
        physical_level = location.level - root_level
        if not 0 <= physical_level <= maximum_physical_level:
            raise ValueError(f"logical location has invalid level: {location}")
        limits = tuple(count << physical_level for count in root_blocks)
        if not all(0 <= coordinate < limit for coordinate, limit in zip(
                (location.lx1, location.lx2, location.lx3), limits)):
            raise ValueError(f"logical location is outside the mesh: {location}")
        if location in seen:
            raise ValueError(f"duplicate logical location: {location}")
        seen.add(location)
    if not all(math.isfinite(cost) and cost > 0.0 for cost in costs):
        raise ValueError("restart contains a non-positive or non-finite MeshBlock cost")

    # Match Mesh::BuildTreeFromRestart: the serialized marker, not a mutable
    # ParameterInput key, determines whether the optional extension is present.
    # On a legacy file the eight-byte probe is rolled back before object state.
    counters: tuple[int, ...] | None = None
    extension_offset = reader.tell()
    magic_struct = struct.Struct(f"{byte_order}Q")
    magic_bytes = reader.read(magic_struct.size)
    magic_present = (
        len(magic_bytes) == magic_struct.size
        and magic_struct.unpack(magic_bytes)[0] == AMR_COUNTER_MAGIC
    )
    if magic_present:
        _validate_metadata_capacity(nmb_total, include_counters=True)
        extension_struct = struct.Struct(f"{byte_order}ii")
        version, count = extension_struct.unpack(
            _read_exact(
                reader, extension_struct.size, "AMR counter extension header"
            )
        )
        if version != AMR_COUNTER_VERSION or count != nmb_total:
            raise ValueError(
                f"invalid AMR counter extension: version={version}, count={count}"
            )
        counter_struct = struct.Struct(f"{byte_order}i")
        counter_bytes = _read_exact(
            reader, nmb_total * counter_struct.size, "AMR cycle counters"
        )
        counters = tuple(
            values[0] for values in counter_struct.iter_unpack(counter_bytes)
        )
        if any(counter < 0 for counter in counters):
            raise ValueError("restart contains a negative AMR cycle counter")
    else:
        reader.seek(extension_offset)

    event_counters: RestartEventCounters | None = None
    event_extension_offset = reader.tell()
    event_magic_bytes = reader.read(magic_struct.size)
    event_magic_present = (
        len(event_magic_bytes) == magic_struct.size
        and magic_struct.unpack(event_magic_bytes)[0] == EVENT_COUNTER_MAGIC
    )
    if event_magic_present:
        _validate_metadata_capacity(
            nmb_total,
            include_counters=counters is not None,
            include_event_counters=True,
        )
        extension_struct = struct.Struct(f"{byte_order}ii")
        version, count = extension_struct.unpack(
            _read_exact(
                reader, extension_struct.size, "event-counter extension header"
            )
        )
        if version != EVENT_COUNTER_VERSION or count != EVENT_SUM_COUNTER_COUNT:
            raise ValueError(
                f"invalid event-counter extension: version={version}, count={count}"
            )
        sum_struct = struct.Struct(
            f"{byte_order}{EVENT_SUM_COUNTER_COUNT}Q"
        )
        sums = sum_struct.unpack(
            _read_exact(reader, sum_struct.size, "event sum counters")
        )
        maxit_struct = struct.Struct(f"{byte_order}i")
        maxit_c2p = maxit_struct.unpack(
            _read_exact(reader, maxit_struct.size, "event maximum C2P iterations")
        )[0]
        if maxit_c2p < 0:
            raise ValueError("restart contains a negative maximum C2P iteration count")
        event_counters = RestartEventCounters(*sums, maxit_c2p=maxit_c2p)
    else:
        reader.seek(event_extension_offset)

    return RestartMetadata(
        source="",
        file_size=file_size,
        parameter_end=parameter_end,
        metadata_end=reader.tell(),
        max_read_offset=reader.max_read_offset,
        byte_order="little" if byte_order == "<" else "big",
        real_bytes=real_bytes,
        nmb_total=nmb_total,
        root_level=root_level,
        time=time,
        last_dt=last_dt,
        cycle=cycle,
        locations=locations,
        costs=costs,
        amr_cycle_counters=counters,
        event_counters=event_counters,
        parameters=parameters,
    )


def read_restart_metadata_stream(
        stream: BinaryIO, file_size: int, source: str = "<stream>") -> RestartMetadata:
    """Read metadata from a seekable single-file restart stream.

    Apart from the eight-byte marker probes used by AthenaK itself for a legacy file,
    the function stops after the optional AMR age-counter and event-counter extensions.
    It never reads per-MeshBlock field arrays.
    """

    parameter_reader = AuditedReader(stream)
    parameter_reader.seek(0)
    parameter_raw, parameter_end = _read_parameter_dump(parameter_reader)
    parameters = _parse_parameters(parameter_raw)
    for block in parameters.values():
        if block.get("file_type") != "rst":
            continue
        text = block.get("single_file_per_rank", "false").strip().lower()
        if text in ("true", "1"):
            raise ValueError("per-rank restart files are not supported by this reader")
        if text not in ("false", "0"):
            raise ValueError(
                "invalid boolean value for restart single_file_per_rank: "
                f"{block.get('single_file_per_rank')}"
            )

    def candidate_reader() -> AuditedReader:
        # Construct after seeking so this candidate's high-water mark cannot inherit
        # bytes touched by a previous, rejected ABI interpretation.
        stream.seek(parameter_end)
        return AuditedReader(stream)

    fixed_failures: list[str] = []

    def fixed_candidates(real_bytes: int) -> list[str]:
        candidates: list[str] = []
        for byte_order in ("<", ">"):
            reader = candidate_reader()
            try:
                _parse_fixed_header(
                    reader, file_size, parameter_end, parameters,
                    byte_order, real_bytes
                )
            except (ValueError, struct.error) as exc:
                fixed_failures.append(f"{byte_order}/Real{real_bytes}: {exc}")
            else:
                candidates.append(byte_order)
        return candidates

    # A failed Real8 interpretation of a valid Real4 file has a longer fixed header
    # and can reach beyond the true metadata boundary.  Discriminate Real4 first using
    # only its short fixed header; inspect Real8 only when no Real4 endian is plausible.
    real_bytes = 4
    byte_orders = fixed_candidates(real_bytes)
    if not byte_orders:
        real_bytes = 8
        byte_orders = fixed_candidates(real_bytes)
    if len(byte_orders) != 1:
        detail = "; ".join(fixed_failures)
        raise ValueError(
            f"could not determine a unique native restart ABI "
            f"({len(byte_orders)} fixed-header candidates); {detail}"
        )

    byte_order = byte_orders[0]
    reader = candidate_reader()
    try:
        metadata = _parse_layout(
            reader, file_size, parameter_end, parameters, byte_order, real_bytes
        )
    except (ValueError, struct.error) as exc:
        raise ValueError(
            f"restart metadata is invalid for the selected native ABI "
            f"({byte_order}/Real{real_bytes}): {exc}"
        ) from exc
    return RestartMetadata(
        **{
            **metadata.__dict__,
            "source": source,
        }
    )


def read_restart_metadata(path: Path) -> RestartMetadata:
    """Read a stable local restart and reject a file that changes during inspection."""

    before = path.stat()
    with path.open("rb") as stream:
        metadata = read_restart_metadata_stream(stream, before.st_size, str(path))
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"restart changed while metadata was being read: {path}")
    return metadata


def load_balance(
        costs: tuple[float, ...] | list[float], nranks: int,
        max_blocks_per_rank: int = 0) -> list[RankPartition]:
    """Reproduce Mesh::LoadBalance's capacity-constrained contiguous partition."""

    block_count = len(costs)
    if nranks < 1:
        raise ValueError("cannot load balance over fewer than one rank")
    if block_count < nranks:
        raise ValueError("cannot keep at least one MeshBlock on every rank")
    if max_blocks_per_rank < 0:
        raise ValueError("maximum MeshBlocks per rank cannot be negative")
    if max_blocks_per_rank and block_count > nranks * max_blocks_per_rank:
        raise ValueError("MeshBlocks exceed the aggregate rank capacity")
    if not all(math.isfinite(cost) and cost > 0.0 for cost in costs):
        raise ValueError("all MeshBlock costs must be positive and finite")

    block_cap = max_blocks_per_rank if max_blocks_per_rank else block_count
    starts = [0] * nranks
    counts = [0] * nranks
    segment_end = block_count
    # Match C++'s explicit left-to-right binary64 += operations.  Python 3.12's sum()
    # uses compensated summation, which can move a closest-cost partition boundary.
    def ordered_sum(begin: int, end: int) -> float:
        total = 0.0
        for index in range(begin, end):
            total += float(costs[index])
        return total

    remaining_cost = ordered_sum(0, block_count)
    for rank in range(nranks - 1, 0, -1):
        minimum_begin = max(rank, segment_end - block_cap)
        maximum_begin = min(segment_end - 1, rank * block_cap)
        if minimum_begin > maximum_begin:
            raise ValueError(f"no feasible capacity boundary for rank {rank}")
        segment_begin = maximum_begin
        segment_cost = ordered_sum(segment_begin, segment_end)
        target_cost = remaining_cost / (rank + 1)
        while segment_begin > minimum_begin:
            candidate_cost = segment_cost + float(costs[segment_begin - 1])
            if abs(candidate_cost - target_cost) > abs(segment_cost - target_cost):
                break
            segment_begin -= 1
            segment_cost = candidate_cost
        starts[rank] = segment_begin
        counts[rank] = segment_end - segment_begin
        segment_end = segment_begin
        remaining_cost = max(0.0, remaining_cost - segment_cost)
    starts[0] = 0
    counts[0] = segment_end

    partitions = []
    for rank, (start, count) in enumerate(zip(starts, counts)):
        if count < 1 or count > block_cap:
            raise ValueError(f"invalid partition count for rank {rank}: {count}")
        partitions.append(RankPartition(
            rank=rank,
            gid_start=start,
            blocks=count,
            cost=ordered_sum(start, start + count),
        ))
    return partitions


def _default_capacity(metadata: RestartMetadata) -> int:
    text = metadata.parameters.get("mesh_refinement", {}).get(
        "max_nmb_per_rank", "0"
    )
    try:
        return int(text)
    except ValueError as exc:
        raise ValueError("invalid <mesh_refinement>/max_nmb_per_rank") from exc


def metadata_dict(
        metadata: RestartMetadata, nranks: int,
        max_blocks_per_rank: int) -> dict[str, object]:
    """Build the stable, compact CLI representation of parsed metadata."""

    partitions = load_balance(metadata.costs, nranks, max_blocks_per_rank)
    expected_costs = tuple(
        math.ldexp(1.0, location.level - metadata.root_level)
        for location in metadata.locations
    )
    return {
        "source": metadata.source,
        "file_size": metadata.file_size,
        "parameter_end": metadata.parameter_end,
        "metadata_end": metadata.metadata_end,
        "max_read_offset": metadata.max_read_offset,
        "field_data_read": False,
        "byte_order": metadata.byte_order,
        "real_bytes": metadata.real_bytes,
        "nmb_total": metadata.nmb_total,
        "root_level": metadata.root_level,
        "time": metadata.time,
        "last_dt": metadata.last_dt,
        "cycle": metadata.cycle,
        "level_subcycling_costs_exact": metadata.costs == expected_costs,
        "physical_levels": metadata.physical_level_rows(),
        "amr_cycle_counters": metadata.counter_rows(),
        "event_counters": (
            metadata.event_counters.as_dict()
            if metadata.event_counters is not None else None
        ),
        "partition": {
            "nranks": nranks,
            "max_blocks_per_rank": max_blocks_per_rank,
            "maximum_assigned_blocks": max(row.blocks for row in partitions),
            "minimum_capacity_headroom": (
                min(max_blocks_per_rank - row.blocks for row in partitions)
                if max_blocks_per_rank else None
            ),
            "ranks": [row.__dict__ for row in partitions],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("restart", type=Path)
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument(
        "--max-blocks-per-rank", type=int,
        help="partition capacity (default: value serialized in the restart)",
    )
    parser.add_argument("--json", action="store_true", help="emit compact JSON")
    args = parser.parse_args()

    try:
        metadata = read_restart_metadata(args.restart)
        capacity = (
            args.max_blocks_per_rank
            if args.max_blocks_per_rank is not None else _default_capacity(metadata)
        )
        payload = metadata_dict(metadata, args.ranks, capacity)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    if args.json:
        print(json.dumps(payload, sort_keys=True, separators=(",", ":")))
    else:
        print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
