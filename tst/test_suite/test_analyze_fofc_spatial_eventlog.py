"""Tests for the strict FOFC spatial event-log analyzer."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import analyze_fofc_spatial_eventlog as analyzer  # noqa: E402
from analyze_fofc_spatial_eventlog import (  # noqa: E402
    FofcTelemetryError,
    analyze_event_log,
    main,
)


HEADER = (
    "#  cycle eos_dfloor eos_efloor eos_tfloor eos_vceil eos_fail c2p_it "
    "fofc cons_adjust mag_adjust c2p_calls fofc_tests"
)
SCHEMA = (
    "# fofc_spatial_v1 kind=schema level_bins=0..31,overflow "
    "stage_bins=other,1,2,3 "
    "reason_bins=unknown,dmp_preflag,scalar,cons_density_floor,"
    "cons_energy_floor,prim_density_floor,prim_temperature_floor,"
    "rho_too_big,rho_too_small,nans_in_cons,mag_too_big,"
    "bracketing_failed,no_solution,invalid_geometry,other_c2p "
    "r_cyl_edges=2,4,8,16,32,64 abs_z_edges=0.5,1,2,4,8,16 "
    "lapse_edges=0.2,0.4,0.6,0.8,1 center1=0 center2=-1.5 center3=2"
)
SUMMARY = (
    "# fofc_spatial_v1 kind=summary cycle=323 count=15 nfofc=15 "
    "unattributed=0"
)
BINS = (
    "# fofc_spatial_v1 kind=bin cycle=323 level_bin=9 stage_bin=1 "
    "reason=dmp_preflag r_cyl_bin=0 abs_z_bin=1 lapse_bin=2 count=7",
    "# fofc_spatial_v1 kind=bin cycle=323 level_bin=9 stage_bin=2 "
    "reason=no_solution r_cyl_bin=0 abs_z_bin=1 lapse_bin=2 count=5",
    "# fofc_spatial_v1 kind=bin cycle=323 level_bin=8 stage_bin=2 "
    "reason=no_solution r_cyl_bin=3 abs_z_bin=4 lapse_bin=0 count=3",
)
EVENT = "323 0 1 2 0 0 14 15 4 5 100 200"
CANONICAL_PREFIX_BIN = (
    "# fofc_spatial_v1 kind=bin cycle=323 level_bin=32 stage_bin=0 "
    "reason=unknown r_cyl_bin=6 abs_z_bin=6 lapse_bin=5 count=4"
)


def _write_log(
    path: Path,
    *,
    schema: tuple[str, ...] = (SCHEMA,),
    summaries: tuple[str, ...] = (SUMMARY,),
    bins: tuple[str, ...] = BINS,
    events: tuple[str, ...] = (EVENT,),
    header: tuple[str, ...] = (HEADER,),
) -> Path:
    lines = ["# Athena event counter data", *header, *schema]
    lines.extend(summaries)
    lines.extend(bins)
    lines.extend(events)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _counts(result: dict[str, object], dimension: str) -> dict[str, int]:
    marginal = result["marginals"][dimension]
    return {row["label"]: row["count"] for row in marginal["bins"]}


def test_valid_cycle_reports_conservative_marginals_and_top_bins(
    tmp_path: Path,
) -> None:
    event_log = _write_log(tmp_path / "events.log")

    result = analyze_event_log(event_log, 323, top_bins=2)

    assert result["valid"]
    assert result["cycle"] == 323
    assert result["policy"] == {"require_unattributed_zero": False}
    assert result["summary"] == {"count": 15, "nfofc": 15, "unattributed": 0}
    assert result["event"]["fofc"] == 15
    assert result["schema"]["center"] == [0.0, -1.5, 2.0]
    assert result["checks"]["bin_sum"] == 15
    assert result["checks"]["marginal_sums"] == {
        "level": 15,
        "stage": 15,
        "reason": 15,
        "r_cyl": 15,
        "abs_z": 15,
        "lapse": 15,
    }
    assert _counts(result, "level") == {"8": 3, "9": 12}
    assert _counts(result, "stage") == {"1": 7, "2": 8}
    assert _counts(result, "reason") == {"dmp_preflag": 7, "no_solution": 8}
    assert _counts(result, "r_cyl") == {"[0,2)": 12, "[8,16)": 3}
    assert _counts(result, "abs_z") == {"[0.5,1)": 12, "[4,8)": 3}
    assert _counts(result, "lapse") == {"[0,0.2)": 3, "[0.4,0.6)": 12}
    assert [row["count"] for row in result["top_bins"]] == [7, 5]
    assert result["top_bins"][0]["reason"] == "dmp_preflag"
    assert result["global_ratios"]["fofc_per_test"]["value"] == 15 / 200
    assert result["global_ratios"]["cons_adjust_per_c2p_call"]["value"] == 4 / 100
    assert result["global_ratios"]["mag_adjust_per_c2p_call"]["value"] == 5 / 100
    assert all(row["defined"] for row in result["global_ratios"].values())
    assert result["input"]["stable_during_read"]
    assert len(result["input"]["sha256"]) == 64
    assert result["checks"]["all_cycles_ordered_and_closed"]
    assert result["records"]["complete_cycle_groups"] == 1


def test_zero_fofc_cycle_is_conservative_without_bin_records(tmp_path: Path) -> None:
    summary = SUMMARY.replace("count=15 nfofc=15", "count=0 nfofc=0")
    event = EVENT.replace("14 15 4", "14 0 4")
    event_log = _write_log(
        tmp_path / "zero.log", summaries=(summary,), bins=(), events=(event,)
    )

    result = analyze_event_log(event_log, 323)

    assert result["checks"]["bin_sum"] == 0
    assert result["top_bins"] == []
    assert all(value["bins"] == [] for value in result["marginals"].values())
    assert all(value["sum"] == 0 for value in result["marginals"].values())


def test_global_ratios_use_null_when_their_denominator_is_zero(
    tmp_path: Path,
) -> None:
    event = "323 0 1 2 0 0 14 15 4 5 0 0"
    event_log = _write_log(tmp_path / "zero_denominators.log", events=(event,))

    result = analyze_event_log(event_log, 323)

    for ratio in result["global_ratios"].values():
        assert ratio["denominator_value"] == 0
        assert ratio["value"] is None
        assert not ratio["defined"]


def test_nonfinite_spatial_overflow_labels_are_not_misreported_as_far_field(
    tmp_path: Path,
) -> None:
    histogram_bin = (
        "# fofc_spatial_v1 kind=bin cycle=323 level_bin=32 stage_bin=1 "
        "reason=invalid_geometry r_cyl_bin=6 abs_z_bin=6 lapse_bin=5 count=15"
    )
    event_log = _write_log(tmp_path / "nonfinite.log", bins=(histogram_bin,))

    result = analyze_event_log(event_log, 323)

    top = result["top_bins"][0]
    assert top["reason"] == "invalid_geometry"
    assert top["r_cyl"] == "other(r_cyl>=64_or_nonfinite)"
    assert top["abs_z"] == "other(abs_z>=16_or_nonfinite)"
    assert top["lapse"] == "other(lapse<0,lapse>=1,or_nonfinite)"
    assert _counts(result, "r_cyl") == {
        "other(r_cyl>=64_or_nonfinite)": 15
    }
    assert _counts(result, "abs_z") == {
        "other(abs_z>=16_or_nonfinite)": 15
    }


@pytest.mark.parametrize(
    ("schemas", "message"),
    [
        ((), "summary record precedes schema"),
        ((SCHEMA, SCHEMA), "duplicate schema record"),
    ],
)
def test_schema_must_be_unique(
    tmp_path: Path, schemas: tuple[str, ...], message: str
) -> None:
    event_log = _write_log(tmp_path / "schema.log", schema=schemas)

    with pytest.raises(FofcTelemetryError, match=message):
        analyze_event_log(event_log, 323)


def test_schema_contract_drift_is_rejected(tmp_path: Path) -> None:
    changed = SCHEMA.replace("lapse_edges=0.2,0.4,0.6,0.8,1", "lapse_edges=0.3,1")
    event_log = _write_log(tmp_path / "schema_drift.log", schema=(changed,))

    with pytest.raises(FofcTelemetryError, match="lapse_edges changed"):
        analyze_event_log(event_log, 323)


@pytest.mark.parametrize(
    ("summaries", "message"),
    [
        ((), "bin for cycle 323 has no active summary"),
        ((SUMMARY, SUMMARY), "missing event row before summary"),
    ],
)
def test_selected_summary_must_be_unique(
    tmp_path: Path, summaries: tuple[str, ...], message: str
) -> None:
    event_log = _write_log(tmp_path / "summary.log", summaries=summaries)

    with pytest.raises(FofcTelemetryError, match=message):
        analyze_event_log(event_log, 323)


@pytest.mark.parametrize(
    ("events", "message"),
    [
        ((), "incomplete telemetry group at EOF"),
        ((EVENT, EVENT), "event row for cycle 323 has no active summary"),
    ],
)
def test_selected_traditional_event_row_must_be_unique(
    tmp_path: Path, events: tuple[str, ...], message: str
) -> None:
    event_log = _write_log(tmp_path / "event.log", events=events)

    with pytest.raises(FofcTelemetryError, match=message):
        analyze_event_log(event_log, 323)


@pytest.mark.parametrize(
    ("summary", "event", "bins", "message"),
    [
        (
            SUMMARY.replace("count=15", "count=16"),
            EVENT,
            BINS,
            "summary.count=16 differs from summary.nfofc=15",
        ),
        (
            SUMMARY.replace("nfofc=15", "nfofc=16"),
            EVENT,
            BINS,
            "summary.count=15 differs from summary.nfofc=16",
        ),
        (
            SUMMARY,
            EVENT.replace("14 15 4", "14 16 4"),
            BINS,
            "summary.nfofc=15 differs from event.fofc=16",
        ),
        (
            SUMMARY,
            EVENT,
            BINS[:-1],
            "histogram sum=12 differs from event.fofc=15",
        ),
    ],
)
def test_authoritative_count_chain_is_fail_closed(
    tmp_path: Path,
    summary: str,
    event: str,
    bins: tuple[str, ...],
    message: str,
) -> None:
    event_log = _write_log(
        tmp_path / "counts.log",
        summaries=(summary,),
        bins=bins,
        events=(event,),
    )

    with pytest.raises(FofcTelemetryError, match=message):
        analyze_event_log(event_log, 323)


def test_first_cycle_canonical_unattributed_prefix_is_accepted_and_reported(
    tmp_path: Path,
) -> None:
    summary = SUMMARY.replace(
        "count=15 nfofc=15 unattributed=0",
        "count=19 nfofc=19 unattributed=4",
    )
    event = EVENT.replace("14 15 4", "14 19 4")
    event_log = _write_log(
        tmp_path / "unattributed.log",
        summaries=(summary,),
        bins=(*BINS, CANONICAL_PREFIX_BIN),
        events=(event,),
    )

    result = analyze_event_log(event_log, 323)

    assert result["valid"]
    assert result["summary"]["unattributed"] == 4
    assert not result["checks"]["unattributed_zero"]
    assert result["checks"]["unattributed_policy_satisfied"]
    assert result["checks"]["unattributed_canonical_bin_exact"]
    assert result["unattributed"] == {
        "count": 4,
        "present": True,
        "selected_cycle_is_first_telemetry_cycle": True,
        "permitted_only_on_first_telemetry_cycle": True,
        "canonical_bin": {
            "level_bin": 32,
            "stage_bin": 0,
            "reason": "unknown",
            "r_cyl_bin": 6,
            "abs_z_bin": 6,
            "lapse_bin": 5,
        },
        "canonical_bin_count": 4,
        "canonical_bin_equals_unattributed": True,
        "policy_satisfied": True,
    }
    assert result["checks"]["bin_sum"] == 19
    assert all(row["sum"] == 19 for row in result["marginals"].values())


def test_unattributed_prefix_outside_canonical_bin_is_rejected(
    tmp_path: Path,
) -> None:
    summary = SUMMARY.replace("unattributed=0", "unattributed=4")
    event_log = _write_log(tmp_path / "noncanonical_unattributed.log", summaries=(summary,))

    with pytest.raises(
        FofcTelemetryError,
        match=(
            "canonical overflow/other/unknown bin count=0 differs from "
            "unattributed=4"
        ),
    ):
        analyze_event_log(event_log, 323)


def test_canonical_prefix_bin_must_equal_not_merely_contain_unattributed(
    tmp_path: Path,
) -> None:
    summary = SUMMARY.replace(
        "count=15 nfofc=15 unattributed=0",
        "count=20 nfofc=20 unattributed=4",
    )
    canonical = CANONICAL_PREFIX_BIN.replace("count=4", "count=5")
    event = EVENT.replace("14 15 4", "14 20 4")
    event_log = _write_log(
        tmp_path / "canonical_excess.log",
        summaries=(summary,),
        bins=(*BINS, canonical),
        events=(event,),
    )

    with pytest.raises(
        FofcTelemetryError,
        match=(
            "canonical overflow/other/unknown bin count=5 differs from "
            "unattributed=4"
        ),
    ):
        analyze_event_log(event_log, 323)


def test_noncanonical_stage_other_bin_is_rejected_even_without_prefix(
    tmp_path: Path,
) -> None:
    invalid = BINS[0].replace("stage_bin=1", "stage_bin=0")
    event_log = _write_log(
        tmp_path / "noncanonical_stage_other.log",
        bins=(invalid, *BINS[1:]),
    )

    with pytest.raises(
        FofcTelemetryError,
        match="stage_bin=0 is reserved for the canonical unattributed prefix bin",
    ):
        analyze_event_log(event_log, 323)


def test_require_unattributed_zero_is_an_explicit_selected_cycle_policy(
    tmp_path: Path,
) -> None:
    summary = SUMMARY.replace(
        "count=15 nfofc=15 unattributed=0",
        "count=19 nfofc=19 unattributed=4",
    )
    event = EVENT.replace("14 15 4", "14 19 4")
    event_log = _write_log(
        tmp_path / "strict_unattributed.log",
        summaries=(summary,),
        bins=(*BINS, CANONICAL_PREFIX_BIN),
        events=(event,),
    )

    with pytest.raises(
        FofcTelemetryError,
        match="violates require_unattributed_zero policy",
    ):
        analyze_event_log(event_log, 323, require_unattributed_zero=True)


def test_duplicate_histogram_key_is_rejected_before_double_counting(
    tmp_path: Path,
) -> None:
    duplicate = BINS[0].replace("count=7", "count=1")
    event_log = _write_log(tmp_path / "duplicate_bin.log", bins=(*BINS, duplicate))

    with pytest.raises(FofcTelemetryError, match="duplicate histogram key"):
        analyze_event_log(event_log, 323)


@pytest.mark.parametrize(
    ("old", "new", "message"),
    [
        ("level_bin=9", "level_bin=33", "level_bin=33 is outside 0..32"),
        ("stage_bin=1", "stage_bin=4", "stage_bin=4 is outside 0..3"),
        ("r_cyl_bin=0", "r_cyl_bin=7", "r_cyl_bin=7 is outside 0..6"),
        ("abs_z_bin=1", "abs_z_bin=7", "abs_z_bin=7 is outside 0..6"),
        ("lapse_bin=2", "lapse_bin=6", "lapse_bin=6 is outside 0..5"),
        ("reason=dmp_preflag", "reason=bad_reason", "unknown reason bin"),
    ],
)
def test_histogram_key_ranges_are_strict(
    tmp_path: Path, old: str, new: str, message: str
) -> None:
    invalid = BINS[0].replace(old, new)
    event_log = _write_log(tmp_path / "range.log", bins=(invalid, *BINS[1:]))

    with pytest.raises(FofcTelemetryError, match=message):
        analyze_event_log(event_log, 323)


def test_duplicate_record_field_is_rejected(tmp_path: Path) -> None:
    summary = SUMMARY + " count=15"
    event_log = _write_log(tmp_path / "duplicate_field.log", summaries=(summary,))

    with pytest.raises(FofcTelemetryError, match="duplicate telemetry key 'count'"):
        analyze_event_log(event_log, 323)


def test_event_header_is_exact_and_unique(tmp_path: Path) -> None:
    event_log = _write_log(tmp_path / "headers.log", header=(HEADER, HEADER))

    with pytest.raises(FofcTelemetryError, match="duplicate event header"):
        analyze_event_log(event_log, 323)


def test_schema_cannot_precede_the_exact_event_header(tmp_path: Path) -> None:
    lines = ["# Athena event counter data", SCHEMA, HEADER, SUMMARY, *BINS, EVENT]
    event_log = tmp_path / "schema_before_header.log"
    event_log.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(
        FofcTelemetryError,
        match="schema record does not follow one exact event header",
    ):
        analyze_event_log(event_log, 323)


def test_event_header_cannot_appear_inside_an_active_group(tmp_path: Path) -> None:
    lines = [
        "# Athena event counter data",
        HEADER,
        SCHEMA,
        SUMMARY,
        HEADER,
        *BINS,
        EVENT,
    ]
    event_log = tmp_path / "header_inside_group.log"
    event_log.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(FofcTelemetryError, match="event header is out of order"):
        analyze_event_log(event_log, 323)


def test_symbolic_link_input_is_rejected(tmp_path: Path) -> None:
    target = _write_log(tmp_path / "real.log")
    link = tmp_path / "link.log"
    link.symlink_to(target)

    with pytest.raises(FofcTelemetryError, match="non-symlink"):
        analyze_event_log(link, 323)


def test_file_change_during_streaming_read_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    event_log = _write_log(tmp_path / "changing.log")
    original_parser = analyzer._parse_bin
    changed = False

    def mutate_after_bin(*args: object, **kwargs: object) -> dict[str, object]:
        nonlocal changed
        result = original_parser(*args, **kwargs)
        if not changed:
            with event_log.open("a", encoding="utf-8") as stream:
                stream.write("# changed while analyzer held its read descriptor\n")
            changed = True
        return result

    monkeypatch.setattr(analyzer, "_parse_bin", mutate_after_bin)

    with pytest.raises(FofcTelemetryError, match="changed while being read"):
        analyze_event_log(event_log, 323)


def test_missing_terminal_newline_is_rejected_as_a_truncated_log(
    tmp_path: Path,
) -> None:
    event_log = _write_log(tmp_path / "truncated.log")
    payload = event_log.read_bytes()
    assert payload.endswith(b"\n")
    event_log.write_bytes(payload[:-1])

    with pytest.raises(FofcTelemetryError, match="not terminated by a newline"):
        analyze_event_log(event_log, 323)


def test_newline_terminated_trailing_half_cycle_is_rejected(
    tmp_path: Path,
) -> None:
    event_log = _write_log(tmp_path / "trailing_half_cycle.log")
    trailing_summary = SUMMARY.replace("cycle=323", "cycle=324")
    trailing_bin = BINS[0].replace("cycle=323", "cycle=324").replace(
        "count=7", "count=15"
    )
    with event_log.open("a", encoding="utf-8") as stream:
        stream.write(trailing_summary + "\n")
        stream.write(trailing_bin + "\n")

    with pytest.raises(
        FofcTelemetryError,
        match="cycle 324: incomplete telemetry group at EOF; missing event row",
    ):
        analyze_event_log(event_log, 323)


def test_complete_later_cycle_can_be_selected_after_full_log_validation(
    tmp_path: Path,
) -> None:
    event_log = _write_log(tmp_path / "two_complete_cycles.log")
    second_summary = (
        "# fofc_spatial_v1 kind=summary cycle=324 count=0 nfofc=0 "
        "unattributed=0"
    )
    second_event = "324 0 0 0 0 0 0 0 0 0 0 0"
    with event_log.open("a", encoding="utf-8") as stream:
        stream.write(second_summary + "\n")
        stream.write(second_event + "\n")

    result = analyze_event_log(event_log, 324, require_unattributed_zero=True)

    assert result["cycle"] == 324
    assert result["records"]["complete_cycle_groups"] == 2
    assert result["records"]["selected_cycle_group_index"] == 1
    assert result["checks"]["all_cycles_ordered_and_closed"]
    assert result["checks"]["cycles_strictly_increasing"]
    assert result["checks"]["cycles_consecutive_under_dcycle_one"]
    assert result["checks"]["bin_sum"] == 0
    assert all(row["sum"] == 0 for row in result["marginals"].values())


def test_complete_groups_must_have_strictly_increasing_cycles(
    tmp_path: Path,
) -> None:
    event_log = _write_log(tmp_path / "decreasing_cycles.log")
    second_summary = (
        "# fofc_spatial_v1 kind=summary cycle=322 count=0 nfofc=0 "
        "unattributed=0"
    )
    second_event = "322 0 0 0 0 0 0 0 0 0 0 0"
    with event_log.open("a", encoding="utf-8") as stream:
        stream.write(second_summary + "\n")
        stream.write(second_event + "\n")

    with pytest.raises(
        FofcTelemetryError,
        match=(
            "telemetry cycles are not strictly increasing: "
            "previous=323, current=322"
        ),
    ):
        analyze_event_log(event_log, 323)


def test_complete_groups_must_be_consecutive_under_dcycle_one_contract(
    tmp_path: Path,
) -> None:
    event_log = _write_log(tmp_path / "cycle_gap.log")
    second_summary = (
        "# fofc_spatial_v1 kind=summary cycle=325 count=0 nfofc=0 "
        "unattributed=0"
    )
    second_event = "325 0 0 0 0 0 0 0 0 0 0 0"
    with event_log.open("a", encoding="utf-8") as stream:
        stream.write(second_summary + "\n")
        stream.write(second_event + "\n")

    with pytest.raises(
        FofcTelemetryError,
        match=(
            "telemetry cycles are not consecutive under the dcycle=1 contract: "
            "expected=324, current=325"
        ),
    ):
        analyze_event_log(event_log, 323)


def test_unattributed_prefix_is_rejected_after_first_complete_cycle(
    tmp_path: Path,
) -> None:
    event_log = _write_log(tmp_path / "later_unattributed.log")
    second_summary = (
        "# fofc_spatial_v1 kind=summary cycle=324 count=4 nfofc=4 "
        "unattributed=4"
    )
    second_bin = CANONICAL_PREFIX_BIN.replace("cycle=323", "cycle=324")
    second_event = "324 0 0 0 0 0 0 4 0 0 0 0"
    with event_log.open("a", encoding="utf-8") as stream:
        stream.write(second_summary + "\n")
        stream.write(second_bin + "\n")
        stream.write(second_event + "\n")

    with pytest.raises(
        FofcTelemetryError,
        match="unattributed=4 is permitted only in the first telemetry cycle",
    ):
        analyze_event_log(event_log, 323)


def test_next_summary_cannot_begin_before_prior_cycle_event(
    tmp_path: Path,
) -> None:
    second_summary = SUMMARY.replace("cycle=323", "cycle=324")
    event_log = _write_log(
        tmp_path / "summary_before_event.log",
        summaries=(SUMMARY, second_summary),
    )

    with pytest.raises(
        FofcTelemetryError,
        match="cycle 323: missing event row before summary for cycle 324",
    ):
        analyze_event_log(event_log, 323)


def test_cli_emits_compact_or_pretty_json(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    event_log = _write_log(tmp_path / "cli.log")

    assert main([str(event_log), "--cycle", "323", "--top-bins", "1", "--compact"]) == 0
    compact = capsys.readouterr()
    assert "\n" not in compact.out.rstrip("\n")
    assert json.loads(compact.out)["top_bins"][0]["count"] == 7

    assert main([str(event_log), "--cycle", "323", "--pretty"]) == 0
    pretty = capsys.readouterr()
    assert "\n  \"checks\"" in pretty.out
    assert json.loads(pretty.out)["valid"]


def test_cli_binds_and_enforces_unattributed_zero_policy(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    clean_log = _write_log(tmp_path / "clean_cli.log")
    assert main(
        [str(clean_log), "--cycle", "323", "--require-unattributed-zero"]
    ) == 0
    clean = json.loads(capsys.readouterr().out)
    assert clean["policy"] == {"require_unattributed_zero": True}
    assert clean["checks"]["unattributed_policy_satisfied"]

    summary = SUMMARY.replace(
        "count=15 nfofc=15 unattributed=0",
        "count=19 nfofc=19 unattributed=4",
    )
    event = EVENT.replace("14 15 4", "14 19 4")
    prefix_log = _write_log(
        tmp_path / "prefix_cli.log",
        summaries=(summary,),
        bins=(*BINS, CANONICAL_PREFIX_BIN),
        events=(event,),
    )
    assert main(
        [str(prefix_log), "--cycle", "323", "--require-unattributed-zero"]
    ) == 2
    rejected = capsys.readouterr()
    assert rejected.out == ""
    assert "require_unattributed_zero policy" in rejected.err


def test_cli_returns_two_for_invalid_telemetry(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    event_log = _write_log(tmp_path / "bad.log", summaries=())

    assert main([str(event_log), "--cycle", "323"]) == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "has no active summary" in captured.err
