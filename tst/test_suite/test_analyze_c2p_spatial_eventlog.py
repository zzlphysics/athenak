from __future__ import annotations

from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import analyze_c2p_spatial_eventlog as analyzer  # noqa: E402
from analyze_c2p_spatial_eventlog import (  # noqa: E402
    C2PTelemetryError,
    analyze_lines,
    analyze_path,
)


SCHEMA = (
    "# c2p_spatial_v1 kind=schema "
    "intervention_bins=cons_adjust,mag_adjust level_bins=0..31,overflow "
    "r_cyl_edges=2,4,8,16,32,64 abs_z_edges=0.5,1,2,4,8,16 "
    "lapse_edges=0.2,0.4,0.6,0.8,1 "
    "density_floor_ratio_edges=1,2,4,16,64,256 "
    "magnetization_limit_ratio_edges=0.01,0.1,0.5,1,2,10 "
    "quantity_invalid_bin=7 stage_bins=other,1,2,3 "
    "geometry_bins=invalid,valid center1=0 center2=0 center3=0"
)
SCHEMA_V2 = (
    "# c2p_spatial_v2 kind=schema "
    "intervention_bins=cons_adjust,mag_adjust level_bins=0..31,overflow "
    "r_cyl_edges=2,4,8,16,32,64 abs_z_edges=0.5,1,2,4,8,16 "
    "lapse_edges=0.2,0.4,0.6,0.8,1 "
    "density_floor_ratio_edges=1,2,4,16,64,256,1e3,1e4,1e5,1e6,"
    "1e7,1e8,1e9,1e10 "
    "magnetization_limit_ratio_edges=0.01,0.1,0.5,1,2,10 "
    "density_floor_ratio_invalid_bin=15 "
    "magnetization_limit_ratio_invalid_bin=7 "
    "stage_bins=other,1,2,3 geometry_bins=invalid,valid "
    "center1=0 center2=0 center3=0"
)
HEADER = (
    "#  cycle eos_dfloor eos_efloor eos_tfloor eos_vceil eos_fail c2p_it "
    "fofc cons_adjust mag_adjust c2p_calls fofc_tests"
)


def _cycle(
    cycle: int, *, cons: int = 3, mag: int = 1,
    cons_unattributed: int = 0, mag_unattributed: int = 0,
) -> list[str]:
    lines = [
        f"# c2p_spatial_v1 kind=summary cycle={cycle} "
        f"intervention=cons_adjust count={cons} authoritative={cons} "
        f"unattributed={cons_unattributed}",
    ]
    if cons:
        if cons_unattributed:
            lines.append(
                f"# c2p_spatial_v1 kind=bin cycle={cycle} "
                "intervention=cons_adjust level_bin=32 r_cyl_bin=6 "
                "abs_z_bin=6 lapse_bin=5 density_floor_ratio_bin=7 "
                "magnetization_limit_ratio_bin=7 "
                f"count={cons_unattributed}"
            )
        attributed = cons-cons_unattributed
        if attributed:
            lines.append(
                f"# c2p_spatial_v1 kind=bin cycle={cycle} "
                "intervention=cons_adjust level_bin=11 r_cyl_bin=3 "
                "abs_z_bin=1 lapse_bin=4 density_floor_ratio_bin=6 "
                "magnetization_limit_ratio_bin=4 "
                f"count={attributed}"
            )
        if cons_unattributed:
            lines.append(
                f"# c2p_spatial_v1 kind=stage cycle={cycle} "
                f"intervention=cons_adjust stage_bin=0 count={cons_unattributed}"
            )
        if attributed:
            lines.append(
                f"# c2p_spatial_v1 kind=stage cycle={cycle} "
                f"intervention=cons_adjust stage_bin=1 count={attributed}"
            )
        if cons_unattributed:
            lines.append(
                f"# c2p_spatial_v1 kind=geometry cycle={cycle} "
                f"intervention=cons_adjust valid=0 count={cons_unattributed}"
            )
        if attributed:
            lines.append(
                f"# c2p_spatial_v1 kind=geometry cycle={cycle} "
                f"intervention=cons_adjust valid=1 count={attributed}"
            )
    lines.append(
        f"# c2p_spatial_v1 kind=summary cycle={cycle} "
        f"intervention=mag_adjust count={mag} authoritative={mag} "
        f"unattributed={mag_unattributed}"
    )
    if mag:
        if mag_unattributed:
            lines.extend([
                f"# c2p_spatial_v1 kind=bin cycle={cycle} "
                "intervention=mag_adjust level_bin=32 r_cyl_bin=6 "
                "abs_z_bin=6 lapse_bin=5 density_floor_ratio_bin=7 "
                "magnetization_limit_ratio_bin=7 "
                f"count={mag_unattributed}",
                f"# c2p_spatial_v1 kind=stage cycle={cycle} "
                f"intervention=mag_adjust stage_bin=0 count={mag_unattributed}",
                f"# c2p_spatial_v1 kind=geometry cycle={cycle} "
                f"intervention=mag_adjust valid=0 count={mag_unattributed}",
            ])
        attributed = mag-mag_unattributed
        if attributed:
            lines.extend([
                f"# c2p_spatial_v1 kind=bin cycle={cycle} "
                "intervention=mag_adjust level_bin=11 r_cyl_bin=3 "
                "abs_z_bin=1 lapse_bin=4 density_floor_ratio_bin=6 "
                "magnetization_limit_ratio_bin=4 "
                f"count={attributed}",
                f"# c2p_spatial_v1 kind=stage cycle={cycle} "
                f"intervention=mag_adjust stage_bin=2 count={attributed}",
                f"# c2p_spatial_v1 kind=geometry cycle={cycle} "
                f"intervention=mag_adjust valid=1 count={attributed}",
            ])
    lines.append(f"{cycle} 0 0 0 0 0 8 0 {cons} {mag} 100 100")
    return lines


def _valid_lines() -> list[str]:
    return [HEADER, SCHEMA, *_cycle(602)]


def test_valid_cycle_is_conservative_and_reports_physical_fractions() -> None:
    report = analyze_lines(_valid_lines(), cycle=602, source="fixture")
    assert report["selected_cycle"] == 602
    assert report["event"]["cons_adjust"] == 3
    cons = report["interventions"]["cons_adjust"]
    mag = report["interventions"]["mag_adjust"]
    assert cons["summary"]["count"] == 3
    assert cons["marginals"]["level_bin"] == {"11": 3}
    assert cons["derived_fractions"]["r_cyl_lt_16"] == 1.0
    assert cons["derived_fractions"]["density_ge_256_times_floor"] == 1.0
    assert mag["derived_fractions"]["magnetization_at_or_above_limit"] == 1.0


def test_v2_extended_density_bins_are_conservative() -> None:
    cycle = [
        line.replace("c2p_spatial_v1", "c2p_spatial_v2").replace(
            "density_floor_ratio_bin=6", "density_floor_ratio_bin=14")
        for line in _cycle(602)
    ]
    report = analyze_lines([HEADER, SCHEMA_V2, *cycle], cycle=602)
    assert report["schema"]["version"] == 2
    cons = report["interventions"]["cons_adjust"]
    assert cons["marginals"]["density_floor_ratio_bin"] == {"14": 3}
    assert cons["derived_fractions"]["density_ge_256_times_floor"] == 1.0


def test_zero_count_intervention_requires_no_fake_bins() -> None:
    report = analyze_lines([HEADER, SCHEMA, *_cycle(10, cons=0, mag=0)])
    for intervention in ("cons_adjust", "mag_adjust"):
        item = report["interventions"][intervention]
        assert item["summary"]["count"] == 0
        assert item["nonzero_joint_bins"] == 0


def test_first_cycle_restart_prefix_is_explicit_and_optionally_rejected() -> None:
    lines = [
        HEADER, SCHEMA,
        *_cycle(20, cons=3, mag=1, cons_unattributed=1, mag_unattributed=1),
        *_cycle(21),
    ]
    report = analyze_lines(lines, cycle=20)
    assert report["interventions"]["cons_adjust"]["summary"]["unattributed"] == 1
    with pytest.raises(C2PTelemetryError, match="unattributed telemetry is nonzero"):
        analyze_lines(lines, cycle=20, require_unattributed_zero=True)


def test_unattributed_prefix_is_forbidden_after_first_cycle() -> None:
    lines = [
        HEADER, SCHEMA, *_cycle(20),
        *_cycle(21, cons_unattributed=1, mag_unattributed=1),
    ]
    with pytest.raises(C2PTelemetryError, match="allowed only in first cycle"):
        analyze_lines(lines)


def test_missing_stage_conservation_is_rejected() -> None:
    lines = [
        line for line in _valid_lines()
        if not ("kind=stage" in line and "intervention=cons_adjust" in line)
    ]
    with pytest.raises(C2PTelemetryError, match="joint/stage/geometry totals"):
        analyze_lines(lines)


def test_duplicate_joint_key_is_rejected() -> None:
    lines = _valid_lines()
    duplicate = next(
        line for line in lines
        if "kind=bin" in line and "intervention=cons_adjust" in line
    )
    lines.insert(-1, duplicate)
    with pytest.raises(C2PTelemetryError, match="duplicate bin key"):
        analyze_lines(lines)


def test_event_counter_mismatch_is_rejected() -> None:
    lines = _valid_lines()
    lines[-1] = "602 0 0 0 0 0 8 0 2 1 100 100"
    with pytest.raises(C2PTelemetryError, match="summary count/authoritative/event"):
        analyze_lines(lines)


def test_streaming_path_reports_identity_and_sha256(tmp_path: Path) -> None:
    path = tmp_path / "event.log"
    path.write_text("\n".join(_valid_lines()) + "\n", encoding="utf-8")
    report = analyze_path(path, cycle=602, require_unattributed_zero=True)
    assert report["source_identity"]["size"] == path.stat().st_size
    assert len(report["source_identity"]["sha256"]) == 64


def test_symlink_input_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "event.log"
    link = tmp_path / "event-link.log"
    path.write_text("\n".join(_valid_lines()) + "\n", encoding="utf-8")
    link.symlink_to(path)
    with pytest.raises(C2PTelemetryError, match="cannot open read-only"):
        analyze_path(link)


def test_cli_returns_two_for_invalid_log(tmp_path: Path, capsys: Any) -> None:
    path = tmp_path / "invalid.log"
    path.write_text("# no telemetry\n", encoding="utf-8")
    assert analyzer.main([str(path)]) == 2
    assert "no supported C2P spatial schema found" in capsys.readouterr().err
