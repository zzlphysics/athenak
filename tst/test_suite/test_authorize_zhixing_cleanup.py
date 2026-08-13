"""Fail-closed tests for the independent cloud-cleanup authorization gate."""

from __future__ import annotations

import base64
import importlib.util
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
SCRIPT = SCRIPTS / "authorize_zhixing_cleanup.py"
SPEC = importlib.util.spec_from_file_location("authorize_zhixing_cleanup", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
AUTHORIZER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = AUTHORIZER
SPEC.loader.exec_module(AUTHORIZER)


HOST_KEY = "SHA256:" + base64.b64encode(bytes(32)).decode("ascii").rstrip("=")


def valid_argv(tmp_path: Path) -> list[str]:
    return [
        "--remote-root", "/data/run01",
        "--remote-manifest", "/data/run01/manifests/segment-1.manifest.ready",
        "--expected-manifest-sha256", "a" * 64,
        "--segment", "segment-1",
        "--destination", str(tmp_path / "nas" / "segment-1"),
        "--expected-mount-source", "server:/volume1/Projects",
        "--expected-mount-fstype", "nfs",
        "--host-key-sha256", HOST_KEY,
        "--delete-file", "state/bin/output.00001.bin",
        "--authorization-dir", str(tmp_path / "authorizations"),
    ]


def test_explicit_binary_candidates_only() -> None:
    assert tuple(
        value.as_posix()
        for value in AUTHORIZER.validate_delete_candidates([
            "bin/global.00001.bin",
            "state/bin/local.00001.bin",
        ])
    ) == ("bin/global.00001.bin", "state/bin/local.00001.bin")


@pytest.mark.parametrize(
    "candidate",
    [
        "rst/run.00001.rst",
        "state/run.hst",
        "state/run.log",
        "evidence/segment.pass.ready",
        "manifests/segment.manifest.ready",
        "bin/*.bin",
        "../bin/escape.bin",
        "/bin/absolute.bin",
    ],
)
def test_protected_or_nonexact_candidate_is_rejected(candidate: str) -> None:
    with pytest.raises(ValueError):
        AUTHORIZER.validate_delete_candidates([candidate])


def test_duplicate_candidate_is_rejected() -> None:
    with pytest.raises(ValueError, match="duplicate"):
        AUTHORIZER.validate_delete_candidates([
            "bin/output.bin",
            "bin/output.bin",
        ])


def test_current_rw_nfs_topology_never_issues_authorization(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    authorization_dir = tmp_path / "authorizations"
    destination = tmp_path / "nas" / "segment-1"
    observed: list[tuple[Path, str, str]] = []

    def fake_verify_mount(path: Path, source: str, fstype: str) -> dict[str, str]:
        observed.append((path, source, fstype))
        return {
            "target": str(tmp_path / "nas"),
            "source": source,
            "fstype": "nfs",
            "options": "rw,hard,vers=3",
        }

    monkeypatch.setattr(AUTHORIZER, "verify_mount", fake_verify_mount)
    with pytest.raises(
        AUTHORIZER.CleanupAuthorizationUnavailable,
        match="multi-client read-write NFS.*transfer ACKs never authorize deletion",
    ):
        AUTHORIZER.main(valid_argv(tmp_path))

    assert observed == [(destination, "server:/volume1/Projects", "nfs")]
    assert not destination.exists()
    assert not authorization_dir.exists()


def test_no_unimplemented_topology_accidentally_opens_success_path() -> None:
    for mount in (
        {"fstype": "nfs4", "options": "rw"},
        {"fstype": "nfs", "options": "ro"},
        {"fstype": "ext4", "options": "rw"},
        {"fstype": "btrfs", "options": "ro"},
    ):
        with pytest.raises(AUTHORIZER.CleanupAuthorizationUnavailable):
            AUTHORIZER.require_supported_cleanup_backend(mount)


def test_gate_exposes_no_force_or_self_asserted_snapshot_escape_hatch() -> None:
    option_strings = {
        option
        for action in AUTHORIZER.build_parser()._actions
        for option in action.option_strings
    }
    assert "--force" not in option_strings
    assert "--assume-single-client" not in option_strings
    assert "--snapshot-id" not in option_strings


@pytest.mark.parametrize("ttl", [0, 29, 901, 3600])
def test_authorization_ttl_is_strictly_bounded(ttl: int) -> None:
    with pytest.raises(ValueError, match="TTL"):
        AUTHORIZER.validate_ttl(ttl)


def test_retention_policy_cannot_be_lowered_from_three_restarts() -> None:
    assert AUTHORIZER.MIN_RETAINED_RESTART_GENERATIONS == 3
    option_strings = {
        option
        for action in AUTHORIZER.build_parser()._actions
        for option in action.option_strings
    }
    assert "--minimum-retained-restarts" not in option_strings
    assert "--delete-restart" not in option_strings
