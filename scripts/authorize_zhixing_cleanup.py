#!/usr/bin/env python3
"""Fail-closed gate for any future Zhixing cloud cleanup authorization.

Transfer acknowledgements are durable receipts only.  They never authorize source
deletion.  The deployed destination is a multi-client read-write NFS mount, and this
repository currently has neither a verifiable storage snapshot/write barrier nor the
transactional per-file executor needed to consume an authorization safely.  Therefore
this program deliberately has no successful authorization path and never creates an
authorization artifact.

The command-line contract records the independently pinned inputs that a future
snapshot-backed implementation must bind.  Supplying those inputs is not proof of a
barrier and cannot override the fail-closed decision.
"""

from __future__ import annotations

import argparse
from pathlib import Path, PurePosixPath
from typing import NoReturn

from zhixing_pull_manifest import (
    DEFAULT_SSH_STATE,
    canonical_remote_absolute,
    normalize_host_key_sha256,
    safe_relative,
    verify_mount,
)


MIN_AUTHORIZATION_TTL_SECONDS = 30
MAX_AUTHORIZATION_TTL_SECONDS = 900
DEFAULT_AUTHORIZATION_TTL_SECONDS = 300
MIN_RETAINED_RESTART_GENERATIONS = 3


class CleanupAuthorizationUnavailable(RuntimeError):
    """The current storage/executor contract cannot safely authorize deletion."""


def validate_segment(value: str) -> str:
    if not value or any(
        char not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_."
        for char in value
    ):
        raise ValueError("segment may contain only letters, digits, '-', '_' and '.'")
    return value


def validate_sha256(value: str, label: str) -> str:
    if len(value) != 64 or any(char not in "0123456789abcdef" for char in value):
        raise ValueError(f"{label} must be 64 lowercase hex digits")
    return value


def validate_ttl(value: int) -> int:
    if not MIN_AUTHORIZATION_TTL_SECONDS <= value <= MAX_AUTHORIZATION_TTL_SECONDS:
        raise ValueError(
            "authorization TTL must be between "
            f"{MIN_AUTHORIZATION_TTL_SECONDS} and "
            f"{MAX_AUTHORIZATION_TTL_SECONDS} seconds"
        )
    return value


def validate_delete_candidates(values: list[str]) -> tuple[PurePosixPath, ...]:
    """Validate an exact, conservative delete set without expanding any path.

    Even a future snapshot-backed gate may initially authorize only binary dumps.
    Restarts, histories, logs, evidence, manifests, ACKs, and ready records remain
    retained.  This is stronger than the minimum-three-restart policy.
    """

    if not values:
        raise ValueError("at least one explicit --delete-file is required")
    candidates: list[PurePosixPath] = []
    seen: set[str] = set()
    for value in values:
        if any(character in value for character in "*?[]"):
            raise ValueError(f"delete candidates cannot contain glob syntax: {value!r}")
        candidate = safe_relative(value)
        text = candidate.as_posix()
        if text in seen:
            raise ValueError(f"duplicate delete candidate: {text}")
        seen.add(text)
        parts = candidate.parts
        in_binary_tree = parts[0] == "bin" or (
            len(parts) >= 2 and parts[0] == "state" and parts[1] == "bin"
        )
        if not in_binary_tree or candidate.suffix != ".bin":
            raise ValueError(
                "cleanup candidates are restricted to explicit bin/**/*.bin files; "
                f"retaining protected output {text!r}"
            )
        candidates.append(candidate)
    return tuple(candidates)


def require_supported_cleanup_backend(mount: dict[str, str]) -> NoReturn:
    """Reject until a verifiable snapshot backend and safe consumer both exist.

    A read-only mode bit, a local ``/proc`` scan, an advisory lock, or a read-only
    remount of the same live NFS export cannot revoke an already-open writable handle
    on another NFS client.  No command-line assertion is accepted as a substitute.
    """

    fstype = mount.get("fstype", "")
    options = set(mount.get("options", "").split(","))
    if fstype in {"nfs", "nfs4"} and "rw" in options:
        detail = (
            "multi-client read-write NFS has no verifiable server-side "
            "snapshot/write barrier"
        )
    else:
        detail = "no verifiable read-only snapshot/write-barrier backend is implemented"
    raise CleanupAuthorizationUnavailable(
        "cleanup authorization disabled: "
        f"{detail}; no transactional per-file authorization consumer is implemented; "
        "transfer ACKs never authorize deletion"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--remote-root", required=True)
    parser.add_argument("--remote-manifest", required=True)
    parser.add_argument("--expected-manifest-sha256", required=True)
    parser.add_argument("--segment", required=True)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--expected-mount-source", required=True)
    parser.add_argument("--expected-mount-fstype", default="nfs")
    parser.add_argument("--ssh-state", type=Path, default=DEFAULT_SSH_STATE)
    parser.add_argument("--host-key-sha256", required=True)
    parser.add_argument(
        "--delete-file",
        action="append",
        required=True,
        help="exact relative binary-dump path; repeat for each candidate",
    )
    parser.add_argument(
        "--authorization-dir",
        required=True,
        type=Path,
        help="future immutable authorization directory (never written by this version)",
    )
    parser.add_argument(
        "--ttl-seconds",
        type=int,
        default=DEFAULT_AUTHORIZATION_TTL_SECONDS,
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        validate_segment(args.segment)
        validate_sha256(
            args.expected_manifest_sha256,
            "--expected-manifest-sha256",
        )
        validate_ttl(args.ttl_seconds)
        validate_delete_candidates(args.delete_file)
        normalize_host_key_sha256(args.host_key_sha256)
        canonical_remote_absolute(args.remote_root, "remote root")
        canonical_remote_absolute(args.remote_manifest, "remote manifest")
    except ValueError as exc:
        parser.error(str(exc))

    # This is deliberately read-only.  In particular, do not create destination,
    # authorization-dir, temporary files, claims, or any remote connection before
    # proving that a supported storage barrier exists.
    destination = args.destination.expanduser().resolve()
    mount = verify_mount(
        destination,
        args.expected_mount_source,
        args.expected_mount_fstype,
    )
    require_supported_cleanup_backend(mount)


if __name__ == "__main__":
    raise SystemExit(main())
