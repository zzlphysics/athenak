#!/usr/bin/env python3
"""One-shot, fail-closed 120-hour renewal for Zhixing instance 641496.

This program deliberately separates the non-idempotent renewal call from all
retries.  ``renew-once`` durably writes an intent before its sole renewal API
call.  Once that intent exists, both ``renew-once`` and ``reconcile`` use only
the read-only instance-list and instance-detail endpoints.  A timeout or lost
HTTP response can therefore never cause an accidental second purchase.

The instance-scoped mutation lock is also the lock contract for subsequent
disk-expansion and stop transactions.  Those transactions must acquire
``INSTANCE_MUTATION_LOCK`` before making any mutating API call.

Only explicitly whitelisted, non-secret fields are persisted.  In particular,
raw API responses are never written because instance responses can contain SSH
credentials and passwords.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, replace
from decimal import Decimal, InvalidOperation
import fcntl
import hashlib
import inspect
import json
import math
import os
from pathlib import Path
import re
import stat
import sys
import time
import types
from typing import Any, Callable, Dict, Iterator, List, Optional, Sequence, Tuple


EXPECTED_INSTANCE_ID = 641496
EXPECTED_GPU_TYPE = "Tesla V100-32GB"
EXPECTED_GPU_NUM = 8
EXPECTED_ADD_DISK_GB = 4500
EXPECTED_ADDING_DISK_GB = 0
EXPECTED_STATUS = 1
EXPECTED_OLD_DUE_TIME = 1787310616  # 2026-08-21 19:10:16 +08:00
PACK_HOURS = 120
EXTENSION_SECONDS = PACK_HOURS * 3600
EXPECTED_NEW_DUE_TIME = 1787742616  # 2026-08-26 19:10:16 +08:00
MAX_COST_CNY = Decimal("2000")
EXPECTED_FINAL_COST_CNY = Decimal("1589.376")
EXPECTED_QUOTE_COMPONENTS = {
    "gpu": {
        "type": "Tesla V100-32GB",
        "contract": "Tesla V100-32GB",
        "title": "租用GPU:8卡",
        "single_hour_price_CNY": Decimal("1.32"),
        "amount": 8,
        "final_discount": Decimal("0.98"),
        "final_cost_CNY": Decimal("1241.856"),
    },
    "spec": {
        "type": "spec",
        "contract": "gcs.g8.xlarge",
        "title": "实例规格gcs.g8.xlarge",
        "single_hour_price_CNY": Decimal("0.2"),
        "amount": 1,
        "final_discount": Decimal("0.98"),
        "final_cost_CNY": Decimal("23.52"),
    },
    "disk": {
        "type": "add_disk",
        "contract": "add_disk:4500GB",
        "title": "添加/扩容数据盘:4500GB",
        "single_hour_price_CNY": Decimal("0.0006"),
        "amount": 4500,
        "final_discount": Decimal("1"),
        "final_cost_CNY": Decimal("324"),
    },
}

API_BASE_SOURCE = Path("/tmp/zhixing_l4_control.py")
API_BASE_SHA256 = (
    "79e0070e8af72bf359c9a770ae901bcdd9fc9e5424d937aa106eaf52edd8bba5"
)
PRIVATE_STATE_DIR = Path(
    "/home/zhangzelin/kx-4t/athenak-zhixing-private/641496"
)
API_BASE = PRIVATE_STATE_DIR / (
    "zhixing_l4_control."
    "79e0070e8af72bf359c9a770ae901bcdd9fc9e5424d937aa106eaf52edd8bba5.py"
)
API_CREDENTIALS_SOURCE = Path("/tmp/zhixing_api.credentials")
API_CREDENTIALS = PRIVATE_STATE_DIR / "zhixing_api.credentials"
INTENT = Path(
    PRIVATE_STATE_DIR
    / "zhixing-v100-641496-renew120-from-1787310616.intent.private.json"
)
RESULT = Path(
    PRIVATE_STATE_DIR
    / "zhixing-v100-641496-renew120-from-1787310616.result.private.json"
)
CALL_MARKER = Path(
    PRIVATE_STATE_DIR
    / "zhixing-v100-641496-renew120-from-1787310616.call.private.json"
)
PRIOR_ANCHOR = PRIVATE_STATE_DIR / "prior-renewal-chain.anchor.private.json"

# This exact path is the cross-operation lock contract.  A future AddDisk or
# stop transaction for instance 641496 must use the same lock path.
INSTANCE_MUTATION_LOCK = PRIVATE_STATE_DIR / "instance-641496.mutation.lock"

QUOTE_ENDPOINT = "/instance/renew_instance_query_price"
RENEW_ENDPOINT = "/instance/renew_instance"
LIST_ENDPOINT = "/instance/get_instance_list"
DETAIL_ENDPOINT = "/instance/get_instance_detail"
EXPECTED_API_BASE_URL = "https://app.ai-galaxy.cn/openapi/v2"
INSTANCE_NAME_PATTERN = re.compile(r"\A[A-Za-z0-9][A-Za-z0-9_.-]{0,255}\Z")

OBSERVE_SECONDS = 120.0
OBSERVE_INTERVAL_SECONDS = 2.0
MAX_PRIVATE_FILE_BYTES = 1 << 20


class ContractViolation(RuntimeError):
    """The local evidence, live instance, or quote violates the contract."""


@dataclass(frozen=True)
class PriorRenewal:
    path: Path
    sha256: str
    old_due_time: int
    new_due_time: int
    extension_seconds: int
    pack_hours: int


@dataclass(frozen=True)
class Config:
    control: Path
    control_sha256: str
    control_source: Path
    credentials: Path
    credentials_source: Path
    intent: Path
    result: Path
    call_marker: Path
    prior_anchor: Path
    mutation_lock: Path
    prior_renewals: Tuple[PriorRenewal, ...]
    instance_id: int = EXPECTED_INSTANCE_ID
    gpu_type: str = EXPECTED_GPU_TYPE
    gpu_num: int = EXPECTED_GPU_NUM
    add_disk_GB: int = EXPECTED_ADD_DISK_GB
    adding_disk_GB: int = EXPECTED_ADDING_DISK_GB
    can_update_disk: bool = True
    status: int = EXPECTED_STATUS
    old_due_time: int = EXPECTED_OLD_DUE_TIME
    pack_hours: int = PACK_HOURS
    new_due_time: int = EXPECTED_NEW_DUE_TIME
    max_cost_CNY: Decimal = MAX_COST_CNY
    observe_seconds: float = OBSERVE_SECONDS
    observe_interval_seconds: float = OBSERVE_INTERVAL_SECONDS
    require_persistent_xfs: bool = True


DEFAULT_CONFIG = Config(
    control=API_BASE,
    control_sha256=API_BASE_SHA256,
    control_source=API_BASE_SOURCE,
    credentials=API_CREDENTIALS,
    credentials_source=API_CREDENTIALS_SOURCE,
    intent=INTENT,
    result=RESULT,
    call_marker=CALL_MARKER,
    prior_anchor=PRIOR_ANCHOR,
    mutation_lock=INSTANCE_MUTATION_LOCK,
    prior_renewals=(
        PriorRenewal(
            path=Path("/tmp/zhixing_v100x8_mid_renew120.private.json"),
            sha256=(
                "b0bf9f291fc53f563fa3ff8b7d88079f012d79ac22969073da5535b9a2c12f77"
            ),
            old_due_time=1786734616,
            new_due_time=1787166616,
            extension_seconds=432000,
            pack_hours=120,
        ),
        PriorRenewal(
            path=Path("/tmp/zhixing_v100x8_mid_renew40.private.json"),
            sha256=(
                "2594df1bec97f1b1e0c40129234aa1a2ad43a943640d49819e34bf45c182f5c6"
            ),
            old_due_time=1787166616,
            new_due_time=1787310616,
            extension_seconds=144000,
            pack_hours=40,
        ),
    ),
)


@dataclass(frozen=True)
class Timers:
    monotonic: Callable[[], float] = time.monotonic
    sleep: Callable[[float], None] = time.sleep
    wall_time: Callable[[], float] = time.time


def _exact_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ContractViolation(f"{label} is not an exact JSON integer")
    return value


def _finite_decimal(value: Any, label: str) -> Decimal:
    if isinstance(value, bool):
        raise ContractViolation(f"{label} is not a decimal number")
    try:
        result = Decimal(str(value))
    except (InvalidOperation, TypeError, ValueError) as error:
        raise ContractViolation(f"{label} is not a decimal number") from error
    if not result.is_finite():
        raise ContractViolation(f"{label} is not finite")
    return result


def _path_occupied(path: Path) -> bool:
    try:
        os.lstat(path)
    except FileNotFoundError:
        return False
    return True


def _unescape_mount_field(value: str) -> str:
    for encoded, decoded in (
        ("\\040", " "),
        ("\\011", "\t"),
        ("\\012", "\n"),
        ("\\134", "\\"),
    ):
        value = value.replace(encoded, decoded)
    return value


def _mount_contract(path: Path) -> Dict[str, Any]:
    resolved = str(path.resolve(strict=True))
    matches: List[Tuple[int, str, str, bool]] = []
    for line in Path("/proc/self/mountinfo").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        try:
            separator = fields.index("-")
            mountpoint = _unescape_mount_field(fields[4])
            options = fields[5].split(",")
            fstype = fields[separator + 1]
            source = _unescape_mount_field(fields[separator + 2])
        except (IndexError, ValueError) as error:
            raise ContractViolation("unrecognized /proc/self/mountinfo record") from error
        prefix = mountpoint.rstrip("/") + "/"
        if resolved == mountpoint or resolved.startswith(prefix):
            matches.append((len(mountpoint), fstype, source, "rw" in options))
    if not matches:
        raise ContractViolation("private state directory has no mount record")
    _, fstype, source, read_write = max(matches)
    return {
        "fstype": fstype,
        "source": source,
        "read_write": read_write,
    }


def validate_transaction_paths(config: Config) -> None:
    paths = (
        config.intent,
        config.result,
        config.call_marker,
        config.prior_anchor,
        config.mutation_lock,
        config.control,
        config.credentials,
    )
    if not all(path.is_absolute() for path in paths):
        raise ContractViolation("transaction evidence paths must be absolute")
    parents = {path.parent for path in paths}
    if len(parents) != 1:
        raise ContractViolation("transaction evidence must share one state directory")
    parent = paths[0].parent
    parent_stat = os.lstat(parent)
    if not stat.S_ISDIR(parent_stat.st_mode) or stat.S_ISLNK(parent_stat.st_mode):
        raise ContractViolation("private state path is not a real directory")
    if parent_stat.st_uid != os.geteuid():
        raise ContractViolation("private state directory has the wrong owner")
    if stat.S_IMODE(parent_stat.st_mode) != 0o700:
        raise ContractViolation("private state directory is not mode 0700")
    if config.require_persistent_xfs:
        if parent != PRIVATE_STATE_DIR:
            raise ContractViolation("live transaction is not on the fixed state path")
        mount = _mount_contract(parent)
        if mount["fstype"] != "xfs" or mount["read_write"] is not True:
            raise ContractViolation("private state directory is not on read-write XFS")
    if not config.control_source.is_absolute() or not config.credentials_source.is_absolute():
        raise ContractViolation("runtime installation sources must be absolute")
    if re.fullmatch(r"[0-9a-f]{64}", config.control_sha256) is None:
        raise ContractViolation("API base controller SHA-256 is invalid")


def _secure_read_bytes(path: Path) -> bytes:
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW | os.O_NONBLOCK
    descriptor = os.open(path, flags)
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ContractViolation(f"private evidence is not a regular file: {path}")
        if before.st_uid != os.geteuid():
            raise ContractViolation(f"private evidence has the wrong owner: {path}")
        if before.st_nlink != 1:
            raise ContractViolation(
                f"private evidence must have exactly one hard link: {path}"
            )
        if stat.S_IMODE(before.st_mode) != 0o600:
            raise ContractViolation(f"private evidence is not mode 0600: {path}")
        if before.st_size <= 0 or before.st_size > MAX_PRIVATE_FILE_BYTES:
            raise ContractViolation(f"private evidence size is unsafe: {path}")
        chunks: List[bytes] = []
        remaining = before.st_size
        while remaining:
            chunk = os.read(descriptor, min(remaining, 65536))
            if not chunk:
                raise ContractViolation(f"private evidence was truncated: {path}")
            chunks.append(chunk)
            remaining -= len(chunk)
        if os.read(descriptor, 1):
            raise ContractViolation(f"private evidence grew while read: {path}")
        after = os.fstat(descriptor)
        stable_fields = (
            "st_dev",
            "st_ino",
            "st_size",
            "st_mtime_ns",
            "st_ctime_ns",
            "st_mode",
            "st_uid",
            "st_gid",
            "st_nlink",
        )
        if any(getattr(before, field) != getattr(after, field)
               for field in stable_fields):
            raise ContractViolation(f"private evidence changed while read: {path}")
        payload = b"".join(chunks)
        if not payload.endswith(b"\n"):
            raise ContractViolation(f"private evidence lacks terminal newline: {path}")
        return payload
    finally:
        os.close(descriptor)


def _secure_read_json(path: Path) -> Tuple[Dict[str, Any], bytes]:
    payload = _secure_read_bytes(path)
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractViolation(f"private evidence is not valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ContractViolation(f"private evidence root is not an object: {path}")
    return value, payload


def write_bytes_exclusive(path: Path, payload: bytes) -> None:
    """Create exactly one durable 0600 file without following symlinks."""

    if not isinstance(payload, bytes) or not payload:
        raise ContractViolation("exclusive file payload must be non-empty bytes")
    descriptor = os.open(
        path,
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | os.O_NOFOLLOW
        | os.O_CLOEXEC,
        0o600,
    )
    try:
        os.fchmod(descriptor, 0o600)
        offset = 0
        while offset < len(payload):
            written = os.write(descriptor, payload[offset:])
            if written <= 0:
                raise OSError("short private evidence write")
            offset += written
        os.fsync(descriptor)
        created = os.fstat(descriptor)
        if not stat.S_ISREG(created.st_mode):
            raise ContractViolation("created private evidence is not a regular file")
        if stat.S_IMODE(created.st_mode) != 0o600:
            raise ContractViolation("created private evidence is not mode 0600")
    finally:
        os.close(descriptor)
    parent = os.open(
        path.parent,
        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC,
    )
    try:
        os.fsync(parent)
    finally:
        os.close(parent)


def write_private_exclusive(path: Path, value: Dict[str, Any]) -> None:
    """Create exactly one durable 0600 JSON file without following symlinks."""

    payload = (json.dumps(
        value, ensure_ascii=False, indent=2, sort_keys=True
    ) + "\n").encode("utf-8")
    write_bytes_exclusive(path, payload)


@contextmanager
def acquire_instance_mutation_lock(path: Path) -> Iterator[int]:
    """Acquire the lock shared by renew, AddDisk, and stop for one instance."""

    descriptor = os.open(
        path,
        os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW | os.O_CLOEXEC,
        0o600,
    )
    try:
        lock_stat = os.fstat(descriptor)
        if not stat.S_ISREG(lock_stat.st_mode):
            raise ContractViolation("instance mutation lock is not a regular file")
        if lock_stat.st_uid != os.geteuid():
            raise ContractViolation("instance mutation lock has the wrong owner")
        os.fchmod(descriptor, 0o600)
        os.fsync(descriptor)
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ContractViolation("instance mutation lock is busy") from error
        yield descriptor
    finally:
        os.close(descriptor)


def _prior_safe_summary(contract: PriorRenewal) -> Dict[str, Any]:
    return {
        "evidence_file": contract.path.name,
        "evidence_sha256": contract.sha256,
        "old_due_time": contract.old_due_time,
        "new_due_time": contract.new_due_time,
        "extension_seconds": contract.extension_seconds,
        "pack_hours": contract.pack_hours,
    }


def validate_prior_renewal_chain(config: Config) -> List[Dict[str, Any]]:
    if not config.prior_renewals:
        raise ContractViolation("prior renewal chain is empty")
    summaries: List[Dict[str, Any]] = []
    previous_new: Optional[int] = None
    for contract in config.prior_renewals:
        value, payload = _secure_read_json(contract.path)
        actual_sha = hashlib.sha256(payload).hexdigest()
        if actual_sha != contract.sha256:
            raise ContractViolation(
                f"prior renewal SHA-256 changed: {contract.path.name}"
            )
        expected = {
            "instance_id": config.instance_id,
            "old_due_time": contract.old_due_time,
            "new_due_time": contract.new_due_time,
            "extension_seconds": contract.extension_seconds,
            "pack_hours": contract.pack_hours,
        }
        for field, wanted in expected.items():
            if _exact_int(value.get(field), f"prior renewal {field}") != wanted:
                raise ContractViolation(
                    f"prior renewal field changed: {contract.path.name}:{field}"
                )
        if contract.new_due_time - contract.old_due_time != contract.extension_seconds:
            raise ContractViolation("prior renewal contract has an invalid delta")
        if contract.extension_seconds != contract.pack_hours * 3600:
            raise ContractViolation("prior renewal contract has invalid hours")
        if previous_new is not None and contract.old_due_time != previous_new:
            raise ContractViolation("prior renewal evidence is not a continuous chain")
        previous_new = contract.new_due_time
        summaries.append(_prior_safe_summary(contract))
    if previous_new != config.old_due_time:
        raise ContractViolation("prior renewal chain does not end at the baseline due time")
    return summaries


def _make_prior_anchor(
    config: Config,
    summaries: Sequence[Dict[str, Any]],
    timers: Timers,
) -> Dict[str, Any]:
    return {
        "schema": "zhixing-renew-prior-anchor-v1",
        "status": "legacy_chain_whitelist_anchored",
        "instance_id": config.instance_id,
        "prior_renewal_chain": list(summaries),
        "terminal_due_time": config.old_due_time,
        "source_policy": "pinned legacy SHA-256; raw API responses not copied",
        "created_unix": timers.wall_time(),
    }


def _validate_prior_anchor(
    config: Config,
) -> Tuple[List[Dict[str, Any]], bytes]:
    anchor, payload = _secure_read_json(config.prior_anchor)
    _validate_exact_keys(anchor, (
        "schema",
        "status",
        "instance_id",
        "prior_renewal_chain",
        "terminal_due_time",
        "source_policy",
        "created_unix",
    ), "prior renewal anchor")
    if anchor.get("schema") != "zhixing-renew-prior-anchor-v1":
        raise ContractViolation("prior renewal anchor schema changed")
    if anchor.get("status") != "legacy_chain_whitelist_anchored":
        raise ContractViolation("prior renewal anchor status changed")
    if _exact_int(
        anchor.get("instance_id"), "prior anchor instance_id"
    ) != config.instance_id:
        raise ContractViolation("prior renewal anchor instance changed")
    if _exact_int(
        anchor.get("terminal_due_time"), "prior anchor terminal_due_time"
    ) != config.old_due_time:
        raise ContractViolation("prior renewal anchor due time changed")
    if anchor.get("source_policy") != (
        "pinned legacy SHA-256; raw API responses not copied"
    ):
        raise ContractViolation("prior renewal anchor source policy changed")
    expected = [_prior_safe_summary(item) for item in config.prior_renewals]
    _validate_intent_prior_chain(anchor.get("prior_renewal_chain"), expected)
    created = anchor.get("created_unix")
    if isinstance(created, bool) or not isinstance(created, (int, float)):
        raise ContractViolation("prior renewal anchor created_unix is invalid")
    if not math.isfinite(float(created)):
        raise ContractViolation("prior renewal anchor created_unix is not finite")
    return expected, payload


def load_or_create_prior_anchor(
    config: Config, timers: Timers
) -> Tuple[List[Dict[str, Any]], bytes]:
    if _path_occupied(config.prior_anchor):
        return _validate_prior_anchor(config)
    summaries = validate_prior_renewal_chain(config)
    anchor = _make_prior_anchor(config, summaries, timers)
    write_private_exclusive(config.prior_anchor, anchor)
    return _validate_prior_anchor(config)


def _validate_credentials_payload(payload: bytes) -> Tuple[str, str]:
    try:
        values = payload.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ContractViolation("API credentials are not UTF-8") from error
    if len(values) != 2:
        raise ContractViolation("API credentials must contain exactly two lines")
    for value in values:
        if re.fullmatch(r"[A-Za-z0-9_-]{8,128}", value) is None:
            raise ContractViolation("API credential line violates the whitelist")
    return values[0], values[1]


def install_runtime(config: Config = DEFAULT_CONFIG) -> Dict[str, Any]:
    """Install pinned controller bytes and credentials on persistent XFS."""

    validate_transaction_paths(config)
    with acquire_instance_mutation_lock(config.mutation_lock):
        if _path_occupied(config.control):
            controller_payload = _secure_read_bytes(config.control)
        else:
            controller_payload = _secure_read_bytes(config.control_source)
            if hashlib.sha256(controller_payload).hexdigest() != config.control_sha256:
                raise ContractViolation("runtime source controller SHA-256 changed")
            write_bytes_exclusive(config.control, controller_payload)
            controller_payload = _secure_read_bytes(config.control)
        if hashlib.sha256(controller_payload).hexdigest() != config.control_sha256:
            raise ContractViolation("persistent controller SHA-256 changed")

        if _path_occupied(config.credentials):
            credential_payload = _secure_read_bytes(config.credentials)
        else:
            credential_payload = _secure_read_bytes(config.credentials_source)
            _validate_credentials_payload(credential_payload)
            write_bytes_exclusive(config.credentials, credential_payload)
            credential_payload = _secure_read_bytes(config.credentials)
        _validate_credentials_payload(credential_payload)
        return {
            "status": "persistent_runtime_installed_and_verified",
            "controller_path": str(config.control),
            "controller_sha256": config.control_sha256,
            "credentials_path": str(config.credentials),
            "credentials_mode": "0600",
            "credentials_hash_persisted": False,
            "network_used": False,
            "mutating_api_call_made": False,
        }


def load_runtime(config: Config) -> Tuple[Any, Any]:
    """Execute only the exact controller bytes that passed the SHA-256 check."""

    source = _secure_read_bytes(config.control)
    actual_sha256 = hashlib.sha256(source).hexdigest()
    if actual_sha256 != config.control_sha256:
        raise ContractViolation("API base controller SHA-256 changed")
    module = types.ModuleType("zhixing_pinned_api_base_79e0070e")
    module.__file__ = str(config.control)
    module.__package__ = None
    try:
        code = compile(source, str(config.control), "exec", dont_inherit=True)
        exec(code, module.__dict__)
    except (Exception, SystemExit) as error:
        raise ContractViolation("pinned API base controller could not load") from error
    if module.__dict__.get("BASE_URL") != EXPECTED_API_BASE_URL:
        raise ContractViolation("pinned API base URL changed")
    endpoints = module.__dict__.get("ENDPOINTS")
    if not isinstance(endpoints, dict):
        raise ContractViolation("pinned API base has no endpoint map")
    if endpoints.get("list") != LIST_ENDPOINT:
        raise ContractViolation("pinned API list endpoint changed")
    if endpoints.get("detail") != DETAIL_ENDPOINT:
        raise ContractViolation("pinned API detail endpoint changed")
    if not callable(module.__dict__.get("call_api")):
        raise ContractViolation("pinned API base has no call_api function")
    if not callable(module.__dict__.get("load_credentials")):
        raise ContractViolation("pinned API base has no load_credentials function")
    credential_payload = _secure_read_bytes(config.credentials)
    credentials = _validate_credentials_payload(credential_payload)

    # The pinned controller's call_api performs a global lookup of
    # load_credentials on every request.  Bind that name to the tuple parsed
    # from the exact bytes securely read above, so replacing the credential
    # path between final preflight and renew cannot switch API accounts.
    def load_verified_credentials(
        verified: Tuple[str, str] = credentials,
    ) -> Tuple[str, str]:
        return verified

    module.CREDENTIALS = config.credentials
    module.load_credentials = load_verified_credentials
    return module, module


def _validated_instance_name(value: Any) -> str:
    if not isinstance(value, str) or INSTANCE_NAME_PATTERN.fullmatch(value) is None:
        raise ContractViolation("live instance name violates the strict whitelist")
    return value


def _safe_live_summary(item: Dict[str, Any]) -> Dict[str, Any]:
    """Whitelist live fields; never persist the raw instance response."""

    return {
        "instance_id": _exact_int(item.get("Id"), "live Id"),
        "status": _exact_int(item.get("Status"), "live Status"),
        "gpu_type": item.get("Gpu_type"),
        "gpu_num": _exact_int(item.get("Gpu_num"), "live Gpu_num"),
        "add_disk_GB": _exact_int(item.get("AddDisk"), "live AddDisk"),
        "adding_disk_GB": _exact_int(
            item.get("AddingDisk"), "live AddingDisk"
        ),
        "can_update_disk": item.get("CanUpdateDisk"),
        "due_time": _exact_int(item.get("Due_time"), "live Due_time"),
    }


def _current_instance(
    config: Config, controller_module: Any, base: Any
) -> Tuple[str, Dict[str, Any], Dict[str, Any]]:
    if controller_module is not base:
        raise ContractViolation("runtime controller is not the pinned API base")
    response = base.call_api(base.ENDPOINTS["list"])
    try:
        items = response["data"]["list"]
    except (KeyError, TypeError) as error:
        raise ContractViolation("unrecognized instance-list response") from error
    if not isinstance(items, list):
        raise ContractViolation("instance-list data is not a list")
    matches = [
        item for item in items
        if isinstance(item, dict)
        and item.get("Id") == config.instance_id
        and not isinstance(item.get("Id"), bool)
    ]
    if len(matches) != 1:
        raise ContractViolation("instance 641496 was not found exactly once")
    list_item = matches[0]
    name = _validated_instance_name(list_item.get("Container_name"))
    detail_response = base.call_api(
        base.ENDPOINTS["detail"], {"instance_name": name}
    )
    try:
        detail_item = detail_response["data"]
    except (KeyError, TypeError) as error:
        raise ContractViolation("unrecognized instance-detail response") from error
    if not isinstance(detail_item, dict):
        raise ContractViolation("instance-detail data is not an object")
    detail_name = _validated_instance_name(detail_item.get("Container_name"))
    if detail_name != name:
        raise ContractViolation("list/detail instance names differ")
    list_safe = _safe_live_summary(list_item)
    detail_safe = _safe_live_summary(detail_item)
    if list_safe != detail_safe:
        raise ContractViolation("list/detail safe instance views differ")
    safe = list_safe
    exact_resources = {
        "instance_id": config.instance_id,
        "status": config.status,
        "gpu_type": config.gpu_type,
        "gpu_num": config.gpu_num,
        "add_disk_GB": config.add_disk_GB,
        "adding_disk_GB": config.adding_disk_GB,
        "can_update_disk": config.can_update_disk,
    }
    for field, wanted in exact_resources.items():
        changed = (
            safe[field] is not wanted
            if isinstance(wanted, bool)
            else safe[field] != wanted
        )
        if changed:
            raise ContractViolation(
                f"live instance field {field} is {safe[field]!r}, expected {wanted!r}"
            )
    return name, list_item, safe


def _require_due(safe: Dict[str, Any], expected: int, phase: str) -> None:
    if safe["due_time"] != expected:
        raise ContractViolation(
            f"{phase} due_time is {safe['due_time']}, expected {expected}"
        )


def _checked_quote_components(sheet: Dict[str, Any], config: Config) -> List[Dict[str, Any]]:
    components = sheet.get("InstanceSubPriceSheet")
    if not isinstance(components, list) or len(components) != 3:
        raise ContractViolation("renewal quote does not contain exactly three components")
    by_resource: Dict[str, Dict[str, Any]] = {}
    for component in components:
        if not isinstance(component, dict):
            raise ContractViolation("renewal quote component is not an object")
        resource_type = component.get("ResourceType")
        if resource_type not in EXPECTED_QUOTE_COMPONENTS:
            raise ContractViolation("renewal quote has an unexpected resource component")
        if resource_type in by_resource:
            raise ContractViolation("renewal quote duplicates a resource component")
        by_resource[resource_type] = component
    safe_components: List[Dict[str, Any]] = []
    component_sum = Decimal("0")
    for resource_type in ("gpu", "spec", "disk"):
        component = by_resource[resource_type]
        expected = EXPECTED_QUOTE_COMPONENTS[resource_type]
        if component.get("Type") != expected["type"]:
            raise ContractViolation(f"renewal {resource_type} quote type changed")
        if component.get("Title") != expected["title"]:
            raise ContractViolation(f"renewal {resource_type} quote title changed")
        if _exact_int(
            component.get("Amount"), f"renewal {resource_type} quote Amount"
        ) != expected["amount"]:
            raise ContractViolation(f"renewal {resource_type} quote amount changed")
        if _exact_int(
            component.get("HourLen"), f"renewal {resource_type} quote HourLen"
        ) != config.pack_hours:
            raise ContractViolation(f"renewal {resource_type} quote hours changed")
        if component.get("PayTypeFirst") != "power":
            raise ContractViolation(f"renewal {resource_type} payment type changed")
        unit_price = _finite_decimal(
            component.get("SingleHourPrice"),
            f"renewal {resource_type} SingleHourPrice",
        )
        final_discount = _finite_decimal(
            component.get("FinalDiscount"),
            f"renewal {resource_type} FinalDiscount",
        )
        final_cost = _finite_decimal(
            component.get("FinalCost"), f"renewal {resource_type} FinalCost"
        )
        if unit_price != expected["single_hour_price_CNY"]:
            raise ContractViolation(f"renewal {resource_type} unit price changed")
        if final_discount != expected["final_discount"]:
            raise ContractViolation(f"renewal {resource_type} discount changed")
        if final_cost != expected["final_cost_CNY"]:
            raise ContractViolation(f"renewal {resource_type} cost changed")
        component_sum += final_cost
        safe_components.append({
            "resource_type": resource_type,
            "type": expected["type"],
            "contract": expected["contract"],
            "amount": expected["amount"],
            "hour_len": config.pack_hours,
            "single_hour_price_CNY": format(unit_price, "f"),
            "final_discount": format(final_discount, "f"),
            "final_cost_CNY": format(final_cost, "f"),
        })
    if component_sum != EXPECTED_FINAL_COST_CNY:
        raise ContractViolation("renewal quote component sum changed")
    return safe_components


def _checked_quote(
    config: Config, controller_module: Any, base: Any
) -> Tuple[str, Dict[str, Any], Dict[str, Any]]:
    name, _, before = _current_instance(config, controller_module, base)
    _require_due(before, config.old_due_time, "pre-renewal")
    payload = {
        "instance_name": name,
        "unit": "hour",
        "pack_time": config.pack_hours,
        "pay_type_first": "power",
    }
    response = base.call_api(QUOTE_ENDPOINT, payload)
    if not isinstance(response, dict):
        raise ContractViolation("renewal quote response is not an object")
    if response.get("success") is not True or str(response.get("code")) != "2000":
        raise ContractViolation("renewal quote response is not successful")
    try:
        data = response["data"]
        sheet = data["PackPriceInfoSheet"]
        item = data["ItemInfoSheet"]
    except (KeyError, TypeError) as error:
        raise ContractViolation("unrecognized renewal quote response") from error
    if not isinstance(sheet, dict) or not isinstance(item, dict):
        raise ContractViolation("renewal quote sheets are not objects")
    if sheet.get("ResourceType") != "instance":
        raise ContractViolation("quote ResourceType is not instance")
    if sheet.get("Type") != "renew_instance":
        raise ContractViolation("quote Type is not renew_instance")
    if _exact_int(sheet.get("HourLen"), "quote HourLen") != config.pack_hours:
        raise ContractViolation("quote HourLen differs from 120 hours")
    if sheet.get("PayTypeFirst") != "power":
        raise ContractViolation("quote payment type changed")
    if _exact_int(item.get("Due_time"), "quote Due_time") != config.old_due_time:
        raise ContractViolation("quote starts from the wrong due time")
    if _exact_int(item.get("NewDuetime"), "quote NewDuetime") != config.new_due_time:
        raise ContractViolation("quote does not end at the exact target due time")
    if item.get("PayTypeFirst") != "power":
        raise ContractViolation("quote item payment type changed")
    safe_components = _checked_quote_components(sheet, config)
    cost = _finite_decimal(sheet.get("FinalCost"), "quote FinalCost")
    if not Decimal("0") < cost <= config.max_cost_CNY:
        raise ContractViolation(
            f"renewal quote {cost} exceeds fail-closed ceiling {config.max_cost_CNY}"
        )
    if cost != EXPECTED_FINAL_COST_CNY:
        raise ContractViolation(
            f"renewal quote total {cost} differs from {EXPECTED_FINAL_COST_CNY}"
        )
    if sum(
        _finite_decimal(component["final_cost_CNY"], "safe quote component")
        for component in safe_components
    ) != cost:
        raise ContractViolation("renewal quote total differs from component sum")
    quote = {
        "resource_type": "instance",
        "operation_type": "renew_instance",
        "hour_len": config.pack_hours,
        "old_due_time": config.old_due_time,
        "new_due_time": config.new_due_time,
        "pay_type_first": "power",
        "final_cost_CNY": format(cost, "f"),
        "components": safe_components,
        "response_success": True,
        "response_code_is_2000": True,
    }
    confirmed_name, _, confirmed = _current_instance(
        config, controller_module, base
    )
    if confirmed_name != name:
        raise ContractViolation("instance identity changed after renewal quote")
    _require_due(confirmed, config.old_due_time, "post-quote pre-intent")
    if confirmed != before:
        raise ContractViolation("instance resource contract changed after renewal quote")
    return name, confirmed, quote


def _resource_contract(config: Config) -> Dict[str, Any]:
    return {
        "status": config.status,
        "gpu_type": config.gpu_type,
        "gpu_num": config.gpu_num,
        "add_disk_GB": config.add_disk_GB,
        "adding_disk_GB": config.adding_disk_GB,
        "can_update_disk": config.can_update_disk,
    }


def _runtime_evidence(config: Config) -> Dict[str, str]:
    return {
        "controller_path": str(config.control),
        "controller_sha256": config.control_sha256,
        "api_base_path": str(config.control),
        "api_base_sha256": config.control_sha256,
        "runtime_load_mode": (
            "secure-read-sha256-compile-exec-same-bytes;"
            "credentials-bound-in-memory"
        ),
        "credentials_path": str(config.credentials),
        "credentials_policy": (
            "persistent mode-0600; verified bytes bound in memory; "
            "no credential hash in evidence"
        ),
    }


def _make_intent(
    config: Config,
    before: Dict[str, Any],
    quote: Dict[str, Any],
    prior: Sequence[Dict[str, Any]],
    prior_anchor_payload: bytes,
    timers: Timers,
) -> Dict[str, Any]:
    return {
        "schema": "zhixing-renew-once-intent-v2",
        "status": "armed_before_non_idempotent_renew_call",
        "instance_id": config.instance_id,
        "old_due_time": config.old_due_time,
        "new_due_time": config.new_due_time,
        "extension_seconds": config.pack_hours * 3600,
        "pack_hours": config.pack_hours,
        "resource_contract": _resource_contract(config),
        "before": before,
        "quote": quote,
        "prior_renewal_chain": list(prior),
        "prior_anchor_file": config.prior_anchor.name,
        "prior_anchor_sha256": hashlib.sha256(prior_anchor_payload).hexdigest(),
        **_runtime_evidence(config),
        "shared_instance_mutation_lock": str(config.mutation_lock),
        "call_marker_file": config.call_marker.name,
        "retry_policy": "never replay renew; reconcile with list+detail only",
        "evidence_policy": "whitelisted fields only; no raw API responses",
        "armed_unix": timers.wall_time(),
    }


def _validate_exact_keys(value: Dict[str, Any], expected: Sequence[str], label: str) -> None:
    if not isinstance(value, dict):
        raise ContractViolation(f"{label} is not an object")
    if set(value) != set(expected):
        raise ContractViolation(f"{label} has non-whitelisted or missing fields")


def _validate_intent_resource_contract(config: Config, value: Any) -> None:
    _validate_exact_keys(value, (
        "status",
        "gpu_type",
        "gpu_num",
        "add_disk_GB",
        "adding_disk_GB",
        "can_update_disk",
    ), "renewal intent resource contract")
    if value.get("gpu_type") != config.gpu_type:
        raise ContractViolation("renewal intent GPU type changed")
    if value.get("can_update_disk") is not config.can_update_disk:
        raise ContractViolation("renewal intent CanUpdateDisk changed")
    expected_integers = {
        "status": config.status,
        "gpu_num": config.gpu_num,
        "add_disk_GB": config.add_disk_GB,
        "adding_disk_GB": config.adding_disk_GB,
    }
    for field, expected in expected_integers.items():
        if _exact_int(
            value.get(field), f"renewal intent resource {field}"
        ) != expected:
            raise ContractViolation(f"renewal intent resource field changed: {field}")


def _validate_safe_snapshot(
    config: Config, value: Any, expected_due_time: int, label: str
) -> None:
    _validate_exact_keys(value, (
        "instance_id",
        "status",
        "gpu_type",
        "gpu_num",
        "add_disk_GB",
        "adding_disk_GB",
        "can_update_disk",
        "due_time",
    ), label)
    if value.get("gpu_type") != config.gpu_type:
        raise ContractViolation(f"{label} GPU type changed")
    if value.get("can_update_disk") is not config.can_update_disk:
        raise ContractViolation(f"{label} CanUpdateDisk changed")
    expected_integers = {
        "instance_id": config.instance_id,
        "status": config.status,
        "gpu_num": config.gpu_num,
        "add_disk_GB": config.add_disk_GB,
        "adding_disk_GB": config.adding_disk_GB,
        "due_time": expected_due_time,
    }
    for field, expected in expected_integers.items():
        if _exact_int(
            value.get(field), f"{label} {field}"
        ) != expected:
            raise ContractViolation(f"{label} field changed: {field}")


def _validate_intent_before(config: Config, value: Any) -> None:
    _validate_safe_snapshot(
        config, value, config.old_due_time, "renewal intent before-summary"
    )


def _validate_intent_prior_chain(
    value: Any, expected: Sequence[Dict[str, Any]]
) -> None:
    if not isinstance(value, list) or len(value) != len(expected):
        raise ContractViolation("renewal intent prior evidence chain changed")
    integer_fields = (
        "old_due_time",
        "new_due_time",
        "extension_seconds",
        "pack_hours",
    )
    for index, (actual, wanted) in enumerate(zip(value, expected)):
        _validate_exact_keys(actual, (
            "evidence_file",
            "evidence_sha256",
            *integer_fields,
        ), f"renewal intent prior evidence {index}")
        for field in ("evidence_file", "evidence_sha256"):
            if not isinstance(actual.get(field), str) or actual[field] != wanted[field]:
                raise ContractViolation(
                    f"renewal intent prior evidence field changed: {index}:{field}"
                )
        for field in integer_fields:
            if _exact_int(
                actual.get(field), f"renewal intent prior {index}:{field}"
            ) != wanted[field]:
                raise ContractViolation(
                    f"renewal intent prior evidence field changed: {index}:{field}"
                )


def _validate_safe_quote(config: Config, quote: Any) -> None:
    _validate_exact_keys(quote, (
        "resource_type",
        "operation_type",
        "hour_len",
        "old_due_time",
        "new_due_time",
        "pay_type_first",
        "final_cost_CNY",
        "components",
        "response_success",
        "response_code_is_2000",
    ), "renewal intent quote")
    expected_integers = {
        "hour_len": config.pack_hours,
        "old_due_time": config.old_due_time,
        "new_due_time": config.new_due_time,
    }
    for field, expected in expected_integers.items():
        if _exact_int(quote.get(field), f"renewal intent quote {field}") != expected:
            raise ContractViolation(f"renewal intent quote field changed: {field}")
    if (
        quote.get("resource_type") != "instance"
        or quote.get("operation_type") != "renew_instance"
        or quote.get("pay_type_first") != "power"
        or quote.get("response_success") is not True
        or quote.get("response_code_is_2000") is not True
    ):
        raise ContractViolation("renewal intent quote contract changed")
    if not isinstance(quote.get("final_cost_CNY"), str):
        raise ContractViolation("renewal intent quote cost is not canonical text")
    cost = _finite_decimal(quote["final_cost_CNY"], "intent quote cost")
    if cost != EXPECTED_FINAL_COST_CNY or cost > config.max_cost_CNY:
        raise ContractViolation("renewal intent quote cost changed")
    components = quote.get("components")
    if not isinstance(components, list) or len(components) != 3:
        raise ContractViolation("renewal intent quote components changed")
    component_sum = Decimal("0")
    for index, resource_type in enumerate(("gpu", "spec", "disk")):
        component = components[index]
        expected = EXPECTED_QUOTE_COMPONENTS[resource_type]
        _validate_exact_keys(component, (
            "resource_type",
            "type",
            "contract",
            "amount",
            "hour_len",
            "single_hour_price_CNY",
            "final_discount",
            "final_cost_CNY",
        ), f"renewal intent quote component {index}")
        if (
            component.get("resource_type") != resource_type
            or component.get("type") != expected["type"]
            or component.get("contract") != expected["contract"]
        ):
            raise ContractViolation("renewal intent quote component identity changed")
        if _exact_int(
            component.get("amount"), f"intent quote component {index} amount"
        ) != expected["amount"]:
            raise ContractViolation("renewal intent quote component amount changed")
        if _exact_int(
            component.get("hour_len"), f"intent quote component {index} hour_len"
        ) != config.pack_hours:
            raise ContractViolation("renewal intent quote component hours changed")
        decimal_fields = {
            "single_hour_price_CNY": expected["single_hour_price_CNY"],
            "final_discount": expected["final_discount"],
            "final_cost_CNY": expected["final_cost_CNY"],
        }
        for field, wanted in decimal_fields.items():
            if not isinstance(component.get(field), str):
                raise ContractViolation("renewal quote component is not canonical text")
            if _finite_decimal(component[field], field) != wanted:
                raise ContractViolation(f"renewal quote component changed: {field}")
        component_sum += expected["final_cost_CNY"]
    if component_sum != cost:
        raise ContractViolation("renewal intent quote component sum changed")


def _validate_existing_intent(
    config: Config,
    prior: Sequence[Dict[str, Any]],
    prior_anchor_payload: bytes,
) -> Tuple[Dict[str, Any], bytes]:
    value, payload = _secure_read_json(config.intent)
    _validate_exact_keys(value, (
        "schema",
        "status",
        "instance_id",
        "old_due_time",
        "new_due_time",
        "extension_seconds",
        "pack_hours",
        "resource_contract",
        "before",
        "quote",
        "prior_renewal_chain",
        "prior_anchor_file",
        "prior_anchor_sha256",
        "controller_path",
        "controller_sha256",
        "api_base_path",
        "api_base_sha256",
        "runtime_load_mode",
        "credentials_path",
        "credentials_policy",
        "shared_instance_mutation_lock",
        "call_marker_file",
        "retry_policy",
        "evidence_policy",
        "armed_unix",
    ), "renewal intent")
    exact_scalars = {
        "schema": "zhixing-renew-once-intent-v2",
        "status": "armed_before_non_idempotent_renew_call",
        "instance_id": config.instance_id,
        "old_due_time": config.old_due_time,
        "new_due_time": config.new_due_time,
        "extension_seconds": config.pack_hours * 3600,
        "pack_hours": config.pack_hours,
        "prior_anchor_file": config.prior_anchor.name,
        "prior_anchor_sha256": hashlib.sha256(prior_anchor_payload).hexdigest(),
        **_runtime_evidence(config),
        "shared_instance_mutation_lock": str(config.mutation_lock),
        "call_marker_file": config.call_marker.name,
        "retry_policy": "never replay renew; reconcile with list+detail only",
        "evidence_policy": "whitelisted fields only; no raw API responses",
    }
    for field, expected in exact_scalars.items():
        if value.get(field) != expected or (
            isinstance(expected, int) and isinstance(value.get(field), bool)
        ):
            raise ContractViolation(f"renewal intent field changed: {field}")
    _validate_intent_resource_contract(config, value.get("resource_contract"))
    _validate_intent_before(config, value.get("before"))
    _validate_safe_quote(config, value.get("quote"))
    _validate_intent_prior_chain(value.get("prior_renewal_chain"), prior)
    armed = value.get("armed_unix")
    if isinstance(armed, bool) or not isinstance(armed, (int, float)):
        raise ContractViolation("renewal intent armed_unix is invalid")
    if not math.isfinite(float(armed)):
        raise ContractViolation("renewal intent armed_unix is not finite")
    return value, payload


def _make_call_marker(
    config: Config,
    intent_payload: bytes,
    final_preflight: Dict[str, Any],
    timers: Timers,
) -> Dict[str, Any]:
    return {
        "schema": "zhixing-renew-call-marker-v1",
        "status": "fsynced_immediately_before_sole_renew_call",
        "instance_id": config.instance_id,
        "intent_sha256": hashlib.sha256(intent_payload).hexdigest(),
        "old_due_time": config.old_due_time,
        "new_due_time": config.new_due_time,
        "pack_hours": config.pack_hours,
        "final_preflight": final_preflight,
        **_runtime_evidence(config),
        "shared_instance_mutation_lock": str(config.mutation_lock),
        "endpoint": RENEW_ENDPOINT,
        "request_summary": {
            "unit": "hour",
            "pack_time": config.pack_hours,
            "pay_type_first": "power",
        },
        "replay_rule": "marker presence means call may have happened; never replay",
        "created_unix": timers.wall_time(),
    }


def _validate_call_marker(
    config: Config, intent_payload: bytes
) -> Tuple[Dict[str, Any], bytes]:
    marker, payload = _secure_read_json(config.call_marker)
    _validate_exact_keys(marker, (
        "schema",
        "status",
        "instance_id",
        "intent_sha256",
        "old_due_time",
        "new_due_time",
        "pack_hours",
        "final_preflight",
        "controller_path",
        "controller_sha256",
        "api_base_path",
        "api_base_sha256",
        "runtime_load_mode",
        "credentials_path",
        "credentials_policy",
        "shared_instance_mutation_lock",
        "endpoint",
        "request_summary",
        "replay_rule",
        "created_unix",
    ), "renewal call marker")
    expected_scalars = {
        "schema": "zhixing-renew-call-marker-v1",
        "status": "fsynced_immediately_before_sole_renew_call",
        "instance_id": config.instance_id,
        "intent_sha256": hashlib.sha256(intent_payload).hexdigest(),
        "old_due_time": config.old_due_time,
        "new_due_time": config.new_due_time,
        "pack_hours": config.pack_hours,
        **_runtime_evidence(config),
        "shared_instance_mutation_lock": str(config.mutation_lock),
        "endpoint": RENEW_ENDPOINT,
        "replay_rule": "marker presence means call may have happened; never replay",
    }
    for field, expected in expected_scalars.items():
        if marker.get(field) != expected or (
            isinstance(expected, int) and isinstance(marker.get(field), bool)
        ):
            raise ContractViolation(f"renewal call marker field changed: {field}")
    _validate_intent_before(config, marker.get("final_preflight"))
    request = marker.get("request_summary")
    _validate_exact_keys(
        request, ("unit", "pack_time", "pay_type_first"), "call request summary"
    )
    if (
        request.get("unit") != "hour"
        or _exact_int(request.get("pack_time"), "call marker pack_time")
        != config.pack_hours
        or request.get("pay_type_first") != "power"
    ):
        raise ContractViolation("renewal call marker request changed")
    created = marker.get("created_unix")
    if isinstance(created, bool) or not isinstance(created, (int, float)):
        raise ContractViolation("renewal call marker created_unix is invalid")
    if not math.isfinite(float(created)):
        raise ContractViolation("renewal call marker created_unix is not finite")
    return marker, payload


def _safe_response_summary(response: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Persist booleans derived from an API response, never response contents."""

    if response is None:
        return None
    if not isinstance(response, dict):
        return {
            "returned": True,
            "response_was_object": False,
            "success_is_true": False,
            "code_is_2000": False,
        }
    return {
        "returned": True,
        "response_was_object": True,
        "success_is_true": response.get("success") is True,
        "code_is_2000": str(response.get("code")) == "2000",
    }


def _reconcile_due(
    config: Config,
    controller_module: Any,
    base: Any,
    timers: Timers,
) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """Use only the read-only list endpoint to prove the exact new due time."""

    deadline = timers.monotonic() + config.observe_seconds
    last_error_class: Optional[str] = None
    while True:
        try:
            _, _, safe = _current_instance(config, controller_module, base)
            last_error_class = None
            if safe["due_time"] == config.new_due_time:
                return safe, None
            if safe["due_time"] != config.old_due_time:
                raise ContractViolation(
                    "observed due time is neither the exact old nor exact target time"
                )
        except ContractViolation:
            raise
        except (Exception, SystemExit) as error:
            last_error_class = type(error).__name__
        remaining = deadline - timers.monotonic()
        if remaining <= 0:
            return None, last_error_class
        timers.sleep(min(config.observe_interval_seconds, remaining))


def _make_result(
    config: Config,
    intent_payload: bytes,
    call_marker_payload: bytes,
    intent: Dict[str, Any],
    observed: Dict[str, Any],
    outcome: str,
    mutation_called_this_invocation: bool,
    mutation_response: Optional[Dict[str, Any]],
    mutation_error_class: Optional[str],
    observe_error_class: Optional[str],
    timers: Timers,
) -> Dict[str, Any]:
    if observed["due_time"] - config.old_due_time != config.pack_hours * 3600:
        raise ContractViolation("observed due time does not prove exactly +120 hours")
    return {
        "schema": "zhixing-renew-once-result-v2",
        "status": "observed_exact_120_hour_extension",
        "outcome": outcome,
        "instance_id": config.instance_id,
        "old_due_time": config.old_due_time,
        "new_due_time": config.new_due_time,
        "extension_seconds": config.pack_hours * 3600,
        "pack_hours": config.pack_hours,
        "resource_contract": _resource_contract(config),
        "quoted_cost_CNY": intent["quote"]["final_cost_CNY"],
        "prior_anchor_file": intent["prior_anchor_file"],
        "prior_anchor_sha256": intent["prior_anchor_sha256"],
        "intent_sha256": hashlib.sha256(intent_payload).hexdigest(),
        "call_marker_file": config.call_marker.name,
        "call_marker_sha256": hashlib.sha256(call_marker_payload).hexdigest(),
        **_runtime_evidence(config),
        "shared_instance_mutation_lock": str(config.mutation_lock),
        "observed": observed,
        "proof": "read-only list+detail observed exact old_due + 120h",
        "mutation_called_this_invocation": mutation_called_this_invocation,
        "mutation_response": _safe_response_summary(mutation_response),
        "mutation_error_class": mutation_error_class,
        "observe_error_class": observe_error_class,
        "renew_must_not_be_replayed": True,
        "completed_unix": timers.wall_time(),
    }


def _validate_optional_error_class(value: Any, label: str) -> None:
    if value is None:
        return
    if not isinstance(value, str) or re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]{0,127}", value) is None:
        raise ContractViolation(f"{label} is not a safe exception class")


def _validate_response_summary(value: Any) -> None:
    if value is None:
        return
    _validate_exact_keys(value, (
        "returned",
        "response_was_object",
        "success_is_true",
        "code_is_2000",
    ), "renewal response summary")
    for field in value:
        if value[field] not in (True, False) or not isinstance(value[field], bool):
            raise ContractViolation("renewal response summary is not boolean-only")


def _validate_result(
    config: Config,
    result: Dict[str, Any],
    intent_payload: bytes,
    call_marker_payload: bytes,
) -> None:
    try:
        validated_intent = json.loads(intent_payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractViolation("validated intent bytes became invalid") from error
    if not isinstance(validated_intent, dict):
        raise ContractViolation("validated intent root changed")
    _validate_exact_keys(result, (
        "schema",
        "status",
        "outcome",
        "instance_id",
        "old_due_time",
        "new_due_time",
        "extension_seconds",
        "pack_hours",
        "resource_contract",
        "quoted_cost_CNY",
        "prior_anchor_file",
        "prior_anchor_sha256",
        "intent_sha256",
        "call_marker_file",
        "call_marker_sha256",
        "controller_path",
        "controller_sha256",
        "api_base_path",
        "api_base_sha256",
        "runtime_load_mode",
        "credentials_path",
        "credentials_policy",
        "shared_instance_mutation_lock",
        "observed",
        "proof",
        "mutation_called_this_invocation",
        "mutation_response",
        "mutation_error_class",
        "observe_error_class",
        "renew_must_not_be_replayed",
        "completed_unix",
    ), "renewal result")
    expected_scalars = {
        "schema": "zhixing-renew-once-result-v2",
        "status": "observed_exact_120_hour_extension",
        "instance_id": config.instance_id,
        "old_due_time": config.old_due_time,
        "new_due_time": config.new_due_time,
        "extension_seconds": config.pack_hours * 3600,
        "pack_hours": config.pack_hours,
        "quoted_cost_CNY": format(EXPECTED_FINAL_COST_CNY, "f"),
        "prior_anchor_file": config.prior_anchor.name,
        "prior_anchor_sha256": validated_intent.get("prior_anchor_sha256"),
        "intent_sha256": hashlib.sha256(intent_payload).hexdigest(),
        "call_marker_file": config.call_marker.name,
        "call_marker_sha256": hashlib.sha256(call_marker_payload).hexdigest(),
        **_runtime_evidence(config),
        "shared_instance_mutation_lock": str(config.mutation_lock),
        "proof": "read-only list+detail observed exact old_due + 120h",
    }
    for field, expected in expected_scalars.items():
        if result.get(field) != expected or (
            isinstance(expected, int) and isinstance(result.get(field), bool)
        ):
            raise ContractViolation(f"renewal result field changed: {field}")
    if result.get("outcome") not in {
        "applied",
        "applied_after_ambiguous_response",
        "applied_reconciled_existing_intent",
    }:
        raise ContractViolation("renewal result outcome is invalid")
    _validate_intent_resource_contract(config, result.get("resource_contract"))
    _validate_safe_snapshot(
        config,
        result.get("observed"),
        config.new_due_time,
        "renewal result observed-summary",
    )
    if not isinstance(result.get("mutation_called_this_invocation"), bool):
        raise ContractViolation("renewal result mutation-called flag is invalid")
    if result.get("renew_must_not_be_replayed") is not True:
        raise ContractViolation("renewal result replay prohibition changed")
    _validate_response_summary(result.get("mutation_response"))
    _validate_optional_error_class(
        result.get("mutation_error_class"), "mutation_error_class"
    )
    _validate_optional_error_class(
        result.get("observe_error_class"), "observe_error_class"
    )
    completed = result.get("completed_unix")
    if isinstance(completed, bool) or not isinstance(completed, (int, float)):
        raise ContractViolation("renewal result completed_unix is invalid")
    if not math.isfinite(float(completed)):
        raise ContractViolation("renewal result completed_unix is not finite")


def validate_persistent_result(
    intent_path: Path,
    result_path: Path,
    mutation_lock: Path,
) -> Tuple[Dict[str, Any], str]:
    """Pure-local validator exported for the storage generation-zero gate."""

    intent_path = Path(intent_path)
    result_path = Path(result_path)
    mutation_lock = Path(mutation_lock)
    if intent_path.name != INTENT.name or result_path.name != RESULT.name:
        raise ContractViolation("renewal evidence basenames are not canonical")
    call_marker = intent_path.parent / CALL_MARKER.name
    config = replace(
        DEFAULT_CONFIG,
        control=intent_path.parent / API_BASE.name,
        credentials=intent_path.parent / API_CREDENTIALS.name,
        intent=intent_path,
        result=result_path,
        call_marker=call_marker,
        prior_anchor=intent_path.parent / PRIOR_ANCHOR.name,
        mutation_lock=mutation_lock,
        require_persistent_xfs=(intent_path.parent == PRIVATE_STATE_DIR),
    )
    validate_transaction_paths(config)
    prior, prior_anchor_payload = _validate_prior_anchor(config)
    _, intent_payload = _validate_existing_intent(
        config, prior, prior_anchor_payload
    )
    _, call_marker_payload = _validate_call_marker(config, intent_payload)
    result, result_payload = _secure_read_json(config.result)
    _validate_result(config, result, intent_payload, call_marker_payload)
    return result, hashlib.sha256(result_payload).hexdigest()


def _uncertain_summary(
    config: Config,
    mutation_called: bool,
    mutation_error_class: Optional[str],
    observe_error_class: Optional[str],
    call_marker_exists: bool,
) -> Dict[str, Any]:
    return {
        "status": "intent_exists_but_exact_extension_not_yet_proven",
        "instance_id": config.instance_id,
        "old_due_time": config.old_due_time,
        "target_due_time": config.new_due_time,
        "mutation_called_this_invocation": mutation_called,
        "mutation_error_class": mutation_error_class,
        "observe_error_class": observe_error_class,
        "call_marker_exists": call_marker_exists,
        "result_written": False,
        "next_action": "run reconcile; never call renew again",
    }


def _write_success_result(
    config: Config,
    result: Dict[str, Any],
    intent_payload: bytes,
    call_marker_payload: bytes,
) -> None:
    if _path_occupied(config.result):
        raise ContractViolation("renewal result path is already occupied")
    write_private_exclusive(config.result, result)
    persisted, _ = _secure_read_json(config.result)
    _validate_result(config, persisted, intent_payload, call_marker_payload)


def quote_only(
    config: Config = DEFAULT_CONFIG,
    runtime_factory: Callable[[Config], Tuple[Any, Any]] = load_runtime,
) -> Dict[str, Any]:
    validate_transaction_paths(config)
    with acquire_instance_mutation_lock(config.mutation_lock):
        if (
            _path_occupied(config.intent)
            or _path_occupied(config.call_marker)
            or _path_occupied(config.result)
        ):
            raise ContractViolation("one-shot renewal transaction is already armed")
        prior, prior_anchor_payload = load_or_create_prior_anchor(
            config, Timers()
        )
        controller_module, base = runtime_factory(config)
        _, before, quote = _checked_quote(config, controller_module, base)
        return {
            "status": "quote_verified",
            "instance": before,
            "quote": quote,
            "prior_renewal_chain": prior,
            "prior_anchor_file": config.prior_anchor.name,
            "prior_anchor_sha256": hashlib.sha256(
                prior_anchor_payload
            ).hexdigest(),
            **_runtime_evidence(config),
            "shared_instance_mutation_lock": str(config.mutation_lock),
            "mutating_call_made": False,
        }


def reconcile_only(
    config: Config = DEFAULT_CONFIG,
    runtime_factory: Callable[[Config], Tuple[Any, Any]] = load_runtime,
    timers: Timers = Timers(),
) -> Tuple[int, Dict[str, Any]]:
    validate_transaction_paths(config)
    with acquire_instance_mutation_lock(config.mutation_lock):
        if _path_occupied(config.result):
            raise ContractViolation("renewal result already exists; transaction is complete")
        if not _path_occupied(config.intent):
            raise ContractViolation("renewal intent does not exist; nothing to reconcile")
        prior, prior_anchor_payload = _validate_prior_anchor(config)
        intent, intent_payload = _validate_existing_intent(
            config, prior, prior_anchor_payload
        )
        marker_payload: Optional[bytes] = None
        if _path_occupied(config.call_marker):
            _, marker_payload = _validate_call_marker(config, intent_payload)
        controller_module, base = runtime_factory(config)
        observed, observe_error_class = _reconcile_due(
            config, controller_module, base, timers
        )
        if observed is None:
            return 2, _uncertain_summary(
                config,
                mutation_called=False,
                mutation_error_class=None,
                observe_error_class=observe_error_class,
                call_marker_exists=marker_payload is not None,
            )
        if marker_payload is None:
            return 3, {
                "status": "target_observed_without_durable_call_marker",
                "instance_id": config.instance_id,
                "target_due_time": config.new_due_time,
                "result_written": False,
                "next_action": "manual provenance review; never replay renew",
            }
        result = _make_result(
            config,
            intent_payload,
            marker_payload,
            intent,
            observed,
            outcome="applied_reconciled_existing_intent",
            mutation_called_this_invocation=False,
            mutation_response=None,
            mutation_error_class=None,
            observe_error_class=observe_error_class,
            timers=timers,
        )
        _write_success_result(config, result, intent_payload, marker_payload)
        return 0, result


def renew_once(
    config: Config = DEFAULT_CONFIG,
    runtime_factory: Callable[[Config], Tuple[Any, Any]] = load_runtime,
    timers: Timers = Timers(),
) -> Tuple[int, Dict[str, Any]]:
    validate_transaction_paths(config)
    with acquire_instance_mutation_lock(config.mutation_lock):
        if _path_occupied(config.result):
            raise ContractViolation("renewal result already exists; refusing replay")
        if _path_occupied(config.intent):
            prior, prior_anchor_payload = _validate_prior_anchor(config)
        else:
            prior, prior_anchor_payload = load_or_create_prior_anchor(
                config, timers
            )
        controller_module, base = runtime_factory(config)

        if _path_occupied(config.intent):
            # SECURITY INVARIANT: an existing intent permanently disables the
            # non-idempotent endpoint.  This branch is read-only reconciliation.
            intent, intent_payload = _validate_existing_intent(
                config, prior, prior_anchor_payload
            )
            marker_payload: Optional[bytes] = None
            if _path_occupied(config.call_marker):
                _, marker_payload = _validate_call_marker(config, intent_payload)
            observed, observe_error_class = _reconcile_due(
                config, controller_module, base, timers
            )
            if observed is None:
                return 2, _uncertain_summary(
                    config,
                    mutation_called=False,
                    mutation_error_class=None,
                    observe_error_class=observe_error_class,
                    call_marker_exists=marker_payload is not None,
                )
            if marker_payload is None:
                return 3, {
                    "status": "target_observed_without_durable_call_marker",
                    "instance_id": config.instance_id,
                    "target_due_time": config.new_due_time,
                    "result_written": False,
                    "next_action": "manual provenance review; never replay renew",
                }
            result = _make_result(
                config,
                intent_payload,
                marker_payload,
                intent,
                observed,
                outcome="applied_reconciled_existing_intent",
                mutation_called_this_invocation=False,
                mutation_response=None,
                mutation_error_class=None,
                observe_error_class=observe_error_class,
                timers=timers,
            )
            _write_success_result(config, result, intent_payload, marker_payload)
            return 0, result

        if _path_occupied(config.call_marker):
            raise ContractViolation("call marker exists without renewal intent")

        name, before, quote = _checked_quote(config, controller_module, base)
        intent = _make_intent(
            config, before, quote, prior, prior_anchor_payload, timers
        )
        write_private_exclusive(config.intent, intent)
        intent_on_disk, intent_payload = _validate_existing_intent(
            config, prior, prior_anchor_payload
        )

        final_name, _, final_preflight = _current_instance(
            config, controller_module, base
        )
        if final_name != name or final_preflight != before:
            raise ContractViolation("final post-intent instance preflight changed")
        _require_due(final_preflight, config.old_due_time, "post-intent pre-call")
        marker = _make_call_marker(
            config, intent_payload, final_preflight, timers
        )
        write_private_exclusive(config.call_marker, marker)
        _, marker_payload = _validate_call_marker(config, intent_payload)

        mutation_response: Optional[Dict[str, Any]] = None
        mutation_error_class: Optional[str] = None
        try:
            # SECURITY INVARIANT: this is the program's sole mutating API call.
            mutation_response = base.call_api(RENEW_ENDPOINT, {
                "instance_name": name,
                "unit": "hour",
                "pack_time": config.pack_hours,
                "pay_type_first": "power",
            })
        except (Exception, SystemExit) as error:
            # The outcome may be ambiguous.  Do not replay; list/reconcile only.
            mutation_error_class = type(error).__name__

        observed, observe_error_class = _reconcile_due(
            config, controller_module, base, timers
        )
        if observed is None:
            return 2, _uncertain_summary(
                config,
                mutation_called=True,
                mutation_error_class=mutation_error_class,
                observe_error_class=observe_error_class,
                call_marker_exists=True,
            )
        outcome = (
            "applied"
            if mutation_error_class is None
            else "applied_after_ambiguous_response"
        )
        result = _make_result(
            config,
            intent_payload,
            marker_payload,
            intent_on_disk,
            observed,
            outcome=outcome,
            mutation_called_this_invocation=True,
            mutation_response=mutation_response,
            mutation_error_class=mutation_error_class,
            observe_error_class=observe_error_class,
            timers=timers,
        )
        _write_success_result(config, result, intent_payload, marker_payload)
        return 0, result


def audit_prior_only(config: Config = DEFAULT_CONFIG) -> Dict[str, Any]:
    prior = validate_prior_renewal_chain(config)
    return {
        "status": "prior_renewal_chain_verified",
        "instance_id": config.instance_id,
        "prior_renewal_chain": prior,
        "current_due_time": config.old_due_time,
        "current_due_time_CST": "2026-08-21 19:10:16 +08:00",
        "target_due_time": config.new_due_time,
        "target_due_time_CST": "2026-08-26 19:10:16 +08:00",
        "network_used": False,
        "mutating_call_made": False,
    }


def self_check(config: Config = DEFAULT_CONFIG) -> Dict[str, Any]:
    source = Path(__file__).read_text(encoding="utf-8")
    renew_source = inspect.getsource(renew_once)
    reconcile_source = inspect.getsource(reconcile_only)
    loader_source = inspect.getsource(load_runtime)
    current_source = inspect.getsource(_current_instance)
    checks = {
        "instance_is_641496": config.instance_id == 641496,
        "exact_8xV100_contract": (
            config.gpu_type == "Tesla V100-32GB" and config.gpu_num == 8
        ),
        "exact_pre_expansion_disk_contract": (
            config.add_disk_GB == 4500 and config.adding_disk_GB == 0
        ),
        "exact_running_status": config.status == 1,
        "exact_old_due_time": config.old_due_time == 1787310616,
        "exact_120h_target": (
            config.new_due_time - config.old_due_time == 120 * 3600
        ),
        "cost_ceiling_2000_CNY": config.max_cost_CNY == Decimal("2000"),
        "exact_three_component_quote": (
            EXPECTED_FINAL_COST_CNY == Decimal("1589.376")
            and sum(
                value["final_cost_CNY"]
                for value in EXPECTED_QUOTE_COMPONENTS.values()
            ) == EXPECTED_FINAL_COST_CNY
        ),
        "intent_and_result_are_distinct": config.intent != config.result,
        "persistent_anchor_and_marker_are_distinct": len({
            config.intent,
            config.result,
            config.call_marker,
            config.prior_anchor,
        }) == 4,
        "shared_lock_is_instance_scoped": (
            config.mutation_lock
            == Path(
                "/home/zhangzelin/kx-4t/athenak-zhixing-private/641496/"
                "instance-641496.mutation.lock"
            )
        ),
        "one_mutating_endpoint_call_site": (
            renew_source.count("base.call_api(RENEW_ENDPOINT") == 1
        ),
        "reconcile_has_no_mutating_endpoint": "RENEW_ENDPOINT" not in reconcile_source,
        "intent_and_marker_precede_mutating_call": (
            source.index("write_private_exclusive(config.intent")
            < source.index("write_private_exclusive(config.call_marker")
            < source.index("base.call_api(RENEW_ENDPOINT")
        ),
        "pinned_persistent_runtime": (
            config.control.parent == PRIVATE_STATE_DIR
            and config.control_sha256 == API_BASE_SHA256
            and config.credentials.parent == PRIVATE_STATE_DIR
            and config.control_source == API_BASE_SOURCE
        ),
        "same_verified_controller_bytes_are_executed": (
            loader_source.index("hashlib.sha256(source).hexdigest()")
            < loader_source.index("compile(source")
        ),
        "credentials_are_bound_in_memory": (
            "module.load_credentials = load_verified_credentials" in loader_source
        ),
        "no_private_current_state_dependency": (
            "controller_module.STATE" not in current_source
            and "load_json" not in current_source
        ),
        "list_and_detail_identity_closure": (
            "list/detail safe instance views differ" in source
        ),
        "exclusive_nofollow_writes": (
            "os.O_EXCL" in inspect.getsource(write_bytes_exclusive)
            and "os.O_NOFOLLOW" in inspect.getsource(write_bytes_exclusive)
            and "os.fsync" in inspect.getsource(write_bytes_exclusive)
        ),
    }
    if not all(checks.values()):
        failed = ",".join(name for name, passed in checks.items() if not passed)
        raise ContractViolation(f"offline self-check failed: {failed}")
    return {
        "status": "offline_self_check_PASS",
        "checks": checks,
        "commands": [
            "self-check",
            "audit-prior",
            "install-runtime",
            "quote",
            "renew-once",
            "reconcile",
        ],
        "network_used": False,
        "mutating_call_made": False,
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    if len(arguments) != 1:
        raise SystemExit(
            "usage: zhixing_v100_mid_renew120_once.py "
            "{self-check|audit-prior|install-runtime|quote|renew-once|reconcile}"
        )
    command = arguments[0]
    try:
        if command == "self-check":
            code, output = 0, self_check()
        elif command == "audit-prior":
            code, output = 0, audit_prior_only()
        elif command == "install-runtime":
            code, output = 0, install_runtime()
        elif command == "quote":
            code, output = 0, quote_only()
        elif command == "renew-once":
            code, output = renew_once()
        elif command == "reconcile":
            code, output = reconcile_only()
        else:
            raise ContractViolation(f"unknown command: {command}")
    except ContractViolation as error:
        raise SystemExit(f"contract violation: {error}") from None
    print(json.dumps(output, ensure_ascii=False, indent=2, sort_keys=True))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
