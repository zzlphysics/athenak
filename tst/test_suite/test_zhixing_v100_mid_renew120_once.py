"""Mock-only tests for the one-shot Zhixing renewal transaction.

No test in this module imports the private cloud controller, opens a network
connection, or invokes SSH.
"""

from __future__ import annotations

from dataclasses import replace
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import sys
from typing import Any, Dict, List, Optional

import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "zhixing_v100_mid_renew120_once.py"
SPEC = importlib.util.spec_from_file_location(
    "zhixing_v100_mid_renew120_once", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
RENEWAL = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = RENEWAL
SPEC.loader.exec_module(RENEWAL)


INSTANCE_NAME = "mock-instance-name"
SECRET = "DO-NOT-PERSIST-FAKE-PASSWORD"


class FakeBase:
    ENDPOINTS = {"list": "mock-list", "detail": "mock-detail"}

    def __init__(
        self,
        states: List[Dict[str, Any]],
        quote: Optional[Dict[str, Any]] = None,
        renew_response: Optional[Dict[str, Any]] = None,
        renew_error: Optional[BaseException] = None,
        intent_path: Optional[Path] = None,
        marker_path: Optional[Path] = None,
        detail_changes: Optional[Dict[str, Any]] = None,
    ) -> None:
        self.states = list(states)
        self.quote_response = quote
        self.renew_response = renew_response or {
            "success": True,
            "code": "2000",
            "data": {
                "Init_passwd": SECRET,
                "RdpPasswd": SECRET,
                "Due_time": RENEWAL.EXPECTED_NEW_DUE_TIME,
            },
        }
        self.renew_error = renew_error
        self.intent_path = intent_path
        self.marker_path = marker_path
        self.detail_changes = detail_changes or {}
        self.calls: List[str] = []
        self.list_calls = 0
        self.detail_calls = 0
        self.quote_calls = 0
        self.renew_calls = 0
        self.last_list_state: Optional[Dict[str, Any]] = None

    def call_api(
        self, endpoint: str, _payload: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        self.calls.append(endpoint)
        if endpoint == self.ENDPOINTS["list"]:
            self.list_calls += 1
            if not self.states:
                raise AssertionError("mock list state was not configured")
            index = min(self.list_calls - 1, len(self.states) - 1)
            self.last_list_state = dict(self.states[index])
            return {
                "success": True,
                "code": "2000",
                "data": {"list": [self.last_list_state]},
                "secret": SECRET,
            }
        if endpoint == self.ENDPOINTS["detail"]:
            self.detail_calls += 1
            if self.last_list_state is None:
                raise AssertionError("detail was called before list")
            if _payload != {"instance_name": INSTANCE_NAME}:
                raise AssertionError("detail used an unexpected instance name")
            detail = dict(self.last_list_state)
            detail.update(self.detail_changes)
            detail["ServerTime"] = 999999
            return {
                "success": True,
                "code": "2000",
                "data": detail,
                "secret": SECRET,
            }
        if endpoint == RENEWAL.QUOTE_ENDPOINT:
            self.quote_calls += 1
            if self.quote_response is None:
                raise AssertionError("unexpected quote API call")
            return self.quote_response
        if endpoint == RENEWAL.RENEW_ENDPOINT:
            self.renew_calls += 1
            if self.intent_path is not None:
                assert self.intent_path.is_file()
                assert stat.S_IMODE(self.intent_path.stat().st_mode) == 0o600
                assert self.intent_path.read_bytes().endswith(b"\n")
            if self.marker_path is not None:
                assert self.marker_path.is_file()
                assert stat.S_IMODE(self.marker_path.stat().st_mode) == 0o600
                assert self.marker_path.read_bytes().endswith(b"\n")
            if self.renew_error is not None:
                raise self.renew_error
            return self.renew_response
        raise AssertionError(f"unexpected mock endpoint: {endpoint}")


class FakeClock:
    def __init__(self) -> None:
        self.now = 10.0

    def monotonic(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.now += seconds

    @staticmethod
    def wall_time() -> float:
        return 1234567890.25

    def timers(self) -> Any:
        return RENEWAL.Timers(
            monotonic=self.monotonic,
            sleep=self.sleep,
            wall_time=self.wall_time,
        )


def _write_private_json(path: Path, value: Dict[str, Any]) -> bytes:
    payload = (json.dumps(
        value, ensure_ascii=False, indent=2, sort_keys=True
    ) + "\n").encode("utf-8")
    path.write_bytes(payload)
    path.chmod(0o600)
    return payload


def _prior_contract(
    path: Path,
    old_due: int,
    new_due: int,
    pack_hours: int,
) -> Any:
    value = {
        "instance_id": RENEWAL.EXPECTED_INSTANCE_ID,
        "old_due_time": old_due,
        "new_due_time": new_due,
        "extension_seconds": pack_hours * 3600,
        "pack_hours": pack_hours,
        # These represent dangerous raw legacy fields.  Validation may read the
        # file but must never copy them into new intent/result evidence.
        "renew_response": {"Init_passwd": SECRET},
        "quote_response": {"UserInfoSheet": {"token": SECRET}},
    }
    payload = _write_private_json(path, value)
    return RENEWAL.PriorRenewal(
        path=path,
        sha256=hashlib.sha256(payload).hexdigest(),
        old_due_time=old_due,
        new_due_time=new_due,
        extension_seconds=pack_hours * 3600,
        pack_hours=pack_hours,
    )


@pytest.fixture
def config(tmp_path: Path) -> Any:
    first = _prior_contract(
        tmp_path / "prior120.private.json",
        1786734616,
        1787166616,
        120,
    )
    second = _prior_contract(
        tmp_path / "prior40.private.json",
        1787166616,
        RENEWAL.EXPECTED_OLD_DUE_TIME,
        40,
    )
    return replace(
        RENEWAL.DEFAULT_CONFIG,
        control=tmp_path / RENEWAL.API_BASE.name,
        control_sha256=RENEWAL.API_BASE_SHA256,
        control_source=tmp_path / "controller.source.py",
        credentials=tmp_path / RENEWAL.API_CREDENTIALS.name,
        credentials_source=tmp_path / "api.credentials.source",
        intent=tmp_path / RENEWAL.INTENT.name,
        result=tmp_path / RENEWAL.RESULT.name,
        call_marker=tmp_path / RENEWAL.CALL_MARKER.name,
        prior_anchor=tmp_path / RENEWAL.PRIOR_ANCHOR.name,
        mutation_lock=tmp_path / "zhixing-instance-641496.mutation.lock",
        prior_renewals=(first, second),
        observe_seconds=0.0,
        observe_interval_seconds=0.0,
        require_persistent_xfs=False,
    )


def live_item(config: Any, due_time: Optional[int] = None, **changes: Any) -> Dict[str, Any]:
    value = {
        "Id": config.instance_id,
        "Container_name": INSTANCE_NAME,
        "Status": config.status,
        "Gpu_type": config.gpu_type,
        "Gpu_num": config.gpu_num,
        "AddDisk": config.add_disk_GB,
        "AddingDisk": config.adding_disk_GB,
        "CanUpdateDisk": config.can_update_disk,
        "Due_time": config.old_due_time if due_time is None else due_time,
        "Init_passwd": SECRET,
        "RdpPasswd": SECRET,
        "Host": "secret-host.example.invalid",
    }
    value.update(changes)
    return value


def quote_response(
    config: Any,
    cost: Any = "1589.376",
    component_changes: Optional[Dict[str, Dict[str, Any]]] = None,
    **changes: Any,
) -> Dict[str, Any]:
    components = []
    for resource_type in ("gpu", "spec", "disk"):
        expected = RENEWAL.EXPECTED_QUOTE_COMPONENTS[resource_type]
        component = {
            "ResourceType": resource_type,
            "Type": expected["type"],
            "Title": expected["title"],
            "SingleHourPrice": str(expected["single_hour_price_CNY"]),
            "Amount": expected["amount"],
            "HourLen": config.pack_hours,
            "FinalDiscount": str(expected["final_discount"]),
            "FinalCost": str(expected["final_cost_CNY"]),
            "PayTypeFirst": "power",
            "Formula": SECRET,
        }
        component.update((component_changes or {}).get(resource_type, {}))
        components.append(component)
    sheet = {
        "ResourceType": "instance",
        "Type": "renew_instance",
        "HourLen": config.pack_hours,
        "PayTypeFirst": "power",
        "FinalCost": cost,
        "InstanceSubPriceSheet": components,
        "Formula": SECRET,
    }
    item = {
        "Due_time": config.old_due_time,
        "NewDuetime": config.new_due_time,
        "PayTypeFirst": "power",
    }
    for key, value in changes.items():
        if key.startswith("item_"):
            item[key[5:]] = value
        else:
            sheet[key] = value
    return {
        "success": True,
        "code": "2000",
        "message": SECRET,
        "data": {
            "PackPriceInfoSheet": sheet,
            "ItemInfoSheet": item,
            "UserInfoSheet": {"credential": SECRET},
        },
    }


def runtime(base: FakeBase):
    return lambda _config: (base, base)


def pinned_base_source(side_effect: str = "") -> bytes:
    return (
        "from pathlib import Path\n"
        + side_effect
        + f"BASE_URL = {RENEWAL.EXPECTED_API_BASE_URL!r}\n"
        + "CREDENTIALS = Path('/tmp/untrusted-original-credentials')\n"
        + "ENDPOINTS = {\n"
        + f"    'list': {RENEWAL.LIST_ENDPOINT!r},\n"
        + f"    'detail': {RENEWAL.DETAIL_ENDPOINT!r},\n"
        + "}\n"
        + "def load_credentials():\n"
        + "    return ('untrusted-access', 'untrusted-secret')\n"
        + "def call_api(endpoint, params=None):\n"
        + "    return {\n"
        + "        'success': True,\n"
        + "        'endpoint': endpoint,\n"
        + "        'params': params,\n"
        + "        'credentials': load_credentials(),\n"
        + "    }\n"
    ).encode("utf-8")


def test_offline_contract_and_prior_chain_are_exact(config: Any) -> None:
    checks = RENEWAL.self_check()
    assert checks["status"] == "offline_self_check_PASS"
    assert all(checks["checks"].values())

    audit = RENEWAL.audit_prior_only(config)
    assert audit["current_due_time"] == 1787310616
    assert audit["current_due_time_CST"] == "2026-08-21 19:10:16 +08:00"
    assert audit["target_due_time"] == 1787742616
    assert audit["target_due_time_CST"] == "2026-08-26 19:10:16 +08:00"
    assert audit["network_used"] is False


def test_pinned_runtime_executes_verified_bytes_not_reopened_path(config: Any) -> None:
    side_effect = "Path(__file__).write_text('tampered-after-read\\n', encoding='utf-8')\n"
    source = pinned_base_source(side_effect)
    config.control.write_bytes(source)
    config.control.chmod(0o600)
    config.credentials.write_text("accesskey\nsecretkey\n", encoding="utf-8")
    config.credentials.chmod(0o600)
    pinned = replace(
        config,
        control_sha256=hashlib.sha256(source).hexdigest(),
    )

    controller, base = RENEWAL.load_runtime(pinned)

    assert controller is base
    assert base.BASE_URL == RENEWAL.EXPECTED_API_BASE_URL
    assert base.CREDENTIALS == config.credentials
    assert base.load_credentials() == ("accesskey", "secretkey")
    assert config.control.read_text(encoding="utf-8") == "tampered-after-read\n"


def test_pinned_runtime_keeps_first_verified_credentials_in_memory(
    config: Any,
) -> None:
    source = pinned_base_source()
    config.control.write_bytes(source)
    config.control.chmod(0o600)
    config.credentials.write_text("firstaccess\nfirstsecret\n", encoding="utf-8")
    config.credentials.chmod(0o600)
    pinned = replace(
        config,
        control_sha256=hashlib.sha256(source).hexdigest(),
    )

    _, base = RENEWAL.load_runtime(pinned)
    first = base.call_api("before-replacement")
    config.credentials.unlink()
    config.credentials.write_text(
        "secondaccess\nsecondsecret\n", encoding="utf-8"
    )
    config.credentials.chmod(0o600)
    second = base.call_api("after-replacement")

    assert first["credentials"] == ("firstaccess", "firstsecret")
    assert second["credentials"] == ("firstaccess", "firstsecret")


def test_secure_read_rejects_hard_linked_private_file(tmp_path: Path) -> None:
    original = tmp_path / "private.json"
    hard_link = tmp_path / "private-hard-link.json"
    original.write_text("{}\n", encoding="utf-8")
    original.chmod(0o600)
    os.link(original, hard_link)

    with pytest.raises(RENEWAL.ContractViolation, match="one hard link"):
        RENEWAL._secure_read_bytes(original)


def test_wrong_controller_sha_fails_before_source_execution(config: Any) -> None:
    marker = config.intent.parent / "source-executed"
    side_effect = f"Path({str(marker)!r}).write_text('executed', encoding='utf-8')\n"
    source = pinned_base_source(side_effect)
    config.control.write_bytes(source)
    config.control.chmod(0o600)
    config.credentials.write_text("accesskey\nsecretkey\n", encoding="utf-8")
    config.credentials.chmod(0o600)

    with pytest.raises(RENEWAL.ContractViolation, match="SHA-256 changed"):
        RENEWAL.load_runtime(config)

    assert not marker.exists()


def test_install_runtime_copies_only_pinned_bytes_to_persistent_paths(
    config: Any,
) -> None:
    source = pinned_base_source()
    config.control_source.write_bytes(source)
    config.control_source.chmod(0o600)
    config.credentials_source.write_text(
        "mockaccess\nmocksecret\n", encoding="utf-8"
    )
    config.credentials_source.chmod(0o600)
    installed = replace(
        config,
        control_sha256=hashlib.sha256(source).hexdigest(),
    )

    first = RENEWAL.install_runtime(installed)
    second = RENEWAL.install_runtime(installed)

    assert first == second
    assert first["network_used"] is False
    assert first["mutating_api_call_made"] is False
    assert first["credentials_hash_persisted"] is False
    assert config.control.read_bytes() == source
    assert stat.S_IMODE(config.control.stat().st_mode) == 0o600
    assert stat.S_IMODE(config.credentials.stat().st_mode) == 0o600
    assert SECRET not in json.dumps(first)


def test_success_writes_intent_before_exactly_one_renew_and_proves_120h(
    config: Any, monkeypatch: pytest.MonkeyPatch,
) -> None:
    base = FakeBase(
        states=[
            live_item(config),
            live_item(config),
            live_item(config),
            live_item(config, due_time=config.new_due_time),
        ],
        quote=quote_response(config),
        intent_path=config.intent,
        marker_path=config.call_marker,
    )

    code, result = RENEWAL.renew_once(
        config, runtime(base), FakeClock().timers()
    )

    assert code == 0
    assert base.renew_calls == 1
    assert base.quote_calls == 1
    assert base.calls == [
        "mock-list",
        "mock-detail",
        RENEWAL.QUOTE_ENDPOINT,
        "mock-list",
        "mock-detail",
        "mock-list",
        "mock-detail",
        RENEWAL.RENEW_ENDPOINT,
        "mock-list",
        "mock-detail",
    ]
    assert result["status"] == "observed_exact_120_hour_extension"
    assert result["new_due_time"] - result["old_due_time"] == 120 * 3600
    assert result["observed"]["due_time"] == config.new_due_time
    assert stat.S_IMODE(config.intent.stat().st_mode) == 0o600
    assert stat.S_IMODE(config.prior_anchor.stat().st_mode) == 0o600
    assert stat.S_IMODE(config.call_marker.stat().st_mode) == 0o600
    assert stat.S_IMODE(config.result.stat().st_mode) == 0o600
    assert config.intent.read_bytes().endswith(b"\n")
    assert config.result.read_bytes().endswith(b"\n")
    combined = (
        config.prior_anchor.read_bytes()
        + config.intent.read_bytes()
        + config.call_marker.read_bytes()
        + config.result.read_bytes()
    )
    assert SECRET.encode() not in combined
    assert b"Init_passwd" not in combined
    assert b"RdpPasswd" not in combined
    intent = json.loads(config.intent.read_text(encoding="utf-8"))
    persisted = json.loads(config.result.read_text(encoding="utf-8"))
    assert intent["schema"] == "zhixing-renew-once-intent-v2"
    assert persisted["schema"] == "zhixing-renew-once-result-v2"
    assert intent["controller_sha256"] == RENEWAL.API_BASE_SHA256
    assert persisted["controller_sha256"] == RENEWAL.API_BASE_SHA256
    assert persisted["quoted_cost_CNY"] == "1589.376"
    monkeypatch.setattr(RENEWAL, "DEFAULT_CONFIG", config)
    validated, result_sha256 = RENEWAL.validate_persistent_result(
        config.intent, config.result, config.mutation_lock
    )
    assert validated == persisted
    assert result_sha256 == hashlib.sha256(config.result.read_bytes()).hexdigest()


def test_ambiguous_http_outcome_can_be_proven_without_replay(config: Any) -> None:
    base = FakeBase(
        states=[
            live_item(config),
            live_item(config),
            live_item(config),
            live_item(config, due_time=config.new_due_time),
        ],
        quote=quote_response(config),
        renew_error=TimeoutError("fake timeout containing " + SECRET),
        intent_path=config.intent,
        marker_path=config.call_marker,
    )

    code, result = RENEWAL.renew_once(
        config, runtime(base), FakeClock().timers()
    )

    assert code == 0
    assert base.renew_calls == 1
    assert result["outcome"] == "applied_after_ambiguous_response"
    assert result["mutation_error_class"] == "TimeoutError"
    assert SECRET.encode() not in config.result.read_bytes()


def test_unproven_http_outcome_leaves_intent_and_rerun_is_read_only(
    config: Any,
) -> None:
    first = FakeBase(
        states=[live_item(config), live_item(config)],
        quote=quote_response(config),
        renew_error=TimeoutError("ambiguous"),
        intent_path=config.intent,
    )
    code, uncertain = RENEWAL.renew_once(
        config, runtime(first), FakeClock().timers()
    )
    assert code == 2
    assert uncertain["result_written"] is False
    assert first.renew_calls == 1
    assert config.intent.is_file()
    assert not config.result.exists()

    second = FakeBase(
        states=[live_item(config, due_time=config.new_due_time)],
        quote=None,
    )
    code, reconciled = RENEWAL.renew_once(
        config, runtime(second), FakeClock().timers()
    )
    assert code == 0
    assert reconciled["outcome"] == "applied_reconciled_existing_intent"
    assert reconciled["mutation_called_this_invocation"] is False
    assert second.calls == ["mock-list", "mock-detail"]
    assert second.renew_calls == 0
    assert second.quote_calls == 0


def test_explicit_reconcile_never_references_quote_or_renew(config: Any) -> None:
    first = FakeBase(
        states=[live_item(config), live_item(config)],
        quote=quote_response(config),
        renew_error=ConnectionError("ambiguous"),
    )
    code, _ = RENEWAL.renew_once(
        config, runtime(first), FakeClock().timers()
    )
    assert code == 2

    observer = FakeBase(
        states=[live_item(config, due_time=config.new_due_time)],
        quote=None,
    )
    code, result = RENEWAL.reconcile_only(
        config, runtime(observer), FakeClock().timers()
    )
    assert code == 0
    assert result["mutation_called_this_invocation"] is False
    assert observer.calls == ["mock-list", "mock-detail"]


def test_persistent_prior_anchor_survives_loss_of_legacy_tmp_evidence(
    config: Any,
) -> None:
    first = FakeBase(
        states=[live_item(config), live_item(config)],
        quote=quote_response(config),
        renew_error=TimeoutError("ambiguous"),
    )
    code, _ = RENEWAL.renew_once(
        config, runtime(first), FakeClock().timers()
    )
    assert code == 2
    assert config.prior_anchor.is_file()
    for prior in config.prior_renewals:
        prior.path.unlink()

    observer = FakeBase(
        states=[live_item(config, due_time=config.new_due_time)],
        quote=None,
    )
    code, result = RENEWAL.reconcile_only(
        config, runtime(observer), FakeClock().timers()
    )

    assert code == 0
    assert result["outcome"] == "applied_reconciled_existing_intent"
    assert observer.renew_calls == 0
    assert observer.quote_calls == 0


def test_list_detail_mismatch_fails_before_quote_intent_or_mutation(
    config: Any,
) -> None:
    base = FakeBase(
        states=[live_item(config)],
        quote=quote_response(config),
        detail_changes={"Due_time": config.old_due_time + 1},
    )

    with pytest.raises(RENEWAL.ContractViolation, match="list/detail"):
        RENEWAL.renew_once(config, runtime(base), FakeClock().timers())

    assert base.calls == ["mock-list", "mock-detail"]
    assert base.quote_calls == 0
    assert base.renew_calls == 0
    assert not config.intent.exists()


def test_invalid_live_instance_name_is_rejected_before_detail(config: Any) -> None:
    base = FakeBase(
        states=[live_item(config, Container_name="../unsafe")],
        quote=quote_response(config),
    )

    with pytest.raises(RENEWAL.ContractViolation, match="name.*whitelist"):
        RENEWAL.renew_once(config, runtime(base), FakeClock().timers())

    assert base.calls == ["mock-list"]
    assert base.renew_calls == 0


def test_post_intent_drift_leaves_no_call_marker_and_can_never_replay(
    config: Any,
) -> None:
    first = FakeBase(
        states=[
            live_item(config),
            live_item(config),
            live_item(config, AddDisk=12500),
        ],
        quote=quote_response(config),
    )
    with pytest.raises(RENEWAL.ContractViolation, match="add_disk_GB"):
        RENEWAL.renew_once(config, runtime(first), FakeClock().timers())
    assert config.intent.is_file()
    assert not config.call_marker.exists()
    assert first.renew_calls == 0

    observer = FakeBase(
        states=[live_item(config, due_time=config.new_due_time)],
        quote=None,
    )
    code, output = RENEWAL.renew_once(
        config, runtime(observer), FakeClock().timers()
    )
    assert code == 3
    assert output["status"] == "target_observed_without_durable_call_marker"
    assert observer.calls == ["mock-list", "mock-detail"]
    assert observer.renew_calls == 0
    assert not config.result.exists()


@pytest.mark.parametrize(
    "changed",
    [
        {"Status": 8},
        {"Gpu_type": "NVIDIA A100-80GB"},
        {"Gpu_num": 4},
        {"AddDisk": 12500},
        {"AddingDisk": 8000},
        {"CanUpdateDisk": False},
        {"Due_time": RENEWAL.EXPECTED_OLD_DUE_TIME + 1},
    ],
)
def test_any_baseline_drift_fails_before_intent_or_mutation(
    config: Any, changed: Dict[str, Any]
) -> None:
    base = FakeBase(
        states=[live_item(config, **changed)],
        quote=quote_response(config),
    )
    with pytest.raises(RENEWAL.ContractViolation):
        RENEWAL.renew_once(config, runtime(base), FakeClock().timers())
    assert base.renew_calls == 0
    assert base.quote_calls == 0
    assert not config.intent.exists()
    assert not config.result.exists()


@pytest.mark.parametrize(
    "quote",
    [
        lambda c: quote_response(c, cost="2000.0001"),
        lambda c: quote_response(c, HourLen=119),
        lambda c: quote_response(c, item_NewDuetime=c.new_due_time + 1),
        lambda c: quote_response(c, ResourceType="disk"),
        lambda c: quote_response(c, Type="add_disk"),
    ],
)
def test_unsafe_quote_fails_without_arming_or_mutating(config: Any, quote) -> None:
    base = FakeBase(states=[live_item(config)], quote=quote(config))
    with pytest.raises(RENEWAL.ContractViolation):
        RENEWAL.renew_once(config, runtime(base), FakeClock().timers())
    assert base.quote_calls == 1
    assert base.renew_calls == 0
    assert not config.intent.exists()


@pytest.mark.parametrize(
    "component_changes",
    [
        {"disk": {"Amount": 1000}},
        {"disk": {"FinalCost": "72"}},
        {"gpu": {"FinalDiscount": "0.94"}},
        {"spec": {"Title": "实例规格gcs.g4.xlarge"}},
        {"gpu": {"Type": "NVIDIA A100-SXM4-80GB"}},
    ],
)
def test_quote_component_drift_requires_manual_review(
    config: Any, component_changes: Dict[str, Dict[str, Any]]
) -> None:
    base = FakeBase(
        states=[live_item(config)],
        quote=quote_response(config, component_changes=component_changes),
    )

    with pytest.raises(RENEWAL.ContractViolation):
        RENEWAL.renew_once(config, runtime(base), FakeClock().timers())

    assert base.quote_calls == 1
    assert base.renew_calls == 0
    assert not config.intent.exists()


def test_exact_component_price_is_allowed_for_read_only_quote(config: Any) -> None:
    base = FakeBase(
        states=[live_item(config), live_item(config)],
        quote=quote_response(config),
    )

    output = RENEWAL.quote_only(config, runtime(base))

    assert output["quote"]["final_cost_CNY"] == "1589.376"
    assert output["mutating_call_made"] is False
    assert base.renew_calls == 0
    assert not config.intent.exists()


def test_post_quote_resource_drift_fails_before_intent_and_mutation(
    config: Any,
) -> None:
    base = FakeBase(
        states=[live_item(config), live_item(config, AddDisk=12500)],
        quote=quote_response(config),
    )

    with pytest.raises(RENEWAL.ContractViolation, match="add_disk_GB"):
        RENEWAL.renew_once(config, runtime(base), FakeClock().timers())

    assert base.quote_calls == 1
    assert base.renew_calls == 0
    assert not config.intent.exists()


def test_nonexact_post_due_is_not_accepted_as_success(config: Any) -> None:
    base = FakeBase(
        states=[
            live_item(config),
            live_item(config),
            live_item(config),
            live_item(config, due_time=config.new_due_time + 1),
        ],
        quote=quote_response(config),
    )
    with pytest.raises(RENEWAL.ContractViolation, match="neither the exact old"):
        RENEWAL.renew_once(config, runtime(base), FakeClock().timers())
    assert base.renew_calls == 1
    assert config.intent.is_file()
    assert not config.result.exists()


def test_existing_intent_rejects_nonwhitelisted_field_before_any_api_call(
    config: Any,
) -> None:
    first = FakeBase(
        states=[live_item(config), live_item(config)],
        quote=quote_response(config),
        renew_error=TimeoutError("ambiguous"),
    )
    code, _ = RENEWAL.renew_once(
        config, runtime(first), FakeClock().timers()
    )
    assert code == 2
    value = json.loads(config.intent.read_text(encoding="utf-8"))
    value["raw_response"] = {"Init_passwd": SECRET}
    _write_private_json(config.intent, value)

    observer = FakeBase(
        states=[live_item(config, due_time=config.new_due_time)],
        quote=None,
    )
    with pytest.raises(RENEWAL.ContractViolation, match="non-whitelisted"):
        RENEWAL.reconcile_only(
            config, runtime(observer), FakeClock().timers()
        )
    assert observer.calls == []


@pytest.mark.parametrize(
    "tamper",
    [
        lambda value: value["resource_contract"].__setitem__("status", True),
        lambda value: value["before"].__setitem__("gpu_num", True),
        lambda value: value["prior_renewal_chain"][0].__setitem__(
            "pack_hours", True
        ),
        lambda value: value["quote"].__setitem__("hour_len", True),
        lambda value: value["quote"].__setitem__("final_cost_CNY", 1500.0),
    ],
)
def test_existing_intent_rejects_json_type_bypasses_before_api(
    config: Any, tamper
) -> None:
    first = FakeBase(
        states=[live_item(config), live_item(config)],
        quote=quote_response(config),
        renew_error=TimeoutError("ambiguous"),
    )
    code, _ = RENEWAL.renew_once(
        config, runtime(first), FakeClock().timers()
    )
    assert code == 2
    value = json.loads(config.intent.read_text(encoding="utf-8"))
    tamper(value)
    _write_private_json(config.intent, value)

    observer = FakeBase(
        states=[live_item(config, due_time=config.new_due_time)],
        quote=None,
    )
    with pytest.raises(RENEWAL.ContractViolation):
        RENEWAL.reconcile_only(
            config, runtime(observer), FakeClock().timers()
        )
    assert observer.calls == []


def test_tampered_call_marker_is_rejected_before_reconciliation_api(
    config: Any,
) -> None:
    first = FakeBase(
        states=[live_item(config), live_item(config)],
        quote=quote_response(config),
        renew_error=TimeoutError("ambiguous"),
    )
    code, _ = RENEWAL.renew_once(
        config, runtime(first), FakeClock().timers()
    )
    assert code == 2
    marker = json.loads(config.call_marker.read_text(encoding="utf-8"))
    marker["request_summary"]["pack_time"] = True
    _write_private_json(config.call_marker, marker)

    observer = FakeBase(
        states=[live_item(config, due_time=config.new_due_time)],
        quote=None,
    )
    with pytest.raises(RENEWAL.ContractViolation):
        RENEWAL.reconcile_only(
            config, runtime(observer), FakeClock().timers()
        )
    assert observer.calls == []


def test_tampered_prior_anchor_is_rejected_without_legacy_fallback_or_api(
    config: Any,
) -> None:
    first = FakeBase(
        states=[live_item(config), live_item(config)],
        quote=quote_response(config),
        renew_error=TimeoutError("ambiguous"),
    )
    code, _ = RENEWAL.renew_once(
        config, runtime(first), FakeClock().timers()
    )
    assert code == 2
    anchor = json.loads(config.prior_anchor.read_text(encoding="utf-8"))
    anchor["raw_response"] = {"Init_passwd": SECRET}
    _write_private_json(config.prior_anchor, anchor)

    observer = FakeBase(
        states=[live_item(config, due_time=config.new_due_time)],
        quote=None,
    )
    with pytest.raises(RENEWAL.ContractViolation, match="non-whitelisted"):
        RENEWAL.reconcile_only(
            config, runtime(observer), FakeClock().timers()
        )
    assert observer.calls == []


def test_persistent_result_validator_rejects_type_tamper(
    config: Any, monkeypatch: pytest.MonkeyPatch,
) -> None:
    base = FakeBase(
        states=[
            live_item(config),
            live_item(config),
            live_item(config),
            live_item(config, due_time=config.new_due_time),
        ],
        quote=quote_response(config),
    )
    code, _ = RENEWAL.renew_once(
        config, runtime(base), FakeClock().timers()
    )
    assert code == 0
    result = json.loads(config.result.read_text(encoding="utf-8"))
    result["observed"]["gpu_num"] = True
    _write_private_json(config.result, result)
    monkeypatch.setattr(RENEWAL, "DEFAULT_CONFIG", config)

    with pytest.raises(RENEWAL.ContractViolation):
        RENEWAL.validate_persistent_result(
            config.intent, config.result, config.mutation_lock
        )


def test_prior_audit_tamper_fails_before_runtime_or_mutation(config: Any) -> None:
    prior = config.prior_renewals[1].path
    prior.write_bytes(prior.read_bytes()[:-1] + b" \n")
    prior.chmod(0o600)
    runtime_loaded = False

    def forbidden_runtime(_config: Any):
        nonlocal runtime_loaded
        runtime_loaded = True
        raise AssertionError("runtime must not load")

    with pytest.raises(RENEWAL.ContractViolation, match="SHA-256"):
        RENEWAL.renew_once(config, forbidden_runtime, FakeClock().timers())
    assert runtime_loaded is False
    assert not config.intent.exists()


def test_shared_instance_lock_blocks_concurrent_mutation(config: Any) -> None:
    runtime_loaded = False

    def forbidden_runtime(_config: Any):
        nonlocal runtime_loaded
        runtime_loaded = True
        raise AssertionError("runtime must not load while lock is busy")

    with RENEWAL.acquire_instance_mutation_lock(config.mutation_lock):
        with pytest.raises(RENEWAL.ContractViolation, match="lock is busy"):
            RENEWAL.quote_only(config, forbidden_runtime)
    assert runtime_loaded is False
    assert stat.S_IMODE(config.mutation_lock.stat().st_mode) == 0o600


def test_exclusive_nofollow_write_refuses_symlink(tmp_path: Path) -> None:
    target = tmp_path / "target"
    target.write_text("unchanged", encoding="utf-8")
    target.chmod(0o600)
    link = tmp_path / "intent.private.json"
    link.symlink_to(target)

    with pytest.raises(FileExistsError):
        RENEWAL.write_private_exclusive(link, {"safe": True})

    assert target.read_text(encoding="utf-8") == "unchanged"


def test_exclusive_write_fsyncs_file_and_parent_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "intent.private.json"
    original_fsync = os.fsync
    synced_kinds: List[str] = []

    def traced_fsync(descriptor: int) -> None:
        mode = os.fstat(descriptor).st_mode
        synced_kinds.append("directory" if stat.S_ISDIR(mode) else "file")
        original_fsync(descriptor)

    monkeypatch.setattr(RENEWAL.os, "fsync", traced_fsync)
    RENEWAL.write_private_exclusive(path, {"schema": "test"})

    assert synced_kinds == ["file", "directory"]
    assert stat.S_IMODE(path.stat().st_mode) == 0o600


def test_result_path_is_exclusive_and_never_overwritten(config: Any) -> None:
    config.result.write_text("sentinel", encoding="utf-8")
    config.result.chmod(0o600)
    base = FakeBase(
        states=[live_item(config)],
        quote=quote_response(config),
    )

    with pytest.raises(RENEWAL.ContractViolation, match="result already exists"):
        RENEWAL.renew_once(config, runtime(base), FakeClock().timers())

    assert config.result.read_text(encoding="utf-8") == "sentinel"
    assert base.calls == []
