"""Security and durability tests for the Zhixing manifest receiver."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import stat
import sys
import time

import paramiko
import pytest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "zhixing_pull_manifest.py"
SPEC = importlib.util.spec_from_file_location("zhixing_pull_manifest", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
PULLER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = PULLER
SPEC.loader.exec_module(PULLER)


def test_host_key_policy_accepts_only_pinned_sha256() -> None:
    trusted = paramiko.RSAKey.generate(1024)
    attacker = paramiko.RSAKey.generate(1024)
    policy = PULLER.PinnedHostKeyPolicy(PULLER.key_sha256(trusted))

    assert policy.missing_host_key(object(), "gateway", trusted) is None
    with pytest.raises(paramiko.SSHException, match="host key mismatch"):
        policy.missing_host_key(object(), "gateway", attacker)


def test_connect_applies_hard_timeout_to_sftp_channel(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    trusted = paramiko.RSAKey.generate(1024)

    class Channel:
        timeout = None

        def settimeout(self, value: float) -> None:
            self.timeout = value

    class Sftp:
        def __init__(self) -> None:
            self.channel = Channel()

        def get_channel(self):
            return self.channel

    class Transport:
        def get_remote_server_key(self):
            return trusted

        def set_keepalive(self, _seconds: int) -> None:
            return None

    class Client:
        def __init__(self) -> None:
            self.sftp = Sftp()

        def set_missing_host_key_policy(self, _policy) -> None:
            return None

        def connect(self, **_kwargs) -> None:
            return None

        def get_transport(self):
            return Transport()

        def open_sftp(self):
            return self.sftp

        def close(self) -> None:
            return None

    monkeypatch.setattr(PULLER.paramiko, "SSHClient", Client)
    config = PULLER.ConnectionConfig(
        host="gateway",
        port=22,
        username="user",
        password="secret",
        host_key_sha256=PULLER.key_sha256(trusted),
    )

    client, sftp = PULLER.connect(config)

    assert client.sftp is sftp
    assert sftp.channel.timeout == PULLER.CHANNEL_IO_TIMEOUT_SECONDS


@pytest.mark.parametrize(
    "fingerprint",
    ["", "SHA256:", "not-base64", "SHA256:YWJj", "SHA256:abc def"],
)
def test_invalid_host_key_fingerprint_is_rejected(fingerprint: str) -> None:
    with pytest.raises(ValueError):
        PULLER.normalize_host_key_sha256(fingerprint)


def test_immutable_write_is_read_only_idempotent_and_no_replace(tmp_path: Path) -> None:
    record = tmp_path / "records" / "segment.ack"
    payload = b'{"schema":1}\n'

    assert PULLER.immutable_write(record, payload)
    original = record.lstat()
    assert original.st_nlink == 1
    assert stat.S_IMODE(original.st_mode) == 0o400
    assert not PULLER.immutable_write(record, payload)
    assert record.lstat().st_ino == original.st_ino
    with pytest.raises(RuntimeError, match="refusing changed"):
        PULLER.immutable_write(record, b"different\n")
    assert record.read_bytes() == payload
    assert record.lstat().st_ino == original.st_ino


def test_immutable_write_rejects_a_writable_existing_record(tmp_path: Path) -> None:
    record = tmp_path / "segment.ack"
    record.write_bytes(b"same\n")
    record.chmod(0o600)
    with pytest.raises(RuntimeError, match="writable"):
        PULLER.immutable_write(record, b"same\n")


def test_segment_lock_recovers_pinned_mode_zero_nfs_inode(tmp_path: Path) -> None:
    destination_path = tmp_path / "destination"
    locks = destination_path / ".locks"
    locks.mkdir(parents=True)
    lock_path = locks / "segment.lock"
    lock_path.write_text("pid=123\n", encoding="ascii")
    lock_path.chmod(0)

    with PULLER.DestinationRoot(destination_path) as destination:
        with PULLER.SegmentLock(destination, "segment"):
            assert stat.S_IMODE(lock_path.stat().st_mode) == 0o600
            assert lock_path.read_text(encoding="ascii") == f"pid={os.getpid()}\n"

    assert stat.S_IMODE(lock_path.stat().st_mode) == 0o600


def test_immutable_write_recovers_linked_temporary_alias(tmp_path: Path) -> None:
    record = tmp_path / "records" / "segment.ack"
    record.parent.mkdir()
    temporary = record.parent / f".{record.name}.123.deadbeef.tmp"
    payload = b"immutable\n"
    temporary.write_bytes(payload)
    temporary.chmod(0o400)
    os.link(temporary, record)
    assert record.lstat().st_nlink == 2

    assert not PULLER.immutable_write(record, payload)
    assert record.read_bytes() == payload
    assert record.lstat().st_nlink == 1
    assert not temporary.exists()


def test_immutable_write_at_recovers_linked_temporary_alias(tmp_path: Path) -> None:
    name = "segment.ack"
    temporary_name = f".{name}.123.deadbeef.tmp"
    payload = b"immutable-at\n"
    temporary = tmp_path / temporary_name
    record = tmp_path / name
    temporary.write_bytes(payload)
    temporary.chmod(0o400)
    os.link(temporary, record)
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        assert not PULLER.immutable_write_at(parent, name, payload)
    finally:
        os.close(parent)
    assert record.read_bytes() == payload
    assert record.lstat().st_nlink == 1
    assert not temporary.exists()


def test_regular_hash_rejects_symlinks_and_hardlinks(tmp_path: Path) -> None:
    source = tmp_path / "source.bin"
    source.write_bytes(b"payload")
    symlink = tmp_path / "symlink.bin"
    symlink.symlink_to(source)
    with pytest.raises(RuntimeError, match="one-link regular"):
        PULLER.sha256_regular(symlink)

    hardlink = tmp_path / "hardlink.bin"
    os.link(source, hardlink)
    with pytest.raises(RuntimeError, match="one-link regular"):
        PULLER.sha256_regular(source)


def test_regular_hash_can_evict_verified_payload_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    payload = tmp_path / "payload.bin"
    payload.write_bytes(b"0123456789")
    monkeypatch.setattr(PULLER, "CHUNK", 4)
    monkeypatch.setattr(PULLER, "HASH_CACHE_WINDOW", 4)
    calls: list[tuple[int, int, int]] = []
    monkeypatch.setattr(
        PULLER.os,
        "posix_fadvise",
        lambda _descriptor, offset, length, advice:
            calls.append((offset, length, advice)),
    )
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        PULLER.sha256_regular_at(parent, payload.name, evict_after=True)
    finally:
        os.close(parent)

    assert calls == [
        (0, 4, PULLER.os.POSIX_FADV_DONTNEED),
        (4, 4, PULLER.os.POSIX_FADV_DONTNEED),
        (0, 0, PULLER.os.POSIX_FADV_DONTNEED),
    ]


def test_secure_destination_rejects_symlinked_payload_parent(tmp_path: Path) -> None:
    destination = tmp_path / "destination"
    outside = tmp_path / "outside"
    destination.mkdir()
    outside.mkdir()
    (destination / "state").symlink_to(outside)
    with pytest.raises(OSError):
        PULLER.secure_relative_parent(
            destination, PULLER.PurePosixPath("state/output.bin")
        )


def test_pinned_destination_parent_survives_path_swap(tmp_path: Path) -> None:
    destination = tmp_path / "destination"
    state = destination / "state"
    destination.mkdir()
    state.mkdir()
    with PULLER.DestinationRoot(destination) as root:
        parent, name = root.open_relative_parent(
            PULLER.PurePosixPath("state/output.bin")
        )
        moved = destination / "state-moved"
        state.rename(moved)
        state.mkdir()
        try:
            stream, size = PULLER.open_append_nofollow_at(parent, name)
            assert size == 0
            with stream:
                stream.write(b"pinned")
        finally:
            os.close(parent)
    assert (moved / "output.bin").read_bytes() == b"pinned"
    assert not (state / "output.bin").exists()


def test_visible_directory_rebind_rejects_renamed_ack_parent(tmp_path: Path) -> None:
    destination = tmp_path / "destination"
    destination.mkdir()
    with PULLER.DestinationRoot(destination) as root:
        ack_parent = root.open_directory((".acks",), create=True)
        moved = destination / ".acks-moved"
        (destination / ".acks").rename(moved)
        (destination / ".acks").mkdir()
        try:
            with pytest.raises(RuntimeError, match="pathname changed"):
                root.assert_directory_visible((".acks",), ack_parent)
        finally:
            os.close(ack_parent)


def test_manifest_validation_reserves_receiver_metadata() -> None:
    manifest = {
        "schema": 1,
        "root": "/remote/root",
        "segment": "segment-1",
        "files": [{"path": ".acks/forged", "size": 0, "sha256": "0" * 64}],
    }
    with pytest.raises(ValueError, match="collides"):
        PULLER.validate_manifest(
            manifest, "segment-1.manifest.ready", "/remote/root", "segment-1"
        )


def test_manifest_validation_rejects_bool_size() -> None:
    manifest = {
        "schema": 1,
        "root": "/remote/root",
        "segment": "segment-1",
        "files": [{"path": "state/a", "size": True, "sha256": "0" * 64}],
    }
    with pytest.raises(ValueError, match="invalid size"):
        PULLER.validate_manifest(
            manifest, "segment-1.manifest.ready", "/remote/root", "segment-1"
        )


@pytest.mark.parametrize(
    "value",
    ["relative/path", "/", "/root/../escape", "/root//file", "/root/./file"],
)
def test_remote_paths_must_be_canonical_absolute(value: str) -> None:
    with pytest.raises(ValueError):
        PULLER.canonical_remote_absolute(value, "remote path")


def _ack_core() -> dict[str, object]:
    return {
        "schema": 1,
        "kind": "athenak_transfer_ack_v2",
        "authorizes_remote_deletion": False,
        "manifest": "segment-1.manifest.ready",
        "manifest_sha256": "a" * 64,
        "segment": "segment-1",
        "destination": "/nas/segment-1",
        "mount": {
            "target": "/nas",
            "source": "server:/nas",
            "fstype": "nfs",
            "options": "rw",
        },
        "file_count": 1,
        "total_bytes": 7,
        "files": [{"path": "state/a", "size": 7, "sha256": "b" * 64}],
    }


def test_production_ack_is_explicitly_non_authorizing() -> None:
    manifest_bytes = b'{"schema":1}\n'
    records = [{
        "path": PULLER.PurePosixPath("state/bin/output.bin"),
        "size": 7,
        "sha256": "b" * 64,
    }]
    core = PULLER.ack_core(
        "segment-1.manifest.ready",
        manifest_bytes,
        "segment-1",
        Path("/nas/segment-1"),
        _ack_core()["mount"],
        records,
    )

    assert core["kind"] == "athenak_transfer_ack_v2"
    assert core["authorizes_remote_deletion"] is False
    encoded = PULLER.choose_ack_bytes(core, None, None, now=123.5)
    assert json.loads(encoded)["authorizes_remote_deletion"] is False

    forged = json.loads(encoded)
    forged["authorizes_remote_deletion"] = True
    forged_bytes = (
        json.dumps(forged, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    with pytest.raises(RuntimeError, match="does not match"):
        PULLER.validate_ack_bytes(forged_bytes, core)

    legacy = json.loads(encoded)
    legacy.pop("kind")
    legacy.pop("authorizes_remote_deletion")
    legacy_bytes = (
        json.dumps(legacy, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    with pytest.raises(RuntimeError, match="unexpected fields"):
        PULLER.validate_ack_bytes(legacy_bytes, core)


def test_ack_reuses_one_sided_transaction_bytes() -> None:
    core = _ack_core()
    first = PULLER.choose_ack_bytes(core, None, None, now=123.5)
    assert json.loads(first)["verified_unix"] == 123.5
    assert PULLER.choose_ack_bytes(core, first, None, now=999.0) == first
    assert PULLER.choose_ack_bytes(core, None, first, now=999.0) == first
    assert PULLER.choose_ack_bytes(core, first, first, now=999.0) == first


def test_ack_refuses_local_remote_divergence() -> None:
    core = _ack_core()
    local = PULLER.choose_ack_bytes(core, None, None, now=123.0)
    remote = PULLER.choose_ack_bytes(core, None, None, now=124.0)
    with pytest.raises(RuntimeError, match="ACKs differ"):
        PULLER.choose_ack_bytes(core, local, remote)


def test_ack_refuses_changed_destination_even_with_valid_json() -> None:
    core = _ack_core()
    encoded = PULLER.choose_ack_bytes(core, None, None, now=123.0)
    changed = dict(core)
    changed["destination"] = "/local-disk"
    with pytest.raises(RuntimeError, match="does not match"):
        PULLER.validate_ack_bytes(encoded, changed)


def test_ack_rejects_duplicate_keys_and_noncanonical_bytes() -> None:
    core = _ack_core()
    encoded = PULLER.choose_ack_bytes(core, None, None, now=123.0)
    ambiguous = encoded.replace(
        b'  "authorizes_remote_deletion": false,\n',
        b'  "authorizes_remote_deletion": true,\n'
        b'  "authorizes_remote_deletion": false,\n',
    )
    with pytest.raises(RuntimeError, match="invalid existing ACK JSON"):
        PULLER.validate_ack_bytes(ambiguous, core)

    compact = json.dumps(json.loads(encoded), sort_keys=True).encode("utf-8")
    with pytest.raises(RuntimeError, match="not canonical"):
        PULLER.validate_ack_bytes(compact, core)


def test_final_verification_hashes_installed_path_and_fsyncs(tmp_path: Path) -> None:
    destination = tmp_path / "destination"
    final = destination / "state" / "closed.bin"
    final.parent.mkdir(parents=True)
    final.write_bytes(b"closed")
    digest = PULLER.hashlib.sha256(b"closed").hexdigest()
    records = [{"path": PULLER.PurePosixPath("state/closed.bin"),
                "size": 6, "sha256": digest}]

    PULLER.verify_final_files(destination, records, seal_read_only=False)
    final.write_bytes(b"damage")
    with pytest.raises(RuntimeError, match="manifest verification failed"):
        PULLER.verify_final_files(destination, records, seal_read_only=False)


def test_final_verification_seals_payload_read_only(tmp_path: Path) -> None:
    destination = tmp_path / "destination"
    final = destination / "state" / "closed.bin"
    final.parent.mkdir(parents=True)
    final.write_bytes(b"closed")
    records = [{
        "path": PULLER.PurePosixPath("state/closed.bin"),
        "size": 6,
        "sha256": PULLER.hashlib.sha256(b"closed").hexdigest(),
    }]

    proofs = PULLER.verify_final_files(destination, records)
    assert stat.S_IMODE(final.stat().st_mode) == 0o400
    assert "state/closed.bin" in proofs
    with pytest.raises(PermissionError):
        final.write_bytes(b"damage")


def test_partial_creation_normalizes_mode_despite_umask(tmp_path: Path) -> None:
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    previous = os.umask(0o777)
    try:
        stream, size = PULLER.open_append_nofollow_at(parent, "payload.part")
    finally:
        os.umask(previous)
    try:
        assert size == 0
        assert stat.S_IMODE(os.fstat(stream.fileno()).st_mode) == 0o600
    finally:
        stream.close()
        os.close(parent)


def test_final_seal_retains_owner_read_for_write_only_payload(
    tmp_path: Path,
) -> None:
    destination = tmp_path / "destination"
    final = destination / "state" / "closed.bin"
    final.parent.mkdir(parents=True)
    final.write_bytes(b"closed")
    records = [{
        "path": PULLER.PurePosixPath("state/closed.bin"),
        "size": 6,
        "sha256": PULLER.hashlib.sha256(b"closed").hexdigest(),
    }]
    # Exercise the production manifest-aware recovery for legacy finals that
    # the previous sealing code accidentally made mode 000.
    os.chmod(final, 0o000)
    PULLER.verify_final_files(destination, records)
    assert stat.S_IMODE(final.stat().st_mode) == 0o400
    assert final.read_bytes() == b"closed"


def test_mode_zero_recovery_refuses_wrong_size_hardlink_and_nonzero_mode(
    tmp_path: Path,
) -> None:
    final = tmp_path / "closed.bin"
    final.write_bytes(b"closed")
    digest = PULLER.hashlib.sha256(b"closed").hexdigest()
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.chmod(final, 0o000)
        with pytest.raises(RuntimeError, match="unsafe legacy"):
            PULLER.hash_manifest_final_at(parent, final.name, 7, digest)
        assert stat.S_IMODE(final.stat().st_mode) == 0o000

        hardlink = tmp_path / "alias.bin"
        os.link(final, hardlink)
        with pytest.raises(RuntimeError, match="one-link"):
            PULLER.hash_manifest_final_at(parent, final.name, 6, digest)
        assert stat.S_IMODE(final.stat().st_mode) == 0o000
        hardlink.unlink()

        os.chmod(final, 0o200)
        with pytest.raises(PermissionError):
            PULLER.hash_manifest_final_at(parent, final.name, 6, digest)
        assert stat.S_IMODE(final.stat().st_mode) == 0o200
    finally:
        os.chmod(final, 0o600)
        os.close(parent)


def test_mode_zero_wrong_digest_is_retained_read_only(tmp_path: Path) -> None:
    final = tmp_path / "closed.bin"
    final.write_bytes(b"closed")
    final.chmod(0o000)
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with pytest.raises(RuntimeError, match="manifest verification failed"):
            PULLER.hash_manifest_final_at(parent, final.name, 6, "0" * 64)
    finally:
        os.close(parent)
    assert final.read_bytes() == b"closed"
    assert stat.S_IMODE(final.stat().st_mode) == 0o400


def test_preopened_writer_can_invalidate_payload_but_not_authorize_cleanup(
    tmp_path: Path,
) -> None:
    """Regress the ACK-publication attack window from the P1 review."""

    destination = tmp_path / "destination"
    final = destination / "state" / "bin" / "closed.bin"
    final.parent.mkdir(parents=True)
    final.write_bytes(b"closed")
    record = {
        "path": PULLER.PurePosixPath("state/bin/closed.bin"),
        "size": 6,
        "sha256": PULLER.hashlib.sha256(b"closed").hexdigest(),
    }
    writer = os.open(final, os.O_RDWR | os.O_CLOEXEC)
    try:
        with PULLER.DestinationRoot(destination) as root:
            proofs = PULLER.verify_final_files(root, [record])
            assert stat.S_IMODE(final.stat().st_mode) == 0o400
            core = PULLER.ack_core(
                "segment-1.manifest.ready",
                b'{"schema":1}\n',
                "segment-1",
                destination,
                {
                    "target": str(destination),
                    "source": "server:/nas",
                    "fstype": "nfs",
                    "options": "rw",
                },
                [record],
            )
            published_receipt = PULLER.choose_ack_bytes(
                core,
                None,
                None,
                now=123.5,
            )
            remote_receipt = tmp_path / "remote" / "segment-1.transfer.ack"
            assert PULLER.immutable_write(remote_receipt, published_receipt)

            # Some test filesystems expose coarse timestamp updates.  Cross one
            # timestamp tick so the existing identity-only postcheck observes the
            # attack deterministically; the changed bytes are the security issue.
            time.sleep(0.02)
            os.pwrite(writer, b"damage", 0)
            os.fsync(writer)
            with pytest.raises(RuntimeError, match="identity changed"):
                PULLER.verify_final_identities(root, proofs)
    finally:
        os.close(writer)

    # The receipt can remain after the failed postcheck, but it cannot be
    # interpreted as a cloud-deletion capability.
    assert remote_receipt.is_file()
    assert stat.S_IMODE(remote_receipt.stat().st_mode) == 0o400
    assert json.loads(published_receipt)["authorizes_remote_deletion"] is False
    assert json.loads(remote_receipt.read_bytes())["authorizes_remote_deletion"] is False


def test_stale_parallel_assembly_is_replaced(tmp_path: Path, monkeypatch) -> None:
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    partial_name = "dump.part"
    monkeypatch.setattr(PULLER, "CHUNK", 4)
    monkeypatch.setattr(PULLER, "ASSEMBLY_CACHE_WINDOW", 4)
    piece0 = tmp_path / f"{partial_name}.piece00"
    piece1 = tmp_path / f"{partial_name}.piece01"
    piece0.write_bytes(b"a" * PULLER.CHUNK)
    piece1.write_bytes(b"b")
    (tmp_path / f"{partial_name}.assembling").write_bytes(b"stale")
    monkeypatch.setattr(PULLER, "PARALLEL_STREAMS", 2)
    evictions: list[tuple[str, int, int]] = []

    def record_eviction(descriptor: int, offset: int, length: int) -> None:
        evictions.append(
            (Path(os.readlink(f"/proc/self/fd/{descriptor}")).name, offset, length)
        )

    def fake_download_verified_piece(
        _config,
        _remote,
        descriptor,
        name,
        start,
        end,
        number,
        _validation,
        _identity,
    ):
        result = PULLER.sha256_regular_at(descriptor, name)
        return PULLER.RangeProof(
            number,
            start,
            end,
            name,
            result[1],
            PULLER._identity(result[0]),
        )

    monkeypatch.setattr(PULLER, "evict_cached_range", record_eviction)
    monkeypatch.setattr(
        PULLER,
        "download_verified_piece",
        fake_download_verified_piece,
    )
    try:
        PULLER.parallel_download(
            object(),
            "/remote/dump",
            parent,
            partial_name,
            PULLER.CHUNK + 1,
        )
    finally:
        os.close(parent)
    assert (tmp_path / partial_name).read_bytes() == b"a" * PULLER.CHUNK + b"b"
    assert not (tmp_path / f"{partial_name}.assembling").exists()
    assert piece0.exists() and piece1.exists()
    assembly_evictions = [
        (offset, length)
        for name, offset, length in evictions
        if name == f"{partial_name}.assembling"
    ]
    assert assembly_evictions == [(0, 4), (4, 1), (0, 0)]


def test_parallel_assembly_retries_without_redownloading(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    partial_name = "dump.part"
    payload = b"aaaab"
    monkeypatch.setattr(PULLER, "CHUNK", 4)
    monkeypatch.setattr(PULLER, "ASSEMBLY_CACHE_WINDOW", 4)
    monkeypatch.setattr(PULLER, "PARALLEL_STREAMS", 2)
    piece0 = tmp_path / f"{partial_name}.piece00"
    piece1 = tmp_path / f"{partial_name}.piece01"
    piece0.write_bytes(payload[:4])
    piece1.write_bytes(payload[4:])
    downloads: list[str] = []

    def prove_existing(
        _config, _remote, descriptor, name, start, end, number,
        _validation, _identity,
    ) -> PULLER.RangeProof:
        downloads.append(name)
        info, digest = PULLER.sha256_regular_at(descriptor, name)
        return PULLER.RangeProof(
            number,
            start,
            end,
            name,
            digest,
            PULLER._identity(info),
        )

    monkeypatch.setattr(PULLER, "download_verified_piece", prove_existing)
    original_write = PULLER._write_parallel_assembly
    assembly_attempts = 0

    def corrupt_first_assembly(descriptor, name, proofs) -> int:
        nonlocal assembly_attempts
        assembly_attempts += 1
        written = original_write(descriptor, name, proofs)
        if assembly_attempts == 1:
            assembly = os.open(
                name,
                os.O_WRONLY | os.O_NOFOLLOW | os.O_CLOEXEC,
                dir_fd=descriptor,
            )
            try:
                assert os.pwrite(assembly, b"x", 0) == 1
                os.fsync(assembly)
            finally:
                os.close(assembly)
        return written

    monkeypatch.setattr(
        PULLER,
        "_write_parallel_assembly",
        corrupt_first_assembly,
    )
    try:
        PULLER.parallel_download(
            object(),
            "/remote/dump",
            parent,
            partial_name,
            len(payload),
            expected_digest=PULLER.hashlib.sha256(payload).hexdigest(),
        )
    finally:
        os.close(parent)

    assert assembly_attempts == 2
    assert sorted(downloads) == [piece0.name, piece1.name]
    assert (tmp_path / partial_name).read_bytes() == payload
    assert piece0.read_bytes() + piece1.read_bytes() == payload
    assert not (tmp_path / f"{partial_name}.assembling").exists()


def test_parallel_assembly_retry_limit_leaves_mismatch_for_quarantine(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    partial_name = "dump.part"
    payload = b"aaaab"
    monkeypatch.setattr(PULLER, "CHUNK", 4)
    monkeypatch.setattr(PULLER, "ASSEMBLY_CACHE_WINDOW", 4)
    monkeypatch.setattr(PULLER, "PARALLEL_STREAMS", 2)
    monkeypatch.setattr(PULLER, "PARALLEL_ASSEMBLY_ATTEMPTS", 2)
    for number, content in enumerate((payload[:4], payload[4:])):
        (tmp_path / f"{partial_name}.piece{number:02d}").write_bytes(content)

    def prove_existing(
        _config, _remote, descriptor, name, start, end, number,
        _validation, _identity,
    ) -> PULLER.RangeProof:
        info, digest = PULLER.sha256_regular_at(descriptor, name)
        return PULLER.RangeProof(
            number,
            start,
            end,
            name,
            digest,
            PULLER._identity(info),
        )

    monkeypatch.setattr(PULLER, "download_verified_piece", prove_existing)
    original_write = PULLER._write_parallel_assembly
    assembly_attempts = 0

    def corrupt_every_assembly(descriptor, name, proofs) -> int:
        nonlocal assembly_attempts
        assembly_attempts += 1
        written = original_write(descriptor, name, proofs)
        assembly = os.open(
            name,
            os.O_WRONLY | os.O_NOFOLLOW | os.O_CLOEXEC,
            dir_fd=descriptor,
        )
        try:
            assert os.pwrite(assembly, b"x", 0) == 1
            os.fsync(assembly)
        finally:
            os.close(assembly)
        return written

    monkeypatch.setattr(
        PULLER,
        "_write_parallel_assembly",
        corrupt_every_assembly,
    )
    try:
        PULLER.parallel_download(
            object(),
            "/remote/dump",
            parent,
            partial_name,
            len(payload),
            expected_digest=PULLER.hashlib.sha256(payload).hexdigest(),
        )
    finally:
        os.close(parent)

    assert assembly_attempts == 2
    assert (tmp_path / partial_name).read_bytes() == b"xaaab"
    assert not (tmp_path / f"{partial_name}.assembling").exists()
    assert len(list(tmp_path.glob(f"{partial_name}.piece*"))) == 2


def test_parallel_ranges_are_ordered_and_uneven(monkeypatch) -> None:
    monkeypatch.setattr(PULLER, "CHUNK", 4)
    monkeypatch.setattr(PULLER, "PARALLEL_STREAMS", 2)

    assert PULLER.parallel_ranges(13, "payload.part") == [
        (0, 0, 8, "payload.part.piece00"),
        (1, 8, 13, "payload.part.piece01"),
    ]


def test_cleanup_parallel_pieces_rejects_same_size_replacement(tmp_path: Path) -> None:
    piece = tmp_path / "payload.part.piece00"
    piece.write_bytes(b"good")
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        info, digest = PULLER.sha256_regular_at(parent, piece.name)
        proof = PULLER.RangeProof(
            0, 0, 4, piece.name, digest, PULLER._identity(info)
        )
        piece.unlink()
        piece.write_bytes(b"evil")
        with pytest.raises(RuntimeError, match="changed before cleanup"):
            PULLER.cleanup_parallel_pieces(parent, [proof])
    finally:
        os.close(parent)
    assert piece.read_bytes() == b"evil"


def test_quarantine_is_read_only_and_idempotent_after_link_crash(
    tmp_path: Path,
) -> None:
    active = tmp_path / "payload.part"
    evidence = tmp_path / "payload.bad"
    active.write_bytes(b"bad-payload")
    digest = PULLER.hashlib.sha256(b"bad-payload").hexdigest()
    active.chmod(0o400)
    os.link(active, evidence)
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        assert PULLER.quarantine_regular_at(
            parent, active.name, evidence.name, len(b"bad-payload"), digest
        ) == evidence.name
    finally:
        os.close(parent)
    assert not active.exists()
    assert evidence.read_bytes() == b"bad-payload"
    assert stat.S_IMODE(evidence.stat().st_mode) == 0o400
    assert evidence.stat().st_nlink == 1


def test_remote_range_hash_uses_exact_pinned_range(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    payload = b"abcdefghijklm"
    attributes = _DownloadAttributes(len(payload))

    class Sftp:
        def lstat(self, _path: str):
            return attributes

    class Session:
        def remote_call(self, _label: str, operation):
            return operation(object(), Sftp())

    captured: list[list[str]] = []

    def capture(_client, argv, **_kwargs):
        captured.append(argv)
        return PULLER.hashlib.sha256(payload[8:13]).hexdigest().encode() + b"  -\n", b""

    monkeypatch.setattr(PULLER, "exec_capture_checked", capture)
    digest = PULLER.remote_range_sha256(
        Session(),
        "/remote/payload with space",
        PULLER._remote_identity(attributes),
        8,
        13,
        "range 1",
    )

    assert digest == PULLER.hashlib.sha256(payload[8:13]).hexdigest()
    assert captured[0][:4] == ["/bin/bash", "-o", "pipefail", "-c"]
    assert "skip=8 count=5" in captured[0][4]
    assert "'/remote/payload with space'" in captured[0][4]


def test_range_hash_mismatch_quarantines_only_bad_piece(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    payload = b"local-range"
    piece = tmp_path / "payload.part.piece01"
    piece.write_bytes(payload)
    validation = _DownloadValidation(
        [_DownloadClient(payload) for _ in range(3)], len(payload)
    )
    monkeypatch.setattr(
        PULLER,
        "remote_range_sha256",
        lambda *_args: PULLER.hashlib.sha256(b"remote-range").hexdigest(),
    )
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with pytest.raises(PULLER.RangeHashMismatch, match="range 1"):
            PULLER.download_verified_piece(
                object(),
                "/remote/payload",
                os.dup(parent),
                piece.name,
                0,
                len(payload),
                1,
                validation,
                _remote_download_identity(len(payload)),
            )
    finally:
        os.close(parent)

    assert not piece.exists()
    evidence = list(tmp_path.glob("range-*.bad"))
    assert len(evidence) == 1
    assert evidence[0].read_bytes() == payload
    assert stat.S_IMODE(evidence[0].stat().st_mode) == 0o400


def test_transient_initial_range_mismatch_requires_two_matching_confirmations(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    payload = b"stable-range"
    expected_digest = PULLER.hashlib.sha256(payload).hexdigest()
    piece = tmp_path / "payload.part.piece00"
    piece.write_bytes(payload)
    validation = _DownloadValidation(
        [_DownloadClient(payload) for _ in range(3)], len(payload)
    )
    real_hash = PULLER.sha256_regular_at
    local_calls = 0
    remote_calls = 0

    def transient_first_local(*args, **kwargs):
        nonlocal local_calls
        info, digest = real_hash(*args, **kwargs)
        local_calls += 1
        if local_calls == 1:
            return info, PULLER.hashlib.sha256(b"transient-read").hexdigest()
        return info, digest

    def stable_remote(*_args):
        nonlocal remote_calls
        remote_calls += 1
        return expected_digest

    monkeypatch.setattr(PULLER, "sha256_regular_at", transient_first_local)
    monkeypatch.setattr(PULLER, "remote_range_sha256", stable_remote)
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        proof = PULLER.download_verified_piece(
            object(),
            "/remote/payload",
            os.dup(parent),
            piece.name,
            0,
            len(payload),
            0,
            validation,
            _remote_download_identity(len(payload)),
        )
    finally:
        os.close(parent)

    assert proof.sha256 == expected_digest
    assert local_calls == 3
    assert remote_calls == 3
    assert piece.read_bytes() == payload
    assert stat.S_IMODE(piece.stat().st_mode) == 0o400
    assert not list(tmp_path.glob("range-*.bad"))


def test_unstable_range_proof_preserves_piece_without_quarantine(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    payload = b"preserve-range"
    piece = tmp_path / "payload.part.piece00"
    piece.write_bytes(payload)
    validation = _DownloadValidation(
        [_DownloadClient(payload) for _ in range(3)], len(payload)
    )
    real_hash = PULLER.sha256_regular_at
    fake_digests = iter(
        PULLER.hashlib.sha256(value).hexdigest()
        for value in (b"first-read", b"second-read", b"third-read")
    )

    def unstable_local(*args, **kwargs):
        info, _digest = real_hash(*args, **kwargs)
        return info, next(fake_digests)

    monkeypatch.setattr(PULLER, "sha256_regular_at", unstable_local)
    monkeypatch.setattr(
        PULLER,
        "remote_range_sha256",
        lambda *_args: PULLER.hashlib.sha256(b"remote-range").hexdigest(),
    )
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        with pytest.raises(PULLER.RangeProofInstability, match="preserving"):
            PULLER.download_verified_piece(
                object(),
                "/remote/payload",
                os.dup(parent),
                piece.name,
                0,
                len(payload),
                0,
                validation,
                _remote_download_identity(len(payload)),
            )
    finally:
        os.close(parent)

    assert piece.read_bytes() == payload
    assert stat.S_IMODE(piece.stat().st_mode) == 0o400
    assert not list(tmp_path.glob("range-*.bad"))


def test_prove_existing_parallel_pieces_only_proves_present_ranges(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(PULLER, "CHUNK", 4)
    monkeypatch.setattr(PULLER, "PARALLEL_STREAMS", 2)
    payload = b"abcdefghijklm"
    piece0 = tmp_path / "payload.part.piece00"
    piece0.write_bytes(payload[:8])
    validation = _DownloadValidation([_DownloadClient(payload)], len(payload))
    monkeypatch.setattr(
        PULLER,
        "remote_range_sha256",
        lambda _session, _path, _identity, start, end, _label:
            PULLER.hashlib.sha256(payload[start:end]).hexdigest(),
    )
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        proofs = PULLER.prove_existing_parallel_pieces(
            "/remote/payload",
            parent,
            "payload.part",
            len(payload),
            validation,
            _remote_download_identity(len(payload)),
        )
    finally:
        os.close(parent)

    assert len(proofs) == 1
    assert proofs[0].number == 0
    assert proofs[0].start == 0 and proofs[0].end == 8
    assert stat.S_IMODE(piece0.stat().st_mode) == 0o400


def test_install_no_replace_refuses_existing_final(tmp_path: Path) -> None:
    partial = tmp_path / "part"
    final = tmp_path / "final"
    partial.write_bytes(b"new")
    final.write_bytes(b"old")
    with pytest.raises(RuntimeError, match="refusing to replace"):
        PULLER.install_no_replace(partial, final)
    assert partial.read_bytes() == b"new"
    assert final.read_bytes() == b"old"


class _LocalSftp:
    def lstat(self, path: str):
        return Path(path).lstat()

    def mkdir(self, path: str, mode: int = 0o777) -> None:
        Path(path).mkdir(mode=mode)

    def open(self, path: str, mode: str):
        return Path(path).open("xb" if mode == "wx" else mode)

    def chmod(self, path: str, mode: int) -> None:
        Path(path).chmod(mode)

    def remove(self, path: str) -> None:
        Path(path).unlink()

    def listdir_attr(self, path: str):
        class Entry:
            def __init__(self, filename: str) -> None:
                self.filename = filename

        values = []
        for child in Path(path).iterdir():
            values.append(Entry(child.name))
        return values


class _DownloadAttributes:
    st_mode = stat.S_IFREG | 0o400
    st_mtime = 123
    st_ino = 456

    def __init__(self, size: int) -> None:
        self.st_size = size


class _DownloadSftp:
    def __init__(self, size: int) -> None:
        self.attributes = _DownloadAttributes(size)
        self.closed = False

    def lstat(self, _path: str):
        return self.attributes

    def close(self) -> None:
        self.closed = True


class _DownloadChannel:
    def __init__(self, payload: bytes, fail_after: int | None) -> None:
        self.payload = payload
        self.position = 0
        self.fail_after = fail_after
        self.timeout = None

    def settimeout(self, timeout: float) -> None:
        self.timeout = timeout

    def recv(self, size: int) -> bytes:
        if self.fail_after is not None and self.position >= self.fail_after:
            raise paramiko.SSHException(
                "Key-exchange timed out waiting for key negotiation"
            )
        limit = len(self.payload)
        if self.fail_after is not None:
            limit = min(limit, self.fail_after)
        end = min(limit, self.position + size)
        result = self.payload[self.position:end]
        self.position = end
        return result

    def recv_exit_status(self) -> int:
        return 0

    def exit_status_ready(self) -> bool:
        return True


class _CaptureChannel:
    def __init__(self, output: bytes) -> None:
        self.output = output
        self.position = 0
        self.timeout = None

    def settimeout(self, timeout: float) -> None:
        self.timeout = timeout

    def recv_ready(self) -> bool:
        return self.position < len(self.output)

    def recv(self, size: int) -> bytes:
        end = min(len(self.output), self.position + size)
        result = self.output[self.position:end]
        self.position = end
        return result

    def recv_stderr_ready(self) -> bool:
        return False

    def recv_stderr(self, _size: int) -> bytes:
        return b""

    def recv_exit_status(self) -> int:
        return 0

    def exit_status_ready(self) -> bool:
        return True


class _DownloadClient:
    def __init__(self, payload: bytes, fail_after: int | None = None) -> None:
        self.payload = payload
        self.fail_after = fail_after
        self.commands: list[str] = []
        self.closed = False

    def exec_command(self, command: str, timeout=None):
        del timeout
        self.commands.append(command)
        if "sha256sum" in command:
            match = PULLER.re.search(r"skip=(\d+) count=(\d+)", command)
            if match is None:
                raise AssertionError(f"malformed range-hash command: {command}")
            start = int(match.group(1))
            count = int(match.group(2))
            digest = PULLER.hashlib.sha256(
                self.payload[start:start + count]
            ).hexdigest().encode("ascii")
            stdout = type(
                "Stdout",
                (),
                {"channel": _CaptureChannel(digest + b"  -\n")},
            )()
            return None, stdout, None
        fields = {
            item.split("=", 1)[0]: item.split("=", 1)[1]
            for item in PULLER.shlex.split(command)
            if "=" in item
        }
        start = int(fields["skip"])
        count = int(fields["count"])
        relative_failure = self.fail_after
        channel = _DownloadChannel(
            self.payload[start:start + count],
            relative_failure,
        )
        stdout = type("Stdout", (), {"channel": channel})()
        return None, stdout, None

    def close(self) -> None:
        self.closed = True


class _DownloadValidation:
    def __init__(self, clients: list[_DownloadClient], payload_size: int) -> None:
        self.clients = clients
        self.payload_size = payload_size
        self.opens = 0

    def open_validated(self):
        if self.opens >= len(self.clients):
            raise AssertionError("unexpected additional reconnect")
        client = self.clients[self.opens]
        self.opens += 1
        return client, _DownloadSftp(self.payload_size)


def _remote_download_identity(size: int) -> tuple[int, int, int, int]:
    return PULLER._remote_identity(_DownloadAttributes(size))


def test_sequential_download_resumes_mid_file_rekey_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"0123456789abcdef"
    clients = [
        _DownloadClient(payload, fail_after=5),
        _DownloadClient(payload),
    ]
    validation = _DownloadValidation(clients, len(payload))
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", len(payload))
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)

    PULLER.download_range(
        object(),
        "/remote/payload",
        parent,
        "payload.part",
        0,
        len(payload),
        0,
        validation,
        _remote_download_identity(len(payload)),
    )

    assert (tmp_path / "payload.part").read_bytes() == payload
    assert validation.opens == 2
    assert "skip=0" in clients[0].commands[0]
    assert "skip=5" in clients[1].commands[0]


def test_download_rotates_transport_before_rekey_threshold(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"0123456789"
    clients = [_DownloadClient(payload) for _ in range(3)]
    validation = _DownloadValidation(clients, len(payload))
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", 4)
    evictions: list[tuple[int, int]] = []
    monkeypatch.setattr(
        PULLER,
        "evict_cached_range",
        lambda _descriptor, offset, length: evictions.append((offset, length)),
    )
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)

    PULLER.download_range(
        object(),
        "/remote/payload",
        parent,
        "payload.part",
        0,
        len(payload),
        0,
        validation,
        _remote_download_identity(len(payload)),
    )

    assert (tmp_path / "payload.part").read_bytes() == payload
    assert validation.opens == 3
    assert evictions == [
        (0, 4), (0, 4),
        (4, 4), (4, 4),
        (8, 2), (8, 2),
    ]
    assert [client.commands[0] for client in clients] == [
        "dd if=/remote/payload iflag=skip_bytes,count_bytes skip=0 count=4 status=none",
        "dd if=/remote/payload iflag=skip_bytes,count_bytes skip=4 count=4 status=none",
        "dd if=/remote/payload iflag=skip_bytes,count_bytes skip=8 count=2 status=none",
    ]


def test_window_hash_mismatch_retries_only_the_bad_window(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"abcdefgh"
    partial = tmp_path / "payload.part"
    partial.write_bytes(payload[:4])
    clients = [_DownloadClient(payload) for _ in range(5)]
    validation = _DownloadValidation(clients, len(payload))
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", 4)
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)
    monkeypatch.setattr(
        PULLER,
        "remote_range_sha256",
        lambda _session, _path, _identity, start, end, _label:
            PULLER.hashlib.sha256(payload[start:end]).hexdigest(),
    )
    real_append = PULLER._append_validated_window
    append_calls = 0

    def corrupt_first_window(
        session,
        remote_path,
        remote_identity,
        parent_descriptor,
        piece_name,
        range_start,
        target_size,
        label,
    ) -> None:
        nonlocal append_calls
        real_append(
            session,
            remote_path,
            remote_identity,
            parent_descriptor,
            piece_name,
            range_start,
            target_size,
            label,
        )
        append_calls += 1
        if append_calls == 1:
            descriptor = os.open(
                piece_name,
                os.O_RDWR | os.O_NOFOLLOW,
                dir_fd=parent_descriptor,
            )
            try:
                os.pwrite(descriptor, b"X", target_size - 1)
                os.fsync(descriptor)
            finally:
                os.close(descriptor)

    monkeypatch.setattr(PULLER, "_append_validated_window", corrupt_first_window)
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    PULLER.download_range(
        object(),
        "/remote/payload",
        parent,
        partial.name,
        0,
        len(payload),
        0,
        validation,
        _remote_download_identity(len(payload)),
    )

    assert partial.read_bytes() == payload
    assert append_calls == 2
    assert validation.opens == 5
    assert "skip=4 count=4" in clients[1].commands[0]
    assert "skip=4 count=4" in clients[4].commands[0]


def test_window_hash_retry_limit_preserves_last_bad_window(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"abcd"
    partial = tmp_path / "payload.part"
    clients = [_DownloadClient(payload) for _ in range(9)]
    validation = _DownloadValidation(clients, len(payload))
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", len(payload))
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)
    monkeypatch.setattr(
        PULLER,
        "remote_range_sha256",
        lambda _session, _path, _identity, start, end, _label:
            PULLER.hashlib.sha256(payload[start:end]).hexdigest(),
    )
    real_append = PULLER._append_validated_window
    append_calls = 0

    def corrupt_every_window(
        session,
        remote_path,
        remote_identity,
        parent_descriptor,
        piece_name,
        range_start,
        target_size,
        label,
    ) -> None:
        nonlocal append_calls
        real_append(
            session,
            remote_path,
            remote_identity,
            parent_descriptor,
            piece_name,
            range_start,
            target_size,
            label,
        )
        append_calls += 1
        descriptor = os.open(
            piece_name,
            os.O_RDWR | os.O_NOFOLLOW,
            dir_fd=parent_descriptor,
        )
        try:
            os.pwrite(descriptor, b"X", target_size - 1)
            os.fsync(descriptor)
        finally:
            os.close(descriptor)

    monkeypatch.setattr(PULLER, "_append_validated_window", corrupt_every_window)
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    with pytest.raises(PULLER.WindowHashMismatch, match="retry limit 3 exhausted"):
        PULLER.download_range(
            object(),
            "/remote/payload",
            parent,
            partial.name,
            0,
            len(payload),
            0,
            validation,
            _remote_download_identity(len(payload)),
        )

    assert append_calls == PULLER.WINDOW_DOWNLOAD_ATTEMPTS
    assert validation.opens == 9
    assert partial.stat().st_size == len(payload)
    assert partial.read_bytes() != payload


def test_resumed_rate_counts_only_bytes_downloaded_this_run(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    window = 1024 * 1024
    payload = b"a" * window + b"b" * window
    partial = tmp_path / "payload.part"
    partial.write_bytes(payload[:window])
    clients = [_DownloadClient(payload) for _ in range(2)]
    validation = _DownloadValidation(clients, len(payload))
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", window)
    monkeypatch.setattr(
        PULLER,
        "remote_range_sha256",
        lambda _session, _path, _identity, start, end, _label:
            PULLER.hashlib.sha256(payload[start:end]).hexdigest(),
    )
    clock = iter((100.0, 100.0, 102.0, 103.0))
    monkeypatch.setattr(PULLER.time, "monotonic", lambda: next(clock))
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    PULLER.download_range(
        object(),
        "/remote/payload",
        parent,
        partial.name,
        0,
        len(payload),
        0,
        validation,
        _remote_download_identity(len(payload)),
    )

    assert partial.read_bytes() == payload
    assert "0.50 MiB/s" in capsys.readouterr().out


def test_data_channel_stall_uses_timeout_and_bounded_reconnects(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"stalled-channel"

    class StalledChannel(_DownloadChannel):
        def recv(self, _size: int) -> bytes:
            raise PULLER.socket.timeout("data channel stalled")

    class StalledClient(_DownloadClient):
        def exec_command(self, command: str, timeout=None):
            del timeout
            self.commands.append(command)
            channel = StalledChannel(payload, None)
            self.channels.append(channel)
            stdout = type("Stdout", (), {"channel": channel})()
            return None, stdout, None

        def __init__(self) -> None:
            super().__init__(payload)
            self.channels: list[StalledChannel] = []

    clients = [
        StalledClient() for _ in range(PULLER.STREAM_RECONNECT_ATTEMPTS + 1)
    ]
    validation = _DownloadValidation(clients, len(payload))
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", len(payload))
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)

    with pytest.raises(RuntimeError, match="reconnect limit .* exhausted"):
        PULLER.download_range(
            object(),
            "/remote/payload",
            parent,
            "payload.part",
            0,
            len(payload),
            0,
            validation,
            _remote_download_identity(len(payload)),
        )

    assert validation.opens == PULLER.STREAM_RECONNECT_ATTEMPTS + 1
    assert all(
        client.channels[0].timeout == PULLER.CHANNEL_IO_TIMEOUT_SECONDS
        for client in clients
    )
    assert (tmp_path / "payload.part").read_bytes() == b""


def test_remote_command_exit_status_has_hard_deadline(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class Channel:
        timeout = None

        def settimeout(self, value: float) -> None:
            self.timeout = value

        def exit_status_ready(self) -> bool:
            return False

        def recv_stderr_ready(self) -> bool:
            return False

        def recv_ready(self) -> bool:
            return False

    channel = Channel()

    class Client:
        def exec_command(self, _command: str, timeout=None):
            assert timeout == PULLER.CHANNEL_IO_TIMEOUT_SECONDS
            stdout = type("Stdout", (), {"channel": channel})()
            return None, stdout, None

    clock = iter((0.0, PULLER.CHANNEL_IO_TIMEOUT_SECONDS + 1.0))
    monkeypatch.setattr(PULLER.time, "monotonic", lambda: next(clock))
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)

    with pytest.raises(PULLER.RetryableRemoteTransport, match="timed out"):
        PULLER.exec_checked(Client(), ["/usr/bin/true"])

    assert channel.timeout == PULLER.CHANNEL_IO_TIMEOUT_SECONDS


def test_parallel_range_reconnect_limit_is_bounded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"abcdefgh"
    clients = [
        _DownloadClient(payload, fail_after=0)
        for _ in range(PULLER.STREAM_RECONNECT_ATTEMPTS + 1)
    ]
    validation = _DownloadValidation(clients, len(payload))
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", len(payload))
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)

    with pytest.raises(RuntimeError, match="reconnect limit .* exhausted"):
        PULLER.download_range(
            object(),
            "/remote/payload",
            parent,
            "range.part",
            0,
            len(payload),
            1,
            validation,
            _remote_download_identity(len(payload)),
        )

    assert validation.opens == PULLER.STREAM_RECONNECT_ATTEMPTS + 1
    assert (tmp_path / "range.part").read_bytes() == b""
    assert not list(tmp_path.glob("*.ack"))


def test_parallel_download_resumes_its_piece_after_rekey(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"parallel-payload"
    clients = [
        _DownloadClient(payload, fail_after=4),
        _DownloadClient(payload),
        _DownloadClient(payload),
    ]
    validation = _DownloadValidation(clients, len(payload))
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", len(payload))
    monkeypatch.setattr(PULLER, "PARALLEL_STREAMS", 1)
    monkeypatch.setattr(
        PULLER,
        "remote_range_sha256",
        lambda _session, _path, _identity, start, end, _label:
            PULLER.hashlib.sha256(payload[start:end]).hexdigest(),
    )
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        PULLER.parallel_download(
            object(),
            "/remote/payload",
            parent,
            "payload.part",
            len(payload),
            validation,
            _remote_download_identity(len(payload)),
        )
    finally:
        os.close(parent)

    assert (tmp_path / "payload.part").read_bytes() == payload
    assert validation.opens == 3
    assert "skip=0" in clients[0].commands[0]
    assert "skip=4" in clients[1].commands[0]
    assert (tmp_path / "payload.part.piece00").read_bytes() == payload


def test_local_fsync_failure_is_not_misclassified_as_remote_rekey(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"local-storage-failure"
    validation = _DownloadValidation([_DownloadClient(payload)], len(payload))
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", len(payload))
    monkeypatch.setattr(
        PULLER.os,
        "fsync",
        lambda _descriptor: (_ for _ in ()).throw(OSError(28, "No space left")),
    )
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)

    with pytest.raises(OSError, match="No space left"):
        PULLER.download_range(
            object(),
            "/remote/payload",
            parent,
            "payload.part",
            0,
            len(payload),
            0,
            validation,
            _remote_download_identity(len(payload)),
        )

    assert validation.opens == 1
    assert (tmp_path / "payload.part").read_bytes() == payload
    assert not list(tmp_path.glob("*.ack"))


def test_rekey_cleanup_fsync_failure_stops_without_reconnect(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"remote-and-local-failure"
    validation = _DownloadValidation(
        [_DownloadClient(payload, fail_after=3), _DownloadClient(payload)],
        len(payload),
    )
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", len(payload))
    monkeypatch.setattr(
        PULLER.os,
        "fsync",
        lambda _descriptor: (_ for _ in ()).throw(OSError(5, "I/O error")),
    )
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)

    with pytest.raises(OSError, match="I/O error") as failure:
        PULLER.download_range(
            object(),
            "/remote/payload",
            parent,
            "payload.part",
            0,
            len(payload),
            0,
            validation,
            _remote_download_identity(len(payload)),
        )

    assert isinstance(failure.value.__cause__, PULLER.RetryableRemoteTransport)
    assert validation.opens == 1
    assert (tmp_path / "payload.part").read_bytes() == payload[:3]
    assert not list(tmp_path.glob("*.ack"))


def test_partial_writer_handles_short_local_writes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"short-write-payload"
    validation = _DownloadValidation([_DownloadClient(payload)], len(payload))
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", len(payload))
    real_open = PULLER.open_append_nofollow_at
    write_calls: list[int] = []

    class ShortWriter:
        def __init__(self, stream) -> None:
            self.stream = stream

        def write(self, data) -> int:
            written = self.stream.write(data[:2])
            write_calls.append(written)
            return written

        def flush(self) -> None:
            self.stream.flush()

        def fileno(self) -> int:
            return self.stream.fileno()

        def close(self) -> None:
            self.stream.close()

    def short_open(parent_descriptor: int, name: str):
        stream, size = real_open(parent_descriptor, name)
        return ShortWriter(stream), size

    monkeypatch.setattr(PULLER, "open_append_nofollow_at", short_open)
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    PULLER.download_range(
        object(),
        "/remote/payload",
        parent,
        "payload.part",
        0,
        len(payload),
        0,
        validation,
        _remote_download_identity(len(payload)),
    )

    assert (tmp_path / "payload.part").read_bytes() == payload
    assert len(write_calls) > 1
    assert all(0 < written <= 2 for written in write_calls)


def test_reconnect_refuses_changed_pinned_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"abcdefgh"
    clients = [_DownloadClient(payload, fail_after=3), _DownloadClient(payload)]
    sftps = [_DownloadSftp(len(payload)), _DownloadSftp(len(payload))]
    manifest_reads = iter((b"manifest-v1", b"manifest-v2"))

    def fake_connect(_config):
        index = fake_connect.calls
        fake_connect.calls += 1
        return clients[index], sftps[index]

    fake_connect.calls = 0
    monkeypatch.setattr(PULLER, "connect", fake_connect)
    monkeypatch.setattr(
        PULLER,
        "read_remote_bytes",
        lambda *_args, **_kwargs: next(manifest_reads),
    )
    monkeypatch.setattr(PULLER, "verify_mount", lambda *_args: {"source": "nas"})
    monkeypatch.setattr(PULLER, "write_probe", lambda *_args: None)
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", len(payload))

    class Destination:
        def assert_visible(self) -> None:
            return None

    validation = PULLER.ReconnectValidation(
        config=object(),
        remote_manifest="/remote/segment.manifest.ready",
        manifest_bytes=b"manifest-v1",
        destination=Destination(),
        destination_path=tmp_path,
        expected_mount_source="nas",
        expected_mount_fstype="nfs",
        mount_snapshot={"source": "nas"},
    )
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    with pytest.raises(PULLER.TrustBoundaryViolation, match="manifest changed"):
        PULLER.download_range(
            object(),
            "/remote/payload",
            parent,
            "payload.part",
            0,
            len(payload),
            0,
            validation,
            _remote_download_identity(len(payload)),
        )

    assert (tmp_path / "payload.part").read_bytes() == payload[:3]
    assert fake_connect.calls == 2
    assert not list(tmp_path.glob("*.ack"))


def test_reconnect_refuses_changed_shared_source_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = b"abcdefgh"
    clients = [_DownloadClient(payload, fail_after=3), _DownloadClient(payload)]
    sftps = [_DownloadSftp(len(payload)), _DownloadSftp(len(payload))]
    sftps[1].attributes.st_mtime = 124

    def fake_connect(_config):
        index = fake_connect.calls
        fake_connect.calls += 1
        return clients[index], sftps[index]

    fake_connect.calls = 0
    monkeypatch.setattr(PULLER, "connect", fake_connect)
    monkeypatch.setattr(PULLER, "read_remote_bytes", lambda *_args, **_kwargs: b"v1")
    monkeypatch.setattr(PULLER, "verify_mount", lambda *_args: {"source": "nas"})
    monkeypatch.setattr(PULLER, "write_probe", lambda *_args: None)
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", len(payload))

    class Destination:
        def assert_visible(self) -> None:
            return None

    validation = PULLER.ReconnectValidation(
        config=object(),
        remote_manifest="/remote/segment.manifest.ready",
        manifest_bytes=b"v1",
        destination=Destination(),
        destination_path=tmp_path,
        expected_mount_source="nas",
        expected_mount_fstype="nfs",
        mount_snapshot={"source": "nas"},
    )
    parent = os.open(tmp_path, os.O_RDONLY | os.O_DIRECTORY)
    with pytest.raises(PULLER.TrustBoundaryViolation, match="source identity changed"):
        PULLER.download_range(
            object(),
            "/remote/payload",
            parent,
            "payload.part",
            0,
            len(payload),
            0,
            validation,
            _remote_download_identity(len(payload)),
        )

    assert (tmp_path / "payload.part").read_bytes() == payload[:3]
    assert fake_connect.calls == 2
    assert not list(tmp_path.glob("*.ack"))


def test_remote_ack_is_create_no_replace_and_read_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def local_exec(_client: object, argv: list[str]) -> None:
        assert argv[:2] == ["/usr/bin/ln", "--"]
        os.link(argv[2], argv[3])

    monkeypatch.setattr(PULLER, "exec_checked", local_exec)
    monkeypatch.setattr(PULLER, "remote_fsync", lambda *_args: None)
    ack = tmp_path / "remote" / "acks" / "segment.ack"
    payload = b'{"schema":1}\n'
    sftp = _LocalSftp()

    assert PULLER.remote_immutable_write(object(), sftp, str(ack), payload)
    original = ack.lstat()
    assert stat.S_IMODE(original.st_mode) == 0o400
    assert original.st_nlink == 1
    assert not PULLER.remote_immutable_write(object(), sftp, str(ack), payload)
    assert ack.lstat().st_ino == original.st_ino
    with pytest.raises(RuntimeError, match="refusing to replace"):
        PULLER.remote_immutable_write(object(), sftp, str(ack), b"different\n")


def test_remote_ack_recovers_linked_temporary_alias(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(PULLER, "remote_fsync", lambda *_args: None)
    ack = tmp_path / "remote" / "acks" / "segment.ack"
    ack.parent.mkdir(parents=True)
    payload = b"remote-immutable\n"
    temporary = ack.parent / f".{ack.name}.123.deadbeef.tmp"
    temporary.write_bytes(payload)
    temporary.chmod(0o400)
    os.link(temporary, ack)
    assert ack.lstat().st_nlink == 2

    assert not PULLER.remote_immutable_write(
        object(),
        _LocalSftp(),
        str(ack),
        payload,
    )
    assert ack.lstat().st_nlink == 1
    assert not temporary.exists()


def test_remote_ack_recovers_when_link_succeeds_but_response_is_lost(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    ack = tmp_path / "remote" / "acks" / "segment.ack"
    ack.parent.mkdir(parents=True)
    payload = b'{"authorizes_remote_deletion":false}\n'
    calls = 0

    def link_then_drop(_client: object, argv: list[str]) -> None:
        nonlocal calls
        calls += 1
        os.link(argv[2], argv[3])
        raise paramiko.SSHException("response lost after successful ln")

    class AckValidation:
        def __init__(self) -> None:
            self.opens = 0

        def open_validated(self):
            self.opens += 1
            return object(), _LocalSftp()

    validation = AckValidation()
    monkeypatch.setattr(PULLER, "exec_checked", link_then_drop)
    monkeypatch.setattr(PULLER, "remote_fsync", lambda *_args: None)
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)

    with PULLER.ValidatedRemoteSession(validation) as session:
        created = session.remote_call(
            "publish ACK",
            lambda client, sftp: PULLER.retryable_remote_io(
                lambda: PULLER.remote_immutable_write(
                    client,
                    sftp,
                    str(ack),
                    payload,
                )
            ),
        )

    assert created is False
    assert validation.opens == 2
    assert calls == 1
    assert ack.read_bytes() == payload
    assert stat.S_IMODE(ack.stat().st_mode) == 0o400
    assert ack.stat().st_nlink == 1


def test_remote_ack_reconnects_when_link_fails_before_creation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    ack = tmp_path / "remote" / "acks" / "segment.ack"
    ack.parent.mkdir(parents=True)
    payload = b'{"authorizes_remote_deletion":false}\n'
    calls = 0

    def drop_then_link(_client: object, argv: list[str]) -> None:
        nonlocal calls
        calls += 1
        if calls == 1:
            raise PULLER.RetryableRemoteTransport("connection lost before ln")
        os.link(argv[2], argv[3])

    class AckValidation:
        def __init__(self) -> None:
            self.opens = 0

        def open_validated(self):
            self.opens += 1
            return object(), _LocalSftp()

    validation = AckValidation()
    monkeypatch.setattr(PULLER, "exec_checked", drop_then_link)
    monkeypatch.setattr(PULLER, "remote_fsync", lambda *_args: None)
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)

    with PULLER.ValidatedRemoteSession(validation) as session:
        created = session.remote_call(
            "publish ACK",
            lambda client, sftp: PULLER.retryable_remote_io(
                lambda: PULLER.remote_immutable_write(
                    client,
                    sftp,
                    str(ack),
                    payload,
                )
            ),
        )

    assert created is True
    assert validation.opens == 2
    assert calls == 2
    assert ack.read_bytes() == payload
    assert ack.stat().st_nlink == 1


def test_ack_reconnect_manifest_change_fails_closed_after_lost_link_response(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    ack = tmp_path / "remote" / "acks" / "segment.ack"
    ack.parent.mkdir(parents=True)
    payload = b'{"authorizes_remote_deletion":false}\n'
    sftps = [_LocalSftp(), _LocalSftp()]
    clients = [object(), object()]
    original_read_remote_bytes = PULLER.read_remote_bytes
    opens = 0
    manifest_reads = iter((b"manifest-v1", b"manifest-v2"))
    link_calls = 0

    def fake_connect(_config):
        nonlocal opens
        result = clients[opens], sftps[opens]
        opens += 1
        return result

    def routed_remote_read(sftp, path: str, **kwargs):
        if path == "/remote/segment.manifest.ready":
            return next(manifest_reads)
        return original_read_remote_bytes(sftp, path, **kwargs)

    def link_then_drop(_client: object, argv: list[str]) -> None:
        nonlocal link_calls
        link_calls += 1
        os.link(argv[2], argv[3])
        raise paramiko.SSHException("response lost after successful ln")

    class Destination:
        def assert_visible(self) -> None:
            return None

    monkeypatch.setattr(PULLER, "connect", fake_connect)
    monkeypatch.setattr(PULLER, "read_remote_bytes", routed_remote_read)
    monkeypatch.setattr(PULLER, "verify_mount", lambda *_args: {"source": "nas"})
    monkeypatch.setattr(PULLER, "write_probe", lambda *_args: None)
    monkeypatch.setattr(PULLER, "exec_checked", link_then_drop)
    monkeypatch.setattr(PULLER, "remote_fsync", lambda *_args: None)
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)
    validation = PULLER.ReconnectValidation(
        config=object(),
        remote_manifest="/remote/segment.manifest.ready",
        manifest_bytes=b"manifest-v1",
        destination=Destination(),
        destination_path=tmp_path,
        expected_mount_source="nas",
        expected_mount_fstype="nfs",
        mount_snapshot={"source": "nas"},
    )

    with pytest.raises(PULLER.TrustBoundaryViolation, match="manifest changed"):
        with PULLER.ValidatedRemoteSession(validation) as session:
            session.remote_call(
                "publish ACK",
                lambda client, sftp: PULLER.retryable_remote_io(
                    lambda: PULLER.remote_immutable_write(
                        client,
                        sftp,
                        str(ack),
                        payload,
                    )
                ),
            )

    assert opens == 2
    assert link_calls == 1
    assert ack.read_bytes() == payload
    assert json.loads(payload)["authorizes_remote_deletion"] is False


def test_ack_publish_reconnect_limit_is_bounded_and_no_ack_commits(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    class AckValidation:
        def __init__(self) -> None:
            self.opens = 0

        def open_validated(self):
            self.opens += 1
            return object(), _LocalSftp()

    validation = AckValidation()
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)
    ack = tmp_path / "remote" / "acks" / "segment.ack"
    attempts = 0

    def always_drop(_client, _sftp):
        nonlocal attempts
        attempts += 1
        raise paramiko.SSHException("ACK channel stalled")

    with pytest.raises(RuntimeError, match="reconnect limit .* exhausted"):
        with PULLER.ValidatedRemoteSession(validation) as session:
            session.remote_call("publish ACK", always_drop)

    assert validation.opens == PULLER.STREAM_RECONNECT_ATTEMPTS + 1
    assert attempts == PULLER.STREAM_RECONNECT_ATTEMPTS + 1
    assert not ack.exists()
