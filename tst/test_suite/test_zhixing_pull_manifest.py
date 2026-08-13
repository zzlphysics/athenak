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
    with pytest.raises(RuntimeError, match="final NAS verification failed"):
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
    piece0 = tmp_path / f"{partial_name}.piece00"
    piece1 = tmp_path / f"{partial_name}.piece01"
    piece0.write_bytes(b"a" * PULLER.CHUNK)
    piece1.write_bytes(b"b")
    (tmp_path / f"{partial_name}.assembling").write_bytes(b"stale")
    monkeypatch.setattr(PULLER, "PARALLEL_STREAMS", 2)
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
            raise paramiko.SSHException("Key-exchange timed out waiting for key negotiation")
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


class _DownloadClient:
    def __init__(self, payload: bytes, fail_after: int | None = None) -> None:
        self.payload = payload
        self.fail_after = fail_after
        self.commands: list[str] = []
        self.closed = False

    def exec_command(self, command: str, timeout=None):
        del timeout
        self.commands.append(command)
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
    assert [client.commands[0] for client in clients] == [
        "dd if=/remote/payload iflag=skip_bytes,count_bytes skip=0 count=4 status=none",
        "dd if=/remote/payload iflag=skip_bytes,count_bytes skip=4 count=4 status=none",
        "dd if=/remote/payload iflag=skip_bytes,count_bytes skip=8 count=2 status=none",
    ]


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
    ]
    validation = _DownloadValidation(clients, len(payload))
    monkeypatch.setattr(PULLER.time, "sleep", lambda _delay: None)
    monkeypatch.setattr(PULLER, "STREAM_WINDOW_BYTES", len(payload))
    monkeypatch.setattr(PULLER, "PARALLEL_STREAMS", 1)
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
    assert validation.opens == 2
    assert "skip=0" in clients[0].commands[0]
    assert "skip=4" in clients[1].commands[0]
    assert not list(tmp_path.glob("*.piece*"))


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
