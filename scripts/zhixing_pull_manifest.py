#!/usr/bin/env python3
"""Resume and verify one immutable AthenaK manifest over Zhixing SFTP.

This program deliberately stops at a durable local+remote acknowledgement.  It never
deletes cloud data, and its ACK never authorizes deletion.  The current independent
cleanup gate fails closed; a future storage-barrier-backed gate and transactional
consumer would have to revalidate all evidence before removing any source file.
"""

from __future__ import annotations

import argparse
import base64
import concurrent.futures
import contextlib
from dataclasses import dataclass
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import posixpath
import secrets
import shlex
import shutil
import socket
import stat
import subprocess
import threading
import time
import urllib.parse

import paramiko


DEFAULT_SSH_STATE = Path("/tmp/zhixing_l4_ssh.private.json")
CHUNK = 8 * 1024 * 1024
PROGRESS_INTERVAL = 256 * 1024 * 1024
FSYNC_INTERVAL = 512 * 1024 * 1024
MIN_RESERVE = 50 * 1024 * 1024 * 1024
PARALLEL_THRESHOLD = 512 * 1024 * 1024
PARALLEL_STREAMS = 2
CONNECT_ATTEMPTS = 5
CONNECT_BACKOFF_SECONDS = 5.0
STREAM_RECONNECT_ATTEMPTS = 4
STREAM_RECONNECT_BACKOFF_SECONDS = 2.0
STREAM_WINDOW_BYTES = 256 * 1024 * 1024
CHANNEL_IO_TIMEOUT_SECONDS = 120.0
CHANNEL_STATUS_POLL_SECONDS = 0.1
RESERVED_DESTINATION_COMPONENTS = frozenset({
    ".acks", ".incoming", ".locks", ".manifests",
})
PRINT_LOCK = threading.Lock()


@dataclass(frozen=True)
class ConnectionConfig:
    host: str
    port: int
    username: str
    password: str
    host_key_sha256: str


def key_sha256(key: paramiko.PKey) -> str:
    encoded = base64.b64encode(hashlib.sha256(key.asbytes()).digest()).decode("ascii")
    return f"SHA256:{encoded.rstrip('=')}"


def normalize_host_key_sha256(value: str) -> str:
    candidate = value.strip()
    if candidate.startswith("SHA256:"):
        candidate = candidate[7:]
    if not candidate or any(char.isspace() for char in candidate):
        raise ValueError("SSH host-key fingerprint must be one SHA256 fingerprint")
    padded = candidate + "=" * (-len(candidate) % 4)
    try:
        decoded = base64.b64decode(padded, validate=True)
    except (ValueError, TypeError) as exc:
        raise ValueError("invalid base64 SSH host-key fingerprint") from exc
    if len(decoded) != hashlib.sha256().digest_size:
        raise ValueError("SSH host-key fingerprint must contain 32 SHA256 bytes")
    return f"SHA256:{base64.b64encode(decoded).decode('ascii').rstrip('=')}"


class PinnedHostKeyPolicy(paramiko.MissingHostKeyPolicy):
    """Accept exactly one out-of-band SHA256 host-key fingerprint."""

    def __init__(self, expected: str) -> None:
        self.expected = normalize_host_key_sha256(expected)

    def missing_host_key(
        self,
        client: paramiko.SSHClient,
        hostname: str,
        key: paramiko.PKey,
    ) -> None:
        del client
        actual = key_sha256(key)
        if actual != self.expected:
            raise paramiko.SSHException(
                f"SSH host key mismatch for {hostname}: {actual} != {self.expected}"
            )


def _identity(info: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_nlink,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )


def _remote_identity(info: paramiko.SFTPAttributes) -> tuple[int, int, int, int]:
    return (
        int(info.st_mode or 0),
        int(info.st_size or 0),
        int(info.st_mtime or 0),
        int(getattr(info, "st_ino", 0) or 0),
    )


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def descriptor_mount_id(descriptor: int) -> int:
    """Read the Linux mount ID for an open descriptor, failing closed."""

    try:
        lines = Path(f"/proc/self/fdinfo/{descriptor}").read_text(
            encoding="ascii", errors="strict"
        ).splitlines()
    except (OSError, UnicodeError) as exc:
        raise RuntimeError(
            f"cannot bind destination mount identity for fd {descriptor}: {exc}"
        ) from exc
    matches = [
        line.split(":", 1)[1].strip()
        for line in lines
        if line.startswith("mnt_id:")
    ]
    if len(matches) != 1 or not matches[0].isdigit():
        raise RuntimeError(f"cannot parse mount ID for fd {descriptor}")
    return int(matches[0])


def safe_directory_component(component: str) -> None:
    if component in ("", ".", "..") or "/" in component or "\x00" in component:
        raise RuntimeError(f"unsafe destination directory component: {component!r}")


class DestinationRoot:
    """Pinned NAS root used for every local pathname operation.

    Children are opened relative to held directory descriptors with `O_NOFOLLOW`.
    Mount IDs, rather than only `st_dev`, reject nested bind mounts, including binds
    of the same filesystem.
    """

    def __init__(self, path: Path) -> None:
        self.path = path
        self.descriptor = -1
        self.identity: tuple[int, int] | None = None
        self.mount_id: int | None = None

    def __enter__(self) -> "DestinationRoot":
        path_info = self.path.lstat()
        if not stat.S_ISDIR(path_info.st_mode):
            raise RuntimeError(f"destination root is not a directory: {self.path}")
        self.descriptor = os.open(
            self.path,
            os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC,
        )
        try:
            descriptor_info = os.fstat(self.descriptor)
            self.identity = (descriptor_info.st_dev, descriptor_info.st_ino)
            self.mount_id = descriptor_mount_id(self.descriptor)
            if self.identity != (path_info.st_dev, path_info.st_ino):
                raise RuntimeError(
                    f"destination root changed while opening: {self.path}"
                )
        except BaseException:
            self.close()
            raise
        return self

    def close(self) -> None:
        if self.descriptor >= 0:
            os.close(self.descriptor)
            self.descriptor = -1

    def __exit__(self, *unused: object) -> None:
        self.close()

    def assert_visible(self) -> None:
        if self.descriptor < 0 or self.identity is None or self.mount_id is None:
            raise RuntimeError("destination root is not open")
        descriptor_info = os.fstat(self.descriptor)
        try:
            path_info = self.path.lstat()
        except OSError as exc:
            raise RuntimeError(
                f"destination root pathname is unavailable: {self.path}: {exc}"
            ) from exc
        if (
            not stat.S_ISDIR(path_info.st_mode)
            or (descriptor_info.st_dev, descriptor_info.st_ino) != self.identity
            or (path_info.st_dev, path_info.st_ino) != self.identity
            or descriptor_mount_id(self.descriptor) != self.mount_id
        ):
            raise RuntimeError(f"destination root identity changed: {self.path}")

    def open_directory(
        self,
        components: tuple[str, ...] = (),
        *,
        create: bool = False,
    ) -> int:
        self.assert_visible()
        assert self.mount_id is not None
        current = os.dup(self.descriptor)
        flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
        try:
            for component in components:
                safe_directory_component(component)
                try:
                    child = os.open(component, flags, dir_fd=current)
                except FileNotFoundError:
                    if not create:
                        raise
                    os.mkdir(component, mode=0o700, dir_fd=current)
                    os.fsync(current)
                    child = os.open(component, flags, dir_fd=current)
                if descriptor_mount_id(child) != self.mount_id:
                    os.close(child)
                    raise RuntimeError(
                        "destination contains a nested mount at component "
                        f"{component!r} below {self.path}"
                    )
                os.close(current)
                current = child
            return current
        except BaseException:
            os.close(current)
            raise

    def assert_directory_visible(
        self,
        components: tuple[str, ...],
        expected_descriptor: int,
    ) -> None:
        expected = os.fstat(expected_descriptor)
        expected_mount = descriptor_mount_id(expected_descriptor)
        fresh = self.open_directory(components)
        try:
            actual = os.fstat(fresh)
            if (
                (actual.st_dev, actual.st_ino) != (expected.st_dev, expected.st_ino)
                or descriptor_mount_id(fresh) != expected_mount
            ):
                joined = "/".join(components)
                raise RuntimeError(
                    f"destination directory pathname changed: {joined}"
                )
        finally:
            os.close(fresh)

    def open_relative_parent(
        self,
        relative: PurePosixPath,
        *,
        prefix: tuple[str, ...] = (),
        create: bool = False,
    ) -> tuple[int, str]:
        return (
            self.open_directory(
                prefix + tuple(relative.parts[:-1]),
                create=create,
            ),
            relative.name,
        )


def secure_directory(
    root: Path,
    components: tuple[str, ...] = (),
    *,
    create: bool = False,
) -> None:
    """Compatibility wrapper; critical operations retain their returned dirfd."""

    with DestinationRoot(root) as destination:
        descriptor = destination.open_directory(components, create=create)
        os.close(descriptor)


def secure_relative_parent(
    destination: Path,
    relative: PurePosixPath,
    *,
    prefix: tuple[str, ...] = (),
    create: bool = False,
) -> None:
    with DestinationRoot(destination) as root:
        descriptor, _ = root.open_relative_parent(
            relative,
            prefix=prefix,
            create=create,
        )
        os.close(descriptor)


def read_regular_bytes(
    path: Path,
    *,
    require_read_only: bool = False,
    require_single_link: bool = True,
) -> bytes:
    before_path = path.lstat()
    if not stat.S_ISREG(before_path.st_mode):
        raise RuntimeError(f"expected regular file: {path}")
    if require_single_link and before_path.st_nlink != 1:
        raise RuntimeError(f"expected one-link regular file: {path}")
    if require_read_only and before_path.st_mode & 0o222:
        raise RuntimeError(f"immutable record is writable: {path}")
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        before = os.fstat(descriptor)
        if _identity(before) != _identity(before_path):
            raise RuntimeError(f"path changed while opening: {path}")
        chunks: list[bytes] = []
        while True:
            chunk = os.read(descriptor, CHUNK)
            if not chunk:
                break
            chunks.append(chunk)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    after_path = path.lstat()
    if _identity(before) != _identity(after) or _identity(after) != _identity(after_path):
        raise RuntimeError(f"file changed while reading: {path}")
    return b"".join(chunks)


def stat_regular_at(
    parent_descriptor: int,
    name: str,
    *,
    require_single_link: bool = True,
    require_read_only: bool = False,
) -> os.stat_result:
    safe_directory_component(name)
    info = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    if not stat.S_ISREG(info.st_mode):
        raise RuntimeError(f"expected regular file: {name}")
    if require_single_link and info.st_nlink != 1:
        raise RuntimeError(f"expected one-link regular file: {name}")
    if require_read_only and info.st_mode & 0o222:
        raise RuntimeError(f"immutable record is writable: {name}")
    return info


def read_regular_bytes_at(
    parent_descriptor: int,
    name: str,
    *,
    require_read_only: bool = False,
    require_single_link: bool = True,
) -> bytes:
    before_path = stat_regular_at(
        parent_descriptor,
        name,
        require_single_link=require_single_link,
        require_read_only=require_read_only,
    )
    descriptor = os.open(
        name,
        os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC,
        dir_fd=parent_descriptor,
    )
    try:
        before = os.fstat(descriptor)
        if _identity(before) != _identity(before_path):
            raise RuntimeError(f"file changed while opening: {name}")
        chunks: list[bytes] = []
        while True:
            chunk = os.read(descriptor, CHUNK)
            if not chunk:
                break
            chunks.append(chunk)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    after_path = stat_regular_at(
        parent_descriptor,
        name,
        require_single_link=require_single_link,
        require_read_only=require_read_only,
    )
    if _identity(before) != _identity(after) or _identity(after) != _identity(after_path):
        raise RuntimeError(f"file changed while reading: {name}")
    return b"".join(chunks)


def sha256_regular_at(
    parent_descriptor: int,
    name: str,
    *,
    sync: bool = False,
    seal_read_only: bool = False,
    require_single_link: bool = True,
) -> tuple[os.stat_result, str]:
    before_path = stat_regular_at(
        parent_descriptor,
        name,
        require_single_link=require_single_link,
    )
    descriptor = os.open(
        name,
        os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC,
        dir_fd=parent_descriptor,
    )
    digest = hashlib.sha256()
    try:
        before = os.fstat(descriptor)
        if _identity(before) != _identity(before_path):
            raise RuntimeError(f"file changed while opening: {name}")
        if seal_read_only and before.st_mode & 0o222:
            os.fchmod(descriptor, stat.S_IMODE(before.st_mode) & ~0o222)
            os.fsync(descriptor)
            before = os.fstat(descriptor)
            sealed_path = stat_regular_at(
                parent_descriptor,
                name,
                require_single_link=require_single_link,
                require_read_only=True,
            )
            if _identity(before) != _identity(sealed_path):
                raise RuntimeError(f"file changed while sealing read-only: {name}")
        while True:
            chunk = os.read(descriptor, CHUNK)
            if not chunk:
                break
            digest.update(chunk)
        if sync:
            os.fsync(descriptor)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    after_path = stat_regular_at(
        parent_descriptor,
        name,
        require_single_link=require_single_link,
        require_read_only=seal_read_only,
    )
    if _identity(before) != _identity(after) or _identity(after) != _identity(after_path):
        raise RuntimeError(f"file changed while hashing: {name}")
    return after, digest.hexdigest()


def recover_local_ready_aliases_at(
    parent_descriptor: int,
    name: str,
    data: bytes,
) -> bool:
    try:
        target = stat_regular_at(
            parent_descriptor,
            name,
            require_single_link=False,
            require_read_only=True,
        )
    except FileNotFoundError:
        return False
    if read_regular_bytes_at(
        parent_descriptor,
        name,
        require_read_only=True,
        require_single_link=False,
    ) != data:
        raise RuntimeError(f"refusing changed immutable record: {name}")
    aliases: list[str] = []
    if target.st_nlink > 1:
        prefix = f".{name}."
        for candidate in os.listdir(parent_descriptor):
            if not candidate.startswith(prefix) or not candidate.endswith(".tmp"):
                continue
            candidate_info = os.stat(
                candidate,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
            if (
                stat.S_ISREG(candidate_info.st_mode)
                and (candidate_info.st_dev, candidate_info.st_ino)
                == (target.st_dev, target.st_ino)
            ):
                aliases.append(candidate)
        if len(aliases) != target.st_nlink - 1:
            raise RuntimeError(
                f"immutable record has unowned hard-link aliases: {name}"
            )
        for alias in aliases:
            current = os.stat(
                alias,
                dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
            if (current.st_dev, current.st_ino) != (target.st_dev, target.st_ino):
                raise RuntimeError(f"immutable-record alias changed: {alias}")
            os.unlink(alias, dir_fd=parent_descriptor)
        os.fsync(parent_descriptor)
    if read_regular_bytes_at(
        parent_descriptor,
        name,
        require_read_only=True,
    ) != data:
        raise RuntimeError(f"immutable-record recovery failed: {name}")
    return True


def immutable_write_at(
    parent_descriptor: int,
    name: str,
    data: bytes,
    mode: int = 0o400,
) -> bool:
    if recover_local_ready_aliases_at(parent_descriptor, name, data):
        return False
    temporary = f".{name}.{os.getpid()}.{secrets.token_hex(8)}.tmp"
    descriptor = os.open(
        temporary,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC,
        0o600,
        dir_fd=parent_descriptor,
    )
    try:
        view = memoryview(data)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError("short immutable-record write")
            view = view[written:]
        os.fsync(descriptor)
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    installed = False
    try:
        try:
            os.link(
                temporary,
                name,
                src_dir_fd=parent_descriptor,
                dst_dir_fd=parent_descriptor,
                follow_symlinks=False,
            )
            installed = True
        except FileExistsError:
            if not recover_local_ready_aliases_at(parent_descriptor, name, data):
                raise RuntimeError(f"immutable-record creation raced: {name}")
    finally:
        try:
            os.unlink(temporary, dir_fd=parent_descriptor)
        except FileNotFoundError:
            pass
    os.fsync(parent_descriptor)
    if read_regular_bytes_at(
        parent_descriptor,
        name,
        require_read_only=True,
    ) != data:
        raise RuntimeError(f"immutable-record readback differs: {name}")
    return installed


def recover_local_ready_aliases(path: Path, data: bytes) -> bool:
    """Finish a link-then-unlink record commit interrupted after its link step."""

    if not path.exists() and not path.is_symlink():
        return False
    if read_regular_bytes(
        path,
        require_read_only=True,
        require_single_link=False,
    ) != data:
        raise RuntimeError(f"refusing changed immutable record: {path}")
    target = path.lstat()
    prefix = f".{path.name}."
    aliases: list[Path] = []
    if target.st_nlink > 1:
        for candidate in path.parent.iterdir():
            if (
                candidate.name.startswith(prefix)
                and candidate.name.endswith(".tmp")
            ):
                candidate_info = candidate.lstat()
                if (
                    stat.S_ISREG(candidate_info.st_mode)
                    and (candidate_info.st_dev, candidate_info.st_ino)
                    == (target.st_dev, target.st_ino)
                ):
                    aliases.append(candidate)
        if len(aliases) != target.st_nlink - 1:
            raise RuntimeError(
                f"immutable record has unowned hard-link aliases: {path}"
            )
        for alias in aliases:
            current = alias.lstat()
            if (current.st_dev, current.st_ino) != (target.st_dev, target.st_ino):
                raise RuntimeError(f"immutable-record alias changed: {alias}")
            alias.unlink()
        fsync_directory(path.parent)
    if read_regular_bytes(path, require_read_only=True) != data:
        raise RuntimeError(f"immutable-record recovery failed: {path}")
    return True


def sha256_regular(path: Path, *, sync: bool = False) -> tuple[int, str]:
    before_path = path.lstat()
    if not stat.S_ISREG(before_path.st_mode) or before_path.st_nlink != 1:
        raise RuntimeError(f"expected one-link regular file: {path}")
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    digest = hashlib.sha256()
    try:
        before = os.fstat(descriptor)
        if _identity(before) != _identity(before_path):
            raise RuntimeError(f"path changed while opening: {path}")
        while True:
            chunk = os.read(descriptor, CHUNK)
            if not chunk:
                break
            digest.update(chunk)
        if sync:
            os.fsync(descriptor)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    after_path = path.lstat()
    if _identity(before) != _identity(after) or _identity(after) != _identity(after_path):
        raise RuntimeError(f"file changed while hashing: {path}")
    return after.st_size, digest.hexdigest()


def immutable_write(path: Path, data: bytes, mode: int = 0o400) -> bool:
    """Create an fsynced, read-only record without ever replacing an existing inode."""

    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    if recover_local_ready_aliases(path, data):
        return False
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{secrets.token_hex(8)}.tmp")
    descriptor = os.open(
        temporary,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
        0o600,
    )
    try:
        view = memoryview(data)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError("short immutable-record write")
            view = view[written:]
        os.fsync(descriptor)
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    installed = False
    try:
        try:
            os.link(temporary, path, follow_symlinks=False)
            installed = True
        except FileExistsError:
            if not recover_local_ready_aliases(path, data):
                raise RuntimeError(f"immutable-record creation raced: {path}")
    finally:
        temporary.unlink(missing_ok=True)
    fsync_directory(path.parent)
    if read_regular_bytes(path, require_read_only=True) != data:
        raise RuntimeError(f"immutable-record readback differs: {path}")
    return installed


def connection_from_state(
    state_path: Path,
    host_key_sha256: str,
) -> ConnectionConfig:
    info = state_path.lstat()
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1 or info.st_mode & 0o077:
        raise SystemExit(
            f"SSH state must be a private one-link regular file: {state_path}"
        )
    if info.st_uid != os.geteuid():
        raise SystemExit(
            f"SSH state owner uid {info.st_uid} != effective uid {os.geteuid()}"
        )
    payload = json.loads(read_regular_bytes(state_path).decode("utf-8"))
    parsed = urllib.parse.urlparse(payload["data"]["ssh_url"])
    query = urllib.parse.parse_qs(parsed.query, strict_parsing=True)
    if parsed.hostname is None:
        raise SystemExit("SSH URL has no gateway hostname")
    return ConnectionConfig(
        host=parsed.hostname,
        port=int(query["port"][0]),
        username=query["username"][0],
        password=base64.b64decode(query["password"][0], validate=True).decode("utf-8"),
        host_key_sha256=normalize_host_key_sha256(host_key_sha256),
    )


def connect(config: ConnectionConfig) -> tuple[paramiko.SSHClient, paramiko.SFTPClient]:
    last_error: BaseException | None = None
    for attempt in range(1, CONNECT_ATTEMPTS + 1):
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(PinnedHostKeyPolicy(config.host_key_sha256))
        try:
            client.connect(
                hostname=config.host,
                port=config.port,
                username=config.username,
                password=config.password,
                timeout=30,
                banner_timeout=30,
                auth_timeout=30,
                look_for_keys=False,
                allow_agent=False,
            )
            transport = client.get_transport()
            if transport is None:
                raise RuntimeError("SSH transport is unavailable")
            actual = key_sha256(transport.get_remote_server_key())
            if actual != config.host_key_sha256:
                raise paramiko.SSHException(
                    "connected SSH host key changed: "
                    f"{actual} != {config.host_key_sha256}"
                )
            transport.set_keepalive(20)
            sftp = client.open_sftp()
            sftp.get_channel().settimeout(CHANNEL_IO_TIMEOUT_SECONDS)
            return client, sftp
        except (OSError, EOFError, paramiko.SSHException) as error:
            last_error = error
            client.close()
            if attempt == CONNECT_ATTEMPTS:
                break
            delay = CONNECT_BACKOFF_SECONDS * attempt
            with PRINT_LOCK:
                print(
                    f"SSH connect attempt {attempt}/{CONNECT_ATTEMPTS} failed; "
                    f"retrying in {delay:.0f}s: {error}",
                    flush=True,
                )
            time.sleep(delay)
    raise RuntimeError(
        f"SSH connection failed after {CONNECT_ATTEMPTS} attempts: {last_error}"
    )


def safe_relative(value: str) -> PurePosixPath:
    if "\x00" in value:
        raise ValueError("manifest path contains NUL")
    path = PurePosixPath(value)
    if path.is_absolute() or not path.parts or any(
        part in ("", ".", "..") for part in path.parts
    ):
        raise ValueError(f"unsafe manifest path: {value!r}")
    if path.parts[0] in RESERVED_DESTINATION_COMPONENTS:
        raise ValueError(f"manifest path collides with receiver metadata: {value!r}")
    return path


def canonical_remote_absolute(value: str, label: str) -> PurePosixPath:
    if "\x00" in value or not value.startswith("/"):
        raise ValueError(f"{label} must be an absolute POSIX path")
    components = value.split("/")[1:]
    if not components or any(
        component in ("", ".", "..") for component in components
    ):
        raise ValueError(
            f"{label} must be canonical without '.', '..', or empty parts"
        )
    path = PurePosixPath(value)
    if path.as_posix() != value:
        raise ValueError(f"{label} is not canonical: {value!r}")
    return path


def existing_ancestor(path: Path) -> Path:
    candidate = path
    while not candidate.exists():
        if candidate.parent == candidate:
            raise SystemExit(f"no existing destination ancestor: {path}")
        candidate = candidate.parent
    return candidate


def verify_mount(destination: Path, source: str, fstype: str) -> dict[str, str]:
    ancestor = existing_ancestor(destination)
    result = subprocess.run(
        [
            "findmnt", "--json", "--target", str(ancestor),
            "--output", "TARGET,SOURCE,FSTYPE,OPTIONS",
        ],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    filesystems = json.loads(result.stdout).get("filesystems", [])
    if len(filesystems) != 1:
        raise SystemExit(f"could not resolve exactly one mount for {ancestor}")
    mount = filesystems[0]
    if mount.get("source") != source:
        raise SystemExit(
            f"destination mount source {mount.get('source')!r} != expected {source!r}"
        )
    if mount.get("fstype") != fstype:
        raise SystemExit(
            f"destination fstype {mount.get('fstype')!r} != expected {fstype!r}"
        )
    options = set(str(mount.get("options", "")).split(","))
    if "rw" not in options or "ro" in options:
        raise SystemExit(f"destination mount is not read-write: {mount.get('options')}")
    return {key: str(mount[key]) for key in ("target", "source", "fstype", "options")}


def write_probe(destination: DestinationRoot | Path) -> None:
    if isinstance(destination, Path):
        destination.mkdir(parents=True, exist_ok=True, mode=0o700)
        with DestinationRoot(destination) as opened:
            write_probe(opened)
        return
    destination.assert_visible()
    probe = f".athenak-nas-probe-{os.getpid()}-{secrets.token_hex(4)}"
    try:
        descriptor = os.open(
            probe,
            os.O_RDWR
            | os.O_CREAT
            | os.O_EXCL
            | os.O_NOFOLLOW
            | os.O_CLOEXEC,
            0o600,
            dir_fd=destination.descriptor,
        )
        try:
            expected = b"athenak-nas-rw\n"
            os.write(descriptor, expected)
            os.fsync(descriptor)
            os.lseek(descriptor, 0, os.SEEK_SET)
            if os.read(descriptor, len(expected) + 1) != expected:
                raise OSError("NAS probe read-back differs")
        finally:
            os.close(descriptor)
        os.fsync(destination.descriptor)
    finally:
        try:
            os.unlink(probe, dir_fd=destination.descriptor)
        except FileNotFoundError:
            pass
        os.fsync(destination.descriptor)
        destination.assert_visible()


def remote_regular_attributes(
    sftp: paramiko.SFTPClient,
    path: str,
    *,
    require_read_only: bool = False,
    require_single_link: bool = True,
) -> paramiko.SFTPAttributes:
    attributes = sftp.lstat(path)
    mode = int(attributes.st_mode or 0)
    if not stat.S_ISREG(mode):
        raise RuntimeError(f"remote path is not a regular file: {path}")
    links = getattr(attributes, "st_nlink", None)
    if require_single_link and links not in (None, 1):
        raise RuntimeError(f"remote path has unexpected link count {links}: {path}")
    if require_read_only and mode & 0o222:
        raise RuntimeError(f"remote immutable record is writable: {path}")
    return attributes


def read_remote_bytes(
    sftp: paramiko.SFTPClient,
    path: str,
    *,
    require_read_only: bool = False,
    require_single_link: bool = True,
) -> bytes:
    before = remote_regular_attributes(
        sftp,
        path,
        require_read_only=require_read_only,
        require_single_link=require_single_link,
    )
    with sftp.open(path, "rb") as stream:
        data = stream.read()
    after = remote_regular_attributes(
        sftp,
        path,
        require_read_only=require_read_only,
        require_single_link=require_single_link,
    )
    if _remote_identity(before) != _remote_identity(after) or len(data) != after.st_size:
        raise RuntimeError(f"remote file changed while reading: {path}")
    return data


class TrustBoundaryViolation(RuntimeError):
    """A reconnect no longer proves the same remote and local transaction."""


class RetryableRemoteTransport(RuntimeError):
    """A failure known to originate in remote SSH/SFTP transport I/O."""


REMOTE_RETRYABLE_ERRORS = (
    EOFError,
    socket.timeout,
    paramiko.SSHException,
    RetryableRemoteTransport,
)


def retryable_remote_io(operation):
    """Classify only an explicitly remote operation as reconnectable."""

    try:
        return operation()
    except TrustBoundaryViolation:
        raise
    except FileNotFoundError as exc:
        raise TrustBoundaryViolation(f"pinned remote path disappeared: {exc}") from exc
    except (EOFError, OSError, paramiko.SSHException) as exc:
        raise RetryableRemoteTransport(str(exc)) from exc


@dataclass(frozen=True)
class ReconnectValidation:
    """Immutable evidence that every replacement SSH session must re-prove."""

    config: ConnectionConfig
    remote_manifest: str
    manifest_bytes: bytes
    destination: DestinationRoot
    destination_path: Path
    expected_mount_source: str
    expected_mount_fstype: str
    mount_snapshot: dict[str, str]

    def open_validated(
        self,
    ) -> tuple[paramiko.SSHClient, paramiko.SFTPClient]:
        client, sftp = connect(self.config)
        try:
            try:
                current_manifest = retryable_remote_io(
                    lambda: read_remote_bytes(
                        sftp,
                        self.remote_manifest,
                        require_read_only=True,
                    )
                )
            except REMOTE_RETRYABLE_ERRORS:
                raise
            except BaseException as exc:
                raise TrustBoundaryViolation(
                    f"cannot revalidate pinned remote manifest: {exc}"
                ) from exc
            if current_manifest != self.manifest_bytes:
                raise TrustBoundaryViolation(
                    "pinned remote manifest changed after reconnect"
                )
            try:
                self.destination.assert_visible()
                mount = verify_mount(
                    self.destination_path,
                    self.expected_mount_source,
                    self.expected_mount_fstype,
                )
                if mount != self.mount_snapshot:
                    raise TrustBoundaryViolation(
                        "NAS mount identity/options changed after reconnect"
                    )
                write_probe(self.destination)
                self.destination.assert_visible()
            except TrustBoundaryViolation:
                raise
            except BaseException as exc:
                raise TrustBoundaryViolation(
                    f"cannot revalidate destination after reconnect: {exc}"
                ) from exc
            return client, sftp
        except BaseException:
            with contextlib.suppress(BaseException):
                sftp.close()
            with contextlib.suppress(BaseException):
                client.close()
            raise


class ValidatedRemoteSession:
    """One replaceable SSH session with bounded fail-closed recovery."""

    def __init__(
        self,
        validation: ReconnectValidation,
        client: paramiko.SSHClient | None = None,
        sftp: paramiko.SFTPClient | None = None,
    ) -> None:
        if (client is None) != (sftp is None):
            raise ValueError("client and SFTP must be supplied together")
        self.validation = validation
        self.client = client
        self.sftp = sftp
        if self.client is None:
            try:
                self.client, self.sftp = validation.open_validated()
            except REMOTE_RETRYABLE_ERRORS as exc:
                self.recover(exc, 0, "validated SSH session")

    def close(self) -> None:
        if self.sftp is not None:
            with contextlib.suppress(BaseException):
                self.sftp.close()
        if self.client is not None:
            with contextlib.suppress(BaseException):
                self.client.close()
        self.sftp = None
        self.client = None

    def recover(
        self,
        error: BaseException,
        attempts_used: int,
        label: str,
    ) -> int:
        """Replace a failed transport, consuming one bounded exponential retry."""

        self.close()
        last_error: BaseException = error
        while attempts_used < STREAM_RECONNECT_ATTEMPTS:
            attempts_used += 1
            delay = STREAM_RECONNECT_BACKOFF_SECONDS * (2 ** (attempts_used - 1))
            with PRINT_LOCK:
                print(
                    f"{label}: transport failed; reconnect "
                    f"{attempts_used}/{STREAM_RECONNECT_ATTEMPTS} in {delay:.0f}s: "
                    f"{last_error}",
                    flush=True,
                )
            time.sleep(delay)
            try:
                self.client, self.sftp = self.validation.open_validated()
                return attempts_used
            except TrustBoundaryViolation:
                self.close()
                raise
            except REMOTE_RETRYABLE_ERRORS as reconnect_error:
                last_error = reconnect_error
            except RuntimeError as reconnect_error:
                # connect() reports its own bounded attempt exhaustion as RuntimeError.
                last_error = reconnect_error
            self.close()
        raise RuntimeError(
            f"{label}: reconnect limit {STREAM_RECONNECT_ATTEMPTS} exhausted: "
            f"{last_error}"
        ) from last_error

    def remote_call(self, label: str, operation):
        attempts = 0
        while True:
            assert self.client is not None and self.sftp is not None
            try:
                return operation(self.client, self.sftp)
            except REMOTE_RETRYABLE_ERRORS as exc:
                attempts = self.recover(exc, attempts, label)

    def __enter__(self) -> "ValidatedRemoteSession":
        return self

    def __exit__(self, *unused: object) -> None:
        self.close()


def remote_mkdirs(sftp: paramiko.SFTPClient, path: str) -> None:
    current = "/"
    for part in PurePosixPath(path).parts[1:]:
        current = posixpath.join(current, part)
        try:
            attributes = sftp.lstat(current)
        except FileNotFoundError:
            sftp.mkdir(current, mode=0o700)
            attributes = sftp.lstat(current)
        if not stat.S_ISDIR(int(attributes.st_mode or 0)):
            raise RuntimeError(f"remote ACK ancestor is not a directory: {current}")


def exec_checked(client: paramiko.SSHClient, argv: list[str]) -> None:
    command = shlex.join(argv)
    _, stdout, _ = client.exec_command(command, timeout=CHANNEL_IO_TIMEOUT_SECONDS)
    channel = stdout.channel
    channel.settimeout(CHANNEL_IO_TIMEOUT_SECONDS)
    deadline = time.monotonic() + CHANNEL_IO_TIMEOUT_SECONDS
    error_chunks: list[bytes] = []
    while not retryable_remote_io(channel.exit_status_ready):
        while retryable_remote_io(channel.recv_stderr_ready):
            error_chunks.append(retryable_remote_io(lambda: channel.recv_stderr(CHUNK)))
        if time.monotonic() >= deadline:
            raise RetryableRemoteTransport(
                f"timed out waiting for remote command: {command}"
            )
        time.sleep(CHANNEL_STATUS_POLL_SECONDS)
    while retryable_remote_io(channel.recv_stderr_ready):
        error_chunks.append(retryable_remote_io(lambda: channel.recv_stderr(CHUNK)))
    status = retryable_remote_io(channel.recv_exit_status)
    error = b"".join(error_chunks).decode("utf-8", errors="replace").strip()
    if status != 0:
        raise RuntimeError(f"remote command failed ({status}): {command}: {error}")


def remote_fsync(
    client: paramiko.SSHClient,
    file_path: str,
    directory_path: str,
) -> None:
    program = (
        "import os,sys; "
        "f=os.open(sys.argv[1],os.O_RDONLY|getattr(os,'O_NOFOLLOW',0)); "
        "os.fsync(f); os.close(f); "
        "d=os.open(sys.argv[2],os.O_RDONLY|os.O_DIRECTORY|"
        "getattr(os,'O_NOFOLLOW',0)); os.fsync(d); os.close(d)"
    )
    exec_checked(client, ["/usr/bin/python3", "-c", program, file_path, directory_path])


def remote_immutable_write(
    client: paramiko.SSHClient,
    sftp: paramiko.SFTPClient,
    path: str,
    data: bytes,
) -> bool:
    """Atomically link a complete ACK into place without replacing an old ACK."""

    parent = posixpath.dirname(path)
    remote_mkdirs(sftp, parent)
    try:
        existing = read_remote_bytes(
            sftp,
            path,
            require_read_only=True,
            require_single_link=False,
        )
    except FileNotFoundError:
        existing = None
    if existing is not None:
        if existing != data:
            raise RuntimeError(f"refusing to replace remote immutable record: {path}")
        prefix = f".{PurePosixPath(path).name}."
        for attributes in sftp.listdir_attr(parent):
            name = str(attributes.filename)
            if not name.startswith(prefix) or not name.endswith(".tmp"):
                continue
            candidate = posixpath.join(parent, name)
            try:
                candidate_data = read_remote_bytes(
                    sftp,
                    candidate,
                    require_read_only=True,
                    require_single_link=False,
                )
            except RetryableRemoteTransport:
                raise
            except (FileNotFoundError, RuntimeError):
                # An incomplete pre-link temporary record is not an alias of the
                # ready ACK.  Leave it untouched for a separate stale-temp audit.
                continue
            if candidate_data == data:
                sftp.remove(candidate)
        remote_fsync(client, path, parent)
        if read_remote_bytes(sftp, path, require_read_only=True) != data:
            raise RuntimeError(f"remote immutable-record recovery failed: {path}")
        return False

    temporary = posixpath.join(
        parent,
        f".{PurePosixPath(path).name}.{os.getpid()}.{secrets.token_hex(8)}.tmp",
    )
    created = False
    try:
        with sftp.open(temporary, "wx") as stream:
            stream.write(data)
            stream.flush()
        sftp.chmod(temporary, 0o400)
        remote_fsync(client, temporary, parent)
        try:
            exec_checked(client, ["/usr/bin/ln", "--", temporary, path])
            created = True
        except RuntimeError as link_error:
            try:
                raced = read_remote_bytes(sftp, path, require_read_only=True)
            except FileNotFoundError:
                # Preserve whether `ln` itself failed deterministically or its
                # response was lost before the link existed.  Re-raising the
                # transport error lets the outer validated transaction reconnect;
                # a deterministic command failure remains fail-closed.
                raise link_error
            if raced != data:
                raise RuntimeError(f"remote immutable-record creation raced: {path}")
    finally:
        try:
            sftp.remove(temporary)
        except FileNotFoundError:
            pass
    remote_fsync(client, path, parent)
    if read_remote_bytes(sftp, path, require_read_only=True) != data:
        raise RuntimeError(f"remote immutable-record readback differs: {path}")
    return created


def validate_manifest(
    manifest: object,
    manifest_name: str,
    remote_root: str,
    expected_segment: str,
) -> list[dict[str, object]]:
    if not isinstance(manifest, dict) or manifest.get("schema") != 1:
        raise ValueError(f"unsupported schema in {manifest_name}")
    if manifest.get("root") != remote_root:
        raise ValueError(f"remote root mismatch in {manifest_name}")
    if manifest.get("segment") != expected_segment:
        raise ValueError(f"segment mismatch in {manifest_name}")
    files = manifest.get("files")
    if not isinstance(files, list) or not files:
        raise ValueError(f"no files in {manifest_name}")
    seen: set[str] = set()
    validated: list[dict[str, object]] = []
    for raw in files:
        if not isinstance(raw, dict) or not isinstance(raw.get("path"), str):
            raise ValueError(f"non-object or invalid path record in {manifest_name}")
        relative = safe_relative(raw["path"])
        size = raw.get("size")
        digest = raw.get("sha256")
        if isinstance(size, bool) or not isinstance(size, int) or size < 0:
            raise ValueError(f"invalid size for {relative}")
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or any(char not in "0123456789abcdef" for char in digest)
        ):
            raise ValueError(f"invalid SHA256 for {relative}")
        if relative.as_posix() in seen:
            raise ValueError(f"duplicate path in manifest: {relative}")
        seen.add(relative.as_posix())
        validated.append({"path": relative, "size": size, "sha256": digest})
    return validated


def open_append_nofollow_at(
    parent_descriptor: int,
    name: str,
) -> tuple[object, int]:
    safe_directory_component(name)
    descriptor = os.open(
        name,
        os.O_WRONLY | os.O_APPEND | os.O_CREAT | os.O_NOFOLLOW | os.O_CLOEXEC,
        0o600,
        dir_fd=parent_descriptor,
    )
    info = os.fstat(descriptor)
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        os.close(descriptor)
        raise RuntimeError(f"partial is not a one-link regular file: {name}")
    # Unbuffered writes make the stat-derived resume offset exactly the bytes that
    # recv() returned, even when an SSH exception tears down the stream mid-window.
    return os.fdopen(descriptor, "ab", buffering=0), info.st_size


def _append_validated_window(
    session: ValidatedRemoteSession,
    remote_path: str,
    expected_remote_identity: tuple[int, int, int, int],
    parent_descriptor: int,
    piece_name: str,
    range_start: int,
    target_size: int,
    label: str,
) -> None:
    """Append through target_size, safely resuming after a torn SSH channel."""

    def operation(client: paramiko.SSHClient, sftp: paramiko.SFTPClient) -> None:
        remote_info = retryable_remote_io(
            lambda: remote_regular_attributes(sftp, remote_path)
        )
        if _remote_identity(remote_info) != expected_remote_identity:
            raise TrustBoundaryViolation(
                f"remote source identity changed after reconnect: {remote_path}"
            )
        try:
            current = stat_regular_at(parent_descriptor, piece_name).st_size
        except FileNotFoundError:
            current = 0
        if current > target_size:
            raise RuntimeError(f"partial grew past transfer window: {piece_name}")
        if current == target_size:
            return
        remaining = target_size - current
        command = (
            f"dd if={shlex.quote(remote_path)} iflag=skip_bytes,count_bytes "
            f"skip={range_start + current} count={remaining} status=none"
        )
        _, stdout, _ = retryable_remote_io(
            lambda: client.exec_command(command, timeout=None)
        )
        channel = stdout.channel
        channel.settimeout(CHANNEL_IO_TIMEOUT_SECONDS)
        stream, opened_size = open_append_nofollow_at(
            parent_descriptor,
            piece_name,
        )
        if opened_size != current:
            stream.close()
            raise RuntimeError(f"partial changed while opening: {piece_name}")
        try:
            while current < target_size:
                chunk = retryable_remote_io(
                    lambda: channel.recv(min(CHUNK, target_size - current))
                )
                if not chunk:
                    raise EOFError(
                        f"early EOF for {piece_name} at {current}/{target_size}"
                    )
                view = memoryview(chunk)
                while view:
                    written = stream.write(view)
                    if written is None or written <= 0:
                        raise OSError("short local partial write")
                    current += written
                    view = view[written:]
            stream.flush()
            os.fsync(stream.fileno())
            os.fsync(parent_descriptor)
        except BaseException as transfer_error:
            # Preserve all complete bytes already returned by recv().  The next
            # validated connection derives its skip offset only from this durable
            # inode.  A local flush/fsync error must override an SSH error so it can
            # never be misclassified as reconnectable.
            try:
                stream.flush()
                os.fsync(stream.fileno())
                os.fsync(parent_descriptor)
            except BaseException as durability_error:
                raise durability_error from transfer_error
            raise
        finally:
            stream.close()
        deadline = time.monotonic() + CHANNEL_IO_TIMEOUT_SECONDS
        while not retryable_remote_io(channel.exit_status_ready):
            if time.monotonic() >= deadline:
                raise RetryableRemoteTransport(
                    f"timed out waiting for remote dd exit status: {piece_name}"
                )
            time.sleep(CHANNEL_STATUS_POLL_SECONDS)
        status = retryable_remote_io(channel.recv_exit_status)
        if status != 0:
            raise paramiko.SSHException(
                f"remote dd failed for {piece_name}: {status}"
            )
        remote_after = retryable_remote_io(
            lambda: remote_regular_attributes(sftp, remote_path)
        )
        if _remote_identity(remote_after) != expected_remote_identity:
            raise TrustBoundaryViolation(
                f"remote source changed while transferring: {remote_path}"
            )

    session.remote_call(label, operation)


def download_range(
    config: ConnectionConfig,
    remote_path: str,
    parent_descriptor: int,
    piece_name: str,
    start: int,
    end: int,
    number: int,
    validation: ReconnectValidation | None = None,
    expected_remote_identity: tuple[int, int, int, int] | None = None,
) -> tuple[int, float]:
    try:
        expected = end - start
        try:
            offset = stat_regular_at(parent_descriptor, piece_name).st_size
        except FileNotFoundError:
            offset = 0
        if offset > expected:
            raise RuntimeError(f"range {number} is too large: {piece_name}")
        if offset == expected:
            return expected, 0.0
        if validation is None or expected_remote_identity is None:
            raise RuntimeError("incomplete reconnect validation for remote range")
        started = time.monotonic()
        while offset < expected:
            target = min(expected, offset + STREAM_WINDOW_BYTES)
            # Proactively use a fresh pinned and revalidated transport for each
            # <=256 MiB window, safely below Paramiko 2.12's 512 MiB rekey point.
            with ValidatedRemoteSession(validation) as session:
                _append_validated_window(
                    session,
                    remote_path,
                    expected_remote_identity,
                    parent_descriptor,
                    piece_name,
                    start,
                    target,
                    f"range {number} {piece_name}",
                )
            offset = stat_regular_at(parent_descriptor, piece_name).st_size
            if offset != target:
                raise RuntimeError(f"range {number} window size mismatch")
            os.fsync(parent_descriptor)
            elapsed = max(time.monotonic() - started, 1e-9)
            with PRINT_LOCK:
                print(
                    f"{piece_name}: {100.0 * offset / expected:.1f}% "
                    f"({offset}/{expected}), "
                    f"{(offset / elapsed) / (1024 * 1024):.2f} MiB/s",
                    flush=True,
                )
        info, _ = sha256_regular_at(parent_descriptor, piece_name)
        if info.st_size != expected:
            raise RuntimeError(f"range {number} size mismatch")
        return expected, time.monotonic() - started
    finally:
        os.close(parent_descriptor)


def parallel_download(
    config: ConnectionConfig,
    remote_path: str,
    parent_descriptor: int,
    partial_name: str,
    expected_size: int,
    validation: ReconnectValidation | None = None,
    expected_remote_identity: tuple[int, int, int, int] | None = None,
) -> None:
    piece_span = (expected_size + PARALLEL_STREAMS - 1) // PARALLEL_STREAMS
    piece_span = ((piece_span + CHUNK - 1) // CHUNK) * CHUNK
    ranges = []
    for number in range(PARALLEL_STREAMS):
        start = number * piece_span
        end = min(expected_size, start + piece_span)
        if start < end:
            ranges.append((number, start, end, f"{partial_name}.piece{number:02d}"))
    try:
        partial_info = stat_regular_at(parent_descriptor, partial_name)
    except FileNotFoundError:
        partial_info = None
    if partial_info is not None and partial_info.st_size:
        _, first_start, first_end, first_piece = ranges[0]
        try:
            stat_regular_at(parent_descriptor, first_piece)
            first_piece_exists = True
        except FileNotFoundError:
            first_piece_exists = False
        if first_piece_exists or partial_info.st_size > first_end - first_start:
            raise RuntimeError("ambiguous or oversized sequential parallel handoff")
        os.replace(
            partial_name,
            first_piece,
            src_dir_fd=parent_descriptor,
            dst_dir_fd=parent_descriptor,
        )
        os.fsync(parent_descriptor)
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(ranges)) as executor:
        futures = [
            executor.submit(
                download_range,
                config,
                remote_path,
                os.dup(parent_descriptor),
                piece_name,
                start,
                end,
                number,
                validation,
                expected_remote_identity,
            )
            for number, start, end, piece_name in ranges
        ]
        for future in concurrent.futures.as_completed(futures):
            future.result()

    assembling_name = f"{partial_name}.assembling"
    try:
        stale = stat_regular_at(parent_descriptor, assembling_name)
    except FileNotFoundError:
        stale = None
    if stale is not None:
        # `.assembling` is never authoritative: only the verified `.part` name
        # can be installed.  Remove this pinned, one-link temporary inode so a
        # crash during assembly is deterministically resumable.
        os.unlink(assembling_name, dir_fd=parent_descriptor)
        os.fsync(parent_descriptor)
    descriptor = os.open(
        assembling_name,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC,
        0o600,
        dir_fd=parent_descriptor,
    )
    try:
        with os.fdopen(descriptor, "wb", closefd=False) as output:
            for number, start, end, piece_name in ranges:
                info, _ = sha256_regular_at(parent_descriptor, piece_name)
                if info.st_size != end - start:
                    raise RuntimeError(f"range {number} changed before assembly")
                source_descriptor = os.open(
                    piece_name,
                    os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC,
                    dir_fd=parent_descriptor,
                )
                with os.fdopen(source_descriptor, "rb") as source:
                    shutil.copyfileobj(source, output, length=CHUNK)
            output.flush()
            os.fsync(output.fileno())
    finally:
        os.close(descriptor)
    if stat_regular_at(parent_descriptor, assembling_name).st_size != expected_size:
        raise RuntimeError("assembled partial size mismatch")
    os.replace(
        assembling_name,
        partial_name,
        src_dir_fd=parent_descriptor,
        dst_dir_fd=parent_descriptor,
    )
    for _, _, _, piece_name in ranges:
        os.unlink(piece_name, dir_fd=parent_descriptor)
    os.fsync(parent_descriptor)


def install_no_replace_at(
    partial_parent: int,
    partial_name: str,
    final_parent: int,
    final_name: str,
) -> None:
    before = stat_regular_at(partial_parent, partial_name)
    try:
        os.link(
            partial_name,
            final_name,
            src_dir_fd=partial_parent,
            dst_dir_fd=final_parent,
            follow_symlinks=False,
        )
    except FileExistsError as exc:
        raise RuntimeError(f"refusing to replace final file: {final_name}") from exc
    linked = stat_regular_at(
        final_parent,
        final_name,
        require_single_link=False,
    )
    if (
        (linked.st_dev, linked.st_ino) != (before.st_dev, before.st_ino)
        or linked.st_nlink != 2
    ):
        raise RuntimeError(f"final link identity mismatch: {final_name}")
    os.unlink(partial_name, dir_fd=partial_parent)
    installed = stat_regular_at(final_parent, final_name)
    if (
        (installed.st_dev, installed.st_ino) != (before.st_dev, before.st_ino)
        or installed.st_nlink != 1
    ):
        raise RuntimeError(f"final identity changed during commit: {final_name}")
    os.fsync(partial_parent)
    if final_parent != partial_parent:
        os.fsync(final_parent)


def install_no_replace(partial: Path, final_path: Path) -> None:
    """Path wrapper retained for unit tests; production uses pinned dirfds."""

    partial_parent = os.open(
        partial.parent,
        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
    )
    final_path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    final_parent = os.open(
        final_path.parent,
        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
    )
    try:
        install_no_replace_at(
            partial_parent,
            partial.name,
            final_parent,
            final_path.name,
        )
    finally:
        os.close(final_parent)
        os.close(partial_parent)


def transfer_file(
    config: ConnectionConfig,
    validation: ReconnectValidation,
    remote_root: str,
    destination: DestinationRoot,
    incoming_segment: str,
    record: dict[str, object],
    index: int,
    count: int,
) -> None:
    relative = record["path"]
    assert isinstance(relative, PurePosixPath)
    expected_size = int(record["size"])
    expected_digest = str(record["sha256"])
    remote_path = posixpath.join(remote_root, relative.as_posix())
    final_parent, final_name = destination.open_relative_parent(
        relative,
        create=True,
    )
    partial_parent, _ = destination.open_relative_parent(
        relative,
        prefix=(".incoming", incoming_segment),
        create=True,
    )
    partial_name = f"{relative.name}.part"
    try:
        try:
            final_info, digest = sha256_regular_at(
                final_parent,
                final_name,
                require_single_link=False,
            )
        except FileNotFoundError:
            final_info = None
        if final_info is not None:
            if final_info.st_size != expected_size or digest != expected_digest:
                raise RuntimeError(
                    f"refusing mismatched final file: {relative}"
                )
            if final_info.st_nlink == 2:
                try:
                    partial_info = stat_regular_at(
                        partial_parent,
                        partial_name,
                        require_single_link=False,
                    )
                except FileNotFoundError as exc:
                    raise RuntimeError(
                        f"final has an unowned hard-link alias: {relative}"
                    ) from exc
                if (partial_info.st_dev, partial_info.st_ino) != (
                    final_info.st_dev,
                    final_info.st_ino,
                ):
                    raise RuntimeError(
                        f"final hard-link alias is not its partial: {relative}"
                    )
                os.unlink(partial_name, dir_fd=partial_parent)
                os.fsync(partial_parent)
                os.fsync(final_parent)
                final_info = stat_regular_at(final_parent, final_name)
            if final_info.st_nlink != 1:
                raise RuntimeError(
                    f"final has unexpected link count {final_info.st_nlink}: {relative}"
                )
            print(f"[{index}/{count}] verified existing {relative}", flush=True)
            return

        with ValidatedRemoteSession(validation) as source_session:
            remote_before = source_session.remote_call(
                f"source stat {relative}",
                lambda _client, current_sftp: retryable_remote_io(
                    lambda: remote_regular_attributes(
                        current_sftp,
                        remote_path,
                    )
                ),
            )
        if remote_before.st_size != expected_size:
            raise RuntimeError(
                f"remote size mismatch for {relative}: "
                f"{remote_before.st_size} != {expected_size}"
            )
        try:
            offset = stat_regular_at(partial_parent, partial_name).st_size
        except FileNotFoundError:
            offset = 0
        if offset > expected_size:
            raise RuntimeError(f"partial file is larger than source: {relative}")
        if expected_size == 0 and offset == 0:
            try:
                stat_regular_at(partial_parent, partial_name)
            except FileNotFoundError:
                stream, opened_size = open_append_nofollow_at(
                    partial_parent,
                    partial_name,
                )
                stream.close()
                if opened_size != 0:
                    raise RuntimeError(f"new empty partial is not empty: {relative}")
        if expected_size >= PARALLEL_THRESHOLD and offset < expected_size:
            parallel_download(
                config,
                remote_path,
                partial_parent,
                partial_name,
                expected_size,
                validation,
                _remote_identity(remote_before),
            )
            offset = expected_size
        if offset < expected_size:
            download_range(
                config,
                remote_path,
                os.dup(partial_parent),
                partial_name,
                0,
                expected_size,
                0,
                validation,
                _remote_identity(remote_before),
            )
            offset = expected_size
        partial_info, digest = sha256_regular_at(
            partial_parent,
            partial_name,
            sync=True,
        )
        if partial_info.st_size != expected_size or digest != expected_digest:
            os.unlink(partial_name, dir_fd=partial_parent)
            os.fsync(partial_parent)
            raise RuntimeError(
                f"local verification failed and corrupt partial was discarded for "
                f"{relative}: size={partial_info.st_size}, sha256={digest}"
            )
        with ValidatedRemoteSession(validation) as source_session:
            remote_after = source_session.remote_call(
                f"final source stat {relative}",
                lambda _client, current_sftp: retryable_remote_io(
                    lambda: remote_regular_attributes(
                        current_sftp,
                        remote_path,
                    )
                ),
            )
        if _remote_identity(remote_before) != _remote_identity(remote_after):
            raise RuntimeError(f"remote source changed while transferring: {relative}")
        install_no_replace_at(
            partial_parent,
            partial_name,
            final_parent,
            final_name,
        )
        print(f"[{index}/{count}] committed {relative}", flush=True)
    finally:
        os.close(partial_parent)
        os.close(final_parent)


def verify_final_files(
    destination: DestinationRoot | Path,
    records: list[dict[str, object]],
    *,
    seal_read_only: bool = True,
) -> dict[str, tuple[int, int, int, int, int, int]]:
    if isinstance(destination, Path):
        with DestinationRoot(destination) as opened:
            return verify_final_files(
                opened,
                records,
                seal_read_only=seal_read_only,
            )
    proofs: dict[str, tuple[int, int, int, int, int, int]] = {}
    for record in records:
        relative = record["path"]
        assert isinstance(relative, PurePosixPath)
        parent, name = destination.open_relative_parent(relative)
        try:
            info, digest = sha256_regular_at(
                parent,
                name,
                sync=True,
                seal_read_only=seal_read_only,
            )
            if info.st_size != record["size"] or digest != record["sha256"]:
                raise RuntimeError(f"final NAS verification failed: {relative}")
            proofs[relative.as_posix()] = _identity(info)
            os.fsync(parent)
        finally:
            os.close(parent)
    os.fsync(destination.descriptor)
    destination.assert_visible()
    return proofs


def verify_final_identities(
    destination: DestinationRoot,
    proofs: dict[str, tuple[int, int, int, int, int, int]],
) -> None:
    """Cheap post-ACK check for any metadata-visible concurrent modification."""

    for relative_text, expected in proofs.items():
        relative = PurePosixPath(relative_text)
        parent, name = destination.open_relative_parent(relative)
        try:
            actual = stat_regular_at(
                parent,
                name,
                require_read_only=True,
            )
            if _identity(actual) != expected:
                raise RuntimeError(
                    f"final NAS identity changed during ACK publication: {relative}"
                )
        finally:
            os.close(parent)
    destination.assert_visible()


def ack_core(
    manifest_name: str,
    manifest_bytes: bytes,
    segment: str,
    destination: Path,
    mount: dict[str, str],
    records: list[dict[str, object]],
) -> dict[str, object]:
    return {
        "schema": 1,
        "kind": "athenak_transfer_ack_v2",
        # A transfer ACK is deliberately never a cloud-deletion capability.
        # The current independent cleanup gate fails closed; only a future
        # storage-barrier-backed verifier and consumer may qualify deletion.
        "authorizes_remote_deletion": False,
        "manifest": manifest_name,
        "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
        "segment": segment,
        "destination": str(destination),
        "mount": mount,
        "file_count": len(records),
        "total_bytes": sum(int(record["size"]) for record in records),
        "files": [
            {
                "path": record["path"].as_posix(),
                "size": record["size"],
                "sha256": record["sha256"],
            }
            for record in records
        ],
    }


class DuplicateJsonKey(ValueError):
    """A security record contains a key with ambiguous duplicate values."""


def reject_duplicate_json_keys(
    pairs: list[tuple[str, object]],
) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise DuplicateJsonKey(f"duplicate JSON key: {key!r}")
        value[key] = item
    return value


def reject_json_constant(value: str) -> object:
    raise ValueError(f"non-finite JSON constant: {value}")


def validate_ack_bytes(data: bytes, expected_core: dict[str, object]) -> None:
    try:
        value = json.loads(
            data,
            object_pairs_hook=reject_duplicate_json_keys,
            parse_constant=reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise RuntimeError("invalid existing ACK JSON") from exc
    if not isinstance(value, dict):
        raise RuntimeError("existing ACK is not an object")
    canonical = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")
    if data != canonical:
        raise RuntimeError("existing ACK bytes are not canonical JSON")
    if set(value) != set(expected_core) | {"verified_unix"}:
        raise RuntimeError("existing ACK has unexpected fields")
    timestamp = value.pop("verified_unix")
    if (
        isinstance(timestamp, bool)
        or not isinstance(timestamp, (int, float))
        or not math.isfinite(float(timestamp))
        or float(timestamp) <= 0
    ):
        raise RuntimeError("existing ACK has invalid verified_unix")
    if value != expected_core:
        raise RuntimeError("existing ACK does not match the verified transfer")


def choose_ack_bytes(
    expected_core: dict[str, object],
    local: bytes | None,
    remote: bytes | None,
    *,
    now: float | None = None,
) -> bytes:
    for candidate in (local, remote):
        if candidate is not None:
            validate_ack_bytes(candidate, expected_core)
    if local is not None and remote is not None and local != remote:
        raise RuntimeError("local and remote immutable ACKs differ")
    if local is not None:
        return local
    if remote is not None:
        return remote
    value = dict(expected_core)
    value["verified_unix"] = time.time() if now is None else now
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def optional_local_record(path: Path) -> bytes | None:
    if not path.exists() and not path.is_symlink():
        return None
    return read_regular_bytes(path, require_read_only=True)


def optional_local_record_at(parent_descriptor: int, name: str) -> bytes | None:
    try:
        return read_regular_bytes_at(
            parent_descriptor,
            name,
            require_read_only=True,
        )
    except FileNotFoundError:
        return None


def optional_remote_record(sftp: paramiko.SFTPClient, path: str) -> bytes | None:
    try:
        return read_remote_bytes(sftp, path, require_read_only=True)
    except FileNotFoundError:
        return None


class SegmentLock:
    def __init__(self, destination: DestinationRoot, segment: str) -> None:
        self.destination = destination
        self.name = f"{segment}.lock"
        self.descriptor: int | None = None
        self.parent_descriptor: int | None = None

    def __enter__(self) -> "SegmentLock":
        self.parent_descriptor = self.destination.open_directory(
            (".locks",),
            create=True,
        )
        try:
            self.descriptor = os.open(
                self.name,
                os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW | os.O_CLOEXEC,
                0o600,
                dir_fd=self.parent_descriptor,
            )
        except BaseException:
            os.close(self.parent_descriptor)
            self.parent_descriptor = None
            raise
        info = os.fstat(self.descriptor)
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            self.__exit__()
            raise RuntimeError(
                f"segment lock is not a one-link regular file: {self.name}"
            )
        try:
            fcntl.flock(self.descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            self.__exit__()
            raise RuntimeError(
                f"another pull owns segment lock: {self.name}"
            ) from exc
        except BaseException:
            self.__exit__()
            raise
        os.ftruncate(self.descriptor, 0)
        os.write(self.descriptor, f"pid={os.getpid()}\n".encode("ascii"))
        os.fsync(self.descriptor)
        os.fsync(self.parent_descriptor)
        return self

    def __exit__(self, *unused: object) -> None:
        if self.descriptor is not None:
            fcntl.flock(self.descriptor, fcntl.LOCK_UN)
            os.close(self.descriptor)
            self.descriptor = None
        if self.parent_descriptor is not None:
            os.close(self.parent_descriptor)
            self.parent_descriptor = None


def main() -> int:
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
    args = parser.parse_args()
    if not args.segment or any(
        char not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_."
        for char in args.segment
    ):
        parser.error("segment may contain only letters, digits, '-', '_' and '.'")
    if (
        len(args.expected_manifest_sha256) != 64
        or any(char not in "0123456789abcdef" for char in args.expected_manifest_sha256)
    ):
        parser.error("--expected-manifest-sha256 must be 64 lowercase hex digits")
    destination = args.destination.expanduser().resolve()
    try:
        remote_root_path = canonical_remote_absolute(args.remote_root, "remote root")
        remote_manifest_path = canonical_remote_absolute(
            args.remote_manifest, "remote manifest"
        )
    except ValueError as exc:
        parser.error(str(exc))
    args.remote_root = remote_root_path.as_posix()
    manifest_name = remote_manifest_path.name
    if not manifest_name.endswith(".manifest.ready"):
        raise SystemExit("remote manifest must end in .manifest.ready")
    config = connection_from_state(args.ssh_state.expanduser(), args.host_key_sha256)
    mount_ancestor = verify_mount(
        destination, args.expected_mount_source, args.expected_mount_fstype
    )
    destination.mkdir(parents=True, exist_ok=True, mode=0o700)
    with DestinationRoot(destination) as local:
        mount_before = verify_mount(
            destination,
            args.expected_mount_source,
            args.expected_mount_fstype,
        )
        if mount_before != mount_ancestor:
            raise RuntimeError(
                "NAS mount identity/options changed while creating destination"
            )
        write_probe(local)
        with SegmentLock(local, args.segment):
            client, sftp = connect(config)
            try:
                manifest_bytes = read_remote_bytes(
                    sftp, args.remote_manifest, require_read_only=True
                )
                actual_manifest_sha = hashlib.sha256(manifest_bytes).hexdigest()
                if actual_manifest_sha != args.expected_manifest_sha256:
                    raise RuntimeError(
                        f"remote manifest SHA256 {actual_manifest_sha} != pinned "
                        f"{args.expected_manifest_sha256}"
                    )
                manifest = json.loads(manifest_bytes)
                records = validate_manifest(
                    manifest,
                    manifest_name,
                    args.remote_root,
                    args.segment,
                )
                validation = ReconnectValidation(
                    config=config,
                    remote_manifest=args.remote_manifest,
                    manifest_bytes=manifest_bytes,
                    destination=local,
                    destination_path=destination,
                    expected_mount_source=args.expected_mount_source,
                    expected_mount_fstype=args.expected_mount_fstype,
                    mount_snapshot=mount_before,
                )
                manifest_parent = local.open_directory(
                    (".manifests",),
                    create=True,
                )
                try:
                    immutable_write_at(
                        manifest_parent,
                        manifest_name,
                        manifest_bytes,
                    )
                finally:
                    os.close(manifest_parent)

                # The bootstrap session has served its only purpose.  Transfers
                # below use bounded <=256 MiB windows on freshly revalidated
                # sessions, so no data channel reaches Paramiko's rekey boundary.
                sftp.close()
                client.close()
                sftp = None
                client = None

                unverified_sizes: list[int] = []
                for record in records:
                    relative = record["path"]
                    assert isinstance(relative, PurePosixPath)
                    parent, name = local.open_relative_parent(
                        relative,
                        create=True,
                    )
                    try:
                        try:
                            info, digest = sha256_regular_at(
                                parent,
                                name,
                                require_single_link=False,
                            )
                        except FileNotFoundError:
                            unverified_sizes.append(int(record["size"]))
                        else:
                            if (
                                info.st_size != record["size"]
                                or digest != record["sha256"]
                            ):
                                raise RuntimeError(
                                    f"mismatched pre-existing final file: {relative}"
                                )
                    finally:
                        os.close(parent)
                extra_assembly = max(
                    (size for size in unverified_sizes if size >= PARALLEL_THRESHOLD),
                    default=0,
                )
                required = sum(unverified_sizes) + extra_assembly + MIN_RESERVE
                free = shutil.disk_usage(existing_ancestor(destination)).free
                if free < required:
                    raise SystemExit(
                        f"insufficient NAS space: free={free}, "
                        f"conservative_required={required}"
                    )

                for index, record in enumerate(records, start=1):
                    transfer_file(
                        config,
                        validation,
                        args.remote_root,
                        local,
                        args.segment,
                        record,
                        index,
                        len(records),
                    )

                with ValidatedRemoteSession(validation) as verify_session:
                    verified_manifest = verify_session.remote_call(
                        "post-transfer remote manifest readback",
                        lambda _client, current_sftp: retryable_remote_io(
                            lambda: read_remote_bytes(
                                current_sftp,
                                args.remote_manifest,
                                require_read_only=True,
                            )
                        ),
                    )
                    if verified_manifest != manifest_bytes:
                        raise RuntimeError("remote manifest changed during transfer")
                mount_mid = verify_mount(
                    destination,
                    args.expected_mount_source,
                    args.expected_mount_fstype,
                )
                if mount_mid != mount_before:
                    raise RuntimeError(
                        "NAS mount identity/options changed during transfer"
                    )
                write_probe(local)
                final_proofs = verify_final_files(local, records)
                mount_final = verify_mount(
                    destination,
                    args.expected_mount_source,
                    args.expected_mount_fstype,
                )
                if mount_final != mount_before:
                    raise RuntimeError("NAS mount identity/options changed before ACK")
                expected_core = ack_core(
                    manifest_name,
                    manifest_bytes,
                    args.segment,
                    destination,
                    mount_final,
                    records,
                )
                local_ack_name = f"{manifest_name}.ack"
                ack_parent = local.open_directory((".acks",), create=True)
                with ValidatedRemoteSession(validation) as ack_session:
                    try:
                        remote_ack = posixpath.join(
                            posixpath.dirname(args.remote_manifest),
                            "acks",
                            f"{manifest_name}.ack",
                        )
                        ack_bytes = choose_ack_bytes(
                            expected_core,
                            optional_local_record_at(ack_parent, local_ack_name),
                            ack_session.remote_call(
                                "read existing remote ACK",
                                lambda _client, current_sftp: retryable_remote_io(
                                    lambda: optional_remote_record(
                                        current_sftp,
                                        remote_ack,
                                    )
                                ),
                            ),
                        )
                        immutable_write_at(ack_parent, local_ack_name, ack_bytes)
                        mount_ack = verify_mount(
                            destination,
                            args.expected_mount_source,
                            args.expected_mount_fstype,
                        )
                        if mount_ack != mount_before:
                            raise RuntimeError(
                                "NAS mount identity/options changed while writing ACK"
                            )
                        if read_regular_bytes_at(
                            ack_parent,
                            local_ack_name,
                            require_read_only=True,
                        ) != ack_bytes:
                            raise RuntimeError(
                                "local ACK readback differs before remote commit"
                            )
                        local.assert_directory_visible((".acks",), ack_parent)
                        visible_ack_parent = local.open_directory((".acks",))
                        try:
                            if read_regular_bytes_at(
                                visible_ack_parent,
                                local_ack_name,
                                require_read_only=True,
                            ) != ack_bytes:
                                raise RuntimeError(
                                    "visible local ACK differs before remote commit"
                                )
                        finally:
                            os.close(visible_ack_parent)
                        # Rebind the local payload immediately before publishing the
                        # non-authorizing receipt. No receipt authorizes cleanup.
                        verify_final_identities(local, final_proofs)
                        ack_session.remote_call(
                            "publish immutable remote ACK",
                            lambda current_client, current_sftp: retryable_remote_io(
                                lambda: remote_immutable_write(
                                    current_client,
                                    current_sftp,
                                    remote_ack,
                                    ack_bytes,
                                )
                            ),
                        )
                        mount_committed = verify_mount(
                            destination,
                            args.expected_mount_source,
                            args.expected_mount_fstype,
                        )
                        if mount_committed != mount_before:
                            raise RuntimeError(
                                "NAS mount identity/options changed after ACK commit"
                            )
                        if read_regular_bytes_at(
                            ack_parent,
                            local_ack_name,
                            require_read_only=True,
                        ) != ack_bytes:
                            raise RuntimeError("local ACK readback differs")
                        local.assert_directory_visible((".acks",), ack_parent)
                        visible_ack_parent = local.open_directory((".acks",))
                        try:
                            if read_regular_bytes_at(
                                visible_ack_parent,
                                local_ack_name,
                                require_read_only=True,
                            ) != ack_bytes:
                                raise RuntimeError("visible local ACK readback differs")
                        finally:
                            os.close(visible_ack_parent)
                        verify_final_identities(local, final_proofs)
                        remote_ack_readback = ack_session.remote_call(
                            "read back immutable remote ACK",
                            lambda _client, current_sftp: retryable_remote_io(
                                lambda: read_remote_bytes(
                                    current_sftp,
                                    remote_ack,
                                    require_read_only=True,
                                )
                            ),
                        )
                        if remote_ack_readback != ack_bytes:
                            raise RuntimeError("remote ACK readback differs")
                        final_manifest_readback = ack_session.remote_call(
                            "final remote manifest readback",
                            lambda _client, current_sftp: retryable_remote_io(
                                lambda: read_remote_bytes(
                                    current_sftp,
                                    args.remote_manifest,
                                    require_read_only=True,
                                )
                            ),
                        )
                        if final_manifest_readback != manifest_bytes:
                            raise RuntimeError(
                                "remote manifest changed while publishing ACK"
                            )
                        print(
                            f"COMPLETE: {len(records)} files, "
                            f"{expected_core['total_bytes']} bytes, durable "
                            "local+remote immutable transfer ACK verified "
                            "(cleanup not authorized)",
                            flush=True,
                        )
                        return 0
                    finally:
                        os.close(ack_parent)
            finally:
                if sftp is not None:
                    sftp.close()
                if client is not None:
                    client.close()


if __name__ == "__main__":
    raise SystemExit(main())
