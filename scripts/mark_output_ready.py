#!/usr/bin/env python3
"""Atomically publish a checksum manifest for a completed AthenaK output segment."""

from __future__ import annotations

import argparse
import errno
import hashlib
import json
import os
from pathlib import Path
import stat
import time


_CHUNK_SIZE = 8 * 1024 * 1024


def _identity(value: os.stat_result) -> tuple[int, int, int, int, int, int]:
    """Fields which must not change while an output is being certified."""

    return (
        value.st_dev,
        value.st_ino,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
        value.st_nlink,
    )


def _require_regular_single_link(value: os.stat_result, display: Path) -> None:
    if not stat.S_ISREG(value.st_mode):
        raise ValueError(f"output is not a regular file: {display}")
    if value.st_nlink != 1:
        raise ValueError(
            f"output must have exactly one hard link: {display} "
            f"(nlink={value.st_nlink})"
        )


def file_sha256_fd(descriptor: int) -> str:
    """Hash an already-open file, never reopening its pathname."""

    digest = hashlib.sha256()
    os.lseek(descriptor, 0, os.SEEK_SET)
    while True:
        chunk = os.read(descriptor, _CHUNK_SIZE)
        if not chunk:
            break
        digest.update(chunk)
    return digest.hexdigest()


def _relative_path(root: Path, value: str) -> Path:
    candidate = Path(value)
    if candidate.is_absolute():
        try:
            relative = candidate.relative_to(root)
        except ValueError as exc:
            raise ValueError(f"output is outside root: {candidate}") from exc
    else:
        relative = candidate

    # Walking with openat(2) below deliberately rejects every symlinked parent.
    # Rejecting '..' up front also makes RESOLVE_BENEATH semantics explicit even
    # on Python versions which do not expose openat2(2).
    if (
        not relative.parts
        or relative.is_absolute()
        or any(part in ("", ".", "..") for part in relative.parts)
    ):
        raise ValueError(f"invalid output path below root: {value}")
    return relative


def _open_parent(root_descriptor: int, relative: Path) -> tuple[int, str]:
    """Open and pin a relative path's parent without following any symlink."""

    nofollow = getattr(os, "O_NOFOLLOW", None)
    if nofollow is None:
        raise RuntimeError("O_NOFOLLOW is required to certify output files")
    directory_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | nofollow
    parent = os.dup(root_descriptor)
    try:
        for component in relative.parts[:-1]:
            child = os.open(component, directory_flags, dir_fd=parent)
            os.close(parent)
            parent = child
        return parent, relative.parts[-1]
    except BaseException:
        os.close(parent)
        raise


def _secure_lstat(root_descriptor: int, relative: Path) -> os.stat_result:
    parent, leaf = _open_parent(root_descriptor, relative)
    try:
        return os.stat(leaf, dir_fd=parent, follow_symlinks=False)
    finally:
        os.close(parent)


def _open_output(
    root_descriptor: int, root: Path, value: str
) -> tuple[int, int, str, Path, os.stat_result]:
    """Open one output through pinned directory descriptors.

    The returned file and parent descriptors remain open through hashing and the
    final pathname check.  This prevents a pathname replacement from changing
    which bytes are certified.
    """

    relative = _relative_path(root, value)
    display = root / relative
    parent, leaf = _open_parent(root_descriptor, relative)
    file_descriptor = -1
    try:
        path_stat = os.stat(leaf, dir_fd=parent, follow_symlinks=False)
        _require_regular_single_link(path_stat, display)
        nofollow = getattr(os, "O_NOFOLLOW", None)
        if nofollow is None:
            raise RuntimeError("O_NOFOLLOW is required to certify output files")
        file_descriptor = os.open(
            leaf,
            os.O_RDONLY | os.O_CLOEXEC | nofollow,
            dir_fd=parent,
        )
        descriptor_stat = os.fstat(file_descriptor)
        _require_regular_single_link(descriptor_stat, display)
        if _identity(path_stat) != _identity(descriptor_stat):
            raise ValueError(f"output pathname changed while opening: {display}")
        return file_descriptor, parent, leaf, relative, descriptor_stat
    except BaseException:
        if file_descriptor >= 0:
            os.close(file_descriptor)
        os.close(parent)
        raise


def _process_exists(process: Path) -> bool:
    try:
        process.stat()
    except FileNotFoundError:
        return False
    except (PermissionError, ProcessLookupError):
        # A permission error is evidence that the process still exists but is not
        # inspectable.  Callers must fail closed in that case.
        return True
    return True


def _proc_mount_options(proc: Path) -> set[str]:
    """Return options for the /proc mount, refusing an uninspectable mount table."""

    if proc != Path("/proc"):
        return set()
    try:
        lines = (proc / "self" / "mountinfo").read_text(
            encoding="utf-8", errors="strict"
        ).splitlines()
    except (OSError, UnicodeError) as exc:
        raise RuntimeError(f"cannot inspect /proc mount visibility: {exc}") from exc

    matches: list[set[str]] = []
    for line in lines:
        fields = line.split()
        try:
            separator = fields.index("-")
        except ValueError:
            continue
        if len(fields) <= separator + 3 or fields[4] != "/proc":
            continue
        options = set(fields[5].split(","))
        options.update(fields[separator + 3].split(","))
        matches.append(options)
    if not matches:
        raise RuntimeError("cannot identify the /proc mount")
    return set().union(*matches)


def _assert_proc_visibility(proc: Path) -> None:
    if not proc.is_dir():
        raise RuntimeError("/proc is unavailable")
    options = _proc_mount_options(proc)
    for option in options:
        if option == "hidepid" or (
            option.startswith("hidepid=")
            and option.split("=", 1)[1] not in ("0", "off", "no")
        ):
            raise RuntimeError(f"/proc process visibility is restricted ({option})")


def _read_process_maps(process: Path) -> list[str] | None:
    """Read maps, returning None only when the process demonstrably vanished."""

    try:
        return (process / "maps").read_text(
            encoding="utf-8", errors="surrogateescape"
        ).splitlines()
    except (FileNotFoundError, ProcessLookupError) as exc:
        if not _process_exists(process):
            return None
        raise RuntimeError(f"cannot inspect {process}/maps: {exc}") from exc
    except PermissionError as exc:
        if not _process_exists(process):
            return None
        raise RuntimeError(f"permission denied inspecting {process}/maps") from exc
    except OSError as exc:
        if exc.errno == errno.ESRCH and not _process_exists(process):
            return None
        raise RuntimeError(f"cannot inspect {process}/maps: {exc}") from exc


def _has_writable_mapping(lines: list[str], target: os.stat_result) -> bool:
    target_device = (os.major(target.st_dev), os.minor(target.st_dev))
    for line in lines:
        fields = line.split(None, 5)
        if len(fields) < 5:
            raise RuntimeError(f"cannot parse /proc maps entry: {line!r}")
        permissions = fields[1]
        device = fields[3]
        inode_text = fields[4]
        try:
            major_text, minor_text = device.split(":", 1)
            mapped_device = (int(major_text, 16), int(minor_text, 16))
            mapped_inode = int(inode_text, 10)
        except (ValueError, TypeError) as exc:
            raise RuntimeError(f"cannot parse /proc maps entry: {line!r}") from exc
        if (
            "w" in permissions
            and mapped_inode == target.st_ino
            and mapped_device == target_device
        ):
            return True
    return False


def processes_holding(
    target: os.stat_result, proc: Path = Path("/proc")
) -> list[tuple[int, str, str]]:
    """Return processes with a descriptor or writable mmap for *target*.

    Any process-visibility or permission gap fails closed.  The only ignored
    inspection race is a process which can be shown to have disappeared.
    """

    _assert_proc_visibility(proc)
    try:
        processes = list(proc.iterdir())
    except (OSError, PermissionError) as exc:
        raise RuntimeError(f"cannot enumerate /proc: {exc}") from exc

    holders: list[tuple[int, str, str]] = []
    for process in processes:
        if not process.name.isdigit() or int(process.name) == os.getpid():
            continue
        pid = int(process.name)
        descriptors = process / "fd"
        try:
            entries = list(descriptors.iterdir())
        except (FileNotFoundError, ProcessLookupError) as exc:
            if not _process_exists(process):
                continue
            raise RuntimeError(f"cannot inspect {descriptors}: {exc}") from exc
        except PermissionError as exc:
            if not _process_exists(process):
                continue
            raise RuntimeError(f"permission denied inspecting {descriptors}") from exc
        except OSError as exc:
            if exc.errno == errno.ESRCH and not _process_exists(process):
                continue
            raise RuntimeError(f"cannot inspect {descriptors}: {exc}") from exc

        held_by_descriptor = False
        for descriptor in entries:
            try:
                descriptor_stat = descriptor.stat()
            except (FileNotFoundError, ProcessLookupError):
                # Descriptors may close independently while the process remains.
                continue
            except PermissionError as exc:
                if not _process_exists(process):
                    break
                raise RuntimeError(
                    f"permission denied inspecting {descriptor}"
                ) from exc
            except OSError as exc:
                if exc.errno in (errno.ENOENT, errno.ESRCH):
                    continue
                raise RuntimeError(f"cannot inspect {descriptor}: {exc}") from exc
            if (descriptor_stat.st_dev, descriptor_stat.st_ino) == (
                target.st_dev,
                target.st_ino,
            ):
                held_by_descriptor = True
                break

        maps = _read_process_maps(process)
        if maps is None:
            continue
        held_by_mapping = _has_writable_mapping(maps, target)
        if not held_by_descriptor and not held_by_mapping:
            continue
        try:
            command = (process / "comm").read_text(encoding="utf-8").strip()
        except (FileNotFoundError, PermissionError, ProcessLookupError, UnicodeError):
            command = "unknown"
        reason = "descriptor+writable-mmap" if (
            held_by_descriptor and held_by_mapping
        ) else ("descriptor" if held_by_descriptor else "writable-mmap")
        holders.append((pid, command, reason))
    return sorted(holders)


def assert_closed(target: os.stat_result, display: Path) -> None:
    try:
        holders = processes_holding(target)
    except RuntimeError as exc:
        raise SystemExit(f"cannot prove output closure for {display}: {exc}") from exc
    if holders:
        detail = ", ".join(
            f"pid={pid}({command},{reason})"
            for pid, command, reason in holders
        )
        raise SystemExit(f"refusing open output file {display}: {detail}")


def _assert_root_identity(root: Path, root_descriptor: int) -> None:
    try:
        pathname_stat = os.stat(root, follow_symlinks=False)
    except (FileNotFoundError, OSError) as exc:
        raise SystemExit(f"output root is no longer reachable: {root}: {exc}") from exc
    descriptor_stat = os.fstat(root_descriptor)
    if not stat.S_ISDIR(pathname_stat.st_mode) or (
        pathname_stat.st_dev,
        pathname_stat.st_ino,
    ) != (descriptor_stat.st_dev, descriptor_stat.st_ino):
        raise SystemExit(f"output root pathname changed during publication: {root}")


def _write_all(descriptor: int, payload: bytes) -> None:
    offset = 0
    while offset < len(payload):
        offset += os.write(descriptor, payload[offset:])


def _publish_manifest(manifest_dir: Path, segment: str, payload: dict) -> Path:
    manifest_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir = manifest_dir.resolve(strict=True)
    nofollow = getattr(os, "O_NOFOLLOW", None)
    if nofollow is None:
        raise SystemExit("O_NOFOLLOW is required to publish a ready manifest")
    directory_descriptor = os.open(
        manifest_dir,
        os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | nofollow,
    )
    destination_name = f"{segment}.manifest.ready"
    temporary_name = f".{segment}.{os.getpid()}.{time.time_ns()}.tmp"
    temporary_created = False
    try:
        try:
            descriptor = os.open(
                temporary_name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | nofollow,
                0o600,
                dir_fd=directory_descriptor,
            )
        except FileExistsError as exc:
            raise SystemExit(
                f"refusing stale temporary manifest: {manifest_dir / temporary_name}"
            ) from exc
        temporary_created = True
        try:
            serialized = (
                json.dumps(payload, indent=2, sort_keys=True) + "\n"
            ).encode("utf-8")
            _write_all(descriptor, serialized)
            if hasattr(os, "fdatasync"):
                os.fdatasync(descriptor)
            else:
                os.fsync(descriptor)
            # The ready inode is non-writable before it becomes visible.  fsync,
            # rather than fdatasync alone, also persists this mode transition.
            os.fchmod(descriptor, 0o444)
            os.fsync(descriptor)
        finally:
            os.close(descriptor)

        try:
            os.link(
                temporary_name,
                destination_name,
                src_dir_fd=directory_descriptor,
                dst_dir_fd=directory_descriptor,
                follow_symlinks=False,
            )
        except FileExistsError as exc:
            raise SystemExit(
                "refusing to replace existing ready manifest: "
                f"{manifest_dir / destination_name}"
            ) from exc
        # Persist the create-no-replace ready name before dropping the private
        # temporary link, then persist removal of that link as a second step.
        os.fsync(directory_descriptor)
        os.unlink(temporary_name, dir_fd=directory_descriptor)
        temporary_created = False
        os.fsync(directory_descriptor)
        return manifest_dir / destination_name
    finally:
        if temporary_created:
            try:
                os.unlink(temporary_name, dir_fd=directory_descriptor)
                os.fsync(directory_descriptor)
            except FileNotFoundError:
                pass
        os.close(directory_descriptor)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manifest-dir", required=True, type=Path)
    parser.add_argument("--segment", required=True)
    parser.add_argument(
        "--min-age",
        type=float,
        default=120.0,
        help="refuse files modified more recently than this many seconds",
    )
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()

    root = args.root.expanduser().resolve(strict=True)
    valid_segment_characters = (
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_."
    )
    if not args.segment or any(
        char not in valid_segment_characters for char in args.segment
    ):
        parser.error("segment may contain only letters, digits, '-', '_' and '.'")

    nofollow = getattr(os, "O_NOFOLLOW", None)
    if nofollow is None:
        raise SystemExit("O_NOFOLLOW is required to certify output files")
    root_descriptor = os.open(
        root,
        os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | nofollow,
    )
    try:
        _assert_root_identity(root, root_descriptor)
        now_ns = time.time_ns()
        records = []
        seen: set[str] = set()
        for value in args.files:
            try:
                descriptor, parent, leaf, relative, stat_before = _open_output(
                    root_descriptor, root, value
                )
            except (OSError, ValueError, RuntimeError) as exc:
                raise SystemExit(str(exc)) from exc
            try:
                relative_text = relative.as_posix()
                if relative_text in seen:
                    continue
                age = (now_ns - stat_before.st_mtime_ns) / 1_000_000_000
                if age < args.min_age:
                    raise SystemExit(
                        f"refusing possibly open file {relative_text}: age {age:.1f}s "
                        f"< {args.min_age}s"
                    )

                display = root / relative
                assert_closed(stat_before, display)
                checksum = file_sha256_fd(descriptor)
                stat_after = os.fstat(descriptor)
                _require_regular_single_link(stat_after, display)
                if _identity(stat_before) != _identity(stat_after):
                    raise SystemExit(f"file changed while hashing: {relative_text}")

                # Check both the pinned original parent entry and a fresh secure
                # traversal from the pinned root.  Either catches a rename/swap.
                pinned_path_stat = os.stat(
                    leaf, dir_fd=parent, follow_symlinks=False
                )
                visible_path_stat = _secure_lstat(root_descriptor, relative)
                for pathname_stat in (pinned_path_stat, visible_path_stat):
                    _require_regular_single_link(pathname_stat, display)
                    if _identity(pathname_stat) != _identity(stat_after):
                        raise SystemExit(
                            f"output pathname changed while hashing: {relative_text}"
                        )

                assert_closed(stat_after, display)
                stat_final = os.fstat(descriptor)
                final_path_stat = _secure_lstat(root_descriptor, relative)
                _require_regular_single_link(stat_final, display)
                _require_regular_single_link(final_path_stat, display)
                if (
                    _identity(stat_after) != _identity(stat_final)
                    or _identity(stat_final) != _identity(final_path_stat)
                ):
                    raise SystemExit(
                        f"file or pathname changed during closure check: {relative_text}"
                    )

                records.append(
                    {
                        "path": relative_text,
                        "size": stat_final.st_size,
                        "sha256": checksum,
                    }
                )
                seen.add(relative_text)
            except ValueError as exc:
                raise SystemExit(str(exc)) from exc
            finally:
                os.close(descriptor)
                os.close(parent)

        _assert_root_identity(root, root_descriptor)
    finally:
        os.close(root_descriptor)

    payload = {
        "schema": 1,
        "segment": args.segment,
        "created_unix": time.time(),
        "root": str(root),
        "files": records,
    }
    destination = _publish_manifest(
        args.manifest_dir.expanduser(), args.segment, payload
    )
    print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
