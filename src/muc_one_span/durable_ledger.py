"""Crash-consistent, transactional ledger management for MucOneUp experiments."""

from __future__ import annotations

import fcntl
import hashlib
import json
import logging
import os
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def compute_sha256(path: Path) -> str:
    """Compute SHA-256 hex digest of a file in 1MB chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fsync_dir(dir_path: Path) -> None:
    """Fsync parent directory to ensure directory entry durability."""
    dir_fd = os.open(str(dir_path), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(dir_fd)
    finally:
        os.close(dir_fd)


@dataclass
class LedgerEntry:
    """Complete sealed ledger entry with full provenance."""

    design_id: int
    design_name: str
    token: str
    split: str
    category: str
    platform: str
    lengths: list[int]
    mutation: str | None
    targets: list[list[int]]
    truth_fa: str
    truth_sha256: str
    reads_file: str
    reads_sha256: str
    usable_records: int
    git_sha: str | None = None
    config_sha256: str | None = None
    status: str = "completed"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class DurableLedger:
    """Thread-safe, crash-consistent authoritative ledger manager with file locking."""

    def __init__(self, output_dir: Path) -> None:
        self.output_dir = output_dir.resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.lock_path = self.output_dir / ".ledger.lock"
        self.sealed_ledger_path = self.output_dir / "ledger_sealed.jsonl"
        self.public_ledger_path = self.output_dir / "ledger_public.jsonl"
        self._entries: dict[tuple[str, str], LedgerEntry] = {}
        self._load_existing()

    @contextmanager
    def file_lock(self) -> Generator[None, None, None]:
        """Acquire exclusive file lock for the ledger directory."""
        with self.lock_path.open("w") as lock_file:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)

    @contextmanager
    def run_lock(self) -> Generator[None, None, None]:
        """Acquire non-blocking run lock for the output directory to prevent overlapping runs."""
        run_lock_path = self.output_dir / ".run.lock"
        with run_lock_path.open("w") as lock_file:
            try:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except (BlockingIOError, OSError) as exc:
                raise RuntimeError(
                    f"Another process is currently running on {self.output_dir}"
                ) from exc
            try:
                yield
            finally:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)

    def _load_existing(self) -> None:
        """Load sealed ledger if present, failing closed on corrupt lines."""
        if not self.sealed_ledger_path.exists():
            return
        for line in self.sealed_ledger_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                data = json.loads(line)
                key = (str(data["design_name"]), str(data["platform"]))
                entry = LedgerEntry(
                    design_id=int(data["design_id"]),
                    design_name=str(data["design_name"]),
                    token=str(data["token"]),
                    split=str(data["split"]),
                    category=str(data["category"]),
                    platform=str(data["platform"]),
                    lengths=list(data["lengths"]),
                    mutation=data.get("mutation"),
                    targets=[list(t) for t in data.get("targets", [])],
                    truth_fa=str(data["truth_fa"]),
                    truth_sha256=str(data["truth_sha256"]),
                    reads_file=str(data["reads_file"]),
                    reads_sha256=str(data["reads_sha256"]),
                    usable_records=int(data["usable_records"]),
                    git_sha=data.get("git_sha"),
                    config_sha256=data.get("config_sha256"),
                    status=str(data.get("status", "completed")),
                )
                self._entries[key] = entry
            except (json.JSONDecodeError, KeyError, ValueError) as exc:
                raise RuntimeError(
                    f"Corrupted ledger entry in {self.sealed_ledger_path}: {line}"
                ) from exc

        # Reconcile public ledger if out of sync
        if self._entries and (
            not self.public_ledger_path.exists()
            or len(self.public_ledger_path.read_text(encoding="utf-8").splitlines())
            != len(self._entries)
        ):
            self._atomic_export()

    def is_verified_complete(self, design_name: str, platform: str) -> bool:
        """Check whether a sample already completed with valid artifacts on disk."""
        key = (design_name, platform)
        if key not in self._entries:
            return False
        entry = self._entries[key]
        if entry.status != "completed":
            return False

        truth_path = self.output_dir / entry.truth_fa
        reads_path = self.output_dir / entry.reads_file

        if not truth_path.exists() or truth_path.stat().st_size == 0:
            return False
        if not reads_path.exists() or reads_path.stat().st_size == 0:
            return False

        # Verify recorded SHA256 if present
        if entry.truth_sha256 and compute_sha256(truth_path) != entry.truth_sha256:
            return False
        return not (entry.reads_sha256 and compute_sha256(reads_path) != entry.reads_sha256)

    def commit_entry(self, entry: LedgerEntry) -> None:
        """Add or update an entry and atomically persist sealed and public ledgers."""
        with self.file_lock():
            self._load_existing()
            key = (entry.design_name, entry.platform)
            self._entries[key] = entry
            self._atomic_export()

    def _atomic_export(self) -> None:
        """Atomically rewrite both sealed and public ledgers on disk."""
        tmp_sealed = self.sealed_ledger_path.with_suffix(".tmp")
        tmp_public = self.public_ledger_path.with_suffix(".tmp")

        sorted_keys = sorted(self._entries.keys())

        # Write sealed ledger
        with tmp_sealed.open("w", encoding="utf-8") as f:
            for k in sorted_keys:
                ent = self._entries[k]
                f.write(json.dumps(ent.to_dict()) + "\n")
            f.flush()
            os.fsync(f.fileno())

        # Write public ledger (blinded token paths only, no truth/design leakage)
        with tmp_public.open("w", encoding="utf-8") as f:
            for k in sorted_keys:
                ent = self._entries[k]
                ext = Path(ent.reads_file).suffix
                pub = {
                    "token": ent.token,
                    "platform": ent.platform,
                    "reads_file": f"blinded_reads/{ent.token}/{ent.platform}/{ent.token}_{ent.platform}{ext}",
                    "usable_records": ent.usable_records,
                }
                f.write(json.dumps(pub) + "\n")
            f.flush()
            os.fsync(f.fileno())

        # Atomic replace
        tmp_sealed.replace(self.sealed_ledger_path)
        tmp_public.replace(self.public_ledger_path)
        fsync_dir(self.output_dir)

    def total_entries(self) -> int:
        return len(self._entries)

    def get_entries(self) -> list[LedgerEntry]:
        return [self._entries[k] for k in sorted(self._entries.keys())]
