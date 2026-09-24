"""Append-only pre-registration of the decision rule, with an audit of test unsealing.

``<out_root>/test/preregistration.jsonl`` holds one JSON line per registration
(``sha256`` of the exact rule text, the text, ``registered_at``). Lines are only
ever appended, under an exclusive lock and fsync. The first scoring of ``test``
writes ``first_evaluation.json`` beside it (exclusive create, never rewritten);
from then on a new pre-registration is refused, so a rule cannot be registered
after the test truth has been seen. The marker is written *before* the scoring
(``evaluate``) or truth reading (``realism``) starts, so an evaluation that later
fails still counts as unsealing: conservative by design.
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

FIRST_EVALUATION = "first_evaluation.json"


def rule_sha256(rule_text: str) -> str:
    """SHA-256 of the exact rule text (UTF-8)."""
    return hashlib.sha256(rule_text.encode("utf-8")).hexdigest()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _marker(path: Path) -> Path:
    return Path(path).parent / FIRST_EVALUATION


def first_evaluation(path: Path) -> str | None:
    """When ``test`` was first scored under the ledger at ``path`` (``None`` if never)."""
    marker = _marker(path)
    if not marker.is_file():
        return None
    return str(json.loads(marker.read_text(encoding="utf-8"))["evaluated_at"])


def mark_first_evaluation(path: Path) -> str:
    """Record the first ``test`` scoring time once (exclusive create); return it."""
    marker = _marker(path)
    marker.parent.mkdir(parents=True, exist_ok=True)
    try:
        fd = os.open(marker, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    except FileExistsError:
        return str(first_evaluation(path))
    stamp = _now()
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        handle.write(json.dumps({"evaluated_at": stamp}) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    return stamp


def preregister(rule_text: str, path: Path) -> str:
    """Append (never rewrite) a pre-registration entry for ``rule_text``; return its sha256.

    Raises:
        PermissionError: If ``test`` has already been scored under this ledger.
    """
    path = Path(path)
    evaluated = first_evaluation(path)
    if evaluated is not None:
        raise PermissionError(f"test was already evaluated at {evaluated}; cannot pre-register")
    digest = rule_sha256(rule_text)
    path.parent.mkdir(parents=True, exist_ok=True)
    entry = {"sha256": digest, "rule_text": rule_text, "registered_at": _now()}
    with path.open("a", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            handle.write(json.dumps(entry, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    return digest


def require_preregistered(path: Path, rule_text: str) -> dict[str, Any]:
    """Return the first matching entry (``sha256``, ``registered_at``) or raise.

    Raises:
        PermissionError: If the ledger is missing or corrupt, or has no entry
            whose stored hash and text both match ``rule_text``.
    """
    digest = rule_sha256(rule_text)
    path = Path(path)
    if not path.is_file():
        raise PermissionError(f"no pre-registration ledger at {path}; run `benchsim preregister`")
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            entry = json.loads(line)
        except ValueError as exc:
            raise PermissionError(f"corrupt pre-registration ledger {path}: {line!r}") from exc
        if entry.get("sha256") == digest and rule_sha256(str(entry.get("rule_text"))) == digest:
            return {"sha256": digest, "registered_at": entry.get("registered_at")}
    raise PermissionError(f"decision rule sha256 {digest} is not pre-registered in {path}")
