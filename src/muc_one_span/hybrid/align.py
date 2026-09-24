"""edlib-based alignment primitives shared by the hybrid engine stages.

edlib is imported lazily: ``rc`` and ``cigar_ops`` (and the synthetic test factory)
work without the optional ``hybrid`` extra, and alignment calls fail with a clear hint.
"""

from __future__ import annotations

import importlib
from dataclasses import dataclass
from types import ModuleType

HYBRID_HINT = "The hybrid engine needs the 'hybrid' extra: pip install 'muc_one_span[hybrid]'"
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _edlib() -> ModuleType:
    try:
        return importlib.import_module("edlib")
    except ImportError as exc:
        raise ImportError(HYBRID_HINT) from exc


def rc(seq: str) -> str:
    """Reverse complement."""
    return seq.translate(_COMP)[::-1]


def cigar_ops(cigar: str) -> list[tuple[int, str]]:
    """Parse an extended CIGAR (=, X, I, D) into (length, op) pairs."""
    ops: list[tuple[int, str]] = []
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            ops.append((int(num), ch))
            num = ""
    return ops


def infix_hit(query: str, target: str, k: int) -> tuple[int, int, int] | None:
    """Best infix location of query in target with at most k edits (start, end_excl, edits)."""
    res = _edlib().align(query, target, mode="HW", task="locations", k=k)
    if res["editDistance"] < 0:
        return None
    start, end = res["locations"][0]
    return int(start), int(end) + 1, int(res["editDistance"])


def edit_distance_infix(query: str, target: str) -> int:
    """Edit distance of query aligned fully inside target (semi-global)."""
    return int(_edlib().align(query, target, mode="HW", task="distance")["editDistance"])


def edit_distance(a: str, b: str) -> int:
    """Global edit distance."""
    return int(_edlib().align(a, b, mode="NW", task="distance")["editDistance"])


@dataclass
class Columns:
    """A read projected onto consensus columns.

    ``cols[i]`` is the read base, ``"-"`` for a deletion, or ``None`` where a partial
    read does not cover column ``i``. ``ins[i]`` holds read bases inserted before
    column ``i`` (``i == len(cons)`` at the end). ``t2q`` maps consensus positions
    (``len(cons) + 1`` entries) to read positions, ``-1`` where uncovered.
    """

    cols: list[str | None]
    ins: dict[int, str]
    t2q: list[int]
    first: int
    last: int  # covered consensus interval [first, last)

    def covers(self, start: int, end: int) -> bool:
        """True when the read covers [start - 1, end] (one anchoring column each side)."""
        return self.first <= max(start - 1, 0) and min(end + 1, len(self.cols)) <= self.last


def project(read: str, cons: str, *, partial: bool = False) -> Columns:
    """Project read onto consensus columns (global, or infix for partial reads)."""
    res = _edlib().align(read, cons, mode="HW" if partial else "NW", task="path")
    n = len(cons)
    first = int(res["locations"][0][0]) if partial else 0
    cols: list[str | None] = [None] * n
    t2q = [-1] * (n + 1)
    ins: dict[int, str] = {}
    ti, qi = first, 0
    for length, op in cigar_ops(res["cigar"]):
        if op in "=X":
            for _ in range(length):
                cols[ti] = read[qi]
                t2q[ti] = qi
                ti += 1
                qi += 1
        elif op == "I":
            ins[ti] = ins.get(ti, "") + read[qi : qi + length]
            qi += length
        else:
            for _ in range(length):
                cols[ti] = "-"
                t2q[ti] = qi
                ti += 1
    t2q[ti] = qi
    return Columns(cols, ins, t2q, first, ti)


def global_columns(read: str, cons: str) -> Columns:
    """Globally align a spanning read to cons and project onto consensus columns."""
    return project(read, cons)
