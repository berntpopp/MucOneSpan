"""``infix_hit`` alignment primitive: normal hits and edlib's None-start quirk.

edlib's fast ``task="locations"`` pass can report a non-negative ``editDistance``
with an unresolved (``None``) start location -- observed in production when a
fragment read's matched anchor consumes the search window down to an empty
target substring (see ``test_spans.py`` for the genomic-fragment regression
that reproduces this with real edlib). These tests pin ``infix_hit``'s handling
of that quirk directly, with a stub edlib so the behaviour does not depend on
the installed edlib version reproducing the same internal edge case.
"""

from __future__ import annotations

from typing import Any

import pytest

from muc_one_span.hybrid import align


class _StubEdlib:
    """Fake ``edlib`` module: records the ``task`` of each call, replays canned results."""

    def __init__(self, locations_result: dict[str, Any], path_result: dict[str, Any]) -> None:
        self.calls: list[str] = []
        self._results = {"locations": locations_result, "path": path_result}

    def align(
        self, query: str, target: str, mode: str = "NW", task: str = "distance", k: int = -1
    ) -> dict[str, Any]:
        self.calls.append(task)
        return self._results[task]


def test_infix_hit_finds_a_normal_match() -> None:
    assert align.infix_hit("ACGT", "NNNACGTNNN", 0) == (3, 7, 0)


def test_infix_hit_returns_none_when_no_match_within_k() -> None:
    assert align.infix_hit("ACGTACGT", "TTTTTTTT", 1) is None


def test_infix_hit_recovers_via_path_retry_when_locations_start_is_none(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A ``locations`` pass with a non-negative distance but ``start=None`` must not be
    trusted as a hit; ``infix_hit`` retries with ``task="path"`` and, when that resolves
    a real location, returns it.
    """
    stub = _StubEdlib(
        locations_result={"editDistance": 2, "locations": [(None, -1)], "cigar": None},
        path_result={"editDistance": 2, "locations": [(3, 8)], "cigar": "6="},
    )
    monkeypatch.setattr(align, "_edlib", lambda: stub)
    assert align.infix_hit("ACGTAC", "NNNACGTACNN", 2) == (3, 9, 2)
    assert stub.calls == ["locations", "path"]


def test_infix_hit_returns_none_when_location_stays_unresolved_after_retry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When the ``path`` retry also cannot resolve a start (the empty-target-window
    case observed in production), there is no valid hit -- ``infix_hit`` must return
    ``None`` rather than raise or return a partial/invalid location.
    """
    stub = _StubEdlib(
        locations_result={"editDistance": 30, "locations": [(None, -1)], "cigar": None},
        path_result={"editDistance": 30, "locations": [(None, -1)], "cigar": None},
    )
    monkeypatch.setattr(align, "_edlib", lambda: stub)
    assert align.infix_hit("ACGTACGTACGTACGTACGTACGTACGTAC", "", 3) is None
    assert stub.calls == ["locations", "path"]


def test_infix_hit_returns_none_when_locations_list_is_empty(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Defensive: an empty ``locations`` list (no entry at all) is also not a hit."""
    stub = _StubEdlib(
        locations_result={"editDistance": 0, "locations": [], "cigar": None},
        path_result={"editDistance": 0, "locations": [], "cigar": None},
    )
    monkeypatch.setattr(align, "_edlib", lambda: stub)
    assert align.infix_hit("ACGT", "ACGT", 0) is None
    assert stub.calls == ["locations", "path"]
