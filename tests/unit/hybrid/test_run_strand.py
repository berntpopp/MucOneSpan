"""Stutter-aware run-site rules of the phase site table (Task 15e), on exact site tables.

The tables are built by hand (no alignment), so every read count is exact: ten peer C
runs of the consensus length show the stutter background, and one target C run
carries a minor length at a chosen share.
"""

from __future__ import annotations

from typing import Any

import pytest

from muc_one_span.hybrid.phase_sites import candidates
from muc_one_span.hybrid.run_strand import run_mixture
from muc_one_span.settings import HybridSettings

S = HybridSettings()
LENGTH = 6  # consensus length of every run in the table (a C6 run)
PEERS = [("run", 100 * k) for k in range(1, 11)]
TARGET = ("run", 2000)
META = dict.fromkeys([*PEERS, TARGET], ("C", LENGTH))
PERIOD = 100  # read-pattern period: shares below are in reads per PERIOD
STUTTER = 5  # per PERIOD reads, a peer run is read one base long, and as many one short
FLOOR = round(S.het_af_min * PERIOD)
# Observed just above het_af_min, but below it once the insertion stutter is removed.
LOW = FLOOR + 1
# Clearly above het_af_min after deconvolution.
HIGH = FLOOR + 2 * STUTTER


def _peer(j: int) -> int:
    k = j % PERIOD
    return LENGTH + 1 if k < STUTTER else LENGTH - 1 if k < 2 * STUTTER else LENGTH


def _table(
    minor_share: dict[str, int], n_reads: dict[str, int]
) -> tuple[list[dict[Any, Any]], list[str]]:
    """Reads per strand; the target shows LENGTH + 1 in ``minor_share`` of each PERIOD."""
    feats, strands = [], []
    for strand, n in n_reads.items():
        for j in range(n):
            f: dict[Any, Any] = {site: _peer(j + 7 * i) for i, site in enumerate(PEERS)}
            k = j % PERIOD
            share = minor_share[strand]
            f[TARGET] = LENGTH + 1 if k < share else LENGTH - 1 if k < share + STUTTER else LENGTH
            feats.append(f)
            strands.append(strand)
    return feats, strands


def _target_kept(feats: list[dict[Any, Any]], strands: list[str]) -> bool:
    return any(c["site"] == TARGET for c in candidates(feats, strands, META, S))


def test_minor_share_is_deconvolved_from_stutter() -> None:
    """Observed above het_af_min, but below it once the insertion stutter is removed."""
    both = {"+": 2 * PERIOD, "-": 2 * PERIOD}
    low = _table({"+": LOW, "-": LOW}, both)
    weight, p = run_mixture(*low, META, TARGET, (LENGTH, LENGTH + 1), S)
    assert weight < S.het_af_min < LOW / PERIOD and p == pytest.approx(1.0)
    assert not _target_kept(*low)
    assert _target_kept(*_table({"+": HIGH, "-": HIGH}, both))


def test_minor_absent_from_a_small_strand_is_rejected() -> None:
    """Too few "-" reads for a significant likelihood-ratio test, but none show the minor."""
    n_minus = S.hp_min_strand_reads
    feats, strands = _table({"+": HIGH, "-": 0}, {"+": 2 * PERIOD, "-": n_minus})
    _weight, p = run_mixture(feats, strands, META, TARGET, (LENGTH, LENGTH + 1), S)
    assert p >= S.phase_strand_bias_alpha  # the likelihood-ratio test alone would keep it
    assert not _target_kept(feats, strands)


def test_run_without_peers_uses_the_column_rules() -> None:
    """A run whose base forms no other run has no stutter model: counts are tested as is."""
    lone = {**META, TARGET: ("A", LENGTH)}
    feats, strands = _table({"+": HIGH, "-": HIGH}, {"+": 2 * PERIOD, "-": 2 * PERIOD})
    assert any(c["site"] == TARGET for c in candidates(feats, strands, lone, S))
