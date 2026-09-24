"""S3/S7 POA backend selection and polishing, including partial-read projection."""

from __future__ import annotations

import random
import sys

import pytest

from muc_one_span.hybrid.poa import get_backend
from muc_one_span.hybrid.polish import draft_consensus, homopolymer_vote, pileup_polish, polish
from muc_one_span.hybrid.spans import Anchors, SpanRead, categorize_reads
from muc_one_span.settings import DEFAULT_SETTINGS, HybridSettings
from tests.unit.hybrid import synth

pytest.importorskip("edlib", reason="edlib (extra 'hybrid') is not installed")
BACKEND = DEFAULT_SETTINGS.hybrid.poa_backend
DUPC_ALLELE = synth.allele(["X"] * 5 + [synth.dupc()] + ["X"] * 6 + ["A", "B"] + ["X"] * 10)


class _RecordingBackend:
    """A stub :class:`~muc_one_span.hybrid.poa.PoaBackend` that records what it was fed."""

    name = "recording"

    def __init__(self) -> None:
        self.seqs: list[str] = []

    def consensus(self, seqs: list[str]) -> str:
        self.seqs = seqs
        return seqs[0]


def _members(seq: str, n: int, seed: int) -> list[SpanRead]:
    return categorize_reads(
        synth.reads(seq, n, err=0.03, seed=seed),
        Anchors.from_dictionary(synth.RD, HybridSettings()),
        HybridSettings(),
    ).spanning


def test_poa_plus_polish_recovers_exact_allele_with_dupc() -> None:
    pytest.importorskip("pyabpoa", reason="pyabpoa (extra 'hybrid') is not installed")
    members = _members(DUPC_ALLELE, 60, 11)
    cons = draft_consensus(members, 40, random.Random(1), get_backend(BACKEND))
    cons, info = polish(cons, [m.seq for m in members], rounds=2, hp_vote=True)
    assert cons == DUPC_ALLELE, info


def test_pyspoa_backend_is_selectable() -> None:
    pytest.importorskip("spoa", reason="pyspoa (extra 'hybrid') is not installed")
    backend = get_backend("pyspoa")
    assert backend.name == "pyspoa"
    assert backend.consensus(["ACGTACGT", "ACGTACGT", "ACGAACGT"]) == "ACGTACGT"


def test_missing_backend_names_the_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "pyabpoa", None)
    with pytest.raises(ImportError, match="hybrid"):
        get_backend("pyabpoa")
    with pytest.raises(ValueError, match="unknown"):
        get_backend("medaka")


def test_partial_reads_do_not_truncate_the_consensus() -> None:
    truth = DUPC_ALLELE
    noisy = truth[:700] + truth[701:1400] + "G" + truth[1400:]  # one deletion, one insertion
    full = [r.seq for r in _members(truth, 20, 12)]
    half = len(truth) // 2
    fragments = [
        r.seq
        for r in synth.reads(truth[:half], 60, err=0.03, seed=13, strand_mix=False, flank_bp=0)
    ]
    cons = noisy
    for _ in range(2):
        cons, _changes = pileup_polish(cons, full, fragments)
        cons, _hp = homopolymer_vote(cons, full, fragments)
    assert cons == truth


def test_homopolymer_vote_uses_median_run_length() -> None:
    cons = "ACGTTGCA" + "C" * 7 + "AGTTGCAT"
    reads = [cons.replace("C" * 7, "C" * 8)] * 6 + [cons] * 3
    new, changes = homopolymer_vote(cons, reads)
    assert new == cons.replace("C" * 7, "C" * 8) and changes == 1


# --- Config-driven tunables (owner directive 2026-09-25): every literal below is a
# validated HybridSettings field; these tests prove the settings actually control
# behaviour, not just that the settings-sourced defaults reproduce the old literals.


def test_homopolymer_vote_min_run_setting_controls_which_runs_are_voted() -> None:
    cons = "ACGT" + "C" * 4 + "TTGCA"
    reads = [cons.replace("C" * 4, "C" * 6)] * 5 + [cons] * 2
    changed, changes = homopolymer_vote(cons, reads)
    assert changed == cons.replace("C" * 4, "C" * 6) and changes == 1
    unchanged, no_changes = homopolymer_vote(cons, reads, min_len=5)
    assert unchanged == cons and no_changes == 0


def test_pileup_polish_insertion_majority_frac_controls_acceptance() -> None:
    cons = "ACGTACGTACGT"
    ins_read = "ACGTACXGTACGT"  # one base inserted after consensus column 6
    full = [ins_read] * 6 + [cons] * 4  # 6/10 vote to insert
    accepted, changes = pileup_polish(cons, full)
    assert accepted == "ACGTACXGTACGT" and changes == 1
    rejected, no_changes = pileup_polish(cons, full, insertion_majority_frac=0.7)
    assert rejected == cons and no_changes == 0


def test_draft_consensus_sample_window_settings_control_selected_members() -> None:
    lengths = [100, 100, 100, 130, 100]  # one 30bp-longer outlier
    members = [
        SpanRead(f"r{i}", "A" * length, 30.0, "+", 0, "motif") for i, length in enumerate(lengths)
    ]
    narrow = _RecordingBackend()
    draft_consensus(
        members,
        10,
        random.Random(0),
        narrow,
        sample_window_floor_bp=5.0,
        sample_window_frac=0.0,
    )
    assert len(narrow.seqs) == 4  # outlier excluded: tolerance 5bp < 30bp distance

    wide = _RecordingBackend()
    draft_consensus(
        members,
        10,
        random.Random(0),
        wide,
        sample_window_floor_bp=40.0,
        sample_window_frac=0.0,
    )
    assert len(wide.seqs) == 5  # outlier included: tolerance 40bp >= 30bp distance


def test_polish_uses_settings_defaults_for_rounds_and_hp_vote() -> None:
    members = _members(DUPC_ALLELE, 20, 21)
    _cons, info = polish(DUPC_ALLELE, [m.seq for m in members])
    assert len(info["rounds"]) == DEFAULT_SETTINGS.hybrid.polish_rounds
