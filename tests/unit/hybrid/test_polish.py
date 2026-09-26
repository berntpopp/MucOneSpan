"""S3/S7 POA backend selection and polishing, including partial-read projection."""

from __future__ import annotations

import random
import sys

import pytest

from muc_one_span.hybrid.poa import get_backend
from muc_one_span.hybrid.polish import (
    consensus_concordance,
    draft_consensus,
    homopolymer_vote,
    pileup_polish,
    polish,
)
from muc_one_span.hybrid.spans import Anchors, SpanRead, categorize_reads
from muc_one_span.settings import DEFAULT_SETTINGS, HybridSettings
from tests.unit.hybrid import synth
from tests.unit.hybrid.synth import LAYOUT

BACKEND = DEFAULT_SETTINGS.hybrid.poa_backend
S = DEFAULT_SETTINGS.hybrid
POLISH = {
    "rounds": S.polish_rounds,
    "hp_vote": S.hp_vote,
    "insertion_majority_frac": S.polish_insertion_majority_frac,
    "hp_min_run": S.hp_vote_min_run,
}
WINDOW = {
    "sample_window_floor_bp": S.poa_sample_window_floor_bp,
    "sample_window_frac": S.poa_sample_window_frac,
}
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
        Anchors.from_dictionary(synth.RD, HybridSettings(), LAYOUT),
        HybridSettings(),
    ).spanning


def test_poa_plus_polish_recovers_exact_allele_with_dupc() -> None:
    members = _members(DUPC_ALLELE, 60, 11)
    cons = draft_consensus(members, S.n_poa, random.Random(S.seed), get_backend(BACKEND), **WINDOW)
    cons, info = polish(cons, [m.seq for m in members], **POLISH)
    assert cons == DUPC_ALLELE, info


def test_polish_corrects_a_substitution_and_a_homopolymer_error_given_realistic_reads() -> None:
    """polish() must not be a no-op: given a deliberately corrupted draft (one
    substitution and one shortened homopolymer run, injected far from the dupC
    template so they cannot be confused with it) and error-realistic reads generated
    from the true allele (the same noise model used throughout this suite), polish()
    must recover the true allele exactly and its round log must show real work done.

    This is the direct test for benchmark "Concerns" item 3 (polish changed only ~3%
    of rounds on real data): it distinguishes a genuinely low change rate driven by
    already-good POA drafts from a silently broken/no-op polish step. ``hp_changes``
    can legitimately stay 0 here -- the ordinary pileup vote (substitution/column and
    insertion-slot voting) already resolves this single-base homopolymer shortening on
    its own at this coverage and error rate; ``homopolymer_vote`` earns its keep on
    harder, more ambiguous runs, not every homopolymer edit.
    """
    truth = DUPC_ALLELE
    sub_pos = 50  # inside the leading PRE repeats, far from the dupC unit (~600-661bp)
    hp_start, hp_end = 559, 563  # a plain "C"*4 run, also well clear of the dupC unit
    assert truth[hp_start:hp_end] == "CCCC"
    new_base = "A" if truth[sub_pos] != "A" else "T"
    corrupted = truth[:sub_pos] + new_base + truth[sub_pos + 1 : hp_start] + "CCC" + truth[hp_end:]
    assert corrupted != truth

    members = _members(truth, 60, 31)
    cons, info = polish(corrupted, [m.seq for m in members], **POLISH)

    assert cons == truth, info
    assert sum(r["changes"] for r in info["rounds"]) > 0  # not a no-op


def test_pyspoa_backend_is_selectable() -> None:
    pytest.importorskip("spoa", reason="pyspoa (extra 'hybrid') is not installed")
    backend = get_backend("pyspoa")
    assert backend.name == "pyspoa"
    assert backend.consensus(["ACGTACGT", "ACGTACGT", "ACGAACGT"]) == "ACGTACGT"


def test_missing_pyspoa_names_the_hybrid_extra(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "spoa", None)
    with pytest.raises(ImportError, match=r"muc_one_span\[hybrid\]"):
        get_backend("pyspoa")


def test_unknown_backend_is_rejected() -> None:
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
        cons, _changes = pileup_polish(
            cons, full, fragments, insertion_majority_frac=S.polish_insertion_majority_frac
        )
        cons, _hp = homopolymer_vote(cons, full, fragments, min_len=S.hp_vote_min_run)
    assert cons == truth


def test_homopolymer_vote_uses_median_run_length() -> None:
    cons = "ACGTTGCA" + "C" * 7 + "AGTTGCAT"
    reads = [cons.replace("C" * 7, "C" * 8)] * 6 + [cons] * 3
    new, changes = homopolymer_vote(cons, reads, min_len=S.hp_vote_min_run)
    assert new == cons.replace("C" * 7, "C" * 8) and changes == 1


# --- Config-driven tunables (owner directive 2026-09-25): every literal below is a
# validated HybridSettings field; these tests prove the settings actually control
# behaviour, not just that the settings-sourced defaults reproduce the old literals.


def test_homopolymer_vote_min_run_setting_controls_which_runs_are_voted() -> None:
    cons = "ACGT" + "C" * 4 + "TTGCA"
    reads = [cons.replace("C" * 4, "C" * 6)] * 5 + [cons] * 2
    changed, changes = homopolymer_vote(cons, reads, min_len=S.hp_vote_min_run)
    assert changed == cons.replace("C" * 4, "C" * 6) and changes == 1
    unchanged, no_changes = homopolymer_vote(cons, reads, min_len=5)
    assert unchanged == cons and no_changes == 0


def test_pileup_polish_insertion_majority_frac_controls_acceptance() -> None:
    cons = "ACGTACGTACGT"
    ins_read = "ACGTACXGTACGT"  # one base inserted after consensus column 6
    full = [ins_read] * 6 + [cons] * 4  # 6/10 vote to insert
    accepted, changes = pileup_polish(
        cons, full, insertion_majority_frac=S.polish_insertion_majority_frac
    )
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


# --- consensus_concordance: real hybrid read-support evidence for the per-allele
# "confidence" concern (the ladder's classify.py confidence is a dictionary-fit
# heuristic fed the same way for both engines and says nothing about how well hybrid's
# own reads support the consensus it built; this is the meaningful alternative). A
# per-position *mean* concordance is used rather than a per-unit unanimity requirement:
# at real hybrid read depths (dozens to thousands of spanning reads), requiring every
# single covering read to agree on every base of a unit collapses to ~0 almost
# everywhere once realistic per-base noise is present, however good the consensus
# actually is -- confirmed on real reconstructed alleles during this fix.


def test_consensus_concordance_is_one_when_every_read_matches_consensus_exactly() -> None:
    cons = "ACGTAC" * 3
    reads = [cons] * 8
    assert consensus_concordance(cons, reads, None) == 1.0


def test_consensus_concordance_is_the_mean_per_base_agreement_rate() -> None:
    unit = 10
    cons = "A" * unit + "C" * unit
    agreeing = "A" * unit + "C" * unit
    disagreeing = "A" * unit + "G" * unit  # mismatches the consensus in the second half
    full = [agreeing] * 5 + [disagreeing] * 1  # 5/6 reads agree wherever they disagree
    # First half: 6/6 agree at every position (score 1.0). Second half: 5/6 agree at
    # every position (score 5/6). Mean over all 20 positions:
    expected = (unit * 1.0 + unit * (5 / 6)) / (2 * unit)
    assert consensus_concordance(cons, full, None) == pytest.approx(expected)


def test_consensus_concordance_treats_an_uncovered_region_as_zero_there() -> None:
    unit = 10
    cons = "A" * unit + "C" * unit
    partial = ["A" * unit]  # an infix/partial read covering only the first half
    # First half: fully covered and unanimous (score 1.0 at each position). Second
    # half: no coverage at all, so it contributes 0 (fail-closed) -- mean is 0.5.
    assert consensus_concordance(cons, [], partial) == 0.5


def test_consensus_concordance_is_zero_with_no_covering_reads() -> None:
    cons = "A" * 10 + "C" * 10
    assert consensus_concordance(cons, [], None) == 0.0


def test_consensus_concordance_fails_closed_for_an_empty_consensus() -> None:
    """An empty consensus has nothing for reads to support; the function's own
    docstring says a position with no covering read fails closed to 0, and an empty
    consensus is that same "no evidence" case in the limit, so it must not instead
    report a misleadingly perfect 1.0.
    """
    assert consensus_concordance("", ["ACGT"], None) == 0.0
    assert consensus_concordance("", [], None) == 0.0


def test_polish_rounds_come_from_the_caller() -> None:
    members = _members(DUPC_ALLELE, 20, 21)
    _cons, info = polish(DUPC_ALLELE, [m.seq for m in members], **POLISH)
    assert len(info["rounds"]) == S.polish_rounds
    _cons, more = polish(DUPC_ALLELE, [m.seq for m in members], **{**POLISH, "rounds": 3})
    assert len(more["rounds"]) == 3


def test_homopolymer_vote_ignores_reads_whose_run_boundary_base_differs() -> None:
    """Task 13b root cause (M3_hifi): a separator base lost or substituted in a read
    merges two runs (unit I ``GCG C5 A`` -> ``G C7 A``); that merged length is not a
    run-length observation and must not vote, or the median lands on a length present
    in neither read population (C6) and a correct consensus is rewritten wrongly.
    """
    unit_i = synth.RD.repeats["I"]
    run5 = "G" + "C" * 5 + "A"
    assert run5 in unit_i  # fixture: dictionary unit I ends GCG C5 A
    cons = synth.allele(["X", "I", "X"])
    stutter = cons.replace("GCG" + "C" * 5 + "A", "GCG" + "C" * 6 + "A")
    merged = cons.replace("GCG" + "C" * 5 + "A", "GCC" + "C" * 5 + "A")  # separator G -> C
    reads = [cons] * 10 + [stutter] * 6 + [merged] * 8
    new, changes = homopolymer_vote(cons, reads, min_len=S.hp_vote_min_run)
    assert new == cons and changes == 0


def test_homopolymer_vote_ignores_reads_that_interrupt_a_merged_draft_run() -> None:
    """Task 13b (D2_hifi): the draft lost unit A's separator G (``GCG CCC G`` drafted
    as ``G C5 G``). Reads that carry the G interrupt the drafted run; their longest C
    stretch (3) is not a length of that run, so they must not vote it down to C3 (the
    pileup round, not the vote, restores the G). The full polish recovers the truth.
    """
    unit_a = synth.RD.repeats["A"]
    truth_tail = "GCG" + "C" * 3 + "GCA"
    assert unit_a.endswith(truth_tail)  # fixture: dictionary unit A ends GCG CCC GCA
    truth = synth.allele(["X", "A", "X"])
    draft = synth.allele(["X", unit_a[: -len(truth_tail)] + "G" + "C" * 5 + "GCA", "X"])
    reads = [truth] * 14 + [draft] * 4
    voted, changes = homopolymer_vote(draft, reads, min_len=S.hp_vote_min_run)
    assert voted == draft and changes == 0
    polished, _info = polish(draft, reads, **POLISH)
    assert polished == truth
