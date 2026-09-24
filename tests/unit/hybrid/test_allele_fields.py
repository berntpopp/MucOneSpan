"""Hybrid selection status and per-allele contract: nothing unresolved is silently resolved."""

from __future__ import annotations

import pytest

from muc_one_span.hybrid.allele_fields import (
    GROUP_SIZE,
    MAX_ALLELES,
    REJECTED_PEAK,
    RESOLVED,
    SINGLE_SITE,
    UNASSIGNED_SPANNING,
    allele_info,
    depth_status,
    selection_detail,
    selection_status,
)
from muc_one_span.hybrid.lengths import GATE_RELEVANT_REJECTIONS, LengthModel
from muc_one_span.hybrid.spans import SpanRead
from muc_one_span.settings import DEFAULT_SETTINGS, HybridSettings
from tests.unit.hybrid import synth

H = DEFAULT_SETTINGS.hybrid
UNIT = synth.RD.repeat_length_bp


def _model(*reasons: str, total: int = 100) -> LengthModel:
    rejected = [{"center_bp": 1.0, "units": 1, "support": 3, "reason": r} for r in reasons]
    return LengthModel([], rejected, [], [], total)


def test_clean_sample_is_resolved() -> None:
    assert selection_status(_model("smear", "noise"), ["none"], 1, 0.0, H) == RESOLVED


@pytest.mark.parametrize("reason", sorted(GATE_RELEVANT_REJECTIONS - {"max_alleles"}))
def test_gate_relevant_rejection_is_unresolved_rejected_peak(reason: str) -> None:
    assert selection_status(_model(reason), ["none"], 1, 0.0, H) == REJECTED_PEAK


def test_third_peak_or_group_is_unresolved_max_alleles() -> None:
    assert selection_status(_model("max_alleles"), ["none", "none"], 2, 0.0, H) == MAX_ALLELES
    assert selection_status(_model(), ["linked_sites", "none"], 3, 0.0, H) == MAX_ALLELES


def test_unconfirmed_splits_have_their_own_status() -> None:
    assert selection_status(_model(), ["unconfirmed_single_site"], 1, 0.0, H) == SINGLE_SITE
    assert selection_status(_model(), ["unconfirmed_group_size"], 1, 0.0, H) == GROUP_SIZE


def test_unassigned_spanning_fraction_uses_the_configured_limit() -> None:
    frac = H.max_unassigned_spanning_fraction
    assert selection_status(_model(), ["none"], 1, frac, H) == RESOLVED
    above = frac + (1 - frac) / 2
    assert selection_status(_model(), ["none"], 1, above, H) == UNASSIGNED_SPANNING
    loose = HybridSettings(max_unassigned_spanning_fraction=1.0)
    assert selection_status(_model(), ["none"], 1, above, loose) == RESOLVED


def test_every_unresolved_status_blocks_a_negative_call() -> None:
    for status in (MAX_ALLELES, SINGLE_SITE, GROUP_SIZE, REJECTED_PEAK, UNASSIGNED_SPANNING):
        assert status.startswith("unresolved")  # clinical_gates.allele_gate_reasons prefix


def test_detail_names_rejected_peaks_and_unassigned_reads() -> None:
    text = selection_detail(_model("smear_ambiguous", "smear"), 1, 7, 0.07)
    assert "smear_ambiguous" in text and "smear," not in text and "7 spanning" in text


def test_depth_status_follows_settings() -> None:
    h = HybridSettings(depth_adequate_spanning=5, depth_low_spanning=2)
    assert [depth_status(n, h) for n in (1, 2, 5)] == ["insufficient", "low", "adequate"]


def test_allele_info_lengths_come_from_unit_and_layout() -> None:
    members = [SpanRead("r", "A" * 10, 20.0, "+", 0, "motif")]
    seq = "A" * (UNIT * 12)
    info = allele_info(
        "allele_1", seq, members, 0, "length", [], (RESOLVED, ""), UNIT, 4, H, concordance=0.75
    )
    assert (info["length"], info["canonical_repeats"]) == (12, 8)
    assert info["depth_status"] == "insufficient" and info["phase_status"] == "phased"
    assert info["allele_genotype_status"] == "not_applicable_read_consensus"
    residual = [{"pos": 1, "ref": "A", "alt": "C", "af": 0.4, "n": 10}]
    info = allele_info(
        "allele_1",
        seq,
        members,
        0,
        "none",
        residual,
        (RESOLVED, ""),
        UNIT,
        4,
        H,
        concordance=1.0,
    )
    assert info["allele_genotype_status"] == "residual_heterogeneity"
    assert info["independent_haplotype_evidence"] is False


def test_allele_info_reports_real_hybrid_evidence_alongside_the_dictionary_confidence() -> None:
    """The ladder's classify.py ``confidence``/``allele_confidence`` is a dictionary-fit
    heuristic fed identically regardless of engine (see pipeline_tail.py), so it is not
    meaningful hybrid evidence. ``allele_info`` additively reports the hybrid engine's own
    read-support evidence (``consensus_concordance_fraction``) plus an explicit marker
    that the ladder confidence does not apply to hybrid reconstruction quality.
    """
    members = [SpanRead("r", "A" * 10, 20.0, "+", 0, "motif")]
    seq = "A" * (UNIT * 3)
    info = allele_info(
        "allele_1", seq, members, 0, "length", [], (RESOLVED, ""), UNIT, 0, H, concordance=0.6
    )
    assert info["consensus_concordance_fraction"] == 0.6
    assert info["classification_confidence_status"] == "not_applicable_dictionary_fit_heuristic"
