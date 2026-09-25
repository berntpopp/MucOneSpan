"""Task 15f fix round 1: an extrapolated stutter profile must keep the alleles apart.

Real ONT "+" reads (PRJEB92208 MP1) show 9.8% deletion stutter at C6 runs and 25.7%
at C7 runs, but only 20.6% at the true C8 run: the trend saturates. Extrapolating it
predicted a C8 allele read as C7 more often than as C8, so the C7 and C8 alleles were
nearly indistinguishable and a 30% wild-type mixture looked like a pure dupC. A strand
whose extrapolated event profile puts more than ``hp_stutter_max_event_confusion`` of
its mass on the no-event run length falls back to the shift model.
"""

from __future__ import annotations

import dataclasses
from typing import Any

import pytest

from muc_one_span.report import compute_clinical_decision
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import test_evidence_stutter as base
from tests.unit.test_clinical_decision import BASE, _gated_summary

pytest.importorskip("edlib", reason="edlib (extra 'hybrid') is not installed")
S = HybridSettings()
SHIFT = dataclasses.replace(S, hp_stutter_model="shift")
# Deletion stutter per run length on real ONT "+" reads (MP1 probe): C6, C7 and the true
# C8 run; shorter runs continue the C6 -> C7 trend downwards.
ONT_PLUS = {6: 0.098, 7: 0.257, 8: 0.206}


def _ont_plus(length: int) -> float:
    if length in ONT_PLUS:
        return ONT_PLUS[length]
    if length > max(ONT_PLUS):
        return ONT_PLUS[max(ONT_PLUS)]
    low = min(ONT_PLUS)
    return ONT_PLUS[low] * (ONT_PLUS[low] / ONT_PLUS[low + 1]) ** (low - length)


def _mix(wild_share: float, seed: int, p_del: Any) -> list[tuple[str, str]]:
    wild = round(base.N_READS * wild_share)
    return base._reads([(base.DUPC, base.N_READS - wild), (base.WILD, wild)], seed, p_del)


def _decision(support: dict[str, Any]) -> str:
    mutation = dict(
        BASE,
        vcf_support=False,
        vcf_support_status="not_applicable_read_consensus",
        read_support=support,
    )
    return str(compute_clinical_decision(_gated_summary([mutation]))["state"])


def _support(reads: list[tuple[str, str]], s: HybridSettings = S) -> dict[str, Any]:
    _result, dupc = base._dupc_support(base.DUPC, reads, s)
    return dict(dupc[0]["read_support"])


def test_decision_fixture_calls_a_supported_event_pathogenic() -> None:
    support = _support(base._reads([(base.DUPC, base.N_READS)], base.SEEDS[0]))
    assert support["status"] == "supported", support
    assert _decision(support) == "PATHOGENIC"


@pytest.mark.parametrize("seed", base.SEEDS)
@pytest.mark.parametrize("wild_share", [0.3, 0.4])
def test_ont_plus_shape_wild_type_mixture_is_never_pathogenic(wild_share: float, seed: int) -> None:
    support = _support(_mix(wild_share, seed, _ont_plus))
    assert support["status"] != "supported", support
    assert _decision(support) != "PATHOGENIC", support


@pytest.mark.parametrize("seed", base.SEEDS)
def test_ont_plus_shape_falls_back_to_the_shift_model(seed: int) -> None:
    reads = _mix(0.3, seed, _ont_plus)
    assert _support(reads) == _support(reads, SHIFT)


def test_confusion_guard_is_its_own_setting() -> None:
    """Without the guard the extrapolated ONT "+" profile hides a 30% wild-type share."""
    reads = _mix(0.3, base.SEEDS[0], _ont_plus)
    unguarded = dataclasses.replace(S, hp_stutter_max_event_confusion=1.0)
    assert _support(reads, unguarded)["status"] == "supported"


def test_confusion_guard_is_validated() -> None:
    with pytest.raises(ValueError, match="hp_stutter_max_event_confusion"):
        dataclasses.replace(S, hp_stutter_max_event_confusion=1.5)
