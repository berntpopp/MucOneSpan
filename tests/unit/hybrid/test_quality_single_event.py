"""A low-accuracy artefact site must not hide an equal-length carrier's event.

Reproduces the shape of the v4 dev HiFi pathogenic cases reported INCONCLUSIVE:
both alleles have the same length (one length peak) and differ at one indel event
(dupA: the carrier X unit reads C7-AA where the wild-type C unit reads C6-AA). A
subset of reads with lower base quality shares a systematic error at another,
unlinked site whose minor allele fraction reaches ``het_af_min``. The peak then has
two candidate events, so the single-event split  is refused and the carrier
event is scored on a merged consensus (read support discordant, INCONCLUSIVE).

With ``phase_quality_single_event`` the 15j low-accuracy test (rank-sum test on mean
base quality and the upper confidence bound of the high-quality minor allele
fraction) is applied to the extra sites of such a peak. When the sites it keeps form
exactly one event, the peak is split on that event under every single-event gate
(share bound, group size, differing drafts). The rule may only turn INCONCLUSIVE into
PATHOGENIC: the selection status stays unresolved, so a negative call stays blocked
whatever the split finds.

The rule is off by default: a wild-type peak with a site-specific
+1 C excess of 0.35-0.40 at one C7 run (beyond the single-event split's validated
range, where it cannot tell an artefact from a real minor) plus a low-accuracy artefact
site was INCONCLUSIVE and becomes PATHOGENIC with the rule on. The stress test below
covers the validated range (excess up to 0.30).
"""

from __future__ import annotations

import dataclasses
import random
from pathlib import Path
from typing import Any

import pytest

from muc_one_span.hybrid.align import rc
from muc_one_span.hybrid.phase import split_by_linked_sites
from muc_one_span.hybrid.spans import PHRED_OFFSET, ReadRecord, categorize_reads
from muc_one_span.settings import DEFAULT_SETTINGS, HybridSettings, RuntimeSettings
from tests.unit.hybrid import synth
from tests.unit.hybrid import test_quality_sites as quality
from tests.unit.hybrid import test_single_event as base
from tests.unit.hybrid import test_single_event_gate as gate

# The heavy synthetic sweeps below carry the safety_sweep marker (their own CI job and
# make test-unit, never skipped); the cheap tests run in the core suite.

S = quality.S  # the 15j test is opt-in since switched on explicitly
NEGATIVE = "NO_PATHOGENIC_VARIANT_DETECTED"
N_PER_ALLELE = base.N_PER_ALLELE
# Share of all reads in the low-accuracy subset carrying the artefact column (the
# 15j panel shape: an allele fraction just above het_af_min).
SUBSET_FRAC = quality.SUBSET_FRAC
SEEDS = (0, 1)
ARTEFACT_UNIT = quality.ARTEFACT_UNIT
# Per-read mean Phred ranges by read quality (ONT-like: disjoint; HiFi-like: overlapping).
PROFILES = {
    "ont": {"good": (18, 26), "poor": (8, 17)},
    "hifi": {"good": quality.GOOD_Q, "poor": quality.POOR_Q},
}
# Mosaic dupC minority (a third haplotype of the wild-type X allele) as a share of all
# reads; 0 is the wild-type counterpart. Low-quality non-carrier reads read the dupC run
# as C8 by insertion stutter at POOR_STUTTER_AT_EVENT (site-specific, the adversarial
# case of test_quality_sites_safety), so the minor set is enriched in poor reads.
MOSAIC_AFS = (0.0, 0.15, 0.22, 0.30)
MOSAIC_SEEDS = (0, 1, 2)
MOSAIC_UNIT = 12  # 0-based inner unit whose C7 run carries the mosaic dupC
POOR_FRAC = 0.4
POOR_STUTTER_AT_EVENT = 0.3


def _settings(**changes: object) -> RuntimeSettings:
    hybrid: HybridSettings = dataclasses.replace(S, **changes)  # type: ignore[arg-type]
    return dataclasses.replace(DEFAULT_SETTINGS, hybrid=hybrid)


# The rule is off by default (see settings_hybrid); these tests
# switch it on.
ON = _settings(phase_quality_single_event=True)
# The split mechanism tests run at the pre-15l single-event share floor (het_af_min):
# at this depth (N_PER_ALLELE per allele) the floor phase_single_event_min_share
# refuses the split. The never-PATHOGENIC tests keep the default floor (ON).
ON_SPLIT = _settings(phase_quality_single_event=True, phase_single_event_min_share=S.het_af_min)


def _read(template: str, strand: str, seed: int, name: str, q: tuple[int, int]) -> ReadRecord:
    rng = random.Random(seed)
    stuttered = base._stutter(template, rng, *base.ONT_LIKE_STUTTER[strand])
    read = synth.reads(stuttered, 1, err=base.ERR, seed=rng.randrange(1 << 30), strand_mix=False)
    seq = read[0].seq if strand == "+" else rc(read[0].seq)
    return ReadRecord(name, seq, chr(rng.randint(*q) + PHRED_OFFSET) * len(seq))


def _with_artefact(template: str) -> str:
    """``template`` with the artefact base at one isolated column of ARTEFACT_UNIT.

    The column is chosen on the wild-type allele and moved by the template's length
    difference (the event sits before ARTEFACT_UNIT), so every allele carries the
    artefact at the same consensus column.
    """
    wild_type = synth.allele(base.WT)
    pos, alt = quality._column(wild_type, ARTEFACT_UNIT)
    pos += len(template) - len(wild_type)
    return template[:pos] + alt + template[pos + 1 :]


def _het(
    seed: int,
    carrier: list[str] = base.MUT,
    profile: str = "hifi",
    per_allele: int = N_PER_ALLELE,
) -> list[ReadRecord]:
    """Equal-length heterozygote (``carrier`` / wild-type C unit) with a low-accuracy
    subset, drawn from both alleles, sharing one artefact column."""
    rng = random.Random(seed)
    q = PROFILES[profile]
    out = []
    for name, inner in (("m", carrier), ("w", base.WT)):
        template = synth.allele(inner)
        for i in range(per_allele):
            strand = "+" if rng.random() < 1 / 2 else "-"
            poor = rng.random() < SUBSET_FRAC
            source = _with_artefact(template) if poor else template
            out.append(
                _read(
                    source,
                    strand,
                    rng.randrange(1 << 30),
                    f"{name}{i}",
                    q["poor" if poor else "good"],
                )
            )
    return out


def _peak(records: list[ReadRecord]) -> Any:
    return categorize_reads(records, base.ANCH, S).spanning


# --- preconditions: the shape of the dev cases -----------------------------------------


@pytest.mark.safety_sweep
def test_rule_is_off_by_default(tmp_path: Path) -> None:
    summary, decision = base._run(tmp_path, _het(SEEDS[0]))
    assert summary["hybrid"]["split_bases"] == ["unconfirmed_single_site"]
    assert decision["state"] == "INCONCLUSIVE"


@pytest.mark.safety_sweep
@pytest.mark.parametrize("seed", SEEDS)
def test_artefact_adds_a_second_event_and_blocks_the_split(tmp_path: Path, seed: int) -> None:
    records = _het(seed)
    res = split_by_linked_sites(synth.allele(base.WT), _peak(records), S, random.Random(S.seed))
    assert res.basis == "unconfirmed_single_site" and len(res.sites) >= 2, res.sites
    summary, decision = base._run(tmp_path, records, _settings(phase_quality_single_event=False))
    assert summary["hybrid"]["split_bases"] == ["unconfirmed_single_site"]
    assert decision["state"] == "INCONCLUSIVE"


# --- the rule ---------------------------------------------------------------------------


@pytest.mark.safety_sweep
@pytest.mark.parametrize("seed", SEEDS)
def test_carrier_event_is_split_and_called_after_the_artefact_is_dropped(
    tmp_path: Path, seed: int
) -> None:
    summary, decision = base._run(tmp_path, _het(seed), ON_SPLIT)
    block = summary["hybrid"]
    assert block["split_bases"] == ["single_event"], block
    assert block["quality_associated_sites"], "the dropped artefact must be recorded"
    assert all(s["af_high_quality"] < S.het_af_min for s in block["quality_associated_sites"])
    mutations = [m for c in summary["classifications"].values() for m in c["mutations"]]
    dupa = [m for m in mutations if m.get("mutation_name") == "dupA"]
    assert len(dupa) == 1 and dupa[0]["read_support"]["status"] == "supported", dupa
    assert decision["state"] == "PATHOGENIC", decision["details"]


@pytest.mark.safety_sweep
@pytest.mark.parametrize("seed", SEEDS)
def test_a_quality_enabled_split_never_resolves_the_selection(tmp_path: Path, seed: int) -> None:
    """Wild-type equal-length heterozygote (X vs C: a run-length difference, no event):
    the split is made, but the sample stays INCONCLUSIVE, never NEGATIVE."""
    summary, decision = base._run(tmp_path, _het(seed, carrier=["X"] * len(base.WT)), ON_SPLIT)
    block = summary["hybrid"]
    assert block["split_bases"] == ["single_event"], block
    assert block["selection_status"] == "unresolved_single_site"
    assert "low-accuracy" in block["selection_detail"]
    assert decision["state"] == "INCONCLUSIVE", decision["details"]


@pytest.mark.safety_sweep
@pytest.mark.parametrize("seed", SEEDS[:1])
def test_a_quality_enabled_split_never_resolves_the_selection_at_the_defaults(
    tmp_path: Path, seed: int
) -> None:
    """The same shape at the default share floor, with the bound's full sample:
    the split is made on the single event and keeps the selection unresolved with the
    low-accuracy reason, so the sample is INCONCLUSIVE, never NEGATIVE."""
    records = _het(
        seed, carrier=["X"] * len(base.WT), per_allele=S.phase_single_event_bound_reads // 2
    )
    summary, decision = base._run(tmp_path, records, ON)
    block = summary["hybrid"]
    assert block["split_bases"] == ["single_event"], block
    assert block["selection_status"] == "unresolved_single_site", block
    assert "low-accuracy" in block["selection_detail"]
    assert decision["state"] == "INCONCLUSIVE", decision["details"]


def test_rule_needs_the_quality_test() -> None:
    """With the 15j test off (phase_quality_alpha 0) the rule would be a silent no-op:
    the configuration is refused, naming both keys."""
    with pytest.raises(ValueError, match=r"phase_quality_single_event.*phase_quality_alpha"):
        _settings(phase_quality_single_event=True, phase_quality_alpha=0.0)


@pytest.mark.safety_sweep
def test_artefact_on_good_reads_keeps_the_peak_unsplit(tmp_path: Path) -> None:
    """A second site carried by high-quality reads is not explained: no split."""
    rng = random.Random(SEEDS[0])
    good = PROFILES["hifi"]["good"]
    records = [
        _read(
            _with_artefact(synth.allele(inner))
            if rng.random() < SUBSET_FRAC
            else synth.allele(inner),
            "+" if rng.random() < 1 / 2 else "-",
            rng.randrange(1 << 30),
            f"{name}{i}",
            good,
        )
        for name, inner in (("m", base.MUT), ("w", base.WT))
        for i in range(N_PER_ALLELE)
    ]
    summary, decision = base._run(tmp_path, records, ON)
    assert summary["hybrid"]["split_bases"] == ["unconfirmed_single_site"]
    assert not summary["hybrid"]["quality_associated_sites"]
    assert decision["state"] == "INCONCLUSIVE"


@pytest.mark.safety_sweep
def test_two_peak_samples_are_unaffected(tmp_path: Path) -> None:
    """The setting acts on single-peak (equal-length) genotypes only; the 15j two-peak
    rule is unchanged."""
    records = quality._sample(SEEDS[0])
    for sub in ("on", "off"):
        (tmp_path / sub).mkdir()
    on, _ = base._run(tmp_path / "on", records, ON)
    off, _ = base._run(tmp_path / "off", records, _settings(phase_quality_single_event=False))
    assert on["hybrid"]["split_bases"] == off["hybrid"]["split_bases"] == ["none", "none"]
    assert on["hybrid"]["quality_associated_sites"] == off["hybrid"]["quality_associated_sites"]


def test_setting_is_validated() -> None:
    with pytest.raises(ValueError, match="phase_quality_single_event"):
        HybridSettings(phase_quality_single_event=1)  # type: ignore[arg-type]


# --- adversarial: a mosaic dupC minority with quality-correlated stutter ---------------


def _mosaic(profile: str, af: float, seed: int) -> list[ReadRecord]:
    """Wild-type equal-length heterozygote (X / C unit) plus a low-accuracy artefact
    subset; ``af`` of all reads are a mosaic X haplotype carrying dupC at MOSAIC_UNIT."""
    rng = random.Random(seed)
    q = PROFILES[profile]
    wt_x = synth.allele(["X"] * len(base.WT))
    inner = ["X"] * len(base.WT)
    inner[MOSAIC_UNIT] = synth.dupc()
    dupc = synth.allele(inner)
    n_total = 2 * N_PER_ALLELE
    n_mosaic = round(af * n_total)
    out = []
    for i in range(n_total):
        strand = "+" if rng.random() < 1 / 2 else "-"
        poor = rng.random() < POOR_FRAC
        if i < n_mosaic:
            template = dupc
        elif i % 2:
            template = synth.allele(base.WT)
        else:
            stutter = poor and rng.random() < POOR_STUTTER_AT_EVENT
            template = dupc if stutter else wt_x
        if poor and rng.random() < SUBSET_FRAC / POOR_FRAC:
            template = _with_artefact(template)
        out.append(
            _read(template, strand, rng.randrange(1 << 30), f"r{i}", q["poor" if poor else "good"])
        )
    return out


@pytest.mark.safety_sweep
@pytest.mark.parametrize("profile", sorted(PROFILES))
@pytest.mark.parametrize("af", MOSAIC_AFS)
@pytest.mark.parametrize("seed", MOSAIC_SEEDS)
def test_mosaic_dupc_is_never_negative_and_wild_type_never_pathogenic(
    tmp_path: Path, profile: str, af: float, seed: int
) -> None:
    summary, decision = base._run(tmp_path, _mosaic(profile, af, seed), ON)
    detail = (summary["hybrid"]["split_bases"], summary["hybrid"]["quality_associated_sites"])
    if af:
        assert decision["state"] != NEGATIVE, detail
    else:
        assert decision["state"] != "PATHOGENIC", detail


# --- adversarial: the single-event FP stress shape plus a low-accuracy artefact -------


def _stress(excess: float, stutter: str, seed: int) -> list[ReadRecord]:
    """test_single_event_gate's wild-type stress reads (a site-specific +1 C excess at
    one C7 run, both strands) with a low-accuracy subset sharing an artefact column:
    without the rule the artefact is a second event that blocks any split."""
    cons = synth.allele(["X"] * len(base.WT))
    start = base._unit_run(cons, gate.STRESS_UNIT)
    edited = cons[:start] + "C" + cons[start:]
    rng = random.Random(seed)
    q = PROFILES["hifi"]
    out = []
    for i in range(gate.STRESS_DEPTH):
        strand = "+" if rng.random() < 1 / 2 else "-"
        source = edited if rng.random() < excess else cons
        poor = rng.random() < SUBSET_FRAC
        template = _with_artefact(source) if poor else source
        stuttered = base._stutter(template, rng, *gate.STUTTER[stutter][strand])
        read = synth.reads(
            stuttered, 1, err=base.ERR, seed=rng.randrange(1 << 30), strand_mix=False
        )
        seq = read[0].seq if strand == "+" else rc(read[0].seq)
        qual = chr(rng.randint(*q["poor" if poor else "good"]) + PHRED_OFFSET) * len(seq)
        out.append(ReadRecord(f"x{seed}_{i}", seq, qual))
    return out


@pytest.mark.safety_sweep
@pytest.mark.parametrize("stutter", sorted(gate.STUTTER))
@pytest.mark.parametrize("excess", gate.EXCESS)
def test_wild_type_stress_with_an_artefact_is_never_pathogenic(
    tmp_path: Path, excess: float, stutter: str
) -> None:
    summary, decision = base._run(tmp_path, _stress(excess, stutter, S.seed), ON)
    assert decision["state"] != "PATHOGENIC", (excess, summary["hybrid"]["split_bases"])
    assert decision["state"] != NEGATIVE or not summary["hybrid"]["quality_associated_sites"]
