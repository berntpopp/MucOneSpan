"""Realism metrics (sim reads vs truth) and comparison with PRJEB92208 targets.

Synthetic sources and reads with known edits; edlib is the only external
dependency (a Python extension, not an external tool). Target shapes follow
the real `targets/prjeb92208_v1.json` keys, not the brief's illustrative shape.
"""

import math
import random
import sys
from pathlib import Path

import pytest

from muc_one_span.benchsim.read_truth import ReadTruth
from muc_one_span.benchsim.realism import (
    aggregate,
    amplicon_sources,
    haplotype_sources,
    read_metrics,
)
from muc_one_span.benchsim.realism_targets import (
    PROFILE_SECTIONS,
    compare,
    js_distance,
    load_targets,
)
from muc_one_span.nomenclature import revcomp

SRC = "ACGTACGTAA" + "C" * 7 + "GATTACAGATTACA" * 4


def _fq(path: Path, reads: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f"@{n} extra\n{s}\n+\n{'I' * len(s)}\n" for n, s in reads))
    return path


def _row(rid: str, hap: int, kind: str, strand: str, start: int, end: int) -> ReadTruth:
    return ReadTruth(rid, hap, 1, kind, strand, start, end, "", "")


def _seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_metrics_per_strand(tmp_path: Path) -> None:
    short_c = SRC.replace("C" * 7, "C" * 6)
    rows = [
        _row("r1", 1, "full", "+", 0, len(SRC)),
        _row("r2", 1, "full", "+", 0, len(SRC)),
        _row("r3", 1, "full", "-", 0, len(SRC)),
    ]
    fq = _fq(tmp_path / "r.fq", [("r1", SRC), ("r2", short_c), ("r3", revcomp(SRC))])
    m = read_metrics(fq, rows, {1: SRC})
    assert m["c7_correct"]["+"] == 0.5 and m["c7_correct"]["-"] == 1.0
    assert m["c7_correct"]["both"] == pytest.approx(2 / 3)
    assert m["error_rate"]["-"] == 0.0 and m["deletion_rate"]["+"] > 0
    assert m["deletion_rate"]["+"] == pytest.approx(1 / (2 * len(SRC)))
    assert m["insertion_rate"]["all"] == 0.0 and m["mismatch_rate"]["all"] == 0.0
    assert m["n_aligned"] == {"+": 2, "-": 1}


def test_mismatch_and_insertion_counted_on_reverse_strand(tmp_path: Path) -> None:
    edited = SRC[:30] + ("A" if SRC[30] != "A" else "C") + SRC[31:40] + "T" + SRC[40:]
    rows = [_row("r1", 1, "full", "-", 0, len(SRC))]
    m = read_metrics(_fq(tmp_path / "r.fq", [("r1", revcomp(edited))]), rows, {1: SRC})
    counts = m["error_counts"]["-"]
    assert counts["ref_bases"] == len(SRC)
    assert (counts["mismatch"], counts["ins"], counts["del"]) == (1, 1, 0)
    assert m["error_counts"]["+"] == {"ref_bases": 0, "mismatch": 0, "ins": 0, "del": 0}
    assert m["error_rate"]["-"] > 0 and m["error_rate"]["+"] is None


def test_fragment_uses_source_interval_and_spanning_fraction(tmp_path: Path) -> None:
    hap = _seq(400, 1)
    rows = [
        _row("f1", 1, "fragment", "+", 50, 350),
        _row("f2", 1, "fragment", "-", 0, 120),
    ]
    fq = _fq(tmp_path / "r.fq", [("f1", hap[50:350]), ("f2", revcomp(hap[0:120]))])
    m = read_metrics(fq, rows, {1: hap}, vntr_span={1: (100, 300)})
    assert m["error_rate"] == {"+": 0.0, "-": 0.0, "all": 0.0}
    assert m["spanning_frac"] == 0.5
    assert m["offtarget_frac"] is None and m["log_ratio_slope"] is None


def _amplicon_case(tmp_path: Path, n_short: int, n_long: int) -> dict:
    short, long_ = _seq(600, 2), _seq(1200, 3)
    rows, reads = [], []
    for i in range(n_short):
        rows.append(_row(f"s{i}", 1, "full", "+" if i % 2 else "-", 0, 600))
        reads.append((f"s{i}", short if i % 2 else revcomp(short)))
    for i in range(n_long):
        rows.append(_row(f"l{i}", 2, "full", "+", 0, 1200))
        reads.append((f"l{i}", long_))
    rows += [
        _row("sm1", 2, "smear", "+", 0, 1200),
        _row("sm2", 1, "smear", "+", 0, 600),
        _row("ch1", 1, "chimera", "+", 0, 600),
        _row("ot1", 1, "offtarget", "+", 10, 300),
    ]
    reads += [
        ("sm1", long_[:300] + long_[700:]),  # 900 bp: between the alleles
        ("sm2", short[:200] + short[500:]),  # 300 bp: below the short allele
        ("ch1", short),  # on the short peak (chimera by truth only)
        ("ot1", short[10:300]),
    ]
    return read_metrics(_fq(tmp_path / f"a{n_short}.fq", reads), rows, {1: short, 2: long_})


def test_amplicon_products_offpeak_and_ratio(tmp_path: Path) -> None:
    m = _amplicon_case(tmp_path, n_short=4, n_long=2)
    n_reads, n_products = 10, 9
    assert m["offtarget_frac"] == pytest.approx(1 / n_reads)
    assert m["smear_frac"] == pytest.approx(2 / n_products)
    assert m["chimera_frac"] == pytest.approx(1 / n_products)
    assert m["span_between_alleles_frac"] == pytest.approx(1 / n_products)
    assert m["span_below_short_frac"] == pytest.approx(1 / n_products)
    assert m["span_off_gt1unit_frac"] == pytest.approx(2 / n_products)
    assert m["allele_ratio"] == {"delta_units": 10.0, "log_ratio": pytest.approx(math.log(2 / 5))}
    assert m["log_ratio_slope"] == pytest.approx(math.log(2 / 5) / 10)
    lt, ge = m["span_offset_hist"]["lt55u"], m["span_offset_hist"]["ge55u"]
    assert sum(lt) == n_products and sum(ge) == 0 and lt[12] == 7


def test_aggregate_pools_counts_and_fits_slope(tmp_path: Path) -> None:
    a = _amplicon_case(tmp_path, n_short=4, n_long=2)
    b = _amplicon_case(tmp_path, n_short=2, n_long=2)
    agg = aggregate([a, b])
    assert agg["n_cases"] == 2
    assert agg["error_counts"]["all"]["ref_bases"] == (
        a["error_counts"]["all"]["ref_bases"] + b["error_counts"]["all"]["ref_bases"]
    )
    frac = agg["offtarget_frac"]
    assert frac["min"] == pytest.approx(1 / 10) and frac["max"] == pytest.approx(1 / 8)
    xs = [10.0, 10.0]
    ys = [math.log(2 / 5), math.log(2 / 3)]
    expected = sum(x * y for x, y in zip(xs, ys, strict=True)) / sum(x * x for x in xs)
    assert agg["log_ratio_slope"] == pytest.approx(expected)
    assert sum(agg["span_offset_hist"]["lt55u"]) == 9 + 7
    assert aggregate([])["log_ratio_slope"] is None


def test_missing_edlib_names_extra(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "edlib", None)
    rows = [_row("r1", 1, "full", "+", 0, len(SRC))]
    with pytest.raises(ImportError, match=r"muc-one-span\[bench\]"):
        read_metrics(_fq(tmp_path / "r.fq", [("r1", SRC)]), rows, {1: SRC})


def test_unknown_read_and_missing_source_rejected(tmp_path: Path) -> None:
    fq = _fq(tmp_path / "r.fq", [("r1", SRC)])
    with pytest.raises(ValueError, match="not in the FASTQ"):
        read_metrics(fq, [_row("zz", 1, "full", "+", 0, len(SRC))], {1: SRC})
    with pytest.raises(ValueError, match="no source sequence"):
        read_metrics(fq, [_row("r1", 2, "full", "+", 0, len(SRC))], {1: SRC})


def test_packaged_targets_are_public_only() -> None:
    targets = load_targets()
    assert set(targets) == {"_meta", "ont_amplicon_PRJEB92208", "ont_wgs_PRJEB92208"}
    assert set(PROFILE_SECTIONS.values()) <= set(targets)
    amp = targets["ont_amplicon_PRJEB92208"]
    assert amp["hp_P_obs_given_true"]["C7|+"]["p_correct"] == pytest.approx(0.5174)


def test_js_distance_bounds() -> None:
    assert js_distance([1, 2, 3], [2, 4, 6]) == pytest.approx(0.0)
    assert js_distance([1, 0], [0, 1]) == pytest.approx(1.0)
    assert js_distance([0, 0], [1, 1]) is None


TARGETS = {
    "ont_amplicon_PRJEB92208": {
        "hp_P_obs_given_true": {"C7|+": {"p_correct": 0.52}, "C7|-": {"p_correct": 0.89}},
        "error_rates_per_ref_base": {
            "all:total": {"min": 0.016, "median": 0.021, "max": 0.033},
            "+:del": {"min": 0.009, "median": 0.011, "max": 0.017},
        },
        "category_frac": {"off_target": {"min": 0.06, "median": 0.30, "max": 0.77}},
        "span_between_alleles_frac": {"min": 0.01, "median": 0.023, "max": 0.042},
        "span_offset_pmf_15bp_bins": {"lt55u": {"p": [0.0] * 12 + [1.0] + [0.0] * 4}},
        "allelic_ratio": {"through_origin_b_per_unit": -0.056},
    }
}


def test_compare_flags_out_of_tolerance() -> None:
    metrics = {
        "c7_correct": {"+": 0.60, "-": 0.88, "both": 0.7},
        "error_rate": {"all": 0.022, "+": None},
        "deletion_rate": {"+": 0.02},
        "offtarget_frac": 0.9,
        "span_between_alleles_frac": 0.02,
        "span_offset_hist": {"lt55u": [0] * 12 + [9] + [0] * 4, "ge55u": [0] * 17},
        "log_ratio_slope": -0.05,
    }
    res = compare(metrics, TARGETS, "ont_amplicon_r10")
    assert res["c7_correct_+"]["pass"] is False and res["c7_correct_-"]["pass"] is True
    assert res["error_rate_all"]["pass"] is True
    assert res["deletion_rate_+"]["pass"] is False
    assert res["offtarget_frac"]["pass"] is False
    assert res["span_between_alleles_frac"]["pass"] is True
    assert res["span_offset_jsd_lt55u"]["pass"] is True
    assert "span_offset_jsd_ge55u" not in res and "error_rate_+" not in res
    assert res["allele_ratio_slope"]["pass"] is True


def test_compare_aggregate_range_and_median() -> None:
    inside = {"span_between_alleles_frac": {"min": 0.015, "median": 0.024, "max": 0.04}}
    shifted = {"span_between_alleles_frac": {"min": 0.015, "median": 0.035, "max": 0.04}}
    wide = {"span_between_alleles_frac": {"min": 0.0, "median": 0.024, "max": 0.04}}
    key = "span_between_alleles_frac"
    assert compare(inside, TARGETS, "ont_amplicon_PRJEB92208")[key]["pass"] is True
    assert compare(shifted, TARGETS, "ont_amplicon_PRJEB92208")[key]["pass"] is False
    assert compare(wide, TARGETS, "ont_amplicon_PRJEB92208")[key]["pass"] is False


def test_compare_unknown_profile() -> None:
    with pytest.raises(KeyError, match="hifi_amplicon"):
        compare({}, TARGETS, "hifi_amplicon")


def test_load_local_targets(tmp_path: Path) -> None:
    path = tmp_path / "local.json"
    path.write_text('{"custom": {"category_frac": {"spanning": {"min": 0.1, "max": 0.2}}}}')
    res = compare({"spanning_frac": 0.15}, load_targets(path), "custom")
    assert res["spanning_frac"]["pass"] is True


def test_haplotype_and_amplicon_sources(tmp_path: Path) -> None:
    fwd, rev = "GGAGAAAAGG", "GCCGTTGTGC"
    body1, body2 = _seq(50, 5), _seq(80, 6)
    hap1 = "tt" + fwd + body1 + revcomp(rev) + "aa"
    hap2 = "c" + fwd + body2 + revcomp(rev)
    fa = tmp_path / "t.simulated.fa"
    fa.write_text(f">haplotype_2 x\n{hap2}\n>haplotype_1\n{hap1}\n")
    haps = haplotype_sources(fa)
    assert haps == {1: hap1.upper(), 2: hap2.upper()}
    amp = amplicon_sources(haps, fwd, rev)
    assert amp[1] == fwd + body1 + revcomp(rev) and len(amp[2]) == 100
    with pytest.raises(ValueError, match="primer sites not unique"):
        amplicon_sources({1: fwd + fwd + revcomp(rev)}, fwd, rev)
