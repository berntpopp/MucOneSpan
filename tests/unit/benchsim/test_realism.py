"""Realism metrics measured on the VNTR interval against per-read truth.

Synthetic sources and reads with known edits. edlib (extra ``bench``) is a
Python extension, not an external tool; tests that align reads skip with a
reason where it has no wheel (e.g. Python 3.14).
"""

import gzip
import json
import math
import random
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, RealismConfig
from muc_one_span.benchsim.read_truth import ReadTruth
from muc_one_span.benchsim.realism import aggregate, case_metrics, read_metrics
from muc_one_span.nomenclature import revcomp

SRC = "ACGTACGTAA" + "C" * 7 + "GATTACAGATTACA" * 4
WHOLE = {1: (0, len(SRC))}
REALISM = DEFAULT_BENCH_CONFIG.realism
HEADER = "read_id\thap\tmolecule\tkind\tstrand\tsrc_start\tsrc_end\tn_hp_edits\thp_edits\tdetail\n"


@pytest.fixture
def edlib() -> object:
    return pytest.importorskip("edlib", reason="edlib (extra 'bench') is not installed")


def _fq(path: Path, reads: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f"@{n} extra\n{s}\n+\n{'I' * len(s)}\n" for n, s in reads))
    return path


def _row(rid: str, hap: int, kind: str, strand: str, start: int, end: int) -> ReadTruth:
    return ReadTruth(rid, hap, 1, kind, strand, start, end, "", "")


def _seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mutate(seq: str, pos: int) -> str:
    return seq[:pos] + ("A" if seq[pos] != "A" else "C") + seq[pos + 1 :]


def test_metrics_per_strand(tmp_path: Path, edlib: object) -> None:
    short_c = SRC.replace("C" * 7, "C" * 6)
    rows = [
        _row("r1", 1, "full", "+", 0, len(SRC)),
        _row("r2", 1, "full", "+", 0, len(SRC)),
        _row("r3", 1, "full", "-", 0, len(SRC)),
    ]
    fq = _fq(tmp_path / "r.fq", [("r1", SRC), ("r2", short_c), ("r3", revcomp(SRC))])
    m = read_metrics(fq, rows, {1: SRC}, WHOLE)
    assert m["c7_correct"]["+"] == 0.5 and m["c7_correct"]["-"] == 1.0
    assert m["c7_correct"]["both"] == pytest.approx(2 / 3)
    assert m["error_rate"]["-"] == 0.0 and m["deletion_rate"]["+"] > 0
    assert m["deletion_rate"]["+"] == pytest.approx(1 / (2 * len(SRC)))
    assert m["insertion_rate"]["all"] == 0.0 and m["mismatch_rate"]["all"] == 0.0
    assert m["n_aligned"] == {"+": 2, "-": 1}


def test_mismatch_and_insertion_counted_on_reverse_strand(tmp_path: Path, edlib: object) -> None:
    edited = _mutate(SRC, 30)[:40] + "T" + SRC[40:]
    rows = [_row("r1", 1, "full", "-", 0, len(SRC))]
    m = read_metrics(_fq(tmp_path / "r.fq", [("r1", revcomp(edited))]), rows, {1: SRC}, WHOLE)
    counts = m["error_counts"]["-"]
    assert counts["ref_bases"] == len(SRC)
    assert (counts["mismatch"], counts["ins"], counts["del"]) == (1, 1, 0)
    assert m["error_counts"]["+"] == {"ref_bases": 0, "mismatch": 0, "ins": 0, "del": 0}
    assert m["error_rate"]["-"] > 0 and m["error_rate"]["+"] is None


def test_only_vntr_of_spanning_fragments_is_scored(tmp_path: Path, edlib: object) -> None:
    hap = _seq(400, 1)
    rows = [
        _row("f1", 1, "fragment", "+", 50, 350),  # spans; flank error at 60
        _row("f2", 1, "fragment", "-", 0, 120),  # overlaps, does not span
        _row("f3", 1, "fragment", "+", 320, 400),  # outside the VNTR
    ]
    reads = [
        ("f1", _mutate(hap, 60)[50:350]),
        ("f2", revcomp(_mutate(hap, 110)[0:120])),
        ("f3", hap[320:400]),
    ]
    m = read_metrics(_fq(tmp_path / "r.fq", reads), rows, {1: hap}, {1: (100, 300)})
    assert m["error_counts"]["all"]["ref_bases"] == 200
    assert m["error_rate"] == {"+": 0.0, "-": None, "all": 0.0}
    assert m["n_aligned"] == {"+": 1, "-": 0}
    assert m["spanning_frac"] == 0.5  # 1 spanning of 2 overlapping the VNTR
    assert m["offtarget_frac"] is None and m["log_ratio_slope"] is None


FLANK_L, FLANK_R = _seq(40, 7), _seq(30, 8)


def _amplicon_case(
    tmp_path: Path, n_short: int, n_long: int, config: RealismConfig = REALISM
) -> dict:
    short, long_ = _seq(600, 2), _seq(1200, 3)
    amp = {1: FLANK_L + short + FLANK_R, 2: FLANK_L + long_ + FLANK_R}
    vntr = {1: (40, 640), 2: (40, 1240)}
    rows, reads = [], []
    for i in range(n_short):
        rows.append(_row(f"s{i}", 1, "full", "+" if i % 2 else "-", 0, 670))
        seq = _mutate(amp[1], 5)  # a primer-flank error, outside the VNTR
        reads.append((f"s{i}", seq if i % 2 else revcomp(seq)))
    for i in range(n_long):
        rows.append(_row(f"l{i}", 2, "full", "+", 0, 1270))
        reads.append((f"l{i}", amp[2]))
    rows += [
        _row("sm1", 2, "smear", "+", 0, 1270),
        _row("sm2", 1, "smear", "+", 0, 670),
        _row("ch1", 1, "chimera", "+", 0, 670),
        _row("ot1", 1, "offtarget", "+", 10, 300),
    ]
    reads += [
        ("sm1", amp[2][:340] + amp[2][740:]),  # VNTR 800 bp: between the alleles
        ("sm2", amp[1][:240] + amp[1][540:]),  # VNTR 300 bp: below the short allele
        ("ch1", amp[1]),  # on the short peak (chimera by truth only)
        ("ot1", short[10:300]),
    ]
    return read_metrics(_fq(tmp_path / f"a{n_short}.fq", reads), rows, amp, vntr, config)


def test_amplicon_products_offpeak_and_ratio(tmp_path: Path, edlib: object) -> None:
    m = _amplicon_case(tmp_path, n_short=4, n_long=2)
    n_reads, n_products = 10, 9
    assert m["error_rate"]["all"] == 0.0  # flank errors are out of scope
    assert m["error_counts"]["all"]["ref_bases"] == 4 * 600 + 2 * 1200
    assert m["offtarget_frac"] == pytest.approx(1 / n_reads)
    assert m["smear_frac"] == pytest.approx(2 / n_products)
    assert m["chimera_frac"] == pytest.approx(1 / n_products)
    assert m["span_between_alleles_frac"] == pytest.approx(1 / n_products)
    assert m["span_below_short_frac"] == pytest.approx(1 / n_products)
    assert m["span_off_gt1unit_frac"] == pytest.approx(2 / n_products)
    assert m["allele_ratio"] == {"delta_units": 10.0, "log_ratio": pytest.approx(math.log(2 / 5))}
    assert m["log_ratio_slope"] == pytest.approx(math.log(2 / 5) / 10)
    lt_key, ge_key = REALISM.size_keys
    lt, ge = m["span_offset_hist"][lt_key], m["span_offset_hist"][ge_key]
    assert sum(lt) == n_products and sum(ge) == 0 and lt[-REALISM.offset_bin_min] == 7
    assert len(lt) == REALISM.n_bins


def test_amplicon_size_split_is_configured(tmp_path: Path, edlib: object) -> None:
    small = replace(REALISM, size_split_units=1)
    m = _amplicon_case(tmp_path, n_short=2, n_long=2, config=small)
    assert set(m["span_offset_hist"]) == set(small.size_keys)
    assert sum(m["span_offset_hist"][small.size_keys[0]]) == 0


def test_aggregate_medians_rates_pools_counts_fits_slope(tmp_path: Path, edlib: object) -> None:
    a = _amplicon_case(tmp_path, n_short=4, n_long=2)
    b = _amplicon_case(tmp_path, n_short=2, n_long=2)
    b["error_rate"] = {"+": 0.02, "-": None, "all": 0.01}
    agg = aggregate([a, b])
    assert agg["n_cases"] == 2
    assert agg["error_rate"] == {"+": 0.01, "-": 0.0, "all": 0.005}
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
    empty = aggregate([])
    assert empty["log_ratio_slope"] is None and empty["error_rate"]["all"] is None


def test_missing_edlib_names_extra(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(sys.modules, "edlib", None)
    rows = [_row("r1", 1, "full", "+", 0, len(SRC))]
    with pytest.raises(ImportError, match=r"muc-one-span\[bench\]"):
        read_metrics(_fq(tmp_path / "r.fq", [("r1", SRC)]), rows, {1: SRC}, WHOLE)


def test_unknown_read_and_missing_source_rejected(tmp_path: Path, edlib: object) -> None:
    fq = _fq(tmp_path / "r.fq", [("r1", SRC)])
    with pytest.raises(ValueError, match="not in the FASTQ"):
        read_metrics(fq, [_row("zz", 1, "full", "+", 0, len(SRC))], {1: SRC}, WHOLE)
    with pytest.raises(ValueError, match="no source sequence"):
        read_metrics(fq, [_row("r1", 2, "full", "+", 0, len(SRC))], {1: SRC}, WHOLE)


def _case_dir(tmp_path: Path, geometry: dict | None) -> Path:
    case_dir = tmp_path / "dev" / "d1"
    (case_dir / "truth").mkdir(parents=True)
    (case_dir / "reads").mkdir()
    hap = "TTTT" + SRC + "GGGG"
    (case_dir / "truth" / "d1.001.simulated.fa").write_text(f">haplotype_1\n{hap}\n")
    with gzip.open(case_dir / "truth" / "d1_read_truth.tsv.gz", "wt") as fh:
        fh.write(HEADER + f"r1\t1\t1\tfull\t+\t0\t{len(hap)}\t0\t\t\n")
    _fq(case_dir / "reads" / "d1_amplicon_ont.fastq", [("r1", hap)])
    case = {"design_id": "d1", "profile": "ont_amplicon_r10"}
    if geometry is not None:
        case["geometry"] = geometry
    (case_dir / "case.json").write_text(json.dumps(case))
    return case_dir


def test_case_metrics_reads_recorded_geometry(tmp_path: Path, edlib: object) -> None:
    n = len(SRC)
    geometry = {
        "frame": "amplicon",
        "primers": {"forward": "TTTT", "reverse": "CCCC"},
        "amplicon": {"1": [0, n + 8]},
        "flank_ext": None,
        "vntr_haplotype": {"1": [4, n + 4]},
        "vntr_source": {"1": [4, n + 4]},
    }
    m = case_metrics(_case_dir(tmp_path, geometry))
    assert m["error_counts"]["all"]["ref_bases"] == n and m["c7_correct"]["+"] == 1.0


def test_case_metrics_legacy_case_needs_config(tmp_path: Path, edlib: object) -> None:
    case_dir = _case_dir(tmp_path, None)
    with pytest.raises(ValueError, match="muconeup_config"):
        case_metrics(case_dir)
    config = tmp_path / "config.json"
    primers = {"forward_primer": "TTTT", "reverse_primer": "CCCC"}
    config.write_text(json.dumps({"amplicon_params": primers}))
    truth = SimpleNamespace(haplotypes=(SimpleNamespace(sequence=SRC),))
    with (
        patch("muc_one_span.benchsim.realism.load_truth", return_value=truth),
        patch(
            "muc_one_span.benchsim.realism.load_repeat_dictionary",
            return_value=SimpleNamespace(flanking_left="TTTT"),
        ),
    ):
        m = case_metrics(case_dir, muconeup_config=config)
    assert m["error_counts"]["all"]["ref_bases"] == len(SRC)
