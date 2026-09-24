"""Task 15d: fast stage calibration of the hybrid length model.

Builds tiny real cases (real bundled repeat dictionary, synthetic noisy spanning
reads via ``tests.unit.hybrid.synth``, a ``load_truth``-loadable truth dir) and
runs the real `hybrid.spans`/`hybrid.lengths`/`hybrid.smear` entry points through
`run_lengths_stage` and `evaluate_lengths_stage` -- no mocking of the length model
itself, only of nothing (the algorithm runs for real, on tiny inputs).
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG
from muc_one_span.benchsim.calibration_grid import base_settings, overlay
from muc_one_span.benchsim.calibration_lengths import (
    LENGTHS_COUNTS,
    LENGTHS_METRICS,
    LENGTHS_RATES,
    _match_alleles,
    evaluate_lengths_stage,
    normalize_lengths_rows,
    run_lengths_stage,
)
from muc_one_span.benchsim.calibration_objective import RATE, load_objective, point_metrics
from muc_one_span.benchsim.generate import FASTQ
from muc_one_span.config import load_repeat_dictionary
from muc_one_span.hybrid.spans import ReadRecord
from muc_one_span.settings import HybridSettings
from tests.unit.hybrid import synth

RD = load_repeat_dictionary()
UNIT = RD.repeat_length_bp
PROFILE = "ont_amplicon_r10"


def _write_config(tmp_path: Path, values: dict[str, object] | None = None) -> Path:
    base = base_settings(None)
    cfg = overlay(base, values or {}, "hybrid")
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "config.json"
    path.write_text(json.dumps(cfg))
    return path


def _write_fastq(path: Path, records: list[ReadRecord]) -> None:
    path.write_text("".join(f"@{r.name}\n{r.seq}\n+\n{r.qual}\n" for r in records))


def _truth(case_dir: Path, hap_units: list[list[str]]) -> None:
    """A ``load_truth``-loadable ``truth/`` dir for haplotypes of given inner units."""
    truth_dir = case_dir / "truth"
    truth_dir.mkdir(parents=True)
    fasta, structures, stats = [], [], []
    for i, units in enumerate(hap_units, start=1):
        name = f"haplotype_{i}"
        structure = synth.PRE + units + synth.POST
        seq = synth.allele(units)
        fasta.append(f">{name}\n{RD.flanking_left + seq + RD.flanking_right}\n")
        structures.append(f"{name}\t{'-'.join(structure)}\n")
        stats.append(
            {
                "repeat_count": len(structure),
                "vntr_length": len(seq),
                "repeat_lengths": [len(RD.repeats[u]) for u in structure],
                "mutant_repeat_count": 0,
                "mutation_details": [],
            }
        )
    (truth_dir / "case.simulated.fa").write_text("".join(fasta))
    (truth_dir / "case.vntr_structure.txt").write_text("".join(structures))
    (truth_dir / "case.simulation_stats.json").write_text(
        json.dumps({"haplotype_statistics": stats, "mutation_info": {}, "provenance": {"seed": 1}})
    )


def _case(
    split_dir: Path,
    design_id: str,
    hap_units: list[list[str]],
    n_reads: int,
    seed: int,
    status: str = "ok",
    bench_set: str = "standard",
    smear_frac: float = 0.0,
    err: float = 0.02,  # matches tests/unit/hybrid/test_lengths.py's proven spans() default
) -> dict[str, object]:
    case_dir = split_dir / design_id
    row: dict[str, object] = {
        "design_id": design_id,
        "profile": PROFILE,
        "status": status,
        "design": {"design_id": design_id, "bench_set": bench_set},
    }
    if status != "ok":
        return row
    _truth(case_dir, hap_units)
    reads_dir = case_dir / "reads"
    reads_dir.mkdir(parents=True)
    records: list[ReadRecord] = []
    for i, units in enumerate(hap_units):
        records += synth.reads(
            synth.allele(units), n_reads, err=err, seed=seed + i, smear_frac=smear_frac
        )
    _write_fastq(reads_dir / FASTQ[PROFILE].format(design_id), records)
    return row


def _manifest(split_dir: Path, rows: list[dict[str, object]]) -> Path:
    split_dir.mkdir(parents=True, exist_ok=True)
    path = split_dir / "manifest.jsonl"
    path.write_text("".join(json.dumps(r) + "\n" for r in rows))
    return path


def _run_and_evaluate(
    tmp_path: Path, rows: list[dict[str, object]], config: Path
) -> tuple[dict, int]:
    split_dir = tmp_path / "data" / "dev"
    manifest = _manifest(split_dir, rows)
    results = tmp_path / "results"
    run_lengths_stage(manifest, results, config, engine="hybrid")
    engine_dir = results / "hybrid"
    ns = argparse.Namespace(
        result_root=engine_dir,
        truth_root=split_dir,
        expected_samples=engine_dir / "inventory_dev.json",
    )
    return evaluate_lengths_stage(ns)


def test_run_lengths_stage_writes_peaks_and_inventory(tmp_path: Path) -> None:
    split_dir = tmp_path / "data" / "dev"
    rows = [_case(split_dir, "c1", [["X"] * 30], 150, seed=1)]
    manifest = _manifest(split_dir, rows)
    config = _write_config(tmp_path)
    results = tmp_path / "results"
    records = run_lengths_stage(manifest, results, config, engine="hybrid")
    assert [r["design_id"] for r in records] == ["c1"]
    data = json.loads((results / "hybrid" / "c1" / "lengths.json").read_text())
    assert data["status"] == "ok"
    assert data["unit_bp"] == UNIT
    (peak,) = data["peaks"]
    assert round(peak["center_bp"] / UNIT) == 39  # 5 pre + 30 inner + 4 post
    assert peak["support"] > 100
    inventory = json.loads((results / "hybrid" / "inventory_dev.json").read_text())
    assert inventory == [
        {
            "sample": "c1",
            "truth_dir": str(split_dir / "c1" / "truth"),
            "result_dir": str(results / "hybrid" / "c1"),
        }
    ]


def test_run_lengths_stage_records_a_manifest_status_that_is_not_ok(tmp_path: Path) -> None:
    split_dir = tmp_path / "data" / "dev"
    rows = [_case(split_dir, "c1", [], 0, seed=1, status="generation_failed")]
    manifest = _manifest(split_dir, rows)
    results = tmp_path / "results"
    run_lengths_stage(manifest, results, _write_config(tmp_path), engine="hybrid")
    data = json.loads((results / "hybrid" / "c1" / "lengths.json").read_text())
    assert data["status"] == "not_attempted"
    assert "generation_failed" in data["error"]


def test_run_lengths_stage_records_a_missing_reads_file(tmp_path: Path) -> None:
    split_dir = tmp_path / "data" / "dev"
    case_dir = split_dir / "c1"
    _truth(case_dir, [["X"] * 30])
    rows = [
        {
            "design_id": "c1",
            "profile": PROFILE,
            "status": "ok",
            "design": {"design_id": "c1", "bench_set": "standard"},
        }
    ]
    manifest = _manifest(split_dir, rows)
    results = tmp_path / "results"
    run_lengths_stage(manifest, results, _write_config(tmp_path), engine="hybrid")
    data = json.loads((results / "hybrid" / "c1" / "lengths.json").read_text())
    assert data["status"] == "not_attempted"
    assert "reads file not found" in data["error"]


def test_run_lengths_stage_records_an_unmapped_profile(tmp_path: Path) -> None:
    split_dir = tmp_path / "data" / "dev"
    rows = [
        {
            "design_id": "c1",
            "profile": "no_such_profile",
            "status": "ok",
            "design": {"design_id": "c1", "bench_set": "standard"},
        }
    ]
    manifest = _manifest(split_dir, rows)
    results = tmp_path / "results"
    run_lengths_stage(manifest, results, _write_config(tmp_path), engine="hybrid")
    data = json.loads((results / "hybrid" / "c1" / "lengths.json").read_text())
    assert data["status"] == "not_attempted"
    assert "no reads filename pattern" in data["error"]


def test_run_lengths_stage_catches_a_malformed_fastq(tmp_path: Path) -> None:
    split_dir = tmp_path / "data" / "dev"
    case_dir = split_dir / "c1"
    reads_dir = case_dir / "reads"
    reads_dir.mkdir(parents=True)
    (reads_dir / FASTQ[PROFILE].format("c1")).write_text("not a fastq record\n")
    rows = [
        {
            "design_id": "c1",
            "profile": PROFILE,
            "status": "ok",
            "design": {"design_id": "c1", "bench_set": "standard"},
        }
    ]
    manifest = _manifest(split_dir, rows)
    results = tmp_path / "results"
    run_lengths_stage(manifest, results, _write_config(tmp_path), engine="hybrid")
    data = json.loads((results / "hybrid" / "c1" / "lengths.json").read_text())
    assert data["status"] == "not_attempted"
    assert "ValueError" in data["error"]


def test_evaluate_records_a_run_that_was_not_attempted(tmp_path: Path) -> None:
    # Same "manifest status is not ok" fixture as
    # test_run_lengths_stage_records_a_manifest_status_that_is_not_ok, carried through
    # evaluate_lengths_stage: the case's own denominator-preserving failure, not a
    # missing-truth failure.
    split_dir = tmp_path / "data" / "dev"
    rows = [_case(split_dir, "c1", [], 0, seed=1, status="generation_failed")]
    report, code = _run_and_evaluate(tmp_path, rows, _write_config(tmp_path))
    assert code == 1
    (row,) = report["rows"]
    assert row["status"] == "not_attempted" and row["failed"]
    assert "generation_failed" in row["error"]


def test_evaluate_records_a_missing_lengths_file(tmp_path: Path) -> None:
    split_dir = tmp_path / "data" / "dev"
    rows = [_case(split_dir, "c1", [["X"] * 30], 20, seed=1)]
    manifest = _manifest(split_dir, rows)
    results = tmp_path / "results"
    run_lengths_stage(manifest, results, _write_config(tmp_path), engine="hybrid")
    (results / "hybrid" / "c1" / "lengths.json").unlink()
    ns = argparse.Namespace(
        result_root=results / "hybrid",
        truth_root=split_dir,
        expected_samples=results / "hybrid" / "inventory_dev.json",
    )
    report, code = evaluate_lengths_stage(ns)
    assert code == 1
    (row,) = report["rows"]
    assert row["status"] == "not_attempted" and row["failed"]
    assert "missing lengths.json" in row["error"]


def test_evaluate_requires_existing_result_and_truth_roots(tmp_path: Path) -> None:
    inventory = tmp_path / "inventory.json"
    inventory.write_text(json.dumps(["c1"]))
    ns = argparse.Namespace(
        result_root=tmp_path / "no-such-results",
        truth_root=tmp_path / "no-such-truth",
        expected_samples=inventory,
    )
    with pytest.raises(ValueError, match="must be existing directories"):
        evaluate_lengths_stage(ns)


def test_match_alleles_skips_a_peak_already_claimed_by_a_closer_truth_length() -> None:
    h = HybridSettings()
    # A single accepted peak sits within both truth lengths' windows; the closer truth
    # length claims it first, and the second (line 153's "already used" branch) is left
    # unmatched rather than double-assigned.
    peaks = [{"center_bp": 1000.0, "support": 50}]
    matched, false_alleles, missed_alleles = _match_alleles([995, 1015], peaks, h, UNIT)
    assert matched == [{"truth_bp": 995, "matched": 1}, {"truth_bp": 1015, "matched": 0}]
    assert false_alleles == 0
    assert missed_alleles == 1


def test_evaluate_scores_a_homozygous_case_as_a_single_exact_peak(tmp_path: Path) -> None:
    rows = [_case(tmp_path / "data" / "dev", "c1", [["X"] * 30], 150, seed=2)]
    report, code = _run_and_evaluate(tmp_path, rows, _write_config(tmp_path))
    assert code == 0
    (row,) = report["rows"]
    assert row["status"] == "ok" and not row["failed"]
    assert row["truth_allele_count"] == 1
    assert row["observed_allele_count"] == 1
    assert row["allele_count_exact"] == 1
    assert row["false_alleles"] == 0 and row["missed_alleles"] == 0
    assert row["case_length_exact"] == 1
    assert row["matched_alleles"] == [{"truth_bp": pytest.approx(39 * UNIT), "matched": 1}]
    assert row["reconstruction_flags"] == []


def test_evaluate_scores_a_heterozygous_case_with_two_matched_peaks(tmp_path: Path) -> None:
    rows = [_case(tmp_path / "data" / "dev", "c1", [["X"] * 30, ["X"] * 70], 150, seed=3)]
    report, code = _run_and_evaluate(tmp_path, rows, _write_config(tmp_path))
    assert code == 0
    (row,) = report["rows"]
    assert row["truth_allele_count"] == 2
    assert row["observed_allele_count"] == 2
    assert row["allele_count_exact"] == 1
    assert row["case_length_exact"] == 1
    assert {m["matched"] for m in row["matched_alleles"]} == {1}


def test_evaluate_records_a_missed_allele_when_a_minor_peak_is_rejected(tmp_path: Path) -> None:
    # Same fixture as tests/unit/hybrid/test_lengths.py::test_rejected_minor_peak_is_reported
    # (spans(30, 200, 6) + spans(80, 5, 7)): 5 reads at a second, distant length are too
    # few to clear the far-peak read-support threshold, so `fit_length_model` keeps only
    # the major peak and the truth's second haplotype is missed.
    split_dir = tmp_path / "data" / "dev"
    rows = [_case(split_dir, "c1", [["X"] * 30], 200, seed=6)]
    case_dir = split_dir / "c1"
    extra = synth.reads(synth.allele(["X"] * 80), 5, err=0.02, seed=7)
    fastq = case_dir / "reads" / FASTQ[PROFILE].format("c1")
    with fastq.open("a") as handle:
        for rec in extra:
            handle.write(f"@{rec.name}\n{rec.seq}\n+\n{rec.qual}\n")
    shutil.rmtree(case_dir / "truth")  # replace the single-haplotype truth _case() wrote
    _truth(case_dir, [["X"] * 30, ["X"] * 80])
    report, code = _run_and_evaluate(tmp_path, rows, _write_config(tmp_path))
    assert code == 0  # the case scored fine; it just did not match perfectly
    (row,) = report["rows"]
    assert row["status"] == "ok" and not row["failed"]
    assert row["truth_allele_count"] == 2
    assert row["observed_allele_count"] == 1
    assert row["missed_alleles"] == 1
    assert row["case_length_exact"] == 0


def test_evaluate_records_smear_ambiguous_as_a_reconstruction_flag(tmp_path: Path) -> None:
    # Same fixture as test_lengths.py::test_borderline_smear_candidate_is_smear_ambiguous
    # (spans(60, 120, 3, smear_frac=0.45)), with the same widened borderline band: a
    # proven scenario where the below-top candidate is "smear_ambiguous", not silent.
    rows = [_case(tmp_path / "data" / "dev", "c1", [["X"] * 60], 120, seed=3, smear_frac=0.45)]
    config = _write_config(tmp_path, {"hybrid.smear_test_borderline_factor": 1e6})
    report, _ = _run_and_evaluate(tmp_path, rows, config)
    (row,) = report["rows"]
    assert row["reconstruction_flags"] == ["smear_ambiguous"]


def test_evaluate_fails_a_case_whose_truth_cannot_be_loaded(tmp_path: Path) -> None:
    split_dir = tmp_path / "data" / "dev"
    rows = [_case(split_dir, "c1", [["X"] * 30], 20, seed=6)]
    (split_dir / "c1" / "truth" / "case.vntr_structure.txt").write_text("garbage\n")
    report, code = _run_and_evaluate(tmp_path, rows, _write_config(tmp_path))
    assert code == 1
    (row,) = report["rows"]
    assert row["status"] == "not_attempted" and row["failed"]


def _evaluate_with_window(tmp_path: Path, name: str, split_dir: Path, base_bp: float) -> dict:
    # ``result_root`` is ``<point_dir>/results/hybrid`` in the real calibration layout
    # (`calibration._run_point`), so `evaluate_lengths_stage` reads the point's own
    # ``config.json`` two levels up from ``result_root`` -- reproduce that layout here.
    point_dir = tmp_path / name
    manifest = split_dir / "manifest.jsonl"
    results = point_dir / "results"
    config = _write_config(
        point_dir, {"hybrid.peak_window_base_bp": base_bp, "hybrid.peak_window_per_unit_bp": 0}
    )
    run_lengths_stage(manifest, results, config, engine="hybrid")
    engine_dir = results / "hybrid"
    ns = argparse.Namespace(
        result_root=engine_dir,
        truth_root=split_dir,
        expected_samples=engine_dir / "inventory_dev.json",
    )
    report, _ = evaluate_lengths_stage(ns)
    (row,) = report["rows"]
    return row


def test_evaluate_uses_the_points_own_window_settings(tmp_path: Path) -> None:
    # Deterministic (zero read error) spanning reads for a 30-inner-unit allele, but a
    # truth declaring 31 inner units: a known, exact ~1-unit (60bp) gap between the
    # fitted peak and the truth length, so matching is decided purely by the window.
    split_dir = tmp_path / "data" / "dev"
    case_dir = split_dir / "c1"
    reads_dir = case_dir / "reads"
    reads_dir.mkdir(parents=True)
    _write_fastq(
        reads_dir / FASTQ[PROFILE].format("c1"),
        synth.reads(synth.allele(["X"] * 30), 20, err=0.0, seed=8),
    )
    _truth(case_dir, [["X"] * 31])
    rows = [
        {
            "design_id": "c1",
            "profile": PROFILE,
            "status": "ok",
            "design": {"design_id": "c1", "bench_set": "standard"},
        }
    ]
    _manifest(split_dir, rows)
    narrow = _evaluate_with_window(tmp_path, "narrow", split_dir, 1.0)
    wide = _evaluate_with_window(tmp_path, "wide", split_dir, 100.0)
    assert narrow["matched_alleles"][0]["matched"] == 0  # ~1bp window << ~60bp gap
    assert wide["matched_alleles"][0]["matched"] == 1  # 100bp window > ~60bp gap


def test_normalize_lengths_rows_joins_bench_set(tmp_path: Path) -> None:
    report = {"rows": [{"sample": "c1"}, {"sample": "c2"}]}
    cases = {"c1": {"design": {"bench_set": "clean"}}, "c2": {"design": {}}}
    rows = normalize_lengths_rows(report, cases, "legacy")
    assert [r["bench_set"] for r in rows] == ["clean", "legacy"]
    with pytest.raises(KeyError):
        normalize_lengths_rows({"rows": [{"sample": "missing"}]}, cases, "legacy")


def test_lengths_registries_feed_point_metrics_and_reason_metrics(tmp_path: Path) -> None:
    assert set(LENGTHS_METRICS) == set(LENGTHS_RATES) | set(LENGTHS_COUNTS)
    assert LENGTHS_METRICS["allele_count_exact"] == RATE
    objective_path = tmp_path / "objective.json"
    objective_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "constraints": {"missed_alleles": {"max": 0}},
                "rank": ["-allele_length_exact"],
                "reason_metrics": {"smear_ambiguous_rate": "smear_ambiguous"},
            }
        )
    )
    objective = load_objective(objective_path, known_metrics=LENGTHS_METRICS)
    rows = [
        {
            "sample": "c1",
            "bench_set": "standard",
            "failed": False,
            "allele_count_exact": 1,
            "case_length_exact": 1,
            "matched_alleles": [{"truth_bp": 100.0, "matched": 1}],
            "false_alleles": 0,
            "missed_alleles": 0,
            "reconstruction_flags": ["smear_ambiguous"],
            "clinical_reasons": [],
        }
    ]
    metrics = point_metrics(
        rows, objective, DEFAULT_BENCH_CONFIG.report, LENGTHS_RATES, LENGTHS_COUNTS
    )
    assert metrics["allele_length_exact"]["value"] == 1.0
    assert metrics["smear_ambiguous_rate"]["value"] == 1.0
    assert metrics["missed_alleles"]["value"] == 0
