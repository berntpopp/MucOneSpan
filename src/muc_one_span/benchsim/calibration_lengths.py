"""Task 15d: ``benchsim calibrate --stage lengths``, a fast length-model-only stage.

Fits only the hybrid length model (`hybrid.spans.categorize_reads` + S1 anchors,
`hybrid.lengths.fit_length_model` S2, `hybrid.smear` significance test) on each
case's spanning reads and scores the accepted peaks against the case truth: no
consensus, phasing or calling runs, so a smear/peak threshold sweep takes seconds
instead of a full-pipeline run per point.

Plugs into the unchanged 15b/15c machinery as a `calibration.RunFn`/`EvaluateFn`
pair (`run_lengths_stage`, `evaluate_lengths_stage`), written to and read from the
same ``<point>/results/<engine>/`` layout `run_cases.run_split` and
`scripts/evaluate.py` use, so `calibration.run_calibration`'s validation, content
addressing and resume logic need no change. `LENGTHS_RATES`/`LENGTHS_COUNTS` are a
second metric registry for the already-generic `calibration_objective.point_metrics`
(constraints, ranking and `calibrate-report`'s report/recommendation files are
unchanged; only the metrics a lengths-stage objective may reference differ).
"""

from __future__ import annotations

import argparse
import itertools
import json
from collections.abc import Callable, Iterator, Sequence
from pathlib import Path
from typing import Any

from muc_one_span.benchsim.calibration import dump, write_atomic
from muc_one_span.benchsim.calibration_grid import OVERLAY_FILE
from muc_one_span.benchsim.calibration_objective import COUNT, RATE
from muc_one_span.benchsim.generate import FASTQ
from muc_one_span.config import load_repeat_dictionary
from muc_one_span.evaluation.artifacts import read_inventory
from muc_one_span.evaluation.truth import TruthValidationError, load_truth
from muc_one_span.hybrid.lengths import fit_length_model, window_bp
from muc_one_span.hybrid.reads_io import read_input
from muc_one_span.hybrid.spans import Anchors, categorize_reads
from muc_one_span.settings import HybridSettings, load_settings

SCHEMA_VERSION = 1
LENGTHS_FILE = "lengths.json"
NOT_ATTEMPTED = "not_attempted"
Row = dict[str, Any]


def _read_manifest(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _fit_one(row: dict[str, Any], split_dir: Path, anchors: Anchors, h: HybridSettings) -> Row:
    """Length-model output for one manifest row, or why it was not attempted."""
    design_id = row["design_id"]
    profile = row.get("profile")
    if row.get("status") != "ok":
        reason = row.get("error") or f"manifest status: {row.get('status')}"
        return {"status": NOT_ATTEMPTED, "error": reason}
    pattern = FASTQ.get(profile) if isinstance(profile, str) else None
    if pattern is None:
        return {
            "status": NOT_ATTEMPTED,
            "error": f"no reads filename pattern for profile {profile!r}",
        }
    fastq = split_dir / design_id / "reads" / pattern.format(design_id)
    if not fastq.is_file():
        return {"status": NOT_ATTEMPTED, "error": f"reads file not found: {fastq}"}
    try:
        cats = categorize_reads(read_input(fastq), anchors, h)
        model = fit_length_model(cats.spanning, h, anchors.unit_bp)
    except (OSError, ValueError) as exc:
        return {"status": NOT_ATTEMPTED, "error": f"{type(exc).__name__}: {exc}"}
    return {
        "status": "ok",
        "unit_bp": anchors.unit_bp,
        "peaks": [{"center_bp": p.center_bp, "support": p.support} for p in model.peaks],
        "rejected": model.rejected,
        "read_categories": cats.counts(),
    }


def run_lengths_stage(
    manifest: Path, results_root: Path, config: Path, *, engine: str
) -> list[dict[str, Any]]:
    """`calibration.RunFn` for ``--stage lengths``: fit the length model, no truth read.

    Mirrors `run_cases.run_split`'s I/O contract: per-case output under
    ``results_root/<engine>/<design_id>/`` plus an ``inventory_<split>.json`` that
    `evaluate_lengths_stage` (which alone sees the truth root) reads.
    """
    manifest = Path(manifest).resolve()
    split_dir = manifest.parent
    split = split_dir.name
    rows = _read_manifest(manifest)
    h = load_settings(config).hybrid
    anchors = Anchors.from_dictionary(load_repeat_dictionary(), h)
    engine_dir = Path(results_root) / engine
    engine_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for row in rows:
        design_id = row["design_id"]
        out_dir = engine_dir / design_id
        out_dir.mkdir(parents=True, exist_ok=True)
        record = _fit_one(row, split_dir, anchors, h)
        write_atomic(out_dir / LENGTHS_FILE, dump(record))
        records.append({"design_id": design_id, **record})
    inventory = [
        {
            "sample": row["design_id"],
            "truth_dir": str(split_dir / row["design_id"] / "truth"),
            "result_dir": str(engine_dir / row["design_id"]),
        }
        for row in rows
    ]
    write_atomic(engine_dir / f"inventory_{split}.json", dump(inventory))
    return records


def _failed_row(name: str, error: str) -> Row:
    return {
        "sample": name,
        "status": NOT_ATTEMPTED,
        "failed": True,
        "error": error,
        "truth_allele_count": None,
        "observed_allele_count": None,
        "allele_count_exact": 0,
        "matched_alleles": [],
        "false_alleles": 0,
        "missed_alleles": 0,
        "case_length_exact": 0,
        "reconstruction_flags": [],
        "clinical_reasons": [],
    }


Assignment = tuple[int | None, ...]  # per truth index, its peak index or None (unmatched)


def _candidate_assignments(n_truth: int, n_peaks: int) -> Iterator[Assignment]:
    """Every truth-index -> peak-index-or-None assignment that uses each peak once.

    Exhaustive: sizes are tiny (at most `PLOIDY` truth lengths and accepted peaks), and
    the search space (``(n_peaks + 1) ** n_truth``, both counts taken from the inputs)
    is not bounded by any literal here.
    """
    for assignment in itertools.product((None, *range(n_peaks)), repeat=n_truth):
        used = [p for p in assignment if p is not None]
        if len(used) == len(set(used)):
            yield assignment


def _assignment_score(
    assignment: Assignment, distance: list[list[float]]
) -> tuple[int, float, tuple[int, ...]]:
    """Sort key picking, in order: fewest unmatched, then least total distance, then a
    fixed deterministic tie-break (the assignment itself, ``None`` sorting first)."""
    unmatched = sum(p is None for p in assignment)
    total_distance = sum(distance[ti][p] for ti, p in enumerate(assignment) if p is not None)
    tie_break = tuple(-1 if p is None else p for p in assignment)
    return unmatched, total_distance, tie_break


def _match_alleles(
    truth_lengths: list[int], peaks: list[dict[str, Any]], h: HybridSettings, unit_bp: int
) -> tuple[list[dict[str, Any]], int, int]:
    """Exact 1:1 match of truth lengths to peaks within `window_bp`: maximum matched
    pairs first, minimum total distance among those, a fixed deterministic tie-break.

    A (truth, peak) pair is a candidate only when their distance is inside the truth
    length's own assignment window (`hybrid.lengths.window_bp`, the same tolerance the
    engine itself uses to assign a read to a peak). A greedy nearest-first match is not
    maximum-cardinality (a truth length reachable through only one peak can be
    stranded by a closer truth length competing for the same peak, when routing it
    through a different peak would have matched both); this instead searches every
    candidate assignment exhaustively (`_candidate_assignments`) and keeps the best by
    `_assignment_score`.
    """
    distance = [[abs(peak["center_bp"] - t) for peak in peaks] for t in truth_lengths]
    windows = [window_bp(t, h, unit_bp) for t in truth_lengths]

    def valid(assignment: Assignment) -> bool:
        return all(p is None or distance[ti][p] <= windows[ti] for ti, p in enumerate(assignment))

    best = min(
        (a for a in _candidate_assignments(len(truth_lengths), len(peaks)) if valid(a)),
        key=lambda a: _assignment_score(a, distance),
    )
    matched_alleles = [
        {"truth_bp": t, "matched": int(best[i] is not None)} for i, t in enumerate(truth_lengths)
    ]
    used_peaks = {p for p in best if p is not None}
    return matched_alleles, len(peaks) - len(used_peaks), len(truth_lengths) - len(used_peaks)


def _score(
    name: str, truth_dir: Path, result_dir: Path, rd: Any, h: HybridSettings, unit_bp: int
) -> Row:
    lengths_path = result_dir / LENGTHS_FILE
    if not lengths_path.is_file():
        return _failed_row(name, f"missing {lengths_path.name}")
    data = json.loads(lengths_path.read_text())
    if data.get("status") != "ok":
        return _failed_row(name, str(data.get("error")))
    try:
        truth = load_truth(truth_dir, rd)
    except TruthValidationError as exc:
        return _failed_row(name, str(exc))
    truth_lengths = sorted({len(hap.sequence) for hap in truth.haplotypes})
    matched_alleles, false_alleles, missed_alleles = _match_alleles(
        truth_lengths, data["peaks"], h, unit_bp
    )
    smear_ambiguous = any(r["reason"] == "smear_ambiguous" for r in data["rejected"])
    return {
        "sample": name,
        "status": "ok",
        "failed": False,
        "truth_allele_count": len(truth_lengths),
        "observed_allele_count": len(data["peaks"]),
        "allele_count_exact": int(len(data["peaks"]) == len(truth_lengths)),
        "matched_alleles": matched_alleles,
        "false_alleles": false_alleles,
        "missed_alleles": missed_alleles,
        "case_length_exact": int(false_alleles == 0 and missed_alleles == 0),
        "reconstruction_flags": ["smear_ambiguous"] if smear_ambiguous else [],
        "clinical_reasons": [],
    }


def evaluate_lengths_stage(ns: argparse.Namespace) -> tuple[dict[str, Any], int]:
    """`calibration.EvaluateFn` for ``--stage lengths``: score peaks against truth.

    Reads the point's own settings from ``<point>/config.json``, two levels above
    ``ns.result_root`` (``<point>/results/<engine>``, the layout
    `calibration._run_point` writes), so the peak/truth matching tolerance
    (`hybrid.lengths.window_bp`) is the point's own, not a hardcoded default.
    """
    inventory = read_inventory(ns.expected_samples)
    if not ns.result_root.is_dir() or not ns.truth_root.is_dir():
        raise ValueError("result root and truth root must be existing directories")
    h = load_settings(ns.result_root.parents[1] / OVERLAY_FILE).hybrid
    rd = load_repeat_dictionary()
    unit_bp = rd.repeat_length_bp
    rows, failed = [], False
    for entry in inventory:
        name = entry["sample"]
        truth_dir = Path(entry.get("truth_dir", ns.truth_root / name / "truth"))
        result_dir = Path(entry.get("result_dir", ns.result_root / name))
        row = _score(name, truth_dir, result_dir, rd, h, unit_bp)
        rows.append(row)
        failed |= bool(row["failed"])
    report = {"schema_version": SCHEMA_VERSION, "rows": rows}
    return report, int(failed)


def _all(rows: Sequence[Row]) -> list[Row]:
    return list(rows)


def _allele_length_rows(rows: Sequence[Row]) -> list[Row]:
    """One row per (case, distinct truth allele length): the `allele_length_exact` cohort."""
    return [
        {"sample": row["sample"], "bench_set": row["bench_set"], **allele}
        for row in rows
        for allele in row["matched_alleles"]
    ]


# Built-in lengths-stage metrics for `calibration_objective.point_metrics`'s injectable
# ``rates``/``counts`` registries (brief: allele count, allele lengths, false/missed
# alleles). `smear_ambiguous_rate` is not built in: an objective declares it as a
# `reason_metrics` entry against `reconstruction_flags`, exactly like the full pipeline.
LENGTHS_RATES: dict[str, tuple[Callable[[Sequence[Row]], list[Row]], str]] = {
    "allele_count_exact": (_all, "allele_count_exact"),
    "allele_length_exact": (_allele_length_rows, "matched"),
    "case_length_exact": (_all, "case_length_exact"),
}
LENGTHS_COUNTS: dict[str, str | None] = {
    "cases": None,
    "false_alleles": "false_alleles",
    "missed_alleles": "missed_alleles",
    "not_completed": "failed",
}
LENGTHS_METRICS = {**dict.fromkeys(LENGTHS_RATES, RATE), **dict.fromkeys(LENGTHS_COUNTS, COUNT)}


def normalize_lengths_rows(
    report: dict[str, Any], cases: dict[str, dict[str, Any]], legacy_set: str
) -> list[Row]:
    """Case rows with ``bench_set`` joined from the manifest (like `report.normalize_rows`).

    Raises:
        KeyError: If a sample has no case (a stratum would silently be lost).
    """
    out = []
    for row in report["rows"]:
        design = cases[row["sample"]].get("design") or {}
        out.append({**row, "bench_set": design.get("bench_set") or legacy_set})
    return out
