"""Denominator-aware sequence, ordered structure and event-annotation metrics."""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import asdict
from typing import Any

from muc_one_span.version import __version__

from .matching import match_alleles
from .models import Event, PredictedAllele, RunObservation, TruthSample

SCHEMA_VERSION = 1
SCORING_POLICY = {
    "diploid_evidence": "identical sequence copies require distinct source IDs and explicit independent genotype evidence for duplicate recovery credit and all-sample exactness",
    "assignment": "maximum cardinality then minimum summed literal global Levenshtein distance; all optima",
    "events": "event_annotation_exact: 1-based total-repeat index, exact parent and name; multiset matching",
    "supported_events": "frameshift AND template_match AND vcf_support are true AND vcf_support_status is exact_sequence_concordance; OR frameshift AND template_match AND read_support_status is supported (hybrid engine)",
    "legacy_supported_events": "historical comparator: frameshift AND template_match AND vcf_support booleans only",
    "accuracy_migration": "sequence_accuracy and structure_accuracy are evidence-aware; literal_*_accuracy retains the former raw assigned-pair numerator",
    "missing_alleles": "missing_alleles is literal output cardinality; independent_missing_alleles additionally excludes unproven duplicate observations",
    "ambiguity": "min/max across sequence-optimal assignments; favorable event labels never break ties",
    "identity_denominator": "max(truth bases, predicted bases); literal ambiguity symbols",
    "biological_event_equivalence": "not_assessable: normalized reference/alternate projection not implemented",
}
ACTIVE_STATUSES = ("completed", "insufficient_evidence", "ambiguous_reconstruction")


def score_events(expected: Sequence[Event], actual: Sequence[Event]) -> tuple[int, int, int]:
    """Return TP,FN,FP, preserving every duplicate and unmatched call."""
    truth, calls = Counter(expected), Counter(actual)
    return (
        sum((truth & calls).values()),
        sum((truth - calls).values()),
        sum((calls - truth).values()),
    )


def _ratio(numerator: int, denominator: int) -> dict[str, Any]:
    return {
        "numerator": numerator,
        "denominator": denominator,
        "value": numerator / denominator if denominator else None,
    }


def _recoverable_prediction_indices(predictions: Sequence[PredictedAllele]) -> set[int]:
    """Choose evidence-only representatives for unproven duplicate observations."""
    recoverable = set(range(len(predictions)))
    sequence_groups: dict[str, list[int]] = {}
    source_groups: dict[str, list[int]] = {}
    for index, prediction in enumerate(predictions):
        sequence_groups.setdefault(prediction.sequence, []).append(index)
        if prediction.sequence_source is not None:
            source_groups.setdefault(prediction.sequence_source, []).append(index)
    for indices in sequence_groups.values():
        sources = [predictions[index].sequence_source for index in indices]
        proven = (
            all(predictions[index].independent_haplotype_evidence for index in indices)
            and all(sources)
            and len(set(sources)) == len(sources)
        )
        if len(indices) > 1 and not proven:
            recoverable.difference_update(indices[1:])
    for indices in source_groups.values():
        if len(indices) > 1:
            recoverable.difference_update(indices[1:])
    return recoverable


def evaluate_sample(truth: TruthSample, observation: RunObservation) -> dict[str, Any]:
    """Score all optimal assignments, including failed samples in truth denominators."""
    predictions = observation.predictions if observation.status in ACTIVE_STATUSES else ()
    assignment = match_alleles(truth.haplotypes, predictions)
    recoverable = _recoverable_prediction_indices(predictions)
    independent = bool(predictions) and len(recoverable) == len(predictions)
    duplicates = len(predictions) - len(recoverable)
    warnings = list(truth.warnings + observation.warnings)
    if duplicates:
        warnings.append(f"Unproven duplicate allele observations: {duplicates}")
    alternatives: list[dict[str, Any]] = []
    total_events = sum(len(h.events) for h in truth.haplotypes)
    for pairs in assignment.alternatives:
        metrics = dict.fromkeys(
            (
                "sequence_exact",
                "structure_exact",
                "independent_sequence_exact",
                "independent_structure_exact",
                "count_exact",
                "count_within1",
                "count_within2",
                "event_tp",
                "event_fn",
                "event_fp",
                "supported_event_tp",
                "supported_event_fn",
                "supported_event_fp",
                "legacy_supported_event_tp",
                "legacy_supported_event_fn",
                "legacy_supported_event_fp",
            ),
            0,
        )
        details = []
        for ti, pi in pairs:
            t, p = truth.haplotypes[ti], predictions[pi]
            distance = assignment.distance_matrix[ti][pi]
            error = p.length - len(t.structure)
            for key, value in [
                ("sequence_exact", t.sequence == p.sequence),
                ("structure_exact", t.structure == p.structure),
                ("count_exact", error == 0),
                ("count_within1", abs(error) <= 1),
                ("count_within2", abs(error) <= 2),
            ]:
                metrics[key] += int(value)
            if pi in recoverable:
                metrics["independent_sequence_exact"] += int(t.sequence == p.sequence)
                metrics["independent_structure_exact"] += int(t.structure == p.structure)
            for prefix, events in [
                ("", p.events),
                ("supported_", tuple(e for e in p.events if e.supported)),
                ("legacy_supported_", tuple(e for e in p.events if e.legacy_supported)),
            ]:
                for suffix, count in zip(
                    ("tp", "fn", "fp"), score_events(t.events, events), strict=True
                ):
                    metrics[f"{prefix}event_{suffix}"] += count
            denominator = max(len(t.sequence), len(p.sequence))
            details.append(
                {
                    "truth": t.name,
                    "prediction": p.name,
                    "repeat_count_error": error,
                    "truth_repeats": len(t.structure),
                    "called_repeats": p.length,
                    "canonical_plus9_matches_reported": p.canonical_repeats + 9 == p.length,
                    "canonical_plus_fixed_matches_reported": p.canonical_repeats
                    + p.fixed_repeat_count
                    == p.length,
                    "fixed_repeat_count": p.fixed_repeat_count,
                    "sequence_edit_distance": distance,
                    "sequence_length_error": len(p.sequence) - len(t.sequence),
                    "sequence_identity": 1 - distance / denominator if denominator else None,
                    "ambiguous_bases": sum(base not in "ACGT" for base in p.sequence),
                    "sequence_exact": t.sequence == p.sequence,
                    "structure_exact": t.structure == p.structure,
                    "truth_structure": list(t.structure),
                    "predicted_structure": list(p.structure),
                }
            )
        missing = [h for i, h in enumerate(truth.haplotypes) if i not in {ti for ti, _ in pairs}]
        extra = [p for i, p in enumerate(predictions) if i not in {pi for _, pi in pairs}]
        for prefix in ("", "supported_", "legacy_supported_"):
            metrics[f"{prefix}event_fn"] += sum(len(t.events) for t in missing)
            if not prefix:
                extra_events = sum(len(p.events) for p in extra)
            elif prefix == "supported_":
                extra_events = sum(sum(e.supported for e in p.events) for p in extra)
            else:
                extra_events = sum(sum(e.legacy_supported for e in p.events) for p in extra)
            metrics[f"{prefix}event_fp"] += extra_events
        cardinality = len(truth.haplotypes) == len(predictions)
        for suffix, key in [
            ("sequences", "sequence_exact"),
            ("structures", "structure_exact"),
            ("counts", "count_exact"),
            ("counts_within1", "count_within1"),
            ("counts_within2", "count_within2"),
        ]:
            metrics[f"all_{suffix}_exact"] = int(
                cardinality and independent and metrics[key] == len(truth.haplotypes)
            )
        alternatives.append(
            {
                "pairs": details,
                "missing_truth": [h.name for h in missing],
                "extra_predictions": [p.name for p in extra],
                "metrics": metrics,
            }
        )
    bounds = {
        key: {
            "min": min(a["metrics"][key] for a in alternatives),
            "max": max(a["metrics"][key] for a in alternatives),
        }
        for key in alternatives[0]["metrics"]
    }
    positive = any(p.events for p in predictions)
    supported_positive = any(e.supported for p in predictions for e in p.events)
    legacy_supported_positive = any(e.legacy_supported for p in predictions for e in p.events)
    # A normal no-call cannot become a TN through absence of candidate events.
    assessed_normal = (
        total_events == 0
        and observation.status == "completed"
        and independent
        and len(predictions) == len(truth.haplotypes)
    )
    return {
        "sample": truth.name,
        "independent_haplotype_evidence": independent,
        "prediction_evidence": {
            p.name: {
                "sequence_source": p.sequence_source,
                "independent_haplotype_evidence": p.independent_haplotype_evidence,
                "phase_status": p.phase_status,
                "genotype_status": p.genotype_status,
                "sequence_identity_status": p.sequence_identity_status,
            }
            for p in predictions
        },
        "support_status_counts": dict(
            Counter(e.support_status for p in predictions for e in p.events)
        ),
        "status": observation.status,
        "truth_status": "valid",
        "truth_haplotypes": len(truth.haplotypes),
        "truth_events": total_events,
        "predicted_alleles": len(predictions),
        "missing_alleles": max(0, len(truth.haplotypes) - len(predictions)),
        "independent_missing_alleles": max(0, len(truth.haplotypes) - len(recoverable)),
        "unproven_duplicate_alleles": duplicates,
        "extra_alleles": max(0, len(predictions) - len(truth.haplotypes)),
        "assignment_cost": assignment.cost,
        "assignment_ties": len(alternatives),
        "distance_matrix": assignment.distance_matrix,
        "alternatives": alternatives,
        "metrics": bounds,
        "sample_alarm": positive,
        "supported_sample_alarm": supported_positive,
        "legacy_supported_sample_alarm": legacy_supported_positive,
        "supported_normal_true_negative": assessed_normal and not supported_positive,
        "supported_normal_false_positive": total_events == 0 and supported_positive,
        "legacy_supported_normal_true_negative": assessed_normal and not legacy_supported_positive,
        "legacy_supported_normal_false_positive": total_events == 0 and legacy_supported_positive,
        "normal_assessable": assessed_normal,
        "normal_true_negative": assessed_normal and not positive,
        "normal_false_positive": total_events == 0 and positive,
        "mutation_positive_truth": total_events > 0,
        "event_annotation_assessable": all(
            e.repeat_index is not None and e.parent is not None and e.name is not None
            for p in predictions
            for e in p.events
        ),
        "run_record": observation.run_record,
        "warnings": warnings,
        "truth_provenance": truth.provenance,
        "truth_hashes": truth.hashes,
        "prediction_events": {p.name: [asdict(e) for e in p.events] for p in predictions},
        "error": observation.error,
    }


def aggregate(samples: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Aggregate conservative endpoints with explicit denominators and null ratios."""
    valid = [s for s in samples if s.get("truth_status") == "valid"]
    totals: dict[str, Any] = {
        "samples": len(samples),
        "valid_truth_samples": len(valid),
        "invalid_truth_samples": len(samples) - len(valid),
        "statuses": dict(Counter(s["status"] for s in samples)),
    }
    for key in (
        "truth_haplotypes",
        "truth_events",
        "predicted_alleles",
        "missing_alleles",
        "independent_missing_alleles",
        "unproven_duplicate_alleles",
        "extra_alleles",
    ):
        totals[key] = sum(s[key] for s in valid)
    metrics = (
        {
            key: {bound: sum(s["metrics"][key][bound] for s in valid) for bound in ("min", "max")}
            for key in valid[0]["metrics"]
        }
        if valid
        else {}
    )
    totals["metrics"] = metrics

    def bound(key: str, endpoint: str = "min") -> int:
        return int(metrics.get(key, {}).get(endpoint, 0))

    for prefix in ("", "supported_", "legacy_supported_"):
        tp, fp = bound(f"{prefix}event_tp"), bound(f"{prefix}event_fp", "max")
        totals[f"{prefix}event_recall"] = _ratio(tp, totals["truth_events"])
        totals[f"{prefix}event_precision"] = _ratio(tp, tp + fp)
    totals["literal_sequence_accuracy"] = _ratio(
        bound("sequence_exact"), totals["truth_haplotypes"]
    )
    totals["literal_structure_accuracy"] = _ratio(
        bound("structure_exact"), totals["truth_haplotypes"]
    )
    totals["sequence_accuracy"] = _ratio(
        bound("independent_sequence_exact"), totals["truth_haplotypes"]
    )
    totals["structure_accuracy"] = _ratio(
        bound("independent_structure_exact"), totals["truth_haplotypes"]
    )
    totals["independently_recovered_sequence_accuracy"] = totals["sequence_accuracy"]
    totals["independently_recovered_structure_accuracy"] = totals["structure_accuracy"]
    totals["exact_sample_reconstruction"] = _ratio(bound("all_sequences_exact"), len(valid))
    totals["call_rate"] = _ratio(sum(s["predicted_alleles"] > 0 for s in valid), len(valid))
    totals["no_call_rate"] = _ratio(sum(s["predicted_alleles"] == 0 for s in valid), len(valid))
    matched = sum(min(s["predicted_alleles"], s["truth_haplotypes"]) for s in valid)
    totals["literal_conditional_sequence_accuracy"] = _ratio(bound("sequence_exact"), matched)
    totals["conditional_sequence_accuracy"] = _ratio(bound("independent_sequence_exact"), matched)
    totals["independently_recovered_conditional_sequence_accuracy"] = totals[
        "conditional_sequence_accuracy"
    ]
    totals["supported_normal_specificity"] = _ratio(
        sum(s["supported_normal_true_negative"] for s in valid),
        sum(
            s["supported_normal_true_negative"] or s["supported_normal_false_positive"]
            for s in valid
        ),
    )
    totals["supported_sample_alarm_sensitivity"] = _ratio(
        sum(s["mutation_positive_truth"] and s["supported_sample_alarm"] for s in valid),
        sum(s["mutation_positive_truth"] for s in valid),
    )
    totals["legacy_supported_normal_specificity"] = _ratio(
        sum(s["legacy_supported_normal_true_negative"] for s in valid),
        sum(
            s["legacy_supported_normal_true_negative"]
            or s["legacy_supported_normal_false_positive"]
            for s in valid
        ),
    )
    totals["legacy_supported_sample_alarm_sensitivity"] = _ratio(
        sum(s["mutation_positive_truth"] and s["legacy_supported_sample_alarm"] for s in valid),
        sum(s["mutation_positive_truth"] for s in valid),
    )
    for prefix in ("", "supported_", "legacy_supported_"):
        alarms = sum(s[f"{prefix}sample_alarm"] for s in valid)
        true_alarms = sum(
            s["mutation_positive_truth"] and s[f"{prefix}sample_alarm"] for s in valid
        )
        totals[f"{prefix}sample_alarm_precision"] = _ratio(true_alarms, alarms)
    totals["normal_specificity"] = _ratio(
        sum(s["normal_true_negative"] for s in valid),
        sum(s["normal_true_negative"] or s["normal_false_positive"] for s in valid),
    )
    normal_count = sum(not s["mutation_positive_truth"] for s in valid)
    totals["normal_confident_negative_rate"] = _ratio(
        sum(s["normal_true_negative"] for s in valid), normal_count
    )
    totals["normal_unresolved_negative_count"] = sum(
        not s["mutation_positive_truth"] and not s["normal_true_negative"] and not s["sample_alarm"]
        for s in valid
    )
    totals["resolved_reconstruction_call_rate"] = _ratio(
        sum(
            s["status"] == "completed"
            and s["missing_alleles"] == 0
            and s["extra_alleles"] == 0
            and s["independent_haplotype_evidence"]
            for s in valid
        ),
        len(valid),
    )
    totals["sample_alarm_sensitivity"] = _ratio(
        sum(s["mutation_positive_truth"] and s["sample_alarm"] for s in valid),
        sum(s["mutation_positive_truth"] for s in valid),
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "package_version": __version__,
        "scoring_policy": SCORING_POLICY,
        "totals": totals,
        "samples": list(samples),
    }
