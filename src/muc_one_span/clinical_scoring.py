"""Source-backed clinical endpoints, separate from simulated haplotype evaluation.

An event name is compared literally: a different alarm cannot recover a confirmed
mutation. Position and allele labels are not inferred from event-only truth. Run
artifacts are revalidated locally by the runner before granting callable credit.
The ledger records provenance assertions; checking that a cited publication really
supports each assertion remains an independent scientific-review responsibility.
"""

import hashlib
import json
from collections import Counter
from copy import deepcopy
from typing import Any
from urllib.parse import urlparse

from muc_one_span.evaluation.matching import match_alleles
from muc_one_span.evaluation.models import PredictedAllele, TruthHaplotype
from muc_one_span.evaluation.scoring import _recoverable_prediction_indices

Json = dict[str, Any]
ARMS = ("primary_amplicon", "secondary_wgs")
STATUSES = (
    "completed",
    "insufficient_evidence",
    "execution_failed",
    "invalid_artifacts",
    "unsupported_input",
    "unattempted",
)


def _text(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field} must be a nonempty string")
    return value


def _list(value: Any, field: str) -> list[Any]:
    if not isinstance(value, list):
        raise ValueError(f"{field} must be a list")
    return value


def _references(value: Any, sources: Json, field: str, *, independent: bool = False) -> None:
    refs = _list(value, field)
    if not refs or any(not isinstance(ref, str) or ref not in sources for ref in refs):
        raise ValueError(f"{field} needs known source IDs")
    if independent and not any(sources[ref]["independent"] for ref in refs):
        raise ValueError(f"{field} needs independent evidence")


def _object(value: Any, field: str) -> Json:
    if not isinstance(value, dict):
        raise ValueError(f"{field} must be an object")
    try:
        json.dumps(value, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must contain finite JSON values") from exc
    return value


def _inventory_runs(inventory: Json) -> Json:
    _object(inventory, "inventory")
    if inventory.get("schema_version") != 1:
        raise ValueError("unsupported inventory schema_version")
    runs = {}
    for run in _list(inventory.get("runs"), "inventory runs"):
        _object(run, "inventory run")
        accession = _text(run.get("run_accession"), "run_accession")
        if accession in runs:
            raise ValueError(f"duplicate inventory run: {accession}")
        if run.get("arm") not in (*ARMS, "excluded"):
            raise ValueError(f"unsupported inventory arm: {accession}")
        runs[accession] = run
    return runs


def _category(event: Json) -> str:
    return str(
        event.get(
            "evidence_category",
            "independent_confirmation"
            if event.get("status") in ("present", "absent")
            else "unknown",
        )
    )


def _sequence_truth(sequence: Any, sources: Json) -> None:
    if sequence is None:
        return
    if not isinstance(sequence, dict):
        raise ValueError("sequence_truth must be an object or null")
    for field in ("version", "coordinates", "boundaries"):
        _text(sequence.get(field), f"sequence {field}")
    _references(sequence.get("source_ids"), sources, "sequence source_ids", independent=True)
    if sequence.get("orientation") != "forward":
        raise ValueError("sequence orientation must be explicitly normalized to forward")
    if sequence.get("normalization") != "uppercase_allele_order":
        raise ValueError("unsupported sequence normalization")
    sequences = _list(sequence.get("sequences"), "sequence sequences")
    if not sequences or any(not _dna(item) for item in sequences):
        raise ValueError("sequence truth needs nonempty DNA sequences")
    hashes = [hashlib.sha256(item.upper().encode("ascii")).hexdigest() for item in sequences]
    if sequence.get("sha256") != hashes:
        raise ValueError("sequence checksum mismatch")
    anchors = sequence.get("boundary_anchors")
    if anchors is not None:
        _object(anchors, "boundary_anchors")
        _references(anchors.get("source_ids"), sources, "anchor source_ids", independent=True)
        if any(not _dna(anchors.get(side)) for side in ("left", "right")):
            raise ValueError("boundary anchors must be nonempty DNA")
        for item in sequences:
            _normalize_sequence(item, anchors)


def _normalize_sequence(sequence: str, anchors: Json | None) -> tuple[str, list[int]]:
    """Select a prespecified, unique, forward interval including both full anchors.

    No approximate matching is allowed. Missing, repeated or reversed anchors
    remain uncomparable rather than allowing convenient crop coordinates.
    """
    value = sequence.upper()
    if anchors is None:
        return value, [0, len(value)]
    positions = []
    for side in ("left", "right"):
        anchor = anchors[side].upper()
        start = value.find(anchor)
        if start < 0 or value.find(anchor, start + 1) >= 0:
            raise ValueError("boundary anchors must occur exactly once (unique exact match)")
        positions.append(start)
    start, right = positions
    if right < start + len(anchors["left"]):
        raise ValueError("boundary anchors overlap or have reversed order")
    end = right + len(anchors["right"])
    return value[start:end], [start, end]


def _sequence_comparison(
    sequence: Json | None, record: Json | None, validated: bool, callable_run: bool
) -> Json:
    result: Json = {
        "sequence_equal": None,
        "literal_sequence_equal": None,
        "sequence_length_equal": None,
        "sequence_length_error_bp": None,
        "sequence_assignments": None,
        "sequence_edit_distance": None,
        "literal_exact_alleles": None,
        "independent_exact_alleles": None,
        "sequence_intervals": None,
        "sequence_comparison_reason": "no independent sequence truth",
    }
    if sequence is None:
        return result
    result["sequence_comparison_reason"] = "run lacks validated completed artifacts"
    if record is None or not validated:
        return result
    predictions = record.get("sequences")
    result["sequence_comparison_reason"] = "missing sequence observations"
    if not isinstance(predictions, list) or not all(
        _dna(item, ambiguous=True) for item in predictions
    ):
        return result
    try:
        normalized = [
            _normalize_sequence(item, sequence.get("boundary_anchors")) for item in predictions
        ]
        expected = [
            _normalize_sequence(item, sequence.get("boundary_anchors"))[0]
            for item in sequence["sequences"]
        ]
    except ValueError as exc:
        result["sequence_comparison_reason"] = str(exc)
        return result
    values = [item[0] for item in normalized]
    flags = record.get("sequence_independent", [])
    sources = record.get("sequence_sources", [])
    alleles = [
        PredictedAllele(
            str(i),
            value,
            (),
            len(value),
            0,
            independent_haplotype_evidence=isinstance(flags, list)
            and i < len(flags)
            and flags[i] is True,
            sequence_source=sources[i] if isinstance(sources, list) and i < len(sources) else None,
        )
        for i, value in enumerate(values)
    ]
    assignment = match_alleles(
        [TruthHaplotype(str(i), value, ()) for i, value in enumerate(expected)], alleles
    )
    recoverable = _recoverable_prediction_indices(alleles)
    alternatives = [
        [
            {
                "truth_allele": ti,
                "predicted_allele": pi,
                "edit_distance": assignment.distance_matrix[ti][pi],
                "length_error_bp": len(values[pi]) - len(expected[ti]),
                "exact": expected[ti] == values[pi],
                "independent_recovery": pi in recoverable,
            }
            for ti, pi in pairs
        ]
        for pairs in assignment.alternatives
    ]
    observed_lengths, expected_lengths = sorted(map(len, values)), sorted(map(len, expected))
    literal = sorted(values) == sorted(expected)
    result.update(
        sequence_equal=literal and len(recoverable) == len(values) and callable_run,
        literal_sequence_equal=literal,
        sequence_length_equal=observed_lengths == expected_lengths,
        sequence_length_error_bp=[
            a - b for a, b in zip(observed_lengths, expected_lengths, strict=True)
        ]
        if len(values) == len(expected)
        else None,
        sequence_assignments=alternatives,
        sequence_edit_distance=assignment.cost,
        literal_exact_alleles=min(sum(pair["exact"] for pair in pairs) for pairs in alternatives),
        independent_exact_alleles=min(
            sum(pair["exact"] and pair["independent_recovery"] for pair in pairs)
            for pairs in alternatives
        )
        if callable_run
        else 0,
        sequence_intervals=[item[1] for item in normalized],
        sequence_comparison_reason=None,
    )
    return result


def _dna(sequence: Any, *, ambiguous: bool = False) -> bool:
    alphabet = "ACGTRYSWKMBDHVN" if ambiguous else "ACGT"
    return isinstance(sequence, str) and bool(sequence) and set(sequence.upper()) <= set(alphabet)


def validate_truth(data: Json, inventory: Json) -> Json:
    """Validate ledger structure and source eligibility without inferring truth.

    Every biological mapping is explicit. Independent confirmations require independent sources. Reported known
    controls are explicitly categorized and scored separately; comparator-only
    observations cannot establish presence or absence. Return
    a copy so normalization or downstream presentation cannot mutate input truth.
    """
    runs = _inventory_runs(inventory)
    _object(data, "truth")
    if data.get("schema_version") != 1:
        raise ValueError("unsupported truth schema_version")
    sources = {}
    for source in _list(data.get("sources"), "sources"):
        _object(source, "source")
        identifier = _text(source.get("id"), "source id")
        if identifier in sources:
            raise ValueError(f"duplicate source id: {identifier}")
        url = urlparse(_text(source.get("url"), "source url"))
        if url.scheme not in ("https", "http") or not url.netloc:
            raise ValueError("source url must be an HTTP(S) URL")
        _text(source.get("location"), "source location")
        if not isinstance(source.get("independent"), bool):
            raise ValueError("source independent must be explicit boolean")
        sources[identifier] = source
    identities: set[str] = set()
    assigned: set[str] = set()
    for sample in _list(data.get("samples"), "samples"):
        _object(sample, "sample")
        identity = _text(sample.get("biological_sample"), "biological_sample")
        if identity in identities:
            raise ValueError(f"duplicate biological sample: {identity}")
        identities.add(identity)
        if sample.get("identity_status", "resolved") not in ("resolved", "unresolved"):
            raise ValueError("unsupported identity_status")
        _references(sample.get("identity_source_ids"), sources, "identity_source_ids")
        accessions = _list(sample.get("run_accessions"), "run_accessions")
        if not accessions:
            raise ValueError("sample needs run_accessions")
        for accession in accessions:
            if not isinstance(accession, str) or accession not in runs:
                raise ValueError("truth run absent from inventory")
            if accession in assigned:
                raise ValueError(f"duplicate truth run mapping: {accession}")
            assigned.add(accession)
        _list(sample.get("relationships"), "relationships")
        _list(sample.get("limitations"), "limitations")
        events: set[str] = set()
        for event in _list(sample.get("event_truth"), "event_truth"):
            _object(event, "event truth")
            name = _text(event.get("event"), "event")
            if name in events:
                raise ValueError(f"duplicate event truth: {name}")
            events.add(name)
            status = event.get("status")
            if status not in ("present", "absent", "unknown"):
                raise ValueError("unsupported event truth status")
            category = _category(event)
            if category not in (
                "independent_confirmation",
                "reported_known_control",
                "comparator_only",
                "unknown",
            ):
                raise ValueError("unsupported evidence category")
            if status != "unknown" and category in ("comparator_only", "unknown"):
                raise ValueError("comparator/unknown category cannot assert confirmed truth")
            if event.get("resolution") != "event_identity":
                raise ValueError("unsupported event truth resolution")
            if any(event.get(field) is not None for field in ("repeat_index", "allele")):
                raise ValueError("event_identity truth cannot assert position or allele")
            refs = _list(event.get("source_ids"), "event source_ids")
            if status != "unknown" or refs:
                _references(
                    refs,
                    sources,
                    "event source_ids",
                    independent=status != "unknown" and category == "independent_confirmation",
                )
            if sample.get("identity_status") == "unresolved" and status != "unknown":
                raise ValueError("unresolved identity cannot carry confirmed truth")
        if "sequence_truth" not in sample:
            raise ValueError("sequence_truth must be explicit, including null")
        _sequence_truth(sample["sequence_truth"], sources)
        if sample.get("identity_status") == "unresolved" and sample["sequence_truth"] is not None:
            raise ValueError("unresolved identity cannot carry sequence truth")
    return deepcopy(data)


def _inventory_matches(run: Json, record: Json) -> bool:
    """Bind scoring to this inventory, beyond the record's self-consistency."""
    try:
        bound = record["provenance"]["run"]
        fields = (
            "run_accession",
            "sample_accession",
            "experiment_accession",
            "study_accession",
            "arm",
            "platform",
        )
        if any(run.get(field) != bound.get(field) for field in fields):
            return False
        expected = [
            {field: item[field] for field in ("url", "bytes", "md5")}
            for item in run.get("files", [])
        ]
        recorded = [
            {field: item[field] for field in ("url", "bytes", "md5")}
            for item in bound.get("files", [])
        ]
        return expected == recorded
    except (KeyError, TypeError, AttributeError):
        return False


def _validated_record(record: Json) -> bool:
    # Deferred import keeps clinical truth usable independently of execution.
    from muc_one_span.clinical_runner import validate_record

    return validate_record(record)


def _ratio(numerator: int, denominator: int, reason: str) -> Json:
    return {
        "numerator": numerator,
        "denominator": denominator,
        "value": numerator / denominator if denominator else None,
        "reason": None if denominator else reason,
    }


def _summary(rows: list[Json]) -> Json:
    counts = Counter(row["status"] for row in rows)
    return {
        "eligible_runs": len(rows),
        "biological_samples": len(
            {row["biological_sample"] for row in rows if row["identity_resolved"]}
        ),
        "status_counts": {status: counts[status] for status in STATUSES},
        "execution_success": _ratio(
            sum(row["execution_completed"] for row in rows), len(rows), "no eligible runs"
        ),
        "callability": _ratio(sum(row["callable"] for row in rows), len(rows), "no eligible runs"),
        "insufficient_evidence": _ratio(
            counts["insufficient_evidence"], len(rows), "no eligible runs"
        ),
    }


def _event_summary(
    rows: list[Json], event: str, category: str = "independent_confirmation"
) -> Json:
    label = (
        "independently confirmed"
        if category == "independent_confirmation"
        else "reported known-control"
    )
    eligible = [row for row in rows if row["event_results"][event]["evidence_category"] == category]
    positives = [row for row in eligible if row["event_results"][event]["truth"] == "present"]
    negatives = [row for row in eligible if row["event_results"][event]["truth"] == "absent"]
    callable_positive = [row for row in positives if row["callable"]]
    callable_negative = [row for row in negatives if row["callable"]]
    recovered = sum(row["event_results"][event]["recovered"] is True for row in positives)
    true_negative = sum(
        row["event_results"][event]["detected"] is False for row in callable_negative
    )
    biological_positive = {
        row["biological_sample"] for row in positives if row["identity_resolved"]
    }
    biological_recovered = {
        row["biological_sample"]
        for row in positives
        if row["identity_resolved"] and row["event_results"][event]["recovered"] is True
    }
    return {
        "biological_recovery_any_run": _ratio(
            len(biological_recovered),
            len(biological_positive),
            f"no {label} biological positives",
        ),
        "evidence_category": category,
        "eligible_positive_runs": len(positives),
        "eligible_negative_runs": len(negatives),
        "confirmed_positive_runs": len(positives) if category == "independent_confirmation" else 0,
        "confirmed_negative_runs": len(negatives) if category == "independent_confirmation" else 0,
        "unknown_truth_runs": len(rows) - len(positives) - len(negatives),
        "all_positive_recovery": _ratio(recovered, len(positives), f"no {label} positive truth"),
        "supported_positive_recovery": _ratio(
            sum(row["event_results"][event]["supported_recovered"] is True for row in positives),
            len(positives),
            f"no {label} positive truth",
        ),
        "callable_positive_recovery": _ratio(
            recovered, len(callable_positive), f"no callable {label} positives"
        ),
        "callable_specificity": _ratio(
            true_negative, len(callable_negative), f"no callable {label} negatives"
        ),
        "uncallable_positive_runs": len(positives) - len(callable_positive),
        "uncallable_negative_runs": len(negatives) - len(callable_negative),
    }


def _row(run: Json, sample: Json | None, record: Json | None, endpoints: list[str]) -> Json:
    status = record["status"] if record is not None else "unattempted"
    execution_completed = (
        status == "completed" and record is not None and record.get("exit_code") == 0
    )
    events = record.get("events") if record else None
    valid_events = isinstance(events, list) and all(
        isinstance(event, dict)
        and "event" in event
        and (event["event"] is None or (isinstance(event["event"], str) and event["event"].strip()))
        for event in events
    )
    if status == "completed" and (
        not execution_completed
        or not valid_events
        or not _inventory_matches(run, record or {})
        or not _validated_record(record or {})
    ):
        status = "invalid_artifacts"
    execution_completed = status == "completed"
    analysis_state = record.get("analysis_state") if record else None
    callable_run = status == "completed" and analysis_state == "completed"
    names = (
        {event["event"] for event in events or [] if event["event"] is not None}
        if callable_run
        else set()
    )
    assertions = {item["event"]: item for item in sample["event_truth"]} if sample else {}
    results = {}
    for endpoint in endpoints:
        assertion = assertions.get(endpoint, {})
        truth_status = assertion.get("status", "unknown")
        detected = endpoint in names if callable_run else None
        results[endpoint] = {
            "truth": truth_status,
            "evidence_category": _category(assertion),
            "source_ids": assertion.get("source_ids", []),
            "resolution": assertion.get("resolution"),
            "detected": detected,
            "recovered": detected is True if truth_status == "present" else None,
            "supported_recovered": any(
                item["event"] == endpoint and item.get("supported") is True for item in events or []
            )
            if callable_run and truth_status == "present"
            else False
            if truth_status == "present"
            else None,
            "wrong_positive_alarm": bool(events) and endpoint not in names
            if callable_run and truth_status == "present"
            else None,
        }
    sequence = sample["sequence_truth"] if sample else None
    return {
        "run_accession": run["run_accession"],
        "arm": run["arm"],
        "biological_sample": sample["biological_sample"] if sample else None,
        "identity_resolved": sample is not None
        and sample.get("identity_status", "resolved") == "resolved",
        "status": status,
        "recorded_status": record["status"] if record else None,
        "execution_completed": execution_completed,
        "callable": callable_run,
        "analysis_state": analysis_state,
        "events": deepcopy(events) if callable_run else None,
        "event_results": results,
        "sequence_truth_available": sequence is not None,
        **_sequence_comparison(sequence, record, status == "completed", callable_run),
        "limitations": sample["limitations"]
        if sample
        else ["No source-backed identity/truth mapping"],
    }


def score_cohort(inventory: Json, truth: Json, records: list[Json]) -> Json:
    """Score every eligible inventory run; retain missing, failed and unknown cases.

    All-positive recovery counts failed/no-call positives as unrecovered. Conditional
    specificity uses only callable, independently confirmed negatives. Replicate
    groups are displayed without treating repeated libraries as new participants.
    No statistical independence or clinical performance generalization is implied.
    """
    runs = _inventory_runs(inventory)
    ledger = validate_truth(truth, inventory)
    indexed = {}
    for record in records:
        _object(record, "record")
        accession = record.get("run_accession")
        if accession not in runs:
            raise ValueError("run record absent from inventory")
        if accession in indexed:
            raise ValueError(f"duplicate run record: {accession}")
        if record.get("schema_version") != 1 or record.get("status") not in STATUSES:
            raise ValueError("unsupported run record schema or status")
        if record.get("arm") != runs[accession]["arm"]:
            raise ValueError("record arm disagrees with inventory")
        indexed[accession] = record
    samples = {
        accession: sample for sample in ledger["samples"] for accession in sample["run_accessions"]
    }
    endpoints = sorted(
        {event["event"] for sample in ledger["samples"] for event in sample["event_truth"]}
    )
    rows = [
        _row(run, samples.get(accession), indexed.get(accession), endpoints)
        for accession, run in runs.items()
        if run["arm"] in ARMS
    ]
    arms = {}
    for arm in ARMS:
        selected = [row for row in rows if row["arm"] == arm]
        arms[arm] = {
            **_summary(selected),
            "events": {event: _event_summary(selected, event) for event in endpoints},
            "reported_controls": {
                event: _event_summary(selected, event, "reported_known_control")
                for event in endpoints
            },
        }
    groups = []
    for sample in ledger["samples"]:
        selected = [row for row in rows if row["biological_sample"] == sample["biological_sample"]]
        if selected:
            groups.append(
                {
                    "biological_sample": sample["biological_sample"],
                    "identity_status": sample.get("identity_status", "resolved"),
                    "run_accessions": [row["run_accession"] for row in selected],
                    "relationships": sample["relationships"],
                    "limitations": sample["limitations"],
                    "callable_runs": sum(row["callable"] for row in selected),
                    "events": {
                        event: {
                            "truth": selected[0]["event_results"][event]["truth"],
                            "evidence_category": selected[0]["event_results"][event][
                                "evidence_category"
                            ],
                            "recovered_runs": sum(
                                row["event_results"][event]["recovered"] is True for row in selected
                            ),
                        }
                        for event in endpoints
                    },
                }
            )
    sequence_rows = [row for row in rows if row["sequence_truth_available"]]
    return {
        "schema_version": 1,
        **_summary(rows),
        "arms": arms,
        "events": {event: _event_summary(rows, event) for event in endpoints},
        "reported_controls": {
            event: _event_summary(rows, event, "reported_known_control") for event in endpoints
        },
        "runs": rows,
        "biological_groups": groups,
        "unmapped_runs": [row["run_accession"] for row in rows if row["biological_sample"] is None],
        "unresolved_identity_runs": [
            row["run_accession"] for row in rows if not row["identity_resolved"]
        ],
        "excluded_runs": [accession for accession, run in runs.items() if run["arm"] == "excluded"],
        "sequence_concordance": _ratio(
            sum(row["sequence_equal"] is True for row in sequence_rows),
            len(sequence_rows),
            "no independently established full sequence truth",
        ),
        "literal_sequence_concordance": _ratio(
            sum(row["literal_sequence_equal"] is True for row in sequence_rows),
            len(sequence_rows),
            "no independently established full sequence truth",
        ),
        "sequence_length_concordance": _ratio(
            sum(row["sequence_length_equal"] is True for row in sequence_rows),
            len(sequence_rows),
            "no independently established full sequence truth",
        ),
        "sequence_comparable_runs": sum(row["sequence_equal"] is not None for row in sequence_rows),
        "sources": ledger["sources"],
        "provenance_boundary": "Callable records revalidated against local hashed artifacts; source interpretation requires independent scientific review.",
    }
