"""Clinical-decision confusion matrix scoring simulated truth against caller calls."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from muc_one_span.config import RepeatDictionary, _apply_mutation
from muc_one_span.evaluation.models import TruthSample
from muc_one_span.report import compute_clinical_decision

_PLACEHOLDER_LENGTH = 60
_FALSE_NEGATIVE_DECISIONS = ("NO_PATHOGENIC_VARIANT_DETECTED", "NO_CALL")


def net_length_change(definition: dict[str, Any]) -> int:
    """Return the net base-count change ``_apply_mutation`` produces for *definition*.

    Applies ``definition["changes"]`` to a neutral 60-base placeholder sequence with
    the same helper the caller's classifier uses (``muc_one_span.config._apply_mutation``),
    so insert/delete/delete_insert length arithmetic never drifts from that ground truth.
    """
    placeholder = "A" * _PLACEHOLDER_LENGTH
    mutated = _apply_mutation(placeholder, definition["changes"])
    return len(mutated) - len(placeholder)


def truth_class(truth: TruthSample, rd: RepeatDictionary) -> str:
    """Classify simulated truth as ``pathogenic``, ``benign``, or ``normal``.

    An event is pathogenic when its net length change (``net_length_change``) is not
    divisible by 3 (frameshift). Truth is pathogenic if any event is, benign if only
    in-frame events exist, and normal when there are no events.
    """
    names = [event.name for hap in truth.haplotypes for event in hap.events if event.name]
    if not names:
        return "normal"
    frameshifts = [net_length_change(rd.mutations[name]) % 3 != 0 for name in names]
    return "pathogenic" if any(frameshifts) else "benign"


def predicted_clinical(result_dir: Path) -> dict[str, Any]:
    """Return the caller's clinical decision and its banner details as ``reasons``.

    ``reasons`` is ``compute_clinical_decision(summary)["details"]`` for every
    decision: the gate reasons of an INCONCLUSIVE banner, the findings and
    quality caveats of a PATHOGENIC one, the allele lines of a negative one.
    An unreadable summary is ``NO_CALL`` without reasons.
    """
    no_call: dict[str, Any] = {"decision": "NO_CALL", "reasons": []}
    try:
        # RecursionError: json.loads uses a recursive parser, so a pathologically
        # nested-but-syntactically-valid document (e.g. thousands of "[") can blow the
        # interpreter's recursion limit instead of raising a JSON decode error.
        summary = json.loads((result_dir / "summary.json").read_text())
    except (OSError, ValueError, RecursionError):
        return no_call
    try:
        # A malformed-but-valid-JSON shape (e.g. a list instead of an object) always
        # surfaces here as one of these three types: AttributeError (no .get on a
        # non-dict), KeyError (missing "state"/expected key), or TypeError (wrong
        # argument shape passed further down) — checked against compute_clinical_decision.
        decision = compute_clinical_decision(summary)
    except (AttributeError, KeyError, TypeError):
        return no_call
    return {"decision": str(decision["state"]), "reasons": [str(r) for r in decision["details"]]}


def predicted_decision(result_dir: Path) -> str:
    """Return the caller's 3-state clinical decision, or ``NO_CALL`` when unreadable."""
    return str(predicted_clinical(result_dir)["decision"])


def confusion(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """Tally truth-class x decision counts and clinical-safety error rates.

    Every row's ``clinical`` sub-dict is counted, including rows with a ``None`` truth
    (invalid or unattempted runs), so denominators are never dropped.
    """
    matrix: dict[Any, Counter[str]] = defaultdict(Counter)
    inconclusive = 0
    critical_false_negative = 0
    false_positive_normal = 0
    false_positive_benign = 0
    for row in rows:
        clinical = row["clinical"]
        truth, decision = clinical["truth"], clinical["decision"]
        matrix[truth][decision] += 1
        if decision == "INCONCLUSIVE":
            inconclusive += 1
        if truth == "pathogenic" and decision in _FALSE_NEGATIVE_DECISIONS:
            critical_false_negative += 1
        if decision == "PATHOGENIC":
            if truth == "normal":
                false_positive_normal += 1
            elif truth == "benign":
                false_positive_benign += 1
    return {
        "matrix": {key: dict(counts) for key, counts in matrix.items()},
        "critical_false_negative": critical_false_negative,
        "false_positive_normal": false_positive_normal,
        "false_positive_benign": false_positive_benign,
        "inconclusive_rate": inconclusive / len(rows) if rows else None,
    }
