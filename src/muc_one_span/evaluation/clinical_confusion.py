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


def predicted_decision(result_dir: Path) -> str:
    """Return the caller's 3-state clinical decision, or ``NO_CALL`` when unreadable."""
    try:
        summary = json.loads((result_dir / "summary.json").read_text())
    except (OSError, ValueError):
        return "NO_CALL"
    try:
        return str(compute_clinical_decision(summary)["state"])
    except (AttributeError, KeyError, TypeError):
        return "NO_CALL"


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
