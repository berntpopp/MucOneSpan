"""Clinical-decision confusion matrix: simulated truth class vs. caller decision."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.evaluation.clinical_confusion import (
    confusion,
    net_length_change,
    predicted_decision,
    truth_class,
)
from muc_one_span.evaluation.models import Event, TruthHaplotype, TruthSample

RD = load_repeat_dictionary()


def _truth(*names: str) -> TruthSample:
    events = tuple(Event(5, "X", n) for n in names)
    return TruthSample(
        "s", (TruthHaplotype("h1", "A", ("X",), events), TruthHaplotype("h2", "A", ("X",)))
    )


def test_truth_class_by_frameshift() -> None:
    assert truth_class(_truth(), RD) == "normal"
    assert truth_class(_truth("dupC"), RD) == "pathogenic"
    assert net_length_change(RD.mutations["dupC"]) == 1


def test_truth_class_benign_in_frame_event() -> None:
    """Every real repeats.json mutation is frameshift (net_length_change % 3 == 1 or 2 for
    all 13 entries, verified against ``_apply_mutation`` directly); the dictionary curates
    only pathogenic ADTKD-MUC1 variants and has no in-frame entry. The benign path is
    therefore proven with a synthetic in-frame definition of the same change schema,
    merged into a copy of the real dictionary, rather than a lookup that would raise
    StopIteration against the shipped data. See task-7-report.md for the ruling.
    """
    benign_def = {"changes": [{"type": "insert", "start": 60, "sequence": "AAA"}]}
    assert net_length_change(benign_def) % 3 == 0
    rd = replace(RD, mutations={**RD.mutations, "synthetic_benign": benign_def})
    assert truth_class(_truth("synthetic_benign"), rd) == "benign"


def test_missing_summary_is_no_call(tmp_path: Path) -> None:
    assert predicted_decision(tmp_path) == "NO_CALL"
    (tmp_path / "summary.json").write_text("{not json")
    assert predicted_decision(tmp_path) == "NO_CALL"


def test_valid_json_wrong_shape_is_no_call(tmp_path: Path) -> None:
    """Valid JSON that is not a summary object still raises inside compute_clinical_decision."""
    (tmp_path / "summary.json").write_text("[1, 2, 3]")
    assert predicted_decision(tmp_path) == "NO_CALL"


def test_confusion_counts_critical_errors() -> None:
    rows = [
        {"clinical": {"truth": "pathogenic", "decision": "NO_CALL"}},
        {"clinical": {"truth": "pathogenic", "decision": "PATHOGENIC"}},
        {"clinical": {"truth": "normal", "decision": "PATHOGENIC"}},
        {"clinical": {"truth": "normal", "decision": "INCONCLUSIVE"}},
    ]
    c = confusion(rows)
    assert c["critical_false_negative"] == 1 and c["false_positive_normal"] == 1
    assert c["inconclusive_rate"] == 0.25
    assert c["matrix"]["pathogenic"]["NO_CALL"] == 1


def test_confusion_counts_false_positive_benign() -> None:
    rows = [{"clinical": {"truth": "benign", "decision": "PATHOGENIC"}}]
    assert confusion(rows)["false_positive_benign"] == 1


def test_confusion_empty_rows_has_null_inconclusive_rate() -> None:
    assert confusion([])["inconclusive_rate"] is None
