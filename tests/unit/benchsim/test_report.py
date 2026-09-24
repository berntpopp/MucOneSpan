"""benchsim.report: row normalization, stratified tables, paired tests, decision rule."""

import json
from pathlib import Path
from typing import Any

import pytest

from muc_one_span.benchsim.report import (
    RULE_TEXT,
    decide,
    normalize_rows,
    paired,
    preregister,
    render_markdown,
    require_preregistered,
    stratified_table,
)


def _rows(exact: list[int], fp: list[int], profile="ont_amplicon_r10"):
    return [
        {
            "sample": f"s{i}",
            "profile": profile,
            "allele_exact": e,
            "normal": True,
            "false_positive": f,
            "critical_false_negative": 0,
        }
        for i, (e, f) in enumerate(zip(exact, fp, strict=True))
    ]


def test_paired_counts_discordant() -> None:
    a, b = _rows([1, 0, 0, 1], [0] * 4), _rows([1, 1, 1, 0], [0] * 4)
    res = paired(a, b, "allele_exact")
    assert (res["b"], res["c"]) == (1, 2)


def test_decide_requires_superiority_and_noninferiority() -> None:
    n = 3000
    base = _rows([0] * n, [0] * n)
    cand = _rows([1] * n, [0] * n)
    assert decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid")["adopt"] is True
    worse = _rows([1] * n, [1] * 50 + [0] * (n - 50))
    assert decide({"ladder": base, "hybrid": worse}, "ladder", "hybrid")["adopt"] is False


def test_sealed_split_requires_preregistration(tmp_path: Path) -> None:
    path = tmp_path / "preregistration.jsonl"
    with pytest.raises(PermissionError):
        require_preregistered(path, "rule v1")
    preregister("rule v1", path)
    require_preregistered(path, "rule v1")
    with pytest.raises(PermissionError):
        require_preregistered(path, "rule v2")


# --- normalize_rows ------------------------------------------------------------------


def _metrics(all_exact: int, per_allele: int) -> dict[str, Any]:
    return {
        "all_sequences_exact": {"min": all_exact, "max": all_exact},
        "independent_sequence_exact": {"min": per_allele, "max": per_allele},
    }


def _sample(
    name: str,
    truth: str | None,
    decision: str,
    status: str = "completed",
    all_exact: int = 1,
    per_allele: int = 2,
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "sample": name,
        "status": status,
        "clinical": {"truth": truth, "decision": decision},
    }
    if truth is not None:
        row.update(
            truth_status="valid", truth_haplotypes=2, metrics=_metrics(all_exact, per_allele)
        )
    else:
        row["truth_status"] = "invalid"
    return row


def _case(name: str, event: str | None, profile: str = "hifi_amplicon") -> dict[str, Any]:
    return {
        "design_id": name,
        "profile": profile,
        "status": "ok",
        "realized_depth": {"1": 30, "2": 25},
        "design": {
            "design_id": name,
            "profile": profile,
            "delta_class": "2",
            "depth": 60,
            "composition": "markov",
            "event": event,
        },
    }


def test_normalize_rows_one_row_per_clinical_class() -> None:
    report = {
        "samples": [
            _sample("p", "pathogenic", "PATHOGENIC"),
            _sample("b", "benign", "PATHOGENIC", all_exact=0, per_allele=1),
            _sample("n", "normal", "NO_PATHOGENIC_VARIANT_DETECTED"),
            _sample("m", "pathogenic", "NO_PATHOGENIC_VARIANT_DETECTED", all_exact=0, per_allele=0),
        ]
    }
    cases = {
        "p": _case("p", "dupC"),
        "b": _case("b", "x"),
        "n": _case("n", None),
        "m": _case("m", "dupC"),
    }
    rows = {r["sample"]: r for r in normalize_rows(report, cases)}
    assert (
        rows["p"]["allele_exact"],
        rows["p"]["false_positive"],
        rows["p"]["critical_false_negative"],
        rows["p"]["normal"],
    ) == (1, 0, 0, False)
    assert rows["b"]["false_positive"] == 1 and rows["b"]["benign"] is True
    assert rows["b"]["allele_exact"] == 0 and rows["b"]["alleles_exact"] == 1
    assert rows["n"]["normal"] is True and rows["n"]["false_positive"] == 0
    assert rows["m"]["critical_false_negative"] == 1
    assert rows["p"]["profile"] == "hifi_amplicon" and rows["p"]["delta_class"] == "2"
    assert rows["p"]["depth"] == 60 and rows["p"]["realized_depth"] == 55
    assert all(r["failed"] is False for r in rows.values())


@pytest.mark.parametrize("status", ["execution_failed", "not_attempted"])
def test_normalize_rows_keeps_failures_as_failures(status: str) -> None:
    # A failed run with a stale PATHOGENIC summary still becomes NO_CALL, exact 0.
    report = {"samples": [_sample("f", "pathogenic", "PATHOGENIC", status=status)]}
    (row,) = normalize_rows(report, {"f": _case("f", "dupC")})
    assert row["decision"] == "NO_CALL" and row["failed"] is True
    assert row["allele_exact"] == 0 and row["alleles_exact"] == 0
    assert row["critical_false_negative"] == 1 and row["false_positive"] == 0


def test_normalize_rows_invalid_truth_uses_design_class() -> None:
    report = {"samples": [_sample("i", None, "NO_CALL", status="invalid_truth")]}
    (row,) = normalize_rows(report, {"i": _case("i", None)})
    assert row["truth"] == "normal" and row["normal"] is True and row["failed"] is True
    assert row["allele_exact"] == 0 and row["decision"] == "NO_CALL"


def test_normalize_rows_requires_every_case() -> None:
    report = {"samples": [_sample("p", "pathogenic", "PATHOGENIC")]}
    with pytest.raises(KeyError):
        normalize_rows(report, {})


# --- stratified_table / paired / decide ----------------------------------------------


def test_stratified_table_counts_and_intervals() -> None:
    rows = _rows([1, 0, 1, 1], [0] * 4) + _rows([0, 0], [0, 0], profile="hifi_amplicon")
    table = stratified_table(rows, ["profile"], "allele_exact")
    by = {t["stratum"]["profile"]: t for t in table}
    assert (by["ont_amplicon_r10"]["k"], by["ont_amplicon_r10"]["n"]) == (3, 4)
    assert by["ont_amplicon_r10"]["rate"] == 0.75
    assert 0 < by["ont_amplicon_r10"]["ci_low"] < 0.75 < by["ont_amplicon_r10"]["ci_high"]
    assert by["hifi_amplicon"]["k"] == 0 and by["hifi_amplicon"]["ci_low"] == 0.0


def test_paired_rejects_mismatched_samples() -> None:
    with pytest.raises(ValueError, match="sample"):
        paired(_rows([1, 0], [0, 0]), _rows([1], [0]), "allele_exact")


def test_decide_rejects_critical_false_negative_increase() -> None:
    n = 3000
    base = _rows([0] * n, [0] * n)
    cand = _rows([1] * n, [0] * n)
    cand[0]["critical_false_negative"] = 1
    result = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid")
    prof = result["profiles"]["ont_amplicon_r10"]
    assert prof["critical_false_negative"]["pass"] is False and result["adopt"] is False


def test_decide_every_profile_must_pass_and_needs_normals() -> None:
    n = 3000
    base = _rows([0] * n, [0] * n) + _rows([0] * 10, [0] * 10, profile="hifi_amplicon")
    cand = _rows([1] * n, [0] * n) + _rows([0] * 10, [0] * 10, profile="hifi_amplicon")
    for i, row in enumerate(base[n:] + cand[n:]):
        row["sample"] = f"h{i % 10}"
    result = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid")
    assert result["profiles"]["ont_amplicon_r10"]["pass"] is True
    assert result["profiles"]["hifi_amplicon"]["pass"] is False and result["adopt"] is False
    for row in base + cand:
        row["normal"] = False
    no_normals = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid")
    fp = no_normals["profiles"]["ont_amplicon_r10"]["false_positive"]
    assert fp["noninferior"] is False and fp["n"] == 0 and no_normals["adopt"] is False


def test_decide_holm_adjusts_primary_family() -> None:
    n = 3000
    result = decide(
        {"ladder": _rows([0] * n, [0] * n), "hybrid": _rows([1] * n, [0] * n)}, "ladder", "hybrid"
    )
    exact = result["profiles"]["ont_amplicon_r10"]["allele_exact"]
    assert exact["p_holm"] >= exact["p"] and set(result["holm_family"]) == {
        "allele_exact",
        "false_positive",
        "critical_false_negative",
    }


def test_decide_empty_is_not_adopted() -> None:
    assert decide({"ladder": [], "hybrid": []}, "ladder", "hybrid")["adopt"] is False


# --- pre-registration and rendering --------------------------------------------------


def test_preregister_is_append_only_and_hashes_rule(tmp_path: Path) -> None:
    path = tmp_path / "test" / "preregistration.jsonl"
    first = preregister(RULE_TEXT, path)
    preregister("other rule", path)
    preregister(RULE_TEXT, path)
    entries = [json.loads(line) for line in path.read_text().splitlines()]
    assert len(entries) == 3 and entries[0]["sha256"] == first == entries[2]["sha256"]
    assert entries[0]["rule_text"] == RULE_TEXT and "registered_at" in entries[0]
    require_preregistered(path, RULE_TEXT)


def test_require_preregistered_rejects_corrupt_ledger(tmp_path: Path) -> None:
    path = tmp_path / "preregistration.jsonl"
    path.write_text("not json\n")
    with pytest.raises(PermissionError):
        require_preregistered(path, RULE_TEXT)


def test_render_markdown_lists_profiles_and_verdict() -> None:
    n = 3000
    result = decide(
        {"ladder": _rows([0] * n, [0] * n), "hybrid": _rows([1] * n, [0] * n)}, "ladder", "hybrid"
    )
    result["tables"] = {
        "allele_exact by profile": stratified_table(
            _rows([1, 0], [0, 0]), ["profile"], "allele_exact"
        )
    }
    text = render_markdown(result)
    assert "ont_amplicon_r10" in text and "ADOPT" in text and "| profile |" in text
    assert "allele_exact by profile" in text
    assert "NOT ADOPTED" in render_markdown({"profiles": {}, "adopt": False})
