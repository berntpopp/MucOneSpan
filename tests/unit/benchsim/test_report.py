"""benchsim.report: row normalization, stratified tables, paired tests, decision rule."""

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

from muc_one_span.benchsim.bench_config import DEFAULT_BENCH_CONFIG, BenchConfig, TargetsConfig
from muc_one_span.benchsim.report import (
    ALLELE_UNIT,
    RULE_TEXT,
    allele_rows,
    decide,
    mark_first_evaluation,
    normalize_rows,
    paired,
    preregister,
    render_markdown,
    require_preregistered,
    rule_sha256,
    rule_text,
    stratified_table,
)


def _rows(
    exact: list[int], fp: list[int], profile: str = "ont_amplicon_r10"
) -> list[dict[str, Any]]:
    """Case rows with one truth allele each (``h1``), exact per ``exact``."""
    return [
        {
            "sample": f"s{i}",
            "profile": profile,
            "bench_set": HEADLINE,
            "alleles": [{"allele": "h1", "allele_exact": e}],
            "normal": True,
            "false_positive": f,
            "critical_false_negative": 0,
        }
        for i, (e, f) in enumerate(zip(exact, fp, strict=True))
    ]


HEADLINE = DEFAULT_BENCH_CONFIG.sets.headline
# Bypasses the task 12e absolute targets (part 2) so pre-12e fixtures that exercise only
# the relative rule (part 1) do not also need pathogenic/decision/inconclusive fields.
NO_TARGETS = BenchConfig(targets=TargetsConfig(by_set={}))


def test_paired_counts_discordant() -> None:
    a, b = _rows([1, 0, 0, 1], [0] * 4), _rows([1, 1, 1, 0], [0] * 4)
    res = paired(allele_rows(a), allele_rows(b), "allele_exact", ALLELE_UNIT)
    assert (res["b"], res["c"]) == (1, 2) and res["unit"] == ["sample", "allele"]


def test_decide_requires_superiority() -> None:
    n = 3000
    base = _rows([0] * n, [0] * n)
    cand = _rows([1] * n, [0] * n)
    assert decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", NO_TARGETS)["adopt"] is True
    not_superior = _rows([0] * n, [0] * n)
    result = decide({"ladder": base, "hybrid": not_superior}, "ladder", "hybrid", NO_TARGETS)
    assert result["adopt"] is False


def test_sealed_split_requires_preregistration(tmp_path: Path) -> None:
    path = tmp_path / "preregistration.jsonl"
    with pytest.raises(PermissionError):
        require_preregistered(path, "rule v1")
    preregister("rule v1", path)
    require_preregistered(path, "rule v1")
    with pytest.raises(PermissionError):
        require_preregistered(path, "rule v2")


# --- normalize_rows ------------------------------------------------------------------


def _pair(truth: str, exact: bool) -> dict[str, Any]:
    return {"truth": truth, "prediction": f"p_{truth}", "sequence_exact": exact}


def _alt(pairs: list[dict[str, Any]], missing: list[str], independent: int) -> dict[str, Any]:
    exact = sum(p["sequence_exact"] for p in pairs)
    return {
        "pairs": pairs,
        "missing_truth": missing,
        "metrics": {"independent_sequence_exact": independent, "sequence_exact": exact},
    }


def _sample(
    name: str,
    truth: str | None,
    decision: str,
    status: str = "completed",
    exact: tuple[bool, bool] = (True, True),
    alternatives: list[dict[str, Any]] | None = None,
    tp: int = 1,
    fp: int = 0,
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "sample": name,
        "status": status,
        "clinical": {"truth": truth, "decision": decision},
    }
    if truth is None:
        return row | {"truth_status": "invalid"}
    alts = alternatives or [_alt([_pair("h1", exact[0]), _pair("h2", exact[1])], [], sum(exact))]
    all_exact = int(all(exact))
    return row | {
        "truth_status": "valid",
        "truth_haplotypes": 2,
        "truth_events": 0 if truth == "normal" else 1,
        "alternatives": alts,
        "metrics": {
            "all_sequences_exact": {"min": all_exact, "max": all_exact},
            "independent_sequence_exact": {"min": sum(exact), "max": sum(exact)},
            "event_tp": {"min": tp, "max": tp},
            "event_fp": {"min": 0, "max": fp},
        },
    }


def _case(name: str, event: str | None, profile: str = "hifi_amplicon") -> dict[str, Any]:
    design = {
        "design_id": name,
        "profile": profile,
        "delta_class": "2",
        "depth": 60,
        "composition": "markov",
        "event": event,
    }
    return {
        "design_id": name,
        "profile": profile,
        "status": "ok",
        "realized_depth": {"1": 30, "2": 25},
        "design": design,
    }


def test_normalize_rows_one_row_per_clinical_class() -> None:
    report = {
        "samples": [
            _sample("p", "pathogenic", "PATHOGENIC"),
            _sample("b", "benign", "PATHOGENIC", exact=(False, True)),
            _sample("n", "normal", "NO_PATHOGENIC_VARIANT_DETECTED", tp=0, fp=1),
            _sample("m", "pathogenic", "NO_PATHOGENIC_VARIANT_DETECTED", exact=(False, False)),
        ]
    }
    cases = {
        "p": _case("p", "dupC"),
        "b": _case("b", "x"),
        "n": _case("n", None),
        "m": _case("m", "dupC"),
    }
    rows = {r["sample"]: r for r in normalize_rows(report, cases)}
    p = rows["p"]
    assert (p["case_exact"], p["false_positive"], p["critical_false_negative"]) == (1, 0, 0)
    assert p["alleles"] == [
        {"allele": "h1", "allele_exact": 1},
        {"allele": "h2", "allele_exact": 1},
    ]
    assert rows["b"]["false_positive"] == 1 and rows["b"]["benign"] is True
    assert rows["b"]["case_exact"] == 0 and rows["b"]["alleles_exact"] == 1
    assert [a["allele_exact"] for a in rows["b"]["alleles"]] == [0, 1]
    assert rows["n"]["normal"] is True and rows["n"]["false_positive"] == 0
    assert (rows["n"]["event_tp"], rows["n"]["event_fp"], rows["n"]["truth_events"]) == (0, 1, 0)
    assert rows["m"]["critical_false_negative"] == 1 and rows["m"]["failure"] == 1
    assert p["profile"] == "hifi_amplicon" and p["delta_class"] == "2"
    assert p["depth"] == 60 and p["realized_depth"] == 55 and p["failure"] == 0
    assert all(r["failed"] is False for r in rows.values())
    flat = allele_rows(list(rows.values()))
    assert len(flat) == 8 and {r["sample"] for r in flat} == {"p", "b", "n", "m"}
    assert flat[0]["profile"] == "hifi_amplicon" and "allele" in flat[0]


@pytest.mark.parametrize("status", ["execution_failed", "not_attempted"])
def test_normalize_rows_keeps_failures_as_failures(status: str) -> None:
    # A failed run with a stale PATHOGENIC summary still becomes NO_CALL, every allele 0.
    alts = [_alt([], ["h1", "h2"], 0)]
    report = {"samples": [_sample("f", "pathogenic", "PATHOGENIC", status, alternatives=alts)]}
    (row,) = normalize_rows(report, {"f": _case("f", "dupC")})
    assert row["decision"] == "NO_CALL" and row["failed"] is True
    assert row["alleles"] == [
        {"allele": "h1", "allele_exact": 0},
        {"allele": "h2", "allele_exact": 0},
    ]
    assert row["case_exact"] == 0 and row["event_tp"] == 0
    assert row["critical_false_negative"] == 1 and row["false_positive"] == 0


def test_normalize_rows_uses_least_favourable_assignment_and_evidence() -> None:
    # Two optimal assignments; the worse one has one exact allele. Of two literal
    # exact pairs only one has independent evidence -> the later truth name is demoted.
    alts = [
        _alt([_pair("h1", True), _pair("h2", True)], [], 2),
        _alt([_pair("h1", True), _pair("h2", True)], [], 1),
    ]
    report = {
        "samples": [_sample("d", "normal", "NO_PATHOGENIC_VARIANT_DETECTED", alternatives=alts)]
    }
    (row,) = normalize_rows(report, {"d": _case("d", None)})
    assert row["alleles"] == [
        {"allele": "h1", "allele_exact": 1},
        {"allele": "h2", "allele_exact": 0},
    ]


def test_normalize_rows_invalid_truth_raises() -> None:
    report = {"samples": [_sample("i", None, "NO_CALL", status="invalid_truth")]}
    with pytest.raises(ValueError, match="no valid truth"):
        normalize_rows(report, {"i": _case("i", None)})


def test_normalize_rows_allele_count_must_match_truth() -> None:
    alts = [_alt([_pair("h1", True)], [], 1)]
    report = {"samples": [_sample("x", "normal", "NO_CALL", alternatives=alts)]}
    with pytest.raises(ValueError, match="truth allele names"):
        normalize_rows(report, {"x": _case("x", None)})


def test_normalize_rows_requires_every_case() -> None:
    report = {"samples": [_sample("p", "pathogenic", "PATHOGENIC")]}
    with pytest.raises(KeyError):
        normalize_rows(report, {})


# --- stratified_table / paired / decide ----------------------------------------------


def test_stratified_table_counts_and_intervals() -> None:
    rows = allele_rows(_rows([1, 0, 1, 1], [0] * 4) + _rows([0, 0], [0, 0], "hifi_amplicon"))
    table = stratified_table(rows, ["profile"], "allele_exact")
    by = {t["stratum"]["profile"]: t for t in table}
    assert (by["ont_amplicon_r10"]["k"], by["ont_amplicon_r10"]["n"]) == (3, 4)
    assert by["ont_amplicon_r10"]["rate"] == 0.75
    assert 0 < by["ont_amplicon_r10"]["ci_low"] < 0.75 < by["ont_amplicon_r10"]["ci_high"]
    assert by["hifi_amplicon"]["k"] == 0 and by["hifi_amplicon"]["ci_low"] == 0.0


def test_paired_rejects_mismatched_units() -> None:
    a = allele_rows(_rows([1, 0], [0, 0]))
    b = allele_rows(_rows([1, 0], [0, 0]))
    b[1]["allele"] = "h2"
    with pytest.raises(ValueError, match="sample/allele"):
        paired(a, b, "allele_exact", ALLELE_UNIT)
    with pytest.raises(ValueError, match="sample"):
        paired(_rows([1, 0], [0, 0]), _rows([1], [0]), "false_positive")


def test_decide_pairs_alleles_not_cases() -> None:
    n = 40
    base, cand = _rows([0] * n, [0] * n), _rows([0] * n, [0] * n)
    for row in base + cand:
        row["alleles"] = [{"allele": "h1", "allele_exact": 0}, {"allele": "h2", "allele_exact": 0}]
    for row in cand:
        row["alleles"][1]["allele_exact"] = 1
    exact = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", NO_TARGETS)["profiles"][
        "ont_amplicon_r10"
    ]["allele_exact"]
    assert (exact["n"], exact["b"], exact["c"]) == (2 * n, 0, n) and exact["superior"] is True


def test_decide_rejects_critical_false_negative_increase() -> None:
    n = 3000
    base = _rows([0] * n, [0] * n)
    cand = _rows([1] * n, [0] * n)
    cand[0]["critical_false_negative"] = 1
    result = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", NO_TARGETS)
    prof = result["profiles"]["ont_amplicon_r10"]
    assert prof["critical_false_negative"]["pass"] is False and result["adopt"] is False


def test_decide_every_profile_must_pass() -> None:
    n = 3000
    base = _rows([0] * n, [0] * n) + _rows([0] * 10, [0] * 10, profile="hifi_amplicon")
    cand = _rows([1] * n, [0] * n) + _rows([0] * 10, [0] * 10, profile="hifi_amplicon")
    for i, row in enumerate(base[n:] + cand[n:]):
        row["sample"] = f"h{i % 10}"
    result = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", NO_TARGETS)
    assert result["profiles"]["ont_amplicon_r10"]["pass"] is True
    assert result["profiles"]["hifi_amplicon"]["pass"] is False and result["adopt"] is False


def test_decide_ignores_false_positives_with_no_normal_or_benign_rows() -> None:
    # Task C1: FP is no longer part of the relative rule. An empty FP cohort
    # (no normal/benign rows) is reported with a null rate/CI, not forced to fail.
    n = 3000
    base = _rows([0] * n, [0] * n)
    cand = _rows([1] * n, [0] * n)
    for row in base + cand:
        row["normal"] = False
    result = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", NO_TARGETS)
    fp = result["profiles"]["ont_amplicon_r10"]["false_positive"]
    assert fp["n"] == 0 and fp["rate"] is None and fp["ci_low"] is None and fp["ci_high"] is None
    assert result["profiles"]["ont_amplicon_r10"]["pass"] is True and result["adopt"] is True


def test_decide_fp_no_longer_gated_by_a_relative_margin_at_planned_test_size() -> None:
    # The exact final-review C1 scenario: 280 normals/profile, 0 FP in both engines.
    # The Newcombe non-inferiority bound could never clear a meaningful margin here
    # (upper bound 0.00957 > a 0.005 margin); the relative rule no longer tests it.
    n_normal = 280
    base = _rows([0] * n_normal, [0] * n_normal)
    cand = _rows([1] * n_normal, [0] * n_normal)
    result = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", NO_TARGETS)
    fp = result["profiles"]["ont_amplicon_r10"]["false_positive"]
    assert fp["fp"] == 0 and fp["fp_baseline"] == 0 and fp["n"] == n_normal
    assert fp["rate"] == 0.0 and "noninferior" not in fp and "margin" not in fp
    assert result["profiles"]["ont_amplicon_r10"]["pass"] is True
    assert result["adopt"] is True


def test_decide_holm_adjusts_primary_family() -> None:
    n = 3000
    result = decide(
        {"ladder": _rows([0] * n, [0] * n), "hybrid": _rows([1] * n, [0] * n)},
        "ladder",
        "hybrid",
        NO_TARGETS,
    )
    exact = result["profiles"]["ont_amplicon_r10"]["allele_exact"]
    assert exact["p_holm"] >= exact["p"]
    assert set(result["holm_family"]) == {
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
    matched = require_preregistered(path, RULE_TEXT)
    assert matched == {"sha256": first, "registered_at": entries[0]["registered_at"]}


def test_changing_a_target_requires_a_new_pre_registration(tmp_path: Path) -> None:
    path = tmp_path / "test" / "preregistration.jsonl"
    preregister(RULE_TEXT, path)
    require_preregistered(path, RULE_TEXT)  # the registered (default-targets) rule seals test
    changed_threshold = DEFAULT_BENCH_CONFIG.targets.by_set["clean"]["pathogenic_rate"]
    changed = replace(
        DEFAULT_BENCH_CONFIG.targets,
        by_set=DEFAULT_BENCH_CONFIG.targets.by_set
        | {
            "clean": DEFAULT_BENCH_CONFIG.targets.by_set["clean"]
            | {"pathogenic_rate": replace(changed_threshold, threshold=0.95)}
        },
    )
    text = rule_text(targets=changed)
    assert rule_sha256(text) != rule_sha256(RULE_TEXT)
    with pytest.raises(PermissionError, match="not pre-registered"):
        require_preregistered(path, text)
    preregister(text, path)
    require_preregistered(path, text)  # now registered too

    # A basis change (point vs. ci_bound) also changes the sha256.
    ci_text = rule_text(targets=replace(DEFAULT_BENCH_CONFIG.targets, basis="ci_bound"))
    assert rule_sha256(ci_text) != rule_sha256(RULE_TEXT) != rule_sha256(ci_text)


def test_preregister_refused_after_test_evaluated(tmp_path: Path) -> None:
    path = tmp_path / "test" / "preregistration.jsonl"
    preregister(RULE_TEXT, path)
    stamp = mark_first_evaluation(path, rule_sha256(RULE_TEXT))
    assert mark_first_evaluation(path, rule_sha256(RULE_TEXT)) == stamp  # written once
    with pytest.raises(PermissionError, match="already evaluated"):
        preregister("a later rule", path)


def test_rule_text_names_the_per_allele_endpoint() -> None:
    assert "per-allele" in RULE_TEXT and "allele pairs" in RULE_TEXT
    assert "cluster-bootstrap" in RULE_TEXT
    assert "Clopper-Pearson" in RULE_TEXT and "information only" in RULE_TEXT


def test_default_rule_text_is_unchanged_by_the_config_refactor() -> None:
    # SHA-256 of the v5 rule text (v4 plus the task C1 owner ruling dropping the
    # relative FP non-inferiority margin; FP is now judged only by the absolute
    # targets and reported informationally). A changed default would silently
    # invalidate existing pre-registrations. No real pre-registration exists yet
    # (task C1), so this is simply re-pinned to the current text.
    pinned = "9249bf485f3c7b11b5415e07f2165bad743df17b476f42a0ece98a9bcbfa96f0"
    assert rule_text(DEFAULT_BENCH_CONFIG.report) == RULE_TEXT
    assert rule_sha256(RULE_TEXT) == pinned


def test_rule_text_and_decision_follow_the_report_config() -> None:
    report = replace(DEFAULT_BENCH_CONFIG.report, alpha=0.1)
    empty = NO_TARGETS.targets
    text = rule_text(report, targets=empty)
    assert "alpha 0.1" in text and "90%" in text
    n = 30
    base, cand = _rows([0] * n, [0] * n), _rows([1] * n, [0] * n)
    result = decide(
        {"ladder": base, "hybrid": cand},
        "ladder",
        "hybrid",
        BenchConfig(report=report, targets=empty),
    )
    assert result["alpha"] == 0.1
    assert result["rule_sha256"] == rule_sha256(text)


def test_require_preregistered_rejects_corrupt_ledger(tmp_path: Path) -> None:
    path = tmp_path / "preregistration.jsonl"
    path.write_text("not json\n")
    with pytest.raises(PermissionError):
        require_preregistered(path, RULE_TEXT)


def test_render_markdown_lists_profiles_and_verdict() -> None:
    n = 3000
    result = decide(
        {"ladder": _rows([0] * n, [0] * n), "hybrid": _rows([1] * n, [0] * n)},
        "ladder",
        "hybrid",
        NO_TARGETS,
    )
    result["tables"] = {
        "allele_exact by profile": stratified_table(
            allele_rows(_rows([1, 0], [0, 0])), ["profile"], "allele_exact"
        )
    }
    text = render_markdown(result)
    assert "ont_amplicon_r10" in text and "ADOPT" in text and "| profile |" in text
    assert "allele_exact by profile" in text and "normal + benign" in text and "no-call" in text
    assert "NOT ADOPTED" in render_markdown({"candidate": "hybrid", "baseline": "ladder"})


def test_render_markdown_without_candidate_is_not_a_verdict() -> None:
    text = render_markdown({"profiles": {}, "adopt": False})
    assert "not evaluated" in text and "ADOPT" not in text


def test_normalize_rows_keeps_reasons_position_and_minimum_depth() -> None:
    sample = _sample("p", "pathogenic", "INCONCLUSIVE")
    sample["clinical"]["reasons"] = ["Allele 1: Genotype phase is unphased or conflicting."]
    sample["reconstruction_flags"] = ["ambiguous_reconstruction", "iupac_bases"]
    case = _case("p", "dupC")
    case["design"]["event_position"] = "first10"
    (row,) = normalize_rows({"samples": [sample]}, {"p": case})
    assert row["clinical_reasons"] == sample["clinical"]["reasons"]
    assert row["reconstruction_flags"] == ["ambiguous_reconstruction", "iupac_bases"]
    assert row["reasons_recorded"] is True and row["event_position"] == "first10"
    assert row["realized_min_allele_depth"] == min(case["realized_depth"].values())
    legacy = _sample("q", "normal", "INCONCLUSIVE")
    (old,) = normalize_rows(
        {"samples": [legacy]}, {"q": _case("q", None) | {"realized_depth": None}}
    )
    assert old["reasons_recorded"] is False and old["clinical_reasons"] == []
    assert old["realized_min_allele_depth"] is None


def test_rule_applies_to_the_headline_set_only() -> None:
    assert f"`{HEADLINE}` benchmark set only" in RULE_TEXT
    other = rule_text(DEFAULT_BENCH_CONFIG.report, "clean")
    assert "`clean` benchmark set only" in other and rule_sha256(other) != rule_sha256(RULE_TEXT)


def test_decide_uses_only_headline_rows() -> None:
    n = 3000
    base, cand = _rows([0] * n, [0] * n), _rows([1] * n, [0] * n)
    stress = [r | {"sample": f"x{i}", "bench_set": "stress"} for i, r in enumerate(base)]
    worse = [r | {"false_positive": 1} for r in stress]
    result = decide(
        {"ladder": base + stress, "hybrid": cand + worse}, "ladder", "hybrid", NO_TARGETS
    )
    assert result["adopt"] is True and result["bench_set"] == HEADLINE
    assert result["profiles"]["ont_amplicon_r10"]["n"] == n


def test_normalize_rows_records_the_set() -> None:
    report = {"samples": [_sample("p", "pathogenic", "PATHOGENIC")]}
    case = _case("p", "dupC")
    (legacy,) = normalize_rows(report, {"p": case})
    assert legacy["bench_set"] == DEFAULT_BENCH_CONFIG.sets.legacy
    (named,) = normalize_rows(
        report, {"p": case | {"design": case["design"] | {"bench_set": "clean"}}}
    )
    assert named["bench_set"] == "clean"
    (other,) = normalize_rows(report, {"p": case}, legacy_set="standard")
    assert other["bench_set"] == "standard"
