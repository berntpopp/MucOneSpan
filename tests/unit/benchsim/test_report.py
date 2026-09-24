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


def test_decide_requires_superiority_and_noninferiority() -> None:
    n = 3000
    base = _rows([0] * n, [0] * n)
    cand = _rows([1] * n, [0] * n)
    assert decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", NO_TARGETS)["adopt"] is True
    worse = _rows([1] * n, [1] * 50 + [0] * (n - 50))
    result = decide({"ladder": base, "hybrid": worse}, "ladder", "hybrid", NO_TARGETS)
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


def test_decide_every_profile_must_pass_and_needs_normals() -> None:
    n = 3000
    base = _rows([0] * n, [0] * n) + _rows([0] * 10, [0] * 10, profile="hifi_amplicon")
    cand = _rows([1] * n, [0] * n) + _rows([0] * 10, [0] * 10, profile="hifi_amplicon")
    for i, row in enumerate(base[n:] + cand[n:]):
        row["sample"] = f"h{i % 10}"
    result = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", NO_TARGETS)
    assert result["profiles"]["ont_amplicon_r10"]["pass"] is True
    assert result["profiles"]["hifi_amplicon"]["pass"] is False and result["adopt"] is False
    for row in base + cand:
        row["normal"] = False
    no_normals = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", NO_TARGETS)
    fp = no_normals["profiles"]["ont_amplicon_r10"]["false_positive"]
    assert fp["noninferior"] is False and fp["n"] == 0 and no_normals["adopt"] is False


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
    margin = DEFAULT_BENCH_CONFIG.report.ni_margin
    assert "cluster-bootstrap" in RULE_TEXT and f"{margin:g}" in RULE_TEXT


def test_default_rule_text_is_unchanged_by_the_config_refactor() -> None:
    # SHA-256 of the v4 rule text (v3 plus the task 12e absolute targets, with an
    # absent targeted set reported as not present); a changed default would
    # silently invalidate existing pre-registrations.
    pinned = "bef9f891331f8fac12f0b8db0c2feecfd0d78042c0ec6e888d9de0f3243bf7e4"
    assert rule_text(DEFAULT_BENCH_CONFIG.report) == RULE_TEXT
    assert rule_sha256(RULE_TEXT) == pinned


def test_rule_text_and_decision_follow_the_report_config() -> None:
    report = replace(DEFAULT_BENCH_CONFIG.report, ni_margin=0.25, alpha=0.1)
    empty = NO_TARGETS.targets
    text = rule_text(report, targets=empty)
    assert "below 0.25" in text and "alpha 0.1" in text and "one-sided 90%" in text
    n = 30
    base, cand = _rows([0] * n, [0] * n), _rows([1] * n, [0] * n)
    result = decide(
        {"ladder": base, "hybrid": cand},
        "ladder",
        "hybrid",
        BenchConfig(report=report, targets=empty),
    )
    assert result["margin"] == 0.25 and result["alpha"] == 0.1
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


# --- absolute targets (task 12e) ------------------------------------------------------


def _target_row(
    sample: str, profile: str, bench_set: str, exact: int, truth: str, decision: str
) -> dict[str, Any]:
    """A normalized case row exercising both the relative rule and the absolute targets."""
    fp = int(decision == "PATHOGENIC" and truth in ("normal", "benign"))
    false_negative = int(truth == "pathogenic" and decision != "PATHOGENIC")
    return {
        "sample": sample,
        "profile": profile,
        "bench_set": bench_set,
        "alleles": [{"allele": "h1", "allele_exact": exact}],
        "normal": truth == "normal",
        "benign": truth == "benign",
        "pathogenic": truth == "pathogenic",
        "decision": decision,
        "inconclusive": int(decision == "INCONCLUSIVE"),
        "false_positive": fp,
        "critical_false_negative": false_negative,
    }


def _standard_fixture(
    n_path: int, n_normal: int
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """(baseline, candidate) ``standard``-set rows: candidate superior (exact=1, fewer
    critical FN) and meeting the default ``standard`` targets (pathogenic_rate >= 0.80,
    inconclusive_rate <= 0.20, false_positive_rate <= 0)."""
    base, cand = [], []
    for i in range(n_path):
        sample = f"p{i}"
        base.append(_target_row(sample, "ont_amplicon_r10", "standard", 0, "pathogenic", "NO_CALL"))
        cand.append(
            _target_row(sample, "ont_amplicon_r10", "standard", 1, "pathogenic", "PATHOGENIC")
        )
    for i in range(n_normal):
        sample, decision = f"n{i}", "NO_PATHOGENIC_VARIANT_DETECTED"
        base.append(_target_row(sample, "ont_amplicon_r10", "standard", 0, "normal", decision))
        cand.append(_target_row(sample, "ont_amplicon_r10", "standard", 1, "normal", decision))
    return base, cand


def _clean_fixture(n_path: int, n_normal: int) -> list[dict[str, Any]]:
    """Candidate-only ``clean``-set rows meeting its targets (>= 0.90 / <= 0.10 / <= 0)."""
    cand = [
        _target_row(f"cp{i}", "ont_amplicon_r10", "clean", 1, "pathogenic", "PATHOGENIC")
        for i in range(n_path)
    ]
    cand += [
        _target_row(
            f"cn{i}", "ont_amplicon_r10", "clean", 1, "normal", "NO_PATHOGENIC_VARIANT_DETECTED"
        )
        for i in range(n_normal)
    ]
    return cand


def test_decide_adopts_only_when_headline_and_clean_targets_also_pass() -> None:
    base, cand_standard = _standard_fixture(2700, 1000)
    clean_ok = _clean_fixture(9, 1)
    result = decide({"ladder": base, "hybrid": cand_standard + clean_ok}, "ladder", "hybrid")
    assert set(result["targets"]) == {"standard", "clean"}
    assert result["targets"]["standard"]["pass"] is True
    assert result["targets"]["clean"]["pass"] is True
    assert result["profiles"]["ont_amplicon_r10"]["pass"] is True  # the relative rule too
    assert result["adopt"] is True

    # Targets are not restricted to the headline set: a `clean`-only regression (one
    # false positive) fails adoption even though the relative rule and the `standard`
    # targets are untouched.
    clean_bad = _clean_fixture(9, 1)
    clean_bad[-1] = clean_bad[-1] | {"decision": "PATHOGENIC", "false_positive": 1}
    result2 = decide({"ladder": base, "hybrid": cand_standard + clean_bad}, "ladder", "hybrid")
    assert result2["targets"]["standard"]["pass"] is True
    assert result2["targets"]["clean"]["pass"] is False
    assert result2["profiles"]["ont_amplicon_r10"]["pass"] is True
    assert result2["adopt"] is False

    # Review I4: a targeted set absent from this root is "not present", never FAIL,
    # and still blocks adoption (fail closed).
    result3 = decide({"ladder": base, "hybrid": cand_standard}, "ladder", "hybrid")
    assert result3["targets"]["clean"]["present"] is False
    assert result3["targets"]["clean"]["pass"] is None
    assert result3["targets_not_present"] == ["clean"]
    assert result3["profiles"]["ont_amplicon_r10"]["pass"] is True
    assert result3["adopt"] is False
    assert "not present" in render_markdown(result3)


def test_decide_targets_use_point_estimate_by_default_and_ci_bound_when_configured() -> None:
    base, cand = _standard_fixture(4, 1)  # 4/4 pathogenic = 1.0, point estimate passes easily
    cand = cand + _clean_fixture(9, 1)
    point = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", NO_TARGETS)
    assert point["targets"] == {}  # NO_TARGETS configures no sets at all

    with_targets = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid")
    assert with_targets["targets"]["standard"]["basis"] == "point"
    ci_cfg = replace(
        DEFAULT_BENCH_CONFIG, targets=replace(DEFAULT_BENCH_CONFIG.targets, basis="ci_bound")
    )
    with_ci = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid", ci_cfg)
    assert with_ci["targets"]["standard"]["basis"] == "ci_bound"
    pooled_path = next(
        r
        for r in with_ci["targets"]["standard"]["table"]
        if r["grouping"] == "pooled" and r["metric"] == "pathogenic_rate"
    )
    assert pooled_path["judged"] == pooled_path["ci_low"]


def test_render_markdown_shows_the_absolute_targets_table() -> None:
    base, cand = _standard_fixture(2700, 1000)
    cand = cand + _clean_fixture(9, 1)
    result = decide({"ladder": base, "hybrid": cand}, "ladder", "hybrid")
    text = render_markdown(result)
    assert "Part 2: absolute targets" in text and "Absolute targets (task 12e)" in text
    assert "| standard | pooled | pathogenic_rate |" in text
    assert "Set verdict: standard=True, clean=True" in text
    assert "**ADOPT**" in text
