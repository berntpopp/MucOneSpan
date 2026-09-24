"""benchsim.report_tables: metric 1-5 tables, cluster-bootstrap pooled CIs, failure atlas."""

from typing import Any

from muc_one_span.benchsim.report import STRATA
from muc_one_span.benchsim.report_tables import (
    build_tables,
    confusion_by_profile,
    event_table,
    pooled_estimates,
    render_engine_tables,
)


def _row(i: int, profile: str, truth: str, decision: str, exact: tuple[int, int]) -> dict[str, Any]:
    return {
        "sample": f"{profile}-{i}",
        "profile": profile,
        "delta_class": "1",
        "depth": 30,
        "composition": "markov",
        "pcr": "calibrated",
        "smear": 0.05,
        "chimera": 0.01,
        "error": "calibrated",
        "event": "dupC" if truth == "pathogenic" else None,
        "truth": truth,
        "decision": decision,
        "normal": truth == "normal",
        "benign": truth == "benign",
        "alleles": [
            {"allele": "h1", "allele_exact": exact[0]},
            {"allele": "h2", "allele_exact": exact[1]},
        ],
        "case_exact": int(all(exact)),
        "truth_events": int(truth != "normal"),
        "event_tp": int(truth != "normal" and decision == "PATHOGENIC"),
        "event_fp": int(truth == "normal" and decision == "PATHOGENIC"),
        "false_positive": int(truth != "pathogenic" and decision == "PATHOGENIC"),
        "critical_false_negative": int(truth == "pathogenic" and decision != "PATHOGENIC"),
        "failure": int(not all(exact)),
        "inconclusive": int(decision == "INCONCLUSIVE"),
        "no_call": int(decision == "NO_CALL"),
    }


ROWS = [
    _row(0, "ont_amplicon_r10", "pathogenic", "PATHOGENIC", (1, 1)),
    _row(1, "ont_amplicon_r10", "pathogenic", "NO_CALL", (0, 0)),
    _row(2, "ont_amplicon_r10", "normal", "PATHOGENIC", (1, 0)),
    _row(3, "hifi_amplicon", "normal", "NO_PATHOGENIC_VARIANT_DETECTED", (1, 1)),
]


def test_pooled_estimates_cluster_bootstrap_per_profile() -> None:
    est = pooled_estimates(ROWS, n_boot=200, seed=1)
    ont = est["ont_amplicon_r10"]
    assert ont["allele_exact"]["n"] == 6 and ont["allele_exact"]["clusters"] == 3
    assert ont["allele_exact"]["point"] == 0.5 and ont["case_exact"]["point"] == 1 / 3
    assert 0 <= ont["allele_exact"]["ci_low"] <= 0.5 <= ont["allele_exact"]["ci_high"] <= 1
    assert est["all"]["allele_exact"]["n"] == 8
    assert pooled_estimates([], n_boot=10)["all"]["allele_exact"]["point"] is None


def test_event_table_separates_dupc() -> None:
    table = {(r["profile"], r["event_class"]): r for r in event_table(ROWS)}
    dupc = table[("ont_amplicon_r10", "dupC")]
    assert (dupc["recall"]["k"], dupc["recall"]["n"]) == (1, 2)
    other = table[("ont_amplicon_r10", "other")]
    assert other["recall"]["rate"] is None and other["precision"]["k"] == 0
    assert other["precision"]["n"] == 1


def test_confusion_by_profile_reuses_clinical_confusion() -> None:
    conf = confusion_by_profile(ROWS)
    assert conf["ont_amplicon_r10"]["critical_false_negative"] == 1
    assert conf["ont_amplicon_r10"]["false_positive_normal"] == 1
    assert conf["all"]["matrix"]["normal"] == {"PATHOGENIC": 1, "NO_PATHOGENIC_VARIANT_DETECTED": 1}


def test_build_tables_has_all_metrics_and_failure_atlas() -> None:
    tables = build_tables(ROWS, n_boot=50)
    names = tables["stratified"]
    assert "allele_exact (metric 1, per allele) by profile" in names
    assert "case_exact (metric 2) by profile x delta_class" in names
    assert "false_positive on normal + benign truths by profile" in names
    assert "no_call on normal + benign truths by profile" in names
    for factor in STRATA:
        assert f"failure atlas: {factor}" in names
    allele = {
        t["stratum"]["profile"]: t for t in names["allele_exact (metric 1, per allele) by profile"]
    }
    assert (allele["ont_amplicon_r10"]["k"], allele["ont_amplicon_r10"]["n"]) == (3, 6)
    text = render_engine_tables(tables)
    assert "cluster bootstrap" in text and "metric 3" in text and "metric 4" in text
    assert "failure atlas: pcr" in text
