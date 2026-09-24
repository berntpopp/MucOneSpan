"""Regressions for report evidence, status and percentage boundaries."""

import re

import pytest

from muc_one_span.report import compute_clinical_decision, generate_report
from muc_one_span.report_igv import create_locus_bed


def summary():
    return {
        "alleles": {"allele_1": {"reads": 100, "length": 20, "contig_name": "contig_11"}},
        "classifications": {"allele_1": {"mutations": []}},
    }


@pytest.mark.parametrize("value", [0, 0.5, 1, 50, 92, 100])
def test_exact_match_uses_percent_units(tmp_path, value):
    html = generate_report(
        summary(),
        tmp_path / "report.html",
        detailed_repeats={"allele_1": {"exact_match_pct": value}},
    ).read_text()
    match = re.search(
        r'aria-valuenow="([^"]+)"[^>]+aria-label="allele_1 exact match percentage"[^>]*>\s*<div[^>]+style="width:([^%]+)%"',
        html,
    )
    assert match
    assert float(match[1]) == float(match[2]) == value
    assert f"{value:.1f}%</span>" in html


@pytest.mark.parametrize("value", [None, "invalid", float("nan"), float("inf"), -1, 101, True])
def test_invalid_exact_match_is_unavailable(tmp_path, value):
    html = generate_report(
        summary(),
        tmp_path / "report.html",
        detailed_repeats={"allele_1": {"exact_match_pct": value}},
    ).read_text()
    assert 'aria-label="allele_1 exact match percentage"' not in html
    assert '<span class="metric-value">—</span>' in html


@pytest.mark.parametrize(
    "status", ["execution_failed", "interrupted", "insufficient_evidence", "running"]
)
def test_known_execution_failure_prevents_negative(status):
    data = summary()
    data["run_status"] = {"status": status}
    assert compute_clinical_decision(data)["state"] == "INCONCLUSIVE"


def test_explicit_status_overrides_stale_summary():
    data = summary()
    data["run_status"] = {"status": "completed"}
    assert (
        compute_clinical_decision(data, execution_status={"status": "execution_failed"})["state"]
        == "INCONCLUSIVE"
    )


@pytest.mark.parametrize("missing", [{}, {"repeat_form": None}, {"repeat_form": ""}])
def test_missing_nomenclature_is_unicode_dash(tmp_path, missing):
    data = summary()
    data["classifications"]["allele_1"]["mutations"] = [
        {"mutation_name": "unrecognized", **missing}
    ]
    html = generate_report(data, tmp_path / "report.html").read_text()
    assert "&amp;mdash;" not in html
    assert "<code>—</code>" in html


def test_requested_contig_must_exist(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">contig_1\nACGT\n")
    with pytest.raises(ValueError, match="contig_11"):
        create_locus_bed(fasta, tmp_path / "locus.bed", ["contig_11"])


@pytest.mark.parametrize("status", ["completed", "analysis_completed"])
def test_successful_explicit_status_overrides_stale_failure(tmp_path, status):
    data = summary()
    data["run_status"] = {"status": "execution_failed"}
    html = generate_report(
        data, tmp_path / "report.html", execution_status={"status": status}
    ).read_text()
    assert ">NEGATIVE</span>" in html
    assert 'class="execution-warning"' not in html


@pytest.mark.parametrize(
    "status", ["execution_failed", "interrupted", "insufficient_evidence", "running"]
)
def test_mutation_evidence_remains_visible_with_execution_warning(tmp_path, status):
    data = summary()
    data["classifications"]["allele_1"]["mutations"] = [
        {
            "mutation_name": "dupC",
            "repeat_index": 8,
            "frameshift": True,
            "template_match": True,
            "vcf_support": True,
            "vcf_support_status": "exact_sequence_concordance",
        }
    ]
    html = generate_report(
        data, tmp_path / "report.html", execution_status={"status": status}
    ).read_text()
    assert ">PATHOGENIC</span>" in html
    assert 'class="execution-warning"' in html
    assert "repeat_8:c.59dupC" in html


def test_legacy_status_never_inferred_from_report_directory(tmp_path):
    (tmp_path / "run_status.json").write_text('{"status":"execution_failed"}')
    html = generate_report(summary(), tmp_path / "report.html").read_text()
    assert ">NEGATIVE</span>" in html
    assert "Execution status unavailable (legacy input)" in html


@pytest.mark.parametrize("status", [{}, {"status": "unrecognized"}])
def test_unknown_explicit_status_is_not_success(status):
    assert compute_clinical_decision(summary(), execution_status=status)["state"] == "INCONCLUSIVE"


def test_missing_exact_match_is_unavailable(tmp_path):
    html = generate_report(
        summary(), tmp_path / "report.html", detailed_repeats={"allele_1": {}}
    ).read_text()
    assert 'aria-label="allele_1 exact match percentage"' not in html


@pytest.mark.parametrize("missing", [{}, {"hgvs_cdna": None}, {"hgvs_cdna": ""}])
def test_template_handles_missing_cdna_without_disabling_escape(tmp_path, monkeypatch, missing):
    # Isolate template fallbacks from enrichment, which usually synthesizes cDNA.
    monkeypatch.setattr(
        "muc_one_span.report._enrich_mutation_nomenclature", lambda mutation: mutation
    )
    data = summary()
    data["classifications"]["allele_1"]["mutations"] = [
        {
            "mutation_name": '<script>alert("x")</script>',
            "repeat_form": '<img src=x onerror="alert(1)">',
            **missing,
        }
    ]
    html = generate_report(data, tmp_path / "report.html").read_text()
    assert "<code>—</code>" in html
    assert "&lt;script&gt;" in html
    assert "&lt;img" in html
    assert "<img src=x" not in html
    assert "<script>alert" not in html


def test_duplicate_allele_contigs_yield_one_bed_locus(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">contig_11\n" + "A" * 2200 + "\n")
    bed = create_locus_bed(fasta, tmp_path / "locus.bed", ["contig_11", "contig_11"])
    assert bed.read_text() == "contig_11\t501\t1700\tMUC1_VNTR\n"


def test_huge_integer_percent_is_unavailable(tmp_path):
    html = generate_report(
        summary(),
        tmp_path / "report.html",
        detailed_repeats={"allele_1": {"exact_match_pct": 10**1000}},
    ).read_text()
    assert 'aria-label="allele_1 exact match percentage"' not in html


@pytest.mark.parametrize("status", [[], {}])
def test_malformed_status_prevents_negative(status):
    assert (
        compute_clinical_decision(summary(), execution_status={"status": status})["state"]
        == "INCONCLUSIVE"
    )
