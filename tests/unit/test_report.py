# tests/unit/test_report.py
"""Tests for HTML report generation."""

from __future__ import annotations

import pytest

try:
    import jinja2  # noqa: F401

    _HAS_JINJA2 = True
except ImportError:
    _HAS_JINJA2 = False

pytestmark = pytest.mark.skipif(not _HAS_JINJA2, reason="jinja2 not installed")

if _HAS_JINJA2:
    from muc_one_span.report import generate_report


@pytest.fixture
def sample_summary():
    return {
        "alleles": {
            "allele_1": {
                "length": 50,
                "reads": 554,
                "canonical_repeats": 41,
                "contig_name": "contig_41",
                "cluster_contigs": ["contig_41", "contig_42"],
            },
            "allele_2": {
                "length": 60,
                "reads": 432,
                "canonical_repeats": 51,
                "contig_name": "contig_51",
                "cluster_contigs": ["contig_51", "contig_52"],
            },
            "homozygous": False,
            "same_length": False,
        },
        "classifications": {
            "allele_1": {
                "structure": "1 2 3 4 5 X X X:dupC A B 6 7 8 9",
                "mutations": [
                    {
                        "repeat_index": 8,
                        "closest_type": "X",
                        "mutation_name": "dupC",
                        "template_match": True,
                        "frameshift": True,
                        "vcf_support": True,
                        "vcf_qual": 23.4,
                        "boundary": False,
                    }
                ],
            },
            "allele_2": {
                "structure": "1 2 3 4 5 X X X X X X 6 7 8 9",
                "mutations": [],
            },
        },
        "tool_versions": {
            "minimap2": "minimap2 2.28-r1209",
            "samtools": "samtools 1.21",
        },
        "pipeline_version": "0.5.0",
    }


class TestGenerateReport:
    """Tests for the generate_report function."""

    def test_creates_html_file(self, tmp_path, sample_summary):
        """Report generates an HTML file containing sample name and project name."""
        out = tmp_path / "report.html"
        result = generate_report(sample_summary, out, sample_name="test_sample")

        assert result == out
        assert out.exists()
        content = out.read_text()
        assert "test_sample" in content
        assert "MucOneSpan" in content

    def test_report_self_contained(self, tmp_path, sample_summary):
        """Report must not reference any external URLs."""
        out = tmp_path / "report.html"
        generate_report(sample_summary, out, sample_name="test_sample")
        content = out.read_text()

        assert "http://" not in content
        assert "https://" not in content

    def test_report_includes_allele_data(self, tmp_path, sample_summary):
        """Report must contain the read counts from the summary."""
        out = tmp_path / "report.html"
        generate_report(sample_summary, out, sample_name="test_sample")
        content = out.read_text()

        assert "554" in content
        assert "432" in content

    def test_report_includes_mutations(self, tmp_path, sample_summary):
        """Report must include HGVS nomenclature, repeat form, ambiguity, and clinical tier."""
        out = tmp_path / "report.html"
        generate_report(sample_summary, out, sample_name="test_sample")
        content = out.read_text()

        assert "dupC" in content
        assert "59dupC" in content
        assert "NM_001204286.1:c.59dupC" in content
        assert "53C[7]&gt;53C[8]" in content or "53C[7]>53C[8]" in content
        assert "(53, 59)" in content
        assert "Tier A" in content

    def test_report_includes_tool_versions(self, tmp_path, sample_summary):
        """Report must display tool version strings."""
        out = tmp_path / "report.html"
        generate_report(sample_summary, out, sample_name="test_sample")
        content = out.read_text()

        assert "minimap2 2.28-r1209" in content
        assert "samtools 1.21" in content

    def test_report_creates_parent_dirs(self, tmp_path, sample_summary):
        """Report auto-creates nested parent directories."""
        out = tmp_path / "deep" / "nested" / "dir" / "report.html"
        generate_report(sample_summary, out, sample_name="test_sample")

        assert out.exists()

    def test_report_size_under_100kb(self, tmp_path, sample_summary):
        """Rendered report must be under 100 KB."""
        out = tmp_path / "report.html"
        generate_report(sample_summary, out, sample_name="test_sample")

        size_kb = out.stat().st_size / 1024
        assert size_kb < 100, f"Report is {size_kb:.1f} KB, expected < 100 KB"

    def test_report_with_tool_versions_override(self, tmp_path, sample_summary):
        """Explicit tool_versions parameter overrides summary values."""
        out = tmp_path / "report.html"
        custom_versions = {"bwa": "bwa-mem 0.7.18"}
        generate_report(
            sample_summary,
            out,
            sample_name="test_sample",
            tool_versions=custom_versions,
        )
        content = out.read_text()

        assert "bwa-mem 0.7.18" in content

    def test_report_with_detailed_repeats(self, tmp_path, sample_summary):
        """Detailed repeats section appears when data is provided."""
        detailed = {
            "allele_1": {
                "repeats": [
                    {
                        "index": 1,
                        "type": "1",
                        "classification": "pre_repeat",
                        "confidence": 0.98,
                        "edit_distance": 0,
                        "identity_pct": 100.0,
                    },
                ],
                "mutations_detected": [],
                "confidence": 0.95,
                "exact_match_pct": 0.92,
            },
        }
        out = tmp_path / "report.html"
        generate_report(
            sample_summary,
            out,
            sample_name="test_sample",
            detailed_repeats=detailed,
        )
        content = out.read_text()

        assert "Detailed Repeats" in content
        assert "Quality Metrics" in content


def test_report_exposes_evidence_without_homozygosity_claim(sample_summary, tmp_path):
    sample_summary["alleles"].update(sequence_identity_status="unresolved", phase_status="unphased")
    html = generate_report(sample_summary, tmp_path / "report.html").read_text()
    assert "Sequence identity: unresolved" in html
    assert "Phase: unphased" in html
    assert "Alignment records" in html


def test_report_distinguishes_unavailable_from_absent_support(sample_summary, tmp_path):
    mutation = sample_summary["classifications"]["allele_1"]["mutations"][0]
    mutation.update(vcf_support=False, vcf_support_status="projection_unavailable", vcf_qual=0.0)
    absent = {**mutation, "repeat_index": 9, "vcf_support_status": "absent"}
    ambiguous = {**mutation, "repeat_index": 10, "vcf_support_status": "localization_ambiguous"}
    sample_summary["classifications"]["allele_1"]["mutations"].append(absent)
    sample_summary["classifications"]["allele_1"]["mutations"].append(ambiguous)

    html = generate_report(sample_summary, tmp_path / "report.html").read_text()
    assert html.count(">Support unavailable<") == 1
    assert html.count(">Support ambiguous<") == 1
    assert html.count(">Unsupported<") == 1
