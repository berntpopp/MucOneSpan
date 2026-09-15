"""Unit tests for offline IGV report assets and generator."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from muc_one_span.report import generate_report
from muc_one_span.report_assets import (
    IGV_ASSET_PATH,
    IGV_GZIP_SHA256,
    IGV_REPORT_LIBRARY_MARKER,
    IGV_SHA256,
    IGV_VERSION,
    REPORT_IGV_EMBEDDED,
    REPORT_IGV_OFF,
    REPORT_IGV_SIDECAR,
    extract_igv_fragments,
    igv_payload,
    igv_provenance,
    js_json_literal,
    verify_asset,
)
from muc_one_span.report_igv import (
    build_igv_context,
    create_locus_bed,
    run_igv_report,
)


def test_verify_asset_success():
    """Vendored asset passes SHA-256 verification."""
    verify_asset(IGV_ASSET_PATH, IGV_GZIP_SHA256)


def test_verify_asset_failure(tmp_path: Path):
    """Corrupted or missing asset raises ValueError."""
    fake_asset = tmp_path / "corrupt.gz"
    fake_asset.write_bytes(b"corrupted data")
    with pytest.raises(ValueError, match="does not match pinned SHA-256"):
        verify_asset(fake_asset, IGV_GZIP_SHA256)

    missing = tmp_path / "missing.gz"
    with pytest.raises(ValueError, match="is missing"):
        verify_asset(missing, IGV_GZIP_SHA256)


def test_igv_payload_modes():
    """igv_payload returns base64 string only for embedded mode."""
    payload = igv_payload(REPORT_IGV_EMBEDDED)
    assert payload is not None
    assert len(payload) > 100000

    assert igv_payload(REPORT_IGV_SIDECAR) is None
    assert igv_payload(REPORT_IGV_OFF) is None

    with pytest.raises(ValueError, match="Unknown --report-igv mode"):
        igv_payload("invalid_mode")


def test_igv_provenance():
    """igv_provenance describes version and location."""
    prov_emb = igv_provenance(REPORT_IGV_EMBEDDED)
    assert IGV_VERSION in prov_emb
    assert IGV_SHA256 in prov_emb
    assert "embedded in this file" in prov_emb

    prov_side = igv_provenance(REPORT_IGV_SIDECAR)
    assert "separate HTML sidecar" in prov_side

    prov_off = igv_provenance(REPORT_IGV_OFF)
    assert prov_off == "not included (--report-igv off)"


def test_js_json_literal():
    """js_json_literal safely escapes '<' characters."""
    assert js_json_literal('{"tag": "<script>alert(1)</script>"}', "{}") == (
        '{"tag":"\\u003cscript>alert(1)\\u003c/script>"}'
    )
    assert js_json_literal("  [1, 2, 3];  ", "[]") == "[1,2,3]"
    assert js_json_literal("", "fallback") == "fallback"
    assert js_json_literal("invalid json", "fallback") == "fallback"


def test_extract_igv_fragments():
    """extract_igv_fragments extracts container, tableJson, and sessionDictionary."""
    html_sample = (
        "<html><body>\n"
        '<div id="container">\n'
        "  <p>IGV Viewer</p>\n"
        "</div>\n"
        "</body>\n"
        '<script>\nconst tableJson = [{"id": 1}];\nconst sessionDictionary = {"0": "sess"};\n</script>'
    )
    content, t_json, s_dict = extract_igv_fragments(html_sample)
    assert '<div id="container">' in content
    assert t_json == '[{"id": 1}];'
    assert s_dict == '{"0": "sess"};'


def test_extract_igv_fragments_missing():
    """extract_igv_fragments handles missing container gracefully."""
    content, t_json, s_dict = extract_igv_fragments("<html>no container</html>")
    assert content == ""
    assert t_json == ""
    assert s_dict == ""


def test_create_locus_bed(tmp_path: Path):
    """create_locus_bed writes a valid 1-interval BED file."""
    fasta = tmp_path / "test.fa"
    fasta.write_text(">chr1_vntr description\nACGTACGT\nACGT\n")
    bed = tmp_path / "locus.bed"

    create_locus_bed(fasta, bed, locus_name="MUC1_REGION")
    assert bed.read_text() == "chr1_vntr\t0\t12\tMUC1_REGION\n"


def test_run_igv_report_invalid_mode(tmp_path: Path):
    """run_igv_report raises on invalid mode."""
    with pytest.raises(ValueError, match="cannot generate an IGV report"):
        run_igv_report(
            tmp_path / "b.bed", tmp_path / "f.fa", tmp_path / "o.html", report_igv="invalid"
        )


def test_run_igv_report_sidecar(tmp_path: Path):
    """run_igv_report sidecar substitutes verified library marker."""
    bed = tmp_path / "locus.bed"
    bed.write_text("c\t0\t10\tL\n")
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">c\nACGTACGTAC\n")
    out_html = tmp_path / "igv_sidecar.html"

    fake_generated = (
        "<html><body>"
        '<div id="container">track</div>'
        "</body>"
        f"<script>{IGV_REPORT_LIBRARY_MARKER}</script>"
        "const tableJson = []\n"
        "const sessionDictionary = {}\n"
        "</html>"
    )

    def fake_run_tool(cmd):
        out_html.write_text(fake_generated, encoding="utf-8")
        return ""

    with patch("muc_one_span.report_igv.run_tool", side_effect=fake_run_tool):
        res = run_igv_report(bed, fasta, out_html, report_igv=REPORT_IGV_SIDECAR)

    assert res == out_html
    text = out_html.read_text(encoding="utf-8")
    assert IGV_REPORT_LIBRARY_MARKER not in text
    assert "igv" in text


def test_build_igv_context_off(tmp_path: Path):
    """build_igv_context returns empty context when report_igv is off."""
    ctx = build_igv_context(tmp_path, tmp_path / "ref.fa", report_igv=REPORT_IGV_OFF)
    assert ctx["mode"] == REPORT_IGV_OFF
    assert ctx["has_igv"] is False


def test_generate_report_with_igv_and_hgvs(tmp_path: Path):
    """generate_report formats HGVS nomenclature and renders embedded IGV."""
    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(">contig_51\n" + "ACGT" * 50 + "\n")

    summary = {
        "alleles": {
            "allele_1": {"length": 51, "reads": 100, "contig_name": "contig_51"},
        },
        "classifications": {
            "allele_1": {
                "structure": "1 2 3 4 5 X 59dupC 6 7 8 9",
                "mutations": [
                    {
                        "repeat_index": 7,
                        "mutation_name": "59dupC",
                        "frameshift": True,
                        "vcf_support": True,
                    }
                ],
            }
        },
        "tool_versions": {"samtools": "1.20"},
    }

    fake_igv_html = (
        "<html><body>"
        '<div id="container"><div id="igvDiv"></div></div>'
        "</body>\n"
        'const tableJson = [{"chrom": "contig_51"}];\n'
        'const sessionDictionary = {"contig_51": "data:application/json;base64,eyJ2ZXJzaW9uIjoiMSJ9"};\n'
        "</html>"
    )

    def fake_run_tool(cmd):
        for i, arg in enumerate(cmd):
            if arg == "--output":
                Path(cmd[i + 1]).write_text(fake_igv_html, encoding="utf-8")
        return ""

    out_report = tmp_path / "report.html"
    with patch("muc_one_span.report_igv.run_tool", side_effect=fake_run_tool):
        generate_report(
            summary,
            out_report,
            sample_name="sample_test",
            report_igv=REPORT_IGV_EMBEDDED,
            fasta_path=ref_fa,
        )

    html = out_report.read_text(encoding="utf-8")
    assert "repeat_7:c.59dupC" in html
    assert "Tier A" in html
    assert "Interactive Alignment Browser" in html
    assert "loadIgv()" in html
    assert "DecompressionStream" in html
    assert "igv.js 3.0.2" in html


@pytest.mark.parametrize("kind", ["fasta", "bed", "bam", "vcf"])
def test_requested_missing_files_are_errors(tmp_path, kind):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">c\nACGT\n")
    arguments = {"fasta_path": fasta, f"{kind}_path": tmp_path / f"missing.{kind}"}
    with pytest.raises(FileNotFoundError, match="Requested"):
        build_igv_context(tmp_path, **arguments)


def test_plural_vcf_precedence_dedup_and_labels(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">c\nACGT\n")
    bed = tmp_path / "locus.bed"
    bed.write_text("c\t0\t4\tlocus\n")
    vcf = tmp_path / "merged.vcf"
    vcf.touch()
    import json

    configs = []

    def capture_config(cmd):
        configs.extend(json.loads(Path(cmd[cmd.index("--track-config") + 1]).read_text()))

    with patch("muc_one_span.report_igv.run_tool", side_effect=capture_config):
        run_igv_report(
            bed,
            fasta,
            tmp_path / "igv.html",
            vcf_file=tmp_path / "ignored.vcf",
            vcf_paths={"allele_2": vcf, "allele_1": vcf},
            report_igv="embedded",
        )
    assert len(configs) == 1
    assert configs[0]["name"] == "Shared variants (Allele 1, Allele 2)"


def test_explicit_bed_does_not_hide_missing_assigned_contig(tmp_path):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">contig_1\nACGT\n")
    bed = tmp_path / "locus.bed"
    bed.write_text("contig_1\t0\t4\tlocus\n")

    def fake_report(cmd):
        Path(cmd[cmd.index("--output") + 1]).write_text("<html></html>")

    with (
        patch("muc_one_span.report_igv.run_tool", side_effect=fake_report),
        pytest.raises(ValueError, match="contig_11"),
    ):
        build_igv_context(tmp_path, fasta, bed_path=bed, contig_names=["contig_11"])


def test_malformed_generated_vcf_is_rejected(tmp_path):
    import base64
    import gzip
    import json

    fasta = tmp_path / "ref.fa"
    fasta.write_text(">c\nACGT\n")
    bed = tmp_path / "locus.bed"
    bed.write_text("c\t0\t4\tlocus\n")
    vcf = tmp_path / "variants.vcf"
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    vcf.write_text("##fileformat=VCFv4.2\n" + header + "\nc\t2\t.\tA\tC\t60\tPASS\t.\tGT\t1/1\n")

    def uri(data):
        return (
            "data:application/gzip;base64,"
            + base64.b64encode(gzip.compress(data.encode())).decode()
        )

    def fake_report(cmd):
        malformed_vcf = header + "c\t2\t.\tA\tC\t60\tPASS\t.\tGT\t1/1\n"
        session = {"tracks": [{"type": "variant", "name": "Variants", "url": uri(malformed_vcf)}]}
        generated = (
            '<div id="container"></div></body>\nconst tableJson = {}\nconst sessionDictionary = '
            + json.dumps({"0": uri(json.dumps(session))})
            + "\n"
        )
        Path(cmd[cmd.index("--output") + 1]).write_text(generated)

    with (
        patch("muc_one_span.report_igv.run_tool", side_effect=fake_report),
        pytest.raises(ValueError, match=r"Malformed VCF.*create_report"),
    ):
        run_igv_report(bed, fasta, tmp_path / "igv.html", vcf_file=vcf, report_igv="embedded")
