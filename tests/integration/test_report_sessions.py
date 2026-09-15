"""Real igv-reports output verifies allele loci and variant payloads."""

import base64
import gzip
import json
import shutil

import pytest

from muc_one_span.report import generate_report
from muc_one_span.report_assets import extract_igv_fragments
from tests.report_fixture_data import fixture_files

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not shutil.which("create_report"), reason="create_report unavailable"),
]


def sessions(html):
    _, table, dictionary = extract_igv_fragments(html)
    raw = json.loads(dictionary.strip().rstrip(";"))
    decoded = [
        json.loads(gzip.decompress(base64.b64decode(uri.split(",", 1)[1]))) for uri in raw.values()
    ]
    return json.loads(table.strip().rstrip(";")), decoded


@pytest.mark.parametrize("mode", ["embedded", "sidecar"])
def test_real_report_uses_allele_contigs(tmp_path, mode):
    summary, fasta, paths = fixture_files(tmp_path)
    output = tmp_path / "report.html"
    generate_report(summary, output, report_igv=mode, fasta_path=fasta, vcf_path=paths["allele_1"])
    html = (tmp_path / "igv_report.html" if mode == "sidecar" else output).read_text()
    table, decoded = sessions(html)
    assert [row[1] for row in table["rows"]] == ["contig_11", "contig_21"]
    assert [session["locus"].split(":")[0] for session in decoded] == ["contig_11", "contig_21"]


@pytest.mark.parametrize("mode", ["embedded", "sidecar"])
def test_real_report_includes_all_vcf_tracks(tmp_path, mode):
    summary, fasta, paths = fixture_files(tmp_path)
    output = tmp_path / "report.html"
    generate_report(
        summary,
        output,
        report_igv=mode,
        fasta_path=fasta,
        vcf_paths=dict(reversed(list(paths.items()))),
    )
    html = (tmp_path / "igv_report.html" if mode == "sidecar" else output).read_text()
    _, decoded = sessions(html)
    for session in decoded:
        variants = [track for track in session["tracks"] if track["type"] == "variant"]
        assert [track["name"] for track in variants] == ["Allele 1 variants", "Allele 2 variants"]
        contig = session["locus"].split(":")[0]
        records = []
        for track in variants:
            payload = gzip.decompress(base64.b64decode(track["url"].split(",", 1)[1])).decode()
            records.extend(line for line in payload.splitlines() if not line.startswith("#"))
        assert records == [f"{contig}\t700\t.\tA\tC\t60\tPASS\t."]


@pytest.mark.parametrize("mode", ["embedded", "sidecar"])
@pytest.mark.parametrize("case", ["single", "shared", "absent", "empty_plural"])
def test_real_report_optional_shared_and_single_inputs(tmp_path, mode, case):
    summary, fasta, paths = fixture_files(tmp_path)
    extra = {}
    expected = ["Variants"]
    if case == "single":
        summary["alleles"].pop("allele_2")
        extra = {"vcf_path": paths["allele_1"]}
    elif case == "shared":
        summary["alleles"]["allele_2"]["contig_name"] = "contig_11"
        extra = {"vcf_paths": dict.fromkeys(paths, paths["allele_1"])}
        expected = ["Shared variants (Allele 1, Allele 2)"]
    elif case == "empty_plural":
        extra = {"vcf_paths": {}, "vcf_path": tmp_path / "ignored_missing.vcf"}
        expected = []
    else:
        expected = []
    output = tmp_path / "report.html"
    generate_report(summary, output, report_igv=mode, fasta_path=fasta, **extra)
    html = (tmp_path / "igv_report.html" if mode == "sidecar" else output).read_text()
    table, decoded = sessions(html)
    assert len(table["rows"]) == (1 if case in {"single", "shared"} else 2)
    for session in decoded:
        variants = [track for track in session["tracks"] if track["type"] == "variant"]
        assert [track["name"] for track in variants] == expected
