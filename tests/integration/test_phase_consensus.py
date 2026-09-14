"""Real bcftools genotype selection and explicit phase uncertainty."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from muc_one_span.consensus import build_consensus, build_consensus_per_allele
from muc_one_span.phasing import annotate_consensus_candidate, phase_evidence
from muc_one_span.tools import run_tool
from muc_one_span.vcf import parse_vcf_variants
from tests.conftest import requires_bcftools, requires_samtools

pytestmark = [pytest.mark.integration, requires_bcftools, requires_samtools]


def phase_fixture(tmp_path: Path, rows: list[tuple[int, str, str, str]]) -> tuple[Path, Path]:
    reference = tmp_path / "ref.fa"
    reference.write_text(">c\nAAAAAAAAAA\n")
    run_tool(["samtools", "faidx", str(reference)])
    vcf = tmp_path / "input.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=c,length=10>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        + "".join(
            f"c\t{pos}\t.\tA\t{alt}\t30\tPASS\t.\tGT:PS\t{gt}:{ps}\n" for pos, alt, gt, ps in rows
        )
    )
    compressed = tmp_path / "input.vcf.gz"
    run_tool(["bcftools", "view", "-Oz", "-o", str(compressed), str(vcf)])
    run_tool(["bcftools", "index", str(compressed)])
    return reference, compressed


def sequence(path: Path) -> str:
    return "".join(line for line in path.read_text().splitlines() if not line.startswith(">"))


@pytest.mark.parametrize(
    "rows,status,expected",
    [
        ([(2, "C", "1|0", "7"), (8, "G", "0|1", "7")], "phased", ["ACAAAAAAAA", "AAAAAAAGAA"]),
        ([(2, "C", "1|0", "7"), (8, "G", "1|0", "7")], "phased", ["ACAAAAAGAA", "AAAAAAAAAA"]),
        (
            [(2, "C", "1|1", "."), (8, "G,T", "1|2", "7")],
            "single_heterozygous_unordered",
            ["ACAAAAAGAA", "ACAAAAATAA"],
        ),
        ([(2, "C", "0/1", "."), (8, "G", "0/1", ".")], "unphased", ["AMAAAAARAA"]),
        ([(2, "C", "1|0", "7"), (8, "G", "0|1", "8")], "disconnected_phase_sets", ["AMAAAAARAA"]),
        ([(2, "C", "1|0", "."), (8, "G", "0|1", ".")], "missing_phase_set", ["AMAAAAARAA"]),
        ([(2, "C", "./.", ".")], "missing_genotype", ["ANAAAAAAAA"]),
        ([(2, "C", "1", ".")], "non_diploid", ["ACAAAAAAAA"]),
        ([], "no_informative_heterozygosity", ["AAAAAAAAAA"]),
        ([(2, "C,G", "2/2", ".")], "no_informative_heterozygosity", ["AGAAAAAAAA"]),
        ([(2, "C", "0/0", ".")], "no_informative_heterozygosity", ["AAAAAAAAAA"]),
        ([(2, "AC", "0/1", ".")], "single_heterozygous_unordered", ["AAAAAAAAAA", "AACAAAAAAAA"]),
    ],
)
def test_genotype_phase_sequences_and_context(tmp_path, rows, status, expected):
    reference, vcf = phase_fixture(tmp_path, rows)
    variants = parse_vcf_variants(vcf)
    evidence = phase_evidence(variants)
    assert evidence["phase_status"] == status
    assert evidence["sequence_identity_status"] == "unresolved"
    alleles = {}
    for i, haplotype in enumerate(evidence["haplotypes"], start=1):
        key = f"allele_{i}"
        alleles[key] = {"contig_name": "c", "length": 1}
        annotate_consensus_candidate(alleles[key], evidence, haplotype, "SAMPLE", str(vcf))
    paths = build_consensus_per_allele(
        reference, dict.fromkeys(alleles, vcf), alleles, tmp_path / "consensus", flank_length=0
    )
    assert [sequence(path) for path in paths.values()] == expected
    for key, info in alleles.items():
        context = info["consensus_context"]
        assert context["vcf_path"] == str(vcf.resolve())
        assert context["haplotype"] == info["consensus_haplotype"]
        assert info["variant_observation_group"] == str(vcf)
        assert info["sequence_source"] == f"{vcf}:GT{info['consensus_haplotype']}"
        assert info["independent_haplotype_evidence"] is (info["consensus_haplotype"] != "I")
        assert context["trim_start"] == 0
        assert context["trim_end"] == len(sequence(paths[key]))
        assert sequence(Path(context["reference_path"])) == "AAAAAAAAAA"
        assert sequence(Path(context["full_consensus_path"])) == sequence(paths[key])
        assert (
            json.loads((tmp_path / "consensus" / f"consensus_{key}_context.json").read_text())
            == context
        )


def test_multisample_consensus_requires_selection(tmp_path):
    reference, vcf = phase_fixture(tmp_path, [(2, "C", "0/1", ".")])
    plain = tmp_path / "input.vcf"
    contents = plain.read_text().replace("\tSAMPLE\n", "\tSAMPLE\tOTHER\n")
    contents = contents.replace("\t0/1:.\n", "\t0/1:.\t1/1:.\n")
    plain.write_text(contents)
    run_tool(["bcftools", "view", "-Oz", "-o", str(vcf), str(plain)])
    run_tool(["bcftools", "index", "-f", str(vcf)])
    with pytest.raises(ValueError, match="sample"):
        build_consensus(reference, vcf, tmp_path / "ambiguous.fa")
    output = build_consensus(reference, vcf, tmp_path / "selected.fa", sample="OTHER")
    assert sequence(output) == "ACAAAAAAAA"


def test_default_and_explicit_iupac_heterozygous_indel(tmp_path):
    reference, vcf = phase_fixture(tmp_path, [(2, "AC", "0/1", ".")])
    old = run_tool(["bcftools", "consensus", "-f", str(reference), str(vcf)])
    explicit = build_consensus(reference, vcf, tmp_path / "mixed.fa")
    assert sequence(explicit) == "AACAAAAAAAA"
    assert "".join(line for line in old.splitlines() if not line.startswith(">")) == sequence(
        explicit
    )


def test_flank_heterozygosity_does_not_establish_vntr_independence(tmp_path):
    reference, vcf = phase_fixture(tmp_path, [(2, "C", "0/1", ".")])
    evidence = phase_evidence(parse_vcf_variants(vcf))
    alleles = {key: {"contig_name": "c", "length": 1} for key in ("allele_1", "allele_2")}
    for i, info in enumerate(alleles.values(), start=1):
        annotate_consensus_candidate(info, evidence, i, "SAMPLE", str(vcf))
    paths = build_consensus_per_allele(
        reference, dict.fromkeys(alleles, vcf), alleles, tmp_path / "consensus", flank_length=3
    )
    assert [sequence(p) for p in paths.values()] == ["AAAA", "AAAA"]
    assert all(not info["independent_haplotype_evidence"] for info in alleles.values())
