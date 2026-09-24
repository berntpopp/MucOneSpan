"""Real normalization must preserve exact mutation-specific sequence concordance."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from muc_one_span.classify import classify_sequence, validate_mutations_against_vcf
from muc_one_span.config import load_repeat_dictionary
from muc_one_span.consensus import build_consensus, trim_flanking
from muc_one_span.tools import run_tool
from muc_one_span.vcf import filter_vcf, parse_vcf_variants

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not shutil.which("bcftools") or not shutil.which("samtools"),
        reason="bcftools and samtools required",
    ),
]


def _named_mutation_case(name: str, genotype: str, tmp_path: Path) -> tuple[str, str, dict]:
    """Build one named mutation VCF, run real filter/consensus/trim, return sequence and result."""
    rd = load_repeat_dictionary()
    mutant, (parent, _) = next(
        (s, label) for s, label in rd.mutated_sequences.items() if label[1] == name
    )
    wild = rd.repeats[parent]
    prefix, suffix = "ACGT" + rd.repeats["X"], rd.repeats["X"] + "TGCA"
    ref, alt = prefix + wild + suffix, prefix + mutant + suffix
    left = 0
    while ref[left] == alt[left]:
        left += 1
    end_ref, end_alt = len(ref), len(alt)
    while end_ref > left and end_alt > left and ref[end_ref - 1] == alt[end_alt - 1]:
        end_ref -= 1
        end_alt -= 1
    # Add the required VCF left anchor for insertions/deletions.
    left -= 1
    rp, vp = tmp_path / "reference.fa", tmp_path / "raw.vcf"
    rp.write_text(">contig_3\n" + ref + "\n")
    run_tool(["samtools", "faidx", str(rp)])
    vp.write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID=contig_3,length={len(ref)}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"contig_3\t{left + 1}\t.\t{ref[left:end_ref]}\t{alt[left:end_alt]}\t30\tPASS\t.\tGT"
        f"\t{genotype}\n"
    )
    filtered = filter_vcf(vp, rp, tmp_path / "normalized")
    full, trimmed = tmp_path / "full.fa", tmp_path / "trimmed.fa"
    build_consensus(rp, filtered, full, sample="S", haplotype="I")
    context = {
        "reference_path": str(rp),
        "full_consensus_path": str(full),
        "chrom": "contig_3",
        "sample": "S",
        "haplotype": "I",
    }
    trim_flanking(full, 4, trimmed, context=context)
    sequence = "".join(
        line for line in trimmed.read_text().splitlines() if not line.startswith(">")
    )
    result = validate_mutations_against_vcf(
        classify_sequence(sequence, rd),
        parse_vcf_variants(filtered),
        sequence=sequence,
        repeat_dict=rd,
        consensus_context=context,
    )
    return sequence, rd.repeats["X"] + mutant + rd.repeats["X"], result


@pytest.mark.parametrize("name", ["dupC", "dupA", "insG", "insCCCC", "del18_31"])
def test_normalized_named_mutation_support(name: str, tmp_path: Path) -> None:
    sequence, expected, result = _named_mutation_case(name, "1/1", tmp_path)
    assert sequence == expected
    matching = [m for m in result["mutations_detected"] if m.get("mutation_name") == name]
    assert len(matching) == 1
    assert matching[0]["vcf_support"] is True


@pytest.mark.parametrize("name", ["dupC", "del18_31"])
def test_iupac_consensus_applies_heterozygous_indel_as_unresolved_event(
    name: str, tmp_path: Path
) -> None:
    """Guards the installed bcftools: -H I applies a 0/1 indel's ALT (verified on 1.17)."""
    sequence, expected, result = _named_mutation_case(name, "0/1", tmp_path)
    assert sequence == expected
    assert result["vcf_projection"]["status"] == "available"
    matching = [m for m in result["mutations_detected"] if m.get("mutation_name") == name]
    assert len(matching) == 1
    assert matching[0]["vcf_support_status"] == "heterozygous_genotype_unresolved"
    assert matching[0]["vcf_support"] is False
