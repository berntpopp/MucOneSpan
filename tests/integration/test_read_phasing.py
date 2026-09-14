"""Actual read-backed phase, primary record identity and exact consensus recovery."""

from __future__ import annotations

import json
import random
import shutil
from collections import defaultdict
from pathlib import Path

import pytest

from muc_one_span.phasing import phase_evidence
from muc_one_span.read_phasing import phase_same_length_reads
from muc_one_span.tools import run_tool
from muc_one_span.vcf import parse_vcf_genotypes
from tests.conftest import requires_bcftools, requires_samtools

pytestmark = [
    pytest.mark.integration,
    requires_bcftools,
    requires_samtools,
    pytest.mark.skipif(shutil.which("whatshap") is None, reason="whatshap is not installed"),
]


def read_fixture(
    directory: Path,
    variants: list[tuple[int, str, str, str, str]],
    haplotypes: list[tuple[int, ...]],
    mode: str = "full",
    collisions: bool = True,
) -> tuple[Path, Path, Path, list[str], dict[int, int]]:
    """Create known primary alignments; source identity comes from record creation."""
    directory.mkdir(parents=True)
    rng = random.Random(94810)
    bases = list("".join(rng.choices("ACGT", k=400)))
    for pos, ref, _alt, _gt, _ps in variants:
        assert len(ref) == 1
        bases[pos - 1] = ref
    reference = "".join(bases)
    reference_path = directory / "ref.fa"
    reference_path.write_text(">c\n" + reference + "\n")
    run_tool(["samtools", "faidx", str(reference_path)])
    input_vcf = directory / "input.vcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=c,length=400>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
        + "".join(
            f"c\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT:PS\t{gt}:{ps}\n"
            for pos, ref, alt, gt, ps in variants
        )
    )
    rows, truth = [], []
    intervals = {
        "full": [(0, 400)],
        "partial": [(75, 225)],
        "disjoint": [(75, 125), (175, 225)],
        "empty": [],
    }[mode]
    for haplotype, indices in enumerate(haplotypes):
        choices = {
            pos - 1: [ref, *alt.split(",")][index]
            for (pos, ref, alt, _gt, _ps), index in zip(variants, indices, strict=True)
        }
        truth.append("".join(choices.get(pos, base) for pos, base in enumerate(reference)))
        for start, end in intervals:
            sequence, cigar, last = "", [], start
            for pos, allele in choices.items():
                if not start <= pos < end:
                    continue
                sequence += reference[last:pos] + allele
                cigar.append(f"{pos - last + 1}M")
                if len(allele) > 1:
                    cigar.append(f"{len(allele) - 1}I")
                last = pos + 1
            sequence += reference[last:end]
            if last < end:
                cigar.append(f"{end - last}M")
            for ordinal in range(10):
                name = "collision" if collisions else f"hap{haplotype}_read{ordinal}_{start}"
                row = (
                    f"{name}\t0\tc\t{start + 1}\t60\t{''.join(cigar)}\t*\t0\t0\t"
                    f"{sequence}\t{'I' * len(sequence)}\tNM:i:0"
                )
                rows.append((start, haplotype, row))
    rows.sort(key=lambda row: row[0])
    sam = directory / "input.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:400\n" + "".join(row[2] + "\n" for row in rows)
    )
    bam = directory / "input.bam"
    run_tool(["samtools", "view", "-b", "-o", str(bam), str(sam)])
    run_tool(["samtools", "index", str(bam)])
    sources = {ordinal: row[1] for ordinal, row in enumerate(rows, start=1)}
    return reference_path, input_vcf, bam, truth, sources


def consensus_pair(reference: Path, vcf: Path, directory: Path) -> list[str]:
    if vcf.suffix != ".gz":
        compressed = directory / "selected.vcf.gz"
        run_tool(["bcftools", "view", "-Oz", "-o", str(compressed), str(vcf)])
        run_tool(["bcftools", "index", str(compressed)])
        vcf = compressed
    sequences = []
    for haplotype in (1, 2):
        fasta = run_tool(
            [
                "bcftools",
                "consensus",
                "-f",
                str(reference),
                "-s",
                "sample",
                "-H",
                str(haplotype),
                str(vcf),
            ]
        )
        sequences.append("".join(fasta.splitlines()[1:]))
    return sequences


@pytest.mark.parametrize(
    "alleles,genotypes,haplotypes,expected_status",
    [
        (["T", "G"], ["0/1", "0/1"], [(0, 0), (1, 1)], "phased"),
        (["T", "G"], ["0/1", "0/1"], [(0, 1), (1, 0)], "phased"),
        (["AT", "CGG"], ["0/1", "0/1"], [(0, 0), (1, 1)], "phased"),
        (["AT", "CGG"], ["0/1", "0/1"], [(0, 1), (1, 0)], "phased"),
        (["AT", "CGG"], ["1/1", "0/1"], [(1, 0), (1, 1)], "not_needed"),
        (["T", "G"], ["1/1", "0/1"], [(1, 0), (1, 1)], "not_needed"),
    ],
    ids=["cis_snps", "trans_snps", "cis_indels", "trans_indels", "shared_indels", "shared_snps"],
)
def test_exact_consensus_with_colliding_independent_records(
    tmp_path, alleles, genotypes, haplotypes, expected_status
):
    variants = [
        (100, "A", alleles[0], genotypes[0], "."),
        (200, "C", alleles[1], genotypes[1], "."),
    ]
    ref, vcf, bam, truth, sources = read_fixture(tmp_path / "collision", variants, haplotypes)
    selected, metadata = phase_same_length_reads(vcf, bam, ref, tmp_path / "phase")
    assert metadata["status"] == expected_status
    assert sorted(consensus_pair(ref, selected, tmp_path)) == sorted(truth)
    assert phase_evidence(parse_vcf_genotypes(selected))["haplotypes"] == [1, 2]
    if expected_status == "phased":
        assert metadata["primary_records"] == 20
        provenance = [
            json.loads(line) for line in Path(metadata["source_map"]).read_text().splitlines()
        ]
        assert len(provenance) == 20
        assert {row["original_qname"] for row in provenance} == {"collision"}
        truth_by_name = {
            row["phase_qname"]: sources[row["input_record_ordinal"]] for row in provenance
        }
        assigned = defaultdict(set)
        for line in Path(metadata["read_list"]).read_text().splitlines():
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            assigned[fields[4]].add(truth_by_name[fields[0]])
        assert {frozenset(values) for values in assigned.values()} == {
            frozenset({0}),
            frozenset({1}),
        }
        assert 2 <= metadata["phasing_assigned_records"] <= 20
    # Name sensitivity may flip haplotype labels, but not the unordered pair.
    unique_ref, unique_vcf, unique_bam, _, _ = read_fixture(
        tmp_path / "unique", variants, haplotypes, collisions=False
    )
    unique_selected, _ = phase_same_length_reads(
        unique_vcf, unique_bam, unique_ref, tmp_path / "unique_phase"
    )
    assert sorted(consensus_pair(unique_ref, unique_selected, tmp_path / "unique")) == sorted(truth)


@pytest.mark.parametrize(
    "mode,genotypes,phase_sets,status,phase_status",
    [
        ("partial", ["0/1", "0/1"], [".", "."], "phased", "phased"),
        ("disjoint", ["0/1", "0/1"], [".", "."], "unresolved", "unphased"),
        ("empty", ["0/1", "0/1"], [".", "."], "insufficient_reads", "unphased"),
        ("full", ["0/.", "0/1"], [".", "."], "not_needed", "missing_genotype"),
        ("full", ["0|1", "1|0"], ["7", "7"], "not_needed", "phased"),
        ("full", ["0|1", "1|0"], ["7", "8"], "phased", "phased"),
        ("disjoint", ["0|1", "1|0"], ["7", "8"], "unresolved", "disconnected_phase_sets"),
        ("full", ["0|1", "1|0"], [".", "."], "phased", "phased"),
    ],
)
def test_only_read_linkage_can_resolve_missing_or_disconnected_phase(
    tmp_path, mode, genotypes, phase_sets, status, phase_status
):
    variants = [
        (100, "A", "T", genotypes[0], phase_sets[0]),
        (200, "C", "G", genotypes[1], phase_sets[1]),
    ]
    ref, vcf, bam, truth, _ = read_fixture(tmp_path / "input", variants, [(0, 1), (1, 0)], mode)
    selected, metadata = phase_same_length_reads(vcf, bam, ref, tmp_path / "phase")
    assert metadata["status"] == status
    evidence = phase_evidence(parse_vcf_genotypes(selected))
    assert evidence["phase_status"] == phase_status
    assert evidence["reference_confidence"] == "unverified"
    if status == "phased":
        assert sorted(consensus_pair(ref, selected, tmp_path)) == sorted(truth)
    else:
        assert selected == vcf


@pytest.mark.parametrize("linked", [False, True])
def test_multiallelic_indices_are_preserved_and_unsupported_phase_stays_unresolved(
    tmp_path, linked
):
    variants = [(100, "A", "T,G", "1/2", ".")]
    haplotypes = [(1,), (2,)]
    if linked:
        variants.append((200, "C", "G", "0/1", "."))
        haplotypes = [(1, 0), (2, 1)]
    ref, vcf, bam, truth, _ = read_fixture(tmp_path / "input", variants, haplotypes)
    selected, metadata = phase_same_length_reads(vcf, bam, ref, tmp_path / "phase")
    parsed = parse_vcf_genotypes(selected)
    assert sorted(parsed[0]["genotype"].replace("|", "/").split("/")) == ["1", "2"]
    if linked and metadata["status"] == "phased":
        # New tool versions may add multiallelic support; they must recover truth.
        assert sorted(consensus_pair(ref, selected, tmp_path)) == sorted(truth)
    elif linked:
        assert metadata["status"] == "unresolved"
        assert selected == vcf
        assert phase_evidence(parsed)["haplotypes"] == ["I"]
    else:
        assert metadata["status"] == "not_needed"
        assert sorted(consensus_pair(ref, selected, tmp_path)) == sorted(truth)


def test_absent_variants_cannot_prove_identical_haplotypes(tmp_path):
    ref, vcf, bam, _, _ = read_fixture(tmp_path / "input", [], [(), ()])
    selected, metadata = phase_same_length_reads(vcf, bam, ref, tmp_path / "phase")
    assert selected == vcf
    assert metadata["status"] == "not_needed"
    evidence = phase_evidence(parse_vcf_genotypes(selected))
    assert evidence["haplotypes"] == ["I"]
    assert evidence["sequence_identity_status"] == "unresolved"
