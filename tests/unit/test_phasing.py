"""Scientific phase decisions, independent of external calling tools."""

import pytest

from muc_one_span import calling


def variants(*genotypes, phase_sets=None):
    return [
        {
            "chrom": "c",
            "pos": i + 1,
            "ref": "A",
            "alt": "C,G",
            "genotype": gt,
            "phase_set": phase_sets[i] if phase_sets else None,
        }
        for i, gt in enumerate(genotypes)
    ]


@pytest.mark.parametrize(
    "records,status,haplotypes",
    [
        ([], "no_informative_heterozygosity", ["I"]),
        (variants("2/2"), "no_informative_heterozygosity", ["I"]),
        (variants("1/2"), "single_heterozygous_unordered", [1, 2]),
        (variants("1|0", "0|1", phase_sets=["8", "8"]), "phased", [1, 2]),
        (variants("1|0", "0|1", phase_sets=["8", "9"]), "disconnected_phase_sets", ["I"]),
        (variants("1|0", "0|1"), "missing_phase_set", ["I"]),
        (variants("0/1", "0/1"), "unphased", ["I"]),
        (variants("0/.", "1|0"), "missing_genotype", ["I"]),
        (variants("1"), "non_diploid", ["I"]),
    ],
)
def test_same_length_retains_only_evidenced_haplotype_candidates(
    tmp_path, monkeypatch, records, status, haplotypes
):
    monkeypatch.setattr(calling, "_extract_and_remap_reads", lambda *a, **kw: tmp_path / "bam")
    monkeypatch.setattr(calling, "run_clair3", lambda *a, **kw: tmp_path / "raw.vcf")
    monkeypatch.setattr(calling, "filter_vcf", lambda *a, **kw: tmp_path / "calls.vcf")
    monkeypatch.setattr(calling, "parse_vcf_genotypes", lambda *a, **kw: records)
    monkeypatch.setattr(
        calling, "phase_same_length_reads", lambda path, *a, **kw: (path, {"status": "unavailable"})
    )
    alleles = {
        "homozygous": True,
        "allele_1": {"contig_name": "c", "cluster_contigs": ["c"]},
        "allele_2": {"contig_name": "c", "cluster_contigs": ["c"]},
    }
    paths = calling.disambiguate_same_length_alleles(
        tmp_path / "bam", tmp_path / "ref", alleles, tmp_path
    )
    assert len(paths) == len(haplotypes)
    assert [alleles[k]["consensus_haplotype"] for k in paths] == haplotypes
    assert all(alleles[k]["phase_status"] == status for k in paths)
    assert alleles["sequence_identity_status"] == "unresolved"
    assert alleles["homozygous"] is False


def test_unequal_length_cluster_retains_residual_heterozygosity_status(tmp_path, monkeypatch):
    monkeypatch.setattr(calling, "_extract_and_remap_reads", lambda *a, **kw: tmp_path / "bam")
    monkeypatch.setattr(calling, "run_clair3", lambda *a, **kw: tmp_path / "raw.vcf")
    monkeypatch.setattr(calling, "filter_vcf", lambda *a, **kw: tmp_path / "calls.vcf")
    monkeypatch.setattr(calling, "parse_vcf_genotypes", lambda *a, **kw: variants("0/1", "0/1"))
    alleles = {"homozygous": False, "allele_1": {"length": 60, "contig_name": "c"}}
    paths = calling.call_variants_per_allele(tmp_path / "bam", tmp_path / "ref", alleles, tmp_path)
    assert alleles["allele_1"]["phase_status"] == "unphased"
    assert alleles["allele_1"]["consensus_haplotype"] == 1
    assert alleles["allele_1"]["vcf_path"] == str(paths["allele_1"])


def test_overlapping_records_cannot_establish_confident_phase(tmp_path, monkeypatch):
    from muc_one_span.phasing import phase_evidence

    records = [
        {"chrom": "c", "pos": 2, "ref": "AA", "alt": "A", "genotype": "1|0", "phase_set": "7"},
        {"chrom": "c", "pos": 3, "ref": "A", "alt": "G", "genotype": "1|0", "phase_set": "7"},
    ]
    evidence = phase_evidence(records)
    assert evidence["phase_status"] == "conflicting_variant_records"
    assert evidence["haplotypes"] == ["I"]
