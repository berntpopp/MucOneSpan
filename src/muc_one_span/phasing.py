"""Conservative phase decisions from selected-sample VCF evidence.

Phase sets establish relationships between retained heterozygous records only;
a sparse VCF cannot establish whole-sequence identity or reference confidence.
"""

from __future__ import annotations

import re


def phase_evidence(variants: list[dict]) -> dict:
    """Describe supported consensus candidates without guessing cross-site phase."""
    state: dict = {
        "phase_status": "no_informative_heterozygosity",
        "genotype_status": "no_retained_variants" if not variants else "no_heterozygosity_observed",
        "haplotypes": ["I"],
        "sequence_identity_status": "unresolved",
        "reference_confidence": "unverified",
    }
    previous_end: dict[str, int] = {}
    for variant in sorted(variants, key=lambda v: (v["chrom"], v["pos"])):
        chrom, start = variant["chrom"], variant["pos"] - 1
        if start < previous_end.get(chrom, 0):
            state.update(phase_status="conflicting_variant_records", genotype_status="conflicting")
            return state
        previous_end[chrom] = start + len(variant["ref"])
    genotypes = [re.split(r"[/|]", v["genotype"]) for v in variants]
    if any("." in gt for gt in genotypes):
        state.update(phase_status="missing_genotype", genotype_status="incomplete")
        return state
    if any(len(gt) != 2 for gt in genotypes):
        state.update(phase_status="non_diploid", genotype_status="non_diploid")
        return state
    hets = [v for v, gt in zip(variants, genotypes, strict=True) if gt[0] != gt[1]]
    if not hets:
        return state
    state["genotype_status"] = "heterozygosity_observed"
    if len(hets) == 1:
        if "|" in hets[0]["genotype"] and hets[0].get("phase_set") not in (None, ".", ""):
            state.update(phase_status="phased", haplotypes=[1, 2])
        else:
            state.update(phase_status="single_heterozygous_unordered", haplotypes=[1, 2])
    elif any("|" not in v["genotype"] for v in hets):
        state["phase_status"] = "unphased"
    elif any(v.get("phase_set") in (None, ".", "") for v in hets):
        state["phase_status"] = "missing_phase_set"
    elif len({(v["chrom"], v["phase_set"]) for v in hets}) != 1:
        state["phase_status"] = "disconnected_phase_sets"
    else:
        state.update(phase_status="phased", haplotypes=[1, 2])
    return state


def annotate_consensus_candidate(
    allele_info: dict, evidence: dict, haplotype: str | int, sample: str | None, vcf_path: str
) -> None:
    """Attach evidence fields; distinguish genotype sources from molecule counts.

    Independent haplotype evidence identifies two genotype-supported sequences,
    not independent read sets: both retain the same variant observation group.
    """
    allele_info.update({k: v for k, v in evidence.items() if k != "haplotypes"})
    allele_info.update(
        consensus_haplotype=haplotype,
        consensus_sample=sample,
        consensus_policy="genotype_iupac_candidate" if haplotype == "I" else "genotype_haplotype",
        vcf_path=vcf_path,
        variant_observation_group=vcf_path,
        sequence_source=f"{vcf_path}:GT{haplotype}",
        independent_haplotype_evidence=(
            haplotype in (1, 2)
            and evidence["phase_status"] in ("phased", "single_heterozygous_unordered")
        ),
        reconstruction_status="candidate_reference_confidence_unverified",
    )
