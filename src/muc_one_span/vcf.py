"""VCF parsing and filtering utilities for muconespan."""

from __future__ import annotations

import math
import re
from pathlib import Path

from muc_one_span.tools import run_tool


def filter_vcf(
    vcf_path: Path,
    reference_path: Path,
    output_dir: Path,
    min_qual: float = 0.0,
    min_dp: int = 0,
    *,
    haploid_majority: bool = False,
    haploid_min_qual: float | None = None,
) -> Path:
    """Normalize and filter a VCF file with bcftools.

    Runs ``bcftools norm -f <reference>`` followed by
    ``bcftools view -f PASS`` (with optional quality filters) and indexes
    the result. When ``haploid_majority=True``, resolves borderline heterozygous
    calls on isolated haploid alignments by setting GT to 1/1 when AF >= 0.5
    and 0/0 otherwise, avoiding spurious IUPAC ambiguity characters in consensus.

    Args:
        vcf_path: Path to input VCF (may be gzipped).
        reference_path: Path to reference FASTA for left-normalisation.
        output_dir: Directory for output files.
        min_qual: Minimum QUAL score to keep a variant (0 = no filter).
        min_dp: Accepted for API compatibility but **not applied**.
            Clair3 HiFi places depth in FORMAT/DP (per-sample), not
            INFO/DP.  Filtering on INFO/DP would crash on Clair3 output.
            QUAL-only filtering is used instead since QUAL already
            integrates depth information.
        haploid_majority: If True, resolve heterozygous GTs on haploid alignments based on AF.

    Returns:
        Path to the filtered, indexed VCF (``variants.vcf.gz``).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    norm_vcf = output_dir / "normalized.vcf.gz"
    filtered = output_dir / "variants.vcf.gz"

    # Left-normalise indels against the reference
    run_tool(
        [
            "bcftools",
            "norm",
            "-f",
            str(reference_path),
            "-o",
            str(norm_vcf),
            "-O",
            "z",
            str(vcf_path),
        ]
    )

    # Skip quality filtering on empty VCFs (no records = no fields to filter).
    # Clair3 produces header-only VCFs for wild-type samples with no variants.
    # Note: a bgzipped header-only VCF may still have non-zero size; this
    # heuristic is conservative (may apply QUAL filter to header-only files,
    # which bcftools handles gracefully). A precise check would use
    # `bcftools view -H | wc -l` but adds an extra subprocess call.
    is_empty = not norm_vcf.exists() or norm_vcf.stat().st_size == 0

    # Build view command with PASS filter and optional quality filters
    view_cmd = [
        "bcftools",
        "view",
        "-f",
        "PASS",
    ]
    # Only add quality filter if VCF has records.
    # Use QUAL only — Clair3 HiFi uses FORMAT/DP not INFO/DP, and
    # QUAL already integrates depth/quality information.
    effective_qual = min_qual
    if haploid_majority and haploid_min_qual is not None:
        effective_qual = haploid_min_qual
    elif haploid_majority and min_qual <= 5.0 and min_qual > 0:
        effective_qual = min(min_qual, 4.0)

    if not is_empty and effective_qual > 0:
        view_cmd.extend(["-i", f"QUAL>={effective_qual}"])

    view_cmd.extend(
        [
            "-o",
            str(filtered),
            "-O",
            "z",
            str(norm_vcf),
        ]
    )
    run_tool(view_cmd)

    if haploid_majority and not is_empty and filtered.exists():
        lines = run_tool(["bcftools", "view", str(filtered)]).splitlines()
        new_lines: list[str] = []
        modified = False
        for line in lines:
            if line.startswith("#"):
                new_lines.append(line)
                continue
            fields = line.split("\t")
            if len(fields) >= 10:
                fmt = fields[8].split(":")
                sample_fields = fields[9].split(":")
                if "GT" in fmt and "AF" in fmt:
                    gt_idx = fmt.index("GT")
                    af_idx = fmt.index("AF")
                    af_str = sample_fields[af_idx]
                    try:
                        af_vals = [float(x) for x in af_str.split(",") if x != "."]
                        if af_vals:
                            max_af = max(af_vals)
                            if max_af >= 0.5:
                                alt_idx = af_vals.index(max_af) + 1
                                target_gt = f"{alt_idx}/{alt_idx}"
                            elif max_af < 0.2:
                                target_gt = "0/0"
                            else:
                                # Borderline heterozygous / mixed evidence: retain original genotype
                                target_gt = sample_fields[gt_idx]
                            if sample_fields[gt_idx] != target_gt:
                                sample_fields[gt_idx] = target_gt
                                fields[9] = ":".join(sample_fields)
                                modified = True
                    except (ValueError, IndexError):
                        pass
            new_lines.append("\t".join(fields))
        if modified:
            tmp_vcf = output_dir / "mod_haploid.vcf"
            tmp_vcf.write_text("\n".join(new_lines) + "\n")
            filtered.unlink(missing_ok=True)
            run_tool(["bcftools", "view", "-o", str(filtered), "-O", "z", str(tmp_vcf)])
            tmp_vcf.unlink(missing_ok=True)

    run_tool(["bcftools", "index", str(filtered)])

    # Remove intermediate normalized VCF
    norm_vcf.unlink(missing_ok=True)

    return filtered


def select_vcf_sample(vcf_path: Path, sample: str | None = None) -> str:
    """Select one VCF sample, rejecting missing or ambiguous selections."""
    samples = run_tool(["bcftools", "query", "-l", str(vcf_path)]).splitlines()
    if sample is not None:
        if sample not in samples:
            raise ValueError(f"VCF sample {sample!r} is not present in {vcf_path}")
        return sample
    if len(samples) != 1:
        raise ValueError(f"Select an explicit sample: VCF has {len(samples)} samples")
    return samples[0]


def parse_vcf_genotypes(vcf_path: Path, sample: str | None = None) -> list[dict]:
    """Return selected-sample variant identity, quality, genotype and phase set.

    Missing QUAL and PS are represented by None. Tool failures and malformed
    records propagate; an empty list means a successful query with no records.
    """
    return parse_vcf_variants(vcf_path, sample=sample)


def parse_vcf_variants(vcf_path: Path, sample: str | None = None) -> list[dict]:
    """Parse CHROM/POS/REF/ALT/QUAL/GT/PS without silently dropping bad records."""
    selected_sample = select_vcf_sample(vcf_path, sample)
    output = run_tool(
        [
            "bcftools",
            "query",
            "-u",
            "-f",
            "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL[\\t%GT\\t%PS]\\n",
            "-s",
            selected_sample,
            str(vcf_path),
        ]
    )
    variants: list[dict] = []
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) != 7:
            raise ValueError(f"Malformed VCF query record: {line!r}")
        chrom, position, ref, alt, quality, genotype, phase_set = fields
        pos = int(position)
        qual = None if quality == "." else float(quality)
        if not chrom or pos < 1 or not ref or not alt:
            raise ValueError(f"Invalid VCF variant identity: {line!r}")
        if qual is not None and (not math.isfinite(qual) or qual < 0):
            raise ValueError(f"Invalid VCF quality: {quality!r}")
        if re.fullmatch(r"(?:[0-9]+|\.)(?:[/|](?:[0-9]+|\.))*", genotype) is None:
            raise ValueError(f"Invalid VCF genotype: {genotype!r}")
        if "/" in genotype and "|" in genotype:
            raise ValueError(f"Mixed genotype separators: {genotype!r}")
        indices = re.split(r"[/|]", genotype)
        if any(int(index) > len(alt.split(",")) for index in indices if index != "."):
            raise ValueError(f"Genotype allele index exceeds ALT count: {genotype!r}")
        variants.append(
            {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "qual": qual,
                "genotype": genotype,
                "phase_set": None if phase_set == "." else phase_set,
                "sample": selected_sample,
            }
        )
    return variants
