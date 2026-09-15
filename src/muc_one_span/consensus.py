"""Per-allele consensus sequence generation with bcftools."""

from __future__ import annotations

import json
import logging
from pathlib import Path

from muc_one_span.config import RepeatDictionary
from muc_one_span.settings import DEFAULT_SETTINGS, ConsensusSettings, ReferenceLayoutSettings
from muc_one_span.tools import run_tool
from muc_one_span.vcf import select_vcf_sample

logger = logging.getLogger(__name__)


def build_consensus(
    reference_path: Path,
    vcf_path: Path,
    output_path: Path,
    *,
    sample: str | None = None,
    haplotype: str | int = "I",
) -> Path:
    """Generate a consensus FASTA by applying VCF variants to a reference.

    Uses an explicit sample and genotype-aware haplotype policy. The default
    ``I`` produces a mixed IUPAC candidate (indel phase remains unresolved).
    Integer 1/2 selects the actual genotype allele index for evidenced haplotypes.
    Multiple VCF samples require explicit selection.

    Args:
        reference_path: Path to the reference FASTA (typically the matching
            ladder contig extracted with ``samtools faidx``).
        vcf_path: Path to the filtered VCF (may be gzipped).
        output_path: Destination path for the consensus FASTA.
        sample: Selected VCF sample; inferred only for a single-sample file.
        haplotype: ``I`` for a genotype-aware mixed candidate, or GT position 1/2.

    Returns:
        Path to the written consensus FASTA.
    """
    if haplotype not in ("I", 1, 2):
        raise ValueError("haplotype must be 'I', 1, or 2")
    selected_sample = select_vcf_sample(vcf_path, sample)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    stdout = run_tool(
        [
            "bcftools",
            "consensus",
            "-s",
            selected_sample,
            "-H",
            str(haplotype),
            "-M",
            "N",
            "-f",
            str(reference_path),
            str(vcf_path),
        ]
    )

    output_path.write_text(stdout)
    return output_path


def _find_anchor(
    sequence: str,
    anchor: str,
    expected_pos: int,
    tolerance: int = DEFAULT_SETTINGS.consensus.anchor_tolerance,
) -> int | None:
    """Find an anchor whose start is near ``expected_pos``.

    Searches for an exact substring match with its start within the inclusive
    interval ``expected_pos +/- tolerance``.
    Returns the position where the anchor ENDS (i.e., the start of the region
    after the anchor).

    Returns None if not found.
    """
    search_start = max(0, expected_pos - tolerance)
    search_end = min(len(sequence), expected_pos + tolerance + len(anchor))
    region = sequence[search_start:search_end]
    idx = region.find(anchor)
    if idx >= 0:
        return search_start + idx + len(anchor)
    return None


def trim_flanking(
    consensus_fasta: Path,
    flank_length: int | None,
    output_path: Path,
    repeat_dict: RepeatDictionary | None = None,
    *,
    context: dict | None = None,
    settings: ConsensusSettings | None = None,
    reference_layout: ReferenceLayoutSettings | None = None,
) -> Path:
    """Remove flanking sequences from a consensus FASTA, keeping only the VNTR.

    The ladder contigs have structure:
    ``[left_flank] [pre-repeats] [N * X] [after-repeats] [right_flank]``

    When *repeat_dict* is provided, uses anchor-based boundary detection
    that is resilient to indels in the flanking region (e.g. from Clair3
    false positives). The flank prefixes match ladder construction. Anchor
    components use up to ``anchor_bases`` available bases; trim coordinates use
    their actual lengths. ``anchor_tolerance`` is the permitted displacement of
    the VNTR boundary, so zero accepts an exact expected boundary. Falls back
    to fixed-position trim when an anchor is absent or no dictionary is supplied.

    Args:
        consensus_fasta: Path to the full-contig consensus FASTA.
        flank_length: Number of flanking bp on each side (must match the
            value used during ladder generation); None uses settings (default 500).
        output_path: Destination path for the trimmed FASTA.
        repeat_dict: Optional repeat dictionary for anchor-based trimming.
        settings: Optional flank and exact-anchor search parameters.
        reference_layout: Selected fixed repeat IDs defining the boundary anchors.
        context: Optional output dictionary receiving actual 0-based half-open
            trim_start/trim_end coordinates on the full consensus.

    Returns:
        Path to the trimmed FASTA containing only the VNTR region.
    """
    settings = settings or DEFAULT_SETTINGS.consensus
    layout = reference_layout or DEFAULT_SETTINGS.reference_layout
    if repeat_dict is not None:
        layout.validate_repeats(repeat_dict.repeats)
    flank_length = settings.flank_length if flank_length is None else flank_length
    if repeat_dict is not None:
        settings.validate_flanks(
            repeat_dict.flanking_left, repeat_dict.flanking_right, flank_length=flank_length
        )
    lines = consensus_fasta.read_text().strip().splitlines()
    header = lines[0] if lines and lines[0].startswith(">") else ">consensus"
    sequence = "".join(line for line in lines if not line.startswith(">"))

    # Validate: if sequence is too short to trim, return it untrimmed
    if len(sequence) < 2 * flank_length:
        logger.warning(
            "Sequence length %d is shorter than 2 * flank_length (%d); "
            "returning full sequence untrimmed.",
            len(sequence),
            2 * flank_length,
        )
        if context is not None:
            context.update(
                trim_start=0,
                trim_end=len(sequence),
                left_trim_method="untrimmed_short_sequence",
                right_trim_method="untrimmed_short_sequence",
            )
        output_path.write_text(f"{header}_vntr\n{sequence}\n")
        return output_path

    # Default: fixed-position trim
    left_trim = flank_length
    right_trim = len(sequence) - flank_length if flank_length > 0 else len(sequence)

    left_method = right_method = "fixed"
    # Try anchor-based trimming if repeat_dict provided
    if repeat_dict is not None and flank_length > 0:
        left_method = right_method = "fixed_anchor_not_found"
        anchor_bases = settings.anchor_bases
        # Use the same flank sequence as ladder.build_contig.
        if settings.proximal_flank:
            left_flank_seq = repeat_dict.flanking_left[-flank_length:]
        else:
            left_flank_seq = repeat_dict.flanking_left[:flank_length]
        left_flank_part = left_flank_seq[-anchor_bases:]
        first_repeat_part = repeat_dict.repeats[layout.left_anchor_id][:anchor_bases]
        left_anchor = left_flank_part + first_repeat_part
        anchor_pos = _find_anchor(
            sequence, left_anchor, flank_length - len(left_flank_part), settings.anchor_tolerance
        )
        if anchor_pos is not None:
            left_trim = anchor_pos - len(first_repeat_part)
            left_method = "exact_anchor"

        last_repeat_part = repeat_dict.repeats[layout.right_anchor_id][-anchor_bases:]
        right_flank_part = repeat_dict.flanking_right[:flank_length][:anchor_bases]
        right_anchor = last_repeat_part + right_flank_part
        right_anchor_pos = _find_anchor(
            sequence,
            right_anchor,
            len(sequence) - flank_length - len(last_repeat_part),
            settings.anchor_tolerance,
        )
        if right_anchor_pos is not None:
            right_trim = right_anchor_pos - len(right_flank_part)
            right_method = "exact_anchor"

    if context is not None:
        context.update(
            trim_start=left_trim,
            trim_end=right_trim,
            left_trim_method=left_method,
            right_trim_method=right_method,
        )
    vntr = sequence[left_trim:right_trim]

    output_path.write_text(f"{header}_vntr\n{vntr}\n")
    return output_path


def build_consensus_per_allele(
    reference_path: Path,
    vcf_paths: dict[str, Path],
    alleles: dict,
    output_dir: Path,
    flank_length: int | None = None,
    repeat_dict: RepeatDictionary | None = None,
    *,
    settings: ConsensusSettings | None = None,
    reference_layout: ReferenceLayoutSettings | None = None,
) -> dict[str, Path]:
    """Build a consensus FASTA for each detected allele and trim flanking.

    For each allele, the function:

    1. Extracts the matching single contig from the ladder reference with
       ``samtools faidx`` (using the peak contig name from allele detection).
    2. Indexes the single-contig FASTA.
    3. Calls :func:`build_consensus` to apply the allele's VCF variants.
    4. Trims flanking sequences to isolate the VNTR region.

    Args:
        reference_path: Path to the full ladder reference FASTA.
        vcf_paths: Mapping from allele key (e.g. ``"allele_1"``) to VCF path.
        alleles: Allele detection result from
            :func:`~muc_one_span.alleles.detect_alleles`.
        output_dir: Base output directory for consensus files.
        flank_length: Explicit flank length; None uses settings (default 500 bp).
        repeat_dict: Optional repeat dictionary for exact boundary anchors.
        settings: Optional flank and anchor parameters.
        reference_layout: Selected fixed repeats defining the boundary anchors.

    Returns:
        Dictionary mapping allele key to the trimmed consensus FASTA path
        (containing only the VNTR region, ready for repeat classification).
    """
    settings = settings or DEFAULT_SETTINGS.consensus
    if repeat_dict is not None:
        settings.validate_flanks(
            repeat_dict.flanking_left, repeat_dict.flanking_right, flank_length=flank_length
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, Path] = {}

    for allele_key, vcf_path in vcf_paths.items():
        allele_info = alleles[allele_key]
        # Use contig_name from allele detection (canonical repeat count),
        # not total length.  Fall back for backwards compatibility.
        contig_name = allele_info.get("contig_name", f"contig_{allele_info['length']}")

        # Extract the single matching contig from the full ladder reference
        contig_fa = output_dir / f"ref_{contig_name}.fa"
        stdout = run_tool(
            [
                "samtools",
                "faidx",
                str(reference_path),
                contig_name,
            ]
        )
        contig_fa.write_text(stdout)

        # Index the single-contig reference so bcftools consensus can use it
        run_tool(["samtools", "faidx", str(contig_fa)])

        # Build full consensus (with flanking)
        full_consensus = output_dir / f"consensus_{allele_key}_full.fa"
        sample = select_vcf_sample(vcf_path, allele_info.get("consensus_sample"))
        haplotype = allele_info.get("consensus_haplotype", "I")
        build_consensus(contig_fa, vcf_path, full_consensus, sample=sample, haplotype=haplotype)
        context = {
            "vcf_path": str(vcf_path.resolve()),
            "reference_path": str(contig_fa.resolve()),
            "full_consensus_path": str(full_consensus.resolve()),
            "chrom": contig_name,
            "sample": sample,
            "haplotype": haplotype,
        }

        # Trim flanking to get VNTR-only sequence for classification
        trimmed = output_dir / f"consensus_{allele_key}.fa"
        trim_flanking(
            full_consensus,
            flank_length,
            trimmed,
            repeat_dict=repeat_dict,
            context=context,
            settings=settings,
            reference_layout=reference_layout,
        )
        allele_info["consensus_context"] = context
        context_path = output_dir / f"consensus_{allele_key}_context.json"
        context_path.write_text(json.dumps(context, indent=2) + "\n")
        results[allele_key] = trimmed

    # Genotype differences outside the retained VNTR cannot establish two
    # independently reconstructed VNTR haplotypes. Preserve both candidates but
    # remove diploid recovery credit for identical trimmed sequences.
    sequences = {
        key: "".join(line for line in path.read_text().splitlines() if not line.startswith(">"))
        for key, path in results.items()
    }
    for key, sequence in sequences.items():
        same = [other for other, value in sequences.items() if value == sequence]
        if len(same) > 1:
            alleles[key]["independent_haplotype_evidence"] = False
            alleles[key]["vntr_phase_status"] = "no_sequence_distinction"
        elif alleles[key].get("independent_haplotype_evidence"):
            alleles[key]["vntr_phase_status"] = "distinct_genotype_candidates"
        else:
            alleles[key]["vntr_phase_status"] = "unresolved"
    return results
