# src/muc_one_span/ladder.py
"""Reference ladder FASTA generation for MUC1 VNTR."""

from __future__ import annotations

import logging
from pathlib import Path

from muc_one_span.config import RepeatDictionary
from muc_one_span.settings import DEFAULT_SETTINGS, ConsensusSettings, ReferenceLayoutSettings

logger = logging.getLogger(__name__)


def build_contig(
    num_repeats: int,
    repeat_dict: RepeatDictionary,
    flank_length: int | None = None,
    *,
    settings: ConsensusSettings | None = None,
    reference_layout: ReferenceLayoutSettings | None = None,
) -> dict[str, str]:
    """Build a single ladder contig for a given repeat count.

    Structure: [left_flank] [pre-repeats 1-5] [N * X] [after-repeats 6-9] [right_flank]

    Args:
        num_repeats: Number of canonical X repeats in the variable region.
        repeat_dict: Loaded repeat dictionary with sequences and flanking.
        flank_length: Explicit flank length; None uses settings (default 500 bp).
        settings: Optional consensus flank defaults.
        reference_layout: Ordered fixed repeat IDs (default 1-5 then 6-9).

    Returns:
        Dict with 'name' and 'sequence' keys.
    """
    settings = settings or DEFAULT_SETTINGS.consensus
    layout = reference_layout or DEFAULT_SETTINGS.reference_layout
    layout.validate_repeats(repeat_dict.repeats)
    flank_length = settings.flank_length if flank_length is None else flank_length
    settings.validate_flanks(
        repeat_dict.flanking_left, repeat_dict.flanking_right, flank_length=flank_length
    )
    parts: list[str] = []

    # Left flanking
    if flank_length > 0 and repeat_dict.flanking_left:
        left = repeat_dict.flanking_left[:flank_length]
        parts.append(left)

    # Selected fixed repeats before the variable region.
    for rid in layout.pre:
        parts.append(repeat_dict.repeats[rid])

    # N canonical X repeats
    x_seq = repeat_dict.repeats[repeat_dict.canonical_repeat]
    parts.append(x_seq * num_repeats)

    # Selected fixed repeats after the variable region.
    for rid in layout.after:
        parts.append(repeat_dict.repeats[rid])

    # Right flanking
    if flank_length > 0 and repeat_dict.flanking_right:
        right = repeat_dict.flanking_right[:flank_length]
        parts.append(right)

    return {
        "name": f"contig_{num_repeats}",
        "sequence": "".join(parts),
    }


def generate_ladder_fasta(
    repeat_dict: RepeatDictionary,
    output_path: Path,
    min_units: int | None = None,
    max_units: int | None = None,
    flank_length: int | None = None,
    line_width: int = 80,
    *,
    settings: ConsensusSettings | None = None,
    reference_layout: ReferenceLayoutSettings | None = None,
) -> Path:
    """Generate a multi-contig FASTA reference ladder.

    Args:
        repeat_dict: Loaded repeat dictionary.
        output_path: Where to write the FASTA file.
        min_units: Minimum number of canonical repeats (default 1).
        max_units: Maximum number of canonical repeats (default 150).
        flank_length: Flanking sequence length per side (default 500bp).
        line_width: FASTA line width (default 80).
        settings: Optional flank defaults; explicit flank_length takes precedence.
        reference_layout: Selected fixed repeats and default ladder range.
            Explicit min_units/max_units take precedence.

    Returns:
        Path to the generated FASTA file.
    """
    settings = settings or DEFAULT_SETTINGS.consensus
    settings.validate_flanks(
        repeat_dict.flanking_left, repeat_dict.flanking_right, flank_length=flank_length
    )
    layout = reference_layout or DEFAULT_SETTINGS.reference_layout
    layout.validate_repeats(repeat_dict.repeats)
    min_units = layout.min_units if min_units is None else min_units
    max_units = layout.max_units if max_units is None else max_units
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(
        "Generating reference ladder (%d-%d repeats) -> %s",
        min_units,
        max_units,
        output_path,
    )

    with output_path.open("w") as f:
        for n in range(min_units, max_units + 1):
            contig = build_contig(
                n, repeat_dict, flank_length, settings=settings, reference_layout=layout
            )
            f.write(f">{contig['name']}\n")
            seq = contig["sequence"]
            for i in range(0, len(seq), line_width):
                f.write(seq[i : i + line_width] + "\n")

    return output_path
