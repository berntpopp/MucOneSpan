# src/muc_one_span/ladder.py
"""Reference ladder FASTA generation for MUC1 VNTR."""

from __future__ import annotations

import logging
from pathlib import Path

from muc_one_span.config import RepeatDictionary

logger = logging.getLogger(__name__)


def build_contig(
    num_repeats: int,
    repeat_dict: RepeatDictionary,
    flank_length: int = 500,
) -> dict[str, str]:
    """Build a single ladder contig for a given repeat count.

    Structure: [left_flank] [pre-repeats 1-5] [N * X] [after-repeats 6-9] [right_flank]

    Args:
        num_repeats: Number of canonical X repeats in the variable region.
        repeat_dict: Loaded repeat dictionary with sequences and flanking.
        flank_length: Length of flanking sequence on each side (bp).

    Returns:
        Dict with 'name' and 'sequence' keys.
    """
    parts: list[str] = []

    # Left flanking
    if flank_length > 0 and repeat_dict.flanking_left:
        left = repeat_dict.flanking_left[:flank_length]
        parts.append(left)

    # Pre-repeats: 1, 2, 3, 4, 5
    pre_ids = ["1", "2", "3", "4", "5"]
    for rid in pre_ids:
        if rid in repeat_dict.repeats:
            parts.append(repeat_dict.repeats[rid])

    # N canonical X repeats
    x_seq = repeat_dict.repeats[repeat_dict.canonical_repeat]
    parts.append(x_seq * num_repeats)

    # After-repeats: 6, 7, 8, 9
    after_ids = ["6", "7", "8", "9"]
    for rid in after_ids:
        if rid in repeat_dict.repeats:
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
    min_units: int = 1,
    max_units: int = 150,
    flank_length: int = 500,
    line_width: int = 80,
) -> Path:
    """Generate a multi-contig FASTA reference ladder.

    Args:
        repeat_dict: Loaded repeat dictionary.
        output_path: Where to write the FASTA file.
        min_units: Minimum number of canonical repeats (default 1).
        max_units: Maximum number of canonical repeats (default 150).
        flank_length: Flanking sequence length per side (default 500bp).
        line_width: FASTA line width (default 80).

    Returns:
        Path to the generated FASTA file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(
        "Generating reference ladder (%d-%d repeats) -> %s",
        min_units,
        max_units,
        output_path,
    )

    with output_path.open("w") as f:
        for n in range(min_units, max_units + 1):
            contig = build_contig(n, repeat_dict, flank_length)
            f.write(f">{contig['name']}\n")
            seq = contig["sequence"]
            for i in range(0, len(seq), line_width):
                f.write(seq[i : i + line_width] + "\n")

    return output_path
