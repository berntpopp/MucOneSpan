#!/usr/bin/env python3
"""Build the bundled reference ladder FASTA."""

from pathlib import Path

from muc_one_span.config import load_repeat_dictionary
from muc_one_span.ladder import generate_ladder_fasta


def main() -> None:
    rd = load_repeat_dictionary()
    output = Path("src/muc_one_span/data/reference/reference_ladder.fa")
    print("Generating reference ladder (1-150 contigs, 500bp flanking)...")
    generate_ladder_fasta(rd, output, min_units=1, max_units=150, flank_length=500)
    print(f"Written to {output}")
    print(f"File size: {output.stat().st_size / 1024:.1f} KB")


if __name__ == "__main__":
    main()
