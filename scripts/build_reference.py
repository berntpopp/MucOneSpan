#!/usr/bin/env python3
"""Build the bundled reference ladder FASTA."""

from pathlib import Path

from open_pacmuci.config import load_repeat_dictionary
from open_pacmuci.ladder import generate_ladder_fasta


def main() -> None:
    rd = load_repeat_dictionary()
    output = Path("src/open_pacmuci/data/reference/reference_ladder.fa")
    print("Generating reference ladder (1-150 contigs, 500bp flanking)...")
    generate_ladder_fasta(rd, output, min_units=1, max_units=150, flank_length=500)
    print(f"Written to {output}")
    print(f"File size: {output.stat().st_size / 1024:.1f} KB")


if __name__ == "__main__":
    main()
