#!/usr/bin/env python3
"""Build and validate the 250-design specification for the 500-dataset initiative.

250 biological designs (500 haplotypes) across 7 clinical challenge categories,
stratified into 150 development, 50 validation, and 50 held-out test datasets.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
STRUCTURES_DIR = ROOT / "examples" / "structures_500"
DESIGN_FILE = ROOT / "examples" / "experiment_500_designs.json"


def generate_identical_structure(length: int) -> str:
    """Generate identical repeat unit chain for homozygous wild-type designs."""
    canonical = length - 9
    if canonical < 1:
        raise ValueError(f"Total length {length} must be >= 10")
    chain = f"1-2-3-4-5-{'-'.join(['X'] * canonical)}-6-7-8-9"
    return f"haplotype_1\t{chain}\nhaplotype_2\t{chain}\n"


def create_structure_files() -> dict[int, str]:
    """Create authored structure files for identical sequence designs."""
    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    lengths = [25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120]
    paths = {}
    for length in lengths:
        file_path = STRUCTURES_DIR / f"identical_{length}_{length}.txt"
        file_path.write_text(generate_identical_structure(length))
        paths[length] = str(file_path.relative_to(ROOT))
    return paths


def build_500_designs() -> list[dict[str, Any]]:
    """Build the 250 biological sample designs (150 dev, 50 val, 50 test)."""
    struct_paths = create_structure_files()
    designs: list[dict[str, Any]] = []

    def add(
        d_id: int,
        name: str,
        split: str,
        category: str,
        lengths: list[int],
        mutation: str | None,
        target: list[int] | None,
        templates: int = 200,
        structure_file: str | None = None,
        platform_mode: str = "all",
    ) -> None:
        bio_seed = 5_000_000 + d_id * 10
        hifi_seed = 6_000_000 + d_id * 10
        ont_seed = 7_000_000 + d_id * 10
        designs.append(
            {
                "design_id": d_id,
                "name": name,
                "split": split,
                "category": category,
                "lengths": lengths,
                "mutation": mutation,
                "targets": [target] if target else [],
                "requested_templates": templates,
                "structure_file": structure_file,
                "bio_seed": bio_seed,
                "hifi_seed": hifi_seed,
                "ont_seed": ont_seed,
                "platform_mode": platform_mode,
            }
        )

    d_id = 1

    # =========================================================================
    # C1: Homozygous Wild-Type (25 designs: 15 dev, 5 val, 5 test)
    # =========================================================================
    c1_lengths = [
        (30, "dev"),
        (40, "dev"),
        (50, "dev"),
        (60, "dev"),
        (70, "dev"),
        (80, "dev"),
        (90, "dev"),
        (100, "dev"),
        (25, "dev"),
        (35, "dev"),
        (45, "dev"),
        (55, "dev"),
        (65, "dev"),
        (75, "dev"),
        (85, "dev"),
        (30, "val"),
        (50, "val"),
        (60, "val"),
        (80, "val"),
        (100, "val"),
        (40, "test"),
        (55, "test"),
        (70, "test"),
        (90, "test"),
        (110, "test"),
    ]
    for length, split in c1_lengths:
        add(
            d_id,
            f"c1_homo_{length}_{length}_{split}",
            split,
            "C1_HOMOZYGOUS_WT",
            [length, length],
            None,
            None,
            structure_file=struct_paths.get(length),
        )
        d_id += 1

    # =========================================================================
    # C2: Asymmetric Length (45 designs: 27 dev, 9 val, 9 test)
    # =========================================================================
    c2_pairs = [
        # 27 dev
        (25, 60),
        (25, 80),
        (25, 100),
        (25, 120),
        (25, 140),
        (30, 70),
        (30, 90),
        (30, 110),
        (30, 130),
        (35, 75),
        (35, 95),
        (35, 115),
        (40, 80),
        (40, 100),
        (40, 120),
        (45, 85),
        (45, 105),
        (45, 125),
        (50, 90),
        (50, 110),
        (50, 130),
        (25, 90),
        (30, 100),
        (35, 110),
        (40, 115),
        (45, 120),
        (50, 125),
        # 9 val
        (25, 70),
        (30, 80),
        (35, 90),
        (40, 95),
        (45, 110),
        (50, 115),
        (25, 110),
        (30, 120),
        (35, 130),
        # 9 test
        (25, 85),
        (30, 95),
        (35, 105),
        (40, 105),
        (45, 115),
        (50, 120),
        (25, 130),
        (30, 140),
        (40, 130),
    ]
    splits_c2 = ["dev"] * 27 + ["val"] * 9 + ["test"] * 9
    for (l1, l2), split in zip(c2_pairs, splits_c2, strict=True):
        add(d_id, f"c2_asym_{l1}_{l2}_{split}", split, "C2_ASYMMETRIC_LENGTH", [l1, l2], None, None)
        d_id += 1

    # =========================================================================
    # C3: Near-Equal Length (40 designs: 24 dev, 8 val, 8 test)
    # =========================================================================
    c3_pairs = [
        # 24 dev (delta 1, 2, 3)
        (40, 41),
        (50, 51),
        (60, 61),
        (70, 71),
        (80, 81),
        (90, 91),
        (100, 101),
        (110, 111),
        (40, 42),
        (50, 52),
        (60, 62),
        (70, 72),
        (80, 82),
        (90, 92),
        (100, 102),
        (110, 112),
        (40, 43),
        (50, 53),
        (60, 63),
        (70, 73),
        (80, 83),
        (90, 93),
        (100, 103),
        (110, 113),
        # 8 val
        (45, 46),
        (55, 56),
        (65, 66),
        (75, 76),
        (45, 47),
        (55, 57),
        (65, 67),
        (75, 77),
        # 8 test
        (35, 36),
        (48, 49),
        (58, 59),
        (68, 69),
        (35, 37),
        (48, 50),
        (58, 60),
        (68, 70),
    ]
    splits_c3 = ["dev"] * 24 + ["val"] * 8 + ["test"] * 8
    for (l1, l2), split in zip(c3_pairs, splits_c3, strict=True):
        add(d_id, f"c3_near_{l1}_{l2}_{split}", split, "C3_NEAR_EQUAL_LENGTH", [l1, l2], None, None)
        d_id += 1

    # =========================================================================
    # C4: Pathogenic 59dupC (50 designs: 30 dev, 10 val, 10 test)
    # =========================================================================
    c4_items = [
        # 30 dev
        ([60, 60], [1, 25], "dev"),
        ([60, 60], [2, 35], "dev"),
        ([50, 50], [1, 20], "dev"),
        ([70, 70], [2, 40], "dev"),
        ([40, 40], [1, 15], "dev"),
        ([80, 80], [2, 50], "dev"),
        ([60, 80], [1, 25], "dev"),
        ([60, 80], [2, 45], "dev"),
        ([50, 70], [1, 20], "dev"),
        ([50, 70], [2, 35], "dev"),
        ([40, 60], [1, 15], "dev"),
        ([40, 60], [2, 30], "dev"),
        ([30, 90], [1, 12], "dev"),
        ([30, 90], [2, 60], "dev"),
        ([25, 100], [1, 10], "dev"),
        ([25, 100], [2, 70], "dev"),
        ([25, 140], [1, 10], "dev"),
        ([25, 140], [2, 90], "dev"),
        ([60, 61], [1, 25], "dev"),
        ([60, 61], [2, 30], "dev"),
        ([50, 51], [1, 20], "dev"),
        ([70, 71], [2, 35], "dev"),
        ([60, 62], [1, 25], "dev"),
        ([60, 62], [2, 35], "dev"),
        ([50, 52], [1, 20], "dev"),
        ([70, 72], [2, 40], "dev"),
        ([60, 63], [1, 25], "dev"),
        ([60, 63], [2, 35], "dev"),
        ([50, 53], [1, 20], "dev"),
        ([70, 73], [2, 40], "dev"),
        # 10 val
        ([60, 60], [1, 30], "val"),
        ([50, 50], [2, 25], "val"),
        ([60, 80], [1, 30], "val"),
        ([50, 70], [2, 40], "val"),
        ([30, 100], [1, 15], "val"),
        ([25, 120], [2, 80], "val"),
        ([60, 61], [1, 28], "val"),
        ([70, 71], [2, 42], "val"),
        ([60, 62], [1, 32], "val"),
        ([50, 52], [2, 28], "val"),
        # 10 test
        ([70, 70], [1, 35], "test"),
        ([80, 80], [2, 45], "test"),
        ([40, 80], [1, 20], "test"),
        ([35, 110], [2, 75], "test"),
        ([25, 130], [1, 12], "test"),
        ([30, 140], [2, 95], "test"),
        ([55, 56], [1, 25], "test"),
        ([65, 66], [2, 38], "test"),
        ([55, 57], [1, 26], "test"),
        ([65, 67], [2, 40], "test"),
    ]
    for lens, tgt, split in c4_items:
        h = f"h{tgt[0]}_r{tgt[1]}"
        add(
            d_id,
            f"c4_dupc_{lens[0]}_{lens[1]}_{h}_{split}",
            split,
            "C4_PATHOGENIC_DUPC",
            lens,
            "dupC",
            tgt,
        )
        d_id += 1

    # =========================================================================
    # C5: Rare Pathogenic Mutations (35 designs: 21 dev, 7 val, 7 test)
    # =========================================================================
    rare_muts = [
        # 21 dev
        ("insG", [60, 60], [1, 20], "dev"),
        ("60dupA", [60, 60], [2, 25], "dev"),
        ("delinsAT", [60, 60], [1, 15], "dev"),
        ("delGCCCA", [60, 60], [2, 5], "dev"),
        ("insCCCC", [60, 60], [1, 30], "dev"),
        ("insG", [60, 80], [1, 25], "dev"),
        ("60dupA", [60, 80], [2, 40], "dev"),
        ("delinsAT", [50, 70], [1, 20], "dev"),
        ("delGCCCA", [50, 70], [2, 6], "dev"),
        ("insCCCC", [40, 60], [1, 18], "dev"),
        ("insG", [30, 90], [2, 50], "dev"),
        ("60dupA", [25, 100], [1, 12], "dev"),
        ("delinsAT", [25, 120], [2, 80], "dev"),
        ("insCCCC", [25, 140], [1, 10], "dev"),
        ("insG", [60, 61], [1, 25], "dev"),
        ("60dupA", [60, 61], [2, 30], "dev"),
        ("delinsAT", [50, 52], [1, 20], "dev"),
        ("delGCCCA", [60, 62], [2, 5], "dev"),
        ("insCCCC", [70, 72], [1, 35], "dev"),
        ("insG", [60, 63], [2, 35], "dev"),
        ("60dupA", [50, 53], [1, 22], "dev"),
        # 7 val
        ("insG", [55, 55], [1, 24], "val"),
        ("60dupA", [65, 65], [2, 32], "val"),
        ("delinsAT", [45, 75], [1, 22], "val"),
        ("insCCCC", [35, 105], [2, 65], "val"),
        ("delGCCCA", [25, 115], [1, 5], "val"),
        ("insG", [55, 56], [2, 28], "val"),
        ("60dupA", [65, 67], [1, 30], "val"),
        # 7 test
        ("insG", [75, 75], [2, 42], "test"),
        ("60dupA", [45, 45], [1, 18], "test"),
        ("delinsAT", [35, 85], [2, 48], "test"),
        ("insCCCC", [25, 125], [1, 10], "test"),
        ("delGCCCA", [30, 135], [2, 6], "test"),
        ("insG", [48, 49], [1, 24], "test"),
        ("60dupA", [58, 60], [2, 32], "test"),
    ]
    for mut, lens, tgt, split in rare_muts:
        h = f"h{tgt[0]}_r{tgt[1]}"
        add(
            d_id,
            f"c5_rare_{mut.lower()}_{lens[0]}_{lens[1]}_{h}_{split}",
            split,
            "C5_RARE_PATHOGENIC",
            lens,
            mut,
            tgt,
        )
        d_id += 1

    # =========================================================================
    # C6: Complex SNVs & Micro-Indels (35 designs: 21 dev, 7 val, 7 test)
    # =========================================================================
    complex_muts = [
        # 21 dev
        ("del18_31", [60, 60], [1, 20], "dev"),
        ("ins16bp", [60, 60], [2, 25], "dev"),
        ("ins25bp", [60, 60], [1, 30], "dev"),
        ("del18_31", [60, 80], [2, 40], "dev"),
        ("ins16bp", [60, 80], [1, 25], "dev"),
        ("ins25bp", [50, 70], [2, 35], "dev"),
        ("del18_31", [50, 70], [1, 22], "dev"),
        ("ins16bp", [40, 60], [2, 30], "dev"),
        ("ins25bp", [40, 60], [1, 18], "dev"),
        ("del18_31", [30, 90], [2, 50], "dev"),
        ("ins16bp", [25, 100], [1, 10], "dev"),
        ("ins25bp", [25, 120], [2, 70], "dev"),
        ("del18_31", [25, 140], [1, 12], "dev"),
        ("ins16bp", [60, 61], [2, 30], "dev"),
        ("ins25bp", [60, 61], [1, 25], "dev"),
        ("del18_31", [50, 52], [2, 28], "dev"),
        ("ins16bp", [60, 62], [1, 25], "dev"),
        ("ins25bp", [70, 72], [2, 40], "dev"),
        ("del18_31", [60, 63], [1, 25], "dev"),
        ("ins16bp", [50, 53], [2, 25], "dev"),
        ("ins25bp", [70, 73], [1, 35], "dev"),
        # 7 val
        ("del18_31", [55, 55], [1, 26], "val"),
        ("ins16bp", [65, 65], [2, 34], "val"),
        ("ins25bp", [45, 75], [1, 22], "val"),
        ("del18_31", [35, 105], [2, 60], "val"),
        ("ins16bp", [25, 115], [1, 10], "val"),
        ("ins25bp", [55, 56], [2, 28], "val"),
        ("del18_31", [65, 67], [1, 32], "val"),
        # 7 test
        ("del18_31", [75, 75], [2, 40], "test"),
        ("ins16bp", [45, 45], [1, 18], "test"),
        ("ins25bp", [35, 85], [2, 45], "test"),
        ("del18_31", [25, 125], [1, 10], "test"),
        ("ins16bp", [30, 135], [2, 80], "test"),
        ("ins25bp", [48, 49], [1, 22], "test"),
        ("del18_31", [58, 60], [2, 30], "test"),
    ]
    for mut, lens, tgt, split in complex_muts:
        h = f"h{tgt[0]}_r{tgt[1]}"
        add(
            d_id,
            f"c6_cmplx_{mut.lower()}_{lens[0]}_{lens[1]}_{h}_{split}",
            split,
            "C6_COMPLEX_SNV_INDEL",
            lens,
            mut,
            tgt,
        )
        d_id += 1

    # =========================================================================
    # C7: Borderline Coverage (20 designs: 12 dev, 4 val, 4 test)
    # =========================================================================
    c7_items = [
        # 12 dev (5x, 8x, 10x, 12x, 15x, 20x)
        ([60, 60], None, None, 20, "dev"),
        ([60, 60], "dupC", [1, 25], 20, "dev"),
        ([60, 80], None, None, 30, "dev"),
        ([60, 80], "dupC", [2, 45], 30, "dev"),
        ([50, 70], None, None, 40, "dev"),
        ([50, 70], "dupC", [1, 20], 40, "dev"),
        ([40, 60], None, None, 50, "dev"),
        ([40, 60], "dupC", [2, 30], 50, "dev"),
        ([25, 100], None, None, 60, "dev"),
        ([25, 100], "dupC", [2, 70], 60, "dev"),
        ([60, 61], None, None, 60, "dev"),
        ([60, 61], "dupC", [1, 25], 60, "dev"),
        # 4 val
        ([55, 55], None, None, 30, "val"),
        ([55, 55], "dupC", [1, 22], 30, "val"),
        ([45, 75], None, None, 45, "val"),
        ([45, 75], "dupC", [2, 40], 45, "val"),
        # 4 test
        ([65, 65], None, None, 25, "test"),
        ([65, 65], "dupC", [2, 30], 25, "test"),
        ([35, 85], None, None, 40, "test"),
        ([35, 85], "dupC", [1, 18], 40, "test"),
    ]
    for lens, c7_mut, c7_tgt, tpl, split in c7_items:
        m_tag = f"_{c7_mut.lower()}" if c7_mut else "_wt"
        add(
            d_id,
            f"c7_cov{tpl}_{lens[0]}_{lens[1]}{m_tag}_{split}",
            split,
            "C7_BORDERLINE_COVERAGE",
            lens,
            c7_mut,
            c7_tgt,
            templates=tpl,
        )
        d_id += 1

    assert len(designs) == 250, f"Expected 250 designs, got {len(designs)}"
    return designs


def main() -> int:
    """Build and write 250 designs to JSON specification."""
    designs = build_500_designs()
    DESIGN_FILE.parent.mkdir(parents=True, exist_ok=True)
    DESIGN_FILE.write_text(json.dumps(designs, indent=2) + "\n", encoding="utf-8")

    splits = {d["split"] for d in designs}
    categories = {d["category"] for d in designs}
    print(f"Successfully generated {len(designs)} designs:")
    for s in sorted(splits):
        count = sum(1 for d in designs if d["split"] == s)
        print(f"  Split '{s}': {count} designs")
    for c in sorted(categories):
        count = sum(1 for d in designs if d["category"] == c)
        print(f"  Category '{c}': {count} designs")
    print(f"Written to: {DESIGN_FILE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
