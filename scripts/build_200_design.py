#!/usr/bin/env python3
"""Build and validate the 100-design specification for the 200-dataset experiment."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
STRUCTURES_DIR = ROOT / "examples" / "structures"
DESIGN_FILE = ROOT / "examples" / "experiment_100_designs.json"


def generate_identical_structure_content(total_length: int) -> str:
    """Generate a valid MucOneUp structure string for identical haplotypes."""
    # Fixed pre repeats: 1-2-3-4-5
    # Canonical variable repeats: X units (total_length - 9)
    # Fixed after repeats: 6-7-8-9
    canonical_count = total_length - 9
    if canonical_count < 1:
        raise ValueError(f"Total length {total_length} must be >= 10")
    variable_part = "-".join(["X"] * canonical_count)
    chain = f"1-2-3-4-5-{variable_part}-6-7-8-9"
    return f"haplotype_1\t{chain}\nhaplotype_2\t{chain}\n"


def create_authored_structures() -> dict[int, str]:
    """Create structure files for identical sequence designs."""
    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    lengths = [30, 40, 50, 60, 70, 80, 100]
    paths = {}
    for length in lengths:
        filename = f"identical_{length}_{length}.txt"
        file_path = STRUCTURES_DIR / filename
        content = generate_identical_structure_content(length)
        file_path.write_text(content)
        paths[length] = str(file_path.relative_to(ROOT))
    return paths


def build_100_designs() -> list[dict[str, Any]]:
    """Build the 100 biological designs (70 dev, 30 final-validation)."""
    identical_paths = create_authored_structures()
    designs: list[dict[str, Any]] = []

    def add_design(
        d_id: int,
        name: str,
        split: str,
        category: str,
        lengths: list[int],
        mutation: str | None,
        target: list[int] | None,
        templates: int,
        structure_file: str | None = None,
        expected_minority: float = 0.5,
    ) -> None:
        bio_seed = 1_000_000 + d_id * 10
        hifi_seed = 2_000_000 + d_id * 10
        ont_seed = 3_000_000 + d_id * 10
        targets = [target] if target else []
        designs.append(
            {
                "design_id": d_id,
                "name": name,
                "split": split,
                "category": category,
                "lengths": lengths,
                "mutation": mutation,
                "targets": targets,
                "requested_templates": templates,
                "structure_file": structure_file,
                "bio_seed": bio_seed,
                "hifi_seed": hifi_seed,
                "ont_seed": ont_seed,
                "expected_minority": expected_minority,
            }
        )

    # -------------------------------------------------------------
    # 1. DEVELOPMENT SPLIT (70 designs: 22 controls, 48 mutation positive)
    # -------------------------------------------------------------
    d_id = 1

    # Dev Controls - Equal Identical (7 designs)
    for length in [30, 40, 50, 60, 70, 80, 100]:
        add_design(
            d_id,
            f"ctrl_ident_{length}_{length}",
            "dev",
            "control_equal_identical",
            [length, length],
            None,
            None,
            200,
            structure_file=identical_paths[length],
        )
        d_id += 1

    # Dev Controls - Equal Variant (5 designs)
    for length, tpl in [(40, 200), (50, 200), (60, 200), (60, 60), (80, 200)]:
        suffix = "_lowcov" if tpl == 60 else ""
        add_design(
            d_id,
            f"ctrl_eqvar_{length}_{length}{suffix}",
            "dev",
            "control_equal_variant",
            [length, length],
            None,
            None,
            tpl,
        )
        d_id += 1

    # Dev Controls - Heterozygous Gaps (10 designs)
    gap_controls_dev = [
        ("ctrl_gap1_60_61", [60, 61], 200),
        ("ctrl_gap1_40_41", [40, 41], 200),
        ("ctrl_gap2_60_62", [60, 62], 200),
        ("ctrl_gap2_50_52", [50, 52], 200),
        ("ctrl_gap3_60_63", [60, 63], 200),
        ("ctrl_typ_40_50", [40, 50], 200),
        ("ctrl_typ_60_80", [60, 80], 200),
        ("ctrl_typ_80_100", [80, 100], 200),
        ("ctrl_asym_25_100", [25, 100], 200),
        ("ctrl_lowcov_60_80", [60, 80], 20),
    ]
    for name, lengths, tpl in gap_controls_dev:
        add_design(d_id, name, "dev", "control_heterozygous", lengths, None, None, tpl)
        d_id += 1

    # Verification: Exactly 22 dev controls created
    assert d_id == 23, f"Expected d_id 23 for dev controls, got {d_id}"

    # Dev Mutation Positive (48 designs)
    # Dev: Equal Lengths with Mutation (7 designs)
    eq_muts_dev = [
        ("mut_eq_dupc_60_60_h1", [60, 60], "dupC", [1, 25], 200),
        ("mut_eq_dupc_60_60_h2", [60, 60], "dupC", [2, 25], 200),
        ("mut_eq_dupa_60_60_h1", [60, 60], "dupA", [1, 30], 200),
        ("mut_eq_insg_60_60_h2", [60, 60], "insG", [2, 20], 200),
        ("mut_eq_inscccc_60_60_h1", [60, 60], "insCCCC", [1, 35], 200),
        ("mut_eq_del18_31_60_60_h2", [60, 60], "del18_31", [2, 25], 200),
        ("mut_eq_delinsat_60_60_h1", [60, 60], "delinsAT", [1, 15], 200),
    ]
    for name, lengths, mut, tgt, tpl in eq_muts_dev:
        add_design(d_id, name, "dev", "equal_length_mutation", lengths, mut, tgt, tpl)
        d_id += 1

    # Dev: Gap 1 Mutation Positive (4 designs)
    gap1_dev = [
        ("mut_gap1_dupc_60_61_h1", [60, 61], "dupC", [1, 25], 200),
        ("mut_gap1_dupa_60_61_h2", [60, 61], "dupA", [2, 25], 200),
        ("mut_gap1_insg_50_51_h1", [50, 51], "insG", [1, 20], 200),
        ("mut_gap1_inscccc_75_76_h2", [75, 76], "insCCCC", [2, 35], 200),
    ]
    for name, lengths, mut, tgt, tpl in gap1_dev:
        add_design(d_id, name, "dev", "gap1_mutation", lengths, mut, tgt, tpl)
        d_id += 1

    # Dev: Gap 2 Mutation Positive (4 designs)
    gap2_dev = [
        ("mut_gap2_dupc_60_62_h1", [60, 62], "dupC", [1, 25], 200),
        ("mut_gap2_dupa_60_62_h2", [60, 62], "dupA", [2, 25], 200),
        ("mut_gap2_insg_40_42_h1", [40, 42], "insG", [1, 18], 200),
        ("mut_gap2_del18_31_80_82_h2", [80, 82], "del18_31", [2, 40], 200),
    ]
    for name, lengths, mut, tgt, tpl in gap2_dev:
        add_design(d_id, name, "dev", "gap2_mutation", lengths, mut, tgt, tpl)
        d_id += 1

    # Dev: Gap 3 Mutation Positive (4 designs)
    gap3_dev = [
        ("mut_gap3_dupc_60_63_h1", [60, 63], "dupC", [1, 25], 200),
        ("mut_gap3_inscccc_60_63_h2", [60, 63], "insCCCC", [2, 30], 200),
        ("mut_gap3_ins16bp_45_48_h1", [45, 48], "ins16bp", [1, 20], 200),
        ("mut_gap3_ins25bp_70_73_h2", [70, 73], "ins25bp", [2, 35], 200),
    ]
    for name, lengths, mut, tgt, tpl in gap3_dev:
        add_design(d_id, name, "dev", "gap3_mutation", lengths, mut, tgt, tpl)
        d_id += 1

    # Dev: Boundary Repeat Mutations (11 designs: terminal canonical 6-8 or L-6..L-4, and delGCCCA on fixed repeat 5/6)
    boundary_dev = [
        ("mut_bnd_dupc_60_80_h1_r6", [60, 80], "dupC", [1, 6], 200),
        ("mut_bnd_dupc_60_80_h2_r74", [60, 80], "dupC", [2, 74], 200),
        ("mut_bnd_dupa_60_80_h1_r7", [60, 80], "dupA", [1, 7], 200),
        ("mut_bnd_dupa_60_80_h2_r75", [60, 80], "dupA", [2, 75], 200),
        ("mut_bnd_insg_50_70_h1_r8", [50, 70], "insG", [1, 8], 200),
        ("mut_bnd_insg_50_70_h2_r65", [50, 70], "insG", [2, 65], 200),
        ("mut_bnd_inscccc_60_80_h1_r54", [60, 80], "insCCCC", [1, 54], 200),
        ("mut_bnd_delgccca_60_80_h1_r5", [60, 80], "delGCCCA", [1, 5], 200),
        ("mut_bnd_delgccca_60_80_h2_r6", [60, 80], "delGCCCA", [2, 6], 200),
        ("mut_bnd_delgccca_40_60_h1_r5", [40, 60], "delGCCCA", [1, 5], 200),
        ("mut_bnd_delgccca_50_70_h2_r6", [50, 70], "delGCCCA", [2, 6], 200),
    ]
    for name, lengths, mut, tgt, tpl in boundary_dev:
        add_design(d_id, name, "dev", "boundary_mutation", lengths, mut, tgt, tpl)
        d_id += 1

    # Dev: Extreme Asymmetry (10 designs with PCR bias minority tiers)
    asym_dev = [
        ("mut_asym_dupc_25_140_h1", [25, 140], "dupC", [1, 15], 200),
        ("mut_asym_dupc_25_140_h2", [25, 140], "dupC", [2, 80], 200),
        ("mut_asym_dupa_25_100_h1", [25, 100], "dupA", [1, 15], 200),
        ("mut_asym_insg_25_100_h2", [25, 100], "insG", [2, 60], 200),
        ("mut_asym_inscccc_30_120_h1", [30, 120], "insCCCC", [1, 20], 200),
        ("mut_asym_del18_31_30_120_h2", [30, 120], "del18_31", [2, 70], 200),
        ("mut_asym_delinsat_20_90_h1", [20, 90], "delinsAT", [1, 12], 200),
        ("mut_asym_ins16bp_20_90_h2", [20, 90], "ins16bp", [2, 50], 200),
        ("mut_asym_tier1_dupc_25_140_low", [25, 140], "dupC", [2, 100], 60),
        ("mut_asym_tier1_dupa_25_140_low", [25, 140], "dupA", [2, 90], 60),
    ]
    for name, lengths, mut, tgt, tpl in asym_dev:
        add_design(d_id, name, "dev", "extreme_asymmetry_mutation", lengths, mut, tgt, tpl)
        d_id += 1

    # Dev: Typical Heterozygous Gaps & Named Mutations Coverage (8 designs)
    # Covering insC_pos23, insG_pos58, insG_pos54, insA_pos54, ins25bp, etc.
    typ_dev = [
        ("mut_typ_insc_pos23_60_80_h1", [60, 80], "insC_pos23", [1, 25], 200),
        ("mut_typ_insc_pos23_60_80_h2", [60, 80], "insC_pos23", [2, 35], 200),
        ("mut_typ_insg_pos58_50_70_h1", [50, 70], "insG_pos58", [1, 20], 200),
        ("mut_typ_insg_pos58_50_70_h2", [50, 70], "insG_pos58", [2, 30], 200),
        ("mut_typ_insg_pos54_60_80_h1", [60, 80], "insG_pos54", [1, 25], 200),
        ("mut_typ_insg_pos54_60_80_h2", [60, 80], "insG_pos54", [2, 40], 200),
        ("mut_typ_insa_pos54_40_60_h1", [40, 60], "insA_pos54", [1, 15], 200),
        ("mut_typ_insa_pos54_40_60_h2", [40, 60], "insA_pos54", [2, 25], 200),
    ]
    for name, lengths, mut, tgt, tpl in typ_dev:
        add_design(d_id, name, "dev", "typical_named_mutation", lengths, mut, tgt, tpl)
        d_id += 1

    # Verification: Exactly 70 dev designs created
    assert d_id == 71, f"Expected d_id 71 after dev, got {d_id}"

    # -------------------------------------------------------------
    # 2. PROTECTED FINAL-VALIDATION SPLIT (30 designs: 10 controls, 20 mutation positive)
    # -------------------------------------------------------------
    # Final Controls - Equal Identical (3 designs)
    for length in [40, 60, 80]:
        add_design(
            d_id,
            f"val_ctrl_ident_{length}_{length}",
            "final",
            "control_equal_identical",
            [length, length],
            None,
            None,
            200,
            structure_file=identical_paths[length],
        )
        d_id += 1

    # Final Controls - Equal Variant (2 designs)
    for length in [50, 70]:
        add_design(
            d_id,
            f"val_ctrl_eqvar_{length}_{length}",
            "final",
            "control_equal_variant",
            [length, length],
            None,
            None,
            200,
        )
        d_id += 1

    # Final Controls - Heterozygous Gaps (5 designs)
    gap_controls_final = [
        ("val_ctrl_gap1_50_51", [50, 51], 200),
        ("val_ctrl_gap2_60_62", [60, 62], 200),
        ("val_ctrl_gap3_70_73", [70, 73], 200),
        ("val_ctrl_typ_60_80", [60, 80], 200),
        ("val_ctrl_asym_25_100", [25, 100], 200),
    ]
    for name, lengths, tpl in gap_controls_final:
        add_design(d_id, name, "final", "control_heterozygous", lengths, None, None, tpl)
        d_id += 1

    # Verification: Exactly 10 final controls created (total 32 controls)
    assert d_id == 81, f"Expected d_id 81 after final controls, got {d_id}"

    # Final Mutation Positive (20 designs)
    # Final: Equal Lengths with Mutation (3 designs)
    eq_muts_val = [
        ("val_mut_eq_dupc_60_60_h1", [60, 60], "dupC", [1, 25], 200),
        ("val_mut_eq_dupa_60_60_h2", [60, 60], "dupA", [2, 25], 200),
        ("val_mut_eq_insg_60_60_h1", [60, 60], "insG", [1, 30], 200),
    ]
    for name, lengths, mut, tgt, tpl in eq_muts_val:
        add_design(d_id, name, "final", "equal_length_mutation", lengths, mut, tgt, tpl)
        d_id += 1

    # Final: Boundary Repeat Mutations (5 designs)
    bnd_muts_val = [
        ("val_mut_bnd_dupc_60_80_h1_r6", [60, 80], "dupC", [1, 6], 200),
        ("val_mut_bnd_dupa_60_80_h2_r74", [60, 80], "dupA", [2, 74], 200),
        ("val_mut_bnd_inscccc_50_70_h1_r7", [50, 70], "insCCCC", [1, 7], 200),
        ("val_mut_bnd_delgccca_60_80_h1_r5", [60, 80], "delGCCCA", [1, 5], 200),
        ("val_mut_bnd_delgccca_60_80_h2_r6", [60, 80], "delGCCCA", [2, 6], 200),
    ]
    for name, lengths, mut, tgt, tpl in bnd_muts_val:
        add_design(d_id, name, "final", "boundary_mutation", lengths, mut, tgt, tpl)
        d_id += 1

    # Final: Gaps 1, 2, 3 Mutation Positive (6 designs)
    gaps_muts_val = [
        ("val_mut_gap1_insg_60_61_h1", [60, 61], "insG", [1, 25], 200),
        ("val_mut_gap1_del18_31_60_61_h2", [60, 61], "del18_31", [2, 25], 200),
        ("val_mut_gap2_dupc_50_52_h1", [50, 52], "dupC", [1, 20], 200),
        ("val_mut_gap2_inscccc_60_62_h2", [60, 62], "insCCCC", [2, 30], 200),
        ("val_mut_gap3_delinsat_60_63_h1", [60, 63], "delinsAT", [1, 25], 200),
        ("val_mut_gap3_ins16bp_60_63_h2", [60, 63], "ins16bp", [2, 30], 200),
    ]
    for name, lengths, mut, tgt, tpl in gaps_muts_val:
        add_design(d_id, name, "final", "gaps_mutation", lengths, mut, tgt, tpl)
        d_id += 1

    # Final: Extreme Asymmetry (5 designs)
    asym_muts_val = [
        ("val_mut_asym_dupc_25_140_h2", [25, 140], "dupC", [2, 80], 200),
        ("val_mut_asym_dupa_25_100_h1", [25, 100], "dupA", [1, 15], 200),
        ("val_mut_asym_insg_30_120_h2", [30, 120], "insG", [2, 70], 200),
        ("val_mut_asym_ins25bp_25_100_h1", [25, 100], "ins25bp", [1, 15], 200),
        ("val_mut_asym_tier1_dupc_25_140_low", [25, 140], "dupC", [2, 100], 60),
    ]
    for name, lengths, mut, tgt, tpl in asym_muts_val:
        add_design(d_id, name, "final", "extreme_asymmetry_mutation", lengths, mut, tgt, tpl)
        d_id += 1

    # Final: Length Extreme (1 design)
    add_design(
        d_id,
        "val_mut_ext_dupc_120_140_h1",
        "final",
        "length_extreme",
        [120, 140],
        "dupC",
        [1, 50],
        200,
    )
    d_id += 1

    assert len(designs) == 100, f"Expected exactly 100 designs, got {len(designs)}"
    return designs


def main() -> None:
    designs = build_100_designs()
    dev_count = sum(1 for d in designs if d["split"] == "dev")
    val_count = sum(1 for d in designs if d["split"] == "final")
    ctrl_count = sum(1 for d in designs if d["mutation"] is None)
    mut_count = sum(1 for d in designs if d["mutation"] is not None)

    print(f"Total designs: {len(designs)}")
    print(f"Development split: {dev_count} (x2 platforms = {dev_count * 2} datasets)")
    print(f"Final-validation split: {val_count} (x2 platforms = {val_count * 2} datasets)")
    print(
        f"Controls: {ctrl_count} (Dev: {sum(1 for d in designs if d['split'] == 'dev' and d['mutation'] is None)}, Final: {sum(1 for d in designs if d['split'] == 'final' and d['mutation'] is None)})"
    )
    print(f"Mutation positive: {mut_count}")

    output_data = {
        "schema_version": 2,
        "description": "100 distinct biological designs for 200-dataset MucOneSpan evaluation",
        "total_designs": len(designs),
        "total_datasets": len(designs) * 2,
        "split_counts": {"dev": dev_count, "final": val_count},
        "designs": designs,
    }
    DESIGN_FILE.write_text(json.dumps(output_data, indent=2) + "\n")
    print(f"Wrote designs to {DESIGN_FILE}")


if __name__ == "__main__":
    main()
