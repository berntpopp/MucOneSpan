"""Frozen experiment sample inventory with tokenization and strict split accounting."""

from __future__ import annotations

import json
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

DEFAULT_TOKEN_SEED = 500_2026_0915
PRIMARY_MODES = ("amplicon_hifi", "genomic_ont")
SUPPLEMENTARY_MODES = ("genomic_pacbio",)
MODE_TO_PLATFORM = {
    "amplicon_hifi": "pacbio",
    "genomic_ont": "ont",
    "genomic_pacbio": "pacbio",
}


@dataclass(frozen=True)
class SampleInventoryRecord:
    """Explicit inventory entry for one experimental dataset."""

    sample_id: str
    token: str
    design_id: int
    design_name: str
    split: str
    category: str
    platform_mode: str
    sequencing_platform: str
    lengths: list[int]
    mutation: str | None
    targets: list[list[int]]
    structure_file: str | None
    bio_seed: int
    platform_seed: int
    requested_templates: int
    is_supplementary: bool = False

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def generate_token_map(
    design_names: list[str],
    seed: int = DEFAULT_TOKEN_SEED,
) -> dict[str, str]:
    """Deterministically map design names to unique sample tokens."""
    rng = random.Random(seed)
    tokens = [f"sample_{i:04d}" for i in range(1, len(design_names) + 1)]
    rng.shuffle(tokens)
    return dict(zip(design_names, tokens, strict=True))


def load_designs(designs_path: Path) -> list[dict[str, Any]]:
    """Load and validate the raw 250-design specification, sorted by design_id."""
    with designs_path.open("r", encoding="utf-8") as f:
        designs = json.load(f)
    if not isinstance(designs, list) or len(designs) != 250:
        raise ValueError(f"Expected 250 designs in {designs_path}, found {len(designs)}")
    designs.sort(key=lambda d: int(d["design_id"]))
    return designs


def build_expected_inventory(
    designs_path: Path,
    token_seed: int = DEFAULT_TOKEN_SEED,
    include_supplementary: bool = False,
) -> list[SampleInventoryRecord]:
    """Construct the frozen sample inventory from biological designs.

    The primary benchmark contains exactly 250 designs x 2 modes (amplicon_hifi,
    genomic_ont) = 500 datasets. If include_supplementary is True, genomic_pacbio
    runs are also added with is_supplementary=True.
    """
    designs = load_designs(designs_path)
    design_names = [str(d["name"]) for d in designs]
    token_map = generate_token_map(design_names, seed=token_seed)

    records: list[SampleInventoryRecord] = []
    modes = PRIMARY_MODES + SUPPLEMENTARY_MODES if include_supplementary else PRIMARY_MODES

    for d in designs:
        d_name = str(d["name"])
        token = token_map[d_name]
        bio_seed = int(d["bio_seed"])
        hifi_seed = int(d["hifi_seed"])
        ont_seed = int(d["ont_seed"])

        for mode in modes:
            is_supp = mode in SUPPLEMENTARY_MODES
            platform = MODE_TO_PLATFORM[mode]
            plat_seed = ont_seed if platform == "ont" else hifi_seed
            sample_id = f"{d_name}_{mode}"

            rec = SampleInventoryRecord(
                sample_id=sample_id,
                token=token,
                design_id=int(d["design_id"]),
                design_name=d_name,
                split=str(d["split"]),
                category=str(d["category"]),
                platform_mode=mode,
                sequencing_platform=platform,
                lengths=list(d["lengths"]),
                mutation=d.get("mutation"),
                targets=[list(t) for t in d.get("targets", [])],
                structure_file=d.get("structure_file"),
                bio_seed=bio_seed,
                platform_seed=plat_seed,
                requested_templates=int(d.get("requested_templates", 200)),
                is_supplementary=is_supp,
            )
            records.append(rec)

    validate_inventory_integrity(records)
    return records


def validate_inventory_integrity(records: list[SampleInventoryRecord]) -> None:
    """Verify inventory invariants: uniqueness, grouping, split counts, and denominators."""
    primary_records = [r for r in records if not r.is_supplementary]
    if len(primary_records) != 500:
        raise ValueError(
            f"Primary inventory must contain exactly 500 records, got {len(primary_records)}"
        )

    sample_ids = [r.sample_id for r in records]
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Duplicate sample_id detected in inventory")

    # Grouping check: all records for a given design must share the same split
    by_design: dict[str, set[str]] = {}
    for r in records:
        by_design.setdefault(r.design_name, set()).add(r.split)
    for d_name, splits in by_design.items():
        if len(splits) != 1:
            raise ValueError(f"Design {d_name} leaked across multiple splits: {splits}")

    # Split breakdown of primary 500 datasets
    split_counts: dict[str, int] = {}
    for r in primary_records:
        split_counts[r.split] = split_counts.get(r.split, 0) + 1

    expected_splits = {"dev": 300, "val": 100, "test": 100}
    if split_counts != expected_splits:
        raise ValueError(f"Expected primary split counts {expected_splits}, got {split_counts}")
