"""Unit tests for frozen experiment inventory and split accounting."""

from pathlib import Path

import pytest

from muc_one_span.inventory import (
    DEFAULT_TOKEN_SEED,
    SampleInventoryRecord,
    build_expected_inventory,
    generate_token_map,
    load_designs,
    validate_inventory_integrity,
)

DESIGNS_PATH = Path("examples/experiment_500_designs.json")


def test_load_designs_loads_exact_250_entries() -> None:
    designs = load_designs(DESIGNS_PATH)
    assert len(designs) == 250
    assert all("name" in d and "split" in d for d in designs)


def test_token_map_reproducibility() -> None:
    designs = load_designs(DESIGNS_PATH)
    names = [str(d["name"]) for d in designs]
    tokens = generate_token_map(names, seed=DEFAULT_TOKEN_SEED)

    assert len(tokens) == 250
    # Spot-check known tokens from existing 500 experiment runs
    assert tokens["c1_homo_30_30_dev"] == "sample_0183"
    assert tokens["c1_homo_40_40_dev"] == "sample_0053"
    assert tokens["c1_homo_50_50_dev"] == "sample_0196"


def test_primary_inventory_has_exact_500_samples() -> None:
    records = build_expected_inventory(DESIGNS_PATH, include_supplementary=False)
    assert len(records) == 500
    assert all(not r.is_supplementary for r in records)

    validate_inventory_integrity(records)

    dev_recs = [r for r in records if r.split == "dev"]
    val_recs = [r for r in records if r.split == "val"]
    test_recs = [r for r in records if r.split == "test"]

    assert len(dev_recs) == 300
    assert len(val_recs) == 100
    assert len(test_recs) == 100


def test_supplementary_inventory_accounting() -> None:
    records = build_expected_inventory(DESIGNS_PATH, include_supplementary=True)
    assert len(records) == 750

    primary = [r for r in records if not r.is_supplementary]
    supp = [r for r in records if r.is_supplementary]

    assert len(primary) == 500
    assert len(supp) == 250
    assert all(r.platform_mode == "genomic_pacbio" for r in supp)

    # Invariants still hold for the primary records
    validate_inventory_integrity(records)


def test_validate_inventory_integrity_detects_leakage() -> None:
    records = build_expected_inventory(DESIGNS_PATH, include_supplementary=False)
    # Corrupt one record to simulate split leakage
    corrupted = list(records)
    r0 = corrupted[0]
    corrupted[0] = SampleInventoryRecord(
        sample_id=r0.sample_id,
        token=r0.token,
        design_id=r0.design_id,
        design_name=r0.design_name,
        split="test",  # Leaked from dev to test
        category=r0.category,
        platform_mode=r0.platform_mode,
        sequencing_platform=r0.sequencing_platform,
        lengths=r0.lengths,
        mutation=r0.mutation,
        targets=r0.targets,
        structure_file=r0.structure_file,
        bio_seed=r0.bio_seed,
        platform_seed=r0.platform_seed,
        requested_templates=r0.requested_templates,
        is_supplementary=r0.is_supplementary,
    )
    with pytest.raises(ValueError, match="leaked across multiple splits"):
        validate_inventory_integrity(corrupted)


def test_validate_inventory_integrity_detects_count_mismatch() -> None:
    records = build_expected_inventory(DESIGNS_PATH, include_supplementary=False)
    with pytest.raises(ValueError, match="Primary inventory must contain exactly 500"):
        validate_inventory_integrity(records[:-1])
