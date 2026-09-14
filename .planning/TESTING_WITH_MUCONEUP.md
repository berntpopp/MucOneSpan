# Testing muconespan with MucOneUp-Simulated Data

**Purpose:** Generate PacBio HiFi amplicon test datasets with known ground truth to validate the muconespan pipeline.

## Prerequisites

- MucOneUp v0.43.2+ installed and configured (https://github.com/berntpopp/muconeup)
- pbsim3 installed (for PacBio HiFi simulation)
- minimap2 installed (for alignment)
- A PacBio HiFi error model file (`.model`) for pbsim3

## Quick Start: Generate a Single Test Sample

### Step 1: Simulate a VNTR Haplotype Pair (Mutant)

```bash
muconeup --config config.json simulate \
  --out-base test_mutant \
  --num-haplotypes 2 \
  --fixed-lengths 60 80 \
  --mutation-name dupC \
  --mutation-targets 1,25 \
  --output-structure \
  --seed 42
```

This generates:
- `test_mutant.001.simulated.fa` (Haplotype 1, 60 repeats, WITH dupC mutation on repeat 25)
- `test_mutant.002.simulated.fa` (Haplotype 2, 80 repeats, normal)
- `test_mutant.simulation_stats.json` (ground truth metadata)

### Step 2: Simulate a Matched Normal Control

```bash
muconeup --config config.json simulate \
  --out-base test_normal \
  --num-haplotypes 2 \
  --fixed-lengths 60 80 \
  --mutation-name normal \
  --output-structure \
  --seed 42
```

Same VNTR structure but no mutation -- serves as negative control.

### Step 3: Generate PacBio HiFi Amplicon Reads

```bash
# Mutant sample
muconeup reads amplicon test_mutant.001.simulated.fa test_mutant.002.simulated.fa \
  --model-file /path/to/pacbio_hifi.model \
  --coverage 200 \
  --seed 42 \
  --platform pacbio

# Normal control
muconeup reads amplicon test_normal.001.simulated.fa test_normal.002.simulated.fa \
  --model-file /path/to/pacbio_hifi.model \
  --coverage 200 \
  --seed 42 \
  --platform pacbio
```

This generates:
- Aligned BAM files (PacBio HiFi reads mapped to hg38)
- PCR bias is modeled (shorter allele gets proportionally more reads)
- Coverage = total template molecules before CCS filtering

### Step 4: Run muconespan

```bash
muconespan run \
  --input test_mutant.aligned.bam \
  --output results/test_mutant/

muconespan run \
  --input test_normal.aligned.bam \
  --output results/test_normal/
```

### Step 5: Validate Against Ground Truth

Compare muconespan output to MucOneUp's `simulation_stats.json`:
- Allele lengths: Should match `--fixed-lengths 60 80`
- Mutation detected: Mutant sample should report dupC; normal should not
- Repeat structure: Should match MucOneUp's output structure

## Comprehensive Test Suite

### Test Set 1: Basic Validation (5 pairs)

Generate 5 matched pairs with varying allele lengths:

```bash
for seed in 100 101 102 103 104; do
  # Mutant
  muconeup --config config.json simulate \
    --out-base pair_${seed}_mut \
    --num-haplotypes 2 \
    --seed ${seed} \
    --mutation-name dupC \
    --mutation-targets 1,25 \
    --output-structure

  # Normal
  muconeup --config config.json simulate \
    --out-base pair_${seed}_normal \
    --num-haplotypes 2 \
    --seed ${seed} \
    --mutation-name normal \
    --output-structure

  # PacBio amplicon reads for both
  muconeup reads amplicon pair_${seed}_mut.*.simulated.fa \
    --model-file /path/to/model \
    --coverage 200 --seed ${seed} --platform pacbio

  muconeup reads amplicon pair_${seed}_normal.*.simulated.fa \
    --model-file /path/to/model \
    --coverage 200 --seed ${seed} --platform pacbio
done
```

**Expected:** 5/5 mutant detected, 0/5 normal false positives, allele lengths within ±2 repeats.

### Test Set 2: Multiple Mutation Types (13 pairs)

One pair per real mutation type:

```bash
mutations=(dupC dupA insG insCCCC insC_pos23 insG_pos58 insG_pos54 insA_pos54 delGCCCA ins25bp ins16bp del18_31 delinsAT)
seed=200

for mut in "${mutations[@]}"; do
  muconeup --config config.json simulate \
    --out-base test_${mut} \
    --num-haplotypes 2 \
    --fixed-lengths 60 80 \
    --mutation-name ${mut} \
    --mutation-targets 1,25 \
    --output-structure \
    --seed ${seed}

  muconeup reads amplicon test_${mut}.*.simulated.fa \
    --model-file /path/to/model \
    --coverage 200 --seed ${seed} --platform pacbio

  seed=$((seed + 1))
done
```

**Expected:** All 13 mutations detected. This validates muconespan handles diverse mutation types.

### Test Set 3: Coverage Sensitivity (1 pair, 5 coverage levels)

```bash
for cov in 500 200 100 50 20; do
  muconeup reads amplicon test_mutant.*.simulated.fa \
    --model-file /path/to/model \
    --coverage ${cov} --seed 42 --platform pacbio
  # Rename output to include coverage level
done
```

**Expected:** Detection at 500x-100x. Possible failure at 20x (near Vrbacka's 10x minimum).

### Test Set 4: Extreme VNTR Lengths (4 pairs)

```bash
# Very short alleles
muconeup --config config.json simulate \
  --out-base test_short --fixed-lengths 25 30 \
  --mutation-name dupC --mutation-targets 1,10 \
  --output-structure --seed 300

# Very long alleles
muconeup --config config.json simulate \
  --out-base test_long --fixed-lengths 120 140 \
  --mutation-name dupC --mutation-targets 1,50 \
  --output-structure --seed 301

# Highly asymmetric
muconeup --config config.json simulate \
  --out-base test_asymmetric --fixed-lengths 25 140 \
  --mutation-name dupC --mutation-targets 1,10 \
  --output-structure --seed 302

# Homozygous length
muconeup --config config.json simulate \
  --out-base test_homozygous --fixed-lengths 60 60 \
  --mutation-name dupC --mutation-targets 1,25 \
  --output-structure --seed 303
```

**Expected:** Tests edge cases for allele length determination algorithm.

### Test Set 5: ONT Amplicon (Cross-Platform)

Same samples as Test Set 1 but with ONT reads. muconespan now supports ONT via `--platform ont`.

**Generate ONT reads:**

```bash
for seed in 100 101 102 103 104; do
  # Mutant
  muconeup reads amplicon pair_${seed}_mut.*.simulated.fa \
    --model-file /path/to/ont_model \
    --coverage 200 --seed ${seed} --platform ont

  # Normal
  muconeup reads amplicon pair_${seed}_normal.*.simulated.fa \
    --model-file /path/to/ont_model \
    --coverage 200 --seed ${seed} --platform ont
done
```

**Run muconespan with ONT mode:**

```bash
for sample_dir in tests/data/generated_ont/*/; do
  bam=$(ls ${sample_dir}*_reads_amplicon_aligned.bam 2>/dev/null | head -1)
  [ -z "$bam" ] && continue
  muconespan run \
    --input "$bam" \
    --output-dir "tests/results_ont/$(basename $sample_dir)" \
    --platform ont \
    --threads 4
done
```

Or use the batch script:

```bash
python scripts/batch_analyze.py tests/data/generated_ont tests/results_ont --platform ont
```

**Minimum ONT test samples** (for `scripts/generate_testdata.py` ONT extension):
- `sample_ont_dupc_60_80` — dupC mutation, ONT reads
- `sample_ont_normal_60_80` — normal control, ONT reads
- `sample_ont_dupa_60_80` — dupA mutation, ONT reads

**Expected:** Validate cross-platform consistency. ONT may show lower sensitivity for 1bp frameshifts due to higher error rate.

## Ground Truth Extraction

MucOneUp's `simulation_stats.json` contains:

```json
{
  "haplotype_1": {
    "total_repeats": 60,
    "mutation_applied": "dupC",
    "mutation_repeat_index": 25
  },
  "haplotype_2": {
    "total_repeats": 80,
    "mutation_applied": null
  }
}
```

Use this to build the ground truth table for sensitivity/specificity calculations.

## Validation Script Template

```python
"""Compare muconespan results to MucOneUp ground truth."""
import json
import csv
from pathlib import Path

def validate_results(results_dir: Path, ground_truth_dir: Path) -> dict:
    """Compare predicted allele lengths and mutations to ground truth."""
    tp = fp = tn = fn = 0
    length_errors = []

    for result_file in results_dir.glob("*/results.json"):
        sample = result_file.parent.name
        gt_file = ground_truth_dir / sample / "simulation_stats.json"

        with open(result_file) as f:
            predicted = json.load(f)
        with open(gt_file) as f:
            truth = json.load(f)

        # Check allele lengths (within tolerance)
        for hap in ["haplotype_1", "haplotype_2"]:
            true_len = truth[hap]["total_repeats"]
            pred_len = predicted.get(f"{hap}_length")
            if pred_len:
                length_errors.append(abs(true_len - pred_len))

        # Check mutation detection
        has_mutation = any(
            truth[h].get("mutation_applied") not in (None, "normal")
            for h in ["haplotype_1", "haplotype_2"]
        )
        detected = predicted.get("mutation_detected", False)

        if has_mutation and detected: tp += 1
        elif has_mutation and not detected: fn += 1
        elif not has_mutation and detected: fp += 1
        else: tn += 1

    return {
        "sensitivity": tp / (tp + fn) if (tp + fn) > 0 else None,
        "specificity": tn / (tn + fp) if (tn + fp) > 0 else None,
        "mean_length_error": sum(length_errors) / len(length_errors) if length_errors else None,
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
    }
```

## Notes

- The `--mutation-targets 1,25` syntax means: mutate haplotype 1, repeat index 25.
  Adjust repeat index based on allele length (must be < total repeats).
- For very short alleles, use a smaller repeat index (e.g., `1,10`).
- The `--platform pacbio` flag in `reads amplicon` generates PacBio HiFi via pbsim3 + CCS.
- The `--platform ont` flag generates ONT reads via single-pass pbsim3.
- Coverage parameter = total template molecules. Actual read count depends on CCS filtering
  (PacBio) or direct output (ONT).
- The PCR bias model means shorter alleles get more reads. This is realistic and muconespan
  should handle the unequal allelic representation.
