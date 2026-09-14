# MucOneSpan Implementation Plan

**Project:** Open-source reconstruction of the workflow described by Vrbacka et al. (2025)
**Date:** 2026-04-06
**Status:** Planning

## Background

The workflow described by Vrbacka et al. (2025) is a bioinformatics pipeline (bioRxiv: 10.1101/2025.09.06.673538) for analyzing PacBio SMRT (CCS/HiFi) sequencing data of PCR amplicons spanning the MUC1 VNTR. The pipeline determines allele lengths, constructs consensus sequences, identifies repeat unit types, and detects frameshift mutations.

**Source code for the workflow described by Vrbacka et al. was not published.** This project reconstructs the pipeline from the methods description as an open-source, containerized tool.

## Pipeline Architecture

```
Input: PacBio HiFi FASTQ/BAM (amplicon reads spanning MUC1 VNTR)
                │
                ▼
┌─────────────────────────────────┐
│ 1. Reference Ladder Generation  │
│    Synthetic contigs with       │
│    20-150 canonical repeats     │
│    (pre-repeats + X*N + after)  │
└──────────────┬──────────────────┘
               │
               ▼
┌─────────────────────────────────┐
│ 2. Read Mapping                 │
│    bwa-mem to ladder reference  │
│    samtools sort + index        │
└──────────────┬──────────────────┘
               │
               ▼
┌─────────────────────────────────┐
│ 3. Allele Length Determination  │
│    samtools idxstats per contig │
│    Peak detection (2 alleles)   │
│    Min 10x coverage filter      │
└──────────────┬──────────────────┘
               │
               ▼
┌─────────────────────────────────┐
│ 4. Read Filtering               │
│    Remove chimeric/incomplete   │
│    Quality filtering            │
└──────────────┬──────────────────┘
               │
               ▼
┌─────────────────────────────────┐
│ 5. Variant Calling              │
│    Clair3 (PacBio-trained)      │
│    bcftools processing          │
└──────────────┬──────────────────┘
               │
               ▼
┌─────────────────────────────────┐
│ 6. Consensus Construction       │
│    bcftools consensus from VCF  │
│    Per-allele consensus FASTA   │
└──────────────┬──────────────────┘
               │
               ▼
┌─────────────────────────────────┐
│ 7. Repeat Unit Classification   │
│    Split into 60bp units        │
│    Classify using nomenclature  │
│    (leverage muc1repeats repo)  │
│    Statistics + mutation report  │
└──────────────┬──────────────────┘
               │
               ▼
Output: Allele lengths, mutations detected, repeat structure, consensus FASTA
```

## Implementation Phases

### Phase 1: Reference Ladder Generator

**Goal:** Build the synthetic multi-contig reference FASTA.

**Details:**
- Each contig represents a possible VNTR allele length (20 to 150 repeat units)
- Contig structure: `[flanking_5'] [pre-repeats 1-5] [X * N] [after-repeats 6-9] [flanking_3']`
- Pre-repeat and after-repeat sequences from MucOneUp config or Vrbacka nomenclature
- Canonical X repeat = 60bp consensus sequence
- Flanking sequences from hg38 MUC1 locus (enough for read anchoring, ~500bp each side)
- Output: `reference_ladder.fa` with 131 contigs (20-150 inclusive)
- Also generate bwa index

**Inputs needed:**
- Canonical repeat sequences (X, 1-9) from MucOneUp config.json or Vrbacka
- hg38 MUC1 flanking coordinates

### Phase 2: Read Mapping and Allele Length Detection

**Goal:** Map amplicon reads and determine the two allele lengths.

**Details:**
- Map with `bwa-mem` (parameters tuned for long reads: `-x pacbio` or custom)
- `samtools sort`, `samtools index`
- `samtools idxstats` to count reads per contig
- Peak detection: find the two contigs with the most mapped reads
- Handle edge cases: homozygous (single peak), ambiguous peaks, low coverage
- Minimum coverage threshold: 10x (from Vrbacka methods)
- Output: allele length estimates + confidence metrics

### Phase 3: Read Filtering and Variant Calling

**Goal:** Call frameshift variants within the VNTR.

**Details:**
- Filter chimeric reads (supplementary alignments), incomplete mappings
- Run Clair3 with PacBio HiFi model on reads mapped to best-matching contigs
- Process VCF with bcftools (normalize, filter)
- Identify frameshift-causing variants (insertions/deletions not divisible by 3)
- Output: filtered VCF with candidate mutations

### Phase 4: Consensus and Repeat Classification

**Goal:** Build per-allele consensus and classify repeat units.

**Details:**
- `bcftools consensus` to generate per-allele FASTA
- Split consensus into 60bp repeat units
- Classify each unit against known repeat types:
  - Use repeat dictionary from https://github.com/pristanna/muc1repeats
  - Vrbacka nomenclature: pre-repeats (1-5), canonical (A-Z, aA-aZ...), after-repeats (6-9)
- Detect mutations: compare each repeat to reference, identify frameshifts
- Output: repeat structure string, mutation report, statistics

### Phase 5: CLI and Packaging

**Goal:** Usable command-line tool with Docker support.

**Details:**
- Click-based CLI: `muconespan run --input reads.bam --output results/`
- Configuration via YAML or JSON
- Docker/Singularity container with all dependencies
- Conda environment.yml
- Comprehensive logging and error handling

### Phase 6: Validation

**Goal:** Verify reconstruction accuracy.

**Details:**
- Test on MucOneUp-simulated PacBio HiFi amplicon data (see TESTING_WITH_MUCONEUP.md)
- Compare results to expected ground truth
- If possible, test on any publicly available real PacBio MUC1 data
- Document concordance with Vrbacka's reported results

## Tool Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| bwa-mem | 0.7.16a+ | Read mapping to ladder reference |
| samtools | 1.9+ | BAM processing, idxstats, sort, index |
| Clair3 | 1.0.10+ | Neural network variant calling (PacBio model) |
| bcftools | 1.10+ | VCF processing, consensus generation |
| Python | 3.10+ | CLI, repeat classification, analysis |
| Click | 8.0+ | CLI framework |

## Repository Structure

```
MucOneSpan/
├── .planning/
│   ├── PLAN.md                    # This file
│   └── TESTING_WITH_MUCONEUP.md   # Test data generation guide
├── src/
│   └── muc_one_span/
│       ├── __init__.py
│       ├── cli.py                 # Click CLI
│       ├── reference.py           # Ladder reference generation
│       ├── mapping.py             # bwa-mem mapping + allele detection
│       ├── calling.py             # Clair3 variant calling
│       ├── consensus.py           # bcftools consensus
│       ├── classify.py            # Repeat unit classification
│       └── config.py              # Configuration management
├── data/
│   ├── repeats/                   # Repeat type definitions
│   └── reference/                 # Generated ladder reference
├── tests/
│   └── ...
├── docker/
│   └── Dockerfile
├── conda/
│   └── environment.yml
├── pyproject.toml
├── README.md
└── LICENSE                        # MIT
```

## Success Criteria

1. Correctly determines allele lengths from simulated PacBio HiFi amplicon data
2. Detects known frameshift mutations (dupC at minimum) with >90% sensitivity
3. Repeat unit classification matches Vrbacka nomenclature
4. Runs in Docker container with no manual setup
5. Processes a sample in <5 minutes on standard hardware
6. Documented and reproducible

## Timeline Estimate

| Phase | Effort | Dependencies |
|-------|--------|--------------|
| Phase 1: Reference ladder | 1-2 days | MUC1 repeat sequences |
| Phase 2: Mapping + allele detection | 2-3 days | Phase 1 |
| Phase 3: Variant calling | 1-2 days | Phase 2 |
| Phase 4: Repeat classification | 2-3 days | Phase 3, muc1repeats repo |
| Phase 5: CLI + Docker | 1-2 days | Phase 4 |
| Phase 6: Validation | 2-3 days | Phase 5, MucOneUp test data |
| **Total** | **~10-15 days** | |

## References

- Vrbacka et al. 2025: bioRxiv 10.1101/2025.09.06.673538
- muc1repeats: https://github.com/pristanna/muc1repeats
- Clair3: https://github.com/HKU-BAL/Clair3
- Wenzel et al. 2018: Scientific Reports 8:4928 (original PacBio MUC1 approach)
- Kirby et al. 2013: Nature Genetics 45:299-303 (MUC1 ADTKD discovery)
