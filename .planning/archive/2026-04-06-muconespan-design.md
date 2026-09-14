# MucOneSpan Design Specification

**Date:** 2026-04-06
**Status:** Approved
**Author:** Bernt Popp + Claude

## Overview

Open-source reconstruction of the workflow described by Vrbacka et al. (2025, bioRxiv: 10.1101/2025.09.06.673538) for analyzing PacBio HiFi sequencing data of PCR amplicons spanning the MUC1 VNTR region. Source code for the workflow described by Vrbacka et al. was not published; this project reconstructs the pipeline from the methods description.

## Approach

Thin wrapper architecture: each pipeline stage is a small Python module that validates inputs, shells out to bioinformatics tools, and parses outputs. The Python layer handles orchestration, configuration, and pure-Python stages (peak detection, repeat classification). Functional style -- modules expose functions, no class hierarchies.

## Documented Deviations from Vrbacka et al.

| Aspect | Published method (Vrbacka et al.) | MucOneSpan | Rationale |
|--------|-------------------|--------------|-----------|
| Read mapper | bwa-mem 0.7.16a | minimap2 (map-hifi preset) | Faster, better HiFi support, standard for PacBio HiFi |
| Ladder range | 1-120 contigs | 1-150 contigs | Cover potential outliers beyond typical 20-125 repeat range |

## Project Structure

```
MucOneSpan/
├── src/muc_one_span/
│   ├── __init__.py
│   ├── version.py
│   ├── cli.py                 # Click CLI with subcommands
│   ├── ladder.py              # Reference ladder (re)generation
│   ├── mapping.py             # minimap2 mapping + samtools
│   ├── alleles.py             # Peak detection / allele length determination
│   ├── calling.py             # Clair3 variant calling + bcftools filtering
│   ├── consensus.py           # bcftools consensus generation
│   ├── classify.py            # Repeat unit classification + mutation detection
│   ├── config.py              # Configuration loading
│   └── tools.py               # Subprocess helpers for external tools
├── data/
│   ├── reference/             # Pre-built ladder FASTA + index (bundled)
│   └── repeats/               # Repeat definitions (from MucOneUp config.json + muc1repeats)
├── tests/
│   ├── conftest.py
│   ├── unit/                  # Pure Python tests (no external tools needed)
│   └── integration/           # Tests requiring conda env tools
├── conda/
│   └── environment.yml        # Single env: minimap2, samtools, bcftools, clair3
├── docker/
│   └── Dockerfile
├── .github/workflows/
│   └── test.yml
├── pyproject.toml             # hatchling build, ruff, mypy, pytest config
├── Makefile
├── .pre-commit-config.yaml    # ruff + mypy + bandit + standard hooks
└── .gitignore
```

## Tooling

Mirrors the MucOneUp project:

- **Package manager:** uv
- **Build system:** hatchling
- **Linter/formatter:** ruff (lint + format)
- **Type checker:** mypy
- **Tests:** pytest with markers (unit, integration, e2e)
- **Pre-commit:** ruff, mypy, bandit, standard hooks
- **CI:** GitHub Actions (quality checks, unit tests across Python 3.10-3.12, integration tests with conda env)
- **Automation:** Makefile for all common tasks

## CLI Design

```
muconespan <subcommand> [options]
```

### Subcommands

| Command | Purpose | Key inputs | Output |
|---------|---------|------------|--------|
| `ladder` | (Re)generate reference ladder FASTA | `--repeats`, `--min-units`, `--max-units`, `--output` | `reference_ladder.fa` + index |
| `map` | Map reads to ladder reference | `--input` (FASTQ/BAM), `--reference` (defaults to bundled), `--output-dir` | sorted + indexed BAM |
| `alleles` | Determine allele lengths from mapping | `--input` (BAM), `--min-coverage 10` | JSON with allele lengths + confidence |
| `call` | Variant calling with Clair3 | `--input` (BAM), `--reference`, `--alleles` (JSON), `--output-dir` | filtered VCF |
| `consensus` | Build per-allele consensus | `--input` (BAM), `--vcf`, `--alleles` (JSON), `--output-dir` | per-allele FASTA |
| `classify` | Classify repeat units | `--input` (FASTA), `--repeats-db`, `--output-dir` | repeat structure JSON + report |
| `run` | Run full pipeline | `--input`, `--output-dir`, `--config` | all outputs |

- Each subcommand writes outputs to `--output-dir` with deterministic filenames
- `run` chains all subcommands, passing intermediate outputs automatically
- `--config` (YAML) allows overriding all defaults in one file
- Bundled ladder is the default for `map --reference`

## Pipeline Data Flow

```
Input: PacBio HiFi FASTQ/BAM (amplicon reads spanning MUC1 VNTR)
         │
         ▼
┌──────────────────────────────────────────────────────┐
│ 1. Mapping (mapping.py)                              │
│    minimap2 -a -x map-hifi + samtools sort/index     │
│    Input: FASTQ/BAM  Output: sorted BAM vs ladder    │
└──────────────────┬───────────────────────────────────┘
                   │
                   ▼
┌──────────────────────────────────────────────────────┐
│ 2. Allele Detection (alleles.py) -- pure Python      │
│    Parse samtools idxstats -> read counts per contig  │
│    Peak detection -> two allele lengths               │
│    Min 10x coverage filter                            │
│    Output: alleles.json                               │
└──────────────────┬───────────────────────────────────┘
                   │
                   ▼
┌──────────────────────────────────────────────────────┐
│ 3. Read Filtering + Variant Calling (calling.py)     │
│    Split reads by allele (by mapped contig)           │
│    Clair3 --platform=hifi per allele                  │
│    bcftools norm + filter                             │
│    Output: filtered VCF per allele                    │
└──────────────────┬───────────────────────────────────┘
                   │
                   ▼
┌──────────────────────────────────────────────────────┐
│ 4. Consensus (consensus.py)                          │
│    bcftools consensus per allele                      │
│    Applies all variants to reconstruct true sequence  │
│    Output: per-allele consensus FASTA                 │
└──────────────────┬───────────────────────────────────┘
                   │
                   ▼
┌──────────────────────────────────────────────────────┐
│ 5. Classification (classify.py) -- pure Python       │
│    Split consensus into 60bp windows                  │
│    Exact match against repeat dictionary              │
│    Unknown repeats: closest match + similarity        │
│    Mutation detection from classification             │
│    Output: repeat structure, mutation report           │
└──────────────────┬───────────────────────────────────┘
                   │
                   ▼
Output: alleles.json, repeats.json, repeats.txt, summary.json
```

## Module Details

### mapping.py

Wraps `minimap2 -a -x map-hifi` + `samtools sort` + `samtools index`. If input is BAM, extracts reads with `samtools fastq` first. Output: sorted, indexed BAM against ladder reference.

### alleles.py (pure Python)

Parses `samtools idxstats` output into contig -> read-count mapping. Peak detection finds two local maxima across contigs 1-150. Handles homozygous case (single dominant peak). Minimum coverage filter (default 10x, configurable).

Output:
```json
{
  "allele_1": {"length": 60, "reads": 245},
  "allele_2": {"length": 80, "reads": 189},
  "homozygous": false
}
```

### calling.py

Extracts reads mapping to the two best-matching contigs. Runs Clair3 with PacBio HiFi model per allele. Filters VCF with `bcftools norm` + `bcftools filter`.

### consensus.py

Runs `bcftools consensus` per allele to apply all called variants back to the reference contig, reconstructing the true allele sequence. Output: per-allele consensus FASTA.

### classify.py (pure Python)

The core analysis module. Classification algorithm:

1. **Split** consensus into 60bp windows
2. **Exact match** each window against the repeat dictionary (from MucOneUp config.json, supplemented by muc1repeats if it contains additional types)
3. **Unknown repeats** (no exact match): find closest known type by edit distance, report:
   - Closest match type and edit distance
   - Percent identity
   - Nature of all differences (each substitution, insertion, deletion with position)
   - Whether differences contain indel(s) -- if so, classify as mutation of closest type
   - Whether substitutions only -- classify as variant (close) or novel repeat (distant)
4. **Frameshift check**: for indel-containing unknowns, check if `sum(indel_lengths) % 3 != 0`
5. **Output**: repeat structure string + per-repeat detail + mutation report

No hard edit distance thresholds for distinguishing mutations from novel repeats -- the tool reports all evidence (closest match, distance, nature of differences) and lets the classification speak for itself. Large indels (e.g., 16bp duplication) are correctly handled as mutations of known types.

### tools.py

Single `run_tool(cmd, ...)` function for subprocess execution with stdout/stderr capture and informative errors. `check_tools(["minimap2", "samtools", ...])` verifies tool availability at CLI startup.

### config.py

Loads repeat definitions from MucOneUp's `config.json`. Optionally supplements with additional repeats from muc1repeats repo. Runtime settings via YAML config or CLI flags: `min_coverage`, `ladder_range`, `clair3_model`, tool paths (auto-detected from PATH by default).

## Output Format

Per-sample output directory:
```
results/<sample>/
  mapping.bam              # Sorted, indexed BAM against ladder
  mapping.bam.bai
  alleles.json             # Allele length determination results
  variants.vcf.gz          # Filtered Clair3 VCF
  consensus_allele1.fa     # Per-allele consensus FASTA
  consensus_allele2.fa
  repeats.json             # Full repeat classification + mutation report
  repeats.txt              # Human-readable repeat structure string
  summary.json             # Combined results
```

### repeats.json schema

```json
{
  "allele_1": {
    "length": 60,
    "structure": "1 2 3 4 5 X D E C F X X A B ... V 6' 7 8 9",
    "repeats": [
      {"index": 1, "type": "1", "match": "exact"},
      {"index": 6, "type": "X", "match": "exact"},
      {
        "index": 25,
        "type": "unknown",
        "closest_match": "X",
        "edit_distance": 1,
        "identity_pct": 98.3,
        "differences": [{"pos": 59, "ref": "", "alt": "C", "type": "insertion"}],
        "classification": "mutation",
        "mutation": "59dupC",
        "frameshift": true
      },
      {
        "index": 12,
        "type": "unknown",
        "closest_match": "A",
        "edit_distance": 5,
        "identity_pct": 91.7,
        "differences": [
          {"pos": 8, "ref": "G", "alt": "A", "type": "substitution"},
          {"pos": 15, "ref": "C", "alt": "T", "type": "substitution"}
        ],
        "classification": "novel_repeat",
        "provisional_label": "U1"
      }
    ],
    "mutations_detected": [
      {"type": "59dupC", "repeat_index": 25, "frameshift": true}
    ]
  },
  "allele_2": { "..." : "..." }
}
```

## Reference Ladder

Pre-built and bundled in `data/reference/`. 150 contigs (1-150 canonical X repeats), each structured as:

```
[flanking_5' ~500bp] [pre-repeats 1-5] [N * canonical X repeat] [after-repeats 6-9] [flanking_3' ~500bp]
```

- Pre-repeat and after-repeat sequences: canonical versions only (1,2,3,4,5 and 6,7,8,9) from the Vrbacka nomenclature. Variants (4', 4'', 5C, 6', 7') appear in real alleles and are identified at classification time, not in the ladder.
- Flanking sequences from hg38 MUC1 locus
- `ladder` subcommand exists for custom regeneration with different parameters

## Repeat Dictionary

**Primary source:** MucOneUp's `config.json` -- contains exact 60bp DNA sequences and protein sequences for all known repeat types (pre-repeats 1-5 with variants, canonical repeats X/A/B/C/D/E/F/G/I/J/N/R/V and lowercase variants, after-repeats 6-9 with variants).

**Supplement:** muc1repeats repo (github.com/pristanna/muc1repeats) -- checked for any additional repeat types not already in MucOneUp's config. If new types are found, they are added to `data/repeats/`.

## External Tool Dependencies

Single conda environment (`conda/environment.yml`):

| Tool | Version | Purpose |
|------|---------|---------|
| minimap2 | >=2.28 | Read mapping (map-hifi preset) |
| samtools | >=1.21 | BAM processing, idxstats, sort, index |
| Clair3 | >=1.0.10 | Neural network variant calling (PacBio HiFi model) |
| bcftools | >=1.10 | VCF processing, consensus generation |

## Testing

### Test structure

```
tests/
  conftest.py              # Shared fixtures, tool availability checks
  unit/
    test_alleles.py         # Peak detection logic
    test_classify.py        # Repeat classification, exact matching, unknown handling
    test_config.py          # Config loading
    test_tools.py           # Tool availability checking (mocked)
  integration/
    test_mapping.py         # minimap2 + samtools pipeline
    test_calling.py         # Clair3 + bcftools pipeline
    test_consensus.py       # bcftools consensus
    test_pipeline.py        # Full end-to-end with MucOneUp test data
```

### Pytest markers

- `@pytest.mark.unit` -- pure Python, no external tools
- `@pytest.mark.integration` -- requires conda env tools
- `@pytest.mark.e2e` -- full pipeline with MucOneUp-generated test data

### CI workflow (GitHub Actions)

- **Quality job:** ruff lint + format check, mypy (Python 3.10)
- **Unit test job:** pytest unit tests across Python 3.10-3.12, coverage
- **Integration test job:** conda env setup, integration + e2e tests with generated test data

### Test data (10 samples, generated on-the-fly)

Generated via `make generate-testdata` using MucOneUp with fixed seeds. PacBio HiFi amplicon mode, 200x coverage. Deterministic and reproducible.

| # | Sample | Allele lengths | Mutation | Purpose |
|---|--------|---------------|----------|---------|
| 1 | `sample_dupc_60_80` | 60, 80 | 59dupC on hap1 | Common mutation, basic case |
| 2 | `sample_dupa_60_80` | 60, 80 | 60dupA on hap1 | Second most common mutation |
| 3 | `sample_insg_60_80` | 60, 80 | 58_59insG on hap1 | Insertion mutation |
| 4 | `sample_dupcccc_60_80` | 60, 80 | 56_59dupCCCC on hap1 | Multi-base duplication |
| 5 | `sample_del_60_80` | 60, 80 | 18_31delGGCCCCGGACACCA on hap1 | Large deletion |
| 6 | `sample_normal_60_80` | 60, 80 | none | Negative control |
| 7 | `sample_homozygous_60_60` | 60, 60 | 59dupC on hap1 | Homozygous length |
| 8 | `sample_asymmetric_25_140` | 25, 140 | 59dupC on hap1 | Extreme length asymmetry |
| 9 | `sample_short_25_30` | 25, 30 | 59dupC on hap1 | Short alleles |
| 10 | `sample_long_120_140` | 120, 140 | 59dupC on hap1 | Long alleles near upper bound |

Ground truth from MucOneUp's `simulation_stats.json` per sample.

### Makefile targets

```makefile
make init              # Install uv + dev dependencies
make dev               # Install package with dev dependencies
make test              # All tests with coverage
make test-fast         # Unit tests only, no coverage
make test-unit         # Unit tests only with coverage
make test-int          # Integration tests only
make lint              # Run ruff linter
make lint-fix          # Run ruff linter with auto-fix
make format            # Format code with ruff
make format-check      # Check formatting
make type-check        # Run mypy
make check             # lint + format-check + type-check + test
make ci-check          # Exact same checks as GitHub Actions
make clean             # Remove build artifacts and caches
make generate-testdata # Generate test data via MucOneUp (requires conda env)
```
