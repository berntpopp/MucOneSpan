# MucOneSpan

**Open-source MUC1 VNTR analysis pipeline for PacBio HiFi and ONT amplicon data**

---

## What is MucOneSpan?

MucOneSpan analyzes MUC1 VNTR long-read amplicon data from PacBio HiFi and ONT sequencing. It draws on some ideas from Vrbacka et al. (2025). Features include **indel-valley allele splitting**, **mutation template matching**, and **per-repeat confidence scoring**.

Since 0.17.0 `muconespan run` uses the read-centric **hybrid engine** by default: it reconstructs each allele from the reads and scores each mutation against the reads of its allele, without minimap2 or Clair3. The ladder engine described under Key Features is deprecated (`--engine ladder`); see the [migration guide](guides/migration.md).

Start with the [installation guide](getting-started/installation.md) or [quickstart](getting-started/quickstart.md).

### Why MucOneSpan?

**Detect frameshift mutations** in the MUC1 VNTR that cause ADTKD-MUC1 kidney disease
**Resolve allele pairs** even when they differ by only 3-9 repeat units (indel-valley splitting)
**Score confidence** per repeat unit and per allele with VCF cross-validation
**Classify repeats** using the Vrbacka nomenclature with a pre-computed mutation catalog
**Fully open source** -- MIT licensed, reproducible, no manual inspection steps

---

## Scientific Context

**Autosomal Dominant Tubulointerstitial Kidney Disease caused by MUC1 mutations (ADTKD-MUC1)** is a rare genetic kidney disease caused by frameshift mutations in the Variable Number Tandem Repeat (VNTR) region of the MUC1 gene. The VNTR consists of 20-125 copies of a degenerate 60-bp repeat unit with extremely high GC content (>80%), making it inaccessible to standard short-read sequencing.

The most common mutation, **59dupC**, duplicates a cytosine in the heptanucleotide C-tract of a canonical X repeat unit, producing a +1 frameshift that leads to expression of the toxic MUC1fs protein. **Long-read sequencing** (PacBio SMRT) combined with specialized bioinformatics is required to resolve the full VNTR structure and locate pathogenic frameshift mutations.

---

## Key Features

### Mutation Catalog

Pre-computed sequence templates for **13 known MUC1 frameshift mutations** enable O(1) lookup for common mutations. Novel mutations are detected via edit-distance fallback.

### Indel-Valley Allele Splitting

Resolves close allele pairs (3-9 repeats apart) that simple peak-finding merges into a single cluster. Analyzes CIGAR indel lengths to find two local minima corresponding to the true allele lengths.

### Confidence Scoring

Per-repeat and per-allele confidence scores with **VCF cross-validation** -- flags low-confidence calls for automated QC.

### VCF Validation

Clair3 variant calls are cross-referenced against repeat classifications to adjust confidence scores and detect compound heterozygotes.

---

## Quick Example

```bash
muconespan run \
  --input reads.fastq \
  --output-dir results/
```

**Output:** Per-allele VNTR structure, mutation calls with exact repeat position, and confidence scores.

---

## Documentation

<div class="grid cards" markdown>

-  **[Getting Started](getting-started/installation.md)**
  Installation, quick start tutorial, and core concepts

-  **[Guides](guides/benchmarking.md)**
  Benchmarking with MucOneUp simulated data

-  **[Reference](reference/cli.md)**
  CLI commands, known mutations, and repeat nomenclature

-  **[About](about/citation.md)**
  Citation guide, license, and changelog

</div>

---

## Acknowledgment

MucOneSpan draws on some ideas from [Vrbacka et al. (2025)](https://doi.org/10.1101/2025.09.06.673538).

---

**Development Status:** Active | **License:** MIT | **Maintained by:** [Bernt Popp](https://github.com/berntpopp)
