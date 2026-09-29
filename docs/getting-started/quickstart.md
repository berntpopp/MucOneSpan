# Quick Start

Get started with MucOneSpan in under 5 minutes. This tutorial walks through the full analysis pipeline.

---

## Prerequisites

- MucOneSpan installed ([Installation Guide](installation.md))
- PacBio HiFi CCS reads or Oxford Nanopore (ONT) Q20+ reads from a MUC1 VNTR PCR amplicon
  (genomic long reads also work; pass `--assay genomic` to record the library type)
- `samtools` on PATH for BAM input; the deprecated ladder engine also needs minimap2,
  bcftools and Clair3

---

## Full Pipeline (Recommended)

Run the default hybrid engine in a single command:

```bash
muconespan run \
  --input reads.fastq \
  --output-dir results/ \
  --report
```

The same command serves HiFi and ONT reads: the hybrid engine takes no platform,
thread, coverage or reference option. Those are ladder-only; given to a hybrid run
they print a warning and are listed in `summary.json["ignored_options"]`. BAM input
is streamed read by read through `samtools fastq`, so subset a WGS BAM to the MUC1
region first (for example `samtools view -b in.bam chr1:155185000-155195000` on hg38,
which covers the MUC1 gene; adjust the coordinates for other builds).

**What happens** (see [Core Concepts](concepts.md#hybrid-engine)):

1. Anchors every read on the conserved motifs and sorts spanning, partial and off-target reads
2. Fits the allele lengths from the spanning reads (smear and PCR-dimer aware)
3. Builds a partial-order-alignment consensus per allele and splits same-length alleles on linked sites
4. Assigns every read, polishes the consensus and classifies each 60bp repeat unit
5. Scores each detected mutation against the reads assigned to its allele

### Ladder engine (deprecated)

The previous five-stage pipeline stays available with `--engine ladder` until a
later release removes it. It prints a deprecation warning and records it in
`summary.json["deprecations"]`:

```bash
muconespan run \
  --engine ladder \
  --input reads.fastq \
  --output-dir results/ \
  --clair3-model /path/to/clair3/models/hifi \
  --threads 8
```

1. Generates a synthetic reference ladder (150 contigs, 1-150 repeat units)
2. Maps reads to the ladder with minimap2
3. Detects two allele lengths from the read count distribution
4. Calls variants per allele with Clair3
5. Builds consensus and classifies each 60bp repeat unit

---

## Step-by-Step Execution

For more control, run each ladder stage individually:

### 1. Generate Reference Ladder

```bash
muconespan ladder --output reference_ladder.fa
```

This creates a FASTA with 150 contigs, each containing 1-150 canonical X repeats plus flanking sequences.

### 2. Map Reads

```bash
muconespan map \
  --input reads.fastq \
  --reference reference_ladder.fa \
  --output-dir results/ \
  --threads 8
```

### 3. Detect Alleles

```bash
muconespan alleles \
  --input results/mapping.bam \
  --output-dir results/
```

### 4. Call Variants

```bash
muconespan call \
  --input results/mapping.bam \
  --reference reference_ladder.fa \
  --alleles-json results/alleles.json \
  --output-dir results/ \
  --clair3-model /path/to/clair3/models/hifi
```

### 5. Build Consensus and Classify

```bash
muconespan consensus \
  --input results/mapping.bam \
  --reference reference_ladder.fa \
  --alleles-json results/alleles.json \
  --output-dir results/

muconespan classify \
  --input results/consensus_allele_1.fa \
  --output-dir results/
```

---

## ONT Data

The default hybrid engine needs no platform option for ONT reads. With the
deprecated ladder engine, add `--platform ont` to `run` or individual subcommands;
the pipeline then auto-selects `minimap2 -x lr:hq` and `Clair3 --platform=ont`:

```bash
muconespan run \
  --engine ladder \
  --input ont_reads.fastq \
  --output-dir results/ \
  --platform ont \
  --threads 8
```

For step-by-step execution, pass `--platform ont` to `map` and `call`:

```bash
muconespan map --input ont_reads.fastq --reference ref.fa --output-dir results/ --platform ont
muconespan call --input results/mapping.bam --reference ref.fa --alleles-json results/alleles.json --output-dir results/ --platform ont
```

You can also override the minimap2 preset explicitly with `--minimap2-preset`:

```bash
muconespan run --input reads.fastq --output-dir results/ --minimap2-preset map-ont
```

---

## Input Format

- **PacBio HiFi CCS reads** or **ONT Q20+ reads** as FASTQ or BAM
- Reads should be from a PCR amplicon spanning the MUC1 VNTR (primers per Wenzel et al. 2018)
- Minimum recommended coverage: **10x per allele** (higher coverage improves confidence)

!!! tip "Coverage recommendation"
    For clinical-grade results, aim for 50-100x per allele. With the default hybrid engine
    an allele needs `hybrid.depth_adequate_spanning` (default 30) spanning reads for an
    adequate depth, and a length peak needs `hybrid.min_peak_reads` reads; see the
    [configuration guide](../guides/configuration.md#hybrid-engine). The `--min-coverage`
    flag (default: 10) applies only to the deprecated ladder engine's allele detection.

---

## Output Files

| File | Description |
|------|-------------|
| `alleles.json` | Detected allele lengths, read counts, contig assignments |
| `allele_1/variants.vcf.gz` | Filtered Clair3 variants for allele 1 |
| `allele_2/variants.vcf.gz` | Filtered Clair3 variants for allele 2 |
| `consensus_allele_1.fa` | Consensus FASTA (VNTR region only) |
| `consensus_allele_2.fa` | Consensus FASTA (VNTR region only) |
| `repeats.json` | Per-repeat classification with confidence scores |
| `repeats.txt` | Human-readable VNTR structure string |
| `summary.json` | Combined allele + classification + mutation report |

### Example Output

```
results/
├── alleles.json
├── allele_1/
│   └── variants.vcf.gz
├── allele_2/
│   └── variants.vcf.gz
├── consensus_allele_1.fa
├── consensus_allele_2.fa
├── repeats.json
├── repeats.txt
└── summary.json
```

### Reading Results

```bash
# View allele lengths
cat results/alleles.json | python -m json.tool

# View VNTR structure
cat results/repeats.txt
# allele_1: 1 2 3 4 5 C F X X X X:dupC B A A B X X X V 6 7 8 9
# allele_2: 1 2 3 4 5 C F X X X X X B A A B X X X V 6 7 8 9

# View confidence scores
python -c "import json; d=json.load(open('results/repeats.json')); print(d['allele_1']['allele_confidence'])"
```

---

## Next Steps

- **[Core Concepts](concepts.md)** -- Understand the pipeline architecture
- **[Differences from the Published Method](deviations.md)** -- What changed and why
- **[CLI Reference](../reference/cli.md)** -- All command options
- **[Benchmarking](../guides/benchmarking.md)** -- Validate with simulated data
