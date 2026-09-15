# Milestone M0: Specification, Architecture, and Acceptance Protocol
**Project:** MucOneSpan Scientific Reliability and Clinical Reporting  
**Date:** 2026-09-15  
**Authors:** Principal Computational Genomicist & Senior Staff Systems Software Engineer; Principal Clinical Bio-UX Designer & Accessibility Architect  

---

## 1. Environment & Preflight Inventory

- **Host & Kernel:** Linux megalodon 7.0.0-31-generic x86_64 GNU/Linux (Ubuntu)
- **Compute Resources:** 32 physical CPU cores, 59 GiB RAM (45 GiB available), 1.9 TiB available NVMe disk space on `/`
- **Python Runtime:** Python 3.12.9 in `.venv`, locked by `uv.lock` with Python 3.10+ compatibility enforced
- **Bioinformatics Tools:**
  - Clair3: `/home/bernt-popp/miniforge3/envs/env_clair3/bin/run_clair3.sh` (using internal Python 3.9 conda environment)
  - Models: `/home/bernt-popp/miniforge3/envs/env_clair3/bin/models/hifi`, `/home/bernt-popp/miniforge3/envs/env_clair3/bin/models/ont`
  - samtools 1.15.1, minimap2 2.28-r1209, bcftools 1.17 in `/home/bernt-popp/miniforge3/envs/env_clair3/bin`
  - MucOneUp: `/home/bernt-popp/miniforge3/bin/muconeup` (v0.44.5)
  - NanoSim: `/home/bernt-popp/miniforge3/envs/env_nanosim/bin/simulator.py` with pre-trained model at `/home/bernt-popp/development/MucOneUp/reference/nanosim/human_giab_hg002_sub1M_kitv14_dorado_v3.2.1`
  - PacBio Simulator: `/home/bernt-popp/miniforge3/envs/env_pacbio/bin/pbsim`, `ccs`
- **Reporting & UI Verification:**
  - Playwright: `/home/bernt-popp/miniforge3/bin/python` with Playwright 1.45.0 and Chromium 145.0.7632.6 verified functional headless
- **Adversarial Reviewer:**
  - Claude CLI: `/home/bernt-popp/.local/bin/claude` verified non-interactive invocation (`claude -p`), responding as `Claude Fable 5.1` (model ID `claude-fable-5-1`)
- **Baseline Quality Gates:**
  - `make dev`: Resolved and synced 76 packages cleanly
  - `make quality`: 131 source/config files verified (all strictly <= 649 lines), Ruff clean, mypy clean, actionlint clean
  - `make test-fast`: 739 passed in 5.06s
  - `make ci-check`: 739 passed in 11.63s, branch-aware aggregate coverage: **91.34%** (exceeds 80% gate)
  - Pre-existing failures: 0

---

## 2. Dataset Matrix Resolution & Inventory Definition

### 2.1 Scope Ambiguity Resolution
- **Biological Truth Manifest:** `examples/experiment_500_designs.json` specifies 250 biological diploid designs across 7 challenge categories:
  - C1_HOMOZYGOUS_WT: 25 designs
  - C2_ASYMMETRIC_LENGTH: 45 designs
  - C3_NEAR_EQUAL_LENGTH: 40 designs
  - C4_PATHOGENIC_DUPC: 50 designs
  - C5_RARE_PATHOGENIC: 35 designs
  - C6_COMPLEX_SNV_INDEL: 35 designs
  - C7_BORDERLINE_COVERAGE: 20 designs
  - Total designs = 250.
- **Split Distribution:**
  - `dev`: 150 designs (60%)
  - `val`: 50 designs (20%)
  - `test`: 50 designs (20%)
- **Primary Benchmark Definition:**
  - Exactly 250 designs × 2 primary acquisition modes (`amplicon_hifi`, `genomic_ont`) = **500 primary datasets**.
  - Split counts for primary benchmark:
    - `dev`: 150 × 2 = 300 datasets
    - `val`: 50 × 2 = 100 datasets
    - `test`: 50 × 2 = 100 datasets
    - **Denominator:** Exactly 500 datasets. Every dataset has a deterministic token (`sample_0001` .. `sample_0500`) and entry in the frozen expected sample inventory.
- **Supplementary Panel Definition:**
  - `genomic_pacbio` simulated across development challenge categories (counted separately, never mixed into the 500-sample denominator).
  - Note: A full 3-mode factorial would be 250 × 3 = 750 datasets.

### 2.2 Leakage & Provenance Guarantees
- Single biological truth generated per design (using deterministic `bio_seed`).
- Platform read simulation uses isolated seeds (`hifi_seed`, `ont_seed`).
- Biological design is the atomic grouping unit: both acquisition modes of a design reside in the exact same split (`dev`, `val`, or `test`).
- Held-out test split (100 datasets) remains sealed and untouched during algorithm tuning.

---

## 3. Architecture Exploration: Alpha, Beta, Gamma

| Dimension | Option Alpha (Conservative Stabilization) | Option Beta (Evidence-Gated Hybrid Reconstruction) | Option Gamma (Repeat-Unit Graph Reconstruction) |
|---|---|---|---|
| **Core Architecture** | Keep minimap2 + allele clustering + haploid remapping + Clair3 + bcftools consensus. Strengthen filtering, adaptive thresholds, and error handling. | Retain Alpha pipeline for clear separation; gate into targeted k-mer / read-phased local graph resolution for hard cases (C3 near-equal, C2 severe asymmetry, ONT homopolymer smear). | Construct global or repeat-level de Bruijn / variation graph across reads; solve path traversal with diploid flow constraints. |
| **Scientific Assumptions** | Length clustering alone separates alleles when $\Delta L \ge 2$; Clair3 + haploid majority call corrects SNVs; unphased identical lengths collapse to 1 call. | Length clustering handles easy cases ($\Delta L \ge 3$); informative repeat k-mers and read linkage provide phase and disambiguation for $\Delta L \in \{0, 1, 2\}$. | Read length and k-mer spectrum accurately constrain graph traversal without cycle explosion in high-error ONT reads. |
| **CPU / RAM Cost** | Low: minimal overhead above baseline (estimated 5–10s / sample, <= 4 threads, < 2GB RAM). | Moderate: additional k-mer / local linkage analysis on difficult subsets only (estimated 8–15s / sample, < 4GB RAM). | High: full graph construction and constraint solving across all reads (estimated 30–60s / sample, high memory spike on ONT). |
| **Dependencies** | Existing locked dependencies (pysam, numpy, scipy, click, jinja2). No new C/Python packages. | Existing locked dependencies. Self-contained pure Python graph/k-mer algorithms. | Likely requires external graph tools or complex solvers (e.g. networkx / ILP / graphviz). |
| **Compatibility & Regressions** | Negligible risk to existing CLI, output schema, or public APIs. | Low risk: gated activation ensures clean fallback to validated baseline. | High risk: shifts underlying coordinate, alignment, and consensus semantics. |
| **Testability & Determinism** | Fully testable with deterministic unit mocks and synthetic FASTQ fixtures. | Highly testable: unit tests can isolate the gating predicate, k-mer counter, and graph traversal. | Complex: hard to isolate edge cases in graph pruning and non-deterministic heuristics. |
| **Primary Failure Modes** | Cannot resolve near-equal lengths ($\Delta L = 1$) or severe coverage imbalance in ONT where peak merges. | Risk of k-mer noise in low-depth ONT (C7) leading to spurious branching or false phase calls. | Graph tangling in repetitive 60-bp units with noisy ONT; catastrophic phase flips; combinatorial explosion. |

### Recommendation
**Adopt Option Alpha first** for concurrency, durable ledger state, inventory tracking, report architecture, and offline IGV (Milestone M1).  
**Evaluate Option Beta** for scientific deconvolution in Milestone M2 (focusing on C3 near-equal length gating and ONT adaptive clustering).  
**Defer Option Gamma** unless Beta fails empirical validation on the development split.

---

## 4. Systems Design: Concurrency, Durable State, and Evaluation Runner

### 4.1 Orchestrator & Worker Architecture
- **Coordinator:** Thread-based concurrent orchestration (`ThreadPoolExecutor` or worker queue) managing external subprocesses (`muconeup`, `run_clair3.sh`, `minimap2`, `samtools`, `bcftools`).
- **Resource Budgeting:**
  - Target: `--parallel-samples 8 --threads 4` (nominal 32 threads).
  - Subprocess containment: Set `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1` on child processes to prevent thread explosion.
- **Truth-Before-Reads Dependency:**
  - Coordinator generates `truth.fa` once per design.
  - Platform read simulations (`amplicon_hifi`, `genomic_ont`) run only after verified completion of truth generation.

### 4.2 Durable State & Crash-Consistent Ledgers
- Single authoritative writer: Workers return immutable structured result dicts to the main thread via thread-safe queue.
- Atomic commit pattern:
  - Writes to `.tmp.jsonl` on the same filesystem.
  - Flush + `os.fsync`.
  - Atomic rename (`replace`) onto authoritative state and exported ledgers (`ledger_sealed.jsonl`, `ledger_public.jsonl`).
- Idempotent Resume:
  - Checks design name, platform, input/output file existence, and non-empty size/hash before skipping.
  - Incomplete or corrupted runs are detected and cleaned prior to reprocessing.

### 4.3 Evaluation Runner Alignment
- Refactor `scripts/run_eval_pipeline.py`:
  - Enforce explicit `--ledger` and `--data-dir` (remove hardcoded `experiment_200` defaults).
  - Support full manifest splits: `dev`, `val`, `test`, `all`, `pilot`.
  - Ingest the frozen 500-sample inventory: missing or errored samples are marked as failure / no-call and remain strictly in the denominator.
  - Strict stage causal diagnosis: `execution_failed`, `insufficient_evidence`, `length_inference`, `consensus`, `variant_calling`, `classification`, `none`.

---

## 5. Clinical UX Design: Report, Accessibility, and Offline IGV

### 5.1 Clinical Decision Hierarchy
- **Header States:**
  - `PATHOGENIC`: Supported pathogenic variant (e.g. `c.530dupC` frameshift or confirmed pathogenic repeat insertion/deletion) with complete sequence and VCF evidence.
  - `BENIGN`: Negative finding with verified adequate coverage (>= 20x usable per allele) and full diploid reconstruction; clearly states testing scope limitation.
  - `INCONCLUSIVE`: Insufficient depth (< 20x), missing allele, ambiguous reconstruction, or unclassified variant.
- Distinct technical status indicator (e.g. `Execution: Complete`, `QC: Pass`).

### 5.2 Accessibility & Print Architecture
- **WCAG 2.1 AAA Target:** Contrast ratio >= 7:1 for normal body text, >= 4.5:1 for large headings across all themes and states.
- Monospace tabular numbers (`font-variant-numeric: tabular-nums lining-nums`).
- Color + icon/text redundant encoding (never rely on color alone).
- Clean print stylesheet (`@media print`):
  - Expand collapsed details.
  - Page-break controls (`break-inside: avoid` on cards/tables).
  - High-contrast monochrome readability.

### 5.3 Offline Sandboxed IGV
- Complete local asset bundling (no CDN, no remote fonts, no external API calls).
- Embed `igv.min.js`, local FASTA reference contig, and BAM/VCF data URIs.
- Explicit `loadDefaultGenomes: false` and strict Content Security Policy (`default-src 'self' 'unsafe-inline' data: blob:`).

---

## 6. Acceptance Protocol & Milestone Roadmap

### 6.1 Acceptance Metrics
- Primary 500-sample accounting: 100% accounted for in `ledger_sealed.jsonl` and `evaluation_report.json`.
- Zero unhandled exceptions or data corruption on worker crash/resume.
- Zero external network requests during Playwright offline report tests.
- Zero console errors in browser stress tests across viewports 320px to 3840px.
- Aggregate test suite maintains >= 80% branch coverage with all files < 650 lines.

### 6.2 Milestone Schedule
- **M0:** Baseline validation, design specification, dataset matrix freeze, Claude Fable 5.1 review.
- **M1:** Parallel scheduler, atomic ledgers, evaluation runner fixes, pilot verification across 3 modes.
- **M2:** Scientific deconvolution (C3 near-equal lengths, C2 asymmetry, ONT smear, HGVS 20.05) on DEV split.
- **M3:** Clinical report redesign, WCAG AAA compliance, offline IGV packaging, Playwright stress suite.
- **M4:** Final candidate selection on VAL split, frozen evaluation on sealed TEST split (500 datasets complete).
- **M5:** Adversarial review, documentation, changelog, final verification diff.
