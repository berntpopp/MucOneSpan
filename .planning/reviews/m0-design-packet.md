Act as the Fable 5.1 Red Team. Review the supplied design for falsifiable genomics, systems, and clinical UX failure modes. Return evidence, counterexamples, and required tests. Do not modify files or execute commands.

# Review Packet: Milestone M0 Design Specification

Review the design document below. Specifically critique:
1. Genomics: identifiability, phase evidence, coordinate errors, false positives, truth leakage across splits, and misleading metrics.
2. Systems: race conditions, partial writes, stale artifacts, resource oversubscription, portability, error propagation, and crash-resilience.
3. UX & Reporting: false reassurance, hidden uncertainty, accessibility failures (WCAG 2.1 AAA), print loss, script injection, and network escape.

Format your findings with severity (High/Medium/Low), precise evidence, counterexamples, reproducible test requirements, and proposed remedies. Distinguish confirmed observations from speculative risks.

---

## Attached Design:

```markdown
# Milestone M0: Specification, Architecture, and Acceptance Protocol
Project: MucOneSpan Scientific Reliability and Clinical Reporting
Date: 2026-09-15

## 1. Environment & Baseline
- 32 CPU cores, 59 GiB RAM (45 GiB available), 1.9 TiB disk.
- Clair3, samtools 1.15.1, minimap2 2.28, bcftools 1.17 in /home/bernt-popp/miniforge3/envs/env_clair3/bin.
- Muconeup v0.44.5, NanoSim 3.2.3, pbsim3 3.0.5, ccs 6.4.0.
- Playwright Chromium 145.0.7632.6 verified.
- Baseline test suite: 739 passed in 11.6s, 91.34% branch coverage, 131 files <= 649 lines.

## 2. Dataset Scope
- 250 biological designs from examples/experiment_500_designs.json:
  - Splits: dev=150 (60%), val=50 (20%), test=50 (20%).
  - Challenge categories: C1 (homozygous WT), C2 (asymmetric 25-140), C3 (near-equal 40-71), C4 (pathogenic dupC), C5 (rare pathogenic), C6 (complex indel/SNV), C7 (borderline 10-20x coverage).
- Primary benchmark: 250 designs x 2 acquisition modes (amplicon_hifi, genomic_ont) = exactly 500 datasets.
- Supplementary panel: genomic_pacbio for dev challenge categories (separately reported).
- Denominator integrity: All 500 datasets are in the expected inventory; missing or failed samples remain in the denominator.
- Leakage protection: Grouping is by biological design; all modes of a design share the same split. Truth generated once per design.

## 3. Architecture Options: Alpha, Beta, Gamma
- Option Alpha (Conservative stabilization): Retain minimap2 + allele clustering + haploid remapping + Clair3 + bcftools consensus. Strengthen filtering, adaptive thresholds, error handling.
- Option Beta (Evidence-gated hybrid reconstruction): Gated activation of k-mer and read-phased local graph for hard cases (C3 near-equal, C2 severe asymmetry, ONT smear), while retaining Alpha for clear cases.
- Option Gamma (Repeat-unit graph reconstruction): Global de Bruijn / variation graph with diploid flow constraints across all reads.
- Decision: Adopt Alpha for M1 (infrastructure/ledgers/runner/report), evaluate Beta for M2 (scientific deconvolution), defer Gamma.

## 4. Systems Concurrency & Durability
- Coordinator: ThreadPoolExecutor managing subprocesses with thread caps (OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1). Budget: 8 parallel samples x 4 threads.
- Durable state: Workers return structured result dicts to main thread. Main thread commits atomically via write to temp file on same filesystem + fsync + atomic rename (replace).
- Evaluation runner: Refactor scripts/run_eval_pipeline.py to accept explicit ledger and data root, support manifest split names (dev, val, test, all, pilot), ingest frozen inventory to preserve failed runs in denominator, and output causal failure diagnosis.

## 5. Clinical Decision & Offline Reporting
- Decision header: PATHOGENIC (supported finding + evidence), BENIGN (negative with >=20x depth and complete reconstruction), INCONCLUSIVE (insufficient depth, missing allele, ambiguous call).
- Accessibility: WCAG 2.1 AAA target (>= 7:1 normal text, >= 4.5:1 large text), tabular lining nums, non-color redundant encoding.
- Print CSS: uncollapse evidence, break-inside: avoid, clean monochrome table layout.
- Offline IGV: fully bundled local assets (no CDN/remote calls), loadDefaultGenomes: false, CSP header.

## 6. Verification & Acceptance
- Zero network requests in headless Playwright tests.
- Zero console errors across 320px to 3840px viewports.
- All files strictly < 650 lines. Aggregated branch coverage >= 80%.
```
