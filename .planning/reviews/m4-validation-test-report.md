# Milestone M4/M5: Comprehensive 500-Dataset Benchmark & Delivery Report

**Date:** 2026-09-15  
**Version:** 0.13.0  
**Repository:** MucOneSpan  
**Branch:** `feat/scientific-reliability-reporting`  
**Ledger:** `tests/data/experiment_500/ledger_sealed.jsonl` (501 entries: 1 pilot + 300 DEV + 100 VAL + 100 TEST)

---

## 1. Executive Summary

The MucOneSpan Scientific Reliability and Clinical Reporting initiative was evaluated across a 500-dataset stratified simulation cohort (250 distinct biological diploid designs across PacBio HiFi amplicons and genomic Oxford Nanopore). The dataset partitioning followed strict pre-registration discipline:
- **DEV (300 datasets):** Optimization, calibration, ablation studies, and adversarial review.
- **VAL (100 datasets):** Unsealed to verify fixed pipeline candidates without tuning.
- **TEST (100 datasets):** Held-out final generalization evaluation unsealed only after validation gate confirmation.

### Benchmark Highlights Across All 500 Datasets

| Cohort / Split | Platform | N | Diploid Seq Exact | Diploid Count Exact | Count $\pm 1$ | Supp. Ev. Precision | Supp. Alarm Precision | Normal Specificity |
|---|---|---|---|---|---|---|---|---|
| **DEV** | Amplicon HiFi | 150 | 42 (28.0%) | 96 (64.0%) | 111 (74.0%) | 95.2% (40/42) | 97.6% (41/42) | 95.0% (19/20) |
| | Genomic ONT | 150 | 0 (0.0%) | 1 (0.7%) | 1 (0.7%) | 66.7% (6/9) | 100.0% (10/10) | 100.0% (18/18) |
| | **DEV Pooled** | 300 | 42 (14.0%) | 97 (32.3%) | 112 (37.3%) | 90.2% (46/51) | 98.1% (51/52) | 97.4% (37/38) |
| **VAL** | Amplicon HiFi | 50 | 15 (30.0%) | 34 (68.0%) | 38 (76.0%) | 100.0% (11/11) | 100.0% (11/11) | 100.0% (6/6) |
| | Genomic ONT | 50 | 1 (2.0%) | 4 (8.0%) | 7 (14.0%) | 100.0% (6/6) | 100.0% (6/6) | 100.0% (7/7) |
| | **VAL Pooled** | 100 | 16 (16.0%) | 38 (38.0%) | 45 (45.0%) | 100.0% (17/17) | 100.0% (17/17) | 100.0% (13/13) |
| **TEST** | Amplicon HiFi | 50 | 11 (22.0%) | 36 (72.0%) | 37 (74.0%) | 91.7% (11/12) | 100.0% (12/12) | 100.0% (4/4) |
| | Genomic ONT | 50 | 1 (2.0%) | 1 (2.0%) | 3 (6.0%) | 75.0% (3/4) | 100.0% (4/4) | 100.0% (4/4) |
| | **TEST Pooled**| 100 | 12 (12.0%) | 37 (37.0%) | 40 (40.0%) | 87.5% (14/16) | 100.0% (16/16) | 100.0% (8/8) |
| **COMBINED**| Amplicon HiFi | 250 | 68 (27.2%) | 166 (66.4%) | 186 (74.4%) | 95.4% (62/65) | 98.5% (64/65) | 96.7% (29/30) |
| | Genomic ONT | 250 | 2 (0.8%) | 6 (2.4%) | 11 (4.4%) | 78.9% (15/19) | 100.0% (20/20) | 100.0% (29/29) |
| | **TOTAL** | **500** | **70 (14.0%)** | **172 (34.4%)**| **197 (39.4%)**| **91.7% (77/84)** | **98.8% (84/85)** | **98.3% (58/59)** |

---

## 2. Key Scientific Findings & Failure Atlas

### A. Elimination of Bimodal +1 Shift (Nearest-Peak Clustering)
In Milestone M2, adversarial audit revealed that integer division midpoint splitting `(c1+c2)//2` assigned intermediate contigs to the lower sub-cluster, causing systematic $+1$ repeat inflation for $\Delta \ge 2$ designs (e.g. 40/42 became 41/42).
Replacing this with nearest-peak assignment (`abs(c - c1) < abs(c - c2)`), equidistant contig exclusion, and peak-anchored selection eliminated this bias:
- HiFi diploid count exactness on held-out TEST reached **72.0% (36/50)**, up from 64.0% on DEV and 68.0% on VAL.
- HiFi count accuracy within $\pm 1$ repeat remained stable at **74.0% – 76.0%** across all splits.

### B. Supported Normal Specificity & Clinical Alarm Precision
- **Normal Control Specificity:** 100.0% on VAL (13/13 true negatives) and 100.0% on TEST (8/8 true negatives). Zero false alarms on normal controls across 21 independent held-out controls.
- **Sample Alarm Precision:** 100.0% on VAL (17/17 true alarms) and 100.0% on TEST (16/16 true alarms).
- **Test Event False Positives (N=2):**
  1. `c5_rare_insg_48_49_h1_r24_test_amplicon_hifi`: Truth design was 48/49 with insG at repeat 24 on H1. Prediction recovered lengths 49 and 48 with insG at repeat 24. Under strict bipartite minimum-distance alignment, predicted allele 1 (len 49) was paired with truth H2 (len 49, WT), causing a cross-haplotype pairing mismatch under repeat-exact matching despite the sample alarm correctly firing (true mutation-positive alarm).
  2. `c6_cmplx_del18_31_58_60_h2_r30_test_genomic_ont`: Complex multiallelic deletion/duplication with non-spanning ONT reads. Sample alarm was correctly positive (true positive at sample level).

### C. Genomic ONT Read-Length Profile
Genomic ONT evaluation across 250 datasets confirms that while caller specificity is pristine (100% normal specificity, 100% alarm precision), diploid count and sequence reconstruction in genomic contexts are physically constrained by the 7.5 kb read-length limit in NanoSim. Without targeted amplicon enrichment, the fraction of molecules spanning the entire VNTR plus unique flanks is insufficient for complete assembly. This is an assay/data limitation, not an algorithmic regression.

---

## 3. Clinical Bio-UX & Accessibility Delivery

### Verified Features
1. **3-State Clinical Decision Banner:**
   - Visual states: `PATHOGENIC` (crimson, alert icon), `NO_PATHOGENIC_VARIANT_DETECTED` (navy, check icon), and `INCONCLUSIVE` (amber, info icon).
   - High-contrast color palette adhering to WCAG 2.1 AA/AAA standards (contrast ratio $\ge 7:1$ against backgrounds).
   - Prominent allele multiplicity caveat warning when single observed lengths are detected in diploid individuals.
2. **HGVS 20.05 Nomenclature Compliance:**
   - Unanchored transcript-level cDNA predictions are safely identified as `"transcript_coordinate_unresolved"` and hidden from report summaries.
   - `repeat_{idx}:c.{edit}` serves as the permanent, unambiguous repeat-relative coordinate.
   - Canonical 55delinsAT net-shift naming verified.
3. **Automated Offline Browser Verification:**
   - Test suite: `tests/browser/test_report_browser.py` executed with Playwright.
   - Verified across headless Chromium: zero external HTTP requests (fully offline/air-gapped), zero console error logs, interactive accessibility focus rings, and proper table formatting.

---

## 4. Verification Gates

```bash
make dev          # Reproduce locked environment: SUCCESS
make quality      # Ruff lint/format, mypy typecheck, file-size gates: SUCCESS
make test-fast    # 760 unit tests: 760 passed in 23.36s (SUCCESS)
make ci-check     # Quality + branch-aware coverage: 90.82% coverage (threshold >=80%): SUCCESS
```

- All authored files strictly satisfy the `< 650` physical line requirement.
- Subprocess execution governed by process-group SIGTERM/SIGKILL timeouts (`run_tool` / `run_tool_iter`).
- Cryptographic ledger integrity verified across all 501 entries.
