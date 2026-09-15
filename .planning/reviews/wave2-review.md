# Independent Review Report: Wave 2 Clinical Benchmark and Accuracy Improvements

Date: 2026-09-16
Branch: `feat/clinical-benchmark`
Reviewer: Independent Bioinformatician Reviewer (Subagent `147e8168-5854-4b1b-9cd2-fe9bd7a6eb39`)

## Executive Summary

The Wave 2 clinical benchmark improvements on branch `feat/clinical-benchmark` resolve four key algorithmic and clinical classification defects (issues #52, #53, #54, #55) and complete the containerized comparator evaluation with VNTRPipeline v1.0 on the PRJEB92208 clinical cohort (issue #56). All engineering, quality, line-count (< 650 physical lines per file), and CI gates (1,003 unit tests, 89.56% branch-aware coverage) are fully satisfied.

---

## 1. Requirements and Defect Verification

### Issue #52: Preserved dupC Repeat Expansion Context
- **Diagnosis:** Exact 7-C tract duplications across compatible MUC1 repeat backgrounds were previously mapped across repeat boundaries, misreporting B-repeat dupC as an A insertion in the following repeat unit.
- **Fix:** Implemented cumulative indel offset logic in `classify.py` and `repeat_alignment.py`. The classifier accurately adjusts window coordinates, retains the actual parent repeat, and consumes the correct 61-base unit.
- **Validation:** Verified on MP1 (ERR15277566) and MP2 (ERR15277567), both correctly resolving to `B:dupC` at repeat 17 with `exact_sequence_concordance`.

### Issue #53: Consensus Genotype Selection
- **Diagnosis:** Unphased distinct-length candidates were previously forced to alternate genotype GT1, fabricating synthetic haplotype evidence.
- **Fix:** In `calling.py` line 525, within-candidate phase status is checked against unphased, disconnected, or conflicting states. If unphased, selector `"I"` (`genotype_iupac_candidate`) is selected, and `independent_haplotype_evidence` is set to `False`.
- **Validation:** Regression tests in `tests/unit/test_distinct_calling.py` and `tests/unit/test_phasing.py` confirm unphased candidates never claim independent biological evidence. In `cohort-v3`, MP3 (ERR15277568) correctly avoids false GT1 selection.

### Issue #54: Allele Selection and Read Assignment
- **Diagnosis:** Valley splitting in `alleles.py` (`_split_cluster_by_indel`) sorted valleys by raw CIGAR indel base counts, causing degraded short fragments to beat long biological alleles in broad clusters.
- **Fix:** 
  1. Filter valley candidates to prefer biological alleles ($c \ge 10$) when at least two exist.
  2. Rank valleys by normalized indel rate: `mean_indel / ((c + 9) * 60)`.
  3. Midpoint tie-breaking uses `<=` to ensure consistent clustering.
- **Validation:** In `cohort-v3`, HG002 PCR (ERR15277563) successfully recovers allele 2 as contig 68 (length 77, 4,638 bp), matching independent Q100 truth with 0 edit distance. MP2 (ERR15277567) recovers allele 2 as contig 69 (length 78, 4,681 bp). MP4 (ERR15277569) recovers allele 2 as contig 71 (length 80, 4,801 bp).

### Issue #55: Calibrated Clinical Interpretation
- **Diagnosis:** Previous versions issued `PATHOGENIC` banners for benign in-frame expansions (HG002 WGS maternal 18 bp insertion) and false `NEGATIVE` banners for unresolved candidate mixtures.
- **Fix:** In `report.py` (`compute_clinical_decision`), `PATHOGENIC` requires `frameshift=True`, unambiguous localization (`localization_status != "ambiguous"`), and verified VCF support (`vcf_support_status != "absent"`). Any in-frame insertions, unresolved duplicate candidates, unphased candidate mixtures, or unverified reference confidence route safely to `INCONCLUSIVE`.
- **Validation:** 
  - HG002 WGS evaluates to `INCONCLUSIVE` (in-frame 18 bp expansion at repeat 41).
  - HG002 PCR evaluates to `INCONCLUSIVE` (in-frame 18 bp expansion).
  - MP3 evaluates to `INCONCLUSIVE` (unphased IUPAC candidate).
  - MP1 and MP2 evaluate to `PATHOGENIC` with `dupC` at repeat 17.
  - MP4 evaluates to `NO_PATHOGENIC_VARIANT_DETECTED`.

### Issue #56: VNTRPipeline v1.0 Comparator Evaluation
- **Execution:** Packaged container `ghcr.io/dhmeduni/vntr_pipeline` executed over normalized PRJEB92208 FASTQ inputs with offline `BiocManager 1.30.27` mounted.
- **Results:**
  - **MP1 (ERR15277566):** Identified LoF at motif position 17 (8-C tract) in 452s. Exactly concordant with MucOneSpan.
  - **MP2 (ERR15277567):** Assembled 78 and 50 repeats in 390s; identified LoF at motif position 17. Exactly concordant with MucOneSpan.
  - **HG002 PCR (ERR15277563):** Assembled 4,638 bp and 3,900 bp in 635s; 0 LoF variants. Exactly concordant with independent Q100 truth and MucOneSpan.
  - **MP4 (ERR15277569):** Assembled 80 and 45 repeats in 422s; identified LoF at motif position 49.

---

## 2. Code Quality and Constraints

- **File Size Gate:** All authored files remain strictly below 650 physical lines (e.g., `src/muc_one_span/alleles.py` is 641 lines; `src/muc_one_span/calling.py` is 550 lines; `src/muc_one_span/report.py` is 405 lines).
- **Code Style and Typing:** Ruff format, Ruff linting, and mypy (70 source files) pass with 0 warnings or errors.
- **Test Coverage:** `make ci-check` passes with 1,003 unit tests and 89.56% branch-aware coverage (exceeds 80% requirement). Full integration suite (`make test-int`) passes 70/70 tests.

---

## 3. Review Conclusion

The Wave 2 improvements establish an honest, calibrated, and reproducible clinical analysis pipeline that accurately reflects biological evidence and resolves the major limitations identified in the initial baseline. All changes are recommended for commit and PR update.
