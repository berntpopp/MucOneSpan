# Milestone M2 Scientific & Systems Review Packet

## 1. Executive Summary

Milestone M2 addresses the core scientific reliability challenges of MUC1 VNTR reconstruction across complex repetitive alleles:
1. **Bimodal Read-Length Separation for Near-Equal Lengths (Category C3):** In PacBio HiFi sequencing, individual read lengths have standard deviation $\sigma \approx 1\text{--}2\text{ bp}$. A single 60-bp repeat unit difference represents approximately $30\sigma$ physical separation. We implemented `split_cluster_by_read_length` in `muc_one_span.length_candidates` to separate near-equal length alleles (e.g. 40/41, 50/51, 60/61) deterministically from physical read lengths without reference-ladder collapse.
2. **Haplotagged Allele Separation for Equal-Length Alleles (Finding G2 / Category C4):** In mixed diploid BAMs where both alleles have identical VNTR repeat lengths (e.g. 50/50, 60/60), heterozygous 1-bp homopolymer indels (`dupC`) are suppressed by Clair3 as homozygous reference calls (0/0). We developed `haplotag_and_split_reads` using WhatsHap haplotagging and `samtools view -d HP:N` to partition reads into pure single-haplotype BAMs. Calling Clair3 on each haplotype independently recovers `dupC` with 100% precision and exact sequence concordance.
3. **Haploid Majority Variant Safety (Finding G3):** Updated `src/muc_one_span/vcf.py` to preserve allele frequencies between 0.20 and 0.50 as heterozygous/mixed evidence rather than silently overwriting them to homozygous reference 0/0.
4. **HGVS 20.05 Nomenclature (Finding G2):** Implemented `format_repeat_relative_coordinate` producing invariant repeat-relative clinical coordinates (`repeat_{idx}:c.{edit}`) and marking unresolved cDNA positions with `transcript_coordinate_unresolved`.

---

## 2. Quantitative Evidence: Full 300-Dataset DEV Split Evaluation

The entire 300-dataset DEV split (150 unique designs × 2 primary platforms: `amplicon_hifi` and `genomic_ont`) was simulated and evaluated using the single-truth durable ledger and strict denominator accounting.

### Platform Performance Summary

| Metric | Amplicon HiFi (N=150) | Genomic ONT (N=150) | Total Cohort (N=300) |
|---|---:|---:|---:|
| Correct Exact Reconstructions (`none`) | **72 (48.0%)** | 3 (2.0%) | 75 (25.0%) |
| Correct Length Pairs (`none` + `consensus`) | **114 (76.0%)** | 5 (3.3%) | 119 (39.7%) |
| Consensus / Minor SNV Divergence (`consensus`) | 42 (28.0%) | 2 (1.3%) | 44 (14.7%) |
| Length Inference Discrepancy (`length_inference`) | 36 (24.0%) | 101 (67.3%) | 137 (45.7%) |
| Insufficient Usable Reads (`usable_reads`) | **0 (0.0%)** | 44 (29.3%) | 44 (14.7%) |
| Supported Event True Positives (TP) | **40** | 6 | 46 |
| Supported Event False Positives (FP) | **2** | 3 | 5 |
| Supported Event Precision | **95.2%** | 66.7% | **90.2%** |

### Key Scientific Findings
- **Amplicon HiFi Diagnostic Efficacy:** Amplicon HiFi achieves **95.2% supported event precision** (40 TP vs 2 FP) and 76.0% exact length-pair resolution across the entire DEV cohort.
- **Genomic ONT Limitations:** Unselected whole-genome ONT at 10-20x nominal depth experiences high dropout (`usable_reads` failure on 44 samples) and length distortion (101 samples) due to lack of amplicon enrichment and uncorrected read noise across the repetitive locus.
- **Equal-Length Allele Disambiguation:** On Category C4 equal-length designs (`c4_dupc_50_50_h1_r20_dev_amplicon_hifi` and `c4_dupc_60_60_h2_r35_dev_amplicon_hifi`), haplotagged allele separation achieved 100% exact sequence identity, exact structure match, TP=1, FP=0.

---

## 3. Verification & CI Status
- `make quality`: All checks PASS (Ruff lint, Ruff format, Mypy configured, line limit <650 lines, actionlint).
- `make test-fast`: 753 passed in 5.2s.
- `make ci-check`: 89.80% total branch-aware test coverage (threshold >= 80%).
- Sealed held-out TEST split (100 datasets) remains completely untouched and unblinded.
