# Empirical Analysis: Sequence vs. Length Discordance in MucOneSpan Consensus Reconstruction

**Date:** 2026-09-15  
**Subject:** Root-Cause Investigation of Sequence Mismatches in Length-Exact Alleles  
**Author:** Antigravity (Pair Programming with User)  

---

## 1. Context & Motivation

During benchmark evaluation of `MucOneSpan` across the 44-sample historical cohort and the 200 simulation dataset catalog (100 diploid designs × 2 platforms), an interesting divergence was observed:
- **Repeat Count Accuracy (`count_exact`):** High (93.2% in standard HiFi, 68.0% across the 200 challenge panel).
- **Diploid Sequence Exactness (`all_sequences_exact`):** 44.3% in standard HiFi, 23.5% across the 200 challenge panel.

Specifically, in **48.9% of all alleles (43 out of 88 in the HiFi cohort)**, the called allele length was **100% exact ($\Delta = 0$ bp error)**, yet the nucleotide sequence failed strict identity (`sequence_exact: false`). The edit distance was consistently small (mean: 35 bp, median: 2 bp).

The user raised three key questions:
1. *Are SNPs inserted into the sequence by MucOneUp?*
2. *Is the wrong repeat unit inserted?*
3. *Are we missing repeat units in our catalogue?*

---

## 2. Empirical Findings

### 2.1 Catalogue Completeness (Question 3)
A programmatic audit of `src/muc_one_span/data/repeats/repeats.json` against all 100 simulation truth structures in `examples/experiment_100_designs.json` and historical fixtures revealed:
- **Catalogue units:** 34 total (`1`, `2`, `3`, `4`, `4p`, `5`, `5C`, `6`, `6p`, `7`, `8`, `9`, `X`, and variants `A` through `W`).
- **Missing units:** **0**. Every repeat unit present in truth designs is represented in the catalogue.

### 2.2 MucOneUp Sequence Generation (Question 1)
MucOneUp generates biological haplotypes by chaining authentic repeat units according to the empirical Vrbacka Markov transition model, with targeted mutations (e.g. `dupC`, `del18_31`) placed at explicit repeat indices. MucOneUp does **not** insert synthetic random background SNPs into repeat units.

### 2.3 The Primary Root Cause: IUPAC Ambiguity Codes from Clair3 `0/1` Calls
Analysis of the 43 length-exact/sequence-mismatched alleles in `hifi-evaluation.json` revealed that **41 out of 43 (95.3%) contain IUPAC ambiguity codes** in the consensus sequence:

| Metric | Value |
|---|---|
| Total evaluated allele pairs | 88 |
| Pairs with `repeat_count_error == 0` | 82 (93.2%) |
| Pairs with `repeat_count_error == 0` AND `sequence_exact == false` | 43 (48.9%) |
| **Length-exact mismatches containing IUPAC ambiguous bases** | **41 / 43 (95.3%)** |

### 2.4 Concrete Case Walkthroughs

#### Case A: `sample_bench_5000` (Allele 2)
- **Truth Length:** 63 repeat units (3,780 bp)
- **Predicted Length:** 63 repeat units (3,780 bp) $\to$ **$\Delta = 0$ bp length error!**
- **Sequence Edit Distance:** 2 bp
- **Truth Structure:** `12345CXAXAABXXXDECXXXXXXXBXAABXXXXXXABXXXXXVVGABXXXXXXBXXXX6789`
- **Pred Structure:** `12345CXAXAABXXXDECXXXXXXXBXAABXXXXXXABXXXXXVVGABXXXXXX?XXXXX6789`
- **Discrepancy:**
  - Base pos 3266 (repeat #54, offset 26): truth=`G`, pred=`S`
  - Base pos 3268 (repeat #54, offset 28): truth=`G`, pred=`S`
- **Root Cause:**
  At repeat #54, the biological truth is unit `B` (carrying two Gs relative to canonical `X`). In Clair3 variant calling on the extracted allele reads, Clair3 observed 14 reads supporting reference `C` and 44 reads supporting `G` (AF = 0.73). Clair3 emitted genotype `0/1`. `bcftools consensus --iupac-codes` converted `0/1` C/G into **`S`**. In strict sequence comparison, `S != G` (edit distance = 2), and `classify.py` labeled the repeat containing `S` as `?` (unknown unit).

#### Case B: `sample_bench_5001` (Allele 2)
- **Truth Length:** 65 repeat units (3,780 bp)
- **Predicted Length:** 65 repeat units (3,780 bp) $\to$ **$\Delta = 0$ bp length error!**
- **Sequence Edit Distance:** 2 bp
- **Discrepancy:**
  - Base pos 1193: truth=`G`, pred=`S` (`contig_56:1694 C→G GT=0/1 AF=0.73`)
  - Base pos 1258: truth=`A`, pred=`M` (`contig_56:1759 C→A GT=0/1 AF=0.77`)
- **Structure:** Truth has `DECF`, consensus was classified as `D?E?XF`.

---

## 3. Secondary Causes (~4.7% of Mismatches)

1. **Clair3 False Negatives (Reference Fill):**
   In low-coverage or high-noise windows, Clair3 may fail to call a true SNP in a variant repeat. The consensus sequence defaults to the reference contig's canonical `X` unit. The repeat count remains exact, but that unit is classified as `X` instead of `A` or `B` (1–2 bp substitution).
2. **ONT 1-bp Homopolymer Indels:**
   In Oxford Nanopore reads, 1-bp indel compression within the 7-C homopolymer tract of unit `X` can lead to slight edit distance shifts.

---

## 4. Architectural Recommendations for Future Releases

1. **Haploid-Forced Consensus Replay:**
   Since reads in `allele_N` have already been partitioned by read dominance to represent a single physical allele, variant calling on an isolated allele should treat dominant variants (e.g. `AF >= 0.60`) as homozygous (`1/1`) for consensus generation.
2. **IUPAC Filtering / Disambiguation:**
   Add a consensus post-processing or VCF filtering step that maps high-frequency heterozygous calls (`0/1` with high ALT depth) to the major allele rather than emitting ambiguous IUPAC symbols, eliminating the `?` repeat-unit classification artifact.
