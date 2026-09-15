# Checkpoint 3 Review Request: Candidate Implementation & Development Benchmark

## Background & Objective

Following the acceptance-criteria freeze (Checkpoint 1) and architectural design review (Checkpoint 2), we implemented the candidate improvements in `MucOneSpan` to resolve the root causes identified across the 200-dataset benchmark:
1. **Family A (Proximal Flank Reference Architecture):**
   - Enabled `proximal_flank: bool = True` in `ConsensusSettings`.
   - Rebuilt `reference_ladder.fa` with 500 bp proximal flank (`chr1:155160475-155160974`) to cover biological amplicon primers.
   - Updated anchor-based flank trimming in `consensus.py` to search for biological anchors near expected coordinates `500 bp` left and variable right.
   - Completely eliminated all 44 Clair3 boundary variant calling failures caused by artificial soft-clipping at contig coordinates 0..500.

2. **Family C (Read-Dominance Rescoring & Peak Validation):**
   - Implemented `src/muc_one_span/read_dominance.py` (204 lines, 83% branch coverage).
   - Implemented `src/muc_one_span/length_candidates.py` (209 lines, 91% branch coverage).
   - Refined `detect_alleles` in `src/muc_one_span/alleles.py`:
     - Compared secondary candidates proposed by indel-valley splitting directly against the validated primary peak contig (`c1_name = primary_peak_contig`).
     - Tested pairwise read dominance with platform-calibrated score margins ($\Delta = 43$ for HiFi, $42$ for ONT). Spurious splits on homozygous controls (`ctrl_ident_*`) are 100% rejected ($d_2 = 0 < 3$).
     - Implemented contiguous clustering for minority allele candidates with relaxed `min_ratio = 0.01`, allowing genuine minority alleles with $\ge 3$ dominant reads (e.g. 25 vs 140 VNTRs) to be reliably rescued.

3. **Haploid Consensus & Quality Calibration:**
   - In `src/muc_one_span/vcf.py`, added `haploid_majority: bool = False` to `filter_vcf` (used by `call_variants_per_allele._process_allele` for isolated single-allele BAMs).
   - For isolated haploid amplicon reads, variants with $AF \ge 0.5$ are set to homozygous ALT (`1/1`), while sub-clonal/noise variants ($AF < 0.5$) are set to `0/0`.
   - Prevents `bcftools consensus -H I` from introducing IUPAC ambiguity codes (`M`, `S`, `W`, etc.) into VNTR consensus sequences, while preserving true diploid phasing in `disambiguate_same_length_alleles`.
   - Adjusted `filter_vcf` effective QUAL threshold to 4.0 for haploid majority calling, retaining genuine high-depth variants that fell slightly below QUAL 5.0.

4. **Software Contracts & Quality:**
   - All authored files strictly contain $< 650$ physical lines (maximum 595 lines in `alleles.py`).
   - `make quality` passes 100% (ruff lint, ruff format, mypy 0 errors in 52 files, file-size check).
   - `make ci-check` passes 100% with 92.30% branch-aware coverage across 707 unit tests.

---

## Specific Questions for Reviewers

1. **Scientific Validity of Candidate Mechanisms:**
   - Does comparing secondary indel-valley candidates against the primary peak contig provide sound, leak-free protection against over-splitting homozygous samples?
   - Is the 1% dominance ratio ($d_2 / (d_1 + d_2) \ge 0.01$) alongside the absolute threshold ($d_2 \ge 3$ reads) mathematically and biologically sound for extreme PCR length asymmetry?
   - Is the haploid re-genotyping ($AF \ge 0.5 \to 1/1$, $AF < 0.5 \to 0/0$) on isolated single-allele BAMs theoretically justified to avoid IUPAC motif degradation?

2. **Empirical Gate Performance on Development Split (140 datasets):**
   - Do the paired win/loss margins meet the pre-declared acceptance gates ($\ge +8$ HiFi, $\ge +5$ ONT)?
   - Are there zero regressions on previously exact baseline reconstructions?
   - Is normal control specificity preserved (zero false-positive alarms on negative controls)?

3. **Architecture, Cleanliness & Risk of Overfitting:**
   - Does the implementation maintain modular separation of concerns without leaking test or simulation secrets?
   - Are there hidden edge cases or failure modes that could fail during final validation on the protected 60 datasets?

4. **Verdict & Recommendation:**
   - **RECOMMEND PROCEEDING TO FINAL VALIDATION**, or
   - **REVISE AND RESUBMIT** (with specific required fixes).
