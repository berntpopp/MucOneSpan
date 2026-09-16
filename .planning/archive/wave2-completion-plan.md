# Wave 2 Clinical Validation and Bug Fix Implementation Plan

> **For Gemini:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Resolve issues #53, #54, and #55 with evidence-consistent implementations and regression tests, complete the #56 VNTRPipeline comparator evaluation with validated artifacts, rerun the 11 clinical libraries into a new versioned directory (`cohort-v3`), update documentation and PR #57, and verify all CI and quality gates.

**Architecture:**
- **#53 (calling/phasing):** Respect the established phase evidence contract in `call_variants_per_allele`. When within-candidate phase is unphased/conflicting/missing (`"I"`), use selector `"I"`, keep unresolved mixture explicit, and do not fabricate `independent_haplotype_evidence`.
- **#54 (allele selection):** Fix candidate valley selection and read assignment in `alleles.py`. Normalize CIGAR indel lengths by contig reference length (indel error rate) so long biological alleles are not penalized relative to degraded short fragments; filter out or deprioritize sub-biological fragment contigs (< 10 canonical units); preserve spanning read evidence; and ensure the long candidate contig is properly refined and mapped.
- **#55 (report/clinical decision):** Make `compute_clinical_decision` require established pathogenic evidence (`frameshift=True`, verified sequence/VCF support, non-ambiguous localization). Do not declare PATHOGENIC for benign in-frame expansions (such as HG002 WGS 18 bp insertion). Do not declare NEGATIVE when candidate reconstruction is unresolved, unphased, or incomplete (such as HG002 PCR fragment candidate).
- **#56 (comparator):** Execute VNTRPipeline v1.0 source over pinned container image with normalized FASTQ headers and provisioned BiocManager runtime dependency. Validate required output artifacts (FASTAs, TRviz plots, LoF spreadsheets).

**Tech Stack:** Python 3.10+, bcftools, samtools, minimap2, Clair3, Docker, pytest, Ruff, mypy.

---

### Task 1: Fix Issue #53 — Consensus Genotype Selection

**Files:**
- Modify: `src/muc_one_span/calling.py:515-535`
- Test: `tests/unit/test_calling.py`

**Step 1: Write the failing test**
In `tests/unit/test_calling.py`, add `TestDistinctLengthGenotypeSelection`:
- Test 1: Unphased candidate with heterozygous variants receives selector `"I"`, `consensus_haplotype="I"`, `consensus_policy="genotype_iupac_candidate"`, and `independent_haplotype_evidence=False`.
- Test 2: Homozygous candidate receives selector `1`, `consensus_haplotype=1`, and `independent_haplotype_evidence=True` (when `len(allele_keys) > 1`).
- Test 3: Phased candidate receives selector `1`, `consensus_haplotype=1`, and `independent_haplotype_evidence=True`.
- Test 4: Single heterozygous site receives selector `1`, and `independent_haplotype_evidence=True`.

**Step 2: Run test to verify it fails**
Run: `pytest tests/unit/test_calling.py -k TestDistinctLengthGenotypeSelection -v`
Expected: FAIL (unphased candidate receives selector 1 and independent_haplotype_evidence=True).

**Step 3: Implement minimal fix in `calling.py`**
Update `_process_allele` in `call_variants_per_allele`:
```python
variants = parse_vcf_genotypes(filtered)
evidence = phase_evidence(variants)
sample = variants[0].get("sample") if variants else None
is_unphased = evidence["phase_status"] in (
    "unphased",
    "missing_phase_set",
    "disconnected_phase_sets",
    "conflicting_variant_records",
    "missing_genotype",
    "non_diploid",
)
haplotype = "I" if is_unphased else 1
annotate_consensus_candidate(allele_info, evidence, haplotype, sample, str(filtered))
if len(allele_keys) > 1 and not is_unphased:
    allele_info["independent_haplotype_evidence"] = True
```

**Step 4: Run test to verify it passes**
Run: `pytest tests/unit/test_calling.py -k TestDistinctLengthGenotypeSelection -v`
Expected: PASS.

---

### Task 2: Fix Issue #55 — Clinical Interpretation Exceeds Evidence

**Files:**
- Modify: `src/muc_one_span/report.py:75-170`
- Test: `tests/unit/test_report.py`, `tests/unit/test_report_wave1.py`

**Step 1: Write the failing tests**
In `tests/unit/test_report.py`, add tests:
- `test_in_frame_ambiguous_expansion_is_not_pathogenic`: HG002 WGS-style 18 bp expansion (`frameshift=False`, `localization_status="ambiguous"`, `vcf_support=False`) evaluates to `INCONCLUSIVE` (not `PATHOGENIC`).
- `test_unresolved_reconstruction_is_not_negative`: HG002 PCR-style allele with unphased/unverified reconstruction (`independent_haplotype_evidence=False`, or `reconstruction_status="candidate_reference_confidence_unverified"` with unphased mixture) evaluates to `INCONCLUSIVE` (not `NEGATIVE`).
- `test_supported_frameshift_is_pathogenic`: Supported `dupC` (`frameshift=True`, `vcf_support=True`, resolved localization) evaluates to `PATHOGENIC`.

**Step 2: Run test to verify it fails**
Run: `pytest tests/unit/test_report.py -k "test_in_frame_ambiguous_expansion_is_not_pathogenic or test_unresolved_reconstruction_is_not_negative" -v`
Expected: FAIL.

**Step 3: Implement minimal fix in `report.py`**
Update `compute_clinical_decision`:
1. Distinguish confirmed pathogenic mutations from ambiguous/uncertain sequence changes:
   - A mutation is pathogenic if `frameshift is True`, `vcf_support is True` (or concordant support status), and `localization_status != "ambiguous"`.
   - Other mutations (in-frame, unsupported, or ambiguous localization) are recorded as uncertain findings.
2. Inconclusive checks:
   - Low coverage (< 30 reads).
   - Ambiguous bases (> 10).
   - Execution warning.
   - Presence of uncertain/unresolved sequence findings.
   - Incomplete or unresolved allele reconstruction (e.g. `independent_haplotype_evidence` is False when diploid, or unphased/unverified reconstruction status).
3. Negative only when:
   - Complete, resolved reconstruction on both alleles without warnings, high quality, and zero pathogenic or uncertain findings.

**Step 4: Run test to verify it passes**
Run: `pytest tests/unit/test_report.py tests/unit/test_report_wave1.py -v`
Expected: PASS.

---

### Task 3: Fix Issue #54 — Allele Selection and Read Assignment

**Files:**
- Modify: `src/muc_one_span/alleles.py:150-200, 280-360, 480-550`
- Test: `tests/unit/test_alleles.py`

**Step 1: Write the failing tests**
In `tests/unit/test_alleles.py`:
- Add synthetic test with two biological alleles (e.g. canonical 41 and canonical 69) and a population of short fragments (canonical 5):
  Verify `_split_cluster_by_indel` and `detect_alleles` selects the two biological alleles (41 and 69), NOT the fragment contig (5).
- Add synthetic test with long full-span reads and unequal allelic representation:
  Verify long reads are assigned to the long candidate, and the short candidate does not steal long spanning reads.

**Step 2: Run test to verify it fails**
Run: `pytest tests/unit/test_alleles.py -k test_synthetic_unequal_alleles_with_fragments -v`
Expected: FAIL.

**Step 3: Implement fix in `alleles.py`**
1. In `_split_cluster_by_indel`:
   - Compute `normalized_indel_rate`: `mean_indel / contig_reference_length` (where reference length is flank (1000 bp) + (9 + c) * 60 bp).
   - Filter out or deprioritize fragment contigs with $c < 10$ canonical units when valid candidate contigs with $c \ge 10$ exist.
   - Choose the two best valleys by normalized indel rate and read support/depth.
2. In `refine_peak_contig`:
   - Ensure the read count threshold does not eliminate a valid second allele mode when sub-clusters are partitioned.
3. Verify read assignment: reads in `sc1` and `sc2` extract to their respective contigs.

**Step 4: Run test to verify it passes**
Run: `pytest tests/unit/test_alleles.py -v`
Expected: PASS.

---

### Task 4: Complete Comparator Evaluation (#56)

**Files:**
- Script: `run_comparator.py`
- Provision: Mount BiocManager to `/opt/conda/envs/python-env/lib/R/library/BiocManager`
- Inputs: `comparator-inputs/` (normalized FASTQ headers)
- Outputs: `comparator-v1/cohort/` (new run directory)

**Step 1: Run comparator across 11 clinical runs**
Execute `python run_comparator.py` with the 11 runs.
Monitor logs and resource usage.

**Step 2: Validate comparator artifacts**
Verify final FASTAs, TRviz plots, and `new_and_lof_seqs.xlsx` for each run.
Score outcomes: generic LoF vs named dupC, HG002 Q100 sequence concordance.

---

### Task 5: Rerun Clinical Cohort (`cohort-v3`) and Validate

**Step 1: Run full 11-library clinical benchmark**
Execute `scripts/clinical_benchmark.py run` into `cohort-v3`.
Execute `scripts/clinical_benchmark.py score` to generate `results-v3.json`.

**Step 2: Comparative analysis**
Compare `cohort-v2` vs `cohort-v3`:
- MP1: supported B:dupC recovery.
- MP2 / MP4: long allele candidate selection, Clair3 calling, consensus, mutation recovery.
- MP3: unphased genotype handling.
- HG002 WGS: in-frame expansion clinical decision (`INCONCLUSIVE` instead of `PATHOGENIC`).
- HG002 PCR: fragment candidate clinical decision (`INCONCLUSIVE` instead of `NEGATIVE`).

**Step 3: Simulation gates**
Run existing simulation and fast regression tests to ensure no regressions on normal controls.

---

### Task 6: Quality Gates, Independent Review, and Delivery

**Step 1: Run all repository gates**
- `make ci-check` (unit tests + branch coverage >= 80% + file sizes < 650 lines + Ruff + mypy).
- `make test-int` (integration tests).
- `make docs-check` (strict MkDocs build).
- `make build-check` (wheel/sdist build).

**Step 2: Independent review**
Conduct independent review of scientific correctness, regressions, and documentation.

**Step 3: Update documentation, changelog, issues, and PR #57**
- Update `docs/guides/clinical-validation-results.md`.
- Update `CHANGELOG.md`.
- Update issue statuses on #53, #54, #55, #56, and PR #57.
