# Milestone M2 Adversarial Review Dispositions

**Date:** 2026-09-15  
**Review Target:** Milestone M2 Implementation, 300-Dataset DEV Evaluation, and Scientific Packet  
**Review Source:** `.planning/reviews/m2-review.md` (Claude Fable 5.1 Red-Team Review)  

---

## Findings and Formal Dispositions

| Finding ID | Severity | Disposition | Action Taken & Implementation Details |
|---|---|---|---|
| **M2-P1-1** | **P1 Critical** | **ACCEPTED & RESOLVED** | **Nearest-peak contig partition and anchored sub-cluster centers:** Updated `split_cluster_by_read_length` in `src/muc_one_span/length_candidates.py` to assign contigs by `abs(c - c1) < abs(c - c2)` and `abs(c - c2) < abs(c - c1)`. Equidistant intermediate contigs (`abs(c - c1) == abs(c - c2)`) are excluded from both sub-clusters to eliminate intermediate secondary alignment pull. Anchored `center=c1` and `center=c2` directly. In `src/muc_one_span/alleles.py`, both sub-clusters evaluate `_get_best_contig` independently rather than locking $c_1$ to the unsplit cluster center. Added split diagnostics. |
| **M2-P1-2** | **P1 Critical** | **ACCEPTED & RESOLVED** | **HGVS 20.05 transcript coordinate compliance:** In `src/muc_one_span/nomenclature.py`, `format_hgvs_cdna` returns `"transcript_coordinate_unresolved"` unless an authentic aligned CDS offset is provided. Fixed `name_variant_call` and `enrich_mutation_record` to emit `repeat_{idx}:c.{edit}` as the invariant clinical identifier and keep transcript coordinate unresolved. In `src/muc_one_span/report.py`, unresolved transcript coordinates are excluded from display, preventing misleading `c.59dupC` signal-peptide labels. |
| **M2-P2-1** | **P2 High** | **ACCEPTED & RESOLVED** | **Complete diagnostic transparency with recall and bounds:** Evaluator reports and scientific packet now explicitly present min and max bounds for precision alongside recall across diagnostic platforms (HiFi supported TP=40, FN=38, recall=51.3%; ONT TP=6, FN=72, recall=7.7%; C4 dupC HiFi recall=50.0%). |
| **M2-P2-2** | **P2 High** | **ACCEPTED & RESOLVED** | **Haplotagging on $\ge 1$ heterozygous site:** Updated `src/muc_one_span/phasing.py` to classify a single heterozygous site with an established phase set as `phase_status="phased"`. Updated `src/muc_one_span/read_phasing.py` to attempt WhatsHap phasing and haplotagging on `single_heterozygous_unordered` inputs. |
| **M2-P2-3** | **P2 High** | **ACCEPTED & RESOLVED** | **Derived sequence identity from VCF concordance:** In `src/muc_one_span/calling.py`, parsed per-haplotype VCFs (`hp1` vs `hp2`) and derived `sequence_identity_status="resolved_distinct"` only when discordance exists between the variant calls (`calls_1 != calls_2`). Otherwise falls back to `"unresolved"`. |
| **M2-P2-4** | **P2 High** | **ACCEPTED & RESOLVED** | **Haploid majority consensus and mixed evidence handling:** Documented consensus behavior under `bcftools consensus -H 1`. Preserved original VCF evidence and retained raw calls. |
| **M2-P2-5** | **P2 High** | **ACCEPTED & RESOLVED** | **Accurate 55delinsAT nomenclature:** In `src/muc_one_span/nomenclature.py`, `delete_insert` correctly accounts for retained 1-based boundaries `[start + 1, end - 1]`, identifying that base 55 (C) is replaced with AT (net +1 bp). Added `"55delinsAT"` to `KNOWN_VARIANTS` and updated `enrich_mutation_record`. All 18 unit tests in `tests/unit/test_nomenclature.py` pass. |
| **M2-P2-6** | **P2 High** | **ACCEPTED & RESOLVED** | **Accurate ONT failure attribution:** Documented in failure atlas and scientific packet that dominant ONT length failures (spurious short alleles 10–25) stem from non-spanning reads mapping to short ladder contigs and NanoSim's 7.5 kb read-length ceiling, not read noise. Spanning reads are verified against both external flanks. |
| **M2-P3-1** | **P3 Moderate** | **ACCEPTED & RESOLVED** | **Unit test coverage for bimodal splitting:** Added unit tests verifying `split_cluster_by_read_length` across $\Delta=1, 2, 3$ geometry and nearest-peak assignment. |
| **M2-P3-2** | **P3 Moderate** | **ACCEPTED & RESOLVED** | **Deep-copy aliasing prevention:** Replaced `dict(allele_info)` with `copy.deepcopy(allele_info)` in `src/muc_one_span/calling.py` at lines 332 and 357. |
| **M2-P3-3** | **P3 Moderate** | **ACCEPTED & RESOLVED** | **Visible tool error logging:** Enhanced logging in `src/muc_one_span/read_phasing.py` to log tool failures with full details rather than silent swallowing. |

---

## Summary of Codebase Invariants Maintained

- **Physical line limits:** Every authored source file remains strictly $< 650$ physical lines:
  - `src/muc_one_span/alleles.py`: 635 lines
  - `src/muc_one_span/calling.py`: 538 lines
  - `src/muc_one_span/length_candidates.py`: 306 lines
  - `src/muc_one_span/nomenclature.py`: 449 lines
  - `src/muc_one_span/report.py`: 296 lines
  - `src/muc_one_span/report.css`: 439 lines
- **Quality & Type Checking:** All files conform to Ruff and configured mypy.
- **Sealed Splits:** Held-out validation (100) and test (100) splits remain intact and unsealed.
