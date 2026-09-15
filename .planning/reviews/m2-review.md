# Milestone M2 Adversarial Review

**Reviewer Persona:** Principal Computational Genomicist & Staff Systems Engineer (Adversarial Red Team)  
**Date:** 2026-09-15  
**Target:** Milestone M2 Implementation, 300-Dataset DEV Evaluation, and Scientific Packet

---

## Executive Summary of Findings

| ID | Severity | Title | Impact |
|---|---|---|---|
| M2-P1-1 | **P1 Critical** | Bimodal split works only for $\Delta=1$; $\Delta \ge 2$ gives systematic $+1$ error on shorter allele | Length inference error on 40/42, 50/52, 60/62, etc. (at least 25 DEV designs) |
| M2-P1-2 | **P1 Critical** | `hgvs_cdna` emits accessioned but wrong transcript coordinate (`c.59dupC` in signal peptide) | Clinically invalid HGVS coordinate in clinical report |
| M2-P2-1 | **P2 High** | Supported-event precision reported with favorable bound; recall omitted | Diagnostic yield transparency (HiFi recall ~51%, ONT ~8%) |
| M2-P2-2 | **P2 High** | Haplotagged separation requires $\ge 2$ heterozygous sites in mixed VCF | Low-heterozygosity founder alleles fail to phase |
| M2-P2-3 | **P2 High** | Distinctness asserted before observed in `calling.py` | Over-claiming resolved haplotypes prior to VCF comparison |
| M2-P2-4 | **P2 High** | AF 0.2–0.5 retention band resolves to REF in `bcftools consensus -H 1` | Retained mixed variant not surfaced in consensus sequence |
| M2-P2-5 | **P2 High** | `delinsAT` nomenclature contradicts truth (`55delinsAT` net +1 vs `54_56delinsAT` net -1) | False negative on rare delinsAT variants |
| M2-P2-6 | **P2 High** | ONT failure mode misattributed (non-spanning reads & 7.5kb simulator cap) | Theoretical vs physical root cause attribution |
| M2-P3-1 | **P3 Moderate** | Missing unit tests for `split_cluster_by_read_length` & geometry assumptions | Test coverage gap for new bimodal logic |
| M2-P3-2 | **P3 Moderate** | Shallow-copy aliasing of `allele_2` in `calling.py` | Mutable dict aliasing risk |
| M2-P3-3 | **P3 Moderate** | Tool exception swallowing in `read_phasing.py` | Violates visible failure contract |

---

## Detailed Critiques

### P1. Bimodal split works only for Δ=1; Δ≥2 gives a systematic +1 error on the shorter allele
C3 HiFi: 0/24 strict exact, 11 correct length pairs, 13 length failures. Every Δ=1 design (40/41, 50/51, 60/61, 70/71, 90/91, 100/101, 110/111) gets correct lengths; every Δ=2 or Δ=3 design fails with the pattern 40/42→[41,42], 50/52→[51,52], 60/62→[61,64], 70/73→[71,71], 90/93→[91,91]. The same signature appears in C4 (60/62→[61,62]), C5, and C6, so this is at least 25 DEV failures.
**Mechanism:** `split_point=(c1+c2)//2` at `length_candidates.py:279` puts the intermediate contig into the lower sub-cluster; that contig collects secondary alignments from both alleles, and the weighted center plus `refine_peak_contig` select it.
**Recommendation:** Assign contigs to the nearest peak rather than a midpoint, count only primary alignments when building sub-cluster centers, record split diagnostics (bins, peaks, delta, which splitter fired, dominance verdict) in `alleles.json`, and add unit tests.

### P1. `hgvs_cdna` emits an accessioned but wrong transcript coordinate
`nomenclature.py:228-243` produces `NM_001204286.1:c.59dupC`. Position c.59 of that transcript lies in the signal-peptide-coding region, not in the VNTR, so the string describes a different variant. M2 added `repeat_relative_coordinate` and `transcript_coordinate_unresolved` but left the wrong string in place, and it feeds the clinical HTML report through `report.py:33`.
**Recommendation:** Stop emitting any `NM_:c.` string unless computed from the reconstructed allele's actual offset; label `repeat_{idx}:c.` explicitly as non-HGVS repeat-relative syntax; emit a genomic `g.` description only if derived from a real alignment.

### P2. Supported-event precision is reported with the favorable bound and recall is omitted
The report gives FP min 5 / max 7 and its own precision as 46/53 = 86.8% (`evaluation_report.json:81-84, 133-137`). The packet's 90.2% uses min FP with min TP; the ablation JSON's per-category FP sums to 7. Recall on the diagnostic platform is absent from the packet (HiFi Supported TP=40, FN=38, Recall=51%; Genomic ONT TP=6, FN=72, Recall=8%; C4 dupC HiFi TP=15, FN=15, Recall=50%).
**Recommendation:** Report min/max bounds as the evaluator does, always alongside recall.

### P2. Haplotagged separation cannot bootstrap on low heterozygosity
The path requires ≥2 PASS heterozygous sites in the mixed VCF (`read_phasing.py:188`) before WhatsHap runs; the 50/50 example had 38 het SNVs because MucOneUp draws random unit compositions per haplotype. A founder dupC on two alleles sharing a background yields 0 or 1 het sites and the dupC stays suppressed. Single-het cases are skipped although WhatsHap haplotag can tag on one phased site.
**Recommendation:** Support phasing and haplotagging on $\ge 1$ heterozygous site; implement homopolymer length tally fallback at the 53-59 C-tract when het-site count is below 1.

### P2. Distinctness is asserted before it is observed
`calling.py:317-319` sets `sequence_identity_status="resolved_distinct"`, `phase_status="phased"`, and `homozygous=False` before parsing the per-haplotype VCFs; per-allele fields still read `unresolved` in the C4 outputs, and `observed_length_candidates`/`allele_multiplicity_status` are never updated.
**Recommendation:** Derive the status from a comparison of the two PASS variant sets, require ≥1 discordant site, and record hp1/hp2/untagged read counts.

### P2. The AF 0.2–0.5 retention band does not surface mixed evidence
`vcf.py:123-130` keeps `0/1` in a supposedly haploid VCF, but consensus uses `bcftools consensus -H 1`, which silently resolves `0/1` to REF.
**Recommendation:** When a retained heterozygous call survives on a haplotype BAM, emit a contamination flag on that allele and degrade its identity status.

### P2. delinsAT nomenclature contradicts the truth
The truth mutated unit is `...CCGCCATCCCCA` (61 bp, net +1, one C replaced by AT), so the correct name is `55delinsAT`. `KNOWN_VARIANTS` at `nomenclature.py:28` and the inclusive reading of `repeats.json:216-228` give `54_56delinsAT` (net −1). The docs diff in this tree already says "net +1", so the code and docs now disagree.
**Recommendation:** Fix the dictionary semantics, the known-variant key, and the concordance projection together, with a test against the mutated-unit FASTA.

### P2. ONT failure mode is misattributed
The dominant ONT length failure is a spurious short allele of 10–25 repeats (25/100→[14,25], 60/60→[13,59], 80/80→[10,34]), not "uncorrected read noise". Non-spanning genomic reads map to short ladder contigs, and NanoSim `max_read_length: 7500` makes any allele above roughly 95 repeats physically unspannable.
**Recommendation:** Require reads to anchor in both external flanks before they count toward a length candidate, report spanning-read counts, and state that the ONT arm is length-limited by design.

### P3. Split algorithm robustness, haplotag hygiene, and nomenclature edge cases
1. Deep-copy `allele_2` in `calling.py` to prevent aliasing.
2. Ensure visible errors in `read_phasing.py` (do not swallow tool failures).
3. Unit tests for `split_cluster_by_read_length`.
4. Fix `read_phase` configuration default so evaluation is reproducible without overrides.
