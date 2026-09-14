I have everything needed. Writing the review now.

**Verdict: REVISE AND RESUBMIT.** The HiFi gains are real and regression-free, but the ONT gate is missed by a wide margin, the dominance test has a hole that already manufactured a false allele on a homozygous dev control, and several claims in the request are contradicted by the reports. I could not run `make ci-check` or timing in this session, so the software-contract and runtime gates are unverified here.

## Gate scorecard (140 dev datasets, paired against v0.11.0)

I rebuilt the paired comparison from both evaluation reports rather than relying on the comparison script, which only checks the four sequence metrics. The baseline failure atlas was produced by an older diagnosis script and disagrees with its own report, so I did not use it.

| Gate (spec §3) | Required | Observed | Status |
|---|---|---|---|
| HiFi paired exact-diploid net | ≥ +8 | +14 (14 wins, 0 losses) | Pass |
| ONT paired exact-diploid net | ≥ +5 | +1 (1 win, 0 losses) | **Fail** |
| HiFi exact diploid recovery gain | ≥ +12 pp | 4/70 to 18/70, +20 pp | Pass |
| Count pairs, wins > losses, no regressions | | 6 wins, 0 losses | Pass |
| Mutation sensitivity ≥ baseline | | sample alarms 72 to 69; exact events 53 to 53 (4 lost, 4 gained) | **Fail** |
| Exact annotation accuracy | ≥ +5 pp | 0 pp | **Fail** |
| Annotation false discoveries decrease | | event FP 62 to 23 | Pass |
| Zero new control false alarms | | 9 to 2, both pre-existing | Pass |
| Absolute control alarm cap | ≤ 1 per platform | HiFi 2, ONT 0 | **Fail** |
| Confident-negative rate | ≥ 65% | 11/44 = 25% | **Fail** |
| Identical-length no false split | ≥ 85% | 13/14, but the one split raised an alarm | Marginal |
| Runtime ≤ +15% | | not measured | Unverified |
| Parameters via typed settings | | five new thresholds hard-coded | **Fail** |

The four lost exact-event true positives are mut_asym_ins16bp_20_90_h2_ont, mut_eq_inscccc_60_60_h1_hifi, mut_gap1_insg_50_51_h1_ont and mut_gap3_dupc_60_63_h1_ont. The ten lost sample alarms all had spurious baseline events, so that drop is mostly cleanup, but the gate text does not distinguish, and the exact-event ledger is flat. The request's claim that all 44 variant-calling failures were eliminated is not supported: exact detections did not increase, only their support status did.

## Answers to the review questions

**Q1a, homozygous protection.** Not leak-free. In read_dominance.py a read with a record on the candidate contig and none on the primary contig is counted as dominant. Mapping runs minimap2 with preset defaults, so each read carries at most five secondaries, and tail reads five or more contigs from the peak never carry a primary-contig record. The test then degenerates to "three records on a distant contig". This is exactly what happened on ctrl_ident_80_80_hifi: the candidate emitted a 75-repeat second allele built from 19 records with 3 primaries, then raised 3 false mutation events. The request's statement that identical controls are "100% rejected" is false for that sample. The design review explicitly asked for pairwise rescoring of all spanning reads against both candidates; the implementation reuses ladder-BAM records instead. Fix by realigning the candidate's reads to the primary contig, or at minimum treating a missing score as ambiguous rather than dominant and requiring primary-record support.

**Q1b, the 1% ratio.** Mathematically harmless but vacuous. With a 3-read floor the ratio only binds when the primary side exceeds 297 decisive reads, which no dev sample reaches. The operative rule is the absolute floor, and 3 records is at the edge of the spec's own "fewer than 3 reads must stay unresolved" clause. The floor is also compared against alignment records, not molecules, which the code elsewhere is careful to distinguish. The ONT floor of 4 is not what the request states.

**Q1c, haploid re-genotyping.** Justified only under an isolation assumption the code never checks. The rewrite in vcf.py forces any site at or above 0.5 allele fraction to homozygous ALT regardless of depth, keeps 0/0 records in the file, and lowers the effective QUAL floor to 4.0 through a hard-coded constant that silently overrides a user's `--min-qual`. The design review recommended keeping QUAL at 5 and validating partition purity; neither was done. The consensus policy label still says "genotype_iupac_candidate" although no IUPAC codes can now be produced, so provenance is misleading.

**Q2.** Covered by the table. ONT fails outright. The dominant ONT failure is a systematic +1 overcall of the longer allele at 80 units and above on at least 14 datasets, unchanged from baseline, and nothing in the candidate addresses it.

**Q3, architecture and overfitting.**
- Dead code: the discovery and validation functions in length_candidates.py are never called from the pipeline, so the quoted 91% coverage measures tests of unused code. The design's read_partition.py was not built.
- Hard-coded tuning: the score margins, read floors, ratio, the record floor of 3 in alleles.py, and the QUAL 4.0 floor are all module constants, violating the config-driven contract and the repo's settings architecture.
- The inner 50/20 tuning and selection split from the spec was not used, and the evaluation runner reads the sealed ledger with design names visible. All dev-set thresholds are in-sample choices. The same runner would unblind the final split unless a public-ledger mode is added.
- The failure atlas labels every equal-length design as a consensus failure because the evaluator cannot credit an unproven duplicate. That is an evaluator artefact and it hides the real unresolved-multiplicity problem.
- The reference ladder and its index stayed byte-compatible only because both flank variants are 500 bp. No test asserts the bundled ladder is reproducible from the dictionary.

**Q4, required before resubmission.**
1. Replace the missing-score rule with explicit realignment or an ambiguous default, add a primary-record requirement, and re-run the identical and equal-variant controls.
2. Address the ONT +1 bias or the ONT gate cannot be met.
3. Move every new threshold into typed settings and remove the QUAL override or make it an explicit setting.
4. Declare which sensitivity and annotation metrics are the gate metrics, and explain the four lost exact events.
5. Wire or delete the unused discovery functions.
6. Report paired runtime, run the inner split or a threshold sensitivity sweep, and add a blinded runner for the final split.