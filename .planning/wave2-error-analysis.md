# Wave 2 scientific error analysis — verified final baseline

Audit covers all eleven completed runs in external `cohort-v2`, scored in `results-v2.json`, using immutable harness commit `1deb5dcfd1dbeede0fbd3ca859cc9a6aab011d47` and the unchanged v0.14.1 caller/scientific settings/model. Independently verified every output hash listed in all eleven final benchmark records and confirmed final scoring JSON equals initial `results-v1.json` exactly. The lead's `cohort-v1-v2-comparison.json` records306 scientific artifacts before/after and zero normalized differences. A separate independent comparison of295 scientific files (196VCFs,66consensus FASTA/context files,33alleles/summary/repeats files) likewise found zero differences after normalizing output paths, volatile VCF command/date headers and JSON ordering; summary comparison retained allele and classification content. Thus the previously identified stage findings are confirmed in the final baseline, not projections from unfinished runs. No scientific code, thresholds, preprocessing, models or settings were changed by this audit.

## Per-run stage observations

Selected lengths are total boundary-inclusive repeat candidates; parenthesized canonical count is obtained by subtracting the ladder's nine fixed units. “Callable” uses the benchmark's final evidence gate, not merely successful process exit. All eleven pipelines and reports completed; six are callable and five retain ambiguous evidence.

| Run / sample | Selected total lengths | Callable | Clinical HTML | Stage observation |
|---|---|---|---|---|
| ERR15277562 / HG001 PCR | 40,72 | yes | NEGATIVE | Lengths agree with published comparator; no independent complete sequence or negative event truth. |
| ERR15277563 / HG002 PCR | 65,15 | no | NEGATIVE | Long Q100 haplotype not reconstructed; second candidate unphased. Literal paternal sequence exact; second900bp candidate mismatches maternal4638bp. |
| ERR15277564 / HG003 PCR | 44,65 | yes | NEGATIVE | Comparator length concordance only; both candidate reference confidences unverified. |
| ERR15277565 / HG004 PCR | 44,10 | yes | NEGATIVE | Long raw physical-span population near4.6kb collapsed within broad candidate cluster; very short second candidate selected. |
| ERR15277566 / MP1 PCR | 39,78 | yes | PATHOGENIC | Raw VCF has dupC-compatible insertion at expected tract; final classification instead describes insertion A, without named dupC or VCF support. |
| ERR15277567 / MP2 PCR | 50,14 | no | NEGATIVE | Long reported control allele not retained; second candidate unphased; filtered VCF contains substitutions only. |
| ERR15277568 / MP3 PCR | 37,133 | no | PATHOGENIC | Raw short-candidate VCF has heterozygous dupC-compatible insertion; GT1 replay selects reference allele. Long candidate has ambiguous dupC-like template label with no VCF support. |
| ERR15277569 / MP4 PCR | 45,15 | no | NEGATIVE | Long reported control allele not retained; second candidate unphased; filtered VCF contains substitutions only. |
| ERR15277570 / MP5 PCR | 61,69 | yes | PATHOGENIC | Generic C insertion has exact VCF/sequence support, but no named dupC; comparator-only novel event remains exploratory. |
| ERR15277552 / HG002 WGS | 65,77 | no | PATHOGENIC | Both literal sequences exactly match independent Q100; biological18bp in-frame expansion is flagged as ambiguous boundary mutation and triggers pathogenic display. |
| ERR15277553 / patient WGS, ENA MP1 alias | 44,32 | yes | NEGATIVE | Few exact full-span inputs, second candidate VCF empty; participant identity unresolved, no reported-control scoring. |

## Aggregate interpretation and denominators

Final v2 execution11/11; primary amplicon callability5/9, secondary WGS1/2. Reported-known-control exact dupC recovery0/4 across MP1–MP4, with three controls uncallable; conditional callable-only recovery0/1. MP1/MP2 are siblings, so four control participants are not four independent families. There are no sample-linked documented orthogonal assay confirmations or independently established negative events in the curated ledger: independently confirmed sensitivity and specificity remain not estimable. WGS patient identity is unresolved and must not become a fifth positive or a technical replicate of MP1 through alias inference.

The broad PATHOGENIC alarm is present for MP1, MP3, MP5 and HG002 WGS. It is not equivalent to exact dupC recovery. HG002 WGS is a concrete case in which exact known sequence reconstruction and an inappropriate mutation interpretation coexist; this can be investigated without manufacturing a global negative clinical truth label for HG002.

## Input geometry, candidate length and support units

All nine PCR libraries contain both sequence orientations and at least1022 reads with both exact conserved anchors; inputs are short primer-bounded products. Numerous fragments persist, including HG004 median422bp, MP4 median394bp and MP5 median577bp. Exact-anchor measurements are conservative because errors or biological boundary variation can break matching. No outcome-based read filter was applied.

HG004 independently demonstrates a long physical-span population:4232 exact-anchor spans round to2640bp,755 to4620bp and143 to4680bp. These are noisy observational bins, not exact haplotypes. The selected major cluster spans contigs19–72 and chooses35; the smaller cluster1–17 chooses1. The long population is represented within a broad major cluster rather than preserved as a separate selected allele. Short fragments plausibly contribute to the tiny candidate, but causation has not been experimentally demonstrated.

HG002 PCR independently demonstrates the practical consequence: one3900bp Q100 haplotype is recovered, while its other candidate is only900bp instead of4638bp. This mismatch is established against actual independent sequence truth, unlike other samples' comparisons with paper repeat counts.

The current ONT `alleles.refine_peak_contig` uses mean CIGAR indel length; the historical #20 issue's alignment-score-maximization explanation is not the current algorithm. Support units remain alignment records, including secondary/supplementary records, not unique molecules. HG003 reports47980+14680 selected records from12485 input reads, and metadata correctly labels `support_basis=alignment_records`, `molecule_count=null`. Such counts cannot measure molecular reconstruction concordance.

## Calling versus downstream interpretation: MP1–MP4

MP1 filtered allele2 VCF contains G>GC at ladder position1512, QUAL35.58, genotype1/1, DP352, AF0.6165. Its reference anchor is immediately followed by seven Cs. With the500bp ladder flank, this lies in total repeat17, using position52 as the VCF anchor. Sliding the inserted C within the homopolymer is compatible with the paper's repeat17/nucleotide59 representation. Yet the final classified mutation is an insertion A at repeat18, no named dupC, `vcf_support=false`, `vcf_support_status=absent`. Thus raw C-insertion evidence is present; final exact event recovery fails downstream of calling, and a generic PATHOGENIC banner cannot rescue it.

MP3 filtered allele1 VCF similarly contains G>GC at position912, immediately before seven Cs, compatible with total repeat7's reported C-tract event. QUAL21.05, genotype0/1, DP18273, AF0.2794. Candidate `consensus_haplotype=1` selects the reference genotype allele at this heterozygous site, so that consensus does not contain the alternate C insertion. Both candidates are unphased. A different, very long133-unit candidate contains an ambiguous dupC template label around repeat80, with `vcf_support=false`/`localization_ambiguous`; its supporting insertion VCF has DP2. This is not defensible exact event recovery or position recovery. Large DP is not necessarily molecular support when alignment records are repeated.

MP2/MP4 selected short second candidates instead of the published long allele observations and contain only substitution records in filtered VCFs. Their final ambiguity prevents callable recovery. The available evidence cannot isolate a model sensitivity defect from the upstream loss or mixing of the correct allele representation.

These are *caller-stage observations*, not independent orthogonal confirmations of the reported patients. Final endpoint recovery remains0/4; raw VCF compatibility is a separate diagnostic finding and must not be silently substituted into the predeclared score.

## Independent Q100 comparison: literal versus evidence-qualified recovery

Source sequences are versioned Q100 v1.1, exact boundaries and minus-genomic/coding-forward normalization documented in `wave2-source-research.md`. Paternal interval yields3900bp; maternal4638bp preserves the biological18bp motif expansion. Allele assignment allows permutation and preserves all sequence differences.

HG002 PCR: literal paternal3900bp sequence is exact (0 edit distance), while its900bp second candidate differs from maternal4638bp by3738 edit operations under the scorer's assignment. Literal exact alleles1/2, pair not exact. The scorer also reports sorted length errors[-3000,-738], a different pairing convention from sequence-optimal assignment; do not conflate those with the assignment-specific [0,-3738] errors. Independent callable-exact credit is0 because the run is uncallable.

HG002 WGS: literal paternal and maternal sequences both exactly equal independently sourced Q100 after boundary normalization, total edit distance0, length errors[0,0], literal exact alleles2/2. Scorer `sequence_equal=false`/independent_exact_alleles0 results from the *run-level callability gate*, not a DNA mismatch. Preserve the strong literal sequence result alongside the failed evidence-qualified endpoint. The per-assignment entries say independent_recovery=true because separate candidate provenance exists; this is different from the run-gated aggregate credit.

The WGS maternal VCF contains an18bp insertion, QUAL28.85, genotype1/1, DP17, AF0.8824. Final classifier interprets a boundary-associated18bp insertion with `frameshift=false`, `localization_status=ambiguous`, no exact VCF support. This makes reconstruction ambiguous and produces a PATHOGENIC banner although the complete sequence exactly matches Q100. The defect is therefore interpretation/localization/display for this result, not failure to reconstruct those bases.

## Consensus certainty and clinical display

All candidate metadata states reference confidence unverified and sequence identity unresolved. Complete segmentation and classification_coverage1.0 often coexist with those limits, and uncertain dictionary motif labels such as `?C`, `?D`, `?A` remain. Zero ambiguous consensus bases does not establish observed molecular support across a reconstructed allele.

`consensus.build_consensus` applies VCF genotypes to ladder reference with `bcftools consensus -H ... -M N`. Missing genotypes represented in VCF can be masked; bases absent from the VCF remain reference. An empty VCF, as in the patient WGS second candidate, can therefore yield a complete reference-derived sequence without independently verified sequence support. Sequenced PCR products do not contain the long artificial reference flanks, so successful trimming cannot verify those flanks were observed.

`report.compute_clinical_decision` enters PATHOGENIC whenever any mutation dictionary is listed, without requiring a named pathogenic event, frameshift=true, exact VCF support or resolved localization. HG002 WGS demonstrates that an in-frame ambiguous boundary change suffices. With no listed mutation, sufficient aggregate support, few ambiguous bases and completed execution, it emits NEGATIVE despite unphased candidates or unverified reference confidence. Actual displayed labels remain part of this frozen baseline; neither branch constitutes benchmark truth or independently established clinical accuracy.

MP5 has one generic C insertion at classified repeat35/base20 with exact sequence/VCF concordance, QUAL79.23, DP402, AF0.8333. The paper describes a comparator-only C insertion at repeat35/base23. Homopolymer representation may shift the nucleotide index, and selected long length69 differs from comparator83. Without independent confirmation or full sequence truth, describe this as exploratory insertion evidence, not a proven diagnostic recovery or positional validation.

## WGS input and identity limitations

HG002 WGS has54 ultralong locus-enriched reads,26 with both exact anchors. Patient WGS has2152 reads, median930bp and only5 exact full spans; it is a materially different input mixture. Neither is a full-genome end-to-end benchmark. The paper names MP4 WGS, ENA alias names MP1 WGS, and Supplement5 has further mismatched WGS rows. Preserve unresolved identity; do not assign the patient WGS empty event output to MP1's recovery denominator.

## Evidence-driven next-wave recommendations

- [#20](https://github.com/berntpopp/MucOneSpan/issues/20): test current CIGAR-indel refinement, broad cluster separation and explicit total/canonical conventions. HG002/HG004 long-allele losses differ from the historical >80-canonical-repeat minus-two claim. Require actual Q100 or independently documented physical spans before calling a length error.
- [#21](https://github.com/berntpopp/MucOneSpan/issues/21): separate candidate loss, genotype selection, indel calling, motif segmentation and exact event validation. MP1/MP3 provide raw C-tract insertions that fail downstream, so indiscriminate lowering of variant thresholds is not evidence-based. MP2/MP4 need upstream allele representation examined before attributing missed events to Clair3.
- [#46](https://github.com/berntpopp/MucOneSpan/issues/46): evaluate unique-read assignment and error-aware molecular concordance, reporting raw versus spanning versus assigned denominators, competing sequences and reference-only regions. Full segmentation and high alignment counts are insufficient. Avoid treating an exact raw-read match requirement as a general ONT accuracy standard.
- [#47](https://github.com/berntpopp/MucOneSpan/issues/47): use HG004's independently visible long span population and HG002's Q100 mismatch as preserved test cases. Evaluate physical span evidence before choosing a tie-breaker; preserve boundaries, short-fragment accounting, ONT errors and in-frame motif expansions. Rounding every length to60bp loses biological information.
- Clinical interpretation/localization follow-up: retain the HG002 WGS exact-sequence/in-frame-PATHOGENIC case and MP1 identity-misclassification case as concrete regression requirements. Separate mutation observations from established pathogenic events and sequence matching from run-level callability without changing this frozen baseline.

No algorithm, new threshold, POA method or clinical tier change is selected by this audit. Completed cohort-v2 results independently confirm the stage findings above; the initial engineering cohort is excluded from final resource aggregates.
