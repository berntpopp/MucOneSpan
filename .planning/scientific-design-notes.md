# Scientific design notes — 2026-09-14

Read-only investigation for the production-validation coordinator. No production source or tests were changed. Sources inspected: root/package/test AGENTS, developer/contributor documentation, classify/repeat_alignment/classify_types/calling/consensus/vcf/alleles, relevant unit tests, and the classification, mapping, consensus and validation audits.

## Bounded repair contracts

### 1. Signed frame accounting

- In `classify_repeat`, define downstream frame change as `_compute_net_indel(differences) % 3 != 0`, not inserted-plus-deleted length modulo three. Record the signed net value separately from total indel burden. Equal insertion/deletion lengths restore the downstream frame; this does not establish absence of an altered translated segment between events.
- For exact mutation templates, calculate the frame change from template length minus its actual parent repeat length. Avoid the current blanket `frameshift=True` assignment; current known templates may all pass that condition, but the contract should remain true for future in-frame templates.
- Document existing `RepeatDifference.pos`: insertion is before the 1-based reference position, with `len(ref)+1` for insertion at the end. The current comment claiming insertion is after that base is wrong. Preserve the existing public coordinate values and traceback tie handling.
- Regression sequences from the audit: `x[:5] + x[7:30] + 'A' + x[30:]` is net -1 and frameshifts; a one-base insertion plus a one-base deletion has net zero; net ±3 does not change the downstream frame. Test difference positions at both ends and deletion endpoints.

### 2. Honest segmentation and completeness

Current `_forward_classify` consumes every best window regardless of edit distance, stopping normally only with <30 bases remaining. Its backward fallback requires >30 bases remaining and therefore cannot normally repair poor forward segmentation. `_apply_bidirectional_fallback` also forces an arbitrarily long gap into a single mutation and does not preserve all backward mutation annotations.

A bounded repair should preserve the fast exact-template path and expose unresolved sequence instead of inventing repeat/mutation interpretations:

- Record 0-based half-open input sequence spans for every accepted repeat and every unresolved interval, plus total/classified/unresolved bases and complete/partial/empty segmentation status.
- Stop accepting forward repeat candidates when the explicitly documented fit rule fails or the candidate contains ambiguous/non-DNA bases. Reusing the backward helper's existing edit-distance limit of three is a plausible conservative candidate rule, but it is a behavior change requiring fixture and matrix evaluation, not a validated universal threshold.
- Recover an exact or accepted suffix from the 3' end if possible. Keep the intervening gap as unresolved sequence; never classify the entire gap as one arbitrary dictionary mutation. Do not assert absolute biological repeat indices for suffix units after a gap of unknown repeat count. Ordered output ordinal and resolved repeat index are different quantities.
- Include every trailing residue, including 1–29 bases, in unresolved coverage. A five-X sequence with 29 Ns must be partial, with 300/329 classified bases. Entirely ambiguous input must not produce mutation calls or a complete structure.
- Preserve legacy `allele_confidence` and `exact_match_pct` if compatibility requires, but label them as dictionary fit among classified windows. Add an explicit coverage fraction and full-sequence status; do not describe the legacy mean as probability of correct reconstruction. If a coverage-adjusted score is added, call it a heuristic and document its formula.
- A mutation from an unresolved segment, ambiguous bases, or an index-ambiguous suffix cannot qualify as an exact localized mutation call.

Tests should include exact dictionary/template controls, trailing residues 1/29/30, internal Ns, an insertion exceeding the probe bound, an unresolved middle with a terminal-repeat suffix, and a suffix mutation whose absolute repeat index cannot be established. Assert complete partition of input spans with no overlap/loss and no extra confident mutations. Existing bidirectional unit tests only check endpoint labels and existence of a mutation; they do not establish that fallback ran or that segmentation is correct.

A fresh local probe using the bundled X sequence with 40 inserted As yielded two poor windows (edit distances 26 and 22) followed by an exact but wrong terminal type. This supports conservative stopping for that particular failure; it is not evidence that a threshold of three solves every segmentation ambiguity.

### 3. Exact VCF event support and coordinate projection

The existing `vcf_support` is proximity-only. `parse_vcf_variants` discards chrom/ref/alt/genotype, and a repeat-wide 90-base window supports a mutation using a SNP in the next repeat. Narrowing this window alone cannot establish mutation identity.

Required interfaces:

1. Parse selected-sample records with CHROM, 1-based POS, REF, each ALT, QUAL, GT, phase separator and phase-set information when present. Preserve missing QUAL as missing, not a discarded variant. Tool errors and malformed records must fail visibly or produce an explicit error state, never an empty successful variant set. Multi-sample files require a selected sample; bracketed `%GT` without a sample delimiter currently concatenates samples.
2. Persist the actual single-contig reference identity/sequence, selected sample and consensus allele policy, and the actual trimmed consensus interval. A nominal flank length is insufficient when flank indels or anchor trimming move the boundary.
3. Represent applied reference edits and their consensus spans from the exact selected alleles. Validate REF against the reference and replay the edits to check agreement with emitted consensus before claiming a coordinate map. A bcftools chain can provide an alternative mapping, but coordinate conversion must still handle insertions/deletions and actual sample/haplotype policy. Ambiguous consensus bases, overlapping incompatible edits, unsupported symbolic records or inconsistent replay give projection-unresolved status.
4. Match the mutation's actual sequence event against the reference-derived event in the projected region, with normalized equivalent indel representation. Require matching chromosome and selected ALT/genotype. Compare complete changed sequence where neighboring SNPs/replacements alter anchors; parent dictionary sequence and actual ladder reference are not interchangeable.
5. `vcf_support=True` only for an exact, fully attributable event. Distinguish exact, absent, insufficient variant identity, and projection-unresolved. Legacy callers providing only pos/qual may retain those values as proximity evidence, but cannot receive exact support. No reference context means no exact reference-event claim.

Useful narrow first implementation: exact normalized simple indels with verified flanking reference correspondence and applied-variant replay; explicitly leave complex/ambiguous or boundary-crossing cases unresolved. This is preferable to silently claiming support for arbitrary replacements. Left normalization can move an event through a homopolymer or across a nominal repeat boundary: local haplotype equivalence, not equality of raw POS alone, determines support. Repeated motifs can make repeat index itself non-identifiable; report ambiguity rather than choosing an unsupported biological location.

Tests: correct insertion/deletion; nearby SNP; different insertion of same length; wrong chromosome; wrong ALT in multi-allelic record; GT0/0 and missing GT; heterozygous unselected ALT; equivalent left-shifted homopolymer indels; earlier VNTR insertion; flank insertion; trimmed boundary change; deletions crossing repeat boundary; multiple mutations; no reference context; missing QUAL; tool failure; malformed query; reference mismatch; ambiguous consensus/projection. Do not label VCF concordance independent experimental validation: the consensus was generated from that same VCF.

### 4. Explicit genotype, phase and consensus states

`disambiguate_same_length_alleles` currently calls no heterozygous records “truly homozygous”, and otherwise places all heterozygous ALTs together. Both are unsound. Query failure currently turns into homozygosity through `[]`.

- Propagate external errors. Successful empty VCF means no retained variants, not proof of homozygosity or coverage.
- Parse allele integers and missing alleles, rather than enumerating only GT0/0 and GT1/1 exceptions. GT2/2 is homozygous-alt; GT1/2 is heterozygous; GT0/. is incomplete; haploid GT1 has no diploid phase claim.
- Add orthogonal metadata: length evidence status, genotype/variant evidence status, phase status, consensus policy, and sequence completeness. Keep legacy booleans for compatibility but make the explicit status authoritative. A single length cluster is unresolved multiplicity; absence of retained heterozygosity is `no_heterozygosity_observed`, not biological identity.
- Generate two haplotypes only when the selected sample has complete genotypes and phase is established across all informative heterozygous loci. Phased GT bars at separate phase sets do not establish common phase. A single heterozygous locus can yield an unordered pair of candidate sequences, while cross-locus unphased variants must remain unresolved. Shared homozygous variants belong in both sequences.
- Specify sample and haplotype explicitly for phase-resolved consensus. The official bcftools manual documents `--samples` and `--haplotype` controls; preserve the input phase/ALT indexing and verify real-tool sequence outputs. Do not adopt blanket `-H A`: the audit demonstrates increased supported false extras and a cis/trans failure.
- For unresolved heterozygosity, retain an explicitly marked mixed/ambiguous consensus artifact or no-call, with no fabricated two-haplotype identity. SNP IUPAC codes do not encode unresolved indel haplotypes, so metadata remains necessary. Unequal-length clusters with residual heterozygosity also need an explicit unresolved state; length separation alone does not prove sequence purity.
- Sparse variant-only VCFs cannot establish whole-VNTR reference confidence or uncovered positions. Preserve that limitation explicitly unless verified read coverage/masking is implemented separately.
- Persist final allele state after calling, and include VCF paths/consensus policies in a manifest usable by standalone consensus. Current pre-call alleles.json can disagree with summary.json, and standalone consensus misses merged/variants.vcf.gz.

Tests: actual bcftools consensus for opposite-phase `1|0` and `0|1`, shared hom-alt, multi-allelic `1|2`, unphased two-locus cis/trans ambiguity, different PS blocks, single het, missing and haploid GT, empty successful VCF, query exception, malformed records, ambiguous indel, and standalone consumption of merged paths. Existing tests assert no-het=>homozygous and unphased het=>two alleles; update those contracts deliberately.

Primary tool documentation checked: [bcftools consensus manual](https://samtools.github.io/bcftools/bcftools.html#consensus). The repository audit's real-tool observations apply specifically to installed bcftools 1.17; current manual version alone cannot substitute for integration tests of supported installed versions.

### 5. Candidate fit versus independent molecule support

- Preserve secondary alignments for candidate-length fit; the fresh and historical primary-only ablations lost every exact length pair. Renaming simulator reads and exact-anchor defaults are not standalone accuracy fixes.
- Describe idxstats and current `reads` as alignment-record counts. Add clearly named alignment-record and independent-read evidence fields. Do not silently reinterpret the old coverage threshold as molecule support while changing counts under it.
- Independent sequence-record IDs must be collision-free before alignment. Existing duplicate QNAMEs label distinct input records; deduplicating retained BAM by QNAME loses biological observations. Without validated identity, molecule count is unknown. Duplicate input record identity is distinct from sequence-identical molecules.
- With valid IDs, count a molecule once per candidate and once in any cluster union, retaining its alternative fits. Assignment to two reported alleles must track shared/ambiguous reads explicitly instead of counting each twice as independent evidence.
- Keep fit summaries (AS, indel burden, alignment counts, candidate reference) separate from support summaries (unique validated records, spanning status, ambiguous assignment). Coverage and terminal spanning are different; a short clipped read may support local variation without supporting repeat length.
- Repair equal-depth plateau splitting so a flat profile does not become two valleys merely from contig spacing. Test equal plateaus and separated meaningful minima without introducing a truth-tuned threshold.
- Make reported inferred length and actual selected-reference length explicitly distinct when current refinement disagrees with cluster center. Do not imply the two are equal. Preserve compatibility unless a validated inference policy is introduced.
- Pool-before-threshold and anchor-based assignment are potentially valuable redesigns, but require held-out low-depth/unequal-coverage/HiFi/ONT evidence. They should not enter this repair as silent default changes.

Tests: one read with primary+secondary records; distinct records sharing original QNAME; one valid read appearing on multiple cluster contigs; molecule shared across two candidate clusters; partial reads; plateau contigs 40/43/46; missing AS; conflicting count/reference lengths. Acceptance includes transparent unavailable molecule counts and unchanged established candidate-fit behavior when only metadata is added.

## Suggested sequencing and acceptance

1. Signed accounting, visible VCF errors, record parsing and additive state/spans fields.
2. Conservative unresolved segmentation and exact-support projection with unsupported complex cases explicit.
3. Phase-resolved consensus paths plus manifests; unresolved paths remain explicit artifacts/no-calls.
4. Evidence accounting and plateau safety without changing unvalidated mapping/anchor/QUAL defaults.

Coordinator should require deterministic fixtures first, real bcftools integration for applied haplotypes/projection, then the preserved error-free dictionary and 106 truth-haplotype controls. Run the full pipeline benchmark on positives and negatives and score exact event name/parent/index, extra calls, no-calls, complete sequence/structure, and missing predictions. A reduced false-positive count obtained by declaring everything unresolved is not improved mutation sensitivity; both outcomes must be reported. No new threshold should be selected using held-out outcomes.

Validation performed in this subtask: read-only source/test/audit review; one locked-environment classification probe for the >30bp insertion. No unit or integration suites rerun and no production repair validated here.
