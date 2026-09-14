I’ll check the proposal against the implementation, frozen experimental spec, and upstream tool semantics. I’m using the academic research review skill to assess the scientific claims and supporting evidence.


Two central claims need revision: selecting reads because they score best against a secondary contig makes the proposed positive-margin check largely circular, and Clair3’s haploid flags do not implement the likelihood model described. I also confirmed that mapping already streams SAM into samtools, so the claimed performance savings need a measured baseline.


The baseline already requires a common phase set and leaves unresolved same-length evidence as a mixed candidate. Its optional read-backed phasing is disabled after a development false-positive result. I’m checking how the proposal addresses that failure, along with the existing depth filters and reference-coordinate contracts.


I found a more fundamental blocker in the saved baseline: Clair3’s pileup log shows `ModuleNotFoundError: No module named 'tensorflow'`, followed by an empty VCF, while the sample is marked completed. I’m checking the extent of this failure. These runs cannot support the proposal’s variant-calling diagnosis or paired accuracy claims.


The failure affects all 140 saved development datasets: all 216 Clair3 pileup logs contain the TensorFlow import error. I also verified that the failure-atlas script reads metric keys absent from the evaluator output, allowing sequence failures to be labelled successful. The review will treat baseline validation as a prerequisite and give concrete revisions for scoring, trimming, haploid consensus, phasing, and regression tests.


**Recommendation: major revision before implementation.** The proximal-flank correction is well motivated, and replacing alignment-record counts with read-level evidence is the right direction. However, the proposed dominance test is circular as defined, Clair3’s haploid semantics are misstated, and a common phase set cannot guarantee correct haplotypes.

More urgently, the saved development baseline is invalid for evaluating variant calling and full reconstruction.

**The baseline must be repaired before drawing architectural conclusions.**

I inspected the checkout at `9535f7ee02033cda3da9b22c5b73752e5473b016`, both specifications, relevant implementation/tests, and saved development artifacts.

All **216 Clair3 pileup logs across all 140 development datasets** contain:

```text
ModuleNotFoundError: No module named 'tensorflow'
```

For example, this [Clair3 log](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/tests/results/baseline_dev_140/ctrl_ident_60_60_hifi/allele_1/clair3/run_clair3.log:78) subsequently reports creation of an empty VCF and exits after the pileup stage. Nevertheless, all 140 pipeline status files report `completed`.

The saved evaluation reports zero exact diploid reconstructions and zero mutation true positives. These results cannot establish that diploid calling causes the observed reconstruction failures: variant inference failed to execute.

There is a second evidence problem. The [failure-atlas implementation](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/scripts/run_eval_pipeline.py:129) reads `length_correct`, `exact_sequence_concordance`, and top-level `event_tp`/`event_fp`. The evaluator instead returns bounded metrics under `metrics`. Missing keys therefore bypass checks or become zero, allowing incorrect reconstructions to be labelled successful.

Before comparing architectures:

1. Validate actual Clair3 inference with its intended interpreter and model, including a known-positive fixture.
2. Detect internal worker failures even when the outer wrapper returns success and writes a syntactically valid empty VCF. An empty VCF alone is neither failure nor evidence of successful calling.
3. Repair the failure atlas to consume the authoritative evaluator schema and conservative metric bounds.
4. Rerun the frozen baseline in the repaired environment, then run candidates under the same environment.
5. Record source provenance beyond `HEAD`: the current worktree includes a scientific change to `config.py`’s `delete_insert` semantics. Its involvement in earlier runs needs an audit.

The existing mapping artifacts can still inform narrowly scoped alignment hypotheses after validation. They cannot validate the proposed variant-calling diagnosis or whole-pipeline runtime estimates.

**1. Family A: the slice correction is appropriate, but reference construction and trimming must change together.**

I confirmed that both the current generator and the bundled ladder use the distal left-flank prefix. This is a **flank-selection error**, not a sequence inversion.

With validated equal flank lengths, switching from prefix to suffix does not itself shift the synthetic reference’s repeat coordinates:

\[
b_L=W,\qquad b_R=|R|-W.
\]

These are boundaries on the **reference**, however—not necessarily on the consensus after applying indels.

The principal integration risks are concrete:

- **Changing `build_contig()` does not change the default runtime reference.** The pipeline loads the bundled FASTA. Regenerate and validate that resource, or explicitly select a newly generated reference.
- **The existing left trimming anchor remains distal.** [consensus.py](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/src/muc_one_span/consensus.py:169) constructs it from `flanking_left[:flank_length]`. Updating the ladder alone will cause anchor failures and silent fixed-position fallback.
- **Old and new references have identical contig names and lengths.** Names and lengths cannot detect an incompatible BAM, VCF, index, or trimming policy. Use sequence hashes and an explicit layout manifest.
- **Alignment outcomes will change.** Seeding, clipping, chaining, primary/secondary selection, MAPQ, and AS distributions can all change. Recalibrate Family C after Family A.
- **Custom dictionaries need an orientation contract.** A suffix is proximal only when the stored left flank has the declared orientation and endpoint.

Prefer one shared reference-layout function that returns the selected flank sequences, fixed-repeat layout, reference boundaries, and anchor components. Generation, trimming, BED generation, and provenance should consume that same object.

The proposed `ConsensusSettings.proximal_flank` follows an existing coupling, but flank orientation is fundamentally a **reference-layout property**. If retained there for compatibility, record it in reference provenance and reject mismatches.

For trimming, distinguish three independent quantities:

| Quantity | Required interpretation |
|---|---|
| Positional tolerance | Allowed displacement from a predicted boundary |
| Sequence tolerance | Allowed substitutions/indels in an anchor |
| Boundary uncertainty | Whether equally plausible alignments imply different cuts |

The current `anchor_tolerance=50` controls positional displacement; matching is still exact. `_find_anchor()` returns the first match, without testing uniqueness.

My recommended policy is:

1. Project reference boundaries through the actual applied edits or a validated reference-to-consensus alignment.
2. Verify each projected boundary with a junction-spanning anchor.
3. Examine all competitive anchor placements; require a unique boundary or agreement among tied placements.
4. Require \(0\le b_L<b_R\le |C|\), correct orientation/order, and preservation of the configured fixed-repeat blocks.
5. If boundary evidence is unresolved, preserve a candidate with explicit uncertainty rather than treating fixed trimming as verified.

For an initial **tuning experiment**, the existing 20 flank bases plus 20 fixed-repeat bases is a reasonable starting anchor. Test exact matching first, then bounded approximate matching—such as one or two edits per 40-base anchor—with explicit uniqueness checks. A narrow search around a correctly projected boundary is preferable to broadly increasing the current ±50 window. These are candidate settings, not validated defaults.

Crucially, a genuine mutation overlapping an anchor must not be “corrected” into dictionary sequence. If an indel crosses the boundary, the specification must define sequence ownership or report ambiguity. Likewise, final VNTR length must not be forced to a multiple of 60: the benchmark deliberately contains length-changing mutations.

**2. Family C: the proposed test does not establish a second allele.**

Under the proposal’s definition,

\[
R_2=\{r:c^*(r)=c_{(2)}\},
\]

every read in \(R_2\) necessarily satisfies

\[
AS(r,c_{(2)})-AS(r,c_{(1)})\ge0.
\]

With unique maxima, all satisfy a strictly positive margin. Thus the 75% positive-margin criterion is largely guaranteed by the selection procedure itself.

The claimed homozygous behavior—reads assigned by maximum AS to a neighboring contig nevertheless having negative margins against the primary contig—contradicts that definition. It could happen if candidate membership came from a different procedure, but that procedure would need to be specified.

The average-margin rule is also weaker than the review request suggests. For example:

\[
\Delta AS=(1,1,1,1000),\qquad \delta=100
\]

passes both the mean threshold and 75% positivity, although only one read exceeds the meaningful margin.

If retaining a margin heuristic, at minimum define:

\[
D_2=\{r:S(r,c_2)-S(r,c_1)>\delta_r\}
\]

and require a minimum number of distinct eligible observations in \(D_2\), alongside support coherence and an independently calibrated false-split criterion. Do not substitute average margin for the number of convincing reads.

Even this is not sufficient: selecting the best secondary candidate among many lengths creates a multiple-search problem, and PCR stutter can produce coherent nonconstitutional length modes.

**Minimap2 does not provide the assumed complete score matrix.**

`-N 10` limits reported secondary alignments; `-p 0.8` suppresses lower-scoring secondary chains and can prevent their subsequent base alignment. Neither ensures that a particular read has scores against both candidate lengths. A missing alignment is missing evidence, not AS=0 or negative infinity. [Minimap2 manual](https://lh3.github.io/minimap2/minimap2.html)

The extreme-asymmetry example is therefore not a valid numerical justification:

- A long read against a shorter reference need not receive a full-length alignment with a large gap penalty; clipping or splitting can intervene.
- Relative to the shorter reference, the excess query sequence is an insertion, not a deletion.
- A large AS difference is not a calibrated probability, much less “100% confidence.”
- Total repeat count 140 corresponds to `contig_131` only for the default nine fixed repeats.

Use minimap2 for candidate discovery, then explicitly obtain comparable scores for the same read interval against each shortlisted candidate. Anchored alignment across the entire VNTR is preferable to comparing local alignments covering different portions of a read. Missing or failed constrained alignment must remain an explicit state.

**A stronger design separates candidate discovery, model selection, and assignment.**

I recommend:

- Discover candidate lengths from anchored read spans and alignment evidence, allowing neighboring lengths to share support.
- Compare one-allele and two-allele models using all eligible reads.
- Assign reads only after deciding which models remain credible.
- Retain an artifact/outlier component for stutter, chimeras, and badly aligned reads.

For example, with a calibrated platform-specific likelihood:

\[
\ell_1=\max_c\sum_i\log p_t(r_i\mid c),
\]

\[
\ell_2=\max_{c_1,c_2,\pi}
\sum_i\log\left[
(1-\pi)p_t(r_i\mid c_1)+\pi p_t(r_i\mid c_2)
\right].
\]

Accept a second candidate using a calibrated improvement threshold plus minimum attributable support. Calibration must repeat the **entire candidate search** under the one-allele null. A standard chi-square likelihood-ratio threshold is inappropriate without justification because mixture identification fails under the null.

Do not silently equate \(\exp(AS)\) with a likelihood. A simpler validated score-based decision can be preferable to an elaborate but uncalibrated probabilistic model.

Also, candidate discovery cannot require an exact-contig peak to exceed \(m_{\min}\): five genuine minority reads distributed across three adjacent contigs may fail that rule despite coherent aggregate support.

**Platform thresholds should reflect measured errors and the decision being made.**

There is no defensible universal numerical AS margin available from the proposal. Raw scores depend on alignment parameters, read length, sequence composition, clipping, and candidate separation. `map-hifi` and `lr:hq` also use different scoring parameters. [Minimap2 v2.28 options](https://raw.githubusercontent.com/lh3/minimap2/v2.28/options.c)

Calibrate margins against homozygous controls, stratified by platform, length, and quality. A possible standardized statistic is:

\[
Z_i=
\frac{\Delta S_i-\mu_{0,t}(L_i,Q_i,\Delta L)}
{\sigma_{0,t}(L_i,Q_i,\Delta L)}.
\]

Estimate its null distribution empirically, including systematic errors.

The proposed floors of three HiFi and four ONT reads are reasonable **experimental discovery floors**, not sufficient evidence for exact sequence reconstruction or confident mutation exclusion. Define separate requirements for:

- length detection;
- allele assignment;
- callable sequence coverage;
- mutation evidence;
- confident negative reporting.

For ONT, stronger requirements for context consistency or effective independent support may be warranted. Simply raising the AS threshold could eliminate the intended minority rescue.

For perspective, under a hypothetical independent binomial sampling model with 100 reads and 5% minority abundance, the probabilities of observing at least three, four, and five minority reads are approximately 88.2%, 74.2%, and 56.4%. The frozen simulator’s allocation must be evaluated separately; this calculation illustrates why unconditional recovery cannot be guaranteed.

Finally, **read records are not necessarily independent molecules**. Preserve distinct source observations despite QNAME collisions, avoid counting secondary/supplementary copies as support, and retain molecule-level dependence information where available.

**Ambiguous support should yield partial evidence, not a forced genotype.**

Recommended fallback behavior:

- Preserve a supported major sequence.
- Retain plausible secondary lengths and their evidence without promoting them to reconstructed alleles.
- Report unresolved multiplicity when a second allele is not established.
- Distinguish length ambiguity, low coverage, assignment ambiguity, and sequence uncertainty.
- Avoid interpreting missing minority evidence as a confident negative.

Equal-length but compositionally different alleles must remain eligible for sequence-based resolution even when length inference finds only one mode.

**3. Family D: haploid calling needs verified partition purity and explicit output support.**

Yes—cross-mapping can create both false alternate calls and missed true variants.

The proposal’s likelihood description is incorrect for the installed Clair3 v1.0.10. Its implementation first interprets variant/genotype predictions, then:

- `--haploid_precise` discards heterozygous predictions;
- `--haploid_sensitive` retains eligible heterozygous predictions;
- haploid output genotypes are rewritten as single allele indices.

This is visible in the installed [CallVariants.py](/home/bernt-popp/miniforge3/envs/env_clair3/bin/clair3/CallVariants.py:1161), and agrees with the documented flag descriptions. [Clair3 documentation](https://github.com/HKU-BAL/Clair3)

Consequently:

- Contamination that produces a heterozygous prediction can become an alternate haploid call in sensitive mode.
- A genuine allele-specific variant predicted as heterozygous because of contamination or model uncertainty can disappear in precise mode.
- Removing IUPAC bases can conceal uncertainty without improving correctness.

Do not choose `precise` versus `sensitive` solely by platform. Validate both against controlled contamination, allele imbalance, and sequence-error conditions.

There is also an immediate compatibility defect: existing [phase_evidence()](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/src/muc_one_span/phasing.py:28) labels single-index genotypes `non_diploid`, while [variant_support.py](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/src/muc_one_span/variant_support.py:38) requires two genotype indices for consensus replay. Adding Clair3 flags alone would break downstream evidence interpretation.

Introduce an explicit haploid consensus/evidence path. Do not represent a haploid observation as evidence for two biological copies.

Before haploid calling, require:

- sufficient support across the full retained sequence;
- credible assignment against competing references;
- checks for residual conflicting variants and contamination;
- retained assignment provenance;
- an unresolved outcome where purity cannot be established.

For close-length cases, assignment can be circular at individual variants: a read carrying an alternate base may preferentially map to the reference containing that same base, inflating the allele fraction later used to validate it. Where practical, test support after masking the candidate site during assignment or using other linked markers.

The current `min_dp` parameter is explicitly **not applied** by `filter_vcf()`. Therefore, the new architecture needs an actual selected-sample depth and callability policy; passing an existing argument will not provide one.

**A single phase set is necessary for unrestricted cross-site haplotype emission, but is not sufficient for correctness.**

A phase set describes a claimed relationship among retained variants. It does not prove correct genotypes, correct repeat placement, two observed haplotypes, or absence of switches. WhatsHap explicitly supports measuring switch errors within phased results. [WhatsHap guide](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-compare-comparing-variant-files)

For example, true haplotypes `00` and `11` can be reconstructed as `01` and `10` within one PS if erroneous or chimeric reads support the wrong link. Connectivity is satisfied; the result is recombinant.

The proposed invariant already exists in the baseline. Optional read-backed phasing is disabled because earlier development validation introduced an extra supported mutation on a collapsed cluster; this is documented in [limitations.md](/home/bernt-popp/development/MucOneSpan/.worktrees/simulation-200/docs/reference/limitations.md:113). The proposal must explain what new evidence requirement addresses that failure.

Strengthen the invariant with:

- multiple independent observations supporting critical phase links;
- explicit support for both haplotypes, not merely one haplotype and its inferred complement;
- rejection or uncertainty for weak bridges and conflicting linkage;
- stability checks after removing a bridge read or resampling observations;
- preserved variant identity and correct multiallelic genotype indices;
- coverage of unresolved sites in the evidence model, rather than making connectivity easier by filtering them away.

Do not require a PS for a single heterozygous site: it supports an unordered pair without cross-site phase.

Refusing to stitch disconnected blocks remains correct. Preserve local block evidence, however. IUPAC can represent SNP ambiguity; it cannot faithfully encode alternative indel haplotypes or every possible relationship among disconnected blocks.

**4. Family B needs explicit protection against reference-driven self-confirmation.**

Two-pass refinement is promising, but AF≥0.70 is insufficient as its principal safety condition.

First, pooled same-length alleles often have genuine heterozygous differences near AF 0.5. Applying only ≥0.70 variants to both references can make them identical, after which competitive alignment assigns nothing. The flowchart’s same-length branch must first establish phased or otherwise supported alternative sequences.

Second, errors incorporated into a reference can recruit the same erroneous reads during pass two. A dominance margin does not independently validate an error that created that margin.

Third, **incorporating a mutation into the tailored reference can erase it from the second-pass VCF**. A now-reference allele may yield no variant record, while current mutation support relies on replaying the VCF against its reference and reverting the event. Without retaining pass-one edits and coordinate mappings, sensitivity or support annotation can regress despite an improved consensus.

Require:

- an immutable original backbone and complete edit provenance;
- composed coordinate mappings across both passes;
- references built from supported haplotypes rather than independent per-site majority choices;
- explicit handling of identical-reference convergence;
- re-estimation of length after substantial indel correction;
- preservation of uncertain reads and coverage masks;
- a defined rollback/acceptance rule when refinement loses evidence.

Internal VCF replay establishes consistency, not independent support for the same variants used to construct the reference.

**5. Families E and F need clearer semantics and measured benefits.**

The ONT proposal’s “chaining gap open 10, gap extend 1” is underspecified. Minimap2’s alignment gap parameters and chaining parameters are different controls. State the exact flags, their ordering after the preset, and the effective configuration. Define ONT chemistry/basecaller/model applicability rather than treating all simplex and duplex data as one error distribution.

QUAL 8 and 6 are tuning candidates, not justified precision guarantees. Evaluate mutation-type-specific sensitivity and control false alarms, especially homopolymer indels.

For performance:

- Initial mapping already streams minimap2 into samtools sort. Allele remapping still materializes FASTQ and SAM, so optimization opportunities exist there.
- Minimap2 v2.28 loads `.mmi` using allocated structures and `fread`; the proposed `mmap` description is inaccurate. The stated startup times require measurement. [Index implementation](https://raw.githubusercontent.com/lh3/minimap2/v2.28/index.c)
- Cache identity must include reference content and effective indexing parameters. Index settings embedded in `.mmi` can override command-line indexing choices. [Minimap2 manual](https://lh3.github.io/minimap2/minimap2.html)
- Use a writable cache, atomic publication, concurrent-build handling, and invalidation. Do not assume installed package resources are writable.
- BED coordinates must come from actual reference boundaries, with explicit padding and 0-based half-open semantics. Hard-coded `450` is incompatible with custom flanks and refined references.
- Excluding flanking variants can remove useful phase links. Interval restriction therefore cannot automatically be classified as scientifically equivalent performance work.
- Budget threads across concurrently running mapping, sorting, calling, and allele jobs.

Streaming also changes a scientific contract if implemented carelessly. `samtools fastq` excludes secondary/supplementary alignments by default and chooses one sequence per QNAME/category. Introducing collation can therefore collapse distinct records with colliding names. Recover source reads through stable observation identities and explicit membership, not merely contig-region extraction. [Samtools 1.21 FASTQ semantics](https://www.htslib.org/doc/1.21/samtools-fasta.html)

Implement pipelines with argument lists through the tool abstraction, drain stderr safely, check every process, and preserve stage failure status.

**6. The module plan should follow evidence responsibilities, not Family labels.**

The projected file counts are not an adequate implementation plan: `alleles.py` currently has 596 lines and `calling.py` 446. Existing behavior cannot simply disappear to meet projected sizes.

A cohesive breakdown would be:

| Component | Responsibility |
|---|---|
| `ladder.py` | Pure reference construction |
| `reference_layout.py` | Shared boundaries, flank selection, layout manifest |
| `alignment_evidence.py` | Observation identities, score completeness, anchored evidence |
| `length_inference.py` | Candidate discovery and one/two-allele decisions |
| `read_assignment.py` | Competitive assignment, ambiguity, purity evidence |
| `alleles.py` | Compatible public interface and result serialization |
| `calling.py` | Calling orchestration |
| `phasing.py`, `read_phasing.py` | Existing phase policy and tool integration |
| `reference_refinement.py` | Tailored references and edit-coordinate provenance |
| `consensus.py`, `trimming.py` | Consensus construction and boundary validation |
| `reference_cache.py`, tool pipeline helper | Index lifecycle and subprocess execution |

Split only where responsibility warrants it; these need not all become large modules. Aim below roughly 450–500 lines to leave maintenance room. Preserve public imports and test patch targets. Keep typed immutable settings, validation, and configuration precedence consistent across CLI and library entry points.

**Required tests should target counterexamples and scientific contracts.**

No unit-test suite can guarantee zero regression on sequencing datasets. Unit tests can establish deterministic invariants; actual-tool integration and paired benchmark execution must establish observed regression behavior.

| Area | Essential tests |
|---|---|
| Execution integrity | Internal Clair3 worker failure plus outer success/empty VCF must become execution failure; genuine successful empty VCF remains valid |
| Evaluation | Authoritative nested metric bounds; missing fields fail explicitly; incorrect sequence never receives “successful reconstruction” diagnosis |
| Reference migration | Distinct prefix/suffix flanks; bundled reference compatibility; stale index/reference hash rejection |
| Trimming | Zero/short/custom flanks; indels before each boundary; mutated anchors; multiple hits; tied alignments with different cuts; crossing/reversed boundaries |
| Dominance | Circular-selection counterexample; mean dominated by one outlier; exact ties; missing AS; missing competitor; clipped versus fully spanning alignment |
| Support identity | Secondary copies add no support; distinct colliding QNAMEs survive; one observation cannot enter two exclusive partitions |
| Homozygous controls | Noisy neighboring-length preferences; stutter modes; identical tailored references; unresolved multiplicity without fabricated duplication |
| Minority alleles | Support spanning 0–8 reads; both long and short minority direction; support spread across adjacent candidates; detection distinct from sequence callability |
| Close lengths | Gaps 1/2/3 under composition variation, imbalance, and boundary mutations |
| Haploid calling | Precise/sensitive semantics; single-index GT parsing, consensus, and replay; contamination causing both false positives and false negatives |
| Phasing | Cis/trans, disconnected blocks, missing PS, one erroneous bridge, one supported haplotype only, conflicting records, multiallelic sites |
| Refinement | AF≈0.5 same-length variants; incorporated mutation remains traceable; rollback; duplicate convergence; insertion/deletion coordinate composition |
| Coverage | Uncovered reference sequence never becomes verified wild type merely because the VCF is empty |
| Performance | Cached/uncached and streamed/materialized scientific equivalence; cache races; upstream/downstream process failures |
| Compatibility | Existing JSON fields, multiplicity semantics, explicit presets, configured layouts, wheel resources, mutation conventions |

Retain existing passing fixtures and explicitly review expected-failure cases when they start passing. Run `make ci-check`, applicable integration checks, and `make build-check` for the bundled-reference change. Aggregate 80% branch-aware coverage is a floor, not proof that these decision boundaries were exercised.

**The experimental acceptance rules also need operational clarification.**

The architecture’s verification section omits several frozen gates, including annotation, confident-negative rate, challenging strata, memory, runtime, and protected validation.

Additional issues:

- Eight net HiFi wins among 70 samples equal **11.43 percentage points**. The separate +12-point requirement requires at least nine net wins.
- A net-win threshold does not establish statistical significance. Five wins and zero losses yield a two-sided exact paired sign/McNemar probability of **0.0625**. Predeclare the test, sidedness, and multiplicity treatment.
- Use exact **sample-level diploid reconstruction** for the primary endpoint, not aggregate allele-level `sequence_accuracy`.
- Keep platform pairing at the biological-design level when reporting combined uncertainty.
- Define FDR behavior when no mutations are called; the present broken baseline has no calls, so its FDR cannot meaningfully establish a decrease.
- Separate confident mutation-negative reporting from independently reconstructed diploid identity. A valid single-sequence negative control with unresolved multiplicity should not require fabricated duplication.
- The supplied runner loads the sealed ledger and uses biological design names. Separate blinded execution from truth-aware evaluation before protected validation.
- The frozen seed text says stride ≥1,000, but its formulas use `design_id * 10`. Resolve the inconsistency and verify actual derived seed collisions.
- Assess Family A independently, then C/D/B and selected interactions. Family F should undergo equivalence checks against an unchanged scientific configuration.

My approval conditions are a valid rerun baseline, a noncircular second-allele decision, complete comparable alignment evidence, explicit haploid/phase uncertainty handling, and reference-to-consensus provenance through refinement.

This was a read-only review; no files were changed and no final-validation truth was opened. I performed source/log inspections and small deterministic checks. `uv run --locked --all-extras --no-sync` could not acquire its lock in the read-only environment; no CI or new sequencing benchmark was run.