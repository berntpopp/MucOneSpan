# Production reconstruction and validation specification

Status: specified before implementation; pre-freeze decisions recorded below.
Baseline d8390b3c244ef8f3240af74b92db12b50dfc77d1.
Date: 2026-09-14. The user's full six-workstream request is the authority.

## Purpose and scope

Turn experimentally verified defects into production changes while testing the
scientific hypotheses behind more ambitious reconstruction methods. Supported
inputs remain PacBio HiFi and ONT amplicon BAM/FASTQ, with existing CLI flags,
reference resources, defaults and version source. This work does not establish
clinical diagnostic performance. No merge, push, publication or release.

Three approaches considered:

1. Globally lower QUAL and force ALT consensus: rejected as a default because
   the exposed ONT normal becomes positive and trans genotypes become cis.
2. Replace ladder inference with exact anchor modes: rejected as a default
   because adjacent lengths are suppressed and most ONT reads fail anchors.
3. Independent exact scoring, explicit evidence contracts, conservative phase
   handling, and gated read-evidence experiments: selected. Each change can be
   reviewed and reverted independently. Scientifically failed experiments remain
   visible; passing a software test does not promote an unvalidated method.

Non-goals: inventing a probability calibration, declaring genotype identity from
absence of VCF records, changing dictionary biological definitions, or masking
accuracy regressions with missing outputs. Broad model retraining and clinical
validation require a different dataset and study.

## Scientific contracts

- `length` is total repeat units, including five pre- and four after-repeat
  units. `canonical_repeats` and `contig_N` refer to variable ladder units;
  total = canonical + 9. Neither names nor nominal 60-base multiplication are
  authoritative truth for mutated lengths.
- Repeat indices are 1-based including pre/after units. Dictionary insertion
  position p inserts BEFORE reference base p; deletion endpoints are inclusive.
  Python intervals are 0-based half-open. VCF POS is 1-based on the actual
  ladder contig, with REF-anchored indels. Consensus spans and VCF positions
  must never be compared without a defined projection.
- Net frameshift means `(inserted_bases - deleted_bases) % 3 != 0` over the
  repeat alignment. It says nothing about transient within-repeat translation.
- Literal Levenshtein substitution costs one even for ambiguity codes; no
  IUPAC wildcard matching or case normalization in the scoring replacement.
  Existing first-dictionary-match ties and traceback order remain unchanged.
- Alignment records describe candidate-reference fit. Primary input records
  approximate independent observed molecules; QNAME is not unique in existing
  simulations. Do not deduplicate distinct records by QNAME. Actual PCR
  molecular independence cannot be inferred without source/UMI evidence.
- Equal lengths, a single observed length, and sequence-identical diploid
  alleles are separate findings. No caller may infer sequence homozygosity
  merely from a lack of detected heterozygous sites.
- Mutation detection, exact name/parent/repeat position on assigned haplotype,
  correct lengths, exact complete nucleotide sequence, and exact ordered repeat
  structure are distinct endpoints.
- Existing confidence values are heuristic similarity/evidence weights, not
  calibrated probabilities. Expose denominator/completeness independently.

## Interfaces and module boundaries

A. Offline evaluation (`evaluation*.py` plus `scripts/evaluate.py`): strict
truth/artifact adapters, pure one-to-one allele assignment and event matching,
and aggregate/report generation. Caller code must not import simulator truth.
Truth adapter validates FASTA records, exact flanks, structure reconstruction,
mutant unit sequence and counts; missing/inconsistent truth is invalid input.
Manifest identifies expected samples explicitly; it is never derived only from
successful output directories. Select BAM/FASTQ/.fq and gzip explicitly; reject
ambiguous selection unless chosen by manifest. Report nominal molecules,
retained records, platforms, seeds, simulator/config/model hashes separately.

Evaluation schema version 1 contains per-sample execution status, truth and
prediction counts, all optimal assignments or ambiguity indication, missing and
extra allele counts, signed repeat-count errors, sequence and structure exact
flags, mutation TP/FN/FP (including positives' extras), and denominators.
`execution_failed`, `not_attempted`, `invalid_truth`, `insufficient_evidence`,
`ambiguous_reconstruction`, and completed predictions are distinct. A zero
accuracy report is valid evaluation; evaluator input/execution failures return
nonzero. Failed callers keep every expected truth denominator. Never score
stale summaries after a failed invocation as successful current predictions.

Assignment minimizes summed exact global sequence distance with one-to-one
matching; unmatched alleles are explicit. Tied assignments are reported, and
haplotype-specific event success cannot be selected by looking at event truth.
Use conservative bounds for metrics affected by ties. Duplicate calls cannot
match one truth event twice. Missing wild-type alleles are reconstruction losses,
not mutation false negatives. Normal no-calls are not true negatives.

B. `repeat_alignment.edit_distance(str, str) -> int`: dependency-free integer
bit-vector exact scoring; preserve the public import from `classify`. Keep
traceback unchanged. Separate controlled timing harness compares both scorer
and entire classification dictionaries on the same inputs in alternating order.

C. Classification keeps legacy keys and adds consumed sequence spans,
`unclassified_regions`, ambiguous-base count, and explicit reconstruction
completeness. A fully consumed uncertain segmentation is not certified correct.
Repair recovery only with deterministic terminal-anchor/known-sequence evidence;
otherwise expose unresolved sequence rather than falsely exact coverage.
VCF queries retain CHROM/POS/REF/ALT/GT/QUAL and propagate tool and malformed-data
failures. `vcf_support` must mean exact compatible variant identity at a verified
reference projection; position-only legacy inputs receive unavailable/unsupported
status, with old positional proximity exposed separately if needed. A SNP cannot
support an indel. Homopolymer-normalized equivalent alleles require equivalent
local sequence on the actual reference, not raw coordinate equality. Where
projection cannot be proven, support is unresolved, not falsely validated.

D. Consensus accepts explicit sample/haplotype selection; reject ambiguous
multi-sample selection. Ordinary per-length calls preserve genotype-aware IUPAC
behavior explicitly. Equal-length shared, cis/trans and multiallelic genotypes
must select real genotype indices. Phased sites are usable only within a proven
common phase block; independently phased blocks cannot be stitched arbitrarily.
Unphased multiple heterozygous sites remain unresolved. No empty-call homozygous
inference. A read-backed phase alternative may be promoted only after known
cis/trans and missing-phase fixtures and read assignment checks. Preserve
candidate sequences where possible with evidence status; do not relabel no-calls
as exact reconstructions. Sidecar/additive metadata distinguishes phase and
sequence identity uncertainty; update persisted alleles after calling.

E. Length/assignment keeps alternative alignment fit and adds actual primary
record support. Flat valleys must not split; diffuse clusters and minority
alleles need explicit evidence. Anchor investigation uses error-tolerant matches,
reports rejected/ambiguous/partial reads and assignment purity where source truth
exists. No forced two peaks or automatic diploid duplication. Existing defaults
remain unless development and fresh acceptance criteria both pass. Unique record
IDs, if required for phasing, retain original ID mapping and every input record;
renaming itself requires a paired sensitivity analysis.

F. Local mutation rescue uses candidate-specific REF/ALT evidence at the assigned
locus and reads, with coverage and alignment ambiguity. No template-anywhere or
universal QUAL reduction. Keep HiFi and ONT results separate. Profile first;
streaming/cache/thread changes require semantic output equality and measured
wall/RSS benefit, with nominal threads distinguished from actual subprocess use.

## Acceptance criteria fixed before final data generation

Software gates: make ci-check (>=80% branch-aware coverage), make test-int,
explicit generated-data e2e selection, make docs-check; build-check for new
packaged modules/resources; security-check if dependencies/security change.
Every authored code/config/template file has <=649 physical lines. New behavior
needs a failing regression test before implementation; no lowered tolerances.

Evaluation: all hand-constructed adversarial examples score exactly; substring,
wrong position/parent, wrong phase, duplicate extras, missing allele, tied
assignment, failed sample with stale output and unknown truth cannot inflate TP.
All 53 exposed full-data samples and 24 perturbations are accounted for.

Exact scoring: zero distance mismatches on >=10,000 deterministic random and
adversarial pairs, zero complete dictionary differences on the cached development
consensus panel when isolating scoring from semantic fixes; median >=2x speedup
on difficult classification over >=3 paired runs. Exact controls must not regress
by >20% in sufficiently batched timings. Report whole-pipeline paired wall/stage
times separately, with no promised 8x pipeline speedup.

Correctness: signed-net fixtures and tool-failure fixtures all pass; every input
base is accounted for; unrelated variants never create exact support. Existing
ideal truth structures retain 106/106 exact results except additive metadata.

Phase: 100% exact sequence recovery for adequately evidenced synthetic cis/trans,
shared and multiallelic fixtures. Zero false phased assertions on unphased,
missing, conflicting or disconnected phase evidence. No negative genotype inferred
from insufficient evidence. Every withheld/unresolved output remains in accuracy
and no-call denominators.

Length/rescue promotion: >=40/44 exact original HiFi count pairs with no lost
previously exact pair, no new normal mutation-positive sample on exposed HiFi or
ONT, and no reduction of exact mutant events from baseline 24/28 and 2/2.
Fresh final comparisons require no lost previously exact sequence/count/event,
no new normal FP, and no reduced call rate for an accuracy-improvement claim.
Report failure even if an aggregate increases. At least one exact reconstruction
or exact event improvement is needed to claim scientific accuracy improvement.
No-call increases are safety-contract changes, never accuracy improvements.

## Development and fresh validation

Development: original 44 HiFi, 3 ONT, exposed seeds 9101–9106 and all existing
depth/name perturbations. The prior six-sample set is no longer held out.
Provenance limitations of original data (null seeds, misleading depth fields,
QNAME collisions, missing source-read map) must be retained in reports.

Fresh panel is generated only after settings/source hashes freeze. Reserve seeds
2026091401–2026091424, twelve designs per platform: normal 60/60, normal 60/61,
normal 60/62, dupC 60/63 H1 repeat25, dupC 60/80 H2 repeat25, insCCCC 60/80 H1
repeat55, dupA 25/30 H2 repeat10, insG 100/120 H1 repeat45, del18_31 60/80 H2
repeat25, dupC 25/140 H2 repeat100, normal 60/80 low depth, dupC 120/140 H1
repeat50. Derive deterministic record-based minority and partial-read subsets
from explicitly mapped source reads where available; otherwise assignment-truth
availability is a disclosed limitation and synthetic fixtures supply that test.
Use exact recorded simulator commands/configurations and actual usable records;
200 requested templates except low-depth design 20. Platform error models and
versions are separate strata. Do not tune on final results; if a bug fix changes
settings after exposure, mark the panel development and reserve new seeds.

Report raw sample-level before/after outcomes, denominator-aware sensitivity,
precision, specificity, call/no-call rates; descriptive sample-level Wilson
intervals only with dependence and design limitations. Alleles within a sample
and perturbations from a common sample are not independent replicates.

## Review, risks and rollback

Fresh installed Claude CLI sessions request claude-fable-5-1 explicitly at spec,
plan, milestones and final diff/report. Save requests, responses, exact source
hash/diff manifests (including untracked files), findings and dispositions.
Model unavailability is reported as unmet, never silently substituted.

Principal risks: phase claims without phase-block evidence, ambiguous repeat
projection, unobserved reference fill, missing molecular truth, simulator bias,
low-count estimates, and unsupported generic large-indel segmentation. Separate
validated fixes from experimental methods. Revert an individual change or retain
it as diagnostic-only if its scientific gate fails. Keep rejected experiments and
unmet workstreams in the final requirement/evidence matrix; archive only completed
plans. No output may claim complete success based solely on tests or reviewer
approval.

## Adversarial review decisions before final freeze

The specification and plan were reviewed in independent actual Claude Fable 5.1
sessions. Final nonempty JSON responses and modelUsage records are preserved.
The review files were empty WHILE those processes were running; they were not
completed empty reviews. Implementations began on fixed independent contracts
while reviews were running; subsequent findings are resolved before final freeze.

- Production segmentation retains candidate windows. The tested strict stop rule
  lost five named events among 106 exposed consensus outputs (three original HiFi,
  two later challenge); it fails the promotion gate. Keep it explicitly experimental
  and expose per-window fit/localization uncertainty and all trailing residues.
  Conservative event metrics compare lower bounds, never choose a favorable tie.
- Pipeline positive rule for legacy comparison is ANY reported mutation candidate;
  exact event rule additionally requires full parent/name/index on the matched
  allele. Supported rule is frameshift + exact template + exact VCF concordance.
  Legacy proximity support is reported separately and cannot be compared as if it
  were the new support rule. Annotation agreement is not arbitrary normalized
  biological event identity. VCF concordance is not independent evidence.
- Consensus writes `consensus_context` containing actual reference_path,
  full_consensus_path, chrom, sample, haplotype and 0-based trim_start/trim_end.
  Support replays selected VCF alleles, verifies the full consensus and trim, and
  reverts an event in whole-sequence context to check the classified parent.
  Homopolymer anchoring equivalence is tested with real bcftools normalization.
- Phase producer in this change is pre-existing genotype/PS evidence only: no
  read-backed phasing tool is being silently invented or represented as complete.
  All multiple heterozygous sites must be phased in one common PS. Single-site
  heterozygosity permits an unordered pair. Read-backed phase remains an unmet
  scientific objective unless a separately validated implementation is delivered.
- `run_status.json` schema1 states completed, insufficient_evidence or
  execution_failed with explicit error details. Insufficient coverage remains a
  nonzero CLI outcome for compatibility, now machine-distinguishable from tool
  failure. The evaluator preserves every expected denominator and ignores stale
  successes. It must not count two copies of one observed sequence as two recovered
  haplotypes without genotype evidence identifying two selections.
- Safety-contract changes (honest unresolved phase/support) are accepted only as
  such, with enumerated lost calls and exact reconstructions. They cannot pass an
  accuracy-improvement claim by increasing no-calls. Accuracy claims require the
  original no-regression/call-rate gates. Equivalence optimization has a separate
  exact-output gate. Unsupported ambitious algorithms remain unmet.
- The shipped ladder uses the distal first500 bases of its left10kb flank, while
  current trim anchors use proximal sequence. This is verified source behavior.
  Preserve the reference in this change to isolate other effects, record failed
  anchor/fixed trim in context, and do not claim flank anchoring or observed flank
  reconstruction. Correcting/regenerating this reference requires its own paired
  reference ablation; this known reference defect remains unresolved here.
- Read coverage must be reported separately from reference-based consensus. Variant
  absence is not covered reference confidence. If masking/coverage integration
  cannot be independently validated, report reference confidence unverified and
  leave coverage-aware reconstruction explicitly unmet, never imply exact sequence
  equality establishes empirical base support.
- Add experimental haploid-prior calling to F hypotheses; the review's claim that
  the diploid prior CAUSED low QUAL is unproven without matched Clair3 reruns.
  Do not promote this causal explanation or new caller mode from one observation.
- Fresh panel tests seed variation of specified simulator designs, not population
  generalization or dictionary biological validity. Add four genuinely new design
  strata per platform (ins25bp, non-X parent mutation, delGCCCA, equal-length H2
  mutation) when generation supports them; freeze actual case manifest before use.
  Read-source capture before merging is preferable to separate haplotype simulation
  because it preserves joint PCR allocation. If unavailable, report assignment
  truth missing rather than guess source from sequence alignment.

### Pre-freeze phase implementation extension

Actual installed WhatsHap1.7 passes linked cis/trans/shared SNP and indel fixtures.
Add optional same-length post-filter read-backed phasing, preserving every primary
record with unique temporary identities and an original-name map. Presence is
checked explicitly; missing optional tool returns unavailable, present tool errors
propagate. Existing usable/empty/missing/conflicting genotypes remain unchanged.
Only unchanged variant/genotype identities accepted by the common-PS policy may
replace the filtered VCF. Multi-site multiallelic phase remains unsupported when
the tool leaves it unresolved. This supersedes the earlier producer-only limitation
only for the validated bounded helper; upstream false length splits and unobserved
reference sequence are still unresolved. No threshold or QUAL change.

Two identical trimmed candidates cannot inherit independent VNTR evidence from
flank variants. Add per-side trim method and per-allele projection failure reason.
Specificity must retain ambiguous positive alarms in FP denominator; report
unresolved negative controls separately, never as TN. Candidate-output call rate
and resolved-reconstruction call rate are distinct measures.

Final sampling is fixed at all 32 main designs plus all three predeclared
perturbations per design (128 paired inputs). Four predetermined main cases,
normal60/61 and dupC120/140 in each platform, receive three alternating paired
repetitions for timing; all other inputs one pair. No selection depends on final
results. Perturbations are dependent sensitivity experiments, not 96 independent
biological replicates. Caller source/settings and evaluator hashes freeze before
reserved seeds are used.

### Final development decisions before reserved seeds

The optional read-backed helper is available only through the Python library's
explicit `read_phase=True`; it is disabled by default and has no new CLI flag.
On five cached difficult samples it recovered four rather than one of ten exact
sequences, but introduced an extra mutation call. Debugging traced this to a
length-assignment failure: all thirteen extracted records came from the short
allele, while the original input also contained five long-allele records. Phasing
split short-allele microvariation 12:1 and could not restore the missing allele.
Six WhatsHap parameter combinations did not correct this error. Minimum support
floors 2/3/5 are not justified: separate known-truth fixtures show downsampling
changes valid 30:2 support to 14:1. The bounded helper passes software fixtures;
the production reconstruction promotion requirement remains unmet.

An independent full 77-input development ablation of existing `--min-coverage 5`
improved HiFi exact count pairs 40/44 to 41/44 and perturbation pairs 9/24 to
15/24, but perturbation extra events increased 11 to 24 and three normal controls
became mutation positive. The default stays 10. No QUAL, ALT-selection, strict
segmentation, length replacement, support-floor, or mutation-rescue change is
promoted from these experiments. These decisions precede final data generation.

For final evaluation, the primary exact sequence/structure accuracy credits only
independently recoverable allele observations. Literal copies are retained in
separate metrics; one sequence duplicated without VNTR genotype evidence cannot
receive two reconstruction successes. Normal ambiguous positive alarms remain
false positives; unresolved negative observations remain outside true negatives.
Both candidate-output and resolved-reconstruction call rates are reported.
