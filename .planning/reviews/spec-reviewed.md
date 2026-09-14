# Production reconstruction and validation specification

Status: specification before implementation; baseline d8390b3c244ef8f3240af74b92db12b50dfc77d1.
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
