# Closed implementation plan; final scientific phase superseded

Closed by the user's explicit direction to finish with the validated work, release
0.11.0 and provide future experiment design. Tasks1–5 are implemented or explicitly
adjudicated; the originally planned fresh128-input phase below is NOTcompleted.
Its unchecked items are retained as unmet/superseded, not converted to passes.
See ../release-scope-0.11.0.md, ../requirement-evidence-matrix.md and
../2026-09-14-production-validation-report.md for the final scope and evidence.

# Production reconstruction implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox syntax for tracking. User authorizes continuous implementation and parallel agents with exclusive ownership.

**Goal:** Deliver independently validated production fixes and explicitly assess all scientific workstreams.
**Architecture:** Offline evaluation is isolated from inference. Exact scoring, classification evidence, consensus phase and read support have separate responsibilities and tests. Algorithms that fail scientific gates remain experiments.
**Tech Stack:** Python >=3.10, uv locked environment, pytest, existing samtools/minimap2/bcftools/Clair3 abstraction, MucOneUp simulations.
**Spec:** `.planning/2026-09-14-production-spec.md`

## Global constraints

- Baseline d8390b3; preserve original checkout and ignored data.
- No push, merge, release, new dependency or global QUAL/ALT change.
- <=649 physical lines per authored code/config/template file.
- Existing flags/imports remain; semantic corrections get additive status metadata and migration documentation.
- Observe failing tests before behavior edits; meaningful scalar oracle for exact scoring.
- Caller never reads generation truth. All 53 exposed full samples and 24 perturbations are development.
- All final seeds generated only after algorithm settings and source hashes freeze.
- Coordinator owns integration, scientific consistency, docs, checks and review dispositions.

## Task 1: Strict offline evaluation (workstream A)

**Owner/files:** evaluation worker exclusively owns new `src/muc_one_span/evaluation/` modules (`__init__.py`, `models.py`, `truth.py`, `artifacts.py`, `matching.py`, `scoring.py`), `scripts/evaluate.py`, and `tests/unit/test_evaluation*.py`. Follow evaluation-design-notes.md.
**Interfaces:** `load_truth(sample_dir, repeat_dict)`, `load_observation(result_dir, run_record)`, `match_alleles(truth,predictions)`, `evaluate_sample`, `aggregate`; thin CLI accepts result root, truth root, explicit expected samples and output. Offline only.
**Dependencies:** repeat dictionary; exact distance via repeat_alignment public API after Task 2 (temporary local scalar use acceptable before then).

- [x] Add tests for swapped alleles, ties with different event scores, missing/extra alleles, invalid truth, failed stale output, normal no-call, wrong site/name/parent, duplicate extras.
  Concrete failure: truth dupC at H1 repeat25 with predictions dupCCCC or dupC repeat999 must have TP=0,FN=1,FP=1. Truth A/C with three predictions A/C/G cannot be exact diploid.
  ```python
  assert score_events([(25, 'X', 'dupC')], [(999, 'X', 'dupC')]) == (0, 1, 1)
  ```
  API may use typed event objects; retain the observable assertion.
- [x] Run `uv run --locked --all-extras pytest tests/unit/test_evaluation*.py --no-cov` and save expected failures.
- [x] Implement strict truth reconstruction using actual mutated FASTA, exact flanks and structure; preserve legacy provenance defects as warnings. Implement all optimum injective assignments and metric bounds.
- [x] CLI writes every expected sample, rejects empty/missing input, distinguishes poor accuracy from execution/validation failure. No sample-name truth guessing.
- [x] Run focused tests; evaluate exposed datasets with explicit inventories, compare endpoints and explain corrections.
- [x] Reviewer checkpoint: adversarial matching/failure/denominator inspection. Coordinator updates benchmarking docs and requirement matrix.

## Task 2: Exact scoring optimization (workstream B)

**Owner/files:** scoring worker exclusively owns `src/muc_one_span/repeat_alignment.py`, new `tests/unit/test_repeat_alignment.py`, `scripts/benchmark_classification.py`.
**Interface:** preserve `edit_distance(str,str)->int` and `characterize_differences`; literal symbol equality, first-hit tie rules unchanged.
**Dependencies:** independent of Task 1/3; task 3 owns classify.py.

- [x] Add independent scalar oracle and >=10,000 seeded pairs including empties, repeats, indels, lowercase, N/IUPAC and lengths around machine-word boundaries. Existing scalar already passes equivalence; assert new optimized entry-point or benchmark threshold separately rather than pretending equivalence tests fail.
  ```python
  assert edit_distance('ACGT', 'ACAGT') == 1
  assert characterize_differences('ACGT', 'ACAGT')[0]['pos'] == 3
  ```
- [x] Run tests against baseline; record behavioral equivalence baseline and measured difficult-case bottleneck.
- [x] Implement exact bit-vector recurrence with correct empty-string handling; leave traceback unchanged except correcting insertion comment.
- [x] Benchmark alternating old/new scorer in the same classifier process and compare complete result dictionaries, then restore scorer safely. Benchmark output records repetitions, medians, equality, Python/platform and source/input hashes.
- [x] Require zero random/dictionary differences and >=2x difficult-case median; batch exact controls to resolve <20% timing effects. Coordinator separately measures whole pipeline.
- [x] Reviewer checkpoint: exact recurrence/oracle independence, tie/traceback semantics and controlled timing.

## Task 3: Classification and variant evidence (workstream C)

**Owner/files:** classification worker exclusively owns `classify.py`, `classify_types.py`, new `variant_support.py`, `tests/unit/test_classification_evidence.py` and necessary `tests/unit/test_classify.py` adaptations. Coordinate VCF parser changes with Task 4 owner.
**Interfaces:** legacy classify_sequence dictionary plus half-open consumed spans, unresolved regions, ambiguous base count, explicit heuristic score semantics. Exact support may consume richer parsed variants and optional actual reference sequence; unavailable projection must remain unavailable.

- [x] Reproduce mixed-indel frameshift, tail omission and unreachable recovery; tests for X[:5]+X[7:30]+'A'+X[30:] net -1, balanced insertion/deletion net0 and X*5+'N'*29.
- [x] Fix signed frameshift using existing `_compute_net_indel(diffs)%3`; include explicit sequence accounting. Test every consumed/unclassified interval partitions input with no overlap/gap.
- [x] Establish justified recovery or unresolved policy for large insertion and terminal exact anchors. Retain full ideal-truth structures; do not claim approximate regions resolved.
- [x] Replace nearby-record support with verified event identity/projection. Tests include unrelated SNP, following repeat indel, upstream insertion, shifted homopolymer equivalent and ambiguous projection. Legacy pos/qual-only fixtures now prove proximity only, not support.
- [x] Run focused tests and ideal truth panel; validate exact event losses separately from intentional support-status corrections.
- [x] Reviewer checkpoint: coordinates, misleading confidence, normalized equivalence, false supports and unresolved segmentation. Coordinator integrates reference context at CLI boundary if needed.

## Task 4: Explicit consensus and phase/failure states (workstream D)

**Owner/files:** consensus worker exclusively owns `calling.py`, `consensus.py`, `vcf.py`, new `phasing.py` if needed, corresponding calling/consensus/VCF tests and new integration phase tests. Coordinator owns cli.py/report changes.
**Interfaces:** `build_consensus(..., sample=None, haplotype=...)`; richer `parse_vcf_variants` dictionaries with chrom,pos,ref,alt,genotype,qual. Query failures propagate. Allele dictionaries carry phase/reconstruction state; return VCF path map stays path-only.

- [x] Add failing query-error/malformed-data tests; remove silent RuntimeError-to-empty conversion only after failure observed.
- [x] Add real bcftools fixtures: reference AAAAAAAAAA with 1|0 at2 A>C and 0|1 at8 A>G in one block => ACAAAAAAAA and AAAAAAAGAA; shared ALT in both; GT1|2 selects different ALTs; missing/unphased/disconnected phase never asserted resolved.
- [x] Make sample/IUPAC behavior explicit. Resolve only proven phase; remove WT/all-ALT split assumptions and empty-call homozygosity assertion. Preserve uncertain candidate output with status wherever possible.
- [x] Investigate read-backed phase using known read fixtures and collision-free record identities. Promote only if tests support all phase scenarios; otherwise mark read-backed phasing unmet and retain honest ambiguity.
- [x] Run focused unit and real-tool tests. Coordinator updates persisted alleles after calling and explicit report status.
- [x] Reviewer checkpoint: genotype indices, phase set connectivity, no-call/negative distinctions, single versus two observed alleles.

## Task 5: Molecule support and length/rescue investigation (workstreams E/F)

**Owner/files:** coordinator or fresh worker after Task 2 owns `alleles.py`, new `read_evidence.py`, `tests/unit/test_read_evidence.py`, `tests/unit/test_allele_evidence.py`, and ignored experiments. No overlap with calling/consensus owner until released.
**Interfaces:** separate candidate alignment count and primary-record support; robust span evidence returns observed span, anchor uncertainty, rejected reason and record identity. No default length substitution without gates.

- [x] Reproduce flat plateau split, diffuse cluster exclusion and misleading read support. Build toy SAM records with colliding QNAME but distinct sequences.
  ```python
  # Flat zero-indel scores at 40,43,46 cannot demonstrate two valleys.
  assert split_flat_profile({40: 0.0, 43: 0.0, 46: 0.0}) is None
  ```
  Adapt test to existing `_split_cluster_by_indels` public-observable behavior.
- [x] Add diagnostic primary-record support without dropping secondary fit records; retain count output unless evidence justifies changed inference. Explain legacy `reads` explicitly.
- [x] Experiment with error-tolerant terminal anchors, gaps0/1/2/3, partial reads and asymmetric support; evaluate candidate set AND assignment. Explicitly document unavailable original molecule truth.
- [x] Test local variant rescue only from actual allele/read evidence, contrasting nearby unrelated variants and normal controls; reject QUAL/global motif shortcuts.
- [x] Profile stage totals, streaming/RSS/thread budget. Promote independent optimization only with full equality and controlled benefit; record rejected/deferred methods with evidence.
- [x] Reviewer checkpoint: inference versus diagnostic evidence, failure/dropout rates and scientific promotion criteria.

## Task 6: Freeze, fresh paired validation, integration and final review

**Owner/files:** coordinator owns `cli.py`, relevant report fields, docs/development.md, docs/reference/limitations.md, docs/guides/benchmarking.md (worktree copy only), CHANGELOG.md, `.planning/` specification/plan/matrix/reviews/report, ignored validation driver/artifacts.
**Dependencies:** Tasks1–5 adjudicated; settings frozen.

- [ ] Record baseline/final source manifests including untracked files, environment/model/input hashes and exact commands. Before final seeds freeze a JSON settings/source manifest.
- [ ] Generate 32 specified seeds and all96 dependent perturbations after freeze using MucOneUp and documented HiFi/ONT models, preserve all truth separately. Produce minority/partial perturbations with clear source limitations.
- [ ] Run old and new code on identical fresh inputs with alternating paired order; retain stages, wall times, thread configuration, peak memory where measured and failures. Never retune exposed final data.
- [ ] Evaluate using frozen matching/scoring, report per-sample regressions and matched-call-rate comparison, platform strata, counts and sample-level uncertainty/dependence.
- [x] Run `make ci-check`, `make test-int`, explicit `pytest tests/integration -m e2e --no-cov` with real data/tools, `make docs-check`, `make build-check`; security only if applicable.
- [x] Final Claude Fable5.1 session receives full tracked/untracked diff manifest, tests, benchmarks, scientific report and requirement matrix; reproduce findings, fix with tests and request fresh scoped review.
- [x] Final diff inspection and preserve all source work. Archive only finished plan components. If scientific methods fail, report unmet requirements without a completion claim.

## Historical milestone disposition update (before the scoped release)

Tasks1–5 implementation/investigation checklists are complete as bounded tasks.
This does not promote failed scientific methods: default read-backed full
haplotypes, joint all-read assignment, reference coverage and local rescue remain
unmet, as listed in requirement-evidence-matrix.md. Actual Claude milestone2
review fixes are implemented with focused regressions and a completed fresh review;
its N1/N2 follow-ups are fixed. User steering added the runtime-configuration
extension before final freeze, specified in 2026-09-14-runtime-configuration-spec.md.
Current gates:626 units,93.59%coverage;47 integration tests;6 generated-data tests
with2 unchanged strict expected failures;docs and build/resource checks pass.
Actual Claude configuration review is running against an immutable source snapshot.
Task6 remains open until fresh paired validation, final review, dispositions and
complete reporting are finished. No reserved final seed has yet been generated.


## Explicit follow-up: repair the two generated-data expected failures

User requested fixing failures after the current gates passed with two strict
expected failures, and requested completing and documenting experiments. Before
any final seed exposure, investigate equal60/60 false valley splitting and
asymmetric25/140 loss. Exclusive read-only prototype ownership: equal-case agent,
asymmetric-case agent, coordinator scientific adjudication; production changes
require a concrete tested design and reviewed evidence first.

Candidate experiments must use all77 exposed regression inputs, preserve adjacent
0/1/2/3-repeat contrasts and minority/partial-read denominators, and obey the
previously fixed no-regression/promotion criteria. Neither passing only the two
failed examples nor changing xfail/tolerance annotations establishes success.
A failed method is recorded with parameters, per-sample regressions and why it
was not promoted. Keep independently validated production fixes. No blanket
QUAL/coverage relaxation, forced two components, or sample/truth branches.

Additional independent audit checks whether existing BAM-to-FASTQ conversion
retains distinct primary records when simulator QNAMEs collide. This is a data
retention contract, not a rationale for deduplicating molecules or inventing
read-source labels. Research evidence is in length-model-research.md and prior
GeneFoundry/PubTator notes. Final settings/source freeze follows adjudication.

## Scoped release closure

Final software gates: 694 unit tests on Python 3.10–3.14, 93.48% primary coverage,
47 integration passes, 6 generated-data passes and 2 unchanged strict xfails;
docs, build, security and quality pass. Actual final Claude Fable 5.1 review found
no blocking defects; eight non-blocking dispositions are committed. The earlier
configuration checkpoint above is historical. Fresh paired validation remains
unperformed under the explicitly changed release scope, not a passed endpoint.
