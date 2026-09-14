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

- [ ] Add tests for swapped alleles, ties with different event scores, missing/extra alleles, invalid truth, failed stale output, normal no-call, wrong site/name/parent, duplicate extras.
  Concrete failure: truth dupC at H1 repeat25 with predictions dupCCCC or dupC repeat999 must have TP=0,FN=1,FP=1. Truth A/C with three predictions A/C/G cannot be exact diploid.
  ```python
  assert score_events([(25, 'X', 'dupC')], [(999, 'X', 'dupC')]) == (0, 1, 1)
  ```
  API may use typed event objects; retain the observable assertion.
- [ ] Run `uv run --locked --all-extras pytest tests/unit/test_evaluation*.py --no-cov` and save expected failures.
- [ ] Implement strict truth reconstruction using actual mutated FASTA, exact flanks and structure; preserve legacy provenance defects as warnings. Implement all optimum injective assignments and metric bounds.
- [ ] CLI writes every expected sample, rejects empty/missing input, distinguishes poor accuracy from execution/validation failure. No sample-name truth guessing.
- [ ] Run focused tests; evaluate exposed datasets with explicit inventories, compare endpoints and explain corrections.
- [ ] Reviewer checkpoint: adversarial matching/failure/denominator inspection. Coordinator updates benchmarking docs and requirement matrix.

## Task 2: Exact scoring optimization (workstream B)

**Owner/files:** scoring worker exclusively owns `src/muc_one_span/repeat_alignment.py`, new `tests/unit/test_repeat_alignment.py`, `scripts/benchmark_classification.py`.
**Interface:** preserve `edit_distance(str,str)->int` and `characterize_differences`; literal symbol equality, first-hit tie rules unchanged.
**Dependencies:** independent of Task 1/3; task 3 owns classify.py.

- [ ] Add independent scalar oracle and >=10,000 seeded pairs including empties, repeats, indels, lowercase, N/IUPAC and lengths around machine-word boundaries. Existing scalar already passes equivalence; assert new optimized entry-point or benchmark threshold separately rather than pretending equivalence tests fail.
  ```python
  assert edit_distance('ACGT', 'ACAGT') == 1
  assert characterize_differences('ACGT', 'ACAGT')[0]['pos'] == 3
  ```
- [ ] Run tests against baseline; record behavioral equivalence baseline and measured difficult-case bottleneck.
- [ ] Implement exact bit-vector recurrence with correct empty-string handling; leave traceback unchanged except correcting insertion comment.
- [ ] Benchmark alternating old/new scorer in the same classifier process and compare complete result dictionaries, then restore scorer safely. Benchmark output records repetitions, medians, equality, Python/platform and source/input hashes.
- [ ] Require zero random/dictionary differences and >=2x difficult-case median; batch exact controls to resolve <20% timing effects. Coordinator separately measures whole pipeline.
- [ ] Reviewer checkpoint: exact recurrence/oracle independence, tie/traceback semantics and controlled timing.

## Task 3: Classification and variant evidence (workstream C)

**Owner/files:** classification worker exclusively owns `classify.py`, `classify_types.py`, new `variant_support.py`, `tests/unit/test_classification_evidence.py` and necessary `tests/unit/test_classify.py` adaptations. Coordinate VCF parser changes with Task 4 owner.
**Interfaces:** legacy classify_sequence dictionary plus half-open consumed spans, unresolved regions, ambiguous base count, explicit heuristic score semantics. Exact support may consume richer parsed variants and optional actual reference sequence; unavailable projection must remain unavailable.

- [ ] Reproduce mixed-indel frameshift, tail omission and unreachable recovery; tests for X[:5]+X[7:30]+'A'+X[30:] net -1, balanced insertion/deletion net0 and X*5+'N'*29.
- [ ] Fix signed frameshift using existing `_compute_net_indel(diffs)%3`; include explicit sequence accounting. Test every consumed/unclassified interval partitions input with no overlap/gap.
- [ ] Establish justified recovery or unresolved policy for large insertion and terminal exact anchors. Retain full ideal-truth structures; do not claim approximate regions resolved.
- [ ] Replace nearby-record support with verified event identity/projection. Tests include unrelated SNP, following repeat indel, upstream insertion, shifted homopolymer equivalent and ambiguous projection. Legacy pos/qual-only fixtures now prove proximity only, not support.
- [ ] Run focused tests and ideal truth panel; validate exact event losses separately from intentional support-status corrections.
- [ ] Reviewer checkpoint: coordinates, misleading confidence, normalized equivalence, false supports and unresolved segmentation. Coordinator integrates reference context at CLI boundary if needed.

## Task 4: Explicit consensus and phase/failure states (workstream D)

**Owner/files:** consensus worker exclusively owns `calling.py`, `consensus.py`, `vcf.py`, new `phasing.py` if needed, corresponding calling/consensus/VCF tests and new integration phase tests. Coordinator owns cli.py/report changes.
**Interfaces:** `build_consensus(..., sample=None, haplotype=...)`; richer `parse_vcf_variants` dictionaries with chrom,pos,ref,alt,genotype,qual. Query failures propagate. Allele dictionaries carry phase/reconstruction state; return VCF path map stays path-only.

- [ ] Add failing query-error/malformed-data tests; remove silent RuntimeError-to-empty conversion only after failure observed.
- [ ] Add real bcftools fixtures: reference AAAAAAAAAA with 1|0 at2 A>C and 0|1 at8 A>G in one block => ACAAAAAAAA and AAAAAAAGAA; shared ALT in both; GT1|2 selects different ALTs; missing/unphased/disconnected phase never asserted resolved.
- [ ] Make sample/IUPAC behavior explicit. Resolve only proven phase; remove WT/all-ALT split assumptions and empty-call homozygosity assertion. Preserve uncertain candidate output with status wherever possible.
- [ ] Investigate read-backed phase using known read fixtures and collision-free record identities. Promote only if tests support all phase scenarios; otherwise mark read-backed phasing unmet and retain honest ambiguity.
- [ ] Run focused unit and real-tool tests. Coordinator updates persisted alleles after calling and explicit report status.
- [ ] Reviewer checkpoint: genotype indices, phase set connectivity, no-call/negative distinctions, single versus two observed alleles.

## Task 5: Molecule support and length/rescue investigation (workstreams E/F)

**Owner/files:** coordinator or fresh worker after Task 2 owns `alleles.py`, new `read_evidence.py`, `tests/unit/test_read_evidence.py`, `tests/unit/test_allele_evidence.py`, and ignored experiments. No overlap with calling/consensus owner until released.
**Interfaces:** separate candidate alignment count and primary-record support; robust span evidence returns observed span, anchor uncertainty, rejected reason and record identity. No default length substitution without gates.

- [ ] Reproduce flat plateau split, diffuse cluster exclusion and misleading read support. Build toy SAM records with colliding QNAME but distinct sequences.
  ```python
  # Flat zero-indel scores at 40,43,46 cannot demonstrate two valleys.
  assert split_flat_profile({40: 0.0, 43: 0.0, 46: 0.0}) is None
  ```
  Adapt test to existing `_split_cluster_by_indels` public-observable behavior.
- [ ] Add diagnostic primary-record support without dropping secondary fit records; retain count output unless evidence justifies changed inference. Explain legacy `reads` explicitly.
- [ ] Experiment with error-tolerant terminal anchors, gaps0/1/2/3, partial reads and asymmetric support; evaluate candidate set AND assignment. Explicitly document unavailable original molecule truth.
- [ ] Test local variant rescue only from actual allele/read evidence, contrasting nearby unrelated variants and normal controls; reject QUAL/global motif shortcuts.
- [ ] Profile stage totals, streaming/RSS/thread budget. Promote independent optimization only with full equality and controlled benefit; record rejected/deferred methods with evidence.
- [ ] Reviewer checkpoint: inference versus diagnostic evidence, failure/dropout rates and scientific promotion criteria.

## Task 6: Freeze, fresh paired validation, integration and final review

**Owner/files:** coordinator owns `cli.py`, relevant report fields, docs/development.md, docs/reference/limitations.md, docs/guides/benchmarking.md (worktree copy only), CHANGELOG.md, `.planning/` specification/plan/matrix/reviews/report, ignored validation driver/artifacts.
**Dependencies:** Tasks1–5 adjudicated; settings frozen.

- [ ] Record baseline/final source manifests including untracked files, environment/model/input hashes and exact commands. Before final seeds freeze a JSON settings/source manifest.
- [ ] Generate 24 specified seeds after freeze using MucOneUp and documented HiFi/ONT models, preserve all truth separately. Produce minority/partial perturbations with clear source limitations.
- [ ] Run old and new code on identical fresh inputs with alternating paired order; retain stages, wall times, thread configuration, peak memory where measured and failures. Never retune exposed final data.
- [ ] Evaluate using frozen matching/scoring, report per-sample regressions and matched-call-rate comparison, platform strata, counts and sample-level uncertainty/dependence.
- [ ] Run `make ci-check`, `make test-int`, explicit `pytest tests/integration -m e2e --no-cov` with real data/tools, `make docs-check`, `make build-check`; security only if applicable.
- [ ] Final Claude Fable5.1 session receives full tracked/untracked diff manifest, tests, benchmarks, scientific report and requirement matrix; reproduce findings, fix with tests and request fresh scoped review.
- [ ] Final diff inspection and preserve all source work. Archive only finished plan components. If scientific methods fail, report unmet requirements without a completion claim.
