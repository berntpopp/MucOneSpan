# MucOneSpan: 200-case simulation and improvement prompt

Copy the prompt below into Gemini/Antigravity. It defines **200 new sequencing
datasets: 100 diploid designs × HiFi and ONT**, with 140 development datasets and
60 protected final-validation datasets.

Research starting points include [NanoRepeat's template fitting and read
clustering](https://pubmed.ncbi.nlm.nih.gov/36262216/) and [TRGT's HiFi repeat
characterization](https://github.com/PacificBiosciences/trgt), without assuming
those approaches transfer directly to MUC1.

---

Act as an exceptionally rigorous senior bioinformatician, Python developer,
statistician, and performance engineer specializing in long-read MUC1 VNTR
reconstruction.

You are the coordinating agent. Use Superpowers, parallel experimental agents,
and actual external Claude Code and Codex CLI sessions. Your objective is to
design, generate, evaluate, investigate, and improve MucOneSpan using 200 new
MucOneUp sequencing datasets. Continue through implementation and validation;
do not stop after writing a plan.

Favor reproducible evidence and correct scientific contracts over plausible
stories, favorable aggregates, or increasingly elaborate heuristics.

## Authorization

I authorize repository inspection, isolated development, simulation, local
experiments, literature research, parallel agents, external model CLI reviews,
implementation, testing, documentation, and local commits.

Make routine decisions autonomously. This explicitly authorizes proceeding
through design and plan checkpoints without repeatedly requesting approval.
Document decisions and use independent review at those checkpoints.

Ask only when an unresolved decision materially changes the objective, requires
new spending beyond existing CLI access, or requires additional authorization.
Do not merge, push, tag, publish, or release.

Preserve unrelated work. Keep generated reads, results, indexes, models,
credentials, and bulky raw artifacts untracked.

## 1. Establish the actual starting point

Read AGENTS.md and scoped instructions, docs/development.md, contributing
guidance, CHANGELOG.md, and these guides:

- docs/guides/simulation-experiments.md
- docs/guides/configuration.md
- docs/guides/validation-results.md
- docs/guides/benchmarking.md
- docs/reference/limitations.md

Read the relevant .planning specifications, validation report,
requirement-evidence matrix, experiment investigations, and review dispositions.

Inspect Git state, tags, worktrees, relevant implementation, tests, and ignored
experimental artifacts. Do not assume main contains the released improvements.

The previous release was v0.11.0 at
9535f7ee02033cda3da9b22c5b73752e5473b016. Verify this rather than trusting it.
Use that immutable release as the comparator. If the working candidate contains
later changes, record and evaluate them separately.

Previously recorded evidence:

- Exact scoring improved difficult classification by about 7.955× locally;
  this was not a whole-pipeline speedup.
- Two known generated-data length failures remained: equal 60/60 and asymmetric
  25/140 alleles.
- Read-backed phasing improved some sequences but introduced a false event,
  so it remained disabled by default.
- Lower coverage/QUAL thresholds and simplified peak suppression had harmful
  tradeoffs.
- No fresh final scientific panel was executed for v0.11.0.

Treat all historical samples, challenge cases, perturbations, reports and
previously examined results as development evidence.

Discover the current MucOneUp release through its authoritative repository.
Compare it with the locally installed version and source revision. Pin the
chosen current release in an isolated environment; do not silently update a
shared environment. Record any inability to use the current release.

## 2. Use Superpowers throughout

Discover and read the installed skills rather than assuming their paths:

- using-superpowers
- brainstorming
- writing-plans
- using-git-worktrees
- systematic-debugging
- test-driven-development
- dispatching-parallel-agents
- subagent-driven-development or executing-plans
- receiving-code-review
- verification-before-completion

Store specifications, executable plans, decisions and evidence in .planning/.
Archive completed plans; preserve unmet requirements explicitly.

Use parallel agents with bounded hypotheses, exclusive file ownership and
separate result directories. Coordinate shared schemas before delegation.
The coordinator owns scientific consistency, integration, promotion decisions,
resource scheduling and final verification.

Do not run multiple unrestricted pipelines that oversubscribe CPU or memory.

## 3. Define the 200-dataset experiment before generation

Interpret 200 cases as:

- 100 distinct diploid biological designs.
- Each design simulated as PacBio HiFi and ONT.
- Exactly 200 new sequencing datasets, excluding historical regressions,
  pilot runs, repeated timings and derivative perturbations.

Generate each biological truth once and use the identical diploid sequences for
its two platform simulations. Use separate recorded sequencing-error seeds.
Verify identical biological truth rather than assuming seed reuse guarantees it.

Split by biological design before generation:

- 70 designs × 2 platforms = 140 development datasets.
- 30 designs × 2 platforms = 60 protected final-validation datasets.

Keep both platforms and every derivative of a design in the same split.
Prevent duplicate full diploid truths from crossing splits. Parameter strata
may overlap; exact truth, source reads and derivatives may not.

Use fresh, explicit seeds. Separate generation, error simulation, sampling,
splitting and algorithm randomness where supported. Never infer molecular
independence from different seeds or QNAMEs alone.

Construct a constrained, stratified design covering:

- At least 30 mutation-negative biological controls, distributed across splits.
- Every supported named mutation type, on either haplotype.
- Mutation positions near VNTR ends and repeat boundaries.
- Equal lengths with identical sequences and equal lengths with different
  repeat structures or variants.
- Repeat-count gaps of 1, 2 and 3, plus wider and strongly asymmetric gaps.
- Short, typical and long VNTRs, including reference-range edge cases.
- Diverse repeat-unit composition, interruption patterns and long X runs.
- Balanced and minority alleles, including substantial allelic imbalance.
- Low, moderate and high realized usable coverage.
- Full-length, partial and imperfect reads, and supported chimera/error cases.
- Appropriate HiFi and ONT error configurations.

Use pairwise/targeted interaction coverage instead of pretending 100 biological
designs exhaust the Cartesian product. Publish the design coverage matrix and
identify combinations that are not covered.

Inspect MucOneUp's actual capabilities before promising these strata. If a
feature requires a deterministic post-generation transformation, specify and
validate it, record its provenance, and label it as such. Do not present an
unsupported transformation as native simulator behavior.

Distinguish requested templates from retained reads, spanning reads, allele
coverage and independent source molecules.

Run a small separate pilot to validate generation, truth reconstruction,
input selection, model compatibility, disk usage and runtime. Pilot cases are
development-only and do not replace any of the 200 intended datasets.

Retain failed intended cases in the inventory. Retry infrastructure failures
with the same declared design where reproducible; never silently replace a
difficult case with an easier seed.

## 4. Protect final validation

The calling process must never receive generation truth, mutation labels,
expected allele lengths or sample names that encode the answer.

Development agents may inspect only development truth and results.
Use neutral caller-facing sample identifiers.

Keep final truth, reads and results outside development-agent contexts until
the candidate and evaluation protocol are frozen. A separate validation
executor may check generation integrity without revealing sample-level outcomes.

Run the immutable baseline and frozen candidate on the 60 final datasets only
after the freeze. Do not tune on those results.

If final results inform another algorithm change, label that panel exposed
development evidence. A new final panel requires new designs/seeds and a
documented extension of scope.

## 5. Require real, independent external model sessions

Use Gemini as coordinator, plus the installed CLIs:

- Claude Code: request claude-fable-5-1 explicitly.
- Codex CLI: request gpt-6-astra explicitly.

Verify:

- claude --version and supported help/flags.
- codex --version and codex exec --help.
- Actual access to each requested model with a minimal invocation.

Do not silently substitute a model or impersonate either reviewer with an
internal subagent. If unavailable, save the concrete error, continue useful
independent work, and mark that review requirement unmet.

Subject to locally verified flags, use fresh read-only sessions such as:

```bash
claude -p \
  --model claude-fable-5-1 \
  --effort high \
  --permission-mode plan \
  --tools "Read,Glob,Grep" \
  --output-format json \
  < review-request.md \
  > claude-response.json \
  2> claude-stderr.log

codex exec \
  --model gpt-6-astra \
  --sandbox read-only \
  --json \
  - \
  < review-request.md \
  > codex-events.jsonl \
  2> codex-stderr.log
```

Use unique paths per session. Verify actual model identity from available
execution metadata; do not equate command-line intent with proven access.

Require independent brainstorming from both models BEFORE showing either the
other's suggestions. Then synthesize disagreements into discriminating tests.

Request both models at:

- Experimental specification and acceptance-criteria freeze.
- Selected architectural designs before production implementation.
- Substantial implementation milestones.
- Final diff and scientific report.

Give reviewers exact revisions, complete diffs including new files, source
manifests, test evidence and relevant experiment artifacts. Never expose sealed
final-validation data during development review.

Require severity, location, failure scenario, evidence/reproducer, and suggested
validation. Reproduce substantive findings before fixing them. Record rejected
findings with evidence and obtain fresh review of substantive fixes.

## 6. Build trustworthy measurement first

Reuse and improve the maintained generation, benchmark and evaluation tooling.
Do not accumulate disconnected one-off scripts as the only reproducible method.

Maintain a manifest and append-only run ledger containing:

- Source revisions and dirty-diff hashes.
- Simulator, caller and external-tool versions.
- Exact commands, effective settings and environment locks.
- Platform/model identities and file hashes.
- Input/truth/reference hashes and seeds.
- Expected, generated, attempted, completed, failed and no-call states.
- Realized usable/spanning reads and allele support where knowable.
- Exit codes, timeouts, stderr, stage timings and resource configuration.

Prevent stale successful output from overriding failed execution.
Use content-aware caching: input, code, reference, model and relevant parameters.
Never reuse an incompatible reference index.

Keep allele assignment truth separate from calling. If simulator source labels
are unavailable, report assignment accuracy as unavailable. Do not manufacture
source labels from ambiguous QNAMEs or truth-based alignment.

Measure:

- Exact full VNTR sequences and complete ordered repeat structures.
- Exact count pairs and signed count errors.
- Exact mutation name, parent, total-repeat position and matched haplotype.
- Missing/extra alleles and false extra calls on mutation-positive samples.
- Assignment ambiguity and phase/connectivity failures.
- Sensitivity, precision, normal-control false-alarm rate, specificity,
  call rate and no-call rate, with explicit denominators.
- Literal sequence equality versus independently supported allele recovery.

Use one-to-one allele matching. Haplotype numbering is permutation-invariant
unless identifiable external phase labels exist. Define deterministic primary
assignment and tie handling before comparing outcomes; do not select a different
favorable assignment for each metric.

Distinguish confident negative, insufficient evidence, ambiguous reconstruction,
execution failure and completed inaccurate reconstruction.

Report platform and difficult-case strata separately. Paired platform datasets
share biological truth and are not independent biological samples.

## 7. Diagnose failures before optimizing

Run the baseline across all 140 development datasets and historical regressions.
Create a sample-level failure atlas and trace errors through:

generation → usable reads → initial mapping → length inference →
read assignment → allele reference → variant calling → phasing →
consensus → trimming → classification → evaluation.

Find the FIRST incorrect stage, then quantify downstream consequences.
Use controlled stage substitutions to test causality.

Truth-derived references, true assignments or true sequences may be used only
in isolated, clearly labelled development oracle experiments to locate the
performance ceiling. Their outputs must never count as deployable improvements.

Profile actual wall time, CPU time, peak memory and stage costs before choosing
performance work. Preserve record identity through BAM/FASTQ transformations.

## 8. Brainstorm and run competing experiments

Build a hypothesis registry. For every proposal record:
mechanism, expected benefit, possible harm, affected strata, configuration,
smallest decisive experiment, cost, falsification criterion, promotion gate
and rollback.

Investigate at least these competing families:

### A. Reference architecture

- Current ladder versus an X-only discovery reference with explicit flanks and
  fixed-repeat accounting.
- Two-pass calling: coarse discovery followed by references containing
  read-supported resolved non-X units.
- A bounded panel of plausible repeat-composition references.
- Correct VNTR-proximal flank construction versus current flank selection.
- A repeat-unit graph or alternative representation if linear-reference bias
  remains the dominant demonstrated failure.

### B. Iterative allele-specific references

- Alternate read assignment, consensus and reference reconstruction.
- Start from multiple hypotheses, not only the initial winning reference.
- Freeze confident blocks while keeping unresolved blocks explicitly unresolved.
- Reconsider previously discarded reads rather than locking in early mistakes.
- Compare one-pass, two-pass and bounded iterative variants.
- Detect cycles, collapse of minority alleles and reinforcement of early errors.
- Use cross-fitting or independent read partitions to assess proposed updates
  where coverage permits.
- Compare candidate fit in a common coordinate/scoring framework; an improved
  alignment score against a changed reference is not automatically evidence
  of a better biological reconstruction.
- Specify convergence, maximum iterations and fallback behavior in configuration.

### C. Length inference and molecule assignment

- Error-tolerant two-flank spanning evidence with partial-read handling.
- Candidate-reference fit separated from independent read support.
- One-versus-two-component models with discrete repeat lengths and an explicit
  unresolved outcome.
- Soft assignment, composition-aware clustering and graph connectivity.
- Selection-aware treatment of minority alleles and per-contig filtering.
- Distinguish one sequence represented twice from two independently supported
  equal-length haplotypes.
- Preserve distinct simulator reads despite QNAME collisions.

### D. Consensus and phasing

- Read-backed SNP/indel linkage, phase-set connectivity and shared variants.
- Local or partial-order consensus on defensibly assigned reads.
- Evidence-aware masking of uncovered reference sequence.
- Haplotype-specific polishing and localized mutation rescue.
- Indel normalization and exact coordinate projection under changing references.
- An explicit ambiguous result when the evidence cannot determine phase.

### E. Technology-specific parameter matrices

Inspect supported versions before choosing parameter ranges. Consider:

- minimap2 presets, seed/chaining/scoring and secondary-alignment settings.
- Clair3 platform/model, supported quality and coverage controls.
- bcftools sample, genotype, haplotype and filtering behavior.
- WhatsHap supported mapping/downsampling settings.
- MucOneSpan thresholds, reference layouts and consensus/classification controls.

Separate technology effects from model-version and coverage effects.
Do not make global QUAL reduction or applying every ALT the default rescue.

Use staged screening, factorial contrasts or space-filling designs followed by
focused interaction experiments. Do not launch an unbounded Cartesian grid.
Use grouped development cross-validation or a development selection split.

### F. Performance

Investigate streaming, eliminating repeated scans, reference/index reuse,
batched exact scoring, subprocess startup, and thread budgets only where
profiles justify them. Include caching disabled and warm-cache conditions.

Propose additional ideas from primary literature and tool documentation.
Use GeneFoundry MCP/PubTator if available; otherwise use primary sources.
Verify publication identity and applicability. A method effective for short STRs
or whole-genome HiFi data may not work for GC-rich MUC1 amplicon VNTRs.

## 9. Promote changes using predeclared criteria

Before candidate selection or final evaluation, write concrete numerical
acceptance criteria for the primary scientific endpoint, safety constraints,
runtime and memory. Use historical evidence and pilot cost estimates to choose
defensible values; do not choose them after seeing favorable final results.

Prioritize complete diploid reconstruction while protecting exact mutation
identity, normal specificity, minority alleles and call rate.

Report absolute changes, paired wins/losses/ties and per-sample regressions.
Compare at matched call rates where meaningful. A normal sample changed into
a no-call is not a specificity improvement.

Use design-level paired uncertainty estimates or resampling that preserves
platform/derivative dependence. Report small stratum denominators and wide
uncertainty honestly. Do not treat reads as independent patient observations
or claim clinical validation from simulations.

For pure performance changes require scientific output equivalence, with an
explicit list of excluded non-scientific fields such as timings and paths.
Benchmark alternating paired repetitions at equal thread and resource budgets.
Report local-stage and whole-pipeline speedups separately.

Promote only a small number of candidates. Perform ablations to determine which
components actually help. If a promising method fails its gate, document and
reject it rather than quietly relaxing the criterion.

## 10. Implement, verify and finish

All scientific and experimental parameters must be configuration-driven and
validated. No branches keyed to sample names, seeds, expected mutations, truth
paths or known failing examples.

Use existing subprocess abstractions, locked environments, typed interfaces
and deterministic tests. Maximum 649 physical lines per authored code,
configuration or template file.

Preserve CLI compatibility where practical. Document additive fields and any
intentional scientific contract correction.

Run:

- make ci-check
- make test-int
- Explicit generated-data end-to-end selection
- make docs-check
- make build-check for packaging/resource/dependency changes
- make security-check for dependency/security changes

Do not weaken coverage, assertions, tolerances or expected failures.
Replace an expected failure only after reproducing corrected behavior and
testing neighboring cases.

After freezing code, settings, models and evaluation:

- Run baseline and candidate on all 60 final-validation datasets.
- Complete the 200-dataset inventory and preserve every failed/no-call outcome.
- Run final scientific analysis and paired performance measurements.
- Obtain independent Claude Fable 5.1 and Codex GPT-6 Astra reviews.
- Commit maintained implementation, tests, documentation and compact planning
  evidence locally; leave raw generated artifacts untracked.

Deliver:

1. Dataset design, manifests and reproducible generation commands.
2. Failure atlas and stage-level causal investigations.
3. Hypothesis registry, parameter matrices, rejected methods and ablations.
4. Tested production improvements with a requirement-to-evidence matrix.
5. Per-sample baseline/candidate tables and uncertainty-aware summary.
6. Exact commands for reproducing and extending the next experiment.
7. Review requests, responses and finding dispositions.
8. A concise TLDR: what improved, what regressed, what remains unresolved,
   actual checks, runtime/resource cost, and the next most informative experiment.

Keep progress updates concise and evidence-based. Do not claim completion
because the dataset was generated, a plan was written, or unit tests passed.

Begin with repository/environment inspection, then produce the experimental
specification and independent brainstorming requests. Continue autonomously
through the full workflow.
