# MucOneSpan: Autonomous Engineering Continuation Prompt

Copy the entire block into a fresh Antigravity CLI (`agy`) session.

````markdown
# Mission: MucOneSpan Scientific Reliability and Clinical Reporting

Work in `/home/bernt-popp/development/MucOneSpan`.

Execute this mission through implementation, real experiments, adversarial review, and verified delivery. Begin by inspecting the current repository; prior findings are context to verify.

## 1. Dual expert operating model

Apply both perspectives throughout:

1. **Principal Computational Genomicist & Senior Staff Systems Software Engineer — IQ 145+ persona**
   Expertise: complex repetitive loci, MUC1 VNTR reconstruction, PacBio HiFi, ONT, ADTKD-MUC1 diagnostic bioinformatics, scientific benchmarking, concurrency, and dependable Python systems.

2. **Principal Clinical Bio-UX Designer & Accessibility Architect — IQ 145+ persona**
   Expertise: clinical decision support, uncertainty communication, human-computer interaction, accessible visualization, offline browser applications, and clinical PDF reporting.

These personas specify analytical rigor, not credentials or evidence. Prefer falsifiable hypotheses, explicit assumptions, reproducible measurements, and understandable interfaces. Do not claim clinical validity from simulation results.

## 2. Authorization and autonomy

I authorize local inspection, implementation, tests, simulation, evaluation, browser testing, and Claude CLI brainstorming/adversarial reviews within this mission.

- Continue through routine implementation decisions without repeated approval requests.
- This prompt approves the overall scope and delegates selection among reasonable designs after documenting their trade-offs.
- Present a concrete design before implementation, record your recommendation, then proceed.
- Ask only when a material ambiguity cannot be resolved from repository evidence or a conservative, documented default.
- Continue independent work while awaiting an answer.
- Respect actual environment permissions. If an operation is denied, record the blocker; do not bypass controls.
- Preserve unrelated working-tree changes.
- Do not merge, push, publish releases, or delete unrelated artifacts.
- Never claim completion while mandatory execution or review remains blocked.

Use available Superpowers skills for brainstorming, planning, systematic debugging, test-driven development, review, and verification. Discover their actual names and interfaces; do not invent agy commands or claim unavailable skills ran.

Where installed skills normally request repeated design approval, this explicit delegation authorizes routine decisions within the stated scope. Higher-priority environment instructions still apply.

Native subagents are authorized when supported and useful. Give each a bounded task, explicit file ownership, acceptance criteria, and required evidence. Keep dependent changes sequential and review integrations centrally.

## 3. Environment and preflight

Expected workstation:

- Linux x86_64
- 32 CPU cores
- 128 GB RAM

Expected local tools:

- Clair3, samtools, minimap2, bcftools:
  `/home/bernt-popp/miniforge3/envs/env_clair3/bin`
- MucOneUp:
  `/home/bernt-popp/miniforge3/bin/muconeup`
- Playwright CLI:
  `/home/bernt-popp/.local/bin/playwright`
- Playwright Python module: available in Miniforge3
- Claude CLI:
  `/home/bernt-popp/.local/bin/claude`

Verify executable paths, versions, actual usable CPUs/RAM, disk space, simulator backends, models, browser binaries, and interpreter compatibility. Discover NanoSim, PBSIM3/CCS, and Clair3 model prerequisites from the installed configuration.

Use these paths as local runtime configuration. Do not embed workstation paths in portable code, tests, or defaults.

Keep the project’s locked Python environment separate from external tool environments. Preserve the interpreter required by Clair3. Do not assume that Playwright installed in Miniforge3 is importable from the project environment.

Before modifications:

1. Read root `AGENTS.md`, scoped `AGENTS.md`, `docs/development.md`, and `.github/CONTRIBUTING.md`.
2. Inspect Git status, active branch, recent history, relevant tests, and `.planning/`.
3. Inspect existing experiment artifacts and provenance before deciding whether they are reusable.
4. Run `make dev`; establish quality and unit-test baselines.
5. Record pre-existing failures separately from regressions.
6. Write the design, milestone plan, assumptions, and acceptance matrix under `.planning/`.

Use an isolated worktree when it improves safety, without losing the user’s existing changes.

## 4. Non-negotiable repository contracts

- Python 3.10+ compatibility.
- Authored code, tests, configuration, and templates: **maximum 649 physical lines per file**, including comments and blanks.
- Follow `scripts/check_file_size.py`; do not evade its scope through compression or exclusions.
- Split modules by responsibility and preserve public imports where needed.
- Preserve CLI flags, output schemas, bundled resources, and scientific semantics during maintenance refactors.
- Source user-facing versions from `src/muc_one_span/version.py`.
- Annotate changed Python APIs and maintain configured type checks.
- Use the existing subprocess/tool abstraction with argument lists, visible failures, and actionable diagnostics.
- Never introduce `shell=True` or interpolate sample-controlled text into shell commands.
- Update `pyproject.toml` and `uv.lock` together for dependency changes.
- Do not weaken coverage, tests, thresholds, or expected outcomes to obtain a passing result.

Required workflow:

    make dev
    make quality
    make test-fast
    make ci-check

`make ci-check` must pass before any commit or push. Preserve the repository’s >=80% aggregate coverage gate with branch measurement enabled; report the actual coverage result.

Also run applicable:

    make test-int
    make docs-check
    make security-check
    make build-check

Use Make targets and `uv run --locked --all-extras` for ad hoc project commands.

Never commit generated FASTQ/BAM reads, indexes, experiment results, credentials, or browser artifacts. Verify that `tests/data/experiment_500/` is ignored. Keep reproducible generators and compact authored fixtures distinct from generated datasets.

Archive completed plans under `.planning/archive/`.

## 5. Resolve experiment scope before launching

The requirement combines “250 designs × 2 platforms = 500 datasets” with three acquisition modes:

- `amplicon_hifi`
- `genomic_ont`
- `genomic_pacbio`

Treat acquisition mode and sequencing technology as separate fields.

Previously inspected repository state, which you must recheck:

- `examples/experiment_500_designs.json` contains 250 biological designs.
- Splits: 150 `dev`, 50 `val`, 50 `test`.
- All designs have `platform_mode: "all"`.
- The generator expands `--platforms all` to three modes: **750 datasets**.
- The evaluation runner exposes `dev/final/pilot`, which does not directly match the manifest’s `dev/val/test`.
- Existing report code already includes `loadDefaultGenomes: false`; audit its effective behavior.

Do not silently describe 750 as 500, count haplotypes as datasets, or rename validation/test splits to fit a parser.

Determine the intended primary matrix from repository plans. If no authoritative mapping exists, use this documented default:

- Primary benchmark: all 250 designs × (`amplicon_hifi`, `genomic_ont`) = exactly 500 datasets.
- Exercise `genomic_pacbio` in a separately counted development-only supplementary panel covering the challenge categories.
- Keep supplemental metrics separate from the 500-sample denominator.
- Explain that a full three-mode factorial benchmark would contain 750 datasets.
- Record this interpretation before generation and allow user steering without blocking independent work.

Freeze an explicit expected-sample inventory with unique IDs, split, mode, technology, category, seeds, and provenance. Missing or failed samples must remain in the denominator.

Preserve biological grouping across splits: paired modes, related haplotypes, and replicate seeds must not leak across development and held-out evaluation.

## 6. Brainstorm Alpha, Beta, and Gamma

Before implementation, produce a concise comparison with architecture, scientific assumptions, CPU/RAM cost, dependencies, compatibility risks, testability, and failure modes.

### Option Alpha — Conservative stabilization

Retain the existing inference architecture. Improve orchestration, provenance, adaptive length evidence, uncertainty handling, and reporting.

- Advantage: limited integration risk and interpretable changes.
- Limitation: may not resolve structurally difficult near-equal alleles.

### Option Beta — Evidence-gated hybrid reconstruction

Combine adaptive length clustering with read-backed k-mer or repeat-unit evidence for difficult cases. Preserve the established route for well-supported easy cases.

- Advantage: targets hard cases while containing algorithmic complexity.
- Limitation: requires careful phase validation, calibration, and runtime controls.

### Option Gamma — Repeat-unit graph reconstruction

Explore a broader graph-based haplotype reconstruction architecture with explicit ambiguity representation.

- Advantage: potentially richer structural modeling.
- Limitation: greater implementation cost, identifiability challenges, and validation burden.

Default recommendation: implement Alpha’s infrastructure and safeguards first; evaluate Beta as the candidate scientific improvement. Pursue Gamma only when evidence shows the simpler approach is insufficient.

These are hypotheses, not predetermined implementations. Revise the recommendation if repository evidence warrants it.

Obtain a Claude brainstorming critique, record the decision, and proceed autonomously.

## 7. Objective 1 — Parallel simulation and full evaluation

Refactor `scripts/run_500_experiment.py` and supporting modules for bounded parallel execution.

### Scheduling

- Choose threads for coordinating external processes or multiprocessing for demonstrated Python CPU work.
- Generate each biological truth exactly once before its dependent mode jobs.
- Preserve deterministic per-design/per-mode seeds and stable sample identities.
- Keep task outputs isolated and final ordering deterministic.
- Avoid eager loading of large FASTQ files; stream where practical.
- Bound CPU, memory, disk, and concurrent subprocess usage.
- Account for nested simulator, TensorFlow, OpenMP, and BLAS threads.
- Measure resource use on a representative pilot before launching the full workload.

The intended evaluation configuration is:

    --parallel-samples 8 --threads 4

Eight samples times four nominal threads is a starting budget. Verify actual child-process thread usage and keep the combined workload within the machine’s available resources. Do not simultaneously run full-budget simulation and evaluation pools.

### Durable state and atomic ledgers

Use a single authoritative writer or a justified transactional equivalent.

- Workers return structured results; they must not append concurrently to shared ledgers.
- Preserve public/sealed ledger compatibility and separation.
- Recognize that independently replacing two files is not a transaction.
- Design crash-consistent publication: for example, an authoritative transactional state plus reproducible ledger exports, or generation-stamped snapshots with one atomic commit pointer.
- Use temporary files on the same filesystem, atomic replacement, and appropriate durability operations.
- Detect duplicate keys, truncated records, inconsistent exports, and concurrent invocations.
- Record pending/running/completed/failed states and retry history.
- Resume only when provenance and output integrity match.
- A nonempty file or existing `summary.json` alone does not establish success.
- Capture failures without losing accounting for remaining samples.
- Propagate meaningful nonzero exit status.

Test worker failure, coordinator interruption, incomplete writes, resume, duplicate prevention, and consistency between public and sealed exports.

### Evaluation runner

Inspect and repair `scripts/run_eval_pipeline.py` as necessary:

- Support the manifest’s actual split names while preserving existing CLI compatibility.
- Specify the experiment ledger and data root explicitly; do not inherit an unrelated experiment’s defaults.
- Build evaluations from the frozen inventory, not only successful ledger entries.
- Preserve execution failure, insufficient evidence, ambiguity, and successful reconstruction as distinct outcomes.
- Reject stale successful artifacts from failed or mismatched invocations.
- Generate `evaluation_report.json` and `failure_atlas.md`, plus machine-readable diagnostics.
- Label inferred failure stages as hypotheses unless supported by intermediate evidence.
- Ensure the atlas agrees with strict evaluator semantics.

Discover current arguments with `--help`. Implement and test missing interfaces before using them. Record every exact command used.

Run the complete primary 500-sample workload with resumable accounting. Keep development, validation, and test results separately identifiable in any aggregate report.

## 8. Objective 2 — Scientific deconvolution and debugging

Use this loop for every scientific change:

1. Reproduce a concrete failure.
2. Trace evidence through the pipeline.
3. Identify the earliest demonstrably incorrect stage.
4. State a falsifiable hypothesis.
5. Add a minimal regression test and negative control.
6. Make one coherent change.
7. Measure paired before/after results and ablations.
8. Accept, revise, or reject based on evidence.

Inspect reads, usable depth, mappings, candidate lengths, allele assignments, VCFs, consensus, trimming, classification, and evaluation. Do not infer correctness from the rendered report.

### Required challenge matrix

Cover:

- VNTR lengths spanning 20–140 repeats.
- Severe asymmetry, especially 25 versus 140.
- Near-equal lengths, especially 45 versus 47 and C3 examples.
- Homozygous and sequence-identical alleles.
- `c.530dupC`-related cases, SNVs, deletions, and complex indels.
- C7 low depth: 10×–15×, including unequal allele support.
- Platform-specific errors and mutation-negative controls.

If required cases are absent from the fixed 250-design manifest, add clearly labeled supplementary fixtures or development experiments. Do not silently alter the primary inventory.

### ONT NanoSim indel smear

Investigate the brittle support threshold of 10 in its actual context.

Evaluate adaptive clustering or density-based length inference using usable depth, uncertainty, mapping quality, and platform error structure. Test for false splitting, merging, outliers, minor-allele loss, and spurious additional alleles.

Do not substitute another arbitrary threshold or tune against held-out truth.

### C3 near-equal allele deconvolution

Compare informative k-mers and repeat-unit graph evidence where justified.

- Require read linkage for phase claims.
- Measure phase consistency and switch errors where truth permits.
- Control noise-induced graph branches and k-mer artifacts.
- Preserve ambiguity when reads cannot distinguish haplotypes.
- Do not force two distinct sequences from unsupported evidence.

### C7 low-depth recovery

Separate requested simulator depth from achieved usable depth and per-allele support.

Evaluate sensitivity, false positives, allele dropout, and no-call rates together. Improve recovery only when evidence supports it; retain explicit insufficient-evidence outcomes.

### Consensus and scientific contracts

Preserve and test:

- 1-based dictionary coordinates versus 0-based Python intervals.
- Insertion anchors and inclusive deletion endpoints.
- Canonical repeat counts versus fixed flanking repeat blocks.
- Explicit sample selection and IUPAC ambiguity in consensus.
- Multi-site haplotype claims only when phase evidence supports them.
- Existing CLI preset overrides and platform-specific tool selection.
- Reference, dictionary, trimming, and consensus provenance.

Do not arbitrarily resolve IUPAC bases. Unphased multi-site ambiguity must not become fabricated haplotypes.

Repository documentation records experimental read-backed phasing and strict segmentation approaches that previously failed development gates. Reassess evidence before reusing them; do not advertise them as established accuracy fixes.

### HGVS 20.05

Implement and validate the specifically requested HGVS **20.05** target using the archived authoritative recommendations. Do not describe it as the latest standard.

- Verify reference accession/version, strand, coordinate system, and sequence context.
- Apply the relevant 3′ normalization and duplication/insertion rules.
- Distinguish legacy dictionary labels from normalized HGVS output.
- Treat `c.530dupC` as a case requiring reference-aware validation, not a universal mapping for every affected repeat.
- Preserve backward-compatible aliases where necessary.
- If reference mapping is unresolved, report that limitation instead of fabricating transcript coordinates.

Keep literal sequence equality, exact event annotation, and normalized biological-event equivalence separate.

### Evaluation integrity

Freeze the evaluator and acceptance criteria before algorithm selection.

- Tune on development data.
- Use validation for candidate selection under a documented protocol.
- Freeze the selected implementation/configuration before final test evaluation.
- Keep all modes of the same biological design in the same split.
- If test data influence subsequent changes, mark them as development-exposed and obtain a new independent holdout for confirmatory claims.
- Generate sealed truth without exposing it to production inference or candidate selection.

Report strict allele recovery, literal sequence equality, ordered structure, repeat-count errors, mutation precision/recall, extra/missing alleles, phase ambiguity, no-calls, execution failures, runtime, and peak memory.

Respect conservative bounds across equally optimal allele assignments. Unsupported duplicated sequences do not establish independent recovery of both alleles. Failed normal samples are not true negatives.

Stratify by split, acquisition mode, category, mutation class, length, asymmetry, and achieved depth. Treat biological design as the paired unit when estimating uncertainty across modes.

Do not invent performance targets. Specify evidence-based acceptance criteria before examining held-out results. Report unresolved limitations honestly.

## 9. Objective 3 — Clinical report redesign

Refactor:

- `src/muc_one_span/templates/report.html.j2`
- `src/muc_one_span/templates/report.css`
- Relevant report-generation, asset, and IGV modules.

Preserve the 649-line limit by splitting templates, styles, and scripts by responsibility.

### Decision hierarchy

Create a prominent header supporting:

- PATHOGENIC
- BENIGN
- INCONCLUSIVE

Define and test an explicit evidence-to-header decision table.

These are scientific assertions:

- PATHOGENIC requires a supported finding and traceable classification evidence.
- BENIGN requires defensible benign interpretation; absence of a detected mutation is not enough.
- INCONCLUSIVE covers insufficient coverage, unresolved reconstruction, unsupported classification, and other material uncertainty.
- Show technical execution status separately.
- Explain negative findings within the tested scope without claiming exclusion of ADTKD-MUC1.
- Synthetic truth annotations must never supply a production clinical classification.

Display sample identity, acquisition mode, analysis version, evidence summary, limitations, and provenance clearly.

Do not present dictionary-fit scores or exact-match percentages as calibrated clinical probabilities.

### Evidence visualization

Include:

- Visual repeat alignment and an accessible diff matrix.
- Explicit repeat indices and allele identities.
- Distinct encodings for substitutions, insertions, deletions, ambiguous bases, and missing evidence.
- Text and symbols alongside color.
- Progressive disclosure from summary to supporting evidence.
- Accessible tables/text alternatives for canvas or IGV content.

### Accessibility and print

- Target WCAG 2.1 AAA text contrast: at least 7:1 for normal text and 4.5:1 for qualifying large text.
- Validate actual foreground/background combinations in all states.
- Include keyboard operation, visible focus, semantic headings, accessible names, suitable tab behavior, reduced motion, and zoom/reflow checks.
- Use `font-variant-numeric: lining-nums tabular-nums`.
- Handle 320px–3840px viewports; confine necessary horizontal scrolling to complex data regions.
- Provide print CSS with repeated table headings, sensible page breaks, readable monochrome output, and visible uncertainty.
- Ensure printed output preserves evidence otherwise hidden in collapsed panels or interactive views.

Contrast checks alone do not establish full WCAG AAA conformance. Record which criteria received automated, visual, and manual assessment.

### Offline, sandboxed IGV

Deliver a self-contained report with:

- Effective `loadDefaultGenomes: false`.
- Embedded/local reference, tracks, indexes, fonts, scripts, styles, and required assets.
- No CDN, remote genome lookup, search service, telemetry, or runtime network dependency.
- A justified iframe sandbox and content-security policy where compatible.
- Escaped untrusted sample names, annotations, and JSON/script payloads.
- Usable fallback content when browser functionality is unavailable.

Disabling default genomes is necessary but does not prove offline behavior. Test the exact packaged report.

“Offline” applies to generated report execution; authorized development downloads and Claude reviews are separate activities.

## 10. Objective 4 — Playwright stress and monkey testing

Build a maintained, reproducible headless suite integrated into documented repository commands.

### Network enforcement

Install interception and event listeners before navigation.

- Block service workers.
- Observe the full browser context, including frames, popups, and workers.
- Record attempted external HTTP/HTTPS requests, then abort them.
- Fail if the attempted external request count is nonzero.
- Do not count only completed requests: blocked attempts still violate the requirement.
- Detect and prohibit unexpected WebSocket connections and other browser network paths.
- Define any loopback test-server allowance narrowly.
- Separately test the distributable offline artifact, preferably via `file://`, where HTTP/HTTPS requests must be zero.
- Do not disable browser security to make the tests pass.

Assert IGV reference/tracks actually load and render while offline. An empty or disabled viewer is not success.

### Deterministic stress testing

Use fixed, recorded random seeds and a bounded action budget.

Exercise:

- Rapid repeated clicks.
- Tab switches and panel expansion.
- Keyboard navigation and focus restoration.
- IGV pan, zoom, and track interactions.
- Viewport resizing across 320px–3840px.
- Representative long-repeat, asymmetric, low-depth, ambiguous, positive, benign, no-call, and failed-run fixtures.
- Print media emulation and PDF generation.
- Reopening reports without cache.

Assert:

- Zero external network attempts.
- Zero console errors.
- Zero uncaught exceptions and unhandled promise rejections.
- No blank primary report content or lost decision header.
- No trapped focus or unusable controls.
- No unexpected whole-page horizontal overflow.
- No sample/evidence mismatch after interactions.

Use readiness assertions and explicit IGV load completion; avoid relying solely on sleeps or `networkidle`.

Save failing seeds, action logs, traces, screenshots, network logs, and relevant report provenance under ignored artifact directories. Convert reproducible failures into focused regression tests.

Perform visual inspection of representative screenshots and printed pages. DOM assertions cannot establish visual legibility or clinical clarity.

## 11. Objective 5 — Claude “Fable 5.1 Red Team”

Use `/home/bernt-popp/.local/bin/claude` strictly for brainstorming and adversarial review through noninteractive `claude -p` invocations.

First verify the installed CLI, authentication, supported model selection, and whether a model identified as “Fable 5.1” is actually available.

“Fable 5.1 Red Team” is the requested reviewer identity, not proof of the underlying model. Record the actual model when verifiable. Do not invent a model identifier or silently claim equivalence. If unavailable, record the requested review as blocked and label any available-model substitute accurately.

For a short static prompt, the requested invocation pattern is:

    echo "Act as the Fable 5.1 Red Team. Review the supplied design for falsifiable genomics, systems, and clinical UX failure modes. Return evidence, counterexamples, and required tests. Do not modify files or execute commands." | /home/bernt-popp/.local/bin/claude -p

For substantive reviews, write a carefully scoped packet and pass it via stdin:

    /home/bernt-popp/.local/bin/claude -p < .planning/reviews/milestone-packet.md > .planning/reviews/milestone-review.md

Use verified CLI controls to enforce review-only operation where available. Do not allow Claude to implement changes or launch experiments. Check exit status and nonempty output; use bounded retries/timeouts.

Do not include credentials, real patient data, raw reads, or held-out truth. Send relevant code, synthetic examples, design choices, and aggregate evidence only.

### Mandatory review milestones

1. Alpha/Beta/Gamma design and acceptance protocol.
2. Parallel scheduler, ledger durability, and resume behavior.
3. Each proposed scientific algorithm change before promotion.
4. Clinical decision semantics, accessibility, and offline architecture.
5. Playwright harness and observed results.
6. Final diff and scientific evidence summary.

Require each review to examine:

- Genomics: identifiability, phase evidence, coordinate errors, false positives, truth leakage, and misleading metrics.
- Systems: races, partial writes, stale artifacts, resource oversubscription, portability, and error propagation.
- UX: false reassurance, hidden uncertainty, accessibility failures, print loss, script injection, and network escape.

Ask for severity, precise evidence, a counterexample, a reproducible test, and a proposed remedy. Require the reviewer to distinguish observations from speculation.

Maintain an issue table with accepted/rejected findings, rationale, test evidence, and resolution status. Investigate findings independently; a model review is not ground truth.

Re-review materially changed areas. After three unproductive cycles on the same issue, change the hypothesis or document a concrete blocker rather than repeating the same prompt.

## 12. Execution milestones

### M0 — Preflight, scope, and design

Deliver environment inventory, repository baseline, explicit dataset matrix, Alpha/Beta/Gamma comparison, acceptance protocol, and milestone plan.

### M1 — Reliable experiment infrastructure

Implement and test concurrency, durable state, inventory accounting, split compatibility, and provenance-aware resume. Validate with a pilot across all three modes before scale-up.

### M2 — Baseline and development fixes

Preserve a reproducible baseline revision. Generate the primary inventory, benchmark development data, build the failure atlas, and investigate the required challenge cases.

Implement focused scientific changes with negative controls and ablations.

### M3 — Report and browser verification

Implement evidence-driven clinical reporting, offline IGV, accessibility, print behavior, and deterministic browser stress tests. This may proceed alongside M2 when file ownership and report contracts are stable.

### M4 — Frozen evaluation

Select using development/validation evidence, freeze code and configuration, and evaluate the held-out test split.

Complete accounting for all 500 primary datasets. Publish split-separated and combined artifacts without presenting tuned development results as independent validation.

### M5 — Final review and delivery

Resolve material review findings, run required checks, inspect the final diff, update documentation/changelog, and deliver an evidence-backed summary.

## 13. Completion criteria and reporting

Maintain a live checklist and concise progress updates. For long-running work, report completed/failed/pending counts, throughput, resource usage, and the next verification step.

Every command outcome must be classified as passed, failed, skipped, blocked, or still running. Include relevant exit codes and test counts. A skipped integration test does not validate its behavior.

Before claiming completion:

- All primary sample IDs are accounted for exactly once.
- All 500 intended simulation/evaluation attempts have durable status; failures remain visible.
- `evaluation_report.json` and `failure_atlas.md` match the frozen inventory and strict evaluator.
- Required scientific edge cases have measured outcomes and regression evidence.
- Clinical decision rules, ambiguity, and negative-result semantics are tested.
- Offline report tests demonstrate zero attempted external requests with working IGV.
- Stress tests demonstrate zero console errors and uncaught exceptions.
- Required quality, coverage, integration, documentation, security, and packaging checks have actually run as applicable.
- Review findings are resolved or explicitly documented.
- Generated data and machine-specific configuration are absent from the staged diff.
- Final changes comply with the maximum 649-line rule.

Distinguish:
1. Work completed and verified.
2. Scientific limitations that remain.
3. Failed acceptance criteria.
4. Environmental or review blockers.

Do not equate experiment completion with perfect reconstruction or clinical validation.

If the session is interrupted, write a continuation record under `.planning/` containing the revision, dirty-state summary, active process IDs, exact commands, artifact paths, inventory/configuration hashes, completed milestones, unresolved findings, and the next action. Resume from verified state rather than rerunning blindly.

Final response:

1. What changed and why.
2. Primary/supplemental dataset accounting.
3. Baseline versus candidate results, separated by split.
4. Actual commands and pass/fail/skip results.
5. Browser, offline, accessibility, and print evidence.
6. Adversarial review outcomes and actual reviewer model.
7. Artifact paths and remaining limitations.

Begin now with repository and environment inspection, then proceed through the milestones.
````

## Why this design works

- **Autonomy with boundaries:** It authorizes routine engineering decisions while requiring concrete evidence before completion.
- **Repository-aware planning:** It addresses the verified 500-versus-750 ambiguity and mismatched split interfaces before expensive execution.
- **Scientific integrity:** It separates development from held-out evaluation, preserves no-calls and failures, and prevents unsupported claims of two-allele recovery.
- **Crash-safe execution:** It distinguishes atomic file replacement from consistency across multiple ledgers.
- **Version-specific nomenclature:** It explicitly targets the requested archived standard. [HGVS versioning](https://archive.hgvs-nomenclature.org/versioning/)
- **Measurable accessibility:** It specifies AAA text-contrast thresholds without treating them as proof of full conformance. [W3C contrast guidance](https://www.w3.org/WAI/WCAG21/Understanding/contrast-enhanced.html)
- **Verified offline behavior:** It treats `loadDefaultGenomes: false` as one configuration setting and requires interception of attempted requests, including service-worker considerations. [IGV configuration](https://igv.org/doc/igvjs/Browser-Creation/), [Playwright networking](https://playwright.dev/python/docs/network)
- **Adversarial review with accountability:** It requires reproducible findings, independently checked remedies, and honest identification of the reviewer model.