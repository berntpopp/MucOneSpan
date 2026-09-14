# Production improvements and scientific validation

Status: all local release software gates pass. Actual Claude Fable 5.1 final
release review found no release-blocking correctness defects; all eight
non-blocking findings have explicit dispositions. The user scoped this release
to verified improvements and a future experiment workflow. Fresh final paired
validation was not performed, and zero reserved seeds have been used. No general
scientific accuracy improvement is claimed. The release tag identifies the
committed candidate; publication is verified through the repository workflows.

## Revisions and reproducibility

Baseline MucOneSpan 0.10.0, `d8390b3c244ef8f3240af74b92db12b50dfc77d1`.
The release candidate is maintained in `.worktrees/production-validation` on
`improve/production-validation`. Original user changes are preserved. Release is now explicitly authorized;
no merge into main is planned. Runtime source is isolated; shared
ignored reads/results are accessed explicitly, never committed.

Locked uv Python 3.12.9; minimap2 2.28-r1209, samtools 1.15.1, bcftools 1.17,
Clair3 v1.0.10, optional WhatsHap 1.7. Model directories are separately fingerprinted
for HiFi and ONT. New simulator MucOneUp 0.44.5 and error-model configuration are
documented in [fresh-generation-design.md](fresh-generation-design.md). Historic
44-sample data used an older simulator; direct paired comparisons always use the
same input bytes for both callers. Provenance evidence: [evidence-inventory.md](evidence-inventory.md).

## Implemented production contracts

- Exact integer bit-vector Levenshtein scoring preserves literal symbol equality,
  first-match ties and unchanged traceback. It changes computational cost only.
- Strict offline evaluation reconstructs actual truth, assigns alleles one-to-one,
  preserves all equally optimal assignments, and counts incorrect/extra event
  names, parents, total-repeat positions and haplotype assignments separately.
- Net frameshift uses signed insertion-minus-deletion. Classification exposes
  unresolved sequence and heuristic confidence; no probability calibration is implied.
- Variant support requires exact selected-VCF sequence concordance at a verified
  reference/consensus projection. This is agreement with the input VCF, not
  independent read evidence. Unavailable projection has an explicit reason.
- Consensus explicitly selects sample and genotype/IUPAC behavior. Absent variants
  do not establish sequence homozygosity. Unknown evidence states remain ambiguous.
- Failure sidecars and explicit invocation records distinguish coverage no-calls,
  tool failures, incomplete/malformed outputs and completed but poor predictions.
- Alignment fit counts, primary records, discarded length candidates, unresolved
  aliases and experimental phase read selection are named separately. QNAME
  collisions never justify collapsing distinct observed records. Selected layout
  IDs determine `fixed_repeat_count`; the evaluator reports
  `canonical_plus_fixed_matches_reported`, retaining the +9 diagnostic only for
  legacy comparison.
- Immutable schema-1 runtime configuration centralizes operational defaults and
  scientific heuristics. Explicit CLI values override file settings, then central
  defaults. Unknown/duplicate fields, wrong types and invalid bounds fail early;
  relative resource paths resolve from the file. Effective settings and
  input/configuration/reference/dictionary hashes are recorded before tools, and
  failed reruns invalidate stale configuration provenance.

See [requirement-evidence-matrix.md](requirement-evidence-matrix.md) for exact test
references and bounded versus unmet workstreams. Public flags remain; additive
metadata clarifies semantics. Benchmark batch output migrates from a bare array
to explicit schema 1 with execution and evaluation status, documented in the guide.

## Configuration evidence before the final freeze

The strict settings suite covers 91 cases; the CLI configuration suite covers
20 cases, including precedence, invalid explicit overrides, failed-rerun
provenance and standalone dictionary forwarding. Named tests include
`test_bad_json_configuration_is_rejected`,
`test_reference_layout_range_roundtrip_and_validation`,
`test_cli_overrides_file_which_overrides_defaults`,
`test_failed_rerun_does_not_retain_old_configuration`, and
`test_standalone_consensus_passes_configured_dictionary` in
`tests/unit/test_runtime_settings.py` and `tests/unit/test_cli_settings.py`.
`tests/unit/test_api_defaults.py` verifies default signatures remain connected
to central settings; `test_canonical_count_consistency_uses_selected_reference_layout`
checks the generic fixed-repeat diagnostic. Classification/reference/consensus
settings have observable nondefault fixtures and bounded default-equivalence
comparisons. Details: [runtime-settings-evidence.md](runtime-settings-evidence.md),
[classification-settings-evidence.md](classification-settings-evidence.md),
[length-consensus-settings-evidence.md](length-consensus-settings-evidence.md),
[api-defaults-evidence.md](api-defaults-evidence.md), and
[runtime-config-dispositions.md](reviews/runtime-config-dispositions.md).

The software left anchor now matches the actual prefix used by ladder generation;
boundary tolerance and oversized component handling were corrected. Among the
154 expected allele slots in the 77 exposed development inputs, 137 available
cached full-consensus files produced identical old/new trimmed FASTA bytes and
trim coordinates. Only `left_trim_method` changed from `fixed_anchor_not_found`
to `exact_anchor`; 17 missing slots received no parity credit. A synthetic
flanking insertion now correctly removes a previously retained flank base, so
universal old/new output identity is not claimed. Exact fixtures:
`tests/unit/test_consensus_settings.py::{test_short_flank_anchor_uses_same_prefix_as_ladder,test_default_prefix_anchor_corrects_one_base_flank_insertion,test_zero_tolerance_finds_anchors_at_expected_boundaries,test_sliced_anchor_components_preserve_short_repeat_boundaries}`.
Evidence: [consensus-anchor-config-fixes.md](consensus-anchor-config-fixes.md).
The biological distal-flank reference choice remains unchanged and unvalidated
against a regenerated proximal reference.

Configuration does not make arbitrary parameter choices scientifically validated.
Confidence remains heuristic, and coordinate/codon/genotype/edit-distance and
evaluation acceptance contracts remain fixed. Thread counts are per-tool settings
with concurrent stages, not a strict process-tree CPU cap. These fixtures and
cached comparisons precede the fresh panel and do not add new accuracy metrics.

## Controlled classification performance

The isolated scorer produced zero complete-output differences in 141 cached
classification cases. At least 10,000 random/adversarial pairs agree with an
independent scalar oracle, including ambiguous symbols and indels. Three paired
alternating difficult-panel repetitions gave median 17.612304 s baseline versus
2.214009 s candidate: **7.955× classification speedup**. Near-exact speedup 5.026×;
batched exact controls differed by less than 1%. This is not a whole-pipeline claim.
Commands, hashes and limits: [task2-evidence.md](task2-evidence.md).

## Development regression evidence, re-scored consistently

The current strict evaluator re-scored all original 53 full inputs and 24 dependent
perturbations, preserving archived predictions and their actual execution records.
It did not rewrite old outputs to imply the final caller had been rerun. Artifacts:
`tests/results/production_validation_20260914/evaluation_current/`; source hashes,
result hashes, explicit inventories and evaluator exit codes are in `provenance.json`.
`refresh_development.py` is the retained reproduction command. Old source outputs
can lack new projection reasons; those are legacy reason-unavailable observations.

| Stratum | Inputs | Independent exact sequences before→after | Exact diploid sequences before→after | Exact annotation TP/truth before→after | Extra events before→after | Samples with output before→after |
| --- | ---: | --- | --- | --- | --- | --- |
| Original HiFi | 44 | 39/88→39/88 | 10/44→10/44 | 24/28→24/28 | 14→14 | 44/44→44/44 |
| Original ONT | 3 | 3/6→3/6 | 0/3→0/3 | 2/2→2/2 | 0→0 | 3/3→3/3 |
| Exposed later HiFi challenges | 6 | 1/12→1/12 | 0/6→0/6 | 2/3→2/3 | 0→0 | 6/6→6/6 |
| Dependent depth/name perturbations | 24 | 6/48→7/48 | 0/24→0/24 | 6/14→7/14 | 12→11 | 17/24→18/24 |

The perturbation change includes repair of an invalid legacy alias output; it is
not independent evidence of a new mutation-reconstruction algorithm. Six low-depth
outcomes remain no-calls: baseline records label them execution failures, candidate
sidecars distinguish insufficient evidence. None is removed from truth denominators.
The original 44 and exposed challenges show no exact sequence/event aggregate gain.

For original HiFi normal controls, 15/16 had no mutation alarm historically, but ten
of those negatives have unresolved reconstruction. The stricter evidence report
therefore records five confident negatives, one false-alarm control and ten
unresolved negatives: conditional specificity 5/(5+1), not 5/5 or 15/16. ONT has one
normal unresolved negative and no assessable TN/FP pair, so conditional specificity
is undefined, not 100%. Explicit confident-negative/all-normal, alarm and no-call
rates must accompany specificity. Mutation-positive samples' extra calls remain
in event precision even when sample-level detection is successful.

## Failed scientific promotions and research

Lowering minimum coverage 10 → 5 on all 77 development inputs improved exact counts
and output rate, but increased extra events 25 → 38 and normal false-alarm samples
4/30→7/30. All three new normal alarms were HiFi perturbations. Complete diploid
exact sequences stayed 10/77. Default 10 remains. Details and per-sample regressions:
[coverage-parameter-experiments.md](coverage-parameter-experiments.md).

Optional read phase on five collapsed/close-length cases increased exact sequences
1/10→4/10 but introduced an extra supported event. The three gains are majority
short-allele sequences after removing IUPAC ambiguity; each new minority candidate
remains twenty repeats shorter than its assigned truth. All thirteen selected
records in the causal example came from the short-read group; five original long
reads had been omitted upstream. Six WhatsHap parameter choices did not correct
the missing allele. Read phase remains experimental and disabled by default in
both Python and CLI use. Python `read_phase=True` and explicit JSON
`calling.read_phase=true` can opt in; there is no dedicated CLI flag. Null WhatsHap
overrides preserve installed-tool defaults, while an invocation records the exact
phase argv and version. Per-phase-set selected counts are not all-read support.
See [read-phase-config-fixes.md](read-phase-config-fixes.md) and
`tests/unit/test_read_phasing.py::test_explicit_tool_overrides_are_recorded_and_forwarded`.

A minimum selected-read floor is not validated: true 30:2 input support becomes
14:1 after WhatsHap selection; applying floor 2 discards an exactly reconstructed
minority with two original records. Forty-five known-sequence real-tool fixtures
establish this counterexample and the cost of floors 2/3/5. Details:
[read-phasing-cached-validation.md](read-phasing-cached-validation.md).

Error-tolerant anchors increase observed spanning evidence but do not reliably
resolve one-repeat gaps or ONT rejection. Strict segmentation stops lost five
named events in 106 cached outputs and remain opt-in experimental. Global QUAL3,
forced ALT and mutation-template-anywhere rescue fail specificity or phase
contracts. No unsupported heuristic was promoted to repair the aggregate.

The requested actual GeneFoundry MCP/PubTator research retrieved full-text methods
for five relevant papers, with exact citations and retrieval records. Sequence-aware
whole-read grouping, candidate consensus and assignment of all eligible original
reads are promising future designs; the papers do not validate a universal MUC1
support floor. Research: [genefoundry-pubtator-research.md](genefoundry-pubtator-research.md),
[phasing-research-alternatives.md](phasing-research-alternatives.md), and
[phasing-parameter-research.md](phasing-parameter-research.md).

## Fresh evaluation

Not performed under the user-directed scoped release; zero reserved seeds were
generated or used. The design is retained in examples/future-validation-experiment.json
for future work. Before generating it, freeze caller, evaluator, models,
simulator configuration, expected 128 primary inputs and preselected timing repeats.
All 32 main simulations and 96 perturbations remain in denominators. The 32 seed
designs are simulated stress cases, not a population-random clinical cohort;
perturbations and paired alleles are dependent. No retuning on exposed final data.

## Checks and independent review

Final local release checks are recorded under
`tests/results/production_validation_20260914/checks/release/`:

- `make ci-check`: 694 passed, 93.48% branch-aware coverage; Ruff, formatting,
  configured mypy and maximum 649-line gate passed.
- `PATH=... make test-int`: 47 passed, 8 generated-data tests deselected.
- `make docs-check`: passed.
- `make build-check`: wheel installation, CLI and bundled resources passed.
- Explicit generated-data selection: 6 passed and two unchanged strict expected failures (equal 60/60 and asymmetric 25/140).
- `make security-check`: passed for release preparation; no known vulnerabilities
  in the locked dependency audit. Packaging now explicitly excludes generated data
  and validates both source and wheel archives.

Actual Claude Code 2.1.268 explicitly used `claude-fable-5-1` for separate
specification, plan and milestone sessions. Requests, immutable tracked/untracked
manifests, JSON responses and finding dispositions are under [reviews/](reviews/).
The actual Claude Fable 5.1 runtime-configuration review completed; all six findings
have documented fixes or scope dispositions in runtime-configuration-final-dispositions.md.
Final release-diff/report review completed without blocking findings; see
[final dispositions](reviews/final-release-dispositions.md). All 694 unit tests
also pass on Python 3.10, 3.11, 3.13 and 3.14; the primary Python 3.12 run supplies
93.48% branch-aware coverage. Fresh paired evaluation was explicitly not performed
under the user-directed release scope.

## Remaining limitations and rollback

Validated joint assignment of all original reads, one-versus-two biological allele
inference, default read-backed full haplotypes, coverage-aware reference masking,
and sensitivity rescue remain unmet. The software trim anchor now matches the
actual selected flank prefix, with the bounded parity and correction evidence
above. The bundled ladder still selects a biologically distal left-flank segment;
changing that reference requires a separate paired regeneration/ablation. Heterozygous IUPAC indels may leave strict
variant support unavailable even when bcftools applies ALT. Multi-tool-version
support is not established beyond executed local versions.

Rollback is by independent component: retain exact scoring/evaluation/failure
contracts while keeping rejected inference methods disabled. No scientific gate,
test tolerance, coverage threshold or expected failure was weakened. Completed
plan components can be archived only after their evidence and review are final;
unmet scientific requirements remain explicitly listed.


## Final length-failure investigations and future experiment workflow

The equal-length case remains 60/63 under current inference. Of 103 unambiguously
joined reference-fit read pairs, every read scores better against the 60-repeat
candidate than the 63-repeat candidate. The false candidate has little full-boundary
support. The asymmetric 25/140 case excludes all five long original records before
calling; its false second short candidate contains one actual short read. Unique-name
controls reproduce the same assignment, so this is not fixed by renaming alone.

Continuous-span Gaussian and Student-t prototypes improved diagnostic count sets
but still falsely split the equal-length control into 59/60. They were not promoted;
mutation/complete-sequence recovery was not validated. Original strict expected
failures remain. See equal-length-failure-investigation.md and
asymmetric-failure-investigation.md for methods, parameters, hashes and limits.

The maintained JSON experiment runner replaces hardcoded sample definitions with
explicit design files and preserves the 29 historical examples. A 32-case future
HiFi/ONT design is provided without claiming it has been executed. New development
seed 19200902 smoke: both platforms generated diploid 25/30 truth with dupC on haplotype 1, repeat 10 and
passed strict truth checks; 12 requested templates produced 10 HiFi and 12 ONT usable FASTQ
records. Read-source assignment truth remains unavailable through the supported
MucOneUp amplicon CLI. Commands/config/model/input/output hashes and all failure states
are inventoried. The two smoke cases are execution validation, not final scientific
or performance evidence. Public guide: docs/guides/simulation-experiments.md.


## Descriptive uncertainty on original HiFi sample endpoints

| Sample endpoint | Observed | 95% Wilson interval |
| --- | --- | --- |
| Exact count pairs | 40/44 | 78.8%–96.4% |
| Exact diploid sequences | 10/44 | 12.8%–37.0% |
| Normal mutation alarms | 1/16 | 1.1%–28.3% |

These are descriptive binomial intervals on sample outcomes, not calibrated
clinical/generalization bounds. The simulations are a selected development design,
not a random population sample; some seeds/provenance are incomplete, and template,
allele and perturbation dependencies preclude treating all records as independent.
The normal-alarm endpoint is separate from confident-negative specificity: ten
normal outputs remain unresolved. No interval is supplied for undefined ONT
conditional specificity. Formula uses z=1.959963984540054 without data-dependent
selection. Fresh-panel uncertainty remains unmeasured because that panel was not run.
