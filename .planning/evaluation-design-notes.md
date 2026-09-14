# Evaluation design notes — 2026-09-14

Read-only investigation of `scripts/benchmark.py`, `scripts/batch_analyze.py`, the
ignored `results/classification-audit/evaluate_pipeline.py` artifact, package result
schemas, repository guidance, and `.planning/2026-09-14-truth-audit.md`. No production
files changed and no scientific performance experiment or tests executed here.

## Observed problems

- The timing benchmark discards classifications, supports BAM only, defaults to
  HiFi, and returns success on a missing input directory.
- Batch scoring uses substring mutation names and sample-wide candidates. It
  counts wrong-site and unsupported calls as TP, wrong events as TP_partial,
  ignores extra calls when one expected call exists, and treats absent truth as
  mutant. Discovery and execution disagree on FASTQ suffix support.
- The artifact evaluator improves on this through global sequence assignment,
  Counter-based event comparisons, and missing/extra accounting. It is a prototype:
  absent truth silently removes samples; malformed files abort the entire run;
  assertions perform validation; tied assignment scores use the first permutation;
  failed executions can still contribute apparent accuracy through leftover
  summaries; all-counts-within3 permits extra predicted alleles. Its tuple equality
  `(repeat_index, closest_type, name)` is event annotation equality, not normalized
  sequence event equivalence or base localization.
- `summary.json` exposes condensed structures and mutation differences; complete
  classifications live in the allele-keyed `repeats.json`. Consensus FASTAs are
  already trimmed by the pipeline;
  do not apply another fixed trim. Truth FASTAs use the full biological flanks.

## Recommended boundaries and interfaces

Keep evaluation reusable and tool-independent in `src/muc_one_span/evaluation/`.
Each module should remain comfortably below the 649-line gate; do not expand the
two scripts into parallel implementations.

| Module | Public responsibilities and illustrative interfaces |
| --- | --- |
| `models.py` | Frozen annotated records: `TruthSample`, `TruthHaplotype`, `TruthEvent`, `PredictedAllele`, `RunObservation`, `Assignment`, `SampleEvaluation`; explicit coordinate/count conventions and enums for execution/truth/metric states. |
| `truth.py` | `load_truth(sample_dir, repeat_dict) -> TruthSample`; strict legacy MucOneUp adapter and optional authoritative manifest adapter; reconstruct and validate actual sequences, count convention, events, flanks, and provenance. Raise contextual `TruthValidationError`. |
| `artifacts.py` | `load_observation(result_dir, run_record=None) -> RunObservation`; parse summary and per-allele artifacts, identify missing/malformed artifacts, distinguish failed/no-call/incomplete/unknown-provenance states. Sample inventory/discovery helpers can live here initially. No subprocess execution. |
| `matching.py` | `sequence_distance(a, b) -> int`, `match_alleles(truth, predictions) -> Assignment`; pure global sequence distance, maximum-cardinality/minimum-cost injective matching, all optimal assignment alternatives or an equivalent ambiguity representation. Explicit unmatched IDs. |
| `scoring.py` | `evaluate_sample(truth, observation, policy) -> SampleEvaluation`; repeat count, sequence, ordered structure, broad sample detection, and one-to-one event annotation metrics. `aggregate(samples) -> EvaluationReport`; denominator and uncertainty rules are centralized here. Split event sequence normalization into its own module only when sufficient inputs are implemented. |
| `__init__.py` | Narrow public entry points; no environment checks or tool discovery on import. |

Add thin `scripts/evaluate_pipeline.py RESULT_ROOT --truth-root ... --output ...
--expected-samples ...`. It performs inventory, per-sample error handling, report
writing, and exit status; computation stays in the package. Preserve the useful
artifact CLI option alias `--truth-dir`. Default output should be a documented
ignored directory rather than the ad hoc classification-audit location.

Keep timing and running concerns in existing scripts. Reuse one input-discovery
helper for BAM, `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz`, and give deterministic
ambiguity errors rather than choosing an arbitrary file. Preserve current script
flags and positional arguments; use the package tool abstraction for commands.
Batch evaluation should call the same scoring API, with accurate statuses; do
not keep substring scoring as the primary advertised metric. If compatibility
requires legacy categories, label them explicitly and separately.

## Required scientific and failure contracts

1. **Inventory and completion:** expected sample manifest is the denominator
   authority. Without one, union discovered truth, observations, and run records;
   report that this is discovery-based. Reject duplicate sample IDs and conflicting
   run records. Every inventoried sample gets a row, including missing truth,
   missing summary, no reads, timeout, nonzero exit, malformed artifacts, and no-call.
   Write a report even when individual samples fail; exit nonzero for execution,
   validation, or input-inventory errors. Biological inaccuracy alone need not make
   the CLI fail unless explicit acceptance thresholds were requested.
2. **Truth:** require unique haplotype IDs agreeing across FASTA/structure/stats,
   nonempty chains, known dictionary symbols, and mutation markers consistent
   with mutation metadata and actual mutated-unit FASTA. Reconstruct each complete
   truth sequence exactly before scoring. Validate total unit count and actual
   VNTR sequence length. Historical stale `repeat_lengths`, null metadata seeds,
   and misleading coverage fields should produce named provenance warnings when
   reconstruction succeeds, not silently override actual truth or reject all
   known valid historical data. Unknown missing mutation sequence is invalid truth.
3. **Counts:** simulator `repeat_count` compares directly to reported `length`;
   verify `canonical_repeats + 9 == length` as a separate invariant. Report exact,
   within-1, and within-2 counts, signed errors, missing alleles, extra alleles,
   and exact diploid reconstruction requiring the right cardinality.
4. **Assignment:** use complete VNTR sequence, not observed mutation labels, to
   select the optimum. Preserve missing and extra alleles. Empty prediction sets
   are legitimate no-call observations; zero matching cost is not exact recovery.
   Deterministic presentation order must not imply biological phase resolution.
   For tied optima, report phase ambiguity and metric ranges where outcomes vary;
   only assert exact phased success if it is established across alternatives.
   Identical biological haplotypes can still have invariant sequence accuracy.
5. **Sequences/structures:** report exactness, edit distance, sequence length error,
   and ambiguous base counts; literal IUPAC symbols must not match A/C/G/T as
   wildcards. Define identity denominator and empty-sequence behavior. Compare
   ordered unit tokens, not histograms. Report dictionary sequence-equivalent
   symbols separately if implemented; never silently collapse nomenclature.
6. **Events:** exact canonical names or an explicit vetted alias table, never
   substrings. One-to-one multiset matching preserves duplicate false calls.
   Unmatched truth is FN and every unmatched prediction is FP, including on
   unmatched alleles and positive samples. Score all candidates and an explicit
   supported rule separately (e.g. frameshift + exact template + VCF support).
   Report the rule in metadata. Wrong event means expected-event FN plus called-event
   FP; broad sample alarm detection remains a separate metric. Keep tuple-based
   site/name agreement named `event_annotation_exact`; only claim biological event
   equivalence after normalized alternate sequence/context and coordinate mapping
   are implemented. Missing precise inputs yield `not_assessable`, not false TP.
7. **Execution versus accuracy:** nonzero exit or incomplete artifacts must not
   produce successful sample/genotype status from stale files. If known valid truth
   exists, missed output counts as failure in end-to-end reconstruction and expected
   events remain unrecovered; conditional accuracy can be shown with its narrower
   denominator. Missing/invalid truth is unassessable and stays visible in completion
   totals. An unexecuted normal sample is not a TN. Report numerator/denominator
   pairs and null for undefined ratios rather than NaN or misleading zero.
8. **Provenance:** version report schema and scoring policy; record package version,
   truth and input hashes, sample inventory, platform, run exit/timing metadata,
   observed software/model versions, available seeds, and known missing provenance.
   Old results without a run manifest have unknown execution provenance; never
   fabricate successful exit status. All-counts success excludes extra alleles.

## Focused deterministic test files

- `test_evaluation_truth.py`: valid tiny normal/mutant fixtures; differing true
  unit lengths; total/canonical convention; duplicate FASTA IDs; missing or multiple
  candidate files; missing haplotype; unknown symbol; marker/metadata discrepancy;
  mutated-unit mismatch; invalid flank; empty flank slicing; stale historical
  lengths warning; actual truth hash changes. No generated sequencing datasets.
- `test_evaluation_matching.py`: swapped alleles; equal length/different sequence;
  identical haplotypes; phase tie with distinct mutation-score outcomes; no
  predictions; one dropout; extra third prediction; no truth error; symmetry and
  insertion/deletion/substitution distance checks against simple trusted DP on
  seeded short strings and unit-token sequences. Keep random oracle checks in tests,
  not CLI startup. Small hand-checked long repeat fixture for bit-vector edge cases.
- `test_evaluation_scoring.py`: wrong allele/site/type; substring collision dupC
  versus dupCCCC; one expected plus one extra; duplicate calls; mutant on haplotype
  2; multiple mutations; unsupported and nonframeshift candidates; unmatched-allele
  events; tied assignment ranges; exact cardinality; normal true negative versus
  failed normal; empty denominators; missing precise event data yields unassessable.
- `test_evaluation_artifacts.py`: missing/malformed summary, missing/multiple-record
  consensus, mismatched allele/classification keys, empty classifications no-call,
  stale summary with failed run record, unknown provenance, duplicate inventory IDs,
  all input suffixes, ambiguous input selection, mixed platform metadata.
- `test_evaluation_cli.py`: valid minimal fixture yields stable versioned JSON;
  missing root/zero eligible samples/failing sample exit nonzero; report includes
  failures and missing truth; one malformed sample does not hide others; output
  parent creation; expected-sample manifest keeps unattempted samples in totals;
  script help and retained option forms; runtime subprocesses mocked at boundary.

Validation for implementation: focused tests first, then `make ci-check`,
`make docs-check` for benchmark documentation, and a representative evaluation of
existing ignored artifacts with expected sample inventory. `make test-int` is
necessary if changing the running scripts' bioinformatics execution behavior;
offline evaluation itself needs no external tool. Existing numerical results
should be compared explicitly, explaining corrected denominators and tie handling
rather than treating an intentional metric change as a regression.
