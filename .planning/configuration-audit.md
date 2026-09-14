# Configuration audit before final freeze

Date: 2026-09-14. Read-only audit of production changes relative to baseline
`d8390b3c244ef8f3240af74b92db12b50dfc77d1`, including untracked new production
modules. No source changes, simulations, caller runs, or settings changes were
made. The coordinator owns runtime configuration design and implementation.

The caller contains no machine-specific installation paths, validation seed or
sample-name branches, or source-truth shortcuts found by this audit. However,
several inherited scientific heuristics and new optional behaviors are only
private/API defaults or literals in routines. A typed runtime configuration
should expose those choices while preserving their current defaults. Fixed
coordinate conventions, genotype contracts, identity checks and algorithmic
recurrences should remain code, not become arbitrary user-adjustable parameters.

## Scope and evidence

Inspected tracked changes in `alleles.py`, `calling.py`, `classify.py`, `cli.py`,
`consensus.py`, `repeat_alignment.py`, `vcf.py` and the report template; new
`read_phasing.py`, `read_evidence.py`, `phasing.py`, `variant_support.py`,
`classification_summary.py`, `benchmarking.py`, `run_status.py`, and the
evaluation package. Read the existing `config.py`, `ladder.py`, `mapping.py`,
`tools.py`, and bundled dictionary to trace configuration and reference semantics.
Compared relevant thresholds directly with `git show d8390b3:<file>`.

Searches included numerical comparisons/defaults, thresholds, coverage, model and
tool selection, `/home/`, `/tmp/`, `/opt/`, development sample names, fresh-data
prefixes, seeds, `source_truth`, and `truth_dir`. Truth directories occur in the
offline benchmarking/evaluation boundary, not scientific caller selection.
The captured WhatsHap 1.7 help in ignored `phasing_debug/whatshap-1.7-help.txt`
was read to distinguish external defaults from wrapper settings.

## Existing configuration that should be retained

| Control | Existing boundary and default | Coverage gap |
| --- | --- | --- |
| Platform / mapping preset | `run`, `map`, `call`: `hifi`; auto `map-hifi` or ONT `lr:hq`; explicit preset override | Centralize the mapping once and retain explicit override precedence. |
| Threads | CLI/API, default 4 | Validate positive integer; retain existing per-allele budget division. |
| Coverage | `run`/`alleles` CLI and `detect_alleles`, default 10 | Preserve 10; the development 5 ablation did not justify a default change. |
| Minimum QUAL | `run`/`call`, default 5; low-level `filter_vcf` default 0 | Preserve both API contracts and pass resolved run value explicitly. |
| Reference | CLI `--reference`, otherwise packaged reference | Record effective path/hash; pair it with correct trimming/reference layout settings. |
| Clair3 model | CLI/API model path, empty means external tool selection/default | Keep explicit configurable model path; no machine-specific model default. |
| Repeat dictionary | JSON-backed `RepeatDictionary`; selected subcommands accept `--repeats-db` | Full `run` currently calls the bundled dictionary directly. |
| Ladder bounds/flanks | `build-ladder`: min 1, max 150 canonical units, flank 500 | Full `run`/consensus needs a coherent reference/flank configuration too. |
| VCF sample selection | Parser API can select a sample explicitly and rejects ambiguity | Clair3 wrapper always names its single sample `sample`. Make that neutral label an operational field if requested. |
| Experimental read phase | Calling API `read_phase=False` | Keep disabled. Configuration may provide explicit opt-in; do not enable during configuration plumbing. |
| Strict segmentation | Classification API `strict_segmentation=False` | Expose as an explicit experimental algorithm option; preserve false. |
| Diagnostic anchors | `read_evidence`: caller-provided anchors, `max_edits=1`; `max_distance` explicitly supplied | Already API-driven; no truth-derived or platform-specific hidden choice. |

`min_dp` is an existing accepted-but-unused compatibility argument. It is not an
effective depth filter, and the current CLI does not offer a `--min-dp` option.
Do not present it as a working runtime knob or silently implement new depth
filtering while moving configuration. Metadata/help should distinguish that
compatibility argument from applied settings.

## Tunable heuristics embedded in routines

These values mostly **predate this branch**. Their age does not make them
scientific constants, and moving them must not alter their comparison operators.

| Location | Current value/behavior | Recommended typed setting |
| --- | --- | --- |
| `alleles._find_clusters` | Private `min_gap=5`; new cluster when difference is **>=** gap | `alleles.cluster_gap=5` |
| `alleles._split_cluster_by_indel` | At least 3 observed indel-series points | `alleles.min_valley_points=3` |
| Same function | Valley separation at least 3 canonical units | `alleles.min_valley_separation=3`; distinct from the preceding setting |
| `alleles._build_allele_info` | Use refined count only within **<=1** of weighted center | `alleles.refined_center_tolerance=1` |
| `classify_sequence` and `_classify_backward` | Duplicated `max_indel_probe=30` | One `classification.max_indel_probe=30` passed through both directions |
| Classification probe/tail logic | Minimum probe/tail span `unit_length // 2` in several paths | One documented `classification.min_unit_fraction=0.5`, with explicit rounding preserving current floor behavior |
| Backward fallback, strict stop, uncertainty labels | Edit distance **>3** stops/marks uncertainty | `classification.max_edit_distance=3`, used consistently in existing paths |
| `classify_repeat` | Substitution-only difference **>2** means novel repeat | `classification.max_variant_edits=2`; do not conflate with segmentation cutoff |
| `consensus._find_anchor` | Search tolerance 50 bases | `consensus.anchor_search_tolerance=50` |
| `consensus.trim_flanking` | 20-base flank and repeat components in both terminal anchors | `consensus.anchor_component_length=20`; derive offsets from it |
| `validate_mutations_against_vcf` | Boundary range 3 repeats, multiplier 0.5 | Existing API fields in `confidence.boundary_repeats=3`, `boundary_penalty=0.5` |
| Same function | Absent exact-support weight 0.3 | `confidence.absent_support_weight=0.3` |
| `_qual_to_confidence` | QUAL below 5: 0.3; 5→20 interpolates 0.5→1.0 | Typed low/high QUAL breakpoints and weights; derive the span 15 and interpolation slope |
| `tools.get_tool_versions` | Version-query timeout 10 seconds | Optional operational `tools.version_timeout_seconds=10`; inherited helper, lower priority than scientific settings |

The confidence curve is a **heuristic weight**, not the Phred probability formula.
Its 5, 20, 0.3 and 0.5 values are tunable policy. The perfect-match/neutral multiplier
1 and percentage conversion 100 are mathematical conventions, not independent
research settings. Display rounding to one/four decimals is presentation behavior.

The private threshold values need API plumbing as well as a loader: merely
reading a configuration file while these routines continue to use their literals
would not satisfy the requirement. Conversely, do not create separate aliases
for the same underlying threshold in forward/backward classification.

## Optional read-phasing tool defaults

The wrapper supplies `--indels`, `--ignore-read-groups`, sample/reference/output
arguments and read-list capture. It does **not** supply coverage or MAPQ knobs.
Captured installed WhatsHap 1.7 help reports `--internal-downsampling` default 15
and `--mapping-quality` default 20. These are external defaults, not hardcoded
MucOneSpan thresholds; option names can vary with tool versions.

The smallest configuration extension is optional typed overrides for supported
internal-downsampling and mapping-quality options, defaulting to `None` to retain
the existing invocation. Record `tool_default` with version/command when omitted;
do not claim a guessed resolved numeric value. Algorithm selection and the
wrapper's deliberate indel/read-group policy can remain documented behavior for
the initial supported adapter, or become enumerated options if explicitly in
scope. Avoid arbitrary free-form tool arguments that bypass genotype/identity
contracts. No configured minimum haplotype-support cutoff is recommended here:
none is validated, and internally selected counts are not all-original-read support.

## Reference-derived values require care

`PRE_AFTER_REPEAT_COUNT=9`, hardcoded terminal IDs `1`/`9`, and the evaluator's
`canonical_plus9_matches_reported` check are tied to the selected ladder layout.
They are **not universal biological constants**. `RepeatDictionary.repeat_length_bp`
already supplies the repeat unit length; mutation templates and their allowed
repeat types already come from data rather than sample-specific code.

Do **not** replace 9 with `len(pre_repeat_ids) + len(after_repeat_ids)`. In the
bundled dictionary these are category lists containing alternatives:

- Pre category: `1,2,3,4,4p,5,5C` (7 entries).
- After category: `6,6p,7,8,9` (5 entries).
- Actual `ladder.build_contig` currently selects `1,2,3,4,5` and `6,7,8,9` (9).

Represent the actual ordered ladder pre/after IDs separately in reference
configuration or reference metadata, with those legacy defaults. Derive the
fixed count and terminal anchor repeat IDs from the actual selected layout.
Do not insert mutually alternative repeat types into one ladder. Custom reference
FASTA, ladder generation, trimming and count diagnostics must share that layout;
otherwise accepting a custom dictionary creates the appearance of configurability
while retaining incompatible assumptions. This adjacent inherited issue should
be explicitly covered or documented as outside the first configuration contract.

## What should remain fixed

- 1-based VCF/repeat coordinates and conversion to Python offsets.
- Codon length 3 in signed net-indel frameshift calculations.
- Diploid two-allele genotype/phase-set contracts and supported haplotype labels.
  Arbitrary ploidy is a different algorithm/schema change, not a threshold edit.
- Myers bit operations, unit-cost Levenshtein recurrence, maximum-cardinality then
  minimum-distance assignment semantics, IUPAC/complement maps and SAM flag definitions.
- Equality/identity checks, zero-division/empty guards, unique-source validation,
  common phase-set requirements, and fail-closed status/evidence contracts.
- The evaluator's named exact/within-one/within-two metrics and all-sample
  denominator policy. Making the benchmark tolerances mutable would change the
  meaning of its metric names; additional metrics need versioned definitions.
- Standard output filenames, schema version identifiers and the source-record
  ordinal naming format, unless a separate output-schema change is requested.

The classifier's `999` edit-distance fallback is a sentinel, not a threshold to
expose. `float("inf")` would communicate that role more clearly in a separate
cleanup. The factor `2 * max_edits` in diagnostic boundary ambiguity follows
opposite endpoint displacement under the same edit bound; it is not a separate
coverage or assignment tolerance to tune independently.

## Smallest maintainable implementation recommendation

Use a standard-library JSON loader and immutable typed dataclasses in a separate
runtime-settings module; keep `config.py`'s repeat-dictionary role clear. Suggested
sections are `run/mapping`, `alleles`, `classification`, `consensus/reference`,
`confidence`, and `read_phasing`. Operational paths/model/sample names belong in
run/reference fields; no environment-specific defaults should be committed.

Maintain one authoritative default object. Existing CLI flags and public API
arguments remain supported. Resolve **explicit CLI override > file setting >
built-in default**; Click's default-valued parameters must not accidentally
override a loaded file. Pass the resolved settings to routines rather than using
mutable globals or re-reading files inside classification loops. Preserve
positional API compatibility and distinguish explicit values from omitted ones.

Validate unknown fields, types (including rejecting bool as int), finite numeric
values, positive threads/lengths, nonnegative edit/coverage settings, probability
weight ranges, ordered QUAL breakpoints, and reference-layout consistency. Resolve
relative data/model paths against the configuration file location and document
the rule. Reject contradictory settings visibly; do not silently clamp them.

Executable selection is already configurable through PATH and environment
activation; Clair3 model selection is already an explicit path. Keep that
mechanism and the existing `run_tool` abstraction. If a `tool_bin` convenience
field is added, scope its environment to the invocation instead of hardcoding
per-tool absolute paths. Resolve/version-check tools using the **same cleaned
PATH** as execution: existing generic `check_tools` checks raw PATH while
`run_tool` cleans virtualenv entries, so recording only the preflight resolution
can misidentify the binary actually used.

Persist the complete effective scientific/operational configuration and its hash,
reference/dictionary/model identities, resolved tool paths/versions, and exact
tool options with each run. Tests should show omitted configuration reproduces
defaults, each exposed field reaches its intended routine, CLI precedence works,
and invalid fields fail. The final validation freeze must be regenerated after
configuration integration; no unseen validation sample should select defaults.
