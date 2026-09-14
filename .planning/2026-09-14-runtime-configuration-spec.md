# Runtime configuration extension before final freeze

User steering: “be sure no hardcoding, everything config driven.” This extends
the current implementation milestone before any final seed is used. Preserve
scientific defaults and the explicit failed-promotion decisions. Configuration
is a reproducibility interface, not license to retune final validation.

## Contract

Central frozen typed settings cover operational defaults and heuristic tunables.
JSON schema1 configuration rejects unknown fields, duplicate keys, wrong types,
nonfinite numbers and invalid ranges before tools execute. CLI values override
file settings, which override central defaults. Relative resource/model paths
resolve relative to the configuration file. Existing CLI flags and Python APIs
remain compatible. Effective settings and input configuration hashes are persisted.
No source branch depends on validation sample name, seed, platform truth or expected
mutation. Tool installation selection continues through PATH and explicit model
paths, without machine-specific defaults or shell interpolation.

Algorithm/format contracts remain fixed: codon size3, half-open/one-based coordinate
conversions, exact unit-cost Levenshtein recurrence/ties, VCF genotype indexing,
one-to-one evaluation, failure semantics and supported-event identity requirements.
These are not scientific tuning knobs. Evaluation endpoints/tolerances are frozen
in the protocol; making them configurable would undermine acceptance criteria.

## Settings sections and default values

- `run`: threads4, platformhifi, min_coverage10, min_qual5, clair3_model empty
  (existing tool default), reference null, minimap2_preset null, reportfalse.
- `allele_selection`: min_gap5, valley_min_points3, valley_min_separation3,
  refinement_max_shift1.
- `classification`: max_indel_probe30, max_fit_edit_distance3,
  novel_repeat_edit_distance2, minimum_unit_fraction0.5,
  early_stop_edit_distance1, strict_segmentationfalse.
- `consensus`: flank_length500, anchor_bases20, anchor_tolerance50.
- `confidence`: qual_low5, qual_high20, weight_below0.3, weight_low0.5,
  weight_high1, absent_weight0.3, boundary_repeats3, boundary_penalty0.5.
- `calling`: sample_name`sample`, read_phasefalse. The latter is experimental
  and remains disabled in the frozen default validation configuration. No
  read-support threshold is invented; optional WhatsHap coverage/MAPQ overrides
  remain null unless explicitly requested through library or JSON settings.
- `read_phasing`: internal_downsampling null, mapping_quality null, preserving
  installed-tool defaults; actual command, version and overrides are recorded.
- `reference_layout`: ordered pre IDs1,2,3,4,5 and after IDs6,7,8,9. Derive the
  fixed count and terminal anchor IDs from actual selected layout, not category
  lists containing alternatives. Validate selected IDs against the dictionary.
  `min_units: 1` and `max_units: 150` configure ladder generation bounds only.
- `repeat_dictionary`: optional resource path; bundled dictionary is default.

All these numeric defaults are named in one typed settings module. Existing
function keyword arguments continue to override their corresponding setting;
inner routines receive immutable settings explicitly. No module-global mutable
runtime state or truth coupling. Diagnostic read-evidence primitives already
accept explicit anchor sequences/edit tolerances; their inputs remain explicit.

## Executable tasks and ownership

1. Add `settings.py` (or split loader by responsibility) and deterministic
   `test_runtime_settings.py`: red unknown/type/range/duplicate/path cases first;
   implement immutable dataclasses and JSON loader; verify defaults and roundtrip.
2. Classification worker: thread classification/confidence settings through
   `classify.py` and `classification_summary.py`, preserving legacy helper aliases
   and explicit boundary arguments. Test at least one observable nondefault per
   category plus complete default-output parity against the preconfiguration
   snapshot. Never weaken old unsupported-confidence assertions; supply verified
   absence fixtures when testing the unchanged penalty.
3. Coordinator: integrate global JSON `--config`, command defaults/explicit CLI
   precedence, configured dictionary/layout, consensus/length/calling settings,
   and effective settings output. Add CLI failure/precedence and stage-forwarding
   tests. Keep every authored file below650lines by splitting helpers/modules.
4. Add optional read-phaser argument overrides with mocked exact command tests;
   no default flag additions or phasing-policy changes. Document experimental
   exposure and implicit tool defaults accurately.
5. Run all gates and complete default classification/length/consensus equivalence
   checks; review the exact tracked/untracked configuration diff with actual
   Claude Fable5.1. Update spec, plan, docs, changelog and final freeze manifests.

Default parity is required for sequence/event/count outputs relative to the
review-fixed candidate before this extension. Additive configuration provenance
is compared separately. Nondefault settings are reproducible user choices and
are not claimed scientifically validated outside the stated fixtures/panel.
Final128-input validation still uses original settings, including disabled phase,
minimum coverage10 and QUAL5; no reserved seed is consumed until this work passes.


## Prefreeze review corrections and compatibility decisions

The independent configuration audit reproduced ignored standalone dictionaries,
invalid explicit stage overrides, stale rerun provenance, incorrect anchor search
coordinates and short-flank anchor mismatches. These are correctness fixes with
failing regressions, not threshold tuning against final data. Standalone consensus
now loads the bundled or configured dictionary and uses the same anchor policy as
the full run. Its prior fixed-only trimming behavior changes explicitly.

Default parity remains required for scoring and length selection. Anchor fixes
must instead be compared on cached development consensuses and reported honestly
if they correct prior outputs; their default differences are not hidden as a
configuration refactor. The existing reference ladder still uses its original
left-flank prefix; changing the biological flank resource is a separate unmet
method. No frozen or fresh evaluation data has been consumed at this decision.
