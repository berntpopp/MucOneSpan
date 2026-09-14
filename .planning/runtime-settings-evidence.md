# Runtime settings validation evidence

Implementation follows `.planning/2026-09-14-runtime-configuration-spec.md` and
coordinator's explicit refinements before the final freeze. Owned changes are
new `src/muc_one_span/settings.py` and `tests/unit/test_runtime_settings.py` only;
CLI and scientific-stage integration belong to the coordinator/other workers.

## Interface

Frozen typed dataclasses: `RunSettings`, `AlleleSelectionSettings`,
`ClassificationSettings`, `ConsensusSettings`, `ConfidenceSettings`,
`CallingSettings`, `ReadPhasingSettings`, `ReferenceLayoutSettings`, and
`RuntimeSettings`. Central `DEFAULT_SETTINGS = RuntimeSettings()` and
`DEFAULT_LAYOUT = DEFAULT_SETTINGS.reference_layout` are immutable.

Runtime sections use the specification's snake_case names. `schema_version` is
integer 1. `repeat_dictionary` is an optional string resource path.
`load_settings(path: Path | None) -> RuntimeSettings` accepts optional sections
but requires an explicit schema_version=1 in a supplied JSON file. Passing None uses
central defaults. `settings_as_dict` produces JSON-compatible effective settings,
including list-encoded layout tuples and selected ladder range.

`ReferenceLayoutSettings` uses ordered `pre` and `after` tuples. It provides
`fixed_repeat_count`, `left_anchor_id` (pre[0]), `right_anchor_id` (after[-1]) and
`validate_repeats(repeats: Mapping[str,str])`. Dictionary membership is checked
when integration loads the selected dictionary; this module loads no resources
or tools. Category lists containing alternative repeats do not determine layout.

The coordinator's approved extension adds `reference_layout.min_units=1` and
`max_units=150`, with max_units >= min_units >=1; no separate LadderSettings.
Explicit existing library/CLI argument precedence is implemented by consumers.

## Validation behavior

The same constructor validation applies to direct library construction and JSON
loading. JSON rejects unknown or duplicated fields at every nesting level,
non-object sections, wrong schema/version types, bool-as-int, numeric strings,
nonfinite/overflowing numeric thresholds and invalid ranges. Layout arrays become
immutable tuples; IDs must be nonempty strings, unique within each side and
disjoint between sides. Direct constructors require tuples, preserving deep
immutability.

Positive minimums: threads1, min_coverage1, min_gap1, valley_min_points3,
valley_min_separation1, anchor_bases1, internal_downsampling1 and ladder minimum1.
Edit/search/tolerance/refinement/boundary-count limits are nonnegative. Unit
fraction is (0,1]; weights and penalties are [0,1]; QUAL high > low >=0.
Mapping-quality override accepts0; absent tool overrides remain None. VCF sample
names reject whitespace and nonprintable/control characters. An empty model path
preserves the existing tool default; other supplied resource/preset paths cannot
be empty. Relative reference/model/dictionary paths resolve against the config
file's directory after raw type validation; filesystem read errors propagate.

## Tests and checks

Tests were written before implementation. Initial collection exposed the absent
settings module. The follow-up stricter stage/range tests produced14 expected
failures before fixes (zero coverage/gap, invalid sample names, missing ladder
range and numeric overflow). Final focused validation on 2026-09-14:

```text
uv run --locked --all-extras pytest -q --no-cov tests/unit/test_runtime_settings.py
91 passed
uv run --locked --all-extras ruff check src/muc_one_span/settings.py tests/unit/test_runtime_settings.py
All checks passed
uv run --locked --all-extras ruff format --check src/muc_one_span/settings.py tests/unit/test_runtime_settings.py
2 files already formatted
uv run --locked --all-extras mypy src/muc_one_span/settings.py
Success: no issues found in 1 source file
```

Tests include all-section nondefault roundtrip, central default roundtrip,
immutability, partial sections, relative and absolute resource paths, dictionary
membership, strict types/ranges, unsafe VCF sample names and ladder range
serialization. No external bioinformatics tools, new dependencies, simulations,
final seeds or fresh validation inputs were used. End-to-end default scientific
parity and global CI remain the coordinator's integration gates.

Final module/test sizes:345/274 physical lines, each below650.

```text
72cfc6fe4148f7ceedd88cc7820d0719a17b132fd690238974c03ec17636d73a  src/muc_one_span/settings.py
1ade8c31378d29d841384b9cd707c8d68c9d978031ad48a23d0ac496b4a6c08e  tests/unit/test_runtime_settings.py
```
