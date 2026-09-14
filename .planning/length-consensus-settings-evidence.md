# Length, consensus, and reference settings integration

Date: 2026-09-14. Worktree: `production-validation`. Scope: `alleles.py`,
`consensus.py`, `ladder.py`, their three owned settings/evidence test files.
The shared settings schema and CLI integration are owned separately.

## Implemented contract

- `detect_alleles(..., *, settings=None, reference_layout=None)` takes
  `AlleleSelectionSettings` and `ReferenceLayoutSettings`. Cluster gap, minimum
  valley points, minimum valley separation, and accepted AS refinement shift
  all resolve from the supplied selection settings. Explicit `min_coverage`
  remains authoritative. The private splitter and allele builder accept the
  corresponding keyword settings so direct library callers receive the same behavior.
- Every built allele records `fixed_repeat_count` from the selected reference
  layout. Both candidate total length and selected reference length use that
  value. Aliased candidates retain it. `PRE_AFTER_REPEAT_COUNT` remains an
  import-compatible alias derived from the default layout, whose count is nine.
- `build_contig` and `generate_ladder_fasta` accept optional keyword
  `settings: ConsensusSettings` and `reference_layout: ReferenceLayoutSettings`.
  `flank_length=None` resolves the configured flank length. Explicit flank length
  wins. The generator's `min_units=None` and `max_units=None` resolve the
  reference layout's configured range; explicit range values win.
- Ladder construction uses the ordered IDs in `reference_layout.pre` and
  `.after`, validating their presence before creating a FASTA output. It does
  not count repeat dictionary categories, which contain alternate repeat IDs.
- `trim_flanking(..., flank_length, output_path, ..., settings=None,
  reference_layout=None)` retains its existing positional slots; callers pass
  `None` in the flank slot to use settings. `build_consensus_per_allele` has an
  optional flank argument defaulting to None and forwards both settings.
  Boundary anchors use the layout's first pre-repeat and last after-repeat;
  component length and search tolerance come from consensus settings.
  Supplied repeat dictionaries are validated even when flank length is zero.

No default scientific thresholds or phase acceptance behavior changed. The
existing left reference flank slicing remains as before; this change does not
repair the separately documented distal flank/reference discrepancy.

## Red to green evidence

Before source edits, the three owned test files produced **11 failed, 7 passed**.
Failures exposed missing settings keywords, absent selected-ID validation, and
missing configurable ladder range. The initial range test referred to a proposed
separate class before coordination settled on the existing reference-layout
section; it was updated to the agreed API before implementation.

After integration, **77 passed, zero skipped** across the existing allele,
ladder, consensus tests and the new settings tests:

```sh
uv run --locked --all-extras pytest -q \
  tests/unit/test_alleles.py tests/unit/test_allele_evidence.py \
  tests/unit/test_ladder.py tests/unit/test_reference_settings.py \
  tests/unit/test_consensus.py tests/unit/test_consensus_settings.py --no-cov
uv run --locked --all-extras ruff check \
  src/muc_one_span/alleles.py src/muc_one_span/ladder.py src/muc_one_span/consensus.py \
  tests/unit/test_reference_settings.py tests/unit/test_consensus_settings.py \
  tests/unit/test_allele_evidence.py
uv run --locked --all-extras mypy \
  src/muc_one_span/alleles.py src/muc_one_span/ladder.py src/muc_one_span/consensus.py
```

Ruff passed; mypy reported no issues in the three source files. Modified Python
files remain below 650 physical lines (alleles 596; ladder 120; consensus 297).
Tests exercise nondefault gaps, valley support/separation, refinement tolerance,
three fixed repeats instead of nine, ordered alternate IDs, missing-ID failure,
configured and explicitly overridden flanks/ranges, two-base anchors with zero
versus two-base tolerance, and both settings/layout propagation through the
per-allele consensus function. External tools are mocked in these unit tests.

## Saved-source default parity

A separate `uv run --locked --all-extras python` check imported each saved source
by absolute file spec from
`tests/results/production_validation_20260914/configuration/before/src/muc_one_span`
and compared it with the configured implementation. Outputs existed only in a
`TemporaryDirectory`, using the bundled repeat dictionary. No fresh reads,
final simulation, calling, or accepted parameter tuning occurred.

Snapshot source SHA-256 values:

| Module | SHA-256 |
|---|---|
| alleles.py | `e26a9d9f62a648fd176981084c26daa6b095649c8aebad1062dd92a90b2c34b7` |
| consensus.py | `d6cca7ca2a4ec7070069b1d0c2b2448947f9790db8c05b8b33d93f6fc6022e0d` |
| ladder.py | `62ba8b4e8f52368086db72e56f6e61f4f87fce390468c4030ea37bb097ad7e6c` |

Parity results:

- 15 contig dictionaries identical: canonical counts 1, 25, 51, 100, 150 crossed
  with flank lengths 0, 100, 500.
- Three entire FASTAs byte-identical: all defaults; range 1–3/flank 0; range
  20–25/flank 100/line width 60.
- Nine trimmed FASTA files and their coordinate/method contexts identical:
  empty sequence with flank 0; four bases with flank 500; AACCCGG with flank 2;
  counts 1/51 crossed with flank 0/100/500 and the bundled dictionary.
- Six allele dictionaries identical after removing only the newly added
  `fixed_repeat_count=9` from each candidate: distributions `{51:30}`,
  `{51:30,80:20}`, `{40:20,44:20}`, `{40:20,45:20}`,
  `{50:9,51:10,80:11,110:10}`, and `{51:30,52:20}`.
- Two mocked-SAM valley outputs identical: points 40/41/43 and 40/41/42, with
  ten insertion bases at 41 and zero at both endpoints.

Full default ladder SHA-256:
`2349a88c8435790aff3deb1c446901888ba9f862d7d4e9f00ed59bd8e8272e34`.
The other two FASTAs have SHA-256
`3e1f0f9afe9d26eeb8d626be1a1bb09c251e6856162496e45fdab642bca6cad6` and
`08cb60d27a3fa171045ad679399dfc44119f29c4c4fa3aa58a2423500935ba5c`.

These checks establish default parity for the exercised pure-Python cases.
They do not establish a biological improvement or substitute for full real-tool
integration and final paired validation. The coordinator owns repository-wide
CI, CLI/evaluator wiring, and freeze preparation.
