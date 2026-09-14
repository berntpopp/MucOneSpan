# Task 2 exact-scoring evidence — 2026-09-14

## Scope and fixed contract

- Baseline: `d8390b3c244ef8f3240af74b92db12b50dfc77d1`.
- Owned implementation: `src/muc_one_span/repeat_alignment.py`,
  `tests/unit/test_repeat_alignment.py`, and
  `scripts/benchmark_classification.py`.
- Public `edit_distance(str, str) -> int` remains available through the existing
  `muc_one_span.classify` import. Equality is literal. The full-matrix traceback
  and its substitution/insertion/deletion tie order are unchanged. Only the
  inaccurate insertion-coordinate comment was corrected.
- No `classify.py` edit belongs to this task. The benchmark replaces
  `muc_one_span.classify.edit_distance` in one loaded process and restores it in
  `finally`, so scalar and bit-vector modes use the same classifier semantics.

## Test-first baseline

Before changing production code, the independent full-matrix scalar oracle was
run on 10,000 pairs in both argument orders. The seed is `20260914`; cases cover
empty and single-symbol strings, repeated symbols, substitutions and indels,
case differences, `N`, IUPAC symbols, and lengths 31/32, 63/64/65, and
127/128/129.

```text
uv run --locked --all-extras pytest \
  tests/unit/test_repeat_alignment.py::test_public_edit_distance_matches_scalar_oracle_on_10k_seeded_pairs \
  tests/unit/test_repeat_alignment.py::test_edit_distance_uses_literal_symbol_equality \
  tests/unit/test_repeat_alignment.py::test_traceback_preserves_first_hit_ties_and_insertion_position \
  --no-cov -q
3 passed in 4.26s
```

The optimization-specific test was then observed failing for the intended
reason:

```text
uv run --locked --all-extras pytest \
  tests/unit/test_repeat_alignment.py::test_bitvector_entry_point_handles_machine_word_boundaries \
  --no-cov -q
FAILED: AttributeError: repeat_alignment has no attribute '_edit_distance_bitvector'
```

The benchmark-harness restoration and cached-panel discovery tests were also
observed failing before their respective script functions existed. The
review-requested near-exact control test failed on the missing
`_near_exact_control` helper before implementation. The later >5 kb symmetry
case passed immediately against the implemented bit-vector recurrence; it
extends adversarial coverage rather than serving as the original red target.

The unmodified scalar classifier baseline, measured before implementation in
the same environment over three calls/batches, was:

```text
difficult_bench5003 median_seconds_per_call 18.06343891800134
exact_X_120 median_seconds_per_call 0.0005805907900139573
```

## Implementation and exactness

`edit_distance` now calls a dependency-free Myers bit-vector recurrence using
Python integers, with explicit identical and empty-string handling. The shorter
argument is selected as the bit pattern, which preserves symmetry and reduces
integer width. A direct 6,000/6,001-base test checks whole-VNTR-scale behavior
and argument symmetry.

The final focused test command after all source changes was:

```text
uv run --locked --all-extras pytest tests/unit/test_repeat_alignment.py \
  tests/unit/test_classify.py --no-cov -q
57 passed in 2.94s
```

After the benchmark-only typing adjustment, the Task 2 test file was rerun:

```text
uv run --locked --all-extras pytest tests/unit/test_repeat_alignment.py --no-cov -q
8 passed in 2.25s
```

## Controlled benchmark

Command:

```text
uv run --locked --all-extras python scripts/benchmark_classification.py
```

Artifact:
`results/production_validation_20260914/scoring/benchmark.json` (ignored through
the existing `results` symlink).

Environment: CPython 3.12.9, Linux 7.0.0-31 x86_64, three paired repetitions.
Each repetition alternates scalar/bit-vector order. The difficult case uses one
classification per observation, the exact control batches 1,000 classifications,
and the scorer-sensitive near-exact control batches 100. Every row compares the
complete classification dictionary.

| Case | Scalar median | Bit-vector median | Speedup | Complete equality |
| --- | ---: | ---: | ---: | --- |
| difficult bench5003 | 17.612304 s | 2.214009 s | 7.955x | yes |
| exact X × 120 | 0.827535 ms | 0.828880 ms | 0.998x | yes |
| near-exact X × 120 | 12.034696 ms | 2.394635 ms | 5.026x | yes |

All 141 available trimmed consensus FASTAs under the cached development root
(HiFi, ONT, former held-out samples, and perturbation outputs) produced equal
complete dictionaries. The harness exited zero with all gates true: zero
dictionary differences, difficult-case speedup at least 2x, and no greater than
20% regression in either batched control.

The artifact records hashes for all 141 FASTAs, the repeat dictionary, loaded
classifier source, scorer source, and benchmark harness. Key hashes from the run:

- difficult FASTA:
  `69d3482c6cd2967705ce136dd5e32f43f867ea6a1becef4dc94e8f790597ae29`
- repeat dictionary:
  `ba901bfdc969bc0214170d0de91697f032963bae6de096252356b9de44f1abb3`
- production scorer:
  `c4a8b89edc304d9234fbf597651a5e74ae64e92b1bebec4fd280b8702e11f36a`
- loaded classifier source:
  `ccc9bca7dac4720b72366f11ebecdbba957708cb0991a07be9a77cef46c3d0c8`
- benchmark harness at timing:
  `3904bb99538064a31b4135d9bf779c5fbb764a9b82459f2c4a2a9cbcc67ef900`

The embedded scalar copy is verified before timing against the normalized AST
hash `2c2eb12bd7c52c46fc62b7e19a20b4bf1db8f7edf4760fc3ec591307b1ce5ac2`
from the baseline module, whose complete file hash is
`4e9fcdfdea4e2ed6a0c0ab7ed69c7e434a29fa050c6f4d36cd51f6575b222ccf`.

After timing, focused mypy required expressing the intentional process-local
module patch through an `Any` cast. This changes no executed operation: both
versions assign the scorer, classify, and restore in the same `try/finally`.
The final harness hash is
`1936231893be02b723cd3d31332a03419c83207e30c4d24b61d7e30cfeed907f`.
The benchmark was not repeated because the coordinator had begun the heavy
53-sample plus perturbation pipeline job; overlapping it would invalidate
controlled timing. The run hash remains unchanged in the artifact to preserve
the exact source provenance of the measured run.

## Focused quality checks

```text
uv run --locked --all-extras ruff check src/muc_one_span/repeat_alignment.py \
  tests/unit/test_repeat_alignment.py scripts/benchmark_classification.py
All checks passed!

uv run --locked --all-extras ruff format --check \
  src/muc_one_span/repeat_alignment.py tests/unit/test_repeat_alignment.py \
  scripts/benchmark_classification.py
3 files already formatted

uv run --locked --all-extras mypy src/muc_one_span/repeat_alignment.py \
  scripts/benchmark_classification.py
Success: no issues found in 2 source files

uv run --locked --all-extras python scripts/benchmark_classification.py --help
exit 0; all benchmark options displayed
```

File lengths are 158, 183, and 377 lines respectively, below the 649-line gate.
`git diff --check` is clean for the owned files. Combined `make ci-check` and
pipeline validation are intentionally left to the coordinator while other
exclusive workstream files are still changing.
