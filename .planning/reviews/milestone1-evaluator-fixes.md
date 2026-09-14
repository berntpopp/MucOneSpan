# Milestone 1 evaluator fixes

Scope: Claude milestone-1 findings 1 and 4 only. Owned implementation files:
`src/muc_one_span/evaluation/artifacts.py` and
`tests/unit/test_evaluation_artifacts.py`. No caller or scoring policy change.

## Findings reproduced before implementation

Focused artifact tests first returned **6 failed, 24 passed**. A real-shape
summary with one observed allele and one explicitly unresolved duplicate alias
was rejected by strict equality of allele/classification keys. Separately,
`running` was rejected during sidecar validation before timeout or nonzero-exit
evidence could be considered. A stale running status without external failure
also lacked an explicit incomplete-run explanation.

## Final behavior

A missing classification/FASTA is accepted only for an allele with both
`reconstruction_status: not_separately_resolved` and a string
`candidate_duplicate_of` pointing to an actually classified allele. A bare alias,
bare unresolved status, missing target, self-reference, and malformed target
remain invalid. Classification extras and arbitrary missing classifications
remain invalid; `repeats.json` must still match actual classifications exactly.
Accepted unresolved aliases receive an explicit warning and make the sample
`ambiguous_reconstruction`. They never create fabricated predictions. The test
scores two truth haplotypes against one retained observed prediction and verifies
one missing allele, one exact sequence, no all-sequences-exact result, and no
normal true negative.

`running` is now a recognized sidecar state. Nonzero process exit, timeout, or
explicit execution failure returns `execution_failed` before opening any stale
summary, preserving the sidecar and external error. Without external failure,
running returns zero predictions and `invalid_artifacts` with the explicit error
`run incomplete: status is running`. This preserves the existing evaluator
nonzero-exit contract without introducing a new status that the evaluator CLI
would accidentally treat as successful. Typed `insufficient_evidence` retains
its previous precedence over the CLI's compatible nonzero no-call exit.

## Validation

- Artifact-only tests after fixes: **30 passed**.
- `uv run --locked --all-extras pytest tests/unit/test_evaluation_artifacts.py
  tests/unit/test_evaluation_scoring.py tests/unit/test_evaluation_cli.py
  --no-cov -q`: **46 passed**.
- Focused Ruff check and format check: passed, both files formatted.
- `uv run --locked --all-extras mypy
  src/muc_one_span/evaluation/artifacts.py`: success, no issues.

The first mypy invocation identified overly narrow inference for the warnings
tuple after adding alias warnings; an explicit `tuple[str, ...]` annotation fixes
that without suppression. Final focused tests and checks were rerun afterward.
No full exposed-panel evaluation or fresh final generation was performed in
this bounded task; those remain coordinator workstreams. No commit was made.
