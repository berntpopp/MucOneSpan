# Milestone2 F5 — fail-closed evidence status derivation

Date: 2026-09-14. Source finding: actual Claude Fable milestone2 response,
`.planning/reviews/milestone2-response.json`, finding F5 only.

## Reproduced defect

The artifact adapter derived `completed` by excluding a fixed list of known
unresolved states. Explicit `unknown`, a new future status, an empty string or
null could therefore pass as completed and become eligible for normal
true-negative credit. Unknown reconstruction status in an allele could also be
ignored when the classification supplied a different status.

Eighteen regression cases failed before implementation: four unknown-value forms
across both phase/reconstruction fields and both allele/classification locations,
plus two cases where a valid status at one location masked an unknown status at
the other. Existing valid producer values and absent legacy metadata remained
passing controls. The red suite reported18 failures and49 passes.

## Repair

Replace the deny-list with an explicit allow-list applied to supplied fields in
both raw metadata locations. This is necessary because the prediction model's
`unknown` defaults cannot distinguish omitted legacy fields from explicitly
unknown evidence.

Known producer phase statuses accepted for completed artifact categorization:

- `phased`
- `single_heterozygous_unordered`
- `no_informative_heterozygosity`

Known producer reconstruction statuses accepted:

- `complete_segmentation` (classifier)
- `candidate_reference_confidence_unverified` (consensus candidate metadata)

Every other explicitly supplied phase/reconstruction status yields
`ambiguous_reconstruction`. The existing malformed-type validation remains in
place. An unknown value in either metadata layer cannot be hidden by a recognized
value in the other. Predictions remain available for literal sequence/event
scoring and truth denominators; they cannot become normal true negatives through
the completed-status predicate.

Legacy artifacts with absent additive fields preserve their historical loader
behavior and model defaults. This compatibility exception applies to absence,
not an explicit `unknown`, empty or null value. Existing unresolved aliases,
ambiguous-base checks and execution-failure precedence are unchanged.

Accepting `complete_segmentation` still describes segmentation, not empirical
reference confidence or complete biological reconstruction; this repair changes
unknown-state handling rather than redefining the accepted scientific statuses.

## Verification

- Artifact suite:67 passed after repair.
- Artifact + evaluation CLI + scoring regression command:
  `uv run --locked --all-extras pytest tests/unit/test_evaluation_artifacts.py tests/unit/test_evaluation_cli.py tests/unit/test_evaluation_scoring.py --no-cov -q`
  —88 passed in0.08s after final formatting.
- Ruff check and format check: pass.
- Configured mypy for artifacts.py: pass.
- Physical lines: artifacts.py242; artifact tests331, both below649.

Only artifacts.py and its unit test file were edited for F5. No models, scoring,
other review findings, heavy datasets or final validation seeds were touched.

Source SHA256:
`8c68d62a28f0d0ac8f820d5a9fcd657b82c73372bf630d2b77648cf361fbd411`.
Test SHA256:
`853085ed204ef89b52fa116a432b866732614c89df8c99d325bc56748594ed43`.
