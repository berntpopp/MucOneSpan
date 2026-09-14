# Milestone 2 N1 projection-unavailable semantics

Date: 2026-09-14. Scope: fresh Claude fix-review finding N1 only.

## Reproduction

Two focused regressions were added before implementation. Both failed:

- `projection_unavailable` with base dictionary-fit confidence 0.8 produced
  repeat and allele confidence 0.24, the same `0.3` multiplier used after a
  completed projection proves exact support absent; and
- the report rendered both `projection_unavailable` and `absent` events with an
  `Unsupported` badge.

The tests preserve `vcf_support=False`, `vcf_qual=0.0`, and the existing status.
They do not infer ALT selection, create support, or change an event call.

## Bounded behavior

VCF confidence weighting now distinguishes evidence states:

- `exact_sequence_concordance` keeps the existing QUAL-derived multiplier;
- only `absent`, after an available projection proves no exact match, keeps the
  existing `0.3` multiplier;
- `projection_unavailable` and `localization_ambiguous` apply no VCF evidence
  multiplier, so the repeat's heuristic dictionary-fit confidence is retained;
  and
- the independent boundary multiplier still applies afterward when the event is
  in the configured terminal repeat interval.

Consequently, for a non-boundary event with base confidence `b`, an unavailable
projection changes the result from `0.3b` to `b`. For a boundary event it changes
from `0.3b * boundary_penalty` to `b * boundary_penalty`. Allele confidence is
then recomputed as the unchanged arithmetic mean of all repeat confidences. The
confidence remains explicitly heuristic dictionary fit, not a probability and
not evidence that the VCF supports the mutation.

The HTML report now renders `projection_unavailable` as `Support unavailable`
and `localization_ambiguous` as `Support ambiguous`. An event with a completed
projection and no exact match remains `Unsupported`; an exact match remains
`Supported`. Historical false booleans without the new status continue to render
`Unsupported` for compatibility.

No sequence, repeat, mutation, genotype, phase, projection, support boolean,
support status, or QUAL value changes. No ALT/projection-selection heuristic was
added.

## Validation

- Red focused run: **2 failed, 16 passed**, on the confidence and badge
  distinctions above.
- Green focused run:
  `uv run --locked --all-extras pytest tests/unit/test_classification_evidence.py
  tests/unit/test_report.py -q --no-cov` — **18 passed**.
- Follow-up clarification reproduced `localization_ambiguous` separately: the
  two focused assertions initially failed because it still received the `0.3`
  multiplier and `Unsupported` badge, before the final status partition.
- The first broader classification/support/report run exposed one pre-N1
  expectation in
  `tests/unit/test_classify.py::TestContinuousQualScoring::test_no_vcf_support_gives_low_confidence`,
  which supplied no projection context but expected the `absent` penalty. The
  coordinator replaced it with an actual replay-verified `absent` fixture rather
  than weakening the assertion. Final broader run: **78 passed**.

Exact new regressions:

- `tests/unit/test_classification_evidence.py::test_projection_unavailable_preserves_fit_confidence_while_absent_penalizes`
- `tests/unit/test_report.py::test_report_distinguishes_unavailable_from_absent_support`

Focused Ruff and formatting checks passed for the three Python files. Configured
mypy reported no issues in `classify.py`; report rendering validates the Jinja
template; `git diff --check` passed. The four authored code/template/test files
are 604, 108, 333 and 200 lines. Full repository gates remain coordinator work.
