# Task 1: Strict offline evaluation evidence

Date: 2026-09-14. Scope: new `src/muc_one_span/evaluation/`,
`scripts/evaluate.py`, and `tests/unit/test_evaluation*.py`.
No caller imports simulator truth. No dependencies or existing running scripts
changed by this worker. Global quality/CI/docs/build checks are coordinated after
integration; the checks below apply to this task's source only.

## Implementation and review dispositions

- Strict actual MucOneUp FASTA/structure/statistics/mutated-unit validation.
  Every complete truth sequence is reconstructed from dictionary units and actual
  mutated units with exact flanks. SNP-bearing truth is explicitly rejected as
  unsupported; mutation coordinates/counts and target markers must agree.
- Historical stale per-unit lengths and missing metadata seeds produce named
  warnings, with actual lengths authoritative. TSV read commands provide nominal
  molecule requests and read seeds separately from misleading Coverage values.
  Retained records/model hashes remain null unless provided by external inventory;
  no molecule counts are fabricated from alignment counts or duplicate QNAMEs.
- Explicit nonempty inventory is mandatory. Objects can provide truth_sample,
  truth_dir, result_dir, input, platform, run_record and additional provenance.
  Every sample receives a row; malformed truth/artifacts remain visible.
- One-to-one maximum-cardinality assignment minimizes complete-sequence global
  distance. Every optimum is retained. Event annotation identity is exact
  1-based repeat index, parent and mutation name, with multiset matching.
  Missing truth events become FN and every unmatched prediction event becomes FP.
- Reports expose conservative min/max endpoints for ties. Direct sequence string
  equality determines sequence exactness independently of the distance scorer.
  Identical consensus copies require distinct source IDs and explicit independent
  genotype evidence before all-sample sequence/structure/count success. Reused
  source IDs cannot establish two haplotypes. Distinct biological sequence matches
  remain individually measurable even if genotype evidence is unavailable.
- `sequence_source`, `independent_haplotype_evidence`, `phase_status`,
  `genotype_status`, and `sequence_identity_status` are accepted from allele
  metadata, falling back to classification metadata. IUPAC consensuses and
  unresolved phase states are `ambiguous_reconstruction`.
- `run_status.json` version 1 supports completed, insufficient_evidence and
  execution_failed. Typed insufficient/failed states score no stale predictions.
  A completed sidecar cannot override a recorded nonzero execution exit.
  Historical missing execution provenance is explicitly warned, never invented.
- All-candidate and supported-event metrics are separate. Supported means all
  three booleans frameshift, template_match and vcf_support are true. Support
  status counts are retained where supplied. Normal no-calls and ambiguous
  reconstructions cannot become true negatives. Rates carry numerator and
  denominator, with null for undefined ratios. Zero biological accuracy is a
  valid report; input/execution failures produce nonzero CLI status.
- Event nucleotide normalization/equivalence is explicitly `not_assessable`.
  The supported rule includes legacy caller boolean evidence; its presence does
  not retrospectively repair old caller projection defects. Algorithm success
  gates must use conservative min endpoints, not max bounds.

## TDD and scoped validation

Observed initial red: three ModuleNotFoundError collection failures before the
package existed. CLI red: three missing-script failures before its implementation.
Subsequent observed behavioral reds covered typed statuses, duplicate consensus
and source evidence, separate supported/no-call rates, metadata molecule semantics,
and malformed evidence metadata. Logs retained below.

Final command: `uv run --locked --all-extras pytest tests/unit/test_evaluation*.py --no-cov -q`:
**50 passed**. No external tools or generated datasets are prerequisites.

Scoped checks passed:

```bash
uv run --locked --all-extras ruff check src/muc_one_span/evaluation scripts/evaluate.py tests/unit/test_evaluation*.py
uv run --locked --all-extras ruff format --check src/muc_one_span/evaluation scripts/evaluate.py tests/unit/test_evaluation*.py
uv run --locked --all-extras mypy src/muc_one_span/evaluation scripts/evaluate.py
uv run --locked --all-extras python scripts/evaluate.py --help
git diff --check -- src/muc_one_span/evaluation scripts/evaluate.py tests/unit/test_evaluation*.py
```

Mypy reports no issues in seven source files. Every authored file is below 649
physical lines. No commit or global formatter was run by this worker.

## Full exposed development evaluation

Artifacts and inventories are retained in ignored
`tests/results/production_validation_20260914/evaluation/`:
`hifi.json`, `ont.json`, `heldout.json`, `perturbations.json`, matching
`*-inventory.json`, CLI `*.log`, `tests.log`, TDD red logs and
`source-sha256.json`. Source hashes include new untracked source/test files.

Inventories were built from all truth directories containing simulation stats,
not successful result directories. Original truth root has 47 samples: 44 HiFi
and three ONT. Exposed-six and perturbation roots contain 6 and 24 respectively.
Each input path was explicitly selected using existing measurement metadata and
its current bytes were hashed. The exposed-six set remains development despite
its historical directory name `heldout`.

Exact commands (the explicit truth root for HiFi/ONT is tests/data/generated):

```bash
uv run --locked --all-extras python scripts/evaluate.py tests/results/deep_validation_20260914/hifi --truth-root tests/data/generated --expected-samples tests/results/production_validation_20260914/evaluation/hifi-inventory.json --output tests/results/production_validation_20260914/evaluation/hifi.json
uv run --locked --all-extras python scripts/evaluate.py tests/results/deep_validation_20260914/ont --truth-root tests/data/generated --expected-samples tests/results/production_validation_20260914/evaluation/ont-inventory.json --output tests/results/production_validation_20260914/evaluation/ont.json
uv run --locked --all-extras python scripts/evaluate.py tests/results/deep_validation_20260914/heldout --truth-root tests/results/deep_validation_20260914/heldout_data --expected-samples tests/results/production_validation_20260914/evaluation/heldout-inventory.json --output tests/results/production_validation_20260914/evaluation/heldout.json
uv run --locked --all-extras python scripts/evaluate.py tests/results/deep_validation_20260914/perturbations --truth-root tests/results/deep_validation_20260914/perturbation_data --expected-samples tests/results/production_validation_20260914/evaluation/perturbations-inventory.json --output tests/results/production_validation_20260914/evaluation/perturbations.json
```

Exit codes are 0, 0, 0, 1 respectively. All **77/77 truth samples** validate.
The nonzero perturbation exit is intentional: six recorded caller execution
failures and one incomplete artifact sample remain in every truth denominator.

| Development stratum | Samples | Exact sequences | Exact count pairs | Candidate event TP/FN/FP | Supported event TP/FN/FP |
| --- | ---: | ---: | ---: | --- | --- |
| Original HiFi | 44 | 39/88 | 40/44 | 24/4/14 | 24/4/0 |
| ONT | 3 | 3/6 | 0/3 | 2/0/0 | 2/0/0 |
| Exposed six | 6 | 1/12 | 2/6 | 2/1/0 | 2/1/0 |
| Existing perturbations | 24 | 6/48 | 9/24 | 6/8/12 | 6/8/1 |

Status counts: HiFi 12 completed +32 ambiguous; ONT three ambiguous; exposed six
all ambiguous; perturbations two completed +15 ambiguous +six execution_failed
+one invalid_artifacts. Ambiguity includes literal IUPAC bases and preserves
candidate sequence/event assessment while withholding normal true-negative claims.

The newly identified incomplete perturbation is
`sample_dupc_60_80__n20_seed1702`: classification keys do not cover the declared
allele keys. The prototype credited its one available allele (and mutation),
producing sequence7/eventTP7. The strict evaluator refuses successful current
predictions for incomplete artifacts, yielding sequence6/eventTP6 while retaining
48 truth haplotypes and 14 expected mutation events. This is a denominator/failure
contract correction, not a claimed biological regression or improvement.

Limitations: old results do not carry phase/source evidence or modern exact VCF
projection statuses; retained independent molecules and model hashes are not
universally available. No confidence intervals are computed by this reusable
scorer; coordinator reporting must keep sample/perturbation dependence explicit.
