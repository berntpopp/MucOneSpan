# Independent integration audit — 2026-09-14

This is an independent Codex audit of the current tracked diff and every new
production Python file visible in the worktree. It is not a Claude review. The
audit was read-only except for this report and did not run a pipeline, caller,
simulator, or other expensive command. One small in-process Python reproduction
was used for evaluator arithmetic; the remaining reproductions follow directly
from the named call paths.

## Actionable defects

### HIGH — duplicate observations still earn duplicate primary accuracy credit

`evaluate_sample` computes `sequence_exact`, `structure_exact`, and count metrics
for every assigned prediction before applying the `independent` evidence gate
(`src/muc_one_span/evaluation/scoring.py:74-143`). The gate only affects the
`all_*_exact` sample flags. `aggregate` then exposes the ungated
`sequence_exact` numerator as the primary `sequence_accuracy` endpoint
(`scoring.py:250`). Consequently, two copies of one observation can earn two
recovered haplotypes even though the same row correctly says that independent
haplotype evidence is false.

Reproduction, executed against the current tree:

```python
truth = TruthSample("s", (
    TruthHaplotype("h1", "A", ("X",)),
    TruthHaplotype("h2", "A", ("X",)),
))
predictions = (
    PredictedAllele("p1", "A", ("X",), 1, -8),
    PredictedAllele("p2", "A", ("X",), 1, -8),
)
row = evaluate_sample(truth, RunObservation("completed", predictions))
assert row["independent_haplotype_evidence"] is False
assert row["metrics"]["sequence_exact"] == {"min": 2, "max": 2}
assert aggregate([row])["totals"]["sequence_accuracy"]["value"] == 1.0
assert row["metrics"]["all_sequences_exact"] == {"min": 0, "max": 0}
```

This conflicts with the governing specification's rule that one observed
sequence cannot be counted as two recovered haplotypes without evidence for two
selections. It can inflate individual sequence, structure, count, and
conditional-sequence endpoints for identical truth haplotypes. Cap evidence-
dependent matched credit per observation/source group, or add evidence-aware
primary numerators and ensure acceptance/reporting uses those. Add a regression
that checks aggregate individual accuracy, not only `all_sequences_exact`.

### HIGH — resolved reconstruction call rate accepts extra alleles

`aggregate` defines a resolved reconstruction as completed, no missing alleles,
and independent evidence (`scoring.py:285-294`). It does not require zero extra
alleles or equal truth/prediction cardinality. Three independent predictions for
diploid truth therefore produce a 100% resolved reconstruction call rate even
though the sample has an explicit extra allele and cannot be a resolved diploid
reconstruction.

Reproduction, executed against the current tree:

```python
truth = TruthSample("s", (
    TruthHaplotype("h1", "A", ("X",)),
    TruthHaplotype("h2", "C", ("X",)),
))
predictions = tuple(
    PredictedAllele(name, seq, ("X",), 1, -8, sequence_source=name)
    for name, seq in (("p1", "A"), ("p2", "C"), ("p3", "G"))
)
row = evaluate_sample(truth, RunObservation("completed", predictions))
assert (row["missing_alleles"], row["extra_alleles"]) == (0, 1)
assert aggregate([row])["totals"]["resolved_reconstruction_call_rate"] == {
    "numerator": 1, "denominator": 1, "value": 1.0,
}
```

Require `extra_alleles == 0` (equivalently exact cardinality) in this rate and
add a three-prediction/diploid regression. Candidate-output `call_rate` can
remain broader because that distinction is documented.

### HIGH — Click validation can leave a stale completed run status scoreable

`@record_run_status` wraps only the `run` callback (`src/muc_one_span/cli.py:350-412`).
Click validates `--input` and an explicit `--reference` with
`click.Path(exists=True)` before invoking that callback. A rerun against an
existing output directory with a missing input therefore exits nonzero without
replacing the old `run_status.json`. If the old sidecar says `completed`, then
`load_observation` accepts the stale summary when evaluation has no external run
record (`src/muc_one_span/evaluation/artifacts.py:124-176`).

Minimal regression:

1. Put a valid `summary.json`, consensus FASTA, and completed schema-1 sidecar in
   an output directory.
2. Invoke `CliRunner().invoke(main, ["run", "--input", missing_path,
   "--output-dir", output_dir])` and assert a nonzero exit.
3. The current sidecar remains `completed`; `load_observation(output_dir)` can
   return completed predictions from the earlier run.

Move existence validation for run inputs inside the status-wrapped callback, or
add a command-level status hook that can resolve `output_dir` before parameter
validation. Preserve the benchmark driver's nonzero run record as an additional
guard, but do not make correctness depend on always running through that driver.

### MEDIUM — the strict “supported” endpoint accepts legacy proximity booleans

`Event.supported` checks only `frameshift`, `template_match`, and the historical
`vcf_support` boolean (`src/muc_one_span/evaluation/models.py:21-24`). The artifact
adapter defaults absent `vcf_support_status` to `unknown` but still passes a true
boolean (`evaluation/artifacts.py:84-104`). Thus historical position-overlap calls
without a projection status are counted in `supported_event_*` and supported
sample metrics as though they had exact sequence concordance.

The one-line reproduction `Event(1, "X", "dupC", True, True, True,
"unknown").supported` returns true. Historical cached summaries contain exactly
this combination: `vcf_support: true` and no status. The governing spec states
that supported means exact template plus exact VCF concordance and that legacy
proximity support must be reported separately. `.planning/task1-evidence.md`
explicitly documents the looser implementation, so this is a documented spec
divergence rather than a hidden caller regression.

Require `support_status == "exact_sequence_concordance"` for the strict endpoint,
and expose a separately named legacy-boolean/proximity endpoint if the baseline
still needs it. Do not compare the old boolean-supported numerator directly with
new projection-supported results.

### MEDIUM — read-backed phase can be accepted with no assigned read evidence

`_assigned_read_count` validates names but permits an empty read list and ignores
the read-list phase/haplotype fields (`src/muc_one_span/read_phasing.py:74-89`).
`phase_same_length_reads` records that count and then accepts any candidate VCF
whose variants pass the common-PS policy (`read_phasing.py:137-151`). A mocked or
future Whatshap result with zero assigned records and phased GT/PS records is
therefore promoted as `status="phased"`, despite the method being described as
read-backed.

Add a unit regression in which the candidate variants are phased but
`assigned_reads.tsv` contains only its header; the current function selects and
compresses the candidate. Require at least one valid assigned record that belongs
to the accepted phase set, and validate the phase/haplotype columns rather than
counting every syntactically valid name. Whether both haplotype groups are
required should follow a stated scientific rule; one linking read can sometimes
determine the complementary phase from a diploid genotype.

Impact is currently bounded because the CLI leaves `read_phase=False`. The helper
is experimental and default-disabled after the recorded extra-event regression,
but its acceptance gate should be sound before any future promotion.

### MEDIUM — mixed-platform inventory can silently reuse one model for both platforms

The benchmark runner infers each sample's platform from simulator metadata, but
resolves `model` independently from one command/environment value
(`src/muc_one_span/benchmarking.py:189-232`). With a mixed HiFi/ONT root and a
HiFi `CLAIR3_MODEL`, the records correctly switch `platform` to `ont` while still
passing the same HiFi model path. Per-entry models work only when supplied in an
explicit inventory and can conflict with an exported global model.

The benchmarking guide mitigates this procedurally by requiring separate
platform directories, but the new discovery interface itself accepts the unsafe
combination. This falls short of the runner requirement to avoid silent HiFi/ONT
mismatch. Reject a discovered mixed-platform inventory when only one model is
available, require explicit per-sample models for mixed runs, or validate model
metadata against the selected platform. Add a two-sample metadata regression
that asserts no caller invocation when one global model would cross platforms.

## Explicitly documented or unmet methods, not introduced defects

- Read-backed phasing is implemented only as an experimental library option and
  is disabled in the CLI/default pipeline because cached validation added a
  supported false event. This is a recorded failed promotion gate, not a missing
  default call accidentally overlooked. Documentation currently conflicts:
  `docs/reference/limitations.md` implies users can enable it, while
  `docs/development.md` says it is not implemented and no CLI flag exists. State
  “library-only, experimental, disabled by default” consistently.
- IUPAC-mode replay deliberately returns `ambiguous_genotype_selection` for
  heterozygous indels even though the recorded bcftools 1.17 experiment shows
  `-H I` selecting ALT for the tested 0/1 indels. This makes support unavailable
  for those alleles, but the limitation and its conservative intent are explicit;
  it should not be described as a consensus byte mismatch.
- Whole-sequence reference confidence, independent molecule counts, robust
  length rescue, the ladder left-flank correction, read-phasing promotion, and
  universal model hashes remain explicitly unvalidated or unavailable. The new
  diagnostic read-span module does not claim to alter calling.

## Scope checked

Reviewed the complete current code in both legacy scripts, `scripts/evaluate.py`,
`scripts/benchmark_classification.py`, all new evaluation modules,
`benchmarking.py`, `classification_summary.py`, `phasing.py`, `read_evidence.py`,
`read_phasing.py`, `run_status.py`, and `variant_support.py`; traced their call
sites through changed allele calling, consensus, VCF parsing, classification,
CLI persistence, report fields, tests, specification, evidence, and review
dispositions. No full test target or external-tool job was run for this bounded
audit.
