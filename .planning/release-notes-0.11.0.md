# MucOneSpan v0.11.0

MucOneSpan 0.11.0 adds configurable, reproducible HiFi/ONT experiments and improves exact scoring, variant evidence, failure reporting and scientific evaluation.

- **Exact classification optimization:** 7.955× faster on the paired difficult classification benchmark, with identical complete cached outputs. This is a stage measurement; whole-pipeline speedup has not been established.
- **Clearer scientific evidence:** signed net-indel frameshifts, variant-specific projected VCF concordance, explicit unresolved sequence and haplotype ambiguity, and visible external-tool failures.
- **Strict evaluation:** one-to-one allele matching, exact mutation identity and position, missing/extra alleles, assignment ambiguity, failed samples and no-calls remain visible.
- **Configuration-driven runs:** strict JSON settings, existing CLI overrides, input/reference/dictionary hashes and explicit consensus sample/haplotype behavior.
- **New MucOneUp experiment designs:** editable HiFi/ONT JSON examples, explicit seeds and model configurations, usable-read counts, truth checks and provenance manifests. The new development smoke generated 10 HiFi and 12 ONT usable reads from 12 requested templates per platform.
- **Packaging protection:** wheel and source distribution validation excludes local generated reads, results, indexes and planning review artifacts. Planning and review records are committed in the repository.

## Start a new experiment

See the [simulation experiment guide](https://github.com/berntpopp/MucOneSpan/blob/v0.11.0/docs/guides/simulation-experiments.md), [configuration guide](https://github.com/berntpopp/MucOneSpan/blob/v0.11.0/docs/guides/configuration.md) and [example designs](https://github.com/berntpopp/MucOneSpan/tree/v0.11.0/examples). Simulator executables and platform/model configurations must be supplied explicitly. Keep generation truth separate from calling, and freeze settings before generating a final validation panel.

## Validation and limits

Local checks: 696 unit tests pass on Python 3.12; the 694-test suite before the final test-fixture fix also passed on Python 3.10–3.14. Primary coverage is 93.48% with branches. Integration: 47 passed. Explicit generated-data end-to-end tests: 6 passed and 2 unchanged strict expected failures. Lint, formatting, configured mypy, file-size/workflow checks, documentation, dependency audit and distribution/install checks pass.

The equal-length 60/60 and asymmetric 25/140 length failures remain. Optional read-backed phasing remains disabled by default after a development false-positive regression. Confidence scores are heuristic; VCF concordance is not independent read support. Stricter evidence fields can distinguish unresolved reconstruction from a confident negative.

The planned fresh final panel was **not executed**. Original samples and later challenges are development evidence. This release does not claim a general increase in complete diploid reconstruction accuracy or clinical performance. See [per-cohort results and rejected experiments](https://github.com/berntpopp/MucOneSpan/blob/v0.11.0/docs/guides/validation-results.md) and the repository's [requirement/evidence matrix](https://github.com/berntpopp/MucOneSpan/blob/v0.11.0/.planning/requirement-evidence-matrix.md).
