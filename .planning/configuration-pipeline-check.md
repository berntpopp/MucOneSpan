# Focused full-pipeline runtime-configuration regression

## Purpose and scope

The ignored driver
`tests/results/production_validation_20260914/configuration/full_pipeline_check.py`
is prepared for a single post-freeze launch. It compares the immutable
preconfiguration source snapshot in `configuration/before` with a one-time copy
of the final candidate `src` tree in `configuration/after`. It runs the existing
`fresh_generation/paired_worker.py` under the candidate worktree's locked uv
Python, the same explicit Clair3 tool environment and platform-specific models,
and four threads.

The driver does not read simulator truth or encode expected biological calls.
It hashes the selected development inputs, both source snapshots, the worker and
driver, locked environment files, external executables, and model files. It
alternates before/after order across cases and executes each pair sequentially.
Any failed or missing worker result remains an explicit execution failure.

## Frozen input proposal

| Case | Platform and role | Repository-relative input | SHA-256 |
| --- | --- | --- | --- |
| `sample_bench_5003` | HiFi difficult classifier | `tests/data/generated/sample_bench_5003/sample_bench_5003_reads_amplicon_aligned.bam` | `ef78136f129d8e837e8a6b131e5217072b1073d5849282d95e0c3c0ddc6c1515` |
| `sample_normal_60_80` | HiFi normal | `tests/data/generated/sample_normal_60_80/sample_normal_60_80_reads_amplicon_aligned.bam` | `d56f364b38bade513d6cc6ce669a3c810fba0c45187f4e8d9f67c43f6aff6626` |
| `sample_homozygous_60_60` | HiFi equal-length | `tests/data/generated/sample_homozygous_60_60/sample_homozygous_60_60_reads_amplicon_aligned.bam` | `f6657358dda7bacc0f030f45d14964484bef12f202278889e0e25e4f0890a0cf` |
| `sample_ont_normal_60_80` | ONT normal | `tests/data/generated/sample_ont_normal_60_80/sample_ont_normal_60_80_reads_amplicon_ont.fastq` | `e013b578471cf8b98f636dda4dca0e12a21a6a22377493d47144e666cabb4338` |
| `sample_ont_dupc_60_80` | ONT mutation | `tests/data/generated/sample_ont_dupc_60_80/sample_ont_dupc_60_80_reads_amplicon_ont.fastq` | `86dab5304372f3c0ce4a6e605c16cc05101f1a27e848ca7eac4e7ae6a2c99dd7` |

These are cached development inputs and do not consume fresh validation seeds.

## Comparison contract

For each paired result, `comparison.json` records exact trimmed-consensus
sequence length and SHA-256, allele count, reported and canonical repeat counts,
classified repeat count, complete repeat structure, and mutation identity. Read
and primary-alignment counts are retained separately. Mutation confidence and
VCF evidence annotations are separated from mutation identity. A difference is
tagged as the known N1 annotation correction only when mutation identities are
unchanged and the current unsupported status is `projection_unavailable` or
`localization_ambiguous`; the driver does not waive sequence, count, structure,
or mutation-call differences.

New `run_configuration.json`, `run_status.json`, and summary schema paths are
reported under provenance/schema sections. These additive records are not
treated as scientific-output differences.

## Prepared launch command

Run only after no other heavy tool jobs are active and the coordinator confirms
the candidate source freeze:

```bash
cd /home/bernt-popp/development/MucOneSpan/.worktrees/production-validation
uv run --locked --all-extras python \
  /home/bernt-popp/development/MucOneSpan/tests/results/production_validation_20260914/configuration/full_pipeline_check.py \
  run \
  --candidate-root /home/bernt-popp/development/MucOneSpan/.worktrees/production-validation \
  --data-root /home/bernt-popp/development/MucOneSpan/tests/data/generated \
  --pv-root /home/bernt-popp/development/MucOneSpan/tests/results/production_validation_20260914 \
  --tool-bin /home/bernt-popp/miniforge3/envs/env_clair3/bin \
  --hifi-model /home/bernt-popp/miniforge3/envs/env_clair3/bin/models/hifi \
  --ont-model /home/bernt-popp/miniforge3/envs/env_clair3/bin/models/r941_prom_hac_g360+g422 \
  --timeout-executable /usr/bin/timeout \
  --timeout-seconds 1800 \
  --output /home/bernt-popp/development/MucOneSpan/tests/results/production_validation_20260914/configuration/full_pipeline_results
```

Candidate snapshot creation and paired jobs are deliberately coupled in `run`;
an existing `configuration/after` or output directory aborts a new invocation
rather than overwriting evidence.

## Completed execution

The coordinator authorized launch after other tool-heavy checks completed. All
ten paired-worker invocations completed with exit code 0. The comparison report
records:

- all five cases have exact equality for trimmed consensus sequence, allele and
  repeat counts, complete classified structure, mutation identity, selection
  fields, mutation annotations, read counts, and primary-alignment counts;
- the difficult HiFi and mutant ONT cases retain the same `dupC` mutation call;
- all other case/allele mutation lists remain empty on both snapshots;
- there are no unexpected-difference samples and no N1 exception was needed;
- the current snapshot adds 63 summary schema paths per case and removes none;
- `run_configuration.json` is present only in current-snapshot results, while
  both snapshots emit `run_status.json`. These are recorded separately from the
  scientific comparison.

Artifacts:

```text
configuration/full_pipeline_results/manifest.json
  SHA-256 a11eb7464bfed1eed69220e23c5c7cf0bba1f68e86ef5d25c21d70b59bb69d69
configuration/full_pipeline_results/comparison.json
  SHA-256 6314aeeb6c439dd19838bd4a262353fe58b7edc9aa3d6f893d062e73a0684dc9
configuration/after/manifest.json
  SHA-256 d4172d15b9f7fc28641e820c77badc79241c10948f0b82d39dd9eca2a9f252bd
```

The first difficult pair completed before the freeze guard stopped on `.fai`
and `.mmi` indexes generated beside the copied bundled reference. Those files
are tool artifacts excluded by the established paired freeze policy. The driver
was repaired to exclude them and resumed from the preserved snapshot; it did not
rerun the completed pair. The manifest records the initial and repaired driver
hashes and resume reason. A later comparison-only repair restricted allele keys
to `allele_<integer>` so additive `allele_*` metadata was not parsed as an allele;
all completed worker records were skipped during that pass.

Wall times ranged from 6.10 to 8.91 seconds. They are descriptive smoke-test
measurements, not a controlled speed comparison: package reference indexes were
cold/warm in different order during the first pair. No fresh validation inputs
or seeds were generated.
