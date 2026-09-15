# Wave 1 scientific preservation panel

## Baseline provenance

- Baseline source: clean `main` commit `1f6c165599fa439bb40cb8d74ec573c241fb88dc`,
  exported with `git archive` to `/tmp/wave1-scientific-baseline/source` before
  implementation. No uncommitted main-checkout code was copied.
- Existing reads are read-only inputs from the main checkout's
  `tests/data/generated/`; generated outputs remain under
  `/tmp/wave1-scientific-baseline/` and are not committed.
- Tools: minimap2 `2.28-r1209`, samtools `1.15.1`, bcftools `1.17`,
  Clair3 `v1.0.10`, from `/home/bernt-popp/miniforge3/envs/env_clair3`.
  Clair3 uses the adjacent `bin/models/hifi` or `bin/models/ont` model.
- Project interpreter/dependencies: Wave 1 worktree `.venv`, invoked through
  `uv run --locked --all-extras`; `PYTHONPATH` selects the archived source.
  The runner places the Conda environment on external-tool PATH.
- All runs use two threads, platform defaults (`map-hifi` / `lr:hq`), normal
  pipeline scientific defaults, and HTML reporting with IGV off.

## Baseline results

All three complete CLI runs exited 0 and produced `run_status.json` with
`completed`, two consensus candidates, and an HTML report. They produced 84
canonical scientific artifact hashes in `baseline/scientific-sha256.json`.

| Fixture | Platform | Selected total repeats / contigs | Consensus lengths (bases) | Mutation identities / total-repeat index | Filtered VCF records |
| --- | --- | --- | --- | --- | --- |
| `sample_normal_60_80` | HiFi | 60/80; `contig_51`/`contig_71` | 3600 / 4800 | None | 59 / 54 |
| `sample_dupc_60_80` | HiFi | 60/80; `contig_51`/`contig_71` | 3601 / 4800 | allele 1 `dupC`, repeat 25 | 61 / 108 |
| `sample_ont_dupc_60_80` | ONT | 60/80; `contig_51`/`contig_71` | 3601 / 4800 | allele 1 `dupC`, repeat 25 | 61 / 108 |

Normal allele 2 retains `heterozygosity_observed` / `unphased`; the other five
candidates retain `no_heterozygosity_observed` /
`no_informative_heterozygosity`. All three retain
`distinct_genotype_candidates` VNTR evidence. This records existing evidence;
it does not assert independently phased diploid reconstruction.

## Reproduction and comparison

From the Wave 1 worktree, substitute equivalent local paths as necessary:

```bash
mkdir -p /tmp/wave1-scientific-baseline/source
git archive 1f6c165599fa439bb40cb8d74ec573c241fb88dc | tar -x -C /tmp/wave1-scientific-baseline/source
uv run --locked --all-extras python .planning/wave1_scientific_panel.py run \
  --source /tmp/wave1-scientific-baseline/source \
  --data /home/bernt-popp/development/MucOneSpan/tests/data/generated \
  --tools /home/bernt-popp/miniforge3/envs/env_clair3 \
  --output /tmp/wave1-scientific-baseline/baseline
uv run --locked --all-extras python .planning/wave1_scientific_panel.py run \
  --source "$PWD" \
  --data /home/bernt-popp/development/MucOneSpan/tests/data/generated \
  --tools /home/bernt-popp/miniforge3/envs/env_clair3 \
  --output /tmp/wave1-scientific-baseline/fixed-final
uv run --locked --all-extras python .planning/wave1_scientific_panel.py compare \
  /tmp/wave1-scientific-baseline/baseline /tmp/wave1-scientific-baseline/fixed-final
```

The runner requires a new output directory to prevent stale artifacts. Each
sample has a 900-second driver timeout, a complete CLI log, and a run record
with arguments, input SHA-256, and actual exit status.

Comparison includes complete `alleles.json`, `repeats.json`, `repeats.txt`, all
full and trimmed consensus FASTAs, all consensus-context JSON, all intermediate
and final compressed/plain VCFs (including genotype FORMAT/sample fields and
Clair3 phased VCFs), and complete summary `alleles`, `classifications`, tool
versions, and pipeline version. JSON key ordering is canonicalized, compressed
VCFs are decompressed, and only checkout/output root paths plus VCF
`##bcftools_*`, `##fileDate=`, and `##cmdline=` provenance headers are normalized.
Record order, alleles, qualities, genotype/phase evidence, sequences, and mutation
identities are retained. Configuration/status/report provenance is intentionally
outside the scientific comparison because this wave changes those contracts.
The runner reports missing/extra/different artifacts and exits nonzero on any
difference. Raw artifacts and metadata remain available for inspection.

## External-report prerequisite observation

An initial baseline normal run with `--report-igv embedded` completed every
scientific stage and failed at report rendering with
`FileNotFoundError: Tool not found: create_report`. The sidecar retained
`execution_failed`, as expected. That separate reproduction is preserved at
`/tmp/wave1-scientific-baseline/before/sample_normal_60_80` with log
`normal-before.log`; it is excluded from the completed scientific panel.
`create_report` was absent from the Conda and project environments at that time.
Real IGV session/browser validation is a separate Wave 1 check.

## Final fixed-run comparison

The authoritative final panel completed against the combined Wave 1 worktree
after the Linux supervisor, source-collision/cleanup/diagnostic fixes, and
terminal-summary synchronization were in place. The previous successful fixed
panel remains preserved separately; it is superseded by this fresh final run.
All three CLI runs exited 0, generated reports, and recorded `completed` in both
the execution sidecar and final summary. Input SHA-256 values match baseline
for all three fixtures. The exact comparison command above exited 0 with:

```json
{
  "before_artifacts": 84,
  "after_artifacts": 84,
  "differences": []
}
```

| Check | Actual result |
| --- | --- |
| HiFi normal full pipeline | Baseline and fixed exit 0; all scientific artifacts identical |
| HiFi dupC full pipeline | Baseline and fixed exit 0; all scientific artifacts identical |
| ONT dupC full pipeline | Baseline and fixed exit 0; all scientific artifacts identical |
| Consensus, complete repeat classifications, mutation identities | Identical after documented normalization |
| Allele selection, genotype/phase evidence, all caller/filtered/phased VCF data | Identical after documented normalization |
| Run status | All six sidecars completed; fixed summaries additionally record completed |
| Input provenance | All three input hashes identical between baseline and fixed runs |

Fixed raw outputs, full logs, exact commands, input hashes and exit codes are
under `/tmp/wave1-scientific-baseline/fixed-final`; the original baseline remains
under `/tmp/wave1-scientific-baseline/baseline`, and the earlier fixed outputs
remain under `/tmp/wave1-scientific-baseline/fixed`. The complete final scientific
comparison passed with no missing, extra, or differing artifacts. This small
panel provides preservation evidence for the selected fixtures; it is not a
500-sample accuracy or speed benchmark and does not replace separate lifecycle
or real IGV/browser verification.

## Final source confirmation after last reviewer correction

The lead repeated all three full runs after the final emergency supervisor reap
reserve fix, output `/tmp/wave1-scientific-baseline/verified`. All exited 0.
`uv run --locked --all-extras python .planning/wave1_scientific_panel.py compare
/tmp/wave1-scientific-baseline/baseline /tmp/wave1-scientific-baseline/verified`
exited 0 and reported 84 before artifacts, 84 after artifacts, differences [].
This is the latest authoritative source confirmation; prior panels are preserved.
