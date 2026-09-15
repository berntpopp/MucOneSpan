# Wave 1 local delivery

Branch: `fix/wave1-mapping-reports`

Worktree: `/home/bernt-popp/development/MucOneSpan-wave1`

Base: `1f6c165599fa439bb40cb8d74ec573c241fb88dc`

At the initial implementation handoff, no commit, push, merge, release, GitHub
comment or issue closure had been performed. The user subsequently authorized
committing this work, pushing the branch and opening a pull request.
The main checkout's unrelated #49 tracked binary diff is identical to its initial
snapshot, and all nine initial untracked file SHA-256 checks passed.

## Proposed PR title

fix: bound mapping process lifecycle and correct report evidence/status

## Proposed PR description

Mapping now keeps minimap2 and samtools streaming under a finite configurable
budget and cleans up launched processes and partial BAM/index outputs on failures,
timeouts and interruptions. An isolated Linux supervisor adopts and reaps orphan
descendants. HiFi/ONT presets and explicit overrides are preserved.

Reports open on actual allele contigs and include all distinct allele VCF files,
with shared files represented once. Existing single-VCF callers remain compatible.
Execution provenance prevents known unfinished, insufficient, failed or interrupted
runs from showing negative results; report failures remain execution failures,
including in portable summaries. Missing nomenclature stays safely escaped and
exact_match_pct is consistently rendered in its established 0–100 percent units.

Validation: 874 unit tests; 91.32% branch-aware coverage; 69 integration tests;
11 real Chrome tests; quality, strict docs and distribution checks passed. All
84 scientific artifact hashes match across three before/after real pipelines.
References #41, #39, #40. Scientific classification policy and thresholds unchanged.

External compatibility: igv-reports 1.13.0 verified. Malformed VCF payloads from
available 1.16.0 are rejected clearly. Non-Linux POSIX retains group cleanup with
orphan reaping delegated to the OS; no non-Linux runtime verification was available.

## Issue-specific completion evidence

| Issue | Reproduced defect | Fix | Regression evidence | Actual verification / remaining limit |
| --- | --- | --- | --- | --- |
| #41 | Raw pipe waits/joins lacked deadlines and group cleanup; launch/peer failure left partial outputs. Review reproduced orphan zombies and input/output collision deletion. | Configurable shared budget, cleaned PATH, streaming pipe, isolated Linux subreaper, bounded termination/reaping/drains, both diagnostics and artifact cleanup; collision rejected before source deletion. | test_mapping_lifecycle.py, test_tool_pipeline.py, test_tool_supervisor.py; ten real handshake scenarios assert every recorded PID disappears, including zombies; presets/override tests. | 78 focused units and ten real process cases passed. Real minimap2/samtools integration and three scientific runs passed. Non-Linux fallback is unit-tested only; uninterruptible kernel tasks surface explicit incomplete cleanup. |
| #39 | Generated session opened contig_1 despite contig_11/21 assignments; only first allele VCF reached report. | summary allele contig_name, all paths through CLI/pipeline/report/IGV, stable labeled/deduplicated tracks, explicit missing inputs/contigs and malformed-payload failures. | test_report_igv.py, test_report_sessions.py and test_report_igv_browser.py cover two/single/shared/absent VCFs, singular/plural precedence, embedded/sidecar. Live variant coordinates/ref/alt checked. | 12 real session tests and 11 browser tests passed; both real standalone CLI modes generated. igv-reports 1.13.0 verified; 1.16.0 malformed output rejected with clear error. |
| #40 | Percent 0.5 displayed 50%, 1 displayed 100%; missing nomenclature emitted literal entity; known failed sidecar ignored; failed rendering left misleading summary state. | Explicit 0–100 boundary, Unicode missing text with autoescape, authoritative lifecycle/sidecar/portable-summary provenance. Valid mutation evidence retained with execution warnings. | test_report_wave1.py, test_wave1_cli.py and test_run_status.py; missing/null/empty/hostile values, all six percentages, invalid values, full/standalone lifecycle, conflicting/legacy status and report failure. | Unit and real browser checks passed. Adjacent biological decision-policy defects are reproduced separately and unchanged; no clinical threshold or tier changes. |

## Final verification commands and results

All commands ran in the isolated worktree on the final source. Logs remain in
`/tmp/wave1-*-final.log`; original red evidence is retained separately.

| Command | Result |
| --- | --- |
| `make test-fast` | 874 passed in 6.13s; zero skips. |
| `make quality` | Ruff/format, configured mypy (63 files), file-size (149 files, max allowed 649), actionlint passed. |
| `make ci-check` | Quality passed; 874 passed in 14.59s; 91.32% branch-aware coverage, exceeding 80%. |
| `make test-int` with external tool PATH below | 69 passed in 16.44s, eight deselected, zero skips. The eight e2e-marker tests are outside this target; three real full pipelines separately tested the representative panel. |
| `make docs-check` | Strict documentation build passed. |
| `uv run --locked --all-extras pytest tests/browser --no-cov` with external tool PATH | 11 passed in 5.04s, zero skips; real Chrome offline loci/tracks/variants, percentages and existing accessibility checks. |
| `make build-check` | Wheel/sdist resource checks, isolated wheel install, CLI smoke passed. |
| `git diff --check` | Passed. |

For the external checks, the executable PATH begins with:

```bash
PATH=/home/bernt-popp/miniforge3/envs/env_clair3/bin:/home/bernt-popp/miniforge3/envs/vntyper/bin:$PATH
```

The first environment supplies minimap2 2.28-r1209, samtools 1.15.1, bcftools 1.17
and Clair3 v1.0.10 with HiFi/ONT models; the second supplies create_report from
igv-reports 1.13.0. These paths describe this machine's validation environment,
not authored runtime configuration. No project dependency or lockfile changed,
so the conditional security gate was not applicable.

## Scientific preservation

After the final supervisor cleanup-reserve correction, the lead ran the full
three-fixture panel again under `/tmp/wave1-scientific-baseline/verified`:
HiFi normal, HiFi dupC and ONT dupC. All three exited 0 and generated reports.
The panel compare command reported 84 before artifacts, 84 after artifacts,
and `differences: []`. Input hashes match. Complete allele/repeat/mutation,
full/trimmed consensus, genotype/phase and intermediate/final VCF evidence is
unchanged. Normal retains no mutations; both mutation fixtures retain dupC at
repeat 25 and 60/80-repeat alleles on contig_51/71.

See `wave1-scientific-evidence.md` and `wave1_scientific_panel.py` for commands,
inputs/models, raw output paths and documented nondeterministic normalization.
Prior baseline/fixed/fixed-final artifacts remain available.

## Independent review and final audit

Independent reviewer approved specification compliance and code quality for Linux
with the verified IGV tool, resolving all five findings and two supervisor followups.
The reviewer additionally exercised an escaped-session descendant and verified its
PID vanished while an unrelated child stayed alive. See `wave1-review.md`.

Final audit covered every numbered task requirement: clean isolation and #49
preservation; red reproductions; finite streaming lifecycle and process cleanup;
all report loci/tracks and real browser evidence; escaping/status/percentage
contracts; scientific artifact invariance; every required gate; docs/changelog;
independent review; archived plan; and this issue table/PR handoff. No essential
Linux/tool/browser verification remains pending.

## Explicit limitations and adjacent observations

- Non-Linux POSIX group cleanup is preserved but orphan reaping is delegated to
  the OS; it was not exercised on a non-Linux host. Kernel-uninterruptible tasks
  cannot be guaranteed to reap by a finite deadline and produce explicit errors.
- Source/output BAM or index collisions are rejected before deleting source data.
- External igv-reports 1.16.0 malformed VCF payloads are rejected; use a working
  create_report implementation (1.13.0 verified). No scientific data are repaired
  or rewritten to hide the external defect.
- Existing clinical decision behavior treats any mutation record as pathogenic
  and can label empty legacy summaries negative. Executed reproductions are
  recorded in `wave1-report-evidence.md`; broader evidence-policy changes remain
  outside this plumbing wave.
- This small panel makes no 500-sample accuracy, independent diploid phase,
  clinical sensitivity or universal speed claim.
