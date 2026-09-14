# Fresh generation driver prepared 2026-09-14

Status: executable prepared, hand-built checks and both actual development
capture smokes pass. No reserved final seed has been generated. Coordinator must
freeze caller/settings/scoring before final-panel launch. Driver resides in ignored
`tests/results/production_validation_20260914/fresh_generation/`. The worktree's
`tests/results` resolves to the original checkout's shared ignored directory.

## Fixed cases and reproducibility

`prepared_v2/prepared_manifest.json` SHA-256:
`67ffa423c2814ebc167cb2490f02573e1a6121bad3300f223da561a1d17df652`.
It contains 107 SHA-256 fingerprints: original simulator Python source,
pyproject/lock, simulator interpreter, generator/capture/audit/check scripts,
copied generation configuration, used tool executables, pbsim Conda metadata,
and both platform error models. Version probes are retained verbatim.

Predefined designs follow the production spec in its listed order. HiFi uses
2026091401–1412; ONT uses 2026091413–1424. Four new designs per platform follow:
ins25bp 60/80 H1 repeat25; non-X insC_pos23 60/80 H2 repeat25;
delGCCCA 60/80 H1 repeat25; dupC 60/60 H2 repeat25. Extension seeds are
2026091425–1428 HiFi and 2026091429–1432 ONT. The insC_pos23 definition permits
only A/E parents, guaranteeing a non-X mutant parent; simulator permissive
conversion selects the allowed parent with its existing seeded semantics.
The configured ins25bp sequence contains exactly 25 inserted bases.

All request 200 initial templates except each platform's normal low-depth
60/80 case requesting 20. The copied config preserves joint PCR preset default,
HiFi ERRHMM-SEQUEL with ten passes/minimum three passes/minimum RQ 0.99,
and ONT QSHMM-ONT-HQ with mean accuracy 0.95 and one pass. Both use two nominal
simulator threads. No human-reference alignment is requested. Metadata's generic
coverage field is not interpreted as achieved coverage.

Recorded read seed offsets follow unchanged source: PCR base seed, pbsim seed
plus 1/2 for H1/H2, CCS seed plus haplotype*100 plus BAM index. Consecutive base
seeds can therefore share a component seed across different designs; these
simulated designs are not claimed to be population-independent replicates.

## Observational capture and caller isolation

`capture_cli.py` imports the original sibling CLI in its original virtualenv,
then wraps only merge boundaries in that child process. It never edits sibling
source or calls separate haplotype simulations. HiFi capture copies each
`hifi_hap{1,2}_NNNN.bam` before the joint BAM merge and runs the same samtools
FASTQ conversion arguments for inspection. ONT capture copies each
`ont_hap{1,2}_NNNN.fq.gz` immediately before the joint FASTQ merge. Both then
call the original merge function unchanged. Failure to observe an expected
merge boundary fails generation rather than inventing source assignment.

`records.py` compares exact sequence/quality AND original-header multisets
between all captured records and the simulator's emitted main FASTQ. It checks
HiFi primary BAM counts against captured usable FASTQ counts. Any lost record
or changed multiset fails the case. It retains all intermediate copies, original
headers, per-source row numbers, exact source haplotypes, and globally unique
source IDs. The main analysis input is the simulator's original-name FASTQ.
Caller commands receive only the selected FASTQ; all source/sequence truth
sidecars are strictly offline evaluator input.

`unique_names.fastq` and `captured_original_names.fastq` are an exact paired
renaming comparison with identical record order and sequence/quality.
The emitted main FASTQ remains separately available; merge ordering can differ.
When identical original name/sequence/quality occurs across haplotypes, its main
record cannot be assigned uniquely after the merge: `main_source_candidates.json`
records null haplotype and all candidates. No arbitrary candidate is promoted to
source truth. Captured source IDs remain exact because they precede the merge.

Each successful case produces three deterministic perturbations, 96 total for
32 successful simulations. Minority H1/H2 keeps all opposite-haplotype records
and target-haplotype records whose SHA256(seed:kind:source_id) is 0 modulo 10.
Partial-read perturbation clips half the sequence/quality from approximately half
the records, choosing prefix/suffix by the same hash. Retained intervals and
actual source identity are recorded per emitted record. These use captured
original names and preserve source sequence/quality slices; they are dependent
within-sample perturbations, not new simulations. Low source support can produce
zero minority records, which remains explicit in counts.

## Launch only after coordinator source/settings freeze

From `/home/bernt-popp/development/MucOneSpan/.worktrees/production-validation`:

```bash
uv run --locked --all-extras python \
  tests/results/production_validation_20260914/fresh_generation/generate.py run \
  --prepared tests/results/production_validation_20260914/fresh_generation/prepared_v2 \
  --prepared-sha256 67ffa423c2814ebc167cb2490f02573e1a6121bad3300f223da561a1d17df652 \
  --freeze COORDINATOR_FREEZE_JSON \
  --freeze-sha256 COORDINATOR_FREEZE_SHA256 \
  --output tests/results/production_validation_20260914/fresh_generation/data
```

Substitute the coordinator's actual freeze path and its explicit immutable SHA.
Both hashes and every prepared dependency are checked before output creation.
An existing output root is rejected, including partial runs. There is no resume
or overwrite mode. The run manifest lists all expected cases before generation,
records exact CLI/capture commands and retained artifact hashes, preserves failed
cases, and exits nonzero if any case fails. Current caller source validation is
the coordinator's freeze responsibility; this driver verifies and archives that
freeze file's bytes but does not interpret its schema or silently create one.

## Executed checks and remaining validation

- `uv run --locked --all-extras python .../check_driver.py`: PASS for 32-case
  seed contract, duplicate QNAME records, genuinely ambiguous source records,
  exact sequence/quality multisets, globally unique IDs, clipping and overwrite
  refusal on hand-built data. It uses fixture seed 7, not a reserved simulation.
- `uv run --locked --all-extras ruff check .../fresh_generation
  --no-respect-gitignore`: all checks passed.
- Matching Ruff format check: four files already formatted. Authored scripts
  range from 62 to 301 lines, all below the 650-line limit.
- Both command parsers' `--help` pass; wrapper `--version` was invoked through
  project `run_tool` under the simulator interpreter and returned MucOneUp 0.44.5.
  Exact command and wrapper log are retained locally.

The first wrapper version probe exposed an incorrect entry-point assumption:
the package's `main()` is a zero-argument launcher, while `click_main.cli`
accepts explicit Click arguments. The corrected import passed the real probe.
The original prepared manifest is superseded and intentionally retained;
`prepared_v2` fingerprints the corrected wrapper. No scientific data existed at
either preparation, so no panel was exposed or reclassified.

## Actual development capture smoke

Coordinator explicitly authorized seed 19200901, lengths 25/30, twelve requested
templates on both platforms, with no calling or settings tuning. A separate
two-case prepared manifest and freeze record identify this as development-only.
No final seed was used. The smoke exercised the existing driver without code or
configuration fixes; final `prepared_v2` and its hash remain current.

Executed command:

```bash
uv run --locked --all-extras python \
  tests/results/production_validation_20260914/fresh_generation/generate.py run \
  --prepared tests/results/production_validation_20260914/fresh_generation/prepared_smoke_19200901 \
  --prepared-sha256 e4025c3d85fadecb2ddbd3bdd4bc62afaca5e9818d3c69ae783b16dae0a787a8 \
  --freeze tests/results/production_validation_20260914/fresh_generation/prepared_smoke_19200901/development_smoke_freeze.json \
  --freeze-sha256 fca2ffb21b567515f79fb96a2540beb64daa5a9d25afbc388a5b1e1cadbf13d3 \
  --output tests/results/production_validation_20260914/fresh_generation/development_smoke_19200901
```

Exit status 0; HiFi 2.377 seconds, ONT 0.334 seconds. Both actual patched merge
boundaries captured two haplotypes and passed the exact sequence/quality and
original-header multiset comparisons against the emitted main input.

| Platform | Requested templates | Usable main records | Source H1/H2 | Duplicate-name records | Ambiguous main source |
| --- | --- | --- | --- | --- | --- |
| HiFi | 12 | 7 | 5/2 | 1 | 0 |
| ONT | 12 | 12 | 7/5 | 5 | 0 |

HiFi retained both actual pre-merge CCS BAMs; primary BAM counts (5/2) equal the
captured FASTQ counts. ONT retained the original haplotype `.fq.gz` files (7/5).
Both retain source maps, unique-name and original-name pairs, main-input maps,
all generation commands/logs, and complete artifact hashes in the generation
manifest. The caller was never invoked.

All six dependent smoke perturbations were produced with exact source truth.
HiFi minority-H1/H2 retained 4/6 records and the partial input 7; ONT retained
7/8 and 12 respectively. Actual source counts remain in each perturbation record;
at these small counts the intended random thinning is not an exact 10% fraction.
The post-smoke audit verified all 107 final prepared dependency hashes still
match. Final data output remains absent. These small fixtures validate execution
and capture mechanics, not final-panel feasibility for every design or scientific
accuracy. No caller or simulator production source was modified.
