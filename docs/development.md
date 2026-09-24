# Development

The repository uses one shared development contract in
[`AGENTS.md`](https://github.com/berntpopp/MucOneSpan/blob/main/AGENTS.md).
Claude Code and Gemini CLI import that file through their small vendor entry
files; Codex reads it directly. Directory-level `AGENTS.md` files cover package,
test, and script details. Change shared policy in its source rather than copying
it into each assistant's entry file.

## Environment and commands

Use Python 3.10 or newer and uv. From the repository root:

```bash
make dev
make ci-check
```

`make dev` synchronizes `uv.lock` with all dependency groups and extras, including development, report,
and documentation dependencies. Normal development uses the lock without
updating dependency resolution. Individual Make targets install only their required
group (`quality`, `test`, or `security`) or extra (`docs` or `report`); developer
tooling is managed with uv dependency groups rather than a published `dev` extra.
The local `dev` group includes all tool groups and pre-commit.

When intentionally changing dependencies, edit
`pyproject.toml`, run `make lock`, review the lock diff, and run `make dev` again.

| Command | Purpose |
| --- | --- |
| `make quality` | Ruff, formatting, configured mypy, file size, and workflow syntax |
| `make test-fast` | Unit tests without coverage or external tools |
| `make test-unit` | Unit tests with the 80% coverage gate |
| `make ci-check` | Quality and unit tests with at least 80% coverage |
| `make test-int` | Tests marked as tool-dependent integration tests |
| `make docs-check` | Strict documentation build |
| `make security-check` | Audit all locked extras against published Python advisories |
| `make build-check` | Distribution build and package validation |
| `make docker-test` | BuildKit runtime checks without image export/load |
| `make docker-build docker-smoke` | Build and test an image through Docker |
| `make format` | Apply Ruff formatting |
| `make lint-fix` | Apply supported Ruff fixes; review the resulting diff |

Ruff covers `src/`, `tests/`, and `scripts/`. Mypy covers the package and scripts
using the settings in `pyproject.toml`, including required function annotations.
`make quality` also validates GitHub Actions syntax and expressions.
Coverage includes branch measurement;
the 80% gate applies to the aggregate unit test result, not each individual file.
Do not weaken checks or swallow failures to make a change pass.

For a focused test:

```bash
uv run --locked --all-extras pytest tests/unit/test_config.py --no-cov
```

Install the repository hooks with `make hooks`. Hooks provide local feedback;
the checked-in Makefile targets and CI remain the common verification contract.

## Architecture

MucOneSpan draws on some ideas from Vrbacka et al. and supports HiFi and ONT
amplicons. The current implementation uses minimap2. The historical bwa-mem
description is not the implementation contract; see
[differences from the published method](getting-started/deviations.md).

| Component | Responsibility |
| --- | --- |
| `cli.py` | Click commands and pipeline orchestration |
| `config.py`, `data/repeats/` | Repeat dictionary, mutation definitions, packaged data |
| `ladder.py`, `data/reference/` | Synthetic reference contigs and bundled ladder |
| `mapping.py` | Read alignment, BAM processing, alignment statistics |
| `alleles.py` | Repeat-count clusters and allele detection |
| `calling.py`, `vcf.py` | Allele read extraction/remapping, Clair3 calls, VCF processing |
| `consensus.py` | Allele consensus sequences using bcftools |
| `classify.py`, `classify_types.py`, `repeat_alignment.py` | Repeat segmentation, nomenclature, mutation interpretation |
| `report.py`, `templates/` | Structured results and optional HTML report |
| `tools.py` | External command execution, environments, errors, tool versions |
| `scripts/` | Reference generation, simulation, benchmark and maintenance helpers |

Package resources live under `src/muc_one_span/data/` and
`src/muc_one_span/templates/`; there is no top-level runtime `data/` directory.
Resource changes need distribution checks, since an editable installation can
hide packaging omissions.

Authored code, configuration, and templates must stay **below 650 physical
lines**, including comments and blank lines. The gate checks tracked and
untracked nonignored files with these suffixes: `.py`, `.sh`, `.yml`, `.yaml`,
`.toml`, `.css`, `.js`, `.ts`, `.html`, and `.j2`, plus Makefiles and Dockerfiles.
Generated/data files, lockfiles, and prose are outside the gate. Split growing
modules into cohesive responsibilities before reaching 650. Retain public
imports where compatibility matters and retain meaningful tests when splitting
test files. Dense formatting and blanket exclusions defeat the purpose.

## Scientific contracts

- MUC1 repeat units are normally 60 bases. Ladder contig names encode canonical
  variable repeats; the pre- and after-repeat blocks are additional sequence.
  The generator defaults to `contig_1` through `contig_150` with 500-base flanks.
- The ladder combines flanks, pre-repeats 1–5, canonical repeats, and after-repeats
  6–9. Dictionary nomenclature and source provenance must remain traceable.
- Dictionary mutation coordinates use 1-based positions; Python sequence slices
  use 0-based positions. Check insertion anchors and inclusive deletion endpoints
  explicitly when changing coordinate conversion or VCF handling.
- Alleles may have equal or nearby repeat counts, and PCR bias can yield very
  unequal coverage. Keep expected behavior explicit for homozygous, close-length,
  asymmetric, low-coverage, and long-repeat cases.
- HiFi defaults to minimap2 `map-hifi`; ONT defaults to `lr:hq`. Platform selection
  also reaches Clair3. Preserve explicit `--minimap2-preset` overrides.
- Changes to thresholds, confidence scores, mutation naming, or allele assignment
  require behavior-specific evidence. Maintenance refactors should preserve them.

See [core concepts](getting-started/concepts.md),
[repeat nomenclature](reference/nomenclature.md), and
[limitations](reference/limitations.md) for domain background. Historical
benchmark observations are snapshots of a stated version and dataset, not
guarantees about a future change.

## Test layers and external tools

Unit tests run without bioinformatics executables or generated sequencing data.
Mock command calls at the boundary and assert returned results, command arguments,
errors, and coordinate behavior with small deterministic fixtures.

Integration tests require the tools exercised by each test: minimap2 and samtools
for alignment, bcftools for consensus/VCF work, and Clair3 with a suitable model
for variant calling. MucOneUp and pbsim3 are generation prerequisites; they are
not prerequisites for every integration test once reads exist. The conda tool
environment is separate from the Python environment created by uv.

```bash
conda env create -f conda/environment.yml
conda activate muconespan-tools
make test-int
```

The optional MucOneUp end-to-end tests use existing reads in
`tests/data/generated/`. With Clair3 and its dependencies on `PATH`, run:

```bash
# Preserve Clair3's Python interpreter on PATH while pytest uses the project env.
.venv/bin/python -m pytest tests/integration --no-cov
```

Clair3 tests discover its adjacent `models/hifi` directory or accept an explicit
`CLAIR3_MODEL` directory. The homozygous 60/60 and asymmetric 25/140 allele cases
have strict expected failures for the [documented limitations](reference/limitations.md).
They keep the original expected counts and tolerance; an unexpected pass requires
reviewing the limitation rather than silently retaining an obsolete expectation.

Inspect skip reasons. A successful invocation in an environment missing tools,
models, or generated data does not demonstrate that those behaviors work.
Report the commands, pass/fail/skip counts, and missing prerequisites in a PR.
See [benchmarking](guides/benchmarking.md) for portable generation and batch runs.

## CI and container efficiency

The workflow design was compared with the sibling `hum-clinical-reporting`
repository. Both benefit from locked dependency layers, small build contexts,
lock-keyed caches, explicit runner versions, and separate routine checks from
expensive tool checks. This pipeline retains Debian for its Conda/TensorFlow
bioinformatics stack; the reporting application's Alpine runtime is specific to
its own dependency set.

- Each Make target selects its required uv group or extra. A fresh Python 3.10
  unit environment installs 16 distributions, compared with 76 when all tooling
  and documentation dependencies were installed in every job.
- All five supported Python versions run the unit suite for runtime changes.
  Python 3.14 alone measures coverage; the other four avoid redundant coverage
  instrumentation and HTML/XML generation.
- PR file filters skip unrelated integration, package, dependency-audit, and
  container work. The always-running CI Gate rejects failures and unexpected
  skips. Branch pushes and manual runs execute all applicable checks.
- The solved integration Conda environment is cached by its specification.
  Python caches depend on `uv.lock`, so a version-only edit does not discard
  third-party packages. The editable project's metadata cache separately tracks
  `version.py`, keeping the installed version accurate after a bump.
- The Docker tools layer depends only on the pinned base and Conda inputs. The
  locked application and report dependencies live in a separate small virtualenv;
  external tool calls use the Conda interpreter. Source edits preserve the large
  tool layer. The smoke test verifies both interpreters, report rendering, and
  real Clair3 inference.
- PR container builds restore cache only and run the smoke test in a BuildKit
  test stage, without exporting/loading a Docker image. This stage adds no
  dependencies to the runtime it tests. Trusted main builds publish that same
  runtime ancestor and export final layers without duplicate builder environments.

Before optimization, the measured GitHub container build step took 424 seconds:
197 seconds were cache export and 119 seconds were image export/loading. The PR
BuildKit test target avoids both export tasks. Cold builds still need
to obtain the bioinformatics tools and models; their size is an inherent cost,
so compare cold and cached builds separately when measuring future changes.

## Changes and review

Keep implementation plans and decisions in `.planning/`, with completed plans
in `.planning/archive/`. Include the behavior being preserved, module boundaries,
relevant checks, and any remaining limitations. Plans are working documents;
developer instructions and published docs are the maintained reference.

Before review, inspect the diff for unrelated changes, generated artifacts,
hardcoded local paths, and compatibility breaks. Run `make ci-check` and any
additional checks relevant to the change. Update the changelog for user-visible
changes. PR descriptions should explain the concrete result and include actual
validation evidence, including skips and follow-up work.

## Evaluation and evidence contracts

Use `scripts/evaluate.py RESULT_ROOT --truth-root TRUTH_ROOT
--expected-samples INVENTORY.json --output REPORT.json` for offline scientific
scoring. The inventory is a JSON array of sample names or objects with `sample`,
`truth_dir`, `result_dir`, `platform`, `input` and optional `run_record` fields.
Explicit inventories retain unattempted samples in denominators. Truth validation
reconstructs each full haplotype from actual repeat and mutated-unit sequences.

The version1 report separates exact event annotation (parent, full mutation name,
1-based total-repeat index on the sequence-matched haplotype), literal full
sequence equality, ordered structure equality, count error, extra/missing alleles,
assignment ambiguity and execution status. Event annotation is not a claim of
normalized arbitrary biological-event equivalence. Every sequence-optimal
assignment is retained; conservative metric bounds govern acceptance. A failed
normal sample is not a true negative, and duplicate sequence observations do not
prove two independently reconstructed haplotypes.

`run_status.json` records `running`, `completed`, `insufficient_evidence` or
`execution_failed`. The latter two keep the existing nonzero CLI outcome; the
status file distinguishes a coverage no-call from a tool failure. Old successful
artifacts cannot override a known failed execution. `alleles.json` is persisted
after calling and consensus so it agrees with final summary evidence.

Consensus uses explicit sample selection and `-H I` for unresolved mixed calls.
Single-site heterozygosity permits an unordered pair; multiple heterozygous loci
require one common phase set before genotype-index haplotypes are emitted. This
implements the [bcftools consensus selectors](https://samtools.github.io/bcftools/bcftools.html#consensus)
with additional phase checks; bcftools selectors alone do not establish phase.
Distinct-length candidates, and the two haplotag partitions of experimental
read-backed phase, are length-partitioned haplotypes. After the
allele-fraction rule, any remaining heterozygous record (single-site, phased or
multi-site) or any conflicting record selects `I`. The allele records
`allele_genotype_status` (`heterozygous_within_length_partition` or
`unresolved_genotype_records`) and `heterozygous_sites`, and it receives no
independent haplotype credit. Only candidates without remaining heterozygosity
use genotype index 1 (`allele_specific_resolved`). `variant_filter` records the
filter applied to that allele's calls, including the QUAL threshold used.
`stage_concordance` compares that Clair3 partition's pileup-stage VCF
(`clair3/pileup.vcf.gz`, split and left-aligned with `bcftools norm`) with the
alleles applied by genotype in the filtered VCF. Its fields are `status`
(`concordant`, `discordant_frameshift` or `not_assessed`), `source`
(`clair3_pileup`), the applied `min_af` and `min_depth`, `records` (each
unapplied pileup frameshift with `chrom`, `pos`, `ref`, `alt`, `af` and `dp`)
and, for `not_assessed`, a `reason` (`pileup_vcf_unavailable` or
`final_vcf_unavailable`). Every allele with its own Clair3 partition carries
the record; an unphased same-length alias shares `allele_1`'s. "Applied" means
after MucOneSpan's own QUAL filter, so a low-QUAL pileup frameshift that the
filter removes is also discordant. A discordant allele blocks NEGATIVE and is a
quality caveat on PATHOGENIC; `not_assessed` is logged and reported as a quality
caveat without blocking NEGATIVE. The gate never creates, removes or edits a call.
Experimental read-backed phase is available through the Python library's explicit
`read_phase=True` or JSON `calling.read_phase: true`. It remains disabled by default
after failing the development false-positive gate; there is no dedicated CLI flag. `consensus_context` records the
actual reference, full consensus, selected sample/haplotype and half-open trim
interval. VCF support verifies replay, then exact event reversion in sequence
context. It describes concordance with the same VCF used to make consensus,
not independent experimental support.
Under `-H I`, bcftools writes IUPAC codes for heterozygous SNVs but applies the
ALT allele of a heterozygous REF/ALT indel. Replay mirrors this, and such an
event gets `vcf_support_status=heterozygous_genotype_unresolved` (not
supported). Other events on the allele keep their own status.
`vcf_projection.unresolved_genotype_edits` counts every heterozygous edit
replayed under `-H I`: IUPAC SNVs and applied indels. Multi-ALT
heterozygous indels remain unprojectable (`ambiguous_genotype_selection`).

Classifier `allele_confidence` and `exact_match_pct` describe dictionary fit among
candidate windows; neither is a calibrated probability. Use reconstruction status,
`unresolved_regions`, ambiguous-base count and coverage together. Candidate-window
scanning remains the default. The opt-in experimental strict segmentation
mode rejected on development evidence must not be advertised as an accuracy fix.

Strict evaluation distinguishes `sequence_accuracy` (independently recoverable
allele observations) from `literal_sequence_accuracy` (raw assigned-pair equality).
An unsupported second copy of one sequence cannot receive independent recovery
credit. `supported_*` metrics require `exact_sequence_concordance` status;
`missing_alleles` reports literal output cardinality; `independent_missing_alleles`
also excludes unproven duplicates, which receive explicit counts and warnings.
`legacy_supported_*` preserves historical boolean-only comparison. A resolved
reconstruction must have neither missing nor extra alleles. Parse-level CLI errors
must be supplied as nonzero run records by evaluation drivers; only invocations
that enter the pipeline callback can update its execution sidecar.


### Runtime configuration

See the [configuration guide](guides/configuration.md) for the strict JSON schema,
central defaults, command precedence and custom reference layout contract. Both
JSON values and explicit stage options use the same range validation. Standalone
consensus now loads the bundled dictionary by default, or `--repeats-db`/configured
`repeat_dictionary`, and applies the same anchor-aware trimming as the full run.
This corrects its prior fixed-trimming-only behavior. CLI flag names remain valid.

`run_configuration.json` records effective values and hashes of the supplied
configuration, input reads, reference and dictionary before tools execute. A new
invocation removes stale configuration provenance before validating its inputs;
failed input validation can therefore leave only the new execution failure status.
`reference_layout.min_units` and `max_units` control ladder generation only; changing
them alone does not require a different reference for a full run. Changing actual
layout IDs, dictionary or flank length requires an explicit compatible reference.
The evaluator records `fixed_repeat_count` and the generic
`canonical_plus_fixed_matches_reported` diagnostic; the older
`canonical_plus9_matches_reported` field remains a literal legacy diagnostic.

## Mapping and report lifecycle

Mapping retains the streaming `minimap2 -> samtools sort` pipe.
`run.mapping_timeout` / `map --mapping-timeout` / `run --mapping-timeout` set a
finite positive total budget (default 3600 seconds). Tool groups are isolated,
stderr drains concurrently, and failure/interruption/timeout terminates the
launched groups and removes partial BAM/index files. On Linux, an isolated
supervisor adopts and reaps orphan descendants without changing the application
process or reaping its unrelated children. Other POSIX platforms retain process
group cleanup and delegate orphan reaping to the operating system.
A cleanup reserve is inside
the deadline; elapsed indexing time shares the mapping budget. Descendants must
remain in their launched process groups; the operating system cannot guarantee
bounded reaping of a process stuck in uninterruptible kernel I/O. Such cleanup
failures are reported explicitly.

Standalone reports retain `--vcf PATH` and accept repeated
`--allele-vcf LABEL PATH`. Explicit plural input overrides singular input.
Without either option, recorded allele VCF paths are used. Library callers can
supply `vcf_paths`; an explicitly empty mapping also overrides `vcf_path`.
Shared files produce one shared track and do not establish independent haplotypes.
IGV loci come from `summary.alleles[key].contig_name`; unavailable requested
contigs and missing requested files produce visible failures rather than silently
opening another contig. IGV still requires the external `create_report` executable. Real session/browser
verification uses igv-reports 1.13.0. The installed 1.16.0 implementation was
observed concatenating the first VCF record into the column header; reports reject
that malformed output explicitly. Select a working `create_report` on PATH;
this does not change the project dependency lock or variant-calling toolchain.

`report --run-status PATH` supplies explicit execution provenance; otherwise the
report command reads `run_status.json` next to its input summary. It never reads
a sidecar from the report destination. Explicit provenance overrides embedded
summary status, so a failed rerun defeats stale successful artifacts. Malformed
sidecar JSON is an error. Legacy inputs without status retain established
decision behavior with an unavailable-status label. Failed, interrupted,
insufficient and running execution cannot produce a reassuring negative banner.
Recorded mutation evidence is retained with execution warnings.

Clinical gates run before the banner is chosen (`clinical_gates.py`). A mutation
supports PATHOGENIC only when all of these hold:

- it is a frameshift;
- it is an exact dictionary template (`template_match` and `mutation_name`);
- its localization is not ambiguous;
- it has explicit support (exact VCF concordance or `read_support.status=supported`);
- its allele's `depth_status` is not `low`.

NEGATIVE additionally requires:

- resolved allele selection (`selection_status`);
- reported length equal to the consensus contig length (`length`/`reference_length`);
- no unresolved length-partition genotype (`allele_genotype_status`);
- adequate per-allele depth.
- no caller-stage discordance (`stage_concordance.status` is not `discordant_frameshift`).

Summaries without per-allele depth fall back to the 30-read total. The gates can
only lower certainty; PATHOGENIC lists remaining problems as quality caveats.

The full pipeline passes an `analysis_completed` rendering context after all
analysis stages succeed. This means analysis completed and the report is being
generated; it is not terminal execution success. The sidecar remains `running`
while rendering and becomes `completed` only after the callback returns.
Summary status becomes completed after successful rendering (or analysis when
no report was requested). The decorator also synchronizes terminal failure,
interruption and insufficient-evidence status into existing readable summaries,
so moving a summary alone preserves known failure provenance. A requested report
failure propagates to the sidecar and nonzero CLI result, including a missing
optional report dependency.

`exact_match_pct` always uses percent units from 0 through 100: `0.5` means
0.5%, and `1` means 1%. Report text, progress width and ARIA use the same value.
Missing or invalid values are unavailable; historical JSON is never guessed to
contain fractions. Unicode em dashes represent absent nomenclature while Jinja
autoescaping continues to protect supplied strings.
