# Offline fresh-read source assignment analysis

Prepared 2026-09-14. Owned implementation is ignored
`tests/results/production_validation_20260914/fresh_generation/assignment.py` and
`assignment_check.py`. This analysis never changes production code, settings,
caller inputs, or source-truth capture. It is prepared for the coordinator to run
after fresh generation, paired calling, and evaluator reports exist.

## Inputs and denominator

Run once per baseline/candidate arm with its `primary_inventory.json`, the
completed `generation_manifest.json`, and the corresponding evaluator report
`evaluation_baseline.json` or `evaluation_candidate.json` produced by the offline
`analyze.py` task. Reports retain the current evaluator schema with `samples`,
per-sample `alternatives`, assigned truth/prediction IDs, and metrics.

The assignment driver constructs expected IDs from all 32 explicitly declared
generation cases and the four predefined input kinds: main, minority_h1,
minority_h2, and partial_50pct. It requires exactly those 128 primary inventory
rows, each once, with repeat number 1. Timing repeats are excluded. A missing or
duplicated inventory row is a visible error, not permission to reduce the cohort.
Failed generation/calling, absent BAMs, and missing evaluator samples remain rows
with unavailable assignment information and their run/evaluator status.

The generation manifest supplies the exact FASTQ path, SHA256 and record count.
Both the inventory input fingerprint and generation fingerprint must agree.
Main input ordinals come from `main_source_candidates.json`, with source IDs and
haplotypes resolved through `source_truth.json`. Perturbation ordinals come from
the corresponding `{kind}.source_truth.json`. Source files must match the
generation manifest's recorded artifact fingerprints. Every input ordinal must
appear exactly once; source IDs and input-record multiplicity must agree.

`--development` allows the smaller development capture manifest. It is used only
for preparation checks and is not part of the final launch.

## Which BAM is measured

The paths were checked against `calling.py` in both the original baseline source
and current worktree:

- Distinct length branches use `allele_1/allele_reads.bam` and
  `allele_2/allele_reads.bam`, respecting the existing homozygous skip rule.
- The same-length branch uses `merged/allele_reads.bam` once.
- `_extract_and_remap_reads` first writes extracted records to that filename,
  then overwrites it with the final sorted remapped BAM. The retained artifact
  therefore describes the final remapped caller observation, not a surviving copy
  of the pre-remapping extraction. This analysis does not claim to recover that
  overwritten intermediate or reads lost by QNAME-based FASTQ conversion.

The tool boundary is `run_tool_iter(["samtools", "view", "-h", bam])`, using the
existing external-tool abstraction. There is no raw shell subprocess and no
calling or phasing command. BAM hashes are recorded in the output.

## Identity, ambiguity and multiplicity

Matching uses exact sequence plus quality, canonicalized with the reverse
complement and reversed quality string. QNAME is never an identity key. Primary
records with misleading names therefore match the same way as uniquely named
records. Secondary/supplementary alignments and unmapped primary records are
counted separately and do not contribute to mapped assignment support.

Soft-clipped reads retain their full sequence and can match. Any hard clipping,
missing sequence/quality, or unmatched sequence/quality remains explicitly
unavailable. No best approximate alignment, header guess, or source-haplotype
preference fills those gaps.

Input records are grouped by canonical sequence/quality identity. Each matched
BAM record retains **all** compatible input ordinals, source IDs and source
haplotypes. If two identical input records originate from different haplotypes,
neither is arbitrarily assigned to a matched record. The input group preserves
its known total source multiplicities. Record-level output omits sequence and
quality themselves but retains the identity hash and candidate identifiers.

For an input group with n records, nH from haplotype H, and m selected records,
the source-count interval under injective matching within that observation group
is `[max(0, m - (n - nH)), min(m, nH)]`. These bounds are summed across exact
identity groups. They are conditional on records not being reused within a
group; identical records cannot prove distinct physical molecules. If m > n,
the excess is recorded as multiplicity/duplication evidence and the aggregate
source-count and purity bounds are unavailable, rather than forcing a mapping.

Purity divides these bounds by the number of exactly matched mapped primary
records. It is **conditional on the matched subset**; `matching_complete` and
the unmatched/clipped counts must accompany it. PCR descendants are records,
not guaranteed independent original templates.

## Omitted inputs and shared observations

Each input identity group records original ordinals/source counts and the number
of matched mapped records in every observed BAM. Union bounds allow the same
input identity to be present in different allele groups. For per-group counts
m1...mk, omitted-input bounds are `[max(0, n - sum(mi)), n - max(mi)]`.
`cross_group_reuse_minimum` sums unavoidable observations beyond unique input
capacity; it is evidence of shared record use, not necessarily an invalid
biological call. Per-haplotype omission intervals are marginal bounds and need
not attain all endpoints simultaneously.

Omission bounds are unavailable if any expected observation BAM is missing,
unmatched/clipped records prevent complete matching, or within-group multiplicity
exceeds input capacity. Zero observed use is not automatically called an omission
when the observation data themselves are missing.

In the same-length branch, both reconstructed allele labels reference **one
shared observation group**. The merged BAM is not counted twice and does not
provide per-haplotype read assignments. Every purity comparison for such a group
is labeled `shared_observation_not_read_level_phase_assignment`. No WhatsHap
selected-read list is substituted for all-read assignment.

## Linking to reconstruction outcomes

The driver retains every sequence-optimal evaluator alternative. For each
assigned pair it reports truth and prediction IDs, the corresponding source
haplotype, observation-group purity interval, repeat-count error, literal pair
sequence/structure exactness and edit distance. It does not select the assignment
that gives the best purity or mutation interpretation. The evaluator's official
per-sample metrics are copied separately; literal pair exactness must not replace
evidence-aware independent-allele recovery credit.

This permits statements such as an allele assigned to H2 being reconstructed
from an observation group containing only H1 reads, while preserving ambiguity
where identical reads or multiple sequence-optimal assignments prevent such a
conclusion. An observed source group is not a validated full-length haplotype.

## Output and launch interface

The output directory must not exist. It contains `report.json` with exactly 128
sample rows and `records.jsonl` with mapped primary observations, including
unavailable identities. `report.json` records input/evaluation/driver hashes,
source-file hashes, BAM hashes, status counts, the matching scope and assumptions.
`status=observed` means the expected BAM artifacts were read; consult
`matching_complete`, caller run status and evaluator status separately.

From the production-validation worktree, set the coordinator's actual locations
and launch each arm after its evaluator report is complete:

```bash
PV=tests/results/production_validation_20260914
PAIRED_ROOT="$PV/fresh_generation/paired_final"
ANALYSIS_ROOT="$PV/fresh_analysis"
GENERATION_MANIFEST="$PV/fresh_generation/data/generation_manifest.json"
PATH=/home/bernt-popp/miniforge3/envs/env_clair3/bin:$PATH \
  uv run --locked --all-extras python "$PV/fresh_generation/assignment.py" \
  --inventory "$PAIRED_ROOT/baseline/primary_inventory.json" \
  --generation "$GENERATION_MANIFEST" \
  --evaluation "$ANALYSIS_ROOT/evaluation_baseline.json" \
  --output "$ANALYSIS_ROOT/assignment_baseline"
PATH=/home/bernt-popp/miniforge3/envs/env_clair3/bin:$PATH \
  uv run --locked --all-extras python "$PV/fresh_generation/assignment.py" \
  --inventory "$PAIRED_ROOT/candidate/primary_inventory.json" \
  --generation "$GENERATION_MANIFEST" \
  --evaluation "$ANALYSIS_ROOT/evaluation_candidate.json" \
  --output "$ANALYSIS_ROOT/assignment_candidate"
```

The directory variables are launch configuration; the actual paired/analysis
locations must agree with the coordinator's launch. The script itself has no
machine-specific caller or tool path. Source truth stays in the offline analysis.

## Preparation checks actually run

`uv run --locked --all-extras python tests/results/production_validation_20260914/fresh_generation/assignment_check.py`
passes the following checks without simulation or calling:

- Both captured development platforms, seed19200901 only, main plus all three
  perturbations: eight input/source manifests; exact input and source hashes,
  record counts, all-read matching, and zero omitted records verified.
- Deterministic fake SAM: reverse complement with reversed quality, soft clipping,
  duplicate/misleading QNAMEs, mixed-source identical records, excluded secondary/
  supplementary/unmapped records, hard clipping, and unequal quality strings.
- Multiplicity capacities, ambiguous source-count intervals, capacity excess,
  omission intervals and cross-group reuse bounds.
- Same-length observation uses only the actual merged BAM path, preserves one
  observation group for both allele labels, and labels purity as shared.
- Sample orchestration with the real development capture truth and a mocked SAM
  tool boundary; no real BAM generation or caller execution.
- A synthetic 32-case inventory with all 128 inputs failed/missing still produces
  128 unavailable rows; a second invocation refuses its existing output directory.

Ruff checks and formatting are run explicitly with `--no-respect-gitignore`.
These checks validate the offline matching/reporting implementation; they do not
validate final data, performance, or biological reconstruction accuracy.
