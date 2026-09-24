# PRJEB92208 clinical-data baseline

This benchmark measures the frozen MucOneSpan v0.14.1 caller on the deposited
MUC1 Oxford Nanopore reads. It separates execution, publication-reported
positive-control recovery, and independently sourced sequence concordance.
It does not establish clinical sensitivity, specificity, or diagnostic safety.

## Dataset and truth

The refreshed [ENA PRJEB92208 inventory](https://www.ebi.ac.uk/ena/browser/view/PRJEB92208)
contains 20 ONT runs: nine MUC1 amplicons, two MUC1 WGS-labelled inputs, and nine
ACAN exclusions. No PacBio run is available. The 11 selected compressed FASTQs
total **198,156,449 bytes**, unchanged from the supplied September 15 snapshot.
Every selected file was checked against ENA byte counts and MD5, with local
SHA-256 recorded. All 20 accessions remain in the curated manifest.

[Madritsch et al.](https://doi.org/10.1038/s41598-025-30441-3) describe MP1–MP4 as
known dupC positive controls. The paper and supplements do not identify their
orthogonal diagnostic assay or map them to independently characterized patients
in cited work. We report **known-control recovery**, separately from independently
confirmed diagnostic sensitivity. MP1 and MP2 are siblings. MP5 is exploratory;
its reported insertion is a comparator observation without independent confirmation.
HG001, HG003 and HG004 are not assumed to be endpoint-specific negatives.
HG003 and HG004 are the parents of HG002.

Independent full-sequence truth is available for HG002 from
[Q100 v1.1](https://github.com/marbl/HG002), using the versioned paternal and maternal
assemblies. Coding-forward, boundary-inclusive MUC1 intervals contain 3,900 and
4,638 bases. The maternal 18-base biological insertion is preserved. Exact
coordinates, strands, anchors, URLs and checksums are in the truth ledger; the
sequence bases are fetched into a local cache. The two HG002 libraries represent
one biological participant.

Three publication inconsistencies limit interpretation:

- ENA labels the patient WGS run MP1; the paper describes MP4 WGS. Its biological
  identity remains unresolved and receives no patient-specific truth score.
- Supplement 5 contains inconsistent WGS labels. Sequence lengths alone do not
  authorize correcting participant identities.
- A Supplement 1 parameter table prints ACAN boundaries beside MUC1 coordinates.
  We use the documented MUC1 motif 1/9 anchors, independently located in Q100.

Published motif positions and allele assignments remain comparator observations.
No independently mapped event-position or allele-assignment denominator exists.
No endpoint-specific independently confirmed negative set was established, so
**specificity is not estimable**. Historical simulation performance is separate
from these clinical-data observations.

## Prespecified processing and endpoints

Preparation validates four-line FASTQ records and decompresses without changing
bases, qualities, IDs, orientation or read selection. No outcome-dependent
filtering or cropping is applied. Amplicon anchor-spanning reads typically retain
only 32–33 upstream and 71–72 downstream bases. Short fragments are retained.
All nine amplicon libraries contain at least 1,022 exact ordered full-anchor spans.

The two WGS-labelled files are analyzed through the existing ONT FASTQ input path
as a separate exploratory arm. HG002 has only 54 deposited reads, with median
length 105,351 bases; the patient-labelled file contains 2,152 reads with median
930 bases and five exact full-anchor spans. These are the deposited inputs, not
complete genomes. Their geometry and small sizes support prior selection; the
original deposition-filter command is unavailable. No end-to-end genome-processing
throughput is claimed.

Execution success requires a completed CLI invocation and current validated
artifacts. **Artifact callability** additionally requires the existing evaluator's
completed observation state. This permits reporting an annotation; it is not a
clinical-quality certification. Ambiguous reconstruction, tool failures, no-calls
and unattempted runs remain in all-sample denominators. Source and artifact hashes
are checked before replay; stale output cannot override a failed invocation.

The event endpoint is the exact named event `dupC`. A different positive alarm is
not recovery. Existing sequence/VCF-supported recovery is reported separately.
Positions and allele assignments are not scored against comparator-only truth.

Sequence comparison requires unique full motif 1 and motif 9 anchors in the
predeclared orientation. It retains every interior base, accepts allele swaps,
and reports optimal sequence assignments, exact equality, edit distance and
length differences. Literal sequence agreement is separate from independently
recoverable haplotype observations. Two copies of one unsupported sequence cannot
prove reconstruction of both haplotypes. The strict `sequence_equal` and
`independent_exact_alleles` fields additionally require **whole-run callability**;
literal equality does not. Thus an ambiguous annotation can remove strict recovery
credit even when the DNA sequences are exact. A length divided by 60 is not a reliable
repeat count when biological motif indels are present.

## Final measured results

The final verification cohort used the reviewed immutable harness. All eleven
records passed input, environment and output-hash validation. Its scored results
exactly match the initial engineering cohort, and all **306 scientific artifact
hashes** agree after documented path/VCF-header normalization. Only the final
cohort contributes to the resource figures below.

| Endpoint | Final result |
| --- | --- |
| Execution completed | **11/11**: amplicon 9/9; WGS 2/2 |
| Artifact callable | **6/11 (54.5%)**: amplicon 5/9; WGS 1/2 |
| Ambiguous reconstruction | **5/11**: four amplicons, one WGS |
| Insufficient-evidence exit, execution failure, invalid artifacts, unsupported, unattempted | **0/11 each** |
| Reported-control exact `dupC` recovery | **0/4** runs and **0/4** participants |
| Reported-control sequence/VCF-supported `dupC` recovery | **0/4** |
| Callable-only reported-control recovery | **0/1**; three positive controls are uncallable |
| Independently confirmed diagnostic sensitivity / specificity | **Not estimable**: no eligible positive / negative event set |
| Independently verified event position / allele assignment | **Not estimable** |

### Per-library execution and resources

All rows completed with CLI exit 0 and generated a report. “Ambiguous” is the
scientific observation state, not an execution failure. RSS is the largest
individual-process high-water value converted from KiB to MiB; see its limits below.

| Run | Library | Reads retained | Artifact state | Wall (s) | RSS lower bound (MiB) |
| --- | --- | ---: | --- | ---: | ---: |
| ERR15277562 | HG001 PCR | 10,197 | Callable | 84.84 | 567.3 |
| ERR15277563 | HG002 PCR | 13,733 | Ambiguous | 147.81 | 763.6 |
| ERR15277564 | HG003 PCR | 12,485 | Callable | 101.31 | 544.3 |
| ERR15277565 | HG004 PCR | 45,325 | Callable | 128.86 | 634.2 |
| ERR15277566 | MP1 PCR | 5,941 | Callable | 51.21 | 441.4 |
| ERR15277567 | MP2 PCR | 5,332 | Ambiguous | 64.65 | 671.1 |
| ERR15277568 | MP3 PCR | 26,672 | Ambiguous | 143.37 | 804.3 |
| ERR15277569 | MP4 PCR | 31,808 | Ambiguous | 68.54 | 498.0 |
| ERR15277570 | MP5 PCR | 8,307 | Callable | 44.90 | 423.9 |
| ERR15277552 | HG002 WGS | 54 | Ambiguous | 6.75 | 412.1 |
| ERR15277553 | Patient WGS; identity unresolved | 2,152 | Callable | 5.13 | 415.1 |

### Biological participants and event-truth eligibility

There are **nine resolved participant identities**, plus one library with unresolved
patient identity. HG002 has two libraries. The unresolved WGS library is not counted
as an independent tenth participant or reassigned to MP1/MP4.

| Participant | Libraries | Relationship | `dupC` truth category | Callable libraries |
| --- | ---: | --- | --- | ---: |
| HG001 | 1 | None documented | Unknown | 1/1 |
| HG002 | 2 | child_of:HG003, child_of:HG004 | Unknown | 0/2 |
| HG003 | 1 | father_of:HG002 | Unknown | 1/1 |
| HG004 | 1 | mother_of:HG002 | Unknown | 1/1 |
| MP1 | 1 | sibling_of:MP2 | Reported positive control | 1/1 |
| MP2 | 1 | sibling_of:MP1 | Reported positive control | 0/1 |
| MP3 | 1 | None documented | Reported positive control | 0/1 |
| MP4 | 1 | None documented | Reported positive control | 0/1 |
| MP5 | 1 | None documented | Exploratory comparator only | 1/1 |
| Unresolved patient WGS | 1 | None documented | Unknown | 1/1 |

### Independent HG002 sequence comparison

| Library | Observed lengths (bp) | Literal exact alleles | Total edit distance | Strict callable exact pair |
| --- | --- | ---: | ---: | --- |
| PCR | 3,900 / 900 | **1/2** | 3,738 | No |
| WGS | 3,900 / 4,638 | **2/2** | **0** | No: ambiguous interpretation |

Across the two libraries of this one participant, literal pair agreement is
**1/2**, literal allele agreement **3/4**, and strict callable sequence recovery
**0/2 pairs, 0/4 alleles**. The strict zero reflects the whole-run evidence gate;
it does not mean the exact WGS sequences were wrong. Both libraries are ambiguous.

For PCR, sequence-optimal assignment gives length errors **0 / −3,738 bp**.
The separate `sequence_length_error_bp` field compares independently sorted length
lists and gives **−3,000 / −738 bp**. These are different pairings; the per-assignment
field is the appropriate one when discussing which Q100 haplotype was reconstructed.
WGS has zero errors under both conventions.

## Subsequent MP1 classifier correction

The frozen results above describe the original caller, before the correction for
[#52](https://github.com/berntpopp/MucOneSpan/issues/52). The current change expands
exact `dupC` templates to repeat backgrounds with the same terminal seven-C tract.
It retains the actual parent repeat and consumes the correct 61-base unit.

Isolated reclassification of all 22 frozen allele consensuses changes exactly one
mutation annotation list: MP1 allele 2 now reports **B:dupC at repeat 17**, with
`exact_sequence_concordance` VCF support, and the spurious downstream A insertion
is absent. The other 21 mutation annotation lists are unchanged. This replay does
not substitute for a newly measured complete cohort; the original 0/4 baseline
remains intact. A separate full MP1 pipeline rerun completes with a callable observation and
exactly this supported `B:dupC` event. Its provenance is recorded in
`mp1-classifier-fix.json`. This is one corrected run, not a full post-fix cohort.

## What the discrepancies mean

These stage observations explain why a successful process or a positive clinical
banner is insufficient to claim exact mutation recovery:

- **Allele selection:** HG002 PCR reconstructs the 3,900-base paternal haplotype,
  but its other candidate is only 900 bases instead of the independent maternal
  4,638 bases. HG004 also has an observed long anchor-spanning read population
  that is merged into a broad candidate cluster; a very short second candidate
  is selected. The HG004 observation is input evidence, not complete haplotype truth.
- **Calling and consensus:** MP1 and MP3 have raw filtered-VCF C insertions in
  seven-C tracts, compatible with the paper's dupC descriptions. MP1's final
  classifier instead reports an A insertion without named dupC support. MP3's
  heterozygous insertion is on the alternate genotype allele, while GT1 consensus
  selects the reference allele. These failures cannot all be described as an
  inability of Clair3 to detect any insertion.
- **Interpretation and display:** HG002 WGS reconstructs both Q100 sequences
  exactly. Its maternal 18-base in-frame expansion is nevertheless interpreted
  as an ambiguous boundary mutation and produces a PATHOGENIC banner. Exact
  sequence reconstruction and correct pathogenic interpretation are separate tasks.
- **Evidence strength:** large alignment-record counts include repeated alignments;
  they are not independent molecules. Complete dictionary segmentation and zero
  ambiguous bases coexist with unverified reference confidence. Empty VCF regions
  can remain reference-derived in consensus.
- **Unknown truth:** MP5's supported generic C insertion is exploratory, and the
  patient WGS identity is unresolved. Neither can repair the four-control denominator.

The following accuracy work should preserve these cases and distinguish allele
loss, genotype selection, indel representation, motif classification and display.
For [#20](https://github.com/berntpopp/MucOneSpan/issues/20) and
[#47](https://github.com/berntpopp/MucOneSpan/issues/47), investigate the current
CIGAR-indel refinement and physical span evidence against Q100; historical issue
assumptions need rechecking. For [#21](https://github.com/berntpopp/MucOneSpan/issues/21),
trace the existing MP1/MP3 insertion evidence through consensus and classification
before changing calling thresholds. For
[#46](https://github.com/berntpopp/MucOneSpan/issues/46), quantify unique-read assignment,
competing haplotypes and reference-only regions. These observations do not select
POA, a new threshold, or a clinical tier change.

## v0.16.0 clinical-safety rerun (PRJEB92208)

The cohort was rerun with the Phase 0 clinical-safety changes (caller commit
`90de0b4`). The decision is `compute_clinical_decision` applied to each
`summary.json`; gate reasons are truncated to 300 characters.

| Library | Sample | Decision | Events (allele:name@repeat:support) | Gate reasons |
| --- | --- | --- | --- | --- |
| ERR15277552 | HG002 WGS | INCONCLUSIVE | allele_2:7@75:localization_ambiguous | Allele 2: Observed sequence variant (Unknown variant at repeat 75) is inconclusive (frameshift not established; event identity not established (no exact dictionary template); localization ambiguous; no explicit sequence-level support (localization_ambiguous); carrying allele is below the per-allele  |
| ERR15277553 | WGS (identity unresolved) | INCONCLUSIVE | - | Allele 1: allele selection unresolved (unresolved_unselected_clusters; secondary mode fraction 0.0).; Allele 1: 17 primary alignments, below the per-allele depth gate (30).; Allele 2: allele selection unresolved (unresolved_secondary_mode; secondary mode fraction 0.25).; Allele 2: 7 primary alignmen |
| ERR15277562 | HG001 | NO_PATHOGENIC_VARIANT_DETECTED | - | Allele 1: 40 repeats (31 canonical units) - 46915 reads; Allele 2: 72 repeats (63 canonical units) - 10166 reads |
| ERR15277563 | HG002 PCR | INCONCLUSIVE | allele_2:7@75:localization_ambiguous | Allele 2: Observed sequence variant (Unknown variant at repeat 75) is inconclusive (frameshift not established; event identity not established (no exact dictionary template); localization ambiguous; no explicit sequence-level support (localization_ambiguous)). |
| ERR15277564 | HG003 | NO_PATHOGENIC_VARIANT_DETECTED | - | Allele 1: 44 repeats (35 canonical units) - 47980 reads; Allele 2: 65 repeats (56 canonical units) - 14680 reads |
| ERR15277565 | HG004 | INCONCLUSIVE | allele_2:7@75:localization_ambiguous | Allele 2: Observed sequence variant (Unknown variant at repeat 75) is inconclusive (frameshift not established; event identity not established (no exact dictionary template); localization ambiguous; no explicit sequence-level support (localization_ambiguous)). |
| ERR15277566 | MP1 | PATHOGENIC | allele_2:dupC@17:exact_sequence_concordance | Allele 2: dupC at repeat unit 17; Quality caveat: Allele 1: reported length 39 differs from the consensus contig length 44. |
| ERR15277567 | MP2 | PATHOGENIC | allele_2:dupC@17:exact_sequence_concordance | Allele 2: dupC at repeat unit 17 |
| ERR15277568 | MP3 | INCONCLUSIVE | allele_1:dupC@7:heterozygous_genotype_unresolved; allele_2:X@35:localization_ambiguous; allele_2:dupC@80:localization_ambiguous | Allele 1: Observed sequence variant (dupC at repeat 7) is inconclusive (heterozygous genotype unresolved within its allele).; Allele 2: Observed sequence variant (Unknown variant at repeat 35) is inconclusive (event identity not established (no exact dictionary template); localization ambiguous; no  |
| ERR15277569 | MP4 | PATHOGENIC | allele_2:dupC@49:exact_sequence_concordance | Allele 2: dupC at repeat unit 49 |
| ERR15277570 | MP5 | INCONCLUSIVE | allele_2:X@35:exact_sequence_concordance | Allele 2: Observed sequence variant (Unknown variant at repeat 35) is inconclusive (event identity not established (no exact dictionary template)).; Allele 1: allele selection unresolved (unresolved_secondary_mode; secondary mode fraction 0.276).; Allele 1: reported length 69 differs from the consen |

The positives are publication-reported known controls from three families; this
rerun claims no sensitivity or specificity. Compared with cohort-v3, only MP4
changes at the scorer level (dupC recovered on allele 2); no control gains an event.

Frozen tools: minimap2 2.28-r1209, samtools 1.15.1, bcftools 1.17 and Clair3
1.0.10 with model `r1041_e82_400bps_sup_v500`, Python 3.12.9.

## Reproduction

The additive command is `scripts/clinical_benchmark.py`. Use Python 3.10 or newer,
uv, minimap2, samtools, bcftools and Clair3 v1.x with a compatible TensorFlow model.
The measured caller used **Python 3.12.9**; use that interpreter for close environment
reproduction. The supported Python range does not imply identical measurements.
The frozen tools were minimap2 2.28-r1209, samtools 1.15.1, bcftools 1.17 and
Clair3 1.0.10. Locked Python packages and model/resource hashes accompany the
benchmark provenance. The clinical model is `r1041_e82_400bps_sup_v500`, selected
before examining predictions from Dorado 0.7.1's contemporary SUP v5.0.0 default.
The study does not disclose the exact original model or sampling rate; this is an
explicit compatibility assumption, not a measured model comparison.

Choose data/output directories outside the checkout and put external tool
executables on `PATH`. Set `CUDA_VISIBLE_DEVICES` empty for the measured CPU path.
The caller remains at commit `c08000a0b97e9ee00e560a369ef320b4d06d3004`;
benchmark helpers can be overlaid from a pinned draft revision without changing
its science. Record `BENCHMARK_REF` with the resulting environment:

```bash
git fetch origin feat/clinical-benchmark
BENCHMARK_REF="$(git rev-parse origin/feat/clinical-benchmark)"
git worktree add --detach ../MucOneSpan-clinical-baseline c08000a0b97e9ee00e560a369ef320b4d06d3004
cd ../MucOneSpan-clinical-baseline
git restore --source "$BENCHMARK_REF" --worktree -- scripts/clinical_benchmark.py \
  src/muc_one_span/clinical_data.py src/muc_one_span/clinical_provenance.py \
  src/muc_one_span/clinical_runner.py src/muc_one_span/clinical_worker.py \
  src/muc_one_span/clinical_scoring.py src/muc_one_span/clinical_truth.py \
  benchmarks/clinical/prjeb92208
uv venv --python 3.12.9
make dev
export CUDA_VISIBLE_DEVICES=''
export DATA_ROOT="$PWD/../clinical-data"
export RESULT_ROOT="$PWD/../clinical-results"
export MODEL_ROOT="$DATA_ROOT/models/r1041_e82_400bps_sup_v500"
mkdir -p "$DATA_ROOT/models" "$RESULT_ROOT"
curl --fail --location --output "$DATA_ROOT/models/model.tar.gz" \
  https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v500.tar.gz
printf '%s  %s\n' 01c05768661bdd7de611e6bae1043c43b7523a54b223e029c683bfac0db7a678 \
  "$DATA_ROOT/models/model.tar.gz" | sha256sum --check
tar -xzf "$DATA_ROOT/models/model.tar.gz" -C "$DATA_ROOT/models"
uv run --locked --all-extras python scripts/clinical_benchmark.py inventory \
  --output "$DATA_ROOT/refreshed-manifest.json"
# For exact accession/checksum replay, use the curated snapshot below.
uv run --locked --all-extras python scripts/clinical_benchmark.py prepare \
  --manifest benchmarks/clinical/prjeb92208/manifest.json --data-root "$DATA_ROOT/reads"
uv run --locked --all-extras python scripts/clinical_benchmark.py truth \
  --ledger benchmarks/clinical/prjeb92208/truth-ledger.json \
  --cache-root "$DATA_ROOT/truth-cache" --output "$DATA_ROOT/truth.json"
uv run --locked --all-extras python scripts/clinical_benchmark.py freeze \
  --checkout . --model "$MODEL_ROOT" --output "$DATA_ROOT/environment.json"
uv run --locked --all-extras python scripts/clinical_benchmark.py run \
  --manifest benchmarks/clinical/prjeb92208/manifest.json \
  --data-root "$DATA_ROOT/reads" --output-root "$RESULT_ROOT" \
  --environment "$DATA_ROOT/environment.json" --threads 2 --timeout 1800
uv run --locked --all-extras python scripts/clinical_benchmark.py score \
  --manifest benchmarks/clinical/prjeb92208/manifest.json --truth "$DATA_ROOT/truth.json" \
  --output-root "$RESULT_ROOT" --output "$RESULT_ROOT/results.json"
```

The original run used harness commit `1deb5dcfd1dbeede0fbd3ca859cc9a6aab011d47`.
The overlay above records the baseline checkout as `harness_commit`; record its
source commit `BENCHMARK_REF` separately and compare helper hashes with the published
environment summary. Tool packages can be recreated from `tool-environment.txt`
with `conda create --name muconespan-clinical-tools --file benchmarks/clinical/prjeb92208/tool-environment.txt`
on Linux x86-64, then activated before freezing. Hardware, operating system and
shared-workstation scheduling still prevent identical runtime measurements.

A nonzero run command retains failures; execute scoring explicitly afterward.
`--run ERR...` selects a wiring smoke case; use a separate smoke output root.
`--resume` reuses only a provenance-validated identical attempt. Changed inputs,
settings or outputs require a new output root. Downloads use verified temporary
files and safely resume or restart partial content. Synthetic unit tests exercise
offline replay, corruption, duplicate identities and source/output collisions.

## Execution and resource interpretation

The final cohort uses one observed execution per library, serially, with two
requested threads and a 1,800-second supervised lifecycle budget. Smoke execution
is excluded from resource aggregates. Report generation is timed separately;
IGV remains at its existing `off` default. The benchmark does not change caller
thresholds, repeat dictionaries, clinical tiers or experimental settings.

Peak memory is the maximum `ru_maxrss` of the isolated worker and its reaped
children, measured in KiB on Linux. This is the largest individual-process
high-water RSS: **a lower bound on concurrent process-tree peak**, not a summed
peak. Timings are single observations on a shared workstation, not medians or a
controlled hardware comparison. Download/preprocessing time is separate.

The measured host was an AMD Ryzen 9 9950X (16 cores / 32 logical CPUs),
with about 59.4 GiB RAM, Linux x86-64, CPU inference and two requested threads.
The filesystem was warm from preliminary execution; other workstation activity
was not controlled. Each library was measured once in the final cohort.

| Arm | Summed wall (s) | Per-library wall range (s) | Timed analysis stages (s) | Report generation (s) | Preparation incl. download (s) |
| --- | ---: | --- | ---: | ---: | ---: |
| 9 amplicons | 835.50 | 44.90–147.81 | 833.56 | 0.199 | 24.57 |
| 2 WGS-labelled inputs | 11.88 | 5.13–6.75 | 11.50 | 0.044 | 1.33 |

Summed final supervised wall time is **847.38 seconds (14.12 minutes)**.
Stage sums exclude startup and supervisory overhead. Preparation measures the
combined acquisition, validation and decompression operation; it does not isolate
network transfer from local work. WGS timings concern only the small deposited
read sets, not whole-genome throughput.

No quantitative runtime, peak-memory or hardware benchmark was found in the
paper, its six supplements, or the pinned VNTRPipeline v1.0/v2.0 README files.
Lab incubation durations are not computational runtimes. A matched speedup or
memory comparison is therefore unavailable.

## Scope and limitations

This small, partly related cohort cannot estimate population-wide clinical
performance. High alignment-record counts are not independent molecule counts.
Completed dictionary classification and zero ambiguous bases can coexist with
unverified reference-base confidence. Neither an empty event list nor a generated
clinical display establishes independently verified event absence.

Raw reads, full participant sequences and generated pipeline outputs stay outside
Git. Only accession metadata, source/sequence provenance and compact aggregate
benchmark records are published. The strict simulation evaluator remains unchanged.

## Issue #44 acceptance and remaining scope

| Requested item | Delivered evidence | Remaining limitation |
| --- | --- | --- |
| Reproducible download/preprocessing | Versioned inventory, checksum validation, deterministic FASTQ preparation and corruption/replay tests | ENA metadata and source URLs may change; curated checksums freeze this snapshot |
| All available PacBio/ONT samples | All nine MUC1 ONT amplicons and both WGS-labelled inputs executed | PRJEB92208 contains no PacBio runs; WGS inputs are selected deposited reads |
| Orthogonal sensitivity/specificity | Separate reported-control recovery and independent HG002 Q100 sequence comparison | No sample-linked orthogonal diagnostic assay or independently confirmed negative set established |
| Runtime and peak-memory comparison | Per-run stages, wall time and explicitly bounded RSS measurement | No published computational measurements available for a matched comparison |
| Documented results | Per-run and participant accounting, provenance, commands and stage-specific errors | This is a frozen baseline with observed failures, not diagnostic validation |

[Issue #44](https://github.com/berntpopp/MucOneSpan/issues/44) remains open for the
unsupported original scope. The benchmark release records an available-data
baseline; it does not resolve the observed ONT accuracy defects.

## Wave 2 clinical improvements and comparator evaluation

Following the baseline measurements above, focused fixes resolved five identified accuracy and interpretation defects:

- [#52](https://github.com/berntpopp/MucOneSpan/issues/52): B-repeat `dupC` template expansion correctly identifies frameshift variants across compatible terminal seven-C tracts without consuming downstream repeats.
- [#53](https://github.com/berntpopp/MucOneSpan/issues/53): Unphased distinct-length candidates receive IUPAC consensus (`"I"`) instead of forcing alternate GT1.
- [#54](https://github.com/berntpopp/MucOneSpan/issues/54): Normalized valley splitting avoids fragment bias, recovering the maternal 4,638 bp allele in HG002 PCR and long alleles in MP2 and MP4.
- [#55](https://github.com/berntpopp/MucOneSpan/issues/55): Calibrated clinical decision calling reports `INCONCLUSIVE` for benign in-frame expansions (HG002 WGS 18 bp insertion) and unphased/unverified reconstructions (HG002 PCR, MP3), and reports `PATHOGENIC` only when frameshift, localization, and VCF support are established.
- [#56](https://github.com/berntpopp/MucOneSpan/issues/56): VNTRPipeline v1.0 executed under isolated container runtime (`ghcr.io/dhmeduni/vntr_pipeline`) with offline BiocManager runtime dependency mounted.

### Measured cohort-v3 vs baseline results

| Metric | Baseline (v1) | Wave 2 (v3) |
| --- | --- | --- |
| Execution completed | 11/11 (100%) | 11/11 (100%) |
| Artifact callable | 6/11 (54.5%) | 7/11 (63.6%) |
| Amplicon callability | 5/9 (55.6%) | 6/9 (66.7%) |
| Reported-control exact dupC recovery | 0/4 (0%) | 2/4 (50%) (MP1, MP2) |
| Supported positive recovery | 0/4 (0%) | 2/4 (50%) (both exact sequence concordance) |
| Callable reported-control recovery | 0/1 (0%) | 2/3 (66.7%) |
| HG002 PCR sequence recovery (Q100) | 1/2 alleles (edit dist 3,738 bp) | 2/2 alleles (edit dist 0 bp) |
| HG002 WGS sequence recovery (Q100) | 2/2 alleles (edit dist 0 bp) | 2/2 alleles (edit dist 0 bp) |
| HG002 clinical decision | PATHOGENIC (false alarm) | INCONCLUSIVE (calibrated in-frame finding) |
| MP3 clinical decision | PATHOGENIC (unphased forced GT1) | INCONCLUSIVE (calibrated unphased IUPAC) |
| MP4 clinical decision | Ambiguous | NO_PATHOGENIC_VARIANT_DETECTED |

### VNTRPipeline v1.0 comparator findings

Under containerized execution with normalized headers and offline BiocManager 1.30.27:
- **MP1 (ERR15277566):** VNTRPipeline identified LoF at motif position 17 with sequence `GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCCCCCCCCA` in 452s. MucOneSpan v3 identified the exact identical motif, repeat position (17), and 61-base sequence with `exact_sequence_concordance`.
- **MP2 (ERR15277567):** VNTRPipeline assembled haplotypes of 4,681 bp (78 repeats) and 3,000 bp (50 repeats) in 390s, identifying LoF at motif position 17 with the identical 61-base sequence. MucOneSpan v3 selected alleles 50 and 78 and called `PATHOGENIC` with `exact_sequence_concordance` at repeat 17.
- **MP4 (ERR15277569):** VNTRPipeline assembled haplotypes of 4,801 bp (80 repeats) and 2,700 bp (45 repeats) in 422s, identifying LoF at motif position 49 with sequence `GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCCA`. MucOneSpan v3 selected allele lengths 49 and 80; Clair3 did not report an insertion call passing variant quality filters, leading to `NO_PATHOGENIC_VARIANT_DETECTED`.
- **HG002 PCR (ERR15277563):** VNTRPipeline assembled haplotypes of 4,638 bp and 3,900 bp in 635s matching independent Q100 truth to the exact base pair, with 0 LoF variants. MucOneSpan v3 similarly reconstructed 4,638 bp and 3,900 bp with 0 edit distance to Q100.
- **Computational performance:** MucOneSpan completes individual amplicon libraries in 45–150 seconds on two CPU threads (~14 minutes for all 11 libraries combined), while VNTRPipeline requires 390–635 seconds per library utilizing 32 minimap2 threads and Canu assembly.

