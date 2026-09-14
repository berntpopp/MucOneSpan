# Mapping, allele assignment, and calling audit — 2026-09-14

Scope: current source plus 44 retained HiFi mapping BAMs from `tests/results/sample_*` and their MucOneUp truth in `tests/data/generated/*/*.001.simulation_stats.json`. These BAMs are historical, not newly generated reads or newly run alignments. Source was rereun against the retained BAMs. Full-pipeline reruns under `tests/results/deep_validation_20260914/hifi` are being evaluated separately. No production edits were made.

## Main conclusions

1. Retain secondary alignments as candidate-length evidence, but distinguish them from independent read support. Removing secondary alignments naively destroyed exact allele-length accuracy in the 44-sample ablation.
2. Exact terminal-repeat anchors offer a strong exploratory alternative for molecule length and read assignment on these HiFi simulations. This is promising evidence, not a validated replacement algorithm.
3. The historical claim that long dupC alleles produce no Clair3 variants is contradicted by the retained raw VCFs. One missing dupC is present at QUAL 3.22 and is removed by the default QUAL 5 filter.
4. Same-length sequence reconstruction lacks haplotype phasing and depends on a false assumption about default consensus behavior in the available bcftools version.
5. Prioritize correctness and evidence accounting before more expensive caller/model changes. Streaming remapping and coherent thread budgeting are independent performance candidates.

## Verified findings

### Alignment records and independent read support are different

`mapping.py:104` invokes minimap2 without disabling secondary alignments. `mapping.py:180` obtains raw idxstats. `alleles.py:67` reads its mapped-record column; `_find_clusters` at `alleles.py:197` thresholds that count, then sums it as `reads`. Both refinement (`alleles.py:119`) and valley splitting (`alleles.py:256`) read unfiltered SAM records. `calling.py:41` extracts the selected contigs and `calling.py:100` runs `samtools fastq`, which excludes secondary/supplementary alignments in this environment.

Measured on historical artifacts, after rerunning current detection:

| Sample | Ladder records / primary mapped | Detected lengths | Reported allele support | Primary reads in those clusters |
|---|---:|---|---|---|
| 60/60 | 1,014 / 169 | 60/63 | 723 / 276 | 164 / 2 |
| 25/140 | 1,044 / 174 | 25/28 | 734 / 278 | 168 / 1 |
| 100/120 dupC | 1,043 / 174 | 100/120 | 657 / 340 | 111 / 56 |
| 120/140 dupC | 1,062 / 177 | 120/140 | 654 / 387 | 110 / 64 |

Thus the false second allele in 60/60 is fed only two primary reads despite a reported 276 reads. The 25/140 data contain five primary long reads, spread over contigs 129–131; its spurious second short allele is fed one primary read. These observations establish a coverage-accounting failure. They do not establish that secondary alignments should simply be removed from length inference.

### Primary-only ablation rejects a tempting simple fix

Script: `tests/results/deep_validation_20260914/mapping_audit/ablate.py`.
Output: `tests/results/deep_validation_20260914/mapping_audit/ablation.json`.

All 44 historical samples with matched stats truth were processed under three methods, with unchanged `min_coverage=10` and current code. Primary means `samtools view -F 2308` (unmapped, secondary, supplementary removed). The third method uses a temporary primary-only indexed BAM for both histogram and AS/indel refinement; temporary BAMs were removed.

| Method | Called | Exact length pairs | Both alleles within ±2 | Mean absolute error per allele |
|---|---:|---:|---:|---:|
| Current all-alignments | 44/44 | 40/44 | 40/44 | 1.398 |
| Primary histogram, all-alignments refinement | 44/44 | 0/44 | 38/44 | 3.023 |
| Primary histogram and primary refinement | 44/44 | 0/44 | 39/44 | 3.034 |

Primary-only inference commonly shifts both lengths down 1–2 units (60/80 becomes 58/78). It improves 60/60 to 59/59 but loses the 40/44 distinction and still misses the asymmetric long allele. Correct candidate contigs frequently receive secondary rather than primary alignments. A sound redesign should separate candidate fit from independent molecule support rather than discard all alternative alignments.

Truth was read from stats, not filenames: e.g. `sample_bench_5000` is 52/63, `sample_bench_5001` is 50/65, and `sample_bench_5004` is 38/47. Filename-based benchmark assumptions would silently mis-score these.

### Exploratory terminal-anchor estimator

Script: `tests/results/deep_validation_20260914/mapping_audit/anchors.py`.
Output: `tests/results/deep_validation_20260914/mapping_audit/anchors.json`.

For each historical primary alignment sequence in reference orientation, find the exact first 20 bases of dictionary repeat `1` and the exact last 20 bases of repeat `9`, then divide the inclusive span by 60 and round. Truth is used only to score the resulting lengths. No ladder contig position, AS score, or truth-derived interval enters the estimate.

- 6,238 of 7,251 primary reads (86.0%) contain both anchors in order.
- 6,222 of those 6,238 lengths (99.74%) round within one repeat of the nearer truth allele.
- Both truth counts have at least three reads rounding exactly to them in 44/44 samples; at least ten in 43/44.
- 50/53: 79 reads at 50 and 74 at 53; 80/83: 66 at 80 and 57 at 83.
- 60/60: 137 at 60, seven at 59; 25/140: 153 at 25, four at 140, one at 24.

These are per-read/histogram observations, not 44/44 validated genotype calls. An unsupervised peak criterion, uncertainty, minority-read minimum, and handling of anchor failures remain to be specified. Exact anchors can introduce ascertainment bias and may fail more often for ONT, mutations in anchor repeats, clipped reads, or nonstandard flanks. The asymmetric long allele still lacks ten independent reads. Do not conceal that uncertainty by lowering a threshold solely for this sample.

### Low coverage is thresholded per contig before pooling

`alleles.py:197` removes every contig below `min_coverage` before clustering. Deterministic probes:

- `detect_alleles({50:9,51:9,52:9},10)` raises despite 27 pooled records.
- `{51:100,110:9,111:9,112:9}` yields duplicated 60/60 calls, losing the diffuse minor cluster.

This can reduce sensitivity for an allele with adequate total support spread across neighboring contigs. Documentation describing the threshold as per allele is misleading. The thresholded contig list also controls downstream extraction, so excluded fringe reads never reach Clair3.

### Flat indel profiles can split arbitrarily

`alleles.py:290` accepts non-strict minima, so every member of a plateau is a valley. The two smallest entries are selected before separation is tested (`alleles.py:298`). Exact zero-indel SAM probes with equal support:

- contigs 40,42,45: no split.
- contigs 40,43,46: splits into 40 and pooled 43/46.

There is no valley prominence requirement. The result changes solely because the first two plateau points are separated by three instead of two. Existing `test_no_valleys_returns_none` uses 40,42,45,48,50, so its passing result does not demonstrate general plateau safety. Near-length alleles one or two repeats apart are explicitly excluded by the minimum separation of three (`alleles.py:306`).

### Read IDs in retained simulator outputs are not unique

Input BAM primary records already contain duplicated `S/N/ccs` names across distinct reads:

| Sample | Primary records | Unique names | Names occurring more than once |
|---|---:|---:|---:|
| 60/60 | 169 | 100 | 69 |
| 100/120 dupC | 174 | 123 | 51 |
| 25/140 | 174 | 169 | 5 |

For 100/120, `S/15/ccs` labels both 7,242- and 6,049-base input reads. This is an upstream simulation artifact and makes unique-QNAME counts invalid molecule counts here. Proposed read-name rescue/phasing must first provide collision-free IDs, keeping an explicit original-ID mapping. Current SAM records must not simply be deduplicated by QNAME. Effects on external phasing and FASTQ conversion are plausible but not isolated in this audit.

### Historical long dupC miss is partly filtering, not absence of a candidate

Retained raw VCF record counts:

| Sample | allele_1 raw records | allele_2 raw records | Relevant indel |
|---|---:|---:|---|
| 60/60 | 61 | 45 | No indels in either raw VCF |
| 100/120 dupC | 85 | 193 | contig_91:3192 G→GC, QUAL 3.22, PASS |
| 120/140 dupC | 183 | 118 | No indels in either raw VCF |

The 100/120 insertion is consistent with the targeted repeat-45 dupC position. `calling.py:318` defaults to QUAL 5; `vcf.py:76` applies `QUAL>=5`, so this candidate is removed. This disproves the broad historical claim of zero variants. It does not show that globally lowering QUAL improves specificity. The 120/140 raw-indel absence remains unresolved: inspect truth-anchored read support, alignment representation, candidate generation, pileup/full-alignment outputs, and model behavior before declaring an inherent Clair3 length limitation.

`vcf.py:28` documents that `min_dp` is ignored, despite a default `min_dp=5` in callers. Low-depth PASS variants can therefore survive: retained 60/60 false second allele contains QUAL ≈5 SNPs at DP=2. Their downstream classification impact needs separate evaluation.

### Same-length reconstruction cannot represent general haplotypes

`calling.py:251` recognizes heterozygous genotypes, then creates a homozygous-alt-only VCF and a VCF containing all variants (`calling.py:268`, `calling.py:287`). Phasing is not preserved as two haplotypes; opposite-phase variants are assigned together. `consensus.py:33` supplies neither sample nor haplotype options.

Real bcftools 1.17 probe (remapping, Clair3, and the normalization wrapper mocked; VCF selection, indexing, genotype query, and consensus real): ten-base `AAAAAAAAAA` reference, A→C at position 2 with GT `1|0`, A→G at position 8 with GT `0|1`. Expected haplotypes are `ACAAAAAAAA` and `AAAAAAAGAA`. Current functions produce:

- allele_1: `AAAAAAAAAA`
- allele_2: `AMAAAAARAA`

Thus available bcftools emits IUPAC ambiguity here, contradicting the code comment that default consensus necessarily produces the intended mutant ALT sequence. Even explicit all-ALT selection alone would create a biologically incorrect cis haplotype in this example. Reconstruct phased haplotypes or report ambiguity; do not claim arbitrary het sites form one mutant allele.

`calling.py:260` declares homozygous when no heterozygous calls remain. Absence of detected heterozygosity is not evidence of biological identity at low coverage. Worse, `vcf.py:112` catches a query RuntimeError and returns an empty list; a deterministic mocked tool failure confirms this silent conversion. Error, insufficient evidence, same length, and homozygous sequence need distinct handling.

### Additional code-level concerns, not isolated end-to-end failures

- `_build_allele_info` (`alleles.py:333`) can keep a count from the center but a reference contig from a distant refinement. Probe center 72/best contig69 reports total length81 while remapping to a 78-repeat reference. Existing tests preserve this behavior. Report/reference mismatch can induce structural indels and confuse length interpretation.
- `cli.py:433` writes alleles.json before calling; same-length calling later mutates homozygous (`calling.py:366`). The persisted alleles file can disagree with summary.json.
- Standalone consensus searches only allele_1/allele_2 VCF directories (`cli.py:277`); a homozygous same-length call returns `merged/variants.vcf.gz`. The standalone command does not discover that path. Its required input BAM is unused and it lacks custom trimming settings exposed by the library.
- `consensus.py:61` uses exact anchors within a fixed window and silently falls back to fixed trimming. Combined with unmasked uncovered reference sequence, consensus may appear complete where reads provide no evidence.
- Neither AS refinement nor the valley metric explicitly checks full VNTR spanning or soft clipping. Mean raw AS compares different read subsets; all-indel totals can conflate sequencing error, repeat composition, and genuine length mismatch. These are hypotheses to test using anchored spans.

## Performance candidates and falsifiable experiments

1. **Anchor-based candidate length and assignment (highest scientific priority).** Freeze a one-/two-component criterion on training seeds, evaluate separate seeds with gaps 0,1,2,3,4,5, lengths 25–150, HiFi/ONT, and skewed coverage. Compare exact pair accuracy, minority-allele no-call rate, assignment purity, and mutation sensitivity/specificity. Preserve ambiguous reads and alignment fallback. Run a no-ladder-map branch to measure whether anchor scanning plus one/two remaps is faster.
2. **Separate secondary score evidence from support (high priority).** Keep alternative contigs for fit; count collision-free primary molecules for coverage. Reject unsupported splits or mark them uncertain. Acceptance: retain ≥40/44 historical exact pairs while resolving 60/60 and gap3 cases; then test held-out seeds. The primary-only ablation is already a rejected naive approach.
3. **Phase same-length alleles (high priority).** Generate two equal-length haplotypes with variants in cis/trans, shared variants, two distinct ALT alleles, missing calls, and a low-depth allele. Test exact reconstructed sequences with supported bcftools versions and a specified sample/haplotype mode. Rename simulator IDs before phasing. Acceptance: exact haplotype recovery when phase is observed, explicit ambiguity otherwise.
4. **Decompose long-indel loss (high priority).** Fix read set/reference/model, sweep QUAL 0/2/3/5 on retained raw VCFs; separately compare raw candidate generation with truth-local read alignment. Measure WT false frameshifts alongside recovered mutations. The observed QUAL 3.22 insertion motivates this experiment but does not justify a global threshold change.
5. **Stream per-allele remapping (low scientific risk).** `_extract_and_remap_reads` (`calling.py:100`, `calling.py:121`) captures whole FASTQ and SAM strings and writes temporary files, whereas top-level mapping already pipes SAM to sort. Reuse a streaming abstraction, compare BAM primary records/CIGAR/VCF/consensus equivalence, peak RSS, bytes written, and wall time at 50/200/1000× and long reads. No speed gain was measured here.
6. **Budget threads across all concurrent stages.** Two allele tasks each pass the full requested threads to remapping (`calling.py:395`), while only Clair3 splits the budget (`calling.py:405`). Top-level mapper and sorter also overlap with each requesting the configured thread count. Benchmark wall time/CPU/RSS at budgets 1/2/4/8 with long alleles; preserve total budget semantics and avoid claiming linear gains.
7. **Avoid repeated scans and small subprocesses.** Valley splitting plus refinement scans the same cluster twice; retain per-contig sufficient statistics in one pass. The extracted BAM is indexed, read once as FASTQ, overwritten with the remapped BAM, and indexed again; measure whether the first index is unnecessary for this path. Cache validated references/indexes across samples. Prove identical scientific outputs before merging any refactor.

## Validation performed

- `uv run --locked --all-extras pytest tests/unit/test_alleles.py tests/unit/test_mapping.py tests/unit/test_calling.py tests/unit/test_consensus.py --no-cov -q`: **88 passed**.
- In-memory deterministic probes for diffuse coverage loss, flat indel profile splitting, count/reference mismatch, and swallowed VCF query failure.
- Real samtools inspection of input/ladder/remapped records and duplicate names; real bcftools 1.17 opposite-phase consensus probe.
- `PATH=/home/bernt-popp/miniforge3/envs/env_clair3/bin:$PATH uv run --locked --all-extras python tests/results/deep_validation_20260914/mapping_audit/ablate.py`: 44 historical BAMs, three methods, all completed.
- Same command pattern for `anchors.py`: 44 historical BAMs, all completed.

Scripts and result JSON remain in ignored generated results. No full pipeline, performance benchmark, clinical validation, or fresh ONT validation was run by this audit subtask. Parent experiment results must be reported separately with their provenance.

## Follow-up: unsupervised anchor modes on original source reads

After the histogram observation, a bounded exploratory mode-selection prototype was run directly on all 44 original HiFi BAMs and three original ONT FASTQs, scanning both orientations. It uses the same exact 20-base terminal-repeat anchors, rounds span/60, repeatedly selects the highest-support bin (ties favor shorter), removes bins within ±1 repeat, and retains **all** peaks with exact-bin support at least 3, 5, or 10. It never forces two peaks. A single detected length remains one length, with no claim of a homozygous genotype. Truth is loaded only after inference.

Scripts/results: `tests/results/deep_validation_20260914/mapping_audit/anchor_modes.py` and `anchor_modes.json`.

| Platform | Threshold | Exact distinct-length sets | Missed truth lengths | Extra lengths | Single-length outputs |
|---|---:|---:|---:|---:|---:|
| HiFi, 44 samples | 3 | 44/44 | 0 | 0 | 1 |
| HiFi, 44 samples | 5 | 43/44 | 1 | 0 | 2 |
| HiFi, 44 samples | 10 | 43/44 | 1 | 0 | 2 |
| ONT, 3 samples | 3 | 3/3 | 0 | 0 | 0 |
| ONT, 3 samples | 5 | 3/3 | 0 | 0 | 0 |
| ONT, 3 samples | 10 | 1/3 | 2 | 0 | 2 |

There were no entirely empty length-call sets. The HiFi loss at thresholds 5/10 is the 140-repeat allele with only four anchored reads in 25/140. At ONT threshold10, the 80-repeat allele is missed in dupC and normal controls. The 60/60 case yields a single 60-repeat peak, and both gap3 examples yield their two correct peaks.

Source-read anchor success: 6,241/7,257 HiFi records (86.0%; reverse orientation), 112/600 ONT reads (18.7%; forward orientation). Failure to find an ordered pair: 1,016 HiFi, 488 ONT. No ambiguous-orientation reads. The source HiFi count exceeds the historical mapped primary count by six because this prototype does not require mapping to the ladder.

This is exploratory validation on previously examined simulations, not a held-out sensitivity estimate. Threshold3's apparent complete recovery does not establish reliability with three independent reads, and duplicated simulated read IDs still require remediation before phasing. Suppression of ±1 repeat makes adjacent-length alleles unresolvable by construction; this is a known prototype limitation. Most ONT reads fail the exact anchor criterion, so fuzzy anchors or a verified fallback are needed. No mutation calling was rerun using these candidate lengths or assignments.

## Fresh mapping ablation confirmation

The same three-method ablation was repeated after all 44 new HiFi mapping BAMs were available under `tests/results/deep_validation_20260914/hifi`, via `ablate_fresh.py` and `ablation_fresh.json`. Every per-sample inferred length pair matched the historical-BAM ablation; aggregate results are identical. These are fresh alignments of historical reads, not independently generated validation samples.

## Held-out generation and fixed-rule validation

Six additional datasets were specified after fixing the prototype rule and generated with sibling MucOneUp **0.44.5** (`../MucOneUp/.venv/bin/muconeup`). All generation commands succeeded. These use independent seeds 9101–9106 and 200 total template molecules per diploid sample, not 200 final HiFi reads per allele. Exact commands, times, requested mutation targets, and exit codes are preserved in `tests/results/deep_validation_20260914/heldout_data/generation_manifest.json`; the driver is `mapping_audit/generate_heldout.py`.

The private copied `heldout_data/generation_config.json` specifies absolute env_pacbio paths for pbsim/samtools/minimap2/ccs, the sibling `reference/pbsim3/ERRHMM-SEQUEL.model`, two CCS threads, and no human-reference remapping. The sibling configuration was not edited. Tool metadata reports CCS 6.4.0, minimap2 2.30-r1287, samtools 1.23. Each generation took approximately 19–22 seconds. Input FASTQs, sequence truth, structures, stats, and logs are retained under `heldout_data/sample_heldout_*` and ignored by Git.

Fixed rule: select exact-anchor length modes with support ≥3 and suppress neighboring ±1 bins, unchanged from the preceding exploratory experiment. Thresholds 5 and 10 were also reported without tuning the mode rule. `anchor_modes_heldout.py` / `anchor_modes_heldout.json` preserve evaluation.

| Held-out sample | Seed | Requested mutation | Truth distinct lengths | Output at 3/5/10 | Anchored / final reads |
|---|---:|---|---|---|---:|
| normal_60_61 | 9101 | none | 60,61 | 60 | 153/178 |
| normal_60_62 | 9102 | none | 60,62 | 60,62 | 156/182 |
| dupc_60_63 | 9103 | dupC, haplotype1 repeat25 | 60,63 | 60,63 | 127/159 |
| normal_60_60 | 9104 | none | 60 | 60 | 151/172 |
| dupc_hap2_60_80 | 9105 | dupC, haplotype2 repeat25 | 60,80 | 60,80 | 143/173 |
| inscccc_terminal_60_80 | 9106 | insCCCC, haplotype1 repeat55 | 60,80 | 60,80 | 155/169 |

At each threshold: **5/6 exact distinct-length sets**, one missed true length, zero extra lengths, zero empty outputs, and two single-length outputs. The 60/61 failure was predicted before generation because the fixed suppression rule cannot preserve adjacent modes. Across all six: 885/1,033 reads anchored (85.7%), 148 anchor failures, all successful reads forward-oriented. A single length in 60/60 does not establish sequence homozygosity. These held-out experiments validate candidate length behavior only; they do not validate mutation calling or phasing from the proposed assignment.

The fresh baseline raw VCF also reproduces the historical QUAL finding exactly: `hifi/sample_dupc_100_120/allele_1/clair3/merge_output.vcf.gz` contains `contig_91:3192 G→GC`, QUAL 3.22 PASS, GT1/1, DP110, AD22,84. The `probe_phase.py` script preserves the real-bcftools same-length consensus reproduction.
