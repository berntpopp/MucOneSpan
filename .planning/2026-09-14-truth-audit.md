# MucOneUp truth and benchmark audit — 2026-09-14

## Scope and evidence status

Read-only audit of existing simulator truth and benchmark code; added three requested ONT datasets in ignored data paths. No production changes or full MucOneSpan pipelines were run by this audit. Commands used the synchronized environment with `uv run --locked --all-extras python`. Shell tools were called through the package `run_tool` for BAM probes. The evidence below is current verification, distinct from historical `.planning/BENCHMARK_RESULTS.md` claims.

## Existing dataset inventory

At audit start there were **44 PacBio datasets, 88 haplotypes**, all with simulated FASTA, structure, statistics, metadata, aligned BAM, and BAM index. Mutant datasets additionally contain the exact mutated repeat FASTA. All 44 statistics files identify MucOneUp **0.44.2**. There were no ONT datasets initially. The checked-in generator defines only 26 of the existing HiFi datasets, plus three ONT definitions; the remaining 18 existing HiFi samples therefore need their original command-line provenance, rather than reconstruction from the current generator alone.

Mutation counts: **22 dupC, 2 dupA, 2 insG, 1 insCCCC, 1 del18_31, 16 normal**. All 28 mutants target haplotype 1 and exactly one repeat. There are no homozygous mutant genotypes, multiple mutation targets, or mutations on haplotype 2 in this set. `sample_homozygous_60_60` is equal-length heterozygosity, with one mutant haplotype, rather than identical haplotype sequences.

| Sample | Total repeat units H1/H2 | Mutation | H1 repeat index (1-based) | CLI seed |
| --- | --- | --- | --- | --- |
| sample_asymmetric_25_140 | 25/140 | dupC | 10 | 1008 |
| sample_bench_5000 | 52/63 | dupC | 20 | 5000 |
| sample_bench_5001 | 50/65 | dupC | 20 | 5001 |
| sample_bench_5002 | 51/58 | dupC | 20 | 5002 |
| sample_bench_5003 | 56/118 | dupC | 20 | 5003 |
| sample_bench_5004 | 38/47 | dupC | 15 | 5004 |
| sample_close_51_58 | 51/58 | dupC | 20 | 5002 |
| sample_del_60_80 | 60/80 | del18_31 | 25 | 1005 |
| sample_dupa_100_120 | 100/120 | dupA | 45 | 1019 |
| sample_dupa_60_80 | 60/80 | dupA | 25 | 1002 |
| sample_dupc_100_120 | 100/120 | dupC | 45 | 1013 |
| sample_dupc_40_50 | 40/50 | dupC | 15 | 1011 |
| sample_dupc_50_55 | 50/55 | dupC | 20 | 1014 |
| sample_dupc_50_57 | 50/57 | dupC | 20 | 1015 |
| sample_dupc_50_60 | 50/60 | dupC | 20 | 1016 |
| sample_dupc_60_80 | 60/80 | dupC | 25 | 1001 |
| sample_dupc_60_80_cov50 | 60/80 | dupC | 25 | 1022 |
| sample_dupc_60_80_s2 | 60/80 | dupC | 25 | 2001 |
| sample_dupc_60_80_s3 | 60/80 | dupC | 25 | 3001 |
| sample_dupc_60_80_s4 | 60/80 | dupC | 25 | 4001 |
| sample_dupc_60_80_s5 | 60/80 | dupC | 25 | 5001 |
| sample_dupc_80_100 | 80/100 | dupC | 35 | 1012 |
| sample_dupcccc_60_80 | 60/80 | insCCCC | 25 | 1004 |
| sample_gap10_40_50 | 40/50 | normal | — | 6107 |
| sample_gap12_55_67 | 55/67 | normal | — | 6108 |
| sample_gap3_50_53 | 50/53 | normal | — | 6100 |
| sample_gap3_80_83 | 80/83 | normal | — | 6111 |
| sample_gap4_40_44 | 40/44 | normal | — | 6101 |
| sample_gap5_55_60 | 55/60 | normal | — | 6102 |
| sample_gap5_70_75 | 70/75 | normal | — | 6109 |
| sample_gap6_30_36 | 30/36 | normal | — | 6110 |
| sample_gap6_45_51 | 45/51 | normal | — | 6103 |
| sample_gap7_60_67 | 60/67 | normal | — | 6104 |
| sample_gap8_35_43 | 35/43 | normal | — | 6105 |
| sample_gap9_50_59 | 50/59 | normal | — | 6106 |
| sample_homozygous_60_60 | 60/60 | dupC | 25 | 1007 |
| sample_insg_100_120 | 100/120 | insG | 45 | 1020 |
| sample_insg_60_80 | 60/80 | insG | 25 | 1003 |
| sample_long_120_140 | 120/140 | dupC | 50 | 1010 |
| sample_normal_100_120 | 100/120 | normal | — | 1021 |
| sample_normal_50_55 | 50/55 | normal | — | 1017 |
| sample_normal_50_60 | 50/60 | normal | — | 1018 |
| sample_normal_60_80 | 60/80 | normal | — | 1006 |
| sample_short_25_30 | 25/30 | dupC | 10 | 1009 |

The historical benchmark report labels `sample_bench_5000–5004` as 60/80; actual stored truth is respectively 52/63, 50/65, 51/58, 56/118, and 38/47. Use stored structure/FASTA rather than that historical shorthand.

## Truth validation and repeat conventions

All **88 existing haplotype sequences reconstruct exactly** from the current sibling configuration's repeat dictionary, the actual stored mutated unit where applicable, and the two 10,000-base hg38 constant flanks. Every reconstructed chain count agrees with `repeat_count`, and every reconstructed sequence length agrees with `vntr_length`. All **28 mutated-unit FASTAs exactly match** the current MucOneUp mutation edit operations. MucOneSpan's complete repeat dictionary is equal to the sibling configuration's dictionary. This establishes no observed sequence drift for these data; it does not substitute for generation-version pinning in future datasets.

MucOneUp `--fixed-lengths N` counts **all repeat units**, including pre-repeats 1–5 and after-repeats 6/6p–9. MucOneSpan `contig_N` counts only the variable canonical units and adds nine fixed units. `detect_alleles(...)[allele]['length']` already adds nine; `canonical_repeats` and contig suffix do not. Compare simulator counts to `length` directly, or compare to `canonical_repeats + 9`. For total 60, the ideal ladder contig is 51. Do not add nine a second time.

Mutation structure indices are 1-based across the entire chain, including pre/after units. Haplotype 1 repeat 25 means variable-region repeat 20. MucOneUp insertion `start=p` inserts **before original base p** (`seq[:p-1] + insertion + seq[p-1:]`), and deletions include both endpoints. For index 25, the unmutated VNTR starts of that unit are 1440 zero-based / 1441 one-based. The dupC/dupA/insCCCC insertion occurs at VNTR interbase offset 1499 (after base 1499); insG is at 1498. del18_31 removes original VNTR one-based bases 1458–1471. Add the actual left flank length for complete haplotype coordinates (10,000 in these truth FASTAs), or 500 for the generated ladder before any alignment transformations. Homopolymer left/right normalization can change VCF anchors while preserving the same alternate sequence; exact event equivalence must be assessed on sequence and aligned context, not raw coordinate equality alone.

## Verified metadata defects and limitations

- All **44** historical `provenance.seed` values are null, despite a valid `--seed` in `command_line`. Recover seeds from those command lines and retain them in a separate authoritative manifest.
- All **28** mutant haplotype `repeat_lengths` arrays still describe unmutated dictionary units. For example dupC 60 repeats reports a sum of 3600 while actual VNTR length is 3601; del18_31 reports 3600 while actual is 3586. Current sibling `simulation_statistics.get_repeat_lengths` still obtains lengths using only `unit.symbol`, ignoring mutated sequence. Use mutated-unit FASTA / reconstructed sequence boundaries for exact lengths.
- All **44** PacBio metadata `Coverage` fields say 30, while command lines specify 200 for 43 datasets and 50 for one. Current sibling source explains the mismatch: amplicon simulation prefers `read_simulation.coverage`, while metadata reads `pacbio_params.coverage`. The current CLI explicitly defines this amplicon parameter as the **total template molecule count before CCS filtering**, split across haplotypes by PCR bias. It is not achieved per-allele sequencing depth. Do not label these datasets simply 200x/50x without measured depth or molecule semantics.
- Primary BAM record counts are 156 (baseline dupC), 44 (low-coverage dupC), 174 (asymmetric), 169 (equal-length), and 177 (long). Thus requested molecules are not retained read counts.
- **Duplicate QNAMEs represent distinct molecules.** Baseline dupC has 156 primary records but 113 distinct names; all 43 repeated names have distinct sequences, with lengths consistent with opposite haplotypes (e.g. `S/1/ccs` lengths 3624 and 4850). Equal-length has 169 records / 100 names; low-coverage 44/33; asymmetric 174/169; long 177/119. `samtools fastq` preserved all 156 baseline and 169 equal-length records in the exact current coordinate-sorted inputs, so loss was **not demonstrated**. Unique-name assumptions, name-collation, or molecule assignment remain unsafe. Future generation should namespace source/haplotype read IDs and retain provenance.
- Existing simulator output does not include a read-source manifest. The current amplicon CLI explicitly rejects `--track-read-source`, so per-read assignment truth cannot be assumed present.

## Existing benchmark correctness gaps

`scripts/benchmark.py` measures stage wall times only. It discards classification results and does not load truth, calculate sensitivity/specificity, compare allele lengths, compare consensus sequences, verify repeat structures, or resolve mutation coordinates. It uses only `*.bam` and silently skips FASTQ-only samples; default mapping/calling is HiFi. Its missing-data branch prints an error and returns exit code zero. Verified command:

```bash
uv run --locked --all-extras python scripts/benchmark.py --data-dir /tmp/muconespan-audit-nonexistent-data
```

Output: `Error: /tmp/muconespan-audit-nonexistent-data not found. Generate test data first.`; exit status **0**.

`scripts/batch_analyze.py` loads structure truth, but scoring is sample-level template-name substring detection. It collects every mutation, without requiring frameshift or VCF support; `has_mutation` and expected haplotype structures are not used for matching. Allele detection is recorded only as a boolean. Repeat identity, mutation position, haplotype assignment, extras, and consensus sequence are not scored. Missing truth (`unknown`) is treated as expected mutant. The CLI does not signal per-sample failures via a nonzero aggregate exit status. Its top-level input discovery excludes `.fq` despite the runner supporting it, and excludes compressed FASTQ; the timing runner excludes all FASTQ.

Executed deterministic probes against the actual current function:

| Expected | Constructed calls | Existing status | Why insufficient |
| --- | --- | --- | --- |
| dupC | dupC, template=True, repeat=999, allele_2, VCF=False | TP | Wrong position/haplotype and no support accepted |
| dupC | dupCCCC, template=True | TP | Substring collision accepted |
| dupC | del18_31, template=True | TP_partial | Wrong mutation is a missed expected event plus an extra event |
| dupC | dupC plus dupA | TP | Additional false call not counted |

Reproduction:

```bash
uv run --locked --all-extras python - <<'PYPROBE'
import runpy
analyze = runpy.run_path('scripts/batch_analyze.py')['analyze_results']
truth = {'mutation': 'dupC', 'haplotypes': {}}
for name, call in [
    ('wrong_site', {'mutation_name': 'dupC', 'template_match': True,
                    'repeat_index': 999, 'vcf_support': False}),
    ('substring', {'mutation_name': 'dupCCCC', 'template_match': True}),
    ('wrong_type', {'mutation_name': 'del18_31', 'template_match': True}),
]:
    result = {'classifications': {'allele_2': {'mutations': [call]}}}
    print(name, analyze(name, truth, result)['status'])
PYPROBE
```

Current integration tests validate six repeat-length cases to **±2 repeats**, with two expected failures, and full pipeline dupC presence plus a normal negative. They do not establish exact length, full exact structure, phased event precision, other mutation sensitivity, or ONT accuracy.

## Proposed metrics and acceptance definitions

Keep detection, reconstruction, and exact resolution distinct:

1. **Execution completeness:** datasets attempted, failed, no-calls, excluded with reason, and completed. Report all denominators. Never convert failed/no-truth runs to success or silently exclude failed positives from sensitivity.
2. **Diploid length reconstruction:** jointly match the two predictions to truth under both permutations using sequence alignment where available; record an ambiguity flag if mapping cannot be uniquely resolved. Report exact pair accuracy, both-within-1 and within-2 repeat accuracy, per-allele signed error/absolute error, and dropout/false split frequency. `same_length` does not imply identical haplotypes. Duplicate predictions from one cluster must not imply two independently recovered sequences.
3. **Consensus sequence reconstruction:** trim verified biological boundaries, globally align each matched VNTR, then report exact diploid-sequence accuracy, per-allele edit distance/identity, inserted/deleted/substituted bases, and length error. Do not assess only the synthetic ladder-length interval. Report flanking sequence reconstruction separately.
4. **Repeat structure reconstruction:** compare ordered repeat units and actual unit sequences using dynamic-programming alignment; report exact structure, edit distance in units, boundaries recovered, and dictionary sequence-equivalence classes where labels are inherently ambiguous. Penalize missing and extra repeats. A correct unordered repeat histogram is not exact structure.
5. **Sample mutation detection:** define an actionable call rule explicitly (e.g. named frameshift with sequence match and VCF evidence), then compute TP/FN on mutant samples and FP/TN on normal samples. Here initial denominators are 28 positive/16 negative. Separately report unfiltered candidate findings. Wrong mutation may count as an alarm for broad mutation-presence sensitivity, but cannot count as the correct genotype.
6. **Exact mutation event resolution:** represent truth by matched haplotype, original repeat index, original base coordinates, normalized ref/alt sequence, and canonical name/explicit aliases. One-to-one match predicted events to truth. A correct call requires equivalent alternate sequence at the corresponding biological locus; a substring name or nearby variant is insufficient. Count unmatched truth as FN and unmatched calls as FP, including extra calls in positive samples. Report event precision, recall, F1; allele assignment accuracy; exact repeat localization; bp localization; and complete diploid genotype accuracy.
7. **Uncertainty and stratification:** report counts plus binomial confidence intervals for sensitivity/specificity, stratified by platform, mutation class, repeat length/gap, achieved per-allele spanning depth, molecule bias, and independently generated seed. Repeated equal biological configurations are not equivalent to broad mutation validation; one deletion/four-C example gives weak per-class evidence.
8. **Runtime:** retain mapping/calling/consensus/classification times, total wall time, CPU/thread budget, peak RSS, input read/base counts, software/model versions, warm/cold cache state, and truth-manifest hash. Compare the same inputs/configuration and report reruns rather than reusing historical artifacts as current results.

An authoritative manifest should include source software version/commit, full command, config hash, truth FASTA hashes, explicit simulator seed, platform/error-model hash, requested template molecules, actual read/molecule counts, source read IDs, actual flanks, repeat sequences/lengths, repeat count convention, mutation events, and any expected ambiguity.

## Fresh ONT datasets generated during this audit

On parent request, generated `sample_ont_dupc_60_80`, `sample_ont_dupa_60_80`, and `sample_ont_normal_60_80` using current sibling **MucOneUp 0.44.5**, seeds 1001/1002/1006, fixed total repeats 60/80, mutation on haplotype 1 repeat 25 where applicable. Output FASTQs are `tests/data/generated/<sample>/<sample>_reads_amplicon_ont.fastq`. All three have **200 FASTQ records** and 131 distinct names; 69 QNAME collisions persist in the current source. Statistics `provenance.seed` remains null.

A copied config at `tests/results/deep_validation_20260914/ont_generation_config.json` supplies explicit installed pbsim3/samtools/minimap2 paths, absolute QSHMM-ONT-HQ model path, and two threads. It omits `read_simulation.human_reference` to retain FASTQ without unnecessary whole-genome alignment. Haplotype simulation used the original sibling configuration. Generation commands and outputs are logged as `<sample>_simulate.log` and `<sample>_reads.log` in that results directory. Source tracked files were not modified. These new ONT data have different simulator provenance from the original 44 HiFi datasets and must be reported separately.

Read generation completed successfully; this audit does **not** claim ONT pipeline performance or accuracy. Parent pipeline validation will supply those results.
