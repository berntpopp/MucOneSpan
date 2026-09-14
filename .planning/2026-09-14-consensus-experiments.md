# Cached-VCF consensus and QUAL experiments — 2026-09-14

## Design and reproducibility

A paired computational experiment holds the freshly called Clair3 raw VCFs and detected alleles fixed across the 44 existing HiFi and three newly generated ONT datasets. No alignment or Clair3 calls are rerun. The driver substitutes existing remapped BAM/raw VCF paths at the calling boundary, copies the original per-contig FASTA into each experiment directory, and runs the existing `call_variants_per_allele` code, including normalization, PASS/QUAL filtering, and same-length genotype disambiguation. The same-length branch can change its genotype result when QUAL filtering changes; this is intended downstream behavior, with raw evidence held fixed.

The paired modes are the existing `bcftools consensus` command and the identical command with explicit **`-H A`**, selecting alternate alleles in heterozygous genotypes. Installed bcftools was verified as **1.17** through the same cleaned tool environment used by the package. QUAL cutoffs are 0, 3, 5, 10, 20; zero still retains the existing PASS-only filter. `min_dp` behavior is unchanged (the current API ignores it).

Classifier acceleration is process-local: only `classify.edit_distance` is replaced by the evaluator's exact bit-vector Levenshtein function. Across all **47 baseline QUAL5/default datasets**, recomputed **complete** VCF-validated classification dictionaries equal original stored `repeats.json`, trimmed consensus FASTAs are byte-identical, and final allele dictionaries equal original `summary.json`. This is stronger than agreement only on positive calls. The evaluator's function was separately validated against ordinary dynamic programming by the classification audit. No production module was edited.

Artifacts:

- Driver: `tests/results/deep_validation_20260914/consensus_experiments.py`
- First paired-run log: `tests/results/deep_validation_20260914/consensus_experiments.log`
- Sweep log: `tests/results/deep_validation_20260914/consensus_sweep.log`
- Conditions: `tests/results/deep_validation_20260914/consensus_experiments/q{QUAL}_{default|alt}/`
- Each condition stores per-sample filtered VCFs, consensus FASTAs, complete classifications, summary, provenance, elapsed time, and `evaluation.json`.

Command (tool environment path reflects this machine only; no tracked config change):

```bash
PATH=/home/bernt-popp/miniforge3/envs/env_clair3/bin:$PATH \
  uv run --locked --all-extras python \
  tests/results/deep_validation_20260914/consensus_experiments.py
```

The evaluator jointly matches predicted sequences to truth by minimum total sequence edit distance under both haplotype assignments; assignment ties, missing truth haplotypes/events, and extra predictions are retained. Exact mutation matching requires **repeat index, parent repeat, and mutation name on the matched haplotype**. This is a much stricter metric than the old sample-name substring scorer, but it does not yet normalize arbitrary equivalent VCF event representations. Both unfiltered mutation candidates and frameshift+VCF-supported candidates are reported.

## QUAL5 paired findings

| Metric | HiFi default | HiFi explicit ALT | ONT default | ONT explicit ALT |
| --- | ---: | ---: | ---: | ---: |
| Samples | 44 | 44 | 3 | 3 |
| Matched truth alleles | 88 | 88 | 6 | 6 |
| Exact VNTR sequences | 39 | 60 | 3 | 3 |
| Exact repeat structures | 39 | 60 | 3 | 3 |
| Complete exact diploid sequences | 10 | 24 | 0 | 0 |
| Ambiguous consensus bases | 1189 | 0 | 107 | 0 |
| Total sequence edit distance | 9437 | 8556 | 476 | 388 |
| Exact mutation events TP / expected | 24 / 28 | 24 / 28 | 2 / 2 | 2 / 2 |
| Additional mutation candidates | 14 | 15 | 0 | 0 |
| Additional supported frameshifts | 7 | 13 | 0 | 0 |
| Normal samples with any extra candidate | 1 / 16 | 1 / 16 | 0 / 1 | 0 / 1 |

All 94 truth haplotypes are matched in both modes; no missing outputs or assignment ties. The 21 gained exact HiFi sequences are not accompanied by any loss of a previously exact sequence. No matched allele's sequence edit distance worsens, and no truth-to-prediction assignment changes at QUAL5.

However, explicit ALT is **not an accuracy solution by itself**: exact mutation sensitivity remains unchanged, and supported extra mutation calls increase. For example `sample_insg_100_120` haplotype 2 improves from 87 to 4 base edits but increases from two to three extra named/candidate mutations. `sample_dupc_50_55` haplotype 2 improves from 17 to seven base edits while supported extras increase from one to three. Several ambiguous candidate sequences previously lacked frameshift status; replacing ambiguity with an alternate sequence can make their erroneous calls more concrete.

The sole normal sample with extra candidate/support calls remains `sample_gap4_40_44`; its haplotype 2 improves from 23 to 11 sequence edits, still has five extra mutation candidates, and supported extras increase from two to three. Thus disappearance of IUPAC ambiguity must not be presented as disappearance of false positives.

The three ONT datasets all retain one incorrect haplotype reconstruction after ambiguity removal. Even though both expected mutation events are correctly resolved, zero of three complete diploid VNTR sequences are exact. ONT support here is therefore preliminary and limited to two mutation-positive samples plus one normal.

## Threshold sweep — completed

All ten conditions completed successfully, each covering **47 samples / 94 matched truth haplotypes**, with zero missing source samples and no missing truth alleles. QUAL0/ALT has one ambiguous haplotype assignment (`sample_gap3_50_53`, equal total edit cost 456 under both permutations); the other nine conditions have no ties. The sweep reused already completed QUAL5 outputs and took 144 seconds wall time; the first full baseline/ALT pair took 38 seconds. These are cached-call experiment times, not end-to-end pipeline runtimes.

HiFi counts (44 samples, 88 haplotypes; 28 expected mutant events; 16 normal samples):

| QUAL | Consensus | Exact sequences / 88 | Exact diploid / 44 | Exact mutant TP / 28 | Extra candidates | Supported extras | Normal FP / 16 | Total sequence edits |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | default | 39 | 10 | 25 | 18 | 11 | 2 | 9435 |
| 0 | alt | 65 | 26 | 26 | 18 | 16 | 2 | 8533 |
| 3 | default | 39 | 10 | 25 | 15 | 8 | 1 | 9431 |
| 3 | alt | 66 | 27 | 26 | 15 | 13 | 1 | 8533 |
| 5 | default | 39 | 10 | 24 | 14 | 7 | 1 | 9437 |
| 5 | alt | 60 | 24 | 24 | 15 | 13 | 1 | 8556 |
| 10 | default | 33 | 7 | 23 | 9 | 5 | 1 | 9442 |
| 10 | alt | 39 | 9 | 23 | 10 | 10 | 1 | 8628 |
| 20 | default | 14 | 1 | 14 | 1 | 0 | 0 | 9851 |
| 20 | alt | 14 | 1 | 14 | 1 | 0 | 0 | 9359 |

ONT counts (three samples, six haplotypes; two expected mutant events; one normal sample):

| QUAL | Consensus | Exact sequences / 6 | Exact mutant TP / 2 | Extra candidates | Supported extras | Normal FP / 1 | Total sequence edits |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | default | 3 | 2 | 1 | 1 | 1 | 479 |
| 0 | alt | 3 | 2 | 1 | 1 | 1 | 391 |
| 3 | default | 3 | 2 | 1 | 1 | 1 | 479 |
| 3 | alt | 3 | 2 | 1 | 1 | 1 | 391 |
| 5 | default | 3 | 2 | 0 | 0 | 0 | 476 |
| 5 | alt | 3 | 2 | 0 | 0 | 0 | 388 |
| 10 | default | 3 | 2 | 0 | 0 | 0 | 469 |
| 10 | alt | 3 | 2 | 0 | 0 | 0 | 397 |
| 20 | default | 0 | 1 | 0 | 0 | 0 | 568 |
| 20 | alt | 0 | 1 | 0 | 0 | 0 | 547 |

### Interpretation and unresolved issues

- HiFi QUAL3/ALT produces the highest exact reconstruction count among the tested settings: **66/88 exact allele sequences and 27/44 exact diploid sequences**, with **26/28 exact mutation events**. This setting retains 15 additional candidates / 13 supported extra frameshifts and one normal false-positive sample. It is an in-sample result, not a validated threshold recommendation.
- Relative to baseline QUAL5/default, QUAL3/ALT recovers `sample_dupc_100_120` at haplotype 1 repeat 45 (dupC, supporting QUAL 3.22), and `sample_dupcccc_60_80` at haplotype 1 repeat 25 (insCCCC, supporting QUAL 4.47). Lowering QUAL alone recovers one; explicit ALT is needed for the second exact classification in this comparison.
- The lower cutoff adds a supported false positive to the sole ONT normal control; its normal false-positive count becomes **1/1** at QUAL0/3. Thus the setting that improves HiFi sensitivity does not preserve observed ONT specificity in this tiny panel.
- QUAL20 substantially suppresses extra calls but also discards true mutation and repeat-identity evidence. It is not a resolution fix: only 14/94 exact sequences remain with ALT over both platforms, and only 15/30 expected mutation events are resolved.
- Calling-derived repeat counts are identical across conditions (85/94 exact, 93/94 within three) because read mapping and detected allele lengths were held fixed. Consensus cannot restore haplotypes lost or length-distorted upstream. Exact reconstruction and mutation precision still need read assignment and consensus improvements.
- The **two residual HiFi event misses** under QUAL3/ALT are `sample_homozygous_60_60` (mutant haplotype edit distance 31) and `sample_long_120_140` (mutant edit distance 40). Both matched mutant predictions have correct total repeat counts but no called mutation. Cached raw VCF inspection finds no insertion near the expected repeat in either case; a consensus-mode/QUAL change cannot restore evidence absent from the caller output. These require upstream read-assignment/calling investigation rather than explanation by the aggregate improvement. The existing panel has only one deletion, one insCCCC, and two ONT positives, so wider held-out mutation and platform validation remains necessary.

No implementation, default threshold, or public behavior was changed by this experiment. Generated artifacts remain ignored; this report records the paired evidence for review.
