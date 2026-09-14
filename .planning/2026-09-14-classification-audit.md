# Classification and consensus-truth audit — 2026-09-14

Scope: investigate repeat/mutation accuracy independently of mapping/calling,
then compare reconstructed consensus sequences against individual truth
haplotypes. Production files were not changed. Evidence scripts and outputs are
ignored under `results/classification-audit/`.

## Reproduction and truth conventions

Run from repository root with the locked environment:

```bash
uv run --locked --all-extras python results/classification-audit/probe.py
uv run --locked --all-extras python results/classification-audit/adversarial.py
uv run --locked --all-extras python results/classification-audit/evaluate_pipeline.py tests/results/deep_validation_20260914/hifi
uv run --locked --all-extras python results/classification-audit/benchmark_distance.py
uv run --locked --all-extras python results/classification-audit/paired_timing.py
uv run --locked --all-extras pytest tests/unit/test_classify.py tests/unit/test_config.py tests/unit/test_vcf.py --no-cov -q
```

The focused existing unit suite passed: **81 passed in 5.95 s**. These tests do
not validate external tools; the parent audit runs those separately. An initial
convenience attempt to import BioPython failed because it is not a dependency;
the evidence scripts use a small FASTA parser and install no packages.

Truth FASTA records are parsed separately. Both complete 10,000-base dictionary
flanks must match exactly before trimming; no fixed guessed boundaries are used.
Expected structure comes from MucOneUp's `vntr_structure.txt`. Mutation labels
such as `Xm` become `X:dupC` using the simulation's mutation metadata. Repeat
indices are 1-based and include the five pre- and four after-repeat units.

The evaluator performs unordered **one-to-one** haplotype assignment minimizing
total exact global nucleotide Levenshtein distance. It reports assignment ties,
missing haplotypes, extra predictions, mutation name/parent/index errors,
supported extra calls, exact nucleotide/structure matches, edit distances, and
the `canonical_repeats + 9 == length` convention. It includes failed entries in
`measurements.json` even when `summary.json` is absent. An optional
`--expected-samples` JSON list includes unattempted inputs. Missing truth
mutations count as false negatives, and sequence/structure denominators include
all expected truth haplotypes. Missing wild-type haplotypes do not add mutation
false negatives. `--truth-dir` aliases `--truth-root` for perturbation datasets.

The evaluator's dependency-free Myers bit-vector distance was checked against
the production scalar DP on 500 deterministic random DNA pairs, lengths 0–100.
This makes whole-consensus comparisons practical without concatenating records.

## Error-free sequence classification

- Original 44 generated datasets: **88/88 haplotypes** match the complete truth
  structure, mutation names, parents, and repeat positions; 349,936 total bases
  classified in 0.0295 s in one run.
- Including three newly generated ONT truth datasets: **94/94 haplotypes** match
  exactly; 0.0325 s in the later run. This checks the ONT source sequences, not
  the error-containing ONT reads or the ONT pipeline.
- Dictionary matrix: **132/132 exact full structures** for all 44 supported
  `(parent repeat, mutation)` templates followed by X, A, or terminal repeat 9.
  Matrix runtime was 0.00228 s in the later run.
- No sequence collisions among the 44 mutation templates or between mutation
  templates and ordinary dictionary repeats were found.

The mutation-template matrix exercises the current dictionary's own coordinate
application. Independent agreement with simulator output comes from the 28
mutated original truth haplotypes; those cover dupC, dupA, insG, insCCCC, and
del18_31. Other dictionary mutations are not independently validated by that
original read dataset matrix. The truth-audit agent separately reconstructed
the generated sequences from current MucOneUp definitions.

## Fresh pipeline results against full truth

Final evaluator outputs are
`tests/results/deep_validation_20260914/{hifi,ont,perturbation}_accuracy.json`.
The additional ONT and perturbation inputs are described in the parent audit.

| Metric | Original HiFi matrix | New ONT matrix | Depth/name perturbations |
| --- | ---: | ---: | ---: |
| Expected samples | 44 | 3 | 24 |
| Expected truth haplotypes | 88 | 6 | 48 |
| Failed samples / missing summaries | 0 | 0 | 6 |
| Missing truth haplotypes | 0 | 0 | 13 |
| Exact sequence and complete structure | 39/88 | 3/6 | 7/48 |
| Exact called repeat count | 82/88 | 3/6 | 27/48 |
| Called count within ±3 repeats | 87/88 | 6/6 | 28/48 |
| Samples with both counts within ±3 | 43/44 | 3/3 | 10/24 |
| Exact mutation parent/name/index TP / FN | 24 / 4 | 2 / 0 | 7 / 7 |
| Extra mutation records | 14 | 0 | 12 |
| VCF-supported extra frameshifts | 7 | 0 | 7 |

All assignments in these three matrices had unique minimum total sequence edit
distance. Uniqueness within these two-reference comparisons does not establish
biological phasing certainty. Correct repeat counts do not imply correct repeat
composition, as the HiFi 82/88 count versus 39/88 complete structure result shows.

The four HiFi mutation false negatives are `sample_dupc_100_120` (dupC, repeat
45), `sample_dupcccc_60_80` (insCCCC, repeat 25), `sample_homozygous_60_60`
(dupC, repeat 25), and `sample_long_120_140` (dupC, repeat 50). The seven
supported extra frameshifts occur in `sample_bench_5003` (one),
`sample_dupc_50_55` (one), `sample_gap4_40_44` (two),
`sample_insg_100_120` (two), and `sample_short_25_30` (one). None has a known
mutation-template name. There are 1,189 ambiguous consensus bases across the
HiFi matrix. These are pipeline findings, not failures of the exact-truth
classifier control.

## Reproducible limitations and defects

1. **Mixed insertions/deletions can receive the wrong frameshift flag.** For X,
   `x[:5] + x[7:30] + 'A' + x[30:]` has length 59 versus the 60-base reference.
   `classify_repeat` finds two deletions and one insertion, edit distance 3,
   but returns `frameshift=False`. `classify.py` checks the sum of inserted and
   deleted lengths modulo 3, rather than their signed net change. The net -1
   clearly changes the downstream frame. The converse, one insertion plus one
   deletion, is marked frameshift despite restoration of the downstream frame;
   any interpretation of the temporarily altered segment should be explicit.
   Evidence: `net-frameshift.json`, `balanced-indels.json`.

2. **The advertised backward rescue cannot currently rescue forward failures.**
   The forward loop stops normally only with fewer than 30 bases remaining.
   `_apply_bidirectional_fallback` requires more than 30 unconsumed bases.
   Otherwise the forward loop consumes even poor matches or raises, so the
   rescue condition cannot be reached through normal public classification.
   A 40-base insertion into the third of five X repeats is classified as
   `X X Em Mm X X`, reports offset -20 instead of the real +40, and creates two
   mutation calls. This insertion exceeds the configured ±30-base probe range;
   that range limitation is expected, but the result does not expose failure
   to segment the region reliably. Evidence: `adversarial.json`.

3. **Unconsumed tails are excluded from confidence.** `X * 5 + 'N' * 29`
   returns exactly five X labels, confidence 1.0, and 100% exact matches without
   reporting the 29 unclassified bases. Confidence describes classified
   windows, not full-sequence completeness. Evidence: `adversarial.json`.

4. **VCF support means a nearby record, not mutation-specific support.** For a
   dupC in repeat 2 and flank length 500, records anywhere from position 560
   through 650 inclusive set `vcf_support=True`. Position 640 is inside the
   next repeat. An unrelated SNP at 585 also suffices; the parser discards REF,
   ALT, and genotype for this validation path. This is consistent with the
   implementation's coarse window but cannot establish support for the named
   indel. Reference/consensus coordinate projection must be handled explicitly
   if replacing it: naive use of cumulative consensus offsets would also be
   wrong against VCF reference coordinates. Evidence: `vcf-support.json`.

5. **Input normalization is a caller responsibility.** A lowercase X is
   classified as `?1` at zero confidence; reverse-complement X yields `3m`,
   confidence 0.533. Standalone classification does not detect orientation or
   normalize case. Mapped consensus is reference-oriented, so these probes do
   not demonstrate a strand failure in the pipeline. Evidence: `adversarial.json`.

6. **Insertion-coordinate comment is misleading.**
   `characterize_differences('ACGT', 'ACAGT')` reports insertion position 3,
   meaning before reference base 3. The inline comment says the position is
   after which insertion occurs. Dictionary application itself uses before
   the 1-based start, and inclusive deletion endpoints, correctly for the
   independent truth cases above.

## Profiling and optimization evidence

Fresh sample_bench_5003 allele 2 has 7,259 consensus bases, versus 7,080 truth
bases, and 286 global sequence edits. The classifier emits 121 repeats for a
called/truth count of 118, and 59 structure edits. Its average confidence is
**0.9828**, although only **51.2%** of repeat matches are exact. The consensus
contains 106 S, 13 M, and 4 R ambiguity symbols. It also produces a VCF-supported
false frameshift at repeat 116 of the wild-type haplotype; the boundary penalty
does not apply because the classifier's erroneous 121-repeat count moves the
boundary window.

`profile_5003.py` measured 33.14 s under cProfile: 80,240 edit-distance calls cost
31.97 s (**96.5%**), and 2,360 difference tracebacks cost 1.09 s. The source is
repeated 34-reference × many-window fallback comparisons, not exact-template
lookups. Profiler timing must not be compared directly with unprofiled timing.

`benchmark_distance.py` changes `classify.edit_distance` only in its Python
process to use the verified bit-vector algorithm. It compares complete returned
dictionaries, including difference traceback, confidence, offsets, mutation
records, and structure. `paired_timing.py` runs three alternating unprofiled
original/replacement comparisons on the same consensus, plus an exact 120-X
control. Other audit jobs may be running concurrently. Timing and full matrix
results are stored in `bitvector-benchmark.json` and `paired-timing.json`.

The three paired unprofiled runs gave:

| Input and scoring | Median | Range |
| --- | ---: | ---: |
| sample_bench_5003 allele 2, original | 19.3504 s | 18.8010–19.6835 s |
| Same consensus, bit-vector | 2.3965 s | 2.3229–2.4000 s |
| Exact 120-X control, original | 0.5815 ms | 0.5787–0.6140 ms |
| Exact 120-X control, bit-vector | 0.5807 ms | 0.5785–0.5843 ms |

Thus the difficult example was **8.07× faster**, with complete dictionary
equality for every run and no meaningful change to the already-fast exact-match
control. A separate sweep of 74 consensus outputs that existed when it started
also had complete returned-dictionary equality; aggregate unprofiled time was
123.62 s original versus 14.56 s bit-vector. The truth-audit agent independently
checked byte-identical consensus and full VCF-validated classification results
for all 47 completed HiFi/ONT samples during its default-consensus ablation.

Potential follow-ups, not implemented:

- Replace scalar edit-distance scoring with a proven exact bit-vector or native
  implementation. Preserve ties and traceback behavior; dependency-free Python
  avoids package changes, while a native library might be faster but adds
  wheel/platform and distribution checks.
- Compute traceback only for the winning window, and cache repeated reference
  masks. This reduces secondary overhead but should follow the measured main
  scoring bottleneck. Cached mutable result dictionaries require copying.
- Make segmentation completeness and unresolved ambiguous bases visible instead
  of treating a high average per-window identity as a calibrated probability.
- Introduce a meaningful stop/recovery policy before backward rescue, with
  synthetic large-indel truth and terminal-boundary tests. Global segmentation
  could improve robustness but costs complexity and requires explicit tie rules.
- Validate named indels using normalized variant alleles and coordinate-aware
  projection. A nearby high-QUAL SNP must not establish named mutation support.
- Correct signed frameshift accounting with focused mixed-indel tests.

The parent audit traced ambiguous consensus output to bcftools 1.17's default
genotype/IUPAC behavior and is testing explicit haplotype selection separately.
That changes scientific behavior and must be judged against mutation and full
structure truth metrics, not runtime alone. No threshold or minimum confidence
has been calibrated by the classifier experiments.
