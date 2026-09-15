# Wave 2 missed-call debugging and isolated MP1 correction

Investigation date: 2026-09-16. Frozen caller `c08000a0b97e9ee00e560a369ef320b4d06d3004`; final benchmark harness `1deb5dcfd1dbeede0fbd3ca859cc9a6aab011d47`. Original outputs remain under external `${DATA_ROOT}/cohort-v2`. Diagnostic scripts/results remain under `${DATA_ROOT}/debug-misses`. No participant sequences or read identifiers are copied into Git. No threshold, caller model, input filter or original benchmark output was changed.

The user subsequently authorized fixing the demonstrated MP1 classifier bug. That correction is isolated on `fix/dupc-repeat-context` in `.worktrees/fix-dupc-repeat-context`; the frozen benchmark checkout's scientific implementation remains untouched. This document distinguishes the frozen result, diagnostic counterfactuals, and the isolated correction.

## Findings with causal evidence

### MP1: retained C insertion becomes an A insertion through wrong repeat segmentation

The existing full consensus exactly replays all retained VCF records using `variant_support._replay`. The target G>GC record at ladder position 1512 is selected as GC by GT1 from genotype 1/1. Therefore neither QUAL filtering nor consensus construction removed this C insertion.

On the actual trimmed candidate, interval `[960,1021)` is **exactly** the existing dictionary B motif with one extra C in its terminal seven-C tract. The frozen bundled mutation definition permits `dupC` only on X (`data/repeats/repeats.json:63`; template construction `config.py:90`). B differs outside this tract, so no exact B:dupC template exists.

The next loss occurs in `classify._forward_classify`:

- Exact-template phase (`classify.py:285`) finds no eligible 61 bp B template.
- The canonical 60 bp fallback (`classify.py:324`) classifies `[960,1020)` as a one-base substitution against B (edit distance 1).
- Instrumented replay records only 60 bp and then 30 bp fallback windows at this position. The 30 bp window has distance 30, but the inherited best distance is already 1; the early-stop condition (`classify.py:350`) exits before evaluating the correct 61 bp window. Additionally, replacement requires strictly lower distance, so merely continuing the loop would not resolve a 60/61 bp distance-1 tie.
- The chosen boundary shifts the terminal A of the true B unit into the next window. Repeat 18 `[1020,1081)` is then described as an A insertion before X. This is a decomposition error in the same unchanged DNA sequence, not a literal C-to-A base conversion during consensus.

**Controlled proofs:**

1. Removing only the retained target C from the full consensus eliminates the downstream classified mutation entirely.
2. A minimal dictionary-derived synthetic sequence `B-with-one-extra-tract-C + X + X` reproduces `?B Xm X`, with the spurious A insertion in repeat 2. Replacing only the motif background with X yields `X:dupC X X`.
3. Evaluating the actual 61 bp interval directly against B identifies the C insertion; classifying the suffix beginning at 1021 produces no mutations.
4. Supplying that correct interval to the existing exact VCF reversion validator produces `exact_sequence_concordance` at repeat 17 for the original G>GC record. The frozen wrong A decomposition receives `absent` support as expected.

These isolate dictionary applicability plus candidate segmentation as the cause. They do not independently prove the participant's complete haplotype or authorize treating comparator position observations as truth.

### MP3: separate-length calling ignores the returned haplotype policy and selects GT1

The retained short-candidate VCF contains the target position-912 G>GC record with genotype 0/1, AF 0.2794, QUAL 21.05. It also has three heterozygous substitution records. `vcf.filter_vcf` retains the original mixed genotype for AF in `[0.2,0.5)` (`vcf.py:116`). `phase_evidence` correctly returns `phase_status=unphased` and `haplotypes=['I']` for the four unphased heterozygous sites (`phasing.py:39`).

However the separate-length branch unconditionally passes haplotype **1** to `annotate_consensus_candidate` (`calling.py:529`), rather than using or preserving the returned policy. `consensus.build_consensus` passes this to `bcftools consensus -H 1` (`consensus.py:48`); genotype position 1 selects reference G at 0/1 and therefore omits the retained alternate C.

**Controlled proofs using existing `tools.run_tool`, each external command bounded to 30 seconds:**

- Real `bcftools consensus -H 1` against the frozen reference and VCF reproduces the original complete consensus exactly.
- Changing only `-H 1` to `-H 2` yields a named, resolved X:dupC at repeat 7. This changes the selection at four unphased sites, so it is not a defensible phased reconstruction.
- A tighter replay changes only the target site's genotype to 1/1 in an external diagnostic VCF, leaving the other records and GT1 selector unchanged. Its full consensus differs from the baseline by exactly one edit and classifies named X:dupC at repeat 7. This proves the target-site selection is sufficient to explain this loss.

Neither diagnostic is a proposed instruction to force all alternate alleles or a replacement for proper phase handling. The original run remains uncallable and its frozen score remains unchanged. Follow-up issue: [#53](https://github.com/berntpopp/MucOneSpan/issues/53).

### MP2 and MP4: target indels are absent before filtering

Inspection of original Clair3 `merge_output.vcf.gz` and final `variants.vcf.gz` gives:

| Run / allele | Raw records | Retained records | Raw indels | Retained indels |
| --- | ---: | ---: | ---: | ---: |
| MP2 allele 1 | 43 | 42 | 0 | 0 |
| MP2 allele 2 | 3 | 3 | 0 | 0 |
| MP4 allele 1 | 37 | 37 | 0 | 0 |
| MP4 allele 2 | 4 | 4 | 0 | 0 |

Thus lowering the final VCF threshold cannot restore an indel missing from these raw candidate calls. Additional read-level tracing now identifies the earlier representation loss, without rerunning a model.

#### Physical long reads survive, but enter the shorter allele's BAM

Using the same exact conserved anchors as the prespecified input audit, a broad diagnostic range of 4200–5400 bp identifies 153 MP2 and 127 MP4 full-span reads. These ranges describe raw physical spans, not independent haplotype or mutation truth. Each read identifier was traced through original ladder BAM flags/contigs and the final per-allele remapped BAM. No input was filtered or modified.

| Measured lineage | MP2 | MP4 |
| --- | ---: | ---: |
| Exact-anchor long spans, 4200–5400 bp | 153 | 127 |
| Long spans with a primary alignment in selected allele 1 cluster | 153 | 127 |
| Long spans actually present and primary-mapped in remapped allele 1 BAM | 153 | 127 |
| Long spans entering selected allele 2 / its remapped BAM | 0 | 0 |
| Original reads shorter than 1200 bp entering remapped allele 2 BAM | 242 | 655 |
| Allele 1 selected reference, total units | 50 | 45 |
| Allele 2 selected reference, total units | 14 | 15 |
| Median long-read insertion bases in allele 1 remap CIGAR | 1685 | 2098 |
| Median long-read soft-clipped bases in allele 1 remap | 75 | 75 |
| Median long-read reference span in allele 1 remap | 3031 bp | 2731 bp |

The long raw observations are **not removed by FASTQ conversion**: all 153/127 survive in the major cluster, then are aligned to the shorter reference as large insertions. Their separately observed length mode has been lost as an allele candidate. The second candidate instead contains hundreds of short fragments and no traced long full-span observations. This establishes a candidate-representation failure; it does not prove which one-base event each noisy long read supports.

Raw exact-span histograms provide independent input-geometry checks: MP2 has 548 observations in the nearest-60-bp 3000 bin and 115 in 4680 (plus 24 in 4620); MP4 has 1001 in 2700 and 95 in 4800 (plus 18 in 4740). Rounding is for diagnostic visualization only and was not applied to caller input or truth.

#### Unchanged selection replay explains why fragments win

Replayed the actual `_find_clusters`, `split_cluster_by_read_length`, `_split_cluster_by_indel`, and `detect_alleles` functions against the frozen BAMs. A local pysam iterator supplies exactly the requested SAM records and flag filters to the existing read-iterator boundary; all computational functions and settings are unchanged. Both final selected contig names, reported lengths, cluster-contig lists and alignment counts exactly match the original stored results.

1. Both runs start with **one** broad passing cluster, canonical contigs 1–72 (MP2) and 1–75 (MP4). The earlier hypothesis that an initial two-cluster branch bypassed splitting is rejected.
2. The physical-read-length splitter returns no split. Its 5 bp bins require at least 15% of eligible primary alignments in an individual bin (`length_candidates.py:247`). MP2 has 4572 eligible records, threshold 685, largest bin only 294; MP4 has 7395, threshold 1109, largest bin only 537. Neither has even one qualifying peak. Separately, this near-length splitter's accepted mode difference is only 45–320 bp (`length_candidates.py:262`), so the observed roughly 1.7/2.1 kb separation is outside its intended regime.
3. The CIGAR-indel valley splitter chooses the two lowest local mean-indel valleys (`alleles.py:311–324`). MP2 chooses canonical 41 (mean 71.9 bp) and **5** (77.5 bp); the long-length canonical 69 valley has 112.5 bp. MP4 chooses 36 (64.4 bp) and **6** (81.8 bp); long-length canonical 71 has 114.7 bp. These are distinct values, not ties resolved by AS.
4. Splitting by nearest selected valley leaves the long and major short populations together. The real selection replay accepts these fragment/major pairs through the existing dominance gate. A fragment's better alignment to a short contig therefore passes without demonstrating an independent full VNTR allele.
5. Final ONT refinement also limits eligible contigs to at least 25% of the largest cluster alignment-record count (`alleles.py:171`). The long-length contig 69 has 701 records below MP2's 730 cutoff; long-length contig 71 has 549 below MP4's 1207 cutoff. These are alignment-record counts, often dominated by secondary records, not molecule counts. The selected shorter contigs have lower mean indel values anyway.

This causally locates why the separate long candidate is missing. Whether a correctly separated long candidate would yield the intended raw one-C variant under the frozen Clair3 model remains untested. A later correction should evaluate complete-span evidence and competing length modes before choosing a strategy, without simply lowering a global threshold or forcing the comparator's reported lengths.

## Isolated MP1 fix: exact tract applicability, without fuzzy matching

Issue [#52](https://github.com/berntpopp/MucOneSpan/issues/52). The smallest supported correction expands only the bundled `dupC.allowed_repeats` from X to the existing 60 bp motifs with **G at base 52, exactly seven Cs at bases 53–59, and a non-C at base 60**: `5, X, B, D, K, W, M, N, O, P, Q, R`.

The mutation remains the same one-C insertion before base 60; exact template lookup remains mandatory. Other mutation definitions, sequence matching, window-search thresholds, VCF support, and clinical decision rules remain unchanged. This allows the existing exact-template phase to consume the correct 61 bp B unit, preserving its actual parent and boundary before fallback can mis-segment it.

Source basis: Madritsch et al. 2025 defines the common MUC1 duplication as a single C insertion within a seven-C homopolymer in one 60 bp VNTR unit (Introduction, PDF p1; Discussion, PDF p7; DOI `10.1038/s41598-025-30441-3`). The original Kirby 2013 citation remains. The added applicability is an exact sequence-context interpretation of that event definition, not independent clinical confirmation for each motif. Motifs with shorter/interrupted tracts are excluded. Fixed precursor motif 5 is inside the boundary-inclusive VNTR and contains this identical tract: "fixed" describes its expected structural placement, not an exception to its nucleotide sequence's insertion identity. This change does not establish patient prevalence or a separate clinical claim for motif 5.

### Red/green and real-data verification

New focused synthetic tests initially failed in 13 cases, including the shifted B boundary and loss of exact VCF support. After the dictionary-only change, all 25 initial new cases and 67 existing classifier/config tests pass (92 total). An additional explicit precursor-5 boundary/normal-control regression also passes (26 new cases). Before that additional test and the lead’s documentation changes, `make ci-check` passed 899 tests with 91.35% branch-aware aggregate coverage; `make build-check` passed wheel/sdist content, installation, CLI and bundled-resource verification. The lead runs final integration checks on the complete revision. They cover every added motif's precise tract, canonical X, unchanged B/X, ineligible shorter/interrupted tracts, wrong-position C insertions, other inserted bases, and homopolymer-equivalent VCF support.

A separate read-only replay reclassified **all 22 frozen cohort consensus sequences** and validated their mutations against the unchanged original VCF/context. Exactly one allele's mutation list changed: MP1 allele 2 now has `B:dupC` at repeat 17, resolved localization, `vcf_support=true`, `vcf_support_status=exact_sequence_concordance`, QUAL 35.58. The wrong downstream A mutation disappears. All 21 other allele mutation lists are unchanged. Original files and final cohort scores were not modified; this is a post-baseline correction replay, not the original frozen benchmark result.

## Reproduction and evidence files

From the frozen benchmark checkout, with the locked environment:

```bash
uv run --locked --all-extras python "$DATA_ROOT/debug-misses/trace_misses.py"
uv run --locked --all-extras python "$DATA_ROOT/debug-misses/isolated_replays.py"
```

`trace-results.json` records VCF filtering comparisons, MP1 complete replay and corrected-boundary support, synthetic B/X controls, MP3 haplotype selection and actual commands. `isolated-results.json` records actual fallback window visits and the one-site MP3 replay. Original and diagnostic full FASTAs/VCFs remain external.

From the isolated fix checkout:

```bash
uv run --locked --all-extras pytest tests/unit/test_dupc_context.py tests/unit/test_classify.py tests/unit/test_config.py --no-cov
uv run --locked --all-extras python "$DATA_ROOT/debug-misses/verify_dupc_context_fix.py"
```

The latter writes `fix-replay.json`, with input FASTA hashes and all before/after mutation annotations. It never writes into `cohort-v2`. Independent review found no blocker and additionally tested all 1156 ordered normal motif pairs without generating a named dupC. The lead owns final-revision integration checks and merging; this file does not claim the fix has merged.


Read lineage and exact selection replay use the already installed external IGV-report environment's pysam without modifying dependencies:

```bash
PYTHONPATH="$BENCHMARK_CHECKOUT/src" "$DATA_ROOT/igv-reports-1.13/bin/python" "$DATA_ROOT/debug-misses/read_lineage.py"
PYTHONPATH="$BENCHMARK_CHECKOUT/src" "$DATA_ROOT/igv-reports-1.13/bin/python" "$DATA_ROOT/debug-misses/replay_selection.py"
```

`read-lineage.json` stores counts, span histograms, contig memberships and CIGAR summaries; `selection-replay.json` stores peak thresholds, valley choices and exact baseline-selection equality. Diagnostic script development encountered an unavailable pysam import in the benchmark uv environment, then used the existing external environment; the first SAM-iterator replay needed hexadecimal flag-mask parsing corrected before completing. Neither setup error modified data or algorithm behavior.
