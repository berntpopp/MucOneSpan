# Ladder 2×2 diagnosis: caller version × read simulation (Task 12d)

Date: 2026-09-24/25. Worktree `MucOneSpan-bench` (feat/benchsim, 0d14ae5; caller
code identical to v0.16.0). Aggregates only. All reads, outputs and scripts are
outside Git in `MucOneSpan-bench-data/diag12d/` (see "Reproduction").

## Owner question

"Mutations found 7/57 → 12/57 is terrible; it used to be much better in the older
open-pacmuci versions, at least on high-coverage simulated data. What happened?"
Is it a caller **regression**, a **simulation** change (MucOneUp 0.45 realistic
reads vs legacy reads), or an **evaluation** change (stricter scoring)?

## Answer

**It is not a caller regression and not a read-simulation change. The gap
comes from the design set (much harder alleles) and the metric (the old number
was an alarm rate, not exact events or a clinical call); the v0.16 gates then
turn most of the remaining mixed-allele detections into INCONCLUSIVE.**

1. **Caller version: no regression.** On identical reads, v0.16.0 is equal to or
   better than open-pacmuci v0.8.0 (the manuscript version) in per-allele exactness
   and exact event recall in every cell, and indistinguishable from v0.14.0 (the
   issue #49 "good" version). The one metric where v0.16.0 is slightly lower is
   manuscript-style alarm specificity (for example old test set 14/17 vs 16/17).
   The extra alarms are "unknown" edits at IUPAC sites; they are INCONCLUSIVE,
   not PATHOGENIC. The evidence:
   - old test set (47 samples; MucOneUp 0.44.2 legacy HiFi, ~165 reads):
     per-allele exact v0.8.0 42/94, v0.14.0 71/94, v0.15.1 69/94, **v0.16.0 71/94**
     (v0.8.0 vs v0.16.0 discordant alleles 1/30, McNemar p = 3e-8);
   - issue #49 baseline, 6 samples (sample_bench_5000, sample_close_51_58,
     sample_del_60_80, sample_dupa_100_120, sample_dupc_100_120,
     sample_normal_100_120; the list is from `MucOneSpan/results/issue49/inventory.json`):
     v0.14.0 and v0.16.0 both 9/12 alleles exact, 5/5 events, PATHOGENIC 5/5.
     The recorded baseline in `results/issue49/baseline/evaluation.json`
     (v0.14.0) was 9/12 and 5/5. Our figures are from `score.py` on
     `runs/oldtest/<caller>/scored.jsonl`, aggregated in
     `diag12d/aggregates/issue49_baseline_subset.json`;
   - manuscript experiment 1, PacBio subset (78 samples, original inputs):
     v0.8.0 77/156, **v0.16.0 111/156** (1/35, p = 1e-9); exact events 29 vs 32/39;
   - v3 clean designs, all four read sets: v0.16.0 ≥ v0.8.0 in every cell
     (e.g. legacy ONT 19 vs 9/60, p = 0.002). v0.14.0 vs v0.16.0: no cell differs
     by more than 3 alleles.
2. **Read simulation: no effect at the same depth** (for the bundle of factors
   that differ between legacy and new reads; see Limitations). Same designs, same caller
   (v0.16.0), legacy vs new clean reads: HiFi 10 vs 9/60 alleles exact, ONT 19 vs
   21/60; exact events 4 vs 6/19 and 6 vs 6/19; no paired difference is
   significant. The same holds for v0.8.0 and v0.14.0 on HiFi. On ONT, v0.8.0
   does better with the new reads (9 → 19/60, p = 0.006), because its
   legacy-ONT depth failure (item 3) disappears. The other read-mode effect is
   the manuscript-style alarm on ONT, which rises with the new reads
   (v0.16.0 15 → 18/19; v0.14.0 5 → 18/19, p = 0.0002) because the new
   calibrated R10 reads suit the R10 Clair3 model better than legacy
   QSHMM-ONT-HQ reads do.
3. **Depth: no degradation for v0.16.0; yes for v0.8.0 on legacy ONT.** v0.16.0
   per-allele exact HiFi 9 (~170 reads) → 10 (~1700) → 9/60 (new, ~1650); ONT
   21 (200) → 19 (4000) → 21/60 (new, 4000). v0.8.0 on legacy ONT drops
   16 → 9/60 from 200 to 4000 reads (7/0 discordant, p = 0.016), with gross
   (≥ 5 unit) length errors in 22/30 cases at 4000 reads. The mechanism is an
   absolute-count threshold:
   - At depth, secondary and supplementary alignments of repeat reads onto short
     ladder contigs form a spurious short cluster. Example: legacy ONT design
     dev-clean-ont_amplicon_r10-0002, 73/73, called 73/31. Its phantom allele_2
     has `reads` 232 and `primary_alignment_records` 4 in
     `diag12d/runs/leg_ont/v0.16.0/dev-clean-ont_amplicon_r10-0002/summary.json`
     (v0.14.0 is identical; v0.8.0 also reports 31 units from 232 reads).
   - This happens with legacy reads at clean depth in every caller version:
     phantom short alleles in 11–17/30 legacy ONT and 4/30 legacy HiFi cases. It
     never happens at 200 templates, and it happens in 1 case with the new reads.
     The per-case list, with called lengths, read and primary-alignment counts
     and the source `summary.json` of each case, is
     `diag12d/aggregates/phantom_alleles.json`. A phantom is a called allele
     shorter than 0.6 × the shorter truth allele.
   - v0.16.0 flags these cases (depth gate, selection unresolved: legacy ONT
     1 → 20/30) instead of reporting a wrong length, so its allele exactness does
     not drop.
4. **Design difficulty is the main effect.** The ladder resolves an allele only
   when the two alleles differ by at least 6 units and the allele is shorter than
   100 units ("easy" stratum). Within each stratum v0.16.0 performs the same
   across every data source:

   | v0.16.0, per-allele exact | easy (Δ ≥ 6, < 100 u) | long (Δ ≥ 6, ≥ 100 u) | close (Δ ≤ 5, incl. equal) | share of alleles in easy |
   | --- | --- | --- | --- | --- |
   | old test set (HiFi) | 58/63 | 5/13 | 8/18 | 67% |
   | manuscript PacBio subset | 102/110 | 3/14 | 6/32 | 71% (full set 72%) |
   | v3 clean HiFi, new reads | 8/9 | 0/7 | 1/44 | **15%** |
   | v3 clean HiFi, legacy reads | 8/9 | 1/7 | 1/44 | 15% |
   | v3 clean ONT, new reads | 13/13 | 3/3 | 5/44 | **22%** |

   The v3 dev designs are stratified by length difference (0_identical,
   0_different, 1, 2, 3–5, 6–20, >20): 22 of 30 designs have Δ ≤ 5, against
   about 20% of alleles in the manuscript set (random Vrbacka lengths) and mostly
   60/80-type pairs in the old test set. With the v3 clean HiFi stratum rates
   reweighted to the manuscript mix, per-allele exact would be ≈ 64%
   (8/9·632 + 0·64 + 1/44·184, over 880 alleles), close to the manuscript-subset
   result (111/156 = 71%).
5. **Evaluation: the metric changed, not the evaluator.** The manuscript's
   0.84 PacBio sensitivity is `parse_open_pacmuci_results`: positive if any
   allele has any mutation record. We reproduce it exactly from the recorded
   v0.8.0 outputs (PacBio 185/220 = 0.841, specificity 208/220 = 0.945;
   ONT 112/220 = 0.509, 220/220). On the same outputs, exact events are 159/220
   and per-allele exact 396/880 (45%). On the v3 clean designs, v0.16.0 reaches
   the manuscript-style sensitivity (HiFi 16/19 = 0.84, ONT 18/19) with
   specificity 10/11 and 11/11. The v0.14.0 evaluator gives the same numbers
   as the current evaluator in every cell, so scoring has not become stricter
   since v0.14.0. The stricter numbers are different metrics: exact event
   (6/19), per-allele exact (9/60) and the clinical decision (PATHOGENIC 1/19).
6. **What v0.16 changes clinically.** Gates move mixed-allele cases from alarms
   to INCONCLUSIVE. On v3 clean HiFi new reads, v0.8.0 raises an alarm in 16/19
   pathogenic and 5/11 normal cases, and 1 normal reaches PATHOGENIC. v0.16.0
   gives PATHOGENIC 1/19, INCONCLUSIVE 26/30 and FP 0/11. 9 of v0.16.0's
   16 alarms name no truth event: they are "unknown" edits at IUPAC sites (old
   evaluator: TP 7, TP_partial 9). Correct exact events blocked by a gate (event
   TP but INCONCLUSIVE): 5 cases in HiFi and 2 in ONT.

## Methods

- **Callers.** Each caller is a detached worktree at its tag with `uv sync --locked
  --all-extras`, at `MucOneSpan-oldcaller-<tag>`:
  - v0.8.0: open-pacmuci, the manuscript version; its algorithm matches the v0.3.0
    benchmark (89% HiFi sensitivity), plus ONT support.
  - v0.14.0: the issue #49 baseline.
  - v0.15.1: the last version before P0.
  - v0.16.0: this worktree.
  All four were installed and run; there was no tool drift. All used the same
  run command, `run --input <fastq> --output-dir --clair3-model --threads 3
  --platform`, with the same Clair3 binaries.
  - Models for v3 designs (bench models): ONT `r1041_e82_400bps_sup_v500`, HiFi
    `models/hifi`.
  - Models for manuscript samples: the manuscript's own, `models/ont` and
    `models/hifi`.
- **Read sets.** Truth is the same haplotype FASTA per design.
  - `new_*`: the existing v3 clean reads. HiFi: v3/clean, 2000 templates, about
    1650 reads. ONT: v3/clean2, 4000 reads, artefact-free.
  - `leg_*`: MucOneUp 0.45.0 legacy mode at the clean depth (`--coverage` =
    requested amount). Legacy mode means `reads amplicon`, no read profile, no
    `--track-read-source`, `--no-align`, and the design's `read_seed`. It uses
    config.json settings: pcr_bias `default`; HiFi pbsim3 ERRHMM-SEQUEL,
    pass 10 + CCS (about 1690 reads); ONT QSHMM-ONT-HQ, accuracy 0.95. Model
    paths are absolute, and threads are 3.
  - The legacy and new read sets differ in several factors at once, not just
    "the read layer":
    - PCR bias: legacy `default` (length-biased allele split; minor template
      share down to 0.06) vs new `no_bias`.
    - Strand: legacy is forward-strand only; new reads mix both strands
      (`forward_frac` 0.5).
    - ONT error model: legacy QSHMM-ONT-HQ vs new calibrated R10.
    - HiFi: both use pbsim3 + CCS; the new profile is not separately calibrated.
    - Per-read truth tracking (new only) does not change the reads.
  - `leg200_*`: the same at 200 templates, the old test-set condition (about
    170 HiFi reads).
  - Old test set: `MucOneSpan/tests/data/generated`, 44 HiFi and 3 ONT samples,
    MucOneUp 0.44.2 legacy, 200 templates. BAM input was converted to FASTQ with
    `samtools fastq -F 0x900`.
  - Manuscript experiment 1: MucOneUp 0.44.4, 500 templates, 220 mutant/normal
    pairs across 13 mutations. The recorded v0.8.0 outputs were scored for all
    440+440 samples. The subset was rerun with v0.8.0 and v0.16.0: 3 random pairs
    per mutation (seed 12), 78 samples per platform, from the exact `.fq` inputs
    that v0.8.0 used. The v0.8.0 rerun reproduced the recorded outputs, with
    0 discordant alleles and events.
- **Metrics.** Every cell was scored by the same current code:
  - manuscript alarm: a verbatim copy of `parse_open_pacmuci_results`, with
    sensitivity on mutant truths and specificity on normal truths;
  - the current `scripts/evaluate.py` plus benchsim `normalize_rows`: per-allele
    exact under the least favourable assignment, exact event recall and
    precision;
  - the current `compute_clinical_decision` on each caller's own `summary.json`.

  Old evaluators:
  - v0.8.0 `batch_analyze.analyze_results` (TP/TP_partial/FN/FP/TN);
  - the v0.14.0 `scripts/evaluate.py`, identical to the current evaluator in
    every cell.

  Paired tests are exact McNemar on alleles or pathogenic samples.
- **Adapters, outside the repository.**
  - Manuscript truth: MucOneUp 0.44.4 stats lack `mutation_info.mutation_targets`.
    They were derived from the same file's per-haplotype `mutation_details`.
    Normal controls (`mutation_name: normal`) were given an empty
    `mutation_info`.
  - No output adapter was needed: v0.8.0 summaries load in the current evaluator.
  - Caveat: the current gates evaluate a v0.8.0 summary only on the fields it has.
    It has no ambiguous-base count, selection, depth or phase fields, so those
    gates never fire on v0.8.0 output, and its "PATHOGENIC" column is an upper
    bound.

## 2×2 per profile (v3 clean designs, 30 per profile: 19 pathogenic, 11 normal)

Values: per-allele exact · exact event recall · PATHOGENIC on pathogenic ·
INCONCLUSIVE · PATHOGENIC on normals · manuscript alarm sensitivity / specificity.

HiFi amplicon (clean depth; legacy about 1690 reads, new about 1650):

| caller | legacy reads | new clean reads |
| --- | --- | --- |
| v0.8.0 | 8/60 · 3/19 · 4/19 · 9/30 · 1/11 · 11/19, 8/11 | 7/60 · 6/19 · 6/19 · 14/30 · 1/11 · 16/19, 6/11 |
| v0.14.0 | 10/60 · 3/19 · 2/19 · 24/30 · 0/11 · 10/19, 11/11 | 10/60 · 4/19 · 1/19 · 26/30 · 0/11 · 13/19, 10/11 |
| v0.16.0 | 10/60 · 4/19 · 2/19 · 25/30 · 0/11 · 12/19, 7/11 | 9/60 · 6/19 · 1/19 · 26/30 · 0/11 · 16/19, 10/11 |

ONT amplicon (4000 reads):

| caller | legacy reads | new clean reads |
| --- | --- | --- |
| v0.8.0 | 9/60 · 2/19 · 3/19 · 10/30 · 0/11 · 13/19, 11/11 | 19/60 · 6/19 · 5/19 · 13/30 · 2/11 · 17/19, 8/11 |
| v0.14.0 | 16/60 · 4/19 · 3/19 · 22/30 · 0/11 · 5/19, 11/11 | 21/60 · 5/19 · 5/19 · 20/30 · 0/11 · 18/19, 11/11 |
| v0.16.0 | 19/60 · 6/19 · 2/19 · 24/30 · 0/11 · 15/19, 11/11 | 21/60 · 6/19 · 4/19 · 21/30 · 0/11 · 18/19, 11/11 |

Old-depth cells (legacy, 200 templates):
- HiFi: v0.8.0 8/60 · 4/19 · 4/19 · 10/30 · 0/11; v0.16.0 9/60 · 5/19 · 0/19 ·
  29/30 · 0/11.
- ONT: v0.8.0 16/60 · 3/19 · 3/19 · 5/30 · 1/11; v0.16.0 21/60 · 6/19 · 2/19 ·
  24/30 · 0/11.

v0.15.1, the last version before P0, on new clean reads:
- HiFi: 10/60 · 6/19 · 1/19 · 26/30 · 0/11 · 16/19, 10/11.
- ONT: 20/60 · 5/19 · 5/19 · 20/30 · 0/11 · 18/19, 11/11.

v0.15.1 vs v0.16.0 discordant alleles: HiFi 1/0, ONT 0/1. The P0 release did
not change reconstruction on these designs.

Paired McNemar (discordant a-only/b-only):
- Caller axis (v0.8.0 vs v0.16.0), alleles:
  - HiFi: legacy 0/2, new 0/2, legacy-200 0/1;
  - ONT: legacy 0/10 (p = 0.002), new 3/5, legacy-200 0/5 (p = 0.06).
- Read axis (legacy vs new), v0.16.0 alleles: HiFi 1/0, ONT 3/5; v0.8.0 ONT 1/11
  (p = 0.006).
- Depth axis, alleles:
  - v0.16.0 HiFi 0/1, ONT 4/2;
  - v0.8.0 ONT 7/0 (p = 0.016).
- None of the PATHOGENIC-decision pairs is significant (smallest p = 0.06, new
  HiFi v0.8.0 5 vs v0.16.0 0 extra).

## Failure modes: present in the old caller?

The table counts cases on v3 clean HiFi new reads / legacy ONT (4000 reads).

| Failure mode | v0.8.0 | v0.14.0 | v0.16.0 | Verdict |
| --- | --- | --- | --- | --- |
| IUPAC consensus (mixed alleles in one length partition) | 28 / 26 | 15 / 4 | 24 / 18 | Present in v0.8.0 in every data set: old test set 35/47, manuscript PacBio 321/440. v0.14.0 hid it (IUPAC in 1/47 old-test cases). By the CHANGELOG this is the forced haploid consensus (`bcftools consensus -H 1`) introduced in v0.13.0; v0.13.0 itself was not run, so this attribution is inferred. Measured: v0.14.0 and v0.16.0 have the same allele exactness (old test set 71/94 both), but v0.14.0 lost ONT alarms (legacy ONT 5/19). v0.15.1 shows IUPAC again (old test set 16/47). |
| Lost/phantom allele (a short second allele built from secondary alignments, 0–4 primary; e.g. 113/113 → 113/19, 73/73 → 73/31, 121/127 → 124/32; per-case counts and source files in `diag12d/aggregates/phantom_alleles.json`) | legacy HiFi 4, legacy ONT 17, new HiFi 1, new ONT 0 | 4 / 11 / 1 / 0 | 4 / 11 / 1 / 0 | Present in all versions, and depth-driven: it never occurs at 200 templates. v0.8.0 reports it as a real allele (a gross length error); v0.16.0 marks it with the depth gate. |
| Length error (≥ 1 unit; ≥ 5 units) | 23 (1) / 26 (22) | 13 (3) / 21 (18) | 13 (3) / 18 (16) | Present in v0.8.0 and more frequent. Gross errors grow with depth on legacy ONT for every version. The genomic 142-vs-77 doubling was not re-run here: no legacy genomic arm, and v0.8.0 has no genomic mode. |
| Selection unresolved (secondary mode) | no such gate | no such gate | 3 / 20 | The gate is new in v0.16. The underlying multi-cluster state exists in all versions and produces the length errors above. |
| Per-allele depth gate (< 30 primary alignments) | min-coverage 10 only | min-coverage only | 4 / 16; legacy HiFi 9; manuscript PacBio 17/78 | The gate is new in v0.16. It fires on real read loss: Δ = 1–2 alleles keep 1–15 primary alignments of about 800 reads, and at 200 templates (about 170 reads in total) the long allele under default PCR bias keeps 9–29. |
| Correct event blocked by a gate ("no explicit sequence-level support (absent)" or "localization ambiguous") | n/a | 1–4 per set | 1–7 per set | Mostly new in v0.16 (v0.14 already had the support gate). It hits both-alleles-exact delinsAT/ins16bp cases. |
| NEGATIVE on a pathogenic truth (critical FN) | new HiFi 3, legacy HiFi 8, legacy-200 HiFi 6, legacy-200 ONT 11, legacy ONT 6 (of 19); manuscript ONT 25/39 | 0–1 per set | 0–1 per set; manuscript ONT 12/39 | Much more frequent in v0.8.0; the later gates removed most of them. The remaining v0.16.0 cases are noisy ONT reads (manuscript, R9 `ont` model) where a single-base insertion is missed while allele QC passes. |

## Candidate ladder fixes (ranked by effect on the v3 failures)

1. **Split mixed length partitions (Δ ≤ 5 and equal-length alleles).** This
   addresses 22/30 v3 designs with 1–9/44 alleles exact in every version: IUPAC
   in 18–24 cases, spurious "unknown" events (HiFi event precision 6/36), and
   phantom second alleles.
   - Cluster reads within a length partition by the heterozygous sites Clair3
     reports: haplotag, then re-call each haplotype haploid. This extends the
     experimental `read_phasing` module from phase reporting to read
     partitioning.
   - For equal lengths with no informative site, report one resolved haplotype
     plus an explicit "second allele not separable" state instead of a 0-read
     phantom.
2. **Assign reads to alleles by per-read span length, not by best ladder-contig
   primary alignment.** Use the anchor-to-anchor repeat count per read and assign
   each read to the nearest allele.
   - This targets allele-split read loss: the Δ = 1–2 minor allele keeps 1–15
     primary reads, which triggers the depth gate.
   - It also targets long alleles (≥ 100 units: 0–1/7 HiFi exact even at Δ ≥ 6)
     and ±2–5 unit length errors on long close pairs.
   - Build length candidates from primary, anchor-to-anchor spanning reads with
     a relative (fractional) support threshold instead of an absolute alignment
     record count. This removes the depth-driven phantom short alleles (legacy
     ONT 11/30 at 4000 reads, 0/30 at 200) and the concatemer sensitivity from
     Task 12c (clean vs clean2).
3. **Read-level support for dictionary-template events.** Count partition reads
   that carry the exact mutant repeat unit versus the wild-type unit at the
   localized repeat, and accept a threshold fraction as explicit support.
   - Today correct, exact calls are blocked by "no explicit sequence-level
     support (absent / localization_ambiguous)": event TP but INCONCLUSIVE in
     1–7 cases per set, including both-alleles-exact delinsAT and ins16bp.
   - The same check gives a read-level veto for the manuscript-ONT critical
     false negatives, where the consensus missed a single-base insertion.

## Limitations

- **The "simulation" arm is a bundle of factors.** Legacy vs new reads change
  the PCR bias (`default` vs `no_bias`), the strand mix (forward-only vs both)
  and, for ONT, the error model, all at once. Conclusion 2 ("no read-simulation
  effect") therefore compares two bundles of factors, not a single factor. A
  null result for the bundle does not exclude single-factor effects that cancel,
  and these factors were not varied one at a time.
- Each arm has n = 30 designs (19 pathogenic). Most clinical-decision
  differences are not significant; the allele-level caller effects are.
- ONT arms use the R10 Clair3 model for legacy (R9-like QSHMM) reads, so the ONT
  read-mode effect mixes the read model with the Clair3 model. The manuscript
  ONT subset used the manuscript's own `ont` model; both callers are poor there
  (8/156 alleles exact).
- The v0.8.0 clinical column is an upper bound; see Methods.
- The genomic arm was not re-run, so the 142-vs-77 genomic length doubling from Task 12c is untested here. Gross amplicon length errors are covered in the failure-mode table.
- Current gates were not relaxed or tuned; nothing in the caller was changed.

## Reproduction (outside Git)

`MucOneSpan-bench-data/diag12d/`:
- `scripts/`: `stage.py`, `gen_legacy.py`, `run_caller.py`, `score.py`,
  `analyze.py`, `strata.py`, `eval14.sh`.
- `legacy_config.json`.
- Drivers: `gen.sh`, `gen2.sh`, `run.sh`, `run2.sh`, `run3.sh`; logs `*.log`.
- `samples/*.json`, and `runs/<set>/<caller>/` holding the caller outputs,
  `scored.jsonl`, `scored_summary.json` and `evaluation_v014.json`.
- `analysis.json` and `analysis.md`: all cells, pairs and per-delta tables.
- `aggregates/phantom_alleles.json`: per-case phantom-allele counts with source
  `summary.json` paths.
- `aggregates/issue49_baseline_subset.json`: the 6-sample issue #49 list,
  per-caller scores and the scoring source.
