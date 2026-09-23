# Spec: realistic simulated benchmark (MucSim-Bench)

Status: draft for review, 2026-09-23. Companion: `2026-09-23-hybrid-engine-spec.md`.
Evidence: `../MucOneSpan-review-20260923/{realprofile,simpanel,heldout}` (local only).

## 1. Problem

Previous simulated panels (the 500-dataset set, the 40-case simpanel, the 40-case
held-out panel) came from MucOneUp defaults. They differ from real data in ways
that change caller behaviour [measured]:

| Property | Real ONT R10 amplicon (PRJEB92208, 9 libs) | MucOneUp default ONT amplicon |
| --- | --- | --- |
| Read strands | both | one strand only (600/600 reads flag 16) |
| Error rate in VNTR | 2.1% median (mm 0.70 / ins 0.67 / del 0.77) | 3.75%, insertion-heavy (44/30/26) |
| C7 read correctly, + / − strand | 0.52 / 0.89 | ~0.85, symmetric |
| Short single-junction "smear" products | 8–52% of spanning reads (median 24%) | none |
| Between-allele chimeras | ~2.3% | none |
| Allelic ratio vs length | ln(long/short) = −0.056 × Δunits (R² 0.85) | model exists, uncalibrated |
| Off-target / short reads | 6–77% (median 30%) | none |
| Read→allele truth | n/a | absent; read names collide |

Genomic ONT (in-house, 5 samples) instead has 85–87% off-VNTR reads in a
100 kb extract, 1.8–9.1% spanning (14–60 reads), 4–7 kb median reads,
long-read-driven allele ratio (~0.75), ~2% smear, and the same strand asymmetry
(C7 correct 0.50 / 0.93).

Because of these gaps, simulated wins may not carry over to real data, and phantom-allele /
dupC-stutter failure modes are under-tested. Current panels are also too small
for statistically meaningful comparisons.

## 2. Goals

1. Generate simulated datasets for three **profiles** whose measurable
   properties match real data within stated tolerances:
   - `ont_amplicon_r10`: PRJEB92208-like PCR amplicons.
   - `ont_genomic_targeted`: in-house-like genomic/targeted R10 reads.
   - `hifi_amplicon`: PacBio HiFi amplicons (no real calibration data; marked
     uncalibrated, tuned to published HiFi ≥Q30).
2. Complete per-read truth: allele of origin, molecule id, strand, product
   kind (full / smear / chimera / concatemer / off-target / fragment), and the
   injected homopolymer edits.
3. A large, stratified, sealed benchmark supporting paired engine comparisons
   with pre-specified statistics.
4. Reproducible from `(design, seed, tool versions, calibration file)`, and
   fast enough to regenerate: target ≤5 s per ONT case and ≤30 s per HiFi case.

Non-goals: basecalling simulation from raw signal; Illumina; somatic mosaicism
beyond an optional two-component mixture; methylation.

## 3. Architecture

```
design row (factors, seed)
  └─ haplotypes: MucOneUp `simulate` (--input-structure or --fixed-lengths,
     --mutation-name/targets, strict mode) → hap1.fa, hap2.fa, structure, stats
  └─ profile generator
       ├─ amplicon profiles: TemplateBuilder
       │     primer-to-primer amplicons → molecule pool:
       │     PCR split (MucOneUp PCRBiasModel, calibrated slope)
       │     strand 50/50, smear (single-junction deletion products),
       │     between-allele chimeras, concatemers, off-target shorts,
       │     per-molecule homopolymer stutter (calibrated table, strand-aware)
       │   → pbsim3 --strategy templ (QSHMM-ONT-HQ | ERRHMM-SEQUEL+ccs),
       │     --id-prefix per molecule, --difference-ratio/--accuracy-mean fitted
       └─ genomic profile: OPEN DECISION (owner). Options: (1) extend MucOneUp
             upstream (per-read truth, strand mix, calibrated stutter, smear/
             chimera, flank length, read-length distribution, skip hg38
             alignment) and call it; (2) MucOneUp NanoSim mode as-is with
             downstream name parsing and post-hoc stutter; (3) local fragment
             sampler + direct pbsim3 templ. Target properties: in-house median
             read length 3.5–6.9 kb, 50/50 strand, ~1:1 alleles, 1.8–9.1%
             spanning reads in a 100 kb extract.
  └─ outputs: reads.fastq.gz, truth/ (haplotypes, structures, mutations with
     1-based repeat index and bp), read_truth.tsv.gz, case.json (factors,
     seeds, versions, calibration hash)
```

MucOneUp is used for haplotype and truth generation only. The read layer is
our own code, because MucOneUp amplicon reads lack strand, truth, chimeras,
smear and calibrated errors. We call pbsim3 directly (as MucOneUp does) and
avoid MucOneUp's hg38 realignment (~12 GB RAM, most of the runtime).

Code location: `src/muc_one_span/benchsim/` (tested library) +
`scripts/benchsim.py` (CLI). External tools (`muconeup`, `pbsim`, `ccs`,
NanoSim `simulator.py`) are run through `tools.run_tool` with argument lists.
Datasets live outside Git in a sibling directory `../MucOneSpan-bench-data/`.

## 4. Calibration

- `benchsim/calibration/*.json`: versioned calibration tables with provenance
  and SHA-256. Only **aggregate** statistics are stored; never reads or
  sequences.
  - `ont_r10_amplicon_v1.json`: derived from public PRJEB92208. May be committed.
  - `ont_r10_genomic_v1.json`: derived from in-house aggregates. Commit only
    with explicit owner approval; otherwise kept local and referenced by hash.
  - `hifi_amplicon_v0.json`: literature defaults, flagged `uncalibrated`.
- Parameters: strand ratio; PCR slope; smear rate and junction-position
  distribution (relative position ~Beta, median 0.2–0.3); chimera and
  concatemer rates; off-target short-read fraction and length distribution;
  homopolymer stutter table P(injected length | true length, base, strand)
  for runs ≥3; pbsim accuracy mean/sd and difference ratio; genomic read
  length median/sd and spanning fraction target.
- **Stutter deconvolution.** The measured table P(obs | true) includes pbsim's
  own homopolymer errors. The injected table is fitted iteratively
  (simulate → measure → adjust) until simulated P(obs | true) for C/G runs
  of 6–8 matches the target within ±0.03 absolute per cell (strand-specific).
- **Realism acceptance** (`benchsim realism-report`): simulate matched designs
  (the real samples' reconstructed structures: HG002 Q100 public; others
  from local consensus, kept local) and compare with real data:
  - error rates by type and strand: within ±20% relative;
  - C7 correct-length fraction per strand: ±0.03;
  - smear fraction and between-allele fraction: within the real
    library range, median within ±25% relative;
  - span-offset pmf (15 bp bins): Jensen–Shannon distance ≤0.1;
  - allele-ratio slope: ±0.01;
  - genomic spanning-read fraction: within the real range.
  Also run both engines on real and matched-simulated data. Per-metric
  sim-vs-real differences (e.g. per-allele exact rate, dupC 8C fraction) are
  reported as the **sim-to-real gap**. Gaps are a benchmark limitation to
  state, not something to hide.

## 5. Experimental design

Factors (stratified; levels sampled by Latin hypercube within strata):

| Factor | Levels |
| --- | --- |
| Profile | ont_amplicon_r10, ont_genomic_targeted, hifi_amplicon |
| Allele A length (units) | 20–130 (uniform over MucOneUp model support) |
| Δ length | 0 identical, 0 different sequence, 1, 2, 3–5, 6–20, >20 |
| Composition | MucOneUp Markov; 20% real-derived structures; 5% with rare/novel units |
| Event | none (normal); dupC; other dictionary frameshifts (insG/insA pos54/58, delGCCCA, delinsAT, dupA, insCCCC, ins25bp, del18_31, ins16bp, insC_pos23); benign in-frame (insCCC/insCCCCCC/delCCC/del3bp benign); combos (event on both alleles, 2 events on one allele) |
| Event allele | shorter, longer, equal-length partner |
| Event position | first 10%, middle, last 10% of repeats |
| Depth | amplicon spanning reads/allele: 5, 10, 20, 30, 60, 150, 500, 2000; genomic: 3, 6, 10, 20, 40, 80 |
| PCR bias | calibrated, strong (2× slope), none |
| Artefact level | smear 0.05 / 0.25 / 0.5; chimera 0.01 / 0.05 |
| Error level | calibrated; +50% ("poor run", Q13-like) |

Splits (design-level, so no haplotype design is shared across splits):

| Split | Cases | Use |
| --- | --- | --- |
| `dev` | 900 (300 / profile) | development; free use |
| `val` | 900 | threshold selection; limited looks, each logged in the ledger |
| `test` | 2,400 (800 / profile) | sealed; one evaluation per release candidate |
| `stress` | ~300 | targeted failure families from §5 edges and past failures |

At least 35% of every split are normals (≥280 per profile in `test`), so zero
false positives give a 95% upper bound ≤1.1% per profile. The
held-out and simpanel sets are kept as small frozen regression sets.

Sample-size rationale (exact tests, simulated power): paired McNemar with
80% power detects +10 pp (vs 2 pp regressions) with ~200 paired cases,
+5 pp with ~400, +3 pp with ~1,200. All-success Clopper–Pearson 95% lower bounds:
n=59 → 0.95, n=300 → 0.99.

Sealing: generate `test` with seeds derived from a secret salt. Record the
design and truth hashes in the existing durable ledger
(`durable_ledger.py`, pattern of `tests/data/experiment_500/ledger_*.jsonl`).
Truth files are encrypted or stored outside the working tree until evaluation.

## 6. Metrics and statistics

Primary, per profile:

1. Per-allele exact sequence (motif 1–motif 9, ACGT only).
2. Case exact (both alleles, with correct zygosity).
3. Event recall and precision on (allele, 1-based repeat index, name); dupC
   reported separately.
4. Clinical decision confusion matrix: truth pathogenic / benign / normal
   against PATHOGENIC / INCONCLUSIVE / NO_PATHOGENIC. False NO_PATHOGENIC on a
   pathogenic truth is the critical error; INCONCLUSIVE rate is reported.
5. False-positive pathogenic calls on normals and on benign events.

Secondary: repeat-count exact and ±1; edit distance; read-assignment
accuracy (hybrid); calibration of reported confidence (reliability curve,
Brier score); runtime and peak RSS.

Statistics:

- Clopper–Pearson 95% CIs per stratum.
- Paired engine comparison: exact McNemar per profile, with the Holm
  correction across primary metrics.
- Pre-registered **non-inferiority** for normal FP rate: margin 0.5 pp,
  one-sided 95%.
- Cluster bootstrap by haplotype design for pooled estimates.
- Failure-mode model: logistic regression of failure on factors, reported with
  odds ratios. Also report a per-factor failure atlas.

Decision rule for adopting an engine or threshold (pre-registered in the
ledger before unsealing `test`):

- superior on metric 1 and non-inferior on metric 5 for every profile; and
- no increase in false NO_PATHOGENIC.

## 7. Outputs and interfaces

- `case.json`: design factors, seeds, tool versions, calibration hash.
- `read_truth.tsv.gz`, columns: `read_id, allele, molecule, strand, kind,
  src_start, src_end, stutter_edits`.
- `truth/`: MucOneUp outputs, compatible with `evaluation/truth.load_truth`.
- `manifest.jsonl` per split.
- Scoring: extend `evaluation/` (not a new scorer) with the clinical-decision
  confusion matrix, read-assignment accuracy, and per-stratum aggregation;
  `scripts/evaluate.py --manifest`.
- Reports: `benchsim report` renders stratified tables and McNemar results
  (JSON + Markdown).

## 8. Compute and storage

| Profile | Generation time per case | Stored per case |
| --- | --- | --- |
| ONT amplicon | ~1–3 s (pbsim, no hg38) | amplicon reads capped by depth factor; ≤10–25 MB gz |
| Genomic | ~5 s | small |
| HiFi | ~25 s (ccs) | small |

About 4,500 cases take ~3–6 CPU-hours and ~40–80 GB outside Git. Runs are
deterministic per seed. `dev` can be regenerated rather than stored.

## 9. Risks

- **Simulator realism ceiling:** pbsim QSHMM is context-free apart from our
  injected homopolymer edits. Mitigation: the realism report and the
  sim-to-real gap are published with every result.
- **Overfitting to the calibration libraries:** 9 public + 5 in-house samples.
  Mitigation: leave-one-library-out realism checks, and include an "error +50%"
  level.
- **HiFi is uncalibrated:** keep HiFi conclusions separate and flagged.
- **MucOneUp quirks** (non-strict mutation rewrite, hard-coded 50/50 6 vs 6p,
  null provenance seed): force strict mode, record structures, and verify each
  planned mutation exists in the output.
