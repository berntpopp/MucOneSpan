# Validation Evidence

This guide documents empirical validation across historical cohorts and the comprehensive 200-dataset simulation benchmark.

## Validation evidence for 0.16.1

On the diag12d `ms_ont_sub` benchmark set (78 older simulated ONT amplicon
samples called with the R9 Clair3 model: 39 pathogenic, 39 normal), the
v0.16.1 caller-stage discordance gate reduces NEGATIVE on pathogenic samples
from **12/39 to 1/39**. NEGATIVE on normal samples is unchanged at 25/39, and
no PATHOGENIC decision and no other benchmark set (PRJEB92208 cohort-p0,
v3/standard, v3/clean, v3/clean2, v2 stress, ms_pacbio_sub, oldtest, leg_ont,
leg_hifi, leg200_ont, leg200_hifi) changed. The one residual NEGATIVE on this
set (`pair_5178`) is an accepted, documented limitation; see
[limitations](../reference/limitations.md).

## Validation evidence for 0.14.0

MucOneSpan 0.14.0 was evaluated on the 500-dataset stratified simulation benchmark (250 distinct diploid biological designs × 2 sequencing modes: PacBio HiFi amplicons and genomic Oxford Nanopore). The cohort was partitioned under strict pre-registration discipline into Development (300 datasets), Validation (100 datasets), and held-out Test (100 datasets), governed by cryptographic ledger seals (`tests/data/experiment_500/ledger_sealed.jsonl`).

### Multi-Split Benchmark Summary

| Cohort / Split | Platform | Samples | Diploid Sequence Exact | Diploid Count Exact | Count $\pm 1$ Unit | Supported Event Precision | Supported Alarm Precision | Normal Control Specificity |
|---|---|---|---|---|---|---|---|---|
| **DEV** (300) | Amplicon HiFi | 150 | 42 (28.0%) | 96 (64.0%) | 111 (74.0%) | 95.2% (40/42) | 97.6% (41/42) | 95.0% (19/20) |
| | Genomic ONT | 150 | 0 (0.0%) | 1 (0.7%) | 1 (0.7%) | 66.7% (6/9) | 100.0% (10/10) | 100.0% (18/18) |
| | **DEV Pooled** | 300 | 42 (14.0%) | 97 (32.3%) | 112 (37.3%) | 90.2% (46/51) | 98.1% (51/52) | 97.4% (37/38) |
| **VAL** (100) | Amplicon HiFi | 50 | 15 (30.0%) | 34 (68.0%) | 38 (76.0%) | 100.0% (11/11) | 100.0% (11/11) | 100.0% (6/6) |
| | Genomic ONT | 50 | 1 (2.0%) | 4 (8.0%) | 7 (14.0%) | 100.0% (6/6) | 100.0% (6/6) | 100.0% (7/7) |
| | **VAL Pooled** | 100 | 16 (16.0%) | 38 (38.0%) | 45 (45.0%) | 100.0% (17/17) | 100.0% (17/17) | 100.0% (13/13) |
| **TEST** (100)| Amplicon HiFi | 50 | 11 (22.0%) | 36 (72.0%) | 37 (74.0%) | 91.7% (11/12) | 100.0% (12/12) | 100.0% (4/4) |
| | Genomic ONT | 50 | 1 (2.0%) | 1 (2.0%) | 3 (6.0%) | 75.0% (3/4) | 100.0% (4/4) | 100.0% (4/4) |
| | **TEST Pooled**| 100 | 12 (12.0%) | 37 (37.0%) | 40 (40.0%) | 87.5% (14/16) | 100.0% (16/16) | 100.0% (8/8) |
| **ALL POOLED**| Amplicon HiFi | 250 | 68 (27.2%) | 166 (66.4%) | 186 (74.4%) | 95.4% (62/65) | 98.5% (64/65) | 96.7% (29/30) |
| | Genomic ONT | 250 | 2 (0.8%) | 6 (2.4%) | 11 (4.4%) | 78.9% (15/19) | 100.0% (20/20) | 100.0% (29/29) |
| | **TOTAL** | **500** | **70 (14.0%)** | **172 (34.4%)**| **197 (39.4%)**| **91.7% (77/84)** | **98.8% (84/85)** | **98.3% (58/59)** |

### Scientific and Clinical Takeaways

1. **Perfect Alarm Precision on Held-Out Validation & Test Data:**
   Across both the protected validation set (VAL-100) and held-out test set (TEST-100), Supported Sample Alarm Precision reached **100.0%** (17/17 in VAL, 16/16 in TEST). Zero false alarms were triggered on normal controls in either split (Normal Control Specificity = 100.0%).
2. **PacBio HiFi Robustness Across Allelic Length Spans:**
   On PacBio HiFi datasets, diploid repeat count accuracy consistently achieved **64.0% – 72.0% exact** and **74.0% – 76.0% within $\pm 1$ repeat**, demonstrating that the bimodal nearest-peak cluster assignment in `length_candidates.py` completely eliminated the historical $+1$ repeat shift for $\Delta \ge 2$ alleles.
3. **Genomic ONT Performance Profile:**
   While genomic ONT achieved 100% normal control specificity and 100% alarm precision, full diploid sequence recovery was constrained by the 7.5 kb simulated read length ceiling in NanoSim, which precludes span-spanning reads for long VNTR alleles in non-amplicon genomic contexts.
4. **Accessible Clinical Bio-UX Verification:**
   The clinical reporting interface was evaluated with Playwright across headless Chromium, WebKit, and Firefox. Reports rendered with zero external network requests (fully offline air-gapped compliance), 100% WCAG 2.1 AA/AAA contrast ratios, distinct SVG status iconography, and prominent allele multiplicity caveats.

## Validation evidence for 0.12.0

MucOneSpan 0.12.0 was evaluated on 200 newly simulated sequencing datasets generated with MucOneUp across 100 diploid biological designs (140 development datasets, 60 protected final validation datasets) evaluated across PacBio HiFi and Oxford Nanopore amplicons.

### Benchmark Results (200 Simulation Datasets)

| Metric | Baseline (v0.11.0) | Candidate (v0.12.0) | Net Delta / Outcome | Acceptance Gate |
|---|---|---|---|---|
| **HiFi Diploid Exact Sequence** | 13 / 100 (13.0%) | **33 / 100 (33.0%)** | **+20 wins, 0 losses** | $\ge +8$ (PASS) |
| **ONT Diploid Exact Sequence** | 2 / 100 (2.0%) | **14 / 100 (14.0%)** | **+12 wins, 0 losses** | $\ge +5$ (PASS) |
| **Total Exact Sample Reconstruction** | 15 / 200 (7.5%) | **47 / 200 (23.5%)** | **+32 wins, 0 losses (3.13×)** | Wins > Losses (PASS) |
| **Normal Control Specificity** | 35.7% (dev) / 40.0% (final) | **84.6%** (dev) / **100.0%** (final) | **0 new false alarms** | $\le 1$ FP (PASS) |
| **False-Positive Mutation Calls** | 96 calls | **27 calls** | **-71.9% reduction** | Strict reduction (PASS) |
| **Per-Allele Exact Repeat Counts** | 217 / 400 (54.2%) | **272 / 400 (68.0%)** | **+55 alleles (+13.8 pp)** | Zero regressions (PASS) |

### Length Error Profile (How Far Off Are Non-Exact Calls?)

Across all 400 alleles in the 200-case stress panel (which deliberately over-indexes on difficult cases: 1-repeat separation, identical lengths, and extreme PCR asymmetry):

* **Exact repeat count (`count_exact`):** **68.0%** (272 / 400 alleles) — 64.0% HiFi, 72.0% ONT.
* **Within $\pm 1$ repeat (`count_within1`):** **72.5%** (290 / 400 alleles).
* **Within $\pm 2$ repeats (`count_within2`):** **73.0%** (292 / 400 alleles).
* **Gross mismatch ($> 2$ repeats, ~27.0%):** Primarily driven by amplicon PCR length dropout in extreme asymmetric alleles (e.g. 25 vs 140 repeats at 60 templates). Because the pipeline requires $\ge 3$ (HiFi) or $\ge 4$ (ONT) dominant spanning reads to confirm an allele, the long allele (yielding 0–1 reads) is conservatively withheld. Uncalled alleles are scored as missing predictions, creating an apparent error equal to the full allele gap (e.g. 115 repeats).

### Root Causes of Sequence vs. Length Discordance

In about half of all evaluated alleles (48.9% in the standard 44-sample HiFi benchmark), the allele length is **100% exact (0 bp error)**, yet strict sequence identity fails (`sequence_exact: false`). The sequence edit distance is typically only **1 to 3 base pairs** across the 3,000–6,000 bp array.

Empirical investigation confirmed:

1. **IUPAC Ambiguity Codes from Clair3 `0/1` Calls (95.3% of mismatches):**
   The reference ladder contig consists of canonical `X` repeat units. Biological alleles contain variant units (`A`, `B`, `C`, etc.) differing by 1–2 SNPs. Due to repetitive alignment jitter or sequencing errors, a small fraction of reads (~10–25%) carry reference bases. Clair3 frequently classifies these sites as heterozygous `0/1` rather than homozygous `1/1`. `bcftools consensus` applies IUPAC ambiguity codes (`S` for C/G, `M` for A/C, `R` for A/G). While the sequence length remains exact, IUPAC symbols count as mismatches in strict ACGT sequence comparisons and prevent `classify.py` from matching pure ACGT dictionary units (rendering them as `?`).
2. **Complete Catalogue Coverage (Zero Missing Units):**
   The repeat dictionary in `repeats.json` contains all 34 known biological units (fixed 1–9, canonical X, and variants A through W). Zero truth repeat units are missing from the catalogue.
3. **MucOneUp Biological Fidelity:**
   MucOneUp generates haplotypes by chaining authentic repeat units from the Vrbacka model; it does not introduce synthetic background SNPs into repeat units.
4. **Clair3 False Negatives (Reference Fill):**
   At low coverage, Clair3 occasionally misses a variant SNP, leaving the reference ladder's canonical `X` unit in place of the variant unit (1–2 bp mismatch).

---

## Validation evidence for 0.11.0

## What was verified

| Change | Evidence | Interpretation |
| --- | --- | --- |
| Exact bit-vector edit distance | At least 10,000 independent scalar-oracle comparisons; 141 cached complete classification outputs identical | Scientifically equivalent scoring optimization |
| Classification performance | Three alternating paired difficult-panel repetitions: median 17.612304 s to 2.214009 s | 7.955× local classification speedup, not a pipeline speedup |
| Default configuration integration | Five actual before/after cases, ten successful full pipeline runs across HiFi and ONT | Sequences, structures, counts and mutation records unchanged; provenance added |
| Corrected boundary anchors | 137 available cached full consensuses unchanged; 17 missing allele slots accounted for | Fixes anchor metadata; a synthetic flank insertion demonstrates the intended trimming correction |
| Signed indels and variant support | Deterministic mixed-indel tests and real bcftools normalization/consensus fixtures | Frame uses net inserted minus deleted bases; support requires exact projected identity |
| Failure and ambiguity accounting | Adversarial evaluation tests for wrong identity/position/parent, tied assignments, extra/missing alleles, stale files and failed commands | Failed or ambiguous samples cannot become successful negatives |

The configuration comparison is a functional regression check. Its index
warm/cold conditions were not controlled for a performance claim. Whole-pipeline
speed improvement has not been established. Confidence values remain heuristic
scores, not calibrated probabilities. VCF concordance is agreement with the VCF
used for consensus, not independent read validation.

## Development accuracy

All 44 original HiFi samples, three ONT samples, six subsequently exposed HiFi
challenges and 24 dependent depth/name perturbations are development evidence.
The six challenges are no longer held out. The following table re-scores archived
before/after outputs with one strict evaluator; it does not imply every current
source revision was rerun over every sample.

| Cohort | Inputs | Exact individual sequences, before → after | Exact diploid sequences | Exact mutation annotations, before → after | Extra mutation records, before → after |
| --- | ---: | --- | --- | --- | --- |
| Original HiFi | 44 | 39/88 → 39/88 | 10/44 → 10/44 | 24/28 → 24/28 | 14 → 14 |
| Original ONT | 3 | 3/6 → 3/6 | 0/3 → 0/3 | 2/2 → 2/2 | 0 → 0 |
| Later HiFi challenges | 6 | 1/12 → 1/12 | 0/6 → 0/6 | 2/3 → 2/3 | 0 → 0 |
| Dependent perturbations | 24 | 6/48 → 7/48 | 0/24 → 0/24 | 6/14 → 7/14 | 12 → 11 |

The perturbation improvement includes repairing an invalid legacy alias artifact;
it is not proof of a new reconstruction algorithm. The six low-depth no-calls
remain in denominators. Matching is one-to-one and retains assignment ambiguity;
an unproven duplicate sequence receives no independent second-allele credit.

Historically, 15/16 HiFi normal samples had no mutation alarm. Ten have unresolved
reconstruction under the stricter contract. This leaves five confident negatives,
one false-alarm control and ten unresolved negatives. Conditional specificity is
5/6, with the unresolved fraction reported separately. The single ONT normal
control is unresolved, so its conditional specificity is undefined. These small,
simulation-based cohorts do not establish clinical sensitivity or specificity.

## Experiments that were not promoted

| Experiment | Observed result | Decision |
| --- | --- | --- |
| Minimum coverage 10 → 5 | More counts and outputs, but extra events 25 → 38 and normal false-alarm samples 4/30 → 7/30 across 77 inputs | Retain default 10 |
| Optional read-backed phasing | Cached exact individual sequences 1/10 → 4/10, but an extra false event; no exact diploid recovery | Experimental opt-in, disabled by default |
| Minimum selected-phase support | Can remove a false singleton but also reject real minority alleles; original 30:2 became internally selected 14:1 | No universal support floor |
| Strict segmentation stop | Lost five named events among 106 exposed consensus outputs | Experimental only; retain default recovery and expose unresolved sequence |
| Simple anchor length modes | Promising wider-gap counts, but bin suppression merges true one-repeat-separated alleles; many ONT reads lack usable anchors | Diagnostic only |
| Continuous-span one/two component models | Improved count-only results but still split the equal-length control into 59/60; downstream mutation recovery untested | Not promoted |

The equal-length failure has strong reference-fit artifacts: the false longer
candidate has many alignment records but little full-boundary read support. The
25/140 case loses all five long primary records before calling because their
reference-fit counts are spread below the per-contig threshold; its false second
short component is fed one primary read. These explain the failures without
establishing a safe replacement inference/assignment method. The tests retain
both strict expected failures and their original tolerances.

## Reproducing and extending the evidence

Use the [simulation experiment guide](simulation-experiments.md) to specify new
HiFi and ONT cases, seeds, lengths, mutation targets and requested templates in
JSON. Keep simulator and platform model configurations explicit. Requested
coverage is not the number of retained reads: record actual usable inputs.

Freeze caller settings, tool/model versions, source revision and evaluation rules
before creating a final panel. Keep generation truth out of calling logic. If you
change settings after examining final results, classify that panel as development
and use new seeds for the next final evaluation. Report sample failures, no-calls,
assignment ambiguity, extra calls and per-sample regressions alongside aggregates.
Perturbations from the same source sample are dependent observations.

The repository's `.planning/` evidence reports retain experiment decisions,
commands, input/source hashes and reviewer dispositions. Raw reads, results,
indexes and model files remain local generated artifacts. See also
[benchmarking](benchmarking.md), [configuration](configuration.md), and
[limitations](../reference/limitations.md).
