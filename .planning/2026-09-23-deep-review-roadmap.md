# Deep review and prioritized roadmap (2026-09-23)

Baseline: `main` @ 10a1440 (v0.15.1). Investigation: code review, prior-art
survey, runtime profiling, a fresh 40-case simulated truth panel, a read-centric
POA prototype, and a read-level homopolymer (dupC) evidence model. Raw
artifacts, scripts and per-agent reports live outside Git in the sibling
`MucOneSpan-review-20260923/` directory (patient-derived files stay local).

Evidence labels: **[M]** measured in this review; **[H]** hypothesis.

## 1. Executive summary

The accuracy ceiling of MucOneSpan is architectural, not a tuning problem.
Two design choices inherited from PacMUC1 cause almost every open accuracy issue:

- **A. Allele selection from alignment-record counts on a 150-contig ladder.**
  83% of idxstats records are secondary alignments on real ONT PCR [M]; gap
  clustering merges alleles into one super-cluster in 9/9 clinical amplicons
  [M]; fragments become "allele 2" while long alleles are lost (#54, #20, #21,
  #58, #42, part of #59).
- **B. Consensus = canonical-X ladder contig edited by diploid Clair3 calls +
  `bcftools consensus`.** Missed/0-1 calls silently leave X bases or IUPAC codes;
  heterozygous indels are applied as ALT by `-H I` [M]; the genotype selector
  logic drops real events (#53, #59, MP4 below; IUPAC class of errors; #48).

A template-free, read-centric prototype (anchor-to-anchor spanning reads →
span-length clustering → abPOA consensus + homopolymer vote → linked-site split
for equal lengths → existing classifier) was built and compared on identical inputs:

| Test set | Metric | v0.15.1 | Prototype |
| --- | --- | --- | --- |
| Sim panel, 40 cases (HiFi+ONT amplicon) [M] | alleles sequence-exact | 49/80 | **80/80** |
| | alleles count-exact | 73/80 | **80/80** |
| | cases fully exact | 18/40 | **40/40** |
| | mutations detected | 18/24 | **24/24** |
| | false-positive cases | 0/40 | 0/40 |
| PRJEB92208 ONT (MP1–MP4 reported dupC) [M] | dupC recovered | 2/4 (cohort-v3) | **4/4**, at published motifs |
| | HG002 PCR vs Q100 (edit distance) | 1/2 exact (3,900/900 bp) | **0 / 0** |
| | dupC in HG001–HG004, MP5 | 0 | 0 |
| Runtime per amplicon [M] | wall | 45–141 s (2 thr) | **3.9–7.6 s (1 thr)** |

Held-out check [M]: a new 40-case panel (new seeds, random lengths 25–130,
five unseen mutation types, dupC near both ends, Δ1 and equal-length designs,
low depth, 4% minor allele), prototype frozen:

| Held-out metric | v0.15.1 | Prototype |
| --- | --- | --- |
| alleles sequence-exact | 45/80 | **76/80** |
| alleles count-exact | 74/80 | **77/80** |
| mutations detected / right allele+repeat | 22/28 / 20/28 | **28/28 / 27/28** |
| false-positive cases | 1/40 | 0/40 |
| close-length / equal-length alleles exact | 0/12 / 0/8 | 11/12 / 8/8 |
| runtime median per case | 9.8 s (4 thr) | 2.4 s (1 thr) |

Prototype failures: (1) ins16bp making alleles 44 bp apart fell under the fixed
45 bp peak separation / 30 bp assignment window, so 46% of spanning reads were
silently dropped, giving a confident false homozygous call with the mutation on
the wrong allele; (2) a 5-read long allele under the 8-read peak minimum was
dropped without a warning; (3) a 51-read allele had 2 polishing errors (one at a
site only 43% of reads carried). Required design changes: length- and
platform-scaled windows (or run the split on all reads within ±1.5 units),
report rejected peaks and the unassigned spanning-read fraction ("possible
unresolved allele"), and a depth flag below ~30–60 spanning reads. Keep both
panels frozen as regression sets. Both panels come from one simulator
(MucOneUp), so real-data validation remains required.
Clinical positives are publication-reported only (three families; MP1/MP2 are
siblings) and controls are one family, so no clinical sensitivity or
specificity is claimed.

An independent read-level dupC evidence channel (C-run length per motif per
read, after splitting by span length; per-sample, per-strand stutter
background; LLR) separates all four positives from all controls by ~10× [M]:
8C fraction 0.60–0.69 at the mutant motif vs ≤0.079 at any site in
HG001–HG004/MP5; ≈30 s for the whole 11-library cohort.

## 2. Confirmed defects in the current release

| ID | Severity | Finding | Location | Evidence |
| --- | --- | --- | --- | --- |
| D1 | **S1 clinical** | MP4 (known dupC) is NEGATIVE: single-het distinct-length path builds consensus from GT allele 1 (REF) and marks `independent_haplotype_evidence=True`. Clair3 PASS `contig_71:3432 G>GC` AD 192,275 (repeat 49) is discarded; `summary.json` mutations `[]`. A unit test asserts this behaviour. | `calling.py:528-539`; `tests/unit/test_distinct_calling.py:157-160` | [M] replay GT2/IUPAC → dupC r49 |
| D2 | S1/S2 | Super-cluster: allele_1 spans contig_1..N in 9/9 clinical amplicons; allele_2 = second-largest cluster (MP3: spurious 133-unit allele). Splitters only run for a single cluster. | `alleles.py:65-79, 219-253, 490, 609-611` | [M] |
| D3 | S2 | Reported repeat count ≠ contig used (MP1 39 vs 44, MP3 37, MP5 69 vs 61). | `alleles.py:391-414` | [M] |
| D4 | S2 | `-H I` applies het indels as ALT (only SNVs become IUPAC); `variant_support` then declares all events unresolvable → `vcf_support=False` (drives #59 INCONCLUSIVE). | `variant_support.py:42-51` | [M] bcftools 1.17 |
| D5 | S2 | Haploid genotype rule uses FORMAT/AF (ALT/DP) with 0.5/0.2 cutoffs; undercounts vs AD fraction (MP4 0.448 vs 0.59). | `vcf.py:113-130` | [M] |
| D6 | S2 | HiFi dupC QUAL is erratic around the filter (15.3/19.9/2.8 at different depths); per-allele cap 500 loses dupC. | Clair3 + `vcf.py` filter | [M] |
| D7 | S3 | `--min-qual` ignored on distinct-length path (hardcoded 4.0); `haploid_*`, `min_dp` settings never read. | `vcf.py:83-87`, `calling.py:515-524` | code |
| D8 | S2 | Clinical gating order: PATHOGENIC decided before coverage/phase; NEGATIVE ignores unsplit super-clusters; coverage counts alignment records. | `report.py:92-221` | code (#55) |
| D9 | S3 | `--report-igv embedded` fails with igv-reports 1.16.0 on PATH (rejected version) → exit 1 after analysis completed. | report IGV path | [M] |
| D10 | S4 | `allele_reads.bam` reused for cluster and remapped BAMs. | `calling.py:52,139` | code |

Sim-panel failure modes of v0.15.1 [M]: close (Δ1–2) and equal-length alleles
merged into one IUPAC-laden allele (0/12 sequence-exact); reference bias at unit
positions 54–59 (A/G/E vs X); delGCCCA missed on both platforms; dupC on a
120-unit allele missed on ONT.

## 3. Prioritized work plan

### P0 — Clinical safety fixes on the current architecture (days)

1. **D1 + D4 + D5**: decide each event from per-allele AD fraction, never by
   genotype index; evaluate VCF support per event, not globally; an event with
   read support but unresolved phase must be INCONCLUSIVE-with-flag, never
   NEGATIVE. Add MP4-shaped synthetic regression; fix the test asserting D1.
2. **D8 (#55)**: gate NEGATIVE on resolved allele selection (no super-cluster,
   no unselected large cluster) and spanning-molecule depth; gate PATHOGENIC on
   event identity + frameshift + support.
3. D3, D7, D9, D10 small correctness fixes.

### P1 — Read-level homopolymer evidence channel (≈1 week)

Promote `homopolymer/hp_extract.py`/`hp_model.py` into a `homopolymer_evidence`
stage (edlib only; no external tools). Per allele×motif: reads, per-strand
counts, 8C fraction, stutter-aware LLR, depth status. Draft call rule (to be
calibrated on simulations): ≥20 allele-specific reads, LLR ≥10, fitted 8C
≥0.5, raw 8C ≥0.30, ≥3× background on both strands; 10–19 reads =
insufficient depth. Clinical decision uses agreement between this channel and
the sequence event (#59, #55). Generalize to 4C tract (MP5 motif 35 4C→5C at
0.96, outside dictionary `insC_pos23` A/E contexts) and other homopolymer
frameshifts (dupA background 3.4–4.1%).

### P2 — Read-centric reconstruction engine for amplicons (2–4 weeks)

Replace ladder selection + Clair3 + bcftools consensus for amplicon mode:

1. Anchor-to-anchor spanning read extraction (edlib; motif 1 / motif 9;
   fallback flank anchors) — promote unused `read_evidence.spanning_evidence`
   (#47). Fragments counted, never form alleles (#54).
2. Span-length peak clustering (minor peak ≥2 units away: 3% support; within
   2 units: 20%), with explicit allelic-dropout reporting (#58 PCR bias).
3. Per-allele abPOA consensus (40 random reads; random, not quality-ranked,
   sampling) + 2 majority-polish rounds + homopolymer median vote (#48).
4. Equal/close-length separation only on ≥2 linked discordant sites;
   single-site candidates reported as unconfirmed (whole-read edit-distance
   clustering fails: ~25 edits noise vs 1-edit difference) [M].
5. Existing classifier unchanged; concordance QC per allele (#46): spanning
   reads, fraction matching consensus, second-most-frequent sequence, and
   "zero PASS variants when reads are remapped to own consensus" (VNTRPipeline QC).
6. Callable gate ≈30 spanning reads per allele (MP4 long-allele dupC lost at
   ~20 reads) [M].

Dependency: `pyabpoa` (MIT; sdist on PyPI, bioconda binaries) or `pyspoa`
(MIT, pip wheels) as an extra behind the tool abstraction; `edlib` (MIT).
Keep ladder+Clair3 as `--engine ladder` fallback and genomic-mode path until
validated. Preserve CLI flags and output schemas; add fields only.

### P2b — Hybrid: assign all reads by known unit structure (measured 2026-09-23)

Degenerate repeat units (X, A, B, C, D, E, F, G, V, ...) are known prior
knowledge, so partial reads inside the VNTR are placeable. From the in-house
allele structures [M]: fragments of ≥12 units (~720 bp) have a unique position
in their allele; ≥25 units (~1.5 kb) are 100% allele-specific, even for 82/83
near-identical alleles. Ground-truth test (fragments cut from spanning reads,
allele known by span; edlib infix edit distance to ladder-flanked POA drafts,
assign if margin ≥3): **3 wrong / ~4,300 decided**; at ≥1.5 kb, 55–100% decided.
Genomic ONT in-house has 57–127 reads with ≥1.5 kb of VNTR vs 14–60 spanning
(LB-level details local only), i.e. ~2–5× usable depth. minimap2 MAPQ discards
most of this for near-identical alleles (70/80 reads MAPQ<20 in the 82/83 case).

Design: ladder → length prior + flanks; spanning reads → POA drafts;
all reads → edit-distance competition against ladder-flanked drafts (undecided
kept separate) → polish + per-motif homopolymer evidence with assigned reads →
self-reference remap QC (Clair3/pileup expects zero variants) → ladder-seeded
iterative polishing when <~10 spanning reads. Validate assignment on simulated
panels with known read origin before use.

### P3 — Speed of the legacy path (if retained)

Measured: minimap2 is ~93% of CPU (ladder map 40–57% of wall, remap 25–35%);
Clair3 has 4–9 s fixed per allele. Options: seeded read cap before ladder
mapping (HG004 55–77 s → ~12–18 s; ONT outputs identical down to ~2k reads)
[M]; split threads across parallel remap workers; stream `minimap2 | samtools
sort` (currently buffered in Python strings); single-pass allele stats; Clair3
`--no_phasing_for_fa`. Do **not** use `--pileup_only` (changes genotypes) [M].
Per-allele caps need dupC re-validation (D6). P2 makes most of this moot.

### P4 — Validation and generalization

- Held-out frozen-parameter simulated panel (running at review time).
- R10-like ONT homopolymer stutter in simulation (pbsim3 ONT models used here
  do not reproduce R10.4.1 7C/8C stutter; Clair3 sim model r941 vs clinical r1041).
- Genomic ONT (12–45 spanning reads): evaluate P2 at low depth (#58 genomic mode).
- Research-only comparators: Medaka Tandem (ONT; best published exact-sequence
  and homopolymer accuracy, but ONT research licence), LongTR (raise max repeat
  length), vamos, TRGT (HiFi, PacBio licence). Finish #56 VNTRPipeline run.
- Record basecaller model/simplex-duplex provenance in reports.
- ERR15277553 (unresolved WGS identity): alleles identical to MP1 (44/78 with
  78:17 dupC from 4 reads) — supports MP1, low depth; not a call.

- In-house genomic ONT (five samples, R10-era, local only; IDs kept out of Git)
  [M]: only 14–60 of 560–1,429 reads span motif 1 to motif 9. v0.15.1 produced
  fragment "alleles" (18–19 units) in 2/5 samples, the #54 pattern. The earlier
  open-pacmuci run called a likely-false frameshift on a 32-unit fragment
  allele. One sample carries dupC at motif 33 of an 82-unit allele that is only
  one unit shorter than the other (83). Per read, the dupC segregates perfectly
  with span length (7/7 vs 0/7 reads), yet the length-mode model merged the
  alleles. At about 6 reads per allele, the POA prototype also produced a
  spurious second dupC on the other allele. Genomic mode therefore needs a
  depth gate, per-read phasing of close lengths, and the unresolved-allele
  report of §3 P2.

### Deprioritized / not recommended

Canu (minutes, collapses same-length heterozygosity); Longshot (SNV only);
DeepConsensus (needs subreads); WhatsHap read phasing default (earlier false
12:1 split); `-H A` genome-wide (more false frameshifts); lowering Clair3
thresholds; #42/#45 policy/geometry refactors before the engine change.

## 4. Issue mapping

| Issue | Root cause | Addressed by |
| --- | --- | --- |
| #54, #20, #21, #47, #58 | A | P2.1–2.2 |
| #53, #59, #48 | B (+ D1, D4, D5) | P0.1, P1, P2.3 |
| #55 | D8 | P0.2, P1 |
| #46 | missing QC | P2.5 |
| #56, #44 | validation | P4 |
| (new) D1 MP4 false negative, D3, D6, D7, D9, D10 | — | P0 |

## 5. Open validation before implementation

- Held-out simulated panel with frozen prototype (gating P2 claims).
- ROC for the homopolymer LLR on simulations (dupC at random motifs, either
  allele, varied depth), plus HiFi.
- Do not restate the prototype as clinical sensitivity; positives are
  publication-reported, three families.
