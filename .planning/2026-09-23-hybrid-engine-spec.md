# Spec: hybrid read-centric reconstruction engine

Status: draft for review, 2026-09-23. Companion: `2026-09-23-realistic-benchmark-spec.md`.
Background: `2026-09-23-deep-review-roadmap.md` (§2 defects, §3 P1/P2/P2b).
Prototype (outside Git): `../MucOneSpan-review-20260923/poa-prototype/{proto.py,hetsplit.py}`,
`homopolymer/{hp_extract.py,hp_model.py}`, `inhouse/hybrid/{hybrid.py,assign_test.py}`.

## 1. Problem and evidence

The ladder engine selects alleles from alignment-record counts and edits a
canonical-X contig with diploid Clair3 calls. This causes phantom or merged
alleles (#54, #20), IUPAC/0-1 artefacts and dropped events (#53, #59, MP4
false NEGATIVE), and ~93% of CPU time in minimap2.

Prototype evidence [measured]:

| Test | Ladder engine v0.15.1 | Prototype |
| --- | --- | --- |
| Dev sim panel, alleles sequence-exact | 49/80 | 80/80 |
| Held-out sim panel, alleles sequence-exact | 45/80 | 76/80 |
| Held-out sim panel, mutations detected | 22/28 | 28/28 |
| PRJEB92208 dupC positives | 2/4 | 4/4 |
| HG002 vs Q100 | 1/2 alleles | 2/2 exact |
| Runtime per case | 10–140 s | 1–8 s |

The prototype's held-out failures come from a fixed 45 bp peak window
(close alleles merged), silent peak rejection, and low-depth polishing errors.
Genomic ONT gives only 14–60 spanning reads. Assigning all reads by
edit-distance competition against the ladder-flanked allele drafts gave 3
wrong out of ~4,300 decided, with 2–5× more usable reads.

## 2. Goals and non-goals

Goals:

1. `--engine hybrid` for `muconespan run` and the benchmark and clinical
   harnesses, producing the existing output contract. New fields may be
   added; nothing is removed.
2. Exact allele sequence reconstruction from amplicon and genomic long reads,
   including equal-length and Δ1-unit heterozygotes where reads show linked
   differences.
3. Explicit evidence and uncertainty: rejected length peaks, unassigned read
   fraction, depth status, and per-event read-level support. Missing evidence
   is never reported as negative.
4. Speed: ≤10 s per amplicon sample and ≤5 s per genomic sample on one thread,
   with no Clair3 in the default path.
5. The ladder engine stays available and stays the default until the
   benchmark decision rule (benchmark spec §6) is met on the sealed test split.

Non-goals: basecalling; replacing classification or dictionary nomenclature;
mosaicism calling (flagged only); WGS whole-genome processing (input is
reads or a region BAM).

## 3. Algorithm

```
reads (FASTQ, or BAM → primary reads)
 S1 anchor + categorize        edlib HW search for motif-1 and motif-9 anchors (≤k edits,
                               both strands), with flank anchors as fallback when a motif
                               anchor is mutated.
                               Categories: spanning / left-anchored / right-anchored /
                               internal / off-target. Orient to motif1→motif9.
 S2 length model               spanning span lengths → peaks. Separation window and
                               assignment window scale with length and platform:
                               w(L) = max(w0, c·sd_platform(L)); sd measured ~0.25 bp/unit + base.
                               Minor peak thresholds: ≥max(n_min, f_far·N) when ≥2 units
                               away, ≥f_near·N within 2 units.
                               Smear-aware: reads shorter than the shortest accepted peak by
                               >1.5 units are "short products", not alleles (amplicon: 8–52%).
                               Rejected peaks and unassigned fractions are recorded.
                               Ladder prior (optional): when spanning < n_span_min, use
                               ladder idxstats / anchored-partial length lower bounds as a
                               prior and mark the length "ladder-assisted".
 S3 draft consensus            per peak: abPOA on ≤n_poa random spanning reads
                               (random, not quality-ranked).
 S4 phase split                within each peak (and across peaks ≤1.5 units apart): align
                               members to the draft; candidate sites (minor AF ≥ af_min,
                               strand-consistent); split only on ≥2 linked sites
                               (φ ≥ phi_min). A single unlinked site →
                               "unconfirmed_candidate" allele with its own consensus,
                               flagged, never a confident call.
 S5 hybrid reference           ref_a = ladder left flank + draft_a + ladder right flank.
 S6 all-read assignment        every non-off-target read (spanning, anchored, internal
                               ≥ min_frag_bp) → edlib HW edit distance to each ref;
                               assign if margin ≥ m_min (default 3), else "undecided".
                               Record per-read assignments.
 S7 polish                     1–2 rounds: align assigned reads to own draft; majority pileup
                               polish; homopolymer median vote per run (strand-balanced).
 S8 residual QC                re-align assigned reads to the final consensus; residual
                               sites with minor AF ≥ qc_af (strand-consistent) →
                               "residual_heterogeneity" (possible unresolved mixture,
                               chimera or mosaicism). Optional Clair3-on-own-consensus
                               (--hybrid-qc clair3) for parity with VNTRPipeline's QC.
 S9 classify                   existing classify_sequence (unchanged) on the ACGT consensus.
 S10 event evidence            per classified event: read-level support on the assigned reads
                               (alt-supporting / ref-supporting / other, per strand).
                               For homopolymer events (C/G/A/T runs), use the stutter-aware
                               LLR (per-sample, per-strand background from same-run
                               non-event sites). Status: supported / insufficient_depth /
                               discordant / not_supported.
 S11 outputs                   alleles.json, consensus_*.fa, *_context.json (engine=hybrid,
                               anchors, reads), repeats.*, summary.json (+ hybrid block),
                               optional per-read assignments TSV.
```

The ladder is used for four things: flanks (S5), repeat numbering and
nomenclature (unchanged classifier), a length prior when spanning reads are
insufficient (S2), and optionally seeding S3 when a peak has fewer than
`n_poa_min` spanning reads. In that last case the ladder contig of the
estimated length is polished iteratively with assigned reads, and the output
is marked `ladder_seeded`.

## 4. Evidence contract and clinical integration

Per allele (added fields):

- `engine: "hybrid"`
- `spanning_reads`, `assigned_reads`, `undecided_reads`
- `depth_status` ∈ {adequate (≥30 spanning or ≥40 assigned), low (10–29), insufficient (<10)}
- `length_basis` ∈ {spanning_peak, ladder_assisted}
- `consensus_basis` ∈ {poa, ladder_seeded}
- `residual_sites`
- `split_basis` ∈ {length, linked_sites, unconfirmed_single_site, none}

Per sample (added fields):

- `rejected_peaks`
- `unassigned_spanning_fraction`
- `short_product_fraction`
- `read_categories`

Per event (added fields): `read_support` = {alt, ref, other, alt_frac,
strand_alt_frac, llr (homopolymer only), status}.

Existing fields, set explicitly:

- `vcf_support: false` and `vcf_support_status: "not_applicable_read_consensus"`.
- `phase_status`: `phased` when split by length or linked sites;
  `no_informative_heterozygosity` when homozygous; `unresolved_single_site`
  otherwise (new).
- `reconstruction_status`: `complete_segmentation` / `ambiguous_reconstruction`
  from classification. New `read_consensus_low_depth` when depth_status ≠ adequate.
- `independent_haplotype_evidence`: true only for length or linked-site splits.
- `sequence_source`: `hybrid:<allele>`.
- `contig_name`: `hybrid_<allele>`, with a FASTA of the hybrid references
  written for IGV.

Required changes elsewhere (the engine plan includes them):

- `report.compute_clinical_decision` requires **explicit** support. PATHOGENIC
  needs a frameshift, a localized event, and either
  `vcf_support_status == "exact_sequence_concordance"` or
  `read_support.status == "supported"`. A missing support field is no longer a
  pass; this also fixes a latent ladder-engine risk.
  NO_PATHOGENIC additionally requires every allele `depth_status == adequate`,
  no `unconfirmed_single_site` allele, `rejected_peaks` empty or all below the
  noise rule, and `unassigned_spanning_fraction ≤ 0.2`. Otherwise the result is
  INCONCLUSIVE with reasons.
- `evaluation/artifacts._statuses_allow_completion` accepts the new statuses
  deliberately. `evaluation.models.Event.supported` accepts
  `read_support.status == "supported"` as support.

## 5. Settings (`HybridSettings`, config file section `hybrid`)

| Setting | Default (prototype-derived; tuned on dev/val only) |
| --- | --- |
| anchor_max_edits | 12 |
| min_span_units / max_span_units | 15 / 160 |
| peak_window_base_bp | 30 |
| peak_window_per_unit_bp | 0.6 |
| min_peak_reads | 8 |
| far_peak_min_frac / near_peak_min_frac | 0.03 / 0.20 |
| n_poa | 40 |
| n_poa_min | 10 |
| polish_rounds | 2 |
| hp_vote | true |
| het_af_min | 0.2 |
| het_min_group | 0.15 |
| link_phi_min | 0.5 |
| min_linked_sites | 2 |
| min_fragment_bp | 1000 |
| assign_margin | 3 |
| qc_residual_af | 0.25 |
| depth_adequate_spanning | 30 |
| depth_low_spanning | 10 |
| hp_llr_min | 10 |
| hp_min_reads | 20 |
| hp_min_alt_frac | 0.30 |
| seed | 1 |

CLI: `--engine {ladder,hybrid}` and `--assay {amplicon,genomic}` (#58). The
assay sets defaults: genomic lowers depth thresholds' meaning and expects no
smear or PCR bias.

## 6. Dependencies

New optional extra `hybrid`: `edlib` (MIT, wheels), `pyabpoa` (MIT; sdist
build, bioconda binary) with `pyspoa` (MIT, wheels) as a fallback backend
behind a small `PoaBackend` protocol. mypy overrides are needed for the untyped
modules. The Dockerfile, conda env and Makefile `UV_TEST` gain the extra.
A missing extra with `--engine hybrid` is a clear error (jinja2 pattern).

## 7. Acceptance

1. The unit suite covers every stage with synthetic reads (strands, smear,
   Δ1, equal length, homopolymer events, low depth, mutated anchor).
   `make ci-check` passes, coverage stays ≥80%, all files are <650 lines.
2. Frozen regressions: the 40-case simpanel and the 40-case held-out panel are
   at or above the prototype (≥76/80 held-out exact), with the 3 known
   prototype failures fixed or explicitly flagged INCONCLUSIVE.
3. PRJEB92208: 4/4 dupC positives PATHOGENIC with read support, HG002 2/2
   exact, no pathogenic calls in HG001–HG004.
4. In-house genomic: no phantom fragment alleles, and depth statuses reported.
   LB-level results stay local.
5. The benchmark decision rule (benchmark spec §6) is met on the sealed test
   split before `hybrid` becomes the default.
6. Docs updated: CLI reference, configuration guide, limitations, and changelog.
