# Hybrid engine IGV report with a detected-mutation track

## Request

Show the location of every detected variant in the IGV view of the HTML report,
not only the alleles. Refactor as needed.

## Findings

- 0.17.0 rejected `--report-igv` for the hybrid engine (no alignment exists).
- `hybrid_references.fa` holds the *draft* alleles (pre-polish) in ladder flanks;
  `repeats.json` coordinates refer to the *polished* consensus. A manual alignment
  to `hybrid_references.fa` therefore cannot carry mutation coordinates exactly.
- `characterize_differences` reports 1-based positions in the reference repeat
  unit; template matches (e.g. dupC) record no differences.

## Decisions

1. With `--report-igv embedded|sidecar`, a hybrid run writes
   `igv/igv_reference.fa`: each polished consensus wrapped in
   `hybrid.assign_flank_bp` ladder flanks (reuses `hybrid_references`), contigs
   `hybrid_<allele>` = the summary `contig_name`.
2. Reads are aligned with the existing `map_reads` (minimap2 + samtools) to that
   reference, display only. New validated setting
   `hybrid.igv_minimap2_preset` (default `lr:hq`, HiFi and ONT Q20+), so the
   hybrid engine stays platform-free. `--threads`/`--mapping-timeout` are used
   for this alignment and no longer listed as ignored when IGV is on.
3. `igv/mutations.bed` (BED8): the feature spans the classified repeat unit;
   thickStart/thickEnd mark the changed bases (window positions from the
   unit-to-window differences with accumulated indel shift; template matches are
   re-characterized against their parent unit). Shown as the
   "Detected mutations" annotation track above the alignments.
4. `igv/loci.bed`: one navigation row per allele and per mutation.
5. Runs without `--report-igv` are unchanged (no new files, same tools).
6. Ladder engine unchanged (deprecated; its consensus coordinates differ from the
   ladder contig after indels, so the projection is hybrid-only).

## Validation

- Unit: coordinate mapping (substitution, isolated deletion, insertion shift,
  dupC on a synthetic allele), BED writing, track config, settings validation,
  ignored-options semantics, end-to-end hybrid run with IGV (mocked mapping).
- Real data: five in-house ONT WGS samples (LB25/LB26, not in the repo) ran
  with `--report-igv embedded` and igv-reports 1.13.0; the dupC thick base is a
  C of the 8-C run, and each embedded session holds both tracks.
- igv-reports 1.16.x fails the existing preflight (VCF header defect); 1.13.0
  works.

## Revision after browser check (same day)

A headless-Chromium check of the real reports showed the variant was invisible:
the reads were aligned to the *mutated* consensus, so carriers matched it. The
reference is now a display reference (mutated units restored to their canonical
parents; built after classification, so the alignment moved into the report
step via an `IgvInputs` callable). The track marks reference bases from the
unit-coordinate differences directly. Mutation navigation rows gained
`hybrid.igv_context_units` (default 1) of context, because `del1G` at a unit
edge was cut off. Names avoid `>` (it was double-escaped in the table).
Screenshots after the fix: LB25-6421 ins42G, LB25-5830 dupC show insertion
columns under the mark; the LB25-5830 del1G calls show a two-base mismatch,
not a clean deletion, in the reads.

## Revision 2 after the owner's browser review

Findings on LB25-5830 (82/83 near-identical alleles): aligning all reads to both
contigs let the alleles swap reads (allele 1 showed ~70x, the caller used 7
spanning reads); reads were not sorted at the variant; the report opened on the
whole allele; the view could not be panned beyond the row.

- Each allele now shows only the reads the engine assigned to it (spanning
  members + assigned partial reads, `HybridResult.read_names`), aligned to its own
  contig (`hybrid/igv_alignment.py`), then merged.
- Mutation rows are the 1-bp site; igv-reports sorts by BASE there (igv.js puts
  reads with an insertion starting at, or a deletion over, the site first).
  Indels are left-normalized within a run (minimap2 placed the LB25-5830 dupC
  insertion before the first C).
- `--window` (`hybrid.igv_context_units`) sets the initial view; `--flanking`
  (`hybrid.igv_padding_units`, default 10 units per side) the embedded span.
- The report opens on the first mutation row (`preferred_rows`).

Verified in headless Chromium on LB25-5830: dupC view centered, 14 carriers sorted
on top, insertion column under the mark; del1G shows a mismatch + 1-bp deletion
in 7/9 allele-2 reads (the reads' real sequence, not a clean single deletion).

## Revision 3: zoom-out and MUC1 gene track

- Zooming out past the embedded reference slice showed reads as all-mismatch
  ("soft-clip"-looking rainbow reads): igv-reports embeds whole reads but only
  the reference slice of row +- flanking/2. Every row now embeds the whole
  contig (`--flanking` = 2 x longest contig); `hybrid.igv_padding_units` was
  removed; `hybrid.igv_max_reads_per_allele` (300, seeded) bounds report size.
- Flanks verified against GRCh38 (MucOneUp reference): left = revcomp of
  chr1:155192240-155202239, right = revcomp of chr1:155178487-155188486, i.e.
  transcript orientation. Options considered for a gene track: (a) map to GRCh38
  (rejected: expanded/contracted alleles become huge indels), (b) project a gene
  model onto the display contig (chosen), (c) repeat-structure only (added as a
  second track). MANE has no MUC1; RefSeq Select NM_001371720.2 is split inside
  the VNTR; Ensembl canonical ENST00000620103.4 (= NM_001204285.2 exons) was
  bundled with provenance. Its two exon-2 blocks around the VNTR are merged.
- Display flanks = transcript extent beyond the VNTR + `igv_gene_margin_bp`
  (1176 / 3163 bp at the default 500).
- Validation: on the bundled flanks with 30/80-unit and mixed alleles the
  projected transcript has GT..AG introns, starts ATG, ends TAG, has no internal
  stop and translates to MUC1 (MTPGTQSPFF...AATSANL*). Browser check of
  LB25-5830: gene, coloured units, dupC mark and sorted carriers; zoom-out to
  7.7 kb shows the whole gene without pseudo-mismatches.
- Note: nomenclature uses NM_001204286.1; the gene track is structural context.

## Revision 4: IGV Desktop output (`--igv-session`)

Owner request: BAMs and real IGV sessions in the output folder, as a tool option.
`run.igv_session` / `--igv-session` keeps `igv/` (reference + .fai, BAM + .bai,
BEDs, `igv_session.xml` with relative paths and regions of interest); without it
the report's intermediate folder is removed. Repeat-unit BED gained
`track itemRgb="On"` (kept by igv-reports, honoured by IGV Desktop). Verified by
loading a session in IGV Desktop 2.19.7 (headless xvfb batch, isolated home):
tracks load, units are coloured, carriers sort at the dupC site.
