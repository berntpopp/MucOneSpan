# Wave 2 input geometry audit

Performed 2026-09-15 against all eleven completed preparation records, before inspecting caller outputs. Owned source: `.planning/wave2-input-audit.md`. External aggregate evidence: `${DATA_ROOT}/input-audit.json`; reproducible audit script: external `sources/audit_inputs.py`. Raw reads and identifiers remain outside Git.

## Decision before outcome inspection

Use validated identity-preserving FASTQ decompression for all nine amplicons. Attempt both WGS-labelled inputs through the existing raw FASTQ ONT path as a separate exploratory arm. Do not crop, filter by mutation or expected haplotype, reverse-complement the input files, restrict to the exact-anchor subset, or discard short fragments. Mapping already handles both orientations. Exact-anchor statistics below characterize inputs, not eligibility of individual reads.

The CLI `muc-one-span run --input INPUT.fastq --platform ont --clair3-model MODEL --output-dir OUTPUT` enters `pipeline.execute_pipeline` and `mapping.map_reads`. The latter accepts FASTQ directly and maps via minimap2 `-a -x lr:hq`, followed by samtools sort/index. No requirement in that input path prohibits locus-derived WGS reads; however the pipeline has MUC1 amplicon-oriented scientific assumptions and these observations do not validate full-genome use or equivalent WGS phasing accuracy. Preserve long WGS flanks and unsupported/insufficient downstream evidence as recorded outcomes.

## Prespecified markers and method

Full file scan, all prepared reads, exact case-insensitive DNA matching in both directions. Full left/right 60bp anchors are from the source audit and pinned VNTRPipeline config; requiring both in correct order defines a conservative observed full-span subset. Counts can underestimate true spanning reads because a sequencing error or biological boundary variant breaks an exact match. Additional20bp anchor seeds, published25bp primers, and25bp markers adjacent to or250/1000bp into each bundled genomic flank test geometry without inspecting a mutation. No reads were selected or modified. Span-tail statistics are conditional on both exact full anchors; a rare read matching both orientations could count twice, though none was expected.

Published primer source: main paper p3 MUC1 section, https://www.nature.com/articles/s41598-025-30441-3.pdf . Remove unspecified barcode positions only from the search patterns, not the input reads: PS2 `GGAGAAAAGGAGACTTCGGCTACCCAG`; PS3 `GCCGTTGTGCACCAGAGTAGAAGCTGA`. Barcodes are unspecified X bases in the paper. Full MUC1 anchors and their direction are documented in `wave2-source-research.md`. The bundled dictionary contains10000bp flanks each (`src/muc_one_span/data/repeats/repeats.json`). Marker sequences and per-run input SHA256 are in the external aggregate JSON.

## Read counts and observed geometry

Both-anchor counts show coding-forward / coding-reverse orientations. Length summaries include every read, not only full spans.

| Run | Label | Reads | Length min / median / max (bp) | Both exact anchors forward / reverse | Both anchors / all reads |
|---|---|---:|---|---|---|
|ERR15277552|HG002 WGS|54|1972 / 105351 / 397962|12 / 14|26 / 54|
|ERR15277553|patient WGS (ENA MP1)|2152|54 / 930 / 51821|4 / 1|5 / 2152|
|ERR15277562|HG001|10197|135 / 2505 / 6822|2159 / 3290|5449 / 10197|
|ERR15277563|HG002|13733|133 / 1605 / 8774|1717 / 2576|4293 / 13733|
|ERR15277564|HG003|12485|189 / 2742 / 8000|2432 / 3574|6006 / 12485|
|ERR15277565|HG004|45325|119 / 422 / 7549|2314 / 3705|6019 / 45325|
|ERR15277566|MP1|5941|140 / 2733 / 5565|665 / 1337|2002 / 5941|
|ERR15277567|MP2|5332|126 / 3085 / 30888|277 / 745|1022 / 5332|
|ERR15277568|MP3|26672|131 / 1718 / 37502|1044 / 3080|4124 / 26672|
|ERR15277569|MP4|31808|123 / 394 / 11240|437 / 1219|1656 / 31808|
|ERR15277570|MP5|8307|133 / 577 / 9065|485 / 960|1445 / 8307|

Every amplicon has both orientations and at least1022 exact full-span observations. Conditional median tails outside the anchors are32–33bp on the coding-upstream side and71–72bp downstream; this strongly supports short primer-bounded products. No amplicon read exactly matched the250bp or1000bp-distance markers on either bundled genomic flank. Thus the assay does not provide the long genomic flanks present in the full ladder reference. This is an assay geometry limitation to keep visible when judging consensus flanks; it is not a reason to synthesize missing sequence or trim informative boundary bases.

Several libraries have many short fragments: HG004 median422bp, MP4 median394bp, MP5 median577bp. These fragments persist in inventory and execution, and read count must not be interpreted as independent full-haplotype evidence. Exact primer sequences are common, consistent with primer-containing deposited reads; exact sequence searches cannot independently establish barcode/adaptor trimming history. ENA headers preserve accession/read identifiers plus original UUID-like identifiers but no detected Dorado/model field. Publication-level Dorado0.7.1 SUP provenance remains the available basecalling evidence; no exact model checksum is recoverable from these headers.

## WGS deposition scope and uncertainty

HG002 WGS (ERR15277552):54 reads,6,553,558 bases;26 exact both-anchor spans, and46 reads contain the exact right-anchor20bp seed. Both orientations present. Conditional median genomic tails36,887bp upstream and73,046bp downstream demonstrate native long locus-spanning molecules.38 and43 reads respectively contain the exact left/right1000bp flank marker. This large enrichment for one locus among only54 ultralong reads strongly supports a locus-selected WGS-derived deposit. It is not a complete genome dataset.

Patient-labelled WGS (ERR15277553):2152 reads,4,122,626 bases;median930bp,5 exact full spans.41 and51 reads respectively match the1000bp left/right flank marker. Full-span conditional median tails7919bp upstream and4821bp downstream distinguish these from the primer-bounded amplicons. The much shorter overall distribution suggests a materially different input mixture from HG002 WGS. Its small total size and MUC1-associated long reads support locus-derived data, but exact selection algorithm, original whole-genome coverage and completeness of region selection cannot be reconstructed from these FASTQs alone. Do not label all2152 reads full-length or claim whole-genome runtime. Paper/ENA/supplement biological-label inconsistency remains unresolved as documented in the source audit.

No outcome-dependent preprocessing is justified. These inputs are syntactically supported by the raw FASTQ path; scientific callability and allele support must be assessed from actual frozen execution, keeping this WGS arm distinct. Neither read count nor this marker audit supplies orthogonal event confirmation.

## Reproduction and evidence integrity

Run external `sources/audit_inputs.py` with Python against prepared files and the frozen bundled dictionary; it writes only aggregate marker counts and statistics to `input-audit.json`. Each row stores prepared input SHA256. This audit read no caller results, consensus, VCFs or clinical classifications. Preparation records report zero excluded reads and equal input/output read counts for all eleven runs. The audit independently iterated the decompressed FASTQs and obtained matching read totals.

## Pre-run Clair3 model choice and compatibility validation

Before any cohort prediction, the initially discovered model was an R9.4.1 HAC model, which does not match the reported Kit14/SUP chemistry. A compatible R10.4.1 SUP model was acquired externally without altering the caller or installed environment.

Chosen model: `r1041_e82_400bps_sup_v500`. This assumes R10.4.1 E8.2 400bps 5kHz and Dorado v5.0.0 SUP. It is a prespecified provenance-based assumption, not a verified original basecalling identifier: the paper says Dorado0.7.1 SUP, but neither exact model version nor sampling frequency is available. Dorado0.7.1 documentation lists v5.0.0 as the latest SUP model for5kHz data and also supports v4.3.0; selecting v5.0.0 matches that contemporary default, without tuning against outcome. If the original study used an explicitly chosen older model, this mismatch remains a limitation.

Sources:

- https://github.com/nanoporetech/dorado/blob/v0.7.1/README.md , available-model table and automatic model selection sections; external copy `sources/dorado-v0.7.1-README.md`.
- https://github.com/nanoporetech/rerio/blob/c7e293b3a60e12734273578523780afa245a9c8c/README.rst , May2024 table mapping this exact model to v5.0.0 SUP / R10.4.1 E8.2 5kHz.
- https://github.com/nanoporetech/rerio/blob/c7e293b3a60e12734273578523780afa245a9c8c/clair3_models/r1041_e82_400bps_sup_v500_model , versioned download URL stub.
- https://github.com/nanoporetech/rerio/blob/198f90f33434437001ad0d023e3ca912cdbca47c/README.rst , explicit compatibility with Clair3 v1.x, not v2.x. Installed Clair3 is1.0.10.

Archive URL: https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v500.tar.gz . Downloaded74,586,126 bytes, SHA256 `01c05768661bdd7de611e6bae1043c43b7523a54b223e029c683bfac0db7a678`. No upstream checksum was published in the inspected stub; this is a local content hash that freezes the retrieved artifact. Safely extracted only regular files/directories within external `models/`; no code execution from archive. Model directory: `${DATA_ROOT}/models/r1041_e82_400bps_sup_v500`. Per-file hashes and selection provenance are in external `models/model-provenance.json`.

Actual compatibility check used installed Clair3 classes `Clair3_P(add_indel_length=False,predict=True)` and `Clair3_F(add_indel_length=True,predict=True)`, matching the respective pipeline stages, loaded both checkpoint prefixes, built correctly shaped zero tensors, asserted all existing model variables matched checkpoint objects, and obtained finite predictions. This tests loading and dimensional compatibility, not scientific accuracy. Details in `models/model-load-validation.json`.

The initial synthetic probe detected a compute-capability12 GPU unsupported by the installed TensorFlow binary and began potentially prolonged JIT initialization. That task-owned probe was terminated. The CPU probe used process-local `CUDA_VISIBLE_DEVICES=''`, one intra/inter-op thread, and completed successfully. An intermediate probe incorrectly enabled indel-length heads for pileup and failed strict matching; matching actual stage options resolved the test setup error. No weights, source, environment installation or study reads were changed. Cohort CPU execution should explicitly disable GPU visibility for reproducible hardware scope.
