# Wave 2 source audit (2026-09-15)

## Sources and preserved evidence

External root: `${DATA_ROOT}/sources`.
No complete participant sequence, raw read, or raw pipeline output is committed; compact provenance and aggregate benchmark records are curated separately. `SHA256.json` records downloaded source-file hashes. Primary paper PDF is `paper.pdf`, extracted `paper.txt`, publisher HTML `paper.html`.

- Paper: https://doi.org/10.1038/s41598-025-30441-3 ; PDF https://www.nature.com/articles/s41598-025-30441-3.pdf . PDF pages below refer to printed page numbers (1–10).
- Supplements share base URL `https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41598-025-30441-3/MediaObjects/` and filenames `41598_2025_30441_MOESM1_ESM.docx`, `41598_2025_30441_MOESM2_ESM.pptx`, `41598_2025_30441_MOESM3_ESM.pptx`, `41598_2025_30441_MOESM4_ESM.xlsx`, `41598_2025_30441_MOESM5_ESM.xlsx`, `41598_2025_30441_MOESM6_ESM.docx`. All downloaded; text extracted from underlying XML into matching `.txt` files. Original graphical supplements retained for inspection.
- Wenzel 2018: https://doi.org/10.1038/s41598-018-22428-0 ; https://www.nature.com/articles/s41598-018-22428-0.pdf (`wenzel2018.pdf`, `wenzel2018.txt`).
- VNTRPipeline cloned outside Git. v2.0/current HEAD `a8d3607d66324a936c64ed5b7e403a32bcb4c0a6` (2026-03-02); v1.0 `f7f594e74ce9cde273deb269f06bb528f0586997`. https://github.com/DHmeduni/VNTRPipeline/blob/a8d3607d66324a936c64ed5b7e403a32bcb4c0a6/README.md .

## Truth resolution and participant identity

Paper p2 Methods Sample collection establishes Coriell identities HG001=NA12878; HG002=NA24385; HG003=NA24149; HG004=NA24143. It calls MP1/MP2 in-house positive controls and MP3/MP4 Cologne positive controls; MP5 was an undiagnosed suspected case after WES. Paper p6 explicitly identifies MP1/MP2 as siblings; p6 identifies HG003/HG004 as father/mother of HG002. PCR and WGS are repeated library evidence, not additional biological participants.

Paper pp6–7 calls MP1–MP4 known dupC controls, specifically a single C insertion into a seven-C homopolymer. This supports *reported pre-existing event-positive status*, with orthogonal assay unspecified. Neither methods nor downloaded supplements supplies an assay result or sample identifier connecting MP1–MP4 to orthogonally characterized patients in cited work. Do not claim Sanger confirmation. Wenzel 2018 Methods p3 and Results p6 establish SNaPshot/probe-extension confirmation for its own F-number families (SMRT first then SNaPshot in F7), but no verified mapping to MP1–MP4 is published in this study. The Supplement 1 notes that a motif equals Wenzel motif L; motif equivalence is not participant identity. Consequently keep `orthogonal_method: null` and distinguish reported known-control recovery from independently documented-assay sensitivity.

MP5 is comparator-only exploratory evidence: reported single-C insertion at nucleotide 23 in motif35/83. Its existence, exact position and causality are not independently confirmed. HG001/3/4 reference labels establish neither absence of every target event nor full sequence truth. HG002 full sequence independently sourced below; any negative endpoint must be explicitly derived/validated against those haplotypes, not inferred from its alias.

Published positions (paper p6 and Fig5) are comparator observations: MP1/MP2 motif17 of78, nucleotide59, longer haplotype; MP3 motif7 of44, nucleotide59, shorter; MP4 motif49 of80, nucleotide59, longer; MP5 motif35 of83, nucleotide23. They must not become independent position truth. Homopolymer insertion representations can shift within the C run; publication coordinates are coding-orientation one-based motif nucleotide positions, not genomic HGVS coordinates. Boundary-inclusive motif numbers cannot directly equal a caller's core-repeat index without explicit mapping.

## Important source inconsistencies

1. Paper p5 Results says HG002 and **MP4** were sequenced by WGS, whereas ENA records the patient WGS alias **MP1_WGS_MUC1**. Preserve ENA alias, mark biological mapping unresolved, and avoid merging that WGS into confirmed sample metrics by unsupported identity inference.
2. Supplement5 spreadsheet sheet1 rows39–40 label MP5 WGS with lengths3900/4638; rows45–46 label HG002 WGS with2640/4681. The first pair agrees with independent HG002 and the second pair with published MP1 PCR. This strongly suggests table labeling errors but does not authorize relabeling participants from length alone.
3. Supplement1 MUC1 parameter table prints ACAN boundary strings, although genomic coordinate annotations correspond to MUC1. Its MUC1 motif visualization (motifs1 and9), and pinned comparator config, supply the appropriate MUC1 anchors. Use actual sequence evidence and pinned config, not the erroneous strings.

## Independent HG002 Q100 v1.1 full sequence acquired

Authoritative project: https://github.com/marbl/HG002 ; v1.1 released July2024, linked GenBank paternal GCA_018852605.3 and maternal GCA_018852615.3. Paper used Q100 through UCSC but does not specify the version; therefore describe this benchmark as independently acquired v1.1 truth, not necessarily the exact paper download. Q100 derives from independent HPRC/GIAB sequencing, not these PRJEB92208 PCR reads.

UCSC hub: https://research.nhgri.nih.gov/CustomTracks/T2T_hubs/T2Tgenomes/hub.txt . Hub genomes file links versioned twoBit sequence files:

- https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/browser/HG002/v1.1/2bit/hg002v1.1.pat.2bit
- https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/browser/HG002/v1.1/2bit/hg002v1.1.mat.2bit

Both boundaries occur exactly once in the queried 10Mb locus neighborhoods. Extract reference plus-strand interval, reverse complement once to MUC1 reading-frame orientation, preserve both 60nt anchors. Exact anchors from pinned comparator config (also Supplement1 motif1/9):

- left `AAGGAGACTTCGGCTACCCAGAGAAGTTCAGTGCCCAGCTCTACTGAGAAGAATGCTGTG`
- right `GGCTCCACCGCCCCTCCAGTCCACAATGTCACCTCGGCCTCAGGCTCTGCATCAGGCTCA`

| Haplotype | Chromosome | Zero-based half-open interval | Orientation | Boundary-inclusive length | SHA256 of uppercase sequence without newline |
|---|---|---|---|---:|---|
| paternal | chr1_PATERNAL | [158194262,158198162) | minus |3900|51f51aae0ffecc6449b753ce20b4e85c2268940570b81788e37c7f065cc4c3d2|
| maternal | chr1_MATERNAL | [150024488,150029126) | minus |4638|675bc2632ae3836189af8135d17f473eeaef2a6fb7bd2ba18f70a33e5d05fba6|

Exact API requests (semicolon-separated UCSC API parameters):

- https://api.genome.ucsc.edu/getData/sequence?hubUrl=https://research.nhgri.nih.gov/CustomTracks/T2T_hubs/T2Tgenomes/hub.txt;genome=HG002v1.1.PAT;chrom=chr1_PATERNAL;start=158194262;end=158198162
- https://api.genome.ucsc.edu/getData/sequence?hubUrl=https://research.nhgri.nih.gov/CustomTracks/T2T_hubs/T2Tgenomes/hub.txt;genome=HG002v1.1.MAT;chrom=chr1_MATERNAL;start=150024488;end=150029126

External FASTA `q100-PAT-muc1-boundary-inclusive.fasta`, `q100-MAT-muc1-boundary-inclusive.fasta`; metadata `q100-muc1-provenance.json` and per-haplotype `q100-{PAT,MAT}-provenance.json`; direct API responses `q100-{PAT,MAT}-muc1-api.json`. Independently repeat exact interval requests and confirm same sequence hashes (performed). Exact agreement permits allele permutation, never removal of biological indels or novel motifs. Maternal 4638bp preserves an 18bp in-frame motif expansion. Repeat counts 65/77 include boundary motifs in the paper; division of4638 by60 is not a valid repeat-count metric. If caller output spans different boundaries, independently locate and record identical anchors before comparing; inability to recover both full spans is uncallable sequence reconstruction, not absence of truth.

## Basecalling, input selection, resource claims

Paper p3 Data analysis: Dorado0.7.1 SUP used for both PCR and WGS. Exact Dorado model identifier/checksum absent, so not reconstructable from version/accuracy tier alone. p3–5 workflow: minimap2 2.28 alignment, samtools1.21 regional filtering, Clair3 1.0.11, WhatsHap2.3, Canu2.2, TRviz. The paper establishes regional selection within its workflow, but does not explicitly state that ENA's deposited WGS FASTQs are those filtered products. Tiny deposited files plus later read-level mapping provide evidence for a locus-selected input; until inspected, do not claim whole-genome throughput or exact original filtering command.

No quantitative analysis runtime, peak RAM, computing hardware, or thread-count benchmark was found in the main paper, all six supplements, or README at v1.0/v2.0. Lab incubation times and Canu genome-size/coverage parameters are not runtime/memory observations. Supplement5 supplies read coverage, not computational resource measurements. Therefore published matched speedup and memory comparison are not estimable; report the locally measured frozen-caller scope only. Version2 README now advertises WGS and PCR/BAM/FASTQ, and test-data packaging; this is software capability documentation, not evidence of a PacBio study arm.

## Lead integration requirements

Use categories for pre-existing reported positives, independently confirmed assay claims, comparator-only positions, independent complete sequences, and unknown negatives. Preserve all four reported controls for clearly labeled known-control recovery, while explicitly stating unknown orthogonal methods. Treat ambiguous patient WGS identity separately. Add actual HG002 sequence endpoint using the external versioned truth above rather than claiming sequence truth inaccessible. Do not commit full participant sequences. Curated ledger can store assembly/coordinates/anchors/hashes and fetch provenance; runtime extraction can remain outside Git. Source errors should remain visible in public limitations.
