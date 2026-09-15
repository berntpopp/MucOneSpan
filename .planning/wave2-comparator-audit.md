# VNTRPipeline comparator audit

## Scope and version identity

Read-only scientific/source review; no comparator execution or production changes by this reviewer. Reviewed full publication-era pipeline/orchestration, length selection, Canu, boundary extraction, modal sequence, TRviz, LoF and QC scripts, plus all substantive v1–v2 diffs, v2 report generation, config and Dockerfile. External checkout: `sources/VNTRPipeline` under the Wave2 data root. Paper/supplement provenance remains in `wave2-source-research.md`.

The paper provides a repository URL but no source commit or container digest. It was received July2025, accepted November2025 and published December2025, with2026 volume numbering. The last source commit before publication is tagged **v1.0**, `f7f594e74ce9cde273deb269f06bb528f0586997` (2025-10-10). This is the closest identifiable publication-era source, not proof of the exact analysis build. **v2.0**, `a8d3607d66324a936c64ed5b7e403a32bcb4c0a6` (2026-03-02), changes length processing and reporting after publication.

Lead independently inspected the cached distributed image `ghcr.io/dhmeduni/vntr_pipeline@sha256:888583f8ef0b69af9c5b386c51c18a2f6a73b0d46dac06424ebd56b01845d5a8`, built2026-03-02: shell/config hashes match v2 source; Canu identifies as master+28 r10516, commit `fb43b3c196a8cbd66a6cd0d3178e7e879eeb918d`. Registry tags v1.0/v2.0 were not available to the lead. This reviewer has not independently executed the image. The selected practical comparator is unmodified **v1 source over this pinned runtime**, explicitly a publication-era-source comparison with newer dependencies, not an exact reproduction of published software/environment. The lead owns execution and freezes actual versions/hashes.

Authoritative source links:

- [v1 source tree](https://github.com/DHmeduni/VNTRPipeline/tree/f7f594e74ce9cde273deb269f06bb528f0586997/VNTR_Pipeline)
- [v1 main entrypoint](https://github.com/DHmeduni/VNTRPipeline/blob/f7f594e74ce9cde273deb269f06bb528f0586997/VNTR_Pipeline/scripts/bin/VNTR_pipeline.sh)
- [v1 mapping/phasing/QC](https://github.com/DHmeduni/VNTRPipeline/blob/f7f594e74ce9cde273deb269f06bb528f0586997/VNTR_Pipeline/scripts/bin/variant_caller.sh)
- [v1 length algorithm](https://github.com/DHmeduni/VNTRPipeline/blob/f7f594e74ce9cde273deb269f06bb528f0586997/VNTR_Pipeline/scripts/bin/find_maxima.py)
- [v1 correction/extraction orchestration](https://github.com/DHmeduni/VNTRPipeline/blob/f7f594e74ce9cde273deb269f06bb528f0586997/VNTR_Pipeline/scripts/bin/whole_pipeline.sh)
- [v1 default config](https://github.com/DHmeduni/VNTRPipeline/blob/f7f594e74ce9cde273deb269f06bb528f0586997/VNTR_Pipeline/scripts/lib/vntr_variables.cfg)
- [paper](https://doi.org/10.1038/s41598-025-30441-3), Methods pp3–5, Further evaluation p5, QC pp6–7, Supplement1/4/5.

## Exact v1 input/default workflow

### Common mapping and reference selection

`VNTR_pipeline -i FILE -o ABSOLUTE_OUTPUT -r t2t -p pcr -v MUC1 DELETE_TMP=N` is the study-style MUC1 PCR entrypoint; replace only `-p pcr` with `-p wgs` for the two deposited WGS inputs. Input can be a single FASTQ, FASTQ.gz or BAM. FASTQ is copied into fresh output temporary work; BAM is converted with samtools fastq. Do not pool biological samples or re-demultiplex ENA per-run files.

`-r t2t` resolves `/references/chm13v2.0.fa`; `hg38` resolves `/references/hg38_p14.fa`. Minimap2 indexes then aligns `-Yax lr:hq -t32 -y --MD`; samtools converts, sorts and region-filters. PCR uses MUC1 coordinates `chr1:154328103-154330802` on CHM13 or `chr1:155188487-155192239` on GRCh38. WGS expands configured start/end by25000bp on either side before samtools regional selection. This retains whole query sequences whose alignments overlap the window; it is not physical cropping of bases to the window. Mapping/chosen genome differs from MucOneSpan's synthetic ladder.

Both source versions explicitly use `/opt/rerio/clair3_models/r1041_e82_400bps_sup_v500` for ONT. Lead reports checkpoint hashes exactly match the model used in the frozen clinical baseline. Therefore model identity is controlled in the local comparison, even though the publication does not disclose its exact basecalling model/sampling rate.

### PCR physical-length partition

Config: assembly_size5k, minimum product2000bp, repeat_size60, MUC1 boundary motifs1 and9. `extract_lengths.sh` rounds size/separation down to tens: threshold2000, minimum peak separation60, filters peak±30bp. `find_maxima.py` reads `query_length` from every BAM alignment record; there is no primary-only filter in that extraction. It constructs10bp bins from lengths strictly greater than threshold. Interior local maxima must exceed their left neighbor, be at least their right neighbor, and reach `max(20,0.025*maximum_histogram_height)`.

The highest-frequency eligible peak is selected first. A longer second peak may be weaker down to the adaptive threshold; a shorter second peak must have frequency at least100% of the first. Highest-frequency qualifying second peak wins. Returned entries are `(frequency,peak)` tuples sorted lexicographically: downstream correctly unpacks them, but haplotype1/2 order is frequency-dependent, not guaranteed shortest/longest despite comments/help. Histogram edge bins are not eligible local maxima. The ±30bp read filter uses SAM sequence length; it does not select expected mutations. Fragments shorter than2000bp cannot define a length peak, unlike the all-input MucOneSpan ladder candidate process.

No raw-read physical flank-to-flank measurement drives this stage: it uses whole read length as an amplicon proxy. Boundary extraction happens after Canu. This distinction matters for primer/barcode tails, concatemer/chimera reads and WGS.

### Clair3 and WhatsHap branches

Paper Methods describes length separation and additionally variant-based phasing. Default v1 source attempts length separation and uses Clair3/WhatsHap as fallback; it does not automatically produce both methods for every successful length-separated sample. When invoked for PCR, Clair3 uses32 threads and `--var_pct_full=1 --ref_pct_full=1 --enable_variant_calling_at_sequence_head_and_tail --var_pct_phasing=1`. These are different from MucOneSpan's frozen per-candidate calling defaults. WGS uses the model/platform with ordinary Clair3 options, without those PCR switches.

WhatsHap `phase --indels --ignore-read-groups`, VCF sort/index, `haplotag --ignore-read-groups --output-haplotag-list`, then `split --discard-unknown-reads` produce haplotype BAMs. This creates read partitions for correction/assembly; it does not create the final MUC1 sequence by VCF replay. Unassigned reads are discarded in this branch and must remain visible in read accounting.

Source orchestration has caveats requiring actual-log validation: fallback tests the first captured output line's word count rather than a typed maxima count; a one-peak output may not trigger phasing as the paper narrative suggests. `WHATSHAP_FORCE=Y` chooses options without the length stage, but the PCR phasing loops inspect a `result_array` populated by that stage; the force path can iterate an empty array. These are static control-flow risks, not measured failures of the selected default execution. Do not silently repair/force branches or claim both methods ran without output evidence.

### Canu and physical boundary extraction

`canu_pcr.sh` converts each haplotype BAM to FASTQ and calls Canu with `genomeSize=5k -nanopore-raw ... -readSamplingCoverage=100 -contigFilter="2 0 1.0 0.5 0"`. The command can execute full Canu stages, but downstream PCR uses corrected/trimmed **`trimmedReads.fasta.gz`**, not Canu's final contig or Clair3 consensus. The paper describes trimmed FASTQ; checked source consumes trimmed FASTA. Record actual generated names.

WGS sets genomeSize50k and Canu `-stopOnLowCoverage=0 -minInputCoverage=2`; downstream uses **`contigs.fasta`**. Neither wrapper explicitly limits Canu maxThreads/maxMemory. The runtime autodetects these; host/cgroup behavior requires measured logs. The paper cites Canu2.2; distributed image's actual master build differs.

Both modes search both orientations for full left/right MUC160bp anchors within edit distance1 and preserve both complete anchors in extracted sequence. Reverse-oriented sequences are reverse-complemented to coding direction. v1 hardcodes1 in both modes; v2 exposes WGS `MISMATCHES` with default1.

`fuzzy_search.pl` strips whitespace from its entire input stream, including FASTA line separators, and scans each approximate start to the first later approximate end. It does not parse sequence records or explicitly remove headers. Consequently this implementation is not formally guaranteed to constrain each hit to a single molecule; malformed/cross-record hits may arise if boundaries are missing. Check actual valid alphabet, anchors and extraction records, rather than assuming “molecule-by-molecule” correctness solely from the paper description.

### Modal sequence and QC evidence

`count_seqs_of_canu_trimmed.R` reverse-complements reverse hits, combines forward/reverse strings, tabulates exact unique corrected strings and sorts by frequency. `*_result.fasta` contains unique sequences with counts; `*_best_hit.fasta` is the first modal string. In WGS, these counts refer to assembled extracted sequences, not raw molecule support. Ties inherit input/factor ordering rather than a demonstrated unique-majority decision. The paper requires reviewing sequence distributions manually; code does not impose a published statistical confidence threshold.

`seqs_distribution.xlsx` reports the unique-sequence counts, sum(other counts)/largest count, and second/largest ratio. These are **not** the modal fraction among all assigned raw reads. Canu correction and subsampling mean corrected-string agreement is not independent molecular replication, but it provides additional sequence-diversity evidence absent from MucOneSpan's current ladder-VCF replay summaries.

For PCR, raw assigned reads are remapped to combined modal consensuses. bcftools mpileup/call and Sniffles produce QC variant files; if either contains records, a `###ERROR...vcf_non_zero...` marker is written. An additional marker warns of same-length sequences produced from length partitioning. These are review artifacts, not necessarily nonzero pipeline exit states. Top-level scripts lack strict shell failure propagation, so validate scientific outputs and per-stage logs even when the command exits0.

## Motifs, LoF and event identity

`trviz_script.py` feeds actual extracted sequences and the25-motif MUC1 list to `TandemRepeatVizWorker.generate_trplot` with package defaults. Its installed TRviz/decomposition/MAFFT versions need runtime capture; source does not pin them. The DNA motif sequences in v1/v2 are identical after uppercasing (independently compared); v2 changes color quoting and adds protein tables/reporting, not a new underlying set of MUC1 DNA sequences.

`analyse_seqs_from_trviz_for_pipeline.R` saves motif dictionary matches and calculates reconstructed sequence frame/stop metrics, but its final whole-sequence LoF spreadsheet write is commented out. The displayed/file LoF classification comes from `trviz_plot_modified.R`: known motifs get their assigned labels; unknown motifs are marked **L** only if motif length is not divisible by3 or a frame1 stop codon is present, otherwise **N**. `NON_CODING=Y` disables LoF classification. `new_and_lof_seqs.xlsx` contains haplotype, N/L class, motif position and actual motif sequence.

The paper explicitly says these rows require **manual investigation** of inserted/deleted base identity and position. There is no standardized exact named `dupC` event output analogous to MucOneSpan's endpoint schema. An L flag or a red tile cannot be scored directly as dupC. For matched evaluation, derive C-tract insertion identity from actual full extracted sequence/motif, document equivalent homopolymer representations, and keep generic LoF alarms separate. Motif positions in v1 spreadsheet originate from the aligned motif representation and may include alignment gaps; v2 changes spreadsheet extraction to the unaligned decomposition input. Neither can be equated blindly with biological core-repeat indices.

The maternal HG00218bp expansion is divisible by3; if correctly decomposed and stop-free it should be **N**, not L. This is a source-derived expectation, not an observed local comparator result. v2 PDF added `UNREMARKABLE` versus `FS/Nonsense detected` using any L row; missing/failed LoF files can leave the default UNREMARKABLE text, so a PDF alone never certifies a negative result.

## Disclosed manual settings and unresolved study choices

Supplement1 exposes LENGTH_1/LENGTH_2, MIN_FREQUENCY and WHATSHAP_FORCE; it does not identify per-participant override values or a complete command manifest. Paper reports evaluating both separation methods and specifically says MUC1 MP4 could not be correctly phased; MP5 showed two similarly frequent corrected sequences21/19 in phasing QC. Published final figures must not be treated as unconditionally automated default-run outputs. No evidence was found that diagnostic labels drove preprocessing, but no per-sample override record establishes that no manual choices were made.

Paper Further evaluation explicitly combines already extracted benchmark/reference FASTAs with best_hit FASTAs and uses `ALL_FIGURES=Y` to make comparisons. Those figures are not independent freshly called comparator outputs. Do not add Q100 truth into the directory used for calling; compare it only afterward.

## Why the comparator may avoid observed MucOneSpan failures (hypotheses)

1. Excluding sub2000bp fragments from length peaks and using a weak-longer-peak rule may retain long alleles that MucOneSpan loses to tiny candidate clusters in HG002/HG004/MP2/MP4. This requires observed comparator spans and counts; it is not proof those defaults succeed on every library.
2. Correcting physical reads within length/phase partitions and taking their modal boundary-contained string avoids constructing the final sequence solely from a selected reference plus VCF genotypes. It can retain an event even when a heterozygous site's reference GT is chosen by MucOneSpan replay. Canu can also introduce consensus bias or remove rare haplotypes; measure both.
3. Generic motif decomposition of actual sequence followed by frame testing may preserve MP1's C-tract insertion rather than shifting it into an A-insertion interpretation, and may leave the HG002 in-frame expansion non-LoF. These are distinct sequence/segmentation/display hypotheses, not evidence that a model upgrade or one threshold fixes all failures.
4. Corrected-sequence distributions and post-consensus remapping provide reviewable uncertainty information. They are not equivalent to a molecular concordance probability or independent truth.

## Concrete local comparison contract for the lead

Primary source overlay: v1.0 commit above, runtime image digest above. Preserve image script executable permissions when overlaying v1 source: Git v1 script modes are100644, while v2 marks them executable. Document executable-mode handling as packaging, not scientific modification. Freeze source and effective config hashes, actual runtime tools (especially Canu/TRviz/Clair3), T2T reference and model checkpoint hashes.

The initial contract used unchanged prepared ENA FASTQ; the observed header incompatibility below requires a documented comment-only adaptation for subsequent attempts. Run one adapted FASTQ per fresh output with `-r t2t -v MUC1 -p pcr` for all nine amplicons and `-p wgs` for the two WGS inputs, `DELETE_TMP=N`. Leave LENGTH_1/LENGTH_2/MIN_FREQUENCY/WHATSHAP_FORCE/VNTR_ALL/NON_CODING/benchmarking unset. No sample-specific length overrides, forced success, threshold tuning, or selection using reported mutation. In v1, config parsing occurs after key=value parsing and can overwrite same-name CLI scientific variables; this makes explicit override experiments especially unsuitable for an unmodified baseline.

Lead's intended resource wrapper is2CPUs/16GiB; source still requests32 threads and Canu autodetects resources. Record this distinction and measured process-tree resource scope. Do not describe it as native two-thread comparator settings or claim a matched speedup. Retain timeout/OOM/missing output states and every eligible input; no container exit0 or PDF suffices to establish usable sequence.

For each attempt save raw invocation/status, logs, mapped counts, length histogram/selected peak windows, partition counts, Canu correction/trim counts, both full best_hit/result sequences and sequence-distribution tables, motif map and LoF workbook, QC error markers. Check fresh files, valid DNA, both anchors, sequence/haplotype counts and uniqueness. Compare final raw sequences to independent Q100 with the same orientation/boundary and allele-permutation rules; score reported-known-control exact C-tract recovery separately from generic LoF. Classifier/LoF display cannot serve as its own truth.

If distributed-v2 source is also run, report it as a separate exploratory comparator: minimum product1800bp, separation50/bands±25, minimum frequency10, longer ratio1%, shorter ratio20%, corrected config precedence, variable WGS anchor tolerance, protein/PDF output and altered LoF-position export. Do not pool these results with v1-source results. No exact published performance/speed comparison is available without the original paper image/dependency/command manifest.


## Observed first-attempt compatibility failures and synthetic confirmation

The initial MP1 attempt in external `comparator-v1/cohort/ERR15277566/container-full.log` produced no usable sequence output despite wrapper exit0. It reached minimap2 alignment of5941 reads, then samtools reported `aux_parse` unrecognized types. Subsequent missing haplotype-directory, R `setwd`, and TRviz `FileNotFoundError` messages are downstream consequences of the failed SAM ingestion. This is an execution/input-format failure, not a negative biological finding. Container-main OOM=false does not exclude child-process OOM; this log establishes a concrete parse failure without establishing a complete memory diagnosis.

The v1 `variant_caller.sh` alignment uses `minimap2 -Yax lr:hq -t 32 -y --MD`. The [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html) documents `-y` as copying FASTA/Q comments to output. ENA headers have an accession first token followed by a UUID-like comment, which is not a SAM auxiliary field. The [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf) defines optional fields as TAG:TYPE:VALUE. Copying an arbitrary UUID verbatim therefore yields invalid SAM optional-field syntax.

An independent four-arm synthetic probe ran inside the exact pinned image, with minimap2 2.30-r1287 and samtools/htslib1.22.1. It used a deterministic4000bp random reference and one matching3000bp read; no participant sequences or identifiers were used. Full script, inputs, commands, versions, SAM/BAM files, logs and hashes are external in `sources/header-comment-probe/`; `result.json` records every result.

| Synthetic input / option | minimap2 exit | samtools exit | Result |
| --- | ---: | ---: | --- |
| UUID-like comment, `-y` | 0 | 1 | `aux_parse` unrecognized type `b`; UUID appended verbatim |
| Comment removed, `-y` | 0 | 0 | Valid SAM/BAM |
| Original comment, no `-y` | 0 | 0 | Valid SAM/BAM |
| Comment represented as valid `XX:Z:` tag, `-y` | 0 | 0 | Valid SAM/BAM |

Every synthetic SAM retained the original first read identifier, bases and quality; the stripped FASTQ retained all non-header lines byte-for-byte. Removing only header comments is therefore a supported input compatibility adaptation while leaving the v1 source and alignment options unchanged. It is **not byte-identical input** or merely decompression. For real adapted inputs, retain source/adapted file hashes, require equal record count/order, first identifiers, bases and qualities, and preserve the source files separately. Apply the same deterministic rule to every eligible input; no read selection, trimming, quality conversion, orientation change or expected-event filtering belongs in this adaptation. Validate actual FASTQ structure rather than stripping whitespace indiscriminately from sequence/quality lines. The synthetic result validates the mechanism; the lead must separately verify the real-input transform.

### Independent R plotting dependency failure

The full log tail also contains `Error in contrib.url(repos, type): trying to use CRAN without setting a mirror`. Unlike the missing-directory errors, this has an independent source explanation. `trviz_plot_modified.R` lines4–5 attempts `install.packages("BiocManager")` whenever that namespace is absent, before inspecting input data. A read-only dependency probe in the pinned image's `/opt/conda/envs/python-env/bin/Rscript` confirmed BiocManager absent, while plotrix, ggplot2, openxlsx, Biostrings and ShortRead are present; `getOption("repos")` returns `@CRAN@`. Earlier scripts set a mirror and issue installation warnings, then load available biological packages, but this plotting script does not set one and halts immediately.

Thus successful SAM ingestion alone cannot establish that the current runtime will produce LoF plotting/workbook outputs. A separately recorded, pinned BiocManager runtime addition can address this startup dependency without altering scientific source; merely setting a mirror under disabled networking does not install the missing package. Preserve the original failed attempt and distinguish any subsequent runtime dependency overlay in provenance. No package installation, source edit or cohort rerun was performed by this audit agent.
