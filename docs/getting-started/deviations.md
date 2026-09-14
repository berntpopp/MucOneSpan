# Differences from the Published Method

MucOneSpan follows the architecture described in Vrbacka et al. 2025 (Figure 1) but makes several deliberate changes to improve accuracy, reproducibility, and usability.

---

## Comparison Table

| Aspect | Published method (Vrbacka et al.) | MucOneSpan | Rationale |
|--------|--------------------------|--------------|-----------|
| **Read mapper** | bwa-mem 0.7.16a | minimap2 2.28+ | minimap2 is the standard for long-read alignment; faster and more accurate for HiFi data |
| **Reference ladder** | 120 contigs (1-120 X repeats) | 150 contigs (1-150 X repeats) | Covers the full observed VNTR length range (Vrbacka reports alleles up to 125 repeats) |
| **Variant caller** | Clair 2.0.7 + Clair3 1.0.10 | Clair3 only (latest HiFi model) | Clair is deprecated; Clair3 subsumes its functionality with better accuracy |
| **Read filtering** | Manual IGV inspection for chimeric/incomplete reads; MAPQ > 4 filter | Automated read extraction per cluster + remapping to peak contig | Reproducible; no manual step required |
| **Allele detection** | `peaks.py` (idxstats peaks, details unpublished) | Gap-based clustering + indel-valley splitting | The indel-valley algorithm resolves close allele pairs (3-9 repeats apart) that simple peak-finding merges into one |
| **Mutation detection** | Custom Python scripts (unpublished) | Pre-computed mutation template catalog (13 known mutations) with exact-match-first probing | O(1) lookup for known mutations; no edit distance needed for common cases |
| **Flanking trim** | Fixed-position trim (assumes flanking is invariant) | Anchor-based trim using flanking/pre-repeat boundary sequences | Resilient to indels in flanking region caused by Clair3 false positive calls |
| **Same-length alleles** | Treated as homozygous (Clair3 skipped) | Disambiguated using Clair3 het genotype detection | Detects compound heterozygotes where both alleles have the same repeat count |
| **Confidence scoring** | Not reported | Per-repeat and per-allele confidence with VCF cross-validation | Enables automated QC; flags low-confidence calls |
| **Source code** | Proprietary (unpublished) | MIT license, fully open source | Reproducibility and community contribution |

---

## Detailed Explanations

### Read Mapper: bwa-mem to minimap2

The workflow described by Vrbacka et al. used **bwa-mem 0.7.16a**, which was designed for short reads. While bwa-mem can align long reads, **minimap2** is the community standard for PacBio HiFi data. It is faster, produces more accurate alignments for reads >1kb, and handles the high GC content of the MUC1 VNTR more reliably.

### Reference Ladder: 120 to 150 Contigs

Vrbacka et al. report alleles with up to 125 repeat units. The original ladder (1-120) would miss these extreme alleles. Extending to 150 contigs provides headroom for rare long alleles without significant computational cost (the ladder FASTA is ~1.5MB).

### Variant Caller: Clair + Clair3 to Clair3 Only

The original pipeline used both **Clair** (deprecated) and **Clair3**. Since Clair3 subsumes Clair's functionality with improved deep neural network models specifically trained for HiFi data, MucOneSpan uses Clair3 exclusively.

### Read Filtering: Manual to Automated

The workflow described by Vrbacka et al. required **manual IGV inspection** to identify chimeric and incomplete reads. This is not reproducible and does not scale. MucOneSpan automates this by extracting reads per allele cluster and remapping them to the peak contig, eliminating the need for visual inspection.

### Allele Detection: peaks.py to Indel-Valley

The original `peaks.py` script (unpublished) detected allele lengths from the idxstats read count distribution. When alleles are far apart (>10 repeats), this works well. For close alleles, the read distributions overlap and simple peak-finding fails.

**Indel-valley splitting** is the key algorithmic innovation in MucOneSpan. It analyzes CIGAR indel lengths to identify two local minima in the indel series, resolving allele pairs as close as 3 repeats apart. See [Core Concepts](concepts.md#indel-valley-splitting) for details.

### Mutation Detection: Scripts to Template Catalog

The original mutation detection used unpublished Python scripts. MucOneSpan pre-computes the exact sequence of each known repeat type after each known mutation is applied. This enables **O(1) exact-match lookup** for the 13 most common mutations, with edit-distance fallback for novel mutations.

### Flanking Trim: Fixed to Anchor-Based

Fixed-position flanking trim assumes the flanking region is invariant. In practice, Clair3 sometimes calls false positive variants near the VNTR boundaries, causing indels in the consensus flanking region. **Anchor-based trim** searches for known boundary sequences (pre-repeat 1 start, after-repeat 9 end) and is resilient to these artifacts.

### Same-Length Alleles: Homozygous Assumption to Het Detection

When both alleles have the same repeat count, the workflow described by Vrbacka et al. treated them as homozygous and skipped Clair3. MucOneSpan instead runs Clair3 and examines **heterozygous genotype calls** to detect compound heterozygotes where both alleles are the same length but carry different mutations.

### Confidence Scoring: None to Per-Repeat Scores

The original pipeline did not report confidence. MucOneSpan assigns a score (0.0-1.0) to each classified repeat unit based on match quality. VCF cross-validation adjusts scores when Clair3 variants confirm or contradict the classification. This enables automated QC filtering.

---

## Next Steps

- **[Core Concepts](concepts.md)** -- Pipeline architecture details
- **[Known Mutations](../reference/mutations.md)** -- Full mutation catalog
- **[Benchmarking](../guides/benchmarking.md)** -- Validation results
