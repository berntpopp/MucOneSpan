# Known Mutations

MucOneSpan includes a catalog of **13 known MUC1 frameshift mutations** with pre-computed sequence templates for exact matching.

---

## Mutation Catalog

| Mutation | Repeat(s) | Effect | Citation |
|----------|-----------|--------|----------|
| **59dupC** | 5, X, B, D, K, W, M, N, O, P, Q, R | +1bp (C duplication in 7C tract) | Kirby et al. 2013 |
| **60dupA** | X | +1bp (A duplication) | Olinger et al. 2020 |
| **58_59insG** | X | +1bp (G insertion) | Olinger et al. 2020 |
| **56_59dupCCCC** | X | +4bp (CCCC duplication) | Vrbacka et al. 2025 |
| **18_31del** | X | -14bp (deletion) | Vrbacka et al. 2025 |
| **42_57dupGGGCTCCACCGCCCCC** | C | +16bp (16bp duplication) | Vrbacka et al. 2025 |
| **31ins25bp** | A, B, J, K, N, S, X | +25bp (25bp insertion) | Saei et al. 2023 |
| **delinsAT** | X | net +1bp (retained boundaries 54..56, delete 55, insert AT) | Olinger et al. 2020 |
| **delGCCCA** | multiple | -5bp (5bp deletion at start) | Saei et al. 2023 |
| **insC_pos23** | A, E | +1bp (C insertion at position 23) | Vrbacka et al. 2025 |
| **insG_pos58** | B, X | +1bp (G insertion at position 58) | Vrbacka et al. 2025 |
| **insG_pos54** | B, J | +1bp (G insertion at position 54) | Vrbacka et al. 2025 |
| **insA_pos54** | A, H | +1bp (A insertion at position 54) | Vrbacka et al. 2025 |

---

## How Mutation Detection Works

### Step 1: Exact Match

Each 60bp unit from the consensus sequence is compared against all 34 known repeat type sequences. If an exact match is found, the repeat is classified without mutation.

### Step 2: Mutation Template Probing

If no exact match, the unit is compared against **pre-computed mutation templates** (across 13 mutations, each applied only to its allowed repeat types). Each template is the exact sequence a repeat type would have after a specific mutation is applied.

```
Example: X repeat (60bp) + dupC mutation = X:dupC template (61bp)
```

The `dupC` templates cover the exact seven-C tract at repeat positions 53–59,
preceded by G and followed by a non-C base, in each listed 60-base repeat. A
single additional C in this tract consumes a 61-base unit and retains its parent
label, for example `B:dupC`. Other insertion positions, bases, and tract contexts
do not receive this exact template match. This sequence-based applicability,
including repeat 5, does not imply a clinically observed carrier for every listed
repeat background. Exact template agreement is separate from VCF support and
diagnostic validation.

!!! note "O(1) lookup"
    For the 13 known mutations, detection is an O(1) hash lookup against the template dictionary. No edit distance computation is needed for common cases.

### Step 3: Edit Distance Fallback

If no template matches, the pipeline computes **Levenshtein edit distance** against all repeat types and reports the closest match as a potential novel mutation.

Novel mutations are still detected but with **lower confidence scores** since they do not match any known template.

---

## Clinical Significance

All 13 cataloged mutations cause **frameshift** in the MUC1 coding sequence, leading to production of the toxic MUC1fs protein. The most common mutation is **59dupC** (also known as dupC), found in the majority of ADTKD-MUC1 families worldwide.

!!! warning "Diagnostic use"
    MucOneSpan is a research tool. Clinical diagnostic use requires validation against established reference standards and is subject to local regulatory requirements.

---

## Adding New Mutations

The mutation catalog is defined in `data/repeats/repeats.json`. To add a new mutation:

1. Define the mutation operation (insert, delete, or substitution)
2. Specify which repeat types it can affect (`allowed_repeats`)
3. The pipeline will automatically pre-compute templates for all affected repeat types

---

## References

- Madritsch et al. Analysis of clinically relevant large tandem repeats using nanopore sequencing. *Sci Rep.* 2025. [doi:10.1038/s41598-025-30441-3](https://doi.org/10.1038/s41598-025-30441-3).
- Kirby A, et al. Mutations causing medullary cystic kidney disease type 1 lie in a large VNTR in MUC1 missed by massively parallel sequencing. *Nat Genet.* 2013;45(3):299-303.
- Olinger E, et al. Clinical and genetic spectra of autosomal dominant tubulointerstitial kidney disease due to mutations in MUC1. *Kidney Int.* 2020;98(2):473-487.
- Saei H, et al. MUC1 VNTR frameshift mutations. *Kidney Int Rep.* 2023.
- Vrbacka A, et al. Long-Read Sequencing of the MUC1 VNTR. *bioRxiv.* 2025. doi:10.1101/2025.09.06.673538.

---

## Next Steps

- **[Repeat Nomenclature](nomenclature.md)** -- Classification system details
- **[Core Concepts](../getting-started/concepts.md)** -- How mutation detection fits in the pipeline
- **[Benchmarking](../guides/benchmarking.md)** -- Validation with simulated mutations
