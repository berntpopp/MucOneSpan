# Repeat Nomenclature

MucOneSpan uses the **Vrbacka nomenclature** (extended from Kirby et al. 2013 and Wenzel et al. 2018) to classify each 60bp repeat unit in the MUC1 VNTR.

---

## Repeat Categories

### Pre-Repeats (N-terminal)

Pre-repeats occur at the beginning of the VNTR array, before the variable-length canonical region.

| Symbol | Description |
|--------|-------------|
| `1` | Pre-repeat 1 (first unit of the VNTR) |
| `2` | Pre-repeat 2 |
| `3` | Pre-repeat 3 |
| `4` | Pre-repeat 4 |
| `4p` | Pre-repeat 4 variant (polymorphic) |
| `5` | Pre-repeat 5 |
| `5C` | Pre-repeat 5 variant |

### Canonical Repeats (Variable Region)

The core VNTR consists of canonical repeat types that differ by 1-3 nucleotide substitutions from the consensus **X** repeat.

| Symbol | Description |
|--------|-------------|
| `X` | Canonical repeat (most common; consensus sequence) |
| `A`, `B`, `C`, `D`, `E`, `F`, `G`, `H`, `I`, `J`, `K`, `L`, `M`, `N`, `O`, `P`, `Q`, `R`, `S`, `V`, `W` | Variant types differing by 1-3 substitutions from X (22 types total including X) |

!!! note "X is the canonical repeat"
    The **X repeat** is the most common type in the variable region. The 21 other canonical repeat types (A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, V, W) are defined by their specific nucleotide differences from X. Note that T, U, Y, and Z are not used.

### After-Repeats (C-terminal)

After-repeats occur at the end of the VNTR array, following the variable-length canonical region.

| Symbol | Description |
|--------|-------------|
| `6` | After-repeat 6 |
| `6p` | After-repeat 6 variant (polymorphic) |
| `7` | After-repeat 7 |
| `8` | After-repeat 8 |
| `9` | After-repeat 9 (final unit of the VNTR) |

---

## Mutation Naming Convention

Mutations are named using the format `RepeatType:MutationDescription`:

```
X:59dupC      -- dupC mutation in an X repeat
C:42_57dup... -- 16bp duplication in a C repeat
A:insC_pos23  -- C insertion at position 23 in an A repeat
```

**Components:**

- **Repeat type** -- which canonical type the mutated unit derives from (e.g., X, C, A)
- **Position** -- nucleotide position(s) within the 60bp repeat unit
- **Operation** -- `dup` (duplication), `ins` (insertion), `del` (deletion), `delins` (deletion + insertion)

---

## Example VNTR Structure

A complete allele structure reads left-to-right (N-terminal to C-terminal):

```
1 2 3 4 5 C F X X X X:dupC B A A B X X X V 6 7 8 9
|---------|   |---------------------------|   |-------|
Pre-repeats   Variable region (canonical)     After-repeats
              with one mutated repeat
```

**Reading this structure:**

- Starts with pre-repeats `1 2 3 4 5`
- Variable region contains canonical types (`C F X X X ... V`)
- One repeat is mutated: `X:dupC` (an X repeat carrying the 59dupC frameshift)
- Ends with after-repeats `6 7 8 9`

---

## Confidence Scores

Each classified repeat receives a confidence score:

| Score Range | Meaning |
|-------------|---------|
| **1.0** | Exact match to known repeat type or mutation template |
| **0.9-0.99** | Very close match (1 substitution) |
| **0.8-0.89** | Close match (2-3 substitutions) |
| **< 0.8** | Low confidence -- possible sequencing error or novel variant |

!!! tip "Interpreting low confidence"
    Repeats with confidence < 0.8 may indicate sequencing errors, novel mutations not in the catalog, or alignment artifacts. Check the VCF for supporting variant calls.

---

## References

- Kirby A, et al. Mutations causing medullary cystic kidney disease type 1 lie in a large VNTR in MUC1 missed by massively parallel sequencing. *Nat Genet.* 2013;45(3):299-303.
- Wenzel A, et al. Single molecule real time sequencing in ADTKD-MUC1 allows complete assembly of the VNTR and detection of sec7A mutation. *Sci Rep.* 2018;8:4170.
- Vrbacka A, et al. Long-Read Sequencing of the MUC1 VNTR. *bioRxiv.* 2025. doi:10.1101/2025.09.06.673538.

---

## Next Steps

- **[Known Mutations](mutations.md)** -- Full mutation catalog
- **[Core Concepts](../getting-started/concepts.md)** -- How classification fits in the pipeline
