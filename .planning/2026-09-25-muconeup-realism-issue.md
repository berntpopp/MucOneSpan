# Draft MucOneUp issue: ONT R10 amplicon span-length shape and strand-specific error rates miss the PRJEB92208 targets

> Draft for the controller to file on berntpopp/MucOneUp. It has not been filed.
> Numbers are aggregates from the MucOneSpan MucSim-Bench dev pilots (Task 12b
> `v2`, Task 12c `v3`), measured with MucOneUp 0.45.0 against the public
> PRJEB92208 aggregate targets packaged in MucOneSpan
> (`src/muc_one_span/benchsim/targets/prjeb92208_v1.json`).

## Summary

Reads from the calibrated ONT R10 profiles (`ont_r10_sup_amplicon_v1`,
`ont_r10_genomic_v1`) reproduce the real **total** VNTR error rate well once
only calibrated-error cases are measured: the amplicon total is +4% of target.
Two gaps remain against the public PRJEB92208 aggregates:

1. **Span-length shape (amplicon).** The span-offset histogram of spanning
   reads has a Jensen-Shannon distance of **0.23 to 0.27** from the real one,
   against the benchmark's 0.1 target. Simulated spans are shifted towards
   shorter lengths and too narrow around the peak. There are too few reads
   0-15 bp longer than the allele, and almost no reads 1-4 units short for
   long alleles.
2. **Strand-specific error rates.** The per-type rates miss the per-strand
   targets although their totals match: deletions on the minus strand are
   **+53%** (deletions overall +29%), and mismatches on the plus strand are
   **-30%**. The genomic profile shows the same pattern, with minus-strand
   insertions +37% and minus-strand mismatches -32%.

The earlier estimate that the calibrated error is "about 25% too high" came
from a run in which half of the cases used a deliberately poor error level
(×1.5). It is **not** confirmed for the calibrated profile alone; see the table
below.

## How it was measured

- Scope, identical for simulation and targets: the VNTR interval (motif-1
  start to motif-9 end) of spanning reads. Error rates are the median of
  per-case rates, compared with the median of the per-library rates.
- Span offset = spanning-read VNTR span minus the true allele length, in 15 bp
  bins from -187 bp to +53 bp (open-ended outer bins), split at 55 units.
- `v3 standard`: 30 ONT amplicon cases at 500/1000/2000 reads and 30 ONT
  genomic cases at 30/60/100 spanning reads. All use calibrated error, PCR
  bias `none` or the profile default, and the profile's own smear 0.24 and
  chimera 0.023. HiFi has no public target.
- `v3 clean`: the same haplotypes at the top depth, with smear, chimera and
  PCR bias set to 0. This run still had the profile's concatemer rate (0.024)
  and off-target fraction (0.3).
- `clean2`: the 30 clean ONT amplicon cases regenerated with **no** molecule
  artefact (smear, chimera, concatemer and off-target all 0, no PCR bias).
- `v2`: the earlier mix (depth 5-2000; error calibrated or ×1.5 "poor"; PCR
  none, calibrated or strong; smear 0.05/0.25; chimera 0.01/0.05).

## Numbers

### ONT R10 amplicon (`ont_r10_sup_amplicon_v1`; target: 9 PRJEB92208 amplicon libraries)

| Metric | Real target (median) | v3 standard | v3 clean (concatemers left) | v2 (half poor error) |
| --- | --- | --- | --- | --- |
| Total error, all strands | 0.02138 | 0.02231 (+4%) | 0.02232 (+4%) | 0.02653 (+24%) |
| Total error + / - | 0.02687 / 0.01719 | 0.02556 / 0.01900 | 0.02576 / 0.01896 | 0.02904 / 0.02148 |
| Mismatch + / - | 0.00811 / 0.00606 | 0.00570 (**-30%**) / 0.00564 | 0.00570 / 0.00563 | 0.00731 / 0.00687 |
| Insertion + / - | 0.00723 / 0.00629 | 0.00762 / 0.00592 | 0.00766 / 0.00592 | 0.00889 / 0.00673 |
| Deletion + / - | 0.01138 / 0.00485 | 0.01233 / 0.00743 (**+53%**) | 0.01235 / 0.00743 | 0.01311 / 0.00783 |
| Deletion, all strands | 0.00767 | 0.00986 (**+29%**) | 0.00988 | 0.01101 |
| C7 correct length + / - | 0.5174 / 0.8876 | 0.5218 / 0.8791 | 0.5217 / 0.8792 | 0.5196 / 0.8758 |
| Span-offset JS distance, < 55 / >= 55 units | 0 (limit 0.1) | **0.232 / 0.270** | 0.462 / 0.450 | 0.209 / 0.247 |
| Off by > 1 unit (median of cases) | 0.2695 | 0.2686 | 0.0244 | 0.2054 |
| Between the alleles (median) | 0.0229 | **0.0** | 0.0 | 0.0 |
| Below the short allele (median) | 0.2361 | 0.2365 | 0.0 | 0.1220 |

The clean sets show that the JS gap is not caused by the smear model alone:
without smear products it grows to 0.45, because the real histogram includes
smear products. In `v3 clean`, the remaining off-by->1-unit share (0.0244) is
about the concatemer rate (0.024), because that run still had concatemers. The
artefact-free `clean2` run has an off-by->1-unit share of 0.0 and a JS distance
of 0.473 / 0.457. Its strand-specific rates are the same as in standard:
deletion + / - 0.01232 / 0.00745, mismatch + / - 0.00568 / 0.00566, total
0.02235. The strand-specific error gap therefore comes from the read error
model and not from molecule artefacts.

Normalized span-offset histogram, `v3 standard` vs real (fraction of spanning reads):

| Bin (bp) | < 55 u sim | < 55 u real | >= 55 u sim | >= 55 u real |
| --- | --- | --- | --- | --- |
| -67 | 0.0011 | 0.0136 | 0.0097 | 0.0224 |
| -52 | 0.0014 | 0.0216 | 0.0567 | 0.0328 |
| -37 | 0.0213 | 0.0686 | 0.1661 | 0.0840 |
| -22 | 0.2350 | 0.1813 | 0.2808 | 0.1801 |
| -7 (peak) | 0.4236 | 0.3315 | 0.1865 | 0.2564 |
| +8 | 0.0238 | 0.0928 | 0.0299 | 0.1302 |
| >= +53 | 0.0437 | 0.0149 | 0.0316 | 0.0095 |
| -172 to -82 (sum) | 0.0177 | 0.0235 | 0.0026 | 0.0430 |

### ONT R10 genomic (`ont_r10_genomic_v1`; target: 2 HG002 WGS libraries, tens of spanning reads)

This target is small and whole-genome, not targeted enrichment, so it is
indicative only.

| Metric | Real target | v3 standard | v2 |
| --- | --- | --- | --- |
| Total error, all strands | 0.00398 | 0.00429 (+8%) | 0.00466 (+17%) |
| Mismatch - | 0.00098 | 0.00067 (**-32%**) | 0.00081 |
| Insertion - | 0.00089 | 0.00122 (**+37%**) | 0.00142 |
| Insertion, all strands | 0.00111 | 0.00133 (+20%) | 0.00142 |
| Deletion + | 0.00229 | 0.00271 (+18%) | 0.00290 |
| Spanning fraction (median) | 0.4195 | **0.1732** | 0.3159 |

The low spanning fraction reflects the generic fragment-length model (median
6 kb, sigma 0.5), which is not fitted to a library.

### Allelic-ratio slope

The v3 standard slope is -0.022 per unit, against a real -0.056. This row is
**not evidence of a simulator defect**. The pooled slope mixes PCR levels, and
half of the standard cases use PCR bias `none` (slope 0) by design. v2, with
PCR levels none, calibrated and strong, gave -0.060. A per-level slope is not
reported yet.

## Suggested fix direction

1. **Per-strand error calibration.** Fit mismatch, insertion and deletion rates
   per read strand (`+` = motif1 to motif9 orientation), not only per type. The
   targets already carry `error_rates_per_ref_base["<strand>:<type>"]`. Today
   the simulated mismatch rate is nearly strand-symmetric (0.0057 / 0.0056),
   while the real one is not (0.0081 / 0.0061). Minus-strand deletions are the
   largest excess.
2. **Span-length noise.** Add a symmetric per-read span-length noise term
   (about ±1 bin, 15 bp) or calibrate the indel-length pmf so that the +8 bp
   bin and the -37 bp bin reach their real weights. Then re-check the JS
   distance against `span_offset_pmf_15bp_bins` (target <= 0.1 per size class).
   The deletion excess above probably explains part of the negative shift for
   alleles >= 55 units.
3. **Smear junction model.** No simulated reads fall between the two alleles
   (real median 2.3%), and long alleles have almost no products 1-4 units
   short. Fit `smear_junction_beta` and `smear_min_keep` to the real
   `offpeak_spanning_products.below.deleted_frac_of_parent` distribution
   (real p50 0.53), including products that keep most of the long allele.
4. **Concatemer or chimera tail.** Spans >= 53 bp longer than the allele are
   3× the real share (4.4% vs 1.5% below 55 units). Check `concatemer_rate`
   (0.024) and chimera placement against the `above` product counts.
5. **Genomic fragments.** Fit `fragments.length_median` / `sigma` to the WGS
   read-length quantiles, or document that the genomic profile's spanning
   fraction is not calibrated.
6. Minor, for batch use: `reads amplicon` and `reads ont --simulator
   pbsim3-fragments` have no `--threads` option. Their thread count comes from
   `pacbio_params.threads` / `ont_amplicon_params.threads`, default 8 (MucOneSpan
   now sets it through the profile `config_overrides`). A CLI flag would make
   the core budget explicit. A HiFi amplicon case with 1000 reads and a
   60-unit length difference needed about 27,000 templates and about 40
   minutes of `ccs` at 8 threads.

## Acceptance

With the calibrated profiles at their own default molecule rates, a
MucSim-Bench `standard` realism run should pass:

- every per-type and per-strand error rate within 20% of the target median;
- the span-offset JS distance at or below 0.1 for both size classes;
- `span_between_alleles_frac` within the real range.
