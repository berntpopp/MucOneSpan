# Phasing parameter research and debugging protocol

2026-09-14; development evidence only, reserved final seeds unused. User requested
root-cause debugging and research of ideas, best practices and parameters before
abandoning the failed phase promotion. No parameter selected from final data.

## Primary sources and version constraints

- [WhatsHap official guide](https://whatshap.readthedocs.io/en/latest/guide.html),
  accessed 2026-09-14. Phase augments genotype calls; reference-backed realignment
  is recommended for noisy reads. Common PS describes connectivity, not evidence
  of two fully recovered VNTRs. Single-sample data may use explicit --sample and
  --ignore-read-groups. Structural differences are outside ordinary phase support.
- [WhatsHap source](https://github.com/whatshap/whatshap) and installed1.7
  `whatshap phase --help` captured in ignored phasing_debug/whatshap-1.7-help.txt.
  This version defaults to indels OFF, mapq20, internal-downsampling15; reference
  and --indels are explicit. Increasing downsampling coverage has exponential
  runtime cost. Read merging is OFF; its error parameters do nothing without
  --merge-reads. Keep it off unless platform-specific evidence justifies it.
- [Clair3 official repository](https://github.com/HKU-BAL/Clair3), accessed
  2026-09-14, cross-checked installed run_clair3.sh1.0.10. Distinguish intermediate
  calling phase from final VCF phase. Installed defaults: minimum mapping quality5;
  indel candidate AF0.08 HiFi/0.15 ONT; phasing variant fraction0.7 except specified
  ONT models0.8. Haploid precise retains only1/1; haploid sensitive also retains0/1.
  These modes do not establish molecular allele identity and must not be applied
  to a mixed diploid cluster as a shortcut.

The current web guide and installed1.7 CLI differ on distrust-genotype behavior:
local1.7 exposes --include-homozygous to allow homozygous-to-heterozygous changes,
whereas guide prose says only heterozygous loci are considered. Test the actual
binary rather than assume newer/stale prose defines local behavior. The existing
production helper deliberately rejects genotype changes; any genotyping experiment
must remain separate until its genotype/likelihood contract is specified.

## Causal observations before parameter changes

In sample_dupc_60_80__n20_seed1701, original20-record input contains14 terminal
spans3557–3615bp, five4762–4804bp, and one partial-anchor record. These are observed
sequence spans, not simulator source labels. The selected cluster contains only
contigs49–53 and Clair3 site depth13. This suggests majority-only extraction,
not a genuine equal-length diploid input; verify retained records by full sequence
and quality, never by QNAME alone. A name collision crosses the two span groups.

If the mutation is already1/1 in a majority-only cluster, genotype-preserving
phasing necessarily puts it in both genotype selections. Raising phasing coverage
cannot recreate an omitted long allele. One PS is insufficient proof of diploid
multiplicity in a potentially single-allele read cluster.

## Bounded development ablations

1. Track input record -> ladder primary/secondary fit -> extraction -> remap ->
   phasing selection. Count exclusions, source identity ambiguity and span groups.
2. Inspect mutation REF/ALT/GT/AD, heterozygous SNPs, phase assignments and exact
   sample-level event changes. Genotype identity and variant phase are distinct.
3. Test default15 vs10/20 internal coverage and mapq20 vs5 on the failure and
   synthetic linked/disconnected controls. Do not increase to full200read coverage.
4. Test distrust-genotypes alone and with include-homozygous as diagnostic
   genotype-revision experiments; preserve old/new calls and likelihood evidence.
5. Test full-input molecule/length-aware assignment if extraction proves causal.
   A targeted fix must retain true equal-length cis/trans/shared variants and
   minority/close-length controls, not suppress the second output to erase FP.

Acceptance remains unchanged: preserve exact TP and prior exact reconstructions,
no extra false calls, no apparent improvement achieved through lower call rate.
Report a failed ablation plainly. Keep no-call reduction, phase block count and
exact reconstruction as separate endpoints. No global QUAL reduction or
mutation-template-anywhere scan is justified by these sources.
