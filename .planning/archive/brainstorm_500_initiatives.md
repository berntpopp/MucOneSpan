You are an expert computational genomicist and peer reviewer specializing in tandem repeat bioinformatics, specifically MUC1 VNTR and ADTKD-MUC1 diagnostics. We are planning the next major architectural milestone for MucOneSpan (a Python package analyzing MUC1 VNTR from PacBio HiFi and ONT sequencing).

We have 5 core initiatives to design:

1. **Adaptive Sampling / Genomic Simulation**:
   Expanding beyond amplicon PCR (which suffered length-dependent dropouts due to PCR bias) by integrating native genomic simulation via MucOneUp (NanoSim for ONT, PBSIM3+CCS for PacBio). Native reads model 1:1 stoichiometric representation, unamplified HMW DNA, and staggered genomic flanks extending 10 kb into unique sequence. How should our benchmark runner support both modes? What are the key simulation parameters (coverage 15x-60x, read length distributions, error models)?

2. **~500-Dataset Stratified Simulation Catalog**:
   Designing ~500 datasets across dev (e.g. 300), validation (e.g. 100), and held-out test (e.g. 100) splits.
   - What strata should be defined (e.g. Amplicon PCR vs Native Genomic, HiFi vs ONT)?
   - What biological challenge profiles should be prioritized (extreme length asymmetry like 25 vs 140 repeats, homozygous identical, 1-RU gaps like 60/61, diverse repeat compositions)?
   - What pathogenic variants (canonical 59dupC, 56_59dupCCCC, 58_59insG, 60dupA, 54_56delinsAT, 1_5delGCCCA, etc.) and benign polymorphisms?
   - What should the dual cryptographic ledger design look like to prevent data leakage and guarantee reproducibility?

3. **IUPAC Ambiguity / Consensus Replay Root Cause**:
   On single-allele partitioned BAMs (where all reads belong to one physical allele), Clair3 called minor-noise/chimeric variants as 0/1 (heterozygous, AF ~0.70-0.85). This caused `bcftools consensus` to inject IUPAC codes (S, M, R, Y), resulting in 95.3% of length-exact sequence mismatches and '?' repeat classifications in classify.py.
   - How should we address this architecturally?
   - E.g., haploid-forced consensus replay, major-allele VCF filtering/sanitizing (mapping AF >= 0.60 to 1/1 and dropping AF < 0.40 noise), or consensus sequence post-processing? What is the cleanest, most defensible approach?

4. **HGVS Nomenclature Alignment & Clinical Reporting**:
   MUC1 is on chromosome 1 on the negative strand (transcription is right-to-left in genomic coordinates).
   - Review HGVS recommendations (specifically the 3'-rule for cDNA NM_001204286.1 vs genomic vs repeat-unit relative shorthand like 59dupC and 53C[7]>53C[8]).
   - Note that HGVS restricts repeat notation `c.seq[N]` in coding regions to multiples of 3, requiring duplications/insertions for frameshifting indels.
   - How should MucOneSpan format, validate, and present these variants in clinical reports (e.g. standard HGVS, repeat form, ambiguity interval, literature reference)? Compare with VNtyper's nomenclature approach.

5. **Interactive IGV Report Integration**:
   Embedding self-contained or sidecar igv.js visualization in HTML reports (referencing ../VNtyper).
   - What are best practices for locus windows, track configuration (BAM, VCF, Reference FASTA), and zero-CDN / offline self-contained HTML packaging?

Please provide your rigorous, creative, structured brainstorm and critical analysis for each of these 5 areas.
