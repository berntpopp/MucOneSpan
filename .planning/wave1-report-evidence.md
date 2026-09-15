# Wave 1 report design and evidence

Scope: issues #39 and #40; report.py, report_igv.py, template and focused tests.
Interfaces agreed with lead: optional plural `vcf_paths: dict[str, Path]` takes
precedence over singular including an empty dict; `execution_status: dict | None`
is explicit authoritative provenance, with summary.run_status fallback only.
Full pipeline passes `analysis_completed` after successful analysis; execution
sidecar stays running until callback/report succeeds. Known noncompleted states
prevent a negative call. Mutation records remain visible with execution warning.
Legacy missing status retains prior decisions with an explicit unavailable label.

Contigs come from sorted summary allele keys, using contig_name, deduplicated.
Explicit missing contigs and requested missing files raise descriptive errors.
Absent optional paths remain optional. With no allele metadata, historical first
reference fallback remains. VCF track configurations supply stable allele labels;
shared resolved file paths form one track explicitly labeled shared, never two
independently resolved haplotypes. Preserve existing locus coordinate arithmetic.

Percent units remain 0–100. Report boundary accepts finite numeric percentages in
range (excluding bool); invalid/missing values render Unicode dash without a
misleading progressbar. Nomenclature placeholders use Unicode with autoescape.

## Red evidence before edits

`uv run --locked --all-extras pytest tests/unit/test_report_wave1.py --no-cov -q`:
17 failed, 5 passed. 0.5 displayed 50%; 1 displayed 100%; invalid metrics crash or
produce misleading bars; failed/interrupted/running status became negative;
explicit status unsupported; unknown mutation repeat form emitted `&amp;mdash;`;
requesting missing contig_11 silently selected contig_1.

Producer evidence: classification_summary.py computes exact_count / repeats *100;
tests/unit/test_classify.py asserts exact_match_pct == 100.0 for exact matches.
Live issue #39/#40 bodies read via gh api on 2026-09-15; both have zero comments.

## Test matrix / execution steps

- [x] Reproduce percent, escaping, status and missing-contig failures.
- [x] Reproduce real session wrong-contig and dropped-track behavior.
- [x] Implement minimal report plumbing and template fixes.
- [x] Test one/two allele contigs, shared VCF, singular compatibility, optional
  paths, missing paths/contigs and plural precedence; real embedded/sidecar sessions.
- [x] Test DOM loci/tracks offline in Chrome, along with escaping and percentages.
- [x] Run focused tests/lint; lead combines lifecycle tests and full gates/review.

## Adjacent evidence policy (not changed)

Current compute_clinical_decision labels every mutation record pathogenic, even
records with frameshift=False/vcf_support=False; it also allows an empty legacy
summary to become negative. These are broader biological interpretation policy,
kept separate from provenance plumbing. Reproduction: compute_clinical_decision({})
returns NO_PATHOGENIC_VARIANT_DETECTED, and a summary classifications entry with
mutations=[{'mutation_name':'unknown','frameshift':False,'vcf_support':False}]
returns PATHOGENIC. No thresholds, clinical tiers or mutation science changed.


## Final worker verification

Real baseline `tests/integration/test_report_sessions.py`: 4 failed before fixes.
Both embedded/sidecar reported `contig_1` instead of `contig_11`, `contig_21`;
plural VCF argument did not exist. Explicit BED missing-contig regression also
failed (DID NOT RAISE), then passed with reference validation.

`PATH=<vntyper-env>/bin:$PATH uv run --locked --all-extras pytest
 tests/unit/test_report.py tests/unit/test_report_wave1.py
 tests/unit/test_report_igv.py tests/integration/test_report_sessions.py
 tests/browser --no-cov -q --tb=short`: **97 passed, 0 skipped** (74 unit,
12 real create_report session checks, 11 Chrome/Playwright browser checks).
The available vntyper environment supplies igv-reports **1.13.0**; external tools
are invoked through the existing run_tool abstraction. Google Chrome is available.
No repository dependencies or bundled JavaScript assets changed.

Actual VCF payloads preserve the exact contig, 1-based position 700, A>C, QUAL 60
record. Browser tests click both visible loci and verify both track labels, the
visible locus search value, and loaded variant features (0-based start 699,
end 700, ref A, alt C). Both modes run offline with zero external requests,
console errors or page exceptions. Six percent cases agree in visible text,
CSS width and ARIA. Three original browser tests also pass.

Ruff lint and formatting passed on all worker-authored Python files. Report.py
and report_igv.py remain below 400 lines; report.html.j2 is 547 lines.
Lead owns combined mypy/file-size/coverage/build checks and independent review.

## Upstream igv-reports limitation, detected rather than masked

Additional payload validation uncovered a real regression in the available base
environment's igv-reports **1.16.0**. VcfReader.slice removes the terminal newline
from the VCF header then appends the first variant. Example malformed payload:
`#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFOcontig_11\t700...`.
Initially labels-only browser tests passed, but decoded-record assertions failed
in both modes (2 failed, 10 passed in the strengthened session suite).

Report generation now compares embedded #CHROM headers with the original VCF
header, including sample columns. A deterministic malformed sample-column test
failed before validation and passed afterward. Real 1.16 reproduction now raises
`ValueError: Malformed VCF track 'Allele 1 variants' from create_report: column
header changed. Use a working igv-reports version (1.13.0 verified); 1.16.0 can
concatenate the first variant into its header.` Both new main HTML and sidecar
remain absent. Scientific VCF data are never rewritten to repair external output.
Use the verified 1.13.0 external tool environment for full browser/session gates.

Malformed status list/dict and an integer percentage 10**1000 additionally had
3 reproduced failures (unhashable status / OverflowError), now covered as
inconclusive status / unavailable metric. Direct execution confirmed both
adjacent biological-policy reproductions above; they remain intentionally
unchanged and require separate scientific policy work.
