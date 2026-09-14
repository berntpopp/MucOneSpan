# Read-only adversarial runtime configuration audit

Scope: configuration integration in cli_settings.py, pipeline.py, cli.py,
calling.py/read_phasing.py, typed settings and scientific-stage consumers.
Inspected while coordinator integration was in progress on 2026-09-14. No source
or test files changed. All executed checks used temporary hand-built files,
Python API calls and mocks; no external tools, real pipelines, fresh data or
reserved seeds were used. Findings describe the inspected state, not a later
coordinator fix. The known min/max-vs-layout comparison issue was excluded from
new findings because the coordinator was already fixing it.

## Findings

### R1 — Configured short left flank can add a flank base to the VNTR

Priority: P1 scientific correctness for accepted nondefault configurations.
Locations: `ladder.py:45` and `consensus.py:161`.

Ladder generation takes `repeat_dict.flanking_left[:flank_length]`, but consensus
anchoring takes the last anchor_bases from the entire dictionary left flank.
When configured flank_length is shorter than that sequence, the anchor does not
match the reference generated with the same configuration. A leading flank indel
then survives fixed-position trimming as a spurious VNTR base.

Deterministic reproduction (uses the real builder/trimmer and no tools):

```python
rd = replace(load_repeat_dictionary(),
    repeats={"first": "GTTG", "X": "ACGT", "last": "CACA"},
    flanking_left="AAAACCCC", flanking_right="GGGGTTTT")
layout = ReferenceLayoutSettings(("first",), ("last",))
settings = ConsensusSettings(flank_length=4, anchor_bases=2, anchor_tolerance=5)
sequence = build_contig(1, rd, settings=settings, reference_layout=layout)["sequence"]
# sequence == "AAAAGTTGACGTCACAGGGG"
consensus.write_text(">c\nT" + sequence + "\n")  # insertion before the original left flank
trim_flanking(consensus, None, output, repeat_dict=rd, settings=settings,
              reference_layout=layout, context=context)
```

Observed VNTR `AGTTGACGTCACA`; expected `GTTGACGTCACA`.
Observed left trim4 / `fixed_anchor_not_found`, right trim17 / `exact_anchor`.
The correct left boundary is5. Proposal: derive the left anchor from the exact
selected flank sequence used by ladder generation, with validated anchor extent;
add this shifted-boundary fixture plus unchanged/default-flank parity controls.
Merely validating anchor_bases <= unit length does not address this defect.

### R2 — Standalone consensus ignores configured dictionary and anchor behavior

Priority: P2 accepted configuration silently ignored.
Location: `cli.py:357` call to build_consensus_per_allele, and
`consensus.py:157` gating exact anchors on repeat_dict being present.

The standalone command forwards consensus/reference_layout settings but never
loads or passes the repeat dictionary. Therefore a custom dictionary and both
anchor parameters are accepted but cannot affect anchor matching. Layout ID
membership is also not checked along this path. Full run does supply the
configured dictionary, so standalone and full-pipeline behavior diverge.

Executed CliRunner with a valid config containing `repeat_dictionary=custom.json`
(where custom.json deliberately contained invalid `{}`) and anchor_bases5. Mocked
only tool availability and build_consensus_per_allele. Exit0; captured kwargs
were `settings=ConsensusSettings(...anchor_bases=5...)` and reference_layout,
with no repeat_dict. The configured dictionary was never parsed.

Proposal: load and validate the configured dictionary before tools and pass it
through, or explicitly reject/declare unsupported anchor/dictionary settings for
this command. The advertised shared configuration favors actual forwarding.

### R3 — Explicit stage CLI values bypass typed settings validation

Priority: P2 boundary validation and failure semantics.
Locations: `cli.py` ladder/map/alleles/call numeric Click options and callbacks;
`cli_settings.py:16` only validates the file before Click applies overrides.

Real run reconstructs RunSettings from effective values, but stage commands pass
explicit integers/floats straight to their implementation. Executed examples:

- `map --threads 0`: exit0 with mocked map_reads; captured threads0.
- `ladder --min-units 20 --max-units 10`: exit0 with mocked generator; captured
  `(20,10,500)`, despite ReferenceLayoutSettings rejecting that range.

Without mocking, reversed ladder bounds produce an empty range/output rather
than a valid ladder. Similarly negative flank length, zero coverage or nonfinite
QUAL can bypass the corresponding constructor guard via explicit CLI values.

Proposal: validate effective per-command RunSettings/ConsensusSettings/layout
using dataclasses.replace, before resource/tool execution. Preserve intentional
explicit CLI precedence while applying the same invariants as JSON. Add
CliRunner regressions asserting rejection and external-tool mocks not called.

### R4 — Failed rerun retains stale effective configuration provenance

Priority: P2 reproducibility and stale-artifact safety.
Locations: `pipeline.py:43` early file validation, `cli_settings.py:69` provenance
written later, and `run_status.py` invocation wrapper.

Prepared an existing output directory with `run_configuration.json` equal to
`{"old":true}`, then invoked `run --input missing.fq --output-dir OUT`.
Observed exit2 and a fresh execution_failed run_status.json, but the previous
run_configuration.json remained byte-for-byte unchanged. Configured-layout or
dictionary validation failures before write_run_configuration have the same risk.

Proposal: replace/invalidate prior configuration provenance at invocation start,
record effective settings plus explicit unavailable hash/error fields as they
become known, and write atomically. A failed status is correctly recorded, but
it must not be paired with an unlabeled configuration from another attempt.
Global JSON parse failures occur before the run callback and should also have a
clearly defined stale-output policy; do not imply those started a pipeline.

### R5 — Zero/narrow anchor tolerance searches the wrong coordinate

Priority: P2 accepted parameter does not implement its documented meaning.
Locations: `consensus.py:82` and calls around164/173.

_find_anchor searches possible anchor START positions around expected_pos.
trim_flanking supplies the VNTR boundary, while the actual anchor starts
anchor_bases earlier. With tolerance0, even an unchanged exact anchor is missed.

Executed fixture: repeats first=GTTG/last=CACA, flank strings AAAA/GGGG,
sequence `AAAAGTTGCACAGGGG`, flank_length4, anchor_bases2, anchor_tolerance0.
Both context methods were `fixed_anchor_not_found`; expected exact anchors at
the unshifted boundaries. Fixed trimming hides the defect on this exact input,
but search tolerance is effectively asymmetric for shifted inputs.

Proposal: define tolerance around an explicitly documented anchor-start or
boundary coordinate and convert consistently at call sites; test tolerance0,
positive/negative shifts at the allowed boundary, and default-output parity.
This is distinct from overlong-anchor validation and R1's selected-flank mismatch.

### R6 — Actual WhatsHap command is not persisted in returned metadata

Priority: P2 provenance requirement/documentation mismatch.
Locations: `read_phasing.py:164` metadata construction and `:212` run_tool call.

Executed the real phase helper with the existing test tool mocks, two unphased
heterozygous input records, phased output records and valid two-record read list.
Configured overrides were internal_downsampling20/mapping_quality5. The helper
returned status phased and correctly recorded version/overrides, but no command
or phase_command field (nor any command string anywhere in serialized metadata).
Captured executed argv began:

```text
whatshap phase --indels --ignore-read-groups --internal-downsampling 20 --mapping-quality 5
```

Reproduction uses tests/unit/test_read_phasing.py helpers loaded with runpy,
pytest.MonkeyPatch.context, input named `input.vcf`, and mock_tools output_records
`records("0|1","0|1",phase_sets=["1","1"])`.

Proposal: construct argv once, store an immutable/copy list in metadata (and an
attempt record if failures need per-attempt provenance), then pass that exact
list to run_tool. Do not rely on optional DEBUG logs to fulfill the persisted
command claim in the spec/guide.

## Additional observations and positive controls

- Explicit CLI threads7/platformhifi/--no-report correctly override configured
  threads2/platformont/reporttrue (real Click parsing, execute_pipeline mocked).
- Experimental read phasing remains false by default, optional None overrides
  add no flags, and explicit overrides are forwarded. Existing phase/record
  identity guards were retained in the inspected helper.
- Default numerical values remain equal to the approved values. However,
  operational defaults are still duplicated as literals in mapping.map_reads
  (threads4), calling._extract_and_remap_reads/disambiguate/call_variants_per_allele
  (threads4, QUAL5, platformhifi), and alleles.detect_alleles (coverage10).
  This is a P3 maintainability gap against the single-source-default requirement;
  replacing literals with central DEFAULT_SETTINGS references preserves behavior.
  Fixed genotype/ploidy/coordinate constants are not part of this concern.
- Thread semantics are per-tool/concurrent-stage values, not a strict global CPU
  cap: two independent alleles with threads1 still receive remap threads1 and
  Clair3 threads1 each; ThreadPoolExecutor has two workers. Mapping also runs
  minimap2 and samtools sort concurrently with the configured value. This behavior
  predates configuration; document it and avoid claiming a global thread budget.
  Any scheduling change should preserve the frozen comparison contract or be
  separately reviewed, not silently bundled with configuration.
- No finding claims scientific accuracy from these small fixtures. They test
  deterministic forwarding, trimming and provenance behavior only.

## Execution summary

Three temporary Python scripts exercised real Click, real settings, real
builder/trimmer, existing read-phasing mock fixtures and actual in-process
calling orchestration with mocked external stages. All reproduced outputs are
reported above; no source edits or production-data inspection occurred. The
coordinator received findings as they were reproduced, before this handoff.
