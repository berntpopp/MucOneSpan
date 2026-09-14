# Central API default references

Bounded implementation changed mapping.py and calling.py default references and
only the min_coverage signature default in alleles.py. Added
`tests/unit/test_api_defaults.py`. No scheduling, scientific filtering, genotype,
reference selection or subprocess behavior was redesigned.

## Change

- Mapping/remapping thread defaults use DEFAULT_SETTINGS.run.threads.
- Calling model/platform/thread/QUAL defaults reference their central RunSettings
  fields. QUAL uses float(DEFAULT_SETTINGS.run.min_qual), preserving the legacy
  float default and default command spelling rather than changing5.0 to5.
- detect_alleles uses DEFAULT_SETTINGS.run.min_coverage for its default.
- mapping.py owns the platform-preset mapping (HiFi map-hifi; ONT lr:hq) and exports
  DEFAULT_MINIMAP2_PRESET, derived from an explicit central preset or the central
  platform. Mapping and calling signatures retain a concrete default string.
- Coordinator will reuse the exported mapping from CLI rather than maintain a
  duplicate route table. Explicit configured/CLI/function presets still win;
  this patch does not add automatic per-call platform inference to old APIs.

All current signature values remain4threads, coverage10, QUAL5.0, model empty,
platformhifi and presetmap-hifi. Runtime JSON choices are still passed explicitly
by integration; no module-global mutable runtime settings were introduced.

## Tests

Initial default-value contract passed on the old code. Two additional tests
executed isolated module definitions with a substituted immutable central
settings object (threads7, coverage12, QUAL8.5, modelreviewed_model, platformont).
Both failed before implementation because map_reads retained its literal4.
After the fix, both passed, including automatic ONT presetlr:hq and an explicit
map-pb preset. Tests do not reload live production modules; monkeypatch restores
all temporary settings/module bindings. No external tools are called.

Final checks on 2026-09-14:

```text
uv run --locked --all-extras pytest -q --no-cov tests/unit/test_api_defaults.py tests/unit/test_mapping.py tests/unit/test_calling.py tests/unit/test_alleles.py
80 passed
uv run --locked --all-extras mypy src/muc_one_span/mapping.py src/muc_one_span/calling.py src/muc_one_span/alleles.py
Success: no issues found in 3 source files
uv run --locked --all-extras ruff check src/muc_one_span/mapping.py src/muc_one_span/calling.py src/muc_one_span/alleles.py tests/unit/test_api_defaults.py
All checks passed
uv run --locked --all-extras ruff format --check src/muc_one_span/mapping.py src/muc_one_span/calling.py src/muc_one_span/alleles.py tests/unit/test_api_defaults.py
4 files already formatted
```

Sizes at handoff: mapping195, calling445, alleles596, API tests104 lines; all below650.
No fresh validation inputs, seeds, simulator, real tools or pipelines were used.
No commits were created. Global CI remains the coordinator's integration gate.

## Remaining semantics and limits

`threads` is a per-tool/concurrent-stage setting, not a strict total CPU cap.
Initial mapping runs minimap2 and samtools sort concurrently with the configured
value; samtools -@ denotes additional threads. Independent allele remappings run
in two executor workers, each receiving the full configured remapping value.
Clair3 divides the requested count across the selected alleles with a minimum1
per caller; therefore threads1 with two alleles still permits two caller workers.
The fixed two-allele scheduling policy and external-tool internal threading are
unchanged. No new executor-count or process-tree-budget parameter was invented.

The legacy min_dp argument remains accepted but unapplied by VCF filtering; it
is not exposed as a functioning JSON depth-filter knob. Existing min_qual remains
the actual configurable QUAL cutoff. Alignment/VCF format flags, diploid allele
cardinality and genotype/coordinate contracts remain fixed behavior rather than
new runtime parameters.
