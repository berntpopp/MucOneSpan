# Platform-Aware Pipeline Design (Issue #16)

**Date:** 2026-04-07
**Issue:** berntpopp/MucOneSpan#16
**Goal:** Support ONT data by making Clair3 platform and minimap2 preset configurable via CLI

---

## Problem

MucOneSpan hardcodes `--platform=hifi` for Clair3 and `-x map-hifi` for minimap2. ONT data requires `--platform=ont` and `-x lr:hq` respectively. With HiFi settings on ONT data, sensitivity drops from 84.6% to 16.7%.

## Design

### CLI Options

1. **`--platform`** on `run` and `call` subcommands
   - Type: `click.Choice(["hifi", "ont"])`, default `"hifi"`
   - Threaded through to `call_variants_per_allele`, `disambiguate_same_length_alleles`, and `run_clair3`

2. **`--minimap2-preset`** on `run` and `map` subcommands
   - Type: `str`, default `None` (auto-selected from platform)
   - Auto-mapping: `hifi` -> `map-hifi`, `ont` -> `lr:hq`
   - If explicitly set, overrides the auto-selection

### Module Changes

**`mapping.py`:**
- `_run_mapping_pipeline()` gets `preset: str` parameter (replaces hardcoded `map-hifi`)
- `map_reads()` gets `preset: str` parameter, passes through
- Replace `-x map-hifi` with `-x {preset}` in minimap2 command

**`calling.py`:**
- `_extract_and_remap_reads()` gets `preset: str = "map-hifi"` parameter (replaces hardcoded `map-hifi` at line 122); passed through to the minimap2 `run_tool` call
- `call_variants_per_allele()` gets `platform: str = "hifi"` and `preset: str = "map-hifi"`, passes `platform` to `run_clair3()` and `preset` to `_extract_and_remap_reads()`
- `disambiguate_same_length_alleles()` gets `platform: str = "hifi"` and `preset: str = "map-hifi"`, passes `platform` to `run_clair3()` and `preset` to `_extract_and_remap_reads()`
- `run_clair3()` already accepts `platform` parameter -- no change needed (but callers must now actually pass it)

**`cli.py`:**
- `run` subcommand: add `--platform` and `--minimap2-preset` options, compute preset from platform if not overridden, pass `preset` to `map_reads()` and both `platform` + `preset` to `call_variants_per_allele()`
- `call` subcommand: add `--platform` and `--minimap2-preset` options, pass both to `call_variants_per_allele()`
- `map` subcommand: add `--minimap2-preset` option, pass to `map_reads()`
- Platform-to-preset mapping: `PLATFORM_PRESETS = {"hifi": "map-hifi", "ont": "lr:hq"}`
- Update main group docstring from "PacBio HiFi amplicon data" to "PacBio HiFi and ONT amplicon data"

### Testing

**Unit tests (no external tools):**
- `test_mapping.py`: Verify `_run_mapping_pipeline` passes correct `-x` preset to minimap2 subprocess mock (`map-hifi` default, `lr:hq` for ONT, custom override)
- `test_calling.py`: Verify `_extract_and_remap_reads` passes correct `-x` preset to minimap2 mock; verify `call_variants_per_allele` and `disambiguate_same_length_alleles` pass both `--platform=ont` to `run_clair3` mock and correct preset to `_extract_and_remap_reads`
- `test_cli.py`: Verify `run` and `call` accept `--platform ont` and `--minimap2-preset`, verify `map` accepts `--minimap2-preset`, verify preset auto-selection logic

**Integration tests (require external tools):**
- Generate ONT simulated test data via MucOneUp with ERRHMM-ONT model
- Run pipeline on ONT samples with `--platform ont`
- Verify allele detection, variant calling, and classification
- `@pytest.mark.integration` marker in `tests/integration/`

**Simulation / test data generation:**
- Update `scripts/generate_testdata.py` to generate ONT samples alongside HiFi
- Minimum ONT samples: `sample_ont_dupc_60_80`, `sample_ont_normal_60_80`, `sample_ont_dupa_60_80`
- Update `scripts/batch_analyze.py` to accept `--platform` flag

### Documentation

- **CLAUDE.md**: Add ONT section to testing strategy, document `--platform`
- **`docs/getting-started/quickstart.md`**: Add ONT usage example
- **`docs/reference/cli.md`**: Auto-generated from Click (picks up new options)
- **`docs/reference/limitations.md`**: Add ONT-specific limitations
- **CHANGELOG.md**: Entry under `[Unreleased]`
- **`.planning/TESTING_WITH_MUCONEUP.md`**: Add ONT data generation commands

### Files Touched

| File | Change |
|------|--------|
| `src/muc_one_span/cli.py` | Add `--platform`, `--minimap2-preset` options; update group docstring for ONT |
| `src/muc_one_span/mapping.py` | Add `preset` parameter to `map_reads` and `_run_mapping_pipeline` |
| `src/muc_one_span/calling.py` | Thread `platform` + `preset` through `_extract_and_remap_reads`, `call_variants_per_allele`, and `disambiguate_same_length_alleles`; fix existing `run_clair3` calls to pass `platform` |
| `tests/unit/test_mapping.py` | Test preset parameter |
| `tests/unit/test_calling.py` | Test platform parameter |
| `tests/unit/test_cli.py` | Test new CLI options |
| `scripts/generate_testdata.py` | Add ONT sample generation |
| `scripts/batch_analyze.py` | Add `--platform` flag |
| `docs/getting-started/quickstart.md` | ONT usage example |
| `docs/reference/limitations.md` | ONT limitations |
| `CLAUDE.md` | ONT testing docs |
| `.planning/TESTING_WITH_MUCONEUP.md` | ONT generation commands |
| `CHANGELOG.md` | Unreleased entry |

### Constraints

- Default behavior unchanged (`hifi` default for all options)
- No auto-detection of platform from BAM -- explicit `--platform` flag required
- Clair3 handles model selection internally based on `--platform`
- `lr:hq` requires minimap2 >= 2.26 (our pin is 2.28 -- satisfied)
