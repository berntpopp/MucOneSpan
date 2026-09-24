# P0 Clinical Safety Fixes (v0.16.0) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the known clinical false negative (MP4 dupC reported NEGATIVE) and the evidence-exceeding decisions of the current ladder engine. Keep CLI flags and output schemas, adding fields only, and release as v0.16.0.

**Architecture:** Every change stays on the existing ladder + Clair3 + bcftools path:

- Genotype and consensus-selector decisions for length-partitioned alleles move into `phasing.length_partition_selection` (D1). They use the allele-specific AD fraction computed in `vcf.haploid_genotype` (D5).
- VCF concordance mirrors `bcftools consensus -H I` and is judged per event (D4).
- Clinical gates move into a new pure module, `clinical_gates.py`. It holds explicit support (hybrid Task 1), event identity, and the depth, selection, genotype and length gates (D8/#55, D3).
- Allele-selection and depth QC is computed once at pipeline time in a new `selection_qc.py` and stored as additive `alleles.json` fields.
- Three small fixes: calling settings are honoured and recorded (D7), an IGV preflight runs before mapping (D9), and the cluster and remapped BAM paths are separated (D10).

**Tech Stack:** Python 3.10+, Click, pytest (unittest.mock at the tool boundary), bcftools (integration only), existing `tools.run_tool` abstraction.

**Spec:** `.planning/2026-09-23-deep-review-roadmap.md` (section 2 defects D1–D10, section 3 P0) and `.planning/2026-09-24-execution-order.md` (Phase 0 table). Explicit-support rule: `.planning/2026-09-23-hybrid-engine-plan.md` Task 1 (moved into this plan as Task 4). Evidence (outside Git): `../MucOneSpan-review-20260923/code-review/REPORT.md`, `.../code-review/mp4_replay/`, `.../code-review/bcf_test/`. Cohort outputs: `../MucOneSpan-wave2-data/cohort-v3/`.

## Global Constraints

- Every authored code, configuration and template file has **fewer than 650 physical lines** (max 649). Split by responsibility; never compress lines or add exclusions. Watch these files:
  - `cli.py` (638): Task 6 keeps its line count unchanged.
  - `alleles.py` (641): this plan does not touch it.
  - `calling.py`: 550 now, 578 after Task 7.
- Run `make ci-check` before every commit. Also run `make test-int`, `make docs-check` and `make build-check` where a task says so. Report the pass/fail/skip counts and any missing prerequisites. A skipped test does not validate its behaviour.
- Preserve CLI flags, output schemas, bundled resources and scientific semantics. Output schemas change by **adding fields only**. User-facing version values come from `src/muc_one_span/version.py`.
- Do not tune thresholds silently. The only thresholds introduced are the settings listed below, each with its default and evidence. Default filter commands stay byte-identical: QUAL>=4.0 on length-partitioned VCFs, and 0.5 / 0.2 haploid cut-offs.
- Never put patient reads, sequences or in-house sample IDs in code, tests, fixtures, commits or PR text. PRJEB92208 public accessions (ERR152775xx) and the names MP1–MP5, HG001–HG004 and HG002 are allowed. Cohort outputs stay outside Git.
- Subprocesses run only through `tools.run_tool`/`check_tools`, with argument lists, no shell and visible failures.
- Unit tests are deterministic and mock external tools. Integration tests carry `@pytest.mark.integration`.
- One commit per issue-sized task, TDD, conventional commit messages that reference the issue, and the trailer `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Do not merge, tag or publish without explicit owner approval.

### Settings introduced (all additive, validated, documented in `docs/guides/configuration.md`)

| Setting | Default | Rationale |
| --- | --- | --- |
| `calling.haploid_alt_fraction` | `0.5` | The existing hardcoded ALT cut-off (`vcf.py:115-130`), now applied to the AD fraction and made configurable. |
| `calling.haploid_ref_fraction` | `0.2` | The existing hardcoded REF cut-off, likewise. The pair is validated so that ref < alt. |
| `calling.haploid_min_qual` | `4.0` (now `float \| null`) | The value is unchanged. `null` makes length-partitioned calls follow `--min-qual` (D7). |
| `allele_selection.secondary_mode_min_fraction` | `0.2` | Measured on cohort-v3 as the largest primary-record count at least `min_gap` (5) units from the cluster peak, divided by the peak count. Controls HG001–HG004: ≤0.045. MP5 fragment mode: 0.276. MP3 merged allele: 0.619. Spurious 133-unit MP3 candidate: 1.0. 0/80 simulated review-panel runs reach 0.2. This matches the roadmap's 20% minor-peak figure. |
| `allele_selection.min_allele_primary_records` | `30` | Roadmap P2.6 per-allele callable gate of about 30 spanning reads (MP4 long-allele dupC was lost at about 20). This replaces the report's old 30-record total, which counted alignment records including secondaries. Primary records serve as the molecule proxy, and the basis is recorded as `depth_basis`. |

The two allele-selection gates can only turn NEGATIVE or PATHOGENIC into INCONCLUSIVE. They never create a call.

## Review Focus

1. **MP4 shape** (`contig_71:3432 G>GC`, `AD 192,275`, `AF 0.4479`, 80-unit allele): Task 2 makes it `1/1`. Without Task 2, Task 1 must give selector `I` and INCONCLUSIVE, never REF with NEGATIVE.
2. **MP3 shape** (`AD 8780,5105`, AD fraction 0.368, in the ambiguous band): the dupC stays `0/1`. Under `I` bcftools still applies it, and support must read `heterozygous_genotype_unresolved`. That is INCONCLUSIVE, not PATHOGENIC and not NEGATIVE.
3. **One unrelated heterozygous indel on an allele** must no longer remove support from a homozygous frameshift on the same allele (Task 3, per event). A multi-ALT heterozygous indel still reports `ambiguous_genotype_selection`.
4. **MP5 shape** (a supported frameshift without a dictionary template, plus a 0.276 fragment mode): becomes INCONCLUSIVE with an explicit reason. It was PATHOGENIC in cohort-v3. The owner must accept this in the PR (Task 8).
5. **Legacy summaries** without the new fields keep their established behaviour, except where the rules below require evidence:
   - a missing `vcf_support` key is no longer support;
   - PATHOGENIC now needs `template_match` and a `mutation_name`;
   - a mismatch between `length` and `reference_length` blocks NEGATIVE.

---

## Defect verification (current code, `main` @ 10a1440)

| ID | Verified | Where (current lines) | Planned in |
| --- | --- | --- | --- |
| D1 | **Real.** `single_heterozygous_unordered` is missing from the unphased tuple, so the selector is GT1 and `independent=True`. Cohort-v3 MP4 `alleles.allele_2` shows exactly this, and the decision is NO_PATHOGENIC. `test_distinct_calling.py:157-160` asserts the bug. | `calling.py:525-540` | Task 1 (#53) |
| D5 | **Real.** The rule uses `FORMAT/AF`. Replaying cohort-v3 with the AD fraction changes 9 records: the MP4 dupC, 7 MP3 records, and 2 records that are already `1/1`. No control changes. | `vcf.py:103-143` (rule 115-130) | Task 2 |
| D4 | **Real.** bcftools 1.17 `-H I` applies `C>CC 0/1` and `CC>C 0/1` as ALT and turns `G>T 0/1` into `K`. `_selected_alt` returns `None`, so the whole allele becomes unavailable (MP3 allele 1). | `variant_support.py:36-51, 54-85, 123-150` | Task 3 |
| Hybrid T1 | **Real.** A missing `vcf_support` key counts as supported. | `report.py:93-96` | Task 4 (#55) |
| D8 | **Real.** PATHOGENIC is decided before coverage and phase are checked, needs no event identity, and never consults allele selection. Coverage counts alignment records, and 83% of those are secondary. | `report.py:75-221` | Task 5 (#55) |
| D3 | **Real.** Cohort-v3 MP1 has 39 vs 44 units, MP3 37 vs 44 and 133 vs 117, MP5 69 vs 61. `reference_length` exists but no status is set and no gate uses it. P0 makes the mismatch explicit and blocks NEGATIVE; choosing the right contig stays with #20/P2. | `alleles.py:391-414` | Task 5 |
| D7 | **Real, narrowed.** `--min-qual` is overridden by `calling.haploid_min_qual` without any message. `vcf.py:86-87` hides a second 4.0 override. `calling.haploid_majority` is never read. `consensus.haploid_*` are dead settings. **Dropped:** `min_dp`, which is documented as accepted-but-not-applied API compatibility. | `calling.py:307, 515-524`; `vcf.py:83-88`; `settings.py:146-147` | Task 6 |
| D9 | **Partly real.** The igv-reports 1.16 rejection is intentional and documented (`docs/development.md`, report lifecycle). The defect is its timing: the check runs only after the whole analysis, and `create_report` is not in the `check_tools` list. A real probe passes with 1.13 and fails with the local 1.16.1. | `pipeline.py:90, 198-214`; `report_igv.py:100-133` | Task 7 |
| D10 | **Real.** Low severity: the path is overwritten while its stale index remains. | `calling.py:52, 102, 139-151` | Task 7 |

D2 (super-cluster splitting) and D6 (HiFi QUAL/cap) are out of P0 scope. Task 5 gates the effect of D2 without fixing it.

All tasks were dry-run in a scratch copy of this worktree. That copy produced 1052 passing unit tests, clean ruff, format and mypy, and passed the file-size gate. The integration suite also passed, with two failures that are pre-existing environment conditions: IGV session tests need igv-reports 1.13 first on `PATH`, and the Clair3 run test fails on baseline too because no model is installed on this machine.

## File Structure

| Path | Change | Responsibility |
| --- | --- | --- |
| `src/muc_one_span/phasing.py` | modify | `length_partition_selection`: selector and genotype status for length-partitioned alleles (Task 1) |
| `src/muc_one_span/calling.py` | modify | Use the selector (T1). Forward the fraction settings (T2). `_haploid_filter` settings and provenance (T6). Separate cluster BAM (T7). |
| `src/muc_one_span/vcf.py` | modify | `haploid_genotype` from the AD fraction (T2). Explicit QUAL rule (T6). |
| `src/muc_one_span/variant_support.py` | modify | Replay mirrors `-H I` for heterozygous REF/ALT indels; support judged per event (T3). |
| `src/muc_one_span/templates/report.html.j2` | modify | "Genotype unresolved" support badge (T3) |
| `src/muc_one_span/clinical_gates.py` | **create** | `mutation_supported` (T4). `mutation_blockers`, `allele_gate_reasons`, `LEGACY_MIN_TOTAL_READS` (T5). |
| `src/muc_one_span/report.py` | modify | Gate order in `compute_clinical_decision` (T4, T5) |
| `src/muc_one_span/selection_qc.py` | **create** | `secondary_mode_fraction`, `assess_allele`, `annotate_selection_qc` (T5) |
| `src/muc_one_span/settings.py` | modify | New fields (T2, T5, T6). Deprecation warning for `consensus.haploid_*` (T6). |
| `src/muc_one_span/pipeline.py` | modify | Call `annotate_selection_qc` (T5). IGV tool check and preflight (T7). |
| `src/muc_one_span/report_igv.py` | modify | `preflight_igv_report` (T7) |
| `src/muc_one_span/cli.py` | modify | `--min-qual` help text only; line count unchanged (T6) |
| `tests/unit/test_distinct_calling.py` | modify | Replace the D1-asserting test; add the MP4 shape and the selector matrix (T1, T2, T6) |
| `tests/unit/test_vcf.py` | modify | AD-fraction and QUAL-rule tests (T2, T6) |
| `tests/unit/test_variant_support.py` | modify | Per-event tests (T3) |
| `tests/integration/test_variant_concordance.py` | modify | Real-bcftools `-H I` guard (T3) |
| `tests/unit/test_clinical_decision.py` | modify | Explicit support, gates, MP-shaped decisions (T1, T4, T5) |
| `tests/unit/test_selection_qc.py` | **create** | Selection, depth and length gates (T5) |
| `tests/unit/test_pipeline_gates.py` | **create** | Mocked pipeline: QC wiring (T5), IGV preflight ordering (T7) |
| `tests/unit/test_runtime_settings.py`, `test_report.py`, `test_report_wave1.py`, `test_report_igv.py`, `test_calling.py` | modify | As listed per task |
| `docs/development.md`, `docs/guides/configuration.md` | modify | Contracts and settings (T1–T6) |
| `CHANGELOG.md`, `src/muc_one_span/version.py` | modify | 0.16.0 (T8) |
| `.planning/2026-09-23-hybrid-engine-plan.md` | modify | Mark Task 1 as done in P0 (T8) |

---

### Task 0: File the issues and post the #53 evidence

**Files:** none in the repository. The drafts live outside it.

- [ ] **Step 1: Create the issues from the drafts** (owner approval required)

The drafts are in the session scratchpad, `p0-issues/D3.md`, `D4.md`, `D5.md`, `D7.md`, `D9.md` and `D10.md`. The first line of each is `TITLE: ...`.

```bash
DRAFTS=/path/to/p0-issues   # the scratchpad directory holding the drafts
for id in D5 D4 D3 D7 D9 D10; do
  title=$(head -1 "$DRAFTS/$id.md" | sed 's/^TITLE: //')
  tail -n +3 "$DRAFTS/$id.md" > "/tmp/$id-body.md"
  gh issue create --title "$title" --body-file "/tmp/$id-body.md"
done
gh issue comment 53 --body-file "$DRAFTS/comment-53.md"
```

- [ ] **Step 2: Record the issue numbers**

In every commit command below, replace each `#65`, `#64`, `#63`, `#66`, `#67` and `#68` with the number that `gh` printed.

---

### Task 1: D1 — never resolve a heterozygous length-partitioned allele to REF (#53)

**Files:**
- Modify: `src/muc_one_span/phasing.py` (append after line 73)
- Modify: `src/muc_one_span/calling.py:11` (import), `calling.py:525-540` (selector block in `_process_allele`)
- Modify: `tests/unit/test_distinct_calling.py` (imports at lines 1-8; delete lines 98-161; append)
- Modify: `tests/unit/test_clinical_decision.py` (append)
- Modify: `docs/development.md:231-235`

**Interfaces:**
- Consumes: `phase_evidence(variants) -> dict` (unchanged) and `annotate_consensus_candidate(allele_info, evidence, haplotype, sample, vcf_path)` (unchanged).
- Produces: `phasing.length_partition_selection(variants: list[dict], evidence: dict) -> tuple[int | str, str, list[dict]]`, returning `(selector, allele_genotype_status, heterozygous_sites)`. The status is one of `allele_specific_resolved`, `heterozygous_within_length_partition` or `unresolved_genotype_records`.
- New additive `alleles.json` fields per distinct-length allele: `allele_genotype_status` (str) and `heterozygous_sites` (a list of `{chrom,pos,ref,alt,genotype}`). Task 5 reads `allele_genotype_status`.

**Decision:** the reads were already partitioned by length, so each candidate is one haplotype. After Task 2's allele-fraction rule, a heterozygous record can remain only in the ambiguous band, or when diploid genotypes were kept by configuration. Such a record is unresolved allele-specific evidence. The rule does **not** pick allele 2. Instead it applies selector `I`, withholds independent haplotype credit, and records an explicit status:

| Case | Selector | Status | Credit |
| --- | --- | --- | --- |
| Homozygous ALT, or no records | `1` | resolved | independent credit |
| Single-site heterozygous | `I` | het | no credit |
| Phased heterozygous | `I` | het | no credit |
| Multi-site unphased | `I` | het | no credit |
| Conflicting or incomplete records | `I` | unresolved | no credit |

Phased heterozygosity inside one length partition means two sequences share one allele slot. Choosing GT1 there would repeat the D1 failure.

- [ ] **Step 1: Replace the test that asserts the bug and add the regression tests**

In `tests/unit/test_distinct_calling.py`, delete `test_distinct_length_phased_and_single_het`, which is the whole function from its `@patch` decorators at line 98 to line 161. Replace the import block (lines 3-7) with:

```python
from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import patch

import pytest

from muc_one_span.calling import call_variants_per_allele
```

Append:

```python
def _variant(
    pos: int, ref: str, alt: str, genotype: str, phase_set: str | None = None
) -> dict[str, Any]:
    return {
        "chrom": "contig_71",
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "qual": 29.93,
        "genotype": genotype,
        "phase_set": phase_set,
        "sample": "S1",
    }


def _run_distinct(
    tmp_path: Path, allele_2_variants: list[dict], **call_kwargs: Any
) -> tuple[dict, list[dict[str, Any]]]:
    """Run distinct-length calling with tools mocked; allele 2 is MP4-shaped (80 units)."""
    alleles = _make_alleles()
    alleles["allele_2"].update(
        length=80, canonical_repeats=71, contig_name="contig_71", cluster_contigs=["contig_71"]
    )
    by_allele: dict[str, list[dict]] = {"allele_1": [], "allele_2": allele_2_variants}
    filter_calls: list[dict[str, Any]] = []

    def fake_filter(vcf: Path, reference: Path, out_dir: Path, **kwargs: Any) -> Path:
        filter_calls.append(kwargs)
        return out_dir / "variants.vcf.gz"

    with (
        patch("muc_one_span.calling._extract_and_remap_reads", return_value=tmp_path / "m.bam"),
        patch("muc_one_span.calling.run_clair3", return_value=tmp_path / "raw.vcf.gz"),
        patch("muc_one_span.calling.filter_vcf", side_effect=fake_filter),
        patch(
            "muc_one_span.calling.parse_vcf_genotypes",
            side_effect=lambda path, sample=None: by_allele[path.parent.name],
        ),
    ):
        call_variants_per_allele(
            tmp_path / "in.bam", tmp_path / "ref.fa", alleles, tmp_path / "out", **call_kwargs
        )
    return alleles, filter_calls


def test_mp4_shaped_single_het_insertion_is_not_resolved_as_reference(tmp_path: Path) -> None:
    """ERR15277569 shape: a 0/1 G>GC left after the allele-fraction rule must not become REF."""
    alleles, _ = _run_distinct(tmp_path, [_variant(3432, "G", "GC", "0/1")])
    a2 = alleles["allele_2"]
    assert a2["phase_status"] == "single_heterozygous_unordered"
    assert a2["consensus_haplotype"] == "I"
    assert a2["consensus_policy"] == "genotype_iupac_candidate"
    assert a2["allele_genotype_status"] == "heterozygous_within_length_partition"
    assert a2["heterozygous_sites"] == [
        {"chrom": "contig_71", "pos": 3432, "ref": "G", "alt": "GC", "genotype": "0/1"}
    ]
    assert a2["independent_haplotype_evidence"] is False
    assert alleles["allele_1"]["consensus_haplotype"] == 1
    assert alleles["allele_1"]["independent_haplotype_evidence"] is True


@pytest.mark.parametrize(
    ("variants", "selector", "status", "independent"),
    [
        pytest.param(
            [_variant(3432, "G", "GC", "1/1")],
            1,
            "allele_specific_resolved",
            True,
            id="homozygous_alt",
        ),
        pytest.param([], 1, "allele_specific_resolved", True, id="no_retained_variants"),
        pytest.param(
            [_variant(3432, "G", "GC", "0/1")],
            "I",
            "heterozygous_within_length_partition",
            False,
            id="single_site_het",
        ),
        pytest.param(
            [_variant(100, "G", "GC", "0|1", "100"), _variant(200, "C", "A", "1|0", "100")],
            "I",
            "heterozygous_within_length_partition",
            False,
            id="phased",
        ),
        pytest.param(
            [_variant(100, "G", "GC", "0/1"), _variant(200, "C", "A", "0/1")],
            "I",
            "heterozygous_within_length_partition",
            False,
            id="multi_site_unphased",
        ),
        pytest.param(
            [_variant(100, "GC", "G", "1/1"), _variant(101, "C", "A", "1/1")],
            "I",
            "unresolved_genotype_records",
            False,
            id="conflicting_records",
        ),
    ],
)
def test_length_partition_selector_follows_allele_specific_genotypes(
    tmp_path: Path, variants: list[dict], selector: int | str, status: str, independent: bool
) -> None:
    alleles, _ = _run_distinct(tmp_path, variants)
    a2 = alleles["allele_2"]
    assert a2["consensus_haplotype"] == selector
    assert a2["allele_genotype_status"] == status
    assert a2["independent_haplotype_evidence"] is independent
```

Append to `tests/unit/test_clinical_decision.py`. This test already passes, because a false `independent_haplotype_evidence` blocks NEGATIVE; it pins that contract:

```python
def test_heterozygous_length_partition_is_never_negative():
    """A length-partitioned allele with an unresolved heterozygous call is not NEGATIVE."""
    summary = _resolved_diploid_summary()
    summary["alleles"]["allele_2"].update(
        phase_status="single_heterozygous_unordered",
        allele_genotype_status="heterozygous_within_length_partition",
        consensus_haplotype="I",
        independent_haplotype_evidence=False,
    )
    assert compute_clinical_decision(summary)["state"] == "INCONCLUSIVE"
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run --locked --all-extras pytest tests/unit/test_distinct_calling.py --no-cov -q`
Expected: 7 tests fail and 1 passes.
- The MP4 test and the `single_site_het`/`phased` cases fail with `assert 1 == 'I'`.
- The remaining parametrized cases fail with `KeyError: 'allele_genotype_status'`.
- `test_distinct_length_unphased_retains_iupac_selector` passes.

- [ ] **Step 3: Implement**

Append to `src/muc_one_span/phasing.py`. The module already imports `re`.

```python
_UNRESOLVED_RECORD_STATUSES = frozenset(
    {"conflicting_variant_records", "missing_genotype", "non_diploid"}
)


def length_partition_selection(
    variants: list[dict], evidence: dict
) -> tuple[int | str, str, list[dict]]:
    """Choose the consensus selector for reads already partitioned by allele length.

    A length partition holds one haplotype, so genotype indices carry no phase
    meaning there. ``filter_vcf`` has already turned clear allele-specific ALT
    support into 1/1 and clear REF support into 0/0; a heterozygous record that
    remains is in the ambiguous allele-fraction band (or diploid genotypes were
    kept by configuration). Selecting GT1 would silently keep REF, so any
    remaining heterozygous or unresolved record selects the IUPAC/unphased
    candidate ``"I"``. Returns ``(selector, allele_genotype_status, sites)``.
    """
    heterozygous = [
        {key: variant[key] for key in ("chrom", "pos", "ref", "alt", "genotype")}
        for variant in variants
        if len(set(re.split(r"[/|]", variant["genotype"]))) > 1
    ]
    if evidence["phase_status"] in _UNRESOLVED_RECORD_STATUSES:
        status = "unresolved_genotype_records"
    elif heterozygous:
        status = "heterozygous_within_length_partition"
    else:
        status = "allele_specific_resolved"
    return (1 if status == "allele_specific_resolved" else "I"), status, heterozygous
```

In `src/muc_one_span/calling.py`, replace line 11 (`from muc_one_span.phasing import annotate_consensus_candidate, phase_evidence`) with:

```python
from muc_one_span.phasing import (
    annotate_consensus_candidate,
    length_partition_selection,
    phase_evidence,
)
```

Replace the block from `variants = parse_vcf_genotypes(filtered)` through `return allele_key, filtered` inside `_process_allele` (currently lines 525-540) with:

```python
        variants = parse_vcf_genotypes(filtered)
        evidence = phase_evidence(variants)
        sample = variants[0].get("sample") if variants else None
        haplotype, genotype_status, heterozygous = length_partition_selection(variants, evidence)
        annotate_consensus_candidate(allele_info, evidence, haplotype, sample, str(filtered))
        allele_info["allele_genotype_status"] = genotype_status
        allele_info["heterozygous_sites"] = heterozygous
        allele_info["independent_haplotype_evidence"] = len(allele_keys) > 1 and haplotype == 1
        return allele_key, filtered
```

In `docs/development.md`, after the sentence ending "bcftools selectors alone do not establish phase." (line 235), insert:

```markdown
Distinct-length candidates are length-partitioned haplotypes. After the
allele-fraction rule, any remaining heterozygous record (single-site, phased or
multi-site) or any conflicting record selects `I`. The allele records
`allele_genotype_status` (`heterozygous_within_length_partition` or
`unresolved_genotype_records`) and `heterozygous_sites`, and it receives no
independent haplotype credit. Only candidates without remaining heterozygosity
use genotype index 1 (`allele_specific_resolved`).
```

- [ ] **Step 4: Run the tests again**

Run: `uv run --locked --all-extras pytest tests/unit/test_distinct_calling.py tests/unit/test_clinical_decision.py tests/unit/test_calling.py tests/unit/test_phasing.py --no-cov -q`
Expected: all pass. Then run `make ci-check`, which must pass.

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/phasing.py src/muc_one_span/calling.py docs/development.md \
  tests/unit/test_distinct_calling.py tests/unit/test_clinical_decision.py
git commit -m "fix(calling): keep heterozygous length-partitioned alleles unresolved instead of REF

A heterozygous genotype left on a length-partitioned (distinct-length) allele
selected GT1, i.e. REF for 0/1, with independent haplotype credit (MP4 dupC,
ERR15277569). Any remaining heterozygous, phased or conflicting record now
selects the IUPAC candidate with an explicit allele_genotype_status and no
independent credit. The removed test asserted the defect.

Fixes #53

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: D5 — haploid genotypes from the allele-specific AD fraction

**Files:**
- Modify: `src/muc_one_span/vcf.py:12-150` (`filter_vcf` signature, docstring, haploid block 103-143) and add `haploid_genotype` before `select_vcf_sample` (line 153)
- Modify: `src/muc_one_span/settings.py:201-218` (`CallingSettings`)
- Modify: `src/muc_one_span/calling.py`: the two `filter_vcf(..., haploid_majority=True, ...)` calls (haplotag path, about line 350; distinct path, about line 563 after Task 1)
- Modify: `tests/unit/test_vcf.py`, `tests/unit/test_runtime_settings.py`, `tests/unit/test_distinct_calling.py` (append)
- Modify: `docs/guides/configuration.md:124-129` (calling table)

**Interfaces:**
- Produces: `vcf.haploid_genotype(fmt: list[str], sample_fields: list[str], *, alt_fraction: float = 0.5, ref_fraction: float = 0.2) -> str | None`.
- Produces: new `filter_vcf` keywords `haploid_alt_fraction: float = 0.5` and `haploid_ref_fraction: float = 0.2`.
- Produces: `CallingSettings.haploid_alt_fraction: float = 0.5` and `CallingSettings.haploid_ref_fraction: float = 0.2`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_vcf.py`:

```python
# Public PRJEB92208 record shapes (cohort-v3 Clair3 output, ladder coordinates).
MP4_DUPC = "contig_71\t3432\t.\tG\tGC\t29.93\tPASS\tF\tGT:GQ:DP:AD:AF\t0/1:29:614:192,275:0.4479"
MP3_DUPC = "contig_35\t912\t.\tG\tGC\t21.05\tPASS\tF\tGT:GQ:DP:AD:AF\t0/1:21:18273:8780,5105:0.2794"


def _haploid_rewrite(tmp_path: Path, records: list[str], **kwargs: float) -> list[str]:
    """Run filter_vcf's haploid rewrite on literal records and return the data lines."""
    vcf = tmp_path / "input.vcf.gz"
    vcf.touch()
    ref = tmp_path / "ref.fa"
    ref.touch()
    out_dir = tmp_path / "out"
    content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        + "".join(record + "\n" for record in records)
    )
    written: list[str] = []

    def side_effect(cmd: list[str]) -> str:
        if cmd[1] == "norm":
            (out_dir / "normalized.vcf.gz").write_bytes(b"data")
        elif cmd[1] == "view" and str(cmd[-1]).endswith("normalized.vcf.gz"):
            (out_dir / "variants.vcf.gz").write_bytes(b"data")
        elif cmd[1] == "view" and str(cmd[-1]).endswith("variants.vcf.gz"):
            return content
        elif cmd[1] == "view" and str(cmd[-1]).endswith("mod_haploid.vcf"):
            written.extend(Path(cmd[-1]).read_text().splitlines())
            (out_dir / "variants.vcf.gz").touch()
        return ""

    with patch("muc_one_span.vcf.run_tool", side_effect=side_effect):
        filter_vcf(vcf, ref, out_dir, min_qual=5.0, haploid_majority=True, **kwargs)
    data = [line for line in written if not line.startswith("#")]
    return data or list(records)


def _genotype(line: str) -> str:
    return line.split("\t")[9].split(":")[0]


def test_haploid_rule_uses_allele_specific_depth_fraction(tmp_path: Path) -> None:
    """MP4 dupC: AF is 0.448 (ALT/DP) but AD gives 275/467 = 0.589, so ALT."""
    assert [_genotype(line) for line in _haploid_rewrite(tmp_path, [MP4_DUPC])] == ["1/1"]


def test_ambiguous_ad_fraction_keeps_heterozygous_genotype(tmp_path: Path) -> None:
    """MP3 dupC: 5105/13885 = 0.368 lies in the ambiguous band and stays 0/1."""
    assert [_genotype(line) for line in _haploid_rewrite(tmp_path, [MP3_DUPC])] == ["0/1"]


def test_low_ad_fraction_becomes_reference(tmp_path: Path) -> None:
    record = "contig_1\t10\t.\tC\tCC\t12.0\tPASS\t.\tGT:AD:AF\t0/1:90,10:0.1"
    assert [_genotype(line) for line in _haploid_rewrite(tmp_path, [record])] == ["0/0"]


def test_multiallelic_ad_selects_supported_alt(tmp_path: Path) -> None:
    record = "contig_1\t10\t.\tC\tCC,CCC\t20.0\tPASS\t.\tGT:AD:AF\t1/2:10,20,70:0.2,0.7"
    assert [_genotype(line) for line in _haploid_rewrite(tmp_path, [record])] == ["2/2"]


def test_af_is_fallback_when_ad_is_missing_or_zero(tmp_path: Path) -> None:
    records = [
        "contig_1\t10\t.\tC\tG\t15.0\tPASS\t.\tGT:AD:AF\t0/1:0,0:0.8",
        "contig_1\t20\t.\tA\tT\t15.0\tPASS\t.\tGT:AD:AF\t0/1:.:0.1",
    ]
    assert [_genotype(line) for line in _haploid_rewrite(tmp_path, records)] == ["1/1", "0/0"]


def test_haploid_fraction_cutoffs_are_configurable(tmp_path: Path) -> None:
    lines = _haploid_rewrite(
        tmp_path, [MP3_DUPC], haploid_alt_fraction=0.35, haploid_ref_fraction=0.1
    )
    assert [_genotype(line) for line in lines] == ["1/1"]
```

Append to `tests/unit/test_runtime_settings.py`:

```python
def test_haploid_fraction_settings_are_validated() -> None:
    calling = CallingSettings()
    assert (calling.haploid_alt_fraction, calling.haploid_ref_fraction) == (0.5, 0.2)
    with pytest.raises(ValueError, match="haploid_ref_fraction must be below"):
        CallingSettings(haploid_alt_fraction=0.3, haploid_ref_fraction=0.3)
    with pytest.raises(ValueError, match="finite number"):
        CallingSettings(haploid_alt_fraction=1.5)
```

In `tests/unit/test_distinct_calling.py`, add `from muc_one_span.settings import CallingSettings` below the `calling` import, then append:

```python
def test_distinct_path_forwards_allele_fraction_settings(tmp_path: Path) -> None:
    settings = CallingSettings(haploid_alt_fraction=0.6, haploid_ref_fraction=0.1)
    _, filter_calls = _run_distinct(tmp_path, [], settings=settings)
    assert len(filter_calls) == 2
    for kwargs in filter_calls:
        assert (kwargs["haploid_alt_fraction"], kwargs["haploid_ref_fraction"]) == (0.6, 0.1)
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run --locked --all-extras pytest tests/unit/test_vcf.py tests/unit/test_runtime_settings.py tests/unit/test_distinct_calling.py --no-cov -q`
Expected: 4 tests fail.
- `test_haploid_rule_uses_allele_specific_depth_fraction` fails with `['0/1'] != ['1/1']`.
- `test_haploid_fraction_cutoffs_are_configurable` fails with a `TypeError` for an unexpected keyword.
- `test_haploid_fraction_settings_are_validated` fails with a `TypeError`.
- `test_distinct_path_forwards_allele_fraction_settings` fails with a `KeyError`.

The existing AF-only test `test_filter_vcf_haploid_majority` still passes.

- [ ] **Step 3: Implement**

`src/muc_one_span/settings.py`, class `CallingSettings`: add two fields after `haploid_min_qual: float = 4.0`, and add these lines after `_number("calling.haploid_min_qual", ...)`:

```python
    haploid_alt_fraction: float = 0.5
    haploid_ref_fraction: float = 0.2
```

```python
        _number("calling.haploid_alt_fraction", self.haploid_alt_fraction, 0.0, 1.0)
        _number("calling.haploid_ref_fraction", self.haploid_ref_fraction, 0.0, 1.0)
        if self.haploid_ref_fraction >= self.haploid_alt_fraction:
            raise ValueError(
                "calling.haploid_ref_fraction must be below calling.haploid_alt_fraction"
            )
```

`src/muc_one_span/vcf.py`:

- Signature: after `haploid_min_qual: float | None = None,` add `haploid_alt_fraction: float = 0.5,` and `haploid_ref_fraction: float = 0.2,`.
- Docstring: replace the sentence "When ``haploid_majority=True``, resolves borderline heterozygous calls on isolated haploid alignments by setting GT to 1/1 when AF >= 0.5 and 0/0 otherwise, avoiding spurious IUPAC ambiguity characters in consensus." with:

```text
    the result. When ``haploid_majority=True``, genotypes on length-partitioned
    (haploid) alignments follow :func:`haploid_genotype`: ALT when the
    allele-specific AD fraction is >= ``haploid_alt_fraction``, REF when it is
    below ``haploid_ref_fraction``, otherwise the original genotype is kept.
```

- Args: replace the `haploid_majority:` line with:

```text
        haploid_majority: If True, resolve genotypes on haploid alignments from AD fractions.
        haploid_alt_fraction: Allele-specific ALT fraction at or above which GT becomes ALT.
        haploid_ref_fraction: ALT fraction below which GT becomes 0/0.
```

- Replace the body of `if haploid_majority and not is_empty and filtered.exists():` up to, but not including, `if modified:` with:

```python
    if haploid_majority and not is_empty and filtered.exists():
        lines = run_tool(["bcftools", "view", str(filtered)]).splitlines()
        new_lines: list[str] = []
        modified = False
        for line in lines:
            fields = line.split("\t")
            if not line.startswith("#") and len(fields) >= 10:
                fmt = fields[8].split(":")
                sample_fields = fields[9].split(":")
                target_gt = haploid_genotype(
                    fmt,
                    sample_fields,
                    alt_fraction=haploid_alt_fraction,
                    ref_fraction=haploid_ref_fraction,
                )
                if "GT" in fmt and target_gt is not None:
                    gt_idx = fmt.index("GT")
                    if sample_fields[gt_idx] != target_gt:
                        sample_fields[gt_idx] = target_gt
                        fields[9] = ":".join(sample_fields)
                        modified = True
            new_lines.append("\t".join(fields))
```

- Insert before `def select_vcf_sample(`:

```python
def haploid_genotype(
    fmt: list[str],
    sample_fields: list[str],
    *,
    alt_fraction: float = 0.5,
    ref_fraction: float = 0.2,
) -> str | None:
    """Return the haploid genotype implied by allele-specific read support.

    Reads were partitioned by allele length, so each site asks REF versus ALT among
    reads supporting either allele: ``ALT_i / sum(AD)``. Reads supporting neither
    (for example other homopolymer lengths) count in DP but not in AD, so FORMAT/AF
    (ALT/DP) undercounts and is used only when AD is absent, malformed or zero.
    Returns None for the ambiguous band [ref_fraction, alt_fraction) or no evidence.
    """
    fractions: list[float] = []
    if "AD" in fmt:
        try:
            depths = [int(value) for value in sample_fields[fmt.index("AD")].split(",")]
        except (ValueError, IndexError):
            depths = []
        if len(depths) > 1 and sum(depths) > 0:
            fractions = [depth / sum(depths) for depth in depths[1:]]
    if not fractions and "AF" in fmt:
        try:
            fractions = [
                float(value) for value in sample_fields[fmt.index("AF")].split(",") if value != "."
            ]
        except (ValueError, IndexError):
            fractions = []
    if not fractions:
        return None
    best = max(fractions)
    if best >= alt_fraction:
        index = fractions.index(best) + 1
        return f"{index}/{index}"
    return "0/0" if best < ref_fraction else None
```

`src/muc_one_span/calling.py`: in both `filter_vcf(...)` calls that pass `haploid_majority=True`, add these two keyword arguments after the existing ones:

```python
                    haploid_alt_fraction=settings.haploid_alt_fraction,
                    haploid_ref_fraction=settings.haploid_ref_fraction,
```

In the distinct path the indentation is one level less. In both functions `settings` is the `CallingSettings` already in scope.

`docs/guides/configuration.md`: after the `calling.read_phase` row, add:

```markdown
| `calling.haploid_alt_fraction` | `0.5` | Length-partitioned calls: the ALT fraction of allele-specific reads (`FORMAT/AD`) at or above which the genotype becomes ALT. Number in [0,1]. |
| `calling.haploid_ref_fraction` | `0.2` | ALT fraction below which the genotype becomes REF (`0/0`). Values in between keep the heterozygous call, which stays unresolved. Must be below `haploid_alt_fraction`. |
```

- [ ] **Step 4: Run the tests again**

Run: `uv run --locked --all-extras pytest tests/unit/test_vcf.py tests/unit/test_runtime_settings.py tests/unit/test_distinct_calling.py tests/unit/test_calling.py --no-cov -q`
Expected: all pass. Then `make ci-check` must pass.

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/vcf.py src/muc_one_span/settings.py src/muc_one_span/calling.py \
  docs/guides/configuration.md tests/unit/test_vcf.py tests/unit/test_runtime_settings.py \
  tests/unit/test_distinct_calling.py
git commit -m "fix(vcf): decide haploid genotypes from allele-specific AD fractions

FORMAT/AF (ALT/DP) undercounts ALT on length-partitioned alleles because DP
includes reads supporting neither allele (MP4 dupC: AF 0.448 vs AD 0.589).
Use ALT/(REF+ALT) from AD with an AF fallback, and expose the existing 0.5/0.2
cut-offs as calling.haploid_alt_fraction / haploid_ref_fraction.

Fixes #65
Refs #59

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: D4 — per-event VCF support under the IUPAC consensus

**Files:**
- Modify: `src/muc_one_span/variant_support.py:36-51` (`_selected_alt`), `:81` (applied edit), `:123` (projection), `:145` (per-event status)
- Modify: `src/muc_one_span/templates/report.html.j2:265-266` (badge)
- Modify: `tests/unit/test_variant_support.py:109-134` (parametrize) and append
- Modify: `tests/unit/test_report.py:194-206`
- Modify: `tests/integration/test_variant_concordance.py` (refactor into a helper and add the heterozygous case)
- Modify: `docs/development.md:242`

**Interfaces:**
- Produces: a new `vcf_support_status` value, `heterozygous_genotype_unresolved`, with `vcf_support=False`. Task 5 reads it.
- Produces: the additive field `classification.vcf_projection.unresolved_genotype_edits: int`.
- Every applied replay edit carries `unresolved_genotype: bool`. This is internal.

- [ ] **Step 1: Write the failing tests**

In `tests/unit/test_variant_support.py`, replace the parametrize entry `("heterozygous_indel", "ambiguous_genotype_selection"),` with `("multiallelic_heterozygous_indel", "ambiguous_genotype_selection"),`. Then replace the matching branch:

```python
    if problem == "multiallelic_heterozygous_indel":
        anchor = variants[0]["ref"]
        variants[0].update(alt=f"{anchor}TT,{anchor}T", genotype="1/2")
```

Append:

```python
def _validate(rd, seq: str, variants: list[dict], context: dict) -> dict:
    return validate_mutations_against_vcf(
        classify_sequence(seq, rd),
        variants,
        sequence=seq,
        repeat_dict=rd,
        consensus_context=context,
    )


def test_unrelated_heterozygous_indel_keeps_per_event_support(tmp_path: Path) -> None:
    """bcftools -H I applies a 0/1 flank insertion; the homozygous VNTR event stays supported."""
    rd, seq, variants, context = fixture(tmp_path)
    reference = Path(context["reference_path"]).read_text().splitlines()[1]
    flank = len(reference) - 4  # right flank "TGCA"; insert A after its G
    assert reference[flank + 1] == "G"
    variants.append(
        {"chrom": "contig_1", "pos": flank + 2, "ref": "G", "alt": "GA", "genotype": "0/1"}
    )
    consensus = Path(context["full_consensus_path"])
    full = consensus.read_text().splitlines()[1]
    consensus.write_text(">contig_1\n" + full[:-2] + "A" + full[-2:] + "\n")
    context["trim_end"] = len(full) + 1 - 5
    result = _validate(rd, seq, variants, context)
    mutation = result["mutations_detected"][0]
    assert result["vcf_projection"]["status"] == "available"
    assert result["vcf_projection"]["unresolved_genotype_edits"] == 1
    assert mutation["vcf_support_status"] == "exact_sequence_concordance"
    assert mutation["vcf_support"] is True


def test_heterozygous_event_under_iupac_is_unresolved_not_supported(tmp_path: Path) -> None:
    rd, seq, variants, context = fixture(tmp_path)
    variants[0]["genotype"] = "0/1"
    result = _validate(rd, seq, variants, context)
    mutation = result["mutations_detected"][0]
    assert result["vcf_projection"]["status"] == "available"
    assert mutation["vcf_support_status"] == "heterozygous_genotype_unresolved"
    assert mutation["vcf_support"] is False


def test_phase_selected_heterozygous_event_remains_supported(tmp_path: Path) -> None:
    rd, seq, variants, context = fixture(tmp_path)
    variants[0]["genotype"] = "0|1"
    context["haplotype"] = 2
    mutation = _validate(rd, seq, variants, context)["mutations_detected"][0]
    assert mutation["vcf_support_status"] == "exact_sequence_concordance"
```

In `tests/unit/test_report.py`, in `test_report_distinguishes_unavailable_from_absent_support`, add after the `ambiguous = ...` line:

```python
    unresolved = {
        **mutation,
        "repeat_index": 11,
        "vcf_support_status": "heterozygous_genotype_unresolved",
    }
```

Then append it with `sample_summary["classifications"]["allele_1"]["mutations"].append(unresolved)`, and add the assertion `assert html.count(">Genotype unresolved<") == 1`.

Replace the body of `tests/integration/test_variant_concordance.py` below `pytestmark`:

```python
def _named_mutation_case(name: str, genotype: str, tmp_path: Path) -> tuple[str, str, dict]:
    """Build one named mutation VCF, run real filter/consensus/trim, return sequence and result."""
    rd = load_repeat_dictionary()
    mutant, (parent, _) = next(
        (s, label) for s, label in rd.mutated_sequences.items() if label[1] == name
    )
    wild = rd.repeats[parent]
    prefix, suffix = "ACGT" + rd.repeats["X"], rd.repeats["X"] + "TGCA"
    ref, alt = prefix + wild + suffix, prefix + mutant + suffix
    left = 0
    while ref[left] == alt[left]:
        left += 1
    end_ref, end_alt = len(ref), len(alt)
    while end_ref > left and end_alt > left and ref[end_ref - 1] == alt[end_alt - 1]:
        end_ref -= 1
        end_alt -= 1
    # Add the required VCF left anchor for insertions/deletions.
    left -= 1
    rp, vp = tmp_path / "reference.fa", tmp_path / "raw.vcf"
    rp.write_text(">contig_3\n" + ref + "\n")
    run_tool(["samtools", "faidx", str(rp)])
    vp.write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID=contig_3,length={len(ref)}>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        f"contig_3\t{left + 1}\t.\t{ref[left:end_ref]}\t{alt[left:end_alt]}\t30\tPASS\t.\tGT"
        f"\t{genotype}\n"
    )
    filtered = filter_vcf(vp, rp, tmp_path / "normalized")
    full, trimmed = tmp_path / "full.fa", tmp_path / "trimmed.fa"
    build_consensus(rp, filtered, full, sample="S", haplotype="I")
    context = {
        "reference_path": str(rp),
        "full_consensus_path": str(full),
        "chrom": "contig_3",
        "sample": "S",
        "haplotype": "I",
    }
    trim_flanking(full, 4, trimmed, context=context)
    sequence = "".join(
        line for line in trimmed.read_text().splitlines() if not line.startswith(">")
    )
    result = validate_mutations_against_vcf(
        classify_sequence(sequence, rd),
        parse_vcf_variants(filtered),
        sequence=sequence,
        repeat_dict=rd,
        consensus_context=context,
    )
    return sequence, rd.repeats["X"] + mutant + rd.repeats["X"], result


@pytest.mark.parametrize("name", ["dupC", "dupA", "insG", "insCCCC", "del18_31"])
def test_normalized_named_mutation_support(name: str, tmp_path: Path) -> None:
    sequence, expected, result = _named_mutation_case(name, "1/1", tmp_path)
    assert sequence == expected
    matching = [m for m in result["mutations_detected"] if m.get("mutation_name") == name]
    assert len(matching) == 1
    assert matching[0]["vcf_support"] is True


@pytest.mark.parametrize("name", ["dupC", "del18_31"])
def test_iupac_consensus_applies_heterozygous_indel_as_unresolved_event(
    name: str, tmp_path: Path
) -> None:
    """Guards the installed bcftools: -H I applies a 0/1 indel's ALT (verified on 1.17)."""
    sequence, expected, result = _named_mutation_case(name, "0/1", tmp_path)
    assert sequence == expected
    assert result["vcf_projection"]["status"] == "available"
    matching = [m for m in result["mutations_detected"] if m.get("mutation_name") == name]
    assert len(matching) == 1
    assert matching[0]["vcf_support_status"] == "heterozygous_genotype_unresolved"
    assert matching[0]["vcf_support"] is False
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run --locked --all-extras pytest tests/unit/test_variant_support.py tests/unit/test_report.py --no-cov -q`
Expected: 3 tests fail.
- `test_unrelated_heterozygous_indel_keeps_per_event_support` fails with `'unavailable' == 'available'`.
- `test_heterozygous_event_under_iupac_is_unresolved_not_supported` fails the same way.
- `test_report_distinguishes_unavailable_from_absent_support` fails on the badge count.

- [ ] **Step 3: Implement**

In `src/muc_one_span/variant_support.py`, `_selected_alt`: replace the final `return None`, the one after the IUPAC branch, with:

```python
    # bcftools consensus -H I applies the ALT of a heterozygous REF/ALT indel
    # (verified with bcftools 1.17; tests/integration guards the installed tool).
    # Multi-ALT heterozygous indels and MNPs have no verified mapping.
    alts = [a for a in selected if a != variant["ref"]]
    if len(alts) == 1 and len(alts[0]) != len(variant["ref"]):
        return alts[0]
    return None
```

In `_replay`, replace `applied.append({**variant, "selected_alt": alt, "query_start": query_pos})` with:

```python
        genotype = set(variant["genotype"].replace("|", "/").split("/"))
        unresolved = context.get("haplotype", "I") == "I" and len(genotype) > 1
        applied.append(
            {
                **variant,
                "selected_alt": alt,
                "query_start": query_pos,
                "unresolved_genotype": unresolved,
            }
        )
```

In `mutation_concordance`, replace `projection.update(status="available", reason="exact_full_consensus_replay")` with:

```python
    projection.update(
        status="available",
        reason="exact_full_consensus_replay",
        unresolved_genotype_edits=sum(edit["unresolved_genotype"] for edit in edits),
    )
```

Also replace `result[index] = ("exact_sequence_concordance" if supporting else "absent", supporting)` with:

```python
        resolved = [edit for edit in supporting if not edit["unresolved_genotype"]]
        if resolved:
            result[index] = ("exact_sequence_concordance", resolved)
        elif supporting:
            result[index] = ("heterozygous_genotype_unresolved", [])
        else:
            result[index] = ("absent", [])
```

`src/muc_one_span/templates/report.html.j2`: after the `localization_ambiguous` badge branch, add:

```jinja
          {% elif m.get("vcf_support_status") == "heterozygous_genotype_unresolved" %}
            <span class="badge badge-info">Genotype unresolved</span>
```

In `docs/development.md` line 242, replace "Unprojectable indels remain unresolved." with:

```markdown
Under `-H I`, bcftools writes IUPAC codes for heterozygous SNVs but applies the
ALT allele of a heterozygous REF/ALT indel. Replay mirrors this, and such an
event gets `vcf_support_status=heterozygous_genotype_unresolved` (not
supported). Other events on the allele keep their own status, and
`vcf_projection.unresolved_genotype_edits` counts these edits. Multi-ALT
heterozygous indels remain unprojectable (`ambiguous_genotype_selection`).
```

- [ ] **Step 4: Run the tests again, including integration**

Run: `uv run --locked --all-extras pytest tests/unit/test_variant_support.py tests/unit/test_report.py tests/unit/test_classify.py --no-cov -q`. Expected: all pass.
Run: `make test-int` with the tool environment active. Expected: `test_variant_concordance.py` has 7 tests passing. If the installed bcftools is not 1.17, record its version in the PR. A failure of `test_iupac_consensus_applies_heterozygous_indel_as_unresolved_event` means that bcftools changed its `-H I` indel behaviour. In that case stop and do not ship.
Then `make ci-check` must pass.

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/variant_support.py src/muc_one_span/templates/report.html.j2 \
  docs/development.md tests/unit/test_variant_support.py tests/unit/test_report.py \
  tests/integration/test_variant_concordance.py
git commit -m "fix(support): evaluate VCF concordance per event under IUPAC consensus

bcftools consensus -H I applies heterozygous REF/ALT indels as ALT, but replay
declared them unresolvable and removed support from every event on the allele.
Mirror bcftools, mark events explained only by such edits as
heterozygous_genotype_unresolved (not supported) and keep exact support for
other events. Add a real-bcftools integration guard.

Fixes #64
Refs #59

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Hybrid-plan Task 1 — explicit support for PATHOGENIC (#55)

**Files:**
- Create: `src/muc_one_span/clinical_gates.py`
- Modify: `src/muc_one_span/report.py:22` (import), `report.py:93-96` (`supp_ok`)
- Modify: `tests/unit/test_clinical_decision.py` (append), `tests/unit/test_report_wave1.py:97-99`

**Interfaces:**
- Produces: `clinical_gates.mutation_supported(mutation: dict[str, Any]) -> bool` and `clinical_gates.SUPPORTED_VCF_STATUSES`.
- The name `report._mutation_supported` exists only in this task. Task 5 drops the alias, and the hybrid plan relies on the behaviour, not on the name.
- The test helpers `_summary(mut)` and `BASE` in `test_clinical_decision.py` are exactly those named in hybrid-plan Task 1/Task 10, so the hybrid Task 10 tests can append to them unchanged.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_clinical_decision.py`. These are the hybrid-plan Task 1 tests, adapted to the real key `state`:

```python
def _summary(mut: dict) -> dict:
    allele = {
        "length": 50,
        "canonical_repeats": 41,
        "reads": 200,
        "phase_status": "phased",
        "independent_haplotype_evidence": True,
    }
    return {
        "alleles": {
            "allele_1": dict(allele),
            "allele_2": dict(allele, length=60, canonical_repeats=51),
            "homozygous": False,
        },
        "classifications": {
            "allele_1": {
                "mutations": [mut],
                "ambiguous_bases": 0,
                "reconstruction_status": "complete_segmentation",
            },
            "allele_2": {
                "mutations": [],
                "ambiguous_bases": 0,
                "reconstruction_status": "complete_segmentation",
            },
        },
    }


BASE = {
    "repeat_index": 20,
    "mutation_name": "dupC",
    "frameshift": True,
    "localization_status": "exact",
    "template_match": True,
}


def test_missing_support_fields_are_not_pathogenic() -> None:
    decision = compute_clinical_decision(_summary(dict(BASE)))
    assert decision["state"] != "PATHOGENIC"


def test_read_support_is_accepted_as_explicit_support() -> None:
    mut = dict(
        BASE,
        vcf_support=False,
        vcf_support_status="not_applicable_read_consensus",
        read_support={"status": "supported", "alt": 180, "ref": 10, "other": 10},
    )
    assert compute_clinical_decision(_summary(mut))["state"] == "PATHOGENIC"


def test_vcf_exact_concordance_still_pathogenic() -> None:
    mut = dict(BASE, vcf_support=True, vcf_support_status="exact_sequence_concordance")
    assert compute_clinical_decision(_summary(mut))["state"] == "PATHOGENIC"


def test_unresolved_heterozygous_event_is_not_pathogenic() -> None:
    mut = dict(BASE, vcf_support=False, vcf_support_status="heterozygous_genotype_unresolved")
    assert compute_clinical_decision(_summary(mut))["state"] == "INCONCLUSIVE"
```

In `tests/unit/test_report_wave1.py`, `test_mutation_evidence_remains_visible_with_execution_warning`, replace the mutation list with a mutation that carries explicit support and identity. The test's intent, "mutation evidence stays visible next to an execution warning", is unchanged:

```python
    data["classifications"]["allele_1"]["mutations"] = [
        {
            "mutation_name": "dupC",
            "repeat_index": 8,
            "frameshift": True,
            "template_match": True,
            "vcf_support": True,
            "vcf_support_status": "exact_sequence_concordance",
        }
    ]
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run --locked --all-extras pytest tests/unit/test_clinical_decision.py tests/unit/test_report_wave1.py --no-cov -q`
Expected: 2 tests fail.
- `test_missing_support_fields_are_not_pathogenic` fails with `'PATHOGENIC' != 'PATHOGENIC'`.
- `test_read_support_is_accepted_as_explicit_support` fails with `'INCONCLUSIVE' == 'PATHOGENIC'`.

- [ ] **Step 3: Implement**

Create `src/muc_one_span/clinical_gates.py`:

```python
"""Evidence gates applied before a clinical decision banner is chosen.

Pure functions over summary dictionaries: they never read files or run tools.
Missing evidence is never treated as support.
"""

from __future__ import annotations

from typing import Any

SUPPORTED_VCF_STATUSES = frozenset({"exact_sequence_concordance"})


def mutation_supported(mutation: dict[str, Any]) -> bool:
    """Return True only for explicit sequence-level support (VCF or read evidence).

    ``read_support.status == "supported"`` is accepted for read-level evidence.
    A legacy ``vcf_support=True`` without a status keeps its historical meaning;
    a missing ``vcf_support`` key is not support.
    """
    read_support = mutation.get("read_support")
    if isinstance(read_support, dict) and read_support.get("status") == "supported":
        return True
    if mutation.get("vcf_support") is True:
        status = mutation.get("vcf_support_status")
        return status is None or status in SUPPORTED_VCF_STATUSES
    return False
```

In `src/muc_one_span/report.py`, add `from muc_one_span.clinical_gates import mutation_supported as _mutation_supported` before the `nomenclature` import. Then replace the four-line `supp_ok = (...)` expression with:

```python
                    supp_ok = _mutation_supported(mut)
```

- [ ] **Step 4: Run the whole clinical and report suites**

Run: `uv run --locked --all-extras pytest tests/unit/test_clinical_decision.py tests/unit/test_report.py tests/unit/test_report_wave1.py tests/unit/test_report_igv.py --no-cov -q`
Expected: all pass. Then `make ci-check` must pass.

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/clinical_gates.py src/muc_one_span/report.py \
  tests/unit/test_clinical_decision.py tests/unit/test_report_wave1.py
git commit -m "fix(report): require explicit VCF or read-level support for pathogenic decisions

A mutation without a vcf_support key counted as supported. Only exact VCF
concordance, legacy vcf_support=True without status, or read_support.status
'supported' is now support. The wave-1 warning test now uses an explicitly
supported mutation; its intent (evidence stays visible) is unchanged.

Refs #55

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: D8 + D3 — gate order, event identity, allele selection, depth and length (#55)

**Files:**
- Create: `src/muc_one_span/selection_qc.py`
- Modify: `src/muc_one_span/clinical_gates.py` (append)
- Modify: `src/muc_one_span/report.py:22` (imports), `report.py:78-221` (`compute_clinical_decision`)
- Modify: `src/muc_one_span/settings.py:83-110` (`AlleleSelectionSettings`)
- Modify: `src/muc_one_span/pipeline.py:39` (import), `pipeline.py:110-111` (call)
- Create: `tests/unit/test_selection_qc.py`, `tests/unit/test_pipeline_gates.py`
- Modify: `tests/unit/test_clinical_decision.py` (fixture at line 38 plus appended tests), `tests/unit/test_report_igv.py:180-185`, `tests/unit/test_runtime_settings.py`
- Modify: `docs/guides/configuration.md:79-82`, `docs/development.md:318`

**Interfaces:**
- Consumes: `alleles[*].fit_metrics[contig].primary_alignment_records`, `primary_alignment_records`, `length`, `reference_length` and `length_selection_evidence.unselected_passing_clusters` (all existing), plus `allele_genotype_status` from Task 1 and `vcf_support_status` from Task 3.
- Produces:
  - `selection_qc.secondary_mode_fraction(fit_metrics: dict[str, dict[str, Any]], min_gap: int) -> float | None`;
  - `selection_qc.assess_allele(info: dict[str, Any], settings: AlleleSelectionSettings) -> dict[str, Any]`;
  - `selection_qc.annotate_selection_qc(alleles: dict[str, Any], settings: AlleleSelectionSettings | None = None) -> dict[str, Any]`.
- New additive per-allele `alleles.json` fields: `selection_status`, `secondary_mode_fraction`, `depth_status`, `depth_basis`, `depth_threshold` and `length_status`.
  - `selection_status` is one of `resolved`, `unresolved_secondary_mode`, `unresolved_unselected_clusters` or `not_assessed`.
  - `depth_status` is `adequate`, `low` or `not_assessed`. Hybrid-plan Task 10 reuses the name and meaning.
- Produces: `clinical_gates.mutation_blockers(mutation) -> list[str]`, `clinical_gates.allele_gate_reasons(info: Any, label: str) -> list[str]` and `clinical_gates.LEGACY_MIN_TOTAL_READS = 30`.
- The decision dict keys are unchanged.

**Gate order in `compute_clinical_decision`:**

1. A mutation is PATHOGENIC-eligible only when all of these hold:
   - `frameshift` is true;
   - event identity is established (`template_match` true and `mutation_name` set);
   - localization is not ambiguous;
   - support is explicit;
   - the carrying allele's `depth_status` is not `low`;
   - for legacy summaries without per-allele depth, the total is not below 30 reads.
2. PATHOGENIC keeps its established behaviour of retaining evidence next to execution warnings. It now lists allele-selection and reconstruction problems as `Quality caveat:` details.
3. NEGATIVE also requires that no allele gate reason exists: resolved selection, reported length equal to the consensus contig (D3), no unresolved length-partition genotype, and depth not low. The existing execution, ambiguity and reconstruction checks still apply.

Legacy summaries without the new fields are gated only on fields they already carry (`length`/`reference_length`, `reads`).

- [ ] **Step 1: Write the failing tests**

Create `tests/unit/test_selection_qc.py`:

```python
"""Allele-selection, depth and length gates from ladder primary-alignment evidence."""

from __future__ import annotations

import pytest

from muc_one_span.selection_qc import annotate_selection_qc, assess_allele, secondary_mode_fraction
from muc_one_span.settings import AlleleSelectionSettings


def _metrics(counts: dict[int, int]) -> dict[str, dict[str, int]]:
    return {f"contig_{unit}": {"primary_alignment_records": n} for unit, n in counts.items()}


# Shapes of PRJEB92208 cohort-v3 clusters (public libraries): a control peak with a
# 1% distant tail, and MP3 allele 1 merging a second allele mode at 0.62 of the peak.
CONTROL = _metrics({6: 36, 27: 284, 28: 2855, 29: 3510, 30: 392})
MERGED = _metrics({32: 2754, 33: 2872, 34: 400, 38: 1778, 39: 1679})


def test_secondary_mode_fraction_separates_control_and_merged_cluster() -> None:
    assert secondary_mode_fraction(CONTROL, 5) == pytest.approx(36 / 3510)
    assert secondary_mode_fraction(MERGED, 5) == pytest.approx(1778 / 2872)
    assert secondary_mode_fraction({}, 5) is None
    assert secondary_mode_fraction(_metrics({10: 0}), 5) is None


def test_assess_allele_reports_all_gates() -> None:
    settings = AlleleSelectionSettings()
    base = {"length": 45, "reference_length": 45, "primary_alignment_records": 632}
    assert assess_allele({**base, "fit_metrics": CONTROL}, settings) == {
        "selection_status": "resolved",
        "secondary_mode_fraction": round(36 / 3510, 4),
        "depth_status": "adequate",
        "depth_basis": "primary_alignment_records",
        "depth_threshold": 30,
        "length_status": "consistent_with_consensus_contig",
    }
    merged = assess_allele({**base, "fit_metrics": MERGED}, settings)
    assert merged["selection_status"] == "unresolved_secondary_mode"
    low = assess_allele({**base, "fit_metrics": CONTROL, "primary_alignment_records": 12}, settings)
    assert low["depth_status"] == "low"
    shifted = assess_allele(
        {**base, "length": 39, "reference_length": 44, "fit_metrics": CONTROL}, settings
    )
    assert shifted["length_status"] == "cluster_center_differs_from_consensus_contig"
    unassessed = assess_allele({"length": 45}, settings)
    assert (
        unassessed["selection_status"],
        unassessed["depth_status"],
        unassessed["length_status"],
    ) == ("not_assessed", "not_assessed", "not_assessed")


def test_unselected_passing_cluster_marks_selection_unresolved() -> None:
    info = {
        "fit_metrics": CONTROL,
        "primary_alignment_records": 632,
        "length_selection_evidence": {"unselected_passing_clusters": 1},
    }
    status = assess_allele(info, AlleleSelectionSettings())["selection_status"]
    assert status == "unresolved_unselected_clusters"


def test_gate_thresholds_come_from_settings() -> None:
    settings = AlleleSelectionSettings(
        secondary_mode_min_fraction=0.7, min_allele_primary_records=700
    )
    result = assess_allele({"fit_metrics": MERGED, "primary_alignment_records": 632}, settings)
    assert (result["selection_status"], result["depth_status"]) == ("resolved", "low")
    assert result["depth_threshold"] == 700


def test_annotate_updates_both_alleles_in_place() -> None:
    alleles = {
        "allele_1": {"fit_metrics": CONTROL, "primary_alignment_records": 632},
        "allele_2": {"fit_metrics": MERGED, "primary_alignment_records": 80},
        "homozygous": False,
    }
    assert annotate_selection_qc(alleles) is alleles
    assert alleles["allele_1"]["selection_status"] == "resolved"
    assert alleles["allele_2"]["selection_status"] == "unresolved_secondary_mode"
    assert alleles["homozygous"] is False
```

Create `tests/unit/test_pipeline_gates.py`. Task 7 reuses its helper:

```python
"""Pipeline wiring of pre-analysis checks and allele-selection gates (all tools mocked)."""

from __future__ import annotations

import json
from collections.abc import Iterable
from contextlib import AbstractContextManager, ExitStack
from pathlib import Path
from typing import Any
from unittest.mock import patch

from click.testing import CliRunner, Result

from muc_one_span.cli import main


def _alleles() -> dict[str, Any]:
    def allele(units: int, primary: int) -> dict[str, Any]:
        return {
            "length": units + 9,
            "reference_length": units + 9,
            "reads": primary * 6,
            "primary_alignment_records": primary,
            "canonical_repeats": units,
            "contig_name": f"contig_{units}",
            "cluster_contigs": [f"contig_{units}"],
            "fit_metrics": {f"contig_{units}": {"primary_alignment_records": primary}},
        }

    return {
        "homozygous": False,
        "same_length": False,
        "allele_1": allele(51, 12),
        "allele_2": allele(71, 90),
    }


def run_mocked_pipeline(
    tmp_path: Path,
    *args: str,
    extra_patches: Iterable[AbstractContextManager[Any]] = (),
) -> tuple[Result, list[str]]:
    """Invoke ``run`` with every tool stage mocked; record check/map call order."""
    input_file = tmp_path / "reads.fastq"
    input_file.touch()
    fasta = tmp_path / "allele.fa"
    fasta.write_text(">allele_vntr\nACGT\n")
    classification = {
        "structure": "1 2 3 X 6 7 8 9",
        "repeats": [],
        "mutations_detected": [],
        "allele_confidence": 1.0,
    }
    calls: list[str] = []

    def check_tools(tools: list[str]) -> bool:
        calls.append("check:" + ",".join(tools))
        return True

    def map_reads(*_args: Any, **_kwargs: Any) -> Path:
        calls.append("map")
        return tmp_path / "mapping.bam"

    with ExitStack() as stack:
        for manager in (
            patch("muc_one_span.tools.check_tools", side_effect=check_tools),
            patch("muc_one_span.tools.get_tool_versions", return_value={}),
            patch("muc_one_span.mapping.map_reads", side_effect=map_reads),
            patch("muc_one_span.mapping.get_idxstats", return_value="contig_51\t3120\t72\t0\n"),
            patch("muc_one_span.alleles.parse_idxstats", return_value={51: 72, 71: 540}),
            patch("muc_one_span.alleles.detect_alleles", side_effect=lambda *a, **k: _alleles()),
            patch(
                "muc_one_span.calling.call_variants_per_allele",
                return_value={
                    "allele_1": tmp_path / "a1.vcf.gz",
                    "allele_2": tmp_path / "a2.vcf.gz",
                },
            ),
            patch(
                "muc_one_span.consensus.build_consensus_per_allele",
                return_value={"allele_1": fasta, "allele_2": fasta},
            ),
            patch("muc_one_span.classify.classify_sequence", return_value=classification),
            patch(
                "muc_one_span.classify.validate_mutations_against_vcf", return_value=classification
            ),
            patch("muc_one_span.vcf.parse_vcf_variants", return_value=[]),
            *extra_patches,
        ):
            stack.enter_context(manager)
        result = CliRunner().invoke(
            main, ["run", "--input", str(input_file), "--output-dir", str(tmp_path / "out"), *args]
        )
    return result, calls


def test_alleles_json_records_selection_and_depth_gates(tmp_path: Path) -> None:
    result, _ = run_mocked_pipeline(tmp_path)
    assert result.exit_code == 0, result.output
    alleles = json.loads((tmp_path / "out" / "alleles.json").read_text())
    assert alleles["allele_1"]["depth_status"] == "low"
    assert alleles["allele_1"]["depth_threshold"] == 30
    assert alleles["allele_2"]["depth_status"] == "adequate"
    assert alleles["allele_2"]["selection_status"] == "resolved"
    assert alleles["allele_2"]["length_status"] == "consistent_with_consensus_contig"
```

Append to `tests/unit/test_runtime_settings.py`:

```python
def test_selection_gate_settings_are_validated() -> None:
    selection = AlleleSelectionSettings()
    assert (selection.secondary_mode_min_fraction, selection.min_allele_primary_records) == (
        0.2,
        30,
    )
    with pytest.raises(ValueError, match="secondary_mode_min_fraction must be > 0"):
        AlleleSelectionSettings(secondary_mode_min_fraction=0)
    with pytest.raises(ValueError, match="min_allele_primary_records"):
        AlleleSelectionSettings(min_allele_primary_records=0)
```

In `tests/unit/test_clinical_decision.py`, `test_supported_frameshift_yields_pathogenic`, add `"template_match": True,` to the mutation dict after `"frameshift": True,`. PATHOGENIC now requires event identity. Then append:

```python
SUPPORTED = dict(BASE, vcf_support=True, vcf_support_status="exact_sequence_concordance")
RESOLVED_GATES = {
    "selection_status": "resolved",
    "secondary_mode_fraction": 0.01,
    "depth_status": "adequate",
    "depth_threshold": 30,
    "primary_alignment_records": 400,
    "length_status": "consistent_with_consensus_contig",
}


def _gated_summary(mutations: list[dict] | None = None) -> dict:
    summary = _summary(dict(BASE))
    summary["classifications"]["allele_1"]["mutations"] = mutations or []
    summary["run_status"] = {"status": "completed"}
    for key in ("allele_1", "allele_2"):
        summary["alleles"][key].update(RESOLVED_GATES)
    return summary


def test_resolved_gates_keep_negative() -> None:
    decision = compute_clinical_decision(_gated_summary())
    assert decision["state"] == "NO_PATHOGENIC_VARIANT_DETECTED"


def test_secondary_mode_blocks_negative() -> None:
    """MP3-shaped: a second primary-record mode at 0.62 of the peak inside one cluster."""
    summary = _gated_summary()
    summary["alleles"]["allele_1"].update(
        selection_status="unresolved_secondary_mode", secondary_mode_fraction=0.619
    )
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"
    assert any("allele selection unresolved" in detail for detail in decision["details"])


def test_low_allele_depth_blocks_negative_and_pathogenic() -> None:
    negative = _gated_summary()
    negative["alleles"]["allele_2"].update(depth_status="low", primary_alignment_records=12)
    assert compute_clinical_decision(negative)["state"] == "INCONCLUSIVE"
    positive = _gated_summary([dict(SUPPORTED)])
    positive["alleles"]["allele_1"].update(depth_status="low", primary_alignment_records=12)
    decision = compute_clinical_decision(positive)
    assert decision["state"] == "INCONCLUSIVE"
    assert any("per-allele depth gate" in detail for detail in decision["details"])


def test_legacy_low_total_reads_blocks_pathogenic() -> None:
    summary = _summary(dict(SUPPORTED))
    summary["alleles"]["allele_1"]["reads"] = 10
    summary["alleles"]["allele_2"]["reads"] = 12
    assert compute_clinical_decision(summary)["state"] == "INCONCLUSIVE"


def test_untemplated_frameshift_is_not_pathogenic() -> None:
    """MP5-shaped: a supported novel frameshift without a dictionary template is inconclusive."""
    novel = {
        "repeat_index": 35,
        "closest_type": "X",
        "differences": [{"type": "insertion", "position": 23}],
        "frameshift": True,
        "vcf_support": True,
        "vcf_support_status": "exact_sequence_concordance",
        "localization_status": "resolved",
    }
    decision = compute_clinical_decision(_gated_summary([novel]))
    assert decision["state"] == "INCONCLUSIVE"
    assert any("event identity not established" in detail for detail in decision["details"])


def test_length_contig_mismatch_blocks_negative() -> None:
    """MP1-shaped allele 1: reported 39 units, consensus from the 44-unit contig."""
    summary = _gated_summary()
    summary["alleles"]["allele_1"].update(length=39, reference_length=44)
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "INCONCLUSIVE"
    assert any(
        "reported length 39 differs from the consensus contig length 44" in detail
        for detail in decision["details"]
    )


def test_pathogenic_lists_quality_caveats() -> None:
    summary = _gated_summary([dict(SUPPORTED)])
    summary["alleles"]["allele_2"].update(
        selection_status="unresolved_secondary_mode", secondary_mode_fraction=0.28
    )
    decision = compute_clinical_decision(summary)
    assert decision["state"] == "PATHOGENIC"
    assert decision["details"][0].startswith("Allele 1: dupC at repeat unit 20")
    assert any(detail.startswith("Quality caveat: Allele 2") for detail in decision["details"])
```

In `tests/unit/test_report_igv.py`, `test_generate_report_with_igv_and_hgvs`, add `"template_match": True,` to the `59dupC` mutation after `"frameshift": True,`. The test asserts the PATHOGENIC banner's HGVS detail.

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run --locked --all-extras pytest tests/unit/test_selection_qc.py tests/unit/test_pipeline_gates.py tests/unit/test_clinical_decision.py tests/unit/test_runtime_settings.py --no-cov -q`
Expected:
- `test_selection_qc.py` fails with `ModuleNotFoundError: muc_one_span.selection_qc`.
- `test_alleles_json_records_selection_and_depth_gates` fails with `KeyError: 'depth_status'`.
- `test_selection_gate_settings_are_validated` fails with a `TypeError`.
- Six decision tests fail: `secondary_mode`, `low_allele_depth`, `legacy_low_total_reads`, `untemplated_frameshift`, `length_contig_mismatch` and `pathogenic_lists_quality_caveats`.

- [ ] **Step 3: Implement**

Create `src/muc_one_span/selection_qc.py`:

```python
"""Allele-selection and depth gates computed from ladder alignment evidence.

Primary alignment records per contig approximate molecules (one primary record
per read). A cluster whose primary records form a second mode at least
``min_gap`` units from its peak may hide another allele or a fragment class, so
its selection is unresolved. These gates only prevent a reassuring negative;
they never create a variant call.
"""

from __future__ import annotations

from typing import Any

from muc_one_span.settings import DEFAULT_SETTINGS, AlleleSelectionSettings


def _contig_units(name: str) -> int:
    return int(name.rsplit("_", 1)[1])


def secondary_mode_fraction(fit_metrics: dict[str, dict[str, Any]], min_gap: int) -> float | None:
    """Largest primary-record count >= ``min_gap`` units from the peak, relative to the peak."""
    counts = {
        _contig_units(contig): int(metrics.get("primary_alignment_records") or 0)
        for contig, metrics in fit_metrics.items()
    }
    if not counts or max(counts.values()) == 0:
        return None
    peak = max(counts, key=lambda contig: (counts[contig], -contig))
    distant = [n for contig, n in counts.items() if abs(contig - peak) >= min_gap]
    return max(distant, default=0) / counts[peak]


def assess_allele(info: dict[str, Any], settings: AlleleSelectionSettings) -> dict[str, Any]:
    """Return additive selection, depth and length-consistency fields for one allele."""
    fraction = secondary_mode_fraction(info.get("fit_metrics") or {}, settings.min_gap)
    unselected = (info.get("length_selection_evidence") or {}).get("unselected_passing_clusters")
    if fraction is None:
        selection = "not_assessed"
    elif fraction >= settings.secondary_mode_min_fraction:
        selection = "unresolved_secondary_mode"
    elif unselected:
        selection = "unresolved_unselected_clusters"
    else:
        selection = "resolved"
    primary = info.get("primary_alignment_records")
    if primary is None:
        depth = "not_assessed"
    else:
        depth = "adequate" if primary >= settings.min_allele_primary_records else "low"
    length, reference_length = info.get("length"), info.get("reference_length")
    if length is None or reference_length is None:
        length_status = "not_assessed"
    elif length == reference_length:
        length_status = "consistent_with_consensus_contig"
    else:
        length_status = "cluster_center_differs_from_consensus_contig"
    return {
        "selection_status": selection,
        "secondary_mode_fraction": None if fraction is None else round(fraction, 4),
        "depth_status": depth,
        "depth_basis": "primary_alignment_records",
        "depth_threshold": settings.min_allele_primary_records,
        "length_status": length_status,
    }


def annotate_selection_qc(
    alleles: dict[str, Any], settings: AlleleSelectionSettings | None = None
) -> dict[str, Any]:
    """Add selection/depth/length gates to ``allele_1``/``allele_2`` in place."""
    settings = settings or DEFAULT_SETTINGS.allele_selection
    for key in ("allele_1", "allele_2"):
        info = alleles.get(key)
        if isinstance(info, dict):
            info.update(assess_allele(info, settings))
    return alleles
```

`src/muc_one_span/settings.py`, class `AlleleSelectionSettings`: add these fields after `min_dominance_ratio: float = 0.01`:

```python
    secondary_mode_min_fraction: float = 0.2
    min_allele_primary_records: int = 30
```

and these validation lines at the end of `__post_init__`:

```python
        _number(
            "allele_selection.secondary_mode_min_fraction",
            self.secondary_mode_min_fraction,
            0.0,
            1.0,
        )
        if self.secondary_mode_min_fraction == 0:
            raise ValueError("allele_selection.secondary_mode_min_fraction must be > 0")
        _integer("allele_selection.min_allele_primary_records", self.min_allele_primary_records, 1)
```

`src/muc_one_span/pipeline.py`: add `from muc_one_span.selection_qc import annotate_selection_qc` to the function-local imports after the `mapping` import. Then, immediately after the `detect_alleles(...)` call and before the first `alleles.json` write, add:

```python
    annotate_selection_qc(alleles_result, settings.allele_selection)
```

Append to `src/muc_one_span/clinical_gates.py`:

```python
LEGACY_MIN_TOTAL_READS = 30

_GENOTYPE_REASONS = {
    "heterozygous_within_length_partition": (
        "heterozygous call left within the length-partitioned allele; "
        "consensus uses unresolved (IUPAC) selection"
    ),
    "unresolved_genotype_records": "conflicting or incomplete genotype records",
}


def mutation_blockers(mutation: dict[str, Any]) -> list[str]:
    """List every reason an observed mutation cannot support a PATHOGENIC banner."""
    blockers: list[str] = []
    if mutation.get("frameshift") is not True:
        blockers.append("frameshift not established")
    if mutation.get("template_match") is not True or not mutation.get("mutation_name"):
        blockers.append("event identity not established (no exact dictionary template)")
    if mutation.get("localization_status") == "ambiguous":
        blockers.append("localization ambiguous")
    if not mutation_supported(mutation):
        status = mutation.get("vcf_support_status")
        blockers.append(
            "heterozygous genotype unresolved within its allele"
            if status == "heterozygous_genotype_unresolved"
            else f"no explicit sequence-level support ({status or 'status unavailable'})"
        )
    return blockers


def allele_gate_reasons(info: Any, label: str) -> list[str]:
    """Reasons an allele's selection, genotype, length or depth prevents a negative call."""
    if not isinstance(info, dict) or not info:
        return []
    reasons: list[str] = []
    selection = info.get("selection_status")
    if isinstance(selection, str) and selection.startswith("unresolved"):
        reasons.append(
            f"{label}: allele selection unresolved ({selection}; secondary mode fraction "
            f"{info.get('secondary_mode_fraction')})."
        )
    length, reference_length = info.get("length"), info.get("reference_length")
    if isinstance(length, int) and isinstance(reference_length, int) and length != reference_length:
        reasons.append(
            f"{label}: reported length {length} differs from the consensus contig length "
            f"{reference_length}."
        )
    if info.get("depth_status") == "low":
        reasons.append(
            f"{label}: {info.get('primary_alignment_records')} primary alignments, below the "
            f"per-allele depth gate ({info.get('depth_threshold')})."
        )
    genotype = info.get("allele_genotype_status")
    if genotype in _GENOTYPE_REASONS:
        reasons.append(f"{label}: {_GENOTYPE_REASONS[genotype]}.")
    return reasons
```

`src/muc_one_span/report.py`: replace the Task 4 import line with:

```python
from muc_one_span.clinical_gates import (
    LEGACY_MIN_TOTAL_READS,
    allele_gate_reasons,
    mutation_blockers,
)
```

In `compute_clinical_decision`, replace everything from the docstring through the old `low_coverage = ...` line (lines 78-105 on `main`, one line later after Task 4) with the code below. This keeps `ambiguous_bases` and the `reconstruction_reasons` block that follow.

```python
    """Derive the 3-state clinical decision banner after all evidence gates.

    Gate order: each observed mutation must pass frameshift, event identity,
    localization, explicit support and carrier-allele depth before it can make
    the banner PATHOGENIC. NEGATIVE additionally requires completed execution,
    resolved allele selection, genotype and length, adequate depth and no
    uncertain observed event. Legacy summaries without per-allele depth fall
    back to the total-read threshold.
    """
    execution = _execution_context(summary, execution_status)
    classifications = summary.get("classifications", {})
    alleles = summary.get("alleles", {})
    a1 = alleles.get("allele_1", {}) if isinstance(alleles, dict) else {}
    a2 = alleles.get("allele_2", {}) if isinstance(alleles, dict) else {}
    carriers = {"allele_1": a1, "allele_2": a2}
    depth_assessed = any(a.get("depth_status") in ("adequate", "low") for a in (a1, a2))
    total_reads = (a1.get("reads", 0) or 0) + (a2.get("reads", 0) or 0)
    low_coverage = (
        not depth_assessed and total_reads < LEGACY_MIN_TOTAL_READS and (bool(a1) or bool(a2))
    )

    pathogenic_mutations: list[dict[str, Any]] = []
    uncertain_mutations: list[dict[str, Any]] = []
    for allele_key, acls in classifications.items():
        if not isinstance(acls, dict):
            continue
        carrier = carriers.get(allele_key) or {}
        for mut in acls.get("mutations", []):
            if not isinstance(mut, dict):
                continue
            mut_copy = dict(mut)
            mut_copy["allele"] = allele_key
            blockers = mutation_blockers(mut)
            if carrier.get("depth_status") == "low":
                blockers.append("carrying allele is below the per-allele depth gate")
            if low_coverage:
                blockers.append("total read depth is below the diagnostic threshold")
            mut_copy["decision_blockers"] = blockers
            (uncertain_mutations if blockers else pathogenic_mutations).append(mut_copy)
```

Then:

- (a) directly after the `reconstruction_reasons` loop and before `if pathogenic_mutations:`, add:

```python
    selection_reasons = allele_gate_reasons(a1, "Allele 1") + allele_gate_reasons(a2, "Allele 2")
```

- (b) in the PATHOGENIC branch, after the loop that appends "Additional uncertain variant ..." details, add:

```python
        details.extend(
            f"Quality caveat: {reason}" for reason in selection_reasons + reconstruction_reasons
        )
```

- (c) in the INCONCLUSIVE condition, add `or bool(selection_reasons)` after `or bool(reconstruction_reasons)`.
- (d) replace the low-coverage reason with:

```python
            reasons.append(
                f"Total read depth ({total_reads} reads) is below diagnostic threshold "
                f"({LEGACY_MIN_TOTAL_READS} reads)."
            )
```

- (e) replace the uncertain-mutation reason, and add the selection reasons after the reconstruction reasons:

```python
            blocker_text = "; ".join(m.get("decision_blockers", []))
            reasons.append(
                f"{allele_name}: Observed sequence variant ({m_name} at repeat {rep_idx}) "
                f"is inconclusive ({blocker_text or 'in-frame or ambiguous localization/support'})."
            )
        if reconstruction_reasons:
            reasons.extend(reconstruction_reasons)
        reasons.extend(selection_reasons)
```

`docs/guides/configuration.md`: after the `allele_selection.refinement_max_shift` row, add:

```markdown
| `allele_selection.secondary_mode_min_fraction` | `0.2` | Clinical gate. Allele selection is unresolved when a cluster's primary-alignment count at least `min_gap` units from its peak reaches this fraction of the peak. Number in (0,1]. Blocks a negative result; never creates a call. |
| `allele_selection.min_allele_primary_records` | `30` | Clinical gate. Minimum primary alignment records (a molecule proxy) per selected allele. Below it, `depth_status` is `low`, which blocks NEGATIVE and PATHOGENIC. Integer >=1. |
```

`docs/development.md`: after line 318 ("Recorded mutation evidence is retained with execution warnings."), add:

```markdown
Clinical gates run before the banner is chosen (`clinical_gates.py`). A mutation
supports PATHOGENIC only when all of these hold:

- it is a frameshift;
- it is an exact dictionary template (`template_match` and `mutation_name`);
- its localization is not ambiguous;
- it has explicit support (exact VCF concordance or `read_support.status=supported`);
- its allele's `depth_status` is not `low`.

NEGATIVE additionally requires:

- resolved allele selection (`selection_status`);
- reported length equal to the consensus contig length (`length`/`reference_length`);
- no unresolved length-partition genotype (`allele_genotype_status`);
- adequate per-allele depth.

Summaries without per-allele depth fall back to the 30-read total. The gates can
only lower certainty; PATHOGENIC lists remaining problems as quality caveats.
```

- [ ] **Step 4: Run the tests again**

Run: `uv run --locked --all-extras pytest tests/unit --no-cov -q`
Expected: all pass. Then `make ci-check` must pass.

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/selection_qc.py src/muc_one_span/clinical_gates.py \
  src/muc_one_span/report.py src/muc_one_span/settings.py src/muc_one_span/pipeline.py \
  docs/guides/configuration.md docs/development.md tests/unit/test_selection_qc.py \
  tests/unit/test_pipeline_gates.py tests/unit/test_clinical_decision.py \
  tests/unit/test_report_igv.py tests/unit/test_runtime_settings.py
git commit -m "fix(report): gate clinical decisions on event identity, allele selection and depth

PATHOGENIC now requires frameshift, exact template identity, resolved
localization, explicit support and adequate carrier-allele depth. NEGATIVE
additionally requires resolved allele selection (no second primary-record mode,
no unselected passing cluster), a reported length equal to the consensus contig,
a resolved length-partition genotype and adequate per-allele primary depth.
Selection/depth/length QC is recorded in alleles.json (additive fields). Two
fixtures gain template_match because identity is now required.

Fixes #55
Fixes #63
Refs #20

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: D7 — honour calling settings and record the applied QUAL threshold

**Files:**
- Modify: `src/muc_one_span/calling.py` (import `Any`; new `_haploid_filter` before `extract_allele_reads`; haplotag path; distinct path)
- Modify: `src/muc_one_span/vcf.py:83-88` (QUAL rule), and the `filter_vcf` docstring args
- Modify: `src/muc_one_span/settings.py` (`CallingSettings.haploid_min_qual`, `ConsensusSettings` docstring, `load_settings`, `import warnings`)
- Modify: `src/muc_one_span/cli.py:256` and `cli.py:480` (help string only, same line count)
- Modify: `tests/unit/test_vcf.py`, `tests/unit/test_distinct_calling.py`, `tests/unit/test_runtime_settings.py` (append)
- Modify: `docs/guides/configuration.md` (calling and consensus rows)

**Interfaces:**
- Produces: `calling._haploid_filter(min_qual: float, settings: CallingSettings) -> tuple[dict[str, Any], dict[str, Any]]`, returning `(filter_vcf kwargs, provenance)`.
- `CallingSettings.haploid_min_qual: float | None = 4.0`.
- New additive per-allele field `variant_filter = {"min_qual", "min_qual_source", "haploid_majority"}`.
- `filter_vcf` QUAL rule: `haploid_min_qual` if given, else `min_qual`.
- Default filter commands are unchanged: QUAL>=4.0 on length-partitioned VCFs, `min_qual` elsewhere.

- [ ] **Step 1: Write the failing tests**

Append to `tests/unit/test_vcf.py`, and add `from typing import Any` to its imports:

```python
def _qual_filter(tmp_path: Path, **kwargs: Any) -> str:
    """Return the -i expression of the bcftools view call for a non-empty VCF."""

    def side_effect(cmd: list[str]) -> str:
        if cmd[:2] == ["bcftools", "norm"]:
            (tmp_path / "normalized.vcf.gz").write_bytes(b"\x00" * 100)
        return ""

    with patch("muc_one_span.vcf.run_tool", side_effect=side_effect) as run:
        filter_vcf(tmp_path / "input.vcf.gz", tmp_path / "ref.fa", tmp_path, **kwargs)
    view = next(c.args[0] for c in run.call_args_list if c.args[0][:2] == ["bcftools", "view"])
    return view[view.index("-i") + 1]


def test_haploid_rewrite_does_not_silently_lower_min_qual(tmp_path: Path) -> None:
    assert _qual_filter(tmp_path, min_qual=5.0, haploid_majority=True) == "QUAL>=5.0"


def test_haploid_min_qual_applies_without_genotype_rewrite(tmp_path: Path) -> None:
    kwargs = {"min_qual": 5.0, "haploid_majority": False, "haploid_min_qual": 4.0}
    assert _qual_filter(tmp_path, **kwargs) == "QUAL>=4.0"
```

Append to `tests/unit/test_distinct_calling.py`:

```python
def test_default_haploid_filter_is_recorded(tmp_path: Path) -> None:
    alleles, filter_calls = _run_distinct(tmp_path, [])
    assert {call["haploid_min_qual"] for call in filter_calls} == {4.0}
    assert {call["haploid_majority"] for call in filter_calls} == {True}
    assert alleles["allele_2"]["variant_filter"] == {
        "min_qual": 4.0,
        "min_qual_source": "calling.haploid_min_qual",
        "haploid_majority": True,
    }


def test_null_haploid_min_qual_follows_min_qual(tmp_path: Path) -> None:
    settings = CallingSettings(haploid_min_qual=None, haploid_majority=False)
    alleles, filter_calls = _run_distinct(tmp_path, [], min_qual=12.0, settings=settings)
    assert {call["haploid_min_qual"] for call in filter_calls} == {12.0}
    assert {call["haploid_majority"] for call in filter_calls} == {False}
    assert alleles["allele_1"]["variant_filter"] == {
        "min_qual": 12.0,
        "min_qual_source": "run.min_qual",
        "haploid_majority": False,
    }


def test_overridden_explicit_min_qual_is_logged(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    _, filter_calls = _run_distinct(tmp_path, [], min_qual=12.0)
    assert {call["haploid_min_qual"] for call in filter_calls} == {4.0}
    assert "not applied to length-partitioned calls" in caplog.text
```

Append to `tests/unit/test_runtime_settings.py`, and add `import warnings` after `import json`:

```python
def test_null_haploid_min_qual_roundtrips(tmp_path: Path) -> None:
    path = tmp_path / "settings.json"
    path.write_text('{"schema_version": 1, "calling": {"haploid_min_qual": null}}')
    settings = load_settings(path)
    assert settings.calling.haploid_min_qual is None
    path.write_text(json.dumps(settings_as_dict(settings)))
    assert load_settings(path) == settings


def test_deprecated_consensus_haploid_settings_warn_only_when_changed(tmp_path: Path) -> None:
    path = tmp_path / "settings.json"
    path.write_text('{"schema_version": 1, "consensus": {"haploid_min_qual": 3.0}}')
    with pytest.warns(DeprecationWarning, match="consensus.haploid"):
        load_settings(path)
    path.write_text(json.dumps(settings_as_dict(DEFAULT_SETTINGS)))
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        load_settings(path)
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run --locked --all-extras pytest tests/unit/test_vcf.py tests/unit/test_distinct_calling.py tests/unit/test_runtime_settings.py --no-cov -q`
Expected: 7 tests fail.
- `test_haploid_rewrite_does_not_silently_lower_min_qual` fails with `QUAL>=4.0 != QUAL>=5.0`.
- `test_haploid_min_qual_applies_without_genotype_rewrite` fails with `QUAL>=5.0`.
- `test_default_haploid_filter_is_recorded` fails with a `KeyError`.
- `test_null_haploid_min_qual_follows_min_qual` fails in validation of `None`.
- `test_overridden_explicit_min_qual_is_logged` fails because nothing is logged.
- `test_null_haploid_min_qual_roundtrips` fails with a `ValueError`.
- `test_deprecated_consensus_haploid_settings_warn_only_when_changed` fails because no warning is raised.

- [ ] **Step 3: Implement**

`src/muc_one_span/vcf.py`: replace the five-line `effective_qual` block, including the implicit `min(min_qual, 4.0)` branch, with:

```python
    effective_qual = min_qual if haploid_min_qual is None else haploid_min_qual
```

Then add this docstring arg line after the `haploid_majority:` line:

```text
        haploid_min_qual: QUAL threshold for length-partitioned calls; None uses min_qual.
```

`src/muc_one_span/settings.py`:

- Add `import warnings` after `import math`.
- In `CallingSettings`, change the field to `haploid_min_qual: float | None = 4.0`, and replace its validation line with:

```python
        if self.haploid_min_qual is not None:
            _number("calling.haploid_min_qual", self.haploid_min_qual, 0.0)
```

- Replace the `ConsensusSettings` docstring with:

```python
    """Flanking reference extent and exact boundary-anchor search parameters.

    ``haploid_majority`` and ``haploid_min_qual`` are deprecated no-ops kept for
    configuration compatibility; ``calling.haploid_*`` controls haploid calling.
    """
```

- In `load_settings`, directly after `settings = RuntimeSettings(**data)`, add:

```python
    if (settings.consensus.haploid_majority, settings.consensus.haploid_min_qual) != (True, 4.0):
        warnings.warn(
            "consensus.haploid_majority and consensus.haploid_min_qual have no effect; "
            "use calling.haploid_majority and calling.haploid_min_qual.",
            DeprecationWarning,
            stacklevel=2,
        )
```

`src/muc_one_span/calling.py`: add `from typing import Any` after `from pathlib import Path`. Insert before `def extract_allele_reads(`:

```python
def _haploid_filter(
    min_qual: float, settings: CallingSettings
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return length-partitioned VCF filter options and their explicit provenance.

    ``calling.haploid_min_qual`` (default 4.0) applies to haploid, length-partitioned
    calls; ``null`` makes them follow ``run.min_qual``/``--min-qual``. An explicit
    non-default ``--min-qual`` that is overridden is logged, never silently dropped.
    """
    if settings.haploid_min_qual is None:
        qual, source = float(min_qual), "run.min_qual"
    else:
        qual, source = float(settings.haploid_min_qual), "calling.haploid_min_qual"
        if min_qual != DEFAULT_SETTINGS.run.min_qual:
            logger.warning(
                "--min-qual %s is not applied to length-partitioned calls; "
                "calling.haploid_min_qual=%s is used. Set it to null to follow --min-qual.",
                min_qual,
                qual,
            )
    options = {
        "haploid_majority": settings.haploid_majority,
        "haploid_min_qual": qual,
        "haploid_alt_fraction": settings.haploid_alt_fraction,
        "haploid_ref_fraction": settings.haploid_ref_fraction,
    }
    provenance = {
        "min_qual": qual,
        "min_qual_source": source,
        "haploid_majority": settings.haploid_majority,
    }
    return options, provenance
```

Haplotag path in `disambiguate_same_length_alleles`: after `per_allele_threads = max(1, threads // 2)` add `hp_filter, _ = _haploid_filter(min_qual, settings)`. In `_call_hp`'s `filter_vcf(...)`, replace the three lines `haploid_majority=True,`, `haploid_alt_fraction=...` and `haploid_ref_fraction=...` with `**hp_filter,`.

Distinct path in `call_variants_per_allele`: directly before `def _process_allele(`, add:

```python
    filter_options, filter_provenance = _haploid_filter(min_qual, settings)
```

Replace the `# Filter VCF` comment, the `hap_min_qual = getattr(...)` line and the whole `filtered = filter_vcf(...)` call with:

```python
        filtered = filter_vcf(
            vcf, contig_ref, allele_dir, min_qual=min_qual, min_dp=min_dp, **filter_options
        )
```

Before `return allele_key, filtered`, add:

```python
        allele_info["variant_filter"] = filter_provenance
```

`src/muc_one_span/cli.py`: replace both occurrences, at lines 256 and 480, of `    help="Minimum QUAL score for VCF filtering (default 5.0).",` with:

```python
    help="Minimum VCF QUAL (default 5.0); see calling.haploid_min_qual for length-split calls.",
```

This is 96 characters and does not change the line count.

`docs/guides/configuration.md`: in the calling table, before the fraction rows added in Task 2, add:

```markdown
| `calling.haploid_majority` | `true` | Apply the allele-fraction genotype rule to length-partitioned (and haplotagged) calls. With `false`, diploid genotypes are kept, and any heterozygous record then leaves the allele unresolved. |
| `calling.haploid_min_qual` | `4.0` | QUAL threshold for length-partitioned calls; number >=0 or `null`. `null` applies `run.min_qual`/`--min-qual`. A non-default `--min-qual` that this value overrides is logged. The applied value is recorded in `alleles.json` as `variant_filter`. |
```

Directly after the "Classification and consensus" table, before the `### Confidence weights` heading, add:

```markdown
`consensus.haploid_majority` and `consensus.haploid_min_qual` are deprecated and
have no effect. Configurations that set non-default values load with a
`DeprecationWarning`; use the `calling.*` fields.
```

- [ ] **Step 4: Run the tests again**

Run: `uv run --locked --all-extras pytest tests/unit/test_vcf.py tests/unit/test_distinct_calling.py tests/unit/test_runtime_settings.py tests/unit/test_calling.py tests/unit/test_api_defaults.py tests/unit/test_cli.py --no-cov -q`
Expected: all pass. Then `make ci-check` must pass; it checks that `cli.py` stays at 638 lines.

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/calling.py src/muc_one_span/vcf.py src/muc_one_span/settings.py \
  src/muc_one_span/cli.py docs/guides/configuration.md tests/unit/test_vcf.py \
  tests/unit/test_distinct_calling.py tests/unit/test_runtime_settings.py
git commit -m "fix(calling): honour haploid calling settings and record the applied QUAL

--min-qual was silently replaced by calling.haploid_min_qual on length-partitioned
calls, filter_vcf lowered QUAL to 4.0 implicitly, calling.haploid_majority was
never read and consensus.haploid_* had no effect. haploid_min_qual may now be
null (follow --min-qual), overrides are logged, the applied threshold is recorded
as alleles[*].variant_filter, and consensus.haploid_* warn as deprecated.
Default filter commands are unchanged.

Fixes #66

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: D9 + D10 — IGV preflight before mapping; separate cluster and remapped BAMs

**Files:**
- Modify: `src/muc_one_span/report_igv.py` (insert `preflight_igv_report` before `build_igv_context`, line 259)
- Modify: `src/muc_one_span/pipeline.py:90` (tool check and preflight)
- Modify: `src/muc_one_span/calling.py:30-67` (`extract_allele_reads` `output_name`), `calling.py:101-107` (cluster BAM)
- Modify: `tests/unit/test_report_igv.py` (import and append), `tests/unit/test_pipeline_gates.py` (append), `tests/unit/test_calling.py` (append)

**Interfaces:**
- Produces: `report_igv.preflight_igv_report(report_igv: str, work_dir: Path) -> None`. It raises `RuntimeError` naming the mode and the verified version.
- Produces: `calling.extract_allele_reads(bam_path, contig_names, output_dir, *, output_name: str = "allele_reads.bam") -> Path`. The public default is unchanged.
- The pipeline adds `create_report` to `check_tools` only when `run.report_igv != "off"`.

- [ ] **Step 1: Write the failing tests**

In `tests/unit/test_report_igv.py`, add `preflight_igv_report,` to the `from muc_one_span.report_igv import (...)` list, between `create_locus_bed` and `run_igv_report`. Then append:

```python
def test_preflight_is_a_noop_when_igv_is_off(tmp_path: Path) -> None:
    with patch("muc_one_span.report_igv.run_igv_report") as run:
        preflight_igv_report(REPORT_IGV_OFF, tmp_path)
    run.assert_not_called()


def test_preflight_runs_create_report_on_a_synthetic_locus(tmp_path: Path) -> None:
    seen: dict[str, str] = {}

    def fake_run(
        bed: Path, fasta: Path, output: Path, *, report_igv: str, vcf_paths: dict[str, Path]
    ) -> Path:
        seen["fasta"] = fasta.read_text()
        seen["vcf"] = vcf_paths["probe"].read_text()
        seen["mode"] = report_igv
        return output

    with patch("muc_one_span.report_igv.run_igv_report", side_effect=fake_run):
        preflight_igv_report(REPORT_IGV_EMBEDDED, tmp_path)
    assert seen["mode"] == REPORT_IGV_EMBEDDED
    assert seen["fasta"].startswith(">probe\n")
    assert "probe\t150\t.\tC\tT" in seen["vcf"]
    assert list(tmp_path.iterdir()) == []


def test_preflight_failure_is_actionable(tmp_path: Path) -> None:
    error = ValueError("Malformed VCF track 'Probe variants' from create_report")
    with (
        patch("muc_one_span.report_igv.run_igv_report", side_effect=error),
        pytest.raises(RuntimeError, match="IGV report preflight failed for --report-igv sidecar"),
    ):
        preflight_igv_report(REPORT_IGV_SIDECAR, tmp_path)
```

Append to `tests/unit/test_pipeline_gates.py`:

```python
def test_failed_igv_preflight_stops_before_mapping(tmp_path: Path) -> None:
    failing = patch(
        "muc_one_span.report_igv.preflight_igv_report",
        side_effect=RuntimeError("IGV report preflight failed for --report-igv embedded"),
    )
    result, calls = run_mocked_pipeline(
        tmp_path, "--report-igv", "embedded", extra_patches=[failing]
    )
    assert result.exit_code != 0
    assert "map" not in calls
    assert calls == ["check:minimap2,samtools,bcftools,run_clair3.sh,create_report"]
    status = json.loads((tmp_path / "out" / "run_status.json").read_text())
    assert status["status"] == "execution_failed"


def test_igv_off_skips_preflight_and_create_report(tmp_path: Path) -> None:
    probe = patch("muc_one_span.report_igv.preflight_igv_report")
    with probe as preflight:
        result, calls = run_mocked_pipeline(tmp_path)
    assert result.exit_code == 0, result.output
    preflight.assert_not_called()
    assert calls[0] == "check:minimap2,samtools,bcftools,run_clair3.sh"
```

Append to `tests/unit/test_calling.py`:

```python
def test_cluster_and_remapped_bams_use_distinct_paths(tmp_path: Path) -> None:
    """The ladder-cluster subset never occupies allele_reads.bam and is removed."""
    out_dir = tmp_path / "out"

    def fake_run_tool(cmd: list[str]) -> str:
        if cmd[:2] in (["samtools", "view"], ["samtools", "sort"]):
            Path(cmd[cmd.index("-o") + 1]).write_bytes(b"bam")
        if cmd[:2] == ["samtools", "index"]:
            Path(cmd[2] + ".bai").write_bytes(b"bai")
        return ""

    with patch("muc_one_span.calling.run_tool", side_effect=fake_run_tool) as run:
        result = _extract_and_remap_reads(
            tmp_path / "mapping.bam",
            ["contig_50", "contig_51"],
            "contig_51",
            tmp_path / "ref.fa",
            out_dir,
            threads=1,
        )
    calls = [c.args[0] for c in run.call_args_list]
    view = next(c for c in calls if c[:2] == ["samtools", "view"])
    fastq = next(c for c in calls if c[:2] == ["samtools", "fastq"])
    sort = next(c for c in calls if c[:2] == ["samtools", "sort"])
    assert view[view.index("-o") + 1] == str(out_dir / "cluster_reads.bam")
    assert fastq[-1] == str(out_dir / "cluster_reads.bam")
    assert sort[sort.index("-o") + 1] == str(out_dir / "allele_reads.bam") == str(result)
    assert not (out_dir / "cluster_reads.bam").exists()
    assert not (out_dir / "cluster_reads.bam.bai").exists()
    assert (out_dir / "allele_reads.bam.bai").exists()
```

- [ ] **Step 2: Run the tests and confirm they fail**

Run: `uv run --locked --all-extras pytest tests/unit/test_report_igv.py tests/unit/test_pipeline_gates.py tests/unit/test_calling.py --no-cov -q`
Expected:
- `test_report_igv.py` fails at import with `ImportError: preflight_igv_report`.
- `test_failed_igv_preflight_stops_before_mapping` and `test_igv_off_skips_preflight_and_create_report` fail with `AttributeError`, because there is nothing to patch.
- `test_cluster_and_remapped_bams_use_distinct_paths` fails because the view output is `allele_reads.bam`.

- [ ] **Step 3: Implement**

`src/muc_one_span/report_igv.py`: insert before `def build_igv_context(`. The module already imports `tempfile` and `REPORT_IGV_OFF`.

```python
def preflight_igv_report(report_igv: str, work_dir: Path) -> None:
    """Fail before analysis when the requested IGV report cannot be produced.

    Runs ``create_report`` once on a synthetic one-record locus in a temporary
    directory, reusing :func:`run_igv_report` validation. A missing executable or
    a version that corrupts VCF tracks (observed with igv-reports 1.16.x) then
    stops the run before mapping instead of after the analysis has completed.
    """
    if report_igv == REPORT_IGV_OFF:
        return
    work_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="igv-preflight-", dir=work_dir) as tmp:
        root = Path(tmp)
        sequence = "ACGT" * 300
        (root / "probe.fa").write_text(f">probe\n{sequence}\n", encoding="utf-8")
        (root / "probe.bed").write_text("probe\t100\t200\tprobe\n", encoding="utf-8")
        (root / "probe.vcf").write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=probe,length=1200>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            f"probe\t150\t.\t{sequence[149]}\tT\t30\tPASS\t.\n",
            encoding="utf-8",
        )
        try:
            run_igv_report(
                root / "probe.bed",
                root / "probe.fa",
                root / "probe.html",
                report_igv=report_igv,
                vcf_paths={"probe": root / "probe.vcf"},
            )
        except (OSError, RuntimeError, ValueError) as exc:
            raise RuntimeError(
                f"IGV report preflight failed for --report-igv {report_igv}: {exc}. "
                "Use a working igv-reports create_report (1.13.0 verified) or --report-igv off."
            ) from exc
```

`src/muc_one_span/pipeline.py`: replace `check_tools(["minimap2", "samtools", "bcftools", "run_clair3.sh"])` with:

```python
    igv_requested = settings.run.report_igv != "off"
    check_tools(
        ["minimap2", "samtools", "bcftools", "run_clair3.sh"]
        + (["create_report"] if igv_requested else [])
    )
    if igv_requested:
        from muc_one_span.report_igv import preflight_igv_report

        preflight_igv_report(settings.run.report_igv, out)
```

`src/muc_one_span/calling.py`, `extract_allele_reads`:

- add `*, output_name: str = "allele_reads.bam",` after `output_dir: Path,`;
- add the docstring arg `output_name: File name of the extracted BAM inside ``output_dir``.`;
- replace `out_bam = output_dir / "allele_reads.bam"` with `out_bam = output_dir / output_name`.

In `_extract_and_remap_reads`, replace steps 1 and 2 with:

```python
    # 1. Extract reads from all cluster contigs (ladder coordinates; intermediate)
    cluster_bam = extract_allele_reads(
        bam_path, cluster_contigs, output_dir, output_name="cluster_reads.bam"
    )

    # 2. Convert to FASTQ, then drop the cluster BAM so allele_reads.bam is only
    #    ever the remapped single-contig BAM.
    fastq_path = output_dir / "cluster_reads.fq"
    stdout = run_tool(["samtools", "fastq", str(cluster_bam)])
    fastq_path.write_text(stdout)
    cluster_bam.unlink(missing_ok=True)
    Path(f"{cluster_bam}.bai").unlink(missing_ok=True)
```

- [ ] **Step 4: Run the tests again, including the real preflight**

Run: `uv run --locked --all-extras pytest tests/unit --no-cov -q`. Expected: all pass.
Run: `make test-int` with the tool environment active and igv-reports 1.13 first on `PATH`. Expected: `test_report_sessions.py` and `test_calling.py::TestExtractAlleleReads` pass. Report every skip and every environment failure, including the Clair3 model test if no model is installed.
Manual real-tool check. Run this once with a 1.13 `create_report` first on `PATH`, where it should print `ok`, and once with a 1.16.x install, where it should raise `RuntimeError: IGV report preflight failed ...`:

```bash
uv run --locked --all-extras python -c "from pathlib import Path; import tempfile; \
from muc_one_span.report_igv import preflight_igv_report; \
preflight_igv_report('embedded', Path(tempfile.mkdtemp())); print('ok')"
```

Then `make ci-check` must pass.

- [ ] **Step 5: Commit**

```bash
git add src/muc_one_span/report_igv.py src/muc_one_span/pipeline.py src/muc_one_span/calling.py \
  tests/unit/test_report_igv.py tests/unit/test_pipeline_gates.py tests/unit/test_calling.py
git commit -m "fix(pipeline): preflight IGV reports and separate cluster/remapped BAM paths

A requested IGV report was validated only after the full analysis; create_report
is now required up front and a synthetic one-record preflight rejects missing or
track-corrupting igv-reports before mapping (the rejection itself is unchanged).
The ladder-cluster subset BAM no longer occupies allele_reads.bam and is removed
after FASTQ conversion.

Fixes #67
Fixes #68

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Validation on PRJEB92208, CHANGELOG and version 0.16.0

**Files:**
- Modify: `CHANGELOG.md` (`[Unreleased]` becomes `[0.16.0]`)
- Modify: `src/muc_one_span/version.py:3` (`0.15.1` becomes `0.16.0`). This is the repository convention: `pyproject.toml` takes the version from here through `[tool.hatch.version]`.
- Modify: `docs/guides/clinical-validation-results.md` (new subsection with the measured table)
- Modify: `.planning/2026-09-23-hybrid-engine-plan.md` (Task 1 status line)
- Outside Git: `../MucOneSpan-wave2-data/environment-p0.json`, `../MucOneSpan-wave2-data/cohort-p0/`, `../MucOneSpan-wave2-data/results-p0.json`

**Interfaces:**
- Consumes: `scripts/clinical_benchmark.py freeze|run|score` (unchanged), `benchmarks/clinical/prjeb92208/manifest.json`, `report.compute_clinical_decision` and `selection_qc.annotate_selection_qc`.

- [ ] **Step 1: Run the full check set**

```bash
make ci-check
make docs-check
make build-check
conda activate env_clair3            # the tool env used for cohort-v3 (bcftools 1.17, Clair3 1.0.10)
export PATH="$PWD/../MucOneSpan-wave2-data/igv-reports-1.13/bin:$PATH"
make test-int
```

Expected: everything passes. Record the pass, fail and skip counts in the PR. Record which integration tests could not run and why, for example a missing Clair3 model or missing generated reads.

- [ ] **Step 2: MP4 smoke run (single library)**

```bash
DATA=../MucOneSpan-wave2-data
export CUDA_VISIBLE_DEVICES=''
uv run --locked --all-extras python scripts/clinical_benchmark.py freeze \
  --checkout . --model "$DATA/models/r1041_e82_400bps_sup_v500" --output "$DATA/environment-p0.json"
uv run --locked --all-extras python scripts/clinical_benchmark.py run \
  --manifest benchmarks/clinical/prjeb92208/manifest.json --data-root "$DATA/reads" \
  --output-root "$DATA/cohort-p0-smoke" --environment "$DATA/environment-p0.json" \
  --threads 2 --timeout 1800 --run ERR15277569
```

Expected results in `$DATA/cohort-p0-smoke/ERR15277569`:
- `allele_2/variants.vcf.gz` has `contig_71:3432 G>GC` with GT `1/1`;
- `alleles.allele_2.allele_genotype_status` is `allele_specific_resolved`;
- `classifications.allele_2.mutations` contains a dupC at repeat 49 with `vcf_support_status=exact_sequence_concordance`.

If any of these fail, stop and investigate. Do not tune thresholds.

- [ ] **Step 3: Full cohort rerun and decision table**

```bash
uv run --locked --all-extras python scripts/clinical_benchmark.py run \
  --manifest benchmarks/clinical/prjeb92208/manifest.json --data-root "$DATA/reads" \
  --output-root "$DATA/cohort-p0" --environment "$DATA/environment-p0.json" \
  --threads 2 --timeout 1800
uv run --locked --all-extras python scripts/clinical_benchmark.py score \
  --manifest benchmarks/clinical/prjeb92208/manifest.json --truth "$DATA/truth.json" \
  --output-root "$DATA/cohort-p0" --output "$DATA/results-p0.json"
uv run --locked --all-extras python - "$DATA" <<'EOF'
import json
import sys
from pathlib import Path

from muc_one_span.report import compute_clinical_decision

data = Path(sys.argv[1])
names = {"52": "HG002 WGS", "53": "WGS (identity unresolved)", "62": "HG001", "63": "HG002 PCR",
         "64": "HG003", "65": "HG004", "66": "MP1", "67": "MP2", "68": "MP3", "69": "MP4",
         "70": "MP5"}
print("| Library | Sample | Decision | Events (allele:name@repeat:support) | Gate reasons |")
print("| --- | --- | --- | --- | --- |")
for suffix, name in names.items():
    summary = json.loads((data / "cohort-p0" / f"ERR152775{suffix}" / "summary.json").read_text())
    decision = compute_clinical_decision(summary)
    events = "; ".join(
        f"{key}:{m.get('mutation_name') or m.get('closest_type')}@{m.get('repeat_index')}:"
        f"{m.get('vcf_support_status')}"
        for key, c in summary["classifications"].items() for m in c.get("mutations", [])
    )
    gates = "; ".join(d for d in decision["details"] if "Allele" in d and ":" in d)[:300]
    print(f"| ERR152775{suffix} | {name} | {decision['state']} | {events or '-'} | {gates or '-'} |")
EOF
```

Acceptance gate, to be met before any release step:

| Library | v0.15.1 (cohort-v3) | Expected after P0 | Must hold |
| --- | --- | --- | --- |
| ERR15277569 MP4 | NO_PATHOGENIC | **PATHOGENIC**, dupC r49, allele 2 | yes |
| ERR15277566 MP1 / 67 MP2 | PATHOGENIC dupC r17 | PATHOGENIC. MP1 adds a caveat: allele 1 length 39 vs contig 44. | stays PATHOGENIC |
| ERR15277568 MP3 | INCONCLUSIVE | INCONCLUSIVE: dupC r7 `heterozygous_genotype_unresolved`, selection unresolved (0.62) | not NEGATIVE |
| ERR15277562 HG001 / 64 HG003 | NO_PATHOGENIC | NO_PATHOGENIC | not PATHOGENIC |
| ERR15277563 HG002 PCR / 65 HG004 / 52 HG002 WGS | INCONCLUSIVE | INCONCLUSIVE | not PATHOGENIC |
| ERR15277570 MP5 | PATHOGENIC (untemplated frameshift at repeat 35) | INCONCLUSIVE: identity not established, 0.28 fragment mode, length 69 vs 61 | not PATHOGENIC. **Owner must accept this intended change.** |
| ERR15277553 (identity unresolved) | NO_PATHOGENIC | INCONCLUSIVE (low per-allele depth, unselected cluster) | no truth; report only |

The expected column for the unchanged libraries was previewed by applying the Task 5 gates to the cohort-v3 summaries. MP4 needs the rerun, because its VCF must pass through Tasks 1–2.

In `results-p0.json`, `all_positive_recovery` must rise by exactly MP4 compared with `$DATA/results-v3.json`, with no new event in any control. Any other deviation must be explained in the PR before release; do not adjust thresholds to hide it.

- [ ] **Step 4: Document the measured result**

In `docs/guides/clinical-validation-results.md`, add a subsection `## v0.16.0 clinical-safety rerun (PRJEB92208)` before `## Reproduction`. It contains:
- the table printed in Step 3, pasted verbatim;
- one sentence stating that positives are publication-reported (three families) and that no sensitivity or specificity is claimed;
- the tool versions from `environment-p0.json`.

Add no patient sequences.

- [ ] **Step 5: CHANGELOG, version and planning status**

Set `__version__ = "0.16.0"` in `src/muc_one_span/version.py`. Replace `## [Unreleased]` in `CHANGELOG.md` with the block below, filling in the release date (`date +%F`) and the issue numbers from Task 0:

```markdown
## [Unreleased]

## [0.16.0] - RELEASE_DATE

### Fixed

- Distinct-length candidates no longer resolve a heterozygous genotype to the
  reference allele. Remaining heterozygous, phased or conflicting records use the
  IUPAC candidate with `allele_genotype_status` and no independent haplotype
  credit. This fixes the PRJEB92208 MP4 dupC false negative (#53).
- Haploid genotypes on length-partitioned alleles use the allele-specific AD
  fraction instead of FORMAT/AF. The 0.5/0.2 cut-offs are now configurable (#65).
- VCF concordance mirrors `bcftools consensus -H I` for heterozygous indels and
  is judged per event. The new status `heterozygous_genotype_unresolved` is
  never support (#64).
- Clinical decisions require explicit support, a frameshift, exact template
  identity, resolved localization and adequate carrier-allele depth for
  PATHOGENIC. NEGATIVE additionally requires resolved allele selection,
  consistent length and adequate per-allele primary depth (#55, #63).
- `--min-qual`, `calling.haploid_majority` and `calling.haploid_min_qual` are
  honoured and the applied threshold is recorded. `consensus.haploid_*` are
  deprecated no-ops (#66).
- A requested IGV report is preflighted before mapping (#67). Cluster and
  remapped BAMs no longer share `allele_reads.bam` (#68).

### Added

- Settings `calling.haploid_alt_fraction`, `calling.haploid_ref_fraction`,
  `allele_selection.secondary_mode_min_fraction` and
  `allele_selection.min_allele_primary_records`.
- Additive output fields:
  - `alleles[*]`: `allele_genotype_status`, `heterozygous_sites`,
    `variant_filter`, `selection_status`, `secondary_mode_fraction`,
    `depth_status`, `depth_basis`, `depth_threshold`, `length_status`;
  - `vcf_projection.unresolved_genotype_edits`.

### Changed

- Untemplated (novel) frameshifts, and events on alleles below the per-allele
  depth gate, are reported INCONCLUSIVE with reasons instead of PATHOGENIC.
  Summaries without `vcf_support` no longer count as supported.
```

In `.planning/2026-09-23-hybrid-engine-plan.md`, directly under the `### Task 1: Explicit-support clinical decision (safety prerequisite)` heading, add the line: `Status: implemented in P0 (v0.16.0, .planning/2026-09-24-p0-clinical-safety-plan.md Task 4); helper lives in clinical_gates.mutation_supported.`

Run `make ci-check` and `make docs-check`. Expected: pass.

- [ ] **Step 6: Commit and hand off for approval**

```bash
git add CHANGELOG.md src/muc_one_span/version.py docs/guides/clinical-validation-results.md \
  .planning/2026-09-23-hybrid-engine-plan.md
git commit -m "chore(release): bump version to 0.16.0 and document changes

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

Push the branch and open a PR only after the owner agrees. The PR body must include:
- the Step 1 check counts;
- the Step 3 table;
- the MP5 decision change for explicit acceptance;
- any integration skips.

End the PR body with the line `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.

Do not merge, tag or create a release without explicit owner approval. After approval: merge, `git tag v0.16.0`, then create the GitHub release from the CHANGELOG section. Once merged, move this plan to `.planning/archive/`.

---

## Self-review

- **Coverage.** Every Phase 0 row of the execution-order table maps to a task:
  - 0.1 is Task 1;
  - 0.2 is Task 2;
  - 0.3 is Task 3;
  - 0.4 is Task 4;
  - 0.5 is Task 5;
  - 0.6 is Tasks 5 (D3), 6 (D7) and 7 (D9, D10).
- **#53 acceptance.**
  - It never picks allele 2 blindly.
  - Unresolved evidence is explicit (`allele_genotype_status`, `heterozygous_sites`, the IUPAC selector) and gets no independent credit.
  - It covers homozygous ALT, phased, single-site and multi-site unphased cases, plus conflicting records.
  - MP4 recovery comes from AD evidence (Task 2). The false-positive check is Task 8.
- **Types.**
  - `length_partition_selection` returns `tuple[int | str, str, list[dict]]`, and `annotate_consensus_candidate` accepts `str | int`.
  - `_haploid_filter` returns kwargs that match the `filter_vcf` keywords added in Task 2.
  - `depth_status` values match the hybrid plan's `adequate`/`low`.
- **Open risks.**
  1. The secondary-mode and depth gates are calibrated on 9 clinical libraries and 80 simulated runs. Phase 1 (`dev` split) must measure the INCONCLUSIVE rate on realistic smear.
  2. The bcftools `-H I` indel behaviour was verified on 1.17 only; the conda environment pins 1.21. The Task 3 integration guard must run on 1.21 before release.
  3. The MP5 downgrade from PATHOGENIC to INCONCLUSIVE is intended by #55. It could hide a true positive, and the owner must decide.
  4. ONT homopolymer stutter can still leave a real dupC in the 0.2–0.5 band, as with MP3. P0 makes that INCONCLUSIVE, not PATHOGENIC. The P1 read-level homopolymer channel (#59) resolves it.
