# Platform-Aware Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make Clair3 `--platform` and minimap2 `-x` preset configurable so ONT data can be analyzed alongside HiFi.

**Architecture:** Thread two new parameters (`platform` for Clair3, `preset` for minimap2) from CLI options through `mapping.py` and `calling.py`. A `PLATFORM_PRESETS` dict auto-maps platform to minimap2 preset, with explicit override via `--minimap2-preset`. All defaults remain `hifi`/`map-hifi` so existing behavior is unchanged.

**Tech Stack:** Python 3.10+, Click 8.0+, pytest, unittest.mock

**Design doc:** `.planning/2026-04-07-platform-support-design.md`

---

### Task 1: Add `preset` parameter to `mapping.py`

**Files:**
- Modify: `src/muc_one_span/mapping.py:77-108` (`_run_mapping_pipeline`)
- Modify: `src/muc_one_span/mapping.py:34-74` (`map_reads`)
- Test: `tests/unit/test_mapping.py`

- [ ] **Step 1: Write failing test for `_run_mapping_pipeline` preset parameter**

Add to `tests/unit/test_mapping.py` inside `TestRunMappingPipeline`:

```python
def test_preset_passed_to_minimap2(self, mocker):
    """_run_mapping_pipeline passes the preset to minimap2 -x flag."""
    from muc_one_span.mapping import _run_mapping_pipeline

    mock_p1 = mocker.MagicMock()
    mock_p1.stdout = mocker.MagicMock()
    mock_p1.returncode = 0
    mock_p1.stderr = mocker.MagicMock()
    mock_p1.stderr.read.return_value = b""

    mock_p2 = mocker.MagicMock()
    mock_p2.returncode = 0
    mock_p2.communicate.return_value = (b"", b"")

    mock_popen = mocker.patch(
        "muc_one_span.mapping.subprocess.Popen",
        side_effect=[mock_p1, mock_p2],
    )

    _run_mapping_pipeline(
        input_path=Path("/tmp/test.fastq"),
        reference_path=Path("/tmp/ref.fa"),
        bam_path=Path("/tmp/out.bam"),
        threads=1,
        preset="lr:hq",
    )

    minimap2_cmd = mock_popen.call_args_list[0][0][0]
    assert "-x" in minimap2_cmd
    x_idx = minimap2_cmd.index("-x")
    assert minimap2_cmd[x_idx + 1] == "lr:hq"

def test_preset_defaults_to_map_hifi(self, mocker):
    """_run_mapping_pipeline defaults to map-hifi when no preset given."""
    from muc_one_span.mapping import _run_mapping_pipeline

    mock_p1 = mocker.MagicMock()
    mock_p1.stdout = mocker.MagicMock()
    mock_p1.returncode = 0
    mock_p1.stderr = mocker.MagicMock()
    mock_p1.stderr.read.return_value = b""

    mock_p2 = mocker.MagicMock()
    mock_p2.returncode = 0
    mock_p2.communicate.return_value = (b"", b"")

    mock_popen = mocker.patch(
        "muc_one_span.mapping.subprocess.Popen",
        side_effect=[mock_p1, mock_p2],
    )

    _run_mapping_pipeline(
        input_path=Path("/tmp/test.fastq"),
        reference_path=Path("/tmp/ref.fa"),
        bam_path=Path("/tmp/out.bam"),
        threads=1,
    )

    minimap2_cmd = mock_popen.call_args_list[0][0][0]
    x_idx = minimap2_cmd.index("-x")
    assert minimap2_cmd[x_idx + 1] == "map-hifi"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_mapping.py::TestRunMappingPipeline::test_preset_passed_to_minimap2 tests/unit/test_mapping.py::TestRunMappingPipeline::test_preset_defaults_to_map_hifi -v --no-cov`

Expected: FAIL — `_run_mapping_pipeline() got an unexpected keyword argument 'preset'`

- [ ] **Step 3: Write failing test for `map_reads` preset parameter**

Add to `tests/unit/test_mapping.py` inside `TestMapReads`:

```python
@patch("muc_one_span.mapping._run_mapping_pipeline")
@patch("muc_one_span.mapping.run_tool")
def test_preset_passed_to_pipeline(self, mock_run_tool, mock_pipeline, tmp_path):
    """The preset parameter is forwarded to _run_mapping_pipeline."""
    mock_run_tool.return_value = ""
    fastq = tmp_path / "reads.fq"
    fastq.touch()
    ref = tmp_path / "ref.fa"
    ref.touch()

    map_reads(fastq, ref, tmp_path, threads=4, preset="lr:hq")

    pipeline_call = mock_pipeline.call_args
    assert pipeline_call[1].get("preset") == "lr:hq" or (
        len(pipeline_call[0]) >= 5 and pipeline_call[0][4] == "lr:hq"
    )
```

- [ ] **Step 4: Run test to verify it fails**

Run: `uv run pytest tests/unit/test_mapping.py::TestMapReads::test_preset_passed_to_pipeline -v --no-cov`

Expected: FAIL — `map_reads() got an unexpected keyword argument 'preset'`

- [ ] **Step 5: Implement `preset` parameter in `_run_mapping_pipeline`**

In `src/muc_one_span/mapping.py`, change the `_run_mapping_pipeline` signature and body:

```python
def _run_mapping_pipeline(
    input_path: Path,
    reference_path: Path,
    bam_path: Path,
    threads: int,
    preset: str = "map-hifi",
) -> None:
```

And replace the hardcoded `"map-hifi"` on line 103 with `preset`:

```python
    minimap2_cmd = [
        "minimap2",
        "-a",
        "-x",
        preset,
        "-t",
        str(threads),
        str(reference_path),
        str(input_path),
    ]
```

- [ ] **Step 6: Implement `preset` parameter in `map_reads`**

In `src/muc_one_span/mapping.py`, change the `map_reads` signature:

```python
def map_reads(
    input_path: Path,
    reference_path: Path,
    output_dir: Path,
    threads: int = 4,
    preset: str = "map-hifi",
) -> Path:
```

And update the `_run_mapping_pipeline` call on line 69:

```python
    _run_mapping_pipeline(actual_input, reference_path, bam_path, threads, preset)
```

Update the docstring to mention the preset parameter:

```
    Args:
        input_path: Path to input FASTQ or BAM file.
        reference_path: Path to reference FASTA.
        output_dir: Directory for output files.
        threads: Number of threads for minimap2/samtools (default 4).
        preset: minimap2 preset (default ``"map-hifi"``). Use ``"lr:hq"`` for ONT.
```

- [ ] **Step 7: Run all mapping tests**

Run: `uv run pytest tests/unit/test_mapping.py -v --no-cov`

Expected: All PASS

- [ ] **Step 8: Commit**

```bash
git add src/muc_one_span/mapping.py tests/unit/test_mapping.py
git commit -m "feat: add preset parameter to mapping.py for configurable minimap2 preset"
```

---

### Task 2: Add `platform` and `preset` parameters to `calling.py`

**Files:**
- Modify: `src/muc_one_span/calling.py:65-150` (`_extract_and_remap_reads`)
- Modify: `src/muc_one_span/calling.py:197-302` (`disambiguate_same_length_alleles`)
- Modify: `src/muc_one_span/calling.py:305-413` (`call_variants_per_allele`)
- Test: `tests/unit/test_calling.py`

- [ ] **Step 1: Write failing test for `_extract_and_remap_reads` preset**

Add to `tests/unit/test_calling.py` inside `TestExtractAndRemapReads`:

```python
@patch("muc_one_span.calling.run_tool")
def test_preset_passed_to_minimap2(self, mock_run_tool, tmp_path):
    """_extract_and_remap_reads passes preset to the minimap2 command."""
    mock_run_tool.return_value = ""
    bam = tmp_path / "mapping.bam"
    bam.touch()
    ref = tmp_path / "ref.fa"
    ref.touch()

    _extract_and_remap_reads(
        bam,
        ["contig_51"],
        "contig_51",
        ref,
        tmp_path / "out",
        threads=1,
        preset="lr:hq",
    )

    minimap2_calls = [
        c[0][0] for c in mock_run_tool.call_args_list if c[0][0][0] == "minimap2"
    ]
    assert minimap2_calls, "No minimap2 call found"
    minimap2_cmd = minimap2_calls[0]
    x_idx = minimap2_cmd.index("-x")
    assert minimap2_cmd[x_idx + 1] == "lr:hq"

@patch("muc_one_span.calling.run_tool")
def test_preset_defaults_to_map_hifi(self, mock_run_tool, tmp_path):
    """_extract_and_remap_reads defaults to map-hifi preset."""
    mock_run_tool.return_value = ""
    bam = tmp_path / "mapping.bam"
    bam.touch()
    ref = tmp_path / "ref.fa"
    ref.touch()

    _extract_and_remap_reads(
        bam,
        ["contig_51"],
        "contig_51",
        ref,
        tmp_path / "out",
        threads=1,
    )

    minimap2_calls = [
        c[0][0] for c in mock_run_tool.call_args_list if c[0][0][0] == "minimap2"
    ]
    minimap2_cmd = minimap2_calls[0]
    x_idx = minimap2_cmd.index("-x")
    assert minimap2_cmd[x_idx + 1] == "map-hifi"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_calling.py::TestExtractAndRemapReads::test_preset_passed_to_minimap2 tests/unit/test_calling.py::TestExtractAndRemapReads::test_preset_defaults_to_map_hifi -v --no-cov`

Expected: FAIL — `_extract_and_remap_reads() got an unexpected keyword argument 'preset'`

- [ ] **Step 3: Write failing test for `call_variants_per_allele` platform and preset**

Add to `tests/unit/test_calling.py` inside `TestCallVariantsPerAllele`:

```python
@patch("muc_one_span.vcf.run_tool")
@patch("muc_one_span.calling.run_tool")
def test_platform_and_preset_threaded_through(self, mock_run_tool, mock_vcf_tool, tmp_path):
    """platform is passed to run_clair3, preset is passed to _extract_and_remap_reads."""
    mock_run_tool.return_value = ""
    mock_vcf_tool.return_value = ""
    bam = tmp_path / "mapping.bam"
    bam.touch()
    ref = tmp_path / "ref.fa"
    ref.touch()
    alleles = self._make_alleles_homozygous()

    call_variants_per_allele(
        bam, ref, alleles, tmp_path,
        platform="ont",
        preset="lr:hq",
    )

    # Check minimap2 was called with lr:hq
    minimap2_calls = [
        c[0][0] for c in mock_run_tool.call_args_list if c[0][0][0] == "minimap2"
    ]
    assert minimap2_calls
    minimap2_cmd = minimap2_calls[0]
    x_idx = minimap2_cmd.index("-x")
    assert minimap2_cmd[x_idx + 1] == "lr:hq"

    # Check run_clair3.sh was called with --platform=ont
    clair3_calls = [
        c[0][0] for c in mock_run_tool.call_args_list if c[0][0][0] == "run_clair3.sh"
    ]
    assert clair3_calls
    assert any("--platform=ont" in arg for arg in clair3_calls[0])
```

- [ ] **Step 4: Write failing test for `disambiguate_same_length_alleles` platform and preset**

Add to `tests/unit/test_calling.py` inside `TestDisambiguateSameLengthAlleles`:

```python
@patch("muc_one_span.vcf.run_tool", return_value="")
@patch("muc_one_span.calling.run_tool", return_value="")
@patch("muc_one_span.calling.parse_vcf_genotypes", return_value=[])
def test_platform_and_preset_threaded_through(
    self, mock_geno, mock_run, mock_vcf_tool, tmp_path
):
    """platform is passed to run_clair3, preset to _extract_and_remap_reads."""
    alleles = {
        "allele_1": {"contig_name": "contig_51", "cluster_contigs": ["contig_51"]},
        "allele_2": {"contig_name": "contig_51", "cluster_contigs": ["contig_51"]},
    }
    disambiguate_same_length_alleles(
        tmp_path / "bam", tmp_path / "ref.fa", alleles, tmp_path,
        platform="ont",
        preset="lr:hq",
    )

    # Check minimap2 was called with lr:hq
    minimap2_calls = [
        c[0][0] for c in mock_run.call_args_list if c[0][0][0] == "minimap2"
    ]
    assert minimap2_calls
    minimap2_cmd = minimap2_calls[0]
    x_idx = minimap2_cmd.index("-x")
    assert minimap2_cmd[x_idx + 1] == "lr:hq"

    # Check run_clair3.sh was called with --platform=ont
    clair3_calls = [
        c[0][0] for c in mock_run.call_args_list if c[0][0][0] == "run_clair3.sh"
    ]
    assert clair3_calls
    assert any("--platform=ont" in arg for arg in clair3_calls[0])
```

- [ ] **Step 5: Run all new calling tests to verify they fail**

Run: `uv run pytest tests/unit/test_calling.py::TestExtractAndRemapReads::test_preset_passed_to_minimap2 tests/unit/test_calling.py::TestCallVariantsPerAllele::test_platform_and_preset_threaded_through tests/unit/test_calling.py::TestDisambiguateSameLengthAlleles::test_platform_and_preset_threaded_through -v --no-cov`

Expected: FAIL

- [ ] **Step 6: Implement `preset` parameter in `_extract_and_remap_reads`**

In `src/muc_one_span/calling.py`, change the signature:

```python
def _extract_and_remap_reads(
    bam_path: Path,
    cluster_contigs: list[str],
    peak_contig: str,
    reference_path: Path,
    output_dir: Path,
    threads: int = 4,
    preset: str = "map-hifi",
) -> Path:
```

And replace the hardcoded `"map-hifi"` on line 122 with `preset`:

```python
    sam_output = run_tool(
        [
            "minimap2",
            "-a",
            "-x",
            preset,
            "-t",
            str(threads),
            str(contig_ref),
            str(fastq_path),
        ]
    )
```

Add to docstring:

```
        preset: minimap2 preset (default ``"map-hifi"``). Use ``"lr:hq"`` for ONT.
```

- [ ] **Step 7: Implement `platform` and `preset` in `disambiguate_same_length_alleles`**

Change the signature:

```python
def disambiguate_same_length_alleles(
    bam_path: Path,
    reference_path: Path,
    alleles: dict,
    output_dir: Path,
    clair3_model: str = "",
    threads: int = 4,
    min_qual: float = 5.0,
    min_dp: int = 5,
    platform: str = "hifi",
    preset: str = "map-hifi",
) -> dict:
```

Pass `preset` to `_extract_and_remap_reads` (around line 222):

```python
    merged_bam = _extract_and_remap_reads(
        bam_path,
        cluster_contigs,
        contig_name,
        reference_path,
        merged_dir,
        threads,
        preset=preset,
    )
```

Pass `platform` to `run_clair3` (around line 234):

```python
    raw_vcf = run_clair3(
        merged_bam,
        contig_ref,
        clair3_dir,
        model_path=clair3_model,
        platform=platform,
        threads=threads,
    )
```

- [ ] **Step 8: Implement `platform` and `preset` in `call_variants_per_allele`**

Change the signature:

```python
def call_variants_per_allele(
    bam_path: Path,
    reference_path: Path,
    alleles: dict,
    output_dir: Path,
    clair3_model: str = "",
    threads: int = 4,
    min_qual: float = 5.0,
    min_dp: int = 5,
    platform: str = "hifi",
    preset: str = "map-hifi",
) -> dict[str, Path]:
```

Pass through to `disambiguate_same_length_alleles` (around line 344):

```python
        disambig = disambiguate_same_length_alleles(
            bam_path,
            reference_path,
            alleles,
            output_dir,
            clair3_model,
            threads,
            min_qual,
            min_dp,
            platform=platform,
            preset=preset,
        )
```

Inside `_process_allele`, pass `preset` to `_extract_and_remap_reads` and `platform` to `run_clair3`:

```python
        allele_bam = _extract_and_remap_reads(
            bam_path,
            cluster_contigs,
            contig_name,
            reference_path,
            allele_dir,
            threads,
            preset=preset,
        )
```

```python
        vcf = run_clair3(
            allele_bam,
            contig_ref,
            clair3_dir,
            model_path=clair3_model,
            platform=platform,
            threads=per_allele_threads,
        )
```

- [ ] **Step 9: Run all calling tests**

Run: `uv run pytest tests/unit/test_calling.py -v --no-cov`

Expected: All PASS

- [ ] **Step 10: Commit**

```bash
git add src/muc_one_span/calling.py tests/unit/test_calling.py
git commit -m "feat: thread platform and preset through calling.py for ONT support"
```

---

### Task 3: Add CLI options and preset auto-selection

**Files:**
- Modify: `src/muc_one_span/cli.py`
- Test: `tests/unit/test_cli.py`

- [ ] **Step 1: Write failing tests for new CLI options**

Add to `tests/unit/test_cli.py`:

```python
class TestPlatformOptions:
    """Tests for --platform and --minimap2-preset CLI options."""

    def test_run_help_shows_platform(self):
        """run --help includes --platform option."""
        runner = CliRunner()
        result = runner.invoke(main, ["run", "--help"])
        assert "--platform" in result.output

    def test_run_help_shows_minimap2_preset(self):
        """run --help includes --minimap2-preset option."""
        runner = CliRunner()
        result = runner.invoke(main, ["run", "--help"])
        assert "--minimap2-preset" in result.output

    def test_call_help_shows_platform(self):
        """call --help includes --platform option."""
        runner = CliRunner()
        result = runner.invoke(main, ["call", "--help"])
        assert "--platform" in result.output

    def test_call_help_shows_minimap2_preset(self):
        """call --help includes --minimap2-preset option."""
        runner = CliRunner()
        result = runner.invoke(main, ["call", "--help"])
        assert "--minimap2-preset" in result.output

    def test_map_help_shows_minimap2_preset(self):
        """map --help includes --minimap2-preset option."""
        runner = CliRunner()
        result = runner.invoke(main, ["map", "--help"])
        assert "--minimap2-preset" in result.output

    def test_run_rejects_invalid_platform(self):
        """run --platform=invalid is rejected by Click."""
        runner = CliRunner()
        result = runner.invoke(main, ["run", "--platform", "invalid", "--help"])
        assert result.exit_code != 0

    def test_run_platform_ont_accepted(self, tmp_path):
        """run --platform ont is accepted."""
        input_file = tmp_path / "reads.fastq"
        input_file.touch()

        with (
            patch("muc_one_span.tools.check_tools"),
            patch("muc_one_span.tools.get_tool_versions", return_value={}),
            patch("muc_one_span.mapping.map_reads", return_value=tmp_path / "m.bam") as mock_map,
            patch("muc_one_span.mapping.get_idxstats", return_value="c51\t3060\t100\t0\n"),
            patch("muc_one_span.alleles.parse_idxstats", return_value={51: 100}),
            patch("muc_one_span.alleles.detect_alleles", return_value={
                "homozygous": True,
                "allele_1": {"length": 60, "reads": 100, "canonical_repeats": 51,
                             "contig_name": "c51", "cluster_contigs": ["c51"]},
                "allele_2": {"length": 60, "reads": 0, "canonical_repeats": 51,
                             "contig_name": "c51", "cluster_contigs": ["c51"]},
            }),
            patch("muc_one_span.calling.call_variants_per_allele",
                  return_value={"allele_1": tmp_path / "a.vcf.gz"}) as mock_call,
            patch("muc_one_span.consensus.build_consensus_per_allele",
                  return_value={"allele_1": tmp_path / "a.fa"}),
            patch("muc_one_span.classify.classify_sequence",
                  return_value={"structure": "X", "repeats": [], "mutations_detected": [],
                                "allele_confidence": 1.0}),
            patch("muc_one_span.classify.validate_mutations_against_vcf",
                  return_value={"structure": "X", "repeats": [], "mutations_detected": [],
                                "allele_confidence": 1.0}),
            patch("muc_one_span.vcf.parse_vcf_variants", return_value=[]),
        ):
            # Write a fake consensus FASTA for classify to read
            fa = tmp_path / "a.fa"
            fa.write_text(">a\nACGT\n")

            runner = CliRunner()
            result = runner.invoke(
                main,
                ["run", "--input", str(input_file), "--output-dir", str(tmp_path),
                 "--platform", "ont"],
            )

            assert result.exit_code == 0, result.output
            # Verify platform=ont was passed to call_variants_per_allele
            call_kwargs = mock_call.call_args
            assert call_kwargs[1].get("platform") == "ont" or "ont" in str(call_kwargs)
            # Verify preset=lr:hq was auto-selected and passed to map_reads
            map_kwargs = mock_map.call_args
            assert "lr:hq" in str(map_kwargs)

    def test_map_minimap2_preset_override(self, tmp_path):
        """map --minimap2-preset overrides default preset."""
        fastq = tmp_path / "reads.fq"
        fastq.touch()
        ref = tmp_path / "ref.fa"
        ref.touch()

        with (
            patch("muc_one_span.tools.check_tools"),
            patch("muc_one_span.mapping._run_mapping_pipeline"),
            patch("muc_one_span.mapping.run_tool", return_value=""),
            patch("muc_one_span.mapping.map_reads", return_value=tmp_path / "m.bam") as mock_map,
        ):
            runner = CliRunner()
            result = runner.invoke(
                main,
                ["map", "--input", str(fastq), "--reference", str(ref),
                 "--minimap2-preset", "map-ont"],
            )

            assert result.exit_code == 0, result.output
            assert "map-ont" in str(mock_map.call_args)

    def test_main_help_mentions_ont(self):
        """Main --help mentions ONT in description."""
        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert "ONT" in result.output or "ont" in result.output.lower()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/unit/test_cli.py::TestPlatformOptions -v --no-cov`

Expected: FAIL — `--platform` not found in help output

- [ ] **Step 3: Add `PLATFORM_PRESETS` constant and update group docstring**

In `src/muc_one_span/cli.py`, add after imports (around line 12):

```python
PLATFORM_PRESETS: dict[str, str] = {"hifi": "map-hifi", "ont": "lr:hq"}
```

Update the `main` group docstring (line 23):

```python
    """muconespan: MUC1 VNTR analysis pipeline for PacBio HiFi and ONT amplicon data."""
```

- [ ] **Step 4: Add `--platform` and `--minimap2-preset` to `run` subcommand**

Add Click options to the `run` command (after the existing `--min-qual` option):

```python
@click.option(
    "--platform",
    type=click.Choice(["hifi", "ont"], case_sensitive=False),
    default="hifi",
    help="Sequencing platform (default: hifi).",
)
@click.option(
    "--minimap2-preset",
    type=str,
    default=None,
    help="minimap2 -x preset (auto-selected from --platform if not set).",
)
```

Update the `run` function signature to include `platform: str` and `minimap2_preset: str | None`.

Add preset resolution at the start of the `run` function body:

```python
    preset = minimap2_preset or PLATFORM_PRESETS[platform]
```

Update `map_reads` call:

```python
    bam = map_reads(Path(input_path), ref, out, threads, preset=preset)
```

Update `call_variants_per_allele` call:

```python
    vcf_paths = call_variants_per_allele(
        bam,
        ref,
        alleles_result,
        out,
        clair3_model,
        threads,
        min_qual=min_qual,
        platform=platform,
        preset=preset,
    )
```

- [ ] **Step 5: Add `--platform` and `--minimap2-preset` to `call` subcommand**

Add Click options to `call` (after `--min-qual`):

```python
@click.option(
    "--platform",
    type=click.Choice(["hifi", "ont"], case_sensitive=False),
    default="hifi",
    help="Sequencing platform (default: hifi).",
)
@click.option(
    "--minimap2-preset",
    type=str,
    default=None,
    help="minimap2 -x preset (auto-selected from --platform if not set).",
)
```

Update the `call` function signature and body:

```python
def call(
    input_path: str,
    reference: str,
    alleles_json: str,
    output_dir: str,
    clair3_model: str,
    threads: int,
    min_qual: float,
    platform: str,
    minimap2_preset: str | None,
) -> None:
```

```python
    preset = minimap2_preset or PLATFORM_PRESETS[platform]
    vcfs = call_variants_per_allele(
        Path(input_path),
        Path(reference),
        alleles_data,
        Path(output_dir),
        clair3_model,
        threads,
        min_qual=min_qual,
        platform=platform,
        preset=preset,
    )
```

- [ ] **Step 6: Add `--minimap2-preset` to `map` subcommand**

Add Click option to `map` (after `--threads`):

```python
@click.option(
    "--minimap2-preset",
    type=str,
    default=None,
    help="minimap2 -x preset (default: map-hifi).",
)
```

Update the `map_cmd` function signature and body:

```python
def map_cmd(
    input_path: str,
    reference: str | None,
    output_dir: str,
    threads: int,
    minimap2_preset: str | None,
) -> None:
```

```python
    preset = minimap2_preset or "map-hifi"
    bam = map_reads(Path(input_path), ref, Path(output_dir), threads, preset=preset)
```

- [ ] **Step 7: Run all CLI tests**

Run: `uv run pytest tests/unit/test_cli.py -v --no-cov`

Expected: All PASS

- [ ] **Step 8: Commit**

```bash
git add src/muc_one_span/cli.py tests/unit/test_cli.py
git commit -m "feat: add --platform and --minimap2-preset CLI options for ONT support"
```

---

### Task 4: Update `scripts/batch_analyze.py`

**Files:**
- Modify: `scripts/batch_analyze.py`

- [ ] **Step 1: Add `--platform` argument support**

In `scripts/batch_analyze.py`, update `run_pipeline` to accept a `platform` parameter:

```python
def run_pipeline(sample_dir: Path, output_dir: Path, platform: str = "hifi") -> dict | None:
```

Add `--platform` to the `cmd` list (after `"--clair3-model"`):

```python
        "--platform",
        platform,
```

- [ ] **Step 2: Update `main()` to accept platform argument**

Update `main()` to parse a `--platform` CLI arg:

```python
def main():
    import argparse

    parser = argparse.ArgumentParser(description="Batch analyze test samples")
    parser.add_argument("samples_dir", nargs="?", default=str(SAMPLES_DIR))
    parser.add_argument("output_dir", nargs="?", default=str(OUTPUT_DIR))
    parser.add_argument("--platform", default="hifi", choices=["hifi", "ont"],
                        help="Sequencing platform (default: hifi)")
    args = parser.parse_args()

    samples_dir = Path(args.samples_dir)
    output_base = Path(args.output_dir)
```

Update the `run_pipeline` call:

```python
        pipeline_result = run_pipeline(sample_dir, out_dir, platform=args.platform)
```

- [ ] **Step 3: Commit**

```bash
git add scripts/batch_analyze.py
git commit -m "feat: add --platform flag to batch_analyze.py"
```

---

### Task 5: Run full test suite and quality checks

**Files:** None (verification only)

- [ ] **Step 1: Run unit tests with coverage**

Run: `uv run pytest tests/unit/ --cov=muc_one_span --cov-fail-under=80 -v`

Expected: All PASS, coverage >= 80%

- [ ] **Step 2: Run linter and formatter**

Run: `make lint && make format`

Expected: No issues

- [ ] **Step 3: Run type checker**

Run: `make type-check`

Expected: No errors

- [ ] **Step 4: Run full CI check**

Run: `make ci-check`

Expected: All checks pass

- [ ] **Step 5: Fix any issues found and commit**

If any CI checks fail, fix the issues and commit the fixes.

---

### Task 6: Documentation updates

**Files:**
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Add changelog entry**

Add under `## [Unreleased]` in `CHANGELOG.md`:

```markdown
### Added
- `--platform` CLI option on `run` and `call` subcommands to select sequencing platform (`hifi` or `ont`)
- `--minimap2-preset` CLI option on `run`, `map`, and `call` subcommands to override minimap2 alignment preset
- Auto-selection of minimap2 preset from platform (`hifi` -> `map-hifi`, `ont` -> `lr:hq`)
- ONT support: Clair3 `--platform=ont` and minimap2 `lr:hq` preset threaded through full pipeline
```

- [ ] **Step 2: Commit**

```bash
git add CHANGELOG.md
git commit -m "docs: add ONT platform support to changelog"
```
