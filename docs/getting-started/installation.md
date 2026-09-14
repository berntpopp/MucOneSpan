# Installation

MucOneSpan requires Python 3.10 or newer; CI tests Python 3.10–3.14. Full pipeline
analysis also needs external bioinformatics tools and a Clair3 model.

## Install the Python package

```bash
pip install 'muc_one_span[report] @ git+https://github.com/berntpopp/MucOneSpan.git@v0.10.0'
muconespan --version
```

The `[report]` extra enables HTML reports through Jinja2. Python installation
includes Click, PyYAML, and packaged reference/repeat data. External alignment
and variant-calling tools are installed separately.

## Install from source

```bash
git clone https://github.com/berntpopp/MucOneSpan.git
cd MucOneSpan
pip install -e ".[report]"
muconespan --version
```

For development with the locked uv environment:

```bash
git clone https://github.com/berntpopp/MucOneSpan.git
cd MucOneSpan
make dev
uv run muconespan --version
make ci-check
```

`make dev` installs all development groups and optional extras. Use `uv run` to
invoke commands from this project environment, or activate `.venv/bin/activate`.
See the [developer guide](../development.md) for checks, hooks, and architecture.

## External tool dependencies

| Tool | Purpose |
| --- | --- |
| minimap2 | Align long reads to the ladder reference |
| samtools | Process and index BAM files |
| bcftools | Filter VCF files and construct consensus sequences |
| Clair3 and a platform-appropriate model | Call variants for HiFi or ONT reads |

The `muconespan run` command uses all four tools. The `ladder` and `classify`
commands can run without external tools. Linux is the primary pipeline platform;
check Clair3's platform requirements before installing on another operating system.

The repository's conda environment supplies pinned minimap2, samtools, bcftools,
and htslib versions:

```bash
conda env create -f conda/environment.yml
conda activate muconespan-tools
```

Install Clair3 and its models following the
[Clair3 installation instructions](https://github.com/HKU-BAL/Clair3).
Make its `run_clair3.sh` available on `PATH` and pass the correct model location
with `--clair3-model`.

Clair3 uses its own Python dependencies. MucOneSpan removes the project virtual
environment's executable path when launching external tools so it does not shadow
that environment. Activate the tool environment or expose its executables:

```bash
export PATH="/path/to/clair3/environment/bin:$PATH"
```

## Docker

```bash
docker pull ghcr.io/berntpopp/muconespan:latest

docker run --rm \
  -v "$(pwd)/data:/data" \
  ghcr.io/berntpopp/muconespan:latest \
  run --input /data/reads.bam --output-dir /data/results/
```

Use a version tag such as `0.10.0` in place of `latest` to select a release.

## Verify installation

```bash
muconespan --version
muconespan --help
muconespan ladder --output test_ladder.fa
```

Then verify external commands and run a real input sample:

```bash
command -v minimap2 samtools bcftools run_clair3.sh
muconespan run \
  --input reads.fastq \
  --output-dir results/ \
  --clair3-model /path/to/clair3/models/hifi
```

For ONT data, add `--platform ont` and select an ONT model.

If `muconespan` is not found after `make dev`, use `uv run muconespan` from the
checkout or activate its `.venv`. If an external command is missing, activate the
appropriate tool environment and check its `PATH` before rerunning the pipeline.

## Next steps

- [Quick Start](quickstart.md): Run your first analysis.
- [Core Concepts](concepts.md): Understand the pipeline architecture.
- [Differences from the Published Method](deviations.md): Compare the implementation.
