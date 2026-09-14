#!/usr/bin/env bash
# Smoke-test the installed image, including Clair3 inference, using Docker only.
set -euo pipefail

container_image=${1:?Usage: bash scripts/check_container.sh IMAGE}

docker run --rm "$container_image" --help
docker run --rm "$container_image" --version

docker run --rm -i --entrypoint bash "$container_image" -s <<'CONTAINER'
set -euo pipefail

# Check the runtime stage, not the builder's source tree or Python environment.
test "$(id -u)" -ne 0
test ! -e /app
minimap2 --version
samtools --version
bcftools --version
run_clair3.sh --version
pypy3 --version

smoke_dir=$(mktemp -d)
cd "$smoke_dir"

python - <<'PY'
from importlib.metadata import version
from pathlib import Path
import random
import shutil
import subprocess

assert subprocess.check_output(["open-pacmuci", "--version"], text=True).strip() == (
    f"open-pacmuci, version {version('open-pacmuci')}"
)
model_root = Path(shutil.which("run_clair3.sh")).parent / "models"
for platform in ("hifi", "ont"):
    for model in ("pileup", "full_alignment"):
        assert (model_root / platform / f"{model}.index").is_file()
        assert list((model_root / platform).glob(f"{model}.data-*"))

# Twenty identical HiFi reads carry one known SNP in a unique reference.
rng = random.Random(42)
reference = "".join(rng.choices("ACGT", k=3000))
alternate_base = next(base for base in "ACGT" if base != reference[1500])
alternate = reference[:1500] + alternate_base + reference[1501:]
Path("ref.fa").write_text(">contig_test\n" + reference + "\n")
Path("reads.fq").write_text(
    "".join(f"@read{i}\n{alternate}\n+\n{'I' * len(alternate)}\n" for i in range(20))
)
Path("expected.tsv").write_text(f"contig_test\t1501\t{reference[1500]}\t{alternate_base}\n")
PY

# Ladder generation also verifies that the installed wheel contains repeat data.
open-pacmuci ladder --min-units 20 --max-units 22 --output ladder.fa
samtools faidx ladder.fa
python - <<'PY'
from pathlib import Path

entries = [line.split("\t") for line in Path("ladder.fa.fai").read_text().splitlines()]
assert [entry[0] for entry in entries] == ["contig_20", "contig_21", "contig_22"]
assert [int(entry[1]) for entry in entries] == [2740, 2800, 2860]
PY

minimap2 -ax map-hifi ref.fa reads.fq | samtools sort -o reads.bam
samtools index reads.bam
samtools faidx ref.fa
clair3_models="$(dirname "$(command -v run_clair3.sh)")/models"
run_clair3.sh \
    --bam_fn="$smoke_dir/reads.bam" \
    --ref_fn="$smoke_dir/ref.fa" \
    --output="$smoke_dir/calls" \
    --threads=2 \
    --platform=hifi \
    --model_path="$clair3_models/hifi" \
    --include_all_ctgs \
    --chunk_size=5000
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT\t%DP]\n' \
    calls/merge_output.vcf.gz > observed.tsv
python - <<'PY'
from pathlib import Path

expected = Path("expected.tsv").read_text().strip().split("\t")
observed = [line.split("\t") for line in Path("observed.tsv").read_text().splitlines()]
assert len(observed) == 1, observed
assert observed[0][:4] == expected, observed
assert observed[0][4:] == ["PASS", "1/1", "20"], observed
print("Container smoke test passed: packaged tools, data, models, and Clair3 SNP call.")
PY
CONTAINER
