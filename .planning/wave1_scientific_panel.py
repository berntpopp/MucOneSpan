"""Reproduce the small Wave 1 real-tool scientific preservation panel.

Run with the project Python, a source checkout/archive, existing read fixtures,
and a Conda tool environment. All generated artifacts remain outside the repo.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

SAMPLES = {
    "sample_normal_60_80": "hifi",
    "sample_dupc_60_80": "hifi",
    "sample_ont_dupc_60_80": "ont",
}


def normalize(text: str, source: Path, output: Path) -> str:
    """Normalize only checkout/output paths and VCF command timestamp headers."""
    return text.replace(str(output.resolve()), "<OUTPUT>").replace(
        str(source.resolve()), "<SOURCE>"
    )


def snapshot(source: Path, output: Path) -> dict[str, str]:
    """Hash complete scientific objects, sequences, and VCF genotype/phase data."""
    artifacts: dict[str, str] = {}
    for sample in SAMPLES:
        root = output / sample
        paths = [root / name for name in ("alleles.json", "repeats.json", "repeats.txt")]
        paths += sorted(root.glob("consensus_*.fa"))
        paths += sorted(root.glob("consensus_*_context.json"))
        paths += sorted(root.rglob("*.vcf")) + sorted(root.rglob("*.vcf.gz"))
        for path in paths:
            text = (
                gzip.decompress(path.read_bytes()).decode()
                if path.suffix == ".gz"
                else path.read_text()
            )
            text = normalize(text, source, output)
            if path.suffix == ".json":
                text = json.dumps(json.loads(text), sort_keys=True, separators=(",", ":"))
            elif ".vcf" in path.name:
                text = "\n".join(
                    line
                    for line in text.splitlines()
                    if not line.startswith(("##bcftools_", "##fileDate=", "##cmdline="))
                )
            artifacts[str(path.relative_to(output))] = hashlib.sha256(text.encode()).hexdigest()
        summary = json.loads(normalize((root / "summary.json").read_text(), source, output))
        # Status/report/configuration provenance intentionally changes in Wave 1.
        scientific = {
            key: summary[key]
            for key in ("alleles", "classifications", "tool_versions", "pipeline_version")
        }
        artifacts[f"{sample}/summary.scientific.json"] = hashlib.sha256(
            json.dumps(scientific, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
    (output / "scientific-sha256.json").write_text(json.dumps(artifacts, indent=2) + "\n")
    return artifacts


def run_panel(source: Path, data: Path, tools: Path, output: Path) -> None:
    """Run three fixed fixtures with unchanged science defaults and two threads."""
    output.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    env["PATH"] = str(tools / "bin") + os.pathsep + env["PATH"]
    env["PYTHONPATH"] = str(source / "src")
    records = []
    for sample, platform in SAMPLES.items():
        inputs = sorted((data / sample).glob("*.bam" if platform == "hifi" else "*.fastq"))
        if len(inputs) != 1:
            raise ValueError(f"Expected exactly one {platform} fixture for {sample}: {inputs}")
        command = [
            sys.executable,
            "-c",
            "from muc_one_span.cli import main; main()",
            "run",
            "--input",
            str(inputs[0]),
            "--output-dir",
            str(output / sample),
            "--clair3-model",
            str(tools / "bin/models" / platform),
            "--platform",
            platform,
            "--threads",
            "2",
            "--report",
        ]
        print(f"Running {sample} ({platform})", flush=True)
        with (output / f"{sample}.log").open("w") as log:
            result = subprocess.run(
                command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=900, check=False
            )
        records.append(
            {
                "sample": sample,
                "platform": platform,
                "command": command,
                "returncode": result.returncode,
                "input_sha256": hashlib.sha256(inputs[0].read_bytes()).hexdigest(),
            }
        )
        (output / "run-records.json").write_text(json.dumps(records, indent=2) + "\n")
        if result.returncode:
            raise RuntimeError(f"{sample} failed with exit {result.returncode}; inspect its log")
    snapshot(source, output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    run = sub.add_parser("run")
    for flag in ("source", "data", "tools", "output"):
        run.add_argument(f"--{flag}", type=Path, required=True)
    snap = sub.add_parser("snapshot")
    snap.add_argument("--source", type=Path, required=True)
    snap.add_argument("--output", type=Path, required=True)
    compare = sub.add_parser("compare")
    compare.add_argument("before", type=Path)
    compare.add_argument("after", type=Path)
    args = parser.parse_args()
    if args.action == "run":
        run_panel(
            args.source.resolve(), args.data.resolve(), args.tools.resolve(), args.output.resolve()
        )
    elif args.action == "snapshot":
        snapshot(args.source.resolve(), args.output.resolve())
    else:
        before = json.loads((args.before / "scientific-sha256.json").read_text())
        after = json.loads((args.after / "scientific-sha256.json").read_text())
        differences = sorted(
            key for key in before.keys() | after.keys() if before.get(key) != after.get(key)
        )
        print(
            json.dumps(
                {
                    "before_artifacts": len(before),
                    "after_artifacts": len(after),
                    "differences": differences,
                },
                indent=2,
            )
        )
        raise SystemExit(bool(differences))


if __name__ == "__main__":
    main()
