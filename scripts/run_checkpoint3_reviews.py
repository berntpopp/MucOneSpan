#!/usr/bin/env python3
"""Run independent external model reviews for Checkpoint 3 (Implementation Review)."""

from __future__ import annotations

import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REVIEWS_DIR = ROOT / ".planning" / "reviews" / "checkpoint3_implementation"
REVIEWS_DIR.mkdir(parents=True, exist_ok=True)

REQUEST_FILE = REVIEWS_DIR / "review-request-implementation.md"
PROMPT_FILE = REVIEWS_DIR / "prompt.txt"


def prepare_prompt() -> None:
    request_text = REQUEST_FILE.read_text()

    # Capture diff against baseline
    git_diff = subprocess.run(
        [
            "git",
            "diff",
            "9535f7ee02033cda3da9b22c5b73752e5473b016..HEAD",
            "--",
            "src/",
            "tests/unit/",
        ],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    ).stdout

    # Capture compare_runs output if available
    compare_script = ROOT / "scripts" / "compare_runs.py"
    compare_out = subprocess.run(
        ["python3", str(compare_script)],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    ).stdout

    prompt = f"""You are acting as an exceptionally rigorous principal bioinformatician, software architect, statistician, and algorithm engineer.

Please perform an independent, in-depth, highly critical review of the candidate implementation and development benchmark results for MucOneSpan VNTR reconstruction.

=== REVIEW REQUEST & QUESTIONS ===
{request_text}

=== DEVELOPMENT BENCHMARK COMPARISON (BASELINE vs CANDIDATE ON 140 DEV DATASETS) ===
{compare_out}

=== GIT DIFF AGAINST BASELINE (v0.11.0 commit 9535f7e) ===
{git_diff[:50000]}

Please provide your rigorous critique, answering each specific question, evaluating all numerical acceptance gates, checking for subtle failure modes or overfitting risks, and issuing a clear final verdict.
"""
    PROMPT_FILE.write_text(prompt)


def run_claude() -> None:
    print("[Reviewer 1/2] Invoking Claude Code (model: claude-fable-5-1)...")
    out_json = REVIEWS_DIR / "claude-response.json"
    err_log = REVIEWS_DIR / "claude-stderr.log"

    cmd = [
        "claude",
        "-p",
        "--model",
        "claude-fable-5-1",
        "--effort",
        "high",
        "--permission-mode",
        "plan",
        "--tools",
        "Read,Glob,Grep",
        "--output-format",
        "json",
    ]
    t0 = time.monotonic()
    with PROMPT_FILE.open("r") as fin, out_json.open("w") as fout, err_log.open("w") as ferr:
        res = subprocess.run(cmd, stdin=fin, stdout=fout, stderr=ferr, cwd=str(ROOT))
    elapsed = time.monotonic() - t0
    print(f"  Claude finished in {elapsed:.1f}s (Exit code: {res.returncode})")


def run_codex() -> None:
    print("[Reviewer 2/2] Invoking Codex CLI (model: gpt-6-astra)...")
    out_jsonl = REVIEWS_DIR / "codex-events.jsonl"
    err_log = REVIEWS_DIR / "codex-stderr.log"

    cmd = [
        "codex",
        "exec",
        "--model",
        "gpt-6-astra",
        "--sandbox",
        "read-only",
        "--json",
        "-",
    ]
    t0 = time.monotonic()
    with PROMPT_FILE.open("r") as fin, out_jsonl.open("w") as fout, err_log.open("w") as ferr:
        res = subprocess.run(cmd, stdin=fin, stdout=fout, stderr=ferr, cwd=str(ROOT))
    elapsed = time.monotonic() - t0
    print(f"  Codex finished in {elapsed:.1f}s (Exit code: {res.returncode})")


def main() -> int:
    prepare_prompt()
    run_claude()
    run_codex()
    print("=== All Checkpoint 3 reviews completed! ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
