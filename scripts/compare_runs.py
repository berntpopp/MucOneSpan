#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

base_path = (
    Path(sys.argv[1])
    if len(sys.argv) > 1
    else Path("tests/results/baseline_dev_140/evaluation_report.json")
)
cand_path = (
    Path(sys.argv[2])
    if len(sys.argv) > 2
    else Path("tests/results/candidate_dev_140/evaluation_report.json")
)

base: dict[str, Any] = (
    json.loads(base_path.read_text())
    if base_path.exists()
    else {"totals": {"metrics": {}}, "samples": []}
)
cand: dict[str, Any] = (
    json.loads(cand_path.read_text())
    if cand_path.exists()
    else {"totals": {"metrics": {}}, "samples": []}
)

print("=== BASELINE METRICS ===")
for k, v in base["totals"]["metrics"].items():
    print(f"{k}: {v}")
print("\n=== CANDIDATE METRICS ===")
for k, v in cand["totals"]["metrics"].items():
    print(f"{k}: {v}")

print("\n=== ACCURACY COMPARISON ===")
for k in [
    "exact_sample_reconstruction",
    "sequence_accuracy",
    "structure_accuracy",
    "normal_specificity",
    "sample_alarm_sensitivity",
]:
    b_val = base["totals"].get(k)
    c_val = cand["totals"].get(k)
    print(f"{k}: base={b_val} cand={c_val}")

b_samples: dict[str, dict[str, Any]] = {s["sample"]: s for s in base["samples"]}
c_samples: dict[str, dict[str, Any]] = {s["sample"]: s for s in cand["samples"]}

for key in ["all_sequences_exact", "all_structures_exact", "sequence_exact", "structure_exact"]:
    hw: list[str] = []
    hl: list[str] = []
    ow: list[str] = []
    ol: list[str] = []
    for s_name, c_row in c_samples.items():
        if s_name not in b_samples:
            continue
        b_row = b_samples[s_name]
        c_val = c_row["metrics"][key]["min"]
        b_val = b_row["metrics"][key]["min"]
        platform = "ont" if s_name.endswith("_ont") else "hifi"
        if c_val > b_val:
            (hw if platform == "hifi" else ow).append(s_name)
        elif c_val < b_val:
            (hl if platform == "hifi" else ol).append(s_name)
    print(f"\nMetric: {key}")
    print(f"  HiFi Wins: {len(hw)} {hw}")
    print(f"  HiFi Losses: {len(hl)} {hl}")
    print(f"  HiFi Net: {len(hw) - len(hl)}")
    print(f"  ONT  Wins: {len(ow)} {ow}")
    print(f"  ONT  Losses: {len(ol)} {ol}")
    print(f"  ONT  Net: {len(ow) - len(ol)}")
