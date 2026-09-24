"""Profile variants: artefact, error and PCR levels layered on a built-in MucOneUp profile.

Rules for transforming base profiles:
- Smear/Chimera (artefact levels): Set molecules.smear_rate and molecules.chimera_rate;
  not applied to ont_genomic_targeted.
- Error levels: For profiles with "errors" dict, error="poor" scales mismatch_rate,
  insertion_rate, and deletion_rate by ``profiles.poor_error_scale`` (stutter
  unchanged); for profiles without "errors" (HiFi), sets
  config_overrides.pacbio_params.accuracy_mean = ``profiles.hifi_poor_accuracy_mean``.
  Calibrated levels leave the base untouched.
- PCR levels: For amplicon profiles (not ont_genomic_targeted), pcr="strong" sets
  amplicon_params.pcr_bias = {"preset": "madritsch2025_r10", "alpha":
  ``strong_pcr_alpha_factor`` x ``r10_pcr_alpha``}; pcr="none" sets
  {"preset": "no_bias"}. Calibrated levels leave the base untouched.
- All numbers come from `bench_config.ProfileConfig`.
- Variant name: Encodes profile base name and all levels (smear/chimera omitted for
  genomic). Provenance records base profile name and SHA256. The file is
  ``<variant name>__<SHA-256 of its content>.json``, so variants with the same
  levels but different settings (or base profile) never overwrite or reuse
  each other.
"""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any

from .bench_config import DEFAULT_BENCH_CONFIG, ProfileConfig
from .design import Design
from .muconeup import BUILTIN_PROFILE

PCR_PRESET_R10 = "madritsch2025_r10"  # MucOneUp pcr_bias preset names
PCR_PRESET_NONE = "no_bias"


def builtin_profile_dir(explicit: Path | None) -> Path:
    if explicit is not None:
        return explicit
    spec = importlib.util.find_spec("muc_one_up")
    if spec is not None and spec.origin:
        candidate = Path(spec.origin).parent / "data" / "read_profiles"
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError(
        "MucOneUp read profiles not found; pass --muconeup-profiles "
        "<MucOneUp checkout>/muc_one_up/data/read_profiles"
    )


def variant_name(design: Design) -> str:
    base = BUILTIN_PROFILE[design.profile]
    parts = [f"pcr-{design.pcr}", f"err-{design.error}"]
    if design.profile != "ont_genomic_targeted":
        parts = [f"s{design.smear}", f"c{design.chimera}", *parts]
    return f"{base}__{'_'.join(parts)}"


def _apply(data: dict[str, Any], design: Design, cfg: ProfileConfig) -> None:
    mol = data.setdefault("molecules", {})
    if design.profile != "ont_genomic_targeted":
        mol["smear_rate"] = design.smear
        mol["chimera_rate"] = design.chimera
    if design.error == "poor":
        if "errors" in data:
            for key in ("mismatch_rate", "insertion_rate", "deletion_rate"):
                data["errors"][key] = data["errors"][key] * cfg.poor_error_scale
        else:
            data.setdefault("config_overrides", {}).setdefault("pacbio_params", {})[
                "accuracy_mean"
            ] = cfg.hifi_poor_accuracy_mean
    if design.profile != "ont_genomic_targeted" and design.pcr != "calibrated":
        pcr = (
            {"preset": PCR_PRESET_NONE}
            if design.pcr == "none"
            else {
                "preset": PCR_PRESET_R10,
                "alpha": cfg.strong_pcr_alpha_factor * cfg.r10_pcr_alpha,
            }
        )
        overrides = data.setdefault("config_overrides", {})
        overrides.setdefault("amplicon_params", {})["pcr_bias"] = pcr


def write_variant(
    base_profile: Path,
    design: Design,
    out_dir: Path,
    config: ProfileConfig = DEFAULT_BENCH_CONFIG.profiles,
) -> tuple[Path, str]:
    raw = base_profile.read_bytes()
    base = json.loads(raw)
    data = copy.deepcopy(base)
    _apply(data, design, config)  # calibrated levels still get a variant file, for provenance
    data["name"] = variant_name(design)
    data.setdefault("provenance", {})["derived_from"] = {
        "name": base["name"],
        "sha256": hashlib.sha256(raw).hexdigest(),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    text = json.dumps(data, indent=2, sort_keys=True) + "\n"
    sha = hashlib.sha256(text.encode()).hexdigest()
    # Content-addressed: equal levels under different settings never share a file.
    path = out_dir / f"{data['name']}__{sha}.json"
    if not path.exists() or path.read_text() != text:
        path.write_text(text)
    return path, sha
