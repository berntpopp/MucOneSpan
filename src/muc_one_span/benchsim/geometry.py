"""Case geometry: amplicon primers and intervals, VNTR bounds in both frames.

Frames:

- **haplotype**: a truth FASTA record, ``flanking_left + VNTR + flanking_right``
  (validated by `evaluation.truth.load_truth`).
- **source**: the frame of the read truth ``src_start``/``src_end``. For
  amplicon profiles it is the primer-to-primer amplicon that MucOneUp extracts
  (exact, unique primer matches, both primers included). For genomic
  fragments it is the flanked source ``--flank-fasta left + haplotype + right``.

The VNTR runs from the motif-1 start to the motif-9 end, i.e. the truth
haplotype's repeat sequence. `generate` records the result as
``case.json["geometry"]``; `source_inputs` turns it back into source
sequences and VNTR bounds for `realism.read_metrics`.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from muc_one_span.evaluation.models import TruthSample
from muc_one_span.evaluation.truth import fasta_records
from muc_one_span.nomenclature import revcomp

GENOMIC = "ont_genomic_targeted"


def primer_pair(muconeup_config: Path) -> tuple[str, str]:
    """Forward and reverse amplicon primers from a MucOneUp config file.

    Raises:
        KeyError: If ``amplicon_params`` lacks either primer.
    """
    params = json.loads(muconeup_config.read_text()).get("amplicon_params") or {}
    if not params.get("forward_primer") or not params.get("reverse_primer"):
        raise KeyError(f"{muconeup_config}: amplicon_params primers are missing")
    return str(params["forward_primer"]).upper(), str(params["reverse_primer"]).upper()


def amplicon_interval(sequence: str, forward: str, reverse: str) -> tuple[int, int]:
    """``[start, end)`` of the primer-to-primer amplicon (MucOneUp extraction).

    Raises:
        ValueError: If a primer site is missing or not unique.
    """
    seq, rev = sequence.upper(), revcomp(reverse.upper())
    f_sites = [m.start() for m in re.finditer(f"(?={forward.upper()})", seq)]
    r_sites = [m.start() for m in re.finditer(f"(?={rev})", seq)]
    if len(f_sites) != 1 or len(r_sites) != 1 or r_sites[0] < f_sites[0]:
        raise ValueError(f"primer sites not unique ({f_sites}, {r_sites})")
    return f_sites[0], r_sites[0] + len(rev)


def haplotype_sequences(truth_fasta: Path) -> dict[int, str]:
    """Haplotype id -> uppercase sequence from a MucOneUp truth FASTA."""
    out = {}
    for index, (name, seq) in enumerate(fasta_records(truth_fasta).items(), 1):
        match = re.search(r"haplotype_(\d+)", name)
        out[int(match.group(1)) if match else index] = seq.upper()
    return out


def flank_sequences(flank_fasta: Path | None) -> tuple[str, str]:
    """The ``left``/``right`` records of a ``--flank-fasta`` (empty without one)."""
    if flank_fasta is None:
        return "", ""
    records = fasta_records(flank_fasta)
    return records.get("left", "").upper(), records.get("right", "").upper()


def vntr_bounds(
    truth: TruthSample, flanking_left: str, haplotypes: dict[int, str]
) -> dict[int, tuple[int, int]]:
    """VNTR ``[start, end)`` per haplotype in the haplotype frame (validated)."""
    out = {}
    for hap, h in enumerate(truth.haplotypes, 1):
        start = len(flanking_left)
        end = start + len(h.sequence)
        if haplotypes[hap][start:end] != h.sequence.upper():
            raise ValueError(f"haplotype {hap}: VNTR does not match the truth FASTA")
        out[hap] = (start, end)
    return out


def case_geometry(
    profile: str,
    haplotypes: dict[int, str],
    vntr: dict[int, tuple[int, int]],
    *,
    primers: tuple[str, str] | None,
    flank_ext: tuple[int, int],
) -> dict[str, Any]:
    """JSON-ready geometry of one case (see the module docstring for frames).

    Args:
        profile: Benchmark profile.
        haplotypes: Haplotype id -> truth FASTA sequence.
        vntr: VNTR bounds in the haplotype frame (`vntr_bounds`).
        primers: Forward/reverse primers (required for amplicon profiles).
        flank_ext: ``--flank-fasta`` left/right lengths (genomic).

    Raises:
        ValueError: If an amplicon profile has no primers or no unique sites.
    """
    genomic = profile == GENOMIC
    out: dict[str, Any] = {
        "frame": "flanked_source" if genomic else "amplicon",
        "primers": None,
        "amplicon": None,
        "flank_ext": list(flank_ext) if genomic else None,
        "vntr_haplotype": {str(h): list(v) for h, v in vntr.items()},
    }
    if genomic:
        shift = dict.fromkeys(vntr, flank_ext[0])
    else:
        if primers is None:
            raise ValueError(f"{profile}: amplicon primers are required")
        amp = {h: amplicon_interval(haplotypes[h], *primers) for h in vntr}
        out["primers"] = {"forward": primers[0], "reverse": primers[1]}
        out["amplicon"] = {str(h): list(a) for h, a in amp.items()}
        shift = {h: -a[0] for h, a in amp.items()}
    out["vntr_source"] = {str(h): [s + shift[h], e + shift[h]] for h, (s, e) in vntr.items()}
    return out


def source_inputs(
    geometry: dict[str, Any],
    haplotypes: dict[int, str],
    flanks: tuple[str, str] = ("", ""),
) -> tuple[dict[int, str], dict[int, tuple[int, int]]]:
    """Source sequences and VNTR bounds in the read-truth (source) frame.

    Args:
        geometry: A `case_geometry` dict (``case.json["geometry"]``).
        haplotypes: Haplotype id -> truth FASTA sequence.
        flanks: The generation ``--flank-fasta`` left/right sequences (genomic).

    Raises:
        ValueError: If the flank lengths differ from the recorded ones.
    """
    vntr = {int(h): (int(v[0]), int(v[1])) for h, v in geometry["vntr_source"].items()}
    if geometry["frame"] == "amplicon":
        amp = geometry["amplicon"]
        return {h: haplotypes[h][amp[str(h)][0] : amp[str(h)][1]] for h in vntr}, vntr
    if [len(flanks[0]), len(flanks[1])] != list(geometry["flank_ext"]):
        raise ValueError(
            f"flank lengths {[len(f) for f in flanks]} != recorded {geometry['flank_ext']}; "
            "pass the --flank-fasta used at generation"
        )
    return {h: flanks[0] + haplotypes[h] + flanks[1] for h in vntr}, vntr
