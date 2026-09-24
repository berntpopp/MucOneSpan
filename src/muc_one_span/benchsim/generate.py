"""Generate one validated benchmark case: MucOneUp truth, reads, read truth, ledger.

Case layout: ``<out_root>/<split>/<design_id>/{truth/, reads/, case.json}``.
``truth/`` holds the ``simulate`` outputs plus copies of the reads'
``*_read_truth.tsv.gz`` and ``*_metadata.tsv``; ``reads/`` holds the FASTQ.

Statuses: ``ok``; ``design_invalid`` when MucOneUp truth fails validation or
does not carry exactly the designed event on the designed haplotype (the case
is recorded, never scored); ``generation_failed`` for tool or output failures.
Only ``ok`` cases enter the sealed ledger; a verified-complete case is reused.

Genomic span interval: MucOneSpan has no motif-1/motif-9 anchor coordinates,
so the span is the VNTR itself inside the flanked source,
``[ext_left + len(flank_left), + len(VNTR))``, where ``flank_left`` is the
MucOneUp flank in the truth FASTA and ``ext_left`` the ``--flank-fasta`` left
record.

``case.json["geometry"]`` (`geometry.case_geometry`) records the amplicon
primer pair (MucOneUp config ``amplicon_params``), each haplotype's amplicon
interval and the VNTR bounds in the haplotype and read-truth source frames.
"""

from __future__ import annotations

import fcntl
import json
import shutil
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from muc_one_span.config import RepeatDictionary, load_repeat_dictionary
from muc_one_span.durable_ledger import DurableLedger, LedgerEntry, compute_sha256
from muc_one_span.evaluation.models import TruthSample
from muc_one_span.evaluation.truth import fasta_records, load_truth
from muc_one_span.tools import run_tool

from .depth import amplicon_templates, capped_minor_share, genomic_reads, pcr_minor_share
from .design import Design
from .geometry import case_geometry, haplotype_sequences, primer_pair, vntr_bounds
from .muconeup import BUILTIN_PROFILE, reads_args, simulate_args
from .profiles import variant_name, write_variant
from .read_truth import composition, load_read_truth, realized_depth
from .structures import prepare_structure, read_chains, scale_targets

FASTQ = {
    "ont_amplicon_r10": "{}_amplicon_ont.fastq",
    "hifi_amplicon": "{}_amplicon_hifi.fastq",
    "ont_genomic_targeted": "{}_ont_fragments.fastq",
}
FRAGMENT_DEFAULT = (5000.0, 0.5)  # MucOneUp FragmentModel defaults


@dataclass(frozen=True)
class GenerateContext:
    """Run-wide inputs shared by every case."""

    executable: str
    config: Path
    out_root: Path
    profile_dir: Path
    flank_fasta: Path | None
    structure_pool: Path | None
    muconeup_version: str


@contextmanager
def _ledger(split_dir: Path) -> Iterator[DurableLedger]:
    """Open the split ledger under a benchsim-wide lock.

    ``DurableLedger.__init__`` may rewrite the public ledger without its own
    lock, so parallel workers serialize ledger access on a separate lock file.
    """
    split_dir.mkdir(parents=True, exist_ok=True)
    with (split_dir / ".benchsim.lock").open("w") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield DurableLedger(split_dir)
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


class DesignInvalidError(ValueError):
    """MucOneUp truth does not realize the design."""


def _write_json(path: Path, data: dict[str, Any]) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def write_manifest(split_dir: Path, cases: list[dict[str, Any]]) -> Path:
    """Write ``manifest.jsonl`` with one case per line, sorted by ``design_id``."""
    split_dir.mkdir(parents=True, exist_ok=True)
    path = split_dir / "manifest.jsonl"
    tmp = path.with_suffix(".tmp")
    rows = sorted(cases, key=lambda c: str(c["design_id"]))
    tmp.write_text("".join(json.dumps(c, sort_keys=True) + "\n" for c in rows))
    tmp.replace(path)
    return path


def validate_events(truth: TruthSample, design: Design) -> list[list[int]]:
    """Return actual ``[hap, repeat]`` targets; raise if they differ from the design."""
    found = [
        (hap, ev.repeat_index, ev.name)
        for hap, h in enumerate(truth.haplotypes, 1)
        for ev in h.events
    ]
    if any(r is None for _, r, _ in found):
        raise DesignInvalidError(f"truth event without a repeat position: {found}")
    if design.event is None:
        if found:
            raise DesignInvalidError(f"normal design but truth has events {found}")
    else:
        hap = design.targets[0][0]
        if [(h, n) for h, _, n in found] != [(hap, design.event)]:
            raise DesignInvalidError(
                f"truth events {found} != design event {design.event} on hap {hap}"
            )
    if design.delta_class == "0_identical":
        # MucOneUp may convert the target unit to an allowed repeat; mask event positions
        masked = {r for _, r, _ in found}
        plain = [
            tuple("*" if i in masked else u.split(":")[0] for i, u in enumerate(h.structure, 1))
            for h in truth.haplotypes
        ]
        if len(set(plain)) != 1:
            raise DesignInvalidError("0_identical design but haplotype structures differ")
    return [[h, r] for h, r, _ in found if r is not None]


def _flank_lengths(flank_fasta: Path | None) -> tuple[int, int]:
    if flank_fasta is None:
        return 0, 0
    records = fasta_records(flank_fasta)
    return len(records.get("left", "")), len(records.get("right", ""))


def _span(
    truth: TruthSample, rd: RepeatDictionary, flank_fasta: Path | None
) -> tuple[dict[int, tuple[int, int]], dict[int, int]]:
    """Per-haplotype VNTR interval and source length in the flanked source."""
    left, right = _flank_lengths(flank_fasta)
    span, source = {}, {}
    for hap, h in enumerate(truth.haplotypes, 1):
        start = left + len(rd.flanking_left)
        span[hap] = (start, start + len(h.sequence))
        source[hap] = start + len(h.sequence) + len(rd.flanking_right) + right
    return span, source


def _amount(
    design: Design,
    base_profile: dict[str, Any],
    truth: TruthSample,
    span: dict[int, tuple[int, int]],
    source: dict[int, int],
) -> tuple[int, bool]:
    """Requested reads (genomic) or templates (amplicon), and whether the amount was capped."""
    if design.profile == "ont_genomic_targeted":
        frag = base_profile.get("fragments") or {}
        median = float(frag.get("length_median", FRAGMENT_DEFAULT[0]))
        sigma = float(frag.get("length_sigma", FRAGMENT_DEFAULT[1]))
        seed = design.read_seed % 2**32
        return max(
            genomic_reads(design.depth, source[h], lo, hi, median, sigma, seed=seed)
            for h, (lo, hi) in span.items()
        ), False
    concatemer = float((base_profile.get("molecules") or {}).get("concatemer_rate", 0.0))
    counts = [len(h.structure) for h in truth.haplotypes]
    share, capped = capped_minor_share(pcr_minor_share((counts[0], counts[-1]), design.pcr))
    rate = design.smear + design.chimera + concatemer
    return amplicon_templates(design.depth, rate, share), capped


def _fastq_records(path: Path) -> int:
    with path.open() as handle:
        lines = sum(1 for line in handle if line.strip())
    if lines % 4:
        raise ValueError(f"{path.name}: truncated FASTQ ({lines} lines)")
    return lines // 4


def _one(root: Path, pattern: str, required: bool = True) -> Path | None:
    found = sorted(root.glob(pattern))
    if len(found) > 1 or (required and not found):
        raise ValueError(f"expected one {pattern} in {root}, found {len(found)}")
    return found[0] if found else None


def _simulate(
    design: Design, ctx: GenerateContext, case_dir: Path, rd: RepeatDictionary, case: dict[str, Any]
) -> tuple[Path, TruthSample]:
    structure, info = prepare_structure(
        design,
        ctx.executable,
        ctx.config,
        ctx.structure_pool,
        case_dir / "structure",
        run_tool,
        known=set(rd.repeats),
    )
    case.update(info)
    targets = design.targets
    if structure is not None and info["structure_source"] == "pool":
        lengths = tuple(len(c) for c in read_chains(structure))
        targets = scale_targets(targets, design.lengths, lengths)
    case["requested_targets"] = [list(t) for t in targets]
    truth_dir = case_dir / "truth"
    lengths_arg = None if structure is not None else design.lengths
    run_tool(
        simulate_args(
            ctx.executable,
            ctx.config,
            truth_dir,
            design.design_id,
            design.bio_seed,
            lengths_arg,
            structure,
            design.event,
            targets,
        )
    )
    return truth_dir, _validated(truth_dir, rd, design, case)


def _validated(
    truth_dir: Path, rd: RepeatDictionary, design: Design, case: dict[str, Any]
) -> TruthSample:
    try:
        truth = load_truth(truth_dir, rd)
        case["actual_targets"] = validate_events(truth, design)
    except ValueError as exc:
        raise DesignInvalidError(str(exc)) from exc
    case["actual_lengths"] = [len(h.structure) for h in truth.haplotypes]
    return truth


def _reads(
    design: Design,
    ctx: GenerateContext,
    case_dir: Path,
    truth_dir: Path,
    truth: TruthSample,
    rd: RepeatDictionary,
    case: dict[str, Any],
) -> tuple[Path, Path, int]:
    base_path = ctx.profile_dir / f"{BUILTIN_PROFILE[design.profile]}.json"
    base_profile = json.loads(base_path.read_text())
    variant, sha = write_variant(base_path, design, ctx.out_root / "profiles")
    case.update(profile_variant=variant_name(design), profile_sha256=sha)
    span, source = _span(truth, rd, ctx.flank_fasta)
    genomic = design.profile == "ont_genomic_targeted"
    amount, capped = _amount(design, base_profile, truth, span, source)
    case["requested_amount"] = amount
    case["amount_capped"] = capped
    case["amount_unit"] = "reads" if genomic else "templates"
    truth_fa = _one(truth_dir, "*.simulated.fa")
    if truth_fa is None:  # _one(required=True) raises first; explicit for type narrowing
        raise ValueError(f"no *.simulated.fa in {truth_dir}")
    haplotypes = haplotype_sequences(truth_fa)
    case["geometry"] = case_geometry(
        design.profile,
        haplotypes,
        vntr_bounds(truth, rd.flanking_left, haplotypes),
        primers=None if genomic else primer_pair(ctx.config),
        flank_ext=_flank_lengths(ctx.flank_fasta),
    )
    reads_dir = case_dir / "reads"
    run_tool(
        reads_args(
            ctx.executable,
            ctx.config,
            design.profile,
            truth_fa,
            reads_dir,
            design.design_id,
            design.read_seed,
            amount=amount,
            profile_path=variant,
            flank_fasta=ctx.flank_fasta,
            pcr_preset=None,
        )
    )
    fastq = reads_dir / FASTQ[design.profile].format(design.design_id)
    manifest = reads_dir / f"{design.design_id}_read_truth.tsv.gz"
    for path in (fastq, manifest):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"MucOneUp did not write {path.name}")
    metadata = _one(reads_dir, "*_metadata.tsv", required=False)
    for copy in (manifest, metadata):
        if copy is not None:
            shutil.copy2(copy, truth_dir / copy.name)
    rows = load_read_truth(manifest)
    n_fastq = _fastq_records(fastq)
    if n_fastq != len(rows):
        raise RuntimeError(f"FASTQ has {n_fastq} records but read truth has {len(rows)} rows")
    if genomic:
        case["span"] = {str(h): list(v) for h, v in span.items()}
    depth = realized_depth(rows, span if genomic else None)
    case["realized_depth"] = {str(h): n for h, n in depth.items()}
    comp = composition(rows)
    comp["reads_per_hap"] = {str(h): n for h, n in comp["reads_per_hap"].items()}
    case["composition"] = comp
    case["hashes"] = {
        "truth_fa": compute_sha256(truth_fa),
        "fastq": compute_sha256(fastq),
        "read_truth": compute_sha256(manifest),
    }
    return truth_fa, fastq, len(rows)


def _commit(
    design: Design, ctx: GenerateContext, case: dict[str, Any], truth_fa: Path, fastq: Path, n: int
) -> None:
    suffix = design.design_id.rsplit("-", 1)[-1]
    with _ledger(ctx.out_root / design.split) as ledger:
        ledger.commit_entry(
            LedgerEntry(
                design_id=int(suffix) if suffix.isdigit() else 0,
                design_name=design.design_id,
                token=design.design_id,
                split=design.split,
                category=f"{design.profile}:{design.delta_class}:{design.event or 'normal'}",
                platform=design.profile,
                lengths=list(case["actual_lengths"]),
                mutation=design.event,
                targets=case["actual_targets"],
                truth_fa=str(truth_fa.relative_to(ledger.output_dir)),
                truth_sha256=case["hashes"]["truth_fa"],
                reads_file=str(fastq.relative_to(ledger.output_dir)),
                reads_sha256=case["hashes"]["fastq"],
                usable_records=n,
            )
        )


def generate_case(design: Design, ctx: GenerateContext) -> dict[str, Any]:
    """Generate, validate and record one case; returns the ``case.json`` dict."""
    split_dir = (ctx.out_root / design.split).resolve()
    case_dir = split_dir / design.design_id
    case_json = case_dir / "case.json"
    if case_json.is_file():
        with _ledger(split_dir) as ledger:
            done = ledger.is_verified_complete(design.design_id, design.profile)
        if done:
            saved: dict[str, Any] = json.loads(case_json.read_text())
            return saved
    for stale in ("truth", "reads", "structure"):
        shutil.rmtree(case_dir / stale, ignore_errors=True)
    case_dir.mkdir(parents=True, exist_ok=True)
    case: dict[str, Any] = {
        "design_id": design.design_id,
        "split": design.split,
        "profile": design.profile,
        "design": design.to_dict(),
        "muconeup_version": ctx.muconeup_version,
        "status": "generation_failed",
        "error": None,
    }
    rd = load_repeat_dictionary()
    try:
        truth_dir, truth = _simulate(design, ctx, case_dir, rd, case)
        truth_fa, fastq, n = _reads(design, ctx, case_dir, truth_dir, truth, rd, case)
        _commit(design, ctx, case, truth_fa, fastq, n)
        case["status"] = "ok"
    except DesignInvalidError as exc:
        case.update(status="design_invalid", error=str(exc))
    except (OSError, RuntimeError, ValueError, KeyError, TimeoutError) as exc:
        case["error"] = f"{type(exc).__name__}: {exc}"
    _write_json(case_json, case)
    result: dict[str, Any] = json.loads(case_json.read_text())  # JSON-normal, as on rerun
    return result
