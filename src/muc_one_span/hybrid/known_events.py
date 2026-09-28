"""Site-table signatures of the repeat dictionary's mutation templates (Task 15k).

The low-accuracy-subset rules (``phase_quality``) must never drop a site whose minor
allele is how a known event looks in the phase site table. That look is derived from
the data, not assumed: every template of the bundled dictionary (insertions,
deletions and delete-inserts) is applied to each of its allowed parent units, and the
parent and the mutated unit are both run through ``phase_sites.features`` inside a
context of canonical units. Every non-column site whose allele differs is a
signature:

* a run site: (base, parent run length, mutated run length). The run can lengthen
  (dupC: a C7 run read as C8) or shorten (insG and delinsAT split the X unit's C7
  run, which the site table sees as C7 -> C6, C5 or C4; the deletions shorten C3/C4
  runs);
* an insertion slot: the inserted string (dupA: "A").

Column sites are not signatures: a column substitution or gap is not specific to an
event. A template/unit pair whose mutation changes no run or insertion slot is
listed in ``KnownEventSites.blind``; the rules cannot protect it. The bundled
dictionary has none since Task 15l: insG_pos54 in unit J inserts a G into the slot
right before a C run, which the site table did not record before. A custom
dictionary can have blind pairs; a run that uses the signatures logs them
(``warn_blind_templates``). The run-minority tier (``run_minor``) reads the run
signatures too: its ``known_events`` scope and its within-run test (insG and delinsAT
split an X unit's C7 run, so a carrier is seen only through its signature).

The signatures are position- and flank-agnostic: a run signature (base, parent
length, mutated length) protects every run of that base and length change anywhere
in the peak, and an insertion signature every slot with that inserted string, not
only the template's own units or neighbouring sequence. That is the conservative
direction: more sites are kept, none is dropped because of where it sits.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

from muc_one_span.config import RepeatDictionary, apply_mutation
from muc_one_span.hybrid.phase_sites import Site, features
from muc_one_span.settings import HybridSettings

RUN, INS = "run", "ins"  # site kinds of phase_sites.features (structural)
logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class KnownEventSites:
    """Protected site signatures of the dictionary's templates (see module docstring).

    ``runs`` holds (base, parent length, mutated length); ``inserted`` the inserted
    strings of insertion slots; ``blind`` the (template, unit) pairs with neither.
    """

    runs: frozenset[tuple[str, int, int]] = frozenset()
    inserted: frozenset[str] = frozenset()
    blind: tuple[tuple[str, str], ...] = ()


def template_signatures(
    parent: str, mutated: str, context: str, settings: HybridSettings
) -> tuple[set[tuple[str, int, int]], set[str]]:
    """Run and insertion-slot signatures of ``parent`` -> ``mutated`` in ``context``."""
    cons = context + parent + context
    (ref, mut), meta = features(cons, [cons, context + mutated + context], settings)
    runs: set[tuple[str, int, int]] = set()
    inserted: set[str] = set()
    for site in set(ref) | set(mut):
        a, b = ref.get(site), mut.get(site)
        if a == b:
            continue
        if site[0] == RUN and a is not None and b is not None:
            runs.add((meta[site][0], int(a), int(b)))
        elif site[0] == INS and b:
            inserted.add(str(b))
    return runs, inserted


def known_event_sites(rd: RepeatDictionary, settings: HybridSettings) -> KnownEventSites:
    """Signatures of every dictionary template in each of its allowed parent units."""
    context = rd.repeats[rd.canonical_repeat]
    runs: set[tuple[str, int, int]] = set()
    inserted: set[str] = set()
    blind: list[tuple[str, str]] = []
    for name, template in rd.mutations.items():
        for unit in template.get("allowed_repeats") or []:
            parent = rd.repeats.get(unit)
            if parent is None:
                continue
            mutated = apply_mutation(parent, template.get("changes") or [])
            r, i = template_signatures(parent, mutated, context, settings)
            if not r and not i:
                blind.append((name, unit))
            runs |= r
            inserted |= i
    return KnownEventSites(frozenset(runs), frozenset(inserted), tuple(blind))


def is_known_event_site(
    site: dict[str, Any], meta: dict[Site, tuple[str, int]], known: KnownEventSites
) -> bool:
    """True when the site's major -> minor change is a known event's signature."""
    kind = site["site"][0]
    if kind == RUN and site["site"] in meta:
        base = meta[site["site"]][0]
        return (base, site["major"], site["minor"]) in known.runs
    if kind == INS and str(site["minor"]).startswith(str(site["major"])):
        return str(site["minor"])[len(str(site["major"])) :] in known.inserted
    return False


def signatures_in_use(settings: HybridSettings) -> bool:
    """True when a setting makes a decision depend on the template signatures.

    The quality rules' guard (``phase_quality_alpha`` > 0), the run-minority tier's
    ``known_events`` scope and its within-run test (``phase_run_minor_in_run`` with the
    tier on) read them; the tier's clean test with scope "all" does not.
    """
    scope = settings.phase_run_minor_scope
    return (
        settings.phase_quality_alpha > 0
        or scope == "known_events"
        or (scope != "off" and settings.phase_run_minor_in_run)
    )


def warn_blind_templates(known: KnownEventSites, settings: HybridSettings) -> None:
    """Log the template/unit pairs no signature covers, when signatures are in use.

    A blind pair can be protected by no quality-rule guard and found by none of the
    run-minority tier's signature tests (``signatures_in_use``).
    """
    if not known.blind or not signatures_in_use(settings):
        return
    pairs = ", ".join(f"{name} in {unit}" for name, unit in known.blind)
    logger.warning(
        "hybrid: these dictionary templates have no phase-site signature, so neither "
        "the low-accuracy rules' guard nor the run-minority tier's signature tests "
        "cover them: %s",
        pairs,
    )
