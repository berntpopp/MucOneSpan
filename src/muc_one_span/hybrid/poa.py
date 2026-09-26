"""Partial-order-alignment consensus backends, selected explicitly by setting.

There is no silent fallback: results must be reproducible from the recorded
``hybrid.poa_backend`` setting. The prototype evidence was produced with pyabpoa.
Both backends are core dependencies since v0.17.0 and are imported at module load.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Protocol

import pyabpoa
import spoa


class PoaBackend(Protocol):
    """Minimal consensus interface used by the engine."""

    name: str

    def consensus(self, seqs: list[str]) -> str: ...


class _Abpoa:
    name = "pyabpoa"

    def __init__(self) -> None:
        self._aligner = pyabpoa.msa_aligner(aln_mode="g")

    def consensus(self, seqs: list[str]) -> str:
        res = self._aligner.msa(seqs, out_cons=True, out_msa=False)
        return str(res.cons_seq[0])


class _Spoa:
    name = "pyspoa"

    def __init__(self) -> None:
        self._poa = spoa.poa

    def consensus(self, seqs: list[str]) -> str:
        cons, _msa = self._poa(seqs, algorithm=1)  # global alignment
        return str(cons)


def get_backend(name: str) -> PoaBackend:
    """Return the named backend; raise ValueError for an unknown name."""
    factories: dict[str, Callable[[], PoaBackend]] = {"pyabpoa": _Abpoa, "pyspoa": _Spoa}
    if name not in factories:
        raise ValueError(f"unknown POA backend {name!r}")
    return factories[name]()
