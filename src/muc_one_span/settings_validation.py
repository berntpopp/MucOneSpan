"""Shared validation helpers for the runtime settings dataclasses.

These small, private-by-convention checks raise ``ValueError`` with a
consistent, field-qualified message so ``__post_init__`` bodies across
``muc_one_span.settings`` (and split-out modules such as
``muc_one_span.settings_hybrid``) stay terse and uniform.
"""

from __future__ import annotations

import math


def _integer(name: str, value: object, minimum: int = 0) -> None:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def _number(name: str, value: object, minimum: float = 0, maximum: float | None = None) -> None:
    try:
        finite = isinstance(value, (int, float)) and math.isfinite(value)
    except OverflowError:
        finite = False
    if (
        not isinstance(value, (int, float))
        or isinstance(value, bool)
        or not finite
        or value < minimum
        or (maximum is not None and value > maximum)
    ):
        limit = f" in [{minimum}, {maximum}]" if maximum is not None else f" >= {minimum}"
        raise ValueError(f"{name} must be a finite number{limit}")


def _boolean(name: str, value: object) -> None:
    if type(value) is not bool:
        raise ValueError(f"{name} must be a boolean")


def _string(name: str, value: object, *, optional: bool = False, empty: bool = False) -> None:
    if optional and value is None:
        return
    if not isinstance(value, str) or (not value.strip() and not (empty and value == "")):
        raise ValueError(f"{name} must be a nonempty string" + (" or null" if optional else ""))


def _open_unit_interval(name: str, value: object) -> None:
    _number(name, value, 0, 1)
    if value in (0, 1):
        raise ValueError(f"{name} must be strictly between 0 and 1")


def _choice(name: str, value: object, allowed: tuple[str, ...]) -> None:
    if not isinstance(value, str) or value not in allowed:
        raise ValueError(f"{name} must be one of {allowed!r}")
