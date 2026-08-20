"""Independent primitive checks used by the record-layout verifier."""
from __future__ import annotations
from fractions import Fraction

from numeric_schema import D_ZERO, Dyadic


def _partition(start: Fraction, end: Fraction, width: Fraction) -> list[tuple[Fraction, Fraction]]:
    cells = []
    left = start
    while left < end:
        right = min(left + width, end)
        cells.append((left, right))
        left = right
    return cells


def _positive_width(record: dict) -> bool:
    try:
        return D_ZERO < Dyadic.from_json(record["width"], "join.width")
    except (KeyError, ValueError):
        return False


__all__ = [name for name in globals() if not name.startswith("__")]
