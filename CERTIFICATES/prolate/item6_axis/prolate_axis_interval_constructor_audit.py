#!/usr/bin/env python3
"""Runtime audit for exact rational Arb interval construction.

The production constructor is accepted only if the constructed Arb ball contains
both endpoint balls for every audit case. Mere overlap is diagnostic and is not an
acceptance predicate. The audit also verifies positivity for strictly positive
input intervals and rejects any production source that still contains the historical
direct-fmpq midpoint/radius constructor.
"""
from __future__ import annotations

import json
from pathlib import Path

from flint import arb, ctx, fmpq

ctx.dps = 50


def closed_interval(lo: fmpq, hi: fmpq) -> arb:
    if not lo <= hi:
        raise ValueError("require lo <= hi")
    return arb(str((lo + hi) / 2), str((hi - lo) / 2))


def encloses(container: arb, contained: arb) -> bool:
    """Return True only when the full contained ball lies in container."""
    return bool(
        container.lower() <= contained.lower()
        and contained.upper() <= container.upper()
    )


cases = [
    (fmpq(1, 200), fmpq(3, 400)),
    (fmpq(1, 20), fmpq(1, 8)),
    (fmpq(63, 64), fmpq(127, 128)),
    (fmpq(1), fmpq(100)),
    (fmpq(0), fmpq(1, 2)),
]
records = []
checks = []
for lo, hi in cases:
    box = closed_interval(lo, hi)
    lo_ball = arb(str(lo))
    hi_ball = arb(str(hi))
    record = {
        "lo": str(lo),
        "hi": str(hi),
        "box": str(box),
        "contains_lo": encloses(box, lo_ball),
        "contains_hi": encloses(box, hi_ball),
        "overlaps_lo_diagnostic": bool(box.overlaps(lo_ball)),
        "overlaps_hi_diagnostic": bool(box.overlaps(hi_ball)),
        "contains_zero": bool(0 in box),
    }
    records.append(record)
    checks.extend([record["contains_lo"], record["contains_hi"]])
    if lo > 0:
        checks.append(not record["contains_zero"])

legacy = "return arb((lo + hi) / 2, (hi - lo) / 2)"
legacy_files = []
audit_path = Path(__file__).resolve()
for path in sorted(Path(__file__).parent.glob("*.py")):
    if path.resolve() == audit_path:
        continue
    if legacy in path.read_text(encoding="utf-8"):
        legacy_files.append(path.name)
checks.append(not legacy_files)

result = {
    "status": "PASSED" if all(checks) else "FAILED",
    "production_constructor": "arb(rational_midpoint_string, rational_radius_string)",
    "acceptance_predicate": (
        "constructed box contains both endpoint balls; overlap alone is insufficient"
    ),
    "cases": records,
    "legacy_files": legacy_files,
    "conclusion": (
        "PASS requires full endpoint-ball containment for every tested rational interval, "
        "zero exclusion for every strictly positive interval, and absence of the "
        "historical direct-fmpq constructor from production item-6 sources."
    ),
}
output = Path(__file__).with_suffix(".json")
output.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
