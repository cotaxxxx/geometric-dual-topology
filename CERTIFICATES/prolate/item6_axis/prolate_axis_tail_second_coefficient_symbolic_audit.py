#!/usr/bin/env python3
"""Exact logarithmic coefficient audit for the scaled axial second derivative."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


w = sp.symbols("w", real=True, finite=True)
pi = sp.pi
H_coefficient = 3 * pi * sp.sqrt(1 - w**2)
Psi_coefficient = sp.simplify(w * H_coefficient)
second_coefficient = sp.simplify(sp.diff(Psi_coefficient, w))
claimed = 3 * pi * (1 - 2 * w**2) / sp.sqrt(1 - w**2)

center_edge = sp.simplify(claimed.subs(w, sp.Rational(1, 20)))
pole_edge = sp.simplify(claimed.subs(w, sp.Rational(3, 4)))

checks = {
    "differentiate_scaled_profile_coefficient": sp.simplify(
        second_coefficient - claimed
    ) == 0,
    "center_coefficient_positive": bool(center_edge > 0),
    "pole_coefficient_negative": bool(pole_edge < 0),
    "sign_change_location": sp.solve(sp.factor(1 - 2 * w**2), w) == [
        -sp.sqrt(2) / 2,
        sp.sqrt(2) / 2,
    ],
}

result = {
    "status": "PASSED" if all(checks.values()) else "FAILED",
    "checks": checks,
    "formulas": {
        "H_log_coefficient": str(H_coefficient),
        "Psi_log_coefficient": str(Psi_coefficient),
        "A_second_over_mu_log_coefficient": str(claimed),
        "center_edge_value": str(center_edge),
        "pole_edge_value": str(pole_edge),
    },
    "conclusion": (
        "The logarithmic coefficient of A_lambda''(w)/mu is positive on "
        "0<=w<1/sqrt(2) and negative on 1/sqrt(2)<w<1. In particular it is "
        "uniformly positive on the center tail w<=1/20 and uniformly negative "
        "on the pole tail w>=3/4."
    ),
    "warning": (
        "The coefficient sign alone is not the tail theorem; a uniform bounded "
        "remainder certificate is still required."
    ),
    "script_sha256": sha256_file(Path(__file__)),
}

output = Path("prolate_axis_tail_second_coefficient_symbolic_audit.json")
output.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
