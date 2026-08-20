#!/usr/bin/env python3
"""Exact symbolic audit for the prolate item 6 axial reduction."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


lam, w, c, t = sp.symbols("lam w c t", positive=True, finite=True)

N = 1 - w * c
R2 = 1 - c**2 + lam**2 * (c - w)**2
S2 = 1 - c**2 + c**2 / lam**2
C = N / sp.sqrt(R2 * S2)

Cw_claim = C * (-c / N + lam**2 * (c - w) / R2)

# Reflection is checked using a separate unconstrained symbol substitution.
reflection = sp.simplify(C.subs({c: -c, w: -w}) - C)
weight_reflection = sp.simplify(N.subs({c: -c, w: -w}) - N)

# Correlation-preserving moving-layer chart used by Region I.
c_left = -1 + (1 + w) * t
c_right = w + (1 - w) * t

checks = {
    "C_w_formula": sp.simplify(sp.diff(C, w) - Cw_claim) == 0,
    "R2_reflection": sp.simplify(R2.subs({c: -c, w: -w}) - R2) == 0,
    "S2_even_in_c": sp.simplify(S2.subs(c, -c) - S2) == 0,
    "C_reflection": reflection == 0,
    "cone_weight_reflection": weight_reflection == 0,
    "left_chart_starts_at_minus_one": sp.simplify(c_left.subs(t, 0) + 1) == 0,
    "left_chart_ends_at_moving_layer": sp.simplify(c_left.subs(t, 1) - w) == 0,
    "right_chart_starts_at_moving_layer": sp.simplify(c_right.subs(t, 0) - w) == 0,
    "right_chart_ends_at_one": sp.simplify(c_right.subs(t, 1) - 1) == 0,
    "left_chart_jacobian": sp.simplify(sp.diff(c_left, t) - (1 + w)) == 0,
    "right_chart_jacobian": sp.simplify(sp.diff(c_right, t) - (1 - w)) == 0,
}

result = {
    "status": "PASSED" if all(checks.values()) else "FAILED",
    "checks": checks,
    "formulas": {
        "N": str(N),
        "R2": str(R2),
        "S2": str(S2),
        "C": str(C),
        "C_w": str(Cw_claim),
        "c_left": str(c_left),
        "dc_left_dt": str(sp.diff(c_left, t)),
        "c_right": str(c_right),
        "dc_right_dt": str(sp.diff(c_right, t)),
    },
    "consequences": {
        "energy_even": "A_lambda(-w)=A_lambda(w), after c -> -c",
        "derivative_odd": "Psi_lambda(-w)=-Psi_lambda(w)",
        "interior_chart": (
            "[-1,w] and [w,1] are each mapped to t in [0,1] while preserving "
            "the shared parameter w inside c-w, 1-wc, and R2"
        ),
    },
    "script_sha256": sha256_file(Path(__file__)),
}

out = Path("prolate_axis_symbolic_audit.json")
out.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
