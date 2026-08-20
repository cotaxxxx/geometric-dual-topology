#!/usr/bin/env python3
"""Exact symbolic audit of the w->1^- axial derivative boundary profile."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


lam, c, w = sp.symbols("lam c w", positive=True, finite=True)
d = 1 - c
N = 1 - w * c
R2 = 1 - c**2 + lam**2 * (c - w) ** 2
S2 = 1 - c**2 + c**2 / lam**2
C2 = N**2 / (R2 * S2)

Rfac = 2 + (lam**2 - 1) * d
S2_boundary = d * (2 - d) + (1 - d) ** 2 / lam**2
C2_boundary = d / (Rfac * S2_boundary)
bracket_boundary = -c - lam**2 * d / Rfac

R2_at_boundary = sp.expand(R2.subs(w, 1))
N_at_boundary = sp.expand(N.subs(w, 1))
weighted_cw_bracket = -c + N * lam**2 * (c - w) / R2

checks = {
    "N_boundary_factor": sp.simplify(N_at_boundary - d) == 0,
    "R2_boundary_factor": sp.simplify(R2_at_boundary - d * Rfac) == 0,
    "S2_boundary_rewrite": sp.simplify(S2 - S2_boundary) == 0,
    "C_squared_boundary_formula": sp.simplify(C2.subs(w, 1) - C2_boundary) == 0,
    "weighted_Cw_boundary_bracket": sp.simplify(
        weighted_cw_bracket.subs(w, 1) - bracket_boundary
    ) == 0,
}

result = {
    "status": "PASSED" if all(checks.values()) else "FAILED",
    "checks": checks,
    "formulas": {
        "d": str(d),
        "R_factor": str(Rfac),
        "S2_boundary": str(S2_boundary),
        "C_boundary": "sqrt(d)/sqrt(R_factor*S2_boundary)",
        "weighted_Cw_bracket": str(bracket_boundary),
        "boundary_integrand": (
            "1/2*(-c*h(C_boundary)+h'(C_boundary)*C_boundary*"
            "weighted_Cw_bracket)"
        ),
    },
    "conclusion": (
        "The one-sided boundary value Psi_lambda(1^-) is represented by a "
        "regular one-dimensional integral with no division by 1-wc."
    ),
    "script_sha256": sha256_file(Path(__file__)),
}

output = Path("prolate_axis_pole_boundary_symbolic_audit.json")
output.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
