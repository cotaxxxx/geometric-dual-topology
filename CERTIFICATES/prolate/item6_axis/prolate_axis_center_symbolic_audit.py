#!/usr/bin/env python3
"""Exact symbolic audit of the item 6 center Hessian kernel."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


lam, w, c, z = sp.symbols("lam w c z", positive=True, finite=True)

L = 1 + (lam**2 - 1) * c**2
W2 = lam**2 * (1 - c**2) + c**2
sqrtL = sp.sqrt(L)
sqrtW = sp.sqrt(W2)

C0 = lam / (sqrtL * sqrtW)
g0 = c * (lam**2 - 1) * (1 - c**2) / L
gw0 = -c**2 - lam**2 / L + 2 * lam**4 * c**2 / L**2

Cw0 = lam * c * (lam**2 - 1) * (1 - c**2) / (L * sqrtL * sqrtW)
Cww0 = (
    lam**3
    * (1 - c**2)
    * (2 * (lam**2 - 1) * c**2 - 1)
    / (L**2 * sqrtL * sqrtW)
)

h = sp.acos(z) ** 2
h1 = -2 * sp.acos(z) / sp.sqrt(1 - z**2)
h2 = 2 / (1 - z**2) - 2 * z * sp.acos(z) / (1 - z**2) ** sp.Rational(3, 2)

h1_at_one = sp.simplify(sp.limit(h1, z, 1, dir="-"))
h2_at_one = sp.simplify(sp.limit(h2, z, 1, dir="-"))

Q_kernel = sp.Rational(1, 2) * (
    -2 * c * h1.subs(z, C0) * Cw0
    + h2.subs(z, C0) * Cw0**2
    + h1.subs(z, C0) * Cww0
)

sphere_Cw0 = sp.simplify(Cw0.subs(lam, 1))
sphere_Cww0 = sp.simplify(Cww0.subs(lam, 1))
sphere_kernel = sp.simplify(
    sp.Rational(1, 2)
    * (
        -2 * c * h1_at_one * sphere_Cw0
        + h2_at_one * sphere_Cw0**2
        + h1_at_one * sphere_Cww0
    )
)
sphere_integral = sp.simplify(sp.integrate(sphere_kernel, (c, -1, 1)))

checks = {
    "h_prime_formula": sp.simplify(sp.diff(h, z) - h1) == 0,
    "h_second_formula": sp.simplify(sp.diff(h, z, 2) - h2) == 0,
    "Cw_center_from_log_derivative": sp.simplify(C0 * g0 - Cw0) == 0,
    "Cww_center_from_log_derivative": sp.simplify(
        C0 * (g0**2 + gw0) - Cww0
    ) == 0,
    "h_prime_removable_limit": sp.simplify(h1_at_one + 2) == 0,
    "h_second_removable_limit": sp.simplify(h2_at_one - sp.Rational(2, 3)) == 0,
    "sphere_center_kernel": sp.simplify(sphere_kernel - (1 - c**2)) == 0,
    "sphere_axial_hessian": sp.simplify(sphere_integral - sp.Rational(4, 3)) == 0,
}

result = {
    "status": "PASSED" if all(checks.values()) else "FAILED",
    "scope": "exact center-kernel identities and sphere endpoint; no positivity proof for all lambda",
    "checks": checks,
    "formulas": {
        "L": str(L),
        "W2": str(W2),
        "C0": str(C0),
        "Cw0": str(Cw0),
        "Cww0": str(Cww0),
        "Q_parallel_kernel": str(Q_kernel),
        "h_prime_at_one": str(h1_at_one),
        "h_second_at_one": str(h2_at_one),
        "sphere_kernel": str(sphere_kernel),
        "sphere_Q_parallel": str(sphere_integral),
    },
    "warning": (
        "The script freezes the exact integrand for Q_parallel(lambda). It does "
        "not certify Q_parallel(lambda)>0 on an interval of lambda."
    ),
    "script_sha256": sha256_file(Path(__file__)),
}

out = Path("prolate_axis_center_symbolic_audit.json")
out.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
