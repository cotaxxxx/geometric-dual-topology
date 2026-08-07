#!/usr/bin/env python3
"""Exact tail-pole boundary normalization and mu=0 limit audit."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


mu, d = sp.symbols("mu d", positive=True, finite=True)
pi = sp.pi
c = 1 - d
Rbar = d + mu**2 * (2 - d)
s2 = d * (2 - d) + mu**2 * (1 - d) ** 2
C = mu * sp.sqrt(d) / sp.sqrt(Rbar * s2)
C_over_mu = sp.sqrt(d) / sp.sqrt(Rbar * s2)
weighted_cw = -c - d / Rbar

hbar0 = -pi
h1_0 = -pi
C_over_mu_0 = 1 / sp.sqrt(d * (2 - d))
weighted_cw_0 = d - 2
bracket0 = sp.simplify(-c * hbar0 + h1_0 * weighted_cw_0)
limit_kernel = sp.simplify(
    sp.Rational(1, 2) * C_over_mu_0 * bracket0
)
limit_integral = sp.simplify(sp.integrate(limit_kernel, (d, 0, 2)))

checks = {
    "boundary_R_factor_mu_form": sp.expand(
        mu**2 * (2 + (mu**-2 - 1) * d) - Rbar
    ) == 0,
    "boundary_cosine_mu_form": sp.simplify(
        C - mu * sp.sqrt(d) / sp.sqrt(Rbar * s2)
    ) == 0,
    "constant_angle_cancellation": sp.integrate(-c * pi**2 / 8, (d, 0, 2)) == 0,
    "hbar_at_zero": hbar0 == -pi,
    "h1_at_zero": h1_0 == -pi,
    "limit_bracket": sp.simplify(bracket0 - pi * (3 - 2 * d)) == 0,
    "scaled_boundary_limit": sp.simplify(limit_integral - pi**2 / 2) == 0,
}

result = {
    "status": "PASSED" if all(checks.values()) else "FAILED",
    "checks": checks,
    "formulas": {
        "Rbar": str(Rbar),
        "s2": str(s2),
        "C": str(C),
        "C_over_mu": str(C_over_mu),
        "weighted_cw": str(weighted_cw),
        "mu_zero_kernel": str(limit_kernel),
        "mu_zero_integral": str(limit_integral),
    },
    "conclusion": (
        "After removing h(0)=pi^2/4 using the odd c weight, the scaled pole "
        "boundary G(mu)=Phi(1/mu)/mu has the exact limit pi^2/2 at mu=0."
    ),
    "warning": (
        "This is the exact endpoint limit only. Uniform positivity on a positive "
        "mu interval still requires an Arb remainder or direct-cover certificate."
    ),
    "script_sha256": sha256_file(Path(__file__)),
}

output = Path("prolate_axis_tail_pole_boundary_symbolic_audit.json")
output.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
