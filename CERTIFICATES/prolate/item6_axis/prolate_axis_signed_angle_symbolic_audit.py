#!/usr/bin/env python3
"""Exact signed-angle reduction audit for the prolate axial profile."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


lam = sp.symbols("lam", positive=True, finite=True)
w, c = sp.symbols("w c", real=True, finite=True)
rho = sp.sqrt(1 - c**2)
N = 1 - w * c
X = rho * (lam * w - (lam - 1 / lam) * c)
r2 = 1 - c**2 + lam**2 * (c - w) ** 2
s2 = 1 - c**2 + c**2 / lam**2
T = X / N
delta = sp.atan(T)
delta_w = sp.simplify(sp.diff(delta, w))
delta_ww = sp.simplify(sp.diff(delta_w, w))
claimed_w = lam * rho / r2
claimed_ww = 2 * lam**3 * rho * (c - w) / r2**2

D = sp.Function("D")
profile_integrand = -sp.Rational(1, 2) * c * D(w) ** 2 + N * D(w) * sp.diff(D(w), w)
second_integrand = sp.simplify(sp.diff(profile_integrand, w))
claimed_second_integrand = (
    -2 * c * D(w) * sp.diff(D(w), w)
    + N * (
        sp.diff(D(w), w) ** 2
        + D(w) * sp.diff(D(w), w, 2)
    )
)

checks = {
    "gram_determinant_identity": sp.simplify(N**2 + X**2 - r2 * s2) == 0,
    "signed_angle_first_derivative": sp.simplify(delta_w - claimed_w) == 0,
    "signed_angle_second_derivative": sp.simplify(delta_ww - claimed_ww) == 0,
    "profile_second_derivative_assembly": sp.simplify(
        second_integrand - claimed_second_integrand
    ) == 0,
}

domain_certificate = {
    "driver_enforces_lambda_positive": True,
    "driver_enforces_0_le_w_lt_1": True,
    "integration_domain_is_minus1_le_c_le_1": True,
    "deduction": (
        "Because 0<=w<1 and -1<=c<=1, w*c<=w<1, hence N=1-w*c>0. "
        "This is an exact order argument enforced by the Arb driver input checks, "
        "not a symbolic polynomial identity."
    ),
}

result = {
    "status": "PASSED" if all(checks.values()) else "FAILED",
    "checks": checks,
    "domain_certificate": domain_certificate,
    "formulas": {
        "N": str(N),
        "signed_cross": str(X),
        "r2": str(r2),
        "s2": str(s2),
        "delta": str(delta),
        "delta_w": str(claimed_w),
        "delta_ww": str(claimed_ww),
        "Psi_integrand": str(profile_integrand),
        "A_second_integrand": str(claimed_second_integrand),
    },
    "branch_statement": (
        "For 0<=w<1 and |c|<=1, N=1-w*c>0. Hence delta=atan(X/N) lies "
        "in (-pi/2,pi/2), cos(delta)=N/sqrt(N^2+X^2)=C, and "
        "acos(C)^2=delta^2."
    ),
    "conclusion": (
        "The original radial-normal angle square and both axial derivatives can "
        "be evaluated using a signed atan and rational-algebraic derivatives, "
        "without hypergeometric angle regularization."
    ),
    "script_sha256": sha256_file(Path(__file__)),
}

output = Path("prolate_axis_signed_angle_symbolic_audit.json")
output.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
