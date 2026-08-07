#!/usr/bin/env python3
"""Exact symbolic audit of the cancellation-regularized tail integrand.

The scaled tail quantity is

    H(mu,w) = Psi_{1/mu}(w)/(mu*w).

Before interval evaluation, the constant h(0)=pi^2/4 is removed from h(C)
using integral_{-1}^1 c dc = 0. This eliminates the apparent 1/mu term away
from the moving layer and leaves the logarithmic layer as the only singular
asymptotic contribution.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


mu, w, c = sp.symbols("mu w c", positive=True, finite=True)
h0 = sp.pi**2 / 4
N = 1 - w * c
P = (c - w) ** 2 + mu**2 * (1 - c**2)
S = 1 - c**2 + mu**2 * c**2
C = mu * N / sp.sqrt(P * S)

h = sp.Function("h")
h1 = sp.Function("h1")
hbar = (h(C) - h0) / C

Cw_claim = C * (-c / N + (c - w) / P)
original_scaled = (
    -c * h(C) + N * h1(C) * Cw_claim
) / (2 * mu * w)
regularized = N / (2 * w * sp.sqrt(P * S)) * (
    -c * hbar + h1(C) * (-c + N * (c - w) / P)
)

difference = sp.simplify(original_scaled - regularized)
expected_difference = -c * h0 / (2 * mu * w)

# Check equivalence with the original lambda geometry after lambda=1/mu.
lam = sp.symbols("lam", positive=True, finite=True)
R2_lam = 1 - c**2 + lam**2 * (c - w) ** 2
S2_lam = 1 - c**2 + c**2 / lam**2
C_lam = N / sp.sqrt(R2_lam * S2_lam)
C_mu_from_lam = sp.simplify(C_lam.subs(lam, 1 / mu))

checks = {
    "mu_geometry_matches_lambda_geometry": sp.simplify(C_mu_from_lam - C) == 0,
    "C_w_formula_in_mu_variables": sp.simplify(sp.diff(C, w) - Cw_claim) == 0,
    "regularized_integrand_difference": sp.simplify(
        difference - expected_difference
    ) == 0,
    "discarded_term_integrates_to_zero": sp.integrate(c, (c, -1, 1)) == 0,
}

result = {
    "status": "PASSED" if all(checks.values()) else "FAILED",
    "checks": checks,
    "formulas": {
        "N": str(N),
        "P": str(P),
        "S": str(S),
        "C": str(C),
        "C_w": str(Cw_claim),
        "h_bar": "(h(C)-pi^2/4)/C",
        "regularized_scaled_integrand": str(regularized),
        "discarded_zero-integral_term": str(expected_difference),
    },
    "conclusion": (
        "The integral of the regularized kernel equals H(mu,w). The apparent "
        "1/mu constant-angle term is removed exactly before interval arithmetic."
    ),
    "limitations": (
        "This audit does not yet subtract the logarithmic main term or bound the "
        "finite remainder Bhat uniformly."
    ),
    "script_sha256": sha256_file(Path(__file__)),
}

output = Path("prolate_axis_tail_regularized_symbolic_audit.json")
output.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
