#!/usr/bin/env python3
"""Exact signed-angle audit for the scaled tail H and M kernels."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


mu, w = sp.symbols("mu w", positive=True, finite=True)
c = sp.symbols("c", real=True, finite=True)
rho = sp.sqrt(1 - c**2)
N = 1 - w * c
P = (c - w) ** 2 + mu**2 * (1 - c**2)
S = 1 - c**2 + mu**2 * c**2
Y = w - c + mu**2 * c
Z = c - w + mu**2 * c
delta = sp.atan(rho * Y / (mu * N))
delta_mu = sp.simplify(sp.diff(delta, mu))
delta_mu_claim = rho * N * Z / (P * S)
angle_floor = sp.pi**2 / 4

psi_integrand = -c * delta**2 / 2 + N * delta * mu * rho / P
H_integrand_raw = sp.simplify(psi_integrand / (mu * w))
H_integrand = -c * (delta**2 - angle_floor) / (2 * mu * w) + N * delta * rho / (w * P)
M_integrand = (
    c / (2 * w) * (
        2 * delta * delta_mu_claim - (delta**2 - angle_floor) / mu
    )
    - mu * N * rho * delta_mu_claim / (w * P)
    + 2 * mu**2 * N * rho**3 * delta / (w * P**2)
)

checks = {
    "mu_derivative_of_signed_angle": sp.simplify(
        delta_mu - delta_mu_claim
    ) == 0,
    "constant_angle_term_integrates_to_zero": sp.integrate(
        -c * angle_floor / (2 * mu * w), (c, -1, 1)
    ) == 0,
    "regularized_H_differs_only_by_zero_integral": sp.simplify(
        H_integrand_raw - H_integrand + c * angle_floor / (2 * mu * w)
    ) == 0,
    "M_is_negative_log_derivative": sp.simplify(
        -mu * sp.diff(H_integrand, mu) - M_integrand
    ) == 0,
    "PS_identity": sp.simplify(
        P * S - (mu**2 * N**2 + (1 - c**2) * Y**2)
    ) == 0,
}

result = {
    "status": "PASSED" if all(checks.values()) else "FAILED",
    "checks": checks,
    "domain_certificate": {
        "mu_positive": True,
        "w_positive_and_less_than_one_enforced_by_driver": True,
        "c_integrated_over_minus1_to_1": True,
        "deduction": (
            "For 0<w<1 and -1<=c<=1, N=1-w*c>0. Therefore the real atan "
            "branch is continuous on every certified positive-mu tail box."
        ),
    },
    "formulas": {
        "P": str(P),
        "S": str(S),
        "delta": str(delta),
        "delta_mu": str(delta_mu_claim),
        "H_integrand": str(H_integrand),
        "M_integrand": str(M_integrand),
    },
    "conclusion": (
        "The scaled tail kernels H=Psi/(mu*w) and M=-mu*dH/dmu are represented "
        "using only a signed atan and rational-algebraic factors. The constant "
        "pi^2/4 term is removed exactly using integral_{-1}^1 c dc=0."
    ),
    "warning": (
        "The formulas are exact. Uniform tail positivity still requires Arb "
        "coverage and a separate mu=0 transfer or remainder theorem."
    ),
    "script_sha256": sha256_file(Path(__file__)),
}

output = Path("prolate_axis_signed_tail_symbolic_audit.json")
output.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
