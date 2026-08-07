#!/usr/bin/env python3
"""Exact audit for the pole-modulus blow-up charts and cancellations."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


u, y, z, r, L = sp.symbols("u y z r L", positive=True, finite=True)
n, radial, s2 = sp.symbols("n radial s2", positive=True, finite=True)
c, q = sp.symbols("c q", real=True, finite=True)
h1, h2 = sp.symbols("h1 h2", finite=True)
root = sp.sqrt(radial * s2)
C = n / root
gbar = -c / n + L * q / radial
gvbar = -c**2 / n**2 - L / radial + 2 * L**2 * q**2 / radial**2

a = L * q / radial
cw_safe = (-c + n * a) / root
cww_safe = (-2 * c * a + n * (3 * a**2 - L / radial)) / root
compact_safe = -2 * c * h1 * cw_safe + n * (h2 * cw_safe**2 + h1 * cww_safe)
compact_raw = -2 * c * h1 * C * gbar + n * (
    h2 * (C * gbar) ** 2 + h1 * C * (gbar**2 + gvbar)
)

d_inner = u**2 * y**2
c_inner = 1 - d_inner
w_inner = 1 - u
N_inner = 1 - w_inner * c_inner
R2_inner = 1 - c_inner**2 + L * (c_inner - w_inner) ** 2
n_inner = 1 + u * y**2 - u**2 * y**2
radial_inner = y**2 * (2 - u**2 * y**2) + L * (1 - u * y**2) ** 2

d_outer = z**2
c_outer = 1 - d_outer
u_outer = r * z
w_outer = 1 - u_outer
N_outer = 1 - w_outer * c_outer
R2_outer = 1 - c_outer**2 + L * (c_outer - w_outer) ** 2
n_outer = r + z - r * z**2
radial_outer = 2 - z**2 + L * (r - z) ** 2

checks = {
    "Cw_product_cancellation": sp.simplify(C * gbar - cw_safe) == 0,
    "Cww_product_cancellation": sp.simplify(
        C * (gbar**2 + gvbar) - cww_safe
    ) == 0,
    "compact_density_assembly": sp.simplify(compact_raw - compact_safe) == 0,
    "inner_N_factor": sp.expand(N_inner - u * n_inner) == 0,
    "inner_R2_factor": sp.expand(R2_inner - u**2 * radial_inner) == 0,
    "inner_c_minus_w_factor": sp.expand(
        c_inner - w_inner - u * (1 - u * y**2)
    ) == 0,
    "inner_half_jacobian": sp.simplify(
        sp.Rational(1, 2) * sp.diff(d_inner, y) - u**2 * y
    ) == 0,
    "outer_N_factor": sp.expand(N_outer - z * n_outer) == 0,
    "outer_R2_factor": sp.expand(R2_outer - z**2 * radial_outer) == 0,
    "outer_c_minus_w_factor": sp.expand(
        c_outer - w_outer - z * (r - z)
    ) == 0,
    "outer_half_jacobian": sp.simplify(
        sp.Rational(1, 2) * sp.diff(d_outer, z) - z
    ) == 0,
}

result = {
    "status": "PASSED" if all(checks.values()) else "FAILED",
    "checks": checks,
    "formulas": {
        "Cw_safe": str(cw_safe),
        "Cww_safe": str(cww_safe),
        "compact_density": str(compact_safe),
        "inner_N": str(u * n_inner),
        "inner_R2": str(u**2 * radial_inner),
        "outer_N": str(z * n_outer),
        "outer_R2": str(z**2 * radial_outer),
    },
    "conclusion": (
        "The identities are checked on N!=0. Their cancellation-safe right-hand "
        "sides extend continuously to the projective corner N=0 and are the "
        "expressions evaluated by the Arb cover. The inner and outer Jacobians "
        "give exactly the transformed densities used by the modulus proof."
    ),
    "branch_note": (
        "On 0<=u<=1/64, 0<=y,r<=1, 0<=z<=sqrt(2), and L>=1, "
        "the geometric radicands are nonnegative; the certificate evaluates "
        "the positive square-root branch inherited from the original geometry."
    ),
    "script_sha256": sha256_file(Path(__file__)),
}

output = Path("prolate_axis_pole_modulus_symbolic_audit.json")
output.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
