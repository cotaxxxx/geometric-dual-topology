#!/usr/bin/env python3
"""Exact symbolic audit for the finite item 6 center-cap kernel."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


lam, w, c, z = sp.symbols("lam w c z", positive=True, finite=True)
H = sp.Function("H")

N = 1 - w * c
R2 = 1 - c**2 + lam**2 * (c - w) ** 2
S2 = 1 - c**2 + c**2 / lam**2
C = N / sp.sqrt(R2 * S2)

g = -c / N + lam**2 * (c - w) / R2
g_w = -c**2 / N**2 - lam**2 / R2 + 2 * lam**4 * (c - w) ** 2 / R2**2
C_w = C * g
C_ww = C * (g**2 + g_w)

h1 = sp.diff(H(z), z).subs(z, C)
h2 = sp.diff(H(z), z, 2).subs(z, C)
second_kernel_claim = sp.Rational(1, 2) * (
    -2 * c * h1 * C_w + N * (h2 * C_w**2 + h1 * C_ww)
)
second_kernel_direct = sp.diff(sp.Rational(1, 2) * N * H(C), w, 2)

checks = {
    "C_w_formula": sp.simplify(sp.diff(C, w) - C_w) == 0,
    "g_w_formula": sp.simplify(sp.diff(g, w) - g_w) == 0,
    "C_ww_formula": sp.simplify(sp.diff(C, w, 2) - C_ww) == 0,
    "second_derivative_kernel": sp.simplify(
        second_kernel_direct - second_kernel_claim
    ) == 0,
}

result = {
    "status": "PASSED" if all(checks.values()) else "FAILED",
    "checks": checks,
    "formulas": {
        "N": str(N),
        "R2": str(R2),
        "S2": str(S2),
        "C": str(C),
        "g": str(g),
        "g_w": str(g_w),
        "C_w": str(C_w),
        "C_ww": str(C_ww),
        "A_second_kernel": str(second_kernel_claim),
    },
    "consequence": (
        "Uniform positivity of this kernel integral on 0<=v<=w0 implies "
        "Psi_lambda(w)>0 on 0<w<=w0 by Psi(w)/w=int_0^1 A''(tw) dt."
    ),
    "script_sha256": sha256_file(Path(__file__)),
}

output = Path("prolate_axis_center_cap_symbolic_audit.json")
output.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
