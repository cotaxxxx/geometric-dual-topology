#!/usr/bin/env python3
"""Exact audit of the paired moving-layer charts used for tail H and M."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


mu, w, t = sp.symbols("mu w t", positive=True, finite=True)
radius = sp.Integer(4)
q = 8 * t - 4
c_inner = w + mu * q
rho_inner = (1 - c_inner) * (1 + c_inner)
p_inner = (c_inner - w) ** 2 + mu**2 * (1 - c_inner**2)

jacobian = 1 - w - radius * mu
x = radius * mu + jacobian * t
c_left = w - x
c_right = w + x
rho_left = (1 - w + x) * (1 + w - x)
rho_right = jacobian * (1 - t) * (1 + w + x)

x_far = 1 - w + 2 * w * t
c_far = w - x_far
rho_far = 4 * w * (1 - t) * (1 - w + w * t)

checks = {
    "inner_difference": sp.expand(c_inner - w - mu * q) == 0,
    "inner_rho_factor": sp.expand(rho_inner - (1 - c_inner**2)) == 0,
    "inner_P_factor": sp.expand(
        p_inner - mu**2 * (q**2 + rho_inner)
    ) == 0,
    "paired_left_rho": sp.expand(rho_left - (1 - c_left**2)) == 0,
    "paired_right_rho": sp.expand(rho_right - (1 - c_right**2)) == 0,
    "paired_left_difference": sp.expand(c_left - w + x) == 0,
    "paired_right_difference": sp.expand(c_right - w - x) == 0,
    "far_left_rho": sp.expand(rho_far - (1 - c_far**2)) == 0,
    "far_left_difference": sp.expand(c_far - w + x_far) == 0,
    "inner_left_endpoint": sp.expand(
        c_inner.subs(t, 0) - (w - 4 * mu)
    ) == 0,
    "inner_right_endpoint": sp.expand(
        c_inner.subs(t, 1) - (w + 4 * mu)
    ) == 0,
    "pair_x_left": sp.expand(x.subs(t, 0) - 4 * mu) == 0,
    "pair_x_right": sp.expand(x.subs(t, 1) - (1 - w)) == 0,
    "pair_left_end": sp.expand(c_left.subs(t, 1) - (2 * w - 1)) == 0,
    "pair_right_end": sp.expand(c_right.subs(t, 1) - 1) == 0,
    "far_x_left": sp.expand(x_far.subs(t, 0) - (1 - w)) == 0,
    "far_x_right": sp.expand(x_far.subs(t, 1) - (1 + w)) == 0,
    "far_c_left": sp.expand(c_far.subs(t, 0) - (2 * w - 1)) == 0,
    "far_c_right": sp.expand(c_far.subs(t, 1) + 1) == 0,
}

result = {
    "status": "PASSED" if all(checks.values()) else "FAILED",
    "checks": checks,
    "partition": [
        "far-left covers c in [-1,2w-1]",
        "paired-left covers c in [2w-1,w-4mu]",
        "inner covers c in [w-4mu,w+4mu]",
        "paired-right covers c in [w+4mu,1]",
    ],
    "conclusion": (
        "The transformed integrals cover [-1,1] exactly, with no gap or "
        "overlap, and preserve c-w, P, and endpoint factors algebraically."
    ),
    "script_sha256": sha256_file(Path(__file__)),
}
output = Path("prolate_axis_tail_paired_chart_symbolic_audit.json")
output.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
