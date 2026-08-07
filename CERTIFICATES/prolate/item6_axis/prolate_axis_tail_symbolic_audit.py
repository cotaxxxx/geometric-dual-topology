#!/usr/bin/env python3
"""Formal outer/Laurent audit for the prolate item 6 aspect-ratio tail.

This verifies the logarithmic coefficient in the large-lambda expansion. It
is an exact symbolic audit of the outer kernel; it is not a uniform remainder
bound and therefore is not by itself an interval certificate.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


w, d, z = sp.symbols("w d z", positive=True, finite=True)
pi = sp.pi

a = 1 - w**2
c = w + d
N = 1 - w * c
G = 1 - c**2

# h(z)=acos(z)^2. The first outer correction uses h'(0)=-pi.
h = sp.acos(z) ** 2
h0 = sp.simplify(h.subs(z, 0))
h1_0 = sp.simplify(sp.diff(h, z).subs(z, 0))

# For mu=1/lambda and fixed d=c-w != 0,
#   C = mu*C1 + O(mu^3),
# with different analytic representatives on the two sides of d=0.
C1_plus = N / (d * sp.sqrt(G))
C1_minus = -N / (d * sp.sqrt(G))

# The corresponding first coefficient of C_w is
#   sqrt(G)/(d*abs(d)).
Cw1_plus = sp.sqrt(G) / d**2
Cw1_minus = -sp.sqrt(G) / d**2

# First outer coefficient K1 in the axial derivative kernel
# K = 1/2(-c h(C) + N h'(C) C_w).
K1_plus = sp.simplify(
    sp.Rational(1, 2) * (-c * h1_0 * C1_plus + N * h1_0 * Cw1_plus)
)
K1_minus = sp.simplify(
    sp.Rational(1, 2) * (-c * h1_0 * C1_minus + N * h1_0 * Cw1_minus)
)

# The d^-2 terms are opposite and cancel under the two-sided matching. The
# d^-1 residues produce the logarithm.
lead_plus = sp.simplify(sp.limit(d**2 * K1_plus, d, 0, dir="+"))
lead_minus = sp.simplify(sp.limit(d**2 * K1_minus, d, 0, dir="-"))
residue_plus = sp.simplify(sp.residue(K1_plus, d, 0))
residue_minus = sp.simplify(sp.residue(K1_minus, d, 0))

expected_residue = sp.Rational(3, 2) * pi * w * sp.sqrt(a)
physical_log_coefficient = sp.simplify(residue_plus - residue_minus)
expected_log_coefficient = 3 * pi * w * sp.sqrt(a)
quotient_coefficient = sp.simplify(physical_log_coefficient / w)
expected_quotient_coefficient = 3 * pi * sp.sqrt(a)
center_tail_coefficient = sp.simplify(
    sp.limit(quotient_coefficient, w, 0, dir="+")
)

checks = {
    "h_at_zero": sp.simplify(h0 - pi**2 / 4) == 0,
    "h_prime_at_zero": sp.simplify(h1_0 + pi) == 0,
    "opposite_d_minus_2_coefficients": sp.simplify(lead_plus + lead_minus) == 0,
    "right_residue": sp.simplify(residue_plus - expected_residue) == 0,
    "left_residue": sp.simplify(residue_minus + expected_residue) == 0,
    "two_sided_log_coefficient": sp.simplify(
        physical_log_coefficient - expected_log_coefficient
    ) == 0,
    "odd_factor_quotient_coefficient": sp.simplify(
        quotient_coefficient - expected_quotient_coefficient
    ) == 0,
    "center_tail_corner": sp.simplify(center_tail_coefficient - 3 * pi) == 0,
}

result = {
    "status": "PASSED" if all(checks.values()) else "FAILED",
    "scope": "formal outer/Laurent coefficient audit; no uniform remainder bound",
    "checks": checks,
    "coefficients": {
        "d_minus_2_right": str(lead_plus),
        "d_minus_2_left": str(lead_minus),
        "d_minus_1_right": str(residue_plus),
        "d_minus_1_left": str(residue_minus),
        "psi_log_coefficient": str(physical_log_coefficient),
        "scaled_odd_quotient_log_coefficient": str(quotient_coefficient),
        "center_tail_limit": str(center_tail_coefficient),
    },
    "asymptotic_form": (
        "Psi_{1/mu}(w) = mu*w*(3*pi*sqrt(1-w^2)*log(1/mu) "
        "+ B(mu,w^2))"
    ),
    "warning": (
        "The audit identifies the logarithmic coefficient exactly but does not "
        "prove that B is uniformly bounded below on a tail strip."
    ),
    "script_sha256": sha256_file(Path(__file__)),
}

out = Path("prolate_axis_tail_symbolic_audit.json")
out.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["status"] == "PASSED" else 1)
