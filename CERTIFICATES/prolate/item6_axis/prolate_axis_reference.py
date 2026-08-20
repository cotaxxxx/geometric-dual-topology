#!/usr/bin/env python3
"""Non-certified high-precision scout for the prolate axial profile.

This script freezes the exact 1D formulas used to design item 6. It uses
mpmath and is not an interval proof.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def angle_data(c: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    """Return h(c)=acos(c)^2 and its regularized first derivative."""
    if c > 1 and c - 1 < 100 * mp.eps:
        c = mp.mpf(1)
    if c < -1 and -1 - c < 100 * mp.eps:
        c = mp.mpf(-1)
    u = mp.acos(c)
    if abs(u) < mp.sqrt(mp.eps):
        # u/sin(u) = 1 + u^2/6 + 7u^4/360 + O(u^6)
        ratio = 1 + u * u / 6 + 7 * u**4 / 360
    else:
        ratio = u / mp.sin(u)
    return u * u, -2 * ratio


def geometry(lam: mp.mpf, w: mp.mpf, c: mp.mpf) -> tuple[mp.mpf, ...]:
    n = 1 - w * c
    r2 = 1 - c * c + lam * lam * (c - w) ** 2
    s2 = 1 - c * c + c * c / (lam * lam)
    cosine = n / mp.sqrt(r2 * s2)
    if cosine > 1 and cosine - 1 < 100 * mp.eps:
        cosine = mp.mpf(1)
    return n, r2, s2, cosine


def energy(lam: mp.mpf, w: mp.mpf) -> mp.mpf:
    def kernel(c: mp.mpf) -> mp.mpf:
        n, _, _, cosine = geometry(lam, w, c)
        h, _ = angle_data(cosine)
        return mp.mpf("0.5") * n * h

    # c=w is the narrow large-aspect-ratio layer. Splitting there improves
    # the non-certified quadrature without changing the exact formula.
    return mp.quad(kernel, [-1, w, 1])


def psi(lam: mp.mpf, w: mp.mpf) -> mp.mpf:
    def kernel(c: mp.mpf) -> mp.mpf:
        n, r2, _, cosine = geometry(lam, w, c)
        h, h1 = angle_data(cosine)
        cosine_w = cosine * (-c / n + lam * lam * (c - w) / r2)
        return mp.mpf("0.5") * (-c * h + n * h1 * cosine_w)

    # Split at the moving layer c=w. This is essential for reliable tail
    # scouting when lambda is large; it remains non-interval quadrature.
    return mp.quad(kernel, [-1, w, 1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--json", default="prolate_axis_reference.json")
    args = parser.parse_args()
    mp.mp.dps = args.dps

    lambdas = [
        mp.mpf("1"),
        mp.mpf("2"),
        mp.mpf("2.0653823"),
        mp.mpf("3.4348684"),
        mp.mpf("4.7243834"),
        mp.mpf("10"),
        mp.mpf("100"),
    ]
    ws = [
        mp.mpf("0.01"),
        mp.mpf("0.1"),
        mp.mpf("0.3"),
        mp.mpf("0.5"),
        mp.mpf("0.7"),
        mp.mpf("0.9"),
        mp.mpf("0.99"),
    ]

    rows = []
    all_positive = True
    for lam in lambdas:
        values = []
        for w in ws:
            value = psi(lam, w)
            all_positive = all_positive and value > 0
            values.append({"w": str(w), "psi": mp.nstr(value, args.dps)})
        rows.append({"lambda": str(lam), "values": values})

    # Independent derivative consistency checks away from the endpoints.
    derivative_checks = []
    for lam, w in [
        (mp.mpf("2"), mp.mpf("0.4")),
        (mp.mpf("4.7243834"), mp.mpf("0.6")),
        (mp.mpf("10"), mp.mpf("0.8")),
    ]:
        direct = psi(lam, w)
        differentiated = mp.diff(lambda z: energy(lam, z), w)
        derivative_checks.append(
            {
                "lambda": str(lam),
                "w": str(w),
                "psi_direct": mp.nstr(direct, args.dps),
                "energy_derivative": mp.nstr(differentiated, args.dps),
                "absolute_difference": mp.nstr(
                    abs(direct - differentiated), args.dps
                ),
            }
        )

    result = {
        "status": "NON_CERTIFIED_REFERENCE",
        "all_sampled_psi_positive": all_positive,
        "decimal_precision": args.dps,
        "lambda_samples": rows,
        "derivative_consistency": derivative_checks,
        "target": "psi(lambda,w)>0 for all lambda>=1 and 0<w<1",
        "warning": "mpmath quadrature is reference evidence, not interval certification",
        "script_sha256": sha256_file(Path(__file__)),
    }

    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
