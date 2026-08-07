#!/usr/bin/env python3
"""Non-certified high-precision scout for the item 6 axial center Hessian."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from prolate_axis_reference import sha256_file


def angle_derivatives(cosine: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    if cosine > 1 and cosine - 1 < 100 * mp.eps:
        cosine = mp.mpf(1)
    u = mp.acos(cosine)
    if abs(u) < mp.sqrt(mp.eps):
        h1 = -2 * (1 + u * u / 6 + 7 * u**4 / 360)
        h2 = mp.mpf(2) / 3 + mp.mpf(4) * u * u / 15
    else:
        sine = mp.sin(u)
        h1 = -2 * u / sine
        h2 = 2 * (sine - u * mp.cos(u)) / sine**3
    return h1, h2


def center_kernel(lam: mp.mpf, c: mp.mpf) -> mp.mpf:
    lam2 = lam * lam
    ell = 1 + (lam2 - 1) * c * c
    w2 = lam2 * (1 - c * c) + c * c
    root = mp.sqrt(w2)
    cosine = lam / (mp.sqrt(ell) * root)
    h1, h2 = angle_derivatives(cosine)

    cosine_w = (
        lam
        * c
        * (lam2 - 1)
        * (1 - c * c)
        / (ell ** mp.mpf("1.5") * root)
    )
    cosine_ww = (
        lam**3
        * (1 - c * c)
        * (2 * (lam2 - 1) * c * c - 1)
        / (ell ** mp.mpf("2.5") * root)
    )
    return mp.mpf("0.5") * (
        -2 * c * h1 * cosine_w + h2 * cosine_w**2 + h1 * cosine_ww
    )


def q_parallel(lam: mp.mpf) -> mp.mpf:
    return 2 * mp.quad(lambda c: center_kernel(lam, c), [0, 1])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=40)
    parser.add_argument("--json", default="prolate_axis_center_reference.json")
    args = parser.parse_args()
    mp.mp.dps = args.dps

    lambdas = [
        mp.mpf("1"),
        mp.mpf("1.1"),
        mp.mpf("1.5"),
        mp.mpf("2"),
        mp.mpf("2.0653823"),
        mp.mpf("3.4348684"),
        mp.mpf("4.7243834"),
        mp.mpf("10"),
        mp.mpf("100"),
        mp.mpf("1e3"),
        mp.mpf("1e4"),
        mp.mpf("1e6"),
    ]

    rows = []
    all_positive = True
    for lam in lambdas:
        value = q_parallel(lam)
        all_positive = all_positive and value > 0
        rows.append({"lambda": str(lam), "Q_parallel": mp.nstr(value, args.dps)})

    lam_lo = mp.mpf("1e5")
    lam_hi = mp.mpf("1e6")
    scaled_lo = lam_lo * q_parallel(lam_lo)
    scaled_hi = lam_hi * q_parallel(lam_hi)
    tail_slope = (scaled_hi - scaled_lo) / mp.log(lam_hi / lam_lo)
    tail_exact = 3 * mp.pi
    tail_relative_error = abs(tail_slope - tail_exact) / tail_exact

    sphere_value = q_parallel(mp.mpf(1))
    sphere_error = abs(sphere_value - mp.mpf(4) / 3)
    result = {
        "status": "NON_CERTIFIED_REFERENCE",
        "all_sampled_Q_parallel_positive": bool(all_positive),
        "sphere_Q_parallel": mp.nstr(sphere_value, args.dps),
        "sphere_absolute_error_from_4_over_3": mp.nstr(sphere_error, 12),
        "lambda_samples": rows,
        "tail_center_decade_slope": mp.nstr(tail_slope, args.dps),
        "tail_center_exact_coefficient": mp.nstr(tail_exact, args.dps),
        "tail_center_relative_error": mp.nstr(tail_relative_error, 12),
        "warning": (
            "The center Hessian values use mpmath quadrature and do not certify "
            "positivity on a lambda interval."
        ),
        "script_sha256": sha256_file(Path(__file__)),
    }

    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
