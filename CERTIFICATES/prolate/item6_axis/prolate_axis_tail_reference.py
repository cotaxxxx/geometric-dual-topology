#!/usr/bin/env python3
"""Non-certified numerical audit of the item 6 logarithmic tail coefficient."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from prolate_axis_reference import psi, sha256_file


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=40)
    parser.add_argument("--lambda-low", default="1e5")
    parser.add_argument("--lambda-high", default="1e6")
    parser.add_argument("--relative-tolerance", default="5e-5")
    parser.add_argument("--json", default="prolate_axis_tail_reference.json")
    args = parser.parse_args()

    mp.mp.dps = args.dps
    lam_lo = mp.mpf(args.lambda_low)
    lam_hi = mp.mpf(args.lambda_high)
    tolerance = mp.mpf(args.relative_tolerance)
    if not (1 < lam_lo < lam_hi):
        raise ValueError("require 1 < lambda-low < lambda-high")

    ws = [
        mp.mpf("0.1"),
        mp.mpf("0.2"),
        mp.mpf("0.35"),
        mp.mpf("0.5"),
        mp.mpf("0.65"),
        mp.mpf("0.8"),
        mp.mpf("0.9"),
    ]

    rows = []
    all_positive = True
    all_within_tolerance = True
    log_ratio = mp.log(lam_hi / lam_lo)

    for w in ws:
        psi_lo = psi(lam_lo, w)
        psi_hi = psi(lam_hi, w)
        all_positive = all_positive and psi_lo > 0 and psi_hi > 0

        # H_lambda(w)=lambda*Psi_lambda(w)/w should be affine in log(lambda)
        # to leading order. Its decade slope isolates the logarithmic
        # coefficient without needing the finite remainder.
        h_lo = lam_lo * psi_lo / w
        h_hi = lam_hi * psi_hi / w
        slope = (h_hi - h_lo) / log_ratio
        exact = 3 * mp.pi * mp.sqrt(1 - w * w)
        relative_error = abs(slope - exact) / exact
        within = relative_error < tolerance
        all_within_tolerance = all_within_tolerance and within

        rows.append(
            {
                "w": str(w),
                "scaled_low": mp.nstr(h_lo, args.dps),
                "scaled_high": mp.nstr(h_hi, args.dps),
                "decade_slope": mp.nstr(slope, args.dps),
                "exact_coefficient": mp.nstr(exact, args.dps),
                "relative_error": mp.nstr(relative_error, 12),
                "within_tolerance": bool(within),
            }
        )

    result = {
        "status": "PASSED_REFERENCE_CHECK"
        if all_positive and all_within_tolerance
        else "FAILED_REFERENCE_CHECK",
        "certification_status": "NON_CERTIFIED_REFERENCE",
        "decimal_precision": args.dps,
        "lambda_low": str(lam_lo),
        "lambda_high": str(lam_hi),
        "relative_tolerance": str(tolerance),
        "all_sampled_psi_positive": bool(all_positive),
        "all_seven_slopes_within_tolerance": bool(all_within_tolerance),
        "rows": rows,
        "asymptotic_target": (
            "Psi_lambda(w) = (w/lambda)*(3*pi*sqrt(1-w^2)*log(lambda) "
            "+ B(1/lambda,w^2))"
        ),
        "warning": (
            "This is a floating-point regression audit. It neither bounds the "
            "remainder uniformly nor certifies positivity on a box."
        ),
        "script_sha256": sha256_file(Path(__file__)),
    }

    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "PASSED_REFERENCE_CHECK" else 1)


if __name__ == "__main__":
    main()
