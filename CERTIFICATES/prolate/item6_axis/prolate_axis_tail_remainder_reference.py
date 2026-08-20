#!/usr/bin/env python3
"""Non-certified scout for a uniform lower target on the tail remainder."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from prolate_axis_reference import psi, sha256_file


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=35)
    parser.add_argument("--json", default="prolate_axis_tail_remainder_reference.json")
    args = parser.parse_args()
    mp.mp.dps = args.dps

    mus = [
        mp.mpf("1e-2"),
        mp.mpf("5e-3"),
        mp.mpf("2e-3"),
        mp.mpf("1e-3"),
        mp.mpf("1e-4"),
    ]
    ws = [
        mp.mpf("0.01"),
        mp.mpf("0.05"),
        mp.mpf("0.1"),
        mp.mpf("0.25"),
        mp.mpf("0.5"),
        mp.mpf("0.75"),
    ]
    target = mp.mpf("-7")
    rows = []
    sampled_min = None
    sampled_argmin = None

    for mu in mus:
        lam = 1 / mu
        for w in ws:
            scaled = psi(lam, w) / (mu * w)
            leading = 3 * mp.pi * mp.sqrt(1 - w * w) * mp.log(1 / mu)
            remainder = scaled - leading
            if sampled_min is None or remainder < sampled_min:
                sampled_min = remainder
                sampled_argmin = (mu, w)
            rows.append({
                "mu": str(mu),
                "lambda": str(lam),
                "w": str(w),
                "scaled_H": mp.nstr(scaled, args.dps),
                "leading_term": mp.nstr(leading, args.dps),
                "Bhat": mp.nstr(remainder, args.dps),
                "above_minus_7": bool(remainder > target),
            })

    mu0 = mp.mpf("1e-2")
    w1 = mp.mpf("0.75")
    minimum_leading_on_target_box = (
        3 * mp.pi * mp.sqrt(1 - w1 * w1) * mp.log(1 / mu0)
    )

    result = {
        "status": "PASSED_REFERENCE_CHECK"
        if sampled_min is not None and sampled_min > target
        else "FAILED_REFERENCE_CHECK",
        "certification_status": "NON_CERTIFIED_REFERENCE",
        "target_box": "0<mu<=1/100, 0<w<=3/4",
        "proposed_uniform_remainder_bound": "Bhat(mu,w^2)>-7",
        "sampled_minimum": mp.nstr(sampled_min, args.dps),
        "sampled_argmin": {
            "mu": str(sampled_argmin[0]),
            "w": str(sampled_argmin[1]),
        },
        "minimum_leading_term_on_closed_outer_corner": mp.nstr(
            minimum_leading_on_target_box, args.dps
        ),
        "implication": (
            "A certified bound Bhat>-7 on the target box would imply H>0, "
            "because the logarithmic leading term is minimized at "
            "mu=1/100, w=3/4 and is greater than 7."
        ),
        "rows": rows,
        "warning": (
            "This grid is floating-point evidence only. It does not prove a "
            "uniform lower bound or control the limits mu->0 and w->0."
        ),
        "script_sha256": sha256_file(Path(__file__)),
    }

    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "PASSED_REFERENCE_CHECK" else 1)


if __name__ == "__main__":
    main()
