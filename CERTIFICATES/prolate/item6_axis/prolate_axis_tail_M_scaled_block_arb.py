#!/usr/bin/env python3
"""Paired, correlation-preserving Arb driver for M=-mu*partial_mu H.

The inner layer is blown up by ``c=w+mu*q``.  The two outer sides are evaluated
at a common distance ``x=|c-w|`` inside one integration kernel, followed by the
unmatched far-left interval.  This exposes the dominant left/right cancellation
before parameter subdivision and preserves all moving-layer correlations.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from flint import acb, arb, ctx, fmpq

import prolate_axis_tail_H_block_arb as cover
from prolate_axis_center_cap_arb import regular_angle_derivatives


LAYER_RADIUS = 4


def paired_tail_M(
    mu_value: arb,
    w_value: arb,
    tolerance: arb,
    integration_depth: int,
    eval_limit: int,
) -> acb:
    mu = acb(mu_value)
    w = acb(w_value)
    mu2 = mu * mu
    radius = acb(LAYER_RADIUS)

    def outer_integrand(
        c: acb,
        difference: acb,
        rho2: acb,
        n: acb,
        analytic: bool,
    ) -> acb:
        p = difference * difference + mu2 * rho2
        s = rho2 + mu2 * c * c
        root = (p * s).sqrt(analytic=analytic)
        cosine = mu * n / root
        h1, h2 = regular_angle_derivatives(cosine)
        hbar = cover.regular_hbar(cosine)

        b = -c + n * difference / p
        q_metric = mu2 * (rho2 / p + c * c / s)
        alpha = 1 - q_metric
        beta = -2 * mu2 * n * difference * rho2 / (p * p)
        f = -c * hbar + h1 * b
        bracket = (
            q_metric * f
            + c * alpha * (h1 - hbar)
            - alpha * h2 * cosine * b
            - h1 * beta
        )
        return n * bracket / (2 * w * root)

    def inner_kernel(t: acb, analytic: bool) -> acb:
        q_coordinate = 2 * radius * t - radius
        c = w + mu * q_coordinate
        rho2 = (1 - c) * (1 + c)
        reduced_p = q_coordinate * q_coordinate + rho2
        s = rho2 + mu2 * c * c
        reduced_root = (reduced_p * s).sqrt(analytic=analytic)
        n = 1 - w * w - mu * w * q_coordinate
        cosine = n / reduced_root
        h1, h2 = regular_angle_derivatives(cosine)
        hbar = cover.regular_hbar(cosine)

        b = -c + n * q_coordinate / (mu * reduced_p)
        q_metric = rho2 / reduced_p + mu2 * c * c / s
        alpha = 1 - q_metric
        beta = -2 * n * q_coordinate * rho2 / (
            mu * reduced_p * reduced_p
        )
        f = -c * hbar + h1 * b
        bracket = (
            q_metric * f
            + c * alpha * (h1 - hbar)
            - alpha * h2 * cosine * b
            - h1 * beta
        )
        return radius * n * bracket / (w * reduced_root)

    pair_jacobian = 1 - w - radius * mu

    def paired_outer_kernel(t: acb, analytic: bool) -> acb:
        x = radius * mu + pair_jacobian * t

        c_left = w - x
        rho2_left = (1 - w + x) * (1 + w - x)
        n_left = 1 - w * w + w * x
        left = outer_integrand(
            c_left, -x, rho2_left, n_left, analytic
        )

        c_right = w + x
        rho2_right = pair_jacobian * (1 - t) * (1 + w + x)
        n_right = 1 - w * w - w * x
        right = outer_integrand(
            c_right, x, rho2_right, n_right, analytic
        )

        return pair_jacobian * (left + right)

    def far_left_kernel(t: acb, analytic: bool) -> acb:
        x = 1 - w + 2 * w * t
        c = w - x
        rho2 = 4 * w * (1 - t) * (1 - w + w * t)
        n = 1 + w - 2 * w * w + 2 * w * w * t
        return 2 * w * outer_integrand(c, -x, rho2, n, analytic)

    inner = cover.rigorous_integral(
        inner_kernel, tolerance, integration_depth, eval_limit
    )
    paired = cover.rigorous_integral(
        paired_outer_kernel, tolerance, integration_depth, eval_limit
    )
    far_left = cover.rigorous_integral(
        far_left_kernel, tolerance, integration_depth, eval_limit
    )
    return inner + paired + far_left


def rename_H_to_M(result: dict) -> dict:
    result["scope"] = result["scope"].replace("H(mu,w)", "M(mu,w)")
    result["derived_statement"] = (
        "M=-mu*partial_mu H>0 on this block, so H increases when mu decreases."
    )
    result["definition"] = "M(mu,w)=-mu*partial_mu H(mu,w)"
    result["conditions"] = {
        key.replace("H > 0", "M > 0"): value
        for key, value in result["conditions"].items()
    }
    for leaf in result.get("leaves", []):
        leaf["M"] = leaf.pop("H")
    for box in result.get("terminal", []):
        if "H" in box:
            box["M"] = box.pop("H")
    worst = result.get("worst_certified_leaf")
    if worst is not None and "H" in worst:
        worst["M"] = worst.pop("H")
    result["integration_chart"] = {
        "type": "paired correlation-preserving moving-layer split",
        "layer_radius": LAYER_RADIUS,
        "pieces": [
            "inner: c=w+mu*(8t-4)",
            "paired: c=w+-x, x=4mu+(1-w-4mu)t",
            "far-left: x=1-w+2wt, c=w-x",
        ],
    }
    result["limitations"] = (
        "This is a compact positive-mu block. A separate endpoint argument is "
        "required to include mu=0."
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--tolerance", default="1e-18")
    parser.add_argument("--integration-depth", type=int, default=22)
    parser.add_argument("--eval-limit", type=int, default=200000)
    parser.add_argument("--max-split-depth", type=int, default=16)
    parser.add_argument("--max-boxes", type=int, default=32768)
    parser.add_argument("--mu-lo-num", type=int, required=True)
    parser.add_argument("--mu-lo-den", type=int, required=True)
    parser.add_argument("--mu-hi-num", type=int, required=True)
    parser.add_argument("--mu-hi-den", type=int, required=True)
    parser.add_argument("--w-lo-num", type=int, required=True)
    parser.add_argument("--w-lo-den", type=int, required=True)
    parser.add_argument("--w-hi-num", type=int, required=True)
    parser.add_argument("--w-hi-den", type=int, required=True)
    parser.add_argument("--json", required=True)
    args = parser.parse_args()

    mu_lo = fmpq(args.mu_lo_num, args.mu_lo_den)
    mu_hi = fmpq(args.mu_hi_num, args.mu_hi_den)
    w_lo = fmpq(args.w_lo_num, args.w_lo_den)
    w_hi = fmpq(args.w_hi_num, args.w_hi_den)
    if not (fmpq(0) < mu_lo < mu_hi):
        raise ValueError("require 0 < mu-lo < mu-hi")
    if not (fmpq(0) < w_lo < w_hi < fmpq(1)):
        raise ValueError("require 0 < w-lo < w-hi < 1")
    if not (w_hi + LAYER_RADIUS * mu_hi < 1):
        raise ValueError("paired outer chart has nonpositive length")

    ctx.dps = args.dps
    cover.scaled_tail_H = paired_tail_M
    result = cover.certify_block(
        args.dps,
        args.tolerance,
        args.integration_depth,
        args.eval_limit,
        args.max_split_depth,
        args.max_boxes,
        mu_lo,
        mu_hi,
        w_lo,
        w_hi,
    )
    result = rename_H_to_M(result)
    result["script_sha256"] = cover.sha256_file(Path(__file__))
    result["formula_audit_sha256"] = cover.sha256_file(
        Path(__file__).with_name("prolate_axis_tail_log_derivative_symbolic_audit.py")
    )
    result["base_cover_sha256"] = cover.sha256_file(
        Path(__file__).with_name("prolate_axis_tail_H_block_arb.py")
    )

    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{result['script_sha256']}  {Path(__file__).name}\n"
        f"{result['formula_audit_sha256']}  prolate_axis_tail_log_derivative_symbolic_audit.py\n"
        f"{result['base_cover_sha256']}  prolate_axis_tail_H_block_arb.py\n"
        f"{cover.sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
