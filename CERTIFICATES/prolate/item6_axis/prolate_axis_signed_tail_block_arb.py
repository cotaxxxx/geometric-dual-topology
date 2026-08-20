#!/usr/bin/env python3
"""Mu-scaled signed-angle Arb certificate for compact tail H or M blocks.

Uses the exact signed angle

    delta = atan(sqrt(1-c^2)*(w-c+mu^2*c)/(mu*(1-w*c)))

and the cancellation-regularized kernels audited in
``prolate_axis_signed_tail_symbolic_audit.py``.  The moving layer is resolved by
three correlated charts

    [-1,w-4mu], [w-4mu,w+4mu], [w+4mu,1].
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from flint import acb, arb, ctx, fmpq

import prolate_axis_tail_H_block_arb as cover


LAYER_RADIUS = 4


def signed_tail_three_chart(
    mu_value: arb,
    w_value: arb,
    quantity: str,
    tolerance: arb,
    integration_depth: int,
    eval_limit: int,
) -> acb:
    mu = acb(mu_value)
    w = acb(w_value)
    mu2 = mu * mu
    radius = acb(LAYER_RADIUS)
    pi = acb(arb.pi())
    angle_floor = pi * pi / 4

    def base_integrand(c: acb, analytic: bool) -> acb:
        c2 = c * c
        rho2 = 1 - c2
        rho = rho2.sqrt(analytic=analytic)
        n = 1 - w * c
        difference = c - w
        p = difference * difference + mu2 * rho2
        s = rho2 + mu2 * c2
        y = w - c + mu2 * c
        z = c - w + mu2 * c
        delta = (rho * y / (mu * n)).atan()

        regularized_square = delta * delta - angle_floor
        if quantity == "H":
            return (
                -c * regularized_square / (2 * mu * w)
                + n * delta * rho / (w * p)
            )

        delta_mu = rho * n * z / (p * s)
        return (
            c / (2 * w) * (
                2 * delta * delta_mu - regularized_square / mu
            )
            - mu * n * rho * delta_mu / (w * p)
            + 2 * mu2 * n * rho * rho2 * delta / (w * p * p)
        )

    left_edge = w - radius * mu
    right_edge = w + radius * mu

    def left_kernel(t: acb, analytic: bool) -> acb:
        jacobian = 1 + left_edge
        c = -1 + jacobian * t
        return jacobian * base_integrand(c, analytic)

    def inner_kernel(t: acb, analytic: bool) -> acb:
        jacobian = 2 * radius * mu
        c = left_edge + jacobian * t
        return jacobian * base_integrand(c, analytic)

    def right_kernel(t: acb, analytic: bool) -> acb:
        jacobian = 1 - right_edge
        c = right_edge + jacobian * t
        return jacobian * base_integrand(c, analytic)

    return (
        cover.rigorous_integral(
            left_kernel, tolerance, integration_depth, eval_limit
        )
        + cover.rigorous_integral(
            inner_kernel, tolerance, integration_depth, eval_limit
        )
        + cover.rigorous_integral(
            right_kernel, tolerance, integration_depth, eval_limit
        )
    )


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
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=45)
    parser.add_argument("--tolerance", default="1e-16")
    parser.add_argument("--integration-depth", type=int, default=20)
    parser.add_argument("--eval-limit", type=int, default=100000)
    parser.add_argument("--max-split-depth", type=int, default=14)
    parser.add_argument("--max-boxes", type=int, default=8192)
    parser.add_argument("--mu-lo-num", type=int, required=True)
    parser.add_argument("--mu-lo-den", type=int, required=True)
    parser.add_argument("--mu-hi-num", type=int, required=True)
    parser.add_argument("--mu-hi-den", type=int, required=True)
    parser.add_argument("--w-lo-num", type=int, required=True)
    parser.add_argument("--w-lo-den", type=int, required=True)
    parser.add_argument("--w-hi-num", type=int, required=True)
    parser.add_argument("--w-hi-den", type=int, required=True)
    parser.add_argument("--quantity", choices=("H", "M"), required=True)
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
    if not (-1 < w_lo - LAYER_RADIUS * mu_hi):
        raise ValueError("left outer chart leaves [-1,1]")
    if not (w_hi + LAYER_RADIUS * mu_hi < 1):
        raise ValueError("right outer chart leaves [-1,1]")

    ctx.dps = args.dps
    cover.scaled_tail_H = lambda mu, w, tol, depth, limit: signed_tail_three_chart(
        mu, w, args.quantity, tol, depth, limit
    )
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
    if args.quantity == "M":
        result = rename_H_to_M(result)

    result["quantity"] = args.quantity
    result["representation"] = (
        "signed atan angle; exact pi^2/4 cancellation; mu-scaled three-chart split"
    )
    result["integration_chart"] = {
        "type": "mu-scaled three-chart moving-layer split",
        "layer_radius": LAYER_RADIUS,
        "pieces": ["[-1,w-4mu]", "[w-4mu,w+4mu]", "[w+4mu,1]"],
    }
    result["script_sha256"] = cover.sha256_file(Path(__file__))
    result["formula_audit_sha256"] = cover.sha256_file(
        Path(__file__).with_name("prolate_axis_signed_tail_symbolic_audit.py")
    )
    result["base_cover_sha256"] = cover.sha256_file(
        Path(__file__).with_name("prolate_axis_tail_H_block_arb.py")
    )

    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{result['script_sha256']}  {Path(__file__).name}\n"
        f"{result['formula_audit_sha256']}  prolate_axis_signed_tail_symbolic_audit.py\n"
        f"{result['base_cover_sha256']}  prolate_axis_tail_H_block_arb.py\n"
        f"{cover.sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
