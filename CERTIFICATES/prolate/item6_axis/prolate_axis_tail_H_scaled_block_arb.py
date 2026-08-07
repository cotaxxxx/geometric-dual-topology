#!/usr/bin/env python3
"""Direct-Psi compact-tail certificate on paired signed-angle charts.

On positive compact tail blocks, ``H=Psi/(mu*w)`` has the same sign as Psi.
The driver therefore certifies Psi directly, before division by ``mu*w``.  This
removes all artificial ``1/mu`` amplification from the inner moving layer.
Endpoint-square caps make ``rho=sqrt(1-c^2)`` analytic at ``c=+/-1``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from flint import acb, arb, ctx, fmpq

import prolate_axis_tail_H_block_arb as base


LAYER_RADIUS = 4
CAP_SIZE_Q = fmpq(1, 256)
CAP_START_Q = 1 - CAP_SIZE_Q
EVALUATION_TRACE: list[dict] = []
TRACE_LIMIT = 32


def mapped_integral(kernel, lo, hi, tolerance, depth, limit) -> acb:
    lo_c = acb(arb(str(lo)))
    jacobian = acb(arb(str(hi - lo)))

    def mapped(t: acb, analytic: bool) -> acb:
        return jacobian * kernel(lo_c + jacobian * t, analytic)

    return base.rigorous_integral(mapped, tolerance, depth, limit)


def signed_psi_density(
    mu: acb,
    w: acb,
    c: acb,
    difference: acb,
    rho2: acb,
    rho: acb,
) -> acb:
    mu2 = mu * mu
    n = 1 - w * c
    p = difference * difference + mu2 * rho2
    y = -difference + mu2 * c
    delta = (rho * y / (mu * n)).atan()
    pi = acb(arb.pi())
    regularized_square = delta * delta - pi * pi / 4
    return -c * regularized_square / 2 + mu * n * delta * rho / p


def paired_tail_Psi(mu_value, w_value, tolerance, depth, limit) -> acb:
    mu = acb(mu_value)
    w = acb(w_value)
    radius = acb(LAYER_RADIUS)
    cap = acb(arb(str(CAP_SIZE_Q)))
    mu2 = mu * mu

    def inner_kernel(t: acb, analytic: bool) -> acb:
        q = 2 * radius * t - radius
        c = w + mu * q
        rho2 = (1 - c) * (1 + c)
        rho = rho2.sqrt(analytic=analytic)
        reduced_p = q * q + rho2
        n = 1 - w * c
        y = -mu * q + mu2 * c
        delta = (rho * y / (mu * n)).atan()
        pi = acb(arb.pi())
        regularized_square = delta * delta - pi * pi / 4
        # dc=8*mu*dt and P=mu^2*reduced_p.
        return (
            -radius * mu * c * regularized_square
            + 2 * radius * n * delta * rho / reduced_p
        )

    pair_jacobian = 1 - w - radius * mu

    def paired_bulk_kernel(t: acb, analytic: bool) -> acb:
        x = radius * mu + pair_jacobian * t
        c_left = w - x
        rho2_left = (1 - w + x) * (1 + w - x)
        rho_left = rho2_left.sqrt(analytic=analytic)
        left = signed_psi_density(mu, w, c_left, -x, rho2_left, rho_left)
        c_right = w + x
        rho2_right = pair_jacobian * (1 - t) * (1 + w + x)
        rho_right = rho2_right.sqrt(analytic=analytic)
        right = signed_psi_density(mu, w, c_right, x, rho2_right, rho_right)
        return pair_jacobian * (left + right)

    def paired_cap_kernel(s: acb, analytic: bool) -> acb:
        u = 1 - s
        t = 1 - cap * u * u
        dt_ds = 2 * cap * u
        x = radius * mu + pair_jacobian * t
        c_left = w - x
        rho2_left = (1 - w + x) * (1 + w - x)
        rho_left = rho2_left.sqrt(analytic=analytic)
        left = signed_psi_density(mu, w, c_left, -x, rho2_left, rho_left)
        c_right = w + x
        right_factor = pair_jacobian * cap * (1 + w + x)
        rho_right = u * right_factor.sqrt(analytic=analytic)
        rho2_right = rho_right * rho_right
        right = signed_psi_density(mu, w, c_right, x, rho2_right, rho_right)
        return dt_ds * pair_jacobian * (left + right)

    def far_bulk_kernel(t: acb, analytic: bool) -> acb:
        x = 1 - w + 2 * w * t
        c = w - x
        rho2 = 4 * w * (1 - t) * (1 - w + w * t)
        rho = rho2.sqrt(analytic=analytic)
        return 2 * w * signed_psi_density(mu, w, c, -x, rho2, rho)

    def far_cap_kernel(s: acb, analytic: bool) -> acb:
        u = 1 - s
        t = 1 - cap * u * u
        dt_ds = 2 * cap * u
        x = 1 - w + 2 * w * t
        c = w - x
        rho_factor = 4 * w * cap * (1 - w + w * t)
        rho = u * rho_factor.sqrt(analytic=analytic)
        rho2 = rho * rho
        return dt_ds * 2 * w * signed_psi_density(mu, w, c, -x, rho2, rho)

    inner = base.rigorous_integral(inner_kernel, tolerance, depth, limit)
    paired_bulk = mapped_integral(
        paired_bulk_kernel, fmpq(0), CAP_START_Q, tolerance, depth, limit
    )
    paired_cap = base.rigorous_integral(paired_cap_kernel, tolerance, depth, limit)
    far_bulk = mapped_integral(
        far_bulk_kernel, fmpq(0), CAP_START_Q, tolerance, depth, limit
    )
    far_cap = base.rigorous_integral(far_cap_kernel, tolerance, depth, limit)
    paired = paired_bulk + paired_cap
    far_left = far_bulk + far_cap
    total = inner + paired + far_left
    if len(EVALUATION_TRACE) < TRACE_LIMIT:
        EVALUATION_TRACE.append({
            "mu_ball": str(mu_value),
            "mu_lower": str(mu_value.lower()),
            "mu_upper": str(mu_value.upper()),
            "w_ball": str(w_value),
            "w_lower": str(w_value.lower()),
            "w_upper": str(w_value.upper()),
            "inner": base.acb_record(inner),
            "paired_bulk": base.acb_record(paired_bulk),
            "paired_cap": base.acb_record(paired_cap),
            "paired_outer": base.acb_record(paired),
            "far_bulk": base.acb_record(far_bulk),
            "far_cap": base.acb_record(far_cap),
            "far_left": base.acb_record(far_left),
            "total": base.acb_record(total),
        })
    return total


def rename_H_to_Psi(result: dict) -> dict:
    result["scope"] = result["scope"].replace("H(mu,w)>0", "Psi_{1/mu}(w)>0")
    result["derived_statement"] = (
        "H(mu,w)>0 on the same block because H=Psi_{1/mu}(w)/(mu*w) and mu*w>0."
    )
    result["conditions"] = {
        key.replace("H > 0", "Psi > 0"): value
        for key, value in result["conditions"].items()
    }
    for leaf in result.get("leaves", []):
        if "H" in leaf:
            leaf["Psi"] = leaf.pop("H")
    for box in result.get("terminal", []):
        if "H" in box:
            box["Psi"] = box.pop("H")
    worst = result.get("worst_certified_leaf")
    if worst is not None and "H" in worst:
        worst["Psi"] = worst.pop("H")
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
    EVALUATION_TRACE.clear()
    base.scaled_tail_H = paired_tail_Psi
    result = base.certify_block(
        args.dps, args.tolerance, args.integration_depth, args.eval_limit,
        args.max_split_depth, args.max_boxes,
        mu_lo, mu_hi, w_lo, w_hi,
    )
    result = rename_H_to_Psi(result)
    result["evaluation_trace"] = EVALUATION_TRACE
    result["quantity"] = "Psi"
    result["representation"] = (
        "direct signed-angle Psi density before division by mu*w; paired charts; "
        "endpoint-square caps"
    )
    result["integration_chart"] = {
        "type": "direct-Psi signed-angle paired split with endpoint-square caps",
        "layer_radius": LAYER_RADIUS,
        "cap_size": str(CAP_SIZE_Q),
        "pieces": [
            "inner: c=w+mu*(8t-4)",
            "paired bulk: t in [0,255/256]",
            "paired endpoint cap: t=1-(1-s)^2/256",
            "far-left bulk: t in [0,255/256]",
            "far-left endpoint cap: t=1-(1-s)^2/256",
        ],
    }
    result["script_sha256"] = base.sha256_file(Path(__file__))
    result["base_driver_sha256"] = base.sha256_file(
        Path(__file__).with_name("prolate_axis_tail_H_block_arb.py")
    )
    result["formula_audit_sha256"] = base.sha256_file(
        Path(__file__).with_name("prolate_axis_signed_tail_symbolic_audit.py")
    )

    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{result['script_sha256']}  {Path(__file__).name}\n"
        f"{result['base_driver_sha256']}  prolate_axis_tail_H_block_arb.py\n"
        f"{result['formula_audit_sha256']}  prolate_axis_signed_tail_symbolic_audit.py\n"
        f"{base.sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
