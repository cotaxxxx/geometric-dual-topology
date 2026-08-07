#!/usr/bin/env python3
"""Adaptive Arb certificate using the exact signed-angle reduction.

Certifies either Psi_lambda(w) or A_lambda''(w) with a requested strict sign on
one exact rational (w,lambda) rectangle.  The c integral is split at c=w and
mapped to two fixed unit intervals.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from collections import deque
from pathlib import Path
from typing import Callable

import flint
from flint import acb, arb, ctx, fmpq


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def closed_interval(lo: fmpq, hi: fmpq) -> arb:
    return arb(str((lo + hi) / 2), str((hi - lo) / 2))


def acb_record(value: acb) -> dict:
    return {
        "real_ball": str(value.real),
        "real_lower": str(value.real.lower()),
        "real_upper": str(value.real.upper()),
        "imag_ball": str(value.imag),
        "imag_contains_zero": bool(0 in value.imag),
    }


def rigorous_integral(
    kernel: Callable[[acb, bool], acb],
    tolerance: arb,
    depth_limit: int,
    eval_limit: int,
) -> acb:
    return acb.integral(
        kernel,
        0,
        1,
        abs_tol=tolerance,
        rel_tol=tolerance,
        depth_limit=depth_limit,
        eval_limit=eval_limit,
    )


def signed_angle_quantities(
    c: acb,
    w: acb,
    lam: acb,
    analytic: bool,
) -> tuple[acb, acb, acb]:
    c2 = c * c
    rho = (1 - c2).sqrt(analytic=analytic)
    n = 1 - w * c
    lam2 = lam * lam
    r2 = 1 - c2 + lam2 * (c - w) ** 2
    cross = rho * (lam * w - (lam - 1 / lam) * c)
    delta = (cross / n).atan()
    delta_w = lam * rho / r2
    delta_ww = 2 * lam * lam2 * rho * (c - w) / (r2 * r2)
    return delta, delta_w, delta_ww


def signed_integral(
    w_value: arb,
    lambda_value: arb,
    quantity: str,
    tolerance: arb,
    integration_depth: int,
    eval_limit: int,
) -> acb:
    w = acb(w_value)
    lam = acb(lambda_value)

    def base(c: acb, analytic: bool) -> acb:
        delta, delta_w, delta_ww = signed_angle_quantities(c, w, lam, analytic)
        n = 1 - w * c
        if quantity == "Psi":
            return -c * delta * delta / 2 + n * delta * delta_w
        return -2 * c * delta * delta_w + n * (
            delta_w * delta_w + delta * delta_ww
        )

    def left_kernel(t: acb, analytic: bool) -> acb:
        jacobian = 1 + w
        c = -1 + jacobian * t
        return jacobian * base(c, analytic)

    def right_kernel(t: acb, analytic: bool) -> acb:
        jacobian = 1 - w
        c = w + jacobian * t
        return jacobian * base(c, analytic)

    return rigorous_integral(
        left_kernel, tolerance, integration_depth, eval_limit
    ) + rigorous_integral(
        right_kernel, tolerance, integration_depth, eval_limit
    )


def certify_rectangle(
    dps: int,
    tolerance_text: str,
    integration_depth: int,
    eval_limit: int,
    max_split_depth: int,
    max_boxes: int,
    w_lo: fmpq,
    w_hi: fmpq,
    lambda_lo: fmpq,
    lambda_hi: fmpq,
    quantity: str,
    sign: str,
) -> dict:
    ctx.dps = dps
    tolerance = arb(tolerance_text)
    w_span = w_hi - w_lo
    lambda_span = lambda_hi - lambda_lo
    queue: deque[tuple[fmpq, fmpq, fmpq, fmpq, int]] = deque(
        [(w_lo, w_hi, lambda_lo, lambda_hi, 0)]
    )
    leaves: list[dict] = []
    terminal: list[dict] = []
    evaluations = 0
    extremal_bound = None
    extremal_leaf = None

    while queue:
        if evaluations >= max_boxes:
            while queue:
                wl, wh, ll, lh, depth = queue.popleft()
                terminal.append({
                    "w_lo": str(wl), "w_hi": str(wh),
                    "lambda_lo": str(ll), "lambda_hi": str(lh),
                    "split_depth": depth,
                    "reason": "max_boxes_exhausted",
                })
            break

        wl, wh, ll, lh, depth = queue.popleft()
        try:
            value = signed_integral(
                closed_interval(wl, wh),
                closed_interval(ll, lh),
                quantity,
                tolerance,
                integration_depth,
                eval_limit,
            )
            error = None
        except Exception as exc:
            value = None
            error = f"{type(exc).__name__}: {exc}"

        evaluations += 1
        accepted = bool(
            value is not None
            and 0 in value.imag
            and ((sign == "positive" and value.real > 0)
                 or (sign == "negative" and value.real < 0))
        )
        if accepted and value is not None:
            record = {
                "w_lo": str(wl), "w_hi": str(wh),
                "lambda_lo": str(ll), "lambda_hi": str(lh),
                "split_depth": depth,
                quantity: acb_record(value),
            }
            leaves.append(record)
            bound = value.real.lower() if sign == "positive" else value.real.upper()
            if (
                extremal_bound is None
                or (sign == "positive" and bound < extremal_bound)
                or (sign == "negative" and bound > extremal_bound)
            ):
                extremal_bound = bound
                extremal_leaf = record
            continue

        if depth < max_split_depth:
            w_relative = (wh - wl) / w_span
            lambda_relative = (lh - ll) / lambda_span
            if w_relative >= lambda_relative:
                midpoint = (wl + wh) / 2
                queue.append((wl, midpoint, ll, lh, depth + 1))
                queue.append((midpoint, wh, ll, lh, depth + 1))
            else:
                midpoint = (ll + lh) / 2
                queue.append((wl, wh, ll, midpoint, depth + 1))
                queue.append((wl, wh, midpoint, lh, depth + 1))
        else:
            record = {
                "w_lo": str(wl), "w_hi": str(wh),
                "lambda_lo": str(ll), "lambda_hi": str(lh),
                "split_depth": depth,
                "reason": "sign_not_certified" if error is None else "evaluation_error",
            }
            if value is not None:
                record[quantity] = acb_record(value)
            if error is not None:
                record["error"] = error
            terminal.append(record)

    leaf_area = sum(
        (
            (fmpq(a["w_hi"]) - fmpq(a["w_lo"]))
            * (fmpq(a["lambda_hi"]) - fmpq(a["lambda_lo"]))
            for a in leaves
        ),
        fmpq(0),
    )
    terminal_area = sum(
        (
            (fmpq(a["w_hi"]) - fmpq(a["w_lo"]))
            * (fmpq(a["lambda_hi"]) - fmpq(a["lambda_lo"]))
            for a in terminal
        ),
        fmpq(0),
    )
    target_area = w_span * lambda_span
    partition_ok = leaf_area + terminal_area == target_area
    exact_coverage = bool(leaves and not terminal and leaf_area == target_area)
    comparison = "> 0" if sign == "positive" else "< 0"
    conditions = {
        f"all accepted leaves have {quantity} {comparison}": bool(leaves),
        "exact binary-partition area invariant": partition_ok,
        "exact rational coverage of block": exact_coverage,
        "zero terminal boxes": not terminal,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"
    return {
        "status": status,
        "certified_sign": sign,
        "quantity": quantity,
        "scope": (
            f"{quantity}{comparison} for {w_lo}<=w<={w_hi}, "
            f"{lambda_lo}<=lambda<={lambda_hi}"
        ),
        "environment": {
            "python": platform.python_version(),
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
        "arithmetic": "python-flint Arb/Acb ball arithmetic",
        "representation": "exact signed atan angle with rational-algebraic derivatives",
        "decimal_precision": dps,
        "integration_tolerance": tolerance_text,
        "integration_depth_limit": integration_depth,
        "integration_eval_limit": eval_limit,
        "parameter_split_depth_limit": max_split_depth,
        "box_evaluation_limit": max_boxes,
        "target_rectangle": {
            "w": f"[{w_lo},{w_hi}]",
            "lambda": f"[{lambda_lo},{lambda_hi}]",
            "exact_area": str(target_area),
        },
        "partition": {
            "leaf_area": str(leaf_area),
            "terminal_area": str(terminal_area),
            "area_invariant": partition_ok,
        },
        "conditions": conditions,
        "counts": {
            "evaluations": evaluations,
            "certified_leaves": len(leaves),
            "terminal_boxes": len(terminal),
        },
        "extremal_certified_leaf": extremal_leaf,
        "leaves": leaves,
        "terminal": terminal,
        "limitations": "One exact finite rectangle; global claims require the strict grid combiner.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=45)
    parser.add_argument("--tolerance", default="1e-16")
    parser.add_argument("--integration-depth", type=int, default=20)
    parser.add_argument("--eval-limit", type=int, default=100000)
    parser.add_argument("--max-split-depth", type=int, default=14)
    parser.add_argument("--max-boxes", type=int, default=8192)
    parser.add_argument("--w-lo-num", type=int, required=True)
    parser.add_argument("--w-lo-den", type=int, required=True)
    parser.add_argument("--w-hi-num", type=int, required=True)
    parser.add_argument("--w-hi-den", type=int, required=True)
    parser.add_argument("--lambda-lo", type=int, required=True)
    parser.add_argument("--lambda-hi", type=int, required=True)
    parser.add_argument("--quantity", choices=("Psi", "A_second"), required=True)
    parser.add_argument("--sign", choices=("positive", "negative"), required=True)
    parser.add_argument("--json", required=True)
    args = parser.parse_args()

    w_lo = fmpq(args.w_lo_num, args.w_lo_den)
    w_hi = fmpq(args.w_hi_num, args.w_hi_den)
    lambda_lo = fmpq(args.lambda_lo)
    lambda_hi = fmpq(args.lambda_hi)
    if not (fmpq(0) <= w_lo < w_hi < fmpq(1)):
        raise ValueError("require 0 <= w-lo < w-hi < 1")
    if not (fmpq(1) <= lambda_lo < lambda_hi):
        raise ValueError("require 1 <= lambda-lo < lambda-hi")

    result = certify_rectangle(
        args.dps,
        args.tolerance,
        args.integration_depth,
        args.eval_limit,
        args.max_split_depth,
        args.max_boxes,
        w_lo,
        w_hi,
        lambda_lo,
        lambda_hi,
        args.quantity,
        args.sign,
    )
    result["script_sha256"] = sha256_file(Path(__file__))
    result["symbolic_audit_sha256"] = sha256_file(
        Path(__file__).with_name("prolate_axis_signed_angle_symbolic_audit.py")
    )
    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{sha256_file(Path(__file__))}  {Path(__file__).name}\n"
        f"{result['symbolic_audit_sha256']}  prolate_axis_signed_angle_symbolic_audit.py\n"
        f"{sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
