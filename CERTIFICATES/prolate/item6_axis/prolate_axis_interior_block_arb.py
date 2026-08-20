#!/usr/bin/env python3
"""Adaptive Arb certificate for one compact-interior parameter block.

The script certifies

    Psi_lambda(w) > 0

on one exact rational rectangle contained in

    0 < w < 1,  1 <= lambda < infinity.

The moving layer c=w is preserved by splitting the c-integral there and mapping
both pieces to t in [0,1]:

    c_left  = -1 + (1+w)t,
    c_right =  w + (1-w)t.

This avoids replacing the correlated quantities c-w, 1-wc, and R^2 by
independent interval balls before integration.
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


def regular_angle_data(cosine: acb) -> tuple[acb, acb]:
    """Return h(c)=acos(c)^2 and h'(c), regular at alignment c=1."""
    one = acb(1)
    z = (one - cosine) / 2
    hyper = z.hypgeom_2f1(one / 2, one / 2, acb(3) / 2)
    h = 4 * z * hyper * hyper
    sinc = (-h / 4).hypgeom_0f1(acb(3) / 2)
    h1 = -2 / sinc
    return h, h1


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


def psi_on_box(
    w_value: arb,
    lambda_value: arb,
    tolerance: arb,
    integration_depth: int,
    eval_limit: int,
) -> acb:
    """Uniform enclosure of Psi_lambda(w) on one real parameter box."""
    w = acb(w_value)
    lam = acb(lambda_value)
    lam2 = lam * lam

    def base_integrand(c: acb, analytic: bool) -> acb:
        c2 = c * c
        n = 1 - w * c
        r2 = 1 - c2 + lam2 * (c - w) ** 2
        s2 = 1 - c2 + c2 / lam2
        r = r2.sqrt(analytic=analytic)
        s = s2.sqrt(analytic=analytic)
        cosine = n / (r * s)
        h, h1 = regular_angle_data(cosine)
        cosine_w = cosine * (-c / n + lam2 * (c - w) / r2)
        return (-c * h + n * h1 * cosine_w) / 2

    def left_kernel(t: acb, analytic: bool) -> acb:
        jacobian = 1 + w
        c = -1 + jacobian * t
        return jacobian * base_integrand(c, analytic)

    def right_kernel(t: acb, analytic: bool) -> acb:
        jacobian = 1 - w
        c = w + jacobian * t
        return jacobian * base_integrand(c, analytic)

    left = rigorous_integral(
        left_kernel, tolerance, integration_depth, eval_limit
    )
    right = rigorous_integral(
        right_kernel, tolerance, integration_depth, eval_limit
    )
    return left + right


def certify_block(
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
    worst_lower = None
    worst_leaf = None

    while queue:
        if evaluations >= max_boxes:
            while queue:
                wl, wh, ll, lh, split_depth = queue.popleft()
                terminal.append({
                    "w_lo": str(wl),
                    "w_hi": str(wh),
                    "lambda_lo": str(ll),
                    "lambda_hi": str(lh),
                    "split_depth": split_depth,
                    "reason": "max_boxes_exhausted",
                })
            break

        wl, wh, ll, lh, split_depth = queue.popleft()
        try:
            value = psi_on_box(
                closed_interval(wl, wh),
                closed_interval(ll, lh),
                tolerance,
                integration_depth,
                eval_limit,
            )
            error = None
        except Exception as exc:
            value = None
            error = f"{type(exc).__name__}: {exc}"

        evaluations += 1
        positive = bool(value is not None and value.real > 0 and 0 in value.imag)

        if positive and value is not None:
            record = {
                "w_lo": str(wl),
                "w_hi": str(wh),
                "lambda_lo": str(ll),
                "lambda_hi": str(lh),
                "split_depth": split_depth,
                "Psi": acb_record(value),
            }
            leaves.append(record)
            lower = value.real.lower()
            if worst_lower is None or lower < worst_lower:
                worst_lower = lower
                worst_leaf = record
            continue

        if split_depth < max_split_depth:
            w_relative = (wh - wl) / w_span
            lambda_relative = (lh - ll) / lambda_span
            if w_relative >= lambda_relative:
                midpoint = (wl + wh) / 2
                queue.append((wl, midpoint, ll, lh, split_depth + 1))
                queue.append((midpoint, wh, ll, lh, split_depth + 1))
            else:
                midpoint = (ll + lh) / 2
                queue.append((wl, wh, ll, midpoint, split_depth + 1))
                queue.append((wl, wh, midpoint, lh, split_depth + 1))
        else:
            record = {
                "w_lo": str(wl),
                "w_hi": str(wh),
                "lambda_lo": str(ll),
                "lambda_hi": str(lh),
                "split_depth": split_depth,
                "reason": "sign_not_certified" if error is None else "evaluation_error",
            }
            if value is not None:
                record["Psi"] = acb_record(value)
            if error is not None:
                record["error"] = error
            terminal.append(record)

    exact_leaf_area = sum(
        (
            (fmpq(leaf["w_hi"]) - fmpq(leaf["w_lo"]))
            * (fmpq(leaf["lambda_hi"]) - fmpq(leaf["lambda_lo"]))
            for leaf in leaves
        ),
        fmpq(0),
    )
    exact_terminal_area = sum(
        (
            (fmpq(box["w_hi"]) - fmpq(box["w_lo"]))
            * (fmpq(box["lambda_hi"]) - fmpq(box["lambda_lo"]))
            for box in terminal
        ),
        fmpq(0),
    )
    target_area = w_span * lambda_span
    partition_area_ok = exact_leaf_area + exact_terminal_area == target_area
    exact_coverage = bool(leaves and not terminal and exact_leaf_area == target_area)

    conditions = {
        "all accepted leaves have Psi > 0": bool(leaves),
        "exact binary-partition area invariant": partition_area_ok,
        "exact rational coverage of block": exact_coverage,
        "zero terminal boxes": not terminal,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"

    return {
        "status": status,
        "scope": (
            f"Psi_lambda(w)>0 for {w_lo}<=w<={w_hi}, "
            f"{lambda_lo}<=lambda<={lambda_hi}"
        ),
        "environment": {
            "python": platform.python_version(),
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
        "arithmetic": "python-flint Arb/Acb ball arithmetic",
        "integration_chart": "moving split c=w mapped to two fixed t intervals",
        "decimal_precision": dps,
        "integration_tolerance": tolerance_text,
        "integration_depth_limit": integration_depth,
        "integration_eval_limit_per_half": eval_limit,
        "parameter_split_depth_limit": max_split_depth,
        "box_evaluation_limit": max_boxes,
        "target_rectangle": {
            "w": f"[{w_lo},{w_hi}]",
            "lambda": f"[{lambda_lo},{lambda_hi}]",
            "exact_area": str(target_area),
        },
        "partition": {
            "construction": "exact recursive bisection from one root rectangle",
            "leaf_area": str(exact_leaf_area),
            "terminal_area": str(exact_terminal_area),
            "area_invariant": partition_area_ok,
        },
        "conditions": conditions,
        "counts": {
            "evaluations": evaluations,
            "certified_leaves": len(leaves),
            "terminal_boxes": len(terminal),
        },
        "worst_certified_leaf": worst_leaf,
        "leaves": leaves,
        "terminal": terminal,
        "limitations": (
            "This certifies one compact-interior block only. It does not cover "
            "the center cap, pole cap, aspect-ratio tail, or full item 6 theorem."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--tolerance", default="1e-18")
    parser.add_argument("--integration-depth", type=int, default=22)
    parser.add_argument("--eval-limit", type=int, default=200000)
    parser.add_argument("--max-split-depth", type=int, default=16)
    parser.add_argument("--max-boxes", type=int, default=32768)
    parser.add_argument("--w-lo-num", type=int, required=True)
    parser.add_argument("--w-lo-den", type=int, required=True)
    parser.add_argument("--w-hi-num", type=int, required=True)
    parser.add_argument("--w-hi-den", type=int, required=True)
    parser.add_argument("--lambda-lo", type=int, required=True)
    parser.add_argument("--lambda-hi", type=int, required=True)
    parser.add_argument("--json", required=True)
    args = parser.parse_args()

    w_lo = fmpq(args.w_lo_num, args.w_lo_den)
    w_hi = fmpq(args.w_hi_num, args.w_hi_den)
    lambda_lo = fmpq(args.lambda_lo)
    lambda_hi = fmpq(args.lambda_hi)
    if not (fmpq(0) < w_lo < w_hi < fmpq(1)):
        raise ValueError("require 0 < w-lo < w-hi < 1")
    if not (fmpq(1) <= lambda_lo < lambda_hi):
        raise ValueError("require 1 <= lambda-lo < lambda-hi")

    result = certify_block(
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
    )
    result["script_sha256"] = sha256_file(Path(__file__))
    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{sha256_file(Path(__file__))}  {Path(__file__).name}\n"
        f"{sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
