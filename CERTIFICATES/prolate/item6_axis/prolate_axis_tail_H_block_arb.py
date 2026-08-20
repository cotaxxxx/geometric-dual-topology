#!/usr/bin/env python3
"""Adaptive Arb certificate for one compact block of the scaled tail H.

For mu=1/lambda and w>0,

    H(mu,w) = Psi_{1/mu}(w)/(mu*w).

The integrand is regularized before interval arithmetic by removing
h(0)=pi^2/4 using integral_{-1}^1 c dc=0. The moving layer c=w is then split
and each side is mapped to t in [0,1].
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
    """Return h(c)=acos(c)^2 and h'(c), regular at c=1."""
    one = acb(1)
    z = (one - cosine) / 2
    hyper = z.hypgeom_2f1(one / 2, one / 2, acb(3) / 2)
    h = 4 * z * hyper * hyper
    sinc = (-h / 4).hypgeom_0f1(acb(3) / 2)
    h1 = -2 / sinc
    return h, h1


def regular_hbar(cosine: acb) -> acb:
    """Return (h(c)-pi^2/4)/c, regular at c=0 and c=1.

    The asin(c)/c hypergeometric form is regular at c=0 but its
    direct 2F1 evaluation is non-finite at the endpoint c=1.  When
    the cosine ball excludes zero, use the endpoint-regular h(c)
    representation instead.  Only balls containing zero use the
    center series.
    """
    pi = acb(arb.pi())
    if 0 not in cosine:
        h, _ = regular_angle_data(cosine)
        return (h - pi * pi / 4) / cosine
    one = acb(1)
    asin_ratio = (cosine * cosine).hypgeom_2f1(
        one / 2, one / 2, acb(3) / 2
    )
    return -pi * asin_ratio + cosine * asin_ratio * asin_ratio


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


def scaled_tail_H(
    mu_value: arb,
    w_value: arb,
    tolerance: arb,
    integration_depth: int,
    eval_limit: int,
) -> acb:
    mu = acb(mu_value)
    w = acb(w_value)
    mu2 = mu * mu

    def base_integrand(c: acb, analytic: bool) -> acb:
        c2 = c * c
        n = 1 - w * c
        delta = c - w
        p = delta * delta + mu2 * (1 - c2)
        s = 1 - c2 + mu2 * c2
        root = (p * s).sqrt(analytic=analytic)
        cosine = mu * n / root
        _, h1 = regular_angle_data(cosine)
        hbar = regular_hbar(cosine)
        bracket = -c + n * delta / p
        return n * (-c * hbar + h1 * bracket) / (2 * w * root)

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
    mu_lo: fmpq,
    mu_hi: fmpq,
    w_lo: fmpq,
    w_hi: fmpq,
) -> dict:
    ctx.dps = dps
    tolerance = arb(tolerance_text)
    mu_span = mu_hi - mu_lo
    w_span = w_hi - w_lo

    queue: deque[tuple[fmpq, fmpq, fmpq, fmpq, int]] = deque(
        [(mu_lo, mu_hi, w_lo, w_hi, 0)]
    )
    leaves: list[dict] = []
    terminal: list[dict] = []
    evaluations = 0
    worst_lower = None
    worst_leaf = None

    while queue:
        if evaluations >= max_boxes:
            while queue:
                ml, mh, wl, wh, depth = queue.popleft()
                terminal.append({
                    "mu_lo": str(ml),
                    "mu_hi": str(mh),
                    "w_lo": str(wl),
                    "w_hi": str(wh),
                    "split_depth": depth,
                    "reason": "max_boxes_exhausted",
                })
            break

        ml, mh, wl, wh, depth = queue.popleft()
        try:
            value = scaled_tail_H(
                closed_interval(ml, mh),
                closed_interval(wl, wh),
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
                "mu_lo": str(ml),
                "mu_hi": str(mh),
                "w_lo": str(wl),
                "w_hi": str(wh),
                "split_depth": depth,
                "H": acb_record(value),
            }
            leaves.append(record)
            lower = value.real.lower()
            if worst_lower is None or lower < worst_lower:
                worst_lower = lower
                worst_leaf = record
            continue

        if depth < max_split_depth:
            mu_relative = (mh - ml) / mu_span
            w_relative = (wh - wl) / w_span
            if w_relative >= mu_relative:
                midpoint = (wl + wh) / 2
                queue.append((ml, mh, wl, midpoint, depth + 1))
                queue.append((ml, mh, midpoint, wh, depth + 1))
            else:
                midpoint = (ml + mh) / 2
                queue.append((ml, midpoint, wl, wh, depth + 1))
                queue.append((midpoint, mh, wl, wh, depth + 1))
        else:
            record = {
                "mu_lo": str(ml),
                "mu_hi": str(mh),
                "w_lo": str(wl),
                "w_hi": str(wh),
                "split_depth": depth,
                "reason": "sign_not_certified" if error is None else "evaluation_error",
            }
            if value is not None:
                record["H"] = acb_record(value)
            if error is not None:
                record["error"] = error
            terminal.append(record)

    leaf_area = sum(
        (
            (fmpq(leaf["mu_hi"]) - fmpq(leaf["mu_lo"]))
            * (fmpq(leaf["w_hi"]) - fmpq(leaf["w_lo"]))
            for leaf in leaves
        ),
        fmpq(0),
    )
    terminal_area = sum(
        (
            (fmpq(box["mu_hi"]) - fmpq(box["mu_lo"]))
            * (fmpq(box["w_hi"]) - fmpq(box["w_lo"]))
            for box in terminal
        ),
        fmpq(0),
    )
    target_area = mu_span * w_span
    partition_ok = leaf_area + terminal_area == target_area
    exact_coverage = bool(leaves and not terminal and leaf_area == target_area)

    conditions = {
        "all accepted leaves have H > 0": bool(leaves),
        "exact binary-partition area invariant": partition_ok,
        "exact rational coverage of block": exact_coverage,
        "zero terminal boxes": not terminal,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"

    return {
        "status": status,
        "scope": (
            f"H(mu,w)>0 for {mu_lo}<=mu<={mu_hi}, "
            f"{w_lo}<=w<={w_hi}"
        ),
        "derived_statement": (
            "Psi_{1/mu}(w)>0 on the same block because mu*w>0."
        ),
        "environment": {
            "python": platform.python_version(),
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
        "arithmetic": "python-flint Arb/Acb ball arithmetic",
        "integrand": "constant-angle cancellation plus moving split c=w",
        "decimal_precision": dps,
        "integration_tolerance": tolerance_text,
        "integration_depth_limit": integration_depth,
        "integration_eval_limit_per_half": eval_limit,
        "parameter_split_depth_limit": max_split_depth,
        "box_evaluation_limit": max_boxes,
        "target_rectangle": {
            "mu": f"[{mu_lo},{mu_hi}]",
            "w": f"[{w_lo},{w_hi}]",
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
        "worst_certified_leaf": worst_leaf,
        "leaves": leaves,
        "terminal": terminal,
        "limitations": (
            "This is one compact tail slab away from w=0 and w=1. It does not "
            "prove the limit mu->0 or the full item 6 theorem."
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

    result = certify_block(
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
    result["script_sha256"] = sha256_file(Path(__file__))
    result["regularization_audit_sha256"] = sha256_file(
        Path(__file__).with_name("prolate_axis_tail_regularized_symbolic_audit.py")
    )
    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{sha256_file(Path(__file__))}  {Path(__file__).name}\n"
        f"{result['regularization_audit_sha256']}  prolate_axis_tail_regularized_symbolic_audit.py\n"
        f"{sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
