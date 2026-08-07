#!/usr/bin/env python3
"""Adaptive Arb certificate for Q_parallel(lambda)>0 on a compact interval.

Default target: 1 <= lambda <= 10. The lambda domain is partitioned with
exact rational endpoints. Each accepted leaf is a uniform Acb integral over
the full lambda interval of that leaf.
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


def regular_angle_derivatives(cosine: acb) -> tuple[acb, acb]:
    """Return h'(c), h''(c) for h(c)=acos(c)^2, regular at c=1."""
    one = acb(1)
    z = (one - cosine) / 2
    H = z.hypgeom_2f1(one / 2, one / 2, acb(3) / 2)
    h = 4 * z * H * H

    # S=sin(u)/u and T=3(sin(u)-u cos(u))/u^3 for u^2=h.
    S = (-h / 4).hypgeom_0f1(acb(3) / 2)
    T = (-h / 4).hypgeom_0f1(acb(5) / 2)
    h1 = -2 / S
    h2 = 2 * T / (3 * S**3)
    return h1, h2


def rigorous_integral(
    kernel: Callable[[acb, bool], acb],
    tol: arb,
    depth: int,
    eval_limit: int,
) -> acb:
    return acb.integral(
        kernel,
        0,
        1,
        abs_tol=tol,
        rel_tol=tol,
        depth_limit=depth,
        eval_limit=eval_limit,
    )


def q_parallel(
    lam_value: arb,
    tol: arb,
    integration_depth: int,
    eval_limit: int,
) -> acb:
    """Uniform enclosure of Q_parallel over one real lambda ball."""
    lam = acb(lam_value)
    lam2 = lam * lam

    def kernel(c: acb, analytic: bool) -> acb:
        c2 = c * c
        one_minus_c2 = 1 - c2
        ell = 1 + (lam2 - 1) * c2
        w2 = lam2 * one_minus_c2 + c2
        sqrt_ell = ell.sqrt(analytic=analytic)
        sqrt_w2 = w2.sqrt(analytic=analytic)
        cosine = lam / (sqrt_ell * sqrt_w2)
        h1, h2 = regular_angle_derivatives(cosine)

        cosine_w = (
            lam
            * c
            * (lam2 - 1)
            * one_minus_c2
            / (ell * sqrt_ell * sqrt_w2)
        )
        cosine_ww = (
            lam**3
            * one_minus_c2
            * (2 * (lam2 - 1) * c2 - 1)
            / (ell**2 * sqrt_ell * sqrt_w2)
        )

        # The center kernel is even in c. Integrating this expression on
        # [0,1] already includes the factor 2 from the full interval.
        return (
            -2 * c * h1 * cosine_w
            + h2 * cosine_w**2
            + h1 * cosine_ww
        )

    return rigorous_integral(kernel, tol, integration_depth, eval_limit)


def certify(
    dps: int,
    tolerance: str,
    integration_depth: int,
    eval_limit: int,
    max_split_depth: int,
    max_boxes: int,
    lambda_hi: int,
) -> dict:
    ctx.dps = dps
    tol = arb(tolerance)

    queue: deque[tuple[fmpq, fmpq, int]] = deque(
        (fmpq(k), fmpq(k + 1), 0) for k in range(1, lambda_hi)
    )
    leaves: list[dict] = []
    leaf_segments: list[tuple[fmpq, fmpq]] = []
    terminal: list[dict] = []
    evaluations = 0
    worst_lower = None
    worst_leaf = None

    while queue:
        if evaluations >= max_boxes:
            while queue:
                lo, hi, split_depth = queue.popleft()
                terminal.append(
                    {
                        "lambda_lo": str(lo),
                        "lambda_hi": str(hi),
                        "split_depth": split_depth,
                        "reason": "max_boxes_exhausted",
                    }
                )
            break

        lo, hi, split_depth = queue.popleft()
        value = q_parallel(
            closed_interval(lo, hi), tol, integration_depth, eval_limit
        )
        evaluations += 1
        positive = bool(value.real > 0 and 0 in value.imag)

        if positive:
            record = {
                "lambda_lo": str(lo),
                "lambda_hi": str(hi),
                "split_depth": split_depth,
                "Q_parallel": acb_record(value),
            }
            leaves.append(record)
            leaf_segments.append((lo, hi))
            lower = value.real.lower()
            if worst_lower is None or lower < worst_lower:
                worst_lower = lower
                worst_leaf = record
            continue

        if split_depth < max_split_depth:
            mid = (lo + hi) / 2
            queue.append((lo, mid, split_depth + 1))
            queue.append((mid, hi, split_depth + 1))
        else:
            terminal.append(
                {
                    "lambda_lo": str(lo),
                    "lambda_hi": str(hi),
                    "split_depth": split_depth,
                    "reason": "sign_not_certified",
                    "Q_parallel": acb_record(value),
                }
            )

    leaf_segments.sort(key=lambda pair: pair[0])
    coverage = bool(
        not terminal
        and leaf_segments
        and leaf_segments[0][0] == fmpq(1)
        and leaf_segments[-1][1] == fmpq(lambda_hi)
        and all(
            leaf_segments[index][1] == leaf_segments[index + 1][0]
            for index in range(len(leaf_segments) - 1)
        )
    )

    conditions = {
        "all accepted leaves have Q_parallel > 0": bool(leaves and not terminal),
        "exact rational coverage of compact lambda interval": coverage,
        "zero terminal boxes": not terminal,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"

    return {
        "status": status,
        "scope": f"Q_parallel(lambda)>0 for 1<=lambda<={lambda_hi}",
        "environment": {
            "python": platform.python_version(),
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
        "arithmetic": "python-flint Arb/Acb ball arithmetic",
        "decimal_precision": dps,
        "integration_tolerance": tolerance,
        "integration_depth_limit": integration_depth,
        "integration_eval_limit": eval_limit,
        "lambda_split_depth_limit": max_split_depth,
        "box_evaluation_limit": max_boxes,
        "conditions": conditions,
        "counts": {
            "evaluations": evaluations,
            "certified_leaves": len(leaves),
            "terminal_boxes": len(terminal),
        },
        "worst_certified_leaf": worst_leaf,
        "leaves": leaves,
        "terminal": terminal,
        "certified_conclusion": (
            f"If status is CERTIFIED, Q_parallel(lambda)>0 uniformly on "
            f"the exact compact interval [1,{lambda_hi}]."
        ),
        "limitations": (
            "This certifies only the center Hessian on a compact lambda "
            "interval. It does not certify a finite-w center cap or the "
            "unbounded aspect-ratio tail."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--tolerance", default="1e-18")
    parser.add_argument("--integration-depth", type=int, default=22)
    parser.add_argument("--eval-limit", type=int, default=200000)
    parser.add_argument("--max-split-depth", type=int, default=8)
    parser.add_argument("--max-boxes", type=int, default=4096)
    parser.add_argument("--lambda-hi", type=int, default=10)
    parser.add_argument(
        "--json", default="prolate_axis_center_compact_arb.json"
    )
    args = parser.parse_args()
    if args.lambda_hi <= 1:
        raise ValueError("lambda-hi must be an integer greater than 1")

    result = certify(
        args.dps,
        args.tolerance,
        args.integration_depth,
        args.eval_limit,
        args.max_split_depth,
        args.max_boxes,
        args.lambda_hi,
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
