#!/usr/bin/env python3
"""Adaptive Arb certificate for a finite axial center cap.

The script certifies A_lambda''(v)>0 on a rational (v,lambda) rectangle.
Since Psi_lambda(0)=0,

    Psi_lambda(w)/w = integral_0^1 A_lambda''(t w) dt,

strict positivity of A'' on 0<=v<=w0 implies Psi_lambda(w)>0 for
0<w<=w0. All parameter boxes have exact rational endpoints.
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
    hyper = z.hypgeom_2f1(one / 2, one / 2, acb(3) / 2)
    h = 4 * z * hyper * hyper

    # If u^2=h, then
    # sin(u)/u = 0F1(3/2;-h/4),
    # 3(sin(u)-u cos(u))/u^3 = 0F1(5/2;-h/4).
    s = (-h / 4).hypgeom_0f1(acb(3) / 2)
    t = (-h / 4).hypgeom_0f1(acb(5) / 2)
    h1 = -2 / s
    h2 = 2 * t / (3 * s**3)
    return h1, h2


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


def axial_second_derivative(
    v_value: arb,
    lambda_value: arb,
    tolerance: arb,
    integration_depth: int,
    eval_limit: int,
) -> acb:
    """Uniform enclosure of A_lambda''(v) on one real parameter box."""
    v = acb(v_value)
    lam = acb(lambda_value)
    lam2 = lam * lam

    def kernel(t: acb, analytic: bool) -> acb:
        # Map t in [0,1] to c in [-1,1]. The Jacobian 2 cancels the
        # factor 1/2 in the exact second-derivative formula.
        c = 2 * t - 1
        c2 = c * c
        n = 1 - v * c
        r2 = 1 - c2 + lam2 * (c - v) ** 2
        s2 = 1 - c2 + c2 / lam2
        r = r2.sqrt(analytic=analytic)
        s = s2.sqrt(analytic=analytic)
        cosine = n / (r * s)
        h1, h2 = regular_angle_derivatives(cosine)

        g = -c / n + lam2 * (c - v) / r2
        g_v = (
            -c2 / n**2
            - lam2 / r2
            + 2 * lam2**2 * (c - v) ** 2 / r2**2
        )
        cosine_v = cosine * g
        cosine_vv = cosine * (g**2 + g_v)

        return (
            -2 * c * h1 * cosine_v
            + n * (h2 * cosine_v**2 + h1 * cosine_vv)
        )

    return rigorous_integral(
        kernel, tolerance, integration_depth, eval_limit
    )


def rectangles_overlap_interior(a: tuple[fmpq, ...], b: tuple[fmpq, ...]) -> bool:
    aw0, aw1, al0, al1 = a
    bw0, bw1, bl0, bl1 = b
    return max(aw0, bw0) < min(aw1, bw1) and max(al0, bl0) < min(al1, bl1)


def certify(
    dps: int,
    tolerance_text: str,
    integration_depth: int,
    eval_limit: int,
    max_split_depth: int,
    max_boxes: int,
    w_num: int,
    w_den: int,
    lambda_hi: int,
) -> dict:
    ctx.dps = dps
    tolerance = arb(tolerance_text)
    w_cap = fmpq(w_num, w_den)
    lambda_span = fmpq(lambda_hi - 1)

    queue: deque[tuple[fmpq, fmpq, fmpq, fmpq, int]] = deque(
        (fmpq(0), w_cap, fmpq(k), fmpq(k + 1), 0)
        for k in range(1, lambda_hi)
    )
    leaves: list[dict] = []
    leaf_rectangles: list[tuple[fmpq, fmpq, fmpq, fmpq]] = []
    terminal: list[dict] = []
    evaluations = 0
    worst_lower = None
    worst_leaf = None

    while queue:
        if evaluations >= max_boxes:
            while queue:
                w_lo, w_hi, lam_lo, lam_hi_box, split_depth = queue.popleft()
                terminal.append(
                    {
                        "w_lo": str(w_lo),
                        "w_hi": str(w_hi),
                        "lambda_lo": str(lam_lo),
                        "lambda_hi": str(lam_hi_box),
                        "split_depth": split_depth,
                        "reason": "max_boxes_exhausted",
                    }
                )
            break

        w_lo, w_hi, lam_lo, lam_hi_box, split_depth = queue.popleft()
        try:
            value = axial_second_derivative(
                closed_interval(w_lo, w_hi),
                closed_interval(lam_lo, lam_hi_box),
                tolerance,
                integration_depth,
                eval_limit,
            )
            error = None
        except Exception as exc:  # preserve the exact failed box for diagnosis
            value = None
            error = f"{type(exc).__name__}: {exc}"

        evaluations += 1
        positive = bool(
            value is not None and value.real > 0 and 0 in value.imag
        )

        if positive and value is not None:
            record = {
                "w_lo": str(w_lo),
                "w_hi": str(w_hi),
                "lambda_lo": str(lam_lo),
                "lambda_hi": str(lam_hi_box),
                "split_depth": split_depth,
                "A_second": acb_record(value),
            }
            leaves.append(record)
            leaf_rectangles.append((w_lo, w_hi, lam_lo, lam_hi_box))
            lower = value.real.lower()
            if worst_lower is None or lower < worst_lower:
                worst_lower = lower
                worst_leaf = record
            continue

        if split_depth < max_split_depth:
            w_relative = (w_hi - w_lo) / w_cap
            lambda_relative = (lam_hi_box - lam_lo) / lambda_span
            if w_relative >= lambda_relative:
                midpoint = (w_lo + w_hi) / 2
                queue.append(
                    (w_lo, midpoint, lam_lo, lam_hi_box, split_depth + 1)
                )
                queue.append(
                    (midpoint, w_hi, lam_lo, lam_hi_box, split_depth + 1)
                )
            else:
                midpoint = (lam_lo + lam_hi_box) / 2
                queue.append(
                    (w_lo, w_hi, lam_lo, midpoint, split_depth + 1)
                )
                queue.append(
                    (w_lo, w_hi, midpoint, lam_hi_box, split_depth + 1)
                )
        else:
            terminal_record = {
                "w_lo": str(w_lo),
                "w_hi": str(w_hi),
                "lambda_lo": str(lam_lo),
                "lambda_hi": str(lam_hi_box),
                "split_depth": split_depth,
                "reason": "sign_not_certified" if error is None else "evaluation_error",
            }
            if value is not None:
                terminal_record["A_second"] = acb_record(value)
            if error is not None:
                terminal_record["error"] = error
            terminal.append(terminal_record)

    inside_domain = all(
        fmpq(0) <= w0 < w1 <= w_cap
        and fmpq(1) <= l0 < l1 <= fmpq(lambda_hi)
        for w0, w1, l0, l1 in leaf_rectangles
    )
    exact_area = sum(
        ((w1 - w0) * (l1 - l0) for w0, w1, l0, l1 in leaf_rectangles),
        fmpq(0),
    )
    target_area = w_cap * fmpq(lambda_hi - 1)
    pairwise_disjoint = all(
        not rectangles_overlap_interior(leaf_rectangles[i], leaf_rectangles[j])
        for i in range(len(leaf_rectangles))
        for j in range(i + 1, len(leaf_rectangles))
    )
    exact_coverage = bool(
        not terminal
        and leaf_rectangles
        and inside_domain
        and pairwise_disjoint
        and exact_area == target_area
    )

    conditions = {
        "all accepted leaves have A_second > 0": bool(leaves and not terminal),
        "all leaves lie inside the target rectangle": inside_domain,
        "leaf interiors are pairwise disjoint": pairwise_disjoint,
        "exact rational leaf-area sum equals target area": exact_area == target_area,
        "exact rational coverage of center-cap rectangle": exact_coverage,
        "zero terminal boxes": not terminal,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"

    return {
        "status": status,
        "scope": (
            f"A_lambda''(v)>0 for 0<=v<={w_cap}, "
            f"1<=lambda<={lambda_hi}"
        ),
        "derived_conclusion": (
            f"If status is CERTIFIED, Psi_lambda(w)>0 for "
            f"1<=lambda<={lambda_hi} and 0<w<={w_cap}."
        ),
        "environment": {
            "python": platform.python_version(),
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
        "arithmetic": "python-flint Arb/Acb ball arithmetic",
        "decimal_precision": dps,
        "integration_tolerance": tolerance_text,
        "integration_depth_limit": integration_depth,
        "integration_eval_limit": eval_limit,
        "parameter_split_depth_limit": max_split_depth,
        "box_evaluation_limit": max_boxes,
        "target_rectangle": {
            "w": f"[0,{w_cap}]",
            "lambda": f"[1,{lambda_hi}]",
            "exact_area": str(target_area),
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
            "This is only Region C on the stated compact parameter range. "
            "It does not cover the compact interior, pole cap, or tail."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--tolerance", default="1e-18")
    parser.add_argument("--integration-depth", type=int, default=22)
    parser.add_argument("--eval-limit", type=int, default=200000)
    parser.add_argument("--max-split-depth", type=int, default=14)
    parser.add_argument("--max-boxes", type=int, default=16384)
    parser.add_argument("--w-num", type=int, default=1)
    parser.add_argument("--w-den", type=int, default=20)
    parser.add_argument("--lambda-hi", type=int, default=10)
    parser.add_argument("--json", default="prolate_axis_center_cap_arb.json")
    args = parser.parse_args()
    if args.w_num <= 0 or args.w_den <= 0 or args.w_num >= args.w_den:
        raise ValueError("require 0 < w-num/w-den < 1")
    if args.lambda_hi <= 1:
        raise ValueError("lambda-hi must be an integer greater than 1")

    result = certify(
        args.dps,
        args.tolerance,
        args.integration_depth,
        args.eval_limit,
        args.max_split_depth,
        args.max_boxes,
        args.w_num,
        args.w_den,
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
