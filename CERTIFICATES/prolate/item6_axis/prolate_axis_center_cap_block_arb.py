#!/usr/bin/env python3
"""Adaptive Arb certificate for one exact lambda block of the axial center cap.

The script certifies A_lambda''(v)>0 on
    0 <= v <= w_cap, lambda_lo <= lambda <= lambda_hi.
It imports the audited kernel from prolate_axis_center_cap_arb.py and preserves
an exact binary-partition invariant for the parameter rectangle.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from collections import deque
from pathlib import Path

import flint
from flint import arb, ctx, fmpq

from prolate_axis_center_cap_arb import (
    acb_record,
    axial_second_derivative,
    closed_interval,
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def certify_block(
    dps: int,
    tolerance_text: str,
    integration_depth: int,
    eval_limit: int,
    max_split_depth: int,
    max_boxes: int,
    w_num: int,
    w_den: int,
    lambda_lo_int: int,
    lambda_hi_int: int,
) -> dict:
    ctx.dps = dps
    tolerance = arb(tolerance_text)
    w_cap = fmpq(w_num, w_den)
    lambda_lo = fmpq(lambda_lo_int)
    lambda_hi = fmpq(lambda_hi_int)
    lambda_span = lambda_hi - lambda_lo

    queue: deque[tuple[fmpq, fmpq, fmpq, fmpq, int]] = deque(
        [(fmpq(0), w_cap, lambda_lo, lambda_hi, 0)]
    )
    leaves: list[dict] = []
    terminal: list[dict] = []
    evaluations = 0
    worst_lower = None
    worst_leaf = None

    while queue:
        if evaluations >= max_boxes:
            while queue:
                v_lo, v_hi, lam_lo, lam_hi, split_depth = queue.popleft()
                terminal.append({
                    "v_lo": str(v_lo),
                    "v_hi": str(v_hi),
                    "lambda_lo": str(lam_lo),
                    "lambda_hi": str(lam_hi),
                    "split_depth": split_depth,
                    "reason": "max_boxes_exhausted",
                })
            break

        v_lo, v_hi, lam_lo, lam_hi, split_depth = queue.popleft()
        try:
            value = axial_second_derivative(
                closed_interval(v_lo, v_hi),
                closed_interval(lam_lo, lam_hi),
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
                "v_lo": str(v_lo),
                "v_hi": str(v_hi),
                "lambda_lo": str(lam_lo),
                "lambda_hi": str(lam_hi),
                "split_depth": split_depth,
                "A_second": acb_record(value),
            }
            leaves.append(record)
            lower = value.real.lower()
            if worst_lower is None or lower < worst_lower:
                worst_lower = lower
                worst_leaf = record
            continue

        if split_depth < max_split_depth:
            v_relative = (v_hi - v_lo) / w_cap
            lambda_relative = (lam_hi - lam_lo) / lambda_span
            if v_relative >= lambda_relative:
                midpoint = (v_lo + v_hi) / 2
                queue.append((v_lo, midpoint, lam_lo, lam_hi, split_depth + 1))
                queue.append((midpoint, v_hi, lam_lo, lam_hi, split_depth + 1))
            else:
                midpoint = (lam_lo + lam_hi) / 2
                queue.append((v_lo, v_hi, lam_lo, midpoint, split_depth + 1))
                queue.append((v_lo, v_hi, midpoint, lam_hi, split_depth + 1))
        else:
            record = {
                "v_lo": str(v_lo),
                "v_hi": str(v_hi),
                "lambda_lo": str(lam_lo),
                "lambda_hi": str(lam_hi),
                "split_depth": split_depth,
                "reason": "sign_not_certified" if error is None else "evaluation_error",
            }
            if value is not None:
                record["A_second"] = acb_record(value)
            if error is not None:
                record["error"] = error
            terminal.append(record)

    exact_leaf_area = sum(
        (
            (fmpq(leaf["v_hi"]) - fmpq(leaf["v_lo"]))
            * (fmpq(leaf["lambda_hi"]) - fmpq(leaf["lambda_lo"]))
            for leaf in leaves
        ),
        fmpq(0),
    )
    exact_terminal_area = sum(
        (
            (fmpq(box["v_hi"]) - fmpq(box["v_lo"]))
            * (fmpq(box["lambda_hi"]) - fmpq(box["lambda_lo"]))
            for box in terminal
        ),
        fmpq(0),
    )
    target_area = w_cap * lambda_span
    partition_area_ok = exact_leaf_area + exact_terminal_area == target_area
    exact_coverage = bool(leaves and not terminal and exact_leaf_area == target_area)

    conditions = {
        "all accepted leaves have A_second > 0": bool(leaves),
        "exact binary-partition area invariant": partition_area_ok,
        "exact rational coverage of block": exact_coverage,
        "zero terminal boxes": not terminal,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"

    return {
        "status": status,
        "scope": (
            f"A_lambda''(v)>0 for 0<=v<={w_cap}, "
            f"{lambda_lo}<=lambda<={lambda_hi}"
        ),
        "derived_conclusion": (
            f"If status is CERTIFIED, Psi_lambda(w)>0 for "
            f"{lambda_lo}<=lambda<={lambda_hi} and 0<w<={w_cap}."
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
            "v": f"[0,{w_cap}]",
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
            "This certifies one Region C lambda block only. "
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
    parser.add_argument("--lambda-lo", type=int, required=True)
    parser.add_argument("--lambda-hi", type=int, required=True)
    parser.add_argument("--json", required=True)
    args = parser.parse_args()

    if args.w_num <= 0 or args.w_den <= 0 or args.w_num >= args.w_den:
        raise ValueError("require 0 < w-num/w-den < 1")
    if args.lambda_lo < 1 or args.lambda_hi <= args.lambda_lo:
        raise ValueError("require 1 <= lambda-lo < lambda-hi")

    result = certify_block(
        args.dps,
        args.tolerance,
        args.integration_depth,
        args.eval_limit,
        args.max_split_depth,
        args.max_boxes,
        args.w_num,
        args.w_den,
        args.lambda_lo,
        args.lambda_hi,
    )
    result["script_sha256"] = sha256_file(Path(__file__))
    result["kernel_script_sha256"] = sha256_file(
        Path(__file__).with_name("prolate_axis_center_cap_arb.py")
    )
    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{sha256_file(Path(__file__))}  {Path(__file__).name}\n"
        f"{result['kernel_script_sha256']}  prolate_axis_center_cap_arb.py\n"
        f"{sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
