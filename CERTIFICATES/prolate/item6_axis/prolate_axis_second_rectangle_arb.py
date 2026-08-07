#!/usr/bin/env python3
"""Adaptive Arb certificate for one signed A_lambda'' rectangle.

The script certifies either A_lambda''(w)>0 or A_lambda''(w)<0 on one exact
rational (w,lambda) rectangle. It reuses the audited second-derivative kernel
from prolate_axis_center_cap_arb.py and preserves an exact recursive
binary-partition invariant.
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
                    "w_lo": str(wl),
                    "w_hi": str(wh),
                    "lambda_lo": str(ll),
                    "lambda_hi": str(lh),
                    "split_depth": depth,
                    "reason": "max_boxes_exhausted",
                })
            break

        wl, wh, ll, lh, depth = queue.popleft()
        try:
            value = axial_second_derivative(
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
        accepted = bool(
            value is not None
            and 0 in value.imag
            and ((sign == "positive" and value.real > 0)
                 or (sign == "negative" and value.real < 0))
        )

        if accepted and value is not None:
            record = {
                "w_lo": str(wl),
                "w_hi": str(wh),
                "lambda_lo": str(ll),
                "lambda_hi": str(lh),
                "split_depth": depth,
                "A_second": acb_record(value),
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
                "w_lo": str(wl),
                "w_hi": str(wh),
                "lambda_lo": str(ll),
                "lambda_hi": str(lh),
                "split_depth": depth,
                "reason": "sign_not_certified" if error is None else "evaluation_error",
            }
            if value is not None:
                record["A_second"] = acb_record(value)
            if error is not None:
                record["error"] = error
            terminal.append(record)

    leaf_area = sum(
        (
            (fmpq(leaf["w_hi"]) - fmpq(leaf["w_lo"]))
            * (fmpq(leaf["lambda_hi"]) - fmpq(leaf["lambda_lo"]))
            for leaf in leaves
        ),
        fmpq(0),
    )
    terminal_area = sum(
        (
            (fmpq(box["w_hi"]) - fmpq(box["w_lo"]))
            * (fmpq(box["lambda_hi"]) - fmpq(box["lambda_lo"]))
            for box in terminal
        ),
        fmpq(0),
    )
    target_area = w_span * lambda_span
    partition_ok = leaf_area + terminal_area == target_area
    exact_coverage = bool(leaves and not terminal and leaf_area == target_area)
    comparison = "> 0" if sign == "positive" else "< 0"

    conditions = {
        f"all accepted leaves have A_second {comparison}": bool(leaves),
        "exact binary-partition area invariant": partition_ok,
        "exact rational coverage of block": exact_coverage,
        "zero terminal boxes": not terminal,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"

    return {
        "status": status,
        "certified_sign": sign,
        "scope": (
            f"A_lambda''(w){comparison} for {w_lo}<=w<={w_hi}, "
            f"{lambda_lo}<=lambda<={lambda_hi}"
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
        "limitations": (
            "This is one signed finite rectangle. Global conclusions require the "
            "strict grid combiner and the relevant center or pole transfer lemma."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=45)
    parser.add_argument("--tolerance", default="1e-16")
    parser.add_argument("--integration-depth", type=int, default=20)
    parser.add_argument("--eval-limit", type=int, default=100000)
    parser.add_argument("--max-split-depth", type=int, default=15)
    parser.add_argument("--max-boxes", type=int, default=16384)
    parser.add_argument("--w-lo-num", type=int, required=True)
    parser.add_argument("--w-lo-den", type=int, required=True)
    parser.add_argument("--w-hi-num", type=int, required=True)
    parser.add_argument("--w-hi-den", type=int, required=True)
    parser.add_argument("--lambda-lo", type=int, required=True)
    parser.add_argument("--lambda-hi", type=int, required=True)
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
        args.sign,
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
