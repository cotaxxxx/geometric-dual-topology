#!/usr/bin/env python3
"""Adaptive Arb certificate for the axial derivative boundary limit.

Certifies the one-sided profile

    Phi(lambda) = lim_{w->1^-} Psi_lambda(w) > 0

on one exact rational lambda block. The boundary formula is a regular 1D
integral obtained after factoring d=1-c from R^2 and combining N*C_w before
setting w=1.
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


def boundary_profile(
    lambda_value: arb,
    tolerance: arb,
    integration_depth: int,
    eval_limit: int,
) -> acb:
    lam = acb(lambda_value)
    lam2 = lam * lam

    def kernel(t: acb, analytic: bool) -> acb:
        c = 2 * t - 1
        d = 1 - c
        r_factor = 2 + (lam2 - 1) * d
        s2 = d * (2 - d) + (1 - d) ** 2 / lam2
        cosine = d.sqrt(analytic=analytic) / (
            r_factor.sqrt(analytic=analytic) * s2.sqrt(analytic=analytic)
        )
        h, h1 = regular_angle_data(cosine)
        weighted_cw = -c - lam2 * d / r_factor
        # The Jacobian dc/dt=2 cancels the factor 1/2 in the profile.
        return -c * h + h1 * cosine * weighted_cw

    return rigorous_integral(
        kernel, tolerance, integration_depth, eval_limit
    )


def certify_block(
    dps: int,
    tolerance_text: str,
    integration_depth: int,
    eval_limit: int,
    max_split_depth: int,
    max_boxes: int,
    lambda_lo: fmpq,
    lambda_hi: fmpq,
) -> dict:
    ctx.dps = dps
    tolerance = arb(tolerance_text)
    queue: deque[tuple[fmpq, fmpq, int]] = deque(
        [(lambda_lo, lambda_hi, 0)]
    )
    leaves: list[dict] = []
    terminal: list[dict] = []
    evaluations = 0
    worst_lower = None
    worst_leaf = None

    while queue:
        if evaluations >= max_boxes:
            while queue:
                lo, hi, depth = queue.popleft()
                terminal.append({
                    "lambda_lo": str(lo),
                    "lambda_hi": str(hi),
                    "split_depth": depth,
                    "reason": "max_boxes_exhausted",
                })
            break

        lo, hi, depth = queue.popleft()
        try:
            value = boundary_profile(
                closed_interval(lo, hi),
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
                "lambda_lo": str(lo),
                "lambda_hi": str(hi),
                "split_depth": depth,
                "Phi": acb_record(value),
            }
            leaves.append(record)
            lower = value.real.lower()
            if worst_lower is None or lower < worst_lower:
                worst_lower = lower
                worst_leaf = record
            continue

        if depth < max_split_depth:
            midpoint = (lo + hi) / 2
            queue.append((lo, midpoint, depth + 1))
            queue.append((midpoint, hi, depth + 1))
        else:
            record = {
                "lambda_lo": str(lo),
                "lambda_hi": str(hi),
                "split_depth": depth,
                "reason": "sign_not_certified" if error is None else "evaluation_error",
            }
            if value is not None:
                record["Phi"] = acb_record(value)
            if error is not None:
                record["error"] = error
            terminal.append(record)

    leaf_length = sum(
        (fmpq(leaf["lambda_hi"]) - fmpq(leaf["lambda_lo"]) for leaf in leaves),
        fmpq(0),
    )
    terminal_length = sum(
        (fmpq(box["lambda_hi"]) - fmpq(box["lambda_lo"]) for box in terminal),
        fmpq(0),
    )
    target_length = lambda_hi - lambda_lo
    partition_ok = leaf_length + terminal_length == target_length
    exact_coverage = bool(leaves and not terminal and leaf_length == target_length)

    conditions = {
        "all accepted leaves have Phi > 0": bool(leaves),
        "exact binary-partition length invariant": partition_ok,
        "exact rational coverage of block": exact_coverage,
        "zero terminal intervals": not terminal,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"

    return {
        "status": status,
        "scope": f"Phi(lambda)>0 for {lambda_lo}<=lambda<={lambda_hi}",
        "definition": "Phi(lambda)=lim_{w->1^-} Psi_lambda(w)",
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
        "interval_evaluation_limit": max_boxes,
        "target_interval": f"[{lambda_lo},{lambda_hi}]",
        "partition": {
            "leaf_length": str(leaf_length),
            "terminal_length": str(terminal_length),
            "target_length": str(target_length),
            "length_invariant": partition_ok,
        },
        "conditions": conditions,
        "counts": {
            "evaluations": evaluations,
            "certified_leaves": len(leaves),
            "terminal_intervals": len(terminal),
        },
        "worst_certified_leaf": worst_leaf,
        "leaves": leaves,
        "terminal": terminal,
        "limitations": (
            "This certifies only the boundary anchor w=1^-; a finite-width "
            "pole-cap transfer is still required."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--tolerance", default="1e-18")
    parser.add_argument("--integration-depth", type=int, default=22)
    parser.add_argument("--eval-limit", type=int, default=200000)
    parser.add_argument("--max-split-depth", type=int, default=12)
    parser.add_argument("--max-boxes", type=int, default=8192)
    parser.add_argument("--lambda-lo", type=int, required=True)
    parser.add_argument("--lambda-hi", type=int, required=True)
    parser.add_argument("--json", required=True)
    args = parser.parse_args()

    if args.lambda_lo < 1 or args.lambda_hi <= args.lambda_lo:
        raise ValueError("require 1 <= lambda-lo < lambda-hi")

    result = certify_block(
        args.dps,
        args.tolerance,
        args.integration_depth,
        args.eval_limit,
        args.max_split_depth,
        args.max_boxes,
        fmpq(args.lambda_lo),
        fmpq(args.lambda_hi),
    )
    result["script_sha256"] = sha256_file(Path(__file__))
    result["symbolic_audit_sha256"] = sha256_file(
        Path(__file__).with_name("prolate_axis_pole_boundary_symbolic_audit.py")
    )
    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{sha256_file(Path(__file__))}  {Path(__file__).name}\n"
        f"{result['symbolic_audit_sha256']}  prolate_axis_pole_boundary_symbolic_audit.py\n"
        f"{sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
