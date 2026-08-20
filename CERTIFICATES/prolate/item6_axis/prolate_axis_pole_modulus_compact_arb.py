#!/usr/bin/env python3
"""Compact high-capacity Arb certificate for the final pole modulus.

This reuses the audited chart densities from ``prolate_axis_pole_modulus_arb``
but stores only exact aggregate volumes, a deterministic digest of accepted
leaves, the worst accepted leaf, and a bounded list of terminal examples.  It
therefore permits a substantially larger adaptive cover without producing a
multi-gigabyte JSON certificate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from collections import Counter, deque
from pathlib import Path
from typing import Callable

import flint
from flint import acb, arb, ctx, fmpq

import prolate_axis_pole_modulus_arb as base


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def box_volume(box: tuple[fmpq, fmpq, fmpq, fmpq, fmpq, fmpq, int]) -> fmpq:
    xl, xh, yl, yh, ll, lh, _ = box
    return (xh - xl) * (yh - yl) * (lh - ll)


def compact_cover_region(
    name: str,
    evaluator: Callable[[arb, arb, arb], acb],
    x_lo: fmpq,
    x_hi: fmpq,
    y_lo: fmpq,
    y_hi: fmpq,
    bound_limit: arb,
    max_split_depth: int,
    max_boxes: int,
    max_terminal_records: int,
) -> dict:
    roots = deque(
        (x_lo, x_hi, y_lo, y_hi, fmpq(lo), fmpq(hi), 0)
        for lo, hi in base.LAMBDA_BANDS
    )
    x_span = x_hi - x_lo
    y_span = y_hi - y_lo
    lambda_span = fmpq(99)
    target_volume = x_span * y_span * lambda_span

    evaluations = 0
    leaf_count = 0
    terminal_count = 0
    leaf_volume = fmpq(0)
    terminal_volume = fmpq(0)
    worst_bound = arb(0)
    worst_leaf = None
    terminal_examples: list[dict] = []
    reason_counts: Counter[str] = Counter()
    leaf_digest = hashlib.sha256()

    def record_terminal(
        box: tuple[fmpq, fmpq, fmpq, fmpq, fmpq, fmpq, int],
        reason: str,
        value: acb | None = None,
        error: str | None = None,
    ) -> None:
        nonlocal terminal_count, terminal_volume
        xl, xh, yl, yh, ll, lh, depth = box
        terminal_count += 1
        terminal_volume += box_volume(box)
        reason_counts[reason] += 1
        if len(terminal_examples) < max_terminal_records:
            record = {
                "x_lo": str(xl), "x_hi": str(xh),
                "y_lo": str(yl), "y_hi": str(yh),
                "lambda_lo": str(ll), "lambda_hi": str(lh),
                "split_depth": depth,
                "reason": reason,
            }
            if value is not None:
                record["density"] = base.acb_record(value)
            if error is not None:
                record["error"] = error
            terminal_examples.append(record)

    while roots:
        if evaluations >= max_boxes:
            while roots:
                record_terminal(roots.popleft(), "max_boxes_exhausted")
            break

        box = roots.popleft()
        xl, xh, yl, yh, ll, lh, depth = box
        try:
            value = evaluator(
                base.closed_interval(xl, xh),
                base.closed_interval(yl, yh),
                base.closed_interval(ll, lh),
            )
            error = None
        except Exception as exc:  # pragma: no cover - recorded in certificate
            value = None
            error = f"{type(exc).__name__}: {exc}"

        evaluations += 1
        accepted = bool(
            value is not None
            and base.finite_ball(value)
            and 0 in value.imag
            and base.abs_upper(value) < bound_limit
        )
        if accepted and value is not None:
            bound = base.abs_upper(value)
            leaf_count += 1
            leaf_volume += box_volume(box)
            payload = {
                "x_lo": str(xl), "x_hi": str(xh),
                "y_lo": str(yl), "y_hi": str(yh),
                "lambda_lo": str(ll), "lambda_hi": str(lh),
                "split_depth": depth,
                "absolute_upper": str(bound),
            }
            leaf_digest.update(
                json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
            )
            if bound > worst_bound:
                worst_bound = bound
                worst_leaf = {**payload, "density": base.acb_record(value)}
            continue

        if depth < max_split_depth:
            relative = [
                ((xh - xl) / x_span if x_span else fmpq(0), "x"),
                ((yh - yl) / y_span if y_span else fmpq(0), "y"),
                ((lh - ll) / lambda_span, "lambda"),
            ]
            _, axis = max(relative, key=lambda item: item[0])
            if axis == "x":
                midpoint = (xl + xh) / 2
                roots.append((xl, midpoint, yl, yh, ll, lh, depth + 1))
                roots.append((midpoint, xh, yl, yh, ll, lh, depth + 1))
            elif axis == "y":
                midpoint = (yl + yh) / 2
                roots.append((xl, xh, yl, midpoint, ll, lh, depth + 1))
                roots.append((xl, xh, midpoint, yh, ll, lh, depth + 1))
            else:
                midpoint = (ll + lh) / 2
                roots.append((xl, xh, yl, yh, ll, midpoint, depth + 1))
                roots.append((xl, xh, yl, yh, midpoint, lh, depth + 1))
        else:
            record_terminal(
                box,
                "bound_not_certified" if error is None else "evaluation_error",
                value,
                error,
            )

    volume_invariant = leaf_volume + terminal_volume == target_volume
    exact_coverage = bool(leaf_count and terminal_count == 0 and leaf_volume == target_volume)
    conditions = {
        "at least one accepted density box": leaf_count > 0,
        "accepted boxes satisfy the requested absolute bound": leaf_count > 0,
        "exact rational coverage": exact_coverage,
        "zero terminal boxes": terminal_count == 0,
        "volume invariant": volume_invariant,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"
    return {
        "name": name,
        "status": status,
        "requested_absolute_bound": str(bound_limit),
        "conditions": conditions,
        "counts": {
            "evaluations": evaluations,
            "certified_leaves": leaf_count,
            "terminal_boxes": terminal_count,
        },
        "target_box": {
            "x": f"[{x_lo},{x_hi}]",
            "y": f"[{y_lo},{y_hi}]",
            "lambda": "[1,100]",
            "exact_volume": str(target_volume),
        },
        "partition": {
            "leaf_volume": str(leaf_volume),
            "terminal_volume": str(terminal_volume),
            "volume_invariant": volume_invariant,
        },
        "worst_absolute_upper": str(worst_bound),
        "worst_certified_leaf": worst_leaf,
        "accepted_leaf_digest_sha256": leaf_digest.hexdigest(),
        "terminal_reason_counts": dict(sorted(reason_counts.items())),
        "terminal_examples": terminal_examples,
        "compact_output": (
            "Accepted leaves are represented by exact aggregate volume, count, "
            "deterministic SHA-256 stream digest, and the worst accepted leaf."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--max-split-depth", type=int, default=24)
    parser.add_argument("--max-boxes", type=int, default=1_000_000)
    parser.add_argument("--inner-density-limit", default="100")
    parser.add_argument("--outer-density-limit", default="144100")
    parser.add_argument("--max-terminal-records", type=int, default=100)
    parser.add_argument("--json", default="prolate_axis_pole_modulus_compact.json")
    args = parser.parse_args()

    ctx.dps = args.dps
    inner_limit = arb(args.inner_density_limit)
    outer_limit = arb(args.outer_density_limit)

    inner = compact_cover_region(
        "inner d=u^2*y^2",
        base.inner_density,
        fmpq(0), base.U_MAX,
        fmpq(0), fmpq(1),
        inner_limit,
        args.max_split_depth,
        args.max_boxes,
        args.max_terminal_records,
    )
    near = compact_cover_region(
        "outer near u=r*sqrt(2)*t",
        base.outer_near_density,
        fmpq(0), base.T_SPLIT,
        fmpq(0), fmpq(1),
        outer_limit,
        args.max_split_depth,
        args.max_boxes,
        args.max_terminal_records,
    )
    far = compact_cover_region(
        "outer far t>=1/128",
        base.outer_far_density,
        base.T_SPLIT, fmpq(1),
        fmpq(0), base.U_MAX,
        outer_limit,
        args.max_split_depth,
        args.max_boxes,
        args.max_terminal_records,
    )

    inner_bound = arb(inner["worst_absolute_upper"])
    near_bound = arb(near["worst_absolute_upper"])
    far_bound = arb(far["worst_absolute_upper"])
    total_bound = (
        inner_bound
        + arb(base.T_SPLIT) * near_bound
        + arb(1 - base.T_SPLIT) * far_bound
    )
    epsilon = fmpq(1, 2**base.DYADIC_EXPONENT)
    modulus_loss = total_bound * arb(epsilon)
    boundary_floor = arb(base.PHI_FLOOR)

    conditions = {
        "inner compact density cover certified": inner["status"] == "CERTIFIED",
        "outer near compact density cover certified": near["status"] == "CERTIFIED",
        "outer far compact density cover certified": far["status"] == "CERTIFIED",
        "modulus loss at 2^-24 is below boundary floor": modulus_loss < boundary_floor,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"
    script_path = Path(__file__)
    base_path = script_path.with_name("prolate_axis_pole_modulus_arb.py")
    audit_path = script_path.with_name("prolate_axis_pole_modulus_symbolic_audit.py")
    result = {
        "status": status,
        "certified_statement": (
            "|A_lambda''(w)|<=C on 1<=lambda<=100, 63/64<=w<1, and "
            "Psi_lambda(1-u)>0 for 0<=u<=2^-24 using Phi>=43/5000."
            if status == "CERTIFIED" else None
        ),
        "constants": {
            "u_max": str(base.U_MAX),
            "t_split": str(base.T_SPLIT),
            "boundary_floor": str(base.PHI_FLOOR),
            "dyadic_exponent": base.DYADIC_EXPONENT,
            "epsilon": str(epsilon),
            "inner_requested_limit": str(inner_limit),
            "outer_requested_limit": str(outer_limit),
            "inner_density_bound": str(inner_bound),
            "outer_near_density_bound": str(near_bound),
            "outer_far_density_bound": str(far_bound),
            "A_second_uniform_bound": str(total_bound),
            "modulus_loss": str(modulus_loss),
            "absolute_safe_ceiling": str(base.PHI_FLOOR * 2**base.DYADIC_EXPONENT),
        },
        "conditions": conditions,
        "regions": [inner, near, far],
        "environment": {
            "python": platform.python_version(),
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
        "arithmetic": "python-flint Arb/Acb ball arithmetic plus exact rational coverage",
        "derivation": (
            "The same three audited compact charts as the full-output driver are "
            "used. The compact driver changes only certificate serialization and "
            "adaptive capacity, not the evaluated density formulas."
        ),
        "script_sha256": sha256_file(script_path),
        "base_driver_sha256": sha256_file(base_path),
        "symbolic_audit_sha256": sha256_file(audit_path),
    }

    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{result['script_sha256']}  {script_path.name}\n"
        f"{result['base_driver_sha256']}  {base_path.name}\n"
        f"{result['symbolic_audit_sha256']}  {audit_path.name}\n"
        f"{sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if status == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
