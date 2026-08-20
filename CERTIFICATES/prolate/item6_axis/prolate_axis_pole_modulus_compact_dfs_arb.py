#!/usr/bin/env python3
"""Memory-bounded depth-first compact cover for the item-6 pole modulus.

This module has the same mathematical acceptance, splitting, exact-volume,
and terminal-record rules as ``prolate_axis_pole_modulus_compact_arb``.  The
only change is traversal order: unresolved siblings are kept on a LIFO stack
instead of a breadth-first queue.  Frontier storage is therefore O(depth)
per lambda root rather than proportional to an entire refinement layer.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from typing import Callable

from flint import acb, arb, fmpq

import prolate_axis_pole_modulus_arb as base


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
    stack = [
        (x_lo, x_hi, y_lo, y_hi, fmpq(lo), fmpq(hi), 0)
        for lo, hi in reversed(base.LAMBDA_BANDS)
    ]
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
    maximum_frontier = len(stack)

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
                "x_lo": str(xl),
                "x_hi": str(xh),
                "y_lo": str(yl),
                "y_hi": str(yh),
                "lambda_lo": str(ll),
                "lambda_hi": str(lh),
                "split_depth": depth,
                "reason": reason,
            }
            if value is not None:
                record["density"] = base.acb_record(value)
            if error is not None:
                record["error"] = error
            terminal_examples.append(record)

    while stack:
        maximum_frontier = max(maximum_frontier, len(stack))
        if evaluations >= max_boxes:
            while stack:
                record_terminal(stack.pop(), "max_boxes_exhausted")
            break

        box = stack.pop()
        xl, xh, yl, yh, ll, lh, depth = box
        try:
            value = evaluator(
                base.closed_interval(xl, xh),
                base.closed_interval(yl, yh),
                base.closed_interval(ll, lh),
            )
            error = None
        except Exception as exc:  # recorded in certificate
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
                "x_lo": str(xl),
                "x_hi": str(xh),
                "y_lo": str(yl),
                "y_hi": str(yh),
                "lambda_lo": str(ll),
                "lambda_hi": str(lh),
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
                right = (midpoint, xh, yl, yh, ll, lh, depth + 1)
                left = (xl, midpoint, yl, yh, ll, lh, depth + 1)
            elif axis == "y":
                midpoint = (yl + yh) / 2
                right = (xl, xh, midpoint, yh, ll, lh, depth + 1)
                left = (xl, xh, yl, midpoint, ll, lh, depth + 1)
            else:
                midpoint = (ll + lh) / 2
                right = (xl, xh, yl, yh, midpoint, lh, depth + 1)
                left = (xl, xh, yl, yh, ll, midpoint, depth + 1)
            stack.append(right)
            stack.append(left)
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
            "maximum_frontier": maximum_frontier,
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
        "traversal": "deterministic depth-first, left child before right child",
        "compact_output": (
            "Accepted leaves are represented by exact aggregate volume, count, "
            "deterministic SHA-256 stream digest, and the worst accepted leaf."
        ),
    }
