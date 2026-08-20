#!/usr/bin/env python3
"""Uniform pole-modulus certificate for item 6.

For u=1-w, split d=1-c at d=u^2.  The inner chart d=u^2 y^2 and
the outer chart d=2 t^2 remove the corner singularity from the A'' integral.
The outer chart is covered by a projective near-corner sector u=r*sqrt(2)*t
and a far sector where t is bounded away from zero.

The resulting compact density bounds imply

    |A_lambda''(1-u)| <= C,  0<=u<=1/64, 1<=lambda<=100,

and therefore

    Psi_lambda(1-u) >= Phi(lambda) - C*u.

The certificate checks that C*2^-24 is below the exact conservative pole
boundary floor 43/5000.
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

from prolate_axis_center_cap_arb import acb_record, closed_interval, regular_angle_derivatives


U_MAX = fmpq(1, 64)
T_SPLIT = fmpq(1, 128)
PHI_FLOOR = fmpq(43, 5000)
DYADIC_EXPONENT = 24
LAMBDA_BANDS = [(1, 2), (2, 3), (3, 5), (5, 10), (10, 20), (20, 40), (40, 70), (70, 100)]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite_ball(value: acb) -> bool:
    text = str(value).lower()
    return not any(token in text for token in ("nan", "inf", "unknown"))


def abs_upper(value: acb) -> arb:
    return abs(value.real).upper()


def compact_density(
    n: acb,
    c: acb,
    radial: acb,
    s2: acb,
    L: acb,
    q: acb,
) -> acb:
    """Cancellation-safe transformed A'' density.

    The products C*gbar and C*(gbar^2+g_vbar) are assembled before
    interval evaluation, eliminating the removable 1/n factors at the
    projective corner n=0.
    """
    root = (radial * s2).sqrt()
    cosine = n / root
    h1, h2 = regular_angle_derivatives(cosine)
    a = L * q / radial
    cwbar = (-c + n * a) / root
    cwwbar = (-2 * c * a + n * (3 * a * a - L / radial)) / root
    return (
        -2 * c * h1 * cwbar
        + n * (h2 * cwbar * cwbar + h1 * cwwbar)
    )


def inner_density(u_value: arb, y_value: arb, lambda_value: arb) -> acb:
    u = acb(u_value)
    y = acb(y_value)
    lam = acb(lambda_value)
    L = lam * lam
    y2 = y * y
    u2 = u * u
    d = u2 * y2
    c = 1 - d
    n = 1 + u * y2 - u2 * y2
    radial = y2 * (2 - u2 * y2) + L * (1 - u * y2) ** 2
    s2 = d * (2 - d) + (1 - d) ** 2 / L
    q = 1 - u * y2
    return u * y * compact_density(n, c, radial, s2, L, q)


def outer_density_from_zr(z: acb, r: acb, lambda_value: arb) -> acb:
    lam = acb(lambda_value)
    L = lam * lam
    z2 = z * z
    c = 1 - z2
    n = r + z - r * z2
    radial = 2 - z2 + L * (r - z) ** 2
    s2 = z2 * (2 - z2) + (1 - z2) ** 2 / L
    q = r - z
    return compact_density(n, c, radial, s2, L, q)


def outer_near_density(t_value: arb, r_value: arb, lambda_value: arb) -> acb:
    t = acb(t_value)
    r = acb(r_value)
    sqrt2 = acb(2).sqrt()
    z = sqrt2 * t
    return sqrt2 * outer_density_from_zr(z, r, lambda_value)


def outer_far_density(t_value: arb, u_value: arb, lambda_value: arb) -> acb:
    t = acb(t_value)
    u = acb(u_value)
    sqrt2 = acb(2).sqrt()
    z = sqrt2 * t
    r = u / z
    return sqrt2 * outer_density_from_zr(z, r, lambda_value)


def cover_region(
    name: str,
    evaluator: Callable[[arb, arb, arb], acb],
    x_lo: fmpq,
    x_hi: fmpq,
    y_lo: fmpq,
    y_hi: fmpq,
    bound_limit: arb,
    max_split_depth: int,
    max_boxes: int,
) -> dict:
    roots = deque(
        (x_lo, x_hi, y_lo, y_hi, fmpq(lo), fmpq(hi), 0)
        for lo, hi in LAMBDA_BANDS
    )
    x_span = x_hi - x_lo
    y_span = y_hi - y_lo
    lambda_span = fmpq(99)
    leaves = []
    terminal = []
    evaluations = 0
    worst_bound = arb(0)
    worst_leaf = None

    while roots:
        if evaluations >= max_boxes:
            while roots:
                xl, xh, yl, yh, ll, lh, depth = roots.popleft()
                terminal.append({
                    "x_lo": str(xl), "x_hi": str(xh),
                    "y_lo": str(yl), "y_hi": str(yh),
                    "lambda_lo": str(ll), "lambda_hi": str(lh),
                    "split_depth": depth,
                    "reason": "max_boxes_exhausted",
                })
            break

        xl, xh, yl, yh, ll, lh, depth = roots.popleft()
        try:
            value = evaluator(
                closed_interval(xl, xh),
                closed_interval(yl, yh),
                closed_interval(ll, lh),
            )
            error = None
        except Exception as exc:
            value = None
            error = f"{type(exc).__name__}: {exc}"

        evaluations += 1
        accepted = bool(
            value is not None
            and finite_ball(value)
            and 0 in value.imag
            and abs_upper(value) < bound_limit
        )
        if accepted and value is not None:
            bound = abs_upper(value)
            record = {
                "x_lo": str(xl), "x_hi": str(xh),
                "y_lo": str(yl), "y_hi": str(yh),
                "lambda_lo": str(ll), "lambda_hi": str(lh),
                "split_depth": depth,
                "density": acb_record(value),
                "absolute_upper": str(bound),
            }
            leaves.append(record)
            if bound > worst_bound:
                worst_bound = bound
                worst_leaf = record
            continue

        if depth < max_split_depth:
            rel = [
                ((xh - xl) / x_span if x_span else fmpq(0), "x"),
                ((yh - yl) / y_span if y_span else fmpq(0), "y"),
                ((lh - ll) / lambda_span, "lambda"),
            ]
            _, axis = max(rel, key=lambda item: item[0])
            if axis == "x":
                mid = (xl + xh) / 2
                roots.append((xl, mid, yl, yh, ll, lh, depth + 1))
                roots.append((mid, xh, yl, yh, ll, lh, depth + 1))
            elif axis == "y":
                mid = (yl + yh) / 2
                roots.append((xl, xh, yl, mid, ll, lh, depth + 1))
                roots.append((xl, xh, mid, yh, ll, lh, depth + 1))
            else:
                mid = (ll + lh) / 2
                roots.append((xl, xh, yl, yh, ll, mid, depth + 1))
                roots.append((xl, xh, yl, yh, mid, lh, depth + 1))
        else:
            record = {
                "x_lo": str(xl), "x_hi": str(xh),
                "y_lo": str(yl), "y_hi": str(yh),
                "lambda_lo": str(ll), "lambda_hi": str(lh),
                "split_depth": depth,
                "reason": "bound_not_certified" if error is None else "evaluation_error",
            }
            if value is not None:
                record["density"] = acb_record(value)
            if error is not None:
                record["error"] = error
            terminal.append(record)

    leaf_volume = sum(
        (
            (fmpq(a["x_hi"]) - fmpq(a["x_lo"]))
            * (fmpq(a["y_hi"]) - fmpq(a["y_lo"]))
            * (fmpq(a["lambda_hi"]) - fmpq(a["lambda_lo"]))
            for a in leaves
        ),
        fmpq(0),
    )
    terminal_volume = sum(
        (
            (fmpq(a["x_hi"]) - fmpq(a["x_lo"]))
            * (fmpq(a["y_hi"]) - fmpq(a["y_lo"]))
            * (fmpq(a["lambda_hi"]) - fmpq(a["lambda_lo"]))
            for a in terminal
        ),
        fmpq(0),
    )
    target_volume = (x_hi - x_lo) * (y_hi - y_lo) * fmpq(99)
    exact_coverage = bool(leaves and not terminal and leaf_volume == target_volume)
    conditions = {
        "all accepted density boxes are finite and real": bool(leaves),
        "all accepted absolute bounds are below the requested limit": bool(leaves),
        "exact rational coverage": exact_coverage,
        "zero terminal boxes": not terminal,
        "volume invariant": leaf_volume + terminal_volume == target_volume,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"
    return {
        "name": name,
        "status": status,
        "conditions": conditions,
        "counts": {
            "evaluations": evaluations,
            "certified_leaves": len(leaves),
            "terminal_boxes": len(terminal),
        },
        "target_box": {
            "x": f"[{x_lo},{x_hi}]",
            "y": f"[{y_lo},{y_hi}]",
            "lambda": "[1,100]",
            "exact_volume": str(target_volume),
        },
        "worst_absolute_upper": str(worst_bound),
        "worst_certified_leaf": worst_leaf,
        "leaves": leaves,
        "terminal": terminal,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--max-split-depth", type=int, default=18)
    parser.add_argument("--max-boxes", type=int, default=65536)
    parser.add_argument("--density-bound-limit", default="140000")
    parser.add_argument("--json", default="prolate_axis_pole_modulus.json")
    args = parser.parse_args()

    ctx.dps = args.dps
    limit = arb(args.density_bound_limit)
    inner = cover_region(
        "inner d=u^2*y^2",
        inner_density,
        fmpq(0), U_MAX,
        fmpq(0), fmpq(1),
        limit,
        args.max_split_depth,
        args.max_boxes,
    )
    near = cover_region(
        "outer near u=r*sqrt(2)*t",
        outer_near_density,
        fmpq(0), T_SPLIT,
        fmpq(0), fmpq(1),
        limit,
        args.max_split_depth,
        args.max_boxes,
    )
    far = cover_region(
        "outer far t>=1/128",
        outer_far_density,
        T_SPLIT, fmpq(1),
        fmpq(0), U_MAX,
        limit,
        args.max_split_depth,
        args.max_boxes,
    )

    inner_bound = arb(inner["worst_absolute_upper"])
    near_bound = arb(near["worst_absolute_upper"])
    far_bound = arb(far["worst_absolute_upper"])
    total_bound = inner_bound + arb(T_SPLIT) * near_bound + arb(1 - T_SPLIT) * far_bound
    epsilon = fmpq(1, 2**DYADIC_EXPONENT)
    modulus_loss = total_bound * arb(epsilon)
    boundary_floor = arb(PHI_FLOOR)

    conditions = {
        "inner compact density cover certified": inner["status"] == "CERTIFIED",
        "outer near compact density cover certified": near["status"] == "CERTIFIED",
        "outer far compact density cover certified": far["status"] == "CERTIFIED",
        "modulus loss at 2^-24 is below boundary floor": modulus_loss < boundary_floor,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"
    result = {
        "status": status,
        "certified_statement": (
            "|A_lambda''(w)|<=C on 1<=lambda<=100, 63/64<=w<1, and "
            "Psi_lambda(1-u)>0 for 0<=u<=2^-24 using Phi>=43/5000."
            if status == "CERTIFIED"
            else None
        ),
        "constants": {
            "u_max": str(U_MAX),
            "t_split": str(T_SPLIT),
            "boundary_floor": str(PHI_FLOOR),
            "dyadic_exponent": DYADIC_EXPONENT,
            "epsilon": str(epsilon),
            "inner_density_bound": str(inner_bound),
            "outer_near_density_bound": str(near_bound),
            "outer_far_density_bound": str(far_bound),
            "A_second_uniform_bound": str(total_bound),
            "modulus_loss": str(modulus_loss),
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
            "A'' is split at d=u^2.  The inner density is integrated over y in "
            "[0,1].  The outer density uses d=2t^2; t in [0,1/128] is resolved "
            "by u=r*sqrt(2)*t, and t in [1/128,1] is nonsingular.  The displayed "
            "density suprema give the stated uniform integral bound."
        ),
        "script_sha256": sha256_file(Path(__file__)),
    }

    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{sha256_file(Path(__file__))}  {Path(__file__).name}\n"
        f"{sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if status == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
