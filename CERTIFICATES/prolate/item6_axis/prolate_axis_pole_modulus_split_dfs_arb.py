#!/usr/bin/env python3
"""Four-region, budget-adaptive pole-modulus certificate for prolate item 6.

The outer-far chart is split into three exact t slabs:

* junction: [1/128, 1/64], evaluated by the audited real signed-angle formula;
* middle:   [1/64, 255/256], evaluated by the factored cosine formula;
* endpoint: [255/256, 1], evaluated by the monotone signed endpoint formula.

The split isolates the only remaining nonfinite enclosures, which occurred at
the outer-near/outer-far interface. All slabs use the same strict residual
uniform far ceiling derived from the certified inner and outer-near bounds.
The final weighted total is checked directly using exact rational slab widths.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import flint
from flint import acb, arb, ctx, fmpq

import prolate_axis_pole_modulus_arb as base
import prolate_axis_pole_modulus_compact_dfs_arb as compact
import prolate_axis_pole_modulus_endpoint_arb as endpoint

JUNCTION_SPLIT = fmpq(1, 64)
ENDPOINT_SPLIT = endpoint.SIGNED_ENDPOINT_T


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def signed_junction_density(
    t_value: arb,
    u_value: arb,
    lambda_value: arb,
) -> acb:
    """Exact signed-angle transformed A'' density near t=1/128.

    On this slab rho is bounded away from zero, so direct nonnegative Arb
    assembly is stable. N=1-(1-u)(1-2t^2)>0 on the full slab, fixing the real
    atan branch without an oriented-angle ambiguity.
    """
    t = arb(t_value)
    u = arb(u_value)
    lam = arb(lambda_value)
    one = arb(1)
    two = arb(2)
    four = arb(4)

    if t.lower() < arb(base.T_SPLIT) or t.upper() > arb(JUNCTION_SPLIT):
        raise ValueError(f"junction evaluator received t outside slab: {t}")
    if lam.lower() <= 0:
        raise ValueError(f"junction evaluator requires lambda>0: {lam}")

    t2 = endpoint.positive_square_hull(t)
    z2 = two * t2
    c = one - z2
    w = one - u
    L = endpoint.positive_square_hull(lam)

    rho2 = endpoint.known_nonnegative_hull(
        four * t2 * (one - t) * (one + t)
    )
    if rho2.lower() <= 0:
        raise ValueError(f"junction rho^2 is not strictly positive: {rho2}")
    rho = rho2.sqrt()

    N = u + z2 - u * z2
    if N.lower() <= 0:
        raise ValueError(f"junction N is not strictly positive: {N}")
    q = u - z2
    q_abs = endpoint.known_nonnegative_hull(abs(q))
    q2 = endpoint.positive_square_hull(q_abs)
    R2 = rho2 + L * q2
    if R2.lower() <= 0:
        raise ValueError(f"junction R^2 is not strictly positive: {R2}")

    cross = rho * (lam * w - (lam - one / lam) * c)
    delta = acb((cross / N).atan())
    delta_w = acb(lam * rho / R2)
    delta_ww = acb(two * lam * L * rho * q / (R2 * R2))
    transformed = (
        -2 * acb(c) * delta * delta_w
        + acb(N) * (delta_w * delta_w + delta * delta_ww)
    )
    return acb(four * t) * transformed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--max-split-depth", type=int, default=24)
    parser.add_argument("--max-boxes", type=int, default=3_000_000)
    parser.add_argument("--inner-density-limit", default="100")
    parser.add_argument("--outer-near-density-limit", default="144100")
    parser.add_argument("--max-terminal-records", type=int, default=100)
    parser.add_argument("--json", default="prolate_axis_pole_modulus_split.json")
    args = parser.parse_args()

    ctx.dps = args.dps
    inner_limit = arb(args.inner_density_limit)
    near_limit = arb(args.outer_near_density_limit)
    safe_ceiling_q = base.PHI_FLOOR * 2**base.DYADIC_EXPONENT
    safe_ceiling = arb(safe_ceiling_q)

    audit = endpoint.endpoint_factorization_audit()
    audit["checks"] = {key: bool(value) for key, value in audit["checks"].items()}
    audit["status"] = "PASSED" if all(audit["checks"].values()) else "FAILED"
    if audit["status"] != "PASSED":
        print(json.dumps(audit, indent=2))
        raise SystemExit(1)

    inner = compact.compact_cover_region(
        "inner d=u^2*y^2",
        base.inner_density,
        fmpq(0), base.U_MAX,
        fmpq(0), fmpq(1),
        inner_limit,
        args.max_split_depth,
        args.max_boxes,
        args.max_terminal_records,
    )
    near = compact.compact_cover_region(
        "outer near u=r*sqrt(2)*t",
        base.outer_near_density,
        fmpq(0), base.T_SPLIT,
        fmpq(0), fmpq(1),
        near_limit,
        args.max_split_depth,
        args.max_boxes,
        args.max_terminal_records,
    )

    inner_bound = arb(inner["worst_absolute_upper"])
    near_bound = arb(near["worst_absolute_upper"])
    residual_far_ceiling = (
        safe_ceiling
        - inner_bound
        - arb(base.T_SPLIT) * near_bound
    ) / arb(1 - base.T_SPLIT)
    far_limit = residual_far_ceiling.lower()
    if far_limit <= 0:
        raise RuntimeError("certified inner and near bounds exhaust the pole budget")

    junction = compact.compact_cover_region(
        "outer junction signed chart 1/128<=t<=1/64",
        signed_junction_density,
        base.T_SPLIT, JUNCTION_SPLIT,
        fmpq(0), base.U_MAX,
        far_limit,
        args.max_split_depth,
        args.max_boxes,
        args.max_terminal_records,
    )
    middle = compact.compact_cover_region(
        "outer middle factored cosine chart 1/64<=t<=255/256",
        endpoint.factored_outer_far_density,
        JUNCTION_SPLIT, ENDPOINT_SPLIT,
        fmpq(0), base.U_MAX,
        far_limit,
        args.max_split_depth,
        args.max_boxes,
        args.max_terminal_records,
    )
    endpoint_region = compact.compact_cover_region(
        "outer endpoint signed chart 255/256<=t<=1",
        endpoint.signed_outer_far_density,
        ENDPOINT_SPLIT, fmpq(1),
        fmpq(0), base.U_MAX,
        far_limit,
        args.max_split_depth,
        args.max_boxes,
        args.max_terminal_records,
    )

    junction_bound = arb(junction["worst_absolute_upper"])
    middle_bound = arb(middle["worst_absolute_upper"])
    endpoint_bound = arb(endpoint_region["worst_absolute_upper"])
    total_bound = (
        inner_bound
        + arb(base.T_SPLIT) * near_bound
        + arb(JUNCTION_SPLIT - base.T_SPLIT) * junction_bound
        + arb(ENDPOINT_SPLIT - JUNCTION_SPLIT) * middle_bound
        + arb(1 - ENDPOINT_SPLIT) * endpoint_bound
    )
    epsilon = fmpq(1, 2**base.DYADIC_EXPONENT)
    modulus_loss = total_bound * arb(epsilon)
    boundary_floor = arb(base.PHI_FLOOR)
    remaining_margin = safe_ceiling - total_bound

    split_checks = {
        "junction begins at outer-near endpoint": base.T_SPLIT == fmpq(1, 128),
        "junction and middle meet exactly": JUNCTION_SPLIT == fmpq(1, 64),
        "middle and endpoint meet exactly": ENDPOINT_SPLIT == fmpq(255, 256),
        "outer slab widths sum to one minus t_split": (
            (JUNCTION_SPLIT - base.T_SPLIT)
            + (ENDPOINT_SPLIT - JUNCTION_SPLIT)
            + (fmpq(1) - ENDPOINT_SPLIT)
            == fmpq(1) - base.T_SPLIT
        ),
        "signed junction N has positive exact lower bound": (
            2 * base.T_SPLIT**2 * (1 - base.U_MAX) > 0
        ),
    }
    conditions = {
        "endpoint and factorization exact audit passed": audit["status"] == "PASSED",
        "exact outer split audit passed": all(split_checks.values()),
        "inner compact density cover certified": inner["status"] == "CERTIFIED",
        "outer near compact density cover certified": near["status"] == "CERTIFIED",
        "outer junction signed cover certified": junction["status"] == "CERTIFIED",
        "outer middle factored cover certified": middle["status"] == "CERTIFIED",
        "outer endpoint signed cover certified": endpoint_region["status"] == "CERTIFIED",
        "all outer slab bounds are below residual uniform ceiling": all(
            bound < far_limit
            for bound in (junction_bound, middle_bound, endpoint_bound)
        ),
        "weighted total is below absolute safe ceiling": total_bound < safe_ceiling,
        "modulus loss at 2^-24 is below boundary floor": modulus_loss < boundary_floor,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"

    script_path = Path(__file__)
    endpoint_path = script_path.with_name("prolate_axis_pole_modulus_endpoint_arb.py")
    dfs_path = script_path.with_name("prolate_axis_pole_modulus_compact_dfs_arb.py")
    base_path = script_path.with_name("prolate_axis_pole_modulus_arb.py")
    result = {
        "status": status,
        "certified_statement": (
            "|A_lambda''(w)| is bounded by the certified four-region weighted "
            "constant on 1<=lambda<=100, 63/64<=w<1, and "
            "Psi_lambda(1-u)>0 for 0<=u<=2^-24 using Phi>=43/5000."
            if status == "CERTIFIED" else None
        ),
        "constants": {
            "u_max": str(base.U_MAX),
            "t_split": str(base.T_SPLIT),
            "junction_split": str(JUNCTION_SPLIT),
            "endpoint_split": str(ENDPOINT_SPLIT),
            "boundary_floor": str(base.PHI_FLOOR),
            "dyadic_exponent": base.DYADIC_EXPONENT,
            "epsilon": str(epsilon),
            "absolute_safe_ceiling": str(safe_ceiling_q),
            "inner_density_bound": str(inner_bound),
            "outer_near_density_bound": str(near_bound),
            "residual_uniform_far_ceiling": str(residual_far_ceiling),
            "outer_far_acceptance_limit": str(far_limit),
            "junction_density_bound": str(junction_bound),
            "middle_density_bound": str(middle_bound),
            "endpoint_density_bound": str(endpoint_bound),
            "A_second_weighted_bound": str(total_bound),
            "remaining_absolute_margin": str(remaining_margin),
            "modulus_loss": str(modulus_loss),
        },
        "conditions": conditions,
        "split_exact_audit": {
            "status": "PASSED" if all(split_checks.values()) else "FAILED",
            "checks": split_checks,
        },
        "regions": [inner, near, junction, middle, endpoint_region],
        "endpoint_factorization_audit": audit,
        "environment": {
            "python": platform.python_version(),
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
        "arithmetic": "python-flint Arb/Acb ball arithmetic plus exact rational coverage",
        "derivation": (
            "The outer integral is partitioned exactly at t=1/128, 1/64, "
            "255/256, and 1. Each slab is bounded independently and multiplied "
            "by its exact rational t width before the final modulus comparison."
        ),
        "script_sha256": sha256_file(script_path),
        "endpoint_wrapper_sha256": sha256_file(endpoint_path),
        "dfs_driver_sha256": sha256_file(dfs_path),
        "base_driver_sha256": sha256_file(base_path),
    }

    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{result['script_sha256']}  {script_path.name}\n"
        f"{result['endpoint_wrapper_sha256']}  {endpoint_path.name}\n"
        f"{result['dfs_driver_sha256']}  {dfs_path.name}\n"
        f"{result['base_driver_sha256']}  {base_path.name}\n"
        f"{sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if status == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
