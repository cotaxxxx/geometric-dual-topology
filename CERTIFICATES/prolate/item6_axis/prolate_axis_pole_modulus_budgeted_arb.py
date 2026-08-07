#!/usr/bin/env python3
"""Budget-adaptive complete pole-modulus certificate for prolate item 6.

The stable inner and outer-near charts are certified first. Their actual
rigorous suprema, rather than coarse requested limits, determine the exact
remaining budget for the unresolved outer-far chart:

    far_ceiling = (safe_ceiling - inner_bound - t_split*near_bound)
                  / (1-t_split).

The endpoint-stable outer-far density receives the lower endpoint of this
remaining Arb budget as its strict acceptance threshold. The final certificate
checks the complete weighted total directly against ``(43/5000)*2^24`` and
requires exact rational coverage with terminal zero in all three regions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import flint
from flint import arb, ctx, fmpq

import prolate_axis_pole_modulus_arb as base
import prolate_axis_pole_modulus_compact_arb as compact
import prolate_axis_pole_modulus_endpoint_arb as endpoint


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--max-split-depth", type=int, default=24)
    parser.add_argument("--max-boxes", type=int, default=1_500_000)
    parser.add_argument("--inner-density-limit", default="100")
    parser.add_argument("--outer-near-density-limit", default="144100")
    parser.add_argument("--max-terminal-records", type=int, default=100)
    parser.add_argument("--json", default="prolate_axis_pole_modulus_budgeted.json")
    args = parser.parse_args()

    ctx.dps = args.dps
    inner_limit = arb(args.inner_density_limit)
    near_limit = arb(args.outer_near_density_limit)
    safe_ceiling_q = base.PHI_FLOOR * 2**base.DYADIC_EXPONENT
    safe_ceiling = arb(safe_ceiling_q)

    audit = endpoint.endpoint_factorization_audit()
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
    remaining_far_budget = (
        safe_ceiling
        - inner_bound
        - arb(base.T_SPLIT) * near_bound
    ) / arb(1 - base.T_SPLIT)
    far_limit = remaining_far_budget.lower()
    if far_limit <= 0:
        raise RuntimeError("certified inner and near bounds exhaust the pole budget")

    far = compact.compact_cover_region(
        "outer far t>=1/128 with dynamic residual budget",
        endpoint.endpoint_stable_outer_far_density,
        base.T_SPLIT, fmpq(1),
        fmpq(0), base.U_MAX,
        far_limit,
        args.max_split_depth,
        args.max_boxes,
        args.max_terminal_records,
    )

    far_bound = arb(far["worst_absolute_upper"])
    total_bound = (
        inner_bound
        + arb(base.T_SPLIT) * near_bound
        + arb(1 - base.T_SPLIT) * far_bound
    )
    epsilon = fmpq(1, 2**base.DYADIC_EXPONENT)
    modulus_loss = total_bound * arb(epsilon)
    boundary_floor = arb(base.PHI_FLOOR)
    remaining_margin = safe_ceiling - total_bound

    conditions = {
        "endpoint factorization exact audit passed": audit["status"] == "PASSED",
        "inner compact density cover certified": inner["status"] == "CERTIFIED",
        "outer near compact density cover certified": near["status"] == "CERTIFIED",
        "outer far dynamic-budget cover certified": far["status"] == "CERTIFIED",
        "outer far observed bound is below dynamic residual ceiling": far_bound < far_limit,
        "weighted total is below absolute safe ceiling": total_bound < safe_ceiling,
        "modulus loss at 2^-24 is below boundary floor": modulus_loss < boundary_floor,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"

    script_path = Path(__file__)
    endpoint_path = script_path.with_name("prolate_axis_pole_modulus_endpoint_arb.py")
    compact_path = script_path.with_name("prolate_axis_pole_modulus_compact_arb.py")
    base_path = script_path.with_name("prolate_axis_pole_modulus_arb.py")
    symbolic_path = script_path.with_name("prolate_axis_pole_modulus_symbolic_audit.py")
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
            "absolute_safe_ceiling": str(safe_ceiling_q),
            "inner_requested_limit": str(inner_limit),
            "outer_near_requested_limit": str(near_limit),
            "inner_density_bound": str(inner_bound),
            "outer_near_density_bound": str(near_bound),
            "derived_outer_far_ceiling": str(remaining_far_budget),
            "outer_far_acceptance_limit": str(far_limit),
            "outer_far_density_bound": str(far_bound),
            "A_second_uniform_bound": str(total_bound),
            "remaining_absolute_margin": str(remaining_margin),
            "modulus_loss": str(modulus_loss),
        },
        "conditions": conditions,
        "regions": [inner, near, far],
        "endpoint_factorization_audit": audit,
        "environment": {
            "python": platform.python_version(),
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
        "arithmetic": "python-flint Arb/Acb ball arithmetic plus exact rational coverage",
        "derivation": (
            "Inner and outer-near are certified first. Their actual rigorous "
            "suprema determine the strict residual outer-far ceiling. The final "
            "weighted total and modulus loss are checked again directly."
        ),
        "script_sha256": sha256_file(script_path),
        "endpoint_wrapper_sha256": sha256_file(endpoint_path),
        "compact_driver_sha256": sha256_file(compact_path),
        "base_driver_sha256": sha256_file(base_path),
        "symbolic_audit_sha256": sha256_file(symbolic_path),
    }

    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{result['script_sha256']}  {script_path.name}\n"
        f"{result['endpoint_wrapper_sha256']}  {endpoint_path.name}\n"
        f"{result['compact_driver_sha256']}  {compact_path.name}\n"
        f"{result['base_driver_sha256']}  {base_path.name}\n"
        f"{result['symbolic_audit_sha256']}  {symbolic_path.name}\n"
        f"{sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if status == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
