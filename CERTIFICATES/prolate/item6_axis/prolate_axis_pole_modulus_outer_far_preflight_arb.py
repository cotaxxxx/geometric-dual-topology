#!/usr/bin/env python3
"""Outer-far viability preflight for the item-6 pole-modulus certificate.

This script runs only the unresolved outer-far chart before the complete
three-region modulus cover.  Its requested density ceiling is derived exactly
from the final pole floor and conservative regional limits for the already
certified inner and outer-near charts:

    inner_limit + t_split*near_limit + (1-t_split)*far_limit
        = (43/5000)*2^24.

A CERTIFIED result therefore proves that the current dyadic scale remains
viable even if the inner and near regions attain their requested limits.  An
INCOMPLETE result is diagnostic only: it distinguishes capacity/dependency
failure from a certified violation and preserves the exact partition volume.
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
    parser.add_argument("--max-boxes", type=int, default=500_000)
    parser.add_argument("--inner-density-limit", default="100")
    parser.add_argument("--outer-near-density-limit", default="144100")
    parser.add_argument("--max-terminal-records", type=int, default=100)
    parser.add_argument(
        "--json",
        default="prolate_axis_pole_modulus_outer_far_preflight.json",
    )
    args = parser.parse_args()

    ctx.dps = args.dps
    inner_limit_q = fmpq(args.inner_density_limit)
    near_limit_q = fmpq(args.outer_near_density_limit)
    safe_ceiling_q = base.PHI_FLOOR * 2**base.DYADIC_EXPONENT
    far_limit_q = (
        safe_ceiling_q
        - inner_limit_q
        - base.T_SPLIT * near_limit_q
    ) / (1 - base.T_SPLIT)
    if far_limit_q <= 0:
        raise ValueError("conservative regional limits leave no positive far ceiling")

    audit = endpoint.endpoint_factorization_audit()
    if audit["status"] != "PASSED":
        print(json.dumps(audit, indent=2))
        raise SystemExit(1)

    far = compact.compact_cover_region(
        "outer far viability preflight t>=1/128",
        endpoint.endpoint_stable_outer_far_density,
        base.T_SPLIT,
        fmpq(1),
        fmpq(0),
        base.U_MAX,
        arb(far_limit_q),
        args.max_split_depth,
        args.max_boxes,
        args.max_terminal_records,
    )

    far_bound = arb(far["worst_absolute_upper"])
    conservative_total = (
        arb(inner_limit_q)
        + arb(base.T_SPLIT) * arb(near_limit_q)
        + arb(1 - base.T_SPLIT) * far_bound
    )
    safe_ceiling = arb(safe_ceiling_q)
    remaining_margin = safe_ceiling - conservative_total

    conditions = {
        "endpoint factorization exact audit passed": audit["status"] == "PASSED",
        "outer far exact compact coverage certified": far["status"] == "CERTIFIED",
        "outer far bound is below conservative far ceiling": far_bound < arb(far_limit_q),
        "conservative total is below absolute safe ceiling": conservative_total < safe_ceiling,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"

    script_path = Path(__file__)
    endpoint_path = script_path.with_name("prolate_axis_pole_modulus_endpoint_arb.py")
    compact_path = script_path.with_name("prolate_axis_pole_modulus_compact_arb.py")
    base_path = script_path.with_name("prolate_axis_pole_modulus_arb.py")
    result = {
        "status": status,
        "certified_statement": (
            "The unresolved outer-far chart fits inside the pole-modulus budget "
            "under conservative inner and outer-near regional limits."
            if status == "CERTIFIED"
            else None
        ),
        "purpose": (
            "Early viability test before the complete endpoint-stable modulus run; "
            "INCOMPLETE does not certify a counterexample."
        ),
        "constants": {
            "boundary_floor": str(base.PHI_FLOOR),
            "dyadic_exponent": base.DYADIC_EXPONENT,
            "absolute_safe_ceiling": str(safe_ceiling_q),
            "inner_conservative_limit": str(inner_limit_q),
            "outer_near_conservative_limit": str(near_limit_q),
            "derived_outer_far_ceiling": str(far_limit_q),
            "outer_far_observed_bound": str(far_bound),
            "conservative_total_bound": str(conservative_total),
            "remaining_absolute_margin": str(remaining_margin),
        },
        "conditions": conditions,
        "outer_far_region": far,
        "endpoint_factorization_audit": audit,
        "environment": {
            "python": platform.python_version(),
            "python_flint": flint.__version__,
            "flint": flint.__FLINT_VERSION__,
        },
        "arithmetic": "python-flint Arb/Acb ball arithmetic plus exact rational coverage",
        "script_sha256": sha256_file(script_path),
        "endpoint_wrapper_sha256": sha256_file(endpoint_path),
        "compact_driver_sha256": sha256_file(compact_path),
        "base_driver_sha256": sha256_file(base_path),
    }

    output = Path(args.json)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    output.with_suffix(output.suffix + ".sha256").write_text(
        f"{result['script_sha256']}  {script_path.name}\n"
        f"{result['endpoint_wrapper_sha256']}  {endpoint_path.name}\n"
        f"{result['compact_driver_sha256']}  {compact_path.name}\n"
        f"{result['base_driver_sha256']}  {base_path.name}\n"
        f"{sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if status == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
