#!/usr/bin/env python3
"""Combine exact adjacent item 6 center-convexity block certificates."""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
from fractions import Fraction
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_interval(text: str) -> tuple[Fraction, Fraction]:
    if not (text.startswith("[") and text.endswith("]")):
        raise ValueError(f"invalid interval: {text}")
    lo, hi = text[1:-1].split(",", 1)
    return Fraction(lo), Fraction(hi)


def lower_fraction(real_lower: str) -> Fraction:
    token = real_lower.split()[0].strip("[")
    return Fraction(token)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pattern", default="prolate_axis_center_cap_block_*.json")
    parser.add_argument("--v-lo", default="0")
    parser.add_argument("--v-hi", default="1/20")
    parser.add_argument("--lambda-lo", default="1")
    parser.add_argument("--lambda-hi", default="10")
    parser.add_argument("--expected-blocks", type=int, default=9)
    parser.add_argument("--json", default="prolate_axis_center_cap_combined.json")
    args = parser.parse_args()

    paths = sorted(Path(path) for path in glob.glob(args.pattern))
    if not paths:
        raise FileNotFoundError(f"no block certificates match {args.pattern}")

    target_v = (Fraction(args.v_lo), Fraction(args.v_hi))
    target_lambda = (Fraction(args.lambda_lo), Fraction(args.lambda_hi))

    blocks = []
    for path in paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        lam_lo, lam_hi = parse_interval(data["target_rectangle"]["lambda"])
        v_lo, v_hi = parse_interval(data["target_rectangle"]["v"])
        blocks.append((lam_lo, lam_hi, v_lo, v_hi, path, data))
    blocks.sort(key=lambda item: item[0])

    common_requested_v = all(
        (block[2], block[3]) == target_v for block in blocks
    )
    adjacent = all(
        blocks[index][1] == blocks[index + 1][0]
        for index in range(len(blocks) - 1)
    )
    requested_endpoints = bool(
        blocks
        and blocks[0][0] == target_lambda[0]
        and blocks[-1][1] == target_lambda[1]
    )
    nonoverlapping = bool(
        blocks
        and all(block[0] < block[1] for block in blocks)
        and all(
            blocks[index][1] <= blocks[index + 1][0]
            for index in range(len(blocks) - 1)
        )
    )
    expected_count = len(blocks) == args.expected_blocks
    all_certified = all(block[5]["status"] == "CERTIFIED" for block in blocks)
    zero_terminal = all(
        block[5]["counts"]["terminal_boxes"] == 0 for block in blocks
    )
    exact_block_coverage = all(
        block[5]["conditions"].get("exact rational coverage of block", False)
        for block in blocks
    )

    worst_candidates = [
        (
            lower_fraction(
                block[5]["worst_certified_leaf"]["A_second"]["real_lower"]
            ),
            block,
        )
        for block in blocks
        if block[5].get("worst_certified_leaf")
    ]
    worst_block = (
        min(worst_candidates, key=lambda item: item[0])[1]
        if worst_candidates
        else None
    )

    conditions = {
        "expected number of lambda blocks present": expected_count,
        "all blocks certified": all_certified,
        "all blocks use the requested v interval": common_requested_v,
        "combined lambda endpoints equal the requested interval": requested_endpoints,
        "lambda blocks are exactly adjacent": adjacent,
        "lambda block interiors do not overlap": nonoverlapping,
        "every block has exact rational coverage": exact_block_coverage,
        "zero terminal boxes across all blocks": zero_terminal,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"

    result = {
        "status": status,
        "scope": (
            f"A_lambda''(v)>0 for {target_v[0]}<=v<={target_v[1]}, "
            f"{target_lambda[0]}<=lambda<={target_lambda[1]}"
        ),
        "derived_conclusion": (
            f"If status is CERTIFIED, Psi_lambda(w)>0 for "
            f"{target_lambda[0]}<=lambda<={target_lambda[1]} and "
            f"0<w<={target_v[1]}."
        ),
        "target_rectangle": {
            "v": f"[{target_v[0]},{target_v[1]}]",
            "lambda": f"[{target_lambda[0]},{target_lambda[1]}]",
        },
        "conditions": conditions,
        "counts": {
            "blocks": len(blocks),
            "evaluations": sum(block[5]["counts"]["evaluations"] for block in blocks),
            "certified_leaves": sum(
                block[5]["counts"]["certified_leaves"] for block in blocks
            ),
            "terminal_boxes": sum(
                block[5]["counts"]["terminal_boxes"] for block in blocks
            ),
        },
        "block_files": [
            {
                "path": block[4].name,
                "sha256": sha256_file(block[4]),
                "v": f"[{block[2]},{block[3]}]",
                "lambda": f"[{block[0]},{block[1]}]",
                "status": block[5]["status"],
                "counts": block[5]["counts"],
            }
            for block in blocks
        ],
        "worst_certified_leaf": (
            worst_block[5]["worst_certified_leaf"] if worst_block else None
        ),
        "combiner_sha256": sha256_file(Path(__file__)),
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
