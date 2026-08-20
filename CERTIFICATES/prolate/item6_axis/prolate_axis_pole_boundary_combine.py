#!/usr/bin/env python3
"""Combine exact adjacent pole-boundary lambda-block certificates."""
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


def lower_fraction(text: str) -> Fraction:
    token = text.split()[0].strip("[")
    return Fraction(token)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pattern", default="prolate_axis_pole_boundary_block_*.json"
    )
    parser.add_argument(
        "--json", default="prolate_axis_pole_boundary_combined.json"
    )
    parser.add_argument("--expected-lo", type=int, default=1)
    parser.add_argument("--expected-hi", type=int, default=100)
    args = parser.parse_args()

    paths = sorted(Path(p) for p in glob.glob(args.pattern))
    if not paths:
        raise FileNotFoundError(f"no block certificates match {args.pattern}")

    blocks = []
    for path in paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        lo, hi = parse_interval(data["target_interval"])
        blocks.append((lo, hi, path, data))
    blocks.sort(key=lambda item: item[0])

    adjacent = all(
        blocks[index][1] == blocks[index + 1][0]
        for index in range(len(blocks) - 1)
    )
    endpoints = bool(
        blocks
        and blocks[0][0] == Fraction(args.expected_lo)
        and blocks[-1][1] == Fraction(args.expected_hi)
    )
    all_certified = all(block[3]["status"] == "CERTIFIED" for block in blocks)
    exact_block_coverage = all(
        block[3]["conditions"].get("exact rational coverage of block", False)
        for block in blocks
    )
    zero_terminal = all(
        block[3]["counts"]["terminal_intervals"] == 0 for block in blocks
    )

    worst_candidates = [
        (
            lower_fraction(block[3]["worst_certified_leaf"]["Phi"]["real_lower"]),
            block,
        )
        for block in blocks
        if block[3].get("worst_certified_leaf")
    ]
    worst_block = (
        min(worst_candidates, key=lambda item: item[0])[1]
        if worst_candidates
        else None
    )

    conditions = {
        "all blocks certified": all_certified,
        "lambda blocks are exactly adjacent": adjacent,
        "combined endpoints equal requested interval": endpoints,
        "every block has exact rational coverage": exact_block_coverage,
        "zero terminal intervals across all blocks": zero_terminal,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"

    result = {
        "status": status,
        "certified_statement": (
            f"Phi(lambda)>0 for {args.expected_lo}<=lambda<={args.expected_hi}"
            if status == "CERTIFIED"
            else None
        ),
        "definition": "Phi(lambda)=lim_{w->1^-} Psi_lambda(w)",
        "conditions": conditions,
        "counts": {
            "blocks": len(blocks),
            "evaluations": sum(
                block[3]["counts"]["evaluations"] for block in blocks
            ),
            "certified_leaves": sum(
                block[3]["counts"]["certified_leaves"] for block in blocks
            ),
            "terminal_intervals": sum(
                block[3]["counts"]["terminal_intervals"] for block in blocks
            ),
        },
        "block_files": [
            {
                "path": block[2].name,
                "sha256": sha256_file(block[2]),
                "lambda": f"[{block[0]},{block[1]}]",
                "status": block[3]["status"],
                "counts": block[3]["counts"],
            }
            for block in blocks
        ],
        "worst_certified_leaf": (
            worst_block[3]["worst_certified_leaf"] if worst_block else None
        ),
        "limitations": (
            "This certifies the boundary anchor only. A finite-width pole-cap "
            "transfer remains necessary for the full item 6 theorem."
        ),
    }
    result["combiner_sha256"] = sha256_file(Path(__file__))

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
