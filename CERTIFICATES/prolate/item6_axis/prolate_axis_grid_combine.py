#!/usr/bin/env python3
"""Strict exact combiner for a two-dimensional grid of block certificates."""
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


def rectangles_overlap_interior(
    a: tuple[Fraction, Fraction, Fraction, Fraction],
    b: tuple[Fraction, Fraction, Fraction, Fraction],
) -> bool:
    ax0, ax1, ay0, ay1 = a
    bx0, bx1, by0, by1 = b
    return max(ax0, bx0) < min(ax1, bx1) and max(ay0, by0) < min(ay1, by1)


def endpoint_fraction(text: str) -> Fraction:
    return Fraction(text.split()[0].strip("["))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pattern", required=True)
    parser.add_argument("--x-key", default="w")
    parser.add_argument("--y-key", default="lambda")
    parser.add_argument("--value-key", required=True)
    parser.add_argument("--sign", choices=("positive", "negative"), required=True)
    parser.add_argument("--x-lo", required=True)
    parser.add_argument("--x-hi", required=True)
    parser.add_argument("--y-lo", required=True)
    parser.add_argument("--y-hi", required=True)
    parser.add_argument("--expected-blocks", type=int, required=True)
    parser.add_argument("--json", required=True)
    args = parser.parse_args()

    paths = sorted(Path(p) for p in glob.glob(args.pattern))
    target = (
        Fraction(args.x_lo),
        Fraction(args.x_hi),
        Fraction(args.y_lo),
        Fraction(args.y_hi),
    )
    target_area = (target[1] - target[0]) * (target[3] - target[2])

    blocks = []
    all_certified = True
    exact_subcoverage = True
    zero_terminal = True
    correct_sign = True
    all_inside = True
    total_evaluations = 0
    total_leaves = 0
    total_terminal = 0
    extremal_value = None
    extremal_record = None

    for path in paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        x0, x1 = parse_interval(data["target_rectangle"][args.x_key])
        y0, y1 = parse_interval(data["target_rectangle"][args.y_key])
        rect = (x0, x1, y0, y1)
        blocks.append((rect, path, data))

        all_certified = all_certified and data.get("status") == "CERTIFIED"
        exact_subcoverage = exact_subcoverage and bool(
            data.get("conditions", {}).get("exact rational coverage of block", False)
        )
        counts = data.get("counts", {})
        terminal_count = int(
            counts.get("terminal_boxes", counts.get("terminal_intervals", 0))
        )
        zero_terminal = zero_terminal and terminal_count == 0
        total_evaluations += int(counts.get("evaluations", 0))
        total_leaves += int(counts.get("certified_leaves", counts.get("leaves", 0)))
        total_terminal += terminal_count

        declared_sign = data.get("certified_sign")
        if declared_sign is not None:
            correct_sign = correct_sign and declared_sign == args.sign

        all_inside = all_inside and (
            target[0] <= x0 < x1 <= target[1]
            and target[2] <= y0 < y1 <= target[3]
        )

        leaf = (
            data.get("extremal_certified_leaf")
            or data.get("worst_certified_leaf")
            or data.get("least_negative_certified_leaf")
        )
        if leaf and args.value_key in leaf:
            record = leaf[args.value_key]
            value = (
                endpoint_fraction(record["real_lower"])
                if args.sign == "positive"
                else endpoint_fraction(record["real_upper"])
            )
            if (
                extremal_value is None
                or (args.sign == "positive" and value < extremal_value)
                or (args.sign == "negative" and value > extremal_value)
            ):
                extremal_value = value
                extremal_record = {"source": str(path), "leaf": leaf}

    rectangles = [item[0] for item in blocks]
    no_overlap = all(
        not rectangles_overlap_interior(rectangles[i], rectangles[j])
        for i in range(len(rectangles))
        for j in range(i + 1, len(rectangles))
    )
    total_area = sum(
        ((rect[1] - rect[0]) * (rect[3] - rect[2]) for rect in rectangles),
        Fraction(0),
    )
    exact_union = all_inside and no_overlap and total_area == target_area
    expected_count = len(paths) == args.expected_blocks

    conditions = {
        "expected number of block files": expected_count,
        "all block certificates are CERTIFIED": all_certified,
        "all blocks have exact rational subcoverage": exact_subcoverage,
        "all blocks have zero terminal boxes": zero_terminal,
        "all declared signs match requested sign": correct_sign,
        "all block rectangles lie inside the target": all_inside,
        "block rectangle interiors are pairwise disjoint": no_overlap,
        "block rectangle areas exactly cover the target": exact_union,
    }
    status = "CERTIFIED" if all(conditions.values()) else "INCOMPLETE"
    comparison = ">0" if args.sign == "positive" else "<0"

    result = {
        "status": status,
        "certified_statement": (
            f"{args.value_key}{comparison} on "
            f"{args.x_lo}<={args.x_key}<={args.x_hi}, "
            f"{args.y_lo}<={args.y_key}<={args.y_hi}"
            if status == "CERTIFIED"
            else None
        ),
        "requested_sign": args.sign,
        "target_rectangle": {
            args.x_key: f"[{args.x_lo},{args.x_hi}]",
            args.y_key: f"[{args.y_lo},{args.y_hi}]",
            "exact_area": str(target_area),
        },
        "conditions": conditions,
        "counts": {
            "block_files": len(paths),
            "evaluations": total_evaluations,
            "certified_leaves": total_leaves,
            "terminal_boxes": total_terminal,
        },
        "extremal_certified_leaf": extremal_record,
        "blocks": [
            {
                "path": str(path),
                "sha256": sha256_file(path),
                args.x_key: f"[{rect[0]},{rect[1]}]",
                args.y_key: f"[{rect[2]},{rect[3]}]",
                "status": data.get("status"),
                "counts": data.get("counts", {}),
                "script_sha256": data.get("script_sha256"),
            }
            for rect, path, data in blocks
        ],
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
