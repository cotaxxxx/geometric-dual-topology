"""Exact candidate and interval helper functions."""
from calibration_context import *
from calibration_security import *

def _nearest_dyadic(value: Fraction, bits: int = 96) -> Dyadic:
    scale = 1 << bits
    numerator = value.numerator * scale
    denominator = value.denominator
    sign = -1 if numerator < 0 else 1
    quotient, remainder = divmod(abs(numerator), denominator)
    if 2 * remainder >= denominator:
        quotient += 1
    return Dyadic.canonical(sign * quotient, bits)

def _candidate_pairs(config: dict[str, Any]) -> list[tuple[Dyadic, Dyadic]]:
    widths = _dyadic_list(config["candidate_lambda_widths"], "candidate_lambda_widths")
    radii = _dyadic_list(config["candidate_tube_radii"], "candidate_tube_radii")
    return [(width, radius) for width in widths for radius in radii]

def _cell_partition(start: Fraction, end: Fraction, width: Fraction, maximum: int):
    cells = []
    left = start
    while left < end:
        right = min(left + width, end)
        if not left < right:
            raise CalibrationError("nonpositive calibration cell")
        cells.append((left, right))
        if len(cells) > maximum:
            raise CalibrationError("maximum cell budget exceeded")
        left = right
    return cells

def _rational_arb(value: Fraction, arb_type):
    return arb_type(value.numerator) / arb_type(value.denominator)

def _dyadic_arb(value: Dyadic, arb_type):
    return arb_type(value.m) / arb_type(1 << value.e)

def _fraction_box(lo: Fraction, hi: Fraction, arb_type):
    midpoint = (lo + hi) / 2
    radius = (hi - lo) / 2
    return _rational_arb(midpoint, arb_type) + _rational_arb(radius, arb_type) * arb_type("+/- 1.0")

def _dyadic_box(interval: DyadicInterval, arb_type):
    midpoint = interval.midpoint()
    radius = (interval.hi - interval.lo) * Dyadic(1, 1)
    return _dyadic_arb(midpoint, arb_type) + _dyadic_arb(radius, arb_type) * arb_type("+/- 1.0")

def _newton_predictor(kernel, arb_type, lam: Fraction, seed: Dyadic, *, iterations: int,
                      tol: str, depth: int, limit: int) -> Dyadic:
    current = seed
    lam_ball = _rational_arb(lam, arb_type)
    for _ in range(iterations):
        point = _dyadic_arb(current, arb_type)
        residual = arb_ball_to_exact_interval(
            kernel.F_arb(point, lam_ball, tol=tol, depth=depth, limit=limit)
        )
        slope = arb_ball_to_exact_interval(
            kernel.dFdr_arb(point, lam_ball, tol=tol, depth=depth, limit=limit)
        )
        slope_mid = slope.midpoint()
        if slope_mid == D_ZERO:
            break
        updated = current.as_fraction() - residual.midpoint().as_fraction() / slope_mid.as_fraction()
        current = _nearest_dyadic(updated)
    return current

def _append_record(records: list[dict[str, Any]], previous: str, body: dict[str, Any]) -> str:
    record = dict(body)
    record["previous_record_sha256"] = previous
    assert_result_namespace(record)
    raw = canonical_json_bytes(record)
    records.append(record)
    return sha256_hex(raw)

__all__ = [name for name in globals() if not name.startswith("__")]
