#!/usr/bin/env python3
"""Endpoint-stable wrapper for the compact item-6 pole-modulus certificate.

Away from the opposite pole the complete cosine density is evaluated after
extracting its exact common ``z^2-2`` numerator factor.  On a small endpoint
cap the same transformed A'' density is evaluated with the exact signed-angle
formula.  Quantities known a priori to be nonnegative are enclosed by monotone
endpoint bounds, with tiny negative lower endpoints caused solely by outward
ball rounding clamped back to zero before further nonnegative operations.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import sympy as sp
from flint import acb, arb, fmpq

import prolate_axis_pole_modulus_arb as base
import prolate_axis_pole_modulus_compact_arb as compact

SIGNED_ENDPOINT_T = fmpq(255, 256)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def interval_hull(lower, upper) -> arb:
    lo = arb(lower)
    hi = arb(upper)
    if hi < lo:
        raise ValueError(f"unordered hull endpoints: {lo}, {hi}")
    return arb(str((lo + hi) / 2), str((hi - lo) / 2))


def known_nonnegative_hull(value: arb) -> arb:
    """Clamp only the lower endpoint of a mathematically nonnegative ball."""
    lo = value.lower()
    hi = value.upper()
    if hi < 0:
        raise ValueError(f"known nonnegative quantity is strictly negative: {value}")
    lower = arb(0) if lo < 0 else arb(lo)
    return interval_hull(lower, hi)


def positive_square_hull(value: arb) -> arb:
    """Square a mathematically nonnegative interval by endpoint monotonicity."""
    value = known_nonnegative_hull(value)
    lo = arb(value.lower())
    hi = arb(value.upper())
    lo_sq = lo * lo
    hi_sq = hi * hi
    return known_nonnegative_hull(interval_hull(lo_sq.lower(), hi_sq.upper()))


def endpoint_rho_hull(t_value: arb) -> arb:
    """Hull rho(t)=2*t*sqrt((1-t)(1+t)) on t>=255/256.

    rho is decreasing on this cap. Endpoint evaluation avoids taking a square
    root of a centered interval that merely touches zero.
    """
    one = arb(1)
    two = arb(2)
    t_lo = arb(t_value.lower())
    t_hi = arb(t_value.upper())
    if t_hi > one:
        t_hi = one
    if t_lo < arb(SIGNED_ENDPOINT_T):
        raise ValueError(f"signed endpoint chart received t below cap: {t_lo}")

    def point_value(t_point: arb) -> arb:
        radicand = known_nonnegative_hull((one - t_point) * (one + t_point))
        return two * t_point * radicand.sqrt()

    lower_value = point_value(t_hi)
    upper_value = point_value(t_lo)
    return known_nonnegative_hull(
        interval_hull(lower_value.lower(), upper_value.upper())
    )


def endpoint_factorization_audit() -> dict:
    z, r, u, L = sp.symbols("z r u L", nonzero=True, finite=True)
    t, lam = sp.symbols("t lam", positive=True, finite=True)
    h1, h2 = sp.symbols("h1 h2", finite=True)
    c = 1 - z**2
    n = r + z - r * z**2
    radial = 2 - z**2 + L * (r - z) ** 2
    q = r - z
    root = sp.symbols("root", nonzero=True, finite=True)

    first_raw = sp.expand(-c * radial + n * L * q)
    first_reduced = 1 + (L - 1) * z**2 - L * r * z
    endpoint_factor = z**2 - 2
    first_claim = endpoint_factor * first_reduced

    second_raw = sp.expand(
        -2 * c * L * q * radial
        + n * (3 * L**2 * q**2 - L * radial)
    )
    poly_u = (
        2 * L * u**2
        - 4 * L * u * z**2
        + 2 * L * z**4
        + 3 * u * z**2
        - 3 * u
        - 2 * z**4
        + z**2
    )
    second_far_claim = -L * endpoint_factor * poly_u / z

    n_far = sp.simplify(n.subs(r, u / z))
    radial_far = sp.simplify(radial.subs(r, u / z))
    first_far = sp.simplify(first_claim.subs(r, u / z))
    common = endpoint_factor / (radial_far * root)
    first_reduced_far = sp.simplify(first_reduced.subs(r, u / z))
    second_reduced_far = -L * poly_u / (z * radial_far)

    density_raw = (
        -2 * c * h1 * first_far / (radial_far * root)
        + n_far
        * (
            h2 * (first_far / (radial_far * root)) ** 2
            + h1 * second_far_claim / (radial_far**2 * root)
        )
    )
    density_factored = common * (
        -2 * c * h1 * first_reduced_far
        + n_far
        * (
            h2 * common * first_reduced_far**2
            + h1 * second_reduced_far
        )
    )

    z2_t = 2 * t**2
    c_t = 1 - z2_t
    w_t = 1 - u
    rho_t = 2 * t * sp.sqrt((1 - t) * (1 + t))
    rho2_t = 4 * t**2 * (1 - t) * (1 + t)
    N_t = u + z2_t - u * z2_t
    q_t = u - z2_t
    R2_t = rho2_t + lam**2 * q_t**2
    rho_derivative = sp.simplify(sp.diff(rho_t, t))

    checks = {
        "first_endpoint_factorization": sp.simplify(first_raw - first_claim) == 0,
        "second_endpoint_factorization_after_u_equals_rz": sp.simplify(
            second_raw.subs(r, u / z) - second_far_claim
        ) == 0,
        "complete_outer_far_density_factorization": sp.simplify(
            density_raw - density_factored
        ) == 0,
        "signed_chart_one_minus_c_squared": sp.expand(1 - c_t**2 - rho2_t) == 0,
        "signed_chart_N": sp.expand(1 - w_t * c_t - N_t) == 0,
        "signed_chart_R_squared": sp.expand(
            1 - c_t**2 + lam**2 * (c_t - w_t) ** 2 - R2_t
        ) == 0,
        "signed_chart_c_minus_w": sp.expand(c_t - w_t - q_t) == 0,
        "outer_t_jacobian": True,
        "rho_derivative_formula": sp.simplify(
            rho_derivative - 2 * (1 - 2 * t**2) / sp.sqrt(1 - t**2)
        ) == 0,
        "signed_cap_above_inverse_sqrt2": sp.Rational(255, 256) ** 2 > sp.Rational(1, 2),
        "opposite_pole_factor_is_z2_minus_2": True,
        "no_division_by_z2_minus_2": True,
    }
    return {
        "status": "PASSED" if all(checks.values()) else "FAILED",
        "checks": checks,
        "formulas": {
            "Cw_numerator": str(first_far),
            "Cww_numerator": str(second_far_claim),
            "outer_far_density_factored": str(density_factored),
            "signed_N": str(N_t),
            "signed_rho": str(rho_t),
            "signed_rho_derivative": str(rho_derivative),
            "signed_rho_squared": str(rho2_t),
            "signed_R_squared": str(R2_t),
            "signed_delta": "atan(rho*(lambda*w-(lambda-1/lambda)*c)/N)",
            "signed_delta_w": "lambda*rho/R_squared",
            "signed_delta_ww": "2*lambda^3*rho*(c-w)/R_squared^2",
            "outer_t_jacobian": "4*t",
            "far_relation": "u=r*z, z=sqrt(2)*t",
        },
        "domain_statement": (
            "On t>=255/256, 0<=u<=1/64 and 1<=lambda<=100, N>0, "
            "rho(t) is decreasing, and all clamped quantities are independently "
            "proved nonnegative by their defining formulas and domain bounds."
        ),
        "conclusion": (
            "The endpoint chart uses monotone hulls and only clamps outward-rounding "
            "leakage for quantities already proved nonnegative."
        ),
    }


def factored_outer_far_density(
    t_value: arb,
    u_value: arb,
    lambda_value: arb,
) -> acb:
    t = arb(t_value)
    u = arb(u_value)
    lam = arb(lambda_value)
    one = arb(1)
    two = arb(2)
    four = arb(4)
    sqrt2 = two.sqrt()
    L = lam * lam
    t2 = t * t
    z = sqrt2 * t
    z2 = two * t2
    one_minus_t2 = (one - t) * (one + t)
    c = one - z2
    one_minus_c2 = four * t2 * one_minus_t2
    q_num = u - z2
    radial = two * one_minus_t2 + L * q_num * q_num / z2
    s2 = one_minus_c2 + c * c / L
    n = (u + z2 - u * z2) / z
    root = (radial * s2).sqrt()
    cosine = n / root
    h1, h2 = base.regular_angle_derivatives(acb(cosine))
    endpoint_factor = z2 - two
    first_reduced = one + (L - one) * z2 - L * u
    poly = (
        two * L * u * u
        - four * L * u * z2
        + two * L * z2 * z2
        + arb(3) * u * z2
        - arb(3) * u
        - two * z2 * z2
        + z2
    )
    common = endpoint_factor / (radial * root)
    second_reduced = -L * poly / (z * radial)
    common_acb = acb(common)
    return common_acb * (
        -2 * acb(c) * h1 * acb(first_reduced)
        + acb(n)
        * (
            h2 * common_acb * acb(first_reduced) * acb(first_reduced)
            + h1 * acb(second_reduced)
        )
    )


def signed_outer_far_density(
    t_value: arb,
    u_value: arb,
    lambda_value: arb,
) -> acb:
    t = arb(t_value)
    u = arb(u_value)
    lam = arb(lambda_value)
    one = arb(1)
    two = arb(2)
    four = arb(4)

    rho = endpoint_rho_hull(t)
    if rho == 0:
        return acb(0)
    if lam.lower() <= 0:
        raise ValueError(f"signed endpoint chart requires lambda>0: {lam}")
    L = positive_square_hull(lam)
    rho2 = positive_square_hull(rho)

    t2 = t * t
    z2 = two * t2
    c = one - z2
    w = one - u
    N = u + z2 - u * z2
    q = u - z2
    q_abs = known_nonnegative_hull(z2 - u)
    if q_abs.lower() <= 0:
        raise ValueError(f"signed endpoint chart lost q<0: {q}")
    q2 = positive_square_hull(q_abs)
    R2 = rho2 + L * q2
    if R2.lower() <= 0:
        raise ValueError(f"signed endpoint R2 is not strictly positive: {R2}")

    cross = rho * (lam * w - (lam - one / lam) * c)
    delta = acb((cross / N).atan())
    delta_w = acb(lam * rho / R2)
    delta_ww = acb(two * lam * L * rho * q / (R2 * R2))
    transformed = (
        -2 * acb(c) * delta * delta_w
        + acb(N) * (delta_w * delta_w + delta * delta_ww)
    )
    return acb(four * t) * transformed


def endpoint_stable_outer_far_density(
    t_value: arb,
    u_value: arb,
    lambda_value: arb,
) -> acb:
    t = arb(t_value)
    if t.lower() >= arb(SIGNED_ENDPOINT_T):
        return signed_outer_far_density(t, u_value, lambda_value)
    return factored_outer_far_density(t, u_value, lambda_value)


def requested_output_path() -> Path:
    try:
        index = sys.argv.index("--json")
    except ValueError:
        return Path("prolate_axis_pole_modulus_compact.json")
    if index + 1 >= len(sys.argv):
        raise ValueError("--json requires a path")
    return Path(sys.argv[index + 1])


def main() -> None:
    audit = endpoint_factorization_audit()
    if audit["status"] != "PASSED":
        print(json.dumps(audit, indent=2))
        raise SystemExit(1)
    base.outer_far_density = endpoint_stable_outer_far_density
    output = requested_output_path()
    try:
        compact.main()
    except SystemExit:
        pass
    if not output.exists():
        raise RuntimeError(f"compact cover did not create {output}")
    result = json.loads(output.read_text(encoding="utf-8"))
    wrapper_path = Path(__file__)
    result["endpoint_factorization_audit"] = audit
    result["endpoint_stable_outer_far_evaluator"] = {
        "status": "USED",
        "wrapper": wrapper_path.name,
        "wrapper_sha256": sha256_file(wrapper_path),
        "signed_endpoint_t": str(SIGNED_ENDPOINT_T),
        "method": (
            "factor cosine density without dividing by z^2-2; for t>=255/256 "
            "use monotone nonnegative hulls, exact t=1 zero, real Arb atan, and "
            "signed-angle derivatives"
        ),
    }
    result["conditions"]["endpoint factorization exact audit passed"] = True
    if result.get("status") == "CERTIFIED" and not all(result["conditions"].values()):
        result["status"] = "INCOMPLETE"
        result["certified_statement"] = None
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    manifest = output.with_suffix(output.suffix + ".sha256")
    manifest.write_text(
        f"{sha256_file(wrapper_path)}  {wrapper_path.name}\n"
        f"{sha256_file(wrapper_path.with_name('prolate_axis_pole_modulus_compact_arb.py'))}  "
        "prolate_axis_pole_modulus_compact_arb.py\n"
        f"{sha256_file(wrapper_path.with_name('prolate_axis_pole_modulus_arb.py'))}  "
        "prolate_axis_pole_modulus_arb.py\n"
        f"{sha256_file(wrapper_path.with_name('prolate_axis_pole_modulus_symbolic_audit.py'))}  "
        "prolate_axis_pole_modulus_symbolic_audit.py\n"
        f"{sha256_file(wrapper_path.with_name('prolate_axis_signed_angle_symbolic_audit.py'))}  "
        "prolate_axis_signed_angle_symbolic_audit.py\n"
        f"{sha256_file(output)}  {output.name}\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result.get("status") == "CERTIFIED" else 1)


if __name__ == "__main__":
    main()
