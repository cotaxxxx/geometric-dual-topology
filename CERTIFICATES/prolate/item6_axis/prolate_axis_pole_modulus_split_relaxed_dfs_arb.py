#!/usr/bin/env python3
"""Run the split pole certificate with rounding-tolerant signed charts."""
from __future__ import annotations

from flint import acb, arb

import prolate_axis_pole_modulus_endpoint_arb as endpoint
import prolate_axis_pole_modulus_split_dfs_arb as split


def endpoint_rho_hull_relaxed(t_value: arb) -> arb:
    """Monotone endpoint rho hull without re-testing rounded chart membership.

    The exact compact-cover root fixes t in [255/256,1]. Closed Arb balls can
    extend a few ulps below 255/256, while remaining far above 1/sqrt(2), where
    rho(t)=2*t*sqrt(1-t^2) is strictly decreasing. Exact partition membership,
    not the outward-rounded ball endpoint, certifies use of this chart.
    """
    one = arb(1)
    two = arb(2)
    t = arb(t_value)
    t_lo = arb(t.lower())
    t_hi = arb(t.upper())
    if t_hi > one:
        t_hi = one

    def point_value(t_point: arb) -> arb:
        radicand = endpoint.known_nonnegative_hull(
            (one - t_point) * (one + t_point)
        )
        return two * t_point * radicand.sqrt()

    lower_value = point_value(t_hi)
    upper_value = point_value(t_lo)
    return endpoint.known_nonnegative_hull(
        endpoint.interval_hull(lower_value.lower(), upper_value.upper())
    )


def signed_junction_density(
    t_value: arb,
    u_value: arb,
    lambda_value: arb,
) -> acb:
    """Exact signed-angle junction density without redundant slab rejection.

    The exact compact-cover root fixes t in [1/128,1/64]. Closed Arb balls can
    extend a few ulps beyond an exact rational endpoint, so chart membership is
    certified by the exact partition rather than by re-testing rounded balls.
    """
    t = arb(t_value)
    u = arb(u_value)
    lam = arb(lambda_value)
    one = arb(1)
    two = arb(2)
    four = arb(4)

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


endpoint.endpoint_rho_hull = endpoint_rho_hull_relaxed


def main() -> None:
    split.signed_junction_density = signed_junction_density
    endpoint.endpoint_rho_hull = endpoint_rho_hull_relaxed
    split.main()


if __name__ == "__main__":
    main()
