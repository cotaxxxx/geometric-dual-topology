#!/usr/bin/env python3
"""Rigorous Arb enclosure of the prolate boundary-entry function B(a)=F(1,a).

Cancellation-preserving N-form
------------------------------
With  s=sin(th), c=cos(th), ell=s^2+a^2 c^2, u=s cos(phi), m=1-u,

    P  := 1 - 2u/ell + 1/ell
    ell*P = (ell-1) + 2m = (a^2-1) c^2 + 2m  =:  N

The identity ell*P = N replaces a difference of O(1) quantities by a sum
of two nonnegative terms, so the enclosure of P no longer collapses near
the contact point (th,phi)=(pi/2,0).  Likewise

    u*N + m^2 = (a^2-1) c^2 u + m(2-m)  =:  K

is written without subtraction of comparable terms.  In these variables

    cos_alpha = C m sqrt(ell) / sqrt(N)
    d/dr [ (1-ru) h(cos_alpha) ] |_{r=1} * s
        = ( -u h(x) - C sqrt(ell) h'(x) m K / N^{3/2} ) * s

The 1/rho blow-up of d(cos_alpha)/dr is cancelled *symbolically* by the
cone weight m, so the integrand is manifestly bounded: since
N >= max((a^2-1)c^2, 2m),

    | m K / N^{3/2} |  <=  sqrt(2 m)

which is the corner bound used for cells that contain N=0.

B(a) = (1/2pi) int_0^{pi/2} int_0^{2pi} G(th,phi,a) dphi dth
lambda_partial is the root of B.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from flint import acb, arb, ctx, fmpq

I = acb(0, 1)


def rball(p: int, q: int = 1) -> arb:
    return arb(fmpq(p, q))


def rinterval(lo: fmpq, hi: fmpq) -> arb:
    return arb((lo + hi) / 2, (hi - lo) / 2)


def acos_analytic(x: acb, analytic: bool) -> acb:
    """Principal acos via logarithms, with branch checking for integration."""
    return -I * (x + I * (1 - x * x).sqrt(analytic=analytic)).log(
        analytic=analytic
    )


def h_and_h1(x: acb, analytic: bool) -> tuple[acb, acb]:
    """h(x)=acos(x)^2 and h'(x); the 0F1 form is regular at x -> 1."""
    beta = acos_analytic(x, analytic)
    S = (-(beta * beta) / 4).hypgeom_0f1(acb(3) / 2)   # sin(beta)/beta
    return beta * beta, -2 / S


def G(theta: acb, phi: acb, a: acb, analytic: bool) -> acb:
    """The r=1 integrand in N-form, including the sin(theta) measure."""
    s = theta.sin()
    c = theta.cos()
    ell = s * s + a * a * (c * c)
    C = a / ((a * a * (s * s) + c * c).sqrt(analytic=analytic)
             * ell.sqrt(analytic=analytic))

    u = s * phi.cos()
    m = 1 - u
    a2m1c2 = (a * a - 1) * (c * c)

    N = a2m1c2 + 2 * m                       # = ell * P, no cancellation
    K = a2m1c2 * u + m * (2 - m)             # = u N + m^2, no cancellation

    rootN = N.sqrt(analytic=analytic)
    rootL = ell.sqrt(analytic=analytic)

    x = C * m * rootL / rootN
    h0, h1 = h_and_h1(x, analytic)

    return (-u * h0 - C * rootL * h1 * m * K / (N * rootN)) * s


def integrate_2d(kernel, tol: arb, depth: int, limit: int) -> acb:
    two_pi = acb(2 * arb.pi())
    half_pi = acb(arb.pi() / 2)

    def outer(theta: acb, an_th: bool) -> acb:
        return acb.integral(
            lambda phi, an_ph: kernel(theta, phi, an_th and an_ph),
            0, two_pi,
            abs_tol=tol, rel_tol=tol,
            depth_limit=depth, eval_limit=limit,
        )

    return acb.integral(
        outer, 0, half_pi,
        abs_tol=tol, rel_tol=tol,
        depth_limit=depth, eval_limit=limit,
    ) / two_pi


def B(a_value, tol, depth, limit) -> acb:
    a = acb(a_value)
    return integrate_2d(lambda th, ph, an: G(th, ph, a, an), tol, depth, limit)


def record(v: acb) -> dict:
    return {
        "ball": str(v.real),
        "lower": str(v.real.lower()),
        "upper": str(v.real.upper()),
        "imag_contains_zero": bool(0 in v.imag),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dps", type=int, default=30)
    ap.add_argument("--tolerance", default="1e-10")
    ap.add_argument("--depth-limit", type=int, default=40)
    ap.add_argument("--eval-limit", type=int, default=300000)
    ap.add_argument("--lo", default="103/50")     # 2.06
    ap.add_argument("--hi", default="207/100")    # 2.07
    ap.add_argument("--json", type=Path,
                    default=Path("lambda_partial_bracket.json"))
    args = ap.parse_args()

    ctx.dps = args.dps
    tol = arb(args.tolerance)
    D, L = args.depth_limit, args.eval_limit

    def to_arb(text: str) -> arb:
        p, _, q = text.partition("/")
        return rball(int(p), int(q) if q else 1)

    b_lo = B(to_arb(args.lo), tol, D, L)
    b_hi = B(to_arb(args.hi), tol, D, L)

    pos = bool(b_lo.real.lower() > 0) and bool(0 in b_lo.imag)
    neg = bool(b_hi.real.upper() < 0) and bool(0 in b_hi.imag)

    out = {
        "target": "sign change of B(a) = F(1,a) bracketing lambda_partial",
        "bracket": f"[{args.lo}, {args.hi}]",
        "B_at_lo": record(b_lo),
        "B_at_hi": record(b_hi),
        "conditions": {"B(lo) > 0": pos, "B(hi) < 0": neg},
        "form": "cancellation-preserving N-form, ell*P = (a^2-1)c^2 + 2m",
        "arithmetic": "python-flint Arb/Acb ball arithmetic",
        "decimal_precision": args.dps,
        "integration_tolerance": args.tolerance,
    }
    out["status"] = ("SIGN_CHANGE_CERTIFIED" if pos and neg
                     else "FAILED_OR_INCONCLUSIVE")

    args.json.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
