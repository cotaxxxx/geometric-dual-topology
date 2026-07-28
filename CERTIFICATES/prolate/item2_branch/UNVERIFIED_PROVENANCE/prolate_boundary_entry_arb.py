#!/usr/bin/env python3
"""Validated Arb enclosure of the prolate boundary-entry function.

    B(lambda) = F(1,lambda) = d/dr E_lambda(r) |_{r=1}

lambda_partial is the root of B: the axis ratio at which the equatorial
stationary circle reaches the boundary r=1.

Contact-centred coordinates
---------------------------
Let chi be the angular distance from the contact direction e_1 and psi
the local azimuth about e_1.  With

    t = sin(chi/2),   u = cos(chi) = 1 - 2t^2,   m = 1 - u = 2t^2,
    q = sin^2(psi),   c^2 = 4 t^2 (1-t^2) q,     d = lambda^2 - 1,

    A = 1 + d (1-t^2) q,
    J = 1 + d (1-2t^2) q,
    W = lambda^2 - d c^2,

the two cancellation-preserving quantities of the (theta,phi) form,

    N = (lambda^2-1) c^2 + 2m,    K = (lambda^2-1) c^2 u + m(2-m),

factor exactly as

    N = 4 t^2 A,                  K = 4 t^2 (1-t^2) J,

so the quotient that carried the 0/0 at the contact point becomes

    m K / N^{3/2} = t (1-t^2) J / A^{3/2},

which is regular at t=0.  Also C sqrt(ell) = lambda / sqrt(W) and

    x = cos_alpha = lambda t / ( sqrt(W) sqrt(A) ).

Hence the fully regularized form

    B(lambda) = int_0^{pi/2} int_0^1 (4t/pi)
                  [ -u h(x) - (lambda/sqrt(W)) h'(x) t(1-t^2) J / A^{3/2} ]
                dt dpsi

Consistency check: applying the same change of variables to the boundary
energy integrand m h(x) gives (4t/pi) m h(x) = (8 t^3/pi) h(x), which is
exactly the kernel weight used by the item-5 boundary_energy routine.

Integration follows the item-5 pattern: four exact psi bands, each
enclosed independently, t inner on [0,1].
"""
from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path
from typing import Callable

import flint
from flint import acb, arb, ctx, fmpq


def rball(p: int, q: int = 1) -> arb:
    return arb(fmpq(p, q))


def rinterval(lo: fmpq, hi: fmpq) -> arb:
    return arb((lo + hi) / 2, (hi - lo) / 2)


def angle_data(x: acb) -> tuple[acb, acb]:
    """h(x)=acos(x)^2 and h'(x), regular at the alignment endpoint x=1."""
    one = acb(1)
    z = (one - x) / 2
    H = z.hypgeom_2f1(one / 2, one / 2, acb(3) / 2)
    h = 4 * z * H * H
    S = (-h / 4).hypgeom_0f1(acb(3) / 2)
    return h, -2 / S


def angle_data2(x: acb) -> tuple[acb, acb, acb]:
    """h, h', h'' via 2F1; h''(1) = 2/3 without cancellation."""
    one = acb(1)
    z = (one - x) / 2
    H = z.hypgeom_2f1(one / 2, one / 2, acb(3) / 2)
    H1 = z.hypgeom_2f1(acb(3) / 2, acb(3) / 2, acb(5) / 2) / 6
    H2 = z.hypgeom_2f1(acb(5) / 2, acb(5) / 2, acb(7) / 2) * 3 / 20
    h = 4 * z * H * H
    h1 = -2 * H * H - 4 * z * H * H1
    h2 = 4 * H * H1 + 2 * z * (H1 * H1 + H * H2)
    return h, h1, h2


def geometry(t: acb, psi: acb, lam, analytic: bool):
    """Shared contact-centred geometry.  lam may be acb or Dual-like."""
    t2 = t * t
    q = psi.sin() ** 2
    d = lam * lam - 1
    u = 1 - 2 * t2
    A = 1 + d * (1 - t2) * q
    J = 1 + d * (1 - 2 * t2) * q
    Wq = lam * lam - d * (4 * t2 * (1 - t2) * q)
    return t2, u, A, J, Wq


def kernel_B(t: acb, psi: acb, lam: acb, analytic: bool) -> acb:
    t2, u, A, J, Wq = geometry(t, psi, lam, analytic)
    rootW = Wq.sqrt(analytic=analytic)
    rootA = A.sqrt(analytic=analytic)
    x = lam * t / (rootW * rootA)
    h, h1 = angle_data(x)
    bracket = -u * h - (lam / rootW) * h1 * t * (1 - t2) * J / (A * rootA)
    return 4 * t * bracket / acb(arb.pi())


def kernel_Bprime(t: acb, psi: acb, lam_interval: arb, analytic: bool) -> acb:
    """d/dlambda of kernel_B, by explicit forward differentiation."""
    lam = acb(lam_interval)
    t2 = t * t
    q = psi.sin() ** 2
    d = lam * lam - 1
    u = 1 - 2 * t2

    A = 1 + d * (1 - t2) * q
    J = 1 + d * (1 - 2 * t2) * q
    c2 = 4 * t2 * (1 - t2) * q
    Wq = lam * lam - d * c2

    dA = 2 * lam * (1 - t2) * q
    dJ = 2 * lam * (1 - 2 * t2) * q
    dW = 2 * lam - 2 * lam * c2

    rootW = Wq.sqrt(analytic=analytic)
    rootA = A.sqrt(analytic=analytic)

    x = lam * t / (rootW * rootA)
    dx = x * (1 / lam - dW / (2 * Wq) - dA / (2 * A))

    h, h1, h2 = angle_data2(x)

    # term1 = -u h(x)
    dterm1 = -u * h1 * dx

    # term2 = -(lam/sqrt(W)) h1(x) t (1-t^2) J A^{-3/2}
    pref = lam / rootW
    dpref = 1 / rootW - lam * dW / (2 * Wq * rootW)
    g = t * (1 - t2) * J / (A * rootA)
    dg = t * (1 - t2) * (dJ / (A * rootA)
                         - acb(3) / 2 * J * dA / (A * A * rootA))
    dterm2 = -(dpref * h1 * g + pref * h2 * dx * g + pref * h1 * dg)

    return 4 * t * (dterm1 + dterm2) / acb(arb.pi())


def psi_partitioned_integral(
    kernel: Callable[[acb, acb, bool], acb], tol: arb, depth: int, limit: int
) -> acb:
    """Four exact psi bands on [0, pi/2]; t inner on [0,1]."""
    nbands = 4
    local_tol = tol / 8
    total = acb(0)
    pi = arb.pi()

    for j in range(nbands):
        lo = acb(pi * j / (2 * nbands))
        hi = acb(pi * (j + 1) / (2 * nbands))

        def outer(psi: acb, an_psi: bool) -> acb:
            return acb.integral(
                lambda t, an_t: kernel(t, psi, an_psi and an_t),
                0, 1,
                abs_tol=local_tol, rel_tol=local_tol,
                depth_limit=depth, eval_limit=limit,
            )

        total += acb.integral(
            outer, lo, hi,
            abs_tol=local_tol, rel_tol=local_tol,
            depth_limit=depth, eval_limit=limit,
        )
    return total


def B(lam_value: arb, tol: arb, depth: int, limit: int) -> acb:
    lam = acb(lam_value)
    return psi_partitioned_integral(
        lambda t, psi, an: kernel_B(t, psi, lam, an), tol, depth, limit
    )


def Bprime(lam_interval: arb, tol: arb, depth: int, limit: int) -> acb:
    return psi_partitioned_integral(
        lambda t, psi, an: kernel_Bprime(t, psi, lam_interval, an),
        tol, depth, limit,
    )


def record(v: acb) -> dict:
    return {
        "ball": str(v.real),
        "lower": str(v.real.lower()),
        "upper": str(v.real.upper()),
        "imag_ball": str(v.imag),
        "imag_contains_zero": bool(0 in v.imag),
    }


def positive(v: acb) -> bool:
    return bool(v.real.lower() > 0) and bool(0 in v.imag)


def negative(v: acb) -> bool:
    return bool(v.real.upper() < 0) and bool(0 in v.imag)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dps", type=int, default=30)
    ap.add_argument("--tolerance", default="1e-10")
    ap.add_argument("--depth-limit", type=int, default=30)
    ap.add_argument("--eval-limit", type=int, default=100000)
    ap.add_argument("--lo", default="103/50")
    ap.add_argument("--hi", default="207/100")
    ap.add_argument("--derivative-cells", type=int, default=10)
    ap.add_argument("--skip-derivative", action="store_true")
    ap.add_argument("--json", type=Path,
                    default=Path("certificate_item2_boundary_entry.json"))
    args = ap.parse_args()

    ctx.dps = args.dps
    tol = arb(args.tolerance)
    D, L = args.depth_limit, args.eval_limit

    def frac(text: str) -> fmpq:
        p, _, q = text.partition("/")
        return fmpq(int(p), int(q) if q else 1)

    lo, hi = frac(args.lo), frac(args.hi)

    b_lo = B(arb(lo), tol, D, L)
    b_hi = B(arb(hi), tol, D, L)

    conditions = {
        f"B({args.lo}) > 0": positive(b_lo),
        f"B({args.hi}) < 0": negative(b_hi),
    }

    deriv = []
    if not args.skip_derivative:
        n = args.derivative_cells
        for k in range(n):
            cl = lo + (hi - lo) * fmpq(k, n)
            cr = lo + (hi - lo) * fmpq(k + 1, n)
            v = Bprime(rinterval(cl, cr), tol, D, L)
            deriv.append({"cell": f"[{cl}, {cr}]", **record(v),
                          "certified_negative": negative(v)})
        conditions["B'(a) < 0 on every cell"] = all(
            item["certified_negative"] for item in deriv
        )

    out = {
        "target": "lambda_partial: unique root of B(lambda) = F(1,lambda)",
        "bracket": f"[{lo}, {hi}]",
        "form": "contact-centred t coordinates; N = 4t^2 A, K = 4t^2(1-t^2) J",
        "B_at_lo": record(b_lo),
        "B_at_hi": record(b_hi),
        "Bprime_cells": deriv,
        "conditions": conditions,
        "arithmetic": "python-flint Arb/Acb ball arithmetic",
        "integration": "four exact psi bands, t inner on [0,1]",
        "decimal_precision": args.dps,
        "integration_tolerance": args.tolerance,
        "environment": {
            "python": platform.python_version(),
            "python_flint": flint.__version__,
        },
    }
    out["status"] = ("CERTIFIED" if all(conditions.values())
                     else "FAILED_OR_INCONCLUSIVE")

    args.json.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({"status": out["status"],
                      "conditions": conditions,
                      "B_lo": out["B_at_lo"]["ball"],
                      "B_hi": out["B_at_hi"]["ball"]}, indent=2))


if __name__ == "__main__":
    main()
