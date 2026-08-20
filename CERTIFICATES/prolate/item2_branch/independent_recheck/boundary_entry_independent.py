"""
Independent rigorous enclosure of the prolate boundary-entry function

    B(lambda) = F(1, lambda) = d/dr E_lambda(r) |_{r=1}

in the contact-centred coordinates verified by
verify_change_of_variables.py:

    B(lambda) = int_0^{pi/2} dpsi int_0^1 dt (4t/pi) * [ -u h(x)
                  - (lambda/sqrt(W)) h'(x) t (1-t^2) J / A^{3/2} ]

    u = 1 - 2 t^2,  q = sin^2(psi),  d = lambda^2 - 1,
    c^2 = 4 t^2 (1-t^2) q,  A = 1 + d (1-t^2) q,
    J = 1 + d (1-2t^2) q,   W = lambda^2 - d c^2,
    x = lambda t / (sqrt(W) sqrt(A)),
    h(x) = acos(x)^2,  h'(x) = -2 acos(x)/sqrt(1-x^2).

Authorship: written from scratch from the published algebraic form.
No file under UNVERIFIED_PROVENANCE/ was consulted.

Angle data uses two representations of h, h' so that the integrand is
analytic on the whole of [0,1]:

  * near x = 1 (z = 1-x^2 small):  h'(x) = -2 * 2F1(1/2,1/2;3/2;z),
                                   h(x)  = z * 2F1(...)^2,
    both analytic at z = 0, which removes the 0/0 at t = 1;
  * away from x = 1: the direct acos form.

Every branch selection is verified on the whole ball; if neither
representation can be certified analytic the routine returns a
non-finite ball, which forces the integrator to subdivide.
"""

import sys
import time
import json
from flint import acb, arb, ctx

HALF = None
THREEHALF = None
NAN = None


def _init_consts():
    global HALF, THREEHALF, NAN
    HALF = acb(1) / 2
    THREEHALF = acb(3) / 2
    NAN = acb(arb('nan'))


def _abs_upper(z):
    """Upper bound for |z| as a Python float."""
    return float(z.abs_upper())


def h_and_hp(x, analytic):
    """Return (h(x), h'(x)) as acb, valid on the whole ball x."""
    z = 1 - x * x
    if _abs_upper(z) < 0.9:
        # 2F1 is analytic on |z| < 1; the ball is verified inside it.
        F = z.hypgeom_2f1(HALF, HALF, THREEHALF)
        return z * F * F, -2 * F
    if analytic:
        # acos has cuts on (-inf,-1] and [1,inf); certify the whole ball
        # avoids them.  If it cannot be certified, return a non-finite
        # ball so that the integrator subdivides.
        r = x.real
        if not (float(r.upper()) < 1.0 and float(r.lower()) > -1.0):
            return NAN, NAN
    b = x.acos()
    s = z.sqrt(analytic=analytic)
    return b * b, -2 * b / s


def integrand(t, psi, lam, analytic):
    d = lam * lam - 1
    sp = psi.sin()
    q = sp * sp
    t2 = t * t
    om = 1 - t2
    u = 1 - 2 * t2
    c2 = 4 * t2 * om * q
    A = 1 + d * om * q
    J = 1 + d * (1 - 2 * t2) * q
    W = lam * lam - d * c2
    sW = W.sqrt(analytic=analytic)
    sA = A.sqrt(analytic=analytic)
    x = lam * t / (sW * sA)
    h, hp = h_and_hp(x, analytic)
    term = -u * h - (lam / sW) * hp * t * om * J / (sA * sA * sA)
    return 4 * t / acb.pi() * term


def inner(psi, lam, rel_tol, eval_limit, depth_limit):
    return acb.integral(
        lambda t, a: integrand(t, psi, lam, True),
        0, 1,
        rel_tol=rel_tol, eval_limit=eval_limit, depth_limit=depth_limit)


def B(lam, bands=4, rel_tol=None, eval_limit=40000, depth_limit=30):
    """Rigorous enclosure of B(lambda). psi is split into `bands` exact bands."""
    lam = acb(lam)
    total = acb(0)
    hi = acb.pi() / 2
    for k in range(bands):
        a = hi * k / bands
        b = hi * (k + 1) / bands
        t0 = time.time()
        part = acb.integral(
            lambda psi, an: inner(psi, lam, rel_tol, eval_limit, depth_limit),
            a, b,
            rel_tol=rel_tol, eval_limit=eval_limit, depth_limit=depth_limit)
        total += part
        print(f"    band {k+1}/{bands}: {part.real.str(12, radius=True)}"
              f"  ({time.time()-t0:.1f} s)", flush=True)
    return total


def report(name, val):
    re = val.real
    lo, up = re.lower(), re.upper()
    sign = "POSITIVE" if arb(lo) > 0 else ("NEGATIVE" if arb(up) < 0 else "UNDECIDED")
    print(f"{name:>28}  = {re.str(12, radius=True)}   -> {sign}")
    return {"ball": re.str(30, radius=True),
            "lower": arb(lo).str(25), "upper": arb(up).str(25),
            "sign": sign}


if __name__ == '__main__':
    ctx.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    _init_consts()
    targets = [("B(103/50)", acb(103) / 50),
               ("B(207/100)", acb(207) / 100)]
    out = {}
    for name, lam in targets:
        t0 = time.time()
        v = B(lam, rel_tol=arb(2) ** -30)
        out[name] = report(name, v)
        out[name]["seconds"] = round(time.time() - t0, 1)
        print(f"{'':>28}    ({out[name]['seconds']} s)")
    print(json.dumps(out, indent=2))
