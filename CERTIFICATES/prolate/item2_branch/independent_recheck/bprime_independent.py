"""
Step 4-3: independent rigorous certification of

    B'(lambda) < 0   for every lambda in [206538/100000, 206539/100000].

Forward differentiation in lambda is done with a Dual number carried over
Arb complex balls (acb): each quantity holds (value, d/dlambda). No finite
differences are used. The angle data reuses exactly the branch switch of
boundary_entry_independent.py (|1-x^2| < 0.9 -> 2F1 form, else acos form),
extended to the second derivative h''(x):

    h(x)   = acos(x)^2
    h'(x)  = -2 acos(x)/sqrt(1-x^2)         = -2 * 2F1(1/2,1/2;3/2; 1-x^2)
    h''(x) = (2 - 2 x acos(x)/sqrt(1-x^2)) / (1-x^2)
                                            = (2x/3) * 2F1(3/2,3/2;5/2; 1-x^2)

Both h' and h'' are regular at x = 1, with h''(1) = 2/3, matching the value
recorded in certificate_item2_boundary_entry.json.

Authorship: written from scratch from the published algebraic form. No file
under UNVERIFIED_PROVENANCE/ was consulted.

Certifying B'(lambda) < 0 for the whole interval at once: lambda is taken as
a single acb ball covering [206538/100000, 206539/100000]. Ball arithmetic
makes integrand_prime(t,psi) enclose the range of d/dlambda of the integrand
over that ball, so the rigorous integral encloses { B'(lambda) : lambda in
the interval }. If its upper bound is < 0, B' < 0 on the whole interval.
"""

import sys
import time
import json
from flint import acb, arb, ctx

HALF = THREEHALF = FIVEHALF = None


def _init_consts():
    global HALF, THREEHALF, FIVEHALF
    HALF = acb(1) / 2
    THREEHALF = acb(3) / 2
    FIVEHALF = acb(5) / 2


def _abs_upper(z):
    return float(z.abs_upper())


class Dual:
    """(value, d/dlambda) over acb, both acb balls."""
    __slots__ = ('v', 'd')

    def __init__(self, v, d):
        self.v = v
        self.d = d

    def __add__(s, o):
        if isinstance(o, Dual):
            return Dual(s.v + o.v, s.d + o.d)
        return Dual(s.v + o, s.d)
    __radd__ = __add__

    def __sub__(s, o):
        if isinstance(o, Dual):
            return Dual(s.v - o.v, s.d - o.d)
        return Dual(s.v - o, s.d)

    def __rsub__(s, o):
        return Dual(o - s.v, -s.d)

    def __mul__(s, o):
        if isinstance(o, Dual):
            return Dual(s.v * o.v, s.v * o.d + s.d * o.v)
        return Dual(s.v * o, s.d * o)
    __rmul__ = __mul__

    def __truediv__(s, o):
        if isinstance(o, Dual):
            return Dual(s.v / o.v, (s.d * o.v - s.v * o.d) / (o.v * o.v))
        return Dual(s.v / o, s.d / o)

    def __rtruediv__(s, o):
        return Dual(o / s.v, -o * s.d / (s.v * s.v))

    def __neg__(s):
        return Dual(-s.v, -s.d)

    def sqrt(s, analytic=True):
        r = s.v.sqrt(analytic=analytic)
        return Dual(r, s.d / (2 * r))


def _h_data(v, analytic):
    """Return (h, h', h'') at acb point v, valid on the whole ball."""
    z = 1 - v * v
    if _abs_upper(z) < 0.9:
        F = z.hypgeom_2f1(HALF, HALF, THREEHALF)
        G = z.hypgeom_2f1(THREEHALF, THREEHALF, FIVEHALF)
        h = z * F * F
        hp = -2 * F
        hpp = (2 * v / 3) * G
        return h, hp, hpp
    r = v.real
    if not (float(r.upper()) < 1.0 and float(r.lower()) > -1.0):
        nan = acb(arb('nan'))
        return nan, nan, nan
    b = v.acos()
    s = z.sqrt(analytic=analytic)
    h = b * b
    hp = -2 * b / s
    hpp = (2 - 2 * v * b / s) / z
    return h, hp, hpp


def h_dual(X, analytic):
    """Given x as a Dual, return (h(x), h'(x)) as Duals."""
    h, hp, hpp = _h_data(X.v, analytic)
    return Dual(h, hp * X.d), Dual(hp, hpp * X.d)


def integrand_prime(t, psi, lamD, analytic):
    """d/dlambda of the (4t/pi)-weighted integrand, as an acb."""
    sp = psi.sin()
    q = sp * sp
    t2 = t * t
    om = 1 - t2
    u = 1 - 2 * t2
    c2 = 4 * t2 * om * q
    d = lamD * lamD - 1
    A = 1 + d * (om * q)
    J = 1 + d * ((1 - 2 * t2) * q)
    W = lamD * lamD - d * c2
    sW = W.sqrt(analytic=analytic)
    sA = A.sqrt(analytic=analytic)
    x = (lamD * t) / (sW * sA)
    H, HP = h_dual(x, analytic)
    term = (-u) * H - (lamD / sW) * HP * (t * om) * J / (sA * sA * sA)
    out = (4 * t / acb.pi()) * term.d
    return out


def inner(psi, lamD, rel_tol, eval_limit, depth_limit):
    return acb.integral(
        lambda t, a: integrand_prime(t, psi, lamD, True),
        0, 1, rel_tol=rel_tol, eval_limit=eval_limit, depth_limit=depth_limit)


def Bprime(lam_ball, bands=4, rel_tol=None, eval_limit=8000, depth_limit=22):
    lamD = Dual(acb(lam_ball), acb(1))
    total = acb(0)
    hi = acb.pi() / 2
    for k in range(bands):
        a = hi * k / bands
        b = hi * (k + 1) / bands
        t0 = time.time()
        part = acb.integral(
            lambda psi, an: inner(psi, lamD, rel_tol, eval_limit, depth_limit),
            a, b, rel_tol=rel_tol, eval_limit=eval_limit, depth_limit=depth_limit)
        total += part
        print(f"    band {k + 1}/{bands}: {part.real.str(12, radius=True)}"
              f"  ({time.time() - t0:.1f}s)", flush=True)
    return total


if __name__ == '__main__':
    ctx.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    _init_consts()

    # --- cross-check h''(1) = 2/3 both branches meet the certificate value.
    h1, hp1, hpp1 = _h_data(acb(1), True)
    hpp1_re = hpp1.real
    target = arb(2) / 3
    ok = (hpp1_re - target).abs_upper() < arb(10) ** -15
    print(f"h''(1) = {hpp1_re.str(20, radius=True)}   "
          f"(target 2/3 = {target.str(20)}) -> {'OK' if ok else 'MISMATCH'}")
    if not ok:
        print("h''(1) MISMATCH -- stopping")
        sys.exit(1)

    # --- interval enclosure of B' over [206538/100000, 206539/100000].
    lo = arb(206538) / arb(100000)
    hi = arb(206539) / arb(100000)
    lam_ball = lo.union(hi)
    print(f"lambda interval ball: {lam_ball.str(20, radius=True)}")
    t0 = time.time()
    v = Bprime(lam_ball, bands=4, rel_tol=arb(2) ** -18)
    re = v.real
    up = arb(re.upper())
    sign = "NEGATIVE" if up < 0 else ("POSITIVE" if arb(re.lower()) > 0 else "UNDECIDED")
    print(f"B'([2.06538,2.06539]) = {re.str(15, radius=True)}  -> {sign}  "
          f"({time.time() - t0:.1f}s)")

    out = {
        "statement": "B'(lambda) < 0 for all lambda in [206538/100000, 206539/100000]",
        "h_second_derivative_at_1": {"value": hpp1_re.str(25, radius=True),
                                     "target": "2/3", "match": bool(ok)},
        "lambda_interval_ball": lam_ball.str(25, radius=True),
        "Bprime_ball": re.str(30, radius=True),
        "Bprime_upper": up.str(25),
        "Bprime_lower": arb(re.lower()).str(25),
        "sign": sign,
        "method": "Dual-over-acb forward differentiation; no finite differences; "
                  "psi split into 4 exact bands; t adaptive rigorous inner",
        "run_parameters": {"dps": ctx.dps, "rel_tol": "2^-18",
                           "eval_limit": 8000, "depth_limit": 22, "bands": 4},
    }
    with open("bprime_result.json", "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))
