#!/usr/bin/env python3
"""Third-order r-jet for the prolate functional.

Factorial-normalized jets  X^[k] = (1/k!) d^k X / dr^k, k = 0..3, so that
products are convolutions:  (XY)^[n] = sum_k X^[k] Y^[n-k].

From the energy jet:
    F     = E'   = 1! * E^[1]
    F_r   = E''  = 2! * E^[2]
    F_rr  = E''' = 3! * E^[3]
    J     = F_rr / r        (continuous extension J(0,a) = H_4(a))

Validation targets:
    F_r(0,a) = Q(a)                      (certified CAP values)
    J(r,a) -> H_4(a) as r -> 0
    F_rr from the jet == central difference of F_r

float64 only.  This is the reference implementation whose structure the
Arb port should follow; it is not a certificate.
"""
import numpy as np

# --------------------------------------------------------------- 0F1

def hyp0f1(b, x, terms=30):
    """0F1(;b;x) by series; |x| <= (pi/2)^2/4 here, so this converges fast."""
    out = np.ones_like(x)
    term = np.ones_like(x)
    for k in range(1, terms):
        term = term * x / ((b + k - 1) * k)
        out = out + term
    return out


def h_derivs(x):
    """h(x)=acos(x)^2 and h',h'',h''' , regular at x -> 1."""
    x = np.clip(x, -1.0, 1.0)
    beta = np.arccos(x)
    z = -(beta * beta) / 4
    S = hyp0f1(1.5, z)
    T = hyp0f1(2.5, z)
    U = hyp0f1(3.5, z)
    h0 = beta * beta
    h1 = -2.0 / S
    h2 = (2.0 / 3.0) * T / S**3
    h3 = (2.0 / 15.0) * U / S**4 - (2.0 / 3.0) * T**2 / S**5
    return h0, h1, h2, h3


# --------------------------------------------------------------- jets

class Jet3:
    """Factorial-normalized 3-jet in r."""
    __slots__ = ('c',)

    def __init__(self, c0, c1=0.0, c2=0.0, c3=0.0):
        self.c = [c0, c1, c2, c3]

    @staticmethod
    def _lift(o):
        return o if isinstance(o, Jet3) else Jet3(o)

    def __add__(s, o):
        o = Jet3._lift(o)
        return Jet3(*[s.c[i] + o.c[i] for i in range(4)])
    __radd__ = __add__

    def __neg__(s):
        return Jet3(*[-v for v in s.c])

    def __sub__(s, o):
        return s + (-Jet3._lift(o))

    def __rsub__(s, o):
        return Jet3._lift(o) + (-s)

    def __mul__(s, o):
        o = Jet3._lift(o)
        a, b = s.c, o.c
        return Jet3(
            a[0]*b[0],
            a[0]*b[1] + a[1]*b[0],
            a[0]*b[2] + a[1]*b[1] + a[2]*b[0],
            a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0],
        )
    __rmul__ = __mul__

    def compose(s, f0, f1, f2, f3):
        """Substitute s into a scalar function with the given derivatives."""
        p1, p2, p3 = s.c[1], s.c[2], s.c[3]
        return Jet3(
            f0,
            f1 * p1,
            f1 * p2 + 0.5 * f2 * p1 * p1,
            f1 * p3 + f2 * p1 * p2 + (f3 / 6.0) * p1**3,
        )

    def pow_alpha(s, alpha):
        p = s.c[0]
        f0 = p**alpha
        f1 = alpha * p**(alpha - 1)
        f2 = alpha * (alpha - 1) * p**(alpha - 2)
        f3 = alpha * (alpha - 1) * (alpha - 2) * p**(alpha - 3)
        return s.compose(f0, f1, f2, f3)


# ---------------------------------------------------------- quadrature

_CACHE = {}

def nodes(nth, nph):
    key = (nth, nph)
    if key not in _CACHE:
        xt, wt = np.polynomial.legendre.leggauss(nth)
        th = 0.25*np.pi*(xt+1.0); wth = 0.25*np.pi*wt
        xp, wp = np.polynomial.legendre.leggauss(nph)
        ph = np.pi*(xp+1.0); wph = np.pi*wp
        TH, PH = np.meshgrid(th, ph, indexing='ij')
        W = np.outer(wth, wph)/(2*np.pi)
        _CACHE[key] = (TH, PH, W)
    return _CACHE[key]


def energy_jet(r, a, nth=160, nph=160):
    """Return E^[0..3] at (r,a)."""
    TH, PH, W = nodes(nth, nph)
    s, c = np.sin(TH), np.cos(TH)
    ell = s*s + a*a*c*c
    C = a/(np.sqrt(a*a*s*s + c*c)*np.sqrt(ell))
    u = s*np.cos(PH)

    R = Jet3(np.full_like(u, float(r)), np.ones_like(u))
    one_ru = 1.0 - R*u
    P = 1.0 - (R*u)*(2.0/ell) + (R*R)*(1.0/ell)
    x = (one_ru * P.pow_alpha(-0.5)) * C
    hx = x.compose(*h_derivs(x.c[0]))
    Wint = (one_ru * hx) * s
    return [float(np.sum(W*Wint.c[k])) for k in range(4)]


def F_and_derivs(r, a, **kw):
    """Return F, F_r, F_rr at (r,a)."""
    e = energy_jet(r, a, **kw)
    return e[1], 2.0*e[2], 6.0*e[3]


def J(r, a, **kw):
    """J = F_rr / r.  Not valid at r = 0; use H4 there."""
    return F_and_derivs(r, a, **kw)[2]/r


if __name__ == '__main__':
    print('--- F_r(0,a) = Q(a) vs certified CAP ---')
    for a, ref in ((4.7, 0.00232242449301), (4.75, -0.00240811894676),
                   (4.72438, 3.2219804e-7), (4.72439, -6.2418318e-7)):
        _, fr, _ = F_and_derivs(0.0, a)
        print('  a=%-9s F_r(0)=%+.12e  cert=%+.12e  diff=%.1e'
              % (a, fr, ref, abs(fr-ref)))

    print()
    print('--- F_rr: 3-jet vs 中心差分 ---')
    for a in (2.0654, 3.0, 4.7):
        for r in (0.05, 0.4, 0.9):
            _, _, frr = F_and_derivs(r, a)
            h = 2e-4
            fd = (F_and_derivs(r+h, a)[1] - F_and_derivs(r-h, a)[1])/(2*h)
            print('  a=%-7s r=%-5s jet=%+.8e  差分=%+.8e  rel=%.1e'
                  % (a, r, frr, fd, abs(frr-fd)/abs(fd)))

    print()
    print('--- J(r,a) -> H_4(a) as r -> 0 ---')
    for a in (2.0654, 3.0, 4.0, 4.72438):
        vals = [J(r, a) for r in (0.2, 0.05, 0.01, 0.002)]
        print('  a=%-9s ' % a + '  '.join('%+.6f' % v for v in vals))
