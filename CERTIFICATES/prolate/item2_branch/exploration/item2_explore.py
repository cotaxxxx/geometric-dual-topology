#!/usr/bin/env python3
"""Exploration kernel for prolate item 2.

E_a(r) = (1/2pi) int_0^{pi/2} int_0^{2pi} (1 - r u) h(cos_alpha(r)) sin(th) dphi dth
  u        = sin(th) cos(phi)
  ell      = sin^2 th + a^2 cos^2 th
  C        = a / ( sqrt(a^2 sin^2 th + cos^2 th) * sqrt(ell) )
  P        = 1 - 2 r u / ell + r^2 / ell
  cos_alpha= C (1 - r u) P^{-1/2}
  h(x)     = arccos(x)^2

Convention fixed by matching Q(a) = d^2/dr^2 E_a(0) against the certified
Arb CAP values.  Float64 Gauss-Legendre; exploration only, not a certificate.
"""
import numpy as np

# ---------- h(x)=acos(x)^2 and derivatives, regular at x->1 ----------

def h_jet(x):
    """Return h, h', h'' at x = cos(beta), beta in [0, pi)."""
    x = np.clip(x, -1.0, 1.0)
    b = np.arccos(x)
    s = np.sin(b)
    small = b < 1e-4
    # S = sin b / b, T = 0F1(5/2; -b^2/4) = 3(sin b - b cos b)/b^3
    S = np.where(small, 1 - b**2/6 + b**4/120, np.divide(s, b, out=np.ones_like(b), where=~small))
    with np.errstate(invalid='ignore', divide='ignore'):
        T_gen = 3*(s - b*np.cos(b))/np.where(small, 1.0, b**3)
    T = np.where(small, 1 - b**2/10 + b**4/280, T_gen)
    h0 = b*b
    h1 = -2.0/S
    h2 = (2.0/3.0)*T/S**3
    return h0, h1, h2


# ---------- 2-jet in r ----------

class J:
    """Truncated Taylor jet (v, d, dd) in the variable r."""
    __slots__ = ('v', 'd', 'dd')

    def __init__(self, v, d=0.0, dd=0.0):
        self.v, self.d, self.dd = v, d, dd

    def __add__(s, o):
        o = o if isinstance(o, J) else J(o)
        return J(s.v+o.v, s.d+o.d, s.dd+o.dd)
    __radd__ = __add__

    def __neg__(s):
        return J(-s.v, -s.d, -s.dd)

    def __sub__(s, o):
        return s + (-(o if isinstance(o, J) else J(o)))

    def __rsub__(s, o):
        return (o if isinstance(o, J) else J(o)) + (-s)

    def __mul__(s, o):
        o = o if isinstance(o, J) else J(o)
        return J(s.v*o.v, s.d*o.v + s.v*o.d, s.dd*o.v + 2*s.d*o.d + s.v*o.dd)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = o if isinstance(o, J) else J(o)
        v = s.v/o.v
        d = (s.d - v*o.d)/o.v
        dd = (s.dd - 2*d*o.d - v*o.dd)/o.v
        return J(v, d, dd)

    def rpow_mhalf(s):
        """s^{-1/2}."""
        v = s.v**-0.5
        d = -0.5*v*s.d/s.v
        dd = -0.5*(d*s.d + v*s.dd)/s.v + 0.5*v*s.d*s.d/(s.v*s.v)
        return J(v, d, dd)


def compose_h(x):
    """h(x(r)) as a jet, given x as a jet."""
    h0, h1, h2 = h_jet(x.v)
    return J(h0, h1*x.d, h2*x.d*x.d + h1*x.dd)


# ---------- quadrature ----------

_CACHE = {}

def nodes(nth, nph):
    key = (nth, nph)
    if key not in _CACHE:
        xt, wt = np.polynomial.legendre.leggauss(nth)
        th = 0.25*np.pi*(xt+1.0)
        wth = 0.25*np.pi*wt
        xp, wp = np.polynomial.legendre.leggauss(nph)
        ph = np.pi*(xp+1.0)
        wph = np.pi*wp
        TH, PH = np.meshgrid(th, ph, indexing='ij')
        W = np.outer(wth, wph)/(2*np.pi)
        _CACHE[key] = (TH, PH, W)
    return _CACHE[key]


def E_jet(r, a, nth=160, nph=160):
    """Return (E, dE/dr, d2E/dr2) at (r,a)."""
    TH, PH, W = nodes(nth, nph)
    s, c = np.sin(TH), np.cos(TH)
    ell = s*s + a*a*c*c
    C = a/(np.sqrt(a*a*s*s + c*c)*np.sqrt(ell))
    u = s*np.cos(PH)

    R = J(r, 1.0, 0.0)
    one_ru = 1.0 - R*u
    P = 1.0 - (R*u)*(2.0/ell) + (R*R)*(1.0/ell)
    x = (one_ru * P.rpow_mhalf()) * C
    Wint = (one_ru * compose_h(x)) * s

    return (float(np.sum(W*Wint.v)),
            float(np.sum(W*Wint.d)),
            float(np.sum(W*Wint.dd)))


def F(r, a, **kw):
    """F(r,a) = dE/dr."""
    return E_jet(r, a, **kw)[1]


def Fr(r, a, h=1e-5, **kw):
    """dF/dr via central difference on the jet's first derivative."""
    return (E_jet(r+h, a, **kw)[1] - E_jet(r-h, a, **kw)[1])/(2*h)


def Q(a, **kw):
    return E_jet(0.0, a, **kw)[2]


if __name__ == '__main__':
    print('--- normalization check: Q(a) vs certified Arb CAP ---')
    for a, ref in ((4.7, 0.00232242449301), (4.75, -0.00240811894676),
                   (4.72438, 3.2219804e-7), (4.72439, -6.2418318e-7)):
        q = Q(a)
        print('  a=%-9s Q=%+.12e   certified=%+.12e   diff=%.2e'
              % (a, q, ref, abs(q-ref)))
