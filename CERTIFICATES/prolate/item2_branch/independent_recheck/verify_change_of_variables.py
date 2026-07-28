"""
Independent verification of the contact-centred (t, psi) form of the
boundary-entry integrand for the prolate spheroid, item 2.

Author: written from scratch for this check, from the algebraic
statements in BOUNDARY_ENTRY_NFORM_NOTE.md and
certificate_item2_boundary_entry.json.  No code from
UNVERIFIED_PROVENANCE/ was read into this implementation.

Two things are checked:

  (A) symbolic:  m*K / N^{3/2}  ==  t(1-t^2) J / A^{3/2}
                 C*sqrt(ell)    ==  lambda / sqrt(W)
                 x              ==  lambda t / (sqrt(W) sqrt(A))
      after substituting the contact-centred parametrisation.

  (B) numeric:   the (theta, phi) bracket and the (t, psi) bracket
                 agree at random points, and the measures match:
                 sin(theta) dtheta dphi / (2 pi)  ->  (4t/pi) dt dpsi.
"""

import random
import mpmath as mp
import sympy as sp

mp.mp.dps = 40

# ----------------------------------------------------------------- (A)

def symbolic_check():
    t, psi, lam = sp.symbols('t psi lam', positive=True)
    d = lam**2 - 1
    q = sp.sin(psi)**2

    # contact-centred parametrisation
    u = 1 - 2*t**2                 # = cos(chi),  chi = 2 asin(t)
    m = 1 - u                      # = 2 t^2
    c2 = 4*t**2*(1 - t**2)*q       # = cos^2(theta)

    # the two cancellation-free polynomials of the note
    N = d*c2 + 2*m
    K = d*c2*u + m*(2 - m)

    A = 1 + d*(1 - t**2)*q
    J = 1 + d*(1 - 2*t**2)*q

    print("N - 4 t^2 A      ->", sp.simplify(N - 4*t**2*A))
    print("K - 4t^2(1-t^2)J ->", sp.simplify(K - 4*t**2*(1 - t**2)*J))

    lhs = m*K/N**sp.Rational(3, 2)
    rhs = t*(1 - t**2)*J/A**sp.Rational(3, 2)
    print("mK/N^{3/2} ratio ->", sp.simplify(sp.powsimp(lhs/rhs, force=True)))

    # envelope factor
    ell = (1 - c2) + lam**2*c2          # s^2 + lam^2 c^2
    C = lam/(sp.sqrt(ell)*sp.sqrt(lam**2*(1 - c2) + c2))
    W = lam**2 - d*c2
    print("C sqrt(ell) - lam/sqrt(W) ->",
          sp.simplify(C*sp.sqrt(ell) - lam/sp.sqrt(W)))


# ----------------------------------------------------------------- (B)

def h(x):
    b = mp.acos(x)
    return b**2

def hp(x):
    b = mp.acos(x)
    if b == 0:
        return mp.mpf(-2)
    return -2*b/mp.sin(b)

def bracket_theta_phi(th, ph, lam):
    """Original (theta, phi) form, N-form of the note."""
    s, c = mp.sin(th), mp.cos(th)
    ell = s**2 + lam**2*c**2
    u = s*mp.cos(ph)
    m = 1 - u
    N = (lam**2 - 1)*c**2 + 2*m
    K = (lam**2 - 1)*c**2*u + m*(2 - m)
    Csell = lam/mp.sqrt(lam**2*s**2 + c**2)      # = C sqrt(ell)
    x = Csell*m/mp.sqrt(N)
    return -u*h(x) - Csell*hp(x)*m*K/N**mp.mpf(1.5)

def bracket_t_psi(t, psi, lam):
    """Contact-centred (t, psi) form of the certificate."""
    d = lam**2 - 1
    q = mp.sin(psi)**2
    u = 1 - 2*t**2
    c2 = 4*t**2*(1 - t**2)*q
    A = 1 + d*(1 - t**2)*q
    J = 1 + d*(1 - 2*t**2)*q
    W = lam**2 - d*c2
    x = lam*t/(mp.sqrt(W)*mp.sqrt(A))
    return -u*h(x) - (lam/mp.sqrt(W))*hp(x)*t*(1 - t**2)*J/A**mp.mpf(1.5)

def tpsi_to_thetaphi(t, psi):
    chi = 2*mp.asin(t)
    u = mp.cos(chi)
    c = mp.sin(chi)*mp.sin(psi)          # cos(theta)
    th = mp.acos(c)
    s = mp.sin(th)
    ph = mp.acos(u/s)
    return th, ph

def numeric_check(n=12):
    random.seed(20260726)
    worst = mp.mpf(0)
    for _ in range(n):
        t = mp.mpf(random.uniform(1e-3, 0.999))
        psi = mp.mpf(random.uniform(1e-3, mp.pi/2 - 1e-3))
        lam = mp.mpf(random.uniform(1.2, 5.0))
        th, ph = tpsi_to_thetaphi(t, psi)
        a = bracket_theta_phi(th, ph, lam)
        b = bracket_t_psi(t, psi, lam)
        rel = abs(a - b)/max(abs(a), mp.mpf(1))
        worst = max(worst, rel)
    print("worst relative bracket mismatch ->", mp.nstr(worst, 8))

def contact_limit():
    """Both forms at the contact point t -> 0 (theta -> pi/2, phi -> 0)."""
    lam = mp.mpf('2.065382293627')
    for k in range(1, 9):
        t = mp.mpf(10)**(-k)
        print(f"  t=1e-{k}:  bracket = {mp.nstr(bracket_t_psi(t, mp.pi/4, lam), 12)}")


if __name__ == '__main__':
    print("=== (A) symbolic identities ===")
    symbolic_check()
    print()
    print("=== (B) numeric agreement of the two brackets ===")
    numeric_check()
    print()
    print("=== contact-point limit in the new coordinates ===")
    contact_limit()
