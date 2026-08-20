# Boundary-entry function `B(a) = F(1,a)` — cancellation-preserving N-form

Date: 2026-07-26
Status: **design note, not a certificate.** No rigorous enclosure of
`B(a)` is claimed here.

## Purpose

Item 2 needs `lambda_partial` rigorously enclosed as the root of

    B(a) := F(1,a) = d/dr E_a(r) |_{r=1}.

A direct Arb implementation fails at the contact point
`(theta,phi) = (pi/2, 0)`, where the base point touches the boundary.
This note gives the algebraic form that removes the failure, together
with the bounds it makes available.

## The failure

With `s=sin(th)`, `c=cos(th)`, `ell = s^2 + a^2 c^2`, `u = s cos(phi)`,
`m = 1-u`, the natural form is

    P  = 1 - 2u/ell + 1/ell
    x  = cos_alpha = C m / sqrt(P)
    x_r = C ( -u/sqrt(P) - m^2/(ell P sqrt(P)) )

Near the contact point, writing `tau = pi/2 - th` and `rho^2 = tau^2 + phi^2`,

    m ~ rho^2/2,        P ~ (a^2 tau^2 + phi^2)/1 ~ rho^2,
    x ~ rho,            x_r ~ -1/rho.

So `x_r` blows up, and it is the cone weight `m` in the term
`m h'(x) x_r` that cancels the blow-up. Ball arithmetic cannot see that
cancellation: `P` is evaluated as a difference of O(1) quantities, its
enclosure contains zero on any cell touching the corner, and `1/sqrt(P)`
becomes infinite.

This is the same failure mode recorded for item 0b, where `P` and
`1 - e_1^2` were overestimated by four orders of magnitude, and it was
cured there by the correlation-preserving `W`-form.

## The identities

Two exact identities remove every subtraction of comparable terms.

**Identity 1.**

    ell * P = (ell - 1) + 2m = (a^2 - 1) c^2 + 2 m  =:  N

*Proof.* `ell*P = ell - 2u + 1 = (ell-1) + 2(1-u)`, and
`ell - 1 = s^2 + a^2 c^2 - 1 = (a^2-1) c^2`.

Both summands are nonnegative for `a > 1`, so `N >= 0` with no
cancellation, and

    N >= (a^2-1) c^2,        N >= 2m.

**Identity 2.**

    u N + m^2 = (a^2 - 1) c^2 u + m (2 - m)  =:  K

*Proof.* Substitute `u = 1 - m` and expand.

## The integrand in N-form

    cos_alpha = x = C m sqrt(ell) / sqrt(N)

    G(th,phi,a) := d/dr[ (1-ru) h(cos_alpha) ]|_{r=1} * sin(th)
                 = ( -u h(x) - C sqrt(ell) h'(x) * m K / N^{3/2} ) * s

    B(a) = (1/2pi) int_0^{pi/2} int_0^{2pi} G dphi dth

Verified against the original form to full double precision at interior
points and at the contact point, and consistent with the `Q` convention
already cross-checked against the certified CAP values.

## Bounds the N-form makes available

**Corner bound.** Since `|u| <= 1` and `0 <= 2-m <= 2`,

    |K| <= (a^2-1) c^2 + 2m = N,

hence

    | m K / N^{3/2} |  <=  m / sqrt(N)  <=  m / sqrt(2m)  =  sqrt(m/2).

Numerically the ratio `|mK/N^{3/2}| / sqrt(m/2)` never exceeds `1.000`
over 4e5 random samples in `a in [1.2,5]`; the bound is sharp.

**Range of the angle.** All factors of `x` are nonnegative, so
`x >= 0`, hence `beta = acos(x) in [0, pi/2]` and

    h(x) = beta^2 <= pi^2/4,        |h'(x)| = 2 beta / sin(beta) <= pi.

**Envelope.** `C sqrt(ell) = a / sqrt(a^2 s^2 + c^2) in [1, a]`, so

    |G| <= |u| pi^2/4 + a pi sqrt(m/2).

Near the contact point `c -> 0` gives `C sqrt(ell) -> 1` and
`sqrt(m/2) <= 1`, so `|G| <= pi^2/4 + pi < 5.62` there. Numerically
`sup |G| = pi^2/4`, attained in the limit at the contact point itself.

**Consequence (corner cap).** A strip `th in [pi/2 - eps, pi/2]`
contributes at most `5.62 * eps` to `B(a)`. Since the sign of `B` must be
resolved at magnitude `~1e-3` for the coarse bracket and `~1.7e-6` for
separation from `206539/100000`, caps of `eps = 1e-5` and `eps = 3e-8`
respectively are already negligible.

## What still blocks a certificate

The N-form fixes the corner but not the parameter dependence. Measured
in this environment (python-flint 0.9.0, dps 25-30):

- The inner `phi` integral is rigorous and fast at any fixed real or
  complex `theta`, including `theta` within `1e-4` of `pi/2`.
- Nested adaptive integration (rigorous `theta` outer over rigorous
  `phi` inner) exhausts the depth limit at the contact corner and
  returns `nan`.
- Replacing the outer integral by an interval Riemann sum over `theta`
  balls fails on dependency: at `theta = 0.9 +/- 1e-3` the inner
  integral encloses to width `7.7e-2`, an amplification of about `77x`
  over the ball radius, and at radius `1e-2` it returns `nan`.
  Extrapolating, a naive outer sweep would need on the order of `1e6`
  cells for `1e-4` accuracy.

The missing ingredient is therefore a mean-value or Taylor form in
`theta` (and later in `a`), i.e. the same machinery as the four-band
phi-partitioned production kernel already used for item 5.

## Recommended next steps

1. Port `B(a)` onto the item-5 phi-partitioned integrator rather than
   naive nested adaptive integration.
2. Certify the coarse bracket first: `B(103/50) > 0`, `B(207/100) < 0`.
   Target magnitude `~1e-3`, three orders easier than the tight case.
3. Add `B'(a) < 0` on that bracket for uniqueness, by the same
   `Dual`-over-`acb` pattern used in `qprime_integrand`.
4. Only then tighten toward separating `lambda_partial` from
   `206539/100000`, which requires resolving `|B| ~ 1.7e-6`.

## Numerical reference values (float64, not certified)

    lambda_partial ~ 2.065382293627
    206539/100000  = 2.06539            (above lambda_partial by 7.7e-6)

    B(2.0653)   = +1.802e-5
    B(2.06535)  = +7.071e-6
    B(2.06539)  = -1.687e-6
    B(2.0654)   = -3.877e-6
    B(2.06)     = +1.181e-3
    B(2.07)     = -1.010e-3
