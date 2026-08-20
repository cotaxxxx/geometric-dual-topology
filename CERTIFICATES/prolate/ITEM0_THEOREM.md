# Item 0 — complete equatorial exclusion theorem

## Theorem

Let `lambda_partial` be the certified boundary-entry parameter. For every

\[
1<\lambda<\lambda_\partial,
\qquad 0<r<1,
\]

one has

\[
\boxed{F(r,\lambda)>0}.
\]

Hence the equatorial plane has no noncentral stationary base point below the
boundary-entry threshold.

## Certified assembly

1. **0c, center band.** On `0 <= r <= 9/20`, `F_r>0`. Since
   `F(0,lambda)=0`, one gets `F>0` for `0<r<=9/20`.
2. **0d, middle band.** On `9/20 <= r <= 3/4`, `F>0` directly.
3. **0b, boundary band.** On `3/4 <= r <=1`, `F_r<0`. Thus
   `F(r,lambda)>=F(1,lambda)`.
4. **0a + Stage 1, right anchor.** For `lambda<lambda_partial`,
   `F(1,lambda)>0`.

The interfaces match at `r=9/20` and `r=3/4`; together they cover the whole
open radial interval.

## Certificate count

- 0c: 1363 leaves
- 0d: 224 leaves
- 0b: 435 leaves
- 0a + Stage 1: 22 + 4 leaves
- total: 2048 leaves

## Reusable certification principles

- Use inf-sup endpoint intervals rather than mid-radius balls when wide
  positive factors lose sign under centered multiplication.
- Preserve algebraic correlation before division; in item 0b the `W` form
  keeps `rho^2 >= b^2 L_2` visible and prevents four-order overestimation.
- Replace broad special-function balls by endpoint values plus monotonicity
  whenever derivative signs are available analytically.
- Permit mixed-label covers: direct `F>0` leaves, `F_r>0` leaves propagating
  from a left anchor, and `F_r<0` leaves propagating from a right anchor may
  interleave arbitrarily. The proof obligation is an exact complete cover
  together with a directed dependency graph in which every transfer leaf is
  reachable from a certified anchor or direct-positive leaf. Disjointness is
  not required; compatible certified overlaps are allowed.

The abstract form of these principles, including the mixed-label anchored
cover theorem, is recorded in
`01_GENERAL_THEORY/CERTIFICATION_ARCHITECTURE.md`.
