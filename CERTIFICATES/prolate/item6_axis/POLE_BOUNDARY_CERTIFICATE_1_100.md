# Item 6 Region P — pole-boundary anchor certificate on `[1,100]`

Status: **CERTIFIED**

## Certified statement

Define the one-sided axial derivative boundary profile

\[
\Phi(\lambda)
=
\lim_{w\to1^-}\Psi_\lambda(w).
\]

The validated Arb/Acb computation proves

\[
\boxed{
\Phi(\lambda)>0
\qquad(1\le\lambda\le100).
}
\]

## Regular boundary formula

With

\[
d=1-c,
\qquad
R_\partial=2+(\lambda^2-1)d,
\]

one has at `w=1`

\[
R^2=dR_\partial,
\]

and

\[
C_\partial
=
\frac{\sqrt d}{
\sqrt{R_\partial}
\sqrt{d(2-d)+(1-d)^2/\lambda^2}}.
\]

Combining the factor `N=1-wc` with `C_w` before taking the boundary limit gives

\[
N\frac{C_w}{C}
\longrightarrow
-c-\frac{\lambda^2d}{R_\partial}.
\]

Thus `Phi(lambda)` is represented by a regular one-dimensional integral and no division by `1-wc` remains in the boundary integrand. The algebraic reduction is independently checked by `prolate_axis_pole_boundary_symbolic_audit.py`.

## Exact block assembly

The parameter range was divided into 16 exact adjacent blocks:

\[
[1,2],[2,3],\ldots,[9,10],
[10,12],[12,16],[16,24],[24,36],
[36,52],[52,72],[72,100].
\]

Every block was certified independently by exact binary bisection. The combiner checked:

- every block status is `CERTIFIED`;
- every block has exact rational coverage;
- the blocks are exactly adjacent from `1` to `100`;
- all terminal-interval counts are zero.

## Completion data

- block certificates: 16
- evaluated intervals: 48
- certified leaves: 32
- terminal intervals: 0
- exact rational coverage: passed on every block
- exact block adjacency: passed

The smallest certified lower endpoint occurs on

\[
\lambda\in[1,9/8]
\]

and is

\[
\Phi([1,9/8])
>
0.0086853328086197058499308881284935494330221451543475.
\]

## Integrity record

- workflow run: `30176841736`
- workflow head: `c53930991475c6ef7d37569f0d7b689a22d3555e`
- combined artifact ID: `8624459353`
- combined artifact digest / downloaded ZIP SHA-256: `9a3fd7f21742ecc8201a36726135347d0da5d13826a3adfe7bf2dfd602c546f4`
- combined JSON SHA-256: `96e522fe849b0b35ff6fc91fc92ceff7552ab6e033685d0fc3ffc8a6931075b2`
- combiner script SHA-256 used by the run: `4a7a335c915347ea57836dc6fe7a3a872ca1b0ff058134eb83433eaf333cf9f8`

The full block JSON files are preserved in the 16 workflow artifacts. The combined artifact records their hashes, counts, exact adjacency, and worst certified leaf.

## Scope

This closes the pole-boundary anchor only. It does not yet prove positivity on a finite-width pole cap. The remaining Region P step is to certify a transfer from `w=1^-` into the interior, naturally by proving

\[
A_\lambda''(w)<0
\]

on a pole strip, which would imply

\[
\Psi_\lambda(w)\ge\Phi(\lambda)>0
\]

there.
