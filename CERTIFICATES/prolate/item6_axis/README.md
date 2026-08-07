# Prolate item 6 — axial profile

Updated: 2026-07-26

Status: **PARTIALLY CERTIFIED / FULL AXIAL THEOREM NOT CERTIFIED**

## Target

For the prolate spheroid

\[
K_\lambda=\left\{(x,y,z):x^2+y^2+\frac{z^2}{\lambda^2}\le1\right\},
\qquad \lambda\ge1,
\]

put the base point on the symmetry axis at

\[
p_{\lambda,w}=(0,0,\lambda w),
\qquad -1<w<1,
\]

and write

\[
A_\lambda(w)=E_{K_\lambda}(p_{\lambda,w}),
\qquad
\Psi_\lambda(w)=A_\lambda'(w).
\]

The item 6 target is

\[
\boxed{\Psi_\lambda(w)>0\quad(\lambda\ge1,\ 0<w<1).}
\]

If completed, this proves that the axial profile is strictly increasing from the center to the pole and that the center is the only stationary base point on the symmetry axis.

## Exact one-dimensional reduction

With `c` the third coordinate on the unit sphere,

\[
A_\lambda(w)
=
\frac12\int_{-1}^{1}
(1-wc)\,h(C_\lambda(c,w))\,dc,
\qquad h(x)=\arccos^2x,
\]

where

\[
C_\lambda(c,w)
=
\frac{1-wc}
{\sqrt{1-c^2+\lambda^2(c-w)^2}
 \sqrt{1-c^2+c^2/\lambda^2}}.
\]

Reflection symmetry gives

\[
A_\lambda(-w)=A_\lambda(w),
\qquad
\Psi_\lambda(-w)=-\Psi_\lambda(w),
\qquad
\Psi_\lambda(0)=0.
\]

## Signed-angle representation

For `0<=w<1`, define

\[
N=1-wc,
\]

\[
X=\sqrt{1-c^2}\left(\lambda w-(\lambda-\lambda^{-1})c\right),
\qquad
\delta=\arctan(X/N).
\]

Since `N>0`,

\[
\arccos(C_\lambda)^2=\delta^2.
\]

The exact derivatives are

\[
\delta_w
=
\frac{\lambda\sqrt{1-c^2}}
{1-c^2+\lambda^2(c-w)^2},
\]

\[
\delta_{ww}
=
\frac{2\lambda^3\sqrt{1-c^2}(c-w)}
{\left(1-c^2+\lambda^2(c-w)^2\right)^2}.
\]

This gives cancellation-safe finite-domain Arb kernels without hypergeometric angle regularization. The exact audit output is `prolate_axis_signed_angle_symbolic_audit.json`.

## Certified nodes

| Node | Domain | Statement | Record |
|---|---|---|---|
| `C-HESSIAN` | `1<=lambda<=100`, `w=0` | `Q_parallel(lambda)>0` | `CENTER_EXTENDED_CERTIFICATE.md` |
| `C-1` | `1<=lambda<=10`, `0<w<=1/20` | `Psi>0` through `A_second>0` | `CENTER_CAP_CERTIFICATE_1_10.md` |
| `P-BOUNDARY` | `1<=lambda<=100`, `w=1^-` | `Phi(lambda)>0` | `POLE_BOUNDARY_CERTIFICATE_1_100.md` |
| `P-MODULUS` | `1<=lambda<=100`, `1-2^-24<w<1` | `Psi>0` from boundary floor minus modulus loss | `POLE_MODULUS_CERTIFICATE.md` |

### Pole modulus result

The boundary certificate gives

\[
\Phi(\lambda)\ge\frac{43}{5000}=0.0086.
\]

The four-region modulus certificate gives

\[
\left|\Psi_\lambda(1-u)-\Phi(\lambda)\right|
\le
0.008565210847426734\ldots
<0.0086
\]

for

\[
1\le\lambda\le100,
\qquad
0\le u\le2^{-24}.
\]

Therefore `Psi_lambda(w)>0` on the final pole layer. The certificate contains 1,614,752 evaluations, 807,396 certified leaves, exact rational coverage, and terminal 0.

Machine record: `prolate_axis_pole_modulus_split.json`.

## Finite-domain assembly

For `1<=lambda<=100`, the remaining exact decomposition is

- `F-CENTER`: `0<w<=1/2`, signed-angle `A_second>0` grid;
- `F-MIDDLE`: `1/2<=w<=3/4`, signed-angle direct `Psi>0` grid;
- `F-POLE`: `3/4<=w<=63/64`, signed-angle direct `Psi>0` grid;
- `P-DYADIC`: `63/64<=w<=1-2^-24`, signed-angle direct `Psi>0` grid;
- `P-MODULUS`: `1-2^-24<w<1`, **CERTIFIED**.

The grid combiners must verify exact requested endpoints, expected file count, pairwise non-overlap, exact rational coverage, every block `CERTIFIED`, and terminal 0.

## Unbounded aspect-ratio tail

Put

\[
\mu=\lambda^{-1},
\qquad
H(\mu,w)=\frac{\Psi_{1/\mu}(w)}{\mu w},
\qquad
M(\mu,w)=-\mu\,\partial_\mu H(\mu,w).
\]

The exact junction is

\[
\lambda=100
\quad\Longleftrightarrow\quad
\mu=\frac1{100}.
\]

The active tail obligations are

- `T-INTERFACE`: certify `H>0` on `1/200<=mu<=1/100`, `1/20<=w<=3/4`;
- `T-MONO`: certify `M>0` on `1/400<=mu<=1/200`, `1/20<=w<=3/4`;
- `T-INTERIOR-0`: extend the monotonicity argument to `mu=0`;
- `T-CENTER`: certify the center-tail overlap;
- `T-POLE`: certify the pole-tail overlap.

Exact endpoint limits and logarithmic coefficients are supporting formulas only until uniform remainder bounds close these nodes.

## Current proof state

The pole modulus workflow succeeded and its artifact is now committed permanently. The ordinary exact audit also succeeded. The finite grid and dyadic pole workflows remain queued or active. Tail diagnostics, pilots, and smoke runs are evidence only and are not used as proof nodes.

See:

- `STATUS.md` — current evidence-state ledger;
- `DEPENDENCY_DAG.md` — exact node and edge assembly;
- `CERTIFICATION_PLAN.md` — implementation plan;
- `POLE_MODULUS_CERTIFICATE.md` — completed final-pole node;
- `SHA256SUMS_POLE_MODULUS.txt` — delivered pole-modulus artifact manifest.

## Remaining obligations

1. certify and archive the combined finite grids `F-CENTER`, `F-MIDDLE`, and `F-POLE`;
2. certify and archive `P-DYADIC`;
3. certify `T-INTERFACE`;
4. certify `T-MONO` and `T-INTERIOR-0`;
5. certify `T-CENTER` and `T-POLE` by uniform remainder bounds;
6. close the dependency DAG at every exact interface, including `lambda=100`.

Only after all finite and tail nodes are joined with exact coverage and terminal 0 may item 6 be marked **CERTIFIED**.
