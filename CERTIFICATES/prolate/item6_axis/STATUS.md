# Item 6 status

Updated: 2026-07-26

## Final target

\[
\Psi_\lambda(w)=A_\lambda'(w)>0
\qquad(\lambda\ge1,\ 0<w<1).
\]

Item 6 as a whole remains **NOT CERTIFIED**. The finite pole modulus node is now certified, but the finite grids and the unbounded-aspect-ratio tail are not yet assembled.

## Certified nodes

### `C-HESSIAN`

\[
Q_\parallel(\lambda)>0
\qquad(1\le\lambda\le100).
\]

- 4135 evaluations
- 2117 exact-rational leaves
- terminal 0
- exact coverage
- worst rigorous lower endpoint:
  `0.00020046868958001537682007201789783354707398439156606`

Records: `CENTER_EXTENDED_CERTIFICATE.md`, `prolate_axis_center_extended_arb_summary.json`.

### `C-1`

\[
\Psi_\lambda(w)>0
\qquad
\left(1\le\lambda\le10,\ 0<w\le\frac1{20}\right).
\]

Certified through `A_lambda''(v)>0` on `0<=v<=1/20`.

- 9 exact adjacent lambda blocks
- 6203 evaluations
- 3106 leaves
- terminal 0
- worst rigorous lower endpoint:
  `0.00012893655173892051421755471961636338423925380545573`

Records: `CENTER_CAP_CERTIFICATE_1_10.md`, `prolate_axis_center_cap_1_10_summary.json`.

### `P-BOUNDARY`

\[
\Phi(\lambda)=\lim_{w\to1^-}\Psi_\lambda(w)>0
\qquad(1\le\lambda\le100).
\]

- 16 exact adjacent blocks
- 48 evaluations
- 32 leaves
- terminal 0
- worst rigorous lower endpoint:
  `0.0086853328086197058499308881284935494330221451543475`

The rational floor used by the final pole layer is

\[
\Phi(\lambda)\ge\frac{43}{5000}=0.0086.
\]

Records: `POLE_BOUNDARY_CERTIFICATE_1_100.md`, `prolate_axis_pole_boundary_1_100_summary.json`.

### `P-MODULUS`

\[
\Psi_\lambda(w)>0
\qquad
\left(1\le\lambda\le100,\ 1-2^{-24}<w<1\right).
\]

The certified weighted modulus estimate is

- `|A_lambda''|` weighted bound:
  `[143700.39247282134784787450023088319074824270803725 +/- 8.47e-45]`
- modulus loss at `epsilon=2^-24`:
  `[0.008565210847426733246318966164045524045720261814431 +/- 6.00e-52]`
- boundary floor: `43/5000 = 0.0086`

Hence the modulus loss is strictly smaller than the certified boundary floor.

The four-region cover reports:

- 1,614,752 evaluations
- 807,396 certified leaves
- exact rational coverage in every region
- terminal 0
- endpoint/factorization audit `PASSED`
- exact split audit `PASSED`

Permanent records:

- `POLE_MODULUS_CERTIFICATE.md`
- `prolate_axis_pole_modulus_split.json`
- `prolate_axis_pole_modulus_split.json.sha256`
- `SHA256SUMS_POLE_MODULUS.txt`

Workflow provenance:

- run `30193776148`
- source head `c2534aec269263a0a585c374ad5f25d71fae9651`
- artifact `8629478702`
- artifact ZIP SHA-256 `b7a9480c6325bf6fa64421d75d3f478122b95a812ab0aaa749295386ca2f655e`

## Exact signed-angle reduction

For `0<=w<1`, `-1<=c<=1`, and `lambda>0`, define

\[
N=1-wc,
\]

\[
X=\sqrt{1-c^2}\left(\lambda w-(\lambda-\lambda^{-1})c\right),
\]

\[
\delta_\lambda(c,w)=\arctan\frac{X}{N}.
\]

Because `N>0`,

\[
\arccos(C_\lambda(c,w))^2=\delta_\lambda(c,w)^2.
\]

The exact derivatives are

\[
\partial_w\delta
=
\frac{\lambda\sqrt{1-c^2}}
{1-c^2+\lambda^2(c-w)^2},
\]

\[
\partial_w^2\delta
=
\frac{2\lambda^3\sqrt{1-c^2}(c-w)}
{\left(1-c^2+\lambda^2(c-w)^2\right)^2}.
\]

The exact audit output is `prolate_axis_signed_angle_symbolic_audit.json`.

## Finite domain `1<=lambda<=100`

The active decomposition is

- `F-CENTER`: `0<w<=1/2`, signed-angle `A_second>0` grid;
- `F-MIDDLE`: `1/2<=w<=3/4`, signed-angle direct `Psi>0` grid;
- `F-POLE`: `3/4<=w<=63/64`, signed-angle direct `Psi>0` grid;
- `P-DYADIC`: `63/64<=w<=1-2^-24`, 144 signed-angle dyadic rectangles;
- `P-MODULUS`: `1-2^-24<w<1`, **CERTIFIED**.

The finite and dyadic combiners must verify requested endpoints, expected file counts, pairwise non-overlap, exact rational coverage, every block `CERTIFIED`, and terminal 0.

## Tail `lambda>=100`

Set

\[
\mu=\lambda^{-1},
\qquad
H(\mu,w)=\frac{\Psi_{1/\mu}(w)}{\mu w},
\qquad
M(\mu,w)=-\mu\,\partial_\mu H(\mu,w).
\]

The exact finite/tail junction is `lambda=100`, equivalently `mu=1/100`.

Active direct targets are

- `T-INTERFACE`: `H>0` on `1/200<=mu<=1/100`, `1/20<=w<=3/4`;
- `T-MONO`: `M>0` on `1/400<=mu<=1/200`, `1/20<=w<=3/4`.

The exact endpoint and logarithmic-coefficient identities are supporting formulas only until uniform remainder bounds close `T-CENTER`, `T-INTERIOR-0`, and `T-POLE`.

## Current workflow state

At source head `c2534aec269263a0a585c374ad5f25d71fae9651`:

- axial formula audit: success;
- pole modulus certificate: success and archived;
- finite grid certificate: queued;
- dyadic pole layer certificate: queued;
- tail H kernel diagnostic: success, diagnostic only;
- direct-Psi tail pilot: success, pilot only;
- first-six tail H smoke: success, smoke only.

A successful diagnostic, pilot, or smoke run is not a proof node unless its exact domain is assembled into a certified combiner output.

## Remaining obligations

1. receive and archive `CERTIFIED` combined outputs for `F-CENTER`, `F-MIDDLE`, and `F-POLE`;
2. receive and archive the `CERTIFIED` combined output for `P-DYADIC`;
3. certify the positive signed-tail `T-INTERFACE` node;
4. certify `T-MONO` and extend the monotonicity argument to `mu=0` in `T-INTERIOR-0`;
5. certify the center and pole tail overlaps `T-CENTER` and `T-POLE` by uniform remainder bounds;
6. assemble the exact dependency DAG with no gaps, no overlaps, exact endpoints, and terminal 0.

The final theorem must remain **NOT CERTIFIED** until every obligation above is closed.
