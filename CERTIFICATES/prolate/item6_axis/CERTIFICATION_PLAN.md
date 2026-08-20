# Item 6 certification plan

The proof will be assembled from four analytic regions and one dependency audit.

The exact finite/tail parameter interface is fixed at

\[
\lambda_0=100,
\qquad
\mu_0=1/100.
\]

The compact finite-parameter proof will not be extended past `lambda=100`; the large-aspect-ratio side is assigned to Region T.

## Region C — center cap

Use odd analyticity

\[
\Psi_\lambda(w)=w\widehat\Psi_\lambda(w^2)
\]

and certify \(\widehat\Psi_\lambda>0\) near \(w=0\).

### C0 — center Hessian anchor

The first coefficient is

\[
Q_\parallel(\lambda)=A_\lambda''(0).
\]

This anchor is now certified on the full compact parameter interval:

\[
Q_\parallel(\lambda)>0
\qquad(1\le\lambda\le100).
\]

The certificate uses 2117 exact-rational leaves, 4135 evaluations, exact coverage, and zero terminal boxes. It is the finite-side anchor for Region C.

### C1 — finite center cap

Use the identity

\[
\frac{\Psi_\lambda(w)}{w}
=
\int_0^1 A_\lambda''(tw)\,dt.
\]

It is sufficient to certify

\[
A_\lambda''(v)>0
\]

on a finite strip `0<=v<=w0`, `1<=lambda<=100`. Start with the already implemented strip `w0=1/20` and split the lambda direction into exact rational blocks. The previous monolithic `[1,10]` job was cancelled and supplies no theorem.

## Region I — compact interior

Validate the one-dimensional integral directly on rational \((w,\lambda)\)-boxes, using regularized angle derivatives and adaptive subdivision.

The exact compact parameter domain is

\[
1\le\lambda\le100.
\]

Region I must overlap Region C in `w` and Region P before `w=1`.

## Region P — pole cap

Introduce a blow-up coordinate near \((w,c)=(1,1)\). Preserve the algebraic correlation among \(N=1-wc\), \(R^2\), and \(1-C^2\) before any division. Use inf-sup arithmetic for wide positive factors.

Region P must overlap both Region I for `lambda<=100` and Region T for `lambda>=100`.

## Region T — aspect-ratio tail

Set

\[
\mu=1/\lambda,
\qquad
s=w^2.
\]

The tail parameter strip is

\[
0<\mu\le\mu_0=1/100.
\]

The unscaled derivative tends to zero as \(\mu\to0\); its leading size is \(\mu\log(1/\mu)\). Remove both positive factors by defining, for \(0<s<1\),

\[
H(\mu,s)
=
\frac{\Psi_{1/\mu}(\sqrt{s})}{\mu\sqrt{s}},
\]

with the \(s=0\) value supplied by odd analyticity. The correct tail decomposition is

\[
H(\mu,s)
=
\log(1/\mu)\,\widehat A(s)
+
\widehat B(\mu,s),
\qquad
\widehat A(s)=3\pi\sqrt{1-s}.
\]

Equivalently,

\[
\Psi_{1/\mu}(w)
=
\mu w\left[
3\pi\sqrt{1-w^2}\log(1/\mu)
+
\widehat B(\mu,w^2)
\right].
\]

The formal outer/Laurent audit proves the coefficient \(\widehat A\): the opposite \((c-w)^{-2}\) terms cancel in the two-sided matching, while the two \((c-w)^{-1}\) residues add to

\[
3\pi w\sqrt{1-w^2}.
\]

The remaining certification target is

\[
\log(1/\mu)\,3\pi\sqrt{1-s}
+
\widehat B(\mu,s)>0
\]

on the exact rational tail strip `0<mu<=1/100`, using a uniform validated lower bound for \(\widehat B\).

At `mu=1/100`, Region T must overlap the finite Region C/I proof. Region T must also join Region P before \(s=1\), where the leading coefficient vanishes.

At the center-tail corner,

\[
\widehat A(0)=3\pi>0,
\]

so the apparent \((w,\mu)\to(0,0)\) obstruction disappears after division by the positive factor \(\mu w\).

## Assembly D — dependency DAG

Every leaf must be one of:

- direct \(\Psi>0\);
- center-anchor transfer;
- pole-anchor transfer;
- parameter-direction transfer.

The certificate is accepted only if exact rational coverage has no gap and every transfer leaf is reachable from a direct or boundary anchor. The interfaces at `lambda=100`, the Region C/I `w` overlap, and the Region I/P/T overlaps must be represented explicitly in the final DAG.
