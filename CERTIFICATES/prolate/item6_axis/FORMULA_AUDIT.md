# Item 6 formula audit

## Exact identities

For

\[
N=1-wc,
\quad
R^2=1-c^2+\lambda^2(c-w)^2,
\quad
S^2=1-c^2+c^2/\lambda^2,
\quad
C=N/(RS),
\]

symbolic differentiation gives

\[
\partial_w C
=C\left(-\frac{c}{N}+\frac{\lambda^2(c-w)}{R^2}\right).
\]

The simultaneous reflection

\[
(c,w)\mapsto(-c,-w)
\]

leaves \(N,R^2,S^2,C\) invariant. After changing the integration variable
\(c\mapsto-c\), this proves

\[
A_\lambda(-w)=A_\lambda(w),
\qquad
\Psi_\lambda(-w)=-\Psi_\lambda(w).
\]

## Independent numerical differentiation

At 30 decimal digits, direct evaluation of the integral formula for
\(\Psi_\lambda(w)\) agrees with numerical differentiation of
\(A_\lambda(w)\) at the following audit points:

| \(\lambda\) | \(w\) | absolute difference |
|---:|---:|---:|
| 2 | 0.4 | 0 at working precision |
| 4.7243834 | 0.6 | 0 at working precision |
| 10 | 0.8 | \(1.87\times10^{-30}\) |

This is an implementation audit only. It does not prove positivity on a box.

## Center Hessian kernel

At \(w=0\), put

\[
L=1+(\lambda^2-1)c^2,
\qquad
W^2=\lambda^2(1-c^2)+c^2,
\qquad
C_0=\frac{\lambda}{\sqrt L\,W}.
\]

The first two \(w\)-derivatives of the cosine are

\[
C_{w,0}
=
\frac{\lambda c(\lambda^2-1)(1-c^2)}{L^{3/2}W},
\]

and

\[
C_{ww,0}
=
\frac{\lambda^3(1-c^2)\left(2(\lambda^2-1)c^2-1\right)}{L^{5/2}W}.
\]

Therefore

\[
Q_\parallel(\lambda)=A_\lambda''(0)
=\int_{-1}^{1}q_\lambda(c)\,dc,
\]

with

\[
q_\lambda(c)
=
\frac12\left[
-2c\,h'(C_0)C_{w,0}
+h''(C_0)C_{w,0}^2
+h'(C_0)C_{ww,0}
\right].
\]

The removable endpoint values are

\[
h'(1)=-2,
\qquad
h''(1)=\frac23.
\]

For the sphere \(\lambda=1\), the kernel reduces exactly to

\[
q_1(c)=1-c^2,
\]

and hence

\[
Q_\parallel(1)=\int_{-1}^{1}(1-c^2)\,dc=\frac43.
\]

`prolate_axis_center_symbolic_audit.py` verifies all of these identities exactly. `prolate_axis_center_reference.py` finds positive values at the recorded samples through \(\lambda=10^6\), and the center-tail decade slope agrees with \(3\pi\) to relative error below \(1.7\times10^{-6}\). These numerical values do not certify \(Q_\parallel(\lambda)>0\) on an interval.

## Aspect-ratio tail coefficient

Put

\[
\mu=1/\lambda,
\qquad
s=w^2,
\qquad
H(\mu,s)=\frac{\Psi_{1/\mu}(\sqrt{s})}{\mu\sqrt{s}}.
\]

The large-aspect-ratio form is

\[
H(\mu,s)
=
3\pi\sqrt{1-s}\log(1/\mu)
+
\widehat B(\mu,s),
\]

or, before removing the positive factor \(\mu w\),

\[
\Psi_{1/\mu}(w)
=
\mu w\left[
3\pi\sqrt{1-w^2}\log(1/\mu)
+
\widehat B(\mu,w^2)
\right].
\]

### Exact outer/Laurent audit

Let \(d=c-w\). The first outer coefficient of the derivative kernel has opposite \(d^{-2}\) terms on the two sides of the moving layer. These cancel in the two-sided matching. The \(d^{-1}\) residues are

\[
\operatorname{Res}_{d=0^+}K_1
=
\frac{3\pi}{2}w\sqrt{1-w^2},
\]

and

\[
\operatorname{Res}_{d=0^-}K_1
=
-\frac{3\pi}{2}w\sqrt{1-w^2}.
\]

After integrating on both sides of the layer, their physical logarithmic contributions add to

\[
3\pi w\sqrt{1-w^2}\log(1/\mu).
\]

Thus the odd-factor quotient coefficient is exactly

\[
\widehat A(s)=3\pi\sqrt{1-s},
\qquad
\widehat A(0)=3\pi>0.
\]

`prolate_axis_tail_symbolic_audit.py` verifies these algebraic identities exactly with SymPy. This is a formal coefficient audit; it does **not** provide a uniform bound for the remainder \(\widehat B\).

### Non-certified seven-point regression audit

`prolate_axis_tail_reference.py` evaluates

\[
\frac{\lambda\Psi_\lambda(w)}{w}
\]

at \(\lambda=10^5\) and \(10^6\). The secant slope with respect to \(\log\lambda\) is compared with \(3\pi\sqrt{1-w^2}\) at

\[
w\in\{0.1,0.2,0.35,0.5,0.65,0.8,0.9\}.
\]

All seven relative errors are below \(1.7\times10^{-6}\), exceeding the requested five-significant-digit agreement. This remains floating-point reference evidence, not interval certification.

## Remaining proof obligations

Region C is not certified until \(Q_\parallel(\lambda)>0\) is interval-certified for all \(\lambda\ge1\) and a finite-\(w\) center-cap remainder is enclosed.

Region T is not certified until a rational tail strip is covered and a uniform validated lower bound for \(\widehat B(\mu,s)\) is proved. Near \(s=1\), where \(\widehat A(s)\) vanishes, the tail strip must overlap the independent pole-cap argument.
