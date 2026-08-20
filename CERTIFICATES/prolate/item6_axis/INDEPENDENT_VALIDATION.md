# Item 6 center Hessian — independent numerical cross-check

Status: **PASSED / NON-CERTIFIED SUPPORTING CHECK**

This note records an independent numerical check supplied after the Arb certificate was produced. It is supporting evidence only and is not part of the interval proof.

## Coordinate convention

The certified quantity is the Hessian in the normalized axial coordinate:

\[
Q_\parallel^{(q)}(\lambda)=A_\lambda''(0),
\qquad
Q_\parallel^{(q)}=\lambda^2Q_\parallel^{(p)}.
\]

At \(\lambda=1\), the two coordinate conventions agree.

## Checks

| Check | Independent result |
|---|---:|
| sphere anchor \(Q_\parallel^{(q)}(1)=4/3\) | `1.33333333`, exact agreement |
| positivity on `[1,10]` | positive throughout; sampled minimum about `1.6075` at `lambda=10` |
| certified worst leaf `[79/8,10]` | pointwise values about `[1.6075,1.6169]`, consistent with rigorous lower bound `0.48955...` |

The `[1,10]` rigorous lower bound is therefore valid but deliberately conservative, with an observed margin loss of about a factor of `3.3` on the worst leaf.

## Large-aspect-ratio structure

The independent values reported were approximately

| \(\lambda\) | \(Q_\parallel^{(q)}(\lambda)\) |
|---:|---:|
| 20 | 1.112 |
| 50 | 0.614 |
| 100 | 0.372 |

They are consistent with

\[
Q_\parallel^{(q)}(\lambda)
\sim
\frac{3\pi\log\lambda}{\lambda}.
\]

At \(\lambda=100\), the leading asymptotic value is about `0.4340`, while the reported numerical value is about `0.3717`; the finite remainder is still negative at that scale.

## Design consequence

The independent check supports the following proof architecture:

- use compact Arb coverage only through the exact junction `lambda_0=100`;
- do not attempt to prove the entire unbounded parameter range by widening finite lambda boxes;
- transfer `lambda>=100` to Region T with `mu<=1/100` and

\[
H(\mu,s)
=
3\pi\sqrt{1-s}\log(1/\mu)+\widehat B(\mu,s).
\]

This note does not certify the asymptotic remainder, the finite center cap, or any full item 6 statement.
