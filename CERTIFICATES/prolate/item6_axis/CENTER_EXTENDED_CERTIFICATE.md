# Item 6 center Hessian — extended compact certificate

Status: **CERTIFIED**

## Certified statement

The validated Arb/Acb computation proves

\[
\boxed{Q_\parallel(\lambda)>0\qquad(1\le\lambda\le100).}
\]

Here

\[
Q_\parallel(\lambda)=A_\lambda''(0)
\]

is the axial Hessian eigenvalue in the normalized axial coordinate \(w\). This is the `q`-coordinate convention; in physical axial displacement coordinates, \(Q_\parallel^{(q)}=\lambda^2Q_\parallel^{(p)}\).

## Exact coverage record

- parameter interval: exact rational interval `[1,100]`
- evaluated boxes: 4135
- certified leaves: 2117
- terminal boxes: 0
- exact rational coverage: passed
- every accepted leaf: strict positive real lower bound
- every imaginary ball: contains zero

The smallest certified lower endpoint occurs on

\[
\lambda\in[231/4,925/16]=[57.75,57.8125]
\]

and is

\[
Q_\parallel([231/4,925/16])
>
0.00020046868958001537682007201789783354707398439156606.
\]

This lower bound is intentionally conservative. The certificate establishes only strict positivity; it is not a sharp enclosure of the pointwise minimum.

## Arithmetic and execution

- Python: 3.13.14
- python-flint: 0.9.0
- FLINT: 3.6.0
- decimal precision: 50 digits
- integration tolerance: `1e-18`
- integration depth limit: 22
- integration evaluation limit: 200000
- lambda split-depth limit: 12
- box evaluation limit: 32768

The driver uses exact rational lambda leaves. Each leaf is evaluated as one real parameter ball by validated Acb integration. A leaf is accepted only when the real lower endpoint is strictly positive and the imaginary enclosure contains zero. Failed leaves are bisected; certification requires exact contiguous coverage and zero terminal leaves.

## Integrity record

- workflow run: `30159579964`
- workflow artifact ID: `8619978829`
- workflow artifact digest / downloaded ZIP SHA-256: `affc28ab8afd169edc057c62d9a7e41a6a1ea868534280b18c336b5d27ab5c80`
- certificate JSON SHA-256: `4726bd601464094f68c1490a510d3b255249b9cbcfbff1df9a95194524deae73`
- SHA-256 manifest file: `eda34c19f6b19fec790c42aa6a91b6c0c0987a660ed6a8c6562b13cf26d72383`
- driver script SHA-256: `7de166c94c196fec299d78e3286f92fc2b45a2b5872e42c61820839e6d5f08c6`

The full 2117-leaf machine-readable JSON and its SHA-256 manifest are preserved in the named GitHub Actions artifact.

## Interface with the tail

This certificate fixes the exact finite/tail junction at

\[
\lambda_0=100,
\qquad
\mu_0=1/100.
\]

The compact center-Hessian strategy is not extended beyond this point. Since

\[
Q_\parallel(\lambda)\sim\frac{3\pi\log\lambda}{\lambda},
\]

the Hessian tends to zero and direct wide-lambda interval evaluation becomes increasingly dependency-limited. The region \(\lambda\ge100\) is therefore assigned to the normalized tail certificate for

\[
H(\mu,s)=\frac{\Psi_{1/\mu}(\sqrt{s})}{\mu\sqrt{s}}
=3\pi\sqrt{1-s}\log(1/\mu)+\widehat B(\mu,s).
\]

## Scope

This closes the center-Hessian statement on `[1,100]`. It does not yet certify:

- positivity of `Psi_lambda(w)/w` on a finite center cap;
- the compact interior region;
- the pole cap;
- a uniform lower bound for the tail remainder `Bhat`;
- the full item 6 theorem.

Accordingly, item 6 as a whole remains **NOT CERTIFIED**.
