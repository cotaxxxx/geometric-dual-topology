# Item 6 center Hessian — compact certificate

Status: **CERTIFIED**

## Certified statement

The validated Arb/Acb computation proves

\[
\boxed{Q_\parallel(\lambda)>0\qquad(1\le\lambda\le10).}
\]

Here

\[
Q_\parallel(\lambda)=A_\lambda''(0)
\]

is the axial Hessian eigenvalue at the center of the prolate spheroid.

## Exact coverage record

- parameter interval: exact rational interval `[1,10]`
- evaluated boxes: 151
- certified leaves: 80
- terminal boxes: 0
- exact rational coverage: passed
- every accepted leaf: strict positive real lower bound
- every imaginary ball: contains zero

The smallest certified lower endpoint occurs on

\[
\lambda\in[79/8,10]
\]

and is

\[
Q_\parallel([79/8,10])
>
0.48955236624670122276879181945444644042425864613135.
\]

## Arithmetic and execution

- Python: 3.13.14
- python-flint: 0.9.0
- FLINT: 3.6.0
- decimal precision: 50 digits
- integration tolerance: `1e-18`
- integration depth limit: 22
- integration evaluation limit: 200000
- lambda split-depth limit: 10
- box evaluation limit: 8192

The driver uses exact rational lambda leaves. Each leaf is evaluated as one real parameter ball by validated Acb integration. A leaf is accepted only when the real lower endpoint is strictly positive and the imaginary enclosure contains zero. Failed leaves are bisected; certification requires exact contiguous coverage and zero terminal leaves.

## Integrity record

- workflow run: `30159239685`
- workflow artifact ID: `8619781926`
- workflow artifact digest: `sha256:0b6292035dde9ce11c67b56564ebd9e4d3162e432ed17c94151bbcc0c1f9da23`
- certificate JSON SHA-256: `72b6cd989b7b9395f857b1f747d529d4fe640b86eb29fb4222505addb926750a`
- driver script SHA-256: `7de166c94c196fec299d78e3286f92fc2b45a2b5872e42c61820839e6d5f08c6`

The full 80-leaf machine-readable JSON and its SHA-256 file are preserved in the named GitHub Actions artifact. This repository record states only data read directly from that successful artifact; it does not reconstruct or alter the leaf certificate.

## Scope

This closes only the compact center-Hessian statement on `[1,10]`. It does not yet certify:

- `Q_parallel(lambda)>0` for `lambda>10`;
- positivity of `Psi_lambda(w)/w` on a finite center cap `0<w<=w0`;
- the compact interior region;
- the pole cap;
- the uniform aspect-ratio-tail remainder.

Accordingly, item 6 as a whole remains **NOT CERTIFIED**.
