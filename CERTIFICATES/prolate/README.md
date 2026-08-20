# Prolate global-classification certificates

長球（prolate spheroid）の大域分類に用いる機械可読証明書を保存する。

## Item 0 — equatorial positivity

For every

\[
1<\lambda<\lambda_\partial,
\qquad 0<r<1,
\]

the four certified stages give

\[
\boxed{F(r,\lambda)>0}.
\]

| Stage | Region | Certified quantity | Status | Leaves |
|---|---|---|---|---:|
| 0c | `0 < r <= 9/20` | `F_r > 0`, hence `F > 0` from `F(0,lambda)=0` | CERTIFIED | 1363 |
| 0d | `9/20 <= r <= 3/4` | `F > 0` directly | CERTIFIED | 224 |
| 0b | `3/4 <= r <= 1` | `F_r < 0`, hence `F(r,lambda) >= F(1,lambda)>0` | CERTIFIED | 435 |
| 0a / Stage 1 | `r = 1` | boundary sign and unique boundary-entry parameter | CERTIFIED | 22 + 4 |

Total certified leaves: **2048**.

The common finite-box range is

\[
1\le\lambda\le\frac{206539}{100000}=2.06539.
\]

The strict final theorem uses the Stage-1 boundary sign for
`lambda < lambda_partial`. Consequently, below boundary entry the equator
contains no stationary point other than the center.

## Directories

- `item0c_center/` — center-band completion record.
- `item0d_interior/` — archived package and independent audit for the middle band.
- `item0b_boundary/` — archived compact certificate and independent audit for the boundary band.

See `ITEM0_THEOREM.md` for the assembled theorem and proof interfaces.
Binary artifacts are immutable. Any replacement must use a new filename and
a new SHA-256 manifest.

## Numbering notes

- `ITEM4_NUMBERING_NOTE.md` — records item 4 as a reserved gap, explains why later item numbers are not renumbered, and separates confirmed facts from unresolved history.
