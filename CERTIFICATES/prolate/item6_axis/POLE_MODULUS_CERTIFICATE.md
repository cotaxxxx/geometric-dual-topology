# Item 6 pole-modulus certificate

Date: 2026-07-26

Status: **CERTIFIED NODE `P-MODULUS`**. Item 6 as a whole remains **NOT CERTIFIED**.

## Statement

For

- `1 <= lambda <= 100`,
- `1 - 2^-24 < w < 1`,
- `Phi(lambda) >= 43/5000`,

the certified four-region weighted modulus bound gives

`|Psi_lambda(w) - Phi(lambda)| <= 0.008565210847426734... < 43/5000 = 0.0086`.

Therefore

`Psi_lambda(w) > 0`

throughout the final pole layer.

## Rigorous bounds

- boundary floor: `43/5000`
- dyadic thickness: `epsilon = 2^-24 = 1/16777216`
- weighted `|A_lambda''|` bound: `[143700.39247282134784787450023088319074824270803725 +/- 8.47e-45]`
- absolute safe ceiling: `90177536/625`
- remaining absolute-bound margin: `[583.66512717865215212549976911680925175729196275 +/- 8.85e-45]`
- modulus loss: `[0.008565210847426733246318966164045524045720261814431 +/- 6.00e-52]`

The modulus loss is strictly smaller than the boundary floor.

## Exact partition

The proof partitions the regularized outer integral at

`t = 1/128`, `1/64`, `255/256`, and `1`.

Each region reports exact rational coverage, zero terminal boxes, and a deterministic accepted-leaf digest.

| Region | Evaluations | Certified leaves | Terminal boxes | Worst absolute upper |
|---|---:|---:|---:|---|
| inner `d=u^2*y^2` | 326872 | 163440 | 0 | `[7.0226482176879782499251090923866396423372845794733 +/- 3.89e-50]` |
| outer near `u=r*sqrt(2)*t` | 146356 | 73182 | 0 | `[139085.32156101032823229814945402151079282411660172 +/- 4.32e-45]` |
| outer junction signed chart `1/128<=t<=1/64` | 57440 | 28724 | 0 | `[142571.62240051379073027342926034431545434687801515 +/- 2.52e-45]` |
| outer middle factored cosine chart `1/64<=t<=255/256` | 996728 | 498368 | 0 | `[144311.48783467518071315177903456029566035893048480 +/- 3.15e-45]` |
| outer endpoint signed chart `255/256<=t<=1` | 87356 | 43682 | 0 | `[5.3406720183296976515360750799885572613717664055493 +/- 3.28e-51]` |

Totals: **1614752 evaluations**, **807396 certified leaves**, **0 terminal boxes**.

## Audits

- endpoint and factorization audit: `PASSED`
- exact outer split audit: `PASSED`
- all four outer slab bounds below the residual uniform ceiling
- weighted total below the absolute safe ceiling
- modulus loss below the boundary floor

## Environment and provenance

- workflow run: `30193776148`
- workflow head: `c2534aec269263a0a585c374ad5f25d71fae9651`
- artifact ID: `8629478702`
- artifact name: `prolate-item6-pole-modulus-split`
- artifact ZIP SHA-256: `b7a9480c6325bf6fa64421d75d3f478122b95a812ab0aaa749295386ca2f655e`
- Python: `3.13.14`
- python-flint: `0.9.0`
- FLINT: `3.6.0`

The machine-readable certificate is `prolate_axis_pole_modulus_split.json`. Script identities are recorded in `prolate_axis_pole_modulus_split.json.sha256` and inside the certificate JSON.

## Scope

This closes only node `P-MODULUS`. The signed dyadic pole node `P-DYADIC`, finite-grid nodes, and the unbounded-`lambda` tail nodes remain separate obligations.
