# Direct-`Psi` smoke comparison — item 6 tail

Date: 2026-07-26

## Artifact

Branch source head:

`d7f08c9507a028c303601a69901880a14d5c1f49`

GitHub Actions run: `30192590525`  
Artifact: `prolate-item6-tail-H-first-six-smoke-summary`  
Artifact ID: `8629119615`  
Artifact digest: `sha256:a3a79c0b2956478724bde63cfce46e45f3a975d62e2fa411bda7f8bf1a618206`

The artifact contains all six requested `w` bands.

## Formal status

- diagnostic status: `COMPLETE_DIAGNOSTIC`;
- six of six band files received;
- each band status: `INCOMPLETE`;
- evaluations per band: `8`;
- certified leaves per band: `0`;
- terminal boxes per band: `9`;
- terminal reason: `max_boxes_exhausted`;
- exact rational area invariant: true.

This is not a `T-INTERFACE` certificate.

## Effect of removing division by `mu*w`

The direct-`Psi` formulation removes the artificial `1/(mu*w)` amplification from the sign test. In the earlier scaled trace, interval radii were on the order of `10^3`. In this direct trace, finite total enclosures have radii on the order of `10^-1` to `1`.

The reduction is substantial and confirms that testing the equivalent quantity

`Psi_{1/mu}(w)>0`

is numerically preferable to enclosing `H=Psi/(mu*w)` directly.

## Wide-box NaN behavior

NaNs are not globally eliminated.

At the full root boxes, bulk chart components can remain nonfinite while endpoint-cap components are finite. The number of finite total enclosures among the first eight evaluations decreases toward the high-`w` bands:

| `w` band | finite totals among first 8 evaluations |
|---|---:|
| `[1/20,1/8]` | 7 |
| `[1/8,1/4]` | 7 |
| `[1/4,3/8]` | 6 |
| `[3/8,1/2]` | 5 |
| `[1/2,5/8]` | 3 |
| `[5/8,3/4]` | 1 |

Thus the earlier statement that one subdivision removes all inspected NaNs applies to the earlier signed-angle diagnostic trace, not uniformly to this direct-`Psi` run. In the highest `w` band, several early boxes remain nonfinite before the first finite total appears.

The endpoint-cap repair is still effective: the cap terms are finite where the wide bulk terms fail.

## Deepest finite enclosures in the bounded trace

The deepest finite record in each band has

- `mu` approximately in `[1/200,1/160]`;
- `w` width approximately `1/32`;
- imaginary interval exactly zero;
- a total real enclosure still crossing zero.

| `w` subinterval | enclosure midpoint | enclosure radius |
|---|---:|---:|
| approximately `[0.05,0.06875]` | `0.01418` | `0.46162` |
| approximately `[0.125,0.15625]` | `0.03364` | `0.51493` |
| approximately `[0.25,0.28125]` | `0.06228` | `0.51302` |
| approximately `[0.375,0.40625]` | `0.08786` | `0.49818` |
| approximately `[0.5,0.53125]` | `0.10861` | `0.47090` |
| approximately `[0.625,0.65625]` | `0.12224` | `0.43269` |

The midpoint progression is consistent with the independent fixed-point values after multiplication by `mu*w`: the expected direct-`Psi` signal is small near the lower-`w` edge and increases toward the upper bands.

These midpoints are not certified point values; they are quoted only to assess scale and contraction.

## Comparison with the subdivision estimate

The observed radii at the deepest bounded records are about `0.43–0.51`, while the expected positive signal scale is approximately `0.01–0.13` across the six bands.

Therefore:

- the direct formulation has removed the dominant scaling amplification;
- the bounded run is still several effective bisections short of sign separation;
- the lowest-`w` band is likely to control the required depth because its positive signal is smallest;
- the high-`w` band may control evaluation robustness because nonfinite bulk enclosures persist longer there.

A further reduction by a factor near `30–40` is needed in the lowest band if the current first-order contraction continues. This is compatible with roughly five to six additional effective halvings after a finite enclosure is reached. Because both parameters must be resolved, a practical leaf count from `10^3` to `10^4` per band remains plausible, but the eight-evaluation trace is too shallow to certify that cost model.

## Next computational step

Use a staged pilot rather than launching the final full cover immediately:

1. increase the per-band box budget enough to observe at least two more complete contraction levels;
2. record finite/nonfinite counts by split depth and the worst finite radius by depth;
3. verify that no chart retains persistent evaluation errors;
4. fit separate contraction rates for the lowest-`w` and highest-`w` bands;
5. only then set the final terminal-zero budget and depth.

Acceptance remains unchanged: exact rational coverage, every accepted leaf strictly positive, no evaluation errors used as leaves, terminal zero, and a reproducible hash manifest.
