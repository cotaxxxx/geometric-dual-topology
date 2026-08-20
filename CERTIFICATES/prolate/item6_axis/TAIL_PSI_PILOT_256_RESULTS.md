# Direct-Psi tail pilot result

Date: 2026-07-26

Run: `30193639845`

## Provenance

- combined output: ID `8629405586`, digest `2ba70fb2acddd0cb4075236ac41f082aa35a53618fc0cd367790902d9be8b653`
- low-w output: ID `8629404048`, digest `8a95a854a5fb6729a3c4b235cba0f1392f26fd7f65bd75c663c2e75c7cbe2686`
- high-w output: ID `8629404250`, digest `ebfcd19aebd3c95c01cce73fbd557a36a1fc217610ff0ada25dc35574f33784f`

The raw JSON hashes agree with the companion manifests.

## Certified high-w subnode

The block

`1/200 <= mu <= 3/400`, `5/8 <= w <= 3/4`

is CERTIFIED for `Psi_{1/mu}(w)>0`.

- evaluations: 211
- certified leaves: 106
- terminal boxes: 0
- exact leaf area: `1/3200`
- exact target area: `1/3200`
- worst rigorous lower endpoint: `0.00029904388710952540427558998498902167`
- worst leaf: `mu in [21/3200,11/1600]`, `w in [21/32,43/64]`, split depth 6

This is a genuine certified subnode. It is not the complete T-INTERFACE node.

## Incomplete low-w subnode

The block

`1/200 <= mu <= 3/400`, `1/20 <= w <= 1/8`

remains INCOMPLETE.

- evaluations: 256
- certified leaves: 0
- terminal boxes: 257
- terminal reason: `max_boxes_exhausted`
- exact terminal and target area: `3/16000`
- finite totals: 255
- nonfinite totals: 1, at the initial wide root box
- best observed radius: `0.06832751363981515`
- maximum observed finite lower endpoint: `-0.05438510869019467`

The remaining obstruction is interval width, not persistent NaN evaluation.

## Cost implication

The low block was evaluated through approximately four bisections in each parameter. Continuing to roughly seven bisections in each parameter gives total depth near 14 and a full-tree scale near 30000 evaluations. This agrees with the independent order-of-magnitude estimate and remains much smaller than the certified pole-modulus computation.

## DAG status

A certified subnode may be recorded for the first mu slab and high-w band. The other first-slab bands, the second mu slab through `mu=1/100`, the monotone extension, center and pole overlaps, and the final exact DAG remain open.

Item 6 remains NOT CERTIFIED.
