# Independent artifact audit — item 6 tail diagnostics

Date: 2026-07-26

## Scope

This record captures an independent verification of three item 6 diagnostic artifacts. It distinguishes:

- byte/hash and exact-rational consistency checks;
- finite fixed-point evaluations;
- incomplete interval-box diagnostics;
- heuristic estimates of the subdivision cost required for a full certificate.

Nothing in this record promotes item 6 or `T-INTERFACE` to `CERTIFIED`.

## Overall conclusion

The reported artifact contents are internally consistent, and the branch-level conclusion

`Item 6 remains NOT CERTIFIED`

is correct.

The phrase “the tail kernel has been finite-ized” must be interpreted narrowly:

- it is correct for the fixed-point kernel diagnostic;
- it is not correct for every wide parameter box;
- in interval evaluation, NaN enclosures remain in wide-`w` bulk charts and disappear after one subdivision in the inspected trace.

The accurate statement is:

> Endpoint-cap repairs succeed. Remaining NaNs are confined to wide bulk parameter boxes in the inspected diagnostic and are removed by subdivision; they have not been globally eliminated from interval evaluation.

## Integrity and exact-coverage checks — PASSED

### Signed first-band artifact

- the main JSON SHA-256 agrees with its recorded value; the independently reported digest begins `b1040cf7...`;
- the four script/driver/audit hashes in the companion `.sha256` file agree with the values declared inside the JSON;
- all nine terminal boxes were parsed as exact rational rectangles;
- their total exact area is

  `3/16000`,

  equal to the target rectangle area;
- the nine terminal rectangles are pairwise interior-disjoint;
- every rectangle lies inside the requested target domain.

Therefore the reported exact area invariant is correct.

### Six-band summary

For every `w` band, the reported `terminal_area` agrees with the exact rational area of that band, including the `1/3200`-scale slabs. The six-band coverage accounting is internally consistent.

## Fixed-point kernel diagnostic — FINITE AND POSITIVE AT THE TEST POINTS

At

`mu = 1/160`,

the three Acb component integrals are finite in all six tested `w` bands, with no NaN values.

The summed direct quantity is clearly positive throughout the six representative points. Representative totals include approximately:

- `+41.4` at `w=7/80`;
- `+30.5` at `w=11/16`.

The values decrease gradually across the tested points and remain on the scale `O(30–40)`.

This is useful reference information, but it is sampled fixed-point evidence rather than a uniform interval theorem.

## Remaining NaNs in interval evaluation

In the signed-angle interval diagnostic, the root box spanning a full `w` band can still produce NaN enclosures in bulk components, including the inspected fields

- `paired_bulk`;
- `paired_outer`;
- `far_bulk`;
- `far_left`.

The endpoint-cap components

- `paired_cap`;
- `far_cap`

are finite already at the root box. This confirms that the `C=±1` endpoint-cap repair is functioning as intended.

After one parameter subdivision, all inspected components become finite. Thus the remaining NaN mechanism in this diagnostic is associated with wide bulk parameter enclosures, not with the endpoint cap itself.

This does not constitute a certified counterexample. It also does not justify saying that NaNs have been removed from all interval evaluations.

## Heuristic subdivision-cost estimate — NOT A PROOF

The first-band trace indicates approximately first-order contraction of the total interval radius under parameter bisection. A rough fit from the inspected data is

`rad ≈ 5×10^5 mu_rad + 3.4×10^4 w_rad`.

For a box with approximately

- `mu_rad = 7.5×10^-4`;
- `w_rad = 2.5×10^-2`,

the observed total radius is about `1.6×10^3`.

Since the fixed-point values are approximately `30–40`, reducing the radius below about `30` appears to require roughly

- five to six additional bisections in `mu`;
- six to seven additional bisections in `w`.

This suggests an order of magnitude near `10^4` leaves per `w` band. Given that the certified pole-modulus node completed with about `1.6×10^6` evaluations, a compact-tail run with a practical box budget appears computationally plausible.

This is an empirical planning estimate only. Certification still requires exact coverage, all accepted leaves strictly positive, terminal zero, and reproducible manifests.

## Consequences for the proof DAG

The following states remain unchanged:

- `P-MODULUS`: `CERTIFIED`;
- fixed-point tail diagnostics: finite and positive sampled evidence;
- bounded six-band tail smoke: `INCOMPLETE`;
- `T-INTERFACE`: `NOT CERTIFIED`;
- `T-MONO`: `NOT CERTIFIED`;
- `T-INTERIOR-0`: `NOT CERTIFIED`;
- `T-CENTER` and `T-POLE`: `NOT CERTIFIED`;
- item 6 as a whole: `NOT CERTIFIED`.

The next numerical check should compare the completed `d7f08c9507a028c303601a69901880a14d5c1f49` direct-`Psi` smoke trace against the contraction-rate estimate above, then raise the box budget only after confirming the expected radius decay and absence of persistent evaluation errors.
