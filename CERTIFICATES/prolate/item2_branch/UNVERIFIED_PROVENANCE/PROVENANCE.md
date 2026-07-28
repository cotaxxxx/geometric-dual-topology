# Unverified provenance — prolate_boundary_entry_arb.py

Recorded: 2026-07-26

## Facts

- SHA-256: `c076ef9d2279e4b55330011a3c11c337e6ba249b0d21f27c5a82174c3cdef1c3`
- Size: 8953 bytes
- mtime / ctime: 2026-07-26 00:58:11 UTC
- Owner uid: 0
- Not present in any branch of `cotaxxxx/geometric-dual-topology`
  (`git log --all -- '*boundary_entry*'` returns nothing)

The file was found pre-existing in the shared working directory. It
appeared approximately 18 minutes after the last file written in this
session. Its authorship was not established.

## Handling

`prolate_boundary_entry_arb.py` was independently cross-checked against
a separately written float64 implementation of the same functional
(`item2_explore.py`, whose normalization is pinned to the certified CAP
values of `Q(a)` to within 5.3e-15). Agreement was to 7 significant
digits at four evaluation points.

It is retained here as reference material only. It is **not** used as
the sole source of any certified conclusion, and it must not be placed
on the execution path of a final proof until its author is identified.

## Status of the lambda_partial result

The lambda_partial enclosure is separable from this provenance question:
every reported enclosure was cross-checked against the independent
implementation before being recorded. Nevertheless, before the result is
archived as CERTIFIED, the enclosures should be re-run from an
implementation of known authorship.

## Addendum — independent re-run from known-author implementation (2026-07-26)

The re-run recommended above was carried out. A from-scratch
implementation of the boundary-entry enclosure and of `B'` was written
without reading `prolate_boundary_entry_arb.py` (it was never opened,
grepped, imported, or placed on any execution path), using only the
published algebraic form in `../BOUNDARY_ENTRY_NFORM_NOTE.md` and
`../certificate_item2_boundary_entry.json`. The implementation and its
run record live in `../independent_recheck/`.

Result: all five certified claims were reproduced with identical signs
and mutually consistent enclosures (see
`../independent_recheck/certificate_item2_independent.json`):

- `B(103/50) > 0`, `B(207/100) < 0`, `B(206538/100000) > 0`,
  `B(206539/100000) < 0`, and `B'(lambda) < 0` on
  `[206538/100000, 206539/100000]`.
- Every enclosure overlaps the original with the same sign; none are
  disjoint. The independent `h''(1)` equals `2/3`.

Consequence for provenance: the `lambda_partial` result no longer depends
on this unknown-provenance file for any of its certified conclusions — an
independent implementation of known authorship now reproduces it. This
file's authorship remains unidentified; it stays here as reference
material only and is still not on any execution path. This addendum does
not delete or amend any statement above; it only records the re-run.
