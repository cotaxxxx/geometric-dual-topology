# Item 6 interval-constructor re-audit gate

Status: **PASSED / BLOCKER CLOSED**

Date: 2026-08-08

## Baseline

The item-6 branch was independently reviewed through source head
`8eff42b26c8644de9e6047715e2bef19075e7605` before this repair sequence.

The historical runtime audit accepted endpoint **overlap** as evidence for the exact-rational
Arb interval constructor. That predicate was weaker than the stated conclusion that the
constructed ball encloses the full rational interval.

## Corrected acceptance rule

The production constructor

```text
arb(str((lo + hi) / 2), str((hi - lo) / 2))
```

is now audited by the strict requirement that the resulting ball contain the complete
`arb(str(lo))` and `arb(str(hi))` endpoint balls. Overlap is retained only as diagnostic
information and cannot produce `PASSED`.

The strengthened audit source is
`prolate_axis_interval_constructor_audit.py`.

## Re-audit result

GitHub Actions run `31228715053` executed the strengthened audit with
`python-flint==0.9.0` at source head
`8c3826175a4c39f6b9413bd6649cc833cae4a0bf`.

The run completed successfully. Every test case reported

```text
contains_lo = true
contains_hi = true
```

and every strictly positive audit interval excluded zero. No historical direct-fmpq
midpoint/radius constructor remained in production item-6 Python sources.

Artifact `9012994554` has digest

```text
sha256:bfe581161ce41655da84659e581cc6f4fd2ef06f83dd1bb5e8e2efb036c9c4f2
```

and the generated JSON reports `status=PASSED`. The tracked
`prolate_axis_interval_constructor_audit.json` has been replaced by that passing result.

Full provenance is recorded in
`INTERVAL_CONSTRUCTOR_REAUDIT_RECEIPT_2026-08-08.md`.

## Certification effect

The interval-constructor blocker identified on 2026-08-08 is closed. This re-audit does
**not** by itself certify any still-open item-6 DAG node and does not change the overall
item-6 theorem status.

`C-HESSIAN`, `C-1`, `P-BOUNDARY`, and `P-MODULUS` remain archived certified nodes. The
full axial theorem remains **NOT CERTIFIED** until the unresolved finite-grid and tail
nodes are closed and the dependency DAG is assembled without gaps.

## State-document rule

The `Current workflow state` paragraph in the older `STATUS.md` is historical because it
is pinned to `c2534aec269263a0a585c374ad5f25d71fae9651`. This re-audit record and later
explicit audit records take precedence for workstream state after that baseline.
