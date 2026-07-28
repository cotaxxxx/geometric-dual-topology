# Item 2 boundary-entry — independent re-derivation and certification

Date: 2026-07-26
Status: **CERTIFIED**

## Purpose and independence

This record discharges the re-run condition formerly recorded in
`../UNVERIFIED_PROVENANCE/PROVENANCE.md`.

The implementation was written from the published algebraic statements in
`../BOUNDARY_ENTRY_NFORM_NOTE.md` and `../certificate_item2_boundary_entry.json`
only. `UNVERIFIED_PROVENANCE/prolate_boundary_entry_arb.py` was not opened,
grepped, imported, or placed on an execution path.

## 1. Exact change of variables

Under

    chi = 2 asin(t),  u = 1 - 2t^2,
    c^2 = 4 t^2 (1-t^2) q,  q = sin^2(psi),

the symbolic audit verifies exactly

    N  = 4 t^2 A,
    K  = 4 t^2 (1-t^2) J,
    m K / N^(3/2) = t (1-t^2) J / A^(3/2),
    C sqrt(ell) = lambda / sqrt(W),

where

    A = 1 + (lambda^2-1)(1-t^2)q,
    J = 1 + (lambda^2-1)(1-2t^2)q,
    W = lambda^2 - (lambda^2-1)c^2.

The measure is independently derived as

    sin(chi) dchi = 4t dt,
    (1/2pi) * 2 * 4t = 4t/pi.

The original N-form bracket and the contact-centred bracket agree to relative
`7.7e-41` at dps 40. The bracket has finite limit `-pi^2/4` as `t -> 0`, so
the weighted integrand is regular and no corner cap is required.

## 2. Non-rigorous high-precision cross-check

At dps 30, independent `mpmath` quadrature gives

| lambda | independent value | original enclosure |
|---|---:|---:|
| `103/50` | `+0.001180508429` | `[0.0011805 +/- 1.91e-8]` |
| `207/100` | `-0.001009560444` | `[-0.0010096 +/- 5.05e-8]` |
| `206538/100000` | `+5.021966056e-7` | `[+5.0009e-7, +5.0430e-7]` |
| `206539/100000` | `-1.687327983e-6` | `[-1.68943e-6, -1.68522e-6]` |

This is a cross-check only; `mp.quad` is not a certificate.

## 3. Rigorous Arb re-run

Environment:

- Python 3.11.15
- python-flint 0.9.0
- SymPy 1.14.0
- mpmath 1.3.0

`boundary_entry_independent.py` uses nested rigorous `acb.integral`, with
`psi` split into four exact bands. Near `x=1`, the analytic hypergeometric
representations of `h`, `h'`, and `h''` are used; elsewhere the direct `acos`
form is used only after branch-cut separation is certified.

All four evaluations completed at dps 30 and `rel_tol = 2^-22`:

| quantity | rigorous enclosure | sign |
|---|---:|---|
| `B(103/50)` | `[0.001181 +/- 7.92e-7]` | positive |
| `B(207/100)` | `[-0.001010 +/- 7.46e-7]` | negative |
| `B(206538/100000)` | lower `2.0032845e-7`, upper `8.0406477e-7` | positive |
| `B(206539/100000)` | lower `-1.9892452e-6`, upper `-1.3854108e-6` | negative |

Only the sum of all four bands is used for a sign decision. The total runtime
for the four evaluations was 7276.2 seconds.

## 4. Rigorous derivative enclosure

`bprime_independent.py` carries a Dual number over `acb` and differentiates
forward in `lambda`; no finite differences are used. The whole interval

    [206538/100000, 206539/100000]

is represented by one ball. The result is

    B'(lambda) in [-0.2189796249644107, -0.2189252946784765]

for every lambda in the interval. Hence `B' < 0` throughout it. The analytic
branch gives `h''(1) = 2/3`, matching the original certificate.

## 5. Certified conclusion and scope

The independent implementation certifies

    lambda_partial in (206538/100000, 206539/100000)

and uniqueness of the root on that interval. Therefore

    lambda_partial < 206539/100000.

The provenance condition is fully discharged. This result covers the
boundary-entry parameter only. Item 2 proper — uniqueness of the interior
stationary branch, requiring the single sign change of `F_r` — remains open.
