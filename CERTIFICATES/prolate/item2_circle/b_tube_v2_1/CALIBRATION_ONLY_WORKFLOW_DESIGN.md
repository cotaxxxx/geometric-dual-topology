# B-TUBE v2.1 — calibration-only workflow design

Status: **IMPLEMENTATION PRESENT; H1/H2 CORRECTION APPLIED; B-LOCAL DEPENDENCY UNPINNED; RUN BLOCKED**

Design commit: `4a1b12a2a1e4f89712c33bc554646b44190f6f5b`

Audited harness source: `CERTIFICATES/prolate/item2_circle/b_tube_v2_1/`

## 1. Purpose and non-purpose

This stage measures candidate numerical operating parameters for a later production
B-TUBE run. It may inspect Krawczyk margins, interval inflation, JOIN widths,
subdivision counts, and evaluation budgets. It does not certify a branch, alter
theorem endpoints, discharge a paper-level dependency, or emit a production
B-TUBE verdict.

The only permitted terminal states are:

- `CALIBRATION_COMPLETE`
- `CALIBRATION_INCOMPLETE`
- `CALIBRATION_FAILED`

Every `CERTIFIED_*` value and every production verdict field is forbidden in
calibration output.

The value `2/1` came from `SELFTEST_ONLY` material and is not a B-LOCAL result. It
is therefore absent from the calibration configuration and is not accepted by any
binding code path.

The current configuration is an explicitly nonbinding diagnostic profile:

- `mode: "DIAGNOSTIC_ONLY"`;
- exact `diagnostic_lambda_start = 21/10`;
- `binding_to_final_lambda_start: false`;
- `blocal_dependency.status: "UNPINNED"`;
- every unavailable B-LOCAL tuple member is exactly `null`.

The diagnostic left endpoint is strictly above the Stage-1 upper bracket
`206539/100000`. It is an engineering probe endpoint only. It is not a replacement,
approximation, rounding, or claim about the final B-LOCAL output.

A diagnostic result must always contain `recommendation: null`,
`state: "CALIBRATION_INCOMPLETE"`, and `coverage_claim: false`, even when every
recorded candidate diagnostic passes. The terminal endpoint remains exactly
`118/25`.

## 2. Source and dependency boundary

The implementation calls the audited v2.1 canonical-byte, dyadic, chain, affine,
Krawczyk, and JOIN primitives. It must not copy or weaken those rules.

The actual imported file bytes of

`CERTIFICATES/prolate/item2_circle/vendor/prolate_circle_F_cleanroom.py`

must hash to

`77e7a93c594ba66ac7d98df29ec3c03107b0c63962a5aa60f8503559082c10ac`.

The resolved imported path must be a regular non-symlink file inside the checkout.
The check uses the actual imported module path and file bytes, not an exported
constant, configured label, manifest name, or wrapper-module identity.

**F and F_r must both be supplied only by this same pinned file. Supplying F_r
from any other module or file is forbidden.** The implementation verifies that
both `F_arb` and `dFdr_arb` are defined by the single loaded module.

B-LOCAL/B-ENTRY is the critical-path dependency. The configuration reserves an
exact dependency tuple with these members:

- artifact ZIP SHA-256;
- certificate SHA-256;
- source-head SHA;
- configuration SHA-256;
- exact final `lambda_start`;
- exact machine conclusion;
- dependency status.

In the present source, the status is `UNPINNED` and all six unavailable values are
`null`. A later separately audited commit must replace them with actual B-LOCAL
artifact values, set the status to `PINNED`, set
`binding_to_final_lambda_start: true`, and implement exact tuple validation before
a binding calibration can execute. Merely changing a status string or inserting a
numerical endpoint is rejected.

The C-G terminal identity tuple remains frozen:

- artifact ZIP SHA-256 `c0f624a955657f906c09c45b016a92f7bcdfa70d26c2508efeb3f06dd7d27381`;
- source head `1e0f671c91798b9c044c04c7a4224a21e1e67830`;
- config SHA-256 `bb6a3655d335240549cbe1f6eec2a9e68e00219eb9c1a2be65796e2e342a0d17`;
- reference-kernel SHA equal to the production-kernel SHA above;
- paper/interface lemma `F_G_FIXED_SLICE_IDENTITY_V1`;
- exact match parameter `118/25`;
- exact bracket `(1/64,11/256)`.

Calibration may record endpoint diagnostics but may not emit the production MATCH
conclusion.

## 3. Immutable configuration

`config.calibration.json` is canonical JSON with no trailing newline, duplicate
keys, floating JSON numbers, BOM, or CR/LF. Its normative fields include:

- explicit diagnostic mode and exact diagnostic endpoint;
- explicit `binding_to_final_lambda_start: false`;
- the complete, currently unpinned B-LOCAL dependency tuple shape;
- exact `lambda_end = 118/25`;
- ordered unique dyadic parameter widths and tube radii;
- predictor refresh cadence;
- Arb working and checker precision with `checker_dps >= dps`;
- maximum cells, subdivisions, and evaluation budget;
- audited-source and design commits;
- production and C-G dependency pins;
- exact affine rule, schema, design version, and chain domain.

Candidate order is normative. The cross-product order is parameter-width order
followed by tube-radius order. Environment variables cannot replace normative
configuration values.

Calibration is fresh-only. Resume files, checkpoints, caches, prior output, and
pre-existing output directories are rejected.

## 4. Execution modes and evaluation protocol

The ordinary command

`python calibration.py run --out <path>`

is the binding path. It fails before kernel evaluation while B-LOCAL remains
unpinned. The workflow calls only this binding path.

The separate command

`python calibration.py run --diagnostic --out <path>`

is an explicit local diagnostic path. The temporary GitHub workflow does not pass
`--diagnostic`; therefore neither a branch push nor an approval-tag workflow can
silently convert the diagnostic profile into an authorized run.

For the diagnostic path, the runner processes `[21/10,118/25]` in deterministic
exact-rational cells. These records describe only the chosen diagnostic interval.
They make no statement about the missing interval beginning at the eventual
B-LOCAL endpoint and cannot support a complete calibration recommendation or a
coverage claim.

Predictor endpoint values are exact dyadics. The only affine rule is
`exact_endpoint_convex_hull_v1`; midpoint substitution for correlated interval
expressions is forbidden.

Each cell record contains exact parameter endpoints, exact predictor endpoints and
tube interval, residual and derivative enclosures, exact preconditioner,
reconstructed Krawczyk image, strict margins or a precise failure reason,
derivative-sign diagnostic, and evaluation/subdivision counts.

Each shared endpoint receives a separate exact JOIN intersection record and width.
A candidate diagnostic passes only when all cells satisfy strict Krawczyk
inclusion, the derivative enclosure is strictly negative, all JOINs have positive
width, and fixed budgets are respected. A passing diagnostic remains nonbinding.

## 5. Independent verification

The binding workflow invokes `calibration.py verify` in fresh Python processes
after the runner and after delivery. Both the standard verifier and the independent
full record-layout verifier:

1. parse configuration and result files through canonical-byte routines;
2. reject duplicate keys, floats, BOM, CR, final JSONL LF, and noncanonical bytes;
3. verify the chain over canonical record-object bytes, excluding JSONL linefeeds;
4. reconstruct width-major/radius-minor candidate order locally with their own
   explicit parsing and nested loops;
5. do not call or import the runner's `_candidate_pairs` helper or a shared
   candidate-pair reconstruction helper;
6. verify candidate completeness and ordered indices;
7. recompute the first passing candidate internally;
8. suppress that candidate in diagnostic mode and require exactly
   `recommendation: null`, `CALIBRATION_INCOMPLETE`, and `coverage_claim: false`;
9. require `machine_conclusion` to equal `{"real_analytic":false}`;
10. reject every `CERTIFIED_*` string and production verdict field.

The pre-delivery binding verifier additionally requires a fully pinned B-LOCAL
tuple. Thus diagnostic records cannot be promoted into a byte-closed binding
delivery.

## 6. In-run receipt byte closure

For a future authorized binding run, delivery is built in a new empty directory.
It copies canonical results, exact replay sources, the pinned production kernel,
requirement lock, design, and workflow; hashes each payload file in sorted order;
builds a deterministic ZIP; hashes the actual ZIP bytes; writes a canonical
receipt; and independently rechecks every referenced digest.

The platform outer artifact ZIP is transport only. No observer may repair or
complete the receipt later.

## 7. Authorization, security, and lifecycle

The temporary workflow has only this trigger:

```yaml
on:
  push:
    tags:
      - "btube-v2-1-calibration-approved-*"
```

No approval tag exists. The workflow invokes the binding command without
`--diagnostic`, so the current unpinned configuration fails closed before numerical
evaluation or artifact creation.

After B-LOCAL design, implementation, certification, tuple pinning, replacement
source audit, and separate run approval, authorization would require the exact tag

`btube-v2-1-calibration-approved-<40-character audited implementation SHA>`

pointing to the same commit. Before checkout, the job requires the tag suffix to
equal `github.sha`; after checkout it requires `git rev-parse HEAD` to equal that
SHA.

The workflow has only `contents: read`; checkout uses
`persist-credentials: false`; actions are commit-pinned; Python-FLINT is
version- and wheel-SHA-pinned and installed with `--require-hashes
--only-binary=:all:`. There is no dispatch or write-capable token path.

The workflow must be deleted before a result-bearing branch is merged to main. The
accepted payload retains exact workflow bytes and SHA for replay.

## 8. Required controls

The implementation contains fail-closed controls for production-kernel mismatch,
path escape, alternate derivative supply, C-G tuple mismatch, missing B-LOCAL pins,
false B-LOCAL promotion, diagnostic endpoint at or below the Stage-1 upper bracket,
any attempt to bind the diagnostic endpoint, terminal endpoint mismatch,
`checker_dps < dps`, duplicate or unordered candidates, noncanonical JSON,
forbidden affine paths, missing records, shared-helper candidate reconstruction,
non-deterministic or diagnostic recommendations, diagnostic coverage claims,
forbidden result vocabulary, stale inputs, workflow authorization mismatch, and
surviving temporary workflow files in a result merge.

Positive controls include precision equality, explicit binding-run failure,
explicit diagnostic-mode acceptance, exact `21/10` endpoint checking, null B-LOCAL
tuple enforcement, and synthetic record-layout fixtures whose verifiers remain
functional when the runner helper is patched to raise.

## 9. Whole-source self-scan and static-audit gate

Every Python file recursively under `b_tube_v2_1/` is self-scanned. The scanner
tokenizes source so comments and string literals do not create false positives and
constructs forbidden spellings by adjacent string fragments.

This correction authorizes static re-audit only. It does not authorize main merge,
creation of an approval tag, a calibration workflow run, a production
configuration, or a production B-TUBE run. The next critical path is B-LOCAL/
B-ENTRY design, implementation, and certification.
