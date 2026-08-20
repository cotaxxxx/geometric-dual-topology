#!/usr/bin/env python3
"""Component and stage diagnostic for the compact-tail H kernel.

This is not a positivity certificate.  It evaluates every algebraic and
special-function stage at exact rational point parameters and fixed t points,
then runs each Acb integral separately.  The purpose is to identify the first
non-finite stage behind a NaN interval integral.
"""
from __future__ import annotations

import json
from pathlib import Path

from flint import acb, arb, ctx, fmpq

import prolate_axis_tail_H_block_arb as base


ctx.dps = 50
RADIUS = acb(4)
TOLERANCE = arb("1e-20")
DEPTH = 18
EVAL_LIMIT = 100000


def point(value: fmpq | int) -> acb:
    return acb(arb(str(value)))


def finite(value: object) -> bool:
    text = str(value).lower()
    return not any(token in text for token in ("nan", "inf", "unknown"))


def value_record(value: object) -> dict:
    record = {"value": str(value), "finite": finite(value)}
    if isinstance(value, acb):
        record.update({
            "real": str(value.real),
            "imag": str(value.imag),
            "imag_contains_zero": bool(0 in value.imag),
        })
    return record


def run_stage(name: str, function, stages: list[dict]):
    try:
        value = function()
        stages.append({"stage": name, **value_record(value)})
        return value
    except Exception as exc:
        stages.append({
            "stage": name,
            "finite": False,
            "error": f"{type(exc).__name__}: {exc}",
        })
        return None


def outer_detail(
    mu: acb,
    w: acb,
    c: acb,
    difference: acb,
    rho2: acb,
    n: acb,
    analytic: bool,
) -> dict:
    stages: list[dict] = []
    mu2 = mu * mu
    p = run_stage("P", lambda: difference * difference + mu2 * rho2, stages)
    s = run_stage("S", lambda: rho2 + mu2 * c * c, stages)
    if p is None or s is None:
        return {"stages": stages}
    root = run_stage("sqrt(P*S)", lambda: (p * s).sqrt(analytic=analytic), stages)
    if root is None:
        return {"stages": stages}
    cosine = run_stage("cosine", lambda: mu * n / root, stages)
    if cosine is None:
        return {"stages": stages}
    angle = run_stage("regular_angle_data", lambda: base.regular_angle_data(cosine), stages)
    hbar = run_stage("regular_hbar", lambda: base.regular_hbar(cosine), stages)
    bracket = run_stage("bracket", lambda: -c + n * difference / p, stages)
    if angle is not None:
        _, h1 = angle
        stages[-3].update({
            "h": str(angle[0]),
            "h1": str(h1),
            "h_finite": finite(angle[0]),
            "h1_finite": finite(h1),
        })
    else:
        h1 = None
    if hbar is not None and bracket is not None and h1 is not None:
        run_stage(
            "density",
            lambda: n * (-c * hbar + h1 * bracket) / (2 * w * root),
            stages,
        )
    return {"stages": stages}


def inner_detail(mu: acb, w: acb, t: acb, analytic: bool) -> dict:
    stages: list[dict] = []
    mu2 = mu * mu
    q = run_stage("q", lambda: 8 * t - 4, stages)
    if q is None:
        return {"stages": stages}
    c = run_stage("c", lambda: w + mu * q, stages)
    if c is None:
        return {"stages": stages}
    rho2 = run_stage("rho2", lambda: (1 - c) * (1 + c), stages)
    reduced_p = run_stage("reduced_P", lambda: q * q + rho2, stages)
    s = run_stage("S", lambda: rho2 + mu2 * c * c, stages)
    if rho2 is None or reduced_p is None or s is None:
        return {"stages": stages}
    root = run_stage(
        "sqrt(reduced_P*S)",
        lambda: (reduced_p * s).sqrt(analytic=analytic),
        stages,
    )
    n = run_stage("N", lambda: 1 - w * w - mu * w * q, stages)
    if root is None or n is None:
        return {"stages": stages}
    cosine = run_stage("cosine", lambda: n / root, stages)
    if cosine is None:
        return {"stages": stages}
    angle = run_stage("regular_angle_data", lambda: base.regular_angle_data(cosine), stages)
    hbar = run_stage("regular_hbar", lambda: base.regular_hbar(cosine), stages)
    bracket = run_stage(
        "bracket", lambda: -c + n * q / (mu * reduced_p), stages
    )
    if angle is not None:
        _, h1 = angle
        stages[-3].update({
            "h": str(angle[0]),
            "h1": str(h1),
            "h_finite": finite(angle[0]),
            "h1_finite": finite(h1),
        })
    else:
        h1 = None
    if hbar is not None and bracket is not None and h1 is not None:
        run_stage(
            "density",
            lambda: RADIUS * n * (-c * hbar + h1 * bracket) / (w * root),
            stages,
        )
    return {"stages": stages}


def build_kernels(mu: acb, w: acb):
    mu2 = mu * mu
    pair_jacobian = 1 - w - RADIUS * mu

    def outer_integrand(c, difference, rho2, n, analytic):
        p = difference * difference + mu2 * rho2
        s = rho2 + mu2 * c * c
        root = (p * s).sqrt(analytic=analytic)
        cosine = mu * n / root
        _, h1 = base.regular_angle_data(cosine)
        hbar = base.regular_hbar(cosine)
        bracket = -c + n * difference / p
        return n * (-c * hbar + h1 * bracket) / (2 * w * root)

    def inner(t, analytic):
        q = 8 * t - 4
        c = w + mu * q
        rho2 = (1 - c) * (1 + c)
        reduced_p = q * q + rho2
        s = rho2 + mu2 * c * c
        reduced_root = (reduced_p * s).sqrt(analytic=analytic)
        n = 1 - w * w - mu * w * q
        cosine = n / reduced_root
        _, h1 = base.regular_angle_data(cosine)
        hbar = base.regular_hbar(cosine)
        bracket = -c + n * q / (mu * reduced_p)
        return RADIUS * n * (-c * hbar + h1 * bracket) / (w * reduced_root)

    def paired(t, analytic):
        x = RADIUS * mu + pair_jacobian * t
        c_left = w - x
        rho2_left = (1 - w + x) * (1 + w - x)
        n_left = 1 - w * w + w * x
        left = outer_integrand(c_left, -x, rho2_left, n_left, analytic)
        c_right = w + x
        rho2_right = pair_jacobian * (1 - t) * (1 + w + x)
        n_right = 1 - w * w - w * x
        right = outer_integrand(c_right, x, rho2_right, n_right, analytic)
        return pair_jacobian * (left + right)

    def far_left(t, analytic):
        x = 1 - w + 2 * w * t
        c = w - x
        rho2 = 4 * w * (1 - t) * (1 - w + w * t)
        n = 1 + w - 2 * w * w + 2 * w * w * t
        return 2 * w * outer_integrand(c, -x, rho2, n, analytic)

    return inner, paired, far_left


def integral_record(kernel) -> dict:
    try:
        value = base.rigorous_integral(kernel, TOLERANCE, DEPTH, EVAL_LIMIT)
        return value_record(value)
    except Exception as exc:
        return {"finite": False, "error": f"{type(exc).__name__}: {exc}"}


def parameter_case(mu_q: fmpq, w_q: fmpq) -> dict:
    mu = point(mu_q)
    w = point(w_q)
    inner, paired, far_left = build_kernels(mu, w)
    t_values = [
        fmpq(0), fmpq(1, 16), fmpq(1, 4), fmpq(1, 2),
        fmpq(3, 4), fmpq(15, 16), fmpq(1),
    ]
    samples = []
    for t_q in t_values:
        t = point(t_q)
        for analytic in (False, True):
            pair_jacobian = 1 - w - RADIUS * mu
            x = RADIUS * mu + pair_jacobian * t
            c_left = w - x
            c_right = w + x
            samples.append({
                "t": str(t_q),
                "analytic": analytic,
                "inner": inner_detail(mu, w, t, analytic),
                "paired_left": outer_detail(
                    mu, w, c_left, -x,
                    (1 - w + x) * (1 + w - x),
                    1 - w * w + w * x,
                    analytic,
                ),
                "paired_right": outer_detail(
                    mu, w, c_right, x,
                    pair_jacobian * (1 - t) * (1 + w + x),
                    1 - w * w - w * x,
                    analytic,
                ),
                "far_left": outer_detail(
                    mu, w,
                    w - (1 - w + 2 * w * t),
                    -(1 - w + 2 * w * t),
                    4 * w * (1 - t) * (1 - w + w * t),
                    1 + w - 2 * w * w + 2 * w * w * t,
                    analytic,
                ),
                "kernel_values": {
                    "inner": value_record(inner(t, analytic)),
                    "paired": value_record(paired(t, analytic)),
                    "far_left": value_record(far_left(t, analytic)),
                },
            })
    integrals = {
        "inner": integral_record(inner),
        "paired_outer": integral_record(paired),
        "far_left": integral_record(far_left),
    }
    return {
        "mu": str(mu_q),
        "w": str(w_q),
        "samples": samples,
        "integrals": integrals,
    }


w_midpoints = [
    fmpq(7, 80), fmpq(3, 16), fmpq(5, 16),
    fmpq(7, 16), fmpq(9, 16), fmpq(11, 16),
]
result = {
    "status": "COMPLETE_DIAGNOSTIC",
    "purpose": "Locate the first non-finite stage in compact-tail H evaluation.",
    "mu": "1/160",
    "cases": [parameter_case(fmpq(1, 160), w) for w in w_midpoints],
}
output = Path("prolate_axis_tail_H_kernel_diagnostic.json")
output.write_text(json.dumps(result, indent=2), encoding="utf-8")
print(json.dumps(result, indent=2))
