"""
Driver for step 4-2: rigorous enclosure of B(lambda) at the four target
rationals, each split into 4 psi-bands, with per-band checkpointing to
RUN_LOG.md and a machine-readable enclosure_partial.json.

This driver only orchestrates the delivered, unmodified kernel
boundary_entry_independent.py (imported as M). It writes no mathematics
of its own; the integrand and angle data come entirely from M.

rel_tol starts at 2^-22 and is tightened (2^-26, 2^-30, 2^-34, 2^-38)
only when the sign of the 4-band total is still undecided. The sign is
NEVER read from a single band: only the sum of all four bands is judged.
"""

import sys
import time
import json
from flint import acb, arb, ctx

import boundary_entry_independent as M

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 30
ctx.dps = DPS
M._init_consts()

LOG = "RUN_LOG.md"
PARTIAL = "enclosure_partial.json"
BANDS = 4
EVAL_LIMIT = 40000
DEPTH_LIMIT = 30
RT_EXPS = [-22, -26, -30, -34, -38]

_partial = {"dps": DPS, "bands": BANDS, "eval_limit": EVAL_LIMIT,
            "depth_limit": DEPTH_LIMIT, "results": {}}


def log(msg):
    with open(LOG, "a") as f:
        f.write(msg + "\n")
    print(msg, flush=True)


def save_partial():
    with open(PARTIAL, "w") as f:
        json.dump(_partial, f, indent=2)


def inner(psi, lam, rel_tol):
    return acb.integral(
        lambda t, a: M.integrand(t, psi, lam, True),
        0, 1, rel_tol=rel_tol, eval_limit=EVAL_LIMIT, depth_limit=DEPTH_LIMIT)


def band_integral(lam, k, rel_tol):
    hi = acb.pi() / 2
    a = hi * k / BANDS
    b = hi * (k + 1) / BANDS
    return acb.integral(
        lambda psi, an: inner(psi, lam, rel_tol),
        a, b, rel_tol=rel_tol, eval_limit=EVAL_LIMIT, depth_limit=DEPTH_LIMIT)


def enclose(name, lam, want):
    lam = acb(lam)
    last = None
    for rt_exp in RT_EXPS:
        rel_tol = arb(2) ** rt_exp
        log("")
        log(f"#### {name}  want={want}  rel_tol=2^{rt_exp}  dps={DPS}")
        t_lam = time.time()
        total = acb(0)
        bands = []
        for k in range(BANDS):
            t0 = time.time()
            part = band_integral(lam, k, rel_tol)
            total += part
            bstr = part.real.str(15, radius=True)
            dt = time.time() - t0
            bands.append({"band": k + 1, "ball": bstr, "seconds": round(dt, 1)})
            log(f"- band {k + 1}/{BANDS}: {bstr}  ({dt:.1f}s)")
        re = total.real
        lo, up = arb(re.lower()), arb(re.upper())
        sign = "POSITIVE" if lo > 0 else ("NEGATIVE" if up < 0 else "UNDECIDED")
        log(f"- TOTAL: {re.str(15, radius=True)}  -> {sign}  "
            f"({time.time() - t_lam:.1f}s)")
        last = {"lambda": name, "want": want, "rel_tol_exp": rt_exp,
                "ball": re.str(30, radius=True),
                "lower": lo.str(25), "upper": up.str(25),
                "sign": sign, "band_detail": bands}
        _partial["results"][name] = last
        save_partial()
        if sign == want:
            return last
        log(f"  (sign still {sign}; tightening rel_tol)")
    return last


TARGETS = [
    ("B(103/50)",        acb(103) / 50,        "POSITIVE"),
    ("B(207/100)",       acb(207) / 100,       "NEGATIVE"),
    ("B(206538/100000)", acb(206538) / 100000, "POSITIVE"),
    ("B(206539/100000)", acb(206539) / 100000, "NEGATIVE"),
]

if __name__ == '__main__':
    log(f"\n## 4-2 enclosure run  ({time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())})")
    log(f"dps={DPS}, bands={BANDS}, eval_limit={EVAL_LIMIT}, "
        f"depth_limit={DEPTH_LIMIT}, rt_exps={RT_EXPS}")
    t_all = time.time()
    for name, lam, want in TARGETS:
        enclose(name, lam, want)
    log(f"\n## 4-2 done in {time.time() - t_all:.1f}s")
    _partial["total_seconds"] = round(time.time() - t_all, 1)
    save_partial()
    log("ALL_TARGETS_ATTEMPTED")
