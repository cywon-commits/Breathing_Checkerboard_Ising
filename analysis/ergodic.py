#!/usr/bin/env python3
"""Ergodicity check: same PT run from different initial conditions (Neel, X-stripe, random), with and without
strip flips (moves=7 removes them).  Reports max |ΔU|/σ and max |Δ<mN^2>| between initial conditions.
Point: K=4, h=0, eps=0.05 (and 0.012), Lb=16."""
import os, subprocess, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from lib import *

def run(eps, init, moves, seed, Tmin, Tmax):
    out = os.path.join(ROOT, "data", "ergodic", f"erg_eps{eps:g}_{init}_m{moves}_s{seed}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if not os.path.exists(out + ".meta.json"):
        subprocess.run([os.path.join(ROOT, "build", "bcising"), "--Lb", "16", "--K", "4", "--h", "0", "--eps", str(eps),
            "--Tmin", str(Tmin), "--Tmax", str(Tmax), "--nrep", "24", "--neq", "3000", "--nmeas", "20000", "--nblocks", "20",
            "--seed", str(seed), "--init", init, "--moves", str(moves), "--out", out], check=True,
            env=dict(os.environ, OMP_NUM_THREADS="1"))
    return load_run(out)

if __name__ == "__main__":
    for eps, (Tmin, Tmax) in {0.05: (1.3, 1.9), 0.012: (0.8, 1.15)}.items():
        res = {}
        for moves, mname in [(15, "all moves"), (7, "no strip flips")]:
            for init in ["neel", "stripe_x", "random"]:
                res[(moves, init)] = run(eps, init, moves, 7, Tmin, Tmax)
        for moves, mname in [(15, "all moves"), (7, "no strip flips")]:
            ref = res[(moves, "neel")]
            Ur, dUr = jackknife(ref["B"], binder_fn("N"))
            Mr = ref["B"].mean(0)[:, FI["mN2"]]
            for init in ["stripe_x", "random"]:
                r = res[(moves, init)]
                U, dU = jackknife(r["B"], binder_fn("N"))
                M = r["B"].mean(0)[:, FI["mN2"]]
                z = np.abs(U - Ur) / np.sqrt(dU ** 2 + dUr ** 2 + 1e-12)
                print(f"eps={eps} {mname:15s} {init:9s} vs neel: max|dU|={np.abs(U-Ur).max():.3f}  max z={z.max():.1f}  "
                      f"max|d<mN^2>|={np.abs(M-Mr).max():.3f}  (lowest-T mN2 {M[0]:.3f} vs {Mr[0]:.3f})")
