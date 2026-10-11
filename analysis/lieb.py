#!/usr/bin/env python3
"""Lieb-ice limit check: eps = 0, J_B > 0 (K = -J_B < 0), J_A = 1.

Every A and B plaquette is then an antiferromagnetic 4-spin cluster and the ground manifold is the six-vertex
square-ice ensemble: S_0 = ln(W)/2 per spin with W = (4/3)^{3/2}  ->  0.2158.
The entropy is obtained by thermodynamic integration of the PT energy data:
    S(T0) = ln 2 - tail - int_{T0}^{Tmax} C(T)/T dT .
usage: python analysis/lieb.py [Lb ...]   (runs the simulations if the data are missing)
"""
import os, sys, subprocess, json
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from lib import load_run, jackknife, FI

S_LIEB = 0.5 * np.log((4 / 3) ** 1.5)


def run(Lb, JB=4.0, nrep=40, Tmin=0.2, Tmax=8.0, neq=5000, nmeas=40000, seed=1, init="random"):
    out = os.path.join(ROOT, "data", "lieb", f"lieb_JB{JB:g}_Lb{Lb}_{init}_s{seed}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    if not os.path.exists(out + ".meta.json"):
        subprocess.run([os.path.join(ROOT, "build", "bcising"), "--Lb", str(Lb), "--K", str(-JB), "--h", "0", "--eps", "0",
                        "--Tmin", str(Tmin), "--Tmax", str(Tmax), "--nrep", str(nrep), "--neq", str(neq), "--nmeas", str(nmeas),
                        "--nblocks", "20", "--seed", str(seed), "--init", init, "--out", out], check=True,
                       env=dict(os.environ, OMP_NUM_THREADS="1"))
    return load_run(out)


def entropy(r):
    T, N = r["T"], 4 * r["meta"]["Lb"] ** 2
    def cv(m): return N * (m[:, FI["e2"]] - m[:, FI["e"]] ** 2) / T ** 2
    def ent(m):
        C = cv(m)
        lnT = np.log(T)
        # cumulative trapezoid of C dlnT from the top
        integ = np.concatenate([[0], np.cumsum(0.5 * (C[1:] + C[:-1]) * np.diff(lnT))])
        tail = C[-1] * T[-1] ** 2 / (2 * T[-1] ** 2)   # C ~ c2/T^2 for T > Tmax
        return np.log(2) - tail - (integ[-1] - integ)
    return jackknife(r["B"], ent)


if __name__ == "__main__":
    Lbs = [int(x) for x in sys.argv[1:]] or [4, 8]
    for Lb in Lbs:
        for init in ["random", "neel"]:
            r = run(Lb, init=init)
            S, dS = entropy(r)
            rv, _ = jackknife(r["B"], lambda m: m[:, FI["rviol"]])
            print(f"Lb={Lb} init={init} swap_acc_min={min(r['meta']['swap_acc']):.2f}")
            for i in range(0, len(r["T"]), 5):
                print(f"   T={r['T'][i]:.3f}  S={S[i]:.4f} +- {dS[i]:.4f}   A-violation density {rv[i]:.2e}")
            print(f"   lowest-T entropy {S[0]:.4f} +- {dS[0]:.4f}   (Lieb: {S_LIEB:.4f})")
