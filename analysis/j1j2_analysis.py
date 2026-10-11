#!/usr/bin/env python3
"""Analysis of the J1-J2 cross-check runs (data/j1j2): Binder crossings, T_c(eps), slope fits,
comparison with the breathing-model S1/S3 results, and the energy histograms for eps<0 (group C).
  python analysis/j1j2_analysis.py     -> results/j1j2.json, results/figs/j1j2_*.png
"""
import glob, json, os, re, sys
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from lib import *
from analyze import crossing

D = os.path.join(ROOT, "data", "j1j2")
LO = float(os.environ.get("BINDER_LO", "-0.6"))   # lower U bound for accepted crossings


def load_group(g, h, eps, Lb):
    runs = [load_run(f[:-len(".meta.json")]) for f in sorted(glob.glob(os.path.join(D, f"j1j2_{g}_h{h:g}_eps{eps:g}_Lb{Lb}_s*.meta.json")))]
    return combine_seeds(runs) if runs else None


def analyse(g, h, eps):
    Lbs = sorted({int(re.search(r"_Lb(\d+)_", f).group(1)) for f in glob.glob(os.path.join(D, f"j1j2_{g}_h{h:g}_eps{eps:g}_Lb*_s*.meta.json"))})
    runs = {L: load_group(g, h, eps, L) for L in Lbs}
    runs = {L: r for L, r in runs.items() if r is not None}
    if len(runs) < 2: return None
    kind = "N" if eps > 0 else "S"
    uf = binder_fn(kind)
    nb = min(r["B"].shape[0] for r in runs.values())
    for L, r in runs.items():
        B = r["B"]; k = B.shape[0] // nb; r["B"] = B[: k * nb].reshape(nb, k, *B.shape[1:]).mean(axis=1)
    U = {L: jackknife(r["B"], uf) for L, r in runs.items()}
    res = dict(Lb=sorted(runs), crossings=[], minswap={L: float(min(r["meta"]["swap_acc"])) for L, r in runs.items()})
    Ls = sorted(runs)
    for L1, L2 in zip(Ls[:-1], Ls[1:]):
        T1, T2 = runs[L1]["T"], runs[L2]["T"]
        x0 = crossing(T1, U[L1][0], T2, U[L2][0], lo=LO)
        xs = []
        for k in range(nb):
            m1 = (runs[L1]["B"].sum(0) - runs[L1]["B"][k]) / (nb - 1); m2 = (runs[L2]["B"].sum(0) - runs[L2]["B"][k]) / (nb - 1)
            xs.append(crossing(T1, uf(m1), T2, uf(m2), lo=LO))
        xs = np.array(xs)
        err = float(np.sqrt((nb - 1) / nb * np.nansum((xs - np.nanmean(xs)) ** 2))) if np.isfinite(x0) else np.nan
        res["crossings"].append(dict(L1=L1, L2=L2, Tx=float(x0), err=err))
    return res, runs, U


if __name__ == "__main__":
    out = {}
    for g, h, epss in [("A", 0.0, [0.003, 0.006, 0.012, 0.025, 0.05]), ("B", 0.5, [0.003, 0.006, 0.012, 0.025]), ("C", 0.0, [-0.012, -0.025]), ("A3", 0.0, [0.003])]:
        for eps in epss:
            r = analyse(g, h, eps)
            if r is None: continue
            res, runs, U = r
            out[f"{g}_eps{eps:g}"] = dict(h=h, eps=eps, g=(1 - 2 * eps) / (2 * (1 + eps)), **res)
            print(g, f"eps={eps:+g} g=J2/J1={(1-2*eps)/(2*(1+eps)):.4f}", [f"({c['L1']},{c['L2']}):{c['Tx']:.4f}+-{c['err']:.4f}" for c in res["crossings"]], "minswap", {k: round(v, 2) for k, v in res["minswap"].items()})
    json.dump(out, open(os.path.join(ROOT, "results", "j1j2.json"), "w"), indent=1)
