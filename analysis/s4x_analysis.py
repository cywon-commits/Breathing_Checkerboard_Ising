#!/usr/bin/env python3
"""Expanded S4 (K = 0.1..0.3, h = 0, 0.5, eps = 0.03, 0.01): T_c(K), B-aligned fraction at T_c.
  python analysis/s4x_analysis.py -> results/s4x.json, results/figs/s4x_*.png
"""
import glob, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from lib import *
from analyze import crossing

D = os.path.join(ROOT, "data", "S4x")
Ks, Hs, Es = [0.1, 0.15, 0.2, 0.3], [0.0, 0.5], [0.03, 0.01]


def runs_for(K, h, eps):
    out = {}
    for Lb in (16, 32, 64):
        rs = find_runs(D, f"S4x_K{K:g}_h{h:g}_eps{eps:g}_Lb{Lb}")
        if rs: out[Lb] = combine_seeds(rs)
    return out


rows = []
for h in Hs:
    for eps in Es:
        for K in Ks:
            runs = runs_for(K, h, eps)
            Lbs = sorted(runs)
            if len(Lbs) < 2: continue
            nb = min(r["B"].shape[0] for r in runs.values())
            for r in runs.values():
                B = r["B"]; k = B.shape[0] // nb; r["B"] = B[: k * nb].reshape(nb, k, *B.shape[1:]).mean(axis=1)
            uf = binder_fn("N")
            U = {L: jackknife(runs[L]["B"], uf) for L in Lbs}
            cr = []
            for L1, L2 in zip(Lbs[:-1], Lbs[1:]):
                T1, T2 = runs[L1]["T"], runs[L2]["T"]
                x0 = crossing(T1, U[L1][0], T2, U[L2][0], lo=0.02, hi=0.66)
                xs = []
                for k in range(nb):
                    m1 = (runs[L1]["B"].sum(0) - runs[L1]["B"][k]) / (nb - 1); m2 = (runs[L2]["B"].sum(0) - runs[L2]["B"][k]) / (nb - 1)
                    xs.append(crossing(T1, uf(m1), T2, uf(m2), lo=0.02, hi=0.66))
                xs = np.array(xs)
                err = float(np.sqrt((nb - 1) / nb * np.nansum((xs - np.nanmean(xs)) ** 2))) if np.isfinite(x0) else np.nan
                cr.append((L1, L2, float(x0), err))
            good = [c for c in cr if np.isfinite(c[2])]
            Lbig = Lbs[-1]
            row = dict(K=K, h=h, eps=eps, Ek=min(8 * K, 4 * (1 - h)), region="B-mediated" if 8 * K < 4 * (1 - h) else "A-violation",
                       crossings=cr, minswap={L: float(min(runs[L]["meta"]["swap_acc"])) for L in Lbs})
            if good:
                Tc = good[-1][2]
                row.update(Tc=Tc, err=good[-1][3], drift=abs(good[-1][2] - good[-2][2]) if len(good) > 1 else None)
                T = runs[Lbig]["T"]; m = runs[Lbig]["B"].mean(0)
                for key in ("rBal", "rviol", "rdiag"):
                    row[key + "_Tc"] = float(CubicSpline(T, m[:, FI[key]])(Tc)) if T[0] < Tc < T[-1] else None
                row["U_Lbig"] = float(CubicSpline(T, U[Lbig][0])(Tc)) if T[0] < Tc < T[-1] else None
                row["Lbig"] = Lbig
            rows.append(row)
            print(f"h={h} eps={eps} K={K}: E_k={row['Ek']:.2f} {row['region']:12s}", f"Tc={row.get('Tc', float('nan')):.4f}+-{row.get('err', float('nan')):.4f}",
                  f"rBal(Tc)={row.get('rBal_Tc')}", "crossings", [(a, b, round(c, 3)) for a, b, c, _ in cr], "minswap", {k: round(v, 2) for k, v in row["minswap"].items()})
json.dump(rows, open(os.path.join(ROOT, "results", "s4x.json"), "w"), indent=1, default=float)

# figures: T_c vs K and rBal vs T for the largest size
fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
for h, mk in zip(Hs, ("o", "s")):
    for eps, c in zip(Es, ("C0", "C1")):
        R = [r for r in rows if r["h"] == h and r["eps"] == eps and "Tc" in r]
        axs[0].errorbar([r["K"] for r in R], [r["Tc"] for r in R], [max(r["err"], 1e-3) for r in R], fmt=mk + "-", color=c, mfc="none" if h else c,
                        label=f"h={h:g}, ε={eps:g}")
axs[0].set_xlabel("K"); axs[0].set_ylabel("T_c (Binder crossing, largest pair)"); axs[0].legend(fontsize=7); axs[0].set_title("T_c(K)", fontsize=9)
for h, mk in zip(Hs, ("o", "s")):
    for eps, c in zip(Es, ("C0", "C1")):
        R = [r for r in rows if r["h"] == h and r["eps"] == eps and r.get("rBal_Tc") is not None]
        axs[1].plot([r["K"] for r in R], [r["rBal_Tc"] for r in R], mk + "-", color=c, mfc="none" if h else c, label=f"h={h:g}, ε={eps:g}")
axs[1].set_xlabel("K"); axs[1].set_ylabel("B-aligned fraction n_B∈{0,4} at T_c"); axs[1].legend(fontsize=7); axs[1].set_title("superspin picture validity", fontsize=9)
for K, c in zip(Ks, ("C0", "C1", "C2", "C3")):
    for h, ls in zip(Hs, ("-", "--")):
        rs = runs_for(K, h, 0.03)
        if not rs: continue
        L = max(rs); r = rs[L]
        axs[2].plot(r["T"] / (rows and [x for x in rows if x["K"] == K and x["h"] == h and x["eps"] == 0.03 and "Tc" in x][0]["Tc"]), r["B"].mean(0)[:, FI["rBal"]], ls, color=c, label=f"K={K:g}, h={h:g}")
axs[2].set_xlabel("T / T_c"); axs[2].set_ylabel("B-aligned fraction"); axs[2].legend(fontsize=6, ncol=2); axs[2].set_title("ε=0.03, largest size", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "results", "figs", "s4x.png"), dpi=130)
