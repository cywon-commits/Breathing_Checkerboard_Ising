#!/usr/bin/env python3
"""Figures for REPORT_task3: results/figs/ice_summary.png (4 panels) and ice_PQ.png (P(Q) histograms from the PT run)."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "results", "ice_summary.json"))); u = np.array(d["u"]); Ls = d["Ls"]
cols = plt.cm.viridis(np.linspace(0, 0.9, len(Ls)))
strip = {1.2: 0.4091, 1.3: 0.4483, 1.4: 0.4929, 1.5: 0.5413, 1.6: 0.5975, 1.7: 0.6567, 1.8: 0.7272}   # task-2 x2(W=28)
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
for L, c in zip(Ls, cols):
    ax[0, 0].plot(u, d["UQ"][str(L)], "o-", ms=3, color=c, label=f"L={L}")
    ax[0, 1].plot(u, d["Q2"][str(L)], "o-", ms=3, color=c)
    g = np.array(d["geff"][str(L)], dtype=float); ax[1, 0].plot(u, g, "o-", ms=3, color=c)
    ax[1, 1].plot(u, np.array(d["xi_par"][str(L)]) / (np.exp(u) / 2), "o-", ms=3, color=c)
ax[1, 0].plot(list(strip), [2 * v for v in strip.values()], "k*", ms=9, label="task 2 strip, 2·x2(W=28)")
ax[0, 0].set_title("Binder cumulant U_Q of the nematic parameter"); ax[0, 0].legend(fontsize=7)
ax[0, 1].set_title("<Q²>"); ax[0, 1].set_yscale("log")
ax[1, 0].set_title("g_eff from <W²> (P(W)∝exp(−π g W²)); g=1 at u_c"); ax[1, 0].set_ylim(0.4, 3.2); ax[1, 0].axhline(1, color="k", ls=":"); ax[1, 0].legend(fontsize=7)
ax[1, 1].set_title("ξ_par / (e^u/2)  (second-moment length along the Q axis)"); ax[1, 1].set_ylim(0.8, 2.2); ax[1, 1].axhline(1, color="k", ls=":")
for a in ax.flat: a.set_xlabel("u"); a.axvline(1.38, color="r", ls=":", lw=0.8)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "results", "figs", "ice_summary.png"), dpi=130)
# P(Q) from the parallel-tempering run (both families visited)
import glob
meta = json.load(open(os.path.join(HERE, "data/pt/ptfull_L32_s1.meta.json"))); uu = meta["u"]
a = np.fromfile(os.path.join(HERE, "data/pt/ptfull_L32_s1.series.f32"), dtype=np.float32).reshape(-1, len(uu), 4)
fig, ax = plt.subplots(1, 5, figsize=(15, 3))
for a_, uv in zip(ax, (1.0, 1.4, 1.8, 2.2, 3.0)):
    Q = a[:, uu.index(uv), 0]; a_.hist(Q, bins=60, range=(-0.4, 0.4), density=True, color="C0"); a_.set_title(f"L=32, u={uv}", fontsize=9); a_.set_xlabel("Q")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "results", "figs", "ice_PQ.png"), dpi=120)
