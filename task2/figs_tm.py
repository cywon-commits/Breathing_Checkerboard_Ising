#!/usr/bin/env python3
"""Summary figure of the task-2 results -> results/summary.png"""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
R = {}
for l in open(os.path.join(HERE, "data", "results.jsonl")):
    r = json.loads(l); R[(r["model"], r["W"], round(r["u"], 6), r["M"])] = r
Ws = sorted({k[1] for k in R})
def eig(m, W, u, M=0):
    r = R.get((m, W, round(u, 6), M)); return None if r is None else np.abs(r["eig_re"])
def x2(m, W, u):
    a, b = eig(m, W, u, 0), eig(m, W, u, 2)
    return None if a is None or b is None else W / (4 * np.pi) * np.log(a[0] / b[0])
def xi(m, W, u):
    a = eig(m, W, u, 0)
    return None if a is None or len(a) < 2 or a[1] >= a[0] else -1 / np.log(a[1] / a[0])
fig, ax = plt.subplots(2, 2, figsize=(11, 8))
us = sorted({k[2] for k in R if k[0] == "stag" and k[3] == 2 and k[1] == 24 and k[2] <= 2.4})
for W, c in zip(Ws, plt.cm.viridis(np.linspace(0, 0.9, len(Ws)))):
    uu = [u for u in us if x2("stag", W, u) is not None]
    ax[0, 0].plot(uu, [x2("stag", W, u) for u in uu], "o-", ms=3, color=c, label=f"W={W}")
ax[0, 0].axhline(0.5, color="k", ls=":"); ax[0, 0].axvline(1.40, color="r", ls=":", lw=0.8)
ax[0, 0].set_xlabel("u"); ax[0, 0].set_ylabel("x2(W,u)"); ax[0, 0].legend(fontsize=7); ax[0, 0].set_title("x2 = (W/4π) ln(Λ0/Λ2): W-independent for 1.4<u<~2", fontsize=9)
ulist = [1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.0, 4.0, 6.0]
for u, c in zip(ulist, plt.cm.plasma(np.linspace(0, 0.9, len(ulist)))):
    pts = [(W, xi("stag", W, u)) for W in Ws if xi("stag", W, u) is not None]
    ax[0, 1].plot([p[0] for p in pts], [p[1] / p[0] for p in pts], "o-", ms=3, color=c, label=f"u={u:g}")
ax[0, 1].set_xlabel("W"); ax[0, 1].set_ylabel("ξ_y/W"); ax[0, 1].legend(fontsize=7, ncol=2); ax[0, 1].set_title("ξ_y/W: flat at u=1.4, rising for u>1.4", fontsize=9)
for u, c in zip((1.6, 2.0, 2.4, 3.0), ("C0", "C1", "C2", "C3")):
    pts = [(W, eig("stag", W, u)) for W in Ws if eig("stag", W, u) is not None and len(eig("stag", W, u)) >= 3]
    ax[1, 0].plot([p[0] for p in pts], [p[0] * np.log(p[1][0] / p[1][1]) for p in pts], "o-", color=c, label=f"u={u:g}: W·g1")
    ax[1, 0].plot([p[0] for p in pts], [p[0] * np.log(p[1][0] / p[1][2]) for p in pts], "s--", color=c, ms=3, label=f"W·g2")
ax[1, 0].set_xlabel("W"); ax[1, 0].set_ylabel("W ln(λ0/λ_i)"); ax[1, 0].legend(fontsize=6, ncol=2); ax[1, 0].set_title("low M=0 spectrum: λ1 closes faster than 1/W", fontsize=9)
for u, c in zip((3.0, 4.0, 5.0, 6.0), ("C0", "C1", "C2", "C3")):
    pts = [(W, eig("stag", W, u)) for W in Ws if eig("stag", W, u) is not None and len(eig("stag", W, u)) >= 2]
    ax[1, 1].plot([p[0] for p in pts], [(p[1][1] / p[1][0] - 0.5) / (p[0] * np.exp(-u)) for p in pts], "o-", color=c, label=f"u={u:g}")
ax[1, 1].set_xlabel("W"); ax[1, 1].set_ylabel("(|λ1/λ0| − 1/2) / (W e^{-u})"); ax[1, 1].legend(fontsize=7); ax[1, 1].set_title("kink-gas signal: → 0.25", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "results", "summary.png"), dpi=130)
