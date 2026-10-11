#!/usr/bin/env python3
"""B-superspin snapshots at 0.7, 1.0, 1.3 x T_c (S1: K=4, h=0, eps=0.025, Lb=64, T_c=1.2)."""
import os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Lb = 64; L = 2 * Lb; Tc = 1.2
raw = np.fromfile(os.path.join(ROOT, "data/snapshots/snap_S1_eps0.025_Lb64.snap.i8"), dtype=np.int8).reshape(3, L, L)
fig, axs = plt.subplots(1, 3, figsize=(11, 3.8))
for ax, s, f in zip(axs, raw, (0.7, 1.0, 1.3)):
    tau = np.zeros((Lb, Lb))
    for Y in range(Lb):
        for X in range(Lb):   # B cell (X,Y): sites (2X+1..2X+2, 2Y+1..2Y+2)
            tau[Y, X] = sum(s[(2 * Y + 1 + dy) % L, (2 * X + 1 + dx) % L] for dx in (0, 1) for dy in (0, 1)) / 4
    ax.imshow(tau, cmap="RdBu", vmin=-1, vmax=1, origin="lower", interpolation="nearest")
    ax.set_title(f"T = {f} T_c = {f*Tc:.2f}", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
fig.suptitle("B superspin τ = S_B/4, S1 (K=4, h=0, ε=0.025), L=128", fontsize=10)
fig.tight_layout()
os.makedirs(os.path.join(ROOT, "results/figs"), exist_ok=True)
fig.savefig(os.path.join(ROOT, "results/figs/snapshots.png"), dpi=130)
