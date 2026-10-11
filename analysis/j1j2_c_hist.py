#!/usr/bin/env python3
"""eps<0 (g = J2/J1 = 0.518, 0.538): energy histograms near the stripe transition, looking for a double peak
(weak first-order transition, revision 8a P5).  Uses the time series stored by j1j2 (series_every=1).
  python analysis/j1j2_c_hist.py -> results/c_hist.json, results/figs/c_hist.png
"""
import glob, json, os, re, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "j1j2")


def load_series(prefix):
    meta = json.load(open(prefix + ".meta.json"))
    nrep = meta["nrep"]
    a = np.fromfile(prefix + ".series.f32", dtype=np.float32).reshape(-1, nrep, 3)
    return meta, a


def peaks(e, nb=120):
    h, edges = np.histogram(e, bins=nb)
    c = 0.5 * (edges[1:] + edges[:-1])
    hs = gaussian_filter1d(h.astype(float), 2.0)
    loc = [i for i in range(1, nb - 1) if hs[i] > hs[i - 1] and hs[i] >= hs[i + 1] and hs[i] > 0.05 * hs.max()]
    return c, hs / hs.sum() / (c[1] - c[0]), loc, hs


out = {}
fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
for ax, eps in zip(axs, (-0.012, -0.025)):
    for Lb, col in zip((16, 32, 64, 128), ("C0", "C1", "C2", "C3")):
        fs = sorted(glob.glob(os.path.join(D, f"j1j2_C_h0_eps{eps:g}_Lb{Lb}_s1.meta.json")))
        if not fs: continue
        meta, a = load_series(fs[0][:-len(".meta.json")])
        e = a[:, :, 0].astype(float)         # (nrec, nrep)
        ms2 = a[:, :, 2].astype(float)
        T = np.array(meta["T"]); N = Lb * Lb
        var = e.var(axis=0) * N / T ** 2      # C_v per spin
        it = int(np.argmax(var))
        c, p, loc, hs = peaks(e[:, it])
        rec = dict(Lb=Lb, T_Cmax=float(T[it]), Cmax=float(var[it]), npeaks=len(loc))
        if len(loc) >= 2:
            i1, i2 = sorted(sorted(loc, key=lambda i: -hs[i])[:2])
            valley = hs[i1:i2 + 1].min()
            rec.update(peak_sep=float(c[i2] - c[i1]), valley_over_lowerpeak=float(valley / min(hs[i1], hs[i2])))
        # energy Binder cumulant V = 1 - <e^4>/(3<e^2>^2) of centred e at the C_max temperature
        d = e[:, it] - e[:, it].mean()
        rec["V"] = float(1 - (d ** 4).mean() / (3 * (d ** 2).mean() ** 2))
        out[f"eps{eps:g}_Lb{Lb}"] = rec
        ax.plot(c - e[:, it].mean(), p, color=col, label=f"L_b={Lb}, T={T[it]:.3f}")
        print(f"eps={eps} Lb={Lb} T(Cmax)={T[it]:.4f} Cmax={var[it]:.2f}", {k: round(v, 4) if isinstance(v, float) else v for k, v in rec.items() if k not in ("Lb", "T_Cmax", "Cmax")})
    ax.set_title(f"ε={eps:g} (g={(1-2*eps)/(2*(1+eps)):.3f}): P(e) at the C_v maximum", fontsize=9); ax.set_xlabel("e − <e>"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "results", "figs", "c_hist.png"), dpi=130)
json.dump(out, open(os.path.join(ROOT, "results", "c_hist.json"), "w"), indent=1)
