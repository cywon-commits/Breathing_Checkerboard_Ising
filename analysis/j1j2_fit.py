#!/usr/bin/env python3
"""J1-J2 results: comparison with breathing-model S1/S3 at identical lattice sizes, T_c/J1, slope fits.
  python analysis/j1j2_fit.py  (needs results/j1j2.json from j1j2_analysis.py and results/fit.json)
"""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = json.load(open(os.path.join(ROOT, "results", "j1j2.json")))
F = json.load(open(os.path.join(ROOT, "results", "fit.json")))


def cross(rec, pair):
    for c in rec["crossings"]:
        if (c["L1"], c["L2"]) == pair and np.isfinite(c["Tx"]):
            return c["Tx"], c["err"]
    return None


print("T_c at identical lattice sizes: breathing model (K=4,8) vs J1-J2 model")
print("| set | ε | pair (B-lattice sizes) | breathing K=4 | breathing K=8 | J1-J2 |")
rows = []
for eps in (0.003, 0.006, 0.012, 0.025, 0.05):
    j = J.get(f"A_eps{eps:g}")
    for pair in ((16, 32), (32, 64)):
        s1 = next((r for r in F["rows"] if r["set"] == "S1" and abs(r["eps"] - eps) < 1e-9), None)
        s2 = next((r for r in F["rows"] if r["set"] == "S2" and abs(r["eps"] - eps) < 1e-9), None)
        v1 = next((c for c in s1["all_crossings"] if tuple(c[:2]) == pair), None) if s1 and s1["Tc"] == s1["Tc"] else None
        v2 = next((c for c in s2["all_crossings"] if tuple(c[:2]) == pair), None) if s2 and s2["Tc"] == s2["Tc"] else None
        vj = cross(j, pair) if j else None
        f = lambda v: "–" if v is None else (f"{v[2]:.3f} ± {v[3]:.3f}" if len(v) == 4 else f"{v[0]:.3f} ± {v[1]:.3f}")
        print(f"| S1/S2 | {eps:g} | {pair} | {f(v1)} | {f(v2)} | {f(vj)} |")
for eps in (0.012, 0.025):
    j = J.get(f"B_eps{eps:g}")
    s3 = next((r for r in F["rows"] if r["set"] == "S3" and abs(r["eps"] - eps) < 1e-9), None)
    for pair in ((16, 32), (32, 64)):
        v3 = next((c for c in s3["all_crossings"] if tuple(c[:2]) == pair), None) if s3 else None
        vj = cross(j, pair) if j else None
        f = lambda v: "–" if v is None else (f"{v[2]:.3f} ± {v[3]:.3f}" if len(v) == 4 else f"{v[0]:.3f} ± {v[1]:.3f}")
        print(f"| S3 (h=0.5) | {eps:g} | {pair} | {f(v3)} | | {f(vj)} |")

# largest-size T_c(eps) of the J1-J2 model (h=0) and slope fit
pts = []
for eps in (0.05, 0.025, 0.012, 0.006):
    j = J[f"A_eps{eps:g}"]
    c = [c for c in j["crossings"] if np.isfinite(c["Tx"])]
    pts.append((eps, c[-1]["Tx"], c[-1]["err"], (c[-1]["L1"], c[-1]["L2"]), c[-1]["Tx"] - c[-2]["Tx"] if len(c) > 1 else np.nan))
print("\nJ1-J2 (h=0), largest-pair crossing:")
for p in pts:
    g = (1 - 2 * p[0]) / (2 * (1 + p[0]))
    print(f"  eps={p[0]:g} g={g:.4f} T_c={p[1]:.3f}±{p[2]:.3f} pair={p[3]} last-step={p[4]:+.3f}  T_c/J1={p[1]/(2*(1+p[0])):.4f}")
e = np.array([p[0] for p in pts]); T = np.array([p[1] for p in pts]); x = np.log(1 / e); y = 1 / T
loc = [(y[i + 1] - y[i]) / (x[i + 1] - x[i]) for i in range(len(x) - 1)]
print("local slopes (eps decreasing):", [round(v, 3) for v in loc])
A = np.polyfit(x, y, 1)
print("straight-line fit over eps=0.006..0.05: slope %.3f (pred 0.25), intercept %.3f (pred ln(4/6)/4 = %.3f)" % (A[0], A[1], np.log(4 / 6) / 4))
fig, ax = plt.subplots(figsize=(5.2, 4))
ax.errorbar(x, y, [p[2] / p[1] ** 2 for p in pts], fmt="o", label="J1-J2 (largest pair)")
s1 = [r for r in F["rows"] if r["set"] == "S1" and r["Tc"] == r["Tc"]]
ax.plot([np.log(1 / r["eps"]) for r in s1], [1 / r["Tc"] for r in s1], "x", color="C1", label="breathing S1 (K=4)")
xx = np.linspace(2.8, 6.0, 5); ax.plot(xx, (xx + np.log(4 / 6)) / 4, "k--", label="prediction (slope 1/4)")
ax.set_xlabel("ln(1/ε)"); ax.set_ylabel("1/T_c"); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(ROOT, "results", "figs", "j1j2_slopes.png"), dpi=130)
