#!/usr/bin/env python3
"""Analysis of the task-2 transfer-matrix data:  calibrations, u_c (KT), xi_y/W, large-u liquid checks.
   python analyze_tm.py    -> results/*.png, results/tm_summary.json, printed tables
"""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
R = {}
for line in open(os.path.join(HERE, "data", "results.jsonl")):
    r = json.loads(line)
    R[(r["model"], r["W"], round(r["u"], 6), r["M"])] = r
Ws = sorted({k[1] for k in R})


def lam(model, W, u, M, i=0):
    r = R.get((model, W, round(u, 6), M))
    if r is None or len(r["eig_re"]) <= i: return None
    return r["eig_re"][i]   # imaginary parts are checked separately


def x2(model, W, u):
    l0, l2 = lam(model, W, u, 0), lam(model, W, u, 2)
    return None if l0 is None or l2 is None else W / (4 * np.pi) * np.log(abs(l0) / abs(l2))


def xi_y(model, W, u):
    l0, l1 = lam(model, W, u, 0, 0), lam(model, W, u, 0, 1)
    if l0 is None or l1 is None or abs(l1) >= abs(l0): return None
    return -1.0 / np.log(abs(l1 / l0))


def entropy(model, W, u):
    l0 = lam(model, W, u, 0)
    return None if l0 is None else (np.log(abs(l0)) - 0.0) / (2 * W)   # per spin, without the e^u weights (u=0 only)


out = dict(Ws=Ws)
maxim = max((abs(v) for r in R.values() for v in r["eig_im"]), default=0.0)
maxleak = max(r["leak"] for r in R.values())
print(f"widths {Ws}; largest |Im lambda| among stored eigenvalues {maxim:.2e}; largest sector leakage {maxleak:.1e}")

# ------------------------------------------------------------ calibrations (C4)
print("\n== Calibration ==")
cal = {}
for name, model, u, target in [("Lieb ice u=0", "F", 0.0, "1/6 = 0.1667 ; S = 0.2158"), ("F model u=0.5", "F", 0.5, "g/2 = 0.308"),
                               ("F model u=ln2 (KT)", "F", np.log(2), "-> 1/2 with 1/ln W correction")]:
    xs = [(W, x2(model, W, u)) for W in Ws if x2(model, W, u) is not None]
    row = dict(target=target, x2={W: v for W, v in xs})
    if u == 0.0:
        row["entropy"] = {W: entropy(model, W, u) for W in Ws if entropy(model, W, u) is not None}
    cal[name] = row
    print(f"{name:20s} target {target}:  " + "  ".join(f"W={W}: {v:.4f}" for W, v in xs))
    if u == 0.0: print("   entropy/spin: " + "  ".join(f"W={W}: {v:.5f}" for W, v in row["entropy"].items()))
# KT log correction at the F-model KT point:  1/(1/2 - x2) = (ln W - ln W0)/c
xs = [(W, x2("F", W, np.log(2))) for W in Ws if x2("F", W, np.log(2)) is not None and W >= 12]
if len(xs) >= 3:
    lw = np.log([a for a, _ in xs]); y = 1 / (0.5 - np.array([b for _, b in xs]))
    p = np.polyfit(lw, y, 1)
    cF = 1 / p[0]; W0F = np.exp(-p[1] / p[0])
    res = y - np.polyval(p, lw)
    print(f"   KT fit  x2 = 1/2 - c/ln(W/W0) over W={[a for a,_ in xs]}:  c = {cF:.3f}, W0 = {W0F:.3f}, max residual of 1/(1/2-x2) = {abs(res).max():.3f}")
    cal["KT_fit_F"] = dict(c=cF, W0=W0F, Ws=[a for a, _ in xs])
out["calibration"] = cal

# ------------------------------------------------------------ u_c
print("\n== x2(W,u) and the crossing x2 = 1/2 ==")
ug = sorted({k[2] for k in R if k[0] == "stag" and k[3] == 2 and k[1] == Ws[-1]})
fig, axs = plt.subplots(1, 3, figsize=(14, 4))
ustar = {}
for model, ugrid in (("stag", [u for u in ug if u <= 2.4]), ("F", sorted({k[2] for k in R if k[0] == "F" and k[1] == Ws[-1]}))):
    ustar[model] = {}
    for W in Ws:
        us = [u for u in ugrid if x2(model, W, u) is not None]
        if len(us) < 4: continue
        xv = np.array([x2(model, W, u) for u in us])
        ax = axs[0] if model == "stag" else axs[1]
        ax.plot(us, xv, "o-", ms=3, label=f"W={W}")
        cs = CubicSpline(us, xv - 0.5)
        uu = np.linspace(us[0], us[-1], 2000); vv = cs(uu)
        sgn = np.nonzero(vv[:-1] * vv[1:] < 0)[0]
        if len(sgn): ustar[model][W] = float(brentq(lambda t: cs(t), uu[sgn[0]], uu[sgn[0] + 1]))
    ax = axs[0] if model == "stag" else axs[1]
    ax.axhline(0.5, color="k", ls=":"); ax.set_xlabel("u"); ax.set_ylabel("x2(W,u)")
    ax.set_title("staggered (A: e^u, B: 1)" if model == "stag" else "uniform F model (KT at u=ln2)", fontsize=9); ax.legend(fontsize=6)
    print(model, "u*(W) with x2=1/2:", {W: round(v, 4) for W, v in ustar[model].items()})
out["ustar"] = ustar
if "KT_fit_F" in cal:
    cF, W0F = cal["KT_fit_F"]["c"], cal["KT_fit_F"]["W0"]
    # calibrated threshold: at the KT point  x2(W) = 1/2 - c/ln(W/W0)  (universal c); u_c(W): x2^stag(W,u)= threshold(W)
    ucal = {}
    for W in Ws:
        us = [u for u in ug if u <= 2.4 and x2("stag", W, u) is not None]
        if len(us) < 4: continue
        thr = 0.5 - 0.5 + (x2("F", W, np.log(2)) if x2("F", W, np.log(2)) is not None else np.nan)   # F-model value at its KT point
        cs = CubicSpline(us, [x2("stag", W, u) - thr for u in us])
        uu = np.linspace(us[0], us[-1], 2000); vv = cs(uu); sg = np.nonzero(vv[:-1] * vv[1:] < 0)[0]
        if len(sg): ucal[W] = float(brentq(lambda t: cs(t), uu[sg[0]], uu[sg[0] + 1]))
    print("calibrated-threshold estimate u_c(W) (threshold = F-model x2 at its KT point u=ln2):", {W: round(v, 4) for W, v in ucal.items()})
    out["u_calibrated"] = ucal
axs[2].set_axis_off()
fig.tight_layout(); fig.savefig(os.path.join(HERE, "results", "x2.png"), dpi=130)

# ------------------------------------------------------------ xi_y / W
print("\n== xi_y(W)/W, M=0 sector (|lambda1/lambda0|) ==")
ugrid = [u for u in sorted({k[2] for k in R if k[0] == "stag" and k[3] == 0}) if 1.5 <= u <= 6.01]
fig2, axs2 = plt.subplots(1, 2, figsize=(11, 4))
tab = {}
for u in ugrid:
    vals = {W: xi_y("stag", W, u) for W in Ws}
    vals = {W: v for W, v in vals.items() if v is not None}
    tab[u] = vals
    if len(vals) >= 3:
        axs2[0].plot(list(vals), [v / W for W, v in vals.items()], "o-", ms=3, label=f"u={u:g}")
        axs2[1].plot(list(vals), list(vals.values()), "o-", ms=3, label=f"u={u:g}")
for u in ugrid:
    if u in (1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 4.0, 5.0, 6.0):
        print(f"u={u:4.1f}: " + "  ".join(f"W={W}: {v:.3f} ({v/W:.3f})" for W, v in tab[u].items()))
axs2[0].set_xlabel("W"); axs2[0].set_ylabel("xi_y/W"); axs2[0].legend(fontsize=6, ncol=2)
axs2[1].set_xlabel("W"); axs2[1].set_ylabel("xi_y (steps of 2 rows)"); axs2[1].set_yscale("log")
fig2.tight_layout(); fig2.savefig(os.path.join(HERE, "results", "xi_y.png"), dpi=130)
out["xi_y"] = {str(u): {str(W): v for W, v in d.items()} for u, d in tab.items()}


# ------------------------------------------------------------ crossings of x2(W,u) curves of different W (phase-interface / KT)
print("\n== crossings u_x(W1,W2) of x2(W1,u) and x2(W2,u) ==")
pairs = list(zip(Ws[:-1], Ws[1:]))
ux = {"stag": {}, "F": {}}
for model in ("stag", "F"):
    ugrid_m = sorted({k[2] for k in R if k[0] == model and k[3] == 2 and k[1] == Ws[-2]})
    for (W1, W2) in pairs:
        us = [u for u in ugrid_m if x2(model, W1, u) is not None and x2(model, W2, u) is not None]
        if len(us) < 4: continue
        d = CubicSpline(us, [x2(model, W2, u) - x2(model, W1, u) for u in us])
        uu = np.linspace(us[0], us[-1], 4000); dv = d(uu)
        sg = np.nonzero(dv[:-1] * dv[1:] < 0)[0]
        # the relevant sign change: x2(W2) < x2(W1) (Coulomb) -> x2(W2) > x2(W1) (beyond u_c); take the last - to + change
        up = [i for i in sg if dv[i] < 0 < dv[i + 1]]
        if up:
            ux[model][(W1, W2)] = float(brentq(lambda t: d(t), uu[up[-1]], uu[up[-1] + 1]))
    print(model, {f"{a}-{b}": round(v, 4) for (a, b), v in ux[model].items()})
out["u_cross_pairs"] = {m: {f"{a}-{b}": v for (a, b), v in d.items()} for m, d in ux.items()}

# ------------------------------------------------------------ large u: sliding liquid checks
print("\n== large u ==")
big = {}
for u in (3.0, 4.0, 5.0, 6.0):
    rows = []
    for W in Ws:
        l0, l1 = lam("stag", W, u, 0, 0), lam("stag", W, u, 0, 1)
        xv = x2("stag", W, u)
        if l0 is None or l1 is None: continue
        r = abs(l1 / l0)
        rows.append(dict(W=W, s_row=float(np.log(abs(l0)) - u * W / 2), ratio=float(r), r_minus_half_over_Wexp=float((r - 0.5) / (W * np.exp(-u))),
                         x2_over_W=None if xv is None else float(xv / W)))
    big[u] = rows
    print(f"u={u}:  " + "  ".join(f"W={r['W']}: s={r['s_row']:.4f}, (r-1/2)/(W e^-u)={r['r_minus_half_over_Wexp']:.3f}, x2/W={r['x2_over_W']:.4f}" for r in rows))
print("(ln 2 = %.4f)" % np.log(2))
out["large_u"] = {str(u): v for u, v in big.items()}
json.dump(out, open(os.path.join(HERE, "results", "tm_summary.json"), "w"), indent=1)
