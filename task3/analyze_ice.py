#!/usr/bin/env python3
"""Analysis of the task-3 production data (data/prod): nematic order parameter Q, Binder crossings, winding fluctuations,
correlation lengths and the anisotropy exponent z.  Jackknife over blocks (seeds concatenated).
   python analyze_ice.py [datadir] -> results/ice_summary.json, results/figs/*.png
"""
import glob, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "data", "prod")
FIG = os.path.join(HERE, "results", "figs"); os.makedirs(FIG, exist_ok=True)
NAMES = ["nc", "ncB", "Q", "Q2", "Q4", "Qabs", "Wx2", "Wy2", "Wx4", "Wy4", "Wpar2", "Wper2", "Spar0", "Spar1", "Sper0", "Sper1", "Sx0", "Sx1", "Sy0", "Sy1", "Wsym2"]
FI = {n: i for i, n in enumerate(NAMES)}


def load(prefix):
    meta = json.load(open(prefix + ".meta.json"))
    raw = np.genfromtxt(prefix + ".blocks.csv", delimiter=",", skip_header=1)
    nb, nu = meta["nblocks"], len(meta["u"])
    B = raw[:, 3:].reshape(nb, nu, len(NAMES))
    # The production runs (data/prod) were made with a driver bug that started the measurement only after 2*neq sweeps, so
    # the last blocks were never filled (they hold garbage ~1e-309): drop them (the filled blocks are correct block averages).
    keep = np.abs(B[:, :, FI["nc"]]).sum(axis=1) > 1e-3   # unfilled blocks hold uninitialised memory (~1e-309), not exact zeros
    return meta, B[keep]


def collect(datadir, L):
    fs = sorted(glob.glob(os.path.join(datadir, f"L{L}_s*.meta.json")))
    if not fs: return None
    parts = [load(f[:-len(".meta.json")]) for f in fs]
    u = np.array(parts[0][0]["u"])
    B = np.concatenate([p[1] for p in parts], axis=0)
    return dict(L=L, u=u, B=B, metas=[p[0] for p in parts], seeds=len(fs))


def jk(B, fn):
    nb = B.shape[0]; tot = B.sum(0)
    full = fn(B.mean(0))
    js = np.array([fn((tot - B[i]) / (nb - 1)) for i in range(nb)])
    return full, np.sqrt((nb - 1) / nb * ((js - js.mean(0)) ** 2).sum(0))


def binderQ(m): return 1 - m[:, FI["Q4"]] / (3 * m[:, FI["Q2"]] ** 2)
def xi2(m, L, par):   # second-moment correlation length of the diagonal field along the parallel / perpendicular direction
    S0, S1 = (m[:, FI["Spar0"]], m[:, FI["Spar1"]]) if par else (m[:, FI["Sper0"]], m[:, FI["Sper1"]])
    return np.sqrt(np.maximum(S0 / S1 - 1, 0) / (4 * np.sin(np.pi / L) ** 2))


def g_from_w2(w2):
    """Coulomb-gas stiffness g with P(W) ~ exp(-pi g W^2) reproducing <W^2> = w2 (isotropic torus)"""
    w = np.arange(-8, 9)
    f = lambda g: (w ** 2 * np.exp(-np.pi * g * w ** 2)).sum() / np.exp(-np.pi * g * w ** 2).sum() - w2
    if w2 <= 1e-12: return np.nan
    try: return brentq(f, 1e-3, 30)
    except ValueError: return np.nan


def crossing(u, U1, U2, lo=-0.1, hi=0.66):
    c1, c2 = CubicSpline(u, U1), CubicSpline(u, U2)
    uu = np.linspace(u[0], u[-1], 3000); d = c1(uu) - c2(uu)
    r = []
    for i in np.nonzero(d[:-1] * d[1:] < 0)[0]:
        x = brentq(lambda t: c1(t) - c2(t), uu[i], uu[i + 1])
        if lo < 0.5 * (c1(x) + c2(x)) < hi: r.append(x)
    return r


def main():
    Ls = [L for L in (16, 24, 32, 48, 64, 96) if collect(DATA, L)]
    D = {L: collect(DATA, L) for L in Ls}
    u = D[Ls[0]]["u"]
    out = dict(Ls=Ls, u=u.tolist(), seeds={L: D[L]["seeds"] for L in Ls})
    # ------- per (L,u) quantities
    Q2 = {L: jk(D[L]["B"], lambda m: m[:, FI["Q2"]]) for L in Ls}
    UQ = {L: jk(D[L]["B"], binderQ) for L in Ls}
    W2 = {L: jk(D[L]["B"], lambda m: m[:, FI["Wsym2"]]) for L in Ls}
    Wpar = {L: jk(D[L]["B"], lambda m: m[:, FI["Wpar2"]]) for L in Ls}
    Wper = {L: jk(D[L]["B"], lambda m: m[:, FI["Wper2"]]) for L in Ls}
    NC = {L: jk(D[L]["B"], lambda m: m[:, FI["nc"]]) for L in Ls}
    Xpar = {L: jk(D[L]["B"], lambda m, L=L: xi2(m, L, True)) for L in Ls}
    Xper = {L: jk(D[L]["B"], lambda m, L=L: xi2(m, L, False)) for L in Ls}
    out["Q2"] = {L: Q2[L][0].tolist() for L in Ls}; out["UQ"] = {L: UQ[L][0].tolist() for L in Ls}
    out["UQ_err"] = {L: UQ[L][1].tolist() for L in Ls}
    out["W2sym"] = {L: W2[L][0].tolist() for L in Ls}
    # ------- Binder crossings
    cross = {}
    for a, b in zip(Ls[:-1], Ls[1:]):
        r = crossing(u, UQ[a][0], UQ[b][0])
        nb = min(D[a]["B"].shape[0], D[b]["B"].shape[0]); xs = []
        for k in range(nb):
            Ba, Bb = D[a]["B"][:nb], D[b]["B"][:nb]
            ua = binderQ((Ba.sum(0) - Ba[k]) / (nb - 1)); ub = binderQ((Bb.sum(0) - Bb[k]) / (nb - 1))
            rr = crossing(u, ua, ub); xs.append(rr[0] if rr else np.nan)
        cross[f"{a}-{b}"] = dict(roots=r, jackknife_err=float(np.sqrt((nb - 1) / nb * np.nansum((np.array(xs) - np.nanmean(xs)) ** 2))) if r else None)
    out["UQ_crossings"] = cross
    # ------- exponents: Q^2 ~ L^-eta, xi ~ L^a
    eta = {}
    for i, uu in enumerate(u):
        q = np.array([Q2[L][0][i] for L in Ls]); lL = np.log(Ls)
        eta[float(uu)] = float(-np.polyfit(lL[-3:], np.log(q[-3:]), 1)[0])
    out["Q2_exponent_last3"] = eta
    zt = {}
    for i, uu in enumerate(u):
        xp = np.array([Xpar[L][0][i] for L in Ls]); xq = np.array([Xper[L][0][i] for L in Ls])
        ok = np.isfinite(xp) & np.isfinite(xq) & (xp > 0) & (xq > 0)
        if ok.sum() >= 3:
            lxp, lxq, lL = np.log(xp[ok][-3:]), np.log(xq[ok][-3:]), np.log(np.array(Ls)[ok][-3:])
            zt[float(uu)] = dict(a_par=float(np.polyfit(lL, lxp, 1)[0]), a_per=float(np.polyfit(lL, lxq, 1)[0]),
                                 z_ratio=float(np.polyfit(lxq, lxp, 1)[0]))
    out["xi_exponents"] = zt
    out["geff"] = {L: [g_from_w2(w) for w in W2[L][0]] for L in Ls}
    out["xi_par"] = {L: Xpar[L][0].tolist() for L in Ls}; out["xi_per"] = {L: Xper[L][0].tolist() for L in Ls}
    out["Wpar2"] = {L: Wpar[L][0].tolist() for L in Ls}; out["Wper2"] = {L: Wper[L][0].tolist() for L in Ls}
    out["nc"] = {L: NC[L][0].tolist() for L in Ls}
    # ------- ergodicity / family comparison
    erg = {}
    for L in Ls:
        sw = np.array([m["Q_sign_switches_per_config"] for m in D[L]["metas"]])   # (seeds, nu) switches per configuration (= per u)
        # |Q|^2 averaged separately for row-start seeds (odd) and column-start seeds (even)
        nbp = D[L]["B"].shape[0] // D[L]["seeds"]
        rows = [D[L]["B"][(s * nbp):((s + 1) * nbp)] for s in range(D[L]["seeds"]) if s % 2 == 0]
        cols_ = [D[L]["B"][(s * nbp):((s + 1) * nbp)] for s in range(D[L]["seeds"]) if s % 2 == 1]
        if not rows or not cols_: continue
        rowB, colB = np.concatenate(rows), np.concatenate(cols_)
        qr = rowB[:, :, FI["Q2"]].mean(0); qc = colB[:, :, FI["Q2"]].mean(0)
        erg[L] = dict(switches=sw.tolist(), Q2_row_start=qr.tolist(), Q2_col_start=qc.tolist(), Q_row_start=rowB[:, :, FI["Q"]].mean(0).tolist(), Q_col_start=colB[:, :, FI["Q"]].mean(0).tolist())
    out["ergodicity"] = erg
    json.dump(out, open(os.path.join(HERE, "results", "ice_summary.json"), "w"), indent=1)
    # ------- figures
    cols = plt.cm.viridis(np.linspace(0, 0.9, len(Ls)))
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for L, c in zip(Ls, cols):
        ax[0, 0].errorbar(u, UQ[L][0], UQ[L][1], fmt="o-", ms=3, color=c, label=f"L={L}")
        ax[0, 1].plot(u, Q2[L][0], "o-", ms=3, color=c, label=f"L={L}")
        ax[0, 2].plot(u, W2[L][0], "o-", ms=3, color=c, label=f"L={L}")
        ax[1, 0].plot(u, Wpar[L][0], "o-", ms=3, color=c); ax[1, 0].plot(u, Wper[L][0], "s--", ms=3, color=c)
        ax[1, 1].plot(u, Xpar[L][0], "o-", ms=3, color=c); ax[1, 1].plot(u, Xper[L][0], "s--", ms=3, color=c)
        ax[1, 2].plot(u, np.array(out["geff"][L]), "o-", ms=3, color=c)
    ax[0, 0].set_title("Binder cumulant of Q"); ax[0, 0].axvline(1.38, color="r", ls=":"); ax[0, 0].set_ylabel("U_Q")
    ax[0, 1].set_title("<Q^2>"); ax[0, 1].set_yscale("log")
    ax[0, 2].set_title("<W^2> (x,y average)"); ax[0, 2].set_yscale("log")
    ax[1, 0].set_title("<W_par^2> (solid), <W_perp^2> (dashed)"); ax[1, 0].set_yscale("log")
    ax[1, 1].set_title("xi_par (solid), xi_perp (dashed)"); ax[1, 1].set_yscale("log")
    ax[1, 2].set_title("g_eff from <W^2>: P(W)~exp(-pi g W^2)"); ax[1, 2].set_ylim(0, 6)
    for a_ in ax.flat: a_.set_xlabel("u"); a_.axvline(1.38, color="r", ls=":", lw=0.7)
    ax[0, 0].legend(fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(FIG, "ice_overview.png"), dpi=120)
    print("crossings of U_Q:", {k: (np.round(v["roots"], 3).tolist(), v["jackknife_err"]) for k, v in cross.items()})


if __name__ == "__main__":
    main()
