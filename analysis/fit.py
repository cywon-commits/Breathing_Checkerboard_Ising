#!/usr/bin/env python3
"""Second analysis stage: T_c table, slope fits (P3), FSS exponents (P4), violation-density scaling,
figures.  Needs results/analysis.json from analyze.py (run that first).

  python analysis/fit.py     -> results/fit.json, results/figs/*.png, prints tables
"""
import json, os, sys, glob
import numpy as np
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from lib import *  # noqa
from analyze import collect, tag, Tc_improved  # noqa

FIG = os.path.join(ROOT, "results", "figs")
os.makedirs(FIG, exist_ok=True)


def wls(x, y, sy):
    """weighted linear fit y = a + b x; returns a, b, errors (scaled by sqrt(chi2/dof) if > 1)."""
    w = 1 / sy ** 2
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    p = cov @ (A.T @ (w * y))
    res = y - A @ p
    dof = max(len(x) - 2, 1)
    chi2 = float((w * res ** 2).sum())
    scale = max(1.0, chi2 / dof) if len(x) > 2 else 1.0
    return p, np.sqrt(np.diag(cov) * scale), chi2, dof


def tc_table(ana):
    """T_c per (set, eps): crossing of the largest pair; systematic = |last - previous crossing|."""
    rows = []
    for sname, S in ana.items():
        Ek = S["cfg"]["Ek"]
        for eps, r in S["eps"].items():
            cr = [c for c in r["crossings"] if np.isfinite(c["Tx"])]
            if not cr:
                rows.append(dict(set=sname, eps=float(eps), Tc=np.nan)); continue
            last = cr[-1]
            drift = abs(cr[-1]["Tx"] - cr[-2]["Tx"]) if len(cr) > 1 else np.nan
            Lmax = r["Lb"][-1]
            rows.append(dict(set=sname, eps=float(eps), Tc=last["Tx"], stat=last["err"], drift=drift,
                             pair=(last["L1"], last["L2"]), Lmax=Lmax, spacing=r["spacing"],
                             small=bool(Lmax < 4 * r["spacing"]), pred=r["Tc_pred"], impr=r["Tc_improved"],
                             Ek=Ek, K=S["cfg"]["K"], h=S["cfg"]["h"], minswap=last["minswap"],
                             all_crossings=[(c["L1"], c["L2"], c["Tx"], c["err"]) for c in cr]))
    return rows


def slope_fits(rows):
    out = {}
    for sname in sorted(set(r["set"] for r in rows)):
        R = [r for r in rows if r["set"] == sname and np.isfinite(r["Tc"]) and r["eps"] > 0 and sname in ("S1", "S2", "S3")]
        if len(R) < 2:
            continue
        R.sort(key=lambda r: r["eps"])
        x = np.array([np.log(1 / r["eps"]) for r in R])
        y = np.array([1 / r["Tc"] for r in R])
        sTc = np.array([np.hypot(r["stat"], 0 if np.isnan(r["drift"]) else r["drift"]) for r in R])
        sy = np.maximum(sTc / np.array([r["Tc"] for r in R]) ** 2, 1e-4)
        Ek = R[0]["Ek"]
        d = dict(Ek=Ek, pred=1 / Ek, eps=[r["eps"] for r in R], local=[])
        p, e, chi2, dof = wls(x, y, sy)
        d["all"] = dict(slope=p[1], err=e[1], intercept=p[0], chi2=chi2, dof=dof, n=len(R))
        good = [i for i, r in enumerate(R) if not r["small"]]
        if len(good) >= 2:
            p, e, chi2, dof = wls(x[good], y[good], sy[good])
            d["valid"] = dict(slope=p[1], err=e[1], intercept=p[0], chi2=chi2, dof=dof, n=len(good))
        for i in range(len(R) - 1):  # local slope between adjacent eps
            s = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
            es = np.hypot(sy[i], sy[i + 1]) / abs(x[i + 1] - x[i])
            d["local"].append(dict(eps_pair=(R[i]["eps"], R[i + 1]["eps"]), slope=s, err=es,
                                   flagged=R[i]["small"] or R[i + 1]["small"]))
        out[sname] = d
    return out


def fss(ana):
    """nu from d U/dT at T_c and gamma/nu from chi(T_c) vs L for each (set, eps>0)."""
    out = {}
    for sname, S in ana.items():
        cfg = yaml.safe_load(open(os.path.join(ROOT, "configs", f"{sname}.yaml")))
        data = collect(sname, cfg)
        for eps, runs in data.items():
            if eps <= 0 or len(runs) < 3:
                continue
            r = S["eps"].get(str(eps))
            if not r:
                continue
            cr = [c for c in r["crossings"] if np.isfinite(c["Tx"])]
            if not cr:
                continue
            Tc = cr[-1]["Tx"]
            Lbs = sorted(runs)
            nb = min(runs[L]["B"].shape[0] for L in Lbs)
            blocks = {}
            for L in Lbs:
                B = runs[L]["B"]; k = B.shape[0] // nb
                blocks[L] = B[: k * nb].reshape(nb, k, *B.shape[1:]).mean(axis=1)

            def obs(L, m):
                T = runs[L]["T"]
                NB = L * L
                U = binder_fn("N")(m)
                chi = chi_fn("N", NB, T)(m)
                cu = CubicSpline(T, U)
                cc = CubicSpline(T, chi)
                return float(cu(Tc, 1)), float(cc(Tc))   # dU/dT and chi at Tc

            vals = {L: [] for L in Lbs}
            for k in range(nb):
                for L in Lbs:
                    m = (blocks[L].sum(0) - blocks[L][k]) / (nb - 1)
                    vals[L].append(obs(L, m))
            res = {}
            full = {L: obs(L, blocks[L].mean(0)) for L in Lbs}
            x = np.log(np.array(Lbs, float))
            def fit(idx):
                ys = np.array([[np.log(abs(vals[L][k][idx])) for L in Lbs] for k in range(nb)])
                y0 = np.array([np.log(abs(full[L][idx])) for L in Lbs])
                sl0 = np.polyfit(x, y0, 1)[0]
                jk = np.array([np.polyfit(x, ys[k], 1)[0] for k in range(nb)])
                err = np.sqrt((nb - 1) / nb * ((jk - jk.mean()) ** 2).sum())
                return sl0, err
            inv_nu, e1 = fit(0)
            gnu, e2 = fit(1)
            out[f"{sname}_eps{eps:g}"] = dict(Tc=Tc, Lb=Lbs, one_over_nu=inv_nu, one_over_nu_err=e1,
                                              nu=1 / inv_nu, gamma_over_nu=gnu, gamma_over_nu_err=e2)
    return out


def violation_scaling(ana):
    """ln(rho_viol) vs 1/T in a +-8% window around T_c (largest size): slope should be -E_k; ratio at T_c."""
    out = {}
    for sname, S in ana.items():
        cfg = yaml.safe_load(open(os.path.join(ROOT, "configs", f"{sname}.yaml")))
        data = collect(sname, cfg)
        for eps, runs in data.items():
            r = S["eps"].get(str(eps))
            if not r or eps <= 0:
                continue
            cr = [c for c in r["crossings"] if np.isfinite(c["Tx"])]
            if not cr:
                continue
            Tc = cr[-1]["Tx"]
            L = max(runs)
            B = runs[L]["B"]
            T = runs[L]["T"]
            rv = B[:, :, FI["rviol"]].mean(0)
            sel = (T > 0.92 * Tc) & (T < 1.08 * Tc)
            if sel.sum() < 3:
                continue
            p = np.polyfit(1 / T[sel], np.log(rv[sel]), 1)
            rv_tc = float(np.exp(np.polyval(p, 1 / Tc)))
            Ek = S["cfg"]["Ek"]
            out[f"{sname}_eps{eps:g}"] = dict(Tc=Tc, Lb=L, slope=float(p[0]), Ek=Ek, rho_viol_Tc=rv_tc,
                                              boltz=float(np.exp(-Ek / Tc)), ratio=rv_tc / float(np.exp(-Ek / Tc)))
    return out


# ------------------------------------------------------------------ figures
def fig_binder(ana):
    for sname, S in ana.items():
        cfg = yaml.safe_load(open(os.path.join(ROOT, "configs", f"{sname}.yaml")))
        data = collect(sname, cfg)
        eps_list = sorted(data, key=lambda e: -abs(e))
        n = len(eps_list)
        fig, axs = plt.subplots(1, n, figsize=(3.6 * n, 3.4), squeeze=False)
        for ax, eps in zip(axs[0], eps_list):
            runs = data[eps]
            kind = "N" if eps > 0 else "S"
            for L in sorted(runs):
                U, dU = jackknife(runs[L]["B"], binder_fn(kind))
                ax.errorbar(runs[L]["T"], U, dU, fmt="o-", ms=3, lw=0.8, label=f"L={2*L}")
            rr = S["eps"].get(str(eps))
            if rr:
                for c in rr["crossings"]:
                    if np.isfinite(c["Tx"]):
                        ax.axvline(c["Tx"], color="gray", lw=0.6, ls=":")
                ax.axvline(rr["Tc_pred"], color="k", lw=0.8, ls="--")
            ax.set_title(f"{sname} ε={eps:g}", fontsize=9)
            ax.set_xlabel("T")
            ax.set_ylabel("Binder U")
            ax.legend(fontsize=6)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, f"binder_{sname}.png"), dpi=130)
        plt.close(fig)


def fig_slopes(rows, fits):
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    colors = {"S1": "C0", "S2": "C1", "S3": "C2"}
    for sname, d in fits.items():
        R = sorted([r for r in rows if r["set"] == sname and np.isfinite(r["Tc"]) and r["eps"] > 0], key=lambda r: r["eps"])
        x = np.array([np.log(1 / r["eps"]) for r in R])
        y = np.array([1 / r["Tc"] for r in R])
        ey = np.array([np.hypot(r["stat"], 0 if np.isnan(r["drift"]) else r["drift"]) / r["Tc"] ** 2 for r in R])
        small = np.array([r["small"] for r in R])
        c = colors[sname]
        ax.errorbar(x[~small], y[~small], ey[~small], fmt="o", color=c, label=f"{sname} (K={R[0]['K']:g}, h={R[0]['h']:g})")
        ax.errorbar(x[small], y[small], ey[small], fmt="o", mfc="none", color=c)
        # predicted slope 1/E_k through the leading-log formula:  1/Tc = (1/Ek) ln(Ek/(6 eps)) = (1/Ek)(ln(1/eps) + ln(Ek/6))
        xx = np.linspace(x.min() - 0.2, x.max() + 0.2, 10)
        ax.plot(xx, (xx + np.log(d["Ek"] / 6)) / d["Ek"], "--", color=c, lw=1)
    ax.set_xlabel("ln(1/ε)")
    ax.set_ylabel("1/T_c")
    ax.set_title("1/T_c vs ln(1/ε); dashed: prediction (slope 1/E_k); open: below size bound", fontsize=8)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "slopes.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    ana = json.load(open(os.path.join(ROOT, "results", "analysis.json")))
    rows = tc_table(ana)
    fits = slope_fits(rows)
    fs = fss(ana)
    viol = violation_scaling(ana)
    for r in rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k != "all_crossings"})
    print(json.dumps(fits, indent=1, default=float))
    print(json.dumps(fs, indent=1, default=float))
    print(json.dumps(viol, indent=1, default=float))
    json.dump(dict(rows=rows, fits=fits, fss=fs, viol=viol), open(os.path.join(ROOT, "results", "fit.json"), "w"), indent=1, default=float)
    fig_binder(ana)
    fig_slopes(rows, fits)
