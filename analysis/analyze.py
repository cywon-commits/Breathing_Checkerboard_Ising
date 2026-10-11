#!/usr/bin/env python3
"""Binder crossings, T_c, slope fits, FSS and violation-density checks for all sets.

  python analysis/analyze.py            -> results/analysis.json and tables printed to stdout
Jackknife over PT blocks (independent seeds concatenated) gives all statistical errors.
"""
import glob, json, os, re, sys
import numpy as np
import yaml
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from lib import *  # noqa

RES = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)


def tag(name, K, h, eps, Lb):
    return f"{name}_K{K:g}_h{h:g}_eps{eps:g}_Lb{Lb}"


def spline_eval(T, y, Tq):
    return CubicSpline(T, y)(Tq)


def crossing(T1, U1, T2, U2, guess=None, lo=0.05, hi=0.64):
    """Crossing of two spline-interpolated Binder curves inside the common T range, restricted to the
    transition region lo < U < hi.  Returns nan if none."""
    a, b = max(T1[0], T2[0]), min(T1[-1], T2[-1])
    Tq = np.linspace(a, b, 600)
    c1, c2 = CubicSpline(T1, U1), CubicSpline(T2, U2)
    d = c1(Tq) - c2(Tq)
    roots = []
    for i in range(len(Tq) - 1):
        if d[i] == 0 or d[i] * d[i + 1] < 0:
            try:
                r = brentq(lambda t: c1(t) - c2(t), Tq[i], Tq[i + 1])
            except ValueError:
                continue
            u = 0.5 * (c1(r) + c2(r))
            if lo < u < hi:
                roots.append(r)
    if not roots:
        return np.nan
    # several crossings can occur in the noisy disordered tail (U ~ 0): take the one on the ordered side
    # of the transition, i.e. the lowest-T root (guess only breaks ties when given)
    return min(roots)


def collect(set_name, cfg):
    """data[(eps)][Lb] = combined run"""
    out = {}
    for eps in cfg["eps"]:
        for Lb in cfg["Lb"]:
            runs = find_runs(os.path.join(ROOT, "data", set_name), tag(set_name, cfg["K"], cfg["h"], eps, Lb))
            if runs:
                out.setdefault(eps, {})[Lb] = combine_seeds(runs)
    return out


def analyse_eps(runs, eps, guess):
    kind = "N" if eps > 0 else "S"
    Lbs = sorted(runs)
    ufun = binder_fn(kind)
    res = dict(Lb=Lbs, crossings=[], kind=kind)
    # jackknife for the crossings (all sizes must have the same number of blocks)
    runs = {k: dict(v) for k, v in runs.items()}
    nb = min(runs[Lb]["B"].shape[0] for Lb in Lbs)
    for Lb in Lbs:  # merge consecutive blocks so that all sizes share nb jackknife blocks
        B = runs[Lb]["B"]
        k = B.shape[0] // nb
        runs[Lb]["B"] = B[: k * nb].reshape(nb, k, *B.shape[1:]).mean(axis=1)
    U = {Lb: jackknife(runs[Lb]["B"], ufun) for Lb in Lbs}
    for L1, L2 in zip(Lbs[:-1], Lbs[1:]):
        T1, T2 = runs[L1]["T"], runs[L2]["T"]
        x0 = crossing(T1, U[L1][0], T2, U[L2][0], guess)
        # jackknife
        xs = []
        B1, B2 = runs[L1]["B"][:nb], runs[L2]["B"][:nb]
        for k in range(nb):
            m1 = (B1.sum(0) - B1[k]) / (nb - 1)
            m2 = (B2.sum(0) - B2[k]) / (nb - 1)
            xs.append(crossing(T1, ufun(m1), T2, ufun(m2), x0 if np.isfinite(x0) else guess))
        xs = np.array(xs)
        err = np.sqrt((nb - 1) / nb * np.nansum((xs - np.nanmean(xs)) ** 2)) if np.isfinite(x0) else np.nan
        sw1, sw2 = min(runs[L1]["meta"]["swap_acc"]), min(runs[L2]["meta"]["swap_acc"])
        res["crossings"].append(dict(L1=L1, L2=L2, Tx=float(x0), err=float(err), minswap=float(min(sw1, sw2))))
    return res, U


def Tc_improved(eps, Ja, h):
    f = lambda T: T * np.log(1 + np.exp(-4 * (Ja - h) / T) + np.exp(-4 * (Ja + h) / T)) - 12 * abs(eps)
    return brentq(f, 1e-3, 50)


def main():
    allres = {}
    for cf in sorted(glob.glob(os.path.join(ROOT, "configs", "S*.yaml"))):
        cfg = yaml.safe_load(open(cf))
        n = cfg["name"]
        data = collect(n, cfg)
        if not data:
            continue
        allres[n] = dict(cfg=dict(K=cfg["K"], h=cfg["h"], Ek=cfg["Ek"], region=cfg["region"]), eps={})
        pil = json.load(open(os.path.join(ROOT, "data", "pilot", "tc_pilot.json"))) if os.path.exists(os.path.join(ROOT, "data", "pilot", "tc_pilot.json")) else {}
        for eps, runs in data.items():
            if len(runs) < 2:
                continue
            guess = pil.get(f"{n}_eps{eps:g}") or cfg["Tc_pred"][cfg["eps"].index(eps)]
            r, U = analyse_eps(runs, eps, guess)
            r["Tc_pred"] = cfg["Tc_pred"][cfg["eps"].index(eps)]
            r["Tc_improved"] = Tc_improved(eps, 1.0, cfg["h"])
            r["spacing"] = r["Tc_pred"] / (6 * abs(eps))  # kink spacing in B units
            allres[n]["eps"][str(eps)] = r
            print(n, f"eps={eps:+.4f}", "sizes", r["Lb"], "crossings",
                  [f"({c['L1']},{c['L2']}):{c['Tx']:.4f}±{c['err']:.4f}" for c in r["crossings"]],
                  f"pred={r['Tc_pred']:.3f} impr={r['Tc_improved']:.3f}")
    json.dump(allres, open(os.path.join(RES, "analysis.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
