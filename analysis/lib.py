"""Shared helpers: loading PT runs, jackknife statistics, Binder cumulants."""
import json, glob, os
import numpy as np

FIELDS = ["e", "e2", "mN", "mN2", "mN4", "mNabs", "ms2", "ms4", "msabs", "rdiag", "rnbr", "rviol", "rBal"]
FI = {f: i for i, f in enumerate(FIELDS)}


def load_run(prefix):
    meta = json.load(open(prefix + ".meta.json"))
    raw = np.genfromtxt(prefix + ".blocks.csv", delimiter=",", skip_header=1)
    nb, nrep = meta["nblocks"], meta["nrep"]
    data = raw[:, 3:].reshape(nb, nrep, len(FIELDS))
    return dict(meta=meta, T=np.array(meta["T"]), B=data, prefix=prefix)


def jackknife(blocks, func):
    """blocks: (nb, ...) array of block means; func maps mean-array -> value(s). Returns (value, error)."""
    nb = blocks.shape[0]
    full = func(blocks.mean(axis=0))
    tot = blocks.sum(axis=0)
    jk = np.array([func((tot - blocks[i]) / (nb - 1)) for i in range(nb)])
    err = np.sqrt((nb - 1) / nb * ((jk - jk.mean(axis=0)) ** 2).sum(axis=0))
    return full, err


def binder_fn(kind, ncomp=None):
    """Binder cumulant from the mean vector (nrep, nfields): kind 'N' (Neel) or 'S' (stripe, 2 components)."""
    if kind == "N":
        return lambda m: 1 - m[:, FI["mN4"]] / (3 * m[:, FI["mN2"]] ** 2)
    return lambda m: 1 - m[:, FI["ms4"]] / (2 * m[:, FI["ms2"]] ** 2)


def chi_fn(kind, NB, T):
    if kind == "N":
        return lambda m: NB * (m[:, FI["mN2"]] - m[:, FI["mNabs"]] ** 2) / T
    return lambda m: NB * (m[:, FI["ms2"]] - m[:, FI["msabs"]] ** 2) / T


def cv_fn(N, T):
    return lambda m: N * (m[:, FI["e2"]] - m[:, FI["e"]] ** 2) / T ** 2


def combine_seeds(runs):
    """Concatenate blocks of independent seeds (same T grid) -> single run dict."""
    B = np.concatenate([r["B"] for r in runs], axis=0)
    return dict(meta=runs[0]["meta"], T=runs[0]["T"], B=B, prefix=runs[0]["prefix"], nseeds=len(runs))


def find_runs(datadir, tag):
    """tag e.g. 'S1_K4_h0_eps0.012_Lb32'; seeds are suffixed _s1,_s2,..."""
    files = sorted(glob.glob(os.path.join(datadir, tag + "_s*.meta.json")))
    return [load_run(f[: -len(".meta.json")]) for f in files]
