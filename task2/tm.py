"""Sector transfer-matrix solver for the staggered six-vertex (ice) model (task 2).  See src/tm.cpp."""
import ctypes, os
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigs

_lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "libtm.so"))
_dp = ctypes.POINTER(ctypes.c_double)
_lib.apply_T.argtypes = [_dp, _dp, ctypes.c_int, ctypes.c_double, ctypes.c_double]
_lib.apply_A.argtypes = [_dp, ctypes.c_int, ctypes.c_double]
_lib.apply_B.argtypes = [_dp, ctypes.c_int, ctypes.c_double]
_lib.set_threads.argtypes = [ctypes.c_int]


def set_threads(n):
    _lib.set_threads(n)


def apply_T(v, W, cA, cB, out=None):
    v = np.ascontiguousarray(v, dtype=np.float64)
    if out is None:
        out = np.empty_like(v)
    _lib.apply_T(v.ctypes.data_as(_dp), out.ctypes.data_as(_dp), W, cA, cB)
    return out


def popcount(W):
    """popcount of all 2^W integers (uint8 array)"""
    pc = np.zeros(1, dtype=np.uint8)
    for _ in range(W):
        pc = np.concatenate([pc, pc + 1])
    return pc


def weights(u, model):
    """model 'stag': A vertices c = e^u, B vertices c = 1;  'F': both e^u (uniform F model)"""
    cA = float(np.exp(u))
    cB = 1.0 if model == "stag" else cA
    return cA, cB


def sector_dense(W, cA, cB, n):
    """dense matrix of T restricted to the sector with n up spins (feasible for W <= 14)"""
    pc = popcount(W)
    idx = np.nonzero(pc == n)[0]
    N = 1 << W
    T = np.empty((len(idx), len(idx)))
    out = np.empty(N)
    for c, b in enumerate(idx):
        e = np.zeros(N); e[b] = 1.0
        T[:, c] = apply_T(e, W, cA, cB, out)[idx]
    return T


def top_eigs(W, u, M, k=1, model="stag", tol=1e-10, seed=1, ncv=None, return_vec=False, dense_max=12):
    """largest-modulus eigenvalues (sorted by modulus) of T in the sector M = sum sigma = 2n - W."""
    n = (W + M) // 2
    cA, cB = weights(u, model)
    if W <= dense_max and not return_vec:
        ev = np.linalg.eigvals(sector_dense(W, cA, cB, n))
        ev = ev[np.argsort(-np.abs(ev))][:k]
        return ev, 0.0
    N = 1 << W
    pc = popcount(W)
    mask = pc == n
    dim = int(mask.sum())
    rng = np.random.default_rng(seed)
    v0 = np.zeros(N)
    v0[mask] = rng.random(dim) + 0.5
    v0 /= np.linalg.norm(v0)
    out = np.empty(N)
    maskf = mask.astype(np.float64)

    def mv(x):   # P T P : projector onto the sector guards against ARPACK's random restart vectors
        y = apply_T(x * maskf, W, cA, cB, out)
        return y * maskf
    op = LinearOperator((N, N), matvec=mv, dtype=np.float64)
    ncv = ncv or min(N, max(2 * k + 6, 30))
    vals, vecs = eigs(op, k=k, which="LM", v0=v0, tol=tol, ncv=ncv, maxiter=5000)
    order = np.argsort(-np.abs(vals))
    vals = vals[order]
    # leakage check: weight outside the sector
    leak = float(np.abs(vecs[~mask, :]).max()) if not mask.all() else 0.0
    if return_vec:
        return vals, leak, vecs[:, order]
    return vals, leak


# ---------------------------------------------------------------- sector-compressed solver (large W)
_slib = None


def _sector_lib():
    global _slib
    if _slib is None:
        _slib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), "libtm_sector.so"))
        _slib.sector_init.argtypes = [ctypes.c_int, ctypes.c_int]
        _slib.sector_init.restype = ctypes.c_int64
        _slib.sector_apply_T.argtypes = [_dp, _dp, ctypes.c_double, ctypes.c_double]
    return _slib


def top_eigs_sector(W, u, M, k=1, model="stag", tol=1e-10, seed=1, ncv=None):
    """Same as top_eigs but on the sector-compressed vector space (dimension C(W, (W+M)/2))."""
    lib = _sector_lib()
    n = (W + M) // 2
    cA, cB = weights(u, model)
    dim = lib.sector_init(W, n)
    rng = np.random.default_rng(seed)
    v0 = rng.random(dim) + 0.5
    v0 /= np.linalg.norm(v0)
    out = np.empty(dim)

    def mv(x):
        x = np.ascontiguousarray(x, dtype=np.float64)
        lib.sector_apply_T(x.ctypes.data_as(_dp), out.ctypes.data_as(_dp), cA, cB)
        return out.copy()
    op = LinearOperator((dim, dim), matvec=mv, dtype=np.float64)
    ncv = ncv or min(dim, max(2 * k + 6, 16))
    vals, _ = eigs(op, k=k, which="LM", v0=v0, tol=tol, ncv=ncv, maxiter=5000)
    order = np.argsort(-np.abs(vals))
    return vals[order], 0.0   # leakage impossible by construction (compressed space)
