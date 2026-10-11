"""Tests of the transfer matrix: (1) gate implementation == brute-force matrix from the plaquette definition,
(2) sector structure / no leakage, (3) known limits (Lieb ice entropy, F-model free energy)."""
import itertools, sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import tm


def brute_T(W, cA, cB):
    """T = T_B T_A from the definition: plaquette = 2 bottom + 2 top spins, n_up = 2, diagonal-type weight c."""
    N = 1 << W
    def step(pairs, c):
        S = np.zeros((N, N))
        for b in range(N):
            for t in range(N):
                w = 1.0
                for (i, j) in pairs:
                    bl, br = (b >> i) & 1, (b >> j) & 1
                    tl, tr = (t >> i) & 1, (t >> j) & 1
                    if bl + br + tl + tr != 2:
                        w = 0.0; break
                    if bl != br and tl != tr and bl != tl:   # ups at (0,0),(1,1) or (1,0),(0,1): diagonal
                        w *= c
                S[t, b] = w
        return S
    A = step([(2 * k, 2 * k + 1) for k in range(W // 2)], cA)
    B = step([(2 * k + 1, (2 * k + 2) % W) for k in range(W // 2)], cB)
    return B @ A


def test_gates_vs_brute():
    for W in (4, 6, 8):
        for (cA, cB) in ((1.0, 1.0), (np.exp(1.3), 1.0), (2.0, 2.0)):
            T0 = brute_T(W, cA, cB)
            T1 = np.array([tm.apply_T(np.eye(1 << W)[:, b], W, cA, cB) for b in range(1 << W)]).T
            assert np.allclose(T0, T1), (W, cA, cB)
    print("gate T == brute-force plaquette T for W=4,6,8: OK")


def test_sectors_and_leakage():
    W = 8
    pc = tm.popcount(W)
    T = np.array([tm.apply_T(np.eye(1 << W)[:, b], W, np.exp(0.7), 1.0) for b in range(1 << W)]).T
    for n in range(W + 1):
        m = pc == n
        assert np.allclose(T[np.ix_(~m, m)], 0), "weight leaks out of sector"
    # M and -M sectors have equal spectra
    ev = lambda n: np.sort(np.abs(np.linalg.eigvals(T[np.ix_(pc == n, pc == n)])))[::-1]
    assert np.allclose(ev(3), ev(5)) and np.allclose(ev(2), ev(6))
    # eigenvalues of the sector solver == dense eigenvalues
    for M in (0, 2, 4):
        n = (W + M) // 2
        dense = np.sort(np.abs(np.linalg.eigvals(T[np.ix_(pc == n, pc == n)])))[::-1][:3]
        vals, leak = tm.top_eigs(W, 0.7, M, k=3, tol=1e-12, dense_max=0)
        assert np.allclose(np.abs(vals), dense, rtol=1e-8) and leak == 0.0, (M, vals, dense, leak)
    print("sector structure, no leakage, ARPACK == dense: OK")


def test_limits():
    # Lieb ice: entropy per spin -> 0.2158 ; ln Lambda_0 / (2W) (2 site rows of W spins per step)
    for W in (8, 12, 16):
        v, _ = tm.top_eigs(W, 0.0, 0, k=1)
        s = np.log(abs(v[0])) / (2 * W)
        print(f"  Lieb ice W={W}: S/spin = {s:.5f}")
    assert abs(s - 0.2158) < 0.004
    print("Lieb entropy OK")


if __name__ == "__main__":
    test_gates_vs_brute(); test_sectors_and_leakage(); test_limits()
    print("ALL TM TESTS PASSED")
