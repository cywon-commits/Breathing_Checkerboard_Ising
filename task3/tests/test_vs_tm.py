"""Independent check of the directed-loop MC against the task-2 transfer matrix: on the 2L x 2L torus
   Z_M(u) = Tr_{sector M} T^L  (T = one step of two rows, width W = 2L), M = 2 W_x.
Compares <nc>/L^2, <W_x^2> and P(W_x = 0) of the MC (ice_pt, a single u) with the exact values."""
import json, os, subprocess, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "task2"))
import tm

def exact(L, u, model="stag"):
    W = 2 * L
    cA, cB = tm.weights(u, model)
    pc = tm.popcount(W)
    Z = {}
    for M in range(-W, W + 1, 2):
        n = (W + M) // 2
        T = tm.sector_dense(W, cA, cB, n)
        Z[M] = np.trace(np.linalg.matrix_power(T, L))
    return Z

def moments(L, u, model="stag"):
    Z = exact(L, u, model); Zt = sum(Z.values())
    wx2 = sum((M / 2) ** 2 * z for M, z in Z.items()) / Zt
    p0 = Z[0] / Zt
    h = 1e-4
    lz = lambda uu: np.log(sum(exact(L, uu, model).values()))
    nc = (lz(u + h) - lz(u - h)) / (2 * h) / L ** 2          # d ln Z / du = <n_c(A)> (+ n_c(B) for the F model)
    return nc, wx2, p0

if __name__ == "__main__":
    tm.set_threads(2)
    ok = True
    for model in ("stag", "F"):
        for L in (4, 6):
            for u in (0.7, 1.6, 2.4):
                nc, wx2, p0 = moments(L, u, model)
                out = f"/tmp/ice_test_{model}_{L}_{u}"
                subprocess.run([os.path.join(ROOT, "build", "ice_pt"), "--L", str(L), "--ulist", str(u), "--model", model, "--neq", "5000", "--nmeas", "400000",
                                "--nblocks", "40", "--series_every", "100000", "--seed", "3", "--init", "rowrand", "--out", out], check=True, capture_output=True)
                B = np.genfromtxt(out + ".blocks.csv", delimiter=",", skip_header=1)
                hdr = open(out + ".blocks.csv").readline().strip().split(",")
                col = {h: i for i, h in enumerate(hdr)}
                def st(name, scale=1.0):
                    v = B[:, col[name]] * scale; return v.mean(), v.std(ddof=1) / np.sqrt(len(v))
                ncm = st("nc") if model == "stag" else (B[:, col["nc"]] + B[:, col["ncB"]]).mean(), None
                if model == "F":
                    v = B[:, col["nc"]] + B[:, col["ncB"]]; ncm = (v.mean(), v.std(ddof=1) / np.sqrt(len(v)))
                else:
                    ncm = st("nc")
                w2 = st("Wx2")
                good = abs(ncm[0] - nc) < 4.5 * ncm[1] + 1e-5 and abs(w2[0] - wx2) < 4.5 * w2[1] + 1e-5
                flag = "" if good else "  <-- MISMATCH"
                ok &= good
                print(f"{model} L={L} u={u}: <nc>/L^2 exact {nc:.5f} MC {ncm[0]:.5f}±{ncm[1]:.5f} | <Wx^2> exact {wx2:.5f} MC {w2[0]:.5f}±{w2[1]:.5f}{flag}")
    print("ALL OK" if ok else "FAILED")
    sys.exit(0 if ok else 1)
