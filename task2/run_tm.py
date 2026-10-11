#!/usr/bin/env python3
"""Run the task-2 transfer-matrix scans; results appended to data/results.jsonl (resumable).
   python run_tm.py [W ...]      default W = 8 12 16 20 24
"""
import json, os, sys, time
import numpy as np
import tm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data", "results.jsonl")
tm.set_threads(int(os.environ.get("TM_THREADS", "4")))

U_C1 = [round(0.8 + 0.1 * i, 1) for i in range(17)]            # 0.8 ... 2.4
U_C2 = [2.6, 2.8, 3.0]                                         # (1.6-3.0 step 0.2 is contained in the union grid)
U_C3 = [4.0, 5.0, 6.0]
U_STAG_M0 = U_C1 + U_C2 + U_C3                                 # top-4 eigenvalues of the M=0 sector
U_STAG_M2 = U_C1 + [3.0, 4.0, 5.0, 6.0]                        # top eigenvalue of M=2 (x_2)
U_F = [0.0, 0.5, 0.6, 0.65, float(np.log(2)), 0.75, 0.8, 0.85, 0.9]   # calibration C4 (uniform F model; u=0 = Lieb ice; KT at ln2)


def jobs(W):
    J = []
    if W >= 28:   # W=28 (40 million states in the M=0 sector): targeted subset only
        us = [1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.0]
        for u in [float(np.log(2)), 0.8]:
            J += [("F", W, u, 0, 1), ("F", W, u, 2, 1)]
        for u in us + [1.5, 1.7]:
            J.append(("stag", W, u, 2, 1))
        for u in us:
            J.append(("stag", W, u, 0, 3))
        return J
    for u in U_F:
        J += [("F", W, u, 0, 1), ("F", W, u, 2, 1)]
    for u in U_STAG_M2:
        J.append(("stag", W, u, 2, 1))
    for u in U_STAG_M0:
        J.append(("stag", W, u, 0, 4))
    return J


def key(j):
    return (j[0], j[1], round(j[2], 6), j[3])


def main():
    Ws = [int(a) for a in sys.argv[1:]] or [8, 12, 16, 20, 24]
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            r = json.loads(line)
            if r["k"] >= 1:
                done.add((r["model"], r["W"], round(r["u"], 6), r["M"], r["k"]))
    for W in Ws:
        for (model, W_, u, M, k) in jobs(W):
            if (model, W_, round(u, 6), M, k) in done:
                continue
            t0 = time.time()
            if W_ >= 16:   # sector-compressed solver (identical eigenvalues to the dense one, much faster and lighter)
                vals, leak = tm.top_eigs_sector(W_, u, M, k=k, model=model, tol=1e-10)
            else:
                vals, leak = tm.top_eigs(W_, u, M, k=k, model=model, tol=1e-10)
            rec = dict(model=model, W=W_, u=u, M=M, k=k, eig_re=[float(v.real) for v in vals],
                       eig_im=[float(v.imag) for v in vals], leak=leak, seconds=time.time() - t0)
            with open(OUT, "a") as f:
                f.write(json.dumps(rec) + "\n")
            print(f"{model} W={W_} u={u:.3f} M={M} k={k}: {rec['eig_re'][0]:.6g} ({rec['seconds']:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
