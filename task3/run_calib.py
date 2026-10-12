#!/usr/bin/env python3
"""Calibration runs: Lieb ice (u=0, staggered model) and the uniform F model (A and B both c-type e^u) at u = 0.3, ln2, 0.9, 1.3."""
import math, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
jobs = []
for L in (16, 32, 64):
    jobs.append(dict(model="stag", u=0.0, L=L, neq=20000 if L < 64 else 5000, nmeas=100000 if L < 64 else 40000))
    for u in (0.3, math.log(2), 0.9, 1.3):
        jobs.append(dict(model="F", u=u, L=L, neq=40000 if L < 64 else 20000, nmeas=100000 if L < 64 else 40000))
jobs.sort(key=lambda j: -j["L"])
env = dict(os.environ, OMP_NUM_THREADS="1")
def work(j):
    out = os.path.join(ROOT, "task3", "data", "calib", f"{j['model']}_u{j['u']:.4f}_L{j['L']}")
    if os.path.exists(out + ".meta.json"): return out
    subprocess.run([os.path.join(ROOT, "build", "ice_pt"), "--L", str(j["L"]), "--ulist", repr(j["u"]), "--model", j["model"], "--neq", str(j["neq"]), "--nmeas", str(j["nmeas"]),
                    "--nblocks", "20", "--series_every", "100", "--seed", "11", "--init", "rowrand", "--out", out], env=env, check=True, capture_output=True)
    return out
with ThreadPoolExecutor(4) as ex:
    for r in ex.map(work, jobs): print("finished", os.path.basename(r), flush=True)
