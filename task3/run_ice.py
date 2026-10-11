#!/usr/bin/env python3
"""Production runs of the directed-loop MC: one process per (L, seed) covering the whole u grid (independent runs at each u,
no replica exchange: its swap acceptance collapses with L, see REPORT_task3.md).  Seeds 1,3 start from row-family ground states,
seeds 2,4 from column-family ground states.
   python run_ice.py [--workers 4] [--dry]
"""
import argparse, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
U = [0.8, 1.0, 1.2] + [round(1.25 + 0.05 * i, 2) for i in range(10)] + [round(1.8 + 0.1 * i, 1) for i in range(17)] + [4.0, 5.0]
LS = {16: 4, 24: 4, 32: 4, 48: 4, 64: 4, 96: 2}
NMEAS = {96: 60000}
INIT = {1: "rowrand", 2: "colrand", 3: "rowrand", 4: "colrand"}


def commit():
    try: return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT).decode().strip()
    except Exception: return "unknown"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=4); ap.add_argument("--dry", action="store_true")
    ap.add_argument("--outdir", default="data/prod"); a = ap.parse_args()
    out = os.path.join(ROOT, "task3", a.outdir); os.makedirs(out, exist_ok=True)
    jobs = []
    for L, ns in LS.items():
        for s in range(1, ns + 1):
            jobs.append(dict(L=L, seed=s, cost=L * L * len(U) * (NMEAS.get(L, 100000) + 20000)))
    jobs.sort(key=lambda j: -j["cost"])
    print(len(jobs), "jobs, est. CPU-h %.1f" % (sum(j["cost"] for j in jobs) * 40e-9 / 3600 * 1.0 / 1.0 * 4), flush=True)
    if a.dry: return
    ch = commit(); env = dict(os.environ, OMP_NUM_THREADS="1")

    def work(j):
        pre = os.path.join(out, f"L{j['L']}_s{j['seed']}")
        if os.path.exists(pre + ".meta.json"): return pre
        cmd = [os.path.join(ROOT, "build", "ice_pt"), "--L", j["L"], "--ulist", ",".join(map(str, U)), "--noswap", 1, "--neq", 20000, "--nmeas", NMEAS.get(j["L"], 100000),
               "--nblocks", 20, "--series_every", 20, "--seed", j["seed"], "--init", INIT[j["seed"]], "--out", pre, "--commit", ch]
        subprocess.run([str(c) for c in cmd], env=env, check=True); return pre
    with ThreadPoolExecutor(a.workers) as ex:
        for r in ex.map(work, jobs): print("finished", os.path.basename(r), time.strftime("%H:%M:%S"), flush=True)

if __name__ == "__main__": main()
