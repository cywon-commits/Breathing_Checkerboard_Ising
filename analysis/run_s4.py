#!/usr/bin/env python3
"""Expanded S4 scan (K = 0.1..0.3, h = 0, 0.5; configs/task1_S4_expanded.yaml): pilot then production with bcising.
  python analysis/run_s4.py pilot | prod [--workers 1]
"""
import argparse, itertools, json, math, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
import yaml
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from lib import load_run, jackknife, binder_fn
BIN = os.path.join(ROOT, "build", "bcising")
CFG = yaml.safe_load(open(os.path.join(ROOT, "configs", "task1_S4_expanded.yaml")))
PIL = os.path.join(ROOT, "data", "S4x_pilot", "tc_pilot.json")


def tag(K, h, eps, Lb, seed): return f"S4x_K{K:g}_h{h:g}_eps{eps:g}_Lb{Lb}_s{seed}"


def level_cross(r, level=0.33):
    U, _ = jackknife(r["B"], binder_fn("N")); T = r["T"]
    for i in range(len(T) - 1):
        if U[i] > level >= U[i + 1]:
            return float(T[i] + (U[i] - level) / (U[i] - U[i + 1]) * (T[i + 1] - T[i]))
    return None


def commit():
    try: return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT).decode().strip()
    except Exception: return "unknown"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("stage"); ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args(); ch = commit(); jobs = []
    pil = json.load(open(PIL)) if os.path.exists(PIL) else {}
    for K, h, eps in itertools.product(CFG["K"], CFG["h"], CFG["eps"]):
        if a.stage == "pilot":
            p = CFG["pilot"]
            jobs.append(dict(out=os.path.join(ROOT, "data", "S4x_pilot", tag(K, h, eps, p["Lb"], 1)), K=K, h=h, eps=eps, Lb=p["Lb"], seed=1,
                             Tmin=p["Tmin"], Tmax=p["Tmax"], nrep=p["nrep"], neq=p["neq"], nmeas=p["nmeas"], nblocks=10, key=f"K{K:g}_h{h:g}_eps{eps:g}"))
        else:
            P = CFG["prod"]; key = f"K{K:g}_h{h:g}_eps{eps:g}"
            tc = pil.get(key)
            if not tc: print("no pilot T_c for", key, "-> skipped"); continue
            w = P["window"]
            Lbs = list(CFG["Lb"]) + ([64] if eps in CFG["Lb64_eps"] else [])
            for Lb in Lbs:
                for seed in range(1, (1 if Lb == 64 else P["seeds"]) + 1):
                    nrep = max(P["nrep_min"], int(math.ceil(math.log((1 + w) / (1 - w)) / (1.2 / Lb))) + 1)
                    jobs.append(dict(out=os.path.join(ROOT, "data", "S4x", tag(K, h, eps, Lb, seed)), K=K, h=h, eps=eps, Lb=Lb, seed=seed,
                                     Tmin=tc * (1 - w), Tmax=tc * (1 + w), nrep=nrep, neq=P["neq"], nmeas=P["nmeas"], nblocks=20, key=key))
    for j in jobs: j["cost"] = j["nrep"] * (j["neq"] + j["nmeas"]) * 4 * j["Lb"] ** 2 * 100e-9
    print(len(jobs), "jobs, ~%.2f CPU-h" % (sum(j["cost"] for j in jobs) / 3600), flush=True)
    jobs.sort(key=lambda j: -j["cost"])
    env = dict(os.environ, OMP_NUM_THREADS="1")

    def work(j):
        os.makedirs(os.path.dirname(j["out"]), exist_ok=True)
        if os.path.exists(j["out"] + ".meta.json"): return j["out"]
        cmd = [BIN, "--Lb", j["Lb"], "--K", j["K"], "--h", j["h"], "--eps", j["eps"], "--Tmin", j["Tmin"], "--Tmax", j["Tmax"], "--nrep", j["nrep"],
               "--neq", j["neq"], "--nmeas", j["nmeas"], "--nblocks", j["nblocks"], "--seed", j["seed"], "--init", "neel", "--out", j["out"], "--commit", ch]
        subprocess.run([str(c) for c in cmd], env=env, check=True); return j["out"]
    with ThreadPoolExecutor(a.workers) as ex:
        for r in ex.map(work, jobs): print("finished", os.path.basename(r), flush=True)
    if a.stage == "pilot":
        res = {}
        for j in jobs:
            tc = level_cross(load_run(j["out"])); res[j["key"]] = tc; print(j["key"], "pilot Tc ~", tc)
        os.makedirs(os.path.dirname(PIL), exist_ok=True); json.dump(res, open(PIL, "w"), indent=1)

if __name__ == "__main__": main()
