#!/usr/bin/env python3
"""Launch the Monte Carlo runs of one or more parameter sets (configs/S*.yaml).

  python analysis/run_set.py pilot S1 S2 S3 S4 S5      # stage 0: locate T_c roughly (Lb=32, window 0.6-2.0 x Tc_pred)
  python analysis/run_set.py prod  S1 S2 S3 S4 S5      # stage 1: production runs, PT window centred on pilot T_c
  options: --profile reduced|full  --workers 4  --dry

Every run writes data/<stage-or-set>/<tag>.{blocks.csv,series.f32,meta.json}; the meta file stores the seed,
parameters, git commit hash, swap / move acceptance rates.
"""
import argparse, json, math, os, subprocess, sys, itertools
from concurrent.futures import ThreadPoolExecutor
import yaml
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "analysis"))
from lib import load_run, jackknife, binder_fn  # noqa: E402

BIN = os.path.join(ROOT, "build", "bcising")


def tag(name, K, h, eps, Lb, seed):
    return f"{name}_K{K:g}_h{h:g}_eps{eps:g}_Lb{Lb}_s{seed}"


def nrep_for(Lb, Tmin, Tmax, minimum=6):
    step = 1.2 / Lb  # ln-T spacing that gave >=25% swap acceptance in the calibration run
    return max(minimum, int(math.ceil(math.log(Tmax / Tmin) / step)) + 1)


def pilot_tc(prefix, eps):
    r = load_run(prefix)
    kind = "N" if eps > 0 else "S"
    U, _ = jackknife(r["B"], binder_fn(kind))
    T = r["T"]
    level = 0.33
    for i in range(len(T) - 1):
        if (U[i] - level) * (U[i + 1] - level) <= 0 and U[i] > U[i + 1]:
            f = (U[i] - level) / (U[i] - U[i + 1])
            return float(T[i] + f * (T[i + 1] - T[i]))
    return None


def commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT).decode().strip()
    except Exception:
        return "unknown"


def run(cmd, env):
    subprocess.run(cmd, env=env, check=True, stderr=subprocess.PIPE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["pilot", "prod"])
    ap.add_argument("sets", nargs="+")
    ap.add_argument("--profile", default="reduced", choices=["reduced", "full"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--maxLb", type=int, default=None)
    ap.add_argument("--init", default="neel")
    a = ap.parse_args()

    ch = commit()
    jobs = []
    pil_file = os.path.join(ROOT, "data", "pilot", "tc_pilot.json")
    pil = json.load(open(pil_file)) if os.path.exists(pil_file) else {}
    for s in a.sets:
        cfg = yaml.safe_load(open(os.path.join(ROOT, "configs", f"{s}.yaml")))
        prof = cfg[a.profile]
        Lbs = prof.get("Lb", cfg["Lb"])
        if s == "S4" and a.stage == "prod":
            Lbs = [l for l in Lbs if l <= 32]  # exploratory set: small sizes only
        if a.maxLb:
            Lbs = [l for l in Lbs if l <= a.maxLb]
        for eps, tcp in zip(cfg["eps"], cfg["Tc_pred"]):
            if a.stage == "pilot":
                Lb, seed = 32, 1
                Tmin, Tmax = tcp * 0.6, tcp * 2.0
                out = os.path.join(ROOT, "data", "pilot", tag(s, cfg["K"], cfg["h"], eps, Lb, seed))
                jobs.append(dict(out=out, Lb=Lb, seed=seed, Tmin=Tmin, Tmax=Tmax, neq=1000, nmeas=4000,
                                 nstrip=4, meas_every=1, set=s, K=cfg["K"], h=cfg["h"], eps=eps, nblocks=10))
            else:
                key = f"{s}_eps{eps:g}"
                tc0 = pil.get(key) or tcp
                w = prof["window"]
                Tmin, Tmax = tc0 * (1 - w), tc0 * (1 + w)
                for Lb in Lbs:
                    for seed in range(1, prof.get("seeds_by_Lb", {}).get(Lb, prof["seeds"]) + 1):
                        out = os.path.join(ROOT, "data", s, tag(s, cfg["K"], cfg["h"], eps, Lb, seed))
                        jobs.append(dict(out=out, Lb=Lb, seed=seed, Tmin=Tmin, Tmax=Tmax, neq=prof["neq"],
                                         nmeas=prof["nmeas"], nstrip=prof["nstrip"],
                                         meas_every=1 if Lb < 32 else 2, set=s, K=cfg["K"], h=cfg["h"],
                                         eps=eps, nblocks=20))
    cost = 0
    for j in jobs:
        j["nrep"] = nrep_for(j["Lb"], j["Tmin"], j["Tmax"])
        j["cost"] = j["nrep"] * (j["neq"] + j["nmeas"]) * 4 * j["Lb"] ** 2 * 90e-9
        cost += j["cost"]
    print(f"{len(jobs)} jobs, estimated {cost/3600:.2f} CPU-hours ({cost/3600/a.workers:.2f} h wall on {a.workers} workers)")
    if a.dry:
        return
    jobs.sort(key=lambda j: -j["cost"])
    env = dict(os.environ, OMP_NUM_THREADS="1")

    def work(j):
        os.makedirs(os.path.dirname(j["out"]), exist_ok=True)
        if os.path.exists(j["out"] + ".meta.json"):
            return j["out"] + " (cached)"
        cmd = [BIN, "--Lb", j["Lb"], "--K", j["K"], "--h", j["h"], "--eps", j["eps"], "--Tmin", j["Tmin"],
               "--Tmax", j["Tmax"], "--nrep", j["nrep"], "--neq", j["neq"], "--nmeas", j["nmeas"],
               "--nblocks", j["nblocks"], "--meas_every", j["meas_every"], "--nstrip", j["nstrip"],
               "--seed", j["seed"], "--init", a.init, "--out", j["out"], "--commit", ch]
        subprocess.run([str(c) for c in cmd], env=env, check=True)
        return j["out"]

    with ThreadPoolExecutor(a.workers) as ex:
        for r in ex.map(work, jobs):
            print("finished", os.path.basename(r), flush=True)
    if a.stage == "pilot":
        res = dict(pil)
        for j in jobs:
            tc = pilot_tc(j["out"], j["eps"])
            res[f"{j['set']}_eps{j['eps']:g}"] = tc
            print(j["set"], j["eps"], "pilot Tc ~", tc)
        json.dump(res, open(pil_file, "w"), indent=1)


if __name__ == "__main__":
    main()
