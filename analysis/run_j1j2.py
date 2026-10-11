#!/usr/bin/env python3
"""Runs of the plain J1-J2 Ising model (src/j1j2.cpp) that corresponds exactly to the breathing checkerboard model
with frozen B superspins (K -> infinity): J1 = 2(J_A+eps), J2 = J_A-2eps, field 4h.   (task-1 revision 8a)

  A: h=0,   eps = 0.003 ... 0.05 (S1/S2 cross-check, sizes up to L_b=128, 256 for eps=0.003)
  B: h=0.5, eps = 0.003 ... 0.025 (S3 cross-check)
  C: h=0,   eps = -0.012, -0.025  (S5: g = J2/J1 = 0.518, 0.538; weak first order?)
  python analysis/run_j1j2.py [A B C] [--workers 4] [--dry]
"""
import argparse, json, math, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "build", "j1j2")

GROUPS = {
    "A": dict(h=0.0, eps={0.003: (0.75, 0.20), 0.006: (0.84, 0.12), 0.012: (0.96, 0.10), 0.025: (1.19, 0.10), 0.05: (1.56, 0.10)},
              Lb=[16, 32, 64, 128], extra={0.003: [256]}),
    "B": dict(h=0.5, eps={0.003: (0.60, 0.22), 0.006: (0.66, 0.15), 0.012: (0.78, 0.12), 0.025: (1.00, 0.12)},
              Lb=[16, 32, 64, 128], extra={}),
    "C": dict(h=0.0, eps={-0.012: (1.16, 0.15), -0.025: (1.33, 0.12)}, Lb=[16, 32, 64, 128], extra={}),
    # round 2 for the large sizes (round-1 L_b = 128/256 had swap acceptance 0.01-0.3 -> not equilibrated; kept in data/j1j2_lowswap)
    "A2": dict(h=0.0, eps={0.003: (0.746, 0.08), 0.006: (0.84, 0.08), 0.012: (0.97, 0.08), 0.025: (1.19, 0.08), 0.05: (1.57, 0.08)},
               Lb=[128], extra={0.003: [256]}, tag="A", step_coef=0.67),
    "B2": dict(h=0.5, eps={0.003: (0.60, 0.10), 0.006: (0.66, 0.08), 0.012: (0.76, 0.08), 0.025: (1.00, 0.08)},
               Lb=[128], extra={}, tag="B", step_coef=0.67),
    # round 3: eps=0.003 transition lies above the round-2 window (U_256 still 0.5 at the upper edge) -> window centred on 0.84
    "A3": dict(h=0.0, eps={0.003: (0.84, 0.10)}, Lb=[32, 64, 128], extra={0.003: [256]}, tag="A3", step_coef=0.67),
    # round 3 for eps<0: narrow windows around the C_v maximum, long runs, series_every=1 (energy histograms, P5 / 8a)
    "C3": dict(h=0.0, eps={-0.012: (1.162, 0.02), -0.025: (1.395, 0.02)}, Lb=[48, 64, 96, 128], extra={}, tag="C3", nrep=12, nmeas=100000),
    "C2": dict(h=0.0, eps={-0.012: (1.16, 0.08), -0.025: (1.33, 0.08)}, Lb=[128], extra={}, tag="C", step_coef=0.67),
}


def commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT).decode().strip()
    except Exception:
        return "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("groups", nargs="*", default=["A", "B", "C"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    jobs = []
    for gname in a.groups:
        G = GROUPS[gname]
        for eps, (tc, w) in G["eps"].items():
            for Lb in G["Lb"] + G["extra"].get(eps, []):
                seeds = 2 if Lb <= 128 and not (gname[0] in "BC" and Lb == 128) else 1
                for seed in range(1, seeds + 1):
                    coef = G.get("step_coef", 2.4)
                    nrep = max(20, int(math.ceil(math.log((1 + w) / (1 - w)) / (coef / Lb))) + 1)
                    nmeas = G.get("nmeas", 10000 if Lb >= 256 else 20000)
                    if "nrep" in G: nrep = G["nrep"]
                    out = os.path.join(ROOT, "data", "j1j2", f"j1j2_{G.get('tag', gname)}_h{G['h']:g}_eps{eps:g}_Lb{Lb}_s{seed}")
                    jobs.append(dict(out=out, Lb=Lb, eps=eps, h=G["h"], seed=seed, Tmin=tc * (1 - w), Tmax=tc * (1 + w), nrep=nrep,
                                     neq=3000, nmeas=nmeas, series_every=1 if gname[0] == "C" else 10,
                                     cost=nrep * (3000 + nmeas) * Lb * Lb * 55e-9))
    jobs = [j for j in jobs]
    tot = sum(j["cost"] for j in jobs)
    print(f"{len(jobs)} jobs, ~{tot/3600:.2f} CPU-h")
    if a.dry:
        return
    jobs.sort(key=lambda j: -j["cost"])
    ch = commit()

    def work(j):
        os.makedirs(os.path.dirname(j["out"]), exist_ok=True)
        if os.path.exists(j["out"] + ".meta.json"):
            return j["out"] + " (cached)"
        cmd = [BIN, "--Lb", j["Lb"], "--eps", j["eps"], "--h", j["h"], "--Tmin", j["Tmin"], "--Tmax", j["Tmax"], "--nrep", j["nrep"],
               "--neq", j["neq"], "--nmeas", j["nmeas"], "--nblocks", 20, "--seed", j["seed"], "--series_every", j["series_every"],
               "--out", j["out"], "--commit", ch]
        subprocess.run([str(c) for c in cmd], check=True)
        return j["out"]

    with ThreadPoolExecutor(a.workers) as ex:
        for r in ex.map(work, jobs):
            print("finished", os.path.basename(r), flush=True)


if __name__ == "__main__":
    main()
