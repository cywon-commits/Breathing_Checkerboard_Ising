# Breathing checkerboard Ising, J_B<0 ordered phase: T_c check

All numbers below were produced by the code in this repository (`src/`, `analysis/`, raw data in `data/`).
Energies are in units of J_A = 1. K = |J_B|.

**Scale of this study.** The commissioning document asks for up to L=256, 24–48 replicas, 2×10⁵ measurement sweeps and
4 seeds, i.e. tens of CPU-hours per set. This session had 4 cores, so the runs are a **reduced version**
(`configs/*.yaml`, block `reduced`):

| | requested | done |
|---|---|---|
| B-lattice sizes L_b | 16, 32, 64, 128 (256 for ε=0.003) | **16, 32, 64** (S4: 16, 32) |
| sweeps per replica | 2×10⁴ equil. + 2×10⁵ meas. | 3×10³ equil. + 2×10⁴ meas. |
| seeds | 4 | 2 (1 for L_b=64 in S3/S5; 2 in S1/S2) |
| PT ladder | 24–48 replicas, T_c ±30 % | ≥20 replicas, window ±10–25 % centred on a pilot T_c (round 2) |

Consequences: the largest lattice (L=128) is **below** the requested size bound (4× kink spacing T_c/6ε in B units) for ε≤0.006 in
S1/S2 (and ε=0.003 in S3); those points are marked and are not used in the "valid" slope fit. Round 1 of the production
data (`data/round1/`) had windows that missed T_c for ε=0.003 and PT ladders too coarse for the steep Binder drop; it was
discarded and redone as round 2 (`data/S*/`, `configs/tc_refined_round2.json`). Statistical errors are jackknife over blocks; the
"drift" column is |T_x(32,64) − T_x(16,32)|, a crude measure of the unfinished finite-size convergence and the
dominant uncertainty. With only two crossings I could not extrapolate to L→∞.

## 0. Correctness checks (section 4 of the commissioning document)

`make test` (all pass, `tests/test_main.cpp`):

- **ε term**: diagonal −8ε, neighbour +4ε, 0 for n≠2 (exact, both signs of ε).
- **Excitation energies** on explicit configurations: B-mediated kink (half flip of a B plaquette) = 8K, single spin in a B plaquette 6K;
  *connected pair* 12K and *bound ring* 16K; Néel wall tension 12|ε| per unit B-length and wall (strip flip gives 2 walls).
  A-violation: the cheapest excitation above the sliding manifold is **8(J_A−|h|)** (checked for h=0, ±0.5 by exhaustive
  enumeration of all 2¹⁶ superspin configurations on a 4×4 B lattice) = a kink–antikink pair of two kinks of
  E_k = 4(J_A−|h|); an isolated kink cannot exist there (Σ_A S_A = 4Σ_B τ is a multiple of 8). *Interpretation note*: the document's
  definitions of "connected pair" and "bound ring" come from an earlier chat I did not have; I used (pair) two neighbouring
  spins from different B plaquettes inside one A plaquette and (ring) the 2×2 block formed by the facing columns of two adjacent B
  plaquettes; the minima over all sliding ground states are exactly 12K and 16K for these.
- **ΔE of each update** (single, half flip, superspin flip, strip flip, with h≠0 and sizes L_b=2,4,6) vs brute-force energy
  difference: 108 000 comparisons, 0 mismatches; incremental energy and B-column/row sums stay in sync.
- **4×4 torus (16 spins) exact enumeration** vs MC with all updates, four parameter points (ε of both signs, h≠0, T=0.9–2.5):
  ⟨e⟩, C, ⟨m_N²⟩, ⟨m_N⁴⟩, ⟨|m_N|⟩ and the m_N distribution agree within 4σ.
- **Detailed balance per update** (each alone, 4×4 torus, 2×10⁷ proposals): transition-frequency ratios
  P(s→s′)/P(s′→s) match e^{−βΔE}; worst deviation 2.3σ (single 1.9σ, half 2.3σ, super 1.6σ, strip 1.1σ).
- **Ergodicity** (`analysis/ergodic.py`, L_b=16, 24 replicas): with all updates, Néel / stripe / random starts agree
  (max 2.8σ in U, |Δ⟨m_N²⟩| ≤ 0.02) for ε=0.05 and 0.012. **Without strip flips** the ε=0.012 run started from stripes disagrees with the
  Néel start (6.4σ in U, ⟨m_N²⟩ 0.57 vs 0.70 at the lowest T): strip flips are needed for small ε. (At ε=0.05 even
  the no-strip run agrees.)
- **Lieb ice limit** (ε=0, J_B=+4): LIEB_PLACEHOLDER

Not done: no per-run autocorrelation analysis beyond the block structure (decimated series are stored in `*.series.f32`
but were not analysed); I did not verify that the autocorrelation time is <1/50 of the measurement window. Swap acceptance is ≥20 %
in all production runs except the two S5 L_b=64 runs (min 0.10–0.11).

## 1. Does a finite-T ordered phase exist for ε≠0? (P1, P2)

**Yes, for ε>0 at every ε, K, h simulated.** The Binder cumulant U_N goes from 2/3 (Néel-ordered) at low T to ≈0 at high T, and the curves for different
L cross (`results/figs/binder_S*.png`). At the lowest T, ⟨m_N²⟩ ≈ 0.97 and the fraction of diagonal-type A plaquettes is ≈1 (P1,
Néel ground state for ε>0; the order parameter is the B superspin staggered magnetisation). Stripe ground states for ε<0
(S5): U_stripe → 1/2 at low T in both S5 sets; a crossing is found for ε=−0.025 (T_c ≈ 1.32 ± 0.02 (stat) ± 0.07 (drift), pair 32/64, with
a (16,32) crossing at 1.40, so not converged) and **no crossing was found for ε=−0.012** (the S5 L_b=64 PT runs also have swap acceptance ≈0.1, i.e.
they are the least reliable). P5 (Ashkin–Teller-like behaviour) was not examined.

## 2. T_c with errors and comparison with the two estimates

Binder crossing of the largest pair (L_b = 32 and 64, L = 64 and 128) with jackknife error; "drift" = difference to the (16,32) crossing.
"leading-log" = E_k/ln(E_k/6|ε|); "improved" = solution of 12|ε| = T ln(1+e^{−4(J_A−h)/T}+e^{−4(J_A+h)/T}).

| set | K | h | ε | L_max | 4× kink spacing (B units) | T_c (Binder, largest pair) | stat ± | drift | leading-log | improved | flag |
|---|---|---|---|---|---|---|---|---|---|---|---|
| S1 | 4 | 0 | +0.003 | 64 | 164 | 0.747 | 0.004 | 0.018 | 0.740 | 1.000 | below size bound |
| S1 | 4 | 0 | +0.006 | 64 | 94 | 0.844 | 0.006 | 0.038 | 0.849 | 1.162 | below size bound |
| S1 | 4 | 0 | +0.012 | 64 | 55 | 0.963 | 0.008 | 0.008 | 0.996 | 1.379 |  |
| S1 | 4 | 0 | +0.025 | 64 | 32 | 1.195 | 0.005 | 0.019 | 1.218 | 1.707 |  |
| S1 | 4 | 0 | +0.05 | 64 | 21 | 1.563 | 0.014 | 0.004 | 1.544 | 2.174 |  |
| S2 | 8 | 0 | +0.003 | 64 | 164 | 0.755 | 0.007 | 0.050 | 0.740 | 1.000 | below size bound |
| S2 | 8 | 0 | +0.006 | 64 | 94 | 0.839 | 0.008 | 0.031 | 0.849 | 1.162 | below size bound |
| S2 | 8 | 0 | +0.012 | 64 | 55 | 0.962 | 0.009 | 0.022 | 0.996 | 1.379 |  |
| S2 | 8 | 0 | +0.025 | 64 | 32 | 1.170 | 0.006 | 0.006 | 1.218 | 1.707 |  |
| S2 | 8 | 0 | +0.05 | 64 | 21 | 1.572 | 0.005 | 0.017 | 1.544 | 2.174 |  |
| S3 | 4 | 0.5 | +0.003 | 64 | 94 | 0.566 | 0.005 | – | 0.425 | 0.685 | below size bound |
| S3 | 4 | 0.5 | +0.006 | 64 | 55 | 0.643 | 0.012 | – | 0.498 | 0.830 |  |
| S3 | 4 | 0.5 | +0.012 | 64 | 33 | 0.774 | 0.014 | 0.021 | 0.602 | 1.038 |  |
| S3 | 4 | 0.5 | +0.025 | 64 | 21 | 1.017 | 0.051 | 0.038 | 0.772 | 1.370 |  |
| S4 | 0.15 | 0.5 | +0.01 | 32 | 27 | 0.355 | 0.003 | – | 0.401 | 0.976 | B-mediated (exploratory) |
| S4 | 0.15 | 0.5 | +0.03 | 32 | 14 | 0.461 | 0.001 | – | 0.633 | 1.478 | B-mediated (exploratory) |
| S5 | | | -0.012 | | | no crossing | | | | | |
| S5 | 4 | 0 | -0.025 | 64 | 32 | 1.324 | 0.017 | 0.073 | 1.218 | 1.707 |  swap<0.2 |

Observations (S1/S2, the A-violation region, E_k=4):

- For ε≥0.012 the measured T_c agrees with the **leading-log formula** to 3–4 % (S1: 0.963/0.996, 1.195/1.218, 1.563/1.544).
  The "improved" independent-kink criterion is much *higher* than the measurement (1.38 vs 0.96 etc.), so in this
  data it is the worse of the two; I have no tested explanation. Note both formulas are asymptotic and T_c is of order J_A.
- For ε=0.003, 0.006 the lattice is below the size bound and the crossings still drift (0.76→0.75, 0.81→0.84).
- **S3 (h=0.5, E_k=2): the formula fails.** Measured T_c is 0.77 (ε=0.012) and 1.02 (ε=0.025) vs predicted 0.60 and 0.77, i.e. 28–32 % too high, and
  T_c ≳ J_A−h = 0.5 so the condition T_c ≪ J_A−h is violated: this set is outside the validity range of the estimate, not
  a clean refutation.

## 3. 1/T_c versus ln(1/ε) (P3) — `results/figs/slopes.png`

Weighted straight-line fits of 1/T_c = a + s·ln(1/ε), errors = stat ⊕ drift; "valid" excludes points below the size bound.

| set | E_k | predicted 1/E_k | fit all ε (n points) | fit valid ε only | local slopes (adjacent ε pairs, ε increasing in pair) |
|---|---|---|---|---|---|
| S1 (K=4, h=0) | 4 | 0.250 | 0.268 ± 0.010 (5), χ²/dof 5.3/3 | 0.280 ± 0.009 (3 pts, 0.012–0.05) | 0.22 ± 0.09*, 0.21 ± 0.08*, 0.275 ± 0.025, 0.284 ± 0.022 |
| S2 (K=8, h=0) | 4 | 0.250 | 0.292 ± 0.018 (5), χ²/dof 9.0/3 | 0.306 ± 0.019 (3 pts) | 0.19 ± 0.14*, 0.22 ± 0.08*, 0.25 ± 0.04, 0.32 ± 0.01 |
| S3 (K=4, h=0.5) | 2 | 0.500 | 0.350 ± 0.021 (4), χ²/dof 1.5/2 | 0.395 ± 0.044 (3 pts) | 0.31 ± 0.05*, 0.38 ± 0.07, 0.42 ± 0.10 |

(* = pair includes a point below the size bound.) The fit uncertainties underestimate the real uncertainty because the
crossings have not converged in L (χ²/dof ≫ 1 for S1/S2).

- **S1/S2 (A-violation kink region):** slopes 0.27–0.31 vs 0.25 predicted: the data are consistent with 1/E_k to within ≈10–20 %, and the local
  slopes scatter around 0.25 (0.19–0.32) without a clear trend toward 0.25 as ε→0 (the small-ε points are the least converged), so the "converges as
  ε→0" criterion can neither be confirmed nor excluded with these sizes.
- **K independence:** K=4 and K=8 slopes 0.27 ± 0.01 vs 0.29 ± 0.02 (valid-only 0.28 vs 0.31) agree within ≈1.5σ; T_c values at the same ε agree within 0.04
  (mostly within the drift). Consistent with E_k being K independent.
- **S3 (h-dependence):** the predicted slope 1/E_k=0.5 (doubling the S1 slope) is **not** reproduced: measured 0.35–0.40, and larger than
  S1's slope (the ratio S3/S1 ≈ 1.3–1.4, not 2). Given the validity violation noted above I read this as "no support, not a refutation".

## 4. Critical exponents (P4)

From the three sizes (L_b=16, 32, 64): 1/ν from d U/dT at T_c and γ/ν from χ(T_c) (log-log slope vs L_b; T_c from the largest crossing; jackknife errors,
no corrections to scaling). Values (full list in `results/fit.json`):

| set, ε | 1/ν | ν | γ/ν |
|---|---|---|---|
| S1, 0.05 | 0.76 ± 0.12 | 1.32 | 1.75 ± 0.09 |
| S2, 0.05 | 1.08 ± 0.11 | 0.92 | 1.73 ± 0.05 |
| S2, 0.025 | 1.06 ± 0.19 | 0.94 | 1.68 ± 0.08 |
| S1, 0.025 | 1.24 ± 0.09 | 0.81 | 1.63 ± 0.04 |
| S3, 0.012 / 0.025 | 1.00 ± 0.24 / 0.94 ± 0.19 | 1.0 / 1.06 | 1.59 ± 0.09 / 1.35 ± 0.05 |
| S1 and S2, ε≤0.012 | 1.4–2.8 (large errors) | – | 1.5–2.2 |

For the largest ε (0.025–0.05) the exponents are roughly consistent with 2D Ising (ν=1, γ/ν=7/4): γ/ν = 1.63–1.75 (errors 0.04–0.09, S1 ε=0.025 is 2.8σ low) and 1/ν = 1.06–1.24 (S2 ε=0.025, 0.05; S1 ε=0.025) but 0.76 ± 0.12 for S1 ε=0.05 (2σ off).
For small ε the exponents are not meaningful (below the size bound, T_c still moving). Only the large-ε results give weak, statistics-limited support for P4 (three sizes, no corrections to scaling, errors are statistical only).

## 5. Mismatches with the predictions and hypotheses

1. **h=0.5 (S3)**: T_c higher than the formula by ≈30 % and slope 0.35–0.40 instead of 0.5 (above). Hypothesis: T_c ≈ J_A−h makes the
   dilute-kink, lowest-order picture inapplicable (kinks of both signs, A violations of n_A=1 and 3 with e^{−4(J_A±h)/T} are no longer a
   small perturbation); a regime T_c ≪ J_A−h would need much smaller ε and much larger L.
2. **Violation-density check (item 5 of the analysis plan) does not verify the kink picture.** Fitting ln ρ_viol (A with n_A≠2, largest L, ±8 % around T_c) against 1/T
   gives slope −9 to −11 for S1/S2 at ε≤0.012 (−16 at the largest ε), not −E_k=−4; and ρ_viol(T_c)/e^{−E_k/T_c} is 2.5 (ε=0.003) down to 0.7
   (ε=0.05). ρ_viol counts all A plaquettes violating the ice rule, which also includes bulk superspin-flip excitations that cost 2E_k=8
   (4 violated A each), so a slope ≈ −8 to −10 is what they would give. A kink-specific observable (A violations adjacent to walls) would be
   needed to isolate e^{−E_k/T}; not done. In S3 the slope is −5 to −8 vs −E_k=−2.
3. **Improved criterion** is above the measured T_c everywhere (see above); no explanation tested.
4. **Finite-size drift**: most crossings still move by 0.01–0.05 between (16,32) and (32,64); results at ε≤0.006 should be regarded as
   preliminary. A proper test needs L_b=128–256 and the requested statistics.
5. **S5**: ε=−0.012 shows no Binder crossing; ε=−0.025 gives a T_c 8–9 % above the formula but is not converged (crossing moves from 1.40 to 1.32).

### S4 (B-mediated kink region, exploratory; K=0.15, h=0.5; formula not valid here)

Round-2 L_b=16, 32 only. U_N has a crossing at T_c = 0.355 ± 0.003 (ε=0.01) and 0.461 ± 0.001 (ε=0.03) (formula: 0.40, 0.63 — "E_k=1.2" is not applicable since T_c ≈ K
and B plaquettes are not rigid superspins). Both are an ordered Néel phase at low T with the same crossing behaviour as the other sets, but with a
single crossing and no convergence information. The ρ_viol slope is −3 to −3.5 (cf. E_k=1.2 expected; again not interpretable). Nothing here is evidence for or against the estimate.

## Figures

- `results/figs/binder_S1.png … binder_S5.png`: U(T) per size with crossings (dotted) and leading-log prediction (dashed).
- `results/figs/slopes.png`: 1/T_c vs ln(1/ε) with predicted slopes (dashed; open symbols = below size bound).
- `results/figs/snapshots.png`: B superspin snapshots at 0.7, 1.0, 1.3 T_c (S1, ε=0.025, L=128, T_c=1.2): perfect Néel order at 0.7 T_c,
  Néel with wall segments / sliding-line defects at T_c, dense line defects (disordered) at 1.3 T_c.

## Reproduce

```
make                                   # build/bcising, build/test_main
make test
python3 analysis/run_set.py pilot S1 S2 S3 S4 S5
python3 analysis/run_set.py prod  S1 S2 S3 S4 S5 --refine configs/tc_refined_round2.json
python3 analysis/analyze.py && python3 analysis/fit.py
python3 analysis/lieb.py 4 8 12 16 ; python3 analysis/ergodic.py ; python3 analysis/snapshots.py
```
Every run's `*.meta.json` stores seed, parameters, git commit hash, swap and per-update acceptance rates. The snapshot run used a binary
with the added `--dump` option (same physics code).
