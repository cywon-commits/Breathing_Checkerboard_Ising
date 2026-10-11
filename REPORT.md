# Breathing checkerboard Ising, J_B<0 ordered phase: T_c check

All numbers below were produced by the code in this repository (`src/`, `analysis/`, raw data in `data/`).
Energies are in units of J_A = 1. K = |J_B|. This report covers the commissioning document including its **revision 8a**
(exact mapping to the J₁–J₂ Ising model); task 2 (transfer-matrix study of the staggered ice model) was **not started**.

**Headline.** (1) A finite-T Néel-ordered phase exists for ε>0 in every case simulated. (2) At large K the model *is* the J₁–J₂ Ising model
(verified exactly, T6, and numerically: T_c of the breathing model at K=4, 8 and of the J₁–J₂ model agree at identical lattice sizes),
so S1–S3 are a reproduction of a known model. (3) The T_c formula follows the measured T_c within 3–4 % at ε≥0.012 (h=0), the straight line
1/T_c vs ln(1/ε) has fitted slope 0.25–0.31 against the predicted 0.25 (J₁–J₂ model, ε=0.006–0.05: 0.252), but T_c is still
drifting upward with system size, the local slope falls (0.29→0.26→0.21) as ε→0 and at ε=0.003 no crossing is reached even at L_b=256, so
the slope test is **not conclusive**. (4) h=0.5 does not follow the 1/E_k=0.5 slope (0.35–0.40). (5) New S4 scan: T_c grows monotonically with K,
and the superspin picture is only partly valid at T_c (B-aligned fraction 0.65–0.95). (6) For ε<0 I find a rapidly growing specific-heat maximum
(∼L^1.5–1.6) but no clean double peak: weak first order is neither confirmed nor excluded.

## Scale of this study

The commissioning document asks for up to L=256, 24–48 replicas, 2×10⁵ measurement sweeps and 4 seeds (tens of CPU-hours per set). This session
had 4 cores, so the runs are a **reduced version** (`configs/*.yaml`):

| | requested | done |
|---|---|---|
| breathing model, B-lattice sizes L_b (L = 2L_b) | 16, 32, 64, 128 (256 for ε=0.003) | **16, 32, 64** (S4: 16, 32; 64 for ε=0.03 with poor swap acceptance) |
| J₁–J₂ model (8a cross-check, `src/j1j2.cpp`) | – | 16, 32, 64, 128, and 256 for ε=0.003 (L_b = J₁–J₂ lattice size = B lattice) |
| sweeps per replica | 2×10⁴ equil. + 2×10⁵ meas. | 3×10³ equil. + 2×10⁴ meas. (10⁵ for the ε<0 histogram runs) |
| seeds | 4 | 2 (1 for some large sizes) |
| PT ladder | 24–48 replicas, T_c ±30 % | ≥20 replicas, window ±8–25 % centred on a pilot/previous T_c |

Revision 8a item 2 (shrink S1–S3): the breathing-model runs of S1–S3 were already finished when 8a arrived; I did not run more of them and used
the J₁–J₂ model for all further/larger S1–S3 work. Round 1 of the breathing-model production data (`data/round1/`) and of the large J₁–J₂ sizes
(`data/j1j2_lowswap/`, swap acceptance 0.01–0.3) had windows/ladders that were inadequate and were discarded and redone; they are kept in the repository only
for transparency. Statistical errors are jackknife over blocks; "drift" = |T_x(largest pair) − T_x(previous pair)| is a crude measure of unfinished
finite-size convergence and the dominant uncertainty. I could not extrapolate to L→∞.

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
- **ΔE of each update** (single, half flip, superspin flip, strip flip, h≠0, L_b=2,4,6) vs brute-force energy difference: 108 000 comparisons, 0 mismatches;
  incremental energy and B-column/row sums stay in sync.
- **4×4 torus (16 spins) exact enumeration** vs MC with all updates, four parameter points (ε of both signs, h≠0, T=0.9–2.5):
  ⟨e⟩, C, ⟨m_N²⟩, ⟨m_N⁴⟩, ⟨|m_N|⟩ and the m_N distribution agree within 4σ.
- **Detailed balance per update** (each alone, 4×4 torus, 2×10⁷ proposals): transition-frequency ratios match e^{−βΔE}; worst deviation 2.3σ.
- **T6 (8a mapping)**: with all B plaquettes aligned the breathing-model energy equals E_J1J2(τ) − 8K·N_B exactly, with J₁=2(J_A+ε), J₂=J_A−2ε, field 4h (4 800 comparisons,
  0 mismatches, L_b=2,4,6, ε=0.05, −0.03, h=0.4); the strip-flip ΔE of both codes is identical; the J₁–J₂ single-spin/strip ΔE match brute force.
- **Ergodicity** (`analysis/ergodic.py`, L_b=16, 24 replicas): with all updates, Néel / stripe / random starts agree
  (max 2.8σ in U, |Δ⟨m_N²⟩| ≤ 0.02) for ε=0.05 and 0.012. **Without strip flips** the ε=0.012 run started from stripes disagrees with the
  Néel start (6.4σ in U, ⟨m_N²⟩ 0.57 vs 0.70 at the lowest T): strip flips are needed for small ε (at ε=0.05 even the no-strip run agrees).
- **Lieb ice limit** (ε=0, J_B=+4, `analysis/lieb.py`; PT ladder 0.2–60, entropy by thermodynamic integration, defect density → 0 at low T:
  A-violation density 0 at T=0.2): lowest-T entropy per spin **0.2341/0.2326 (L_b=4), 0.2192/0.2196 (8), 0.2187/0.2183 (12), 0.2168/0.2181 (16)** (random / Néel start)
  vs 0.2158 (Lieb): within 1 % at L_b≥12, decreasing with size (finite-torus excess). An earlier version of this test used T_max=8 and gave 0.242–0.26;
  the error was the high-T tail correction (checked against the exact 4×4 result) and was fixed by integrating up to T=60.

Not done: per-run autocorrelation analysis (decimated series are stored in `*.series.f32` but were not analysed); I did not verify that the autocorrelation time
is <1/50 of the measurement window, so statistical errors may be underestimated. Swap acceptance is ≥20 % in all production runs except the S4x L_b=64 runs
(0.02–0.08), the S5 L_b=64 runs (≈0.1) and the ε<0 J₁–J₂ runs where a first-order-like transition limits it (0.06–0.17).

## 1. Does a finite-T ordered phase exist for ε≠0? (P1, P2)

**Yes, for ε>0 at every ε, K, h simulated.** U_N goes from 2/3 (Néel-ordered) at low T to ≈0 (or negative, see below) at high T, and the curves for different
L cross (`results/figs/binder_S*.png`). At the lowest T, ⟨m_N²⟩ ≈ 0.97 and the fraction of diagonal-type A plaquettes is ≈1 (P1; Néel ground state for ε>0).
For ε<0 the low-T state is a stripe state (U_stripe → 1/2 at low T, S5). Binder cumulants of the Néel parameter dip **below zero** on the high-T side for small ε
(e.g. ε=0.003, J₁–J₂: min U = −0.47, −0.93, −0.96 for L_b=64, 128, 256): a non-Gaussian, bimodal-like distribution of m_N that is also
seen at first-order transitions, but which can also arise from the line-defect (sliding) liquid just above T_c; I did not resolve which.

**ε<0 (P5, revised in 8a: Ashkin–Teller prediction withdrawn; literature expects a weak first-order transition for 1/2<g<g*≈0.67).** g=J₂/J₁=0.518 (ε=−0.012) and
0.538 (ε=−0.025). Binder (stripe, two-component) crossings exist at 1.15–1.16 (ε=−0.012: (16,32) 1.163, (32,64) 1.155, (64,128) 1.154 ± 0.016; low swap acceptance 0.06–0.17 at the larger sizes) and
for ε=−0.025 1.399, 1.198 ± 0.033, 1.377 ± 0.005, i.e. not converged. Energy histograms (`analysis/j1j2_c_hist.py`, narrow window ±2 % around C_v max, 10⁵ measurements × 2 seeds,
L_b=48–128, `results/figs/c_hist.png`): the maximum of the specific heat per spin grows as

| ε | C_max (L_b=48, 64, 96, 128) | effective exponent C_max ∝ L^x (48→128) |
|---|---|---|
| −0.012 | 6.9, 11.3, 23.3, 33.9 | 1.6 |
| −0.025 | 9.8, 13.4, 26.3, 40.8 | 1.45 |

(x=2 expected for a first-order transition, ≈0 (log) for an Ising-type transition.) A double-peaked P(e) is seen only in one case (ε=−0.025, L_b=128: two peaks of separation 0.20, valley/peak 0.40, but V=−0.06 and few tunnelling events);
in the other cases the histograms are single- or noisy multi-peaked. **So the data are compatible with a weak first-order transition (growing C_max, exponent between 1 and 2) but do not
establish a double peak; longer runs (10⁶–10⁷ sweeps) and L_b≥128–256 would be needed.** I did not check the histograms for ε<0 in the breathing model itself.

## 2. T_c with errors and comparison with the two estimates (K=4, 8; breathing model)

Binder crossing of the largest pair (L_b = 32 and 64, L = 64 and 128) with jackknife error; "drift" = difference to the (16,32) crossing.
"leading-log" = E_k/ln(E_k/6|ε|); "improved" = solution of 12|ε| = T ln(1+e^{−4(J_A−h)/T}+e^{−4(J_A+h)/T}).
(Rows for S4 are in section 6; those for the J₁–J₂ model in section 3.)

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
| S5 | | | -0.012 | | | no crossing | | | | | |
| S5 | 4 | 0 | -0.025 | 64 | 32 | 1.324 | 0.017 | 0.073 | 1.218 | 1.707 |  swap<0.2 |

Observations (S1/S2, A-violation region, E_k=4):

- For ε≥0.012 the measured T_c agrees with the **leading-log formula** to 3–4 % (S1: 0.963/0.996, 1.195/1.218, 1.563/1.544), but see section 3: at larger J₁–J₂ sizes T_c at these ε is 3–4 % *higher*
  still (1.219, 1.607), so the agreement is partly accidental finite-size cancellation.
- The "improved" independent-kink criterion is much **higher** than the measurement (1.38 vs 0.96 etc.). This is in line with the remark in 8a that it is an over-estimate
  (8a: at g=0.4955 it gives T_c/J₁≈0.50, whereas here T_c/J₁ ≈ 0.37–0.42 for ε=0.003, and 0.4299 at ε=0.006, 0.49 at ε=0.012 from the J₁–J₂ runs; the 8a text quotes literature finite-size C peaks 0.380→0.320 at g=0.5; I could not access
  arXiv, so no direct comparison with Kalz et al.'s Fig. 3 was possible).
- ε=0.003 and 0.006 are below the size bound (4×kink spacing = 164, 94 B-units vs L_b ≤ 64 in the breathing model).
- **S3 (h=0.5, E_k=2): the formula fails.** Measured T_c is 0.77 (ε=0.012) and 1.02 (ε=0.025) vs predicted 0.60 and 0.77, 28–32 % too high, and
  T_c ≳ J_A−h = 0.5, so T_c ≪ J_A−h is violated; this set is outside the validity range, not a clean refutation.

## 3. J₁–J₂ cross-check (8a items 1, 2) and 1/T_c vs ln(1/ε) (P3)

### 3a. The breathing model at large K *is* the J₁–J₂ model

T_c (Binder crossing) at **identical** lattice sizes (B lattice = J₁–J₂ lattice):

| set | ε | pair (L_b) | breathing K=4 | breathing K=8 | J₁–J₂ |
|---|---|---|---|---|---|
| S1/S2 | 0.003 | (16,32) | 0.765 ± 0.008 | 0.705 ± 0.014 | 0.735 ± 0.004 |
| S1/S2 | 0.003 | (32,64) | 0.747 ± 0.004 | 0.755 ± 0.007 | 0.746 ± 0.007 |
| S1/S2 | 0.006 | (16,32) | 0.806 ± 0.004 | 0.808 ± 0.005 | 0.807 ± 0.004 |
| S1/S2 | 0.006 | (32,64) | 0.844 ± 0.006 | 0.839 ± 0.008 | 0.865 ± 0.011 |
| S1/S2 | 0.012 | (16,32) | 0.955 ± 0.003 | 0.940 ± 0.003 | 0.951 ± 0.004 |
| S1/S2 | 0.012 | (32,64) | 0.963 ± 0.008 | 0.962 ± 0.009 | 0.967 ± 0.005 |
| S1/S2 | 0.025 | (16,32) | 1.176 ± 0.004 | 1.164 ± 0.006 | 1.179 ± 0.004 |
| S1/S2 | 0.025 | (32,64) | 1.195 ± 0.005 | 1.170 ± 0.006 | 1.167 ± 0.004 |
| S1/S2 | 0.05 | (16,32) | 1.559 ± 0.006 | 1.556 ± 0.005 | 1.556 ± 0.007 |
| S1/S2 | 0.05 | (32,64) | 1.563 ± 0.014 | 1.572 ± 0.005 | 1.554 ± 0.004 |
| S3 (h=0.5) | 0.012 | (16,32) | 0.753 ± 0.014 | | 0.754 ± 0.004 |
| S3 (h=0.5) | 0.012 | (32,64) | 0.774 ± 0.014 | | 0.766 ± 0.014 |
| S3 (h=0.5) | 0.025 | (16,32) | 0.978 ± 0.006 | | 0.970 ± 0.006 |
| S3 (h=0.5) | 0.025 | (32,64) | 1.017 ± 0.051 | | 0.990 ± 0.010 |

All but a few entries agree within ~2σ (largest differences: S1 ε=0.025 (32,64) 1.195 vs 1.167, S1/S2 ε=0.006 (32,64), ≈3–4σ — plausibly statistics/autocorrelation not captured by the jackknife).
K=4 and K=8 agree, so K-independence (S2) holds as expected from the mapping (B excitations ∼e^{−6K/T}, ≈10⁻⁷ at T≈1.5, K=4). S1–S3 are therefore reproductions of the J₁–J₂ model, not independent physics.

### 3b. Larger sizes, T_c and slope (J₁–J₂ model, h=0; `results/figs/j1j2_slopes.png`)

Largest-pair crossings (swap acceptance ≥0.48 at all these sizes, round-2 ladders; the L_b=128 ladders of round 1 were discarded):

| ε | g = J₂/J₁ | T_c (largest pair) | pair | step vs previous pair | T_c/J₁ |
|---|---|---|---|---|---|
| 0.05 | 0.4286 | 1.607 ± 0.006 | (64,128) | +0.053 | 0.765 |
| 0.025 | 0.4634 | 1.219 ± 0.005 | (64,128) | +0.052 | 0.595 |
| 0.012 | 0.4822 | 0.991 ± 0.015 | (64,128) | +0.024 | 0.490 |
| 0.006 | 0.4911 | 0.865 ± 0.006 | (64,128) | 0.000 | 0.430 |
| 0.003 | 0.4955 | ≥0.83 (no crossing found; (64,128): 0.826 ± 0.011 from a wide window run, (32,64): 0.759 ± 0.002) | – | still rising | ≳0.41 |

For ε=0.003 the (128,256) curves do not cross inside the window 0.76–0.92: U₂₅₆>U₁₂₈ at all T<0.82 and both go strongly negative (U_min −0.93, −0.96), so T_c(ε=0.003) is above ≈0.83 and the
crossing was not located even at L_b=256 (the size bound 4×41=164 is nominally satisfied). The T_c(L) from crossings is still increasing with size for ε≥0.012, so all values above are lower bounds of T_c(∞) at the few-% level.

**Slope of 1/T_c vs ln(1/ε)** (largest-pair crossings, ε=0.006…0.05, J₁–J₂): straight-line fit slope **0.252** (prediction 1/E_k=0.25; intercept −0.122 vs −0.101 from the formula).
Local slopes between adjacent ε (ε decreasing): **0.286, 0.257, 0.212**; including the ε=0.003 lower bound gives a last-pair slope ≤0.08.
So the *functional form* 1/T_c ≈ a + 0.25 ln(1/ε) is a good description over 0.006≤ε≤0.05, but the local slope does not tend to 0.25 as ε→0 — it falls — while the smallest ε are the least converged
(crossings still drifting upward, which lowers 1/T_c at small ε and would lower the slope further). Not including the drift, I cannot give a trustworthy error bar; statistical errors on the slope alone are ≈0.01.

Breathing-model fits for completeness (errors = stat ⊕ drift; "valid" excludes ε below the size bound):

| set | E_k | predicted 1/E_k | fit all ε (n points) | fit valid ε only | local slopes (adjacent ε pairs, ε increasing in pair) |
|---|---|---|---|---|---|
| S1 (K=4, h=0) | 4 | 0.250 | 0.268 ± 0.010 (5), χ²/dof 5.3/3 | 0.280 ± 0.009 (3 pts, 0.012–0.05) | 0.22 ± 0.09*, 0.21 ± 0.08*, 0.275 ± 0.025, 0.284 ± 0.022 |
| S2 (K=8, h=0) | 4 | 0.250 | 0.292 ± 0.018 (5), χ²/dof 9.0/3 | 0.306 ± 0.019 (3 pts) | 0.19 ± 0.14*, 0.22 ± 0.08*, 0.25 ± 0.04, 0.32 ± 0.01 |
| S3 (K=4, h=0.5) | 2 | 0.500 | 0.350 ± 0.021 (4), χ²/dof 1.5/2 | 0.395 ± 0.044 (3 pts) | 0.31 ± 0.05*, 0.38 ± 0.07, 0.42 ± 0.10 |

(* = pair includes a point below the size bound; `results/figs/slopes.png`.) These uncertainties underestimate the real one (χ²/dof ≫ 1, no L-convergence).

- **K independence:** K=4 and K=8 slopes 0.27 ± 0.01 vs 0.29 ± 0.02 agree within ≈1.5σ (and by the exact mapping they must).
- **S3 (h-dependence):** the predicted slope 0.5 (twice S1's) is **not** reproduced: measured 0.35–0.40 (S3/S1 ≈ 1.3–1.4, not 2); with T_c ≳ J_A−h the formula's premise fails. Per 8a, the judgment is
  limited to the 1/ln functional form and the slope; the form holds for h=0, and for h=0.5 only with a different (smaller) slope.

## 4. Critical exponents (P4)

8a: the Ising universality of the Néel side is already established in the literature (ν=0.99(1) at g=0.3), so this is a consistency check only.
From the three sizes (L_b=16, 32, 64, breathing model): 1/ν from dU/dT at T_c and γ/ν from χ(T_c) (log-log slope vs L_b; T_c from the largest crossing; jackknife errors,
no corrections to scaling):

| set, ε | 1/ν | ν | γ/ν |
|---|---|---|---|
| S1, 0.05 | 0.76 ± 0.12 | 1.32 | 1.75 ± 0.09 |
| S2, 0.05 | 1.08 ± 0.11 | 0.92 | 1.73 ± 0.05 |
| S2, 0.025 | 1.06 ± 0.19 | 0.94 | 1.68 ± 0.08 |
| S1, 0.025 | 1.24 ± 0.09 | 0.81 | 1.63 ± 0.04 |
| S3, 0.012 / 0.025 | 1.00 ± 0.24 / 0.94 ± 0.19 | 1.0 / 1.06 | 1.59 ± 0.09 / 1.35 ± 0.05 |
| S1 and S2, ε≤0.012 | 1.4–2.8 (large errors) | – | 1.5–2.2 |

For the largest ε (0.025–0.05) the exponents are roughly consistent with 2D Ising (ν=1, γ/ν=7/4): γ/ν = 1.63–1.75 and 1/ν = 1.06–1.24, except S1 ε=0.05 (1/ν 2σ off).
For small ε they are not meaningful (below the size bound, T_c still moving, and the negative Binder dips indicate non-Ising-like distributions).

## 5. Mismatches with the predictions and hypotheses

1. **h=0.5 (S3)**: T_c higher than the formula by ≈30 % and slope 0.35–0.40 instead of 0.5. Hypothesis: T_c ≈ J_A−h makes the dilute-kink picture inapplicable
   (kinks of both signs, A violations of n_A=1 and 3 with e^{−4(J_A±h)/T} are no longer a small perturbation).
2. **Local slope falls as ε→0** (0.29, 0.26, 0.21, ≤0.08 in J₁–J₂): at odds with "converges to 1/E_k"; may be finite-size (the crossings are still drifting) or may be a genuine deviation from
   the leading-log form (8a already says the coefficient can be off by 30–50 %). Not decided.
3. **Violation-density check (item 5 of the analysis plan) does not verify the kink picture.** ln ρ_viol (A with n_A≠2, largest L, ±8 % around T_c) vs 1/T has slope −9 to −11 for S1/S2 at
   ε≤0.012 (−16 at the largest ε), not −E_k=−4; ρ_viol(T_c)/e^{−E_k/T_c} is 2.5 (ε=0.003) down to 0.7 (ε=0.05). ρ_viol includes bulk superspin-flip excitations that cost 2E_k=8 (4 violated A each), so a slope ≈−8 to −10 is
   what they would give. A kink-specific observable would be needed; not done. In S3 the slope is −5 to −8 vs −E_k=−2.
4. **Improved criterion** over-estimates T_c everywhere (consistent with 8a).
5. **Finite-size drift**: T_c(L) still increases with L for ε≥0.012 (J₁–J₂: +0.02…+0.05 per size doubling at L_b=64→128), so all T_c above are finite-size values; negative Binder dips above T_c
   (bimodal-like m_N distributions) complicate the crossing analysis for small ε.
6. **S5 / ε<0**: Binder crossings not converged (ε=−0.025) or hard to define (ε=−0.012); see section 1 for the histogram analysis. In the breathing model itself the S5 L_b=64 runs have swap acceptance ≈0.1.

## 6. S4 (new physics: finite K, B-mediated kink competes with the A-violation kink; 8a item 3)

`configs/task1_S4_expanded.yaml`, `analysis/run_s4.py`, `analysis/s4x_analysis.py`; figure `results/figs/s4x.png`. K = 0.1, 0.15, 0.2, 0.3; h = 0, 0.5; ε = 0.03, 0.01; L_b = 16, 32 (+64 for ε=0.03, 1 seed, poor swap acceptance 0.02–0.08, so
the (32,64) crossings† are unreliable); 2 seeds, 2×10⁴ sweeps. This region does **not** reduce to the J₁–J₂ model and the leading-log formula (E_k=min(8K, 4(J_A−h))) is outside its validity range (T_c is of order K to a few K), so
nothing below is a test of it. For h=0 all K here have 8K<4(J_A−h) (B-mediated); for h=0.5, K=0.1–0.2 are B-mediated and K=0.3 (E_k=2, 8K=2.4) on the A-violation side.

| h | ε | K | E_k = min(8K, 4(1−h)) | kink type | T_c (16,32) | T_c (32,64)† | T_c/K | B-aligned fraction at T_c | min swap acc. at largest L_b |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.03 | 0.1 | 0.8 | B-mediated | 0.396 ± 0.002 | 0.390 ± 0.004 | 3.90 | 0.67 | 0.04 (L_b=64) |
| 0 | 0.03 | 0.15 | 1.2 | B-mediated | 0.488 ± 0.002 | 0.469 ± 0.003 | 3.12 | 0.83 | 0.08 (L_b=64) |
| 0 | 0.03 | 0.2 | 1.6 | B-mediated | 0.570 ± 0.002 | 0.581 ± 0.030 | 2.91 | 0.80 | 0.05 (L_b=64) |
| 0 | 0.03 | 0.3 | 2.4 | B-mediated | 0.717 ± 0.008 | 0.729 ± 0.007 | 2.43 | 0.87 | 0.06 (L_b=64) |
| 0 | 0.01 | 0.1 | 0.8 | B-mediated | – | – | – | – | 0.45 (L_b=32) |
| 0 | 0.01 | 0.15 | 1.2 | B-mediated | – | – | – | – | 0.38 (L_b=32) |
| 0 | 0.01 | 0.2 | 1.6 | B-mediated | 0.457 ± 0.003 | – | 2.28 | 0.87 | 0.45 (L_b=32) |
| 0 | 0.01 | 0.3 | 2.4 | B-mediated | 0.569 ± 0.002 | – | 1.90 | 0.94 | 0.38 (L_b=32) |
| 0.5 | 0.03 | 0.1 | 0.8 | B-mediated | 0.364 ± 0.004 | 0.386 ± 0.001 | 3.86 | 0.65 | 0.04 (L_b=64) |
| 0.5 | 0.03 | 0.15 | 1.2 | B-mediated | 0.456 ± 0.003 | 0.482 ± 0.005 | 3.21 | 0.73 | 0.06 (L_b=64) |
| 0.5 | 0.03 | 0.2 | 1.6 | B-mediated | 0.550 ± 0.002 | 0.533 ± 0.029 | 2.67 | 0.88 | 0.02 (L_b=64) |
| 0.5 | 0.03 | 0.3 | 2.0 | A-violation | 0.681 ± 0.004 | 0.676 ± 0.002 | 2.25 | 0.90 | 0.05 (L_b=64) |
| 0.5 | 0.01 | 0.1 | 0.8 | B-mediated | 0.275 ± 0.004 | – | 2.75 | 0.79 | 0.42 (L_b=32) |
| 0.5 | 0.01 | 0.15 | 1.2 | B-mediated | 0.348 ± 0.002 | – | 2.32 | 0.88 | 0.39 (L_b=32) |
| 0.5 | 0.01 | 0.2 | 1.6 | B-mediated | 0.429 ± 0.001 | – | 2.14 | 0.89 | 0.39 (L_b=32) |
| 0.5 | 0.01 | 0.3 | 2.0 | A-violation | 0.528 ± 0.003 | – | 1.76 | 0.95 | 0.50 (L_b=32) |

† swap acceptance 0.02–0.08 at L_b=64. "–" for ε=0.01 means no Binder crossing in the window or only the (16,32) pair available (no L_b=64 for ε=0.01); for h=0, ε=0.01, K=0.1, 0.15 no crossing was found at all.

What happens:
- **T_c(K) rises monotonically with K** in all four (h, ε) families (h=0, ε=0.03: 0.40, 0.49, 0.57, 0.72 for K=0.1…0.3 from the (16,32) crossings), roughly linearly with T_c/K falling from 3.9 to 2.4 (h=0, ε=0.03; 2.75→1.76 at ε=0.01, h=0.5). For comparison,
  K→∞ (J₁–J₂, h=0, ε=0.03) is ≈1.3 (interpolating the J₁–J₂ values), so K=0.3 reaches only ≈55 % of it. h=0.5 lowers T_c by 3–8 % relative to h=0 at the same K ((16,32) crossings, ε=0.03: 0.364/0.456/0.550/0.681 vs 0.396/0.488/0.570/0.717 for K=0.1…0.3; ε=0.01: 0.429 vs 0.457 at K=0.2, 0.528 vs 0.569 at K=0.3).
  No visible discontinuity at the K≈0.25 boundary between B-mediated and A-violation kinks at h=0.5, but with errors ±0.03 at K=0.2 and unreliable L_b=64 data I cannot resolve one.
- **Superspin picture** (B-aligned fraction n_B∈{0,4} at T_c): 0.65–0.67 at K=0.1, 0.73–0.83 at K=0.15, 0.80–0.89 at K=0.2, 0.87–0.95 at K=0.3 (an uncorrelated random state would give 0.125).
  It falls steadily with T (right panel of the figure: from ≈0.9 at 0.9 T_c to 0.5–0.7 at 1.2 T_c), so for K≲0.15 about a third of the B plaquettes are not rigid at T_c: the picture of frozen superspins (hence the J₁–J₂ mapping) is
  broken there; at K=0.3 it holds to ≈10 %. I did not define a sharp "breakdown point"; the data show a smooth crossover in K.
- Ordered Néel phase exists at low T with a Binder crossing at all (K, h) where a crossing could be located; the order is of the same Néel type.

## Figures

- `results/figs/binder_S1.png … binder_S5.png`: U(T) per size with crossings (dotted) and leading-log prediction (dashed).
- `results/figs/slopes.png`, `results/figs/j1j2_slopes.png`: 1/T_c vs ln(1/ε) with the predicted slope (dashed; open symbols = below size bound).
- `results/figs/snapshots.png`: B superspin snapshots at 0.7, 1.0, 1.3 T_c (S1, ε=0.025, L=128, T_c=1.2): perfect Néel order at 0.7 T_c,
  Néel with wall segments / line defects at T_c, dense line defects (disordered) at 1.3 T_c.
- `results/figs/s4x.png`: T_c(K), B-aligned fraction at T_c, B-aligned fraction vs T/T_c. `results/figs/c_hist.png`: ε<0 energy histograms.

## Reproduce

```
make                                   # build/bcising, build/j1j2, build/test_main
make test
python3 analysis/run_set.py pilot S1 S2 S3 S4 S5
python3 analysis/run_set.py prod  S1 S2 S3 S4 S5 --refine configs/tc_refined_round2.json
python3 analysis/analyze.py && python3 analysis/fit.py
python3 analysis/run_j1j2.py A B C A2 B2 C2 A3 C3   && python3 analysis/j1j2_analysis.py && python3 analysis/j1j2_fit.py && python3 analysis/j1j2_c_hist.py
python3 analysis/run_s4.py pilot && python3 analysis/run_s4.py prod --workers 2 && python3 analysis/s4x_analysis.py
python3 analysis/lieb.py 4 8 12 16 ; python3 analysis/ergodic.py ; python3 analysis/snapshots.py
```
Every run's `*.meta.json` stores seed, parameters, git commit hash, swap and per-update acceptance rates. The snapshot run used a binary
with the added `--dump` option (same physics code). `data/round1/` and `data/j1j2_lowswap/` hold discarded first-round runs. Task 2 (transfer-matrix study) is not included.
