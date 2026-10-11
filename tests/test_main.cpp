// Unit tests for the breathing checkerboard Ising model.
//   T1  eps term: diagonal -8 eps, neighbour +4 eps, zero for n != 2
//   T2  explicit-configuration excitation energies (see comments for the interpretation used)
//   T3  incremental dE of every update == brute-force energy difference (random states, several sizes)
//   T4  4x4 torus (16 spins) exact enumeration vs. Monte Carlo with all updates
//   T5  detailed balance per update: measured transition-frequency ratios vs Boltzmann ratios
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

#include "../src/model.hpp"

static int nfail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { nfail++; printf("  FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)
static bool near(double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol * (1 + std::fabs(b)); }

// ---------------------------------------------------------------- T1
static void test_eps_term() {
    printf("T1 eps term\n");
    for (double eps : {0.05, -0.03}) {
        Model M(2, 1.0, eps, 4.0, 0.0);
        for (int c = 0; c < 16; c++) {
            int s[4]; int n = 0;
            for (int p = 0; p < 4; p++) { s[p] = (c >> p & 1) ? 1 : -1; n += s[p] == 1; }
            int S = s[0] + s[1] + s[2] + s[3];
            double rest = M.JA * (S * S - 4) / 2.0;  // J_A part depends on S only
            double extra = M.EA4(s[0], s[1], s[2], s[3]) - rest;
            double expect = 0;
            if (n == 2) expect = (s[0] == s[3]) ? -8 * eps : 4 * eps;  // diagonal-type : neighbour-type
            CHECK(near(extra, expect), "eps=%g config %d: got %g expected %g", eps, c, extra, expect);
        }
    }
}

// ---------------------------------------------------------------- helpers
static void random_state(const Model& M, Replica& r, uint64_t seed) { M.init(r, 3, seed, 7); }
static double brute_dE_flip(const Model& M, Replica& r, const std::vector<int>& sites) {
    double e0 = M.total_energy(r.s.data());
    for (int i : sites) r.s[i] = -r.s[i];
    double e1 = M.total_energy(r.s.data());
    for (int i : sites) r.s[i] = -r.s[i];
    return e1 - e0;
}
// sliding ground states: B superspin configuration (+-1 per B cell) -> spins
static void set_superspins(const Model& M, Replica& r, const std::vector<int>& tau) {
    r.s.assign(M.N, 1);
    for (int b = 0; b < M.NB; b++)
        for (int p = 0; p < 4; p++) r.s[M.Bsite[4 * b + p]] = (int8_t)tau[b];
    M.rebuild(r);
}

// ---------------------------------------------------------------- T2
static void test_excitations() {
    printf("T2 excitation energies\n");
    // (a) B-mediated kink = half flip of a ferromagnetic B plaquette = 8K ; single spin in B = 6K
    {
        double K = 3.0;
        Model M(4, 1.0, 0.0, K, 0.0);
        Replica r; std::vector<int> tau(M.NB, 1); set_superspins(M, r, tau);
        int st[2] = {M.Bsite[0], M.Bsite[1]};
        double dEB = 0;  // B part only: compute from EB
        { double e0 = M.EB(r.s.data(), 0); r.s[st[0]] = -1; r.s[st[1]] = -1; dEB = M.EB(r.s.data(), 0) - e0; r.s[st[0]] = 1; r.s[st[1]] = 1; }
        CHECK(near(dEB, 8 * K), "B-mediated kink B-energy %g != 8K=%g", dEB, 8 * K);
        double e0 = M.EB(r.s.data(), 0); r.s[st[0]] = -1; double d1 = M.EB(r.s.data(), 0) - e0; r.s[st[0]] = 1;
        CHECK(near(d1, 6 * K), "single spin in B: %g != 6K", d1);
    }
    // (b) A-violation kink: lowest excitation above the sliding manifold for K large is 4(J_A -+ h)
    //     (exhaustive over all 2^16 B-superspin configurations of a 4x4 B lattice).
    for (double h : {0.0, 0.5, -0.5}) {
        Model M(4, 1.0, 0.0, 50.0, h);
        Replica r; std::vector<int> tau(M.NB);
        double Emin = 1e30;
        std::vector<double> E(1 << 16);
        for (int c = 0; c < (1 << 16); c++) {
            for (int b = 0; b < M.NB; b++) tau[b] = (c >> b & 1) ? 1 : -1;
            set_superspins(M, r, tau); E[c] = r.E; Emin = std::min(Emin, E[c]);
        }
        // sector with zero net violations at h=0: Emin; the manifold with all n_A=2 has magnetisation 0.
        // find lowest energy among configs having at least one n_A != 2 (relative to the n_A=2 manifold)
        double Eground = 1e30, Eexc = 1e30;
        for (int c = 0; c < (1 << 16); c++) {
            for (int b = 0; b < M.NB; b++) tau[b] = (c >> b & 1) ? 1 : -1;
            set_superspins(M, r, tau);
            Obs o = M.measure(r);
            if (o.rviol == 0) Eground = std::min(Eground, E[c]); else Eexc = std::min(Eexc, E[c]);
        }
        double gap = Eexc - Eground;
        // A single superspin flip violates its four A plaquettes (n_A = 2 -> 1 or 3, 2J_A -+ 2h each):
        // 8(J_A - |h|) for the favourable sign = a kink-antikink pair of two kinks of E_k = 4(J_A - |h|).
        // (A lone kink cannot exist in this 4x4 enumeration: sum_A S_A = 4 sum_B tau is a multiple of 8.)
        double expect = 8 * (M.JA - std::fabs(h));
        CHECK(near(gap, expect), "h=%g: lowest A-violation excitation %g, expected 2*E_k=8(J_A-|h|)=%g", h, gap, expect);
        printf("  h=%+.1f: lowest A-violation excitation %.6f  (= 2 E_k = 8(J_A-|h|) = %.6f)\n", h, gap, expect);
    }
    // (c) 'connected pair' (12K) and 'bound ring' (16K): minimum over placements and sliding ground
    //     states (eps=0, h=0) of the energy change when flipping
    //       pair: two neighbouring spins that sit in different B plaquettes (same A plaquette),
    //       ring: the 2x2 block made of the two facing columns of two adjacent B plaquettes.
    {
        double K = 20.0;  // large so that the pure-A part is a small correction to the K multiple
        Model M(4, 1.0, 0.0, K, 0.0);
        Replica r; std::vector<int> tau(M.NB);
        double minPair = 1e30, minRing = 1e30;
        for (int c = 0; c < (1 << 16); c++) {
            for (int b = 0; b < M.NB; b++) tau[b] = (c >> b & 1) ? 1 : -1;
            set_superspins(M, r, tau);
            if (M.measure(r).rviol != 0) continue;  // sliding ground states only
            for (int Y = 0; Y < M.Lb; Y++)
                for (int X = 0; X < M.Lb; X++) {
                    int a = Y * M.Lb + X;
                    // pair across a B boundary: A positions (0,1) horizontal and (0,2) vertical
                    for (auto pr : {std::pair<int,int>{0, 1}, {0, 2}, {2, 3}, {1, 3}}) {
                        std::vector<int> st = {M.Asite[4 * a + pr.first], M.Asite[4 * a + pr.second]};
                        minPair = std::min(minPair, brute_dE_flip(M, r, st));
                    }
                    // ring: two A cells stacked / adjacent share the 2x2 block -> use A-cell origin shifted
                    // block of 4 spins = sites (2X+2..2X+3, 2Y+1..2Y+2) (spans two A cells, two B cells) horizontally
                    {
                        std::vector<int> st = {M.idx(2 * X + 2, 2 * Y + 1), M.idx(2 * X + 3, 2 * Y + 1),
                                               M.idx(2 * X + 2, 2 * Y + 2), M.idx(2 * X + 3, 2 * Y + 2)};
                        minRing = std::min(minRing, brute_dE_flip(M, r, st));
                    }
                    {
                        std::vector<int> st = {M.idx(2 * X + 1, 2 * Y + 2), M.idx(2 * X + 2, 2 * Y + 2),
                                               M.idx(2 * X + 1, 2 * Y + 3), M.idx(2 * X + 2, 2 * Y + 3)};
                        minRing = std::min(minRing, brute_dE_flip(M, r, st));
                    }
                }
        }
        // B part is exactly 12K / 16K; the A part is >= 0 and vanishes if the A plaquettes keep n_A=2
        printf("  min pair dE = %.4f (12K = %.1f), min ring dE = %.4f (16K = %.1f)\n", minPair, 12 * K, minRing, 16 * K);
        CHECK(near(minPair, 12 * K, 1e-9), "connected pair minimum %g != 12K", minPair);
        CHECK(near(minRing, 16 * K, 1e-9), "bound ring minimum %g != 16K", minRing);
    }
    // (d) Neel stripe defect: tension 12|eps| per unit B-length and wall (two walls for a strip flip)
    for (double eps : {0.0125, -0.0125}) {
        int Lb = 8;
        Model M(Lb, 1.0, eps, 8.0, 0.0);
        Replica r;
        // ground state: Neel (eps>0) or X stripes (eps<0)
        M.init(r, eps > 0 ? 0 : 1, 3, 0);
        double dE = M.dE_strip(r, eps > 0 ? 0 : 1, 2, 3);  // columns 2..4 ; for stripes flip rows (parallel to stripes -> wall costs)
        // for eps>0 (Neel): flipping a column strip creates two straight walls
        // for eps<0 (X stripes, i.e. alternating columns): flipping a *row* strip creates 2 walls of the other kind
        double per_len = dE / (2.0 * Lb);
        printf("  eps=%+.4f: strip flip dE = %.6f -> %.6f per unit length per wall (expect 12|eps| = %.6f)\n", eps, dE, per_len, 12 * std::fabs(eps));
        if (eps > 0) CHECK(near(per_len, 12 * std::fabs(eps)), "Neel wall tension %g != 12 eps", per_len);
        // brute-force check of the strip dE in both cases
        Replica r2 = r; int ori = eps > 0 ? 0 : 1; int x0 = 2, w = 3;
        std::vector<int> sites;
        for (int k = 0; k < w; k++) for (int t = 0; t < Lb; t++) {
            int b = ori == 0 ? t * Lb + (x0 + k) : (x0 + k) * Lb + t;
            for (int p = 0; p < 4; p++) sites.push_back(M.Bsite[4 * b + p]);
        }
        CHECK(near(brute_dE_flip(M, r2, sites), dE), "strip dE != brute force");
    }
    // Neel wall as a straight wall of a *row* strip in the stripe state would be a defect of other type:
    // verified implicitly by the brute-force comparison in T3.
}

// ---------------------------------------------------------------- T3
static void test_dE_vs_bruteforce() {
    printf("T3 dE of every update vs brute force\n");
    int bad = 0, total = 0;
    for (int Lb : {2, 4, 6}) {
        for (double h : {0.0, 0.37}) {
            Model M(Lb, 1.0, 0.07, 1.3, h);
            Replica r; random_state(M, r, 11 + Lb);
            Rng g(5, Lb);
            for (int trial = 0; trial < 3000; trial++) {
                // random walk the state a bit
                int i = g.below(M.N);
                // single spin
                { double d = M.dE_single(r, i); double b = brute_dE_flip(M, r, {i}); total++; if (!near(d, b)) bad++; }
                // half flip
                { int b = g.below(M.NB), e = g.below(4);
                  int st[2] = {M.Bsite[4 * b + Model::HALF_P[e][0]], M.Bsite[4 * b + Model::HALF_P[e][1]]};
                  double d = M.dE_flipset(r, st, 2); double bb = brute_dE_flip(M, r, {st[0], st[1]}); total++; if (!near(d, bb)) bad++; }
                // superspin flip
                { int b = g.below(M.NB); const int* st = &M.Bsite[4 * b];
                  double d = M.dE_flipset(r, st, 4); double bb = brute_dE_flip(M, r, {st[0], st[1], st[2], st[3]}); total++; if (!near(d, bb)) bad++; }
                // strip (all w, both orientations)
                { int ori = g.below(2), x0 = g.below(Lb), w = 1 + g.below(Lb - 1);
                  double d = M.dE_strip(r, ori, x0, w);
                  std::vector<int> sites;
                  for (int k = 0; k < w; k++) for (int t = 0; t < Lb; t++) {
                      int b = ori == 0 ? t * Lb + (x0 + k) % Lb : ((x0 + k) % Lb) * Lb + t;
                      for (int p = 0; p < 4; p++) sites.push_back(M.Bsite[4 * b + p]);
                  }
                  double bb = brute_dE_flip(M, r, sites); total++; if (!near(d, bb)) bad++; }
                // advance state with accepted random moves (keeps colB/rowB bookkeeping exercised)
                { double d = M.dE_single(r, i); M.flip(r, i); r.E += d; }
                if (trial % 500 == 0) {
                    Replica chk = r; M.rebuild(chk);
                    CHECK(near(chk.E, r.E, 1e-8), "tracked energy drift");
                    CHECK(chk.colB == r.colB && chk.rowB == r.rowB, "column/row sums out of sync");
                }
            }
        }
    }
    CHECK(bad == 0, "%d of %d dE mismatches", bad, total);
    printf("  %d comparisons, %d mismatches\n", total, bad);
}

// ---------------------------------------------------------------- T4
static void test_enumeration() {
    printf("T4 4x4 torus exact enumeration vs MC\n");
    struct Cfg { double eps, K, h, T; } cfgs[] = {{0.15, 0.8, 0.2, 1.0}, {-0.1, 1.5, 0.0, 0.9}, {0.3, 0.6, 0.3, 2.5}, {-0.2, 0.4, 0.1, 2.0}};
    for (auto c : cfgs) {
        Model M(2, 1.0, c.eps, c.K, c.h);
        double beta = 1.0 / c.T;
        Replica r;
        // exact sums
        double Z = 0, sE = 0, sE2 = 0, sm2 = 0, sm4 = 0, sabs = 0;
        std::vector<double> hist(9, 0.0);  // distribution of mN in 9 bins of width 0.25 over [-1.125,1.125]
        std::vector<int> st(16);
        // first pass for Emin
        double Emin = 1e30;
        for (int cfg = 0; cfg < 65536; cfg++) {
            r.s.assign(16, 1); for (int i = 0; i < 16; i++) r.s[i] = (cfg >> i & 1) ? 1 : -1;
            Emin = std::min(Emin, M.total_energy(r.s.data()));
        }
        for (int cfg = 0; cfg < 65536; cfg++) {
            r.s.assign(16, 1); for (int i = 0; i < 16; i++) r.s[i] = (cfg >> i & 1) ? 1 : -1;
            M.rebuild(r);
            double w = std::exp(-beta * (r.E - Emin));
            Obs o = M.measure(r);
            Z += w; sE += w * o.E; sE2 += w * o.E * o.E; sm2 += w * o.mN * o.mN; sm4 += w * std::pow(o.mN, 4); sabs += w * std::fabs(o.mN);
            int bin = (int)std::floor((o.mN + 1.125) / 0.25); hist[bin] += w;
        }
        double eE = sE / Z, eE2 = sE2 / Z, em2 = sm2 / Z, em4 = sm4 / Z, eabs = sabs / Z;
        double eC = beta * beta * 16 * (eE2 - eE * eE);  // specific heat per spin: beta^2 N var(e)
        for (auto& x : hist) x /= Z;
        // Monte Carlo: nb independent blocks
        const int nb = 40; const long sweeps = 40000;
        std::vector<double> bE(nb), bE2(nb), bm2(nb), bm4(nb), babs(nb); std::vector<std::vector<double>> bh(nb, std::vector<double>(9));
        M.init(r, 3, 17, 3);
        long tr[N_MOVES] = {0}, ac[N_MOVES] = {0};
        for (int i = 0; i < 2000; i++) M.sweep(r, beta, 2, tr, ac);
        for (int b = 0; b < nb; b++) {
            for (long k = 0; k < sweeps; k++) {
                M.sweep(r, beta, 2, tr, ac);
                Obs o = M.measure(r);
                bE[b] += o.E; bE2[b] += o.E * o.E; bm2[b] += o.mN * o.mN; bm4[b] += std::pow(o.mN, 4); babs[b] += std::fabs(o.mN);
                bh[b][(int)std::floor((o.mN + 1.125) / 0.25)] += 1;
            }
            for (double* p : {&bE[b], &bE2[b], &bm2[b], &bm4[b], &babs[b]}) *p /= sweeps;
            for (auto& x : bh[b]) x /= sweeps;
        }
        auto stat = [&](std::vector<double>& v, double& mean, double& err) {
            mean = 0; for (double x : v) mean += x; mean /= v.size();
            double s2 = 0; for (double x : v) s2 += (x - mean) * (x - mean);
            err = std::sqrt(s2 / (v.size() * (v.size() - 1)));
        };
        printf("  eps=%+.2f K=%.1f h=%.1f T=%.1f   (exact | MC +- err)\n", c.eps, c.K, c.h, c.T);
        auto cmp = [&](const char* nm, double ex, std::vector<double>& v) {
            double m, e; stat(v, m, e);
            bool ok = std::fabs(m - ex) < 4 * e + 1e-12;
            printf("    %-8s %.6f | %.6f +- %.6f %s\n", nm, ex, m, e, ok ? "ok" : "MISMATCH");
            CHECK(ok, "%s: exact %g MC %g +- %g", nm, ex, m, e);
        };
        cmp("<e>", eE, bE); cmp("<mN^2>", em2, bm2); cmp("<mN^4>", em4, bm4); cmp("<|mN|>", eabs, babs);
        {   // specific heat from block variance of e: use pooled estimate  C = beta^2 N (<e^2>-<e>^2)
            std::vector<double> bC(nb);
            for (int b = 0; b < nb; b++) bC[b] = beta * beta * 16 * (bE2[b] - bE[b] * bE[b]);
            cmp("C/spin", eC, bC);
        }
        for (int k = 0; k < 9; k++) {
            if (hist[k] < 1e-3) continue;
            std::vector<double> v(nb); for (int b = 0; b < nb; b++) v[b] = bh[b][k];
            char nm[32]; snprintf(nm, 32, "P(mN bin%d)", k);
            cmp(nm, hist[k], v);
        }
    }
}

// ---------------------------------------------------------------- T5
static void test_detailed_balance() {
    printf("T5 detailed balance (transition-frequency ratios) per update\n");
    // Each update alone acts on the 4x4 torus; we count transitions s -> s' (s != s').
    // Detailed balance with symmetric proposals: [n(s->s')/n(s)] / [n(s'->s)/n(s')] = exp(-beta dE).
    struct Mv { const char* name; unsigned bit; } mvs[] = {{"single", 1}, {"half", 2}, {"super", 4}, {"strip", 8}};
    Model M(2, 1.0, 0.1, 0.5, 0.15);
    double beta = 1.0 / 2.0;
    for (auto mv : mvs) {
        Replica r; M.init(r, 3, 99, mv.bit);
        std::map<std::pair<uint32_t, uint32_t>, long> trans; std::map<uint32_t, long> visits;
        auto code = [&](const Replica& q) { uint32_t c = 0; for (int i = 0; i < 16; i++) c |= (q.s[i] > 0) << i; return c; };
        long tr[N_MOVES] = {0}, ac[N_MOVES] = {0};
        uint32_t cur = code(r);
        const long steps = 20000000;
        for (long k = 0; k < steps; k++) {
            M.sweep(r, beta, 4, tr, ac, mv.bit);  // a "sweep" of only this move type
            // (the sweep applies many proposals; track state before/after is too coarse) -> use single-proposal mode below
            (void)k;
            break;
        }
        // single-proposal mode: emulate one proposal of the chosen move at a time
        Rng g(4242 + mv.bit);
        cur = code(r);
        for (long k = 0; k < steps; k++) {
            visits[cur]++;
            double dE = 0; bool acc;
            int st[4], n = 0;
            if (mv.bit == 1) { int i = g.below(16); st[0] = i; n = 1; dE = M.dE_single(r, i); }
            else if (mv.bit == 2) { int b = g.below(4), e = g.below(4); st[0] = M.Bsite[4 * b + Model::HALF_P[e][0]]; st[1] = M.Bsite[4 * b + Model::HALF_P[e][1]]; n = 2; dE = M.dE_flipset(r, st, 2); }
            else if (mv.bit == 4) { int b = g.below(4); for (int p = 0; p < 4; p++) st[p] = M.Bsite[4 * b + p]; n = 4; dE = M.dE_flipset(r, st, 4); }
            acc = Model::accept(dE, beta, g);
            if (mv.bit == 8) {
                int ori = g.below(2), x0 = g.below(2), w = 1 + g.below(1);
                double d = M.dE_strip(r, ori, x0, w);
                if (Model::accept(d, beta, g)) { M.apply_strip(r, ori, x0, w, d); uint32_t nx = code(r); trans[{cur, nx}]++; cur = nx; }
                continue;
            }
            if (acc) {
                for (int j = 0; j < n; j++) M.flip(r, st[j]);
                uint32_t nx = code(r); trans[{cur, nx}]++; cur = nx;
            }
        }
        int checked = 0, bad = 0; double worst = 0;
        for (auto& kv : trans) {
            uint32_t s = kv.first.first, t = kv.first.second;
            auto rev = trans.find({t, s});
            if (rev == trans.end() || s > t) continue;
            long n1 = kv.second, n2 = rev->second;
            if (n1 < 1000 || n2 < 1000) continue;
            Replica a, b; a.s.assign(16, 1); b.s.assign(16, 1);
            for (int i = 0; i < 16; i++) { a.s[i] = (s >> i & 1) ? 1 : -1; b.s[i] = (t >> i & 1) ? 1 : -1; }
            double Ea = M.total_energy(a.s.data()), Eb = M.total_energy(b.s.data());
            double ratio = (double)n1 / visits[s] / ((double)n2 / visits[t]);   // P(s->t)/P(t->s)
            double expect = std::exp(-beta * (Eb - Ea));
            double rel_err = 1.0 / std::sqrt((double)n1) + 1.0 / std::sqrt((double)n2);
            double dev = std::fabs(ratio / expect - 1) / rel_err;
            worst = std::max(worst, dev); checked++;
            if (dev > 5) bad++;
        }
        printf("  %-6s: %d state pairs checked, worst deviation %.2f sigma, %d beyond 5 sigma\n", mv.name, checked, worst, bad);
        CHECK(checked > 0 && bad == 0, "%s: detailed balance violated (%d bad pairs of %d)", mv.name, bad, checked);
    }
}

int main(int argc, char** argv) {
    std::string only = argc > 1 ? argv[1] : "";
    if (only.empty() || only == "T1") test_eps_term();
    if (only.empty() || only == "T2") test_excitations();
    if (only.empty() || only == "T3") test_dE_vs_bruteforce();
    if (only.empty() || only == "T4") test_enumeration();
    if (only.empty() || only == "T5") test_detailed_balance();
    printf(nfail ? "\n%d FAILURES\n" : "\nALL TESTS PASSED\n", nfail);
    return nfail != 0;
}
