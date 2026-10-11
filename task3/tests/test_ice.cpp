// Tests of the directed-loop code:  T1 scattering tables (normalisation, detailed balance), T2 ground-state Q values (+-1/2),
// T3 ice rule / flux conservation under loops, T4 exact enumeration of the 4x4-site torus (L=2) vs Monte Carlo (stag and uniform F).
#include <cstdio>
#include <cstdlib>
#include <map>
#include "../src/ice.hpp"

static int nfail = 0;
#define CHECK(c, ...) do { if (!(c)) { nfail++; printf("  FAIL: "); printf(__VA_ARGS__); printf("\n"); } } while (0)

static void t1_tables() {
    printf("T1 scattering tables\n");
    for (double c : {1.0, 1.3, 1.9999, 2.0, 2.5, 7.0, 150.0}) {
        Table t; t.build(c);
        double maxbounce = 0; int bad = 0;
        for (int S = 0; S < 16; S++) {
            if (__builtin_popcount(S) != 2) continue;
            double W = vertex_weight(S, c);
            for (int pe = 0; pe < 4; pe++) {
                double sum = 0;
                for (int q = 0; q < 4; q++) {
                    double p = t.P[S][pe][q];
                    if (p < -1e-14) bad++;
                    sum += p;
                    if (q != pe) {
                        int S2 = S ^ (1 << pe) ^ (1 << q);
                        // detailed balance: W(S) P(S;pe->q) = W(S') P(S';q->pe)
                        double lhs = W * p, rhs = vertex_weight(S2, c) * t.P[S2][q][pe];
                        if (std::fabs(lhs - rhs) > 1e-12) bad++;
                    } else maxbounce = std::max(maxbounce, p);
                }
                if (std::fabs(sum - 1) > 1e-12) bad++;
            }
        }
        printf("  c=%-8g  max bounce probability %.4f   violations %d\n", c, maxbounce, bad);
        CHECK(bad == 0, "table c=%g invalid", c);
    }
}

static void t2_ground_states() {
    printf("T2 ground-state nematic order parameter\n");
    for (int L : {4, 8}) {
        Ice M(L); Rng g(5, L); std::vector<int8_t> s;
        double q[4];
        for (int k = 0; k < 4; k++) {
            M.init_state(s, k, g);
            CHECK(M.is_ice(s.data()), "init kind %d not ice", k);
            auto o = M.measure(s.data()); q[k] = o.Q;
            CHECK(o.nc == M.NV, "init kind %d: not all A vertices c type", k);
        }
        printf("  L=%d: Q(row alt) = %+.3f, Q(col alt) = %+.3f, Q(random row phases) = %+.3f, Q(random col phases) = %+.3f\n", L, q[0], q[1], q[2], q[3]);
        CHECK(std::fabs(q[0] - 0.5) < 1e-12 && std::fabs(q[1] + 0.5) < 1e-12, "alternating ground states must give Q=+-1/2");
    }
}

static void run_mc(const Ice& M, double cA, double cB, int nsweeps, int nblocks, uint64_t seed, std::vector<std::vector<double>>& blk, int kind = 2, bool bothc = false) {
    Table tA, tB; tA.build(cA); tB.build(cB);
    Rng g(seed); std::vector<int8_t> s; M.init_state(s, kind, g);
    std::vector<int> flips;
    double avglen = 4; int nloops = std::max(1, (int)(M.N / avglen));
    long steps = 0, loops = 0;
    blk.assign(nblocks, std::vector<double>());
    const int NO = 4 + 4;  // nc, Q2, Wx2, Wy2, then 4 histogram-like extras
    for (int b = 0; b < nblocks; b++) {
        std::vector<double> acc(64, 0.0);
        int per = nsweeps / nblocks;
        for (int it = 0; it < per; it++) {
            for (int l = 0; l < nloops; l++) { int st = M.loop_update(s.data(), g, tA, tB, flips, 1000000, cB); steps += st; loops++; }
            auto o = M.measure(s.data());
            int ncount = M.count_c(s.data(), bothc);
            acc[0] += ncount; acc[1] += o.Q * o.Q; acc[2] += (double)o.Wx * o.Wx; acc[3] += (double)o.Wy * o.Wy;
            acc[4 + 0] += (o.Wx == 0 && o.Wy == 0); acc[4 + 1] += (std::abs(o.Wx) == 1 && o.Wy == 0); acc[4 + 2] += (o.Wx == 0 && std::abs(o.Wy) == 1);
            acc[4 + 3] += o.Q;
            if (b == 0 && it == 200) { avglen = (double)steps / loops; nloops = std::max(1, (int)(M.N / avglen)); }
        }
        for (auto& x : acc) x /= per;
        blk[b] = acc;
        CHECK(M.is_ice(s.data()) && M.fluxes_consistent(s.data()), "ice rule or flux conservation violated");
    }
    (void)NO;
}

static void t3_t4_exact() {
    printf("T3/T4 exact enumeration (L=2, 4x4 sites) vs directed-loop MC; ice rule and flux conservation checked at every block\n");
    Ice M(2);
    for (int model = 0; model < 2; model++) {
        for (double u : {0.0, 0.7, 1.6, 3.0}) {
            double cA = std::exp(u), cB = model == 0 ? 1.0 : std::exp(u);
            // exact sums over all ice configurations
            double Z = 0, nc = 0, q2 = 0, wx2 = 0, wy2 = 0, p00 = 0, p10 = 0, p01 = 0, q1 = 0; long nice = 0;
            std::vector<int8_t> s(16);
            for (int cfg = 0; cfg < 65536; cfg++) {
                for (int i = 0; i < 16; i++) s[i] = (cfg >> i & 1) ? 1 : -1;
                if (!M.is_ice(s.data())) continue;
                nice++;
                auto o = M.measure(s.data());
                int n = M.count_c(s.data(), model == 1);
                double w = std::exp(u * n);
                Z += w; nc += w * n; q2 += w * o.Q * o.Q; q1 += w * o.Q; wx2 += w * o.Wx * o.Wx; wy2 += w * o.Wy * o.Wy;
                p00 += w * (o.Wx == 0 && o.Wy == 0); p10 += w * (std::abs(o.Wx) == 1 && o.Wy == 0); p01 += w * (o.Wx == 0 && std::abs(o.Wy) == 1);
            }
            double ex[8] = {nc / Z, q2 / Z, wx2 / Z, wy2 / Z, p00 / Z, p10 / Z, p01 / Z, q1 / Z};
            std::vector<std::vector<double>> blk;
            run_mc(M, cA, cB, 400000, 40, 12345 + (uint64_t)(u * 100) + model, blk, 2, model == 1);
            printf("  %s u=%.1f (%ld ice states): ", model == 0 ? "stag" : "F   ", u, nice);
            const char* nm[8] = {"<nc>", "<Q^2>", "<Wx^2>", "<Wy^2>", "P(0,0)", "P(+-1,0)", "P(0,+-1)", "<Q>"};
            for (int k = 0; k < 8; k++) {
                double m = 0, e = 0; int nb = blk.size();
                for (auto& b : blk) m += b[k]; m /= nb;
                for (auto& b : blk) e += (b[k] - m) * (b[k] - m); e = std::sqrt(e / (nb * (nb - 1)));
                bool ok = std::fabs(m - ex[k]) < 4.5 * e + 1e-9;
                if (k == 7 && std::fabs(ex[k]) < 1e-9) ok = std::fabs(m) < 4.5 * e + 1e-9;
                if (!ok || k < 4) printf("%s %.4f|%.4f(%.4f)%s ", nm[k], ex[k], m, e, ok ? "" : " MISMATCH");
                CHECK(ok, "model %d u=%g %s: exact %g MC %g +- %g", model, u, nm[k], ex[k], m, e);
            }
            printf("\n");
        }
    }
}

int main() {
    t1_tables(); t2_ground_states(); t3_t4_exact();
    printf(nfail ? "\n%d FAILURES\n" : "\nALL TESTS PASSED\n", nfail);
    return nfail != 0;
}
