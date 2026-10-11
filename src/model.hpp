// Breathing checkerboard Ising model: geometry, energies, Monte Carlo updates.
//
// H = sum_A [J1 sum_edges s_i s_j + J2 sum_diag s_i s_j] + sum_B (J_B/2) S_B^2 - h sum_i s_i
// with J1 = J_A + eps, J2 = J_A - 2 eps, J_B = -K.
//
// Geometry (L x L torus, L % 4 == 0, Lb = L/2 even):
//   B cell (X,Y) has origin site (2X+1, 2Y+1), A cell (X,Y) has origin site (2X+2, 2Y+2).
//   Local position p = dx + 2*dy of a site inside a plaquette (dx,dy in {0,1}).
//   The four spins of A cell (X,Y) belong to the four B cells (X+dx, Y+dy) (same p).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

struct Rng {
    uint64_t s[4];
    static uint64_t splitmix(uint64_t& x) {
        uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
    explicit Rng(uint64_t seed = 1, uint64_t stream = 0) {
        uint64_t x = seed * 0x2545F4914F6CDD1DULL + stream * 0x9E3779B97F4A7C15ULL + 12345;
        for (auto& v : s) v = splitmix(x);
    }
    static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    inline uint64_t next() {
        uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
        return r;
    }
    inline double uniform() { return (next() >> 11) * (1.0 / 9007199254740992.0); }
    inline int below(int n) { return (int)(((unsigned __int128)next() * (unsigned)n) >> 64); }
};

enum Move { MV_SINGLE = 0, MV_HALF = 1, MV_SUPER = 2, MV_STRIP = 3, N_MOVES = 4 };

struct Replica {
    std::vector<int8_t> s;
    double E = 0;
    std::vector<int> colB, rowB;  // sum of spins in B columns / B rows (for the strip field energy)
    Rng rng;
};

struct Obs {
    double E, mN, mx, my, rdiag, rnbr, rviol, rBal;
};

struct Model {
    int L, Lb, N, NB, NA;
    double JA, eps, K, h, J1, J2, JB;
    std::vector<int> Asite, Bsite, Aof, Bof;

    Model(int Lb_, double JA_, double eps_, double K_, double h_)
        : Lb(Lb_), JA(JA_), eps(eps_), K(K_), h(h_) {
        L = 2 * Lb; N = L * L; NB = NA = Lb * Lb;
        J1 = JA + eps; J2 = JA - 2 * eps; JB = -K;
        Asite.assign(4 * NA, 0); Bsite.assign(4 * NB, 0); Aof.assign(N, 0); Bof.assign(N, 0);
        for (int Y = 0; Y < Lb; Y++)
            for (int X = 0; X < Lb; X++) {
                int c = Y * Lb + X;
                for (int p = 0; p < 4; p++) {
                    int dx = p & 1, dy = p >> 1;
                    int ia = idx(2 * X + 2 + dx, 2 * Y + 2 + dy);
                    int ib = idx(2 * X + 1 + dx, 2 * Y + 1 + dy);
                    Asite[4 * c + p] = ia; Aof[ia] = c;
                    Bsite[4 * c + p] = ib; Bof[ib] = c;
                }
            }
    }
    inline int idx(int x, int y) const { return ((y % L + L) % L) * L + ((x % L + L) % L); }

    // ---------------- plaquette energies ----------------
    inline double EA(const int8_t* s, int a) const {
        const int* q = &Asite[4 * a];
        int s0 = s[q[0]], s1 = s[q[1]], s2 = s[q[2]], s3 = s[q[3]];
        return J1 * (s0 * s1 + s2 * s3 + s0 * s2 + s1 * s3) + J2 * (s0 * s3 + s1 * s2);
    }
    inline double EB(const int8_t* s, int b) const {
        const int* q = &Bsite[4 * b];
        int S = s[q[0]] + s[q[1]] + s[q[2]] + s[q[3]];
        return 0.5 * JB * S * S;
    }
    // energy of an isolated A plaquette with spins s0..s3 (positions 0..3); used by unit tests
    double EA4(int s0, int s1, int s2, int s3) const {
        return J1 * (s0 * s1 + s2 * s3 + s0 * s2 + s1 * s3) + J2 * (s0 * s3 + s1 * s2);
    }
    double total_energy(const int8_t* s) const {
        double e = 0; long m = 0;
        for (int a = 0; a < NA; a++) e += EA(s, a);
        for (int b = 0; b < NB; b++) e += EB(s, b);
        for (int i = 0; i < N; i++) m += s[i];
        return e - h * m;
    }

    // ---------------- replica handling ----------------
    inline void flip(Replica& r, int i) const {
        r.s[i] = -r.s[i];
        int b = Bof[i];
        int X = b % Lb, Y = b / Lb;
        r.colB[X] += 2 * r.s[i];
        r.rowB[Y] += 2 * r.s[i];
    }
    void rebuild(Replica& r) const {
        r.colB.assign(Lb, 0); r.rowB.assign(Lb, 0);
        for (int b = 0; b < NB; b++) {
            int X = b % Lb, Y = b / Lb, S = 0;
            for (int p = 0; p < 4; p++) S += r.s[Bsite[4 * b + p]];
            r.colB[X] += S; r.rowB[Y] += S;
        }
        r.E = total_energy(r.s.data());
    }
    // kind: 0 neel, 1 stripe in X (m_x), 2 stripe in Y (m_y), 3 random, 4 ferro
    void init(Replica& r, int kind, uint64_t seed, uint64_t stream) const {
        r.rng = Rng(seed, stream);
        r.s.assign(N, 1);
        for (int i = 0; i < N; i++) {
            int b = Bof[i], X = b % Lb, Y = b / Lb;
            int sg = 1;
            if (kind == 0) sg = ((X + Y) & 1) ? -1 : 1;
            else if (kind == 1) sg = (X & 1) ? -1 : 1;
            else if (kind == 2) sg = (Y & 1) ? -1 : 1;
            else if (kind == 3) sg = (r.rng.next() >> 63) ? 1 : -1;
            r.s[i] = (int8_t)sg;
        }
        rebuild(r);
    }

    // ---------------- energy changes (non-mutating) ----------------
    inline double dE_single(const Replica& r, int i) const {
        const int8_t* s = r.s.data();
        int a = Aof[i], b = Bof[i];
        const int* qa = &Asite[4 * a];
        const int* qb = &Bsite[4 * b];
        int p = 0;
        for (; p < 4; p++) if (qa[p] == i) break;
        int sg = s[i];
        double eA = sg * (J1 * (s[qa[p ^ 1]] + s[qa[p ^ 2]]) + J2 * s[qa[p ^ 3]]);
        int S = s[qb[0]] + s[qb[1]] + s[qb[2]] + s[qb[3]];
        return -2 * eA + 2 * JB * (1 - sg * S) + 2 * h * sg;
    }
    // flip an arbitrary small set of spins (all in distinct A plaquettes or not): returns dE
    double dE_flipset(Replica& r, const int* sites, int n) const {
        int as[4], bs[4], na = 0, nb = 0;
        for (int k = 0; k < n; k++) {
            int a = Aof[sites[k]], b = Bof[sites[k]];
            bool fa = false, fb = false;
            for (int j = 0; j < na; j++) fa |= (as[j] == a);
            for (int j = 0; j < nb; j++) fb |= (bs[j] == b);
            if (!fa) as[na++] = a;
            if (!fb) bs[nb++] = b;
        }
        int8_t* s = r.s.data();
        double e0 = 0, f0 = 0;
        for (int j = 0; j < na; j++) e0 += EA(s, as[j]);
        for (int j = 0; j < nb; j++) e0 += EB(s, bs[j]);
        for (int k = 0; k < n; k++) f0 += s[sites[k]];
        for (int k = 0; k < n; k++) s[sites[k]] = -s[sites[k]];
        double e1 = 0;
        for (int j = 0; j < na; j++) e1 += EA(s, as[j]);
        for (int j = 0; j < nb; j++) e1 += EB(s, bs[j]);
        for (int k = 0; k < n; k++) s[sites[k]] = -s[sites[k]];
        return e1 - e0 + 2 * h * f0;
    }
    inline void apply_flipset(Replica& r, const int* sites, int n, double dE) const {
        for (int k = 0; k < n; k++) flip(r, sites[k]);
        r.E += dE;
    }
    // energy change of flipping the superspins of B columns (orient 0) / rows (orient 1)
    // x0 .. x0+w-1 (mod Lb), 1 <= w <= Lb-1.  Only the boundary A plaquettes change (plus field).
    double dE_strip(const Replica& r, int orient, int x0, int w) const {
        const int8_t* s = r.s.data();
        double d = 0;
        int cA, cB, mA, mB;
        if (orient == 0) { cA = (x0 - 1 + Lb) % Lb; cB = (x0 + w - 1) % Lb; mA = 0b1010; mB = 0b0101; }
        else             { cA = (x0 - 1 + Lb) % Lb; cB = (x0 + w - 1) % Lb; mA = 0b1100; mB = 0b0011; }
        for (int t = 0; t < Lb; t++) {
            int a1 = orient == 0 ? t * Lb + cA : cA * Lb + t;
            int a2 = orient == 0 ? t * Lb + cB : cB * Lb + t;
            d += EA_masked(s, a1, mA) - EA(s, a1);
            d += EA_masked(s, a2, mB) - EA(s, a2);
        }
        long sum = 0;
        const std::vector<int>& line = orient == 0 ? r.colB : r.rowB;
        for (int k = 0; k < w; k++) sum += line[(x0 + k) % Lb];
        return d + 2 * h * sum;
    }
    inline double EA_masked(const int8_t* s, int a, int mask) const {
        const int* q = &Asite[4 * a];
        int v[4];
        for (int p = 0; p < 4; p++) v[p] = ((mask >> p) & 1) ? -s[q[p]] : s[q[p]];
        return J1 * (v[0] * v[1] + v[2] * v[3] + v[0] * v[2] + v[1] * v[3]) + J2 * (v[0] * v[3] + v[1] * v[2]);
    }
    void apply_strip(Replica& r, int orient, int x0, int w, double dE) const {
        for (int k = 0; k < w; k++) {
            int c = (x0 + k) % Lb;
            for (int t = 0; t < Lb; t++) {
                int b = orient == 0 ? t * Lb + c : c * Lb + t;
                for (int p = 0; p < 4; p++) flip(r, Bsite[4 * b + p]);
            }
        }
        r.E += dE;
    }

    // ---------------- Monte Carlo sweep ----------------
    static inline bool accept(double dE, double beta, Rng& g) {
        return dE <= 0 || g.uniform() < std::exp(-beta * dE);
    }
    static constexpr int HALF_P[4][2] = {{0, 1}, {2, 3}, {0, 2}, {1, 3}};

    // one sweep = N single-spin, NB half flips, NB superspin flips, n_strip strip flips.
    // tried/acc (length N_MOVES) accumulate the statistics.
    void sweep(Replica& r, double beta, int n_strip, long* tried, long* acc, unsigned moves = 15) const {
        Rng& g = r.rng;
        if (moves & 1)
            for (int k = 0; k < N; k++) {
                int i = g.below(N);
                double dE = dE_single(r, i);
                tried[MV_SINGLE]++;
                if (accept(dE, beta, g)) { flip(r, i); r.E += dE; acc[MV_SINGLE]++; }
            }
        if (moves & 2)
            for (int k = 0; k < NB; k++) {
                int b = g.below(NB), e = g.below(4);
                int st[2] = {Bsite[4 * b + HALF_P[e][0]], Bsite[4 * b + HALF_P[e][1]]};
                double dE = dE_flipset(r, st, 2);
                tried[MV_HALF]++;
                if (accept(dE, beta, g)) { apply_flipset(r, st, 2, dE); acc[MV_HALF]++; }
            }
        if (moves & 4)
            for (int k = 0; k < NB; k++) {
                int b = g.below(NB);
                const int* st = &Bsite[4 * b];
                double dE = dE_flipset(r, st, 4);
                tried[MV_SUPER]++;
                if (accept(dE, beta, g)) { apply_flipset(r, st, 4, dE); acc[MV_SUPER]++; }
            }
        if (moves & 8)
            for (int k = 0; k < n_strip; k++) {
                int orient = g.below(2), x0 = g.below(Lb), w = 1 + g.below(Lb - 1);
                double dE = dE_strip(r, orient, x0, w);
                tried[MV_STRIP]++;
                if (accept(dE, beta, g)) { apply_strip(r, orient, x0, w, dE); acc[MV_STRIP]++; }
            }
    }

    // ---------------- observables ----------------
    Obs measure(const Replica& r) const {
        const int8_t* s = r.s.data();
        double mN = 0, mx = 0, my = 0, E = 0;
        long nBal = 0, nd = 0, nn = 0, nv = 0, mag = 0;
        for (int b = 0; b < NB; b++) {
            int X = b % Lb, Y = b / Lb;
            const int* q = &Bsite[4 * b];
            int S = s[q[0]] + s[q[1]] + s[q[2]] + s[q[3]];
            double tau = S * 0.25;
            mN += ((X + Y) & 1) ? -tau : tau;
            mx += (X & 1) ? -tau : tau;
            my += (Y & 1) ? -tau : tau;
            if (S == 4 || S == -4) nBal++;
            E += 0.5 * JB * S * S;
            mag += S;
        }
        for (int a = 0; a < NA; a++) {
            const int* q = &Asite[4 * a];
            int s0 = s[q[0]], s1 = s[q[1]], s2 = s[q[2]], s3 = s[q[3]];
            int S = s0 + s1 + s2 + s3;
            if (S == 0) { if (s0 == s3) nd++; else nn++; } else nv++;
            E += J1 * (s0 * s1 + s2 * s3 + s0 * s2 + s1 * s3) + J2 * (s0 * s3 + s1 * s2);
        }
        Obs o;
        o.E = (E - h * mag) / N;
        o.mN = mN / NB; o.mx = mx / NB; o.my = my / NB;
        o.rdiag = (double)nd / NA; o.rnbr = (double)nn / NA; o.rviol = (double)nv / NA;
        o.rBal = (double)nBal / NB;
        return o;
    }
};
