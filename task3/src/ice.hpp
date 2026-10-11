// Staggered six-vertex (ice) model on an L x L A-lattice torus (site lattice 2L x 2L): directed-loop (worm) Monte Carlo.
//
// Geometry as in the task-1 code: A plaquette (X,Y) has origin site (2X+2, 2Y+2), B plaquette (X,Y) has origin (2X+1, 2Y+1);
// position p = dx + 2 dy inside a plaquette (0 = SW, 1 = SE, 2 = NW, 3 = NE).  Every site (spin) belongs to exactly one A and one B plaquette.
// Ice rule: every plaquette has two up and two down spins.  Vertex model: vertices = plaquettes, edges = spins; an edge points
// from A to B if sigma = +1.  'out' positions of an A vertex = up spins, of a B vertex = down spins (the same arrow is out of A and into B).
// Weights: A vertex of c type (the two up spins on a diagonal) e^u, B vertex all 1 (staggered model) or the same c (uniform F model).
//
// Directed-loop update (see REPORT_task3.md): pick a random edge e0, flip it, the head enters its A vertex; at every vertex the exit edge
// is drawn from the table P(entry -> exit | pre-state) = m(S,S')/W(S) (symmetric m, bounce allowed), flipped, and the head moves on; the
// loop closes when the head reaches the tail vertex (the B vertex of e0, weight ratio 1).  Detailed balance follows from m(S,S') = m(S',S).
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
        uint64_t x = seed * 0x2545F4914F6CDD1DULL + stream * 0x9E3779B97F4A7C15ULL + 777;
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

// ---------------------------------------------------------------- scattering tables
// state mask S: bit p set = position p is an 'out' edge of the vertex (2 bits set for an ice state)
static inline bool is_diag(int S) { return S == 0b1001 || S == 0b0110; }
static inline double vertex_weight(int S, double c) { return is_diag(S) ? c : 1.0; }

// symmetric directed-loop weights m(S,S') for neighbouring ice states S' = S - {p} + {q}
static inline double m_weight(int S, int S2, double c) {
    bool dS = is_diag(S), dS2 = is_diag(S2);
    if (dS && dS2) return 0.0;                       // never neighbours
    double x = c >= 2.0 ? 1.0 : 0.5 * c, y = c >= 2.0 ? 0.0 : 1.0 - 0.5 * c;
    return (dS || dS2) ? x : y;
}

struct Table {
    // cum[S][pe][q] : cumulative probability of exit position q (q == pe is the bounce) for entry position pe at a vertex in pre-state S
    double P[16][4][4];
    double cum[16][4][4];
    void build(double c) {
        std::memset(P, 0, sizeof P); std::memset(cum, 0, sizeof cum);
        for (int S = 0; S < 16; S++) {
            if (__builtin_popcount(S) != 2) continue;
            double W = vertex_weight(S, c);
            for (int pe = 0; pe < 4; pe++) {
                bool pe_out = (S >> pe) & 1;
                double tot = 0;
                for (int q = 0; q < 4; q++) {
                    if (q == pe) continue;
                    bool q_out = (S >> q) & 1;
                    if (q_out == pe_out) continue;            // exit must have the opposite polarity
                    int S2 = S ^ (1 << pe) ^ (1 << q);
                    P[S][pe][q] = m_weight(S, S2, c) / W;
                    tot += P[S][pe][q];
                }
                P[S][pe][pe] = 1.0 - tot;                      // bounce
                double a = 0;
                for (int q = 0; q < 4; q++) { a += P[S][pe][q]; cum[S][pe][q] = a; }
            }
        }
    }
};

struct Ice {
    int L, Ls, N, NV;                                  // A lattice size, site lattice size 2L, #spins, #A (= #B) vertices (L*L)
    std::vector<int> Asite, Bsite, Aof, Bof;
    std::vector<int8_t> Apos, Bpos;
    Ice(int L_) : L(L_), Ls(2 * L_), N(4 * L_ * L_), NV(L_ * L_) {
        Asite.assign(4 * NV, 0); Bsite.assign(4 * NV, 0); Aof.assign(N, 0); Bof.assign(N, 0); Apos.assign(N, 0); Bpos.assign(N, 0);
        for (int Y = 0; Y < L; Y++)
            for (int X = 0; X < L; X++) {
                int c = Y * L + X;
                for (int p = 0; p < 4; p++) {
                    int dx = p & 1, dy = p >> 1;
                    int ia = idx(2 * X + 2 + dx, 2 * Y + 2 + dy), ib = idx(2 * X + 1 + dx, 2 * Y + 1 + dy);
                    Asite[4 * c + p] = ia; Aof[ia] = c; Apos[ia] = (int8_t)p;
                    Bsite[4 * c + p] = ib; Bof[ib] = c; Bpos[ib] = (int8_t)p;
                }
            }
    }
    inline int idx(int x, int y) const { return ((y % Ls + Ls) % Ls) * Ls + ((x % Ls + Ls) % Ls); }

    inline int stateA(const int8_t* s, int a) const {
        const int* q = &Asite[4 * a];
        return (s[q[0]] > 0) | ((s[q[1]] > 0) << 1) | ((s[q[2]] > 0) << 2) | ((s[q[3]] > 0) << 3);
    }
    inline int stateB(const int8_t* s, int b) const {
        const int* q = &Bsite[4 * b];
        return (s[q[0]] < 0) | ((s[q[1]] < 0) << 1) | ((s[q[2]] < 0) << 2) | ((s[q[3]] < 0) << 3);
    }
    bool is_ice(const int8_t* s) const {
        for (int v = 0; v < NV; v++)
            if (__builtin_popcount(stateA(s, v)) != 2 || __builtin_popcount(stateB(s, v)) != 2) return false;
        return true;
    }
    int count_c(const int8_t* s, bool both) const {
        int n = 0;
        for (int v = 0; v < NV; v++) {
            n += is_diag(stateA(s, v));
            if (both) n += is_diag(stateB(s, v));
        }
        return n;
    }

    // ---------------- initial ice states: all A vertices c type, sigma(SW)=s, SE=NW=-s, NE=s in A block (X,Y)
    // kind: 0 'row alternating' s=(-1)^Y (Q=+1/2), 1 'column alternating' s=(-1)^X (Q=-1/2), 2 random row phases s=b(Y), 3 random column phases s=a(X)
    void init_state(std::vector<int8_t>& s, int kind, Rng& g) const {
        s.assign(N, 1);
        std::vector<int> ph(L);
        for (int i = 0; i < L; i++) ph[i] = (kind == 0 || kind == 1) ? ((i & 1) ? -1 : 1) : ((g.next() >> 63) ? 1 : -1);
        for (int Y = 0; Y < L; Y++)
            for (int X = 0; X < L; X++) {
                int sg = (kind == 0 || kind == 2) ? ph[Y] : ph[X];
                const int* q = &Asite[4 * (Y * L + X)];
                s[q[0]] = sg; s[q[1]] = -sg; s[q[2]] = -sg; s[q[3]] = sg;
            }
    }

    // ---------------- one directed-loop update; returns the number of edge traversals (0 if aborted and undone)
    // cB != 1 (uniform F model): the tail vertex is a weighted B vertex, so the closed loop is accepted with min(1, W(final)/W(initial)) there.
    int loop_update(int8_t* s, Rng& g, const Table& tA, const Table& tB, std::vector<int>& flips, int maxsteps, double cB = 1.0) const {
        int e0 = g.below(N);
        int tail = Bof[e0];
        int Stail0 = cB != 1.0 ? stateB(s, tail) : 0;
        flips.clear();
        int type = 0, v = Aof[e0], pe = Apos[e0];
        int Spre = stateA(s, v);
        s[e0] = -s[e0]; flips.push_back(e0);
        int steps = 1;
        for (;;) {
            const Table& t = type == 0 ? tA : tB;
            double u = g.uniform();
            const double* cm = t.cum[Spre][pe];
            int q = u < cm[0] ? 0 : u < cm[1] ? 1 : u < cm[2] ? 2 : 3;
            int e = type == 0 ? Asite[4 * v + q] : Bsite[4 * v + q];
            // move across edge e to the other vertex
            int ntype = 1 - type, nv, npe, nS;
            if (type == 0) { nv = Bof[e]; npe = Bpos[e]; nS = ntype == 1 ? stateB(s, nv) : 0; }
            else { nv = Aof[e]; npe = Apos[e]; nS = stateA(s, nv); }
            s[e] = -s[e]; flips.push_back(e); steps++;
            if (ntype == 1 && nv == tail) {                      // head reaches the tail vertex: loop closed
                if (cB != 1.0) {
                    double r = vertex_weight(stateB(s, tail), cB) / vertex_weight(Stail0, cB);
                    if (r < 1.0 && g.uniform() >= r) { for (int f : flips) s[f] = -s[f]; return 0; }
                }
                return steps;
            }
            type = ntype; v = nv; pe = npe; Spre = nS;
            if (steps > maxsteps) {                              // symmetric abort (depends only on the path length): undo
                for (int f : flips) s[f] = -s[f];
                return 0;
            }
        }
    }

    // ---------------- observables
    struct Obs {
        int nc, ncB, Nh, Nv, Wx, Wy;
        double Q;
        double Sx0, Sx1, Sy0, Sy1;          // structure factors of the diagonal field d along x / along y at k = 0, 2 pi / L
    };
    Obs measure(const int8_t* s) const {
        Obs o;
        std::vector<int8_t> d(NV, 0);
        o.nc = 0; o.ncB = 0;
        for (int v = 0; v < NV; v++) {
            int S = stateA(s, v);
            if (S == 0b1001) { d[v] = 1; o.nc++; } else if (S == 0b0110) { d[v] = -1; o.nc++; }
            if (is_diag(stateB(s, v))) o.ncB++;
        }
        o.Nh = o.Nv = 0;
        for (int Y = 0; Y < L; Y++)
            for (int X = 0; X < L; X++) {
                int a = d[Y * L + X];
                if (!a) continue;
                int r = d[Y * L + (X + 1) % L], u = d[((Y + 1) % L) * L + X];
                if (r == a) o.Nh++;
                if (u == a) o.Nv++;
            }
        o.Q = (double)(o.Nh - o.Nv) / (2.0 * NV);
        // structure factors
        double c1 = std::cos(2 * M_PI / L), s1 = std::sin(2 * M_PI / L);
        double sx0 = 0, sx1 = 0, sy0 = 0, sy1 = 0;
        for (int Y = 0; Y < L; Y++) {
            double sum = 0, re = 0, im = 0, cr = 1, ci = 0;
            for (int X = 0; X < L; X++) {
                double dd = d[Y * L + X]; sum += dd; re += dd * cr; im += dd * ci;
                double t = cr * c1 - ci * s1; ci = cr * s1 + ci * c1; cr = t;
            }
            sx0 += sum * sum / L; sx1 += (re * re + im * im) / L;
        }
        for (int X = 0; X < L; X++) {
            double sum = 0, re = 0, im = 0, cr = 1, ci = 0;
            for (int Y = 0; Y < L; Y++) {
                double dd = d[Y * L + X]; sum += dd; re += dd * cr; im += dd * ci;
                double t = cr * c1 - ci * s1; ci = cr * s1 + ci * c1; cr = t;
            }
            sy0 += sum * sum / L; sy1 += (re * re + im * im) / L;
        }
        o.Sx0 = sx0 / L; o.Sx1 = sx1 / L; o.Sy0 = sy0 / L; o.Sy1 = sy1 / L;
        // windings = conserved fluxes: row flux sum_x sigma(x,y) (-1)^y / 2 and column flux sum_y sigma(x,y) (-1)^x / 2
        int fr = 0, fc = 0;
        for (int x = 0; x < Ls; x++) fr += s[idx(x, 0)];
        for (int y = 0; y < Ls; y++) fc += s[idx(0, y)];
        o.Wx = fr / 2; o.Wy = fc / 2;
        return o;
    }
    // consistency check of the conserved fluxes (used in tests)
    bool fluxes_consistent(const int8_t* s) const {
        int fr0 = 0, fc0 = 0;
        for (int x = 0; x < Ls; x++) fr0 += s[idx(x, 0)];
        for (int y = 0; y < Ls; y++) fc0 += s[idx(0, y)];
        for (int y = 0; y < Ls; y++) { int f = 0; for (int x = 0; x < Ls; x++) f += s[idx(x, y)]; if (f * ((y & 1) ? -1 : 1) != fr0) return false; }
        for (int x = 0; x < Ls; x++) { int f = 0; for (int y = 0; y < Ls; y++) f += s[idx(x, y)]; if (f * ((x & 1) ? -1 : 1) != fc0) return false; }
        return true;
    }
};
