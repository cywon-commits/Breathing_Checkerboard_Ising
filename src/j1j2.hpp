// Plain J1-J2 Ising model on an Lb x Lb torus (energy, single-spin and strip-flip energy changes).
#pragma once
#include <cstdint>
#include <vector>

struct J1J2 {
    int Lb, N;
    double J1, J2, hf;
    J1J2(int Lb_, double J1_, double J2_, double hf_) : Lb(Lb_), N(Lb_ * Lb_), J1(J1_), J2(J2_), hf(hf_) {}
    inline int id(int x, int y) const { return ((y + Lb) % Lb) * Lb + ((x + Lb) % Lb); }
    double energy(const std::vector<int8_t>& t) const {
        double e = 0; long m = 0;
        for (int y = 0; y < Lb; y++)
            for (int x = 0; x < Lb; x++) {
                int s = t[id(x, y)];
                e += J1 * s * (t[id(x + 1, y)] + t[id(x, y + 1)]) + J2 * s * (t[id(x + 1, y + 1)] + t[id(x + 1, y - 1)]);
                m += s;
            }
        return e - hf * m;
    }
    inline double dE_single(const std::vector<int8_t>& t, int x, int y) const {
        int s = t[id(x, y)];
        int nn = t[id(x + 1, y)] + t[id(x - 1, y)] + t[id(x, y + 1)] + t[id(x, y - 1)];
        int nnn = t[id(x + 1, y + 1)] + t[id(x - 1, y - 1)] + t[id(x + 1, y - 1)] + t[id(x - 1, y + 1)];
        return -2.0 * s * (J1 * nn + J2 * nnn) + 2.0 * hf * s;
    }
    // flip columns x0..x0+w-1 (orient 0) or rows (orient 1): dE = -2 sum over bonds with exactly one end inside
    // + field part.  Evaluated from the two boundary lines only.
    double dE_strip(const std::vector<int8_t>& t, int orient, int x0, int w) const {
        auto inside = [&](int c) { return ((c - x0) % Lb + Lb) % Lb < w; };
        double d = 0; long sum = 0;
        int lines[2] = {x0 % Lb, (x0 + w - 1) % Lb};
        int nl = (w == 1) ? 1 : 2;
        for (int k = 0; k < nl; k++) {
            int c = lines[k];
            for (int u = 0; u < Lb; u++) {
                int x = orient == 0 ? c : u, y = orient == 0 ? u : c;
                int s = t[id(x, y)];
                for (int dx = -1; dx <= 1; dx++)
                    for (int dy = -1; dy <= 1; dy++) {
                        if (!dx && !dy) continue;
                        int xx = x + dx, yy = y + dy;
                        int cc = orient == 0 ? ((xx % Lb) + Lb) % Lb : ((yy % Lb) + Lb) % Lb;
                        if (inside(cc)) continue;
                        double J = (dx && dy) ? J2 : J1;
                        d += -2.0 * J * s * t[id(xx, yy)];
                    }
            }
        }
        for (int k = 0; k < w; k++) {
            int c = (x0 + k) % Lb;
            for (int u = 0; u < Lb; u++) sum += orient == 0 ? t[id(c, u)] : t[id(u, c)];
        }
        return d + 2.0 * hf * sum;
    }
    void apply_strip(std::vector<int8_t>& t, int orient, int x0, int w) const {
        for (int k = 0; k < w; k++) {
            int c = (x0 + k) % Lb;
            for (int u = 0; u < Lb; u++) { int8_t& s = orient == 0 ? t[id(c, u)] : t[id(u, c)]; s = -s; }
        }
    }
};

