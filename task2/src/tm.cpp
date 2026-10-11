// Sector transfer matrix of the staggered six-vertex (ice) model on a strip of W spin columns (task 2).
//
// Row state = W spins (bit i = column i, 1 = up).  One transfer step goes up two site rows: an A-plaquette row
// (plaquettes on column pairs (2k,2k+1), weight c_A for the "diagonal" vertex, 1 otherwise) followed by a
// B-plaquette row (column pairs (2k+1,2k+2 mod W), weight c_B).  Every plaquette has two up and two down spins
// (ice rule), built into the local gate.
//
// Arrow (staggered) representation: the spins of each row are multiplied by (-1)^row, so that the ice rule
// becomes "number of up spins in the bottom pair = number of up spins in the top pair".  The local gate then
// conserves the up-spin number of the pair -> the full step conserves M = sum sigma exactly, and no weight can leak
// between sectors (the 'false degeneracy' warned about in the task):
//     pair state p = bit_l + 2 bit_r :  (0,0)->(0,0): 1 ; (1,1)->(1,1): 1 ; (1,0)->(0,1): 1 ; (0,1)->(1,0): 1 ;
//                                       (1,0)->(1,0): c ; (0,1)->(0,1): c      (c = diagonal-type plaquette)
// (in the original spin picture "c" is the plaquette whose two up spins sit on a diagonal).
#include <cstdint>
#include <cstring>
#include <omp.h>

static inline void gate(double* v, int W, int i, int j, double c) {
    const uint64_t mi = 1ULL << i, mj = 1ULL << j;
    int lo = i < j ? i : j, hi = i < j ? j : i;
    const uint64_t N4 = 1ULL << (W - 2);
    const uint64_t lowmask = (1ULL << lo) - 1, midmask = (1ULL << (hi - lo - 1)) - 1;
#pragma omp parallel for schedule(static)
    for (int64_t x = 0; x < (int64_t)N4; x++) {
        uint64_t idx = (x & lowmask) | (((x >> lo) & midmask) << (lo + 1)) | ((x >> (hi - 1)) << (hi + 1));
        double vb = v[idx | mi], vc = v[idx | mj];   // p=1 (bit i set), p=2 (bit j set)
        v[idx | mi] = c * vb + vc;
        v[idx | mj] = vb + c * vc;
    }
}

extern "C" {
// w = T v, T = T_B T_A ; work buffer not needed (in-place gates on a copy)
void apply_T(const double* in, double* out, int W, double cA, double cB) {
    const uint64_t N = 1ULL << W;
    if (out != in) {
#pragma omp parallel for schedule(static)
        for (int64_t k = 0; k < (int64_t)N; k++) out[k] = in[k];
    }
    for (int k = 0; k < W / 2; k++) gate(out, W, 2 * k, 2 * k + 1, cA);
    for (int k = 0; k < W / 2; k++) gate(out, W, 2 * k + 1, (2 * k + 2) % W, cB);
}
// same for the half steps (used by tests)
void apply_A(double* v, int W, double cA) { for (int k = 0; k < W / 2; k++) gate(v, W, 2 * k, 2 * k + 1, cA); }
void apply_B(double* v, int W, double cB) { for (int k = 0; k < W / 2; k++) gate(v, W, 2 * k + 1, (2 * k + 2) % W, cB); }
void set_threads(int n) { omp_set_num_threads(n); }
}
