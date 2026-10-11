// Sector-compressed transfer matrix (fixed number n of up spins, arrow representation: every gate conserves n).
// State s = (hi, lo) with lo = low H1 bits, hi = high W-H1 bits; index = base[a] + rankHi[hi]*cntLo[n-a] + rankLo[lo],
// a = popcount(hi).  Used for W = 28 (M=0 sector: 40 million states).
#include <cstdint>
#include <vector>
#include <omp.h>
#include <algorithm>

static int W_, n_, H1_, H2_;              // total width, up spins, low bits, high bits
static std::vector<int32_t> rankLo, rankHi;   // rank of a word among words of equal popcount
static std::vector<std::vector<uint32_t>> wordsLo, wordsHi;   // words by popcount
static std::vector<int64_t> base_;        // block offsets (by popcount a of the high part)
static std::vector<int64_t> cntLoPop, cntHiPop;
static int64_t dim_ = 0;

static inline int64_t idx_of(uint64_t s) {
    uint32_t lo = (uint32_t)(s & ((1ULL << H1_) - 1)), hi = (uint32_t)(s >> H1_);
    int a = __builtin_popcount(hi);
    return base_[a] + (int64_t)rankHi[hi] * cntLoPop[n_ - a] + rankLo[lo];
}

extern "C" {
int64_t sector_init(int W, int n) {
    W_ = W; n_ = n; H1_ = W / 2; H2_ = W - H1_;
    auto build = [](int bits, std::vector<int32_t>& rank, std::vector<std::vector<uint32_t>>& words) {
        rank.assign(1u << bits, 0); words.assign(bits + 1, {});
        for (uint32_t w = 0; w < (1u << bits); w++) { int p = __builtin_popcount(w); rank[w] = (int32_t)words[p].size(); words[p].push_back(w); }
    };
    build(H1_, rankLo, wordsLo); build(H2_, rankHi, wordsHi);
    cntLoPop.assign(H1_ + 1, 0); cntHiPop.assign(H2_ + 1, 0);
    for (int p = 0; p <= H1_; p++) cntLoPop[p] = wordsLo[p].size();
    for (int p = 0; p <= H2_; p++) cntHiPop[p] = wordsHi[p].size();
    base_.assign(H2_ + 2, 0);
    int64_t acc = 0;
    for (int a = 0; a <= H2_; a++) {
        base_[a] = acc;
        if (n - a >= 0 && n - a <= H1_) acc += cntHiPop[a] * cntLoPop[n - a];
    }
    base_[H2_ + 1] = acc; dim_ = acc;
    return dim_;
}

// enumerate the states: states[k] = s of index k (needed once by the Python side for the start vector, optional)
void sector_states(uint64_t* out) {
    #pragma omp parallel for schedule(dynamic)
    for (int a = 0; a <= H2_; a++) {
        if (n_ - a < 0 || n_ - a > H1_) continue;
        int64_t k = base_[a];
        for (uint32_t h : wordsHi[a]) for (uint32_t l : wordsLo[n_ - a]) out[k++] = ((uint64_t)h << H1_) | l;
    }
}

static void gate_sector(double* v, int i, int j, double c) {
    const uint64_t mi = 1ULL << i, mj = 1ULL << j, flip = mi | mj;
    #pragma omp parallel for schedule(dynamic)
    for (int a = 0; a <= H2_; a++) {
        if (n_ - a < 0 || n_ - a > H1_) continue;
        int64_t k = base_[a];
        for (uint32_t h : wordsHi[a]) {
            for (uint32_t l : wordsLo[n_ - a]) {
                uint64_t s = ((uint64_t)h << H1_) | l;
                if ((s & mi) && !(s & mj)) {           // p = 1 state; partner p = 2
                    int64_t k2 = idx_of(s ^ flip);
                    double vb = v[k], vc = v[k2];
                    v[k] = c * vb + vc; v[k2] = vb + c * vc;
                }
                k++;
            }
        }
    }
}

void sector_apply_T(const double* in, double* out, double cA, double cB) {
    if (out != in) {
        #pragma omp parallel for schedule(static)
        for (int64_t k = 0; k < dim_; k++) out[k] = in[k];
    }
    for (int k = 0; k < W_ / 2; k++) gate_sector(out, 2 * k, 2 * k + 1, cA);
    for (int k = 0; k < W_ / 2; k++) gate_sector(out, 2 * k + 1, (2 * k + 2) % W_, cB);
}
}
