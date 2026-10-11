// Parallel-tempering (in u) directed-loop Monte Carlo of the staggered six-vertex model (task 3).
//   ice_pt --L 32 --ulist 0.8,1.0,... --neq 20000 --nmeas 100000 --nblocks 20 --seed 1 --init mixed --out prefix
// outputs: <out>.blocks.csv  block averages per u;  <out>.series.f32 (Q, nc/L^2, Wx, Wy per replica, float32);  <out>.meta.json
#include <chrono>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>

#include "ice.hpp"

int main(int argc, char** argv) {
    std::map<std::string, std::string> a;
    for (int i = 1; i + 1 < argc; i += 2) a[argv[i] + 2] = argv[i + 1];
    auto get = [&](const char* k, const char* d) { return a.count(k) ? a[k] : std::string(d); };
    int L = std::stoi(get("L", "16"));
    std::string model = get("model", "stag"), initname = get("init", "mixed"), out = get("out", "icerun"), commit = get("commit", "unknown");
    long neq = std::stol(get("neq", "20000")), nmeas = std::stol(get("nmeas", "100000"));
    int nblocks = std::stoi(get("nblocks", "20")), meas_every = std::stoi(get("meas_every", "1")), series_every = std::stoi(get("series_every", "10"));
    uint64_t seed = std::stoull(get("seed", "1"));
    std::vector<double> U;
    {
        std::string ul = get("ulist", "1.0");
        size_t p = 0;
        while (p < ul.size()) { size_t q = ul.find(',', p); if (q == std::string::npos) q = ul.size(); U.push_back(std::stod(ul.substr(p, q - p))); p = q + 1; }
    }
    const int nrep = (int)U.size();
    const bool isF = model == "F";
    Ice M(L);

    std::vector<std::vector<int8_t>> spins(nrep);
    std::vector<Rng> rngs;
    for (int c = 0; c < nrep; c++) rngs.emplace_back(seed, 100 + c);
    for (int c = 0; c < nrep; c++) {
        int kind = initname == "row" ? 0 : initname == "col" ? 1 : initname == "rowrand" ? 2 : initname == "colrand" ? 3 : (c % 4);
        M.init_state(spins[c], kind, rngs[c]);
    }
    std::vector<Table> tA(nrep), tB(nrep);
    for (int t = 0; t < nrep; t++) { double cA = std::exp(U[t]); tA[t].build(cA); tB[t].build(isF ? cA : 1.0); }
    std::vector<int> perm(nrep);
    for (int t = 0; t < nrep; t++) perm[t] = t;
    Rng swapRng(seed, 424242);
    std::vector<int> nloops(nrep, std::max(1, M.N / 4));
    std::vector<double> lensum(nrep, 0.0); std::vector<long> loopcnt(nrep, 0), abortcnt(nrep, 0);
    std::vector<long> swt(nrep, 0), swa(nrep, 0);
    std::vector<int> sgn(nrep, 0), switches(nrep, 0);     // per configuration: last |Q|>0.1 sign and number of sign changes

    const int NF = 21;
    const char* names[NF] = {"nc", "ncB", "Q", "Q2", "Q4", "Qabs", "Wx2", "Wy2", "Wx4", "Wy4", "Wpar2", "Wper2", "Spar0", "Spar1", "Sper0", "Sper1", "Sx0", "Sx1", "Sy0", "Sy1", "Wpar2sym"};
    long nrec = nmeas / meas_every, blen = std::max(1L, nrec / nblocks);
    int nbl = (int)(nrec / blen);
    std::vector<double> acc(nrep * NF, 0.0), blocks;
    std::vector<float> series;
    long rec = 0, inblock = 0, warm = std::min(neq / 2, 1000L);
    auto t0 = std::chrono::steady_clock::now();
    const int maxsteps = 200 * M.N;

    for (long it = 0; it < neq + nmeas; it++) {
        bool measuring = it >= neq && ((it - neq) % meas_every == 0) && rec < (long)nbl * blen;
        std::vector<Ice::Obs> obs(nrep);
        std::vector<int> ncnow(nrep);
#pragma omp parallel for schedule(static)
        for (int t = 0; t < nrep; t++) {
            int c = perm[t];
            std::vector<int> flips; flips.reserve(1024);
            long steps = 0, lps = 0, ab = 0;
            for (int l = 0; l < nloops[t]; l++) {
                int st = M.loop_update(spins[c].data(), rngs[c], tA[t], tB[t], flips, maxsteps, isF ? std::exp(U[t]) : 1.0);
                if (st == 0) ab++; else steps += st;
                lps++;
            }
            lensum[t] += steps; loopcnt[t] += lps; abortcnt[t] += ab;
            ncnow[t] = M.count_c(spins[c].data(), isF);
            if (measuring) obs[t] = M.measure(spins[c].data());
        }
        if (it == warm && neq > 0) {   // fix the number of loops per sweep (mean traversed edges ~ N) after a short warm-up
            for (int t = 0; t < nrep; t++) {
                double mean = lensum[t] / std::max(1L, loopcnt[t]);
                nloops[t] = std::max(1, (int)std::ceil(M.N / std::max(1.0, mean)));
            }
            std::fill(lensum.begin(), lensum.end(), 0.0); std::fill(loopcnt.begin(), loopcnt.end(), 0L); std::fill(abortcnt.begin(), abortcnt.end(), 0L);
        }
        if (measuring) {
            for (int t = 0; t < nrep; t++) {
                const Ice::Obs& o = obs[t];
                double Q = o.Q, W2x = (double)o.Wx * o.Wx, W2y = (double)o.Wy * o.Wy;
                bool xpar = Q >= 0;
                double Spar0 = xpar ? o.Sx0 : o.Sy0, Spar1 = xpar ? o.Sx1 : o.Sy1, Sper0 = xpar ? o.Sy0 : o.Sx0, Sper1 = xpar ? o.Sy1 : o.Sx1;
                double v[NF] = {(double)o.nc / M.NV, (double)o.ncB / M.NV, Q, Q * Q, Q * Q * Q * Q, std::fabs(Q), W2x, W2y, W2x * W2x, W2y * W2y,
                                xpar ? W2x : W2y, xpar ? W2y : W2x, Spar0, Spar1, Sper0, Sper1, o.Sx0, o.Sx1, o.Sy0, o.Sy1, 0.5 * (W2x + W2y)};
                for (int f = 0; f < NF; f++) acc[t * NF + f] += v[f];
                int c = perm[t];
                int sg = Q > 0.1 ? 1 : Q < -0.1 ? -1 : 0;
                if (sg != 0) { if (sgn[c] != 0 && sg != sgn[c]) switches[c]++; sgn[c] = sg; }
            }
            if (rec % series_every == 0)
                for (int t = 0; t < nrep; t++) {
                    const Ice::Obs& o = obs[t];
                    series.push_back((float)o.Q); series.push_back((float)o.nc / M.NV); series.push_back((float)o.Wx); series.push_back((float)o.Wy);
                }
            rec++; inblock++;
            if (inblock == blen) { for (int i = 0; i < nrep * NF; i++) { blocks.push_back(acc[i] / blen); acc[i] = 0; } inblock = 0; }
        }
        for (int t = (int)(it & 1); t + 1 < nrep; t += 2) {   // replica exchange in u
            double d = (U[t] - U[t + 1]) * (double)(ncnow[t + 1] - ncnow[t]);
            swt[t]++;
            if (d >= 0 || swapRng.uniform() < std::exp(d)) { std::swap(perm[t], perm[t + 1]); swa[t]++; }
        }
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    {
        std::ofstream f(out + ".blocks.csv"); f.precision(9);
        f << "block,t,u";
        for (int k = 0; k < NF; k++) f << "," << names[k];
        f << "\n";
        for (int b = 0; b < nbl; b++) for (int t = 0; t < nrep; t++) {
            f << b << "," << t << "," << U[t];
            for (int k = 0; k < NF; k++) f << "," << blocks[((size_t)b * nrep + t) * NF + k];
            f << "\n";
        }
    }
    { std::ofstream f(out + ".series.f32", std::ios::binary); f.write((const char*)series.data(), series.size() * sizeof(float)); }
    {
        std::ofstream f(out + ".meta.json"); f.precision(9);
        f << "{\n \"model\": \"" << model << "\", \"L\": " << L << ", \"neq\": " << neq << ", \"nmeas\": " << nmeas << ", \"meas_every\": " << meas_every
          << ", \"series_every\": " << series_every << ", \"nblocks\": " << nbl << ", \"blocklen\": " << blen << ", \"seed\": " << seed
          << ", \"init\": \"" << initname << "\", \"commit\": \"" << commit << "\", \"seconds\": " << secs << ",\n \"u\": [";
        for (int t = 0; t < nrep; t++) f << (t ? "," : "") << U[t];
        f << "],\n \"swap_acc\": [";
        for (int t = 0; t + 1 < nrep; t++) f << (t ? "," : "") << (swt[t] ? (double)swa[t] / swt[t] : 0.0);
        f << "],\n \"loops_per_sweep\": [";
        for (int t = 0; t < nrep; t++) f << (t ? "," : "") << nloops[t];
        f << "],\n \"mean_loop_len\": [";
        for (int t = 0; t < nrep; t++) f << (t ? "," : "") << (loopcnt[t] ? lensum[t] / loopcnt[t] : 0.0);
        f << "],\n \"aborted_fraction\": [";
        for (int t = 0; t < nrep; t++) f << (t ? "," : "") << (loopcnt[t] ? (double)abortcnt[t] / loopcnt[t] : 0.0);
        f << "],\n \"Q_sign_switches_per_config\": [";
        for (int c = 0; c < nrep; c++) f << (c ? "," : "") << switches[c];
        f << "]\n}\n";
    }
    double minsw = 1;
    for (int t = 0; t + 1 < nrep; t++) if (swt[t]) minsw = std::min(minsw, (double)swa[t] / swt[t]);
    fprintf(stderr, "done %s: %.1fs, min swap acc %.3f\n", out.c_str(), secs, minsw);
    return 0;
}
