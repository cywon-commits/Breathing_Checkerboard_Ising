// Parallel-tempering Monte Carlo driver for the breathing checkerboard Ising model.
//
// usage: bcising --Lb 16 --K 4 --h 0 --eps 0.012 --Tmin 0.7 --Tmax 1.4 --nrep 24
//                --neq 20000 --nmeas 200000 --nblocks 50 --seed 1 --init neel --out prefix
// Writes  <out>.blocks.csv   block averages of all observables per temperature
//         <out>.series.f32   decimated time series (float32, shape [nrec][nrep][2]: mN^2, ms^2)
//         <out>.meta.json    parameters, provenance, acceptance rates
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "model.hpp"

int main(int argc, char** argv) {
    std::map<std::string, std::string> a;
    for (int i = 1; i + 1 < argc; i += 2) a[argv[i] + 2] = argv[i + 1];
    auto get = [&](const char* k, const char* d) { return a.count(k) ? a[k] : std::string(d); };
    int Lb = std::stoi(get("Lb", "16"));
    double JA = std::stod(get("JA", "1")), K = std::stod(get("K", "4")), h = std::stod(get("h", "0"));
    double eps = std::stod(get("eps", "0.01"));
    double Tmin = std::stod(get("Tmin", "0.7")), Tmax = std::stod(get("Tmax", "1.4"));
    int nrep = std::stoi(get("nrep", "24"));
    long neq = std::stol(get("neq", "20000")), nmeas = std::stol(get("nmeas", "200000"));
    int nblocks = std::stoi(get("nblocks", "50"));
    int meas_every = std::stoi(get("meas_every", "1"));
    int series_every = std::stoi(get("series_every", "10"));
    int n_strip = std::stoi(get("nstrip", "4"));
    unsigned moves = (unsigned)std::stoi(get("moves", "15"));
    uint64_t seed = std::stoull(get("seed", "1"));
    std::string initname = get("init", "neel"), out = get("out", "run"), commit = get("commit", "unknown");
    std::string tlist = get("Tlist", "");

    if (Lb % 2) { fprintf(stderr, "Lb must be even (L multiple of 4)\n"); return 1; }
    Model M(Lb, JA, eps, K, h);

    std::vector<double> T(nrep);
    if (!tlist.empty()) {  // comma separated explicit ladder
        nrep = 0; T.clear();
        size_t p = 0;
        while (p < tlist.size()) {
            size_t q = tlist.find(',', p); if (q == std::string::npos) q = tlist.size();
            T.push_back(std::stod(tlist.substr(p, q - p))); p = q + 1;
        }
        nrep = (int)T.size();
    } else {
        for (int t = 0; t < nrep; t++)
            T[t] = nrep == 1 ? Tmin : Tmin * std::pow(Tmax / Tmin, (double)t / (nrep - 1));
    }

    std::vector<Replica> reps(nrep);
    for (int c = 0; c < nrep; c++) {
        int kind = 0;
        if (initname == "neel") kind = 0;
        else if (initname == "stripe_x") kind = 1;
        else if (initname == "stripe_y") kind = 2;
        else if (initname == "random") kind = 3;
        else if (initname == "ferro") kind = 4;
        else if (initname == "mixed") kind = c % 4;
        M.init(reps[c], kind, seed, 1000 + c);
    }
    std::vector<int> perm(nrep);
    for (int t = 0; t < nrep; t++) perm[t] = t;
    Rng swapRng(seed, 999999);

    const int NF = 13;  // number of recorded observables
    long nrec = nmeas / meas_every;
    long blen = std::max(1L, nrec / nblocks);
    int nbl = (int)(nrec / blen);
    std::vector<double> acc(nrep * NF, 0.0);
    std::vector<double> blocks;  // nbl * nrep * NF
    std::vector<float> series;
    std::vector<long> mtried(nrep * N_MOVES, 0), macc(nrep * N_MOVES, 0), swt(nrep, 0), swa(nrep, 0);
    double maxdrift = 0;
    long rec = 0, inblock = 0;
    auto t0 = std::chrono::steady_clock::now();

    long total = neq + nmeas;
    for (long it = 0; it < total; it++) {
        bool measuring = it >= neq && ((it - neq) % meas_every == 0) && rec < (long)nbl * blen;
        std::vector<Obs> obs(nrep);
#pragma omp parallel for schedule(static)
        for (int t = 0; t < nrep; t++) {
            Replica& r = reps[perm[t]];
            M.sweep(r, 1.0 / T[t], n_strip, &mtried[t * N_MOVES], &macc[t * N_MOVES], moves);
            if (measuring) obs[t] = M.measure(r);
        }
        if (measuring) {
            for (int t = 0; t < nrep; t++) {
                Replica& r = reps[perm[t]];
                double Eex = obs[t].E * M.N;
                maxdrift = std::max(maxdrift, std::fabs(Eex - r.E));
                r.E = Eex;
                const Obs& o = obs[t];
                double ms2 = o.mx * o.mx + o.my * o.my;
                double v[NF] = {o.E, o.E * o.E, o.mN, o.mN * o.mN, o.mN * o.mN * o.mN * o.mN, std::fabs(o.mN),
                                ms2, ms2 * ms2, std::sqrt(ms2), o.rdiag, o.rnbr, o.rviol, o.rBal};
                for (int f = 0; f < NF; f++) acc[t * NF + f] += v[f];
            }
            if (rec % series_every == 0)
                for (int t = 0; t < nrep; t++) {
                    const Obs& o = obs[t];
                    series.push_back((float)(o.mN * o.mN));
                    series.push_back((float)(o.mx * o.mx + o.my * o.my));
                }
            rec++; inblock++;
            if (inblock == blen) {
                for (int i = 0; i < nrep * NF; i++) { blocks.push_back(acc[i] / blen); acc[i] = 0; }
                inblock = 0;
            }
        }
        // replica exchange (alternate even / odd pairs)
        for (int t = (int)(it & 1); t + 1 < nrep; t += 2) {
            double d = (1.0 / T[t] - 1.0 / T[t + 1]) * (reps[perm[t]].E - reps[perm[t + 1]].E);
            swt[t]++;
            if (d >= 0 || swapRng.uniform() < std::exp(d)) { std::swap(perm[t], perm[t + 1]); swa[t]++; }
        }
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    const char* names[NF] = {"e", "e2", "mN", "mN2", "mN4", "mNabs", "ms2", "ms4", "msabs",
                             "rdiag", "rnbr", "rviol", "rBal"};
    {
        std::ofstream f(out + ".blocks.csv");
        f.precision(10);
        f << "block,t,T";
        for (int k = 0; k < NF; k++) f << "," << names[k];
        f << "\n";
        for (int b = 0; b < nbl; b++)
            for (int t = 0; t < nrep; t++) {
                f << b << "," << t << "," << T[t];
                for (int k = 0; k < NF; k++) f << "," << blocks[((size_t)b * nrep + t) * NF + k];
                f << "\n";
            }
    }
    {
        std::ofstream f(out + ".series.f32", std::ios::binary);
        f.write((const char*)series.data(), series.size() * sizeof(float));
    }
    {
        std::ofstream f(out + ".meta.json");
        f.precision(10);
        f << "{\n \"Lb\": " << Lb << ", \"L\": " << M.L << ", \"JA\": " << JA << ", \"K\": " << K << ", \"h\": " << h
          << ", \"eps\": " << eps << ",\n \"nrep\": " << nrep << ", \"neq\": " << neq << ", \"nmeas\": " << nmeas
          << ", \"meas_every\": " << meas_every << ", \"nblocks\": " << nbl << ", \"blocklen\": " << blen
          << ", \"series_every\": " << series_every << ", \"nstrip\": " << n_strip << ", \"moves\": " << moves
          << ",\n \"seed\": " << seed << ", \"init\": \"" << initname << "\", \"commit\": \"" << commit
          << "\", \"seconds\": " << secs << ", \"max_energy_drift\": " << maxdrift << ",\n \"T\": [";
        for (int t = 0; t < nrep; t++) f << (t ? "," : "") << T[t];
        f << "],\n \"swap_acc\": [";
        for (int t = 0; t + 1 < nrep; t++) f << (t ? "," : "") << (swt[t] ? (double)swa[t] / swt[t] : 0.0);
        f << "],\n \"move_acc\": {";
        const char* mn[N_MOVES] = {"single", "half", "super", "strip"};
        for (int m = 0; m < N_MOVES; m++) {
            f << (m ? "," : "") << "\"" << mn[m] << "\": [";
            for (int t = 0; t < nrep; t++)
                f << (t ? "," : "") << (mtried[t * N_MOVES + m] ? (double)macc[t * N_MOVES + m] / mtried[t * N_MOVES + m] : 0.0);
            f << "]";
        }
        f << "}\n}\n";
    }
    double minsw = 1;
    for (int t = 0; t + 1 < nrep; t++) if (swt[t]) minsw = std::min(minsw, (double)swa[t] / swt[t]);
    fprintf(stderr, "done %s: %.1fs, min swap acc %.3f, max E drift %.2e\n", out.c_str(), secs, minsw, maxdrift);
    return 0;
}
