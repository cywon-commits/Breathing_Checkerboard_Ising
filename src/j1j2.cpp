// Plain square-lattice J1-J2 Ising model with parallel tempering and strip flips (cross-check for task 1, section 8a).
//
//   H = J1 sum_<nn> t_i t_j + J2 sum_<nnn> t_i t_j - hf sum_i t_i        (Lb x Lb torus, Lb even)
//
// When the B plaquettes of the breathing checkerboard model are frozen into superspins t = +-1 the model maps
// exactly onto this one with  J1 = 2 (J_A + eps),  J2 = J_A - 2 eps,  hf = 4 h  (tests/test_main.cpp, T6).
// Updates: single spin Metropolis, strip flips (a block of w consecutive rows/columns, symmetric proposal), PT.
// Output files have the same layout as bcising (blocks.csv with the same columns, meta.json, series.f32 with
// 3 floats [e, mN^2, ms^2] per replica every `series_every` measurements).
//
//   j1j2 --Lb 32 --eps 0.05 [--h 0]            (breathing-model mapping)    or   --J1 1 --J2 0.45 [--hf 0]
#include <chrono>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>

#include "j1j2.hpp"
#include "model.hpp"  // Rng only

struct Rep { std::vector<int8_t> t; double E; Rng rng; };

int main(int argc, char** argv) {
    std::map<std::string, std::string> a;
    for (int i = 1; i + 1 < argc; i += 2) a[argv[i] + 2] = argv[i + 1];
    auto get = [&](const char* k, const char* d) { return a.count(k) ? a[k] : std::string(d); };
    int Lb = std::stoi(get("Lb", "32"));
    double JA = std::stod(get("JA", "1")), eps = std::stod(get("eps", "0.05")), h = std::stod(get("h", "0"));
    double J1 = 2 * (JA + eps), J2 = JA - 2 * eps, hf = 4 * h;
    if (a.count("J1")) J1 = std::stod(a["J1"]);
    if (a.count("J2")) J2 = std::stod(a["J2"]);
    if (a.count("hf")) hf = std::stod(a["hf"]);
    double Tmin = std::stod(get("Tmin", "1")), Tmax = std::stod(get("Tmax", "2"));
    int nrep = std::stoi(get("nrep", "24"));
    long neq = std::stol(get("neq", "3000")), nmeas = std::stol(get("nmeas", "20000"));
    int nblocks = std::stoi(get("nblocks", "20")), meas_every = std::stoi(get("meas_every", "1"));
    int series_every = std::stoi(get("series_every", "10")), n_strip = std::stoi(get("nstrip", "4"));
    uint64_t seed = std::stoull(get("seed", "1"));
    std::string out = get("out", "j1j2run"), initname = get("init", "neel"), commit = get("commit", "unknown");
    if (Lb % 2) { fprintf(stderr, "Lb must be even\n"); return 1; }
    J1J2 M(Lb, J1, J2, hf);

    std::vector<double> T(nrep);
    for (int i = 0; i < nrep; i++) T[i] = nrep == 1 ? Tmin : Tmin * std::pow(Tmax / Tmin, (double)i / (nrep - 1));
    std::vector<Rep> reps(nrep);
    for (int c = 0; c < nrep; c++) {
        Rep& r = reps[c];
        r.rng = Rng(seed, 5000 + c);
        r.t.assign(M.N, 1);
        int kind = initname == "neel" ? 0 : initname == "stripe_x" ? 1 : initname == "stripe_y" ? 2 : initname == "random" ? 3 : 4;
        if (initname == "mixed") kind = c % 4;
        for (int y = 0; y < Lb; y++)
            for (int x = 0; x < Lb; x++) {
                int s = 1;
                if (kind == 0) s = ((x + y) & 1) ? -1 : 1;
                else if (kind == 1) s = (x & 1) ? -1 : 1;
                else if (kind == 2) s = (y & 1) ? -1 : 1;
                else if (kind == 3) s = (r.rng.next() >> 63) ? 1 : -1;
                r.t[M.id(x, y)] = (int8_t)s;
            }
        r.E = M.energy(r.t);
    }
    std::vector<int> perm(nrep);
    for (int i = 0; i < nrep; i++) perm[i] = i;
    Rng swapRng(seed, 888888);

    const int NF = 13;
    long nrec = nmeas / meas_every, blen = std::max(1L, nrec / nblocks);
    int nbl = (int)(nrec / blen);
    std::vector<double> acc(nrep * NF, 0.0), blocks;
    std::vector<float> series;
    std::vector<long> swt(nrep, 0), swa(nrep, 0), tried(nrep * 2, 0), accd(nrep * 2, 0);
    double maxdrift = 0;
    long rec = 0, inblock = 0;
    auto t0 = std::chrono::steady_clock::now();

    struct O { double e, mN, mx, my; };
    for (long it = 0; it < neq + nmeas; it++) {
        bool measuring = it >= neq && ((it - neq) % meas_every == 0) && rec < (long)nbl * blen;
        std::vector<O> obs(nrep);
        for (int i = 0; i < nrep; i++) {
            Rep& r = reps[perm[i]];
            double beta = 1.0 / T[i];
            Rng& g = r.rng;
            for (int k = 0; k < M.N; k++) {
                int x = g.below(Lb), y = g.below(Lb);
                double dE = M.dE_single(r.t, x, y);
                tried[i * 2]++;
                if (dE <= 0 || g.uniform() < std::exp(-beta * dE)) { r.t[M.id(x, y)] = -r.t[M.id(x, y)]; r.E += dE; accd[i * 2]++; }
            }
            for (int k = 0; k < n_strip; k++) {
                int orient = g.below(2), x0 = g.below(Lb), w = 1 + g.below(Lb - 1);
                double dE = M.dE_strip(r.t, orient, x0, w);
                tried[i * 2 + 1]++;
                if (dE <= 0 || g.uniform() < std::exp(-beta * dE)) { M.apply_strip(r.t, orient, x0, w); r.E += dE; accd[i * 2 + 1]++; }
            }
            if (measuring) {
                double mN = 0, mx = 0, my = 0;
                for (int y = 0; y < Lb; y++)
                    for (int x = 0; x < Lb; x++) {
                        int s = r.t[M.id(x, y)];
                        mN += ((x + y) & 1) ? -s : s; mx += (x & 1) ? -s : s; my += (y & 1) ? -s : s;
                    }
                double Eex = M.energy(r.t);
                maxdrift = std::max(maxdrift, std::fabs(Eex - r.E));
                r.E = Eex;
                obs[i] = {Eex / M.N, mN / M.N, mx / M.N, my / M.N};
            }
        }
        if (measuring) {
            for (int i = 0; i < nrep; i++) {
                const O& o = obs[i];
                double ms2 = o.mx * o.mx + o.my * o.my;
                double v[NF] = {o.e, o.e * o.e, o.mN, o.mN * o.mN, o.mN * o.mN * o.mN * o.mN, std::fabs(o.mN), ms2, ms2 * ms2,
                                std::sqrt(ms2), 0, 0, 0, 0};
                for (int f = 0; f < NF; f++) acc[i * NF + f] += v[f];
            }
            if (rec % series_every == 0)
                for (int i = 0; i < nrep; i++) {
                    const O& o = obs[i];
                    series.push_back((float)o.e); series.push_back((float)(o.mN * o.mN));
                    series.push_back((float)(o.mx * o.mx + o.my * o.my));
                }
            rec++; inblock++;
            if (inblock == blen) {
                for (int i = 0; i < nrep * NF; i++) { blocks.push_back(acc[i] / blen); acc[i] = 0; }
                inblock = 0;
            }
        }
        for (int i = (int)(it & 1); i + 1 < nrep; i += 2) {
            double d = (1.0 / T[i] - 1.0 / T[i + 1]) * (reps[perm[i]].E - reps[perm[i + 1]].E);
            swt[i]++;
            if (d >= 0 || swapRng.uniform() < std::exp(d)) { std::swap(perm[i], perm[i + 1]); swa[i]++; }
        }
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const char* names[NF] = {"e", "e2", "mN", "mN2", "mN4", "mNabs", "ms2", "ms4", "msabs", "rdiag", "rnbr", "rviol", "rBal"};
    {
        std::ofstream f(out + ".blocks.csv"); f.precision(10);
        f << "block,t,T";
        for (int k = 0; k < NF; k++) f << "," << names[k];
        f << "\n";
        for (int b = 0; b < nbl; b++)
            for (int i = 0; i < nrep; i++) {
                f << b << "," << i << "," << T[i];
                for (int k = 0; k < NF; k++) f << "," << blocks[((size_t)b * nrep + i) * NF + k];
                f << "\n";
            }
    }
    { std::ofstream f(out + ".series.f32", std::ios::binary); f.write((const char*)series.data(), series.size() * sizeof(float)); }
    {
        std::ofstream f(out + ".meta.json"); f.precision(10);
        f << "{\n \"model\": \"j1j2\", \"Lb\": " << Lb << ", \"J1\": " << J1 << ", \"J2\": " << J2 << ", \"hf\": " << hf
          << ", \"g\": " << J2 / J1 << ", \"eps\": " << eps << ", \"h\": " << h << ",\n \"nrep\": " << nrep << ", \"neq\": " << neq
          << ", \"nmeas\": " << nmeas << ", \"meas_every\": " << meas_every << ", \"nblocks\": " << nbl << ", \"blocklen\": " << blen
          << ", \"series_every\": " << series_every << ", \"series_fields\": [\"e\",\"mN2\",\"ms2\"], \"nstrip\": " << n_strip
          << ",\n \"seed\": " << seed << ", \"init\": \"" << initname << "\", \"commit\": \"" << commit << "\", \"seconds\": " << secs
          << ", \"max_energy_drift\": " << maxdrift << ",\n \"T\": [";
        for (int i = 0; i < nrep; i++) f << (i ? "," : "") << T[i];
        f << "],\n \"swap_acc\": [";
        for (int i = 0; i + 1 < nrep; i++) f << (i ? "," : "") << (swt[i] ? (double)swa[i] / swt[i] : 0.0);
        f << "],\n \"move_acc\": {\"single\": [";
        for (int i = 0; i < nrep; i++) f << (i ? "," : "") << (tried[i * 2] ? (double)accd[i * 2] / tried[i * 2] : 0.0);
        f << "], \"strip\": [";
        for (int i = 0; i < nrep; i++) f << (i ? "," : "") << (tried[i * 2 + 1] ? (double)accd[i * 2 + 1] / tried[i * 2 + 1] : 0.0);
        f << "]}\n}\n";
    }
    double minsw = 1;
    for (int i = 0; i + 1 < nrep; i++) if (swt[i]) minsw = std::min(minsw, (double)swa[i] / swt[i]);
    fprintf(stderr, "done %s: %.1fs, min swap acc %.3f, max E drift %.2e\n", out.c_str(), secs, minsw, maxdrift);
    return 0;
}
