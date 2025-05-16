#pragma once

#include <string>
#include <unordered_map>
#include <fstream>

struct Config {
    int limit;             // stop coarsening when graph size <= limit
    int spectral_runs;     // number of spectral starts
    int spectral_topk;     // number of best spectral starts to select
    int spectral_iters;    // iterations for power method
    int random_restarts;
    int fm_passes;         // FM refinement passes per level
    int sa_steps;          // total SA steps
    double sa_T0;          // initial temperature
    double sa_cooling;     // cooling rate
};

bool LoadConfig(const std::string &path, Config &cfg);
