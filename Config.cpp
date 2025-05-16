#include "Config.hpp"

static std::string Trim(const std::string &s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return std::string();
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

bool LoadConfig(const std::string &path, Config &cfg) {
    std::ifstream in(path);
    if (!in) return false;
    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        size_t pos = line.find('=');
        if (pos == std::string::npos) continue;
        std::string key = Trim(line.substr(0, pos));
        std::string val = Trim(line.substr(pos + 1));
        kv[key] = val;
    }
    try {
        cfg.limit = std::stoi(kv.at("limit"));
        cfg.spectral_runs = std::stoi(kv.at("spectral_runs"));
        cfg.spectral_topk = std::stoi(kv.at("spectral_topk"));
        cfg.spectral_iters = std::stoi(kv.at("spectral_iters"));
        cfg.random_restarts = std::stoi(kv.at("random_restarts"));
        cfg.fm_passes = std::stoi(kv.at("fm_passes"));
        cfg.sa_steps = std::stoi(kv.at("sa_steps"));
        cfg.sa_T0 = std::stod(kv.at("sa_T0"));
        cfg.sa_cooling = std::stod(kv.at("sa_cooling"));
    } catch (...) {
        return false;
    }
    return true;
}
