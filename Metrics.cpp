#include "Metrics.hpp"

#include <numeric>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <queue>
#include <iostream>

// Compute cut size (# edges crossing partition)
static int CutSize(const Graph &g, const std::vector<int> &part) {
    int n = g.NumVertices(), cut = 0;
    for (int u = 0; u < n; ++u) {
        for (int v: g.Neighbors(u)) {
            if (u < v && part[u] != part[v]) {
                ++cut;
            }
        }
    }
    return cut;
}

// Heavy-edge matching for coarsening
static std::vector<int> HeavyEdgeMatching(const Graph &g) {
    int n = g.NumVertices();
    std::vector<int> mate(n, -1), order(n);
    std::iota(order.begin(), order.end(), 0);
    std::mt19937 rng(std::random_device{}());
    std::shuffle(order.begin(), order.end(), rng);
    for (int u: order) {
        if (mate[u] == -1) {
            int best = -1, best_deg = -1;
            for (int v: g.Neighbors(u)) {
                if (mate[v] == -1) {
                    int deg = static_cast<int>(g.Neighbors(v).size());
                    if (deg > best_deg) {
                        best_deg = deg;
                        best = v;
                    }
                }
            }
            mate[u] = (best != -1 ? best : u);
            if (best != -1) mate[best] = u;
        }
    }
    return mate;
}

// Build coarse graph by contracting matched pairs
static std::pair<Graph, std::vector<std::vector<int>>> BuildCoarse(
        const Graph &g,
        const std::vector<int> &mate) {
    int n = g.NumVertices();
    std::vector<int> cid(n, -1);
    int cnt = 0;
    for (int u = 0; u < n; ++u) {
        if (u <= mate[u]) {
            cid[u] = cnt;
            if (mate[u] != u) cid[mate[u]] = cnt;
            ++cnt;
        }
    }
    Graph cg(cnt);
    std::vector<std::vector<int>> cmap(cnt);
    std::unordered_set<long long> seen;
    seen.reserve(static_cast<size_t>(n) * 2);
    auto key = [&](int a, int b) {
        return (static_cast<long long>(a) << 32) | static_cast<unsigned>(b);
    };
    for (int u = 0; u < n; ++u) {
        int cu = cid[u];
        cmap[cu].push_back(u);
        for (int v: g.Neighbors(u)) {
            int cv = cid[v];
            if (cu == cv) continue;
            int a = std::min(cu, cv);
            int b = std::max(cu, cv);
            long long k = key(a, b);
            if (seen.insert(k).second) {
                cg.AddEdge(a, b);
            }
        }
    }
    return {std::move(cg), std::move(cmap)};
}

// Naive Kernighan-Lin refinement (FM pass)
static void RefineNaive(
        const Graph &g,
        std::vector<int> &part,
        int passes) {
    int n = g.NumVertices();
    auto Gain = [&](int u) {
        int same = 0, diff = 0;
        for (int v: g.Neighbors(u)) {
            if (part[v] == part[u]) ++same;
            else ++diff;
        }
        return diff - same;
    };
    for (int it = 0; it < passes; ++it) {
        std::vector<char> locked(n, 0);
        int moved = 0;
        while (moved < n) {
            int bestu = -1;
            int bestg = std::numeric_limits<int>::min();
            for (int u = 0; u < n; ++u) {
                if (locked[u]) continue;
                int g0 = Gain(u);
                int to = part[u] ^ 1;
                int cto = static_cast<int>(
                        std::count(part.begin(), part.end(), to));
                int cfrom = n - cto;
                if (cto + 1 > (n + 1) / 2 || cfrom - 1 < n / 2) continue;
                if (g0 > bestg) {
                    bestg = g0;
                    bestu = u;
                }
            }
            if (bestu < 0 || bestg <= 0) break;
            locked[bestu] = 1;
            part[bestu] ^= 1;
            ++moved;
        }
    }
}

// Spectral initialization (power iteration)
static std::vector<double> SpectralVector(const Graph &g, int iters) {
    int n = g.NumVertices();
    std::vector<double> x(n), y(n);
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (double &xi: x) xi = dist(rng);
    for (int it = 0; it < iters; ++it) {
        for (int u = 0; u < n; ++u) {
            double sum = g.Neighbors(u).size() * x[u];
            for (int v: g.Neighbors(u)) sum -= x[v];
            y[u] = sum;
        }
        double mean = std::accumulate(y.begin(), y.end(), 0.0) / n;
        double norm = 0;
        for (double &yi: y) {
            yi -= mean;
            norm += yi * yi;
        }
        norm = std::sqrt(norm);
        for (int i = 0; i < n; ++i) x[i] = (norm > 1e-14 ? y[i] / norm : y[i]);
    }
    return x;
}


// Simulated Annealing for final refinement
static void SimulatedAnnealing(
        const Graph &g,
        std::vector<int> &part,
        const Config &cfg) {
    int n = g.NumVertices();
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    double T = cfg.sa_T0;
    int count0 = static_cast<int>(std::count(part.begin(), part.end(), 0));
    int count1 = n - count0;
    for (int step = 0; step < cfg.sa_steps; ++step) {
        int u = rng() % n;
        int from = part[u], to = from ^ 1;
        if ((to == 1 && count1 + 1 > (n + 1) / 2) || (to == 0 && count0 + 1 > n / 2)) continue;
        int same = 0, diff = 0;
        for (int v: g.Neighbors(u)) (part[v] == from ? ++same : ++diff);
        int delta = diff - same;  // maximize gain for min-cut is diff-same <0, but we keep this
        if (delta > 0 || ud(rng) < std::exp(-delta / T)) {
            part[u] = to;
            if (to == 1) {
                ++count1;
                --count0;
            } else {
                ++count0;
                --count1;
            }
        }
        T *= cfg.sa_cooling;
    }
}

// Multilevel bisection (min-cut only)
static std::vector<int> MultilevelBisection(
        const Graph &g_input, const Config &cfg) {
    Graph g = g_input;
    std::vector<std::pair<std::vector<std::vector<int>>, Graph>> stk;
    // Coarsen
    while (g.NumVertices() > cfg.limit) {
        auto mate = HeavyEdgeMatching(g);
        auto coarse = BuildCoarse(g, mate);
        stk.emplace_back(std::move(coarse.second), std::move(g));
        g = std::move(coarse.first);
    }
    int cn = g.NumVertices();
    // Spectral + random candidates
    std::vector<std::pair<int, std::vector<int>>> cand;
    cand.reserve(cfg.spectral_runs + cfg.random_restarts);
    // spectral starts
    for (int sr = 0; sr < cfg.spectral_runs; ++sr) {
        auto vec = SpectralVector(g, cfg.spectral_iters);
        std::vector<int> idx(cn);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return vec[a] < vec[b]; });
        std::vector<int> part(cn, 0);
        for (int i = cn / 2; i < cn; ++i) part[idx[i]] = 1;
        RefineNaive(g, part, cfg.fm_passes);
        cand.emplace_back(CutSize(g, part), part);
    }
    // random restarts
    {
        std::mt19937 rng(std::random_device{}());
        for (int rr = 0; rr < cfg.random_restarts; ++rr) {
            std::vector<int> idx(cn);
            std::iota(idx.begin(), idx.end(), 0);
            std::shuffle(idx.begin(), idx.end(), rng);
            std::vector<int> part(cn, 0);
            for (int i = 0; i < cn / 2; ++i) part[idx[i]] = 1;
            RefineNaive(g, part, cfg.fm_passes);
            cand.emplace_back(CutSize(g, part), part);
        }
    }
    // pick best
    std::sort(cand.begin(), cand.end(), [](auto &a, auto &b) { return a.first < b.first; });
    int topk = std::min(cfg.spectral_topk, static_cast<int>(cand.size()));
    int bestc = std::numeric_limits<int>::max();
    std::vector<int> bestp;
    for (int i = 0; i < topk; ++i) {
        int c0 = cand[i].first;
        if (c0 < bestc) {
            bestc = c0;
            bestp = cand[i].second;
        }
    }
    // uncoarsen + refine + SA
    std::vector<int> part = bestp;
    while (!stk.empty()) {
        auto pr = std::move(stk.back());
        stk.pop_back();
        auto &cmap = pr.first;
        Graph fine = std::move(pr.second);
        int fn = fine.NumVertices();
        std::vector<int> fpart(fn);
        for (size_t cv = 0; cv < cmap.size(); ++cv) {
            for (int v: cmap[cv]) fpart[v] = part[cv];
        }
        RefineNaive(fine, fpart, cfg.fm_passes);
        part.swap(fpart);
    }
    SimulatedAnnealing(g_input, part, cfg);
    return part;
}

Metrics ComputeMetrics(const Graph &g, const Config &cfg) {
    int n = g.NumVertices();
    // Degrees
    std::vector<int> degree(n);
    for (int u = 0; u < n; ++u) {
        degree[u] = static_cast<int>(g.Neighbors(u).size());
    }
    int min_deg = *std::min_element(degree.begin(), degree.end());
    int max_deg = *std::max_element(degree.begin(), degree.end());

    // Packing density: N_total/N_compute * sum(eccentricity[u] * degree[u])
    // accumulate pack_sum over ALL vertices, not just compute-nodes
    double pack_sum = 0.0;
    auto compute_nodes = g.ComputeNodes();
    int n_compute = static_cast<int>(std::count(compute_nodes.begin(), compute_nodes.end(), true));

    // For path metrics
    long long total_dist = 0;
    long long total_pairs = 0;
    double total_paths = 0;

    int global_diameter = 0;

    // For each source, run BFS to get distances and path counts
    for (int src = 0; src < n; ++src) {
        // BFS
        std::vector<int> dist(n, std::numeric_limits<int>::max());
        std::vector<uint64_t> path_count(n, 0);
        std::queue<int> q;

        dist[src] = 0;
        path_count[src] = 1;
        q.push(src);

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            int du = dist[u];
            for (int v: g.Neighbors(u)) {
                if (dist[v] == std::numeric_limits<int>::max()) {
                    // first visit
                    dist[v] = du + 1;
                    path_count[v] = path_count[u];
                    q.push(v);
                } else if (dist[v] == du + 1) {
                    // another shortest path
                    path_count[v] += path_count[u];
                }
            }
        }

        // Eccentricity = max finite dist
        int ecc = 0;
        for (int v = 0; v < n; ++v) {
            if (dist[v] < std::numeric_limits<int>::max()) {
                ecc = std::max(ecc, dist[v]);
            }
            else
            {
                std::cerr << "Graph is unconnected\n";
                exit(1);
            }
        }
        global_diameter = std::max(global_diameter, ecc);

        // accumulate cost for all nodes
        if (degree[src] == 1 && compute_nodes[src] && !compute_nodes[g.Neighbors(src)[0]])
            pack_sum += static_cast<double>(ecc) * degree[g.Neighbors(src)[0]];
        else
            pack_sum += static_cast<double>(ecc) * degree[src];

        // Accumulate distances and path counts for pairs (src, v)
        for (int v = 0; v < n; ++v) {
            if (v == src) continue;
            if (dist[v] < std::numeric_limits<int>::max()) {
                total_dist += dist[v];
                total_paths += path_count[v];
                ++total_pairs;
            }
        }
    }

    Metrics m;
    m.nodes_count = n;
    m.compute_nodes_count = n_compute;
    m.diameter = global_diameter;
    // average over ordered pairs: divide by total_pairs
    m.average_path_length = total_pairs > 0
                            ? static_cast<double>(total_dist) / total_pairs
                            : 0.0;
    m.min_degree = min_deg;
    m.max_degree = max_deg;

    // global packing density
    m.global_packing_density = static_cast<double>(m.nodes_count) / (global_diameter * max_deg);
    // Packing density normalized by sum of (ecc * degree)
    if (n_compute > 0 && pack_sum != 0.0) {
        // Packing density based on Moore bound: (D_min * Delta) / sum_ecc_deg
        int Delta = max_deg;
        int D_min = 1;
        if (Delta > 1) {
            unsigned long long moore = 1ULL;   // 1 + Δ * Σ(Δ-1)^i
            unsigned long long power = 1ULL;   // (Δ-1)^0
            for (D_min = 1;; ++D_min) {
                if (D_min == 1) {
                    moore = 1ULL + static_cast<unsigned long long>(Delta) * power;
                } else {
                    power *= static_cast<unsigned long long>(Delta - 1);
                    moore += static_cast<unsigned long long>(Delta) * power;
                }
                if (moore >= static_cast<unsigned long long>(n)) break;
            }
        }
        if (pack_sum > 0.0) {
            // normalization includes compute node fraction
            m.norm_local_packing_density =
                    (static_cast<double>(D_min) * Delta / pack_sum) * (static_cast<double>(n_compute));
        } else {
            m.norm_local_packing_density = 0.0;
        };
    } else {
        m.norm_local_packing_density = 0.0;
    }
    m.average_num_shortest_paths = total_pairs > 0
                                   ? static_cast<double>(total_paths) / total_pairs
                                   : 0.0;

    m.norm_average_num_shortest_paths = total_pairs > 0
                                        ? static_cast<double>(total_paths) / (static_cast<double>(total_pairs) * global_diameter)
                                        : 0.0;

    auto part = MultilevelBisection(g, cfg);
    m.bisection_width = CutSize(g, part);
    return m;
}
