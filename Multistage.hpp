#pragma once

#include "Graph.hpp"

Graph GenerateDelta(int num_cores) {
    if (num_cores <= 0 || (num_cores & (num_cores - 1)) != 0) {
        return Graph(0);
    }

    int stages = 0;
    while ((1 << stages) < num_cores) {
        ++stages;
    }
    int switches_per_stage = num_cores / 2;
    int total_switches = switches_per_stage * stages;
    int total_vertices = num_cores + total_switches;

    Graph graph(total_vertices);
    std::vector<bool> compute_flags(total_vertices, true);

    for (int stage = 0; stage < stages; ++stage) {
        for (int row = 0; row < switches_per_stage; ++row) {
            int sw = num_cores
                     + stage * switches_per_stage
                     + row;
            compute_flags[sw] = false;

            if (stage == 0) {
                for (int p = 0; p < 2; ++p) {
                    int pe = row * 2 + p;
                    graph.AddEdge(pe, sw);
                }
            } else {
                int prev1 = num_cores
                            + (stage - 1) * switches_per_stage
                            + row;
                graph.AddEdge(sw, prev1);

                int bitpos = (stages - 1 - stage);
                int m = row ^ (1 << bitpos);
                int prev2 = num_cores
                            + (stage - 1) * switches_per_stage
                            + m;
                graph.AddEdge(sw, prev2);
            }

            if (stage == stages - 1) {
                for (int p = 0; p < 2; ++p) {
                    int pe = row * 2 + p;
                    graph.AddEdge(sw, pe);
                }
            }
        }
    }

    graph.SetCompute(compute_flags);
    return graph;
}

// 1) Banyan: R(i,j) → R(i+1, 2*j % S) и R(i+1, (2*j+1) % S)
//    классический banyan с наращиванием по веткам
//    C. Franklin, “VLSI Performance of Banyan…,” Fig. 2(b) :contentReference[oaicite:0]{index=0}:contentReference[oaicite:1]{index=1}
Graph GenerateBanyan(int num_cores) {
    if (num_cores <= 0 || (num_cores & (num_cores - 1)) != 0) return Graph(0);
    int stages = 0;
    while ((1 << stages) < num_cores) ++stages;
    int switches_per_stage = num_cores / 2;
    int total_switches = switches_per_stage * stages;
    int total_vertices = num_cores + total_switches;

    Graph graph(total_vertices);
    std::vector<bool> compute_flags(total_vertices, true);

    // PE ↔ первый and последний stage
    for (int stage = 0; stage < stages; ++stage) {
        for (int row = 0; row < switches_per_stage; ++row) {
            int sw = num_cores + stage * switches_per_stage + row;
            compute_flags[sw] = false;
            if (stage == 0) {
                // входы PE → switch
                for (int p = 0; p < 2; ++p) {
                    graph.AddEdge(row * 2 + p, sw);
                }
            }
            if (stage < stages - 1) {
                // Banyan-соединения
                int child0 = (2 * row) % switches_per_stage;
                int child1 = (2 * row + 1) % switches_per_stage;
                int sw0 = num_cores + (stage + 1) * switches_per_stage + child0;
                int sw1 = num_cores + (stage + 1) * switches_per_stage + child1;
                graph.AddEdge(sw, sw0);
                graph.AddEdge(sw, sw1);
            }
            if (stage == stages - 1) {
                // выходы switch → PE
                for (int p = 0; p < 2; ++p) {
                    graph.AddEdge(sw, row * 2 + p);
                }
            }
        }
    }

    graph.SetCompute(compute_flags);
    return graph;
}

// 2) Butterfly: R(i,j) → R(i+1, j) и R(i+1, j XOR (1<<i))
//    классический butterfly, выбор по i-му биту адреса :contentReference[oaicite:2]{index=2}:contentReference[oaicite:3]{index=3}
Graph GenerateButterfly(int num_cores) {
    if (num_cores <= 0 || (num_cores & (num_cores - 1)) != 0) return Graph(0);
    int stages = 0;
    while ((1 << stages) < num_cores) ++stages;
    int switches_per_stage = num_cores / 2;
    int total_switches = switches_per_stage * stages;
    int total_vertices = num_cores + total_switches;

    Graph graph(total_vertices);
    std::vector<bool> compute_flags(total_vertices, true);

    for (int stage = 0; stage < stages; ++stage) {
        for (int row = 0; row < switches_per_stage; ++row) {
            int sw = num_cores + stage * switches_per_stage + row;
            compute_flags[sw] = false;
            if (stage == 0) {
                for (int p = 0; p < 2; ++p)
                    graph.AddEdge(row * 2 + p, sw);
            }
            if (stage < stages - 1) {
                int sw0 = num_cores + (stage + 1) * switches_per_stage + row;
                int sw1 = num_cores + (stage + 1) * switches_per_stage
                          + (row ^ (1 << stage));
                graph.AddEdge(sw, sw0);
                graph.AddEdge(sw, sw1);
            }
            if (stage == stages - 1) {
                for (int p = 0; p < 2; ++p)
                    graph.AddEdge(sw, row * 2 + p);
            }
        }
    }

    graph.SetCompute(compute_flags);
    return graph;
}

// 3) Omega: межэтапный «perfect shuffle» + swap-1 LSB
//    shuffle: left-rotate индекс j на 1 бит :contentReference[oaicite:4]{index=4}:contentReference[oaicite:5]{index=5}
Graph GenerateOmega(int num_cores) {
    if (num_cores <= 0 || (num_cores & (num_cores - 1)) != 0) return Graph(0);
    int stages = 0;
    while ((1 << stages) < num_cores) ++stages;
    int switches_per_stage = num_cores / 2;
    int total_switches = switches_per_stage * stages;
    int total_vertices = num_cores + total_switches;
    int k = stages;  // log2(num_cores)

    Graph graph(total_vertices);
    std::vector<bool> compute_flags(total_vertices, true);

    for (int stage = 0; stage < stages; ++stage) {
        for (int row = 0; row < switches_per_stage; ++row) {
            int sw = num_cores + stage * switches_per_stage + row;
            compute_flags[sw] = false;
            if (stage == 0) {
                for (int p = 0; p < 2; ++p)
                    graph.AddEdge(row * 2 + p, sw);
            }
            if (stage < stages - 1) {
                // perfect shuffle
                int mask = (1 << k) - 1;
                int shuffled = ((row << 1) & mask) | (row >> (k - 1));
                int sw0 = num_cores + (stage + 1) * switches_per_stage + (shuffled % switches_per_stage);
                int sw1 = num_cores + (stage + 1) * switches_per_stage + ((shuffled ^ 1) % switches_per_stage);
                graph.AddEdge(sw, sw0);
                graph.AddEdge(sw, sw1);
            }
            if (stage == stages - 1) {
                for (int p = 0; p < 2; ++p)
                    graph.AddEdge(sw, row * 2 + p);
            }
        }
    }

    graph.SetCompute(compute_flags);
    return graph;
}

Graph GenerateFlip(int num_cores) {
    if (num_cores <= 0 || (num_cores & (num_cores - 1)) != 0)
        return Graph(0);

    /* основные параметры */
    int stages = 0;
    while ((1 << stages) < num_cores) ++stages;   // k = log2 N
    int switches_per_stage = num_cores / 2;       // N / 2
    int total_switches = switches_per_stage * stages;
    int total_vertices = num_cores + total_switches;

    Graph graph(total_vertices);
    std::vector<bool> compute(total_vertices, true);   // PE = true

    auto sw_idx = [&](int stage, int row) {
        return num_cores + stage * switches_per_stage + row;
    };

    /* помечаем коммутаторы */
    for (int st = 0; st < stages; ++st)
        for (int row = 0; row < switches_per_stage; ++row)
            compute[sw_idx(st, row)] = false;

    /* подключаем PE ↔ первый / последний stage */
    for (int row = 0; row < switches_per_stage; ++row) {
        int pe0 = row * 2;
        int pe1 = row * 2 + 1;
        int sw_first = sw_idx(0, row);
        int sw_last = sw_idx(stages - 1, row);
        graph.AddEdge(pe0, sw_first);
        graph.AddEdge(pe1, sw_first);
        graph.AddEdge(sw_last, pe0);
        graph.AddEdge(sw_last, pe1);
    }

    /* inverse perfect shuffle между стадиями */
    int k = stages;                 // число бит полного адреса линии
    int row_bits = k - 1;           // столько бит кодирует row (0..N/2−1)
    int msb_mask = 1 << (row_bits - 1);

    for (int st = 0; st < stages - 1; ++st) {
        for (int row = 0; row < switches_per_stage; ++row) {
            int cur = sw_idx(st, row);

            for (int port = 0; port < 2; ++port) {
                /* формула  next_row = port·2^{k-2}  |  (row >> 1) */
                int next_row = ((row >> 1) & (msb_mask - 1))
                               | (port << (row_bits - 1));
                int nxt = sw_idx(st + 1, next_row);
                graph.AddEdge(cur, nxt);
            }
        }
    }

    graph.SetCompute(compute);
    return graph;
}


Graph GenerateBenes(int num_terminals) {
    auto is_pow2 = [](int x) { return x > 0 && (x & (x - 1)) == 0; };
    if (!is_pow2(num_terminals) || num_terminals < 2) return Graph(0);

    int m = std::log2(num_terminals);      // log2 N
    int stages = 2 * m - 1;                     // S
    int sw_per = num_terminals / 2;             // 2×2 switches per stage
    int total_sw = stages * sw_per;
    int first_sw = num_terminals;                 // индекс первого переключателя
    int total_v = num_terminals + total_sw;      // терминалы + свитчи

    Graph g(total_v);
    std::vector<bool> is_comp(total_v, true);    // PE = true

    auto rot_left = [&](int w) { return ((w << 1) & (num_terminals - 1)) | (w >> (m - 1)); };
    auto rot_right = [&](int w) { return (w >> 1) | ((w & 1) << (m - 1)); };

    /* 1) входные терминалы → стадия 0 */
    for (int t = 0; t < num_terminals; ++t) {
        int sw = first_sw + t / 2;
        g.AddEdge(t, sw);
    }

    /* 2) соединения между стадиями */
    for (int s = 0; s < stages - 1; ++s) {
        int base_cur = first_sw + s * sw_per;
        int base_next = first_sw + (s + 1) * sw_per;

        bool left_half = (s < m - 1);
        bool middle = (s == m - 1);

        for (int sw = 0; sw < sw_per; ++sw) {
            int self = base_cur + sw;
            is_comp[self] = false;

            for (int port = 0; port < 2; ++port) {
                int wire = sw * 2 + port;               // глобальный индекс провода
                int wire_next = middle ? wire
                                       : left_half ? rot_left(wire)
                                                   : rot_right(wire);

                int sw_next = wire_next / 2;
                int neigh = base_next + sw_next;
                g.AddEdge(self, neigh);
            }
        }
    }

    /* 3) последняя стадия → выходные терминалы (идентичны входным) */
    int last_base = first_sw + (stages - 1) * sw_per;
    for (int sw = 0; sw < sw_per; ++sw) {
        int self = last_base + sw;
        is_comp[self] = false;

        for (int port = 0; port < 2; ++port) {
            int term = sw * 2 + port;            // тот же индекс, что вход
            g.AddEdge(self, term);
        }
    }

    /* 4) помечаем все переключатели как non-compute */
    for (int v = first_sw; v < total_v; ++v) is_comp[v] = false;
    g.SetCompute(is_comp);
    return g;
}


Graph GenerateFlattenedButterfly(int radix, int layers) {
    if (radix < 2 || layers < 2) {
        return Graph(0);
    }

    /* ---------- вычисляем размеры ---------- */
    long long num_pe = 1;
    for (int i = 0; i < layers; ++i) {
        num_pe *= radix;                        // k^n
        if (num_pe > (1LL << 30)) {         // грубая защита от переполнения
            return Graph(0);
        }
    }
    long long num_rt = 1;
    for (int i = 0; i < layers - 1; ++i) {
        num_rt *= radix;                        // k^(n-1)
    }
    long long total_vertices = num_pe + num_rt;

    Graph graph(static_cast<int>(total_vertices));
    std::vector<bool> compute_flags(total_vertices, true);      // PE-узлы = true

    /* ---------- pre-compute strides for k-ичное разложение ---------- */
    std::vector<long long> stride(layers - 1);
    stride[layers - 2] = 1;
    for (int d = layers - 3; d >= 0; --d) {
        stride[d] = stride[d + 1] * radix;
    }

    /* ---------- маршрутизаторы и их связи ---------- */
    for (long long r_idx = 0; r_idx < num_rt; ++r_idx) {
        int router_v = static_cast<int>(num_pe + r_idx);
        compute_flags[router_v] = false;                    // router — не compute

        /* ----- связь router ↔ local PE узлы ----- */
        long long pe_base = r_idx * radix;                      // x0 = 0
        for (int m = 0; m < radix; ++m) {
            int pe_v = static_cast<int>(pe_base + m);
            graph.AddEdge(router_v, pe_v);
        }

        /* ----- восстановить координаты маршрутизатора ----- */
        long long rem = r_idx;
        std::vector<int> coord(layers - 1);
        for (int d = 0; d < layers - 1; ++d) {
            coord[d] = static_cast<int>(rem / stride[d]);
            rem %= stride[d];
        }

        /* ----- меж-роутерные связи ----- */
        for (int dim = 0; dim < layers - 1; ++dim) {
            int old_val = coord[dim];
            for (int val = 0; val < radix; ++val) {
                if (val == old_val) continue;
                long long neighbor_idx = r_idx + (val - old_val) * stride[dim];
                int neighbor_v = static_cast<int>(num_pe + neighbor_idx);
                if (router_v < neighbor_v) {
                    graph.AddEdge(router_v, neighbor_v);
                }
            }
        }
    }

    graph.SetCompute(compute_flags);
    return graph;
}

// Butterfly-Fat-Tree
//  ─ лист (lvl-0):  к каждому роутеру подключены radix PE
//  ─ lvl-0 → lvl-1: каждый роутер имеет radix/2 uplink-ов подряд
//  ─ lvl ≥ 1:   child r соединён ровно с radix/2 роутерами lvl+1,
//               индексы которых имеют тот же остаток (r mod shrink),
//               где shrink = radix / 2
//
// radix  (чётное ≥ 4)  — число портов роутера
// levels (≥ 2)         — число уровней роутеров
Graph GenerateBFT(int radix, int levels)
{
    if (radix < 2 || (radix & 1) || levels < 2)
        return Graph(0);

    const int shrink = radix / 2;              // k/2 uplink-ов

    /* ---------- 1. PE-узлы ---------- */
    long long pe_cnt = 1;
    for (int i = 0; i < levels; ++i) {
        pe_cnt *= radix;                       // radix^levels
        if (pe_cnt > (1LL << 30)) return Graph(0);
    }

    /* ---------- 2. routers per level ---------- */
    std::vector<long long> rt(levels);
    rt[0] = 1;
    for (int i = 0; i < levels - 1; ++i) rt[0] *= radix;   // radix^(L-1)
    for (int l = 1; l < levels; ++l) rt[l] = rt[l - 1] / shrink;

    long long rt_total = 0;
    for (auto v : rt) rt_total += v;
    long long total_v = pe_cnt + rt_total;

    Graph g((int)total_v);
    std::vector<bool> is_comp(total_v, true);

    /* ---------- 3. base-индексы уровней ---------- */
    std::vector<long long> base(levels);
    base[0] = pe_cnt;
    for (int l = 1; l < levels; ++l) base[l] = base[l - 1] + rt[l - 1];

    /* ---------- 4. PE → leaf-роутеры ---------- */
    for (long long pe = 0; pe < pe_cnt; ++pe) {
        long long leaf = base[0] + pe / radix;
        g.AddEdge((int)pe, (int)leaf);
    }

    /* ---------- 5. lvl-0 → lvl-1 (как было) ---------- */
    {
        long long R0 = rt[0], R1 = rt[1];
        for (long long r = 0; r < R0; ++r) {
            int child = (int)(base[0] + r);
            long long block = r / shrink;
            long long start = block * shrink;
            for (int h = 0; h < shrink; ++h) {
                int parent = (int)(base[1] + (start + h) % R1);
                g.AddEdge(child, parent);
            }
        }
    }

    /* ---------- 6. остальные уровни (исправлено) ---------- */
    for (int l = 1; l + 1 < levels; ++l) {
        long long Rl   = rt[l];
        long long Rlp1 = rt[l + 1];
        long long bl   = base[l];
        long long blp1 = base[l + 1];
        for (long long r = 0; r < Rl; ++r) {
            int child  = (int)(bl + r);
            long long offset = r % Rlp1;           // позиция внутри «кольца» родителей

            for (int t = 0; t < shrink; ++t) {
                long long p_idx = (offset * shrink + t) % Rlp1;
                int parent      = (int)(blp1 + p_idx);
                g.AddEdge(child, parent);
            }
        }
    }

    /* ---------- 7. пометка compute-узлов ---------- */
    for (int l = 0; l < levels; ++l)
        for (long long r = 0; r < rt[l]; ++r)
            is_comp[(int)(base[l] + r)] = false;

    g.SetCompute(is_comp);
    return g;
}



// SMBFT(radix, depth)
//  • каждый роутер lvl ℓ ≥ 1:
//      – sibling-линки внутри группы size=radix
//      – ровно shrink=radix/2 uplink-ов mod‑shrink
Graph GenerateSMBFT(int radix, int depth) {
    if (radix < 2 || (radix & 1) || depth < 1) return Graph(0);

    int pe_count = 1;
    for (int i = 0; i < depth; ++i) pe_count *= radix;
    int levels = depth - 1;
    int shrink = radix / 2;

    // routers per level ℓ=1..levels
    std::vector<int> R(levels + 1);
    for (int l = 1; l <= levels; ++l) {
        R[l] = pe_count / static_cast<int>(std::pow(radix, l));
    }

    // base indices
    std::vector<int> base(levels + 1);
    base[1] = pe_count;
    for (int l = 2; l <= levels; ++l) {
        base[l] = base[l - 1] + R[l - 1];
    }

    int total = pe_count + std::accumulate(R.begin() + 1, R.end(), 0);
    Graph g(total);
    std::vector<bool> is_compute(total, true);

    // 1) PE → lvl1
    for (int pe = 0; pe < pe_count; ++pe) {
        int leaf = base[1] + pe / radix;
        g.AddEdge(pe, leaf);
    }

    // 2) для каждого lvl ℓ=1..levels:
    for (int l = 1; l <= levels; ++l) {
        int cnt = R[l];
        int b = base[l];

        // a) sibling-линки, если в уровне ≥ radix
        if (cnt >= radix) {
            for (int i = 0; i < cnt; ++i) {
                int self = b + i;
                is_compute[self] = false;
                int group = (i / radix) * radix;
                for (int off = 1; off < radix; ++off) {
                    int sib = b + group + (i + off) % radix;
                    if (self < sib) g.AddEdge(self, sib);
                }
            }
        }

        // b) uplink-ы lvl ℓ → ℓ+1, если не последний уровень
        if (l < levels) {
            int next_cnt = R[l + 1];
            for (int i = 0; i < cnt; ++i) {
                int child = b + i;
                int parent = base[l + 1] + (i % shrink) * (next_cnt / shrink)
                             + (i / shrink) % (next_cnt / shrink);
                // эквивалент мод‑shrink: parent_block=i/shrink, then shrink parents
                g.AddEdge(child, parent);
                is_compute[child] = false;
            }
        }
    }

    g.SetCompute(is_compute);
    return g;
}


// Hybrid-SMBFT (radix even ≥4, depth ≥2)
Graph GenerateHSMBFT(int radix, int depth) {
    if (radix < 4 || (radix & 1) || depth < 2) return Graph(0);

    /* ---------- sizes ---------- */
    int pe_cnt = 1;
    for (int i = 0; i < depth; ++i) pe_cnt *= radix;

    int levels = depth - 1;                 // router layers 1..levels
    std::vector<int> R(levels + 1);         // routers per level
    for (int l = 1; l <= levels; ++l)
        R[l] = pe_cnt / static_cast<int>(std::pow(radix, l));

    /* ---------- base indices ---------- */
    std::vector<int> base(levels + 1);
    base[1] = pe_cnt;
    for (int l = 2; l <= levels; ++l)
        base[l] = base[l - 1] + R[l - 1];

    int total = pe_cnt + std::accumulate(R.begin() + 1, R.end(), 0);
    Graph g(total);
    std::vector<bool> is_compute(total, true);

    /* ---------- PE → level-1 ---------- */
    for (int pe = 0; pe < pe_cnt; ++pe) {
        int leaf = base[1] + pe / radix;
        g.AddEdge(pe, leaf);
    }

    /* ---------- level-1: sibling + single uplink ---------- */
    {
        int cnt = R[1], b = base[1], R2 = R[2];
        for (int i = 0; i < cnt; ++i) {
            int leaf = b + i;
            is_compute[leaf] = false;

            /* sibling-links */
            int grp = (i / radix) * radix;
            for (int off = 1; off < radix; ++off) {
                int sib = b + grp + (i + off) % radix;
                if (leaf < sib) g.AddEdge(leaf, sib);
            }

            /* single uplink */
            int parent = base[2] + (i % R2);
            g.AddEdge(leaf, parent);
        }
    }

    /* ---------- levels ℓ = 2 … levels-1 : ровно shrink uplinks, равномерно разбросанных ---------- */
    const int shrink = radix / 2;
    for (int l = 2; l < levels; ++l) {
        int cnt  = R[l];         // routers on this level
        int next = R[l + 1];     // routers on level l+1
        int b    = base[l];
        int bn   = base[l + 1];

        for (int i = 0; i < cnt; ++i) {
            int self   = b + i;
            is_compute[self] = false;

            // offset внутри кольца из next родителей
            int offset = i % next;
            // каждый роутер даёт ровно shrink uplink’ов
            for (int t = 0; t < shrink; ++t) {
                int p_idx = (offset * shrink + t) % next;
                int parent = bn + p_idx;
                g.AddEdge(self, parent);
            }
        }
    }

    /* ---------- top level routers ---------- */
    for (int i = 0; i < R[levels]; ++i)
        is_compute[base[levels] + i] = false;

    g.SetCompute(is_compute);
    return g;
}


Graph GenerateFBFT(int radix_fly, int stages_fly) {
    if (radix_fly < 2 || stages_fly < 2) {
        return Graph(0);                                   // invalid request
    }

    /* ---------- basic sizes ---------- */
    int dimensions_flat = stages_fly - 1;                  // n-1
    long long pe_count = 1;
    for (int i = 0; i < stages_fly; ++i) {
        pe_count *= radix_fly;                             // k^n
    }
    long long router_count = 1;
    for (int i = 0; i < dimensions_flat; ++i) {
        router_count *= radix_fly;                         // k^(n-1)
    }

    long long total_vertices = pe_count + router_count;
    Graph graph(static_cast<int>(total_vertices));
    std::vector<bool> compute_flags(total_vertices, true); // PE = true by default

    /* ---------- strides for radix-k positional system ---------- */
    std::vector<long long> stride(dimensions_flat);
    stride[dimensions_flat - 1] = 1;
    for (int d = dimensions_flat - 2; d >= 0; --d) {
        stride[d] = stride[d + 1] * radix_fly;
    }

    /* ---------- build routers ---------- */
    long long router_base = pe_count;                      // first router id

    for (long long r_idx = 0; r_idx < router_count; ++r_idx) {
        int router_vertex = static_cast<int>(router_base + r_idx);
        compute_flags[router_vertex] = false;              // mark as router

        /* ----- attach k local PE nodes ----- */
        for (int local_id = 0; local_id < radix_fly; ++local_id) {
            int pe_vertex = static_cast<int>(r_idx * radix_fly + local_id);
            graph.AddEdge(router_vertex, pe_vertex);
        }

        /* ----- router-to-router links (flattened dimensions) ----- */
        long long rem = r_idx;
        std::vector<int> coord(dimensions_flat);
        for (int d = 0; d < dimensions_flat; ++d) {        // decode index → coordinates
            coord[d] = static_cast<int>(rem / stride[d]);
            rem = rem % stride[d];
        }

        for (int d = 0; d < dimensions_flat; ++d) {
            int original = coord[d];
            for (int v = 0; v < radix_fly; ++v) {
                if (v == original) continue;               // skip self
                long long neighbour_idx = r_idx + (v - original) * stride[d];
                int neighbour_vertex = static_cast<int>(router_base + neighbour_idx);
                if (router_vertex < neighbour_vertex) {    // add once
                    graph.AddEdge(router_vertex, neighbour_vertex);
                }
            }
        }
    }

    graph.SetCompute(compute_flags);
    return graph;
}

// Dragonfly Topology
// p   – terminals per router
// a   – routers per group
// g   – total number of groups
// (каждый роутер имеет (g–1) глобальных каналов)
//
// V = p·(a·g) + (a·g)
//
// connections:
//  1) PE → свой router
//  2) intra-group mesh (a-clique)
//  3) global: каждый router i в группе grp соединяется с router i
//     в каждой другой группе, но ребро добавляем только если self<neighbor
Graph GenerateDragonfly(int p, int a, int g) {
    if (p < 1 || a < 1 || g < 2) return Graph(0);

    int num_groups = g;
    int num_routers = a * num_groups;
    int num_terminals = p * num_routers;
    int total = num_terminals + num_routers;

    Graph graph(total);
    std::vector<bool> is_compute(total, true);

    // 1) PE → router
    for (int r = 0; r < num_routers; ++r) {
        int router_id = num_terminals + r;
        int t_base = r * p;
        for (int t = 0; t < p; ++t) {
            graph.AddEdge(t_base + t, router_id);
        }
    }

    // 2) intra-group: полный mesh в каждой группе
    for (int grp = 0; grp < num_groups; ++grp) {
        int base_r = num_terminals + grp * a;
        for (int i = 0; i < a; ++i) {
            int r1 = base_r + i;
            is_compute[r1] = false;
            for (int j = i + 1; j < a; ++j) {
                int r2 = base_r + j;
                graph.AddEdge(r1, r2);
            }
        }
    }

    // 3) global links: роутер i в grp соединяется с роутером i в других группах
    for (int grp = 0; grp < num_groups; ++grp) {
        int base_r = num_terminals + grp * a;
        for (int i = 0; i < a; ++i) {
            int self = base_r + i;
            // добавить связь только один раз: self < neighbor
            for (int other = 0; other < num_groups; ++other) {
                if (other == grp) continue;
                int neighbor = num_terminals + other * a + i;
                if (self < neighbor) {
                    graph.AddEdge(self, neighbor);
                }
            }
        }
    }

    // отметить роутеры
    for (int v = num_terminals; v < total; ++v) {
        is_compute[v] = false;
    }
    graph.SetCompute(is_compute);
    return graph;
}


Graph GenerateClos(int num_terminals, int radix) {
    if (num_terminals % radix != 0 || radix < 2) return Graph(0);

    int in = num_terminals / radix;
    int cen = in;
    int out = in;
    int baseIM = num_terminals;
    int baseCM = baseIM + in;
    int baseOM = baseCM + cen;
    int total = num_terminals + in + cen + out;

    Graph g(total);
    std::vector<bool> is_compute(total, true);

    // 1) PE → its IM
    for (int t = 0; t < num_terminals; ++t) {
        int im = baseIM + (t / radix);
        g.AddEdge(t, im);
    }

    // 2) IM ↔ CM (full bipartite)
    for (int i = 0; i < in; ++i) {
        int im = baseIM + i;
        is_compute[im] = false;
        for (int c = 0; c < cen; ++c) {
            g.AddEdge(im, baseCM + c);
        }
    }

    // 3) CM ↔ OM (full bipartite)
    for (int c = 0; c < cen; ++c) {
        int cm = baseCM + c;
        is_compute[cm] = false;
        for (int o = 0; o < out; ++o) {
            g.AddEdge(cm, baseOM + o);
        }
    }

    // 4) OM → its PEs
    for (int t = 0; t < num_terminals; ++t) {
        int om = baseOM + (t / radix);
        g.AddEdge(om, t);
    }
    for (int o = 0; o < out; ++o) {
        is_compute[baseOM + o] = false;
    }

    g.SetCompute(is_compute);
    return g;
}

Graph GenerateCnoc(int groups, int terminals_per_group, int radix) {
    if (groups < 2 || radix < 2 || (radix & 1) ||
        terminals_per_group % radix != 0) {
        return Graph(0);
    }

    int G = groups;
    int Np = terminals_per_group;
    int B = radix;
    int N = G * Np;

    int im_cnt = Np / B;      // IM и OM на группу
    int cm_cnt = B;           // CM на группу
    int om_cnt = im_cnt;

    int basePE = 0;
    int baseIM = N;
    int baseCM = baseIM + G * im_cnt;
    int baseOM = baseCM + G * im_cnt;
    int total = baseOM + G * om_cnt;

    Graph g(total);
    std::vector<bool> is_compute(total, true);

    for (int gi = 0; gi < G; ++gi) {
        int pe_group_base = gi * Np;
        int im_group_base = baseIM + gi * im_cnt;
        for (int t = 0; t < Np; ++t) {
            int pe = pe_group_base + t;
            int im_i = t / B;
            g.AddEdge(pe, im_group_base + im_i);
        }
    }

    for (int gi = 0; gi < G; ++gi) {
        int im_base = baseIM + gi * im_cnt;
        int cm_base = baseCM + gi * im_cnt;
        for (int i = 0; i < im_cnt; ++i) {
            int vIM = im_base + i;
            is_compute[vIM] = false;
            for (int j = 0; j < im_cnt; ++j) {
                g.AddEdge(vIM, cm_base + j);
            }
        }
    }

    // 3) CM ↔ OM (полносвязно внутри группы)
    for (int gi = 0; gi < G; ++gi) {
        int cm_base = baseCM;
        int om_base = baseOM;
        for (int j = 0; j < om_cnt; ++j) {
            int vCM = cm_base + j + gi * im_cnt;
            is_compute[vCM] = false;
            for (int k = 0; k < G; ++k) {
                g.AddEdge(vCM, om_base + j + k * im_cnt);
            }
        }
    }

    // 5) пометка всех роутеров
    for (int v = baseIM; v < total; ++v) {
        is_compute[v] = false;
    }
    g.SetCompute(is_compute);
    return g;
}
