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


Graph GenerateBanyan(int num_cores) {
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

                for (int p = 0; p < 2; ++p) {
                    graph.AddEdge(row * 2 + p, sw);
                }
            }
            if (stage < stages - 1) {

                int child0 = (2 * row) % switches_per_stage;
                int child1 = (2 * row + 1) % switches_per_stage;
                int sw0 = num_cores + (stage + 1) * switches_per_stage + child0;
                int sw1 = num_cores + (stage + 1) * switches_per_stage + child1;
                graph.AddEdge(sw, sw0);
                graph.AddEdge(sw, sw1);
            }
            if (stage == stages - 1) {

                for (int p = 0; p < 2; ++p) {
                    graph.AddEdge(sw, row * 2 + p);
                }
            }
        }
    }

    graph.SetCompute(compute_flags);
    return graph;
}


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


Graph GenerateOmega(int num_cores) {
    if (num_cores <= 0 || (num_cores & (num_cores - 1)) != 0) return Graph(0);
    int stages = 0;
    while ((1 << stages) < num_cores) ++stages;
    int switches_per_stage = num_cores / 2;
    int total_switches = switches_per_stage * stages;
    int total_vertices = num_cores + total_switches;
    int k = stages;

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


    int stages = 0;
    while ((1 << stages) < num_cores) ++stages;
    int switches_per_stage = num_cores / 2;
    int total_switches = switches_per_stage * stages;
    int total_vertices = num_cores + total_switches;

    Graph graph(total_vertices);
    std::vector<bool> compute(total_vertices, true);

    auto sw_idx = [&](int stage, int row) {
        return num_cores + stage * switches_per_stage + row;
    };


    for (int st = 0; st < stages; ++st)
        for (int row = 0; row < switches_per_stage; ++row)
            compute[sw_idx(st, row)] = false;


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


    int k = stages;
    int row_bits = k - 1;
    int msb_mask = 1 << (row_bits - 1);

    for (int st = 0; st < stages - 1; ++st) {
        for (int row = 0; row < switches_per_stage; ++row) {
            int cur = sw_idx(st, row);

            for (int port = 0; port < 2; ++port) {

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

    int m = std::log2(num_terminals);
    int stages = 2 * m - 1;
    int sw_per = num_terminals / 2;
    int total_sw = stages * sw_per;
    int first_sw = num_terminals;
    int total_v = num_terminals + total_sw;

    Graph g(total_v);
    std::vector<bool> is_comp(total_v, true);

    auto rot_left = [&](int w) { return ((w << 1) & (num_terminals - 1)) | (w >> (m - 1)); };
    auto rot_right = [&](int w) { return (w >> 1) | ((w & 1) << (m - 1)); };


    for (int t = 0; t < num_terminals; ++t) {
        int sw = first_sw + t / 2;
        g.AddEdge(t, sw);
    }


    for (int s = 0; s < stages - 1; ++s) {
        int base_cur = first_sw + s * sw_per;
        int base_next = first_sw + (s + 1) * sw_per;

        bool left_half = (s < m - 1);
        bool middle = (s == m - 1);

        for (int sw = 0; sw < sw_per; ++sw) {
            int self = base_cur + sw;
            is_comp[self] = false;

            for (int port = 0; port < 2; ++port) {
                int wire = sw * 2 + port;
                int wire_next = middle ? wire
                                       : left_half ? rot_left(wire)
                                                   : rot_right(wire);

                int sw_next = wire_next / 2;
                int neigh = base_next + sw_next;
                g.AddEdge(self, neigh);
            }
        }
    }


    int last_base = first_sw + (stages - 1) * sw_per;
    for (int sw = 0; sw < sw_per; ++sw) {
        int self = last_base + sw;
        is_comp[self] = false;

        for (int port = 0; port < 2; ++port) {
            int term = sw * 2 + port;
            g.AddEdge(self, term);
        }
    }


    for (int v = first_sw; v < total_v; ++v) is_comp[v] = false;
    g.SetCompute(is_comp);
    return g;
}


Graph GenerateFlattenedButterfly(int radix, int layers) {
    if (radix < 2 || layers < 2) {
        return Graph(0);
    }


    long long num_pe = 1;
    for (int i = 0; i < layers; ++i) {
        num_pe *= radix;
        if (num_pe > (1LL << 30)) {
            return Graph(0);
        }
    }
    long long num_rt = 1;
    for (int i = 0; i < layers - 1; ++i) {
        num_rt *= radix;
    }
    long long total_vertices = num_pe + num_rt;

    Graph graph(static_cast<int>(total_vertices));
    std::vector<bool> compute_flags(total_vertices, true);


    std::vector<long long> stride(layers - 1);
    stride[layers - 2] = 1;
    for (int d = layers - 3; d >= 0; --d) {
        stride[d] = stride[d + 1] * radix;
    }


    for (long long r_idx = 0; r_idx < num_rt; ++r_idx) {
        int router_v = static_cast<int>(num_pe + r_idx);
        compute_flags[router_v] = false;


        long long pe_base = r_idx * radix;
        for (int m = 0; m < radix; ++m) {
            int pe_v = static_cast<int>(pe_base + m);
            graph.AddEdge(router_v, pe_v);
        }


        long long rem = r_idx;
        std::vector<int> coord(layers - 1);
        for (int d = 0; d < layers - 1; ++d) {
            coord[d] = static_cast<int>(rem / stride[d]);
            rem %= stride[d];
        }


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


Graph GenerateBFT(int radix, int levels) {
    if (radix < 2 || (radix & 1) || levels < 2)
        return Graph(0);

    const int shrink = radix / 2;


    long long pe_cnt = 1;
    for (int i = 0; i < levels; ++i) {
        pe_cnt *= radix;
        if (pe_cnt > (1LL << 30)) return Graph(0);
    }


    std::vector<long long> rt(levels);
    rt[0] = 1;
    for (int i = 0; i < levels - 1; ++i) rt[0] *= radix;
    for (int l = 1; l < levels; ++l) rt[l] = rt[l - 1] / shrink;

    long long rt_total = 0;
    for (auto v: rt) rt_total += v;
    long long total_v = pe_cnt + rt_total;

    Graph g((int) total_v);
    std::vector<bool> is_comp(total_v, true);


    std::vector<long long> base(levels);
    base[0] = pe_cnt;
    for (int l = 1; l < levels; ++l) base[l] = base[l - 1] + rt[l - 1];


    for (long long pe = 0; pe < pe_cnt; ++pe) {
        long long leaf = base[0] + pe / radix;
        g.AddEdge((int) pe, (int) leaf);
    }


    {
        long long R0 = rt[0], R1 = rt[1];
        for (long long r = 0; r < R0; ++r) {
            int child = (int) (base[0] + r);
            long long block = r / shrink;
            long long start = block * shrink;
            for (int h = 0; h < shrink; ++h) {
                int parent = (int) (base[1] + (start + h) % R1);
                g.AddEdge(child, parent);
            }
        }
    }


    for (int l = 1; l + 1 < levels; ++l) {
        long long Rl = rt[l];
        long long Rlp1 = rt[l + 1];
        long long bl = base[l];
        long long blp1 = base[l + 1];
        for (long long r = 0; r < Rl; ++r) {
            int child = (int) (bl + r);
            long long offset = r % Rlp1;

            for (int t = 0; t < shrink; ++t) {
                long long p_idx = (offset * shrink + t) % Rlp1;
                int parent = (int) (blp1 + p_idx);
                g.AddEdge(child, parent);
            }
        }
    }


    for (int l = 0; l < levels; ++l)
        for (long long r = 0; r < rt[l]; ++r)
            is_comp[(int) (base[l] + r)] = false;

    g.SetCompute(is_comp);
    return g;
}


Graph GenerateSMBFT(int radix, int depth) {
    if (radix < 2 || (radix & 1) || depth < 1) return Graph(0);

    int pe_count = 1;
    for (int i = 0; i < depth; ++i) pe_count *= radix;
    int levels = depth - 1;
    int shrink = radix / 2;


    std::vector<int> R(levels + 1);
    for (int l = 1; l <= levels; ++l) {
        R[l] = pe_count / static_cast<int>(std::pow(radix, l));
    }


    std::vector<int> base(levels + 1);
    base[1] = pe_count;
    for (int l = 2; l <= levels; ++l) {
        base[l] = base[l - 1] + R[l - 1];
    }

    int total = pe_count + std::accumulate(R.begin() + 1, R.end(), 0);
    Graph g(total);
    std::vector<bool> is_compute(total, true);


    for (int pe = 0; pe < pe_count; ++pe) {
        int leaf = base[1] + pe / radix;
        g.AddEdge(pe, leaf);
    }


    for (int l = 1; l <= levels; ++l) {
        int cnt = R[l];
        int b = base[l];


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


        if (l < levels) {
            int next_cnt = R[l + 1];
            for (int i = 0; i < cnt; ++i) {
                int child = b + i;
                int parent = base[l + 1] + (i % shrink) * (next_cnt / shrink)
                             + (i / shrink) % (next_cnt / shrink);

                g.AddEdge(child, parent);
                is_compute[child] = false;
            }
        }
    }

    g.SetCompute(is_compute);
    return g;
}


Graph GenerateHSMBFT(int radix, int depth) {
    if (radix < 4 || (radix & 1) || depth < 2) return Graph(0);


    int pe_cnt = 1;
    for (int i = 0; i < depth; ++i) pe_cnt *= radix;

    int levels = depth - 1;
    std::vector<int> R(levels + 1);
    for (int l = 1; l <= levels; ++l)
        R[l] = pe_cnt / static_cast<int>(std::pow(radix, l));


    std::vector<int> base(levels + 1);
    base[1] = pe_cnt;
    for (int l = 2; l <= levels; ++l)
        base[l] = base[l - 1] + R[l - 1];

    int total = pe_cnt + std::accumulate(R.begin() + 1, R.end(), 0);
    Graph g(total);
    std::vector<bool> is_compute(total, true);


    for (int pe = 0; pe < pe_cnt; ++pe) {
        int leaf = base[1] + pe / radix;
        g.AddEdge(pe, leaf);
    }


    {
        int cnt = R[1], b = base[1], R2 = R[2];
        for (int i = 0; i < cnt; ++i) {
            int leaf = b + i;
            is_compute[leaf] = false;


            int grp = (i / radix) * radix;
            for (int off = 1; off < radix; ++off) {
                int sib = b + grp + (i + off) % radix;
                if (leaf < sib) g.AddEdge(leaf, sib);
            }


            int parent = base[2] + (i % R2);
            g.AddEdge(leaf, parent);
        }
    }


    const int shrink = radix / 2;
    for (int l = 2; l < levels; ++l) {
        int cnt = R[l];
        int next = R[l + 1];
        int b = base[l];
        int bn = base[l + 1];

        for (int i = 0; i < cnt; ++i) {
            int self = b + i;
            is_compute[self] = false;


            int offset = i % next;

            for (int t = 0; t < shrink; ++t) {
                int p_idx = (offset * shrink + t) % next;
                int parent = bn + p_idx;
                g.AddEdge(self, parent);
            }
        }
    }


    for (int i = 0; i < R[levels]; ++i)
        is_compute[base[levels] + i] = false;

    g.SetCompute(is_compute);
    return g;
}


Graph GenerateFBFT(int radix_fly, int stages_fly) {
    if (radix_fly < 2 || stages_fly < 2) {
        return Graph(0);
    }


    int dimensions_flat = stages_fly - 1;
    long long pe_count = 1;
    for (int i = 0; i < stages_fly; ++i) {
        pe_count *= radix_fly;
    }
    long long router_count = 1;
    for (int i = 0; i < dimensions_flat; ++i) {
        router_count *= radix_fly;
    }

    long long total_vertices = pe_count + router_count;
    Graph graph(static_cast<int>(total_vertices));
    std::vector<bool> compute_flags(total_vertices, true);


    std::vector<long long> stride(dimensions_flat);
    stride[dimensions_flat - 1] = 1;
    for (int d = dimensions_flat - 2; d >= 0; --d) {
        stride[d] = stride[d + 1] * radix_fly;
    }


    long long router_base = pe_count;

    for (long long r_idx = 0; r_idx < router_count; ++r_idx) {
        int router_vertex = static_cast<int>(router_base + r_idx);
        compute_flags[router_vertex] = false;


        for (int local_id = 0; local_id < radix_fly; ++local_id) {
            int pe_vertex = static_cast<int>(r_idx * radix_fly + local_id);
            graph.AddEdge(router_vertex, pe_vertex);
        }


        long long rem = r_idx;
        std::vector<int> coord(dimensions_flat);
        for (int d = 0; d < dimensions_flat; ++d) {
            coord[d] = static_cast<int>(rem / stride[d]);
            rem = rem % stride[d];
        }

        for (int d = 0; d < dimensions_flat; ++d) {
            int original = coord[d];
            for (int v = 0; v < radix_fly; ++v) {
                if (v == original) continue;
                long long neighbour_idx = r_idx + (v - original) * stride[d];
                int neighbour_vertex = static_cast<int>(router_base + neighbour_idx);
                if (router_vertex < neighbour_vertex) {
                    graph.AddEdge(router_vertex, neighbour_vertex);
                }
            }
        }
    }

    graph.SetCompute(compute_flags);
    return graph;
}


Graph GenerateDragonfly(int p, int a, int g) {
    if (p < 1 || a < 1 || g < 2) return Graph(0);

    int num_groups = g;
    int num_routers = a * num_groups;
    int num_terminals = p * num_routers;
    int total = num_terminals + num_routers;

    Graph graph(total);
    std::vector<bool> is_compute(total, true);


    for (int r = 0; r < num_routers; ++r) {
        int router_id = num_terminals + r;
        int t_base = r * p;
        for (int t = 0; t < p; ++t) {
            graph.AddEdge(t_base + t, router_id);
        }
    }


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


    for (int grp = 0; grp < num_groups; ++grp) {
        int base_r = num_terminals + grp * a;
        for (int i = 0; i < a; ++i) {
            int self = base_r + i;

            for (int other = 0; other < num_groups; ++other) {
                if (other == grp) continue;
                int neighbor = num_terminals + other * a + i;
                if (self < neighbor) {
                    graph.AddEdge(self, neighbor);
                }
            }
        }
    }


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


    for (int t = 0; t < num_terminals; ++t) {
        int im = baseIM + (t / radix);
        g.AddEdge(t, im);
    }


    for (int i = 0; i < in; ++i) {
        int im = baseIM + i;
        is_compute[im] = false;
        for (int c = 0; c < cen; ++c) {
            g.AddEdge(im, baseCM + c);
        }
    }


    for (int c = 0; c < cen; ++c) {
        int cm = baseCM + c;
        is_compute[cm] = false;
        for (int o = 0; o < out; ++o) {
            g.AddEdge(cm, baseOM + o);
        }
    }


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

    int im_cnt = Np / B;
    int cm_cnt = B;
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


    for (int v = baseIM; v < total; ++v) {
        is_compute[v] = false;
    }
    g.SetCompute(is_compute);
    return g;
}
