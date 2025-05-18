#pragma once

#include "Graph.hpp"

Graph GenerateFT(int radix, int levels) {

    if (radix < 1 || levels < 1) return Graph(0);

    std::vector<int> count(levels), start(levels);
    count[0] = 1;
    for (int l = 1; l < levels; ++l)
        count[l] = count[l - 1] * radix;
    start[0] = 0;
    for (int l = 1; l < levels; ++l)
        start[l] = start[l - 1] + count[l - 1];
    int N = start[levels - 1] + count[levels - 1];


    Graph g(N);
    for (int l = 0; l < levels - 1; ++l) {
        for (int i = 0; i < count[l]; ++i) {
            int parent = start[l] + i;
            int child_base = start[l + 1] + i * radix;
            for (int j = 0; j < radix; ++j) {
                g.AddEdge(parent, child_base + j);
            }
        }
    }


    std::vector<bool> compute(N, false);
    int leaf_start = start[levels - 1];
    int leaf_count = count[levels - 1];
    for (int i = 0; i < leaf_count; ++i) {
        compute[leaf_start + i] = true;
    }
    g.SetCompute(compute);

    return g;
}

Graph GenerateXTree(int radix, int levels) {

    if (radix < 1 || levels < 1) return Graph(0);


    std::vector<int> count(levels), start(levels);
    count[0] = 1;
    for (int l = 1; l < levels; ++l)
        count[l] = count[l - 1] * radix;
    start[0] = 0;
    for (int l = 1; l < levels; ++l)
        start[l] = start[l - 1] + count[l - 1];
    int N = start[levels - 1] + count[levels - 1];


    Graph g(N);


    for (int l = 0; l + 1 < levels; ++l) {
        for (int i = 0; i < count[l]; ++i) {
            int parent = start[l] + i;
            int child_base = start[l + 1] + i * radix;
            for (int j = 0; j < radix; ++j) {
                g.AddEdge(parent, child_base + j);
            }
        }
    }


    for (int l = 0; l < levels; ++l) {
        for (int i = 0; i + 1 < count[l]; ++i) {
            int u = start[l] + i;
            int v = start[l] + i + 1;
            g.AddEdge(u, v);
        }
    }


    return g;
}


Graph GenerateGFT(int h, int m, int w) {


    if (h < 0 || m < 1 || w < 1) return Graph(0);


    Graph g_current(1);
    std::vector<int> leaf_ids = {0};
    std::vector<int> root_ids = {0};


    for (int level = 1; level <= h; ++level) {
        int N_sub = g_current.NumVertices();
        int R_sub = (int) root_ids.size();
        int N_next = N_sub * m + R_sub * w;
        Graph g_next(N_next);


        for (int rep = 0; rep < m; ++rep) {
            int off = rep * N_sub;
            for (int u = 0; u < N_sub; ++u)
                for (int v: g_current.Neighbors(u))
                    if (u < v)
                        g_next.AddEdge(off + u, off + v);
        }


        std::vector<int> new_roots;
        new_roots.reserve(R_sub * w);
        int base = N_sub * m;
        for (int i = 0; i < R_sub; ++i)
            for (int k = 0; k < w; ++k)
                new_roots.push_back(base + i * w + k);


        for (int rep = 0; rep < m; ++rep) {
            int off = rep * N_sub;
            for (int i = 0; i < R_sub; ++i) {
                int parent = off + root_ids[i];
                for (int k = 0; k < w; ++k) {
                    int child = base + i * w + k;
                    if (parent < child) g_next.AddEdge(parent, child);
                    else g_next.AddEdge(child, parent);
                }
            }
        }


        std::vector<int> new_leaves;
        new_leaves.reserve(leaf_ids.size() * m);
        for (int rep = 0; rep < m; ++rep) {
            int off = rep * N_sub;
            for (int id: leaf_ids)
                new_leaves.push_back(off + id);
        }


        g_current = std::move(g_next);
        root_ids = std::move(new_roots);
        leaf_ids = std::move(new_leaves);
    }


    std::vector<bool> compute(g_current.NumVertices(), false);
    for (int id: leaf_ids) compute[id] = true;
    g_current.SetCompute(compute);

    return g_current;
}

Graph GeneratePyramid(int radix, int levels) {

    if (radix < 1 || levels < 1) return Graph(0);


    std::vector<int> side(levels);
    side[0] = 1;
    for (int l = 1; l < levels; ++l)
        side[l] = side[l - 1] * radix;


    std::vector<int> start(levels);
    start[0] = 0;
    for (int l = 1; l < levels; ++l)
        start[l] = start[l - 1] + side[l - 1] * side[l - 1];


    int N = start[levels - 1] + side[levels - 1] * side[levels - 1];
    Graph g(N);


    for (int l = 0; l < levels; ++l) {
        int s = side[l], base = start[l];
        for (int r = 0; r < s; ++r) {
            for (int c = 0; c < s; ++c) {
                int u = base + r * s + c;
                if (c + 1 < s) g.AddEdge(u, base + r * s + (c + 1));
                if (r + 1 < s) g.AddEdge(u, base + (r + 1) * s + c);
            }
        }
    }


    for (int l = 0; l + 1 < levels; ++l) {
        int s = side[l], sc = side[l + 1];
        int baseP = start[l], baseC = start[l + 1];
        for (int r = 0; r < s; ++r) {
            for (int c = 0; c < s; ++c) {
                int parent = baseP + r * s + c;
                int br = r * radix, bc = c * radix;
                for (int dr = 0; dr < radix; ++dr) {
                    for (int dc = 0; dc < radix; ++dc) {
                        int child = baseC + (br + dr) * sc + (bc + dc);
                        g.AddEdge(parent, child);
                    }
                }
            }
        }
    }


    std::vector<bool> compute(N, false);
    int leafStart = start[levels - 1];
    for (int i = leafStart; i < N; ++i)
        compute[i] = true;
    g.SetCompute(compute);

    return g;
}

Graph GenerateHierarchicalClique(int radix, int levels) {

    if (radix < 1 || levels < 1) return Graph(0);


    std::vector<int> count(levels), start(levels);
    count[0] = radix;
    for (int l = 1; l < levels; ++l)
        count[l] = count[l - 1] * radix;
    start[0] = 0;
    for (int l = 1; l < levels; ++l)
        start[l] = start[l - 1] + count[l - 1];

    int N = start[levels - 1] + count[levels - 1];
    Graph g(N);


    for (int a = 0; a < count[0]; ++a)
        for (int b = a + 1; b < count[0]; ++b)
            g.AddEdge(start[0] + a, start[0] + b);


    for (int l = 1; l < levels; ++l) {
        int r = radix;
        int base = start[l];
        int prevBase = start[l - 1];
        for (int p = 0; p < count[l - 1]; ++p) {
            int parent = prevBase + p;
            int childBase = base + p * r;

            for (int j = 0; j < r; ++j)
                g.AddEdge(parent, childBase + j);

            for (int a = 0; a < r; ++a)
                for (int b = a + 1; b < r; ++b)
                    g.AddEdge(childBase + a, childBase + b);
        }
    }


    std::vector<bool> compute(N, false);
    int leafStart = start[levels - 1];
    for (int i = 0; i < count[levels - 1]; ++i)
        compute[leafStart + i] = true;
    g.SetCompute(compute);

    return g;
}

Graph GenerateCubeTreeHybrid(int levels) {

    if (levels < 2) return Graph(0);


    int numLeaves = 1 << (levels + 2);


    std::vector<int> count(levels), start(levels);
    count[0] = numLeaves;
    count[1] = count[0] >> 2;
    for (int l = 2; l < levels; ++l)
        count[l] = count[l - 1] >> 1;


    start[0] = 0;
    for (int l = 1; l < levels; ++l)
        start[l] = start[l - 1] + count[l - 1];

    int N = start.back() + count.back();
    Graph g(N);


    std::vector<bool> compute(N, false);
    for (int i = 0; i < count[0]; ++i)
        compute[i] = true;
    g.SetCompute(compute);


    for (int l = 1; l < levels; ++l) {
        int D = (l == 1 ? 4 : 2);
        for (int i = 0; i < count[l]; ++i) {
            int parent = start[l] + i;
            int childBase = start[l - 1] + i * D;
            for (int j = 0; j < D; ++j)
                g.AddEdge(parent, childBase + j);
        }
    }


    for (int l = 2; l < levels; ++l) {
        int groups = count[l] / 4;
        for (int gi = 0; gi < groups; ++gi) {
            int baseIdx = start[l] + gi * 4;
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b)
                    g.AddEdge(baseIdx + a, baseIdx + b);
        }
    }


    int topBase = start[levels - 1];
    int topCount = count[levels - 1];
    for (int a = 0; a < topCount; ++a)
        for (int b = a + 1; b < topCount; ++b)
            g.AddEdge(topBase + a, topBase + b);

    return g;
}
