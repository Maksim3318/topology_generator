#pragma once

#include "Graph.hpp"

Graph GenerateFT(int radix, int levels) {
    // radix: branching factor, levels: number of tree levels (≥1)
    if (radix < 1 || levels < 1) return Graph(0);
    // 1) compute number of nodes and start index per level
    std::vector<int> count(levels), start(levels);
    count[0] = 1;
    for (int l = 1; l < levels; ++l)
        count[l] = count[l-1] * radix;
    start[0] = 0;
    for (int l = 1; l < levels; ++l)
        start[l] = start[l-1] + count[l-1];
    int N = start[levels-1] + count[levels-1];

    // 2) build graph
    Graph g(N);
    for (int l = 0; l < levels - 1; ++l) {
        for (int i = 0; i < count[l]; ++i) {
            int parent = start[l] + i;
            int child_base = start[l+1] + i * radix;
            for (int j = 0; j < radix; ++j) {
                g.AddEdge(parent, child_base + j);
            }
        }
    }

    // 3) mark compute nodes: only leaves (last level)
    std::vector<bool> compute(N, false);
    int leaf_start = start[levels-1];
    int leaf_count = count[levels-1];
    for (int i = 0; i < leaf_count; ++i) {
        compute[leaf_start + i] = true;
    }
    g.SetCompute(compute);

    return g;
}

Graph GenerateXTree(int radix, int levels) {
    // radix: branching factor, levels: tree height (≥1)
    if (radix < 1 || levels < 1) return Graph(0);

    // 1) compute counts and start indices per level
    std::vector<int> count(levels), start(levels);
    count[0] = 1;
    for (int l = 1; l < levels; ++l)
        count[l] = count[l-1] * radix;
    start[0] = 0;
    for (int l = 1; l < levels; ++l)
        start[l] = start[l-1] + count[l-1];
    int N = start[levels-1] + count[levels-1];

    // 2) build X-Tree
    Graph g(N);

    // 2.1) parent→children links
    for (int l = 0; l + 1 < levels; ++l) {
        for (int i = 0; i < count[l]; ++i) {
            int parent = start[l] + i;
            int child_base = start[l+1] + i * radix;
            for (int j = 0; j < radix; ++j) {
                g.AddEdge(parent, child_base + j);
            }
        }
    }

    // 2.2) chain links on each level
    for (int l = 0; l < levels; ++l) {
        for (int i = 0; i + 1 < count[l]; ++i) {
            int u = start[l] + i;
            int v = start[l] + i + 1;
            g.AddEdge(u, v);
        }
    }

    // compute flags are true by default in the constructor

    return g;
}


Graph GenerateGFT(int h, int m, int w) {
    // h: число уровней (≥0)
    // m: фактор репликации подграфа (≥1)
    // w: fan-out новых корней на каждом уровне (≥1)
    if (h < 0 || m < 1 || w < 1) return Graph(0);

    // 1) уровень 0: один лист (вычислительный узел)
    Graph g_current(1);
    std::vector<int> leaf_ids = {0};
    std::vector<int> root_ids = {0};

    // 2) наращиваем уровни 1…h
    for (int level = 1; level <= h; ++level) {
        int N_sub   = g_current.NumVertices();
        int R_sub   = (int)root_ids.size();
        int N_next  = N_sub * m + R_sub * w;
        Graph g_next(N_next);

        // 2.1 Репликация подграфа m раз
        for (int rep = 0; rep < m; ++rep) {
            int off = rep * N_sub;
            for (int u = 0; u < N_sub; ++u)
                for (int v : g_current.Neighbors(u))
                    if (u < v)
                        g_next.AddEdge(off + u, off + v);
        }

        // 2.2 Создание новых корней
        std::vector<int> new_roots;
        new_roots.reserve(R_sub * w);
        int base = N_sub * m;
        for (int i = 0; i < R_sub; ++i)
            for (int k = 0; k < w; ++k)
                new_roots.push_back(base + i * w + k);

        // 2.3 Связывание старых корней с новыми
        for (int rep = 0; rep < m; ++rep) {
            int off = rep * N_sub;
            for (int i = 0; i < R_sub; ++i) {
                int parent = off + root_ids[i];
                for (int k = 0; k < w; ++k) {
                    int child = base + i * w + k;
                    if (parent < child) g_next.AddEdge(parent, child);
                    else               g_next.AddEdge(child, parent);
                }
            }
        }

        // 2.4 Обновляем листовые индексы (только реплики старых листьев)
        std::vector<int> new_leaves;
        new_leaves.reserve(leaf_ids.size() * m);
        for (int rep = 0; rep < m; ++rep) {
            int off = rep * N_sub;
            for (int id : leaf_ids)
                new_leaves.push_back(off + id);
        }

        // 2.5 Переходим на следующий уровень
        g_current = std::move(g_next);
        root_ids  = std::move(new_roots);
        leaf_ids  = std::move(new_leaves);
    }

    // 3) Помечаем вычислительные узлы — только листья
    std::vector<bool> compute(g_current.NumVertices(), false);
    for (int id : leaf_ids) compute[id] = true;
    g_current.SetCompute(compute);

    return g_current;
}

Graph GeneratePyramid(int radix, int levels) {
    // radix ≥ 1, levels ≥ 1
    if (radix < 1 || levels < 1) return Graph(0);

    // 1) side[l] = grid side length at level l (top l=0 → 1, bottom l=levels-1 → radix^(levels-1))
    std::vector<int> side(levels);
    side[0] = 1;
    for (int l = 1; l < levels; ++l)
        side[l] = side[l-1] * radix;

    // 2) compute starting index of each level
    std::vector<int> start(levels);
    start[0] = 0;
    for (int l = 1; l < levels; ++l)
        start[l] = start[l-1] + side[l-1] * side[l-1];

    // 3) total vertices
    int N = start[levels-1] + side[levels-1] * side[levels-1];
    Graph g(N);

    // 4) intra‐level grid edges (4‐neighbor mesh, no wrap)
    for (int l = 0; l < levels; ++l) {
        int s = side[l], base = start[l];
        for (int r = 0; r < s; ++r) {
            for (int c = 0; c < s; ++c) {
                int u = base + r*s + c;
                if (c + 1 < s) g.AddEdge(u, base + r*s + (c+1));
                if (r + 1 < s) g.AddEdge(u, base + (r+1)*s + c);
            }
        }
    }

    // 5) inter‐level “pyramid” edges (each parent → radix×radix children)
    for (int l = 0; l + 1 < levels; ++l) {
        int s   = side[l],   sc = side[l+1];
        int baseP = start[l], baseC = start[l+1];
        for (int r = 0; r < s; ++r) {
            for (int c = 0; c < s; ++c) {
                int parent = baseP + r*s + c;
                int br = r * radix, bc = c * radix;
                for (int dr = 0; dr < radix; ++dr) {
                    for (int dc = 0; dc < radix; ++dc) {
                        int child = baseC + (br+dr)*sc + (bc+dc);
                        g.AddEdge(parent, child);
                    }
                }
            }
        }
    }

    // 6) mark compute only on bottom level
    std::vector<bool> compute(N, false);
    int leafStart = start[levels-1];
    for (int i = leafStart; i < N; ++i)
        compute[i] = true;
    g.SetCompute(compute);

    return g;
}

Graph GenerateHierarchicalClique(int radix, int levels) {
    // radix ≥ 1, levels ≥ 1
    if (radix < 1 || levels < 1) return Graph(0);

    // 1) compute counts and start indices
    std::vector<int> count(levels), start(levels);
    count[0] = radix;
    for (int l = 1; l < levels; ++l)
        count[l] = count[l-1] * radix;
    start[0] = 0;
    for (int l = 1; l < levels; ++l)
        start[l] = start[l-1] + count[l-1];

    int N = start[levels-1] + count[levels-1];
    Graph g(N);

    // 2) level 0 clique
    for (int a = 0; a < count[0]; ++a)
        for (int b = a+1; b < count[0]; ++b)
            g.AddEdge(start[0]+a, start[0]+b);

    // 3) for each subsequent level: parent→children + children clique
    for (int l = 1; l < levels; ++l) {
        int r = radix;
        int base = start[l];
        int prevBase = start[l-1];
        for (int p = 0; p < count[l-1]; ++p) {
            int parent = prevBase + p;
            int childBase = base + p * r;
            // parent→children
            for (int j = 0; j < r; ++j)
                g.AddEdge(parent, childBase + j);
            // children clique
            for (int a = 0; a < r; ++a)
                for (int b = a+1; b < r; ++b)
                    g.AddEdge(childBase + a, childBase + b);
        }
    }

    // 4) mark compute only on leaves (last level)
    std::vector<bool> compute(N, false);
    int leafStart = start[levels-1];
    for (int i = 0; i < count[levels-1]; ++i)
        compute[leafStart + i] = true;
    g.SetCompute(compute);

    return g;
}

Graph GenerateCubeTreeHybrid(int levels) {
    // levels: число уровней иерархии ≥2
    if (levels < 2) return Graph(0);

    // 1) вычисляем число листов: 2^(levels+2)
    int numLeaves = 1 << (levels + 2);

    // 2) считаем количество узлов на каждом уровне l = 0..levels-1
    std::vector<int> count(levels), start(levels);
    count[0] = numLeaves;
    count[1] = count[0] >> 2;
    for (int l = 2; l < levels; ++l)
        count[l] = count[l-1] >> 1;

    // 3) вычисляем индексы начала каждого уровня
    start[0] = 0;
    for (int l = 1; l < levels; ++l)
        start[l] = start[l-1] + count[l-1];

    int N = start.back() + count.back();
    Graph g(N);

    // 4) помечаем вычислительные узлы: только на уровне 0 (листья)
    std::vector<bool> compute(N, false);
    for (int i = 0; i < count[0]; ++i)
        compute[i] = true;
    g.SetCompute(compute);

    // 5) «деревянные» связи: каждый узел уровня l соединён с D детьми уровня l−1
    for (int l = 1; l < levels; ++l) {
        int D = (l == 1 ? 4 : 2);
        for (int i = 0; i < count[l]; ++i) {
            int parent    = start[l]   + i;
            int childBase = start[l-1] + i * D;
            for (int j = 0; j < D; ++j)
                g.AddEdge(parent, childBase + j);
        }
    }

    // 6) «кубовые» клики на уровнях ≥2: каждая группа по 4 узла образует K4
    for (int l = 2; l < levels; ++l) {
        int groups = count[l] / 4;
        for (int gi = 0; gi < groups; ++gi) {
            int baseIdx = start[l] + gi * 4;
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b)
                    g.AddEdge(baseIdx + a, baseIdx + b);
        }
    }

    // 7) подключаем все узлы на самом верхнем уровне в один кластер (clique),
    //    чтобы гарантировать связность
    int topBase = start[levels-1];
    int topCount = count[levels-1];
    for (int a = 0; a < topCount; ++a)
        for (int b = a + 1; b < topCount; ++b)
            g.AddEdge(topBase + a, topBase + b);

    return g;
}
