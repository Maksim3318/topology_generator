#pragma once

#include "Graph.hpp"

Graph GenerateRicobit(int K) {
    // K: number of concentric rings (≥1)
    if (K < 1) return Graph(0);

    // 1) compute ring sizes 2^L and start indices
    std::vector<int> ringSize(K + 1), start(K + 1);
    for (int L = 1; L <= K; ++L)
        ringSize[L] = 1 << L;
    start[1] = 0;
    for (int L = 2; L <= K; ++L)
        start[L] = start[L - 1] + ringSize[L - 1];

    int N = start[K] + ringSize[K];
    Graph g(N);
    // all nodes are compute by default

    // 2) intra-ring edges (each ring is a cycle)
    for (int L = 1; L <= K; ++L) {
        int sz = ringSize[L], base = start[L];
        for (int i = 0; i < sz; ++i) {
            int j = (i + 1) % sz;
            if (i < j)
                g.AddEdge(base + i, base + j);
            else
                g.AddEdge(base + j, base + i);
        }
    }

    // 3) inter-ring “radial” edges
    //    node n in ring L connects to 2n,2n+1 in ring L+1
    for (int L = 1; L < K; ++L) {
        int base = start[L], nextBase = start[L + 1];
        int sz = ringSize[L];
        for (int n = 0; n < sz; ++n) {
            int u = base + n;
            int v0 = nextBase + 2 * n;
            int v1 = nextBase + 2 * n + 1;
            g.AddEdge(u, v0);
            g.AddEdge(u, v1);
        }
    }

    return g;
}


inline Graph GenerateTreeMesh(int mesh_k) {
    // mesh_k must be a power of two (≥2)
    if (mesh_k < 2 || (mesh_k & (mesh_k - 1)) != 0) return Graph(0);

    // 1) tree depth = log2(mesh_k)
    int levels = 0;
    while ((1 << levels) < mesh_k) ++levels;

    int meshN = mesh_k * mesh_k;

    // 2) count internal tree nodes: sum_{t=0..levels-1} 4^t
    std::vector<int> treeStart(levels);
    int internalN = 0;
    for (int t = 0; t < levels; ++t) {
        treeStart[t] = meshN + internalN;
        internalN += (1 << (2 * t));
    }
    int N = meshN + internalN;
    Graph g(N);

    // 3) mesh edges among level-0 nodes
    for (int r = 0; r < mesh_k; ++r) {
        for (int c = 0; c < mesh_k; ++c) {
            int u = r * mesh_k + c;
            if (c + 1 < mesh_k) g.AddEdge(u, u + 1);
            if (r + 1 < mesh_k) g.AddEdge(u, u + mesh_k);
        }
    }

    // 4) tree edges
    for (int t = 0; t < levels; ++t) {
        int blocks = 1 << t;
        int childBlocks = blocks << 1;
        int base = treeStart[t];
        for (int rr = 0; rr < blocks; ++rr) {
            for (int cc = 0; cc < blocks; ++cc) {
                int parent = base + rr * blocks + cc;
                if (t + 1 == levels) {
                    // children are mesh leaves
                    int size = mesh_k >> t;
                    int r0 = rr * size, c0 = cc * size;
                    for (int dr = 0; dr < 2; ++dr) {
                        for (int dc = 0; dc < 2; ++dc) {
                            int leaf = (r0 + dr * (size - 1)) * mesh_k + (c0 + dc * (size - 1));
                            g.AddEdge(parent, leaf);
                        }
                    }
                } else {
                    int childBase = treeStart[t + 1];
                    for (int dr = 0; dr < 2; ++dr) {
                        for (int dc = 0; dc < 2; ++dc) {
                            int child = childBase + (rr * 2 + dr) * childBlocks + (cc * 2 + dc);
                            g.AddEdge(parent, child);
                        }
                    }
                }
            }
        }
    }

    // 5) mark compute only on mesh nodes (level 0)
    std::vector<bool> compute(N, false);
    for (int i = 0; i < meshN; ++i)
        compute[i] = true;
    g.SetCompute(compute);

    return g;
}


Graph GenerateMeshOfTrees(int mesh_k) {
    // mesh_k: side length of mesh, must be power of two ≥2
    if (mesh_k < 2 || (mesh_k & (mesh_k - 1)) != 0) return Graph(0);

    int h = 0;
    while ((1 << h) < mesh_k) ++h;

    int leafCount = mesh_k * mesh_k;
    int internalPerTree = mesh_k - 1;
    int totalInternal = 2 * mesh_k * internalPerTree;
    int N = leafCount + totalInternal;
    Graph g(N);

    int baseRow = leafCount;
    int baseCol = leafCount + mesh_k * internalPerTree;
    int maxHeap = 2 * mesh_k - 1;

    // row trees
    for (int r = 0; r < mesh_k; ++r) {
        int offset = baseRow + r * internalPerTree;
        for (int i = 1; i <= internalPerTree; ++i) {
            int parent = offset + (i - 1);
            int left = 2 * i, right = 2 * i + 1;
            if (left <= maxHeap) {
                int child = (left <= internalPerTree
                             ? offset + (left - 1)
                             : r * mesh_k + (left - mesh_k));
                g.AddEdge(parent, child);
            }
            if (right <= maxHeap) {
                int child = (right <= internalPerTree
                             ? offset + (right - 1)
                             : r * mesh_k + (right - mesh_k));
                g.AddEdge(parent, child);
            }
        }
    }

    // column trees
    for (int c = 0; c < mesh_k; ++c) {
        int offset = baseCol + c * internalPerTree;
        for (int i = 1; i <= internalPerTree; ++i) {
            int parent = offset + (i - 1);
            int left = 2 * i, right = 2 * i + 1;
            if (left <= maxHeap) {
                int child = (left <= internalPerTree
                             ? offset + (left - 1)
                             : (left - mesh_k) * mesh_k + c);
                g.AddEdge(parent, child);
            }
            if (right <= maxHeap) {
                int child = (right <= internalPerTree
                             ? offset + (right - 1)
                             : (right - mesh_k) * mesh_k + c);
                g.AddEdge(parent, child);
            }
        }
    }

    // mark compute only on mesh leaves
    std::vector<bool> compute(N, false);
    for (int i = 0; i < leafCount; ++i)
        compute[i] = true;
    g.SetCompute(compute);

    return g;
}

inline Graph GenerateFatHTree(int n) {
    // n: log2 of grid side length (≥1)
    if (n < 1) return Graph(0);

    int Nside = 1 << n;
    int coreCount = Nside * Nside;

    // 1) compute counts and start indices for red ranks 0..n
    std::vector<int> count(n+1), start_red(n+1);
    for (int lvl = 0; lvl <= n; ++lvl) {
        int dim = Nside >> lvl;
        count[lvl] = dim * dim;
        start_red[lvl] = (lvl == 0 ? 0 : start_red[lvl-1] + count[lvl-1]);
    }
    // 2) compute start indices for black ranks 1..n
    std::vector<int> start_black(n+1);
    start_black[1] = start_red[n] + count[n];
    for (int lvl = 2; lvl <= n; ++lvl) {
        start_black[lvl] = start_black[lvl-1] + count[lvl-1];
    }
    int N = start_black[n] + count[n];
    Graph g(N);

    // 3) mark only rank-0 (cores) as compute nodes
    std::vector<bool> compute(N, false);
    for (int i = 0; i < coreCount; ++i) compute[i] = true;
    g.SetCompute(compute);

    // 4) red tree: for each rank lvl = 0..n-1, connect to lvl+1
    for (int lvl = 0; lvl < n; ++lvl) {
        int dim    = Nside >> lvl;
        int nextDim= Nside >> (lvl+1);
        int base   = start_red[lvl];
        int baseP  = start_red[lvl+1];
        for (int y = 0; y < dim; ++y) {
            for (int x = 0; x < dim; ++x) {
                int u = base + y*dim + x;
                int vx = x >> 1, vy = y >> 1;
                int v = baseP + vy*nextDim + vx;
                if (u < v) g.AddEdge(u, v);
                else       g.AddEdge(v, u);
            }
        }
    }

    // 5) black tree: level 0→1 from cores
    for (int y = 0; y < Nside; ++y) {
        for (int x = 0; x < Nside; ++x) {
            int X0 = (x + Nside - 1) & (Nside - 1);
            int Y0 = (y + Nside - 1) & (Nside - 1);
            int u  = y * Nside + x;
            int base1 = start_black[1];
            int dim1  = Nside >> 1;
            int vx = X0 >> 1, vy = Y0 >> 1;
            int v  = base1 + vy*dim1 + vx;
            if (u < v) g.AddEdge(u, v);
            else       g.AddEdge(v, u);
        }
    }
    // black tree: levels 1..n-1
    for (int lvl = 1; lvl < n; ++lvl) {
        int dim     = Nside >> lvl;
        int nextDim = Nside >> (lvl+1);
        int base    = start_black[lvl];
        int baseP   = start_black[lvl+1];
        for (int y = 0; y < dim; ++y) {
            for (int x = 0; x < dim; ++x) {
                int u = base + y*dim + x;
                int vx = x >> 1, vy = y >> 1;
                int v = baseP + vy*nextDim + vx;
                if (u < v) g.AddEdge(u, v);
                else       g.AddEdge(v, u);
            }
        }
    }

    return g;
}

Graph GenerateMeshOfStars(int leavesPerHub, int meshSide) {
    int numHubs    = meshSide * meshSide;
    int numCompute = numHubs * leavesPerHub;
    Graph g(numCompute + numHubs);

    // 1) помечаем первые numCompute вершин как вычислительные
    std::vector<bool> compute(numCompute + numHubs, false);
    for (int id = 0; id < numCompute; ++id) {
        compute[id] = true;
    }
    g.SetCompute(compute);

    // 2) подключаем вычислительные узлы (leaves) к своим концентраторам (hubs)
    for (int h = 0; h < numHubs; ++h) {
        int hubId    = numCompute + h;
        int baseLeaf = h * leavesPerHub;
        for (int k = 0; k < leavesPerHub; ++k) {
            g.AddEdge(hubId, baseLeaf + k);
        }
    }

    // 3) строим mesh-связи между концентраторами
    auto hubIndex = [&](int i, int j) {
        return numCompute + (i * meshSide + j);
    };
    for (int i = 0; i < meshSide; ++i) {
        for (int j = 0; j < meshSide; ++j) {
            int u = hubIndex(i, j);
            if (j + 1 < meshSide) {
                int v = hubIndex(i, j + 1);
                g.AddEdge(u, v);
            }
            if (i + 1 < meshSide) {
                int v = hubIndex(i + 1, j);
                g.AddEdge(u, v);
            }
        }
    }

    return g;
}

Graph GenerateMeshOfSpidergon(int n, int meshSide) {
    int clusterSize = 4 * n;
    int numClusters = meshSide * meshSide;
    int N = clusterSize * numClusters;
    Graph g(N);

    // 1) Внутрикластерная топология Spidergon:
    //    – кольцо из 4n вершин
    //    – диагональные «противоположные» связи k<->k+2n
    for (int c = 0; c < numClusters; ++c) {
        int base = c * clusterSize;
        // кольцо
        for (int k = 0; k < clusterSize; ++k) {
            int u = base + k;
            int v = base + ((k + 1) % clusterSize);
            g.AddEdge(u, v);
        }
        // противоположные связи
        for (int k = 0; k < 2 * n; ++k) {
            int u = base + k;
            int v = base + (k + 2 * n);
            g.AddEdge(u, v);
        }
    }

    // 2) Межкластерные mesh-связи через специальные узлы:
    //    – «вправо»: узел n  к узлу 3n у соседа справа
    //    – «вниз»:   узел 2n к узлу 0  у соседа внизу
    auto clusterBase = [&](int i, int j) {
        return (i * meshSide + j) * clusterSize;
    };
    for (int i = 0; i < meshSide; ++i) {
        for (int j = 0; j < meshSide; ++j) {
            int base = clusterBase(i, j);
            if (j + 1 < meshSide) {
                int rightBase = clusterBase(i, j + 1);
                g.AddEdge(base +     n, rightBase + 3 * n);
            }
            if (i + 1 < meshSide) {
                int downBase = clusterBase(i + 1, j);
                g.AddEdge(base + 2 * n, downBase + 0);
            }
        }
    }

    return g;
}
