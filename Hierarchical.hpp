#pragma once

#include "Graph.hpp"

Graph GenerateCubeConnectedCycles(int d) {
    // число вершин гиперкуба
    int M = 1 << d;
    // итоговое число вершин в CCC = M * d
    int N = M * d;
    Graph g(N);

    // Для каждого гиперкубного адреса A и позиции i на цикле:
    for (int A = 0; A < M; ++A) {
        for (int i = 0; i < d; ++i) {
            int u = A * d + i;

            // 1) ребро цикла внутри одной «копии» A:
            int v_cycle = A * d + ((i + 1) % d);
            g.AddEdge(u, v_cycle);

            // 2) соединение по гиперкубному измерению i:
            //    переходим к соседу A' = A ^ (1<<i), остаёмся в той же позиции i
            int A2 = A ^ (1 << i);
            int v_cube = A2 * d + i;
            // убедимся, что добавляем каждое ребро один раз
            if (u < v_cube) {
                g.AddEdge(u, v_cube);
            }
        }
    }

    return g;
}

Graph GenerateCubeConnectedCirculants(int hyperDim,
                                      int circBase,
                                      int circExp) {
    // 1) Параметры
    int d = hyperDim;                  // размерность гиперкуба
    int n = circBase;                  // основание «мультипликативного» циркулянта
    int r = circExp;                   // экспонента: циркулянт размера n^r
    int M = 1 << d;                    // число вершин гиперкуба = 2^d
    int C = static_cast<int>(std::pow(n, r)); // число вершин в каждом циркулянте
    int N = M * C;                     // итоговое число вершин

    Graph g(N);

    auto idx = [&](int A, int x) {
        return A * C + x;
    };

    // 2) Вычисляем сами «мультипликативные» генераторы ±n^j mod C
    std::vector<int> gens;
    gens.reserve(2 * r);
    for (int j = 0; j < r; ++j) {
        int step = static_cast<int>(std::pow(n, j)) % C;
        gens.push_back(step);              // +n^j
        gens.push_back((C - step) % C);    // −n^j ≡ C−n^j mod C
    }

    // 3) Внутренние циркулянтные рёбра
    for (int A = 0; A < M; ++A) {
        for (int x = 0; x < C; ++x) {
            int u = idx(A, x);
            for (int step: gens) {
                int y = (x + step) % C;
                int v = idx(A, y);
                if (u < v) g.AddEdge(u, v);
            }
        }
    }

    // 4) Кубовые рёбра: «сдвиг измерения на d·x» как в статье
    for (int A = 0; A < M; ++A) {
        for (int x = 0; x < C; ++x) {
            int u = idx(A, x);
            for (int i = 0; i < d; ++i) {
                // направление «i + d·x» mod d
                int di = (i + d * x) % d;
                int A2 = A ^ (1 << di);
                int v = idx(A2, x);
                if (u < v) g.AddEdge(u, v);
            }
        }
    }

    return g;
}


enum class SubnetTopo {
    Mesh, RCTM
};
enum class HubTopo {
    Ring, MeshExternal
};

// ————— Вспомогательные сборщики подсетей (фикc 4×4) —————

inline void buildSubnetMesh4x4(Graph &g, int base) {
    const int W = 4, H = 4;
    auto idx = [&](int x, int y) { return base + y * W + x; };
    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
            int u = idx(x, y);
            if (x + 1 < W) g.AddEdge(u, idx(x + 1, y));
            if (y + 1 < H) g.AddEdge(u, idx(x, y + 1));
        }
}

inline void buildSubnetRCTM4x4(Graph &g, int base) {
    const int W = 4, H = 4;
    auto idx = [&](int x, int y) { return base + y * W + x; };
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            int u = idx(x, y);
            // 1) горизонталь/вертикаль
            if (x + 1 < W) g.AddEdge(u, idx(x + 1, y));
            if (y + 1 < H) g.AddEdge(u, idx(x, y + 1));
            // 2) обе диагонали, когда они внутри поля
            if (x + 1 < W && y + 1 < H)  // юго-восточная
                g.AddEdge(u, idx(x + 1, y + 1));
            if (x - 1 >= 0 && y + 1 < H) // юго-западная
                g.AddEdge(u, idx(x - 1, y + 1));
        }
    }
}



// ————— Вспомогательные сборщики «проводных» хабов —————

inline void buildHubRing(Graph &g, int K, int base) {
    for (int i = 0; i < K; ++i) {
        int u = base + i;
        int v = base + (i + 1) % K;
        if (u < v)
            g.AddEdge(u, v);
    }
}

inline void buildHubMesh(Graph &g, int S, int base) {
    auto idx = [&](int r, int c) { return base + r * S + c; };
    for (int r = 0; r < S; ++r)
        for (int c = 0; c < S; ++c) {
            int u = idx(r, c);
            int v1 = idx(r, c + 1);
            int v2 = idx(r + 1, c);
            if (c + 1 < S && u < v1)
                g.AddEdge(u, v1);
            if (r + 1 < S && u < v2)
                g.AddEdge(u, v2);
        }
}

// ————— Выбор беспроводных хабов —————

// 1) Для кольца из n^2 хабов: номера i, где i % n == 0
static std::vector<int> selectWirelessRing(int K) {
    int n = int(std::round(std::sqrt(K)));
    std::vector<int> w;
    for (int i = 0; i < K; i += n)
        w.push_back(i);
    return w;
}

// 2) Для квадратной S×S сетки: решение задачи N-ферзей для S
static bool solveQueen(int col, int S,
                       std::vector<int> &pos,
                       std::vector<bool> &usedCol,
                       std::vector<bool> &usedDiag1,
                       std::vector<bool> &usedDiag2) {
    if (col == S) return true;
    for (int c = 0; c < S; ++c) {
        int d1 = col + c, d2 = col - c + (S - 1);
        if (!usedCol[c] && !usedDiag1[d1] && !usedDiag2[d2]) {
            pos[col] = c;
            usedCol[c] = usedDiag1[d1] = usedDiag2[d2] = true;
            if (solveQueen(col + 1, S, pos, usedCol, usedDiag1, usedDiag2))
                return true;
            usedCol[c] = usedDiag1[d1] = usedDiag2[d2] = false;
        }
    }
    return false;
}

static std::vector<int> selectWirelessMesh(int S) {
    if (S == 2)
        return {0,3};
    // pos[r] = column of queen in row r
    std::vector<int> pos(S);
    std::vector<bool> usedCol(S, false),
            usedDiag1(2 * S - 1, false),
            usedDiag2(2 * S - 1, false);
    bool ok = solveQueen(0, S, pos, usedCol, usedDiag1, usedDiag2);
    if (!ok)
        exit(2);
    std::vector<int> w;
    for (int r = 0; r < S; ++r)
        w.push_back(r * S + pos[r]);
    return w;
}

// ————— Общая функция иерархии со связями хабов —————

Graph GenerateHierarchical(
        SubnetTopo subnet,
        HubTopo hub,
        int extSize
) {
    const int W = 4, H = 4;
    int numSubnets = (hub == HubTopo::Ring ? extSize : extSize * extSize);
    int totalPEs = W * H * numSubnets;
    int baseHub = totalPEs;
    Graph g(totalPEs + numSubnets);

    // Массивы «супер‑узлов» для Mesh и RCTM
    static const int superMesh[4] = {1, 7, 8, 14};
    static const int superRCTM[4] = {5, 6, 9, 10};

    // 1) строим подсети и подключаем их к хабам по «супер‑узлам»
    for (int s = 0; s < numSubnets; ++s) {
        int basePE = s * (W * H);
        if (subnet == SubnetTopo::Mesh) buildSubnetMesh4x4(g, basePE);
        else buildSubnetRCTM4x4(g, basePE);

        int hubIdx = baseHub + s;
        // выбираем правильный массив офсетов
        const int *superOff = (subnet == SubnetTopo::Mesh
                               ? superMesh
                               : superRCTM);
        for (int i = 0; i < 4; ++i) {
            g.AddEdge(hubIdx, basePE + superOff[i]);
        }
    }

    // 2) проводные связи между хабами
    if (hub == HubTopo::Ring) buildHubRing(g, extSize, baseHub);
    else buildHubMesh(g, extSize, baseHub);

    // 3) беспроводные связи
    std::vector<int> wireless = (hub == HubTopo::Ring
                                 ? selectWirelessRing(extSize)
                                 : selectWirelessMesh(extSize));
    for (size_t i = 0; i < wireless.size(); ++i) {
        for (size_t j = i + 1; j < wireless.size(); ++j) {
            g.AddEdge(baseHub + wireless[i],
                      baseHub + wireless[j]);
        }
    }

    std::vector<bool> compute(g.NumVertices(), true);
    for (int i = totalPEs; i < g.NumVertices(); i++)
        compute[i] = false;

    g.SetCompute(compute);

    return g;
}

// ————— Удобные обёртки —————

Graph GenerateRingMesh(int ringSize) {
    return GenerateHierarchical(
            SubnetTopo::Mesh,
            HubTopo::Ring,
            ringSize
    );
}

Graph GenerateMeshMesh(int meshSide) {
    return GenerateHierarchical(
            SubnetTopo::Mesh,
            HubTopo::MeshExternal,
            meshSide
    );
}

Graph GenerateRingRCTM(int ringSize) {
    return GenerateHierarchical(
            SubnetTopo::RCTM,
            HubTopo::Ring,
            ringSize
    );
}

Graph GenerateMeshRCTM(int meshSide) {
    return GenerateHierarchical(
            SubnetTopo::RCTM,
            HubTopo::MeshExternal,
            meshSide
    );
}


Graph GenerateTESH(const std::vector<int> &sizes) {
    int L = sizes.size();  // sizes[0]=baseSize, sizes[1..]=torus sizes
    long long N = 1;
    for (int S: sizes) N *= 1LL * S * S;
    Graph g((int) N);

    int dims = 2 * L;
    std::vector<int> c(dims), n(dims);

    for (int u = 0; u < N; ++u) {
        long long tmp = u;
        // decode mixed-radix coords
        for (int lvl = 0; lvl < L; ++lvl) {
            int S = sizes[lvl];
            c[2 * lvl + 0] = tmp % S;
            tmp /= S;
            c[2 * lvl + 1] = tmp % S;
            tmp /= S;
        }

        int base = sizes[0];
        // 1) local mesh (level 1) sparse: only BL and BR get global channels,
        //    but all nodes get the 4-connectivity within baseSize×baseSize
        for (int d = 0; d < 2; ++d) {
            for (int delta: {-1, +1}) {
                int vCoord = c[d] + delta;
                if (0 <= vCoord && vCoord < base) {
                    n = c;
                    n[d] = vCoord;
                    int v = 0, m = 1;
                    for (int lvl = 0; lvl < L; ++lvl) {
                        int S = sizes[lvl];
                        v += n[2 * lvl + 0] * m;
                        m *= S;
                        v += n[2 * lvl + 1] * m;
                        m *= S;
                    }
                    if (u < v)
                        g.AddEdge(u, v);
                }
            }
        }

        bool isBL = (c[0] == 0) && (c[1] == base - 1);
        bool isBR = (c[0] == base - 1) && (c[1] == base - 1);

        // 2) sparse torus channels levels ≥2
        if (isBL || isBR) {
            for (int lvl = 1; lvl < L; ++lvl) {
                int S = sizes[lvl];
                int dx = 2 * lvl, dy = dx + 1;
                if (isBL) {
                    // vertical both ways on y_lv
                    for (int delta: {-1, +1}) {
                        n = c;
                        n[dy] = (c[dy] + delta + S) % S;
                        int v = 0, m = 1;
                        for (int ll = 0; ll < L; ++ll) {
                            int SS = sizes[ll];
                            v += n[2 * ll + 0] * m;
                            m *= SS;
                            v += n[2 * ll + 1] * m;
                            m *= SS;
                        }
                        if (u < v)
                            g.AddEdge(u, v);
                    }
                }
                if (isBR) {
                    // horizontal both ways on x_lv
                    for (int delta: {-1, +1}) {
                        n = c;
                        n[dx] = (c[dx] + delta + S) % S;
                        int v = 0, m = 1;
                        for (int ll = 0; ll < L; ++ll) {
                            int SS = sizes[ll];
                            v += n[2 * ll + 0] * m;
                            m *= SS;
                            v += n[2 * ll + 1] * m;
                            m *= SS;
                        }
                        if (u < v)
                            g.AddEdge(u, v);
                    }
                }
            }
        }
    }
    return g;
}

Graph GenerateMTESH(const std::vector<int> &sizes) {
    int L = sizes.size();
    long long N = 1;
    for (int S: sizes) N *= 1LL * S * S;
    Graph g((int) N);

    int dims = 2 * L;
    std::vector<int> c(dims), n(dims);

    for (int u = 0; u < N; ++u) {
        long long tmp = u;
        for (int lvl = 0; lvl < L; ++lvl) {
            int S = sizes[lvl];
            c[2 * lvl + 0] = tmp % S;
            tmp /= S;
            c[2 * lvl + 1] = tmp % S;
            tmp /= S;
        }

        // 1) local torus (level 1): full wrap-around on x0,y0
        int base = sizes[0];
        for (int d = 0; d < 2; ++d) {
            for (int delta: {-1, +1}) {
                n = c;
                n[d] = (c[d] + delta + base) % base;
                int v = 0, m = 1;
                for (int lvl = 0; lvl < L; ++lvl) {
                    int S = sizes[lvl];
                    v += n[2 * lvl + 0] * m;
                    m *= S;
                    v += n[2 * lvl + 1] * m;
                    m *= S;
                }
                if (u < v)
                    g.AddEdge(u, v);
            }
        }

        bool isBL = (c[0] == 0) && (c[1] == base - 1);
        bool isBR = (c[0] == base - 1) && (c[1] == base - 1);

        // 2) sparse torus channels levels ≥2
        if (isBL || isBR) {
            for (int lvl = 1; lvl < L; ++lvl) {
                int S = sizes[lvl];
                int dx = 2 * lvl, dy = dx + 1;
                if (isBL) {
                    for (int delta: {-1, +1}) {
                        n = c;
                        n[dy] = (c[dy] + delta + S) % S;
                        int v = 0, m = 1;
                        for (int ll = 0; ll < L; ++ll) {
                            int SS = sizes[ll];
                            v += n[2 * ll + 0] * m;
                            m *= SS;
                            v += n[2 * ll + 1] * m;
                            m *= SS;
                        }
                        if (u < v)
                            g.AddEdge(u, v);
                    }
                }
                if (isBR) {
                    for (int delta: {-1, +1}) {
                        n = c;
                        n[dx] = (c[dx] + delta + S) % S;
                        int v = 0, m = 1;
                        for (int ll = 0; ll < L; ++ll) {
                            int SS = sizes[ll];
                            v += n[2 * ll + 0] * m;
                            m *= SS;
                            v += n[2 * ll + 1] * m;
                            m *= SS;
                        }
                        if (u < v)
                            g.AddEdge(u, v);
                    }
                }
            }
        }
    }
    return g;
}

Graph GenerateRDT(int rows, int cols, int R, int n = 2) {
    int N = rows * cols;
    Graph g(N);

    // 1) Вычисляем векторы направлений для каждого ранга
    //    xr[r], yr[r] — это смещения по (i,j) для rank = r
    std::vector<std::pair<int, int>> xr(R + 1), yr(R + 1);
    xr[0] = {1, 0};
    yr[0] = {0, 1};
    for (int r = 1; r <= R; ++r) {
        int xrp = xr[r - 1].first, yrp = xr[r - 1].second;
        int xyp = yr[r - 1].first, yyp = yr[r - 1].second;
        // рекуррентные формулы из Definition 2
        xr[r].first = n * xrp + n * xyp;
        xr[r].second = n * yrp + n * yyp;
        yr[r].first = -n * xrp + n * xyp;
        yr[r].second = -n * yrp + n * yyp;
    }

    // 2) Для каждого узла (i,j) добавляем рёбра:
    //    — базовый тор: смещения xr[0], yr[0]
    //    — для каждого rank = 1..R: смещения ±xr[r], ±yr[r]
    auto mod = [&](int v, int m) {
        int t = v % m;
        return t < 0 ? t + m : t;
    };

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int v = i * cols + j;
            // базовый тор
            {
                auto [di, dj] = xr[0];
                int ni = mod(i + di, rows), nj = mod(j + dj, cols);
                if (v < ni * cols + nj)
                    g.AddEdge(v, ni * cols + nj);
                auto [ei, ej] = yr[0];
                ni = mod(i + ei, rows);
                nj = mod(j + ej, cols);
                if (v < ni * cols + nj)
                    g.AddEdge(v, ni * cols + nj);
            }
            // верхние ранги
            for (int r = 1; r <= R; ++r) {
                std::array<std::pair<int, int>, 4> offs = {
                        xr[r], yr[r],
                        std::make_pair(-xr[r].first, -xr[r].second),
                        std::make_pair(-yr[r].first, -yr[r].second)
                };
                for (auto [di, dj]: offs) {
                    int ni = mod(i + di, rows);
                    int nj = mod(j + dj, cols);
                    int u = ni * cols + nj;
                    // чтобы не дублировать ребра, добавляем только
                    // когда v < u
                    if (v < u) {
                        g.AddEdge(v, u);
                    }
                }
            }
        }
    }

    return g;
}


Graph GenerateMTN(int m, int level, int q) {
    int S = 1 << m;  // сторона базового торуса

    // --- Level 1: простой тор S×S ---
    if (level == 1) {
        long long N = 1LL * S * S;
        Graph g(N);
        auto modS = [&](int x) {
            int t = x % S;
            return t < 0 ? t + S : t;
        };
        for (int i = 0; i < S; ++i) {
            for (int j = 0; j < S; ++j) {
                long long v = 1LL * i * S + j;
                // вправо
                g.AddEdge(v, 1LL * i * S + modS(j + 1));
                // вниз
                g.AddEdge(v, 1LL * modS(i + 1) * S + j);
            }
        }
        return g;
    }

    // --- Level >1: рекурсивно строим предыдущий уровень ---
    Graph g_prev = GenerateMTN(m, level - 1, q);
    long long N_prev = g_prev.NumVertices();
    int M = S * S;               // число модулей
    long long N = N_prev * M;    // общее число вершин
    Graph g(N);

    // 1) Встраиваем S×S копий g_prev
    for (int mi = 0; mi < S; ++mi) {
        for (int mj = 0; mj < S; ++mj) {
            long long base = (1LL * mi * S + mj) * N_prev;
            for (long long v = 0; v < N_prev; ++v) {
                for (long long u: g_prev.Neighbors(v)) {
                    if (u > v) {
                        g.AddEdge(base + v, base + u);
                    }
                }
            }
        }
    }

    // 2) Вычисляем "смещённые" порты по Midimew
    std::vector<std::pair<int, int>> horiz, vert, diag;
    int ports = 1 << q;
    for (int t = 0; t < ports; ++t) {
        int gap = S / ports;
        int offset = gap / 2 + t * gap;
        // верх/низ
        horiz.emplace_back(0, offset);
        horiz.emplace_back(S - 1, offset);
        // лево/право
        vert.emplace_back(offset, 0);
        vert.emplace_back(offset, S - 1);
        // диагонали
        diag.emplace_back(offset, offset);
        diag.emplace_back(offset, S - 1 - offset);
    }
    auto mod = [&](int x) {
        int t = x % S;
        return t < 0 ? t + S : t;
    };

    // 3) Межмодульные связи по Midimew
    for (int mi = 0; mi < S; ++mi) {
        for (int mj = 0; mj < S; ++mj) {
            long long base = (1LL * mi * S + mj) * N_prev;

            // — горизонтальные (право)
            for (auto [i0, j0]: horiz) {
                long long u = base + (1LL * i0 * S + j0);
                int ni = mi, nj = (mj + 1) % S;
                long long v = (1LL * ni * S + nj) * N_prev + (1LL * i0 * S + j0);
                g.AddEdge(u, v);
            }

            // — вертикальные (вниз)
            for (auto [i0, j0]: vert) {
                long long u = base + (1LL * i0 * S + j0);
                int ni = (mi + 1) % S, nj = mj;
                long long v = (1LL * ni * S + nj) * N_prev + (1LL * i0 * S + j0);
                g.AddEdge(u, v);
            }

            // — диагональные (вниз-вправо и вниз-влево)
            for (auto [i0, j0]: diag) {
                long long u = base + (1LL * i0 * S + j0);

                // вниз-вправо
                {
                    int ni = (mi + 1) % S, nj = (mj + 1) % S;
                    long long v = (1LL * ni * S + nj) * N_prev + (1LL * i0 * S + j0);
                    g.AddEdge(u, v);
                }
                // вниз-влево
                {
                    int ni = (mi + 1) % S, nj = mod(mj - 1);
                    long long v = (1LL * ni * S + nj) * N_prev + (1LL * i0 * S + j0);
                    g.AddEdge(u, v);
                }
            }
        }
    }

    return g;
}


Graph GenerateCCN_HIN(const std::vector<int> &sizes) {
    // --- 1) Строим базовый полный граф K_{sizes[0]} ---
    int N0 = sizes[0];
    Graph g(1LL * N0);
    for (int i = 0; i < N0; ++i) {
        for (int j = i + 1; j < N0; ++j) {
            g.AddEdge(i, j);
        }
    }

    // --- 2) Для каждого последующего уровня i = 1..L-1:
    //     разворачиваем текущий g в sizes[i] копий и
    //     «полносвязываем» их между собой по репрезентантам.
    for (size_t lvl = 1; lvl < sizes.size(); ++lvl) {
        int C = sizes[lvl];          // число копий предыдущего графа
        long long N_prev = g.NumVertices();
        long long N_new  = N_prev * C;
        Graph g_new(N_new);

        // 2.1) Встраиваем C копий текущего g
        for (int c = 0; c < C; ++c) {
            long long base = 1LL * c * N_prev;
            // идём по каждой вершине v и её соседям u > v, чтобы не дублировать
            for (long long v = 0; v < N_prev; ++v) {
                for (long long u : g.Neighbors(v)) {
                    if (u > v) {
                        g_new.AddEdge(base + v, base + u);
                    }
                }
            }
        }

        // 2.2) Межкопийные (межкластерные) связи:
        //      для каждой пары кластеров i<j соединяем их «репрезентанты»:
        //      вершину 0 в каждой копии.
        for (int i = 0; i < C; ++i) {
            for (int j = i + 1; j < C; ++j) {
                long long u = 1LL * i * N_prev + 0;
                long long v = 1LL * j * N_prev + 0;
                g_new.AddEdge(u, v);
            }
        }

        // переходим к следующему уровню
        g = std::move(g_new);
    }

    return g;
}

Graph GenerateH3DMesh(int m, int n, int L) {
    // --- Level 1: 3D-torus m×m×m ---
    if (L == 1) {
        long long N = 1LL * m * m * m;
        Graph g(N);
        auto mod = [&](int x, int M) {
            int t = x % M;
            return t < 0 ? t + M : t;
        };
        for (int x = 0; x < m; ++x)
            for (int y = 0; y < m; ++y)
                for (int z = 0; z < m; ++z) {
                    long long v = (x * m + y) * m + z;
                    g.AddEdge(v, (mod(x+1,m)*m + y)*m + z);
                    g.AddEdge(v, (x*m + mod(y+1,m))*m + z);
                    g.AddEdge(v, (x*m + y)*m + mod(z+1,m));
                }
        return g;
    }

    // --- Level >1: рекурсивно строим уровень L-1 ---
    Graph g_prev = GenerateH3DMesh(m, n, L-1);
    long long N_prev = g_prev.NumVertices();
    int M = n * n;                        // число модулей на этом уровне
    long long N = N_prev * M;
    Graph g(N);

    // 1) встраиваем n×n копий g_prev
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            long long base = (1LL * i * n + j) * N_prev;
            for (long long v = 0; v < N_prev; ++v) {
                for (long long u : g_prev.Neighbors(v)) {
                    if (u > v) {
                        g.AddEdge(base + v, base + u);
                    }
                }
            }
        }
    }

    // 2) собираем список «свободных» портов на границе базового модуля (контуры XY-плоскости)
    std::vector<std::tuple<int,int,int>> freePorts;
    // y = 0 и y = m-1 (сев/юг)
    for (int x = 0; x < m; ++x)
        for (int z = 0; z < m; ++z) {
            freePorts.emplace_back(x, 0,   z);
            freePorts.emplace_back(x, m-1, z);
        }
    // x = 0 и x = m-1 (зап/вост), чтобы не дублировать углы, y от 1 до m-2
    for (int y = 1; y < m-1; ++y)
        for (int z = 0; z < m; ++z) {
            freePorts.emplace_back(0,   y, z);
            freePorts.emplace_back(m-1, y, z);
        }

    // 3) межмодульная 2D-mesh-связность (без wrap-around)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            long long base = (1LL * i * n + j) * N_prev;
            // право (j+1)
            if (j + 1 < n) {
                long long baseE = (1LL * i * n + (j+1)) * N_prev;
                for (auto [x,y,z] : freePorts) {
                    long long u = base + ((x*m + y)*m + z);
                    long long v = baseE + ((x*m + y)*m + z);
                    if (u < v)
                        g.AddEdge(u, v);
                }
            }
            // вниз (i+1)
            if (i + 1 < n) {
                long long baseS = (1LL * (i+1) * n + j) * N_prev;
                for (auto [x,y,z] : freePorts) {
                    long long u = base + ((x*m + y)*m + z);
                    long long v = baseS + ((x*m + y)*m + z);
                    if (u < v)
                        g.AddEdge(u, v);
                }
            }
        }
    }

    return g;
}

Graph GenerateH3DTorus(int m, int n, int L) {
    // --- Level 1: простой 3D-mesh (без wrap-around) ---
    if (L == 1) {
        long long N = 1LL * m * m * m;
        Graph g(N);
        auto mod = [&](int x, int M) {
            int t = x % M;
            return t < 0 ? t + M : t;
        };
        for (int x = 0; x < m; ++x) {
            for (int y = 0; y < m; ++y) {
                for (int z = 0; z < m; ++z) {
                    long long v = (x * m + y) * m + z;
                    // +x
                    if (x + 1 < m)
                        g.AddEdge(v, ((x+1)*m + y)*m + z);
                    // +y
                    if (y + 1 < m)
                        g.AddEdge(v, (x*m + (y+1))*m + z);
                    // +z
                    if (z + 1 < m)
                        g.AddEdge(v, (x*m + y)*m + (z+1));
                }
            }
        }
        return g;
    }

    // --- Level >1: рекурсивно строим уровень L−1 ---
    Graph g_prev = GenerateH3DTorus(m, n, L-1);
    long long N_prev = g_prev.NumVertices();
    // число модулей = n^3
    long long N = N_prev * n * n * n;
    Graph g(N);

    // 1) встраиваем n×n×n копий g_prev
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                long long base = ((1LL*i*n + j)*n + k) * N_prev;
                for (long long v = 0; v < N_prev; ++v) {
                    for (long long u : g_prev.Neighbors(v)) {
                        if (u > v) {
                            g.AddEdge(base + v, base + u);
                        }
                    }
                }
            }
        }
    }

    // 2) собираем «свободные» порты: узлы на 4 вертикальных ребрах BM:
    //    (0,0,z), (0,m-1,z), (m-1,0,z), (m-1,m-1,z), z=0..m-1
    std::vector<std::tuple<int,int,int>> freePorts;
    for (int z = 0; z < m; ++z) {
        freePorts.emplace_back(0,    0,    z);
        freePorts.emplace_back(0,    m-1,  z);
        freePorts.emplace_back(m-1,  0,    z);
        freePorts.emplace_back(m-1,  m-1,  z);
    }

    auto modn = [&](int x) {
        int t = x % n;
        return t < 0 ? t + n : t;
    };

    // 3) межмодульные связи: 3D-torus по осям i, j, k (wrap-around)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                long long base = ((1LL*i*n + j)*n + k) * N_prev;

                // сосед по +X-торусу
                {
                    int ni = modn(i+1), nj = j, nk = k;
                    long long base2 = ((1LL*ni*n + nj)*n + nk) * N_prev;
                    for (auto [x,y,z] : freePorts) {
                        long long u = base  + ((x*m + y)*m + z);
                        long long v = base2 + ((x*m + y)*m + z);
                        g.AddEdge(u, v);
                    }
                }
                // сосед по +Y-торусу
                {
                    int ni = i, nj = modn(j+1), nk = k;
                    long long base2 = ((1LL*ni*n + nj)*n + nk) * N_prev;
                    for (auto [x,y,z] : freePorts) {
                        long long u = base  + ((x*m + y)*m + z);
                        long long v = base2 + ((x*m + y)*m + z);
                        g.AddEdge(u, v);
                    }
                }
                // сосед по +Z-торусу
                {
                    int ni = i, nj = j, nk = modn(k+1);
                    long long base2 = ((1LL*ni*n + nj)*n + nk) * N_prev;
                    for (auto [x,y,z] : freePorts) {
                        long long u = base  + ((x*m + y)*m + z);
                        long long v = base2 + ((x*m + y)*m + z);
                        g.AddEdge(u, v);
                    }
                }
            }
        }
    }

    return g;
}


Graph GenerateHTN(int m, int n, int L) {
    // --- Level 1: полный 3D-torus m×m×m ---
    if (L == 1) {
        long long N = 1LL * m * m * m;
        Graph g(N);
        auto mod = [&](int x, int M) {
            int t = x % M; return t < 0 ? t + M : t;
        };
        for (int x = 0; x < m; ++x) {
            for (int y = 0; y < m; ++y) {
                for (int z = 0; z < m; ++z) {
                    long long v = (1LL*x*m + y)*m + z;
                    // +x, −x
                    g.AddEdge(v, (1LL*mod(x+1,m)*m + y)*m + z);
                    // +y, −y
                    g.AddEdge(v, (1LL*x*m + mod(y+1,m))*m + z);
                    // +z, −z
                    g.AddEdge(v, (1LL*x*m + y)*m + mod(z+1,m));
                }
            }
        }
        return g;
    }

    // --- Level >1: рекурсивно строим уровень L−1 ---
    Graph g_prev = GenerateHTN(m, n, L-1);
    long long N_prev = g_prev.NumVertices();
    int M = n * n;
    long long N = N_prev * M;
    Graph g(N);

    // 1) встраиваем n×n копий g_prev
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            long long base = (1LL*i*n + j) * N_prev;
            for (long long v = 0; v < N_prev; ++v) {
                for (long long u : g_prev.Neighbors(v)) {
                    if (u > v) {
                        g.AddEdge(base + v, base + u);
                    }
                }
            }
        }
    }

    // 2) собираем «свободные» порты на контуре базового модуля (PEs с y=0,y=m−1 или x=0,x=m−1)
    std::vector<std::tuple<int,int,int>> freePorts;
    // север/юг (y=0, y=m−1)
    for (int x = 0; x < m; ++x)
        for (int z = 0; z < m; ++z) {
            freePorts.emplace_back(x, 0,   z);
            freePorts.emplace_back(x, m-1, z);
        }
    // запад/восток (x=0, x=m−1), y=1..m−2 чтобы не дублировать углы
    for (int y = 1; y < m-1; ++y)
        for (int z = 0; z < m; ++z) {
            freePorts.emplace_back(0,   y, z);
            freePorts.emplace_back(m-1, y, z);
        }

    // 3) межмодульная 2D-torus-связь (wrap-around по oсям i,j)
    auto modn = [&](int x) {
        int t = x % n; return t < 0 ? t + n : t;
    };
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            long long base = (1LL*i*n + j) * N_prev;
            // сосед по +X-торусу (i+1 mod n, j)
            {
                long long base2 = (1LL*modn(i+1)*n + j) * N_prev;
                for (auto [x,y,z] : freePorts) {
                    long long u = base  + ((1LL*x*m + y)*m + z);
                    long long v = base2 + ((1LL*x*m + y)*m + z);
                    g.AddEdge(u, v);
                }
            }
            // сосед по +Y-торусу (i, j+1 mod n)
            {
                long long base2 = (1LL*i*n + modn(j+1)) * N_prev;
                for (auto [x,y,z] : freePorts) {
                    long long u = base  + ((1LL*x*m + y)*m + z);
                    long long v = base2 + ((1LL*x*m + y)*m + z);
                    g.AddEdge(u, v);
                }
            }
        }
    }

    return g;
}

Graph GenerateMH3DT(int m, int n, int L, int q) {
    // Level 1: полный 3D-torus m×m×m (wrap-around по всем трём осям)
    if (L == 1) {
        long long N = 1LL * m * m * m;
        Graph g(N);
        auto mod = [&](int x, int M) {
            int t = x % M; return t < 0 ? t + M : t;
        };
        for (int x = 0; x < m; ++x) {
            for (int y = 0; y < m; ++y) {
                for (int z = 0; z < m; ++z) {
                    long long v = (1LL*x*m + y)*m + z;
                    // три направления +x,+y,+z достаточно, AddEdge — неориентированный
                    g.AddEdge(v, (mod(x+1,m)*m + y)*m + z);
                    g.AddEdge(v, (x*m + mod(y+1,m))*m + z);
                    g.AddEdge(v, (x*m + y)*m + mod(z+1,m));
                }
            }
        }
        return g;
    }

    // Level >1: рекурсивно строим L-1
    Graph g_prev = GenerateMH3DT(m, n, L - 1, q);
    long long N_prev = g_prev.NumVertices();
    int M = n * n * n;
    long long N = N_prev * M;
    Graph g(N);

    // 1) Встраиваем n^3 копий g_prev
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                long long base = ((1LL*i*n + j)*n + k) * N_prev;
                for (long long v = 0; v < N_prev; ++v) {
                    for (long long u : g_prev.Neighbors(v)) {
                        if (u > v) {
                            g.AddEdge(base + v, base + u);
                        }
                    }
                }
            }
        }
    }

    // 2) Собираем free-ports на вертикальных ребрах BM (четыре ребра по (x,y) = {0,m-1}^2)
    //    и выбираем по q узлов вдоль z на каждом ребре равномерно.
    std::vector<std::tuple<int,int,int>> freePorts;
    int ports = std::min(q, 2);  // q=0..2
    for (int t = 0; t < ports; ++t) {
        int gap = m / (ports + 1);
        int offset = gap/2 + t * gap;
        for (auto [x,y] : std::vector<std::pair<int,int>>{{0,0},{0,m-1},{m-1,0},{m-1,m-1}}) {
            freePorts.emplace_back(x, y, offset);
        }
    }

    auto modn = [&](int x) {
        int t = x % n; return t < 0 ? t + n : t;
    };

    // 3) Межмодульная 3D-torus-связность по осям i,j,k с wrap-around
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                long long base = ((1LL*i*n + j)*n + k) * N_prev;

                // +i
                {
                    long long base2 = ((1LL*modn(i+1)*n + j)*n + k) * N_prev;
                    for (auto [x,y,z] : freePorts) {
                        long long u = base  + ((x*m + y)*m + z);
                        long long v = base2 + ((x*m + y)*m + z);
                        g.AddEdge(u, v);
                    }
                }
                // +j
                {
                    long long base2 = ((1LL*i*n + modn(j+1))*n + k) * N_prev;
                    for (auto [x,y,z] : freePorts) {
                        long long u = base  + ((x*m + y)*m + z);
                        long long v = base2 + ((x*m + y)*m + z);
                        g.AddEdge(u, v);
                    }
                }
                // +k
                {
                    long long base2 = ((1LL*i*n + j)*n + modn(k+1)) * N_prev;
                    for (auto [x,y,z] : freePorts) {
                        long long u = base  + ((x*m + y)*m + z);
                        long long v = base2 + ((x*m + y)*m + z);
                        g.AddEdge(u, v);
                    }
                }
            }
        }
    }

    return g;
}

Graph GenerateStarMesh(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    for (int i = 1; i < rows; i += 3) {
        for (int j = 1; j < cols; j += 3) {
            g.AddEdge(i * cols + j, (i - 1) * cols + j - 1);
            g.AddEdge(i * cols + j, (i - 1) * cols + j + 1);
            g.AddEdge(i * cols + j, (i + 1) * cols + j - 1);
            g.AddEdge(i * cols + j, (i + 1) * cols + j + 1);

            if (i + 3 < rows) {
                g.AddEdge(i * cols + j, (i + 3) * cols + j);
            }
            if (j + 3 < cols) {
                g.AddEdge(i * cols + j, i * cols + j + 3);
            }
        }
    }
    return g;
}

Graph GenerateL2Star(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    for (uint16_t i = 1; i < rows; i += 3) {
        for (uint16_t j = 1; j < cols; j += 3) {
            if (i + 3 < rows) {
                g.AddEdge(i * cols + j, (i + 3) * cols + j);
            }
            if (j + 3 < cols) {
                g.AddEdge(i * cols + j, i * cols + j + 3);
            }
        }
    }

    for (uint16_t i = 2; i < rows; i += 3) {
        for (uint16_t j = 0; j < cols; j += 3) {
            if (i + 3 < rows) {
                g.AddEdge(i * cols + j, (i + 3) * cols + j);
            }
            if (j + 3 < cols) {
                g.AddEdge(i * cols + j, i * cols + j + 3);
            }
        }
    }
    return g;
}
