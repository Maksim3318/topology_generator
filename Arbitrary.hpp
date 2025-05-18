#pragma once

#include <unordered_map>
#include <tuple>
#include <climits>

#include "Graph.hpp"

Graph GenerateHoneycombMesh(int t) {

    std::vector<std::tuple<int, int, int>> coords;
    coords.reserve(6 * t * t);
    for (int x = -t + 1; x <= t; ++x)
        for (int y = -t + 1; y <= t; ++y)
            for (int z = -t + 1; z <= t; ++z)
                if (1 <= x + y + z && x + y + z <= 2)
                    coords.emplace_back(x, y, z);
    int N = coords.size();
    Graph g(N);
    auto pack = [&](int x, int y, int z) {
        return ((long long) x << 40) ^ ((long long) y << 20) ^ (long long) z;
    };
    std::unordered_map<long long, int> id;
    id.reserve(N);
    for (int i = 0; i < N; ++i)
        id[pack(std::get<0>(coords[i]), std::get<1>(coords[i]), std::get<2>(coords[i]))] = i;
    auto inside = [&](int x, int y, int z) {
        return x >= -t + 1 && x <= t
               && y >= -t + 1 && y <= t
               && z >= -t + 1 && z <= t
               && 1 <= x + y + z && x + y + z <= 2;
    };
    int D[3][3] = {{1, 0, 0},
                   {0, 1, 0},
                   {0, 0, 1}};
    for (int v = 0; v < N; ++v) {
        int x = std::get<0>(coords[v]);
        int y = std::get<1>(coords[v]);
        int z = std::get<2>(coords[v]);
        for (int a = 0; a < 3; ++a)
            for (int dir: {-1, 1}) {
                int nx = x + dir * D[a][0];
                int ny = y + dir * D[a][1];
                int nz = z + dir * D[a][2];
                if (!inside(nx, ny, nz)) continue;
                int u = id[pack(nx, ny, nz)];
                if (v < u) g.AddEdge(v, u);
            }
    }
    return g;
}

Graph GenerateHoneycombTorus(int t) {

    std::vector<std::tuple<int, int, int>> coords;
    coords.reserve(6 * t * t);
    for (int x = -t + 1; x <= t; ++x)
        for (int y = -t + 1; y <= t; ++y)
            for (int z = -t + 1; z <= t; ++z)
                if (1 <= x + y + z && x + y + z <= 2)
                    coords.emplace_back(x, y, z);
    int N = coords.size();
    Graph g(N);
    auto pack = [&](int x, int y, int z) {
        return ((long long) x << 40) ^ ((long long) y << 20) ^ (long long) z;
    };
    std::unordered_map<long long, int> id;
    id.reserve(N);
    for (int i = 0; i < N; ++i)
        id[pack(std::get<0>(coords[i]), std::get<1>(coords[i]), std::get<2>(coords[i]))] = i;
    auto inside = [&](int x, int y, int z) {
        return x >= -t + 1 && x <= t
               && y >= -t + 1 && y <= t
               && z >= -t + 1 && z <= t
               && 1 <= x + y + z && x + y + z <= 2;
    };
    int D[3][3] = {{1, 0, 0},
                   {0, 1, 0},
                   {0, 0, 1}};
    int k = 2 * t - 1;
    auto wrap = [&](int axis, int dir) -> std::tuple<int, int, int> {
        if (axis == 0) return {dir > 0 ? k : -k, dir > 0 ? -t : t, dir > 0 ? -t : t};
        if (axis == 1) return {dir > 0 ? -t : t, dir > 0 ? k : -k, dir > 0 ? -t : t};
        return {dir > 0 ? -t : t, dir > 0 ? -t : t, dir > 0 ? k : -k};
    };
    for (int v = 0; v < N; ++v) {
        int x = std::get<0>(coords[v]);
        int y = std::get<1>(coords[v]);
        int z = std::get<2>(coords[v]);
        for (int a = 0; a < 3; ++a)
            for (int dir: {-1, 1}) {
                int nx = x + dir * D[a][0];
                int ny = y + dir * D[a][1];
                int nz = z + dir * D[a][2];
                if (!inside(nx, ny, nz)) {
                    auto [wx, wy, wz] = wrap(a, dir);
                    nx = x + wx;
                    ny = y + wy;
                    nz = z + wz;
                }
                int u = id[pack(nx, ny, nz)];
                if (v < u) g.AddEdge(v, u);
            }
    }
    return g;
}


Graph GenerateHexStarMesh(int rows, int cols, int d) {

    auto cubeDist = [](int x, int y, int z) {
        return (std::abs(x) + std::abs(y) + std::abs(z)) / 2;
    };
    auto pack = [](int x, int y, int z) {
        return ((long long) x << 40) ^ ((long long) y << 20) ^ (long long) z;
    };


    std::vector<std::array<int, 3>> base;
    for (int x = -2; x <= 2; ++x)
        for (int y = -2; y <= 2; ++y) {
            int z = -x - y;
            if (std::abs(z) > 2) continue;
            int r = cubeDist(x, y, z);
            if (r == 1 || r == 2) base.push_back({x, y, z});
        }


    int Rlayer = d - 1;
    std::vector<std::array<int, 3>> local;
    for (int q = -Rlayer; q <= Rlayer; ++q)
        for (int r = -Rlayer; r <= Rlayer; ++r) {
            int s = -q - r;
            if (std::max({std::abs(q), std::abs(r), std::abs(s)}) > Rlayer) continue;
            int cx = 3 * q, cy = 3 * r, cz = 3 * s;
            for (auto &v: base) local.push_back({v[0] + cx, v[1] + cy, v[2] + cz});
        }
    std::sort(local.begin(), local.end());
    local.erase(std::unique(local.begin(), local.end()), local.end());
    int nLocal = static_cast<int>(local.size());


    constexpr int dir[8][3] = {
            {1,  -1, 0},
            {1,  0,  -1},
            {0,  1,  -1},
            {-1, 1,  0},
            {-1, 0,  1},
            {0,  -1, 1},
            {2,  -1, -1},
            {-2, 1,  1}
    };
    std::unordered_map<long long, int> l2i;
    for (int i = 0; i < nLocal; ++i)
        l2i[pack(local[i][0], local[i][1], local[i][2])] = i;
    std::vector<std::pair<int, int>> localEdge;
    for (int v = 0; v < nLocal; ++v) {
        int x = local[v][0], y = local[v][1], z = local[v][2];
        for (auto dxyz: dir) {
            int nx = x + dxyz[0], ny = y + dxyz[1], nz = z + dxyz[2];
            auto it = l2i.find(pack(nx, ny, nz));
            if (it != l2i.end() && v < it->second)
                localEdge.emplace_back(v, it->second);
        }
    }


    int Rhex = 3 * d - 1;
    int idxE = l2i[pack(Rhex, -Rhex, 0)];
    int idxW = l2i[pack(-Rhex, Rhex, 0)];
    int idxSE = l2i[pack(0, Rhex, -Rhex)];
    int idxNW = l2i[pack(0, -Rhex, Rhex)];


    int spacing = 2 * Rhex;
    const int east[3] = {spacing, -spacing, 0};
    const int se[3] = {0, spacing, -spacing};


    std::unordered_map<long long, int> g2i;
    std::vector<std::array<int, 3>> global;
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int cx = col * east[0] + row * se[0];
            int cy = col * east[1] + row * se[1];
            int cz = col * east[2] + row * se[2];
            for (auto &v: local) {
                long long key = pack(v[0] + cx, v[1] + cy, v[2] + cz);
                if (g2i.find(key) == g2i.end()) {
                    int id = static_cast<int>(global.size());
                    g2i[key] = id;
                    global.push_back({v[0] + cx, v[1] + cy, v[2] + cz});
                }
            }
        }
    }
    Graph g(static_cast<int>(global.size()));


    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int cx = col * east[0] + row * se[0];
            int cy = col * east[1] + row * se[1];
            int cz = col * east[2] + row * se[2];

            auto toGlobal = [&](int localIdx) {
                auto &v = local[localIdx];
                return g2i[pack(v[0] + cx, v[1] + cy, v[2] + cz)];
            };
            for (auto [a, b]: localEdge) {
                int u = toGlobal(a), v = toGlobal(b);
                if (u < v) g.AddEdge(u, v);
            }
        }
    }


    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int cx = col * east[0] + row * se[0];
            int cy = col * east[1] + row * se[1];
            int cz = col * east[2] + row * se[2];

            int centerE_x = cx + east[0], centerE_y = cy + east[1], centerE_z = cz + east[2];
            int centerSE_x = cx + se[0], centerSE_y = cy + se[1], centerSE_z = cz + se[2];


            if (col + 1 < cols) {
                int vE = g2i[pack(local[idxE][0] + cx, local[idxE][1] + cy, local[idxE][2] + cz)];
                int uW = g2i[pack(local[idxW][0] + centerE_x, local[idxW][1] + centerE_y, local[idxW][2] + centerE_z)];
                if (vE < uW) g.AddEdge(vE, uW);
            }

            if (row + 1 < rows) {
                int vSE = g2i[pack(local[idxSE][0] + cx, local[idxSE][1] + cy, local[idxSE][2] + cz)];
                int uNW = g2i[pack(local[idxNW][0] + centerSE_x, local[idxNW][1] + centerSE_y,
                                   local[idxNW][2] + centerSE_z)];
                if (vSE < uNW) g.AddEdge(vSE, uNW);
            }
        }
    }
    return g;
}


Graph GenerateDeBruijn(int d, int k) {

    if (d < 2 || k < 1) return Graph(0);


    int N = 1;
    for (int i = 0; i < k; ++i) {
        if (N > INT_MAX / d) return Graph(0);
        N *= d;
    }
    Graph g(N);

    const int step = N / d;

    for (int v = 0; v < N; ++v) {
        int v_shiftR = v / d;


        int baseL = v * d % N;
        for (int r = 0; r < d; ++r) {
            int u = baseL + r;
            if (v < u) g.AddEdge(v, u);
        }


        for (int r = 0; r < d; ++r) {
            int u = v_shiftR + r * step;
            if (v < u) g.AddEdge(v, u);
        }
    }
    return g;
}

Graph Generate3DDeBruijn(int d, int k, int layers) {
    if (d < 2 || k < 1 || layers < 2) return Graph(0);


    long long n_plane = 1;
    for (int i = 0; i < k; ++i) {
        n_plane *= d;
        if (n_plane > INT_MAX) return Graph(0);
    }
    int Nplane = static_cast<int>(n_plane);
    int Ntotal = Nplane * layers;
    Graph g(Ntotal);

    const int step = Nplane / d;

    for (int z = 0; z < layers; ++z) {
        int base = z * Nplane;
        for (int v = 0; v < Nplane; ++v) {
            int vR = v / d;
            int baseL = (v * d) % Nplane;

            for (int r = 0; r < d; ++r) {
                int u = baseL + r;
                if (v < u) g.AddEdge(base + v, base + u);
            }
            for (int r = 0; r < d; ++r) {
                int u = vR + r * step;
                if (v < u) g.AddEdge(base + v, base + u);
            }
        }
    }

    for (int z = 0; z < layers - 1; ++z) {
        int offsetZ = z * Nplane;
        int offsetZ1 = (z + 1) * Nplane;
        for (int v = 0; v < Nplane; ++v)
            g.AddEdge(offsetZ + v, offsetZ1 + v);
    }
    return g;
}

Graph GenerateDIMB(int d, int k) {
    if (d < 2 || k < 1) return Graph(0);

    int N1 = 1;
    for (int i = 0; i < k; ++i) {
        N1 *= d;
    }
    int N = N1 * N1;
    Graph g(N);
    int step = N1 / d;

    for (int y = 0; y < N1; ++y) {
        for (int x = 0; x < N1; ++x) {
            int src = y * N1 + x;
            std::vector<int> nbrs;
            nbrs.reserve(4 * d);

            int baseL = (x * d) % N1;
            int xR = x / d;
            for (int r = 0; r < d; ++r) {
                nbrs.push_back(y * N1 + (baseL + r));
                nbrs.push_back(y * N1 + (xR + r * step));
            }

            int baseC = (y * d) % N1;
            int yR = y / d;
            for (int r = 0; r < d; ++r) {
                nbrs.push_back((baseC + r) * N1 + x);
                nbrs.push_back((yR + r * step) * N1 + x);
            }

            std::sort(nbrs.begin(), nbrs.end());
            nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
            for (int u: nbrs) {
                if (u != src) {
                    g.AddEdge(src, u, true);
                }
            }
        }
    }
    return g;
}

Graph GenerateHERT(int k, int m) {
    int N = k * m * 4;
    Graph g(N);

    auto id = [&](int i, int j, int l) {
        return ((i * m + j) * 4 + l);
    };

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int l = 0; l < 4; ++l) {
                int v = id(i, j, l);


                int l1 = (l + 1) % 4;
                int l2 = (l + 3) % 4;
                int u1 = id(i, j, l1);
                int u2 = id(i, j, l2);
                if (v < u1) g.AddEdge(v, u1);
                if (v < u2) g.AddEdge(v, u2);

                if ((l & 1) == 0) {
                    int jn = (j + 1) % m;
                    int jp = (j - 1 + m) % m;
                    int uc = id(i, jn, l);
                    int ud = id(i, jp, l);
                    if (v < uc) g.AddEdge(v, uc);
                    if (v < ud) g.AddEdge(v, ud);
                } else {
                    int in = (i + 1) % k;
                    int ip = (i - 1 + k) % k;
                    int ur = id(in, j, l);
                    int us = id(ip, j, l);
                    if (v < ur) g.AddEdge(v, ur);
                    if (v < us) g.AddEdge(v, us);
                }
            }
        }
    }

    return g;
}

Graph GenerateSEM(int R) {
    if (R < 2 || (R & (R - 1)) != 0) return Graph(0);
    int n = 0;
    while ((1 << n) < R) ++n;
    int N = R * R;
    Graph g(N);

    auto shuffle = [&](int x) {
        return ((x << 1) & (R - 1)) | (x >> (n - 1));
    };

    for (int r = 0; r < R; ++r) {
        for (int c = 0; c < R; ++c) {
            int src = r * R + c;

            int dst = r * R + shuffle(c);
            if (dst != src) g.AddEdge(src, dst, true);
            dst = r * R + (c ^ 1);
            if (dst != src) g.AddEdge(src, dst, true);

            dst = shuffle(r) * R + c;
            if (dst != src) g.AddEdge(src, dst, true);
            dst = (r ^ 1) * R + c;
            if (dst != src) g.AddEdge(src, dst, true);
        }
    }

    return g;
}

Graph GenerateHyperX(const std::vector<int> &dims) {
    int m = dims.size();
    long long N = 1;
    for (int d: dims) N *= d;
    Graph g((int) N);

    std::vector<long long> mult(m);
    mult[0] = 1;
    for (int i = 1; i < m; ++i) {
        mult[i] = mult[i - 1] * dims[i - 1];
    }

    for (int d = 0; d < m; ++d) {
        long long stride = mult[d];
        long long blockLen = stride * dims[d];
        for (long long base = 0; base < N; base += blockLen) {
            for (long long offset = 0; offset < stride; ++offset) {
                for (int i = 0; i < dims[d]; ++i) {
                    for (int j = i + 1; j < dims[d]; ++j) {
                        int u = (int) (base + offset + i * stride);
                        int v = (int) (base + offset + j * stride);
                        g.AddEdge(u, v);
                    }
                }
            }
        }
    }

    return g;
}

