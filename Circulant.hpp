#pragma once

#include <set>
#include "Graph.hpp"

Graph GenerateCirculant(int num_vertices, const std::vector<int> &generators) {
    Graph g(num_vertices);
    for (int i = 0; i < num_vertices; i++) {
        for (auto x: generators) {
            g.AddEdge(i, (i + x) % num_vertices);
        }
    }
    return g;
}

bool IsPrime(uint16_t x) {
    if (x < 2) return false;
    for (uint16_t i = 2; i * i <= x; ++i)
        if (x % i == 0) return false;
    return true;
}

Graph GeneratePaley(int num_vertices, int degree) {
    if (!IsPrime(num_vertices))
        return Graph(0);

    switch (degree) {
        case 2:
            if (num_vertices % 4 != 1)
                return Graph(0);
            break;
        case 3:
            if (num_vertices % 3 != 1)
                return Graph(0);
            break;
        case 4:
            if (num_vertices % 8 != 1)
                return Graph(0);
            break;
        default:
            return Graph(0);
    }

    std::set<uint64_t> residues;
    for (uint64_t i = 1; i < num_vertices; ++i) {
        uint64_t prod = 1;
        for (int j = 0; j < degree; j++)
            prod *= i;
        uint64_t res = prod % num_vertices;
        residues.insert(res);
    }

    std::vector<int> generators;
    for (uint64_t s: residues) {
        if (s != 0 && s <= num_vertices / 2 && residues.count(num_vertices - s)) {
            generators.push_back(static_cast<int>(s));
        }
        if (residues.count(num_vertices - s) == 0) {
            std::cout << "!!!" << s << "\n";
        }
    }

    std::sort(generators.begin(), generators.end());
    return GenerateCirculant(num_vertices, generators);
}

Graph GenerateMultiplicativeCirculant(int base, int degree) {
    std::vector<int> generators(degree);
    int prod = 1;
    for (int i = 0; i < degree; i++) {
        generators[i] = prod;
        prod *= base;
    }
    return GenerateCirculant(prod * base, generators);
}

Graph GenerateRing(int num_vertices) {
    Graph g(num_vertices);
    for (int i = 0; i < num_vertices; i++) {
        g.AddEdge(i, (i + 1) % num_vertices);
    }
    return g;
}

Graph GenerateSpidergon(int num_vertices) {
    Graph g = GenerateRing(num_vertices);
    for (int i = 0; i < num_vertices; i++) {
        if (i < num_vertices / 2 == 0) {
            g.AddEdge(i, (i + num_vertices / 2) % num_vertices);
        }
    }
    return g;
}

Graph GenerateHyperCirculant(int num_vertices, const std::vector<std::vector<int>> &genx) {
    Graph g(num_vertices);
    int degree = genx.size();
    for (int i = 0; i < degree; ++i) {
        for (int node = i; node < num_vertices; node += degree) {
            for (int gen: genx[i]) {
                g.AddEdge(node, (node + gen + num_vertices) % num_vertices, true);
                g.AddEdge(node, (node - gen + num_vertices) % num_vertices, true);
            }
        }
    }
    return g;
}

Graph GenerateMidimew(int N, int b) {
    return GenerateCirculant(N, {b, b - 1});
}


inline long long CoordKey(int x, int y) {
    return (static_cast<long long>(x) << 32) |
           static_cast<unsigned int>(y);
}

Graph GenerateLNetwork(int a, int b, int c, int d) {

    long long det = static_cast<long long>(a) * d
                    - static_cast<long long>(b) * c;
    const int N = std::abs(det);
    Graph g(N);


    std::vector<std::pair<int, int>> rep;
    rep.reserve(N);


    int xs[4] = {0, a, b, a + b};
    int ys[4] = {0, c, d, c + d};
    int x0 = *std::min_element(xs, xs + 4);
    int x1 = *std::max_element(xs, xs + 4);
    int y0 = *std::min_element(ys, ys + 4);
    int y1 = *std::max_element(ys, ys + 4);


    for (int x = x0; x < x1; ++x) {
        for (int y = y0; y < y1; ++y) {


            long long s_num = static_cast<long long>(d) * x
                              - static_cast<long long>(b) * y;
            long long t_num = -static_cast<long long>(c) * x
                              + static_cast<long long>(a) * y;
            if (det > 0) {
                if (0 <= s_num && s_num < det
                    && 0 <= t_num && t_num < det)
                    rep.emplace_back(x, y);
            } else {

                if (det < s_num && s_num <= 0
                    && det < t_num && t_num <= 0)
                    rep.emplace_back(x, y);
            }
        }
    }


    std::unordered_map<long long, int> pos2id;
    pos2id.reserve(N * 2);
    for (int id = 0; id < N; ++id) {
        pos2id[CoordKey(rep[id].first, rep[id].second)] = id;
    }


    const int stepX[4] = {1, -1, 0, 0};
    const int stepY[4] = {0, 0, 1, -1};


    for (int v = 0; v < N; ++v) {
        int x = rep[v].first;
        int y = rep[v].second;

        for (int dir = 0; dir < 4; ++dir) {
            int nx = x + stepX[dir];
            int ny = y + stepY[dir];


            bool found = false;
            for (int p = -1; p <= 1 && !found; ++p) {
                for (int q = -1; q <= 1; ++q) {
                    int rx = nx - (a * p + b * q);
                    int ry = ny - (c * p + d * q);

                    auto it = pos2id.find(CoordKey(rx, ry));
                    if (it != pos2id.end()) {
                        int u = it->second;

                        if (v < u) g.AddEdge(v, u);
                        found = true;
                        break;
                    }
                }
            }
        }
    }

    return g;
}


Graph GenerateGaussian(int a, int b) {
    return GenerateLNetwork(a, -b, b, a);
}

Graph GenerateTwistedTorus(int width, int height, int twist) {
    return GenerateLNetwork(width, twist, 0, height);
}

Graph GenerateDiagonalToroidalMesh(int k1, int k2) {
    if ((k1 % 2 == 1) && (k2 % 2 == 1)) {
        int p = (k1 + k2) / 2;
        int q = (k1 - k2) / 2;
        return GenerateLNetwork(p, q, q, p);
    }
    int p = k1 / 2;
    int q = k2 / 2;
    return GenerateLNetwork(p, -q, q, p);
}


