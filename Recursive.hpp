#pragma once

#include <climits>
#include "Graph.hpp"


Graph GenerateWKRecursive(int k, int L) {
    if (k < 2 || L < 1) return Graph(0);

    /* ---- total number of vertices ---- */
    long long N64 = 1;
    for (int i = 0; i < L; ++i) N64 *= k;
    if (N64 > INT_MAX) return Graph(0);
    int N = static_cast<int>(N64);
    Graph g(N);                           // все вершины вычислительные

    std::vector<int> d1(L), d2(L), pref;
    // Для каждого j = 0..L-1 будем добавлять «j-flipping» и «substituting» рёбра.
    for (int j = 0; j < L; ++j) {
        int prefixLen = L - (j + 1);  // число старших цифр
        int suffixLen = j;            // число младших цифр

        // Всего вариантов старшей части — k^prefixLen
        int numPrefixes = 1;
        for (int i = 0; i < prefixLen; ++i)
            numPrefixes *= k;
        pref.resize(prefixLen);

        // Перебираем все возможные старшие части
        for (int p = 0; p < numPrefixes; ++p) {
            // Раскладываем p в систему счисления с основанием k
            int tmp = p;
            for (int i = 0; i < prefixLen; ++i) {
                pref[i] = tmp % k;
                tmp /= k;
            }

            // Перебираем все пары значений b < c из {0,1,...,k-1}
            for (int b = 0; b < k; ++b) {
                for (int c = b + 1; c < k; ++c) {
                    // Собираем первую вершину d1 и вторую d2
                    // Сначала младшие suffixLen цифр
                    for (int i = 0; i < suffixLen; ++i) {
                        d1[i] = c;
                        d2[i] = b;
                    }
                    // Затем цифра на позиции j
                    d1[suffixLen] = b;
                    d2[suffixLen] = c;
                    // Старшие prefixLen цифр
                    for (int i = suffixLen + 1; i < L; ++i) {
                        d1[i] = pref[i - (suffixLen + 1)];
                        d2[i] = pref[i - (suffixLen + 1)];
                    }

                    // Переводим «цифровой» вектор в числовой индекс
                    int idx1 = 0, idx2 = 0, power = 1;
                    for (int i = 0; i < L; ++i) {
                        idx1 += d1[i] * power;
                        idx2 += d2[i] * power;
                        power *= k;
                    }

                    // Добавляем ребро в обе стороны (для неориентированного графа)
                    g.AddEdge(idx1, idx2);
                }
            }
        }
    }
    return g;
}

Graph GenerateTHIN(int L) {
    return GenerateWKRecursive(3, L);
}

Graph GenerateRNT(int k, int L, int layers) {
    // N = k^L — количество вершин в одном слое
    int N = 1;
    for (int i = 0; i < L; ++i) N *= k;

    // Всего вершин = N * layers
    Graph g(N * layers);

    // Вспомог. контейнеры для цифровых представлений вершин
    std::vector<int> d1(L), d2(L), pref;

    // 1) Для каждого слоя строим WK-Recursive «на месте»
    for (int z = 0; z < layers; ++z) {
        int offset = z * N;  // смещение индексов текущего слоя

        // Проходим по всем уровням j = 0..L-1
        for (int j = 0; j < L; ++j) {
            int prefixLen = L - (j + 1);
            int suffixLen = j;

            // Количество вариантов старшей части
            int numPrefixes = 1;
            for (int t = 0; t < prefixLen; ++t) numPrefixes *= k;
            pref.resize(prefixLen);

            // Перебираем все старшие части
            for (int p = 0; p < numPrefixes; ++p) {
                // Раскладываем p в k-ичную старшую часть
                int tmp = p;
                for (int i = 0; i < prefixLen; ++i) {
                    pref[i] = tmp % k;
                    tmp /= k;
                }

                // Перебираем пары b < c из {0..k-1}
                for (int b = 0; b < k; ++b) {
                    for (int c = b + 1; c < k; ++c) {
                        // Заполняем младшие suffixLen цифр
                        for (int i = 0; i < suffixLen; ++i) {
                            d1[i] = c;
                            d2[i] = b;
                        }
                        // Цифра на позиции j
                        d1[suffixLen] = b;
                        d2[suffixLen] = c;
                        // Старшие prefixLen цифр
                        for (int i = suffixLen + 1; i < L; ++i) {
                            d1[i] = pref[i - (suffixLen + 1)];
                            d2[i] = pref[i - (suffixLen + 1)];
                        }

                        // Переводим в индексы и добавляем offset
                        int idx1 = 0, idx2 = 0, power = 1;
                        for (int i = 0; i < L; ++i) {
                            idx1 += d1[i] * power;
                            idx2 += d2[i] * power;
                            power *= k;
                        }
                        idx1 += offset;
                        idx2 += offset;

                        // Добавляем неориентированное ребро
                        g.AddEdge(idx1, idx2);
                    }
                }
            }
        }
    }

    int step = 1;
    for (int i = 0; i < L - 1; ++i) step *= k;  // step = k^(L-1)
    int A = 0;
    int B = step - 1;
    int C = N - step;
    int D = N - 1;

    for (int z = 0; z + 1 < layers; ++z) {
        int off1 = z * N;
        int off2 = (z + 1) * N;

        // Соединяем соответствующие углы
        int corners1[4] = {off1 + A, off1 + B, off1 + C, off1 + D};
        int corners2[4] = {off2 + A, off2 + B, off2 + C, off2 + D};
        for (int i = 0; i < 4; ++i) {
            g.AddEdge(corners1[i], corners2[i]);
        }
    }

    return g;
}

Graph GenerateRCR(int k, int r, int j) {
    // p = общее число «кубовых» измерений
    int p = k + j;
    // число адресов A = 2^p
    int M = 1 << p;
    // всего вершин
    int N = M * r;

    Graph g(N);

    // для каждого адреса A и кольцевого индекса b
    for (int A = 0; A < M; ++A) {
        for (int b = 0; b < r; ++b) {
            int u = A * r + b;

            // 1) кольцевые ребра: (b) —> (b+1)%r
            int v_ring = A * r + ((b + 1) % r);
            if (u < v_ring) {
                g.AddEdge(u, v_ring);
            }

            // 2) «кубовые» ребра: флип каждого из p битов A
            for (int bit = 0; bit < p; ++bit) {
                int A2 = A ^ (1 << bit);
                int v_cube = A2 * r + b;
                // чтобы не дублировать AddEdge
                if (u < v_cube) {
                    g.AddEdge(u, v_cube);
                }
            }
        }
    }

    return g;
}