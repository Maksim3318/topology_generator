#pragma once

#include <climits>
#include "Graph.hpp"


Graph GenerateWKRecursive(int k, int L) {
    if (k < 2 || L < 1) return Graph(0);


    long long N64 = 1;
    for (int i = 0; i < L; ++i) N64 *= k;
    if (N64 > INT_MAX) return Graph(0);
    int N = static_cast<int>(N64);
    Graph g(N);

    std::vector<int> d1(L), d2(L), pref;

    for (int j = 0; j < L; ++j) {
        int prefixLen = L - (j + 1);
        int suffixLen = j;


        int numPrefixes = 1;
        for (int i = 0; i < prefixLen; ++i)
            numPrefixes *= k;
        pref.resize(prefixLen);


        for (int p = 0; p < numPrefixes; ++p) {

            int tmp = p;
            for (int i = 0; i < prefixLen; ++i) {
                pref[i] = tmp % k;
                tmp /= k;
            }


            for (int b = 0; b < k; ++b) {
                for (int c = b + 1; c < k; ++c) {


                    for (int i = 0; i < suffixLen; ++i) {
                        d1[i] = c;
                        d2[i] = b;
                    }

                    d1[suffixLen] = b;
                    d2[suffixLen] = c;

                    for (int i = suffixLen + 1; i < L; ++i) {
                        d1[i] = pref[i - (suffixLen + 1)];
                        d2[i] = pref[i - (suffixLen + 1)];
                    }


                    int idx1 = 0, idx2 = 0, power = 1;
                    for (int i = 0; i < L; ++i) {
                        idx1 += d1[i] * power;
                        idx2 += d2[i] * power;
                        power *= k;
                    }


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

    int N = 1;
    for (int i = 0; i < L; ++i) N *= k;


    Graph g(N * layers);


    std::vector<int> d1(L), d2(L), pref;


    for (int z = 0; z < layers; ++z) {
        int offset = z * N;


        for (int j = 0; j < L; ++j) {
            int prefixLen = L - (j + 1);
            int suffixLen = j;


            int numPrefixes = 1;
            for (int t = 0; t < prefixLen; ++t) numPrefixes *= k;
            pref.resize(prefixLen);


            for (int p = 0; p < numPrefixes; ++p) {

                int tmp = p;
                for (int i = 0; i < prefixLen; ++i) {
                    pref[i] = tmp % k;
                    tmp /= k;
                }


                for (int b = 0; b < k; ++b) {
                    for (int c = b + 1; c < k; ++c) {

                        for (int i = 0; i < suffixLen; ++i) {
                            d1[i] = c;
                            d2[i] = b;
                        }

                        d1[suffixLen] = b;
                        d2[suffixLen] = c;

                        for (int i = suffixLen + 1; i < L; ++i) {
                            d1[i] = pref[i - (suffixLen + 1)];
                            d2[i] = pref[i - (suffixLen + 1)];
                        }


                        int idx1 = 0, idx2 = 0, power = 1;
                        for (int i = 0; i < L; ++i) {
                            idx1 += d1[i] * power;
                            idx2 += d2[i] * power;
                            power *= k;
                        }
                        idx1 += offset;
                        idx2 += offset;


                        g.AddEdge(idx1, idx2);
                    }
                }
            }
        }
    }

    int step = 1;
    for (int i = 0; i < L - 1; ++i) step *= k;
    int A = 0;
    int B = step - 1;
    int C = N - step;
    int D = N - 1;

    for (int z = 0; z + 1 < layers; ++z) {
        int off1 = z * N;
        int off2 = (z + 1) * N;


        int corners1[4] = {off1 + A, off1 + B, off1 + C, off1 + D};
        int corners2[4] = {off2 + A, off2 + B, off2 + C, off2 + D};
        for (int i = 0; i < 4; ++i) {
            g.AddEdge(corners1[i], corners2[i]);
        }
    }

    return g;
}

Graph GenerateRCR(int k, int r, int j) {

    int p = k + j;

    int M = 1 << p;

    int N = M * r;

    Graph g(N);


    for (int A = 0; A < M; ++A) {
        for (int b = 0; b < r; ++b) {
            int u = A * r + b;


            int v_ring = A * r + ((b + 1) % r);
            if (u < v_ring) {
                g.AddEdge(u, v_ring);
            }


            for (int bit = 0; bit < p; ++bit) {
                int A2 = A ^ (1 << bit);
                int v_cube = A2 * r + b;

                if (u < v_cube) {
                    g.AddEdge(u, v_cube);
                }
            }
        }
    }

    return g;
}