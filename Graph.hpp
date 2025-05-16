#pragma once

#include <vector>


class Graph {
public:
    explicit Graph(int n);

    void AddEdge(int node1, int node2, bool is_directed = false);


    int NumVertices() const;

    const std::vector<int> &Neighbors(int u) const;

    const std::vector<bool> &ComputeNodes() const;

    void SetCompute(const std::vector<bool> &compute) { is_compute = compute; }
private:
    int n_;
    std::vector<std::vector<int>> adj_;
    std::vector<bool> is_compute;
};
