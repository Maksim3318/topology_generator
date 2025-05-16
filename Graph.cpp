#include "Graph.hpp"

Graph::Graph(int n) :
        n_(n), adj_(n), is_compute(n, true) {}

void Graph::AddEdge(int node1, int node2, bool is_directed) {
    if (node1 < n_ && node2 < n_) {
        adj_[node1].push_back(node2);
        if (!is_directed)
            adj_[node2].push_back(node1);
    }
}

int Graph::NumVertices() const {
    return n_;
}

const std::vector<int> &Graph::Neighbors(int u) const {
    return adj_[u];
}

const std::vector<bool> &Graph::ComputeNodes() const {
    return is_compute;
}
