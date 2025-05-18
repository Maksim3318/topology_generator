#pragma once

#include "Graph.hpp"

struct Delay {
    int router_to_node;
    int node_to_router;
    int router_to_router;
};

void CreateAnynetConfig(const std::string &filename, const Graph &g, Delay delay) {
    std::ofstream file(filename);
    if (!file.is_open())
        std::cerr << "Can't create file\n";

    auto compute = g.ComputeNodes();
    for (int i = 0; i < g.NumVertices(); i++)
    {
        file << "router " << i << " ";
        for (auto x: g.Neighbors(i)) {
            file << "router " << x << " ";
            if (delay.router_to_router != 1)
                file << delay.router_to_router << " ";
        }

        if (compute[i]) {
            file << "node " << i;
            if (delay.router_to_node != 1)
                file << " " << delay.router_to_node;

            if (delay.node_to_router != 1)
                file << "\nnode " << i << " router " << i << " " << delay.node_to_router;
        }

        file << "\n";
    }

    file.close();
    std::cout << "file " << filename << " successfully saved\n";
}
