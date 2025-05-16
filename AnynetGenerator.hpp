#pragma once

#include "Graph.hpp"

void CreateAnynetConfig(const std::string &filename, const Graph &g) {
    std::ofstream file(filename);
    if (!file.is_open())
        std::cerr << "Can't create file\n";

    auto compute = g.ComputeNodes();
    for (int i = 0; i < g.NumVertices(); i++)
    {
        file << "router " << i << " ";
        for (auto x: g.Neighbors(i))
            file << "router " << x << " ";

        if (compute[i])
            file << "node " << i;

        file << "\n";
    }

    file.close();
    std::cout << "file " << filename << "successfully saved\n";
}
