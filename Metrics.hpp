#pragma once

#include "Graph.hpp"
#include "Config.hpp"

struct Metrics {
    int nodes_count;
    int compute_nodes_count;
    int diameter;
    double average_path_length;
    int min_degree;
    int max_degree;
    int bisection_width;
    double global_packing_density;
    double norm_local_packing_density;
    double average_num_shortest_paths;
    double norm_average_num_shortest_paths;
};

Metrics ComputeMetrics(const Graph &g, const Config &cfg);


