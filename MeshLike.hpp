#pragma once

#include "Graph.hpp"

Graph GenerateMesh(int rows, int cols, int n = 0) {
    if (n == 0)
        n = rows * cols;
    Graph g(n);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int current_vertex = i * cols + j;
            if (i != 0) {
                g.AddEdge(current_vertex, current_vertex - cols);
            }
            if (j != 0) {
                g.AddEdge(current_vertex, current_vertex - 1);
            }
        }
    }
    return g;
}

Graph GenerateTMesh(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    g.AddEdge(0, cols - 1);
    g.AddEdge(0, (rows - 1) * cols);
    g.AddEdge((rows - 1) * cols, rows * cols - 1);
    g.AddEdge(cols - 1, rows * cols - 1);
    return g;
}

Graph GenerateSDMesh(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int current_vertex = i * cols + j;
            if (i != rows - 1 && j != cols - 1) {
                g.AddEdge(current_vertex, current_vertex + cols + 1);
            }
        }
    }
    return g;
}

Graph GenerateZMesh(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int current_vertex = i * cols + j;
            if (i != rows - 1 && j != cols - 1 && i % 2 == 1) {
                g.AddEdge(current_vertex, current_vertex + cols + 1);
            }
            if (i != rows - 1 && j != 0 && i % 2 == 0) {
                g.AddEdge(current_vertex, current_vertex + cols - 1);
            }
        }
    }
    return g;
}

Graph GenerateDMesh(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            int current_vertex = i * cols + j;
            if (i != rows - 1 && j != cols - 1) {
                g.AddEdge(current_vertex, current_vertex + cols + 1);
            }
            if (i != rows - 1 && j != 0) {
                g.AddEdge(current_vertex, current_vertex + cols - 1);
            }
        }
    }
    return g;
}

Graph GenerateDCM(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    for (int i = 0; i < rows - 1; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (i % 2 == j % 2) {
                if (j != cols - 1) {
                    g.AddEdge(cols * i + j, cols * (i + 1) + j + 1);
                }
            } else if (j != 0) {
                g.AddEdge(cols * i + j, cols * (i + 1) + j - 1);
            }
        }
    }
    return g;
}

Graph GenerateMCC(int rows, int cols) {
    Graph g(rows * cols + cols / 2 * rows / 2 + ((cols - 1) / 2) * ((rows - 1) / 2));
    int crossbar_idx = rows * cols;
    std::vector<bool> nodes(crossbar_idx, true);
    for (int i = crossbar_idx; i < g.NumVertices(); i++) {
        nodes.push_back(false);
    }
    g.SetCompute(nodes);

    for (int i = 0; i < rows - 1; ++i) {
        for (int j = 0; j < cols - 1; ++j) {
            if (i % 2 == j % 2) {
                g.AddEdge(crossbar_idx, i * cols + j);
                g.AddEdge(crossbar_idx, (i + 1) * cols + j);
                g.AddEdge(crossbar_idx, i * cols + j + 1);
                g.AddEdge(crossbar_idx, (i + 1) * cols + j + 1);
                crossbar_idx++;
            }
        }
    }

    for (int i = 0; i < cols - 1; ++i) {
        if (i % 2 != 0) {
            g.AddEdge(i, i + 1);
        }
        if (i % 2 != rows % 2) {
            g.AddEdge((rows - 1) * cols + i, (rows - 1) * cols + i + 1);
        }
    }

    for (int i = 0; i < rows - 1; ++i) {
        if (i % 2 != 0) {
            g.AddEdge(i * cols, (i + 1) * cols);
        }
        if (i % 2 != cols % 2) {
            g.AddEdge((i + 1) * cols - 1, (i + 2) * cols - 1);
        }
    }
    return g;
}

Graph GenerateDiametricalMesh(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    for (int i = 0; i < rows - 1; i++) {
        for (int j = 0; j < cols; j++) {
            if (j < cols - 3) {
                g.AddEdge(i * cols + j, (i + 3) * cols + j + 3);
            }
            if (j > 2) {
                g.AddEdge(i * cols + j, (i + 3) * cols + j - 3);
            }
        }
    }
    return g;
}

Graph GenerateCBPMesh(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    for (int i = 0; i < rows - 1; i += 2) {
        for (int j = 0; j < cols; j += 2) {
            if (j < cols - 2) {
                g.AddEdge(i * cols + j, (i + 2) * cols + j + 2);
            }
            if (j > 1) {
                g.AddEdge(i * cols + j, (i + 2) * cols + j - 2);
            }
        }
    }
    return g;
}

Graph GenerateOCBPMesh(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    for (int i = 0; i < rows - 1; i++) {
        for (int j = 0; j < cols; j++) {
            if (j < cols - 2) {
                g.AddEdge(i * cols + j, (i + 2) * cols + j + 2);
            }
            if (j > 1) {
                g.AddEdge(i * cols + j, (i + 2) * cols + j - 2);
            }
        }
    }
    return g;
}

Graph GenerateC2Mesh(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    for (int i = 0; i < rows - 1; i += 2) {
        for (int j = 0; j < cols; j += 2) {
            if (((j % 4) ^ (i % 4)) != 0)
                continue;
            if (j < cols - 2) {
                g.AddEdge(i * cols + j, (i + 2) * cols + j + 2);
            }
            if (j > 1) {
                g.AddEdge(i * cols + j, (i + 2) * cols + j - 2);
            }
        }
    }
    return g;
}

Graph GenerateXDMesh(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    for (int i = 1; i < rows; i += 2) {
        for (int j = 1; j < cols; j += 2) {
            g.AddEdge(i * cols + j, (i - 1) * cols + j - 1);
            g.AddEdge(i * cols + j, (i - 1) * cols + j + 1);
            g.AddEdge(i * cols + j, (i + 1) * cols + j - 1);
            g.AddEdge(i * cols + j, (i + 1) * cols + j + 1);
        }
    }
    return g;
}

Graph GenerateXNetwork(int rows, int cols) {
    int n = rows * cols + (cols - 1) * (rows - 1);
    Graph g = GenerateMesh(rows, cols, n);
    int router_idx = rows * cols;
    std::vector<bool> nodes(rows * cols + (cols - 1) * (rows - 1), true);
    for (int i = router_idx; i < n; i++) {
        nodes[i] = false;
    }
    g.SetCompute(nodes);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (i != rows - 1 && j != cols - 1) { // se
                g.AddEdge(router_idx + i * (cols - 1) + j, i * cols + j);
            }
            if (i != 0 && j != cols - 1) { // ne
                g.AddEdge(router_idx + i * (cols - 1) + j - cols + 1, i * cols + j);
            }
            if (i != rows - 1 && j != 0) { // sw
                g.AddEdge(router_idx + i * (cols - 1) + j - 1, i * cols + j);
            }
            if (i != 0 && j != 0) { // nw
                g.AddEdge(router_idx + i * (cols - 1) + j - cols, i * cols + j);
            }
        }
    }
    return g;
}

Graph GenerateTorus(int rows, int cols) {
    Graph g = GenerateMesh(rows, cols);
    for (int i = 0; i < rows; ++i) {
        g.AddEdge(i * cols, i * cols + cols - 1);
    }
    for (int i = 0; i < cols; ++i) {
        g.AddEdge(i, (rows - 1) * cols + i);
    }
    return g;
}

Graph GenerateSDTorus(int rows, int cols) {
    Graph g = GenerateSDMesh(rows, cols);
    for (int i = 0; i < rows; ++i) {
        g.AddEdge(i * cols, i * cols + cols - 1);
    }
    for (int i = 0; i < cols; ++i) {
        g.AddEdge(i, (rows - 1) * cols + i);
    }
    return g;
}

Graph GenerateDTorus(int rows, int cols) {
    Graph g = GenerateDMesh(rows, cols);
    for (int i = 0; i < rows; ++i) {
        g.AddEdge(i * cols, i * cols + cols - 1);
    }
    for (int i = 0; i < cols; ++i) {
        g.AddEdge(i, (rows - 1) * cols + i);
    }
    return g;
}

Graph GenerateCBPTorus(int rows, int cols) {
    Graph g = GenerateCBPMesh(rows, cols);
    for (int i = 0; i < rows; ++i) {
        g.AddEdge(i * cols, i * cols + cols - 1);
    }
    for (int i = 0; i < cols; ++i) {
        g.AddEdge(i, (rows - 1) * cols + i);
    }
    return g;
}

Graph GenerateOCBPTorus(int rows, int cols) {
    Graph g = GenerateOCBPMesh(rows, cols);
    for (int i = 0; i < rows; ++i) {
        g.AddEdge(i * cols, i * cols + cols - 1);
    }
    for (int i = 0; i < cols; ++i) {
        g.AddEdge(i, (rows - 1) * cols + i);
    }
    return g;
}

Graph GenerateC2Torus(int rows, int cols) {
    Graph g = GenerateC2Mesh(rows, cols);
    for (int i = 0; i < rows; ++i) {
        g.AddEdge(i * cols, i * cols + cols - 1);
    }
    for (int i = 0; i < cols; ++i) {
        g.AddEdge(i, (rows - 1) * cols + i);
    }
    return g;
}

Graph GenerateXDTorus(int rows, int cols) {
    Graph g = GenerateXDMesh(rows, cols);
    for (int i = 0; i < rows; ++i) {
        g.AddEdge(i * cols, i * cols + cols - 1);
    }
    for (int i = 0; i < cols; ++i) {
        g.AddEdge(i, (rows - 1) * cols + i);
    }
    return g;
}

Graph Generate3DMesh(int x_dim, int y_dim, int z_dim) {
    int num_vertices = x_dim * y_dim * z_dim;
    Graph g(num_vertices);

    auto get_index = [&](int x, int y, int z) {
        return x * y_dim * z_dim
               + y * z_dim
               + z;
    };

    for (int x = 0; x < x_dim; ++x) {
        for (int y = 0; y < y_dim; ++y) {
            for (int z = 0; z < z_dim; ++z) {
                int current_index = get_index(x, y, z);

                int neighbor_x = x + 1;
                if (neighbor_x < x_dim) {
                    int neighbor_index = get_index(neighbor_x, y, z);
                    g.AddEdge(current_index, neighbor_index);
                }

                int neighbor_y = y + 1;
                if (neighbor_y < y_dim) {
                    int neighbor_index = get_index(x, neighbor_y, z);
                    g.AddEdge(current_index, neighbor_index);
                }

                int neighbor_z = z + 1;
                if (neighbor_z < z_dim) {
                    int neighbor_index = get_index(x, y, neighbor_z);
                    g.AddEdge(current_index, neighbor_index);
                }
            }
        }
    }

    return g;
}

Graph GenerateSMesh(int x_dim, int y_dim, int z_dim) {
    int num_vertices = x_dim * y_dim * z_dim;
    Graph graph(num_vertices);

    auto get_index = [&](int x, int y, int z) {
        return x * y_dim * z_dim
               + y * z_dim
               + z;
    };

    for (int x = 0; x < x_dim; ++x) {
        for (int y = 0; y < y_dim; ++y) {
            for (int z = 0; z < z_dim; ++z) {
                int current = get_index(x, y, z);

                int next_x = x + 1;
                if (next_x < x_dim) {
                    int neighbor = get_index(next_x, y, z);
                    graph.AddEdge(current, neighbor);
                }

                int next_y = y + 1;
                if (next_y < y_dim) {
                    int neighbor = get_index(x, next_y, z);
                    graph.AddEdge(current, neighbor);
                }
            }
        }
    }

    for (int x = 0; x < x_dim; ++x) {
        for (int y = 0; y < y_dim; ++y) {
            std::vector<int> plane_nodes;
            for (int z = 0; z < z_dim; ++z) {
                int node = get_index(x, y, z);
                plane_nodes.push_back(node);
            }
            int count = plane_nodes.size();
            for (int i = 0; i < count; ++i) {
                for (int j = i + 1; j < count; ++j) {
                    graph.AddEdge(plane_nodes[i], plane_nodes[j]);
                }
            }
        }
    }

    return graph;
}

Graph Generate3DTorus(int x_dim, int y_dim, int z_dim) {
    int num_vertices = x_dim * y_dim * z_dim;
    Graph g(num_vertices);

    auto get_index = [&](int x, int y, int z) {
        return x * y_dim * z_dim
               + y * z_dim
               + z;
    };

    for (int x = 0; x < x_dim; ++x) {
        for (int y = 0; y < y_dim; ++y) {
            for (int z = 0; z < z_dim; ++z) {
                int current_index = get_index(x, y, z);

                int next_x = (x + 1) % x_dim;
                int neighbor_x = get_index(next_x, y, z);
                g.AddEdge(current_index, neighbor_x);

                int next_y = (y + 1) % y_dim;
                int neighbor_y = get_index(x, next_y, z);
                g.AddEdge(current_index, neighbor_y);

                int next_z = (z + 1) % z_dim;
                int neighbor_z = get_index(x, y, next_z);
                g.AddEdge(current_index, neighbor_z);
            }
        }
    }

    return g;
}


Graph GenerateHypercube(int dimensions) {
    int num_vertices = 1 << dimensions;
    Graph g(num_vertices);

    for (int vertex = 0; vertex < num_vertices; ++vertex) {
        for (int bit = 0; bit < dimensions; ++bit) {
            int neighbor = vertex ^ (1 << bit);
            if (vertex < neighbor) {
                g.AddEdge(vertex, neighbor);
            }
        }
    }

    return g;
}

// Fibonacci-cube of rank n (n ≥ 1)
Graph GenerateFibonacciCube(int n) {
    if (n < 1 || n > 30) return Graph(0);

    const int full = 1 << n;
    std::vector<int> id(full, -1);
    std::vector<int> list;
    list.reserve(full);

    for (int v = 0; v < full; ++v) {
        if ((v & (v >> 1)) == 0) {
            id[v] = (int) list.size();
            list.push_back(v);
        }
    }
    int V = (int) list.size();

    Graph g(V);
    std::vector<bool> comp(V, true);

    for (int idx = 0; idx < V; ++idx) {
        int bits = list[idx];
        for (int b = 0; b < n; ++b) {
            int nbrBits = bits ^ (1 << b);
            int nbrId = id[nbrBits];
            if (nbrId > idx) g.AddEdge(idx, nbrId);
        }
    }
    g.SetCompute(comp);
    return g;
}


Graph GenerateHypermesh(const std::vector<int> &dimensions) {
    int dimension_count = dimensions.size();
    int num_vertices = 1;
    for (int size: dimensions) {
        num_vertices = num_vertices * size;
    }

    Graph g(num_vertices);

    std::vector<int> strides(dimension_count);
    strides[dimension_count - 1] = 1;
    for (int i = dimension_count - 2; i >= 0; --i) {
        strides[i] = strides[i + 1] * dimensions[i + 1];
    }

    for (int current_index = 0; current_index < num_vertices; ++current_index) {
        std::vector<int> coords(dimension_count);
        int remainder = current_index;
        for (int i = 0; i < dimension_count; ++i) {
            coords[i] = remainder / strides[i];
            remainder = remainder % strides[i];
        }

        for (int i = 0; i < dimension_count; ++i) {
            int next_coord = coords[i] + 1;
            if (next_coord < dimensions[i]) {
                int neighbor_index = current_index + strides[i];
                g.AddEdge(current_index, neighbor_index);
            }
        }
    }

    return g;
}

Graph GenerateHypertorus(const std::vector<int> &dimensions) {
    int dimension_count = dimensions.size();

    Graph graph = GenerateHypermesh(dimensions);

    std::vector<int> strides(dimension_count);
    strides[dimension_count - 1] = 1;
    for (int i = dimension_count - 2; i >= 0; --i) {
        strides[i] = strides[i + 1] * dimensions[i + 1];
    }

    int vertex_count = graph.NumVertices();

    for (int current_index = 0; current_index < vertex_count; ++current_index) {
        std::vector<int> coords(dimension_count);
        int remainder = current_index;
        for (int i = 0; i < dimension_count; ++i) {
            coords[i] = remainder / strides[i];
            remainder = remainder % strides[i];
        }

        for (int dim = 0; dim < dimension_count; ++dim) {
            bool is_last_layer = (coords[dim] + 1 == dimensions[dim]);
            if (is_last_layer) {
                std::vector<int> neighbor_coords = coords;
                neighbor_coords[dim] = 0;

                int neighbor_index = 0;
                for (int k = 0; k < dimension_count; ++k) {
                    neighbor_index += neighbor_coords[k] * strides[k];
                }

                graph.AddEdge(current_index, neighbor_index);
            }
        }
    }

    return graph;
}


Graph GenerateKAryNCubes(int size_per_dim,
                         int num_dimensions) {
    std::vector<int> dims;
    dims.reserve(num_dimensions);
    for (int i = 0; i < num_dimensions; ++i) {
        dims.push_back(size_per_dim);
    }

    return GenerateHypermesh(dims);
}


// Metacube MC(k,m): k — число бит class-ID, m — число бит node-ID
Graph GenerateMetacube(int k, int m) {
    if (k < 0 || m <= 0) {
        return Graph(0);
    }

    int num_classes = 1 << k;
    int total_bits = k + m * num_classes;
    int vertex_count = 1 << total_bits;

    Graph graph(vertex_count);

    int class_offset = m * num_classes;

    for (int vertex = 0; vertex < vertex_count; ++vertex) {
        int class_id = (vertex >> class_offset) & ((1 << k) - 1);
        int node_offset = class_id * m;

        for (int bit_index = 0; bit_index < m; ++bit_index) {
            int neighbor = vertex ^ (1 << (node_offset + bit_index));
            if (neighbor > vertex) {
                graph.AddEdge(vertex, neighbor);
            }
        }

        for (int bit_index = 0; bit_index < k; ++bit_index) {
            int neighbor = vertex ^ (1 << (class_offset + bit_index));
            if (neighbor > vertex) {
                graph.AddEdge(vertex, neighbor);
            }
        }
    }

    return graph;
}