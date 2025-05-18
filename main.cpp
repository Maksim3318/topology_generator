#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <functional>
#include <limits>

#include "Graph.hpp"
#include "Config.hpp"
#include "Metrics.hpp"
#include "GenerateGraphs.hpp"
#include "AnynetGenerator.hpp"

using namespace std;


void printUsage(const char *prog) {
    cerr << "Usage: " << prog << " <mode> <topology> [parameters...]\n"
         << "  mode: params | net\n"
         << "  [filename] - for net mode\n"
         << "  [delays router_to_node node_to_router router_to_router] - for net mode"
         << "  topology: one of [list of Generate functions]\n"
         << "Examples:\n"
         << "  " << prog << " params BFT 4 2\n"
         << "  " << prog << " net out.file 1 2 3 Circulant 1024 1 144 258 276\n";
}

// Generator function type and registry
using GeneratorFn = function<Graph(const vector<string> &)>;
struct GeneratorEntry {
    GeneratorFn fn;
    int minArgs;
    int maxArgs;
};
static unordered_map<string, GeneratorEntry> generators;

void InitGenerators() {
    // utility to split strings
    auto split = [](const string &s, char delim) {
        vector<string> res;
        string token;
        istringstream iss(s);
        while (getline(iss, token, delim)) { if (!token.empty()) res.push_back(token); }
        return res;
    };
    generators.clear();
    generators["3DDeBruijn"] = {[](const vector<string> &args) -> Graph {
        int d = stoi(args[0]);
        int k = stoi(args[1]);
        int layers = stoi(args[2]);
        return Generate3DDeBruijn(d, k, layers);
    }, 3, 3};
    generators["3DMesh"] = {[](const vector<string> &args) -> Graph {
        int x_dim = stoi(args[0]);
        int y_dim = stoi(args[1]);
        int z_dim = stoi(args[2]);
        return Generate3DMesh(x_dim, y_dim, z_dim);
    }, 3, 3};
    generators["3DTorus"] = {[](const vector<string> &args) -> Graph {
        int x_dim = stoi(args[0]);
        int y_dim = stoi(args[1]);
        int z_dim = stoi(args[2]);
        return Generate3DTorus(x_dim, y_dim, z_dim);
    }, 3, 3};
    generators["BFT"] = {[](const vector<string> &args) -> Graph {
        int radix = stoi(args[0]);
        int levels = stoi(args[1]);
        return GenerateBFT(radix, levels);
    }, 2, 2};
    generators["Banyan"] = {[](const vector<string> &args) -> Graph {
        int num_cores = stoi(args[0]);
        return GenerateBanyan(num_cores);
    }, 1, 1};
    generators["Benes"] = {[](const vector<string> &args) -> Graph {
        int num_terminals = stoi(args[0]);
        return GenerateBenes(num_terminals);
    }, 1, 1};
    generators["Butterfly"] = {[](const vector<string> &args) -> Graph {
        int num_cores = stoi(args[0]);
        return GenerateButterfly(num_cores);
    }, 1, 1};
    generators["C2Mesh"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateC2Mesh(rows, cols);
    }, 2, 2};
    generators["C2Torus"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateC2Torus(rows, cols);
    }, 2, 2};
    generators["CBPMesh"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateCBPMesh(rows, cols);
    }, 2, 2};
    generators["CBPTorus"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateCBPTorus(rows, cols);
    }, 2, 2};
    generators["CCN_HIN"] = {[](const vector<string> &args) -> Graph {
        vector<int> sizes;
        for (size_t i = 0; i < args.size(); ++i) sizes.push_back(stoi(args[i]));
        return GenerateCCN_HIN(sizes);
    }, 1, numeric_limits<int>::max()};
    generators["Circulant"] = {[](const vector<string> &args) -> Graph {
        int num_vertices = stoi(args[0]);
        vector<int> generators;
        for (size_t i = 1; i < args.size(); ++i) generators.push_back(stoi(args[i]));
        return GenerateCirculant(num_vertices, generators);
    }, 2, numeric_limits<int>::max()};
    generators["Clos"] = {[](const vector<string> &args) -> Graph {
        int num_terminals = stoi(args[0]);
        int radix = stoi(args[1]);
        return GenerateClos(num_terminals, radix);
    }, 2, 2};
    generators["Cnoc"] = {[](const vector<string> &args) -> Graph {
        int groups = stoi(args[0]);
        int terminals_per_group = stoi(args[1]);
        int radix = stoi(args[2]);
        return GenerateCnoc(groups, terminals_per_group, radix);
    }, 3, 3};
    generators["CubeConnectedCirculants"] = {[](const vector<string> &args) -> Graph {
        int hyperDim = stoi(args[0]);
        int circBase = stoi(args[1]);
        int circExp = stoi(args[2]);
        return GenerateCubeConnectedCirculants(hyperDim, circBase, circExp);
    }, 3, 3};
    generators["CubeConnectedCycles"] = {[](const vector<string> &args) -> Graph {
        int d = stoi(args[0]);
        return GenerateCubeConnectedCycles(d);
    }, 1, 1};
    generators["CubeTreeHybrid"] = {[](const vector<string> &args) -> Graph {
        int levels = stoi(args[0]);
        return GenerateCubeTreeHybrid(levels);
    }, 1, 1};
    generators["DCM"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateDCM(rows, cols);
    }, 2, 2};
    generators["DIMB"] = {[](const vector<string> &args) -> Graph {
        int d = stoi(args[0]);
        int k = stoi(args[1]);
        return GenerateDIMB(d, k);
    }, 2, 2};
    generators["DMesh"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateDMesh(rows, cols);
    }, 2, 2};
    generators["DTorus"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateDTorus(rows, cols);
    }, 2, 2};
    generators["DeBruijn"] = {[](const vector<string> &args) -> Graph {
        int d = stoi(args[0]);
        int k = stoi(args[1]);
        return GenerateDeBruijn(d, k);
    }, 2, 2};
    generators["Delta"] = {[](const vector<string> &args) -> Graph {
        int num_cores = stoi(args[0]);
        return GenerateDelta(num_cores);
    }, 1, 1};
    generators["DiagonalToroidalMesh"] = {[](const vector<string> &args) -> Graph {
        int k1 = stoi(args[0]);
        int k2 = stoi(args[1]);
        return GenerateDiagonalToroidalMesh(k1, k2);
    }, 2, 2};
    generators["DiametricalMesh"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateDiametricalMesh(rows, cols);
    }, 2, 2};
    generators["Dragonfly"] = {[](const vector<string> &args) -> Graph {
        int p = stoi(args[0]);
        int a = stoi(args[1]);
        int g = stoi(args[2]);
        return GenerateDragonfly(p, a, g);
    }, 3, 3};
    generators["FBFT"] = {[](const vector<string> &args) -> Graph {
        int radix_fly = stoi(args[0]);
        int stages_fly = stoi(args[1]);
        return GenerateFBFT(radix_fly, stages_fly);
    }, 2, 2};
    generators["FT"] = {[](const vector<string> &args) -> Graph {
        int radix = stoi(args[0]);
        int levels = stoi(args[1]);
        return GenerateFT(radix, levels);
    }, 2, 2};
    generators["FibonacciCube"] = {[](const vector<string> &args) -> Graph {
        int n = stoi(args[0]);
        return GenerateFibonacciCube(n);
    }, 1, 1};
    generators["FlattenedButterfly"] = {[](const vector<string> &args) -> Graph {
        int radix = stoi(args[0]);
        int layers = stoi(args[1]);
        return GenerateFlattenedButterfly(radix, layers);
    }, 2, 2};
    generators["Flip"] = {[](const vector<string> &args) -> Graph {
        int num_cores = stoi(args[0]);
        return GenerateFlip(num_cores);
    }, 1, 1};
    generators["GFT"] = {[](const vector<string> &args) -> Graph {
        int h = stoi(args[0]);
        int m = stoi(args[1]);
        int w = stoi(args[2]);
        return GenerateGFT(h, m, w);
    }, 3, 3};
    generators["Gaussian"] = {[](const vector<string> &args) -> Graph {
        int a = stoi(args[0]);
        int b = stoi(args[1]);
        return GenerateGaussian(a, b);
    }, 2, 2};
    generators["H3DMesh"] = {[](const vector<string> &args) -> Graph {
        int m = stoi(args[0]);
        int n = stoi(args[1]);
        int L = stoi(args[2]);
        return GenerateH3DMesh(m, n, L);
    }, 3, 3};
    generators["H3DTorus"] = {[](const vector<string> &args) -> Graph {
        int m = stoi(args[0]);
        int n = stoi(args[1]);
        int L = stoi(args[2]);
        return GenerateH3DTorus(m, n, L);
    }, 3, 3};
    generators["HERT"] = {[](const vector<string> &args) -> Graph {
        int k = stoi(args[0]);
        int m = stoi(args[1]);
        return GenerateHERT(k, m);
    }, 2, 2};
    generators["HSMBFT"] = {[](const vector<string> &args) -> Graph {
        int radix = stoi(args[0]);
        int depth = stoi(args[1]);
        return GenerateHSMBFT(radix, depth);
    }, 2, 2};
    generators["HTN"] = {[](const vector<string> &args) -> Graph {
        int m = stoi(args[0]);
        int n = stoi(args[1]);
        int L = stoi(args[2]);
        return GenerateHTN(m, n, L);
    }, 3, 3};
    generators["HexStarMesh"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        int d = stoi(args[2]);
        return GenerateHexStarMesh(rows, cols, d);
    }, 3, 3};
    generators["HierarchicalClique"] = {[](const vector<string> &args) -> Graph {
        int radix = stoi(args[0]);
        int levels = stoi(args[1]);
        return GenerateHierarchicalClique(radix, levels);
    }, 2, 2};
    generators["HoneycombMesh"] = {[](const vector<string> &args) -> Graph {
        int t = stoi(args[0]);
        return GenerateHoneycombMesh(t);
    }, 1, 1};
    generators["HoneycombTorus"] = {[](const vector<string> &args) -> Graph {
        int t = stoi(args[0]);
        return GenerateHoneycombTorus(t);
    }, 1, 1};
    generators["HyperCirculant"] = {[&split](const vector<string> &args) -> Graph {
        int num_vertices = stoi(args[0]);
        vector<vector<int>> genx;
        // args[1] should be group list separated by ';', e.g. "1,2;3,4"
        auto group_strs = split(args[1], ';');
        for (auto &gs: group_strs) {
            auto elems = split(gs, ',');
            vector<int> grp;
            for (auto &es: elems) grp.push_back(stoi(es));
            genx.push_back(grp);
        }
        return GenerateHyperCirculant(num_vertices, genx);
    }, 2, numeric_limits<int>::max()};
    generators["HyperX"] = {[](const vector<string> &args) -> Graph {
        vector<int> dims;
        for (size_t i = 0; i < args.size(); ++i) dims.push_back(stoi(args[i]));
        return GenerateHyperX(dims);
    }, 1, numeric_limits<int>::max()};
    generators["Hypercube"] = {[](const vector<string> &args) -> Graph {
        int dimensions = stoi(args[0]);
        return GenerateHypercube(dimensions);
    }, 1, 1};
    generators["Hypermesh"] = {[](const vector<string> &args) -> Graph {
        vector<int> dimensions;
        for (size_t i = 0; i < args.size(); ++i) dimensions.push_back(stoi(args[i]));
        return GenerateHypermesh(dimensions);
    }, 1, numeric_limits<int>::max()};
    generators["Hypertorus"] = {[](const vector<string> &args) -> Graph {
        vector<int> dimensions;
        for (size_t i = 0; i < args.size(); ++i) dimensions.push_back(stoi(args[i]));
        return GenerateHypertorus(dimensions);
    }, 1, numeric_limits<int>::max()};
    generators["KAryNCubes"] = {[](const vector<string> &args) -> Graph {
        int size_per_dim = stoi(args[0]);
        int num_dimensions = stoi(args[1]);
        return GenerateKAryNCubes(size_per_dim, num_dimensions);
    }, 2, 2};
    generators["L2Star"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateL2Star(rows, cols);
    }, 2, 2};
    generators["LNetwork"] = {[](const vector<string> &args) -> Graph {
        int a = stoi(args[0]);
        int b = stoi(args[1]);
        int c = stoi(args[2]);
        int d = stoi(args[3]);
        return GenerateLNetwork(a, b, c, d);
    }, 4, 4};
    generators["MCC"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateMCC(rows, cols);
    }, 2, 2};
    generators["MH3DT"] = {[](const vector<string> &args) -> Graph {
        int m = stoi(args[0]);
        int n = stoi(args[1]);
        int L = stoi(args[2]);
        int q = stoi(args[3]);
        return GenerateMH3DT(m, n, L, q);
    }, 4, 4};
    generators["MTESH"] = {[](const vector<string> &args) -> Graph {
        vector<int> sizes;
        for (size_t i = 0; i < args.size(); ++i) sizes.push_back(stoi(args[i]));
        return GenerateMTESH(sizes);
    }, 1, numeric_limits<int>::max()};
    generators["MTN"] = {[](const vector<string> &args) -> Graph {
        int m = stoi(args[0]);
        int level = stoi(args[1]);
        int q = stoi(args[2]);
        return GenerateMTN(m, level, q);
    }, 3, 3};
    generators["Mesh"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateMesh(rows, cols);
    }, 2, 2};
    generators["MeshMesh"] = {[](const vector<string> &args) -> Graph {
        int meshSide = stoi(args[0]);
        return GenerateMeshMesh(meshSide);
    }, 1, 1};
    generators["MeshOfSpidergon"] = {[](const vector<string> &args) -> Graph {
        int n = stoi(args[0]);
        int meshSide = stoi(args[1]);
        return GenerateMeshOfSpidergon(n, meshSide);
    }, 2, 2};
    generators["MeshOfStars"] = {[](const vector<string> &args) -> Graph {
        int leavesPerHub = stoi(args[0]);
        int meshSide = stoi(args[1]);
        return GenerateMeshOfStars(leavesPerHub, meshSide);
    }, 2, 2};
    generators["MeshOfTrees"] = {[](const vector<string> &args) -> Graph {
        int mesh_k = stoi(args[0]);
        return GenerateMeshOfTrees(mesh_k);
    }, 1, 1};
    generators["MeshRCTM"] = {[](const vector<string> &args) -> Graph {
        int meshSide = stoi(args[0]);
        return GenerateMeshRCTM(meshSide);
    }, 1, 1};
    generators["Metacube"] = {[](const vector<string> &args) -> Graph {
        int k = stoi(args[0]);
        int m = stoi(args[1]);
        return GenerateMetacube(k, m);
    }, 2, 2};
    generators["Midimew"] = {[](const vector<string> &args) -> Graph {
        int N = stoi(args[0]);
        int b = stoi(args[1]);
        return GenerateMidimew(N, b);
    }, 2, 2};
    generators["MultiplicativeCirculant"] = {[](const vector<string> &args) -> Graph {
        int base = stoi(args[0]);
        int degree = stoi(args[1]);
        return GenerateMultiplicativeCirculant(base, degree);
    }, 2, 2};
    generators["OCBPMesh"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateOCBPMesh(rows, cols);
    }, 2, 2};
    generators["OCBPTorus"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateOCBPTorus(rows, cols);
    }, 2, 2};
    generators["Omega"] = {[](const vector<string> &args) -> Graph {
        int num_cores = stoi(args[0]);
        return GenerateOmega(num_cores);
    }, 1, 1};
    generators["Paley"] = {[](const vector<string> &args) -> Graph {
        int num_vertices = stoi(args[0]);
        int degree = stoi(args[1]);
        return GeneratePaley(num_vertices, degree);
    }, 2, 2};
    generators["Pyramid"] = {[](const vector<string> &args) -> Graph {
        int radix = stoi(args[0]);
        int levels = stoi(args[1]);
        return GeneratePyramid(radix, levels);
    }, 2, 2};
    generators["RCR"] = {[](const vector<string> &args) -> Graph {
        int k = stoi(args[0]);
        int r = stoi(args[1]);
        int j = stoi(args[2]);
        return GenerateRCR(k, r, j);
    }, 3, 3};
    generators["RDT"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        int R = stoi(args[2]);
        int n = stoi(args[3]);
        return GenerateRDT(rows, cols, R, n);
    }, 4, 4};
    generators["RNT"] = {[](const vector<string> &args) -> Graph {
        int k = stoi(args[0]);
        int L = stoi(args[1]);
        int layers = stoi(args[2]);
        return GenerateRNT(k, L, layers);
    }, 3, 3};
    generators["Ricobit"] = {[](const vector<string> &args) -> Graph {
        int K = stoi(args[0]);
        return GenerateRicobit(K);
    }, 1, 1};
    generators["Ring"] = {[](const vector<string> &args) -> Graph {
        int num_vertices = stoi(args[0]);
        return GenerateRing(num_vertices);
    }, 1, 1};
    generators["RingMesh"] = {[](const vector<string> &args) -> Graph {
        int ringSize = stoi(args[0]);
        return GenerateRingMesh(ringSize);
    }, 1, 1};
    generators["RingRCTM"] = {[](const vector<string> &args) -> Graph {
        int ringSize = stoi(args[0]);
        return GenerateRingRCTM(ringSize);
    }, 1, 1};
    generators["SDMesh"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateSDMesh(rows, cols);
    }, 2, 2};
    generators["SDTorus"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateSDTorus(rows, cols);
    }, 2, 2};
    generators["SEM"] = {[](const vector<string> &args) -> Graph {
        int R = stoi(args[0]);
        return GenerateSEM(R);
    }, 1, 1};
    generators["SMBFT"] = {[](const vector<string> &args) -> Graph {
        int radix = stoi(args[0]);
        int depth = stoi(args[1]);
        return GenerateSMBFT(radix, depth);
    }, 2, 2};
    generators["SMesh"] = {[](const vector<string> &args) -> Graph {
        int x_dim = stoi(args[0]);
        int y_dim = stoi(args[1]);
        int z_dim = stoi(args[2]);
        return GenerateSMesh(x_dim, y_dim, z_dim);
    }, 3, 3};
    generators["Spidergon"] = {[](const vector<string> &args) -> Graph {
        int num_vertices = stoi(args[0]);
        return GenerateSpidergon(num_vertices);
    }, 1, 1};
    generators["StarMesh"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateStarMesh(rows, cols);
    }, 2, 2};
    generators["TESH"] = {[](const vector<string> &args) -> Graph {
        vector<int> sizes;
        for (size_t i = 0; i < args.size(); ++i) sizes.push_back(stoi(args[i]));
        return GenerateTESH(sizes);
    }, 1, numeric_limits<int>::max()};
    generators["THIN"] = {[](const vector<string> &args) -> Graph {
        int L = stoi(args[0]);
        return GenerateTHIN(L);
    }, 1, 1};
    generators["TMesh"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateTMesh(rows, cols);
    }, 2, 2};
    generators["Torus"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateTorus(rows, cols);
    }, 2, 2};
    generators["TwistedTorus"] = {[](const vector<string> &args) -> Graph {
        int width = stoi(args[0]);
        int height = stoi(args[1]);
        int twist = stoi(args[2]);
        return GenerateTwistedTorus(width, height, twist);
    }, 3, 3};
    generators["WKRecursive"] = {[](const vector<string> &args) -> Graph {
        int k = stoi(args[0]);
        int L = stoi(args[1]);
        return GenerateWKRecursive(k, L);
    }, 2, 2};
    generators["XDMesh"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateXDMesh(rows, cols);
    }, 2, 2};
    generators["XDTorus"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateXDTorus(rows, cols);
    }, 2, 2};
    generators["XNetwork"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateXNetwork(rows, cols);
    }, 2, 2};
    generators["XTree"] = {[](const vector<string> &args) -> Graph {
        int radix = stoi(args[0]);
        int levels = stoi(args[1]);
        return GenerateXTree(radix, levels);
    }, 2, 2};
    generators["ZMesh"] = {[](const vector<string> &args) -> Graph {
        int rows = stoi(args[0]);
        int cols = stoi(args[1]);
        return GenerateZMesh(rows, cols);
    }, 2, 2};
    generators["FatHTree"] = {[](const vector<string> &args) -> Graph {
        int n = stoi(args[0]);
        return GenerateFatHTree(n);
    }, 1, 1};
    generators["TreeMesh"] = {[](const vector<string> &args) -> Graph {
        int n = stoi(args[0]);
        return GenerateTreeMesh(n);
    }, 1, 1};
}

int main(int argc, char *argv[]) {
    Config cfg;
    if (!LoadConfig("config.conf", cfg)) {
        cerr << "Failed to load config.conf" << endl;
        return 1;
    }
    InitGenerators();
    if (argc < 3) {
        printUsage(argv[0]);
        return 1;
    }
    string mode = argv[1];

    if (mode == "net") {
        if (argc < 4) {
            printUsage(argv[0]);
            return 1;
        }
        string filename = argv[2];
        Delay delay;
        delay.router_to_node = atoi(argv[3]);
        delay.node_to_router = atoi(argv[4]);
        delay.router_to_router = atoi(argv[5]);
        string topo_net = argv[6];
        auto it_net = generators.find(topo_net);
        if (it_net == generators.end()) {
            cerr << "Unknown topology: " << topo_net << endl;
            printUsage(argv[0]);
            return 1;
        }
        vector<string> args_net;
        for (int i = 7; i < argc; ++i)
            args_net.emplace_back(argv[i]);
        Graph g_net = it_net->second.fn(args_net);
        if (g_net.NumVertices() == 0) {
            cerr << "Error: generated empty graph" << endl;
            return 1;
        }
        CreateAnynetConfig(filename, g_net, delay);
        return 0;
    } else if (mode == "params") {
        string topo = argv[2];
        auto it = generators.find(topo);
        if (it == generators.end()) {
            cerr << "Unknown topology: " << topo << endl;
            printUsage(argv[0]);
            return 1;
        }
        vector<string> args;
        for (int i = 3; i < argc; ++i) args.emplace_back(argv[i]);
        int argcnt = args.size();
        int minA = it->second.minArgs;
        int maxA = it->second.maxArgs;
        if (argcnt < minA || argcnt > maxA) {
            printUsage(argv[0]);
            return 1;
        }
        Graph g = it->second.fn(args);


        if (g.NumVertices() == 0) {
            cerr << "Graph is empty, validate parameters" << endl;
            return 1;
        }
        Metrics metrics;
        metrics = ComputeMetrics(g, cfg);
        cout << "Nodes count: " << metrics.nodes_count << endl;
        cout << "Compute nodes count: " << metrics.compute_nodes_count << endl;
        cout << "Diameter: " << metrics.diameter << endl;
        cout << "Avg path length: " << metrics.average_path_length << endl;
        cout << "Min degree: " << metrics.min_degree << endl;
        cout << "Max degree: " << metrics.max_degree << endl;
        cout << "Bisection width: " << metrics.bisection_width << endl;
        cout << "Global Packing density: " << metrics.global_packing_density << endl;
        cout << "Norm local Packing density: " << metrics.norm_local_packing_density << endl;
        cout << "Avg #shortest paths: " << metrics.average_num_shortest_paths << endl;
        cout << "Norm Avg #shortest paths: " << metrics.norm_average_num_shortest_paths << endl;
    } else {
        cerr << "Unknown mode: " << mode << endl;
        return 1;
    }

    return 0;
}