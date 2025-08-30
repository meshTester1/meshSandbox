// gmshReader.cpp - Gmsh file reader implementation
#include "gmshReader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <map>
#include <array>
#include <algorithm>

namespace GmshReader {

    static bool starts_with(const std::string& s, const char* tag) {
        size_t n = std::strlen(tag);
        return s.size() >= n && std::equal(tag, tag + n, s.begin());
    }

    bool readGmshV1Ascii(const std::string& path, GmshData& out) {
        std::ifstream in(path);
        if (!in) {
            std::cerr << "Cannot open " << path << std::endl;
            return false;
        }

        std::string line;
        std::map<int, int> id2idx;

        auto read_nodes = [&](int n) -> bool {
            out.nodes.resize(n);
            id2idx.clear();
            for (int i = 0; i < n; ++i) {
                if (!std::getline(in, line)) {
                    std::cerr << "Failed to read node " << i << std::endl;
                    return false;
                }
                std::istringstream iss(line);
                int id;
                double x, y, z;
                if (!(iss >> id >> x >> y >> z)) {
                    std::cerr << "Failed to parse node line: " << line << std::endl;
                    return false;
                }
                int idx = (int)id2idx.size();
                id2idx[id] = idx;
                out.nodes[idx] = { x, y, z };
            }
            return true;
            };

        auto read_elements = [&](int m) -> bool {
            for (int i = 0; i < m; ++i) {
                if (!std::getline(in, line)) {
                    std::cerr << "Failed to read element " << i << std::endl;
                    return false;
                }
                std::istringstream iss(line);
                std::vector<long long> t;
                long long v;
                while (iss >> v) t.push_back(v);

                if (t.size() < 3) continue;
                int type = (int)t[1]; // element type

                // Triangle (type 2)
                if (type == 2 && t.size() >= 6) {
                    int a = (int)t[t.size() - 3], b = (int)t[t.size() - 2], c = (int)t[t.size() - 1];
                    Triangle T;
                    T.v[0] = id2idx.count(a) ? id2idx[a] : a - 1;
                    T.v[1] = id2idx.count(b) ? id2idx[b] : b - 1;
                    T.v[2] = id2idx.count(c) ? id2idx[c] : c - 1;

                    // Validate indices
                    if (T.v[0] >= 0 && T.v[0] < (int)out.nodes.size() &&
                        T.v[1] >= 0 && T.v[1] < (int)out.nodes.size() &&
                        T.v[2] >= 0 && T.v[2] < (int)out.nodes.size()) {
                        out.tris.push_back(T);
                    }
                }
                // Tetrahedron (type 4)
                else if (type == 4 && t.size() >= 7) {
                    int a = (int)t[t.size() - 4], b = (int)t[t.size() - 3],
                        c = (int)t[t.size() - 2], d = (int)t[t.size() - 1];
                    Tet T;
                    T.v[0] = id2idx.count(a) ? id2idx[a] : a - 1;
                    T.v[1] = id2idx.count(b) ? id2idx[b] : b - 1;
                    T.v[2] = id2idx.count(c) ? id2idx[c] : c - 1;
                    T.v[3] = id2idx.count(d) ? id2idx[d] : d - 1;

                    // Validate indices
                    if (T.v[0] >= 0 && T.v[0] < (int)out.nodes.size() &&
                        T.v[1] >= 0 && T.v[1] < (int)out.nodes.size() &&
                        T.v[2] >= 0 && T.v[2] < (int)out.nodes.size() &&
                        T.v[3] >= 0 && T.v[3] < (int)out.nodes.size()) {
                        out.tets.push_back(T);
                    }
                }
            }
            return true;
            };

        // Parse Gmsh v1 format only
        while (std::getline(in, line)) {
            if (starts_with(line, "$NOD")) {
                if (!std::getline(in, line)) {
                    std::cerr << "Failed to read node count" << std::endl;
                    return false;
                }
                int n = 0;
                { std::istringstream iss(line); iss >> n; }
                if (n <= 0) {
                    std::cerr << "Invalid node count: " << n << std::endl;
                    return false;
                }
                if (!read_nodes(n)) return false;
                if (!std::getline(in, line)) return false; // $ENDNOD
            }
            else if (starts_with(line, "$ELM")) {
                if (!std::getline(in, line)) {
                    std::cerr << "Failed to read element count" << std::endl;
                    return false;
                }
                int m = 0;
                { std::istringstream iss(line); iss >> m; }
                if (m <= 0) {
                    std::cerr << "Invalid element count: " << m << std::endl;
                    return false;
                }
                if (!read_elements(m)) return false;
                if (!std::getline(in, line)) return false; // $ENDELM
            }
        }

        if (out.nodes.empty() || out.tets.empty()) {
            std::cerr << "Parsed nodes=" << out.nodes.size()
                << ", tets=" << out.tets.size()
                << ", tris=" << out.tris.size() << std::endl;
            std::cerr << "ERROR: Need both nodes and tetrahedra for mesh conversion" << std::endl;
            return false;
        }

        // Derive boundary triangles if not present
        if (out.tris.empty()) {
            std::cout << "No boundary triangles found, deriving from tetrahedra..." << std::endl;
            auto key3 = [&](int a, int b, int c) {
                std::array<int, 3> s{ a, b, c };
                std::sort(s.begin(), s.end());
                return s;
                };
            std::map<std::array<int, 3>, int> cnt;
            auto add = [&](int a, int b, int c) { cnt[key3(a, b, c)]++; };

            for (const auto& T : out.tets) {
                add(T.v[0], T.v[1], T.v[2]);
                add(T.v[0], T.v[1], T.v[3]);
                add(T.v[0], T.v[2], T.v[3]);
                add(T.v[1], T.v[2], T.v[3]);
            }

            for (const auto& kv : cnt) {
                if (kv.second == 1) { // boundary face
                    out.tris.push_back(Triangle{ {kv.first[0], kv.first[1], kv.first[2]} });
                }
            }
            std::cout << "Derived " << out.tris.size() << " boundary triangles" << std::endl;
        }
        return true;
    }

    void convertTetsToPolyCells(const GmshData& g, std::vector<PolyCell>& cells) {
        cells.resize(g.tets.size());

        for (size_t t = 0; t < g.tets.size(); ++t) {
            const Tet& tet = g.tets[t];
            PolyCell& cell = cells[t];

            // Add the 4 vertices of the tetrahedron
            cell.verts.resize(4);
            cell.verts[0] = g.nodes[tet.v[0]];
            cell.verts[1] = g.nodes[tet.v[1]];
            cell.verts[2] = g.nodes[tet.v[2]];
            cell.verts[3] = g.nodes[tet.v[3]];

            // Add the 4 triangular faces of the tetrahedron
            // Face 0: vertices 0,1,2 (opposite to vertex 3)
            // Face 1: vertices 0,1,3 (opposite to vertex 2)
            // Face 2: vertices 0,2,3 (opposite to vertex 1)
            // Face 3: vertices 1,2,3 (opposite to vertex 0)
            cell.faces.resize(4);
            cell.faces[0] = { 0, 1, 2 }; // face opposite to vertex 3
            cell.faces[1] = { 0, 1, 3 }; // face opposite to vertex 2
            cell.faces[2] = { 0, 2, 3 }; // face opposite to vertex 1
            cell.faces[3] = { 1, 2, 3 }; // face opposite to vertex 0
        }

        std::cout << "Converted " << g.tets.size() << " tetrahedra to polyhedral cells" << std::endl;
    }

}