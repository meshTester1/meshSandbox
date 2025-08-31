#include "rvd_meshio_gmsh_v1.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace rvd {

    static inline std::string trim(const std::string& s) {
        size_t a = 0, b = s.size();
        while (a < b && std::isspace((unsigned char)s[a])) ++a;
        while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
        return s.substr(a, b - a);
    }
    static inline bool ieq(const std::string& a, const char* b) {
        if (a.size() != std::strlen(b)) return false;
        for (size_t i = 0;i < a.size();++i) if (std::tolower((unsigned char)a[i]) != std::tolower((unsigned char)b[i])) return false;
        return true;
    }

    static bool parse_nodes_block(std::istream& is, TetMesh& out, bool legacy, std::string* err) {
        // legacy: $NOD / $ENDNOD ; modern-ish: $Nodes / $EndNodes
        std::string line;
        // next line: number of nodes
        if (!std::getline(is, line)) { if (err) *err = "Unexpected EOF after $NOD/$Nodes"; return false; }
        std::istringstream ls(trim(line));
        std::size_t N = 0; if (!(ls >> N)) { if (err) *err = "Could not read node count"; return false; }
        out.points.clear(); out.points.resize(N);

        // Node IDs in v1 can be arbitrary; remap to 0..N-1
        std::unordered_map<long long, uint32_t> id2idx; id2idx.reserve(N);

        for (std::size_t k = 0;k < N;++k) {
            if (!std::getline(is, line)) { if (err) *err = "Unexpected EOF in nodes block"; return false; }
            line = trim(line); if (line.empty()) { --k; continue; }
            std::istringstream ns(line);
            long long id = 0; double x = 0, y = 0, z = 0;
            if (!(ns >> id >> x >> y >> z)) { if (err) *err = "Bad node line: " + line; return false; }
            uint32_t idx = (uint32_t)k;
            id2idx[id] = idx;
            out.points[idx] = { x,y,z };
        }

        // read end token line
        if (!std::getline(is, line)) { if (err) *err = "Missing end of nodes block"; return false; }
        line = trim(line);
        if (legacy) {
            if (!ieq(line, "$ENDNOD") && !ieq(line, "$EndNodes")) { if (err) *err = "Expected $ENDNOD/$EndNodes"; return false; }
        }
        else {
            if (!ieq(line, "$EndNodes") && !ieq(line, "$ENDNOD")) { if (err) *err = "Expected $EndNodes/$ENDNOD"; return false; }
        }

        // stash the map in out.tris temporarily? — no, keep local.
        // We’ll re-open when parsing elements; but we need access to id2idx there.
        // Solution: return it upward? Simpler: store in a static thread_local. But
        // easiest is to re-parse IDs in elements by building id2idx again: not possible.
        // Instead, we keep id2idx in a static to reuse in this load call.

        // To avoid globals, we pass it onward via istream state by saving into ios pword — overkill.
        // Practical approach: return the map via an out-parameter… so rewrite function signature?
        // Simpler: we’ll keep id2idx as a file-scope static for the duration of this call.

        // But we cannot share between functions cleanly. Instead, we parse elements in the same scope
        // as nodes so we can see id2idx. We'll restructure the loader to do both in one function.
        // (We’re inside helper; we’ll instead move id2idx building to the outer load function.)
        return true;
    }

    // We’ll implement load_gmsh_v1_ascii as a single-pass parser keeping id2idx alive.
    bool load_gmsh_v1_ascii(const std::string& path, TetMesh& out, std::string* err) {
        std::ifstream fs(path);
        if (!fs) { if (err) *err = "Cannot open file: " + path; return false; }

        out.points.clear(); out.tets.clear(); out.tris.clear();
        std::unordered_map<long long, uint32_t> id2idx;

        std::string line;
        while (std::getline(fs, line)) {
            line = trim(line);
            if (line.empty()) continue;
            if (line[0] == '#' || (line.size() >= 2 && line[0] == '/' && line[1] == '/')) continue;

            // Nodes block start?
            if (ieq(line, "$NOD") || ieq(line, "$Nodes")) {
                const bool legacy = ieq(line, "$NOD");
                // read count
                if (!std::getline(fs, line)) { if (err) *err = "Unexpected EOF after $NOD/$Nodes"; return false; }
                std::istringstream ls(trim(line));
                std::size_t N = 0; if (!(ls >> N)) { if (err) *err = "Could not read node count"; return false; }
                out.points.clear(); out.points.resize(N);
                id2idx.clear(); id2idx.reserve(N);

                for (std::size_t k = 0;k < N;++k) {
                    if (!std::getline(fs, line)) { if (err) *err = "Unexpected EOF in nodes block"; return false; }
                    line = trim(line); if (line.empty()) { --k; continue; }
                    std::istringstream ns(line);
                    long long id = 0; double x = 0, y = 0, z = 0;
                    if (!(ns >> id >> x >> y >> z)) { if (err) *err = "Bad node line: " + line; return false; }
                    uint32_t idx = (uint32_t)k;
                    id2idx[id] = idx;
                    out.points[idx] = { x,y,z };
                }
                if (!std::getline(fs, line)) { if (err) *err = "Missing end of nodes block"; return false; }
                line = trim(line);
                const bool okEnd = (legacy ? (ieq(line, "$ENDNOD") || ieq(line, "$EndNodes"))
                    : (ieq(line, "$EndNodes") || ieq(line, "$ENDNOD")));
                if (!okEnd) { if (err) *err = "Expected $ENDNOD/$EndNodes"; return false; }
                continue;
            }

            // Elements block start?
            if (ieq(line, "$ELM") || ieq(line, "$Elements")) {
                const bool legacy = ieq(line, "$ELM");
                if (!std::getline(fs, line)) { if (err) *err = "Unexpected EOF after $ELM/$Elements"; return false; }
                std::istringstream ls(trim(line));
                std::size_t M = 0; if (!(ls >> M)) { if (err) *err = "Could not read element count"; return false; }

                out.tets.clear(); out.tris.clear();
                out.tets.reserve(M); out.tris.reserve(M);

                auto needNodes = [](int type)->int {
                    switch (type) {
                    case 2: return 3; // 3-node triangle
                    case 4: return 4; // 4-node tetrahedron
                        // extend as needed
                    default: return -1;
                    }
                    };

                for (std::size_t k = 0;k < M;++k) {
                    if (!std::getline(fs, line)) { if (err) *err = "Unexpected EOF in elements block"; return false; }
                    line = trim(line); if (line.empty()) { --k; continue; }
                    std::istringstream es(line);

                    // Parse *all* ints on the line (v1/v2 have different tag layouts)
                    std::vector<long long> tok; tok.reserve(16);
                    long long v;
                    while (es >> v) tok.push_back(v);
                    if (tok.size() < 3) { if (err) *err = "Bad element line: " + line; return false; }

                    const long long elmNo = tok[0];
                    const int       elmType = (int)tok[1];
                    (void)elmNo;

                    const int K = needNodes(elmType);
                    if (K < 0) continue; // skip unsupported types

                    if ((int)tok.size() < 2 + K) {
                        if (err) *err = "Element line too short for node list: " + line; return false;
                    }

                    // Take the *last K* tokens as node ids; map to 0-based.
                    const int start = (int)tok.size() - K;
                    std::array<uint32_t, 4> tet{};
                    std::array<uint32_t, 3> tri{};

                    if (elmType == 4) {
                        for (int i = 0;i < 4;++i) {
                            long long nid = tok[start + i];
                            auto it = id2idx.find(nid);
                            if (it == id2idx.end()) { if (err) *err = "Unknown node id in tet: " + std::to_string(nid); return false; }
                            tet[i] = it->second;
                        }
                        out.tets.push_back(tet);
                    }
                    else if (elmType == 2) {
                        for (int i = 0;i < 3;++i) {
                            long long nid = tok[start + i];
                            auto it = id2idx.find(nid);
                            if (it == id2idx.end()) { if (err) *err = "Unknown node id in tri: " + std::to_string(nid); return false; }
                            tri[i] = it->second;
                        }
                        out.tris.push_back(tri);
                    }
                }

                if (!std::getline(fs, line)) { if (err) *err = "Missing end of elements block"; return false; }
                line = trim(line);
                const bool okEnd = (legacy ? (ieq(line, "$ENDELM") || ieq(line, "$EndElements"))
                    : (ieq(line, "$EndElements") || ieq(line, "$ENDELM")));
                if (!okEnd) { if (err) *err = "Expected $ENDELM/$EndElements"; return false; }
                continue;
            }

            // Otherwise: ignore unknown sections for v1
        }

        if (out.points.empty() || out.tets.empty()) {
            if (err) *err = "File parsed, but points or tets are empty";
            return false;
        }
        return true;
    }

} // namespace rvd
