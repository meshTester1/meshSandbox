// main.cpp — Orchestrator with standalone RVD path (default) and legacy tet path.
// Usage:
//   vorpoliteGeo <mesh.msh> [outDir] [--mode=rvd|tet] [--quiet]
//
// Examples:
//   vorpoliteGeo mediumCube.msh
//   vorpoliteGeo block.msh case/constant/polyMesh --mode=rvd
//   vorpoliteGeo block.msh polyMesh --mode=tet

#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>
#include <type_traits>

#include "meshData.h"
#include "gmshReader.h"
#include "foamWriter.h"   // used by the legacy --mode=tet path
#include "convertRVD.h"   // new standalone RVD module (no geogram)

// ---------- small helpers ----------
static std::string getArg(int argc, char** argv, int i, const std::string& dflt) {
    return (i < argc ? std::string(argv[i]) : dflt);
}

// Detect common tet layouts in your repo so we can read 4 indices safely.
template <typename T, typename = void> struct has_index_op : std::false_type {};
template <typename T>
struct has_index_op<T, std::void_t<decltype(std::declval<const T&>()[0])>> : std::true_type {};

template <typename T, typename = void> struct has_v_array : std::false_type {};
template <typename T>
struct has_v_array<T, std::void_t<decltype(std::declval<const T&>().v[0])>> : std::true_type {};

template <typename T, typename = void> struct has_abcd : std::false_type {};
template <typename T>
struct has_abcd<T, std::void_t<
    decltype(std::declval<const T&>().a),
    decltype(std::declval<const T&>().b),
    decltype(std::declval<const T&>().c),
    decltype(std::declval<const T&>().d)>> : std::true_type {};

template <typename T, typename = void> struct has_i0123 : std::false_type {};
template <typename T>
struct has_i0123<T, std::void_t<
    decltype(std::declval<const T&>().i0),
    decltype(std::declval<const T&>().i1),
    decltype(std::declval<const T&>().i2),
    decltype(std::declval<const T&>().i3)>> : std::true_type {};

template <typename T, typename = void> struct has_n0123 : std::false_type {};
template <typename T>
struct has_n0123<T, std::void_t<
    decltype(std::declval<const T&>().n0),
    decltype(std::declval<const T&>().n1),
    decltype(std::declval<const T&>().n2),
    decltype(std::declval<const T&>().n3)>> : std::true_type {};

// Generic accessor: fills out[4] with 0-based vertex ids from a tet record
template <typename Tet>
inline void getTet4(const Tet& t, uint32_t out[4]) {
    if constexpr (has_index_op<Tet>::value) {
        out[0] = uint32_t(t[0]); out[1] = uint32_t(t[1]);
        out[2] = uint32_t(t[2]); out[3] = uint32_t(t[3]);
    }
    else if constexpr (has_v_array<Tet>::value) {
        out[0] = uint32_t(t.v[0]); out[1] = uint32_t(t.v[1]);
        out[2] = uint32_t(t.v[2]); out[3] = uint32_t(t.v[3]);
    }
    else if constexpr (has_abcd<Tet>::value) {
        out[0] = uint32_t(t.a); out[1] = uint32_t(t.b);
        out[2] = uint32_t(t.c); out[3] = uint32_t(t.d);
    }
    else if constexpr (has_i0123<Tet>::value) {
        out[0] = uint32_t(t.i0); out[1] = uint32_t(t.i1);
        out[2] = uint32_t(t.i2); out[3] = uint32_t(t.i3);
    }
    else if constexpr (has_n0123<Tet>::value) {
        out[0] = uint32_t(t.n0); out[1] = uint32_t(t.n1);
        out[2] = uint32_t(t.n2); out[3] = uint32_t(t.n3);
    }
    else {
        static_assert(sizeof(Tet) == 0, "Unsupported tet record type: add a branch to getTet4()");
    }
}

int main(int argc, char** argv) {
    const std::string mshPath = getArg(argc, argv, 1, "mediumCube.msh");
    std::string outDir = getArg(argc, argv, 2, "constant/polyMesh");

    // flags
    std::string mode = "rvd";   // "rvd" (default) or "tet"
    bool verbose = true;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a.rfind("--mode=", 0) == 0) {
            mode = a.substr(std::string("--mode=").size());
            std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);
        }
        else if (a == "--quiet") {
            verbose = false;
        }
    }

    std::cout << "Tetrahedral → OpenFOAM converter\n";
    std::cout << "Input     : " << mshPath << "\n";
    std::cout << "Output dir: " << outDir << "\n";
    std::cout << "Mode      : " << mode << "  (rvd|tet)\n\n";

    // 1) Read Gmsh v1 (ASCII) — your existing reader
    GmshData gmshData;
    if (!GmshReader::readGmshV1Ascii(mshPath, gmshData)) {
        std::cerr << "ERROR: Failed to parse Gmsh mesh\n";
        return 1;
    }
    if (verbose) {
        std::cout << "[gmsh] nodes=" << gmshData.nodes.size()
            << " tets=" << gmshData.tets.size()
            << " tris=" << gmshData.tris.size() << " (boundary)\n";
    }

    // Branch A: NEW — Standalone RVD → OpenFOAM (default)
    if (mode == "rvd") {
        // Pack nodes into contiguous double[3*N]
        std::vector<double> pts; pts.reserve(gmshData.nodes.size() * 3);
        for (const auto& p : gmshData.nodes) {
            // assumes a Vec3-like {x,y,z}; adjust if your type differs
            pts.push_back(p.x); pts.push_back(p.y); pts.push_back(p.z);
        }

        // Pack tets into uint32[4*M] via robust accessor
        std::vector<uint32_t> tetIdx; tetIdx.reserve(gmshData.tets.size() * 4);
        for (const auto& t : gmshData.tets) {
            uint32_t v[4]; getTet4(t, v);
            tetIdx.push_back(v[0]); tetIdx.push_back(v[1]);
            tetIdx.push_back(v[2]); tetIdx.push_back(v[3]);
        }

        // Build thin views for the RVD module
        rvd2foam::TetMeshView tv{
            pts.data(),
            gmshData.nodes.size(),
            tetIdx.data(),
            gmshData.tets.size()
        };

        // Seeds: empty ⇒ one seed per tet (tet centroid) inside the module
        rvd2foam::SeedCloud seeds{ nullptr, 0 };

        if (verbose) std::cout << "[rvd] converting …\n";
        const bool ok = rvd2foam::convertTetToOpenFOAM(tv, seeds, outDir, verbose);
        if (!ok) {
            std::cerr << "ERROR: RVD conversion failed.\n";
            return 2;
        }
        if (verbose) std::cout << "[rvd] done. Run `checkMesh` in the case directory.\n";
        return 0;
    }

    // Branch B: LEGACY — direct tet → poly cells → foamWriter (your existing path)
    if (mode == "tet") {
        std::vector<PolyCell> cells;
        GmshReader::convertTetsToPolyCells(gmshData, cells);

        // Basic sanity
        size_t valid = 0;
        for (const auto& c : cells) {
            if (!c.faces.empty() && !c.verts.empty()) ++valid;
        }
        if (valid == 0) {
            std::cerr << "ERROR: No valid polyhedral cells generated from tets.\n";
            return 3;
        }
        if (verbose) {
            std::cout << "[tet] valid cells: " << valid << " / " << cells.size() << "\n";
        }

        PolyMeshOut foamMesh;
        FoamWriter::buildOpenFoamTopology(cells, foamMesh);

        if (!FoamWriter::writePolyMesh(outDir, foamMesh)) {
            std::cerr << "ERROR: Failed writing OpenFOAM polyMesh files.\n";
            return 4;
        }

        if (verbose) {
            std::cout << "[tet] wrote polyMesh to: " << outDir << "\n";
            std::cout << "      points=" << foamMesh.points.size()
                << " faces=" << foamMesh.faces.size()
                << " (internal=" << foamMesh.nInternalFaces
                << ", boundary=" << (foamMesh.faces.size() - foamMesh.nInternalFaces)
                << ") cells=" << foamMesh.cellFaces.size() << "\n";
        }
        return 0;
    }

    std::cerr << "ERROR: Unknown mode '" << mode << "'. Use --mode=rvd or --mode=tet.\n";
    return 5;
}
