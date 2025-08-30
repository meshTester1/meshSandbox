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

#include "meshData.h"
#include "gmshReader.h"
#include "foamWriter.h"     // used by the legacy --mode=tet path
#include "convertRVD.h"     // new standalone RVD module (no geogram)

static std::string getArg(int argc, char** argv, int i, const std::string& dflt) {
    return (i < argc ? std::string(argv[i]) : dflt);
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

        // Pack tets into uint32[4*M]
        std::vector<uint32_t> tetIdx; tetIdx.reserve(gmshData.tets.size() * 4);
        for (const auto& t : gmshData.tets) {
            // assumes an indexable type (e.g., std::array<int,4>) with 0-based ids
            tetIdx.push_back(static_cast<uint32_t>(t[0]));
            tetIdx.push_back(static_cast<uint32_t>(t[1]));
            tetIdx.push_back(static_cast<uint32_t>(t[2]));
            tetIdx.push_back(static_cast<uint32_t>(t[3]));
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
