// main.cpp — Orchestrator with standalone RVD path (default) and legacy tet path.

#include <iostream>
#include <vector>
#include <string>
#include <cstdint>
#include <algorithm>
#include <cctype>

#include "meshData.h"
#include "gmshReader.h"
//#include "foamWriter.h"   // for --mode=tet

// ---------- helpers ----------
static std::string getArg(int argc, char** argv, int i, const std::string& dflt) {
    return (i < argc ? std::string(argv[i]) : dflt);
}

// Simple tolower function to avoid any compiler issues
static std::string toLowerSimple(const std::string& str) {
    std::string result = str;
    for (size_t i = 0; i < result.length(); ++i) {
        if (result[i] >= 'A' && result[i] <= 'Z') {
            result[i] = result[i] - 'A' + 'a';
        }
    }
    return result;
}

// Simple function to extract tet vertices - avoid complex templates
static void extractTetVertices(const Tet& t, uint32_t out[4]) {
    out[0] = static_cast<uint32_t>(t.v[0]);
    out[1] = static_cast<uint32_t>(t.v[1]);
    out[2] = static_cast<uint32_t>(t.v[2]);
    out[3] = static_cast<uint32_t>(t.v[3]);
}

int main(int argc, char** argv) {
    const std::string mshPath = getArg(argc, argv, 1, "mediumCube.msh");
    std::string outDir = getArg(argc, argv, 2, "constant/polyMesh");

    // Parse flags
    std::string mode = "rvd";
    bool verbose = true;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a.find("--mode=") == 0) {
            mode = a.substr(7); // skip "--mode="
            mode = toLowerSimple(mode);
        }
        else if (a == "--quiet") {
            verbose = false;
        }
    }

    std::cout << "Tetrahedral → OpenFOAM converter\n";
    std::cout << "Input     : " << mshPath << "\n";
    std::cout << "Output dir: " << outDir << "\n";
    std::cout << "Mode      : " << mode << "  (rvd|tet)\n\n";

    if (mshPath.empty()) {
        std::cerr << "ERROR: No input mesh file specified\n";
        return 1;
    }

    // 1) Read Gmsh v1 (ASCII)
    GmshData gmshData;
    if (!GmshReader::readGmshV1Ascii(mshPath, gmshData)) {
        std::cerr << "ERROR: Failed to parse Gmsh mesh: " << mshPath << "\n";
        return 1;
    }
    if (verbose) {
        std::cout << "[gmsh] nodes=" << gmshData.nodes.size()
            << " tets=" << gmshData.tets.size()
            << " tris=" << gmshData.tris.size() << " (boundary)\n";
    }

    if (gmshData.nodes.empty()) {
        std::cerr << "ERROR: No nodes found in mesh\n";
        return 1;
    }
    if (gmshData.tets.empty()) {
        std::cerr << "ERROR: No tetrahedra found in mesh\n";
        return 1;
    }

    // Branch A: RVD → OpenFOAM (default)
    if (mode == "rvd") {
        // Pack nodes into contiguous array
        std::vector<double> pts;
        pts.reserve(gmshData.nodes.size() * 3);
        for (size_t i = 0; i < gmshData.nodes.size(); ++i) {
            const Vec3& p = gmshData.nodes[i];
            pts.push_back(p.x);
            pts.push_back(p.y);
            pts.push_back(p.z);
        }

        // Pack tets into uint32 array with validation
        std::vector<uint32_t> tetIdx;
        tetIdx.reserve(gmshData.tets.size() * 4);

        for (size_t i = 0; i < gmshData.tets.size(); ++i) {
            const Tet& t = gmshData.tets[i];
            uint32_t v[4];
            extractTetVertices(t, v);

            // Validate indices
            bool validTet = true;
            for (int j = 0; j < 4; ++j) {
                if (v[j] >= gmshData.nodes.size()) {
                    if (verbose) {
                        std::cerr << "WARNING: Invalid tet vertex index: " << v[j]
                            << " >= " << gmshData.nodes.size() << "\n";
                    }
                    validTet = false;
                    break;
                }
            }

            if (validTet) {
                tetIdx.push_back(v[0]);
                tetIdx.push_back(v[1]);
                tetIdx.push_back(v[2]);
                tetIdx.push_back(v[3]);
            }
        }

        if (tetIdx.empty()) {
            std::cerr << "ERROR: No valid tetrahedra after validation\n";
            return 1;
        }

    // Branch B: Legacy tet → poly cells → foamWriter
    else if (mode == "tet") {
        std::vector<PolyCell> cells;

        GmshReader::convertTetsToPolyCells(gmshData, cells);

        // Count valid cells
        size_t valid = 0;
        for (size_t i = 0; i < cells.size(); ++i) {
            const PolyCell& c = cells[i];
            if (!c.faces.empty() && !c.verts.empty()) {
                ++valid;
            }
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
            std::cerr << "ERROR: Failed writing OpenFOAM polyMesh files to: " << outDir << "\n";
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