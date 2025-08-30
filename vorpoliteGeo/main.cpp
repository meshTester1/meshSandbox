// main.cpp - Main orchestrator for tetrahedral to OpenFOAM converter
#include <iostream>
#include <vector>
#include <string>

#include "meshData.h"
#include "gmshReader.h"
#include "foamWriter.h"

int main(int argc, char** argv) {
    const std::string mshPath = (argc >= 2 ? argv[1] : "mediumCube.msh");

    std::cout << "Tetrahedral mesh to OpenFOAM polyMesh converter" << std::endl;
    std::cout << "Input file: " << mshPath << std::endl;

    // 1) Read Gmsh v1 format only
    GmshData gmshData;
    if (!GmshReader::readGmshV1Ascii(mshPath, gmshData)) {
        std::cerr << "Failed to parse Gmsh mesh" << std::endl;
        return 1;
    }
    std::cout << "Read: " << gmshData.nodes.size() << " nodes, " << gmshData.tets.size()
        << " tets, " << gmshData.tris.size() << " boundary tris" << std::endl;

    // 2) Convert tetrahedra directly to polyhedral cells
    std::vector<PolyCell> cells;
    GmshReader::convertTetsToPolyCells(gmshData, cells);

    // 3) Validate cells
    size_t valid = 0;
    for (const auto& c : cells) {
        if (!c.faces.empty() && !c.verts.empty()) ++valid;
    }

    if (valid == 0) {
        std::cerr << "ERROR: No valid polyhedral cells generated" << std::endl;
        std::cerr << "Check that input mesh contains valid tetrahedra" << std::endl;
        return 1;
    }
    std::cout << "Valid cells: " << valid << " / " << cells.size() << std::endl;

    // 4) Build OpenFOAM topology
    PolyMeshOut foamMesh;
    FoamWriter::buildOpenFoamTopology(cells, foamMesh);

    // 5) Write polyMesh files
    const std::string polyDir = "polyMesh";
    if (!FoamWriter::writePolyMesh(polyDir, foamMesh)) {
        std::cerr << "Failed writing OpenFOAM polyMesh files" << std::endl;
        return 1;
    }

    std::cout << "Tetrahedral mesh to OpenFOAM conversion completed successfully" << std::endl;
    std::cout << "Output polyMesh: points=" << foamMesh.points.size()
        << ", faces=" << foamMesh.faces.size()
        << " (internal=" << foamMesh.nInternalFaces
        << ", boundary=" << (foamMesh.faces.size() - foamMesh.nInternalFaces)
        << "), cells=" << foamMesh.cellFaces.size() << std::endl;

    return 0;
}