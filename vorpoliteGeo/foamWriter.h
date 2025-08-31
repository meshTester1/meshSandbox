//// foamWriter.h - OpenFOAM mesh writer interface
//#pragma once
//
//#include "meshData.h"
//#include <string>
//#include <vector>
//
//namespace FoamWriter {
//
//    // Build OpenFOAM topology from polyhedral cells
//    void buildOpenFoamTopology(const std::vector<PolyCell>& cells, PolyMeshOut& out);
//
//    // Write complete polyMesh to directory
//    bool writePolyMesh(const std::string& dir, const PolyMeshOut& mesh);
//
//}