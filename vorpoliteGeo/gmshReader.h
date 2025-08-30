// gmshReader.h - Gmsh file reader interface
#pragma once

#include "meshData.h"
#include <string>
#include <vector>

namespace GmshReader {

    // Read Gmsh v1 ASCII format
    bool readGmshV1Ascii(const std::string& path, GmshData& out);

    // Convert tetrahedral mesh to polyhedral cells
    void convertTetsToPolyCells(const GmshData& g, std::vector<PolyCell>& cells);

}