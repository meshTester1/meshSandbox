#pragma once
#include <vector>
#include <array>
#include <string>
#include <cstdint>

namespace rvd2foam {

    // ---------- Inputs ----------
    struct TetMeshView {
        // points: x0,y0,z0, x1,y1,z1, ...
        const double* points = nullptr;   // required
        std::size_t   numPoints = 0;      // number of 3D points (not *3)
        // tets: v0,v1,v2,v3 (0-based indices)
        const uint32_t* tets = nullptr;   // required
        std::size_t     numTets = 0;      // number of tetrahedra
    };

    struct SeedCloud {
        // seeds: x0,y0,z0, x1,y1,z1, ... (optional)
        const double* xyz = nullptr;
        std::size_t   count = 0;
    };

    // ---------- Output ----------
    struct OpenFOAMPolyMesh {
        std::vector<std::array<double, 3>> points;        // global points
        std::vector<std::vector<int>>     faces;         // polygonal faces (vertex indices)
        std::vector<int>                  owner;         // size = faces.size()
        std::vector<int>                  neighbour;     // size = #internal faces (prefix)
        // Single default boundary patch (we can split later)
        std::string boundaryName = "defaultFaces";
        std::string boundaryType = "patch";
    };

    // Build RVD polyhedra and convert to OpenFOAM lists (no external deps).
    // - If 'seeds' is empty, uses one seed per tet (tet centroid).
    // - Returns true on success. 'verbose' prints progress.
    bool convertTetToRVDPolyMesh(
        const TetMeshView& tet,
        const SeedCloud& seeds,
        OpenFOAMPolyMesh& outMesh,
        bool verbose = false);

    // Minimal ASCII writer for OpenFOAM polyMesh (use if you don't call your foamWriter).
    bool writeOpenFOAMPolyMesh(
        const OpenFOAMPolyMesh& m,
        const std::string& casePolyMeshDir,   // e.g. ".../constant/polyMesh"
        bool verbose = false);

    // One-shot convenience: compute + write polyMesh.
    bool convertTetToOpenFOAM(
        const TetMeshView& tet,
        const SeedCloud& seeds,
        const std::string& casePolyMeshDir,
        bool verbose = false);

} // namespace rvd2foam
