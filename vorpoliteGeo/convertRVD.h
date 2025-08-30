#pragma once

#include <vector>
#include <array>
#include <string>
#include <cstdint>

namespace rvd2foam {

    // Input tetrahedral mesh view (shallow wrapper around raw arrays)
    struct TetMeshView {
        const double* points;      // xyz coordinates: points[3*i+0], points[3*i+1], points[3*i+2]
        std::size_t numPoints;
        const std::uint32_t* tets; // tetrahedra: tets[4*j+0], tets[4*j+1], tets[4*j+2], tets[4*j+3] (0-based indices)
        std::size_t numTets;
    };

    // Optional seed points for RVD (if nullptr, uses tet centroids)
    struct SeedCloud {
        const double* xyz;         // xyz coordinates: xyz[3*i+0], xyz[3*i+1], xyz[3*i+2]
        std::size_t count;
    };

    // Output OpenFOAM polyMesh structure
    struct OpenFOAMPolyMesh {
        std::vector<std::array<double, 3>> points;      // vertex coordinates
        std::vector<std::vector<int>> faces;            // faces as vertex indices
        std::vector<int> owner;                         // face owner cell
        std::vector<int> neighbour;                     // face neighbour cell (internal faces only)

        // Boundary patch info (single patch by default)
        std::string boundaryName = "walls";
        std::string boundaryType = "wall";
    };

    // Main conversion function: tetrahedral mesh → RVD polyhedral cells
    bool convertTetToRVDPolyMesh(
        const TetMeshView& tet,
        const SeedCloud& seeds,
        OpenFOAMPolyMesh& out,
        bool verbose = false
    );

    // Write OpenFOAM polyMesh files (points, faces, owner, neighbour, boundary)
    bool writeOpenFOAMPolyMesh(
        const OpenFOAMPolyMesh& mesh,
        const std::string& polyMeshDir,
        bool verbose = false
    );

    // High-level convenience function: tet mesh → OpenFOAM files
    bool convertTetToOpenFOAM(
        const TetMeshView& tet,
        const SeedCloud& seeds,
        const std::string& casePolyMeshDir,
        bool verbose = false
    );

} // namespace rvd2foam