// main.cpp — Clean C++17 implementation avoiding all compiler ambiguities
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <cstring>

#include "rvd_types.h"
#include "rvd_geometry.h"
#include "rvd_meshio_gmsh_v1.h"

namespace {
    // Helper functions in anonymous namespace
    std::string getArgument(int argc, char** argv, int index, const std::string& defaultValue) {
        if (index < argc) {
            return std::string(argv[index]);
        }
        return defaultValue;
    }

    bool hasFlag(int argc, char** argv, const std::string& flag) {
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == flag) {
                return true;
            }
        }
        return false;
    }

    // Explicit min/max functions to avoid macro issues
    template<typename T>
    constexpr const T& minimum(const T& a, const T& b) {
        return (a < b) ? a : b;
    }

    template<typename T>
    constexpr const T& maximum(const T& a, const T& b) {
        return (a > b) ? a : b;
    }
}

int main(int argc, char** argv) {
    // Parse command line arguments
    const std::string meshPath = getArgument(argc, argv, 1, "mediumCube.msh");
    const bool verbose = !hasFlag(argc, argv, "--quiet");

    if (meshPath.empty()) {
        std::cerr << "ERROR: No input mesh file specified\n";
        return 1;
    }

    // Step 1: Load the tetrahedral mesh from Gmsh file
    rvd::TetMesh mesh;
    std::string errorMessage;

    const bool loadSuccess = rvd::load_gmsh_v1_ascii(meshPath, mesh, &errorMessage);
    if (!loadSuccess) {
        std::cerr << "ERROR: Failed to parse Gmsh v1 mesh: " << meshPath << "\n"
            << "       " << errorMessage << "\n";
        return 1;
    }

    // Validate mesh data
    if (mesh.points.empty() || mesh.tets.empty()) {
        std::cerr << "ERROR: mesh has no points or no tetrahedra\n";
        return 1;
    }

    if (verbose) {
        std::cout << "[gmsh] points=" << mesh.points.size()
            << " tets=" << mesh.tets.size()
            << " tris=" << mesh.tris.size() << " (boundary candidates)\n";
    }

    // Step 2: Convert to packed arrays for processing
    const std::size_t numPoints = mesh.points.size();
    const std::size_t numTets = mesh.tets.size();

    // Create packed point array (x,y,z,x,y,z,...)
    std::vector<double> packedPoints;
    packedPoints.reserve(numPoints * 3);

    for (std::size_t i = 0; i < numPoints; ++i) {
        const std::array<double, 3>& point = mesh.points[i];
        packedPoints.push_back(point[0]); // x
        packedPoints.push_back(point[1]); // y  
        packedPoints.push_back(point[2]); // z
    }

    // Create packed tetrahedron index array
    std::vector<std::uint32_t> packedTets;
    packedTets.reserve(numTets * 4);

    for (std::size_t t = 0; t < numTets; ++t) {
        const std::array<std::uint32_t, 4>& tet = mesh.tets[t];
        packedTets.push_back(tet[0]);
        packedTets.push_back(tet[1]);
        packedTets.push_back(tet[2]);
        packedTets.push_back(tet[3]);
    }

    // Create mesh view
    rvd::TetMeshView meshView;
    meshView.points = packedPoints.data();
    meshView.numPoints = numPoints;
    meshView.tets = packedTets.data();
    meshView.numTets = numTets;

    std::cout << "Loaded tetrahedral mesh OK.\n";

    // Step 3: Geometry smoke test - clip first tetrahedron
    if (meshView.numTets > 0) {
        const std::uint32_t* firstTet = meshView.tets;

        // Create polyhedron from first tetrahedron
        rvd::Poly3 tetrahedronPoly = rvd::make_tet_poly(meshView.points, firstTet);
        rvd::V3 tetCentroid = rvd::tet_centroid(meshView.points, firstTet);

        // Calculate bounding box for epsilon calculation
        double minX = 1e300, minY = 1e300, minZ = 1e300;
        double maxX = -1e300, maxY = -1e300, maxZ = -1e300;

        for (std::size_t i = 0; i < meshView.numPoints; ++i) {
            const double x = meshView.points[3 * i + 0];
            const double y = meshView.points[3 * i + 1];
            const double z = meshView.points[3 * i + 2];

            minX = minimum(minX, x);
            maxX = maximum(maxX, x);
            minY = minimum(minY, y);
            maxY = maximum(maxY, y);
            minZ = minimum(minZ, z);
            maxZ = maximum(maxZ, z);
        }

        const double deltaX = maxX - minX;
        const double deltaY = maxY - minY;
        const double deltaZ = maxZ - minZ;
        const double diagonal = std::sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
        const double epsilon = (diagonal > 0.0) ? (1e-9 * diagonal) : 1e-12;

        // Create offset centroid for bisector plane
        rvd::V3 offsetCentroid;
        offsetCentroid.x = tetCentroid.x + 0.05 * diagonal;
        offsetCentroid.y = tetCentroid.y;
        offsetCentroid.z = tetCentroid.z;

        // Create power bisector plane
        rvd::Plane bisectorPlane = rvd::power_bisector(tetCentroid, 0.0, offsetCentroid, 0.0);

        // Perform clipping
        std::vector<rvd::V3> cutPolygon;
        rvd::Poly3 clippedPoly = rvd::clip_polyhedron(tetrahedronPoly, bisectorPlane, epsilon, &cutPolygon);

        // Report results
        if (clippedPoly.V.empty() || clippedPoly.F.empty()) {
            std::cout << "[smoke] bisector removed the entire tet (unexpected for tiny offset)\n";
        }
        else if (cutPolygon.size() >= 3) {
            std::cout << "[smoke] first tet clipped by bisector: cut polygon has "
                << cutPolygon.size() << " vertices\n";
        }
        else {
            std::cout << "[smoke] first tet clipped: no visible intersection polygon\n";
        }
    }

    return 0;
}