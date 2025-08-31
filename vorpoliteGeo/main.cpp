// main.cpp — Bulletproof implementation avoiding ALL Windows macro issues
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <cstdint>
#include <cmath>

#include "rvd_types.h"
#include "rvd_geometry.h"
#include "rvd_meshio_gmsh_v1.h"

// Completely custom utility functions to avoid any macro conflicts
namespace utils {

    std::string get_arg(int argc, char** argv, int index, const std::string& default_value) {
        if (index >= 0 && index < argc) {
            return std::string(argv[index]);
        }
        return default_value;
    }

    bool has_flag(int argc, char** argv, const std::string& flag) {
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == flag) {
                return true;
            }
        }
        return false;
    }

    double get_minimum(double a, double b) {
        return (a < b) ? a : b;
    }

    double get_maximum(double a, double b) {
        return (a > b) ? a : b;
    }
}

int main(int argc, char** argv) {
    // Get command line arguments
    std::string mesh_file = utils::get_arg(argc, argv, 1, "mediumCube.msh");
    bool verbose_output = !utils::has_flag(argc, argv, "--quiet");

    if (mesh_file.empty()) {
        std::cerr << "ERROR: No input mesh file specified\n";
        return 1;
    }

    // Load mesh from file
    rvd::TetMesh tetmesh;
    std::string error_msg;

    bool load_ok = rvd::load_gmsh_v1_ascii(mesh_file, tetmesh, &error_msg);
    if (!load_ok) {
        std::cerr << "ERROR: Failed to parse Gmsh v1 mesh: " << mesh_file << "\n";
        std::cerr << "       " << error_msg << "\n";
        return 1;
    }

    // Check if mesh is valid
    if (tetmesh.points.empty() || tetmesh.tets.empty()) {
        std::cerr << "ERROR: mesh has no points or no tetrahedra\n";
        return 1;
    }

    if (verbose_output) {
        std::cout << "[gmsh] points=" << tetmesh.points.size()
            << " tets=" << tetmesh.tets.size()
            << " tris=" << tetmesh.tris.size() << " (boundary candidates)\n";
    }

    // Convert to packed format
    std::size_t point_count = tetmesh.points.size();
    std::size_t tet_count = tetmesh.tets.size();

    // Pack points into flat array
    std::vector<double> point_data;
    point_data.resize(point_count * 3);

    for (std::size_t i = 0; i < point_count; ++i) {
        std::array<double, 3> point = tetmesh.points[i];
        point_data[3 * i + 0] = point[0];
        point_data[3 * i + 1] = point[1];
        point_data[3 * i + 2] = point[2];
    }

    // Pack tets into flat array
    std::vector<std::uint32_t> tet_data;
    tet_data.resize(tet_count * 4);

    for (std::size_t i = 0; i < tet_count; ++i) {
        std::array<std::uint32_t, 4> tet = tetmesh.tets[i];
        tet_data[4 * i + 0] = tet[0];
        tet_data[4 * i + 1] = tet[1];
        tet_data[4 * i + 2] = tet[2];
        tet_data[4 * i + 3] = tet[3];
    }

    // Create view
    rvd::TetMeshView mesh_view;
    mesh_view.points = point_data.data();
    mesh_view.numPoints = point_count;
    mesh_view.tets = tet_data.data();
    mesh_view.numTets = tet_count;

    std::cout << "Loaded tetrahedral mesh OK.\n";

    // Test geometry operations on first tet
    if (mesh_view.numTets > 0) {
        const std::uint32_t* first_tet = mesh_view.tets;

        // Create tet polyhedron
        rvd::Poly3 tet_poly = rvd::make_tet_poly(mesh_view.points, first_tet);
        rvd::V3 tet_center = rvd::tet_centroid(mesh_view.points, first_tet);

        // Find bounding box
        double min_x = 1e300;
        double max_x = -1e300;
        double min_y = 1e300;
        double max_y = -1e300;
        double min_z = 1e300;
        double max_z = -1e300;

        for (std::size_t i = 0; i < mesh_view.numPoints; ++i) {
            double x = mesh_view.points[3 * i + 0];
            double y = mesh_view.points[3 * i + 1];
            double z = mesh_view.points[3 * i + 2];

            min_x = utils::get_minimum(min_x, x);
            max_x = utils::get_maximum(max_x, x);
            min_y = utils::get_minimum(min_y, y);
            max_y = utils::get_maximum(max_y, y);
            min_z = utils::get_minimum(min_z, z);
            max_z = utils::get_maximum(max_z, z);
        }

        double size_x = max_x - min_x;
        double size_y = max_y - min_y;
        double size_z = max_z - min_z;
        double diagonal = std::sqrt(size_x * size_x + size_y * size_y + size_z * size_z);
        double epsilon = (diagonal > 0.0) ? (1e-9 * diagonal) : 1e-12;

        // Create second point for bisector
        rvd::V3 offset_center;
        offset_center.x = tet_center.x + 0.05 * diagonal;
        offset_center.y = tet_center.y;
        offset_center.z = tet_center.z;

        // Create bisector plane
        rvd::Plane bisector = rvd::power_bisector(tet_center, 0.0, offset_center, 0.0);

        // Clip the tetrahedron
        std::vector<rvd::V3> cut_vertices;
        rvd::Poly3 clipped_tet = rvd::clip_polyhedron(tet_poly, bisector, epsilon, &cut_vertices);

        // Show results
        if (clipped_tet.V.empty() || clipped_tet.F.empty()) {
            std::cout << "[smoke] bisector removed the entire tet (unexpected for tiny offset)\n";
        }
        else if (cut_vertices.size() >= 3) {
            std::cout << "[smoke] first tet clipped by bisector: cut polygon has "
                << cut_vertices.size() << " vertices\n";
        }
        else {
            std::cout << "[smoke] first tet clipped: no visible intersection polygon\n";
        }
    }

    return 0;
}