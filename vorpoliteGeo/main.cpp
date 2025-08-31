// main.cpp — Minimal RVD harness using new modules (no foamWriter, no legacy deps)
// Build: C++17
// Usage: vorpoliteGeo <mesh.msh> [--quiet]
//
// This loads a Gmsh v1 ASCII tetra mesh via rvd_meshio_gmsh_v1,
// packs it into a TetMeshView, prints counts, and runs a tiny clip
// smoke-test on the first tet using rvd_geometry (modules 1–2).

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>

#include "rvd_types.h"
#include "rvd_geometry.h"
#include "rvd_meshio_gmsh_v1.h"

// ---------- helpers ----------
static std::string getArg(int argc, char** argv, int i, const std::string& dflt) {
    return (i < argc ? std::string(argv[i]) : dflt);
}

static bool hasFlag(int argc, char** argv, const std::string& flag) {
    for (int i = 1; i < argc; ++i) if (flag == argv[i]) return true;
    return false;
}

int main(int argc, char** argv) {
    const std::string mshPath = getArg(argc, argv, 1, "mediumCube.msh");
    const bool verbose = !hasFlag(argc, argv, "--quiet");

    if (mshPath.empty()) {
        std::cerr << "ERROR: No input mesh file specified\n";
        return 1;
    }

    // 1) Read Gmsh v1 ASCII tetra mesh (modules: rvd_types + rvd_meshio_gmsh_v1)
    rvd::TetMesh mesh;
    std::string err;
    if (!rvd::load_gmsh_v1_ascii(mshPath, mesh, &err)) {
        std::cerr << "ERROR: Failed to parse Gmsh v1 mesh: " << mshPath << "\n"
            << "       " << err << "\n";
        return 1;
    }

    if (verbose) {
        std::cout << "[gmsh] points=" << mesh.points.size()
            << " tets=" << mesh.tets.size()
            << " tris=" << mesh.tris.size() << " (boundary candidates)\n";
    }

    if (mesh.points.empty() || mesh.tets.empty()) {
        std::cerr << "ERROR: mesh has no points or no tetrahedra\n";
        return 1;
    }

    // 2) Pack into contiguous arrays for a TetMeshView (future passes use this)
    std::vector<double> pts;
    pts.resize(mesh.points.size() * 3);
    for (size_t i = 0; i < mesh.points.size(); ++i) {
        pts[3 * i + 0] = mesh.points[i][0];
        pts[3 * i + 1] = mesh.points[i][1];
        pts[3 * i + 2] = mesh.points[i][2];
    }

    std::vector<uint32_t> tets;
    tets.resize(mesh.tets.size() * 4);
    for (size_t e = 0; e < mesh.tets.size(); ++e) {
        tets[4 * e + 0] = mesh.tets[e][0];
        tets[4 * e + 1] = mesh.tets[e][1];
        tets[4 * e + 2] = mesh.tets[e][2];
        tets[4 * e + 3] = mesh.tets[e][3];
    }

    rvd::TetMeshView tv{ pts.data(), pts.size() / 3, tets.data(), tets.size() / 4 };

    std::cout << "Loaded tetrahedral mesh OK.\n";

    // 3) Quick geometry smoke-test (modules: rvd_geometry) — optional but handy now
    //    We clip the first tet by a simple Voronoi bisector and report the cut polygon size.
    if (tv.numTets > 0) {
        const uint32_t* T0 = tv.tets + 0;
        // Build the tet as a convex poly and compute its centroid
        rvd::Poly3 poly = rvd::make_tet_poly(tv.points, T0);
        rvd::V3 c0 = rvd::tet_centroid(tv.points, T0);

        // Create a second "site" slightly offset from c0 to define a bisector plane
        // Scale offset by bbox diagonal for stability
        double xmin = +1e300, ymin = +1e300, zmin = +1e300;
        double xmax = -1e300, ymax = -1e300, zmax = -1e300;
        for (size_t i = 0; i < tv.numPoints; ++i) {
            xmin = (std::min)(xmin, tv.points[3 * i + 0]);
            ymin = (std::min)(ymin, tv.points[3 * i + 1]);
            zmin = (std::min)(zmin, tv.points[3 * i + 2]);
            xmax = (std::max)(xmax, tv.points[3 * i + 0]);
            ymax = (std::max)(ymax, tv.points[3 * i + 1]);
            zmax = (std::max)(zmax, tv.points[3 * i + 2]);
        }

        const double dx = xmax - xmin, dy = ymax - ymin, dz = zmax - zmin;
        const double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double eps = (diag > 0 ? 1e-9 * diag : 1e-12);

        rvd::V3 c1{ c0.x + 0.05 * diag, c0.y, c0.z }; // slight +x offset
        rvd::Plane H = rvd::power_bisector(c0, 0.0, c1, 0.0);

        std::vector<rvd::V3> cutPoly;
        rvd::Poly3 clipped = rvd::clip_polyhedron(poly, H, eps, &cutPoly);

        if (clipped.V.empty() || clipped.F.empty()) {
            std::cout << "[smoke] bisector removed the entire tet (unexpected for tiny offset)\n";
        }
        else if (cutPoly.size() >= 3) {
            std::cout << "[smoke] first tet clipped by bisector: cut polygon has "
                << cutPoly.size() << " vertices\n";
        }
        else {
            std::cout << "[smoke] first tet clipped: no visible intersection polygon\n";
        }
    }

    // That's it for modules 1–2 + loader. Next steps will plug in seeds/exterior/power/builder.
    return 0;
}
