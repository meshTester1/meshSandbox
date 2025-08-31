#pragma once
// rvd_types.h — shared PODs & keys for the Voronoi/RVD pipeline (no OpenFOAM)

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace rvd {

    // Thin view (pointer-based) over a tet mesh: packed xyz, packed 4-ids.
    struct TetMeshView {
        const double* points = nullptr;   // length = 3 * numPoints
        std::size_t     numPoints = 0;
        const uint32_t* tets = nullptr;   // length = 4 * numTets (0-based)
        std::size_t     numTets = 0;
    };

    // Concrete container for a tet mesh loaded from file.
    struct TetMesh {
        std::vector<std::array<double, 3>>  points; // xyz
        std::vector<std::array<uint32_t, 4>> tets;  // 0-based indices
        std::vector<std::array<uint32_t, 3>> tris;  // optional boundary tris (0-based)
    };

    // Input sites (Voronoi/Laguerre). If xyz==nullptr, we’ll synthesize one per tet later.
    struct SeedCloud {
        const double* xyz = nullptr;  // length = 3*count (or nullptr)
        const double* w = nullptr;  // optional weights (nullptr => zeros)
        std::size_t   count = 0;
    };

    // Internal face record (placeholder for later steps).
    struct FaceRec {
        std::vector<int> verts;  // vertex ids into a global point array
        int owner = -1;
        int neighbour = -1;      // -1 => boundary
    };

    // Sorted triangle key (for outer-surface detection / dedup).
    struct TriKey {
        uint32_t a, b, c; // sorted
        bool operator==(const TriKey& o) const noexcept { return a == o.a && b == o.b && c == o.c; }
    };
    struct TriKeyHash {
        std::size_t operator()(const TriKey& k) const noexcept {
            std::size_t h = 1469598103934665603ull;
            auto mix = [&](uint32_t v) { h ^= (std::size_t)v + 0x9e3779b97f4a7c15ull; h *= 1099511628211ull; };
            mix(k.a); mix(k.b); mix(k.c);
            return h;
        }
    };

} // namespace rvd
