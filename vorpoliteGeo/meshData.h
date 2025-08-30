// meshData.h - All shared data structures for the mesh converter
#pragma once

#include <vector>
#include <cmath>

// ---------- tiny math / hashing ----------
struct Vec3 {
    double x, y, z;
};

struct Vec3Key {
    long long xi, yi, zi;
    bool operator==(const Vec3Key& o) const {
        return xi == o.xi && yi == o.yi && zi == o.zi;
    }
};

struct Vec3KeyHash {
    size_t operator()(const Vec3Key& k) const {
        size_t h = 1469598103934665603ull;
        auto mix = [&](long long v) { h ^= (size_t)v; h *= 1099511628211ull; };
        mix(k.xi); mix(k.yi); mix(k.zi);
        return h;
    }
};

static inline Vec3Key quantize(const Vec3& p, double q = 1e-10) {
    return Vec3Key{
        (long long)std::llround(p.x / q),
        (long long)std::llround(p.y / q),
        (long long)std::llround(p.z / q)
    };
}

// ---------- Mesh Elements ----------
struct Triangle {
    int v[3];
};

struct Tet {
    int v[4];
};

struct GmshData {
    std::vector<Vec3> nodes;
    std::vector<Triangle> tris; // optional boundary
    std::vector<Tet> tets;      // volume tets (required)
};

// ---------- Per-cell polyhedron (skeleton for future polyhedral support) ----------
struct PolyCell {
    std::vector<Vec3> verts;
    std::vector<std::vector<int>> faces; // indices into verts
};

struct PolyMeshOut {
    std::vector<Vec3> points;
    std::vector<std::vector<int>> faces;     // global point indices
    std::vector<int> owner;                  // faces.size()
    std::vector<int> neighbour;              // nInternalFaces
    std::vector<std::vector<int>> cellFaces; // per cell -> list of global face ids
    int nInternalFaces = 0;
};