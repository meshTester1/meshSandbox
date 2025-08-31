#pragma once
// rvd_geometry.h — math & clipping primitives (pure; no IO)

#include <cstdint>
#include <vector>

namespace rvd {

    struct V3 { double x, y, z; };
    struct Plane { V3 n; double c; };   // half-space: n·x <= c

    // vector ops
    V3     v_add(V3 a, V3 b);
    V3     v_sub(V3 a, V3 b);
    V3     v_scale(V3 a, double s);
    double v_dot(V3 a, V3 b);
    V3     v_cross(V3 a, V3 b);
    double v_norm(V3 a);

    // polygon helpers
    V3 newell_normal(const std::vector<V3>& P);
    V3 poly_centroid(const std::vector<V3>& P);

    // basic geometry
    V3    tet_centroid(const double* pts, const uint32_t* T);  // T -> 4 vertex ids
    V3    tri_normal(V3 a, V3 b, V3 c);
    Plane oriented_tri_plane_for_tet(V3 a, V3 b, V3 c, V3 tetCentroid);

    // Voronoi/Laguerre bisector: weights w0,w1 (set both 0 for unweighted Voronoi)
    Plane power_bisector(V3 s0, double w0, V3 s1, double w1);

    // convex polyhedron container
    struct Poly3 {
        std::vector<V3>               V;   // vertices
        std::vector<std::vector<int>> F;   // faces: CCW vertex loops
    };

    // make a tet poly from packed buffers
    Poly3 make_tet_poly(const double* pts, const uint32_t* T);

    // clip polygon by plane n·x <= c (Sutherland–Hodgman)
    std::vector<V3> clip_polygon(const std::vector<V3>& poly, const Plane& P, double eps);

    // clip convex polyhedron; optionally returns the intersection polygon with the plane
    Poly3 clip_polyhedron(const Poly3& Pin, const Plane& P, double eps, std::vector<V3>* cutPoly);

} // namespace rvd
