#include "convertRVD.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>  
#include <cstdlib>   
#include <iostream>

namespace rvd2foam {

    // =========================== small geometry ===========================

    struct V3 { double x, y, z; };

    static inline V3 v_add(V3 a, V3 b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
    static inline V3 v_sub(V3 a, V3 b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
    static inline V3 v_scale(V3 a, double s) { return { a.x * s, a.y * s, a.z * s }; }
    static inline double v_dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
    static inline V3 v_cross(V3 a, V3 b) {
        return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
    }
    static inline double v_norm(V3 a) { return std::sqrt(v_dot(a, a)); }

    struct Plane { V3 n; double c; };

    static inline double signed_dist(const Plane& P, V3 x) {
        return v_dot(P.n, x) - P.c;
    }

    static inline V3 lerp(V3 a, V3 b, double t) {
        return { a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, a.z + (b.z - a.z) * t };
    }

    static inline V3 tri_normal(V3 a, V3 b, V3 c) {
        return v_cross(v_sub(b, a), v_sub(c, a));
    }

    static inline V3 poly_centroid(const std::vector<V3>& poly) {
        V3 c = { 0,0,0 };
        if (poly.empty()) return c;
        for (size_t i = 0; i < poly.size(); ++i) {
            c = v_add(c, poly[i]);
        }
        return v_scale(c, 1.0 / double(poly.size()));
    }

    static inline V3 newell_normal(const std::vector<V3>& P) {
        double nx = 0, ny = 0, nz = 0;
        size_t n = P.size();
        if (n < 3) return { 0,0,0 };
        for (size_t i = 0; i < n; i++) {
            const V3& a = P[i];
            const V3& b = P[(i + 1) % n];
            nx += (a.y - b.y) * (a.z + b.z);
            ny += (a.z - b.z) * (a.x + b.x);
            nz += (a.x - b.x) * (a.y + b.y);
        }
        return { nx,ny,nz };
    }

    // =========================== polyhedron struct ===========================

    struct Poly3 {
        std::vector<V3> V;
        std::vector<std::vector<int>> F;
    };

    static Poly3 make_tet_poly(const V3& p0, const V3& p1, const V3& p2, const V3& p3) {
        Poly3 P;
        P.V.push_back(p0);
        P.V.push_back(p1);
        P.V.push_back(p2);
        P.V.push_back(p3);

        std::vector<int> f0, f1, f2, f3;
        f0.push_back(0); f0.push_back(2); f0.push_back(1);
        f1.push_back(0); f1.push_back(1); f1.push_back(3);
        f2.push_back(1); f2.push_back(2); f2.push_back(3);
        f3.push_back(2); f3.push_back(0); f3.push_back(3);

        P.F.push_back(f0);
        P.F.push_back(f1);
        P.F.push_back(f2);
        P.F.push_back(f3);
        return P;
    }

    // Simple vertex key using only basic types
    struct VKey {
        long long x, y, z;
        bool operator<(const VKey& other) const {
            if (x != other.x) return x < other.x;
            if (y != other.y) return y < other.y;
            return z < other.z;
        }
    };

    static VKey make_vertex_key(const V3& a, double eps) {
        const double s = 1.0 / eps;
        VKey key;
        key.x = static_cast<long long>(std::llround(a.x * s));
        key.y = static_cast<long long>(std::llround(a.y * s));
        key.z = static_cast<long long>(std::llround(a.z * s));
        return key;
    }

    static bool point_inside_plane(const V3& q, const Plane& P, double eps) {
        return signed_dist(P, q) <= eps;
    }

    static bool face_is_degenerate(const std::vector<int>& f, const std::vector<V3>& vertices) {
        if (f.size() < 3) return true;
        V3 a = vertices[f[0]];
        V3 b = vertices[f[1]];
        V3 c = vertices[f[2]];
        return v_norm(tri_normal(a, b, c)) < 1e-18;
    }

    // Clip a polygon by plane n·x <= c
    static std::vector<V3> clip_polygon(const std::vector<V3>& poly, const Plane& P, double eps) {
        std::vector<V3> out;
        if (poly.empty()) return out;

        V3 S = poly.back();
        bool Sin = point_inside_plane(S, P, eps);

        for (size_t i = 0; i < poly.size(); ++i) {
            const V3& E = poly[i];
            bool Ein = point_inside_plane(E, P, eps);
            if (Ein) {
                if (!Sin) {
                    V3 d = v_sub(E, S);
                    double denom = v_dot(P.n, d);
                    double t = (P.c - v_dot(P.n, S)) / (denom == 0 ? 1e-30 : denom);
                    t = std::max(0.0, std::min(1.0, t));
                    out.push_back(lerp(S, E, t));
                }
                out.push_back(E);
            }
            else if (Sin) {
                V3 d = v_sub(E, S);
                double denom = v_dot(P.n, d);
                double t = (P.c - v_dot(P.n, S)) / (denom == 0 ? 1e-30 : denom);
                t = std::max(0.0, std::min(1.0, t));
                out.push_back(lerp(S, E, t));
            }
            S = E; Sin = Ein;
        }

        std::vector<V3> comp;
        for (size_t i = 0; i < out.size(); ++i) {
            const V3& q = out[i];
            if (comp.empty() || v_norm(v_sub(q, comp.back())) > eps) {
                comp.push_back(q);
            }
        }
        if (comp.size() >= 3 && v_norm(v_sub(comp.front(), comp.back())) <= eps) {
            comp.back() = comp.front();
            comp.pop_back();
        }
        return comp;
    }

    static Poly3 clip_polyhedron(const Poly3& Pin, const Plane& P, double eps, std::vector<V3>* outCutPolygon = nullptr) {
        Poly3 Q;
        std::vector<V3> cutLoopAccum;

        for (size_t fi = 0; fi < Pin.F.size(); ++fi) {
            const std::vector<int>& face = Pin.F[fi];
            std::vector<V3> poly;
            poly.reserve(face.size());
            for (size_t i = 0; i < face.size(); ++i) {
                poly.push_back(Pin.V[face[i]]);
            }

            std::vector<V3> clipped = clip_polygon(poly, P, eps);
            if (clipped.size() >= 3) {
                int base = static_cast<int>(Q.V.size());
                for (size_t i = 0; i < clipped.size(); ++i) {
                    Q.V.push_back(clipped[i]);
                }
                std::vector<int> idx;
                for (size_t k = 0; k < clipped.size(); ++k) {
                    idx.push_back(base + static_cast<int>(k));
                }
                Q.F.push_back(idx);
            }
            for (size_t i = 0; i < clipped.size(); ++i) {
                if (std::fabs(signed_dist(P, clipped[i])) <= 2 * eps) {
                    cutLoopAccum.push_back(clipped[i]);
                }
            }
        }

        // Weld vertices using simple map
        std::map<VKey, int> vmap;
        std::vector<V3> WV;
        WV.reserve(Q.V.size());
        std::vector<int> remap(Q.V.size(), -1);

        for (size_t i = 0; i < Q.V.size(); ++i) {
            VKey K = make_vertex_key(Q.V[i], eps);
            std::map<VKey, int>::iterator it = vmap.find(K);
            if (it == vmap.end()) {
                int nid = static_cast<int>(WV.size());
                vmap[K] = nid;
                WV.push_back(Q.V[i]);
                remap[i] = nid;
            }
            else {
                remap[i] = it->second;
            }
        }

        for (size_t fi = 0; fi < Q.F.size(); ++fi) {
            std::vector<int>& f = Q.F[fi];
            for (size_t i = 0; i < f.size(); ++i) {
                f[i] = remap[f[i]];
            }

            std::vector<int> g;
            g.reserve(f.size());
            for (size_t k = 0; k < f.size(); ++k) {
                if (g.empty() ||
                    WV[g.back()].x != WV[f[k]].x ||
                    WV[g.back()].y != WV[f[k]].y ||
                    WV[g.back()].z != WV[f[k]].z) {
                    g.push_back(f[k]);
                }
            }
            f = g;
        }
        Q.V = WV;

        std::vector<std::vector<int>> Fclean;
        Fclean.reserve(Q.F.size());
        for (size_t i = 0; i < Q.F.size(); ++i) {
            if (!face_is_degenerate(Q.F[i], Q.V)) {
                Fclean.push_back(Q.F[i]);
            }
        }
        Q.F = Fclean;

        // Build cut polygon if requested
        if (outCutPolygon) {
            std::vector<V3> Cw;
            std::set<VKey> seen;
            for (size_t i = 0; i < cutLoopAccum.size(); ++i) {
                VKey K = make_vertex_key(cutLoopAccum[i], eps);
                if (seen.find(K) == seen.end()) {
                    seen.insert(K);
                    Cw.push_back(cutLoopAccum[i]);
                }
            }
            if (Cw.size() >= 3) {
                V3 n = P.n;
                double nn = v_norm(n);
                if (nn < 1e-30) {
                    n.x = 0; n.y = 0; n.z = 1;
                }
                else {
                    n = v_scale(n, 1.0 / nn);
                }

                V3 tmp = (std::fabs(n.x) < 0.9) ? V3{ 1,0,0 } : V3{ 0,1,0 };
                V3 u_raw = v_sub(tmp, v_scale(n, v_dot(n, tmp)));
                double ulen = v_norm(u_raw);
                V3 u = (ulen < 1e-30) ? V3{ 1,0,0 } : v_scale(u_raw, 1.0 / ulen);
                V3 v = v_cross(n, u);

                V3 C = { 0,0,0 };
                for (size_t i = 0; i < Cw.size(); ++i) {
                    C = v_add(C, Cw[i]);
                }
                C = v_scale(C, 1.0 / double(Cw.size()));

                struct Node { V3 q; double ang; };
                std::vector<Node> nodes;
                nodes.reserve(Cw.size());
                for (size_t i = 0; i < Cw.size(); ++i) {
                    V3 d = v_sub(Cw[i], C);
                    double x = v_dot(d, u), y = v_dot(d, v);
                    Node node;
                    node.q = Cw[i];
                    node.ang = std::atan2(y, x);
                    nodes.push_back(node);
                }

                // Simple sort function
                std::sort(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
                    return a.ang < b.ang;
                    });

                std::vector<V3> poly;
                poly.reserve(nodes.size());
                for (size_t i = 0; i < nodes.size(); ++i) {
                    poly.push_back(nodes[i].q);
                }
                *outCutPolygon = poly;
            }
            else {
                outCutPolygon->clear();
            }
        }
        return Q;
    }

    // =========================== Simple face key using basic types ===========================

    struct FaceKey {
        std::vector<int> vertices;

        bool operator<(const FaceKey& other) const {
            if (vertices.size() != other.vertices.size()) {
                return vertices.size() < other.vertices.size();
            }
            for (size_t i = 0; i < vertices.size(); ++i) {
                if (vertices[i] != other.vertices[i]) {
                    return vertices[i] < other.vertices[i];
                }
            }
            return false;
        }
    };

    static int add_global_point(const V3& p, std::vector<std::array<double, 3>>& P,
        std::map<VKey, int>& map, double eps)
    {
        VKey key = make_vertex_key(p, eps);
        std::map<VKey, int>::iterator it = map.find(key);
        if (it != map.end()) return it->second;
        int id = static_cast<int>(P.size());
        std::array<double, 3> arr = { { p.x, p.y, p.z } };
        P.push_back(arr);
        map[key] = id;
        return id;
    }

    // =========================== RVD core ===========================

    static inline Plane bisector_plane(V3 s0, V3 si) {
        V3 n = v_sub(si, s0);
        double c = 0.5 * (v_dot(si, si) - v_dot(s0, s0));
        Plane plane;
        plane.n = n;
        plane.c = c;
        return plane;
    }

    static V3 tet_centroid(const double* pts, const uint32_t* T) {
        V3 a = { pts[3 * T[0] + 0], pts[3 * T[0] + 1], pts[3 * T[0] + 2] };
        V3 b = { pts[3 * T[1] + 0], pts[3 * T[1] + 1], pts[3 * T[1] + 2] };
        V3 c = { pts[3 * T[2] + 0], pts[3 * T[2] + 1], pts[3 * T[2] + 2] };
        V3 d = { pts[3 * T[3] + 0], pts[3 * T[3] + 1], pts[3 * T[3] + 2] };
        V3 s = v_add(v_add(a, b), v_add(c, d));
        return v_scale(s, 0.25);
    }

    static V3 seed_at(const double* S, std::size_t i) {
        return { S[3 * i + 0], S[3 * i + 1], S[3 * i + 2] };
    }

    bool convertTetToRVDPolyMesh(
        const TetMeshView& tet,
        const SeedCloud& seedsIn,
        OpenFOAMPolyMesh& out,
        bool verbose)
    {
        if (!tet.points || tet.numPoints == 0 || !tet.tets || tet.numTets == 0) return false;

        // Calculate epsilon
        double xmin = std::numeric_limits<double>::infinity();
        double ymin = xmin, zmin = xmin, xmax = -xmin, ymax = -xmin, zmax = -xmin;
        for (std::size_t i = 0; i < tet.numPoints; i++) {
            double x = tet.points[3 * i + 0], y = tet.points[3 * i + 1], z = tet.points[3 * i + 2];
            xmin = std::min(xmin, x); xmax = std::max(xmax, x);
            ymin = std::min(ymin, y); ymax = std::max(ymax, y);
            zmin = std::min(zmin, z); zmax = std::max(zmax, z);
        }
        double dx = xmax - xmin, dy = ymax - ymin, dz = zmax - zmin;
        double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double EPS = (diag > 0 ? 1e-9 * diag : 1e-12);

        // Seeds setup
        std::vector<double> seeds;
        const double* S = nullptr;
        std::size_t NS = 0;
        if (seedsIn.xyz && seedsIn.count) {
            S = seedsIn.xyz;
            NS = seedsIn.count;
        }
        else {
            seeds.resize(tet.numTets * 3);
            for (std::size_t e = 0; e < tet.numTets; ++e) {
                const uint32_t* T = tet.tets + 4 * e;
                V3 c = tet_centroid(tet.points, T);
                seeds[3 * e + 0] = c.x;
                seeds[3 * e + 1] = c.y;
                seeds[3 * e + 2] = c.z;
            }
            S = seeds.data();
            NS = tet.numTets;
        }

        out.points.clear();
        out.faces.clear();
        out.owner.clear();
        out.neighbour.clear();

        std::map<VKey, int> pointMap;
        std::map<FaceKey, int> faceMap;
        std::vector<std::vector<int>> cellVerts(NS);

        std::size_t totalCuts = 0;
        for (std::size_t e = 0; e < tet.numTets; ++e) {
            const uint32_t* T = tet.tets + 4 * e;
            V3 p0 = { tet.points[3 * T[0] + 0], tet.points[3 * T[0] + 1], tet.points[3 * T[0] + 2] };
            V3 p1 = { tet.points[3 * T[1] + 0], tet.points[3 * T[1] + 1], tet.points[3 * T[1] + 2] };
            V3 p2 = { tet.points[3 * T[2] + 0], tet.points[3 * T[2] + 1], tet.points[3 * T[2] + 2] };
            V3 p3 = { tet.points[3 * T[3] + 0], tet.points[3 * T[3] + 1], tet.points[3 * T[3] + 2] };

            std::size_t ownerSeed = (seedsIn.xyz && seedsIn.count) ? (e % NS) : e;
            V3 s0 = seed_at(S, ownerSeed);

            Poly3 cell = make_tet_poly(p0, p1, p2, p3);

            for (std::size_t j = 0; j < NS; ++j) {
                if (j == ownerSeed) continue;
                V3 si = seed_at(S, j);
                Plane H = bisector_plane(s0, si);

                std::vector<V3> cutPoly;
                Poly3 clipped = clip_polyhedron(cell, H, EPS, &cutPoly);

                if (clipped.V.empty() || clipped.F.empty()) {
                    cell.V.clear();
                    cell.F.clear();
                    break;
                }

                if (cutPoly.size() >= 3) {
                    std::vector<int> vids;
                    vids.reserve(cutPoly.size());
                    for (size_t k = 0; k < cutPoly.size(); ++k) {
                        int gi = add_global_point(cutPoly[k], out.points, pointMap, EPS);
                        vids.push_back(gi);
                        cellVerts[ownerSeed].push_back(gi);
                    }

                    // Register face
                    FaceKey key;
                    key.vertices = vids;
                    std::sort(key.vertices.begin(), key.vertices.end());

                    std::map<FaceKey, int>::iterator it = faceMap.find(key);
                    if (it == faceMap.end()) {
                        int idx = static_cast<int>(out.faces.size());
                        faceMap[key] = idx;
                        out.faces.push_back(vids);
                        out.owner.push_back(static_cast<int>(ownerSeed));
                        out.neighbour.push_back(static_cast<int>(j));
                    }

                    totalCuts++;
                }

                cell = clipped;
            }

            // Domain boundary faces
            for (size_t fi = 0; fi < cell.F.size(); ++fi) {
                const std::vector<int>& f = cell.F[fi];
                std::vector<int> vids;
                vids.reserve(f.size());
                for (size_t i = 0; i < f.size(); ++i) {
                    const V3& q = cell.V[f[i]];
                    int gi = add_global_point(q, out.points, pointMap, EPS);
                    vids.push_back(gi);
                    cellVerts[ownerSeed].push_back(gi);
                }

                if (vids.size() < 3) continue;

                FaceKey key;
                key.vertices = vids;
                std::sort(key.vertices.begin(), key.vertices.end());

                if (faceMap.find(key) == faceMap.end()) {
                    int idx = static_cast<int>(out.faces.size());
                    faceMap[key] = idx;
                    out.faces.push_back(vids);
                    out.owner.push_back(static_cast<int>(ownerSeed));
                }
            }
        }

        if (verbose) {
            std::cout << "[RVD] seeds=" << NS << "  tets=" << tet.numTets
                << "  internalCuts=" << totalCuts
                << "  faces(total)=" << out.faces.size()
                << "  points=" << out.points.size() << "\n";
        }

        // Orientation fix
        std::vector<V3> ownerC(NS);
        for (std::size_t cid = 0; cid < NS; ++cid) {
            ownerC[cid].x = ownerC[cid].y = ownerC[cid].z = 0.0;

            std::vector<int>& vv = cellVerts[cid];
            std::sort(vv.begin(), vv.end());
            vv.erase(std::unique(vv.begin(), vv.end()), vv.end());

            V3 acc = { 0,0,0 };
            for (size_t i = 0; i < vv.size(); ++i) {
                int gi = vv[i];
                if (gi >= 0 && gi < static_cast<int>(out.points.size())) {
                    acc.x += out.points[gi][0];
                    acc.y += out.points[gi][1];
                    acc.z += out.points[gi][2];
                }
            }
            if (!vv.empty()) {
                double inv = 1.0 / double(vv.size());
                ownerC[cid].x = acc.x * inv;
                ownerC[cid].y = acc.y * inv;
                ownerC[cid].z = acc.z * inv;
            }
        }

        const int n_internal = static_cast<int>(out.neighbour.size());

        std::vector<std::vector<int>> faces_internal;
        std::vector<int> owners_internal, neigh_internal;

        std::vector<std::vector<int>> faces_boundary;
        std::vector<int> owners_boundary;

        // Process internal faces
        for (int i = 0; i < n_internal; i++) {
            std::vector<int> vids = out.faces[i];
            int own = out.owner[i];
            int nei = out.neighbour[i];

            if (own < 0 || own >= static_cast<int>(NS) || nei < 0 || nei >= static_cast<int>(NS)) continue;

            std::vector<V3> P;
            P.reserve(vids.size());
            for (size_t j = 0; j < vids.size(); ++j) {
                int gi = vids[j];
                if (gi >= 0 && gi < static_cast<int>(out.points.size())) {
                    V3 point;
                    point.x = out.points[gi][0];
                    point.y = out.points[gi][1];
                    point.z = out.points[gi][2];
                    P.push_back(point);
                }
            }

            if (P.size() < 3) continue;

            V3 Cf = poly_centroid(P);
            V3 n = newell_normal(P);
            V3 Co = ownerC[own];
            if (v_dot(n, v_sub(Cf, Co)) < 0.0) {
                std::reverse(vids.begin(), vids.end());
            }
            faces_internal.push_back(vids);
            owners_internal.push_back(own);
            neigh_internal.push_back(nei);
        }

        // Process boundary faces
        for (size_t i = n_internal; i < out.faces.size(); ++i) {
            std::vector<int> vids = out.faces[i];
            int own = out.owner[i];

            if (own < 0 || own >= static_cast<int>(NS)) continue;

            std::vector<V3> P;
            P.reserve(vids.size());
            for (size_t j = 0; j < vids.size(); ++j) {
                int gi = vids[j];
                if (gi >= 0 && gi < static_cast<int>(out.points.size())) {
                    V3 point;
                    point.x = out.points[gi][0];
                    point.y = out.points[gi][1];
                    point.z = out.points[gi][2];
                    P.push_back(point);
                }
            }

            if (P.size() < 3) continue;

            V3 Cf = poly_centroid(P);
            V3 n = newell_normal(P);
            V3 Co = ownerC[own];
            if (v_dot(n, v_sub(Cf, Co)) < 0.0) {
                std::reverse(vids.begin(), vids.end());
            }
            faces_boundary.push_back(vids);
            owners_boundary.push_back(own);
        }

        // Rebuild output
        out.faces.clear();
        out.owner.clear();

        for (size_t i = 0; i < faces_internal.size(); ++i) {
            out.faces.push_back(faces_internal[i]);
            out.owner.push_back(owners_internal[i]);
        }
        out.neighbour = neigh_internal;

        for (size_t i = 0; i < faces_boundary.size(); ++i) {
            out.faces.push_back(faces_boundary[i]);
            out.owner.push_back(owners_boundary[i]);
        }

        if (verbose) {
            std::cout << "[RVD] faces: internal=" << out.neighbour.size()
                << " boundary=" << (out.faces.size() - out.neighbour.size())
                << " total=" << out.faces.size() << "\n";
        }
        return true;
    }

    // =========================== OpenFOAM writer ===========================

    static void writeHeader(std::ostream& os, const std::string& cls, const std::string& obj) {
        os <<
            "FoamFile\n"
            "{\n"
            "    version     2.0;\n"
            "    format      ascii;\n"
            "    class       " << cls << ";\n"
            "    location    \"constant/polyMesh\";\n"
            "    object      " << obj << ";\n"
            "}\n\n";
    }

    static bool ensureDir(const std::string& path) {
#if defined(_WIN32)
        std::string cmd = "mkdir \"" + path + "\" >nul 2>nul";
#else
        std::string cmd = "mkdir -p \"" + path + "\"";
#endif
        return std::system(cmd.c_str()) >= 0;
    }

    bool writeOpenFOAMPolyMesh(
        const OpenFOAMPolyMesh& m,
        const std::string& polyMeshDir,
        bool verbose)
    {
        if (!ensureDir(polyMeshDir)) return false;

        if (m.points.empty() || m.faces.empty() || m.owner.empty()) {
            std::cerr << "ERROR: Invalid mesh data\n";
            return false;
        }

        // points
        {
            std::ofstream os(polyMeshDir + "/points");
            if (!os) return false;
            writeHeader(os, "vectorField", "points");
            os << m.points.size() << "\n(\n";
            for (size_t i = 0; i < m.points.size(); ++i) {
                const std::array<double, 3>& p = m.points[i];
                os << "(" << std::scientific << std::setprecision(15)
                    << p[0] << " " << p[1] << " " << p[2] << ")\n";
            }
            os << ")\n";
        }

        // faces
        {
            std::ofstream os(polyMeshDir + "/faces");
            if (!os) return false;
            writeHeader(os, "faceList", "faces");
            os << m.faces.size() << "\n(\n";
            for (size_t i = 0; i < m.faces.size(); ++i) {
                const std::vector<int>& f = m.faces[i];
                if (f.size() >= 3) {
                    os << f.size() << "(";
                    for (std::size_t j = 0; j < f.size(); ++j) {
                        os << f[j] << (j + 1 < f.size() ? " " : "");
                    }
                    os << ")\n";
                }
            }
            os << ")\n";
        }

        // owner
        {
            std::ofstream os(polyMeshDir + "/owner");
            if (!os) return false;
            writeHeader(os, "labelList", "owner");
            os << m.owner.size() << "\n(\n";
            for (size_t i = 0; i < m.owner.size(); ++i) {
                os << m.owner[i] << "\n";
            }
            os << ")\n";
        }

        // neighbour
        {
            std::ofstream os(polyMeshDir + "/neighbour");
            if (!os) return false;
            writeHeader(os, "labelList", "neighbour");
            os << m.neighbour.size() << "\n(\n";
            for (size_t i = 0; i < m.neighbour.size(); ++i) {
                os << m.neighbour[i] << "\n";
            }
            os << ")\n";
        }

        // boundary
        {
            const int startFace = static_cast<int>(m.neighbour.size());
            const int nBoundary = static_cast<int>(m.faces.size()) - startFace;
            std::ofstream os(polyMeshDir + "/boundary");
            if (!os) return false;
            writeHeader(os, "polyBoundaryMesh", "boundary");
            os << "1\n(\n";
            os << "    " << m.boundaryName << "\n";
            os << "    {\n";
            os << "        type            " << m.boundaryType << ";\n";
            os << "        nFaces          " << nBoundary << ";\n";
            os << "        startFace       " << startFace << ";\n";
            os << "    }\n";
            os << ")\n";
        }

        if (verbose) {
            std::cout << "[RVD->OpenFOAM] Wrote polyMesh to: " << polyMeshDir << "\n";
        }
        return true;
    }

    bool convertTetToOpenFOAM(
        const TetMeshView& tet,
        const SeedCloud& seeds,
        const std::string& casePolyMeshDir,
        bool verbose)
    {
        OpenFOAMPolyMesh M;
        if (!convertTetToRVDPolyMesh(tet, seeds, M, verbose)) return false;
        return writeOpenFOAMPolyMesh(M, casePolyMeshDir, verbose);
    }

} // namespace rvd2foam