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
#include <unordered_map>
#include <utility>
#include <iostream>
#include <vector>

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

    struct Plane { V3 n; double c; }; // half-space: n·x <= c

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
        V3 c{ 0,0,0 }; if (poly.empty()) return c;
        for (auto& p : poly) { c = v_add(c, p); }
        return v_scale(c, 1.0 / double(poly.size()));
    }

    static inline V3 newell_normal(const std::vector<V3>& P) {
        double nx = 0, ny = 0, nz = 0; size_t n = P.size(); if (n < 3) return { 0,0,0 };
        for (size_t i = 0;i < n;i++) {
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
        std::vector<V3> V;                // vertices
        std::vector<std::vector<int>> F;  // faces as CCW vertex indices (arbitrary start)
        // (We will fix orientation later against owner centroids)
    };

    static Poly3 make_tet_poly(const V3& p0, const V3& p1, const V3& p2, const V3& p3) {
        Poly3 P;
        P.V = { p0,p1,p2,p3 };
        // 4 faces: opposite each vertex
        // We'll use consistent vertex cycles; orientation will be corrected later.
        P.F.push_back({ 0,2,1 }); // opp v3
        P.F.push_back({ 0,1,3 }); // opp v2
        P.F.push_back({ 1,2,3 }); // opp v0
        P.F.push_back({ 2,0,3 }); // opp v1
        return P;
    }

    // Clip a polygon by plane n·x <= c (3D; polygon lies in some plane).
    // Sutherland–Hodgman generalization in 3D (in/out test by signed distance).
    static std::vector<V3> clip_polygon(const std::vector<V3>& poly, const Plane& P, double eps) {
        std::vector<V3> out; if (poly.empty()) return out;
        auto inside = [&](const V3& q) { return signed_dist(P, q) <= eps; };
        V3 S = poly.back();
        bool Sin = inside(S);
        for (const V3& E : poly) {
            bool Ein = inside(E);
            if (Ein) {
                if (!Sin) {
                    // S->E crosses plane: add intersection
                    V3 d = v_sub(E, S);
                    double denom = v_dot(P.n, d);
                    double t = (P.c - v_dot(P.n, S)) / (denom == 0 ? 1e-30 : denom);
                    t = std::max(0.0, std::min(1.0, t));
                    out.push_back(lerp(S, E, t));
                }
                out.push_back(E);
            }
            else if (Sin) {
                // leaving: add intersection
                V3 d = v_sub(E, S);
                double denom = v_dot(P.n, d);
                double t = (P.c - v_dot(P.n, S)) / (denom == 0 ? 1e-30 : denom);
                t = std::max(0.0, std::min(1.0, t));
                out.push_back(lerp(S, E, t));
            }
            S = E; Sin = Ein;
        }
        // collapse nearly-duplicate consecutive vertices
        std::vector<V3> comp;
        for (const auto& q : out) {
            if (comp.empty() || v_norm(v_sub(q, comp.back())) > eps) comp.push_back(q);
        }
        if (comp.size() >= 3 && v_norm(v_sub(comp.front(), comp.back())) <= eps) {
            comp.back() = comp.front();
            comp.pop_back();
        }
        return comp;
    }

    // Clip a convex polyhedron P by plane n·x <= c.
    // Returns new convex poly + the new "cut face" polygon (if clipping occurred).
    static Poly3 clip_polyhedron(const Poly3& Pin, const Plane& P, double eps, std::vector<V3>* outCutPolygon = nullptr) {
        Poly3 Q;
        std::vector<V3> cutLoopAccum;

        // Clip each face polygon
        for (const auto& face : Pin.F) {
            std::vector<V3> poly; poly.reserve(face.size());
            for (int vi : face) poly.push_back(Pin.V[std::size_t(vi)]);
            auto clipped = clip_polygon(poly, P, eps);
            if (clipped.size() >= 3) {
                // store new face; indices will be remapped after vertex weld
                int base = int(Q.V.size());
                for (const auto& v : clipped) Q.V.push_back(v);
                std::vector<int> idx(clipped.size());
                for (size_t k = 0;k < clipped.size();++k) idx[k] = base + int(k);
                Q.F.push_back(std::move(idx));
            }

            // Collect intersection segments endpoints to reconstruct the new cut face.
            // We approximate by collecting all points that lie (within eps) on the plane.
            for (const auto& v : clipped) {
                if (std::fabs(signed_dist(P, v)) <= 2 * eps) cutLoopAccum.push_back(v);
            }
        }

        // Weld vertices: build mapping old->new unique by rounding
        auto keyOf = [&](const V3& a) {
            const double s = 1.0 / eps;
            long long ix = (long long)std::llround(a.x * s);
            long long iy = (long long)std::llround(a.y * s);
            long long iz = (long long)std::llround(a.z * s);
            return std::tuple<long long, long long, long long>(ix, iy, iz);
            };
        std::unordered_map<std::tuple<long long, long long, long long>, int> vmap;
        std::vector<V3> WV; WV.reserve(Q.V.size());
        std::vector<int> remap(Q.V.size(), -1);
        for (size_t i = 0;i < Q.V.size();++i) {
            auto K = keyOf(Q.V[i]);
            auto it = vmap.find(K);
            if (it == vmap.end()) {
                int nid = int(WV.size());
                vmap.emplace(K, nid);
                WV.push_back(Q.V[i]);
                remap[i] = nid;
            }
            else {
                remap[i] = it->second;
            }
        }
        for (auto& f : Q.F) {
            for (int& vi : f) vi = remap[std::size_t(vi)];
            // remove duplicate consecutive indices
            std::vector<int> g; g.reserve(f.size());
            for (int k = 0;k < (int)f.size();++k) {
                if (g.empty() || WV[std::size_t(g.back())].x != WV[std::size_t(f[k])].x
                    || WV[std::size_t(g.back())].y != WV[std::size_t(f[k])].y
                    || WV[std::size_t(g.back())].z != WV[std::size_t(f[k])].z) {
                    g.push_back(f[k]);
                }
            }
            f.swap(g);
        }
        Q.V.swap(WV);

        // Remove degenerate faces (<3 verts or collinear)
        auto is_degenerate = [&](const std::vector<int>& f)->bool {
            if (f.size() < 3) return true;
            V3 a = Q.V[std::size_t(f[0])];
            V3 b = Q.V[std::size_t(f[1])];
            V3 c = Q.V[std::size_t(f[2])];
            return v_norm(tri_normal(a, b, c)) < 1e-18;
            };
        std::vector<std::vector<int>> Fclean; Fclean.reserve(Q.F.size());
        for (const auto& f : Q.F) if (!is_degenerate(f)) Fclean.push_back(f);
        Q.F.swap(Fclean);

        // Build the new cut face (if we truly cut off something).
        if (outCutPolygon) {
            // Weld and sort the on-plane points to a cycle
            // 1) Weld
            std::vector<V3> Cw;
            std::set<std::tuple<long long, long long, long long>> seen;
            for (const auto& q : cutLoopAccum) {
                auto K = keyOf(q);
                if (seen.insert(K).second) Cw.push_back(q);
            }
            // 2) If enough points, sort them around plane normal to form a polygon
            if (Cw.size() >= 3) {
                // Compute local 2D basis on plane
                V3 n = P.n;
                double nn = v_norm(n);
                if (nn < 1e-30) n = { 0,0,1 }; else n = v_scale(n, 1.0 / nn);
                V3 tmp = (std::fabs(n.x) < 0.9 ? V3{ 1,0,0 } : V3{ 0,1,0 });
                V3 u = v_scale(v_sub(tmp, v_scale(n, v_dot(n, tmp))), 1.0 / std::max(1e-30, v_norm(v_sub(tmp, v_scale(n, v_dot(n, tmp))))));
                V3 v = v_cross(n, u);
                // centroid
                V3 C{ 0,0,0 }; for (auto& q : Cw) C = v_add(C, q);
                C = v_scale(C, 1.0 / double(Cw.size()));
                // angle sort
                struct Node { V3 q; double ang; };
                std::vector<Node> nodes; nodes.reserve(Cw.size());
                for (auto& q : Cw) {
                    V3 d = v_sub(q, C);
                    double x = v_dot(d, u), y = v_dot(d, v);
                    nodes.push_back({ q, std::atan2(y,x) });
                }
                std::sort(nodes.begin(), nodes.end(), [](const Node& A, const Node& B) { return A.ang < B.ang; });
                std::vector<V3> poly; poly.reserve(nodes.size());
                for (auto& nd : nodes) poly.push_back(nd.q);
                *outCutPolygon = poly;
            }
            else {
                outCutPolygon->clear();
            }
        }
        return Q;
    }

    // ====================== face & vertex global bookkeeping ======================

    struct FaceRec { std::vector<int> vid; int owner = -1; int neighbour = -1; };

    struct VHash {
        std::size_t operator()(const std::tuple<long long, long long, long long>& k) const noexcept {
            auto [a, b, c] = k;
            std::size_t h = 1469598103934665603ull;
            h ^= std::hash<long long>{}(a); h *= 1099511628211ull;
            h ^= std::hash<long long>{}(b); h *= 1099511628211ull;
            h ^= std::hash<long long>{}(c); h *= 1099511628211ull;
            return h;
        }
    };

    static int add_global_point(const V3& p, std::vector<std::array<double, 3>>& P,
        std::unordered_map<std::tuple<long long, long long, long long>, int, VHash>& map,
        double eps)
    {
        const double s = 1.0 / eps;
        auto key = std::make_tuple(std::llround(p.x * s), std::llround(p.y * s), std::llround(p.z * s));
        auto it = map.find(key);
        if (it != map.end()) return it->second;
        int id = (int)P.size();
        P.push_back({ p.x,p.y,p.z });
        map.emplace(key, id);
        return id;
    }

    static std::size_t face_key_hash(const std::vector<int>& v) {
        std::vector<int> k = v; std::sort(k.begin(), k.end());
        std::size_t h = 1469598103934665603ull;
        for (int x : k) { h ^= std::size_t(x + 0x9e3779b9); h *= 1099511628211ull; }
        return h;
    }
    struct FKey { std::vector<int> k; };
    struct FKeyHasher {
        std::size_t operator()(const FKey& f) const noexcept {
            return face_key_hash(f.k);
        }
    };
    struct FKeyEq {
        bool operator()(const FKey& a, const FKey& b) const noexcept {
            if (a.k.size() != b.k.size()) return false;
            for (size_t i = 0;i < a.k.size();++i) if (a.k[i] != b.k[i]) return false;
            return true;
        }
    };

    // =========================== RVD core ===========================

    // Build plane for bisector between s0 and si: { x | ||x-s0|| <= ||x-si|| }.
    // Inequality: 2(si-s0)·x <= ||si||^2 - ||s0||^2  -> n·x <= c
    static inline Plane bisector_plane(V3 s0, V3 si) {
        V3 n = v_sub(si, s0);
        double c = 0.5 * (v_dot(si, si) - v_dot(s0, s0));
        return { n,c };
    }

    static V3 tet_centroid(const double* pts, const uint32_t* T) {
        V3 a{ pts[3 * T[0] + 0], pts[3 * T[0] + 1], pts[3 * T[0] + 2] };
        V3 b{ pts[3 * T[1] + 0], pts[3 * T[1] + 1], pts[3 * T[1] + 2] };
        V3 c{ pts[3 * T[2] + 0], pts[3 * T[2] + 1], pts[3 * T[2] + 2] };
        V3 d{ pts[3 * T[3] + 0], pts[3 * T[3] + 1], pts[3 * T[3] + 2] };
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

        // Epsilon (merge/clip tolerance) ~ 1e-9 of bbox diagonal
        double xmin = +std::numeric_limits<double>::infinity();
        double ymin = xmin, zmin = xmin, xmax = -xmin, ymax = -xmin, zmax = -xmin;
        for (std::size_t i = 0;i < tet.numPoints;i++) {
            double x = tet.points[3 * i + 0], y = tet.points[3 * i + 1], z = tet.points[3 * i + 2];
            xmin = std::min(xmin, x); xmax = std::max(xmax, x);
            ymin = std::min(ymin, y); ymax = std::max(ymax, y);
            zmin = std::min(zmin, z); zmax = std::max(zmax, z);
        }
        double diag = std::sqrt((xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin) + (zmax - zmin) * (zmax - zmin));
        const double EPS = std::max(1e-12, 1e-9 * diag);

        // Seeds: use provided or per-tet centroids
        std::vector<double> seeds;
        const double* S = nullptr;
        std::size_t   NS = 0;
        if (seedsIn.xyz && seedsIn.count) {
            S = seedsIn.xyz; NS = seedsIn.count;
        }
        else {
            seeds.resize(tet.numTets * 3);
            for (std::size_t e = 0;e < tet.numTets;++e) {
                const uint32_t* T = tet.tets + 4 * e;
                V3 c = tet_centroid(tet.points, T);
                seeds[3 * e + 0] = c.x; seeds[3 * e + 1] = c.y; seeds[3 * e + 2] = c.z;
            }
            S = seeds.data(); NS = tet.numTets;
        }

        // Global vertex map
        out.points.clear(); out.faces.clear(); out.owner.clear(); out.neighbour.clear();
        std::unordered_map<std::tuple<long long, long long, long long>, int, VHash> pointMap;
        pointMap.reserve(tet.numTets * 8);

        // Global face map: sorted vertex indices -> face index in out.faces
        std::unordered_map<FKey, int, FKeyHasher, FKeyEq> faceMap;
        faceMap.reserve(tet.numTets * 16);

        // Accumulate per-cell vertices to estimate centroids for orientation fix
        std::vector<std::vector<int>> cellVerts(NS);

        auto register_face = [&](const std::vector<int>& vtx, int owner, int neighbour) {
            // canonical key
            FKey key; key.k = vtx; std::sort(key.k.begin(), key.k.end());
            auto it = faceMap.find(key);
            if (it == faceMap.end()) {
                int idx = (int)out.faces.size();
                faceMap.emplace(key, idx);
                out.faces.push_back(vtx);
                out.owner.push_back(owner);
                if (neighbour >= 0) out.neighbour.push_back(neighbour);
            }
            else {
                int idx = it->second;
                // set neighbour only if not yet set (neighbour is stored in the internal prefix)
                // If this face was previously added as boundary, we cannot rewrite neighbour inline;
                // but by construction we only register internal "cut" faces here before boundary.
                int internalCount = (int)out.neighbour.size();
                if (idx < internalCount) {
                    // already internal; nothing to do
                }
                else {
                    // This is unexpected in this simple pass; ignore.
                }
            }
            };

        // RVD by per-tet clipping
        // NOTE: For simplicity, we consider all other seeds when clipping.
        // We’ll optimize later by culling distant seeds.
        std::size_t totalCuts = 0;
        for (std::size_t e = 0;e < tet.numTets;++e) {
            const uint32_t* T = tet.tets + 4 * e;
            V3 p0{ tet.points[3 * T[0] + 0], tet.points[3 * T[0] + 1], tet.points[3 * T[0] + 2] };
            V3 p1{ tet.points[3 * T[1] + 0], tet.points[3 * T[1] + 1], tet.points[3 * T[1] + 2] };
            V3 p2{ tet.points[3 * T[2] + 0], tet.points[3 * T[2] + 1], tet.points[3 * T[2] + 2] };
            V3 p3{ tet.points[3 * T[3] + 0], tet.points[3 * T[3] + 1], tet.points[3 * T[3] + 2] };

            // If no explicit seeds: seed 0 for this tet is its centroid (index = e)
            // Otherwise, you can choose a policy; here we still use the centroid index as the "owner id".
            // You can map tet->seed id differently later.
            std::size_t ownerSeed = (seedsIn.xyz && seedsIn.count) ? e % NS : e;

            V3 s0 = seed_at(S, ownerSeed);

            // Start with the tet as Ωᵉ
            Poly3 cell = make_tet_poly(p0, p1, p2, p3);

            // For all other seeds si, clip by the perpendicular bisector (keep closer to s0)
            for (std::size_t j = 0;j < NS;++j) {
                if (j == ownerSeed) continue;
                V3 si = seed_at(S, j);
                Plane H = bisector_plane(s0, si);

                std::vector<V3> cutPoly;
                Poly3 clipped = clip_polyhedron(cell, H, EPS, &cutPoly);

                if (clipped.V.empty() || clipped.F.empty()) {
                    cell.V.clear(); cell.F.clear(); break; // nothing remains for s0 in this tet
                }

                // If this plane actually cut the cell, we got a cut face polygon
                if (cutPoly.size() >= 3) {
                    // Register cut face as an internal face between (ownerSeed, neighbour=j)
                    std::vector<int> vids; vids.reserve(cutPoly.size());
                    for (const auto& q : cutPoly) {
                        int gi = add_global_point(q, out.points, pointMap, EPS);
                        vids.push_back(gi);
                        cellVerts[ownerSeed].push_back(gi);
                    }
                    // Face orientation will be fixed globally later.
                    register_face(vids, (int)ownerSeed, (int)j);
                    totalCuts++;
                }
                std::swap(cell, clipped);
            }

            // Any remaining faces lying on the original tet boundary are domain boundary faces.
            // Add them now as boundary faces (neighbour = -1).
            for (const auto& f : cell.F) {
                // Build polygon geometry
                std::vector<V3> poly; poly.reserve(f.size());
                for (int vi : f) poly.push_back(cell.V[std::size_t(vi)]);
                // Convert to global indices
                std::vector<int> vids; vids.reserve(f.size());
                for (const auto& q : poly) {
                    int gi = add_global_point(q, out.points, pointMap, EPS);
                    vids.push_back(gi);
                    cellVerts[ownerSeed].push_back(gi);
                }
                // Register as boundary face
                // (May be duplicated across adjacent tets/owners; dedup map will coalesce.)
                FKey key; key.k = vids; std::sort(key.k.begin(), key.k.end());
                if (faceMap.find(key) == faceMap.end()) {
                    int idx = (int)out.faces.size();
                    faceMap.emplace(key, idx);
                    out.faces.push_back(vids);
                    out.owner.push_back((int)ownerSeed);
                    // boundary faces are appended *after* internal; neighbour omitted
                }
            }
        }

        if (verbose) {
            std::cout << "[RVD] seeds=" << NS << "  tets=" << tet.numTets
                << "  internalCuts=" << totalCuts
                << "  faces(total)=" << out.faces.size()
                << "  points=" << out.points.size() << "\n";
        }

        // ---------- Orientation fix (OpenFOAM expects owner-facing outward) ----------
        // Build crude owner centroids from unique vertex averages.
        std::vector<V3> ownerC(NS, V3{ 0,0,0 });
        std::vector<int> ownerN(NS, 0);
        for (std::size_t cid = 0; cid < NS; ++cid) {
            auto& vv = cellVerts[cid];
            std::sort(vv.begin(), vv.end());
            vv.erase(std::unique(vv.begin(), vv.end()), vv.end());
            V3 acc{ 0,0,0 };
            for (int gi : vv) {
                acc.x += out.points[std::size_t(gi)][0];
                acc.y += out.points[std::size_t(gi)][1];
                acc.z += out.points[std::size_t(gi)][2];
            }
            if (!vv.empty()) {
                double inv = 1.0 / double(vv.size());
                ownerC[cid] = { acc.x * inv, acc.y * inv, acc.z * inv };
                ownerN[cid] = (int)vv.size();
            }
        }

        // Split into internal and boundary lists as OpenFOAM expects:
        // neighbour applies only to internal faces and must be a contiguous prefix.
        std::vector<std::vector<int>> faces_internal;
        std::vector<int> owners_internal, neigh_internal;
        faces_internal.reserve(out.faces.size());
        owners_internal.reserve(out.faces.size());
        neigh_internal.reserve(out.faces.size());

        std::vector<std::vector<int>> faces_boundary;
        std::vector<int> owners_boundary;
        faces_boundary.reserve(out.faces.size());
        owners_boundary.reserve(out.faces.size());

        // First pass: detect which faces are internal (have neighbour registered earlier).
        // In our simple assembler, internal faces got registered first via 'register_face'
        // and boundary faces were added afterwards without neighbour entries.
        // So we use a simple rule: a face is internal if it appears (owner, neighbour) pair existed
        // in the first pass. Here we recompute: if both owner & neighbour cells saw this face,
        // it will appear twice; dedup earlier kept the first as owner and we pushed 'neighbour'.
        // We can detect internal by finding any other face with the exact same sorted-vertex key
        // in the prefix where neighbour exists; but we already constructed the neighbour list
        // as we went. So we take 'neighbour.size()' as internal count.
        const int n_internal = (int)out.neighbour.size();

        for (int i = 0;i < n_internal;i++) {
            auto vids = out.faces[std::size_t(i)];
            int own = out.owner[std::size_t(i)];
            int nei = out.neighbour[std::size_t(i)];
            // Orientation: outward from owner
            std::vector<V3> P; P.reserve(vids.size());
            for (int gi : vids) P.push_back(V3{ out.points[std::size_t(gi)][0],
                                                out.points[std::size_t(gi)][1],
                                                out.points[std::size_t(gi)][2] });
            V3 Cf = poly_centroid(P);
            V3 n = newell_normal(P);
            V3 Co = ownerC[std::size_t(own)];
            if (v_dot(n, v_sub(Cf, Co)) < 0.0) {
                std::reverse(vids.begin(), vids.end());
            }
            faces_internal.push_back(std::move(vids));
            owners_internal.push_back(own);
            neigh_internal.push_back(nei);
        }

        for (std::size_t i = n_internal; i < out.faces.size(); ++i) {
            auto vids = out.faces[i];
            int own = out.owner[i];
            std::vector<V3> P; P.reserve(vids.size());
            for (int gi : vids) P.push_back(V3{ out.points[std::size_t(gi)][0],
                                                out.points[std::size_t(gi)][1],
                                                out.points[std::size_t(gi)][2] });
            V3 Cf = poly_centroid(P);
            V3 n = newell_normal(P);
            V3 Co = ownerC[std::size_t(own)];
            if (v_dot(n, v_sub(Cf, Co)) < 0.0) {
                std::reverse(vids.begin(), vids.end());
            }
            faces_boundary.push_back(std::move(vids));
            owners_boundary.push_back(own);
        }

        // Repack in OpenFOAM order: internal faces first (have neighbour),
        // boundary faces after (no neighbour entries).
        out.faces.clear();
        out.owner.clear();
        // keep neighbour content
        std::vector<int> N = std::move(out.neighbour);

        for (size_t i = 0;i < faces_internal.size();++i) {
            out.faces.push_back(std::move(faces_internal[i]));
            out.owner.push_back(owners_internal[i]);
        }
        // neighbour stays intact (same count)
        out.neighbour = std::move(neigh_internal);

        for (size_t i = 0;i < faces_boundary.size();++i) {
            out.faces.push_back(std::move(faces_boundary[i]));
            out.owner.push_back(owners_boundary[i]);
        }

        if (verbose) {
            std::cout << "[RVD] faces: internal=" << out.neighbour.size()
                << " boundary=" << (out.faces.size() - out.neighbour.size())
                << " total=" << out.faces.size() << "\n";
        }
        return true;
    }

    // =========================== OpenFOAM writer (ASCII) ===========================

    static void writeHeader(std::ostream& os, const std::string& cls, const std::string& obj) {
        os <<
            "FoamFile\n"
            "{\n"
            "    version     2.0;\n"
            "    format      ascii;\n"
            "    class       " << cls << ";\n"
            "    location    \"polyMesh\";\n"
            "    object      " << obj << ";\n"
            "}\n\n";
    }

    // tiny mkdir wrapper
    static bool ensureDir(const std::string& path) {
#if defined(_WIN32)
        std::string cmd = "mkdir \"" + path + "\" >nul 2>nul";
#else
        std::string cmd = "mkdir -p \"" + path + "\"";
#endif
        return std::system(cmd.c_str()) == 0;
    }

    bool writeOpenFOAMPolyMesh(
        const OpenFOAMPolyMesh& m,
        const std::string& polyMeshDir,
        bool verbose)
    {
        if (!ensureDir(polyMeshDir)) return false;

        // points
        {
            std::ofstream os(polyMeshDir + "/points");
            writeHeader(os, "vectorField", "points");
            os << m.points.size() << "\n(\n";
            for (const auto& p : m.points) {
                os << "(" << std::setprecision(17) << p[0] << " " << p[1] << " " << p[2] << ")\n";
            }
            os << ")\n";
        }
        // faces
        {
            std::ofstream os(polyMeshDir + "/faces");
            writeHeader(os, "faceList", "faces");
            os << m.faces.size() << "\n(\n";
            for (const auto& f : m.faces) {
                os << f.size() << "(";
                for (std::size_t i = 0;i < f.size();++i) {
                    os << f[i] << (i + 1 < f.size() ? " " : "");
                }
                os << ")\n";
            }
            os << ")\n";
        }
        // owner
        {
            std::ofstream os(polyMeshDir + "/owner");
            writeHeader(os, "labelList", "owner");
            os << m.owner.size() << "\n(\n";
            for (int c : m.owner) os << c << "\n";
            os << ")\n";
        }
        // neighbour (internal faces only)
        {
            std::ofstream os(polyMeshDir + "/neighbour");
            writeHeader(os, "labelList", "neighbour");
            os << m.neighbour.size() << "\n(\n";
            for (int c : m.neighbour) os << c << "\n";
            os << ")\n";
        }
        // boundary (single patch for now)
        {
            const int startFace = (int)m.neighbour.size();
            const int nBoundary = (int)m.faces.size() - startFace;
            std::ofstream os(polyMeshDir + "/boundary");
            writeHeader(os, "polyBoundaryMesh", "boundary");
            os << 1 << "\n(\n";
            os << "    " << m.boundaryName << "\n    {\n";
            os << "        type            " << m.boundaryType << ";\n";
            os << "        nFaces          " << nBoundary << ";\n";
            os << "        startFace       " << startFace << ";\n";
            os << "    }\n";
            os << ")\n";
        }

        if (verbose) {
            std::cout << "[RVD->OpenFOAM] Wrote polyMesh to: " << polyMeshDir << "\n";
            std::cout << "  points=" << m.points.size()
                << " faces=" << m.faces.size()
                << " owner=" << m.owner.size()
                << " neighbour=" << m.neighbour.size()
                << " patches=1\n";
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
