// rvd_geometry.cpp — implementation for math & clipping

#if defined(_WIN32) && !defined(NOMINMAX)
#define NOMINMAX
#endif

#include "rvd_geometry.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <tuple>
#include <unordered_map>

namespace rvd {

    V3 v_add(V3 a, V3 b) { return { a.x + b.x,a.y + b.y,a.z + b.z }; }
    V3 v_sub(V3 a, V3 b) { return { a.x - b.x,a.y - b.y,a.z - b.z }; }
    V3 v_scale(V3 a, double s) { return { a.x * s,a.y * s,a.z * s }; }
    double v_dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
    V3 v_cross(V3 a, V3 b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }
    double v_norm(V3 a) { return std::sqrt(v_dot(a, a)); }

    V3 newell_normal(const std::vector<V3>& P) {
        if (P.size() < 3) return { 0,0,0 };
        double nx = 0, ny = 0, nz = 0;
        for (size_t i = 0;i < P.size();++i) {
            const V3& a = P[i]; const V3& b = P[(i + 1) % P.size()];
            nx += (a.y - b.y) * (a.z + b.z);
            ny += (a.z - b.z) * (a.x + b.x);
            nz += (a.x - b.x) * (a.y + b.y);
        }
        return { nx,ny,nz };
    }
    V3 poly_centroid(const std::vector<V3>& P) {
        V3 c{ 0,0,0 }; if (P.empty()) return c;
        for (const auto& p : P) c = v_add(c, p);
        return v_scale(c, 1.0 / double(P.size()));
    }

    V3 tet_centroid(const double* pts, const uint32_t* T) {
        V3 a{ pts[3 * T[0] + 0],pts[3 * T[0] + 1],pts[3 * T[0] + 2] };
        V3 b{ pts[3 * T[1] + 0],pts[3 * T[1] + 1],pts[3 * T[1] + 2] };
        V3 c{ pts[3 * T[2] + 0],pts[3 * T[2] + 1],pts[3 * T[2] + 2] };
        V3 d{ pts[3 * T[3] + 0],pts[3 * T[3] + 1],pts[3 * T[3] + 2] };
        V3 s = v_add(v_add(a, b), v_add(c, d));
        return v_scale(s, 0.25);
    }

    V3 tri_normal(V3 a, V3 b, V3 c) { return v_cross(v_sub(b, a), v_sub(c, a)); }

    Plane oriented_tri_plane_for_tet(V3 a, V3 b, V3 c, V3 tetCentroid) {
        V3 n = tri_normal(a, b, c);
        double cst = v_dot(n, a);
        // Orient so tetCentroid is inside: n·tetC <= c
        if (v_dot(n, tetCentroid) > cst) { n = v_scale(n, -1.0); cst = -cst; }
        return { n,cst };
    }

    Plane power_bisector(V3 s0, double w0, V3 s1, double w1) {
        // 2 (s1 - s0) · x <= (||s1||^2 - w1) - (||s0||^2 - w0)
        V3 n = v_sub(s1, s0);
        double c = 0.5 * ((v_dot(s1, s1) - w1) - (v_dot(s0, s0) - w0));
        return { n, c };
    }

    static inline double sd(const Plane& P, V3 x) { return v_dot(P.n, x) - P.c; }
    static inline V3 lerp(V3 a, V3 b, double t) { return { a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, a.z + (b.z - a.z) * t }; }

    rvd::Poly3 make_tet_poly(const double* pts, const uint32_t* T) {
        Poly3 P;
        P.V = {
            {pts[3 * T[0] + 0],pts[3 * T[0] + 1],pts[3 * T[0] + 2]},
            {pts[3 * T[1] + 0],pts[3 * T[1] + 1],pts[3 * T[1] + 2]},
            {pts[3 * T[2] + 0],pts[3 * T[2] + 1],pts[3 * T[2] + 2]},
            {pts[3 * T[3] + 0],pts[3 * T[3] + 1],pts[3 * T[3] + 2]}
        };
        // faces (each opposite the missing vertex), CCW as viewed from outside
        P.F = { {0,2,1}, {0,1,3}, {1,2,3}, {2,0,3} };
        return P;
    }

    std::vector<V3> clip_polygon(const std::vector<V3>& poly, const Plane& P, double eps) {
        std::vector<V3> out; if (poly.empty()) return out;
        auto inside = [&](const V3& q) { return sd(P, q) <= eps; };

        V3 S = poly.back(); bool Sin = inside(S);
        for (const V3& E : poly) {
            bool Ein = inside(E);
            if (Ein) {
                if (!Sin) {
                    V3 d = v_sub(E, S);
                    double denom = v_dot(P.n, d);
                    double t = (P.c - v_dot(P.n, S)) / (denom == 0 ? 1e-30 : denom);
                    t = (std::max)(0.0, (std::min)(1.0, t));
                    out.push_back(lerp(S, E, t));
                }
                out.push_back(E);
            }
            else if (Sin) {
                V3 d = v_sub(E, S);
                double denom = v_dot(P.n, d);
                double t = (P.c - v_dot(P.n, S)) / (denom == 0 ? 1e-30 : denom);
                t = (std::max)(0.0, (std::min)(1.0, t));
                out.push_back(lerp(S, E, t));
            }
            S = E; Sin = Ein;
        }
        // compress near-duplicates
        std::vector<V3> comp;
        for (const auto& q : out) {
            if (comp.empty() || v_norm({ q.x - comp.back().x,q.y - comp.back().y,q.z - comp.back().z }) > eps)
                comp.push_back(q);
        }
        if (comp.size() >= 3 &&
            v_norm({ comp.front().x - comp.back().x,comp.front().y - comp.back().y,comp.front().z - comp.back().z }) <= eps) {
            comp.back() = comp.front(); comp.pop_back();
        }
        return comp;
    }

    Poly3 clip_polyhedron(const Poly3& Pin, const Plane& P, double eps, std::vector<V3>* cutPoly) {
        Poly3 Q;
        std::vector<V3> cutPts;

        for (const auto& face : Pin.F) {
            std::vector<V3> ring; ring.reserve(face.size());
            for (int vi : face) ring.push_back(Pin.V[std::size_t(vi)]);
            auto clipped = clip_polygon(ring, P, eps);
            if (clipped.size() >= 3) {
                int base = (int)Q.V.size();
                for (const auto& v : clipped) Q.V.push_back(v);
                std::vector<int> idx(clipped.size());
                for (size_t k = 0;k < clipped.size();++k) idx[k] = base + (int)k;
                Q.F.push_back(std::move(idx));
            }
            for (const auto& v : clipped) if (std::fabs(sd(P, v)) <= 2 * eps) cutPts.push_back(v);
        }

        // weld Q.V on grid 1/eps
        auto keyOf = [&](V3 a) {
            double s = 1.0 / eps;
            long long ix = llround(a.x * s), iy = llround(a.y * s), iz = llround(a.z * s);
            return std::tuple<long long, long long, long long>(ix, iy, iz);
            };
        std::unordered_map<std::tuple<long long, long long, long long>, int> vm;
        std::vector<V3> WV; WV.reserve(Q.V.size());
        std::vector<int> remap(Q.V.size(), -1);
        for (size_t i = 0;i < Q.V.size();++i) {
            auto K = keyOf(Q.V[i]); auto it = vm.find(K);
            if (it == vm.end()) { int id = (int)WV.size(); vm.emplace(K, id); WV.push_back(Q.V[i]); remap[i] = id; }
            else remap[i] = it->second;
        }
        for (auto& f : Q.F) for (int& vi : f) vi = remap[std::size_t(vi)];
        Q.V.swap(WV);

        auto deg = [&](const std::vector<int>& f)->bool {
            if (f.size() < 3) return true;
            V3 a = Q.V[f[0]], b = Q.V[f[1]], c = Q.V[f[2]];
            return v_norm(tri_normal(a, b, c)) < 1e-18;
            };
        std::vector<std::vector<int>> F2; F2.reserve(Q.F.size());
        for (const auto& f : Q.F) if (!deg(f)) F2.push_back(f);
        Q.F.swap(F2);

        if (cutPoly) {
            std::set<std::tuple<long long, long long, long long>> seen;
            std::vector<V3> Cw;
            for (const auto& q : cutPts) {
                auto K = keyOf(q);
                if (seen.insert(K).second) Cw.push_back(q);
            }
            if (Cw.size() >= 3) {
                V3 n = P.n; double nn = v_norm(n); if (nn < 1e-30) { n = { 0,0,1 }; nn = 1.0; }
                n = v_scale(n, 1.0 / nn);
                V3 t = (std::fabs(n.x) < 0.9) ? V3{ 1,0,0 } : V3{ 0,1,0 };
                V3 u = v_sub(t, v_scale(n, v_dot(n, t)));
                double ul = v_norm(u); if (ul < 1e-30) u = { 1,0,0 }; else u = v_scale(u, 1.0 / ul);
                V3 v = v_cross(n, u);

                V3 C{ 0,0,0 }; for (auto& q : Cw) C = v_add(C, q);
                C = v_scale(C, 1.0 / double(Cw.size()));

                struct Node { V3 q; double ang; };
                std::vector<Node> nodes; nodes.reserve(Cw.size());
                for (auto& q : Cw) {
                    V3 d = v_sub(q, C); double x = v_dot(d, u), y = v_dot(d, v);
                    nodes.push_back({ q,std::atan2(y,x) });
                }
                std::sort(nodes.begin(), nodes.end(), [](const Node& A, const Node& B) { return A.ang < B.ang; });
                std::vector<V3> poly; poly.reserve(nodes.size());
                for (auto& nd : nodes) poly.push_back(nd.q);
                *cutPoly = std::move(poly);
            }
            else {
                cutPoly->clear();
            }
        }

        return Q;
    }

} // namespace rvd
