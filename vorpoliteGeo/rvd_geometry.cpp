#include "rvd_geometry.h"

#include <cmath>
#include <vector>
#include <unordered_set>
#include <unordered_map>

namespace rvd {

    // Internal utility functions - completely self-contained
    namespace internal {

        static double clamp_value(double value, double min_val, double max_val) {
            if (value < min_val) return min_val;
            if (value > max_val) return max_val;
            return value;
        }

        static double get_minimum(double a, double b) {
            return (a < b) ? a : b;
        }

        static double get_maximum(double a, double b) {
            return (a > b) ? a : b;
        }

        // Custom round function to avoid llround macro issues
        static long long safe_round(double x) {
            return (long long)(x + (x >= 0.0 ? 0.5 : -0.5));
        }

        struct GridKey {
            long long x, y, z;

            bool operator==(const GridKey& other) const {
                return x == other.x && y == other.y && z == other.z;
            }
        };

        struct GridKeyHash {
            std::size_t operator()(const GridKey& key) const {
                std::size_t h1 = std::hash<long long>{}(key.x);
                std::size_t h2 = std::hash<long long>{}(key.y);
                std::size_t h3 = std::hash<long long>{}(key.z);
                return h1 ^ (h2 << 1) ^ (h3 << 2);
            }
        };
    }

    V3 v_add(V3 a, V3 b) {
        return { a.x + b.x, a.y + b.y, a.z + b.z };
    }

    V3 v_sub(V3 a, V3 b) {
        return { a.x - b.x, a.y - b.y, a.z - b.z };
    }

    V3 v_scale(V3 a, double s) {
        return { a.x * s, a.y * s, a.z * s };
    }

    double v_dot(V3 a, V3 b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    V3 v_cross(V3 a, V3 b) {
        return {
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x
        };
    }

    double v_norm(V3 a) {
        return std::sqrt(v_dot(a, a));
    }

    V3 newell_normal(const std::vector<V3>& P) {
        if (P.size() < 3) return { 0.0, 0.0, 0.0 };

        double nx = 0.0, ny = 0.0, nz = 0.0;
        for (size_t i = 0; i < P.size(); ++i) {
            const V3& a = P[i];
            const V3& b = P[(i + 1) % P.size()];
            nx += (a.y - b.y) * (a.z + b.z);
            ny += (a.z - b.z) * (a.x + b.x);
            nz += (a.x - b.x) * (a.y + b.y);
        }
        return { nx, ny, nz };
    }

    V3 poly_centroid(const std::vector<V3>& P) {
        V3 center = { 0.0, 0.0, 0.0 };
        if (P.empty()) return center;

        for (const auto& point : P) {
            center = v_add(center, point);
        }
        return v_scale(center, 1.0 / double(P.size()));
    }

    V3 tet_centroid(const double* pts, const uint32_t* T) {
        V3 a = { pts[3 * T[0] + 0], pts[3 * T[0] + 1], pts[3 * T[0] + 2] };
        V3 b = { pts[3 * T[1] + 0], pts[3 * T[1] + 1], pts[3 * T[1] + 2] };
        V3 c = { pts[3 * T[2] + 0], pts[3 * T[2] + 1], pts[3 * T[2] + 2] };
        V3 d = { pts[3 * T[3] + 0], pts[3 * T[3] + 1], pts[3 * T[3] + 2] };

        V3 sum = v_add(v_add(a, b), v_add(c, d));
        return v_scale(sum, 0.25);
    }

    V3 tri_normal(V3 a, V3 b, V3 c) {
        return v_cross(v_sub(b, a), v_sub(c, a));
    }

    Plane oriented_tri_plane_for_tet(V3 a, V3 b, V3 c, V3 tetCentroid) {
        V3 normal = tri_normal(a, b, c);
        double constant = v_dot(normal, a);

        // Orient so tetCentroid is inside: n·tetC <= c
        if (v_dot(normal, tetCentroid) > constant) {
            normal = v_scale(normal, -1.0);
            constant = -constant;
        }
        return { normal, constant };
    }

    Plane power_bisector(V3 s0, double w0, V3 s1, double w1) {
        V3 normal = v_sub(s1, s0);
        double constant = 0.5 * ((v_dot(s1, s1) - w1) - (v_dot(s0, s0) - w0));
        return { normal, constant };
    }

    Poly3 make_tet_poly(const double* pts, const uint32_t* T) {
        Poly3 poly;
        poly.V = {
            {pts[3 * T[0] + 0], pts[3 * T[0] + 1], pts[3 * T[0] + 2]},
            {pts[3 * T[1] + 0], pts[3 * T[1] + 1], pts[3 * T[1] + 2]},
            {pts[3 * T[2] + 0], pts[3 * T[2] + 1], pts[3 * T[2] + 2]},
            {pts[3 * T[3] + 0], pts[3 * T[3] + 1], pts[3 * T[3] + 2]}
        };
        // faces (each opposite the missing vertex), CCW as viewed from outside
        poly.F = { {0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3} };
        return poly;
    }

    std::vector<V3> clip_polygon(const std::vector<V3>& poly, const Plane& P, double eps) {
        std::vector<V3> result;
        if (poly.empty()) return result;

        auto is_inside = [&](const V3& q) {
            return (v_dot(P.n, q) - P.c) <= eps;
            };

        auto interpolate = [](V3 a, V3 b, double t) {
            return V3{
                a.x + (b.x - a.x) * t,
                a.y + (b.y - a.y) * t,
                a.z + (b.z - a.z) * t
            };
            };

        V3 prev_vertex = poly.back();
        bool prev_inside = is_inside(prev_vertex);

        for (const V3& curr_vertex : poly) {
            bool curr_inside = is_inside(curr_vertex);

            if (curr_inside) {
                if (!prev_inside) {
                    // Entering - find intersection
                    V3 direction = v_sub(curr_vertex, prev_vertex);
                    double denominator = v_dot(P.n, direction);
                    double t = (P.c - v_dot(P.n, prev_vertex)) /
                        (denominator == 0.0 ? 1e-30 : denominator);
                    t = internal::clamp_value(t, 0.0, 1.0);
                    result.push_back(interpolate(prev_vertex, curr_vertex, t));
                }
                result.push_back(curr_vertex);
            }
            else if (prev_inside) {
                // Exiting - find intersection
                V3 direction = v_sub(curr_vertex, prev_vertex);
                double denominator = v_dot(P.n, direction);
                double t = (P.c - v_dot(P.n, prev_vertex)) /
                    (denominator == 0.0 ? 1e-30 : denominator);
                t = internal::clamp_value(t, 0.0, 1.0);
                result.push_back(interpolate(prev_vertex, curr_vertex, t));
            }

            prev_vertex = curr_vertex;
            prev_inside = curr_inside;
        }

        // Remove near-duplicate vertices
        std::vector<V3> cleaned;
        for (const auto& vertex : result) {
            if (cleaned.empty()) {
                cleaned.push_back(vertex);
            }
            else {
                V3 diff = v_sub(vertex, cleaned.back());
                if (v_norm(diff) > eps) {
                    cleaned.push_back(vertex);
                }
            }
        }

        // Check if first and last are too close
        if (cleaned.size() >= 3) {
            V3 diff = v_sub(cleaned.front(), cleaned.back());
            if (v_norm(diff) <= eps) {
                cleaned.pop_back();
            }
        }

        return cleaned;
    }

    Poly3 clip_polyhedron(const Poly3& Pin, const Plane& P, double eps, std::vector<V3>* cutPoly) {
        Poly3 result;
        std::vector<V3> cut_points;

        // Process each face
        for (const auto& face : Pin.F) {
            std::vector<V3> face_vertices;
            face_vertices.reserve(face.size());

            for (int vertex_index : face) {
                face_vertices.push_back(Pin.V[std::size_t(vertex_index)]);
            }

            auto clipped_face = clip_polygon(face_vertices, P, eps);

            if (clipped_face.size() >= 3) {
                int base_index = (int)result.V.size();

                for (const auto& vertex : clipped_face) {
                    result.V.push_back(vertex);
                }

                std::vector<int> face_indices(clipped_face.size());
                for (size_t k = 0; k < clipped_face.size(); ++k) {
                    face_indices[k] = base_index + (int)k;
                }
                result.F.push_back(std::move(face_indices));
            }

            // Collect potential cut points
            for (const auto& vertex : clipped_face) {
                double distance = std::fabs(v_dot(P.n, vertex) - P.c);
                if (distance <= 2.0 * eps) {
                    cut_points.push_back(vertex);
                }
            }
        }

        // Vertex welding using simple grid
        std::unordered_map<internal::GridKey, int, internal::GridKeyHash> vertex_map;
        std::vector<V3> welded_vertices;
        welded_vertices.reserve(result.V.size());
        std::vector<int> vertex_remap(result.V.size(), -1);

        double grid_scale = 1.0 / eps;

        for (size_t i = 0; i < result.V.size(); ++i) {
            internal::GridKey key = {
                internal::safe_round(result.V[i].x * grid_scale),
                internal::safe_round(result.V[i].y * grid_scale),
                internal::safe_round(result.V[i].z * grid_scale)
            };

            auto it = vertex_map.find(key);
            if (it == vertex_map.end()) {
                int new_id = (int)welded_vertices.size();
                vertex_map[key] = new_id;
                welded_vertices.push_back(result.V[i]);
                vertex_remap[i] = new_id;
            }
            else {
                vertex_remap[i] = it->second;
            }
        }

        // Update face indices
        for (auto& face : result.F) {
            for (int& vertex_index : face) {
                vertex_index = vertex_remap[std::size_t(vertex_index)];
            }
        }
        result.V = std::move(welded_vertices);

        // Remove degenerate faces
        auto is_degenerate = [&](const std::vector<int>& face) -> bool {
            if (face.size() < 3) return true;
            V3 a = result.V[face[0]];
            V3 b = result.V[face[1]];
            V3 c = result.V[face[2]];
            return v_norm(tri_normal(a, b, c)) < 1e-18;
            };

        std::vector<std::vector<int>> valid_faces;
        valid_faces.reserve(result.F.size());

        for (const auto& face : result.F) {
            if (!is_degenerate(face)) {
                valid_faces.push_back(face);
            }
        }
        result.F = std::move(valid_faces);

        // Create cut polygon if requested - using simple vector instead of unordered_set
        if (cutPoly) {
            std::vector<V3> unique_cut_points;
            double grid_scale = 1.0 / eps;

            // Simple deduplication using manual comparison
            for (const auto& point : cut_points) {
                bool already_exists = false;
                for (const auto& existing : unique_cut_points) {
                    V3 diff = v_sub(point, existing);
                    if (v_norm(diff) <= eps) {
                        already_exists = true;
                        break;
                    }
                }
                if (!already_exists) {
                    unique_cut_points.push_back(point);
                }
            }

            if (unique_cut_points.size() >= 3) {
                // Create local coordinate system
                V3 normal = P.n;
                double normal_length = v_norm(normal);
                if (normal_length < 1e-30) {
                    normal = V3{ 0.0, 0.0, 1.0 };
                    normal_length = 1.0;
                }
                normal = v_scale(normal, 1.0 / normal_length);

                V3 tangent = (std::fabs(normal.x) < 0.9) ? V3{ 1.0, 0.0, 0.0 } : V3{ 0.0, 1.0, 0.0 };
                V3 u = v_sub(tangent, v_scale(normal, v_dot(normal, tangent)));
                double u_length = v_norm(u);
                if (u_length < 1e-30) {
                    u = V3{ 1.0, 0.0, 0.0 };
                }
                else {
                    u = v_scale(u, 1.0 / u_length);
                }
                V3 v = v_cross(normal, u);

                // Calculate centroid
                V3 centroid = { 0.0, 0.0, 0.0 };
                for (const auto& point : unique_cut_points) {
                    centroid = v_add(centroid, point);
                }
                centroid = v_scale(centroid, 1.0 / double(unique_cut_points.size()));

                // Sort by angle
                struct AnglePoint {
                    V3 point;
                    double angle;
                };

                std::vector<AnglePoint> angle_points;
                angle_points.reserve(unique_cut_points.size());

                for (const auto& point : unique_cut_points) {
                    V3 offset = v_sub(point, centroid);
                    double x = v_dot(offset, u);
                    double y = v_dot(offset, v);
                    angle_points.push_back({ point, std::atan2(y, x) });
                }

                // Simple insertion sort
                for (size_t i = 1; i < angle_points.size(); ++i) {
                    AnglePoint key = angle_points[i];
                    size_t j = i;
                    while (j > 0 && angle_points[j - 1].angle > key.angle) {
                        angle_points[j] = angle_points[j - 1];
                        --j;
                    }
                    angle_points[j] = key;
                }

                std::vector<V3> sorted_points;
                sorted_points.reserve(angle_points.size());
                for (const auto& ap : angle_points) {
                    sorted_points.push_back(ap.point);
                }
                *cutPoly = std::move(sorted_points);
            }
            else {
                cutPoly->clear();
            }
        }

        return result;
    }

} // namespace rvd
