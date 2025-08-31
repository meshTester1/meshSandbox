//// foamWriter.cpp - OpenFOAM mesh writer implementation
//#include "foamWriter.h"
//#include <fstream>
//#include <iostream>
//#include <unordered_map>
//#include <algorithm>
//
//#ifdef _WIN32
//#  include <direct.h>
//#else
//#  include <sys/stat.h>
//#  include <unistd.h>
//#endif
//
//namespace FoamWriter {
//
//    static void writeFoamHeader(std::ofstream& f, const std::string& cls, const std::string& obj) {
//        f << "FoamFile\n{\n"
//            << "    version     2.0;\n"
//            << "    format      ascii;\n"
//            << "    class       " << cls << ";\n"
//            << "    location    \"polyMesh\";\n"
//            << "    object      " << obj << ";\n"
//            << "}\n\n";
//    }
//
//    static bool writePoints(const std::string& dir, const std::vector<Vec3>& pts) {
//        std::ofstream f(dir + "/points");
//        if (!f) {
//            std::cerr << "Failed to create points file" << std::endl;
//            return false;
//        }
//        writeFoamHeader(f, "vectorField", "points");
//        f << pts.size() << "\n(\n";
//        f.setf(std::ios::scientific);
//        for (const auto& p : pts) {
//            f << " (" << p.x << " " << p.y << " " << p.z << ")\n";
//        }
//        f << ")\n";
//        return true;
//    }
//
//    static bool writeFaces(const std::string& dir, const std::vector<std::vector<int>>& faces) {
//        std::ofstream f(dir + "/faces");
//        if (!f) {
//            std::cerr << "Failed to create faces file" << std::endl;
//            return false;
//        }
//        writeFoamHeader(f, "faceList", "faces");
//        f << faces.size() << "\n(\n";
//        for (const auto& fa : faces) {
//            f << " " << fa.size() << "(";
//            for (size_t i = 0; i < fa.size(); ++i) {
//                f << (i ? " " : "") << fa[i];
//            }
//            f << ")\n";
//        }
//        f << ")\n";
//        return true;
//    }
//
//    static bool writeOwner(const std::string& dir, const std::vector<int>& owner) {
//        std::ofstream f(dir + "/owner");
//        if (!f) {
//            std::cerr << "Failed to create owner file" << std::endl;
//            return false;
//        }
//        writeFoamHeader(f, "labelList", "owner");
//        f << owner.size() << "\n(\n";
//        for (int v : owner) {
//            f << " " << v << "\n";
//        }
//        f << ")\n";
//        return true;
//    }
//
//    static bool writeNeighbour(const std::string& dir, const std::vector<int>& neigh) {
//        std::ofstream f(dir + "/neighbour");
//        if (!f) {
//            std::cerr << "Failed to create neighbour file" << std::endl;
//            return false;
//        }
//        writeFoamHeader(f, "labelList", "neighbour");
//        f << neigh.size() << "\n(\n";
//        for (int v : neigh) {
//            f << " " << v << "\n";
//        }
//        f << ")\n";
//        return true;
//    }
//
//    static bool writeCells(const std::string& dir, const std::vector<std::vector<int>>& cellFaces) {
//        std::ofstream f(dir + "/cells");
//        if (!f) {
//            std::cerr << "Failed to create cells file" << std::endl;
//            return false;
//        }
//        writeFoamHeader(f, "cellList", "cells");
//        f << cellFaces.size() << "\n(\n";
//        for (const auto& fl : cellFaces) {
//            f << " " << fl.size() << "(";
//            for (size_t i = 0; i < fl.size(); ++i) {
//                f << (i ? " " : "") << fl[i];
//            }
//            f << ")\n";
//        }
//        f << ")\n";
//        return true;
//    }
//
//    static bool writeBoundarySinglePatch(const std::string& dir, int nFacesTotal, int nInternalFaces) {
//        std::ofstream f(dir + "/boundary");
//        if (!f) {
//            std::cerr << "Failed to create boundary file" << std::endl;
//            return false;
//        }
//        writeFoamHeader(f, "polyBoundaryMesh", "boundary");
//        int nBoundaryFaces = nFacesTotal - nInternalFaces;
//        f << "1\n(\n"
//            << " boundary\n{\n"
//            << "     type            wall;\n"
//            << "     nFaces          " << nBoundaryFaces << ";\n"
//            << "     startFace       " << nInternalFaces << ";\n"
//            << "}\n)\n";
//        return true;
//    }
//
//    void buildOpenFoamTopology(const std::vector<PolyCell>& cells, PolyMeshOut& out) {
//        out = PolyMeshOut();
//
//        std::unordered_map<Vec3Key, int, Vec3KeyHash> gpmap;
//        auto add_point = [&](const Vec3& v)->int {
//            auto k = quantize(v);
//            auto it = gpmap.find(k);
//            if (it != gpmap.end()) return it->second;
//            int id = (int)out.points.size();
//            out.points.push_back(v);
//            gpmap.emplace(k, id);
//            return id;
//            };
//
//        struct FaceKeyHash {
//            size_t operator()(const std::vector<int>& f) const {
//                size_t h = 1469598103934665603ull;
//                for (int v : f) {
//                    h ^= (size_t)v;
//                    h *= 1099511628211ull;
//                }
//                return h;
//            }
//        };
//
//        struct FaceKeyEq {
//            bool operator()(const std::vector<int>& a, const std::vector<int>& b) const {
//                if (a.size() != b.size()) return false;
//                for (size_t i = 0; i < a.size(); ++i) {
//                    if (a[i] != b[i]) return false;
//                }
//                return true;
//            }
//        };
//
//        struct PendingFace {
//            std::vector<int> verts;
//            int ownerCell = -1;
//        };
//
//        std::unordered_map<std::vector<int>, PendingFace, FaceKeyHash, FaceKeyEq> pending;
//
//        out.cellFaces.resize(cells.size());
//        auto canonical_sorted = [](std::vector<int> f) {
//            std::sort(f.begin(), f.end());
//            return f;
//            };
//
//        std::vector<std::vector<int>> internalFaces;
//        std::vector<int> internalOwner, internalNeigh;
//
//        // Process all cells and their faces
//        for (int c = 0; c < (int)cells.size(); ++c) {
//            const auto& cell = cells[c];
//            if (cell.faces.empty() || cell.verts.empty()) {
//                continue; // Skip empty cells
//            }
//
//            for (const auto& fLocal : cell.faces) {
//                std::vector<int> fGlob;
//                fGlob.reserve(fLocal.size());
//
//                for (int li : fLocal) {
//                    if (li >= 0 && li < (int)cell.verts.size()) {
//                        fGlob.push_back(add_point(cell.verts[li]));
//                    }
//                }
//
//                if ((int)fGlob.size() < 3) continue; // Skip degenerate faces
//
//                auto key = canonical_sorted(fGlob);
//                auto it = pending.find(key);
//
//                if (it == pending.end()) {
//                    // First time seeing this face
//                    pending.emplace(key, PendingFace{ fGlob, c });
//                }
//                else {
//                    // Second time - this is an internal face
//                    internalFaces.push_back(it->second.verts);
//                    internalOwner.push_back(it->second.ownerCell);
//                    internalNeigh.push_back(c);
//
//                    int faceId = (int)internalFaces.size() - 1;
//                    out.cellFaces[it->second.ownerCell].push_back(faceId);
//                    out.cellFaces[c].push_back(faceId);
//
//                    pending.erase(it);
//                }
//            }
//        }
//
//        // Set up internal faces
//        out.nInternalFaces = (int)internalFaces.size();
//        out.faces = std::move(internalFaces);
//        out.owner = std::move(internalOwner);
//        out.neighbour = std::move(internalNeigh);
//
//        // Add boundary faces
//        int nextFaceId = out.nInternalFaces;
//        for (auto& kv : pending) {
//            auto& pf = kv.second;
//            int faceId = nextFaceId++;
//            out.faces.push_back(pf.verts);
//            out.owner.push_back(pf.ownerCell);
//            out.cellFaces[pf.ownerCell].push_back(faceId);
//        }
//
//        std::cout << "OpenFOAM topology: " << out.points.size() << " points, "
//            << out.faces.size() << " faces (" << out.nInternalFaces << " internal, "
//            << (out.faces.size() - out.nInternalFaces) << " boundary)" << std::endl;
//    }
//
//    bool writePolyMesh(const std::string& dir, const PolyMeshOut& mesh) {
//        // Create directory
//#ifdef _WIN32
//        _mkdir(dir.c_str());
//#else
//        mkdir(dir.c_str(), 0755);
//#endif
//
//        std::cout << "Writing OpenFOAM polyMesh files..." << std::endl;
//
//        bool ok = true;
//        ok = ok && writePoints(dir, mesh.points);
//        ok = ok && writeFaces(dir, mesh.faces);
//        ok = ok && writeOwner(dir, mesh.owner);
//        ok = ok && writeNeighbour(dir, mesh.neighbour);
//        ok = ok && writeCells(dir, mesh.cellFaces);
//        ok = ok && writeBoundarySinglePatch(dir, (int)mesh.faces.size(), mesh.nInternalFaces);
//
//        if (!ok) {
//            std::cerr << "Failed writing OpenFOAM polyMesh files" << std::endl;
//            return false;
//        }
//
//        return true;
//    }
//
//}