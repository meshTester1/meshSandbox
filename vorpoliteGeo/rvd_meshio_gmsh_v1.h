#pragma once
// rvd_meshio_gmsh_v1.h — read Gmsh MSH v1 ASCII (robust to $NOD/$ELM and $Nodes/$Elements)
#include "rvd_types.h"
#include <string>
namespace rvd {
	// Loads a tetrahedral mesh from a Gmsh v1 ASCII .msh file.
	// - Accepts both legacy blocks: $NOD/$ENDNOD, $ELM/$ENDELM
	//   and the newer $Nodes/$EndNodes, $Elements/$EndElements tokens.
	// - Extracts 3-node triangles (type 2) and 4-node tets (type 4).
	// - Handles extra tags by taking the *last K* integers on each element line as node ids.
	// Returns true on success; on failure, sets err (if provided).
	bool load_gmsh_v1_ascii(const std::string& path, TetMesh& out, std::string* err = nullptr);
} // namespace rvd