#pragma once
// rvd_meshio_gmsh_v1.h — read Gmsh MSH v1 ASCII ($NOD/$ELM and $Nodes/$Elements)

#include "rvd_types.h"
#include <string>

namespace rvd {

	// Loads tetrahedra (type=4) and triangles (type=2) from a Gmsh v1 ASCII .msh.
	// Accepts both $NOD/$ENDNOD, $ELM/$ENDELM and $Nodes/$EndNodes, $Elements/$EndElements.
	// Extra tags on element lines are ignored by taking the last K ints as node IDs.
	// Returns true on success; on failure returns false and fills `err` (if provided).
	bool load_gmsh_v1_ascii(const std::string& path, TetMesh& out, std::string* err = nullptr);

} // namespace rvd
