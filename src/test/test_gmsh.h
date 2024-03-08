#pragma once

#include <string>
//#include <gmsh.h>
//
//void auto_mesh_surface() {
//	
//	std::string file_path = "../test_data/box.step";
//
//	gmsh::initialize();
//	gmsh::option::setNumber("General.Terminal", 1);
//
//	gmsh::vectorpair outDimTags;
//	gmsh::model::occ::importShapes(file_path, outDimTags);
//	gmsh::model::occ::synchronize();
//	gmsh::model::mesh::generate();
//
//	gmsh::write("../test_data/mesh.msh");
//	// ¹Ø±Õ gmsh
//	gmsh::finalize();
//
//
//}