#pragma once

#include <string>

#include <gmsh.h>

void auto_mesh_surface() {
	
	std::string file_path = "../test_data/sphere.brep";
    gmsh::initialize();
 
    std::vector<std::pair<int, int> > v;
    try {
        gmsh::model::occ::importShapes(file_path, v);
    }
    catch (...) {
        std::cout << "gmsh handle error! use occ mesher to mesh shape!" << std::endl;
        gmsh::finalize();
        return;
    }
    gmsh::model::occ::synchronize();
    gmsh::model::mesh::generate(2);

    // 获取所有元素
  /*  std::vector<std::size_t> faceNodes;
    gmsh::model::mesh::getElementFaceNodes(elementType, 3, faceNodes);*/


    int elementType = gmsh::model::mesh::getElementType("Triangle", 1);

    std::vector<std::size_t> nodeTags;
    std::vector<double> coord;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodesByElementType(elementType, nodeTags, coord, parametricCoord);
    
    std::vector<std::size_t> elementTags2;
    std::vector<std::size_t> nodeTags2;
    gmsh::model::mesh::getElementsByType(elementType, elementTags2, nodeTags2);

    std::cout << "elementType  = " << elementType << std::endl;

    std::cout << "nodeTags size = " << nodeTags.size() << std::endl;
    std::cout << "coord size = " << coord.size() << std::endl;
    std::cout << "parametricCoord size = " << parametricCoord.size() << std::endl;

    std::cout << "elementTags2 size = " << elementTags2.size() << std::endl;
    std::cout << "nodeTags2 size = " << nodeTags2.size() << std::endl;

    // 打印每个三角形的顶点坐标
    for (int i = 0; i < nodeTags.size(); i += 3) { // 每个三角形有3个顶点
        std::cout << "三角形 " << (i / 3) + 1 << " 的顶点坐标: ";
        for (int j = 0; j < 3; ++j) {
            double x = coord[3 * i + 3 * j];
            double y = coord[3 * i + 3 * j + 1];
            double z = coord[3 * i + 3 * j + 2];
            std::cout << "(" << x << ", " << y << ", " << z << ") ";
        }
        std::cout << std::endl;
    }


    gmsh::finalize();



} 