#include <iostream>
#include <TopoDS.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>


int main() {
    std::cout << "start ..." << std::endl;

    BRepPrimAPI_MakeBox box(1.0, 1.0, 1.0);
    const TopoDS_Shape& shape = box.Shape();
    BRepTools::Write(shape, "../box.brep");

}