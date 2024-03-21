

#include "test/test.h"
#include "hex/PrepareGridData.h"
#include "test/osgwindows.h"
#include "test/osgwidgetperformance.h"
#include "test/osggeometry.h"

#include <iostream>
#include <fstream>
#include <filesystem>
#include <OSD_OpenFile.hxx>
#include "test/test_gmsh.h"
#include "test/test_gtest.h"


int main(int argc, char** argv) {
	test_google_start(argc, argv);
    //test_start();

	//auto_mesh_surface();
   /* std::string data_dir = "D:/work/myocc-dev2/log";
    HW_PGD::gen_fdtd_prepare_data(data_dir);*/

    // 假设的命令行参数
  
    //osggeometry(argc, argv);

	/*std::string prepare_data_dir = "C:/Users/hw/Downloads/测试box/file";


	std::locale::global(std::locale(""));
	std::filesystem::path  config_file_path = prepare_data_dir + "./CONFIG.json";
	std::ifstream ifs;
	OSD_OpenStream(ifs, config_file_path.c_str(), std::ios_base::in | std::ios_base::binary);


	if (!ifs.is_open())
	{
		std::cout << "[HW_PGD::read_fdtd_prepare_data] Error opening CONFIG file. prepare data fail!" << std::endl;
		exit(-1);
	}
	else {
		std::cout << "open file!" << std::endl;
	}*/

}