#ifdef ENABLE_POLT_TEST
#include "PlotTestData.h"
#endif //  ENABLE_POLT_TEST
#include "PrepareGridData.h"

namespace HW_PGD
{

	void gen_fdtd_prepare_data(const std::string& prepare_data_dir)
	{

		configor::json material;
		std::vector<double> XOY_LINES;
		std::vector<double> YOZ_LINES;
		std::vector<double> ZOX_LINES;
		std::unordered_map<int, std::string> objId2matId;
		TopTools_DataMapOfIntegerShape occ_shape_list;

		RTCDevice device = rtcNewDevice(NULL);
		RTCScene scene = rtcNewScene(device);
		std::unordered_map<int, int> geoId2objId;

		bool success = false;
		// read data
		double min_Len; // unit: [mm]
		std::array<int, 2> boundary_level_x, boundary_level_y, boundary_level_z;
		success = read_fdtd_prepare_data(occ_shape_list, scene, device, geoId2objId, objId2matId, material, XOY_LINES, YOZ_LINES, ZOX_LINES, prepare_data_dir, boundary_level_x, boundary_level_y, boundary_level_z, min_Len);
		rtcCommitScene(scene);

		std::unordered_map<std::string, configor::json> materials; // matId => material
		for (auto it = material.begin(); it != material.end(); it++)
		{
			auto mk = it.key();
			auto v = it.value();

			materials[mk] = v;
		}
		// 获取场景中的第一个几何对象

		/*std::cout << "XOY_LINES size = " << XOY_LINES.size() << std::endl;
		std::cout << "YOZ_LINES size = " << YOZ_LINES.size() << std::endl;
		std::cout << "ZOX_LINES size = " << ZOX_LINES.size() << std::endl;
		std::cout << "objid_2_mkid size = " << objid_2_mkid.size() << std::endl;
		std::cout << "occ_shape_list size = " << occ_shape_list.Size() << std::endl;
		std::cout << "geoid_2_objid size = " << geoid_2_objid.size() << std::endl;

		materials.dump(std::cout);*/

		std::vector<FACE> XOY_FACES, YOZ_FACES, ZOX_FACES;

		gen_mesh_faces(XOY_FACES, YOZ_FACES, ZOX_FACES, XOY_LINES, YOZ_LINES, ZOX_LINES);

		std::cout << "XOY_FACES size = " << XOY_FACES.size() << std::endl;
		std::cout << "YOZ_FACES size = " << YOZ_FACES.size() << std::endl;
		std::cout << "ZOX_FACES size = " << ZOX_FACES.size() << std::endl;

		handle_faces_material(occ_shape_list, scene, XOY_FACES, YOZ_FACES, ZOX_FACES, geoId2objId, XOY_LINES, YOZ_LINES, ZOX_LINES, false);
		std::cout << "handle_faces_material  " << std::endl;
#ifdef ENABLE_POLT_TEST

		plot::plot_mesh_faces2(XOY_FACES, YOZ_FACES, ZOX_FACES, XOY_LINES, YOZ_LINES, ZOX_LINES, 0);

#endif //  ENABLE_POLT_TEST
		Mesh _mesh;
		gen_mesh_info(_mesh, XOY_LINES.size(), YOZ_LINES.size(), ZOX_LINES.size(), XOY_FACES, YOZ_FACES, ZOX_FACES, boundary_level_x, boundary_level_y, boundary_level_z, (double)min_Len * 1e-3);
		std::unordered_map<int, int> objId2meshmatlId = gen_objId2meshmatlId(materials, objId2matId, _mesh);
		std::vector<int> PEC_Obj_id = gen_PEC_obj_id(objId2matId, materials);

#pragma omp parallel sections
		{
#pragma omp section
			{
				refine_grid(false, occ_shape_list, scene, geoId2objId, objId2meshmatlId, PEC_Obj_id, XOY_LINES, YOZ_LINES, ZOX_LINES, XOY_FACES, 1);
				std::cout << "XOY_FACES refine_grid" << std::endl;
			}
#pragma omp section
			{
				refine_grid(false, occ_shape_list, scene, geoId2objId, objId2meshmatlId, PEC_Obj_id, XOY_LINES, YOZ_LINES, ZOX_LINES, YOZ_FACES, 2);
				std::cout << "YOZ_FACES refine_grid" << std::endl;
			}
#pragma omp section
			{
				refine_grid(false, occ_shape_list, scene, geoId2objId, objId2meshmatlId, PEC_Obj_id, XOY_LINES, YOZ_LINES, ZOX_LINES, ZOX_FACES, 3);
				std::cout << "ZOX_FACES refine_grid" << std::endl;
			}
		}
		//refine_grid(false, occ_shape_list, scene, geoId2objId, objId2meshmatlId, PEC_Obj_id, XOY_LINES, YOZ_LINES, ZOX_LINES, XOY_FACES, 1);
		//std::cout << "XOY_FACES refine_grid" << std::endl;
		//refine_grid(false, occ_shape_list, scene, geoId2objId, objId2meshmatlId, PEC_Obj_id, XOY_LINES, YOZ_LINES, ZOX_LINES, YOZ_FACES, 2);
		//std::cout << "YOZ_FACES refine_grid" << std::endl;
		//refine_grid(false, occ_shape_list, scene, geoId2objId, objId2meshmatlId, PEC_Obj_id, XOY_LINES, YOZ_LINES, ZOX_LINES, ZOX_FACES, 3);
		//std::cout << "ZOX_FACES refine_grid" << std::endl;

		gen_face_material(false, objId2meshmatlId, XOY_FACES, _mesh);
		gen_face_material(false, objId2meshmatlId, YOZ_FACES, _mesh);
		gen_face_material(false, objId2meshmatlId, ZOX_FACES, _mesh);

		//////////////////////////unreal lines/////////////////////////////
		std::vector<double> UNREAL_XOY_LINES, UNREAL_YOZ_LINES, UNREAL_ZOX_LINES;
		UNREAL_XOY_LINES = gen_unreal_lines(XOY_LINES);
		UNREAL_YOZ_LINES = gen_unreal_lines(YOZ_LINES);
		UNREAL_ZOX_LINES = gen_unreal_lines(ZOX_LINES);

		std::vector<FACE> UNREAL_XOY_FACES, UNREAL_YOZ_FACES, UNREAL_ZOX_FACES;
		gen_unreal_mesh_faces(UNREAL_XOY_FACES, UNREAL_YOZ_FACES, UNREAL_ZOX_FACES, UNREAL_XOY_LINES, UNREAL_YOZ_LINES, UNREAL_ZOX_LINES);

		std::cout << "XOY_FACES size = " << UNREAL_XOY_FACES.size() << std::endl;
		std::cout << "YOZ_FACES size = " << UNREAL_YOZ_FACES.size() << std::endl;
		std::cout << "ZOX_FACES size = " << UNREAL_ZOX_FACES.size() << std::endl;

		handle_faces_material(occ_shape_list, scene, UNREAL_XOY_FACES, UNREAL_YOZ_FACES, UNREAL_ZOX_FACES, geoId2objId, UNREAL_XOY_LINES, UNREAL_YOZ_LINES, UNREAL_ZOX_LINES, true);
		std::cout << "handle_faces_material  " << std::endl;
#ifdef ENABLE_POLT_TEST

		//plot::plot_mesh_faces2(UNREAL_XOY_FACES, UNREAL_YOZ_FACES, UNREAL_ZOX_FACES, XOY_LINES, YOZ_LINES, ZOX_LINES, -1);

#endif //  ENABLE_POLT_TEST

#pragma omp parallel sections
		{
#pragma omp section
			{
				refine_grid(true, occ_shape_list, scene, geoId2objId, objId2meshmatlId, PEC_Obj_id, UNREAL_XOY_LINES, UNREAL_YOZ_LINES, UNREAL_ZOX_LINES, UNREAL_XOY_FACES, 1);
				std::cout << "XOY_FACES refine_grid" << std::endl;
			}
#pragma omp section
			{
				refine_grid(true, occ_shape_list, scene, geoId2objId, objId2meshmatlId, PEC_Obj_id, UNREAL_XOY_LINES, UNREAL_YOZ_LINES, UNREAL_ZOX_LINES, UNREAL_YOZ_FACES, 2);
				std::cout << "YOZ_FACES refine_grid" << std::endl;
			}
#pragma omp section
			{
				refine_grid(true, occ_shape_list, scene, geoId2objId, objId2meshmatlId, PEC_Obj_id, UNREAL_XOY_LINES, UNREAL_YOZ_LINES, UNREAL_ZOX_LINES, UNREAL_ZOX_FACES, 3);
				std::cout << "ZOX_FACES refine_grid" << std::endl;
			}
		}
		// refine_grid(true, occ_shape_list, scene, geoId2objId, objId2meshmatlId, PEC_Obj_id, UNREAL_XOY_LINES, UNREAL_YOZ_LINES, UNREAL_ZOX_LINES, UNREAL_XOY_FACES, 1);
		// std::cout << "XOY_FACES refine_grid" << std::endl;
		// refine_grid(true, occ_shape_list, scene, geoId2objId, objId2meshmatlId, PEC_Obj_id, UNREAL_XOY_LINES, UNREAL_YOZ_LINES, UNREAL_ZOX_LINES, UNREAL_YOZ_FACES, 2);
		// std::cout << "YOZ_FACES refine_grid" << std::endl;
		// refine_grid(true, occ_shape_list, scene, geoId2objId, objId2meshmatlId, PEC_Obj_id, UNREAL_XOY_LINES, UNREAL_YOZ_LINES, UNREAL_ZOX_LINES, UNREAL_ZOX_FACES, 3);
		// std::cout << "ZOX_FACES refine_grid" << std::endl;

		gen_face_material(true, objId2meshmatlId, UNREAL_XOY_FACES, _mesh);
		gen_face_material(true, objId2meshmatlId, UNREAL_YOZ_FACES, _mesh);
		gen_face_material(true, objId2meshmatlId, UNREAL_ZOX_FACES, _mesh);

		int MPI_ip = 0;
		PARALLEL_INF _PARALLEL_INF;
		Field _field;

		_PARALLEL_INF.MPI_Lx = 1;
		_PARALLEL_INF.MPI_Ly = 1;
		_PARALLEL_INF.MPI_Lz = 1;
		_PARALLEL_INF.N_x = _mesh.NX;
		_PARALLEL_INF.N_y = _mesh.NY;
		_PARALLEL_INF.N_z = _mesh.NZ;
		_PARALLEL_INF.N_start_x = 0;
		_PARALLEL_INF.N_start_y = 0;
		_PARALLEL_INF.N_start_z = 0;

		gen_field_index_parallel(_PARALLEL_INF, XOY_FACES, YOZ_FACES, ZOX_FACES, UNREAL_XOY_FACES, UNREAL_YOZ_FACES, UNREAL_ZOX_FACES, _mesh, _field, MPI_ip);

		// std::string filepath = prepare_data_dir + "fdtd_cad_data.bin";
		// write_bin(_mesh, _field, filepath);
		return;
	}

	// zzz

	bool read_fdtd_prepare_data(TopTools_DataMapOfIntegerShape& occ_shape_list, RTCScene& scene, RTCDevice& device, std::unordered_map<int, int>& geoid_2_objid, std::unordered_map<int, std::string>& objid_2_mkid, configor::json& materials, std::vector<double>& XOY_LINES, std::vector<double>& YOZ_LINES, std::vector<double>& ZOX_LINES, const std::string& prepare_data_dir, std::array<int, 2>& boundary_layer_x, std::array<int, 2>& boundary_layer_y, std::array<int, 2>& boundary_layer_z, double& min_Len)
	{
		bool success = false;

		std::string config_file_path = prepare_data_dir + "./CONFIG.json";

		std::ifstream ifs(config_file_path);

		if (!ifs.is_open())
		{
			std::cout << "[HW_PGD::read_fdtd_prepare_data] Error opening CONFIG file. prepare data fail!" << std::endl;
			return false;
		}

		// 将文件内容读取到字符串
		std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
		configor::json jsonData;
		jsonData = configor::json::parse(content);

		// get material
		if (jsonData.find("materials") != jsonData.end())
		{
			materials = jsonData["materials"];
		}
		else
		{
			std::cout << "[HW_PGD::read_fdtd_prepare_data] can not find the materials, prepare data fail!" << std::endl;
			return false;
		}
		// get XOY_LINES
		if (jsonData.find("xoy_lines") != jsonData.end())
		{
			auto xoy_arr = jsonData["xoy_lines"];
			for (size_t i = 0; i < xoy_arr.size(); i++)
			{
				double x = xoy_arr[i];
				XOY_LINES.push_back(x);
			}
		}
		else
		{
			std::cout << "[HW_PGD::read_fdtd_prepare_data] can not find the xoy_lines, prepare data fail!" << std::endl;
			return false;
		}
		// get YOZ_LINES
		if (jsonData.find("yoz_lines") != jsonData.end())
		{
			auto yoz_arr = jsonData["yoz_lines"];
			for (size_t i = 0; i < yoz_arr.size(); i++)
			{
				double x = yoz_arr[i];
				YOZ_LINES.push_back(x);
			}
		}
		else
		{
			std::cout << "[HW_PGD::read_fdtd_prepare_data] can not find the yoz_lines, prepare data fail!" << std::endl;
			return false;
		}
		// get ZOX_LINES
		if (jsonData.find("zox_lines") != jsonData.end())
		{
			auto zox_arr = jsonData["zox_lines"];
			for (size_t i = 0; i < zox_arr.size(); i++)
			{
				double x = zox_arr[i];
				ZOX_LINES.push_back(x);
			}
		}
		else
		{
			std::cout << "[HW_PGD::read_fdtd_prepare_data] can not find the zox_lines, prepare data fail!" << std::endl;
			return false;
		}

		// get airbox
		if (jsonData.find("airbox") != jsonData.end())
		{
			auto airbox = jsonData["airbox"];
			int id = airbox["id"];
			std::string mk = airbox["mk"];
			objid_2_mkid[id] = mk;
		}
		else
		{
			std::cout << "[HW_PGD::read_fdtd_prepare_data] can not find the airbox, prepare data fail!" << std::endl;
			return false;
		}

		// get models
		if (jsonData.find("models") != jsonData.end())
		{
			auto models = jsonData["models"];

			for (size_t i = 0; i < models.size(); i++)
			{
				auto model = models[i];
				std::string file_name = model["filename"];
				std::string mk = model["mk"];
				int id = model["id"];
				int type = model["type"];
				std::string file_path = prepare_data_dir + "./" + file_name;

				objid_2_mkid[id] = mk;

				if (type == 0)
				{ // occ shape
					TopoDS_Shape shape;
					BRep_Builder builder;

					BRepTools::Read(shape, file_path.c_str(), builder);
					occ_shape_list.Bind(id, shape);
				}
				else
				{ // mesh shape
					RTCGeometry geom = read_mesh_data(device, file_path);

					rtcCommitGeometry(geom);

					unsigned int geo_id = rtcAttachGeometry(scene, geom);

					geoid_2_objid[geo_id] = id;
				}
			}
		}
		else
		{
			std::cout << "[HW_PGD::read_fdtd_prepare_data] can not find the models, prepare data fail!" << std::endl;
			return false;
		}

		// get boundary layers x
		if (jsonData.find("boundary_level_x") != jsonData.end())
		{
			auto boundary_level_x = jsonData["boundary_level_x"];
			boundary_layer_x[0] = (double)boundary_level_x[0];
			boundary_layer_x[1] = (double)boundary_level_x[1];
		}
		else
		{
			std::cout << "[HW_PGD::read_fdtd_prepare_data] can not find the boundary_level_x, prepare data fail!" << std::endl;
			return false;
		}
		// get boundary layers y
		if (jsonData.find("boundary_level_y") != jsonData.end())
		{
			auto boundary_level_y = jsonData["boundary_level_y"];
			boundary_layer_y[0] = (double)boundary_level_y[0];
			boundary_layer_y[1] = (double)boundary_level_y[1];
		}
		else
		{
			std::cout << "[HW_PGD::read_fdtd_prepare_data] can not find the boundary_level_y, prepare data fail!" << std::endl;
			return false;
		}
		// get boundary layers z
		if (jsonData.find("boundary_level_z") != jsonData.end())
		{
			auto boundary_level_z = jsonData["boundary_level_z"];
			boundary_layer_z[0] = (double)boundary_level_z[0];
			boundary_layer_z[1] = (double)boundary_level_z[1];
		}
		else
		{
			std::cout << "[HW_PGD::read_fdtd_prepare_data] can not find the boundary_level_z, prepare data fail!" << std::endl;
			return false;
		}

		// get cell min_Len
		if (jsonData.find("min_Len") != jsonData.end())
		{
			min_Len = (double)jsonData["min_Len"];
		}
		else
		{
			std::cout << "[HW_PGD::read_fdtd_prepare_data] can not find the min_Len, prepare data fail!" << std::endl;
			return false;
		}
		return true;
	}

	RTCGeometry read_mesh_data(RTCDevice& device, const std::string& file_path)
	{

		RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
		std::ifstream file(file_path);

		if (!file.is_open())
		{
			std::cout << "[PGD::read_mesh_data] can not open mesh file!" << std::endl;
			return false;
		}

		double Xmin = INFINITY, Ymin = INFINITY, Zmin = INFINITY, Xmax = -INFINITY, Ymax = -INFINITY, Zmax = -INFINITY;
		std::vector<float> V;
		std::vector<int> F;
		size_t vertices_count;
		size_t indices_count;

		std::string line;
		while (std::getline(file, line))
		{

			std::istringstream iss(line);

			std::string first_word;
			iss >> first_word;

			if (first_word == "VERTICES")
			{

				if (iss >> vertices_count)
				{
					for (int i = 0; i < vertices_count; i++)
					{
						double x, y, z;
						if (std::getline(file, line) && sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z) == 3)
						{

							V.push_back(x);
							V.push_back(y);
							V.push_back(z);

							if (x < Xmin)
							{
								Xmin = x;
							}

							if (y < Ymin)
							{
								Ymin = y;
							}

							if (z < Zmin)
							{
								Zmin = z;
							}

							if (x > Xmax)
							{
								Xmax = x;
							}

							if (y > Ymax)
							{
								Ymax = y;
							}

							if (z > Zmax)
							{
								Zmax = z;
							}
						}
						else
						{
							std::cout << "Error reading vertex data." << std::endl;
							return false;
						}
					}
				}
				else
				{
					std::cerr << "[PGD::read_mesh_data] Error reading vertices count." << std::endl;
				}
			}
			else if (first_word == "FACES")
			{

				if (iss >> indices_count)
				{
					for (int i = 0; i < indices_count; i++)
					{
						size_t index0, index1, index2;

						if (std::getline(file, line) && sscanf(line.c_str(), "%d %d %d", &index0, &index1, &index2) == 3)
						{

							F.push_back(index0);
							F.push_back(index1);
							F.push_back(index2);
						}
						else
						{
							std::cout << "Error reading indices data." << std::endl;
							return false;
						}
					}
				}
				else
				{
					std::cerr << "[PGD::read_mesh_data] Error reading indices count." << std::endl;
				}
			}
		}

		unsigned int* indices = (unsigned int*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(unsigned int) * 3, indices_count);
		float* vertices = (float*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(float) * 3, vertices_count);


		for (int i = 0; i < V.size(); i += 3)
		{
			vertices[i] = V[i];
			vertices[i + 1] = V[i + 1];
			vertices[i + 2] = V[i + 2];
		}

		for (int i = 0; i < F.size(); i += 3)
		{
			indices[i] = F[i];
			indices[i + 1] = F[i + 1];
			indices[i + 2] = F[i + 2];
		}

		HW_PGD::GEO_USER_DATA ud;

		ud.Xmin = Xmin;
		ud.Xmax = Xmax;
		ud.Ymin = Ymin;
		ud.Ymax = Ymax;
		ud.Zmin = Zmin;
		ud.Zmax = Zmax;

		rtcSetGeometryUserData(geom, &ud);

		file.close();

		return geom;
	}

	std::vector<double> gen_unreal_lines(const std::vector<double>& origin_lines)
	{
		std::vector<double> UNREAL_LINES;

		for (int i = 0; i < origin_lines.size(); i++)
		{
			double center;
			if (i == 0)
			{
				center = (origin_lines[i] + origin_lines[i + 1]) / 2;
				double diff = center - origin_lines[i];
				double fistLine = origin_lines[i] - diff;
				UNREAL_LINES.push_back(fistLine);
				UNREAL_LINES.push_back(center);
			}
			else if (i == origin_lines.size() - 1)
			{
				center = (origin_lines[i] + origin_lines[i - 1]) / 2;

				double diff = origin_lines[i] - center;
				double lastLine = origin_lines[i] + diff;
				UNREAL_LINES.push_back(lastLine);
			}
			else
			{
				center = (origin_lines[i] + origin_lines[i + 1]) / 2;

				UNREAL_LINES.push_back(center);
			}
		}

		return UNREAL_LINES;
	}

	bool gen_mesh_faces(std::vector<FACE>& XOY_FACES, std::vector<FACE>& YOZ_FACES, std::vector<FACE>& ZOX_FACES, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES)
	{

		XOY_FACES.clear();
		YOZ_FACES.clear();
		ZOX_FACES.clear();

#pragma omp parallel sections
		{
#pragma omp section
			{
				// calc all faces that parallel to XOY
				for (int i = 0; i < ZOX_LINES.size(); i++)
				{

					for (int j = 0; j < (YOZ_LINES.size() - 1); j++)
					{

						for (int k = 0; k < (XOY_LINES.size() - 1); k++)
						{

							FACE f;

							f.line_indices[0] = k;
							f.line_indices[1] = j;
							f.line_indices[2] = i;

							double h = Abs(XOY_LINES[k + 1] - XOY_LINES[k]);
							double v = Abs(YOZ_LINES[j + 1] - YOZ_LINES[j]);
							f.h = h;
							f.v = v;

							XOY_FACES.push_back(f);
						}
					}
				}
			}
			// calc all faces that parallel to YOZ
#pragma omp section
			{
				for (int i = 0; i < XOY_LINES.size(); i++)
				{

					for (int j = 0; j < (ZOX_LINES.size() - 1); j++)
					{

						for (int k = 0; k < (YOZ_LINES.size() - 1); k++)
						{

							FACE f;

							f.line_indices[0] = k;
							f.line_indices[1] = j;
							f.line_indices[2] = i;

							double h = Abs(YOZ_LINES[k + 1] - YOZ_LINES[k]);
							double v = Abs(ZOX_LINES[j + 1] - ZOX_LINES[j]);
							f.h = h;
							f.v = v;

							YOZ_FACES.push_back(f);
						}
					}
				}
			}
			// calc all faces that parallel to ZOX
#pragma omp section
			{
				for (int i = 0; i < YOZ_LINES.size(); i++)
				{

					for (int j = 0; j < (XOY_LINES.size() - 1); j++)
					{

						for (int k = 0; k < (ZOX_LINES.size() - 1); k++)
						{

							FACE f;

							f.line_indices[0] = k;
							f.line_indices[1] = j;
							f.line_indices[2] = i;

							double h = Abs(ZOX_LINES[k + 1] - ZOX_LINES[k]);
							double v = Abs(XOY_LINES[j + 1] - XOY_LINES[j]);
							f.h = h;
							f.v = v;

							ZOX_FACES.push_back(f);
						}
					}
				}
			}
		}

		return true;
	}

	bool gen_unreal_mesh_faces(std::vector<FACE>& XOY_FACES, std::vector<FACE>& YOZ_FACES, std::vector<FACE>& ZOX_FACES, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES)
	{

		XOY_FACES.clear();
		YOZ_FACES.clear();
		ZOX_FACES.clear();

#pragma omp parallel sections
		{
#pragma omp section
			{
				// calc all faces that parallel to XOY
				for (int i = 1; i < ZOX_LINES.size() - 1; i++)
				{

					for (int j = 0; j < (YOZ_LINES.size() - 1); j++)
					{

						for (int k = 0; k < (XOY_LINES.size() - 1); k++)
						{

							FACE f;

							f.line_indices[0] = k;
							f.line_indices[1] = j;
							f.line_indices[2] = i;

							double h = Abs(XOY_LINES[k + 1] - XOY_LINES[k]);
							double v = Abs(YOZ_LINES[j + 1] - YOZ_LINES[j]);
							f.h = h;
							f.v = v;

							XOY_FACES.push_back(f);
						}
					}
				}
			}
			// calc all faces that parallel to YOZ
#pragma omp section
			{
				for (int i = 1; i < XOY_LINES.size() - 1; i++)
				{

					for (int j = 0; j < (ZOX_LINES.size() - 1); j++)
					{

						for (int k = 0; k < (YOZ_LINES.size() - 1); k++)
						{

							FACE f;

							f.line_indices[0] = k;
							f.line_indices[1] = j;
							f.line_indices[2] = i;

							double h = Abs(YOZ_LINES[k + 1] - YOZ_LINES[k]);
							double v = Abs(ZOX_LINES[j + 1] - ZOX_LINES[j]);
							f.h = h;
							f.v = v;

							YOZ_FACES.push_back(f);
						}
					}
				}
			}
			// calc all faces that parallel to ZOX
#pragma omp section
			{
				for (int i = 1; i < YOZ_LINES.size() - 1; i++)
				{

					for (int j = 0; j < (XOY_LINES.size() - 1); j++)
					{

						for (int k = 0; k < (ZOX_LINES.size() - 1); k++)
						{

							FACE f;

							f.line_indices[0] = k;
							f.line_indices[1] = j;
							f.line_indices[2] = i;

							double h = Abs(ZOX_LINES[k + 1] - ZOX_LINES[k]);
							double v = Abs(XOY_LINES[j + 1] - XOY_LINES[j]);
							f.h = h;
							f.v = v;

							ZOX_FACES.push_back(f);
						}
					}
				}
			}
		}

		return true;
	}

	gp_Pnt get_xoy_face_center(const FACE& face, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES)
	{

		size_t line0 = face.line_indices[0]; // x line
		size_t line1 = face.line_indices[1]; // y line
		size_t line2 = face.line_indices[2]; // z line

		double x = XOY_LINES[line0] + (XOY_LINES[line0 + 1] - XOY_LINES[line0]) / 2;
		double y = YOZ_LINES[line1] + (YOZ_LINES[line1 + 1] - YOZ_LINES[line1]) / 2;
		double z = ZOX_LINES[line2];

		return gp_Pnt(x, y, z);
	}

	gp_Pnt get_yoz_face_center(const FACE& face, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES)
	{

		size_t line0 = face.line_indices[0]; // y line
		size_t line1 = face.line_indices[1]; // z line
		size_t line2 = face.line_indices[2]; // x line

		double y = YOZ_LINES[line0] + (YOZ_LINES[line0 + 1] - YOZ_LINES[line0]) / 2;
		double z = ZOX_LINES[line1] + (ZOX_LINES[line1 + 1] - ZOX_LINES[line1]) / 2;
		double x = XOY_LINES[line2];

		return gp_Pnt(x, y, z);
	}

	gp_Pnt get_zox_face_center(const FACE& face, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES)
	{

		size_t line0 = face.line_indices[0]; // y line
		size_t line1 = face.line_indices[1]; // z line
		size_t line2 = face.line_indices[2]; // x line

		double z = ZOX_LINES[line0] + (ZOX_LINES[line0 + 1] - ZOX_LINES[line0]) / 2;
		double x = XOY_LINES[line1] + (XOY_LINES[line1 + 1] - XOY_LINES[line1]) / 2;
		double y = YOZ_LINES[line2];

		return gp_Pnt(x, y, z);
	}

	std::vector<gp_Pnt> get_occ_shape_ray_intersection_points(const TopoDS_Shape& shape, const gp_Pnt& originPnt, const gp_Dir& dir)
	{

		std::vector<gp_Pnt> pnts;

		// create ray
		Handle(Geom_Curve) ray = new Geom_Line(originPnt, dir);

		GeomAPI_IntCS intersector;

		TopExp_Explorer exp(shape, TopAbs_FACE); // iterator
		while (exp.More())
		{
			const TopoDS_Face& face = TopoDS::Face(exp.Current());
			const Handle(Geom_Surface)& surface = BRep_Tool::Surface(face);
			intersector.Perform(ray, surface);

			if (intersector.IsDone() && intersector.NbPoints() > 0)
			{
				for (int i = 1; i <= intersector.NbPoints(); i++)
				{
					const gp_Pnt& pnt = intersector.Point(i); // 获取交点的坐标
					pnts.push_back(pnt);
				}
			}

			exp.Next();
		}

		return pnts;
	}

	std::unordered_map<int, std::vector<gp_Pnt>> get_occ_shape_ray_intersection_points(const TopTools_DataMapOfIntegerShape& occ_shape_list, const gp_Pnt& originPnt, const gp_Dir& dir)
	{
		double tol = 4 * 1e-6;

		std::unordered_map<int, std::vector<gp_Pnt>> intersectPoints;

		gp_Lin ray(originPnt, dir);

		Handle(Geom_Line) aGeomLine = new Geom_Line(ray);
		aGeomLine->SetLocation(originPnt);
		aGeomLine->SetDirection(dir);

		gp_XYZ c = dir.XYZ();

		TopTools_DataMapIteratorOfDataMapOfIntegerShape iter(occ_shape_list);
		// 遍历数据映射中的所有元素
		for (; iter.More(); iter.Next())
		{

			int objId = iter.Key();
			const TopoDS_Shape& shape = iter.Value();
			// ray geo obj
			std::vector<gp_Pnt> pnts = get_occ_shape_ray_intersection_points(shape, originPnt, dir);
			std::vector<gp_Pnt> resultPnts;

			if (pnts.size() > 1)
			{

				// sort
				std::sort(pnts.begin(), pnts.end(), [c](const gp_Pnt& p1, const gp_Pnt& p2)
				{
					gp_XYZ c1 = p1.XYZ();
					gp_XYZ c2 = p2.XYZ();
					return (c1.Dot(c) < c2.Dot(c)); });

				// del the repeat value
				for (int i = 0; i < pnts.size() - 1; i++)
				{

					if (i == 0)
					{
						resultPnts.push_back(pnts[0]);
					}
					gp_XYZ c1 = resultPnts[resultPnts.size() - 1].XYZ();
					gp_XYZ c2 = pnts[i + 1].XYZ();
					if (Abs(c1.Dot(c) - c2.Dot(c)) >= tol)
					{

						// pnts[i + 1].DumpJson(std::cout);
						resultPnts.push_back(pnts[i + 1]);
					}
				}
			}

			intersectPoints[objId] = resultPnts;
		}

		return intersectPoints;
	}

	std::unordered_map<int, std::vector<gp_Pnt>> get_mesh_ray_intersection_points(const RTCScene& scene, const gp_Pnt& originPnt, const gp_Dir& dir)
	{
		double tol = 4 * 1e-6;

		std::unordered_map<int, std::vector<gp_Pnt>> intersectPoints;

		// 创建一个 RTCRayHit 对象，表示一条射线和一个交点
		RTCRayHit rayhit;

		rayhit.ray.org_x = originPnt.X();
		rayhit.ray.org_y = originPnt.Y();
		rayhit.ray.org_z = originPnt.Z();
		rayhit.ray.dir_x = dir.X();
		rayhit.ray.dir_y = dir.Y();
		rayhit.ray.dir_z = dir.Z();
		rayhit.ray.tnear = tol;									  // ray min distance
		rayhit.ray.tfar = std::numeric_limits<float>::infinity(); // ray max distance
		rayhit.ray.mask = 0xFFFFFFFF;							  // ray code
		rayhit.ray.flags = 0;									  // ray flat

		rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID; // geo ID
		rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID; // tri ID
		rayhit.hit.Ng_x = 0.0f;						 // intersecting point normal x
		rayhit.hit.Ng_y = 0.0f;						 // intersecting point normal y
		rayhit.hit.Ng_z = 0.0f;						 // intersecting point normal z

		rtcIntersect1(scene, &rayhit);

		while (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
		{ // ray next
			int geoId = rayhit.hit.geomID;
			// calc intersecting point
			double hitX = rayhit.ray.org_x + rayhit.ray.tfar * rayhit.ray.dir_x;
			double hitY = rayhit.ray.org_y + rayhit.ray.tfar * rayhit.ray.dir_y;
			double hitZ = rayhit.ray.org_z + rayhit.ray.tfar * rayhit.ray.dir_z;

			if (intersectPoints.find(geoId) == intersectPoints.end())
			{
				std::vector<gp_Pnt> pnts;
				pnts.push_back(gp_Pnt(hitX, hitY, hitZ));
				intersectPoints[geoId] = pnts;
			}
			else
			{
				intersectPoints[geoId].push_back(gp_Pnt(hitX, hitY, hitZ));
			}

			rayhit.ray.org_x = hitX + tol * dir.X();
			rayhit.ray.org_y = hitY + tol * dir.Y();
			rayhit.ray.org_z = hitZ + tol * dir.Z();
			rayhit.ray.dir_x = dir.X();
			rayhit.ray.dir_y = dir.Y();
			rayhit.ray.dir_z = dir.Z();
			rayhit.ray.tnear = tol;									  // ray min distance
			rayhit.ray.tfar = std::numeric_limits<float>::infinity(); // ray max distance
			rayhit.ray.mask = 0xFFFFFFFF;							  // ray code
			rayhit.ray.flags = 0;									  // ray flat

			rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID; // geo ID
			rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID; // tri ID
			rayhit.hit.Ng_x = 0.0f;						 // intersecting point normal x
			rayhit.hit.Ng_y = 0.0f;						 // intersecting point normal y
			rayhit.hit.Ng_z = 0.0f;						 // intersecting point normal z

			rtcIntersect1(scene, &rayhit);
		}

		return intersectPoints;
	}

	void set_face_belong_to_obj(FACE& face, int obj_id)
	{

		if (obj_id != -1)
		{
			face.Obj_id.push_back(obj_id);
			if (face.Obj_id.size() > 1)
			{
				face.flag = 0;
			}
			else
			{
				face.flag = -1;
			}
		}
		else // AirBox obj_id = -1
		{
			if (face.Obj_id.empty())
			{ // if had been set other obj_id, can not set air_box id
				face.Obj_id = std::vector<int>{ obj_id };
				face.flag = 1;
			}
		}
	}

	bool handle_faces_material(const TopTools_DataMapOfIntegerShape& occ_shape_list, const RTCScene& scene, std::vector<FACE>& XOY_FACES, std::vector<FACE>& YOZ_FACES, std::vector<FACE>& ZOX_FACES, const std::unordered_map<int, int>& geoId2objId, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES, bool is_unreal)
	{

		// handle XOY face material
#pragma omp parallel sections
		{
			// 1. handle XOY  ray
#pragma omp section
			{
				size_t ZOX_COUNT = (is_unreal == false) ? ZOX_LINES.size() : (ZOX_LINES.size() - 2);
				size_t ONE_XOY_FACES_COUNT = (XOY_LINES.size() - 1) * (YOZ_LINES.size() - 1);

				for (int i = 0; i < ONE_XOY_FACES_COUNT; i++)
				{

					const FACE& face = XOY_FACES[i];
					const gp_Pnt& face_center = get_xoy_face_center(face, XOY_LINES, YOZ_LINES, ZOX_LINES);

					gp_Pnt startPoint(face_center.X(), face_center.Y(), face_center.Z() - 1);

					gp_Dir direction(0, 0, 1);
					///////////////////////////////////////////// ray occ obj

					std::unordered_map<int, std::vector<gp_Pnt>> objId2Points = get_occ_shape_ray_intersection_points(occ_shape_list, startPoint, direction);
					const std::unordered_map<int, std::vector<gp_Pnt>>& geoId2Points = get_mesh_ray_intersection_points(scene, startPoint, direction);

					for (auto it = geoId2Points.begin(); it != geoId2Points.end(); it++)
					{
						int geoId = it->first;
						int objId = geoId2objId.at(geoId);

						const std::vector<gp_Pnt>& pnts = it->second;
						objId2Points[objId] = pnts;
					}

					for (auto it = objId2Points.begin(); it != objId2Points.end(); it++)
					{
						int objId = it->first;

						const std::vector<gp_Pnt>& pnts = it->second;

						if (pnts.size() <= 1)
						{ //
							// for (int n = 0; n < ZOX_LINES.size(); n++) {
							for (int n = 0; n < ZOX_COUNT; n++)
							{

								size_t next_face_index = i + n * ONE_XOY_FACES_COUNT;
								// set face material air_box
								int air_box_obj_id = -1;
								set_face_belong_to_obj(XOY_FACES[next_face_index], air_box_obj_id);
							}

							continue;
						}

						// update face material
						for (int k = 0; k < pnts.size(); k += 2)
						{

							if (k >= pnts.size() || k + 1 >= pnts.size())
							{
								break;
							}
							const gp_Pnt& pnt0 = pnts[k];
							const gp_Pnt& pnt1 = pnts[k + 1];
							double minZ = std::min(pnt0.Z(), pnt1.Z());
							double maxZ = std::max(pnt0.Z(), pnt1.Z());

							// update cell material
							// for (int n = 0; n < ZOX_LINES.size(); n++) {
							for (int n = 0; n < ZOX_COUNT; n++)
							{

								size_t next_face_index = i + n * ONE_XOY_FACES_COUNT;

								const gp_Pnt& face_center = get_xoy_face_center(XOY_FACES[next_face_index], XOY_LINES, YOZ_LINES, ZOX_LINES);

								double z = face_center.Z();

								double tol = 1e-7;
								if (n < ZOX_COUNT - 1)
								{
									tol = (ZOX_LINES[n + 1] - ZOX_LINES[n]) / 2;
								}
								else
								{
									tol = (ZOX_LINES[n] - ZOX_LINES[n - 1]) / 2;
								}

								if (z >= (minZ + tol) && z <= (maxZ - tol))
								{ // inner obj
									set_face_belong_to_obj(XOY_FACES[next_face_index], objId);
								}
								else if (Abs(z - minZ) < tol || Abs(z - maxZ) < tol)
								{ // middle
									int air_box_obj_id = -1;
									set_face_belong_to_obj(XOY_FACES[next_face_index], air_box_obj_id);
									set_face_belong_to_obj(XOY_FACES[next_face_index], objId);
								}
								else
								{ // outer
									// set face material air_box
									int air_box_obj_id = -1;
									set_face_belong_to_obj(XOY_FACES[next_face_index], air_box_obj_id);
								}
							}
						}
					}

					if (objId2Points.empty())
					{

						for (int n = 0; n < ZOX_COUNT; n++)
						{

							size_t next_face_index = i + n * ONE_XOY_FACES_COUNT;

							// set face material air_box
							int air_box_obj_id = -1;
							set_face_belong_to_obj(XOY_FACES[next_face_index], air_box_obj_id);
						}
					}
				}
			}

			// 2. hanle YOZ ray
#pragma omp section
			{
				size_t XOY_COUNT = (is_unreal == false) ? XOY_LINES.size() : (XOY_LINES.size() - 2);
				size_t ONE_YOZ_FACES_COUNT = (YOZ_LINES.size() - 1) * (ZOX_LINES.size() - 1);

				for (int i = 0; i < ONE_YOZ_FACES_COUNT; i++)
				{

					const FACE& face = YOZ_FACES[i];
					const gp_Pnt& face_center = get_yoz_face_center(face, XOY_LINES, YOZ_LINES, ZOX_LINES);

					gp_Pnt startPoint(face_center.X() - 1, face_center.Y(), face_center.Z());

					gp_Dir direction(1, 0, 0);
					///////////////////////////////////////////// ray occ obj

					std::unordered_map<int, std::vector<gp_Pnt>> objId2Points = get_occ_shape_ray_intersection_points(occ_shape_list, startPoint, direction);
					const std::unordered_map<int, std::vector<gp_Pnt>>& geoId2Points = get_mesh_ray_intersection_points(scene, startPoint, direction);

					for (auto it = geoId2Points.begin(); it != geoId2Points.end(); it++)
					{
						int geoId = it->first;
						int objId = geoId2objId.at(geoId);

						const std::vector<gp_Pnt>& pnts = it->second;
						objId2Points[objId] = pnts;
					}

					for (auto it = objId2Points.begin(); it != objId2Points.end(); it++)
					{
						int objId = it->first;

						const std::vector<gp_Pnt>& pnts = it->second;

						if (pnts.size() <= 1)
						{ //
							for (int n = 0; n < XOY_COUNT; n++)
							{

								size_t next_face_index = i + n * ONE_YOZ_FACES_COUNT;

								// set face material air_box
								int air_box_obj_id = -1;
								set_face_belong_to_obj(YOZ_FACES[next_face_index], air_box_obj_id);
							}

							continue;
						}

						// update face material
						for (int k = 0; k < pnts.size(); k += 2)
						{

							if (k >= pnts.size() || k + 1 >= pnts.size())
							{
								break;
							}
							const gp_Pnt& pnt0 = pnts[k];
							const gp_Pnt& pnt1 = pnts[k + 1];
							double minX = std::min(pnt0.X(), pnt1.X());
							double maxX = std::max(pnt0.X(), pnt1.X());

							// update cell material
							for (int n = 0; n < XOY_COUNT; n++)
							{

								size_t next_face_index = i + n * ONE_YOZ_FACES_COUNT;
								FACE next_face = YOZ_FACES[next_face_index];

								const gp_Pnt& face_center = get_yoz_face_center(next_face, XOY_LINES, YOZ_LINES, ZOX_LINES);

								double x = face_center.X();

								double tol = 1e-7;
								if (n < XOY_COUNT - 1)
								{
									tol = (XOY_LINES[n + 1] - XOY_LINES[n]) / 2;
								}
								else
								{
									tol = (XOY_LINES[n] - XOY_LINES[n - 1]) / 2;
								}

								if (x >= (minX + tol) && x <= (maxX - tol))
								{ // inner obj
									set_face_belong_to_obj(YOZ_FACES[next_face_index], objId);
								}
								else if (Abs(x - minX) < tol || Abs(x - maxX) < tol)
								{ // middle
									int air_box_obj_id = -1;
									set_face_belong_to_obj(YOZ_FACES[next_face_index], air_box_obj_id);
									set_face_belong_to_obj(YOZ_FACES[next_face_index], objId);
								}
								else
								{
									// set face material air_box
									int air_box_obj_id = -1;
									set_face_belong_to_obj(YOZ_FACES[next_face_index], air_box_obj_id);
								}
							}
						}
					}

					if (objId2Points.empty())
					{
						for (int n = 0; n < XOY_COUNT; n++)
						{

							size_t next_face_index = i + n * ONE_YOZ_FACES_COUNT;

							// set face material air_box
							int air_box_obj_id = -1;
							set_face_belong_to_obj(YOZ_FACES[next_face_index], air_box_obj_id);
						}
					}
				}
			}

			// 3. handle ZOX
#pragma omp section
			{
				size_t YOZ_COUNT = (is_unreal == false) ? YOZ_LINES.size() : (YOZ_LINES.size() - 2);
				size_t ONE_ZOX_FACES_COUNT = (ZOX_LINES.size() - 1) * (XOY_LINES.size() - 1);

				for (int i = 0; i < ONE_ZOX_FACES_COUNT; i++)
				{

					const FACE& face = ZOX_FACES[i];
					const gp_Pnt& face_center = get_zox_face_center(face, XOY_LINES, YOZ_LINES, ZOX_LINES);

					gp_Pnt startPoint(face_center.X(), face_center.Y(), face_center.Z() - 1);

					gp_Dir direction(0, 1, 0);
					///////////////////////////////////////////// ray occ obj

					std::unordered_map<int, std::vector<gp_Pnt>> objId2Points = get_occ_shape_ray_intersection_points(occ_shape_list, startPoint, direction);
					//	///////////////////////////////////////////// ray mesh obj
					const std::unordered_map<int, std::vector<gp_Pnt>>& geoId2Points = get_mesh_ray_intersection_points(scene, startPoint, direction);

					for (auto it = geoId2Points.begin(); it != geoId2Points.end(); it++)
					{
						int geoId = it->first;
						int objId = geoId2objId.at(geoId);

						const std::vector<gp_Pnt>& pnts = it->second;
						objId2Points[objId] = pnts;
					}

					for (auto it = objId2Points.begin(); it != objId2Points.end(); it++)
					{
						int objId = it->first;

						const std::vector<gp_Pnt>& pnts = it->second;

						if (pnts.size() <= 1)
						{ //
							for (int n = 0; n < YOZ_COUNT; n++)
							{

								size_t next_face_index = i + n * ONE_ZOX_FACES_COUNT;

								// set face material air_box
								int air_box_obj_id = -1;
								set_face_belong_to_obj(ZOX_FACES[next_face_index], air_box_obj_id);
							}

							continue;
						}

						// update face material
						for (int k = 0; k < pnts.size(); k += 2)
						{

							if (k >= pnts.size() || k + 1 >= pnts.size())
							{
								break;
							}
							const gp_Pnt& pnt0 = pnts[k];
							const gp_Pnt& pnt1 = pnts[k + 1];
							double minY = std::min(pnt0.Y(), pnt1.Y());
							double maxY = std::max(pnt0.Y(), pnt1.Y());

							// update cell material
							for (int n = 0; n < YOZ_COUNT; n++)
							{

								size_t next_face_index = i + n * ONE_ZOX_FACES_COUNT;

								const gp_Pnt& face_center = get_zox_face_center(ZOX_FACES[next_face_index], XOY_LINES, YOZ_LINES, ZOX_LINES);

								double y = face_center.Y();

								double tol = 1e-7;
								if (n < YOZ_COUNT - 1)
								{
									tol = (YOZ_LINES[n + 1] - YOZ_LINES[n]) / 2;
								}
								else
								{
									tol = (YOZ_LINES[n] - YOZ_LINES[n - 1]) / 2;
								}

								if (y >= (minY + tol) && y <= (maxY - tol))
								{ // inner obj
									set_face_belong_to_obj(ZOX_FACES[next_face_index], objId);
								}
								else if (Abs(y - minY) < tol || Abs(y - maxY) < tol)
								{ // middle
									int air_box_obj_id = -1;
									set_face_belong_to_obj(ZOX_FACES[next_face_index], air_box_obj_id);
									set_face_belong_to_obj(ZOX_FACES[next_face_index], objId);
								}
								else
								{
									// set face material air_box
									int air_box_obj_id = -1;
									set_face_belong_to_obj(ZOX_FACES[next_face_index], air_box_obj_id);
								}
							}
						}
					}

					if (objId2Points.empty())
					{
						for (int n = 0; n < YOZ_COUNT; n++)
						{

							size_t next_face_index = i + n * ONE_ZOX_FACES_COUNT;

							// set face material air_box
							int air_box_obj_id = -1;
							set_face_belong_to_obj(ZOX_FACES[next_face_index], air_box_obj_id);
						}
					}
				}
			}
		}

		return true;
	}

	// zxm

	std::vector<int> gen_PEC_obj_id(std::unordered_map<int, std::string> objId2matId, std::unordered_map<std::string, configor::json>& _Material)
	{
		std::vector<int> PEC_Obj_id;

		for (auto iter = _Material.begin(); iter != _Material.end(); ++iter)
		{
			if (iter->second["type"] == 1)
			{
				for (auto iter2 = objId2matId.begin(); iter2 != objId2matId.end(); ++iter2)
				{
					if (iter2->second == iter->first)
					{
						PEC_Obj_id.emplace_back(iter2->first);
					}
				}
			}
		}
		return PEC_Obj_id;
	}

	void gen_mesh_info(Mesh& _mesh, int XOY_LINES_COUNT, int YOZ_LINES_COUNT, int ZOX_LINES_COUNT, std::vector<FACE> xoy_face, std::vector<FACE> yoz_face, std::vector<FACE> zox_face, std::array<int, 2> boundary_level_x, std::array<int, 2> boundary_level_y, std::array<int, 2> boundary_level_z, double min_cell_size)
	{
		double PI = 3.141592653589793;
		double Mu = 4 * PI * 1e-7;
		double DIE = 8.854187817 * 1e-12;
		double cs = 1.0 / std::sqrt(DIE * Mu);

		//_mesh.layer_for_PML = 6;

		_mesh.boundary_layer[0] = boundary_level_x[0];
		_mesh.boundary_layer[1] = boundary_level_x[1];
		_mesh.boundary_layer[2] = boundary_level_y[0];
		_mesh.boundary_layer[3] = boundary_level_y[1];
		_mesh.boundary_layer[4] = boundary_level_z[0];
		_mesh.boundary_layer[5] = boundary_level_z[1];

		_mesh.NX = XOY_LINES_COUNT - 1;
		_mesh.NY = YOZ_LINES_COUNT - 1;
		_mesh.NZ = ZOX_LINES_COUNT - 1;

		_mesh.dx.resize(_mesh.NX);
		_mesh.dy.resize(_mesh.NY);
		_mesh.dz.resize(_mesh.NZ);

		for (int i = 0; i < _mesh.NX; i++)
		{
			_mesh.dx[i] = xoy_face[i].h * 1e-3;
		}
		for (int j = 0; j < _mesh.NY; j++)
		{
			_mesh.dy[j] = yoz_face[j].h * 1e-3;
		}
		for (int k = 0; k < _mesh.NZ; k++)
		{
			_mesh.dz[k] = zox_face[k].h * 1e-3;
		}

		_mesh.dt = min_cell_size / (cs * std::sqrt(3));

		H_inf _h_inf;
		_h_inf.CP = 0.0;
		_h_inf.CQ_S1 = 0.0;
		_h_inf.l[0] = 0.0;
		_h_inf.l[1] = 0.0;
		_h_inf.l[2] = 0.0;
		_h_inf.l[3] = 0.0;
		_mesh.h_inf.emplace_back(_h_inf);
	}

	std::unordered_map<int, int> gen_objId2meshmatlId(std::unordered_map<std::string, configor::json> materials, std::unordered_map<int, std::string> objId2matId, Mesh& _mesh)
	{
		double PI = 3.141592653589793;
		double Mu = 4 * PI * 1e-7;
		double DIE = 8.854187817 * 1e-12;
		double z0 = sqrt(Mu / DIE);
		// new
		std::unordered_map<std::string, int> matId2meshmatId;
		std::unordered_map<int, int> objId2meshmatlId;
		for (auto iter = materials.begin(); iter != materials.end(); ++iter)
		{
			if (iter->second["dielectric"]["type"] == 0)
			{
				_mesh.epsilon.emplace_back(iter->second["dielectric"]["epsilon"]);
				_mesh.sigma.emplace_back(iter->second["dielectric"]["sigma"]);
			}
			if (iter->second["permeability"]["type"] == 0)
			{
				_mesh.mu.emplace_back(iter->second["permeability"]["mu"]);
				_mesh.sigma_mu.emplace_back(iter->second["permeability"]["sigma_m"]);
			}
			matId2meshmatId[iter->first] = _mesh.epsilon.size() - 1;

			double eps = _mesh.epsilon.back() * DIE;
			double miu = _mesh.mu.back() * Mu;

			double temp1 = _mesh.sigma.back() * _mesh.dt / eps / 2.0;
			_mesh.CA.emplace_back((-2 * temp1) / (1 + temp1));
			_mesh.CB.emplace_back((_mesh.dt / eps / z0) / (1 + temp1));
			double temp2 = _mesh.sigma_mu.back() * _mesh.dt / miu / 2.0;
			_mesh.CP.emplace_back((-2 * temp2) / (1 + temp2));
			_mesh.CQ.emplace_back((_mesh.dt * z0 / miu) / (1 + temp2));
		}
		for (auto iter2 = objId2matId.begin(); iter2 != objId2matId.end(); ++iter2)
		{
			for (auto iter = matId2meshmatId.begin(); iter != matId2meshmatId.end(); ++iter)
			{
				if (iter->first == iter2->second)
				{
					objId2meshmatlId[iter2->first] = iter->second;
				}
			}
		}

		return objId2meshmatlId;
	}

	void gen_rayhit(RTCRayHit& rayhit, float ox, float oy, float oz, gp_Dir direction, float t_near)
	{
		rayhit.ray.org_x = ox;
		rayhit.ray.org_y = oy;
		rayhit.ray.org_z = oz;
		rayhit.ray.dir_x = direction.X();
		rayhit.ray.dir_y = direction.Y();
		rayhit.ray.dir_z = direction.Z();
		rayhit.ray.tnear = t_near;
		rayhit.ray.tfar = std::numeric_limits<float>::infinity();
		rayhit.ray.mask = -1;
		rayhit.ray.flags = 0;

		rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
		rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
		rayhit.hit.instPrimID[0] = RTC_INVALID_GEOMETRY_ID;
	}

	std::vector<gp_Pnt> get_occ_shape_ray_intersection_points_new(const TopoDS_Shape& shape, const gp_Pnt& originPnt, const gp_Dir& dir)
	{
		double tol = 4e-6;
		std::vector<gp_Pnt> pnts;

		// create ray
		Handle(Geom_Curve) ray = new Geom_Line(originPnt, dir);

		GeomAPI_IntCS intersector;

		TopExp_Explorer exp(shape, TopAbs_FACE); // iterator

		gp_XYZ c = dir.XYZ();
		gp_XYZ startp = originPnt.XYZ();

		while (exp.More())
		{
			const TopoDS_Face& face = TopoDS::Face(exp.Current());
			const Handle(Geom_Surface)& surface = BRep_Tool::Surface(face);
			intersector.Perform(ray, surface);

			if (intersector.IsDone() && intersector.NbPoints() > 0)
			{
				for (int i = 1; i <= intersector.NbPoints(); i++)
				{
					const gp_Pnt& pnt = intersector.Point(i);	// 获取交点的坐标
					if ((pnt.XYZ().Dot(c) - startp.Dot(c)) > 0) // 与dir同向
					{
						pnts.push_back(pnt);
					}
				}
			}

			exp.Next();
		}

		std::vector<gp_Pnt> resultPnts;

		if (pnts.size() > 1)
		{

			// sort
			std::sort(pnts.begin(), pnts.end(), [c](const gp_Pnt& p1, const gp_Pnt& p2)
			{
				gp_XYZ c1 = p1.XYZ();
				gp_XYZ c2 = p2.XYZ();
				return (c1.Dot(c) < c2.Dot(c)); });

			// del the repeat value
			for (int i = 0; i < pnts.size() - 1; i++)
			{

				if (i == 0)
				{
					resultPnts.push_back(pnts[0]);
				}
				gp_XYZ c1 = resultPnts[resultPnts.size() - 1].XYZ();
				gp_XYZ c2 = pnts[i + 1].XYZ();
				if (Abs(c1.Dot(c) - c2.Dot(c)) >= tol)
				{

					// pnts[i + 1].DumpJson(std::cout);
					resultPnts.push_back(pnts[i + 1]);
				}
			}
		}
		return resultPnts;
	}

	std::vector<std::pair<int, gp_Pnt>> get_occ_shape_ray_intersection_points_new(TopTools_DataMapOfIntegerShape& occ_shape_list, const gp_Pnt& originPnt, const gp_Dir& dir)
	{
		double tol = 4 * 1e-6;

		std::vector<std::pair<int, gp_Pnt>> intersectPoints;

		gp_Lin ray(originPnt, dir);

		Handle(Geom_Line) aGeomLine = new Geom_Line(ray);
		aGeomLine->SetLocation(originPnt);
		aGeomLine->SetDirection(dir);

		TopTools_DataMapIteratorOfDataMapOfIntegerShape iter(occ_shape_list);
		// 遍历数据映射中的所有元素
		for (; iter.More(); iter.Next())
		{

			int objId = iter.Key();
			const TopoDS_Shape& shape = iter.Value();
			// ray geo obj
			std::vector<gp_Pnt> resultPnts = get_occ_shape_ray_intersection_points_new(shape, originPnt, dir);
			for (int k = 0; k < resultPnts.size(); k++)
			{
				intersectPoints.emplace_back(std::make_pair(objId, resultPnts[k]));
			}
		}
		return intersectPoints;
	}

	std::vector<std::pair<int, gp_Pnt>> sort_ray_intersection_points(std::vector<std::pair<int, gp_Pnt>> occ_objId2Points, std::vector<std::pair<int, gp_Pnt>> mesh_objId2Points, SEGMENT& _segment, int& total_cross, int& seg_cross)
	{
		std::vector<std::pair<int, gp_Pnt>> total_obj2point = occ_objId2Points;
		total_obj2point.insert(total_obj2point.end(), mesh_objId2Points.begin(), mesh_objId2Points.end());

		gp_XYZ c = _segment.direction.XYZ();

		// sort
		if (total_obj2point.size() > 1)
		{
			// sort
			std::sort(total_obj2point.begin(), total_obj2point.end(), [c](const std::pair<int, gp_Pnt>& p1, const std::pair<int, gp_Pnt>& p2)
			{
				gp_XYZ c1 = p1.second.XYZ();
				gp_XYZ c2 = p2.second.XYZ();
				return (c1.Dot(c) < c2.Dot(c)); });
		}

		// generate seg_cross
		double total_distance = 0.0;
		total_cross = total_obj2point.size();

		for (int i = 0; i < total_cross; i++)
		{
			double distance = 0;
			if (i == 0)
			{
				distance = _segment.startPoint.Distance(total_obj2point[i].second);
			}
			else
			{
				distance = total_obj2point[i - 1].second.Distance(total_obj2point[i].second);
			}
			total_distance += distance;
			// std::cout << "total_distance: " << total_distance << std::endl;
			if (total_distance <= _segment.length)
			{
				_segment.cross_dir.emplace_back(distance);
				seg_cross++;
			}
		}
		// std::cout << "total cross: " << total_cross <<" seg_cross: "<< seg_cross << std::endl;
		return total_obj2point;
	}

	void gen_conformal_status(bool is_unreal, std::vector<int> obj_id1, std::vector<int> obj_id2, std::vector<int> PEC_id, std::vector<int>& obj_id, int& conformal_status, int& case22_status)
	{
		int size1 = obj_id1.size();
		int size2 = obj_id2.size();
		if (!(((size1 == 1) && (size2 == 2)) || ((size1 == 2) && (size2 == 1)) || ((size1 == 2) && (size2 == 2))))
		{
			//Logger::info("obj num is not match!");
			conformal_status = -1;
			return;
		}
		else if (((size1 == 1) && (size2 == 2)) || ((size1 == 2) && (size2 == 1)))
		{
			if (obj_id2[0] != obj_id1[0])
			{
				std::cout << "case 1&3 obj id is not match!" << std::endl;
				//Logger::info("case 1&3 obj id is not match!");
				return;
			}
			else
			{
				if ((is_unreal == false) && (std::find(PEC_id.begin(), PEC_id.end(), obj_id1[0]) != PEC_id.end()))
				{
					if (size1 == 1)
					{
						obj_id = obj_id2;
					}
					else if (size2 == 1)
					{
						obj_id = obj_id1;
					}
					conformal_status = 1;
				}
				else
				{
					if (size1 == 1)
					{
						obj_id = { -1, obj_id1[0] };
					}
					else if (size2 == 1)
					{
						obj_id = obj_id1;
					}
					conformal_status = 2;
				}
			}
		}
		else if ((size1 == 2) && (size2 == 2))
		{
			if ((obj_id1[0] != obj_id2[1]) || (obj_id1[1] != obj_id2[0]))
			{
				std::cout << "case 2&2 obj id is not match!" << std::endl;
				//Logger::info("case 2&2 obj id is not match!");
				return;
			}
			else
			{
				bool obj1 = false, obj2 = false;
				if (is_unreal == false)
				{
					obj1 = std::find(PEC_id.begin(), PEC_id.end(), obj_id1[0]) != PEC_id.end();
					obj2 = std::find(PEC_id.begin(), PEC_id.end(), obj_id2[0]) != PEC_id.end();
				}
				if ((obj1 == true) && (obj2 == true))
				{
					conformal_status = 1;
					case22_status = 0; // 0:multi pec conformal
					obj_id = obj_id1;
				}
				else if (obj1 == true)
				{
					conformal_status = 1;
					case22_status = 1;
					obj_id = obj_id1;
				}
				else if (obj2 == true)
				{
					conformal_status = 1;
					case22_status = 2;
					obj_id = obj_id2;
				}
				else
				{
					conformal_status = 2;
					case22_status = -1;
					obj_id = obj_id1; // medium - Conformal: default obj_id = obj_id1
				}
			}
		}
	}

	void gen_conformal_status_2points(bool is_unreal, std::vector<int> obj_id1, std::vector<int> obj_id2, std::vector<int> PEC_id, std::vector<int>& obj_id, int& conformal_status, int& case22_status)
	{
		int size1 = obj_id1.size();
		int size2 = obj_id2.size();
		if (!(((size1 == 1) && (size2 == 2)) || ((size1 == 2) && (size2 == 1)) || ((size1 == 2) && (size2 == 2))))
		{
			//Logger::info("obj num is not match!");
			conformal_status = -1;
			return;
		}
		else if (((size1 == 1) && (size2 == 2)) || ((size1 == 2) && (size2 == 1))) //[air -> obj] & [obj -> air]
		{
			if (obj_id2[0] != obj_id1[0])
			{
				std::cout << "case 1&3 obj id is not match!" << std::endl;
				//Logger::info("case 1&3 obj id is not match!");
				abort;
			}
			else
			{
				if ((is_unreal == false) && (std::find(PEC_id.begin(), PEC_id.end(), obj_id1[0]) != PEC_id.end()))
				{
					if (size1 == 1)
					{
						obj_id = obj_id2;
					}
					else if (size2 == 1)
					{
						obj_id = obj_id1;
					}
					conformal_status = 1;
				}
				else
				{
					if (size1 == 1)
					{
						obj_id = { -1, obj_id1[0] };
					}
					else if (size2 == 1)
					{
						obj_id = obj_id1;
					}
					conformal_status = 2;
				}
			}
		}
		else if ((size1 == 2) && (size2 == 2))
		{
			if ((obj_id1[0] != obj_id2[1]) || (obj_id1[1] != obj_id2[0]))
			{
				std::cout << "case 2&2 obj id is not match!" << std::endl;
				//Logger::info("case 2&2 obj id is not match!");
				abort;
			}
			else
			{
				bool obj1 = false, obj2 = false;
				if (is_unreal == false)
				{
					obj1 = std::find(PEC_id.begin(), PEC_id.end(), obj_id1[0]) != PEC_id.end();
					obj2 = std::find(PEC_id.begin(), PEC_id.end(), obj_id2[0]) != PEC_id.end();
				}
				if ((obj1 == true) && (obj2 == true))
				{
					conformal_status = 1;
					case22_status = 0; // 0:multi pec conformal
					obj_id = obj_id1;
				}
				else if (obj1 == true)
				{
					conformal_status = 1;
					case22_status = 1;
					obj_id = obj_id1;
				}
				else if (obj2 == true)
				{
					conformal_status = 1;
					case22_status = 2;
					obj_id = obj_id2;
				}
				else
				{
					conformal_status = 2;
					case22_status = -1;
					obj_id = obj_id1; // medium - Conformal: default obj_id = obj_id1
				}
			}
		}
	}

	void check_conformal_restraint(FACE& _face, std::vector<int> Obj_id, std::vector<double> dir, double area, double face_area)
	{
		std::vector<double> seg(4);
		seg[0] = dir[0] / _face.h / 1e-3;
		seg[1] = dir[1] / _face.v / 1e-3;
		seg[2] = dir[2] / _face.h / 1e-3;
		seg[3] = dir[3] / _face.v / 1e-3;
		double max_seg = *std::max_element(seg.begin(), seg.end());
		double area_rate = area / face_area;
		if ((area_rate <= 0.05) || ((max_seg / area_rate) >= 12))
		{
			_face.Conformal_status = -1;
			int obj = Obj_id[1];
			Obj_id = std::vector<int>{ obj };
		}
	}

	void gen_face_conformal_status_2(FACE& _face, bool is_unreal, std::array<SEGMENT, 4>& segment, std::vector<int> PEC_id)
	{
		double face_v = _face.v * 1e-3;
		double face_h = _face.h * 1e-3;
		double face_area = _face.v * _face.h;

		std::array<double, 4> seg_dir = { segment[0].cross_dir[0] * 1e-3, segment[1].cross_dir[0] * 1e-3, segment[2].cross_dir[0] * 1e-3, segment[3].cross_dir[0] * 1e-3 };
		std::array<int, 4> seg_status; // crosspoint status: -1 no cross, 0 regular crosss, 1 corner point
		std::array<int, 4> dir_status; // direction status: -1 no cross, 0 case3, 1 case1, 2 case2
		for (int seg = 0; seg < 4; seg++)
		{
			dir_status[seg] = segment[seg].status[0];
			if (seg_dir[seg] == 0)
			{
				seg_status[seg] = -1;
			}
			else if (seg_dir[seg] == segment[seg].length)
			{
				seg_status[seg] = 1;
			}
			else
			{
				seg_status[seg] = 0;
			}
		}
		_face.Obj_id.resize(2);

		int case22_status = -1;
		// dir:-x +x -y +y / -y +y -z +z / -z +z -x +x
		if ((seg_status[3] != -1) && (seg_status[1] != -1)) // diagnal - four cases
		{
			int case_status = -1;
			double area = 0.0;
			// PEC - Conformal or Medium - Conformal
			gen_conformal_status(is_unreal, segment[3].Obj_ID, segment[1].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
			if ((seg_status[3] == 1) && (seg_status[1] == 1))
			{
				area = face_area / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
					{
						_face.dir = { 0.0, face_v, face_h, 0.0 };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
					else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
					{
						_face.dir = { face_v, 0.0, 0.0, face_h };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { area, area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
			else if ((seg_status[3] == 1) && (seg_status[1] == 0))
			{
				area = face_v * seg_dir[1] / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
					{
						_face.dir = { 0.0, face_v, face_h, face_h - seg_dir[1] };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
					else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
					{
						_face.dir = { face_v, 0.0, 0.0, seg_dir[1] };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { face_area - area, area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
			else if ((seg_status[1] == 1) && (seg_status[3] == 0))
			{
				double area = face_v * (face_h - seg_dir[3]) / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
					{
						_face.dir = { 0.0, face_v, seg_dir[3], 0.0 };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
					else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
					{
						_face.dir = { face_v, 0.0, face_h - seg_dir[3], face_h };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { area, face_area - area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
			else if ((seg_status[3] == 0) && (seg_status[1] == 0))
			{
				double area = face_v * (face_h - seg_dir[3] + seg_dir[1]) / 2.0; // Down seg area
				if (_face.Conformal_status == 1)								 // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - Down Corner
					{
						_face.dir = { 0.0, face_v, seg_dir[3], face_h - seg_dir[1] };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
					else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - Up Corner
					{
						_face.dir = { face_v, 0.0, face_h - seg_dir[3], seg_dir[1] };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { face_area - area, area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
		}
		else if ((seg_status[2] != -1) && (seg_status[0] != -1)) // diagnal - four cases
		{
			gen_conformal_status(is_unreal, segment[2].Obj_ID, segment[0].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
			if ((seg_status[2] == 1) && (seg_status[0] == 1))
			{
				double area = face_v * face_h / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
					{
						_face.dir = { 0.0, face_v, 0.0, face_h };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
					else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
					{
						_face.dir = { face_v, 0.0, face_h, 0.0 };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { face_area - area, area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
			else if ((seg_status[2] == 1) && (seg_status[0] == 0))
			{
				double area = seg_dir[0] * face_h / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
					{
						_face.dir = { face_v - seg_dir[0], face_v, 0.0, face_h };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
					else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
					{
						_face.dir = { seg_dir[0], 0.0, face_h, 0.0 };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { face_area - area, area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
			else if ((seg_status[0] == 1) && (seg_status[2] == 0))
			{
				double area = seg_dir[2] * face_h / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
					{
						_face.dir = { 0.0, seg_dir[2], 0.0, face_h };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
					else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
					{
						_face.dir = { face_v, face_v - seg_dir[2], face_h, 0.0 };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { area, face_area - area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
			else if ((seg_status[2] == 0) && (seg_status[0] == 0))
			{
				double area = face_h * (face_v - seg_dir[2] + seg_dir[0]) / 2.0; // 左边面积
				if (_face.Conformal_status == 1)								 // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - Left Corner
					{
						_face.dir = { face_v - seg_dir[0], seg_dir[2], 0.0, face_h };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
					else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - Right Corner
					{
						_face.dir = { seg_dir[0], face_v - seg_dir[2], face_h, 0.0 };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { face_area - area, area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
		}
		else if ((seg_status[1] != -1) && (seg_status[0] != -1)) // Adjacent edges - two cases
		{
			gen_conformal_status(is_unreal, segment[1].Obj_ID, segment[0].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
			if ((seg_status[1] == 1) && (seg_status[0] == 0))
			{
				double area = (face_v - seg_dir[0]) * face_h / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[1] == 0) || (case22_status == 2)) // PEC - LeftUp Corner
					{
						_face.dir = { face_v - seg_dir[0], 0.0, 0.0, face_h };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
					else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightDown Corner
					{
						_face.dir = { seg_dir[0], face_v, face_h, 0.0 };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { area, face_area - area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
			else if ((seg_status[0] == 0) && (seg_status[1] == 0))
			{
				double area = (face_v - seg_dir[0]) * seg_dir[1] / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[1] == 0) || (case22_status == 2)) // PEC - LeftUp Corner
					{
						_face.dir = { face_v - seg_dir[0], 0.0, 0.0, seg_dir[1] };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
					else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightDown Corner
					{
						_face.dir = { seg_dir[0], face_v, face_h, face_h - seg_dir[1] };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { area, face_area - area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
		}
		else if ((seg_status[2] != -1) && (seg_status[1] != -1)) // Adjacent edges - two cases
		{
			gen_conformal_status(is_unreal, segment[2].Obj_ID, segment[1].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
			if ((seg_status[2] == 1) && (seg_status[1] == 0))
			{
				double area = (face_h - seg_dir[1]) * face_v / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
					{
						_face.dir = { 0.0, face_v, 0.0, face_h - seg_dir[1] };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
					else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - RightUp Corner
					{
						_face.dir = { face_v, 0.0, face_h, seg_dir[1] };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { area, face_area - area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
			else if ((seg_status[1] == 0) && (seg_status[2] == 0))
			{
				double area = seg_dir[2] * (face_h - seg_dir[1]) / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
					{
						_face.dir = { 0.0, seg_dir[2], 0.0, face_h - seg_dir[1] };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
					else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - RightUp Corner
					{
						_face.dir = { face_v, face_v - seg_dir[2], face_h, seg_dir[1] };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { area, face_area - area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
		}
		else if ((seg_status[3] != -1) && (seg_status[2] != -1)) // Adjacent edges - two cases
		{
			gen_conformal_status(is_unreal, segment[3].Obj_ID, segment[2].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
			if ((seg_status[3] == 1) && (seg_status[2] == 0))
			{
				double area = face_h * (face_v - seg_dir[2]) / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
					{
						_face.dir = { 0.0, face_v - seg_dir[2], face_h, 0.0 };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
					else if ((dir_status[2] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
					{
						_face.dir = { face_v, seg_dir[2], 0.0, face_h };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { area, face_area - area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
			else if ((seg_status[2] == 0) && (seg_status[3] == 0))
			{
				double area = seg_dir[3] * (face_v - seg_dir[2]) / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
					{
						_face.dir = { 0.0, face_v - seg_dir[2], seg_dir[3], 0.0 };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
					else if ((dir_status[2] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
					{
						_face.dir = { face_v, seg_dir[2], face_h - seg_dir[3], face_h };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { area, face_area - area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
		}
		else if ((seg_status[3] != -1) && (seg_status[0] != -1)) // Adjacent edges - two cases
		{
			gen_conformal_status(is_unreal, segment[3].Obj_ID, segment[0].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
			if ((seg_status[0] == 1) && (seg_status[3] == 0))
			{
				double area = (face_h - seg_dir[3]) * face_v / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
					{
						_face.dir = { 0.0, face_v, seg_dir[3], face_h };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
					else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
					{
						_face.dir = { face_v, 0.0, face_h - seg_dir[3], 0.0 };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { face_area - area, area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
			else if ((seg_status[3] == 0) && (seg_status[0] == 0))
			{
				double area = seg_dir[0] * (face_h - seg_dir[3]) / 2.0;
				if (_face.Conformal_status == 1) // PEC - Conformal
				{
					if (case22_status == 0)
					{
						_face.dir = { 0.0, 0.0, 0.0, 0.0 };
					}
					else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
					{
						_face.dir = { face_v - seg_dir[0], face_v, seg_dir[3], face_h };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
					}
					else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
					{
						_face.dir = { seg_dir[0], 0.0, face_h - seg_dir[3], 0.0 };
						check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
					}
				}
				else if (_face.Conformal_status == 2) // Medium - Conformal
				{
					_face.sec_area = { face_area - area, area };
				}
				else if (_face.Conformal_status == -1) // obj num is not match to be determined
				{
					_face.Obj_id = std::vector<int>{ -1 };
				}
			}
		}
	}

	void refine_grid(bool is_unreal, TopTools_DataMapOfIntegerShape& occ_shape_list, RTCScene& scene, const std::unordered_map<int, int>& geoId2objId, const std::unordered_map<int, int>& objId2meshmatlId, std::vector<int> PEC_id, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES, std::vector<FACE>& face, int plane)
	{
		double threshold = 4e-6;

		// unit: ray tracing [mm] filling refine grid [m]
		for (int face_id = 0; face_id < face.size(); face_id++)
		{
			FACE& _face = face[face_id];

			double face_v = _face.v * 1e-3;
			double face_h = _face.h * 1e-3;

			if (_face.flag == 0)
			{
				//std::cout << "face_id: " << face_id << " start" << std::endl;
				double face_area = face_v * face_h;

				size_t line0 = _face.line_indices[0]; // x line
				size_t line1 = _face.line_indices[1]; // y line
				size_t line2 = _face.line_indices[2]; // z line

				std::array<SEGMENT, 4> segment;
				segment[0].length = _face.h;
				segment[1].length = _face.v;
				segment[2].length = _face.h;
				segment[3].length = _face.v;

				if (plane == 1)
				{
					segment[0].startPoint = { XOY_LINES[line0], YOZ_LINES[line1], ZOX_LINES[line2] };
					segment[1].startPoint = { XOY_LINES[line0 + 1], YOZ_LINES[line1], ZOX_LINES[line2] };
					segment[2].startPoint = { XOY_LINES[line0 + 1], YOZ_LINES[line1 + 1], ZOX_LINES[line2] };
					segment[3].startPoint = { XOY_LINES[line0], YOZ_LINES[line1 + 1], ZOX_LINES[line2] };

					segment[0].direction = { 1, 0, 0 };
					segment[1].direction = { 0, 1, 0 };
					segment[2].direction = { -1, 0, 0 };
					segment[3].direction = { 0, -1, 0 };
				}
				else if (plane == 2)
				{
					segment[0].startPoint = { XOY_LINES[line2], YOZ_LINES[line0], ZOX_LINES[line1] };
					segment[1].startPoint = { XOY_LINES[line2], YOZ_LINES[line0 + 1], ZOX_LINES[line1] };
					segment[2].startPoint = { XOY_LINES[line2], YOZ_LINES[line0 + 1], ZOX_LINES[line1 + 1] };
					segment[3].startPoint = { XOY_LINES[line2], YOZ_LINES[line0], ZOX_LINES[line1 + 1] };

					segment[0].direction = { 0, 1, 0 };
					segment[1].direction = { 0, 0, 1 };
					segment[2].direction = { 0, -1, 0 };
					segment[3].direction = { 0, 0, -1 };
				}
				else if (plane == 3)
				{
					segment[0].startPoint = { XOY_LINES[line1], YOZ_LINES[line2], ZOX_LINES[line0] };
					segment[1].startPoint = { XOY_LINES[line1], YOZ_LINES[line2], ZOX_LINES[line0 + 1] };
					segment[2].startPoint = { XOY_LINES[line1 + 1], YOZ_LINES[line2], ZOX_LINES[line0 + 1] };
					segment[3].startPoint = { XOY_LINES[line1 + 1], YOZ_LINES[line2], ZOX_LINES[line0] };

					segment[0].direction = { 0, 0, 1 };
					segment[1].direction = { 1, 0, 0 };
					segment[2].direction = { 0, 0, -1 };
					segment[3].direction = { -1, 0, 0 };
				}

				// 非边界情况（射线入射与某条几何边完全重合导致交点数异常）
				for (int i = 0; i < 4; i++)
				{
					// occ ray tracing

					std::vector<std::pair<int, gp_Pnt>> occ_objId2Points = get_occ_shape_ray_intersection_points_new(occ_shape_list, segment[i].startPoint, segment[i].direction);

					// mesh ray tracing
					std::vector<std::pair<int, gp_Pnt>> mesh_obj2point;
					RTCRayHit rayhit;
					bool intersect = true;
					double distance = 0.0;
					float ox = segment[i].startPoint.X();
					float oy = segment[i].startPoint.Y();
					float oz = segment[i].startPoint.Z();
					while (intersect == true)
					{
						gen_rayhit(rayhit, ox, oy, oz, segment[i].direction);
						rtcIntersect1(scene, &rayhit);
						if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
						{
							intersect = true;
							distance = rayhit.ray.tfar + threshold;
							ox = rayhit.ray.org_x + distance * rayhit.ray.dir_x;
							oy = rayhit.ray.org_y + distance * rayhit.ray.dir_y;
							oz = rayhit.ray.org_z + distance * rayhit.ray.dir_z;

							int obj_id;
							if (rayhit.hit.geomID == -1)
							{
								obj_id = -1;
							}
							for (auto iter = geoId2objId.begin(); iter != geoId2objId.end(); ++iter)
							{
								if (iter->first == rayhit.hit.geomID)
								{
									obj_id = iter->second;
								}
							}
							mesh_obj2point.emplace_back(std::make_pair(obj_id, gp_Pnt{ ox, oy, oz }));
						}
						else
						{
							int obj_id;
							if (rayhit.hit.geomID == -1)
							{
								obj_id = -1;
							}
							intersect = false;
						}
					}

					// gen total ray seg_points info
					int seg_cross = 0;	 // 当前线段上的交点数
					int total_cross = 0; // 射线与几何体的总交点数
					std::vector<std::pair<int, gp_Pnt>> total_obj2point;
					total_obj2point = sort_ray_intersection_points(occ_objId2Points, mesh_obj2point, segment[i], total_cross, seg_cross);

					// gen segment crosspoint status
					if (seg_cross == 0)
					{
						segment[i].status.emplace_back(-1);
						segment[i].cross_num = 0;
						segment[i].cross_dir = std::vector<double>{ 0.0 };
					}
					else if (seg_cross == 1)
					{
						segment[i].cross_num = 1;
						//_face.cross_vtx.emplace_back(total_obj2point[0].second);
						if ((total_cross == seg_cross)) // case1. single cross [inside obj -> airbox]
						{
							segment[i].Obj_ID.emplace_back(total_obj2point[0].first);
							segment[i].Obj_ID.emplace_back(-1);
							segment[i].status.emplace_back(1);
						}
						else if ((total_obj2point[0].first != total_obj2point[1].first)) // case2. multi cross [obj1 -> obj2]
						{
							segment[i].Obj_ID.emplace_back(total_obj2point[0].first);
							segment[i].Obj_ID.emplace_back(total_obj2point[1].first);
							segment[i].status.emplace_back(2);
						}
						else if ((total_obj2point[0].first == total_obj2point[1].first)) // case3. multi cross [airbox -> obj]
						{
							segment[i].Obj_ID.emplace_back(total_obj2point[0].first);
							segment[i].status.emplace_back(0);
						}
					}
					else if (seg_cross == 2)
					{
						segment[i].status.emplace_back(-2);
						segment[i].cross_num = 0;
						segment[i].cross_dir = std::vector<double>{ 0.0 };
						// count2++;
						//  To be determined
					}
					else // seg_cross > 2
					{
						segment[i].status.emplace_back(-3);
						segment[i].cross_num = 0;
						segment[i].cross_dir = std::vector<double>{ 0.0 };
						// countn++;
						//  To be determined
					}
					// std::cout << "seg: " << face_id << " i: " << i << std::endl;
				}
				// std::cout << "face_id: " << face_id << " mid" << std::endl;
				//  generate face crosspoint info
				int face_cross_num = segment[0].cross_num + segment[1].cross_num + segment[2].cross_num + segment[3].cross_num;
				// calculating material segment area
				if ((face_cross_num == 0) || (face_cross_num == 1))
				{
					_face.Conformal_status = -1;
					_face.Obj_id.emplace_back(-1);
				}
				else if (face_cross_num == 2) // two segment areas
				{
					// gen_face_conformal_status_2(_face, is_unreal, segment, PEC_id);
					std::array<double, 4> seg_dir = { segment[0].cross_dir[0] * 1e-3, segment[1].cross_dir[0] * 1e-3, segment[2].cross_dir[0] * 1e-3, segment[3].cross_dir[0] * 1e-3 };
					std::array<int, 4> seg_status; // crosspoint status: -1 no cross, 0 regular crosss, 1 corner point
					std::array<int, 4> dir_status; // direction status: -1 no cross, 0 case3, 1 case1, 2 case2
					for (int seg = 0; seg < 4; seg++)
					{
						dir_status[seg] = segment[seg].status[0];
						if (seg_dir[seg] == 0)
						{
							seg_status[seg] = -1;
						}
						else if (seg_dir[seg] == segment[seg].length)
						{
							seg_status[seg] = 1;
						}
						else
						{
							seg_status[seg] = 0;
						}
					}
					_face.Obj_id.resize(face_cross_num);

					int case22_status = -1;
					// dir:-x +x -y +y / -y +y -z +z / -z +z -x +x
					if ((seg_status[3] != -1) && (seg_status[1] != -1)) // diagnal - four cases
					{
						int case_status = -1;
						double area = 0.0;
						// PEC - Conformal or Medium - Conformal
						gen_conformal_status(is_unreal, segment[3].Obj_ID, segment[1].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
						if ((seg_status[3] == 1) && (seg_status[1] == 1))
						{
							area = face_area / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
								{
									_face.dir = { 0.0, face_v, face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v, 0.0, 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[3] == 1) && (seg_status[1] == 0))
						{
							area = face_v * seg_dir[1] / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
								{
									_face.dir = { 0.0, face_v, face_h, face_h - seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
								else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v, 0.0, 0.0, seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[1] == 1) && (seg_status[3] == 0))
						{
							double area = face_v * (face_h - seg_dir[3]) / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
								{
									_face.dir = { 0.0, face_v, seg_dir[3], 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v, 0.0, face_h - seg_dir[3], face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[3] == 0) && (seg_status[1] == 0))
						{
							double area = face_v * (face_h - seg_dir[3] + seg_dir[1]) / 2.0; // Down seg area
							if (_face.Conformal_status == 1)								 // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - Down Corner
								{
									_face.dir = { 0.0, face_v, seg_dir[3], face_h - seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
								else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - Up Corner
								{
									_face.dir = { face_v, 0.0, face_h - seg_dir[3], seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
					}
					else if ((seg_status[2] != -1) && (seg_status[0] != -1)) // diagnal - four cases
					{
						gen_conformal_status(is_unreal, segment[2].Obj_ID, segment[0].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
						if ((seg_status[2] == 1) && (seg_status[0] == 1))
						{
							double area = face_v * face_h / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { 0.0, face_v, 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { face_v, 0.0, face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[2] == 1) && (seg_status[0] == 0))
						{
							double area = seg_dir[0] * face_h / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { face_v - seg_dir[0], face_v, 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { seg_dir[0], 0.0, face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[0] == 1) && (seg_status[2] == 0))
						{
							double area = seg_dir[2] * face_h / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { 0.0, seg_dir[2], 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { face_v, face_v - seg_dir[2], face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[2] == 0) && (seg_status[0] == 0))
						{
							double area = face_h * (face_v - seg_dir[2] + seg_dir[0]) / 2.0; // 左边面积
							if (_face.Conformal_status == 1)								 // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - Left Corner
								{
									_face.dir = { face_v - seg_dir[0], seg_dir[2], 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - Right Corner
								{
									_face.dir = { seg_dir[0], face_v - seg_dir[2], face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
					}
					else if ((seg_status[1] != -1) && (seg_status[0] != -1)) // Adjacent edges - two cases
					{
						gen_conformal_status(is_unreal, segment[1].Obj_ID, segment[0].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
						if ((seg_status[1] == 1) && (seg_status[0] == 0))
						{
							double area = (face_v - seg_dir[0]) * face_h / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[1] == 0) || (case22_status == 2)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v - seg_dir[0], 0.0, 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightDown Corner
								{
									_face.dir = { seg_dir[0], face_v, face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[0] == 0) && (seg_status[1] == 0))
						{
							double area = (face_v - seg_dir[0]) * seg_dir[1] / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[1] == 0) || (case22_status == 2)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v - seg_dir[0], 0.0, 0.0, seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightDown Corner
								{
									_face.dir = { seg_dir[0], face_v, face_h, face_h - seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
					}
					else if ((seg_status[2] != -1) && (seg_status[1] != -1)) // Adjacent edges - two cases
					{
						gen_conformal_status(is_unreal, segment[2].Obj_ID, segment[1].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
						if ((seg_status[2] == 1) && (seg_status[1] == 0))
						{
							double area = (face_h - seg_dir[1]) * face_v / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { 0.0, face_v, 0.0, face_h - seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { face_v, 0.0, face_h, seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[1] == 0) && (seg_status[2] == 0))
						{
							double area = seg_dir[2] * (face_h - seg_dir[1]) / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { 0.0, seg_dir[2], 0.0, face_h - seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { face_v, face_v - seg_dir[2], face_h, seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
					}
					else if ((seg_status[3] != -1) && (seg_status[2] != -1)) // Adjacent edges - two cases
					{
						gen_conformal_status(is_unreal, segment[3].Obj_ID, segment[2].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
						if ((seg_status[3] == 1) && (seg_status[2] == 0))
						{
							double area = face_h * (face_v - seg_dir[2]) / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
								{
									_face.dir = { 0.0, face_v - seg_dir[2], face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[2] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v, seg_dir[2], 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[2] == 0) && (seg_status[3] == 0))
						{
							double area = seg_dir[3] * (face_v - seg_dir[2]) / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
								{
									_face.dir = { 0.0, face_v - seg_dir[2], seg_dir[3], 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[2] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v, seg_dir[2], face_h - seg_dir[3], face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
					}
					else if ((seg_status[3] != -1) && (seg_status[0] != -1)) // Adjacent edges - two cases
					{
						gen_conformal_status(is_unreal, segment[3].Obj_ID, segment[0].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
						if ((seg_status[0] == 1) && (seg_status[3] == 0))
						{
							double area = (face_h - seg_dir[3]) * face_v / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { 0.0, face_v, seg_dir[3], face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { face_v, 0.0, face_h - seg_dir[3], 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[3] == 0) && (seg_status[0] == 0))
						{
							double area = seg_dir[0] * (face_h - seg_dir[3]) / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { face_v - seg_dir[0], face_v, seg_dir[3], face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { seg_dir[0], 0.0, face_h - seg_dir[3], 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
					}
				}
				else
				{
					_face.Conformal_status = -1;
					_face.Obj_id.emplace_back(-1);
					// To be determined
				}
				// std::cout << "face_id: " << face_id << " finish" << std::endl;
			}
			else
			{
				_face.Conformal_status = -1;
				// normal cell
			}
		}
	}

	void refine_grid_new(bool is_unreal, TopTools_DataMapOfIntegerShape& occ_shape_list, RTCScene& scene, const std::unordered_map<int, int>& geoId2objId, const std::unordered_map<int, int>& objId2meshmatlId, std::vector<int> PEC_id, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES, std::vector<FACE>& face, int plane)
	{
		double threshold = 4e-6;
		int f = 0, f0 = 0, f2 = 0;
		// unit: ray tracing [mm] filling refine grid [m]
		for (int face_id = 0; face_id < face.size(); face_id++)
		{
			FACE& _face = face[face_id];

			double face_v = _face.v * 1e-3;
			double face_h = _face.h * 1e-3;

			if (_face.flag == 0)
			{
				f++;
				double face_area = face_v * face_h;

				size_t line0 = _face.line_indices[0]; // x line
				size_t line1 = _face.line_indices[1]; // y line
				size_t line2 = _face.line_indices[2]; // z line

				std::array<SEGMENT, 4> segment;
				segment[0].length = _face.h;
				segment[1].length = _face.v;
				segment[2].length = _face.h;
				segment[3].length = _face.v;

				if (plane == 1)
				{
					segment[0].startPoint = { XOY_LINES[line0], YOZ_LINES[line1], ZOX_LINES[line2] };
					segment[1].startPoint = { XOY_LINES[line0 + 1], YOZ_LINES[line1], ZOX_LINES[line2] };
					segment[2].startPoint = { XOY_LINES[line0 + 1], YOZ_LINES[line1 + 1], ZOX_LINES[line2] };
					segment[3].startPoint = { XOY_LINES[line0], YOZ_LINES[line1 + 1], ZOX_LINES[line2] };

					segment[0].direction = { 1, 0, 0 };
					segment[1].direction = { 0, 1, 0 };
					segment[2].direction = { -1, 0, 0 };
					segment[3].direction = { 0, -1, 0 };
				}
				else if (plane == 2)
				{
					segment[0].startPoint = { XOY_LINES[line2], YOZ_LINES[line0], ZOX_LINES[line1] };
					segment[1].startPoint = { XOY_LINES[line2], YOZ_LINES[line0 + 1], ZOX_LINES[line1] };
					segment[2].startPoint = { XOY_LINES[line2], YOZ_LINES[line0 + 1], ZOX_LINES[line1 + 1] };
					segment[3].startPoint = { XOY_LINES[line2], YOZ_LINES[line0], ZOX_LINES[line1 + 1] };

					segment[0].direction = { 0, 1, 0 };
					segment[1].direction = { 0, 0, 1 };
					segment[2].direction = { 0, -1, 0 };
					segment[3].direction = { 0, 0, -1 };
				}
				else if (plane == 3)
				{
					segment[0].startPoint = { XOY_LINES[line1], YOZ_LINES[line2], ZOX_LINES[line0] };
					segment[1].startPoint = { XOY_LINES[line1], YOZ_LINES[line2], ZOX_LINES[line0 + 1] };
					segment[2].startPoint = { XOY_LINES[line1 + 1], YOZ_LINES[line2], ZOX_LINES[line0 + 1] };
					segment[3].startPoint = { XOY_LINES[line1 + 1], YOZ_LINES[line2], ZOX_LINES[line0] };

					segment[0].direction = { 0, 0, 1 };
					segment[1].direction = { 1, 0, 0 };
					segment[2].direction = { 0, 0, -1 };
					segment[3].direction = { -1, 0, 0 };
				}

				// 非边界情况（射线入射与某条几何边完全重合导致交点数异常）
				for (int i = 0; i < 4; i++)
				{
					// occ ray tracing

					std::vector<std::pair<int, gp_Pnt>> occ_objId2Points = get_occ_shape_ray_intersection_points_new(occ_shape_list, segment[i].startPoint, segment[i].direction);

					// mesh ray tracing
					std::vector<std::pair<int, gp_Pnt>> mesh_obj2point;
					RTCRayHit rayhit;
					bool intersect = true;
					double distance = 0.0;
					float ox = segment[i].startPoint.X();
					float oy = segment[i].startPoint.Y();
					float oz = segment[i].startPoint.Z();
					while (intersect == true)
					{
						gen_rayhit(rayhit, ox, oy, oz, segment[i].direction);
						rtcIntersect1(scene, &rayhit);
						if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
						{
							intersect = true;
							distance = rayhit.ray.tfar + threshold;
							ox = rayhit.ray.org_x + distance * rayhit.ray.dir_x;
							oy = rayhit.ray.org_y + distance * rayhit.ray.dir_y;
							oz = rayhit.ray.org_z + distance * rayhit.ray.dir_z;

							int obj_id;
							if (rayhit.hit.geomID == -1)
							{
								obj_id = -1;
							}
							for (auto iter = geoId2objId.begin(); iter != geoId2objId.end(); ++iter)
							{
								if (iter->first == rayhit.hit.geomID)
								{
									obj_id = iter->second;
								}
							}
							mesh_obj2point.emplace_back(std::make_pair(obj_id, gp_Pnt{ ox, oy, oz }));
						}
						else
						{
							intersect = false;
						}
					}

					// gen total ray seg_points info
					int seg_cross = 0;	 // 当前线段上的交点数
					int total_cross = 0; // 射线与几何体的总交点数
					std::vector<std::pair<int, gp_Pnt>> total_obj2point;
					total_obj2point = sort_ray_intersection_points(occ_objId2Points, mesh_obj2point, segment[i], total_cross, seg_cross);

					// gen segment crosspoint status
					segment[i].cross_num = seg_cross;
					if (seg_cross == 0)
					{
						segment[i].status.emplace_back(-1);
						segment[i].cross_dir = std::vector<double>{ 0.0 };
					}
					else if (seg_cross == 1)
					{
						if ((total_cross == seg_cross)) // case1. single cross [inside obj -> airbox]
						{
							segment[i].Obj_ID.emplace_back(total_obj2point[0].first);
							segment[i].Obj_ID.emplace_back(-1);
							segment[i].status.emplace_back(1);
						}
						else if ((total_obj2point[0].first != total_obj2point[1].first)) // case2. multi cross
						{
							if ((total_cross > 2) && (total_obj2point[1].first == total_obj2point[2].first)) // [obj1->air->obj2]
							{
								segment[i].Obj_ID.emplace_back(total_obj2point[0].first);
								segment[i].Obj_ID.emplace_back(-1);
								segment[i].status.emplace_back(1);
							}
							else if (total_cross == 2) // [obj1 -> obj2]
							{
								segment[i].Obj_ID.emplace_back(total_obj2point[0].first);
								segment[i].Obj_ID.emplace_back(total_obj2point[1].first);
								segment[i].status.emplace_back(2);
							}
						}
						else if ((total_obj2point[0].first == total_obj2point[1].first)) // case3. multi cross [airbox -> obj]
						{
							segment[i].Obj_ID.emplace_back(total_obj2point[0].first);
							segment[i].status.emplace_back(0);
						}
					}
					else if (seg_cross == 2)
					{
						if ((total_cross == seg_cross)) // case1. double cross [air -> obj -> air]
						{
							if ((total_obj2point[0].first != total_obj2point[1].first))
							{
								segment[i].status.emplace_back(-1);
								segment[i].cross_dir = std::vector<double>{ 0.0 };
								std::cout << "status 3 abnormal in seg " << i << " face " << face_id << " plane " << plane << std::endl;
							}
							else
							{
								segment[i].Obj_ID.emplace_back(total_obj2point[0].first);
								segment[i].status.emplace_back(3);
								//segment[i].status.emplace_back(0);
								//segment[i].status.emplace_back(1);
							}
						}
						else if ((total_obj2point[0].first == total_obj2point[1].first)) // air hole in obj [obj -> air -> obj]
						{
							segment[i].status.emplace_back(4);
							//segment[i].status.emplace_back(1);
							//segment[i].status.emplace_back(0);
							segment[i].Obj_ID.emplace_back(total_obj2point[0].first);
						}
						else if ((total_obj2point[0].first != total_obj2point[1].first))// air hole between objs [obj->air->obj2]
						{
							segment[i].status.emplace_back(5);
							segment[i].Obj_ID.emplace_back(total_obj2point[0].first);
							segment[i].Obj_ID.emplace_back(total_obj2point[1].first);
						}

						// count2++;
						//  To be determined
					}
					else // seg_cross > 2
					{
						segment[i].status.emplace_back(-1 * seg_cross);

						// countn++;
						//  To be determined
					}
				}
				// std::cout << "face_id: " << face_id << " mid" << std::endl;
				//  generate face crosspoint info
				int face_cross_num = segment[0].cross_num + segment[1].cross_num + segment[2].cross_num + segment[3].cross_num;
				// calculating material segment area
				if ((face_cross_num == 0) || (face_cross_num == 1))
				{
					f0++;
					_face.Conformal_status = -1;
					_face.Obj_id.emplace_back(-1);
				}
				else if (face_cross_num == 2) // two segment areas
				{
					f2++;
					// gen_face_conformal_status_2(_face, is_unreal, segment, PEC_id);
					std::array<double, 4> seg_dir = { segment[0].cross_dir[0] * 1e-3, segment[1].cross_dir[0] * 1e-3, segment[2].cross_dir[0] * 1e-3, segment[3].cross_dir[0] * 1e-3 };
					std::array<int, 4> seg_status; // crosspoint status: -1 no cross, 0 regular crosss, 1 corner point
					std::array<int, 4> dir_status; // direction status: -1 no cross, 0 case3, 1 case1, 2 case2
					for (int seg = 0; seg < 4; seg++)
					{
						dir_status[seg] = segment[seg].status[0];
						if (seg_dir[seg] == 0)
						{
							seg_status[seg] = -1;
						}
						else if (seg_dir[seg] == segment[seg].length)
						{
							seg_status[seg] = 1;
						}
						else
						{
							seg_status[seg] = 0;
						}
					}
					_face.Obj_id.resize(face_cross_num);

					int case22_status = -1;
					// dir:-x +x -y +y / -y +y -z +z / -z +z -x +x
					if ((seg_status[3] != -1) && (seg_status[1] != -1)) // diagnal - four cases
					{
						int case_status = -1;
						double area = 0.0;
						// PEC - Conformal or Medium - Conformal
						gen_conformal_status_2points(is_unreal, segment[3].Obj_ID, segment[1].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
						if ((seg_status[3] == 1) && (seg_status[1] == 1))
						{
							area = face_area / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
								{
									_face.dir = { 0.0, face_v, face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v, 0.0, 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[3] == 1) && (seg_status[1] == 0))
						{
							area = face_v * seg_dir[1] / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
								{
									_face.dir = { 0.0, face_v, face_h, face_h - seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
								else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v, 0.0, 0.0, seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[1] == 1) && (seg_status[3] == 0))
						{
							double area = face_v * (face_h - seg_dir[3]) / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
								{
									_face.dir = { 0.0, face_v, seg_dir[3], 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v, 0.0, face_h - seg_dir[3], face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[3] == 0) && (seg_status[1] == 0))
						{
							double area = face_v * (face_h - seg_dir[3] + seg_dir[1]) / 2.0; // Down seg area
							if (_face.Conformal_status == 1)								 // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - Down Corner
								{
									_face.dir = { 0.0, face_v, seg_dir[3], face_h - seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
								else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - Up Corner
								{
									_face.dir = { face_v, 0.0, face_h - seg_dir[3], seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
					}
					else if ((seg_status[2] != -1) && (seg_status[0] != -1)) // diagnal - four cases
					{
						gen_conformal_status_2points(is_unreal, segment[2].Obj_ID, segment[0].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
						if ((seg_status[2] == 1) && (seg_status[0] == 1))
						{
							double area = face_v * face_h / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { 0.0, face_v, 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { face_v, 0.0, face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[2] == 1) && (seg_status[0] == 0))
						{
							double area = seg_dir[0] * face_h / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { face_v - seg_dir[0], face_v, 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { seg_dir[0], 0.0, face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[0] == 1) && (seg_status[2] == 0))
						{
							double area = seg_dir[2] * face_h / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { 0.0, seg_dir[2], 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { face_v, face_v - seg_dir[2], face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[2] == 0) && (seg_status[0] == 0))
						{
							double area = face_h * (face_v - seg_dir[2] + seg_dir[0]) / 2.0; // 左边面积
							if (_face.Conformal_status == 1)								 // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - Left Corner
								{
									_face.dir = { face_v - seg_dir[0], seg_dir[2], 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - Right Corner
								{
									_face.dir = { seg_dir[0], face_v - seg_dir[2], face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
					}
					else if ((seg_status[1] != -1) && (seg_status[0] != -1)) // Adjacent edges - two cases
					{
						gen_conformal_status_2points(is_unreal, segment[1].Obj_ID, segment[0].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
						if ((seg_status[1] == 1) && (seg_status[0] == 0))
						{
							double area = (face_v - seg_dir[0]) * face_h / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[1] == 0) || (case22_status == 2)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v - seg_dir[0], 0.0, 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightDown Corner
								{
									_face.dir = { seg_dir[0], face_v, face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[0] == 0) && (seg_status[1] == 0))
						{
							double area = (face_v - seg_dir[0]) * seg_dir[1] / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[1] == 0) || (case22_status == 2)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v - seg_dir[0], 0.0, 0.0, seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightDown Corner
								{
									_face.dir = { seg_dir[0], face_v, face_h, face_h - seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
					}
					else if ((seg_status[2] != -1) && (seg_status[1] != -1)) // Adjacent edges - two cases
					{
						gen_conformal_status_2points(is_unreal, segment[2].Obj_ID, segment[1].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
						if ((seg_status[2] == 1) && (seg_status[1] == 0))
						{
							double area = (face_h - seg_dir[1]) * face_v / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { 0.0, face_v, 0.0, face_h - seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { face_v, 0.0, face_h, seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[1] == 0) && (seg_status[2] == 0))
						{
							double area = seg_dir[2] * (face_h - seg_dir[1]) / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[2] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { 0.0, seg_dir[2], 0.0, face_h - seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[1] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { face_v, face_v - seg_dir[2], face_h, seg_dir[1] };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
					}
					else if ((seg_status[3] != -1) && (seg_status[2] != -1)) // Adjacent edges - two cases
					{
						gen_conformal_status_2points(is_unreal, segment[3].Obj_ID, segment[2].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
						if ((seg_status[3] == 1) && (seg_status[2] == 0))
						{
							double area = face_h * (face_v - seg_dir[2]) / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
								{
									_face.dir = { 0.0, face_v - seg_dir[2], face_h, 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[2] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v, seg_dir[2], 0.0, face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[2] == 0) && (seg_status[3] == 0))
						{
							double area = seg_dir[3] * (face_v - seg_dir[2]) / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - RightDown Corner
								{
									_face.dir = { 0.0, face_v - seg_dir[2], seg_dir[3], 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
								else if ((dir_status[2] == 0) || (case22_status == 1)) // PEC - LeftUp Corner
								{
									_face.dir = { face_v, seg_dir[2], face_h - seg_dir[3], face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { area, face_area - area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
					}
					else if ((seg_status[3] != -1) && (seg_status[0] != -1)) // Adjacent edges - two cases
					{
						gen_conformal_status_2points(is_unreal, segment[3].Obj_ID, segment[0].Obj_ID, PEC_id, _face.Obj_id, _face.Conformal_status, case22_status);
						if ((seg_status[0] == 1) && (seg_status[3] == 0))
						{
							double area = (face_h - seg_dir[3]) * face_v / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { 0.0, face_v, seg_dir[3], face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { face_v, 0.0, face_h - seg_dir[3], 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
						else if ((seg_status[3] == 0) && (seg_status[0] == 0))
						{
							double area = seg_dir[0] * (face_h - seg_dir[3]) / 2.0;
							if (_face.Conformal_status == 1) // PEC - Conformal
							{
								if (case22_status == 0)
								{
									_face.dir = { 0.0, 0.0, 0.0, 0.0 };
								}
								else if ((dir_status[3] == 0) || (case22_status == 2)) // PEC - LeftDown Corner
								{
									_face.dir = { face_v - seg_dir[0], face_v, seg_dir[3], face_h };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, face_area - area, face_area);
								}
								else if ((dir_status[0] == 0) || (case22_status == 1)) // PEC - RightUp Corner
								{
									_face.dir = { seg_dir[0], 0.0, face_h - seg_dir[3], 0.0 };
									check_conformal_restraint(_face, _face.Obj_id, _face.dir, area, face_area);
								}
							}
							else if (_face.Conformal_status == 2) // Medium - Conformal
							{
								_face.sec_area = { face_area - area, area };
							}
							else if (_face.Conformal_status == -1) // obj num is not match to be determined
							{
								_face.Obj_id = std::vector<int>{ -1 };
							}
						}
					}
				}
				else
				{
					_face.Conformal_status = -1;
					_face.Obj_id.emplace_back(-1);
					// To be determined
				}
				// std::cout << "face_id: " << face_id << " finish" << std::endl;
			}
			else
			{
				_face.Conformal_status = -1;
				// normal cell
			}
		}
	}

	void gen_face_material(bool is_unreal, std::unordered_map<int, int> objId2meshmatlId, std::vector<FACE>& face, Mesh _mesh)
	{
		double PI = 3.141592653589793;
		double Mu = 4 * PI * 1e-7;
		double DIE = 8.854187817 * 1e-12;
		double z0 = sqrt(Mu / DIE);

		for (int m = 0; m < face.size(); m++)
		{
			FACE& _face = face[m];
			if (_face.Conformal_status == -1)
			{
				for (auto iter = objId2meshmatlId.begin(); iter != objId2meshmatlId.end(); ++iter)
				{
					if (iter->first == _face.Obj_id[0])
					{
						_face.material_idx = iter->second;
					}
				}
			}
			else if (_face.Conformal_status == 1)
			{
				double epsilon, mu, sigma_mu;
				H_inf uni_Hinf;
				for (int a = 0; a < 4; a++)
				{
					uni_Hinf.l[a] = _face.dir[a];
				}
				for (auto iter = objId2meshmatlId.begin(); iter != objId2meshmatlId.end(); ++iter)
				{
					if (iter->first == _face.Obj_id[1])
					{
						epsilon = _mesh.epsilon[iter->second] * DIE;
						mu = _mesh.mu[iter->second] * Mu;
						sigma_mu = _mesh.sigma_mu[iter->second];
					}
				}

				double temp = sigma_mu * _mesh.dt / mu / 2.0;
				uni_Hinf.CP = (-2 * temp) / (1 + temp);
				uni_Hinf.CQ_S1 = (_mesh.dt * z0 / mu) / (1 + temp);
				_mesh.h_inf.emplace_back(uni_Hinf);
				_face.material_idx = (-1) * (_mesh.h_inf.size() - 1);
			}
			else if (_face.Conformal_status == 2)
			{
				double epsilon = 0.0, sigma = 0.0, mu = 0.0, sigma_mu = 0.0;
				double face_area = _face.v * _face.h * 1e-6;
				for (int a = 0; a < _face.Obj_id.size(); a++)
				{
					for (auto iter = objId2meshmatlId.begin(); iter != objId2meshmatlId.end(); ++iter)
					{
						if (iter->first == _face.Obj_id[a])
						{
							epsilon += _face.sec_area[a] * _mesh.epsilon[iter->second];
							sigma += _face.sec_area[a] * _mesh.sigma[iter->second];
							mu += _face.sec_area[a] * _mesh.mu[iter->second];
							sigma_mu += _face.sec_area[a] * _mesh.sigma_mu[iter->second];
						}
					}
				}
				epsilon = epsilon / face_area;
				sigma = sigma / face_area;
				mu = mu / face_area;
				sigma_mu = sigma_mu / face_area;

				_mesh.epsilon.emplace_back(epsilon);
				_mesh.sigma.emplace_back(sigma);
				_mesh.mu.emplace_back(mu);
				_mesh.sigma_mu.emplace_back(sigma_mu);

				epsilon = epsilon * DIE;
				mu = mu * Mu;

				if (is_unreal == true)
				{
					double temp1 = sigma * _mesh.dt / epsilon / 2.0;
					_mesh.CA.emplace_back((-2 * temp1) / (1 + temp1));
					_mesh.CB.emplace_back((_mesh.dt / epsilon / z0) / (1 + temp1));
					_face.material_idx = _mesh.CA.size() - 1;
				}
				else if (is_unreal == false)
				{
					double temp2 = sigma_mu * _mesh.dt / mu / 2.0;
					_mesh.CP.emplace_back((-2 * temp2) / (1 + temp2));
					_mesh.CQ.emplace_back((_mesh.dt * z0 / mu) / (1 + temp2));
					_face.material_idx = _mesh.CP.size() - 1;
				}
			}
		}
	}

	void gen_field_index_parallel(PARALLEL_INF& _PARALLEL_INF, std::vector<FACE>& XOY_FACES, std::vector<FACE>& YOZ_FACES, std::vector<FACE>& ZOX_FACES, std::vector<FACE>& UNREAL_XOY_FACES, std::vector<FACE>& UNREAL_YOZ_FACES, std::vector<FACE>& UNREAL_ZOX_FACES, Mesh& _mesh, Field& _Field, int MPI_ip)
	{
		int MPI_size = _PARALLEL_INF.MPI_Lx * _PARALLEL_INF.MPI_Ly * _PARALLEL_INF.MPI_Lz;
		int x_size, y_size, z_size;

		for (int ip = 0; ip < MPI_size; ip++)
		{
			if (MPI_ip == ip)
			{
				_Field.resize(_PARALLEL_INF.N_x, _PARALLEL_INF.N_y, _PARALLEL_INF.N_z);

				for (int select_EH = 0; select_EH < 2; select_EH++)
				{
					for (int plane = 1; plane < 4; plane++)
					{
						int material_idx;
						if (select_EH == 0)
						{
							if (plane == 1) //_face Ez (NX+1)*(NY+1)*NZ
							{
								x_size = _PARALLEL_INF.N_x;
								y_size = _PARALLEL_INF.N_y;
								z_size = _PARALLEL_INF.N_z + 1;
							}
							else if (plane == 2) //_face Ex NX*(NY+1)*(NZ+1)
							{
								x_size = _PARALLEL_INF.N_x + 1;
								y_size = _PARALLEL_INF.N_y;
								z_size = _PARALLEL_INF.N_z;
							}
							else if (plane == 3) //_face Ey (NX+1)*NY*(NZ+1)
							{
								x_size = _PARALLEL_INF.N_x;
								y_size = _PARALLEL_INF.N_y + 1;
								z_size = _PARALLEL_INF.N_z;
							}
						}
						else if (select_EH == 1)
						{
							if (plane == 1) //_face Ez (NX+1)*(NY+1)*NZ
							{
								x_size = _PARALLEL_INF.N_x + 1;
								y_size = _PARALLEL_INF.N_y + 1;
								z_size = _PARALLEL_INF.N_z;
							}
							else if (plane == 2) //_face Ex NX*(NY+1)*(NZ+1)
							{
								x_size = _PARALLEL_INF.N_x;
								y_size = _PARALLEL_INF.N_y + 1;
								z_size = _PARALLEL_INF.N_z + 1;
							}
							else if (plane == 3) //_face Ey (NX+1)*NY*(NZ+1)
							{
								x_size = _PARALLEL_INF.N_x + 1;
								y_size = _PARALLEL_INF.N_y;
								z_size = _PARALLEL_INF.N_z + 1;
							}
						}

						for (int i = 0; i < x_size; i++)
						{
							int id_x = i + _PARALLEL_INF.N_start_x;
							for (int j = 0; j < y_size; j++)
							{
								int id_y = j + _PARALLEL_INF.N_start_y;
								for (int k = 0; k < z_size; k++)
								{
									int id_z = k + _PARALLEL_INF.N_start_z;
									int id;
									if (plane == 1)
									{
										id = id_x + id_y * x_size + id_z * (x_size * y_size);
									}
									else if (plane == 2)
									{
										id = id_y + id_z * y_size + id_x * (y_size * z_size);
									}
									else if (plane == 3)
									{
										id = id_z + id_x * z_size + id_y * (z_size * x_size);
									}

									if (select_EH == 0)
									{
										if (plane == 1)
										{
											material_idx = XOY_FACES[id].material_idx;
											_Field.index_Hz[i + 1][j + 1][k] = material_idx;
										}
										else if (plane == 2)
										{
											material_idx = YOZ_FACES[id].material_idx;
											_Field.index_Hx[i][j + 1][k + 1] = material_idx;
										}
										else if (plane == 3)
										{
											material_idx = ZOX_FACES[id].material_idx;
											_Field.index_Hy[i + 1][j][k + 1] = material_idx;
										}
									}
									else if (select_EH == 1)
									{
										if (plane == 1)
										{
											material_idx = UNREAL_XOY_FACES[id].material_idx;
											_Field.index_Ez[i][j][k] = material_idx;
										}
										else if (plane == 2)
										{
											material_idx = UNREAL_YOZ_FACES[id].material_idx;
											_Field.index_Ex[i][j][k] = material_idx;
										}
										else if (plane == 3)
										{
											material_idx = UNREAL_ZOX_FACES[id].material_idx;
											_Field.index_Ey[i][j][k] = material_idx;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	void write_bin(Mesh& _mesh, Field& _field, std::string file_path)
	{
		// open bin file & write data
		std::ofstream outFile(file_path, std::ios::binary);
		// std::ofstream outFile(file_path, std::ios::out);
		//  check openfile status
		if (!outFile)
		{
			//Logger::error("Cannot open file for writing!");
			return;
		}
		//_mesh.layer_for_PML = 6;
		// outFile.write((char*)&_mesh.layer_for_PML, 4);

		for (int i = 0; i < 6; i++)
		{
			outFile.write((char*)&_mesh.boundary_layer[i], 4);
		}
		int size = _mesh.epsilon.size();
		outFile.write((char*)&size, 4);
		for (int i = 0; i < size; i++)
		{
			outFile.write((char*)&_mesh.epsilon[i], 8);
			outFile.write((char*)&_mesh.sigma[i], 8);
			outFile.write((char*)&_mesh.mu[i], 8);
			outFile.write((char*)&_mesh.sigma_mu[i], 8);
		}
		size = _mesh.CA.size();
		outFile.write((char*)&size, 4);
		for (int i = 0; i < size; i++)
		{
			outFile.write((char*)&_mesh.CA[i], 8);
			outFile.write((char*)&_mesh.CB[i], 8);
		}
		size = _mesh.CP.size();
		outFile.write((char*)&size, 4);
		for (int i = 0; i < size; i++)
		{
			outFile.write((char*)&_mesh.CP[i], 8);
			outFile.write((char*)&_mesh.CQ[i], 8);
		}
		outFile.write((char*)&_mesh.NX, 4);
		outFile.write((char*)&_mesh.NY, 4);
		outFile.write((char*)&_mesh.NZ, 4);
		for (int i = 0; i < _mesh.NX; i++)
		{
			outFile.write((char*)&_mesh.dx[i], 8);
		}
		for (int j = 0; j < _mesh.NY; j++)
		{
			outFile.write((char*)&_mesh.dy[j], 8);
		}
		for (int k = 0; k < _mesh.NZ; k++)
		{
			outFile.write((char*)&_mesh.dz[k], 8);
		}
		int hinf_size = _mesh.h_inf.size();
		outFile.write((char*)&hinf_size, 4);
		for (int i = 0; i < hinf_size; i++)
		{
			outFile.write((char*)&_mesh.h_inf[i].CP, 8);
			outFile.write((char*)&_mesh.h_inf[i].CQ_S1, 8);
			for (int j = 0; j < 4; j++)
			{
				outFile.write((char*)&_mesh.h_inf[i].l[j], 8);
			}
		}
		for (int i = 0; i < _field.index_Ex.size(); i++)
		{
			for (int j = 0; j < _field.index_Ex[i].size(); j++)
			{
				for (int k = 0; k < _field.index_Ex[i][j].size(); k++)
				{
					outFile.write((char*)&_field.index_Ex[i][j][k], 8);
				}
			}
		}

		for (int i = 0; i < _field.index_Ey.size(); i++)
		{
			for (int j = 0; j < _field.index_Ey[i].size(); j++)
			{
				for (int k = 0; k < _field.index_Ey[i][j].size(); k++)
				{
					outFile.write((char*)&_field.index_Ey[i][j][k], 8);
				}
			}
		}

		for (int i = 0; i < _field.index_Ez.size(); i++)
		{
			for (int j = 0; j < _field.index_Ez[i].size(); j++)
			{
				for (int k = 0; k < _field.index_Ez[i][j].size(); k++)
				{
					outFile.write((char*)&_field.index_Ez[i][j][k], 8);
				}
			}
		}
		for (int i = 0; i < _field.index_Hx.size(); i++)
		{
			for (int j = 0; j < _field.index_Hx[i].size(); j++)
			{
				for (int k = 0; k < _field.index_Hx[i][j].size(); k++)
				{
					outFile.write((char*)&_field.index_Hx[i][j][k], 8);
				}
			}
		}

		for (int i = 0; i < _field.index_Hy.size(); i++)
		{
			for (int j = 0; j < _field.index_Hy[i].size(); j++)
			{
				for (int k = 0; k < _field.index_Hy[i][j].size(); k++)
				{
					outFile.write((char*)&_field.index_Hy[i][j][k], 8);
				}
			}
		}

		for (int i = 0; i < _field.index_Hz.size(); i++)
		{
			for (int j = 0; j < _field.index_Hz[i].size(); j++)
			{
				for (int k = 0; k < _field.index_Hz[i][j].size(); k++)
				{
					outFile.write((char*)&_field.index_Hz[i][j][k], 8);
				}
			}
		}
		outFile.close();
	}

}