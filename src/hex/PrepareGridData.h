#ifndef __PREPAREGRIDDATA__
#define __PREPAREGRIDDATA__

#include <fstream>
#include <string>
#include <vector>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <embree4/rtcore.h>
#include <embree4/rtcore_ray.h>
#include <BRepExtrema_DistShapeShape.hxx>
#include <gp_Lin.hxx>
#include <gp_Pnt.hxx>
#include <vector>
#include <algorithm>
#include <gp_Pnt.hxx>
#include <omp.h>

#include "configor/json.hpp"

#include <TopTools_DataMapOfIntegerShape.hxx>

#include <TopoDS_Shape.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <gp_Lin.hxx>
#include <gp_Pnt.hxx>
#include <Geom_Line.hxx>
#include <GeomAPI_IntCS.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>

#include <iostream>
#include <embree4/rtcore.h>
#include <embree4/rtcore_ray.h>
#include <math.h>

typedef double AX;

namespace HW_PGD {

	const int INIT_FLAG = 99999;

	struct FACE
	{
		std::array<size_t, 3> line_indices;
		double v;
		double h;
		std::vector<int> Obj_id;
		int flag; // 1 outer, 0 middle, -1 inner

		int Conformal_status = -2;
		std::vector<double> dir;
		std::vector<double> sec_area;
		int material_idx;
	};

	struct GEO_USER_DATA
	{
		double Xmin;
		double Xmax;
		double Ymin;
		double Ymax;
		double Zmin;
		double Zmax;
	};

	struct SEGMENT
	{
		gp_Pnt startPoint;
		gp_Dir direction; // (1,0,0) of (0,1,0) or (0,0,1)
		double length;

		int cross_num;
		std::vector<int> status;
		//0: multi&1 [airbox -> obj],1: 1&1 [inside obj -> airbox],2: multi&1 [obj1 -> obj2]
		//3: 2&2 [air -> obj -> air],4: multi&2 [obj -> air -> obj], 5: multi&2 [obj1->air->obj2]
		std::vector<int> Obj_ID;
		//std::vector<gp_Pnt> cross_vtx; 
		std::vector<double> cross_dir;
	};

	// PEC conformal
	struct H_inf
	{
		AX CP;
		AX CQ_S1;
		AX l[4];
	};

	struct PARALLEL_INF
	{
		//节点网格中x,y,z方向的长度
		int MPI_Lx;
		int MPI_Ly;
		int MPI_Lz;

		//本节点MPI_ip所在的位置
		int MPI_x;
		int MPI_y;
		int MPI_z;

		//本节点分配的网格在x/y/z维度的数量
		long long N_x;
		long long N_y;
		long long N_z;

		//本节点分配的网格在x/y/z维度的起始位置
		int N_start_x;
		int N_start_y;
		int N_start_z;

		//其他参数
		int N_iter;
		AX res;//单位线性值
	};

	struct Mesh
	{
		int layer_for_PML;
		//int layer_for_PML[6];
		int boundary_layer[6]; // 边界层数 -x +x -y +y -z +z
		std::vector<AX> epsilon;
		std::vector<AX> sigma;
		std::vector<AX> mu;
		std::vector<AX> sigma_mu;
		std::vector<AX> CA;
		std::vector<AX> CB;
		std::vector<AX> CP;
		std::vector<AX> CQ;
		std::vector<AX> dx;
		std::vector<AX> dy;
		std::vector<AX> dz;

		int NX;
		int NY;
		int NZ;
		AX dt;

		std::vector<H_inf> h_inf;
	};

	struct Field
	{
		std::vector<std::vector<std::vector<long long>>> index_Ex;
		std::vector<std::vector<std::vector<long long>>> index_Ey;
		std::vector<std::vector<std::vector<long long>>> index_Ez;
		std::vector<std::vector<std::vector<long long>>> index_Hx;
		std::vector<std::vector<std::vector<long long>>> index_Hy;
		std::vector<std::vector<std::vector<long long>>> index_Hz;

		void resize(int Nx, int Ny, int Nz)
		{
			index_Ex.resize(Nx);
			for (int i = 0; i < Nx; i++)
			{
				index_Ex[i].resize(Ny + 1);
				for (int j = 0; j < Ny + 1; j++)
				{
					index_Ex[i][j].resize(Nz + 1);
				}
			}

			index_Ey.resize(Nx + 1);
			for (int i = 0; i < Nx + 1; i++)
			{
				index_Ey[i].resize(Ny);
				for (int j = 0; j < Ny; j++)
				{
					index_Ey[i][j].resize(Nz + 1);
				}
			}

			index_Ez.resize(Nx + 1);
			for (int i = 0; i < Nx + 1; i++)
			{
				index_Ez[i].resize(Ny + 1);
				for (int j = 0; j < Ny + 1; j++)
				{
					index_Ez[i][j].resize(Nz);
				}
			}

			//---------------------------------------------

			index_Hx.resize(Nx + 1);
			for (int i = 0; i < Nx + 1; i++)
			{
				index_Hx[i].resize(Ny + 1);
				for (int j = 0; j < Ny + 1; j++)
				{
					index_Hx[i][j].resize(Nz + 1);
				}
			}

			index_Hy.resize(Nx + 1);
			for (int i = 0; i < Nx + 1; i++)
			{
				index_Hy[i].resize(Ny + 1);
				for (int j = 0; j < Ny + 1; j++)
				{
					index_Hy[i][j].resize(Nz + 1);
				}
			}

			index_Hz.resize(Nx + 1);
			for (int i = 0; i < Nx + 1; i++)
			{
				index_Hz[i].resize(Ny + 1);
				for (int j = 0; j < Ny + 1; j++)
				{
					index_Hz[i][j].resize(Nz + 1);
				}
			}
		}
	};


	// gen fdtd prepare data
	void gen_fdtd_prepare_data(const std::string& prepare_data_dir);

	// write fdtd prepare data
	bool read_fdtd_prepare_data(TopTools_DataMapOfIntegerShape& occ_shape_list, RTCScene& scene, RTCDevice& device, std::unordered_map<int, int>& geoid_2_objid, std::unordered_map<int, std::string>& objid_2_mkid, configor::json& materials, std::vector<double>& XOY_LINES, std::vector<double>& YOZ_LINES, std::vector<double>& ZOX_LINES, const std::string& prepare_data_dir, std::array<int, 2>& boundary_level_x, std::array<int, 2>& boundary_level_y, std::array<int, 2>& boundary_level_z, double& min_Len);

	// write mesh 
	RTCGeometry read_mesh_data(RTCDevice& device, const std::string& file_path);

	// unreal lines
	std::vector<double> gen_unreal_lines(const std::vector<double>& origin_lines);

	// gen face
	bool gen_mesh_faces(std::vector<FACE>& XOY_FACES, std::vector<FACE>& YOZ_FACES, std::vector<FACE>& ZOX_FACES, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES);

	// get xoy face center
	gp_Pnt get_xoy_face_center(const FACE& face, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES);
	// get yoz face center
	gp_Pnt get_yoz_face_center(const FACE& face, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES);
	// get zox face center
	gp_Pnt get_zox_face_center(const FACE& face, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES);

	// get occ shape ray intersection points
	std::vector<gp_Pnt> get_occ_shape_ray_intersection_points(const TopoDS_Shape& shape, const gp_Pnt& originPnt, const gp_Dir& dir);
	std::unordered_map<int, std::vector<gp_Pnt>> get_occ_shape_ray_intersection_points(const TopTools_DataMapOfIntegerShape& occ_shape_list, const gp_Pnt& originPnt, const gp_Dir& dir);

	// get mesh ray intersection points
	std::unordered_map<int, std::vector<gp_Pnt>> get_mesh_ray_intersection_points(const RTCScene& scene, const gp_Pnt& originPnt, const gp_Dir& dir);

	// handle face material
	bool handle_faces_material(const TopTools_DataMapOfIntegerShape& occ_shape_list, const RTCScene& scene, std::vector<FACE>& XOY_FACES, std::vector<FACE>& YOZ_FACES, std::vector<FACE>& ZOX_FACES, const std::unordered_map<int, int>& geoId2objId, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES, bool is_unreal);

	// set face belong to the obj
	void set_face_belong_to_obj(FACE& face, int obj_id);

	// zxm
	// generate obj_id with PEC
	std::vector<int> gen_PEC_obj_id(std::unordered_map<int, std::string> objId2matId, std::unordered_map<std::string, configor::json>& materials);
	// generate mesh basic info
	void gen_mesh_info(Mesh& _mesh, int XOY_LINES_COUNT, int YOZ_LINES_COUNT, int ZOX_LINES_COUNT, std::vector<FACE> xoy_face, std::vector<FACE> yoz_face, std::vector<FACE> zox_face, std::array<int, 2> boundary_level_x, std::array<int, 2> boundary_level_y, std::array<int, 2> boundary_level_z, double min_cell_size);
	// generate basic objId2meshmatlId for Mesh
	std::unordered_map<int, int> gen_objId2meshmatlId(std::unordered_map<std::string, configor::json> materials, std::unordered_map<int, std::string> objId2matId, Mesh& _mesh);
	// refine mesh XOY, YOZ, ZOX faces
	void refine_grid(bool is_unreal, TopTools_DataMapOfIntegerShape& occ_shape_list, RTCScene& scene, const std::unordered_map<int, int>& geoId2objId, const std::unordered_map<int, int>& objId2meshmatlId, std::vector<int> PEC_Obj_id, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES, std::vector<FACE>& _face, int plane); // plane:1-XOY 2-YOZ 3-ZOX
	void refine_grid_new(bool is_unreal, TopTools_DataMapOfIntegerShape& occ_shape_list, RTCScene& scene, const std::unordered_map<int, int>& geoId2objId, const std::unordered_map<int, int>& objId2meshmatlId, std::vector<int> PEC_Obj_id, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES, std::vector<FACE>& _face, int plane); // plane:1-XOY 2-YOZ 3-ZOX
	void gen_face_conformal_status_2(FACE& _face, bool is_unreal, std::array<SEGMENT, 4>& segment, std::vector<int> PEC_id);

	// subfunction of refine_grid: generate single ray for embree
	void gen_rayhit(RTCRayHit& rayhit, float ox, float oy, float oz, gp_Dir direction, float t_near = 0);
	// get occ_shape_ray_intersection_points
	std::vector<gp_Pnt> get_occ_shape_ray_intersection_points_new(const TopoDS_Shape& shape, const gp_Pnt& originPnt, const gp_Dir& dir);
	std::vector<std::pair<int, gp_Pnt>> get_occ_shape_ray_intersection_points_new(TopTools_DataMapOfIntegerShape& occ_shape_list, const gp_Pnt& originPnt, const gp_Dir& dir);
	// sort intersection_points from embree & occ
	std::vector<std::pair<int, gp_Pnt>> sort_ray_intersection_points(std::vector<std::pair<int, gp_Pnt>> occ_objId2Points, std::vector<std::pair<int, gp_Pnt>> mesh_objId2Points, SEGMENT& _segment, int& total_cross, int& seg_cross);
	// generate face conformal status
	void gen_conformal_status(bool is_unreal, std::vector<int> obj_id1, std::vector<int> obj_id2, std::vector<int> PEC_id, std::vector<int>& obj_id, int& conformal_status, int& case22_status);
	void gen_conformal_status_2points(bool is_unreal, std::vector<int> obj_id1, std::vector<int> obj_id2, std::vector<int> PEC_id, std::vector<int>& obj_id, int& conformal_status, int& case22_status);

	// check_conformal_restraint
	void check_conformal_restraint(FACE& _face, std::vector<int> Obj_id, std::vector<double> dir, double area, double face_area);
	// generate face material_id
	void gen_face_material(bool is_unreal, std::unordered_map<int, int> objId2meshmatlId, std::vector<FACE>& _face, Mesh _mesh);
	// generate unreal face (remove the first and last plane)
	bool gen_unreal_mesh_faces(std::vector<FACE>& XOY_FACES, std::vector<FACE>& YOZ_FACES, std::vector<FACE>& ZOX_FACES, const std::vector<double>& XOY_LINES, const std::vector<double>& YOZ_LINES, const std::vector<double>& ZOX_LINES);

	// generate field index
	void gen_field_index_parallel(PARALLEL_INF& _PARALLEL_INF, std::vector<FACE>& XOY_FACES, std::vector<FACE>& YOZ_FACES, std::vector<FACE>& ZOX_FACES, std::vector<FACE>& UNREAL_XOY_FACES, std::vector<FACE>& UNREAL_YOZ_FACES, std::vector<FACE>& UNREAL_ZOX_FACES, Mesh& _mesh, Field& _Field, int MPI_ip);

	// test write bin
	void write_bin(Mesh& _mesh, Field& _field, std::string file_path);
}

#endif