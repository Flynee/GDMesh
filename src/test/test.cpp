#define NOMINMAX


#include "test.h"

#include <iostream>
#include <string>
#include <string_view>

#include <eigen/Eigen>
#include <igl/readOFF.h>
#include <igl/readSTL.h>
#include <igl/decimate.h>
#include <igl/opengl/glfw/Viewer.h>

#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/decimate.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/parallel_for.h>
#include <igl/read_triangle_mesh.h>


#include <OSD_OpenFile.hxx>
#include <gp_Vec.hxx>
#include <gp_Pnt.hxx>
#include <STEPControl_Reader.hxx>
#include <TopoDS_Shape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <BRep_Tool.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <STEPControl_Writer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Lin.hxx>
#include <GeomAPI_IntCS.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Line.hxx>
#include <Geom_Surface.hxx>
#include <BRep_Tool.hxx>


#include <chrono>
#include <array>


#include <embree4/rtcore.h>
#include <embree4/rtcore_ray.h>

#define TEST_DATA_PATH "../test_data"



#include <iostream>
#include <cmath>
#include <vector>




void test_start() {
	std::cout << "test start..." << std::endl;

	double L = 19;
	double dmax = 2;
	double max_step = 1;
	double min_step = 0.5;
	double R = 1.5;

	const std::vector<double> steps = dA_to_dM_to_dB(L, dmax, min_step, max_step, R);


	double sum = 0;
	for (size_t i = 0; i < steps.size(); i++)
	{
		std::cout << "step = " << steps[i] << std::endl;
		sum += steps[i];

	}

	std::cout << "sum = " << sum << std::endl;

}

std::vector<double> small_to_big(double L, double dmax, double max_step, double min_step, double R) {


	std::vector<double> steps;


	int NA = Floor(Abs(1 + (log(dmax / min_step) / log(R))));
	double LA = min_step * (1 - Pow(R, NA)) / (1 - R);

	std::cout << "NA : " << NA << std::endl;
	std::cout << "LA : " << LA << std::endl;

	double gap = L - LA - max_step;

	std::cout << "gap : " << gap << std::endl;

	int GB = gap / max_step;
	std::cout << "GB : " << GB << std::endl;

	int split = 5;

	if (GB > split) { // general gen
		int NB = GB + 1;
		double LB = NB * max_step;
		std::cout << "LB : " << LB << std::endl;

		double diff = L - LA - LB;
		std::cout << "diff : " << diff << std::endl;

		for (int i = 0; i < NA; i++) {
			double step = min_step * Pow(R, i);
			steps.push_back(step);
		}

		if (diff >= (1 / R) * max_step){

			for (int i = 0; i < NB + 1; i++) {
				double step = max_step;

				if (i == 0) {
					step = diff;
				}

				steps.push_back(step);

			}
		}
		else {

			diff = max_step - diff;
			for (int i = 0; i < NB + 1; i++) {
				double step = max_step;

				if (i < split) {
					step -= diff / split;
				}

				steps.push_back(step);

			}
		}
		
	}
	else { // fast gen
		R = Min(R, (max_step - min_step) / L + 1);
		std::cout << "R : " << R << std::endl;

		NA = Floor(Abs((log(max_step / min_step) / log(R))));
		LA = min_step * (1 - Pow(R, NA)) / (1 - R);
		std::cout << "NA : " << NA << std::endl;
		std::cout << "LA : " << LA << std::endl;

		double diff = (L - LA) / NA;

		std::cout << "diff : " << diff << std::endl;

		for (int i = 0; i < NA; i++) {
			double step = min_step * Pow(R, i);
			step += diff;

			steps.push_back(step);

		}
		

	}



	return steps;

}

std::vector<double> dA_to_dM_to_dB(double L, double dmax, double dA, double dB, double R) {

	std::vector<double> steps;

	int NA = Floor(Abs(1 + (log(dmax / dA) / log(R))));
	double LA = dA * (1 - Pow(R, NA)) / (1 - R);

	int NB = Floor(Abs(1 + (log(dmax / dB) / log(R))));
	double LB = dB * (1 - Pow(R, NB)) / (1 - R);

	double LM = L - LA - LB;

	int NM = LM / dmax;

	int split = 5;

	if (NM < split) {
		double min_step = Min(dA, dB);
		double max_step = Max(dA, dB);

		R = Min(R, (max_step - min_step) / L + 1);
		std::cout << "R : " << R << std::endl;

		NA = Floor(Abs((log(max_step / min_step) / log(R))));
		LA = min_step * (1 - Pow(R, NA)) / (1 - R);
		std::cout << "NA : " << NA << std::endl;
		std::cout << "LA : " << LA << std::endl;

		double diff = (L - LA) / NA;

		std::cout << "diff : " << diff << std::endl;

		for (int i = 0; i < NA; i++) {
			double step = min_step * Pow(R, i);
			step += diff;

			steps.push_back(step);

		}
		if (dA > dB) {
			std::reverse(steps.begin(), steps.end());
		}
	}
	else {

		// LA
		for (int i = 0; i < NA; i++) {
			double step = dA * Pow(R, i);
			steps.push_back(step);
		}

		// LM
		double diff = LM - (NM * dmax);
		std::cout << "diff : " << diff << std::endl;

		if (diff >= (1 / R) * dmax) {

			for (int i = 0; i < NM + 1; i++) {
				double step = dmax;

				if (i == 0) {
					step = diff;
				}

				steps.push_back(step);

			}
		}
		else {

			diff = dmax - diff;
			for (int i = 0; i < NM + 1; i++) {
				double step = dmax;

				if (i < split) {
					step -= diff / split;
				}

				steps.push_back(step);

			}
		}

		// LB
		for (int i = NB - 1; i >= 0; i--) {
			double step = dB * Pow(R, i);
			steps.push_back(step);
		}

	}

	return steps;

}

void sort_pnt_byZ() {
	double tol = 1e-7;

	gp_Dir dir(0, 0, 1);
	gp_XYZ c = dir.XYZ();
	std::vector<gp_Pnt> pnts = { gp_Pnt(1,0, 3), gp_Pnt(3, 0, 1), gp_Pnt(2,0, 2), gp_Pnt(1, 0, 2), };
	std::sort(pnts.begin(), pnts.end(), [c](const gp_Pnt& p1, const gp_Pnt& p2) {
		gp_XYZ c1 = p1.XYZ();
		gp_XYZ c2 = p2.XYZ();
		return (c1.Dot(c) < c2.Dot(c));
	});

	std::vector<gp_Pnt> resultPnts;
	
	for (int i = 0; i < pnts.size() - 1; i++) {
	
		if (i == 0) {
			resultPnts.push_back(pnts[0]);
		}
		gp_XYZ c1 = resultPnts[resultPnts.size() - 1].XYZ();
		gp_XYZ c2 = pnts[i+1].XYZ();
		if (Abs(c1.Dot(c) - c2.Dot(c)) > tol) {
			resultPnts.push_back(pnts[i+1]);
		}
	}

	for (int i = 0; i < resultPnts.size(); i++) {
		resultPnts[i].DumpJson(std::cout);

	}
}

void dA_dB_dmax() {

	double L = 9;
	double dmax = 2;
	double dB = 2;
	double dA = 2;
	double R = 1.5;

	double line = 0;
	std::vector<double> lines;
	double end = L;


	int N = Floor(L / dA);
	double TL = N * dA;
	double tol = 1e-7;

	if ((L - TL) > tol) {

		N += 1;
		if (N % 2 == 0) { // 偶数个
			std::cout << "偶数" << std::endl;
			double cd = (dA * 2 + (L - TL)) / 3;

			for (int i = 0; i < N; i++) {
				if (i == N / 2 || i == (N / 2 - 1) || i == (N / 2 + 1)) {
					line += cd;
				}
				else {
					line += dA;
				}

				lines.push_back(line);

			}

		}
		else { // 奇数个
			std::cout << "奇数" << std::endl;
			double cd = (dA * 2 + (L - TL)) / 3;

			for (int i = 0; i < N; i++) {
				if (i == (N - 1) / 2 || i == ((N - 1) / 2 - 1) || i == ((N - 1) / 2 + 1)) {
					line += cd;
				}
				else {
					line += dA;
				}

				lines.push_back(line);

			}
		}

	}
	else {

		double d = L / N;

		for (int i = 0; i < N; i++) {
			line += d;

			if (i == N - 1) {
				line = end;
			}
			lines.push_back(line);
		}
	}


	for (size_t i = 0; i < lines.size(); i++)
	{
		std::cout << "line = " << lines[i] << std::endl;

	}

}


void dA_dM_dB() {
	double L = 7.;
	double dmax = 2;
	double dB = 0.5;
	double dA = 2;
	double R = 1.5;

	double line = 0;
	std::vector<double> lines;
	double end = L;


	// double L = step;
	// double R = gRmax;
	
	std::cout << "dM > 0" << std::endl;
	int NA = Floor(Abs(1 + (log(dmax / dA) / log(R))));
	double LA = dA * (1 - Pow(R, NA)) / (1 - R);
	double LastDA = dA * Pow(R, NA - 1);

	int NB = Floor(Abs(1 + (log(dmax / dB) / log(R))));
	double LB = dB * (1 - Pow(R, NB)) / (1 - R);
	double LastDB = dB * Pow(R, NB - 1);

	double LM = L - (LA + LB);


	if (LM > 0) { // LM > 0
	
		std::cout << "LM = " << LM << std::endl;

		if (LM < dA / R || LM < dB / R) {
			std::cout << "LM < dA / R || LM < dB / R" << std::endl;

			double dM = (LastDA + LastDB + LM) / 3;

			// LA
			for (int i = 0; i < NA; i++) {
				double d = dA * Pow(R, i);

				if (i == NA - 1) {
					d = dM;
				}

				line += d;
				lines.push_back(line);
			}

			// LM
			line += dM;
			lines.push_back(line);

			// LB
			for (int i = NB - 1; i >= 0; i--) {
				double d = dB * Pow(R, i);

				if (i == NB - 1) {
					d = dM;
				}
				line += d;
				
				if (i == 0) {
					line = end;
				}

				lines.push_back(line);
			}
		}
		else if (LM <= dmax) { 
			std::cout << "LM <= dmax"<< std::endl;

			double dM = LM;
			// LA
			for (int i = 0; i < NA; i++) {
				double d = dA * Pow(R, i);

				line += d;
				lines.push_back(line);
			}

			// LM
			line += dM;
			lines.push_back(line);

			// LB
			for (int i = NB - 1; i >= 0; i--) {
				double d = dB * Pow(R, i);

				line += d;

				if (i == 0) {
					line = end;
				}

				lines.push_back(line);
			}
		}
		else if (LM <= 2*dmax) {
			std::cout << "LM <= 2*dmax" << std::endl;

			int k = 3 + Floor(LM / LastDA);
			double dM = (LastDA + LastDB + LM) / k;

			// LA
			for (int i = 0; i < NA; i++) {
				double d = dA * Pow(R, i);

				if (i == NA - 1) {
					d = dM;
				}

				line += d;
				lines.push_back(line);
			}

			// LM
			for (int i = 0; i < k - 2; i++) {
				line += dM;
				lines.push_back(line);
			}

			// LB
			for (int i = NB - 1; i >= 0; i--) {
				double d = dB * Pow(R, i);

				if (i == NB - 1) {
					d = dM;
				}
				line += d;

				if (i == 0) {
					line = end;
				}

				lines.push_back(line);
			}
		}
		else {

			// LA
			for (int i = 0; i < NA; i++) {
				double d = dA * Pow(R, i);

				line += d;
				lines.push_back(line);
			}

			// LM
			double dM = dmax;
			int NM = Floor(LM / dM);
			double TL = NM * dM;
			double tol = 1e-7;

			if ((LM - TL) > tol) {

				NM += 1;
				if (NM % 2 == 0) { // 偶数个
					std::cout << "偶数" << std::endl;
					double cd = (dM * 2 + (LM - TL)) / 3;

					for (int i = 0; i < NM; i++) {
						if (i == NM / 2 || i == (NM / 2 - 1) || i == (NM / 2 + 1)) {
							line += cd;
						}
						else {
							line += dM;
						}

						lines.push_back(line);

					}

				}
				else { // 奇数个
					std::cout << "奇数" << std::endl;
					double cd = (dM * 2 + (LM - TL)) / 3;

					for (int i = 0; i < NM; i++) {
						if (i == (NM - 1) / 2 || i == ((NM - 1) / 2 - 1) || i == ((NM - 1) / 2 + 1)) {
							line += cd;
						}
						else {
							line += dM;
						}

						lines.push_back(line);

					}
				}

			}
			else {

				double d = LM / NM;

				for (int i = 0; i < NM; i++) {
					line += d;
					lines.push_back(line);
				}
			}


			// LB
			for (int i = NB - 1; i >= 0; i--) {
				double d = dB * Pow(R, i);

				line += d;

				if (i == 0) {
					line = end;
				}

				lines.push_back(line);
			}
		}


	}
	else { // LM <= 0
		std::cout << "LM = "<< LM << std::endl;
		

		std::cout << "dA = " << dA << std::endl;
		std::cout << "dB = " << dB << std::endl;

		double diff = LM;
		double avg = diff / (NA - 1 + NB - 1);

		for (int i = 0; i < NA; i++) {
			double d = dA * Pow(R, i);

			if (i > 0) {
				d += avg;
			}

			line += d;
			lines.push_back(line);
		}

		// LB
		for (int i = NB - 1; i >= 0; i--) {
			double d = dB * Pow(R, i);

			d += avg;
			line += d;

			if (i == 0) {
				line = end;
			}

			lines.push_back(line);
		}
	}

	double sum = 0;
	for (size_t n = 0; n < lines.size(); n++)
	{
		//std::cout << "line = " << lines[n] << std::endl;

		if (n < lines.size() - 1) {
			double step = lines[n + 1] - lines[n];
			std::cout << "step = " << step << std::endl;

			sum += step;
		}

	}

	std::cout << "sum = " << sum << std::endl;

}




void draw_mesh() {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	// Load a mesh in OFF format
	std::string fileName = TEST_DATA_PATH"/bunny.off";

	std::cout << "fileName: " << fileName << std::endl;
	igl::readOFF(fileName, V, F);

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	viewer.launch();
}


void draw_box() {
	// 网格数据
	Eigen::Matrix<double, 8, 3> V{
		{0.5, 0.5, 0.5 },
		{-0.5, 0.5, 0.5},
		{-0.5, -0.5, 0.5},
		{0.5, -0.5, 0.5},
		{0.5, -0.5, -0.5},
		{0.5, 0.5, -0.5},
		{-0.5, 0.5, -0.5},
		{-0.5, -0.5, -0.5},
	};
	Eigen::Matrix<int, 12, 3> F{
		{0, 1, 2},
		{ 0, 2, 3}, // 前面
		{0, 3, 4},
		{ 0, 4, 5}, // 右面
		{0, 5, 6},
		{ 0, 6, 1}, // 上面
		{1, 6, 7},
		{ 1, 7, 2}, // 左面
		{7, 4, 3},
		{ 7, 3, 2}, // 下面
		{4, 7, 6},
		{ 4, 6, 5} // 后面

	};


	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);

	// 根据面的法线来着色
	viewer.data().set_face_based(true);


	viewer.launch();
}

void read_stl() {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd N;
	// Load a mesh in OFF format
	std::string fileName = TEST_DATA_PATH"/Inlet_binary_3000w.stl";

	//igl::readSTL(fileName, V, F, N);
	FILE* fp = fopen(fileName.c_str(), "rb");

	igl::readSTL(fp, V, F, N);


	std::cout << "V size " << V.size() * sizeof(V(0, 0)) / 1024 / 1024 << " M " << std::endl;
	std::cout << "F size " << F.size() * sizeof(F(0, 0)) / 1024 / 1024 << " M " << std::endl;
	std::cout << "N size " << N.size() * sizeof(N(0, 0)) / 1024 / 1024 << " M " << std::endl;

}

void mesh_decimate() {

	std::string fileName = TEST_DATA_PATH"/Inlet_binary_3000w.stl";
	//igl::readSTL(fileName, V, F, N);
	FILE* fp = fopen(fileName.c_str(), "rb");

	Eigen::MatrixXd V, OV, ON;
	Eigen::MatrixXi F, OF;
	igl::readSTL(fp, OV, OF, ON);

	igl::opengl::glfw::Viewer viewer;

	// Prepare array-based edge data structures and priority queue
	Eigen::VectorXi EMAP;
	Eigen::MatrixXi E, EF, EI;
	igl::min_heap< std::tuple<double, int, int> > Q;
	Eigen::VectorXi EQ;
	// If an edge were collapsed, we'd collapse it to these points:
	Eigen::MatrixXd C;
	int num_collapsed;

	// Function to reset original mesh and data structures
	/*const auto& reset = [&]()
	{*/
	F = OF;
	V = OV;
	igl::edge_flaps(F, E, EMAP, EF, EI);
	C.resize(E.rows(), V.cols());
	Eigen::VectorXd costs(E.rows());
	// https://stackoverflow.com/questions/2852140/priority-queue-clear-method
	// Q.clear();
	Q = {};
	EQ = Eigen::VectorXi::Zero(E.rows());
	{
		Eigen::VectorXd costs(E.rows());
		igl::parallel_for(E.rows(), [&](const int e)
		{
			double cost = e;
			Eigen::RowVectorXd p(1, 3);
			igl::shortest_edge_and_midpoint(e, V, F, E, EMAP, EF, EI, cost, p);
			C.row(e) = p;
			costs(e) = cost;
		}, 1);
		for (int e = 0; e < E.rows(); e++)
		{
			Q.emplace(costs(e), e, 0);
		}
	}


	num_collapsed = 0;
	viewer.data().clear();
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
	//};
}


void read_step() {
	// Load a mesh in OFF format
	std::string fileName = TEST_DATA_PATH"/union-box.step";
	std::string outFileName = TEST_DATA_PATH"/daoru.brep";


	// 读取 step 模型
	STEPControl_Reader aReader;
	IFSelect_ReturnStatus status = aReader.ReadFile(fileName.c_str());
	if (status == IFSelect_RetDone) {
		bool failsonly = false;
		//aReader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity);
		int nbr = aReader.NbRootsForTransfer();
		//aReader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity);
		std::cout << "nbr " << nbr << std::endl;

		for (Standard_Integer n = 1; n <= nbr; n++) {
			bool ok = aReader.TransferRoot(n);
			int nbs = aReader.NbShapes();

			std::cout << "nbs " << nbs << std::endl;

			if (nbs > 0) {
				TopoDS_Shape aShape = aReader.OneShape();
				// 转换为 brep 模型
				BRep_Builder aBuilder;
				TopoDS_Compound aCompound;
				aBuilder.MakeCompound(aCompound);
				aBuilder.Add(aCompound, aShape);
				// 写入本地文件
				BRepTools::Write(aCompound, outFileName.c_str());
			}
		}
	}


}

void write_step() {
	std::string fileName = TEST_DATA_PATH"/union-box.brep";
	std::string outFileName = TEST_DATA_PATH"/daoru.step";

	// 创建一个空的形状
	TopoDS_Shape shape_Brep;
	// 创建一个构建器
	BRep_Builder builder_Brep;
	// 从文件中读取形状
	BRepTools::Read(shape_Brep, fileName.c_str(), builder_Brep);

	// 创建一个写入器
	STEPControl_Writer aWriter;

	// 将 shape 转换并添加到写入器中
	IFSelect_ReturnStatus status = aWriter.Transfer(shape_Brep, STEPControl_AsIs);
	// 检查转换是否成功
	if (status == IFSelect_RetDone) {
		// 将写入器中的数据写入到文件中
		IFSelect_ReturnStatus status3 = aWriter.Write(outFileName.c_str());
		// 检查写入是否成功
		if (status3 != IFSelect_RetDone) {
			std::cout << "Writing error" << std::endl;
		}
	}

}


void ray_mesh() {

	auto startRead = std::chrono::high_resolution_clock::now();
	// include embree header file
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd N;
	// Load a mesh in OFF format
	std::string fileName = TEST_DATA_PATH"/sphere_r10.stl";

	//igl::readSTL(fileName, V, F, N);
	FILE* fp = fopen(fileName.c_str(), "rb");

	igl::readSTL(fp, V, F, N);
	auto endRead = std::chrono::high_resolution_clock::now();

	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(endRead - startRead);
	std::cout << "Read Time: " << elapsed.count() << " ms.\n";



	auto starthandleData = std::chrono::high_resolution_clock::now();


	int rowsV = V.rows();
	int colsV = V.cols();

	// 动态分配内存
	float* vertices = new float[rowsV * colsV];

	// 将顶点数据从MatrixXd复制到float数组
	for (int i = 0; i < rowsV; ++i) {
		for (int j = 0; j < colsV; ++j) {
			vertices[i * colsV + j] = static_cast<float>(V(i, j));
		}
	}

	// 获取矩阵的行数和列数
	int rowsF = F.rows();
	int colsF = F.cols();

	// 创建一个unsigned int数组来存储三角形索引数据
	unsigned int* triangles = new unsigned int[rowsF * colsF];

	// 将三角形索引数据从MatrixXi复制到unsigned int数组
	for (int i = 0; i < rowsF; ++i) {
		for (int j = 0; j < colsF; ++j) {
			triangles[i * colsF + j] = static_cast<unsigned int>(F(i, j));
		}
	}


	auto endhandleData = std::chrono::high_resolution_clock::now();
	auto elapsed1 = std::chrono::duration_cast<std::chrono::milliseconds>(endhandleData - starthandleData);
	std::cout << "Handle Data Time: " << elapsed1.count() << " ms.\n";


	auto startBuildBVH = std::chrono::high_resolution_clock::now();

	// create a new embree device
	RTCDevice device = rtcNewDevice(NULL);
	// create a new embree scene
	RTCScene scene = rtcNewScene(device);
	// create a new embree geometry for the triangle mesh
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
	// set the vertex buffer of the geometry
	rtcSetSharedGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4, vertices, 0, sizeof(float) * 4, 4);
	// set the index buffer of the geometry
	rtcSetSharedGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, triangles, 0, sizeof(unsigned int) * 3, 2);
	// commit the geometry
	rtcCommitGeometry(geom);
	// attach the geometry to the scene
	unsigned int geomID = rtcAttachGeometry(scene, geom);

	// commit the scene
	rtcCommitScene(scene);

	auto endBuildBVH = std::chrono::high_resolution_clock::now();
	auto elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(endBuildBVH - startBuildBVH);
	std::cout << "Build BVH Time: " << elapsed2.count() << " ms.\n";



	/////////////////////////////////////// ray test ///////////////////////////////////////
	// 
	// 在球上随机生成1000个点，并将这些点作为射线的终点

	auto startRayIntersect = std::chrono::high_resolution_clock::now();

	const int numRays = 1;
	std::vector<RTCRayHit> rayHits(numRays);
	// 生成球上的点
	for (int i = 0; i < numRays; ++i) {
		RTCHit hit;
		hit.geomID = geomID;
		hit.primID = i;

		RTCRay ray;



		// 极坐标表示，生成球面上的点
		float theta = static_cast<float>(rand()) / RAND_MAX * 2.0f * M_PI; // 角度范围 [0, 2π]
		float phi = static_cast<float>(rand()) / RAND_MAX * M_PI; // 极角范围 [0, π]

		// 球坐标转换为直角坐标
		float x = 5000.0f * sin(phi) * cos(theta);
		float y = 5000.0f * sin(phi) * sin(theta);
		float z = 5000.0f * cos(phi);

		// 设置射线起点
		ray.org_x = 0;
		ray.org_y = 0;
		ray.org_z = -30;

		// 设置射线方向（向球上的点发射）
		ray.dir_x = 0.;
		ray.dir_y = 0.;
		ray.dir_z = 1.;

		ray.tnear = 0.0f;
		ray.tfar = 10000.;
		ray.mask = 0xFFFFFFFF;
		ray.id = i;
		ray.flags = 0;


		RTCRayHit rayhit;

		rayhit.hit = hit;
		rayhit.ray = ray;
		rtcIntersect1(scene, &rayhit);
		rayHits.push_back(rayhit);
	}

	auto endRayIntersect = std::chrono::high_resolution_clock::now();
	auto elapsed3 = std::chrono::duration_cast<std::chrono::milliseconds>(endRayIntersect - startRayIntersect);
	std::cout << "Ray intersect Time: " << elapsed3.count() << " ms.\n";

	// 处理射线拾取的结果
	for (int i = 0; i < numRays; ++i) {
		if (rayHits[i].hit.geomID != RTC_INVALID_GEOMETRY_ID) {
			std::cout << "Ray " << i << " hit triangle " << rayHits[i].hit.primID << " at distance " << rayHits[i].ray.tfar << std::endl;
		}
		else {
			std::cout << "Ray " << i << " hit nothing." << std::endl;
		}
	}

	/////////////////////////////////////// ray test ///////////////////////////////////////

	// release the embree objects
	rtcReleaseGeometry(geom);
	rtcReleaseScene(scene);
	rtcReleaseDevice(device);
	delete[] vertices;
	delete[] triangles;
}



void ray_mesh2() {

	auto startRead = std::chrono::high_resolution_clock::now();
	// include embree header file
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd N;
	// Load a mesh in OFF format
	std::string fileName = TEST_DATA_PATH"/sphere_r10.stl";

	//igl::readSTL(fileName, V, F, N);
	FILE* fp = fopen(fileName.c_str(), "rb");

	igl::readSTL(fp, V, F, N);
	auto endRead = std::chrono::high_resolution_clock::now();

	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(endRead - startRead);
	std::cout << "Read Time: " << elapsed.count() << " ms.\n";



	auto starthandleData = std::chrono::high_resolution_clock::now();


	int rowsV = V.rows();
	int colsV = V.cols();

	// 动态分配内存
	//float* vertices = new float[rowsV * colsV];

	//std::unique_ptr<float[]> vertices = std::make_unique<float[]>(rowsV * colsV);

	//// 将顶点数据从MatrixXd复制到float数组
	//for (int i = 0; i < rowsV; ++i) {
	//	for (int j = 0; j < colsV; ++j) {
	//		vertices[i * colsV + j] = static_cast<float>(V(i, j));
	//	}
	//}

	// 获取矩阵的行数和列数
	int rowsF = F.rows();
	int colsF = F.cols();

	// 创建一个unsigned int数组来存储三角形索引数据
	//unsigned int* triangles = new unsigned int[rowsF * colsF];
	//std::unique_ptr<float[]> triangles = std::make_unique<float[]>(rowsF * colsF);

	//// 将三角形索引数据从MatrixXi复制到unsigned int数组
	//for (int i = 0; i < rowsF; ++i) {
	//	for (int j = 0; j < colsF; ++j) {
	//		triangles[i * colsF + j] = static_cast<unsigned int>(F(i, j));
	//	}
	//}


	auto endhandleData = std::chrono::high_resolution_clock::now();
	auto elapsed1 = std::chrono::duration_cast<std::chrono::milliseconds>(endhandleData - starthandleData);
	std::cout << "Handle Data Time: " << elapsed1.count() << " ms.\n";


	auto startBuildBVH = std::chrono::high_resolution_clock::now();

	// create a new embree device
	RTCDevice device = rtcNewDevice(NULL);
	// create a new embree scene
	RTCScene scene = rtcNewScene(device);
	// create a new embree geometry for the triangle mesh
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

	float* vertices = (float*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(float) * 3, rowsV);

	for (int i = 0; i < rowsV; ++i) {
		for (int j = 0; j < colsV; ++j) {
			vertices[i * colsV + j] = static_cast<float>(V(i, j));
		}
	}

	unsigned int* indices = (unsigned int*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(unsigned int) * 3, rowsF);
	for (int i = 0; i < rowsF; ++i) {
		for (int j = 0; j < colsF; ++j) {
			indices[i * colsF + j] = static_cast<unsigned int>(F(i, j));
		}
	}

	// set the vertex buffer of the geometry
	//rtcSetSharedGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4, vertices.get(), 0, sizeof(float) * 4, rowsV);
	// set the index buffer of the geometry
	//rtcSetSharedGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, triangles.get(), 0, sizeof(unsigned int) * 3, rowsF);
	// commit the geometry
	rtcCommitGeometry(geom);
	// attach the geometry to the scene
	unsigned int geomID = rtcAttachGeometry(scene, geom);


	std::cout << "geomID: " << geomID << std::endl;
	// commit the scene
	rtcCommitScene(scene);

	auto endBuildBVH = std::chrono::high_resolution_clock::now();
	auto elapsed2 = std::chrono::duration_cast<std::chrono::milliseconds>(endBuildBVH - startBuildBVH);
	std::cout << "Build BVH Time: " << elapsed2.count() << " ms.\n";



	/////////////////////////////////////// ray test ///////////////////////////////////////
	

	auto startRayIntersect = std::chrono::high_resolution_clock::now();

	std::vector<RTCRayHit> rayHits;

	// 创建一个 RTCRayHit 对象，表示一条射线和一个交点
	RTCRayHit rayhit;


	rayhit.ray.org_x = 0.5;
	rayhit.ray.org_y = 0.5;
	rayhit.ray.org_z = -30;
	rayhit.ray.dir_x = 0;
	rayhit.ray.dir_y = 0;
	rayhit.ray.dir_z = 1;
	rayhit.ray.tnear = 1E-7; // ray min distance
	rayhit.ray.tfar = 100; // ray max distance
	rayhit.ray.mask = 0xFFFFFFFF; // ray code
	rayhit.ray.flags = 0; // ray flat

	rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID; // geo ID
	rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID; // tri ID
	rayhit.hit.Ng_x = 0.0f; // intersecting point normal x
	rayhit.hit.Ng_y = 0.0f; // intersecting point normal y
	rayhit.hit.Ng_z = 0.0f; // intersecting point normal z


	rtcIntersect1(scene, &rayhit);


	// Plot the stl mesh & ray
	igl::opengl::glfw::Viewer viewer;


	while (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID) { // ray next
		int geoId = rayhit.hit.geomID;
		// calc intersecting point
		double hitX = rayhit.ray.org_x + rayhit.ray.tfar * rayhit.ray.dir_x;
		double hitY = rayhit.ray.org_y + rayhit.ray.tfar * rayhit.ray.dir_y;
		double hitZ = rayhit.ray.org_z + rayhit.ray.tfar * rayhit.ray.dir_z;

			Eigen::MatrixXd hitP(1, 3);
		hitP << hitX, hitY, hitZ;

		Eigen::MatrixXd hitC(1, 3);
		hitC << .0, 0.0, 1.0;

		viewer.data().add_points(hitP, hitC);



		rayhit.ray.org_x = hitX;
		rayhit.ray.org_y = hitY;
		rayhit.ray.org_z = hitZ + 1E-6;
		rayhit.ray.dir_x = 0;
		rayhit.ray.dir_y = 0;
		rayhit.ray.dir_z = 1;
		rayhit.ray.tnear = 1e-7; // ray min distance
		rayhit.ray.tfar = std::numeric_limits<float>::infinity(); // ray max distance
		rayhit.ray.mask = 0xFFFFFFFF; // ray code
		rayhit.ray.flags = 0; // ray flat

		rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID; // geo ID
		rayhit.hit.primID = RTC_INVALID_GEOMETRY_ID; // tri ID
		rayhit.hit.Ng_x = 0.0f; // intersecting point normal x
		rayhit.hit.Ng_y = 0.0f; // intersecting point normal y
		rayhit.hit.Ng_z = 0.0f; // intersecting point normal z

		rtcIntersect1(scene, &rayhit);

	}


	//double hitX, hitY, hitZ;

	//// 检查是否有交点
	//if (rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	//{
	//	// first point
	//	rayHits.push_back(rayhit);
	//	std::cout << "Ray hit geoId = " << rayhit.hit.geomID << std::endl;

	//	// 如果有相交点，打印出交点坐标
	//	hitX = rayhit.ray.org_x + rayhit.ray.tfar * rayhit.ray.dir_x;
	//	hitY = rayhit.ray.org_y + rayhit.ray.tfar * rayhit.ray.dir_y;
	//	hitZ = rayhit.ray.org_z + rayhit.ray.tfar * rayhit.ray.dir_z;

	//	std::cout << "Ray hit tirP = " << hitX << " , " << hitY << " , " << hitZ << std::endl;


	//	Eigen::MatrixXd hitP(1, 3);
	//	hitP << hitX, hitY, hitZ;

	//	Eigen::MatrixXd hitC(1, 3);
	//	hitC << .0, 0.0, 1.0;

	//	viewer.data().add_points(hitP, hitC);

	//	// second point 
	//	rayhit.ray.org_x = hitX;
	//	rayhit.ray.org_y = hitY;
	//	rayhit.ray.org_z = hitZ;
	//	rayhit.ray.tfar = 100;


	//	rtcIntersect1(scene, &rayhit);


	//	hitX = rayhit.ray.org_x + rayhit.ray.tfar * rayhit.ray.dir_x;
	//	hitY = rayhit.ray.org_y + rayhit.ray.tfar * rayhit.ray.dir_y;
	//	hitZ = rayhit.ray.org_z + rayhit.ray.tfar * rayhit.ray.dir_z;

	//	std::cout << "Ray hit tirP2 = " << hitX << " , " << hitY << " , " << hitZ << std::endl;

	//	Eigen::MatrixXd hitP2(1, 3);
	//	hitP2 << hitX, hitY, hitZ;

	//	Eigen::MatrixXd hitC2(1, 3);
	//	hitC2 << .0, 0.0, 1.0;

	//	viewer.data().add_points(hitP2, hitC2);


	//}
	//else {
	//	std::cout << "Ray hit nothing. " << std::endl;
	//}


	/////////////////////////////////////// ray test ///////////////////////////////////////

	
	
	// mesh
	viewer.data().set_mesh(V, F);

	// ray
	Eigen::MatrixXd P1(1, 3);
	P1 << rayhit.ray.org_x, rayhit.ray.org_y, rayhit.ray.org_z;

	Eigen::MatrixXd P2(1, 3);
	P2 << rayhit.ray.org_x + rayhit.ray.tfar * rayhit.ray.dir_x, rayhit.ray.org_y + rayhit.ray.tfar * rayhit.ray.dir_y, rayhit.ray.org_z + rayhit.ray.tfar * rayhit.ray.dir_z;

	Eigen::MatrixXd C(1, 3);
	C << 1.0, 0.0, 0.0;

	viewer.data().add_points(P1, C);
	viewer.data().add_edges(P1, P2, C);

	// hit point
	

	
	viewer.launch();


	// release the embree objects
	rtcReleaseGeometry(geom);
	rtcReleaseScene(scene);
	rtcReleaseDevice(device);


}

void fake_uniform_hex() {
	
	// read F-35
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd N;
	// Load a mesh in OFF format
	std::string fileName = TEST_DATA_PATH"/F35.stl";
	std::ifstream input(fileName);

	if (input.is_open()) {
		igl::readSTL(input, V, F, N);
	}
	else {
		std::cout << "Fail to open file " << fileName << std::endl;
		return;
	}

	// 1.1 generate geo
	int rowsV = V.rows();
	int colsV = V.cols();

	// 动态分配内存
	//float* vertices = new float[rowsV * colsV];

	std::unique_ptr<float[]> vertices = std::make_unique<float[]>(rowsV * colsV);

	// 将顶点数据从MatrixXd复制到float数组
	for (int i = 0; i < rowsV; ++i) {
		for (int j = 0; j < colsV; ++j) {
			vertices[i * colsV + j] = static_cast<float>(V(i, j));
		}
	}

	// 获取矩阵的行数和列数
	int rowsF = F.rows();
	int colsF = F.cols();

	// 创建一个unsigned int数组来存储三角形索引数据
	unsigned int* triangles = new unsigned int[rowsF * colsF];

	// 将三角形索引数据从MatrixXi复制到unsigned int数组
	for (int i = 0; i < rowsF; ++i) {
		for (int j = 0; j < colsF; ++j) {
			triangles[i * colsF + j] = static_cast<unsigned int>(F(i, j));
		}
	}

	// create a new embree device
	RTCDevice device = rtcNewDevice(NULL);
	// create a new embree scene
	RTCScene scene = rtcNewScene(device);
	// create a new embree geometry for the triangle mesh
	RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
	// set the vertex buffer of the geometry
	rtcSetSharedGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4, vertices.get(), 0, sizeof(float) * 4, 4);
	// set the index buffer of the geometry
	rtcSetSharedGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, triangles, 0, sizeof(unsigned int) * 3, 2);
	// commit the geometry
	rtcCommitGeometry(geom);
	// attach the geometry to the scene
	unsigned int geomID = rtcAttachGeometry(scene, geom);
	// commit the scene
	rtcCommitScene(scene);

	// 1.0 get scene bbox (air_box)
	RTCBounds AirBox;
	rtcGetSceneBounds(scene, &AirBox);
	float x_center = (AirBox.lower_x + AirBox.upper_x) / 2.0f;
	float y_center = (AirBox.lower_y + AirBox.upper_y) / 2.0f;
	float z_center = (AirBox.lower_z + AirBox.upper_z) / 2.0f;
	AirBox.lower_x = 1.2 * (AirBox.lower_x - x_center) + x_center;
	AirBox.lower_y = 1.2 * (AirBox.lower_y - y_center) + y_center;
	AirBox.lower_z = 1.2 * (AirBox.lower_z - z_center) + z_center;
	AirBox.upper_x = 1.2 * (AirBox.upper_x - x_center) + x_center;
	AirBox.upper_y = 1.2 * (AirBox.upper_y - y_center) + y_center;
	AirBox.upper_z = 1.2 * (AirBox.upper_z - z_center) + z_center;


	// 1.1 get geo bbox
	RTCBoundsFunctionArguments GeoBoxArgs;
	RTCBounds GeoBox;
	GeoBoxArgs.geometryUserPtr = geom;
	GeoBoxArgs.bounds_o = &GeoBox;
	GeoBoxArgs.primID = geomID;

	
	rtcSetGeometryBoundsFunction(geom, (RTCBoundsFunction)(&GeoBoxArgs), NULL);
	
	const double D_MAX = 0.3;
	std::vector<double> XOY_LINES, YOZ_LINES, ZOX_LINES;
	double X0 = AirBox.lower_x;
	double Y0 = AirBox.lower_y;
	double Z0 = AirBox.lower_z;

	double XL = AirBox.upper_x - AirBox.lower_x;
	double YL = AirBox.upper_y - AirBox.lower_y;
	double ZL = AirBox.upper_z - AirBox.lower_z;
	double tol = 10e-8;

	for (int L = 0; L < XL; L += D_MAX) {
		XOY_LINES.push_back(X0 + L);
	}

	for (int L = 0; L < YL; L += D_MAX) {
		YOZ_LINES.push_back(Y0 + L);
	}

	for (int L = 0; L < ZL; L += D_MAX) {
		ZOX_LINES.push_back(Y0 + L);
	}
	









	// Plot the mesh
	/*igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V, F);
	viewer.launch();*/

}


void occ_ray() {


	// 创建一个盒子形状
	BRepPrimAPI_MakeBox box(10.0, 10.0, 10.0);
	TopoDS_Shape boxShape = box.Shape();

	// 创建一条射线
	Handle(Geom_Curve) ray = new Geom_Line(gp_Pnt(5.0, 5.0, 5.0), gp_Dir(.0, .0, 9.0));
	// 创建一个曲线和曲面的交点计算器
	GeomAPI_IntCS intersector;

	// 假设已经定义了一个 TopoDS_Shell 对象 shell
	TopExp_Explorer exp(boxShape, TopAbs_FACE); // 创建一个拓扑遍历器，指定遍历的类型为面
	while (exp.More()) // 遍历 shell 中的所有面
	{
		const TopoDS_Face& face = TopoDS::Face(exp.Current()); // 获取当前的面
		const Handle(Geom_Surface)& surface = BRep_Tool::Surface(face); // 将面转为几何曲面
		intersector.Perform(ray, surface); // 指定曲线和曲面，进行求交计算



		// 遍历交点
		if (intersector.IsDone() && intersector.NbPoints() > 0) // 如果计算成功并且有交点
		{
			for (int i = 1; i <= intersector.NbPoints(); i++) // 遍历每个交点
			{
				gp_Pnt pnt = intersector.Point(i); // 获取交点的坐标
				Standard_Real u, v, w;
				intersector.Parameters(i, u, v, w); // 获取交点在曲线和曲面上的参数
				std::cout << "Intersection point " << i << ": " << pnt.X() << ", " << pnt.Y() << ", " << pnt.Z() << std::endl; // 打印交点的信息
			}
		}
		else // 如果计算失败或者没有交点
		{
			std::cout << "No intersection found." << std::endl; // 打印提示信息
		}


		// 对几何曲面进行一些操作
		exp.Next(); // 移动到下一个面
	}
}