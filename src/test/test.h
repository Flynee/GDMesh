#ifndef __GD__TESTH__
#define __GD__TESTH__

#include <vector>

void test_start();


std::vector<double> small_to_big(double L, double dmax, double max_step, double min_step, double R);
std::vector<double> avg_to_avg(double L, double step);
std::vector<double> dA_to_dM_to_dB(double L, double dmax, double dA, double dB, double R);

// draw mesh
void draw_mesh();

// draw box
void draw_box();

// read stl
void read_stl();

// decimate mesh
void mesh_decimate();

// read step
void read_step();

// write to step
void write_step();

// ray insection & mesh
void ray_mesh();

void ray_mesh2();

// non-uniform hex
void fake_uniform_hex();

// occ ray
void occ_ray();

void occ_ray2();

void occ_ray3();


// dA = dB = dmax
void dA_dB_dmax();

// dA increase to dM, dM > 0 && dM = dmax, dM descrease to dB
void dA_dM_dB();


// sort pnts by Z
void sort_pnt_byZ();

#endif