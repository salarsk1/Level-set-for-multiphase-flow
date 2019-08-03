//Written by Salar Safarkhani

#pragma once
#include<iostream>
#include "datag.h"
#include "parallel.h"
#include<fstream>
class init{
public:

	void init_m_init();
	void init_geometry();
	double G_init_value(const vector<double> &xyz);
	double G_inject(const vector<double> &, shape_t *shape);
	double G_notched_circle(const vector<double> &xyz);
	double G_circle(const vector<double> &xyz);
	double G_circle_xz(const vector<double> &xyz);
	double G_column(const vector<double> &xyz);
	double G_column_cap_x(const vector<double> &xyz);
	double G_column_cap_y(const vector<double> &xyz);
	double G_column_cap_z(const vector<double> &xyz);
	double G_column_pool(const vector<double> &xyz);
	double G_sphere_pool(const vector<double> &xyz);
	double G_randomCircle(const vector<double> &xyz);
	double G_randomJet(const vector<double> &xyz);
	double G_randomDrops(const vector<double> &xyz);
	double G_ring(const vector<double> &xyz);
	double G_rod(const vector<double> &xyz);
	double G_from_sphere(const vector<double>&xyz,const std::array<double,10>&data);
	double G_sphere(const vector<double> &xyz);
	double G_sphere2(const vector<double> &xyz);
	double G_sphere_cyl(const vector<double> &xyz);
	double G_disc_cyl(const vector<double> &xyz);
	double G_plane(const vector<double> &xyz);
	double G_sine(const vector<double> &xyz);
	double G_siney(const vector<double> &xyz);
	double G_sine_column(const vector<double> &xyz);
	double G_siney_column(const vector<double> &xyz);
	double G_cosine3D(const vector<double> &xyz);
	double G_RT3D(const vector<double> &xyz);
	double G_cosine(const vector<double> &xyz);
	double G_cosiney(const vector<double> &xyz);
	double G_sheet(const vector<double> &xyz);
	double G_randomSheet(const vector<double> &xyz);
	double G_deformed_column(const vector<double> &xyz);
	double G_ellipse(const vector<double> &xyz);
	double G_deformed_sphere(const vector<double> &xyz);
	double G_milk_crown(const vector<double> &xyz);
	double G_bursting_bubble(const vector<double> &xyz);
	double G_bursting_bubble_3D(const vector<double> &xyz);
	double G_bursting_bubble_2D(const vector<double> &xyz);
	double G_bursting_bubble_rim_3D(const vector<double> &xyz);
	double G_rayleigh_axi(const vector<double> &xyz);
	double G_rayleigh(const vector<double> &xyz);
	double G_bursting_bubble_rim_2D(const vector<double> &xyz);
	double G_ligament_2D(const vector<double> &xyz);
	double G_droplens(const vector<double> &xyz);
	double G_Dam_break_2D(const vector<double> &xyz);
	double G_function2D(const vector<double> &xyz);
	double arctan(double ddx, double ddy);
private:
	double **xd;
	std::unique_ptr<double[]> dd;
};

