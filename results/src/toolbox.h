//Written by Salar Safarkhani

#pragma once
#include<mpi.h>
#include"datag.h"
#include"parallel.h"
#include"bl.h"
#include"bound.h"
#include"tricubic.h"
#include"litBuffer.h"
#include"timing.h"
#include"init.h"
#include"monitor.h"
class toolbox{
public:
	void toolbox_m_init();
	double delta(const double G);
	double delta_L(const double G);
	bool cell_contains_front(const Gnode_t &Gnh, const int ic, block_t *bb);
	double volume_delft(const double G, const std::vector<double> &gradG);
	double volume_delft_cyl(const double,const double,const double,\
									const std::vector<double> &gradG);
	double cell_volume_fraction(const Gnode_t&,const int,block_t*,const int indic=0);
	bool is_indicator(const Gnode_t&, const int ic, block_t*); // ic checked
	void findClosestFrontPointDirect(size_t, size_t, double **); //checked
	void curvature(size_t, double*);
	void solve_linear_system(std::vector<std::vector<double>>&,\
									 std::vector<double>&,std::vector<double>&, const int n);
	double heaviside_st(const double);
	void curvature_indicator(size_t, double*, size_t, double*); // it is related to getR2Buffer; //checked
	void direct_curvature_indicator(size_t,double *,size_t, double*); // checked
	double triLinear(const vector<double>&, const vector<double>&);
	void calcCurvature(size_t,double*,size_t,double*);
	void normal(size_t, size_t, double**); // it is related to getR2Buffer; // checked
	double volume();
	double volumeTband();
	void interfaceHolder();
	double calcAmplitude_quarter();
	void diag_RT3D(vector<double>&,vector<double>&,vector<double>&,\
						vector<double>&,vector<double>&);
	void GaussPivotLarge(vector<vector<double>>&,vector<double>&,const int,\
								vector<double>&);
	void diag_RMI(vector<double>&);
	void diag_ligament(vector<double>&);
	void diag_drop(vector<double>&);
	void diagnostics_init();
	void diagnostics();
	void checkBandLayerGmin();
	void rotating_cyl_mfix_init();
	void rotating_cyl_mfix_adjust(double inVol=-1.0);
private:
	double width_H;
	double r2pi, pi_wh, r2wh;
	static double width_d;
	static double pi_wd;
	static double r2wd;
	static double r13;
};
