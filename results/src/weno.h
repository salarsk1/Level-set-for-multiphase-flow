//Written by Salar Safarkhani

#pragma once
#include"datag.h"
#include"parallel.h"
#include"litBuffer.h"
#include"bound.h"
#include"param.h"
class weno{
	public:
//	void prepare_WENO_5th(double**,double**,double**,double**,bool TbandOnly=false); not used anywhere
	void prepare_WENO_ghost(size_t, size_t, double**,size_t, size_t, double**);
	double G_WENO_5th(const double, const double, const double, const double);
	inline double wWENO3(const double a, const double b){
		double r;
		r = (WENOepsilon+b*b)/(WENOepsilon+a*a);
		return 0.5/(1.0+2.0*r*r);
	};
//	void prepare_WENO_3rd(double**, double**); This function is not used anywhere
	private:
	const double WENOepsilon = 1.0e-6;
	const double r13 = 1.0/3.0;
	const double r16 = 1.0/6.0;
};

