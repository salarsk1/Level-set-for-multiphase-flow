//Written by Salar Safarkhani

#pragma once
#include"datag.h"
#include<array>
#include<memory>
#include<vector>
#include "matrix.h"
using std::array;
using std::vector;
class triCubic{
public:
	void triCubicInit();
	void triCubicGetCoeff(const vector<double> &f,const vector<double> &d3fdx,\
								 const vector<vector<double>> &dfdx,\
								 const vector<vector<double>> &d2fdx,\
								 vector<double> &coeff);
	void triCubicGetCoeff2(const vector<double> &xx, vector<double> &coeff);
	double triCubicEval(const vector<double> &coeff,const vector<double> &xyz);
	double triCubicDerEval(const vector<double> &coeff,\
								  const vector<double> &xyz,\
								  const vector<int> &ddir);
private:
	static array<array<double,64>, 64> A;
	static bool arraySet;
};

