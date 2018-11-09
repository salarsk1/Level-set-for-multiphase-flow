#pragma once
#include<string>
#include<vector>
using std::vector;
class advection{
public:
	advection(std::string s);
	void advect_sol(const double, const double, const double, \
						 const vector<vector<double>>&,const vector<vector<double>>&,\
						 vector<vector<double>> &);

private:
	std::string advection_method;
};
