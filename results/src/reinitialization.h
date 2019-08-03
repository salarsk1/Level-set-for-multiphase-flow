#pragma once
#include<vector>
#include<string>
using std::vector;
class reinitialize{
public:
	vector<vector<double>> sign_func(const double, const double, \
												const vector<vector<double>> &);
	void pde_reinit(double dt,const double dx,const double dy,\
						 vector<vector<double>> &);
	bool trigger_reinitialize(const double, const double, const int, \
									  vector<vector<double>>&);
private:
	double trigger_min = 0.0001;
	double trigger_max = 2.0;
	bool trigger;
	const std::string boundary = "periodic";

};
