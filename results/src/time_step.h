#pragma once
#include<vector>
#include<algorithm>
using std::vector;
class time_step{
public:
	void simulate(double, const double, const double, double &,
					  vector<vector<double>>&, vector<vector<double>> &,\
					  vector<vector<double>> &);

	void stable_time(vector<vector<double>>&,vector<vector<double>>&, \
							 double, double, const double);
private:
	double dt_sub;
	int n_sub;
	const double CFL_adv = 0.25;
	double lev_CFL;
};
