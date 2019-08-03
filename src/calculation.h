#pragma once
#include<vector>
#include<cmath>
#include<algorithm>
using std::vector;
class calculation{
public:

	void calcTimeStep(const vector<vector<double>>&,const vector<vector<double>>&, \
							const int, const int, const std::vector<std::vector<double>> &,\
							const std::vector<std::vector<double>> &, double &);

	void predict(const vector<vector<double>>&,const vector<vector<double>> &,\
					 const vector<vector<double>> &, \
					 const double, const int, const int,\
					 std::vector<std::vector<double>> &,\
					 std::vector<std::vector<double>>&,\
					 std::vector<std::vector<double>>&,\
					 std::vector<std::vector<double>>&, \
					 std::vector<std::vector<double>>&, \
					 std::vector<std::vector<double>>&, \
					 const std::vector<std::vector<double>> &, 
					 const std::vector<std::vector<double>> &);

	void poisson(const double, const int, const int, double dt, \
					 const std::vector<std::vector<double>>&,\
					 std::vector<std::vector<double>>&,\
					 std::vector<std::vector<double>>&, \
					 std::vector<std::vector<double>>&);

	void correct(const vector<vector<double>>&,const double dt, const int, const int,\
					 std::vector<std::vector<double>>&, \
					 std::vector<std::vector<double>>&,\
				  	 std::vector<std::vector<double>>&, \
					 std::vector<std::vector<double>>&,\
				 	 std::vector<std::vector<double>>&);

	void gs2d(const int,const int,const std::vector<std::vector<double>>&,\
				 std::vector<std::vector<double>>&, const int);

private:

	double pi = 4.0*atan(1.0);

};
