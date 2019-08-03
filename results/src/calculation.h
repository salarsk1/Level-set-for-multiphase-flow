#pragma once
#include<vector>
#include<cmath>
#include<algorithm>
class calculation{
public:

	void calcTimeStep(const double, const int, const int,const std::vector<std::vector<double>> &,\
							const std::vector<std::vector<double>> &, double &);

	void predict(const double, const double, const int, const int,\
					 std::vector<std::vector<double>> &,std::vector<std::vector<double>>&,\
					 std::vector<std::vector<double>>&,std::vector<std::vector<double>>&);

	void poisson(const int, const int, double dt, const std::vector<std::vector<double>>&,\
					std::vector<std::vector<double>>&, int &, const double);

	void correct(const double dt, const int, const int,\
					 std::vector<std::vector<double>>&,std::vector<std::vector<double>>&,\
				  	 std::vector<std::vector<double>>&,std::vector<std::vector<double>>&,\
				 	 std::vector<std::vector<double>>&);

	void gs2d(const int,const int,const std::vector<std::vector<double>>&,\
				 std::vector<std::vector<double>>&, const int);

private:

	double pi = 4.0*atan(1.0);

};
