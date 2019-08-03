#pragma once
#include<vector>
using std::vector;
class utility{
public:
	double heaviside(const double phi, const double ds);
	void curvature(double dx, double dy, vector<vector<double>> &, \
						vector<vector<double>> &);
	void alpha_cal(vector<vector<double>>&, const vector<vector<double>>&,\
						const int M, const int N, const double ds);

	inline double dHeav_dwid(const double phi){
		const double pi = 4.0 * atan(1.0);
		if(phi <= width && phi >= -width){
			return 0.5*(1.0 + cos(pi * phi / width))/width;
		}
		else return 0.0;
	}
private:
	const double width = 0.008;
};
