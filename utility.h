#pragma once
#include<vector>
using std::vector;
class utility{
public:
	double heaviside(const double phi);
	void curvature(double dx, double dy, vector<vector<double>> &, \
						vector<vector<double>> &);
	void tridiag(vector<double>&, vector<double>&, vector<double>&, vector<double>&,\
					 vector<double>&);
inline double deltamo(const double phi){
   const double width = 0.12;
   const double pi = 4.0*atan(1.0);
   if(std::abs(phi) < width){
      return 0.5 * (1.0 + cos(pi * phi/width)) / width;
   }
   else return 0.0;
}
};
