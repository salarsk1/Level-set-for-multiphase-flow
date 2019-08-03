#include<cmath>
#include"utility.h"
#include<iostream>
double utility::heaviside(const double phi){
	double H;
	const double width = 0.12;
	const double pi = 4.0 * atan(1.0);
	if(phi <= -width) H = 0.0;
	else if(phi >= width) H = 1.0;
	else H = 0.5 + 0.5 * phi / width + 0.5 * sin(pi * phi / width) / 	pi;
//	if(std::abs(phi)<1.0) H = phi;
//	else if(phi<=-1.0) H = -1.0;
//	else H = 1.0;
	return H;
}
void utility::curvature(double dx, double dy, vector<vector<double>> &kappa,\
							   vector<vector<double>> &phi){
	auto size1 = phi.size();
	auto size2 = phi[0].size();
	double phi_x, phi_y, phi_xx, phi_yy, phi_xy, temp;
	for(auto i=2; i<size1-2; i++){
		for(auto j=2; j<size2-2; j++){
			phi_x  = (phi[i+1][j] - phi[i-1][j]) / (2.0*dx);
			phi_xx = (phi[i+1][j] - 2.0*phi[i][j] + phi[i-1][j]) / (pow(dx,2));
			phi_y  = (phi[i][j+1] - phi[i][j-1]) / (2.0*dy);
			phi_yy = (phi[i][j+1] - 2.0*phi[i][j] + phi[i][j-1]) / (pow(dy,2));
			phi_xy = (phi[i+1][j+1] - phi[i-1][j+1])/(4.0*dx*dy) - \
						(phi[i+1][j-1] - phi[i-1][j-1])/(4.0*dx*dy);
			temp = pow(pow(phi_x,2) + pow(phi_y,2) + 10.0e-7, 1.5);
			kappa[i-2][j-2] = (phi_xx*pow(phi_y,2) + phi_yy*pow(phi_x,2) -
							   2.0*phi_x*phi_y*phi_xy)/temp;
		}
	}
}
void utility::tridiag(vector<double>&l, vector<double>&d, vector<double>&u,\
							 vector<double>&y, vector<double>&x){
	auto N = d.size();
	vector<double> d_bar(N, 0.0);
	for(auto j=0; j<N; j++){
		d_bar[j] = d[j];
	}
	//forward elimination
	for(auto j=1; j<N; j++){
		d_bar[j] += -l[j] / d_bar[j-1] * u[j-1];
		y[j] += -l[j] / d_bar[j-1] * y[j-1];
	}
	//backward substitution
	x[N-1] = y[N-1] / d_bar[N-1];
	for(int j=N-2; j>-1; j--){
		x[j] = (y[j] - u[j]*x[j+1]) / d_bar[j];
	}
}

