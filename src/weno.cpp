#include"weno.h"
#include<cmath>
double weno::phi_weno_5th(const double a, const double b, const double c, const double d){
	const double weno_eps = 1.0e-6;
	const double r13 = 1.0/3.0;
	const double r16 = 1.0/6.0;
	double w0, w2, a0, a1, a2, IS0, IS1, IS2, rh;
	IS0 = 13.0 * pow(a-b, 2) + 3.0 * pow(a-3.0*b, 2);
	IS1 = 13.0 * pow(b-c, 2) + 3.0 * pow(b+c, 2);
	IS2 = 13.0 * pow(c-d, 2) + 3.0 * pow(3.0*c-d, 2);
	a0 = pow(weno_eps + IS0, -2);
	a1 = 6.0 * pow(weno_eps + IS1, -2);
	a2 = 3.0 * pow(weno_eps + IS2, -2);
	rh = 1.0/(a0 + a1 + a2);
	w0 = a0 * rh;
	w2 = a2 * rh;
	return (r13 * w0 * (a-2.0*b+c)+r16*(w2-0.5)*(b-2.0*c+d));
}
