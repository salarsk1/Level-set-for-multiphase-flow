#include"init.h"
#include<cmath>
#include<iostream>
init::init(std::string s, const std::vector<double>&xyz){
	shape = s;
	shape_info = xyz;
}
void init::init_geom(std::string s){
	shape = s;
}
double init::phi_init_value(std::string s, const std::vector<double>&xyz){
	init_geom(s);
	double GV;
	if(shape == "circle"){
		GV = phi_circle(xyz);
	}
	else if(shape == "plane"){
		GV = phi_plane(xyz);
	}
	else if(shape == "zalesak"){
		GV = phi_zalesak(xyz);
	}
	return GV;
}
double init::phi_circle(const std::vector<double> &xyz){
	double res;
	res = shape_info[2] - sqrt(pow(xyz[0]-shape_info[0],2)+pow(xyz[1]-shape_info[1],2));
	return res;
}

double init::phi_plane(const std::vector<double> &xyz){
	double res;
	res = shape_info[0] - xyz[0];
	return res;
}
double init::phi_zalesak(const std::vector<double> &xyz){
	double c, bb, b1, b2, h1, h2;
	c=shape_info[2]-sqrt(pow(xyz[0]-shape_info[0],2)+pow(xyz[1]-shape_info[1],2));
	b1 = shape_info[0] - 0.5*shape_info[3];
	b2 = shape_info[0] + 0.5*shape_info[3];
	h1 = shape_info[1] - shape_info[2]*cos(asin(0.5*shape_info[3]/shape_info[2]));
	h2 = shape_info[1] - shape_info[2]+shape_info[4];

	if(0.0 <= c && xyz[0] <= b1 && xyz[1] <= h2){
		bb = b1 - xyz[0];
		return std::min(c,bb);
	}
	else if(0.0 <= c && b2 <= xyz[0] && h2 >= xyz[1]){
		bb = xyz[0] - b2;
		return std::min(c,bb);
	}
	else if(0.0 <= c && b1 <= xyz[0] && xyz[0] <= b2 && h2 <= xyz[1]){
		bb = xyz[1] - h2;
		return std::min(c,bb);
	}
	else if(0.0 <= c && xyz[0] <= b1 && h2 <= xyz[1]){
		bb = sqrt(pow(xyz[0]-b1,2)+pow(xyz[1]-h2,2));
		return std::min(c,bb);
	}
	else if(0.0 <= c && b2 <= xyz[0] && h2 <= xyz[1]){
		bb = sqrt(pow(xyz[0]-b2,2)+pow(xyz[1]-h2,2));
		return std::min(c,bb);
	}
	else if(b1 <= xyz[0] && xyz[0] <= b2 && xyz[1] <= h2 && h1 <= xyz[1]){
		auto d1 = std::abs(xyz[0]-b1);
		auto d2 = std::abs(xyz[0]-b2);
		auto d3 = std::abs(xyz[1]-h2);
		return -std::min(std::min(d1,d2),d3);
	}
	else if(b1 <= xyz[0] && xyz[0]<=b2 && xyz[1] <= h1){
		auto d1 = sqrt(pow(xyz[0]-b1,2)+pow(xyz[1]-h1,2));
		auto d2 = sqrt(pow(xyz[0]-b2,2)+pow(xyz[1]-h1,2));
		return -std::min(d1,d2);
	}
	else{
		return c;
	}

}
