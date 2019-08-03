#pragma once
#include<string>
#include<vector>
class init{
public:
	init(std::string, const std::vector<double>&);
	void init_geom(std::string);
	double phi_init_value(std::string, const std::vector<double> &);
	double phi_circle(const std::vector<double> &);
	double phi_plane(const std::vector<double> &);
	double phi_zalesak(const std::vector<double> &);
private:
	std::string shape;
	std::vector<double> shape_info;
};
