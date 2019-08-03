#pragma once
class weno{
public:
	void prepare_weno_ghost();
	double phi_weno_5th(const double, const double, const double, const double);



};
