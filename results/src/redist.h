//Written by Salar Safarkhani

#pragma once
#include"datag.h"
#include"parallel.h"
#include"bound.h"
#include"weno.h"
#include"toolbox.h"
#include"litBuffer.h"
#include <cmath>
class redist{
public:
	void lit_redist_pde(size_t, size_t, double**); 
	inline void redist_m_init(){
		std::unique_ptr<global_variable> pgv(new global_variable);
		max_iter_redist = ceil(pgv->nVelFilter * 2.0);
	}
private:
	const int WENO_redist = 5;
	const int RK_redist = 1;
	static int max_iter_redist;
	const double r23 = 2.0/3.0;
	const double r112 = 1.0/12.0;
};

