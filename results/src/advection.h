//Written by Salar Safarkhani

#pragma once
#include<memory>
#include"redist.h"
#include"datag.h"
#include"parallel.h"
#include"bl.h"
#include"bound.h"
#include"timing.h"
#include"litBuffer.h"
#include"reinit.h"
#include"param.h"
class advection {
public:
	void lit_advection(const double);
	void filterVelocity();
	void advection_m_init();
private:
	const int WENO_advection = 5;
	const int RK_advection = 3;
	const double r23  = 2.0/3.0;
	const double r112 = 1.0/12.0;
	static bool velocityFiltered;
};

