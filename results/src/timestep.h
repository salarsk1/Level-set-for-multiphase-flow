//Written by Salar Safarkhani

#pragma once
#include "datag.h"
#include "sg.h"
#include "bound.h"
#include "advection.h"
#include "reinit.h"
#include "parallel.h"
#include "toolbox.h"
#include "monitor.h"
#include "param.h"
class timestep{
public:
	void lit_timestep();
	void lit_stable_dt(double& dt, int& n_subcycle);
	void timestep_init();
private:
	double litCFL;
	int n_wo_reinit, n_wo_reband; // do not need to be defined as static
};
