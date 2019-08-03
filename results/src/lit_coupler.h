//Written by Salar Safarkhani

#pragma once
#include "datag.h"
#include "timestep.h"
#include "io.h"
#include "parallel.h"
#include "toolbox.h"
#include "advection.h"
#include "monitor.h"
#include "timing.h"
#include "band.h"
#include "bound.h"
#include "param.h"
#include "redist.h"
#include "init.h"
#include "timestep.h"
#include "sg.h"
#include "litBuffer.h"
#include "reinit.h"
#include<boost/filesystem.hpp>
class lit_coupler{
public:
	void litRunIteration(const double, const bool, int arg=1);
	void lit_initialize(int argc, char *argv[]);
	void lit_initialize2(int argc, char*argv[], std::string);
	void lit_finalize();
};
