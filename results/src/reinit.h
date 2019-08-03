//Written by Salar Safarkhani

#pragma once
#include"datag.h"
#include"parallel.h"
#include"bl.h"
#include"bound.h"
#include"weno.h"
#include"toolbox.h"
#include"litBuffer.h"
#include"param.h"
#include"monitor.h"
#include"timing.h"
#include"param.h"
class reinit{
public:
	void lit_reinit_pde();
	bool trigger_reinit(const int);
	void reinit_m_init();
	void lit_clip_G();
	void calcSignFunction(double*);
	void calcHOCR2rijk(double *);
	void calcHOCR2fijk(double*, double*);
	static double triggerGradMin, triggerGradMax;
	static int min_iter_reinit;
private:
	static int RK_reinit;
	static int max_iter_reinit;
	//static int min_iter_reinit;
	const double r23 = 2.0/3.0;
	const double r112 = 1.0/12.0;
	const double hocr2w = 0.5;
	//static double triggerGradMin, triggerGradMax;
	static vector<vector<double>> alphaRK;
};

