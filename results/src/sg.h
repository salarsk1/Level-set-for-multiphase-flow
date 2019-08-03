//Written by Salar Safarkhani

#pragma once
#include"datag.h"
#include"monitor.h"
#include"timing.h"
#include"parallel.h"
#include"bl.h"
#include"band.h"
#include"bound.h"
class sg{
public:
	void sg_m_init();
	void sg_init();
	void sg_band_update();
	void sg_calc_active();
	void sg_load_balance(bool force, const char);
	void sg_load_balance_metis(bool force);
	void sg_load_balance_parmetis(bool force);
	void sg_load_balance_parmetis_fix(bool force);
	void sg_load_balance_lit();
	void executeBlockExchange(const int nex,const vector<vector<int>> &exScript);
	std::vector<int> xyz_to_ijk_sg(const std::array<double, 3> &xyz);
	void sort_pick(vector<vector<int>> &arr);
	void monitor_load(double &imbalance);
};

