//Written by Salar Safarkhani

#pragma once
#include "datag.h"
#include"monitor.h"
#include "parallel.h"
#include"param.h"
class lit_timer_t;
class timing{
public:
	double get_time(bool no_block);
	double get_time();
	void lit_timing_find_timer(const std::string &name, int &itimer);
	void lit_timing_pre_init();
	void lit_timing_post_init();
	void lit_timing_create(const std::string &name);
	void lit_timing_start(const std::string &name);
	void lit_timing_start(const std::string &name, bool no_block);
	void lit_timing_stop(const std::string &name, bool no_block);
	void lit_timing_stop(const std::string &name);
	void lit_timing_start_full();
	void lit_timing_monitor();
private:
	const int lit_ntimers_max = 32;
	static int lit_ntimers;
	static int lit_nstarted;
	static std::unique_ptr<lit_timer_t []> lit_timers;
	static std::unique_ptr<int []> lit_started; //pointer to array; checked
	static double lit_full_time_in;
	static double  *lit_mval; //pointer to array ; checked
};

class lit_timer_t{
public:
	std::string name;
	double time;
	double time_in;
	bool started;
};

