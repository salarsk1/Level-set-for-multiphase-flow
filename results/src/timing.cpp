//Written by Salar Safarkhani

#include "timing.h"
#include<mpi.h>
#include<iostream>
#include "monitor.h"
#include"float.h"
double timing::get_time(bool no_block){
return MPI_Wtime();
};
double timing::get_time(){
MPI_Barrier(MPI_COMM_WORLD);
return MPI_Wtime();
};
void timing::lit_timing_find_timer(const std::string &name, int &itimer){
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<parallel> pp(new parallel);
for(itimer=0; itimer<lit_ntimers; itimer++){
	if(ppar->trim(lit_timers[itimer].name) == ppar->trim(name)) break;
};
if(itimer+1 > lit_ntimers){
	std::cout<<"Timer : " << ppar->trim(name)<<std::endl;
	pp->litError("timing_find_timer","unknown timer");
};
};

void timing::lit_timing_pre_init(){
std::unique_ptr<parallel> pp(new parallel);
lit_ntimers  = 0;
lit_nstarted = 0;
lit_full_time_in = get_time();
if(pp->ierr == 0){
	lit_timers.reset(new lit_timer_t[lit_ntimers_max]);
	lit_started.reset(new int[lit_ntimers_max]);
}
else{
	pp->litError("lit_timing_pre_init","allocation error for lit_timers,lit_started");
};
};

void timing::lit_timing_post_init(){
double tick;
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<litparam> ppar(new litparam);
if(pgv->verbose && pp->myrank == 0){
	tick = MPI_Wtick();
	std::cout<<"Timing resolution is "<<tick<<" s"<<std::endl;
};
if(pp->ierr==0){
	lit_mval = new double[2*lit_ntimers+4];
}
else{
	pp->litError("lit_timing_post_init", "allocation error for for lit_mval");
};
pmo->lit_monitor_create_file_step("lit_timing", 2*lit_ntimers+4);
pmo->lit_monitor_set_header(1, "Total [s]" , 'r');
pmo->lit_monitor_set_header(2, "TIme/points [us]", 'r');
for(auto i=0; i<lit_ntimers; i++){
	pmo->lit_monitor_set_header((i+1)*2+1, ppar->trim(lit_timers[i].name)+"[s]", 'r'); // checked
	pmo->lit_monitor_set_header((i+1)*2+2, ppar->trim(lit_timers[i].name)+"[%]", 'r'); // checked
};
pmo->lit_monitor_set_header(2*lit_ntimers+3, "Rest [s]", 'r');
pmo->lit_monitor_set_header(2*lit_ntimers+4, "Rest [%]", 'r');
};

void timing::lit_timing_create(const std::string &name){
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<parallel> pp(new parallel);
lit_ntimers = lit_ntimers + 1;
if(lit_ntimers > lit_ntimers_max){
	pp->litError("timing_create", "too many timers");
};
lit_timers[lit_ntimers-1].name = ppar->trim(name);
lit_timers[lit_ntimers-1].time = 0.0;
lit_timers[lit_ntimers-1].time_in = 0.0;
lit_timers[lit_ntimers-1].started = false;
};
void timing::lit_timing_start(const std::string &name){
double time_in;
int itimer;
std::unique_ptr<parallel> pp(new parallel);
lit_timing_find_timer(name, itimer);
if(lit_timers[itimer].started){
	std::cout<<"Timer : "<<name<<std::endl;
	pp->litError("timing_start","timer already started");
};
time_in = get_time();
lit_timers[itimer].time_in = time_in;
lit_timers[itimer].started = true;
if(lit_nstarted != 0){
	auto d1 = lit_started[lit_nstarted-1]-1;
   lit_timers[d1].time=lit_timers[d1].time+time_in-lit_timers[d1].time_in;
	lit_timers[d1].time_in = 0.0;
};
lit_nstarted = lit_nstarted + 1;
lit_started[lit_nstarted-1] = itimer;
};

void timing::lit_timing_start(const std::string &name, bool no_block){
double time_in;
int itimer = 0;
std::unique_ptr<parallel> pp(new parallel);
lit_timing_find_timer(name, itimer);
if(lit_timers[itimer].started){
	std::cout<<"Timer : "<<name<<std::endl;
	pp->litError("timing_start","timer already started");
};
time_in = get_time(no_block);
lit_timers[itimer].time_in = time_in;
lit_timers[itimer].started = true;
if(lit_nstarted != 0){
	auto d1 = lit_started[lit_nstarted-1];
   lit_timers[d1].time=lit_timers[d1].time+time_in-lit_timers[d1].time_in;
	lit_timers[d1].time_in = 0.0;
};
lit_nstarted = lit_nstarted + 1;
lit_started[lit_nstarted-1] = itimer;
};

void timing::lit_timing_stop(const std::string &name){
int itimer = 0;
double time_out;
std::unique_ptr<parallel> pp(new parallel);
lit_timing_find_timer(name, itimer);
if(!lit_timers[itimer].started){
	std::cout<<"Timer : "<<name<<std::endl;
	pp->litError("timing_stop","timer already stopped");
};
time_out = get_time();
auto d1 = itimer;
lit_timers[d1].time=lit_timers[d1].time+time_out-lit_timers[d1].time_in;
lit_timers[d1].time_in = 0.0;
lit_timers[d1].started = false;
lit_nstarted = lit_nstarted - 1;
if(lit_nstarted != 0){
	auto d1 = lit_started[lit_nstarted-1]-1;
	lit_timers[d1].time_in = time_out;
};
};
void timing::lit_timing_stop(const std::string &name, bool no_block){
int itimer = 0;
double time_out;
std::unique_ptr<parallel> pp(new parallel);
lit_timing_find_timer(name, itimer);
if(!lit_timers[itimer].started){
	std::cout<<"Timer : "<<name<<std::endl;
	pp->litError("timing_stop","timer already stopped");
};
time_out = get_time(no_block);
auto d1 = itimer;
lit_timers[d1].time=lit_timers[d1].time+time_out-lit_timers[d1].time_in;
lit_timers[d1].time_in = 0.0;
lit_timers[d1].started = false;
lit_nstarted = lit_nstarted - 1;
if(lit_nstarted != 0){
	auto d1 = lit_started[lit_nstarted-1]-1;
	lit_timers[d1].time_in = time_out;
};
};

void timing::lit_timing_start_full(){
lit_full_time_in = get_time();
};

void timing::lit_timing_monitor(){
double full_time_out, rest_time;
int n;
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<monitor> pmo(new monitor);
full_time_out = get_time();
rest_time = 0.0;
lit_mval[0] = full_time_out-lit_full_time_in+DBL_EPSILON;
n = pp->nodesInBandAll('N');
lit_mval[1] = 1.0e6 * lit_mval[0] / static_cast<double>(n);
for(auto i=0; i<lit_ntimers; i++){
	rest_time = rest_time+lit_timers[i].time;
	lit_mval[2*(i+1)] = lit_timers[i].time; // checked
	lit_mval[2*(i+1)+1] = 100.0 * lit_timers[i].time/lit_mval[0];
};
lit_mval[2*lit_ntimers+2] = lit_mval[0] - rest_time; // checked
lit_mval[2*lit_ntimers+3] = 100.0 *(lit_mval[0]-rest_time)/lit_mval[0]; // checked
pmo->lit_monitor_select_file("lit_timing");
pmo->lit_monitor_set_array_values(lit_mval); // checked
lit_full_time_in = full_time_out;
for(auto i=0; i<lit_ntimers; i++){
	lit_timers[i].time = 0.0;
	lit_timers[i].time_in = 0.0;

};
};

int timing::lit_ntimers = 0;
int timing::lit_nstarted = 0;
std::unique_ptr<lit_timer_t []> timing::lit_timers(nullptr);
std::unique_ptr<int []> timing::lit_started(nullptr);
double timing::lit_full_time_in = 0.0;
double* timing::lit_mval = nullptr;

