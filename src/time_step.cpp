#include<memory>
#include<cmath>
#include<iostream>
#include"time_step.h"
#include"alg.h"
#include"reinitialization.h"
#include"advection.h"
void time_step::simulate(double dt, double dx, double dy, double &CFL, \
								 vector<vector<double>>& u, vector<vector<double>> &v, \
							    vector<vector<double>> & phi){
	std::unique_ptr<advection> padv(new advection("WENO-5"));
	std::unique_ptr<reinitialize> prein(new reinitialize);
	stable_time(u,v,dx,dy,CFL,dt);
	for(auto i=0; i<n_sub; i++){
		std::cout<<"subiteration  "<<i+1<<"  of  "<<n_sub<<std::endl;
		padv->advect_sol(dt_sub,dx,dy,u,v,phi);
		if(prein->trigger_reinitialize(dx,dy,0,phi)){
			prein->pde_reinit(dt_sub, dx, dy, phi);
		}
	}
}
void time_step::stable_time(vector<vector<double>>&u,vector<vector<double>>&v,\
									 double dx, double dy, double &CFL,\
									 const double dt){
	alg algor;
	auto size1_u = u.size();
	auto size1_v = v.size();
	auto size2_u = u[0].size();
	auto size2_v = v[0].size();
	double max_vel;
	max_vel = std::max(algor.maxval(size1_u,u), algor.maxval(size1_v, v));
	dt_sub  = CFL * std::min(dx,dy) / max_vel;
	n_sub = ceil(dt/dt_sub);
	dt_sub = dt / static_cast<double> (n_sub);
//	CFL = max_vel / (std::min(dx,dy)/dt_sub);
	std::cout<<"max_vel = "<<max_vel<<'\t'<<"dt_sub = "<<dt_sub<<'\t';
	std::cout<<"n_sub = "<<n_sub<<'\t'<<"CFL = "<<CFL<<std::endl;
}
