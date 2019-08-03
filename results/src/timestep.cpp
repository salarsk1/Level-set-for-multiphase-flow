//Written by Salar Safarkhani

#include<cmath>
#include "timestep.h"
void timestep::lit_timestep(){
double dt, imbalance;
int n_subcycles, dn_reinit, dn_reband, i_cycle;
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<advection> padv(new advection);
bool force;
if(pgv->verbose) pp->parallel_barrier();
if(pgv->verbose && pp->myrank == 0){
	std::cout<<pgv->clit<<"***Starting step "<<pgv->fs_step<<" time = "<<pgv->time;
	std::cout<<" dt ="<< pgv->fs_dt<<std::endl;
}
dt = 0.0;
n_subcycles = 0;
lit_stable_dt(dt, n_subcycles);
std::cout<<dt<<'\t'<<n_subcycles<<std::endl;
if(n_subcycles > pgv->maxn_wo_reinit){
	auto temp1 = ceil(static_cast<double>(n_subcycles)/static_cast<double>(pgv->maxn_wo_reinit));
	auto temp2 = ceil(static_cast<double>(n_subcycles)/static_cast<double>(temp1));
	dn_reinit = std::min(static_cast<double>(pgv->maxn_wo_reinit), static_cast<double>(temp2));
}
else dn_reinit = pgv->maxn_wo_reinit;
dn_reband = dn_reinit;
pmo->monitor_nSubCycles = n_subcycles;
for(auto i=0; i<n_subcycles; i++){
	ptool->interfaceHolder();
	padv->lit_advection(dt);
	ptool->interfaceHolder();
	psg->sg_band_update();

	n_wo_reinit += 1;
	n_wo_reband += 1;
	if(n_wo_reband >= dn_reband){
		psg->sg_band_update();
		n_wo_reband = 0;
	}
	if(prein->trigger_reinit(0)){
		if(n_wo_reband != 0){
			psg->sg_band_update();
				n_wo_reband = 0;
		}
		if(ppar->trim(pgv->reinit_solver) == "PDE")	prein->lit_reinit_pde();
		n_wo_reinit = 0;
	}
	else{
		pmo->lit_monitor_select_file("lit_reinit");
		pmo->lit_monitor_set_single_values(pmo->monitor_gradG_min, pmo->monitor_gradG_max);
		pmo->lit_monitor_dump_values_iter(2,0);
		prein->lit_clip_G();
	}
	pgv->time += dt;
}
if(n_wo_reband > 0){
	psg->sg_band_update();
	n_wo_reband = 0;
	if(prein->trigger_reinit(0)){
		if(ppar->trim(pgv->reinit_solver) == "PDE"){
			if(n_wo_reband != 0){
				pp->litError("timestep.cpp","Must update band before reinitializing");
			}
			prein->lit_reinit_pde();
		}
		n_wo_reinit = 0;
	}
	else{
		pmo->lit_monitor_select_file("lit_reinit");
		pmo->lit_monitor_set_single_values(pmo->monitor_gradG_min, pmo->monitor_gradG_max);
		pmo->lit_monitor_dump_values_iter(2,0);
		prein->lit_clip_G();
	}
}
psg->sg_load_balance(force, 'n');

//monitor
//monitor
ptool->diagnostics();
if(pgv->verbose) pp->parallel_barrier();
if(pgv->verbose && pp->myrank == 0){
	std::cout<<pgv->clit<<'\t'<<"*** Finished step "<<pgv->fs_step<<'\t';
	std::cout<<" time = "<<pgv->time<<std::endl;
}
};

void timestep::lit_stable_dt(double& dt, int& n_subcycles){
double my_maxvel, maxvel, my_dt, min_required_band_size,current_band_size;
double dxyz_l, maxvel_l;
double my_velRdxyz, velRdxyz;
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
my_maxvel = 1.0e-77;
my_dt = 1.0e77;
my_velRdxyz = 1.0e-77;
if(pgv->cylindrical && pgv->ijkm_gl[2]>1){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinN;ic++){
			auto j = pgv->Gn[ic].ijk[1]-1;
			double temp = std::max(fabs(pgv->Gn[ic].V[0]), fabs(pgv->Gn[ic].V[1]));
			maxvel_l = std::max(temp, static_cast<double>(fabs(pgv->Gn[ic].V[2])));
			double temp1 = std::min(pgv->dxyz[0], pgv->dxyz[1]);
			double temp2 = pgv->yc[j+pgv->nghost] * pgv->dxyz[2];
			dxyz_l = std::min(temp1, temp2);
			temp = my_dt;
			temp1 = pgv->CFL_advection*dxyz_l/maxvel_l;
			my_dt = std::min(temp, temp1);
			temp1 = my_maxvel;
			my_maxvel = std::max(temp1, maxvel_l);
			temp1 = my_velRdxyz;
			temp2 = maxvel_l/dxyz_l;
			my_velRdxyz = std::max(temp1, temp2);
		}
	}
	pp->parallel_all_min(my_dt, dt);
	pp->parallel_all_max(my_maxvel, maxvel);
	pp->parallel_all_max(my_velRdxyz, velRdxyz);
	n_subcycles = ceil(pgv->fs_dt/dt);
	dt = pgv->fs_dt/static_cast<double>(n_subcycles);
	litCFL = dt*velRdxyz;
}
else{
	double abs1, abs2, abs3;
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinN; ic++){
			abs1 = fabs(pgv->Gn[ic].V[0]);
			abs2 = fabs(pgv->Gn[ic].V[1]);
			abs3 = fabs(pgv->Gn[ic].V[2]);
			double temp1 = std::max(abs1, abs2);
			double temp2 = std::max(temp1, abs3);
			my_maxvel = std::max(my_maxvel,temp2);
		}
	}
	my_dt = pgv->CFL_advection * pgv->dxyz_min/my_maxvel;
	pp->parallel_all_min(my_dt, dt);
	pp->parallel_all_max(my_maxvel, maxvel);
	n_subcycles = ceil(pgv->fs_dt/dt);
	dt = pgv->fs_dt/static_cast<double>(n_subcycles);
	litCFL = maxvel / (pgv->dxyz_min/dt);
}
if(pgv->verbose && pp->myrank == 0){
	std::cout<<pgv->clit<<"max vel, sub-iteration dt, number of subsycles, litCFL"<<'\t';
	std::cout<<maxvel<<'\t'<<dt<<'\t'<<n_subcycles<<'\t'<<litCFL<<std::endl;
}
min_required_band_size = pgv->fs_dt*maxvel;
//int temp = std::accumulate(pgv->band_size.begin()+2, pgv->band_size.end(), 0);
int temp = pgv->band_size[2]+pgv->band_size[3]+pgv->band_size[4];
current_band_size = pgv->dxyz_min*static_cast<double>(temp);
if(!pgv->cylindrical){
	if(min_required_band_size > current_band_size){
		std::cout<<pgv->clit<< "####################################"<<std::endl;
		std::cout<<pgv->clit<<"WARNING: ["<<pp->myrank<<"] band sizes might not be sufficient."<<std::endl;
		auto temp = ceil(min_required_band_size/pgv->dxyz_min);
		temp += (-pgv->band_size[2]-pgv->band_size[3]);
		std::cout<<pgv->clit<<"  flow solver timestep requires X-band size of "<<temp<<std::endl;


		std::cout<<pgv->clit<<"     but only have "<<pgv->band_size[4]<<std::endl;
		std::cout<<"fs_dt, dxyz= "<< pgv->fs_dt <<" "<<pgv->dxyz[0]<<'\t'<<pgv->dxyz[1]<<'\t';
		std::cout<<pgv->dxyz[2]<<std::endl;
		std::cout<<pgv->clit<<"##########################"<<std::endl;
	}
}
};

void timestep::timestep_init(){
std::unique_ptr<monitor> pmo(new monitor);
n_wo_reinit = 0;
n_wo_reband = 0;
pmo->lit_monitor_create_file_step("lit_timestep",3);
pmo->lit_monitor_set_header(1, "CFL", 'r');
pmo->lit_monitor_set_header(2, "subcycles", 'i');
pmo->lit_monitor_set_header(3, "nreinit", 'i');
pmo->lit_monitor_select_file("lit_timestep");
pmo->lit_monitor_set_single_values(0.0, 0.0, 0.0);
};


