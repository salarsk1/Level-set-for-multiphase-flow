//Written by Salar Safarkhani

#include "lit_coupler.h"
#include<exception>
void lit_coupler::litRunIteration(const double deltaT, const bool triggerSolutionDump, int arg){
	std::unique_ptr<global_variable> pgv(new global_variable);
	std::unique_ptr<parallel>pp(new parallel);
	std::unique_ptr<timing>pti(new timing);
	std::unique_ptr<monitor>pmo(new monitor);
	std::unique_ptr<advection>padv(new advection);
	std::unique_ptr<timestep>pts(new timestep);
	std::unique_ptr<toolbox>ptool(new toolbox);
	std::unique_ptr<io>pio(new io);
	std::unique_ptr<band>pband(new band);
	bool a1, a2;
	if(pgv->verbose) pp->parallel_barrier();
	if(pgv->verbose && pp->myrank == 0) std::cout<<pgv->clit<<"starting lit_runIteration"<<std::endl;
	pti->lit_timing_start_full();
	pmo->monitor_reinit_steps = 0;
	pgv->fs_dt = static_cast<double>(deltaT);
//	padv->filterVelocity();
	pgv->fs_step += 1;
	pts->lit_timestep();
	ptool->interfaceHolder();
	if(arg == 2){
//		if(triggerSolutionDump) pio->dumpSolution(a1,a2,0);
	}
//	else pio->dumpSolution(a1,a2,0);
	pband->band_monitor();
	pti->lit_timing_monitor();
	pmo->lit_monitor_dump_values_step();
	if(pgv->verbose) pp->parallel_barrier();
	if(pgv->verbose && pp->myrank == 0) std::cout<<pgv->clit<<"Finished lit_runIteration"<<std::endl;
};

void lit_coupler::lit_initialize(int argc, char*argv[]){
	std::unique_ptr<global_variable> pgv(new global_variable);
	std::unique_ptr<parallel>pp(new parallel);
	std::unique_ptr<monitor>pmo(new monitor);
	std::unique_ptr<litparam> ppar(new litparam);
	bool dir_is_there;
	pgv->verbose = true;
	pgv->fs_step = 0;
	pgv->time = 0.0;
	pp->parallel_start(argc, argv);
	pmo->lit_monitor_pre_init();

#ifdef DEBUG_MODE
	if(pp->myrank==0){
		std::cout<<pgv->clit<<"+++++++++++++++++++++++++"<<std::endl;
		std::cout<<pgv->clit<<"+++++++++++++++++++++++++"<<std::endl;
		std::cout<<pgv->clit<<"WARNING: Additional runtime consistency ";
		std::cout<<"will slow"<<std::endl;
		std::cout<<pgv->clit<<"down computation. Undefine DEBUG_MODE in ";
		std::cout<<"Makefile"<<std::endl;
		std::cout<<pgv->clit<<"if these should not be performed"<<std::endl;
		std::cout<<pgv->clit<<"+++++++++++++++++++++++++"<<std::endl;
		std::cout<<pgv->clit<<"+++++++++++++++++++++++++"<<std::endl;
	}
#endif
	if(pp->myrank == 0)
		boost::filesystem::create_directories(ppar->trim(pgv->dump_dir));
	pp->parallel_barrier();
}

void lit_coupler::lit_initialize2(int argc, char *argv[], std::string input_filename){
	bool l;
	std::unique_ptr<parallel>pp(new parallel);
	std::unique_ptr<timing>pti(new timing);
	std::unique_ptr<io>pio(new io);
	std::unique_ptr<litBuffer>plbuf(new litBuffer);
	std::unique_ptr<init>pinit(new init);
	std::unique_ptr<band>pband(new band);
	std::unique_ptr<advection>padv(new advection);
	std::unique_ptr<reinit>prein(new reinit);
	std::unique_ptr<redist>prd(new redist);
	std::unique_ptr<timestep>pts(new timestep);
	std::unique_ptr<toolbox>ptool(new toolbox);
	std::unique_ptr<sg>psg(new sg);
	std::unique_ptr<monitor>pmo(new monitor);
	pti->lit_timing_pre_init();
	pio->io_init();
	std::unique_ptr<bound>pbou(new bound);
//	pio->read_input(argc, argv, input_filename, 3);
	plbuf->litBuffer_m_init();
	pbou->bound_m_init();
	pinit->init_m_init();
	pinit->init_geometry();
	pband->band_m_init();
	padv->advection_m_init();
	prein->reinit_m_init();
	prd->redist_m_init();
	pts->timestep_init();
	ptool->toolbox_m_init();
	psg->sg_m_init();
	pp->parallel_buffer_attach();
	psg->sg_init();
//	pio->dumpSolution(true, 1);
	l = prein->trigger_reinit(10000);
	ptool->diagnostics_init();
	ptool->diagnostics();
	ptool->rotating_cyl_mfix_init();
	pti->lit_timing_post_init();
	pmo->lit_monitor_post_init();
	pband->band_monitor();
	pti->lit_timing_monitor();
	pmo->lit_monitor_dump_values_step();
	pmo->lit_monitor_log("SIMULATION INITIALIZATION");
}

void lit_coupler::lit_finalize(){
	bool skipDumph;
	std::unique_ptr<global_variable> pgv(new global_variable);
	std::unique_ptr<parallel>pp(new parallel);
	std::unique_ptr<band>pband(new band);
	std::unique_ptr<monitor>pmo(new monitor);
	if(pgv->verbose && pp->myrank == 0) std::cout<<"start finalize"<<std::endl;
	skipDumph = false;
	pband->band_m_cleanup();
	if(pgv->verbose && pp->myrank == 0) std::cout<<"end finalize"<<std::endl;
	pmo->lit_monitor_log("Simulation ended.");
	pmo->lit_monitor_finalize();
}
