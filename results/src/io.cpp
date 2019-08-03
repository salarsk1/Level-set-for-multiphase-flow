//Written by Salar Safarkhani

#include<sstream>
#include<iomanip>
#include<fstream>
#include<ios>
#include "io.h"
void io::io_init(){
std::unique_ptr<timing>pti(new timing);
pti->lit_timing_create("io_ensight");
pti->lit_timing_create("io_restart");
};

void io::read_input(int argc, char* argv[], std::string input_filename, int arg){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<bl> pbl(new bl);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<band> pba(new band);
std::string line, chelp, fname, namep, name2;
const bool notnecessary = true;
const bool necessary = false;
int nGhostMin;
std::string::size_type sz, sz1, sz2, sz3, sz4, sz5, sz6, sz7, sz8;
param_t* param;
bool lhelp;
ppar->init_param(argc, argv, input_filename);
pgv->verbose = ppar->get_logical_param("LIT.VERBOSE", false);
//NOHDF is used;
pgv->do_restart = ppar->get_logical_param("LIT.RESTART", false);
pgv->cylindrical = ppar->get_logical_param("LIT.CYLINDRICAL", false);
pgv->rotating_cyl_mfix = ppar->get_logical_param("LIT.ROTATING_CYL_MFIX", false);
pgv->diagnostics_type = ppar->get_string_param("LIT.DIAGNOSTICS", "none", 2);
if(ppar->trim(pgv->diagnostics_type)=="none"|| ppar->trim(pgv->diagnostics_type)=="NONE"|| \
	ppar->trim(pgv->diagnostics_type)=="None" || ppar->trim(pgv->diagnostics_type)=="0"|| \
	ppar->trim(pgv->diagnostics_type)=="F" || ppar->trim(pgv->diagnostics_type)=="False" || \
	ppar->trim(pgv->diagnostics_type)=="false" || ppar->trim(pgv->diagnostics_type)=="FALSE")
		pgv->do_diagnostics = false;
else pgv->do_diagnostics = true;
std::vector<double> tempr;
std::vector<int> tempi;
pgv->xyzs_sg = ppar->get_real3_param("LIT.XYZS_SG", tempr);

pgv->xyze_sg = ppar->get_real3_param("LIT.XYZE_SG", tempr);
pgv->ijkm_gl = ppar->get_integer3_param("LIT.IJKM_GL", tempi);
pgv->ijkm_bl = ppar->get_integer3_param("LIT.IJKM_BL", tempi);
pgv->max_bl = ppar->get_integer_param("LIT.MAX_BL", 512, 2);
pgv->loadBalancer = ppar->get_string_param("LIT.LOAD_BALANCER", "PARMETIS", 2);
lhelp = ppar->get_logical_param("LIT.OLD_LOAD_BAlANCE", false, 2);
if(lhelp) pgv->loadBalancer = "LIT";
chelp = ppar->get_string_param("LIT.LOAD_BALANCE_BAND", "N", 2);
if(ppar->trim(chelp) == "A") pgv->loadBalanceband = 1;
else if(ppar->trim(chelp) == "T") pgv->loadBalanceband = 2;
else if(ppar->trim(chelp) == "N") pgv->loadBalanceband = 3;
else if(ppar->trim(chelp) == "W") pgv->loadBalanceband = 4;
else if(ppar->trim(chelp) == "X") pgv->loadBalanceband = 5;
else{
	pp->litError("read_input", "unknown LIT.LOAD_BALANCE_BAND == ");
	std::cout<<ppar->trim(chelp)<<std::endl;
}
pp->max_load_imbalance = ppar->get_real_param("LIT.MAX_LOAD_IMBALANCE", 1.25, 2);
pgv->schemeReinit = ppar->get_string_param("LIT.SCHEME_REINIT", "WENO-5", 2);
pgv->schemeAdvect = ppar->get_string_param("LIT.SCHEME_ADVECT", "WENO-5", 2);
nGhostMin = 0;
if(ppar->trim(pgv->schemeReinit) == "UPWIND-1") nGhostMin = 1;
else if(ppar->trim(pgv->schemeReinit) == "WENO-3") nGhostMin = 2;
else if(ppar->trim(pgv->schemeReinit) == "WENO-5") nGhostMin = 3;
else{
	pp->litError("read_input", "unknown LIT.SCHEME_REINIT =");
	std::cout<<pgv->schemeReinit<<std::endl;
}
if(ppar->trim(pgv->schemeAdvect) == "UPWIND-1") nGhostMin = std::max(nGhostMin, 1);
else if(ppar->trim(pgv->schemeAdvect) == "WENO-3" || ppar->trim(pgv->schemeAdvect) == "UC-3")
	nGhostMin = std::max(nGhostMin,2);
else if(ppar->trim(pgv->schemeAdvect) == "WENO-5" || ppar->trim(pgv->schemeAdvect) == "UC-5")
	nGhostMin = std::max(nGhostMin, 3);
else if(ppar->trim(pgv->schemeAdvect) == "UC-7") nGhostMin = std::max(nGhostMin, 4);
else if(ppar->trim(pgv->schemeAdvect) == "UC-9") nGhostMin = std::max(nGhostMin, 5);
else if(ppar->trim(pgv->schemeAdvect) == "UC-11") nGhostMin = std::max(nGhostMin, 6);
else {
	pp->litError("read_input","unknown LIT.SCHEME_ADVECT =");
	std::cout<<ppar->trim(pgv->schemeAdvect)<<std::endl;
}
pgv->nghost = ppar->get_integer_param("LIT.NGHOST", nGhostMin, 2);
if(pgv->nghost < nGhostMin){
	std::cout<<"Error: must have at least "<<nGhostMin<<" ghost cells for requested \
					advection/reinitialization schemes"<<std::endl;
	pp->litError("read_input", "insufficient number of requested ghost cells");
}
pgv->sg_periodic[0] = ppar->get_logical_param("LIT.PERIODIC_X", false, 2);
pgv->sg_periodic[1] = ppar->get_logical_param("LIT.PERIODIC_Y", false, 2);
pgv->sg_periodic[2] = ppar->get_logical_param("LIT.PERIODIC_Z", false, 2);
if(pgv->cylindrical && pgv->ijkm_gl[2] > 1) pgv->sg_periodic[2] = true;
if(pgv->verbose && pp->myrank == 0){
	std::cout<<pgv->clit<<" global periodic bcs =";
	for(auto i=0; i<3; i++) std::cout<<pgv->sg_periodic[i]<<'\t';
	std::cout<<std::endl;
}
pgv->delta_width_factor = ppar->get_real_param("LIT.DELTA_WIDTH_FACTOR", 1.5, 2);
pgv->calcCurvatureType  = ppar->get_integer_param("LIT.CALC_CURVATURE_TYPE", pgv->ccDirect, 2);
if(pgv->verbose && pp->myrank == 0){
	if(pgv->calcCurvatureType == pgv->ccNodes)
		std::cout<<pgv->clit<<" curvature calculation at nodes"<<std::endl;
	else if(pgv->calcCurvatureType == pgv->ccDirect)
		std::cout<<pgv->clit<<" curvature calculation by direct surface method"<<std::endl;
	else if(pgv->calcCurvatureType == pgv->ccTriCubic)
		std::cout<<pgv->clit<<" curvature calculation by Chopps method"<<std::endl;
}
pgv->fs_max_drop_vol = ppar->get_real_param("LIT.FS_MAX_DROP_VOL", -1.0, 2);
if(pgv->fs_max_drop_vol < 0.0){
	pgv->fs_max_drop_vol = ppar->get_real_param("LIT.FS_MAX_DROP_R", -1.0, 2);
	if(pgv->fs_max_drop_vol < 0.0)
		pgv->fs_max_drop_vol = 0.0;
	else
		pgv->fs_max_drop_vol = pgv->pi/6.0 * pow(2.0*pgv->fs_max_drop_vol,3);
}
pgv->slender_size = ppar->get_integer_param("LIT.SLENDER_SIZE", 0, 2);
pgv->reinit_solver = ppar->get_string_param("LIT.REINIT_SOLVER", "PDE", 2);
param = ppar->get_param("LIT.HOLDER");
if(param != nullptr){
	if(pp->ierr==0) pgv->holder_name = param->line.substr(param->argv[1]-1); // line 213
	else{
		std::cout<<"Error: interpreting LIT.HOLDER:  "<<ppar->trim(param->line);
		pp->parallel_kill(0);
	}
	pgv->holder = true;

	if(ppar->trim(pgv->holder_name) == "plane" || ppar->trim(pgv->holder_name) == "Plane"||\
		ppar->trim(pgv->holder_name) == "PLANe"){
		if(pp->ierr == 0){
			// line 222; ask Dr Herrmann how it is read in fortran
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->holder_xyzs[0] = std::stod(strtemp, &sz);
			pgv->holder_xyzs[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->holder_xyzs[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->holder_xyze[0] = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->holder_xyze[1] = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->holder_xyze[2] = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4));
		}
		else{
			std::cout<<"Error: interpreting interface holder plane: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.HOLDER plane [start 0 to 2] [end 0 to 2]"<<std::endl;
			pp->parallel_kill(0);
		}
		
	}
	else if(ppar->trim(pgv->holder_name)=="nozzle"||ppar->trim(pgv->holder_name)=="Nozzle"||
			  ppar->trim(pgv->holder_name)=="NOZZLE"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->holder_xyzs[0] = std::stod(strtemp, &sz);
			pgv->holder_xyze[0] = std::stod(strtemp.substr(sz), &sz1);
			pgv->holder_center[1] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->holder_center[2] = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->holder_radius = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->holder_dr = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4));

		}
		else{
			std::cout<<"Error: iterpreting interface holder nozzle:  "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.HOLDER nozzle [xstart] [xend] [centery] [centerx] radius radialextend"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->holder_name)=="nozzle_z"||ppar->trim(pgv->holder_name)=="Nozzle_z"||
			  ppar->trim(pgv->holder_name)=="NOZZLE_Z"){

		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->holder_xyzs[2] = std::stod(strtemp, &sz);
			pgv->holder_xyze[2] = std::stod(strtemp.substr(sz), &sz1);
			pgv->holder_center[0] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->holder_center[1] = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->holder_radius = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->holder_dr = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4));
		}
		else{
			std::cout<<"Error: iterpreting interface holder nozzle z:  "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.HOLDER nozzle [zstart] [zend] [centerx] [centery] radius radialextend"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->holder_name)=="nozzle_y"||ppar->trim(pgv->holder_name)=="Nozzle_y"||
			  ppar->trim(pgv->holder_name)=="NOZZLE_Y"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->holder_xyzs[1] = std::stod(strtemp, &sz);
			pgv->holder_xyze[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->holder_center[0] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->holder_center[2] = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->holder_radius = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->holder_dr = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4));
		}
		else{
			std::cout<<"Error: iterpreting interface holder nozzle y:  "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.HOLDER nozzle [zstart] [zend] [centerx] [centery] radius radialextend"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->holder_name)=="sphere_cap"||ppar->trim(pgv->holder_name)=="Sphere_Cap"||
			  ppar->trim(pgv->holder_name)=="SPHERE_CAP"){
		if(pp->ierr == 0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->holder_xyzs[2] = std::stod(strtemp, &sz);
			pgv->holder_xyze[2] = std::stod(strtemp.substr(sz), &sz1);
			pgv->holder_dt = std::stod(strtemp.substr(sz+sz1));
		}
		else{
			std::cout<<"Error: iterpreting interface holder sphere_cap:  "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.HOLDER nozzle [zstart] [zend] releasetime"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->holder_name)=="sphere_cap_cyl"||ppar->trim(pgv->holder_name)=="Sphere_Cap_Cyl"||
			  ppar->trim(pgv->holder_name)=="SPHERE_CAP_CYL"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->holder_xyzs[0] = std::stod(strtemp, &sz);
			pgv->holder_xyze[0] = std::stod(strtemp.substr(sz), &sz1);
			pgv->holder_dt = std::stod(strtemp.substr(sz+sz1));		
		}
		else{
			std::cout<<"Error: iterpreting interface holder sphere_cap:  "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.HOLDER nozzle [zstart] [zend] releasetime"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->holder_name)=="circle_cap"||ppar->trim(pgv->holder_name)=="Circle_Cap"||
			  ppar->trim(pgv->holder_name)=="CIRCLE_CAP"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->holder_xyzs[1] = std::stod(strtemp, &sz);
			pgv->holder_xyze[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->holder_dt = std::stod(strtemp.substr(sz+sz1));				
		}
		else{
			std::cout<<"Error: iterpreting interface holder sphere_cap:  "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.HOLDER nozzle [zstart] [zend] releasetime"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else pgv->holder = false;
}
std::vector<double> xyz_temp(3);
for(auto i=0; i<3; i++) xyz_temp[i] = pgv->xyzs_sg[i];
pgv->xyzs_init_bb = ppar->get_real3_param("LIT.XYZS_INIT_BB", xyz_temp, 2);
for(auto i=0; i<3; i++) xyz_temp[i] = pgv->xyze_sg[i];
pgv->xyze_init_bb = ppar->get_real3_param("LIT>XYZE_INIT_BB", xyz_temp, 2);
param = ppar->get_param("LIT.INIT_SHAPE");
if(param != nullptr){
	if(pp->ierr==0) pgv->init_shape = param->line.substr(param->argv[1]-1);
	else{
		std::cout<<"Error: interpreting LIT.INIT_SHAPE:  "<<ppar->trim(param->line);
		pp->parallel_kill(0);
	}
	if(ppar->trim(pgv->init_shape) == "notched_circle"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_width = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->init_height = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4));
		}
		else{
			std::cout<<"Error: interpeting init shape notched_circle: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE notched_circle center[0 to 2] [radius] [width] [height]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="circle"||ppar->trim(pgv->init_shape)=="column"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2));
		}
		else{
			std::cout<<"Error: interpreting init shape circle: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE [circle,column] centrer[0 to 2] [radius]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape) == "circle_xz"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2));
		}
		else{
			std::cout<<"Error: interpreting init shape circle_xz: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE [circle_xz] center[0 to 3] [radius]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="column_cap_z"||ppar->trim(pgv->init_shape)=="column_cap_y"|| \
			  ppar->trim(pgv->init_shape)=="column_cap_x"){

		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_length = std::stod(strtemp.substr(sz+sz1+sz2+sz3));
		}
		else{
			std::cout<<"Error: interpreting init shape column_cap_z: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE column_cap_[x,y,z] center[0 to 2] [radius] [length]"<<std::endl;
			pp->parallel_kill(0);
		}

	}
	else if(ppar->trim(pgv->init_shape)=="column_pool"){
		if(pp->ierr == 0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_height = std::stod(strtemp.substr(sz+sz1+sz2+sz3));
		}
		else{
			std::cout<<"Error: interpreting init shape column_pool: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE column_pool center[0 to 2] [radius] [pool height from y-bottom]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="sphere_pool"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_height = std::stod(strtemp.substr(sz+sz1+sz2+sz3));
		}
		else{
			std::cout<<"Error: interpreting init shape sphere_pool: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE sphere_pool center[0 to 2] [radius] [pool height from y-bottom]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="random_circle"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_amplitude = std::stod(strtemp.substr(sz+sz1+sz2+sz3));		
		}
		else{
			std::cout<<"Error: interpreting init shape random_circle: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE random_circle center[0 to 2] [radius] [amplitude]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="random_jet"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			fname = strtemp.substr(sz+sz1+sz2+sz3);
		}
		else{
			std::cout<<"Error: interpreting init shape random jet: "<<ppar->trim(param->line);
			std::cout<<"syntax is LIT.INIT_SHAPE random_jet center[0 to 2] [radius] [init_random_filename]"<<std::endl;
			pp->parallel_kill(0);
		}
		std::ifstream myfile; // line 356 of fortran ; check the format
		myfile.open(ppar->trim(fname));
		if(!myfile.is_open()) std::cout<<"FILE IS NOT OPEN line: 367 of io.cpp"<<std::endl;
		if(myfile.good()) myfile>>pgv->n1_random>>pgv->n2_random;
		if(pgv->a1_random != nullptr){
			if(pp->ierr==nullptr){
				delete [] pgv->a1_random;
				delete [] pgv->a2_random;
				delete [] pgv->phi1_random;
				delete [] pgv->phi2_random;
			}
			else{
				pp->litError("read_input","deallocation error for a1_random,a2_random,phi1_random,phi2_random");
			}
		}
		pgv->a1_random = new double[pgv->n1_random];
		pgv->phi1_random = new double[pgv->n1_random];
		pgv->a2_random = new double[pgv->n2_random];
		pgv->phi2_random = new double[pgv->n2_random];
		for(auto ig=0; ig<pgv->n1_random; ig++){
			myfile>>pgv->a1_random[ig]>>pgv->phi1_random[ig];
		}
		for(auto ig=0; ig<pgv->n2_random; ig++){
			myfile>>pgv->a2_random[ig]>>pgv->phi2_random[ig];
		}
		myfile.close();
	}
	else if(ppar->trim(pgv->init_shape)=="random_drops"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_mod = std::stoi(strtemp, &sz);
			pgv->init_number = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_xyz_min[0] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_xyz_max[0] = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_xyz_min[1] = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->init_xyz_max[1] = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4), &sz5);
			pgv->init_xyz_min[2] = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4+sz5), &sz6);
			pgv->init_xyz_max[2] = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4+sz5+sz6), &sz7);
			pgv->init_radius_min = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4+sz5+sz6+sz7), &sz8);
			pgv->init_radius_max = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4+sz5+sz6+sz7+sz8));
		}
		else{
			std::cout<<"Error: interpreting init shape random_drops: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE random_drops [read file(1) or not(0)]";
			std::cout<<"[number] [xmin,xmax,ymin,ymax,zmin,zmax] [minimum radius] [maximum radius]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="ring"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_width = std::stod(strtemp.substr(sz+sz1+sz2+sz3));
		}
		else{
			std::cout<<"Error: interpreting init shape ring: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE ring center[0 to 2] [radius] [width]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="rod"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_length = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_width = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->init_angle = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4));
		}
		else{
			std::cout<<"Error: interpreting init shape rod: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE rod center[0 to 2] [length] [width] [angle in deg]"<<std::endl;
			pp->parallel_kill(0);
		}
		pgv->init_angle  = pgv->init_angle/180.0 * pgv->pi;
	}
	else if(ppar->trim(pgv->init_shape) == "sphere"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2));
		}
		else{
			std::cout<<"Error: interpreting init shape sphere: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE sphere center[0 to 2] [radius]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape) == "sphere2"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_center2[0] = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->init_center2[1] = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4), &sz5);
			pgv->init_center2[2] = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4+sz5), &sz6);
			pgv->init_radius2 = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4+sz5+sz6));
		}
		else{
			std::cout<<"Error: interpreting init shape sphere2: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE sphere2 center[0 to 2] [radius]";
			std::cout<<"center2[0 to 2] [radius2]"<<std::endl;
			pp->parallel_kill(0);		
		}
	}
	else if(ppar->trim(pgv->init_shape) == "sphere_cyl"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2));
		}
		else{
			std::cout<<"Error: interpreting init shape sphere_cyl: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE sphere_cyl center[0 to 2] [radius]"<<std::endl;
			pp->parallel_kill(0);		
		}
	}
	else if(ppar->trim(pgv->init_shape)=="disc_cyl"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_length = std::stod(strtemp.substr(sz+sz1+sz2+sz3));
		}
		else{
			std::cout<<"Error: interpreting init shape disc_cyl: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE disc_cyl center[0 to 2] [radius] [length]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="plane"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_mod = std::stoi(strtemp.substr(sz+sz1+sz2));
		}
		else{
			std::cout<<"Error: interpreting init shape plane: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE plane center[0 to 2] [mod]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="sine"||ppar->trim(pgv->init_shape)=="cosine" || \
			  ppar->trim(pgv->init_shape)=="siney"||ppar->trim(pgv->init_shape)=="cosiney" || \
			  ppar->trim(pgv->init_shape)=="cosine3D"||ppar->trim(pgv->init_shape)=="sine_column"|| \
			  ppar->trim(pgv->init_shape)=="siney_column"||ppar->trim(pgv->init_shape)=="RT3D"){

		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_amplitude = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_wavelength = std::stod(strtemp.substr(sz+sz1+sz2+sz3));
		}
		else{
			std::cout<<"Error: interpreting init shape sin,cosine,cosine3D,siney,cosiney";
			std::cout<<"sine_column,siney_column,RT3D: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE sine center[0 to 2] [amplitude] [wavelength]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape) == "RMI"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_amplitude = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_wavelength = std::stod(strtemp.substr(sz+sz1+sz1+sz3));
		}
		else{
			std::cout<<"Error: interpreting init shape RMI: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE RMI center[0 to 2] [amplitude*k] [wavelength]"<<std::endl;
			pp->parallel_kill(0);
		}
		pgv->init_amplitude = pgv->init_amplitude*pgv->init_wavelength/(2.0*pgv->pi);
	}
	else if(ppar->trim(pgv->init_shape)=="sheet"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_normal[0] = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_normal[1] = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->init_normal[2] = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4), &sz5);
			pgv->init_width = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4+sz5));
		}
		else{
			std::cout<<"Error: interpreting init shape sheet: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE RMI center[0 to 2] normal[0 to 2] [width]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="random_sheet"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_normal[0] = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_normal[1] = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->init_normal[2] = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4), &sz5);
			pgv->init_width = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4+sz5), &sz6);
			pgv->init_amplitude = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4+sz5+sz6));
		}
		else{
			std::cout<<"Error: interpreting init shape random_sheet: "<<ppar->trim(param->line)<<std::cout;
			std::cout<<"syntax is LIT.INIT_SHAPE random_sheet center[0 to 2] normal[0 to 2] [width] [amplitude]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="deformed_column"||ppar->trim(pgv->init_shape)=="deformed_sphere"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_amplitude = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->init_mod = std::stoi(strtemp.substr(sz+sz1+sz2+sz3+sz4));
		}
		else{
			std::cout<<"Error: interpreting init shape deformed_column,deformed_shape: ";
			std::cout<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE deformed_column center[0 to 2] [radius] [amplitude] [mod]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="ellipse"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_a = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_b = std::stod(strtemp.substr(sz+sz1+sz2+sz3));
		}
		else{
			std::cout<<"Error: interpreting init shape ellipse (x/a)^2+(y/b)^2=1: ";
			std::cout<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE ellipse center[0 to 2] [a] [b]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="milk_crown"|| \
			  ppar->trim(pgv->init_shape)=="bursting_bubble"|| \
			  ppar->trim(pgv->init_shape)=="bursting_bubble_3D"|| \
			  ppar->trim(pgv->init_shape)=="bursting_bubble_2D"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_height = std::stod(strtemp.substr(sz+sz1+sz2+sz3));
		}
		else{
			std::cout<<"Error: interpreting init shape milk_crown, bursting_bubble,";
			std::cout<<"bursting_bubble_3D,bursting_bubble_2D: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE milk_crown center[0 to 2] [radius] [height]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="bursting_bubble_rim_3D"|| \
			  ppar->trim(pgv->init_shape)=="bursting_bubble_rim_2D"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_height = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->init_hole_radius = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4), &sz5);
			pgv->init_rim_radius = std::stod(strtemp.substr(sz+sz1+sz2+sz3+sz4+sz5));
		}
		else{
			std::cout<<"Error: interpreting init shape bursting_bubble_rim_3D,";
			std::cout<<"bursting_bubble_rim_2: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE bursting_bubble_rim_3D center[0 to 2] [radius]";
			std::cout<<"[height] [hole_radius] [rim_radius]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="droplens"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_radius = std::stod(strtemp, &sz);
			pgv->init_height = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_radius2 = std::stod(strtemp.substr(sz+sz1));
		}
		else{
			std::cout<<"Error: interpreting init droplens: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE droplens [hole redius] [hole height] [drop radius]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="ligament_2D"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_length = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_width = std::stod(strtemp.substr(sz+sz1+sz2+sz3));
		}
		else{
			std::cout<<"Error: interpreting init ligament_2D: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE ligament_2D center[0 to 2] [length] [width]"<<std::endl;
			pp->parallel_kill(0);	
		}
	}
	else if(ppar->trim(pgv->init_shape)=="rayleigh"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_amplitude = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->init_mod = std::stoi(strtemp.substr(sz+sz1+sz2+sz3+sz4));
		}
		else{
			std::cout<<"Error: interpreting init rayleigh: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE rayleigh center[0 to 2] [radius] [amplitude] [mod]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
	else if(ppar->trim(pgv->init_shape)=="rayleigh_axi"){
		if(pp->ierr==0){
			string strtemp = param->line.substr(param->argv[2]-1);
			pgv->init_center[0] = std::stod(strtemp, &sz);
			pgv->init_center[1] = std::stod(strtemp.substr(sz), &sz1);
			pgv->init_center[2] = std::stod(strtemp.substr(sz+sz1), &sz2);
			pgv->init_radius = std::stod(strtemp.substr(sz+sz1+sz2), &sz3);
			pgv->init_amplitude = std::stod(strtemp.substr(sz+sz1+sz2+sz3), &sz4);
			pgv->init_mod = std::stoi(strtemp.substr(sz+sz1+sz2+sz3+sz4));
		}
		else{
			std::cout<<"Error: interpreting init rayleigh: "<<ppar->trim(param->line)<<std::endl;
			std::cout<<"syntax is LIT.INIT_SHAPE rayleigh_axi center[0 to 2] [radius] [amplitude] [mod]"<<std::endl;
			pp->parallel_kill(0);
		}
	}
}
param = ppar->get_param("WRITE_LIT_STEP");
pgv->dump_dstep_fs = -1;
while(param != nullptr){
	if(pp->ierr==0){
		string strtemp = param->line.substr(param->argv[2]-1);
		std::vector<string> vectortemp = split(strtemp);
		namep = vectortemp[0];
		pgv->case_name = vectortemp[1];
		pgv->dump_dstep_fs = std::stoi(vectortemp[2]);
		name2 = vectortemp[3];
	}
	else{
		std::cout<<pgv->clit<<'\t'<<"Error: interpreting WRITE_LIT_STEP: "<<ppar->trim(param->line)<<std::endl;
		pp->parallel_kill(0);
	}
	if(name2=="T-BAND"||name2=="T-band"||name2=="T-Band"||name2=="t-band"){
		pgv->dump_band = "T";
	}
	if(name2=="N-BAND"||name2=="N-band"||name2=="N-Bnad"||name2=="n-band"){
		pgv->dump_band = "N";
	}
	if(name2=="X-BAND"||name2=="X-band"||name2=="X-Band"||name2=="x-band"){
		pgv->dump_band = "X";
	}
	if(name2=="Z-BAND"||name2=="Z-band"||name2=="Z-Band"||name2=="Z-band"){
		pgv->dump_band = "Z";
	}
	else pgv->dump_band = "A";
	nVarDump = param->argc-5;
	auto sizediff = 0;
	for(auto i=5; i<param->argc; i++){
		if(pp->ierr==0){
			sizediff = param->argv[i+1]-param->argv[i];
			varDumpName[i-5] = param->line.substr(param->line[param->argv[i]-1], sizediff);
		}
	}
	if(namep=="TECPLOT") pgv->lit_die("Tecplot output not implemented in LIT");
	ppar->set_next_param(param, "WRITE_LIT_STEP");
}
param = ppar->get_param("WRITE_LIT_NGA");
pgv->dump_dtime_fs = -1.0e20;
while(param != nullptr){
	if(pp->ierr==0){
		string strtemp = param->line.substr(param->argv[2]-1);
		std::vector<string> vectortemp = split(strtemp);
		namep = vectortemp[0];
		pgv->case_name = vectortemp[1];
		pgv->dump_dtime_fs = std::stod(vectortemp[2]);
		name2 = vectortemp[3];
	}
	else{
		std::cout<<pgv->clit<<'\t'<<"Error: interpreting WRITE_LIT_NGA: ";
		std::cout<<ppar->trim(param->line)<<std::endl;
		std::cout<<pgv->clit<<'\t'<<"   must be: WRITE_LIT_NGA(ENSIGHT/TECPLOT) nameprefix ";
		std::cout<<"band-name varnames"<<std::endl;
		pp->parallel_kill(0);
	}
	if(name2=="T-BAND"||name2=="T-band"||name2=="T-Band"||name2=="t-band")
		pgv->dump_band="T";
	if(name2=="N-BAND"||name2=="N-band"||name2=="N-Band"||name2=="n-band")
		pgv->dump_band="N";
	if(name2=="X-BAND"||name2=="X-band"||name2=="X-Band"||name2=="x-band")
		pgv->dump_band="X";
	if(name2=="Z-BAND"||name2=="Z-band"||name2=="Z-Band"||name2=="z-band")
		pgv->dump_band="Z";
	else pgv->dump_band = "A";
	nVarDump = param->argc-4;
	auto sizediff=0;
	for(auto i=4; i<param->argc; i++){
		if(pp->ierr==0){
			sizediff = param->argv[i+1]-param->argv[i];
			varDumpName[i-4] = param->line.substr(param->line[param->argv[i]-1], sizediff);
		}
	}
	if(namep=="TECPLOT") pgv->lit_die("Tecplot output not implemented in LIT");
	ppar->set_next_param(param, "WRITE_LIT_NGA");
}
param = ppar->get_param("WRITE_LIT_CELL_STEP");
pgv->dump_dstep_cell_fs = -1;
while(param != nullptr){
	if(pp->ierr==0){
		string strtemp = param->line.substr(param->argv[2]-1);
		std::vector<string> vectortemp = split(strtemp);
		namep = vectortemp[0];
		pgv->case_name = vectortemp[1];
		pgv->dump_dstep_cell_fs = std::stoi(vectortemp[2]);
		name2 = vectortemp[3];
	}
	else{
		std::cout<<pgv->clit<<'\t'<<"Error: interpreting WRITE_LIT_CELL_STEP: ";
		std::cout<<ppar->trim(param->line)<<std::endl;
		pp->parallel_kill(0);
	}
	if(name2=="T-BAND"||name2=="T-band"||name2=="T-Band"||name2=="t-band")
		pgv->dump_band_cell="T";
	if(name2=="N-BAND"||name2=="N-band"||name2=="N-Band"||name2=="n-band")
		pgv->dump_band_cell="N";
	if(name2=="X-BAND"||name2=="X-band"||name2=="X-Band"||name2=="x-band")
		pgv->dump_band_cell="X";
	if(name2=="Z-BAND"||name2=="Z-band"||name2=="Z-Band"||name2=="z-band")
		pgv->dump_band_cell="Z";
	else pgv->dump_band_cell = "A";
	nVarDump = param->argc-5;
	auto sizediff = 0;
	for(auto i=5; i<param->argc; i++){
		if(pp->ierr==0){
			sizediff = param->argv[i+1]-param->argv[i];
			varCellDumpName[i-5] = param->line.substr(param->line[param->argv[i]-1], sizediff);
		}
	}
	if(namep == "TECPLOT") pgv->lit_die("Tecplot output not implemented in LIT");
	ppar->set_next_param(param, "WRITE_LIT_STEP");
}
param = ppar->get_param("WRITE_LIT_CELL_NAME");
pgv->dump_dtime_cell_fs == -1.0e-20;
while(param != nullptr){
	pgv->lit_die("Write_LIT_CELL_TIME");
	if(pp->ierr==0){
		string strtemp = param->line.substr(param->argv[2]-1);
		std::vector<string> vectortemp = split(strtemp);
		namep = vectortemp[0];
		pgv->case_name = vectortemp[1];
		pgv->dump_dtime_cell_fs = std::stoi(vectortemp[2]);
		name2 = vectortemp[3];
	}
	else{
		std::cout<<pgv->clit<<'\t'<<"Error: interpreting WRITE_CELL_TIME: "<<ppar->trim(param->line)<<std::endl;
		pp->parallel_kill(0);
	}
	if(name2=="T-BAND"||name2=="T-band"||name2=="T-Band"||name2=="t-band")
		pgv->dump_band_cell="T";
	if(name2=="N-BAND"||name2=="N-band"||name2=="N-Band"||name2=="n-band")
		pgv->dump_band_cell="N";
	if(name2=="X-BAND"||name2=="X-band"||name2=="X-Band"||name2=="x-band")
		pgv->dump_band_cell="X";
	if(name2=="Z-BAND"||name2=="Z-band"||name2=="Z-Band"||name2=="z-band")
		pgv->dump_band_cell="Z";
	else pgv->dump_band_cell = "A";
	nVarDump = param->argc-5;
	auto sizediff = 0;
	for(auto i=5; i<param->argc; i++){
		if(pp->ierr==0){
			sizediff = param->argv[i+1]-param->argv[i];
			varCellDumpName[i-5] = param->line.substr(param->line[param->argv[i]-1], sizediff);
		}
	}
	if(namep == "TECPLOT") pgv->lit_die("Tecplot output not implemented in LIT");
	ppar->set_next_param(param, "WRITE_LIT_CELL_TIME");
}
};

void io::dumpSolution(bool newv, bool force, int arg){
bool newh, forceh, timeTrigger;
double twriteNext;
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<bl> pbl(new bl);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<band> pba(new band);
if(arg==1 || arg == 2) newh = newv;
else if(arg==3 || arg == 0) newh = false;
if(arg == 2 || arg == 3) forceh = force;
else if(arg == 1 || arg == 0) forceh = false;
if(pp->myrank == 0){
	twriteNext = static_cast<double>(round(pgv->time/pgv->dump_dtime_fs))*pgv->dump_dtime_fs;
	if(abs(pgv->time-twriteNext) <= 1.0e-6*pgv->fs_dt) timeTrigger=true;
	else timeTrigger = false;
}
pp->parallel_BCast(timeTrigger);
nVarDump = 2;
if(newh||forceh||(pgv->dump_dstep_fs>0 && pgv->fs_step%pgv->dump_dstep_fs==0)||timeTrigger){
	if(nVarDump>0){
		ptool->interfaceHolder();
		pmo->lit_monitor_log("Dumping Ensight data");
		if(arg == 1|| arg == 2){
			dumpEnsight(ppar->trim(pgv->case_name), pgv->dump_band, newv, 3);
		}
		else{
			dumpEnsight(ppar->trim(pgv->case_name), pgv->dump_band, newv);
		}
		
	}
}
if(pp->myrank==0){
twriteNext = static_cast<double>(pgv->time/pgv->dump_dtime_cell_fs)/pgv->dump_dtime_cell_fs;
	if(abs(pgv->time - twriteNext) <= 1.0e-6*pgv->fs_dt){
		timeTrigger = true;
	}
	else{
		timeTrigger = false;
	}
}
pp->parallel_BCast(timeTrigger);
if(newh || forceh || (pgv->dump_dstep_cell_fs>0 && pgv->fs_step%pgv->dump_dstep_cell_fs==0)|| timeTrigger){
	if(nVarCellDump > 0){
		pmo->lit_monitor_log("Dumping Ensight cell data");
		if(arg ==1 || arg==2)
			dumpEnsightCell(ppar->trim(pgv->case_name), pgv->dump_band, newv, 3);
		else{
			dumpEnsightCell(ppar->trim(pgv->case_name),pgv->dump_band_cell, newv);
		}
	}
}
};






void io::dumpEnsight(const std::string name, const std::string band1, const bool newv, int arg){
int nnodes, nelem, nelem_l, inode, which_Band, inode_offset, nnodes_l, inth;
int ielem,_offset, info, strlen, koff, ielem_offset;
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<bl> pbl(new bl);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<band> pba(new band);
MPI_File ifile;
std::vector<int> nhelp(pp->nprocs);
bool exists, append, fileExists;
std::string errsrtr, filename, fname2, chelp6;
std::string line(80, ' '); // was null
double * coord;
int **elem;
int * nodeIc;
bool *nodeFlag;
double *kappa, *marker, *S;
double theta;
pti->lit_timing_start("io_ensight");
if(pgv->verbose && pp->myrank == 0){
	std::cout<<pgv->clit<<"Starting dump_Gnodes_Ensight of file "<<ppar->trim(name)<<std::endl;
}
size_t nodeFlag_S;
nodeFlag = plbuf->getL1Buffer(nodeFlag_S);
for(auto i=0; i<nodeFlag_S; i++) nodeFlag[i] = false;
size_t nodeIc_S;
nodeIc   = plbuf->getI1Buffer(nodeIc_S);
for(auto i=0; i<nodeIc_S; i++) nodeIc[i] = 0;
if(ppar->trim(band1) == "A") 		 which_Band = 1;
else if(ppar->trim(band1) == "T") which_Band = 2;
else if(ppar->trim(band1) == "N") which_Band = 3;
else if(ppar->trim(band1) == "W") which_Band = 4;
else if(ppar->trim(band1) == "X") which_Band = 5;
else{
	std::cout<<pgv->clit<<"unknown band name for dump_Gnodes_Ensight: "<<ppar->trim(band1)<<std::endl;
	pp->parallel_kill(1);
}
if(pp->myrank == 0){
	if(arg == 1 || arg == 2){
		if(newv){
			write_new_case_file(name, 1);
			append = true;
		}
		else{
			append = true;
		}
	}
	else{
		append = true;
	}
	if(append){
		exists = exists_file(ppar->trim(pgv->dump_dir)+ppar->trim(name)+".case");
		if(!exists) write_new_case_file(name, 1);
		else append_case_file(name, 1);
	}
}
pp->ierr = MPI_Bcast(&dump_nfiles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
#ifdef WITH_CDP
if(pgv->dump_dtime_fs > 0.0){
	chelp6 = std::to_string(dump_nfiles-1); // correct the format
}
else{
	chelp6 = std::to_string(pgv->fs_step); // corect the format
}
#else
chelp6 = std::to_string(dump_nfiles-1); // correct the format
#endif
std::stringstream ssss;
ssss<<std::setw(6)<<std::setfill('0')<<chelp6;
std::string sss = ssss.str();
filename = ppar->trim(pgv->dump_dir)+"sol_"+ppar->trim(name)+"_"+sss+".geo";
if(pp->myrank == 0){
	fileExists = exists_file(ppar->trim(filename));
	if(fileExists){
		if(pgv->verbose) std::cout<<pgv->clit<<"deleting geo"<<std::endl;
		std::string temp = ppar->trim(filename);
		char* temp_mpi = const_cast<char*>(temp.c_str());
		pp->ierr=MPI_File_delete(temp_mpi, pp->MPI_IO_LIT_FILE_HINT);
	}
	if(pgv->verbose) std::cout<<pgv->clit<<"writing geo "<<ppar->trim(filename)<<std::endl;
}
filename = ppar->trim(pp->file_prefix)+ppar->trim(filename);
char* temp_mpi = const_cast<char*>(ppar->trim(filename).c_str());
pp->ierr=MPI_File_open(MPI_COMM_WORLD, temp_mpi, MPI_MODE_WRONLY+MPI_MODE_CREATE, pp->MPI_IO_LIT_FILE_HINT, &ifile);
if(pp->myrank == 0){
	line = "C Binary";
	line.append(72, ' ');
	pp->ierr=MPI_File_write(ifile, const_cast<char *>(line.data()), 80, MPI_CHARACTER, &pp->status);

	line = ppar->trim(band1)+"-Band nodes";
	line.append(68, ' ');
	pp->ierr=MPI_File_write(ifile, const_cast<char *>(line.data()), 80, MPI_CHARACTER, &pp->status);

	line = "";
	line.append(80,' ');
	pp->ierr=MPI_File_write(ifile, const_cast<char *>(line.data()), 80, MPI_CHARACTER, &pp->status);

	line = "node id off";
	line.append(69, ' ');
	pp->ierr=MPI_File_write(ifile, const_cast<char *>(line.data()), 80, MPI_CHARACTER, &pp->status);

	line = "element id off";
	line.append(66, ' ');
	pp->ierr=MPI_File_write(ifile, const_cast<char *>(line.data()), 80, MPI_CHARACTER, &pp->status);

	line = "part";
	line.append(76, ' ');
	pp->ierr=MPI_File_write(ifile, const_cast<char *>(line.data()), 80, MPI_CHARACTER, &pp->status);

	inth = 1;
	pp->ierr=MPI_File_write(ifile,&inth,1, MPI_INTEGER, &pp->status);

	line = ppar->trim(band1)+"-Band coordinates";
	line.append(62, ' ');
	pp->ierr=MPI_File_write(ifile, const_cast<char *>(line.data()), 80, MPI_CHARACTER, &pp->status);

	line = "coordinates";
	line.append(69, ' ');
	pp->ierr=MPI_File_write(ifile, const_cast<char *>(line.data()), 80, MPI_CHARACTER, &pp->status);
}
pp->disp_header = 8*80*1+1*4;
if(pgv->cylindrical){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinBand[which_Band-1]; ic++){
			auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
			auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
			auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
			for(auto ii=i-1; ii<i+2; ii++){
				for(auto jj=j-1; jj<j+2; jj++){
					for(auto kk=k-1; kk<k+2; kk++){
						if(pgv->i2c[ii][jj][kk] > 0){
							if(ii+pgv->b->imino_<1 && pgv->sg_periodic[0]) continue;
							if(jj+pgv->b->jmino_<1 && pgv->sg_periodic[1]) continue;
							if(ii+pgv->b->imino_>pgv->ijkm_gl[0]&&pgv->sg_periodic[0])continue;
							if(jj+pgv->b->jmino_>pgv->ijkm_gl[1]&&pgv->sg_periodic[1])continue;
							if(kk+pgv->b->kmino_>pgv->ijkm_gl[2]&&pgv->sg_periodic[2])continue;
							nodeIc[pgv->ioff+pgv->i2c[ii][jj][kk]-1]=1;
						}
					}
				}
			}
		}
	}
}
else{
	if(pgv->ijkm_gl[2] == 1) koff=0;
	else koff = 1;
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinBand[which_Band-1]; ic++){
			auto i=pgv->Gn[ic].ijk[0]-pgv->b->imino_;	
			auto j=pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
			auto k=pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
			for(auto ii=i-1; ii<i+2; ii++){
				for(auto jj=j-1; jj<j+2; jj++){
					for(auto kk=k-koff; kk<k+koff+1; kk++){
						if(pgv->i2c[ii][jj][kk] > 0){
							if(ii+pgv->b->imino_<1 && pgv->sg_periodic[0]) continue;
							if(jj+pgv->b->jmino_<1 && pgv->sg_periodic[1]) continue;
							if(kk+pgv->b->kmino_<1 && pgv->sg_periodic[2]) continue;
							if(ii+pgv->b->imino_>pgv->ijkm_gl[0]&&pgv->sg_periodic[0]) continue;
							if(jj+pgv->b->jmino_>pgv->ijkm_gl[1]&&pgv->sg_periodic[1]) continue;
							if(kk+pgv->b->kmino_>pgv->ijkm_gl[2]&&pgv->sg_periodic[2]) continue;
							nodeIc[pgv->ioff+pgv->i2c[ii][jj][kk]-1]=1;
						}
					}
				}
			}
		}
	}
}
nnodes_l = 0;
for(auto i_temp=0; i_temp<nodeIc_S; i_temp++){
	if(nodeIc[i_temp]==1) nnodes_l++;
}
pp->parallel_all_sum(nnodes_l, nnodes);
if(pp->myrank == 0)
	pp->ierr=MPI_File_write(ifile, &nnodes, 1, MPI_INTEGER, &pp->status);
pp->disp_header += 4;
//pp->parallel_all_gather(nnodes_l, nhelp);
if(pp->myrank > 0)
	inode_offset = nnodes_l; //std::accumulate(nhelp.begin(), nhelp.begin()+pp->myrank, 0);
else inode_offset = 0;
if(pp->ierr==0) coord = new double[nnodes_l];
else pp->litError("dumpEnsight","allocation error for coord");

for(auto id=0; id<3; id++){
	inode = -1;
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinZ; ic++){
			if(nodeIc[pgv->ioff+ic]>0){
				inode++;
				nodeFlag[pgv->ioff+ic] = true;
				if(id==0){
						nodeIc[pgv->ioff+ic]=inode+inode_offset+1;
						coord[inode] = pgv->xc[pgv->Gn[ic].ijk[0]-(1-pgv->nghost)];
				}
				else if(id == 1){
						if(pgv->cylindrical)
							coord[inode] = pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)]*\
												cos(pgv->zc[pgv->Gn[ic].ijk[2]-(1-pgv->nghost)]);
						else
							coord[inode] = pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)];
				}
				else if(id == 2){
						if(pgv->cylindrical)
							coord[inode]=pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)]*\
											 sin(pgv->zc[pgv->Gn[ic].ijk[2]-(1-pgv->nghost)]);
						else coord[inode] = pgv->zc[pgv->Gn[ic].ijk[2]-(1-pgv->nghost)];
		
				}
			}
		}
	}
	pp->disp=pp->disp_header+(static_cast<int>(id)*static_cast<int>(nnodes)+static_cast<int>(inode_offset))*4;

	pp->ierr==MPI_File_write_at(ifile, pp->disp, coord, nnodes_l, MPI_DOUBLE, &pp->status);


	pp->MPI_ierr_Handler(1);
}
if(pp->ierr==0&& coord != nullptr) delete [] coord;
else pp->litError("dumpEnsight","deallocation error for coord");
pbou->updateGhostI1(nodeIc);
if(pp->myrank==0){
	pp->disp = pp->disp_header + nnodes*3*4;
	if(pgv->ijkm_gl[2]==1){
		line = "quad4";
		line.append(75, ' ');
	}
	else{
		line = "hexa8";
		line.append(75, ' ');
	}
	pp->ierr = MPI_File_write_at(ifile, pp->disp, const_cast<char *>(line.data()), 80, MPI_CHARACTER, &pp->status);
	pp->MPI_ierr_Handler(2);
}
pp->disp_header += nnodes*3*4+80;
nelem_l = 0;
size_t elem_S1, elem_S2;
if(pgv->ijkm_gl[2]==1){
	if(pp->ierr==0){
		elem = new int*[4];
		for(auto iiii=0; iiii<4; iiii++){
			elem[iiii] = new int[nnodes_l];
		}
		elem_S1 = 4;
		elem_S2 = nnodes_l;
	}
	else pp->litError("dumpEnsight","allocation error for elem");
}
else{
	if(pp->ierr==0){
		elem = new int*[8];
		for(auto iiii=0; iiii<8; iiii++) elem[iiii] = new int[nnodes_l];
		elem_S1 = 8;
		elem_S2 = nnodes_l;
	}
	else pp->litError("dumpEnsight","allocation error for elem");
}
if(pgv->cylindrical){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		for(auto ic=0; ic<pgv->b->NinZ; ic++){
			int i = pgv->Gn[ic].ijk[0] - pgv->b->imino_;
			int j = pgv->Gn[ic].ijk[1] - pgv->b->jmino_;
			int k = pgv->Gn[ic].ijk[2] - pgv->b->kmino_;
			if(i<pgv->b->imin1_ - pgv->b->imino_) continue;
			if(j<pgv->b->jmin1_ - pgv->b->jmino_) continue;
			if(k<pgv->b->kmin1_ - pgv->b->kmino_) continue;
			if(i+pgv->b->imino_<1 && pgv->sg_periodic[0]) continue;
			if(j+pgv->b->jmino_<1 && pgv->sg_periodic[1]) continue;
			if(i==pgv->ijkm_gl[0]-pgv->b->imino_ && pgv->sg_periodic[0]) continue;
			if(j==pgv->ijkm_gl[1]-pgv->b->jmino_ && pgv->sg_periodic[1]) continue;
			if(k==pgv->ijkm_gl[2]-pgv->b->kmino_ && pgv->sg_periodic[2]) continue;
			if(i>pgv->b->imax_-pgv->b->imino_) continue;
			if(j>pgv->b->jmax_-pgv->b->jmino_) continue;
			if(k>pgv->b->kmax_-pgv->b->kmino_) continue;
			if(pgv->i2c[i][j][k] <= 0) continue;
			if(pgv->i2c[i+1][j][k] <= 0) continue;
			if(pgv->i2c[i+1][j+1][k] <= 0) continue;
			if(pgv->i2c[i][j+1][k] <= 0) continue;
			if(pgv->i2c[i][j][k+1] <= 0) continue;
			if(pgv->i2c[i+1][j][k+1] <= 0) continue;
			if(pgv->i2c[i+1][j+1][k+1] <= 0) continue;
			if(pgv->i2c[i][j+1][k+1] <= 0) continue;
			if(nodeIc[pgv->ioff+pgv->i2c[i][j][k]-1] <= 0) continue;
			if(nodeIc[pgv->ioff+pgv->i2c[i+1][j][k]-1]<=0) continue;
			if(nodeIc[pgv->ioff+pgv->i2c[i+1][j+1][k]-1]<=0) continue;
			if(nodeIc[pgv->ioff+pgv->i2c[i][j+1][k]-1] <= 0) continue;
			if(nodeIc[pgv->ioff+pgv->i2c[i][j][k+1]-1]<=0) continue;
			if(nodeIc[pgv->ioff+pgv->i2c[i+1][j][k+1]-1]<=0) continue;
			if(nodeIc[pgv->ioff+pgv->i2c[i+1][j+1][k+1]-1]<=0) continue;
			if(nodeIc[pgv->ioff+pgv->i2c[i][j+1][k+1]-1] <= 0) continue;
			nelem_l += 1;
			elem[0][nelem_l-1] = nodeIc[pgv->ioff+pgv->i2c[i][j][k]-1];
			elem[1][nelem_l-1] = nodeIc[pgv->ioff+pgv->i2c[i][j][k+1]-1];
			elem[2][nelem_l-1] = nodeIc[pgv->ioff+pgv->i2c[i][j+1][k+1]-1];
			elem[3][nelem_l-1] = nodeIc[pgv->ioff+pgv->i2c[i][j+1][k]-1];
			elem[4][nelem_l-1] = nodeIc[pgv->ioff+pgv->i2c[i+1][j][k]-1];
			elem[5][nelem_l-1] = nodeIc[pgv->ioff+pgv->i2c[i+1][j][k+1]-1];
			elem[6][nelem_l-1] = nodeIc[pgv->ioff+pgv->i2c[i+1][j+1][k+1]-1];
			elem[7][nelem_l-1] = nodeIc[pgv->ioff+pgv->i2c[i+1][j+1][k]-1];
		}
	}
}
else{
	if(pgv->ijkm_gl[2]==1){
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinZ; ic++){
				if(pgv->Gn[ic].ijk[2] != 1) continue;
				int i = pgv->Gn[ic].ijk[0] - pgv->b->imino_;
				int j = pgv->Gn[ic].ijk[1] - pgv->b->jmino_;
				if(i<pgv->b->imin1_-pgv->b->imino_) continue;
				if(j<pgv->b->jmin1_-pgv->b->jmino_) continue;
				if(i==pgv->ijkm_gl[0]-pgv->b->imino_ && pgv->sg_periodic[0]) continue;
				if(j==pgv->ijkm_gl[1]-pgv->b->jmino_ && pgv->sg_periodic[1]) continue;
				if(i>pgv->b->imax_-pgv->b->imino_) continue;
				if(j>pgv->b->jmax_-pgv->b->jmino_) continue;
				if(pgv->i2c[i][j][3] <= 0) continue;
				if(pgv->i2c[i+1][j][3] <= 0) continue;
				if(pgv->i2c[i+1][j+1][3] <= 0) continue;
				if(pgv->i2c[i][j+1][3] <= 0) continue;
				if(nodeIc[pgv->ioff+pgv->i2c[i][j][3]-1] <= 0) continue;
				if(nodeIc[pgv->ioff+pgv->i2c[i+1][j][3]-1] <= 0) continue;
				if(nodeIc[pgv->ioff+pgv->i2c[i+1][j+1][3]-1] <= 0) continue;
				if(nodeIc[pgv->ioff+pgv->i2c[i][j+1][3]-1] <= 0) continue;
				nelem_l += 1;

				elem[0][nelem_l-1] = nodeIc[pgv->ioff+pgv->i2c[i][j][1-pgv->b->kmino_]-1];
				elem[1][nelem_l-1] = nodeIc[pgv->ioff+pgv->i2c[i+1][j][1-pgv->b->kmino_]-1];
				elem[2][nelem_l-1] = nodeIc[pgv->ioff+pgv->i2c[i+1][j+1][1-pgv->b->kmino_]-1];
				elem[3][nelem_l-1] = nodeIc[pgv->ioff+pgv->i2c[i][j+1][1-pgv->b->kmino_]-1];

		}
		}
	}
	else{
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinZ; ic++){
				int i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
				int j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
				int k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
				if(i<pgv->b->imin1_-pgv->b->imino_) continue;	
				if(j<pgv->b->jmin1_-pgv->b->jmino_) continue;
				if(k<pgv->b->kmin1_-pgv->b->kmino_) continue;
				if(i<1-pgv->b->imino_ && pgv->sg_periodic[0]) continue;
				if(j<1-pgv->b->jmino_ && pgv->sg_periodic[1]) continue;
				if(k<1-pgv->b->kmino_ && pgv->sg_periodic[2]) continue;
				if(i==pgv->ijkm_gl[0]-pgv->b->imino_ && pgv->sg_periodic[0]) continue;
				if(j==pgv->ijkm_gl[1]-pgv->b->jmino_ && pgv->sg_periodic[1]) continue;
				if(k==pgv->ijkm_gl[2]-pgv->b->kmino_ && pgv->sg_periodic[2]) continue;
				if(i>pgv->b->imax_-pgv->b->imino_) continue;
				if(j>pgv->b->jmax_-pgv->b->jmino_) continue;
				if(k>pgv->b->kmax_-pgv->b->kmino_) continue;
				if(pgv->i2c[i][j][k] <= 0) continue;
				if(pgv->i2c[i+1][j][k] <= 0) continue;				
				if(pgv->i2c[i+1][j+1][k] <= 0) continue;
				if(pgv->i2c[i][j+1][k] <= 0) continue;
				if(pgv->i2c[i][j][k+1] <= 0) continue;
				if(pgv->i2c[i+1][j][k+1] <= 0) continue;
				if(pgv->i2c[i+1][j+1][k+1] <= 0) continue;
				if(pgv->i2c[i][j+1][k+1] <= 0) continue;
				if(nodeIc[pgv->ioff+pgv->i2c[i][j][k]-1]<=0)continue;
				if(nodeIc[pgv->ioff+pgv->i2c[i+1][j][k]-1]<=0)continue;
				if(nodeIc[pgv->ioff+pgv->i2c[i+1][j+1][k]-1]<=0)continue;
				if(nodeIc[pgv->ioff+pgv->i2c[i][j+1][k]-1]<=0)continue;
				if(nodeIc[pgv->ioff+pgv->i2c[i][j][k+1]-1]<=0)continue;
				if(nodeIc[pgv->ioff+pgv->i2c[i+1][j][k+1]-1]<=0)continue;
				if(nodeIc[pgv->ioff+pgv->i2c[i+1][j+1][k+1]-1]<=0)continue;
				if(nodeIc[pgv->ioff+pgv->i2c[i][j+1][k+1]-1]<=0)continue;
				nelem_l += 1;
				elem[0][nelem-1] = nodeIc[pgv->ioff+pgv->i2c[i][j][k]-1];
				elem[1][nelem-1] = nodeIc[pgv->ioff+pgv->i2c[i+1][j][k]-1];
				elem[2][nelem-1] = nodeIc[pgv->ioff+pgv->i2c[i+1][j+1][k]-1];
				elem[3][nelem-1] = nodeIc[pgv->ioff+pgv->i2c[i][j+1][k]-1];
				elem[4][nelem-1] = nodeIc[pgv->ioff+pgv->i2c[i][j][k+1]-1];
				elem[5][nelem-1] = nodeIc[pgv->ioff+pgv->i2c[i+1][j][k+1]-1];
				elem[6][nelem-1] = nodeIc[pgv->ioff+pgv->i2c[i+1][j+1][k+1]-1];
				elem[7][nelem-1] = nodeIc[pgv->ioff+pgv->i2c[i][j+1][k+1]-1];
			}
		}
	}
}

pp->parallel_sum(nelem_l, nelem);
std::cout<<nelem_l<<std::endl;


if(pp->myrank == 0){
	pp->disp = pp->disp_header;
	pp->ierr=MPI_File_write_at(ifile,pp->disp,&nelem,1,MPI_INTEGER,&pp->status);
	pp->MPI_ierr_Handler(3);
}
pp->disp_header += 4;
//pp->parallel_all_gather(nelem_l, nhelp);
if(pp->myrank == 0) ielem_offset = 0;
else ielem_offset = nelem_l; //std::accumulate(nhelp.begin(), nhelp.begin()+pp->myrank, 0);
;
if(pgv->ijkm_gl[2]==1){
	pp->disp = pp->disp_header+ielem_offset*4*4;
	pp->ierr=MPI_File_write_at(ifile, pp->disp, elem, 4*nelem_l,MPI_INTEGER,&pp->status);
}
else{
	pp->disp = pp->disp_header+ielem_offset*8*4;
	pp->ierr=MPI_File_write_at(ifile,pp->disp, elem, 8*nelem_l, MPI_INTEGER,&pp->status);
}
pp->MPI_ierr_Handler(4);
pp->ierr=MPI_File_close(&ifile);
if(pp->myrank==0){
	for(auto i=0; i<elem_S1; i++) 
		delete [] elem[i];
	delete [] elem;
}
else pp->litError("dumpEnsight","deallocation error for elem");
if(pp->ierr==0) coord = new double[nnodes_l];
else pp->litError("dumpEnsight","allocation error for coord");
for(auto iv=5; iv<nVarDump+5; iv++){ // need to be read correctly in input file
	if(ppar->trim(varDumpName[iv])=="G"){
		writeEnsightHeaderScalar("G",filename,name,ifile,line,chelp6,inth,fileExists, "node");
		inode = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinZ; ic++){
				if(nodeFlag[ic]){
					inode += 1;
					coord[inode-1] = static_cast<double>(pgv->Gn[ic-pgv->ibl2buf[ibl]].G);
				}
			}
		}
		writeEnsightScalar(nnodes_l, coord, ifile, inode_offset);


	}
	if(ppar->trim(varDumpName[iv])=="VOF"){
		writeEnsightHeaderScalar("VOF",filename,name,ifile,line,chelp6,inth,fileExists, "node");
		inode = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinZ; ic++){
				if(nodeFlag[ic]){
					inode += 1;
					coord[inode-1]=ptool->cell_volume_fraction(pgv->Gn[ic-pgv->ioff],ic-pgv->ioff+1,pgv->b,2);
				}
			}
		}
		writeEnsightScalar(nnodes_l, coord, ifile, inode_offset);
	}
	else if(ppar->trim(varDumpName[iv])=="S"){
		writeEnsightHeaderScalar("S",filename,name,ifile,line,chelp6,inth,fileExists, "node");
		inode = 0;
		size_t S_S;
		S = plbuf->getR1Buffer(S_S);
		prein->calcSignFunction(S);
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinZ; ic++){
				if(nodeFlag[ic]){
					inode += 1;
					coord[inode-1] = S[ic];
				}
			}
		}
		plbuf->freeR1Buffer(S);
		writeEnsightScalar(nnodes_l, coord, ifile, inode_offset);
	}
	else if(ppar->trim(varDumpName[iv])=="PE"){
		writeEnsightHeaderScalar("PE",filename,name,ifile,line,chelp6,inth,fileExists, "node");
		inode = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinZ; ic++){
				if(nodeFlag[ic]){
					inode += 1;
					coord[inode-1] = static_cast<double>(pp->myrank+1);
				}
			}
		}
		writeEnsightScalar(nnodes_l,coord, ifile, inode_offset);
	}
	else if(ppar->trim(varDumpName[iv])=="LAYER"){
		writeEnsightHeaderScalar("LAYER",filename,name,ifile,line,chelp6,inth,fileExists, "node");
		inode = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			auto ic = pgv->ioff;
			auto temp_co = std::accumulate(pgv->band_size.begin(), pgv->band_size.begin()+which_Band, 0);
			for(auto ib=0; ib<temp_co; ib++){
				for(auto il=0; il<pgv->b->NinBandLayer[ib]; il++){
					ic += 1;
					if(nodeFlag[ic-1]){
						inode += 1;
						coord[inode] = static_cast<double>(ib);
					}
				}
			}
		}
		writeEnsightScalar(nnodes_l, coord, ifile, inode_offset);
	}
	else if(ppar->trim(varDumpName[iv])=="BAND_ID"){
		writeEnsightHeaderScalar("BAND_ID",filename,name,ifile,line,chelp6,inth,fileExists, "node");
		inode = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinA; ic++){
				if(nodeFlag[pgv->ioff+ic]){
					inode += 1;
					coord[inode-1] = 1.0;
				}
			}
			for(auto ic=pgv->b->NinA; ic<pgv->b->NinT; ic++){
				if(nodeFlag[pgv->ioff+ic]){
					inode += 1;
					coord[inode-1] = 2.0;
				}
			}
			for(auto ic=pgv->b->NinT; ic<pgv->b->NinN; ic++){
				if(nodeFlag[pgv->ioff+ic]){
					inode += 1;
					coord[inode-1] = 3.0;
				}
			}
			for(auto ic=pgv->b->NinN; ic<pgv->b->NinW; ic++){
				if(nodeFlag[pgv->ioff+ic]){
					inode += 1;
					coord[inode-1] = 4.0;
				}
			}
			for(auto ic=pgv->b->NinW; ic<pgv->b->NinX; ic++){
				if(nodeFlag[pgv->ioff+ic]){
					inode += 1;
					coord[inode-1] = 5.0;
				}
			}
			for(auto ic=pgv->b->NinX; ic<pgv->b->NinZ; ic++){
				if(nodeFlag[pgv->ioff+ic]){
					inode += 1;
					coord[inode-1] = 6.0;

				}
			}
		}
		writeEnsightScalar(nnodes_l, coord, ifile, inode_offset);
	}
	else if(ppar->trim(varDumpName[iv])=="KAPPA"){
		writeEnsightScalar(nnodes_l, coord, ifile, inode_offset);
		size_t kappa_S;
		kappa  = plbuf->getR1Buffer(kappa_S);
		size_t marker_S;
		marker = plbuf->getR1Buffer(marker_S);
		ptool->calcCurvature(kappa_S,kappa,marker_S,marker);
		inode = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinZ; ic++){
				if(nodeFlag[ic]){
					inode += 1;
					coord[inode-1] = static_cast<double>(kappa[ic]*marker[ic]);
				}
			}
		}
		writeEnsightScalar(nnodes_l, coord, ifile, inode_offset);
		plbuf->freeR1Buffer(kappa);
		plbuf->freeR1Buffer(marker);
	}
	else if(ppar->trim(varDumpName[iv]) == "DIST_DIJK"){
		writeEnsightHeaderScalar("DIST_DIJK",filename,name,ifile,line,chelp6,inth,fileExists, "node");
		inode = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinZ; ic++){
				if(nodeFlag[ic]){
					inode += 1;
					auto temp_sum=pow(pgv->Gn[ic-pgv->ibl2buf[ibl]].dijk[0],2)+\
									  pow(pgv->Gn[ic-pgv->ibl2buf[ibl]].dijk[1],2)+\
									  pow(pgv->Gn[ic-pgv->ibl2buf[ibl]].dijk[2],2);
					coord[inode-1] = static_cast<double>(sqrt(static_cast<double>(temp_sum)));
				}
			}
		}
		writeEnsightScalar(nnodes_l, coord, ifile, inode_offset);
	}
	else if(ppar->trim(varDumpName[iv]) == "BLOCK"){
		writeEnsightHeaderScalar("BLOCK",filename,name,ifile,line,chelp6,inth,fileExists, "node");
		inode = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinZ; ic++){
				if(nodeFlag[ic]){
					inode += 1;
					coord[inode-1] = static_cast<double>(ibl);
				}
			}
		}
		writeEnsightScalar(nnodes_l, coord, ifile, inode_offset);
	}
	else if(ppar->trim(varDumpName[iv])=="V"){
		writeEnsightHeaderVector("V",filename,name,ifile,line,chelp6,inth,fileExists, "node");
		if(pgv->cylindrical){
			for(auto id=0; id<3; id++){
				inode = 0;
				for(auto ibl=0; ibl<pgv->nbl; ibl++){
					pgv->setBlockPointers(ibl);
					for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinZ; ibl++){
						if(nodeFlag[ic]){
							inode += 1;
							theta = pgv->zc[pgv->Gn[ic-pgv->ioff].ijk[2]-(1-pgv->nghost)];
							if(id == 1){
								coord[inode-1]=static_cast<double>(pgv->Gn[ic-pgv->ioff].V[1]*\
													cos(theta)-pgv->Gn[ic-pgv->ioff].V[2]*sin(theta));
							}
							else if(id == 2){
								coord[inode-1]=static_cast<double>(pgv->Gn[ic-pgv->ioff].V[1]*\
													sin(theta)+pgv->Gn[ic-pgv->ioff].V[2]*cos(theta));
							}
							else{
								coord[inode-1]=static_cast<double>(pgv->Gn[ic-pgv->ioff].V[id]);
							}
						}
					}
				}
				writeEnsightVector(nnodes_l, coord, id, ifile, inode_offset, nnodes);
			}
		}
		else{
			for(auto id=0; id<3; id++){
				inode = 0;
				for(auto ibl=0; ibl<pgv->nbl; ibl++){
					pgv->setBlockPointers(ibl);
					for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinZ; ic++){
						if(nodeFlag[ic]){
							inode += 1;
							coord[inode-1]=static_cast<double>(pgv->Gn[ic-pgv->ioff].V[id]);
						}
					}
				}
				writeEnsightVector(nnodes_l, coord, id, ifile, inode_offset, nnodes);
			}
		}
	}
	else if(ppar->trim(varDumpName[iv])=="DIJK"){
		writeEnsightHeaderVector("DIJK",filename,name,ifile,line,chelp6,inth,fileExists,"node");
		if(pgv->cylindrical){
			for(auto id=0; id<3; id++){
				inode = 0;
				for(auto ibl=0; ibl<pgv->nbl; ibl++){
					pgv->setBlockPointers(ibl);
					for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinZ; ic++){
						if(nodeFlag[ic]){
							inode += 1;
							theta = pgv->zc[pgv->Gn[ic-pgv->ioff].ijk[2]-(1-pgv->nghost)];
							if(id == 1){
								coord[inode-1]=static_cast<double>(pgv->Gn[ic-pgv->ioff].dijk[1]*cos(theta)-pgv->Gn[ic-pgv->ioff].dijk[2]*sin(theta));
							}
							else if(id==2){
								coord[inode-1]=static_cast<double>(pgv->Gn[ic-pgv->ioff].dijk[1]*sin(theta)+pgv->Gn[ic-pgv->ioff].dijk[2]*cos(theta));
							}
							else{
								coord[inode-1]=static_cast<double>(pgv->Gn[ic-pgv->ioff].dijk[id]);
							}
						}
					}
				}
				writeEnsightVector(nnodes_l,coord, id, ifile, inode_offset, nnodes);
			}
		}
		else{
			for(auto id=0; id<3; id++){
				inode = 0;
				for(auto ibl=0; ibl<pgv->nbl; ibl++){
					pgv->setBlockPointers(ibl);
					for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinZ; ic++){
						if(nodeFlag[ic]){
							inode += 1;
							coord[inode-1]=static_cast<double>(pgv->Gn[ic-pgv->ioff].dijk[id]);
						}
					}
					writeEnsightVector(nnodes_l,coord, id, ifile, inode_offset, nnodes);
				}
			}
		}
	}
}
if(pp->ierr==0) delete [] coord;
else pp->litError("dumpEnsight","deallocation error, for coord");
plbuf->freeL1Buffer(nodeFlag);
plbuf->freeI1Buffer(nodeIc);
if(pgv->verbose && pp->myrank==0){
	std::cout<<pgv->clit<<"Starting dump_Gnodes_Ensight of file "<<ppar->trim(name)<<"... Done"<<std::endl;
}
pti->lit_timing_stop("io_ensight");
};

void io::writeEnsightScalar(size_t var_S, double *var, MPI_File & ifile, int &inode_offset){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<bl> pbl(new bl);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<band> pba(new band);
pp->disp = pp->disp_header + static_cast<int>(inode_offset)*4;
pp->ierr=MPI_File_write_at(ifile,pp->disp,&var,var_S,MPI_DOUBLE,&pp->status);
pp->MPI_ierr_Handler(5);
pp->ierr = MPI_File_close(&ifile);
pp->MPI_ierr_Handler(51);
};

void io::writeEnsightVector(size_t var_S, double *var, int id, MPI_File &ifile,\
									 int &inode_offset, int &nnodes){

std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<bl> pbl(new bl);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<band> pba(new band);

pp->disp=pp->disp_header+static_cast<int>(inode_offset)*4+(static_cast<int>(id)-1)*\
																				 static_cast<int>(nnodes)*4;
pp->ierr=MPI_File_write_at(ifile,pp->disp,&var,var_S,MPI_DOUBLE,&pp->status);

pp->MPI_ierr_Handler(5);
if(id == 2){
	pp->ierr = MPI_File_close(&ifile);
	pp->MPI_ierr_Handler(51);
}
};

void io::writeEnsightHeaderScalar(const std::string varname, std::string &filename, \
											 const std::string name, MPI_File &ifile,std::string &line, \
											 std::string &chelp6,int &inth,bool &fileExists,\
											 const std::string type){

std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<bl> pbl(new bl);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<band> pba(new band);
std::stringstream sss;
sss<<std::setw(6)<<std::setfill('0')<<chelp6;
std::string ssss = sss.str();
if(type=="node"){
	filename = ppar->trim(pgv->dump_dir)+"sol_"+ppar->trim(name)+"_"+ssss+"_"+ppar->trim(varname)+".scalar";
}
else if(type == "elem"){
	filename=ppar->trim(pgv->dump_dir)+"sol_cell_"+ppar->trim(name)+"_"+ssss+"_"+ppar->trim(varname)+".scalar";
}
if(pp->myrank == 0){
	fileExists = exists_file(ppar->trim(filename));
	if(fileExists){
		if(pgv->verbose) std::cout<<pgv->clit<<"deleting"<<'\t'<<ppar->trim(varname)<<".scalar"<<std::endl;
		pp->ierr=MPI_File_delete(const_cast<char *>(filename.data()), pp->MPI_IO_LIT_FILE_HINT);
	}
	if(pgv->verbose) std::cout<<pgv->clit<<"writing"<<'\t'<<ppar->trim(varname)<<".scalar"<<std::endl;
}
pp->parallel_barrier();
//filename = ppar->trim(pp->file_prefix)+ppar->trim(filename);
pp->ierr=MPI_File_open(MPI_COMM_WORLD,const_cast<char*>(filename.data()),\
							  MPI_MODE_WRONLY+MPI_MODE_CREATE, pp->MPI_IO_LIT_FILE_HINT, &ifile);
std::cout<<pp->ierr<<std::endl;
if(pp->myrank == 0){
	line = ppar->trim(varname);
	line.append(80-line.size(), ' ');
	pp->ierr = MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
	line = "part";
	line.append(80-line.size(), ' ');
	pp->ierr = MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
	inth = 1;
	pp->ierr = MPI_File_write(ifile,&inth,1,MPI_INTEGER,&pp->status);
	if(type=="node"){
		line = "coordinates";
		line.append(80-line.size(), ' ');
	}
	else if(type=="elem"){
		line=="hexa8";
		line.append(80-line.size(), ' ');
	}
	pp->ierr = MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
}
pp->disp_header = 3*80*1+1*4;
};

void io::writeEnsightHeaderVector(const std::string varname, std::string &filename, \
											 const std::string name, MPI_File &ifile,std::string &line, \
											 std::string &chelp6,int &inth,bool &fileExists,\
											 const std::string type){

std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<bl> pbl(new bl);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<band> pba(new band);
std::stringstream sss;
sss<<std::setw(6)<<std::setfill('0')<<chelp6;
std::string ssss = sss.str();
if(type=="node"){
	filename=ppar->trim(pgv->dump_dir)+"sol_"+ppar->trim(name)+"_"+ssss+"_"+ppar->trim(varname)+".vector";
}
else if(type=="elem"){
	filename=ppar->trim(pgv->dump_dir)+"sol_cell_"+ppar->trim(name)+"_"+ssss+"_"+ppar->trim(varname)+".vector";
}
if(pp->myrank == 0){
	fileExists = exists_file(ppar->trim(filename));
	if(fileExists){
		if(pgv->verbose) std::cout<<pgv->clit<<"deleting"<<'\t'<<ppar->trim(varname)<<".vector"<<std::endl;
		pp->ierr=MPI_File_delete(const_cast<char*>(filename.data()), pp->MPI_IO_LIT_FILE_HINT);
	}
	if(pgv->verbose) std::cout<<pgv->clit<<"writing"<<'\t'<<ppar->trim(varname)<<".vector"<<std::endl;
}
//filename = ppar->trim(pp->file_prefix)+ppar->trim(filename);
MPI_File_open(MPI_COMM_WORLD,const_cast<char*>(filename.data()),MPI_MODE_WRONLY+MPI_MODE_CREATE, \
				  pp->MPI_IO_LIT_FILE_HINT,&ifile);
if(pp->myrank == 0){
	line = ppar->trim(varname);
	pp->ierr=MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
	line = "part";
	pp->ierr=MPI_File_write(ifile,const_cast<char*>(line.data()), 80, MPI_CHARACTER,&pp->status);
	inth = 1;
	pp->ierr=MPI_File_write(ifile,&inth, 1, MPI_INTEGER,&pp->status);
	if(type=="node") line = "coordinates";
	else if(type=="elem") line="hexa8";
	MPI_File_write(ifile, const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
}
pp->disp_header = 3*80*1+1*4;
};

void io::dumpEnsightCell(const std::string name, const std::string band1, bool newv, int arg){

std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<bl> pbl(new bl);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<band> pba(new band);

int inth, nCorners_l, nCorners, iCorner, iCorner_offset,nelem,nelem_l;
int ielem, which_band, ielem_offset, info, strlen;
std::vector<int> nhelp(pp->nprocs);
MPI_File ifile;
bool exists, append, fileExists;
std::string line(80, ' '); // was null
std::string errstr(1024, ' '), filename(1024, ' '),fname2(1024, ' '); // was null
std::string chelp6(6, ' '); // was null
std::vector<double> coord;
std::vector<int> elem, cMarker;
std::vector<bool> lMarker;
double *kappa, *marker, *S;
pti->lit_timing_start("io_ensight");
if(pgv->verbose && pp->myrank == 0){
	std::cout<<pgv->clit<<"Starting dumpEnsightCell of file "<<ppar->trim(name)<<std::endl;
}
if(ppar->trim(band1) == "A") which_band = 1;
else if(ppar->trim(band1)=="T") which_band = 2;
else if(ppar->trim(band1)=="N") which_band = 3;
else if(ppar->trim(band1)=="W") which_band = 4;
else if(ppar->trim(band1)=="X") which_band = 5;
else{
	std::cout<<pgv->clit<<"unknown band name for dumpEnsightCell: "<<ppar->trim(band1)<<std::endl;
	pp->parallel_kill(1);
}
if(pp->myrank == 0){
	if(arg==4){
		if(newv){
			write_new_case_file(name,2);
			append = false;
		}
		else append = true;
	}
	else append = true;
	if(append){
		exists=exists_file(ppar->trim(pgv->dump_dir)+ppar->trim(name)+".case");
		if(!exists) write_new_case_file(name,2);
		else append_case_file(name,2);
	}
}
pp->ierr = MPI_Bcast(&dump_nfiles, 1, MPI_INTEGER,0,MPI_COMM_WORLD);
#ifdef WITH_CDP
if(pgv->dump_dtime_fs > 0.0){
	chelp6 = std::to_string(dump_nfiles-1); // correct format
}
else{
	chelp6 = std::to_string(pgv->fs_step); // correct format
}
#else
	chelp6 = std::to_string(dump_nfiles-1); // correct format
#endif
std::stringstream sss;
sss<<std::setw(6)<<std::setfill('0')<<chelp6;
std::string ssss = sss.str();
filename=ppar->trim(pgv->dump_dir)+"sol_cell_"+ppar->trim(name)+"_"+ssss+".geo";
if(pp->myrank==0){
	fileExists = exists_file(ppar->trim(filename));
	if(fileExists){
		if(pgv->verbose) std::cout<<pgv->clit<<"deleting geo"<<std::endl;
		pp->ierr=MPI_File_delete(const_cast<char*>(filename.data()),pp->MPI_IO_LIT_FILE_HINT);
	}
	if(pgv->verbose) std::cout<<pgv->clit<<"writing geo"<<ppar->trim(filename)<<std::endl;
}
pp->parallel_barrier();
filename=ppar->trim(pp->file_prefix)+ppar->trim(filename);
pp->ierr=MPI_File_open(MPI_COMM_WORLD,const_cast<char*>(filename.data()), \
			MPI_MODE_WRONLY+MPI_MODE_CREATE,pp->MPI_IO_LIT_FILE_HINT, &ifile);
if(pp->myrank==0){
	line = "C Binary ";
	pp->ierr=MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
	line = ppar->trim(band1)+"-Band nodes";
	pp->ierr=MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
	line= " ";
	pp->ierr=MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
	line = "node id off";
	pp->ierr=MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
	line = "element id off";
	pp->ierr=MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
	line = "part";
	pp->ierr=MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
	inth = 1;
	pp->ierr==MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_INTEGER,&pp->status);
	line = ppar->trim(band1)+"-Bnad coordinates";
	pp->ierr=MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
	line = "coordinates";
	pp->ierr=MPI_File_write(ifile,const_cast<char*>(line.data()),80,MPI_CHARACTER,&pp->status);
}
pp->disp_header = 8*80*1 +1*4;
auto s1 = pgv->ijkm_bl[0]+1;
auto s2 = pgv->ijkm_bl[1]+1;
auto s3 = pgv->ijkm_bl[2]+1;
if(pp->ierr==0){
	cMarker.resize(s1*s2*s3*pgv->nbl,-1);
	lMarker.resize(s1*s2*s3*pgv->nbl,false);
}
else pp->litError("dumpEnsightCell","cannot allocate cMarker, lMarker");
//this is a different algorithm from fortran; line 1899; // checked but 4D array is fine
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinBand[which_band-1]; ic++){
		auto i = pgv->Gn[ic].ijk[0]-pgv->b->ijk0[0]-1;
		auto j = pgv->Gn[ic].ijk[1]-pgv->b->ijk0[1]-1;
		auto k = pgv->Gn[ic].ijk[2]-pgv->b->ijk0[2]-1;
		for(auto ii=i; ii<=i+1; ii++){
			for(auto jj=j; jj<=j+1; jj++){
				for(auto kk=k; kk<=k+1; kk++){
					lMarker[ii+s2*(jj+s3*(kk+pgv->nbl*ibl))]=true; // 4d array is fine; checked
				}
			}
		}
	}
}
nCorners_l = 0;
for(auto i:lMarker){
	if(i==true) nCorners_l++;

pp->parallel_all_sum(nCorners_l, nCorners);
if(pp->myrank == 0) 
	iCorner_offset = std::accumulate(nhelp.begin(),nhelp.begin()+pp->myrank, 0);
else iCorner_offset = 0;
if(pp->ierr==0) coord.resize(nCorners_l);
else pp->litError("dumpEnsightCell","allocation error for coord");
for(auto id=0; id<3; id++){
	iCorner = 0;
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto i=pgv->b->imin_-1; i<pgv->b->imax_+1; i++){
			for(auto j=pgv->b->jmin_-1; j<pgv->b->jmax_+1; j++){
				for(auto k=pgv->b->kmin_-1; k<pgv->b->kmax_+1; k++){
					auto t1 = i-pgv->b->ijk0[0];
					auto t2 = j-pgv->b->ijk0[1];
					auto t3 = k-pgv->b->ijk0[2];
					if(lMarker[t1+s2*(t2+s3*(t3+pgv->nbl*ibl))]){
						iCorner += 1;
						cMarker[t1+s2*(t2+s3*(t3+pgv->nbl*ibl))] = iCorner+iCorner_offset;
						switch(id){
							case 1:
								coord[iCorner-1] = pgv->lxf[i+pgv->nghost];
								break;
							case 2:
								coord[iCorner-1] = pgv->lyf[j+pgv->nghost];
								break;
							case 3:
								coord[iCorner-1] = pgv->lzf[k+pgv->nghost];
								break;
						}
					}
				}
			}
		}
	}
	pp->disp = pp->disp_header+(id*nCorners+iCorner_offset)*4;
	pp->ierr==MPI_File_write_at(ifile,pp->disp,&coord,nCorners_l,MPI_DOUBLE,&pp->status);
	pp->MPI_ierr_Handler(1);
}
if(pp->ierr==0){
	coord.clear();
	lMarker.clear();
}
else pp->litError("dumpEnsightCell","cannot deallocate coord_lMarker");
if(pp->myrank==0){
	pp->disp = pp->disp_header+nCorners*3*4;
	line = "hexa8";
	pp->ierr=MPI_File_write_at(ifile,pp->disp,const_cast<char*>(line.data()),\
				80,MPI_CHARACTER,&pp->status);
	pp->MPI_ierr_Handler(2);
}
pp->disp_header += nCorners*3*4+80;
nelem_l = 0;
if(pp->ierr==0) elem.resize(8*nCorners_l);
else pp->litError("dumpEnsightCell","cannot allocate elem");
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinBand[which_band-1]; ic++){
		auto i = pgv->Gn[ic].ijk[0]-pgv->b->ijk0[0]-1;
		auto j = pgv->Gn[ic].ijk[1]-pgv->b->ijk0[1]-1;
		auto k = pgv->Gn[ic].ijk[2]-pgv->b->ijk0[2]-1;
		nelem_l += 1;
		elem[0+nCorners_l*(nelem-1)] = cMarker[i+s2*(j+s3*(k+pgv->nbl*ibl))];
		elem[1+nCorners_l*(nelem-1)] = cMarker[(i+1)+s2*(j+s3*(k+pgv->nbl*ibl))];
		elem[2+nCorners_l*(nelem-1)] = cMarker[(i+1)+s2*((j+1)+s3*(k+pgv->nbl*ibl))];
		elem[3+nCorners_l*(nelem-1)] = cMarker[i+s2*((j+1)+s3*(k+pgv->nbl*ibl))];
		elem[4+nCorners_l*(nelem-1)] = cMarker[i+s2*(j+s3*((k+1)+pgv->nbl*ibl))];
		elem[5+nCorners_l*(nelem-1)] = cMarker[(i+1)+s2*(j+s3*((k+1)+pgv->nbl*ibl))];
		elem[6+nCorners_l*(nelem-1)] = cMarker[(i+1)+s2*((j+1)+s3*((k+1)+pgv->nbl*ibl))];
		elem[7+nCorners_l*(nelem-1)] = cMarker[i+s2*((j+1)+s3*((k+1)+pgv->nbl*ibl))];
	}
}
if(pp->ierr==0) cMarker.clear();
else pp->litError("dumpEnsightCell","cannot deallocate cMarker");
pp->parallel_all_sum(nelem_l, nelem);
if(pp->myrank == 0){
	pp->disp = pp->disp_header;
	pp->ierr=MPI_File_write_at(ifile,pp->disp,&nelem,1,MPI_INTEGER,&pp->status);
	pp->MPI_ierr_Handler(3);
}
pp->disp_header += 4;
pp->parallel_all_gather(nelem_l,nhelp);
if(pp->myrank > 0) ielem_offset = std::accumulate(nhelp.begin(),nhelp.begin()+pp->myrank, 0);
else ielem_offset = 0;
pp->disp = pp->disp_header + ielem_offset*8*4;
pp->ierr==MPI_File_write_at(ifile,pp->disp,&elem,8*nelem_l,MPI_INTEGER,&pp->status);
pp->MPI_ierr_Handler(4);
pp->ierr=MPI_File_close(&ifile);
if(pp->ierr==0) elem.clear();
else pp->litError("dumpEnsightCell","deallocation error for elem");
if(pp->ierr==0) coord.resize(nelem_l);
else pp->litError("dumpEnsight","allocation error for coord");
for(auto iv=0; iv<nVarCellDump; iv++){
	if(ppar->trim(varCellDumpName[iv])== "G"){
		ielem = 0;
		for(auto ibl=0; ibl<pgv->nbl;ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinBand[which_band-1]; ic++){
				ielem += 1;
				coord[ielem-1] = static_cast<double>(pgv->Gn[ic].G);
			}
		}
		writeEnsightScalar(coord.size(), coord.data(), ifile, ielem_offset);
	}
	else if(ppar->trim(varCellDumpName[iv]) == "VOF"){
		writeEnsightHeaderScalar("VOF",filename,name,ifile,line,ssss,inth,fileExists,"elem");
		ielem = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinBand[which_band-1]; ic++){
				ielem += 1;
				coord[ielem-1] = static_cast<double>(ptool->cell_volume_fraction(pgv->Gn[ic],ic+1,pgv->b,3));
			}
		}
		writeEnsightScalar(coord.size(), coord.data(),ifile,ielem_offset);
	}
	else if(ppar->trim(varCellDumpName[iv])=="S"){
		writeEnsightHeaderScalar("S",filename,name,ifile,line,ssss,inth,fileExists,"elem");
		ielem = 0;
		size_t S_S;
		S = plbuf->getR1Buffer(S_S);
		prein->calcSignFunction(S);
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinBand[which_band-1]; ic++){
				ielem += 1;
				coord[ielem-1] = static_cast<double>(S[pgv->ioff+ic]);
			}
		}
		plbuf->freeR1Buffer(S);
		writeEnsightScalar(coord.size(), coord.data(), ifile,ielem_offset);
	}
	else if(ppar->trim(varCellDumpName[iv])=="PE"){
		writeEnsightHeaderScalar("PE",filename,name,ifile,line,ssss,inth,fileExists,"elem");
		ielem = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinBand[which_band-1]; ic++){
				ielem += 1;
				coord[ielem-1] = static_cast<double>(pp->myrank+1);
			}
		}
		writeEnsightScalar(coord.size(), coord.data(), ifile, ielem_offset);
	}
	else if(ppar->trim(varCellDumpName[iv])=="LAYER"){
		writeEnsightHeaderScalar("LAYER",filename,name,ifile,line,ssss,inth,fileExists,"elem");
		ielem = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			auto ul = std::accumulate(pgv->band_size.begin(), pgv->band_size.begin()+which_band, 0);
			for(auto ib=0; ib<ul; ib++){
				for(auto il=0; il<pgv->b->NinBandLayer[ib]; il++){
					ielem += 1;
					coord[ielem-1] = static_cast<double>(ib+1);
				}
			}
		}
		writeEnsightScalar(coord.size(), coord.data(), ifile, ielem_offset);
	}
	else if(ppar->trim(varCellDumpName[iv]) == "BAND_ID"){
		writeEnsightHeaderScalar("BAND_ID",filename,name,ifile,line,ssss,inth,fileExists,"elem");
		ielem = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinA; ic++){
				ielem += 1;
				coord[ielem-1] = 1.0;
			}
			auto ul = std::min(pgv->b->NinBand[which_band-1],pgv->b->NinT);
			for(auto ic=pgv->b->NinA; ic<ul; ic++){
				ielem += 1;
				coord[ielem-1] = 2.0;
			}
			ul = std::min(pgv->b->NinBand[which_band-1],pgv->b->NinN);
			for(auto ic=pgv->b->NinT; ic<ul; ic++){
				ielem += 1;
				coord[ielem-1]  =3.0;
			}
			ul = std::min(pgv->b->NinBand[which_band-1], pgv->b->NinW);
			for(auto ic=pgv->b->NinN; ic<ul; ic++){
				ielem += 1;
				coord[ielem-1] = 4.0;
			}
			ul = std::min(pgv->b->NinBand[which_band-1], pgv->b->NinX);
			for(auto ic=pgv->b->NinW; ic<ul; ic++){
				ielem += 1;
				coord[ielem-1] = 5.0;
			}
			ul = std::min(pgv->b->NinBand[which_band-1], pgv->b->NinZ);
			for(auto ic=pgv->b->NinX; ic<ul; ic++){
				ielem += 1;
				coord[ielem-1] = 6.0;
			}
		}
		writeEnsightScalar(coord.size(), coord.data(), ifile, ielem_offset);
	}
	else if(ppar->trim(varCellDumpName[iv])=="KAPPA"){
		writeEnsightHeaderScalar("KAPPA",filename,name,ifile,line,ssss,inth,fileExists,"elem");
		size_t kappa_S;
		kappa = plbuf->getR1Buffer(kappa_S);
		size_t marker_S;
		marker = plbuf->getR1Buffer(marker_S);
		ptool->calcCurvature(kappa_S, kappa,marker_S, marker);
		ielem = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinBand[which_band-1]; ic++){
				ielem += 1;
				auto sum = pow(pgv->Gn[ic].dijk[0],2)+pow(pgv->Gn[ic].dijk[1],2)+pow(pgv->Gn[ic].dijk[2],2);
				coord[ielem-1] = static_cast<double>(sqrt(sum));
			}
		}
		writeEnsightScalar(coord.size(), coord.data(), ifile, ielem_offset);
	}
	else if(ppar->trim(varCellDumpName[iv])=="DIST_DIJK"){
		writeEnsightHeaderScalar("DIST_DIJK",filename,name,ifile,line,ssss,inth,fileExists,"elem");
		ielem = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinBand[which_band-1]; ic++){
				ielem += 1;
				auto sum = pow(pgv->Gn[ic].dijk[0],2)+pow(pgv->Gn[ic].dijk[1],2)+pow(pgv->Gn[ic].dijk[2],2);
				coord[ielem-1] = static_cast<double>(sqrt(sum));
			}
		}
		writeEnsightScalar(coord.size(), coord.data(), ifile, ielem_offset);
	}
	else if(ppar->trim(varCellDumpName[iv]) == "BLOCK"){
		ielem = 0;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinBand[which_band-1]; ic++){
				ielem += 1;
				coord[ielem-1] = static_cast<double>(ibl);
			}
		}
		writeEnsightScalar(coord.size(), coord.data(), ifile, ielem_offset);
	}
	else if(ppar->trim(varCellDumpName[iv]) == "V"){
		writeEnsightHeaderVector("V",filename,name,ifile,line,ssss,inth,fileExists,"elem");
		for(auto id=0; id<3; id++){
			ielem = 0;
			for(auto ibl=0; ibl<pgv->nbl; ibl++){
				pgv->setBlockPointers(ibl);
				for(auto ic=0; ic<pgv->b->NinBand[which_band-1]; ic++){
					ielem += 1;
					coord[ielem-1] = static_cast<double>(pgv->Gn[ic].V[id]);
				}
			}
			writeEnsightVector(coord.size(), coord.data(), id,ifile,ielem_offset,nelem);
		}
	}
	else if(ppar->trim(varCellDumpName[iv]) == "DIJK"){
		writeEnsightHeaderVector("DIJK",filename,name,ifile,line,ssss,inth,fileExists,"elem");
		for(auto id=0; id<3; id++){
			ielem = 0;
			for(auto ibl=0; ibl<pgv->nbl; ibl++){
				pgv->setBlockPointers(ibl);
				for(auto ic=0; ic<pgv->b->NinBand[which_band-1]; ic++){
					ielem += 1;
					coord[ielem-1] = static_cast<double>(pgv->Gn[ic].dijk[id]);
				}
			}
		}
	}
}
// line 2211; checked
if(pgv->verbose && pp->myrank==0){
	std::cout<<pgv->clit<<"Starting dump_Gnodes_Ensight of file "<<ppar->trim(name)<<" ... Done"<<std::endl;
}
};
}
void io::writeCaseFileHeader(const std::string name){

std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<bl> pbl(new bl);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<band> pba(new band);

std::string chelp6(6, ' '); // was null
int filenameIncr;
#ifdef WITH_CDP
if(pgv->dump_dstep_fs > 0) filenameIncr = pgv->dump_dstep_fs;
else filenameIncr = 1;
#else
filenameIncr = 1;
#endif
chelp6 = std::to_string(filenameIncr); // check the format
global_file = ppar->trim(pgv->dump_dir)+ppar->trim(name)+".case";
std::fstream myfile(global_file, std::ios::out);
if(!myfile.is_open()) std::cout<<"cannot open file line 2194 io.cpp"<<std::endl;
myfile<< "FORMAT"<<std::endl; // format
myfile<< "type:	ensight gold"<<std::endl; // format
myfile<< "GEOMETRY" << std::endl; // format
myfile<< "model:	1	sol_"<<ppar->trim(name)<<"_******.geo"<<std::endl; // format
myfile<< "VARIABLE"<<std::endl; // format
for(auto iv=5; iv<nVarDump+5; iv++){ // this is temperory
	if(isVector(ppar->trim(varDumpName[iv]))){
		myfile<<"vector per node:	1	"<<ppar->trim(varDumpName[iv])<<"	sol_"<<ppar->trim(name);
		myfile<<"_******_"<<ppar->trim(varDumpName[iv])<<".vector"<<std::endl; // format
	}
	else{
		myfile<<"scalar per node:	1	"<<ppar->trim(varDumpName[iv])<<"	sol_"<<ppar->trim(name);
		myfile<<"_******_"<<ppar->trim(varDumpName[iv])<<".scalar"<<std::endl;
	}
}
myfile<<"TIME"<<std::endl; // format
myfile<<"time set:	1"<<std::endl; // format
myfile<<"number of steps: "<<dump_nfiles<<std::endl; // format
myfile<<"filename start number: 000000"<<std::endl; // format
std::stringstream sss;
sss<<std::setw(6)<<std::setfill('0')<<chelp6;
std::string ssss = sss.str();
myfile<<"filename increment: "<<ssss<<std::endl; // format
myfile<<"time values: "<<std::endl; // format
myfile.close(); // we should close it; checked;
};

void io::writeCaseFileCellHeader(const std::string name){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<bl> pbl(new bl);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<band> pba(new band);

std::string chelp6(6, ' '); // was null
int filenameIncr;
#ifdef WITH_CDP
if(pgv->dump_dstep_fs > 0) filenameIncr = pgv->dump_dstep_cell_fs;
else filenameIncr = 1;
#else
filenameIncr = 1;
#endif
chelp6 = filenameIncr;
global_file = ppar->trim(pgv->dump_dir)+ppar->trim(name)+"_cell.case"; //checked
std::fstream myfile(global_file, std::ios::out); // it is new file; checked
if(!myfile.is_open()) std::cout<<"cannot open the file line 2078 io.cpp"<<std::endl;
myfile<< "FORMAT" <<std::endl; // format
myfile<<	"type:	ensight gold"<<std::endl; // format
myfile<< "GEOMETRY"<<std::endl; // format
myfile<< "model:	1	sol_cell_"<<ppar->trim(name)<<"_******.geo"<<std::endl; // format
myfile<< "VARIABLE" <<std::endl; // format
for(auto iv=0; iv<nVarCellDump; iv++){
	if(isVector(ppar->trim(varCellDumpName[iv]))){
		myfile<<"vector per element:	1	"<<ppar->trim(varCellDumpName[iv])<<"	sol_cell_";//format
		myfile<<ppar->trim(name)<<"_******_"<<ppar->trim(varCellDumpName[iv])<<".vector"<<std::endl;
	}
	else{
		myfile<<"scalar per element:	1	"<<ppar->trim(varCellDumpName[iv])<<"	sol_cel_";//format
		myfile<<ppar->trim(name)<<"_******_"<<ppar->trim(varCellDumpName[iv])<<".scalar"<<std::endl;
	}
}
myfile<<"TIME"<<std::endl; // format
myfile<<"time set:	1"<<std::endl; // format
myfile<<"number of steps: "<<dump_nfiles<<std::endl; // format
myfile<<"filename start number: 000000"<<std::endl; // format
std::stringstream sss;
sss<<std::setw(6)<<std::setfill('0')<<chelp6;
std::string ssss = sss.str();
myfile<<"filename increment: "<<ssss<<std::endl; // format
myfile<<"time values: "<<std::endl; // format
myfile.close(); // we need this; checked
}

void io::write_new_case_file(const std::string name, const int typ){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<bl> pbl(new bl);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<band> pba(new band);

std::string chelp6(6, ' '); // was null
dump_nfiles = 0;
if(typ == 1) writeCaseFileHeader(name);
else writeCaseFileCellHeader(name);
std::fstream myfile(global_file, std::ios::out|std::ios::app); // it is new file; checked;
if(myfile.good()) myfile<<pgv->time<<std::endl;
else{
	std::cout<<"The file is not open in write_new_case_file io.cpp"<<std::endl;
	std::exit(1);
}
myfile.close();
};

void io::append_case_file(std::string name, const int typ){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<reinit> prein(new reinit);
std::unique_ptr<bl> pbl(new bl);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<sg> psg(new sg);
std::unique_ptr<band> pba(new band);

std::string chelp6(6, ' ');
std::string line, garbage;
std::vector<double> th;
std::fstream myfile;
string temp_file;
if(typ == 1){
	temp_file = ppar->trim(pgv->dump_dir)+ppar->trim(name)+".case";
//	myfile.open(temp_file,std::ios::app); // append; checked
	myfile.open(temp_file);
}
else{
	temp_file = ppar->trim(pgv->dump_dir)+ppar->trim(name)+"_cell.case";
	myfile.open(temp_file,std::ios::app); // append; // checked
//	myfile.open(temp_file);
}
if(!myfile.is_open()) std::cout<<"cannot open file line 2337 or 2341 io.cpp"<<std::endl;
getline(myfile,garbage);
getline(myfile,garbage);
getline(myfile,garbage);
getline(myfile,garbage);
getline(myfile,garbage);
if(typ == 1){
	for(auto iv=0; iv<nVarDump; iv++) getline(myfile,garbage);
}
else{
	for(auto iv=0; iv<nVarCellDump; iv++) getline(myfile,garbage);
}
getline(myfile,garbage);
getline(myfile,garbage);
getline(myfile,line);
auto index = line.find(':')+1;
std::string temp_string = line.substr(index, line.size()-index);
dump_nfiles = std::stoi(temp_string);
//std::cout<<dump_nfiles<<std::endl; std::exit(1);
std::cout<<dump_nfiles<<std::endl;
getline(myfile,garbage);
getline(myfile,garbage);
getline(myfile,garbage);
if(pp->ierr == 0) th.resize(dump_nfiles);
else pp->litError("append_case_file","allocation error for th");
for(auto ic=0; ic<dump_nfiles; ic++){
	myfile>>th[ic]; // checked
}
myfile.close();
int ich = 0;
for(auto ic=0; ic<dump_nfiles; ic++){
	ich += 1;
	if(th[ic] > pgv->time){
		ich -= 1;
		break;
	}
}
dump_nfiles = ich;
dump_nfiles += 1;
if(typ == 1) writeCaseFileHeader(name);
else writeCaseFileCellHeader(name);
std::ofstream myfile1;
myfile1.open(global_file, std::ios::app);
myfile1.precision(15);
for(auto ic=0; ic<dump_nfiles-1; ic++){
	myfile1<<th[ic]<<std::endl;
}
myfile1<<pgv->time<<std::endl;
myfile1.close();
};

string io::global_file; // check if we need to initialize it;
int io::dump_nfiles = 0;
int io::nVarDump = 0;
int io::nVarCellDump = 0;
std::vector<string> io::varDumpName(99);
std::vector<string> io::varCellDumpName(99);
// check line 2139 of C++ io.cpp
