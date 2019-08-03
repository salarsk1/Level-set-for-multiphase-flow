//Written by Salar Safarkhani

#include<iostream>
#include<algorithm>
#include "datag.h"
void global_variable::setBlockPointers(const int ibl) {
	b = ibl2bl[ibl].p;
	Gn =  b -> Gnodes;
	i2c = b -> ijk2ic;
	ioff = ibl2buf[ibl];
};
void global_variable::shift_ijk_sg(const std::vector<int> &ijk_sg, \
											  const std::vector<int> &dir,\
											  std::vector<int> &ijk_sg_shift,\
											  std::vector<int> &ijk_shift) {
	
	for(auto itr=0; itr<ijk_sg_shift.size(); itr++) ijk_sg_shift[itr] = ijk_sg[itr] + dir[itr];
	std::fill_n(ijk_shift.begin(), ijk_shift.size(), 0);
	for(auto i=0; i<sg_periodic.size(); i++){
		if(sg_periodic[i] && ijk_sg_shift[i] == 0) {
			ijk_sg_shift[i] = ijkm_sg[i];
			ijk_shift[i] = ijkm_gl[i];
		}
	}
	for(auto i=0; i<sg_periodic.size(); i++){
		if(sg_periodic[i] && ijk_sg_shift[i]==ijkm_sg[i]+1){
			ijk_sg_shift[i] = 1;
			ijk_shift[i] = - ijkm_gl[i];
		}
	}
};
void global_variable::shift_ijk_sg_cyl(const std::vector<int> &ijk_sg, \
											  		const std::vector<int> &dir,\
													std::vector<int> &ijk_sg_shift,\
													std::vector<int> &ijk_shift) {
for(auto i=0;i<ijk_sg.size();i++) ijk_sg_shift[i] = ijk_sg[i]+dir[i];
std::fill_n(ijk_shift.begin(), ijk_shift.size(), 0);
if(ijk_sg_shift[1] ==0){
	if((ijkm_sg[2]>1 || dir[0] !=0) && !(ijkm_sg[2]==2 && dir[2]!=0 && dir[0]!=0)){
		ijk_sg_shift[1] = ijk_sg[1];
		ijk_shift[1] = 0;
		if(ijk_sg[2]<=ijkm_sg[2]/2){
			ijk_sg_shift[2] = ijk_sg[2] + ijkm_sg[2]/2;
			ijk_shift[2] = ijkm_gl[2]/2;
			if(dir[2]!=0){
				ijk_sg_shift[2] = ijk_sg_shift[2] - dir[2];
				ijk_shift[2] = dir[2] * ijkm_gl[2]/2;
				if(ijk_sg_shift[2]==ijkm_sg[2]+1){
					ijk_sg_shift[2] = 1;
					ijk_shift[2] = -ijkm_gl[2]/2;
				}
				if(ijk_sg[2]<ijkm_sg[2]/2 && dir[2] == -1){
					ijk_shift[2] = ijkm_gl[2]/2;
				};
			}
		}
		else{
			ijk_sg_shift[2] = ijk_sg[2] - ijkm_sg[2]/2;
			ijk_shift[2] = -ijkm_gl[2]/2;
			if(dir[2]!=0){
				ijk_sg_shift[2] = ijk_sg_shift[2] - dir[2];
				ijk_shift[2] = dir[2] * ijkm_gl[2]/2;
				if(ijk_sg_shift[2]==0){
					ijk_sg_shift[2] = ijkm_sg[2];
					ijk_shift[2] = ijkm_gl[2]/2;
				}
				if(ijk_sg[2]>ijkm_sg[2]/2+1 && dir[2]==1){
					ijk_shift[2] = -ijkm_gl[2]/2;
				}
			}
		};
		if(ijkm_sg[2]>1){
			if(ijk_sg_shift[0]==0){
				if(!(ijkm_sg[2]==2 && dir[2] != 0)){
					ijk_sg_shift[0] = 1;
				}
			}
			if(ijk_sg_shift[0]==ijkm_sg[0]+1){
				if(!(ijkm_sg[2]==2 && dir[2]!=0)){
					ijk_sg_shift[0] = ijkm_sg[0];
				}
			}
		}
	}
}
for(auto i=0; i < sg_periodic.size(); i++){
	if(sg_periodic[i]!=0 && ijk_sg_shift[i] ==0){
		ijk_sg_shift[i] = ijkm_sg[i];
		ijk_shift[i] = ijkm_gl[i];
	}
};
for(auto i=0; i < sg_periodic.size(); i++){
	if(sg_periodic[i]!=0 && ijk_sg_shift[i] == ijkm_sg[i]+1){
		ijk_sg_shift[i] = 1;
		ijk_shift[i] = -ijkm_gl[i];
	}
};
};

void global_variable::lit_die(std::string str) {
	std::cout<< this -> clit << ',' << "*******************************" << std::endl;
	std::cout<< this -> clit << ',' << "will kill job to error: " << str << std::endl;
	std::exit(1);
};

int global_variable::lit_iopen(){
	int icall = 1;//it might be global variable; //checked
	if(icall == 1){
		lit_fileindex = 1;
		icall = 0;
	}
	lit_iunits[lit_fileindex - 1] = 0;
	int i;
	for(i=0; i<lit_fileindex; i++){
		if(lit_iunits[i] == 0) break;
	}
	if(i == lit_fileindex-1){
		lit_fileindex += 1;
		if(lit_fileindex >= 128) lit_die("iopen: maximum units number exceeded");
	}
	lit_iunits[i] = 1;
	return i+128;
}

int global_variable::lit_iclose(int iu){// the return checked
	iu -= 128;
	if(iu >0 && iu < lit_fileindex){
		lit_iunits[iu-1] = 0;
		return iu + 128;
	}
	else return -1;
}

bool global_variable::verbose = false;
std::array<int, global_variable::nbands> global_variable::band_size;
int global_variable::nBandLayers = 0;
double global_variable::G_max = 0.0;
double global_variable::G_min = 0.0;
int*** global_variable::sg_rank_block = nullptr;
bool*** global_variable::sg_active = nullptr;
block_t* global_variable::bl = nullptr;
block_p_t* global_variable::ibl2bl = nullptr;
int global_variable::nbl_gl = 0;
int global_variable::nbl = 0;
int global_variable::nghost = 0;
std::array<int, 3> global_variable::ijkm_gl;
std::array<int, 3> global_variable::ijkm_sg;
std::array<int, 3> global_variable::ijkm_bl;
std::array<bool, 3> global_variable::sg_periodic;
int global_variable::max_bl = 512;
std::unique_ptr<double []> global_variable::xc;
std::unique_ptr<double []> global_variable::yc;
std::unique_ptr<double []> global_variable::zc;
std::unique_ptr<double []> global_variable::lxf;
std::unique_ptr<double []> global_variable::lyf;
std::unique_ptr<double []> global_variable::lzf;
std::unique_ptr<double []> global_variable::ldx;
std::unique_ptr<double []> global_variable::ldy;
std::unique_ptr<double []> global_variable::ldz;
std::unique_ptr<double []> global_variable::lrdx;
std::unique_ptr<double []> global_variable::lrdy;
std::unique_ptr<double []> global_variable::lrdz;
std::unique_ptr<double []> global_variable::ryc;
std::array<double, 3> global_variable::xyzs_sg;
std::array<double, 3> global_variable::xyze_sg;
std::array<double, 3> global_variable::xyzs_init_bb;
std::array<double, 3> global_variable::xyze_init_bb;
std::array<double, 3> global_variable::dxyz_sg;
std::array<double, 3> global_variable::rdxyz_sg;
std::array<double, 3> global_variable::dxyz;
std::array<double, 3> global_variable::rdxyz;
std::array<double, 3> global_variable::rdxyz2;
bool global_variable::cylindrical = false;
double global_variable::dxyz_min = 0.0;
double global_variable::dxyz_min_2 = 0.0;
double global_variable::cell_volume = 0.0;
double global_variable::delta_width_factor = 0.0;
double global_variable::kappa_min = 0.0;
double global_variable::kappa_max = 0.0;
int global_variable::nVelFilter = 0;
int global_variable::calcCurvatureType = 0;
std::array<int, 26> global_variable::ibs;
std::array<int, 26> global_variable::ibe;
std::array<int, 26> global_variable::jbs;
std::array<int, 26> global_variable::jbe;
std::array<int, 26> global_variable::kbs;
std::array<int, 26> global_variable::kbe;
std::array<int, 26> global_variable::igs;
std::array<int, 26> global_variable::ige;
std::array<int, 26> global_variable::jgs;
std::array<int, 26> global_variable::jge;
std::array<int, 26> global_variable::kgs;
std::array<int, 26> global_variable::kge;
std::array<int, 26> global_variable::gc1_start;
int global_variable::nmax_gc = 0;
std::string global_variable::schemeReinit;
std::string global_variable::schemeAdvect;
std::array<std::array<double, 3>, 4> global_variable::ucWeight3;
std::array<std::array<double, 3>, 6> global_variable::ucWeight5;
std::array<std::array<double, 3>, 8> global_variable::ucWeight7;
std::array<std::array<double, 3>, 10> global_variable::ucWeight9;
std::array<std::array<double, 3>, 12> global_variable::ucWeight11;
int global_variable::fs_step = 0;
int global_variable::dump_dstep_cell_fs = 0;
int global_variable::maxn_wo_reinit = 0;
int global_variable::dump_dstep_fs = 0;
double global_variable::time = 0.0;
double global_variable::fs_dt = 0.0;
double global_variable::dump_dtime_fs = 0.0;
double global_variable::dump_dtime_cell_fs = 0.0;
bool global_variable::do_diagnostics = false;
std::string global_variable::reinit_solver;
std::string global_variable::diagnostics_type;
int global_variable::loadBalanceband = 0;
std::string global_variable::loadBalancer;
std::array<double, 3> global_variable::fs_dxyzmax;
double global_variable::fs_lengthScale = -1.0;
double global_variable::fs_max_drop_vol = 0.0;
std::string global_variable::case_name;
std::string global_variable::dump_band;
std::string global_variable::dump_band_cell;
std::string global_variable::dump_dir = "dump/";
std::array<double, 3> global_variable::holder_xyzs;
std::array<double, 3> global_variable::holder_xyze;
std::array<double, 3> global_variable::holder_center;
double global_variable::holder_radius = 0.0;
double global_variable::holder_dr = 0.0;
double global_variable::holder_dt = 0.0;
std::string global_variable::holder_name;
std::string global_variable::init_shape;
std::array<double, 3> global_variable::init_center;
std::array<double, 3> global_variable::init_center2;
std::array<double, 3> global_variable::init_normal;
std::array<double, 3> global_variable::init_xyz_min;
std::array<double, 3> global_variable::init_xyz_max;
double global_variable::init_radius = 0.0;
double global_variable::init_radius2 = 0.0;
double global_variable::init_hole_radius = 0.0;
double global_variable::init_rim_radius = 0.0;
double global_variable::init_width = 0.0;
double global_variable::init_height = 0.0;
double global_variable::init_amplitude = 0.0;
double global_variable::init_wavelength = 0.0;
double global_variable::init_angle = 0.0;
double global_variable::init_length = 0.0;
double global_variable::init_a = 0.0;
double global_variable::init_b = 0.0;
double global_variable::init_radius_min = 0.0;
double global_variable::init_radius_max = 0.0;
double global_variable::init_volume = 0.0;
bool global_variable::rotating_cyl_mfix = false;
int global_variable::init_mod = 0;
int global_variable::init_number = 0;
double* global_variable::a1_random = nullptr;
double* global_variable::phi1_random = nullptr;
double* global_variable::a2_random = nullptr;
double* global_variable::phi2_random = nullptr;
int global_variable::n1_random = 0;
int global_variable::n2_random = 0;
int* global_variable::ibl2buf = nullptr;
block_t* global_variable::b = nullptr;
int*** global_variable::i2c = nullptr;
Gnode_t* global_variable::Gn = nullptr;
int global_variable::ioff = 0;
int global_variable::lit_fileindex = 0;
std::array<int, 128> global_variable::lit_iunits;
int global_variable::IIII = 0;
block_t::block_t(){
	std::cout<<"block_t constructure was called"<<std::endl;
//	ijk2ic = nullptr;
//	person = nullptr;
};
block_t::block_t(const block_t& other)
{
/*	std::cout<<"block_t copy constructure was called"<<std::endl;
	global_variable GV;
	imino_ = other.imino_; imaxo_ = other.imaxo_;
	jmino_ = other.jmino_; jmaxo_ = other.jmaxo_;
	kmino_ = other.kmino_; kmaxo_ = other.kmaxo_;
	ijkm = other.ijkm;
	ijk_sg = other.ijk_sg;
	ijk0 = other.ijk0;
	ijk_b2g = other.ijk_b2g;
	boundType = other.boundType;
	imin_ = other.imin_; imax_ = other.imax_;
	jmin_ = other.jmin_; jmax_ = other.jmax_;
	kmin_ = other.kmin_; kmax_ = other.kmax_;
	nG = other.nG; nGmax = other.nGmax;
	NinA = other.NinA;
	NinT = other.NinT;
	NinW = other.NinW;
	NinX = other.NinX;
	NinZ = other.NinZ;
	NinBand = other.NinBand;

	GnodeGhost_s = other.GnodeGhost_s;
	nGnodeb = other.nGnodeb;
	next = other.next;
	imin1_ = other.imin1_; imax1_ = other.imax1_;
	jmin1_ = other.jmin1_; jmax1_ = other.jmax1_;
	kmin1_ = other.kmin1_; kmax1_ = other.kmax1_;
	Gnodes.reset(new Gnode_t[nGmax]);
	for(auto i=0; i<nGmax; i++) Gnodes[i] = other.Gnodes[i];
	GnodesOld.reset(new Gnode_t[nGmax]);
	for(auto i=0; i<nGmax; i++) GnodesOld[i] = other.GnodesOld[i];
	ic0_bl.reset(new int[GV.nBandLayers]);
	for(auto i=0; i<GV.nBandLayers; i++) ic0_bl[i] = other.ic0_bl[i];
	NinBandLayer.reset(new int[GV.nBandLayers]);
	for(auto i=0; i<GV.nBandLayers; i++) NinBandLayer[i] = other.NinBandLayer[i];
	auto d1 = imaxo_ - imino_ + 1;
	auto d2 = jmaxo_ - jmino_ + 1;
	auto d3 = kmaxo_ - kmino_ + 1;
	std::cout<<d1<<'\t'<<d2<<'\t'<<d3<<std::endl;
	for(auto i=0; i<d1; i++){
		for(auto j=0; j<d2; j++){
			for(auto k=0; k<d3; k++){
				ijk2ic[i][j][k] = other.ijk2ic[i][j][k];
			}
		}
	}
	for(auto i=0; i<d1; i++){
		for(auto j=0; j<d2; j++){
			for(auto k=0; k<d3; k++){
				person[i][j][k] = other.person[i][j][k];
			}
		}
	}*/
};
block_t& block_t::operator=(const block_t& other){
/*	std::cout<<"block_t copy assignment was called"<<std::endl;
	global_variable GV;
	imino_ = other.imino_; imaxo_ = other.imaxo_;
	jmino_ = other.jmino_; jmaxo_ = other.jmaxo_;
	kmino_ = other.kmino_; kmaxo_ = other.kmaxo_;
	ijkm = other.ijkm;
	ijk_sg = other.ijk_sg;
	ijk0 = other.ijk0;
	ijk_b2g = other.ijk_b2g;
	boundType = other.boundType;
	imin_ = other.imin_; imax_ = other.imax_;
	jmin_ = other.jmin_; jmax_ = other.jmax_;
	kmin_ = other.kmin_; kmax_ = other.kmax_;
	nG = other.nG; nGmax = other.nGmax;
	for(auto i=0; i<nGmax; i++) Gnodes[i] = other.Gnodes[i];
	NinA = other.NinA;
	NinT = other.NinT;
	NinW = other.NinW;
	NinX = other.NinX;
	NinZ = other.NinZ;
	NinBand = other.NinBand;
	GnodeGhost_s = other.GnodeGhost_s;
	nGnodeb = other.nGnodeb;
	next = other.next;
	imin1_ = other.imin1_; imax1_ = other.imax1_;
	jmin1_ = other.jmin1_; jmax1_ = other.jmax1_;
	kmin1_ = other.kmin1_; kmax1_ = other.kmax1_;
	Gnodes.reset(new Gnode_t[nGmax]);
	for(auto i=0; i<nGmax; i++) Gnodes[i] = other.Gnodes[i];
	GnodesOld.reset(new Gnode_t[nGmax]);
	for(auto i=0; i<nGmax; i++) GnodesOld[i] = other.GnodesOld[i];
	ic0_bl.reset(new int[GV.nBandLayers]);
	for(auto i=0; i<GV.nBandLayers; i++) ic0_bl[i] = other.ic0_bl[i];
	NinBandLayer.reset(new int[GV.nBandLayers]);
	for(auto i=0; i<GV.nBandLayers; i++) NinBandLayer[i] = other.NinBandLayer[i];
	auto d1 = imaxo_ - imino_ + 1;
	auto d2 = jmaxo_ - jmino_ + 1;
	auto d3 = kmaxo_ - kmino_ + 1;
	std::cout<<d1<<'\t'<<d2<<'\t'<<d3<<std::endl;
	for(auto i=0; i<d1; i++){
		for(auto j=0; j<d2; j++){
			for(auto k=0; k<d3; k++){
				ijk2ic[i][j][k] = other.ijk2ic[i][j][k];
			}
		}
	}
	for(auto i=0; i<d1; i++){
		for(auto j=0; j<d2; j++){
			for(auto k=0; k<d3; k++){
				person[i][j][k] = other.person[i][j][k];
			}
		}
	}

	return *this;*/

}
block_t::~block_t(){
/*	std::cout<<"block_t destructure was called"<<std::endl;
   size_t d1 = imaxo_ - imino_ + 1;
   size_t d2 = jmaxo_ - jmino_ + 1;
	if(ijk2ic != nullptr){
	   for(auto i=0;i<d1;i++){
	      for(auto j=0;j<d2;j++){
	         delete [] ijk2ic[i][j];
	      }
	      delete [] ijk2ic[i];
	   };
	   delete [] ijk2ic;
		ijk2ic = nullptr;
	}
	if(person != nullptr){
		for(auto i=0;i<d1;i++){
	      for(auto j=0;j<d2;j++){
	         delete [] person[i][j];
	      }
	      delete [] person[i];
	   };
	   delete [] person;
		person = nullptr;
	}*/
};

