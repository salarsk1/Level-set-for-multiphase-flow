//Written by Salar Safarkhani

#pragma once
#include<vector>
#include<string>
#include<memory>
#include<array>
#include<algorithm>
#include<iostream>
class block_t; 
class block_p_t;
class Gnode_t{
	public:
	double G;
	std::array<double, 3> V;
	std::array<int, 3> ijk;
	std::array<int, 3> dijk;
};

class Gnode_short_t {
	public:
	double G;
	std::array<int,3> ijk;
};


class shape_t {
	public:
	int code;
	std::array<double, 10> rdata;
};

class global_variable {
public:
	static int IIII;
	int lit_iopen();
	int lit_iclose(int);
	static bool verbose;
	const double pi = 3.141592654;
	static const int nbands = 5;
	static std::array<int, nbands> band_size;
	static int nBandLayers; 
	static double G_max, G_min;
	double G_Tband, G_Bband;
	double G_band_Gm3B, G_band_rGmB3;
	const int sg_notset = 47114711;
	const int bCenterLine = 3;
	const int bNeumann = 2;
	const int bConnect = 1;
	static int*** sg_rank_block;
	static bool*** sg_active;
	static block_t *bl;
	static block_p_t *ibl2bl;
	static int nbl_gl;
	static int nbl;
	static int nghost;
	static std::array<int, 3> ijkm_gl, ijkm_sg, ijkm_bl;
	static std::array<bool, 3> sg_periodic;
	static int max_bl;
	static std::unique_ptr<double []> xc, yc, zc, lxf, lyf, lzf;
	static std::unique_ptr<double []> ldx, ldy, ldz, lrdx, lrdy, lrdz, ryc;
	static std::array<double, 3> xyzs_sg, xyze_sg, xyzs_init_bb, xyze_init_bb,\
										  dxyz_sg, rdxyz_sg, dxyz, rdxyz, rdxyz2;
	static bool cylindrical;
	static double dxyz_min, dxyz_min_2, cell_volume, delta_width_factor;
	static double kappa_min, kappa_max;
	static int nVelFilter;
	static int calcCurvatureType;
	const int ccNodes = 1;
	const int ccDirect = 2;
	const int ccTriLinear = 3;
	const int ccTriCubic = 4;
	const int nneighbors = 6;
	const int dijk_neighbor[3][6] = \
	{{-1, 0, 0, 1, 0, 0}, {0, -1, 0, 0, 1, 0}, {0, 0, -1, 0, 0, 1}};

	int dijk_nc[26][3] = \
	{
		{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1},{-1,0,-1},{1,0,-1},{ 0,-1,-1},{0,1,-1},\
		{-1,0,1},{1,0,1},{0,-1,1},{0,1,1},{-1,-1,0},{1,-1,0},{-1,1,0},{1,1,0},{-1,-1,-1},{1,-1,-1},\
		{-1,1,-1},{1,1,-1},{-1,-1,1},{1,-1,1},{-1,1,1},{1,1,1}
	};

	static std::array<int, 26> ibs,ibe,jbs,jbe,kbs,kbe,igs,ige,jgs,jge,kgs,kge;
	const std::vector<int> ibf2 = {2, 1, 4, 3, 6, 5, 12, 11, 14, 13, 8, 7,\
											 10, 9, 18, 17, 16, 15, 26, 25, 24, 23, 22,\
											 21, 20, 19};
	static std::array<int, 26> gc1_start;
	static int nmax_gc;
	const std::string clit = "LIT: ";
	static std::string schemeReinit, schemeAdvect;
	static std::array<std::array<double,3>,  4> ucWeight3;
	static std::array<std::array<double,3>,  6> ucWeight5;
	static std::array<std::array<double,3>,  8> ucWeight7;
	static std::array<std::array<double,3>, 10> ucWeight9;
	static std::array<std::array<double,3>, 12> ucWeight11;
	
	
	static int fs_step, dump_dstep_cell_fs, maxn_wo_reinit, dump_dstep_fs;
	int max_wo_reinit; 
	
	static double time, fs_dt, dump_dtime_fs, dump_dtime_cell_fs;
	const double CFL_advection = 1.0;
	const double CFL_reinit = 0.5;
	bool do_restart;
	static bool do_diagnostics;
	static std::string reinit_solver, diagnostics_type;
	static int loadBalanceband;
	static std::string loadBalancer;
	static std::array<double, 3> fs_dxyzmax;
	static double fs_lengthScale;
	static double fs_max_drop_vol;
	int slender_size;
	static std::string case_name;
	static std::string dump_band, dump_band_cell;
	static std::string dump_dir;
	static std::array<double, 3> holder_xyzs, holder_xyze, holder_center;
	static double holder_radius, holder_dr, holder_dt;
	static std::string holder_name;
	bool holder;
	static std::string init_shape;
	static std::array<double, 3> init_center, init_xyz_min, init_xyz_max, init_center2, init_normal;
	
	static double init_radius, init_radius2, init_hole_radius,\
		    init_rim_radius, init_width, init_height, init_amplitude,\
			 init_wavelength, init_angle, init_length, init_a, init_b,\
			 init_radius_min, init_radius_max, init_volume;
	
	static bool rotating_cyl_mfix;
	static int init_mod, init_number;
	static double *a1_random, *phi1_random, *a2_random, *phi2_random;
	static int n1_random, n2_random;
	static int* ibl2buf;
	static block_t *b;
	static int ***i2c;
	static Gnode_t *Gn;
	static int ioff;
	static int lit_fileindex;
	static std::array<int, 128> lit_iunits;
	void setBlockPointers(const int ibl);
	void shift_ijk_sg(const std::vector<int> &, const std::vector<int> &, std::vector<int>&,\
							std::vector<int>&);
	void shift_ijk_sg_cyl(const std::vector<int>&, const std::vector<int>&,\
								 std::vector<int>&, std::vector<int>&);
	void lit_die(std::string);
};

class block_t {
	public:
	block_t();
	block_t(const block_t& other);
	block_t &operator=(const block_t&);
	~block_t();
	std::array<int, 3> ijkm;
	std::vector<int> ijk_sg;
	std::array<int, 3> ijk0;
	std::array<std::array<int, 26>, 3> ijk_b2g;
//	std::array<26,std::array<>int ijk_b2g[26][3];
	std::array<int, 26> boundType;
	int imin_, imax_, imin1_, imax1_;
	int jmin_, jmax_, jmin1_, jmax1_;
	int kmin_, kmax_, kmin1_, kmax1_;
	int imino_, imaxo_, jmino_, jmaxo_, kmino_, kmaxo_;
	int nG, nGmax;
	Gnode_t *Gnodes;
	Gnode_t *GnodesOld;
	int ***ijk2ic;
	int NinA = 0;
	int NinT = 0;
	int NinN = 0;
	int NinW = 0;
	int NinX = 0;
	int NinZ = 0;
	std::array<int, global_variable::nbands> NinBand;
	std::unique_ptr<int[]> ic0_bl;
	bool ***person;
	std::unique_ptr<int[]> NinBandLayer;
	std::array<int, 26> GnodeGhost_s;
	std::array<int, 26> nGnodeb;
	block_t *next;
};
class block_p_t{
public:
	block_t *p;
};


