//Written by Salar Safarkhani

#pragma once
#include<vector>
#include<memory>
#include<cmath>
#include<limits>
#include"datag.h"
#include"parallel.h"
#include"param.h"
#include"timing.h"
#include"monitor.h"
#include"bl.h"
#include"bound.h"
#include"gnodes.h"
#include"init.h"
class ghostcloth_bl_t;
class band{
	public:
	void band_m_init();
	void band_init();
	void band_inject(const int nshapes, shape_t* shapes);
	void band_seed_initG();
	void band_seed();
	void band_grow(bool init, char ch);
	bool is_ghostcloth(block_t* bb, const Gnode_t &Gnh, const int ibl);
	void band_push_ghostcloth(const int ib);
	void calc_partner_block_cyl(block_t*,const int,int &,int &, int &, bool &);
	void calc_partner_block(block_t*,const int,int &,int &, int&, bool &);
	void band_glue_cloth(double*,int**,const int,const int,const int);
	void band_set_block_data();
	void band_regenerate_ijk2ic();
	void band_cleanup();
	void band_m_cleanup();
	void band_monitor();
private:
	static std::vector<int> skin_s, skin_e, cloth_s, cloth_e, nskin, ncloth, ncount;
	static std::vector<ghostcloth_bl_t> ghostcloth_bl;
	static int bandLayerCounter;
	static std::vector<std::vector<int>> NinBand;
	void is_ghostcloth_bl();
};

class ghostcloth_bl_t{
public:
//	ghostcloth_bl_t(){};
//	~ghostcloth_bl_t(){
//		delete [] G_gc;
//		for(auto i=0; i<3; i++){
//			delete [] ijk_gc[i];
//			delete [] dijk_gc[i];
//		}
//		delete [] ijk_gc;
//		delete [] dijk_gc;
//	};
	std::array<int, 26> nin_bf;
	double* G_gc  = nullptr;
	int** ijk_gc   = nullptr;
	int** dijk_gc = nullptr;
private:
	size_t G_gc_S;
	size_t ijk_gc_S1, ijk_gc_S2;
	size_t dijk_gc_S1, dijk_gc_S2;
};
