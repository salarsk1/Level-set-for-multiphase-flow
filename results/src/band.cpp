//Written by Salar Safarkhani

#include"band.h"
void band::band_m_init(){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<litparam>ppar(new litparam);
std::unique_ptr<monitor>pmo(new monitor);

pgv->band_size[0] = 1;
pgv->band_size[1] = 5+pgv->nghost;
pgv->band_size[2] = std::max(3, pgv->nghost);
pgv->band_size[3] = pgv->nghost;
pgv->band_size[1] = ppar->get_integer_param("LIT.BAND_SIZE_T",pgv->band_size[1],2);
pgv->band_size[2] = ppar->get_integer_param("LIT.BAND_SIZE_N",pgv->band_size[2],2);
pgv->band_size[4] = ppar->get_integer_param("LIT.BAND_SIZE_X",-1,2);
if(pgv->band_size[4]<0){
	if(pgv->cylindrical){
		if(pgv->fs_lengthScale<0.0){
			auto m1 = std::max(pgv->fs_dxyzmax[0]/pgv->dxyz[0],\
									 pgv->fs_dxyzmax[1]/pgv->dxyz[1]);
			m1 = std::max(m1,pgv->fs_dxyzmax[2]/pgv->dxyz[2]);
			m1 = ceil(m1);
			auto c1 = m1-(pgv->band_size[2]+pgv->band_size[3]);
			pgv->band_size[4]=std::max(0.0, c1);
		}
		else{
			auto m1 = std::min(pgv->dxyz[0], pgv->dxyz[1]);
			m1 = pgv->fs_lengthScale/m1;
			m1 = ceil(m1);
			int c1 = m1-(pgv->band_size[2]+pgv->band_size[3]);
			pgv->band_size[4] = std::max(0, c1);
		}
	}
	else{
		if(pgv->fs_lengthScale<0.0){
			auto m1 = std::max(pgv->fs_dxyzmax[0]/pgv->dxyz[0],pgv->fs_dxyzmax[1]/pgv->dxyz[1]);
			m1 = std::max(m1,pgv->fs_dxyzmax[2]/pgv->dxyz[2]);
			m1 = ceil(m1);
			int c1 = m1-(pgv->band_size[2]+pgv->band_size[3]);
			pgv->band_size[4]=std::max(0, c1);
		}
		else{
			auto m1 = std::min(pgv->dxyz[0], pgv->dxyz[1]);
			m1 = std::min(m1, pgv->dxyz[2]);
			m1 = pgv->fs_lengthScale/m1;
			m1 = ceil(m1);
			int c1 = m1-(pgv->band_size[2]+pgv->band_size[3]);
			pgv->band_size[4] = std::max(0, c1);
		}
	}
}
pgv->nBandLayers = 0;
for(auto i=0; i<5; i++){
	pgv->nBandLayers += pgv->band_size[i];
}
if(pgv->cylindrical){
 if(std::max(pgv->fs_dxyzmax[0]/pgv->dxyz[0],pgv->fs_dxyzmax[1]/pgv->dxyz[1])>1.1){
	pgv->nVelFilter = static_cast<int>(std::max(pgv->fs_dxyzmax[0]/pgv->dxyz[0],\
												           pgv->fs_dxyzmax[1]/pgv->dxyz[1]));
 }
 else{
	pgv->nVelFilter = 0;
 }
}
else{
	auto m1=std::max(pgv->fs_dxyzmax[0]/pgv->dxyz[0], pgv->fs_dxyzmax[1]/pgv->dxyz[1]);
	m1 = std::max(m1, pgv->fs_dxyzmax[2]/pgv->dxyz[2]);
	if(m1 > 1.1)
		pgv->nVelFilter = static_cast<int>(m1);
	else
		pgv->nVelFilter = 0;
}
pgv->nVelFilter = ppar->get_integer_param("LIT.NVELFILTER ",pgv->nVelFilter,2);
if(pgv->verbose && pp->myrank==0){
	std::cout<<pgv->clit<<"setting nVelFilter"<<pgv->nVelFilter<<std::endl;
}
if(pgv->verbose && pp->myrank==0){
	std::cout<<pgv->clit<<pp->cmy<<"setting N-band size to "<< pgv->band_size[2]<<std::endl;
	std::cout<<pgv->clit<<pp->cmy<<"setting X-band size to "<< pgv->band_size[4]<<std::endl;
}
if(pgv->cylindrical){
	if(pgv->ijkm_gl[2]>1){
		auto m1 = std::max(pgv->dxyz[0], pgv->dxyz[1]);
		m1 = std::max(m1, pgv->yc[pgv->ijkm_gl[1]/2-(1-pgv->nghost)]*pgv->dxyz[2]);
		pgv->G_Bband = static_cast<double	>(pgv->band_size[1]-3)*m1;
		pgv->G_Tband = static_cast<double>(pgv->band_size[1])*m1;
		pgv->G_max = (static_cast<double>(pgv->band_size[1]+pgv->band_size[2])-2.0)*m1;
	}
	else{
		auto m1 = std::max(pgv->dxyz[0], pgv->dxyz[1]);
		pgv->G_Bband = static_cast<double>(pgv->band_size[1]-3)*m1;
		pgv->G_Tband = static_cast<double>(pgv->band_size[1])*m1;
		pgv->G_max = (static_cast<double>(pgv->band_size[1]+pgv->band_size[2])-2.0)*m1;
	}
}
else{
		auto m1 = std::max(pgv->dxyz[0], pgv->dxyz[1]);
		m1 = std::max(m1, pgv->dxyz[2]);
		pgv->G_Bband = static_cast<double>(pgv->band_size[1]-3)*m1;
		pgv->G_Tband = static_cast<double>(pgv->band_size[1])*m1;
		pgv->G_max = (static_cast<double>(pgv->band_size[1]+pgv->band_size[2])-2.0)*m1;
}
pgv->G_min = -pgv->G_max;
pgv->G_band_Gm3B = pgv->G_Tband-3.0*pgv->G_Bband;
pgv->G_band_rGmB3 = 1.0/pow(pgv->G_Tband-pgv->G_Bband,3);
pgv->maxn_wo_reinit = std::max(1, pgv->band_size[1]-6);
pgv->gc1_start[0]  = 0;
pgv->gc1_start[1]  = pgv->gc1_start[0] + pgv->ijkm_bl[1]*pgv->ijkm_bl[2];
pgv->gc1_start[2]  = pgv->gc1_start[1] + pgv->ijkm_bl[1]*pgv->ijkm_bl[2];
pgv->gc1_start[3]  = pgv->gc1_start[2] + pgv->ijkm_bl[0]*pgv->ijkm_bl[2];
pgv->gc1_start[4]  = pgv->gc1_start[3] + pgv->ijkm_bl[0]*pgv->ijkm_bl[2];
pgv->gc1_start[5]  = pgv->gc1_start[4] + pgv->ijkm_bl[1]*pgv->ijkm_bl[0];
pgv->gc1_start[6]  = pgv->gc1_start[5] + pgv->ijkm_bl[1]*pgv->ijkm_bl[0];
pgv->gc1_start[7]  = pgv->gc1_start[6] + pgv->ijkm_bl[1];
pgv->gc1_start[8]  = pgv->gc1_start[7] + pgv->ijkm_bl[1];
pgv->gc1_start[9]  = pgv->gc1_start[8] + pgv->ijkm_bl[0];
pgv->gc1_start[10] = pgv->gc1_start[9] + pgv->ijkm_bl[0];
pgv->gc1_start[11] = pgv->gc1_start[10] + pgv->ijkm_bl[1];
pgv->gc1_start[12] = pgv->gc1_start[11] + pgv->ijkm_bl[1];
pgv->gc1_start[13] = pgv->gc1_start[12] + pgv->ijkm_bl[0];
pgv->gc1_start[14] = pgv->gc1_start[13] + pgv->ijkm_bl[0];
pgv->gc1_start[15] = pgv->gc1_start[14] + pgv->ijkm_bl[2];
pgv->gc1_start[16] = pgv->gc1_start[15] + pgv->ijkm_bl[2];
pgv->gc1_start[17] = pgv->gc1_start[16] + pgv->ijkm_bl[2];
pgv->gc1_start[18] = pgv->gc1_start[17] + pgv->ijkm_bl[2];
for(auto i=19; i<26; i++){
	pgv->gc1_start[i] = pgv->gc1_start[i-1] + 1;
};

pgv->ibs[0]  = 1; pgv->ibe[0]  = pgv->nghost;
pgv->jbs[0]  = 1; pgv->jbe[0]  = pgv->ijkm_bl[1];
pgv->kbs[0]  = 1; pgv->kbe[0]  = pgv->ijkm_bl[2];

pgv->ibs[1]  = pgv->ijkm_bl[0]+1-pgv->nghost; pgv->ibe[1]  = pgv->ijkm_bl[0];
pgv->jbs[1]  = 1; pgv->jbe[1]  = pgv->ijkm_bl[1];
pgv->kbs[1]  = 1; pgv->kbe[1]  = pgv->ijkm_bl[2];

pgv->ibs[2]  = 1; pgv->ibe[2]  = pgv->ijkm_bl[0];
pgv->jbs[2]  = 1; pgv->jbe[2]  = pgv->nghost;
pgv->kbs[2]  = 1; pgv->kbe[2]  = pgv->ijkm_bl[2];

pgv->ibs[3]  = 1; pgv->ibe[3]  = pgv->ijkm_bl[0];
pgv->jbs[3]  = pgv->ijkm_bl[1]+1-pgv->nghost; pgv->jbe[3] = pgv->ijkm_bl[1];
pgv->kbs[3]  = 1; pgv->kbe[3]  = pgv->ijkm_bl[2];

pgv->ibs[4]  = 1; pgv->ibe[4]  = pgv->ijkm_bl[0];
pgv->jbs[4]  = 1; pgv->jbe[4]  = pgv->ijkm_bl[1];
pgv->kbs[4]  = 1; pgv->kbe[4]  = pgv->nghost;

pgv->ibs[5]  = 1; pgv->ibe[5]  = pgv->ijkm_bl[0];
pgv->jbs[5]  = 1; pgv->jbe[5]  = pgv->ijkm_bl[1];
pgv->kbs[5]  = pgv->ijkm_bl[2]+1-pgv->nghost; pgv->kbe[5]  = pgv->ijkm_bl[2];

pgv->ibs[6]  = 1; pgv->ibe[6]  = pgv->nghost;
pgv->jbs[6]  = 1; pgv->jbe[6]  = pgv->ijkm_bl[1];
pgv->kbs[6]  = 1; pgv->kbe[6]  = pgv->nghost;

pgv->ibs[7]  = pgv->ijkm_bl[0]+1-pgv->nghost; pgv->ibe[7]=pgv->ijkm_bl[0];
pgv->jbs[7]  = 1; pgv->jbe[7]  = pgv->ijkm_bl[1];
pgv->kbs[7]  = 1; pgv->kbe[7]  = pgv->nghost;

pgv->ibs[8]  = 1; pgv->ibe[8]  = pgv->ijkm_bl[0];
pgv->jbs[8]  = 1; pgv->jbe[8]  = pgv->nghost;
pgv->kbs[8]  = 1; pgv->kbe[8]  = pgv->nghost;

pgv->ibs[9]  = 1; pgv->ibe[9]  = pgv->ijkm_bl[0];
pgv->jbs[9]  = pgv->ijkm_bl[1]+1-pgv->nghost; pgv->jbe[9]  = pgv->ijkm_bl[1];
pgv->kbs[9]  = 1; pgv->kbe[9]  = pgv->nghost;

pgv->ibs[10] = 1; pgv->ibe[10] = pgv->nghost;
pgv->jbs[10] = 1; pgv->jbe[10] = pgv->ijkm_bl[1];
pgv->kbs[10] = pgv->ijkm_bl[2]+1-pgv->nghost; pgv->kbe[10] = pgv->ijkm_bl[2];

pgv->ibs[11] = pgv->ijkm_bl[0]+1-pgv->nghost; pgv->ibe[11] = pgv->ijkm_bl[0];
pgv->jbs[11] = 1; pgv->jbe[11] = pgv->ijkm_bl[1];
pgv->kbs[11] = pgv->ijkm_bl[2]+1-pgv->nghost; pgv->kbe[11] = pgv->ijkm_bl[2];

pgv->ibs[12] = 1; pgv->ibe[12] = pgv->ijkm_bl[0];
pgv->jbs[12] = 1; pgv->jbe[12] = pgv->nghost;
pgv->kbs[12] = pgv->ijkm_bl[2]+1-pgv->nghost; pgv->kbe[12] = pgv->ijkm_bl[2];

pgv->ibs[13] = 1; pgv->ibe[13] = pgv->ijkm_bl[0];
pgv->jbs[13] = pgv->ijkm_bl[1]+1-pgv->nghost; pgv->jbe[13] = pgv->ijkm_bl[1];
pgv->kbs[13] = pgv->ijkm_bl[2]+1-pgv->nghost; pgv->kbe[13] = pgv->ijkm_bl[2];

pgv->ibs[14] = 1; pgv->ibe[14] = pgv->nghost;
pgv->jbs[14] = 1; pgv->jbe[14] = pgv->nghost;
pgv->kbs[14] = 1; pgv->kbe[14] = pgv->ijkm_bl[2];

pgv->ibs[15] = pgv->ijkm_bl[0]+1-pgv->nghost; pgv->ibe[15] = pgv->ijkm_bl[0];
pgv->jbs[15] = 1; pgv->jbe[15] = pgv->nghost;
pgv->kbs[15] = 1; pgv->kbe[15] = pgv->ijkm_bl[2];

pgv->ibs[16] = 1; pgv->ibe[16] = pgv->nghost;
pgv->jbs[16] = pgv->ijkm_bl[1]+1-pgv->nghost; pgv->jbe[16] = pgv->ijkm_bl[1];
pgv->kbs[16] = 1; pgv->kbe[16] = pgv->ijkm_bl[2];

pgv->ibs[17] = pgv->ijkm_bl[0]+1-pgv->nghost; pgv->ibe[17] = pgv->ijkm_bl[0];
pgv->jbs[17] = pgv->ijkm_bl[1]+1-pgv->nghost; pgv->jbe[17] = pgv->ijkm_bl[1];
pgv->kbs[17] = 1; pgv->kbe[17] = pgv->ijkm_bl[2];

pgv->ibs[18] = 1; pgv->ibe[18] = pgv->nghost;
pgv->jbs[18] = 1; pgv->jbe[18] = pgv->nghost;
pgv->kbs[18] = 1; pgv->kbe[18] = pgv->nghost;

pgv->ibs[19] = pgv->ijkm_bl[0]+1-pgv->nghost; pgv->ibe[19] = pgv->ijkm_bl[0];
pgv->jbs[19] = 1; pgv->jbe[19] = pgv->nghost;
pgv->kbs[19] = 1; pgv->kbe[19] = pgv->nghost;

pgv->ibs[20] = 1; pgv->ibe[20] = pgv->nghost;
pgv->jbs[20] = pgv->ijkm_bl[1]+1-pgv->nghost; pgv->jbe[20] = pgv->ijkm_bl[1];
pgv->kbs[20] = 1; pgv->kbe[20] = pgv->nghost;

pgv->ibs[21] = pgv->ijkm_bl[0]+1-pgv->nghost; pgv->ibe[21] = pgv->ijkm_bl[0];
pgv->jbs[21] = pgv->ijkm_bl[1]+1-pgv->nghost; pgv->jbe[21] = pgv->ijkm_bl[1];
pgv->kbs[21] = 1; pgv->kbe[21] = pgv->nghost;

pgv->ibs[22] = 1; pgv->ibe[22] = pgv->nghost;
pgv->jbs[22] = 1; pgv->jbe[22] = pgv->nghost;
pgv->kbs[22] = pgv->ijkm_bl[2]+1-pgv->nghost; pgv->kbe[22] = pgv->ijkm_bl[2];

pgv->ibs[23] = pgv->ijkm_bl[0]+1-pgv->nghost; pgv->ibe[23] = pgv->ijkm_bl[0];
pgv->jbs[23] = 1; pgv->jbe[23] = pgv->nghost;
pgv->kbs[23] = pgv->ijkm_bl[2]+1-pgv->nghost; pgv->kbe[23] = pgv->ijkm_bl[2];

pgv->ibs[24] = 1; pgv->ibe[24] = pgv->nghost;
pgv->jbs[24] = pgv->ijkm_bl[1]+1-pgv->nghost; pgv->jbe[24] = pgv->ijkm_bl[1];
pgv->kbs[24] = pgv->ijkm_bl[2]+1-pgv->nghost; pgv->kbe[24] = pgv->ijkm_bl[2];

pgv->ibs[25] = pgv->ijkm_bl[0]+1-pgv->nghost; pgv->ibe[25] = pgv->ijkm_bl[0];
pgv->jbs[25] = pgv->ijkm_bl[1]+1-pgv->nghost; pgv->jbe[25] = pgv->ijkm_bl[1];
pgv->kbs[25] = pgv->ijkm_bl[2]+1-pgv->nghost; pgv->kbe[25] = pgv->ijkm_bl[2];

for(auto i=0; i<26; i++){
	pgv->igs[i] = pgv->ibs[i];
	pgv->ige[i] = pgv->ibe[i];
	pgv->jgs[i] = pgv->jbs[i];
	pgv->jge[i] = pgv->jbe[i];
	pgv->kgs[i] = pgv->kbs[i];
	pgv->kge[i] = pgv->kbe[i];
}
pgv->igs[0]  = 1-pgv->nghost; pgv->ige[0]  = 0;
pgv->igs[1]  = pgv->ijkm_bl[0]+1; pgv->ige[1]  = pgv->ijkm_bl[1]+pgv->nghost;
pgv->jgs[2]  = 1-pgv->nghost; pgv->jge[2]  = 0;
pgv->jgs[3]  = pgv->ijkm_bl[1]+1; pgv->jge[3]  = pgv->ijkm_bl[1]+pgv->nghost;
pgv->kgs[4]  = 1-pgv->nghost; pgv->kge[4]  = 0;
pgv->kgs[5]  = pgv->ijkm_bl[2]+1; pgv->kge[5]  = pgv->ijkm_bl[2]+pgv->nghost;
pgv->igs[6]  = 1-pgv->nghost; pgv->ige[6]  = 0;
pgv->kgs[6]  = 1-pgv->nghost; pgv->kge[6]  = 0;
pgv->igs[7]  = pgv->ijkm_bl[0]+1; pgv->ige[7]  = pgv->ijkm_bl[0]+pgv->nghost;
pgv->kgs[7]  = 1-pgv->nghost; pgv->kge[7]  = 0;
pgv->jgs[8]  = 1-pgv->nghost; pgv->jge[8]  = 0;
pgv->kgs[8]  = 1-pgv->nghost; pgv->kge[8]  = 0;
pgv->jgs[9]  = pgv->ijkm_bl[1]+1; pgv->jge[9]  = pgv->ijkm_bl[1]+pgv->nghost;
pgv->kgs[9]  = 1-pgv->nghost; pgv->kge[9]  = 0;
pgv->igs[10] = 1-pgv->nghost; pgv->ige[10] = 0;
pgv->kgs[10] = pgv->ijkm_bl[2]+1; pgv->kge[10] = pgv->ijkm_bl[2]+pgv->nghost;
pgv->igs[11] = pgv->ijkm_bl[0]+1; pgv->ige[11] = pgv->ijkm_bl[0]+pgv->nghost;
pgv->kgs[11] = pgv->ijkm_bl[2]+1; pgv->kge[11] = pgv->ijkm_bl[2]+pgv->nghost;
pgv->jgs[12] = 1-pgv->nghost; pgv->jge[12] = 0;
pgv->kgs[12] = pgv->ijkm_bl[2]+1; pgv->kge[12] = pgv->ijkm_bl[2]+pgv->nghost;
pgv->jgs[13] = pgv->ijkm_bl[1]+1; pgv->jge[13] = pgv->ijkm_bl[1]+pgv->nghost;
pgv->kgs[13] = pgv->ijkm_bl[2]+1; pgv->kge[13] = pgv->ijkm_bl[2]+pgv->nghost;
pgv->igs[14] = 1-pgv->nghost; pgv->ige[14] = 0;
pgv->jgs[14] = 1-pgv->nghost; pgv->jge[14] = 0;
pgv->kgs[14] = 1; pgv->kge[14] = pgv->ijkm_bl[2];
pgv->igs[15] = pgv->ijkm_bl[0]+1; pgv->ige[15] = pgv->ijkm_bl[0]+pgv->nghost;
pgv->jgs[15] = 1-pgv->nghost; pgv->jge[16] = 0;
pgv->igs[16] = 1-pgv->nghost; pgv->ige[16] = 0;
pgv->jgs[16] = pgv->ijkm_bl[1]+1; pgv->jge[16] = pgv->ijkm_bl[1]+pgv->nghost;
pgv->igs[17] = pgv->ijkm_bl[0]+1; pgv->ige[17] = pgv->ijkm_bl[0]+pgv->nghost;
pgv->jgs[17] = pgv->ijkm_bl[1]+1; pgv->jge[17] = pgv->ijkm_bl[1]+pgv->nghost;
pgv->igs[18] = 1-pgv->nghost; pgv->ige[18] = 0;
pgv->jgs[18] = 1-pgv->nghost; pgv->jge[18] = 0;
pgv->kgs[18] = 1-pgv->nghost; pgv->kge[18] = 0;
pgv->igs[19] = pgv->ijkm_bl[0]+1; pgv->ige[19] = pgv->ijkm_bl[0]+pgv->nghost;
pgv->jgs[19] = 1-pgv->nghost; pgv->jge[19] = 0;
pgv->kgs[19] = 1-pgv->nghost; pgv->kge[19] = 0;
pgv->igs[20] = 1-pgv->nghost; pgv->ige[20] = 0;
pgv->jgs[20] = pgv->ijkm_bl[1]+1; pgv->jge[20] = pgv->ijkm_bl[1]+pgv->nghost;
pgv->kgs[20] = 1-pgv->nghost; pgv->kge[20] = 0;
pgv->igs[21] = pgv->ijkm_bl[0]+1; pgv->ige[21] = pgv->ijkm_bl[0]+pgv->nghost;
pgv->jgs[21] = pgv->ijkm_bl[1]+1; pgv->jge[21] = pgv->ijkm_bl[1]+pgv->nghost;
pgv->kgs[21] = 1-pgv->nghost; pgv->kge[21] = 0;
pgv->igs[22] = 1-pgv->nghost; pgv->ige[22] = 0;
pgv->jgs[22] = 1-pgv->nghost; pgv->jge[22] = 0;
pgv->kgs[22] = pgv->ijkm_bl[2]+1; pgv->kge[22] = pgv->ijkm_bl[2]+pgv->nghost;
pgv->igs[23] = pgv->ijkm_bl[0]+1; pgv->ige[23] = pgv->ijkm_bl[0]+pgv->nghost;
pgv->jgs[23] = 1-pgv->nghost; pgv->jge[23] = 0;
pgv->kgs[23] = pgv->ijkm_bl[2]+1; pgv->kge[23] = pgv->ijkm_bl[2]+pgv->nghost;
pgv->igs[24] = 1-pgv->nghost; pgv->ige[24] = 0;
pgv->jgs[24] = pgv->ijkm_bl[1]+1; pgv->jge[24] = pgv->ijkm_bl[1]+pgv->nghost;
pgv->kgs[24] = pgv->ijkm_bl[2]+1; pgv->kge[24] = pgv->ijkm_bl[2]+pgv->nghost;
pgv->igs[25] = pgv->ijkm_bl[0]+1; pgv->ige[25] = pgv->ijkm_bl[0]+pgv->nghost;
pgv->jgs[25] = pgv->ijkm_bl[1]+1; pgv->jge[25] = pgv->ijkm_bl[1]+pgv->nghost;
pgv->kgs[25] = pgv->ijkm_bl[2]+1; pgv->kge[25] = pgv->ijkm_bl[2]+pgv->nghost;
pgv->nmax_gc = pgv->gc1_start[25] + 1;
if(pp->ierr==0){
	NinBand.resize(pgv->nbands);
	for(auto i=0; i<pgv->nbands; i++) NinBand[i].resize(pgv->max_bl);
	for(auto i=0; i<pgv->nbands; i++){
		for(auto j=0; j<pgv->max_bl; j++){
			NinBand[i][j] = 0;
		}
	}
}
else
	pp->litError("band_m_init","allocation error for NinBand");

pmo->lit_monitor_create_file_step("lit_band",5);
pmo->lit_monitor_set_header(1, "NinA",'i');
pmo->lit_monitor_set_header(2, "NinT",'i');
pmo->lit_monitor_set_header(3, "NinN",'i');
pmo->lit_monitor_set_header(4, "NinX",'i');
pmo->lit_monitor_set_header(5, "NinZ",'i');
pmo->lit_monitor_create_gnuplot(1,2,5);
};

void band::band_init(){
	std::unique_ptr<global_variable>pgv(new global_variable);
	std::unique_ptr<parallel>pp(new parallel);
	std::unique_ptr<bound>pbou(new bound);
	if(pgv->verbose && pp->myrank==0)
		std::cout<<pgv->clit<<"Starting band_init"<<std::endl;
	band_seed_initG();
	band_grow(true, 'y');
	band_set_block_data();
	band_regenerate_ijk2ic();
	band_cleanup();
	pbou->prepareGhostNodes();
	pbou->updateGhostNodes();
	if(pgv->verbose && pp->myrank ==0)
		std::cout<<pgv->clit<<"Starting band_init ... Done."<<std::endl;
};

void band::band_inject(const int nshapes, shape_t* shapes){
std::array<double, 3> bb_min, bb_max;
std::array<int,    3> sg_min, sg_max, bl_min, bl_max;
int nGhelp;
Gnode_t* Gnhelp;
double Gnew;
block_t* bb;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<bl>pbl(new bl);
std::unique_ptr<gnodes>pgn(new gnodes);
std::unique_ptr<init>pinit(new init);
if(pgv->verbose && pp->myrank==0)
	std::cout<<pgv->clit<<"starting band_inject"<<std::endl;
for(auto is=0; is<nshapes; is++){
	switch(shapes[is].code){
		case 25:
			for(auto i=0; i<3; i++){
				bb_min[i]=shapes[is].rdata[i]-shapes[is].rdata[3]-pgv->G_max;
				bb_max[i]=shapes[is].rdata[i]+shapes[is].rdata[3]+pgv->G_max;
			}
		break;
		default:
			std::cout<<"unknown injected shape code:"<<shapes[is].code<<std::endl;
			pp->parallel_kill(0);
			break;
	};
	for(auto i=0; i<3; i++){
		auto m1 = static_cast<int>((bb_min[i]-pgv->xyzs_sg[i])*pgv->rdxyz_sg[i]);
		auto m2 = static_cast<int>((bb_max[i]-pgv->xyzs_sg[i])*pgv->rdxyz_sg[i]);
		sg_min[i] = std::max(1, m1+1);
		sg_max[i] = std::max(m2+1, pgv->ijkm_sg[i]);
	}
	for(auto isg=sg_min[0]; isg<=sg_max[0]; isg++){
		for(auto jsg=sg_min[1]; jsg<=sg_max[1]; jsg++){
			for(auto ksg=sg_min[2]; ksg<=sg_max[2]; ksg++){
				if(pgv->sg_rank_block[isg-1][jsg-1][ksg-1] == pp->myrank){
					if(pp->ierr==0){
						bb = new block_t;
					}
					else{
						pp->litError("band_inject","allocation error for bb");
					}
					std::vector<int> v = {isg, jsg, ksg};
					pbl->bl_activate_new(bb, v);
					isg = v[0]; jsg = v[1]; ksg = v[2];
					if(pgv->nbl>0) pgv->ibl2bl[pgv->nbl-1].p->next = bb;
					else pgv->bl = bb;
					pgv->nbl += 1;
					auto ibl = pgv->nbl;
					pgv->ibl2bl[ibl-1].p = bb;
					pgv->sg_active[isg-1][jsg-1][ksg-1] = true;
					pgv->sg_rank_block[isg-1][jsg-1][ksg-1] = -ibl;
					for(auto i=0; i<NinBand.size(); i++) NinBand[i][ibl-1] = 0;
					bb->NinA = 0;
					bb->NinT = 0;
					bb->NinN = 0;
					bb->NinW = 0;
					bb->NinX = 0;
					bb->NinZ = 0;
					bb->nG   = 0;
				}
				if(pgv->sg_rank_block[isg-1][jsg-1][ksg-1]<0){
					auto ibl = -pgv->sg_rank_block[isg-1][jsg-1][ksg-1];
					pgv->setBlockPointers(ibl);
					auto f = floor((bb_min[0]-pgv->xyzs_sg[0])*pgv->rdxyz[0]);
					bl_min[0] = std::max(static_cast<double>(pgv->b->imin_), f);

					f = floor((bb_min[1]-pgv->xyzs_sg[1])*pgv->rdxyz[1]);
					bl_min[1] = std::max(static_cast<double>(pgv->b->jmin_), f);

					f = floor((bb_min[2]-pgv->xyzs_sg[2])*pgv->rdxyz[2]);
					bl_min[2] = std::max(static_cast<double>(pgv->b->kmin_), f);

					f = ceil((bb_max[0]-pgv->xyzs_sg[0])*pgv->rdxyz[0]);
					bl_max[0] = std::min(static_cast<double>(pgv->b->imax_), f);

					f = ceil((bb_max[1]-pgv->xyzs_sg[1])*pgv->rdxyz[1]);
					bl_max[1] = std::min(static_cast<double>(pgv->b->jmax_), f);

					f = ceil((bb_max[2]-pgv->xyzs_sg[2])*pgv->rdxyz[2]);
					bl_max[2] = std::min(static_cast<double>(pgv->b->kmax_), f);
					
					nGhelp = 0;

					for(auto i=bl_min[0]-1;i<bl_max[0]; i++){
						for(auto j=bl_min[1]-1; j<bl_max[1]; j++){
							for(auto k=bl_min[2]-1; k<bl_max[2]; k++){
								vector<double> v = {pgv->xc[i+pgv->nghost], pgv->yc[j+pgv->nghost],\
														  pgv->zc[k+pgv->nghost]};
								Gnew = pinit->G_inject(v, shapes+is);
								pgv->xc[i+pgv->nghost] = v[0]; pgv->yc[j+pgv->nghost] = v[1];
								pgv->zc[k+pgv->nghost] = v[2];
								if(Gnew>-pgv->G_max+std::numeric_limits<double>::epsilon()){
									if(pgv->i2c[i+1-pgv->b->imino_][j+1-pgv->b->jmino_]\
										[k+1-pgv->b->kmino_]>0){
										auto ic = pgv->i2c[i+1-pgv->b->imino_][j+1-pgv->b->jmino_]\
																[k+1-pgv->b->kmino_];
										if(Gnew>0.0 && pgv->Gn[ic-1].G>0.0){
											pgv->Gn[ic-1].G = std::min(Gnew, pgv->Gn[ic-1].G);
										}
										else if(Gnew>0.0 && pgv->Gn[ic-1].G<=0.0){
											pgv->Gn[ic-1].G = Gnew;
										}
										else if(Gnew<=0.0 && pgv->Gn[ic-1].G>0.0){
											pgv->Gn[ic-1].G = pgv->Gn[ic-1].G;
										}
										else{
											pgv->Gn[ic-1].G = std::max(Gnew, pgv->Gn[ic-1].G);
										};
									}
									else{
										if(Gnhelp!=nullptr){
											if(pp->ierr==0){
												auto d = pgv->b->ijkm[0]*pgv->b->ijkm[1]*pgv->b->ijkm[2];
												Gnhelp = new Gnode_t[d];
											}
											else{
												pp->litError("band_inject","cannot allocate Gnhelp");
											};
										}
											nGhelp += 1;
											Gnhelp[nGhelp-1].G = Gnew;
											Gnhelp[nGhelp-1].V[0] = 0.0;
											Gnhelp[nGhelp-1].V[1] = 0.0;
											Gnhelp[nGhelp-1].V[2] = 0.0;
											Gnhelp[nGhelp-1].ijk[0] = i+1;
											Gnhelp[nGhelp-1].ijk[1] = j+1;
											Gnhelp[nGhelp-1].ijk[2] = k+1;
											Gnhelp[nGhelp-1].dijk[0] = 0;
											Gnhelp[nGhelp-1].dijk[1] = 0;
											Gnhelp[nGhelp-1].dijk[2] = 0;
									}
								}
							}
						}
					}
					if(nGhelp>0){
						pgn->ensureSizeGnodes(pgv->b, pgv->b->NinZ+nGhelp);
						for(auto ic=pgv->b->NinZ+nGhelp; ic>=pgv->b->NinT+1+nGhelp; ic--)
							pgv->Gn[ic-1] = pgv->Gn[ic-1-nGhelp];
						for(auto ic=pgv->b->NinT+1; ic<=pgv->b->NinT+nGhelp; ic++)
							pgv->Gn[ic-1] = Gnhelp[ic - 1 - pgv->b->NinT];
						pgv->b->NinT = pgv->b->NinT + nGhelp;
						pgv->b->NinN = pgv->b->NinN + nGhelp;
						pgv->b->NinW = pgv->b->NinW + nGhelp;
						pgv->b->NinX = pgv->b->NinX + nGhelp;
						pgv->b->NinZ = pgv->b->NinZ + nGhelp;
						for(auto i=1; i<pgv->nbands; i++){
							pgv->b->NinBand[i] += nGhelp;
						}
						pgv->b->nG += nGhelp;
						for(auto ic=pgv->b->NinT-nGhelp; ic<pgv->b->NinZ; ic++){
							auto d1 = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
							auto d2 = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
							auto d3 = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
							pgv->i2c[d1][d2][d3] = ic+1;
						}
						if(pp->ierr==0) delete [] Gnhelp;
						else pp->litError("band_inject","cannot deallocate Gnhelp");
						for(auto ic=0; ic<pgv->b->NinZ; ic++){
							if(pgv->Gn[ic].ijk[0]<pgv->b->imino_ ||\
							   pgv->Gn[ic].ijk[0]>pgv->b->imaxo_ ||\
							 	pgv->Gn[ic].ijk[1]<pgv->b->jmino_ ||\
								pgv->Gn[ic].ijk[1]>pgv->b->jmaxo_ ||\
								pgv->Gn[ic].ijk[2]<pgv->b->kmino_ ||\
								pgv->Gn[ic].ijk[2]>pgv->b->kmaxo_){
								std::cout<<"post ERROR! i,j,k,ic =",pgv->Gn[ic].ijk[0]<<'\t';
							std::cout<<pgv->Gn[ic].ijk[1]<<'\t'<<pgv->Gn[ic].ijk[2]<<'\t'<<ic<<std::endl;
							std::cout<<"pgv->b->ijkmin,max="<<pgv->b->imino_<<'\t'<<pgv->b->imaxo_<<'\t';
								std::cout<<"NinATNWXYZ="<<pgv->b->NinA<<'\t'<<pgv->b->NinT<<'\t';
								std::cout<<pgv->b->NinN<<'\t'<<pgv->b->NinW<<'\t'<<pgv->b->NinX<<'\t';
								std::cout<<pgv->b->NinZ<<std::endl;
								std::exit(1);
							}
							auto d1 = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
							auto d2 = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
							auto d3 = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
							if(pgv->i2c[d1][d2][d3] != ic+1){
								std::cout<<"ERROR! i,j,k,ic ="<<pgv->Gn[ic].ijk[0]<<'\t'<<pgv->Gn[ic].ijk[1]<<'\t';
								std::cout<<pgv->Gn[ic].ijk[2]<<'\t'<<ic<<std::endl;
								std::exit(1);
							}
						}
						for(auto i=pgv->b->imino_-1; i<pgv->b->imaxo_;i++){
							for(auto j=pgv->b->jmino_-1; j<pgv->b->jmaxo_;j++){
								for(auto k=pgv->b->kmino_-1;k<pgv->b->kmaxo_;k++){
									if(pgv->i2c[i+1-pgv->b->imino_][j+1-pgv->b->jmino_]\
										[k+1-pgv->b->kmino_]>0){
										auto d1 = pgv->i2c[i+1-pgv->b->imino_][j+1-pgv->b->jmino_]\
																[k+1-pgv->b->kmino_];
										if(i+1 != pgv->Gn[d1].ijk[0] ||\
											j+1 != pgv->Gn[d1].ijk[1] ||\
											k+1 != pgv->Gn[d1].ijk[2]){
											std::cout<<"ERROR! i,j,k,c ="<<pgv->Gn[d1].ijk[0]<<'\t';
											std::cout<<pgv->Gn[d1].ijk[1]<<'\t';
											std::cout<<pgv->Gn[d1].ijk[2]<<std::endl;
											std::exit(1);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
};
if(pgv->verbose && pp->myrank==0)
	std::cout<<pgv->clit<<"Starting band_inject ... Done."<<std::endl;
};

void band::band_seed_initG(){
double*** G;
int n;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<gnodes>pgn(new gnodes);
std::unique_ptr<init>pinit(new init);
Gnode_t Gnh;
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	if(pp->myrank==0 && pgv->verbose){
		std::cout<<pgv->clit<<"band init, block # "<<ibl+1<<" of "<<pgv->nbl<<std::endl;
	}
	for(auto i=0; i<pgv->nbands; i++){
		pgv->b->NinBandLayer[i] = 0;
	}
	pgv->b->nG = 0;
	n = 0;
	auto d1  = pgv->b->imax1_-pgv->b->imin1_+1;
	auto d2  = pgv->b->jmax1_-pgv->b->jmin1_+1;
	auto d3  = pgv->b->kmax1_-pgv->b->kmin1_+1;
	auto n11 = pgv->b->imin_-pgv->b->imin1_;
	auto n12 = pgv->b->imax_-pgv->b->imin1_+1;
	auto n21 = pgv->b->jmin_-pgv->b->jmin1_;
	auto n22 = pgv->b->jmax_-pgv->b->jmin1_+1;
	auto n31 = pgv->b->kmin_-pgv->b->kmin1_;
	auto n32 = pgv->b->kmax_-pgv->b->kmin1_+1;
	if(pp->ierr==0){
		G = new double**[d1];
		for(auto i=0; i<d1;i++){
			G[i] = new double* [d2];
			for(auto j=0; j<d2;j++){
				G[i][j] = new double[d3];
			}
		}
	}
	else
		pp->litError("band_seed_initG","cannot allocate G.");

	for(auto i=0; i<d1; i++){
		for(auto j=0; j<d2; j++){
			for(auto k=0; k<d3; k++){
				vector<double> v = {pgv->xc[i+pgv->b->imin1_-(1-pgv->nghost)],\
										  pgv->yc[j+pgv->b->jmin1_-(1-pgv->nghost)],\
										  pgv->zc[k+pgv->b->kmin1_-(1-pgv->nghost)]};
				G[i][j][k] = pinit->G_init_value(v);
			}
		}
	}
	auto n=0;
	for(auto i=n11; i<n12; i++){
		for(auto j=n21; j<n22; j++){
			for(auto k=n31; k<n32; k++){
				if( fabs(G[i][j][k]) < pgv->G_max ) n += 1;
			}
		}
	}
	pgn->ensureSizeGnodes(pgv->b, n);
	for(auto i=0; i<NinBand.size(); i++) NinBand[i][ibl] = 0;
	for(auto i=n11; i<n12; i++){
		for(auto j=n21; j<n22; j++){
			for(auto k=n31; k<n32; k++){
				if(G[i][j][k]*G[i-1][j][k]<=0.0 || G[i][j][k]*G[i+1][j][k]<=0.0 ||\
					G[i][j][k]*G[i][j-1][k]<=0.0 || G[i][j][k]*G[i][j+1][k]<=0.0 ||\
					G[i][j][k]*G[i][j][k-1]<=0.0 || G[i][j][k]*G[i][j][k+1]<=0.0){
					NinBand[0][ibl] += 1;
					Gnh.G = G[i][j][k];
					Gnh.V[0] = 0.0; Gnh.V[1] = 0.0; Gnh.V[2] = 0.0;
					Gnh.ijk[0] = i+pgv->b->imin1_; Gnh.ijk[1] = j+pgv->b->jmin1_;
					Gnh.ijk[2] = k+pgv->b->kmin1_;
					Gnh.dijk[0] = 0; Gnh.dijk[1] = 0; Gnh.dijk[2] = 0;
					pgn->addGnode(pgv->b, Gnh);
				}
			}
		}
	}
	if(pp->ierr==0){
		for(auto i=0; i<d1; i++){
			for(auto j=0; j<d2; j++){
				delete [] G[i][j];
			}
			delete [] G[i];
		}
	delete [] G;
	G = nullptr;
	}
	else pp->litError("band_seed_initG","cannot deallocate G.");
	pgv->b->NinBandLayer[0] = NinBand[0][ibl];
}
bandLayerCounter = 1;
};

void band::band_seed(){
Gnode_t Goh;
bool found;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<gnodes>pgn(new gnodes);
std::unique_ptr<parallel>pp(new parallel);
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	pgv->b->GnodesOld = pgv->b->Gnodes;
	pgv->b->Gnodes = nullptr;
	pgv->b->nG = 0;
	for(auto i=0; i<NinBand.size(); i++) NinBand[i][ibl] = 0;
	for(auto i=0; i<pgv->nbands; i++) pgv->b->NinBandLayer[i] = 0;
	for(auto ic=0; ic<pgv->b->NinZ; ic++){
		pgv->b->GnodesOld[ic].dijk[0]=0;
		pgv->b->GnodesOld[ic].dijk[1]=0;
		pgv->b->GnodesOld[ic].dijk[2]=0;
	}
	for(auto ic = 0; ic < pgv->b->NinT; ic++){
		Goh = pgv->b->GnodesOld[ic];
		auto i = Goh.ijk[0]-pgv->b->imino_; 
		auto j = Goh.ijk[1]-pgv->b->jmino_; 
		auto k = Goh.ijk[2]-pgv->b->kmino_;
		found =false;
		if(pgv->i2c[i+1][j][k]>0){
			if(Goh.G*pgv->b->GnodesOld[pgv->i2c[i+1][j][k]-1].G<=0.0){
				found = true;
			}
		}
		if(!found){
			if(pgv->i2c[i-1][j][k]>0){
				if(Goh.G*pgv->b->GnodesOld[pgv->i2c[i-1][j][k]-1].G<=0.0){
					found = true;
				}
			}
			if(!found){
				if(pgv->i2c[i][j+1][k]>0){
					if(Goh.G*pgv->b->GnodesOld[pgv->i2c[i][j+1][k]-1].G<=0.0){
						found = true;
					}
				}
				if(!found){
					if(pgv->i2c[i][j-1][k]>0){
						if(Goh.G*pgv->b->GnodesOld[pgv->i2c[i][j-1][k]-1].G<=0.0){
							found = true;
						}
					}
					if(!found){
						if(pgv->i2c[i][j][k+1]>0){
							if(Goh.G*pgv->b->GnodesOld[pgv->i2c[i][j][k+1]-1].G<=0.0){
								found = true;
							}
						}
						if(!found){
							if(pgv->i2c[i][j][k-1]>0){
								if(Goh.G*pgv->b->GnodesOld[pgv->i2c[i][j][k-1]-1].G<=0.0){
									found = true;
								}
							}
						}
					}
				}
			}
		}
		if(found){
			NinBand[0][ibl] += 1;
			if(NinBand[0][ibl]==1){
				if(pp->ierr==0) pgv->b->Gnodes = new Gnode_t[pgv->b->nGmax];
				else pp->litError("band_seed","Allocation error for pgv->b->Gnodes in resizeGnodes");
				pgv->Gn = pgv->b->Gnodes;
			}
			pgn->addGnode(pgv->b, Goh);
		}
	}
	pgv->b->NinBandLayer[0] = NinBand[0][ibl];
	if(pgv->b->Gnodes == nullptr){
		if(pp->ierr==0){
			pgv->b->Gnodes = new Gnode_t[pgv->b->nGmax];
		}
		else{
			pp->litError("band_seed","allocation error for pgv->b->Gnodes in resizeGnodes.");
		}
		pgv->Gn = pgv->b->Gnodes;
	}
}
bandLayerCounter = 1;
};

void band::band_grow(bool init1, char ch){
double rband2;
int irband, kkn;
std::array<int, 3> ijk;
Gnode_t Gnh;
bool do_init;
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<gnodes>pgn(new gnodes);
std::unique_ptr<init>pinit(new init);
bool*** person;
if(ch == 'y'){
	do_init = init1;
}
else{
	do_init = false;
}
if(pp->ierr==0){
	skin_s.resize(pgv->max_bl);
	skin_e.resize(pgv->max_bl);
	cloth_s.resize(pgv->max_bl);
	cloth_e.resize(pgv->max_bl);
	nskin.resize(pgv->max_bl);
	ncloth.resize(pgv->max_bl);
	ncount.resize(pgv->max_bl);
	ghostcloth_bl.resize(pgv->max_bl);
}
else
	pp->parallel_die("Allocation error of ghostcloth_bl... in band_grow");
for(auto i=0;i<pgv->max_bl;i++){
	cloth_s[i] = -1;
	cloth_e[i] = -1;
}
for(auto ibl=0; ibl<pgv->max_bl; ibl++){
	for(auto j=0; j<26; j++){
		ghostcloth_bl[ibl].nin_bf[j] = 0;
	}
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	auto d1 = pgv->b->imax1_-pgv->b->imin1_+1;
	auto d2 = pgv->b->jmax1_-pgv->b->jmin1_+1;
	auto d3 = pgv->b->kmax1_-pgv->b->kmin1_+1;
	if(pp->ierr==0){
		pgv->b->person = new bool**[d1];
		for(auto i=0; i<d1; i++){
			pgv->b->person[i] = new bool*[d2];
			for(auto j=0; j<d2; j++){
				pgv->b->person[i][j] = new bool[d3];
			}
		}
	}
	else{
		pp->litError("band_grow","Allocation error of pgv->b->person in band_grow");
	}
	person = pgv->b->person;
	for(auto i=0; i<d1;i++){
		for(auto j=0; j<d2; j++){
			for(auto k=0; k<d3; k++){
				person[i][j][k] = false;
			}
		}
	}
	for(auto ic=0; ic<NinBand[0][ibl]; ic++){
			person[pgv->Gn[ic].ijk[0]-pgv->b->imin1_][pgv->Gn[ic].ijk[1]-pgv->b->jmin1_]\
					[pgv->Gn[ic].ijk[2]-pgv->b->kmin1_] = true;
	}
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	skin_s[ibl] = 1;
	skin_e[ibl] = NinBand[0][ibl];
}
irband = 0;
for(auto ib=1; ib<pgv->nbands; ib++){
	for(auto i=0; i<ncount.size(); i++){
		ncount[i] = 0;
	}
	for(auto is=0; is<pgv->band_size[ib]; is++){
		bandLayerCounter += 1;
		for(auto i=0; i<ncloth.size(); i++){
			ncloth[i] = 0;
		}
		irband += 1;
		rband2 = (static_cast<double>(irband)+0.5)*(static_cast<double>(irband)+0.5);
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			person  = pgv->b->person;
			for(auto isk=skin_s[ibl]-1; isk<skin_e[ibl]; isk++){
				for(auto dk=-1; dk<2; dk++){
					auto kn = pgv->Gn[isk].ijk[2]+dk;
					auto kr = pgv->Gn[isk].dijk[2]+dk;
					auto kr2 = kr*kr;
					for(auto dj=-1; dj<2; dj++){
						auto jn = pgv->Gn[isk].ijk[1]+dj;
						auto jr = pgv->Gn[isk].dijk[1]+dj;
						auto jr2 = jr*jr;
						for(auto di=-1; di<2; di++){
							if(di==0 && dj==0 && dk==0) continue;
							auto in = pgv->Gn[isk].ijk[0]+di;
							auto ir = pgv->Gn[isk].dijk[0]+di;
							if(person[in-pgv->b->imin1_][jn-pgv->b->jmin1_][kn-pgv->b->kmin1_])continue;
							if(kr2+jr2+ir*ir>rband2) continue;
							if(do_init){
								Gnh.ijk[0] = in; Gnh.ijk[1] = jn; Gnh.ijk[2] = kn;
								Gnh.dijk[0] = ir; Gnh.dijk[1] = jr; Gnh.dijk[2] = kr;
								vector<double> v={pgv->xc[in-(1-pgv->nghost)],\
														pgv->yc[jn-(1-pgv->nghost)],\
														pgv->zc[kn-(1-pgv->nghost)]};
								Gnh.G = pinit->G_init_value(v);
								Gnh.V[0] = 0.0; Gnh.V[1] = 0.0; Gnh.V[2] = 0.0;
							}
							else{
								if(pgv->i2c[in-pgv->b->imino_][jn-pgv->b->jmino_][kn-pgv->b->kmino_]>0){
									Gnh = pgv->b->GnodesOld[pgv->i2c[in-pgv->b->imino_][jn-pgv->b->jmino_]\
																  [kn-pgv->b->kmino_]-1];
									Gnh.dijk[0]= ir; Gnh.dijk[1]=jr; Gnh.dijk[2]=kr;
								}	
								else{
									Gnh.ijk[0]=in; Gnh.ijk[1]=jn; Gnh.ijk[2]=kn;
									Gnh.dijk[0]=ir; Gnh.dijk[1]=jr; Gnh.dijk[2]=kr;
									if(pgv->cylindrical && pgv->ijkm_gl[2]>1){
										auto s = (pgv->Gn[isk].G>=0.0) ? 1.0 : -1.0;
										Gnh.G = pgv->Gn[isk].G + s*(abs(di)*pgv->dxyz[0]+\
												  abs(dj)*pgv->dxyz[1]+abs(dk)*pgv->yc[jn-(1-pgv->nghost)]*\
												  pgv->dxyz[2]);
									}
									else{
										double s = (pgv->Gn[isk].G>=0) ? pgv->G_max : -pgv->G_max;
										Gnh.G = s;
									}
									Gnh.V[0] = 0.0; Gnh.V[1] = 0.0; Gnh.V[2] = 0.0;
								}
							}
							person[in-pgv->b->imin1_][jn-pgv->b->jmin1_][kn-pgv->b->kmin1_] = true;
							if(!is_ghostcloth(pgv->b, Gnh, ibl)){
								ncloth[ibl] += 1;
								pgn->addGnode(pgv->b, Gnh);
								if(ncloth[ibl] == 1) cloth_s[ibl] = pgv->b->nG;
							}
							if(pgv->cylindrical && pgv->ijkm_gl[2]>1 && jn==0){
								for(auto kkn=pgv->b->kmin1_-1; kkn<pgv->b->kmax1_; kkn++){
									if(person[in-pgv->b->imin1_][jn-pgv->b->jmin1_][kkn-pgv->b->kmin1_])
										continue;
									if(do_init){
										Gnh.ijk[0]=in; Gnh.ijk[1]=jn; Gnh.ijk[2]=kkn+1;
										Gnh.dijk[0]=ir; Gnh.dijk[1]=jr; Gnh.dijk[2]=kr;
										vector<double> v={pgv->xc[in-(1-pgv->nghost)],\
																pgv->yc[jn-(1-pgv->nghost)],\
												   			pgv->zc[kkn+pgv->nghost]};
										Gnh.G=pinit->G_init_value(v);
										Gnh.V[0]=0.0; Gnh.V[1]=0.0; Gnh.V[2]=0.0;
									}
									else{
										if(pgv->i2c[in-pgv->b->imino_][jn-pgv->b->jmino_]\
													  [kkn+1-pgv->b->kmino_]>0){
											Gnh = pgv->b->GnodesOld[pgv->i2c[in-pgv->b->imino_]\
																		[jn-pgv->b->jmino_][kkn+1-pgv->b->kmino_]-1];
											Gnh.dijk[0]=ir; Gnh.ijk[1]=jr; Gnh.ijk[2]=kr;
										}
										else{
											Gnh.ijk[0]=in; Gnh.ijk[1]=jn; Gnh.ijk[2]=kkn+1;
											Gnh.dijk[0]=ir; Gnh.dijk[1]=jr; Gnh.dijk[2]=kr;
											auto s=(pgv->Gn[isk].G>=0) ? 1.0 : -1.0;
											Gnh.G=pgv->Gn[isk].G+s*(abs(di)*pgv->dxyz[0]+
													abs(dj)*pgv->dxyz[1]+abs(dk)*pgv->yc[jn-(1-pgv->nghost)]*\
													pgv->dxyz[2]);
											Gnh.V[0]=0.0; Gnh.V[1]=0.0; Gnh.V[2]=0.0;
										}
									}
									person[in-pgv->b->imin1_][jn-pgv->b->jmin1_][kn-pgv->b->kmin1_] = true;
									if(!is_ghostcloth(pgv->b, Gnh, ibl)){
										ncloth[ibl] += 1;
										pgn->addGnode(pgv->b, Gnh);
										if(ncloth[ibl]==1) cloth_s[ibl]=pgv->b->nG;
									}
								}
							}
						}
					}
				}
			}
		}
		band_push_ghostcloth(ib);
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->b = pgv->ibl2bl[ibl].p;
			cloth_e[ibl]=pgv->b->nG;
			if(cloth_s[ibl]>0){
				skin_s[ibl] = cloth_s[ibl];
				skin_e[ibl] = cloth_e[ibl];
			}
			ncount[ibl] += ncloth[ibl];
			pgv->b->NinBandLayer[bandLayerCounter-1] = ncloth[ibl];
		}
	}
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		NinBand[ib][ibl] = NinBand[ib-1][ibl]+ncount[ibl];
	}
}
if(pp->ierr==0){
	skin_s.clear(); //skin_s.shrink_to_fit();
	skin_e.clear(); //skin_e.shrink_to_fit();
	cloth_s.clear(); //cloth_s.shrink_to_fit();
	cloth_e.clear(); //cloth_e.shrink_to_fit();
	nskin.clear(); //nskin.shrink_to_fit();
	ncloth.clear(); //ncloth.shrink_to_fit();
	ncount.clear(); //ncount.shrink_to_fit();
	ghostcloth_bl.clear(); //ghostcloth_bl.shrink_to_fit();
}
else
	pp->litError("band_grow","Allocation error of ghostcloth_bl in band_grow");
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->b = pgv->ibl2bl[ibl].p;
	if(pgv->b->GnodesOld != nullptr){
		if(pp->ierr==0){
			delete [] pgv->b->GnodesOld;
			pgv->b->GnodesOld = nullptr;
		}
		else pp->litError("band_grow", "Allocation error of pgv->b->GnodesOld");
	}
	if(pgv->b->person != nullptr){
		auto d1 = pgv->b->imax1_-pgv->b->imin1_+1;
		auto d2 = pgv->b->jmax1_-pgv->b->jmin1_+1;
		auto d3 = pgv->b->kmax1_-pgv->b->kmin1_+1;
		for(auto i=0; i<d1; i++){
			for(auto j=0; j<d2; j++){
				delete [] pgv->b->person[i][j];
			}
			delete[] pgv->b->person[i];
		}
		delete [] pgv->b->person;
		pgv->b->person = nullptr;
	}
	else
		pp->litError("band_grow","Allocation error of pgv->b->person in band_grow");
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->b = pgv->ibl2bl[ibl].p;
	if(pgv->b->nG < static_cast<int>(0.75*static_cast<double>(pgv->b->nGmax)))
		pgn->shortenGnodes(pgv->b);
}
};

bool band::is_ghostcloth(block_t* bb, const Gnode_t &Gnh, const int ibl){
bool add;
int ibf;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
if(Gnh.ijk[0]>=bb->imin_ && Gnh.ijk[0]<=bb->imax_ &&\
	Gnh.ijk[1]>=bb->jmin_ && Gnh.ijk[1]<=bb->jmax_ &&\
	Gnh.ijk[2]>=bb->kmin_ && Gnh.ijk[2]<=bb->kmax_){
		return false;
}
else{
	auto i=Gnh.ijk[0]; auto j=Gnh.ijk[1]; auto k=Gnh.ijk[2];
	return true;
	if(i==bb->imax1_ && j==bb->jmax1_ && k==bb->kmax1_){
		ibf=26;
		goto mylable;
	}
	if(i==bb->imin1_ && j==bb->jmax1_ && k==bb->kmax1_){
		ibf=25;
		goto mylable;
	}
	if(i==bb->imax1_ && j==bb->jmin1_ && k==bb->kmax1_){
		ibf=24;
		goto mylable;
	}
	if(i==bb->imin1_ && j==bb->jmin1_ && k==bb->kmax1_){
		ibf=23;
		goto mylable;
	}
	if(i==bb->imax1_ && j==bb->jmax1_ && k==bb->kmin1_){
		ibf=22;
		goto mylable;
	}
	if(i==bb->imin1_ && j==bb->jmax1_ && k==bb->kmin1_){
		ibf=21;
		goto mylable;
	}
	if(i==bb->imax1_ && j==bb->jmin1_ && k==bb->kmin1_){
		ibf=20;
		goto mylable;
	}
	if(i==bb->imin1_ && j==bb->jmin1_ && k==bb->kmin1_){
		ibf=19;
		goto mylable;
	}
	if(i==bb->imax1_ && j==bb->jmax1_                 ){
		ibf=18;
		goto mylable;
	}
	if(i==bb->imin1_ && j==bb->jmax1_                 ){
		ibf=17;
		goto mylable;
	}
	if(i==bb->imax1_ && j==bb->jmin1_                 ){
		ibf=16;
		goto mylable;
	}
	if(i==bb->imin1_ && j==bb->jmin1_                 ){
		ibf=15;
		goto mylable;
	}
	if(                 j==bb->jmax1_ && k==bb->kmax1_){
		ibf=14;
		goto mylable;
	}
	if(						 j==bb->jmin1_ && k==bb->kmax1_){
		ibf=13;
		goto mylable;
	}
	if(i==bb->imax1_						&& k==bb->kmax1_){
		ibf=12;
		goto mylable;
	}
	if(i==bb->imin1_                  && k==bb->kmax1_){
		ibf=11;
		goto mylable;
	}
	if(						 j==bb->jmax1_ && k==bb->kmin1_){
		ibf=10;
		goto mylable;
	}
	if(						 j==bb->jmin1_ && k==bb->kmin1_){
		ibf=9;
		goto mylable;
	}
	if(i==bb->imax1_                  && k==bb->kmin1_){
		ibf=8;
		goto mylable;
	}
	if(i==bb->imin1_                  && k==bb->kmin1_){
		ibf=7;
		goto mylable;
	}
	if(												k==bb->kmax1_){
		ibf=6;
		goto mylable;
	}
	if(												k==bb->kmin1_){
		ibf=5;
		goto mylable;
	}
	if(						 j==bb->jmax1_						 ){
		ibf=4;
		goto mylable;
	}
	if(						 j==bb->jmin1_   					 ){
		ibf=3;
		goto mylable;
	}
	if(i==bb->imax1_											 ){
		ibf=2;
		goto mylable;
	}
	if(i==bb->imin1_											 ){
		ibf=1;
		goto mylable;
	}
	return true;
//	pp->litError("is_ghostcloth","Should not be here.");
mylable:
	{
		if(ghostcloth_bl[ibl].G_gc==nullptr){
			if(pp->ierr==0){
				ghostcloth_bl[ibl].G_gc = new double[pgv->nmax_gc];
				ghostcloth_bl[ibl].ijk_gc = new int*[3];
				for(auto q1=0; q1<3; q1++){
					ghostcloth_bl[ibl].ijk_gc[q1] = new int[pgv->nmax_gc];
				}
				ghostcloth_bl[ibl].dijk_gc = new int*[3];
				for(auto q1=0; q1<3; q1++){
					ghostcloth_bl[ibl].dijk_gc[q1] = new int[pgv->nmax_gc];
				}
			}
			else
				pp->litError("is_ghostcloth","Allocation error of ghost_cloth G_gc");
		}
		ghostcloth_bl[ibl].nin_bf[ibf-1] += 1;
	 	ghostcloth_bl[ibl].G_gc[pgv->gc1_start[ibf-1]+\
		ghostcloth_bl[ibl].nin_bf[ibf-1]]=Gnh.G;
		for(auto q=0; q<3; q++){
			ghostcloth_bl[ibl].ijk_gc[q][pgv->gc1_start[ibf-1]+ghostcloth_bl[ibl].\
			nin_bf[ibf-1]]=Gnh.ijk[q];
			ghostcloth_bl[ibl].dijk_gc[q][pgv->gc1_start[ibf-1]+ghostcloth_bl[ibl].\
			nin_bf[ibf-1]]=Gnh.dijk[q];
		}
	}
}
};

void band::band_push_ghostcloth(const int ib){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::array<int, 3> ijk_sg, ijk_sg_n, ijk_sg_skip, ijk_sg_s, ijk_sg_e;
int ibl_primary, irank_primary, ibl_secondary, ibf_secondary, irank_secondary;
int i_sg, j_sg, k_sg, id, ipe, nn, proc, procTarget;
int nProcsRecv, nProcsSend, nMessagesSend, sizeComm, selfj, ic, selfSize;
std::vector<int> nNodesPe, iNodePe, sizeMessage, myBuf;
std::vector<MPI_Request> requestsInt, requests;
bool selfComm, skip;
double *sendBuf, *recvBuf;
int **sendBufInt, **recvBufInt;
if(pp->ierr==0){
	nNodesPe.resize(pp->nprocs+1, 0);
	iNodePe.resize(pp->nprocs,0);
	sizeMessage.resize(pp->nprocs, 1);
	myBuf.resize(pp->nprocs, 0);
	requests.resize(pp->nprocs);
	requestsInt.resize(pp->nprocs);
}
else
	pp->parallel_die("allocation error of vNodePe ... in band_push_ghostcloth");
if(pgv->cylindrical && pgv->ijkm_gl[2]>1){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ibf=0; ibf<26; ibf++){
			if(ghostcloth_bl[ibl].nin_bf[ibf]>0){
				calc_partner_block_cyl(pgv->b, ibf, i_sg, j_sg, k_sg, skip);
				if(skip) continue;
				if(pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1]<0){
					myBuf[pp->myrank]==1;
					nNodesPe[pp->myrank+1] += ghostcloth_bl[ibl].nin_bf[ibf];
				}
				else{
					myBuf[pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1]]=1;
					nNodesPe[pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1]+1] += \
					ghostcloth_bl[ibl].nin_bf[ibf];
				}		
			}
		}
	}
}
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ibf=0; ibf<26; ibf++){
			if(ghostcloth_bl[ibl].nin_bf[ibf]>0){

				calc_partner_block(pgv->b, ibf, i_sg, j_sg, k_sg, skip);
				if(skip) continue;
				if(pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1]<0){
					myBuf[pp->myrank]==1;
					nNodesPe[pp->myrank+1] += ghostcloth_bl[ibl].nin_bf[ibf];
				}
				else{
					myBuf[pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1]]=1;
					nNodesPe[pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1]+1] += \
					ghostcloth_bl[ibl].nin_bf[ibf];
				}		
			}
		}
	}
}
pp->ierr = MPI_Reduce_scatter(&myBuf.front(),&nProcsRecv,&sizeMessage.front(),\
									   MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
nProcsSend=0;
for(auto i=0; i<pp->nprocs; i++) nProcsSend += myBuf[i];
for(auto ipe=1; ipe<pp->nprocs+1; ipe++) nNodesPe[ipe] += nNodesPe[ipe-1];
if(pp->ierr==0){
	sendBuf = new double[nNodesPe[pp->nprocs]];
	sendBufInt = new int*[6];
	for(auto i=0; i<6; i++) sendBufInt[i]=new int[nNodesPe[pp->nprocs]];
}
else pp->parallel_die("allocation error of sendBuf in band_push_ghostcloth");
if(pgv->cylindrical && pgv->ijkm_gl[2]>1){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ibf=0; ibf<26; ibf++){
			if(ghostcloth_bl[ibl].nin_bf[ibf]>0){
				calc_partner_block_cyl(pgv->b, ibf, i_sg, j_sg, k_sg, skip);
				if(skip) continue;
				if(pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1]<0) ipe=pp->myrank;
				else ipe = pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1];
				for(auto ic=0; ghostcloth_bl[ibl].nin_bf[ibf]; ic++){
					iNodePe[ipe] += 1;
					sendBuf[iNodePe[ipe]+nNodesPe[ipe]-1]=ghostcloth_bl[ibl].G_gc\
																	  [pgv->gc1_start[ibf]+ic];
					for(auto i=0; i<3; i++){
						sendBufInt[i][iNodePe[ipe]+nNodesPe[ipe]-1]=ghostcloth_bl[ibl].\
										  				ijk_gc[i][pgv->gc1_start[ibf]+ic];
					}
					for(auto i=3; i<6; i++){
						sendBufInt[i][iNodePe[ipe]+nNodesPe[ipe]-1]=ghostcloth_bl[ibl].\
										  				dijk_gc[i-3][pgv->gc1_start[ibf]+ic];
					}
				}
			}
		}
		for(auto i=0; i<26; i++){
			ghostcloth_bl[ibl].nin_bf[i] = 0;
		}
		if(ghostcloth_bl[ibl].G_gc!=nullptr){
			if(pp->ierr==0){
				delete [] ghostcloth_bl[ibl].G_gc;
				if(ghostcloth_bl[ibl].ijk_gc != nullptr){
					for(auto i=0; i<3; i++){
						delete [] ghostcloth_bl[ibl].ijk_gc[i];
					}
					delete [] ghostcloth_bl[ibl].ijk_gc;
				}
				if(ghostcloth_bl[ibl].dijk_gc != nullptr){
					for(auto i=0; i<3; i++){
						delete [] ghostcloth_bl[ibl].dijk_gc[i];
					}
					delete [] ghostcloth_bl[ibl].dijk_gc;
				}
			}
			else
				pp->parallel_die("deallocation error of (ghostcloth_bl) ... in band_push_ghostcloth");
		}
	}
}
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ibf=0; ibf<26; ibf++){
			if(ghostcloth_bl[ibl].nin_bf[ibf]>0){
				calc_partner_block(pgv->b, ibf, i_sg, j_sg, k_sg, skip);
				if(skip) continue;
				if(pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1]<0) ipe=pp->myrank;
				else ipe = pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1];
				for(auto ic=0; ghostcloth_bl[ibl].nin_bf[ibf]; ic++){
					iNodePe[ipe] += 1;
					sendBuf[iNodePe[ipe]+nNodesPe[ipe]-1]=ghostcloth_bl[ibl].G_gc\
																	  [pgv->gc1_start[ibf]+ic];
					for(auto i=0; i<3; i++){
						sendBufInt[i][iNodePe[ipe]+nNodesPe[ipe]-1]=ghostcloth_bl[ibl].\
										  				ijk_gc[i][pgv->gc1_start[ibf]+ic];
					}
					for(auto i=3; i<6; i++){
						sendBufInt[i][iNodePe[ipe]+nNodesPe[ipe]-1]=ghostcloth_bl[ibl].\
										  				dijk_gc[i-3][pgv->gc1_start[ibf]+ic];
					}
				}
			}
		}
		for(auto i=0; i<26; i++){
			ghostcloth_bl[ibl].nin_bf[i] = 0;
		}
		if(ghostcloth_bl[ibl].G_gc!=nullptr){
			if(pp->ierr==0){
				delete [] ghostcloth_bl[ibl].G_gc;
				if(ghostcloth_bl[ibl].ijk_gc != nullptr){
					for(auto i=0; i<3; i++){
						delete [] ghostcloth_bl[ibl].ijk_gc[i];
					}
					delete [] ghostcloth_bl[ibl].ijk_gc;
				}
				if(ghostcloth_bl[ibl].dijk_gc != nullptr){
					for(auto i=0; i<3; i++){
						delete [] ghostcloth_bl[ibl].dijk_gc[i];
					}
					delete [] ghostcloth_bl[ibl].dijk_gc;
				}
			}
			else
				pp->parallel_die("deallocation error of (ghostcloth_bl)");
		}		
	}
}
nMessagesSend = 0;
selfComm = false;
for(auto ipe=0; ipe<pp->nprocs; ipe++){
	if(nNodesPe[ipe+1]>nNodesPe[ipe]){
		sizeComm = nNodesPe[ipe+1]-nNodesPe[ipe];
		auto j = nNodesPe[ipe]+1;
		if(ipe==pp->myrank){
			selfComm = true;
			selfj = j;
			selfSize = sizeComm;
		}
		else{
			nMessagesSend += 1;
			procTarget = ipe;
			pp->ierr = MPI_Isend(&sendBuf[j-1], sizeComm, MPI_DOUBLE, ipe, ipe+1,\
									   MPI_COMM_WORLD, &requestsInt[nMessagesSend]);
			pp->ierr = MPI_Isend(&sendBufInt[0][j-1], 6*sizeComm, MPI_INTEGER, ipe, \
							pp->nprocs+ipe+1, MPI_COMM_WORLD, &requestsInt[nMessagesSend]);
		}
	}
}
auto ii = 0;
if(selfComm){
	auto j = selfj;
	band_glue_cloth(sendBuf, sendBufInt, j, selfSize, ib);
	ii += 1;
}
for(auto nn=ii; nn<nProcsRecv; nn++){
	pp->ierr=MPI_Probe(MPI_ANY_SOURCE, pp->myrank+1, MPI_COMM_WORLD, &pp->status);
	proc = pp->status.MPI_SOURCE;
	pp->ierr=MPI_Get_count(&pp->status, MPI_DOUBLE, &sizeComm);
	if(pp->ierr==0){
		recvBuf = new double[sizeComm];
		recvBufInt = new int*[6];
		for(auto i=0; i<6; i++) recvBufInt[i] = new int[sizeComm];
	}
	else{
		pp->parallel_die("allocation error of recBuf, recvBufInt in band_push_ghost");
	}
	pp->ierr=MPI_Recv(recvBuf,sizeComm,MPI_DOUBLE,proc,pp->myrank+1,MPI_COMM_WORLD,\
							&pp->status);
	pp->ierr=MPI_Recv(&recvBufInt[0][0],6*sizeComm,MPI_INTEGER,proc,pp->myrank+1+\
							pp->nprocs, MPI_COMM_WORLD, &pp->status);
	auto j = 1;
	band_glue_cloth(recvBuf, recvBufInt, j, sizeComm, ib);
	if(pp->ierr==0){
		delete [] recvBuf;	
		for(auto i=0; i<6; i++) delete[]recvBufInt[i];
		delete[] recvBufInt;
	}
	else{
		pp->parallel_die("deallocation error of recvBuf, recvBufInt in band_push_gho");
	}
}
sizeComm = nMessagesSend;
for(auto nn=0; nn<nMessagesSend; nn++){
	pp->ierr=MPI_Waitany(sizeComm,requests.data(),&proc,&pp->status);
	pp->ierr=MPI_Waitany(sizeComm, requestsInt.data(), &proc, &pp->status);
}
if(pp->ierr==0){
	delete [] sendBuf;
	for(auto i=0; i<6; i++)	delete [] sendBufInt[i];
	delete [] sendBufInt;
	nNodesPe.clear();
	iNodePe.clear();
	sizeMessage.clear();
	myBuf.clear();
	requests.clear();
	requestsInt.clear();
}
else
	pp->parallel_die("deallocation error at 1 of sendBuf in band_push_ghostcloth");
	pp->parallel_sync_sg_active();
};

void band::calc_partner_block_cyl(block_t *bb,const int ibf,int &i_sg, int &j_sg,\
											 int &k_sg, bool & skip){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
skip = true;
i_sg = bb->ijk_sg[0]+pgv->dijk_nc[ibf][0];
if(i_sg<1){
	if(pgv->sg_periodic[0]) i_sg += pgv->ijkm_sg[0];
	return;
}
if(i_sg>pgv->ijkm_sg[0]){
	if(pgv->sg_periodic[0]) i_sg -= pgv->ijkm_sg[0];
	else return;
}
j_sg = bb->ijk_sg[1]+pgv->dijk_nc[ibf][1];
if(j_sg<1){
	j_sg = 1;
	if(bb->ijk_sg[2]<=pgv->ijkm_sg[2]/2){
		k_sg = bb->ijk_sg[2]+pgv->ijkm_sg[2]/2+pgv->dijk_nc[ibf][2];
		if(k_sg==pgv->ijkm_sg[2]+1) k_sg=1;
	}
	else{
		k_sg=bb->ijk_sg[2] - pgv->ijkm_sg[2]/2 + pgv->dijk_nc[ibf][2];
		if(k_sg==0) k_sg=pgv->ijkm_sg[2];
	}
	if(pgv->ijkm_sg[2]==1) k_sg = 1;
}
if(j_sg>pgv->ijkm_sg[1]){
	if(pgv->sg_periodic[1]) j_sg -= pgv->ijkm_sg[1];
	else return;
}
if(bb->ijk_sg[1]+pgv->dijk_nc[ibf][1]!= 0){
	k_sg = bb->ijk_sg[2]+pgv->dijk_nc[2][ibf];
	if(k_sg<1){
		if(pgv->sg_periodic[2]) k_sg += pgv->ijkm_sg[2];
		else return;
	}
	if(k_sg>pgv->ijkm_sg[2]){
		if(pgv->sg_periodic[2]) k_sg -= pgv->ijkm_sg[2];
		else return;
	}
}
skip = false;
};

void band::calc_partner_block(block_t* bb, const int ibf, int& i_sg, int& j_sg,\
										int& k_sg, bool& skip){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);

skip = true;
i_sg = bb->ijk_sg[0]+pgv->dijk_nc[ibf][0];
if(i_sg<1){
	if(pgv->sg_periodic[0]) i_sg += pgv->ijkm_sg[0];
	return;
}
if(i_sg>pgv->ijkm_sg[0]){
	if(pgv->sg_periodic[0]) i_sg -= pgv->ijkm_sg[0];
	else return;
}
j_sg = bb->ijk_sg[1] + pgv->dijk_nc[ibf][1];
if(j_sg<1){
	if(pgv->sg_periodic[1]) j_sg += pgv->ijkm_sg[1];
	else return;
}
if(j_sg>pgv->ijkm_sg[1]){
	if(pgv->sg_periodic[1]) j_sg -= pgv->ijkm_sg[1];
	else return;
}
k_sg = bb->ijk_sg[2] + pgv->dijk_nc[ibf][2];
if(k_sg<1){
	if(pgv->sg_periodic[2]) k_sg += pgv->ijkm_sg[2];
	else return;
}
if(k_sg>pgv->ijkm_sg[2]){
	if(pgv->sg_periodic[2]) k_sg -= pgv->ijkm_sg[2];
	else return;
}
skip = false;
};

void band::band_glue_cloth(double *buf, int **bufInt,const int jstart,\
									const int size,const int ib){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<gnodes>pgn(new gnodes);
std::unique_ptr<bl>pbl(new bl);
block_t* bb;
Gnode_t crc;
std::array<int, 3> ijk, ijk0;
bool ***person;
bool new_block;
auto j = jstart;
int ibl;
for(auto ic=j-1; ic<j+size-1; ic++){
	auto i_sg = (bufInt[0][ic]-1)/pgv->ijkm_bl[0]+1;
	auto j_sg = (bufInt[1][ic]-1)/pgv->ijkm_bl[1]+1;
	auto k_sg = (bufInt[2][ic]-1)/pgv->ijkm_bl[2]+1;
	for(auto i=0; i<3; i++) ijk0[i] = bufInt[i][ic];
	if(pgv->sg_periodic[0]){
		if(bufInt[0][ic]){
			i_sg = pgv->ijkm_sg[0];
			bufInt[0][ic] = pgv->ijkm_gl[0];
		}
		else if(bufInt[0][ic]==pgv->ijkm_gl[0]+1){
			i_sg = 1;
			bufInt[0][ic] = 1;
		}
	}
	if(pgv->cylindrical && pgv->ijkm_gl[2]>1){
		if(bufInt[1][ic]==0){
			j_sg = 1;
			bufInt[1][ic] = 1;
			if(bufInt[2][ic]<=pgv->ijkm_gl[2]/2) bufInt[2][ic] += pgv->ijkm_gl[2]/2;
			else bufInt[2][ic] -= pgv->ijkm_gl[2]/2;
			k_sg = (bufInt[2][ic]-1)/pgv->ijkm_bl[2]+1;
			bufInt[4][ic] = -bufInt[4][ic];
		}
	}
	else{
		if(pgv->sg_periodic[1]){
			if(bufInt[1][ic]==0){
				j_sg = pgv->ijkm_sg[1];
				bufInt[1][ic] = pgv->ijkm_gl[1];
			}
			else if(bufInt[1][ic]==pgv->ijkm_gl[1]+1){
				j_sg = 1;
				bufInt[1][ic] = 1;
			}
		}
	}
	if(pgv->sg_periodic[2]){
		if(bufInt[2][ic]==0){
			k_sg = pgv->ijkm_sg[2];
			bufInt[2][ic] = pgv->ijkm_gl[2];
		}
		else if(bufInt[2][ic] == pgv->ijkm_gl[2]+1){
			k_sg = 1;
			bufInt[2][ic] = 1;
		}
	}
	if(!pgv->sg_active[i_sg-1][j_sg-1][k_sg-1]){
		if(pp->ierr==0) bb = new block_t;
		else pp->litError("band_glue_cloth","allocation error for bb");
		std::vector<int> v = {i_sg, j_sg, k_sg};
		pbl->bl_activate_new(bb,v);
		i_sg = v[0]; j_sg = v[1]; k_sg = v[2];
		if(pgv->nbl>0) pgv->ibl2bl[pgv->nbl-1].p->next = bb;
		else pgv->bl = bb;
		auto d1 = bb->imax1_ - bb->imin1_+1;
		auto d2 = bb->jmax1_ - bb->jmin1_+1;
		auto d3 = bb->kmax1_ - bb->kmin1_+1;
		if(pp->ierr==0){
			bb->person = new bool**[d1];
			for(auto i=0; i<d1; i++){
				bb->person[i] = new bool*[d2];
				for(auto j=0; j<d2; j++){
					bb->person[i][j] = new bool[d3];
				}
			}
		}
		else{
			pp->litError("band_glue_cloth","deallocation error for bb->person");
		}
		person = bb->person;
		for(auto i=0; i<d1; i++){
			for(auto j=0; j<d2; j++){
				for(auto k=0; k<d3; k++){
					person[i][j][k] = false;
				}
			}
		}
		pgv->nbl += 1;
		ibl = pgv->nbl-1;
		pgv->ibl2bl[ibl].p = bb;
		pgv->sg_active[i_sg-1][j_sg-1][k_sg-1] = true;
		pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1] = -pgv->nbl;
		skin_s[ibl] = 0;
		skin_e[ibl] = 0;
		for(auto i=0;i<NinBand.size();i++) NinBand[i][ibl]=0;
		new_block = true;
	}
	else{
		ibl = -pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1]-1;
		bb = pgv->ibl2bl[ibl].p;
		person = bb->person;
		bool check = true;
		for(auto i=0; i<NinBand.size(); i++){
			if(NinBand[i][ibl] != 0){
				check = false;
				break;
			}
		}
		if(check && ncloth[ibl] == 0) new_block = true;
		else new_block = false;
	}
	pgv->setBlockPointers(ibl);
	ijk[0] = bufInt[0][ic]; ijk[1] = bufInt[1][ic]; ijk[2] = bufInt[2][ic];
	if(person[ijk[0]-bb->imin1_][ijk[1]-bb->jmin1_][ijk[2]-bb->kmin1_]) continue;
	person[ijk[0]-bb->imin1_][ijk[1]-bb->jmin1_][ijk[2]-bb->kmin1_] = true;
	if(pgv->i2c[ijk[0]-pgv->b->imino_][ijk[1]-pgv->b->jmino_][ijk[2]-pgv->b->kmino_]){
		crc = pgv->b->GnodesOld[pgv->i2c[ijk[0]-pgv->b->imino_][ijk[1]-pgv->b->jmino_]\
									  [ijk[2]-pgv->b->kmino_]-1];
		crc.dijk[0] = bufInt[3][ic];
		crc.dijk[2] = bufInt[4][ic];
		crc.dijk[3] = bufInt[5][ic];
	}
	else{
		crc.G = buf[ic];
		crc.V[0] = 0.0; crc.V[1] = 0.0; crc.V[2] = 0.0;
		crc.ijk[0] = ijk[0]; crc.ijk[1] = ijk[1]; crc.ijk[2] = ijk[2];
		for(auto i=0; i<3; i++) crc.dijk[i] = bufInt[i-3][ic];
	}
	pgn->addGnode(pgv->b, crc);
	if(new_block){
		cloth_s[ibl] = pgv->b->nG;
		cloth_e[ibl] = pgv->b->nG;
		new_block = false;
	}
	ncloth[ibl] += 1;
	if(ncloth[ibl] == 1) cloth_s[ibl] = pgv->b->nG;
}
};

void band::band_set_block_data(){
	std::unique_ptr<global_variable> pgv(new global_variable);	
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto i=0;i<pgv->nbands; i++){
			pgv->b->NinBand[i]=NinBand[i][ibl];
		}
		pgv->b->NinA = NinBand[0][ibl];
		pgv->b->NinT = NinBand[1][ibl];
		pgv->b->NinN = NinBand[2][ibl];
		pgv->b->NinW = NinBand[3][ibl];
		pgv->b->NinX = NinBand[4][ibl];
		pgv->b->NinZ = pgv->b->NinX;
	}
};

void band::band_regenerate_ijk2ic(){
	std::unique_ptr<global_variable> pgv(new global_variable);
	size_t d1, d2, d3;
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);

		d1 = pgv->b->imaxo_ - pgv->b->imino_ + 1;
		d2 = pgv->b->jmaxo_ - pgv->b->jmino_ + 1;
		d3 = pgv->b->kmaxo_ - pgv->b->kmino_ + 1;
		for(auto i=0; i<d1; i++){
			for(auto j=0; j<d2; j++){
				for(auto k=0; k<d3; k++){
					pgv->i2c[i][j][k] = 0;
				}
			}
		}

		for(auto ic=0; ic<pgv->b->NinX; ic++){
			pgv->i2c[pgv->Gn[ic].ijk[0]-pgv->b->imino_][pgv->Gn[ic].ijk[1]-pgv->b->jmino_]\
					  	[pgv->Gn[ic].ijk[2]-pgv->b->kmino_] = ic+1;

		}
	}
};

void band::band_cleanup(){
	for(auto i=0; i<NinBand.size(); i++){
		for(auto j=0; j<NinBand[i].size(); j++){
	 		NinBand[i][j] = 0;
		}
	}
};

void band::band_m_cleanup(){
	std::unique_ptr<parallel>pp(new parallel);
	if(pp->ierr==0){
		for(auto i=0; i<NinBand.size(); i++){
			NinBand[i].clear();
			NinBand[i].shrink_to_fit();
		}
		NinBand.clear();
		NinBand.shrink_to_fit();
	}
	else pp->litError("band_m_cleanup","deallocation error of Ninband");
};

void band::band_monitor(){
	std::unique_ptr<parallel> pp(new parallel);
	std::unique_ptr<monitor>pmo(new monitor);
	auto NinA = pp->nodesInBand('A');
	auto NinT = pp->nodesInBand('T');
	auto NinN = pp->nodesInBand('N');
	auto NinX = pp->nodesInBand('X');
	auto NinZ = pp->nodesInBand('Z');
	pmo->lit_monitor_select_file("lit_band");
	auto v1 = static_cast<double>(NinA);
	auto v2 = static_cast<double>(NinT);
	auto v3 = static_cast<double>(NinN);
	auto v4 = static_cast<double>(NinX);
	auto v5 = static_cast<double>(NinZ);
	pmo->lit_monitor_set_single_values(v1, v2, v3, v4, v5);
};

std::vector<std::vector<int>> band::NinBand;
int band::bandLayerCounter = 0;
std::vector<ghostcloth_bl_t> band::ghostcloth_bl;
std::vector<int> band::skin_s;
std::vector<int> band::skin_e;
std::vector<int> band::cloth_s;
std::vector<int> band::cloth_e;
std::vector<int> band::nskin;
std::vector<int> band::ncloth;
std::vector<int> band::ncount;

