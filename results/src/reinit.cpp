//Written by Salar Safarkhani

#include"reinit.h"
#include<limits>
#include<cmath>
#include<cfloat>
void reinit::lit_reinit_pde(){
double *S, *fijk, *rijk; // checked
double ** LG, **dpG, **dmdpG;
vector<vector<double>> dmG, dmdpmG, dmdppG;
std::array<double, 3> grad_G_c, grad_G_m, grad_G_p;
double base, rm, rp, vm, vp, cylrad;
bool is_done;
int iter;
vector<int> inR, inA;
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<monitor>pmo(new monitor);
std::unique_ptr<litBuffer>plbuf(new litBuffer);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<weno>pweno(new weno);
std::unique_ptr<litparam>ppar(new litparam);
std::unique_ptr<timing>pti(new timing);
triggerGradMin = ppar->get_real_param("LIT.REINIT_TRIGGER_MIN",1.0e-4,2);
triggerGradMax = ppar->get_real_param("LIT.REINIT_TRIGGER_MAX",2.0,2);
bool fixG0Newton = ppar->get_logical_param("LIT_FIX_GO_NEWTON",false,2);
bool hocr2 = ppar->get_logical_param("LIT_REINIT_HOCR2",false,2);
if(max_iter_reinit == 0){
	if(pp->myrank==0){
		std::cout<<"WARNING:Skipping reinit but still clipping"<<std::endl;
		lit_clip_G();
		pmo->lit_monitor_select_file("lit_reinit");
		pmo->lit_monitor_set_single_values(pmo->monitor_gradG_min,\
													  pmo->monitor_gradG_max);
		
		pmo->lit_monitor_dump_values_iter(2,0);
	}
	return;
}
if(pgv->verbose && pp->myrank==0)
	std::cout<<pgv->clit<<"Starting lit_reinit"<<std::endl; // sth should be added // checked
pti->lit_timing_start("reinit");
if(fixG0Newton)
	pp->litError("lit_reinit_pde","removed fixG0Newton code. out it back");
if(!pgv->cylindrical){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic = pgv->b->NinN; ic<pgv->b->NinX; ic++){
			auto sign = (pgv->Gn[ic].G>=0) ? 1.0: -1.0;
			pgv->Gn[ic].G= sign*pgv->G_max;
		}
	}
}
size_t S_S;
S = plbuf->getR1Buffer(S_S);
size_t fijk_S;
fijk = plbuf->getR1Buffer(fijk_S);
for(auto i=0; i<fijk_S; i++) fijk[i] = 0.0;
size_t LG_S1, LG_S2;
LG = plbuf->getR2Buffer(LG_S1, LG_S2);
size_t dpG_S1, dpG_S2, dmdpG_S1, dmdpG_S2;
if(ppar->trim(pgv->schemeReinit) == "WENO-3"||ppar->trim(pgv->schemeReinit)=="WENO-5"){
	dpG = plbuf->getR2Buffer(dpG_S1, dpG_S2);
	dmdpG = plbuf->getR2Buffer(dmdpG_S1, dmdpG_S2);
}

calcSignFunction(S);
pmo->lit_monitor_select_file("lit_reinit");
pmo->lit_monitor_set_single_values(pmo->monitor_gradG_min,pmo->monitor_gradG_max);
pmo->lit_monitor_dump_values_iter(2,0);
if(hocr2){
	//allocate memory; allocated; checked
	size_t rijk_S;
	rijk = plbuf->getR1Buffer(rijk_S);
	calcHOCR2rijk(rijk);
}
for(iter=0; iter<max_iter_reinit; iter++){
	for(auto irk=0; irk<RK_reinit; irk++){
		if(pgv->schemeReinit == "WENO-3" || pgv->schemeReinit == "WENO-5")
			pweno->prepare_WENO_ghost(dpG_S1, dpG_S2, dpG,dmdpG_S1, dmdpG_S2, dmdpG);
		if(hocr2) calcHOCR2fijk(fijk, rijk);
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			if(pgv->schemeReinit == "UPWIND-1"){
				if(pgv->cylindrical){
					for(auto ic=0; ic<pgv->b->NinN; ic++){
						auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
						auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
						auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
						grad_G_p[0]=(pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[ic].G)*\
										 pgv->rdxyz[0];
						grad_G_p[1]=(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[ic].G)*\
										 pgv->rdxyz[1];
						grad_G_p[2]=(pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[ic].G)*\
										 pgv->rdxyz[2]*pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
						grad_G_m[0]=(pgv->Gn[pgv->i2c[i-1][j][k]-1].G-pgv->Gn[ic].G)*\
										 pgv->rdxyz[0];
						grad_G_m[1]=(pgv->Gn[pgv->i2c[i][j-1][k]-1].G-pgv->Gn[ic].G)*\
										 pgv->rdxyz[1];
						grad_G_m[2]=(pgv->Gn[pgv->i2c[i][j][k-1]-1].G-pgv->Gn[ic].G)*\
										 pgv->rdxyz[2]*pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
//						LG[irk][pgv->ioff+ic]=-(std::max(S[pgv->ioff+ic],0.0)*(sqrt(\
						std::max(pow(std::max(grad_G_m[0],0.0),2),pow(std::min(grad_G_p[1],\
						0.0),2))+std::max(pow(std::max(grad_G_m[1],0.0),2),pow(std::min(\
						grad_G_p[1],0.0),2))+std::max(pow(std::max(grad_G_m[2],0.0),2),\
						pow(std::min(grad_G_p[2],0.0),2)))-1.0)+std::min(S[pgv->ioff+ic],\
						0.0)*(sqrt(std::max(pow(std::min(grad_G_m[0],0.0),2),pow(std::max(\
						grad_G_p[0],0.0),2))+std::max(pow(std::min(grad_G_m[1],0.0),2),\
						pow(std::max(grad_G_p[1],0.0),2))+std::max(pow(std::min(grad_G_m\
						[2],0.0),2),pow(std::max(grad_G_p[2],0.0),2))-1.0)));
						double v1 = std::max(S[pgv->ioff+ic], 0.0);
						double v2 = std::max(grad_G_m[0], 0.0);
						double v3 = std::min(grad_G_p[0], 0.0);
						double v4 = std::max(pow(v2,2), pow(v3,2));
						double v5 = std::max(grad_G_m[1], 0.0);
						double v6 = std::min(grad_G_p[1], 0.0);
						double v7 = std::max(pow(v5,2), pow(v6,2));
						double v8 = std::max(grad_G_m[2], 0.0);
						double v9 = std::min(grad_G_p[2], 0.0);
						double v10 = std::max(pow(v8,2), pow(v9,2));

						double v11 = std::min(S[pgv->ioff+ic], 0.0);
						double v12 = std::min(grad_G_m[0], 0.0);
						double v13 = std::max(grad_G_p[0], 0.0);
						double v14 = std::max(pow(v12, 2), pow(v13, 2));
						double v15 = std::min(grad_G_m[1], 0.0);
						double v16 = std::max(grad_G_p[1], 0.0);
						double v17 = std::max(pow(v15, 2), pow(v16, 2));
						double v18 = std::min(grad_G_m[2], 0.0);
						double v19 = std::max(grad_G_p[2], 0.0);
						double v20 = std::max(pow(v18, 2), pow(v19, 2));

						LG[irk][pgv->ioff+ic]=-(v1*(sqrt(v4+v7+v10)-1.0)+ v11*sqrt(v14+v17+v20)-1.0);
					}
				}
				else{
					for(auto ic=0; ic<pgv->b->NinN; ic++){
						auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
						auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
						auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
						grad_G_p[0]=(pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[ic].G)*\
										 pgv->rdxyz[0];
						grad_G_p[1]=(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[ic].G)*\
										 pgv->rdxyz[1];
						grad_G_p[2]=(pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[ic].G)*\
										 pgv->rdxyz[2];
						grad_G_m[0]=(pgv->Gn[pgv->i2c[i-1][j][k]-1].G-pgv->Gn[ic].G)*\
										 pgv->rdxyz[0];
						grad_G_m[1]=(pgv->Gn[pgv->i2c[i][j-1][k]-1].G-pgv->Gn[ic].G)*\
										 pgv->rdxyz[1];
						grad_G_m[2]=(pgv->Gn[pgv->i2c[i][j][k-1]-1].G-pgv->Gn[ic].G)*\
										 pgv->rdxyz[2];
						LG[irk][pgv->ioff+ic]=-(std::max(S[pgv->ioff+ic],0.0)*(sqrt(\
						std::max(pow(std::max(grad_G_m[0],0.0),2),pow(std::min(grad_G_p[1],\
						0.0),2))+std::max(pow(std::max(grad_G_m[1],0.0),2),pow(std::min(\
						grad_G_p[1],0.0),2))+std::max(pow(std::max(grad_G_m[2],0.0),2),\
						pow(std::min(grad_G_p[2],0.0),2)))-1.0)+std::min(S[pgv->ioff+ic],\
						0.0)*(sqrt(std::max(pow(std::min(grad_G_m[0],0.0),2),pow(std::max(\
						grad_G_p[0],0.0),2))+std::max(pow(std::min(grad_G_m[1],0.0),2),\
						pow(std::max(grad_G_p[1],0.0),2))+std::max(pow(std::min(grad_G_m\
						[2],0.0),2),pow(std::max(grad_G_p[2],0.0),2))-1.0)));
					}
				}
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					for(auto irki=0;irki<=irk; irki++){
						pgv->Gn[ic].G = pgv->Gn[ic].G+alphaRK[irki][irk]*LG[irki][pgv->ioff+ic];
					}
				}
			}
			else if(pgv->schemeReinit=="WENO-3"){
				if(pgv->cylindrical){
					for(auto ic = 0; ic<pgv->b->NinN; ic++){
						auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
						auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
						auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
						if(ic+1>pgv->b->NinT){
							grad_G_p[0] = dpG[0][pgv->ioff+ic];
							grad_G_p[1] = dpG[1][pgv->ioff+ic];
							grad_G_p[2] = dpG[2][pgv->ioff+ic]*pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
							grad_G_m[0] = dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1];
							grad_G_m[1] = dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1];
							grad_G_m[2] = dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]*\
							pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
						}
						else{
							base = 0.5*(dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]+dpG[0][\
							pgv->ioff+ic]);
							grad_G_p[0] = base-pweno->wWENO3(dmdpG[0][pgv->ioff+ic], dmdpG[0][\
							pgv->ioff+pgv->i2c[i+1][j][k]-1])*(dmdpG[0][pgv->ioff+pgv->i2c\
							[i+1][j][k]-1]-dmdpG[0][pgv->ioff+ic]);
							grad_G_m[0] = base-pweno->wWENO3(dmdpG[0][pgv->ioff+ic], dmdpG[0][\
							pgv->ioff+pgv->i2c[i-1][j][k]-1])*(dmdpG[0][pgv->ioff+ic]-dmdpG[0]\
							[pgv->ioff+pgv->i2c[i-1][j][k]-1]);
							base = 0.5*(dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]+dpG[1][\
							pgv->ioff+ic]);
							grad_G_p[1] = base-pweno->wWENO3(dmdpG[1][pgv->ioff+ic],dmdpG[1][\
							pgv->ioff+pgv->i2c[i][j+1][k]-1])*(dmdpG[1][pgv->ioff+pgv->i2c\
							[i][j+1][k]-1]-dmdpG[1][pgv->ioff+ic]);
							grad_G_m[1] = base-pweno->wWENO3(dmdpG[1][pgv->ioff+ic], dmdpG[1][\
							pgv->ioff+pgv->i2c[i][j-1][k]-1])*(dmdpG[0][pgv->ioff+ic]-dmdpG[\
							1][pgv->ioff+pgv->i2c[i][j-1][k]-1]);
							base = 0.5*(dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]+dpG[2][\
							pgv->ioff+ic]);
							grad_G_p[2] = (base-pweno->wWENO3(dmdpG[2][pgv->ioff+ic], dmdpG[2][\
							pgv->ioff+pgv->i2c[i][j][k+1]-1])*(dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1]-\
							dmdpG[2][pgv->ioff+ic]))*pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
							grad_G_m[2] = (base-pweno->wWENO3(dmdpG[2][pgv->ioff+ic], dmdpG[2][\
							pgv->ioff+pgv->i2c[i][j][k-1]-1])*(dmdpG[2][pgv->ioff+ic]-dmdpG[2]\
							[pgv->ioff+pgv->i2c[i][j][k-2]-1]))*pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
						}
						LG[irk][pgv->ioff+ic]=-(std::max(S[pgv->ioff+ic],0.0)*(sqrt(\
						std::max(pow(std::max(grad_G_m[0],0.0),2),pow(std::min(grad_G_p[1],\
						0.0),2))+std::max(pow(std::max(grad_G_m[1],0.0),2),pow(std::min(\
						grad_G_p[1],0.0),2))+std::max(pow(std::max(grad_G_m[2],0.0),2),\
						pow(std::min(grad_G_p[2],0.0),2)))-1.0)+std::min(S[pgv->ioff+ic],\
						0.0)*(sqrt(std::max(pow(std::min(grad_G_m[0],0.0),2),pow(std::max(\
						grad_G_p[0],0.0),2))+std::max(pow(std::min(grad_G_m[1],0.0),2),\
						pow(std::max(grad_G_p[1],0.0),2))+std::max(pow(std::min(grad_G_m\
						[2],0.0),2),pow(std::max(grad_G_p[2],0.0),2))-1.0)));
					}
				}
				else{
					for(auto ic=0; ic<pgv->b->NinN; ic++){
						auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
						auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
						auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
						if(ic+1>pgv->b->NinT){
							for(auto q=0; q<3; q++) grad_G_p[q]=dpG[q][pgv->ioff+ic];
							grad_G_m[0] = dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1];
							grad_G_m[1] = dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1];
							grad_G_m[2] = dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1];
						}
						else{
							base = 0.5*(dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]+dpG[0][pgv->ioff+ic]);
							grad_G_p[0] = base-pweno->wWENO3(dmdpG[0][pgv->ioff+ic], dmdpG[0][\
							pgv->ioff+pgv->i2c[i+1][j][k]-1])*(dmdpG[0][pgv->ioff+pgv->i2c\
							[i+1][j][k]-1]-dmdpG[0][pgv->ioff+ic]);
							grad_G_m[0] = base-pweno->wWENO3(dmdpG[0][pgv->ioff+ic], dmdpG[0][\
							pgv->ioff+pgv->i2c[i-1][j][k]-1])*(dmdpG[0][pgv->ioff+ic]-dmdpG[\
							0][pgv->ioff+pgv->i2c[i-1][j][k]-1]);
							base = 0.5*(dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]+dpG[1][\
							pgv->ioff+ic]);
							grad_G_p[1] = base-pweno->wWENO3(dmdpG[1][pgv->ioff+ic], dmdpG[1][\
							pgv->ioff+pgv->i2c[i][j+1][k]-1])*(dmdpG[1][pgv->ioff+pgv->i2c\
							[i][j+1][k]-1]-dmdpG[1][pgv->ioff+ic]);
							grad_G_m[1] = base-pweno->wWENO3(dmdpG[1][pgv->ioff+ic], dmdpG[1][\
							pgv->ioff+pgv->i2c[i][j-1][k]-1])*(dmdpG[0][pgv->ioff+ic]-dmdpG[\
							1][pgv->ioff+pgv->i2c[i][j-1][k]-1]);
							base = 0.5*(dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]+dpG[2][\
							pgv->ioff+ic]);
							grad_G_p[2] = base-pweno->wWENO3(dmdpG[2][pgv->ioff+ic], dmdpG[2][\
							pgv->ioff+pgv->i2c[i][j][k+1]-1])*(dmdpG[2][pgv->ioff+pgv->i2c\
							[i][j][k+1]-1]-dmdpG[2][pgv->ioff+ic]);
							grad_G_m[2] = base-pweno->wWENO3(dmdpG[2][pgv->ioff+ic], dmdpG[2][\
							pgv->ioff+pgv->i2c[i][j][k-1]-1])*(dmdpG[2][pgv->ioff+ic]-dmdpG[\
							2][pgv->ioff+pgv->i2c[i][j][k-1]-1]);
						}
						LG[irk][pgv->ioff+ic]=-(std::max(S[pgv->ioff+ic],0.0)*(sqrt(\
						std::max(pow(std::max(grad_G_m[0],0.0),2),pow(std::min(grad_G_p[1],\
						0.0),2))+std::max(pow(std::max(grad_G_m[1],0.0),2),pow(std::min(\
						grad_G_p[1],0.0),2))+std::max(pow(std::max(grad_G_m[2],0.0),2),\
						pow(std::min(grad_G_p[2],0.0),2)))-1.0)+std::min(S[pgv->ioff+ic],\
						0.0)*(sqrt(std::max(pow(std::min(grad_G_m[0],0.0),2),pow(std::max(\
						grad_G_p[0],0.0),2))+std::max(pow(std::min(grad_G_m[1],0.0),2),\
						pow(std::max(grad_G_p[1],0.0),2))+std::max(pow(std::min(grad_G_m\
						[2],0.0),2),pow(std::max(grad_G_p[2],0.0),2))-1.0)));
					}
				}
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					for(auto irki=0; irki<=irk; irki++){
						pgv->Gn[ic].G += alphaRK[irki][irk]*LG[irki][pgv->ioff+ic];
					}
				}
			}
			else if(pgv->schemeReinit=="WENO-5"){
				if(pgv->cylindrical){
					for(auto ic=0; ic<pgv->b->NinN; ic++){
						auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
						auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
						auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
						if(ic+1>pgv->b->NinT){
							grad_G_p[0] = dpG[0][pgv->ioff+ic];
							grad_G_p[1] = dpG[1][pgv->ioff+ic];
							grad_G_p[2] = dpG[2][pgv->ioff+ic]*pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
							grad_G_m[0] = dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1];
							grad_G_m[1] = dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1];
							grad_G_m[2] = dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]*\
											  pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
						}
						else{
							base = r112*(-dpG[0][pgv->ioff+pgv->i2c[i-2][j][k]-1]+7.0*dpG\
 						   [0][pgv->ioff+pgv->i2c[i-1][j][k]-1]+7.0*dpG[0][pgv->\
							ioff+ic]-dpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1]);
							grad_G_p[0]=base+pweno->G_WENO_5th(dmdpG[0][pgv->ioff+pgv->i2c\
							[i+1][j][k]-1],dmdpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1],dmdpG\
							[0][pgv->ioff+ic],dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]);
							grad_G_m[0]=base-pweno->G_WENO_5th(dmdpG[0][pgv->ioff+pgv->i2c\
							[i-2][j][k]-1],dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1],dmdpG\
							[0][pgv->ioff+ic],dmdpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1]);
							base = r112*(-dpG[1][pgv->ioff+pgv->i2c[i][j-2][k]-1]+7.0*dpG\
 						   [1][pgv->ioff+pgv->i2c[i][j-1][k]-1]+7.0*dpG[1][pgv->\
							ioff+ic]-dpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1]);
							grad_G_p[1]=base+pweno->G_WENO_5th(dmdpG[1][pgv->ioff+pgv->i2c\
							[i][j+1][k]-1],dmdpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1],dmdpG\
							[1][pgv->ioff+ic],dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]);
							grad_G_m[1]=base-pweno->G_WENO_5th(dmdpG[1][pgv->ioff+pgv->i2c\
							[i][j-2][k]-1],dmdpG[1][pgv->ioff+pgv->i2c[i][j-2][k]-1],dmdpG\
							[1][pgv->ioff+ic],dmdpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1]);

							base = r112*(-dpG[2][pgv->ioff+pgv->i2c[i][j][k-2]-1]+7.0*dpG\
 						   [2][pgv->ioff+pgv->i2c[i][j][k-1]-1]+7.0*dpG[2][pgv->\
							ioff+ic]-dpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1]);
							grad_G_p[2]=(base+pweno->G_WENO_5th(dmdpG[2][pgv->ioff+pgv->i2c\
							[i][j][k+1]-1],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1],dmdpG\
							[2][pgv->ioff+ic],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]))*\
							pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
							grad_G_m[2]=(base-pweno->G_WENO_5th(dmdpG[2][pgv->ioff+pgv->i2c\
							[i][j][k-2]-1],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-2]-1],dmdpG\
							[2][pgv->ioff+ic],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1]))*\
							pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
						}
						LG[irk][pgv->ioff+ic]=-(std::max(S[pgv->ioff+ic],0.0)*(sqrt(\
						std::max(pow(std::max(grad_G_m[0],0.0),2),pow(std::min(grad_G_p[1],\
						0.0),2))+std::max(pow(std::max(grad_G_m[1],0.0),2),pow(std::min(\
						grad_G_p[1],0.0),2))+std::max(pow(std::max(grad_G_m[2],0.0),2),\
						pow(std::min(grad_G_p[2],0.0),2)))-1.0)+std::min(S[pgv->ioff+ic],\
						0.0)*(sqrt(std::max(pow(std::min(grad_G_m[0],0.0),2),pow(std::max(\
						grad_G_p[0],0.0),2))+std::max(pow(std::min(grad_G_m[1],0.0),2),\
						pow(std::max(grad_G_p[1],0.0),2))+std::max(pow(std::min(grad_G_m\
						[2],0.0),2),pow(std::max(grad_G_p[2],0.0),2))-1.0)))-hocr2w*\
						fijk[pgv->ioff+ic];
					}
				}
				else{
					for(auto ic=0; ic<pgv->b->NinN; ic++){
						auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
						auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
						auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
						if(ic+1>pgv->b->NinT){
							for(auto q=0; q<3; q++) grad_G_p[q]=dpG[q][pgv->ioff+ic];
							grad_G_m[0] = dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1];
							grad_G_m[1] = dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1];
							grad_G_m[2] = dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1];
						}
						else{
							base = r112*(-dpG[0][pgv->ioff+pgv->i2c[i-2][j][k]-1]+7.0*dpG\
 						   [0][pgv->ioff+pgv->i2c[i-1][j][k]-1]+7.0*dpG[0][pgv->\
							ioff+ic]-dpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1]);
							
							grad_G_p[0]=base+pweno->G_WENO_5th(dmdpG[0][pgv->ioff+pgv->i2c\
							[i+2][j][k]-1],dmdpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1],dmdpG\
							[0][pgv->ioff+ic],dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]);

							grad_G_m[0]=base-pweno->G_WENO_5th(dmdpG[0][pgv->ioff+pgv->i2c\
							[i-2][j][k]-1],dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1],dmdpG\
							[0][pgv->ioff+ic],dmdpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1]);
							

							base = r112*(-dpG[1][pgv->ioff+pgv->i2c[i][j-2][k]-1]+7.0*dpG\
 						   [1][pgv->ioff+pgv->i2c[i][j-1][k]-1]+7.0*dpG[1][pgv->\
							ioff+ic]-dpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1]);

							grad_G_p[1]=base+pweno->G_WENO_5th(dmdpG[1][pgv->ioff+pgv->i2c\
							[i][j+2][k]-1],dmdpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1],dmdpG\
							[1][pgv->ioff+ic],dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]);

							grad_G_m[1]=base-pweno->G_WENO_5th(dmdpG[1][pgv->ioff+pgv->i2c\
							[i][j-2][k]-1],dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1],dmdpG\
							[1][pgv->ioff+ic],dmdpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1]);


							base = r112*(-dpG[2][pgv->ioff+pgv->i2c[i][j][k-2]-1]+7.0*dpG\
 						   [2][pgv->ioff+pgv->i2c[i][j][k-1]-1]+7.0*dpG[2][pgv->\
							ioff+ic]-dpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1]);

							grad_G_p[2]=base+pweno->G_WENO_5th(dmdpG[2][pgv->ioff+pgv->i2c\
							[i][j][k+2]-1],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1],dmdpG\
							[2][pgv->ioff+ic],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]);

							grad_G_m[2]=base-pweno->G_WENO_5th(dmdpG[2][pgv->ioff+pgv->i2c\
							[i][j][k-2]-1],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1],dmdpG\
							[2][pgv->ioff+ic],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1]);
						}

						double v1 = std::max(S[pgv->ioff+ic], 0.0);
						double v2 = std::max(grad_G_m[0], 0.0);
						double v3 = std::min(grad_G_p[0], 0.0);
						double v4 = std::max(pow(v2,2), pow(v3,2));
						double v5 = std::max(grad_G_m[1], 0.0);
						double v6 = std::min(grad_G_p[1], 0.0);
						double v7 = std::max(pow(v5,2), pow(v6,2));
						double v8 = std::max(grad_G_m[2], 0.0);
						double v9 = std::min(grad_G_p[2], 0.0);
						double v10 = std::max(pow(v8,2), pow(v9,2));

						double v11 = std::min(S[pgv->ioff+ic], 0.0);
						double v12 = std::min(grad_G_m[0], 0.0);
						double v13 = std::max(grad_G_p[0], 0.0);
						double v14 = std::max(pow(v12, 2), pow(v13, 2));
						double v15 = std::min(grad_G_m[1], 0.0);
						double v16 = std::max(grad_G_p[1], 0.0);
						double v17 = std::max(pow(v15, 2), pow(v16, 2));
						double v18 = std::min(grad_G_m[2], 0.0);
						double v19 = std::max(grad_G_p[2], 0.0);
						double v20 = std::max(pow(v18, 2), pow(v19, 2));

						LG[irk][pgv->ioff+ic]=-(v1*(sqrt(v4+v7+v10)-1.0)+ v11*(sqrt(v14+v17+v20)-1.0))-\
														hocr2w*fijk[pgv->ioff+ic];
					}
				}
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					for(auto irki=0; irki<=irk; irki++){
						pgv->Gn[ic].G += alphaRK[irki][irk]*LG[irki][pgv->ioff+ic];
					}
				}
			}
		}
		if(!(iter==max_iter_reinit && irk==RK_reinit)){
			lit_clip_G();
			pbou->updateGhostNodes();
		}
	}
	is_done = !(trigger_reinit(iter));
	pmo->lit_monitor_select_file("lit_reinit");
	pmo->lit_monitor_set_single_values(pmo->monitor_gradG_min,pmo->monitor_gradG_max);
	pmo->lit_monitor_dump_values_iter(2,iter);
	if(is_done)	break;
}
lit_clip_G();
pbou->updateGhostNodes();
plbuf->freeR1Buffer(S);
plbuf->freeR2Buffer(LG_S1, LG_S2, LG);
plbuf->freeR1Buffer(fijk);
if(hocr2) plbuf->freeR1Buffer(rijk);
if(ppar->trim(pgv->schemeReinit) == "WENO-3" || ppar->trim(pgv->schemeReinit)=="WENO-5"){
	plbuf->freeR2Buffer(dpG_S1, dpG_S2, dpG);
	plbuf->freeR2Buffer(dmdpG_S1, dmdpG_S2, dmdpG);
}
pmo->monitor_reinit_steps += std::min(iter+1,max_iter_reinit);
pp->ierr = MPI_Barrier(MPI_COMM_WORLD);
pti->lit_timing_stop("reinit");
if(pgv->verbose && pp->myrank==0){
	std::cout<<pgv->clit<<"Starting lit_reinit ... Done in # steps ";
	std::cout<<std::min(iter,max_iter_reinit)<<std::endl;
}
};

bool reinit::trigger_reinit(const int iter){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<litparam> ppar(new litparam);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<monitor> pmo(new monitor);
double myGradMax, myGradMin, gradMin, gradMax, gradG;
bool trigger;
min_iter_reinit = ppar->get_integer_param("LIT.REINIT_STEPS_MIN",0, 2);
if(iter>=min_iter_reinit){
	myGradMin = DBL_MAX;
	myGradMax = -DBL_MAX;
	if(pgv->cylindrical){
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinT; ic++){
				auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
				auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
				auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
				gradG=pow((pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[pgv->i2c[i-1][j][k]].\
						G)*pgv->rdxyz[0],2)+pow((pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn\
						[pgv->i2c[i][j-1][k]].G)*pgv->rdxyz[1],2)+pow((pgv->Gn[pgv->i2c\
						[i][j][k+1]-1].G-pgv->Gn[pgv->i2c[i][j][k-1]-1].G)*pgv->rdxyz[2]*\
						pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)],2);
				myGradMin = std::min(myGradMin, 0.25*gradG);
				myGradMax = std::max(myGradMax, 0.25*gradG);
			}
		}
	}
	else{
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinT; ic++){
				auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
				auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
				auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
				gradG=pow((pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[pgv->i2c[i-1][j][k]-1].G)*\
							  pgv->rdxyz[0],2)+\
						pow((pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[pgv->i2c[i][j-1][k]-1].G)*\
							  pgv->rdxyz[1],2)+\
						pow((pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[pgv->i2c[i][j][k-1]-1].G)*\
							  pgv->rdxyz[2],2);
				myGradMin = std::min(myGradMin, 0.25*gradG);
				myGradMax = std::max(myGradMax, 0.25*gradG);
			}
		}
	}
	pp->parallel_all_min(myGradMin, gradMin);
	pp->parallel_all_max(myGradMax, gradMax);
	if(gradMax < -1.0){
		pmo->monitor_gradG_max = -0.01;
		pmo->monitor_gradG_min = -0.01;
		return false;
	}
	if(sqrt(gradMin)<triggerGradMin || sqrt(gradMax)>triggerGradMax) trigger=true;
	else trigger = false;
	if(pgv->verbose && pp->myrank == 0){
		std::cout<<pgv->clit<<"|grad G| min,max = "<<sqrt(gradMin)<<'\t'<<sqrt(gradMax)<<std::endl;
	}
	pmo->monitor_gradG_max = sqrt(gradMax);
	pmo->monitor_gradG_min = sqrt(gradMin);
}
else trigger = true;
return trigger;
};

void reinit::reinit_m_init(){
auto i0 = 0;
auto i1 = 1;
auto i2 = 2;
auto i3 = 3;
std::unique_ptr<timing>pti(new timing);
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<litparam>ppar(new litparam);
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<monitor>pmo(new monitor);
pti->lit_timing_create("reinit");
pti->lit_timing_start("reinit");
auto dt_reinit = pgv->CFL_reinit*pgv->dxyz_min;
auto s = pgv->band_size[0]+pgv->band_size[1]+pgv->band_size[2];
auto i = ceil(static_cast<double>(s)/pgv->CFL_reinit);
max_iter_reinit = ppar->get_integer_param("LIT.REINIT_STEPS_MAX",i,2);
min_iter_reinit = ppar->get_integer_param("LIT.REINIZ_STEPS_MIN",0,2);
if(pgv->verbose && pp->myrank==0){
	if(i != max_iter_reinit){
		std::cout<<pgv->clit<<"WARNING: setting number of reinitialization steps from, to";
	}
	else{
		std::cout<<pgv->clit<<"Setting number of reinitialization steps to ";
		std::cout<<max_iter_reinit<<std::endl;
	}
}
RK_reinit = ppar->get_integer_param("LIT.REINIT_RK",3,2);
int WENO_reinit = ppar->get_integer_param("LIT_REINIT_WENO",5,2);
if(pp->ierr==0){
	alphaRK.resize(RK_reinit);
	for(auto q=0; q<RK_reinit; q++) alphaRK[q].resize(RK_reinit);
}
else pp->litError("reinit_m_init","allocation error for alphaRK");
if(RK_reinit==1){
	alphaRK[0][0] = 1.0;
}
else if(RK_reinit == 3){
	alphaRK[i0][i1-1]=  1.0;	alphaRK[i1][i1-1]=0.0;    alphaRK[i2][i1-1]=0.0;
	alphaRK[i0][i2-1]= -0.75;  alphaRK[i1][i2-1]=0.25;   alphaRK[i2][i2-1]=0.0;
	alphaRK[i0][i3-1]= -r112;  alphaRK[i1][i3-1]=-r112;  alphaRK[i2][i3-1]=r23;
}
else{
	std::cout<<pgv->clit<<"ERROR unknown RK-order in reinit"<<std::endl;
	pp->parallel_kill(0);
}
for(auto m=0; m<alphaRK.size(); m++){
	for(auto n=1; n<=RK_reinit; n++){
		alphaRK[m][n-1] *= dt_reinit;
	}
}
triggerGradMin = ppar->get_real_param("LIT.REINIT_TRIGGER_MIN",1.0e-4,2);
triggerGradMax = ppar->get_real_param("LIT.REINIT_TRIGGER_MAX",2.0,2);
pmo->lit_monitor_create_file_iter("lit_reinit",2,2);
pmo->lit_monitor_set_header(1,"min |grad G|",'r');
pmo->lit_monitor_set_header(2,"max |grad G|",'r');
pmo->lit_monitor_create_gnuplot(1,2);
pti->lit_timing_stop("reinit");
};

void reinit::lit_clip_G(){
std::unique_ptr<global_variable>pgv(new global_variable);
double m;
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinW; ic++){
		m = std::max(pgv->G_min, pgv->Gn[ic].G);
		pgv->Gn[ic].G = std::min(pgv->G_max, m);
	}
	for(auto ic=pgv->b->NinW; ic<pgv->b->NinX; ic++){
		m = std::max(pgv->G_min, pgv->Gn[ic].G);
		pgv->Gn[ic].G = std::min(pgv->G_max,m);
	}
}
};

void reinit::calcSignFunction(double *S){
std::unique_ptr<global_variable>pgv(new global_variable);
std::array<double, 3> grad_G_c;
if(pgv->cylindrical){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic =0; ic<pgv->b->NinT; ic++){
			auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
			auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
			auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
			grad_G_c[0]=0.5*(pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[pgv->i2c[i-1][j]\
																			[k]-1].G)*pgv->rdxyz[0];
			grad_G_c[1]=0.5*(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[pgv->i2c[i][j-1]\
																			[k]-1].G)*pgv->rdxyz[1];
			grad_G_c[2]=0.5*(pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[pgv->i2c[i][j]\
										[k+1]-1].G)*pgv->rdxyz[2]*pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
			S[pgv->ioff+ic]=pgv->Gn[ic].G/sqrt(pow(pgv->Gn[ic].G,2)+(pow(grad_G_c[0],2)\
			+pow(grad_G_c[1],2)+pow(grad_G_c[2],2))*pgv->dxyz_min_2);
		}
		for(auto ic=pgv->b->NinT; ic<pgv->b->NinN; ic++){
			auto sign = (pgv->Gn[ic].G >= 0.0) ? 1.0: -1.0;
			S[pgv->ioff+ic] = sign;
		}
	}
}
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinT; ic++){
			auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
			auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
			auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
				grad_G_c[0]=0.5*(pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[pgv->i2c[i-1][j]\
																				[k]-1].G)*pgv->rdxyz[0];
				grad_G_c[1]=0.5*(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[pgv->i2c[i][j-1]\
																				[k]-1].G)*pgv->rdxyz[1];
				grad_G_c[2]=0.5*(pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[pgv->i2c[i][j]\
															[k-1]-1].G)*pgv->rdxyz[2];
			S[pgv->ioff+ic]=pgv->Gn[ic].G/sqrt(pow(pgv->Gn[ic].G,2)+(pow(grad_G_c[0],2)\
			+pow(grad_G_c[1],2)+pow(grad_G_c[2],2))*pgv->dxyz_min_2);
		}
		for(auto ic=pgv->b->NinT; ic<pgv->b->NinN; ic++){
			double sign = (pgv->Gn[ic].G >= 0.0) ? 1.0 : -1.0;
			S[pgv->ioff+ic] = sign;
		}
	}
}
};

void reinit::calcHOCR2rijk(double *rijk){
double phi_sum;
std::unique_ptr<global_variable> pgv(new global_variable);
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinA; ic++){
		auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
		phi_sum = 0.0;
		if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i-1][j][k]-1].G < 0.0)
			phi_sum += pgv->Gn[pgv->i2c[i-1][j][k]-1].G;
		if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i+1][j][k]-1].G < 0.0)
			phi_sum += pgv->Gn[pgv->i2c[i+1][j][k]-1].G;
		if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i][j-1][k]-1].G < 0.0)
			phi_sum += pgv->Gn[pgv->i2c[i][j-1][k]-1].G;
		if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i][j+1][k]-1].G < 0.0)
			phi_sum += pgv->Gn[pgv->i2c[i][j+1][k]-1].G;
		if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i][j][k-1]-1].G < 0.0)
			phi_sum += pgv->Gn[pgv->i2c[i][j][k-1]-1].G;
		if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i][j][k+1]-1].G < 0.0)
			phi_sum += pgv->Gn[pgv->i2c[i][j][k+1]-1].G;
		rijk[pgv->ioff+ic] = pgv->Gn[ic].G/phi_sum;
	}
}
};

void reinit::calcHOCR2fijk(double *fijk, double *rijk){
double phi_sum;
int hasNeighbor;
std::unique_ptr<global_variable>pgv(new global_variable);
for(auto ibl=0; ibl< pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinA; ic++){
		phi_sum = 0.0;
		hasNeighbor = 0;
		int i = pgv->Gn[ic].G-pgv->b->imino_;
		int j = pgv->Gn[ic].G-pgv->b->jmino_;
		int k = pgv->Gn[ic].G-pgv->b->kmino_;
		if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i-1][j][k]-1].G < 0.0){
			phi_sum += pgv->Gn[pgv->i2c[i-1][j][k]-1].G;
			hasNeighbor = 1;
		}
		if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i+1][j][k]-1].G < 0.0){
			phi_sum += pgv->Gn[pgv->i2c[i+1][j][k]-1].G;
			hasNeighbor = 1;
		}
		if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i][j-1][k]-1].G < 0.0){
			phi_sum += pgv->Gn[pgv->i2c[i][j-1][k]-1].G;
			hasNeighbor = 1;
		}
		if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i][j+1][k]-1].G < 0.0){
			phi_sum += pgv->Gn[pgv->i2c[i][j+1][k]-1].G;
			hasNeighbor = 1;
		}
		if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i][j][k-1]-1].G < 0.0){
			phi_sum += pgv->Gn[pgv->i2c[i][j][k-1]-1].G;
			hasNeighbor = 1;
		}
		if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i][j][k+1]-1].G < 0.0){
			phi_sum += pgv->Gn[pgv->i2c[i][j][k+1]-1].G;
			hasNeighbor = 1;
		}
		if(hasNeighbor == 1)
			fijk[pgv->ioff+ic]=(rijk[pgv->ioff+ic]*phi_sum-pgv->Gn[ic].G)/pgv->dxyz_min;
		else fijk[pgv->ioff+ic] = 0.0;
	}
}
};
vector<vector<double>> reinit::alphaRK(0.0); // check it
int reinit::RK_reinit = 3;
int reinit::max_iter_reinit = 0;
double reinit::triggerGradMin = 0.0;
double reinit::triggerGradMax = 0.0;
int reinit::min_iter_reinit = 0;
