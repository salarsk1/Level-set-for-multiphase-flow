//Written by Salar Safarkhani

#include"weno.h"
#include<cmath>
/*void weno::prepare_WENO_5th(double **dpG, double** dmdpG,double** dmdpmG,
									 double **dmdppG, bool TbandOnly){
bool onlyTband;
onlyTband = TbandOnly;
if(onlyTband && pgv->band_size[2]<3)
	pp->litError("prepare_WENO_5th","Nband size must be at least 3 for Tband only..");
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinN; ic++){
		auto i = pgv->Gn[ic].ijk[0]-1;
		auto j = pgv->Gn[ic].ijk[1]-1;
		auto k = pgv->Gn[ic].ijk[2]-1;
		dpG[0][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i+1][j][k]].G-pgv->Gn[ic].G)*\
										  pgv->rdxyz[0];
		dpG[1][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i][j][k+1]].G-pgv->Gn[ic].G)*\
										  pgv->rdxyz[1];
		dpG[2][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i][j][k+1]].G-pgv->Gn[ic].G)*\
										  pgv->rdxyz[2];
	}
	for(auto ic = pgv->b->NinN; ic<pgv->b->NinW; ic++){
		auto i = pgv->Gn[ic].ijk[0]-1;
		auto j = pgv->Gn[ic].ijk[1]-1;
		auto k = pgv->Gn[ic].ijk[2]-1;
		if(pgv->i2c[i+1][j][k]>0){
			dpG[0][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i+1][j][k]].G-pgv->Gn[ic].G*\
											  pgv->rdxyz[0]);
		}
		else dpG[0][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i][j+1][k]>0){
			dpG[1][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i][j+1][k]].G-pgv->Gn[ic].G*\
											  pgv->rdxyz[1]);
		}
		else dpG[1][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i][j][k+1]>0){
			dpG[2][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i][j][k+1]].G-pgv->Gn[ic].G*\
											  pgv->rdxyz[2]);
		}
		else dpG[2][pgv->ioff+ic] = 0.0;
	}
	// I do not know what is the start point of dpG;
	for(auto i=0; i<sizeof(dpG)/sizeof(dpG[0]); i++){
		for(auto j=pgv->ioff+pgv->b->NinW; j<pgv->ioff+pgv->b->NinZ;j++){
			dpG[i][j]=0.0;
		}
	}
}
pbou->updateGhostR2(dpG, pbou->S_DP);
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinN; ic++){
		auto i = pgv->Gn[ic].ijk[0]-1;
		auto j = pgv->Gn[ic].ijk[1]-1;
		auto k = pgv->Gn[ic].ijk[2]-1;
		dmdpG[0][pgv->ioff+ic] = dpG[0][pgv->ioff+ic]-dpG[0][pgv->ioff+pgv->i2c\
																			 [i-1][j][k]];
		dmdpG[1][pgv->ioff+ic] = dpG[1][pgv->ioff+ic]-dpG[1][pgv->ioff+pgv->i2c\
																			 [i][j-1][k]];
		dmdpG[2][pgv->ioff+ic] = dpG[2][pgv->ioff+ic]-dpG[2][pgv->ioff+pgv->i2c\
																			 [i][j][k-1]];
	}
	for(auto ic=pgv->b->NinN;ic<pgv->b->NinW;ic++){
		auto i = pgv->Gn[ic].ijk[0]-1;
		auto j = pgv->Gn[ic].ijk[1]-1;
		auto k = pgv->Gn[ic].ijk[2]-1;
		if(pgv->i2c[i-1][j][k]>0){
			dmdpG[0][pgv->ioff+ic]=dpG[0][pgv->ioff+ic]-dpG[0][pgv->ioff+pgv->i2c\
																			  [i-1][j][k]];
		}
		else dmdpG[0][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i][j-1][k]){
			dmdpG[1][pgv->ioff+ic]=dpG[1][pgv->ioff+ic]-dpG[1][pgv->ioff+pgv->i2c\
																			  [i][j-1][k]];
		}
		else dmdpG[1][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i][j][k-1]){
			dmdpG[2][pgv->ioff+ic]=dpG[2][pgv->ioff+ic]-dpG[2][pgv->ioff+pgv->i2c\
																			  [i][j][k-1]];
		}
		else dmdpG[2][pgv->ioff+ic] = 0.0;
	}
//the followinf for loop could be written in an other way
	for(auto i=0; i<sizeof(dmdpG)/sizeof(dmdpG[0]); i++){
		for(auto j=pgv->ioff+pgv->b->NinW; j<pgv->ioff+pgv->b->NinZ;j++){
			dmdpG[i][j]=0.0;
		}
	}
}
pbou->updateGhostR2(dmdpG);
for(auto ibl=0; pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinN; ic++){
		auto i = pgv->Gn[ic].ijk[0]-1;
		auto j = pgv->Gn[ic].ijk[1]-1;
		auto k = pgv->Gn[ic].ijk[2]-1;
		dmdpmG[0][pgv->ioff+ic] = dmdpG[0][pgv->ioff+ic]-dmdpG[0][pgv->ioff+pgv->i2c\
																					[i-1][j][k]];
		dmdpmG[1][pgv->ioff+ic] = dmdpG[1][pgv->ioff+ic]-dmdpG[1][pgv->ioff+pgv->i2c\
																					[i][j-1][k]];
		dmdpmG[2][pgv->ioff+ic] = dmdpG[2][pgv->ioff+ic]-dmdpG[2][pgv->ioff+pgv->i2c\
																					[i][j][k-1]];
		dmdppG[0][pgv->ioff+ic] = dmdpG[0][pgv->ioff+pgv->i2c[i+1][j][k]]-dmdpG\
																				 [0][pgv->ioff+ic];
		dmdppG[1][pgv->ioff+ic] = dmdpG[1][pgv->ioff+pgv->i2c[i][j+1][k]]-dmdpG\
																				 [1][pgv->ioff+ic];
		dmdppG[2][pgv->ioff+ic] = dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+1]]-dmdpG\
																				 [2][pgv->ioff+ic];
	}
	for(auto ic=pgv->b->NinN; ic<pgv->b->NinW; ic++){
		auto i = pgv->Gn[ic].ijk[0]-1;
		auto j = pgv->Gn[ic].ijk[1]-1;
		auto k = pgv->Gn[ic].ijk[2]-1;
		if(pgv->i2c[i-1][j][k]>0){
			dmdpmG[0][pgv->ioff+ic] = dmdpG[0][pgv->ioff+ic]-dmdpG[0][pgv->ioff+\
											  pgv->i2c[i-1][j][k]];
		}
		else dmdpmG[0][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i][j-1][k]>0){
			dmdpmG[1][pgv->ioff+ic] = dmdpG[1][pgv->ioff+ic]-dmdpG[1][pgv->ioff+\
											  pgv->i2c[i][j-1][k]];
		}
		else dmdpmG[1][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i][j][k-1]>0){
			dmdpmG[2][pgv->ioff+ic] = dmdpG[2][pgv->ioff+ic]-dmdpG[2][pgv->ioff+\
											  pgv->i2c[i][j][k-1]];
		}
		else dmdpmG[2][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i+1][j][k]){
			dmdppG[0][pgv->ioff+ic] = dmdpG[0][pgv->ioff+pgv->i2c[i+1][j][k]]-\
											  dmdpG[0][pgv->ioff+ic];
		}
		else dmdppG[0][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i][j+1][k]){
			dmdppG[1][pgv->ioff+ic] = dmdpG[1][pgv->ioff+pgv->i2c[i][j+1][k]]-\
											  dmdpG[1][pgv->ioff+ic];
		}
		else dmdppG[1][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i][j][k+1]){
			dmdppG[2][pgv->ioff+ic] = dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+1]]-\
											  dmdpG[2][pgv->ioff+ic];
		}
		else dmdppG[2][pgv->ioff+ic] = 0.0;
	}
	for(auto i=0; i<sizeof(dmdpmG)/sizeof(dmdpmG[0]); i++){
		for(auto j=pgv->b->NinW; j<pgv->ioff+pgv->b->NinZ; j++){
			dmdpmG[i][j] = 0.0;
		}
	}
	for(auto i=0; i<sizeof(dmdppG)/sizeof(dmdppG[0]); i++){
		for(auto j=pgv->b->NinW; j<pgv->ioff+pgv->b->NinZ; j++){
			dmdppG[i][j] = 0.0;
		}
	}
}
pbou->updateGhostR2(dmdpmG, pbou->S_DM);
pbou->updateGhostR2(dmdppG, pbou->S_DP);
};*/

void weno::prepare_WENO_ghost(size_t dpG_S1,size_t dpG_S2,double **dpG,size_t dmdpG_S1, size_t dmdpG_S2, double **dmdpG){
std::unique_ptr<global_variable> pgv(new global_variable);
for(auto i=0; i<dpG_S1; i++){
	for(auto j=0; j<dpG_S2; j++){
		dpG[i][j] = 0.0;
	}
}
for(auto i=0; i<dmdpG_S1; i++){
	for(auto j=0; j<dmdpG_S2; j++){
		dmdpG[i][j] = 0.0;
	}
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinN; ic++){
		auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
		dpG[0][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[0];
		dpG[1][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[1];
		dpG[2][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[2];
	}
	for(auto ic=pgv->b->NinN; ic<pgv->b->NinW; ic++){
		auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
		if(pgv->i2c[i+1][j][k]>0)
			dpG[0][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[0];
		if(pgv->i2c[i][j+1][k]>0)
			dpG[1][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[1];	
		if(pgv->i2c[i][j][k+1]>0)
			dpG[2][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[2];
	}
	for(auto ic=pgv->b->NinX; ic<pgv->b->NinZ; ic++){
		auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;

		if(i+pgv->b->imino_ < pgv->b->imaxo_){
			if(pgv->i2c[i+1][j][k]>0){
				dpG[0][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[0];
			}
		}
		if(j+pgv->b->jmino_ < pgv->b->jmaxo_){
			if(pgv->i2c[i][j+1][k]>0)
				dpG[1][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[1];
		}
		if(k+pgv->b->kmino_ < pgv->b->kmaxo_){
			if(pgv->i2c[i][j][k+1]>0)
				dpG[2][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[2];
		}
	}
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinN; ic++){
		auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
		dmdpG[0][pgv->ioff+ic]=dpG[0][pgv->ioff+ic]-dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1];
		dmdpG[1][pgv->ioff+ic]=dpG[1][pgv->ioff+ic]-dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1];
		dmdpG[2][pgv->ioff+ic]=dpG[2][pgv->ioff+ic]-dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1];
	}
	for(auto ic=pgv->b->NinN; ic<pgv->b->NinW; ic++){
		auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
		if(pgv->i2c[i-1][j][k]>0)
			dmdpG[0][pgv->ioff+ic]=dpG[0][pgv->ioff+ic]-dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1];
		if(pgv->i2c[i][j-1][k]>0)
			dmdpG[1][pgv->ioff+ic]=dpG[1][pgv->ioff+ic]-dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1];
		if(pgv->i2c[i][j][k-1]>0)
			dmdpG[2][pgv->ioff+ic]=dpG[2][pgv->ioff+ic]-dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1];
	}
	for(auto ic=pgv->b->NinX; ic<pgv->b->NinZ; ic++){
		auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
		if(i>0){
			if(pgv->i2c[i-1][j][k]>0)
				dmdpG[0][pgv->ioff+ic]=dpG[0][pgv->ioff+ic]-dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1];
		}
		if(j>0){
			if(pgv->i2c[i][j-1][k]>0)
				dmdpG[1][pgv->ioff+ic]=dpG[1][pgv->ioff+ic]-dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1];
		}
		if(k>0){
			if(pgv->i2c[i][j][k-1]>0)
				dmdpG[2][pgv->ioff+ic]=dpG[2][pgv->ioff+ic]-dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1];
		}
	}
}
};

double weno::G_WENO_5th(const double a,const double b,const double c,const double d){
double w0, w2, a0, a1, a2, IS0, IS1, IS2, rh;
IS0 = 13.0 * pow(a-b, 2) + 3.0 * pow(a-3.0*b, 2);
IS1 = 13.0 * pow(b-c, 2) + 3.0 * pow(b+c, 2);
IS2 = 13.0 * pow(c-d, 2) + 3.0 * pow(3.0*c-d, 2);
a0 = pow(WENOepsilon+IS0, -2);
a1 = 6.0 * pow(WENOepsilon+IS1, -2);
a2 = 3.0 * pow(WENOepsilon+IS2, -2);
rh = 1.0/(a0+a1+a2);
w0 = a0 * rh;
w2 = a2 * rh;
return (r13 * w0 * (a-2.0*b+c)+r16*(w2-0.5)*(b-2.0*c+d));
};

/*void weno::prepare_WENO_3rd(double **dpG,double**dmG){
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinN; ic++){
		auto i = pgv->Gn[ic].ijk[0]-1;
		auto j = pgv->Gn[ic].ijk[1]-1;
		auto k = pgv->Gn[ic].ijk[2]-1;
		dpG[0][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i+1][j][k]].G-pgv->Gn[ic].G)*\
																			pgv->rdxyz[0];
		dpG[1][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i][j+1][k]].G-pgv->Gn[ic].G)*\
																			pgv->rdxyz[1];
		dpG[2][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i][j][k+1]].G-pgv->Gn[ic].G)*\
																			pgv->rdxyz[2];
		dmG[0][pgv->ioff+ic] = (pgv->Gn[ic].G-pgv->Gn[pgv->i2c[i-1][j][k]].G)*\
																			pgv->rdxyz[0];
		dmG[1][pgv->ioff+ic] = (pgv->Gn[ic].G-pgv->Gn[pgv->i2c[i][j-1][k]].G)*\
																			pgv->rdxyz[1];
		dmG[2][pgv->ioff+ic] = (pgv->Gn[ic].G-pgv->Gn[pgv->i2c[i][j][k-1]].G)*\
																			pgv->rdxyz[2];
	}
	for(auto ic=pgv->b->NinN; ic<pgv->b->NinW; ic++){
		auto i = pgv->Gn[ic].ijk[0]-1;
		auto j = pgv->Gn[ic].ijk[1]-1;
		auto k = pgv->Gn[ic].ijk[2]-1;
		if(pgv->i2c[i+1][j][k]>0)
			dpG[0][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i+1][j][k]].G-pgv->Gn[ic].G)*\
																				pgv->rdxyz[0];
		else dpG[0][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i][j+1][k]>0)
			dpG[1][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i][j+1][k]].G-pgv->Gn[ic].G)*\
																				pgv->rdxyz[1];
		else dpG[1][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i][j][k+1]>0)
			dpG[2][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i][j][k+1]].G-pgv->Gn[ic].G)*\
																				pgv->rdxyz[2];
		else dpG[2][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i-1][j][k]>0)
			dmG[0][pgv->ioff+ic]=(pgv->Gn[ic].G-pgv->Gn[pgv->i2c[i-1][j][k]].G)*\
																				pgv->rdxyz[0];
		else dmG[0][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i][j-1][k]>0)
			dmG[1][pgv->ioff+ic]=(pgv->Gn[ic].G-pgv->Gn[pgv->i2c[i][j-1][k]].G)*\
																				pgv->rdxyz[1];
		else dmG[1][pgv->ioff+ic] = 0.0;
		if(pgv->i2c[i][j][k-1]>0)
			dmG[2][pgv->ioff+ic]=(pgv->Gn[ic].G-pgv->Gn[pgv->i2c[i][j][k-1]].G)*\
																				pgv->rdxyz[2];
		else dmG[2][pgv->ioff+ic] = 0.0;
	}
	for(auto i=0; i<sizeof(dpG)/sizeof(dpG[0]); i++){
		for(auto j=pgv->ioff+pgv->b->NinW; j<sizeof(dpG[i])/sizeof(double); j++){
			dpG[i][j] = 0.0;
		}
	}
	for(auto i=0; i<sizeof(dmG)/sizeof(dmG[0]); i++){
		for(auto j=pgv->ioff+pgv->b->NinW; j<sizeof(dmG[i])/sizeof(double); j++){
			dmG[i][j] = 0.0;
		}
	}
}
pbou->updateGhostR2(dpG, pbou->S_DP);
pbou->updateGhostR2(dmG, pbou->S_DM);
};*/








