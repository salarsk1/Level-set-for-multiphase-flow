//Written by Salar Safarkhani


#include"advection.h"

void advection::lit_advection(const double dt){
	std::unique_ptr<global_variable> pgv;
	std::unique_ptr<parallel> pp;
	std::unique_ptr<weno> pweno;
	std::unique_ptr<bound> pbou;
	std::unique_ptr<timing> pti;
	std::unique_ptr<litBuffer> plbuf;
	std::unique_ptr<redist> plred;
	std::unique_ptr<reinit> prein;
	std::unique_ptr<litparam> ppar;
std::array<std::array<double, 3>, 3> alphaRK;
std::array<double, 3> grad_G;
//std::ofstream myfile;
double **dpG, **dmdpG, **LG;
double * Gold;
double G_WENO, base, a;
pti.reset(new timing);
pgv.reset(new global_variable);
pp.reset(new parallel);
pweno.reset(new weno);
pbou.reset(new bound);
pti.reset(new timing);
plbuf.reset(new litBuffer);
prein.reset(new reinit);
ppar.reset(new litparam);
pti->lit_timing_start("advection");
if(pgv->verbose && pp->myrank==0)
	std::cout<<pgv->clit<<"Starting lit_advection"<<std::endl; //correct it
if(ppar->trim(pgv->schemeAdvect) == "UC-5" || ppar->trim(pgv->schemeAdvect) =="UC-7" ||\
	ppar->trim(pgv->schemeAdvect) == "UC-9" || ppar->trim(pgv->schemeAdvect) =="UC-11"){
	alphaRK[0][0]=dt; 		alphaRK[1][0]=0.0; 		alphaRK[2][0]=0.0;
	alphaRK[0][1]=0.25*dt;  alphaRK[1][1]=0.25*dt;  alphaRK[2][1]=0.0;
	alphaRK[0][2]=dt/6.0; 	alphaRK[1][2]=dt/6.0; 	alphaRK[2][2]=2.0*dt/3.0;
}
else{
	alphaRK[0][0]=dt; 		 alphaRK[1][0]=0.0; 		  alphaRK[2][0]=0.0;
	alphaRK[0][1]=-0.75*dt;  alphaRK[1][1]=0.25*dt;   alphaRK[2][1]=0.0;
	alphaRK[0][2]=-r112*dt;  alphaRK[1][2]=-r112*dt;  alphaRK[2][2]=r23*dt;
}
size_t LG_S1, LG_S2;
LG_S1 = 0;
LG_S2 = 0;
LG = plbuf->getR2Buffer(LG_S1, LG_S2);
size_t Gold_S;
Gold = plbuf->getR1Buffer(Gold_S);
size_t dpG_S1, dpG_S2, dmdpG_S1, dmdpG_S2;
if(ppar->trim(pgv->schemeAdvect) == "WENO-3"||ppar->trim(pgv->schemeAdvect)=="WENO-5"){
	dpG   = plbuf->getR2Buffer(dpG_S1, dpG_S2);
	dmdpG = plbuf->getR2Buffer(dmdpG_S1, dmdpG_S2);
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinN; ic++){
		Gold[pgv->ioff+ic] = pgv->Gn[ic].G;
	}
}

for(int irk=0; irk<RK_advection; irk++){
	if(ppar->trim(pgv->schemeAdvect)=="WENO-3" || ppar->trim(pgv->schemeAdvect)=="WENO-5")
		pweno->prepare_WENO_ghost(dpG_S1, dpG_S2, dpG,dmdpG_S1, dmdpG_S2, dmdpG);
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		if(ppar->trim(pgv->schemeAdvect)=="UPWIND-1"){
			if(pgv->cylindrical){
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=(pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[0];
					}
					else{
						grad_G[0]=(pgv->Gn[ic].G-pgv->Gn[pgv->i2c[i-1][j][k]-1].G)*pgv->rdxyz[0];
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[1];
					}
					else{
						grad_G[1]=(pgv->Gn[ic].G-pgv->Gn[pgv->i2c[i][j-1][k]-1].G)*pgv->rdxyz[1];
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=(pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[ic].G)*\
									  pgv->rdxyz[2]*pgv->ryc[j+pgv->nghost];
					}
					else{
						grad_G[2]=(pgv->Gn[ic].G-pgv->Gn[pgv->i2c[i][j][k-1]-1].G)*\
									  pgv->rdxyz[2]*pgv->ryc[j+pgv->nghost];
					}
					LG[irk][pgv->ioff+ic] = -(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].\
													  V[1]+grad_G[2]*pgv->Gn[ic].V[2]);
				}
			}
			else{
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;



					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0] = (pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[0];
					}
					else{
						grad_G[0] = (pgv->Gn[ic].G-pgv->Gn[pgv->i2c[i-1][j][k]-1].G)*pgv->rdxyz[0];
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1] = (pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[1];
					}
					else{
						grad_G[1] = (pgv->Gn[ic].G-pgv->Gn[pgv->i2c[i][j-1][k]-1].G)*pgv->rdxyz[1];
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2] = (pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[ic].G)*pgv->rdxyz[2];
					}
					else{
						grad_G[2] = (pgv->Gn[ic].G-pgv->Gn[pgv->i2c[i][j][k-1]-1].G)*pgv->rdxyz[2];
					}
					LG[irk][pgv->ioff+ic] = -(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].V[1]\
													 +grad_G[2]*pgv->Gn[ic].V[2]);
				}
			}
			for(auto ic=0; ic<pgv->b->NinN; ic++){
				auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
				auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
				auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
				for(auto irki=0; irki<=irk; irki++){
					pgv->Gn[ic].G += alphaRK[irki][irk]*LG[irki][pgv->ioff+ic];
				}
			}
		}
		else if(ppar->trim(pgv->schemeAdvect) == "UC-5"){
			if(pgv->cylindrical){
				for(auto ic=0; ic<pgv->b->NinN;ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=(-pgv->ucWeight5[0][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].\
									G-pgv->ucWeight5[1][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G-\
									pgv->ucWeight5[2][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G-\
									pgv->ucWeight5[3][0]*pgv->Gn[ic].G-pgv->ucWeight5[4][0]*\
									pgv->Gn[pgv->i2c[i-1][j][k]-1].G-pgv->ucWeight5[5][0]*\
									pgv->Gn[pgv->i2c[i-2][j][k]-1].G);
					}
					else{
						grad_G[0]=(pgv->ucWeight5[0][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G+\
									  pgv->ucWeight5[1][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G+\
									  pgv->ucWeight5[2][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G+\
									  pgv->ucWeight5[3][0]*pgv->Gn[ic].G+pgv->ucWeight5[4][0]\
									  *pgv->Gn[pgv->i2c[i+1][j][k]-1].G+pgv->ucWeight5[5][0]*\
									  pgv->Gn[pgv->i2c[i+2][j][k]-1].G);
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=(-pgv->ucWeight5[0][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].\
									G-pgv->ucWeight5[1][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G-\
									pgv->ucWeight5[2][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G-\
									pgv->ucWeight5[3][1]*pgv->Gn[ic].G-pgv->ucWeight5[4][1]*\
									pgv->Gn[pgv->i2c[i][j-1][k]-1].G-pgv->ucWeight5[5][1]*\
									pgv->Gn[pgv->i2c[i][j-2][k]-1].G);
					}
					else{
						grad_G[1]=(pgv->ucWeight5[0][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G+\
									  pgv->ucWeight5[1][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G+\
									  pgv->ucWeight5[2][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G+\
									  pgv->ucWeight5[3][1]*pgv->Gn[ic].G+pgv->ucWeight5[4][1]\
									  *pgv->Gn[pgv->i2c[i][j+1][k]-1].G+pgv->ucWeight5[5][1]*\
									  pgv->Gn[pgv->i2c[i][j+2][k]-1].G);
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=(-pgv->ucWeight5[0][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].\
									G-pgv->ucWeight5[1][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G-\
									pgv->ucWeight5[2][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G-\
									pgv->ucWeight5[3][2]*pgv->Gn[ic].G-pgv->ucWeight5[4][2]*\
									pgv->Gn[pgv->i2c[i][j][k-1]-1].G-pgv->ucWeight5[5][2]*\
									pgv->Gn[pgv->i2c[i][j][k-2]-1].G)*pgv->ryc[j+pgv->b->jmino_-\
									(1-pgv->nghost)];
					}
					else{
						grad_G[2]=(pgv->ucWeight5[0][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G+\
									  pgv->ucWeight5[1][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G+\
									  pgv->ucWeight5[2][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G+\
									  pgv->ucWeight5[3][2]*pgv->Gn[ic].G+pgv->ucWeight5[4][2]\
									  *pgv->Gn[pgv->i2c[i][j][k+1]-1].G+pgv->ucWeight5[5][2]*\
									  pgv->Gn[pgv->i2c[i][j][k+2]-1].G)*pgv->ryc[j+pgv->b->jmino_-\
									  (1-pgv->nghost)];
					}
					LG[irk][pgv->ioff+ic]-(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].\
											  V[1]+grad_G[2]*pgv->Gn[ic].V[2]);
				}
			}
			else{
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;	
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=(-pgv->ucWeight5[0][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G-\
										pgv->ucWeight5[1][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G-\
										pgv->ucWeight5[2][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G-\
										pgv->ucWeight5[3][0]*pgv->Gn[ic].G-
										pgv->ucWeight5[4][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G-\
										pgv->ucWeight5[5][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G);
					}
					else{
						grad_G[0]=(pgv->ucWeight5[0][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G-\
									  pgv->ucWeight5[1][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G-\
									  pgv->ucWeight5[2][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G-\
									  pgv->ucWeight5[3][0]*pgv->Gn[ic].G-\
									  pgv->ucWeight5[4][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G-\
									  pgv->ucWeight5[5][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G);
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=(-pgv->ucWeight5[0][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G-\
										pgv->ucWeight5[1][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G-\
										pgv->ucWeight5[2][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G-\
										pgv->ucWeight5[3][1]*pgv->Gn[ic].G-
										pgv->ucWeight5[4][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G-\
										pgv->ucWeight5[5][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G);
					}
					else{
						grad_G[1]=(pgv->ucWeight5[0][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G-\
									  pgv->ucWeight5[1][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G-\
									  pgv->ucWeight5[2][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G-\
									  pgv->ucWeight5[3][1]*pgv->Gn[ic].G-\
									  pgv->ucWeight5[4][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G-\
									  pgv->ucWeight5[5][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G);
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=(-pgv->ucWeight5[0][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G-\
										pgv->ucWeight5[1][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G-\
										pgv->ucWeight5[2][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G-\
										pgv->ucWeight5[3][2]*pgv->Gn[ic].G-
										pgv->ucWeight5[4][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G-\
										pgv->ucWeight5[5][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G);
					}
					else{
						grad_G[2]=(pgv->ucWeight5[0][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G-\
									  pgv->ucWeight5[1][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G-\
									  pgv->ucWeight5[2][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G-\
									  pgv->ucWeight5[3][2]*pgv->Gn[ic].G-\
									  pgv->ucWeight5[4][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G-\
									  pgv->ucWeight5[5][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G);
					}
					LG[irk][pgv->ioff+ic]=-(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].\
											  V[1]+grad_G[2]*pgv->Gn[ic].V[2]);
				}
			}
			for(auto ic=0; ic<pgv->b->NinN; ic++){
				pgv->Gn[ic].G=Gold[pgv->ioff+ic];

			}
			for(auto ic=0; ic<pgv->b->NinN; ic++){
				for(auto irki=0; irki<=irk; irki++){
					pgv->Gn[ic].G += alphaRK[irki][irk]*LG[irki][pgv->ioff+ic];
				}
			}
		}
		else if(ppar->trim(pgv->schemeAdvect) == "UC-7"){
			if(pgv->cylindrical){
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=(-pgv->ucWeight7[0][0]*pgv->Gn[pgv->i2c[i+4][j][k]-1].G-\
										pgv->ucWeight7[1][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G-\
										pgv->ucWeight7[2][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G-\
										pgv->ucWeight7[3][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G-\
										pgv->ucWeight7[4][0]*pgv->Gn[ic].G-\
										pgv->ucWeight7[5][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G-\
										pgv->ucWeight7[6][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G-\
										pgv->ucWeight7[7][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G);
					}
					else{
						grad_G[0]=(pgv->ucWeight7[0][0]*pgv->Gn[pgv->i2c[i-4][j][k]-1].G+\
									  pgv->ucWeight7[1][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G+\
									  pgv->ucWeight7[2][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G+\
									  pgv->ucWeight7[3][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G+\
									  pgv->ucWeight7[4][0]*pgv->Gn[ic].G+\
									  pgv->ucWeight7[5][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G+\
									  pgv->ucWeight7[6][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G+\
									  pgv->ucWeight7[7][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G);
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=(-pgv->ucWeight7[0][1]*pgv->Gn[pgv->i2c[i][j+4][k]-1].G\
										-pgv->ucWeight7[1][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G\
										-pgv->ucWeight7[2][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G\
										-pgv->ucWeight7[3][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G\
										-pgv->ucWeight7[4][1]*pgv->Gn[ic].G-pgv->ucWeight7[5][1]*\
										 pgv->Gn[pgv->i2c[i][j-1][k]-1].G-pgv->ucWeight7[6][1]*\
										 pgv->Gn[pgv->i2c[i][j-2][k]-1].G-pgv->ucWeight7[7][1]*\
										 pgv->Gn[pgv->i2c[i][j-3][k]-1].G);
					}
					else{
						grad_G[1]=(pgv->ucWeight7[0][1]*pgv->Gn[pgv->i2c[i][j-4][k]-1].G+\
									  pgv->ucWeight7[1][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G+\
									  pgv->ucWeight7[2][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G+\
									  pgv->ucWeight7[3][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G+\
									  pgv->ucWeight7[4][1]*pgv->Gn[ic].G+\
									  pgv->ucWeight7[5][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G+\
									  pgv->ucWeight7[6][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G+\
									  pgv->ucWeight7[7][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G);
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=(-pgv->ucWeight7[0][2]*pgv->Gn[pgv->i2c[i][j][k+4]-1].G-\
										pgv->ucWeight7[1][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G-\
										pgv->ucWeight7[2][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G-\
										pgv->ucWeight7[3][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G-\
										pgv->ucWeight7[4][2]*pgv->Gn[ic].G-\
										pgv->ucWeight7[5][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G-\
										pgv->ucWeight7[6][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G-\
										pgv->ucWeight7[7][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G)*\
										pgv->ryc[(j+pgv->b->jmino_)-(1-pgv->nghost)];
					}
					else{
						grad_G[2]=(pgv->ucWeight7[0][2]*pgv->Gn[pgv->i2c[i][j][k-4]-1].G+\
									  pgv->ucWeight7[1][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G+\
									  pgv->ucWeight7[2][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G+\
									  pgv->ucWeight7[3][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G+\
									  pgv->ucWeight7[4][2]*pgv->Gn[ic].G+\
									  pgv->ucWeight7[5][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G+\
									  pgv->ucWeight7[6][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G+\
									  pgv->ucWeight7[7][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G)*\
									  pgv->ryc[(j+pgv->b->jmino_)-(1-pgv->nghost)];
					}
					LG[irk][pgv->ioff+ic] = -(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].\
											  V[1]+grad_G[2]*pgv->Gn[ic].V[2]);
				}
			}
			else{
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=(-pgv->ucWeight7[0][0]*pgv->Gn[pgv->i2c[i+4][j][k]-1].G-\
										pgv->ucWeight7[1][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G-\
										pgv->ucWeight7[2][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G-\
										pgv->ucWeight7[3][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G-\
										pgv->ucWeight7[4][0]*pgv->Gn[ic].G-\
										pgv->ucWeight7[5][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G-\
										pgv->ucWeight7[6][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G-\
										pgv->ucWeight7[7][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G);
					}
					else{
						grad_G[0]=(pgv->ucWeight7[0][0]*pgv->Gn[pgv->i2c[i-4][j][k]-1].G+\
									  pgv->ucWeight7[1][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G+\
									  pgv->ucWeight7[2][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G+\
									  pgv->ucWeight7[3][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G+\
									  pgv->ucWeight7[4][0]*pgv->Gn[ic].G+\
									  pgv->ucWeight7[5][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G+\
									  pgv->ucWeight7[6][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G+\
									  pgv->ucWeight7[7][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G);
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=(-pgv->ucWeight7[0][1]*pgv->Gn[pgv->i2c[i][j+4][k]-1].G-\
										pgv->ucWeight7[1][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G-\
										pgv->ucWeight7[2][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G-\
										pgv->ucWeight7[3][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G-\
										pgv->ucWeight7[4][1]*pgv->Gn[ic].G-\
										pgv->ucWeight7[5][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G-\
										pgv->ucWeight7[6][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G-\
										pgv->ucWeight7[7][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G);
					}
					else{
						grad_G[1]=(pgv->ucWeight7[0][1]*pgv->Gn[pgv->i2c[i][j-4][k]-1].G+\
									  pgv->ucWeight7[1][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G+\
									  pgv->ucWeight7[2][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G+\
									  pgv->ucWeight7[3][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G+\
									  pgv->ucWeight7[4][1]*pgv->Gn[ic].G+\
									  pgv->ucWeight7[5][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G+\
									  pgv->ucWeight7[6][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G+\
									  pgv->ucWeight7[7][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G);
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=(-pgv->ucWeight7[0][2]*pgv->Gn[pgv->i2c[i][j][k+4]-1].G-\
										pgv->ucWeight7[1][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G-\
										pgv->ucWeight7[2][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G-\
										pgv->ucWeight7[3][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G-\
										pgv->ucWeight7[4][2]*pgv->Gn[ic].G-\
										pgv->ucWeight7[5][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G-\
										pgv->ucWeight7[6][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G-\
										pgv->ucWeight7[7][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G)*\
										pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
					}
					else{
						grad_G[2]=(pgv->ucWeight7[0][2]*pgv->Gn[pgv->i2c[i][j][k-4]-1].G+\
									  pgv->ucWeight7[1][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G+\
									  pgv->ucWeight7[2][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G+\
									  pgv->ucWeight7[3][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G+\
									  pgv->ucWeight7[4][2]*pgv->Gn[ic].G+\
									  pgv->ucWeight7[5][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G+\
									  pgv->ucWeight7[6][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G+\
									  pgv->ucWeight7[7][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G)*\
									  pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
					}
					LG[irk][pgv->ioff+ic]=-(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].\
											  V[1]+grad_G[2]*pgv->Gn[ic].V[2]);
				}
			}
			for(auto ic=0; ic<pgv->b->NinN; ic++){
				pgv->Gn[ic].G = Gold[pgv->ioff+ic];
			}
			for(auto ic=0; ic<pgv->b->NinN; ic++){
				for(auto irki=0; irki<=irk; irki++){
					pgv->Gn[ic].G += alphaRK[irki][irk]*LG[irki][pgv->ioff+ic];
				}
			}
		}
		else if(ppar->trim(pgv->schemeAdvect) == "UC-9"){
			if(pgv->cylindrical){
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=(-pgv->ucWeight9[0][0]*pgv->Gn[pgv->i2c[i+5][j][k]-1].G-\
										pgv->ucWeight9[1][0]*pgv->Gn[pgv->i2c[i+4][j][k]-1].G-\
										pgv->ucWeight9[2][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G-\
										pgv->ucWeight9[3][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G-\
										pgv->ucWeight9[4][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G-\
										pgv->ucWeight9[5][0]*pgv->Gn[ic].G-\
										pgv->ucWeight9[6][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G-\
										pgv->ucWeight9[7][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G-\
										pgv->ucWeight9[8][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G-\
										pgv->ucWeight9[9][0]*pgv->Gn[pgv->i2c[i-4][j][k]-1].G);
					}
					else{
						grad_G[0]=(pgv->ucWeight9[0][0]*pgv->Gn[pgv->i2c[i-5][j][k]-1].G+\
									  pgv->ucWeight9[1][0]*pgv->Gn[pgv->i2c[i-4][j][k]-1].G+\
									  pgv->ucWeight9[2][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G+\
									  pgv->ucWeight9[3][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G+\
									  pgv->ucWeight9[4][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G+\
									  pgv->ucWeight9[5][0]*pgv->Gn[ic].G-\
									  pgv->ucWeight9[6][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G+\
									  pgv->ucWeight9[7][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G+\
									  pgv->ucWeight9[8][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G+\
									  pgv->ucWeight9[9][0]*pgv->Gn[pgv->i2c[i+4][j][k]-1].G);
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=(-pgv->ucWeight9[0][1]*pgv->Gn[pgv->i2c[i][j+5][k]-1].G-\
										pgv->ucWeight9[1][1]*pgv->Gn[pgv->i2c[i][j+4][k]-1].G-\
										pgv->ucWeight9[2][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G-\
										pgv->ucWeight9[3][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G-\
										pgv->ucWeight9[4][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G-\
										pgv->ucWeight9[5][1]*pgv->Gn[ic].G-\
										pgv->ucWeight9[6][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G-\
										pgv->ucWeight9[7][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G-\
										pgv->ucWeight9[8][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G-\
										pgv->ucWeight9[9][1]*pgv->Gn[pgv->i2c[i][j-4][k]-1].G);
					}
					else{
						grad_G[1]=(pgv->ucWeight9[0][1]*pgv->Gn[pgv->i2c[i][j-5][k]-1].G+\
									  pgv->ucWeight9[1][1]*pgv->Gn[pgv->i2c[i][j-4][k]-1].G+\
									  pgv->ucWeight9[2][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G+\
									  pgv->ucWeight9[3][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G+\
									  pgv->ucWeight9[4][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G+\
									  pgv->ucWeight9[5][1]*pgv->Gn[ic].G-\
									  pgv->ucWeight9[6][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G+\
									  pgv->ucWeight9[7][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G+\
									  pgv->ucWeight9[8][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G+\
									  pgv->ucWeight9[9][1]*pgv->Gn[pgv->i2c[i][j+4][k]-1].G);
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=(-pgv->ucWeight9[0][2]*pgv->Gn[pgv->i2c[i][j][k+5]-1].G-\
										pgv->ucWeight9[1][2]*pgv->Gn[pgv->i2c[i][j][k+4]-1].G-\
										pgv->ucWeight9[2][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G-\
										pgv->ucWeight9[3][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G-\
										pgv->ucWeight9[4][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G-\
										pgv->ucWeight9[5][2]*pgv->Gn[ic].G-\
										pgv->ucWeight9[6][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G-\
										pgv->ucWeight9[7][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G-\
										pgv->ucWeight9[8][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G-\
										pgv->ucWeight9[9][2]*pgv->Gn[pgv->i2c[i][j][k-4]-1].G)*\
										pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
					}
					else{
						grad_G[2]=(pgv->ucWeight9[0][2]*pgv->Gn[pgv->i2c[i][j][k-5]-1].G+\
									  pgv->ucWeight9[1][2]*pgv->Gn[pgv->i2c[i][j][k-4]-1].G+\
									  pgv->ucWeight9[2][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G+\
									  pgv->ucWeight9[3][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G+\
									  pgv->ucWeight9[4][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G+\
									  pgv->ucWeight9[5][2]*pgv->Gn[ic].G-\
									  pgv->ucWeight9[6][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G+\
									  pgv->ucWeight9[7][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G+\
									  pgv->ucWeight9[8][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G+\
									  pgv->ucWeight9[9][2]*pgv->Gn[pgv->i2c[i][j][k+4]-1].G)*\
									  pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
					}
					LG[irk][pgv->ioff+ic]=-(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].\
											  V[1]+grad_G[2]*pgv->Gn[ic].V[2]);
				}
			}
			else{
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=(-pgv->ucWeight9[0][0]*pgv->Gn[pgv->i2c[i+5][j][k]-1].G-\
										pgv->ucWeight9[1][0]*pgv->Gn[pgv->i2c[i+4][j][k]-1].G-\
										pgv->ucWeight9[2][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G-\
										pgv->ucWeight9[3][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G-\
										pgv->ucWeight9[4][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G-\
										pgv->ucWeight9[5][0]*pgv->Gn[ic].G-\
										pgv->ucWeight9[6][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G-\
										pgv->ucWeight9[7][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G-\
										pgv->ucWeight9[8][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G-\
										pgv->ucWeight9[9][0]*pgv->Gn[pgv->i2c[i-4][j][k]-1].G);
					}
					else{
						grad_G[0]=(pgv->ucWeight9[0][0]*pgv->Gn[pgv->i2c[i-5][j][k]-1].G+\
									  pgv->ucWeight9[1][0]*pgv->Gn[pgv->i2c[i-4][j][k]-1].G+\
									  pgv->ucWeight9[2][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G+\
									  pgv->ucWeight9[3][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G+\
									  pgv->ucWeight9[4][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G+\
									  pgv->ucWeight9[5][0]*pgv->Gn[ic].G-\
									  pgv->ucWeight9[6][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G+\
									  pgv->ucWeight9[7][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G+\
									  pgv->ucWeight9[8][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G+\
									  pgv->ucWeight9[9][0]*pgv->Gn[pgv->i2c[i+4][j][k]-1].G);
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=(-pgv->ucWeight9[0][1]*pgv->Gn[pgv->i2c[i][j+5][k]-1].G-\
										pgv->ucWeight9[1][1]*pgv->Gn[pgv->i2c[i][j+4][k]-1].G-\
										pgv->ucWeight9[2][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G-\
										pgv->ucWeight9[3][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G-\
										pgv->ucWeight9[4][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G-\
										pgv->ucWeight9[5][1]*pgv->Gn[ic].G-\
										pgv->ucWeight9[6][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G-\
										pgv->ucWeight9[7][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G-\
										pgv->ucWeight9[8][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G-\
										pgv->ucWeight9[9][1]*pgv->Gn[pgv->i2c[i][j-4][k]-1].G);
					}
					else{
						grad_G[1]=(pgv->ucWeight9[0][1]*pgv->Gn[pgv->i2c[i][j-5][k]-1].G+\
									  pgv->ucWeight9[1][1]*pgv->Gn[pgv->i2c[i][j-4][k]-1].G+\
									  pgv->ucWeight9[2][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G+\
									  pgv->ucWeight9[3][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G+\
									  pgv->ucWeight9[4][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G+\
									  pgv->ucWeight9[5][1]*pgv->Gn[ic].G-\
									  pgv->ucWeight9[6][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G+\
									  pgv->ucWeight9[7][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G+\
									  pgv->ucWeight9[8][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G+\
									  pgv->ucWeight9[9][1]*pgv->Gn[pgv->i2c[i][j+4][k]-1].G);
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=(-pgv->ucWeight9[0][2]*pgv->Gn[pgv->i2c[i][j][k+5]-1].G-\
										pgv->ucWeight9[1][2]*pgv->Gn[pgv->i2c[i][j][k+4]-1].G-\
										pgv->ucWeight9[2][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G-\
										pgv->ucWeight9[3][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G-\
										pgv->ucWeight9[4][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G-\
										pgv->ucWeight9[5][2]*pgv->Gn[ic].G-\
										pgv->ucWeight9[6][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G-\
										pgv->ucWeight9[7][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G-\
										pgv->ucWeight9[8][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G-\
										pgv->ucWeight9[9][2]*pgv->Gn[pgv->i2c[i][j][k-4]-1].G);
					}
					else{
						grad_G[2]=(pgv->ucWeight9[0][2]*pgv->Gn[pgv->i2c[i][j][k-5]-1].G+\
									  pgv->ucWeight9[1][2]*pgv->Gn[pgv->i2c[i][j][k-4]-1].G+\
									  pgv->ucWeight9[2][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G+\
									  pgv->ucWeight9[3][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G+\
									  pgv->ucWeight9[4][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G+\
									  pgv->ucWeight9[5][2]*pgv->Gn[ic].G-\
									  pgv->ucWeight9[6][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G+\
									  pgv->ucWeight9[7][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G+\
									  pgv->ucWeight9[8][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G+\
									  pgv->ucWeight9[9][2]*pgv->Gn[pgv->i2c[i][j][k+4]-1].G);
					}
					LG[irk][pgv->ioff+ic]=-(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].\
											  V[1]+grad_G[2]*pgv->Gn[ic].V[2]);

				}
			}
			for(auto ic=0; ic<pgv->b->NinN; ic++){
				pgv->Gn[ic].G = Gold[pgv->ioff+ic];
			}
			for(auto ic=0; ic<pgv->b->NinN; ic++){
				for(auto irki=0; irki<=irk; irki++){
					pgv->Gn[ic].G += alphaRK[irki][irk]*LG[irki][pgv->ioff+ic];
				}
			}
		}
		else if(ppar->trim(pgv->schemeAdvect) == "UC-11"){
			if(pgv->cylindrical){
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=(-pgv->ucWeight11[0][0]*pgv->Gn[pgv->i2c[i+6][j][k]-1].G-\
										pgv->ucWeight11[1][0]*pgv->Gn[pgv->i2c[i+5][j][k]-1].G-\
										pgv->ucWeight11[2][0]*pgv->Gn[pgv->i2c[i+4][j][k]-1].G-\
										pgv->ucWeight11[3][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G-\
										pgv->ucWeight11[4][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G-\
										pgv->ucWeight11[5][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G-\
										pgv->ucWeight11[6][0]*pgv->Gn[ic].G-\
										pgv->ucWeight11[7][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G-\
										pgv->ucWeight11[8][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G-\
										pgv->ucWeight11[9][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G-\
										pgv->ucWeight11[10][0]*pgv->Gn[pgv->i2c[i-4][j][k]-1].G-\
									   pgv->ucWeight11[11][0]*pgv->Gn[pgv->i2c[i-5][j][k]-1].G);
					}
					else{
						grad_G[0]=(pgv->ucWeight11[0][0]*pgv->Gn[pgv->i2c[i-6][j][k]-1].G+\
									  pgv->ucWeight11[1][0]*pgv->Gn[pgv->i2c[i-5][j][k]-1].G+\
									  pgv->ucWeight11[2][0]*pgv->Gn[pgv->i2c[i-4][j][k]-1].G+\
									  pgv->ucWeight11[3][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G+\
									  pgv->ucWeight11[4][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G+\
									  pgv->ucWeight11[5][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G+\
									  pgv->ucWeight11[6][0]*pgv->Gn[ic].G-\
								 	  pgv->ucWeight11[7][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G+\
								  	  pgv->ucWeight11[8][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G+\
								  	  pgv->ucWeight11[9][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G+\
								     pgv->ucWeight11[10][0]*pgv->Gn[pgv->i2c[i+4][j][k]-1].G+\
								     pgv->ucWeight11[11][0]*pgv->Gn[pgv->i2c[i+5][j][k]-1].G);
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=(-pgv->ucWeight11[0][1]*pgv->Gn[pgv->i2c[i][j+6][k]-1].G-\
										pgv->ucWeight11[1][1]*pgv->Gn[pgv->i2c[i][j+5][k]-1].G-\
										pgv->ucWeight11[2][1]*pgv->Gn[pgv->i2c[i][j+4][k]-1].G-\
										pgv->ucWeight11[3][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G-\
										pgv->ucWeight11[4][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G-\
										pgv->ucWeight11[5][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G-\
										pgv->ucWeight11[6][1]*pgv->Gn[ic].G-\
										pgv->ucWeight11[7][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G-\
										pgv->ucWeight11[8][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G-\
										pgv->ucWeight11[9][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G-\
						  			  pgv->ucWeight11[10][1]*pgv->Gn[pgv->i2c[i][j-4][k]-1].G-\
									  pgv->ucWeight11[11][1]*pgv->Gn[pgv->i2c[i][j-5][k]-1].G);
					}
					else{
						grad_G[1]=(pgv->ucWeight11[0][1]*pgv->Gn[pgv->i2c[i][j-6][k]-1].G+\
									  pgv->ucWeight11[1][1]*pgv->Gn[pgv->i2c[i][j-5][k]-1].G+\
									  pgv->ucWeight11[2][1]*pgv->Gn[pgv->i2c[i][j-4][k]-1].G+\
									  pgv->ucWeight11[3][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G+\
									  pgv->ucWeight11[4][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G+\
									  pgv->ucWeight11[5][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G+\
									  pgv->ucWeight11[6][1]*pgv->Gn[ic].G-\
									  pgv->ucWeight11[7][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G+\
									  pgv->ucWeight11[8][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G+\
									  pgv->ucWeight11[9][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G+\
									  pgv->ucWeight11[10][1]*pgv->Gn[pgv->i2c[i][j+4][k]-1].G+\
									  pgv->ucWeight11[11][1]*pgv->Gn[pgv->i2c[i][j+5][k]-1].G);
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=(-pgv->ucWeight11[0][2]*pgv->Gn[pgv->i2c[i][j][k+6]-1].G-\
										pgv->ucWeight11[1][2]*pgv->Gn[pgv->i2c[i][j][k+5]-1].G-\
										pgv->ucWeight11[2][2]*pgv->Gn[pgv->i2c[i][j][k+4]-1].G-\
										pgv->ucWeight11[3][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G-\
										pgv->ucWeight11[4][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G-\
										pgv->ucWeight11[5][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G-\
										pgv->ucWeight11[6][2]*pgv->Gn[ic].G-\
										pgv->ucWeight11[7][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G-\
										pgv->ucWeight11[8][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G-\
										pgv->ucWeight11[9][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G-\
									  	pgv->ucWeight11[10][2]*pgv->Gn[pgv->i2c[i][j][k-4]-1].G-\
									  	pgv->ucWeight11[11][2]*pgv->Gn[pgv->i2c[i][j][k-5]-1].G)*\
										pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
					}
					else{
						grad_G[2]=(pgv->ucWeight11[0][2]*pgv->Gn[pgv->i2c[i][j][k-6]-1].G+\
									  pgv->ucWeight11[1][2]*pgv->Gn[pgv->i2c[i][j][k-5]-1].G+\
									  pgv->ucWeight11[2][2]*pgv->Gn[pgv->i2c[i][j][k-4]-1].G+\
									  pgv->ucWeight11[3][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G+\
									  pgv->ucWeight11[4][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G+\
									  pgv->ucWeight11[5][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G+\
								  	  pgv->ucWeight11[6][2]*pgv->Gn[ic].G-\
								  	  pgv->ucWeight11[7][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G+\
								  	  pgv->ucWeight11[8][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G+\
								     pgv->ucWeight11[9][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G+\
								     pgv->ucWeight11[10][2]*pgv->Gn[pgv->i2c[i][j][k+4]-1].G+\
								     pgv->ucWeight11[11][2]*pgv->Gn[pgv->i2c[i][j][k+5]-1].G)*\
									  pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
					}
					LG[irk][pgv->ioff+ic] = -(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].\
												  V[1]+grad_G[2]*pgv->Gn[ic].V[2]);
				}
			}
			else{
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=(-pgv->ucWeight11[0][0]*pgv->Gn[pgv->i2c[i+6][j][k]-1].G-\
										pgv->ucWeight11[1][0]*pgv->Gn[pgv->i2c[i+5][j][k]-1].G-\
										pgv->ucWeight11[2][0]*pgv->Gn[pgv->i2c[i+4][j][k]-1].G-\
										pgv->ucWeight11[3][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G-\
										pgv->ucWeight11[4][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G-\
										pgv->ucWeight11[5][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G-\
										pgv->ucWeight11[6][0]*pgv->Gn[ic].G-\
										pgv->ucWeight11[7][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G-\
										pgv->ucWeight11[8][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G-\
										pgv->ucWeight11[9][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G-\
					 				   pgv->ucWeight11[10][0]*pgv->Gn[pgv->i2c[i-4][j][k]-1].G-\
									   pgv->ucWeight11[11][0]*pgv->Gn[pgv->i2c[i-5][j][k]-1].G);
					}
					else{
						grad_G[0]=(pgv->ucWeight11[0][0]*pgv->Gn[pgv->i2c[i-6][j][k]-1].G+\
									  pgv->ucWeight11[1][0]*pgv->Gn[pgv->i2c[i-5][j][k]-1].G+\
									  pgv->ucWeight11[2][0]*pgv->Gn[pgv->i2c[i-4][j][k]-1].G+\
									  pgv->ucWeight11[3][0]*pgv->Gn[pgv->i2c[i-3][j][k]-1].G+\
									  pgv->ucWeight11[4][0]*pgv->Gn[pgv->i2c[i-2][j][k]-1].G+\
									  pgv->ucWeight11[5][0]*pgv->Gn[pgv->i2c[i-1][j][k]-1].G+\
									  pgv->ucWeight11[6][0]*pgv->Gn[ic].G-\
									  pgv->ucWeight11[7][0]*pgv->Gn[pgv->i2c[i+1][j][k]-1].G+\
									  pgv->ucWeight11[8][0]*pgv->Gn[pgv->i2c[i+2][j][k]-1].G+\
									  pgv->ucWeight11[9][0]*pgv->Gn[pgv->i2c[i+3][j][k]-1].G+\
								     pgv->ucWeight11[10][0]*pgv->Gn[pgv->i2c[i+4][j][k]-1].G+\
								     pgv->ucWeight11[11][0]*pgv->Gn[pgv->i2c[i+5][j][k]-1].G);
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=(-pgv->ucWeight11[0][1]*pgv->Gn[pgv->i2c[i][j+6][k]-1].G-\
										pgv->ucWeight11[1][1]*pgv->Gn[pgv->i2c[i][j+5][k]-1].G-\
										pgv->ucWeight11[2][1]*pgv->Gn[pgv->i2c[i][j+4][k]-1].G-\
										pgv->ucWeight11[3][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G-\
										pgv->ucWeight11[4][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G-\
										pgv->ucWeight11[5][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G-\
										pgv->ucWeight11[6][1]*pgv->Gn[ic].G-\
										pgv->ucWeight11[7][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G-\
										pgv->ucWeight11[8][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G-\
										pgv->ucWeight11[9][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G-\
									   pgv->ucWeight11[10][1]*pgv->Gn[pgv->i2c[i][j-4][k]-1].G-\
									   pgv->ucWeight11[11][1]*pgv->Gn[pgv->i2c[i][j-5][k]-1].G);
					}
					else{
						grad_G[1]=(pgv->ucWeight11[0][1]*pgv->Gn[pgv->i2c[i][j-6][k]-1].G+\
									  pgv->ucWeight11[1][1]*pgv->Gn[pgv->i2c[i][j-5][k]-1].G+\
									  pgv->ucWeight11[2][1]*pgv->Gn[pgv->i2c[i][j-4][k]-1].G+\
									  pgv->ucWeight11[3][1]*pgv->Gn[pgv->i2c[i][j-3][k]-1].G+\
									  pgv->ucWeight11[4][1]*pgv->Gn[pgv->i2c[i][j-2][k]-1].G+\
									  pgv->ucWeight11[5][1]*pgv->Gn[pgv->i2c[i][j-1][k]-1].G+\
									  pgv->ucWeight11[6][1]*pgv->Gn[ic].G-\
									  pgv->ucWeight11[7][1]*pgv->Gn[pgv->i2c[i][j+1][k]-1].G+\
									  pgv->ucWeight11[8][1]*pgv->Gn[pgv->i2c[i][j+2][k]-1].G+\
									  pgv->ucWeight11[9][1]*pgv->Gn[pgv->i2c[i][j+3][k]-1].G+\
								     pgv->ucWeight11[10][1]*pgv->Gn[pgv->i2c[i][j+4][k]-1].G+\
								     pgv->ucWeight11[11][1]*pgv->Gn[pgv->i2c[i][j+5][k]-1].G);
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=(-pgv->ucWeight11[0][2]*pgv->Gn[pgv->i2c[i][j][k+6]-1].G-\
										pgv->ucWeight11[1][2]*pgv->Gn[pgv->i2c[i][j][k+5]-1].G-\
										pgv->ucWeight11[2][2]*pgv->Gn[pgv->i2c[i][j][k+4]-1].G-\
										pgv->ucWeight11[3][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G-\
										pgv->ucWeight11[4][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G-\
										pgv->ucWeight11[5][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G-\
										pgv->ucWeight11[6][2]*pgv->Gn[ic].G-\
										pgv->ucWeight11[7][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G-\
										pgv->ucWeight11[8][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G-\
										pgv->ucWeight11[9][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G-\
									   pgv->ucWeight11[10][2]*pgv->Gn[pgv->i2c[i][j][k-4]-1].G-\
									   pgv->ucWeight11[11][2]*pgv->Gn[pgv->i2c[i][j][k-5]-1].G);
					}
					else{
						grad_G[2]=(pgv->ucWeight11[0][2]*pgv->Gn[pgv->i2c[i][j][k-6]-1].G+\
									  pgv->ucWeight11[1][2]*pgv->Gn[pgv->i2c[i][j][k-5]-1].G+\
									  pgv->ucWeight11[2][2]*pgv->Gn[pgv->i2c[i][j][k-4]-1].G+\
									  pgv->ucWeight11[3][2]*pgv->Gn[pgv->i2c[i][j][k-3]-1].G+\
									  pgv->ucWeight11[4][2]*pgv->Gn[pgv->i2c[i][j][k-2]-1].G+\
									  pgv->ucWeight11[5][2]*pgv->Gn[pgv->i2c[i][j][k-1]-1].G+\
									  pgv->ucWeight11[6][2]*pgv->Gn[ic].G-\
									  pgv->ucWeight11[7][2]*pgv->Gn[pgv->i2c[i][j][k+1]-1].G+\
									  pgv->ucWeight11[8][2]*pgv->Gn[pgv->i2c[i][j][k+2]-1].G+\
									  pgv->ucWeight11[9][2]*pgv->Gn[pgv->i2c[i][j][k+3]-1].G+\
								     pgv->ucWeight11[10][2]*pgv->Gn[pgv->i2c[i][j][k+4]-1].G+\
								     pgv->ucWeight11[11][2]*pgv->Gn[pgv->i2c[i][j][k+5]-1].G);
					}
					LG[irk][pgv->ioff+ic] = -(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].\
													  V[1]+grad_G[2]*pgv->Gn[ic].V[2]);
				}
			}
			for(auto ic=0; ic<pgv->b->NinN; ic++){
				pgv->Gn[ic].G = Gold[pgv->ioff+ic];
			}
			for(auto ic=0; ic<pgv->b->NinN; ic++){
				for(auto irki=0; irki<=irk; irki++){
					pgv->Gn[ic].G += alphaRK[irki][irk]*LG[irki][pgv->ioff+ic];
				}
			}
		}
		else if(ppar->trim(pgv->schemeAdvect)=="WENO-3"){
			if(pgv->cylindrical){
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=0.5*(dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]+dpG[0][pgv->\
									 ioff+ic])-pweno->wWENO3(dmdpG[0][pgv->ioff+ic], dmdpG\
									 [0][pgv->ioff+pgv->i2c[i+1][j][k]-1])*(dmdpG[0][pgv->\
									 ioff+ic]-dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]);
					}
					else{
						grad_G[0]=0.5*(dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]+dpG[0]\
									 [pgv->ioff+ic])-pweno->wWENO3(dmdpG[0][pgv->ioff+ic], \
									 dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1])*(dmdpG[0]\
								   [pgv->ioff+ic]-dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]);
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=0.5*(dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]+dpG[1][pgv->\
									 ioff+ic])-pweno->wWENO3(dmdpG[1][pgv->ioff+ic], dmdpG\
									 [1][pgv->ioff+pgv->i2c[i][j+1][k]-1])*(dmdpG[1][pgv->\
									 ioff+ic]-dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]);
					}
					else{
						grad_G[1]=0.5*(dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]+dpG[1]\
									 [pgv->ioff+ic])-pweno->wWENO3(dmdpG[1][pgv->ioff+ic], \
									 dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1])*(dmdpG[1]\
								   [pgv->ioff+ic]-dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]);
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=(0.5*(dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]+dpG[2][pgv->\
									 ioff+ic])-pweno->wWENO3(dmdpG[2][pgv->ioff+ic], dmdpG\
									 [2][pgv->ioff+pgv->i2c[i][j][k+1]-1])*(dmdpG[0][pgv->\
									 ioff+ic]-dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]))*\
									 pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
					}
					else{
						grad_G[2]=(0.5*(dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]+dpG[2]\
									 [pgv->ioff+ic])-pweno->wWENO3(dmdpG[2][pgv->ioff+ic], \
									 dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1])*(dmdpG[2]\
								   [pgv->ioff+ic]-dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]))*\
									pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
					}
					LG[irk][pgv->ioff+ic]=-(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].\
											  V[1]+grad_G[2]*pgv->Gn[ic].V[2]);
					for(auto irki=0; irki<=irk; irki++){
						pgv->Gn[ic].G += alphaRK[irki][irk]*LG[irki][pgv->ioff+ic];
					}
				}
			}
			else{
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=0.5*(dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]+dpG[0][pgv->\
									 ioff+ic])-pweno->wWENO3(dmdpG[0][pgv->ioff+ic], dmdpG\
									 [0][pgv->ioff+pgv->i2c[i+1][j][k]-1])*(dmdpG[0][pgv->\
									 ioff+ic]-dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]);
					}
					else{
						grad_G[0]=0.5*(dpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]+dpG[0]\
									 [pgv->ioff+ic])-pweno->wWENO3(dmdpG[0][pgv->ioff+ic], \
									 dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1])*(dmdpG[0]\
								   [pgv->ioff+ic]-dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]);
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=0.5*(dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]+dpG[1][pgv->\
									 ioff+ic])-pweno->wWENO3(dmdpG[1][pgv->ioff+ic], dmdpG\
									 [1][pgv->ioff+pgv->i2c[i][j+1][k]-1])*(dmdpG[1][pgv->\
									 ioff+ic]-dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]);
					}
					else{
						grad_G[1]=0.5*(dpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]+dpG[1]\
									 [pgv->ioff+ic])-pweno->wWENO3(dmdpG[1][pgv->ioff+ic], \
									 dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1])*(dmdpG[1]\
								   [pgv->ioff+ic]-dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]);
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=0.5*(dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]+dpG[2][pgv->\
									 ioff+ic])-pweno->wWENO3(dmdpG[2][pgv->ioff+ic], dmdpG\
									 [2][pgv->ioff+pgv->i2c[i][j][k+1]-1])*(dmdpG[0][pgv->\
									 ioff+ic]-dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]);
					}
					else{
						grad_G[2]=0.5*(dpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]+dpG[2]\
									 [pgv->ioff+ic])-pweno->wWENO3(dmdpG[2][pgv->ioff+ic], \
									 dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1])*(dmdpG[2]\
								   [pgv->ioff+ic]-dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]);
					}
					LG[irk][pgv->ioff+ic]=-(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].\
											  V[1]+grad_G[2]*pgv->Gn[ic].V[2]);

					for(auto irki=0; irki<=irk; irki++){
						pgv->Gn[ic].G += alphaRK[irki][irk]*LG[irki][pgv->ioff+ic];
					}
				}
			}
		}
		else if(ppar->trim(pgv->schemeAdvect)=="WENO-5"){
			if(pgv->cylindrical){
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=r112*(-dpG[0][pgv->ioff+pgv->i2c[i-2][j][k]-1]+7.0*dpG\
									 [0][pgv->ioff+pgv->i2c[i-1][j][k]-1]+7.0*dpG[0][pgv->\
									 ioff+ic]-dpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1])+\
									 pweno->G_WENO_5th(dmdpG[0][pgv->ioff+pgv->i2c[i+2][j][k]\
									 -1],dmdpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1],dmdpG[0]\
									[pgv->ioff+ic],dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]);
					}
					else{
						grad_G[0]=r112*(-dpG[0][pgv->ioff+pgv->i2c[i-2][j][k]-1]+7.0*dpG\
									 [0][pgv->ioff+pgv->i2c[i-1][j][k]-1]+7.0*dpG[0][pgv->\
									 ioff+ic]-dpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1])-\
									 pweno->G_WENO_5th(dmdpG[0][pgv->ioff+pgv->i2c[i-2][j][k]\
									 -1],dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1],dmdpG[0]\
									[pgv->ioff+ic],dmdpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1]);
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=r112*(-dpG[1][pgv->ioff+pgv->i2c[i][j-2][k]-1]+7.0*dpG\
									 [1][pgv->ioff+pgv->i2c[i][j-1][k]-1]+7.0*dpG[1][pgv->\
									 ioff+ic]-dpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1])+\
									 pweno->G_WENO_5th(dmdpG[1][pgv->ioff+pgv->i2c[i][j+2][k]\
									 -1],dmdpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1],dmdpG[1]\
									[pgv->ioff+ic],dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]);
					}
					else{
						grad_G[1]=r112*(-dpG[1][pgv->ioff+pgv->i2c[i][j-2][k]-1]+7.0*dpG\
									 [1][pgv->ioff+pgv->i2c[i][j-1][k]-1]+7.0*dpG[1][pgv->\
									 ioff+ic]-dpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1])-\
									 pweno->G_WENO_5th(dmdpG[1][pgv->ioff+pgv->i2c[i][j-2][k]\
									 -1],dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1],dmdpG[1]\
									[pgv->ioff+ic],dmdpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1]);
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=r112*(-dpG[2][pgv->ioff+pgv->i2c[i][j][k-2]-1]+7.0*dpG\
									 [2][pgv->ioff+pgv->i2c[i][j][k-1]-1]+7.0*dpG[2][pgv->\
									 ioff+ic]-dpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1])+\
									 pweno->G_WENO_5th(dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+2]\
									 -1],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1],dmdpG[2]\
								  [pgv->ioff+ic],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1])*\
									pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
					}
					else{
						grad_G[2]=r112*(-dpG[2][pgv->ioff+pgv->i2c[i][j][k-2]-1]+7.0*dpG\
									 [2][pgv->ioff+pgv->i2c[i][j][k-1]-1]+7.0*dpG[2][pgv->\
									 ioff+ic]-dpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1])-\
									 pweno->G_WENO_5th(dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-2]\
									 -1],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1],dmdpG[2]\
								 [pgv->ioff+ic],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1])*\
								 pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
					}
					LG[irk][pgv->ioff+ic]=-(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].\
											  V[1]+grad_G[2]*pgv->Gn[ic].V[2]);
					for(auto irki=0; irki<irk; irki++){
						pgv->Gn[ic].G += alphaRK[irki][irk]*LG[irki+1][pgv->ioff+ic];
					}
				}
			}
			else{
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
					auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
					auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
					if(pgv->Gn[ic].V[0] <= 0.0){
						grad_G[0]=r112*(-dpG[0][pgv->ioff+pgv->i2c[i-2][j][k]-1]+7.0*dpG\
									 [0][pgv->ioff+pgv->i2c[i-1][j][k]-1]+7.0*dpG[0][pgv->\
									 ioff+ic]-dpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1])+\
									 pweno->G_WENO_5th(dmdpG[0][pgv->ioff+pgv->i2c[i+2][j][k]-1],\
									 dmdpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1],dmdpG[0]\
									 [pgv->ioff+ic],dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1]);
					}
					else{
						grad_G[0]=r112*(-dpG[0][pgv->ioff+pgv->i2c[i-2][j][k]-1]+7.0*dpG\
									 [0][pgv->ioff+pgv->i2c[i-1][j][k]-1]+7.0*dpG[0][pgv->\
									 ioff+ic]-dpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1])-\
									 pweno->G_WENO_5th(dmdpG[0][pgv->ioff+pgv->i2c[i-2][j][k]-1],\
									 dmdpG[0][pgv->ioff+pgv->i2c[i-1][j][k]-1],dmdpG[0]\
									 [pgv->ioff+ic],dmdpG[0][pgv->ioff+pgv->i2c[i+1][j][k]-1]);
					}
					if(pgv->Gn[ic].V[1] <= 0.0){
						grad_G[1]=r112*(-dpG[1][pgv->ioff+pgv->i2c[i][j-2][k]-1]+7.0*dpG\
									 [1][pgv->ioff+pgv->i2c[i][j-1][k]-1]+7.0*dpG[1][pgv->\
									 ioff+ic]-dpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1])+\
									 pweno->G_WENO_5th(dmdpG[1][pgv->ioff+pgv->i2c[i][j+2][k]-1],\
									 dmdpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1],dmdpG[1]\
									 [pgv->ioff+ic],dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1]);
					}
					else{
						grad_G[1]=r112*(-dpG[1][pgv->ioff+pgv->i2c[i][j-2][k]-1]+7.0*dpG\
									 [1][pgv->ioff+pgv->i2c[i][j-1][k]-1]+7.0*dpG[1][pgv->\
									 ioff+ic]-dpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1])-\
									 pweno->G_WENO_5th(dmdpG[1][pgv->ioff+pgv->i2c[i][j-2][k]-1],\
									 dmdpG[1][pgv->ioff+pgv->i2c[i][j-1][k]-1],dmdpG[1]\
									 [pgv->ioff+ic],dmdpG[1][pgv->ioff+pgv->i2c[i][j+1][k]-1]);
					}
					if(pgv->Gn[ic].V[2] <= 0.0){
						grad_G[2]=r112*(-dpG[2][pgv->ioff+pgv->i2c[i][j][k-2]-1]+7.0*dpG\
									 [2][pgv->ioff+pgv->i2c[i][j][k-1]-1]+7.0*dpG[2][pgv->\
									 ioff+ic]-dpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1])+\
									 pweno->G_WENO_5th(dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+2]-1],\
									 dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1],dmdpG[2]\
								    [pgv->ioff+ic],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1]);
					}
					else{
						grad_G[2]=r112*(-dpG[2][pgv->ioff+pgv->i2c[i][j][k-2]-1]+7.0*dpG\
									 [2][pgv->ioff+pgv->i2c[i][j][k-1]-1]+7.0*dpG[2][pgv->\
									 ioff+ic]-dpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1])-\
									 pweno->G_WENO_5th(dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-2]-1],\
									 dmdpG[2][pgv->ioff+pgv->i2c[i][j][k-1]-1],dmdpG[2]\
									 [pgv->ioff+ic],dmdpG[2][pgv->ioff+pgv->i2c[i][j][k+1]-1]);
					}
 
					LG[irk][pgv->ioff+ic]=-(grad_G[0]*pgv->Gn[ic].V[0]+grad_G[1]*pgv->Gn[ic].V[1]+\
													grad_G[2]*pgv->Gn[ic].V[2]);

					for(auto irki=0; irki<=irk; irki++){
						pgv->Gn[ic].G += alphaRK[irki][irk] * LG[irki][pgv->ioff+ic];
					}
				}
			}
		}
	}
	pbou->updateGhostNodes();
}
if(ppar->trim(pgv->schemeAdvect)=="WENO-3"||ppar->trim(pgv->schemeAdvect)=="WENO-5"){
	plbuf->freeR2Buffer(dpG_S1, dpG_S2, dpG);
	plbuf->freeR2Buffer(dmdpG_S1, dmdpG_S2, dmdpG);
}
plbuf->freeR2Buffer(LG_S1, LG_S2, LG);
plbuf->freeR1Buffer(Gold);
pti->lit_timing_stop("advection");
if(pgv->verbose && pp->myrank==0)
	std::cout<<pgv->clit<<"Starting lit_advection ... Done"<<std::endl;
};

void advection::filterVelocity(){
	std::unique_ptr<global_variable> pgv;
	std::unique_ptr<parallel> pp;
	std::unique_ptr<weno> pweno;
	std::unique_ptr<bound> pbou;
	std::unique_ptr<timing> pti;
	std::unique_ptr<litBuffer> plbuf;
	std::unique_ptr<redist> plred;
	std::unique_ptr<reinit> prein;
	std::unique_ptr<litparam> ppar;
	double **V;
	double**** Vnode;
	bool*** marker;
	double weight;
	std::array<double, 3> Vh;
	int iFilter;
	plred.reset(new redist);
	pgv.reset(new global_variable);
	pp.reset(new parallel);
	pbou.reset(new bound);
	plbuf.reset(new litBuffer);
	plred.reset(new redist);
	prein.reset(new reinit);
	if(velocityFiltered) return;
	velocityFiltered = true;
	if(pgv->nVelFilter <= 0) return;
	if(pgv->verbose && pp->myrank==0){
		std::cout<<pgv->clit<<"Starting filterVelocity, nFilter=";
		std::cout<<pgv->nVelFilter<<std::endl;
	}
size_t V_S1, V_S2;
V = plbuf->getR2Buffer(V_S1, V_S2);
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinZ; ic++){
		for(auto i=0; i<3; i++){
			V[i][pgv->ioff+ic] = pgv->Gn[ic].V[i];
		}
	}
}
plred->lit_redist_pde(V_S1, V_S2, V);
pbou->updateGhostR2(V, V_S1, V_S2, pbou->S_GRADIENT);
for(auto iFilter=0; iFilter<pgv->nVelFilter; iFilter++){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		if(pp->ierr==0){
			Vnode = new double***[3];
			for(auto i=0; i<3; i++){
				auto m1 = pgv->b->imax_ - pgv->b->imino_+1;
				Vnode[i] = new double**[m1];
				for(auto j=0; j<m1; j++){
					auto m2 = pgv->b->jmax_ - pgv->b->jmino_+1;
					Vnode[i][j] = new double*[m2];
					for(auto k=0; k<m2; k++){
						auto m3 = pgv->b->kmax_ - pgv->b->kmino_ +1;
						Vnode[i][j][k] = new double[m3];
					}
				}
			}
			auto m1 = pgv->b->imax_ - pgv->b->imino_ + 1;
			marker = new bool**[m1];
			for(auto i=0; i<m1; i++){
				auto m2 = pgv->b->jmax_ - pgv->b->jmino_ + 1;
				marker[i] = new bool*[m2];
				for(auto j=0; j<m2; j++){
					auto m3 = pgv->b->kmax_ - pgv->b->kmino_ + 1;
					marker[i][j] = new bool[m3];
				}
			}
			// it is already zero; line 877
			// it is already false; line 888
		}
		else pp->litError("filterVelocity","allocation error Vnode, marker");

		for(auto i=0; i<pgv->b->imax_ - pgv->b->imino_ +1; i++){
			for(auto j=0; j<pgv->b->jmax_ - pgv->b->jmino_ +1; j++){
				for(auto k=0; k<pgv->b->kmax_ - pgv->b->kmino_ +1; k++){
					weight = 0.0;
					if(pgv->i2c[i][j][k] >0){
						for(auto q=0; q<3; q++){
							Vnode[q][i][j][k] += V[q][pgv->ioff+pgv->i2c[i][j][k]-1];
							weight += 1.0;
						}
					}
					if(pgv->i2c[i+1][j][k] >0){
						for(auto q=0; q<3; q++){
							Vnode[q][i][j][k] += V[q][pgv->ioff+pgv->i2c[i+1][j][k]-1];
							weight += 1.0;
						}
					}
					if(pgv->i2c[i+1][j+1][k] > 0){
						for(auto q=0; q<3; q++){
							Vnode[q][i][j][k] += V[q][pgv->ioff+pgv->i2c[i+1][j+1][k]-1];
							weight += 1.0;
						}
					}
					if(pgv->i2c[i][j+1][k] > 0){
						for(auto q=0; q<3; q++){
							Vnode[q][i][j][k] += V[q][pgv->ioff+pgv->i2c[i][j+1][k]-1];
							weight += 1.0;
						}
					}
					if(pgv->i2c[i][j][k+1] > 0){
						for(auto q=0; q<3; q++){
							Vnode[q][i][j][k] += V[q][pgv->ioff+pgv->i2c[i][j][k+1]-1];
							weight += 1.0;
						}
					}
					if(pgv->i2c[i+1][j][k+1] > 0){
						for(auto q=0; q<3; q++){
							Vnode[q][i][j][k] += V[q][pgv->ioff+pgv->i2c[i+1][j][k+1]-1];
							weight += 1.0;
						}
					}
					if(pgv->i2c[i+1][j+1][k+1] > 0){
						for(auto q=0; q<3; q++){
							Vnode[q][i][j][k] += V[q][pgv->ioff+pgv->i2c[i+1][j+1][k+1]-1];
							weight += 1.0;
						}
					}
					if(pgv->i2c[i][j+1][k+1] > 0){
						for(auto q=0; q<3; q++){
							Vnode[q][i][j][k] += V[q][pgv->ioff+pgv->i2c[i][j+1][k+1]-1];
							weight += 1.0;
						}
					}
					if(weight >0.1){
						marker[i][j][k] = true;
						for(auto q=0; q<3; q++){
							Vnode[q][i][j][k] /= weight;
						}
					}
					else{
						marker[i][j][k] = true;
					}
				}
			}
		}
		for(auto ic=0; ic<pgv->b->NinX;ic++){
			auto i = pgv->Gn[ic].ijk[0];
			auto j = pgv->Gn[ic].ijk[1];
			auto k = pgv->Gn[ic].ijk[2];
			weight = 0.0;
			Vh[0] = 0.0; Vh[1] = 0.0; Vh[2] = 0.0;
	  	   if(marker[i-pgv->b->imino_][j-pgv->b->jmino_][k-pgv->b->kmino_]){
		    for(auto q=0; q<3; q++){
		    Vh[q]+=Vnode[q][i-1-pgv->b->imino_][j-1-pgv->b->jmino_][k-1-pgv->b->kmino_];
			  weight += 1.0;
		    }
			}
	   	if(marker[i-pgv->b->imino_][j-1-pgv->b->jmino_][k-1-pgv->b->kmino_]){
		    for(auto q=0; q<3; q++){
		     Vh[q]+=Vnode[q][i-pgv->b->imino_][j-1-pgv->b->jmino_][k-1-pgv->b->kmino_];
			  weight += 1.0;
		    }
			}
	   	if(marker[i-pgv->b->imino_][j-pgv->b->jmino_][k-1-pgv->b->kmino_]){
		 	 for(auto q=0; q<3; q++){
		     Vh[q]+=Vnode[q][i-pgv->b->imino_][j-pgv->b->jmino_][k-1-pgv->b->kmino_];
			  weight += 1.0;
		    }
			}
	   	if(marker[i-1-pgv->b->imino_][j-pgv->b->jmino_][k-1-pgv->b->kmino_]){
		 	 for(auto q=0; q<3; q++){
		     Vh[q]+=Vnode[q][i-1-pgv->b->imino_][j-pgv->b->jmino_][k-1-pgv->b->kmino_];
			  weight += 1.0;
		    }
		   }
	   	if(marker[i-1-pgv->b->imino_][j-1-pgv->b->jmino_][k-pgv->b->kmino_]){
		 	 for(auto q=0; q<3; q++){
		     Vh[q]+=Vnode[q][i-1-pgv->b->imino_][j-1-pgv->b->jmino_][k-pgv->b->kmino_];
			  weight += 1.0;
		 	 }
			}
	   	if(marker[i-pgv->b->imino_][j-1-pgv->b->jmino_][k-pgv->b->kmino_]){
		    for(auto q=0; q<3; q++){
		     Vh[q]+=Vnode[q][i-pgv->b->imino_][j-1-pgv->b->jmino_][k-pgv->b->kmino_];
			  weight += 1.0;
		 	 }
			}
	   	if(marker[i-pgv->b->imino_][j-pgv->b->jmino_][k-pgv->b->kmino_]){
		 	 for(auto q=0; q<3; q++){
		     Vh[q]+=Vnode[q][i-pgv->b->imino_][j-pgv->b->jmino_][k-pgv->b->kmino_];
			  weight += 1.0;
		 	 }
			}
	   	if(marker[i-1-pgv->b->imino_][j-pgv->b->jmino_][k-pgv->b->kmino_]){
		 	 for(auto q=0; q<3; q++){
		     Vh[q]+=Vnode[q][i-1-pgv->b->imino_][j-pgv->b->jmino_][k-pgv->b->kmino_];
			  weight += 1.0;
		 	 }
			}
			if(weight > 0.1){
				for(auto q=0; q<3; q++){
					V[q][pgv->ioff+ic] = Vh[q]/weight;
				}
			}
			else{
				for(auto q=0; q<3; q++){
					V[q][pgv->ioff+ic] = 0.0;
				}
			}
		}
		if(pp->ierr==0){
			auto m1 = pgv->b->imax_ -pgv->b->imino_ + 1;
			auto m2 = pgv->b->jmax_ -pgv->b->jmino_ + 1;
			auto m3 = pgv->b->kmax_ -pgv->b->kmino_ + 1;
			for(auto i=0; i<3; i++){
				for(auto j=0; j<m1; j++){
					for(auto k=0; k<m2; k++){
						delete [] Vnode[i][j][k];
					}
					delete [] Vnode[i][j];
				}
				delete [] Vnode[i];
			}
			delete [] Vnode;
			for(auto i=0; i<m1; i++){
				for(auto j=0; j<m2; j++){
					delete [] marker[i][j];
				}
				delete [] marker[i];
			}
			delete[] marker;
		}
		else
			pp->litError("filterVelocity","deallocation error of Vnode, marker");
	}
	pbou->updateGhostR2(V, V_S1, V_S2, pbou->S_GRADIENT);
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinZ; ic++){
		for(auto i=0; i<3; i++) pgv->Gn[ic].V[i] = V[i][pgv->ioff+ic];
	}
}
plbuf->freeR2Buffer(V_S1, V_S2, V);
if(pgv->verbose && pp->myrank==0){
	std::cout<<pgv->clit<<"Starting filterVelocity, nFilter="<<pgv->nVelFilter;
	std::cout<<" ... Done";
}
};

void advection::advection_m_init(){
	std::unique_ptr<timing> pti(new timing);
	pti->lit_timing_create("advection");
};

bool advection::velocityFiltered = false;
