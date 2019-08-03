//Written by Salar Safarkhani

#include<limits>
#include<cmath>
#include"redist.h"
#include<float.h>
void redist::lit_redist_pde(size_t psi_S1, size_t psi_S2, double** psi){
double **n, **psiOld;
double* S;
double myBuf, dtRedist;
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<toolbox> ptool(new toolbox);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<bound> pbou(new bound);
std::array<double, 3> psiMax, psiMin;
std::array<double, 6> buf2, myBuf2;
if(pgv->verbose && pp->myrank==0){
	std::cout<<pgv->clit<<"Starting lit_redist"<<std::endl;
}
auto npsi = psi_S1;
size_t n_S1, n_S2;
n = plbuf->getR2Buffer(n_S1, n_S2);
ptool->normal(n_S1, n_S2, n);
size_t S_S;
S = plbuf->getR1Buffer(S_S);
size_t psiOld_S1, psiOld_S2;
psiOld = plbuf->getR2Buffer(psiOld_S1, psiOld_S2);
for(auto i=0; i<npsi; i++){
	psiMin[i] =  DBL_MAX;
	psiMax[i] = -DBL_MAX;
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinN; ic++){
		S[pgv->ioff+ic]=pgv->Gn[ic].G/sqrt(pow(pgv->Gn[ic].G,2)+pgv->dxyz_min_2);
	}
	for(auto ic=0; ic<pgv->b->NinA+(pgv->b->NinT-pgv->b->NinA)/2; ic++){
		for(auto i=0; i<3; i++){
			psiMin[i] = std::min(psiMin[i], psi[i][pgv->ioff+ic]);
			psiMax[i] = std::max(psiMax[i], psi[i][pgv->ioff+ic]);
		}
	}
}
for(auto ic=0; ic<npsi; ic++){
	myBuf2[2*ic] = psiMin[ic];
	myBuf2[2*ic+1] = -psiMax[ic];
}
pp->ierr=MPI_Allreduce(&myBuf2[0],&buf2[0],6,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
for(auto ic=0; ic<3; ic++){
	psiMin[ic] = buf2[2*ic];
	psiMax[ic] = -buf2[2*ic+1];
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinW; ic++){
		for(auto i=0; i<3; i++){
			auto m1 = std::min(psiMax[i], psi[i][ic]);
			psi[i][ic]=std::max(psiMin[i], m1);
			n[i][ic] *= S[ic];
		}
	}
}
dtRedist = 0.5*pgv->dxyz_min;
for(auto iter=0; iter<max_iter_redist; iter++){
	psiOld = psi;
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=pgv->ioff; ic<pgv->ioff+pgv->b->NinN; ic++){
			auto i = pgv->b->Gnodes[ic-pgv->ioff].ijk[0]-pgv->b->imino_;
			auto j = pgv->b->Gnodes[ic-pgv->ioff].ijk[1]-pgv->b->jmino_;
			auto k = pgv->b->Gnodes[ic-pgv->ioff].ijk[2]-pgv->b->kmino_;
			if(pgv->cylindrical){
				for(auto q=0; i<npsi; i++){
					psi[q][ic]=psiOld[q][ic]-dtRedist*(pgv->rdxyz[0]*(std::max(0.0,\
					n[0][ic])*(psiOld[q][ic]-psiOld[q][pgv->ioff+pgv->i2c[i-1][j][k]-1])+\
					std::min(n[0][ic],0.0)*(psiOld[q][pgv->ioff+pgv->i2c[i+1][j][k]-1]-\
					psiOld[q][ic]))+pgv->rdxyz[2]*(std::max(n[2][ic],0.0)*(psiOld[q][ic]\
					-psiOld[q][pgv->ioff+pgv->i2c[i][j][k-1]-1]))+std::min(n[2][ic],0.0)*\
					(psiOld[q][pgv->ioff+pgv->i2c[i][j][k+1]-1]-psiOld[q][ic]));
				}
			}
			else{
				for(auto q=0; i<npsi; i++){
					psi[q][ic]=psiOld[q][ic]-dtRedist*(pgv->rdxyz[0]*(std::max(0.0,\
					n[0][ic])*(psiOld[q][ic]-psiOld[q][pgv->ioff+pgv->i2c[i-1][j][k]-1])+\
					std::min(n[0][ic],0.0)*(psiOld[q][pgv->ioff+pgv->i2c[i+1][j][k]-1]-\
					psiOld[q][ic]))+pgv->rdxyz[1]*(std::max(n[1][ic],0.0)*(psiOld[q][ic]\
					-psiOld[q][pgv->ioff+pgv->i2c[i][j-1][k]-1]))+std::min(n[1][ic],0.0)*\
					(psiOld[q][pgv->ioff+pgv->i2c[i][j+1][k]-1]-psiOld[q][ic])+\
					pgv->rdxyz[2]*(std::max(n[2][ic],0.0)*(psiOld[q][ic]-psiOld[q]\
					[pgv->ioff+pgv->i2c[i][j][k-1]-1])+std::min(n[2][ic],0.0)*(psiOld[q]\
					[pgv->ioff+pgv->i2c[i][j][k+1]-1]-psiOld[q][ic])));
				}
			}
		}
	}
	pbou->updateGhostR2(psi, psi_S1, psi_S2);
}
plbuf->freeR2Buffer(n_S1, n_S2, n);
plbuf->freeR1Buffer(S);
plbuf->freeR2Buffer(psiOld_S1, psiOld_S2, psiOld);
if(pgv->verbose && pp->myrank==0)
	std::cout<<pgv->clit<<"Starting lit_redist ... Done"<<std::endl;
};

int redist::max_iter_redist = 0;

