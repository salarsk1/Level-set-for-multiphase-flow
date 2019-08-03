//Written by Salar Safarkhani

#include "bl.h"
#include"arrayND.h"
void bl::bl_init(block_t *bb, std::vector<int> &ijk_sg){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<init>pinit(new init);
bool plus_there, minus_there;
vector<double> xyz(3);
bl_activate_new(bb, ijk_sg);
plus_there = false;
minus_there = false;
for(auto k=bb->kmino_-1; k<bb->kmaxo_; k++){
	for(auto j=bb->jmino_-1; j<bb->jmaxo_; j++){
		for(auto i=bb->imino_-1; i<bb->imaxo_; i++){
			xyz[0] = pgv->xc[i+pgv->nghost];
			xyz[1] = pgv->yc[j+pgv->nghost];
			xyz[2] = pgv->zc[k+pgv->nghost];
			if(0.0<=pinit->G_init_value(xyz)){
				plus_there = true;
			}
			else{
				minus_there = true;
			};
			if(plus_there && minus_there) break;
		}
		if(plus_there && minus_there) break;
	}
	if(plus_there && minus_there) break;
};
if(plus_there && minus_there){
	bb->NinN = 1; 
	bb->NinW = 1;
	bb->NinX = 1;
	bb->NinZ = 1;
}
else{
	bb->NinN = 0; 
	bb->NinW = 0;
	bb->NinX = 0;
	bb->NinZ = 0;
	size_t d1 = bb->imaxo_ - bb->imino_ + 1;
	size_t d2 = bb->jmaxo_ - bb->jmino_ + 1;
	if(pp->ierr==0){
		for(auto i=0;i<d1;i++){
			for(auto j=0;j<d2;j++){
				delete [] bb->ijk2ic[i][j];
			}
			delete [] bb->ijk2ic[i];
		};
		delete [] bb->ijk2ic;
		delete [] bb->Gnodes;
		bb->ijk2ic = nullptr;
		bb->Gnodes = nullptr;
	}
	else
		pp->litError("bl_init","Deallocataion error of pgv->b->ijk2ic,pgv->b->Gnodes in bl_init");
}
};

void bl::bl_clear_all_ghosts(bool old){
std::unique_ptr<global_variable> pgv(new global_variable);
if(old){
	for(auto ibl=0;ibl<pgv->nbl;ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=pgv->b->NinX;ic<pgv->b->NinZ; ic++){
			auto d1 = pgv->b->GnodesOld[ic].ijk[0]-pgv->b->imino_;
			auto d2 = pgv->b->GnodesOld[ic].ijk[1]-pgv->b->jmino_;
			auto d3 = pgv->b->GnodesOld[ic].ijk[2]-pgv->b->kmino_;
			pgv->i2c[d1][d2][d3] = 0;
		}
		pgv->b->NinZ = pgv->b->NinX;
		for(auto i=0; i<26; i++) pgv->b->GnodeGhost_s[i] = 0;
	}
}
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=pgv->b->NinX;ic<pgv->b->NinZ; ic++){
			auto d1 = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
			auto d2 = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
			auto d3 = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
			pgv->i2c[d1][d2][d3] = 0;
		}
		pgv->b->nG = pgv->b->nG - (pgv->b->NinZ-pgv->b->NinX);
		pgv->b->NinZ = pgv->b->NinX;
		for(auto i=0; i<26; i++) pgv->b->GnodeGhost_s[i] = 0;
	}
};
};

void bl::bl_remove_from_list(block_t **bb, block_t **bm, bool first){
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<global_variable> pgv(new global_variable);
if(first){
	pgv->bl = (*bb)->next;
}
else{
	(*bm)->next = (*bb)->next;
};
bl_free(*bb);
if(pp->ierr == 0){
	delete *bb;
	*bb = nullptr;
}
else
	pp->litError("bl_remove_from_list","Deallocation error of bb");
if(first){
	*bb = pgv->bl;
	*bm = pgv->bl;
}
else{
	*bb = (*bm)->next;
};
};
void bl::bl_free(block_t *bb){
std::unique_ptr<parallel>pp(new parallel);
if(bb->Gnodes != nullptr){
	if(pp->ierr==0){
		delete [] bb->Gnodes;
		bb->Gnodes = nullptr;
	}
	else{
		pp->litError("b_free","Deallocation error of bb->Gnodes");
	};
};
if(bb->GnodesOld != nullptr){
	if(pp->ierr == 0){
		delete [] bb->GnodesOld;
		bb->GnodesOld = nullptr;
	}
	else{
		pp->litError("bl_free", "Deallocation error of bb->GnodesOld");
	};
};
if(bb->ijk2ic != 0){
	if(pp->ierr==0){
		size_t d1 = bb->imaxo_ - bb->imino_ + 1;
		size_t d2 = bb->jmaxo_ - bb->jmino_ + 1;
	 	for(auto i=0;i<d1;i++){
	   	for(auto j=0;j<d2;j++){
	      	delete [] bb->ijk2ic[i][j];
	    	}
	   	delete [] bb->ijk2ic[i];
	 	};
	 	delete [] bb->ijk2ic;
		bb->ijk2ic = nullptr;
	}
	else{
		pp->litError("bl_free","Deallocation error of bb->ijk2ic");
	};
}
if(pp->ierr == 0){
	bb->NinBandLayer.release();
	bb->ic0_bl.release();
}
else{
	pp->litError("bl_free","Deallocation error of bb->NinBandLayer,bb->ic0_bl");
};
};
void bl::bl_activate_new(block_t *bb, std::vector<int> &ijk_sg){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
bb->ijkm = pgv->ijkm_bl;
bb->ijk_sg = ijk_sg;
for(auto i=0; i<3; i++){
	bb->ijk0[i] = (ijk_sg[i]-1)*pgv->ijkm_bl[i];
};
bb->imin_ = (ijk_sg[0]-1) * pgv->ijkm_bl[0]+1;
bb->imax_ =  ijk_sg[0] * pgv->ijkm_bl[0];
bb->jmin_ = (ijk_sg[1]-1) * pgv->ijkm_bl[1]+1;
bb->jmax_ =  ijk_sg[1] * pgv->ijkm_bl[1];
bb->kmin_ = (ijk_sg[2]-1) * pgv->ijkm_bl[2]+1;
bb->kmax_ =  ijk_sg[2] * pgv->ijkm_bl[2];
bb->imino_ = bb->imin_ - pgv->nghost;
bb->imaxo_ = bb->imax_ + pgv->nghost;
bb->jmino_ = bb->jmin_ - pgv->nghost;
bb->jmaxo_ = bb->jmax_ + pgv->nghost;
bb->kmino_ = bb->kmin_ - pgv->nghost;
bb->kmaxo_ = bb->kmax_ + pgv->nghost;
bb->imin1_ = bb->imin_ - 1;
bb->imax1_ = bb->imax_ + 1;
bb->jmin1_ = bb->jmin_ - 1;
bb->jmax1_ = bb->jmax_ + 1;
bb->kmin1_ = bb->kmin_ - 1;
bb->kmax1_ = bb->kmax_ + 1;
for(auto i=0; i<3; i++){
	for(auto j=0; j<26; j++){
		bb->ijk_b2g[i][j] = 2*(bb->ijk0[i]+(pgv->dijk_nc[j][i]+1)/2*pgv->ijkm_bl[i]);
	};
};
bb->NinBand = {0};
auto d1 = bb->imaxo_ - bb->imino_ + 1;
auto d2 = bb->jmaxo_ - bb->jmino_ + 1;
auto d3 = bb->kmaxo_ - bb->kmino_ + 1;
if(pp->ierr == 0){
  bb->ijk2ic = new int**[d1];
  for(size_t i=0; i<d1; i++) {
     bb->ijk2ic[i] = new int *[d2];
      for(size_t j=0; j<d2; j++) {
         bb->ijk2ic[i][j] = new int [d3];
      }
   }
}
else{
   pp->litError("bl_activate_new","Allocation error of bb->ijk2ic");
};
for(auto i=0; i<d1; i++){
	for(auto j=0; j<d2; j++){
		for(auto k=0; k<d3; k++){
			bb->ijk2ic[i][j][k] = 0;
		}
	}
}
if(pp->ierr ==0){
	bb->ic0_bl.reset(new int[pgv->nBandLayers]);
	bb->NinBandLayer.reset(new int[pgv->nBandLayers]);
}
else{
	pp->litError("bl_activate_new","Allocation error of bb->ic0_bl");
};
bb->nG = 0;
bb->nGmax = 100;
if(pp->ierr == 0){
	bb->Gnodes = new Gnode_t[bb->nGmax];
}
else{
	pp->litError("bl_activate_new","Allocation error of bb->Gnodes in bl_activate_new");
};
};


