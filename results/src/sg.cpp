//Written by Salar Safarkhani

#include"sg.h"
#include"mpi.h"
#include<algorithm>
#include<iterator>
void sg::sg_m_init(){
std::unique_ptr<timing> pti(new timing);
std::unique_ptr<monitor> pmo(new monitor);
pti->lit_timing_create("load_balance");
pti->lit_timing_create("band");
pmo->lit_monitor_create_file_step("lit_balance",4);
pmo->lit_monitor_set_header(1, "imbalance", 'r');
pmo->lit_monitor_set_header(2, "comm/comp", 'r');
pmo->lit_monitor_set_header(3, "min blocks", 'i');
pmo->lit_monitor_set_header(4, "max blocks", 'i');
};

void sg::sg_init(){
std::array<int, 3> ijks_sg_bb, ijke_sg_bb;
bool lskipped;
double imbalance;
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<timing>pti(new timing);
std::unique_ptr<band> pband(new band);
std::unique_ptr<bound>pbou(new bound);
std::unique_ptr<bl>pbl(new bl);
if(pgv->verbose && pp->myrank == 0) 
	std::cout<<pgv->clit<<"Starting sg_init"<<std::endl;
pti->lit_timing_start("band");
if(pp->ierr==0){
	pgv->sg_rank_block = new int**[pgv->ijkm_sg[0]];
	pgv->sg_active = new bool**[pgv->ijkm_sg[0]];
	for(auto i=0; i<pgv->ijkm_sg[0]; i++){
		pgv->sg_rank_block[i] = new int*[pgv->ijkm_sg[1]];
		pgv->sg_active[i] = new bool*[pgv->ijkm_sg[1]];
		for(auto j=0; j<pgv->ijkm_sg[1]; j++){
			pgv->sg_rank_block[i][j] = new int[pgv->ijkm_sg[2]];
			pgv->sg_active[i][j] = new bool[pgv->ijkm_sg[2]];
		}
	}
}
else
	pp->litError("sg_init","allocation error for sg_rank_block, sg_active");
for(auto i=0; i<pgv->ijkm_sg[0]; i++){
	for(auto j=0; j<pgv->ijkm_sg[1]; j++){
		for(auto k=0; k<pgv->ijkm_sg[2]; k++){
			pgv->sg_rank_block[i][j][k] = pgv->sg_notset;
			pgv->sg_active[i][j][k] = false;
		}
	}
}
if(pp->nprocs <= pgv->ijkm_sg[2] && pgv->ijkm_sg[2]%pp->nprocs ==0){
	auto irk = 0;
	for(auto k_sg=0; k_sg<pgv->ijkm_sg[2]; k_sg++){
		for(auto i=0; i<pgv->ijkm_sg[0]; i++){
			for(auto j=0; j<pgv->ijkm_sg[1]; j++){
				pgv->sg_rank_block[i][j][k_sg] = irk;
			}
		}
		irk++;
		if(irk == pp->nprocs) irk = 0;
	}
}
else{
	auto irk = 0;
	auto lskipped = false;
	auto irk_double = 0;
	for(auto i_sg=0; i_sg<pgv->ijkm_sg[0]; i_sg++){
		for(auto j_sg=0; j_sg<pgv->ijkm_sg[1]; j_sg++){
			for(auto k_sg=0; k_sg<pgv->ijkm_sg[2]; k_sg++){
				pgv->sg_rank_block[i_sg][j_sg][k_sg] = irk;
				if(irk==irk_double && !lskipped){
					lskipped = true;
					irk--;
				}
				irk++;
				if(irk==pp->nprocs){
					irk = 0;
					lskipped = false;
					irk_double++;
					if(irk_double == pp->nprocs) irk_double = 0;
				}
			}
		}
	}
}
std::vector<int> arr1(3);
std::vector<int> arr2(3);
arr1 = xyz_to_ijk_sg(pgv->xyzs_init_bb);
arr2 = xyz_to_ijk_sg(pgv->xyze_init_bb);
for(auto i=0; i<3; i++){
	ijks_sg_bb[i] = std::max(1, arr1[i]);
	ijke_sg_bb[i] = std::min(pgv->ijkm_sg[i], arr2[i]);
}
if(pp->ierr==0) pgv->ibl2bl = new block_p_t[pgv->max_bl*2];
else  pp->litError("sg_init","allocation error for ibl2bl");
auto nblh = 0;
if(pp->myrank==0 && pgv->verbose){
	for(auto i=0; i<pgv->ijkm_sg[0]; i++){
		for(auto j=0; j<pgv->ijkm_sg[1]; j++){
			for(auto k=0; k<pgv->ijkm_sg[2]; k++){
				if(pgv->sg_rank_block[i][j][k] == 0) nblh += 1;
			}
		}
	}
}
auto irk = 0;
auto ibl = 0;
for(auto k_sg=ijks_sg_bb[2]; k_sg<=ijke_sg_bb[2]; k_sg++){
	for(auto j_sg=ijks_sg_bb[1]; j_sg<=ijke_sg_bb[1]; j_sg++){
		for(auto i_sg=ijks_sg_bb[0]; i_sg<=ijke_sg_bb[0]; i_sg++){
			if(pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1] == pp->myrank){
				ibl += 1;
				pgv->sg_rank_block[i_sg-1][j_sg-1][k_sg-1] = -ibl;
				if(ibl == 1){
					if(pp->ierr==0) pgv->bl = new block_t;
					else pp->litError("sg_init","allocation error for bl");
					pgv->ibl2bl[ibl-1].p = pgv->bl;
				}
				else{
					if(pp->ierr==0)pgv->ibl2bl[ibl-1].p = new block_t;
					else pp->litError("sg_init","allocation error for ibl2bl[ibl].p");
					pgv->ibl2bl[ibl-2].p->next = pgv->ibl2bl[ibl-1].p;
				}
				std::vector<int> vec_temp{i_sg,j_sg, k_sg};
				pbl->bl_init(pgv->ibl2bl[ibl-1].p, vec_temp);
				i_sg = vec_temp[0]; j_sg = vec_temp[1]; k_sg = vec_temp[2];
			}
		}
	}
}
pgv->nbl = ibl;
sg_calc_active();
pband->band_init();
sg_calc_active();
pti->lit_timing_stop("band");
sg_load_balance(true,'y');
pbou->prepareGhostNodes();
pbou->updateGhostNodes();
if(pgv->verbose && pp->myrank == 0){
	std::cout<<pgv->clit<<'\t'<<"Starting sg_init ... Done"<<std::endl;
}
};

void sg::sg_band_update(){
std::unique_ptr<bl>pbl(new bl);
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<band>pband(new band);
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<timing>pti(new timing);
std::unique_ptr<bound> pbou(new bound);
if(pgv->verbose && pp->myrank==0)
	std::cout<<pgv->clit<<"Starting band update"<<std::endl;
	pti->lit_timing_start("band");
	pband->band_seed();
	pbl->bl_clear_all_ghosts(true);
	pband->band_grow(false, 'n');
	pband->band_set_block_data();
	pband->band_regenerate_ijk2ic();
	sg_calc_active();
	pband->band_cleanup();
	pbou->prepareGhostNodes();
	pbou->updateGhostNodes();
	pti->lit_timing_stop("band");
	if(pgv->verbose && pp->myrank==0)
		std::cout<<pgv->clit<<"Starting band update ... Done"<<std::endl;
};

void sg::sg_calc_active(){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<bl>pbl(new bl);
bool first, change;
block_t *bb, *bm;

if(pgv->verbose && pp->myrank==0) std::cout<<pgv->clit<<"Starting sg_calc_active"<<std::endl;
for(auto i=0; i<pgv->ijkm_sg[0]; i++){
	for(auto j=0; j<pgv->ijkm_sg[1]; j++){
		for(auto k=0; k<pgv->ijkm_sg[2]; k++){
			pgv->sg_active[i][j][k] = false;
		}
	}
}

bb = pgv->bl;
bm = pgv->bl;

first = true;
change = false;
auto count = 0;
while(bb != nullptr){
	if(bb->NinX > 0){
		pgv->sg_active[bb->ijk_sg[0]-1][bb->ijk_sg[1]-1][bb->ijk_sg[2]-1] = true;
		bb = bb->next;
		if(!first) bm = bm->next;
		else first = false;
	}
	else{
	 change = true;
	 pgv->sg_rank_block[bb->ijk_sg[0]-1][bb->ijk_sg[1]-1][bb->ijk_sg[2]-1]=pp->myrank;
	 pbl->bl_remove_from_list(&bb,&bm,first);
	}
}
if(change){
	for(auto ibl=0; ibl<pgv->nbl; ibl++) pgv->ibl2bl[ibl].p = nullptr;
	bb = pgv->bl;
	auto ibl = 0;
	while(bb != nullptr){
		ibl++;
		pgv->ibl2bl[ibl-1].p = bb;
		pgv->sg_rank_block[bb->ijk_sg[0]-1][bb->ijk_sg[1]-1][bb->ijk_sg[2]-1]=-ibl;
		bb = bb->next;
	}
	pgv->nbl = ibl;
}
pp->parallel_sync_sg_active();
pp->parallel_all_sum(pgv->nbl,pgv->nbl_gl);
if(pgv->verbose && pp->myrank==0){
	std::cout<<pgv->clit<<"number of active blocks = "<<pgv->nbl_gl<<std::endl;
	std::cout<<pgv->clit<<"Starting sg_calc_active ... Done"<<std::endl;
}
};

void sg::sg_load_balance(bool force, const char yn){
double imbalance;
bool lforce;
std::unique_ptr<timing>pti(new timing);
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<global_variable>pgv(new global_variable);
pti->lit_timing_start("load_balance");
monitor_load(imbalance);
if(pp->nprocs == 1){
	pti->lit_timing_stop("load_balance");
	return;
}
if(yn == 'y') lforce = force;
else lforce = false;
if(!lforce){
	if(imbalance < pp->max_load_imbalance) {
		pti->lit_timing_stop("load_balance");
		if(pp->myrank==0 && pgv->verbose){
			std::cout<<pgv->clit<<"Ending sg_load_balance ...";
			std::cout<<"imbalance,imbalanceTrigger"<<imbalance;
			std::cout<<pp->max_load_imbalance<<std::endl;
			return;
		}
	}
	if(pp->myrank==0 && pgv->verbose){
		std::cout<<pgv->clit<<"Starting sg_load_balance ...,imbalance";
		std::cout<<"imbalanceTrigger="<<imbalance<<pp->max_load_imbalance<<std::endl;
	}
	if(pgv->loadBalancer=="lit"||pgv->loadBalancer=="<Lit"||pgv->loadBalancer=="LIT"){
		sg_load_balance_lit();
	}
	else if(pgv->loadBalancer=="parmetis" || pgv->loadBalancer=="Parmetis" ||\
			  pgv->loadBalancer=="PARMETIS"){
		sg_load_balance_parmetis(lforce);
	}
	else if(pgv->loadBalancer=="parmetis_fix"||pgv->loadBalancer=="Parmetis_fix"||\
			  pgv->loadBalancer=="PARMETIS_FIX"){
		sg_load_balance_parmetis_fix(lforce);
	}
	else if(pgv->loadBalancer=="metis"||pgv->loadBalancer=="Metis"||\
			  pgv->loadBalancer=="METIS"){
		sg_load_balance_metis(lforce);
	}
	else
		pp->litError("sg_load_balance","inknown load balancer");
	monitor_load(imbalance);
	if(pgv->verbose && pp->myrank==0){
		std::cout<<pgv->clit<<"post-load balance node-imbalance";
		std::cout<<imbalance<<std::endl;
		std::cout<<pgv->clit<<"Starting sg_laod_balance ... Done"<<std::endl;
	}
}
};

void sg::sg_load_balance_metis(bool force){
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<global_variable> pgv(new global_variable);
bool do_bl_list_update;
int nbln, ibls, iedge, nedge, edgecut, wgtflag, iblNeighbor, ipe, nex;
int myNodeCount, cch1, cch2, ch3, cch4, myCommCount;
vector<int> nbl_pe(pp->nprocs), nbl_pe_true(pp->nprocs), iblOffset(pp->nprocs);
vector<bool> remove_bl(pgv->max_bl);
vector<int> myNodeC(pp->nprocs), NodeC(pp->nprocs), myCommC(pp->nprocs);
vector<int> CommC(pp->nprocs);
vector<int> myWeight, weight, xadj, adjwgt, adjncy, vwgt, part;
vector<vector<int>> bijksg, mybijksg, myEdgeCount, edgeCount, exScript;
int*** blgl;
std::array<int, 5> options;
int ncon, nblgl;
std::vector<int> ijk_sg_n(3), ijk(3), ijk_shift(3), ijkNeighbor(3), ihelp(3);
vector<double> tpwgts, ubvec;
block_t* bb;
if(pp->myrank==0 && pgv->verbose)
	std::cout<<pgv->clit<<"Starting metis load balancer"<<std::endl;
pp->parallel_all_sum(pgv->nbl, nblgl);
if(nblgl < pp->nprocs){
	sg_load_balance_lit();
	return;
}
pp->ierr = MPI_Gather(&pgv->nbl, 1, MPI_INTEGER,&nbl_pe[0], 1, MPI_INTEGER,\
							 0,MPI_COMM_WORLD);
if(pp->myrank==0){
	iblOffset[0] = 0;
	for(auto irk=1; irk<pp->nprocs; irk++){
		iblOffset[irk] = iblOffset[irk-1] + nbl_pe[irk-1];
	}
	nblgl = iblOffset[pp->nprocs-1]+nbl_pe[pp->nprocs-1];
}
else{
nblgl = 1;
}
if(pp->ierr==0){
	mybijksg.resize(3);
	bijksg.resize(3);
	for(auto i=0; i<3; i++){
		mybijksg[i].resize(pgv->nbl);
		bijksg[i].resize(nblgl);
	}
}
else pp->litError("sg_load_balance","allocation error for mybijksg,bijksg");
for(auto i=0; i<mybijksg.size(); i++){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		mybijksg[i][ibl] = pgv->ibl2bl[ibl].p->ijk_sg[i];
	}
}
auto recvcounts = nbl_pe;
auto displs     = iblOffset;
for(auto &i:recvcounts) i *= 3;
for(auto &i:displs) i *= 3;
pp->ierr=MPI_Gatherv(mybijksg.data(),3*pgv->nbl,MPI_INTEGER,bijksg.data(),\
							recvcounts.data(),displs.data(),MPI_INTEGER,0,MPI_COMM_WORLD);
if(pp->ierr==0){
	mybijksg.clear();
	mybijksg.shrink_to_fit();
}
else pp->litError("sg_load_balance","allocation error for mybijksg");
if(pp->myrank==0){
	if(pp->ierr==0){
		blgl = new int**[pgv->ijkm_sg[0]];
		for(auto i=0; i<pgv->ijkm_sg[0]; i++){
			blgl[i] = new int*[pgv->ijkm_sg[1]];
			for(auto j=0; j<pgv->ijkm_sg[1]; j++){
				blgl[i][j] = new int[pgv->ijkm_sg[2]];
			}
		}
	}
	else pp->litError("sg_load_balance","allocation error for blgl");
	for(auto ibl=0; ibl<nblgl; ibl++){
		blgl[bijksg[0][ibl]-1][bijksg[1][ibl]-1][bijksg[2][ibl]-1] = ibl+1;
	}
}
if(pp->ierr==0){
	myWeight.resize(pgv->nbl);
	weight.resize(nblgl);
	xadj.resize(nblgl+1);
	myEdgeCount.resize(26);
	edgeCount.resize(26);
	for(auto i=0; i<26; i++){
		myEdgeCount[i].resize(pgv->nbl);
		edgeCount[i].resize(nblgl);
	}
}
else pp->litError("sg_load_balance","allocation error for myWeight,myEdgeCount...");
// line 515, 516; checked
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	switch(pgv->loadBalanceband){
		case 1:
			myWeight[ibl] = pgv->b->NinA;
			break;
		case 2:
			myWeight[ibl] = pgv->b->NinT;
			break;
		case 3:
			myWeight[ibl] = pgv->b->NinN;
			break;
		case 4:
			myWeight[ibl] = pgv->b->NinW;
			break;
		case 5:
			myWeight[ibl] = pgv->b->NinX;
			break;
	}
	for(auto i=0; i<myEdgeCount.size(); i++) 
		myEdgeCount[i][ibl]=pgv->b->nGnodeb[i];
}
pp->ierr=MPI_Gatherv(myWeight.data(),pgv->nbl,MPI_INTEGER,weight.data(),nbl_pe.data(),\
							iblOffset.data(),MPI_INTEGER,0,MPI_COMM_WORLD);

recvcounts = nbl_pe;
displs     = iblOffset;
for(auto &i:recvcounts) i *= 26;
for(auto &i:displs)     i *= 26;
pp->ierr=MPI_Gatherv(myEdgeCount.data(),26*pgv->nbl,MPI_INTEGER,edgeCount.data(),\
							recvcounts.data(),displs.data(),MPI_INTEGER,0,MPI_COMM_WORLD);
if(pp->myrank==0){
	for(auto i=0; i<weight.size(); i++){
		weight[i] = std::max(1, weight[i]); // checked
	}
	nedge = 0;
	xadj[0] = 0;
	for(auto ibl=0; ibl<nblgl; ibl++){
		for(auto ibf=0; ibf<26; ibf++){
			vector<int> v1, v2;
			for(auto i=0; i<bijksg.size(); i++) v1.push_back(bijksg[i][ibl]);
			for(auto i=0;i<3; i++) v2.push_back(pgv->dijk_nc[ibf][i]);
			pgv->shift_ijk_sg(v1,v2,ijkNeighbor,ihelp);
			bool check = false;
			for(auto i=0; i<ijkNeighbor.size(); i++){
				if(ijkNeighbor[i] == 0) check = true;
				break;
			}
			if(check) continue;
			if(ijkNeighbor[0]>pgv->ijkm_sg[0]) continue;
			if(ijkNeighbor[1]>pgv->ijkm_sg[1]) continue;
			if(ijkNeighbor[2]>pgv->ijkm_sg[2]) continue;
			if(!(pgv->sg_active[ijkNeighbor[0]-1][ijkNeighbor[1]-1][ijkNeighbor[2]-1])){
				continue;
			}
			iblNeighbor = blgl[ijkNeighbor[0]-1][ijkNeighbor[1]-1][ijkNeighbor[2]-1];
			if(edgeCount[ibf][ibl]+edgeCount[pgv->ibf2[ibf]-1][iblNeighbor-1]>0) nedge++;
		}
		xadj[ibl+1] = nedge;
	}
	if(pp->ierr==0){
		adjwgt.resize(nedge);
		adjncy.resize(nedge);
	}
	else pp->litError("sg_load_balance", "allocation error for adjwgt,adjncy");
	iedge = 0;
	for(auto ibl=0; ibl<nblgl; ibl++){
		for(auto ibf=0; ibf<26; ibf++){
			vector<int> v1, v2;
			for(auto i=0; i<bijksg.size(); i++) v1.push_back(bijksg[i][ibl]);
			for(auto i=0;i<3; i++) v2.push_back(pgv->dijk_nc[ibf][i]);
			pgv->shift_ijk_sg(v1,v2,ijkNeighbor,ihelp);
			bool check = false;
			for(auto i=0; i<ijkNeighbor.size(); i++){
				if(ijkNeighbor[i] == 0) check = true;
				break;
			}
			if(check) continue;
			if(ijkNeighbor[0]>pgv->ijkm_sg[0]) continue;
			if(ijkNeighbor[1]>pgv->ijkm_sg[1]) continue;
			if(ijkNeighbor[2]>pgv->ijkm_sg[2]) continue;
			if(!(pgv->sg_active[ijkNeighbor[0]-1][ijkNeighbor[1]-1][ijkNeighbor[2]-1])){
				continue;
			}
			iblNeighbor = blgl[ijkNeighbor[0]-1][ijkNeighbor[1]-1][ijkNeighbor[2]-1];
			if(edgeCount[ibf][ibl]+edgeCount[pgv->ibf2[ibf]-1][iblNeighbor-1]>0){
				iedge++;
				adjncy[iedge-1] = iblNeighbor-1;
				adjwgt[iedge-1] = edgeCount[ibf][ibl]+edgeCount[pgv->ibf2[ibf]-1]\
										[iblNeighbor-1];
			}
		}
	}
	if(pp->ierr==0){
		for(auto i=0; i<pgv->ijkm_sg[0]; i++){
			for(auto j=0; j<pgv->ijkm_sg[1]; j++){
				delete [] blgl[i][j];
			}
			delete [] blgl[i];
		}
		delete [] blgl;
	}
	else pp->litError("sg_load_balance","deallocation error for blgl");
	wgtflag = 3;
	if(pp->ierr==0) part.resize(nblgl);
	else pp->litError("sg_load_balance","allocation error for part");
	//The c++ version of the following metis function needs 13 varivals;
	//In fortran only 11 is provided. The last 2 arguments is nullptr;
	//It should be corrected;
//	METIS_PartGraphKway(&nblgl, &xadj[0], &adjncy[0],&weight[0],&adjwgt[0],\
							  &wgtflag,nullptr, &(pp->nprocs), reinterpret_cast<float*>\
							  (&options[0]), reinterpret_cast<float*>(&edgecut),&part[0],\
							  nullptr, nullptr);
	if(pp->ierr==0){
		adjwgt.clear(); adjwgt.shrink_to_fit();
		adjncy.clear(); adjncy.shrink_to_fit();
		xadj.clear(); xadj.shrink_to_fit();
	}
	else pp->litError("sg_load_balance_metis","deallocation error for adjwgt");
	if(pp->ierr==0){
		exScript.resize(6);
		for(auto i=0; i<6; i++) exScript[i].resize(nblgl);
	}
	else pp->litError("sg_load_balance_metis","allocation error for exScript");
	auto ipe = 0;
	auto nex = 0;
	for(auto ibl=0; ibl<nblgl; ibl++){
		while(ibl+1>iblOffset[ipe]+nbl_pe[ipe]) ipe++; // checked
		if(ipe != part[ibl]){
			nex++;
			exScript[0][nex-1] = ipe;
			exScript[1][nex-1] = part[ibl];
			exScript[2][nex-1] = ibl+1-iblOffset[ipe];
			exScript[3][nex-1] = bijksg[0][ibl];
			exScript[4][nex-1] = bijksg[1][ibl];
			exScript[5][nex-1] = bijksg[2][ibl];
		}
	}
	if(pp->ierr==0){
		part.clear(); part.shrink_to_fit();
	}
	else pp->litError("sg_load_balance_metis","deallocation error of part");
}
pp->ierr = MPI_Bcast(&nex, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
if(pp->myrank==0){
	if(pp->ierr==0){
		exScript.resize(6);
		for(auto i=0; i<6; i++) exScript[i].resize(nex);
	}
	else pp->litError("sg_load_balance_metis","allocation error for exScript");
}
pp->ierr==MPI_Bcast(&exScript[0][0],6*nex, MPI_INTEGER, 0, MPI_COMM_WORLD);
executeBlockExchange(nex, exScript);
};

void sg::sg_load_balance_parmetis(bool force){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bl> pbl(new bl);
bool do_bl_list_update;
vector<bool> remove_bl(pgv->max_bl);
int nbln, iedge, nedge, edgecut, wgtflag;
vector<int> nbl_pe(pp->nprocs), nbl_pe_true(pp->nprocs), vtxdist(pp->nprocs+1);
int myNodeCount, cch1, cch2, cch3, cch4, myCommCount;
vector<int> myNodeC(pp->nprocs), NodeC(pp->nprocs), myCommC(pp->nprocs);
vector<int> CommC(pp->nprocs);
int ***blgl, ***myblgl;
vector<int> partDest, partSrc, xadj, adjwgt, adjncy, vwgt, vwgt2, mypart;
vector<vector<int>> bijksg, mybijksg, nGnoder;
std::array<int, 5> options;
int ncon;
std::array<int, 3> ijk_sg_n, ijk, ijk_shift;
vector<double> tpwgts, ubvec;
block_t* bb;
if(pp->myrank=0 && pgv->verbose)
	std::cout<<pgv->clit<<"Starting parmetis load balancer"<<std::endl;
pp->ierr = MPI_Allgather(&(pgv->nbl),1,MPI_INTEGER,nbl_pe.data(),1,MPI_INTEGER,MPI_COMM_WORLD);
bool check=false;
for(auto i=0; i<nbl_pe.size(); i++){
	if(nbl_pe[i] == 0){
		check = true;
		break;
	}
}
if(check){
	sg_load_balance_metis(force);
	return;
}
vtxdist[0] = 1;
for(auto irk=0; irk<pp->nprocs; irk++){
	vtxdist[irk+1] = vtxdist[irk] + nbl_pe[irk];
}
myblgl = new int**[pgv->ijkm_sg[0]]; blgl = new int**[pgv->ijkm_sg[0]];
if(pp->ierr==0){
	for(auto i=0; i<pgv->ijkm_sg[0]; i++){
		myblgl[i] = new int*[pgv->ijkm_sg[1]];
		blgl[i]   = new int*[pgv->ijkm_sg[1]];
		for(auto j=0; j<pgv->ijkm_sg[1]; j++){
			myblgl[i][j] = new int[pgv->ijkm_sg[2]];
			blgl[i][j]   = new int[pgv->ijkm_sg[2]];
		}
	}
}
else pp->litError("sg_load_balance","allocation error for myblgl,blgl");
if(pp->ierr==0){
	mybijksg.resize(3);
	bijksg.resize(3);
	for(auto i=0; i<3; i++){
		mybijksg[i].resize(pgv->nbl);
		bijksg[i].resize(vtxdist[pp->nprocs]);
	}
}
else pp->litError("allocation error","mybijksg");
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	myblgl[pgv->b->ijk_sg[0]-1][pgv->b->ijk_sg[1]-1][pgv->b->ijk_sg[2]-1] = \
	vtxdist[pp->myrank]+ibl;
	for(auto i=0; i<mybijksg.size(); i++) mybijksg[i][ibl] = pgv->b->ijk_sg[i];
}
auto d=pgv->b->ijk_sg[0]*pgv->b->ijk_sg[1]*pgv->b->ijk_sg[2];
pp->ierr=MPI_Allreduce(myblgl,blgl,d,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
if(pp->ierr==0){
	for(auto i=0; i<pgv->b->ijk_sg[0]; i++){
		for(auto j=0; j<pgv->b->ijk_sg[1]; j++){
			delete [] myblgl[i][j];
		}
		delete [] myblgl[i];
	}
	delete [] myblgl;
}
else pp->litError("sg_load_balance","deallocation error for myblgl");
if(pp->ierr=0){
	xadj.resize(pgv->nbl+1);
	vwgt.resize(pgv->nbl);
	nGnoder.resize(26);
	for(auto i=0; i<26; i++){
		nGnoder[i].resize(pgv->nbl);
	}
}
else pp->litError("sg_load_balance","allocation error for xadj");
iedge = 0;
myCommCount = 0;
myNodeCount = 0;
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	switch(pgv->loadBalanceband){
		case 1:
			vwgt[ibl] = pgv->b->NinA;
			break;
		case 2:
			vwgt[ibl] = pgv->b->NinT;
			break;
		case 3:
			vwgt[ibl] = pgv->b->NinN;
			break;
		case 4:
			vwgt[ibl] = pgv->b->NinW;
			break;
		case 5:
			vwgt[ibl] = pgv->b->NinX;
			break;
	}
	myNodeCount += static_cast<int>(vwgt[ibl]);
	xadj[ibl] = iedge+1;
	for(auto ibf=0; ibf<26; ibf++){
		if(!(pgv->b->boundType[ibf]==pgv->bNeumann)){
			auto counter = 0;
			for(auto i=pgv->b->ijk0[0]+pgv->igs[ibf]-pgv->b->imino_;\
				 i<=pgv->b->ijk0[0]+pgv->ige[ibf]-pgv->b->imino_; i++){
				for(auto j=pgv->b->ijk0[1]+pgv->jgs[ibf]-pgv->b->jmino_;\
					 j<=pgv->b->ijk0[1]+pgv->jge[ibf]-pgv->b->jmino_; j++){
					for(auto k=pgv->b->ijk0[2]+pgv->kgs[ibf]-pgv->b->kmino_;\
						 k<=pgv->b->ijk0[2]+pgv->kge[ibf]-pgv->b->kmino_; k++){
						if(pgv->b->ijk2ic[i][j][k]>0) counter++; 

					}
				}
			}
			nGnoder[ibf][ibl] = counter;
		}
		else nGnoder[ibf][ibl] = 0;
		if(pgv->b->nGnodeb[ibf]+nGnoder[ibf][ibl]>0){
			iedge++;
			if(pgv->sg_rank_block[pgv->b->ijk_sg[0]+pgv->dijk_nc[ibf][0]-1]\
				[pgv->b->ijk_sg[1]+pgv->dijk_nc[ibf][1]-1]\
				[pgv->b->ijk_sg[2]+pgv->dijk_nc[ibf][2]-1] >= 0){
				myCommCount += static_cast<int>(pgv->b->nGnodeb[ibf]);
			}
		}
	}
}
xadj[pgv->nbl] = iedge+1;
nedge = iedge;
if(pp->ierr==0){
	adjwgt.resize(nedge);
	adjncy.resize(nedge);
}
else pp->litError("sg_load_balance","allocation error for adjwgt,adjncy");
iedge = 0;
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ibf=0; ibf<26; ibf++){
		if(pgv->b->nGnodeb[ibf]+nGnoder[ibf][ibl] > 0){
			iedge++;
			adjwgt[iedge-1] = pgv->b->nGnodeb[ibf]+nGnoder[ibf][ibl];
			adjncy[iedge-1] = blgl[pgv->b->ijk_sg[0]+pgv->dijk_nc[ibf][0]-1]\
									  [pgv->b->ijk_sg[1]+pgv->dijk_nc[ibf][1]-1]\
									  [pgv->b->ijk_sg[2]+pgv->dijk_nc[ibf][2]-1];
		}
	}
}
ncon=2;
wgtflag=3;
if(pp->ierr==0){
	tpwgts.resize(ncon*pp->nprocs);
	ubvec.resize(ncon);
	vwgt2.resize(pgv->nbl*ncon);
	mypart.resize(pgv->nbl);
}
else pp->litError("sg_load_balance","allocation error for tpwgts,ubvec");
std::fill(tpwgts.begin(),tpwgts.begin()+tpwgts.size(),1.0/static_cast<double>(ncon*pp->nprocs));
std::fill(ubvec.begin(),ubvec.begin()+ubvec.size(),1.05);
// line 780 checked
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	vwgt2[2*ibl]=1;
	vwgt2[2*ibl+1] = std::max(vwgt[ibl]/1000,1);
}
for(auto i=0; i<adjwgt.size(); i++) adjwgt[i] = std::max(adjwgt[i]/1000,1);
//int* ad;//check it;
//ParMETIS_V3_PartKway(&vtxdist[0],&xadj[0],&adjncy[0],&vwgt2[0],&adjwgt[0],&wgtflag,\
							ad,&ncon,&(pp->nprocs),reinterpret_cast<float*>(&tpwgts[0]),\
							reinterpret_cast<float*>(&ubvec[0]),&options[0],&edgecut,\
							&mypart[0],&(pp->MY_G_WORLD));
if(pp->ierr==0){
	xadj.clear(); xadj.shrink_to_fit();
	adjwgt.clear(); adjwgt.shrink_to_fit();
	adjncy.clear(); adjncy.shrink_to_fit();
	vwgt.clear(); vwgt.shrink_to_fit();
	vwgt2.clear(); vwgt2.shrink_to_fit();
	tpwgts.clear(); tpwgts.shrink_to_fit();
	ubvec.clear(); ubvec.shrink_to_fit();
}
else pp->litError("sg_load_balance","deallocation error for xadj, ...");
if(pp->ierr==0){
	partDest.resize(vtxdist[pp->nprocs]);
	partSrc.resize(vtxdist[pp->nprocs]);
}
else pp->litError("sg_load_balance","allocation error for partDest, partSrc");
for(auto i=0; i<pp->nprocs; i++){
	for(auto j=vtxdist[i]-1;j<vtxdist[i+1]-1; j++){
		partSrc[j] = i;
	}
}
for(auto i=0; i<vtxdist.size(); i++) vtxdist[i]--;
for(auto i=0; i<mypart.size(); i++) mypart[i]--;
pp->ierr=MPI_Allgatherv(mypart.data(),pgv->nbl,MPI_INTEGER,partDest.data(),\
								nbl_pe.data(),vtxdist.data(),MPI_INTEGER, MPI_COMM_WORLD);
for(auto i=0; i<vtxdist.size(); i++) vtxdist[i] *= 3;
for(auto i=0; i<nbl_pe.size(); i++) nbl_pe[i] *= 3;
pp->ierr=MPI_Allgatherv(mybijksg.data(),3*pgv->nbl,MPI_INTEGER,bijksg.data(),\
								nbl_pe.data(),vtxdist.data(),MPI_INTEGER, MPI_COMM_WORLD);
for(auto i=0; i<vtxdist.size(); i++) vtxdist[i] = vtxdist[i]/3 + 1;
for(auto i=0; i<nbl_pe.size(); i++) nbl_pe[i] /= 3;
pp->ierr=MPI_Reduce(&myNodeCount,&cch1,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD);
pp->ierr=MPI_Reduce(&myNodeCount,&cch2,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD);
pp->ierr=MPI_Reduce(&myNodeCount,&cch3,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD);
if(pp->myrank==0 && pgv->verbose){
	auto d = static_cast<double>(cch1)/static_cast<double>(pp->nprocs);
	d = static_cast<double>(cch2)/d;
	std::cout<<pgv->clit<<"pre-load balance node-imbalance:"<<d<<std::endl;
	std::cout<<pgv->clit<<"pre-load balance # of communication nodes: ";
	std::cout<<cch3/2<<std::endl;
	std::cout<<pgv->clit<<"pre-load balance ratio communication/computenodes: ";
	std::cout<<static_cast<double>(0.5*cch3)/static_cast<double>(cch1);
}
std::fill(myNodeC.begin(), myNodeC.begin()+myNodeC.size(), 0);
std::fill(myCommC.begin(), myCommC.begin()+myNodeC.size(), 0);
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	switch(pgv->loadBalanceband){
		case 1:
			myNodeC[mypart[ibl]] += static_cast<int>(pgv->b->NinA);
			break;
		case 2:
			myNodeC[mypart[ibl]] += static_cast<int>(pgv->b->NinT);
			break;
		case 3:
			myNodeC[mypart[ibl]] += static_cast<int>(pgv->b->NinN);
			break;
		case 4:
			myNodeC[mypart[ibl]] += static_cast<int>(pgv->b->NinW);
			break;
		case 5:
			myNodeC[mypart[ibl]] += static_cast<int>(pgv->b->NinX);
			break;
	}
	for(auto ibf=0; ibf<26; ibf++){
		if(pgv->b->nGnodeb[ibf]+nGnoder[ibf][ibl]>0){
			auto i = pgv->b->ijk_sg[0]+pgv->dijk_nc[ibf][0] - 1;
			auto j = pgv->b->ijk_sg[1]+pgv->dijk_nc[ibf][1] - 1;
			auto k = pgv->b->ijk_sg[2]+pgv->dijk_nc[ibf][2] - 1;
			if(mypart[ibl] != partDest[blgl[i][j][k]]){
				myCommC[mypart[ibl]] += static_cast<int>(pgv->b->nGnodeb[ibf]);
			}
		}
	}
}
if(pp->ierr==0){
	nGnoder.clear();
	mypart.clear();
	nGnoder.shrink_to_fit();
	mypart.shrink_to_fit();
}
else pp->litError("sg_load_balance","deallocation error for nGnoder,mypart");
pp->ierr=MPI_Reduce(myNodeC.data(),NodeC.data(),pp->nprocs, MPI_INTEGER, MPI_SUM,\
						  0, MPI_COMM_WORLD);
pp->ierr=MPI_Reduce(myCommC.data(),CommC.data(),pp->nprocs, MPI_INTEGER, MPI_SUM,\
						  0, MPI_COMM_WORLD);
if(pp->myrank==0 && pgv->verbose){
	auto sum = 0;
	for (auto i:NodeC) sum += i;
	auto d1 = static_cast<double>(sum)/static_cast<double>(pp->nprocs);
	auto d2 = std::max_element(NodeC.begin(),NodeC.begin()+NodeC.size());
	auto d3 = *d2;
	std::cout<<pgv->clit<<"post-load balance node-imbalance: "<<d3/d1<<std::endl;
	sum = 0;
	for(auto i:CommC) sum += i/2;
	std::cout<<pgv->clit<<"post-load balance # of communication nodes: ";
	std::cout<<sum<<std::endl;
	sum = 0;
	for(auto i:CommC) sum += i;
	auto sum2 = 0;
	for(auto i:NodeC) sum2 += i;
	auto d = static_cast<double>(sum)/static_cast<double>(sum2);
	std::cout<<pgv->clit<<"post-load balance ratio communication/compute nodes";
	std::cout<<0.5*d<<std::endl;
}
std::fill(remove_bl.begin(), remove_bl.begin()+remove_bl.size(),false);
do_bl_list_update = false;
for(auto i=0; i<vtxdist[pp->nprocs]-1; i++){
	if(partSrc[i] == pp->myrank){
		if(partDest[i] != pp->myrank){
			auto ibl = i+1-vtxdist[pp->myrank];
			bb = pgv->ibl2bl[ibl].p;
			pp->parallel_block_send(bb, partDest[i]);
			pgv->sg_rank_block[bb->ijk_sg[0]-1][bb->ijk_sg[1]-1][bb->ijk_sg[2]-1]=\
									 partDest[i];
			pbl->bl_free(bb);
			remove_bl[ibl] = true;
			do_bl_list_update = true;
		}
	}
	else if(partDest[i] == pp->myrank){
		pgv->nbl++;
		auto ibl = pgv->nbl-1;
#ifdef DEBUG_MODE
		if(pgv->ibl2bl[ibl].p != nullptr){
			std::cout<<pp->cmy<<"trying to allocate block already exist"<<std::endl;
			pp->parallel_kill(1);
		}
#endif
		if(pp->ierr==0) pgv->ibl2bl[ibl].p = new block_t;
		else pp->litError("sg_load_balance","allocation error for ibl2bl[ibl].p");
		bb = pgv->ibl2bl[ibl].p;
		pp->parallel_block_recv(bb, partSrc[i]);
		do_bl_list_update = true;
		pgv->sg_rank_block[bb->ijk_sg[0]-1][bb->ijk_sg[1]-1][bb->ijk_sg[2]-1]=-ibl-1;
#ifdef DEBUG_MODE
		if(!pgv->sg_active[bb->ijk_sg[0]-1][bb->ijk_sg[1]-1][bb->ijk_sg[2]-1]){
			std::cout<<pgv->clit<<"["<<pp->myrank<<"] error:will receive block locally is marked as not a active"<<std::endl;
			pp->parallel_kill(23);
		}
#endif
	}
	else{
		auto sum = 0;
		for(auto q=0; q<3; q++) sum+=bijksg[q][i];
		if(sum > 0)
			pgv->sg_rank_block[bijksg[0][i]-1][bijksg[1][i]-1][bijksg[2][i]-1] = partDest[i];
	}
}
if(do_bl_list_update){
	for(auto ibl=0; ibl<pgv->nbl-1; ibl++){
		if(remove_bl[ibl]){
			if(pp->ierr==0) pgv->ibl2bl[ibl].p = nullptr;
			else pp->litError("sg_load_balance","deallocation error for ibl2bl[ibl].p");
		}
	}
	auto ibl = 0;
	while(remove_bl[ibl] && ibl<pgv->nbl) ibl++;
	if(remove_bl[ibl]){
		pgv->bl = nullptr;
		std::cout<<pp->cmy<<"What am I doing here?"<<std::endl;
		nbln = 0;
	}
	else{
		pgv->bl = pgv->ibl2bl[ibl].p;
		bb = pgv->bl;
		bb->next = nullptr;
		nbln = 0;
		auto ibls = ibl+1;
		for(auto ibl=ibls; ibl<pgv->nbl; ibl++){
			if(!remove_bl[ibl]){
				nbln++;
				bb->next = pgv->ibl2bl[ibl].p;
				bb = bb->next;
				bb->next = nullptr;
			}
		}
	}
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->ibl2bl[ibl].p = nullptr;
	}
	pgv->nbl = nbln;
	bb = pgv->bl;
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->ibl2bl[ibl].p = bb;
		pgv->sg_rank_block[bb->ijk_sg[0]-1][bb->ijk_sg[1]-1][bb->ijk_sg[2]-1]=-ibl-1;
		bb = bb->next;
	}
}
sg_band_update();
if(pp->ierr==0){
	for(auto i=0; i<pgv->ijkm_sg[0]; i++){
		for(auto j=0; j<pgv->ijkm_sg[1]; j++){
			delete [] blgl[i][j];
		}
		delete [] blgl[i];
	}
	delete [] blgl;
}
else pp->litError("sg_load_balance","deallocation error for blgl");
};

void sg::sg_load_balance_parmetis_fix(bool force){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bl> pbl(new bl);

bool do_bl_list_update;
vector<bool> remove_bl(pgv->max_bl);
int nbln, iedge, nedge, edgecut, wgtflag, nblh;
vector<int> nbl_pe(pp->nprocs), nbl_pe_true(pp->nprocs);
vector<int> vtxdist(pp->nprocs+1);
int myNodeCount, cch1, cch2, cch3, cch4, myCommCount;
vector<int> myNodeC(pp->nprocs),NodeC(pp->nprocs);
vector<int> myCommC(pp->nprocs),CommC(pp->nprocs);
int ***blgl, ***myblgl;
vector<int> partDest, partSrc, xadj, adjwgt, adjncy, vwgt, vwgt2, mypart;
vector<vector<int>> bijksg, mybijksg, nGnoder;
int ncon;
std::array<int, 5> options;
std::array<int, 3> ijk_sg_n, ijk, ijk_shift;
vector<double> tpwgts, ubvec;
block_t* bb;
if(pp->myrank == 0)
	std::cout<<pgv->clit<<"Starting parmetis fix load balancer"<<std::endl;
if(pgv->nbl==0) nblh = 2;
else nblh = pgv->nbl;
pp->ierr=MPI_Allgather(&pgv->nbl, 1,MPI_INTEGER,&nbl_pe_true[0],1,MPI_INTEGER,\
							  MPI_COMM_WORLD);
pp->ierr=MPI_Allgather(&nblh,1,MPI_INTEGER, &nbl_pe[0],1,MPI_INTEGER,\
							  MPI_COMM_WORLD);
vtxdist[0] = 1;
for(auto irk=0; irk<pp->nprocs; irk++){
	vtxdist[irk+1] = vtxdist[irk] + nbl_pe[irk];
}
myblgl = new int**[pgv->ijkm_sg[0]];
blgl = new int**[pgv->ijkm_sg[0]];
for(auto i=0; i<pgv->ijkm_sg[0]; i++){
	myblgl[i] = new int*[pgv->ijkm_sg[1]];
	blgl[i] = new int*[pgv->ijkm_sg[1]];
	for(auto j=0; j<pgv->ijkm_sg[1]; j++){
		myblgl[i][j] = new int[pgv->ijkm_sg[2]];
		blgl[i][j] = new int[pgv->ijkm_sg[2]];
	}
}
mybijksg.resize(3); bijksg.resize(3);
for(auto i=0; i<3; i++){
	mybijksg[i].resize(nblh);
	bijksg[i].resize(vtxdist[pp->nprocs]);
}
//1087 & 1088 checked
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	myblgl[pgv->b->ijk_sg[0]-1][pgv->b->ijk_sg[1]-1][pgv->b->ijk_sg[2]-1]=\
	vtxdist[pp->myrank]+ibl;
	for(auto i=0; i<mybijksg.size(); i++) mybijksg[i][ibl] = pgv->b->ijk_sg[i];
}
auto d=pgv->ijkm_sg[0]*pgv->ijkm_sg[1]*pgv->ijkm_sg[2];
pp->ierr=MPI_Allreduce(myblgl,blgl,d,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
if(pp->ierr==0){
	for(auto i=0; i<pgv->ijkm_sg[0]; i++){
		for(auto j=0; j<pgv->ijkm_sg[1]; j++){
			delete [] myblgl[i][j];
		}
		delete [] myblgl[i];
	}
	delete [] myblgl;
}
else pp->litError("sg_load_balance","deallocation error for myblgl");
xadj.resize(nblh+1);
vwgt.resize(nblh);
nGnoder.resize(26);
for(auto i=0; i<26; i++) nGnoder[i].resize(nblh);
iedge = 0;
myCommCount = 0;
myNodeCount = 0;
if(pgv->nbl==0){	
	std::fill(vwgt.begin(), vwgt.begin()+vwgt.size(),1);
	xadj[0] = 1;
	xadj[1] = 3;
	xadj[2] = 5;
	nedge = 4;
}
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		switch(pgv->loadBalanceband){
			case 1:
				vwgt[ibl]=pgv->b->NinA;
				break;
			case 2:
				vwgt[ibl]=pgv->b->NinT;
				break;
			case 3:
				vwgt[ibl]=pgv->b->NinN;
				break;
			case 4:
				vwgt[ibl]=pgv->b->NinW;
				break;
			case 5:
				vwgt[ibl]=pgv->b->NinX;
				break;
		}
		myNodeCount += static_cast<int>(vwgt[ibl]);
		xadj[ibl] = iedge+1;
		if(vtxdist[pp->myrank] == 1 && ibl == 0){
			for(auto i=0; i<pp->nprocs; i++){
				if(nbl_pe_true[i]==0) iedge+=2;
			}
		}
		for(auto ibf=0; ibf<26; ibf++){
			if(!(pgv->b->boundType[ibf]==pgv->bNeumann)){
				auto counter = 0;
				for(auto i=pgv->b->ijk0[0]+pgv->igs[ibf]-pgv->b->imino_;i<=pgv->b->ijk0[0]+\
					 pgv->ige[ibf]-pgv->b->imino_; i++){
					for(auto j=pgv->b->ijk0[1]+pgv->jgs[ibf]-pgv->b->jmino_; j<=pgv->b->ijk0[1]+\
						 pgv->jge[ibf]-pgv->b->jmino_; j++){
						for(auto k=pgv->b->ijk0[2]+pgv->kgs[ibf]-pgv->b->kmino_; k<=pgv->b->ijk0[2]+\
							 pgv->kge[ibf]-pgv->b->kmino_; k++){
							if(pgv->b->ijk2ic[i][j][k]>0) counter++;
						}
					}
				}
				nGnoder[ibf][ibl] = counter;
			}
			else nGnoder[ibf][ibl] = 0;
			if(pgv->b->nGnodeb[ibf]+nGnoder[ibf][ibl]>0){
				iedge++;
				if(pgv->sg_rank_block[pgv->b->ijk_sg[0]+pgv->dijk_nc[ibf][0]-1]\
											[pgv->b->ijk_sg[1]+pgv->dijk_nc[ibf][1]-1]\
											[pgv->b->ijk_sg[2]+pgv->dijk_nc[ibf][2]-1]>=0){
					myCommCount += static_cast<int>(pgv->b->nGnodeb[ibf]);
				}
			}
		}
	}
	xadj[pgv->nbl] = iedge+1;
	nedge = iedge;
}
adjwgt.resize(nedge);
adjncy.resize(nedge);
if(pgv->nbl==0){
	adjncy[0]=1;
	adjncy[1]=vtxdist[pp->myrank];
	adjncy[2]=1;
	adjncy[3]=vtxdist[pp->myrank]+1;
	std::fill_n(adjwgt.begin(),adjwgt.size(),1);
	iedge = 4;
}
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		if(vtxdist[pp->myrank]==1 && ibl==0){
			for(auto i=0; i<pp->nprocs; i++){
				if(nbl_pe_true[i]==0){
					iedge++;
					adjncy[iedge-1] = vtxdist[i];
					adjwgt[iedge-1] = 1;
					iedge++;
					adjncy[iedge-1] = vtxdist[i]+1;
					adjwgt[iedge-1] = 1;
				}
			}
		}
		for(auto ibf=0; ibf<26; ibf++){
			if(pgv->b->nGnodeb[ibf]+nGnoder[ibf][ibl]>0){
				iedge++;
				adjwgt[iedge-1]=pgv->b->nGnodeb[ibf]+nGnoder[ibf][ibl];
				adjncy[iedge-1]=blgl[pgv->b->ijk_sg[0]+pgv->dijk_nc[ibf][0]-1]\
										  [pgv->b->ijk_sg[1]+pgv->dijk_nc[ibf][1]-1]\
										  [pgv->b->ijk_sg[2]+pgv->dijk_nc[ibf][2]-1];
			}
		}
	}
}
ncon = 2;
wgtflag = 3;
tpwgts.resize(ncon*pp->nprocs);
ubvec.resize(ncon);
vwgt2.resize(nblh*ncon);
mypart.resize(nblh);
std::fill_n(tpwgts.begin(),tpwgts.size(),1.0/static_cast<double>(ncon*pp->nprocs));
std::fill_n(ubvec.begin(), ubvec.size(), 1.05);
std::fill_n(options.begin(), options.size(), 0); //  checked
for(auto ibl=0; ibl<nblh; ibl++){
	vwgt2[2*ibl]=1;
	vwgt2[2*ibl+1]=std::max(vwgt[ibl]/1000, 1);
}
for(auto i=0; i<adjwgt.size(); i++) adjwgt[i]=std::max(adjwgt[i]/1000,1);
auto numflag = 1;
//ParMETIS_V3_PartKway(&vtxdist[0],&xadj[0],&adjncy[0],&vwgt2[0],&adjwgt[0],\
	&wgtflag,&numflag,&ncon,&pp->nprocs,reinterpret_cast<float*>(&tpwgts[0]),\
	reinterpret_cast<float*>(&ubvec[0]),&options[0],&edgecut,&mypart[0],&(pp->MY_G_WORLD));
xadj.clear(); xadj.shrink_to_fit();
adjwgt.clear(); adjwgt.shrink_to_fit();
adjncy.clear(); adjncy.shrink_to_fit();
vwgt.clear(); vwgt.shrink_to_fit();
vwgt2.clear(); vwgt2.shrink_to_fit();
tpwgts.clear(); tpwgts.shrink_to_fit();
ubvec.clear(); ubvec.shrink_to_fit();
if(pgv->nbl==0) std::fill_n(mypart.begin(),mypart.size(),pp->myrank+1);
partDest.resize(vtxdist[pp->nprocs]);
partSrc.resize(vtxdist[pp->nprocs]);
for(auto i=0; i<pp->nprocs;i++){
	for(auto j=vtxdist[i]-1; j<vtxdist[i+1]-1; j++){
		partSrc[j] = i;
	}
}
for(auto i=0; i<vtxdist.size(); i++) vtxdist[i]--;
for(auto i=0; i<mypart.size(); i++) mypart[i]--;
pp->ierr=MPI_Allgatherv(mypart.data(),nblh,MPI_INTEGER,partDest.data(),nbl_pe.data(),\
								vtxdist.data(),MPI_INTEGER,MPI_COMM_WORLD);
for(auto i=0; i<vtxdist.size(); i++) vtxdist[i] *= 3;
for(auto i=0; i<nbl_pe.size(); i++) nbl_pe[i] *= 3;
pp->ierr=MPI_Allgatherv(mybijksg.data(),3*nblh,MPI_INTEGER,bijksg.data(),\
								nbl_pe.data(),vtxdist.data(), MPI_INTEGER, MPI_COMM_WORLD);
for(auto i=0; i<vtxdist.size(); i++) vtxdist[i] = vtxdist[i]/3+1;
for(auto i=0; i<nbl_pe.size(); i++) nbl_pe[i] /= 3;
for(auto i=0; i<vtxdist.size(); i++) vtxdist[i]++;
pp->ierr=MPI_Reduce(&myNodeCount,&cch1,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD);
pp->ierr=MPI_Reduce(&myNodeCount,&cch2,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD);
pp->ierr=MPI_Reduce(&myNodeCount,&cch3,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD);
if(pp->myrank==0 && pgv->verbose){
	auto d1 = static_cast<double>(cch1)/static_cast<double>(pp->nprocs);
	auto d  = static_cast<double>(cch2)/d1;
	std::cout<<pgv->clit<<"pre-load balance node-imbalance"<<d<<std::endl;
	std::cout<<pgv->clit<<"pre-load balance # of communication nodes: ";
	std::cout<<cch3/2<<std::endl;
	d = static_cast<double>(0.5*cch3)/static_cast<double>(cch1);
	std::cout<<pgv->clit<<"pre-load balance ratio communication/compit nodes";
	std::cout<<d<<std::endl;
}
std::fill_n(myNodeC.begin(),myNodeC.size(),0);
std::fill_n(myCommC.begin(),myCommC.size(),0);
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	switch(pgv->loadBalanceband){
		case 1:
			myNodeC[mypart[ibl]] += static_cast<int>(pgv->b->NinA);
			break;
		case 2:
			myNodeC[mypart[ibl]] += static_cast<int>(pgv->b->NinT);
			break;
		case 3:
			myNodeC[mypart[ibl]] += static_cast<int>(pgv->b->NinN);
			break;
		case 4:
			myNodeC[mypart[ibl]] += static_cast<int>(pgv->b->NinW);
			break;
		case 5:
			myNodeC[mypart[ibl]] += static_cast<int>(pgv->b->NinX);
			break;
	}
	for(auto ibf=0; ibf<26; ibf++){
		if(pgv->b->nGnodeb[ibf]+nGnoder[ibf][ibl]>0){
			auto i = pgv->b->ijk_sg[0]+pgv->dijk_nc[ibf][0]-1;
			auto j = pgv->b->ijk_sg[1]+pgv->dijk_nc[ibf][1]-1;
			auto k = pgv->b->ijk_sg[2]+pgv->dijk_nc[ibf][2]-1;
			if(mypart[ibl] != partDest[blgl[i][j][k]-1]){
				myCommC[mypart[ibl]] += static_cast<int>(pgv->b->nGnodeb[ibf]);
			}
		}
	}
}
nGnoder.clear(); nGnoder.shrink_to_fit();
mypart.clear(); mypart.shrink_to_fit();
pp->ierr=MPI_Reduce(myNodeC.data(),&NodeC,pp->nprocs,MPI_INTEGER,MPI_SUM,0,\
						  MPI_COMM_WORLD);
pp->ierr=MPI_Reduce(myCommC.data(),&CommC,pp->nprocs,MPI_INTEGER,MPI_SUM,0,\
						  MPI_COMM_WORLD);
if(pp->myrank == 0 && pgv->verbose){
	auto m = *std::max_element(NodeC.begin(),NodeC.begin()+pp->nprocs);
	auto sum1 = 0;
	for(auto i=0; i<pp->nprocs; i++) sum1 += NodeC[i];
	auto d = static_cast<double>(sum1)/static_cast<double>(pp->nprocs);
	d = static_cast<double>(m)/d;
	std::cout<<pgv->clit<<"post-load balance node-imbalance: "<<d<<std::endl;
	auto sum2 = 0;
	for(auto i=0; i<pp->nprocs; i++) sum2 += CommC[i];
	std::cout<<pgv->clit<<"post-load balance # of communication nodes: ";
	std::cout<<sum2/2<<std::endl;
	sum1 = static_cast<double>(sum1);
	sum2 = static_cast<double>(sum2);
	std::cout<<pgv->clit<<"post-load balance ratio communication/compute nodes";
	std::cout<<0.5*sum2/sum1<<std::endl;
}
std::fill(remove_bl.begin(),remove_bl.end(),false);
do_bl_list_update = false;
for(auto i=0; i<vtxdist[pp->nprocs]-1; i++){
	if(partSrc[i] == pp->myrank){
		if(partDest[i] != pp->myrank){
			auto ibl = i-vtxdist[pp->myrank]+1;
			bb = pgv->ibl2bl[ibl].p;
			pp->parallel_block_send(bb, partDest[i]);
			pgv->sg_rank_block[bb->ijk_sg[0]-1][bb->ijk_sg[1]-1][bb->ijk_sg[2]-1]=\
			partDest[i];
			pbl->bl_free(bb);
			remove_bl[ibl] = true;
			do_bl_list_update = true;
		}
	}
	else if(partDest[i] == pp->myrank){
		pgv->nbl++;
		auto ibl = pgv->nbl-1;
#ifdef DEBUG_MODE
		if(pgv->ibl2bl[ibl].p != nullptr){
			std::cout<<pp->cmy<<"trying to allocate block already exists"<<std::endl;
			pp->parallel_kill(1);
		}
#endif
		if(pp->ierr==0) pgv->ibl2bl[ibl].p = new block_t;
		else pp->litError("sg_load_balance","allocation error fot ibl2bl[ibl].p");
		bb = pgv->ibl2bl[ibl].p;
		pp->parallel_block_recv(bb, partSrc[i]);
		do_bl_list_update = true;
		pgv->sg_rank_block[bb->ijk_sg[0]-1][bb->ijk_sg[1]-1][bb->ijk_sg[2]-1]=-ibl-1;
#ifdef DEBUG_MODE
		if(!pgv->sg_active[bb->ijk_sg[0]-1][bb->ijk_sg[1]-1][bb->ijk_sg[2]-1]){
			std::cout<<pgv->clit<<pp->myrank;
			std::cout<<"error:will receive block marked  as not active"<<std::endl;
			pp->parallel_kill(23);
		}
#endif
	}
	else{
		if(bijksg[0][i]+bijksg[1][i]+bijksg[2][i]>0)
			pgv->sg_rank_block[bijksg[0][i]-1][bijksg[1][i]-1][bijksg[2][i]-1]=
			partDest[i];
	}
}
if(do_bl_list_update){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		if(remove_bl[ibl]){
			if(pp->ierr==0) pgv->ibl2bl[ibl].p=nullptr;
			else pp->litError("sg_load_balance","deallocation error for ibl2bl[ibl].p");
		}
	}
	auto ibl = 0;
	while(remove_bl[ibl] && ibl+1<pgv->nbl){
		ibl++;
	}
	if(remove_bl[ibl]){
		pgv->bl = nullptr;
		std::cout<<pp->cmy<<"What am I doing here"<<std::endl;
		nbln = 0;
	}
	else{
		pgv->bl = pgv->ibl2bl[ibl].p;
		bb = pgv->bl;
		bb->next = nullptr;
		nbln = 1;
		auto ibls = ibl+1;
		for(auto ibl=ibls; ibl<pgv->nbl; ibl++){
			if(!remove_bl[ibl]){
				nbln++;
				bb->next = pgv->ibl2bl[ibl].p;
				bb = bb->next;
				bb->next = nullptr;
			}
		}
	}
	for(auto ibl=0; ibl<pgv->nbl; ibl++) pgv->ibl2bl[ibl].p = nullptr;
	pgv->nbl = nbln;
	bb = pgv->bl;
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->ibl2bl[ibl].p = bb;
		pgv->sg_rank_block[bb->ijk_sg[0]-1][bb->ijk_sg[1]-1][bb->ijk_sg[2]-1]=-ibl-1;
		bb = bb->next;
	}
}
sg_band_update();
if(pp->ierr==0){
	for(auto i=0; i<pgv->ijkm_sg[0]; i++){
		for(auto j=0; j<pgv->ijkm_sg[1]; j++){
			delete [] blgl[i][j];
		}
		delete [] blgl[i];
	}
	delete [] blgl;
}
else pp->litError("sg_load_balance","deallocation error for blgl");
};

void sg::sg_load_balance_lit(){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
bool do_bl_list_update, done, found;
vector<bool> remove_bl(pgv->max_bl);
int nin_pe_source, nex_pe, nin_l, nin_pe_target,nbln,ihh,pe_source,pe_target;
double nin_pe_av;
vector<int> nbl_pe(pp->nprocs), nbl_pe5(pp->nprocs), displ(pp->nprocs);
vector<int> displ5(pp->nprocs), nin_pe(pp->nprocs);
vector<vector<int>> nin_bl_pe, nin_bl, ex_script;
vector<bool> bl_remove;
block_t* bb;
if(pp->myrank==0 && pgv->verbose){
	std::cout<<pgv->clit<<"Starting lit load balancer"<<std::endl;
}
pp->parallel_gather(pgv->nbl, nbl_pe);
if(pp->myrank==0){
	auto sum=0;
	for(auto i=0; i<nbl_pe.size(); i++) sum += nbl_pe[i];
	pgv->nbl_gl = sum;
	displ[0]=0;
	for(auto irk=1; irk<pp->nprocs; irk++) displ[irk]=displ[irk-1]+nbl_pe[irk-1];
}
nin_bl.resize(5);
for(auto i=0; i<5; i++) nin_bl[i].resize(pgv->nbl);
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	bb = pgv->ibl2bl[ibl].p;
	nin_bl[0][ibl] = bb->NinN;
	nin_bl[1][ibl] = ibl+1;
	nin_bl[2][ibl] = bb->ijk_sg[0];
	nin_bl[3][ibl] = bb->ijk_sg[1];
	nin_bl[4][ibl] = bb->ijk_sg[2];
}
nin_l = 0;
for(auto i=0; i<pgv->nbl; i++) nin_l += nin_bl[0][i];
sort_pick(nin_bl);
if(pp->myrank==0){
	nin_bl_pe.resize(5);
	for(auto i=0; i<5; i++) nin_bl_pe[i].resize(pgv->nbl_gl);
}
else{
	nin_bl_pe.resize(5);
	for(auto i=0; i<5; i++) nin_bl_pe[i].resize(1);
}
for(auto i=0; i<nbl_pe5.size(); i++){
	nbl_pe5[i] = 5*nbl_pe5[i];
	displ5[i]  = 5*displ[i];
}
pp->ierr=MPI_Gatherv(nin_bl.data(),5*pgv->nbl,MPI_INTEGER,nin_bl_pe.data(),\
							nbl_pe5.data(),displ5.data(),MPI_INTEGER,0,MPI_COMM_WORLD);
pp->parallel_gather(nin_l,nin_pe);
if(pp->myrank==0){
	bl_remove.resize(pgv->nbl_gl);
	ex_script.resize(6);
	for(auto i=0; i<6; i++) ex_script[i].resize(pgv->nbl_gl);
	std::fill(bl_remove.begin(),bl_remove.end(), false);
	auto sum = 0;
	for(auto i=0; i<nin_pe.size(); i++) sum += nin_pe[i];
	nin_pe_av = static_cast<double>(sum)/static_cast<double>(pp->nprocs);
	done = false;
	nex_pe = 0;
	while(!done){
		auto res = std::max_element(nin_pe.begin(),nin_pe.end());
		ihh = std::distance(nin_pe.begin(),res);
		pe_source = ihh;
		nin_pe_source = nin_pe[pe_source];
		if(nbl_pe[pe_source]==0) done=true;
		else{
			auto temp = static_cast<double>(nin_pe_source)/static_cast<double>\
																		  (nin_pe_av);
			if(temp>pp->max_load_imbalance){
				auto res = std::min_element(nin_pe.begin(), nin_pe.end());
				ihh = std::distance(nin_pe.begin(), res);
				pe_target = ihh;
				auto ic = displ[pe_source];
				found = false;
				while(!found && nin_bl_pe[0][ic]>0 && ic+1<=displ[pe_source]+\
																		  nbl_pe[pe_source]){
					if(!bl_remove[ic]){
						nin_pe_target = nin_pe[pe_target]+nin_bl_pe[0][ic];
						if(nin_pe_target < nin_pe[pe_source]) found=true;
						if(!found) ic++;
					}
					else ic++;
					if(ic+1>displ[pe_source]+nbl_pe[pe_source]) break;
				}
				if(!found) done = true;
				else{
					nex_pe++;
					ex_script[0][nex_pe-1] = pe_source-1;
					ex_script[1][nex_pe-1] = pe_target-1;
					ex_script[2][nex_pe-1] = nin_bl_pe[1][ic];
					for(auto i=3; i<6; i++){
						ex_script[i][nex_pe-1] = nin_bl_pe[i-1][ic];
					}
					nin_pe[pe_source] -= nin_bl_pe[0][ic];
					nin_pe[pe_target] = nin_pe_target;
					bl_remove[ic] = true;
				}
			}
			else{
				done = true;
			}
		}
	}
	bl_remove.clear(); bl_remove.shrink_to_fit();
}
pp->ierr=MPI_Bcast(&nex_pe, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
if(pp->myrank==0){
	ex_script.resize(6);
	for(auto i=0; i<6; i++) ex_script[i].resize(nex_pe);
}
pp->ierr=MPI_Bcast(ex_script.data(),6*nex_pe,MPI_INTEGER,0,MPI_COMM_WORLD);
executeBlockExchange(nex_pe, ex_script);
// 1666 checked
};

void sg::executeBlockExchange(const int nex,const vector<vector<int>> &exScript){
bool do_bl_list_update = false;
int nbln, ibls;
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bl> pbl(new bl);
vector<bool> remove_bl(pgv->max_bl);
block_t* bb;
std::fill(remove_bl.begin(), remove_bl.end(), false);
for(auto iex=0; iex<nex; iex++){
	if(exScript[0][iex] == pp->myrank){
		bb = pgv->ibl2bl[exScript[2][iex]-1].p;
		pp->parallel_block_send(bb,exScript[0][iex]);
		pbl->bl_free(bb);
		pgv->sg_rank_block[exScript[3][iex]-1][exScript[4][iex]-1][exScript[5][iex]-1]\
		= exScript[1][iex];
		remove_bl[exScript[2][iex]] = true;
		do_bl_list_update = true;
	}
	else if(exScript[1][iex] == pp->myrank){
		pgv->nbl++;
		auto ibl = pgv->nbl-1;
#ifdef DEBUG_MODE
		if(pgv->ibl2bl[ibl].p != nullptr)
			pp->litError("executeBlockExchange","trying to allocate block that exist");
#endif
		if(pp->ierr==0) pgv->ibl2bl[ibl].p = new block_t;
		else pp->litError("executeBlockExchange","allocation error for ibl2bl...");
		bb = pgv->ibl2bl[ibl].p;
		pp->parallel_block_recv(bb,exScript[0][iex]);
		do_bl_list_update = true;
		pgv->sg_rank_block[exScript[3][iex]-1][exScript[4][iex]-1][exScript[5][iex]-1]=-ibl-1;
#ifdef DEBUG_MODE
	  if(!pgv->sg_active[exScript[3][iex]-1][exScript[4][iex]-1][exScript[5][iex]-1])
			pp->litError("executeBlockExchange","will receive block marked not active");
#endif
	}
	else
		pgv->sg_rank_block[exScript[3][iex]-1][exScript[4][iex]-1][exScript[5][iex]-1]\
		= exScript[1][iex];
}
if(do_bl_list_update){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		if(remove_bl[ibl]){
			if(pp->ierr==0) pgv->ibl2bl[ibl].p = nullptr;
			else pp->litError("executeBlockExchange","deallocation error for ibl2bl[ibl].p");
		}
	}
	auto ibl=0;
	while(remove_bl[ibl] && ibl+1<pgv->nbl) ibl++;
	if(remove_bl[ibl]){
		pgv->bl = nullptr;
		nbln = 0;
	}
	else{
		pgv->bl = pgv->ibl2bl[ibl].p;
		bb = pgv->bl;
		bb->next = nullptr;
		nbln = 1;
		ibls = ibl+1;
		for(ibl=ibls; ibl<pgv->nbl; ibl++){
			if(!remove_bl[ibl]){
				nbln++;
				bb->next = pgv->ibl2bl[ibl].p;
				bb = bb->next;
				bb->next = nullptr;
			}
		}
	}
	for(auto ibl=0; ibl<pgv->nbl; ibl++) pgv->ibl2bl[ibl].p = nullptr;
	pgv->nbl = nbln;
	if(pgv->nbl>0){
		bb = pgv->bl;
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->ibl2bl[ibl].p = bb;
			pgv->sg_rank_block[bb->ijk_sg[0]-1][bb->ijk_sg[1]-1][bb->ijk_sg[2]-1] = \
			-ibl-1;
			bb = bb->next;
		}
	}
}
sg_band_update();
};
std::vector<int> sg::xyz_to_ijk_sg(const std::array<double, 3> &xyz){
std::unique_ptr<global_variable> pgv(new global_variable);
std::vector<int> arr(3);
for(auto i=0; i<3; i++){
	arr[i] = static_cast<int>((xyz[i]-pgv->xyzs_sg[i])*pgv->rdxyz_sg[i])+1;
}
return arr;
};

void sg::sort_pick(vector<vector<int>>& arr){
std::array<int, 5> a;
auto n = arr[0].size();
for(auto j=0; j<n; j++){
	for(auto q=0; q<5; q++){
		a[q] = arr[q][j];
	}
	for(auto i=j-1; i>=0; i--){
		if(arr[0][i] >= a[0]) break;
		for(auto q=0; q<arr.size(); q++) arr[q][i+1] = arr[q][i];
	}
	for(auto q=0; q<arr.size(); q++) arr[q][1] = a[q]; //checked
}
};

void sg::monitor_load(double &imbalance){
int blockMin, blockMax;
int myN, myC, nMax, nSum, cSum;
double ratio;
std::unique_ptr<monitor> pmo(new monitor);
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
pp->parallel_min(pgv->nbl, blockMin);
pp->parallel_max(pgv->nbl, blockMax);
myN = 0; myC = 0;
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	auto sum = 0;
	for(auto i=0; i<pgv->ibl2bl[ibl].p->nGnodeb.size(); i++){
		sum += static_cast<int>(pgv->ibl2bl[ibl].p->nGnodeb[i]);
	}
	myN += static_cast<int>(pgv->ibl2bl[ibl].p->NinN);
	myC += sum; 
}
pp->ierr=MPI_Allreduce(&myN,&myC,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD);
pp->ierr=MPI_Allreduce(&myN,&nSum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
pp->ierr=MPI_Allreduce(&myN,&cSum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
imbalance = static_cast<double>(pp->nprocs)*static_cast<double>(nMax);
imbalance /= static_cast<double>(nSum);
ratio = 0.5*static_cast<double>(cSum)/static_cast<double>(nSum);
pmo->lit_monitor_select_file("lit_balance");
auto v3 = static_cast<double>(blockMin);
auto v4 = static_cast<double>(blockMax);
pmo->lit_monitor_set_single_values(imbalance, ratio, v3, v4);
};

