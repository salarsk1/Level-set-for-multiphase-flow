//Written by Salar Safarkhani

#include<array>
#include<algorithm>
#include<utility>
#include<vector>
#include<mpi.h>
#include"bound.h"
#include"timing.h"
#include"monitor.h"
#include"parallel.h"
#include"bl.h"
#include"init.h"
#include"gnodes.h"
void bound::bound_m_init(){
	std::unique_ptr<global_variable>pgv(new global_variable);
	std::unique_ptr<parallel>pp(new parallel);
	std::unique_ptr<timing>pti(new timing);
	std::unique_ptr<monitor>pmo(new monitor);

	if(pp->ierr==0){
		bound_nNodesPe.resize(pp->nprocs+1);
		ghost_nNodesPe.resize(pp->nprocs+1);
		ghost_nNodesNeumann.resize(pgv->max_bl+1);
	}
	else{
		pp->litError("bound_m_init","bound_nNodesPe, ghost_nNodesNeumann has problem");
	};
	for(auto i=0; i<pgv->max_bl+1; i++) ghost_nNodesNeumann[i] = 0;
	pti->lit_timing_create("bound_setup");
	pti->lit_timing_create("bound_g");
	pmo->lit_monitor_create_file_step("lit_bound",8);
	pmo->lit_monitor_set_header(1, "connect send max", 'i');
	pmo->lit_monitor_set_header(2, "connect recv max", 'i');
	pmo->lit_monitor_set_header(3, "connect local max", 'i');
	pmo->lit_monitor_set_header(4, "external max", 'i');
	pmo->lit_monitor_set_header(5, "connect send avg", 'i');
	pmo->lit_monitor_set_header(6, "connect recv avg", 'i');
	pmo->lit_monitor_set_header(7, "connect local avg", 'i');
	pmo->lit_monitor_set_header(8, "exteral avg", 'i');
};

void bound::prepareGhostNodes(){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<timing>pti(new timing);
std::unique_ptr<monitor>pmo(new monitor);
std::unique_ptr<gnodes>pgn(new gnodes);
std::unique_ptr<bl>pbl(new bl);
std::vector<int> ijk_sg_n(3, 0);
std::vector<int> ijk_shift(3, 0);
std::vector<int> ijk(3, 0);
block_t* bb;
Gnode_t Gnh;
int ibl, ibf, irank, nrecvs, nsends, is, ibln, ibfn, sizeComm, n;
int proc, ics, is2, ncells, irecv, ic, isn, index, itag, ipe, i, j, k;
int nMessagesSend, selfj, selfSize, procTarget, ii, jj, kk, nn;
std::array<int, 4> myEdgeCount, edgeCountSum, edgeCountMax;
bool selfComm;
std::vector<std::vector<int>> sendBufInt, recvBufInt;
std::vector<std::vector<int>> ncells_recv(26);
std::vector<std::vector<int>> ncells_send(26);

for(auto i=0; i<26; i++){
	ncells_recv[i].resize(pgv->nbl);
	ncells_send[i].resize(pgv->nbl);
}
std::vector<int> recv_requests(26*pgv->nbl), send_requests(26*pgv->nbl);
std::vector<int> irecv2iscript(26*pgv->nbl);
size_t sss = pp->nprocs;
int *iNodePe;
iNodePe = new int[pp->nprocs];
for(auto i=0; i<pp->nprocs; i++) iNodePe[i] = 0;
std::vector<int> iNodeNeumann(pgv->nbl, 0);
int bf2sg[3][26][pgv->nbl];
int bf2sh[3][26][pgv->nbl];
vector<int> ncount(pgv->nbl);
MPI_Request requestsInt[pp->nprocs];
double ttt;
pti->lit_timing_start("bound_setup");
pbl->bl_clear_all_ghosts(false);
if(bf2pe.size() != 0){
	bf2pe.clear();
	bf2pe.shrink_to_fit();
}
if(pp->ierr==0){
	bf2pe.resize(26);
	for(auto i=0;i<26;i++) bf2pe[i].resize(pgv->nbl);
}
else{
	pp->litError("prepareGhostNodes","Error");
};
for(auto i=0; i<26; i++){
	for(auto j=0; j<pgv->nbl; j++){
		bf2pe[i][j] = -1;
	}
}
for(auto i=0; i<3; i++){
	for(auto j=0; j<26; j++){
		for(auto k=0; k<pgv->nbl; k++){
			bf2sg[i][j][k] = -1;
			bf2sh[i][j][k] = 0;
		}
	}
}
if(bound_Gnodes.size()!=0){
	if(pp->ierr==0){
		bound_Gnodes.clear();
		bound_Gnodes.shrink_to_fit();
	}
	else 
		pp->litError("prepareGhostNodes","Deallocation bound_Gnodes err");
};
if(ghost_Neumann.size()!=0){
	if(pp->ierr==0){
		ghost_Neumann.clear();
		ghost_Neumann.shrink_to_fit();
	}
	else 
		pp->litError("prepareGhostNodes","Deallocation ghost_Neumann err");
};

if(ghost_index.size()!=0){
	if(pp->ierr==0){
		ghost_index.clear();
		ghost_index.shrink_to_fit();
	}
	else 
		pp->litError("prepareGhostNodes","Deallocation ghost_index err");
};
std::fill_n(bound_nNodesPe.begin(), bound_nNodesPe.size(), 0);
std::fill_n(ghost_nNodesPe.begin(), ghost_nNodesPe.size(), 0);
for(auto i=0; i<pgv->nbl+1; i++) ghost_nNodesNeumann[i] = 0;

n = 0;
if(pgv->cylindrical && pgv->ijkm_gl[2]>1){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ibf=0;ibf<26; ibf++){
			vector<int> v_temp(3);
			for(auto i=0;i<3;i++){
				v_temp[i] = pgv->dijk_nc[ibf][i];
			}
				pgv->shift_ijk_sg_cyl(pgv->b->ijk_sg, v_temp, ijk_sg_n,ijk_shift);
				bool check1=false;
				for(auto i=0;i<ijk_sg_n.size();i++){
					if(ijk_sg_n[i]>=1 && ijk_sg_n[i]<=pgv->ijkm_sg[i]) check1=true;
				}
				if(check1){
					auto d1 = ijk_sg_n[0];
					auto d2 = ijk_sg_n[1];
					auto d3 = ijk_sg_n[2];
					if(pgv->sg_active[d1-1][d2-1][d3-1]){
						if(pgv->sg_rank_block[d1-1][d2-1][d3-1]){
							ipe = pp->myrank+1;
						}
						else{
							ipe = pgv->sg_rank_block[d1-1][d2-1][d3-1] + 1;
						};
						bf2pe[ibf][ibl] = ipe;
						for(auto i=0;i<3;i++) bf2sg[i][ibf][ibl] = ijk_sg_n[i];
						for(auto i=0;i<3;i++) bf2sh[i][ibf][ibl] = ijk_shift[i];
						if(pgv->dijk_nc[ibf][1]==-1 && pgv->b->ijk_sg[1]==1){
							if(pgv->ijkm_sg[2]==1){
								auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
								auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
								auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
								auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
								auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]+\
 					         pgv->dijk_nc[ibf][2]*(pgv->nghost-pgv->ijkm_gl[2]/2)-pgv->b->kmino_;
								auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]+\
 					         pgv->dijk_nc[ibf][2]*(pgv->nghost-pgv->ijkm_gl[2]/2)-pgv->b->kmino_;
								for(auto i=d11; i<=d12; i++){
									for(auto j=d21; j<=d22; j++){
										for(auto k=d31; k<=d32; k++){
											if(pgv->b->ijk2ic[i][j][k] > 0) n+=1;
										}
									}
								}
							}
							else{
								auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
								auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
								auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
								auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
								auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]+\
								pgv->dijk_nc[ibf][2]*(pgv->nghost-pgv->ijkm_bl[2])-pgv->b->kmino_;
								auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]+\
								pgv->dijk_nc[ibf][2]*(pgv->nghost-pgv->ijkm_bl[2])-pgv->b->kmino_;
								for(auto i=d11; i<=d12; i++){
									for(auto j=d21; j<=d22; j++){
										for(auto k=d31; k<=d32; k++){
											if(pgv->b->ijk2ic[i][j][k] > 0) n+=1;
										}
									}
								}
							};
						}
						else{
							auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
							auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
							auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
							auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
							auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]-pgv->b->kmino_;
							auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]-pgv->b->kmino_;
							for(auto i=d11; i<=d12; i++){
								for(auto j=d21; j<=d22; j++){
									for(auto k=d31; k<=d32; k++){
										if(pgv->b->ijk2ic[i][j][d31] > 0) n+=1;
									}
								}
							}
						};
						bound_nNodesPe[ipe] += n;
						pgv->b->boundType[ibf] = pgv->bConnect;
						pgv->b->nGnodeb[ibf] = n;
					}
					else{
						pgv->b->nGnodeb[ibf] = 0;
					};
				}
				else{
					n = 0;
					if(pgv->dijk_nc[ibf][1]=-1 && pgv->b->ijk_sg[1]==1){
						if(pgv->ijkm_sg[2]==1){
							auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
							auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
							auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
							auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
							auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]+\
						   pgv->dijk_nc[ibf][2]*(pgv->nghost-pgv->ijkm_gl[2]/2)-pgv->b->kmino_;
							auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]+\
						   pgv->dijk_nc[ibf][2]*(pgv->nghost-pgv->ijkm_gl[2]/2)-pgv->b->kmino_;
							for(auto i=d11; i<=d12; i++){
								for(auto j=d21; j<=d22; j++){
									for(auto k=d31; k<=d32; k++){
										if(pgv->b->ijk2ic[i][j][k] > 0) n+=1;
									}
								}
							}
						}
						else{
							auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
							auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
							auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
							auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
							auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]+\
				         pgv->dijk_nc[ibf][2]*(pgv->nghost-pgv->ijkm_bl[2])-pgv->b->kmino_;
							auto d32= pgv->b->ijk0[2]+pgv->kbe[ibf]+\
				         pgv->dijk_nc[ibf][2]*(pgv->nghost-pgv->ijkm_bl[2])-pgv->b->kmino_;
							for(auto i=d11; i<=d12; i++){
								for(auto j=d21; j<=d22; j++){
									for(auto k=d31; k<=d32; k++){
										if(pgv->b->ijk2ic[i][j][k] > 0) n+=1;
									}
								}
							}
						};
						if(n>0){
							pgv->b->boundType[ibf] = pgv->bCenterLine;
							ghost_nNodesNeumann[ibl+1] += n;
						};
						pgv->b->nGnodeb[ibf] = 0;
					}
					else{
						if(pgv->ijkm_gl[2] == 1){
							if(pgv->dijk_nc[ibf][2] != 0){
								auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
								auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
								auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
								auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
								for(auto i=d11; i<=d12; i++){
									for(auto j=d21; j<=d22; j++){
										if(pgv->b->ijk2ic[i][j][1-pgv->b->kmino_]>0) n+=1;
									}
								}
								n*=pgv->nghost;
							}
							else{
								auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
								auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
								auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
								auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
								for(auto i=d11; i<=d12; i++){
									for(auto j=d21; j<=d22; j++){
										if(pgv->b->ijk2ic[i][j][1-pgv->b->kmino_]>0) n+=1;
									}
								}
							};
						}
						else{
							auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
							auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
							auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
							auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
							auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]-pgv->b->kmino_;
							auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]-pgv->b->kmino_;
							for(auto i=d11; i<=d12; i++){
								for(auto j=d21; j<=d22; j++){
									for(auto k=d31; k<=d32; k++){
										if(pgv->b->ijk2ic[i][j][k]>0) n+=1;
									}
								}
							}
						};
						if(n>0){
							pgv->b->boundType[ibf] = pgv->bNeumann;
							ghost_nNodesNeumann[ibl+1] += n;
						};
						pgv->b->nGnodeb[ibf] = 0;
					};
				};
		};
	};
}
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ibf=0; ibf<26; ibf++){
			n = 0;
			std::vector<int> v_temp;
			for(auto t1=0; t1<3; t1++) v_temp.push_back(pgv->dijk_nc[ibf][t1]);
			pgv->shift_ijk_sg(pgv->b->ijk_sg, v_temp, ijk_sg_n, ijk_shift);
			bool check=true;
			for(auto i=0; i<ijk_sg_n.size(); i++){
				if(ijk_sg_n[i] <1 || ijk_sg_n[i]>pgv->ijkm_sg[i]){
					check=false;
					break;
				}
			}
			if(check){
				n = 0;
				if(pgv->sg_active[ijk_sg_n[0]-1][ijk_sg_n[1]-1][ijk_sg_n[2]-1]){
					if(pgv->sg_rank_block[ijk_sg_n[0]-1][ijk_sg_n[1]-1][ijk_sg_n[2]-1]<0){
						ipe = pp->myrank + 1;
					}
					else{
						ipe = pgv->sg_rank_block[ijk_sg_n[0]-1][ijk_sg_n[1]-1][ijk_sg_n[2]-1]+1;
					};
					bf2pe[ibf][ibl] = ipe;
					for(auto i=0; i<3; i++){
						bf2sg[i][ibf][ibl] = ijk_sg_n[i];
						bf2sh[i][ibf][ibl] = ijk_shift[i];
					};
					auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
					auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
					auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
					auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
					auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]-pgv->b->kmino_;
					auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]-pgv->b->kmino_;
					for(auto i=d11; i<=d12; i++){
						for(auto j=d21; j<=d22; j++){
							for(auto k=d31; k<=d32; k++){
								if(pgv->b->ijk2ic[i][j][k]>0) n+=1;
							}
						}
					}
					bound_nNodesPe[ipe] +=  n;
					pgv->b->boundType[ibf]  = pgv->bConnect;
					pgv->b->nGnodeb[ibf]  = n;
				}
				else{
					pgv->b->nGnodeb[ibf] = 0;
				};
			}
			else{
				n = 0;
				if(pgv->ijkm_gl[2] == 1){
					if(pgv->dijk_nc[ibf][2] != 0){
						auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
						auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
						auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
						auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
						for(auto i=d11; i<=d12; i++){
							for(auto j=d21; j<=d22; j++){
								if(pgv->b->ijk2ic[i][j][1-pgv->b->kmino_]>0) n+=1;
							}
						}
						n *= pgv->nghost;
					}
					else{
						auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
						auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
						auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
						auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
						for(auto i=d11; i<=d12; i++){
							for(auto j=d21; j<=d22; j++){
								if(pgv->b->ijk2ic[i][j][1-pgv->b->kmino_]>0) n+=1;
							}
						}
					};
				}
				else{
					n = 0;
					auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
					auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
					auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
					auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
					auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]-pgv->b->kmino_;
					auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]-pgv->b->kmino_;
					for(auto i=d11; i<=d12; i++){
						for(auto j=d21; j<=d22; j++){
							for(auto k=d31; k<=d32; k++){
								if(pgv->b->ijk2ic[i][j][k]>0) n+=1;
							}
						}
					}
				};
				if(n>0){
					pgv->b->boundType[ibf] = pgv->bNeumann;
					ghost_nNodesNeumann[ibl+1] += n;
				};
				pgv->b->nGnodeb[ibf] = 0;
			};
		};
	};
};
pp->ierr = MPI_Alltoall(&bound_nNodesPe[1], 1, MPI_INTEGER, &ghost_nNodesPe[1], 1, \
								 MPI_INTEGER,MPI_COMM_WORLD);
ghost_nNodesRecv = 0;
ghost_nProcsRecv = 0;
ghost_nProcsSend = 0;
for(auto i=1; i<pp->nprocs+1; i++){
	ghost_nNodesRecv += ghost_nNodesPe[i];
	if(ghost_nNodesPe[i]>0) ghost_nProcsRecv += 1;
	if(bound_nNodesPe[i]>0) ghost_nProcsSend += 1;
};
for (ipe=1; ipe<pp->nprocs+1; ipe++){
	bound_nNodesPe[ipe] += bound_nNodesPe[ipe-1];
	ghost_nNodesPe[ipe] = ghost_nNodesPe[ipe] + ghost_nNodesPe[ipe-1];
};
for(auto ibl=1; ibl<pgv->nbl+1; ibl++){
	ghost_nNodesNeumann[ibl] += ghost_nNodesNeumann[ibl-1];
};
if(pp->ierr==0){
	bound_Gnodes.resize(3);
	ghost_Neumann.resize(5);
	for(auto i=0; i<3; i++) bound_Gnodes[i].resize(bound_nNodesPe[pp->nprocs]+1);
	for(auto i=0; i<5; i++) ghost_Neumann[i].resize(ghost_nNodesNeumann[pgv->nbl]);
}
else{
	pp->litError("prepareGhostNodes","Allocation error of bound_Gnodes.");
};
for(auto i=0; i<5; i++){
	for(auto j=0; j<ghost_nNodesNeumann[pgv->nbl]; j++){
		ghost_Neumann[i][j] = 0;
	}
}

for(auto tmp=0; tmp<3; tmp++) bound_Gnodes[tmp][0] = -1.0;

if(pp->ierr==0){
	sendBufInt.resize(6);
	for(auto i=0; i<6; i++) sendBufInt[i].resize(bound_nNodesPe[pp->nprocs]);
}
else{
	pp->litError("prepareGhostNodes","Allocation error of sendBufInt.");
};
if(pgv->cylindrical && pgv->ijkm_gl[2]>1){ 
	for(auto ibl=0; ibl<26; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ibf=0; ibf<26; ibf++){
			auto d11 = pgv->b->ijk0[0] + pgv->ibs[ibf]-pgv->b->imino_;
			auto d12 = pgv->b->ijk0[0] + pgv->ibe[ibf]-pgv->b->imino_;
			auto d21 = pgv->b->ijk0[1] + pgv->jbs[ibf]-pgv->b->jmino_;
			auto d22 = pgv->b->ijk0[1] + pgv->jbe[ibf]-pgv->b->jmino_;
			auto d31 = pgv->b->ijk0[2] + pgv->kbs[ibf]-pgv->b->kmino_;
			auto d32 = pgv->b->ijk0[2] + pgv->kbe[ibf]-pgv->b->kmino_;
			if(pgv->b->boundType[ibf] == pgv->bNeumann){
				if(pgv->ijkm_gl[2] == 1 && pgv->dijk_nc[ibf][2] !=0){
					for(auto i=d11; i<=d12; i++){
						for(auto j=d21; j<=d22; j++){
							for(auto k=d31; k<=d32; k++){
								if(pgv->b->ijk2ic[i][j][1-pgv->b->kmino_]>0){
									iNodeNeumann[ibl] = iNodeNeumann[ibl]+1;
									ghost_Neumann[0][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1] = ibf;
									ghost_Neumann[1][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1] = pgv->i2c[i][j][1-pgv->b->kmino_];
									ghost_Neumann[2][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1] = (i+pgv->b->imino_)+abs(pgv->dijk_nc[ibf]\
									[0])*(pgv->b->ijk_b2g[0][ibf]-2*(i+pgv->b->imino_)+1);
									ghost_Neumann[3][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1] = (j+pgv->b->jmino_)+abs(pgv->dijk_nc[ibf]\
									[1])*(pgv->b->ijk_b2g[1][ibf]-2*(j+pgv->b->jmino_)+1);
									ghost_Neumann[4][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1] = (k+pgv->b->kmino_)+abs(pgv->dijk_nc[ibf]\
									[2])*(pgv->b->ijk_b2g[2][ibf]-2*(k+pgv->b->kmino_)+1);
									for(auto q=0; q<3; q++){
										Gnh.ijk[q] = ghost_Neumann[q+2]\
										[ghost_nNodesNeumann[ibl]+iNodeNeumann[ibl]];
										Gnh.V[q] = 0.0;
										Gnh.dijk[q] = 0;
									}
									pgn->addGnode(pgv->b, Gnh);
									pgv->i2c[Gnh.ijk[0]-pgv->b->imino_][Gnh.ijk[1]-pgv->b->jmino_]\
											  [Gnh.ijk[2]-pgv->b->kmino_]=pgv->b->nG;
								}
							}
						}
					}
				}
				else{
					for(auto i=d11; i<=d12; i++){
						for(auto j=d21; j<=d22; j++){
							for(auto k=d31; k<=d32; k++){
								if(pgv->b->ijk2ic[i][j][k]>0){
									iNodeNeumann[ibl] = iNodeNeumann[ibl]+1;
									ghost_Neumann[0][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1] = ibf;
									ghost_Neumann[1][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1] = pgv->i2c[i][j][k];
									ghost_Neumann[2][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1] = (i+pgv->b->imino_)+abs(pgv->dijk_nc[ibf]\
									[0])*(pgv->b->ijk_b2g[0][ibf]-2*(i+pgv->b->imino_)+1);
									ghost_Neumann[3][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1] = (j+pgv->b->jmino_)+abs(pgv->dijk_nc[ibf]\
									[1])*(pgv->b->ijk_b2g[1][ibf]-2*(j+pgv->b->jmino_)+1);
									ghost_Neumann[4][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1] = (k+pgv->b->kmino_)+abs(pgv->dijk_nc[ibf]\
									[2])*(pgv->b->ijk_b2g[2][ibf]-2*(k+pgv->b->kmino_)+1);
									for(auto q=0; q<3; q++){
										Gnh.ijk[q] = ghost_Neumann[q+2]\
										[ghost_nNodesNeumann[ibl]+iNodeNeumann[ibl]];
										Gnh.V[q] = 0.0;
										Gnh.dijk[q] = 0;
									}
									pgn->addGnode(pgv->b, Gnh);
									pgv->i2c[Gnh.ijk[0]-pgv->b->imino_][Gnh.ijk[1]-pgv->b->jmino_]\
											  [Gnh.ijk[2]-pgv->b->kmino_]=pgv->b->nG;
								}
							}
						}
					}
				};
			}
			else if(pgv->b->boundType[ibf] == pgv->bCenterLine){
				auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
				auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
				auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
				auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
				auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]+pgv->dijk_nc[ibf][2]*\
							  (pgv->nghost-pgv->ijkm_gl[2]/2)-pgv->b->imino_;
				auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]+pgv->dijk_nc[ibf][2]*\
							  (pgv->nghost-pgv->ijkm_gl[2]/2)-pgv->b->imino_;
				for(auto i=d11;i<=d12;i++){
					for(auto j=d21; j<=d22; j++){
						for(auto k=d31;k<=d32; k++){
							if(pgv->b->ijk2ic[i][j][k]>0){
								iNodeNeumann[ibl] = iNodeNeumann[ibl]+1;
								ghost_Neumann[0][ghost_nNodesNeumann[ibl]+\
													  iNodeNeumann[ibl]-1] = ibf;
								ghost_Neumann[1][ghost_nNodesNeumann[ibl]+\
											  iNodeNeumann[ibl]-1] = pgv->i2c[i][j][k];
								if(pgv->ijkm_sg[2]==1 && pgv->dijk_nc[ibf][2]!=0){
									if((k+pgv->b->kmino_)-pgv->b->ijk0[2]<=pgv->ijkm_bl[2]/2){
										ghost_Neumann[2][ghost_nNodesNeumann[ibl]+\
															  iNodeNeumann[ibl]-1] = \
					 (i+pgv->b->imino_)+abs(pgv->dijk_nc[ibf][0])*pgv->b->ijk_b2g[0][ibf]-2*(i+pgv->b->imino_)+1;
										ghost_Neumann[3][ghost_nNodesNeumann[ibl]+\
															  iNodeNeumann[ibl]-1] = \
										1-(j+pgv->b->jmino_);
										ghost_Neumann[4][ghost_nNodesNeumann[ibl]+iNodeNeumann[ibl]-1]= \
										(k+pgv->b->kmino_)-pgv->ijkm_gl[2]/2;
									}
									else{
											ghost_Neumann[2][ghost_nNodesNeumann[ibl]+\
											iNodeNeumann[ibl]-1] = (i+pgv->b->imino_)+abs(pgv->dijk_nc[ibf][0])*\
											(pgv->b->ijk_b2g[0][ibf])-2*(i+pgv->b->imino_)+1;
											ghost_Neumann[3][ghost_nNodesNeumann[ibl]+\
											iNodeNeumann[ibl]-1] = 1-(j+pgv->b->jmino_);
											ghost_Neumann[4][ghost_nNodesNeumann[ibl]+\
											iNodeNeumann[ibl]-1] = (k+pgv->b->kmino_)+pgv->ijkm_gl[2]/2;
									}
								}
								else{
									if((k+pgv->b->kmino_)-pgv->b->ijk0[2]<=pgv->ijkm_bl[2]/2){
										ghost_Neumann[2][ghost_nNodesNeumann[ibl]+\
															  iNodeNeumann[ibl]-1]=
									 (i+pgv->b->imino_)+abs(pgv->dijk_nc[ibf][0])*pgv->b->ijk_b2g[0][ibf]-2*(i+pgv->b->imino_)+1;
										ghost_Neumann[3][ghost_nNodesNeumann[ibl]+\
															  iNodeNeumann[ibl]-1]=1-(j+pgv->b->jmino_);
										ghost_Neumann[4][ghost_nNodesNeumann[ibl]+\
														  iNodeNeumann[ibl]-1]=(k+pgv->b->kmino_)+pgv->ijkm_gl[2]/2;
									}
									else{
										ghost_Neumann[2][ghost_nNodesNeumann[ibl]+\
															  iNodeNeumann[ibl]-1]=
									 (i+pgv->b->imino_)+abs(pgv->dijk_nc[ibf][0])*pgv->b->ijk_b2g[0][ibf]-2*(i+pgv->b->imino_)+1;
										ghost_Neumann[3][ghost_nNodesNeumann[ibl]+\
															  iNodeNeumann[ibl]-1]=1-(j+pgv->b->jmino_);
										ghost_Neumann[4][ghost_nNodesNeumann[ibl]+\
														  iNodeNeumann[ibl]-1]=(k+pgv->b->kmino_)-pgv->ijkm_gl[2]/2;
									};
									nn = (i+pgv->b->imino_)*pgv->ijkm_gl[0]+\
										  (j+pgv->b->jmino_)*pgv->ijkm_gl[1]+\
										  (k+pgv->b->kmino_)*pgv->ijkm_gl[2];
									ii = k;
									for(auto q=0; q<3; q++){
										Gnh.ijk[q]=ghost_Neumann[q+2]\
													 [ghost_nNodesNeumann[ibl]+iNodeNeumann[ibl]-1];
										Gnh.V[q] = 0.0;
										pgn->addGnode(pgv->b, Gnh);
										pgv->i2c[Gnh.ijk[0]-pgv->b->imino_][Gnh.ijk[1]-pgv->b->jmino_]\
												  [Gnh.ijk[2]-pgv->b->kmino_]	= pgv->b->nG;		
									};
								};
							}
						}
					}
				}
			}
			else if(pgv->b->ijk_sg[1]==1 && pgv->dijk_nc[ibf][1]==-1){
				if(pgv->ijkm_sg[2]==1){
					if(bf2pe[ibf][ibl]>=0){
						ipe = bf2pe[ibf][ibl];
						auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
						auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
						auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
						auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
						auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]+pgv->dijk_nc[ibf][2]*\
									  (pgv->nghost-pgv->ijkm_gl[2]/2)-pgv->b->kmino_;
						auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]+pgv->dijk_nc[ibf][2]*\
									  (pgv->nghost-pgv->ijkm_gl[2]/2)-pgv->b->kmino_;
						for(auto i=d11; i<=d12; i++){
							for(auto j=d21; j<=d22; j++){
								for(auto k=d31; k<=d32; k++){
									if(pgv->b->ijk2ic[i][j][k]>0){
										iNodePe[ipe-1] +=1;
										bound_Gnodes[0][bound_nNodesPe[ipe-1]+\
															 iNodePe[ipe-1]]=ibl;
										bound_Gnodes[1][bound_nNodesPe[ipe-1]+\
															 iNodePe[ipe-1]]=pgv->i2c[i][j][k];
										bound_Gnodes[2][bound_nNodesPe[ipe-1]+\
															 iNodePe[ipe-1]]=ibf;
										if((k+pgv->b->kmino_)-pgv->b->ijk0[2]<=pgv->ijkm_bl[2]/2){
										 sendBufInt[0][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
										 (i+pgv->b->imino_)+bf2sh[0][ibf][ibl];
										 sendBufInt[1][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
										 1-(j+pgv->b->jmino_);
										 sendBufInt[2][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
										 (k+pgv->b->kmino_)+pgv->ijkm_gl[2]/2;
										}
										else{
										 sendBufInt[0][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
										 (i+pgv->b->imino_)+bf2sh[0][ibf][ibl];
										 sendBufInt[1][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
										 1-(j-pgv->b->jmino_);
										 sendBufInt[2][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
										 (k+pgv->b->kmino_)-pgv->ijkm_gl[2]/2;
										};
										for(auto q=3; q<6; q++){
										 sendBufInt[q][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
										 bf2sg[q-3][ibf][ibl];
										}
									}
								}
							}
						}
					}
				}
				else{
					if(bf2pe[ibf][ibl]>=0){
						ipe = bf2pe[ibf][ibl];
						auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
						auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
						auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
						auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
						auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]+pgv->dijk_nc[ibf][2]*\
									  (pgv->nghost-pgv->ijkm_bl[2])-pgv->b->kmino_;
						auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]+pgv->dijk_nc[ibf][2]*\
									  (pgv->nghost-pgv->ijkm_bl[2])-pgv->b->kmino_;
						for(auto i=d11;i<=d12; i++){
							for(auto j=d21;j<=d22; j++){
								for(auto k=d31; k<=d32; k++){
									if(pgv->b->ijk2ic[i][j][k]>0){
										iNodePe[ipe-1] += 1;
									 bound_Gnodes[0][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									 ibl;
									 bound_Gnodes[1][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									 pgv->i2c[i][j][k];
									 bound_Gnodes[2][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									 ibf;
										sendBufInt[0][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
										(i+pgv->b->imino_)+bf2sh[0][ibf][ibl];
										sendBufInt[1][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
										1-(j+pgv->b->jmino_);
										sendBufInt[2][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
										(k+pgv->b->kmino_)+bf2sh[2][ibf][ibl];
										for(auto q=3; q<6; q++){
									 	sendBufInt[q][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
										bf2sg[q-3][ibf][ibl];
										}
									}
								}
							}
						}
					}
				}
			}
			else{
				if(bf2pe[ibf][ibl]>=0){
					ipe = bf2pe[ibf][ibl];
					auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
					auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
					auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
					auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
					auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]-pgv->b->kmino_;
					auto d32 = pgv->b->ijk0[2]+pgv->kbs[ibf]-pgv->b->kmino_;
					for(auto i=d11; i<=d12; i++){
						for(auto j=d21; j<=d22; j++){
							for(auto k=d31; k<=d32; k++){
								if(pgv->i2c[i][j][k]>0){
									iNodePe[ipe-1] += 1;
									bound_Gnodes[0][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									ibl;
									bound_Gnodes[1][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									pgv->i2c[i][j][k];
									bound_Gnodes[2][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									ibf;
									sendBufInt[0][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									(i+pgv->b->imino_)+bf2sh[0][ibf][ibl];
									sendBufInt[1][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									(j+pgv->b->jmino_)+bf2sh[1][ibf][ibl];
									sendBufInt[2][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									(k+pgv->b->kmino_)+bf2sh[2][ibf][ibl];
									for(auto q=3; q<6; q++){
									sendBufInt[q][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									bf2sg[q-3][ibf][ibl];
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
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ibf=0; ibf<26; ibf++){

			if(pgv->b->boundType[ibf]==pgv->bNeumann){
				if(pgv->ijkm_gl[2]==1 && pgv->dijk_nc[ibf][2] != 0){

					auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
					auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
					auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
					auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
					auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]-pgv->b->kmino_;
					auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]-pgv->b->kmino_;
					for(auto i=d11; i<=d12; i++){
						for(auto j=d21; j<=d22; j++){
							for(auto k=d31; k<=d32; k++){
								if(pgv->b->ijk2ic[i][j][1-pgv->b->kmino_]>0){

									iNodeNeumann[ibl] += 1;
					
									ghost_Neumann[0][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1]=ibf+1;
					
									ghost_Neumann[1][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1]=pgv->i2c[i][j][1-pgv->b->kmino_];
					
									ghost_Neumann[2][ghost_nNodesNeumann[ibl]+iNodeNeumann[ibl]-1]=\
									i+pgv->b->imino_+abs(pgv->dijk_nc[ibf][0])*(pgv->b->ijk_b2g[0][ibf]\
									-2*(i+pgv->b->imino_)+1);

									ghost_Neumann[3][ghost_nNodesNeumann[ibl]+iNodeNeumann[ibl]-1]=\
									j+pgv->b->jmino_+abs(pgv->dijk_nc[ibf][1])*(pgv->b->ijk_b2g[1][ibf]\
									-2*(j+pgv->b->jmino_)+1);

									ghost_Neumann[4][ghost_nNodesNeumann[ibl]+iNodeNeumann[ibl]-1]=\
									k+pgv->b->kmino_+abs(pgv->dijk_nc[ibf][2])*(pgv->b->ijk_b2g[2][ibf]\
									-2*(k+pgv->b->kmino_)+1);

									for(int q=0; q<3; q++){
										Gnh.ijk[q] = ghost_Neumann[q+2][ghost_nNodesNeumann[ibl]+\
														 iNodeNeumann[ibl]-1];
										Gnh.V[q] = 0.0;
										Gnh.dijk[q] = 0;
									}
									pgn->addGnode(pgv->b, Gnh);

									pgv->i2c[Gnh.ijk[0]-pgv->b->imino_][Gnh.ijk[1]-pgv->b->jmino_]\
									[Gnh.ijk[2]-pgv->b->kmino_] = pgv->b->nG;

								}
							}
						}
					}
				}
				else{
					auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
					auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
					auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
					auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
					auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]-pgv->b->kmino_;
					auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]-pgv->b->kmino_;
					for(auto i=d11; i<=d12; i++){
						for(auto j=d21; j<=d22; j++){
							for(auto k=d31; k<=d32; k++){
								if(pgv->b->ijk2ic[i][j][k]>0){
									iNodeNeumann[ibl] += 1;
									ghost_Neumann[0][ghost_nNodesNeumann[ibl]+iNodeNeumann[ibl]-1]=ibf+1;

									ghost_Neumann[1][ghost_nNodesNeumann[ibl]+\
									iNodeNeumann[ibl]-1]=pgv->i2c[i][j][k];

									ghost_Neumann[2][ghost_nNodesNeumann[ibl]+iNodeNeumann[ibl]-1]=\
									(i+pgv->b->imino_)+abs(pgv->dijk_nc[ibf][0])*\
									(pgv->b->ijk_b2g[0][ibf]-2*(i+pgv->b->imino_)+1);

									ghost_Neumann[3][ghost_nNodesNeumann[ibl]+iNodeNeumann[ibl]-1]=\
									(j+pgv->b->jmino_)+abs(pgv->dijk_nc[ibf][1])*\
									(pgv->b->ijk_b2g[1][ibf]-2*(j+pgv->b->jmino_)+1);

									ghost_Neumann[4][ghost_nNodesNeumann[ibl]+iNodeNeumann[ibl]-1]=\
									(k+pgv->b->kmino_)+abs(pgv->dijk_nc[ibf][2])*\
									(pgv->b->ijk_b2g[2][ibf]-2*(k+pgv->b->kmino_)+1);

									for(auto q=0; q<3; q++){
										Gnh.ijk[q] = ghost_Neumann[q+2][ghost_nNodesNeumann[ibl]+\
														 iNodeNeumann[ibl]-1];
										Gnh.V[q] = 0.0;
										Gnh.dijk[q] = 0;
									}

									pgn->addGnode(pgv->b, Gnh);
									pgv->i2c[Gnh.ijk[0]-pgv->b->imino_][Gnh.ijk[1]-pgv->b->jmino_]\
											  [Gnh.ijk[2]-pgv->b->kmino_] = pgv->b->nG;
								}
							}
						}
					}				
				};
			}
			else{
				if(bf2pe[ibf][ibl]>=0){
					ipe = bf2pe[ibf][ibl];
					auto d11 = pgv->b->ijk0[0]+pgv->ibs[ibf]-pgv->b->imino_;
					auto d12 = pgv->b->ijk0[0]+pgv->ibe[ibf]-pgv->b->imino_;
					auto d21 = pgv->b->ijk0[1]+pgv->jbs[ibf]-pgv->b->jmino_;
					auto d22 = pgv->b->ijk0[1]+pgv->jbe[ibf]-pgv->b->jmino_;
					auto d31 = pgv->b->ijk0[2]+pgv->kbs[ibf]-pgv->b->kmino_;
					auto d32 = pgv->b->ijk0[2]+pgv->kbe[ibf]-pgv->b->kmino_;
					for(auto i=d11; i<=d12; i++){
						for(auto j=d21; j<=d22; j++){
							for(auto k=d31; k<=d32; k++){
								if(pgv->i2c[i][j][k]>0){
									iNodePe[ipe-1] += 1;

									bound_Gnodes[0][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									ibl+1;
	
									bound_Gnodes[1][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									pgv->i2c[i][j][k];

									bound_Gnodes[2][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]]=\
									ibf+1;

									sendBufInt[0][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
									(i+pgv->b->imino_)+bf2sh[0][ibf][ibl];

									sendBufInt[1][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
									(j+pgv->b->jmino_)+bf2sh[1][ibf][ibl];

									sendBufInt[2][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
									(k+pgv->b->kmino_)+bf2sh[2][ibf][ibl];

									for(auto q=3; q<6; q++){
										sendBufInt[q][bound_nNodesPe[ipe-1]+iNodePe[ipe-1]-1]=\
										bf2sg[q-3][ibf][ibl];
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
nMessagesSend = 0;
selfComm = false;
for(auto ipe=1; ipe<=pp->nprocs; ipe++){
	if(bound_nNodesPe[ipe]>bound_nNodesPe[ipe-1]){
		sizeComm = (bound_nNodesPe[ipe]-bound_nNodesPe[ipe-1]);
		j = bound_nNodesPe[ipe-1]+1;
		if(ipe==pp->myrank+1){
			selfComm = true;
			selfj    = j;
			selfSize = sizeComm;
		}
		else{
			nMessagesSend += 1;
			procTarget     = ipe-1;
			MPI_Issend(&sendBufInt[0][j-1], 6*sizeComm, MPI_INTEGER, ipe-1, \
						  ipe, MPI_COMM_WORLD, &requestsInt[nMessagesSend-1]);
		}
	}
}
ii=1;
ic = 0;
if(pp->ierr == 0){
	ghost_index.resize(6);
	for(auto i=0; i<6; i++) ghost_index[i].resize(ghost_nNodesRecv);
}
else{
	pp->litError("prepareGhostNodes","Allocation error of ghost_index");
}
if(selfComm){
	j = selfj;
	i = 0;
	for(auto q=ghost_nNodesPe[pp->myrank]; q<ghost_nNodesPe[pp->myrank]+selfSize; q++){
		ghost_index[0][q] = sendBufInt[0][q-(ghost_nNodesPe[pp->myrank]+1-selfj)];
	}
	for(auto q=ghost_nNodesPe[pp->myrank]; q<ghost_nNodesPe[pp->myrank]+selfSize; q++){
		ghost_index[1][q] = sendBufInt[1][q-(ghost_nNodesPe[pp->myrank]+1-selfj)];
	}
	for(auto q=ghost_nNodesPe[pp->myrank]; q<ghost_nNodesPe[pp->myrank]+selfSize; q++){
		ghost_index[2][q] = sendBufInt[2][q-(ghost_nNodesPe[pp->myrank]+1-selfj)];
	}
	for(auto q=ghost_nNodesPe[pp->myrank]; q<ghost_nNodesPe[pp->myrank]+selfSize; q++){
		ghost_index[3][q] = sendBufInt[3][q-(ghost_nNodesPe[pp->myrank]+1-selfj)];
	}
	for(auto q=ghost_nNodesPe[pp->myrank]; q<ghost_nNodesPe[pp->myrank]+selfSize; q++){
		ghost_index[4][q] = sendBufInt[4][q-(ghost_nNodesPe[pp->myrank]+1-selfj)];
	}
	for(auto q=ghost_nNodesPe[pp->myrank]; q<ghost_nNodesPe[pp->myrank]+selfSize; q++){
		ghost_index[5][q] = sendBufInt[5][q-(ghost_nNodesPe[pp->myrank]+1-selfj)];
	}
	ic  = selfSize;
	ii += 1;
}
for(nn=ii; nn<=ghost_nProcsRecv; nn++){
	pp->ierr = MPI_Probe(MPI_ANY_SOURCE, pp->myrank+1, MPI_COMM_WORLD, &(pp->status));
	proc = pp->status.MPI_SOURCE;
	MPI_Get_count(&(pp->status), MPI_INTEGER, &sizeComm);
	sizeComm /= 6;
	MPI_Recv(&ghost_index[0][ghost_nNodesPe[proc]], 6*sizeComm, MPI_INTEGER,\
				proc, pp->myrank+1, MPI_COMM_WORLD, &(pp->status));
	ic += sizeComm;
}
sizeComm = nMessagesSend;

for(nn=0; nn<nMessagesSend; nn++){
  MPI_Waitany(sizeComm, requestsInt, &proc, &(pp->status));
}
if(pp->ierr==0){
	sendBufInt.clear();
	sendBufInt.shrink_to_fit();
}
else{
 pp->litError("prepareGhostNodes", "deallocation error at 1 of sendBuf.");
}
std::fill_n(ncount.begin(), pgv->nbl, 0);
for(ic=0; ic<ghost_nNodesRecv; ic++){
	ibl=-pgv->sg_rank_block[ghost_index[3][ic]-1][ghost_index[4][ic]-1][ghost_index[5][ic]-1];
	ncount[ibl-1] += 1;
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	bb = pgv->ibl2bl[ibl].p;
	if(bb->nG+ncount[ibl]>bb->nGmax){
		pgn->enlargenGnodes(bb, ncount[ibl]);
	}
}
for(ic=0; ic<ghost_nNodesRecv; ic++){
	ibl = -(pgv->sg_rank_block[ghost_index[3][ic]-1][ghost_index[4][ic]-1][ghost_index[5][ic]-1]);
	bb = pgv->ibl2bl[ibl-1].p;
	for(auto i=0; i<3; i++){
		Gnh.ijk[i] = ghost_index[i][ic];
		Gnh.V[i] = 0.0;
		Gnh.dijk[i] = 0;
		Gnh.G = 0.0;
	}
	if(ghost_index[0][ic]<bb->imino_ || ghost_index[0][ic]>bb->imaxo_ ||\
		ghost_index[1][ic]<bb->jmino_ || ghost_index[1][ic]>bb->jmaxo_ ||\
		ghost_index[2][ic]<bb->kmino_ || ghost_index[2][ic]>bb->kmaxo_){
		std::cout<<"error: ibl, ic, nG, ijk0, i, j, k="<<ibl<<'\t'<<ic<<'\t'<<\
						'\t'<<bb->ijk0[0]<<'\t'<<bb->ijk0[1]<<'\t'<<bb->ijk0[2]<<\
						'\t'<<ghost_index[0][ic]<<'\t'<<ghost_index[1][ic]<<'\t'<<\
						ghost_index[2][ic]<<std::endl;
		std::cout<<"nneumann, nrecv"<<(ghost_nNodesNeumann[ibl+1]-ghost_nNodesNeumann\
												 [ibl])<<'\t'<<ghost_nNodesRecv<<std::endl;
		pp->litError("prepareGhostNodes","consistency error");
	}
	pgn->addGnode(bb, Gnh);
	bb->ijk2ic[ghost_index[0][ic]-bb->imino_][ghost_index[1][ic]-bb->jmino_][ghost_index[2][ic]-bb->kmino_]=bb->nG;
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	pgv->b->NinZ = pgv->b->NinX+(ghost_nNodesNeumann[ibl+1]-ghost_nNodesNeumann[ibl])+ncount[ibl];
}
myEdgeCount[0]=bound_nNodesPe[pp->nprocs]-(bound_nNodesPe[pp->myrank+1]-bound_nNodesPe[pp->myrank]);
myEdgeCount[1]=ghost_nNodesRecv-(bound_nNodesPe[pp->myrank+1]-bound_nNodesPe[pp->myrank]);
myEdgeCount[2] = (bound_nNodesPe[pp->myrank+1]-bound_nNodesPe[pp->myrank]);
myEdgeCount[3] = ghost_nNodesNeumann[pgv->nbl];
pp->ierr = MPI_Reduce(&myEdgeCount[0], &edgeCountSum[0], 4, MPI_INTEGER,MPI_SUM,0, MPI_COMM_WORLD);
pp->ierr = MPI_Reduce(&myEdgeCount[0], &edgeCountMax[0], 4, MPI_INTEGER,MPI_MAX,0, MPI_COMM_WORLD);
if(pp->myrank ==0){
	pmo->lit_monitor_select_file("lit_bound");
	auto v1 = static_cast<double>(edgeCountMax[0]);
	auto v2= static_cast<double>(edgeCountMax[1]);
	auto v3= static_cast<double>(edgeCountMax[2]);
	auto v4= static_cast<double>(edgeCountMax[3]);
	auto v5= static_cast<double>(edgeCountSum[0])/static_cast<double>(pp->nprocs);
	auto v6= static_cast<double>(edgeCountSum[1])/static_cast<double>(pp->nprocs);
	auto v7= static_cast<double>(edgeCountSum[2])/static_cast<double>(pp->nprocs);
	auto v8= static_cast<double>(edgeCountSum[3])/static_cast<double>(std::max(pgv->nbl,1));
	pmo->lit_monitor_set_single_values(v1, v2, v3, v4, v5, v6, v7, v8);
}

pti->lit_timing_stop("bound_setup");
};
void bound::updateGhostNodes(){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<timing>pti(new timing);
std::unique_ptr<parallel>pp(new parallel);
block_t* bb = nullptr;
int ibl, ibf, id, ic, sizeComm, n, proc, ii, nMessagesSend, ipe;
int j, selfj, procTarget, i, nn, k, iic, jjc, kkc, nnc, selfSize;
std::vector<MPI_Request> requests(pp->nprocs,0); 
std::vector<MPI_Request> requestsInt(pp->nprocs,0);
bool selfComm;
std::vector<double> xyz(3,0.0);
std::vector<double> sendBuf, recvBuf;
pti->lit_timing_start("bound_g");
if(pp->ierr==0)
	sendBuf.resize(bound_nNodesPe[pp->nprocs]);
else
	pp->litError("updateGhostNodes","Allocation error of sendBuf");
for(auto ic=0; ic<bound_nNodesPe[pp->nprocs]; ic++){
	sendBuf[ic] = pgv->ibl2bl[bound_Gnodes[0][ic+1]-1].p->Gnodes[bound_Gnodes[1][ic+1]-1].G;
}
nMessagesSend = 0;
selfComm = false;
for(auto ipe=1; ipe<=pp->nprocs; ipe++){
	if(bound_nNodesPe[ipe]>bound_nNodesPe[ipe-1]){
		sizeComm = (bound_nNodesPe[ipe]-bound_nNodesPe[ipe-1]);
		j = bound_nNodesPe[ipe-1]+1;
		if(ipe==pp->myrank+1){
			selfComm = true;	
			selfj = j;
			selfSize = sizeComm;
		}
		else{
			nMessagesSend += 1;
			procTarget = ipe-1;
			MPI_Issend(&sendBuf[j-1], sizeComm, MPI_DOUBLE, ipe-1,\
						  pp->nprocs+ipe, MPI_COMM_WORLD, &requests[nMessagesSend-1]);
		}
	}
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	if(ghost_nNodesNeumann[ibl+1]>ghost_nNodesNeumann[ibl]){
		pgv->setBlockPointers(ibl);
		for(auto i=ghost_nNodesNeumann[ibl]; i<ghost_nNodesNeumann[ibl+1]; i++){


			pgv->Gn[pgv->i2c[ghost_Neumann[2][i]-pgv->b->imino_][ghost_Neumann[3][i]-pgv->b->jmino_]\
			[ghost_Neumann[4][i]-pgv->b->kmino_]-1].G=pgv->Gn[ghost_Neumann[1][i]-1].G;

			pgv->Gn[pgv->i2c[ghost_Neumann[2][i]-pgv->b->imino_][ghost_Neumann[3][i]-pgv->b->jmino_]\
						  [ghost_Neumann[4][i]-pgv->b->kmino_]-1].V=pgv->Gn[ghost_Neumann[1][i]-1].V;

			pgv->Gn[pgv->i2c[ghost_Neumann[2][i]-pgv->b->imino_][ghost_Neumann[3][i]-pgv->b->jmino_]\
						  [ghost_Neumann[4][i]-pgv->b->kmino_]-1].dijk=\
						  pgv->Gn[ghost_Neumann[1][i]-1].dijk;

		}
	}
}
ii = 1;
if(selfComm){
	ic = ghost_nNodesPe[pp->myrank];
	j  = selfj;
	for(j=selfj; j<selfj+selfSize; j++){
		ic += 1;
	ibl=-pgv->sg_rank_block[ghost_index[3][ic-1]-1][ghost_index[4][ic-1]-1][ghost_index[5][ic-1]-1];
		bb = pgv->ibl2bl[ibl-1].p;
		bb->Gnodes[bb->ijk2ic[ghost_index[0][ic-1]-bb->imino_][ghost_index[1][ic-1]-bb->jmino_]\
									 [ghost_index[2][ic-1]-bb->kmino_]-1].G = sendBuf[j-1];
	}
	ii++;
}
for(nn=ii-1; nn<ghost_nProcsRecv; nn++){
	MPI_Probe(MPI_ANY_SOURCE,pp->nprocs+pp->myrank+1,MPI_COMM_WORLD,&(pp->status));
	proc = pp->status.MPI_SOURCE;
	MPI_Get_count(&(pp->status), MPI_DOUBLE, &sizeComm);
	if(pp->ierr==0)
		recvBuf.resize(sizeComm);
	else
		pp->litError("updateGhostNodes","allocation error of recvBuf");
	MPI_Recv(recvBuf.data(), sizeComm, MPI_DOUBLE, proc, pp->nprocs+pp->myrank+1,\
				MPI_COMM_WORLD, &(pp->status));
	ic = ghost_nNodesPe[proc]-1;
	for(auto j=0; j<sizeComm; j++){
		ic += 1;
		ibl = -pgv->sg_rank_block[ghost_index[3][ic]-1][ghost_index[4][ic]-1][ghost_index[5][ic]-1];
		bb = pgv->ibl2bl[ibl-1].p;
		bb->Gnodes[bb->ijk2ic[ghost_index[0][ic]-bb->imino_][ghost_index[1][ic]-bb->jmino_]\
									[ghost_index[2][ic]-bb->kmino_]-1].G = recvBuf[j];
	}
	if(pp->ierr==0){
		recvBuf.clear();
		recvBuf.shrink_to_fit();
	}
	else
		pp->litError("updateGhostNodes","deallocation error of recvBuf");
}
sizeComm = nMessagesSend;
for(auto nn=0; nn<nMessagesSend; nn++){
	MPI_Waitany(sizeComm, &requests[nn], &proc, &(pp->status));
}
if(pp->ierr==0){
	sendBuf.clear();
	sendBuf.shrink_to_fit();
}
else
	pp->litError("updateGhostNodes","deallocation at 1 of sendBuf");
pti->lit_timing_stop("bound_g");
};

void bound::updateGhostR1(double* r1, int stype){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
block_t* bb = nullptr;
int ibl, ibf, id, ic, sizeComm, n, nMessagesSend, ipe, j, selfj;
int selfSize, procTarget, s_type, i, nn, proc;
MPI_Request requests[pp->nprocs], requestsInt[pp->nprocs];
bool selfComm;
double* sendBuf;
double* recvBuf;

s_type = S_SCALAR;
if(stype!=0) s_type = stype;
if(pp->ierr==0) sendBuf= new double[bound_nNodesPe[pp->nprocs]];
else pp->litError("updateGhostNodes","Allocation error of sendBuf.");
for(auto ic=0; ic<bound_nNodesPe[pp->nprocs]; ic++){
	ibl = bound_Gnodes[0][ic+1];
	sendBuf[ic] = r1[pgv->ibl2buf[ibl-1]+bound_Gnodes[1][ic+1]-1];
}
nMessagesSend = 0;
selfComm = false;
for(auto ipe=0;ipe<pp->nprocs; ipe++){
	if(bound_nNodesPe[ipe+1]>bound_nNodesPe[ipe]){
		sizeComm = bound_nNodesPe[ipe+1]-bound_nNodesPe[ipe];
		j = bound_nNodesPe[ipe] +1;
		if(ipe == pp->myrank){
			selfComm = true;
			selfj = j;
			selfSize = sizeComm;
		}
		else{
			nMessagesSend += 1;
			procTarget = ipe;
			MPI_Issend(&sendBuf[j-1], sizeComm, MPI_DOUBLE, ipe,\
						  3*pp->nprocs+ipe+1, MPI_COMM_WORLD,\
						  &requests[nMessagesSend-1]);
		}
	}
}
if(s_type==S_SCALAR){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		if(ghost_nNodesNeumann[ibl+1]>ghost_nNodesNeumann[ibl]){
			pgv->setBlockPointers(ibl);
			for(auto i=ghost_nNodesNeumann[ibl];i<ghost_nNodesNeumann[ibl+1]; i++){
				r1[pgv->ioff+pgv->i2c[ghost_Neumann[2][i]-pgv->b->imino_][ghost_Neumann\
				   [3][i]-pgv->b->jmino_][ghost_Neumann[4][i]-pgv->b->kmino_]-1] =\
			   r1[pgv->ioff+ghost_Neumann[1][i]];
			}
		}
	}
}
else if(s_type = S_GRADIENT){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		if(ghost_nNodesNeumann[ibl+1]>ghost_nNodesNeumann[ibl]){
			pgv->setBlockPointers(ibl);
			for(auto i=ghost_nNodesNeumann[ibl];i<ghost_nNodesNeumann\
															  [ibl+1]; i++){
				r1[pgv->ioff+pgv->i2c[ghost_Neumann[2][i]-pgv->b->imino_][ghost_Neumann\
				   [3][i]-pgv->b->jmino_][ghost_Neumann[4][i]-pgv->b->kmino_]-1] =\
			  -r1[pgv->ioff+ghost_Neumann[1][i]];
			}
		}
	}
}
else if(s_type == S_ZERO){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		if(ghost_nNodesNeumann[ibl+1]>ghost_nNodesNeumann[ibl]){
			pgv->setBlockPointers(ibl);
			for(auto i=ghost_nNodesNeumann[ibl];i<ghost_nNodesNeumann\
															  [ibl+1]; i++){
				r1[pgv->ioff+pgv->i2c[ghost_Neumann[2][i]-pgv->b->imino_][ghost_Neumann\
				   [3][i]-pgv->b->jmino_][ghost_Neumann[4][i]-pgv->b->kmino_]-1] = 0.0;
			}
		}
	}
};
auto ii=1;
if(selfComm){
	j = selfj;
	ic = ghost_nNodesPe[pp->myrank] - 1;
	for(auto j=selfj-1; j<selfj-selfSize-1; j++){
		ic += 1;
		ibl=-pgv->sg_rank_block[ghost_index[3][ic]-1][ghost_index[4][ic]-1][ghost_index[5][ic]-1]-1;
		bb = pgv->ibl2bl[ibl].p;
		r1[pgv->ibl2buf[ibl]+bb->ijk2ic[ghost_index[0][ic]-bb->imino_][ghost_index[1][ic]-bb->jmino_]\
												 [ghost_index[2][ic]-bb->kmino_]-1] = sendBuf[j];
	}
	ii += 1;
}
for(auto nn=ii-1; nn<ghost_nProcsRecv; nn++){
	MPI_Probe(MPI_ANY_SOURCE, 3*pp->nprocs+pp->myrank+1, MPI_COMM_WORLD,\
				 &(pp->status));
	proc = pp->status.MPI_SOURCE;
	MPI_Get_count(&(pp->status), MPI_DOUBLE, &sizeComm);
	if(pp->ierr==0){
		recvBuf = new double[sizeComm];
	}
	else{
		pp->litError("updateGhostR1", "allocation error of recvBuf");
	}
	MPI_Recv(recvBuf, sizeComm, MPI_DOUBLE, proc, 3*pp->nprocs+pp->myrank+1,\
				MPI_COMM_WORLD, &(pp->status));
	ic = ghost_nNodesPe[proc]-1;
	for(auto j=0; j<sizeComm; j++){
		ic += 1;
		ibl = -pgv->sg_rank_block[ghost_index[3][ic]-1][ghost_index[4][ic]-1][ghost_index[5][ic]-1]-1;
		bb = pgv->ibl2bl[ibl].p;
		r1[pgv->ibl2buf[ibl]+bb->ijk2ic[ghost_index[0][ic]-bb->imino_][ghost_index[1][ic]-bb->jmino_]\
												 [ghost_index[2][ic]-bb->kmino_]-1] = recvBuf[j];
	}
	if(pp->ierr==0) delete [] recvBuf;
	else pp->parallel_die("deallocation error of rwecvBuf in updateGhostNodes");
}
sizeComm = nMessagesSend;
for(auto nn=0; nn<nMessagesSend; nn++){
	MPI_Waitany(sizeComm, requests, &proc, &(pp->status));
}
if(pp->ierr==0) delete [] sendBuf;
else pp->litError("updateGhostR1", "deallocatipn error at sendBuf.");
};

void bound::updateGhostI1(int* i1, int stype){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
block_t* bb;
int sizeComm, n, proc, nMessagesSend, selfj, selfSize, procTarget;
int s_type, ii;
std::vector<MPI_Request> requestsInt(pp->nprocs, 0);
bool selfComm;
int *sendBufInt;
int *recvBufInt;
s_type = S_SCALAR;
if(s_type != 0) s_type = stype;
if(pp->ierr==0)
	sendBufInt = new int[bound_nNodesPe[pp->nprocs]];
else
	pp->litError("updateGhostNodes","Allocation error of sendBufInt.");
for(auto ic=0; ic<bound_nNodesPe[pp->nprocs];ic++){
	auto ibl = bound_Gnodes[0][ic+1]-1;
	sendBufInt[ic] = i1[pgv->ibl2buf[ibl]+bound_Gnodes[1][ic+1]-1];
}
nMessagesSend = 0;
selfComm = false;
for(auto ipe=0; ipe<pp->nprocs; ipe++){
	if(bound_nNodesPe[ipe+1]>bound_nNodesPe[ipe]){
		sizeComm = bound_nNodesPe[ipe+1]-bound_nNodesPe[ipe];
		auto j = bound_nNodesPe[ipe] + 1;
		if(ipe == pp->myrank){
			selfComm = false;
			selfj = j;
			selfSize = sizeComm;
		}
		else{
			nMessagesSend += 1;
			procTarget = ipe;
			MPI_Issend(&sendBufInt[j],sizeComm,MPI_INTEGER,ipe,5*pp->nprocs+ipe+1,\
						  MPI_COMM_WORLD, &requestsInt[nMessagesSend-1]);
		}
	}	
}
if(s_type != S_KEEP){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		if(ghost_nNodesNeumann[ibl+1]>ghost_nNodesNeumann[ibl]){
			pgv->b = pgv->ibl2bl[ibl].p;
			for(auto i=ghost_nNodesNeumann[ibl];i<ghost_nNodesNeumann[ibl+1];i++){
				i1[pgv->ibl2buf[ibl]+pgv->b->ijk2ic[ghost_Neumann[2][i]-pgv->b->imino_]\
						 [ghost_Neumann[3][i]-pgv->b->jmino_][ghost_Neumann[4][i]-pgv->b->kmino_]-1]=\
				i1[pgv->ibl2buf[ibl]+ghost_Neumann[1][i]-1];
			}
		}
	}
}
ii = 0;
if(selfComm){
	auto ic = ghost_nNodesPe[pp->myrank]-1;
	for(auto j=selfj-1; j<selfj+selfSize-1; j++){
		ic +=1;
	auto ibl=-pgv->sg_rank_block[ghost_index[3][ic]-1][ghost_index[4][ic]-1][ghost_index[5][ic]-1]-1;
		bb = pgv->ibl2bl[ibl].p;
		i1[pgv->ibl2buf[ibl]+bb->ijk2ic[ghost_index[0][ic]-bb->imino_][ghost_index[1][ic]-bb->jmino_]\
												 [ghost_index[2][ic]-bb->kmino_]-1] = sendBufInt[j];
	}
ii += 1;
}
for(auto nn=ii; nn<ghost_nProcsRecv; nn++){
	MPI_Probe(MPI_ANY_SOURCE,5*pp->nprocs+pp->myrank+1,MPI_COMM_WORLD,&(pp->status));
	proc = pp->status.MPI_SOURCE;
	pp->ierr = MPI_Get_count(&(pp->status), MPI_INTEGER, &sizeComm);
	if(pp->ierr==0)
		recvBufInt= new int[sizeComm];
	else
		pp->litError("updateGhostI1","allocation error of recvBufInt");
	pp->ierr = MPI_Recv(recvBufInt, sizeComm, MPI_INTEGER, proc,\
							   5*pp->nprocs+pp->myrank+1, MPI_COMM_WORLD, &(pp->status));
	auto ic = ghost_nNodesPe[proc]-1;
	for(auto j=0; j<sizeComm; j++){
		ic += 1;
	auto ibl=-pgv->sg_rank_block[ghost_index[3][ic]-1][ghost_index[4][ic]-1][ghost_index[5][ic]-1]-1;
		bb = pgv->ibl2bl[ibl].p;
		i1[pgv->ibl2buf[ibl]+bb->ijk2ic[ghost_index[0][ic]-bb->imino_][ghost_index[1][ic]-bb->jmino_]\
												 [ghost_index[2][ic]-bb->kmino_]-1] = recvBufInt[j];
	}
	if(pp->ierr==0)
		delete [] recvBufInt;
	else
		pp->parallel_die("deallocation error of recvBufInt in updateGhostNodes");
}
sizeComm = nMessagesSend;
for(auto nn=0; nn<nMessagesSend;nn++){
	pp->ierr = MPI_Waitany(sizeComm, requestsInt.data(), &proc, &(pp->status));
}
if(pp->ierr==0) delete [] sendBufInt;
else pp->litError("updateGhostI1","deallocation error at sendBuf.");
};
void bound::updateGhostR2(double **r2, size_t r2_S1, size_t r2_S2, int stype){
block_t* bb=nullptr;
int sizeComm, proc, nMessagesSend, selfj, selfSize, procTarget, ndirs, s_type;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
MPI_Request requests[pp->nprocs];
MPI_Request requestsInt[pp->nprocs];
std::vector<int> dirs(3,0);
bool selfComm;
double **sendBuf;
double **recvBuf;
s_type= S_SCALAR;
if(stype!=0) s_type = stype;
sendBuf = new double*[3];
if(pp->ierr==0)
	for(auto i=0; i<3; i++) sendBuf[i] = new double[bound_nNodesPe[pp->nprocs]];
else
	pp->litError("updateGhostNodes","Allocation error of sendBuf");
for(auto ic=0; ic<bound_nNodesPe[pp->nprocs]; ic++){
	auto ibl = bound_Gnodes[0][ic+1]-1;
	for(auto i=0; i<3; i++) sendBuf[i][ic] = r2[i][pgv->ibl2buf[ibl]+bound_Gnodes[1][ic+1]-1];
}
nMessagesSend = 0;
selfComm = false;
for(auto ipe=0; ipe<pp->nprocs;ipe++){
	if(bound_nNodesPe[ipe+1]>bound_nNodesPe[ipe]){
		sizeComm = bound_nNodesPe[ipe+1] - bound_nNodesPe[ipe];
		auto j = bound_nNodesPe[ipe]+1;
		if(ipe==pp->myrank){
			selfComm = true;
			selfj = j;
			selfSize = sizeComm;
		}
		else{
			nMessagesSend += 1;
			procTarget = ipe;
			pp->ierr = MPI_Issend(&sendBuf[0][j-1], 3*sizeComm, MPI_DOUBLE, ipe,\
						  6*pp->nprocs+ipe+1, MPI_COMM_WORLD, &requests[nMessagesSend-1]);
		}
	}
}
if(pgv->cylindrical){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		if(ghost_nNodesNeumann[ibl+1]>ghost_nNodesNeumann[ibl]){
			pgv->setBlockPointers(ibl);
			for(auto i=ghost_nNodesNeumann[ibl];i<ghost_nNodesNeumann[ibl+1]; i++){
				auto ibf = ghost_Neumann[0][i];
				auto icg = pgv->ioff+pgv->i2c[ghost_Neumann[2][i]-pgv->b->imino_]\
							  [ghost_Neumann[3][i]-pgv->b->kmino_][ghost_Neumann[4][i]-1];
				auto ic = ghost_Neumann[1][i];
				auto icb = pgv->ioff+ic;
				for(auto q=0; q<3; q++) r2[q][icg-1] = r2[q][icb-1];
				if(pgv->b->boundType[ibf-1]==pgv->bCenterLine && s_type==S_GRADIENT)
					for(auto q=1; q<3; q++) r2[q][icg-1] = -r2[q][icb-1];
				else{
					get_bound_normals(ibf, ndirs, dirs);
					for(auto idir=0;idir<ndirs;idir++){
						if(s_type==S_GRADIENT)
							r2[dirs[idir]-1][icg-1] = -r2[dirs[idir]-1][icb-1];
						else if(s_type==S_DM){
							if(ibf==2 || ibf ==4 || ibf==6) r2[dirs[idir]-1][icg-1] = 0.0;
							else{
								auto ih=pgv->Gn[ic-1].ijk[0]-pgv->dijk_nc[0][ibf-1];
								auto jh=pgv->Gn[ic-1].ijk[1]-pgv->dijk_nc[1][ibf-1];
								auto kh=pgv->Gn[ic-1].ijk[2]-pgv->dijk_nc[2][ibf-1];
								if(pgv->i2c[ih-pgv->b->imino_][jh-pgv->b->jmino_][kh-pgv->b->kmino_]>0){
									r2[dirs[idir]-1][icg-1] = \
								  -r2[dirs[idir]-1][pgv->ioff+pgv->i2c[ih-pgv->b->imino_]\
									  [jh-pgv->b->jmino_][kh-pgv->b->kmino_]-1];
								}
								else{
									r2[dirs[idir]-1][icg-1] = 0.0;
								}
							}
						}
						else if(s_type = S_DP){
							if(ibf==1 || ibf==3 || ibf==5) r2[dirs[idir]-1][icg-1];
							else{
								auto ih=pgv->Gn[ic-1].ijk[0]-pgv->dijk_nc[0][ibf-1];
								auto jh=pgv->Gn[ic-1].ijk[1]-pgv->dijk_nc[1][ibf-1];
								auto kh=pgv->Gn[ic-1].ijk[2]-pgv->dijk_nc[2][ibf-1];
								if(pgv->i2c[ih-pgv->b->imino_][jh-pgv->b->jmino_][kh-pgv->b->kmino_]>0){
									r2[dirs[idir]-1][icg-1] = \
								  -r2[dirs[idir]-1][pgv->ioff+pgv->i2c[ih-pgv->b->imino_]\
									  [jh-pgv->b->jmino_][kh-pgv->b->kmino_]-1];
								}
								else{
									r2[dirs[idir]-1][icg-1] = 0.0;
								}
							}

						}
						else if(s_type==S_ZERO){
							r2[dirs[idir]-1][icg-1] = 0.0;
						};
					}
				};
			}
		}
	}
}
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		if(ghost_nNodesNeumann[ibl+1]>ghost_nNodesNeumann[ibl]){
			pgv->setBlockPointers(ibl);
			for(size_t i=ghost_nNodesNeumann[ibl]; i<ghost_nNodesNeumann[ibl+1]; i++){
				auto ibf = ghost_Neumann[0][i];
				auto icg = pgv->ioff+pgv->i2c[ghost_Neumann[2][i]-pgv->b->imino_]\
							  [ghost_Neumann[3][i]-pgv->b->jmino_][ghost_Neumann[4][i]-pgv->b->kmino_];
				auto ic = ghost_Neumann[1][i];
				auto icb = pgv->ioff+ic;
//				std::cout<<r2_S1<<'\t'<<r2_S2<<'\t'<<icg-1<<'\t'<<icb-1<<std::endl;
//				for(auto q=0; q<3; q++) std::cout<<r2[q][icg-1]<<'\t'<<r2[q][icb-1]<<std::endl;
				for(auto q=0;q<3;q++) r2[q][icg-1] = r2[q][icb-1];
				get_bound_normals(ibf, ndirs, dirs);
				for(auto idir=0;idir<ndirs;idir++){
					if(s_type==S_GRADIENT)
						r2[dirs[idir]-1][icg-1] = -r2[dirs[idir]-1][icb-1];
					else if(s_type==S_DM){
						if(ibf==2 || ibf ==4 || ibf==6) r2[dirs[idir]-1][icg-1] = 0.0;
						else{
							auto ih=pgv->Gn[ic-1].ijk[0]-pgv->dijk_nc[ibf-1][0];
							auto jh=pgv->Gn[ic-1].ijk[1]-pgv->dijk_nc[ibf-1][1];
							auto kh=pgv->Gn[ic-1].ijk[2]-pgv->dijk_nc[ibf-1][2];
							if(pgv->i2c[ih-pgv->b->imino_][jh-pgv->b->jmino_][kh-pgv->b->kmino_]>0){
								r2[dirs[idir]-1][icg-1] = \
							  -r2[dirs[idir]-1][pgv->ioff+pgv->i2c[ih-pgv->b->imino_][jh-pgv->b->jmino_]\
																			  [kh-pgv->b->kmino_]-1];
							}
							else{
								r2[dirs[idir]-1][icg-1] = 0.0;
							}
						}
					}
					else if(s_type = S_DP){
						if(ibf==1 || ibf==3 || ibf==5) r2[dirs[idir]-1][icg-1];
						else{
							auto ih=pgv->Gn[ic-1].ijk[0]-pgv->dijk_nc[ibf-1][0];
							auto jh=pgv->Gn[ic-1].ijk[1]-pgv->dijk_nc[ibf-1][1];
							auto kh=pgv->Gn[ic-1].ijk[2]-pgv->dijk_nc[ibf-1][2];
							if(pgv->i2c[ih-pgv->b->imino_][jh-pgv->b->jmino_][kh-pgv->b->kmino_]>0){
								r2[dirs[idir]-1][icg-1] = \
							  -r2[dirs[idir]-1][pgv->ioff+pgv->i2c[ih-pgv->b->imino_][jh-pgv->b->jmino_]\
													 [kh-pgv->b->kmino_]-1];
							}
							else{
								r2[dirs[idir]-1][icg-1] = 0.0;
							}
						}
					}
					else if(s_type==S_ZERO){
						r2[dirs[idir]-1][icg-1] = 0.0;
					};
				}
			}
		}
	}
}
auto ii=0;
if(selfComm){
	auto ic = ghost_nNodesPe[pp->myrank]-1;
	if(pgv->cylindrical && s_type == S_GRADIENT){
		for(auto j=selfj-1; j<selfj+selfSize-1; j++){
			ic += 1;
			auto ibl = -pgv->sg_rank_block[ghost_index[3][ic]-1][ghost_index[4][ic]-1]\
													[ghost_index[5][ic]-1]-1;
			bb = pgv->ibl2bl[ibl].p;
			if(ghost_index[1][ic]<1){
				r2[0][pgv->ibl2buf[ibl]+bb->ijk2ic[ghost_index[0][ic]-bb->imino_]\
					  [ghost_index[1][ic]-bb->jmino_][ghost_index[2][ic]-bb->kmino_]-1] = sendBuf[0][j];
				for(auto q=1; q<3; q++){
					r2[q][pgv->ibl2buf[ibl]+bb->ijk2ic[ghost_index[0][ic]-bb->imino_]\
						  [ghost_index[1][ic]-bb->jmino_][ghost_index[2][ic]-bb->kmino_]-1] = \
						  -sendBuf[q][j];
				}
			}
			else{
				for(auto q=0; q<3; q++){
					r2[q][pgv->ibl2buf[ibl]+bb->ijk2ic[ghost_index[0][ic]-bb->imino_]\
						  [ghost_index[1][ic]-bb->jmino_][ghost_index[2][ic]-bb->kmino_]-1] =\
						  sendBuf[q][j];
				}
			}
		}
	}
	else{
		for(auto j=selfj-1; j<selfj+selfSize-1; j++){
			ic += 1;
			auto ibl = -pgv->sg_rank_block[ghost_index[3][ic]-1][ghost_index[4][ic]-1]\
													[ghost_index[5][ic]-1]-1;
			bb = pgv->ibl2bl[ibl].p;
			for(auto q=0; q<3; q++){
				r2[q][pgv->ibl2buf[ibl]+bb->ijk2ic[ghost_index[0][ic]-bb->imino_]\
					  [ghost_index[1][ic]-bb->jmino_][ghost_index[2][ic]-bb->kmino_]-1]=sendBuf[q][j];
			}
		}
	}
	ii += 1;
}
for(auto nn=ii; nn<ghost_nProcsRecv; nn++){
	pp->ierr = MPI_Probe(MPI_ANY_SOURCE,6*pp->nprocs+pp->myrank+1,\
				 				MPI_COMM_WORLD,&(pp->status));
	proc = pp->status.MPI_SOURCE;
	pp->ierr = MPI_Get_count(&(pp->status), MPI_DOUBLE, &sizeComm);
	sizeComm /= 3;
	if(pp->ierr==0){
		recvBuf = new double* [3];
		for(auto i=0; i<3; i++){
			recvBuf[i] = new double[sizeComm];
		}
	}
	else{
		pp->litError("updateGhostR2","allocation error of recvBuf");
	}
	pp->ierr = MPI_Recv(recvBuf, 3*sizeComm, MPI_DOUBLE, proc,\
							   6*pp->nprocs+pp->myrank+1, MPI_COMM_WORLD, &(pp->status));
	auto ic = ghost_nNodesPe[proc]-1;
	if(pgv->cylindrical && s_type==S_GRADIENT){
		for(auto j=0; j<sizeComm; j++){
			ic += 1;
			auto ibl=-pgv->sg_rank_block[ghost_index[3][ic]-1][ghost_index[4][ic]-1]\
												 [ghost_index[5][ic]-1]-1;
			bb = pgv->ibl2bl[ibl].p;
			if(ghost_index[1][ic]<1){
				r2[0][pgv->ibl2buf[ibl]+bb->ijk2ic[ghost_index[0][ic]-bb->imino_]\
					  [ghost_index[1][ic]-bb->jmino_][ghost_index[2][ic]-bb->kmino_]-1] = recvBuf[1][j];
				for(auto q=1; q<3; q++){
				r2[q][pgv->ibl2buf[ibl]+bb->ijk2ic[ghost_index[0][ic]-bb->imino_]\
					  [ghost_index[1][ic]][ghost_index[2][ic]]-1] = -recvBuf[q][j];
				}
			}
			else{
				for(auto q=0; q<3; q++){
					r2[q][pgv->ibl2buf[ibl]+bb->ijk2ic[ghost_index[0][ic]-bb->imino_]\
					  [ghost_index[1][ic]][ghost_index[2][ic]]-1] = recvBuf[q][j];
				}
			}
		}
	}
	else{
		for(auto j=0; j<sizeComm; j++){
			ic += 1;
			auto ibl = -pgv->sg_rank_block[ghost_index[3][ic]-1][ghost_index[4][ic]-1]\
											 [ghost_index[5][ic]-1]-1;
			bb = pgv->ibl2bl[ibl].p;
			for(auto q=0; q<3; q++){
				r2[q][pgv->ibl2buf[ibl]+bb->ijk2ic[ghost_index[0][ic]-bb->imino_]\
				  [ghost_index[1][ic]][ghost_index[2][ic]]-1] = recvBuf[q][j];
			}
		}
	}
	if(pp->ierr==0){
		for(auto i=0; i<3; i++)	delete [] recvBuf[i];
		delete [] recvBuf;
		recvBuf = nullptr;
	}
	else{
		pp->parallel_die("deallocation error of recvBuf in updateGhostNodes");
	}
}
sizeComm = nMessagesSend;
for(auto nn=0; nn<nMessagesSend; nn++){
	pp->ierr = MPI_Waitany(sizeComm, requests, &proc, &(pp->status));
}
if(pp->ierr==0){
	for(auto i=0; i<3; i++) delete [] sendBuf[i];
	delete [] sendBuf;
}
else{
	pp->litError("updateGhostR2", "deallocation error at sendBuf");
}
};

void bound::get_bound_normals(const int ibf, int &ndirs, vector<int> &dirs){
if(ibf>=1 && ibf<=6){
	ndirs = 1;
	dirs[0] = (ibf+1)/2;
}
else if(ibf==7 || ibf==8 || ibf==11 || ibf==12){
	ndirs = 2;
	dirs[0] = 1;
	dirs[1] = 3;
}
else if(ibf==9 || ibf==10 || ibf==13 || ibf==14){
	ndirs = 2;
	dirs[0] = 2;
	dirs[1] = 3;
}
else if(ibf>=15 && ibf<=18){
	ndirs = 2;
	dirs[0] = 1;
	dirs[1] = 2;
}
else if(ibf>=19 && ibf<=26){
	ndirs = 3;
	dirs[0] = 1;
	dirs[1] = 2;
	dirs[2] = 3;
};
};
vector<int> bound::bound_nNodesPe;
vector<int> bound::ghost_nNodesPe; 
vector<int> bound::ghost_nNodesNeumann;
vector<vector<int>> bound::ghost_Neumann;
vector<vector<int>> bound::ghost_index;
vector<vector<int>> bound::bound_Gnodes;
int bound::ghost_nProcsSend = 0;
int bound::ghost_nProcsRecv = 0;
int bound::ghost_nNodesRecv = 0;
vector<vector<int>> bound::bf2pe;
