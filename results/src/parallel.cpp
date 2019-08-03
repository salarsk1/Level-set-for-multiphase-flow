//Written by Salar Safarkhani

#include "parallel.h"
#include<mpi.h>
#include<iostream>
#include "datag.h"
#include<memory>
#include<algorithm>
#include<array>
using std::vector;
void parallel::parallel_start(int argc, char*argv[], MPI_Comm communicator) {
#ifdef MPIDEFS2
	Gnode_short_t Gsh[3];
	MPI_Datatype types[3], types2[6];
	int lengths[3], length2[6];
	MPI_Aint base, displacement[3], displacement2[6];
#else
	int blockcounts[2];
	MPI_Aint offsets[2];
	MPI_Datatype oldtypes[2];
	MPI_Aint extent;
	MPI_Request request;	//checked;
#endif
	if(communicator)
		MY_G_WORLD = communicator;
	else {
		MPI_Init(&argc, &argv);
		MY_G_WORLD = MPI_COMM_WORLD;
	};
	ierr = MPI_Comm_size(MY_G_WORLD, &nprocs);
	ierr = MPI_Comm_rank(MY_G_WORLD, &myrank);
	cmyrank = std::to_string(myrank);
	cmy = '['+std::to_string(myrank)+']';
#ifdef MPIDEFS2
	ierr=MPI_Barrier(MY_G_WORLD);
	types[0] = MPI_DOUBLE_PRECISION;
	types[1] = MPI_INTEGER;
	types[2] = MPI_UB;
	lengths[0] = 1;
	lengths[1] = 3;
	lengths[2] = 2;
	ierr=MPI_Get_address(&Gsh[0], &base);
	ierr=MPI_Get_address(&(Gsh[0].G), &displacement[0]);
	ierr=MPI_Get_address(&(Gsh[0].ijk[0]), &displacement[1]);
	ierr=MPI_Get_address(&Gsh[1], &displacement[2]);
	for(auto &i:displacement)
		i = i - base;


	ierr=MPI_Type_create_struct(3, lengths, displacement, types, &MPI_GNODE_SHORT);
	ierr=MPI_Type_commit(&MPI_GNODE_SHORT); 
	if(!ierr != 0) parallel_kill(99);
#else
	offsets    [0] = 0;
	oldtypes   [0] = MPI_DOUBLE;
	blockcounts[0] = 1;
	ierr=MPI_Type_extent(MPI_DOUBLE, &extent);
	offsets    [1] = blockcounts[0] * extent;
	oldtypes   [1] = MPI_INTEGER;
	blockcounts[1] = 3;
	ierr=MPI_Type_struct(2, blockcounts, offsets, oldtypes, &MPI_GNODE_SHORT);
	ierr=MPI_Type_commit(&MPI_GNODE_SHORT);
#endif
#ifdef PANFS_CUNCURRENT_WRITE
	ierr=MPI_Info_create(&MPI_IO_LIT_FILE_HINT);
	char temp = '1';
	ierr=MPI_Info_set(MPI_IO_LIT_FILE_HINT, "pasnfs_concurrent_write", &temp);
	file_prefix = "panfs:";
#else
	MPI_IO_LIT_FILE_HINT = MPI_INFO_NULL;
	file_prefix = "ufs:";
#endif
};
void parallel::parallel_kill(int icode) {
	std::unique_ptr<global_variable>pgv(new global_variable);
	std::cout << pgv -> clit << "*******************************" << std::endl;
	std::cout<< pgv-> clit<<"pe#"<<myrank<<" will kill job due to error code "<<icode<< std::endl;
	ierr=MPI_Abort(MPI_COMM_WORLD, icode);
	exit(1);
};
void parallel::parallel_die(std::string str) {
	std::unique_ptr<global_variable> pgv1(new global_variable);
	std::cout << pgv1 -> clit <<"******************************"<<std::endl;
	std::cout << pgv1 -> clit <<"pe#"<<myrank<<"will kill job due to error"<<str<<std::endl;
	exit(1);
};
void parallel::parallel_barrier() {
	ierr=MPI_Barrier(MPI_COMM_WORLD);
};
void parallel::parallel_end() {
	ierr=MPI_Barrier(MPI_COMM_WORLD);
	ierr=MPI_Finalize();
};
void parallel::MPI_ierr_Handler(int from) {
	std::unique_ptr<global_variable>pgv(new global_variable);
	char str[MPI_MAX_ERROR_STRING];
	int lenstr;
	int eclass;//different from fortran; checked
	if(ierr != 0) {
		MPI_Error_class(ierr, &eclass);
		ierr = MPI_Error_string(ierr, str, &lenstr);
		std::cout<<pgv -> clit<<"["<<myrank<<"]"<<"MPI_error is"<<str<<std::endl;
		std::cout<<pgv -> clit<<"["<<myrank<<"]"<<"MPI_error occured at position = "<<from<<std::endl;
		parallel_kill(0);
	};
};
void parallel::parallel_buffer_attach() {
	std::unique_ptr<global_variable>pgv(new global_variable);
	if(pgv -> verbose && myrank == 0) {
		std::cout<<pgv->clit<<"WARNING: not attaching any MPI Buffer"<<std::endl;
		std::cout<<pgv->clit<<" if FMM must be used, attach buffer buffer again"<<std::endl;
	};
};
void parallel::parallel_buffer_detach() {
	int NewMessage;
	int itag, ipe;
	std::unique_ptr<global_variable>pgv(new global_variable);
	ierr=MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &NewMessage, &mpi_status);
	MPI_ierr_Handler();
	if(NewMessage) {
		itag = mpi_status.MPI_TAG; // checked
		ipe  = mpi_status.MPI_SOURCE; //checked
		std::cout<<pgv->clit<<cmy<<"strange, still have a message in the buffer befor detach";
		std::cout<<std::endl;
	};
};
// The following functions will change to template functions

void parallel::parallel_sum(int i, int &isum){
	ierr=MPI_Reduce(&i, &isum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
};
void parallel::parallel_sum(double r, double &rsum) {
	ierr=MPI_Reduce(&r, &rsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
};
void parallel::parallel_sum(vector<double> &rv1,vector<double> &rv1sum) {
	ierr=MPI_Reduce(&rv1[0],&rv1sum[0],rv1.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
};

void parallel::parallel_all_sum(int i, int &isum) {
	ierr=MPI_Allreduce(&i, &isum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
};
void parallel::parallel_all_sum(vector<vector<int>> &iv2,vector<vector<int>> &iv2sum){
	ierr=MPI_Allreduce(&iv2[0][0],&iv2sum[0][0],iv2.size()*iv2[0].size(),\
							 MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
};
void parallel::parallel_all_sum(double &r, double &rsum){
	ierr=MPI_Allreduce(&r,&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
};

void parallel::parallel_min(int i, int &imin) {
	ierr=MPI_Reduce(&i, &imin, 1, MPI_INTEGER, MPI_MIN, 0, MPI_COMM_WORLD);
};
void parallel::parallel_min(double r, double &rmin) {
	ierr=MPI_Reduce(&r, &rmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
};
void parallel::parallel_min(vector<double> &rv1,vector<double> &rv1min) {
	ierr=MPI_Reduce(&rv1[0],&rv1min[0],rv1.size(),MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
};

void parallel::parallel_all_min(double r, double &rmin){
	ierr=MPI_Allreduce(&r, &rmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
};
void parallel::parallel_all_min(int i, int &imin){
	ierr=MPI_Allreduce(&i, &imin, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD);
};
void parallel::parallel_all_min(vector<double> &rv1, vector<double> &rv1min){
	ierr=MPI_Allreduce(&rv1[0], &rv1min[0], rv1.size(), MPI_INTEGER, MPI_MIN,MPI_COMM_WORLD);
};

void parallel::parallel_max(int i, int &imax){
	ierr=MPI_Reduce(&i, &imax, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD);
};
void parallel::parallel_max(double r, double &rmax){
	ierr=MPI_Reduce(&r, &rmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
};
void parallel::parallel_all_max(int i, int &imax){
	ierr=MPI_Allreduce(&i, &imax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);
};
void parallel::parallel_all_max(double r, double &rmax){
	ierr=MPI_Allreduce(&r, &rmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
};
void parallel::parallel_all_or(bool &my_l, bool &l){
	ierr=MPI_Allreduce(&my_l, &l, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD);
};
void parallel::parallel_gather(int &i, vector<int> &ivec){
	ierr=MPI_Gather(&i, 1, MPI_INTEGER, &ivec, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
};
void parallel::parallel_all_gather(int &i, vector<int> &ivec){
	ierr=MPI_Allgather(&i, 1, MPI_INTEGER, &ivec, 1, MPI_INTEGER, MPI_COMM_WORLD);
};
void parallel::parallel_sync_sg_active(){
	std::unique_ptr<global_variable>pgv(new global_variable);
	const size_t d1 = pgv -> ijkm_sg[0];
	const size_t d2 = pgv -> ijkm_sg[1];
	const size_t d3 = pgv -> ijkm_sg[2];
	bool ***l_result;
	if(ierr == 0){
		l_result = new bool**[d1];
		for(auto i=0; i<d1; i++) {
			l_result[i] = new bool *[d2];
			for(auto j=0; j<d2; j++) {
				l_result[i][j] = new bool [d3];
			}
		}
	}
	else{
		litError("Parallel_sync_sg_active","allocation error for l_result");
	};
	size_t NNN = pgv->ijkm_sg[0]*pgv->ijkm_sg[1]*pgv->ijkm_sg[2];
	ierr = MPI_Allreduce(&pgv->sg_active[0][0][0],&l_result[0][0][0],NNN, MPI_LOGICAL,\
								 MPI_LOR, MPI_COMM_WORLD);
	for(int i=0;i<d1;++i){
		for(int j=0;j<d2;++j){
			for(int k=0;k<d3;++k){
				pgv->sg_active[i][j][k] = l_result[i][j][k];
			}
		}
	}
	if(ierr==0){
		for(int i=0; i<d1; ++i){
			for(int j=0; j<d2; ++j){
				delete [] l_result[i][j];
			};
			delete [] l_result[i];
		}
		delete [] l_result;
	}
	else litError("parallel_sync_sg_active", "deallocation error for l_result");
}
void parallel::parallel_BCast(bool &l){
	ierr = MPI_Bcast(&l, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
};
void parallel::parallel_BCast(int &i){
	ierr = MPI_Bcast(&i, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
};
void parallel::parallel_BCast(double &r){
	ierr = MPI_Bcast(&r, 1, MPI_DOUBLE, 0 ,MPI_COMM_WORLD);
};
void parallel::parallel_BCast(char &c1){
	int ii;
	if(myrank ==0)	ii = static_cast<int>(c1);
	ierr = MPI_Bcast(&ii, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
	if(myrank != 0) c1 = static_cast<char>(ii);
};
void parallel::parallel_BCast(string &c50){
int ic50[50];
int ii;
if(myrank ==0) {
	for(ii=0; ii<50; ii++){
		ic50[ii] = static_cast<int>(c50[ii]);
	};
};
ierr = MPI_Bcast(ic50, 50, MPI_INTEGER, 0, MPI_COMM_WORLD);
if(myrank != 0){
	for(ii=0; ii<50; ii++){
		c50[ii] == static_cast<char>(ic50[ii]);
	};
};
};
void parallel::parallel_BCast(vector<bool> &lv){
ierr = MPI_Bcast(&lv, lv.size(), MPI_LOGICAL, 0, MPI_COMM_WORLD);
};
void parallel::parallel_BCast(vector<int> &iv){
ierr = MPI_Bcast(&iv, iv.size(), MPI_INTEGER, 0, MPI_COMM_WORLD);
};
void parallel::parallel_BCast(vector<double> &rv){
ierr = MPI_Bcast(&rv, rv.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
};
void parallel::parallel_BCast(vector<vector<double>> &rv2){
ierr = MPI_Bcast(&rv2, rv2.size() * rv2[0].size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
};
void parallel::parallel_BCast(vector<string> &c50v){
int ii, id;
std::array<int, 50> ic50;
for(id=0; c50v.size(); id++){
	if(myrank == 0){
		for(ii=0; ii<50; ii++){
			ic50[ii] = static_cast<int> (c50v[id][ii]);
		};
	};
};
};
void parallel::parallel_block_send(block_t* bb, int irk_target){
std::array<int, 12> dat;
int ic;
std::copy(std::begin(bb->ijkm), std::end(bb->ijkm), std::begin(dat));
std::copy(std::begin(bb->ijk_sg),std::end(bb->ijk_sg),std::begin(dat)+3);
dat[6] = bb->NinA; dat[7]  = bb->NinT; dat[8]  = bb->NinN;
dat[9] = bb->NinW; dat[10] = bb->NinX; dat[11] = bb->NinZ;
#ifdef DEBUG_MODE
ierr = MPI_Ssend(&dat, dat.size(), MPI_INTEGER, irk_target, MPI_TAG_BLOCK1, MPI_COMM_WORLD);
#else
ierr = MPI_Ssend(&dat, dat.size(), MPI_INTEGER, irk_target, MPI_TAG_BLOCK1, MPI_COMM_WORLD);
#endif
if(ierr != 0){
	vector<double> sendBuf(bb->NinZ);
	vector<vector<double>> sendBufInt(3, vector<double>(bb->NinZ));
for(ic=0; ic < bb->NinZ; ic++){
	sendBuf      [ic] = (bb->Gnodes)[ic].G;
	sendBufInt[0][ic] = (bb->Gnodes)[ic].ijk[0];
	sendBufInt[1][ic] = (bb->Gnodes)[ic].ijk[1];
	sendBufInt[2][ic] = (bb->Gnodes)[ic].ijk[2];
};
#ifdef DEBUG_MODE
ierr=MPI_Ssend(sendBuf.data(), bb->NinZ, MPI_DOUBLE, irk_target, MPI_TAG_BLOCK2,MPI_COMM_WORLD);
ierr=MPI_Ssend(sendBufInt.data(),3*bb->NinZ,MPI_INTEGER,irk_target,MPI_TAG_BLOCK3,MPI_COMM_WORLD);
#else
ierr=MPI_Send(sendBuf.data(), bb->NinZ, MPI_DOUBLE, irk_target, MPI_TAG_BLOCK2, MPI_COMM_WORLD);
ierr=MPI_Send(sendBufInt.data(),3*bb->NinZ,MPI_INTEGER,irk_target,MPI_TAG_BLOCK3,MPI_COMM_WORLD);
#endif
sendBuf.clear();
sendBufInt.clear();
sendBuf.shrink_to_fit();
sendBufInt.shrink_to_fit();
}
else
	litError("parallel_block_send","cannot allocate sendBuf, sendBufInt");
};
void parallel::parallel_block_recv(block_t *bb, int irk_source){
vector<double> recvBuf;
vector<vector<int>> recvBufInt;
std::array<int, 12> dat;
int ibf;
std::unique_ptr<global_variable>pgv(new global_variable);
ierr=MPI_Recv(&dat, dat.size(), MPI_INTEGER, irk_source, MPI_TAG_BLOCK1,MPI_COMM_WORLD, &mpi_status);
std::copy(std::begin(dat), std::begin(dat)+3, std::begin(bb->ijkm));
std::copy(std::begin(dat)+3, std::begin(dat)+6, std::begin(bb->ijk_sg));
bb->NinA = dat[6];
bb->NinT = dat[7];
bb->NinN = dat[8];
bb->NinW = dat[9];
bb->NinX = dat[10];
bb->NinZ = dat[11];
std::copy(std::begin(dat)+6, std::begin(dat)+9, std::begin(bb->NinBand));
bb->NinBand[3] = dat[9];
bb->NinBand[4] = dat[10];
std::fill(bb->nGnodeb.begin(), bb->nGnodeb.end(), 0);
for(auto i=0; i<3; i++) {
	bb->ijk0[i] = (bb->ijk_sg[i] -1)*(pgv->ijkm_bl[i]);
};
bb->imin_ = (bb->ijk_sg[0] - 1)*pgv->ijkm_bl[0] + 1;
bb->imax_ = bb->ijk_sg[0]*pgv->ijkm_bl[0];
bb->jmin_ = (bb->ijk_sg[1] - 1)*pgv->ijkm_bl[1] + 1;
bb->imax_ = bb->ijk_sg[1]*pgv->ijkm_bl[1];
bb->kmin_ = (bb->ijk_sg[2] - 1)*pgv->ijkm_bl[2] + 1;
bb->kmax_ = bb->ijk_sg[2]*pgv->ijkm_bl[2];
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
for (auto j=0; j<3; j++){
	for(ibf=0; ibf<26; ibf++){
		bb->ijk_b2g[j][ibf]=2*((bb->ijk0)[1]+((pgv->dijk_nc)[ibf][j]+1)/2 * pgv->ijkm_bl[j]);
	}
};
auto d1 = bb->imaxo_ - bb->imino_+1;
auto d2 = bb->jmaxo_ - bb->jmino_+1;
auto d3 = bb->kmaxo_ - bb->kmino_+1;
if(ierr==0){
	bb->ijk2ic = new int**[d1];
	for(size_t i=0; i<d1; i++) {
		(bb->ijk2ic)[i] = new int *[d2];
		for(size_t j=0; i<d2; i++) {
			(bb->ijk2ic)[i][j] = new int [d3];
		}
	}
}
else
	litError("bl_activate_new","Allocation error of bb->ijk2ic");
if(ierr==0){
	bb->ic0_bl.reset(new int[pgv->nBandLayers]);
	bb->NinBandLayer.reset(new int[pgv->nBandLayers]);
};
if(ierr==0){
	recvBuf.resize(bb->NinZ);
	recvBufInt.resize(3);
	for(auto i=0; i<bb->NinZ; i++){
		recvBufInt[i].resize(bb->NinZ);
	};
}
else litError("parallel_block_recv","cannot allocate recvBuf, recvBufInt.");
ierr=MPI_Recv(recvBuf.data(), bb->NinZ, MPI_DOUBLE, irk_source, MPI_TAG_BLOCK2,\
				  MPI_COMM_WORLD, &mpi_status);
ierr=MPI_Recv(recvBufInt.data(), 3*bb->NinZ, MPI_DOUBLE, irk_source,\
		   MPI_TAG_BLOCK3, MPI_COMM_WORLD, &mpi_status);
bb->nGmax = static_cast<int>(bb->NinZ*1.1);

if(bb->Gnodes != nullptr){
	if(ierr==0){
		delete [] bb->Gnodes;
	}
	else
		litError("parallel_block_recv","cannot deallocate bb->Gnodes.");
};
if(ierr==0){
	bb->Gnodes = new Gnode_t[bb->nGmax];
}
else{
	litError("parallel_block_recv","cannot allocate bb->Gnodes.");
};
for(auto ic=0; ic < bb->NinZ; ic++){
	bb->Gnodes[ic].G = recvBuf.at(ic);
	bb->Gnodes[ic].V = {0};
	bb->Gnodes[ic].ijk[0] = recvBufInt[0][ic];
	bb->Gnodes[ic].ijk[1] = recvBufInt[1][ic];
	bb->Gnodes[ic].ijk[2] = recvBufInt[2][ic];
	bb->Gnodes[ic].dijk[0] = 0;
	bb->Gnodes[ic].dijk[1] = 0;
	bb->Gnodes[ic].dijk[2] = 0;
	auto d1 = recvBufInt[0][ic]-pgv->b->imino_;
	auto d2 = recvBufInt[1][ic]-pgv->b->jmino_;
	auto d3 = recvBufInt[2][ic]-pgv->b->kmino_;
	bb->ijk2ic[d1][d2][d3] = ic+1;
};
if(ierr==0){
	recvBuf.clear();
	recvBufInt.clear();
	recvBuf.shrink_to_fit();
	recvBufInt.shrink_to_fit();
}
else{
	litError("parallel_block_recv","cannot deallocate recvBuf,recvBufInt");
};
};

const int parallel::nodesInBand(char band){
int myN = 0;
int n;
std::unique_ptr<global_variable>pgv(new global_variable);
switch(band){
case 'A':
	for(auto ibl=0; ibl<pgv->nbl;ibl++){
		myN += static_cast<int>(pgv->ibl2bl[ibl].p->NinA);
	};
break;
case 'T':
	for(auto ibl=0; ibl<pgv->nbl;ibl++){
		myN += static_cast<int>(pgv->ibl2bl[ibl].p->NinT);
	};
break;
case 'N':
	for(auto ibl=0; ibl<pgv->nbl;ibl++){
		myN += static_cast<int>(pgv->ibl2bl[ibl].p->NinN);
	};
break;
case 'X':
	for(auto ibl=0; ibl<pgv->nbl;ibl++){
		myN += static_cast<int>(pgv->ibl2bl[ibl].p->NinX);
	};
break;
case 'Z':
	for(auto ibl=0; ibl<pgv->nbl;ibl++){
		myN += static_cast<int>(pgv->ibl2bl[ibl].p->NinZ);
	};
break;
};
ierr = MPI_Reduce(&myN, &n, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
return n;
};
const int parallel::nodesInBandAll(char band){
int n;
int myN = 0;
std::unique_ptr<global_variable>pgv(new global_variable);
switch(band){
case 'A':
	for(auto ibl=0;ibl<pgv->nbl;ibl++){
		myN = myN + pgv->ibl2bl[ibl].p->NinA;
	};
break;
case 'T':
	for(auto ibl=0;ibl<pgv->nbl;ibl++){
		myN = myN + pgv->ibl2bl[ibl].p->NinT;
	};
break;
case 'N':
	for(auto ibl=0;ibl<pgv->nbl;ibl++){
		myN = myN + pgv->ibl2bl[ibl].p->NinN;
	};
break;
case 'Z':
	for(auto ibl=0;ibl<pgv->nbl;ibl++){
		myN = myN + pgv->ibl2bl[ibl].p->NinZ;
	};
break;
};
parallel_all_sum(myN, n);
return(n);
};
const int parallel::nodesInBandMaxAll(char band){
int n;
int myN = 0;
std::unique_ptr<global_variable>pgv(new global_variable);
switch(band){
case 'A':
	for(auto ibl=0;ibl<pgv->nbl;ibl++){
		myN = myN + pgv->ibl2bl[ibl].p->NinA;
	};
break;
case 'T':
	for(auto ibl=0;ibl<pgv->nbl;ibl++){
		myN = myN + pgv->ibl2bl[ibl].p->NinT;
	};
break;
case 'N':
	for(auto ibl=0;ibl<pgv->nbl;ibl++){
		myN = myN + pgv->ibl2bl[ibl].p->NinN;
	};
break;
case 'Z':
	for(auto ibl=0;ibl<pgv->nbl;ibl++){
		myN = myN + pgv->ibl2bl[ibl].p->NinZ;
	};
break;
};
parallel_all_sum(myN, n);
return n;
};
void parallel::Gnode_imbalance(double &imbalance){
int nSum, nMax;
int myN = 0;
std::unique_ptr<global_variable>pgv(new global_variable);
for(auto ibl=0; ibl<pgv->nbl; ibl++) myN = myN + pgv->ibl2bl[ibl].p->NinN;
ierr = MPI_Reduce(&myN, &nMax, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD);
ierr = MPI_Reduce(&myN, &nSum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);
if(myrank==0) imbalance = static_cast<double>(nprocs*nMax/nSum);
ierr = MPI_Bcast(&imbalance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
};
void parallel::litError(string errorLocation, string message){
	std::unique_ptr<global_variable>pgv(new global_variable);
	std::cout<<pgv->clit<<"Error! Aborting lit! myrank = "<<myrank<<std::endl;
	std::cout<<pgv->clit<<"Error occured in subroutine "<<errorLocation<< std::endl;
	std::cout<<pgv->clit<<"Error message: "<<message<<std::endl;
parallel_kill(0);
};
int parallel::ierr = 0;
int parallel::nprocs = 0;
int parallel::myrank = 0;
string parallel::cmy;
double parallel::max_load_imbalance = 0.0;
string parallel::file_prefix;
int parallel::disp_header = 0;
MPI_Status parallel::status;
MPI_Comm parallel::MY_G_WORLD = 0;
MPI_Info parallel::MPI_IO_LIT_FILE_HINT = 0;
MPI_Offset parallel::disp = 0;
