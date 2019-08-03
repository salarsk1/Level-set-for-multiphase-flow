//Written by Salar Safarkhani

#include "litBuffer.h"
#include<algorithm>
#include<iostream>
void litBuffer::litBuffer_m_init(){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
if(pp->ierr == 0){
	pgv->ibl2buf = new int[pgv->max_bl];
}
else{
	pp->litError("litBuffer_m_init","Cannot allocate ibl2buf");
};
};

bool* litBuffer::getL1Buffer(size_t &n, char band, int* i2b){
int bandID;
int *i2bh;
bool* l1;
block_t *bb;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
if(band == 'A') bandID = 1;
else if(band == 'T') bandID = 2;
else if(band == 'N') bandID = 3;
else if(band == 'W') bandID = 4;
else if(band == 'X') bandID = 5;
else if(band == 'Z') bandID = 6;
else if(band == 'F') bandID = 6;
else pp->litError("getL1Buffer","unknown band name");
if (i2b != nullptr) i2bh = i2b;
else{
	if(bandID==6) i2bh = pgv->ibl2buf;
	else pp->litError("getBuffer","must provide ibl2buf when using band\
							 other than Z.");
};
n = 0;
for(auto ibl = 0; ibl<pgv->nbl; ibl++){
	bb = pgv->ibl2bl[ibl].p;
	i2bh[ibl] = n;
	if(bandID == 6) n += bb->NinZ;
	else n += bb->NinBand[bandID-1];
};
if(pp->ierr == 0)
	l1 = new bool[n];
else
	pp->litError("getL1Buffer","Memory allocation failur for l1buffer");
for(auto i=0; i<n; i++) l1[i] = false;
return l1;
};

int* litBuffer::getI1Buffer(size_t &n, char band, int* i2b){
int ibl, bandID;
int *i2bh;
int* i1;
block_t *bb;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
if(band == 'A') bandID = 1;
else if(band == 'T') bandID = 2;
else if(band == 'N') bandID = 3;
else if(band == 'W') bandID = 4;
else if(band == 'X') bandID = 5;
else if(band == 'Z') bandID = 6;
else if(band == 'F') bandID = 6;
else pp->litError("getL1Buffer","unknown band name");
if (i2b != nullptr) i2bh = i2b;
else{
	if(bandID==6) i2bh = pgv->ibl2buf;
	else pp->litError("getBuffer","must provide ibl2buf when using band\
							 other than Z.");
}
n = 0;
for(auto ibl = 0; ibl<pgv->nbl; ibl++){
	bb = pgv->ibl2bl[ibl].p;
	i2bh[ibl] = n;
	if(bandID == 6) n += bb->NinZ;
	else n += bb->NinBand[bandID-1];
};
if(pp->ierr == 0)
	i1 = new int[n];
else
	pp->litError("getI1Buffer","Memory allocation failur for i1buffer");
for(auto i=0; i<n; i++) i1[i] = 0;
return i1;
};

double* litBuffer::getR1Buffer(size_t &n, char band, int* i2b){
int ibl, bandID;
int *i2bh;
double* r1;
block_t *bb;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
if(band == 'A') bandID = 1;
else if(band == 'T') bandID = 2;
else if(band == 'N') bandID = 3;
else if(band == 'W') bandID = 4;
else if(band == 'X') bandID = 5;
else if(band == 'Z') bandID = 6;
else if(band == 'F') bandID = 6;
else pp->litError("getL1Buffer","unknown band name");
if (i2b != nullptr) i2bh = i2b;
else{
	if(bandID==6) i2bh = pgv->ibl2buf;
	else pp->litError("getBuffer","must provide ibl2buf when using band other than Z.");
}
n = 0;
for(auto ibl = 0; ibl<pgv->nbl; ibl++){
	bb = pgv->ibl2bl[ibl].p;
	i2bh[ibl] = n;
	if(bandID == 6) n += bb->NinZ;
	else n += bb->NinBand[bandID-1];
};
if(pp->ierr == 0) r1 = new double[n];
else pp->litError("getI1Buffer","Memory allocation failur for r1buffer");
for(auto i=0; i<n; i++) r1[i] = 0.0;
return r1;
};

double** litBuffer::getR2Buffer(size_t&n2, size_t&n,char band, int nr2, int* i2b){
int* i2bh;
int bandID;
block_t* bb;
double** r2;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
if(band == 'A') bandID = 1;
else if(band == 'T') bandID = 2;
else if(band == 'N') bandID = 3;
else if(band == 'W') bandID = 4;
else if(band == 'X') bandID = 5;
else if(band == 'Z') bandID = 6;
else if(band == 'F') bandID = 6;
else pp->litError("getL1Buffer","unknown band name");
if (i2b != nullptr) i2bh = i2b;
else{
	if(bandID==6) i2bh = pgv->ibl2buf;
	else pp->litError("getBuffer","must provide ibl2buf when using band\
							 other than Z.");
}
if(nr2 != 0) n2 = nr2;
else n2 = 3;
n = 0;
for(auto ibl = 0; ibl<pgv->nbl; ibl++){
	bb = pgv->ibl2bl[ibl].p;
	i2bh[ibl] = n;
	if(bandID == 6) n += bb->NinZ;
	else n += bb->NinBand[bandID-1];
};
if(pp->ierr == 0){
	r2 = new double*[n2];
	for(auto i=0; i<n2; i++){
		r2[i] = new double[n];
	}
}
else{
	pp->litError("getI1Buffer","Memory allocation failur for r1buffer");
};
for(auto i=0; i<n2; i++){
	for(auto j=0; j<n; j++){
		r2[i][j] = 0.0;
	}
}
return r2;
};

void litBuffer::freeL1Buffer(bool *l1){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
if(l1 != nullptr){
	if(pp->ierr==0){
		delete [] l1;
		l1 = nullptr;
	}
	else{
 		pp->litError("freeL1Buffer","Memory deallocation\
						  failure for l1buffer");
	}
}
else{
	std::cout<<pgv->clit<<"Warning: not deallocation l1buffer";
	std::cout<<"pointer is not associated"<<std::endl;
};
};

void litBuffer::freeI1Buffer(int *i1){
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<global_variable>pgv(new global_variable);
if(i1 != nullptr){
	if(pp->ierr==0){
		delete [] i1;
		i1 = nullptr;
	}
	else{
 		pp->litError("freeI1Buffer","Memory deallocation\
						  failure for i1buffer");
	}
}
else{
	std::cout<<pgv->clit<<"Warning: not deallocation i1buffer";
	std::cout<<"since pointer is not associated"<<std::endl;
};
};

void litBuffer::freeR1Buffer(double *r1){
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<global_variable>pgv(new global_variable);
if(r1 != nullptr){
	if(pp->ierr==0){
		delete [] r1;
		r1 = nullptr;
	}
	else{
 		pp->litError("freeR1Buffer","Memory deallocation\
						  failure for r1buffer");
	}
}
else{
	std::cout<<pgv->clit<<"Warning: not deallocation r1buffer";
	std::cout<<"since pointer is not associated"<<std::endl;
};
};

void litBuffer::freeR2Buffer(size_t r2_S1, size_t r2_S2, double **r2){
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<global_variable>pgv(new global_variable);
if(r2 != nullptr){
	if(pp->ierr==0){
		for(auto i=0; i<r2_S1; i++){
			delete [] r2[i];
		}
		delete [] r2;
		r2 = nullptr;
	}
	else{
 		pp->litError("freeR2Buffer","Memory deallocation failure for r2buffer");
	}
}
else{
	std::cout<<pgv->clit<<"Warning: not deallocation r2buffer";
	std::cout<<"since pointer is not associated"<<std::endl;
};
};
