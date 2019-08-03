//Written by Salar Safarkhani

#include"gnodes.h"
#include<algorithm>
#include<numeric>
#include<iostream>
void gnodes::enlargenGnodes(block_t* bb, int n){
	std::unique_ptr<global_variable>pgv(new global_variable);
	std::unique_ptr<parallel>pp(new parallel);
	Gnode_t* Gh;
	bool shortcut;
	if(n !=0){
		bb->nGmax = bb->nGmax + n;
	}
	else{
		bb->nGmax = static_cast<int>(static_cast<double>(bb->nGmax)*1.1);
	};
	if(pgv->Gn == bb->Gnodes){
		shortcut = true;
		pgv->Gn = nullptr;
	}
	else shortcut = false;
	if(pp->ierr == 0)	Gh = new Gnode_t[bb->nGmax];
	else 	pp->litError("enlargenGnodes","Allocation error for Gh");

	for(auto i=0; i<bb->nG; i++)
		Gh[i] = bb->Gnodes[i];
	if(bb->Gnodes != nullptr){
		if(pp->ierr==0){
			delete [] bb->Gnodes;
			bb->Gnodes = nullptr;
		}
		else pp->litError("enlargenGnodes","Deallocation error for bb->Gnodes");
	}
	if(pp->ierr==0)
		bb->Gnodes = new Gnode_t[bb->nGmax];
	else
		pp->litError("enlargen Gnodes","Allocation Error for bb->Gnodes in resizeGnodes");
	for(auto j=0; j<pgv->b->nG; j++) bb->Gnodes[j] = Gh[j];
	if(pp->ierr==0 && Gh != nullptr){
		delete [] Gh;
		Gh = nullptr;
	}
	else
		pp->litError("enlargenGnodes","Deallocation error for Gh in resizeGnodes");
	if(shortcut)
		pgv->Gn = bb->Gnodes;
};

void gnodes::shortenGnodes(block_t* bb){
	std::unique_ptr<global_variable>pgv(new global_variable);
	std::unique_ptr<parallel>pp(new parallel);
	Gnode_t* Gh;
	bool shortcut;
	bb->nGmax = std::max(100, static_cast<int>(static_cast<double>(bb->nG)*1.1));
	if(pgv->Gn == bb->Gnodes){
		shortcut = true;
		pgv->Gn = nullptr;
	}
	else{
		shortcut = false;
	};
	if(pp->ierr==0){
		Gh = new Gnode_t[bb->nG];
	}
	else{
		pp->litError("shortenGnodes", "Allocation error for Gh");
	};
	for(auto i=0; i < bb->nG; i++) Gh[i] = bb->Gnodes[i];
	if(bb->Gnodes != nullptr){
		if(pp->ierr==0){
			delete [] bb->Gnodes;
			bb->Gnodes = nullptr;
		}
		else
			pp->litError("shortenGnodes","Deallocation error for bb->Gnodes");
	};
	if(pp->ierr==0) bb->Gnodes = new Gnode_t[bb->nGmax];
	else pp->litError("shortenGnodes", "Allocation error for bb->Gnodes");
	for(auto ic=0; ic<bb->nG; ic++) bb->Gnodes[ic] = Gh[ic];
	if(pp->ierr==0) delete[] Gh;
	else pp->litError("ShortenGnodes","Deallaocation error for Gh");
	if(shortcut) pgv->Gn = bb->Gnodes;
};

void gnodes::addGnode(block_t* bb, const Gnode_t &Gnh){
	if(bb->nG+1>bb->nGmax) enlargenGnodes(bb);
	bb->nG += 1;
	if(bb->Gnodes == nullptr) std::cout<<"huh?"<<std::endl;
	bb->Gnodes[bb->nG-1] = Gnh;
};
	
void gnodes::ensureSizeGnodes(block_t* bb, int n){
	if(n>bb->nGmax)
		enlargenGnodes(bb, n-bb->nGmax+1);
};

