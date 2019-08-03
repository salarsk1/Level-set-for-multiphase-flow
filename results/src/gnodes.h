//Written by Salar Safarkhani

#ifndef GNODES_H
#define GNODES_H
#include"datag.h"
#include"parallel.h"
class gnodes{
	public:
	void enlargenGnodes(block_t* bb, int n=0);
	void shortenGnodes(block_t* bb);
	void addGnode(block_t* bb, const Gnode_t &Gnh);
	void ensureSizeGnodes(block_t *bb, int n);
	private:
};
#endif
