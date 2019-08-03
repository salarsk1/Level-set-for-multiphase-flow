//Written by Salar Safarkhani

#pragma once
#include "datag.h"
#include<mpi.h>
#include<memory>
using std::vector;
using std::string;
class parallel {
public:
	void parallel_start(int argc, char *argv[], MPI_Comm communicator = 0);
	void parallel_kill(int icode);
	void parallel_die(string str);
	void parallel_barrier();
	void parallel_end();
	void MPI_ierr_Handler(int from = 0);
	void parallel_buffer_attach();
	void parallel_buffer_detach();
	void parallel_sum(int i, int &isum);
	void parallel_sum(double r,double &rsum);
	void parallel_sum(vector<double> &rv1, vector<double> &rv1sum);
	void parallel_all_sum(int i, int &isum);
	void parallel_all_sum(vector<vector<int>> &iv2,vector<vector<int>> &iv2sum);
	void parallel_all_sum(double &r, double &rsum);
	void parallel_min(int i, int &imin);
	void parallel_min(double r, double &rmin);
	void parallel_min(vector<double> &rv1, vector<double> &rv1min);
	void parallel_all_min(double r, double &rmin);
	void parallel_all_min(int i, int &imin);
	void parallel_all_min(vector<double> &r, vector<double> &rv1min);
	void parallel_max(int i, int &imax);
	void parallel_max(double r, double &rmax);
	void parallel_all_max(int i, int &imax);
	void parallel_all_max(double r, double &rmax);
	void parallel_all_or(bool &my_l, bool &l);
	void parallel_gather(int &i, vector<int> &ivec);
	void parallel_all_gather(int &i, vector<int> &ivec);
	void parallel_sync_sg_active();
	void parallel_BCast(bool &l);
	void parallel_BCast(int &i);
	void parallel_BCast(double &r);
	void parallel_BCast(char &c1);
	void parallel_BCast(string &c50);
	void parallel_BCast(vector<bool> &lv);
	void parallel_BCast(vector<int> &iv);
	void parallel_BCast(vector<double> &rv);
	void parallel_BCast(vector<vector<double>> &rv2);
	void parallel_BCast(vector<string> &c50v);
	void parallel_block_send(block_t* bb, int irk_target);
	void parallel_block_recv(block_t *bb, int irk_source);
	const int nodesInBand(char band);
	const int nodesInBandAll(char band);
	const int nodesInBandMaxAll(char band);
	void Gnode_imbalance(double &imbalance);
	void litError(string errorLocation, string message);
	static int ierr;
	static int myrank;
	static int nprocs;
	static string cmy;
	static double max_load_imbalance;
	static MPI_Comm MY_G_WORLD; //  This is true only when MPI_G_WORLD = MPI_COMM_WORLD
	static MPI_Info MPI_IO_LIT_FILE_HINT; // checked
	static MPI_Offset disp; // checked 
	static string file_prefix;
	static int disp_header;
	static MPI_Status status;
private:
	MPI_Datatype MPI_GNODE_SHORT; // do not need to be static
	MPI_Status mpi_status; // do not need to be static
	string cmyrank;
	const int MPI_TAG_GHOSTCLOTH1 = 1;
	const int MPI_TAG_GHOSTCLOTH2 = 2;
	const int MPI_TAG_GHOSTCLOTH3 = 3;
	const int MPI_TAG_BOUND1 = 4;
	const int MPI_TAG_BOUND2 = 5;
	const int MPI_TAG_BOUND3 = 6;
	const int MPI_TAG_BOUND4 = 7;
	const int MPI_TAG_BOUND5 = 8;
	const int MPI_TAG_BOUND6 = 9;
	const int MPI_TAG_BOUND7 = 20;
	const int MPI_TAG_BLOCK1 = 10;
	const int MPI_TAG_BLOCK2 = 11;
	const int MPI_TAG_BLOCK3 = 12;
	const int MPI_TAG_BLOCK4 = 13;
	const int MPI_TAG_BLOCK5 = 14;
	const int MPI_TAG_BLOCK6 = 15;
	const int MPI_TAG_BLOCK7 = 16;
};

