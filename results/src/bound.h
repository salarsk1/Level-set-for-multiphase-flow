//Written by Salar Safarkhani

#pragma once
#include"datag.h"
#include"parallel.h"
#include"timing.h"
#include"monitor.h"
#include"gnodes.h"
#include"bl.h"
#include<vector>
using std::vector;
class bound{
public:
	void bound_m_init();
	void prepareGhostNodes();
	void updateGhostNodes();
	void updateGhostR1(double r1[], int stype=0);
	void updateGhostI1(int i1[], int stype=0);
	void updateGhostR2(double **, size_t, size_t, int stype=0);
//	void updateGhostR2_v2(double **, int stype=0);
	void get_bound_normals(const int ibf, int &ndirs, vector<int> &dirs);
	const int S_SCALAR   = 1;
	const int S_GRADIENT = 2;
	const int S_ZERO     = 3;
	const int S_DM       = 4;
	const int S_DP       = 5;
	const int S_KEEP     = 6;
private:
	static int ghost_nProcsSend;
	static int ghost_nProcsRecv;
	static int ghost_nNodesRecv;
	static vector<int>  bound_nNodesPe;
	static vector<int> ghost_nNodesPe;
	static vector<int> ghost_nNodesNeumann;
	static vector<vector<int>> ghost_Neumann;
	static vector<vector<int>> ghost_index;
	static vector<vector<int>> bound_Gnodes;
	static vector<vector<int>> bf2pe;
	const int LOCAL  = 1;
	const int REMOTE = 2;
	const int NEUMAN = 3;
};
