#include"flow.h"
#include"calculation.h"
#include"alg.h"
#include<memory>
void flow::flow_solver(int M,int N,double t_final,const double rho_min, \
							  const vector<vector<double>> &rho_n,\
							  const vector<vector<double>> &rho,\
							  const vector<vector<double>> &mu,\
							  vector<vector<double>>&uin,\
							  vector<vector<double>>&vin, \
							  vector<vector<double>>&poo, \
							  vector<vector<double>>&po,\
							  vector<vector<double>>&pn,\
							  const vector<vector<double>>& kappa ,\
							  const vector<vector<double>> &lev){
	std::unique_ptr<calculation> pcal(new calculation);
	std::vector<std::vector<double>> u, v, un, vn;
	double dt;
	alg algo;
	int num_itr;
	double dx = 1.0 / M;
	double dy = 1.5 / N;
	u.resize(M+1);
	un.resize(M+1);
	for(auto i=0; i<M+1; i++){
		u[i].resize(N+2);
		un[i].resize(N+2);
	}
	v.resize(M+2);
	vn.resize(M+2);
	for(auto i=0; i<M+2; i++){
		v[i].resize(N+1);
		vn[i].resize(N+1);
	}
	for(auto i=0; i<M+1; i++){
		for(auto j=0; j<N+2; j++){
			u[i][j]  = 0.0;
			un[i][j] = 0.0;
		}
	}
	for(auto i=0; i<M+2; i++){
		for(auto j=0; j<N+1; j++){
			v[i][j]  = 0.0;
			vn[i][j] = 0.0;
		}
	}
	for(auto i=0; i<M+1; i++){
		for(auto j=0; j<N+2; j++){
			u[i][j] = uin[i][j];
		}
	}
	for(auto i=0; i<M+2; i++){
		for(auto j=0; j<N+1; j++){
			v[i][j] = vin[i][j];
		}
	}
	int k = 0;
	double time = 0.0;
//		pcal->calcTimeStep(rho, mu, M, N, u, v, dt);
		time += dt;
		if(time > t_final){
			time -= dt;
			dt = t_final - time;
			time += dt;
		}
//		pcal->predict(rho_n, rho, mu, t_final, M, N, u, v, un, vn, poo, po,  kappa, lev);

//		u = un; v = vn;
		un = u; vn = v;
		pcal->poisson(rho_min, M, N, t_final, un, vn, po, pn);
		//Impose Neumann Boundary Condition for poo, po
/*		for(auto i=0; i<M+2; i++){
			po[i][0] = po[i][1];
			pn[i][0] = pn[i][1];
			po[i][N+1] = po[i][N];
			pn[i][N+1] = pn[i][N];
		}
		for(auto j=0; j<N+2; j++){
			po[0][j] = po[1][j];
			po[M+1][j] = po[M][j];
			pn[0][j] = pn[1][j];
			pn[M+1][j] = pn[M][j];
		}*/
		// set poo as po
		for(auto i=0; i<M+2; i++){
			for(auto j=0; j<N+2; j++){
				poo[i][j] = po[i][j];
			}
		}
		for(auto i=0; i<M+2; i++){
			for(auto j=0; j<N+2; j++){
				po[i][j] = pn[i][j];
			}
		}
		for(auto i=0; i<M+1; i++){
			for(auto j=0; j<N+2; j++) uin[i][j] = un[i][j];
		}
		for(auto i=0; i<M+2; i++){
			for(auto j=0; j<N+1; j++) vin[i][j] = vn[i][j];
		}
}
