#include"calculation.h"
#include"alg.h"
#include<iostream>
#include<limits>
void calculation::calcTimeStep(const double Re, const int M, const int N,const std::vector<std::vector<double>> &u,\
                  const std::vector<std::vector<double>> &v, double & dt){
	alg algo;
	double a, b, Re1, Re2;
	double dx = 1.0/M;
	double dy = 1.0/N;
	double dt1 = pow(dx,2)/(4.0/Re);
	std::vector<std::vector<double>> dt2, dt3;
	dt2.resize(M-1);
	dt3.resize(N-1);
	for(auto i=0; i<M-1; i++){
		dt2[i].resize(N-1);
		dt3[i].resize(N-1);
	}
	for(auto i=1; i<M; i++){
		for(auto j=1; j<N; j++){
			a = u[i+1][j] + u[i-1][j];
			b = v[i][j+1] + v[i][j-1];
			if(std::abs(a + 0.5*b) < dx) dt2[i-1][j-1] = dx/(a+0.5*b);
			else dt2[i-1][j-1] = 1.0;
			if(std::abs(0.5*a + b) < dy) dt3[i-1][j-1] = dy/(0.5*a+b);
			else dt3[i-1][j-1] = 1.0;
			Re1 = std::abs(a + 0.5*b) * dx / (1.0/Re);
			Re2 = std::abs(0.5*a + b) * dy / (1.0/Re);
		}
	}

	auto s1 = dt2.size();
	auto s2 = dt3.size();
	double w1 = algo.minval(s1, dt2);
	double w2 = algo.minval(s2, dt3);
	w1 = std::min(w1, w2);
	dt = std::min(dt1, w1);
}

void calculation::predict(const double Re, const double dt, const int M, const int N,\
          std::vector<std::vector<double>> &u,std::vector<std::vector<double>>&v,\
       	 std::vector<std::vector<double>>&us,std::vector<std::vector<double>>&vs){
	double dx = 1.0/M;
	double dy = 1.0/N;
	double ud = 0.0; // b.c.
	double ur = 0.0; // b.c.
	double uu = 1.0; // b.c.
	double ul = 0.0; // b.c.
	for(auto i=1; i<M; i++){
		for(auto j=0; j<N; j++){
			us[i][j] = u[i][j+1] + dt * (-(pow(0.5*(u[i+1][j+1]+u[i][j+1]),2)-\
						  pow(0.5*(u[i][j+1]+u[i-1][j+1]),2))/dx\
						  -(0.5*(u[i][j+2]+u[i][j+1])*0.5*(v[i+1][j+1]+v[i][j+1])-\
					     0.5*(u[i][j]+u[i][j+1])*0.5*(v[i][j]+v[i+1][j]))/dy+\
						  ((u[i+1][j+1]-2.0*u[i][j+1]+u[i-1][j+1])/pow(dx,2)+\
						  (u[i][j+2]-2.0*u[i][j+1]+u[i][j])/pow(dy,2))/Re);
		}
	}
	for(auto i=0; i<M; i++){
		for(auto j=1; j<N; j++){
			vs[i][j] = v[i+1][j] + dt * (-(0.5*(v[i+2][j]+v[i+1][j])*0.5*\
						  (u[i+1][j+1]+u[i+1][j])-(0.5*(v[i+1][j]+v[i][j])*0.5*\
						  (u[i][j+1]+u[i][j])))/dx-(pow(0.5*(v[i+1][j]+v[i+1][j+1]),2)-\
						  pow(0.5*(v[i+1][j]+v[i+1][j-1]),2))/dy+(((v[i+1][j+1]-2.0*\
						  v[i+1][j]+v[i+1][j-1])/pow(dy,2)+(v[i+2][j]-2.0*\
						  v[i+1][j]+v[i][j])/pow(dx,2))/Re));
		}
	}
	// set the boundary conditions for u* and v* which is set ti the values of the
	// u and v at the walls.
	for(auto i=0; i<N; i++) us[0][i] = u[0][i+1];
	for(auto i=0; i<N; i++) us[M][i] = u[M][i+1];
	for(auto i=0; i<M; i++) vs[i][0] = v[i+1][0];
	for(auto i=0; i<M; i++) vs[i][N] = v[i+1][N];
}
void calculation::poisson(const int M, const int N, double dt, \
				 const std::vector<std::vector<double>>&rhs,\
         std::vector<std::vector<double>>& phi,int &num_itr,const double vst_diff){
	alg algo;
	double dx = 1.0/M;
	double dy = 1.0/N;
	double limit = 1.0;
	num_itr = 0;
	int k = 0;
	std::vector<std::vector<double>> res;
	res.resize(M);
	for(auto i=0; i<M; i++) res[i].resize(N);
	for(auto i=0; i<M; i++){
		for(auto j=0; j<N; j++){
			res[i][j] = 0.0;
		}
	}
	while(1){
		for(auto i=1; i<M+1; i++){
			for(auto j=1; j<N+1; j++){
				phi[i][j] = (pow(dy,2)*(phi[i-1][j]+phi[i+1][j])+pow(dx,2)*(phi[i][j-1]\
								 +phi[i][j+1])-pow(dx,2)*pow(dy,2)*rhs[i-1][j-1])\
								/(2.0*(pow(dx,2)+pow(dy,2)));

				phi[i][0]   = phi[i][1];
				phi[i][N+1] = phi[i][N];
				phi[0][j]	= phi[1][j];
				phi[M+1][j] = phi[M][j];
			}
		}
		// convergence criterion is checked every 10 iteration
		if(k%10 == 0){
			//these loops go over x and y index to calculate the residuals
			for(auto i=1; i<M+1; i++){
				for(auto j=1; j<N+1; j++){
					res[i-1][j-1] = rhs[i-1][j-1] - \
					((phi[i+1][j]-2.0*phi[i][j]+phi[i-1][j])/pow(dx,2)+\
					(phi[i][j+1]-2.0*phi[i][j]+phi[i][j-1])/pow(dy,2));
				}
			}
			// limit is the maximum value of residual
			auto s1 = res.size();
			limit = algo.maxval(s1, res);
			num_itr += 1;
			double temp = std::max(0.8 * vst_diff/dt, pow(10.0,-7));
//			std::cout<<temp<<'\t'<<limit<<std::endl;
			if(limit < temp) break;
		}
		k += 1;
	}
}
void calculation::correct(const double dt, const int M, const int N,\
             std::vector<std::vector<double>>& us,\
				 std::vector<std::vector<double>>& vs,\
             std::vector<std::vector<double>>& un,\
				 std::vector<std::vector<double>>& vn,\
             std::vector<std::vector<double>>& phi){
	double dx = 1.0/M;
	double dy = 1.0/N;
	// correrction of u
	for(auto i=1; i<M; i++){
		for(auto j=0; j<N; j++){
			un[i][j] = us[i][j] - dt * (phi[i+1][j+1]-phi[i][j+1]) / dx;
		}
	}
	// correction of v	
	for(auto i=0; i<M; i++){
		for(auto j=1; j<N; j++){
			vn[i][j] = vs[i][j] - dt * (phi[i+1][j+1]-phi[i+1][j]) / dy;
		}
	}
	auto s1 = un[0].size();
	auto s2 = vn.size();
	for(auto i=0; i<s1; i++){
		un[0][i] = 0.0; 
		un[M][i] = 0.0;
	}
	for(auto i=0; i<s2; i++){
		vn[i][0] = 0.0;
		vn[i][N] = 0.0;
	}
}

void calculation::gs2d(const int M,const int N,const std::vector<std::vector<double>>& rhs,\
          std::vector<std::vector<double>>& phi, const int num_itr){
	double dx = 1.0/M;
	double dy = 1.0/N;
	for(auto k=0; k<num_itr; k++){
		for(auto i=1; i<M; i++){
			for(auto j=1; j<N; j++){
				phi[i][j] = (pow(dy,2)*(phi[i-1][j]+phi[i+1][j])+pow(dx,2)*\
								(phi[i][j-1]+phi[i][j+1])-pow(dx,2)*pow(dy,2)*rhs[i][j])\
								/(2.0*(pow(dx,2)+pow(dy,2)));
			}
		}
	}
}





