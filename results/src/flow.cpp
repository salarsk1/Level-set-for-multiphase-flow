#include"flow.h"
#include"calculation.h"
#include"alg.h"
#include<memory>
void flow::flow_solver(int M, int N, double t_final,double Re, vector<vector<double>>&uin,vector<vector<double>>&vin){
	std::unique_ptr<calculation> pcal(new calculation);
	double dt;
	alg algo;
	int num_itr;
	std::vector<std::vector<double>> u, v, us, vs, un, vn, rhs, phi;
	std::vector<std::vector<double>> uo, vo, u_steady,  v_steady, omega, psi;
	double dx = 1.0 / M;
	double dy = 1.0 / N;
	u.resize(M+1);
	un.resize(M+1);
	us.resize(M+1);
	uo.resize(M+1);
	u_steady.resize(M+1);
	omega.resize(M+1);
	psi.resize(M+1);
	for(auto i=0; i<M+1; i++){
		u[i].resize(N+2);
		us[i].resize(N);
		un[i].resize(N);
		uo[i].resize(N);
		u_steady[i].resize(N);
		omega[i].resize(N+1);
		psi[i].resize(N+1);
	}
	v.resize(M+2);
	phi.resize(M+2);
	for(auto i=0; i<M+2; i++){
		v[i].resize(N+1);
		phi[i].resize(N+2);
	}
	vs.resize(M);
	vn.resize(M);
	rhs.resize(M);
	vo.resize(M);
	v_steady.resize(M);
	for(auto i=0; i<M; i++){
		vs[i].resize(N+1);
		vn[i].resize(N+1);
		rhs[i].resize(N);
		vo[i].resize(N+1);
		v_steady[i].resize(N+1);
	}
	for(auto i=0; i<M+2; i++){
		for(auto j=0; j<N+2; j++){
			phi[i][j] = 0.0;
		}
	}
	for(auto i=0; i<M+1; i++){
		for(auto j=0; j<N; j++){
			uo[i][j] = 0.0;
		}
	}
	for(auto i=0; i<M; i++){
		for(auto j=0; j<N+1; j++){
			vo[i][j] = 0.0;
		}
	}
	for(auto i=0; i<M+1; i++){
		for(auto j=0; j<N; j++){
			un[i][j] = 0.0;
		}
	}
	for(auto i=0; i<M; i++){
		for(auto j=0; j<N+1; j++){
			vn[i][j] = 0.0;
		}
	}
	for(auto i=0; i<M+1; i++){
		for(auto j=0; j<N+2; j++){
			u[i][j] = 0.0;
		}
	}
	for(auto i=0; i<M+2; i++){
		for(auto j=0; j<N+1; j++){
			v[i][j] = 0.0;
		}
	}
	double vst_diff = 0.0;
	double vst_diff_old = 1.0;
	double ud = 0.0;
	double ur = 0.0;
	double uu = 1.0;
	double ul = 0.0;
	for(auto i=0; i<M+1; i++){
		u[i][N+1] = 2.0 * uu - u[i][N];
		u[i][0]   = 2.0 * ud - u[i][1];
	}
	for(auto j=0; j<N+1; j++){
		v[0][j] = 2.0 * ul - v[1][j];
		v[M+1][j] = 2.0 * ur -v[M][j];
	}
//	std::ifstream input;
//	input.open("u_velocity.dat");
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
	while(time < t_final){
		pcal->calcTimeStep(Re, M, N, u, v, dt);
		time += dt;
		if(time > t_final){
			time -= dt;
			dt = t_final - time;
			time += dt;
		}
		pcal->predict(Re, dt, M, N, u, v, us, vs);
		std::cout<<dt<<'\t'<<time<<std::endl;
		for(auto i=0; i<M; i++){
			for(auto j=0; j<N; j++){
				rhs[i][j] = ((us[i+1][j]-us[i][j])/dx+(vs[i][j+1]-vs[i][j])/dy)/dt;
			}
		}
		pcal->poisson(M, N, dt, rhs, phi, num_itr, pow(10.0,-6));
		pcal->correct(dt, M, N, us, vs, un, vn, phi);
		for(auto i=0; i<M+1; i++){
			for(auto j=1; j<N+1; j++){
				u[i][j] = un[i][j-1];
			}
		}
		for(auto i=1; i<M+1; i++){
			for(auto j=0; j<N+1; j++){
				v[i][j] = vn[i-1][j];
			}
		}

		//Impose Boundary Conditions
		for(auto i=0; i<M+1; i++){
			u[i][N+1] = 2.0 * uu - u[i][N];
			u[i][0]   = 2.0 * ud - u[i][1];
		}
		for(auto j=0; j<N+1; j++){
			v[0][j] = 2.0 * ul - v[1][j];
			v[M+1][j] = 2.0 * ur - v[M][j];
		}
	}
	for(auto i=0; i<M+1; i++){
		for(auto j=0; j<N+2; j++) uin[i][j] = u[i][j];
	}
	for(auto i=0; i<M+2; i++){
		for(auto j=0; j<N+1; j++) vin[i][j] = v[i][j];
	}
}
