#include"calculation.h"
#include"utility.h"
#include"alg.h"
#include<omp.h>
#include<memory>
#include<iostream>
#include<limits>
#include<cfloat>
void calculation::calcTimeStep(const vector<vector<double>> &rho, \
										 const vector<vector<double>> &mu, \
										 const int M,const int N, \
										 const std::vector<std::vector<double>> &u,\
				             const std::vector<std::vector<double>> &v, double & dt){
	alg algo;
	double a, b, Re1, Re2;
	double dx = 1.0/M;
	double dy = 1.5/N;
	auto size1 = rho.size();
	auto size2 = rho[0].size();
	double max_mu_rho = -DBL_MAX;
	for(auto i=0; i<size1; i++){
		for(auto j=0; j<size2; j++){
			max_mu_rho = std::max(mu[i][j]/rho[i][j], max_mu_rho);
		}
	}
	double dt1 = pow(dx,2)/(4.0 * max_mu_rho);
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
			Re1 = std::abs(a + 0.5*b) * dx / (mu[i+3][j+3] / rho[i+3][j+3]);
			Re2 = std::abs(0.5*a + b) * dy / (mu[i+3][j+3] / rho[i+3][j+3]);
		}
	}
	auto s1 = dt2.size();
	auto s2 = dt3.size();
	double w1 = algo.minval(s1, dt2);
	double w2 = algo.minval(s2, dt3);
	w1 = std::min(w1, w2);
	dt = std::min(dt1, w1);
	dt *= 0.25;
//	std::cout<<"Delta t for Flow Solver = "<<dt<<std::endl;
}

void calculation::predict(const vector<vector<double>>& rho_n,\
								  const vector<vector<double>> &rho, \
								  const vector<vector<double>> &mu, const double dt, \
								  const int M, const int N, \
								  std::vector<std::vector<double>> &u, \
								  std::vector<std::vector<double>>&v,\
					       	  std::vector<std::vector<double>>&un, \
								  std::vector<std::vector<double>>&vn, \
								  std::vector<std::vector<double>>&poo, \
								  std::vector<std::vector<double>>&po, \
								  const std::vector<std::vector<double>> &kappa, \
								  const std::vector<std::vector<double>> &phi){
	std::unique_ptr<utility> putil(new utility);
	std::unique_ptr<alg> palg(new alg);
	double gx = 0.0;
	double gy = -1.81;
	double dx = 1.0/M;
	double dy = 1.5/N;
	double mux1, mux2, muy1, muy2;
	double mduu_xx, mduu_yy, mdvv_xy;
	double mdvv_yy, mdvv_xx, mduu_xy;
	double grad_phi, phi_m, kappa_m, x_surf, y_surf, deltamo_val, rho_m, rhon_m;
	double poo_m, po_m, u_m, v_m;
//	std::vector<double> L1((M-1)*N,0.0), D1((M-1)*N,0.0);
//	std::vector<double> U1((M-1)*N,0.0), Y1((M-1)*N,0.0), X1((M-1)*N,0.0);

//	std::vector<double> L2(M*(N-1),0.0), D2(M*(N-1),0.0);
//	std::vector<double> U2(M*(N-1),0.0), Y2(M*(N-1),0.0), X2(M*(N-1),0.0);

	std::vector<double> L1(M-1,0.0),D1(M-1,0.0),U1(M-1, 0.0);
	std::vector<double> Y1(M-1,0.0),X1(M-1,0.0);
	std::vector<double> L2(N,0.0), D2(N,0.0);
	std::vector<double> U2(N,0.0), Y2(N,0.0), X2(N,0.0);

	std::vector<std::vector<double>> u_ms, v_ms;
	double rhs1, rhs2;
	std::vector<std::vector<double>> res, rhs;
	res.resize(M-1); rhs.resize(M-1);
	for(auto i=0; i<M-1; i++){
		res[i].resize(N);
		rhs[i].resize(N);
	}
	u_ms.resize(M+1); v_ms.resize(M+2);

	for(auto i=0; i<M+1; i++) u_ms[i].resize(N+2);
	for(auto i=0; i<M+2; i++) v_ms[i].resize(N+1);
	int k = 0;
	auto s1 = res.size();
	for(auto j=0; j<N+2; j++){
		un[0][j] = 0.0;
		un[M][j]  = 0.0;
	}
	for(auto i=0; i<M+2; i++){
		vn[i][0] = 0.0;
		vn[i][N] = 0.0;
	}
	for(auto i=1; i<M; i++){
		for(auto j=1; j<N+1; j++){
			un[0][j] = 0.0;
			un[M][j] = 0.0;
			un[i][0] = -un[i][1];
			un[i][N+1] = -un[i][N];
			grad_phi = (phi[i+3][j+2] - phi[i+2][j+2]) / dx;
			kappa_m  = 0.5 * (kappa[i+1][j] + kappa[i][j]);
			phi_m    = 0.5 * (phi[i+3][j+2] + phi[i+2][j+2]);
			deltamo_val = putil->deltamo(phi_m);
			rho_m    = 0.5 * (rho[i+3][j+2] + rho[i+2][j+2]);

			x_surf = 0.05 * (kappa_m * deltamo_val * grad_phi);

			mux1 = 0.25*(mu[i+3][j+2]+mu[i+3][j+3]+mu[i+2][j+3]+mu[i+2][j+2]);

			mux2 = 0.25*(mu[i+3][j+2]+mu[i+3][j+1]+mu[i+2][j+1]+mu[i+2][j+2]);

			poo_m = 0.25*(poo[i+1][j+1]+poo[i][j+1]-poo[i+1][j-1]-poo[i][j-1]);

			po_m  = 0.25*(po[i+1][j+1]+po[i][j+1]-po[i+1][j-1]-po[i][j-1]);

			v_m   = 0.25*(v[i][j] + v[i+1][j] + v[i][j-1] + v[i+1][j-1]);

			rhs1 = (poo[i+1][j]-poo[i][j])/dx+poo_m/dx-2.0*((po[i+1][j]-\
					  po[i][j])/dx+po_m/dx)+x_surf+rho_m*u[i][j]/dt;

			rhs2=(mu[i+3][j+2]*(u[i+1][j]-u[i][j])+mu[i+2][j+2]*(u[i][j]-u[i-1][j]))/pow(dx,2)+\
			(mux1*(u[i][j+1]-u[i][j])-mux2*(u[i][j]-u[i][j-1]))/pow(dy,2) - \
			u[i][j]*rho_m*0.5*(u[i+1][j]-u[i-1][j])/dx-\
			v_m*(u[i][j+1]-u[i][j-1])*0.5*rho_m/dx;

			rhs[i-1][j-1] = rhs1 + rhs2;

			un[i][j]= (dt/rho_m)*(rhs[i-1][j-1]);
			un[0][j] = 0.0;
			un[M][j] = 0.0;
			un[i][0] = -un[i][1];
			un[i][N+1] = 0.0-un[i][N];
		}
	}
	for(auto j=1; j<N+1; j++) un[0][j] = 0.0;
	for(auto j=1; j<N+1; j++) un[M][j] = 0.0;
	for(auto i=0; i<M+1; i++){
		un[i][0] = -un[i][1];
		un[i][N+1] = -un[i][N];
	}
	res.clear();
	rhs.clear();
	res.resize(M); rhs.resize(M);
	for(auto i=0; i<M; i++){
		res[i].resize(N-1);
		rhs[i].resize(N-1);
	}
	for(auto i=1; i<M+1; i++){
		for(auto j=1; j<N; j++){
			vn[i][0] = 0.0;
			vn[i][N] = 0.0;
			grad_phi = (phi[i+2][j+3] - phi[i+2][j+2]) / dy;
			kappa_m  = 0.5 * (kappa[i][j+1] + kappa[i][j]);
			phi_m    = 0.5 * (phi[i+2][j+3] + phi[i+2][j+2]);
			deltamo_val = putil->deltamo(phi_m);
			rho_m    = 0.5 * (rho[i+2][j+3] + rho[i+2][j+2]);
			y_surf = 0.05 * (kappa_m * deltamo_val * grad_phi) / rho_m;

			muy1 = 0.25*(mu[i+3][j+2]+mu[i+3][j+3]+mu[i+2][j+3]+mu[i+2][j+2]);

			muy2 = 0.25*(mu[i+2][j+2]+mu[i+2][j+3]+mu[i+1][j+3]+mu[i+1][j+2]);

			poo_m = 0.25*(poo[i+1][j]+poo[i+1][j+1]-poo[i][j-1]-poo[i-1][j+1]);
			po_m  = 0.25*(po[i+1][j]+po[i+1][j+1]-po[i][j-1]-po[i-1][j+1]);
			u_m   = 0.25*(u[i][j] + u[i][j+1] + u[i-1][j] + u[i-1][j+1]);

 			rhs1 = (poo[i][j+1]-poo[i][j])/dy+poo_m/dy-2.0*((po[i][j+1]-\
					  po[i][j])/dy+po_m/dy)+y_surf+rho_m*v[i][j]/dt - rho_m*gy;

			rhs2=(mu[i+2][j+3]*(v[i][j+1]-v[i][j])-mu[i+2][j+2]*(v[i][j]-v[i][j-1]))/pow(dy,2)+\
			(muy1*(v[i+1][j]-v[i][j])-muy2*(v[i][j]-v[i][j-1]))/pow(dx,2)-\
			v[i][j]*rho_m*0.5*(v[i][j+1]-v[i][j-1])/dy-\
			u_m*(v[i+1][j]-v[i-1][j])*0.5*rho_m/dx;
			rhs[i-1][j-1] = rhs1 + rhs2;
 			vn[i][j]= (dt/rho_m)*(rhs[i-1][j-1]);
		}
	}
	for(auto i=1; i<M+1; i++) vn[i][0] = 0.0;
	for(auto i=1; i<M+1; i++) vn[i][N] = 0.0;
	//set neumann boundary condition
	for(auto j=0; j<N+1; j++){
		vn[0][j] = -vn[1][j]; 
		vn[M+1][j] = -vn[M][j];
	}
	v = vn;
	u = un;
}
while(1){
	for(auto i=1; i<M; i++){
		for(auto j=1; j<N+1; j++){
			un[0][j] = 0.0;
			un[M][j] = 0.0;
			un[i][0] = -un[i][1];
			un[i][N+1] = 0.01-un[i][N];
			grad_phi = (phi[i+3][j+2] - phi[i+2][j+2]) / dx;
			kappa_m  = 0.5 * (kappa[i+1][j] + kappa[i][j]);
			phi_m    = 0.5 * (phi[i+3][j+2] + phi[i+2][j+2]);
			deltamo_val = putil->deltamo(phi_m);
			rho_m    = 0.5 * (rho[i+3][j+2] + rho[i+2][j+2]);

			x_surf = 0.0;//0.05 * (kappa_m * deltamo_val * grad_phi);

			mux1 = 0.25*(mu[i+3][j+2]+mu[i+3][j+3]+mu[i+2][j+3]+mu[i+2][j+2]);

			mux2 = 0.25*(mu[i+3][j+2]+mu[i+3][j+1]+mu[i+2][j+1]+mu[i+2][j+2]);

			poo_m = 0.25*(poo[i+1][j+1]+poo[i][j+1]-poo[i+1][j-1]-poo[i][j-1]);

			po_m  = 0.25*(po[i+1][j+1]+po[i][j+1]-po[i+1][j-1]-po[i][j-1]);

			v_m   = 0.25*(v[i][j] + v[i+1][j] + v[i][j-1] + v[i+1][j-1]);

			rhs1 = (poo[i+1][j]-poo[i][j])/dx+poo_m/dx-2.0*((po[i+1][j]-\
					  po[i][j])/dx+po_m/dx)+x_surf+rho_m*u[i][j]/dt;

			rhs2 = (mu[i+3][j+2]*un[i+1][j]+mu[i+2][j+2]*un[i-1][j])/pow(dx,2)+\
					 (mux1*un[i][j+1]+mux2*un[i][j-1])/pow(dy,2)-\
					  u[i][j]*rho_m*0.5*(un[i+1][j]-un[i-1][j])/dx-v_m*0.5*\
					 (un[i][j+1]-un[i][j-1])*rho_m/dy;
			std::cout<<rhs1<<'\t'<<rhs2<<'\t'<<i<<'\t'<<j<<std::endl;

			rhs[i-1][j-1] = rhs1 + rhs2;

			un[i][j]= (1.0/(rho_m/dt+(mu[i+3][j+2]+mu[i+2][j+2])/pow(dx,2)+\
						 (mux1+mux2)/pow(dy,2)))*(rhs[i-1][j-1]);
			un[0][j] = 0.0;
			un[M][j] = 0.0;
			un[i][0] = -un[i][1];
			un[i][N+1] = 0.01-un[i][N];
		}
	}
	for(auto j=1; j<N+1; j++) un[0][j] = 0.0;
	for(auto j=1; j<N+1; j++) un[M][j] = 0.0;
	for(auto i=0; i<M+1; i++){
		un[i][0] = -un[i][1];
		un[i][N+1] = -un[i][N];
	}
	if(k%10 == 0){
		for(auto i=1; i<M; i++){
			for(auto j=1; j<N+1; j++){
				rhs1 = (poo[i+1][j]-poo[i][j])/dx+poo_m/dx-2.0*((po[i+1][j]-\
				po[i][j])/dx+po_m/dx)+x_surf+rho_m*u[i][j]/dt;
				res[i-1][j-1]=\
				rho_m*un[i][j]/dt - (mu[i+3][j+2]*(un[i+1][j]-un[i][j])-\
				mu[i+2][j+2]*(un[i][j]-un[i-1][j]))/pow(dx,2)-(mux1*(un[i][j+1]-un[i][j])-\
				mux2*(un[i][j]-un[i][j-1]))/pow(dy,2)+u[i][j]*rho_m*0.5*(un[i+1][j]-\
				un[i-1][j])/dy+0.5*v_m*rho_m*(un[i][j+1]-un[i][j-1])/dy - rhs1;
			}
		}
		limit = palg->maxval(s1,res);
		if(k%1000==0) std::cout<<"U velocity residual = "<<limit<<std::endl;
		if(limit<1.0e-7||k>200000) break;// || std::abs(limit-limit_old)<1.0e-15) break;
		limit_old = limit;
	}
	k+=1;
}
res.clear();
rhs.clear();
res.resize(M); rhs.resize(M);
for(auto i=0; i<M; i++){
	res[i].resize(N-1);
	rhs[i].resize(N-1);
}
k=1;
limit = 1.0;
limit_old = 0.0;
while(1){
	for(auto i=1; i<M+1; i++){
		for(auto j=1; j<N; j++){
			vn[i][0] = 0.0;
			vn[i][N] = 0.0;
			grad_phi = (phi[i+2][j+3] - phi[i+2][j+2]) / dy;
			kappa_m  = 0.5 * (kappa[i][j+1] + kappa[i][j]);
			phi_m    = 0.5 * (phi[i+2][j+3] + phi[i+2][j+2]);
			deltamo_val = putil->deltamo(phi_m);
			rho_m    = 0.5 * (rho[i+2][j+3] + rho[i+2][j+2]);
			y_surf = 0.0;//0.05 * (kappa_m * deltamo_val * grad_phi) / rho_m;

			muy1 = 0.25*(mu[i+3][j+2]+mu[i+3][j+3]+mu[i+2][j+3]+mu[i+2][j+2]);

			muy2 = 0.25*(mu[i+2][j+2]+mu[i+2][j+3]+mu[i+1][j+3]+mu[i+1][j+2]);

			poo_m = 0.25*(poo[i+1][j]+poo[i+1][j+1]-poo[i][j-1]-poo[i-1][j+1]);
			po_m  = 0.25*(po[i+1][j]+po[i+1][j+1]-po[i][j-1]-po[i-1][j+1]);
			u_m   = 0.25*(u[i][j] + u[i][j+1] + u[i-1][j] + u[i-1][j+1]);

 			rhs1 = (poo[i][j+1]-poo[i][j])/dy+poo_m/dy-2.0*((po[i][j+1]-\
					  po[i][j])/dy+po_m/dy)+y_surf+rho_m*v[i][j]/dt - rho_m*gy;

			rhs2 = (mu[i+2][j+3]*vn[i][j+1]+mu[i+2][j+2]*vn[i][j-1])/pow(dy,2)+\
					 (muy1*vn[i+1][j]+muy2*vn[i-1][j])/pow(dx,2)-\
					  v[i][j]*rho_m*0.5*(vn[i][j+1]-un[i][j-1])/dy-u_m*0.5*\
					 (vn[i+1][j]-un[i-1][j])*rho_m/dx;
			rhs[i-1][j-1] = rhs1 + rhs2;
 			vn[i][j]= (1.0/(rho_m/dt+(mu[i+2][j+3]+mu[i+2][j+2])/pow(dy,2)+\
						 (muy1+muy2)/pow(dx,2)))*(rhs[i-1][j-1]);
		}
	}
	for(auto i=1; i<M+1; i++) vn[i][0] = 0.0;
	for(auto i=1; i<M+1; i++) vn[i][N] = 0.0;
	//set neumann boundary condition
	for(auto j=0; j<N+1; j++){
		vn[0][j] = -vn[1][j]; 
		vn[M+1][j] = -vn[M][j];
	}
	if(k%10 == 0){
		for(auto i=1; i<M+1; i++){
			for(auto j=1; j<N; j++){
	 			rhs1 = (poo[i][j+1]-poo[i][j])/dy+poo_m/dy-2.0*((po[i][j+1]-\
						  po[i][j])/dy+po_m/dy)+y_surf+rho_m*v[i][j]/dt-rho_m*gy;
				res[i-1][j-1]=\
				rho_m*vn[i][j]/dt - (mu[i+2][j+3]*(vn[i][j+1]-vn[i][j])-\
				mu[i+2][j+2]*(vn[i][j]-vn[i][j-1]))/pow(dy,2)-(mux1*(vn[i+1][j]-vn[i][j])-\
				mux2*(vn[i][j]-un[i-1][j]))/pow(dx,2)+v[i][j]*rho_m*0.5*(vn[i][j+1]-\
				vn[i][j-1])/dy+0.5*v_m*rho_m*(vn[i+1][j]-un[i-1][j])/dx - rhs1;
			}
		}
		limit = palg->maxval(s1,res);
		if(k%1000==0) std::cout<<"V velocity residual = "<<limit<<std::endl;
		if(limit<1.0e-7 ||k>200000) break;// std::abs(limit-limit_old)<1.0e-15) break;
		limit_old = limit;
	}
	k+=1;
v = vn;
u = un;
}
double new_term;
// first step of ADI for u
	for(auto j=1; j<N+1; j++){
		for(auto i=1; i<M; i++){
			u_ms[0][j] = 0.0;
			u_ms[M][j] = 0.0;		
			grad_phi = (phi[i+3][j+2] - phi[i+2][j+2]) / dx;
			kappa_m  = 0.5 * (kappa[i+1][j] + kappa[i][j]);
			phi_m    = 0.5 * (phi[i+3][j+2] + phi[i+2][j+2]);
			deltamo_val = putil->deltamo(phi_m);
			rho_m    = 0.5 * (rho[i+3][j+2] + rho[i+2][j+2]);
			rhon_m   = 0.5 * (rho_n[i+3][j+2] + rho_n[i+2][j+2]);
			new_term=0.5*dt*(rhon_m-rho_m)+0.5*(0.5*rho_m*(u[i+1][j]-u[i-1][j])/dx+\
			rho_m*0.5*(v[i][j]+v[i+1][j]-v[i][j-1]-v[i+1][j-1])/dy);
			
			x_surf = 0.05 * (kappa_m * deltamo_val * grad_phi);

			mux1 = 0.25*(mu[i+3][j+2]+mu[i+3][j+3]+mu[i+2][j+3]+mu[i+2][j+2]);

			mux2 = 0.25*(mu[i+3][j+2]+mu[i+3][j+1]+mu[i+2][j+1]+mu[i+2][j+2]);

			poo_m = 0.25*(poo[i+1][j+1]+poo[i][j+1]-poo[i+1][j-1]-poo[i][j-1]);
			po_m  = 0.25*(po[i+1][j+1]+po[i][j+1]-po[i+1][j-1]-po[i][j-1]);
			v_m   = 0.25*(v[i][j] + v[i+1][j] + v[i][j-1] + v[i+1][j-1]);

			Y1[i-1] = (poo[i+1][j] - poo[i][j]) / dx + poo_m / dy- \
						  2.0*(po[i+1][j] - po[i][j]) / dx -2.0*(po_m / dx) - \
							v_m*0.5*(u[i][j+1]-u[i][j-1])*rho_m/dy+\
				(mux1*(u[i][j+1]-u[i][j])-mux2*(u[i][j]-u[i][j-1]))/pow(dy,2)\
				 + x_surf + rho_m*u[i][j]/(0.5*dt);

			L1[i-1] = -mu[i+2][j+2]/pow(dx,2)-u[i][j]*0.5*rho_m/dx;

			U1[i-1] = -mu[i+3][j+2]/pow(dx,2)+u[i][j]*0.5*rho_m/dx;

			D1[i-1] = rho_m/(0.5*dt)+(mu[i+3][j+2]+mu[i+2][j+2])/pow(dx,2)+new_term;
		}
		putil->tridiag(L1, D1, U1, Y1, X1);
		for(auto i=1; i<M; i++){
			u_ms[i][j] = X1[i-1];
		}
//		u_ms[0][j] = 0.0;
//		u_ms[M][j] = 0.0;
		for(auto i=0; i<M+1; i++){
			u_ms[i][0] = -u_ms[i][1];
			u_ms[i][N+1] = - u_ms[i][N];
		}
	}

	//second step of the ADI for u
	for(auto i=1; i<M; i++){
		for(auto j=1; j<N+1; j++){
			un[0][j] = 0.0;
			un[M][j] = 0.0;

			grad_phi = (phi[i+3][j+2] - phi[i+2][j+2]) / dx;
			kappa_m  = 0.5 * (kappa[i+1][j] + kappa[i][j]);
			phi_m    = 0.5 * (phi[i+3][j+2] + phi[i+2][j+2]);
			deltamo_val = putil->deltamo(phi_m);
			rho_m    = 0.5 * (rho[i+3][j+2] + rho[i+2][j+2]);
			rhon_m   = 0.5 * (rho_n[i+3][j+2] + rho_n[i+2][j+2]);
			x_surf = 0.05 * (kappa_m * deltamo_val * grad_phi);

			new_term=0.5*dt*(rhon_m-rho_m)+0.5*(0.5*rho_m*(u[i+1][j]-u[i-1][j])/dx+\
			rho_m*0.5*(v[i][j]+v[i+1][j]-v[i][j-1]-v[i+1][j-1])/dy);

			mux1 = 0.25*(mu[i+3][j+2]+mu[i+3][j+3]+mu[i+2][j+3]+mu[i+2][j+2]);

			mux2 = 0.25*(mu[i+3][j+2]+mu[i+3][j+1]+mu[i+2][j+1]+mu[i+2][j+2]);

			poo_m = 0.25*(poo[i+1][j+1]+poo[i][j+1]-poo[i+1][j-1]-poo[i][j-1]);
			po_m  = 0.25*(po[i+1][j+1]+po[i][j+1]-po[i+1][j-1]-po[i][j-1]);
			v_m   = 0.25*(v[i][j] + v[i+1][j] + v[i][j-1] + v[i+1][j-1]);

			Y2[j-1] = (poo[i+1][j]-poo[i][j])/dx+poo_m / dx - \
						  2.0*(po[i+1][j] - po[i][j]) / dx -2.0*(po_m / dx) + \
						  (mu[i+3][j+2]*(u_ms[i+1][j]-u_ms[i][j])-mu[i+2][j+2]*\
						  (u_ms[i][j]-u_ms[i-1][j]))/pow(dx,2)-\
						  (u[i][j]*0.5*(u_ms[i+1][j]-u_ms[i-1][j])*rho_m/dx)\
						  + x_surf + rho_m*u_ms[i][j]/(0.5*dt);
			if(j==1){
				L2[j-1] = 0.0;
				D2[j-1] = -(-mux2/pow(dy,2)-v_m*0.5*rho_m/dy)+\
							 rho_m/(0.5*dt)+((mux1+mux2)/pow(dy,2))+new_term;
			}
			else if(j==N){
				U2[j-1] = 0.0;
				D2[j-1] = rho_m/(0.5*dt)+((mux1+mux2)/pow(dy,2))+new_term\
							 -(-mux1/pow(dy,2)+v_m*0.5*rho_m/dy);
			}
			else{
				L2[j-1] = (-mux2/pow(dy,2)-v_m*0.5*rho_m/dy);
				U2[j-1] = (-mux1/pow(dy,2)+v_m*0.5*rho_m/dy);
				D2[j-1] = rho_m/(0.5*dt)+((mux1+mux2)/pow(dy,2))+new_term;
			}
//			un[0][j] = 0.0;
//			un[M][j] = 0.0;
		}
		putil->tridiag(L2, D2, U2, Y2, X2);
		for(auto j=1; j<N+1; j++){
			un[i][j] = X2[j-1];
		}
		un[i][0] = -un[i][1];
		un[i][N+1] = -un[i][N];
	}
//first step of the ADI for v velocity
L1.clear(); D1.clear(); U1.clear(); Y1.clear(); X1.clear();
L2.clear(); D2.clear(); U2.clear(); Y2.clear(); X2.clear();

L1.resize(M); D1.resize(M); U1.resize(M);
Y1.resize(M); X1.resize(M);

L2.resize(N-1); D2.resize(N-1); U2.resize(N-1);
Y2.resize(N-1); X2.resize(N-1);

	for(auto j=1; j<N; j++){
		for(auto i=1; i<M+1; i++){
			v_ms[i][0] = 0.0;
			v_ms[i][N] = 0.0;
			grad_phi = (phi[i+2][j+3] - phi[i+2][j+2]) / dy;
			kappa_m  = 0.5 * (kappa[i][j+1] + kappa[i][j]);
			phi_m    = 0.5 * (phi[i+2][j+3] + phi[i+2][j+2]);
			deltamo_val = putil->deltamo(phi_m);
			rho_m    = 0.5 * (rho[i+2][j+3] + rho[i+2][j+2]);

			rhon_m   = 0.5 * (rho_n[i+2][j+2] + rho_n[i+2][j+3]);

			new_term=0.5*dt*(rhon_m-rho_m)+0.5*(0.5*rho_m*(v[i][j+1]-v[i][j-1])/dy+\
			rho_m*0.5*(u[i][j+1]+u[i-1][j+1]-u[i][j]-u[i-1][j])/dx);

			y_surf = 0.05 * (kappa_m * deltamo_val * grad_phi);


			muy1 = 0.25*(mu[i+3][j+2]+mu[i+3][j+3]+mu[i+2][j+3]+mu[i+2][j+2]);

			muy2 = 0.25*(mu[i+2][j+2]+mu[i+2][j+3]+mu[i+1][j+3]+mu[i+1][j+2]);

			poo_m = 0.25*(poo[i+1][j]+poo[i+1][j+1]-poo[i][j-1]-poo[i-1][j+1]);
			po_m  = 0.25*(po[i+1][j]+po[i+1][j+1]-po[i][j-1]-po[i-1][j+1]);
			u_m   = 0.25*(u[i][j] + u[i][j+1] + u[i-1][j] + u[i-1][j+1]);
			
			Y2[i-1] = ((mu[i+2][j+3]*(v[i][j+1]-v[i][j])-\
			mu[i+2][j+2]*(v[i][j]-v[i][j-1]))/pow(dy,2)\
			-v[i][j]*0.5*(v[i][j+1]-v[i][j-1])*rho_m/dy+(poo[i][j+1]-poo[i][j])/dy+ \
			poo_m/dy-2.0*(po[i][j+1]-po[i][j])/dy-2.0*po_m/dy)-rho_m*gy\
			+y_surf+rho_m*v[i][j]/(0.5*dt);
			if(i==1){
				L1[i-1] = 0.0;
				D1[i-1] = (muy1+muy2)/pow(dx,2) + rho_m/(0.5*dt)+new_term\
							 -(-muy2/pow(dx,2)-0.5*u_m*rho_m/dx);
			}
			else if(i==M){
				U1[i-1] = 0.0;
				D1[i-1] = ((muy1+muy2)/pow(dx,2)) + rho_m/(0.5*dt)+new_term\
							 -(-muy1/pow(dx,2)+u_m*0.5*rho_m/dx);
			}
			else{
				L1[i-1] = (-muy2/pow(dx,2)-0.5*u_m*rho_m/dx);
				U1[i-1] = (-muy1/pow(dx,2)+u_m*0.5*rho_m/dx);
				D1[i-1] = ((muy1+muy2)/pow(dx,2)) + rho_m/(0.5*dt)+new_term;
			}
//			v_ms[i][0] = 0.0;
//			v_ms[i][N] = 0.0;
		}
		putil->tridiag(L1, D1, U1, Y1, X1);
		for(auto i=1; i<M+1; i++){
			v_ms[i][j] = X1[i-1];
		}
		//set neumann boundary condition
		for(auto j=0; j<N+1; j++){
			v_ms[0][j] = -v_ms[1][j];
			v_ms[M+1][j] = -v_ms[M][j];
		}
	}

//second step of ADI
	for(auto i=1; i<M+1; i++){
		for(auto j=1; j<N; j++){
			vn[i][0] = 0.0;
			vn[i][N] = 0.0;
			grad_phi = (phi[i+2][j+3] - phi[i+2][j+2]) / dy;
			kappa_m  = 0.5 * (kappa[i][j+1] + kappa[i][j]);
			phi_m    = 0.5 * (phi[i+2][j+3] + phi[i+2][j+2]);
			deltamo_val = putil->deltamo(phi_m);
			rho_m    = 0.5 * (rho[i+2][j+3] + rho[i+2][j+2]);
			y_surf = 0.05 * (kappa_m * deltamo_val * grad_phi) / rho_m;

			rhon_m   = 0.5 * (rho_n[i+2][j+2] + rho_n[i+2][j+3]);

			new_term=0.5*dt*(rhon_m-rho_m)+0.5*(0.5*rho_m*(v[i][j+1]-v[i][j-1])/dy+\
			rho_m*0.5*(u[i][j+1]+u[i-1][j+1]-u[i][j]-u[i-1][j])/dx);

			muy1 = 0.25*(mu[i+3][j+2]+mu[i+3][j+3]+mu[i+2][j+3]+mu[i+2][j+2]);

			muy2 = 0.25*(mu[i+2][j+2]+mu[i+2][j+3]+mu[i+1][j+3]+mu[i+1][j+2]);

			poo_m = 0.25*(poo[i+1][j]+poo[i+1][j+1]-poo[i][j-1]-poo[i-1][j+1]);
			po_m  = 0.25*(po[i+1][j]+po[i+1][j+1]-po[i][j-1]-po[i-1][j+1]);
			u_m   = 0.25*(u[i][j] + u[i][j+1] + u[i-1][j] + u[i-1][j+1]);

			Y2[j-1] = (muy1*(v_ms[i+1][j]-v_ms[i][j])-\
			muy2*(v_ms[i][j]-v_ms[i-1][j]))/pow(dx,2)-0.5*rho_m*u_m*\
			(v_ms[i+1][j]-v_ms[i-1][j])/dx+(poo[i][j+1]-poo[i][j])/dy+\
			poo_m/dy-2.0*(po[i][j+1]-po[i][j])/dy-2.0*po_m/dy
			+ y_surf+rho_m*v_ms[i][j]/(0.5*dt)-rho_m*gy;

			L2[j-1] = (-mu[i+2][j+2]/pow(dy,2)-0.5*v[i][j]*rho_m/dy);

			U2[j-1] = (-mu[i+2][j+3]/pow(dy,2)+0.5*v[i][j]*rho_m/dy);

			D2[j-1] = ((mu[i+2][j+3]+mu[i+2][j+2])/pow(dy,2))+rho_m/(0.5*dt)+new_term;
		}
		putil->tridiag(L2, D2, U2, Y2, X2);
		for(auto j=1; j<N; j++){
			vn[i][j] = X2[j-1];
		}
//		vn[i][0] = 0.0;
//		vn[i][N] = 0.0;
		for(auto j=0; j<N+1; j++){
			vn[0][j] = -vn[1][j]; 
			vn[M+1][j] = -vn[M][j];
		}
	}
	u = un;
	v = vn;
void calculation::poisson(const double rho, const int M, \
								  const int N,double dtm,const vector<vector<double>>&u,\
         					  std::vector<std::vector<double>>& v, \
								  std::vector<std::vector<double>>& po,\
								  std::vector<std::vector<double>>& pn){
	alg algo;
	std::unique_ptr<utility> putil(new utility);
	double dx = 1.0/M;
	double dy = 1.5/N;
	double limit = 1.0;
	double limit_old = 0.0;
	double limit_diff;
	auto num_itr = 0;
	int k = 0;
	double rhs1, rhs2;
	std::vector<std::vector<double>> res, rhs;
	res.resize(M); rhs.resize(M);
	for(auto i=0; i<M; i++) res[i].resize(N);
	for(auto i=0; i<M; i++) rhs[i].resize(N);
	for(auto i=0; i<M; i++){
		for(auto j=0; j<N; j++){
			res[i][j] = 0.0;
			rhs[i][j] = 0.0;
		}
	}
	std::vector<std::vector<double>> p_ms;
	p_ms.resize(M+2);
	for(auto i=0; i<M+2; i++) p_ms[i].resize(N+2);

	std::vector<double> L1(M,0.0), D1(M,0.0);
	std::vector<double> U1(M,0.0), Y1(M,0.0), X1(M,0.0);

	std::vector<double> L2(N,0.0), D2(N,0.0);
	std::vector<double> U2(N,0.0), Y2(N,0.0), X2(N,0.0);

	auto s1 = res.size();
	auto ss1 = pn.size();
	auto ss2 = pn[0].size();
	
	for(auto i=0; i<ss1; i++){
		for(auto j=0; j<ss2; j++) pn[i][j] = 0.0;
	}
	std::vector<std::vector<double>> p_i = po;
auto dt = dtm/1.0;
for(auto manual=0; manual<1;  manual++){
//first step of ADI in x direction
	for(auto j=1; j<N+1; j++){
		for(auto i=1; i<M+1; i++){
			rhs1 = rho*((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy) / dt;
			rhs2 = (po[i+1][j]+po[i-1][j]+po[i][j+1]+po[i][j-1]-4.0*po[i][j])/pow(dx,2);
			Y1[i-1] = (rhs1+rhs2-(p_i[i][j+1]-2.0*p_i[i][j]+p_i[i][j-1])/\
									 pow(dy,2))*pow(dx,2);

			if(i==1){
				D1[i-1] = -1.0;
				L1[i-1] = 0.0;
			}
			else if(i==M){
				D1[i-1] = -1.0;
				U1[i-1] = 0.0;
			}
			else{
				D1[i-1] = -2.0;
				L1[i-1] = 1.0;
				U1[i-1] = 1.0;
			}
		}
		putil->tridiag(L1, D1, U1, Y1, X1);
		for(auto i=1; i<M+1; i++){
			p_ms[i][j] = X1[i-1];
		}
		p_ms[0][j] = p_ms[1][j];
		p_ms[M+1][j] = p_ms[M][j];
	}
//second step of ADI in y direction
	for(auto i=1; i<M+1; i++){
		for(auto j=1; j<N+1; j++){
			rhs1 = rho*((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy) / dt;
			rhs2 = (po[i+1][j]+po[i-1][j]+po[i][j+1]+po[i][j-1]-4.0*po[i][j])/pow(dx,2);
			Y2[j-1] = (rhs1+rhs2-(p_ms[i+1][j]-2.0*p_ms[i][j]+p_ms[i-1][j])/\
									 pow(dx,2))*pow(dy,2);
			if(j==1){
				D2[j-1] = -1.0;
				L2[j-1] = 0.0;
			}
			else if(j==N){
				D2[j-1] = -1.0;
				U2[j-1] = 0.0;
			}
			else{
				D2[j-1] = -2.0;
				L2[j-1] = 1.0;
				U2[j-1] = 1.0;
			}
		}
		putil->tridiag(L2, D2, U2, Y2, X2);
		for(auto j=1; j<N+1; j++){
			pn[i][j] = X2[j-1];
		}
		pn[i][0] = pn[i][1];
		pn[i][N+1] = pn[i][N];
	}
	p_i = pn;
}
std::cout<<algo.maxval(pn.size(), pn)<<"******"<<std::endl;
}
	while(1){
		for(auto i=1; i<M+1; i++){
			for(auto j=1; j<N+1; j++){
				rhs1 = (po[i+1][j]-2.0*po[i][j]+po[i-1][j])/pow(dx,2)+\
						 (po[i][j+1]-2.0*po[i][j]+po[i][j-1])/pow(dy,2);
				rhs2 = rho*((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy)/dt;
				rhs[i-1][j-1]  = rhs1 + rhs2;
				pn[i][j] = (pow(dy,2)*(pn[i+1][j]+pn[i-1][j])+\
								pow(dx,2)*(pn[i][j+1]+pn[i][j-1])-rhs[i-1][j-1]*pow(dx*dy,2))/ \
							  (2.0*(pow(dx,2)+pow(dy,2)));
				
				pn[i][0  ] = pn[i][1];
				pn[i][N+1] = pn[i][N];
				pn[0  ][j] = pn[1][j];
				pn[M+1][j] = pn[M][j];
			}
		}
		// convergence criterion is checked every 10 iteration
		if(k%10 == 0){
			//these loops go over x and y index to calculate the residuals
			for(auto i=1; i<M+1; i++){
				for(auto j=1; j<N+1; j++){

					res[i-1][j-1]=rhs[i-1][j-1]-((pn[i+1][j]-2.0*pn[i][j]+pn[i-1][j])/pow(dx,2)+\
									  (pn[i][j+1]-2.0*pn[i][j]+pn[i][j-1])/pow(dy,2));
				}
			}
//			limit is the maximum value of residual
			limit = algo.maxval(s1, res);
			num_itr += 1;
			double temp = 10.0e-15;//std::max(0.8 * vst_diff/dt, pow(10.0,-15));
			if(k%500==0) std::cout<<"pressure poisson residual = "<<limit<<std::endl;
			if(limit < temp || std::abs(limit-limit_old)<1.0e-17) break;
			limit_old = limit;
		}
		k += 1;
	}
std::cout<<algo.maxval(pn.size(), pn)<<"******"<<std::endl;

void calculation::correct(const vector<vector<double>> &rho,const double dt,\
								 const int M, const int N,\
				             std::vector<std::vector<double>>& us,\
								 std::vector<std::vector<double>>& vs,\
   				          std::vector<std::vector<double>>& un,\
								 std::vector<std::vector<double>>& vn,\
				             std::vector<std::vector<double>>& phi){
	double dx = 1.0/M;
	double dy = 1.5/N;
	double rhoa, rhoc;
	// correrction of u
	for(auto i=1; i<M; i++){
		for(auto j=0; j<N; j++){
			rhoa = 0.5*(rho[i+4][j+3] + rho[i+3][j+3]);
//			un[i][j] = us[i][j] - dt * (phi[i+1][j+1]-phi[i][j+1]) / (rho[i][j]*dx);
			un[i][j] = us[i][j] - dt * (phi[i+1][j+1]-phi[i][j+1])/(rhoa*dx);
		}
	}
	// correction of v
	for(auto i=0; i<M; i++){
		for(auto j=1; j<N; j++){
			rhoc = 0.5*(rho[i+3][j+3] + rho[i+3][j+4]);
//			vn[i][j] = vs[i][j] - dt * (phi[i+1][j+1]-phi[i+1][j]) / (rho[i][j]*dy);
			vn[i][j] = vs[i][j]-dt*(phi[i+1][j+1]-phi[i+1][j])/(rhoc*dy);
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





