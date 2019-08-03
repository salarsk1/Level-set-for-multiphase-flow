#include"reinitialization.h"
#include<omp.h>
#include<array>
#include<cmath>
#include<algorithm>
#include"weno.h"
#include<memory>
#include<cfloat>
#include<iostream>
vector<vector<double>> reinitialize::sign_func(const double dx, const double dy, \
							  const vector<vector<double>> & phi){
	double rdx = 1.0/dx;
	double rdy = 1.0/dy;
	std::array<double, 2> grad_phi;
	auto size1 = phi.size();
	auto size2 = phi[0].size();
	vector<vector<double>> ret;
	ret.resize(size1-6);
	double d_min = std::min(dx, dy);
	for(auto i=0; i<size1-6; i++) ret[i].resize(size2-6);
	for(auto i=3; i<size1-3; i++){
		for(auto j=3; j<size2-3; j++){
			grad_phi[0] = 	0.5*(phi[i+1][j] - phi[i-1][j]) * rdx;
			grad_phi[1] = 	0.5*(phi[i][j+1] - phi[i][j-1]) * rdy;
			ret[i-3][j-3] = (phi[i][j] / sqrt(pow(phi[i][j],2)+(pow(grad_phi[0],2)+ \
								 pow(grad_phi[1],2))*pow(d_min,2))>=0.0) ? 1.0 : -1.0;
		}
	}
	return ret;
}

void reinitialize::pde_reinit(double dt, const double dx, const double dy, \
										vector<vector<double>> &phi){
	std::cout<<"Starting reinitialization"<<std::endl;
	std::unique_ptr<weno> pweno(new weno);
	bool is_done;
	auto size1 = phi.size();
	auto size2 = phi[0].size();
	std::vector<std::vector<double>> S;
	S = sign_func(dx, dy, phi);
	double rdx = 1.0/dx;
	double rdy = 1.0/dy;
   const double r23  = 2.0/3.0;
   const double r112 = 1.0/12.0;

   std::vector<std::vector<double>> alphaRK;
   alphaRK.resize(3);
   for(auto i=0; i<3; i++) alphaRK[i].resize(3);
   alphaRK[0][0]=dt;        alphaRK[1][0]=0.0;       alphaRK[2][0]=0.0;
   alphaRK[0][1]=-0.75*dt;  alphaRK[1][1]=0.25*dt;   alphaRK[2][1]=0.0;
   alphaRK[0][2]=-r112*dt;  alphaRK[1][2]=-r112*dt;  alphaRK[2][2]=r23*dt;

	vector<double> d_phix(6,0.0);
	vector<double> d_phiy(6,0.0);
	vector<double> grad_phi_p(2,0.0);
	vector<double> grad_phi_m(2,0.0);
	double v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15;
	double*** v_dot_grad;
   v_dot_grad =new double**[3];
   for(auto i=0; i<3; i++){
      v_dot_grad[i] = new double*[size1-6];
   }
   for(auto i=0; i<3; i++){
      for(auto j=0; j<size1-6; j++){
         v_dot_grad[i][j] = new double[size2-6];
      }
   }

	double base;
   for(auto j=3; j<size2-3; j++) phi[0][j] = phi[5][j];
   for(auto j=3; j<size2-3; j++) phi[1][j] = phi[4][j];
   for(auto j=3; j<size2-3; j++) phi[2][j] = phi[3][j];
   for(auto j=3; j<size2-3; j++) phi[size1-1][j] = phi[size1-6][j];
   for(auto j=3; j<size2-3; j++) phi[size1-2][j] = phi[size1-5][j];
   for(auto j=3; j<size2-3; j++) phi[size1-3][j] = phi[size1-4][j];

   for(auto i=3; i<size1-3; i++) phi[i][0] = phi[i][5];
   for(auto i=3; i<size1-3; i++) phi[i][1] = phi[i][4];
   for(auto i=3; i<size1-3; i++) phi[i][2] = phi[i][3];
   for(auto i=3; i<size1-3; i++) phi[i][size2-1] = phi[i][size2-6];
   for(auto i=3; i<size1-3; i++) phi[i][size2-2] = phi[i][size2-5];
   for(auto i=3; i<size1-3; i++) phi[i][size2-3] = phi[i][size2-4];
	for(auto iter=0; iter<=10; iter++){
		// Time integeration
	   for(auto irk=0; irk<3; irk++){
   	   for(auto i=3; i<size1-3; i++){
	         for(auto j=3; j<size2-3; j++){
		      	d_phix[0] = (phi[i-2][j] - phi[i-3][j]) * rdx;
            	d_phix[1] = (phi[i-1][j] - phi[i-2][j]) * rdx;
            	d_phix[2] = (phi[i  ][j] - phi[i-1][j]) * rdx;
            	d_phix[3] = (phi[i+1][j] - phi[i  ][j]) * rdx;
            	d_phix[4] = (phi[i+2][j] - phi[i+1][j]) * rdx;
            	d_phix[5] = (phi[i+3][j] - phi[i+2][j]) * rdx;

            	d_phiy[0] = (phi[i][j-2] - phi[i][j-3]) * rdy;
            	d_phiy[1] = (phi[i][j-1] - phi[i][j-2]) * rdy;
	            d_phiy[2] = (phi[i][j  ] - phi[i][j-1]) * rdy;
	            d_phiy[3] = (phi[i][j+1] - phi[i][j  ]) * rdy;
	            d_phiy[4] = (phi[i][j+2] - phi[i][j+1]) * rdy;
	            d_phiy[5] = (phi[i][j+3] - phi[i][j+2]) * rdy;

					grad_phi_p[0] = \
	            r112*(-d_phix[1]+7.0*d_phix[2]+7.0*d_phix[3]-d_phix[4]) + \
	            pweno->phi_weno_5th(d_phix[5]-d_phix[4],d_phix[4]-d_phix[3], \
	                                d_phix[3]-d_phix[2],d_phix[2]-d_phix[1]);

	            grad_phi_m[0] = \
	            r112*(-d_phix[1]+7.0*d_phix[2]+7.0*d_phix[3]-d_phix[4]) - \
	            pweno->phi_weno_5th(d_phix[1]-d_phix[0],d_phix[2]-d_phix[1], \
	                                d_phix[3]-d_phix[2],d_phix[4]-d_phix[3]);

	            grad_phi_p[1] = \
	            r112*(-d_phiy[1]+7.0*d_phiy[2]+7.0*d_phiy[3]-d_phiy[4]) + \
	            pweno->phi_weno_5th(d_phiy[5]-d_phiy[4],d_phiy[4]-d_phiy[3], \
	                                d_phiy[3]-d_phiy[2],d_phiy[2]-d_phiy[1]);

   	         grad_phi_m[1] = \
	            r112*(-d_phiy[1]+7.0*d_phiy[2]+7.0*d_phiy[3]-d_phiy[4]) - \
	            pweno->phi_weno_5th(d_phiy[1]-d_phiy[0],d_phiy[2]-d_phiy[1], \
	                                d_phiy[3]-d_phiy[2],d_phiy[4]-d_phiy[3]);

	            double v1 = std::max(S[i-3][j-3], 0.0);
	            double v2 = std::max(grad_phi_m[0], 0.0);
	            double v3 = std::min(grad_phi_p[0], 0.0);
	            double v4 = std::max(pow(v2,2), pow(v3,2));
	            double v5 = std::max(grad_phi_m[1], 0.0);
	            double v6 = std::min(grad_phi_p[1], 0.0);
	            double v7 = std::max(pow(v5,2), pow(v6,2));

	            double v11 = std::min(S[i-3][j-3], 0.0);
	            double v12 = std::min(grad_phi_m[0], 0.0);
	            double v13 = std::max(grad_phi_p[0], 0.0);
	            double v14 = std::max(pow(v12, 2), pow(v13, 2));
	            double v15 = std::min(grad_phi_m[1], 0.0);
	            double v16 = std::max(grad_phi_p[1], 0.0);
	            double v17 = std::max(pow(v15, 2), pow(v16, 2));
	            v_dot_grad[irk][i-3][j-3] = -(v1*(sqrt(v4+v7)-1.0)+\
															v11*(sqrt(v14+v17)-1.0));
	            for(auto irki=0; irki<=irk; irki++){
	               phi[i][j] += alphaRK[irki][irk] * v_dot_grad[irki][i-3][j-3];
	            }
				}
			}
	      for(auto j=3; j<size2-3; j++) phi[0      ][j] = phi[5      ][j];
			for(auto j=3; j<size2-3; j++) phi[1      ][j] = phi[4      ][j];
			for(auto j=3; j<size2-3; j++) phi[2      ][j] = phi[3      ][j];
			for(auto j=3; j<size2-3; j++) phi[size1-1][j] = phi[size1-6][j];
			for(auto j=3; j<size2-3; j++) phi[size1-2][j] = phi[size1-5][j];
			for(auto j=3; j<size2-3; j++) phi[size1-3][j] = phi[size1-4][j];
			for(auto i=3; i<size1-3; i++) phi[i][0      ] = phi[i][5      ];
			for(auto i=3; i<size1-3; i++) phi[i][1      ] = phi[i][4      ];
			for(auto i=3; i<size1-3; i++) phi[i][2      ] = phi[i][3      ];
			for(auto i=3; i<size1-3; i++) phi[i][size2-1] = phi[i][size2-6];
			for(auto i=3; i<size1-3; i++) phi[i][size2-2] = phi[i][size2-5];
			for(auto i=3; i<size1-3; i++) phi[i][size2-3] = phi[i][size2-4];
		}

		is_done = !(trigger_reinitialize(dx,dy,iter,phi));
		if(is_done) break;
	}
   if(v_dot_grad != nullptr){
      for(auto i=0; i<3; i++){
         for(auto j=0; j<size1-6; j++){
            delete [] v_dot_grad[i][j];
         }
         delete [] v_dot_grad[i];
      }
      delete [] v_dot_grad;
   }
	std::cout<<"End reinitialization"<<std::endl;
}

bool reinitialize::trigger_reinitialize(const double dx, const double dy, \
						 const int itr, vector<vector<double>>& phi){
	double my_min_grad, my_max_grad, min_grad, max_grad, grad_phi;
	double rdx = 1.0/dx;
	double rdy = 1.0/dy;
	auto size1 = phi.size();
	auto size2 = phi[0].size();
	if(itr>=0){
		my_min_grad =  DBL_MAX;
		my_max_grad = -DBL_MAX;
		for(auto i=3; i<size1-3; i++){
			for(auto j=3; j<size2-3; j++){
				grad_phi = \
				pow((phi[i+1][j]-phi[i-1][j])*rdx,2)+ \
				pow((phi[i][j+1]-phi[i][j-1])*rdy,2);
				my_min_grad = std::min(my_min_grad, 0.25*grad_phi);
				my_max_grad = std::max(my_max_grad, 0.25*grad_phi);

			}
		}
		min_grad = my_min_grad;
		max_grad = my_max_grad;
		if(max_grad < -1.0) return false;
		if(sqrt(min_grad)<trigger_min || sqrt(max_grad)>trigger_max) trigger=true;
		else trigger = false;
		std::cout<<"|grad phi|, min, max = "<<sqrt(min_grad)<<'\t';
		std::cout<<sqrt(max_grad)<<std::endl;
	}
	else trigger = true;
	return trigger;
}
