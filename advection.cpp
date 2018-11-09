#include"advection.h"
#include"weno.h"
#include<iostream>
#include<memory>
#include<omp.h>
advection::advection(std::string s){
	advection_method = s;
}

void advection::advect_sol(const double dt, const double dx, const double dy, \
									const vector<vector<double>>& vel_x, \
									const vector<vector<double>>& vel_y, \
									vector<vector<double>>& phi){
	std::unique_ptr<weno> pweno(new weno);
	std::cout<<"Starting Advection"<<std::endl;
	const int RK_advection   = 3;
	auto size1 = phi.size();
	auto size2 = phi[0].size();
	double rdx = 1.0/dx;
	double rdy = 1.5/dy; // check this always
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
	const double r23  = 2.0/3.0;
	const double r112 = 1.0/12.0;
	std::vector<std::vector<double>> alphaRK;
	alphaRK.resize(3);
	for(auto i=0; i<3; i++) alphaRK[i].resize(3);
	// set the coefficients of the TVD-RK3 

   alphaRK[0][0]=dt;        alphaRK[1][0]=0.0;       alphaRK[2][0]=0.0;
   alphaRK[0][1]=-0.75*dt;  alphaRK[1][1]=0.25*dt;   alphaRK[2][1]=0.0;
   alphaRK[0][2]=-r112*dt;  alphaRK[1][2]=-r112*dt;  alphaRK[2][2]=r23*dt;
	
	//set ghost nodes
	for(auto j=3; j<size2-3; j++) phi[0][j] = phi[5][j];
	for(auto j=3; j<size2-3; j++) phi[1][j] = phi[4][j];
	for(auto j=3; j<size2-3; j++) phi[2][j] = phi[3][j];
	for(auto j=3; j<size2-3; j++) phi[size1-1][j] = phi[size1-6][j];
	for(auto j=3; j<size2-3; j++) phi[size1-2][j] = phi[size1-5][j];
	for(auto j=3; j<size2-3; j++)	phi[size1-3][j] = phi[size1-4][j];

	for(auto i=3; i<size1-3; i++) phi[i][0] = phi[i][5];
	for(auto i=3; i<size1-3; i++) phi[i][1] = phi[i][4];
	for(auto i=3; i<size1-3; i++) phi[i][2] = phi[i][3];
	for(auto i=3; i<size1-3; i++) phi[i][size2-1] = phi[i][size2-6];
	for(auto i=3; i<size1-3; i++) phi[i][size2-2] = phi[i][size2-5];
	for(auto i=3; i<size1-3; i++) phi[i][size2-3] = phi[i][size2-4];
	vector<double> grad_phi(2,0.0);
	vector<double> d_phix (6, 0.0);
	vector<double> d_phiy (6, 0.0);
	double d2_phix;
	double d2_phiy;
	const int WENO_advection = 5;
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

				if(vel_x[i-3][j-3] <= 0.0){
					grad_phi[0] = \
					r112*(-d_phix[1]+7.0*d_phix[2]+7.0*d_phix[3]-d_phix[4]) + \
					pweno->phi_weno_5th(d_phix[5]-d_phix[4],d_phix[4]-d_phix[3], \
											  d_phix[3]-d_phix[2],d_phix[2]-d_phix[1]);
				}
				else{
					grad_phi[0] = \
					r112*(-d_phix[1]+7.0*d_phix[2]+7.0*d_phix[3]-d_phix[4]) - \
					pweno->phi_weno_5th(d_phix[1]-d_phix[0],d_phix[2]-d_phix[1], \
											  d_phix[3]-d_phix[2],d_phix[4]-d_phix[3]);
				}
				if(vel_y[i-3][j-3] <= 0.0){
					grad_phi[1] = \
					r112*(-d_phiy[1]+7.0*d_phiy[2]+7.0*d_phiy[3]-d_phiy[4]) + \
					pweno->phi_weno_5th(d_phiy[5]-d_phiy[4],d_phiy[4]-d_phiy[3], \
											  d_phiy[3]-d_phiy[2],d_phiy[2]-d_phiy[1]);
				}
				else{
					grad_phi[1] = \
					r112*(-d_phiy[1]+7.0*d_phiy[2]+7.0*d_phiy[3]-d_phiy[4]) - \
					pweno->phi_weno_5th(d_phiy[1]-d_phiy[0],d_phiy[2]-d_phiy[1], \
											  d_phiy[3]-d_phiy[2],d_phiy[4]-d_phiy[3]);
				}

				v_dot_grad[irk][i-3][j-3] = -(grad_phi[0] * vel_x[i-3][j-3] + \
														grad_phi[1] * vel_y[i-3][j-3]);
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
	if(v_dot_grad != nullptr){
		for(auto i=0; i<3; i++){
			for(auto j=0; j<size1-6; j++){
				delete [] v_dot_grad[i][j];
			}
			delete [] v_dot_grad[i];
		}
		delete [] v_dot_grad;
	}
	std::cout<<"End Advection"<<std::endl;
}

void advection::advect_sol(const double dt, const double dx, const double dy, \
									const vector<vector<double>>& vel_x, \
									const vector<vector<double>>& vel_y, \
									vector<vector<double>>& phi){
	std::unique_ptr<weno> pweno(new weno);
	std::cout<<"Starting Advection"<<std::endl;
	const int RK_advection   = 3;
	auto size1 = phi.size();
	auto size2 = phi[0].size();
	double rdx = 1.0/dx;
	double rdy = 1.0/dy;
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
	const double r23  = 2.0/3.0;
	const double r112 = 1.0/12.0;
	std::vector<std::vector<double>> alphaRK;
	alphaRK.resize(3);
	for(auto i=0; i<3; i++) alphaRK[i].resize(3);
	// set the coefficients of the TVD-RK3 

   alphaRK[0][0]=dt;        alphaRK[1][0]=0.0;       alphaRK[2][0]=0.0;
   alphaRK[0][1]=-0.75*dt;  alphaRK[1][1]=0.25*dt;   alphaRK[2][1]=0.0;
   alphaRK[0][2]=-r112*dt;  alphaRK[1][2]=-r112*dt;  alphaRK[2][2]=r23*dt;
	
	//set ghost nodes
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
	vector<double> grad_phi(2,0.0);
	vector<double> d_phix (6, 0.0);
	vector<double> d_phiy (6, 0.0);
	double d2_phix;
	double d2_phiy;
	const int WENO_advection = 5;
	double mean1, mean2;
	for(auto irk=0; irk<3; irk++){
#pragma opm parallel for
		for(auto i=3; i<size1-3; i++){
			for(auto j=3; j<size2-3; j++){
				mean1 = 0.5 * (vel_x[i-1][j-2] + vel_x[i-2][j-2]);
				mean2 = 0.5 * (vel_x[i-2][j-2] + vel_x[i-3][j-2]);
				d_phix[0] = (mean1*phi[i-2][j] - mean2*phi[i-3][j]) * rdx;

				mean1 = 0.5 * (vel_x[i  ][j-2] + vel_x[i-1][j-2]);
				mean2 = 0.5 * (vel_x[i-1][j-2] + vel_x[i-2][j-2]);
				d_phix[1] = (mean1*phi[i-1][j] - mean2*phi[i-2][j]) * rdx;

				mean1 = 0.5 * (vel_x[i+1][j-2] + vel_x[i  ][j-2]);
				mean2 = 0.5 * (vel_x[i  ][j-2] + vel_x[i-1][j-2]);
				d_phix[2] = (mean1*phi[i  ][j] - mean2*phi[i-1][j]) * rdx;

				mean1 = 0.5 * (vel_x[i+2][j-2] + vel_x[i+1][j-2]);
				mean2 = 0.5 * (vel_x[i+1][j-2] + vel_x[i][j-2]);
				d_phix[3] = (mean1*phi[i+1][j] - mean2*phi[i  ][j]) * rdx;

				mean1 = 0.5 * (vel_x[i+3][j-2] + vel_x[i+2][j-2]);
				mean2 = 0.5 * (vel_x[i+2][j-2] + vel_x[i+1][j-2]);
				d_phix[4] = (mean1*phi[i+2][j] - mean2*phi[i+1][j]) * rdx;

				mean1 = 0.5 * (vel_x[i+4][j-2] + vel_x[i+3][j-2]);
				mean2 = 0.5 * (vel_x[i+3][j-2] + vel_x[i+2][j-2]);
				d_phix[5] = (mean1*phi[i+3][j] - mean2*phi[i+2][j]) * rdx;

				mean1 = 0.5 * (vel_y[i-2][j-1] + vel_y[i-2][j-2]);
				mean2 = 0.5 * (vel_y[i-2][j-2] + vel_y[i-2][j-3]);
				d_phiy[0] = (mean1*phi[i][j-2] - mean2*phi[i][j-3]) * rdy;

				mean1 = 0.5 * (vel_y[i-2][j  ] + vel_y[i-2][j-1]);
				mean2 = 0.5 * (vel_y[i-2][j-1] + vel_y[i-2][j-2]);
				d_phiy[1] = (mean1*phi[i][j-1] - mean2*phi[i][j-2]) * rdy;

				mean1 = 0.5 * (vel_y[i-2][j+1] + vel_y[i-2][j  ]);
				mean2 = 0.5 * (vel_y[i-2][j  ] + vel_y[i-2][j-1]);
				d_phiy[2] = (mean1*phi[i][j  ] - mean2*phi[i][j-1]) * rdy;

				mean1 = 0.5 * (vel_y[i-2][j+2] + vel_y[i-2][j+1]);
				mean2 = 0.5 * (vel_y[i-2][j+1] + vel_y[i-2][j]);
				d_phiy[3] = (mean1*phi[i][j+1] - mean2*phi[i][j  ]) * rdy;

				mean1 = 0.5 * (vel_y[i-2][j+3] + vel_y[i-2][j+2]);
				mean2 = 0.5 * (vel_y[i-2][j+2] + vel_y[i-2][j+1]);
				d_phiy[4] = (mean1*phi[i][j+2] - mean2*phi[i][j+1]) * rdy;

				mean1 = 0.5 * (vel_y[i-2][j+4] + vel_y[i-2][j+3]);
				mean2 = 0.5 * (vel_y[i-2][j+3] + vel_y[i-2][j+2]);
				d_phiy[5] = (mean1*phi[i][j+3] - mean2*phi[i][j+2]) * rdy;

				if(vel_x[i-3][j-3] <= 0.0){
					grad_phi[0] = \
					r112*(-d_phix[1]+7.0*d_phix[2]+7.0*d_phix[3]-d_phix[4]) + \
					pweno->phi_weno_5th(d_phix[5]-d_phix[4],d_phix[4]-d_phix[3], \
											  d_phix[3]-d_phix[2],d_phix[2]-d_phix[1]);
				}
				else{
					grad_phi[0] = \
					r112*(-d_phix[1]+7.0*d_phix[2]+7.0*d_phix[3]-d_phix[4]) - \
					pweno->phi_weno_5th(d_phix[1]-d_phix[0],d_phix[2]-d_phix[1], \
											  d_phix[3]-d_phix[2],d_phix[4]-d_phix[3]);
				}
				if(vel_y[i-3][j-3] <= 0.0){
					grad_phi[1] = \
					r112*(-d_phiy[1]+7.0*d_phiy[2]+7.0*d_phiy[3]-d_phiy[4]) + \
					pweno->phi_weno_5th(d_phiy[5]-d_phiy[4],d_phiy[4]-d_phiy[3], \
											  d_phiy[3]-d_phiy[2],d_phiy[2]-d_phiy[1]);
				}
				else{
					grad_phi[1] = \
					r112*(-d_phiy[1]+7.0*d_phiy[2]+7.0*d_phiy[3]-d_phiy[4]) - \
					pweno->phi_weno_5th(d_phiy[1]-d_phiy[0],d_phiy[2]-d_phiy[1], \
											  d_phiy[3]-d_phiy[2],d_phiy[4]-d_phiy[3]);
				}

				v_dot_grad[irk][i-3][j-3] = -(grad_phi[0] + grad_phi[1]);
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
	if(v_dot_grad != nullptr){
		for(auto i=0; i<3; i++){
			for(auto j=0; j<size1-6; j++){
				delete [] v_dot_grad[i][j];
			}
			delete [] v_dot_grad[i];
		}
		delete [] v_dot_grad;
	}
	std::cout<<"End Advection"<<std::endl;
}
void advection::advect_sol(const double dt, const double dx, const double dy, \
									const vector<vector<double>>& vel_x, \
									const vector<vector<double>>& vel_y, \
									vector<vector<double>>& phi){
	std::unique_ptr<weno> pweno(new weno);
	std::cout<<"Starting Advection"<<std::endl;
	const int RK_advection   = 3;
	auto size1 = phi.size();
	auto size2 = phi[0].size();
	double rdx = 1.0/dx;
	double rdy = 1.0/dy;
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
	const double r23  = 2.0/3.0;
	const double r112 = 1.0/12.0;
	std::vector<std::vector<double>> alphaRK;
	alphaRK.resize(3);
	for(auto i=0; i<3; i++) alphaRK[i].resize(3);
	// set the coefficients of the TVD-RK3 

   alphaRK[0][0]=dt;        alphaRK[1][0]=0.0;       alphaRK[2][0]=0.0;
   alphaRK[0][1]=-0.75*dt;  alphaRK[1][1]=0.25*dt;   alphaRK[2][1]=0.0;
   alphaRK[0][2]=-r112*dt;  alphaRK[1][2]=-r112*dt;  alphaRK[2][2]=r23*dt;
	
	//set ghost nodes
	for(auto j=3; j<size2-3; j++) phi[0][j] = phi[5][j];
	for(auto j=3; j<size2-3; j++) phi[1][j] = phi[4][j];
	for(auto j=3; j<size2-3; j++) phi[2][j] = phi[3][j];
	for(auto j=3; j<size2-3; j++) phi[size1-1][j] = phi[size1-6][j];
	for(auto j=3; j<size2-3; j++) phi[size1-2][j] = phi[size1-5][j];
	for(auto j=3; j<size2-3; j++)	phi[size1-3][j] = phi[size1-4][j];

	for(auto i=3; i<size1-3; i++) phi[i][0] = phi[i][5];
	for(auto i=3; i<size1-3; i++) phi[i][1] = phi[i][4];
	for(auto i=3; i<size1-3; i++) phi[i][2] = phi[i][3];
	for(auto i=3; i<size1-3; i++) phi[i][size2-1] = phi[i][size2-6];
	for(auto i=3; i<size1-3; i++) phi[i][size2-2] = phi[i][size2-5];
	for(auto i=3; i<size1-3; i++) phi[i][size2-3] = phi[i][size2-4];
	vector<double> grad_phi(2,0.0);
	vector<double> d_phix (6, 0.0);
	vector<double> d_phiy (6, 0.0);
	double d2_phix;
	double d2_phiy;
	const int WENO_advection = 5;
	for(auto irk=0; irk<3; irk++){
		for(auto i=3; i<size1-3; i++){
			for(auto j=3; j<size2-3; j++){
				d_phix[0] = (phi[i-1][j] - phi[i-2][j]) * rdx;
				d_phix[1] = (phi[i][j] - phi[i-1][j]) * rdx;
				d_phix[2] = (phi[i+1][j] - phi[i][j]) * rdx;
				d_phix[3] = (phi[i+2][j] - phi[i+1  ][j]) * rdx;
				d_phix[4] = (phi[i+3][j] - phi[i+2][j]) * rdx;
//				d_phix[5] = (phi[i+4][j] - phi[i+3][j]) * rdx;

				d_phiy[0] = (phi[i][j-1] - phi[i][j-2]) * rdy;
				d_phiy[1] = (phi[i][j] - phi[i][j-1]) * rdy;
				d_phiy[2] = (phi[i][j+1] - phi[i][j]) * rdy;
				d_phiy[3] = (phi[i][j+2] - phi[i][j+1]) * rdy;
				d_phiy[4] = (phi[i][j+3] - phi[i][j+2]) * rdy;
//				d_phiy[5] = (phi[i][j+4] - phi[i][j+3]) * rdy;

				if(vel_x[i-3][j-3] <= 0.0){
					grad_phi[0] = \
					r112*(-d_phix[1]+7.0*d_phix[2]+7.0*d_phix[3]-d_phix[4]) + \
					pweno->phi_weno_5th(d_phix[4]-d_phix[3],d_phix[3]-d_phix[2], \
											  d_phix[2]-d_phix[1],d_phix[1]-d_phix[0]);
				}
				else{
					grad_phi[0] = \
					r112*(-d_phix[1]+7.0*d_phix[2]+7.0*d_phix[3]-d_phix[4]) - \
					pweno->phi_weno_5th(d_phix[1]-d_phix[0],d_phix[2]-d_phix[1], \
											  d_phix[3]-d_phix[2],d_phix[4]-d_phix[3]);
				}
				if(vel_y[i-3][j-3] <= 0.0){
					grad_phi[1] = \
					r112*(-d_phiy[1]+7.0*d_phiy[2]+7.0*d_phiy[3]-d_phiy[4]) + \
					pweno->phi_weno_5th(d_phiy[4]-d_phiy[3],d_phiy[3]-d_phiy[2], \
											  d_phiy[2]-d_phiy[1],d_phiy[1]-d_phiy[0]);
				}
				else{
					grad_phi[1] = \
					r112*(-d_phiy[1]+7.0*d_phiy[2]+7.0*d_phiy[3]-d_phiy[4]) - \
					pweno->phi_weno_5th(d_phiy[1]-d_phiy[0],d_phiy[2]-d_phiy[1], \
											  d_phiy[3]-d_phiy[2],d_phiy[4]-d_phiy[3]);
				}

				v_dot_grad[irk][i-3][j-3] = -(grad_phi[0] * vel_x[i-3][j-3] + \
														grad_phi[1] * vel_y[i-3][j-3]);
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
	if(v_dot_grad != nullptr){
		for(auto i=0; i<3; i++){
			for(auto j=0; j<size1-6; j++){
				delete [] v_dot_grad[i][j];
			}
			delete [] v_dot_grad[i];
		}
		delete [] v_dot_grad;
	}
	std::cout<<"End Advection"<<std::endl;
}

