#include"matrix.h"
#include<iostream>
#include"advection.h"
#include"reinitialization.h"
#include"init.h"
#include<fstream>
#include<memory>
#include"utility.h"
#include"flow.h"
#include"time_step.h"
#include<omp.h>
#include<algorithm>
int main(){
	std::unique_ptr<advection> padv(new advection("WENO-5"));
	std::unique_ptr<reinitialize> prein(new reinitialize);
	std::unique_ptr<utility> putil(new utility);
	std::unique_ptr<flow> pflow(new flow);
	std::unique_ptr<time_step> ptime(new time_step);
	double rho_l= 995.65;
	double mu_l = 0.0007977;
	double rho_g= 1.161;
	double mu_g = 0.0000186;
	double H_value;
	double rho_min = 0.5*std::min(rho_l, rho_g);
	std::vector<double> x(3);
	x[0] = 0.5; x[1] = 0.5; x[2] = 0.25;
	init pp1("circle", x);
//	std::vector<double> x(1);
//	x[0] = 0.5;
//	init pp1("plane", x);
	double CFL = 0.25;
	double dt;
	int M = 80;
	int N = 120;
	double dx = 1.0/M;
	double dy = 1.5/N; // check this always
	std::vector<double> y(2);
	std::vector<std::vector<double>> phi, vel_x, vel_y, u_flow, v_flow;
	std::vector<std::vector<double>> vvx, vvy;
	vvx.resize(M+7); vvy.resize(M+2);
	for(auto i=0; i<M+7; i++) vvx[i].resize(N+2);
	for(auto i=0; i<M+2; i++) vvy[i].resize(N+7);
	std::vector<std::vector<double>> rho, mu, kappa, rho_n;
	std::vector<std::vector<double>> po, poo, pn;
	po.resize(M+2); poo.resize(M+2); pn.resize(M+2);
	for(auto i=0; i<M+2; i++){
		poo[i].resize(N+2);
		po[i].resize(N+2);
		pn[i].resize(N+2);
	}
	phi.resize(M+6);
	for(auto i=0; i<M+6; i++) phi[i].resize(N+6);
//*****************************
	vel_x.resize(M);
	vel_y.resize(M);
	kappa.resize(M+2);
	rho.resize(M+6);
	rho_n.resize(M+6);
	mu.resize(M+6);
	u_flow.resize(M+1);
	v_flow.resize(M+2);
	for(auto i=0; i<M; i++){
		vel_x[i].resize(N);
		vel_y[i].resize(N);
	}
	for(auto i=0; i<M+2; i++) kappa[i].resize(N+2);
	for(auto i=0; i<M+6; i++){
		mu[i].resize(N+6);
		rho[i].resize(N+6);
		rho_n[i].resize(N+6);
	}
	for(auto i=0; i<M+1; i++){
		u_flow[i].resize(N+2);
	}
	for(auto i=0; i<M+2; i++){
		v_flow[i].resize(N+1);
	}
	for(auto i=0; i<M+1; i++){
		for(auto j=0; j<N+2; j++){
			u_flow[i][j] = 0.0;
		}
	}
	for(auto i=0; i<M+2; i++){
		for(auto j=0; j<N+1; j++){
			v_flow[i][j] = 0.0;
		}
	}
	for(auto i=0; i<M; i++){
		for(auto j=0; j<N; j++){
			vel_x[i][j] = 0.0;
			vel_y[i][j] = 0.0;
		}
	}
//*****************************
	for(auto i=3; i<M+3; i++){
		for(auto j=3; j<N+3; j++){
			y[0] = (i-3.0+0.5)*dx;
			y[1] = (j-3.0+0.5)*dy;
			phi[i][j] = pp1.phi_init_value("circle",y);
		}
	}
	for(auto i=0; i<M+6; i++){
		for(auto j=0; j<N+6; j++){
			H_value = putil->heaviside(phi[i][j]);
//			rho[i][j] = rho_g * H_value + (1.0 - H_value) * rho_l;
//			rho[i][j] = rho_l +(rho_g-rho_l) * H_value;
//			mu[i][j] = mu_l +(mu_g-mu_l) * H_value;
			rho[i][j] = 0.5*(rho_g-rho_l)*H_value + 0.5*(rho_g+rho_l);
//			mu[i][j]  = mu_g * H_value + (1.0 - H_value) * mu_l;
			mu[i][j] =  0.5*(mu_g-mu_l) * H_value + 0.5*(mu_g + mu_l);
		}
	}
	std::ofstream myfile;
	myfile.open("out0.dat");
	for(auto i=3; i<M+3; i++){
		for(auto j=3; j<N+3; j++){
			myfile<<phi[i][j]<<'\t';
		}
		myfile<<std::endl;
	}
	myfile.close();
	dt = 0.01;
	for(auto i=0; i<M+1; i++) u_flow[i][N+1] = -u_flow[i][N];
	double x_p, y_p;
	std::vector<std::vector<double>> uuu,vvv;
	uuu.resize(M+1); vvv.resize(M+2);
	for(auto i=0; i<M+1; i++) uuu[i].resize(N+2);
	for(auto i=0; i<M+2; i++) vvv[i].resize(N+1);
// make them zero
	for(auto i=0; i<M+1; i++){
		for(auto j=0; j<N+2; j++) uuu[i][j] = 0.0;
	}

	for(auto i=0; i<M+2; i++){
		for(auto j=0; j<N+1; j++) vvv[i][j] = 0.0;
	}

	double pi = 4.0 * atan(1.0);
	for(auto idt=0; idt<300; idt++){
		for(auto i=0; i<M+1; i++){
			for(auto j=1; j<N+1; j++){
				x_p = static_cast<double>(i)*dx;
				y_p = static_cast<double>(j)*dy;	
				uuu[i][j]  = pi*sin(2.0*pi*y_p)*pow(sin(pi*x_p),2)*sin(idt*dt);
			}
		}
		for(auto i=1; i<M+1; i++){
			for(auto j=0; j<N+1; j++){
				x_p = static_cast<double>(i)*dx;
				y_p = static_cast<double>(j)*dy;	
				vvv[i][j]  =-pi*sin(2.0*pi*x_p)*pow(sin(pi*y_p),2)*sin(idt*dt);
			}
		}
		std::cout<<"time step  "<<idt+1<<"  of  "<<"..."<<std::endl;
		putil->curvature(dx, dy, kappa, phi);
		for(auto i=0; i<M+6; i++){
			for(auto j=0; j<N+6; j++){
				H_value = putil->heaviside(phi[i][j]);
//				rho_n[i][j] = rho_g * H_value + (1.0 - H_value) * rho_l;
//				rho_n[i][j] = rho_l +(rho_g-rho_l) * H_value;
				rho_n[i][j] = 0.5*(rho_g-rho_l)*H_value + 0.5*(rho_g+rho_l);
			}
		}
//		pflow->flow_solver(M, N, dt,rho_min, rho_n,rho, mu, uuu, vvv,\
								 poo, po, pn,kappa, phi);
//		for(auto i=0; i<M+1; i++){
//			for(auto j=0; j<N+2; j++) u_flow[i][j] = 1.0;
//		}

//#######################

//		for(auto i=0; i<M; i++){
//			for(auto j=1; j<N+1; j++){
//				vel_x[i][j-1] = 0.5 * (u_flow[i][j] + u_flow[i+1][j]);
//			}
//		}
//		for(auto i=1; i<M+1; i++){
//			for(auto j=0; j<N; j++){
//				vel_y[i-1][j] = 0.5 * (v_flow[i][j] + v_flow[i][j+1]);
//			}
//		}
		ptime->simulate(dt, dx, dy, CFL, uuu, vvv, phi);
		for(auto i=3; i<M+4; i++){
			for(auto j=0; j<N+2; j++) vvx[i][j] = u_flow[i-3][j];
		}
		for(auto j=0; j<N+2; j++) vvx[0][j] = -vvx[6][j];
		for(auto j=0; j<N+2; j++) vvx[1][j] = -vvx[5][j];
		for(auto j=0; j<N+2; j++) vvx[2][j] = -vvx[4][j];

		for(auto i=0; i<M+2; i++){
			for(auto j=3; j<N+4; j++) vvy[i][j] = v_flow[i][j-3];
		}
		for(auto i=0; i<M+2; i++) vvy[i][0] = -vvy[i][6];
		for(auto i=0; i<M+2; i++) vvy[i][1] = -vvy[i][5];
		for(auto i=0; i<M+2; i++) vvy[i][2] = -vvy[i][4];
		for(auto i=0; i<M+6; i++){
			for(auto j=0; j<N+6; j++){
				H_value = putil->heaviside(phi[i][j]);
//				rho[i][j] = rho_g * H_value + (1.0 - H_value) * rho_l;
//				rho[i][j] = rho_l +(rho_g-rho_l) * H_value;
//				mu[i][j] = mu_l +(mu_g-mu_l) * H_value;
				rho[i][j] = 0.5*(rho_g-rho_l)*H_value + 0.5*(rho_g+rho_l);
//				mu[i][j]  = mu_g * H_value + (1.0 - H_value) * mu_l;
				mu[i][j] =  0.5*(mu_g-mu_l) * H_value + 0.5*(mu_g + mu_l);
			}
		}
//		ptime->simulate(dt, dx, dy, CFL, u_flow, v_flow, phi);

		for(auto i=3; i<M+4; i++){
			for(auto j=0; j<N+2; j++) u_flow[i-3][j] = vvx[i][j];
		}

		for(auto i=0; i<M+2; i++){
			for(auto j=3; j<N+4; j++) v_flow[i][j-3] =vvy[i][j];
		}

		if(idt ==27){
			myfile.open("out1.dat");
			for(auto i=3; i<M+3; i++){
				for(auto j=3; j<N+3; j++){
					myfile<<phi[i][j]<<'\t';
				}
				myfile<<std::endl;
			}
			myfile.close();
		}
		if(idt == 57){
			myfile.open("out2.dat");
			for(auto i=3; i<M+3; i++){
				for(auto j=3; j<N+3; j++){
					myfile<<phi[i][j]<<'\t';
				}
				myfile<<std::endl;
			}
			myfile.close();
		}
		if(idt == 87){
			myfile.open("out3.dat");
			for(auto i=3; i<M+3; i++){
				for(auto j=3; j<N+3; j++){
					myfile<<phi[i][j]<<'\t';
				}
				myfile<<std::endl;
			}
			myfile.close();
		}
		if(idt == 117){
			myfile.open("out4.dat");
			for(auto i=3; i<M+3; i++){
				for(auto j=3; j<N+3; j++){
					myfile<<phi[i][j]<<'\t';
				}
				myfile<<std::endl;
			}
			myfile.close();
		}
		if(idt == 49){
			myfile.open("out5.dat");
			for(auto i=3; i<M+3; i++){
				for(auto j=3; j<N+3; j++){
					myfile<<phi[i][j]<<'\t';
				}
				myfile<<std::endl;
			}
			myfile.close();
		}
		if(idt == 150){
			myfile.open("out6.dat");
			for(auto i=3; i<M+3; i++){
				for(auto j=3; j<N+3; j++){
					myfile<<phi[i][j]<<'\t';
				}
				myfile<<std::endl;
			}
			myfile.close();
		}
		if(idt == 250){
			myfile.open("out7.dat");
			for(auto i=3; i<M+3; i++){
				for(auto j=3; j<N+3; j++){
					myfile<<phi[i][j]<<'\t';
				}
				myfile<<std::endl;
			}
			myfile.close();
		}
		if(idt == 300){
			myfile.open("out8.dat");
			for(auto i=3; i<M+3; i++){
				for(auto j=3; j<N+3; j++){
					myfile<<phi[i][j]<<'\t';
				}
				myfile<<std::endl;
			}
			myfile.close();
		}
		if(idt == 89){
			myfile.open("out9.dat");
			for(auto i=3; i<M+3; i++){
				for(auto j=3; j<N+3; j++){
					myfile<<phi[i][j]<<'\t';
				}
				myfile<<std::endl;
			}
			myfile.close();
		}
		if(idt == 99){
			myfile.open("out10.dat");
			for(auto i=3; i<M+3; i++){
				for(auto j=3; j<N+3; j++){
					myfile<<phi[i][j]<<'\t';
				}
				myfile<<std::endl;
			}
			myfile.close();
		}
	}
	myfile.open("velocity.dat");
	for(auto i=0; i<M+1; i++){
		for(auto j=0; j<N+2; j++){
			myfile<<u_flow[i][j]<<'\t';
		}
		myfile<<std::endl;
	}
	myfile.close();
//	myfile.open("out.dat");
//	for(auto i=3; i<M+3; i++){
//		for(auto j=3; j<N+3; j++){
//			myfile<<phi[i][j]<<'\t';
//		}
//		myfile<<std::endl;
//	}
//	myfile.close();
	return 0;
}
