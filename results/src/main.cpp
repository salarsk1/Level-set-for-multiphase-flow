//Written by Salar Safarkhani

#include "datag.h"
#include "parallel.h"
#include "lit_coupler.h"
#include "io.h"
#include "reinit.h"
#include "bound.h"
#include <memory>
#include <fstream>
#include <iostream>
//#include "/home/salar/developments/cfd_561/final/flow.h"
void setVelocity(const int);
void calcShapeError();
int main(int argc, char *argv[]){
   int is;
   double delta_t, fs_dx_max, end_time, dt_dump;
   std::vector<double> r1_fs(2);
   std::string solvername;
   bool triggerSol;
   std::unique_ptr<global_variable> pgv(new global_variable);
   std::unique_ptr<parallel> pp(new parallel);
   std::unique_ptr<lit_coupler> plco(new lit_coupler);
   std::unique_ptr<io> pio(new io);
   std::unique_ptr<reinit> prein(new reinit);
   std::unique_ptr<toolbox> ptool(new toolbox);
   for(auto i=0; i<3; i++) pgv->fs_dxyzmax[i] = 1.0/128.0;
   pgv->time = 0.0;
   pgv->fs_max_drop_vol = 0.0;
   plco->lit_initialize(argc, argv);
		
//   pgv->init_shape = "notched_circle";
//   pgv->xyzs_sg[0] = 0.0; pgv->xyzs_sg[1] = 0.0; pgv->xyzs_sg[2] = 0.0;
//   pgv->xyze_sg[0] = 1.0; pgv->xyze_sg[1] = 1.0; pgv->xyze_sg[2] = 0.01;
//   pgv->ijkm_gl[0] = 128; pgv->ijkm_gl[1] = 128; pgv->ijkm_gl[2] = 1;
//   pgv->ijkm_bl[0] = 128; pgv->ijkm_bl[1] = 128;  pgv->ijkm_bl[2] = 1;
//   pgv->init_center[0] = 0.5; pgv->init_center[1] = 0.75; pgv->init_center[2] = 0.0;
//   pgv->init_radius = 0.15; pgv->init_width = 0.05; pgv->init_height = 0.25;

//   pgv->init_shape = "circle";
//   pgv->xyzs_sg[0] = 0.0; pgv->xyzs_sg[1] = 0.0; pgv->xyzs_sg[2] = 0.0;
//   pgv->xyze_sg[0] = 1.0; pgv->xyze_sg[1] = 1.0; pgv->xyze_sg[2] = 0.0078125;
//   pgv->ijkm_gl[0] = 100; pgv->ijkm_gl[1] = 100; pgv->ijkm_gl[2] = 1;
//	pgv->ijkm_bl[0] = 100; pgv->ijkm_bl[1] = 100; pgv->ijkm_bl[2] = 1;
//   pgv->init_center[0] = 0.5; pgv->init_center[1] = 0.5; pgv->init_center[2] = 0.0;
//	pgv->init_radius = 0.2;

   pgv->init_shape = "plane";
   pgv->xyzs_sg[0] = 0.0; pgv->xyzs_sg[1] = 0.0; pgv->xyzs_sg[2] = 0.0;
   pgv->xyze_sg[0] = 1.0; pgv->xyze_sg[1] = 1.0; pgv->xyze_sg[2] = 0.0078125;
   pgv->ijkm_gl[0] = 32; pgv->ijkm_gl[1] = 32; pgv->ijkm_gl[2] = 1;
   pgv->ijkm_bl[0] = 32; pgv->ijkm_bl[1] = 32; pgv->ijkm_bl[2] = 1;
	pgv->init_center[0] = 0.5; pgv->init_center[1] = 0.5; pgv->init_center[2] = 0.0;
	pgv->init_mod = 1;
//   pgv->init_shape = "ellipse";
//   pgv->xyzs_sg[0] = 0.0; pgv->xyzs_sg[1] = 0.0; pgv->xyzs_sg[2] = 0.0;
//   pgv->xyze_sg[0] = 1.0; pgv->xyze_sg[1] = 1.0; pgv->xyze_sg[2] = 0.0078125;
//   pgv->ijkm_gl[0] = 128; pgv->ijkm_gl[1] = 128; pgv->ijkm_gl[2] = 1;
//	pgv->ijkm_bl[0] = 128; pgv->ijkm_bl[1] = 128; pgv->ijkm_bl[2] = 1;
//   pgv->init_center[0] = 0.5; pgv->init_center[1] = 0.5; pgv->init_center[2] = 0.0;
//	pgv->init_a = 0.1;
//	pgv->init_b = 0.3;
//   pgv->init_shape = "ring";
//   pgv->xyzs_sg[0] = 0.0; pgv->xyzs_sg[1] = 0.0; pgv->xyzs_sg[2] = 0.0;
//   pgv->xyze_sg[0] = 1.0; pgv->xyze_sg[1] = 1.0; pgv->xyze_sg[2] = 0.0078125;
//   pgv->ijkm_gl[0] = 128; pgv->ijkm_gl[1] = 128; pgv->ijkm_gl[2] = 1;
//	pgv->ijkm_bl[0] = 128; pgv->ijkm_bl[1] = 128; pgv->ijkm_bl[2] = 1;
//   pgv->init_center[0] = 0.5; pgv->init_center[1] = 0.5; pgv->init_center[2] = 0.0;
//	pgv->init_radius = 0.3;
//	pgv->init_width  = 0.1;

	pgv->nghost = 3;
	pgv->xyzs_init_bb[0] = pgv->xyzs_sg[0];
	pgv->xyzs_init_bb[1] = pgv->xyzs_sg[1];
	pgv->xyzs_init_bb[2] = pgv->xyzs_sg[2];

	pgv->xyze_init_bb[0] = pgv->xyze_sg[0];
	pgv->xyze_init_bb[1] = pgv->xyze_sg[1];
	pgv->xyze_init_bb[2] = pgv->xyze_sg[2];
	pgv->sg_periodic[0] = false;
	pgv->sg_periodic[1] = false;
	pgv->sg_periodic[2] = false;
   prein->triggerGradMin = 1.0e-4;
   prein->triggerGradMax = 2.0;
   prein->min_iter_reinit = 1;
   pgv->nVelFilter = 0;

	pio->varDumpName[0] = "WRITE_LIT_STEP";
	pio->varDumpName[1] = "ENSIGHT";
	pio->varDumpName[2] = "dc";
	pio->varDumpName[3] = "512";
	pio->varDumpName[4] = "T-BAND";
	pio->varDumpName[5] = "G";
//	pio->varDumpName[6] = "V";
	pgv->dump_band = 'T';
	pgv->schemeAdvect = "WENO-5";
	pgv->schemeReinit = "WENO-5";
	pgv->reinit_solver = "PDE";
	pgv->case_name = "lit";
	pgv->dump_dtime_fs = 0.5;
   plco->lit_initialize2(argc, argv, "input");

	if(pp->myrank == 0){
		auto dd1 = pgv->b->NinA;
		std::ofstream myfile;
	   std::cout<<"mr="<<pp->myrank<<std::endl;
	   setVelocity(0);
	   delta_t = 0.001;//2.0*pgv->pi/628.0;
		myfile.open("A1.dat");
		for(auto i=0; i<dd1; i++){
			myfile<<pgv->xc[pgv->Gn[i].ijk[0]-(1-pgv->nghost)]<<'\t';
			myfile<<pgv->yc[pgv->Gn[i].ijk[1]-(1-pgv->nghost)]<<'\t';
			myfile<<pgv->Gn[i].G<<std::endl;
		}
		myfile.close();
		dd1 = pgv->b->NinT;
		myfile.open("T1.dat");
		for(auto i=0; i<dd1; i++){
			myfile<<pgv->xc[pgv->Gn[i].ijk[0]-(1-pgv->nghost)]<<'\t';
			myfile<<pgv->yc[pgv->Gn[i].ijk[1]-(1-pgv->nghost)]<<'\t';
			myfile<<pgv->Gn[i].G<<std::endl;
		}
		myfile.close();
		dd1 = pgv->b->NinX;
		myfile.open("X1.dat");
		for(auto i=0; i<dd1; i++){
			myfile<<pgv->xc[pgv->Gn[i].ijk[0]-(1-pgv->nghost)]<<'\t';
			myfile<<pgv->yc[pgv->Gn[i].ijk[1]-(1-pgv->nghost)]<<'\t';
			myfile<<pgv->Gn[i].G<<std::endl;
		}
		myfile.close();
		myfile.close();
	}
	flow myflow;
	std::vector<std::vector<double>> uin1, uin2, vin1, vin2;
	uin1.resize(33);
	uin2.resize(33);
	for(auto i=0 ;i<33; i++){
		uin1[i].resize(34);
		uin2[i].resize(34);
	}
	vin1.resize(34);
	vin2.resize(34);
	for(auto i=0; i<34; i++){
		vin1[i].resize(33);
		vin2[i].resize(33);
	}
	double xpi, ypi;
   for(auto is=0; is<1000; is++){
		pgv->time = static_cast<double>(is)*delta_t;
      if(pp->myrank == 0) std::cout<<"Step "<<is<<" of 600, t ="<<pgv->time<<std::endl;
//		myflow.flow_solver(32,32,delta_t,100,uin1, vin1);
//		myflow.flow_solver(32,32,delta_t,100,uin2, vin2);
		setVelocity(is);
	   for(auto ibl=0; ibl<pgv->nbl; ibl++){
   	   pgv->setBlockPointers(ibl);
	      for(auto ic=0; ic<pgv->b->NinT; ic++){
				xpi = pgv->Gn[ic].ijk[0]-(1-pgv->nghost)-3;
				ypi = pgv->Gn[ic].ijk[1]-(1-pgv->nghost)-3;
				if(pgv->Gn[ic].G >= 0.0){
					pgv->Gn[ic].V[0] = 1.0;//uin1[xpi][ypi];
					pgv->Gn[ic].V[1] = 0.0;//vin1[xpi][ypi];
				}
				else{
					pgv->Gn[ic].V[0] = 1.0;//uin2[xpi][ypi];
					pgv->Gn[ic].V[1] = 0.0;//vin2[xpi][ypi];
				}
				pgv->Gn[ic].V[2] = 0.0;
			}
		}
      if((is+1)%256 == 0) triggerSol = true;
      else triggerSol = false;
      plco->litRunIteration(delta_t, triggerSol, 2);
		int dd1;
		int dd2 = pgv->b->NinT;
		double *kappa;
		kappa = new double[dd2];
//		ptool->curvature(dd2, kappa);
//		for(auto nn=0; nn<dd2; nn++) std::cout<<kappa[nn]<<std::endl;
		std::ofstream myfile;
		if(is==49||is==99||is==149||is==199||is==249||is==299||is==349||is==399||is==449|| \
			is==499 || is==549 || is==599){
			dd1 = pgv->b->NinA;
			if		 (is==49)  myfile.open("A2.dat");
			else if(is==99)  myfile.open("A3.dat");
			else if(is==149) myfile.open("A4.dat");
			else if(is==199) myfile.open("A5.dat");
			else if(is==249) myfile.open("A6.dat");
			else if(is==299) myfile.open("A7.dat");
			else if(is==349) myfile.open("A8.dat");
			else if(is==399) myfile.open("A9.dat");
			else if(is==449) myfile.open("A10.dat");
			else if(is==499) myfile.open("A11.dat");
			else if(is==549) myfile.open("A12.dat");
			else if(is==599) myfile.open("A13.dat");
			for(auto i=0; i<dd1; i++){
				myfile<<pgv->xc[pgv->Gn[i].ijk[0]-(1-pgv->nghost)]<<'\t';
				myfile<<pgv->yc[pgv->Gn[i].ijk[1]-(1-pgv->nghost)]<<'\t';
				myfile<<pgv->Gn[i].G<<std::endl;
			}
			myfile.close();
		}
/*		if(is==156||is==313||is==470||is==627){
			dd1 = pgv->b->NinT;
			if(is==156) 		myfile.open("T2.dat");
			else if(is==313)  myfile.open("T3.dat");
			else if(is==470)  myfile.open("T4.dat");
			else if(is==627) myfile.open("T5.dat");
			for(auto i=0; i<dd1; i++){
				myfile<<pgv->xc[pgv->Gn[i].ijk[0]-(1-pgv->nghost)]<<'\t';
				myfile<<pgv->yc[pgv->Gn[i].ijk[1]-(1-pgv->nghost)]<<'\t';
				myfile<<pgv->Gn[i].G<<std::endl;
			}
			myfile.close();
		}
		if(is==156||is==313||is==470||is==627){
			dd1 = pgv->b->NinX;
			if(is==156) 		myfile.open("X2.dat");
			else if(is==313)  myfile.open("X3.dat");
			else if(is==470)  myfile.open("X4.dat");
			else if(is==627) myfile.open("X5.dat");
			for(auto i=0; i<dd1; i++){
				myfile<<pgv->xc[pgv->Gn[i].ijk[0]-(1-pgv->nghost)]<<'\t';
				myfile<<pgv->yc[pgv->Gn[i].ijk[1]-(1-pgv->nghost)]<<'\t';
				myfile<<pgv->Gn[i].G<<std::endl;
			}
			myfile.close();
		}*/

   }

//	calcShapeError();
   plco->lit_finalize();
   pp->parallel_end();
   return 0;
}
void setVelocity(const int is){
   std::unique_ptr<global_variable> pgv(new global_variable);
	double periodT, tt, xpi, ypi, delta_t;
	std::vector<std::vector<double>> u_temp, v_temp, u2_temp, v2_temp;
	u_temp.resize(200);
	v_temp.resize(200);
	u2_temp.resize(200);
	v2_temp.resize(200);
	for(auto i=0; i<200; i++){
		u_temp[i].resize(200);
		v_temp[i].resize(200);
		u2_temp[i].resize(200);
		v2_temp[i].resize(200);
	}
	std::ifstream myfile1, myfile2, myfile3, myfile4;
	myfile1.open("u_velocity.dat");
	myfile2.open("v_velocity.dat");
	myfile3.open("u1_velocity.dat");
	myfile4.open("v1_velocity.dat");
	for(auto i=0; i<200; i++){
		for(auto j=0; j<200; j++){
			myfile1>>u_temp[i][j];
			myfile2>>v_temp[i][j];
			myfile3>>u2_temp[i][j];
			myfile4>>v2_temp[i][j];
		}
	}
	myfile1.close();
	myfile2.close();
	myfile3.close();
	myfile4.close();
//	delta_t = 3.90625e-3;

//	periodT = 8.0;
//	tt = (8.0/2048.0)*(static_cast<double>(is)+0.5);
   for(auto ibl=0; ibl<pgv->nbl; ibl++){
      pgv->setBlockPointers(ibl);
      for(auto ic=0; ic<pgv->b->NinT; ic++){
			xpi = pgv->Gn[ic].ijk[0]-(1-pgv->nghost)-3;
			ypi = pgv->Gn[ic].ijk[1]-(1-pgv->nghost)-3;
			if(pgv->Gn[ic].G >= 0.0){
				pgv->Gn[ic].V[0] = 1.0;//u_temp[xpi][ypi];
				pgv->Gn[ic].V[1] = 0.0;//v_temp[xpi][ypi];
			}
			else{
				pgv->Gn[ic].V[0] = 1.0;//u2_temp[xpi][ypi];
				pgv->Gn[ic].V[1] = 0.0;//v2_temp[xpi][ypi];
			}
			pgv->Gn[ic].V[2] = 0.0;
			xpi = pgv->pi*pgv->xc[pgv->Gn[ic].ijk[0]-(1-pgv->nghost)];
			ypi = pgv->pi*pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)];
//       pgv->Gn[ic].V[0] =  1.0;//-2.0*pow(sin(xpi),2)*sin(ypi)*cos(ypi)*cos(pgv->pi*tt/periodT);
//       pgv->Gn[ic].V[1] =  0.0;// 2.0*pow(sin(ypi),2)*sin(xpi)*cos(xpi)*cos(pgv->pi*tt/periodT);
//      pgv->Gn[ic].V[2] =  	0.0;
//	pgv->Gn[ic].V[0] =  0.5-(pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)]-pgv->xyzs_sg[1])/ \
										  (pgv->xyze_sg[1]-pgv->xyzs_sg[1]);
//	pgv->Gn[ic].V[1] = -0.5+(pgv->xc[pgv->Gn[ic].ijk[0]-(1-pgv->nghost)]-pgv->xyzs_sg[0])/ \
											(pgv->xyze_sg[0]-pgv->xyzs_sg[0]);
//			pgv->Gn[ic].V[2] = 0.0;
      }
   }
}

void calcShapeError(){
	std::unique_ptr<global_variable> pgv(new global_variable);
	std::unique_ptr<init> pinit(new init);
	std::unique_ptr<parallel> pp(new parallel);
	int n = 1000;
	double error, G00, G01, G10, G11, G00ex, G10ex, G01ex, G11ex, ax, ay, ddh, ddx, ddy;
	std::vector<double>xyz(3);
	double errorsum, G0, G1, G0ex, GG, Gex, G1ex;
	if(pp->myrank == 0)
		std::cout<<"*** Calculating shape error(this might take a while)***"<<std::endl;
	error = 0.0;
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinN; ic++){
			auto i = pgv->Gn[ic].ijk[0];
			auto j = pgv->Gn[ic].ijk[1];
			auto k = pgv->Gn[ic].ijk[2];
			G00 = pgv->Gn[ic].G;
			G10 = pgv->Gn[pgv->i2c[i+1-pgv->b->imino_][j-pgv->b->jmino_][k-pgv->b->kmino_]-1].G;
			G01 = pgv->Gn[pgv->i2c[i-pgv->b->imino_][j+1-pgv->b->jmino_][k-pgv->b->kmino_]-1].G;
			if(pgv->i2c[i+1-pgv->b->imino_][j+1-pgv->b->jmino_][k-pgv->b->kmino_]>0){
				G11 = pgv->Gn[pgv->i2c[i+1-pgv->b->imino_][j+1-pgv->b->jmino_][k-pgv->b->kmino_]-1].G;
			}
			else
				G11=0.5*(2.0*G10-pgv->Gn[pgv->i2c[i+1-pgv->b->imino_][j-1-pgv->b->jmino_]\
												[k-pgv->b->kmino_]-1].G);
			G00 = pgv->Gn[ic].G;
			G10 = pgv->Gn[pgv->i2c[i+1-pgv->b->imino_][j-pgv->b->jmino_][k-pgv->b->kmino_]-1].G;
			G01 = pgv->Gn[pgv->i2c[i-pgv->b->imino_][j+1-pgv->b->jmino_][k-pgv->b->kmino_]-1].G;
			if(pgv->i2c[i+1-pgv->b->imino_][j+1-pgv->b->jmino_][k-pgv->b->kmino_]>0){
				G11 = pgv->Gn[pgv->i2c[i+1-pgv->b->imino_][j+1-pgv->b->jmino_][k-pgv->b->kmino_]-1].G;
			}
			else
			G11=0.5*(2.0*G10-pgv->Gn[pgv->i2c[i+1-pgv->b->imino_][j-1-pgv->b->jmino_][k-pgv->b->kmino_]-1].G+\
						2.0*G01-pgv->Gn[pgv->i2c[i-1-pgv->b->imino_][j+1-pgv->b->jmino_][k-pgv->b->kmino_]-1].G);
			xyz[0] = pgv->xc[i-(1-pgv->nghost)];
			xyz[1] = pgv->yc[j-(1-pgv->nghost)];
			xyz[2] = 0.0;
			G00ex  = pinit->G_init_value(xyz);
			xyz[0]+= pgv->dxyz[0];
			G10ex  = pinit->G_init_value(xyz);
			xyz[1]+= pgv->dxyz[1];
			G11ex  = pinit->G_init_value(xyz);
			xyz[0] = pinit->G_init_value(xyz);
			G01ex  = pinit->G_init_value(xyz);
			ddh    = 1.0/static_cast<double>(n);
			ax     = -ddh;
			ddx    = ddh*pgv->dxyz[0];
			ddy    = ddh*pgv->dxyz[1];
			xyz[0] = pgv->xc[i-(1-pgv->nghost)] - ddx;
			xyz[1] = pgv->yc[j-(1-pgv->nghost)] - ddy;
			xyz[2] = 0.0;
			for(auto ii=0; ii<n; ii++){
				ax    += ddh;
				xyz[0]+= ddx;
				ay     = -ddh;
				xyz[1] = pgv->yc[j-(1-pgv->nghost)] - ddy;
				G0     = G00+ax*(G10-G00);
				G1     = G01+ax*(G11-G01);
				G0ex   = G00ex+ax*(G10ex-G00ex);
				G1ex   = G01ex+ax*(G11ex-G01ex);
				for(auto jj=0; jj<n; jj++){
					ay    += ddh;
					xyz[1]+= ddy;
					GG     = G0+ay*(G1-G0);
					Gex    = G0ex+ay*(G1ex-G0ex);
					if(GG*Gex < 0) error += 1.0;
				}
			}
		}
	}
	pp->parallel_sum(error, errorsum);
	if(pp->myrank == 0){
		errorsum = errorsum/(static_cast<double>(pow(n,2)))*(pgv->dxyz[0]*pgv->dxyz[1]);
		std::cout<<"======================================================="<<std::endl;
		std::cout<<"expected shape error    =   1.647020886230496E-002"<<std::endl;
		std::cout<<"calculated shape error  =   "<<errorsum<<std::endl;
		std::cout<<"======================================================="<<std::endl;
	}
}


