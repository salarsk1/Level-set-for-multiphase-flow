//Written by Salar Safarkhani

#include<iostream>
#include "parallel.h"
#include "init.h"
#include<cmath>
#include<limits>
#include<iostream>
#include "datag.h"
#include<cfloat>
void init::init_m_init(){
};
void init::init_geometry(){
double x0;
int i, j, k;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
for(auto i=0; i<pgv->xyze_sg.size(); i++){
	pgv->dxyz[i]   = (pgv->xyze_sg[i]-pgv->xyzs_sg[i])/static_cast<double>(pgv->ijkm_gl[i]);
	pgv->rdxyz[i]  = 1.0/pgv->dxyz[i];
	pgv->rdxyz2[i] = pow(pgv->rdxyz[i],2);
};
pgv->cell_volume = pgv->dxyz[0]*pgv->dxyz[1]*pgv->dxyz[2];
if(pgv->cylindrical){
	auto d1 = pgv->dxyz[0];
	auto d2 = pgv->dxyz[1];
	auto d3 = 0.5 * pgv->dxyz[1] * pgv->dxyz[2];
	pgv->dxyz_min = std::min(std::min(d1, d2), d3);
	pgv->kappa_max = 2.0/pgv->dxyz_min;
}
else{
	auto d1 = pgv->dxyz[0];
	auto d2 = pgv->dxyz[1];
	auto d3 = pgv->dxyz[2];
	pgv->dxyz_min  = std::min(std::min(d1, d2), d3);
	pgv->kappa_max = 1.0/pgv->dxyz_min; 
};
pgv->dxyz_min_2 = pow(pgv->dxyz_min,2);
pgv->kappa_min  = pgv->kappa_max;
for(auto i=0; i<pgv->ijkm_sg.size(); i++){
	pgv->ijkm_sg.at(i) = pgv->ijkm_gl.at(i)/pgv->ijkm_bl.at(i);
};
if((pgv->cylindrical && pgv->ijkm_sg[2] !=0) ||pgv->ijkm_sg[2]%2 == 0){
	pp->litError("init_geometry","need an even number of supergrid cells\
	  				  in theta.");
};
for(auto i=0; i<pgv->dxyz_sg.size(); i++){
	pgv->dxyz_sg[i]  = (pgv->xyze_sg[i]-pgv->xyzs_sg[i])/static_cast<double>(pgv->ijkm_sg[i]);
	pgv->rdxyz_sg[i] = 1.0/pgv->dxyz_sg[i];
};
if(pp->ierr == 0){
	pgv->xc.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[0]]);
	pgv->yc.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[1]]);
	pgv->zc.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[2]]);
	pgv->lxf.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[0]]);
	pgv->lyf.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[1]]);
	pgv->lzf.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[2]]);
	pgv->ldx.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[0]]);
	pgv->ldy.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[1]]);
	pgv->ldz.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[2]]);
	pgv->lrdx.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[0]]);
	pgv->lrdy.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[1]]);
	pgv->lrdz.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[2]]);
	pgv->ryc.reset(new double[2 * pgv->nghost + pgv->ijkm_gl[1]]);
};

if(pp->ierr != 0){
	pp->litError("init_geometry","Allocation error of x,y,z ... .");
};
auto d1 = 2 * pgv->nghost + pgv->ijkm_gl[0];
for(auto i=0; i< d1; i++){
	pgv->xc[i]  =  pgv->xyzs_sg[0] + (static_cast<double>(i)-0.5+static_cast<double>(1-pgv->nghost))*pgv->dxyz[0];
	pgv->	lxf[i] = pgv->xyzs_sg[0] + (static_cast<double>(i)-1.0+static_cast<double>(1-pgv->nghost))*pgv->dxyz[0];
	pgv->ldx[i]  = pgv->dxyz[0];	
	pgv->lrdx[i] = pgv->rdxyz[0];
};
d1 = 2 * pgv->nghost + pgv->ijkm_gl[1];
for(auto i=0; i<d1; i++){
	pgv->yc[i]  = pgv->xyzs_sg[1] + (static_cast<double>(i)-0.5+static_cast<double>(1-pgv->nghost))*pgv->dxyz[1];
	pgv->ryc[i] = 1.0/pgv->yc[i];
	pgv->lyf[i] = pgv->xyzs_sg[1] + (static_cast<double>(i)-1.0+static_cast<double>(1-pgv->nghost))*pgv->dxyz[1];
	pgv->ldy[i]  = pgv->dxyz[1];
	pgv->lrdy[i] = pgv->rdxyz[1];
};
d1 = 2 * pgv->nghost + pgv->ijkm_gl[2];
for(auto i=0; i<d1; i++){
	pgv->zc[i]  = pgv->xyzs_sg[2] + (static_cast<double>(i)-0.5+static_cast<double>(1-pgv->nghost))*pgv->dxyz[2];
	pgv->lzf[i] = pgv->xyzs_sg[2] + (static_cast<double>(i)-1.0+static_cast<double>(1-pgv->nghost))*pgv->dxyz[2];
	pgv->ldz[i]  = pgv->dxyz[2];
	pgv->lrdz[i] = pgv->rdxyz[2];
};
pgv->ucWeight11[0][0]  =  1.0/2772.0;
pgv->ucWeight11[1][0]  = -1.0/210.0;
pgv->ucWeight11[2][0]  =  5.0/168.0;
pgv->ucWeight11[3][0]  = -5.0/42.0;
pgv->ucWeight11[4][0]  =  5.0/14.0;
pgv->ucWeight11[5][0]  = -1.0;
pgv->ucWeight11[6][0]  =  1.0/6.0;
pgv->ucWeight11[7][0]  =  5.0/7.0;
pgv->ucWeight11[8][0]  = -5.0/28.0;
pgv->ucWeight11[9][0]  =  5.0/126.0;
pgv->ucWeight11[10][0] = -47.0/9240.0;
pgv->ucWeight11[11][0] = -2.0/2310.0;
for(auto i=0;i<12;i++){
	pgv->ucWeight11[i][1] = pgv->ucWeight11[i][0] * pgv->rdxyz[1];
	pgv->ucWeight11[i][2] = pgv->ucWeight11[i][0] * pgv->rdxyz[2];
	pgv->ucWeight11[i][0] = pgv->ucWeight11[i][0] * pgv->rdxyz[0];
};
pgv->ucWeight9[0][0]  =  1.0/630.0;
pgv->ucWeight9[1][0]  =  1.0/56.0;
pgv->ucWeight9[2][0]  = -2.0/21.0;
pgv->ucWeight9[3][0]  =  1.0/3.0;
pgv->ucWeight9[4][0]  = -1.0;
pgv->ucWeight9[5][0]  =  1.0/5.0;
pgv->ucWeight9[6][0]  =  2.0/3.0;
pgv->ucWeight9[7][0]  = -1.0/7.0;
pgv->ucWeight9[8][0]  =  1.0/42.0;
pgv->ucWeight9[9][0]  =  1.0/504;
for(auto i=0;i<10;i++){
	pgv->ucWeight9[i][1] = pgv->ucWeight9[i][0] * pgv->rdxyz[1];
	pgv->ucWeight9[i][2] = pgv->ucWeight9[i][0] * pgv->rdxyz[2];
	pgv->ucWeight9[i][0] = pgv->ucWeight9[i][0] * pgv->rdxyz[0];
};
pgv->ucWeight7[0][0]  =    3.0/420.0;
pgv->ucWeight7[1][0]  =  -28.0/420.0;
pgv->ucWeight7[2][0]  =  126.0/420.0;
pgv->ucWeight7[3][0]  = -420.0/420.0;
pgv->ucWeight7[4][0]  =  105.0/420.0;
pgv->ucWeight7[5][0]  =  252.0/420.0;
pgv->ucWeight7[6][0]  =  -42.0/420.0;
pgv->ucWeight7[7][0]  =    4.0/420.0;
for(auto i=0;i<8;i++){
	pgv->ucWeight7[i][1] = pgv->ucWeight7[i][0] * pgv->rdxyz[1];
	pgv->ucWeight7[i][2] = pgv->ucWeight7[i][0] * pgv->rdxyz[2];
	pgv->ucWeight7[i][0] = pgv->ucWeight7[i][0] * pgv->rdxyz[0];
};
pgv->ucWeight5[0][0]  =  -2.0/60.0;
pgv->ucWeight5[1][0]  =  15.0/60.0;
pgv->ucWeight5[2][0]  = -60.0/60.0;
pgv->ucWeight5[3][0]  =  20.0/60.0;
pgv->ucWeight5[4][0]  =  30.0/60.0;
pgv->ucWeight5[5][0]  =  -3.0/60.0;
for(auto i=0;i<6;i++){
	pgv->ucWeight5[i][1] = pgv->ucWeight5[i][0] * pgv->rdxyz[1];
	pgv->ucWeight5[i][2] = pgv->ucWeight5[i][0] * pgv->rdxyz[2];
	pgv->ucWeight5[i][0] = pgv->ucWeight5[i][0] * pgv->rdxyz[0];
};
pgv->ucWeight3[0][0]  =  1.0/6.0;
pgv->ucWeight3[1][0]  = -1.0;
pgv->ucWeight3[2][0]  =  0.5;
pgv->ucWeight3[3][0]  =  1.0/3.0;
for(auto i=0;i<4;i++){
	pgv->ucWeight3[i][1] = pgv->ucWeight3[i][0] * pgv->rdxyz[1];
	pgv->ucWeight3[i][2] = pgv->ucWeight3[i][0] * pgv->rdxyz[2];
	pgv->ucWeight3[i][0] = pgv->ucWeight3[i][0] * pgv->rdxyz[0];
};
};
double init::G_init_value(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
double GV;
if(pgv->init_shape == "plane")
	GV = G_plane(xyz);
else if(pgv->init_shape == "sheet")
	GV = G_sheet(xyz);
else if(pgv->init_shape == "random_sheet")
	GV = G_randomSheet(xyz);
else if(pgv->init_shape == "sine_column")
	GV = G_sine_column(xyz);
else if(pgv->init_shape == "sine")
	GV = G_sine(xyz);
else if(pgv->init_shape == "siney")
	GV = G_siney(xyz);
else if(pgv->init_shape == "cosine" || pgv->init_shape == "RMI")
	GV = G_cosine(xyz);
else if(pgv->init_shape == "cosine3D")
	GV = G_cosine3D(xyz);
else if(pgv->init_shape == "RT3D")
	GV = G_RT3D(xyz);
else if(pgv->init_shape == "cosiney")
	GV = G_cosiney(xyz);
else if(pgv->init_shape == "notched_circle")
	GV = G_notched_circle(xyz);
else if(pgv->init_shape == "circle")
	GV = G_circle(xyz);
else if(pgv->init_shape == "circle_xz")
	GV = G_circle(xyz);
else if(pgv->init_shape == "column")
	GV = G_column(xyz);
else if(pgv->init_shape == "column_cap_x")
	GV = G_column_cap_x(xyz);
else if(pgv->init_shape == "column_cap_y")
	GV = G_column_cap_y(xyz);
else if(pgv->init_shape == "column_cap_z")
	GV = G_column_cap_z(xyz);
else if(pgv->init_shape == "random_circle")
	GV = G_randomCircle(xyz);
else if(pgv->init_shape == "random_jet")
	GV = G_randomJet(xyz);
else if(pgv->init_shape == "random_drops")
	GV = G_randomDrops(xyz);
else if(pgv->init_shape == "function")
	GV = G_function2D(xyz);
else if(pgv->init_shape == "ring")
	GV = G_ring(xyz);
else if(pgv->init_shape == "rod")
	GV = G_rod(xyz);
else if(pgv->init_shape == "sphere")
	GV = G_sphere(xyz);
else if(pgv->init_shape == "sphere2")
	GV = G_sphere2(xyz);
else if(pgv->init_shape == "sphere_cyl")
	GV = G_sphere_cyl(xyz);
else if(pgv->init_shape == "disc_cyl")
	GV = G_disc_cyl(xyz);
else if(pgv->init_shape == "deformed_column")
	GV = G_deformed_column(xyz);
else if(pgv->init_shape == "deformed_sphere")
	GV = G_deformed_sphere(xyz);
else if(pgv->init_shape == "ellipse")
	GV = G_ellipse(xyz);
else if(pgv->init_shape == "milk_crown")
	GV = G_milk_crown(xyz);
else if(pgv->init_shape == "bursting_bubble")
	GV = G_bursting_bubble(xyz);
else if(pgv->init_shape == "column_pool")
	GV = G_column_pool(xyz);
else if(pgv->init_shape == "sphere_pool")
	GV = G_sphere_pool(xyz);
else if(pgv->init_shape == "bursting_bubble_3D")
	GV = G_bursting_bubble_3D(xyz);
else if(pgv->init_shape == "bursting_bubble_2D")
	GV = G_bursting_bubble_2D(xyz);
else if(pgv->init_shape == "bursting_bubble_rim_3D")
	GV = G_bursting_bubble_rim_3D(xyz);
else if(pgv->init_shape == "bursting_bubble_rim_2D")
	GV = G_bursting_bubble_rim_2D(xyz);
else if(pgv->init_shape == "ligament_2D")
	GV = G_ligament_2D(xyz);
else if(pgv->init_shape == "droplens")
	GV = G_droplens(xyz);
else if(pgv->init_shape == "Dam_break_2D")
	GV = G_Dam_break_2D(xyz);
else if(pgv->init_shape == "rayleigh")
	GV = G_rayleigh(xyz);
else if(pgv->init_shape == "rayleigh_axi")
	GV = G_rayleigh_axi(xyz);
else{
	std::cout<<"unknown initial G shape:"<< pgv->init_shape<<std::endl;
	pp->parallel_kill(0);
};
double m = std::min(pgv->G_max, GV);
return std::max(pgv->G_min, m);
};

double init::G_inject(const vector<double> &xyz, shape_t *shape){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
double GI;
switch(shape->code){
case(25):
	GI = G_from_sphere(xyz, shape->rdata);
break;
default:
	std::cout<<"unknown inject G shape code:"<<shape->code<<std::endl;
	pp->parallel_kill(0);
};
double m = std::min(pgv->G_max, GI);
return std::max(pgv->G_min, m);
};

double init::G_notched_circle(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
double c, bb, b1, b2, h1, h2;
c=pgv->init_radius-sqrt(pow(xyz.at(0)-pgv->init_center.at(0),2)+\
								pow(xyz.at(1)-pgv->init_center.at(1),2));

b1 = pgv->init_center.at(0) - 0.5 * pgv->init_width;
b2 = pgv->init_center.at(0) + 0.5 * pgv->init_width;
h1 = pgv->init_center.at(1) - pgv->init_radius * \
	  cos(asin(0.5*pgv->init_width/pgv->init_radius));
h2 = pgv->init_center.at(1) - pgv->init_radius+pgv->init_height;
if(0.0 <= c && xyz[0] <= b1 && xyz[1] <= h2){
	bb = b1 - xyz[0];
	return std::min(c, bb);
}
else if(0.0 <= c && b2 <= xyz[0] && h2 >= xyz[1]){
	bb = xyz[0] - b2;
	return std::min(c, bb);
}
else if(0.0 <= c && b1 <= xyz[0] && xyz[0] <= b2 && h2 <= xyz[1]){
	bb = xyz[1] - h2;
	return std::min(c, bb);
}
else if(0.0 <= c && xyz[0] <= b1 && h2 <= xyz[1]){
	bb = sqrt(pow(xyz.at(0)-b1,2)+pow(xyz.at(1)-h2,2));
	return std::min(c, bb);
}
else if(0.0 <= c && b2 <= xyz[0] && h2 <= xyz[1]){
	bb = sqrt(pow(xyz.at(0)-b2,2)+pow(xyz.at(1)-h2,2));
	return std::min(c, bb);
}
else if(b1<=xyz[0] && xyz[0]<=b2 && xyz[1]<=h2 && h1<=xyz[1]){
	auto d1 = fabs(xyz[0] - b1);
	auto d2 = fabs(xyz[0] - b2);
	auto d3 = fabs(xyz[1] - h2);
	return (-std::min(std::min(d1,d2),d3));
}
else if(b1<=xyz[0] && xyz[0]<=b2 && xyz[1] <=h1){
	auto d1 = sqrt(pow(xyz[0]-b1,2)+pow(xyz[1]-h1,2));
	auto d2 = sqrt(pow(xyz[0]-b2,2)+pow(xyz[1]-h1,2));
	return -std::min(d1, d2);
}
else{
	return c;
};
};

double init::G_circle(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
return pgv->init_radius-sqrt(pow(xyz[0]-pgv->init_center[0],2)+\
									  pow(xyz[1]-pgv->init_center[1],2));
};

double init::G_circle_xz(const vector<double> &xyz){
std::array<double, 3> xz;
std::unique_ptr<global_variable>pgv(new global_variable);
for(auto i=0; i<xz.size(); i++) xz[i] = xyz[i];
return pgv->init_radius-sqrt(pow(xz[0]-pgv->init_center[0],2)+\
									  pow(xz[1]-pgv->init_center[1],2)+\
									  pow(xz[2]-pgv->init_center[2],2));
};

double init::G_column(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
return pgv->init_radius-sqrt(pow(xyz[1]-pgv->init_center[1],2)+\
									  pow(xyz[2]-pgv->init_center[2],2));
};

double init::G_column_cap_x(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
if(xyz[0]-pgv->init_center[0] < pgv->init_length){
	return pgv->init_radius-sqrt(pow(xyz[1]-pgv->init_center[1],2)+\
										  pow(xyz[2]-pgv->init_center[2],2));
}
else{
	return pgv->init_radius-sqrt(pow(xyz[1]-pgv->init_center[1],2)+\
										  pow(xyz[2]-pgv->init_center[2],2)+
 						   pow(xyz[0]-pgv->init_center[0]-pgv->init_length,2));
};
};

double init::G_column_cap_y(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
if(xyz[1] - pgv->init_center[1] < pgv->init_length){
	return pgv->init_radius-sqrt(pow(xyz[0]-pgv->init_center[0],2)+\
										  pow(xyz[2]-pgv->init_center[2],2));
}
else{
	return pgv->init_radius-sqrt(pow(xyz[0]-pgv->init_center[0],2)+\
										  pow(xyz[2]-pgv->init_center[2],2)+\
						   pow(xyz[1]-pgv->init_center[1]-pgv->init_length,2));
}
};
double init::G_column_cap_z(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
if(xyz[2]-pgv->init_center[2]<pgv->init_length){
	return pgv->init_radius-sqrt(pow(xyz[0]-pgv->init_center[0],2)+\
										  pow(xyz[1]-pgv->init_center[1],2));
}
else{
	return pgv->init_radius-sqrt(pow(xyz[0]-pgv->init_center[0],2)+\
										  pow(xyz[1]-pgv->init_center[1],2)+\
					    pow(xyz[2]-pgv->init_center[2]-pgv->init_length,2));
};
};
double init::G_column_pool(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
double G1, G2;
G1 = pgv->init_radius-sqrt(pow(xyz[0]-pgv->init_center[0],2)+\
									pow(xyz[1]-pgv->init_center[1],2));
G2 = pgv->init_height-(xyz[1]-pgv->xyzs_sg[1]);
if(0.0 <= G2) 
	return G2;
else if(0.0 < G1) 
	return G1;
else
	return std::max(G1, G2);
};

double init::G_sphere_pool(const vector<double> &xyz){
double G1, G2;
std::unique_ptr<global_variable>pgv(new global_variable);
G1 = pgv->init_radius-sqrt(pow(xyz[0]-pgv->init_center[0],2)+\
									pow(xyz[1]-pgv->init_center[1],2)+\
									pow(xyz[2]-pgv->init_center[2],2));
G2 = pgv->init_height-(xyz[2]-pgv->xyzs_sg[2]);
if(0.0 <= G2)
	return G2;
else if(0.0 < G1)
	return G1;
else
	return std::max(G1, G2);
};
double init::G_randomCircle(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
double rand1 = static_cast<double>(rand())/RAND_MAX;
return pgv->init_radius-sqrt(pow(xyz[0]-pgv->init_center[0],2)+\
									  pow(xyz[1]-pgv->init_center[1],2))+\
									  2.0*(rand1-0.5)*pgv->init_amplitude;
};

double init::G_randomJet(const vector<double> &xyz){
double BoxLength;
double z_perturb, the_perturb, theta;
std::unique_ptr<global_variable>pgv(new global_variable);
BoxLength = pgv->xyze_sg[2] - pgv->xyzs_sg[2];
z_perturb = 0.0;
the_perturb = 0.0;
theta = arctan(xyz[0], xyz[1]);
for(auto i=0; i<pgv->n1_random; i++){
	z_perturb = z_perturb + pgv->a1_random[i]*cos(2.0*pgv->pi*\
									 static_cast<double>(i)/BoxLength*xyz[2]+\
									 +pgv->phi1_random[i]);
									 
};
for(auto i=0; i<pgv->n2_random; i++){
	z_perturb = z_perturb + pgv->a2_random[i]*cos(static_cast<double>(i)*\
									theta+pgv->phi2_random[i]);
};
return (pgv->init_radius-sqrt(pow(xyz[0]-pgv->init_center[0],2)+\
										pow(xyz[1]-pgv->init_center[1],2)))+\
								z_perturb*the_perturb*pgv->init_radius*2.0;
};
double init::G_randomDrops(const vector<double> &xyz){
static int nd;
static bool firstCall;
double rand1, rand2, rand3, gg, g, x, y, z, v;
std::array<int, 3> i, ijksg;
bool done;
char c1;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
std::ifstream input;
std::ofstream output;
if(firstCall){
	firstCall = false;
	if(pgv->init_mod == 1){
		if(pp->myrank == 0){
			input.open("random_drops.kal");
			input>>nd;
			if(nd != pgv->init_number){
				std::cout<<"WARNING: random_drops uses"<<nd<<"drops as read \
					from random_drops.kal, not"<<pgv->init_number<<std::endl;
			};
			if(pp->ierr==0){
				xd = new double*[3];
				for(auto i=0; i<3; i++){
					xd[i] = new double[nd];
				};
				dd.reset(new double[nd]);
				for(auto i=0; i<nd; i++){
					input>>xd[0][i]>>c1>>xd[1][i]>>c1>>xd[2][i]>>c1>>v>>c1\
						  >>dd[i];
				};
				input.close();
			};
			pp->ierr = MPI_Bcast(&nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
		}
		else{
			pp->ierr = MPI_Bcast(&nd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
			if(pp->ierr==0){
				xd = new double*[3];
				for(auto i=0;i<3;i++){
					xd[i] = new double[nd];
				};
				dd.reset(new double[nd]);
			};
		};
		pp->ierr = MPI_Bcast(xd, 2*nd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		pp->ierr = MPI_Bcast(dd.get(),   nd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	else{
		nd = pgv->init_number;
		if(pp->ierr==0){
			xd = new double*[3];
			for(auto i=0; i<3; i++){
				xd[i] = new double[nd];
			};
			dd.reset(new double[nd]);
		};
		if(pp->myrank == 0){
			for(auto i=0;i<nd;i++){
				rand1 = rand()/RAND_MAX;
				rand2 = rand()/RAND_MAX;
				rand3 = rand()/RAND_MAX;
				xd[0][i] = pgv->init_xyz_min[0]+rand1*(pgv->init_xyz_max[0]-\
							  										pgv->init_xyz_min[0]);
				xd[1][i] = pgv->init_xyz_min[1]+rand1*(pgv->init_xyz_max[1]-\
							  										pgv->init_xyz_min[1]);
				xd[2][i] = pgv->init_xyz_min[2]+rand1*(pgv->init_xyz_max[2]-\
							  										pgv->init_xyz_min[2]);
				rand1 = rand()/RAND_MAX;
				dd[i] = pgv->init_radius_min+rand1*(pgv->init_radius_max-\
																pgv->init_radius_min);
			};
			output.open("random_drops.kal");
			output<<nd;
			for(auto i=0;i<nd;i++){
				output<<xd[0][i]<<'\t'<<xd[1][i]<<'\t'<<xd[2][i]<<'\t'<<\
										4.0/3.0*pgv->pi*pow(dd[i],3)<<'\t'<<dd[i];
			};
			output.close();
		};
		pp->ierr = MPI_Bcast(&xd, 3*nd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		pp->ierr = MPI_Bcast(&xd,   nd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	};
};
gg = -DBL_MAX;
for(auto i=0;i<nd; i++){
	g = dd[i]-sqrt(pow(xyz[0]-xd[0][i],2)+pow(xyz[0]-xd[0][i],2)+\
													  pow(xyz[0]-xd[0][i],2));
	if(abs(g) < abs(gg)){
		if(!(0.0 < gg && g < 0.0)) gg = g;
	}
	else if(gg < 0.0 && 0.0 < g){
		 gg = g;
	};
};
return gg;
};
	
double init::G_ring(const vector<double> &xyz){
double G1, G2, Rx;
std::unique_ptr<global_variable>pgv(new global_variable);
Rx = sqrt(pow(xyz[0]-pgv->init_center[0],2)+pow(xyz[1]-pgv->init_center[1],2));
if(Rx < pgv->init_radius-0.5*pgv->init_width){
	return (Rx - pgv->init_radius+0.5*pgv->init_width);
}
else if(Rx>pgv->init_radius+0.5*pgv->init_width-Rx){
	return (pgv->init_radius+0.5*pgv->init_width-Rx);
}
else{
	auto d1 = Rx-pgv->init_radius+0.5*pgv->init_width;
	auto d2 = pgv->init_radius+0.5*pgv->init_width-Rx;
	return (std::min(d1, d2));
};
};

double init::G_rod(const vector<double> &xyz){
double G1, G2, Rx, alphalt, ln, lt;
std::array<double, 3> tang, norm, center;
std::unique_ptr<global_variable> pgv(new global_variable);
tang[0] = cos(pgv->init_angle);
tang[1] = sin(pgv->init_angle);
tang[2] = 0.0;
norm[0] = -tang[1];
norm[1] =  tang[0];
norm[2] =  0.0;
lt = (xyz[0]-pgv->init_center[0])*tang[0]+\
	  (xyz[1]-pgv->init_center[1])*tang[1];
ln = (xyz[0]-pgv->init_center[0])*norm[0]+\
	  (xyz[1]-pgv->init_center[1])*norm[1];
if(abs(lt) <= 0.5*pgv->init_length){
	return (0.5*pgv->init_width-abs(ln));
}
else{
	if(0.0 < lt){
		for(auto i=0; i<3; i++){
			center[i] = pgv->init_center[i]+0.5*pgv->init_length*tang[i];
		};
	}
	else{
		for(auto i=0; i<3; i++){
			center[i] = pgv->init_center[i]-0.5*pgv->init_length*tang[i];
		};
	};
	return(0.5*pgv->init_width-sqrt(pow(xyz[0]-center[0],2)+\
											  pow(xyz[1]-center[1],2)));
};
};

inline double init::G_from_sphere(const vector<double> &xyz, const std::array<double,10> &rdata){
	return(0.5*rdata[3]-sqrt(pow(xyz[0]-rdata[0],2)+pow(xyz[1]-rdata[1],2)+\
									 pow(xyz[2]-rdata[2],2)));
};

inline double init::G_sphere(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
double d1 = 0;
for(auto i=0; i<3; i++){
	d1 = d1 + pow(xyz[i] - pgv->init_center[i],2);
};
return(pgv->init_radius-sqrt(d1));
};

double init::G_sphere2(const vector<double> &xyz){
double G1, G2;
double d1, d2;
std::unique_ptr<global_variable>pgv(new global_variable);
for(auto i=0; i<3; i++){
	d1 = d1 + pow(xyz[i] - pgv->init_center[i],2);
};
for(auto i=0; i<3; i++){
	d2 = d2 + pow(xyz[i] - pgv->init_center2[i],2);
};
G1 = pgv->init_radius-sqrt(d1);
G2 = pgv->init_radius-sqrt(d2);
if(G1 < 0.0 && G2 < 0.0){
	return(-std::min(abs(G1), abs(G2)));
}
else if(0.0 <= G1 && 0.0 <= G2){
	return(std::min(G1, G2));
}
else{
	return(std::max(G1, G2));
};
};

inline double init::G_sphere_cyl(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
return( pgv->init_radius-sqrt(pow(xyz[0]-pgv->init_center[0],2)+\
										pow(xyz[1]*cos(xyz[2])-pgv->init_center[1],2)+\
										pow(xyz[1]*sin(xyz[2])-pgv->init_center[2],2)));
};

inline double init::G_disc_cyl(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
double G1;
G1 = 0.5*pgv->init_length-sqrt(pow(xyz[0]-pgv->init_center[0],2));
if(0.0 < G1){
	double d1 = xyz[1]*cos(xyz[2])-pgv->init_center[1];
	double d2 = xyz[1]*sin(xyz[2])-pgv->init_center[2];
	return 0.5*pgv->init_length*(pgv->init_radius-sqrt(pow(d1,2)+pow(d2,2)));
}
else{
	return G1;
};
};

inline double init::G_plane(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
return(pgv->init_center[pgv->init_mod-1]-xyz[pgv->init_mod-1]);
};

double inline init::G_sine(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
return(pgv->init_center[1]-xyz[1]+pgv->init_amplitude*\
    sin(2.0*pgv->pi*(xyz[0]-pgv->init_center[0])/pgv->init_wavelength));
};

inline double init::G_siney(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
return((pgv->init_center[0]-xyz[0])+pgv->init_amplitude*\
    sin(2.0*pgv->pi*(xyz[1]-pgv->init_center[1])/pgv->init_wavelength));
};

inline double init::G_sine_column(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
return(-abs(pgv->init_center[1]-xyz[1])+pgv->init_radius+pgv->\
init_amplitude*sin(2.0*pgv->pi*(xyz[0]-pgv->init_center[0])/pgv->\
						 init_wavelength));
};

inline double init::G_siney_column(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
return(-abs(pgv->init_center[0]-xyz[0])+pgv->init_radius+pgv->\
init_amplitude*sin(2.0*pgv->pi*(xyz[1]-pgv->init_center[1])/pgv->\
						 init_wavelength));
}; 

inline double init::G_cosine3D(const vector<double> &xyz){
double r;
std::unique_ptr<global_variable>pgv(new global_variable);
auto d1 = sqrt(pow(xyz[0]-pgv->init_center[0],2)+pow(xyz[2]-pgv->\
					init_center[2],2)/pgv->init_wavelength);
r = std::min(0.5, d1);
return (-(pgv->init_center[1]-xyz[1]-pgv->init_amplitude*cos(2.0*pgv->pi*r)));
};
inline double init::G_RT3D(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
double k = 2.0*pgv->pi/pgv->init_wavelength;
return(xyz[1]-pgv->init_center[1]-pgv->init_amplitude*(cos(k*xyz[0])+cos(k*xyz[2])));
};

inline double init::G_cosine(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
return -(pgv->init_center[1]-xyz[1]+pgv->init_amplitude*\
		  cos(2.0*pgv->pi*(xyz[0]-pgv->init_center[0])/pgv->init_wavelength));
};

inline double init::G_cosiney(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
return -(pgv->init_center[0]-xyz[0]+pgv->init_amplitude*\
		   cos(2.0*pgv->pi*(xyz[1]-pgv->init_center[1])/pgv->init_wavelength));
};

inline double init::G_sheet(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
double d1 = 0.0;
for(auto i=0; i<3; i++){
	d1= d1 + (xyz[i]-pgv->init_center[i])*pgv->init_normal[i];
};
return (0.5*pgv->init_width-abs(d1));
};

double init::G_randomSheet(const vector<double> &xyz){
double dn, dA;
std::unique_ptr<global_variable>pgv(new global_variable);
dn = 0.0;
for(auto i=0; i<3; i++){
	dn = dn + (xyz[i]-pgv->init_center[i])*pgv->init_normal[i];
};
dA = rand()/RAND_MAX;
dA = pgv->init_amplitude*dA;
return 0.5*(pgv->init_width+dA)-abs(dn);
};

double init::G_deformed_column(const vector<double> &xyz){
double angle, costheta, theta0;
std::array<double, 3> xyzn;
std::unique_ptr<global_variable>pgv(new global_variable);
angle = 0.0;
theta0 = angle/180.0*pgv->pi;
xyzn[0]  =  (xyz[0]-pgv->init_center[0])*cos(theta0)+(xyz[1]-pgv->\
			    init_center[1])*sin(theta0)+pgv->init_center[0];
xyzn[1]  = -(xyz[0]-pgv->init_center[0])*sin(theta0)+(xyz[1]-pgv->\
			    init_center[1])*cos(theta0)+pgv->init_center[1];
costheta = acos(xyzn[0]-pgv->init_center[0])/sqrt(pow(xyzn[0]-pgv->\
					 init_center[0],2)+pow(xyzn[1]-pgv->init_center[1],2));
return pgv->init_radius-pgv->init_amplitude*cos(static_cast<double>(pgv->\
		 init_mod)*costheta)-sqrt(pow(xyzn[0]-pgv->init_center[0],2)+\
		 pow(xyzn[1]-pgv->init_center[1],2));
};

double init::G_ellipse(const vector<double> &xyz){
std::array<double, 3> xyzn;
std::unique_ptr<global_variable>pgv(new global_variable);
for(auto i=0; i<2; i++) xyzn[i] = xyz[i] -pgv->init_center[i];
return 1.0-(pow(xyzn[0]/pgv->init_a,2)+pow(xyzn[1]/pgv->init_b,2));
};

double init::G_deformed_sphere(const vector<double> &xyz){
double phi, theta, Sn, angle, x;
std::array<double, 3> xyzn;
std::unique_ptr<global_variable>pgv(new global_variable);
for(auto i=0; i<3; i++) xyzn[i] = xyz[i] - pgv->init_center[i];
phi = arctan(xyzn[1],xyzn[2]);
theta = arctan(xyzn[0], sqrt(pow(xyzn[1],2)+pow(xyzn[2],2)));
if(pgv->init_mod==2){
	x = cos(theta);
	Sn = 0.5*(3.0*pow(x,2)-1.0);
};
Sn = Sn*pgv->init_amplitude;
double d1 = 0.0;
for(auto i:xyzn) d1 = d1 + pow(i,2);
return pgv->init_radius+Sn-sqrt(d1);
};

inline double init::G_milk_crown(const vector<double> &xyz){
double G1, G2;
std::unique_ptr<global_variable>pgv(new global_variable);
G1 = 0.0;
for(auto i=0; i<3; i++){
	G1 = G1 + pow(xyz[i]-pgv->init_center[i],2);
};
G1 = pgv->init_radius - sqrt(G1);
G2 = pgv->init_height - xyz[2];
return std::max(G1, G2);
};

double init::G_bursting_bubble(const vector<double> &xyz){
double G1, G2;
std::unique_ptr<global_variable>pgv(new global_variable);
G1 = 0.0;
for(auto i=0; i<2; i++){
	G1 = G1 + pow(xyz[i]-pgv->init_center[i],2);
};
G1 = -pgv->init_radius+sqrt(G1);
G2 = pgv->init_height-xyz[1];
return std::min(G1, G2);
};

double init::G_bursting_bubble_3D(const vector<double> &xyz){
double G1, G2;
std::unique_ptr<global_variable>pgv(new global_variable);
if(pgv->init_height < xyz[2]){
	return pgv->init_height-xyz[2];
}
else{
	G1 = 0.0;
	for(auto i=0; i<3; i++){
		G1 = G1 + pow(xyz[i]-pgv->init_center[i],2);
	};
	G1 = -pgv->init_radius+sqrt(G1);
	G2 = pgv->init_height - xyz[2];
	return std::min(G1, G2);
};
};

double init::G_bursting_bubble_2D(const vector<double> &xyz){
double G1, G2;
std::unique_ptr<global_variable>pgv(new global_variable);
if(pgv->init_height < xyz[1]){
	return pgv->init_height-xyz[1];
}
else{
	G1 = 0.0;
	for(auto i=0; i<2; i++){
		G1 = G1 + pow(xyz[i]-pgv->init_center[i],2);
	};
	G1 = -pgv->init_radius + sqrt(G1);
	G2 =  pgv->init_height - xyz[1];
	return std::min(G1, G2);
};
};

double init::G_bursting_bubble_rim_3D(const vector<double> &xyz){
std::array<double, 3> bubbleCenter, remCenter;
double r, zz, GG, rCr, rCz, bCr, bCz, rr, alpha, beta, G;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
r=0.0;
for(auto i=0;i<2; i++){
	r = r + pow(xyz[i]-pgv->init_center[i],2);
};
r = sqrt(r);
zz = xyz[2];
G = DBL_MAX;
GG = pgv->init_height - zz;
if(r > pgv->init_hole_radius){
	if(0.0<=GG){
		G = std::min(abs(GG), abs(G));
	}
	else{
		G = -std::min(abs(GG), abs(G));
	};
};
rCr = pgv->init_hole_radius;
rCz = pgv->init_height - pgv->init_rim_radius;
rr  = sqrt(pow(r-rCr,2)+pow(zz-rCz,2));
GG  = pgv->init_rim_radius - rr;
if(r<=rCr){
	if(zz>rCz){
		if(abs(GG) < abs(G)){
			G = GG;
		};		
	}
	else{
		alpha = acos((rCz-zz)/rr);
		beta  = asin(pgv->init_hole_radius/(pgv->init_hole_radius+\
														pgv->init_rim_radius));
		if(alpha>beta){
			if(abs(GG) < abs(G)) G = GG;
		};
	};
};
bCr = 0.0;
beta=asin(pgv->init_hole_radius/(pgv->init_radius+pgv->init_rim_radius));
bCz = pgv->init_height-pgv->init_rim_radius-cos(beta)*(pgv->init_radius+\
																	pgv->init_rim_radius);
rr = sqrt(pow(r-bCr,2)+pow(zz-bCz,2));
GG = rr - pgv->init_radius;
if(zz <= bCz){
	if(abs(GG) < abs(G)) G = GG;
}
else{
	alpha = acos((zz-bCz)/rr);
 beta=asin(pgv->init_hole_radius/(pgv->init_radius+pgv->init_rim_radius));
	if(beta <= alpha){
		if(abs(GG) < abs(G)) G = GG;
	};
};
if(G>1.0e20){
	std::cout<< pgv->clit<<"ERROR setting init G for bursting bubble with\
									rim in 3D at xyz="<<xyz[0]<<'\t'<<xyz[1]<<\
									'\t'<<xyz[2];
	pp->parallel_kill(34);
};
};
inline double init::G_rayleigh_axi(const vector<double> &xyz){
double R;
std::unique_ptr<global_variable>pgv(new global_variable);
R=pgv->init_amplitude*cos((xyz[0]-pgv->init_center[0])/(pgv->xyze_sg[0]-\
								pgv->xyzs_sg[0])*2.0*pgv->pi*static_cast<double>\
								(pgv->init_mod));
return R+pgv->init_radius-sqrt(pow(xyz[1]-pgv->init_center[1],2));
};

inline double init::G_rayleigh(const vector<double> &xyz){
double R;
std::unique_ptr<global_variable>pgv(new global_variable);
R=pgv->init_amplitude*cos((xyz[0]-pgv->init_center[0])/(pgv->xyze_sg[0]-\
								pgv->xyzs_sg[0])*2.0*pgv->pi*static_cast<double>\
								(pgv->init_mod));
return R+pgv->init_radius-sqrt(pow(xyz[1]-pgv->init_center[1],2)+\
		 pow(xyz[2]-pgv->init_center[2],2));
};
double init::G_bursting_bubble_rim_2D(const vector<double> &xyz){
std::array<double, 3> bubbleCenter, rimCenter;
double r, zz, GG, rCr, rCz, bCr, bCz, rr, alpha, beta, G;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
r  = xyz[0];
zz = xyz[2];
G = DBL_MAX;
GG = pgv->init_height - zz;
if(r>pgv->init_hole_radius){
	if(0 <= GG){
		G = std::min(abs(GG), abs(G));
	}
	else{
		G = -std::min(abs(GG), abs(G));
	};
};
rCr = pgv->init_hole_radius;
rCz = pgv->init_height-pgv->init_rim_radius;
rr  = sqrt(pow(r-rCr,2)+pow(zz-rCz,2));
GG = pgv->init_rim_radius - rr;
if(r <= rCr){
	if(zz > rCz){
		if(abs(GG) < abs(G)) G = GG;
	}
	else{
		alpha = acos((rCz-zz)/rr);
 beta=asin(pgv->init_hole_radius/(pgv->init_radius+pgv->init_rim_radius));
		if(alpha>beta){
			if(abs(GG) < abs(G)) G = GG;
		};
	};
};
bCr = 0.0;
beta=asin(pgv->init_hole_radius/(pgv->init_radius+pgv->init_rim_radius));
bCz=pgv->init_height-pgv->init_rim_radius-cos(beta)*(pgv->init_radius+\
			pgv->init_rim_radius);
rr = sqrt(pow(r-bCr,2)+pow(zz-bCz,2));
GG = rr - pgv->init_radius;
if(zz <= bCz){
	if(abs(GG) < abs(G)) G = GG;
}
else{
	alpha=acos((zz-bCz)/rr);
 beta=asin(pgv->init_hole_radius/(pgv->init_radius+pgv->init_rim_radius));
	if(beta <= alpha){
		if(abs(GG) < abs(G)) G = GG;
	};
};
if(G > 1.0e20){
	std::cout<<pgv->clit<<"ERROR setting init G for bursting bubble with\
								  rim in 3D at xyz="<<xyz[0]<<'\t'<<xyz[1]<<\
								  '\t'<<xyz[2]<<std::endl;
	pp->parallel_kill(34);
};
return G;
};

double init::G_ligament_2D(const vector<double> &xyz){
std::array<double, 3> bubbleCenter, rimCenter;
std::unique_ptr<global_variable>pgv(new global_variable);
if(xyz[0]<-pgv->init_length/2.0){
	return 0.5*pgv->init_width-sqrt(pow(xyz[0]+0.5*pgv->init_length,2)+\
											  pow(xyz[2],2));
}
else if(xyz[0] < pgv->init_length/2.0){
	return 0.5*pgv->init_width-abs(xyz[1]);
}
else{
 	return 0.5*pgv->init_width-sqrt(pow(xyz[0]-0.5*pgv->init_length,2)+\
											  pow(xyz[2],2));
};
};
double init::G_droplens(const vector<double> &xyz){
double x1, x2, bR;
std::unique_ptr<global_variable>pgv(new global_variable);
x1 = 0.5*pgv->init_height+sqrt(std::max(0.0, pow(pgv->init_radius2,2)-\
	  													   pow(pgv->init_radius,2)));
x2 = x1 - pgv->init_height;
if(abs(xyz[0]) < 0.5*pgv->init_height){
	return pgv->init_radius-sqrt(pow(xyz[1],2)+pow(xyz[2],2));
}
else if(xyz[0] > 0.5*pgv->init_height){
	return pgv->init_radius2-sqrt(pow(xyz[0]-x1,2)+pow(xyz[1],2)+\
											pow(xyz[2],2));
}
else{
	return pgv->init_radius2-sqrt(pow(xyz[0]-x2,2)+pow(xyz[1],2)+\
											pow(xyz[2],2));
};
};

double init::G_Dam_break_2D(const vector<double> &xyz){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
if(xyz[0] < pgv->init_width && xyz[1] < pgv->init_height){
	return  std::min(pgv->init_width-xyz[0], pgv->init_height-xyz[1]);
}
else if(xyz[0] > pgv->init_width && xyz[1] >pgv->init_height){
	return -sqrt(pow(xyz[0]-pgv->init_width,2)+pow(xyz[1]-pgv->init_height,2));
}
else if(xyz[0] > pgv->init_width){
	return -xyz[0] + pgv->init_width;
}
else{
	return -xyz[1] + pgv->init_height;
};
};

inline double init::G_function2D(const vector<double> &xyz){
return -(pow(xyz[0]/3.0,2)+pow(xyz[1]/2.0,2)-1.0);
};

double init::arctan(double ddx, double ddy){
double alpha;
std::unique_ptr<global_variable>pgv(new global_variable);
if(abs(ddx)+abs(ddy) < 1.0e-9){
	alpha = 0.0;
}
else{
	alpha = atan(ddy/ddx);
};
if(ddx <= 0.0){
	alpha = pgv->pi+alpha;
}
else if(ddy <=0.0 && 0.0 < ddx){
	alpha = 2.0*pgv->pi+alpha;
};
};
