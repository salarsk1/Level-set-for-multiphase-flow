//Written by Salar Safarkhani

#include<cmath>
#include<algorithm>
#include<fstream>
#include<cfloat>
#include"toolbox.h"
void toolbox::toolbox_m_init(){
std::unique_ptr<global_variable> pgv(new global_variable);
width_d = pgv->delta_width_factor*(*std::min_element(pgv->dxyz.begin(),pgv->dxyz.end()));
width_H = width_d;
r2pi = 0.5/pgv->pi;
pi_wh = pgv->pi/width_H;
r2wh = 0.5/width_H;
pi_wd = pgv->pi/width_d;
r2wd = 0.5/width_d;
r13 = 1.0/3.0;
}
double toolbox::delta(const double G){
if(fabs(G) >= width_d) return 0.0;
else return r2wd*(1.0+cos(pi_wd*G));
}
double toolbox::delta_L(const double G){
if(fabs(G) > width_d) return 0.0;
else return (1.0-fabs(G/width_d))/width_d;
}
bool toolbox::cell_contains_front(const Gnode_t &Gnh, const int ic, block_t *bb){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);
if(ic <= bb->NinT){
	auto ii = Gnh.ijk[0];
	auto jj = Gnh.ijk[1];
	auto kk = Gnh.ijk[2];
	if(kk==1 && ii==1) return true;
	else if(kk == 1 && ii == bb->ijkm[0]) return true;
	else if(kk == 1 && jj == 1) return true;
	else if(kk == 1 && jj == bb->ijkm[1]) return true;
	else if(kk == bb->ijkm[2] && ii == 1) return true;
	else if(kk == bb->ijkm[2] && ii == bb->ijkm[0]) return true;
	else if(kk == bb->ijkm[2] && jj == 1) return true;
	else if(kk == bb->ijkm[2] && jj == bb->ijkm[1]) return true;
	else if(ii == 1 && jj == 1) return true;
	else if(ii == bb->ijkm[0] && jj == 1) return true;
	else if(ii == bb->ijkm[0] && jj == bb->ijkm[1]) return true;
	else if(ii == 1 && jj == bb->ijkm[1]) return true;
	std::exit(1);
	if(Gnh.G >= 0.0){
		for(auto i=Gnh.ijk[0]-bb->imino_-1; i<Gnh.ijk[0]-bb->imino_+2; i++){
			for(auto j=Gnh.ijk[1]-bb->jmino_-1; j<Gnh.ijk[1]-bb->jmino_+2; j++){
				for(auto k=Gnh.ijk[2]-bb->kmino_-1; j<Gnh.ijk[2]-bb->kmino_+2; k++){
					if(bb->Gnodes[bb->ijk2ic[i][j][k]-1].G <= 0.0) return true;
				}
			}
		}
		return false;
	}
	else{
		for(auto i=Gnh.ijk[0]-1-pgv->b->imino_; i<=Gnh.ijk[0]+1-pgv->b->imino_; i++){
			for(auto j=Gnh.ijk[1]-1-pgv->b->jmino_; j<=Gnh.ijk[1]+1-pgv->b->jmino_; j++){
				for(auto k=Gnh.ijk[2]-1-pgv->b->kmino_; j<=Gnh.ijk[2]+1-pgv->b->kmino_; k++){
					if(bb->Gnodes[bb->ijk2ic[i][j][k]-1].G <= 0.0) return true;
				}
			}
		}
		return false;
	}		
}
else return false;
};
double toolbox::volume_delft(const double G, const std::vector<double> &gradG){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

const double COMPUTE_PSI_EPS = 1.0e-8;
double Gm, Dxi, Deta, Dzeta, A, B, C, D, E;
double Psi;
std::array<double, 3> abs_gradG;
Gm = -fabs(G);
for(auto i=0; i<3; i++) abs_gradG[i] = fabs(gradG[i]);
Dxi = *std::max_element(abs_gradG.begin(), abs_gradG.end());
Dzeta = *std::min_element(abs_gradG.begin(), abs_gradG.end());
Deta = std::accumulate(abs_gradG.begin(), abs_gradG.end(), 0.0) - Dxi - Dzeta;
A = std::max(Gm + 0.5*( Dxi + Deta + Dzeta),0.0);
B = std::max(Gm + 0.5*( Dxi + Deta - Dzeta),0.0);
C = std::max(Gm + 0.5*( Dxi - Deta + Dzeta),0.0);
D = std::max(Gm + 0.5*(-Dxi + Deta + Dzeta),0.0);
E = std::max(Gm + 0.5*( Dxi - Deta - Dzeta),0.0);
if(Dxi>COMPUTE_PSI_EPS){
	if(Deta>COMPUTE_PSI_EPS){
		if(Dzeta>COMPUTE_PSI_EPS)
			Psi = (pow(A,3)-pow(B,3)-pow(C,3)-pow(D,3)+pow(E,3))/(6.0*Dxi*Deta*Dzeta);
		else
			Psi = (pow(A,2) - pow(C,2))/(2.0*Dxi*Deta);
	}
	else{
		Psi = A/Dxi;
	}
}
else{
	if(Gm != 0.0) Psi = 0.0;
	else Psi = 0.5;
}
if(G > 0.0) Psi = 1.0 - Psi;
return Psi;
};

double toolbox::volume_delft_cyl(const double r,const double dr,const double G,const std::vector<double> &gradG){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double psi;
const double COMPUTE_PSI_EPS = 1.0e-8;
double A, B, C, D;
if(fabs(gradG[0]<COMPUTE_PSI_EPS)){
	A = G + 0.5*gradG[1];
	B = G - 0.5*gradG[1];
	psi = 0.5+0.5*(fabs(A)-fabs(B)*(1.0/gradG[1]-dr*G/(pow(r*gradG[1],2))))+\
			0.25*dr/(r*pow(gradG[1],2))*(A*fabs(A)-B*fabs(B));
}
else if(fabs(gradG[1])<COMPUTE_PSI_EPS){
	A = G + 0.5*gradG[0];
	B = G - 0.5*gradG[0];
	psi = 0.5*(1.0+(fabs(A)-fabs(B))/gradG[0]);
}
else{
	A = G + 0.5*( gradG[0]+gradG[1]);
	B = G + 0.5*( gradG[0]-gradG[1]);
	C = G + 0.5*(-gradG[0]+gradG[1]);
	D = G + 0.5*(-gradG[0]-gradG[1]);
	psi = 0.5*dr/(gradG[0]*pow(gradG[1],2)*r)*(r13*(pow(A,2)*fabs(A)-pow(B,2)\
			*fabs(B)-pow(C,2)*fabs(C)+pow(D,2)*fabs(D)))-\
		   0.5*(G+0.5*gradG[0])*(A*fabs(A)-B*fabs(B))+\
			0.5*(G-0.5*gradG[0])*(C*fabs(C)-D*fabs(D))+0.5+\
			0.25/(gradG[0]*gradG[1])*(A*fabs(A)-B*fabs(B)-C*fabs(C)-D*fabs(D));
}
psi = std::max(psi, 0.0);
psi = std::max(psi, 1.0);
return psi;
};
double toolbox::cell_volume_fraction(const Gnode_t &Gnh,const int ic,\
												block_t* bb , const int indic){
/*std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double V;
std::vector<double> Gx(3);
if(cell_contains_front(Gnh, ic, bb)){
	auto i = Gnh.ijk[0]-1;
	auto j = Gnh.ijk[1]-1;
	auto k = Gnh.ijk[2]-1;
	auto ii = i+1-pgv->b->imino_;
	auto jj = j+1-pgv->b->jmino_;
	auto kk = k+1-pgv->b->kmino_;
	Gx[0] = (bb->Gnodes[bb->ijk2ic[ii+1][jj][kk]-1].G-\
				bb->Gnodes[bb->ijk2ic[ii-1][jj][kk]-1].G)*0.5;
	Gx[1] = (bb->Gnodes[bb->ijk2ic[ii][jj+1][kk]-1].G-\
				bb->Gnodes[bb->ijk2ic[ii][jj-1][kk]-1].G)*0.5;
	Gx[2] = (bb->Gnodes[bb->ijk2ic[ii][jj][kk+1]-1].G-\
				bb->Gnodes[bb->ijk2ic[ii][jj][kk-1]-1].G)*0.5;
	if(indic != 0)	V = volume_delft(Gnh.G, Gx);
	else{
		if(pgv->cylindrical){
			V = volume_delft_cyl(pgv->yc[j+pgv->nghost], pgv->dxyz[1],Gnh.G, Gx);
		}
		else{
			V = volume_delft(Gnh.G, Gx);
		}
	}
}
else{
	if(Gnh.G > 0.0) V = 1.0;
	else V = 0.0;
}
return V;*/

};
bool toolbox::is_indicator(const Gnode_t &Gnh, const int ic, block_t* bb){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double V;
if(pgv->cylindrical) V = cell_volume_fraction(Gnh, ic, bb, 1);
else V = cell_volume_fraction(Gnh, ic, bb);
if(V>0.0 && V < 1.0) return true;
else return false;
};

void toolbox::findClosestFrontPointDirect(size_t xBase_S1, size_t xBase_S2, double **xBase){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

std::array<double, 3> Gx;
double scale, R, R0;
for(auto i=0; i<xBase_S1; i++)
	for(auto j=0; j<xBase_S2; j++) 
		xBase[i][j]=std::numeric_limits<double>::max();
if(pgv->cylindrical){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinT; ic++){
			auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
			auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
			auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
			Gx[0]=(pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[pgv->i2c[i-1][j][k]-1].G)*\
					 0.5*pgv->rdxyz[0];
			Gx[1]=(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[pgv->i2c[i][j-1][k]-1].G)*\
					 0.5*pgv->rdxyz[1];
			Gx[2]=(pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[pgv->i2c[i][j][k-1]-1].G)*\
					 0.5*pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)]*pgv->rdxyz[2];
			auto temp = std::accumulate(Gx.begin(), Gx.end(), 0.0);
			temp = pow(temp,2);
			temp = std::max(1.0e-15, temp);
			xBase[0][pgv->ioff+ic]=pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)]-pgv->Gn[ic].G*Gx[0]/temp;
			xBase[1][pgv->ioff+ic]=pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]-pgv->Gn[ic].G*Gx[1]/temp;
			xBase[2][pgv->ioff+ic]=pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)]-pgv->Gn[ic].G*Gx[2]/temp;
		}
	}
}
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinT; ic++){
			auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
			auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
			auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
 			Gx[0]=(pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[pgv->i2c[i-1][j][k]-1].G)*\
					 0.5*pgv->rdxyz[0];
			Gx[1]=(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[pgv->i2c[i][j-1][k]-1].G)*\
					 0.5*pgv->rdxyz[1];
			Gx[2]=(pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[pgv->i2c[i][j][k-1]-1].G)*\
					 0.5*pgv->rdxyz[2];
			auto temp = std::accumulate(Gx.begin(), Gx.end(), 0.0);
			temp = pow(temp,2);
			temp = std::max(1.0e-15, temp);
			xBase[0][pgv->ioff+ic]=pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)]-pgv->Gn[ic].G*Gx[0]/temp;
			xBase[1][pgv->ioff+ic]=pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]-pgv->Gn[ic].G*Gx[1]/temp;
			xBase[2][pgv->ioff+ic]=pgv->zc[i+pgv->b->kmino_-(1-pgv->nghost)]-pgv->Gn[ic].G*Gx[2]/temp;
		}
	}
}
};
void toolbox::curvature(size_t kappa_S, double *kappa){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double **S, **Gx;
std::array<double, 3> Gxx;
double Gxy, Gyz, Gxz, r12;
double xloc, yloc, zloc, x0, y0, z0, r, theta;
std::array<std::array<double, 10>, 27> LS_A;
std::vector<std::vector<double>> A;
A.resize(10);
for(auto i=0; i<10; i++) A[i].resize(10);
std::vector<double> sol(10), bb(10);
std::vector<double> bb1(27);
r12 = 1.0/12.0;
size_t S_S1, S_S2, Gx_S1, Gx_S2;
S = plbuf->getR2Buffer(S_S1, S_S2);
Gx = plbuf->getR2Buffer(Gx_S1, Gx_S2);
if(pgv->cylindrical){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinN; ic++){
			auto i = pgv->Gn->ijk[0]-pgv->b->imino_;
			auto j = pgv->Gn->ijk[1]-pgv->b->jmino_;
			auto k = pgv->Gn->ijk[2]-pgv->b->kmino_;
			Gx[0][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[pgv->i2c[i-1][j][k]-1].G)*\
										  0.5*pgv->rdxyz[0];
			Gx[1][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[pgv->i2c[i][j-1][k]-1].G)*\
										  0.5*pgv->rdxyz[1];
			Gx[2][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[pgv->i2c[i][j][k-1]-1].G)*\
										  0.5*pgv->rdxyz[2];
			Gxx[0] = (pgv->Gn[pgv->i2c[i+1][j][k]-1].G-2.0*pgv->Gn[ic].G+\
						 pgv->Gn[pgv->i2c[i-1][j][k]-1].G)*pgv->rdxyz2[0];
			Gxx[1] = (pgv->Gn[pgv->i2c[i][j+1][k]-1].G-2.0*pgv->Gn[ic].G+\
						 pgv->Gn[pgv->i2c[i][j-1][k]-1].G)*pgv->rdxyz2[1];
			Gxx[2] = (pgv->Gn[pgv->i2c[i][j][k+1]-1].G-2.0*pgv->Gn[ic].G+\
						 pgv->Gn[pgv->i2c[i][j][k-1]-1].G)*pgv->rdxyz2[2];
			S[0][pgv->ioff+ic] = Gxx[0]*(pow(Gx[1][pgv->ioff+ic],2)+pow(Gx[2][pgv->ioff+ic]*\
										pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)],2))+\
										Gxx[1]*(pow(Gx[0][pgv->ioff+ic],2)+pow(Gx[2][pgv->ioff+ic]*\
										pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)],2))+\
										Gxx[2]*(pow(Gx[0][pgv->ioff+ic],2)+pow(Gx[1][pgv->ioff+ic],2)*\
										pow(pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)],2))+\
										Gx[1][pgv->ioff+ic]*(pow(Gx[0][pgv->ioff+ic],2)+\
										pow(Gx[1][pgv->ioff+ic],2)+2.0*pow(Gx[2][pgv->ioff+ic]*\
										pgv->ryc[j+pgv->b->jmino_-\
										(1-pgv->nghost)],2))*pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
			auto temp = sqrt(pow(Gx[0][pgv->ioff+ic],2)+pow(Gx[1][pgv->ioff+ic],2)+\
										pow(Gx[2][pgv->ioff+ic]*pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)],2));
			S[1][pgv->ioff+ic] = std::max(temp, 1.0e-10);
		}
	}
}
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinT; ic++){
			auto i = pgv->Gn[ic].ijk[0] - pgv->b->imino_;
			auto j = pgv->Gn[ic].ijk[1] - pgv->b->jmino_;
			auto k = pgv->Gn[ic].ijk[2] - pgv->b->kmino_;
			Gx[0][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[pgv->i2c[i-1][j][k]-1].G)*\
										0.5*pgv->rdxyz[0];
			Gx[1][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[pgv->i2c[i][j-1][k]-1].G)*\
										0.5*pgv->rdxyz[1];
			Gx[2][pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[pgv->i2c[i][j][k-1]-1].G)*\
										0.5*pgv->rdxyz[2];
 			Gxx[0] = (pgv->Gn[pgv->i2c[i+1][j][k]-1].G-2.0*pgv->Gn[ic].G\
						+pgv->Gn[pgv->i2c[i-1][j][k]-1].G)*pgv->rdxyz2[0];
			Gxx[1] = (pgv->Gn[pgv->i2c[i][j+1][k]-1].G-2.0*pgv->Gn[ic].G+\
						 pgv->Gn[pgv->i2c[i][j-1][k]-1].G)*pgv->rdxyz2[1];
			Gxx[2] = (pgv->Gn[pgv->i2c[i][j][k+1]-1].G-2.0*pgv->Gn[ic].G+\
						 pgv->Gn[pgv->i2c[i][j][k-1]-1].G)*pgv->rdxyz2[2];
			S[0][pgv->ioff+ic] = Gxx[0]*(pow(Gx[1][pgv->ioff+ic],2)+pow(Gx[2][pgv->ioff+ic],2))+\
										Gxx[1]*(pow(Gx[0][pgv->ioff+ic],2)+pow(Gx[2][pgv->ioff+ic],2))+\
										Gxx[2]*(pow(Gx[0][pgv->ioff+ic],2)+pow(Gx[1][pgv->ioff+ic],2));
			auto temp=pow(Gx[0][pgv->ioff+ic],2)+pow(Gx[1][pgv->ioff+ic],2)+pow(Gx[2][pgv->ioff+ic],2);
			S[1][pgv->ioff+ic]=std::max(sqrt(temp),1.0e-10);
		}
	}
}

pbou->updateGhostR2(Gx, Gx_S1, Gx_S2, pbou->S_GRADIENT);
for(auto i=0; i<kappa_S; i++) kappa[i] = 0.0; // checked
if(pgv->cylindrical){
	if(pgv->ijkm_gl[2]>1){
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinT; ic++){
				auto i = pgv->Gn[ic].ijk[0]-pgv->b->imino_;
				auto j = pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
				auto k = pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
				if(j+pgv->b->jmino_==1){
					auto count = -1;
					xloc = pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)];
					yloc = pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]*\
							 cos(pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)]);
					zloc = pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]*\
							 cos(pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)]);
					for(auto ii=i-1; ii<i+1; ii++){
						for(auto jj=j-1; jj<j+1; jj++){
							for(auto kk=k-1; kk<k+1; kk++){
								count += 1;
								LS_A[0][count]=1.0;
								LS_A[1][count]=pgv->xc[ii+pgv->b->imino_-(1-pgv->nghost)]-xloc;
								LS_A[2][count]=pgv->yc[jj+pgv->b->jmino_-(1-pgv->nghost)]*\
													cos(pgv->zc[kk+pgv->b->kmino_-(1-pgv->nghost)])-yloc;
								LS_A[3][count]=pgv->yc[jj+pgv->b->jmino_-(1-pgv->nghost)]*\
													sin(pgv->zc[kk+pgv->b->kmino_-(1-pgv->nghost)])-zloc;
								LS_A[4][count]=0.5*pow(pgv->xc[ii+pgv->b->imino_-(1-pgv->nghost)]-xloc,2);
								LS_A[5][count]=0.5*pow(pgv->yc[jj+pgv->b->jmino_-(1-pgv->nghost)]*\
														cos(pgv->zc[kk+pgv->b->kmino_-(1-pgv->nghost)])-yloc,2);
								LS_A[6][count]=0.5*pow(pgv->yc[jj+pgv->b->jmino_-(1-pgv->nghost)]*\
														sin(pgv->zc[kk+pgv->b->kmino_-(1-pgv->nghost)])-zloc,2);
								LS_A[7][count]=(pgv->xc[ii+pgv->b->imino_-(1-pgv->nghost)]-xloc)*\
													(pgv->yc[jj+pgv->b->jmino_-(1-pgv->nghost)]*\
													 cos(pgv->zc[kk+pgv->b->kmino_-(1-pgv->nghost)]-yloc));
								LS_A[8][count]=(pgv->xc[ii+pgv->b->imino_-(1-pgv->nghost)]-xloc)*\
													(pgv->yc[jj+pgv->b->jmino_-(1-pgv->nghost)]\
													 *sin(pgv->zc[kk+pgv->b->kmino_-(1-pgv->nghost)])-zloc);
								LS_A[9][count]=(pgv->yc[jj+pgv->b->jmino_-(1-pgv->nghost)]*\
													cos(pgv->zc[kk+pgv->b->kmino_-(1-pgv->nghost)])-yloc)*\
												(pgv->yc[jj+pgv->b->jmino_-(1-pgv->nghost)]*\
												sin(pgv->zc[kk+pgv->b->kmino_-(1-pgv->nghost)])-zloc);
							}
						}
					}
					for(auto ii=0; ii<10; ii++){
						for(auto jj=0; jj<10; jj++){
							auto temp=0.0;
							for(auto q=0; q<LS_A[0].size(); q++){
								temp+=LS_A[ii][q]*LS_A[jj][q];
							}
							A[ii][jj] = temp;
						}
					}
					count = -1;
					for(auto ii=i-1; ii<i+1; ii++){
						for(auto jj=j-1; jj<j+1; jj++){
							for(auto kk=k-1; kk<k+1; kk++){
								bb1[count]=heaviside_st(pgv->Gn[pgv->i2c[ii][jj][kk]-1].G);
							}
						}
					}
					for(auto ii=0; ii<10; i++){
						auto temp = 0.0;
						for(auto q=0; bb1.size(); q++){
							temp += LS_A[ii][q]*bb1[q];
						}
						bb[ii] = temp;
					}
					auto n=10;
					solve_linear_system(A,bb,sol,n);
					kappa[pgv->ioff+ic]=\
					-(pow(sol[1],2)*sol[5]+pow(sol[1],2)*sol[6]+\
					  pow(sol[2],2)*sol[6]+pow(sol[2],2)*sol[4]+\
					  pow(sol[3],2)*sol[4]+pow(sol[3],2)*sol[5]-\
					  2.0*sol[1]*sol[2]*sol[7]-\
					  2.0*sol[1]*sol[3]*sol[8]-\
					  2.0*sol[2]*sol[3]*sol[9])/(pow(pow(sol[1],2)+pow(sol[2],2)+pow(sol[3],2),1.5)+\
														  std::numeric_limits<double>::epsilon());
				}
				else{
					Gxy=(Gx[0][pgv->ioff+pgv->i2c[i][j+1][k]]-Gx[0][pgv->ioff+pgv->i2c[i][j-1][k]])*\
						  0.5*pgv->rdxyz[1];
					Gxz=(Gx[0][pgv->ioff+pgv->i2c[i][j][k+1]]-Gx[0][pgv->ioff+pgv->i2c[i][j][k-1]])*\
						  0.5*pgv->rdxyz[2];
					Gyz=(Gx[1][pgv->ioff+pgv->i2c[i][j][k+1]]-Gx[1][pgv->ioff+pgv->i2c[i][j][k-1]])*\
						  0.5*pgv->rdxyz[2];
					auto value=Gxy*Gx[0][pgv->ioff+ic]*Gx[1][pgv->ioff+ic]+\
								 (Gxz*Gx[0][pgv->ioff+ic]*Gx[2][pgv->ioff+ic]+\
								  Gyz*Gx[1][pgv->ioff+ic]*Gx[2][pgv->ioff+ic])*\
								  pow(pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)],2);
					kappa[pgv->ioff+ic]=-(S[0][pgv->ioff+ic]-2.0*value)/ \
												 pow(S[1][pgv->ioff+ic],3);
				}
				if(kappa[pgv->ioff+ic]<0.0){
					kappa[pgv->ioff+ic]=std::max(kappa[pgv->ioff+ic],pgv->kappa_min);
				}
				else{
					kappa[pgv->ioff+ic]=std::min(kappa[pgv->ioff+ic],pgv->kappa_max);
				}
			}
		}
	}
	else{
		for(auto ibl=0; ibl<pgv->nbl; ibl++){
			pgv->setBlockPointers(ibl);
			for(auto ic=0; ic<pgv->b->NinT; ic++){
				auto i=pgv->Gn[ic].ijk[0]-pgv->b->imino_;
				auto j=pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
				auto k=pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
				Gxy=(Gx[0][pgv->ioff+pgv->i2c[i][j+1][k]]-\
					  Gx[0][pgv->ioff+pgv->i2c[i][j-1][k]])*0.5*pgv->rdxyz[1];
				Gxz=(Gx[0][pgv->ioff+pgv->i2c[i][j][k+1]]-\
					  Gx[0][pgv->ioff+pgv->i2c[i][j][k-1]])*0.5*pgv->rdxyz[2];
				Gyz=(Gx[1][pgv->ioff+pgv->i2c[i][j][k+1]]-\
					  Gx[1][pgv->ioff+pgv->i2c[i][j][k-1]])*0.5*pgv->rdxyz[2];
				auto value=Gxy*Gx[0][pgv->ioff+ic]*Gx[1][pgv->ioff+ic]+\
							 (Gxz*Gx[0][pgv->ioff+ic]*Gx[2][pgv->ioff+ic]+\
							  Gyz*Gx[1][pgv->ioff+ic]*Gx[2][pgv->ioff+ic])*\
							  pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
				kappa[pgv->ioff+ic]=-(S[0][pgv->ioff+ic]-2.0*value)/ \
											 pow(S[1][pgv->ioff+ic],3);
				if(kappa[pgv->ioff+ic]<0.0){
					kappa[pgv->ioff+ic]=std::max(kappa[pgv->ioff+ic],pgv->kappa_min);
				}
				else{
					kappa[pgv->ioff+ic]=std::min(kappa[pgv->ioff+ic],pgv->kappa_max);
				}
			}
		}
	}
}
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinT; ic++){
			auto i=pgv->Gn[ic].ijk[0]-pgv->b->imino_;
			auto j=pgv->Gn[ic].ijk[1]-pgv->b->jmino_;
			auto k=pgv->Gn[ic].ijk[2]-pgv->b->kmino_;
			Gxy=(Gx[0][pgv->ioff+pgv->i2c[i][j+1][k]]-\
				  Gx[0][pgv->ioff+pgv->i2c[i][j-1][k]])*0.5*pgv->rdxyz[1];
			Gxz=(Gx[0][pgv->ioff+pgv->i2c[i][j][k+1]]-\
				  Gx[0][pgv->ioff+pgv->i2c[i][j][k-1]])*0.5*pgv->rdxyz[2];
			Gyz=(Gx[1][pgv->ioff+pgv->i2c[i][j][k+1]]-\
				  Gx[1][pgv->ioff+pgv->i2c[i][j][k-1]])*0.5*pgv->rdxyz[2];
			auto value=Gxy*Gx[0][pgv->ioff+ic]*Gx[1][pgv->ioff+ic]+\
						 (Gxz*Gx[0][pgv->ioff+ic]*Gx[2][pgv->ioff+ic]+\
						  Gyz*Gx[1][pgv->ioff+ic]*Gx[2][pgv->ioff+ic]);
			kappa[pgv->ioff+ic]=-(S[0][pgv->ioff+ic]-2.0*value)/ \
										 pow(S[1][pgv->ioff+ic],3);
			if(kappa[pgv->ioff+ic]<0.0){
				kappa[pgv->ioff+ic]=std::max(kappa[pgv->ioff+ic],pgv->kappa_min);
			}
			else{
				kappa[pgv->ioff+ic]=std::min(kappa[pgv->ioff+ic],pgv->kappa_max);
			}
		}
	}
}
plbuf->freeR2Buffer(S_S1, S_S2, S);
plbuf->freeR2Buffer(Gx_S1, Gx_S2, Gx);
//pbou->updateGhostR1(kappa);
};
void toolbox::solve_linear_system(std::vector<std::vector<double>> &A,\
					std::vector<double> &B,std::vector<double> &X, const int n){
for(auto i=0; i<n; i++){
	B[i] /= A[i][i];
	for(auto j=0; j<n; j++) A[i][j] /= A[i][i];
	for(auto l=i+1; l<n; l++){
		B[l] = B[l]-A[l][i]*B[i];
		for(auto j=0; j<n; j++) A[l][j] = A[l][j] - A[l][i]*A[i][j];
	}
}
for(auto i=n-1; i>-1; i--){
	X[i] = B[i];
	for(auto l=i+1; l<n; l++) X[i] = A[i][l]*X[l];
}
};
double toolbox::heaviside_st(const double G){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double H;
const double heaviside_width=0.12;
const double Pi=3.14159;
if(G<=-heaviside_width) H = 0.0;
else if(G >= heaviside_width) H = 1.0;
else H = 0.5 + 0.5/heaviside_width * G+0.5/Pi * sin(Pi/heaviside_width*G);
H = H - 0.5;
return H;
}
void toolbox::curvature_indicator(size_t kappa_S,double*kappa,size_t indicator_S,double*indicator){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double v;
curvature(kappa_S, kappa);
for(auto i=0; i<indicator_S; i++) indicator[i] = 0.0;
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	if(pgv->cylindrical){
		for(auto ic=0; ic<pgv->b->NinT; ic++){
			if(is_indicator(pgv->Gn[ic], ic, pgv->b)) indicator[pgv->ioff+ic] = 1.0;
		}
	}
	else{
		for(auto ic=0; ic<pgv->b->NinT; ic++){
			if(is_indicator(pgv->Gn[ic], ic, pgv->b)) indicator[pgv->ioff+ic] = 1.0;
		}
	}
}
pbou->updateGhostR1(indicator);
};
void toolbox::direct_curvature_indicator(size_t kappa_S,double*kappa,size_t indicator_S,double *indicator){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double *kappa_cv;
double **xBase;
std::vector<double> curv(8), xk(3), v;
size_t xBase_S1, xBase_S2;
xBase = plbuf->getR2Buffer(xBase_S1, xBase_S2);
findClosestFrontPointDirect(xBase_S1, xBase_S2, xBase);
size_t kappa_cv_S;
kappa_cv = plbuf->getR1Buffer(kappa_cv_S);
for(auto i=0; i<indicator_S; i++) indicator[i] = 0.0;
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinN; ic++){
		if(is_indicator(pgv->Gn[ic], ic, pgv->b)){
			auto ip=std::min(pgv->b->imax_,\
					static_cast<int>(floor((xBase[0][pgv->ioff+ic]-pgv->xc[pgv->nghost])*pgv->rdxyz[0])+1));
			auto jp=std::min(pgv->b->jmax_,\
					static_cast<int>(floor((xBase[1][pgv->ioff+ic]-pgv->yc[pgv->nghost])*pgv->rdxyz[1])+1));
			auto kp=std::min((pgv->b->kmax_),\
					static_cast<int>(floor((xBase[2][pgv->ioff+ic]-pgv->zc[pgv->nghost])*pgv->rdxyz[2])+1));
			auto i = std::max(pgv->b->imino_,ip)-pgv->b->imino_;
			auto j = std::max(pgv->b->jmino_,jp)-pgv->b->jmino_;
			auto k = std::max(pgv->b->kmino_,kp)-pgv->b->kmino_;
			xk[0] = xBase[0][pgv->ioff+ic]-pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)]*pgv->rdxyz[0];
			xk[1] = xBase[1][pgv->ioff+ic]-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]*pgv->rdxyz[1];
			xk[2] = xBase[2][pgv->ioff+ic]-pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)]*pgv->rdxyz[2];
			curv[0]=kappa_cv[pgv->ioff+pgv->i2c[i][j][k]-1];
			curv[1]=kappa_cv[pgv->ioff+pgv->i2c[i+1][j][k]-1];
			curv[2]=kappa_cv[pgv->ioff+pgv->i2c[i][j+1][k]-1];
			curv[3]=kappa_cv[pgv->ioff+pgv->i2c[i+1][j+1][k]-1];
			curv[4]=kappa_cv[pgv->ioff+pgv->i2c[i][j+1][k]-1];
			curv[5]=kappa_cv[pgv->ioff+pgv->i2c[i][j][k+1]-1];
			curv[6]=kappa_cv[pgv->ioff+pgv->i2c[i+1][j][k+1]-1];
			curv[7]=kappa_cv[pgv->ioff+pgv->i2c[i+1][j+1][k+1]-1];
			kappa[pgv->ioff+ic] = triLinear(curv, xk);
			indicator[pgv->ioff+ic] = 1.0;
		}
		else{
			kappa[pgv->ioff+ic] = 0.0;
			indicator[pgv->ioff+ic] = 0.0;
		}
	}
}
plbuf->freeR2Buffer(xBase_S1, xBase_S2, xBase);
plbuf->freeR1Buffer(kappa_cv);
pbou->updateGhostR1(kappa);
pbou->updateGhostR1(indicator);
};
double toolbox::triLinear(const std::vector<double>& p, const std::vector<double> &x){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

std::array<double, 6> pp;
for(auto i=0; i<4; i++){
	pp[i] = p[i]+(p[i+4]-p[i])*x[2];
}
pp[4] = pp[0]+(pp[2]-pp[0])*x[1];
pp[5] = pp[1]+(pp[3]-pp[1])*x[1];
return pp[4]+(pp[5]-pp[4])*x[0];
};

void toolbox::calcCurvature(size_t kappa_S,double*kappa,size_t indicator_S,double* indicator){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

if(pgv->calcCurvatureType == pgv->ccNodes)
	curvature_indicator(kappa_S,kappa,indicator_S,indicator);
else if(pgv->calcCurvatureType == pgv->ccDirect)
	direct_curvature_indicator(kappa_S,kappa,indicator_S,indicator);
else{
	std::cout<<pgv->clit<<"Error: implement this curvature calc type="<<\
	pgv->calcCurvatureType<<std::endl;
	pp->parallel_kill(87);
}
}

void toolbox::normal(size_t n_S1, size_t n_S2, double **n){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double denom;
for(auto i=0; i<n_S1; i++){
	for(auto j=0; j<n_S2; j++){
		n[i][j] = 0.0;
	}
}
if(pgv->cylindrical){
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinT; ic++){
			auto i = pgv->Gn[ic].ijk[0] - pgv->b->imino_;
			auto j = pgv->Gn[ic].ijk[1] - pgv->b->jmino_;
			auto k = pgv->Gn[ic].ijk[2] - pgv->b->kmino_;
			n[0][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[pgv->i2c[i-1][j][k]-1].G)*\
										 0.5*pgv->rdxyz[0];
			n[1][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[pgv->i2c[i][j-1][k]-1].G)*\
										 0.5*pgv->rdxyz[1];
			n[2][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[pgv->i2c[i][j][k-1]-1].G)*\
										 0.5*pgv->rdxyz[2]*pgv->ryc[j+pgv->b->jmino_-(1-pgv->nghost)];
			auto temp = 0.0;
			for(auto i=0; i<3; i++){
				temp += pow(n[i][pgv->ioff+ic],2);
			}
			denom = std::max(sqrt(temp), 1.0e-10);
			for(auto i=0; i<3; i++){
				n[i][pgv->ioff+ic] /= denom;
			}
		}
	}
}
else{
	for (auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinT; ic++){
			auto i = pgv->Gn[ic].ijk[0] - pgv->b->imino_;
			auto j = pgv->Gn[ic].ijk[1] - pgv->b->jmino_;
			auto k = pgv->Gn[ic].ijk[2] - pgv->b->kmino_;
			n[0][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[pgv->i2c[i-1][j][k]-1].G)*\
										 0.5*pgv->rdxyz[0];
			n[1][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[pgv->i2c[i][j-1][k]-1].G)*\
										 0.5*pgv->rdxyz[1];
			n[2][pgv->ioff+ic] = (pgv->Gn[pgv->i2c[i][j][k+1]-1].G-pgv->Gn[pgv->i2c[i][j][k-1]-1].G)*\
										 0.5*pgv->rdxyz[2];
			auto temp = 0.0;
			for(auto i=0; i<3; i++){
				temp += pow(n[i][pgv->ioff+ic],2);
			}
			denom = std::max(sqrt(temp), 1.0e-10);
			for(auto i=0; i<3; i++){
				n[i][pgv->ioff+ic] /= denom;
			}
		}
	}
}
pbou->updateGhostR2(n, n_S1, n_S2, pbou->S_GRADIENT);
};
double toolbox::volume(){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double V, Vl;
std::array<double, 3> Gx;
Vl = 0.0;
if(pgv->cylindrical){
	auto temp = pgv->dxyz[0]*pgv->dxyz[1]*pgv->dxyz[2];
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinX; ic++){
			auto i = pgv->Gn[ic].ijk[0] - 1;
			auto j = pgv->Gn[ic].ijk[1] - 1;
			auto k = pgv->Gn[ic].ijk[2] - 1;
			Vl += cell_volume_fraction(pgv->Gn[ic], ic+1, pgv->b)*pgv->yc[j+pgv->nghost]*temp;
		}
	}
}
else{
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinX; ic++){
			Vl += Vl + cell_volume_fraction(pgv->Gn[ic], ic+1, pgv->b);
		}
	}
	Vl *= pgv->cell_volume;
}
pp->ierr = MPI_Allreduce(&Vl, &V, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
return V;
}
double toolbox::volumeTband(){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double V;
double Vl = 0.0;
std::array<double, 3> Gx;
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic = 0; ic<pgv->b->NinT; ic++){
		Vl += cell_volume_fraction(pgv->Gn[ic], ic+1, pgv->b);
	}
}
Vl *= pgv->cell_volume;
pp->ierr = MPI_Allreduce(&Vl, &V, 1, MPI_DOUBLE, MPI_SUM, pp->MY_G_WORLD);
return V;
}
void toolbox::interfaceHolder(){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);
block_t *bb = nullptr;
MPI_Request *requests, *requestsInt;
requests = new MPI_Request[pp->nprocs];  // delete the pointer; // checked
requestsInt = new MPI_Request[pp->nprocs]; // delete the pointer; // checked
bool selfComm;
double r, vin, vout;
std::vector<double> xyz(3);
if(pgv->holder){
	vin = volumeTband();
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		if(pgv->holder_name == "plane"||pgv->holder_name=="Plane"||pgv->holder_name=="PLANE"){
			for(auto ic=0; ic<pgv->b->NinT; ic++){
				xyz[0]=pgv->xc[pgv->Gn[ic].ijk[0]-(1-pgv->nghost)];
				xyz[1]=pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)];
				xyz[2]=pgv->zc[pgv->Gn[ic].ijk[2]-(1-pgv->nghost)];
				if(xyz[0] >= pgv->holder_xyzs[0] && xyz[0] <= pgv->holder_xyze[0] && \
					xyz[1] >= pgv->holder_xyzs[1] && xyz[1] <= pgv->holder_xyze[1] && \
					xyz[2] >= pgv->holder_xyzs[2] && xyz[2] <= pgv->holder_xyze[2]){
					pgv->Gn[ic].G = pinit->G_init_value(xyz);
				}
			}
		}
		else if(pgv->holder_name=="nozzle"||pgv->holder_name=="Nozzle"||pgv->holder_name=="NOZZLE"){
			for(auto ic=0; ic<pgv->b->NinT; ibl++){
				xyz[0]=pgv->xc[pgv->Gn[ic].ijk[0]-(1-pgv->nghost)];
				xyz[1]=pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)];
				xyz[2]=pgv->zc[pgv->Gn[ic].ijk[2]-(1-pgv->nghost)];
				if(xyz[0] >= pgv->holder_xyzs[0] && xyz[0] <= pgv->holder_xyze[0]){
					auto temp=pow(xyz[1]-pgv->holder_center[1],2)+pow(xyz[2]-pgv->holder_center[2],2);
					auto r=pgv->holder_radius-sqrt(temp);
					if(fabs(r) < pgv->holder_dr) pgv->Gn[ic].G = r;
				}
			}
		}
		else if(pgv->holder_name=="nozzle_z"||pgv->holder_name=="Nozzle_z"||\
				  pgv->holder_name=="NOZZLE_Z"){
			for(auto ic=0; ic<pgv->b->NinT; ibl++){
				xyz[0]=pgv->xc[pgv->Gn[ic].ijk[0]-(1-pgv->nghost)];
				xyz[1]=pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)];
				xyz[2]=pgv->zc[pgv->Gn[ic].ijk[2]-(1-pgv->nghost)];
				if(xyz[2] >= pgv->holder_xyzs[2] && xyz[2] <= pgv->holder_xyze[2]){
					auto temp=pow(xyz[0]-pgv->holder_center[0],2)+pow(xyz[1]-pgv->holder_center[1],2);
					auto r=pgv->holder_radius-sqrt(temp);
					if(fabs(r) < pgv->holder_dr) pgv->Gn[ic].G = r;
				}
			}
		}
		else if(pgv->holder_name=="nozzle_y"||pgv->holder_name=="Nozzle-y"||\
				  pgv->holder_name=="NOZZLE_Y"){
			for(auto ic=0; ic<pgv->b->NinT; ibl++){
				xyz[0]=pgv->xc[pgv->Gn[ic].ijk[0]-(1-pgv->nghost)];
				xyz[1]=pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)];
				xyz[2]=pgv->zc[pgv->Gn[ic].ijk[2]-(1-pgv->nghost)];
				if(xyz[1] >= pgv->holder_xyzs[1] && xyz[1] <= pgv->holder_xyze[1]){
					auto temp=pow(xyz[0]-pgv->holder_center[0],2)+pow(xyz[2]-pgv->holder_center[2],2);
					auto r = pgv->holder_radius - sqrt(temp);
					if(fabs(r)<pgv->holder_dr)
					pgv->Gn[ic].G = pinit->G_init_value(xyz);	
				}
			}
		}
		else if(pgv->holder_name=="sphere_cap"||pgv->holder_name=="Sphere_Cap"||\
				  pgv->holder_name=="SPHERE_CAP"){
			if(plt->time < pgv->holder_dt){
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					xyz[0]=pgv->xc[pgv->Gn[ic].ijk[0]-(1-pgv->nghost)];
					xyz[1]=pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)];
					xyz[2]=pgv->zc[pgv->Gn[ic].ijk[2]-(1-pgv->nghost)];
					if(xyz[2] > pgv->holder_xyzs[2] && xyz[2] < pgv->holder_xyze[2]){
						pgv->Gn[ic].G = pinit->G_init_value(xyz);
					}
				}
			}
		}
		else if(pgv->holder_name=="sphere_cap_cyl"||pgv->holder_name=="Sphere_Cap_Cyl"||\
				  pgv->holder_name=="SPHERE_CAP_CYL"){
			if(plt->time < pgv->holder_dt){
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					xyz[0]=pgv->xc[pgv->Gn[ic].ijk[0]-(1-pgv->nghost)];
					xyz[1]=pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)];
					xyz[2]=pgv->zc[pgv->Gn[ic].ijk[2]-(1-pgv->nghost)];
					if(xyz[0] > pgv->holder_xyzs[0] && xyz[0] < pgv->holder_xyze[0]){
						pgv->Gn[ic].G = pinit->G_init_value(xyz);
					}
				}
			}
		}
		else if(pgv->holder_name=="circle_cap"||pgv->holder_name=="Circle_Cap"||\
				  pgv->holder_name=="CIRCLE_CAP"){
			if(plt->time < pgv->holder_dt){
				for(auto ic=0; ic<pgv->b->NinN; ic++){
					xyz[0]=pgv->xc[pgv->Gn[ic].ijk[0]-(1-pgv->nghost)];
					xyz[1]=pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)];
					xyz[2]=pgv->zc[pgv->Gn[ic].ijk[2]-(1-pgv->nghost)];
					if(xyz[1] > pgv->holder_xyzs[1] && xyz[1] < pgv->holder_xyze[1]){
						pgv->Gn[ic].G = pinit->G_init_value(xyz);
					}
				}
			}
		}
	}
//	vout = volumeTband();
	if(pp->myrank==0 && pgv->verbose)
		std::cout<<pgv->clit<<"volume change in interface holder: "<<vout-vin<<std::endl;
}
delete [] requests;
delete [] requestsInt;
}
double toolbox::calcAmplitude_quarter(){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double A, Ah;
Ah = std::numeric_limits<double>::max();
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	if(pgv->b->ijk0[0]+1<pgv->ijkm_gl[0]/4 && pgv->b->ijk0[0]+pgv->b->ijkm[0]>=pgv->ijkm_gl[0]/4){
		for(auto ic=0; ic<pgv->b->NinT; ic++){
			auto i = pgv->Gn[ic].ijk[0] - pgv->b->imino_;
			auto j = pgv->Gn[ic].ijk[1] - pgv->b->jmino_;
			auto k = pgv->Gn[ic].ijk[2] - pgv->b->kmino_;
			if(i+pgv->b->imino_ == pgv->ijkm_gl[0]/4 && k+pgv->b->kmino_ == 1){
				if(pgv->Gn[ic].G*pgv->Gn[pgv->i2c[i][j+1][k]-1].G <= 0.0){
					if(fabs(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[ic].G<1.0e-10)){
						Ah = 0.5*(pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]+\
									 pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]);
					}
					else{
						Ah = pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]-pgv->Gn[ic].G*\
							  pgv->dxyz[1]/(pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[ic].G);
					}
				}
			}
		}
	}
}
pp->ierr = MPI_Reduce(&Ah, &A, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
return A;
}

void toolbox::diag_RT3D(vector<double> &val, vector<double> &val_spike, vector<double> &val_min,\
								vector<double> &val_bubble, vector<double> & val_max){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);


std::vector<double> velXminh(3), velXmaxh(3);
double xMinh, xMaxh, kappaXMinh, kappaXMaxh, sumIndicatorMinh, sumIndicatorMaxh;
double sumIndicatorMin, sumIndicatorMax, areah, area, xh, yh, zh;
double *kappa_cv;
vector<double> valh_min(7), valh_max(7), valh(7,0);
vector<vector<double>> valh_spike(7, vector<double>(4,0));
vector<vector<double>> valh_bubble(7, vector<double>(4,0));
vector<vector<double>> bufh;
double ***bufhh;
std::fill(val.begin(), val.end(), 0.0);
std::fill(val_min.begin(), val_min.end(), 0.0);
std::fill(val_max.begin(), val_max.end(), 0.0);
std::fill(val_spike.begin(), val_spike.end(), 0.0);
std::fill(val_bubble.begin(), val_bubble.end(), 0.0);
valh_min[1] = std::numeric_limits<double>::max();
valh_max[1] = -std::numeric_limits<double>::max();
for(auto i=0; i<4; i++){
	valh_spike[1][i] =  std::numeric_limits<double>::max();
	valh_bubble[1][i]= -std::numeric_limits<double>::max();
}
size_t kappa_cv_S;
kappa_cv = plbuf->getR1Buffer(kappa_cv_S);
curvature(kappa_cv_S,kappa_cv);
for(auto ibl=0; ibl< pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinT; ic++){
		auto i = pgv->Gn[ic].ijk[0] - pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1] - pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2] - pgv->b->kmino_;
		valh[1] += delta_L(pgv->Gn[ic].G);
		auto ic2 = pgv->i2c[i][j+1][k];
		if(pgv->Gn[ic].G*pgv->Gn[ic2-1].G < 0.0){
			if(fabs(pgv->Gn[ic2-1].G-pgv->Gn[ic].G) < 1.0e-10){
				yh = 0.5*(pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]+\
							 pgv->yc[j+1+pgv->b->jmino_-(1-pgv->nghost)]);
			}
			else{
				yh = pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]-pgv->Gn[ic].G*\
					  pgv->dxyz[1]/(pgv->Gn[ic2-1].G-pgv->Gn[ic].G);
			}
			if(yh < valh_min[1]){
				valh_min[0] = pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)];
				valh_min[1] = yh;
				valh_min[2] = pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)];
				for(auto q=0; q<3; q++){
					valh[q+3]=pgv->Gn[ic].V[q]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
								 pgv->rdxyz[1]*(pgv->Gn[ic2-1].V[q]-pgv->Gn[ic].V[q]);

				}
				valh_min[6]=kappa_cv[ic]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*pgv->rdxyz[1]*\
								(kappa_cv[ic2-1]-kappa_cv[ic]);

			}
			if(yh>valh_max[1]){
				valh_max[0] = pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)];
				valh_max[1] = yh;
				valh_max[2] = pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)];
				for(auto q=0; q<3; q++){
					valh_max[q+3]=pgv->Gn[ic].V[q]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
									  pgv->rdxyz[1]*(pgv->Gn[ic2-1].V[q]-pgv->Gn[ic].V[q]);
				}
				valh_max[6]=kappa_cv[ic]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*pgv->rdxyz[1]*\
							  (kappa_cv[ic2-1]-kappa_cv[ic]);
			}
			if(i+pgv->b->imino_==pgv->ijkm_gl[0]/2 && k+pgv->b->kmino_==pgv->ijkm_gl[2]/2){
				if(yh<valh_spike[1][0]){
					valh_spike[0][0] = pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)];
					valh_spike[1][0] = yh;
					valh_spike[2][0] = pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)];
					for(auto q=0; q<3; q++){
					 valh_spike[q+3][0]=pgv->Gn[ic].V[q]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
							pgv->rdxyz[1]*(pgv->Gn[ic2-1].V[q]-pgv->Gn[ic].V[q]);
					}
					valh_spike[6][0]=kappa_cv[ic]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
										  pgv->rdxyz[1]*(kappa_cv[ic2-1]-kappa_cv[ic]);
				}
			}
			if(i+pgv->b->imino_==pgv->ijkm_gl[0]/2 && k+pgv->b->kmino_==pgv->ijkm_gl[2]/2+1){
				if(yh<valh_spike[1][1]){
					valh_spike[0][1] = pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)];
					valh_spike[1][1] = yh;
					valh_spike[2][1] = pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)];
					for(auto q=0; q<3; q++){
					 valh_spike[q+3][1]=pgv->Gn[ic].V[q]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
							pgv->rdxyz[1]*(pgv->Gn[ic2-1].V[q]-pgv->Gn[ic].V[q]);
					}
					valh_spike[6][1]=kappa_cv[ic]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
										  pgv->rdxyz[1]*(kappa_cv[ic2-1]-kappa_cv[ic]);
				}
			}
			if(i+pgv->b->imino_==pgv->ijkm_gl[0]/2+1 && k+pgv->b->kmino_==pgv->ijkm_gl[2]/2){
				if(yh<valh_spike[1][2]){
					valh_spike[0][2] = pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)];
					valh_spike[1][2] = yh;
					valh_spike[2][2] = pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)];
					for(auto q=0; q<3; q++){
					 valh_spike[q+3][2]=pgv->Gn[ic].V[q]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
							pgv->rdxyz[1]*(pgv->Gn[ic2-1].V[q]-pgv->Gn[ic].V[q]);
					}
					valh_spike[6][1]=kappa_cv[ic]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
										  pgv->rdxyz[1]*(kappa_cv[ic2-1]-kappa_cv[ic]);
				}
			}
			if(i+pgv->b->imino_==pgv->ijkm_gl[0]/2+1 && k+pgv->b->kmino_==pgv->ijkm_gl[2]/2){
				if(yh<valh_spike[1][3]){
					valh_spike[0][3] = pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)];
					valh_spike[1][3] = yh;
					valh_spike[2][3] = pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)];
					for(auto q=0; q<3; q++){
					 valh_spike[q+3][3]=pgv->Gn[ic].V[q]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
							pgv->rdxyz[1]*(pgv->Gn[ic2].V[q]-pgv->Gn[ic].V[q]);
					}
					valh_spike[6][3]=kappa_cv[ic]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
										  pgv->rdxyz[1]*(kappa_cv[ic2-1]-kappa_cv[ic]);
				}
			}
			if(i+pgv->b->imino_== 1 && k+pgv->b->kmino_== 1){
				if(yh>valh_bubble[1][0]){
					valh_bubble[0][0] = pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)];
					valh_bubble[1][0] = yh;
					valh_bubble[2][0] = pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)];
					for(auto q=0; q<3; q++){
					 valh_bubble[q+3][0]=pgv->Gn[ic].V[q]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])\
					 *pgv->rdxyz[1]*(pgv->Gn[ic2-1].V[q]-pgv->Gn[ic].V[q]);
					}
					valh_bubble[6][0]=kappa_cv[ic]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
											pgv->rdxyz[1]*(kappa_cv[ic2-1]-kappa_cv[ic]);
				}
			}
			if(i+pgv->b->imino_ == 1 && k+pgv->b->kmino_ == pgv->ijkm_gl[2]){
				if(yh>valh_bubble[1][1]){
					valh_bubble[0][1] = pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)];
					valh_bubble[1][1] = yh;
					valh_bubble[2][1] = pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)];
					for(auto q=0; q<3; q++){
					valh_bubble[q+3][1]=pgv->Gn[ic].V[q]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
					pgv->rdxyz[1]*(pgv->Gn[ic2-1].V[q]-pgv->Gn[ic].V[q]);
					}
					valh_bubble[6][1]=kappa_cv[ic]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
											pgv->rdxyz[1]*(kappa_cv[ic2-1]-kappa_cv[ic]);				
				}

			}
			if(i+pgv->b->imino_==pgv->ijkm_gl[0] && k+pgv->b->kmino_ == 1){
				if(yh>valh_bubble[1][2]){
					valh_bubble[0][2] = pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)];
					valh_bubble[1][2] = yh;
					valh_bubble[2][2] = pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)];
					for(auto q=0; q<3; q++){
					valh_bubble[q+3][2]=pgv->Gn[ic].V[q]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])\
					*pgv->rdxyz[1]*(pgv->Gn[ic2-1].V[q]-pgv->Gn[ic].V[q]);
					}
					valh_bubble[6][2]=kappa_cv[ic]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
											pgv->rdxyz[1]*(kappa_cv[ic2-1]-kappa_cv[ic]);
				}
			}
			if(i+pgv->b->imino_==pgv->ijkm_gl[0] && k+pgv->b->kmino_ == pgv->ijkm_gl[2]){
				if(yh>valh_bubble[1][3]){
					valh_bubble[0][3] = pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)];
					valh_bubble[1][3] = yh;
					valh_bubble[2][1] = pgv->zc[k+pgv->b->kmino_-(1-pgv->nghost)];
					for(auto q=0; q<3; q++){
					valh_bubble[q+3][3]=pgv->Gn[ic].V[q]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])\
					*pgv->rdxyz[1]*(pgv->Gn[ic2-1].V[q]-pgv->Gn[ic].V[q]);
					}
					valh_bubble[6][3]=kappa_cv[ic]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
					pgv->rdxyz[1]*(kappa_cv[ic2-1]-kappa_cv[ic]);
				}
			}
		}
	}
}
valh[1] *= pgv->cell_volume;
bufh.resize(7);
for(auto i=0; i<7; i++) bufh[i].resize(pp->nprocs);
bufhh = new double**[7];
for(auto i=0; i<7; i++){
	bufhh[i] = new double*[4];
	for(auto j=0; j<4; j++){
		bufhh[i][j] = new double[pp->nprocs];
	}
}
pp->ierr=MPI_Gather(valh_min.data(), valh_min.size(),MPI_DOUBLE,bufh.data(),\
						  valh_min.size(), MPI_DOUBLE, 0, pp->MY_G_WORLD);
if(pp->myrank==0){
	val_min[1] = std::numeric_limits<double>::max();
	for(auto i=0; i<pp->nprocs; i++){
		if(bufh[1][i] < val_min[1]){
			for(auto q=0; q<val_min.size(); q++) val_min[q]=bufh[q][i];
		}
	}
}
pp->ierr=MPI_Gather(valh_max.data(), valh_max.size(), MPI_DOUBLE,bufh.data(),\
							valh_max.size(), MPI_DOUBLE, 0, pp->MY_G_WORLD);
if(pp->myrank==0){
	val_max[1] = -DBL_MAX;
	for(auto i=0; i<pp->nprocs; i++){
		if(bufh[1][i] > val_max[1]){
			for(auto q=0; q<val_max.size(); q++) val_max[q]=bufh[q][i];
		}
	}
}
auto si = 7*4*pp->nprocs;
pp->ierr=MPI_Gather(valh_bubble.data(),28,MPI_DOUBLE,bufhh,si,MPI_DOUBLE,0,pp->MY_G_WORLD);
if(pp->myrank==0){
	for(auto i=0; i<7; i++){
		for(auto j=0; j<4; j++) valh_bubble[i][j] = 0.0;
	}
	for(auto q=0; q<valh_bubble[q].size(); q++){
		valh_bubble[1][q] = -std::numeric_limits<double>::max();
	}
	for(auto i=0; i<pp->nprocs; i++){
		for(auto j=0; j<4; j++){
			if(bufhh[1][j][i]>valh_bubble[1][j]){
				for(auto q=0; q<valh_bubble[q].size(); q++){
					valh_bubble[q][j] = bufhh[q][j][i];
				}
			}
		}
	}
	for(auto i=0; i<7; i++){
		val_bubble[i]=std::accumulate(valh_bubble[i].begin(),valh_bubble[i].end(),0.0)/4.0;
	}
}
pp->ierr=MPI_Gather(valh_spike.data(),28,MPI_DOUBLE,bufhh,28,MPI_DOUBLE,0,pp->MY_G_WORLD);
if(pp->myrank==0){
	for(auto i=0; i<7; i++){
		for(auto j=0; j<4; j++) valh_spike[i][j] = 0.0;
	}
	for(auto q=0; valh_spike[0].size(); q++){
		valh_spike[1][q] = std::numeric_limits<double>::max();
	}
	for(auto i=0; i<pp->nprocs; i++){
		for(auto j=0; j<4; j++){
			if(bufhh[1][j][i]<valh_spike[1][j]){
				for(auto q=0; q<valh_spike.size(); q++){
					valh_spike[q][j] = bufhh[q][j][i];
				}
			}
		}
	}
	for(auto i=0; i<7; i++){
		val_spike[i]=std::accumulate(valh_spike[i].begin(),valh_spike[i].end(),0.0)/4.0;
	}
}
pp->ierr=MPI_Reduce(&valh[0],&val[0],valh.size(),MPI_DOUBLE,MPI_SUM,0,pp->MY_G_WORLD);
plbuf->freeR1Buffer(kappa_cv);
bufh.clear();
for(auto i=0; i<7; i++){
	for(auto j=0; j<4; j++){
		delete [] bufhh[i][j];
	}
	delete[] bufhh[i];
}
delete [] bufhh;

};
void toolbox::GaussPivotLarge(vector<vector<double>>&A,vector<double>&b,const int n,\
									   vector<double>&x){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

vector<vector<double>> Ab;
vector<double> rh(n+1);
double M;
Ab = A;
Ab.push_back(b);
for(auto j=0; j<n-1; j++){
	auto kh = fabs(Ab[j][j]);
	for(auto q=0; q<n+1; q++){
		if(kh<fabs(Ab[j+q][j])) kh = fabs(Ab[j+q][j]);
	}
	auto k = kh -1 + j;
	rh = Ab[k];
	Ab[j] = Ab[k];
	Ab[k] = rh;
	for(auto i=j+1; i<n; i++){ // checked
		M = Ab[i][j]/Ab[j][j];
		for(auto q=j; q<Ab[0].size(); q++){
			Ab[i][j+q] -= M*Ab[j][j+q];
		}
	}
}
std::fill(x.begin(), x.end(), 0.0);
x[n] = Ab[n-1][n]/Ab[n-1][n-1];
for(auto i=n-2; i >= 0; i--){
	x[i]=	(Ab[i][n]-std::inner_product(&Ab[i][i+1], &Ab[i][i+1]+n, &x[i+1],0.0))/Ab[i][i];
}
}
void toolbox::diag_RMI(vector<double> &val){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double dir, R, R1, R2, xwidth, xmid;
std::array<double, 2> xyc;
double xMinh, xMaxh, kappaXMinh, kappaXMaxh, sumIndicatorMinh, sumIndicatorMaxh;
double sumIndicatorMin, sumIndicatorMax, areah, area, xh, yh, zh;
double *kappa_cv, *v_cv;
int ncs, ncb;
std::vector<double> rhs, aa;
std::array<std::array<double, 3>, 4> valh;
double ***bufhh;
int indx;
std::vector<std::vector<double>> yyh, yy, xy, Bh, Ah;
std::fill(val.begin(), val.end(), 0.0);
for(auto i=0; i<4; i++) 	std::fill(valh[i].begin(), valh[i].end(), 0.0);
valh[0][0] = -std::numeric_limits<double>::max();
valh[0][1] = -std::numeric_limits<double>::max();
valh[0][2] =  std::numeric_limits<double>::max();
valh[0][3] =  std::numeric_limits<double>::max();
ncs = 0;
ncb = 0;
size_t kappa_cv_S;
kappa_cv = plbuf->getR1Buffer(kappa_cv_S);
size_t v_cv_S;
v_cv = plbuf->getR1Buffer(v_cv_S);
curvature(kappa_cv_S,kappa_cv);
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinN; ibl++) v_cv[pgv->ioff+ic]=pgv->Gn[ic].V[1];
}
pbou->updateGhostR1(v_cv);
yyh.resize(2); yy.resize(2);
for(auto i=0; i<2; i++){
	yyh[i].resize(pgv->ijkm_gl[0]);
	yy.resize(pgv->ijkm_gl[0]);
}
std::fill(yyh[0].begin(), yyh[0].end(), std::numeric_limits<double>::max());
std::fill(yy[1].begin(), yy[1].end(), -std::numeric_limits<double>::max());
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinT; ic++){
		auto i = pgv->Gn[ic].ijk[0] - pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1] - pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2] - pgv->b->kmino_;
		auto ic2 = pgv->i2c[i][j+1][k] - 1;
		if(pgv->Gn[ic].G*pgv->Gn[ic2].G< 0.0){
			if(fabs(pgv->Gn[ic2].G-pgv->Gn[ic].G)<1.0e-10){
				for(auto q=0; q<yyh.size(); q++){
					yyh[q][i] = 0.5*(pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]+\
										  pgv->yc[j+1+pgv->b->jmino_-(1-pgv->nghost)]);
				}
			}
			else{
				for(auto q=0; q<yyh.size(); q++){
					yyh[q][i] = pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]-pgv->Gn[ic].G*\
									pgv->dxyz[1]/(pgv->Gn[ic2].G-pgv->Gn[ic].G);
				}
			}
		}
	}
}
pp->ierr=MPI_Reduce(yyh.data(),yy.data(),yyh[0].size(),MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
pp->ierr=MPI_Reduce(&yyh[1][0],&yy[1][0],yyh[0].size(),MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
if(pp->myrank==0){
	xy.resize(2);
	for(auto q=0; q<2; q++){
		xy[q].resize(pgv->ijkm_gl[0]);
	}
	xmid = (pgv->xyze_sg[0]-pgv->xyzs_sg[0])*0.5+pgv->xyzs_sg[0];
	xwidth = 2.0/64.0*(pgv->xyze_sg[0]-pgv->xyzs_sg[0])*0.5;
	auto nx = -1;
	for(auto i=0; i<pgv->ijkm_gl[0]; i++){
		if(pgv->xc[i+pgv->nghost]>=xmid-xwidth && pgv->xc[i+pgv->nghost] <= xmid+xwidth){
			nx += 1;
			xy[0][nx] = pgv->xc[i+pgv->nghost]-xmid;
			xy[1][nx] = yy[0][i];
		}
	}
	auto ppp = xy.data();
	auto temp = std::accumulate(&xy[0][0], &xy[0][0]+nx+1,0.0);
	if(fabs(temp)>0.01*pgv->dxyz[0]){
		std::cout<<"Error!point data for bubble polynomial fit in RMI is not centered"<<std::endl;
		std::cout<<"xi=";
//		for(auto q=0; q<nx; q++) std::cout<<xy[0][q]<<" ";
		pp->parallel_kill(11);
	}
	Bh.resize(nx+1); Ah.resize(5);rhs.resize(5); aa.resize(5);
	for(auto q=0; q<nx+1; q++) Bh[q].resize(5);
	for(auto i=0; i<5; i++){
		for(auto j=0; j<5; j++){
			Bh[i][j] = pow(xy[0][i],j);
		}
	}
	for(auto q=0; q<5; q++) Ah.resize(5);
	for(auto i=0; i<5; i++){
		for(auto j=0; j<5; j++){
			auto temp = 0.0;
			for(auto q=0; q<Bh.size(); q++) temp+= Bh[q][i]*Bh[q][j];
			Ah[i][j] = temp;
		}
		temp = 0.0;
		for(auto q=0; q<nx+1; q++) temp += Bh[q][i]*xy[1][q];
		rhs[i] = temp;
	}
	GaussPivotLarge(Ah, rhs,5,aa);
	val[6] = aa[2];
	for(auto q=13; q<18; q++) val[q] = aa[q-13];
	val[18] = 0.0;
	for(auto i=0; i<nx+1; i++){
		val[18]+=pow(xy[1][i]-(aa[0]+aa[1]*xy[0][i]+aa[2]*pow(xy[0][i],2)+aa[3]*pow(xy[0][i],3)+aa[4]*pow(xy[0][i],4)),2);
	}
	val[18] /= static_cast<double>(nx+1);
	Bh.clear(); Ah.clear(); rhs.clear(); aa.clear();
	nx = -1;
	for(auto i=0; i<pgv->ijkm_gl[0]; i++){
		if(pgv->xc[i+pgv->nghost] <= pgv->xyzs_sg[0]+xwidth){
			nx += 1;
			xy[0][nx] = pgv->xc[i+pgv->nghost];
			xy[1][nx] = yy[0][i];
		}
		else if(pgv->xc[i+pgv->nghost]>=pgv->xyze_sg[0]-xwidth){
			nx += 1;
			xy[0][nx] = pgv->xc[i+pgv->nghost]-(pgv->xyze_sg[0]-pgv->xyzs_sg[0]);
			xy[1][nx] = yy[0][i];
		}
	}
	temp = std::accumulate(&xy[0][0],&xy[0][0]+nx+1, 0.0);
	if(fabs(temp)> 0.01 * pgv->dxyz[0]){
		std::cout<<"Error! point data for spike polynomial fit in RMI is not centered";
//		for(auto q=0; q<nx+1; q++) std::cout<<"xi="<<xy[0][q]<<" "<<std::endl;
		pp->parallel_kill(11);
	}
	Bh.resize(nx+1); Ah.resize(5); rhs.resize(5); aa.resize(5);
	for(auto q=0; q<nx+1; q++) Bh[q].resize(5);
	for(auto q=0; q<5; q++) Ah[q].resize(5);
	for(auto i=0; i<5; i++){
		for(auto j=0; j<5; j++)	Bh[i][j] = pow(xy[0][i],j);
	}
	for(auto i=0; i<5; i++){
		for(auto j=0; j<5; j++){
			auto temp = 0.0;
			for(auto q=0; q<Bh.size(); q++) temp+= Bh[q][i]*Bh[q][j];
			Ah[i][j] = temp;
		}
		temp = 0.0;
		for(auto q=0; q<nx+1; q++) temp += Bh[q][i]*xy[1][q];
		rhs[i] = temp;
	}
	GaussPivotLarge(Ah, rhs, 5, aa);
	val[3] = aa[2];
	for(auto q=7; q<12; q++) val[1] = aa[q-7];
	val[12] = 0.0;
	for(auto i=0; i<nx+1; i++){
		val[13]+=pow(xy[1][i]-(aa[0]+aa[1]*xy[0][i]+aa[2]*pow(xy[0][i],2)+aa[3]*pow(xy[0][i],3)+aa[4]*pow(xy[0][i],4)),2);
	}
	val[12] /= static_cast<double>(nx);
	std::ofstream myfile;
	myfile.open("circle_bubble");
	if(myfile.is_open()){
		for(auto i=0; i<nx; i++){
			myfile<<xy[0][i]<<'\t'<<xy[1][i]<<'\t'<<aa[0]+aa[1]*xy[0][i]+aa[2]*pow(xy[0][i],2)+\
					  aa[3]*pow(xy[0][i],3)+aa[4]*pow(xy[0][i],4)<<std::endl; // correct format
		}
	}
	else std::cout<<"the file circle_bubble is toolbox cannot be opened"<<std::endl;
	Bh.clear(); Ah.clear(); rhs.clear(); aa.clear(); xy.clear();
}
yyh.clear(); yy.clear();
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinT; ic++){
		auto i = pgv->Gn[ic].ijk[0] - pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1] - pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2] - pgv->b->kmino_;
		if((i+pgv->b->imino_==1 || i+pgv->b->imino_==pgv->ijkm_gl[0])||((pgv->ijkm_gl[0]%2==0) &&\
			(i+pgv->b->imino_==pgv->ijkm_gl[0]/2 || i+pgv->b->imino_==pgv->ijkm_gl[0]/2+1)) || \
			((pgv->ijkm_gl[0]%2 != 0)&&(i+pgv->b->imino_==(pgv->ijkm_gl[0]+1)/2))){
			auto ic2 = pgv->i2c[i][j+1][k] - 1;
			if(pgv->Gn[ic].G*pgv->Gn[ic2].G < 0.0){
				if(fabs(pgv->Gn[ic2].G-pgv->Gn[ic].G)<1.0e-10){
					yh = 0.5*(pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]+\
								 pgv->yc[j+1+pgv->b->jmino_-(1-pgv->nghost)]);
				}
				else{
					yh=pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]-pgv->Gn[ic].G*\
						pgv->dxyz[1]/(pgv->Gn[ic2].G-pgv->Gn[ic].G);
				}
				if(i+pgv->b->imino_==1){
					indx = 0;
					dir = 1.0;
				}
				else if(i+pgv->b->imino_==pgv->ijkm_gl[0]){
					indx = 1;
					dir = 1.0;
				}
				else if((pgv->ijkm_gl[0]%2==0 && i+pgv->b->imino_==pgv->ijkm_gl[0]/2) || \
							(pgv->ijkm_gl[0]/2!=0 && i+pgv->b->imino_==(pgv->ijkm_gl[0]+1)/2)){
					indx = 2;
					dir = -1.0;
				}
				else if(pgv->ijkm_gl[0]%2==0 && i+pgv->b->imino_==pgv->ijkm_gl[0]/2+1){
					indx = 3;
					dir = -1.0;
				}
				if(dir*yh > dir*valh[0][indx]){
					valh[0][indx] = yh;
					valh[1][indx] = v_cv[pgv->ioff+ic]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])*\
										 pgv->rdxyz[1]*(v_cv[pgv->ioff+ic2]-v_cv[pgv->ioff+ic]);
					valh[2][indx]=kappa_cv[pgv->ioff+ic]+(yh-pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)])\
									  *pgv->rdxyz[1]*(kappa_cv[pgv->ioff+ic2]-kappa_cv[pgv->ioff+ic]);
				}
			}
		}
	}
}
bufhh = new double**[3];
for(auto i=0; i<3; i++){
	bufhh[i] = new double*[4];
	for(auto j=0; j<4; j++){
		bufhh[i][j] = new double[pp->nprocs];
	}
}
auto s1 = valh.size()*valh[0].size();
pp->ierr=MPI_Gather(valh.data(),s1,MPI_DOUBLE,bufhh,s1,MPI_DOUBLE,0,MPI_COMM_WORLD);
if(pp->myrank==0){
	for(auto q=0; q<valh.size(); q++) std::fill(valh[q].begin(),valh[q].end(),0.0);
	valh[0][0] = -std::numeric_limits<double>::max();
	valh[0][1] = -std::numeric_limits<double>::max();
	valh[0][2] =  std::numeric_limits<double>::max();
	valh[0][3] =  std::numeric_limits<double>::max();
	for(auto i=0; i<pp->nprocs; i++){
		for(auto j=0; j<2; j++){
			if(bufhh[0][j][i] > valh[0][j]){
				for(auto q=0; q<pp->nprocs; q++){
					valh[q][j] = bufhh[q][j][i];
				}
			}
		}

	}	
	val[0] = 0.25*(valh[0][0]+valh[0][1]-valh[0][2]-valh[0][3]);
	val[1] = 0.5*(valh[1][0]+valh[1][1]);
	val[2] = 0.5*(valh[2][0]+valh[2][1]);
	val[4] = 0.5*(valh[1][2]+valh[1][3]);
	val[5] = 0.5*(valh[2][2]+valh[2][3]);
	val[1] *= 2.0*pgv->pi/pgv->init_wavelength;
	for(auto q=2; q<4; q++){
		val[q]=fabs(val[q]*pgv->init_wavelength/(2.0*pgv->pi));
	}
	for(auto q=5; q<7; q++){
		val[q]=fabs(val[q]*pgv->init_wavelength/(2.0*pgv->pi));
	}
}
plbuf->freeR1Buffer(v_cv);
plbuf->freeR1Buffer(kappa_cv);
for(auto i=0; i<3; i++){
	for(auto j=0; j<4; j++){
		delete [] bufhh[i][j];
	}
	delete [] bufhh[i];
}
delete [] bufhh;
};
void toolbox::diag_ligament(vector<double>&val){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

std::array<double, 7> myval;
// myval is already zero; // checked
myval[0] =  std::numeric_limits<double>::max();
myval[1] = -std::numeric_limits<double>::max();
myval[2] =  std::numeric_limits<double>::max();
myval[3] = -std::numeric_limits<double>::max();
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinT; ic++){
		auto i = pgv->Gn[ic].ijk[0] - pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1] - pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2] - pgv->b->kmino_;
		if(j+pgv->b->jmino_==pgv->ijkm_gl[1]/2){
			auto ic2 = pgv->i2c[i+1][j][k] - 1;
			if(pgv->Gn[ic].G*pgv->Gn[ic2].G <= 0.0){
				if(fabs(pgv->Gn[ic2].G-pgv->Gn[ic].G < 1.0e-10)){
					myval[0] = std::min(myval[0], 0.5*(pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)]+\
											  pgv->xc[i+1+pgv->b->imino_-(1-pgv->nghost)]));
					myval[1] = std::max(myval[1], 0.5*(pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)]+\
											  pgv->xc[i+1+pgv->b->imino_-(1-pgv->nghost)]));
				}
				else{
					auto temp=pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)]-pgv->Gn[ic].G*\
											pgv->dxyz[0]/(pgv->Gn[ic2].G-pgv->Gn[ic].G);
					myval[0] = std::min(myval[0], temp);
					myval[1] = std::max(myval[1], temp);
				}
			}
		}
		if(i+pgv->b->imino_==pgv->ijkm_gl[0]/2){
			auto ic2 = pgv->i2c[i][j+1][k] - 1;
			if(pgv->Gn[ic].G-pgv->Gn[ic].G <= 0.0){
				if(fabs(pgv->Gn[ic2].G-pgv->Gn[ic].G)<1.0e-10){
					myval[2] = std::min(myval[2],0.5*(pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]+\
											  pgv->yc[j+1+pgv->b->jmino_-(1-pgv->nghost)]));
					myval[3] = std::min(myval[3],0.5*(pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]+\
											  pgv->yc[j+1+pgv->b->jmino_-(1-pgv->nghost)]));
				}
				else{
					auto temp=pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]-pgv->Gn[ic].G*\
								 pgv->dxyz[1]/(pgv->Gn[ic2].G-pgv->Gn[ic].G);
					myval[2] = std::min(myval[2],temp);
					myval[3] = std::max(myval[3],temp);
				}
			}
		}
	}
}
myval[0] = -myval[0];
myval[2] = -myval[2];
pp->ierr=MPI_Reduce(&myval[0],&val[0],myval.size(),MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
if(pp->myrank==0){
	val[0] = -val[0];
	val[2] = -val[2];
}
}

void toolbox::diag_drop(std::vector<double> &val){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

std::array<double, 7> myval;
// myval is already zero;
myval[0] =  std::numeric_limits<double>::max();
myval[1] = -std::numeric_limits<double>::max();
myval[2] =  std::numeric_limits<double>::max();
myval[3] = -std::numeric_limits<double>::max();
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinT; ic++){
		auto i = pgv->Gn[ic].ijk[0] - pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1] - pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2] - pgv->b->kmino_;
		if(j+pgv->b->jmino_==pgv->ijkm_gl[1]/2){
			auto ic2 = pgv->i2c[i+1][j][k] - 1;
			if(pgv->Gn[ic].G*pgv->Gn[ic2].G <= 0.0){
				if(fabs(pgv->Gn[ic2].G-pgv->Gn[ic].G)<1.0e-10){
					myval[0]=std::min(myval[0],0.5*(pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)]+\
								pgv->xc[i+1+pgv->b->imino_-(1-pgv->nghost)]));
					myval[1]=std::max(myval[1],0.5*(pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)]+\
								pgv->xc[i+1+pgv->b->imino_-(1-pgv->nghost)]));
				}
				else{
					auto temp=pgv->xc[i+pgv->b->imino_-(1-pgv->nghost)]-pgv->Gn[ic].G*\
											pgv->dxyz[0]/(pgv->Gn[ic2].G-pgv->Gn[ic].G);
					myval[0]=std::min(myval[0],temp);
					myval[1]=std::max(myval[1],temp);
				}
			}
		}
		if(i+pgv->b->imino_==pgv->ijkm_gl[0]/2){
			auto ic2 = pgv->i2c[i][j+1][k] - 1;
			if(pgv->Gn[ic].G*pgv->Gn[ic2].G <= 0.0){
				if(fabs(pgv->Gn[ic2].G-pgv->Gn[ic].G) < 1.0e-10){
					myval[2]=std::min(myval[2],0.5*(pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]+\
								pgv->yc[j+1+pgv->b->jmino_-(1-pgv->nghost)]));
					myval[3]=std::min(myval[3],0.5*(pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]+\
								pgv->yc[j+1+pgv->b->jmino_-(1-pgv->nghost)]));
				}
				else{
 					auto temp=pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]-pgv->Gn[ic].G*\
								 pgv->dxyz[1]/(pgv->Gn[ic2].G-pgv->Gn[ic].G);
					myval[2]=std::min(myval[2],temp);
					myval[3]=std::max(myval[3],temp);
				}
			}
		}
	}
}
myval[0] = -myval[0];
myval[2] = -myval[2];
pp->ierr=MPI_Reduce(&myval[0],&val[0],myval.size(),MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
if(pp->myrank==0){
	val[0] = -val[0];
	val[2] = -val[2];
}
val[4] = volume();
}
void toolbox::diagnostics_init(){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);
if(!pgv->do_diagnostics) return;
if(pgv->diagnostics_type=="RMI"){
	pmo->lit_monitor_create_file_step("lit_RMI",19);
	pmo->lit_monitor_set_header(1,"a",'r');
	pmo->lit_monitor_set_header(2,"v_b",'r');
	pmo->lit_monitor_set_header(3,"kappa_b",'r');
	pmo->lit_monitor_set_header(4,"zeta1_b/k",'r');
	pmo->lit_monitor_set_header(5,"v_s",'r');
	pmo->lit_monitor_set_header(6,"kappa_s",'r');
	pmo->lit_monitor_set_header(7,"zeta1_s/k",'r');
	pmo->lit_monitor_set_header(8,"a0_b",'r');
	pmo->lit_monitor_set_header(9,"a1_b",'r');
	pmo->lit_monitor_set_header(10,"a2_b",'r');
	pmo->lit_monitor_set_header(11,"a3_b",'r');
	pmo->lit_monitor_set_header(12,"a4_b",'r');
	pmo->lit_monitor_set_header(13,"E_b^2./n_b",'r');
	pmo->lit_monitor_set_header(14,"a0_s",'r');
	pmo->lit_monitor_set_header(15,"a1_s",'r');
	pmo->lit_monitor_set_header(16,"a2_s",'r');
	pmo->lit_monitor_set_header(17,"a3_s",'r');
	pmo->lit_monitor_set_header(18,"a4_s",'r');
	pmo->lit_monitor_set_header(19,"E/s^2/n_s",'r');
}
else if(pgv->diagnostics_type=="RT3D"){
	pmo->lit_monitor_create_file_step("lit_RT3D",3);
	pmo->lit_monitor_set_header(1,"A spike/bubble",'r');
	pmo->lit_monitor_set_header(2,"A min/max",'r');
	pmo->lit_monitor_set_header(3,"Surface_area",'r');
	pmo->lit_monitor_create_gnuplot(1);
	pmo->lit_monitor_create_file_step("lit_RT3D_spike",12);
	pmo->lit_monitor_set_header(1,"y_spike",'r');
	pmo->lit_monitor_set_header(2,"U_spike",'r');
	pmo->lit_monitor_set_header(3,"V_spike",'r');
	pmo->lit_monitor_set_header(4,"W_spike",'r');
	pmo->lit_monitor_set_header(5,"kappa_spike",'r');
	pmo->lit_monitor_set_header(6,"x_min",'r');
	pmo->lit_monitor_set_header(7,"y_min",'r');
	pmo->lit_monitor_set_header(8,"z_min",'r');
	pmo->lit_monitor_set_header(9,"U_min",'r');
	pmo->lit_monitor_set_header(10,"V_min",'r');
	pmo->lit_monitor_set_header(11,"W_min",'r');
	pmo->lit_monitor_set_header(12,"kappa_min",'r');
	pmo->lit_monitor_create_file_step("lit_RT3D_bubble",12);
	pmo->lit_monitor_set_header(1,"y_bubble",'r');
	pmo->lit_monitor_set_header(2,"U_bubble",'r');
	pmo->lit_monitor_set_header(3,"V_bubble",'r');
	pmo->lit_monitor_set_header(4,"W_bubble",'r');
	pmo->lit_monitor_set_header(5,"kappa_bubble",'r');
	pmo->lit_monitor_set_header(6,"x_max",'r');
	pmo->lit_monitor_set_header(7,"y_max",'r');
	pmo->lit_monitor_set_header(8,"z_max",'r');
	pmo->lit_monitor_set_header(9,"U_max",'r');
	pmo->lit_monitor_set_header(10,"V_max",'r');
	pmo->lit_monitor_set_header(11,"W_max",'r');
	pmo->lit_monitor_set_header(12,"kappa_max",'r');
}
else if(pgv->diagnostics_type=="ligament"){
	pmo->lit_monitor_create_file_step("lit_ligamenet",6);
	pmo->lit_monitor_set_header(1,"xmin",'r');
	pmo->lit_monitor_set_header(2,"xmax",'r');
	pmo->lit_monitor_set_header(3,"ymin",'r');
	pmo->lit_monitor_set_header(4,"ymax",'r');
	pmo->lit_monitor_set_header(5,"a",'r');
	pmo->lit_monitor_set_header(6,"b",'r');
}
else if(pgv->diagnostics_type=="DSW"){
	pmo->lit_monitor_create_file_step("lit_dsw",1);
	pmo->lit_monitor_set_header(1,"A",'r');
}
else if(pgv->diagnostics_type=="drop"){
	pmo->lit_monitor_create_file_step("lit_drop",5);
	pmo->lit_monitor_set_header(1,"xmin",'r');
	pmo->lit_monitor_set_header(2,"xmax",'r');
	pmo->lit_monitor_set_header(3,"ymin",'r');
	pmo->lit_monitor_set_header(4,"ymax",'r');
	pmo->lit_monitor_set_header(5,"V",'r');
}
else{
	pp->parallel_die("Unknown diagnostics type: ");
}
}

void toolbox::diagnostics(){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

bool first_call = true;
std::array<double, 3> velXmin, velXmax;
std::vector<double> val_spike(19),val_bubble(19),val_min(19),val_max(19),val(19);
double xmin, xmax, kappaXmin, kappaXmax, area, A;
bool exists;
if(!pgv->do_diagnostics) return;
if(pgv->diagnostics_type=="RMI"){
	diag_RMI(val);
	pmo->lit_monitor_select_file("lit_RMI");
	pmo->lit_monitor_set_single_values(val[0],val[1],val[2],val[3],val[4],val[5],val[6],\
												  val[7],val[8],val[9],val[10],val[11],val[12],val[13],\
												  val[14],val[15],val[16],val[17],val[18]);
}
else if(pgv->diagnostics_type=="RT3D"){
	diag_RT3D(val,val_spike,val_min,val_bubble,val_max);
	pmo->lit_monitor_select_file("lit_RT3D");
	pmo->lit_monitor_set_single_values(0.5*(val_bubble[1]-val_spike[1]),\
										0.5*(val_max[1]-val_min[1]),val[1]);
	pmo->lit_monitor_select_file("lit_RT3D_spike");
	pmo->lit_monitor_set_single_values(val_spike[1],val_spike[3],val_spike[4],val_spike[5],\
												  val_spike[6],val_min[0],val_min[1],val_min[2],\
												  val_min[3],val_min[4],val_min[5],val_min[6]);
	pmo->lit_monitor_select_file("lit_RT3D_bubble");
	pmo->lit_monitor_set_single_values(val_bubble[1],val_bubble[3],val_bubble[4],\
												  val_bubble[5],val_bubble[6],val_max[0],val_max[1],\
												  val_max[2],val_max[3],val_max[4],val_max[5],val_max[6]);
}
else if(pgv->diagnostics_type=="ligament"){
	diag_ligament(val);
	pmo->lit_monitor_select_file("lit_ligament");
	pmo->lit_monitor_set_single_values(val[0],val[1],val[2],val[3],0.5*(val[1]-val[0]),\
												  0.5*(val[3]-val[2]));
}
else if(pgv->diagnostics_type=="DSW"){
	A = calcAmplitude_quarter();
	pmo->lit_monitor_select_file("lit_dsw");
	pmo->lit_monitor_set_single_values(A);
}
else if(pgv->diagnostics_type=="drop"){
	diag_drop(val);
	pmo->lit_monitor_set_single_values(val[0],val[1],val[2],val[3],val[4]);
}
else{
pp->parallel_die("Unknown diagnostics type: ");
}
}

void toolbox::checkBandLayerGmin(){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double bandGmin;
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	auto si=pgv->band_size[0]+pgv->band_size[1]+pgv->band_size[2];
	for(auto il=0; il<si; il++){
		auto temp=static_cast<double>(il-2.0);
		temp *= *std::min_element(pgv->dxyz.begin(),pgv->dxyz.end())/sqrt(3.0);
		bandGmin=std::min(pgv->G_max, temp);
		for(auto ic=0; ic<pgv->b->NinBandLayer[il]; ic++){
			if(fabs(pgv->Gn[ic].G<bandGmin-1.0e-10)){
				std::cout<<pgv->clit<<'\t'<<pp->myrank<<"WARNING: G smaller than minimum sensible value in band layer"<<std::endl; // format needs to be corrected
				std::cout<<pgv->clit<<'\t'<<pp->myrank<<"    layer#,x,y,z,G,sensible Gmin="<<il<<'\t';
				std::cout<<pgv->xc[pgv->Gn[ic].ijk[0]-(1-pgv->nghost)]<<'\t';
				std::cout<<pgv->yc[pgv->Gn[ic].ijk[1]-(1-pgv->nghost)]<<'\t';
				std::cout<<pgv->zc[pgv->Gn[ic].ijk[2]-(1-pgv->nghost)]<<std::endl; // correct the format
				std::cout<<pgv->clit<<'\t'<<pp->myrank<<"    ijk="<<pgv->Gn[ic].ijk[0]<<'\t';
				std::cout<<pgv->Gn[ic].ijk[1]<<'\t'<<pgv->Gn[ic].ijk[2]<<std::endl; //correct format
				std::cout<<pgv->clit<<'\t'<<pp->myrank<<"    ibl,ijk_sg="<<ibl<<'\t';
				std::cout<<pgv->b->ijk_sg[0]<<'\t'<<pgv->b->ijk_sg[1]<<'\t'<<pgv->b->ijk_sg[2]<<std::endl; // correct format
				std::cout<<pgv->clit<<'\t'<<pp->myrank<<"    il,band_size(0 to 4)"<<'\t';
				std::cout<<il<<'\t'<<pgv->band_size[0]<<'\t'<<pgv->band_size[1]<<'\t'<<pgv->band_size[2]<<'\t';
				std::cout<<pgv->band_size[3]<<'\t'<<pgv->band_size[4]<<std::endl; //correct format
				std::cout<<pgv->clit<<'\t'<<pp->myrank<<"    NinBandLayer(0 to 5)=";pgv->b->NinBandLayer[0]<<'\t';
				std::cout<<pgv->b->NinBandLayer[1]<<'\t'<<pgv->b->NinBandLayer[2]<<'\t';
				std::cout<<pgv->b->NinBandLayer[3]<<'\t'<<pgv->b->NinBandLayer[4]<<'\t';
				std::cout<<pgv->b->NinBandLayer[5]<<std::endl; // correct format
				std::cout<<pgv->clit<<'\t'<<pp->myrank<<"    ic,NinA,NinT,NinN,ic_max="<<ic<<'\t';
				std::cout<<pgv->b->NinA<<'\t'<<pgv->b->NinT<<'\t'<<pgv->b->NinN<<'\t';
				std::cout<<pgv->b->NinBandLayer[il]<<std::endl; // correct format
			}
		}
	}
}
}
void toolbox::rotating_cyl_mfix_init(){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

if(!pgv->rotating_cyl_mfix) return;
// should be init_mod-1; checked
pgv->init_volume=pgv->pi*pow(pgv->xyze_sg[1]-pgv->xyzs_sg[1],2)*(pgv->\
					  init_center[pgv->init_mod-1]-pgv->xyzs_sg[pgv->init_mod-1]);
pmo->lit_monitor_create_file_step("lit_mfix",4);
pmo->lit_monitor_set_header(1,"Init_Vol",'r');
pmo->lit_monitor_set_header(2,"Curr_Vol",'r');
pmo->lit_monitor_set_header(3,"Mass_Fix",'r');
pmo->lit_monitor_set_header(4,"G_shift",'r');
pmo->lit_monitor_select_file("lit_mfix");
pmo->lit_monitor_set_single_values(pgv->init_volume,pgv->init_volume,0.0,0.0);
}
void toolbox::rotating_cyl_mfix_adjust(double inVol){
std::unique_ptr<global_variable> pgv(new global_variable);
std::unique_ptr<parallel> pp(new parallel);
std::unique_ptr<bound> pbou(new bound);
std::unique_ptr<litBuffer> plbuf(new litBuffer);
std::unique_ptr<lit_timer_t> plt(new lit_timer_t);
std::unique_ptr<init> pinit(new init);
std::unique_ptr<monitor> pmo(new monitor);

double myVol, tVol, VolErr, inactive_cell_volume, inactive_cell_VOF;
double sgVol, G_shift, dr, ny, denom;
double *nx;
if(!pgv->rotating_cyl_mfix) return;
size_t nx_S;
nx = plbuf->getR1Buffer(nx_S);
if(inVol==-1.0){
	myVol = 0.0;
	for(auto ibl=0; ibl<pgv->nbl; ibl++){
		pgv->setBlockPointers(ibl);
		for(auto ic=0; ic<pgv->b->NinT; ic++){
			auto i = pgv->Gn[ic].ijk[0] - 1;
			auto j = pgv->Gn[ic].ijk[1] - 1;
			auto k = pgv->Gn[ic].ijk[2] - 1;
			myVol += cell_volume_fraction(pgv->Gn[ic],ic+1,pgv->b)*pgv->yc[j+pgv->nghost]*\
					 	pgv->cell_volume;
		}
		for(auto ic=pgv->b->NinT; ic>pgv->b->NinX; ic--){
			auto i = pgv->Gn[ic].ijk[0] - 1;
			auto j = pgv->Gn[ic].ijk[1] - 1;
			auto k = pgv->Gn[ic].ijk[2] - 1;
			auto temp = 0.5;
			if(pgv->Gn[ic].G<0) temp =-0.5;
			myVol += (0.5+temp)*pgv->yc[j+pgv->nghost]*pgv->cell_volume;
		}
		auto k = 1-pgv->b->kmino_;
		inactive_cell_VOF = -1.0;
		inactive_cell_volume = 0.0;
		for(auto i=pgv->b->imin_-pgv->b->imino_; i<pgv->b->imax_-pgv->b->imino_; i++){
			for(auto j=pgv->b->jmin_-pgv->b->jmino_; j<pgv->b->jmax_-pgv->b->jmino_; j++){
				if(pgv->i2c[i][j][k] > 0){
					auto temp = 0.5;
					if(pgv->Gn[pgv->i2c[i][j][k]-1].G<0) temp = -0.5;
					inactive_cell_VOF=0.5+temp;
					myVol += inactive_cell_volume*inactive_cell_VOF;
					inactive_cell_volume = 0.0;
					break;
				}
				else{
				 inactive_cell_volume=inactive_cell_volume+pgv->yc[j+pgv->b->jmino_-(1-pgv->nghost)]*\
											 pgv->cell_volume;
				}
			}
		}
		myVol += inactive_cell_volume*inactive_cell_VOF;
		inactive_cell_volume = 0.0;
	}
	if(pp->myrank==0){
		 auto k_sg = 0;
		for(auto i_sg = 0; i_sg<pgv->ijkm_sg[0]; i_sg++){
			for(auto j_sg=0; j_sg<pgv->ijkm_sg[1]; j_sg++){
				if(pgv->sg_active[i_sg][j_sg][k_sg]) break;
				auto t1 = static_cast<double>(j_sg);
				auto t2 = pgv->dxyz_sg[0]*pgv->dxyz_sg[1]*pgv->dxyz_sg[2];
				sgVol = (pgv->xyzs_sg[1]+0.5*t1*pgv->dxyz_sg[1])*t2;
				myVol += sgVol;
			}
		}
	}
	pp->parallel_all_sum(myVol, tVol);
}
else pp->parallel_all_sum(inVol, tVol);
VolErr = tVol-pgv->init_volume;
dr = pgv->xyze_sg[1]-pgv->xyzs_sg[1];
G_shift = VolErr/(pgv->pi*pow(dr,2));
for(auto i=0; i<nx_S; i++) nx[i] = 0.0; //3137 checked
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinT; ic++){
		auto i = pgv->Gn[ic].ijk[0] - pgv->b->imino_;
		auto j = pgv->Gn[ic].ijk[1] - pgv->b->jmino_;
		auto k = pgv->Gn[ic].ijk[2] - pgv->b->kmino_;
		nx[pgv->ioff+ic]=(pgv->Gn[pgv->i2c[i+1][j][k]-1].G-pgv->Gn[pgv->i2c[i-1][j][k]-1].G)*\
				  0.5*pgv->rdxyz[0];
		ny = (pgv->Gn[pgv->i2c[i][j+1][k]-1].G-pgv->Gn[pgv->i2c[i][j-1][k]-1].G)*\
				pgv->rdxyz[1];
		auto temp = sqrt(pow(nx[ic],2)+pow(ny,2));
		denom = std::max(temp, 1.0e-10);
		nx[pgv->ioff+ic] = denom/nx[pgv->ioff+ic];
	}
}
for(auto ibl=0; ibl<pgv->nbl; ibl++){
	pgv->setBlockPointers(ibl);
	for(auto ic=0; ic<pgv->b->NinN; ic++){
		pgv->Gn[ic].G -=G_shift;
	}
}
pbou->updateGhostNodes();
pmo->lit_monitor_select_file("lit_mfix");
pmo->lit_monitor_set_single_values(pgv->init_volume, tVol, VolErr, G_shift);
plbuf->freeR1Buffer(nx);
}
double toolbox::width_d = 0.0;
double toolbox::pi_wd = 0.0;
double toolbox::r2wd = 0.0;
double toolbox::r13 = 0.0;








