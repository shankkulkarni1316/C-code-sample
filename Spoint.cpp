
//Change it from heat transfer problem into mechanical problem. But dynamic

#include"Spoint.h"
//#include "Smatrix.h"

Point::Point(int tmpn){
	point_number=tmpn;
	x=0.0;y=0.0;
	x0=0.0;y0=0.0;
	dx=0.0;dy=0.0;
	rp=0.0;
	density=0.0;
	volume=0.0;
	mass=0.0;
	E_point=0.0;
	Bx=0.0; By=0.0;     //@Shank
	K_shape.Setdim(2,2); // set dimention of shape tensor K
	K_shape_inv.Setdim(2,2);
	F_dg.Setdim(2,2); // set dimention of deformation gradient F
	F_dg_1.Setdim(2,2);
	point_strain.Setdim(2,2);
	
	D_point.Setdim(2,2);
	//type=0;
	
	fx=0.0;fy=0.0;
	vx=0.0;vy=0.0;
	ax=0.0;ay=0.0;
	flag_fix_x=0;flag_fix_y=0;//defalt displacement is free
	flag_fixfx=0;flag_fixfy=0;//defalt force is free
	flag_fixvx=0;flag_fixvy=0;//defalt velocity is free
	flag_fixax=0;flag_fixay=0;//defalt acceleration is free
	
	//value for calculate forcestate
	theta=0.0;B=0.0;q=0.0;alpha=0.0;
	point_omega=0.0;

	Sxx= 0.0; Syy= 0.0; Sxy=0.0;
	Exx= 0.0; Eyy= 0.0; Exy=0.0;
	Exx_c= 0.0; Eyy_c= 0.0; Exy_c=0.0;
	broken_parameter=0.0;
}

Point::Point(){
	point_number=0;
	x=0.0;y=0.0;
	x0=0.0;y0=0.0;
	rp=0.0;
	density=0.0;
	volume=0.0;
	mass=0.0;
	E_point=0.0;
	Bx=0.0; By=0.0;      //@Shank
	K_shape.Setdim(2,2); // set dimention of shape tensor K
	K_shape_inv.Setdim(2,2);
	F_dg.Setdim(2,2); // set dimention of deformation gradient F
	F_dg_1.Setdim(2,2);
	point_strain.Setdim(2,2);
	
	D_point.Setdim(2,2);
	//type=0;
	
	fx=0.0;fy=0.0;
	vx=0.0;vy=0.0;
	dx=0.0;dy=0.0;
	ax=0.0;ay=0.0;
	flag_fix_x=0;flag_fix_y=0;//defalt displacement is free
	flag_fixfx=0;flag_fixfy=0;//defalt force is free
	flag_fixvx=0;flag_fixvy=0;//defalt velocity is free
	flag_fixax=0;flag_fixay=0;//defalt acceleration is free
	
	//value for calculate forcestate
	theta=0.0;B=0.0;q=0.0;alpha=0.0;


	point_omega=0.0;
	Sxx= 0.0; Syy= 0.0; Sxy=0.0;
	Exx= 0.0; Eyy= 0.0; Exy=0.0;
	Exx_c= 0.0; Eyy_c= 0.0; Exy_c=0.0;
	broken_parameter=0.0;
}

//settle the current coordinate of the atom
//it automatically update the disp
void Point::SetX(int n,double temp){
	if(flag_fix_x==0&&n==1){x=temp;dx=x-x0;}
	if(flag_fix_y==0&&n==2){y=temp;dy=y-y0;}
}
void Point::SetX(double tmpx,double tmpy){
	if (flag_fix_x==0){x=tmpx;dx=x-x0;}
	if(flag_fix_y==0){y=tmpy;dy=y-y0;}
}

//settle the displacement of the atom
//it automatically update the current coordinate
void Point::SetD(int n,double temp){
	if(flag_fix_x==0&&n==1){dx=temp;x=x0+dx;}
	if(flag_fix_y==0&&n==2){dy=temp;y=y0+dy;}
}
void Point::SetD(double tmpx,double tmpy){
	if (flag_fix_x==0){dx=tmpx;x=x0+dx;}
	if(flag_fix_y==0){dy=tmpy;y=y0+dy;}
}

//set the nodal acceleration, n=1 for x, 2 for y
void Point::Seta(int n,double tempa){
	if(flag_fixax==0&&n==1){ax=tempa;}
	if(flag_fixay==0&&n==2){ay=tempa;}
}
void Point::Seta(double tmpx,double tmpy){
	if(flag_fixax==0){ax=tmpx;}
	if(flag_fixay==0){ay=tmpy;}
}

//set the force of the node, only one direction
//n=1 for x 2 for y
void Point::SetF(int n,double temp){
	if(flag_fixfx==0&&n==1){fx=temp;}
	if(flag_fixfy==0&&n==2){fy=temp;}
}
void Point::SetF(double tmpx,double tmpy){
	if(flag_fixfx==0){fx=tmpx;}
	if(flag_fixfy==0){fy=tmpy;}
}

//set the body force density of the node, only one direction       //@Shank
//n=1 for x 2 for y
void Point::SetBF(int n,double temp){
	if(flag_fixBx==0&&n==1){Bx=temp;}
	if(flag_fixBy==0&&n==2){By=temp;}
}
void Point::SetBF(double tmpx,double tmpy){
	if(flag_fixBx==0){Bx=tmpx;}
	if(flag_fixBy==0){By=tmpy;}                                    //@Shank
}

//accumulate atom force. n=1 for x and 2 for y direction
//the force of a point is determined by many neighbors, so it is "accumulated"
void Point::AccumulateF(int n,double temp){
	if (flag_fixfx==0&&n==1){fx=fx+temp;}
	if (flag_fixfy==0&&n==2){fy=fy+temp;}
}
void Point::AccumulateF(double tmpx,double tmpy){
	if (flag_fixfx==0){fx=fx+tmpx;}
	if (flag_fixfy==0){fy=fy+tmpy;}
}

//set the point velocity, n=1 for x 2 for y direction
void Point::SetV(int n,double tempa){
	if(flag_fixvx==0&&n==1){vx=tempa;}
	if(flag_fixvy==0&&n==2){vy=tempa;}
}
void Point::SetV(double tmpx,double tmpy){
	if(flag_fixvx==0){vx=tmpx;}
	if(flag_fixvy==0){vy=tmpy;}
	//vx=tmpx;
	//vy=tmpy;
}

//settle the flag of t
/*void Atom::SetFlagT(int temp){
	flag_fixt=temp;
}*/

//settle the flag of displacement. fixed for 1, free for 0, n=1 for x, 2 for y direction
void Point::SetFlagX(int n,int temp){
	if(n==1){flag_fix_x=temp;}
	if(n==2){flag_fix_y=temp;}
}
/*void Point::SetFlagX(int tmpx,int tmpy){
	flag_fix_x=tmpx;
	flag_fix_y=tmpy;
}*/

//settle the flag of f, tempf can only be 0 or 1, n=1 for x, 2 for y direction
void Point::SetFlagF(int n,int tempf){
	if(n==1){flag_fixfx=tempf;}
	if(n==2){flag_fixfy=tempf;}
}

void Point::SetFlagBF(int n,int tempf){       //@Shank
	if(n==1){flag_fixBx=tempf;}
	if(n==2){flag_fixBy=tempf;}
}
/*void Point::SetFlagF(int tmpx,int tmpy){
	flag_fixfx=tmpx;
	flag_fixfy=tmpy;
}*/

//settle the flag of v
void Point::SetFlagV(int n,int temp){
	if(n==1){flag_fixvx=temp;}
	if(n==2){flag_fixvy=temp;}
}
/*void Point::SetFlagV(int tmpx,int tmpy){
	flag_fixvx=tempf;
	flag_fixvy=tempf;
}*/

//settle the flag of a
//n=1 for x, 2 for y direction
void Point::SetFlagA(int n,int temp){
	if(n==1){flag_fixax=temp;}
	if(n==2){flag_fixay=temp;}
}
/*void Pointm::SetFlagA(int tmpx,int tmpy){
	flag_fixax=tmpx;
	flag_fixay=tmpy;
}*/

//settle nodal mass
void Point::SetMass(){
	//mass=rp*rp*density
	mass=density*rp;
	mass=mass*rp;
}

void Point::SetE_point(double temp){
	E_point = temp;
}
double Point::GetE_point(int n){
	if (n==1) {return E_point;}
}

//settle atom type
/*void Point::SetType(int temp){
	type=temp;
}*/

//get the displacement of the point
double Point::GetD(int n){
	double d;
	if(n==1){d=x-x0;return d;}
	if(n==2){d=y-y0;return d;}
}

//get current coordinate of the point
//n =1 for x, 2 for y direction
double Point::GetX(int n){
	if(n==1){return x;}
	if(n==2){return y;}
}

//get originate coordinate of the point
//n =1 for x, 2 for y direction
double Point::GetX0(int n){
	if(n==1){return x0;}
	if(n==2){return y0;}
}

//get force of the point
//n =1 for x, 2 for y direction
double Point::GetF(int n){
	if(n==1){return fx;}
	if(n==2){return fy;}
}

//get acceleration of the point
//n =1 for x, 2 for y direction
double Point::Geta(int n){
	if(n==1){return ax;}
	if(n==2){return ay;}
}

//get velocity of the point
//n =1 for x, 2 for y direction
double Point::GetV(int n){
	if(n==1){return vx;}
	if(n==2){return vy;}
}


//get Body Force of the point             //@Shank
//n =1 for x, 2 for y direction
double Point::GetBF(int n){
	if(n==1){return Bx;}
	if(n==2){return By;}
}                                         //@Shank

//get stress of the point             //@Shank
//n =1 for x, 2 for y direction
double Point::GetSigma(int n){
	if(n==0){return Sxx;}
	if(n==1){return Syy;}
	if(n==2){return Sxy;}
}                                         //@Shank

void Point::SetSigma(int n,double temp){
	if(n==0){Sxx=temp;}
	if(n==1){Syy=temp;}
	if(n==2){Sxy=temp;}
}

double Point::GetStrain(int n){
	if(n==0){return Exx;}
	if(n==1){return Eyy;}
	if(n==2){return Exy;}
}                                         //@Shank

void Point::SetStrain(int n,double temp){
	if(n==0){Exx=temp;}
	if(n==1){Eyy=temp;}
	if(n==2){Exy=temp;}
}

double Point::GetStrain_c(int n){
	if(n==0){return Exx_c;}
	if(n==1){return Eyy_c;}
	if(n==2){return Exy_c;}
}                                         //@Shank

void Point::SetStrain_c(int n,double temp){
	if(n==0){Exx_c=temp;}
	if(n==1){Eyy_c=temp;}
	if(n==2){Exy_c=temp;}
}

void Point::SetK_INVmatrix(int n,double temp){
	if(n==0){K_shape_inv.matrix[0][0]=temp;}
	if(n==1){K_shape_inv.matrix[0][1]=temp;}
	if(n==2){K_shape_inv.matrix[1][0]=temp;}
	if(n==3){K_shape_inv.matrix[1][1]=temp;}
}

double Point::GetK_INVmatrix(int n){
	if(n==0){return K_shape_inv.matrix[0][0];}
	if(n==1){return K_shape_inv.matrix[0][1];}
	if(n==2){return K_shape_inv.matrix[1][0];}
	if(n==3){return K_shape_inv.matrix[1][1];}
}

void Point::Setomega(double temp){
	point_omega=temp;
}
double Point::Getomega(){
	return point_omega;
}