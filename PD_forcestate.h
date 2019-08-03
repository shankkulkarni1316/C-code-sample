
//2D paradynamic force state, plane stress
//based on Q.V.Le etl, 'A two-dimensional ordinary, state-based peridynamic model for linearly elastic solids', Int, J. Number. Meth, 2014
//updated the function to calculate the force state

#ifndef HForceState_Elastic
#define HForceState_Elastic

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>

#include "Smatrix.h"
#include "Spoint.h"
#include "Step.h"

using namespace std;

class ForceState_Elastic {
private:
double E;//Young's moudulus
double k;//bulk modulus
double miu;//shear modulus
double v;//depth of the potential well
double kp;//kprime
double delta;//horizon size
double rp;//uniform node size
double row;//density

public:
ForceState_Elastic(){
	E=0.0;
	k=0.0;
	miu=0.0;
	v=0.0;
	kp=0.0;
	rp=0.0;
	delta=0.0;
	row=0.0;
}

double Get_E(){return E;}
double Get_k(){return k;}
double Get_miu(){return miu;}
double Get_kp(){return kp;}
double Get_delta(){return delta;}
double Get_rp(){return rp;}


void compute(vector <Point> &pointgroup,vector <Point> &pointgroup1);
void ReadParameter(char parameterfile[100]);
void assignparameter(double,double,double,double,double);
void FindNeighbor(vector <Point> &pointgroup);
void FindNeighbor_1(vector <Point> &pointgroup, int star);
void ForceState_compute(Point &centralpoint, vector <Point> &pointgroup);
void PointValueCal(Point & centralpoint, vector <Point> &pointgroup);
void ForceStateCal(Point &pointx, Point &pointxp,double &forcestatex,double &forcestatey, double w);
double cal_w(double r0,double r,double x2,double y2,Point &centralpoint,int nei_position,Point &neighborpoint);
};
#endif
