
//only static problem, 1D bar. nonlocal FEM
#ifndef HGlobal
#define HGlobal

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <iomanip>
#include <cmath>

#include "Spoint.h"
#include "Step.h"
#include "Smatrix.h"
#include "PD_forcestate.h"

//using namespace std;

class Global {
private:
int n_points;//number of nodes
int n_dof;//degree of freedom


double row;//density
double rp;//node's size, width or length of the node
double mass;//uniform mass of point
double delta;//uniform horizon size

int n_steps;//total number of steps
double dt;//delta t. time interval of the step
double BF_value;//body force value

int stepini, stepend;
int step_creep;    //--------@Shank
int stepinterval;//how many steps interval which you dump the result
double outputmagx;//enlarge the output displacement
double outputmagy;//enlarge the output displacement

vector<Step> GlobalSteps;//stall information of steps. For now try to only keep 3 steps
                
//store which displacement boundary condition are known
//first row stores the DOF
//second row stores the fixed "displacement"
Matrix BC_d_PD;//fixed end
Matrix BC_BF;  // for body force  //@Shank

//Matrix BC_v_PD;//loading

ForceState_Elastic PDforce;//paradynamic elastic potential

public:

double E;//Young's modulus. Should be uniform
double v;//poisson ratio
Matrix Creep_parameters; // for next time step

void ReadMaterial(char filename[100]);
void ReadPoints(char filename[100],vector<Point> &GlobalPoint);
void Read_BC_PD(char filename[100]);
void Read_BF(char filename[100]);       //@Shank
//void Read_LC_PD(char filename[100]);
void settle_BC_PD(vector <Point> &GlobalPoint,int Flag_creep_analysis);
void settle_BF(vector <Point> &GlobalPoint);    //@Shank
void settle_mass(vector <Point> &GlobalPoint, double E);


void Initialize_central_PD(char pointfile[100],char BCPDfile[100],char material_file[100]);
void Integration_central_PD(int step0,int step1,int step2);
void Solve();
void Shape_tensor(Point & centralpoint, vector <Point> &pointgroup);
void Damage_CreepStrain(Point & centralpoint, vector <Point> &pointgroup);

//void Initialize_creep(double dt);
//void find_force(Point & centralpoint, vector <Point> &pointgroup);
//void find_newKmatrix(Point & centralpoint, vector <Point> &pointgroup);
void Calculate_stress(int step2);
void Creep_evoluation(int step2);
void Add_force(int step2);
void Add_Damage(int step2);
void Add_Damage2(int step2);
void Add_Damage3(int step2);
void Add_Damage4(int step2);
void Add_Damage5(int step2);
void Crack_Propa(int step2);

//when initialize the Global, we have to tell the element, nodes and BC condition.
Global(char pointfile[100],char BCPDfile[100],char Material_file[100]){
	Initialize_central_PD(pointfile,BCPDfile,Material_file);


}

};

#endif
