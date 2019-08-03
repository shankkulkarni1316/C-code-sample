
//left and right end of the panel is fixed
//now the BC_PD can represent the initial velocity
//considered the half mass in the top and bottom

#include<iostream>
#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

#include "Global.h"
#include "Spoint.h"
#define Maxline 100000 // in the bc file the max number of line is 100000. 


//Read material from material file
//It is now run before point is read
void Global::ReadMaterial(char filename[100]){
	printf("Now read material\n");
	FILE * materialfile;
	materialfile=fopen(filename,"r");
	char buf[Maxline];
	fgets(buf,Maxline,materialfile);//stip the first line

	//double E,v,row,rp;
	fscanf(materialfile,"E: %lf\n",&E);//young's modulus
	fscanf(materialfile,"v: %lf\n",&v);//poisson's ratio
	fscanf(materialfile,"row: %lf\n",&row);//young's modulus
	fscanf(materialfile,"rp: %lf\n",&rp);//point's size, it is uniform
	fscanf(materialfile,"delta: %lf\n",&delta);//uniform horizon size
	fclose (materialfile);

	mass=rp*rp*row;//update uniform mass of point
	printf("finished reading material coefficients\n");
}

void Global::settle_mass(vector<Point> &GlobalPoint, double E){
	if (mass==0){printf("error: point mass hasn't assigned\n");exit(1);}
	double halfmass=mass/2.0;
	double ymin=100.00;double ymax=-100.00;
	for (int i=0;i<n_points;i++){
		double y=GlobalPoint[i].GetX0(2);
		if (y>=ymax){
			ymax=y;
		}
		if (y<=ymin){
			ymin=y;
		}
	}
	
	for (int i=0;i<n_points;i++){
		double y=GlobalPoint[i].GetX0(2);
		if (y==ymin||y==ymax){
			GlobalPoint[i].SetMass(halfmass);
		}
		else{GlobalPoint[i].SetMass(mass);}
		GlobalPoint[i].SetE_point(E);
	}
}

//read information of points from a file and assign it to a point vector
//it should be used after material is read
void Global::ReadPoints(char filename[100],vector<Point> &GlobalPoint){
	printf("Now read points\n");
	
	FILE * pointfile;
	pointfile=fopen(filename,"r");
	char buf[Maxline];
	fscanf(pointfile, "Points Number: %d\n", &n_points);//total number of nodes
	
	n_dof=2*n_points;//2D problem
	
	GlobalPoint.resize(n_points);//settle the space of the vector
		
	fgets(buf,Maxline,pointfile);//stip the second line
	
	//int nnode;double x, y, dx, dy, fx, fy;//number of the node, coordinate of the node
	int npoint;double x,y;//number of the node, coordinate of the node
	
	for (int i=0;i<n_points;i++)
	{
		fscanf(pointfile,"%lf %lf\n",&x,&y);//nodes must put in sequence
		npoint=i+1;
		//printf("for point %d, x=%lf, y=%lf\n",npoint,x,y);
		
		Point tempoint(npoint);//set up a node whose number is nnode
		tempoint.SetX0(x,y);//assign initial coordinate of the current node.
		tempoint.SetX(x,y);//update the current coordinate
		//I use PD_forcestate to instore mass and rp. Just to save space
		//tempoint.Setdensity(row);//density
		//tempoint.Setrp(rp);//point's size, automatically update volumn
		//tempoint.SetMass();//update mass of point

		GlobalPoint[i]=tempoint;
		//printf("Now settled %d points\n",i+1);
	}


	fclose (pointfile);
	
	//update the neighbor list
	PDforce.FindNeighbor(GlobalPoint);
	
}

//read BC informatin from file
//assign BC_d and BC_v value which stores all the BC informatin for the whole running
//also get number of total steps.
void Global::Read_BC_PD(char filename[100]){
	printf("now read BC_PD\n");

	FILE * bcfile;
	bcfile=fopen(filename,"r");
	char buf[Maxline];
	fgets(buf,Maxline,bcfile);//stip the first line. It is a head line.

	//initial assign the disp
	int n_fixed_points;
	fscanf(bcfile, "steprange: %d %d\n", &stepini,&stepend);//initial and last step number
	fscanf(bcfile, "outputinterval: %d\n", &stepinterval);//how many steps between each dump
	//fscanf(bcfile, "creepStep: %d\n", &step_creep);//how many steps between each dump-------------------@Shank
	//fscanf(bcfile, "timestep: %lf\n", &BF_value);
	fscanf(bcfile, "timestep: %lf\n", &dt);
	fscanf(bcfile, "n_fixed_Point: %d\n", &n_fixed_points);//the num of dof that are fixed
	fgets(buf,Maxline,bcfile);//stip the sixth line. It is a head line.
	printf("there are %d fixed points\n", n_fixed_points);

	if (n_fixed_points>0){
		BC_d_PD.Setdim(4,n_fixed_points);//settle dimension of BC matrix
	}
	int fixed_point_number;double dx,vx;//position where to put the disp in the d vector
	int direction;//1 for x 2 for y
	for (int i=0;i<n_fixed_points;i++){
		fscanf(bcfile,"%d %d %lf %lf\n",&fixed_point_number,&direction,&dx,&vx);
	    BC_d_PD.matrix[0][i]=fixed_point_number;//node number which is fixed
	    BC_d_PD.matrix[1][i]=direction;
		BC_d_PD.matrix[2][i]=dx;//corresponding disp value
		BC_d_PD.matrix[3][i]=vx;//corresponding disp value
		printf("%d Point is fixed as %e disp\n",fixed_point_number,dx);
	}
	fscanf(bcfile, "creepStep: %d\n", &step_creep);
	fclose (bcfile);

	n_steps=stepend-stepini+1;//how many steps totally

	printf("finished reading BC_PD.\n\n");
}


void Global::settle_BC_PD(vector<Point> &GlobalPoint, int Flag_creep_analysis){
	printf("now settle boundary condition\n");
	//assign the fixed displacement
	int n_fixed_points=BC_d_PD.Getdim(2);//how many points are fixed for displacements
	for (int i=0;i<n_fixed_points;i++){
		//find the node and dof
		double ipointdouble;int ipoint;double value_d,value_v;int direction;double directiondouble;
		ipointdouble=BC_d_PD.matrix[0][i];//which point is fixed
		ipoint=ipointdouble;
		directiondouble=BC_d_PD.matrix[1][i];//1 for x and 2 for y
		direction=directiondouble;
		value_d=BC_d_PD.matrix[2][i];//fixed displacement value
		value_v=BC_d_PD.matrix[3][i];//value of initial velocity
		//debug
		//printf("ipoint: %d, direction: %d, dx: %e, vx: %e\n",ipoint,direction,value_d,value_v);

		//assign the flag and displacement value of the atom
		GlobalPoint[ipoint-1].SetFlagX(direction,0);//free the flag of freedom first
		GlobalPoint[ipoint-1].SetD(direction,value_d);//assign the disp value, also update the current position
		GlobalPoint[ipoint-1].SetV(direction,value_v);
	}
	printf("debug1 finished\n");

	//freeze the left end of the panel
	double xmin=100.00;
	double xmax=-100.00;
	for (int i=0;i<n_points;i++){
		double x=GlobalPoint[i].GetX0(1);
		if(x>xmax){xmax=x;}
		if (x<xmin){xmin=x;}
	}

	xmin=xmin+1.0;xmax=xmax-1.0;
	int direction=1;
//	for (int i=0;i<n_points;i++){
//		double xpoint=GlobalPoint[i].GetX0(1);
//		double ypoint=GlobalPoint[i].GetX0(2);
//		if ((xpoint==0.0) && (ypoint==20.0)){    //xpoint<=xmin||
//			GlobalPoint[i].SetFlagX(direction,1);//freeze the disp
//			GlobalPoint[i].SetFlagV(direction,1);//freeze the velocity
//			GlobalPoint[i].SetFlagA(direction,1);//freeze the acceleration
//			GlobalPoint[i].SetFlagF(direction,1);//freeze the force
//		}
//	}

	int direction2=2;
	if (Flag_creep_analysis == 0) {
		BF_value = -0.27777;
		for (int i=0;i<n_points;i++){
			double xpoint=GlobalPoint[i].GetX0(1);
			// setting initial damage to point as 0
			GlobalPoint[i].Setomega(0.0);
			if (xpoint<=xmin){    
			GlobalPoint[i].SetFlagX(direction,1);//freeze the disp
            GlobalPoint[i].SetFlagX(direction2,1);
			GlobalPoint[i].SetFlagV(direction,1);//freeze the velocity
            GlobalPoint[i].SetFlagV(direction2,1);
			GlobalPoint[i].SetFlagA(direction,1);//freeze the acceleration
            GlobalPoint[i].SetFlagA(direction2,1);
			GlobalPoint[i].SetFlagF(direction,1);//freeze the force
            GlobalPoint[i].SetFlagF(direction2,1);
			}
			if (xpoint>=xmax){    
				GlobalPoint[i].SetFlagBF(2,0);//freeze the disp
				GlobalPoint[i].SetBF(2,BF_value);
			}
		}
	}
	


	printf("finished settling boundary condition\n");
}


//read Body force informatin from file                                                        //@Shank
//assign BFx and BFy value which stores all the BF informatin for the whole running
void Global::Read_BF(char filename[100]){
	printf("now read Body force\n");

	FILE * bffile;
	bffile=fopen(filename,"r");
	char buf[Maxline];
	fgets(buf,Maxline,bffile);//stip the first line. It is a head line.

	int BF_value;

	fscanf(bffile, "body_force: %d\n", &BF_value);

	fclose (bffile);
	printf("Body force is %d.\n\n",&BF_value);
}


//initialize PD
void Global::Initialize_central_PD(char pointfile[100],char BCPDfile[100],char material_file[100]){
	printf("Now initialize PD part\n");
	
	//read input files
	//PDforce.ReadParameter(material_file);//read material into PDforce.
	ReadMaterial(material_file);
	PDforce.assignparameter(E,v,row,rp,delta);
	
	vector<Point> GlobalPoint;//this point group is used as basic model for all the rest steps of points to save time
	ReadPoints(pointfile,GlobalPoint);//read coordinate and update neighbor list
	settle_mass(GlobalPoint,E);//assigh mass to all the points
	
	Read_BC_PD(BCPDfile);
	int Flag_creep_analysis = 0;
	settle_BC_PD(GlobalPoint,Flag_creep_analysis);

	//Read_BF(BF_PDfile);   //@Shank

	//settle space for initial three steps
	GlobalSteps.clear();
	Step tmpstep;
	GlobalSteps.push_back(tmpstep);
	GlobalSteps.push_back(tmpstep);
	GlobalSteps.push_back(tmpstep);

	Step * Istepm1=&GlobalSteps[0];
	Step * Istep0=&GlobalSteps[1];
	Step * Istep1=&GlobalSteps[2];

	//settle step 0
	printf("Now build step0\n");
	Istep0->settle_istep(stepini);//assign the first step number
	Istep0->settle_dt(dt);//time interval
	Istep0->settle_time(0.0);//suppose time begins at 0.0.
	Istep0->settle_n_points(n_points);//number of points
	Istep0->copy_Points(GlobalPoint);//assign points. copy the initial position, mass, and neighbor list
	printf("I am here\n");
	PDforce.compute(Istep0->GlobalPoint,Istep0->GlobalPoint);//update force and acceleration for step 0. update force and acceleration once positions are settled
	printf("I am here too!\n");
	outputmagx=1.0;outputmagy=1.0;//this is the defalt value of output magnitude
	Istep0->Write_vtk_PD("dumpinitial.vtk","w",outputmagx, outputmagy);

	printf("Now build step -1\n");
	//build step -1
	Istepm1->settle_istep(stepini-1);
	Istepm1->settle_dt(dt);//time interval
	Istepm1->settle_time(-1.0*dt);//suppose time begins at 0.0.
	Istepm1->settle_n_points(n_points);//number of points
	Istepm1->copy_Points(GlobalPoint);//assign points

	//settle step -1 point position
	for (int ipoint=0;ipoint<n_points;ipoint++){
		double x0=Istep0->GlobalPoint[ipoint].GetX(1);
		double y0=Istep0->GlobalPoint[ipoint].GetX(2);
		double a0x=Istep0->GlobalPoint[ipoint].Geta(1);
		double a0y=Istep0->GlobalPoint[ipoint].Geta(2);
		double v0x=Istep0->GlobalPoint[ipoint].GetV(1);
		double v0y=Istep0->GlobalPoint[ipoint].GetV(2);
		//x-1=x0-dt*v0+1/2*dt^2*a0
		double xm1=x0-dt*v0x+1/2*dt*dt*a0x;
		double ym1=y0-dt*v0y+1/2*dt*dt*a0y;
		Istepm1->GlobalPoint[ipoint].SetX(xm1,ym1);

		if(ipoint==98){
		printf("point 99, x0=%lf, v0x=%lf, a0x=%lf, xm1=%lf, dt=%e\n",x0,v0x,a0x,xm1,dt);
		}
	}

	printf("Now build initial step 1\n");
	//build initial step
	Istep1->settle_istep(stepini+1);
	Istep1->settle_dt(dt);//time interval
	Istep1->settle_time(1*dt);//suppose time begins at 0.0.
	Istep1->settle_n_points(n_points);//number of points

	Istep1->copy_Points(Istep0->GlobalPoint);//copy points from last step

	//settle step 1 point position
	for (int ipoint=0;ipoint<n_points;ipoint++){
		double x0=Istep0->GlobalPoint[ipoint].GetX(1);
		double y0=Istep0->GlobalPoint[ipoint].GetX(2);
		double a0x=Istep0->GlobalPoint[ipoint].Geta(1);
		double a0y=Istep0->GlobalPoint[ipoint].Geta(2);
		double xm1=Istepm1->GlobalPoint[ipoint].GetX(1);
		double ym1=Istepm1->GlobalPoint[ipoint].GetX(2);
		//x1=2*x0-x-1+dt^2*a0;
		double x1=2*x0-xm1+dt*dt*a0x;
		double y1=2*y0-ym1+dt*dt*a0y;
		Istep1->GlobalPoint[ipoint].SetX(x1,y1);

		if (ipoint==98){
			printf("point 99, xm1=%e, x0=%e, x1=%e, a0=%e,dt=%e",xm1,x0,x1,a0x,dt);
			printf("\n");
		}

	}
	PDforce.compute(Istep1->GlobalPoint,Istep0->GlobalPoint);//update force and acceleration of step 1, once new positins are settled

	printf("finished building of step -1, 0, 1 of PD\n");
}

void Global::Integration_central_PD(int step0, int step1, int step2){
	printf("Now integrate PD step %d, %d and %d\n", step0, step1, step2);
	//find the content of the assigned step number
	Step * Step0, *Step1;
	Step * Step2;
	int count=0;int size=GlobalSteps.size();
	for (int i=0;i<size;i++){
		int nstep=GlobalSteps[i].GetN();
		if (nstep==step0){Step0=&GlobalSteps[i];count++;}
		if (nstep==step1){Step1=&GlobalSteps[i];count++;}
		if (nstep==step2){Step2=&GlobalSteps[i];count++;}
		//printf("now i is %d, nstep=%d, count=%d\n",i, nstep, count);
	}
	if (count<2){
		printf("these steps haven't be settle by the FEM part\n");
		exit(1);}

	//if the third step hasn't developed
	if (count==2){
		Step istep2(step2);
		GlobalSteps.push_back(istep2);
		Step0=&GlobalSteps[1];
		Step1=&GlobalSteps[2];
		Step2=&GlobalSteps[3];
	}

	double t1=Step1->Get_time();double t2=t1+dt;
	Step2->settle_dt(dt);
	Step2->settle_time(t2);
	Step2->settle_n_points(n_points);//number of points
	Step2->copy_Points(Step1->GlobalPoint);//copy points, use latest displacement

	//settle step 2 point position
	for (int ipoint=0;ipoint<n_points;ipoint++){
		double x1=Step1->GlobalPoint[ipoint].GetX(1);
		double y1=Step1->GlobalPoint[ipoint].GetX(2);
		double a1x=Step1->GlobalPoint[ipoint].Geta(1);
		double a1y=Step1->GlobalPoint[ipoint].Geta(2);
		double x0=Step0->GlobalPoint[ipoint].GetX(1);
		double y0=Step0->GlobalPoint[ipoint].GetX(2);
		double f_brok = Step1->GlobalPoint[ipoint].broken_parameter;

		//x2=2*x1-x0+dt^2*a1;

		double x2=2*x1-x0+dt*dt*a1x;
		double y2=2*y1-y0+dt*dt*a1y;

		Step2->GlobalPoint[ipoint].SetX(x2,y2);//points position of step2   

		//v1=(x2-x0)/(2*dt)
		//double v1x=(x1-x0)/(dt);                  //@Shank
		//double v1y=(y1-y0)/(dt);
		double v1x=(x2-x0)/(2*dt);                  //@Shank
		double v1y=(y2-y0)/(2*dt);

		if (f_brok == 1.0) {
			v1x = 0.0;
			v1y = 0.0;
		}
		Step1->GlobalPoint[ipoint].SetV(v1x,v1y);//velocity of step1


	}


	PDforce.compute(Step2->GlobalPoint,Step1->GlobalPoint);//update force and acceleration of step 2
	printf("	finished integration of step %d, %d and %d.\n",step0, step1, step2);
}

//
void Global::Solve(){
	int step0, step1, step2;
	//int step_creep=5;
	

	Creep_parameters.Setdim(1,n_points);

	for (int i=0;i+2<n_steps;i++){
		step0=i;step1=i+1;step2=i+2;
		Integration_central_PD(step0, step1, step2);

		//find position of step1 in the global vector of Steps
		int nstep,position;int count=0;int size=GlobalSteps.size();
		for (int j=0;j<size;j++){
			nstep=GlobalSteps[j].GetN();
			if (nstep==step1){position=j;count++;break;}
	    }
		
		//write dump files.
		//build filename
        double intpart,fracpart,temp,idouble;
        idouble=i;
        temp=idouble/stepinterval;
        fracpart=modf(temp,&intpart);
        if (fracpart==0){
        	//char filename_FEM[20];
			char steptime[20];
			char filename_PD[20];
        	sprintf(steptime,"%d",i);//changethe step time from integer to char

        	strcpy(filename_PD,"dumpPD");
        	strcat(filename_PD,steptime);
        	strcat(filename_PD,".vtk");
			
        	printf("write vtk file for step %d\n",GlobalSteps[position].GetN());
			GlobalSteps[position].Write_vtk_PD(filename_PD,"w",outputmagx,outputmagy);
        }
		//printf("fraction is %e\n",step_creep);

		if (i+3==step_creep) {
			Calculate_stress(step2);
			//GlobalSteps[position+1].Write_vtk_PD("stressAndStrain01.vtk","w",outputmagx,outputmagy);
			//Creep_evoluation(step2);
			//Add_force(step2);
            //Add_Damage(step2);
			GlobalSteps[position+1].Write_vtk_Damage("Damage01.vtk","w",outputmagx,outputmagy);
            //GlobalSteps[position+1].Write_Damage_forNextStep("Omega_01.txt","w",outputmagx,outputmagy);
            //GlobalSteps[position+1].Write_Force_forNextStep("Force_01.txt","w",outputmagx,outputmagy);
		}
		if (i+3==2*step_creep) {
			//Calculate_stress(step2);
			//GlobalSteps[position+1].Write_vtk_PD("stressAndStrain02.vtk","w",outputmagx,outputmagy);
			//Creep_evoluation(step2);
			//Add_force(step2);
            //Add_Damage2(step2);
			GlobalSteps[position+1].Write_vtk_Damage("Damage02.vtk","w",outputmagx,outputmagy);
		}
		if (i+3==3*step_creep) {
			//Calculate_stress(step2);
			//GlobalSteps[position+1].Write_vtk_PD("stressAndStrain03.vtk","w",outputmagx,outputmagy);
			//Creep_evoluation(step2);
			//Add_force(step2);
            //Add_Damage3(step2);
			GlobalSteps[position+1].Write_vtk_Damage("Damage03.vtk","w",outputmagx,outputmagy);
		}
		if (i+3==4*step_creep) {
			//Calculate_stress(step2);
			//GlobalSteps[position+1].Write_vtk_PD("stressAndStrain04.vtk","w",outputmagx,outputmagy);
			//Creep_evoluation(step2);
			//Add_force(step2);
            //Add_Damage4(step2);
			GlobalSteps[position+1].Write_vtk_Damage("Damage04.vtk","w",outputmagx,outputmagy);
		}
		if (i+3==5*step_creep) {
			//Calculate_stress(step2);
			//GlobalSteps[position+1].Write_vtk_PD("stressAndStrain05.vtk","w",outputmagx,outputmagy);
			//Creep_evoluation(step2);
			//Add_force(step2);
            //Add_Damage5(step2);
			GlobalSteps[position+1].Write_vtk_Damage("Damage05.vtk","w",outputmagx,outputmagy);
		}
		/*
		double intpart1,fracpart1,temp1,idouble1,stepdouble;
		idouble1=i;
		stepdouble = step_creep;
		temp1 = (i+3)/step_creep;
		fracpart1=modf(temp1,&intpart1);
		if (i+3==step_creep) {
			char steptime[20];
			char filename_PD[20];
        	sprintf(steptime,"%d",i);//changethe step time from integer to char

        	strcpy(filename_PD,"Damage");
        	strcat(filename_PD,steptime);
        	strcat(filename_PD,".vtk");
						
			//Creep_evoluation(step2);
			//Add_force(step2);
			//Crack_Propa(step2);
			//GlobalSteps[position+1].Write_vtk_Damage(filename_PD,"w",outputmagx,outputmagy);
		}
		*/
        //delete the very first step in the global steps vector to save some space
        //GlobalSteps.erase(GlobalSteps.begin());
		Step * Step0=new Step;
		Step * Step1=new Step;
		Step * Step2=new Step;
		*Step0=GlobalSteps[1];
		*Step1=GlobalSteps[2];
		*Step2=GlobalSteps[3];
		
		//destruct GlobalSteps
		{vector <Step> tmp1;
		tmp1.swap(GlobalSteps);}
		
		GlobalSteps.resize(3);
		GlobalSteps[0]=*Step0;
		GlobalSteps[1]=*Step1;
		GlobalSteps[2]=*Step2;
		
		delete Step0;
		delete Step1;
		delete Step2;
		
		
	}
	

}

void Global::Add_Damage(int step2){
	Step * Step2;
	int count=0;int size=GlobalSteps.size();
	for (int i=0;i<size;i++){
		int nstep=GlobalSteps[i].GetN();
		if (nstep==step2){Step2=&GlobalSteps[i];count++;}
	}
	vector <Point> &pointgroup = Step2->GlobalPoint;
	for (int i=0;i<n_points;i++){
		Point * Ipoint=&pointgroup[i];//current atom
		
		double Exx_c,Eyy_c,Exy_c;

		Point & centralpoint = *Ipoint;


		double x01,y01;
		x01=centralpoint.GetX0(1);y01=centralpoint.GetX0(2);
        double omega_present;
        omega_present = 0.75;
        if ((x01 == 0.0)&&(y01 == 0.0)) {
            centralpoint.Setomega(omega_present);
            
        }
        if ((x01 == 0.2)&&(y01 == 0.0)) {
            centralpoint.Setomega(omega_present);
            
        }
        if ((x01 == -0.2)&&(y01 == 0.0)) {
            centralpoint.Setomega(omega_present);
            
        }
        if ((x01 == 0.0)&&(y01 == 0.2)) {
            centralpoint.Setomega(omega_present);
            
        }
        if ((x01 == 0.2)&&(y01 == 0.2)) {
            centralpoint.Setomega(omega_present);
            
        }
        if ((x01 == -0.2)&&(y01 == 0.2)) {
            centralpoint.Setomega(omega_present);
            
        }        

		
	}
}



void Global::Crack_Propa(int step2){
	printf("now propagate crack.\n");
	Step * Step2;
	int count=0;int size=GlobalSteps.size();
	for (int i=0;i<size;i++){
		int nstep=GlobalSteps[i].GetN();
		if (nstep==step2){Step2=&GlobalSteps[i];count++;}
	}
	vector <Point> &pointgroup = Step2->GlobalPoint;
	PDforce.FindNeighbor_1(pointgroup, 1);
}

void Global::Add_force(int step2){
	Step * Step2;
	int count=0;int size=GlobalSteps.size();
	for (int i=0;i<size;i++){
		int nstep=GlobalSteps[i].GetN();
		if (nstep==step2){Step2=&GlobalSteps[i];count++;}
	}
	vector <Point> &pointgroup = Step2->GlobalPoint;
	for (int i=0;i<n_points;i++){
		Point * Ipoint=&pointgroup[i];//current atom
		
		double Exx_c,Eyy_c,Exy_c;

		Point & centralpoint = *Ipoint;

		Exx_c = centralpoint.GetStrain_c(0);
		Eyy_c = centralpoint.GetStrain_c(1);
		Exy_c = centralpoint.GetStrain_c(2);
		//omega = centralpoint.Getomega();

		double x01,y01,x02,y02,x1,y1,x2,y2;
		x01=centralpoint.GetX0(1);y01=centralpoint.GetX0(2);
		x1=centralpoint.GetX(1);y1=centralpoint.GetX(2);//deformed coordinate of the subject point

		//double E_temp = centralpoint.GetE_point(1);
		//printf(" E is %e\n",E_temp);

		int n_neighbors=centralpoint.neighbor_list.size();//how many neighbors it has		
		
		for (int j=0;j<n_neighbors;j++){
			int position=centralpoint.neighbor_list[j];//number of the current neighbor point
			//Point neighborpoint=pointgroup[position-1];
			Point * Ipoint1=&pointgroup[position-1];
			Point & neighborpoint = *Ipoint1;

			x02=neighborpoint.GetX0(1);y02=neighborpoint.GetX0(2);//coordinate of the neighbor point
			x2=neighborpoint.GetX(1);y2=neighborpoint.GetX(2);//deformed coordinate of the neighbor point
			double r0=sqrt((x02-x01)*(x02-x01)+(y02-y01)*(y02-y01));//bond length
			double r=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));// deformed bond length
			double ebar=r-r0;//bond enlongation
			//double w=1.0;
			double w=centralpoint.neighbor_w[j];
			
			double bfx,bfy;
			bfx = neighborpoint.GetBF(1);
			bfy = neighborpoint.GetBF(2);

			//volumn defraction
			double Rmin=delta-(rp/2.0);
			double Rmax=delta+(rp/2.0);
			double vf;
			if (r0<=Rmin){vf=1.0;}
			else if (r0<Rmax){vf=0.5+((delta-r0)/rp);}
			else {vf=0.0;}
			
			double K11,K12,K22;
			K11 = centralpoint.GetK_INVmatrix(0);
			K12 = centralpoint.GetK_INVmatrix(1);
			K22 = centralpoint.GetK_INVmatrix(3);

			double sigma_xx = w*(E/(1-(v*v)))*((Exx_c)+(v*(Eyy_c)));
			double sigma_yy = w*(E/(1-(v*v)))*((Eyy_c)+(v*(Exx_c)));
			double sigma_xy = w*((0.25*E)/(1+v))*((Exy_c)+(v*(Exy_c)));

			double bfx_n1, bfy_n1;
			double l_p = 0.1114;

			double bfx_n1_temp = ((((sigma_xx*K11)+(sigma_xy*K12))*(x02-x01)) + (((sigma_xx*K12)+(sigma_xy*K22))*(y02-y01)))*(vf)*(exp((-(r0*r0))/(l_p*l_p)));
			double bfy_n1_temp = ((((sigma_xy*K11)+(sigma_yy*K12))*(x02-x01)) + (((sigma_xy*K12)+(sigma_yy*K22))*(y02-y01)))*(vf)*(exp((-(r0*r0))/(l_p*l_p)));

			if(bfx_n1_temp>=3.0){
				bfx_n1_temp = 3.0;
			}

			bfx_n1 = bfx + bfx_n1_temp;
			bfy_n1 = bfy + bfy_n1_temp;


			neighborpoint.SetFlagBF(1,0);
			//bfx_n1 = 0.0;
			neighborpoint.SetBF(1,bfx_n1);
			//bfy_n1 = 0.0;
			neighborpoint.SetFlagBF(2,0);
			neighborpoint.SetBF(2,bfy_n1);
		}
		
	}
}


void Global::Calculate_stress(int step2){
	printf("Now calculating Stress.\n");
	Step * Step2;
	int count=0;int size=GlobalSteps.size();
	for (int i=0;i<size;i++){
		int nstep=GlobalSteps[i].GetN();
		if (nstep==step2){Step2=&GlobalSteps[i];count++;}
	}

	// Writing file for next time step
	ofstream nextTime1;
	nextTime1.open("nextTime1.txt");

	vector <Point> &pointgroup = Step2->GlobalPoint;

	//FILE * pointfile;
	//pointfile=fopen("Points.txt","r");
	//char buf1[Maxline];
	//char buf2[Maxline];

	//fgets(buf1,Maxline,pointfile);
	//fgets(buf2,Maxline,pointfile);

	//FILE * pointfile2;//-------------------- test
	//pointfile2=fopen("Points_FEM.txt","r");//-------------------- test
	//char buf3[Maxline];//-------------------- test
	//char buf4[Maxline];//-------------------- test

	//fgets(buf3,Maxline,pointfile2);//-------------------- test
	//fgets(buf4,Maxline,pointfile2);//-------------------- test

	//for (int i=0;i<n_points;i++){//-------------------- test
		//Point * Ipoint=&pointgroup[i];//current atom//-------------------- test

		//double q,r;//-------------------- test
		//fscanf(pointfile2,"%lf %lf\n",&q,&r);//-------------------- test
	
		//Point & centralpoint = *Ipoint;//-------------------- test
	
		//centralpoint.SetX(q,r);//-------------------- test
		
	//}
	//fclose (pointfile2); //-------------------- test
	for (int i=0;i<n_points;i++){
		Point * Ipoint=&pointgroup[i];//current atom
		//double x,y;

		//fscanf(pointfile,"%lf %lf\n",&x,&y);
		Point & centralpoint = *Ipoint;
		//centralpoint.SetX0(x,y);

		Shape_tensor(*Ipoint,pointgroup);
		
	}

	//fclose (pointfile);

	
	nextTime1.close();

	// Write text file for next time increment
	double omega;
	ofstream nextTime;
	nextTime.open("nextTime.txt");

	// write points position for next time
	ofstream points_nextTime;
	points_nextTime.open("points_nextTime.txt");
	points_nextTime<<"Points Number: "<<n_points<<"\n";
	points_nextTime<<"#x y \n";

	for (int i=0;i<n_points;i++){
		Point * Ipoint=&pointgroup[i];//current atom
		Point & centralpoint = *Ipoint;
		double x01,y01,x1,y1;
		x01=centralpoint.GetX0(1);y01=centralpoint.GetX0(2);
		x1=centralpoint.GetX(1);y1=centralpoint.GetX(2);

		double Exx, Eyy, Exy, Sxx, Syy, Sxy;
		Exx = centralpoint.GetStrain(0);
		Eyy = centralpoint.GetStrain(1);
		Exy = centralpoint.GetStrain(2);

		Sxx = centralpoint.GetSigma(0);
		Syy = centralpoint.GetSigma(1);
		Sxy = centralpoint.GetSigma(2);

		//double omega;

		omega = centralpoint.Getomega();

		nextTime<<x1<<"\t"<<y1<<"\t"<<Exx<<"\t"<<Eyy<<"\t"<<Exy<<"\t"<<Sxx<<"\t"<<Syy<<"\t"<<Sxy<<"\t"<<omega<<"\n";
		points_nextTime<<x1<<"\t"<<y1<<"\n";
	}
	nextTime.close();
	points_nextTime.close();
}

void Global::Shape_tensor(Point & centralpoint, vector <Point> &pointgroup){
	double x01,y01,x02,y02,x1,y1,x2,y2;
	//x01=centralpoint.GetX0(1);y01=centralpoint.GetX0(2);//coordinate of the subject point
	// reading co odinates of original points as stress must be calculated according to original positions

	x01=centralpoint.GetX0(1);y01=centralpoint.GetX0(2);
	x1=centralpoint.GetX(1);y1=centralpoint.GetX(2);//deformed coordinate of the subject point

	int n_neighbors=centralpoint.neighbor_list.size();//how many neighbors it has
	
	
	double omega;
	omega = centralpoint.Getomega();
	// initialize all to zero
	centralpoint.K_shape.matrix[0][0] = 0.0;
	centralpoint.K_shape.matrix[0][1] = 0.0;
	centralpoint.K_shape.matrix[1][0] = 0.0;
	centralpoint.K_shape.matrix[1][1] = 0.0;
	centralpoint.F_dg_1.matrix[0][0] = 0.0;
	centralpoint.F_dg_1.matrix[0][1] = 0.0;
	centralpoint.F_dg_1.matrix[1][0] = 0.0;
	centralpoint.F_dg_1.matrix[1][1] = 0.0;

	for (int i=0;i<n_neighbors;i++){
		int position=centralpoint.neighbor_list[i];//number of the current neighbor point
		Point neighborpoint=pointgroup[position-1];
		x02=neighborpoint.GetX0(1);y02=neighborpoint.GetX0(2);//coordinate of the neighbor point
		x2=neighborpoint.GetX(1);y2=neighborpoint.GetX(2);//deformed coordinate of the neighbor point
		double r0=sqrt((x02-x01)*(x02-x01)+(y02-y01)*(y02-y01));//bond length
		double r=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));// deformed bond length
		double ebar=r-r0;//bond enlongation
		//double w=1.0;

		//volumn defraction
		double Rmin=delta-(rp/2.0);
		double Rmax=delta+(rp/2.0);
		double vf;
		if (r0<=Rmin){vf=1.0;}
		else if (r0<Rmax){vf=0.5+((delta-r0)/rp);}
		else {vf=0.0;}
		double w=centralpoint.neighbor_w[i];

		// Calculation of Shape tensor K                                              @Shank
		centralpoint.K_shape.matrix[0][0] = centralpoint.K_shape.matrix[0][0] + w*(x02-x01)*(x02-x01)*vf;
		centralpoint.K_shape.matrix[0][1] = centralpoint.K_shape.matrix[0][1] + w*(x02-x01)*(y02-y01)*vf;
		centralpoint.K_shape.matrix[1][0] = centralpoint.K_shape.matrix[1][0] + w*(y02-y01)*(x02-x01)*vf;
		centralpoint.K_shape.matrix[1][1] = centralpoint.K_shape.matrix[1][1] + w*(y02-y01)*(y02-y01)*vf;

		// Calculation of deformation gradient F                                              @Shank
		centralpoint.F_dg_1.matrix[0][0] = centralpoint.F_dg_1.matrix[0][0] + w*(x2-x1)*(x02-x01)*vf;
		centralpoint.F_dg_1.matrix[0][1] = centralpoint.F_dg_1.matrix[0][1] + w*(x2-x1)*(y02-y01)*vf;
		centralpoint.F_dg_1.matrix[1][0] = centralpoint.F_dg_1.matrix[1][0] + w*(y2-y1)*(x02-x01)*vf;
		centralpoint.F_dg_1.matrix[1][1] = centralpoint.F_dg_1.matrix[1][1] + w*(y2-y1)*(y02-y01)*vf;

    //if (x01 == 0.0 && y01 == 0.3){
		//printf("negibor x02 and y02 %f\t%f\t%f\n",x02,y02,vf);
		//printf("negibor x2 and y2 %f\t%f\t%f\n",x2,y2,vf);
		//printf("negibor y02 %f\n",y02);
		//printf("K12 %f\n",centralpoint.K_shape.matrix[0][1]);
	//}
			

	}



	double det_K = (centralpoint.K_shape.matrix[0][0]*centralpoint.K_shape.matrix[1][1])-(centralpoint.K_shape.matrix[0][1]*centralpoint.K_shape.matrix[1][0]);
	if (det_K == 0.0) {
		centralpoint.K_shape_inv.matrix[0][0] = 0.0;
		centralpoint.K_shape_inv.matrix[1][1] = 0.0;
		centralpoint.K_shape_inv.matrix[0][1] = 0.0;
		centralpoint.K_shape_inv.matrix[1][0] = 0.0;
	}
	else {
		centralpoint.K_shape_inv.matrix[0][0] = (1/det_K)*(centralpoint.K_shape.matrix[1][1]);
		centralpoint.K_shape_inv.matrix[1][1] = (1/det_K)*(centralpoint.K_shape.matrix[0][0]);
		centralpoint.K_shape_inv.matrix[0][1] = (-1.0)*(1/det_K)*(centralpoint.K_shape.matrix[1][0]);
		centralpoint.K_shape_inv.matrix[1][0] = (-1.0)*(1/det_K)*(centralpoint.K_shape.matrix[0][1]);
	}


	centralpoint.F_dg = centralpoint.F_dg_1*centralpoint.K_shape_inv;

	if (det_K == 0.0) {
		centralpoint.point_strain.matrix[0][0] = ((0.5)*((centralpoint.F_dg.matrix[0][0])+(centralpoint.F_dg.matrix[0][0])));
		centralpoint.point_strain.matrix[0][1] = ((0.5)*((centralpoint.F_dg.matrix[1][0])+(centralpoint.F_dg.matrix[0][1])));
		centralpoint.point_strain.matrix[1][0] = ((0.5)*((centralpoint.F_dg.matrix[0][1])+(centralpoint.F_dg.matrix[1][0])));
		centralpoint.point_strain.matrix[1][1] = ((0.5)*((centralpoint.F_dg.matrix[1][1])+(centralpoint.F_dg.matrix[1][1])));
	}
	else {
		centralpoint.point_strain.matrix[0][0] = ((0.5)*((centralpoint.F_dg.matrix[0][0])+(centralpoint.F_dg.matrix[0][0]))) - 1.0;
		centralpoint.point_strain.matrix[0][1] = ((0.5)*((centralpoint.F_dg.matrix[1][0])+(centralpoint.F_dg.matrix[0][1])));
		centralpoint.point_strain.matrix[1][0] = ((0.5)*((centralpoint.F_dg.matrix[0][1])+(centralpoint.F_dg.matrix[1][0])));
		centralpoint.point_strain.matrix[1][1] = ((0.5)*((centralpoint.F_dg.matrix[1][1])+(centralpoint.F_dg.matrix[1][1]))) - 1.0;
	}


	double sigma_xx = (1-omega)*(E/(1-(v*v)))*((centralpoint.point_strain.matrix[0][0])+(v*(centralpoint.point_strain.matrix[1][1])));
	double sigma_yy = (1-omega)*(E/(1-(v*v)))*((centralpoint.point_strain.matrix[1][1])+(v*(centralpoint.point_strain.matrix[0][0])));
	double sigma_xy = (1-omega)*((0.25*E)/(1+v))*((centralpoint.point_strain.matrix[0][1])+(v*(centralpoint.point_strain.matrix[1][0])));
	double Exx, Eyy, Exy;

	Exx = centralpoint.point_strain.matrix[0][0];
	Eyy = centralpoint.point_strain.matrix[1][1];
	Exy = centralpoint.point_strain.matrix[0][1];

	//centralpoint.point_stress.push_back(sigma_xx);
	//centralpoint.point_stress.push_back(sigma_yy);
	//centralpoint.point_stress.push_back(sigma_xy);

	centralpoint.SetStrain(0,Exx);
	centralpoint.SetStrain(1,Eyy);
	centralpoint.SetStrain(2,Exy);

	// HARD CODE to 240,0,0
	//sigma_xx = 240.0;
	//sigma_yy = 0.0;
	//sigma_xy = 0.0;
	//sigma_xx = sigma_xx/3;
	//sigma_yy = sigma_yy/3;
	//sigma_xy = sigma_xy/3;
	centralpoint.SetSigma(0,sigma_xx);
	centralpoint.SetSigma(1,sigma_yy);
	centralpoint.SetSigma(2,sigma_xy);

}

void Global::Creep_evoluation(int step2){
	printf("Now evoluting creep\n");

	Step * Step2;
	int count=0;int size=GlobalSteps.size();
	for (int i=0;i<size;i++){
		int nstep=GlobalSteps[i].GetN();
		if (nstep==step2){Step2=&GlobalSteps[i];count++;}
	}
	vector <Point> &pointgroup = Step2->GlobalPoint;
	for (int i=0;i<n_points;i++){
		Point * Ipoint=&pointgroup[i];//current atom
		Point & centralpoint = *Ipoint;
		Damage_CreepStrain(*Ipoint,pointgroup);
	}

}


void Global::Damage_CreepStrain(Point & centralpoint, vector <Point> &pointgroup){
	double alpha_creep = 0.47845;   //  
	double C_creep = 1.47*1e-29;
	double n_creep = 10.147;
	double D_creep = 2.73*1e-30;
	double p_creep = 10.949;
	double q_creep = 6.35;
	double omega_past, omega_present;
	double dt_c = 1.0;

	double Sxx, Syy, Sxy, Exx_c,Eyy_c,Exy_c, omega;
	Sxx = centralpoint.GetSigma(0);
	Syy = centralpoint.GetSigma(1);
	Sxy = centralpoint.GetSigma(2);
	Exx_c = centralpoint.GetStrain_c(0);
	Eyy_c = centralpoint.GetStrain_c(1);
	Exy_c = centralpoint.GetStrain_c(2);
	omega = centralpoint.Getomega();

	double S1,S2,von_mises,sigma_1,Sigma_r,d_omega,d_creep_Exx,d_creep_Eyy,d_creep_Exy;
	S1 = (0.5*(Sxx+Syy))+(sqrt((0.25*(Sxx-Syy)*(Sxx-Syy))+(Sxy*Sxy)));
	S2 = (0.5*(Sxx+Syy))-(sqrt((0.25*(Sxx-Syy)*(Sxx-Syy))+(Sxy*Sxy)));
	von_mises = sqrt((S1*S1)+(S2*S2)-(S1*S2));
	sigma_1 = S1; 
	Sigma_r= (alpha_creep*sigma_1) + ((1-alpha_creep)*von_mises);
	d_omega = dt_c*D_creep*(pow(Sigma_r,p_creep))*(exp(q_creep*omega))*(((1.0)-(exp(-q_creep)))/(q_creep));
	//d_omega = dt_c*D_creep*(pow(Sxx,p_creep))*(exp(q_creep*omega))*(((1.0)-(exp(-q_creep)))/(q_creep));
	if (von_mises == 0.0) {
		d_creep_Exx = 0.0;
		d_creep_Eyy = 0.0;
		d_creep_Exy = 0.0;
	}
	else {
		d_creep_Exx=dt_c*(1.5)*(C_creep)*(pow(von_mises,n_creep))*((((Sxx*0.5))-(Syy*0.5))/(von_mises))*(exp(((2*(n_creep+1))/(3.141592654*(sqrt(1+(3/n_creep)))))*(((sigma_1*sigma_1)/(von_mises*von_mises)))*(pow(omega,1.5))));
		d_creep_Eyy=dt_c*(1.5)*(C_creep)*(pow(von_mises,n_creep))*((((Syy*0.5))-(Sxx*0.5))/(von_mises))*(exp(((2*(n_creep+1))/(3.141592654*(sqrt(1+(3/n_creep)))))*(((sigma_1*sigma_1)/(von_mises*von_mises)))*(pow(omega,1.5))));
		d_creep_Exy=dt_c*(1.5)*(C_creep)*(pow(von_mises,n_creep))*((Sxy)/(von_mises))*(exp(((2*(n_creep+1))/(3.141592654*(sqrt(1+(3/n_creep)))))*(((sigma_1*sigma_1)/(von_mises*von_mises)))*(pow(omega,1.5))));
		//d_creep_Exx=dt_c*(1.5)*(C_creep)*(pow(von_mises,n_creep))*((((Sxx*0.5))-(Syy*0.5))/(von_mises))*(exp(((2*(n_creep+1))/(3.141592654*(sqrt(1+(3/n_creep)))))*(((sigma_1*sigma_1)/(von_mises*von_mises)))*(pow(omega,1.5))));
		//d_creep_Eyy=dt_c*(1.5)*(C_creep)*(pow(von_mises,n_creep))*((((Syy*0.5))-(Sxx*0.5))/(von_mises))*(exp(((2*(n_creep+1))/(3.141592654*(sqrt(1+(3/n_creep)))))*(((sigma_1*sigma_1)/(von_mises*von_mises)))*(pow(omega,1.5))));
		//d_creep_Exy=dt_c*(1.5)*(C_creep)*(pow(von_mises,n_creep))*((Sxy)/(von_mises))*(exp(((2*(n_creep+1))/(3.141592654*(sqrt(1+(3/n_creep)))))*(((sigma_1*sigma_1)/(von_mises*von_mises)))*(pow(omega,1.5))));
	}

	//if (d_omega>=5.0){
		
	//	printf("point number and value %e\n",d_omega);
	//	d_omega = 0.75;
	//}
	//else if (d_omega>=0.75){
	//	d_omega = 0.5;
		//printf("point number and value %e\n",d_omega);
	//}
	omega_present = omega + d_omega;
	Exx_c = Exx_c + d_creep_Exx;
	Eyy_c = Eyy_c + d_creep_Eyy;
	Exy_c = Exy_c + d_creep_Exy;

	centralpoint.Setomega(omega_present);
	//if (omega_present>=0.5){
	//	//d_omega = 0.5;
	//	printf("point number and value of new omega %e\n",omega_present);
	//}	
	centralpoint.SetStrain_c(0,Exx_c);
	centralpoint.SetStrain_c(1,Eyy_c);
	centralpoint.SetStrain_c(2,Exy_c);

	//update damage youngs modulus
	//double E_dam = centralpoint.GetE_point(1);
	double E_dam = (1-omega_present)*(148000.0);
	//centralpoint.SetE_point(E_dam);
	
}

/*

void Global::Initialize_creep(double dt_c){
	printf("Initializing creep.\n");

	// Initilizing core code
	ReadMaterial("Material.txt");
	PDforce.assignparameter(E,v,row,rp,delta);
	vector<Point> GlobalPoint;
	ReadPoints("points_nextTime.txt",GlobalPoint);
	settle_mass(GlobalPoint,E);//assigh mass to all the points

	FILE * creepfile;
	creepfile=fopen("nextTime.txt","r");

	//read disp, strain, stress and omega
	Creep_parameters.Setdim(9,n_points);//settle dimension of creep matrix

	double x1, y1, Exx, Eyy, Exy, Sxx, Syy, Sxy, omega;
	double S1,S2,von_mises,sigma_1,Sigma_r,d_omega,d_creep_Exx,d_creep_Eyy,d_creep_Exy;
	double alpha_creep = 0.47845;   //  
	double C_creep = 1.47*1e-29;
	double n_creep = 10.147;
	double D_creep = 2.73*1e-30;
	double p_creep = 10.949;
	double q_creep = 6.35;
	double omega_past, omega_present;

	for (int i=0;i<n_points;i++){
		fscanf(creepfile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&x1,&y1,&Exx,&Eyy,&Exy,&Sxx,&Syy,&Sxy,&omega);
	    //Creep_parameters.matrix[0][i]=x1;
	    //Creep_parameters.matrix[1][i]=y1;
		Creep_parameters.matrix[2][i]=Exx;
		Creep_parameters.matrix[3][i]=Eyy;
		Creep_parameters.matrix[4][i]=Exy;
		Creep_parameters.matrix[5][i]=Sxx;
		Creep_parameters.matrix[6][i]=Syy;
		Creep_parameters.matrix[7][i]=Sxy;
		Creep_parameters.matrix[8][i]=omega;

		//calculating increments using Liu-Murakami model
		S1 = (0.5*(Sxx+Syy))+(sqrt((0.25*(Sxx-Syy)*(Sxx-Syy))+(Sxy*Sxy)));
		S2 = (0.5*(Sxx+Syy))-(sqrt((0.25*(Sxx-Syy)*(Sxx-Syy))+(Sxy*Sxy)));
		von_mises = sqrt((S1*S1)+(S2*S2)-(S1*S2));
		sigma_1 = S1; 
		Sigma_r= (alpha_creep*sigma_1) + ((1-alpha_creep)*von_mises);
		d_omega = dt_c*D_creep*(pow(Sigma_r,p_creep))*(exp(q_creep*omega))*(((1.0)-(exp(-q_creep)))/(q_creep));

		if (von_mises == 0.0) {
			d_creep_Exx = 0.0;
			d_creep_Eyy = 0.0;
			d_creep_Exy = 0.0;
		}
		else {
			d_creep_Exx=dt_c*(1.5)*(C_creep)*(pow(von_mises,n_creep))*((((Sxx*0.5))-(Syy*0.5))/(von_mises))*(exp(((2*(n_creep+1))/(3.141592654*(sqrt(1+(3/n_creep)))))*(((sigma_1*sigma_1)/(von_mises*von_mises)))*(pow(omega,1.5))));
			d_creep_Eyy=dt_c*(1.5)*(C_creep)*(pow(von_mises,n_creep))*((((Syy*0.5))-(Sxx*0.5))/(von_mises))*(exp(((2*(n_creep+1))/(3.141592654*(sqrt(1+(3/n_creep)))))*(((sigma_1*sigma_1)/(von_mises*von_mises)))*(pow(omega,1.5))));;
			d_creep_Exy=dt_c*(1.5)*(C_creep)*(pow(von_mises,n_creep))*((Sxy)/(von_mises))*(exp(((2*(n_creep+1))/(3.141592654*(sqrt(1+(3/n_creep)))))*(((sigma_1*sigma_1)/(von_mises*von_mises)))*(pow(omega,1.5))));;
		}

		// Updating the damage omega
		//omega_past = GlobalPoint[i].Getomega();
		omega_present = omega + d_omega;
		GlobalPoint[i].Setomega(omega_present);

		// Find increment in stress using increment in creep strain
		double d_sigma_xx = (1-omega_present)*(E/(1-(v*v)))*((d_creep_Exx)+(v*(d_creep_Eyy)));
		double d_sigma_yy = (1-omega_present)*(E/(1-(v*v)))*((d_creep_Eyy)+(v*(d_creep_Exx)));
		double d_sigma_xy = (1-omega_present)*((0.25*E)/(1+v))*((d_creep_Exy)+(v*(d_creep_Exy)));

		//Sxx = Sxx + d_sigma_xx;
		//Syy = Syy + d_sigma_yy;
		//Sxy = Sxy + d_sigma_xy;

		//Exx = Exx + d_creep_Exx;
		//Eyy = Eyy + d_creep_Eyy;
		//Exy = Exy + d_creep_Exy;

		GlobalPoint[i].point_stress.push_back(d_sigma_xx);
		GlobalPoint[i].point_stress.push_back(d_sigma_yy);
		GlobalPoint[i].point_stress.push_back(d_sigma_xy);

		GlobalPoint[i].SetSigma(0,Sxx);
		GlobalPoint[i].SetSigma(1,Syy);
		GlobalPoint[i].SetSigma(2,Sxy);
		//GlobalPoint[i].SetSigma(0,Sxx);
		//GlobalPoint[i].SetSigma(1,Syy);
		//GlobalPoint[i].SetSigma(2,Sxy);

		//GlobalPoint[i].SetStrain(0,Exx);
		//GlobalPoint[i].SetStrain(1,Eyy);
		//GlobalPoint[i].SetStrain(2,Exy);

		//printf("stress values are: %lf \t %lf \t %lf \n",d_sigma_xx, d_sigma_yy, d_sigma_xy);
		if (i==555){
			printf("x1 %f\n",x1);
			printf("y1 %f\n",y1);
			printf("Sxx %f\n",Sxx);
			printf("Syy %f\n",Syy);
			printf("Sxy %f\n",Sxy);
			printf("Exx %f\n",Exx);
			printf("Eyy %f\n",Eyy);
			printf("Exy %f\n",Exy);
			printf("d_sigma_xx %f\n",d_sigma_xx);
			printf("d_sigma_yy %f\n",d_sigma_yy);
			printf("d_sigma_xy %f\n",d_sigma_xy);
			printf("d_creep_Exx %f\n",d_creep_Exx);
			printf("d_creep_Eyy %f\n",d_creep_Eyy);
			printf("d_creep_Exy %f\n",d_creep_Exy);
		}

	}
	fclose (creepfile);
	

	Read_BC_PD("BC_PD.txt");
	int Flag_creep_analysis = 1;
	settle_BC_PD(GlobalPoint,Flag_creep_analysis);

	printf("Calculate and apply body force due to change in creep strain.\n");
	vector <Point> &pointgroup = GlobalPoint;
	ofstream bodyforce;
	bodyforce.open("bodyforce.txt");

	for (int i=0;i<n_points;i++){
		Point * Ipoint=&pointgroup[i];
		find_newKmatrix(*Ipoint,pointgroup);
	}
	for (int i=0;i<n_points;i++){
		Point * Ipoint=&pointgroup[i];
		find_force(*Ipoint,pointgroup);
		Point & centralpoint=*Ipoint;
		double force_x, force_y;
		double tem_x, tem_y;
		tem_x = centralpoint.GetX0(1);
		tem_y = centralpoint.GetX0(2);
		force_x = centralpoint.GetBF(1);
		force_y = centralpoint.GetBF(2);
		bodyforce<<tem_x<<"\t"<<tem_y<<"\t"<<force_x<<"\t"<<force_y<<"\n";
		if (i==555){
			printf("body force values are: %lf \t %lf \n",force_x, force_y);

			printf("co cordinates are: %lf \t %lf \n",tem_x, tem_y);
		}
	}
	bodyforce.close();

	//settle space for initial three steps
	GlobalSteps.clear();
	Step tmpstep;
	GlobalSteps.push_back(tmpstep);
	GlobalSteps.push_back(tmpstep);
	GlobalSteps.push_back(tmpstep);

	Step * Istepm1=&GlobalSteps[0];
	Step * Istep0=&GlobalSteps[1];
	Step * Istep1=&GlobalSteps[2];

	//settle step 0
	printf("Now build step0\n");
	Istep0->settle_istep(stepini);//assign the first step number
	Istep0->settle_dt(dt);//time interval
	Istep0->settle_time(0.0);//suppose time begins at 0.0.
	Istep0->settle_n_points(n_points);//number of points
	Istep0->copy_Points(GlobalPoint);//assign points. copy the initial position, mass, and neighbor list

	PDforce.compute(Istep0->GlobalPoint,Istep0->GlobalPoint);//update force and acceleration for step 0. update force and acceleration once positions are settled

	outputmagx=1.0;outputmagy=1.0;//this is the defalt value of output magnitude
	Istep0->Write_vtk_PD("dumpinitial.vtk","w",outputmagx, outputmagy);

	printf("Now build step -1\n");
	//build step -1
	Istepm1->settle_istep(stepini-1);
	Istepm1->settle_dt(dt);//time interval
	Istepm1->settle_time(-1.0*dt);//suppose time begins at 0.0.
	Istepm1->settle_n_points(n_points);//number of points
	Istepm1->copy_Points(GlobalPoint);//assign points

	//settle step -1 point position
	for (int ipoint=0;ipoint<n_points;ipoint++){
		double x0=Istep0->GlobalPoint[ipoint].GetX(1);
		double y0=Istep0->GlobalPoint[ipoint].GetX(2);
		double a0x=Istep0->GlobalPoint[ipoint].Geta(1);
		double a0y=Istep0->GlobalPoint[ipoint].Geta(2);
		double v0x=Istep0->GlobalPoint[ipoint].GetV(1);
		double v0y=Istep0->GlobalPoint[ipoint].GetV(2);
		//x-1=x0-dt*v0+1/2*dt^2*a0
		double xm1=x0-dt*v0x+1/2*dt*dt*a0x;
		double ym1=y0-dt*v0y+1/2*dt*dt*a0y;
		Istepm1->GlobalPoint[ipoint].SetX(xm1,ym1);

		if(ipoint==98){
		printf("point 99, x0=%lf, v0x=%lf, a0x=%lf, xm1=%lf, dt=%e\n",x0,v0x,a0x,xm1,dt);
		}
	}

	printf("Now build initial step 1\n");
	//build initial step
	Istep1->settle_istep(stepini+1);
	Istep1->settle_dt(dt);//time interval
	Istep1->settle_time(1*dt);//suppose time begins at 0.0.
	Istep1->settle_n_points(n_points);//number of points

	Istep1->copy_Points(Istep0->GlobalPoint);//copy points from last step

	//settle step 1 point position
	for (int ipoint=0;ipoint<n_points;ipoint++){
		double x0=Istep0->GlobalPoint[ipoint].GetX(1);
		double y0=Istep0->GlobalPoint[ipoint].GetX(2);
		double a0x=Istep0->GlobalPoint[ipoint].Geta(1);
		double a0y=Istep0->GlobalPoint[ipoint].Geta(2);
		double xm1=Istepm1->GlobalPoint[ipoint].GetX(1);
		double ym1=Istepm1->GlobalPoint[ipoint].GetX(2);
		//x1=2*x0-x-1+dt^2*a0;
		double x1=2*x0-xm1+dt*dt*a0x;
		double y1=2*y0-ym1+dt*dt*a0y;
		Istep1->GlobalPoint[ipoint].SetX(x1,y1);

		if (ipoint==98){
			printf("point 99, xm1=%e, x0=%e, x1=%e, a0=%e,dt=%e",xm1,x0,x1,a0x,dt);
			printf("\n");
		}

	}
	PDforce.compute(Istep1->GlobalPoint,Istep0->GlobalPoint);//update force and acceleration of step 1, once new positins are settled

	printf("finished building of step -1, 0, 1 of PD\n");


}

void Global::find_force(Point & centralpoint,vector <Point> &pointgroup){
	
	double x01,y01,x02,y02,x1,y1,x2,y2;
	double w=1.0;//doesn't consider fracture
	int n_neighbors=centralpoint.neighbor_list.size();
	x01=centralpoint.GetX0(1);y01=centralpoint.GetX0(2);

	/*
	centralpoint.K_shape.matrix[0][0] = 0.0;
	centralpoint.K_shape.matrix[0][1] = 0.0;
	centralpoint.K_shape.matrix[1][0] = 0.0;
	centralpoint.K_shape.matrix[1][1] = 0.0;

	for (int i=0;i<n_neighbors;i++){
		int position=centralpoint.neighbor_list[i];//number of the current neighbor point
		Point neighborpoint=pointgroup[position-1];
		x02=neighborpoint.GetX0(1);y02=neighborpoint.GetX0(2);//coordinate of the neighbor point
		//x2=neighborpoint.GetX(1);y2=neighborpoint.GetX(2);//deformed coordinate of the neighbor point
		double r0=sqrt((x02-x01)*(x02-x01)+(y02-y01)*(y02-y01));//bond length
		//double r=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));// deformed bond length
		//double ebar=r-r0;//bond enlongation

		//volumn defraction
		double Rmin=delta-rp/2.0;
		double Rmax=delta+rp/2.0;
		double vf;
		if (r0<=Rmin){vf=1.0;}
		else if (r0<Rmax){vf=1.0/2.0+(delta-r0)/rp;}
		else {vf=0.0;}

		// Calculation of Shape tensor K                                              @Shank
		centralpoint.K_shape.matrix[0][0] = centralpoint.K_shape.matrix[0][0] + w*(x02-x01)*(x02-x01)*vf;
		centralpoint.K_shape.matrix[0][1] = centralpoint.K_shape.matrix[0][1] + w*(x02-x01)*(y02-y01)*vf;
		centralpoint.K_shape.matrix[1][0] = centralpoint.K_shape.matrix[1][0] + w*(y02-y01)*(x02-x01)*vf;
		centralpoint.K_shape.matrix[1][1] = centralpoint.K_shape.matrix[1][1] + w*(y02-y01)*(y02-y01)*vf;
	}

	double det_K = (centralpoint.K_shape.matrix[0][0]*centralpoint.K_shape.matrix[1][1])-(centralpoint.K_shape.matrix[0][1]*centralpoint.K_shape.matrix[1][0]);
	centralpoint.K_shape_inv.matrix[0][0] = (1/det_K)*(centralpoint.K_shape.matrix[1][1]);
	centralpoint.K_shape_inv.matrix[1][1] = (1/det_K)*(centralpoint.K_shape.matrix[0][0]);
	centralpoint.K_shape_inv.matrix[0][1] = (-1.0)*(1/det_K)*(centralpoint.K_shape.matrix[1][0]);
	centralpoint.K_shape_inv.matrix[1][0] = (-1.0)*(1/det_K)*(centralpoint.K_shape.matrix[0][1]);
	

	double K11,K12,K21,K22;
	K11 = centralpoint.GetK_INVmatrix(0);
	K12 = centralpoint.GetK_INVmatrix(1);
	K21 = centralpoint.GetK_INVmatrix(2);
	K22 = centralpoint.GetK_INVmatrix(3);
	
	double Sxx,Syy,Sxy;
	Sxx = centralpoint.GetSigma(0);
	Syy = centralpoint.GetSigma(1);
	Sxy = centralpoint.GetSigma(2);

	double Fx,Fy;
	Fx=0.0;
	Fy=0.0;
	for (int i=0;i<n_neighbors;i++){
		int position=centralpoint.neighbor_list[i];
		Point neighborpoint=pointgroup[position-1];

		x02=neighborpoint.GetX0(1);y02=neighborpoint.GetX0(2);
		double r0=sqrt((x02-x01)*(x02-x01)+(y02-y01)*(y02-y01));//bond length
		//volumn defraction
		double Rmin=delta-rp/2.0;
		double Rmax=delta+rp/2.0;
		double vf;
		if (r0<=Rmin){vf=1.0;}
		else if (r0<Rmax){vf=1.0/2.0+(delta-r0)/rp;}
		else {vf=0.0;}

		double Kn11,Kn12,Kn21,Kn22;
		Kn11 = neighborpoint.GetK_INVmatrix(0);
		Kn12 = neighborpoint.GetK_INVmatrix(1);
		Kn21 = neighborpoint.GetK_INVmatrix(2);
		Kn22 = neighborpoint.GetK_INVmatrix(3);
		double Snxx,Snyy,Snxy;
		Snxx = neighborpoint.GetSigma(0);
		Snyy = neighborpoint.GetSigma(1);
		Snxy = neighborpoint.GetSigma(2);

		Fx = Fx + ((w*((((Sxx*(K11))+(Sxy*(K12)))*(x02-x01))+(((Sxx*(K12))+(Sxy*(K22)))*(y02-y01))))-((w*((((Snxx*(Kn11))+(Snxy*(Kn12)))*(x01-x02))+(((Snxx*(Kn12))+(Snxy*(Kn22)))*(y01-y02))))))*(vf);
		Fy = Fy + ((w*((((Sxy*(K11))+(Syy*(K12)))*(x02-x01))+(((Sxy*(K12))+(Syy*(K22)))*(y02-y01))))-((w*((((Snxy*(Kn11))+(Snyy*(Kn12)))*(x01-x02))+(((Snxy*(Kn12))+(Snyy*(Kn22)))*(y01-y02))))))*(vf);
	}
	//printf("stress, Sxx=%lf,Syy=%lf, Sxy=%lf\n",Sxx,Syy,Sxy);
	//printf("force, Fx=%lf, Fy=%lf\n",Fx,Fy);
	centralpoint.SetFlagBF(1,0);
	centralpoint.SetBF(1,Fx);
	centralpoint.SetFlagBF(2,0);
	centralpoint.SetBF(2,Fy);
}


void Global::find_newKmatrix(Point & centralpoint,vector <Point> &pointgroup){

	double x01,y01,x02,y02,x1,y1,x2,y2;
	double w=1.0;//doesn't consider fracture
	int n_neighbors=centralpoint.neighbor_list.size();
	x01=centralpoint.GetX0(1);y01=centralpoint.GetX0(2);

	centralpoint.K_shape.matrix[0][0] = 0.0;
	centralpoint.K_shape.matrix[0][1] = 0.0;
	centralpoint.K_shape.matrix[1][0] = 0.0;
	centralpoint.K_shape.matrix[1][1] = 0.0;

	for (int i=0;i<n_neighbors;i++){
		int position=centralpoint.neighbor_list[i];//number of the current neighbor point
		Point neighborpoint=pointgroup[position-1];
		x02=neighborpoint.GetX0(1);y02=neighborpoint.GetX0(2);//coordinate of the neighbor point
		//x2=neighborpoint.GetX(1);y2=neighborpoint.GetX(2);//deformed coordinate of the neighbor point
		double r0=sqrt((x02-x01)*(x02-x01)+(y02-y01)*(y02-y01));//bond length
		//double r=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));// deformed bond length
		//double ebar=r-r0;//bond enlongation

		//volumn defraction
		double Rmin=delta-rp/2.0;
		double Rmax=delta+rp/2.0;
		double vf;
		if (r0<=Rmin){vf=1.0;}
		else if (r0<Rmax){vf=1.0/2.0+(delta-r0)/rp;}
		else {vf=0.0;}

		// Calculation of Shape tensor K                                              @Shank
		centralpoint.K_shape.matrix[0][0] = centralpoint.K_shape.matrix[0][0] + w*(x02-x01)*(x02-x01)*vf;
		centralpoint.K_shape.matrix[0][1] = centralpoint.K_shape.matrix[0][1] + w*(x02-x01)*(y02-y01)*vf;
		centralpoint.K_shape.matrix[1][0] = centralpoint.K_shape.matrix[1][0] + w*(y02-y01)*(x02-x01)*vf;
		centralpoint.K_shape.matrix[1][1] = centralpoint.K_shape.matrix[1][1] + w*(y02-y01)*(y02-y01)*vf;
	}

	double det_K = (centralpoint.K_shape.matrix[0][0]*centralpoint.K_shape.matrix[1][1])-(centralpoint.K_shape.matrix[0][1]*centralpoint.K_shape.matrix[1][0]);
	centralpoint.K_shape_inv.matrix[0][0] = (1/det_K)*(centralpoint.K_shape.matrix[1][1]);
	centralpoint.K_shape_inv.matrix[1][1] = (1/det_K)*(centralpoint.K_shape.matrix[0][0]);
	centralpoint.K_shape_inv.matrix[0][1] = (-1.0)*(1/det_K)*(centralpoint.K_shape.matrix[1][0]);
	centralpoint.K_shape_inv.matrix[1][0] = (-1.0)*(1/det_K)*(centralpoint.K_shape.matrix[0][1]);

	double K11,K12,K21,K22;
	K11 = centralpoint.K_shape_inv.matrix[0][0];
	K12 = centralpoint.K_shape_inv.matrix[0][1];
	K21 = centralpoint.K_shape_inv.matrix[1][0];
	K22 = centralpoint.K_shape_inv.matrix[1][1];

	centralpoint.SetK_INVmatrix(0,K11);
	centralpoint.SetK_INVmatrix(1,K12);
	centralpoint.SetK_INVmatrix(2,K21);
	centralpoint.SetK_INVmatrix(3,K22);
}
*/