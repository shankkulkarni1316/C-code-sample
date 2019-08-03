
//a0 is initially given
#include<iostream>
#include<stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>

#include "Step.h"

#define Maxline 100000 // in the bc file the max number of line is 100000. 


//write vtk file for PD
void Step::Write_vtk_PD(char filename[50],char option[1], double outputmagx, double outputmagy){
	FILE * vtkfile;
	vtkfile=fopen(filename,option);
	fprintf(vtkfile,"# vtk DataFile Version 1.0\n");
	fprintf(vtkfile,"PD result Serina copyright\n");
	fprintf(vtkfile,"ASCII\n\n");

	//nodes position
	fprintf(vtkfile,"DATASET POLYDATA\n");
	fprintf(vtkfile,"POINTS %d float\n",n_points);
    for (int i=0;i<n_points;i++){
    	double x,y,z;
    	//x=GlobalNode[i].GetX()+outputmagy*GlobalNode[i].GetDy();//final position
    	//y=GlobalNode[i].GetY()+outputmagx*GlobalNode[i].GetDx();//final position
    	x=GlobalPoint[i].GetX(1);//original coordinate
		y=GlobalPoint[i].GetX(2);
    	//z=100.00*GlobalPoint[i].GetD(1);;//now use z to represent dx
		z=0.0;

    	//fprintf(vtkfile,"%lf %lf %lf\n",x,outputmag*y,z);
    	fprintf(vtkfile,"%e %e %e\n",x,y,z);
		if (i==50){
			ofstream myfile;
			myfile.open("disp.txt", std::ios_base::app);
			myfile<<x<<"\t"<<y;
			myfile << endl << endl;
			myfile.close();
		}
    }
	
	//Vertices
	int temp=n_points;
    int n_numbers=2*temp;
    fprintf(vtkfile,"\nVERTICES %d %d\n",temp,n_numbers);
    for (int i=0;i<temp;i++){
    	fprintf(vtkfile,"1 %d\n",i);
    }

    //Point data
    fprintf(vtkfile,"\nPOINT_DATA %d\n",n_points);
    //points coordinate
	fprintf(vtkfile,"VECTORS c float\n");
    for (int i=0;i<n_points;i++){
		double x=GlobalPoint[i].GetX(1);
		double y=GlobalPoint[i].GetX(2);
		double z=0.0;
    	fprintf(vtkfile,"%e %e %e\n",x,y,z);
    }
	
	//points' displacement
	fprintf(vtkfile,"VECTORS d float\n");
    for (int i=0;i<n_points;i++){

		double dx=GlobalPoint[i].GetD(1);
		double dy=GlobalPoint[i].GetD(2);
		double dz=0.0;

    	fprintf(vtkfile,"%e %e %e\n",dx,dy,dz);
    }

	//points' mass
	fprintf(vtkfile,"VECTORS m float\n");
    for (int i=0;i<n_points;i++){

		double mx=GlobalPoint[i].GetMass();
		double my=GlobalPoint[i].GetMass();
		double mz=GlobalPoint[i].GetMass();

    	fprintf(vtkfile,"%e %e %e\n",mx,my,mz);
    }

    //nodal velocity
    fprintf(vtkfile,"\nVECTORS v float\n");
    for (int i=0;i<n_points;i++){
		double vx,vy,vz;
		vx=GlobalPoint[i].GetV(1);
		vy=GlobalPoint[i].GetV(2);
		vz=0.0;
    	fprintf(vtkfile,"%e %e %e\n",vx,vy,vz);
    }

    //nodal acceleration
    fprintf(vtkfile,"\nVECTORS a float\n");
    for (int i=0;i<n_points;i++){
		double ax,ay,az;
		ax=GlobalPoint[i].Geta(1);
		ay=GlobalPoint[i].Geta(2);
		az=0.0;
    	fprintf(vtkfile,"%e %e %e\n",ax,ay,az);
    }

    //nodal force
    fprintf(vtkfile,"\nVECTORS f float\n");
    for (int i=0;i<n_points;i++){
		double fx,fy,fz;
		fx=GlobalPoint[i].GetF(1);
		fy=GlobalPoint[i].GetF(2);
		fz=0.0;
    	fprintf(vtkfile,"%e %e %e\n",fx,fy,fz);
    }

    //nodal force
    fprintf(vtkfile,"\nVECTORS s float\n");
    for (int i=0;i<n_points;i++){
		double sx,sy,sz;
		sx=GlobalPoint[i].GetSigma(0);
		sy=GlobalPoint[i].GetSigma(1);
		sz=GlobalPoint[i].GetSigma(2);
    	fprintf(vtkfile,"%e %e %e\n",sx,sy,sz);
    }

    //nodal force
    fprintf(vtkfile,"\nVECTORS e float\n");
    for (int i=0;i<n_points;i++){
		double ex,ey,ez;
		ex=GlobalPoint[i].GetStrain(0);
		ey=GlobalPoint[i].GetStrain(1);
		ez=GlobalPoint[i].GetStrain(2);
    	fprintf(vtkfile,"%e %e %e\n",ex,ey,ez);
    }

   //points' mass, and broken coefficient
   fprintf(vtkfile,"\nVECTORS m&bro float\n");
   for (int i=0;i<n_points;i++){
    double mx=GlobalPoint[i].GetMass();
    double brok=GlobalPoint[i].broken_parameter;
   	fprintf(vtkfile,"%e %e %e\n",mx,brok,0.0);
   }

    fclose(vtkfile);
}


void Step::Write_vtk_Damage(char filename[50],char option[1], double outputmagx, double outputmagy){
	FILE * vtkfile;
	vtkfile=fopen(filename,option);
	fprintf(vtkfile,"# vtk DataFile Version 1.0\n");
	fprintf(vtkfile,"PD result Serina copyright\n");
	fprintf(vtkfile,"ASCII\n\n");

	//nodes position
	fprintf(vtkfile,"DATASET POLYDATA\n");
	fprintf(vtkfile,"POINTS %d float\n",n_points);
    for (int i=0;i<n_points;i++){
    	double x,y,z;
    	//x=GlobalNode[i].GetX()+outputmagy*GlobalNode[i].GetDy();//final position
    	//y=GlobalNode[i].GetY()+outputmagx*GlobalNode[i].GetDx();//final position
    	x=GlobalPoint[i].GetX(1);//original coordinate
		y=GlobalPoint[i].GetX(2);
    	//z=100.00*GlobalPoint[i].GetD(1);;//now use z to represent dx
		z=0.0;

    	//fprintf(vtkfile,"%lf %lf %lf\n",x,outputmag*y,z);
    	fprintf(vtkfile,"%e %e %e\n",x,y,z);

    }
	
	//Vertices
	int temp=n_points;
    int n_numbers=2*temp;
    fprintf(vtkfile,"\nVERTICES %d %d\n",temp,n_numbers);
    for (int i=0;i<temp;i++){
    	fprintf(vtkfile,"1 %d\n",i);
    }

    //Point data
    fprintf(vtkfile,"\nPOINT_DATA %d\n",n_points);
    //points coordinate
	fprintf(vtkfile,"VECTORS c float\n");
    for (int i=0;i<n_points;i++){
		double x=GlobalPoint[i].GetX(1);
		double y=GlobalPoint[i].GetX(2);
		double z=0.0;
    	fprintf(vtkfile,"%e %e %e\n",x,y,z);
    }
	
	//points' displacement
	fprintf(vtkfile,"VECTORS d float\n");
    for (int i=0;i<n_points;i++){

		double dx=GlobalPoint[i].GetD(1);
		double dy=GlobalPoint[i].GetD(2);
		double dz=0.0;

    	fprintf(vtkfile,"%e %e %e\n",dx,dy,dz);
    }



    //nodal damage
    fprintf(vtkfile,"\nVECTORS w float\n");
    for (int i=0;i<n_points;i++){
		double wx,wy,wz;
		wx=GlobalPoint[i].Getomega();
		wy = 0.0;
		wz = 0.0;
    	fprintf(vtkfile,"%e %e %e\n",wx,wy,wz);
    }

	//creep strain
    fprintf(vtkfile,"\nVECTORS ec float\n");
    for (int i=0;i<n_points;i++){
		double ex_c,ey_c,ez_c;
		ex_c=GlobalPoint[i].GetStrain_c(0);
		ey_c = GlobalPoint[i].GetStrain_c(1);
		ez_c = GlobalPoint[i].GetStrain_c(2);
    	fprintf(vtkfile,"%e %e %e\n",ex_c,ey_c,ez_c);
    }
	//stress
    fprintf(vtkfile,"\nVECTORS st float\n");
    for (int i=0;i<n_points;i++){
		double Sxx,Syy,Sxy;
		Sxx = GlobalPoint[i].GetSigma(0);
		Syy = GlobalPoint[i].GetSigma(1);
		Sxy = GlobalPoint[i].GetSigma(2);        
    	fprintf(vtkfile,"%e %e %e\n",Sxx,Syy,Sxy);
    }
	//body force
    fprintf(vtkfile,"\nVECTORS bf float\n");
    for (int i=0;i<n_points;i++){
		double fbx,fby,fbz;
		fbx=GlobalPoint[i].GetBF(1);
		fby=GlobalPoint[i].GetBF(2);
		fbz=0.0;
    	fprintf(vtkfile,"%e %e %e\n",fbx,fby,fbz);
		//printf("strain value %e %e %e \n",fbx,fby,fbz);
    }

	//damage youngs modulus
    fprintf(vtkfile,"\nVECTORS dy float\n");
    for (int i=0;i<n_points;i++){
		double dyx,dyy,dyz;
		dyx=GlobalPoint[i].GetE_point(1);
		dyy=0.0;
		dyz=0.0;
    	fprintf(vtkfile,"%e %e %e\n",dyx,dyy,dyz);
		//printf("strain value %e %e %e \n",fbx,fby,fbz);
    }

    fclose(vtkfile);
}

void Step::Write_Damage_forNextStep(char filename[50],char option[1], double outputmagx, double outputmagy){
	FILE * vtkfile;
	vtkfile=fopen(filename,option);
	//fprintf(vtkfile,"# vtk DataFile Version 1.0\n");
	//fprintf(vtkfile,"PD result Serina copyright\n");
	//fprintf(vtkfile,"ASCII\n\n");

	//nodes position
	//fprintf(vtkfile,"DATASET POLYDATA\n");
	//fprintf(vtkfile,"POINTS %d float\n",n_points);

    //nodal damage
    //fprintf(vtkfile,"\nVECTORS w float\n");
    for (int i=0;i<n_points;i++){
		double wx;
		wx=GlobalPoint[i].Getomega();
    	fprintf(vtkfile,"%e\n",wx);
    }

    fclose(vtkfile);
}

void Step::Write_Force_forNextStep(char filename[50],char option[1], double outputmagx, double outputmagy){
	FILE * vtkfile;
	vtkfile=fopen(filename,option);

    for (int i=0;i<n_points;i++){
		double wx;
		wx=GlobalPoint[i].GetBF(1);
    	fprintf(vtkfile,"%e\n",wx);
    }

    fclose(vtkfile);
}

//settle number of points
void Step:: settle_n_points(int tmpn){
    /*printf("n_atoms is %d\n",n_atoms);
	if (n_atoms!=0){
		printf("warning: the amount of atoms is not zero now, are you sure to resettle it again? Step.cpp (line: 179)\n");
	}*/
	n_points=tmpn;
	GlobalPoint.resize(n_points);//settle space for GlobalNode vector
}

//copy points
void Step::copy_Points(vector <Point> &tmppoint){
	//settle the size of the global node vector
	if (n_points==0){
		n_points=tmppoint.size();
		GlobalPoint.resize(n_points);
	}
	//check if there is a conflict about the total atom number
	if (n_points!=0&&n_points!=tmppoint.size()){
		printf("there is a conflict about the total atom number. Step.cpp (line 606)\n");
		exit(1);
	}

	for (int i=0;i<n_points;i++)
		{
			GlobalPoint[i]=tmppoint[i];
		}
	
}




