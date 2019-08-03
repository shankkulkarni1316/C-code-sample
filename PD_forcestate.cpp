

#include "PD_forcestate.h"

#define Maxline 100000 // in the bc file the max number of line is 100000.

//given position of a group of points, accumulate all their force based on elastic
//also settle all the acceleration of points
//neighbor list should already be settled before
void ForceState_Elastic::compute(vector <Point> &pointgroup,vector <Point> &pointgroup1){
	int n_points=pointgroup.size();//how many atoms are there
	
	//clear all the residual atom force 
	for (int ipoint=0;ipoint<n_points;ipoint++){
		pointgroup[ipoint].SetF(0.0,0.0);
	}
	//printf("I am here!");
	//calculate all the values for forcestate calculation
	for (int i=0;i<n_points;i++){
		Point * Ipoint=&pointgroup[i];//current atom
		//printf("I am here too!");
		PointValueCal(*Ipoint,pointgroup);
		//printf("I am here too3!");
	}
	//printf("I am here too4!");
	//accumulate point force. 
	for (int ipoint=0;ipoint<n_points;ipoint++){
		Point * Ipoint=&pointgroup[ipoint];//current atom
		ForceState_compute(*Ipoint,pointgroup);
	}
	
	//calculate the acceleration of all points
	for (int ipoint=0;ipoint<n_points;ipoint++){
		Point * Ipoint=&pointgroup[ipoint];
		Point * Ipoint1=&pointgroup1[ipoint];
		double fx=Ipoint->GetF(1);
		double fy=Ipoint->GetF(2);

		double Bx=Ipoint->GetBF(1);   //@Shank
		double By=Ipoint->GetBF(2);

		double f_brok = Ipoint->broken_parameter;

		//double vx=Ipoint->GetV(1);
		double vx=Ipoint1->GetV(1);
		double vy=Ipoint1->GetV(2);
		double m=row*rp*rp;
		double ax=(fx+Bx-10.0*vx)/m;          //@Shank
		double ay=(fy+By-10.0*vy)/m;

		if (f_brok == 1.0) {
			ax = 0.0;
			ay = 0.0;
		}

		Ipoint->Seta(ax,ay);
		//printf("velocity %f\n",vx);
	}
}

//read parameter for LJ
void ForceState_Elastic::ReadParameter(char parameterfile[100]){
	//read parameters from file
	printf("now read harmonic parameter\n");
	FILE * parafile;
	parafile=fopen(parameterfile,"r");
	char buf[Maxline];
	fgets(buf,Maxline,parafile);//skip first lines
	fscanf(parafile,"E: %lf\n",&E);
	fscanf(parafile,"v: %lf\n",&v);
	fscanf(parafile,"row: %lf\n",&row);
	fscanf(parafile,"rp: %lf\n",&rp);//node size
	fscanf(parafile,"delta: %lf\n",&delta);//horizon size
	fclose (parafile);

	//derive the rest parameters
	//k=E/(3(1-2v))  bulk modulus
	double temp=2.0*v;
	temp=1.0-temp;
	temp=3.0*temp;
	k=E/temp;

	//G=E/(2(1+v)) shear modulus
	temp=1.0+v;
	temp=2.0*temp;
	miu=E/temp;

	//kp=k+miu/9*(v+1)^2/(2v-1)^2
	double temp1,temp2,temp3;
	temp1=2.0*v;temp1=temp1-1.0;temp1=temp1*temp1;
	temp2=v+1.0;temp2=temp2*temp2;
	temp3=miu/9.0;
	temp3=temp3*temp2;
    temp3=temp3/temp1;
	kp=k+temp3;

	printf("Now finished reading harmonic parameter\n");
}

//assign parameter by given values
void ForceState_Elastic::assignparameter(double tmpE,double tmpv,double tmprow,double tmprp,double tmpdelta){
	E=tmpE;
	v=tmpv;
	row=tmprow;
	rp=tmprp;
	delta=tmpdelta;

	//derive the rest parameters
	//k=E/(3(1-2v))  bulk modulus
	double temp=2.0*v;
	temp=1.0-temp;
	temp=3.0*temp;
	k=E/temp;

	//G=E/(2(1+v)) shear modulus
	temp=1.0+v;
	temp=2.0*temp;
	miu=E/temp;

	//kp=k+miu/9*(v+1)^2/(2v-1)^2
	double temp1,temp2,temp3;
	temp1=2.0*v;temp1=temp1-1.0;temp1=temp1*temp1;
	temp2=v+1.0;temp2=temp2*temp2;
	temp3=miu/9.0;
	temp3=temp3*temp2;
    temp3=temp3/temp1;
	kp=k+temp3;
}

//find neighbor list for each point
//small deformation supposed.
//the neighbor is found based on initial positions
void ForceState_Elastic::FindNeighbor(vector <Point> &pointgroup){
	printf("now find neighbor\n");
	int n_points=pointgroup.size();//how many atoms are there
	for (int ipoint=0;ipoint<n_points;ipoint++){
		Point * pointi=&pointgroup[ipoint];
		//double xi=atomi->GetX();
		double xi=pointi->GetX0(1);
		double yi=pointi->GetX0(2);
		//only search atom larger
		for (int jpoint=ipoint+1;jpoint<n_points;jpoint++){
			Point * pointj=& pointgroup[jpoint];
			//double xj=atomj->GetX();
			double xj=pointj->GetX0(1);
			double yj=pointj->GetX0(2);
			//double yj=atomj->GetY();
			double rij=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj));
			double Rmax=delta+rp/2.0;
			if (rij<= Rmax){

					pointi->neighbor_list.push_back(pointj->GetN());
					pointj->neighbor_list.push_back(pointi->GetN());
					pointi->neighbor_w.push_back(1.0);
					pointj->neighbor_w.push_back(1.0);
				
			

			}
		}
	}
}

void ForceState_Elastic::FindNeighbor_1(vector <Point> &pointgroup, int star){
	printf("now find neighbor\n");
	int n_points=pointgroup.size();//how many atoms are there
	for (int ipoint=0;ipoint<n_points;ipoint++){
		Point * pointi=&pointgroup[ipoint];
		pointi->neighbor_list.clear();
	}
	for (int ipoint=0;ipoint<n_points;ipoint++){
		Point * pointi=&pointgroup[ipoint];
		//double xi=atomi->GetX();
		double xi=pointi->GetX0(1);
		double yi=pointi->GetX0(2);
		//only search atom larger
		for (int jpoint=ipoint+1;jpoint<n_points;jpoint++){
			Point * pointj=& pointgroup[jpoint];
			//double xj=atomj->GetX();
			double xj=pointj->GetX0(1);
			double yj=pointj->GetX0(2);
			//double yj=atomj->GetY();
			double rij=sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj));
			double Rmax=delta+rp/2.0;
			if (rij<= Rmax){
				if ((xi*xj<=0.0) && (yi<0.8||yj<0.8)) {
					continue;
				}
				else {
					pointi->neighbor_list.push_back(pointj->GetN());
					pointj->neighbor_list.push_back(pointi->GetN());
				}
			
				//pointi->neighbor_w.push_back(1.0);
				//pointj->neighbor_w.push_back(1.0);
			}
		}
	}
}

//calculate theta, ed(*)x=B and q value for the current point
void ForceState_Elastic::PointValueCal(Point & centralpoint, vector <Point> &pointgroup){
	double x01,y01,x02,y02,x1,y1,x2,y2;
	x01=centralpoint.GetX0(1);y01=centralpoint.GetX0(2);//coordinate of the subject point
	x1=centralpoint.GetX(1);y1=centralpoint.GetX(2);//deformed coordinate of the subject point

	double E_temp = centralpoint.GetE_point(1);     //--------------------------------Shank _ variable E
	//derive the rest parameters
	//k=E/(3(1-2v))  bulk modulus
	double temp=2.0*v;
	temp=1.0-temp;
	temp=3.0*temp;
	k=E_temp/temp;

	//G=E/(2(1+v)) shear modulus
	temp=1.0+v;
	temp=2.0*temp;
	miu=E_temp/temp;

	//kp=k+miu/9*(v+1)^2/(2v-1)^2
	double temp1,temp2,temp3;
	temp1=2.0*v;temp1=temp1-1.0;temp1=temp1*temp1;
	temp2=v+1.0;temp2=temp2*temp2;
	temp3=miu/9.0;
	temp3=temp3*temp2;
    temp3=temp3/temp1;
	kp=k+temp3;


	int n_neighbors=centralpoint.neighbor_list.size();//how many neighbors it has
	vector <double> X; vector <double> Y; //store the deformed coordinate of all bond vectors
	vector <double> R0; //bond length of all the bond vectors
	vector <double> R; //deformed bond length
	vector <double> ee; //bond enlongation of all the bond vectors
	vector <double> Vf;//volumn fraction of neighbor points
	vector <double> Ed;//delatation
	double w=1.0;//doesn't consider fracture
	for (int i=0;i<n_neighbors;i++){
		int position=centralpoint.neighbor_list[i];//number of the current neighbor point
		Point neighborpoint=pointgroup[position-1];
		x02=neighborpoint.GetX0(1);y02=neighborpoint.GetX0(2);//coordinate of the neighbor point
		x2=neighborpoint.GetX(1);y2=neighborpoint.GetX(2);//deformed coordinate of the neighbor point
		double r0=sqrt((x02-x01)*(x02-x01)+(y02-y01)*(y02-y01));//bond length
		double r=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));// deformed bond length
		double ebar=r-r0;//bond enlongation
		double x,y;
		x=x2-x1;y=y2-y1;

		//volumn defraction
		double Rmin=delta-rp/2.0;
		double Rmax=delta+rp/2.0;
		double vf;
		if (r0<=Rmin){vf=1.0;}
		else if (r0<Rmax){vf=1.0/2.0+(delta-r0)/rp;}      // changed r to r0       
		else {vf=0.0;}

		R0.push_back(r0);//bond length
		R.push_back(r);
		X.push_back(x);//deformed bond vector
		Y.push_back(y);//deformed bond vector
		ee.push_back(ebar);//bond enlongation 
		Vf.push_back(vf);//volumn fraction

	}
	
	for (int i=0;i<n_neighbors;i++){
		int position=centralpoint.neighbor_list[i];//number of the current neighbor point
		Point neighborpoint=pointgroup[position-1];
		x02=neighborpoint.GetX0(1);y02=neighborpoint.GetX0(2);//coordinate of the neighbor point
		x2=neighborpoint.GetX(1);y2=neighborpoint.GetX(2);//deformed coordinate of the neighbor point
		double r0=sqrt((x02-x01)*(x02-x01)+(y02-y01)*(y02-y01));//bond length
		double r=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));// deformed bond length
		double ebar=r-r0;//bond enlongation

		//volumn defraction
		double Rmin=delta-rp/2.0;
		double Rmax=delta+rp/2.0;
		double vf;
		if (r0<=Rmin){vf=1.0;}
		else if (r0<Rmax){vf=1.0/2.0+(delta-r0)/rp;}         // changed r to r0
		else {vf=0.0;}
		//calculate w 
        cal_w(r0,r,x02,y02,centralpoint,i,neighborpoint);

	}
	//printf("I am here Shank22!\n");
	//calculate p and q, two integration
	//q=wx(*)x, p=wx(*)e
	//don't consider fracture, w always equals 1.0
	double vol=rp*rp;//volumn of each point. uniform point size
	double p,q;
	//double p;//updated by Serina Dec14/2015
	p=0.0;q=0.0;
	for (int i=0;i<n_neighbors;i++){
		double r0=R0[i];
		double vf=Vf[i];
		double ebar=ee[i];
		//double w=centralpoint.neighbor_w[i];
		q=q+w*r0*r0*vf*vol;
		p=p+w*r0*ebar*vf*vol;
	}
	
	double alpha=8.0*miu/q;
	//alpha=8.0*miu/q;

	//theta=2(2v-1)/(v-1)*p/q
	double scalar;
	double theta;
	scalar=2.0*(2.0*v-1.0)/(v-1.0);
	theta=p/q;theta=scalar*theta;
	//printf("I am here Shank22!\n");

	//calculate ed=e-theta*x/3
	for (int i=0;i<n_neighbors;i++){
		double r0=R0[i];
		double ebar=ee[i];
		double ed=ebar-theta*r0/3.0;
		Ed.push_back(ed);
	}
	

	//printf("I am here Shank22!\n");
	//B=w*ed(*)x=p-theta/3*q
	double B=0.0;

	//B=0.0;
	for (int i=0;i<n_neighbors;i++){
		//double w=1.0;//doesn't consider fracture
		double ed=Ed[i];
		double r0=R0[i];
		double vf=Vf[i];
		//double w=centralpoint.neighbor_w[i];

		B=B+w*ed*r0*vf*vol;

	}
	
	//assign to the point value. it will be later used for cal the force state
	centralpoint.theta=theta;
	centralpoint.B=B;
	centralpoint.alpha=alpha;
	centralpoint.q=q;


	//if (n_neighbors == 0){
	//	int tiktok = centralpoint.GetN();
	//	printf("values are: %d \t %d\n",tiktok,n_neighbors);
	//}
}

//calculate the forcestate. it acts on point x, it comes from point xp
void ForceState_Elastic::ForceStateCal(Point &pointx, Point &pointxp,double &forcestatex,double &forcestatey, double w){
	double x01,y01,x02,y02,x1,y1,x2,y2;
	x01=pointx.GetX0(1);y01=pointx.GetX0(2);//coordinate of the subject point
	x1=pointx.GetX(1);y1=pointx.GetX(2);//deformed coordinate of the subject point
	
	x02=pointxp.GetX0(1);y02=pointxp.GetX0(2);//coordinate of the neighbor point
	x2=pointxp.GetX(1);y2=pointxp.GetX(2);//deformed coordinate of the neighbor point

	double E_temp = pointx.GetE_point(1);     //--------------------------------Shank _ variable E
	//derive the rest parameters
	//k=E/(3(1-2v))  bulk modulus
	double temp=2.0*v;
	temp=1.0-temp;
	temp=3.0*temp;
	k=E_temp/temp;

	//G=E/(2(1+v)) shear modulus
	temp=1.0+v;
	temp=2.0*temp;
	miu=E_temp/temp;

	//kp=k+miu/9*(v+1)^2/(2v-1)^2
	double temp1,temp2,temp3;
	temp1=2.0*v;temp1=temp1-1.0;temp1=temp1*temp1;
	temp2=v+1.0;temp2=temp2*temp2;
	temp3=miu/9.0;
	temp3=temp3*temp2;
    temp3=temp3/temp1;
	kp=k+temp3;

	double r0=sqrt((x02-x01)*(x02-x01)+(y02-y01)*(y02-y01));//bond length
	double r=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));// deformed bond length
	double ebar=r-r0;//bond enlongation
	double x,y;
	x=x2-x1;y=y2-y1;
	//double w=1.0;
    
	//double omega_check = pointx.Getomega();
	//if (omega_check >= 0.99) {
	//	w = 0.0;
	//}
	
	//t=2(2v-1)/(v-1)*(k'*theta-alpha/3*B)*wx/q+alpha*w*ed
	double scalar1;
	double theta=pointx.theta;
	double alpha=pointx.alpha;
	double q=pointx.q;
	double B=pointx.B;
	double ed=ebar-theta*r0/3.0;
	
	scalar1=2.0*(2.0*v-1.0)/(v-1.0)*(kp*theta-alpha/3.0*B)/q;
	//if (scalar1<1.0e-10){printf("\n");exit(1);}

	double forcescalar;
	forcescalar=scalar1*w*r0+alpha*w*ed;

	//double fx,fy;//deformed bond vector
	forcestatex=forcescalar*x/r;
	forcestatey=forcescalar*y/r;
}

//calculate the force state of a central point
void ForceState_Elastic::ForceState_compute(Point &centralpoint, vector <Point> &pointgroup){
	double fx,fy;//point force density of central point
	fx=0.0;fy=0.0;
	
	double Rmin=delta-rp/2.0;
	double Rmax=delta+rp/2.0;
	
	double x01,y01,x02,y02;
	x01=centralpoint.GetX0(1);y01=centralpoint.GetX0(2);//coordinate of the subject point

	int n_broken_bond=0;

	int n_neighbors=centralpoint.neighbor_list.size();//how many neighbors it has
	for (int i=0;i<n_neighbors;i++){
		int position=centralpoint.neighbor_list[i];//number of the current neighbor point
		Point neighborpoint=pointgroup[position-1];
		
		//volumn defraction
		x02=neighborpoint.GetX0(1);y02=neighborpoint.GetX0(2);//coordinate of the neighbor point
		double r0=sqrt((x02-x01)*(x02-x01)+(y02-y01)*(y02-y01));//bond length
		//double w=1.0;
		
		double vf;double vol;
		if (r0<=Rmin){vf=1.0;}
		else if (r0<Rmax){vf=1.0/2.0+(delta-r0)/rp;}
		else {vf=0.0;}
		vol=rp*rp;
		
		double fx_xp_x,fx_xp_y;//force state, it acts on x, it comes from xp
		double fxp_x_x,fxp_x_y;//force state, it acts on xp, it comes from x
		
		double w=centralpoint.neighbor_w[i];
        if (w==0.0){
			n_broken_bond++;
        }
		ForceStateCal(centralpoint,neighborpoint,fx_xp_x,fx_xp_y,w);
		ForceStateCal(neighborpoint,centralpoint,fxp_x_x,fxp_x_y,w);
		
		double integrandfx=fx_xp_x-fxp_x_x;
		double integrandfy=fx_xp_y-fxp_x_y;
		
		fx=fx+integrandfx*vol*vf;
		fy=fy+integrandfy*vol*vf;
	}

    double den_value=n_neighbors;
    double num_value=n_broken_bond;
	if (n_neighbors == 0){
		centralpoint.broken_parameter=0.0;
	}
	else{
		double f_brok=num_value/den_value;
		centralpoint.broken_parameter=f_brok;
	}


	//assign the point force 
	double pointforcex,pointforcey;
	pointforcex=fx*rp*rp;
	pointforcey=fy*rp*rp;
	centralpoint.SetF(pointforcex,pointforcey);
}


//crack propagation mechanism. prevent the crack heal.
//nei_position is the position of the neighbor point int he neighbor_list of the  central point
//this is used when the neighbor position is known
double ForceState_Elastic::cal_w(double r0,double r,double x02,double y02,Point &centralpoint,int nei_position,Point &neighborpoint){
  double w=centralpoint.neighbor_w[nei_position];
  double omega_D1 = centralpoint.Getomega();
  double omega_D2 = neighborpoint.Getomega();
  //no healing
  if(w==0.0){
    return w;
  }

  else{ 
    //double w=1.0;
    //double delta=centralpoint.Getdelta();
    //double x01=centralpoint.GetX0(1);
    //double y01=centralpoint.GetX0(2);
    //normal calculation
    //double s0=sqrt(5.0*G0/(9.0*k*delta));//critical stretch
	double s0=1.0;
    double boundstretch=(r-r0)/r0;
    if ((boundstretch>=s0)||(omega_D1>=0.75)||(omega_D2>=0.75)){
      w=0.0;
    }else {
      w=1.0;
    }

    //cut the initial edge crack
    //if(x01*x02<0.0){
    //  if (y01<0||y02<0){
    //    w=0.0;
    // }
    //} 
    centralpoint.neighbor_w[nei_position]=w;//update
    return w;
  }
}

//calculate the force state of a central point
/*void ForceState_Elastic::ForceState_compute(Point &centralpoint, vector <Point> &pointgroup){
	double x01,y01,x02,y02,x1,y1,x2,y2;
	x01=centralpoint.GetX0(1);y01=centralpoint.GetX0(2);//coordinate of the subject point
	x1=centralpoint.GetX(1);y1=centralpoint.GetX(2);//deformed coordinate of the subject point

	int n_neighbors=centralpoint.neighbor_list.size();//how many neighbors it has
	vector <double> X; vector <double> Y; //store the deformed coordinate of all bond vectors
	vector <double> R0; //bond length of all the bond vectors
	vector <double> R; //deformed bond length
	vector <double> ee; //bond enlongation of all the bond vectors
	vector <double> Vf;//volumn fraction of neighbor points
	vector <double> Ed;//delatation
	for (int i=0;i<n_neighbors;i++){
		int position=centralpoint.neighbor_list[i];//number of the current neighbor point
		Point neighborpoint=pointgroup[position-1];
		x02=neighborpoint.GetX0(1);y02=neighborpoint.GetX0(2);//coordinate of the neighbor point
		x2=neighborpoint.GetX(1);y2=neighborpoint.GetX(2);//deformed coordinate of the neighbor point
		double r0=sqrt((x02-x01)*(x02-x01)+(y02-y01)*(y02-y01));//bond length
		double r=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));// deformed bond length
		double ebar=r-r0;//bond enlongation
		double x,y;
		x=x2-x1;y=y2-y1;

		//volumn defraction
		double Rmin=delta-rp/2.0;
		double Rmax=delta+rp/2.0;
		double vf;
		if (r0<=Rmin){vf=1.0;}
		else if (r0<Rmax){vf=1.0/2.0+(delta-r)/rp;}
		else {vf=0.0;}

		R0.push_back(r0);//bond length
		R.push_back(r);
		X.push_back(x);//deformed bond vector
		Y.push_back(y);//deformed bond vector
		ee.push_back(ebar);//bond enlongation 
		Vf.push_back(vf);//volumn fraction

	}

	//calculate p and q, two integration
	//q=wx(*)x, p=wx(*)e
	//don't consider fracture, w always equals 1.0
	double vol=rp*rp;//volumn of each point. uniform point size
	double p,q;
	p=0.0;q=0.0;
	for (int i=0;i<n_neighbors;i++){
		double r0=R0[i];
		double vf=Vf[i];
		double ebar=ee[i];
		q=q+1.0*r0*r0*vf*vol;
		p=p+1.0*r0*ebar*vf*vol;
	}
	double alpha=8.0*miu/q;

	//theta=2(2v-1)/(v-1)*p/q
	double scalar;double theta;
	scalar=2.0*(2.0*v-1.0)/(v-1.0);
	theta=p/q;theta=scalar*theta;


	//calculate ed=e-theta*x/3
	for (int i=0;i<n_neighbors;i++){
		double r0=R0[i];
		double ebar=ee[i];
		double ed=ebar-theta*r0/3.0;
		Ed.push_back(ed);
	}
	


	//B=w*ed(*)x=p-theta/3*q
	double B=0.0;
	for (int i=0;i<n_neighbors;i++){
		double w=1.0;//doesn't consider fracture
		double ed=Ed[i];
		double r0=R0[i];
		double vf=Vf[i];
		B=B+w*ed*r0*vf*vol;
	}

	//accumulate force state from all neighbors
	//t=2(2v-1)/(v-1)*(k'*theta-alpha/3*B)*wx/q+alpha*w*ed
	double scalar1;
	scalar1=2.0*(2.0*v-1.0)/(v-1.0)*(kp*theta-alpha/3.0*B)/q;

	double forcescalar;
	double forcestatex=0.0;//initially settle them 0
	double forcestatey=0.0;
	for (int i=0;i<n_neighbors;i++){
		double w=1.0;
		double r0=R0[i];
		double ed=Ed[i];
		forcescalar=scalar1*w*r0+alpha*w*ed;

		double x,y,r,fx,fy;//deformed bond vector
		x=X[i];y=Y[i];r=R[i];
		fx=forcescalar*x/r;
		fy=forcescalar*y/r;
		forcestatex=forcestatex+fx;
		forcestatey=forcestatey+fy;

		//int pointnumber=centralpoint.GetN();
		//int neighbornumber=centralpoint.neighbor_list[i];

	}

	//assign the point force 
	centralpoint.SetF(forcestatex,forcestatey);

}*/


