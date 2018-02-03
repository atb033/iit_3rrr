# include "inertiamodel_3rrr.h"
# include "../SWZ/swz_len.h"
# include <float.h>


using namespace std;

int main()
{
	// Object
	InertiaModel test;


	// link lengths
	// double b = 0.577, lcr = 0.5, r = 0.5, a = 0.086, alphamid = 68;
	// double b = 0.733, lcr = 0.375, r = 0.5, a = 0.217, alphamid = -17;
	// double b = 0.844, lcr = 0.360, r = 0.485, a = 0.3, alphamid = -8.82; // Torque minimization
	double b = 0.466, lcr = 0.460, r = 0.499, a = 0.299, alphamid = 115.55;
	double l[] = {b, lcr, r, a};
  	double delalpha = 5;

	test.setlength(l);

	// SWZ Calculation

	// center point of the workspace
	double X0 = l[0]*sqrt(3)/2;
	double Y0 = l[0]*0.5;
	double scanrad = 0.4;
	double swzrad;

	swz_len(l[0], l[1], l[2], l[3], alphamid, delalpha, X0, Y0, scanrad, swzrad);
	// cout<<zrange[0]*180/M_PI<<endl<<zrange[1]*180/M_PI<<endl<<zrange[2]*180/M_PI<<endl;
	// cout<<"SWZ radius :"<<swzrad<<endl;
	// cin.ignore();

	// Defining Payload
	double mpl = 1.5;
	double ipl = 1;
	double xypl[] = {0, 0};

	test.setpayload(mpl, ipl, xypl);
	// test.output();

	double par[] = {1,1}, torque;
	torque = test.torquescan(alphamid, delalpha, 250*0.95, par);

	cout<<"SWZ radius: "<<swzrad<<endl<<"Torque: "<<torque<<endl;


	 // test.output();

	// Dynamic index calculation
//	double index[2];
//	test.dynScan(alphamid, delalpha, swzrad*0.9, index);
//	cout<<"SWZ rad: "<<swzrad<<endl;
//	cout<<"Index 1: "<<index[0]<<endl<<"Index 2: "<<index[1]<<endl; 
	// // cout<<"SWZ: "<<zrange[0]*180/M_PI<<'\t'<<zrange[1]*180/M_PI<<'\t'<<zrange[2]*180/M_PI<<'\n';
	// cout<<"Dynamic Indices: "<<index[0]<<", "<<index[1]<<endl;
	return 0;
}
