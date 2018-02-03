// Torque test 3RRR

#include "inertiamodel_3rrr.h"
#include <float.h>
# include <fstream>
using namespace std;

void InertiaModel::torquetest(double xmid, double ymid, double rad, double alpmid, double velmax, int leg)
{
	int N = 200; // Number of points in the tragectory
	
	double angvelmax = velmax/rad;

	double T = 3*M_PI/angvelmax;
	double delt = T/N;
	double t = 0;
	double theta, dtheta, ddtheta;

	double torque[3];

	double task[3], dtask[3], ddtask[3];

	std::ofstream myfile;
	myfile.open ("torqueplot.txt");	

	for(int i=0; i<N; i++)
	{
		theta = ((2*angvelmax*angvelmax)/(3*M_PI))*t*t - ((4*angvelmax*angvelmax*angvelmax)/(27*M_PI*M_PI))*t*t*t;
		dtheta = ((4*angvelmax*angvelmax)/(3*M_PI))*t - ((4*angvelmax*angvelmax*angvelmax)/(9*M_PI*M_PI))*t*t;
		ddtheta = ((4*angvelmax*angvelmax)/(3*M_PI)) - ((8*angvelmax*angvelmax*angvelmax)/(9*M_PI*M_PI))*t;

		// cout<<i<<'\t'<<theta<<'\t'<<dtheta<<'\t'<<ddtheta<<endl;
		// cin.ignore();

		task[0] = X0 + xmid + rad * cos(theta);
		task[1] = Y0 + ymid + rad * sin(theta);
		task[2] = alpmid;

		dtask[0] = - rad * sin(theta) * dtheta;
		dtask[1] = rad * cos(theta) *dtheta;
		dtask[2] = 0;

		ddtask[0] = - rad * sin(theta) * ddtheta - rad * cos(theta) * dtheta;
		ddtask[1] = rad * cos(theta) *ddtheta - rad * sin(theta) * dtheta;
		ddtask[2] = 0;

		// cout<<task[0]<<'\t'<<task[1]<<'\t'<<task[2]<<endl;
		// cout<<dtask[0]<<'\t'<<dtask[1]<<'\t'<<dtask[2]<<endl;	
		// cout<<ddtask[0]<<'\t'<<ddtask[1]<<'\t'<<ddtask[2]<<endl;			
		// cin.ignore();

		torquecalc(task, dtask, ddtask, torque);
		// cout<<theta<<'\t'<<torque[leg]<<endl;

	    // Storing Torque values in text file
	    myfile <<theta<<'\t'<<torque[leg]<<'\n';	

		t+=delt;
	}
    myfile.close();
}