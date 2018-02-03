#include <iostream>	
#include <math.h>

using namespace std;

//////////////////////////////////
#define    PI         acos(-1)
#define    EPS        1.0e-10
#define    sign(n)    ((n) < 0) ? -1 : 1
#define    max(A,B)   ((A) > (B)) ? A : B
#define    min(A,B)   ((A) < (B)) ? A : B
//////////////////////////////////
double scanner(double (*function)(double, double, double*), double x0, double y0, double minrad, double circreso, double *parameters)
{
	double r, r_stepsize = circreso;
	double theta, theta_stepsize=PI;

	int initialSign;
	double fnValue;
	double x,y;

	x= x0;
	y= y0;

	fnValue=function(x,y,parameters);
	initialSign= sign(fnValue);

	long fn_samples=0;

	for(r=r_stepsize; r<= minrad; r+=r_stepsize)
	{
		if(theta_stepsize >= PI/180)
			theta_stepsize=circreso/r;
		for(theta=0; theta < 2*PI; theta+=theta_stepsize)
		{
			x= x0 + r*cos(theta);
			y= y0 + r*sin(theta);

			fnValue=function(x,y,parameters);
			fn_samples++;
			// std::cout<<x<<'\t'<<y<<'\t'<<fnValue<<std::endl;
			if( (sign(fnValue)) != initialSign || isnan(fnValue) ) 
			{
				return r - r_stepsize;
			}
		}
	}
	return r-r_stepsize;
}
