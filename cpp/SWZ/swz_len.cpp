/** @file 3rrr_Main.cpp 
The program is used to find designs of 3rrr which satisfy the kinematic design constraints using the concept of "safe working zone" by searching through the design space of its link lengths.
*/

#include <vector>
#include <complex.h>
#include <iostream>
//#include <time.h>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <fstream>

#include "swz_3rrr.h"
#include "scanner.h"
#include "swz_len.h"

#define PAR // Comment out if not parallel version
using namespace std;


//Main function
//----------------------------------------------------------

// The structure of this program involves looping over a range of the design variables and finding the swz data for each of them. The output is written in radData with all the design variables. 

void swz_len(double des_rb, double des_lcr, double des_lst, double des_rt, double des_alpmid, double des_delalp, double X0, double Y0, double scanrad, double &swzrad)
{
  const int fn_num = 4; // The total number of functions to be iterated for the S1-S4 functions.
  const int param_num = 5; // Number of Design parameters.

  const double alpReso = 1; // search resolution for alpha (1 degree)
  const double circReso = 1e-3; // Circumferal resolution of scanning (10 mm) 

  double parameters[param_num];

  int samples = (des_delalp*2)/alpReso;

  double printradii[samples+1], mTOmm = 1000, minRadii;
  double radii[fn_num];
  int i2;

  #ifdef PAR
  #pragma omp parallel for private(parameters,minRadii,radii, i2) shared(printradii)//, zmin, zmax, cSamplesX, cSamplesY, SAMPLESZ)
  #endif
  for(int i=0; i<=samples; i++)
  {
    parameters[0] = des_rb; //crank position
    parameters[1] = des_lcr; //l crank
    parameters[2] = des_lst; //l strut
    parameters[3] = des_rt;  // l top platform
       
    parameters[4] = (des_alpmid - des_delalp + (2*des_delalp) * i/samples) * M_PI/180; //ith z value

    minRadii = scanrad;
    // cout<<(des_alpmid - des_delalp + (2*des_delalp) * i/samples)<<'\t';
    for(i2 = 0; i2<fn_num; i2++)
    {
    	radii[i2]=scanner(fn_ptr_selector(i2), X0, Y0, minRadii, circReso, parameters);
    	// cout<<radii[i2]<<'\t';
      if ( minRadii > radii[i2] ) 
      {
    	  minRadii = radii[i2];
    	}
    	printradii[i] = minRadii*mTOmm;
    }
    // cout<<endl;

  }//END ONE DESIGN	

  swzrad = 1000000;
  for(int i=0; i<=samples; i++)
  {
    if(printradii[i]< swzrad)
      swzrad = printradii[i];
  }
}
