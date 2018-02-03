#include <iostream>
#include <cmath>
#include "swz_len.h"
#include "swz_3rrr.h"
#include "scanner.h"

using namespace std; 
int main()
{

  double b = 0.794, lcr = 0.346, r = 0.5, a = 0.281;
  double alphamid = -11, delalpha = 5;
  double l[] = {b, lcr, r, a};

  // SWZ
  // center point of the workspace
  double X0 = l[0]*sqrt(3)/2;
  double Y0 = l[0]*0.5;

  double scanrad = 0.40; // Scanning radius

  double swzrad;
  swz_len(l[0], l[1], l[2], l[3], alphamid, delalpha, X0, Y0, scanrad, swzrad);

  cout<<"SWZ Radius: "<<swzrad<<endl;  

	return 0;
}
