/* Test problem definitions */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"
# include "../SWZ/swz_len.h"
# include "../Dynamics/inertiamodel_3rrr.h"

void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
  double l[4];

  l[0] = xreal[0];// crank position
  l[1] = xreal[1];//lCr
  l[2] = xreal[2];//lSt
  l[3] = xreal[3];//rt

  double alphamid = xreal[4]; // alpha mid
  double delalpha = 5;

  // SWZ

  double X0 = l[0]*sqrt(3)/2;
  double Y0 = l[0]*0.5;
  double scanrad = 0.3;
  double swzrad;

  swz_len(l[0], l[1], l[2], l[3], alphamid, delalpha, X0, Y0, scanrad, swzrad);

  // Defining test object
  InertiaModel simulator;

  simulator.setlength(l);

  // Defining Payload
  double mpl = 1.5;
  double ipl = 1;
  double xypl[] = {0, 0};
  simulator.setpayload(mpl, ipl, xypl);

  //Torque scanning zone 1
  double param[2] = {1,1};
  double torqueZone1 = simulator.torquescan(alphamid, delalpha, swzrad*0.2, param);

  //Torque scanning zone 2
  param[0] = xreal[5];
  param[1] = xreal[6];
  double torqueZone2 = simulator.torquescan(alphamid, delalpha, swzrad*0.9, param);

  // Objectives
  obj[0] = -xreal[5];
  obj[1] = -xreal[6];

  // Constraints
  double rreq = 250;
  double torquemax = 3;
  constr[0] = swzrad - rreq; 
  constr[1] = torquemax - torqueZone1;
  constr[2] = torquemax - torqueZone2;

  return;
}
