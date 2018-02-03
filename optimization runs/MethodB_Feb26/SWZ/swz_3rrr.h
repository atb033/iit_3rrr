/** @file swz_3rrr.h 
The program is used to find the "safe working zone" of a 3rrr-parallel manipulator.
*/
#ifndef SWZ_3rrr_H
#define SWZ_3rrr_H 1
#include <cmath>
#include <iostream>

typedef double (*fptr)(double alpha, double beta, double* params);
fptr fn_ptr_selector(int fn_no);

inline void swz_3rrr_rrdyad(double l1, double l2, double x1, double y1, double x2, double y2, double soln[2]);
inline void swz_3rrr_rrdyad1(double l1, double l2, double x1, double y1, double x2, double y2, double &soln);

int swz_3rrr_invkin_theta1(double , double , double *, double &);
int swz_3rrr_invkin_theta2(double , double , double *, double &);
int swz_3rrr_invkin_theta3(double , double , double *, double &);

int swz_3rrr_invkin1(double , double , double *, double *);
int swz_3rrr_invkin2(double , double , double *, double *);
int swz_3rrr_invkin3(double , double , double *, double *);

double swz_3rrr_S11(double , double , double *);
double swz_3rrr_S12(double , double , double *);
double swz_3rrr_S13(double , double , double *);

double swz_3rrr_S2(double alpha, double beta, double* params);


/** swz_3rrr_rrdyad: Inverse Kinematics for 3rrr for finding the angles of the RRDyad of each leg

@param l1 The length of link1 in rrdyad
@param l2 The length of link2 in rrdyad
@param x1 x-position of rrdyad base
@param x2 x-position of rrdyad end
@param y1 y-position of rrdyad base
@param y2 y-position of rrdyad end
@param soln variable to store the solution in the order, theta and gamma

Variables used:
Variable | Description
-------- |------------
temp_i   | Temporary variable used for expressions evaluated repeatedly

@returns: void
 */
inline void swz_3rrr_rrdyad(double l1, double l2, double x1, double y1, double x2, double y2, double soln[2])
{
  double temp_1, temp_2, temp_3;
  temp_1 = x2-x1;
  temp_2 = y2-y1;
  temp_3 = temp_1*temp_1+temp_2*temp_2;
  soln[0] = -acos((l1*l1+temp_3-l2*l2)/(2*sqrt(temp_3)*l1)) + atan2(temp_2,temp_1);
  soln[1] = atan2(temp_2-l1*sin(soln[0]),temp_1-l1*cos(soln[0]));

}

/** swz_3rrr_rrdyad1: Inverse Kinematics for 3rrr for finding only theta of the RRDyad of each leg

@param l1 The length of link1 in rrdyad
@param l2 The length of link2 in rrdyad
@param x1 x-position of rrdyad base
@param x2 x-position of rrdyad end
@param y1 y-position of rrdyad base
@param y2 y-position of rrdyad end
@param soln variable to store the solution for theta

Variables used:
Variable | Description
@returns: void
 */
inline void swz_3rrr_rrdyad1(double l1, double l2, double x1, double y1, double x2, double y2, double &soln)
{
  double temp_1, temp_2, temp_3;
  temp_1 = x2-x1;
  temp_2 = y2-y1;
  temp_3 = temp_1*temp_1+temp_2*temp_2;
  soln = -acos((l1*l1+temp_3-l2*l2)/(2*sqrt(temp_3)*l1)) + atan2(temp_2,temp_1);
}


#endif//SWZ_3rrr_H
