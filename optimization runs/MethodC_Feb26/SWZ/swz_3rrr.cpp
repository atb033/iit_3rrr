
/** @file swz_3rrr.cpp 
The program is used to find the "safe working zone" of the parallel manipulator 3rrr
*/

//#include "global.h"
#include "swz_3rrr.h"
#include <math.h>

const double SQRT3 = sqrt(3);
using namespace std;

/** fn_ptr_selector: Function to select a desired function for the SWZ of 3rrr

@param fn_no Number of the desired function

Variables used:
Variable | Description
-------- |------------
fn       | An array to hold the pointers to the functions

@returns: pointer to the desired function
 */
fptr fn_ptr_selector(int fn_no)
{

  double (*fn[4])(double,double,double*);

  fn[0] = swz_3rrr_S11;
  fn[1] = swz_3rrr_S12;
  fn[2] = swz_3rrr_S13;
  fn[3] = swz_3rrr_S2;

  return fn[fn_no];
}


/** swz_3rrr_invkin_theta1: Inverse Kinematics to find theta1, crank angle of 1st leg

@param xpos The x-coordinate of centroid of the platform
@param ypos The y-coordinate of centroid of the platform
@param theta1 The variable to store the solution

Variables used:
Variable | Description
-------- |-------------------------------------
rb       | Base Radius 
ra       | Radius of end effector plate
lcr      | Length of the crank
lst      | Length of the strut 
-----------------------------------------------


@returns: 0, if solved properly

*/

int swz_3rrr_invkin_theta1(double xpos, double ypos, double *params, double &theta1)
{
  double lcr, lst, ra, alpha;
  
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  alpha = params[4];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);

  double b1[2] = {0, 0};      
  double q1[2]={xpos - (SQRT3*ra*c_alpha)/2. + (ra*s_alpha)/2.,ypos - (ra*c_alpha)/2. - (SQRT3*ra*s_alpha)/2.};

  swz_3rrr_rrdyad1(lcr, lst, b1[0], b1[1], q1[0], q1[1], theta1);
  return 0;
}

/** swz_3rrr_invkin_theta2: Inverse Kinematics to find theta1, crank angle of 2nd leg

@param xpos The x-coordinate of centroid of the platform
@param ypos The y-coordinate of centroid of the platform
@param theta2 The variable to store the solution

Variables used:
Variable | Description
-------- |-------------------------------------
rb       | Base Radius 
ra       | Radius of end effector plate
lcr      | Length of the crank
lst      | Length of the strut 
-----------------------------------------------


@returns: 0, if solved properly

*/
int swz_3rrr_invkin_theta2(double xpos, double ypos, double *params, double &theta2)
{
  double rb, lcr, lst, ra, alpha;
  
  rb = params[0]; 
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  alpha = params[4];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);

  double b2[2] = {SQRT3*rb, 0};
  double q2[2]={(2*xpos + SQRT3*ra*c_alpha + ra*s_alpha)/2.,ypos - (ra*c_alpha)/2. + (SQRT3*ra*s_alpha)/2.};
  
  swz_3rrr_rrdyad1(lcr, lst, b2[0], b2[1], q2[0], q2[1], theta2);

  return 0;
}

/** swz_3rrr_invkin_theta2: Inverse Kinematics to find theta1, crank angle of 3rd leg

@param xpos The x-coordinate of centroid of the platform
@param ypos The y-coordinate of centroid of the platform
@param theta3 The variable to store the solution

Variables used:
Variable | Description
-------- |-------------------------------------
rb       | Base Radius 
ra       | Radius of end effector plate
lcr      | Length of the crank
lst      | Length of the strut 
-----------------------------------------------


@returns: 0, if solved properly

*/
int swz_3rrr_invkin_theta3(double xpos, double ypos, double *params, double &theta3)
{
  double rb, lcr, lst, ra, alpha;
  
  rb = params[0]; 
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  alpha = params[4];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);

  double b3[2] = { rb*SQRT3/2, rb * 1.5};      
  double q3[2]={xpos - ra*s_alpha,ypos + ra*c_alpha};
  
  swz_3rrr_rrdyad1(lcr, lst, b3[0], b3[1], q3[0], q3[1], theta3);

  return 0;
}

/** swz_3rrr_invkin1: Inverse Kinematics to find both crank angle and strut angle of first leg

@param xpos The x-coordinate of centroid of the platform
@param ypos The y-coordinate of centroid of the platform
@param sol1 The variable to store the solution

Variables used:
Variable | Description
-------- |-------------------------------------
rb       | Base Radius 
ra       | Radius of end effector plate
lcr      | Length of the crank
lst      | Length of the strut 
-----------------------------------------------


@returns: 0, if solved properly

*/
int swz_3rrr_invkin1(double xpos, double ypos, double *params, double sol1[2])
{
  double lcr, lst, ra, alpha;
  
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  alpha = params[4];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);

  double b1[2] = {0, 0};      
  double q1[2]={xpos - (SQRT3*ra*c_alpha)/2. + (ra*s_alpha)/2.,ypos - (ra*c_alpha)/2. - (SQRT3*ra*s_alpha)/2.};

  swz_3rrr_rrdyad(lcr, lst, b1[0], b1[1], q1[0], q1[1], sol1);

  return 0;
}

/** swz_3rrr_invkin2: Inverse Kinematics to find both crank angle and strut angle of second leg

@param xpos The x-coordinate of centroid of the platform
@param ypos The y-coordinate of centroid of the platform
@param sol2 The variable to store the solution

Variables used:
Variable | Description
-------- |-------------------------------------
rb       | Base Radius 
ra       | Radius of end effector plate
lcr      | Length of the crank
lst      | Length of the strut 
-----------------------------------------------


@returns: 0, if solved properly

*/
int swz_3rrr_invkin2(double xpos, double ypos, double *params, double sol2[2])
{
  double rb, lcr, lst, ra, alpha;
  
  rb = params[0]; 
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  alpha = params[4];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);

  double b2[2] = {SQRT3*rb, 0};
  double q2[2] ={(2*xpos + SQRT3*ra*c_alpha + ra*s_alpha)/2.,ypos - (ra*c_alpha)/2. + (SQRT3*ra*s_alpha)/2.};
  
  swz_3rrr_rrdyad(lcr, lst, b2[0], b2[1], q2[0], q2[1], sol2);

  return 0;
}

/** swz_3rrr_invkin3: Inverse Kinematics to find both crank angle and strut angle of third leg

@param xpos The x-coordinate of centroid of the platform
@param ypos The y-coordinate of centroid of the platform
@param sol2 The variable to store the solution

Variables used:
Variable | Description
-------- |-------------------------------------
rb       | Base Radius 
ra       | Radius of end effector plate
lcr      | Length of the crank
lst      | Length of the strut 
-----------------------------------------------


@returns: 0, if solved properly

*/
int swz_3rrr_invkin3(double xpos, double ypos, double *params, double sol3[2])
{
  double rb, lcr, lst, ra, alpha;
  
  rb = params[0]; 
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  alpha = params[4];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);

  double b3[2] = { rb*SQRT3/2, rb * 1.5};      
  double q3[2] ={xpos - ra*s_alpha,ypos + ra*c_alpha};
  
  swz_3rrr_rrdyad(lcr, lst, b3[0], b3[1], q3[0], q3[1], sol3);

  return 0;
}

/** S1 function for 1st leg: S11

@param xpos The x-coordinate of centroid of the platform
@param ypos The y-coordinate of centroid of the platform
@param value The variable to store the solution

Variables used:
Variable | Description
-------- |-------------------------------------
rb       | Base Radius 
ra       | Radius of end effector plate
lcr      | Length of the crank
lst      | Length of the strut 
-----------------------------------------------


@returns: value of s function, if solved properly

*/

double swz_3rrr_S11(double xpos, double ypos, double *params)
{
  double lcr, lst, ra, alpha;

  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  alpha = params[4];

  double theta1;
  if(swz_3rrr_invkin_theta1(xpos, ypos, params, theta1))
    std::cout<<"Error: Check S11";   

  double c_alpha = cos(alpha), s_alpha = sin(alpha);
  double c_theta1 = cos(theta1), s_theta1 = sin(theta1);

  double value = -2*ypos*c_theta1 + ra*c_alpha*c_theta1 + SQRT3*ra*c_theta1*s_alpha + lst*s_theta1 + 2*lcr*c_theta1*s_theta1;

  return value;

}

/** S1 function for 2nd leg: S12

@param xpos The x-coordinate of centroid of the platform
@param ypos The y-coordinate of centroid of the platform
@param value The variable to store the solution

Variables used:
Variable | Description
-------- |-------------------------------------
rb       | Base Radius 
ra       | Radius of end effector plate
lcr      | Length of the crank
lst      | Length of the strut 
-----------------------------------------------


@returns: value of s function, if solved properly

*/

double swz_3rrr_S12(double xpos, double ypos, double *params)
{
  double lcr, lst, ra, alpha;

  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  alpha = params[4];

  double theta2;
  if(swz_3rrr_invkin_theta2(xpos, ypos, params, theta2))
    std::cout<<"Error: Check S12";   

  double c_alpha = cos(alpha), s_alpha = sin(alpha);
  double c_theta2 = cos(theta2), s_theta2 = sin(theta2);

  double value = -2*ypos*c_theta2 + ra*c_alpha*c_theta2 - SQRT3*ra*c_theta2*s_alpha + lst*s_theta2 + 2*lcr*c_theta2*s_theta2;

//  cout<<value;

  return value;
}

/** S1 function for 3rd leg: S13

@param xpos The x-coordinate of centroid of the platform
@param ypos The y-coordinate of centroid of the platform
@param value The variable to store the solution

Variables used:
Variable | Description
-------- |-------------------------------------
rb       | Base Radius 
ra       | Radius of end effector plate
lcr      | Length of the crank
lst      | Length of the strut 
-----------------------------------------------


@returns: value of s function, if solved properly

*/

double swz_3rrr_S13(double xpos, double ypos, double *params)
{
  double rb, lcr, lst, ra, alpha;

  rb = params[0]; 
  lcr = params[1]; 
  lst = params[2];
  ra = params[3]; 
  alpha = params[4];

  double theta3;
  if(swz_3rrr_invkin_theta3(xpos, ypos, params, theta3))
    std::cout<<"Error: Check S13";   

  double c_alpha = cos(alpha);
  double c_theta3 = cos(theta3), s_theta3 = sin(theta3);

  double value = 3*rb*c_theta3 - 2*ypos*c_theta3 - 2*ra*c_alpha*c_theta3 + lst*s_theta3 + 2*lcr*c_theta3*s_theta3;

  return value;
}

/** S2 function : S2

@param xpos The x-coordinate of centroid of the platform
@param ypos The y-coordinate of centroid of the platform
@param  The variable to store the solution

Variables used:
Variable | Description
-------- |-------------------------------------
rb       | Base Radius 
ra       | Radius of end effector plate
lcr      | Length of the crank
lst      | Length of the strut 
-----------------------------------------------


@returns: value of s function, if solved properly

*/

double swz_3rrr_S2(double xpos, double ypos, double *params)
{
  double rb, lcr, lst;

  rb = params[0]; 
  lcr = params[1]; 
  lst = params[2];

  double sol1[2], sol2[2], sol3[2];

  if(swz_3rrr_invkin1(xpos, ypos, params, sol1))
    std::cout<<"ERROR: Check inside S2";
  if(swz_3rrr_invkin2(xpos, ypos, params, sol2))
    std::cout<<"ERROR: Check inside S2";
  if(swz_3rrr_invkin3(xpos, ypos, params, sol3))
    std::cout<<"ERROR: Check inside S2";      

  double theta1 = sol1[0], theta2 = sol2[0], theta3 = sol3[0];
  double phi1 = sol1[1], phi2 = sol2[1], phi3 = sol3[1];

  double value = (lcr*sin(theta1 - phi2) - lcr*sin(theta2 - phi2) + lst*sin(phi1 - \
phi2) + SQRT3*rb*sin(phi2))*(-3*rb*cos(phi1) + 2*lcr*sin(theta1 - \
phi1) - 2*lcr*sin(theta3 - phi1) + SQRT3*rb*sin(phi1) + \
2*lst*sin(phi1 - phi3))*(-3*rb*cos(phi3) + 2*lcr*sin(theta2 - phi3) - \
2*lcr*sin(theta3 - phi3) + 2*lst*sin(phi2 - phi3) - \
SQRT3*rb*sin(phi3)) - (lcr*sin(theta1 - phi1) - lcr*sin(theta2 - \
phi1) + SQRT3*rb*sin(phi1) + lst*sin(phi1 - phi2))*(-3*rb*cos(phi2) + \
2*lcr*sin(theta2 - phi2) - 2*lcr*sin(theta3 - phi2) - \
SQRT3*rb*sin(phi2) + 2*lst*sin(phi2 - phi3))*(-3*rb*cos(phi3) + \
2*lcr*sin(theta1 - phi3) - 2*lcr*sin(theta3 - phi3) + 2*lst*sin(phi1 \
- phi3) + SQRT3*rb*sin(phi3));

  return value;
}