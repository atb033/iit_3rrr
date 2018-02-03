// Torque calculation functions

#include "inertiamodel_3rrr.h"
#include <float.h>
using namespace std;
#define PAR

/** torquescan Function
Input     |
------------------------------------------------------
zrange    | range of alpha to perform scan
rad       | radius of scanning about neutral position
par       | homogenization parameters
------------------------------------------------------
Output    |
maxTorque | Maximum torque in the given region
--------------------------------------------------------
**/
double InertiaModel::torquescan(double alphamid, double delalpha, double rad, double par[2])
{
//Find z-location of minimum torque

  
  double linvel = 0.2;
  double linacc = 0.3*g;
  double angvel = 10 * M_PI/180;
  double angacc = 100 * M_PI/180;

  double task[] = {X0, Y0, 0};
  double dtask[] = {linvel, 0, angvel};
  double ddtask[] = {par[0]*linacc, 0, par[1]*angacc};

  int alpiter = ceil(2*delalpha);
  double step = 2.5;

  double torque;

  if (rad == 0)
    return DBL_MAX;

  int i;
  double maxTorque[alpiter];

  double j, ang, th;

  #ifdef PAR
  #pragma omp parallel for private(i, j, ang, th,task,ddtask, torque) shared(maxTorque) //, zmin, zmax, cSamplesX, cSamplesY, SAMPLESZ)
  #endif
  for(i = 0; i<=alpiter; i++)
  {
    maxTorque[i]=0;
    task[2] = (alphamid -delalpha +  2*delalpha *  i/alpiter)*M_PI/180;
    for(j = step; j<rad; j+= step)
    {
      for(ang = 0; ang < 2*M_PI; ang = ang + step/j)
      {
        task[0] = X0 + j*cos(ang)/1000;
        task[1] = Y0 + j*sin(ang)/1000;
        torquecalcmax(task, dtask, ddtask, torque);

        for(th = 0; th<=2*M_PI; th+=M_PI/12)
        {
          ddtask[0] = cos(th)*linacc*par[0];
          ddtask[1] = sin(th)*linacc*par[0];
          ddtask[2] = par[1]*angacc;
          // cout<<"x: "<<task[0]<<'\t'<<"y: "<<task[1]<<'\t'<<"alpha:   "<<task[2]<<endl;          

          if(std::isnan(torque))
          {
            cout<<"torque nan";
            cin.ignore();
          }
          if(torque > maxTorque[i])
            maxTorque[i] = torque;            
            // cout<<torque<<'\t'<<maxTorque[i]<<endl;
        }
      }
    }
    // cin.ignore();
  }
  double tor = -DBL_MAX;
  // cin.ignore();
  for (int k = 0; k<=alpiter; k++)
  {
    // cout<<k<<'\t'<<maxTorque[k]<<endl;
    if(tor<maxTorque[k])
      tor = maxTorque[k];
  }

  return tor;
}


/** torquecalcmax Function
Input     |
------------------------------------------------------
task      | task space variables
dtask     | first derivatives of task space variables
ddtask    | second derivatives of task space variables
------------------------------------------------------
Output    |
torque    | Maximum of three actuator torques 
--------------------------------------------------------
**/

void InertiaModel::torquecalcmax(double task[3], double dtask[3], double ddtask[3], double &torque)
{
  // Defining variables and their derivatives
  VectorXd q9(9);
  VectorXd dq9(9);
  VectorXd ddq9(9);

  // Defining Mass, Coriolis, and Gravity Matrices
  MatrixXd MqMat(9,9);
  MatrixXd CqMat(9,9);

  // Torque vector and derivatives of active variables
  VectorXd Toutheta(3);

  VectorXd dtht(3);
  VectorXd ddtht(3);

  MatrixXd J_q_tht(9,3);
  MatrixXd J_q_tht_trans(3,9);
  MatrixXd J_q_tht_dot(9,3);
  

  /********Torque calculation********/

  invk1(task, dtask, ddtask,q9, dq9, ddq9);
  Mq(q9, MqMat);
  Cq(q9, dq9, CqMat);

  Jqtht(q9, J_q_tht);
  Jqthtdot(q9, dq9, J_q_tht_dot);

  J_q_tht_trans = J_q_tht.transpose();
  
  dtht << dq9(0), dq9(1), dq9(2);
  ddtht << ddq9(0), ddq9(1), ddq9(2);
  
  Toutheta = J_q_tht_trans*MqMat*J_q_tht*ddtht + J_q_tht_trans*(MqMat*J_q_tht_dot +  CqMat*J_q_tht)*dtht;

  torque = Toutheta.lpNorm<Infinity>();
}


/** torquecalc Function
Input     |
------------------------------------------------------
task      | task space variables
dtask     | first derivatives of task space variables
ddtask    | second derivatives of task space variables
------------------------------------------------------
Output    |
torque    | Actuator torques 
--------------------------------------------------------
**/

void InertiaModel::torquecalc(double task[3], double dtask[3], double ddtask[3], double torque[3])
{

  // Defining variables and their derivatives
  VectorXd q9(9);
  VectorXd dq9(9);
  VectorXd ddq9(9);

  // Defining Mass, Coriolis, and Gravity Matrices
  MatrixXd MqMat(9,9);
  MatrixXd CqMat(9,9);

  // Torque vector and derivatives of active variables
  VectorXd Toutheta(3);

  VectorXd dtht(3);
  VectorXd ddtht(3);

  MatrixXd J_q_tht(9,3);
  MatrixXd J_q_tht_trans(3,9);
  MatrixXd J_q_tht_dot(9,3);


  /********Torque calculation********/

  invk1(task, dtask, ddtask,q9, dq9, ddq9);

  Mq(q9, MqMat);
  Cq(q9, dq9, CqMat);

  Jqtht(q9, J_q_tht);
  Jqthtdot(q9, dq9, J_q_tht_dot);
  
  J_q_tht_trans = J_q_tht.transpose();
  
  dtht << dq9(0), dq9(1), dq9(2);
  ddtht << ddq9(0), ddq9(1), ddq9(2);
  
  Toutheta = J_q_tht_trans*MqMat*J_q_tht*ddtht + J_q_tht_trans*(MqMat*J_q_tht_dot +  CqMat*J_q_tht)*dtht;

  for(int i = 0; i<3; i++)
  {
    torque[i] = Toutheta(i);
  }

}
