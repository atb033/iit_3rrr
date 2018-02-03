#ifndef INERTIAMODEL_3rrr_H
#define INERTIAMODEL_3rrr_H 1
#include <iostream>
#include <cmath>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Eigenvalues>

using Eigen::MatrixXd;
using Eigen::Matrix3d;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::Infinity;
class InertiaModel
{
  double g, Eyng, density, P;
  double hmin, mindef, b1, b2, t, f;

  double rb, lcr, lst, ra;
  double X0, Y0;
  double tTp, dt;
  double mass[4];
  double inertia[4];
  double xPl, yPl, zPl;
  double SQRT3;


  void Mq(VectorXd &q9, MatrixXd &MqMat);
  void Cq(VectorXd &q9, VectorXd &dq9, MatrixXd &CqMat);
  void Gq(VectorXd &q9, VectorXd &GqVec);
  void Jqtht(VectorXd &q9,  MatrixXd &J_q_tht);  
  void Jqthtdot(VectorXd &q9, VectorXd &dq9, MatrixXd &J_q_tht_dot);

  bool rrdyad(double, double, double, double, double, double, double *);
  void calc();

  bool invk(double *, double *, double *);
  bool invk1(double *, double *, double *, VectorXd &q9, VectorXd &dq9, VectorXd &ddq9);  

  public:
  
  InertiaModel();

  void output();
  void setlength(double *);
  void setpayload(double, double, double *);
  void torquecalc(double *, double *, double *, double *);
  void torquecalcmax(double *, double *, double *, double &);
  void MqTerm(double *, double *, double *, double *, double &);
  void GqTerm(double *, double *, double *, double *, double &);

  double torquescan(double, double, double, double*);
  void plot(double *, double *, double *, double *, double *);

  void calcMthetaEigen(VectorXd &q9, double *);
  void dynindices(double *, double *);
  bool dynScan(double, double, double, double *);


};
#endif /*INERTIAMODEL_3rrr_H */
