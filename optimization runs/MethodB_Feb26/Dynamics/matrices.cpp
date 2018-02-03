#include "inertiamodel_3rrr.h"
using namespace std;
void InertiaModel::Mq(VectorXd &q9, MatrixXd &MqMat)
{

  double tht1 = q9(0);
  double tht2 = q9(1);
  double tht3 = q9(2);
  

  double phi1 = q9(3);
  double phi2 = q9(4);
  double phi3 = q9(5);
  
  double alpha = q9(8);

  double ICrz = inertia[0];
  double IStz = inertia[1];
  double ITpz = inertia[2];
  double IPlz = inertia[3];

  double mCr = mass[0];
  double mSt = mass[1];
  double mEe = mass[2];
  double mPl = mass[3];

  double lcrsqr = lcr*lcr;
  double lstsqr = lst*lst;

  double cthphi1 = cos(tht1 - phi1);
  double cthphi2 = cos(tht2 - phi2);
  double cthphi3 = cos(tht3 - phi3);

  double calp = cos(alpha), salp = sin(alpha);

  double xPlsqr = xPl*xPl;
  double yPlsqr = yPl*yPl;

  MqMat <<   
        ICrz + (lcrsqr*(mCr + 4*mSt))/4.,0,0,(lcr*lst*mSt*cthphi1)/2.,0,0,0,0,0,
        0,ICrz + (lcrsqr*(mCr + 4*mSt))/4.,0,0,(lcr*lst*mSt*cthphi2)/2.,0,0,0,0,
        0,0,ICrz + (lcrsqr*(mCr + 4*mSt))/4.,0,0,(lcr*lst*mSt*cthphi3)/2.,0,0,0,
        (lcr*lst*mSt*cthphi1)/2.,0,0,IStz + (lstsqr*mSt)/4.,0,0,0,0,0,
        0,(lcr*lst*mSt*cthphi2)/2.,0,0,IStz + (lstsqr*mSt)/4.,0,0,0,0,
        0,0,(lcr*lst*mSt*cthphi3)/2.,0,0,IStz + (lstsqr*mSt)/4.,0,0,0,
        0,0,0,0,0,0,mEe + mPl,0,-(mPl*(yPl*calp + xPl*salp)),
        0,0,0,0,0,0,0,mEe + mPl,mPl*(xPl*calp - yPl*salp),
        0,0,0,0,0,0,-(mPl*(yPl*calp + xPl*salp)),mPl*(xPl*calp - yPl*salp),IPlz + ITpz + mPl*(xPlsqr + yPlsqr);

}

void InertiaModel::Cq(VectorXd &q9, VectorXd &dq9, MatrixXd &CqMat)
{

//State variables
  double tht1 = q9(0);
  double tht2 = q9(1);
  double tht3 = q9(2);
  
  double phi1 = q9(3);
  double phi2 = q9(4);
  double phi3 = q9(5);
  
  double alpha = q9(8);

// Derivatives of state variables
  double dtht1 = dq9(0);
  double dtht2 = dq9(1);
  double dtht3 = dq9(2);
  
  double dphi1 = dq9(3);
  double dphi2 = dq9(4);
  double dphi3 = dq9(5);

  double dalpha = dq9(8);

//Mass parameters
 
  double mSt = mass[1];
  double mPl = mass[3];

  double calp = cos(alpha);
  double salp = sin(alpha);

  double sthphi1 = sin(tht1 - phi1);
  double sthphi2 = sin(tht2 - phi2);
  double sthphi3 = sin(tht3 - phi3);

  CqMat << 
      0,0,0,(dphi1*lcr*lst*mSt*sthphi1)/2.,0,0,0,0,0,
      0,0,0,0,(dphi2*lcr*lst*mSt*sthphi2)/2.,0,0,0,0,
      0,0,0,0,0,(dphi3*lcr*lst*mSt*sthphi3)/2.,0,0,0,
      -(dtht1*lcr*lst*mSt*sthphi1)/2.,0,0,0,0,0,0,0,0,
      0,-(dtht2*lcr*lst*mSt*sthphi2)/2.,0,0,0,0,0,0,0,
      0,0,-(dtht3*lcr*lst*mSt*sthphi3)/2.,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,-(dalpha*mPl*(xPl*calp - yPl*salp)),
      0,0,0,0,0,0,0,0,dalpha*mPl*(-(yPl*calp) - xPl*salp),
      0,0,0,0,0,0,0,0,0;

}

void InertiaModel::Jqtht(VectorXd &q9, MatrixXd &J_q_tht)
{

  double tht1 = q9(0);
  double tht2 = q9(1);
  double tht3 = q9(2);
  
  double phi1 = q9(3);
  double phi2 = q9(4);
  double phi3 = q9(5);
  
  double alpha = q9(8);
  
  double ctht1 = cos(tht1);
  double ctht2 = cos(tht2);
  double ctht3 = cos(tht3); 

  double cphi1 = cos(phi1);
  double cphi2 = cos(phi2);
  double cphi3 = cos(phi3); 

  double stht1 = sin(tht1);
  double stht2 = sin(tht2);
  double stht3 = sin(tht3); 

  double sphi1 = sin(phi1);
  double sphi2 = sin(phi2);
  double sphi3 = sin(phi3); 

  double calp = cos(alpha);
  double salp = sin(alpha);

  MatrixXd J_eta_tht(6,3);
  MatrixXd J_eta_phi(6,6);

  J_eta_tht << 
  -(lcr*stht1),0,0,
  lcr*ctht1,0,0,
  0,-(lcr*stht2),0,
  0,lcr*ctht2,0,
  0,0,-(lcr*stht3),
  0,0,lcr*ctht3;

  
  J_eta_phi << 
  -(lst*sphi1),0,0,-1,0,-(ra*calp)/2. - (SQRT3*ra*salp)/2.,
  lst*cphi1,0,0,0,-1,(SQRT3*ra*calp)/2. - (ra*salp)/2.,
  0,-(lst*sphi2),0,-1,0,(-(ra*calp) + SQRT3*ra*salp)/2.,
  0,lst*cphi2,0,0,-1,-(SQRT3*ra*calp)/2. - (ra*salp)/2.,
  0,0,-(lst*sphi3),-1,0,ra*calp,
  0,0,lst*cphi3,0,-1,ra*salp;


  MatrixXd J_eta_phi_inv = J_eta_phi.inverse();
  MatrixXd J_phi_tht = (-J_eta_phi_inv*J_eta_tht);

  J_q_tht<<  Matrix3d::Identity(), J_phi_tht;

}


void InertiaModel::Jqthtdot(VectorXd &q9, VectorXd &dq9, MatrixXd &J_q_tht_dot)
{

//State variables
  double tht1 = q9(0);
  double tht2 = q9(1);
  double tht3 = q9(2);
  
  double phi1 = q9(3);
  double phi2 = q9(4);
  double phi3 = q9(5);
  
  double alp = q9(8);

// Derivatives of state variables
  double dtht1 = dq9(0);
  double dtht2 = dq9(1);
  double dtht3 = dq9(2);
  
  double dphi1 = dq9(3);
  double dphi2 = dq9(4);
  double dphi3 = dq9(5);

  double dalp = dq9(8);

  double a = ra;

  J_q_tht_dot << 

  0, 0, 0,
  0, 0, 0,
  0, 0, 0,
  (lcr*((3*(dalp + dphi1 - dphi2 - dphi3)*cos(alp + phi1 - phi2 - phi3) \
- 3*(dalp - dphi1 + dphi2 - dphi3)*cos(alp - phi1 + phi2 - phi3) - \
(dalp + dphi1 - dphi2 - dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - \
(dalp - dphi1 + dphi2 - dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + \
2*(dalp - dphi1 - dphi2 + dphi3)*SQRT3*sin(alp - phi1 - phi2 + \
phi3))*(SQRT3*cos(alp + phi2 - phi3 - tht1) - 2*SQRT3*cos(alp - phi2 \
+ phi3 - tht1) + SQRT3*cos(alp - phi2 - phi3 + tht1) - 3*sin(alp + \
phi2 - phi3 - tht1) + 3*sin(alp - phi2 - phi3 + tht1)) - \
(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + phi2 - \
phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + phi1 - \
phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))*(-3*(dalp + dphi2 - \
dphi3 - dtht1)*cos(alp + phi2 - phi3 - tht1) + 3*(dalp - dphi2 - \
dphi3 + dtht1)*cos(alp - phi2 - phi3 + tht1) - (dalp + dphi2 - dphi3 \
- dtht1)*SQRT3*sin(alp + phi2 - phi3 - tht1) + 2*(dalp - dphi2 + \
dphi3 - dtht1)*SQRT3*sin(alp - phi2 + phi3 - tht1) - (dalp - dphi2 - \
dphi3 + dtht1)*SQRT3*sin(alp - phi2 - phi3 + \
tht1))))/(lst*pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp \
- phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3),2)),(2*lcr*((-dphi2 + dtht2)*cos(phi2 - tht2)*(3*cos(alp - \
phi3) + SQRT3*sin(alp - phi3))*(SQRT3*cos(alp + phi1 - phi2 - phi3) + \
SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + \
phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3)) - (dalp - dphi3)*(SQRT3*cos(alp - phi3) - 3*sin(alp - \
phi3))*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + \
phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))*sin(phi2 - \
tht2) + (3*cos(alp - phi3) + SQRT3*sin(alp - phi3))*(3*(dalp + dphi1 \
- dphi2 - dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - dphi1 + \
dphi2 - dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 - dphi2 \
- dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 + dphi2 \
- dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - dphi1 - \
dphi2 + dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))*sin(phi2 - \
tht2)))/(lst*pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp \
- phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3),2)),(4*lcr*SQRT3*(-((-dphi3 + dtht3)*cos(phi3 - tht3)*sin(alp - \
phi2)*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + \
phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))) + (dalp - \
dphi2)*cos(alp - phi2)*(SQRT3*cos(alp + phi1 - phi2 - phi3) + \
SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + \
phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3))*sin(phi3 - tht3) - sin(alp - phi2)*(3*(dalp + dphi1 - dphi2 - \
dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - dphi1 + dphi2 - \
dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 - dphi2 - \
dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 + dphi2 - \
dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - dphi1 - dphi2 \
+ dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))*sin(phi3 - \
tht3)))/(lst*pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp \
- phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3),2)),



(2*lcr*((-dphi1 + dtht1)*cos(phi1 - tht1)*(-3*cos(alp - phi3) + \
SQRT3*sin(alp - phi3))*(SQRT3*cos(alp + phi1 - phi2 - phi3) + \
SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + \
phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3)) - (dalp - dphi3)*(SQRT3*cos(alp - phi3) + 3*sin(alp - \
phi3))*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + \
phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))*sin(phi1 - \
tht1) + (-3*cos(alp - phi3) + SQRT3*sin(alp - phi3))*(3*(dalp + dphi1 \
- dphi2 - dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - dphi1 + \
dphi2 - dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 - dphi2 \
- dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 + dphi2 \
- dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - dphi1 - \
dphi2 + dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))*sin(phi1 - \
tht1)))/(lst*pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp \
- phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3),2)),(lcr*((3*(dalp + dphi1 - dphi2 - dphi3)*cos(alp + phi1 - \
phi2 - phi3) - 3*(dalp - dphi1 + dphi2 - dphi3)*cos(alp - phi1 + phi2 \
- phi3) - (dalp + dphi1 - dphi2 - dphi3)*SQRT3*sin(alp + phi1 - phi2 \
- phi3) - (dalp - dphi1 + dphi2 - dphi3)*SQRT3*sin(alp - phi1 + phi2 \
- phi3) + 2*(dalp - dphi1 - dphi2 + dphi3)*SQRT3*sin(alp - phi1 - \
phi2 + phi3))*(SQRT3*cos(alp + phi1 - phi3 - tht2) - 2*SQRT3*cos(alp \
- phi1 + phi3 - tht2) + SQRT3*cos(alp - phi1 - phi3 + tht2) + \
3*sin(alp + phi1 - phi3 - tht2) - 3*sin(alp - phi1 - phi3 + tht2)) - \
(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + phi2 - \
phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + phi1 - \
phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))*(3*(dalp + dphi1 - \
dphi3 - dtht2)*cos(alp + phi1 - phi3 - tht2) - 3*(dalp - dphi1 - \
dphi3 + dtht2)*cos(alp - phi1 - phi3 + tht2) - (dalp + dphi1 - dphi3 \
- dtht2)*SQRT3*sin(alp + phi1 - phi3 - tht2) + 2*(dalp - dphi1 + \
dphi3 - dtht2)*SQRT3*sin(alp - phi1 + phi3 - tht2) - (dalp - dphi1 - \
dphi3 + dtht2)*SQRT3*sin(alp - phi1 - phi3 + \
tht2))))/(lst*pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp \
- phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3),2)),(4*lcr*SQRT3*(-((-dphi3 + dtht3)*cos(phi3 - tht3)*sin(alp - \
phi1)*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + \
phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))) + (dalp - \
dphi1)*cos(alp - phi1)*(SQRT3*cos(alp + phi1 - phi2 - phi3) + \
SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + \
phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3))*sin(phi3 - tht3) - sin(alp - phi1)*(3*(dalp + dphi1 - dphi2 - \
dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - dphi1 + dphi2 - \
dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 - dphi2 - \
dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 + dphi2 - \
dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - dphi1 - dphi2 \
+ dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))*sin(phi3 - \
tht3)))/(lst*pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp \
- phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3),2)),


(2*lcr*((-dphi1 + dtht1)*cos(phi1 - tht1)*(-3*cos(alp - phi2) + \
SQRT3*sin(alp - phi2))*(SQRT3*cos(alp + phi1 - phi2 - phi3) + \
SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + \
phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3)) - (dalp - dphi2)*(SQRT3*cos(alp - phi2) + 3*sin(alp - \
phi2))*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + \
phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))*sin(phi1 - \
tht1) + (-3*cos(alp - phi2) + SQRT3*sin(alp - phi2))*(3*(dalp + dphi1 \
- dphi2 - dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - dphi1 + \
dphi2 - dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 - dphi2 \
- dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 + dphi2 \
- dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - dphi1 - \
dphi2 + dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))*sin(phi1 - \
tht1)))/(lst*pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp \
- phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3),2)),(2*lcr*((-dphi2 + dtht2)*cos(phi2 - tht2)*(3*cos(alp - \
phi1) + SQRT3*sin(alp - phi1))*(SQRT3*cos(alp + phi1 - phi2 - phi3) + \
SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + \
phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3)) - (dalp - dphi1)*(SQRT3*cos(alp - phi1) - 3*sin(alp - \
phi1))*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + \
phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))*sin(phi2 - \
tht2) + (3*cos(alp - phi1) + SQRT3*sin(alp - phi1))*(3*(dalp + dphi1 \
- dphi2 - dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - dphi1 + \
dphi2 - dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 - dphi2 \
- dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 + dphi2 \
- dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - dphi1 - \
dphi2 + dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))*sin(phi2 - \
tht2)))/(lst*pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp \
- phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3),2)),(lcr*(-((2*SQRT3*cos(alp - phi1 - phi2 + tht3) - 2*cos(alp \
- tht3)*(SQRT3*cos(phi1 - phi2) + 3*sin(phi1 - phi2)))*(3*(dalp + \
dphi1 - dphi2 - dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - \
dphi1 + dphi2 - dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 \
- dphi2 - dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 \
+ dphi2 - dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - \
dphi1 - dphi2 + dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))) + \
(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + phi2 - \
phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + phi1 - \
phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))*(-2*(dphi1 - \
dphi2)*cos(alp - tht3)*(3*cos(phi1 - phi2) - SQRT3*sin(phi1 - phi2)) \
+ 2*(dalp - dtht3)*(SQRT3*cos(phi1 - phi2) + 3*sin(phi1 - \
phi2))*sin(alp - tht3) - 2*(dalp - dphi1 - dphi2 + \
dtht3)*SQRT3*sin(alp - phi1 - phi2 + tht3))))/(lst*pow(SQRT3*cos(alp \
+ phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + phi2 - phi3) - \
2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + phi1 - phi2 - \
phi3) - 3*sin(alp - phi1 + phi2 - phi3),2)),



(2*lcr*((-dphi1 + dtht1)*cos(phi1 - tht1)*(SQRT3*cos(alp + phi1 - \
phi2 - phi3) + SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp \
- phi1 - phi2 + phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - \
phi1 + phi2 - phi3))*(sin(alp)*(-(SQRT3*cos(phi2)) + \
3*sin(phi2))*sin(phi3) + cos(alp)*(2*cos(phi3)*sin(phi2) + (cos(phi2) \
+ SQRT3*sin(phi2))*sin(phi3))) - (SQRT3*cos(alp + phi1 - phi2 - phi3) \
+ SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 \
+ phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3))*(dphi3*(cos(phi3)*sin(alp)*(-(SQRT3*cos(phi2)) + 3*sin(phi2)) \
+ cos(alp)*(cos(phi2)*cos(phi3) + sin(phi2)*(SQRT3*cos(phi3) - \
2*sin(phi3)))) - dalp*(2*cos(phi3)*sin(alp)*sin(phi2) + \
(SQRT3*cos(alp - phi2) + cos(phi2)*sin(alp) - \
3*cos(alp)*sin(phi2))*sin(phi3)) + dphi2*(sin(alp)*(3*cos(phi2) + \
SQRT3*sin(phi2))*sin(phi3) + cos(alp)*(-(sin(phi2)*sin(phi3)) + \
cos(phi2)*(2*cos(phi3) + SQRT3*sin(phi3)))))*sin(phi1 - tht1) + \
(sin(alp)*(-(SQRT3*cos(phi2)) + 3*sin(phi2))*sin(phi3) + \
cos(alp)*(2*cos(phi3)*sin(phi2) + (cos(phi2) + \
SQRT3*sin(phi2))*sin(phi3)))*(3*(dalp + dphi1 - dphi2 - \
dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - dphi1 + dphi2 - \
dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 - dphi2 - \
dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 + dphi2 - \
dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - dphi1 - dphi2 \
+ dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))*sin(phi1 - \
tht1)))/pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - \
phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3),2),(2*lcr*(-((-dphi2 + dtht2)*cos(phi2 - tht2)*(SQRT3*cos(alp + \
phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + phi2 - phi3) - \
2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + phi1 - phi2 - \
phi3) - 3*sin(alp - phi1 + phi2 - phi3))*(sin(alp)*(SQRT3*cos(phi1) + \
3*sin(phi1))*sin(phi3) + cos(alp)*(2*cos(phi3)*sin(phi1) + (cos(phi1) \
- SQRT3*sin(phi1))*sin(phi3)))) + (SQRT3*cos(alp + phi1 - phi2 - \
phi3) + SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 \
- phi2 + phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + \
phi2 - phi3))*(dalp*(-2*cos(phi3)*sin(alp)*sin(phi1) + (SQRT3*cos(alp \
- phi1) - cos(phi1)*sin(alp) + 3*cos(alp)*sin(phi1))*sin(phi3)) + \
dphi3*(cos(phi3)*sin(alp)*(SQRT3*cos(phi1) + 3*sin(phi1)) + \
cos(alp)*(cos(phi1)*cos(phi3) - sin(phi1)*(SQRT3*cos(phi3) + \
2*sin(phi3)))) + dphi1*(sin(alp)*(3*cos(phi1) - \
SQRT3*sin(phi1))*sin(phi3) - cos(alp)*(sin(phi1)*sin(phi3) + \
cos(phi1)*(-2*cos(phi3) + SQRT3*sin(phi3)))))*sin(phi2 - tht2) - \
(sin(alp)*(SQRT3*cos(phi1) + 3*sin(phi1))*sin(phi3) + \
cos(alp)*(2*cos(phi3)*sin(phi1) + (cos(phi1) - \
SQRT3*sin(phi1))*sin(phi3)))*(3*(dalp + dphi1 - dphi2 - \
dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - dphi1 + dphi2 - \
dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 - dphi2 - \
dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 + dphi2 - \
dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - dphi1 - dphi2 \
+ dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))*sin(phi2 - \
tht2)))/pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - \
phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3),2),(2*lcr*(-((-dphi3 + dtht3)*cos(phi3 - \
tht3)*(cos(alp)*(cos(phi2)*sin(phi1) + (-cos(phi1) + \
2*SQRT3*sin(phi1))*sin(phi2)) - SQRT3*sin(alp)*sin(phi1 + \
phi2))*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + \
phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))) + (-((dphi1 + \
dphi2)*SQRT3*cos(phi1 + phi2)*sin(alp)) - \
dalp*sin(alp)*(cos(phi2)*sin(phi1) + (-cos(phi1) + \
2*SQRT3*sin(phi1))*sin(phi2)) + \
cos(alp)*(-(dphi2*(cos(phi1)*cos(phi2) + \
sin(phi1)*(-2*SQRT3*cos(phi2) + sin(phi2)))) + \
dphi1*(sin(phi1)*sin(phi2) + cos(phi1)*(cos(phi2) + \
2*SQRT3*sin(phi2)))) - dalp*SQRT3*cos(alp)*sin(phi1 + \
phi2))*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + \
phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))*sin(phi3 - \
tht3) - (cos(alp)*(cos(phi2)*sin(phi1) + (-cos(phi1) + \
2*SQRT3*sin(phi1))*sin(phi2)) - SQRT3*sin(alp)*sin(phi1 + \
phi2))*(3*(dalp + dphi1 - dphi2 - dphi3)*cos(alp + phi1 - phi2 - \
phi3) - 3*(dalp - dphi1 + dphi2 - dphi3)*cos(alp - phi1 + phi2 - \
phi3) - (dalp + dphi1 - dphi2 - dphi3)*SQRT3*sin(alp + phi1 - phi2 - \
phi3) - (dalp - dphi1 + dphi2 - dphi3)*SQRT3*sin(alp - phi1 + phi2 - \
phi3) + 2*(dalp - dphi1 - dphi2 + dphi3)*SQRT3*sin(alp - phi1 - phi2 \
+ phi3))*sin(phi3 - tht3)))/pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + \
SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + \
phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3),2),


(2*lcr*(-((-dphi1 + dtht1)*cos(phi1 - tht1)*(SQRT3*cos(alp + phi1 - \
phi2 - phi3) + SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp \
- phi1 - phi2 + phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - \
phi1 + phi2 - phi3))*(cos(alp)*cos(phi3)*(3*cos(phi2) + \
SQRT3*sin(phi2)) + sin(alp)*(cos(phi3)*sin(phi2) + \
cos(phi2)*(-(SQRT3*cos(phi3)) + 2*sin(phi3))))) + (SQRT3*cos(alp + \
phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + phi2 - phi3) - \
2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + phi1 - phi2 - \
phi3) - 3*sin(alp - phi1 + phi2 - \
phi3))*(-(dalp*(cos(phi3)*sin(alp)*(3*cos(phi2) + SQRT3*sin(phi2)) + \
cos(alp)*(-(cos(phi3)*sin(phi2)) + cos(phi2)*(SQRT3*cos(phi3) - \
2*sin(phi3))))) + dphi2*(cos(alp)*cos(phi3)*(SQRT3*cos(phi2) - \
3*sin(phi2)) + sin(alp)*(cos(phi2)*cos(phi3) + \
sin(phi2)*(SQRT3*cos(phi3) - 2*sin(phi3)))) + \
dphi3*(-((SQRT3*cos(alp) + sin(alp))*sin(phi2)*sin(phi3)) + \
cos(phi2)*(2*cos(phi3)*sin(alp) + (-3*cos(alp) + \
SQRT3*sin(alp))*sin(phi3))))*sin(phi1 - tht1) - \
(cos(alp)*cos(phi3)*(3*cos(phi2) + SQRT3*sin(phi2)) + \
sin(alp)*(cos(phi3)*sin(phi2) + cos(phi2)*(-(SQRT3*cos(phi3)) + \
2*sin(phi3))))*(3*(dalp + dphi1 - dphi2 - dphi3)*cos(alp + phi1 - \
phi2 - phi3) - 3*(dalp - dphi1 + dphi2 - dphi3)*cos(alp - phi1 + phi2 \
- phi3) - (dalp + dphi1 - dphi2 - dphi3)*SQRT3*sin(alp + phi1 - phi2 \
- phi3) - (dalp - dphi1 + dphi2 - dphi3)*SQRT3*sin(alp - phi1 + phi2 \
- phi3) + 2*(dalp - dphi1 - dphi2 + dphi3)*SQRT3*sin(alp - phi1 - \
phi2 + phi3))*sin(phi1 - tht1)))/pow(SQRT3*cos(alp + phi1 - phi2 - \
phi3) + SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 \
- phi2 + phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + \
phi2 - phi3),2),(2*lcr*((-dphi2 + dtht2)*cos(phi2 - \
tht2)*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + \
phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3))*(cos(alp)*cos(phi3)*(3*cos(phi1) - SQRT3*sin(phi1)) + \
sin(alp)*(cos(phi3)*sin(phi1) + cos(phi1)*(SQRT3*cos(phi3) + \
2*sin(phi3)))) - (SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp \
- phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3))*(dphi3*((SQRT3*cos(alp) - sin(alp))*sin(phi1)*sin(phi3) + \
cos(phi1)*(2*cos(phi3)*sin(alp) - (3*cos(alp) + \
SQRT3*sin(alp))*sin(phi3))) + dalp*(cos(phi3)*sin(alp)*(-3*cos(phi1) \
+ SQRT3*sin(phi1)) + cos(alp)*(cos(phi3)*sin(phi1) + \
cos(phi1)*(SQRT3*cos(phi3) + 2*sin(phi3)))) - \
dphi1*(cos(alp)*cos(phi3)*(SQRT3*cos(phi1) + 3*sin(phi1)) + \
sin(alp)*(-(cos(phi1)*cos(phi3)) + sin(phi1)*(SQRT3*cos(phi3) + \
2*sin(phi3)))))*sin(phi2 - tht2) + (cos(alp)*cos(phi3)*(3*cos(phi1) - \
SQRT3*sin(phi1)) + sin(alp)*(cos(phi3)*sin(phi1) + \
cos(phi1)*(SQRT3*cos(phi3) + 2*sin(phi3))))*(3*(dalp + dphi1 - dphi2 \
- dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - dphi1 + dphi2 - \
dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 - dphi2 - \
dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 + dphi2 - \
dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - dphi1 - dphi2 \
+ dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))*sin(phi2 - \
tht2)))/pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - \
phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3),2),(2*lcr*((-dphi3 + dtht3)*cos(phi3 - \
tht3)*(cos(phi2)*(SQRT3*cos(alp) - sin(alp))*sin(phi1) + \
cos(phi1)*(-2*SQRT3*cos(phi2)*sin(alp) + (SQRT3*cos(alp) + \
sin(alp))*sin(phi2)))*(SQRT3*cos(alp + phi1 - phi2 - phi3) + \
SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + \
phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3)) - (dphi1*(SQRT3*cos(alp)*cos(phi1 + phi2) - \
sin(alp)*(cos(phi1)*cos(phi2) + sin(phi1)*(-2*SQRT3*cos(phi2) + \
sin(phi2)))) + dphi2*(SQRT3*cos(alp)*cos(phi1 + phi2) + \
sin(alp)*(sin(phi1)*sin(phi2) + cos(phi1)*(cos(phi2) + \
2*SQRT3*sin(phi2)))) - dalp*(cos(alp)*(cos(phi2)*sin(phi1) + \
cos(phi1)*(2*SQRT3*cos(phi2) - sin(phi2))) + SQRT3*sin(alp)*sin(phi1 \
+ phi2)))*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 \
+ phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))*sin(phi3 - \
tht3) + (cos(phi2)*(SQRT3*cos(alp) - sin(alp))*sin(phi1) + \
cos(phi1)*(-2*SQRT3*cos(phi2)*sin(alp) + (SQRT3*cos(alp) + \
sin(alp))*sin(phi2)))*(3*(dalp + dphi1 - dphi2 - dphi3)*cos(alp + \
phi1 - phi2 - phi3) - 3*(dalp - dphi1 + dphi2 - dphi3)*cos(alp - phi1 \
+ phi2 - phi3) - (dalp + dphi1 - dphi2 - dphi3)*SQRT3*sin(alp + phi1 \
- phi2 - phi3) - (dalp - dphi1 + dphi2 - dphi3)*SQRT3*sin(alp - phi1 \
+ phi2 - phi3) + 2*(dalp - dphi1 - dphi2 + dphi3)*SQRT3*sin(alp - \
phi1 - phi2 + phi3))*sin(phi3 - tht3)))/pow(SQRT3*cos(alp + phi1 - \
phi2 - phi3) + SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp \
- phi1 - phi2 + phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - \
phi1 + phi2 - phi3),2),

(4*lcr*((-dphi1 + dtht1)*cos(phi1 - tht1)*sin(phi2 - \
phi3)*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + \
phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3)) - (dphi2 - \
dphi3)*cos(phi2 - phi3)*(SQRT3*cos(alp + phi1 - phi2 - phi3) + \
SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + \
phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3))*sin(phi1 - tht1) + sin(phi2 - phi3)*(3*(dalp + dphi1 - dphi2 - \
dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - dphi1 + dphi2 - \
dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 - dphi2 - \
dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 + dphi2 - \
dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - dphi1 - dphi2 \
+ dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))*sin(phi1 - \
tht1)))/(a*pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - \
phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3),2)),(4*lcr*(-((-dphi2 + dtht2)*cos(phi2 - tht2)*sin(phi1 - \
phi3)*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + \
phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3))) + (dphi1 - \
dphi3)*cos(phi1 - phi3)*(SQRT3*cos(alp + phi1 - phi2 - phi3) + \
SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + \
phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3))*sin(phi2 - tht2) - sin(phi1 - phi3)*(3*(dalp + dphi1 - dphi2 - \
dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - dphi1 + dphi2 - \
dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 - dphi2 - \
dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 + dphi2 - \
dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - dphi1 - dphi2 \
+ dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))*sin(phi2 - \
tht2)))/(a*pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - \
phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3),2)),(4*lcr*((-dphi3 + dtht3)*cos(phi3 - tht3)*sin(phi1 - \
phi2)*(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - phi1 + \
phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + 3*sin(alp + \
phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3)) - (dphi1 - \
dphi2)*cos(phi1 - phi2)*(SQRT3*cos(alp + phi1 - phi2 - phi3) + \
SQRT3*cos(alp - phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + \
phi3) + 3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - \
phi3))*sin(phi3 - tht3) + sin(phi1 - phi2)*(3*(dalp + dphi1 - dphi2 - \
dphi3)*cos(alp + phi1 - phi2 - phi3) - 3*(dalp - dphi1 + dphi2 - \
dphi3)*cos(alp - phi1 + phi2 - phi3) - (dalp + dphi1 - dphi2 - \
dphi3)*SQRT3*sin(alp + phi1 - phi2 - phi3) - (dalp - dphi1 + dphi2 - \
dphi3)*SQRT3*sin(alp - phi1 + phi2 - phi3) + 2*(dalp - dphi1 - dphi2 \
+ dphi3)*SQRT3*sin(alp - phi1 - phi2 + phi3))*sin(phi3 - \
tht3)))/(a*pow(SQRT3*cos(alp + phi1 - phi2 - phi3) + SQRT3*cos(alp - \
phi1 + phi2 - phi3) - 2*SQRT3*cos(alp - phi1 - phi2 + phi3) + \
3*sin(alp + phi1 - phi2 - phi3) - 3*sin(alp - phi1 + phi2 - phi3),2));


}