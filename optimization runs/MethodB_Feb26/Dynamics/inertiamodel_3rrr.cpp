#include "inertiamodel_3rrr.h"
using namespace std;
InertiaModel::InertiaModel()
{
  g = 9.812;
  Eyng = 69*1e9; // in Pa
  density = 2700; // in kg/m^3

  mindef = 0.01*1e-3; // minimum deflection
  hmin = 30*1e-3; // minimum height

  b1 = 30*1e-3; //Minimum value of b
  b2 = 10*1e-3;

  t = 5*1e-3; // width of link
  f = 50*1e-3; // space for joint
  tTp = 5*1e-3; // Top platform thickness

  SQRT3 = sqrt(3);
}
void InertiaModel::setlength(double l[])
{
  rb = l[0];
  lcr = l[1];
  lst = l[2];
  ra = l[3];

  X0 = SQRT3*rb*0.5;
  Y0 = 0.5*rb;

}

void InertiaModel::setpayload(double mpl, double ipl, double xypl[2])
{
  // This function accepts payloads and calculates total mass, center of mass and moment of inertia 
  
  xPl = xypl[0];
  yPl = xypl[1];
  
  // Inertia of Payload 
  inertia[3] = ipl;
  mass[3] = mpl; 

  P = mpl * 1.2 * g; // payload mass + link masses (assuming link masses to be 20 % of that of payload)
  calc();
}

void InertiaModel::calc()
{
  
  double hCr, hSt;
  
 // Active link
  hCr = lcr*std::pow(2.0*P/(t*Eyng*mindef),1/3.);  
  if(hCr<hmin)
    hCr = hmin;

  mass[0] = ((lcr+f)*hCr*b1 - (lcr-f)*hCr*(b1-2*t)) * density; 
  inertia[0] =  density/12 *  ( ((lcr+f)*hCr*b1)*((lcr+f)*(lcr+f) + b1 * b1)  -  (lcr-f)*hCr*(b1-2*t) * ((lcr-f)*(lcr-f) + (b1-2*t)*(b1-2*t)));


// Passive link
  hSt = lst*std::pow(2.0*P/(b2*Eyng*mindef),1/3.);  
  if(hSt<hmin)
    hSt = hmin;

  mass[1] = b2*hSt*lst*density; 
  inertia[1] = mass[1]*((lst*lst + b2*b2)/12);

//Top Plate

  mass[2] = 3 * SQRT3*ra * b2 * tTp*density; 
  inertia[2] = mass[2]*(ra*ra - (ra - b2)* (ra - b2))/4;

}

bool InertiaModel::rrdyad(double l1, double l2, double x1, double y1, double x2, double y2, double soln[2])
{
  double temp_1, temp_2, temp_3;
  temp_1 = x2-x1;
  temp_2 = y2-y1;
  temp_3 = temp_1*temp_1+temp_2*temp_2;
  
  soln[0] = -acos((l1*l1+temp_3-l2*l2)/(2*sqrt(temp_3)*l1)) + atan2(temp_2,temp_1);
  soln[1] = atan2(temp_2-l1*sin(soln[0]),temp_1-l1*cos(soln[0]));

  return 0;
}

bool InertiaModel::invk(double task[3], double theta[3], double phi[3])
{
 
  double xpos = task[0];
  double ypos = task[1];
  double alpha = task[2];

  double c_alpha = cos(alpha), s_alpha = sin(alpha);

  double sol1[2], sol2[2], sol3[2];

  double b1[2] = {0, 0};      
  double q1[2]={xpos - (SQRT3*ra*c_alpha)/2. + (ra*s_alpha)/2.,ypos - (ra*c_alpha)/2. - (SQRT3*ra*s_alpha)/2.};

  double b2[2] = {SQRT3*rb, 0};
  double q2[2] ={(2*xpos + SQRT3*ra*c_alpha + ra*s_alpha)/2.,ypos - (ra*c_alpha)/2. + (SQRT3*ra*s_alpha)/2.};

  double b3[2] = { rb*SQRT3/2, rb * 1.5};      
  double q3[2] ={xpos - ra*s_alpha,ypos + ra*c_alpha};
  
  rrdyad(lcr, lst, b1[0], b1[1], q1[0], q1[1], sol1);
  rrdyad(lcr, lst, b2[0], b2[1], q2[0], q2[1], sol2);
  rrdyad(lcr, lst, b3[0], b3[1], q3[0], q3[1], sol3);  
  
  theta[0] = sol1[0]; theta[1] = sol2[0]; theta[2] = sol3[0];
  phi[0] = sol1[1]; phi[1] = sol2[1]; phi[2] = sol3[1];

  return 0;
  
}
bool InertiaModel::invk1(double task[3], double dtask[3], double ddtask[3], VectorXd &q9, VectorXd &dq9, VectorXd &ddq9)
{
  double theta[3], phi[3];

  if(invk(task, theta, phi))
    {
      std::cout << "Error: Check SWZ"<< std::endl;
      return 1;
    }
  double x = task[0], y = task[1], alpha = task[2];
  double dx = dtask[0], dy = dtask[1], dalpha = dtask[2];
  double ddx = ddtask[0], ddy = ddtask[1], ddalpha = ddtask[2];
  
  double tht1 = theta[0];
  double tht2 = theta[1];
  double tht3 = theta[2];

  double phi1 = phi[0];
  double phi2 = phi[1];
  double phi3 = phi[2];
  
  double stht1, stht2, stht3, ctht1, ctht2, ctht3, sphi1, sphi2, sphi3, cphi1, cphi2, cphi3;
  
  stht1 = sin(tht1);
  stht2 = sin(tht2);
  stht3 = sin(tht3);

  ctht1 = cos(tht1);
  ctht2 = cos(tht2);
  ctht3 = cos(tht3);

  sphi1 = sin(phi1);
  sphi2 = sin(phi2);
  sphi3 = sin(phi3);

  cphi1 = cos(phi1);
  cphi2 = cos(phi2);
  cphi3 = cos(phi3);
 
  double calp = cos(alpha), salp = sin(alpha);
  double cscthphi1 = 1/sin(tht1 - phi1), cscthphi2 = 1/sin(tht2 - phi2), cscthphi3 = 1/sin(tht3 - phi3);


  MatrixXd InvJetaBeta(6,6);

  InvJetaBeta << -((cphi1*cscthphi1)/lcr),-((cscthphi1*sphi1)/lcr),0,0,0,0,
              0,0,-((cphi2*cscthphi2)/lcr),-((cscthphi2*sphi2)/lcr),0,0, 
              0,0,0,0,-((cphi3*cscthphi3)/lcr),-((cscthphi3*sphi3)/lcr), 
              (ctht1*cscthphi1)/lst,(cscthphi1*stht1)/lst,0,0,0,0, 
              0,0,(ctht2*cscthphi2)/lst,(cscthphi2*stht2)/lst,0,0, 
              0,0,0,0,(ctht3*cscthphi3)/lst,(cscthphi3*stht3)/lst;
  
  

  MatrixXd JetaX(6,3);
   
  JetaX <<  -1,0,-(ra*(calp + SQRT3*salp))/2., 
            0,-1,(ra*(SQRT3*calp - salp))/2., 
            -1,0,(ra*(-calp + SQRT3*salp))/2., 
            0,-1,-(ra*(SQRT3*calp + salp))/2., 
            -1,0,ra*calp, 
            0,-1,ra*salp;

   Vector3d dX;
   dX <<  dx, dy, dalpha;
   
   VectorXd dBeta(6);
   dBeta = -InvJetaBeta*(JetaX*dX);

   //   VectorXd q9(9);
   q9 << theta[0], theta[1], theta[2], phi[0], phi[1], phi[2],  x, y, alpha;

   // VectorXd dq9(9);
   dq9 << dBeta, dx, dy, dalpha;
   

   double dtht1, dtht2, dtht3, dphi1, dphi2, dphi3;

   dtht1 = dq9(0);
   dtht2 = dq9(1);
   dtht3 = dq9(2);

   dphi1 = dq9(3);
   dphi2 = dq9(4);
   dphi3 = dq9(5);
   
   MatrixXd dJetaBeta(6,6);

   dJetaBeta << -(dtht1*lcr*ctht1),0,0,-(dphi1*lst*cphi1),0,0, 
                -(dtht1*lcr*stht1),0,0,-(dphi1*lst*sphi1),0,0, 
                0,-(dtht2*lcr*ctht2),0,0,-(dphi2*lst*cphi2),0, 
                0,-(dtht2*lcr*stht2),0,0,-(dphi2*lst*sphi2),0, 
                0,0,-(dtht3*lcr*ctht3),0,0,-(dphi3*lst*cphi3), 
                0,0,-(dtht3*lcr*stht3),0,0,-(dphi3*lst*sphi3);

   Vector3d ddX;
   ddX <<  ddx, ddy, ddalpha;

   MatrixXd dJetaX(6,3);

    dJetaX <<   0,0,-(ra*(SQRT3*dalpha*calp - dalpha*salp))/2., 
                0,0,(ra*(-(dalpha*calp) - SQRT3*dalpha*salp))/2., 
                0,0,(ra*(SQRT3*dalpha*calp + dalpha*salp))/2.,
                0,0,-(ra*(dalpha*calp - SQRT3*dalpha*salp))/2.,
                0,0,-(dalpha*ra*salp),
                0,0,dalpha*ra*calp;


   VectorXd ddBeta(6);
   ddBeta = -InvJetaBeta*(JetaX*ddX + dJetaX*dX + dJetaBeta*dBeta);


   // VectorXd ddq9(9);
   ddq9 <<  ddBeta, ddx, ddy, ddalpha;

   return 0;

}

void InertiaModel::output()
{
  cout<<"rb: "<<rb<<endl;
  cout<<"lcr: "<<lcr<<endl;
  cout<<"lst: "<<lst<<endl;
  cout<<"ra: "<<ra<<endl;

  cout<<"mcr: "<<mass[0]<<endl;
  cout<<"mst: "<<mass[1]<<endl;
  cout<<"mEe: "<<mass[2]<<endl;
  cout<<"mPl: "<<mass[3]<<endl;

  cout<<"Icr: "<<inertia[0]<<endl;
  cout<<"Ist: "<<inertia[1]<<endl;
  cout<<"IEe: "<<inertia[2]<<endl;
  cout<<"IPl: "<<inertia[3]<<endl;
}