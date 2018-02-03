# include "inertiamodel_3rrr.h"
# include <float.h>

using namespace std;
/*
	dynScan function
	Input: z-range of scanning
	Output: Global dynamic indices
*/
bool InertiaModel::dynScan(double alphamid, double delalpha, double rad, double indices[2])
{
	double eigenval[3];

	int alpiter = floor(2*delalpha);
	double step = 2.5;

	double index1 = 0;
	double index2 = -DBL_MAX;
	int count = 0;

	double temp[2];

	VectorXd q9(9);
	
	double theta[3];
	double phi[3];

	double task[3];


	for(int i = 0; i<alpiter; i+=1)
	{
		task[2] = (alphamid - delalpha + 2*delalpha * i /alpiter)*M_PI/180;	
		for(double j=step; j<=rad; j+= step) // Scanning
		{
			for(double ang=0; ang<=2*M_PI; ang+=step/j)
			{
		    		task[0] = X0 + j*cos(ang)/1000;
				task[1] = Y0 + j*sin(ang)/1000;
				if(  invk(task, theta, phi))
				{
				  std::cout << "Error: Check SWZ (inside dynscan() function)"<< std::endl;
				}
				q9 << theta[0], theta[1], theta[2], phi[0], phi[1], phi[2], task[0], task[1], task[2];

				calcMthetaEigen(q9, eigenval);
				dynindices(eigenval, temp);
				if(!std::isnan(temp[0])){				
					index1 += temp[0];
					count++;
				}
				if(temp[1] > index2)
					index2 = temp[1];
			}	
		}

	}
	indices[0] = index1/count;
	indices[1] = index2;

	return 0;
}



/*
	dynindices function

	Input: Eigenvalues of Mtheta matrix
	Output: Local dynamic indices
*/
void InertiaModel::dynindices(double eigenval[3], double index[2])
{
	// Calculating Index 1
	double inv1 = eigenval[0]+eigenval[1]+eigenval[2];
	double inv2 =  eigenval[0]*eigenval[1] + eigenval[0]*eigenval[2]  + eigenval[1]*eigenval[2];
	double inv3 = eigenval[0]*eigenval[1]*eigenval[2];

	index[0] = 9.*inv3/(inv2*inv1); 
	
	// Calculating Index 2

	if(eigenval[0]>eigenval[1])
	index[1] = eigenval[0];
	else
	index[1] = eigenval[1];

	if(eigenval[2]>index[1])
	index[1] = eigenval[2];

	// std::cout<<index[0]<<'\t'<<index[1]<<std::endl;	
}

/*
	calcMthetaEigen funtion

	Input: Vector q9: configuration space variables
	Output: Eigenvalues of Mtheta matrix 
*/
void InertiaModel::calcMthetaEigen(VectorXd &q9, double eigenval[3])
{

	// Defining variables and matrices

	MatrixXd MqMat(9,9);
	MatrixXd j_q_tht(9,3);

	Mq(q9,MqMat);
	Jqtht(q9, j_q_tht);

	Matrix3d mtheta;
	mtheta = j_q_tht.transpose()*MqMat*j_q_tht;		

	double p1, eig1, eig2, eig3;
	 
	p1 = mtheta(0,1)*mtheta(0,1) + mtheta(0,2)*mtheta(0,2) + mtheta(1,2)*mtheta(1,2);
	double phi;

	if (p1 == 0)
	{
		// A is diagonal.
		eig1 = mtheta(0,0);
		eig2 = mtheta(1,1);
		eig3 = mtheta(2,2);
	}
	else
	{
		double q = (mtheta(0,0)+mtheta(1,1)+mtheta(2,2))/3.;
		double p2 = (mtheta(1,1) - q)*(mtheta(1,1) - q) + (mtheta(2,2) - q)*(mtheta(2,2) - q) + (mtheta(0,0) - q)*(mtheta(0,0) - q) + 2. * p1;
		double p = sqrt(p2 / 6.);
		Matrix3d B = (1 / p) * (mtheta - q * Matrix3d::Identity());
		double r = B.determinant() / 2.;

		// In exact arithmetic for a symmetric matrix  -1 <= r <= 1
		// but computation error can leave it slightly outside this range.
		if (r <= -1)
		{
		phi = M_PI / 3.;
		}
		else if (r >= 1)
		 phi = 0;
		else
		{
		phi = acos(r) / 3.;
		}

		// the eigenvalues satisfy eig3 <= eig2 <= eig1
		eig1 = q + 2. * p * cos(phi);
		eig3 = q + 2. * p * cos(phi + (2.*M_PI/3.));
		eig2 = 3. * q - eig1 - eig3;     //since trace(A) = eig1 + eig2 + eig3
	}

	eigenval[0] =eig1;
	eigenval[1]= eig2;
	eigenval[2] = eig3;	
}
