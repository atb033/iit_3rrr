# include "inertiamodel_3rrr.h"
# include "../SWZ/swz_len.h"
# include <float.h>


using namespace std;

int main()
{
	// Object
	InertiaModel test;

	double b = 0.365, lcr = 0.386, r = 0.414, a = 0.199, alphamid = 128.32; // Torque minimization	

	double l[] = {b, lcr, r, a};

	test.setlength(l);

	// Defining Payload
	double mpl = 1.5;
	double ipl = 1;
	double xypl[] = {0, 0};

	test.setpayload(mpl, ipl, xypl);

	double x0 = 0, y0 = 0;
	double radius = 0.04;
	double velmax = 0.3;

	test.torquetest(x0, y0, radius, alphamid, velmax, 0);
}		