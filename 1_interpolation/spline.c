#include <stdio.h>
#include <stdlib.h>
#include <math.h>



// Implementation of quadratic spline, "qspline"

// define the qspline structure
typedef struct
{
	int n;
	double x[], y[], b[], c[];
} qspline;



// Implementation of linear spline, "lspline"
double lspline(int n, double x[], double y[], double z)
{
	// check that the values make sense
	assert(n > 1 && x[0]<= z && z <= x[n-1]);

	// binary search for the correct coordinate
	int i = 0, j = n/2, k;
	while(j-i > 1)
	{
		k = (i+j)/2;
		if (z > x[k])
			i = m;
		else
			j = m;
	}

	return y[i] + (y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);
}
