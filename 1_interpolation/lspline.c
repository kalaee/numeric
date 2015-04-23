#include <assert.h>
#include <stdio.h>

double lspline(int n, double x[], double y[], double z)
{
	assert(n > 1 && x[0]<= z && z <= x[n-1]);

	// binary search for the correct coordinate
	int i = 0, j = n, k;
	while(j-i > 1)
	{
		k = (i+j)/2;
		if (z > x[k])
			i = k;
		else
			j = k;
	}

	return y[i] + (y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);
}
