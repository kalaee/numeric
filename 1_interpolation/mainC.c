#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "cspline.h"

int main(void)
{
	// generating source data
	int length = 10, i;
	double x[length], y[length];
	FILE * f = fopen("C.dat","w");
	for(i=0; i<length; i++)
	{
		x[i] = i+0.5 * sin(i);
		y[i] = i + cos(i*i);
		fprintf(f,"%g\t%g\n",x[i],y[i]);
	}

	int steps = 1000;
	double dz = (x[length-1]-x[0])/steps, z;
	
	// finding creating the cspline with second derivatives 3 and -2 at the end points respectively
	cspline * s = cspline_alloc(length, x, y,3,-2);
	fprintf(f,"\n\n");
	for(i=0; i<steps; i++)
	{
		z = x[0]+i*dz;
		fprintf(f,"%g\t%g\t%g\t%g\n",z,cspline_eval(s,z),cspline_derivative(s,z),cspline_integral(s,z));
	}
	
	fclose(f);
	cspline_free(s);

	return 0;
}
