#include <stdio.h>
#include <math.h>

#include "qspline.h"

int main(void)
{
	// generating source data
	int length = 10, i;
	double x[10], y[10];
	FILE * f = fopen("B.dat","w");
	for(i=0; i<10; i++)
	{
		x[i] = i+0.5 * sin(i);
		y[i] = i + cos(i*i);
		fprintf(f,"%g\t%g\n",x[i],y[i]);
	}

	// generate z-coordinates
	int steps = 500;
	double dz = (x[length-1]-x[0])/steps, z[steps];
	for(i=0; i<steps; i++)
	{
		z[i] = x[0]+i*dz;
	}

	// finding the interpolation
	qspline * s = qspline_alloc(length, x, y);
	fprintf(f,"\n\n");
	for(i=0; i<steps; i++)
	{
		fprintf(f,"%g\t%g\n",z[i],qspline_eval(s,z[i]));
	}

	// finding the derivative
	fprintf(f,"\n\n");
	for(i=0; i<steps; i++)
	{
		fprintf(f,"%g\t%g\n",z[i],qspline_derivative(s,z[i]));
	}

	// finding the integral
	fprintf(f,"\n\n");
	for(i=0; i<steps; i++)
	{
		fprintf(f,"%g\t%g\n",z[i],qspline_integral(s,z[i]));
	}

	fclose(f);
	qspline_free(s);

	return 0;
}
