#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "lspline.h"
#include "qspline.h"

int main(void)
{
	// generating source data
	int length = 10, i;
	double x[10], y[10];
	FILE * f = fopen("A.dat","w");
	for(i=0; i<10; i++)
	{
		x[i] = i+0.5 * sin(i);
		y[i] = i + cos(i*i);
		fprintf(f,"%g\t%g\n",x[i],y[i]);
	}

	// applying linear spline
	int steps = 500;
	double dz = (x[length-1]-x[0])/steps, z[steps];
	fprintf(f,"\n\n");
	for(i=0; i<steps; i++)
	{
		z[i] = x[0]+i*dz;
		fprintf(f,"%g\t%g\n",z[i],lspline(length,x,y,z[i]));
	}
	
	// applying quadratic spline
	qspline * s = qspline_alloc(length, x, y);
	fprintf(f,"\n\n");
	for(i=0; i<steps; i++)
	{
		fprintf(f,"%g\t%g\n",z[i],qspline_eval(s,z[i]));
	}
	fclose(f);
	qspline_free(s);

	return 0;
}

