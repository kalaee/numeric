#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_integration.h>

#include "cspline.h"

double fun(double z, void * params)
{
	cspline* s = (cspline*) params;
	return cspline_eval(s,z);
}

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
	// creating the cspline with second derivatives 3 and -2 at the end points respectively
	cspline * s = cspline_alloc(length, x, y,3,-2);
	fprintf(f,"\n\n");
	for(i=0; i<steps; i++)
	{
		z = x[0]+i*dz;
		fprintf(f,"%g\t%g\t%g\t%g\n",z,cspline_eval(s,z),cspline_derivative(s,z),cspline_integral(s,z));
	}

	// GSL integration
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	double result, error;
	gsl_function F;
	F.function = &fun;
	F.params = s;
	fprintf(f,"\n\n");
	dz = (x[length-1]-x[0])/10;
	for(i=0; i<10; i++)
	{
		z = x[0]+i*dz;
		gsl_integration_qag(&F,x[0],z,1e-7,1e-7,1000,6,w,&result,&error);
		fprintf(f,"%g\t%g\n",z,result);
	}
	fclose(f);
	cspline_free(s);
	gsl_integration_workspace_free(w);

	return 0;
}
