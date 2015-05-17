#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdio.h>
#include "ode_rkf45.h"
#include "ode_driver.h"

#define EPS	1e-12
#define ACC	1e-12


// ODE for sine and cosine, y'' = -y
void f1(double t, gsl_vector* y, gsl_vector* dydt)
{
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-gsl_vector_get(y,0));
}

// Bernoulli's ODE: y' = 2y/t - t^2 y^2
void f2(double t, gsl_vector* y0, gsl_vector* dydt)
{
	double y = gsl_vector_get(y0,0);
	gsl_vector_set(dydt,0,-t*t*y*y+2*y/t);
}

int main(void)
{
	// allocate workspace and system
	ode_workspace* W = ode_workspace_alloc(2,ODE_RKF45_ALLOC);
	int status;
	double t = 0, h;
	gsl_vector* y = gsl_vector_alloc(2);
	gsl_vector_set(y,0,0);
	gsl_vector_set(y,1,1);

	// solve ODE for sine, y(0) = 0, y'(0) = 1
	h = 1;
	fprintf(stdout,"%g\t%g\t%g\n",t,gsl_vector_get(y,0),gsl_vector_get(y,1));
	while(t < 10)
	{
		status = ode_driver(f1,ODE_RKF45,&t,10,&h,y,ACC,EPS,W);
		if (status != ODE_SUCCES)
		{
			break;
		}
		fprintf(stdout,"%g\t%g\t%g\n",t,gsl_vector_get(y,0),gsl_vector_get(y,1));
	}
	// output coordinates for known solution
	fprintf(stdout,"\n\n");
	for(t=0; t<10; t+=1)
	{
		fprintf(stdout,"%g\t%g\t%g\n",t,sin(t),cos(t));
	}

	// Solve Bernoulli ODE, y(0.1) = 0.0125, ignore second coordinate of y vector
	fprintf(stdout,"\n\n");
	gsl_vector_set(y,1,0);
	gsl_vector_set(y,0,0.0125);
	t = 0.1;
	h = 1;
	fprintf(stdout,"%g\t%g\n",t,gsl_vector_get(y,0));
	while(t < 10)
	{
		status = ode_driver(f2,ODE_RKF45,&t,10,&h,y,ACC,EPS,W);
		if (status != ODE_SUCCES)
		{
			break;
		}
		fprintf(stdout,"%g\t%g\n",t,gsl_vector_get(y,0));
	}
	// output known coordinates for known solution
	fprintf(stdout,"\n\n");
	for(t=0.5; t<10; t+=1)
	{
		fprintf(stdout,"%g\t%g\n",t,5*t*t/(pow(t,5)+3.99999));
	}

	// free allocated memory
	gsl_vector_free(y);
	ode_workspace_free(W,ODE_RKF45_FREE);
	return 0;
}
