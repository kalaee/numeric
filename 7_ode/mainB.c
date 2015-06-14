#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdio.h>
#include "ode_rk5.h"
#include "ode_rkf45.h"
#include "ode_driver.h"

#define EPS	1e-9
#define ACC	1e-9

int counter;

// ODE for sine and cosine, y'' = -y
void f1(double t, gsl_vector* y, gsl_vector* dydt)
{
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-gsl_vector_get(y,0));
	counter++;
}

// Bernoulli's ODE: y' = 2y/t - t^2 y^2
void f2(double t, gsl_vector* y0, gsl_vector* dydt)
{
	double y = gsl_vector_get(y0,0);
	gsl_vector_set(dydt,0,-t*t*y*y+2*y/t);
	counter++;
}

int main(void)
{
	// allocate workspace and system
	ode_workspace* Wrk5 = ode_workspace_alloc(2,ODE_RK5_ALLOC);
	ode_workspace* Wrkf45 = ode_workspace_alloc(2,ODE_RKF45_ALLOC);
	int i, top;
	double t;
	gsl_matrix* y = gsl_matrix_alloc(1000,3);
	gsl_vector* y0 = gsl_vector_alloc(2);
	gsl_vector_set(y0,0,0);
	gsl_vector_set(y0,1,1);

	// solve ODE for sine, y(0) = 0, y'(0) = 1 using RK5
	counter = 0;
	top = ode_driver(f1,ODE_RK5,0,10,0.1,y0,y,ACC,EPS,Wrk5);
	for(i=0; i<=top; i++)
	{
		fprintf(stdout,"%g\t%g\t%g\n",gsl_matrix_get(y,i,0),
										gsl_matrix_get(y,i,1),
										gsl_matrix_get(y,i,2));
	}

	fprintf(stderr,"Trigonometric system, y'' = -y:\nRK5:\t%d\n",counter);
	
	// solve ODE for sine, y(0) = 0, y'(0) = 1 using RKF45
	gsl_vector_set(y0,0,0);
	gsl_vector_set(y0,1,1);
	counter = 0;
	ode_driver(f1,ODE_RKF45,0,10,0.1,y0,y,ACC,EPS,Wrkf45);
	fprintf(stderr,"RKF45:\t%d\n",counter);

	// output coordinates for known solution
	fprintf(stdout,"\n\n");
	for(t=0; t<10; t+=1)
	{
		fprintf(stdout,"%g\t%g\t%g\n",t,sin(t),cos(t));
	}

	// Solve Bernoulli ODE, y(0.1) = 0.0125 using RK5
	// ignore second coordinate of y vector
	fprintf(stdout,"\n\n");
	gsl_vector_set(y0,1,0);
	gsl_vector_set(y0,0,0.0125);
	counter = 0;
	top = ode_driver(f2,ODE_RK5,0.1,10,0.1,y0,y,ACC,EPS,Wrk5);
	for(i=0; i<=top; i++)
	{
		fprintf(stdout,"%g\t%g\n",gsl_matrix_get(y,i,0),
										gsl_matrix_get(y,i,1));
	}

	fprintf(stderr,"Bernoulli's ODE:\nRK5:\t%d\n",counter);

	// Solve Bernoulli ODE, y(0.1) = 0.0125 using RKF45
	// ignore second coordinate of y vector
	gsl_vector_set(y0,1,0);
	gsl_vector_set(y0,0,0.0125);
	counter = 0;
	ode_driver(f2,ODE_RKF45,0.1,10,0.1,y0,y,ACC,EPS,Wrkf45);
	fprintf(stderr,"RKF45:\t%d\n",counter);

	// output known coordinates for known solution
	fprintf(stdout,"\n\n");
	for(t=0.5; t<10; t+=1)
	{
		fprintf(stdout,"%g\t%g\n",t,5*t*t/(pow(t,5)+3.99999));
	}

	// free allocated memory
	gsl_vector_free(y0);
	gsl_matrix_free(y);
	ode_workspace_free(Wrk5,ODE_RK5_FREE);
	ode_workspace_free(Wrkf45,ODE_RKF45_FREE);
	return 0;
}
