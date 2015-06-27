#include <gsl/gsl_vector.h>
#include "../libs/ode_driver.h"
#include "../libs/ode_rkf45.h"

// solve integrals by rewriting the expression as an ODE system an solve it using ode_evolve
// f is ode function, a is beginning coordinate, b end coordinate, ACC absolute tolerance, EPS relative tolerance
// W ode_workspace with the option RKF45
double ode_integrate(void f(double t, gsl_vector* y, gsl_vector* dy), double a, double b, double ACC, double EPS, ode_workspace* W)
{
	int status;
	// assume initial step size of one hundredth the interval
	double h = (b-a)/100.;
	// ode_evolve uses gsl_vector
	gsl_vector* y = gsl_vector_calloc(1);
	// while loop, terminate when ode_evolve has reached b from a
	while (a < b)
	{
		// take next step
		status = ode_evolve(f,ODE_RKF45,&a,b,&h,y,ACC,EPS,W);
		// if step failed, break loop
		if (status != ODE_SUCCESS)
		{
			break;
		}
	}
	// output end-value of y
	return gsl_vector_get(y,0);
}
