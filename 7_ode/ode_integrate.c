#include <gsl/gsl_vector.h>
#include "ode_driver.h"
#include "ode_rkf45.h"

double ode_integrate(void f(double t, gsl_vector* y, gsl_vector* dy), double a, double b, double ACC, double EPS, ode_workspace* W)
{
	int status;
	double h = (b-a)/100.;
	gsl_vector* y = gsl_vector_calloc(1);
	while (a < b)
	{
		status = ode_driver(f,ODE_RKF45,&a,b,&h,y,ACC,EPS,W);
		if (status != ODE_SUCCESS)
		{
			break;
		}
	}
	return gsl_vector_get(y,0);
}
