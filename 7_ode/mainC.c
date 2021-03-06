#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdio.h>
#include "ode_integrate.h"

#define ACC	1e-9
#define EPS	1e-9

// integrand f(x) = log(1+tan(t))
void f1(double t, gsl_vector* y0, gsl_vector* dydt)
{
	gsl_vector_set(dydt,0,log(1+tan(t)));
}

// integrand f(x) = log(cos(x))
void f2(double t, gsl_vector* y0, gsl_vector* dydt)
{
	gsl_vector_set(dydt,0,log(cos(t)));
}

// integrand f(x) = sqrt(1-x^2)
void f3(double t, gsl_vector* y0, gsl_vector* dydt)
{
	gsl_vector_set(dydt,0,sqrt(1-t*t));
}

int main(void)
{
	double estimate, trueval;
	ode_workspace* W = ode_workspace_alloc(1,ODE_RKF45_ALLOC);

	// integrate function f1
	fprintf(stdout,"int_0^{pi/4} log(1+tan(t))\n");
	estimate = ode_integrate(f1,0,M_PI/4.,ACC,EPS,W);
	trueval =M_PI/8.*log(2);
	fprintf(stdout,"Numerical estimate (RKF45): %g\n",estimate);
	fprintf(stdout,"True value (Schaum's): log(2)*pi/8 =\t%g\n",trueval);
	fprintf(stdout,"Relative deviation: %g\n",estimate/trueval-1);

	// integrate function f2
	fprintf(stdout,"\nint_0^{pi/2} log(cos(t))\n");
	estimate = ode_integrate(f2,0,M_PI/2.,ACC,EPS,W);
	trueval = -M_PI/2.*log(2);
	fprintf(stdout,"Numerical estimate (RKF45): %g\n",estimate);
	fprintf(stdout,"True value (Schaum's): -log(2)*pi/2 = \t%g\n",trueval);
	fprintf(stdout,"Relative deviation: %g\n",estimate/trueval-1);

	// integrate function f3
	fprintf(stdout,"\nint_0^{1} sqrt(1-x^2)\n");
	estimate = ode_integrate(f3,0,1,ACC,EPS,W);
	trueval = M_PI/4;
	fprintf(stdout,"Numerical estimate (RKF45): %g\n",estimate);
	fprintf(stdout,"True value (Schaum's): pi/4 = \t%g\n",trueval);
	fprintf(stdout,"Relative deviation: %g\n",estimate/trueval-1);

	ode_workspace_free(W,ODE_RKF45_FREE);
	return 0;
}
