#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdio.h>
#include "../libs/ode_integrate.h"
#include "qaro.h"
#include "qarc.h"
#include "qasc.h"
#include "qaso.h"

#define ACC	1e-9
#define EPS	1e-9
int counter;
// integrand f(x) = cos(x)*exp(1+sin(x))
void f1_ode(double t, gsl_vector* y0, gsl_vector* dydt)
{
	gsl_vector_set(dydt,0,cos(t)*exp(1+sin(t)));
	counter++;
	return;
}

double f1_adapt(double t)
{
	counter++;
	return cos(t)*exp(1+sin(t));
}
int main(void)
{
	counter = 0;
	double estimate, err;
	ode_workspace* W = ode_workspace_alloc(1,ODE_RKF45_ALLOC);

	// integrate function f1
	fprintf(stdout,"True value:\n");
	fprintf(stdout,"int_0^{2*pi} cos(x)*exp(1+sin(x)) = 0\n\n");
	fprintf(stdout,"Routine\tEstimate\tCalls\n");
	
	// RKF45, an ODE routine
	counter = 0;
	estimate = ode_integrate(f1_ode,0,10*M_PI,ACC,EPS,W);
	fprintf(stdout,"RKF45:\t%g\t%d\n",estimate,counter);
	
	// QARO
	counter = 0;
	estimate = qaro(f1_adapt,0,10*M_PI,ACC,EPS,&err);
	fprintf(stdout,"QARO:\t%g\t%d\n",estimate,counter);
	
	// QARC
	counter = 0;
	estimate = qarc(f1_adapt,0,10*M_PI,ACC,EPS,&err);
	fprintf(stdout,"QARC:\t%g\t%d\n",estimate,counter);
	
	// QASO
	counter = 0;
	estimate = qaso(f1_adapt,0,10*M_PI,ACC,EPS,&err);
	fprintf(stdout,"QASO:\t%g\t%d\n",estimate,counter);
	
	// QASC
	counter = 0;
	estimate = qasc(f1_adapt,0,10*M_PI,ACC,EPS,&err);
	fprintf(stdout,"QASC:\t%g\t%d\n",estimate,counter);

	ode_workspace_free(W,ODE_RKF45_FREE);
	return 0;
}
