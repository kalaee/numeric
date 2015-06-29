#include <stdio.h>
#include <math.h>
#include "monte_carlo.h"
#include "adapt_2d.h"
#include <gsl/gsl_rng.h>

#define ACC	0.001
#define EPS 0.001

int counter;

double f1_adapt(double x, double y)
{
	counter++;
	return sin(x)*y*exp(-x*y);
}

double f1_mc(double* t)
{
	counter++;
	return sin(t[0])*t[1]*exp(-t[0]*t[1]);
}

int main(void)
{
	// prepare GSL random number generator
	const gsl_rng_type* T;
	gsl_rng* r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	
	double exact, estim, err, a[2], b[2];
	exact = 1.92405;
	a[0] = 0; a[1] = 0;
	b[0] = M_PI; b[1] = 2*M_PI;

	fprintf(stdout,"Tolerances\nACC:\t%g\nEPS:\t%g\n\n",ACC,EPS);

	fprintf(stdout,"Adaptive integrator w/ open intervals in 2D\n");
	counter = 0;
	estim = adapt_2d(f1_adapt,a[0],a[1],b[0],b[1],ACC,EPS,&err);
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n\n",counter);

	fprintf(stdout,"Plain Monte Carlo w/ tolerance\n");
	counter = 0;
	estim = monte_carlo_tolerance(f1_mc,a,b,2,ACC,EPS,100,&err,r);
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n\n",counter);

	gsl_rng_free(r);
}
