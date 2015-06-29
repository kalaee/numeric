#include <stdio.h>
#include <math.h>
#include "monte_carlo.h"
#include "adapt_2d.h"
#include <gsl/gsl_rng.h>

#define ACC	0.0001
#define EPS 0.0001

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

double f2_adapt(double x, double y)
{
	counter++;
	double r = sqrt(x*x+y*y);
	if(r < 1)
	{
		return exp(-r);
	}
	else
	{
		return 0;
	}
}

double f2_mc(double* t)
{
	counter++;
	double r = sqrt(t[0]*t[0]+t[1]*t[1]);
	if(r < 1)
	{
		return exp(-r);
	}
	else
	{
		return 0;
	}
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
	exact = 1.924052409955306;
	a[0] = 0; a[1] = 0;
	b[0] = M_PI; b[1] = 2*M_PI;

	fprintf(stdout,"Tolerances\nACC:\t%g\nEPS:\t%g\n\n",ACC,EPS);

	fprintf(stdout,"Integral: int sin(x)*y*exp(-x*y), x from 0 to pi, y from 0 to 2*pi\n");
	fprintf(stdout,"Exact value (Mathematica):\t%g\n\n",exact);

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

	a[0] = -1.1; a[1] = -1.1;
	b[0] = 1.1; b[1] = 1.1;
	exact = 1.66028;
	fprintf(stdout,"Integral: exp(-r) for r = sqrt(x^2+y^2) < 1, zero otherwise\n");
	fprintf(stdout,"x from -1.1 to 1.1, y from -1.1 to 1.1\n");

	fprintf(stdout,"Adaptive integrator w/ open intervals in 2D\n");
	counter = 0;
	estim = adapt_2d(f2_adapt,a[0],a[1],b[0],b[1],ACC,EPS,&err);
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n\n",counter);

	fprintf(stdout,"Plain Monte Carlo w/ tolerance\n");
	counter = 0;
	estim = monte_carlo_tolerance(f2_mc,a,b,2,ACC,EPS,100,&err,r);
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n\n",counter);


	gsl_rng_free(r);
}
