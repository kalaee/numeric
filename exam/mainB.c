#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "monte_carlo_speclim.h"


#define ACC	1e-4
#define EPS	1e-4

int counter;

double f1(double* t)
{
	counter++;
	return exp(t[1]-t[0]);
}

double d1(double* t, int p)
{
	if (p == 0)
	{
		return 0;
	}
	else
	{
		return sin(t[0])-0.5;
	}
}

double u1(double* t, int p)
{
	if (p == 0)
	{
		return 2*M_PI;
	}
	else
	{
		return cos(0.5*t[0]*t[0])+0.5;
	}
}


double f2(double* t)
{
	counter++;
	return sin(t[0]+t[1]);
}

double d2(double* t, int p)
{
	if (p == 0)
	{
		return -M_PI/4;
	}
	else
	{
		return 1/cos(t[0])-1;
	}
}

double u2(double* t, int p)
{
	if (p == 0)
	{
		return M_PI/4;
	}
	else if(t[0] == 0)
	{
		return 1;
	}
	else
	{
		return t[0]/tan(t[0]);
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

	double exact, estim, err;
	fprintf(stdout,"Tolerances\nACC:\t%g\nEPS:\t%g\n\n",ACC,EPS);


	fprintf(stdout,"Integrand: exp(y-x),\nx from 0 to 2*pi, y from d(x) to u(x),\n");
	fprintf(stdout,"where d(x) = sin(x)-0.5 and u(x) = cos(0.5*x^2)+0.5\n");
	exact = 2.55398;
	fprintf(stdout,"Actual value (WolframAlpha):\t%g\n",exact);
	counter = 0;
	estim = monte_carlo_speclim(f1,d1,u1,2,ACC,EPS,100,&err,r);
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n\n",counter);


	fprintf(stdout,"Integrand: sin(x+y),\nx from -pi/4 to pi/4, y from d(x) to u(x),\n");
	fprintf(stdout,"where d(x) = 1/cos(x)-1 and u(x) = x/tan(x)\n");
	exact = 0.55951;
	fprintf(stdout,"Actual value (WolframAlpha):\t%g\n",exact);
	counter = 0;
	estim = monte_carlo_speclim(f2,d2,u2,2,ACC,EPS,100,&err,r);
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n\n",counter);

	gsl_rng_free(r);

	return 0;
}
