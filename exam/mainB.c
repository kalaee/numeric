#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "monte_carlo_speclim.h"
#include "adapt_nd_speclim.h"

#define ACC	1e-4
#define EPS	1e-4

unsigned int counter;

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
	return cos(t[0])+sin(t[1]+t[2])+t[3];
}
double d2(double* t, int p)
{
	if (p == 0)
	{
		return 0;
	}
	else if (p == 1)
	{
		return -t[0];
	}
	else if (p == 2)
	{
		return 4;
	}
	else
	{
		return -t[2];
	}
}
double u2(double* t, int p)
{
	if (p == 0)
	{
		return 3;
	}
	else if (p == 1)
	{
		return t[0];
	}
	else if (p == 2)
	{
		return t[1] + 1;
	}
	else
	{
		return t[0]+t[1];
	}
}

double f3(double* t)
{
	counter++;
	return sin(t[0]+t[1]+t[2]);
}
double d3(double* t, int p)
{
	if (p == 0)
	{
		return 0;
	}
	else if (p == 1)
	{
		return t[0];
	}
	else
	{
		return t[1]+t[0];
	}
}
double u3(double*t, int p)
{
	if (p == 0)
	{
		return 1;
	}
	else if (p ==1)
	{
		return 2*t[0]*t[0];
	}
	else
	{
		return t[1] - t[0];
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

	fprintf(stdout,"============================================================\n");
	fprintf(stdout,"Integrand: exp(y-x),\nx from 0 to 2*pi, y from d(x) to u(x),\n");
	fprintf(stdout,"where d(x) = sin(x)-0.5 and u(x) = cos(0.5*x^2)+0.5\n");
	exact = 2.55398;
	fprintf(stdout,"============================================================\n");
	fprintf(stdout,"Actual value (Mathematica):\t%g\n\n",exact);
	counter = 0;
	estim = adapt_nd_speclim(f1,d1,u1,2,ACC,EPS,&err);
	fprintf(stdout,"Adaptive n-dim iterated integral routine:\n");
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n",counter);
	fprintf(stdout,"------------------------------------------------------------\n");
	fprintf(stdout,"Monte Carlo n-dim iterated integral:\n");
	counter = 0;
	estim = monte_carlo_speclim(f1,d1,u1,2,ACC,EPS,100,&err,r);
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n\n",counter);

	fprintf(stdout,"============================================================\n");
	fprintf(stdout,"Integrand: cos(x)+sin(y+z)+q,\n");
	fprintf(stdout,"x: d = 0, u = 3\n");
	fprintf(stdout,"y: d(x) = -x, u(x) = x\n");
	fprintf(stdout,"z: d(x,y) = 4, u(x,y) = y+1\n");
	fprintf(stdout,"q: d(x,y,z) = -z, u(x,y,z) = y+x\n");
	exact = 74.325525793728760077;
	fprintf(stdout,"============================================================\n");
	fprintf(stdout,"Actual value (Mathematica):\t%g\n",exact);
	counter = 0;
	estim = adapt_nd_speclim(f2,d2,u2,4,ACC,EPS,&err);
	fprintf(stdout,"Adaptive n-dim iterated integral routine:\n");
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n",counter);
	fprintf(stdout,"------------------------------------------------------------\n");
	fprintf(stdout,"Monte Carlo n-dim iterated integral:\n");
	counter = 0;
	estim = monte_carlo_speclim(f2,d2,u2,4,ACC,EPS,100,&err,r);
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n\n",counter);

	fprintf(stdout,"============================================================\n");
	fprintf(stdout,"Integrand: sin(x+y+z)\n");
	fprintf(stdout,"x: d = 0,\tu = 1\n");
	fprintf(stdout,"y: d(x) = x,\tu(x) = 2*x^2\n");
	fprintf(stdout,"z: d(x,y) = y+x,\tu(x,y) = y - x\n");
	exact = 0.03366664173724119;
	fprintf(stdout,"============================================================\n");
	fprintf(stdout,"Actual value (Mathematica):\t%g\n",exact);
	counter = 0;
	estim = adapt_nd_speclim(f3,d3,u3,3,ACC,EPS,&err);
	fprintf(stdout,"Adaptive n-dim iterated integral routine:\n");
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n",counter);
	fprintf(stdout,"------------------------------------------------------------\n");
	fprintf(stdout,"Monte Carlo n-dim iterated integral:\n");
	counter = 0;
	estim = monte_carlo_speclim(f3,d3,u3,3,ACC,EPS,100,&err,r);
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n\n",counter);
	gsl_rng_free(r);

	return 0;
}
