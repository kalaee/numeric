#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <gsl/gsl_rng.h>


double monte_carlo_speclim(double f(double* x), double d(double* x, int p), double u(double* x, int p), int dim, double ACC, double EPS, double THRESHOLD, double *err, gsl_rng* r)
{
	unsigned int i, n = 0;
	double x[dim], a[dim], h[dim], sum, sum2, fx, avg, I;

	a[0] = d(x,0);
	h[0] = u(x,0) - a[0];
	sum = 0;
	sum2 = 0;

	do
	{
		n++;
		x[0] = a[0] + gsl_rng_uniform(r)*h[0];
		for (i=1; i<dim; i++)
		{
			a[i] = d(x,i);
			h[i] = u(x,i) - a[i];
			x[i] = a[i] + gsl_rng_uniform(r)*h[i];
		}
		fx = f(x);
		for(i=1; i<dim; i++)
		{
			fx *=h[i];
		}
		sum += fx;
		sum2 += fx*fx;
		avg = sum/n;
		*err = h[0] * sqrt((sum2/n - avg*avg)/n);
		I = h[0] * avg;
		if (n > UINT_MAX-1)
		{
			fprintf(stderr,"Attention, monte_carlo_plain_speclim reached UINT_MAX samples!\n");
			fprintf(stderr,"Integrator will return a result violating allowed tolerance.\n");
			break;
		}
	} while (*err > ACC + EPS*I || n < THRESHOLD);
	return I;
}
