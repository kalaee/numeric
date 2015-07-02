#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <gsl/gsl_rng.h>

// routine for integrating n-dimensional iterated integrals using monte carlo integration
// note that because the limits are function of previous coordinates, integration uses
// importance sampling to weight the samples accordingly
// the integrand is function f which takes the dim length array x as argument
// d i lower limits in the integral, thus d(x,p) gives the lower limit for the p'th coordinate, x[p]
// which is a function of coordinates 0, ..., p-1. Likewise u(x,p) givet the upper limit for x[p]
// ACC and EPS are absolute and relative tolerance, THRESHOLD minimum number of samples
// before monte carlo is allowed to terminate and *err pointer to where the estimated error is stored
double monte_carlo_speclim(double f(double* x), double d(double* x, int p), double u(double* x, int p), int dim, double ACC, double EPS, double THRESHOLD, double *err, gsl_rng* r)
{
	int i;
	// unsigned int to increase possible number of samples
	unsigned int n = 0;
	double x[dim], a[dim], h, hi, hx, sum, sum2, fx, avg, I;

	// width of the zeroth coordinate
	a[0] = d(x,0);
	hx = u(x,0) - a[0];

	sum = 0;
	sum2 = 0;

	// sample until desired precision achieved or
	// upper limit UINT_MAX samples reached
	do
	{
		n++;
		// random zeroth coordinate
		x[0] = a[0] + gsl_rng_uniform(r)*hx;
		// for i'th coordinate, find the limits based in subsequent coordinate
		// choose random i'th coordinate and store width h[i]
		// repeat overall coordinates
		// note that the uneven widths mean non-uniform probability distribution
		// in our approach, hence the need for importance sampling
		h = 1;
		for (i=1; i<dim; i++)
		{
			a[i] = d(x,i);
			hi = u(x,i) - a[i];
			h *= hi;
			x[i] = a[i] + gsl_rng_uniform(r)*hi;
		}
		// store integrand value at this point and scale with widths
		// this is where the importance sampling is included
		fx = h*f(x);
		// add weighted integrand values linearly and quadratically
		sum += fx; 
		sum2 += fx*fx;
		avg = sum/n;	// <h_1*...*h_(dim-1)*f(x)>
		*err = hx * sqrt((sum2/n - avg*avg)/n);	// h_0*sqrt( (<(h...hf)^2> - <h...hf>^2) / n)
		I = hx * avg;
		// if number of samples exceeded our ability to count them, terminate monte carlo
		// and only return the unfinished result
		if (n > UINT_MAX-1)
		{
			fprintf(stderr,"Attention, monte_carlo_plain_speclim reached UINT_MAX samples!\n");
			fprintf(stderr,"Integrator will return a result violating allowed tolerance.\n");
			break;
		}
	// if the desired precision is achieved and minimum number of
	// sample points estimated, end integration
	} while (*err > ACC + EPS*I || n < THRESHOLD);
	return I;
}
