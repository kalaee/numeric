#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <limits.h>

// This function creates a random sample vector in the volume between a and b
// x is pointer to where the vector is stored, a the beginning of the interval,
// h = b - a, the difference between the two vectors, dim the dimensions of the system
// and r the GSL random number generator
// it works by, for each coordinate, choosing a random point between a and b
// using the uniform distribution in [0;1) provided by gsl_rng_uniform
void monte_carlo_sample_vector(double* x, double* a, double* h, int dim, gsl_rng* r)
{
	int i;
	for(i=0; i<dim; i++)
	{
		x[i] = a[i] + gsl_rng_uniform(r)*h[i];
	}
	return;
}

// plain monte carlo integration of function f with start coordinates a and end corodinates b
// of dimension dim. N is number of vector samples, *err is pointer to where the error is stored
// anf *r is GSL random number generator which the routine need for pseudorandom sampling
double monte_carlo_plain(double f(double* t), double* a, double* b, int dim, int N, double *err, gsl_rng* r)
{
	// prepare routine
	int i;
	double V, h[dim], x[dim], fx, sum, sum2;
	V = 1;
	sum = 0;
	sum2 = 0;

	// estimate the difference vector h
	for(i=0; i<dim; i++)
	{
		h[i] = b[i] - a[i];
		V *= h[i];
	}
	
	// repeat sampling N times, add the function values f(x) in sum
	// and the squared values in sum2
	for(i=0; i<N; i++)
	{
		monte_carlo_sample_vector(x,a,h,dim,r);
		fx = f(x);
		sum += fx;
		sum2 += fx*fx;
	}

	// estimate avg fx
	sum /= N;
	// estimate error
	*err = V*sqrt((sum2/N - sum*sum)/N);
	
	// return integral estimate
	return V*sum;
}

// plain monte carlo integration where the routine collects samples until the desired precision is achieved.
// THRESHOLD is minimum number of attempts needed before allowing the routine to terminate
double monte_carlo_tolerance(double f(double* t), double* a, double* b, int dim, double acc, double eps, int THRESHOLD, double *err, gsl_rng* r)
{
	// prepare routine
	int i;
	double I, V, h[dim], x[dim], fx, sum, sum2, avg;
	V = 1;
	sum = 0;
	sum2 = 0;

	// estimate the difference vector h
	for(i=0; i<dim; i++)
	{
		h[i] = b[i] - a[i];
		V *= h[i];
	}
	
	i = 0;
	// do-while loop continues until desired precision is reached
	// it does, however, require at least THRESHOLD samples before allowing
	// the integration to be considered accomplished
	do
	{
		i++;
		monte_carlo_sample_vector(x,a,h,dim,r);
		fx = f(x);
		sum += fx;
		sum2 += fx*fx;
		avg = sum/i;
		*err = V*sqrt((sum2/i-avg*avg)/i);
		I = V*avg;
		// if we are sampling more points than we can count, send error message
		// and terminate the integration routine
		if (!(i+1 < INT_MAX))
		{
			fprintf(stderr,"Error, monte_carlo_tolerance reached INT_MAX samples!\n");
			break;
		}
	} while (*err > acc + eps*fabs(I) || i < THRESHOLD);
	return I;
}
