#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include "monte_carlo_speclim.h"
#include "adapt_nd_speclim.h"


#define ACC	1e-3
#define EPS	1e-3

unsigned long long int counter;
int dim;

// for volume of n-ball
double f (double* t)
{
	counter++;
	return 1;
}
// for exponential function in n-ball
double fexp (double* t)
{
	counter++;
	double x = 0;
	int i;
	for(i=0; i<dim; i++)
	{
		x += t[i];
	}
	return exp(x);
}
// lower limits
double d(double* t, int p)
{
	int i;
	double x2 = 0;
	for(i = 0; i<p; i++)
	{
		x2 += t[i]*t[i];
	}
	return -sqrt(1-x2);
}
// upper limits
double u(double* t, int p)
{
	return -d(t,p);
}

// volume of unit ball in n-dim space
double volume(int dim)
{
	return pow(M_PI,dim/2.)/gsl_sf_gamma(dim/2.+1);
}

int main(void)
{
	// prepare GSL random number generator
	const gsl_rng_type* T;
	gsl_rng* r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	double estim, err;

	fprintf(stdout,"# Volume of n-dim ball:\n");
	fprintf(stdout,"# DIM\tActual Vol.\tADAPT(Estim.)\tADAPT(counts)\tMC(Estim.)\tMC(Counts)\n");
	for(dim=1; dim < 7; dim++)
	{
		fprintf(stdout,"%d\t%g",dim,volume(dim));
		counter = 0;
		estim = adapt_nd_speclim(f,d,u,dim,ACC,EPS,&err);
		fprintf(stdout,"\t%g\t%llu",estim,counter);
		counter = 0;
		estim = monte_carlo_speclim(f,d,u,dim,ACC,EPS,100,&err,r);
		fprintf(stdout,"\t%g\t%llu\n",estim,counter);
	}

	fprintf(stdout,"\n\n# Integral of exp(sum_i x_i) in n-ball\n");
	fprintf(stdout,"# Dim\tADAPT\tMonte Carlo\n");

	for(dim=1; dim < 7; dim++)
	{
		fprintf(stdout,"%d",dim);
		counter = 0;
		estim = adapt_nd_speclim(fexp,d,u,dim,ACC,EPS,&err);
		fprintf(stdout,"\t%llu",counter);
		counter = 0;
		estim = monte_carlo_speclim(fexp,d,u,dim,ACC,EPS,100,&err,r);
		fprintf(stdout,"\t%llu\n",counter);
	}

	gsl_rng_free(r);

	return 0;
}
