#include <stdio.h>
#include <math.h>
#include "monte_carlo.h"

// f1 is spherical volume integral r^2*sin(theta), where
// the coordinates (r, phi, theta) are to be integrated
// from 0 -> 1, 0 -> 2*pi, 0 -> pi
// the volume is then 4/3*pi
double f1(double* t)
{
	return t[0]*t[0]*sin(t[2]);
}

// integration of cos(r) over the unit circle
// due to symmetry the integral is exactly zero
// the vector is t = (r,theta)
double f2(double* t)
{
	return cos(t[0])*t[0]*sin(t[1]);
}

// f3, integral from assignment
double f3(double *t)
{
	return 1./(1.-cos(t[0])*cos(t[1])*cos(t[2]));
}
int main(void)
{
	// prepare GSL random number generator
	const gsl_rng_type* T;
	gsl_rng* r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	int N;
	double a[3], b[3], xct, estim, err;

	// f1, Volume of sphere
	a[0] = 0; a[1] = 0, a[2] = 0;
	b[0] = 1; b[1] = 2*M_PI; b[2] = M_PI;
	xct = 4./3.*M_PI;
	fprintf(stdout,"----Volume of sphere----\nExact volume:\t%g\n\n",xct);
	for (N=10; N<=1e8; N *= 10)
	{
		estim = monte_carlo_plain(f1,a,b,3,N,&err,r);
		fprintf(stdout,"N = %d\n",N);
		fprintf(stdout,"Estimate:\t%g\n",estim);
		fprintf(stdout,"Estim. error:\t%g\n",err);
		fprintf(stdout,"True error:\t%g\n\n",estim-xct);
	}

	// f3, 2D integral of cos(r) over unit circle
	// note that the true value is zero
	// the estimate is thus also the deviation from correct value
	// The integration of f2 is used to test the power law of the error
	fprintf(stdout,"\n----Integral of cos(r) over unit circle----\nExact value:\t0\n\n");
	for(N=10; N<=1e9; N*= 10)
	{
		estim = monte_carlo_plain(f2,a,b,2,N,&err,r);
		fprintf(stdout,"N = %d\n",N);
		fprintf(stdout,"Estimate:\t%g\n",estim);
		fprintf(stdout,"Estim. error:\t%g\n",err);
		fprintf(stderr,"%d\t%g\n",N,fabs(estim));
	}

	// f3, integral from assignment
	b[0] = M_PI; b[1] = M_PI; b[2] = M_PI;
	N = 1e9;
	xct = 1.3932039296856768591842462603255;
	fprintf(stdout,"\n----Assignment----\nTrue results:\t%g\n\n",xct);
	fprintf(stdout,"N = %d\n",N);
	estim = monte_carlo_plain(f3,a,b,3,N,&err,r)/M_PI/M_PI/M_PI;
	fprintf(stdout,"Estimate:\t%g\n",estim);
	fprintf(stdout,"Estim. error:\t%g\n",err);
	fprintf(stdout,"True error:\t%g\n\n",estim-xct);

	gsl_rng_free(r);
	return 0;
}
