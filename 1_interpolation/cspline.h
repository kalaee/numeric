#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

typedef struct
{
	int n;
	double *x, *y, *b, *c, *d, *h;
} cspline;

// binary search algorithm
int cspline_binsearch(cspline * s, double z)
{
	int i = 0, j = s->n-1, k;
	while(j-i > 1)
	{
		k = (i+j)/2;
		if (z > s->x[k])
			i = k;
		else
			j = k;
	}

	return i;
}

// construct the quadratic spline
cspline* cspline_alloc(int n, double x[], double y[], double ds0, double dds0, double dsn, double ddsn)
{
	assert(n>1);

	// allocating memory
	cspline* s = (cspline*) malloc(n*sizeof(cspline));
	s->n = n;
	s->x = (double*) malloc( n*sizeof(double) );
	s->y = (double*) malloc( n*sizeof(double) );
	s->b = (double*) malloc( (n-1)*sizeof(double) );
	s->c = (double*) malloc( (n-1)*sizeof(double) );
	s->d = (double*) malloc( (n-1)*sizeof(double) );
	s->h = (double*) malloc( (n-1)*sizeof(double) );
	
	// transfering x and y into s
	int i;
	for(i=0; i<n; i++)
	{
		s->x[i] = x[i];
		s->y[i] = y[i];
	}

	// estimating p[i] and h[i]
	double p[n-1];
	for(i=0; i<n-1; i++)
	{
		s->h[i] = x[i+1] - x[i];
		p[i] = (y[i+1]-y[i])/s->h[i];
	}

	// calculates b[i]
	s->b[0] = ds0;
	s->b[1] = 3.0*p[0] - dds0*s->h[0]/2.0 - 2.0*s->b[0];
	for(i = 2; i < n-1; i++)
	{
		s->b[i] = s->h[i-1]/s->h[i-2]*(3.0*(p[i-2]+p[i-1]*s->h[i-2]/s->h[i-1])-2.0*(s->h[i-2]/s->h[i-1]+1)*s->b[i-1]-s->b[i-2]);
	}

	// calculate c[i]
	s->c[0] = dds0/2.0;
	s->c[n-2] = (2.0*dsn-2.0*s->b[n-2]+s->h[n-2]*ddsn)/(2.0*s->h[n-2]);
	for(i = 1; i < n-2; i++)
	{
		s->c[i] = (-2.0*s->b[i]-s->b[i+1]+3.0*p[i])/s->h[i];
	}

	// calculate d[i]
	s->d[n-2] = (s->b[n-2]+s->h[n-2]*ddsn - dsn)/(3.0*pow(s->h[n-2],2));
	for(i = 0; i < n-2; i++)
	{
		s->d[i] = (s->b[i]+s->b[i+1]-2.0*p[i])*pow(s->h[i],-2);
	}

	return s;
}

// evaluates the quadratic spline
double cspline_eval(cspline * s, double z)
{
	assert(s->x[0] <= z && z <= s->x[s->n-1]);
	
	// find the correct coordinate
	int i = cspline_binsearch(s,z);
	double h = z-s->x[i];
	return s->y[i]+s->b[i]*h+s->c[i]*pow(h,2)+s->d[i]*pow(h,3);
}

// frees the memory
void cspline_free(cspline * s)
{
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s->d);
	free(s->h);
	free(s);
}

// evaluate the derivative of the quadratic spline
double cspline_derivative(cspline * s, double z)
{
	assert(s->x[0] <= z && z <= s->x[s->n-1]);
	
	// search for the correct coordinate
	int i = cspline_binsearch(s,z);
	double h = z - s->x[i];
	return s->b[i]+2.0*s->c[i]*h+3.0*s->d[i]*pow(h,2);
}

// evaluate the integral from x[0] to z
double cspline_integral(cspline *s, double z)
{
	assert(s->x[0] <= z && z <= s->x[s->n-1]);
	
	// binary search for the correct coordinate, i
	int i = cspline_binsearch(s,z);

	// the integral for x[0] to x[i]
	double h, zint = 0;
	int j;
	for(j=0; j<i; j++)
	{
		zint += s->y[j]*s->h[j] + s->b[j]*pow(s->h[j],2)/2.0
				+ s->c[j]*pow(s->h[j],3)/3.0 + s->d[j]*pow(s->h[j],4)/4.0;
	}
	
	// the integral for x[i] to z
	h = z - s->x[i];
	zint += s->y[i]*h + s->b[i]*pow(h,2)/2.0
			+ s->c[i]*pow(h,3)/3.0 + s->d[i]*pow(h,4)/4.0;

	return zint;
}
