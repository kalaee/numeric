#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

typedef struct
{
	int n;
	double *x, *y, *b, *c, *d;
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
cspline* cspline_alloc(int n, double x[], double y[])
{
	cspline * s = (cspline*) malloc(sizeof(cspline));
	s->n = n;
	s->x = (double*) malloc(n*(sizeof(double)));
	s->y = (double*) malloc(n*(sizeof(double)));
	s->b = (double*) malloc(n*(sizeof(double)));
	s->c = (double*) malloc(n*(sizeof(double)));
	s->d = (double*) malloc(n*(sizeof(double)));

	int i;
	for(i=0; i<n; i++)
	{
		s->x[i] = x[i];
		s->y[i] = y[i];
	}

	double h[n-1], p[n-1];
	for (i = 0; i<n; i++)
	{
		 h[i] = x[i+1]-x[i];
		 p[i] = (y[i+1]-y[i])/h[i];
	}

	double D[n], Q[n-1], B[n];
	D[0]=2;
	Q[0]=1;
	B[0]=3*p[0];
	for(i = 1; i<n-1; i++)
	{
		D[i] = 2*h[i-1]/h[i]+2;
		Q[i] = h[i-1]/h[i];
		B[i] = 3*(p[i-1]+p[i]*h[i-1]/h[i]);
	}
	D[n-1] = 2;
	B[n-1] = 3*p[n-1];

	for (i=1; i<n; i++)
	{
		D[i] -= Q[i-1]/D[i-1];
		B[i] -= B[i-1]/D[i-1];
	}
	s->b[n-1] = B[n-1]/D[n-1];
	for(i = n-2; i >= 0; i--)
	{
		s->b[i] =(B[i]-Q[i]*s->b[i+1]/D[i]);
	}

	for (i = 0; i<n-1; i++)
	{
		s->c[i] = (-2*s->b[i]-s->b[i+1]+3*p[i])/h[i];
		s->d[i] = (s->b[i]+s->b[i+1]-2*p[i])/h[i]/h[i];
	}

	return s;
}

// evaluates the quadratic spline
double cspline_eval(cspline * s, double z)
{
	assert(s->x[0] <= z && z <= s->x[s->n-1]);
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
	free(s);
}

// evaluate the derivative of the quadratic spline
double cspline_derivative(cspline * s, double z)
{
	assert(s->x[0] <= z && z <= s->x[s->n-1]);
	int i = cspline_binsearch(s,z);
	double h = z - s->x[i];
	return s->b[i]+2*s->c[i]*h+3*s->d[i]*pow(h,2);
}

// evaluate the integral from x[0] to z
double cspline_integral(cspline *s, double z)
{
	assert(s->x[0] <= z && z <= s->x[s->n-1]);
	int i = cspline_binsearch(s,z);

	// the integral for x[0] to x[i]
	double h, zint = 0;
	int j;
	for(j=0; j<i; j++)
	{
		h = s->x[j+1] - s->x[j];
		zint += s->y[j]*h + s->b[j]*pow(h,2)/2.0
				+ s->c[j]*pow(h,3)/3.0 + s->d[j]*pow(h,4)/4.0;
	}
	
	// the integral for x[i] to z
	h = z - s->x[i];
	zint += s->y[i]*h + s->b[i]*pow(h,2)/2.0
			+ s->c[i]*pow(h,3)/3.0 + s->d[i]*pow(h,4)/4.0;

	return zint;
}
