#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef struct
{
	int n;
	double *x, *y, *b, *c;
} qspline;

// binary search algorithm
int qspline_binsearch(qspline * s, double z)
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
qspline* qspline_alloc(int n, double x[], double y[])
{
	assert(n>1);

	// allocating memory
	qspline* s = (qspline*) malloc(n*sizeof(qspline));
	s->n = n;
	s->x = (double*) malloc( n*sizeof(double) );
	s->y = (double*) malloc( n*sizeof(double) );
	s->b = (double*) malloc( (n-1)*sizeof(double) );
	s->c = (double*) malloc( (n-1)*sizeof(double) );
	
	// transfering x and y into s
	int i;
	for(i=0; i<n; i++)
	{
		s->x[i] = x[i];
		s->y[i] = y[i];
	}

	// estimating p[i]
	double p[n-1], h[n-1];
	for(i=0; i<n-1; i++)
	{
		h[i] = x[i+1] - x[i];
		p[i] = (y[i+1]-y[i])/h[i];
	}

	// forward recursion with c[0] = 0
	s->c[0] = 0;
	for(i = 1; i<n-1; i++)
	{
		s->c[i] = (p[i]-p[i-1]-s->c[i-1]*h[i-1])/h[i];
	}

	// backward recursion from 1/2*c[n-2]
	s->c[n-2] *= 0.5;
	for(i=n-2; i>0; i--)
	{
		s->c[i-1] = (p[i]-p[i-1]-s->c[i]*h[i])/h[i-1];
	}

	// estimating b[i]
	for(i=0; i<n-1; i++)
	{
		s->b[i] = p[i] - s->c[i]*h[i];
	}

	return s;
}

// evaluates the quadratic spline
double qspline_eval(qspline * s, double z)
{
	assert(s->x[0] <= z && z <= s->x[s->n-1]);
	
	// find the correct coordinate
	int i = qspline_binsearch(s,z);

	return s->y[i] + s->b[i]*(z - s->x[i]) + s->c[i]*(z - s->x[i])*(z - s->x[i]);
}

// frees the memory
void qspline_free(qspline * s)
{
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s);
}

// evaluate the derivative of the quadratic spline
double qspline_derivative(qspline * s, double z)
{
	assert(s->x[0] <= z && z <= s->x[s->n-1]);
	
	// search for the correct coordinate
	int i = qspline_binsearch(s,z);

	return s->b[i] + 2.0*s->c[i]*(z - s->x[i]);
}

// evaluate the integral from x[0] to z
double qspline_integral(qspline *s, double z)
{
	assert(s->x[0] <= z && z <= s->x[s->n-1]);
	
	// binary search for the correct coordinate, i
	int i = qspline_binsearch(s,z);

	// the integral for x[0] to x[i]
	double h, zint = 0;
	int j;
	for(j=0; j<i; j++)
	{
		h = s->x[j+1] - s->x[j];
		zint += s->y[j]*h + s->b[j]*h*h/2.0 + s->c[j]*h*h*h/3.0;
	}
	
	// the integral for x[i] to z
	h = z - s->x[i];
	zint += s->y[i]*h + s->b[i]*h*h/2.0 + s->c[i]*h*h*h/3.0;

	return zint;
}
